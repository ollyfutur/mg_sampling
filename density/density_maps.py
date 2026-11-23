"""
Forward density maps of one replicate of a simulated system, aligned onto the reference
structure, as a cryomap `GridEnsemble`:
- rigidly (whole RNA), written to `density/global/<name>`;
- locally, with `--align-radius` and `--out-radius`: each residue is fitted on the RNA atoms
  within `align-radius` of it and carries the atoms within `out-radius` of it (see cryomap's
  `ForwardDensity`), written to `density/local/align<A>_out<O>/<name>`.

Examples:
    python density_maps.py plain_md/7ez0/mMg_tip3p 3     # replicate0/rep03
    python density_maps.py replex/7ez0/mMg_tip3p 1       # replicate1 (rep08-15)
    python density_maps.py plain_md/7ez0/mMg_tip3p 3 --align-radius 10 --out-radius 10
    python density_maps.py replex/7ez0/mMg_tip3p 0 --test 20 --out-dir /tmp/test

Existing outputs are skipped (use --overwrite to recompute).
"""

import argparse
import os
import shutil
import time
from pathlib import Path

import MDAnalysis as mda
from cryomap import ForwardDensity, masks, read_mrc
from cryomap.transformations import unwrap_consecutive, wrap_around

from mg_sampling.replex import ParallelReader
from mg_sampling.structure import (
    check_alignable,
    reorder_residue_atoms,
    standardize_names,
)
from pipeline import PROJECT_DIR
from pipeline.systems import SYSTEMS, is_replex, replica_dirs

REF_MAP = PROJECT_DIR / "density" / "emd_42499.mrc"
REF_PDB = PROJECT_DIR / "structures" / "9CBU_fixed.pdb"
RSHELL = 5.0  # Å (solvent shell around RNA)
REF_ALIGN_SEL = "nucleic and not element H and not (resid 22 and name P OP1 OP2)"
MOB_ALIGN_SEL = "nucleic and not element H"
OUT_SEL = "all"
BACKEND = "dask"  # MDAnalysis parallel backend
DEVICE = "cuda"  # torch device
CHUNK_SIZE = 4096  # atoms per GPU chunk

# Frame subsampling per system (production window: pipeline.systems.SYSTEMS)
STEP = {
    "replex/7ez0/mMg_tip3p": 1,
    "replex/7ez0_noion/mMg_tip3p": 1,
    "plain_md/7ez0/mMg_tip3p": 8,
    "plain_md/7ez0/nMg_tip3p": 8,
    "plain_md/7ez0/nMg_tip3p_soft": 8,
    "plain_md/7ez0_noion/mMg_tip3p": 8,
    "plain_md/7ez0_noion/nMg_tip3p": 8,
    "plain_md/7ez0_nomg/mMg_tip3p": 5,
}


def inputs(system, k):
    """Topology, trajectories and bias index file (replex only) of replicate k."""
    reps = replica_dirs(system, k)
    bias_index = reps[0] / "bias_index.xvg" if is_replex(system) else None
    return reps[0] / "run.tpr", [rep / "run.xtc" for rep in reps], bias_index


def frames(system, dt):
    """Frame range (start, stop, step) of the production window, for a trajectory saved every `dt` ps."""
    window = SYSTEMS[system]
    return round(window.start * 1000 / dt), round(window.end * 1000 / dt), STEP[system]


def output_name(system, k):
    return "_".join(system.split("/")) + f"_{k}"


def default_out_dir(align_radius, out_radius):
    if align_radius is None:
        return PROJECT_DIR / "density" / "global"
    return PROJECT_DIR / "density" / "local" / f"align{align_radius:g}_out{out_radius:g}"


def reference():
    """Reference alignment atoms and sparse grid parameters (RNA + solvent shell)."""
    ref_u = mda.Universe(REF_PDB, to_guess=("elements",))
    standardize_names(ref_u)
    ref_rna = ref_u.select_atoms("nucleic")
    ref_align = reorder_residue_atoms(ref_u.select_atoms(REF_ALIGN_SEL))
    params = read_mrc(REF_MAP).params
    params = params.sparse(
        mask=masks.within(params, ref_rna.positions, masks.vdw_radii(ref_rna) + RSHELL)
    )
    return ref_align, params


def compute(system, k, out_dir, args):
    name = output_name(system, k)
    outpath = out_dir / name
    if outpath.exists() and not args.overwrite:
        print(f"[{name}] exists, skipping: {outpath}", flush=True)
        return
    topfile, trajfiles, bias_index = inputs(system, k)
    missing = [f for f in [topfile, *trajfiles, bias_index] if f and not f.exists()]
    if missing:
        raise SystemExit(f"[{name}] missing input files: {[str(f) for f in missing]}")

    t0 = time.time()
    ref_align, params = reference()
    replex = is_replex(system)
    u = mda.Universe(
        topfile,
        trajfiles,
        format=ParallelReader if replex else None,
        to_guess=("elements",),
    )
    if replex:
        # This is in fact "re-multiplexing" continuous trajs to discontinuous, but constant bias trajs
        u.trajectory.demux(bias_index)
    standardize_names(u)
    mob_align = reorder_residue_atoms(u.select_atoms(MOB_ALIGN_SEL))
    check_alignable(mob_align, ref_align)
    mob_align_light = mob_align.select_atoms("name P")[::5]
    mob_out = u.select_atoms(OUT_SEL)

    # Whole RNA, other atoms wrapped around it; the fit itself is done by ForwardDensity
    u.trajectory.add_transformations(
        unwrap_consecutive(mob_align),
        wrap_around(mob_out, mob_align_light),
    )

    start, stop, step = frames(system, u.trajectory.dt)
    if args.test is not None:
        stop = min(stop, start + args.test * step)
    n_frames = len(range(start, min(stop, len(u.trajectory)), step))
    mode = "rigid" if args.align_radius is None else f"local (align {args.align_radius:g} Å, out {args.out_radius:g} Å)"
    print(
        f"[{name}] {n_frames} frames ({start}:{stop}:{step}), {mode} alignment, n_workers={args.n_workers}",
        flush=True,
    )

    ana = ForwardDensity(
        mob_out,
        params,
        align=mob_align,
        reference=ref_align,
        align_radius=args.align_radius,
        out_radius=args.out_radius,
        ensemble=True,
        backend="torch",
        device=DEVICE,
        chunk_size=CHUNK_SIZE,
    )
    ana.run(
        start=start, stop=stop, step=step, backend=BACKEND, n_workers=args.n_workers
    )

    ensemble = ana.results.ensemble
    if replex:
        replica_index = u.trajectory.replica_temp[ana.frames, 0]
        ensemble = ensemble.with_frame_info(replica_index=replica_index)

    # Write to a temporary directory first, so that an interrupted job (e.g. time
    # limit) never leaves a partial output that would be skipped on the next run
    tmp = outpath.with_name(outpath.name + ".tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    ensemble.write(tmp)
    shutil.rmtree(outpath, ignore_errors=True)
    tmp.rename(outpath)
    print(f"[{name}] wrote {outpath} in {time.time() - t0:.0f} s", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "system",
        choices=SYSTEMS,
        metavar="system",
        help=f"One of: {', '.join(SYSTEMS)}",
    )
    parser.add_argument(
        "replicate", type=int, help="Replicate number (replex: 0-1, plain MD: 0-7)"
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")),
        help="Number of dask worker processes (default: $SLURM_CPUS_PER_TASK, or 4)",
    )
    parser.add_argument(
        "--align-radius",
        type=float,
        help="Local alignment: fit each residue on the RNA atoms within this radius (Å)",
    )
    parser.add_argument(
        "--out-radius",
        type=float,
        help="Local alignment: each residue carries the atoms within this radius (Å)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Output directory (default: density/global, or density/local/align<A>_out<O>)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Recompute existing output"
    )
    parser.add_argument(
        "--test", type=int, metavar="N", help="Only compute the first N frames"
    )
    args = parser.parse_args()
    if (args.align_radius is None) != (args.out_radius is None):
        parser.error("--align-radius and --out-radius go together")

    args.out_dir = args.out_dir or default_out_dir(args.align_radius, args.out_radius)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    compute(args.system, args.replicate, args.out_dir, args)


if __name__ == "__main__":
    main()
