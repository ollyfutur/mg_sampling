"""
Forward density maps of one replicate of a simulated system, globally aligned onto the
reference structure, written to `density/global/<name>` (cryomap `GridEnsemble`).

Examples:
    python density_maps.py plain_md/7ez0/mMg_tip3p 3     # replicate0/rep03
    python density_maps.py replex/7ez0/mMg_tip3p 1       # replicate1 (rep08-15)
    python density_maps.py replex/7ez0/mMg_tip3p 0 --test 20 --out-dir /tmp/test

Existing outputs are skipped (use --overwrite to recompute).
"""

import argparse
import os
import shutil
import time
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.transformations import fit_rot_trans

from cryomap import ForwardDensity, masks, read_mrc
from cryomap.transformations import unwrap_consecutive, wrap_around
from mg_sampling import PROJECT_DIR
from mg_sampling.density.align import (
    check_alignable,
    reorder_residue_atoms,
    standardize_names,
)
from mg_sampling.parallel import ParallelReader

REF_MAP = PROJECT_DIR / "density" / "emd_42499.mrc"
REF_PDB = PROJECT_DIR / "structures" / "9CBU_fixed.pdb"
RSHELL = 5.0  # Å (solvent shell around RNA)
REF_ALIGN_SEL = "nucleic and not element H and not (resid 22 and name P OP1 OP2)"
MOB_ALIGN_SEL = "nucleic and not element H"
OUT_SEL = "all"
BACKEND = "dask"  # MDAnalysis parallel backend
DEVICE = "cuda"   # torch device
CHUNK_SIZE = 4096  # atoms per GPU chunk

# Frames (start, stop, step) per system, after discarding equilibration
FRAMES = {
    "replex/7ez0/mMg_tip3p": (750, 3000, 1),  # 75-300 ns
    "replex/7ez0_noion/mMg_tip3p": (750, 3000, 1),
    "plain_md/7ez0/mMg_tip3p": (1500, 6000, 8),  # 150-600 ns
    "plain_md/7ez0/nMg_tip3p": (1500, 6000, 8),
    "plain_md/7ez0/nMg_tip3p_soft": (1500, 6000, 8),
    "plain_md/7ez0_noion/mMg_tip3p": (1500, 6000, 8),
    "plain_md/7ez0_noion/nMg_tip3p": (1500, 6000, 8),
    "plain_md/7ez0_nomg/mMg_tip3p": (1000, 4000, 5),  # 100-400 ns
}
N_REPLICAS = 8  # replicas per replex replicate


def inputs(system, k):
    """Topology, trajectories and bias index file (replex only) of replicate k."""
    path = PROJECT_DIR / system
    if system.startswith("replex"):
        reps = [path / f"replicate{k}" / f"rep{i:02d}" for i in range(N_REPLICAS * k, N_REPLICAS * (k + 1))]
        return reps[0] / "run.tpr", [rep / "run.xtc" for rep in reps], reps[0] / "bias_index.xvg"
    rep = path / "replicate0" / f"rep{k:02d}"
    return rep / "run.tpr", [rep / "run.xtc"], None


def output_name(system, k):
    return "_".join(system.split("/")) + f"_{k}"


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
    is_replex = bias_index is not None
    u = mda.Universe(
        topfile,
        trajfiles,
        format=ParallelReader if is_replex else None,
        to_guess=("elements",),
    )
    if is_replex:
        # This is in fact "re-multiplexing" continuous trajs to discontinuous, but constant bias trajs
        u.trajectory.demux(bias_index)
    standardize_names(u)
    mob_align = reorder_residue_atoms(u.select_atoms(MOB_ALIGN_SEL))
    check_alignable(mob_align, ref_align)
    mob_align_light = mob_align.select_atoms("name P")[::5]
    mob_out = u.select_atoms(OUT_SEL)

    u.trajectory.add_transformations(
        unwrap_consecutive(mob_align),
        wrap_around(mob_out, mob_align_light),
        fit_rot_trans(mob_align, ref_align),
    )

    start, stop, step = FRAMES[system]
    if args.test is not None:
        stop = min(stop, start + args.test * step)
    n_frames = len(range(start, min(stop, len(u.trajectory)), step))
    print(f"[{name}] {n_frames} frames ({start}:{stop}:{step}), n_workers={args.n_workers}", flush=True)

    ana = ForwardDensity(
        mob_out, params, ensemble=True, backend="torch", device=DEVICE, chunk_size=CHUNK_SIZE
    )
    ana.run(start=start, stop=stop, step=step, backend=BACKEND, n_workers=args.n_workers)

    ensemble = ana.results.ensemble
    if is_replex:
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
    parser.add_argument("system", choices=FRAMES, metavar="system", help=f"One of: {', '.join(FRAMES)}")
    parser.add_argument("replicate", type=int, help="Replicate number (replex: 0-1, plain MD: 0-7)")
    parser.add_argument(
        "--n-workers", type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)),
        help="Number of dask worker processes (default: $SLURM_CPUS_PER_TASK, or 4)",
    )
    parser.add_argument("--out-dir", type=Path, default=PROJECT_DIR / "density" / "global",
                        help="Output directory (default: density/global)")
    parser.add_argument("--overwrite", action="store_true", help="Recompute existing output")
    parser.add_argument("--test", type=int, metavar="N", help="Only compute the first N frames")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    compute(args.system, args.replicate, args.out_dir, args)


if __name__ == "__main__":
    main()
