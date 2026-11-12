"""
Forward density maps (one ensemble per setup), globally aligned onto the reference
structure, written to `density/global/<setup>` (cryomap `GridEnsemble` directory).

Examples:
    python density_maps.py --list                       # list setups (index and name)
    python density_maps.py replex_7ez0_mMg_tip3p_0      # one setup, by name
    python density_maps.py 3 --n-workers 8              # one setup, by index
    python density_maps.py                              # all setups, sequentially
    python density_maps.py 0 --test 20 --out-dir /tmp/test   # quick check on 20 frames

Setups whose output already exists are skipped (use --overwrite to recompute).
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

SETUPS = {
    # Replex (one per replicate)
    f"replex_{sys}_mMg_tip3p_{rep}": {
        "topfile": PROJECT_DIR / "replex" / sys / "mMg_tip3p" / f"replicate{rep}" / f"rep{8*rep:02d}" / "run.tpr",
        "trajfiles": [
            PROJECT_DIR / "replex" / sys / "mMg_tip3p" / f"replicate{rep}" / f"rep{i+8*rep:02d}" / "run.xtc"
            for i in range(8)
        ],
        "bias_index": PROJECT_DIR / "replex" / sys / "mMg_tip3p" / f"replicate{rep}" / f"rep{8*rep:02d}" / "bias_index.xvg",
        "start": 750,
        "stop": 3000,
        "step": 1,
    } for sys in ["7ez0", "7ez0_noion"] for rep in range(2)
}
SETUPS.update({
    # Plain MD (one per replicate)
    f"plain_md_{sys}_{ff}_{rep}": {
        "topfile": PROJECT_DIR / "plain_md" / sys / ff / "replicate0" / f"rep{rep:02d}" / "run.tpr",
        "trajfiles": [PROJECT_DIR / "plain_md" / sys / ff / "replicate0" / f"rep{rep:02d}" / "run.xtc"],
        "start": 1500,
        "stop": 6000,
        "step": 8,
    } for sys in ["7ez0", "7ez0_noion"] for ff in ["mMg_tip3p", "nMg_tip3p"] for rep in range(8)
})
SETUPS.update({
    f"plain_md_7ez0_nMg_tip3p_soft_{rep}": {
        "topfile": PROJECT_DIR / "plain_md" / "7ez0" / "nMg_tip3p_soft" / "replicate0" / f"rep{rep:02d}" / "run.tpr",
        "trajfiles": [PROJECT_DIR / "plain_md" / "7ez0" / "nMg_tip3p_soft" / "replicate0" / f"rep{rep:02d}" / "run.xtc"],
        "start": 1500,
        "stop": 6000,
        "step": 8,
    } for rep in range(8)
})
SETUPS.update({
    f"plain_md_7ez0_nomg_mMg_tip3p_{rep}": {
        "topfile": PROJECT_DIR / "plain_md" / "7ez0_nomg" / "mMg_tip3p" / "replicate0" / f"rep{rep:02d}" / "run.tpr",
        "trajfiles": [PROJECT_DIR / "plain_md" / "7ez0_nomg" / "mMg_tip3p" / "replicate0" / f"rep{rep:02d}" / "run.xtc"],
        "start": 1000,
        "stop": 4000,
        "step": 5,
    } for rep in range(8)
})
SETUP_NAMES = list(SETUPS)


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


def compute_setup(setup, ref_align, params, out_dir, args):
    config = SETUPS[setup]
    outpath = out_dir / setup
    if outpath.exists() and not args.overwrite:
        print(f"[{setup}] exists, skipping: {outpath}", flush=True)
        return
    missing = [f for f in [config["topfile"], *config["trajfiles"], config.get("bias_index")] if f and not Path(f).exists()]
    if missing:
        print(f"[{setup}] missing input files, skipping: {[str(f) for f in missing]}", flush=True)
        return

    t0 = time.time()
    is_replex = setup.startswith("replex")
    u = mda.Universe(
        config["topfile"],
        config["trajfiles"],
        format=ParallelReader if is_replex else None,
        to_guess=("elements",),
    )
    if is_replex:
        # This is in fact "re-multiplexing" continuous trajs to discontinuous, but constant bias trajs
        u.trajectory.demux(config["bias_index"])
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

    start, stop, step = config["start"], config["stop"], config["step"]
    if args.test is not None:
        stop = min(stop, start + args.test * step)
    n_frames = len(range(start, min(stop, len(u.trajectory)), step))
    print(f"[{setup}] {n_frames} frames ({start}:{stop}:{step}), backend={args.backend}, "
          f"n_workers={args.n_workers}", flush=True)

    ana = ForwardDensity(
        mob_out, params, ensemble=True, backend="torch", device=args.device, chunk_size=args.chunk_size
    )
    run_kwargs = {} if args.backend == "serial" else {"n_workers": args.n_workers}
    ana.run(start=start, stop=stop, step=step, backend=args.backend, **run_kwargs)

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
    print(f"[{setup}] wrote {outpath} in {time.time() - t0:.0f} s", flush=True)


def parse_setups(items):
    names = []
    for item in items:
        if item.isdigit():
            names.append(SETUP_NAMES[int(item)])
        elif item in SETUPS:
            names.append(item)
        else:
            raise SystemExit(f"Unknown setup '{item}'. Use --list to see available setups.")
    return names


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("setups", nargs="*", help="Setup names or indices (default: all)")
    parser.add_argument("--list", action="store_true", help="List setups and exit")
    parser.add_argument(
        "--n-workers", type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)),
        help="Number of dask worker processes (default: $SLURM_CPUS_PER_TASK, or 4)",
    )
    parser.add_argument("--backend", choices=["dask", "serial"], default="dask",
                        help="MDAnalysis parallelization backend (default: dask)")
    parser.add_argument("--device", default="cuda", help="Torch device (default: cuda)")
    parser.add_argument("--chunk-size", type=int, default=4096, help="Atoms per GPU chunk (default: 4096)")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_DIR / "density" / "global",
                        help="Output directory (default: density/global)")
    parser.add_argument("--overwrite", action="store_true", help="Recompute existing outputs")
    parser.add_argument("--test", type=int, metavar="N", help="Only compute the first N frames of each setup")
    args = parser.parse_args()

    if args.list:
        for i, name in enumerate(SETUP_NAMES):
            done = "done" if (args.out_dir / name).exists() else ""
            print(f"{i:3d}  {name:40s} {done}")
        return

    setups = parse_setups(args.setups) if args.setups else SETUP_NAMES
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ref_align, params = reference()
    for setup in setups:
        compute_setup(setup, ref_align, params, args.out_dir, args)


if __name__ == "__main__":
    main()
