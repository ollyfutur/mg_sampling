import os
import itertools

import MDAnalysis as mda
from MDAnalysis.analysis.align import alignto

from dask import compute, delayed
from tqdm import tqdm
from mg_sampling import PROJECT_DIR
from mg_sampling.density.align import (
    check_alignable,
    reorder_residue_atoms,
    standardize_names,
    unwrap_sequential,
    wrap_around_atomgroup,
)
from mg_sampling.density.forward_density import atomgroup_to_map, voxels_around
from mg_sampling.density.grid import Grid


def run(
    top_file,
    traj_file,
    ref_file,
    out_file,
    grid_params,
    device="cuda",
    start=None,
    end=None,
    step=None,
    chunk_size=None,
):
    print("Computing:", top_file, traj_file)

    u = mda.Universe(
        top_file,
        traj_file,
    )
    ref = mda.Universe(ref_file)

    standardize_names(u)
    standardize_names(ref)

    mob_align = reorder_residue_atoms(u.select_atoms("not type H and nucleic"))
    # For efficiently wrapping solvent around
    mob_align_light = mob_align.select_atoms("name P")[::5]
    ref_align = reorder_residue_atoms(
        ref.select_atoms("not type H and nucleic and not (resid 22 and name P OP1 OP2)")
    )

    vox_indices = voxels_around(ref_align, 6, grid_params)

    check_alignable(mob_align, ref_align)

    mob_out = u.select_atoms("not type H")
    mob_out = mob_out.select_atoms(
        "group mob_align or around 10 group mob_align",
        mob_align=mob_align,
        updating=True,
    )

    tot_map = 0
    n_frames = len(u.trajectory[start:end:step])
    for ts in tqdm(u.trajectory[start:end:step], miniters=10):
        unwrap_sequential(mob_align)
        wrap_around_atomgroup(mob_out, mob_align_light)
        alignto(mob_align, ref_align, subselection=mob_out)
        map_ = atomgroup_to_map(
            mob_out,
            shape=grid_params["shape"],
            delta=grid_params["delta"],
            origin=grid_params["origin"],
            vox_indices=vox_indices,
            device=device,
            chunk_size=chunk_size,
        )
        tot_map += map_
    tot_map /= n_frames

    print("Saving output map:", out_file, flush=True)
    tot_map.save(out_file)


if __name__ == "__main__":
    device = "cuda"
    chunk_size = 16384
    num_workers = int(os.environ["SLURM_CPUS_PER_TASK"])

    # Reference structure to align to
    ref_file = PROJECT_DIR / "density/9CBU.pdb"
    # Reference CryoEM map (for map geometry)
    grid_params = Grid(
        PROJECT_DIR / "density/emd_42499.mrc", file_format="mrc"
    ).grid_params

    systems = ["7ez0", "7ez0_noion", "7ez0_nomg"]
    forcefields = ["mMg_tip3p", "nMg_tip3p", "nMg_tip3p_soft"]
    replicates = range(8)
    start = 1500  # after 150 ns equilibration

    tasks = []
    for sys, ff, rep in itertools.product(systems, forcefields, replicates):
        path = PROJECT_DIR / f"plain_md/{sys}/{ff}/rep{rep:02d}"

        if not path.exists():
            print("Skip:", path)
            continue

        outpath = PROJECT_DIR / f"density/plain_md_{sys}_{ff}"
        outpath.mkdir(exist_ok=True, parents=True)

        top_file = path / "input.gro"
        traj_file = path / "run.xtc"
        out_file = outpath / f"rep{rep:02d}.pkl"

        tasks.append(
            delayed(run)(
                top_file,
                traj_file,
                ref_file,
                out_file,
                grid_params,
                start=start,
                device=device,
                chunk_size=chunk_size,
            )
        )

    print(f"Computing {len(tasks)} tasks…", flush=True)
    compute(*tasks, scheduler="processes", num_workers=num_workers)
