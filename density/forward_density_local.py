import itertools
import numpy as np
import MDAnalysis as mda

from dask import delayed, compute
from dask.distributed import Client

from mg_sampling import PROJECT_DIR
from mg_sampling.tools import read_plumed
from mg_sampling.density.align import (
    TransformationFromFunction,
    local_align,
    check_alignable,
    reorder_residue_atoms,
    standardize_names,
    unwrap_sequential,
    wrap_around_atomgroup,
)
from mg_sampling.density.forward_density import voxels_around, atomgroup_to_map
from mg_sampling.density.grid import Grid

if __name__ == "__main__":
    device = "cuda"
    chunk_size = 8192
    n_workers = 6

    # Reference CryoEM map (for map geometry)
    grid_params = Grid(
        PROJECT_DIR / "density/emd_42499.mrc", file_format="mrc"
    ).grid_params

    # Reference structure to align to
    ref_file = PROJECT_DIR / "structures/9CBU_fixed.pdb"
    ref = mda.Universe(ref_file)
    ref_align = reorder_residue_atoms(
        ref.select_atoms("not type H and nucleic and not (resid 22 and name P OP1 OP2)")
    )
    vox_indices = voxels_around(ref_align, 6, grid_params)

    # Iterate
    # methods = ["plain_md", "replex"]
    # systems = ["7ez0", "7ez0_noion", "7ez0_nomg"]
    # forcefields = ["mMg_tip3p", "nMg_tip3p", "nMg_tip3p_soft"]
    methods = ["plain_md"]
    systems = ["7ez0_noion"]
    forcefields = ["mMg_tip3p"]
    replicates = range(2)
    repids = range(8)
    out_sels = {'all': "not type H"}

    def run(method, sys, ff, rep, i, out):
        path = PROJECT_DIR / f"{method}/{sys}/{ff}/replicate{rep}/rep{i+rep*8:02d}"

        if not path.is_dir():
            print("Ignoring", path)
            return

        outpath = PROJECT_DIR / f"density/{method}_{sys}_{ff}/local/{out}"
        outpath.mkdir(exist_ok=True, parents=True)

        top_file = path / "input.gro"
        traj_file = path / "run.xtc"
        out_file = outpath / f"replicate{rep}_rep{i:02d}.pkl"

        if out_file.exists():
            print("Exists already:", out_file)
            return

        u = mda.Universe(top_file, traj_file)

        standardize_names(u)
        mob_align = reorder_residue_atoms(u.select_atoms("not type H and nucleic"))
        check_alignable(mob_align, ref_align)

        # For efficiently wrapping solvent around
        mob_align_light = mob_align.select_atoms("name P")[::5]

        mob_out = u.select_atoms(out_sels[out])
        mob_out_or_align = mob_out | mob_align
        mob_out_sub = mob_out_or_align.select_atoms(
            "group mob_out and (group mob_align or around 10 group mob_align)",
            mob_align=mob_align,
            mob_out=mob_out,
            updating=True,
        )

        workflow = [
            TransformationFromFunction(unwrap_sequential, mob_align),
            TransformationFromFunction(wrap_around_atomgroup, mob_out_sub, mob_align_light),
        ]
        u.trajectory.add_transformations(*workflow)

        if method == 'replex':
            reptemp = read_plumed(path / f'bias.COLVAR.{i}')['bias.index'].astype(int)
            frames = (reptemp == 0).to_numpy()[::10]
            assert len(frames) == len(u.trajectory)
            sel = np.zeros_like(frames, dtype=bool)
            # Skip 1/4 = 50 ns as equilibration
            # Total frames per map: 1500 * 2 replicates ~ 3000
            sel[500:2000:1] = True  
            frames = frames & sel
        elif method == 'plain_md':
            sel = np.zeros(len(u.trajectory), dtype=bool)
            if sys == '7ez0_nomg':
                # Skip first 1/4 = 100 ns
                # Subsample every 8 frames
                # Total frames per map: 3000 ÷ 8 * 8 ~ 3000
                sel[1000:4000:8] = True
            else:
                # Skip first 1/4 = 150 ns
                # Subsample every 12 frames
                # Total frames per map: 4500 ÷ 12 * 8 ~ 3000
                sel[1000:4000:8] = True
            frames = sel
        print(f"Selected {int(frames.sum())} frames from {len(u.trajectory)}")

        tot_map = 0
        for ts in u.trajectory[frames]:
            new_u = local_align(
                mob_out_or_align,
                ref_align,
                align_sel='nucleic',
                align_radius=10,
                out_radius=10,
                reorder=True,
            )
            out_file_ts = out_file.with_stem(out_file.stem + f"_{ts.frame}")
            map_ = atomgroup_to_map(
                new_u.atoms,
                shape=grid_params['shape'],
                delta=grid_params['delta'],
                origin=grid_params['origin'],
                vox_indices=vox_indices,
                cutoff=4.0,
                device=device,
                chunk_size=chunk_size,
            )
            tot_map += map_
            if ts.frame % 100 == 0:
                with mda.Writer(out_file_ts.with_suffix(".pdb")) as w:
                    w.write(new_u.atoms)
                    del w.obj
                map_.save(out_file_ts.with_suffix(".pkl"))

        tot_map /= len(u.trajectory[frames])

        print("Saving output map:", out_file, flush=True)
        tot_map.save(out_file)

        # Save total number of frames as weight for later analysis
        weight_file = out_file.with_stem(out_file.stem + "_weight").with_suffix(".txt")
        with open(weight_file, "w") as f:
            f.write(f"{frames.sum()}")

    tasks = []
    for method, sys, ff, rep, i, out in itertools.product(methods, systems, forcefields, replicates, repids, out_sels):
        tasks.append(
            delayed(run)(method, sys, ff, rep, i, out)
        )

    with Client(threads_per_worker=1, n_workers=n_workers) as client:
        print(f"Computing {len(tasks)} tasks…", flush=True)
        print("Dashboard:", client.dashboard_link)
        compute(*tasks)
