import itertools
import numpy as np
import MDAnalysis as mda
from MDAnalysis.transformations import fit_rot_trans
from MDAnalysis.transformations.base import TransformationBase

from mg_sampling import PROJECT_DIR
from mg_sampling.tools import read_plumed
from mg_sampling.density.align import (
    check_alignable,
    reorder_residue_atoms,
    standardize_names,
    unwrap_sequential,
    wrap_around_atomgroup,
)
from mg_sampling.density.forward_density import voxels_around, ForwardDensity
from mg_sampling.density.grid import Grid

class TransformationFromFunction(TransformationBase):
    def __init__(self, func, *args, max_threads=1, parallelizable=True, **kwargs):
        super().__init__(max_threads=max_threads,
                          parallelizable=parallelizable)
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def _transform(self, ts):
        self.func(*self.args, **self.kwargs)
        return ts

if __name__ == "__main__":
    device = "cuda"
    chunk_size = 8192
    n_workers = None
    backend = 'serial'

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
    #systems = ["7ez0", "7ez0_noion"]
    systems = ["7ez0"]
    forcefields = ["mMg_tip3p"]
    replicates = range(2)
    # replicates = [0]
    repids = range(8)
    # repids = [0]
    start = 500  # after 50 ns equilibration
    stop = 2000
    step = 1  # 1500 frames * 2 replicates ~ 3000 frames

    out_sel = {'all': "not type H", 'mg': "name MG"}

    tasks = []
    for sys, ff, rep, i, out in itertools.product(systems, forcefields, replicates, repids, out_sel):
        path = PROJECT_DIR / f"replex/{sys}/{ff}/replicate{rep}/rep{i+rep*8:02d}"

        outpath = PROJECT_DIR / f"density/replex_{sys}_{ff}/{out}"
        outpath.mkdir(exist_ok=True, parents=True)

        top_file = path / "input.gro"
        traj_file = path / "run.xtc"
        out_file = outpath / f"replicate{rep}_rep{i:02d}.pkl"

        u = mda.Universe(top_file, traj_file)

        standardize_names(u)
        mob_align = reorder_residue_atoms(u.select_atoms("not type H and nucleic"))
        check_alignable(mob_align, ref_align)

        # For efficiently wrapping solvent around
        mob_align_light = mob_align.select_atoms("name P")[::5]

        mob_out = u.select_atoms(out_sel[out])
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
            fit_rot_trans(mob_align, ref_align),
        ]
        u.trajectory.add_transformations(*workflow)

        reptemp = read_plumed(path / f'bias.COLVAR.{i}')['bias.index'].astype(int)
        frames = (reptemp == 0).to_numpy()[::10]
        assert len(frames) == len(u.trajectory)
        sel = np.zeros_like(frames, dtype=bool)
        sel[start:stop:step] = True
        frames = frames & sel
        print(f"Selected {int(frames.sum())} frames from {len(u.trajectory)}")

        anal = ForwardDensity(
            mob_out_sub,
            shape=grid_params['shape'],
            delta=grid_params['delta'],
            origin=grid_params['origin'],
            out_sel="around 10 (not type H and nucleic)",
            vox_indices=vox_indices,
            cutoff=4.0,
            device=device,
            chunk_size=chunk_size,
        )

        verbose = True if backend == 'serial' else False
        anal.run(frames=frames, n_workers=n_workers, backend=backend, verbose=verbose)

        print("Saving output map:", out_file, flush=True)
        anal.results.map.save(out_file)

        weight_file = out_file.with_stem(out_file.stem + "_weight").with_suffix(".txt")
        with open(weight_file, "w") as f:
            f.write(f"{frames.sum()}")


