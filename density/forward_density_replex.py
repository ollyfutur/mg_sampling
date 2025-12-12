import itertools

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from matplotlib.colors import CenteredNorm
from MDAnalysis.analysis.align import alignto
from tqdm import tqdm

from dask import compute, delayed

from mg_sampling import PROJECT_DIR
from mg_sampling.density.align import (
    check_alignable,
    reorder_residue_atoms,
    standardize_names,
    unwrap_sequential,
    wrap_around_atomgroup,
)
from mg_sampling.density.forward_density import atomgroup_to_map, voxels_around
from mg_sampling.density.grid import Grid, SparseGrid

out_radius = 10
device = 'cuda'
chunk_size = 16384

# Reference structure to align to
ref = mda.Universe(PROJECT_DIR / 'density/9CBU.pdb')
# Reference CryoEM map (for map geometry)
refmap = Grid(PROJECT_DIR / 'density/emd_42499.mrc', file_format='mrc')

methods = ['replex', 'plain_md']
systems = ['7ez0', '7ez0_noion', '7ez0_nomg']
forcefields = ['mMg_tip3p', 'nMg_tip3p', 'nMg_tip3p_soft']
replicates = range(2)
repids = range(8)
start = 500  # after 50 ns equilibration

def run_plain_md(params):
    meth, sys, ff, rep = params

    path = PROJECT_DIR / f'{meth}/{sys}/{ff}/replicate{rep}'

    if not path.exists():
        print("Skip:", path)
        return

    outpath = PROJECT_DIR / f'density/{meth}_{sys}_{ff}'
    outpath.mkdir(exist_ok=True, parents=True)

    # Load all replicas
    paths = [path / f'rep{i+8*rep:02d}' for i in repids]
    u = mda.Universe(
        paths[0] / 'input.gro',
        *(p / 'run.xtc' for p in paths),
    )

    standardize_names(u)

    mob_align = reorder_residue_atoms(u.select_atoms('not type H and nucleic'))
    mob_align_light = mob_align.select_atoms('name P')[::5]  # For efficiently wrapping solvent around
    ref_align = reorder_residue_atoms(ref.select_atoms('not type H and nucleic and not (resid 22 and name P OP1 OP2)'))

    vox_indices = voxels_around(ref_align, 6, refmap.grid_params)

    check_alignable(mob_align, ref_align)

    mob_out = u.select_atoms("not type H")
    mob_out = mob_out.select_atoms(f"group mob_align or around {out_radius} group mob_align", mob_align=mob_align, updating=True)

    tot_map = 0
    for ts in tqdm(u.trajectory[start:]):
        unwrap_sequential(mob_align)
        wrap_around_atomgroup(mob_out, mob_align_light)
        alignto(mob_align, ref_align, subselection=mob_out)
        map_ = atomgroup_to_map(
            mob_out,
            shape=refmap.shape,
            delta=refmap.delta,
            origin=refmap.origin,
            vox_indices=vox_indices,
            device=device,
            chunk_size=chunk_size,
        )
        tot_map += map_
    tot_map /= len(u.trajectory[start:])
    tot_map.save(outpath / f'global_rep{i:02d}.pkl')

if __name__ == "__main__":
    tasks = [delayed(run_plain_md)(params) for params in itertools.product(methods, systems, forcefields, replicates)]
    compute(*tasks, scheduler='processes', num_workers=8)
