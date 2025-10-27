import argparse

import numpy as np
from scipy.stats import linregress

import mrcfile
import MDAnalysis as mda

from .forward_density import get_map_params, voxels_around

def main():
    parser = argparse.ArgumentParser(
        description="Compare density maps.",
    )
    parser.add_argument('mapA', type=str, help='Input MRC map')
    parser.add_argument("mapB", type=str, help="Input MRC map")
    parser.add_argument('--ref', type=str, help="Reference structure used for voxel selection")
    parser.add_argument('--vox-indices', type=str, help="Voxel indices file name")
    parser.add_argument('--vox-around', nargs=2, metavar=('R', 'SEL'), help="Select voxels around R the SEL in reference structure", action='append')
    parser.add_argument('--vox-not-around', nargs=2, metavar=('R', 'SEL'), help="Exclude voxels around R the SEL in reference structure", action='append')
    parser.add_argument('--density-threshold', type=float, default=-np.inf)
    parser.add_argument("-v", "--verbose", action="store_true", help="Show progress")
    args = parser.parse_args()

    with mrcfile.open(args.mapA, permissive=True) as mrcA, mrcfile.open(args.mapB, mode="r+", permissive=True) as mrcB:
        mapA_params = get_map_params(mrcA)
        print("Input Cryo-EM map A:")
        print("  Filename:", args.mapA)
        print("  Parameters:")
        print("    crs:", mapA_params['filemap'])
        print("    x0:", mapA_params['x0'])
        print("    dx:", mapA_params['dx'])
        print("    nbin:", mapA_params['nbin'])
        print("    nbin_tot:", mapA_params['nbin_tot'])

        mapB_params = get_map_params(mrcB)
        print("Input Cryo-EM map B:")
        print("  Filename:", args.mapB)
        print("  Parameters:")
        print("    crs:", mapB_params['filemap'])
        print("    x0:", mapB_params['x0'])
        print("    dx:", mapB_params['dx'])
        print("    nbin:", mapB_params['nbin'])
        print("    nbin_tot:", mapB_params['nbin_tot'])

        if any(mapA_params[key] != mapB_params[key] for key in ('x0', 'dx', 'nbin')):
            raise ValueError("Maps must have identical parameters", mapA_params, mapB_params)
        map_params = mapA_params

        if args.vox_indices is not None:
            vox_indices = np.load(args.vox_indices)
            print("Voxel index file:", args.vox_indices)
        else:
            vox_indices = np.array([], dtype=int)

        if args.vox_around is not None:
            if args.ref is None:
                parser.error("--vox-around requires setting --ref")
            vox_ref = mda.Universe(args.ref)

            for vox_r, vox_sel in args.vox_around:
                vox_r = float(vox_r)
                vox_ag = vox_ref.select_atoms(vox_sel)
                _inds = voxels_around(vox_ag, vox_r, map_params)
                vox_indices = np.union1d(vox_indices, _inds)
                print(f"Selecting {len(_inds)} voxels around:")
                print("  Reference:", args.ref)
                print("  Selection:", vox_sel)
                print("  Radius:", vox_r)

        if args.vox_not_around is not None:
            if args.ref is None:
                parser.error("--vox-not-around requires setting --ref")
            vox_ref = mda.Universe(args.ref)

            for vox_r, vox_sel in args.vox_not_around:
                vox_r = float(vox_r)
                vox_ag = vox_ref.select_atoms(vox_sel)
                _inds = voxels_around(vox_ag, vox_r, map_params)
                assert _inds.dtype == int
                vox_indices = np.setdiff1d(vox_indices, _inds)
                print(f"Excluding {len(_inds)} voxels around:")
                print("  Reference:", args.ref)
                print("  Selection:", vox_sel)
                print("  Radius:", vox_r)

        orderA = np.argsort(mapA_params['filemap'])
        dataA = mrcA.data.transpose(orderA).ravel()[vox_indices]
        orderB = np.argsort(mapB_params['filemap'])
        dataB = mrcB.data.transpose(orderB).ravel()[vox_indices]

        thresh_mask = (dataA > args.density_threshold) * (dataB > args.density_threshold)
        if (~thresh_mask).sum() > 0:
            print(f"Excluding {thresh_mask.sum()} indices where density is below threshold {args.density_threshold} in A or B.")
            dataA = dataA[thresh_mask]
            dataB = dataB[thresh_mask]
            vox_indices = vox_indices[thresh_mask]

        assert vox_indices.dtype == int
        print("Final voxel number:", len(vox_indices))
        np.savetxt('vox_indices.dat', vox_indices)

        def summarize_data(name, data):
            print(f"Map {name} ({data.size} voxels) stats:")
            print(f"  Mean: {data.mean()}")
            print(f"  Min: {data.min()}")
            print(f"  Max: {data.max()}")
            print(f"  Std: {data.std()}")

        summarize_data("A", dataA)
        summarize_data("B", dataB)

        print("Linear fitting of map A to map B:")
        slope, intercept, r_value, p_value, std_err = linregress(dataA, dataB)
        print(f"  {slope = :.5f}")
        print(f"  {intercept = :.5f}")
        print(f"  {r_value = :.5f}")
        print(f"  {p_value = :.5f}")
        print(f"  {std_err = :.5f}")

        summarize_data("A_fitted", slope * dataA + intercept)

if __name__ == '__main__':
    main()

