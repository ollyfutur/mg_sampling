"""
Backmapping of density map from (an ensemble of) structures
Inspired from code by Max Bonomi (mbonomi@pasteur.fr)
"""

import argparse
import time
from collections import Counter
from collections.abc import Sequence
from functools import partial
from numbers import Number
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import psutil
import torch
import torch.multiprocessing as mp
from MDAnalysis import AtomGroup
from tqdm import tqdm

from mg_sampling.density.grid import Grid, SparseGrid
from mg_sampling.density.utils import voxels_around
from mg_sampling.tools import batched, vectorize_dict_lookup

π = np.pi

## vdw radii [Ang] dictionary
## from ./modules/cctbx_project/cctbx/eltbx/van_der_waals_radii.py
## based on this http://www.phenix-online.org/pipermail/phenixbb/2014-August/021077.html
VDW_RADII = {
    "C": 1.775,
    "O": 1.45,
    "N": 1.50,
    "S": 1.80,
    "P": 1.9,
    "F": 1.47,
    "NA": 2.27,
    "MG": 1.73,
    "K": 2.75,
    "CL": 1.75,
}
_get_vdw_radii = vectorize_dict_lookup(VDW_RADII)

## RPROBE and RSHRINK from Jiang and Brunger
RPROBE = 1.0
RSHRINK = 1.1

## 5-Gaussians parameters
# map of atom types to A and B coefficients of scattering factor
# f(s) = A * exp(-B*s**2)
# B is in Angstrom squared
# Elastic atomic scattering factors of electrons for neutral atoms
# and s up to 6.0 A^-1: as implemented in PLUMED
### SIGMA
# B_={}
# B_["C"]=np.array([0.1140, 1.0825, 5.4281, 17.8811, 51.1341])
# B_["O"]=np.array([0.0652, 0.6184, 2.9449,  9.6298, 28.2194])
# B_["N"]=np.array([0.0541, 0.5165, 2.8207, 10.6297, 34.3764])
# B_["S"]=np.array([0.0838, 0.7788, 4.3462, 15.5846, 44.63655])
### WEIGHT
# A_={}
# A_["C"]=np.array([0.0489, 0.2091, 0.7537, 1.1420, 0.3555])
# A_["O"]=np.array([0.0365, 0.1729, 0.5805, 0.8814, 0.3121])
# A_["N"]=np.array([0.0267, 0.1328, 0.5301, 1.1020, 0.4215])
# A_["S"]=np.array([0.0915, 0.4312, 1.0847, 2.4671, 1.0852])
# Elastic atomic scattering factors of electrons for neutral atoms
# and s up to 2.0 A^-1: as implemented in phenix
### SIGMA
SF_WEIGHT = {
    "C": np.array([0.0893, 0.2563, 0.7570, 1.0487, 0.3575]),
    "O": np.array([0.0974, 0.2921, 0.6910, 0.6990, 0.2039]),
    "N": np.array([0.1022, 0.3219, 0.7982, 0.8197, 0.1715]),
    "S": np.array([0.2497, 0.5628, 1.3899, 2.1865, 0.7715]),
    "P": np.array([0.2548, 0.6106, 1.4541, 2.3204, 0.8477]),
    "F": np.array([0.1083, 0.3175, 0.6487, 0.5846, 0.1421]),
    "MG": np.array([0.2314, 0.6866, 0.9677, 2.1882, 1.1339]),
    "K": np.array([0.4115, 1.4031, 2.2784, 2.6742, 2.2162]),
    "NA": np.array([0.2142, 0.6853, 0.7692, 1.6589, 1.4482]),
    "CL": np.array([0.2443, 0.5397, 1.3919, 2.0197, 0.6621]),
}
SF_SIGMA = {
    "C": np.array([0.2465, 1.7100, 6.4094, 18.6113, 50.2523]),
    "O": np.array([0.2067, 1.3815, 4.6943, 12.7105, 32.4726]),
    "N": np.array([0.2451, 1.7481, 6.1925, 17.3894, 48.1431]),
    "S": np.array([0.2681, 1.6711, 7.0267, 19.5377, 50.3888]),
    "P": np.array([0.2908, 1.8740, 8.5176, 24.3434, 63.2996]),
    "F": np.array([0.2057, 1.3439, 4.2788, 11.3932, 28.7881]),
    "MG": np.array([0.3278, 2.2720, 10.9241, 39.2898, 101.9748]),
    "K": np.array([0.3703, 3.3874, 13.1029, 68.9592, 194.4329]),
    "NA": np.array([0.3334, 2.3446, 10.0830, 48.3037, 138.2700]),
    "CL": np.array([0.2468, 1.5242, 6.1537, 16.6687, 42.3086]),
}
_get_sf_weight = vectorize_dict_lookup(SF_WEIGHT)
_get_sf_sigma = vectorize_dict_lookup(SF_SIGMA)


## the Gaussian in density (real) space is the FT of scattering factor
## f(r) = A * (pi/B)**1.5 * exp(-pi**2/B*r**2)
## B must be augmented by bfactors (bfact)
## B -> B+bfact/4.0
## we define the following 2 auxiliary variables:
## invs2 = 2.0 * pi**2 / ( B + bfact/4.0 )
## pref  = A * (pi/( B + bfact/4.0 ))**1.5
## and we write the forward model as:
## f(r) = pref * exp(-0.5 * invs2 * r**2)


def get_fmod_param(ag: AtomGroup):
    tempfactors = ag.tempfactors if hasattr(ag, "tempfactors") else np.zeros(len(ag))
    pref = (
        _get_sf_weight(ag.types)
        * (π / (_get_sf_sigma(ag.types) + tempfactors[:, None] / 4)) ** 1.5
    )
    invs2 = 2 * π**2 / (_get_sf_sigma(ag.types) + tempfactors[:, None] / 4)
    return pref, invs2


def get_vdw_radii(ag: AtomGroup):
    return _get_vdw_radii(ag.types)


def atomgroup_to_map(
    ag: AtomGroup,
    shape: tuple[int],
    delta: tuple[float] | float,
    origin: tuple[float] = (0, 0, 0),
    vox_indices: Sequence[int] | None = None,
    cutoff: float = 4.0,
    device: str = "cpu",
    chunk_size: int | None = None,
) -> Grid:
    """
    Return map from the atoms in `ag` on a grid provided with `map_params`.
    """
    if isinstance(delta, Number):
        delta = (delta, delta, delta)

    if chunk_size is None:
        chunk_size = len(ag)

    delta = torch.tensor(delta, dtype=torch.float, device=device)  # (3,)
    origin = torch.tensor(origin, dtype=torch.float, device=device)  # (3,)
    grid = torch.zeros(shape, dtype=torch.float, device=device)  # (X, Y, Z)

    if vox_indices is not None:
        vox_indices = torch.tensor(vox_indices, dtype=torch.int64, device=device)

    a_xyz = torch.tensor(ag.positions, dtype=torch.float, device=device)  # (N, 3)
    if hasattr(ag, "occupancies"):
        a_occ = torch.tensor(ag.occupancies, dtype=torch.float, device=device)  # (N,)
    else:
        a_occ = torch.ones(len(ag), dtype=torch.float, device=device)  # (N,)
    a_pref, a_invs2 = get_fmod_param(ag)
    a_pref = torch.tensor(a_pref, dtype=torch.float, device=device)  # (N, 5)
    a_invs2 = torch.tensor(a_invs2, dtype=torch.float, device=device)  # (N, 5)

    # Cutoff -> box of offset indices
    oi, oj, ok = torch.ceil(cutoff / delta).to(torch.int)
    offsets = [torch.arange(-o, o + 1, device=device) for o in (oi, oj, ok)]
    offset_kernel = torch.cartesian_prod(*offsets)  # 3 * (Δ, 3)

    for _a_xyz, _a_occ, _a_pref, _a_invs2 in batched(
        a_xyz, a_occ, a_pref, a_invs2, batch_size=chunk_size
    ):
        a_ijk = torch.floor((_a_xyz - origin) / delta).to(torch.int)  # (N, 3)

        av_ijk = a_ijk[:, None] + offset_kernel[None, :]  # (N, Δ, 3)
        del a_ijk

        av_xyz = origin + av_ijk * delta  # (N, Δ, 3)

        av_d2 = torch.sum((av_xyz - _a_xyz[:, None]) ** 2, dim=-1)  # (N, Δ)
        del av_xyz

        av_dens = (
            (av_d2 <= cutoff**2)
            * torch.sum(
                _a_pref[:, None, :]
                * torch.exp(
                    -0.5 * _a_invs2[:, None, :] * av_d2[:, :, None]
                ),  # (N, Δ, 5)
                dim=-1,
            )
            * _a_occ[:, None]
        )  # (N, Δ)
        del av_d2

        av_i, av_j, av_k = av_ijk.moveaxis(-1, 0)
        grid.index_put_((av_i, av_j, av_k), av_dens, accumulate=True)

    # metadata = dict(info=f"Density map from {ag} at {ag.ts}") <- This takes ages

    # TODO: It would be better to just return on-device but `Grid` only accepts `numpy` arrays
    if vox_indices is None:
        return Grid(
            grid=grid.cpu().numpy(),
            origin=origin.cpu().numpy(),
            delta=delta.cpu().numpy(),
        )
    else:
        data = grid.reshape(-1)[vox_indices]
        return SparseGrid(
            data=data.cpu().numpy(),
            vox_indices=vox_indices.cpu().numpy(),
            shape=grid.shape,
            origin=origin.cpu().numpy(),
            delta=delta.cpu().numpy(),
        )


def file_to_map(
    file,
    shape: tuple[int],
    delta: tuple[float] | float,
    origin: tuple[float] = (0, 0, 0),
    vox_indices=None,
    sel="all",
    cutoff=4.0,
    tempfactor=None,
    device="cpu",
    verbose=False,
):
    """
    Computes density from `file`.
    Returns either a 3-dimensional map (`outshape='map'`) or a flat array
    corresponding to `vox_indices`.
    """
    if verbose:
        pid = f"[{mp.current_process().name}]" if mp.current_process().daemon else ""
        print(f"{pid} Computing density from:", file, flush=True)

    t0 = time.time()
    u = mda.Universe(file)
    ag = u.select_atoms(sel)
    if tempfactor is not None:
        ag.tempfactors = tempfactor
    atom_count = dict(Counter(ag.types))
    t1 = time.time()
    map_ = atomgroup_to_map(
        ag,
        shape=shape,
        delta=delta,
        origin=origin,
        cutoff=cutoff,
        device=device,
    )
    t2 = time.time()

    if verbose:
        print(f"{pid} Done computing map.")
        print("  Atoms:", atom_count)
        print(f"  Loading: {t1 - t0:.2f} s")
        print(f"  Map: {t2 - t1:.2f} s", flush=True)

    if vox_indices is not None:
        return map_.to_sparse(vox_indices=vox_indices)
    else:
        return map_


def main():
    parser = argparse.ArgumentParser(
        prog="python forward_density.py",
        description="Forward mapping and validation of structural models and trajectories against cryo-EM map",
    )
    parser.add_argument("input", nargs="+", type=str, help="Input structure(s)")
    parser.add_argument(
        "--sel",
        type=str,
        default="not type H",
        help="('all') MDAnalysis selection string for atoms that will be backmapped.",
    )
    parser.add_argument("--map", type=str, help="Input MRC map")
    parser.add_argument(
        "-o",
        "--average-out",
        type=str,
        default="map_out.mrc",
        help="(map_out.mrc) Output average MRC map (*.mrc)",
    )
    parser.add_argument(
        "--frame-out",
        type=str,
        help="Output sparse map for each input structure/frame as a numpy (*.npy) flat array corresponding to selected voxels.",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=4.0,
        help="(4.0) Cutoff distance from atoms for Gaussians (nm)",
    )
    parser.add_argument(
        "--device",
        type=str,
        choices=["auto", "cpu", "gpu"],
        default="auto",
        help="(auto) Device",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        help="(None) Size of chunks (number of atoms) to divide computation. If None will try to guess from available memory. Try a smaller number if you get CUDA out-of-memory errors.",
    )
    parser.add_argument(
        "--ref", type=str, help="Reference structure used for voxel selection"
    )
    parser.add_argument("--vox-index", type=str, help="Input voxel index NumPy file")
    parser.add_argument(
        "--vox-index-out",
        type=str,
        default="vox_indices.npy",
        help="Output voxel index NumPy file",
    )
    parser.add_argument(
        "--vox-around",
        nargs=2,
        metavar=("R", "SEL"),
        help="Select voxels around R the SEL in reference structure",
        action="append",
    )
    parser.add_argument(
        "--vox-not-around",
        nargs=2,
        metavar=("R", "SEL"),
        help="Exclude voxels around R the SEL in reference structure",
        action="append",
    )
    parser.add_argument(
        "--tempfactor", type=float, help="Enforce this temperature factor (B-factor)"
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=0,
        help="Number of processes for multiprocessing.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="More info")
    parser.add_argument("--progress", action="store_true", help="Display progress bar")

    args = parser.parse_args()

    if args.device == "auto":
        if torch.cuda.is_available():
            device = "cuda"
        else:
            device = "cpu"
    elif args.device == "gpu":
        device = "cuda"
    else:
        device = "cpu"
    print(f"Device: {device}")

    if args.map is not None:
        refmap = Grid(args.map)
        map_params = {
            "x0": refmap.origin,
            "dx": refmap.delta,
            "shape": refmap.shape,
            "size": refmap.size,
        }
        print("Input Cryo-EM map:")
        print("  Filename:", args.map)
    else:
        raise NotImplementedError("Reference map is needed")

    print("  Parameters:")
    print("    x0:", map_params["origin"])
    print("    dx:", map_params["delta"])
    print("    shape:", map_params["shape"])
    print("    size:", map_params["size"])

    if args.vox_indices is None and args.vox_around is None:
        vox_indices = np.arange(map_params["size"])
    else:
        vox_indices = np.array([], dtype=int)

    if args.vox_indices is not None:
        vox_indices = np.load(args.vox_indices)
        print(f"Selecting {len(vox_indices)} from voxel index file:", args.vox_indices)

    if args.vox_around is not None:
        if args.ref is None:
            parser.error("`--vox-around` requires setting `--ref`")
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
            parser.error("`--vox-not-around` requires setting `--ref`")
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

    print("Final voxel number:", len(vox_indices))
    np.save(args.vox_indices_out, vox_indices)
    print("Saved voxel index file:", args.vox_indices_out)

    if args.chunk_size is None:
        mem = (
            torch.cuda.memory.mem_get_info()[0]
            if device == "cuda"
            else psutil.virtual_memory().available
        )
        if args.verbose:
            print(f"Detected total of free memory: {mem * 1e-9:.2f} Gb")
        # Heuristic
        kernel_size = np.prod(
            2 * np.ceil(args.cutoff / np.asarray(map_params["delta"])) + 1
        )
        chunk_size = 2 ** int(
            np.log2(mem)
            - np.log2(kernel_size)
            - np.log2(torch.float.itemsize)
            - np.log2(max(1, args.processes))
            - np.log2(18)
        )
    else:
        chunk_size = args.chunk_size
    print(f"Using maximum chunk size of {chunk_size}")

    if len(vox_indices) == map_params["size"]:
        vox_indices = None

    _get_model_map = partial(
        file_to_map,
        sel=args.sel,
        tempfactor=args.tempfactor,
        vox_indices=vox_indices,
        device=device,
        cutoff=args.cutoff,
        chunk_size=chunk_size,
        verbose=args.verbose,
        **map_params,
    )

    def run(iter):
        tot_map = None
        if args.progress:
            iter = tqdm(iter)
        for i, map_ in enumerate(iter):
            if args.frame_out is not None:
                outfile = Path(args.frame_out)
                outfile = outfile.with_stem(f"{outfile.stem}_{i}").with_suffix(".pkl")
                map_.save(outfile)
                print(f"Saving frame {i}: {outfile}")
            tot_map = map_ if tot_map is None else tot_map + map_
            if args.verbose:
                print(f"Accumulating model map {i}/{len(args.input)}", flush=True)
        return tot_map

    if args.processes > 0:
        mp.set_start_method("spawn")
        with mp.Pool(args.processes) as pool:
            iter = pool.imap(_get_model_map, args.input)
            tot_map = run(iter)
    else:
        iter = args.input
        tot_map = run(iter)

    print("Averaging", len(args.input), "models.")
    tot_map /= len(args.input)

    tot_map.to_dense().export(args.average_out, file_format="mrc")
    print("Output file:", args.average_out)


if __name__ == "__main__":
    main()
