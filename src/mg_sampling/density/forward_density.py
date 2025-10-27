"""
Backmapping of density map from (an ensemble of) structures
"""

import argparse
from collections import Counter

import numpy as np

import mrcfile

import MDAnalysis as mda
from MDAnalysis import AtomGroup

import torch

from tqdm import trange, tqdm

from mg_sampling.density.utils import get_map_params, save_map, get_voxels, voxels_around

π = np.pi

## you need the following libraries:
##  conda install -c conda-forge mrcfile mdanalysis biopython
##  conda install -c anaconda scipy
##
## Adapted from code by Max Bonomi (mbonomi@pasteur.fr)
##
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
WEIGHT = {
    "C": [0.0893, 0.2563, 0.7570, 1.0487, 0.3575],
    "O": [0.0974, 0.2921, 0.6910, 0.6990, 0.2039],
    "N": [0.1022, 0.3219, 0.7982, 0.8197, 0.1715],
    "S": [0.2497, 0.5628, 1.3899, 2.1865, 0.7715],
    "P": [0.2548, 0.6106, 1.4541, 2.3204, 0.8477],
    "F": [0.1083, 0.3175, 0.6487, 0.5846, 0.1421],
    "MG": [0.2314, 0.6866, 0.9677, 2.1882, 1.1339],
    "K": [0.4115, 1.4031, 2.2784, 2.6742, 2.2162],
    "NA": [0.2142, 0.6853, 0.7692, 1.6589, 1.4482],
    "CL": [0.2443, 0.5397, 1.3919, 2.0197, 0.6621],
}
SIGMA = {
    "C": [0.2465, 1.7100, 6.4094, 18.6113, 50.2523],
    "O": [0.2067, 1.3815, 4.6943, 12.7105, 32.4726],
    "N": [0.2451, 1.7481, 6.1925, 17.3894, 48.1431],
    "S": [0.2681, 1.6711, 7.0267, 19.5377, 50.3888],
    "P": [0.2908, 1.8740, 8.5176, 24.3434, 63.2996],
    "F": [0.2057, 1.3439, 4.2788, 11.3932, 28.7881],
    "MG": [0.3278, 2.2720, 10.9241, 39.2898, 101.9748],
    "K": [0.3703, 3.3874, 13.1029, 68.9592, 194.4329],
    "NA": [0.3334, 2.3446, 10.0830, 48.3037, 138.2700],
    "CL": [0.2468, 1.5242, 6.1537, 16.6687, 42.3086],
}


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
    pref = [
        WEIGHT[a.type] * (π / (SIGMA[a.type] + a.tempfactor / 4)) ** 1.5 for a in ag
    ]
    invs2 = [2 * π**2 / (SIGMA[a.type] + a.tempfactor / 4) for a in ag]
    return np.array(pref), np.array(invs2)


# get vdw radii
def get_vdw_radii(ag: AtomGroup):
    vdw = [VDW_RADII[a.type] for a in ag]
    return np.array(vdw)


def get_model_map(
    ag: AtomGroup,
    device: str = "cpu",
    vox_indices=None,
    cutoff: float = 4.0,
    chunk_size: int = 100,
    **map_params,
):
    """
    Returns (flattened) map from the atoms in `ag`.
    """

    # convert to tensor and copy to device [nat, 3]
    pos_g = torch.tensor(ag.positions).to(torch.float).to(device)
    occ_g = torch.tensor(ag.occupancies).to(torch.float).to(device)
    vox_g = torch.tensor(get_voxels(map_params, vox_indices)).to(torch.float).to(device)

    pref, invs2 = get_fmod_param(ag)
    pref, invs2 = (
        torch.tensor(pref).to(torch.float).to(device),
        torch.tensor(invs2).to(torch.float).to(device),
    )

    nvox = len(vox_g)
    niter = max(1, int(np.ceil(nvox / chunk_size)))

    for i in trange(niter):
        i0 = i * chunk_size
        i1 = i0 + chunk_size if i != niter - 1 else nvox
        if i1 - i0 <= 1:
            continue
        # neighbour list index atoms, index voxels, temporary
        nl_ia_t, nl_iv_t = torch.where(torch.cdist(pos_g, vox_g[i0:i1]) < cutoff)
        if i == 0:
            pos_nl = torch.index_select(pos_g, 0, nl_ia_t)
            vox_nl = torch.index_select(vox_g, 0, nl_iv_t)
            occ_nl = torch.index_select(occ_g, 0, nl_ia_t)
            nl_ia = nl_ia_t
            nl_iv = nl_iv_t
        else:
            temp_pos = torch.index_select(pos_g, 0, nl_ia_t)
            temp_vox = torch.index_select(vox_g, 0, nl_iv_t + i0)
            temp_occ = torch.index_select(occ_g, 0, nl_ia_t)

            pos_nl = torch.cat((pos_nl, temp_pos), 0)
            vox_nl = torch.cat((vox_nl, temp_vox), 0)
            occ_nl = torch.cat((occ_nl, temp_occ), 0)
            nl_ia = torch.cat((nl_ia, nl_ia_t), 0)
            nl_iv = torch.cat((nl_iv, nl_iv_t + i0), 0)
    # calculate distance squared positions / voxels in neighbor list [nvox_nl]
    dist2 = torch.sum(torch.pow(pos_nl - vox_nl, 2), 1)
    # expand distances cutoff [nvox_nl, 5]
    dist2 = dist2[:, None].expand(-1, 5)
    # calculate map [nvox_nl]
    map_g = (
        torch.sum(
            torch.index_select(pref, 0, nl_ia)
            * torch.exp(-0.5 * torch.index_select(invs2, 0, nl_ia) * dist2),
            1,
        )
        * occ_nl
    )
    # initialize sum map on device [nvox]
    map_sum = torch.zeros(nvox, dtype=torch.float, device=device)
    # sum entries of map_g [nvox]
    map_sum.index_add_(0, nl_iv, map_g)
    # and return map [on GPU]
    return map_sum

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
        help="('all') MDAnalysis selection string",
    )
    parser.add_argument("--map", type=str, help="Input MRC map")
    parser.add_argument("-o", "--out", type=str, help="(map_out.mrc) Output MRC map")
    # parser.add_argument(
    # "--threshold", type=float, metavar="TH", default=-10000.0, help="map thresold"
    # )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=4.0,
        help="(4.0) Cutoff distance for Gaussians (nm)",
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
        default=1000,
        help="(1000) Size of chunks to divide voxel computation",
    )
    parser.add_argument(
        "--ref", type=str, help="Reference structure used for voxel selection"
    )
    parser.add_argument("--vox-indices", type=str, help="Voxel indices file name")
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
        "--tempfactor", type=float, help="Force this temperature factor (B-factor)"
    )

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
        with mrcfile.open(args.map, permissive=True) as mrc:
            map_params = get_map_params(mrc)
        print("Input Cryo-EM map:")
        print("  Filename:", args.map)
        print("  Parameters:")
        print("    x0:", map_params["x0"])
        print("    dx:", map_params["dx"])
        print("    nbin:", map_params["nbin"])
        print("    nbin_tot:", map_params["nbin_tot"])

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

    assert vox_indices.dtype == int
    print("Final voxel number:", len(vox_indices))

    tot_map = None
    for topfile in tqdm(args.input):
        u = mda.Universe(topfile)
        ag = u.select_atoms(args.sel)
        if args.tempfactor is not None:
            ag.tempfactors = args.tempfactor
        model_map = get_model_map(
            ag,
            vox_indices=vox_indices,
            device=device,
            cutoff=args.cutoff,
            chunk_size=args.chunk_size,
            **map_params,
        )
        print("Computing density from:", topfile)
        print("  Atoms:", dict(Counter(ag.types)))
        if tot_map is None:
            tot_map = model_map
        else:
            tot_map += model_map
    print("Averaging", len(args.input), "models.")
    tot_map /= len(args.input)

    if args.out is not None:
        _map = np.zeros(map_params["nbin_tot"], dtype=np.float32)
        _map[vox_indices] = tot_map.cpu().numpy()
        _map = _map.reshape(map_params["nbin"])
        save_map(args.out, _map, map_params)
        print("Output file:", args.out)


if __name__ == "__main__":
    main()
