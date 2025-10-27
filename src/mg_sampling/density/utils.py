"""
Tools for density maps
"""

import numpy as np

import mrcfile

import MDAnalysis as mda
from MDAnalysis import AtomGroup

def get_map_params(mrc) -> dict:
    mrc_params = {
        "map": [
            int(mrc.header.mapc),
            int(mrc.header.mapr),
            int(mrc.header.maps),
        ],  # data organization
        "nbin": [int(mrc.header.nx), int(mrc.header.ny), int(mrc.header.nz)],
        "nbin_tot": int(mrc.header.nx * mrc.header.ny * mrc.header.nz),
        "x0": [
            float(mrc.header.origin.x + float(mrc.header.nxstart) * mrc.voxel_size.x),
            float(mrc.header.origin.y + float(mrc.header.nystart) * mrc.voxel_size.y),
            float(mrc.header.origin.z + float(mrc.header.nzstart) * mrc.voxel_size.z),
        ],
        # dimension of one voxel the in x, y, z directions
        "dx": [
            float(mrc.voxel_size.x),
            float(mrc.voxel_size.y),
            float(mrc.voxel_size.z),
        ],
    }
    # internally, reorder so that is always xyz format
    axes_order = np.argsort(mrc_params["map"])
    mrc_params["filemap"] = mrc_params["map"]
    mrc_params["map"] = [mrc_params["map"][k] for k in axes_order]  # should be [0, 1, 2]
    mrc_params["nbin"] = [mrc_params["nbin"][k] for k in axes_order]
    mrc_params["x0"] = [mrc_params["x0"][k] for k in axes_order]
    mrc_params["dx"] = [mrc_params["dx"][k] for k in axes_order]
    # return dictionary
    return mrc_params


def save_map(path, map, map_params):
    with mrcfile.new(path, overwrite=True) as mrc:
        mrc.set_data(map)
        mrc.voxel_size = map_params["dx"]


def get_voxels(map_params, vox_indices=None):
    """
    Returns a list of voxel coordinates from `map_params`, possibly restricted to
    the indices in `vox_indices`.
    """
    if vox_indices is None:
        vox_indices = np.arange(map_params["nbin_tot"])

    i, j, k = np.unravel_index(vox_indices, map_params["nbin"], order="F")
    x, y, z = (
        map_params["x0"][0] + i * map_params["dx"][0],
        map_params["x0"][1] + j * map_params["dx"][1],
        map_params["x0"][2] + k * map_params["dx"][2],
    )
    vox_xyz = np.stack([x, y, z], axis=1)
    return vox_xyz


def voxels_around(ag: AtomGroup, radius: float, map_params):
    # TODO: take into account vdW radii?
    x0 = np.asarray(map_params["x0"])
    dx = np.asarray(map_params["dx"])
    n = (np.asarray(radius) // dx).astype(int)

    # Fast bounding box preselection
    vox_indices = []
    for a in ag:
        ijk = ((a.position - x0) // dx).astype(int)
        index_ranges = [np.arange(low, high + 1) for low, high in zip(ijk - n, ijk + n)]
        bbox_ijk = tuple(map(np.ravel, np.meshgrid(*index_ranges, indexing="ij")))
        bbox_indices = np.ravel_multi_index(bbox_ijk, map_params["nbin"], order="F")
        # Distance check
        bbox_xyz = get_voxels(map_params, vox_indices=bbox_indices)
        mask = np.sum((bbox_xyz - a.position) ** 2, axis=-1) <= radius**2
        vox_indices.append(bbox_indices[mask])
    vox_indices = np.unique(np.concatenate(vox_indices))

    return vox_indices