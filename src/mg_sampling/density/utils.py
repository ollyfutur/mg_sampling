"""
Tools for density maps
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.interpolate import interp1d

import mrcfile

from MDAnalysis import AtomGroup


def get_grid_params(mrc) -> dict:
    mrc_params = {
        "map": [
            int(mrc.header.mapc),
            int(mrc.header.mapr),
            int(mrc.header.maps),
        ],  # data organization, 1 for x, 2 for y, 3 for z
        "shape": [int(mrc.header.nx), int(mrc.header.ny), int(mrc.header.nz)],
        "size": int(mrc.header.nx * mrc.header.ny * mrc.header.nz),
        "origin": [
            float(mrc.header.origin.x + float(mrc.header.nxstart) * mrc.voxel_size.x),
            float(mrc.header.origin.y + float(mrc.header.nystart) * mrc.voxel_size.y),
            float(mrc.header.origin.z + float(mrc.header.nzstart) * mrc.voxel_size.z),
        ],
        # dimension of one voxel the in x, y, z directions
        "delta": [
            float(mrc.voxel_size.x),
            float(mrc.voxel_size.y),
            float(mrc.voxel_size.z),
        ],
    }

    return mrc_params


def save_map(path, map, grid_params):
    # mrcfile expects numpy array with SRC convention
    # i.e. map[k, j, i] (kji <-> zyx)
    map = map.transpose((2, 1, 0))
    with mrcfile.new(path, overwrite=True) as mrc:
        mrc.set_data(map)
        mrc.voxel_size = grid_params["delta"]


def get_voxel_coords(grid_params, vox_indices=None):
    """
    Returns a list of voxel coordinates from `grid_params`, possibly restricted to
    the indices in `vox_indices`.
    """
    if vox_indices is None:
        vox_indices = np.arange(grid_params["size"])

    i, j, k = np.unravel_index(vox_indices, grid_params["shape"])
    x, y, z = (
        grid_params["origin"][0] + i * grid_params["delta"][0],
        grid_params["origin"][1] + j * grid_params["delta"][1],
        grid_params["origin"][2] + k * grid_params["delta"][2],
    )
    vox_xyz = np.stack([x, y, z], axis=1)
    return vox_xyz


def voxels_around(ag: AtomGroup, radius: float, grid_params):
    # TODO: take into account vdW radii?
    origin = np.asarray(grid_params["origin"])
    delta = np.asarray(grid_params["delta"])
    n = (np.asarray(radius) // delta).astype(int)

    # Fast bounding box preselection
    vox_indices = []
    for a in ag:
        ijk = ((a.position - origin) // delta).astype(int)
        index_ranges = [np.arange(low, high + 1) for low, high in zip(ijk - n, ijk + n)]
        bbox_ijk = tuple(map(np.ravel, np.meshgrid(*index_ranges, indexing="ij")))
        bbox_indices = np.ravel_multi_index(bbox_ijk, grid_params["shape"])
        # Distance check
        bbox_xyz = get_voxel_coords(grid_params, vox_indices=bbox_indices)
        mask = np.sum((bbox_xyz - a.position) ** 2, axis=-1) <= radius**2
        vox_indices.append(bbox_indices[mask])
    vox_indices = np.unique(np.concatenate(vox_indices))

    return vox_indices


def indices_to_mask(vox_indices, grid_params):
    i, j, k = np.unravel_index(vox_indices, grid_params["shape"])
    mask = np.zeros(grid_params["shape"], dtype=bool)
    mask[i, j, k] = True
    return mask


def flatmap_to_dense(flatmap, vox_indices, **grid_params):
    map = np.zeros(grid_params["size"], dtype=flatmap.dtype)
    map[vox_indices] = flatmap
    return map.reshape(grid_params["shape"])


class QuantileTransformation:
    """
    Given two PDFs
        x ~ p
        y ~ q
    Find f such that
        f(y) ~ p
    """

    def __init__(self, p, q, x_range, y_range=None):
        if y_range is None:
            y_range = x_range

        self.x_range = np.asarray(x_range)
        self.y_range = np.asarray(y_range)
        self.p = np.asarray(p) / trapezoid(p, x_range)
        self.q = np.asarray(q) / trapezoid(q, y_range)

        self.invP = interp1d(
            cumulative_trapezoid(p, x_range, initial=0),
            x_range,
            fill_value="extrapolate",
        )
        self.Q = interp1d(
            y_range,
            cumulative_trapezoid(q, y_range, initial=0),
            fill_value="extrapolate",
        )

    def __call__(self, x):
        return self.invP(self.Q(x))

