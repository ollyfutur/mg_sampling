"""
Extension of `MDAnalysis.GridDataFormats` for manipulating maps.
"""

from math import prod
import operator
from numbers import Number
import pickle

import gridData

import numpy as np

import mrcfile


class Grid(gridData.Grid):
    """
    Overload `gridData.Grid` to add functionalities.
    """

    def __init__(
        self,
        grid=None,
        edges=None,
        origin=None,
        delta=None,
        metadata=None,
        interpolation_spline_order=3,
        file_format=None,
    ):
        super().__init__(
            grid,
            edges,
            origin,
            delta,
            metadata,
            interpolation_spline_order,
            file_format,
        )

        self._exporters.update({"MRC": self._export_mrc})

    def _export_mrc(self, filename, overwrite=True, **kwargs):
        with mrcfile.new(filename, overwrite=overwrite) as mrc:
            mrc.set_data(
                self.grid.transpose((2, 1, 0)),  # Default C-style ordering Z, Y, X
            )
            mrc.voxel_size = tuple(self.delta)
            mrc.header.origin = tuple(self.origin)

    def export(self, filename, file_format=None, type=None, typequote='"', **kwargs):
        filename = str(filename)
        exporter = self._get_exporter(filename, file_format=file_format)
        exporter(filename, type=type, typequote=typequote, **kwargs)

    def translate(self, vector, inplace=False):
        """
        Translates the data with respect to the frame by a given `vector` in coordinate space, interpolating values.
        To translate the frame together with the data (simple change of coordinates) use `translate_origin`
        """
        midpoints = (x - offset for x, offset in zip(self.midpoints, vector))
        coordinates = np.meshgrid(*midpoints, indexing="ij")
        newgrid = self.interpolated(*coordinates).astype(self.grid.dtype)
        if inplace:
            self.grid = newgrid
        else:
            return self.__class__(newgrid, self.edges)

    def translate_origin(self, vector, inplace=False):
        newedges = [x + offset for x, offset in zip(self.edges, vector)]
        if inplace:
            self.edges = newedges
            self._update()
        else:
            return self.__class__(self.grid, newedges)

    def rotate(self, R, point=(0, 0, 0)):
        """
        Rotates the data about pivot `point` by applying rotation matrix `R` and interpolating values.
        """
        R = np.asarray(R)
        point = np.asarray(point)
        coordinates = np.meshgrid(*self.midpoints, indexing="ij")
        coordinates = np.stack([np.ravel(xx) for xx in coordinates], axis=-1)
        coordinates -= point
        coordinates = np.dot(coordinates, R)
        coordinates += point
        coordinates = np.unstack(coordinates, axis=-1)
        coordinates = [np.reshape(xx, self.shape) for xx in coordinates]
        newgrid = self.interpolated(*coordinates).astype(self.grid.dtype)
        return self.__class__(newgrid, self.edges)

    def apply_tempfactor(self, B, inplace=False):
        """
        Convolve the map with a gaussian corresponding to the given B-factor value (Å²).
        Negative values result in sharpening.
        """
        fft_grid = np.fft.fftn(self.grid)
        freqs = [
            np.fft.fftfreq(n, d=dx) for n, dx in zip(self.shape, self.delta)
        ]  # Å⁻¹
        kx, ky, kz = np.meshgrid(*freqs, indexing="ij")
        k2 = kx**2 + ky**2 + kz**2  # Å⁻²
        kernel = np.exp(-B * k2 / 4)
        fft_grid = fft_grid * kernel
        newgrid = np.real(np.fft.ifftn(fft_grid)).astype(self.grid.dtype)
        if inplace:
            self.grid = newgrid
        else:
            return self.__class__(newgrid, self.edges)

    @property
    def shape(self):
        return self.grid.shape

    @property
    def size(self):
        return self.grid.size

    @property
    def grid_params(self):
        return dict(shape=self.shape, edges=self.edges, origin=self.origin, delta=self.delta)

    def to_sparse(self, *, vox_indices=None, mask=None):
        if vox_indices is None and mask is None:
            raise ValueError("Must provide either `vox_indices` or a boolean `mask`")
        elif vox_indices is None:
            vox_indices = np.flatnonzero(mask)
        data = self.grid.reshape(-1)[vox_indices]
        return SparseGrid(data, vox_indices, self.edges, metadata=self.metadata)

    def to_dense(self):
        return self


class SparseGrid(object):
    """
    A lightweight representation of a grid object only storing
    values for a set of voxel indices.
    """

    def __init__(
        self,
        data,
        vox_indices,
        edges=None,
        shape=None,
        origin=None,
        delta=None,
        metadata=None,
    ):
        if edges is None:
            if shape is None or origin is None or delta is None:
                raise ValueError(
                    "To set up SparseGrid either `edges` or `shape`/`origin`/`delta` must be provided."
                )
            edges = [
                x0 + (np.arange(n + 1) - 0.5) * dx
                for x0, n, dx in zip(origin, shape, delta)
            ]
        self.data = np.asarray(data)
        self.vox_indices = np.asarray(vox_indices)
        self.edges = edges
        self.metadata = metadata if metadata is not None else {}

    @property
    def shape(self):
        return tuple(len(edges) - 1 for edges in self.edges)

    @property
    def size(self):
        return prod(self.shape)

    def to_dense(self, fill_value=0):
        grid = np.full(self.size, fill_value, dtype=self.data.dtype)
        grid[self.vox_indices] = self.data
        grid = grid.reshape(self.shape)
        return Grid(grid=grid, edges=self.edges, metadata=self.metadata)

    def check_compatible(self, other):
        """
        Warning: to avoid comparing all indices only a simple sanity check
        is performed (number of voxels and first and last indices).
        """
        if isinstance(other, self.__class__):
            is_compatible = (
                all(
                    np.allclose(other_edge, self_edge)
                    for other_edge, self_edge in zip(other.edges, self.edges)
                )
                and len(other.vox_indices) == len(self.vox_indices)
                and other.vox_indices[0] == self.vox_indices[0]
                and other.vox_indices[-1] == self.vox_indices[-1]
            )
        else:
            is_compatible = False
        if not is_compatible:
            raise TypeError(
                "Can only perform binary operations between `SparseGrid`s "
                "with identical `edges` and `vox_indices`."
            )
        return True

    def _operate(self, other, op):
        if isinstance(other, self.__class__):
            self.check_compatible(other)
            return self.__class__(
                op(self.data, other.data), self.vox_indices, self.edges
            )
        elif isinstance(other, Number):
            return self.__class__(op(self.data, other), self.vox_indices, self.edges)
        elif isinstance(other, Grid):
            return op(self.to_dense(), other)
        else:
            raise TypeError("Cannot perform binary operation between")

    def save(self, filename):
        data = dict(
            data=self.data,
            vox_indices=self.vox_indices,
            edges=self.edges,
            metadata=self.metadata,
        )
        with open(filename, "wb") as f:
            pickle.dump(data, f)

    @classmethod
    def load(cls, filename):
        with open(filename, "rb") as f:
            data = pickle.load(f)
        data["metadata"].update(filename=filename)
        return cls(**data)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (
                np.allclose(self.data, other.data)
                and np.array_equal(self.vox_indices, other.vox_indices)
                and all(map(np.allclose, self.edges, other.edges))
            )


def _get_binary_op_method(opname):
    op = getattr(operator, opname)

    def method(self, other):
        return self._operate(other, op)

    return method


for opname in ["add", "sub", "mul", "truediv", "floordiv", "pow"]:
    method = _get_binary_op_method(opname)
    setattr(SparseGrid, f"__{opname}__", method)
    setattr(SparseGrid, f"__r{opname}__", method)
