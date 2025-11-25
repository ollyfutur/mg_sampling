import gridData

import numpy as np

import mrcfile

class Grid(gridData.Grid):
    """
    Overload `gridData.Grid` to add functionalities.
    """

    def __init__(self, grid=None, edges=None, origin=None, delta=None,
                 metadata=None, interpolation_spline_order=3,
                 file_format=None):
        super().__init__(grid, edges, origin, delta, metadata, interpolation_spline_order, file_format)

        self._exporters.update({
            'MRC': self._export_mrc
        })
    
    def _export_mrc(self, filename, overwrite=False, **kwargs):
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
        coordinates = np.meshgrid(*midpoints, indexing='ij')
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
        coordinates = np.meshgrid(*self.midpoints, indexing='ij')
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
        freqs = [np.fft.fftfreq(n, d=dx) for n, dx in zip(self.shape, self.delta)]  # Å⁻¹
        kx, ky, kz = np.meshgrid(*freqs, indexing='ij')
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