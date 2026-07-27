import numpy as np

import matplotlib.pyplot as plt

from mg_sampling.density.grid import Grid


def plot_map_slice(
    grid: Grid,
    func: str = "imshow",
    plane="xy",
    depth=0.5,
    depth_mode="fraction",
    thickness=0,
    ax=None,
    **kwargs,
):
    if ax is None:
        ax = plt.gca()

    if isinstance(func, str):
        func = getattr(ax, func)

    x_axis, y_axis = ["xyz".index(dim) for dim in plane]
    depth_axis = ({0, 1, 2} - set((x_axis, y_axis))).pop()
    midpoints = grid.midpoints[depth_axis]
    dx = grid.delta[depth_axis]
    x_edges, y_edges = grid.edges[x_axis], grid.edges[y_axis]
    extent = [*x_edges[[0, -1]], *y_edges[[0, -1]]]
    z_edges = grid.edges[depth_axis]

    if depth_mode == "fraction":
        k = int(depth * len(midpoints))
        dk = int(thickness * len(midpoints))
    elif depth_mode == "index":
        k = int(depth)
        dk = int(thickness)
    elif depth_mode == "coord":
        k = int(np.searchsorted(midpoints, depth, side="left"))
        dk = int(np.searchsorted(midpoints - midpoints[0], thickness, side="left"))

    if thickness == 0:
        s = tuple(slice(None) if d in (x_axis, y_axis) else k for d in (0, 1, 2))
        img = grid.grid[s].T
    else:
        s = tuple(slice(None) if d in (x_axis, y_axis) else slice(k-dk, k+dk) for d in (0, 1, 2))
        img = grid.grid[s].mean(axis=depth_axis).T

    out = func(img, origin="lower", extent=extent, **kwargs)

    ax.set_xlabel("xyz"[x_axis])
    ax.set_ylabel("xyz"[y_axis])

    return out