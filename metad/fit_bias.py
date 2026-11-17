"""Fit the barrier-flattening bias parameters (E, sX, sY) on a metadynamics run."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from pipeline.bias import ASSOCIATIVE_BOUNDS, ASSOCIATIVE_P_0, PARAMS_NAMES, bias_fn
from pipeline.constants import RT
from mg_sampling.io import read_plumed, read_plumed_fes, to_wsv
from mg_sampling.stats import edges_to_range


class MetadSim:
    """CVs and metadynamics weights from `static.COLVAR` (CV names from `fes.dat`)."""

    def __init__(self, path):
        colvar = read_plumed(path / "static.COLVAR")
        (x_name, y_name), _, _ = read_plumed_fes(path / "fes.dat")
        self.x = colvar[x_name].to_numpy()
        self.y = colvar[y_name].to_numpy()
        self.X = self.x - self.y + 6
        self.weights = np.exp(colvar["metad.bias"].to_numpy() / RT)


def associative_mask(x, y):
    X = x - y + 6
    Y = x + y - 6
    return (X > 0) * (X < 2) * (Y > 0) * (Y < 2)


def associative_loss_fn(sim: MetadSim, bias_fn, params):
    """Barrier height along X in the associative square, after removing the trend between the two minima."""
    X_edges = np.linspace(0, 2, 2 * 20 + 1)
    X_range = edges_to_range(X_edges)
    mask = associative_mask(sim.x, sim.y)
    bias_weights = np.exp(-bias_fn(sim.x, sim.y, *params) / RT) * mask
    ρX, _ = np.histogram(sim.X, bins=X_edges, density=True, weights=sim.weights * bias_weights)
    FX = -RT * np.log(ρX)
    imid = len(FX) // 2
    i1, i2 = np.argmin(FX[:imid]), np.argmin(FX[imid:]) + imid
    FX -= X_range * (FX[i2] - FX[i1]) / (X_range[i2] - X_range[i1])
    return np.max(FX) - np.min(FX)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="Metadynamics folder (with static.COLVAR and fes.dat)")
    parser.add_argument("-o", "--output", type=Path, help="Output file (default: <path>/bias_params.wsv)")
    args = parser.parse_args()

    output = args.output or args.path / "bias_params.wsv"
    if output.exists():
        raise SystemExit(f"{output} exists, not overwriting.")

    sim = MetadSim(args.path)
    with np.errstate(divide="ignore", invalid="ignore"):
        res = minimize(
            lambda params: associative_loss_fn(sim, bias_fn, params),
            ASSOCIATIVE_P_0,
            method="Nelder-Mead",
            bounds=ASSOCIATIVE_BOUNDS,
        )
    print(res)
    to_wsv(pd.DataFrame([res.x], columns=PARAMS_NAMES), output)


if __name__ == "__main__":
    main()
