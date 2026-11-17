"""Barrier-flattening bias on the (n_OP, n_OW) coordination plane, and its fit settings."""

import numpy as np

from pipeline.constants import π

PARAMS_NAMES = ('E', 'sX', 'sY')
ASSOCIATIVE_P_0 = (1, 0.1, 0.2)
ASSOCIATIVE_BOUNDS = ((0, np.inf), (-1/π, 1/π), (0, 1/2))

def bias_fn(x, y, E, sX, sY):
    """
    Bias function based on a cosine envelope for the associative pathway.
    """
    X = x - y + 6
    Y = x + y - 6
    _X = X + sX * (1 - np.cos(X*π))
    _Y = Y - sY * (1 - np.cos(X*π))

    return -E * (_X > 0) * (_Y > 0) * (_Y < 1/2) * (1 - np.cos(_X*π)) * (1 - np.cos(_Y*π*4)) / 4
