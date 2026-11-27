"""Generic helpers."""

import numpy as np


def asunique(seq):
    """
    Ensure all elements of `seq` are equal and return the unique element.
    """
    it = iter(seq)
    x0 = next(it)
    for x in it:
        if x != x0:
            raise ValueError(f"Expected all elements to be equal ({x} != {x0})")
    return x0


def vectorize_dict_lookup(d: dict):
    """
    Return a function that accepts numpy arrays and maps
    them according to dictionary `d` of array-like values.
    All values of `d` must have the same shape.
    """
    dtype, shape = asunique(
        (np.asarray(v).dtype, np.asanyarray(v).shape) for v in d.values()
    )

    def vectorized(keys):
        keys = np.asarray(keys)
        out = np.empty((len(keys), *shape), dtype=dtype)
        mapped = np.zeros(len(keys), dtype=bool)
        for key in d:
            out[keys == key] = d[key]
            mapped[keys == key] = True
        if not np.all(mapped):
            raise KeyError(keys[~mapped])
        return out

    return vectorized
