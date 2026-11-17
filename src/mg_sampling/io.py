"""Readers and writers for PLUMED outputs and whitespace-separated tables."""

import numpy as np
import pandas as pd


def to_wsv(df: pd.DataFrame, path):
    with open(path, "w") as f:
        f.write(df.to_string(index=False))


def read_wsv(path):
    df = pd.read_csv(path, delimiter=r"\s+")
    return df


def read_plumed(path, drop_duplicates=True):
    with open(path, "r") as f:
        header = next(f)
        columns = header.split()[2:]
    df = pd.read_csv(path, delimiter=r"\s+", comment="#", header=None)
    df.columns = columns
    if drop_duplicates:
        df = df.drop_duplicates(subset="time", keep="last")
    return df


def read_plumed_fes(path):
    """
    Reads output of `sum_hills` (fes.dat).
    Returns
    -------
    cv_names : tuple of str
    cv_ranges : tuple of 1d ndarrays
    fes : ndarray
    """
    df = read_plumed(path, drop_duplicates=False)
    d = (df.shape[1] - 1) // 2
    df = df.set_index(list(df.columns[:d])).sort_index()
    fes_shape = tuple([len(x) for x in df.index.levels])
    fes = np.reshape(df["file.free"], fes_shape)
    cv_names = tuple([idx.name for idx in df.index.levels])
    cv_ranges = tuple([np.array(idx) for idx in df.index.levels])
    return cv_names, cv_ranges, fes
