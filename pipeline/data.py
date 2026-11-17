"""Loaders of analysis outputs of the simulated systems."""

import numpy as np
import pandas as pd

from mg_sampling.io import read_plumed
from mg_sampling.motifs import compress_motifs
from pipeline.systems import in_window, is_replex, replica_dirs, replicas


def bias_state(path, i, time_ps):
    """λ index of replex replica `i` (folder `path`) at times `time_ps`, from `bias.COLVAR.<i>`."""
    colvar = read_plumed(path / f"bias.COLVAR.{i}")
    state = colvar.set_index(colvar["time"].round().astype(int))["bias.index"]
    return state.reindex(np.round(time_ps).astype(int)).to_numpy().astype(int)


def load_occupancy(system, k, ion="MG", lambda_index=0):
    """
    `compress_motifs` of each continuous trajectory of replicate `k` of `system` (`SYSTEMS`),
    restricted to the production window and, for replex, to the frames at `lambda_index`.
    """
    blocks = []
    for i, path in enumerate(replica_dirs(system, k)):
        df = pd.read_parquet(path / f"{ion}_motifs.parquet")
        keep = in_window(system, df["time"].to_numpy() / 1000)
        if is_replex(system):
            keep &= bias_state(path, i, df["time"]) == lambda_index
        blocks.append(compress_motifs(df[keep]))
    return blocks


def _load(system, filename, columns=None, stride_ps=None):
    """Concatenate a per-replica table of `system`, with `replicate`, `replica`, `state` (λ index)."""
    tables = []
    for k, i, path in replicas(system):
        if not (path / filename).exists():
            continue
        df = pd.read_parquet(path / filename, columns=columns)
        if stride_ps is not None:
            df = df[np.round(df["time"]).astype(int) % stride_ps == 0]
        state = bias_state(path, i, df["time"]) if is_replex(system) else 0
        tables.append(df.assign(replicate=k, replica=i, state=state))
    return pd.concat(tables, ignore_index=True) if tables else None


def load_motifs(system, ion="MG", columns=None, stride_ps=1000):
    """`{ion}_motifs.parquet` of all replicas of `system` every `stride_ps` (None if absent); time in ns."""
    df = _load(system, f"{ion}_motifs.parquet", columns, stride_ps)
    if df is not None:
        df["time"] = df["time"] / 1000
    return df


def load_ion_counts(system, filename="ion_counts.parquet"):
    """Ion counts (`analysis/ion_counts.py`) of all replicas of `system`; time in ns."""
    df = _load(system, filename)
    if df is not None:
        df["time"] = df["time"] / 1000
    return df
