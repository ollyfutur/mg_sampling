"""Loaders of analysis outputs of the simulated systems."""

import pandas as pd

from mg_sampling.io import read_plumed
from mg_sampling.motifs import compress_motifs
from pipeline.systems import SYSTEMS, is_replex, replica_dirs


def load_occupancy(system, k, ion="MG", lambda_index=0):
    """
    `compress_motifs` of each continuous trajectory of replicate `k` of `system` (`SYSTEMS`),
    restricted to the production window and, for replex, to the frames at `lambda_index`.
    """
    window = SYSTEMS[system]
    blocks = []
    for i, path in enumerate(replica_dirs(system, k)):
        df = pd.read_parquet(path / f"{ion}_motifs.parquet")
        time_ns = df["time"].to_numpy() / 1000
        keep = (time_ns >= window.start) & (time_ns < window.end)
        if is_replex(system):
            colvar = read_plumed(path / f"bias.COLVAR.{i}")
            state = colvar.set_index(colvar["time"].round().astype(int))["bias.index"]
            keep &= (
                state.reindex(df["time"].round().astype(int)).to_numpy() == lambda_index
            )
        blocks.append(compress_motifs(df[keep]))
    return blocks
