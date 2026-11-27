"""
Mean map of the frames where each Mg²⁺ binding motif is occupied, per continuous trajectory
(block: plain MD replica, or replex replica within a replicate), from the per-frame maps of
`density_maps.py`. Only the voxels within `MOTIF_RADIUS + BLUR_MARGIN` of the motif atoms
in 9CBU are kept (see `pipeline.density`). Written to `density/motifs/<family>/<sys>_<ff>/`:
- `motifs.parquet`: the motif table, with the total numbers of occupied frames and voxels;
- `<motif>/`: a `GridEnsemble` with one frame per block where the motif is occupied,
  weighted by its number of occupied frames (frame info: replicate, replica, n_present,
  n_frames).
The mean map of the frames where a motif is not occupied follows from the block means M of
`block_maps.py`: (n_frames·M − n_present·P) / (n_frames − n_present).

Examples:
    python motif_maps.py                                  # all families, systems with Mg
    python motif_maps.py --family local plain_md/7ez0/mMg_tip3p

Existing outputs are skipped (use --overwrite to recompute).
"""

import argparse
import shutil
import time
from collections import defaultdict

import numpy as np
import pandas as pd
from cryomap import GridEnsemble, read_ensemble

from pipeline.data import load_occupancy
from pipeline.density import FAMILIES, family_dir, motif_coords, motif_maps_path, motif_voxels
from pipeline.systems import SYSTEMS

DEFAULT_SYSTEMS = [s for s in SYSTEMS if "nomg" not in s]  # 7ez0_nomg has no Mg²⁺
CHUNK_FRAMES = 64


def occupancy_table(system, k, times, replicas):
    """Motif table of replicate `k`, and occupied (bool, frames x motifs) at each map frame
    (`times` in ps, continuous trajectory `replicas`)."""
    blocks = load_occupancy(system, k)
    info = pd.concat([i for _, i in blocks])
    info = info[~info.index.duplicated()]
    present = np.zeros((len(times), len(info)), dtype=bool)
    for i, (occupancy, _) in enumerate(blocks):
        rows = np.flatnonzero(replicas == i)
        occupancy = occupancy.set_axis(np.round(occupancy.index).astype(int))
        occupancy = occupancy.reindex(columns=info.index, fill_value=0).reindex(times[rows])
        if occupancy.isna().any(axis=None):
            raise ValueError(f"{system} {k} replica {i}: map frames without motif data")
        present[rows] = occupancy.to_numpy() > 0
    return info, present


def compute(family, system, overwrite=False):
    output = motif_maps_path(family, system)
    if output.exists() and not overwrite:
        print(f"[{family} {system}] exists, skipping", flush=True)
        return
    name = system.replace("/", "_")
    inputs = [family_dir(family) / f"{name}_{k}" for k in range(SYSTEMS[system].n_replicates)]
    missing = [str(p) for p in inputs if not p.exists()]
    if missing:
        print(f"[{family} {system}] missing replicates, skipping: {missing}", flush=True)
        return

    t0 = time.time()
    params = None
    infos, voxels = [], {}  # motif: flat indices of its voxels, None if not in 9CBU
    sums = defaultdict(dict)  # motif: {(replicate, replica): sum of the occupied frames}
    n_present, n_frames = defaultdict(dict), {}
    for k, path in enumerate(inputs):
        ensemble = read_ensemble(path, mmap=True)
        params = params or ensemble.params
        if not np.array_equal(ensemble.params.voxel_indices, params.voxel_indices):
            raise ValueError(f"{path}: other voxels than replicate 0")
        stored = params.voxel_indices
        times = np.round(ensemble.frame_info["time"]).astype(int)
        replicas = ensemble.frame_info.get("replica_index", np.zeros(len(times), int))
        info, present = occupancy_table(system, k, times, replicas)
        infos.append(info)
        motifs = list(info.index)
        for motif in motifs:
            if motif not in voxels:
                coords = motif_coords(info.loc[motif])
                # Positions of its voxels among the stored ones
                voxels[motif] = None if coords is None else np.searchsorted(stored, motif_voxels(params, coords))
        for i in np.unique(replicas):
            n_frames[k, i] = int(np.sum(replicas == i))

        for start in range(0, len(times), CHUNK_FRAMES):
            sl = slice(start, start + CHUNK_FRAMES)
            chunk = np.asarray(ensemble.data[sl])
            for j in np.flatnonzero(present[sl].any(axis=0)):
                motif = motifs[j]
                if voxels[motif] is None:
                    continue
                for i in np.unique(replicas[sl][present[sl, j]]):
                    rows = present[sl, j] & (replicas[sl] == i)
                    values = chunk[rows][:, voxels[motif]].sum(axis=0, dtype=np.float64)
                    sums[motif][k, i] = sums[motif].get((k, i), 0) + values
                    n_present[motif][k, i] = n_present[motif].get((k, i), 0) + int(rows.sum())
        print(f"[{family} {system}] replicate {k}: {len(times)} frames, {len(motifs)} motifs", flush=True)

    info = pd.concat(infos)
    info = info[~info.index.duplicated()]
    skipped = [m for m in info.index if voxels[m] is None]
    info = info.loc[[m for m in info.index if m in sums]]
    info["n_present"] = [sum(n_present[m].values()) for m in info.index]
    info["n_voxels"] = [len(voxels[m]) for m in info.index]

    tmp = output.with_name(output.name + ".tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    for motif in info.index:
        if "/" in motif:
            raise ValueError(f"Motif name {motif!r} is not a valid directory name")
        keys = sorted(sums[motif])
        counts = np.array([n_present[motif][key] for key in keys])
        data = np.stack([sums[motif][key] for key in keys]) / counts[:, None]
        GridEnsemble(
            data.astype(np.float32),
            params.sparse(voxel_indices=params.voxel_indices[voxels[motif]]),
            log_weights=np.log(counts),
            frame_info={
                "replicate": np.array([k for k, _ in keys]),
                "replica": np.array([i for _, i in keys]),
                "n_present": counts,
                "n_frames": np.array([n_frames[key] for key in keys]),
            },
        ).write(tmp / motif)
    info.to_parquet(tmp / "motifs.parquet")
    shutil.rmtree(output, ignore_errors=True)
    tmp.rename(output)
    print(
        f"[{family} {system}] {len(info)} motifs written, {len(skipped)} skipped (atoms not in 9CBU), "
        f"in {time.time() - t0:.0f} s",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("systems", nargs="*", metavar="system",
                        help=f"Systems (default: all with Mg²⁺): {', '.join(DEFAULT_SYSTEMS)}")
    parser.add_argument("--family", choices=FAMILIES, action="append",
                        help="Map family, repeatable (default: all)")
    parser.add_argument("--overwrite", action="store_true", help="Recompute existing outputs")
    args = parser.parse_args()
    unknown = set(args.systems) - set(DEFAULT_SYSTEMS)
    if unknown:
        parser.error(f"unknown systems or without Mg²⁺: {', '.join(sorted(unknown))}")
    for family in args.family or FAMILIES:
        for system in args.systems or DEFAULT_SYSTEMS:
            compute(family, system, args.overwrite)


if __name__ == "__main__":
    main()
