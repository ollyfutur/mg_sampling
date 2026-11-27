"""
Mean density map of each continuous trajectory (block) of each system, from the per-frame
maps of `density_maps.py`: one plain MD replica, or one replex replica within a replicate.
Written to `density/blocks/<family>/<sys>_<ff>` as a `GridEnsemble` with one frame per block,
weighted by its number of frames (frame info: replicate, replica, n_frames).

Examples:
    python block_maps.py                                  # all families and systems
    python block_maps.py --family local plain_md/7ez0/mMg_tip3p

Existing outputs are skipped (use --overwrite to recompute).
"""

import argparse
import shutil
import time

import numpy as np
from cryomap import GridEnsemble, read_ensemble

from pipeline.density import FAMILIES, blocks_path, family_dir, group_means
from pipeline.systems import SYSTEMS, is_replex


def compute(family, system, overwrite=False):
    output = blocks_path(family, system)
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
    means, replicates, replicas, counts = [], [], [], []
    for k, path in enumerate(inputs):
        ensemble = read_ensemble(path, mmap=True)
        labels = ensemble.frame_info["replica_index"] if is_replex(system) else np.zeros(ensemble.n_frames, int)
        groups, group_mean, group_count = group_means(ensemble, labels)
        means.append(group_mean)
        replicates += [k] * len(groups)
        replicas += list(groups)
        counts += list(group_count)
    counts = np.array(counts)
    blocks = GridEnsemble(
        np.concatenate(means).astype(np.float32),
        ensemble.params,
        log_weights=np.log(counts),
        frame_info={"replicate": np.array(replicates), "replica": np.array(replicas), "n_frames": counts},
    )

    tmp = output.with_name(output.name + ".tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    blocks.write(tmp)
    shutil.rmtree(output, ignore_errors=True)
    tmp.rename(output)
    print(f"[{family} {system}] {len(counts)} blocks, {counts.sum()} frames, in {time.time() - t0:.0f} s", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("systems", nargs="*", metavar="system",
                        help=f"Systems (default: all): {', '.join(SYSTEMS)}")
    parser.add_argument("--family", choices=FAMILIES, action="append",
                        help="Map family, repeatable (default: all)")
    parser.add_argument("--overwrite", action="store_true", help="Recompute existing outputs")
    args = parser.parse_args()
    unknown = set(args.systems) - set(SYSTEMS)
    if unknown:
        parser.error(f"unknown systems: {', '.join(sorted(unknown))}")
    for family in args.family or FAMILIES:
        for system in args.systems or SYSTEMS:
            compute(family, system, args.overwrite)


if __name__ == "__main__":
    main()
