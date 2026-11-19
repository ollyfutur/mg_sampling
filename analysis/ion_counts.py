"""
Total and bulk counts of water oxygens and ions along each trajectory of one replicate of a
system, for several bulk cutoffs. Bulk: farther than `radius` from every RNA phosphorus.
Written to `<replica folder>/ion_counts.parquet` (one row per time, name and radius; time in ps).

Examples:
    python ion_counts.py plain_md/7ez0/mMg_tip3p 3
    python ion_counts.py plain_md/7ez0/mMg_tip3p 3 --radii 20 30 40 50 60
    python ion_counts.py replex/7ez0/mMg_tip3p 1 --stride 100   # quick check

Existing outputs are skipped (use --overwrite to recompute).
"""

import argparse
import time

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from pipeline.systems import SYSTEMS, replica_dirs

RADII = [30, 40, 50, 60, 70]  # Å
NAMES = ["OW", "MG", "NA", "CL"]
OUTPUT = "ion_counts.parquet"


def count(path, stride, radii):
    u = mda.Universe(path / "input.gro", path / "run.xtc")
    phosphorus = u.select_atoms("name P")
    groups = [u.select_atoms(f"name {name}") for name in NAMES]
    atoms = sum(groups[1:], groups[0])
    bounds = np.cumsum([0] + [len(g) for g in groups])
    rows = []
    for ts in u.trajectory[::stride]:
        # Distance of each atom to the closest phosphorus (inf beyond the largest radius)
        pairs, distances = capped_distance(atoms.positions, phosphorus.positions, max(radii), box=ts.dimensions)
        d_min = np.full(len(atoms), np.inf)
        np.minimum.at(d_min, pairs[:, 0], distances)
        for name, start, stop in zip(NAMES, bounds[:-1], bounds[1:]):
            for radius in radii:
                rows.append({
                    "time": ts.time,
                    "name": name,
                    "radius": radius,  # Å
                    "n_bulk": int(np.sum(d_min[start:stop] > radius)),
                    "n_tot": int(stop - start),
                    "volume": ts.volume,  # Å³
                })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("system", choices=SYSTEMS, metavar="system", help=f"One of: {', '.join(SYSTEMS)}")
    parser.add_argument("replicate", type=int, help="Replicate number (replex: 0-1, plain MD: 0-7)")
    parser.add_argument("--radii", type=float, nargs="+", default=RADII,
                        help=f"Bulk cutoffs in Å (default: {' '.join(map(str, RADII))})")
    parser.add_argument("--stride", type=int, default=10, help="Every N trajectory frames (default: 10, i.e. 1 ns)")
    parser.add_argument("--output", default=OUTPUT, help=f"Output file name in each replica folder (default: {OUTPUT})")
    parser.add_argument("--overwrite", action="store_true", help="Recompute existing outputs")
    args = parser.parse_args()

    for path in replica_dirs(args.system, args.replicate):
        output = path / args.output
        if output.exists() and not args.overwrite:
            print(f"[{path}] exists, skipping", flush=True)
            continue
        t0 = time.time()
        df = count(path, args.stride, args.radii)
        tmp = output.with_name(output.name + ".tmp")
        df.to_parquet(tmp)
        tmp.rename(output)
        print(f"[{path}] {df['time'].nunique()} frames in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
