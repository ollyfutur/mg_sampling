"""
Total and bulk counts of water oxygens and ions along each trajectory of one replicate of a
system. Bulk: farther than 40 Å from every RNA phosphorus. Written to
`<replica folder>/ion_counts.parquet` (time in ps).

Examples:
    python ion_counts.py plain_md/7ez0/mMg_tip3p 3
    python ion_counts.py replex/7ez0/mMg_tip3p 1 --stride 100   # quick check

Existing outputs are skipped (use --overwrite to recompute).
"""

import argparse
import time

import MDAnalysis as mda
import pandas as pd

from pipeline.systems import SYSTEMS, replica_dirs

RADIUS = 40  # Å
NAMES = ["OW", "MG", "NA", "CL"]
OUTPUT = "ion_counts.parquet"


def count(path, stride):
    u = mda.Universe(path / "input.gro", path / "run.xtc")
    total = {name: u.select_atoms(f"name {name}") for name in NAMES}
    bulk = {
        name: total[name].select_atoms(f"not around {RADIUS} global name P", updating=True)
        for name in NAMES
    }
    rows = []
    for ts in u.trajectory[::stride]:
        for name in NAMES:
            rows.append({
                "time": ts.time,
                "name": name,
                "n_bulk": len(bulk[name]),
                "n_tot": len(total[name]),
                "volume": ts.volume,  # Å³
            })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("system", choices=SYSTEMS, metavar="system", help=f"One of: {', '.join(SYSTEMS)}")
    parser.add_argument("replicate", type=int, help="Replicate number (replex: 0-1, plain MD: 0-7)")
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
        df = count(path, args.stride)
        tmp = output.with_name(output.name + ".tmp")
        df.to_parquet(tmp)
        tmp.rename(output)
        print(f"[{path}] {df['time'].nunique()} frames in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
