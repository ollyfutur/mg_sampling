#!/usr/bin/env python3

import argparse
import numpy as np

from mg_sampling import PROJECT_DIR
from mg_sampling.analysis import IonNucleicAnalysis
import MDAnalysis as mda

def main():
    parser = argparse.ArgumentParser(
        prog="ion_motif_analysis",
        description="Anaylize ion-RNA contacts and classify them according to MgRNA nomenclature",
    )
    parser.add_argument('topfile', type=str, help="MDAnalysis-readable file for topology")
    parser.add_argument('trajfiles', nargs='*', type=str, help="MDAnalysis-readable file(s) for trajectory")
    parser.add_argument('--sel', type=str, default='name MG NA CL',
        help="('name MG NA CL') MDAnalysis selection string for ions")
    parser.add_argument('--ligsel', type=str, default='(nucleic and not type H) or name MG NA CL',
        help="('(nucleic and not type H) or name MG NA CL') MDAnalysis selection string for binding sites")
    parser.add_argument('--inner_cutoff', type=float, default=2.5, help="(2.5) Inner-sphere cutoff radius in Å")
    parser.add_argument('--outer_cutoff', type=float, default=4.8, help="(4.8) Outer-sphere cutoff radius in Å; set to zero or negative to skip outer-sphere analysis")
    parser.add_argument('-o', '--outfile', type=str, default='ion_analysis.csv', help="(ion_analysis.csv) Output CSV file")
    parser.add_argument('--start', type=int, default=0, help="(0) Starting frame")
    parser.add_argument('--stop', type=int, default=None, help="(None) Ending frame")
    parser.add_argument('--step', type=int, default=1, help="(1) Subsample every N frame")
    parser.add_argument('-n', '--n_workers', type=int, default=1, help="(1) Number of parallel processes")

    args = parser.parse_args()
    if args.outer_cutoff <= 0:
        args.outer_cutoff = None

    u = mda.Universe(args.topfile, args.trajfiles)
    anal = IonNucleicAnalysis(u, args.sel, args.inner_cutoff, ligsel=args.ligsel, outer_cutoff=args.outer_cutoff)
    print(f"Running analysis with {args.n_workers} processes...")
    anal.run(start=args.start, stop=args.stop, step=args.step, backend='dask', n_workers=args.n_workers)
    df = anal.to_pandas()
    df.to_csv(args.outfile)

if __name__ == "__main__":
    main()

