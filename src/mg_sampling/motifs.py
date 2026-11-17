"""Ion binding motifs (MgRNA nomenclature) and the `ion_motif_analysis` command."""

import argparse
from collections import defaultdict, Counter

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis import Universe
from MDAnalysis.core.groups import Atom, AtomGroup
from MDAnalysis.lib.distances import calc_angles, distance_array
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.results import ResultsGroup

from mg_sampling.stats import block_bootstrap
from mg_sampling.structure import standardize_names

MOTIF_INFO = ['n_inner', 'inner_indices', 'inner_names', 'inner_resids', 'inner_resnames', 'stereo', 'mgrna']

MGRNA_LIGTYPES = ['Oph', 'Ob', 'Or', 'Nb', 'X', 'Ow']
MGRNA_LIGTYPES_DICT = defaultdict(
    lambda: 'X',
    {
        **{name: 'Oph' for name in ['OP1', 'OP2', 'O1P', 'O2P']},
        **{name: 'Ob' for name in ['O2', 'O4', 'O6']},
        **{name: 'Ow' for name in ['OW']},
        **{name: 'Nb'for name in ['N1', 'N3', 'N7']},
        **{name: 'Or' for name in ["O2'", "O4'", "O3'", "O5'"]},
    },
)

def get_mgrna_classification(stereo: str, atomnames: list[str]) -> str:
    ligtypes = [MGRNA_LIGTYPES_DICT[name] for name in atomnames]
    counter = Counter(ligtypes)
    if stereo:
        stereo += '-'
    label = stereo + ''.join(f"{counter[lig] if counter[lig] > 1 else ''}{lig}" for lig in MGRNA_LIGTYPES if counter[lig])
    return label

def get_stereochemistry(a: Atom, ag: AtomGroup) -> str:
    """Gets the stereochemistry (cis-, trans-, fac-, mer-, Cis-, Trans-) designation of 
    a group of ligands `ag` around a central atom `a`, assuming an octahedral geometry."""
 
    match tuple(ag):
        case (b, c):
            box = a.universe.dimensions
            if calc_angles(b.position, a.position, c.position, box) <= 2.3562:  # 135° = (90 + 180) / 2
                return 'cis'
            else:
                return 'trans'
        case (b, c, d):
            box = a.universe.dimensions
            angle_sum = (
                calc_angles(b.position, a.position, c.position, box)
                + calc_angles(b.position, a.position, d.position, box)
                + calc_angles(c.position, a.position, d.position, box)
            )
            if angle_sum <= 5.4978:  # 315° = (90+90+90 + 90+90+180) / 2
                return 'fac'
            else:
                return 'mer'
        case (b, c, d, e):
            box = a.universe.dimensions
            angle_sum = (
                calc_angles(b.position, a.position, c.position, box)
                + calc_angles(b.position, a.position, d.position, box)
                + calc_angles(b.position, a.position, e.position, box)
                + calc_angles(c.position, a.position, d.position, box)
                + calc_angles(c.position, a.position, e.position, box)
                + calc_angles(d.position, a.position, e.position, box)
            )
            if angle_sum <= 11.7810:  # 675° = ((5×90+180) + (4×90+2×180)) / 2
                return 'Cis'
            else:
                return 'Trans'
        case _:
            return ''

class IonNucleicAnalysis(AnalysisBase):
    """Analyze ion binding to RNA (though it could be applied to anything really).

    Classifies the binding motifs according to MgRNA nomenclature, e.g.:
        - 'Oph': Bound to one non-bridging phosphate oxygen atom
        - 'fac-3Oph': Bound to 3 non-bridging phosphate oxygen atoms with `fac` stereoisomery
        - etc.

    RNA ligand types:
        - Oph: non-bridging phosphate oxygen atoms (O1P, O2P)
        - Ob: base oxygen
        - Or: ribose oxygen
        - Nb: base nitrogen
        - Ow: water oxygen
        - X: other
    
    """

    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(cls):
        return ('serial', 'multiprocessing', 'dask',)

    def _get_aggregator(self):
        return ResultsGroup(lookup={'timeseries': ResultsGroup.flatten_sequence})

    def __init__(
            self,
            universe: Universe,
            ionsel: str = 'name MG',
            ligsel: str = 'nucleic',
            inner_cutoff: float = 2.5,
            outer_cutoff: float | None = None,
            output_distances: bool = False,
            keep_empty: bool = True,
            verbose: bool = True,
            ):
        """
        Parameters
        ----------
        universe : Universe
        ionsel: str, optional ('name MG')
            selection of "ions", central particles whose inner (and outer) shell will be analyzed
        ligsel: str, optional ('nucleic')
            selection whose atoms will be treated as ligands for `ionsel`
        inner_cutoff: float, optional (2.5 Angstroms)
            cutoff radius for the inner (first) solvation shell
        outer_cutoff: float | None, optional (None)
            cutoff radius for the outer (second) solvation shell
        output_distances: bool, optional (False)
            save distances to ligands
        keep_empty: bool, optional (True)
            whether to keep entries where no ligand is detected
        verbose: bool, optional (True)
            show progress
        """
        super(IonNucleicAnalysis, self).__init__(universe.trajectory, verbose=verbose)

        self.u = universe
        self.ionsel = ionsel
        self.ligsel = ligsel

        if (inner_cutoff <= 0) or ((outer_cutoff is not None) and (outer_cutoff <= inner_cutoff)):
            raise ValueError("Must have 0 < `inner_cutoff` < `outer_cutoff`")
        self.inner_cutoff = inner_cutoff
        self.outer_cutoff = outer_cutoff
        self.output_distances = output_distances
        self.keep_empty = keep_empty
    
    def _prepare(self):
        self._ions = self.u.select_atoms(self.ionsel)
        if not self._ions:
            raise ValueError("Ion selection must contain at least one atom.")

        cutoff = self.outer_cutoff if self.outer_cutoff is not None else self.inner_cutoff

        # This selection is costly and is updated every frame
        self._ag = self.u.select_atoms(f'(({self.ligsel}) and around {cutoff} ({self.ionsel})) or ({self.ionsel})', updating=True)

        # These subselections are less costly and also updated every frame
        self._ags_lig_inner = [
            self._ag.select_atoms(f'({self.ligsel}) and around {self.inner_cutoff} (index {a.index})',
                                  updating=True)
            for a in self._ions
        ]

        if self.outer_cutoff is not None:
            self._ags_lig_outer = [
                self._ag.select_atoms(
                    f'({self.ligsel}) and sphlayer {self.inner_cutoff} {self.outer_cutoff} (index {a.index})',
                    updating=True)
                for a in self._ions
            ]
        self.results.timeseries = []

    def _single_frame(self):
        frame = []
        for i, a in enumerate(self._ions):
            ag_lig_inner = self._ags_lig_inner[i]
            ag_lig_outer = self._ags_lig_outer[i] if self.outer_cutoff is not None else None

            if self.keep_empty or ag_lig_inner or ag_lig_outer:
                stereo = get_stereochemistry(a, ag_lig_inner)
                data = {
                    'frame': self._ts.frame,
                    'time': self._ts.time,
                    'index': a.index,
                    'name': a.name,
                    'n_inner': len(ag_lig_inner),
                    'inner_indices': ag_lig_inner.indices.tolist(),
                    'inner_names': ag_lig_inner.names.tolist(),
                    'inner_resids': ag_lig_inner.resids.tolist(),
                    'inner_resnames': ag_lig_inner.resnames.tolist(),
                    'stereo': stereo,
                    'mgrna': get_mgrna_classification(stereo, ag_lig_inner.names.tolist()),
                }

                if ag_lig_outer is not None:
                    data.update({
                        'n_outer': len(ag_lig_outer),
                        'outer_indices': ag_lig_outer.indices.tolist(),
                        'outer_names': ag_lig_outer.names.tolist(),
                        'outer_resids': ag_lig_outer.resids.tolist(),
                        'outer_resnames': ag_lig_outer.resnames.tolist(),
                    })

                if self.output_distances:
                    inner_dist = distance_array(a.position, ag_lig_inner.positions, box=self._ts.dimensions).reshape(-1)
                    outer_dist = distance_array(a.position, ag_lig_outer.positions, box=self._ts.dimensions).reshape(-1)
                    data.update({
                        'inner_dist': inner_dist.tolist(),
                        'outer_dist': outer_dist.tolist(),
                    })

                frame.append(data)
        self.results.timeseries.append(frame)

    def to_pandas(self):
        df = pd.DataFrame([data for frame in self.results.timeseries for data in frame])
        return df


def motif_names(df):
    """Motif of each row of an `IonNucleicAnalysis` table: "{resname}{resid}{name},…-{mgrna}"."""
    keys = pd.Series(list(zip(
        df['inner_resnames'].map(tuple), df['inner_resids'].map(tuple), df['inner_names'].map(tuple), df['mgrna'],
    )))
    codes, uniques = pd.factorize(keys)
    names = []
    for resnames, resids, atoms, mgrna in uniques:
        name = ",".join(f"{resname}{resid}{atom}" for resname, resid, atom in zip(resnames, resids, atoms))
        names.append(f"{name}-{mgrna}" if mgrna else name)
    return pd.Series(np.asarray(names, dtype=object)[codes], index=df.index)


def compress_motifs(df):
    """
    Motif occupancy along one trajectory and the table of motifs, from an `IonNucleicAnalysis` table.

    occupancy: frames (index `time`) x motifs, number of ions in each motif (inner sphere).
    info: one row per motif with its inner ligands, stereochemistry, MgRNA class and the ions seen in it.
    """
    bound = df[df['n_inner'] > 0].assign(motif=lambda x: motif_names(x))
    occupancy = (
        bound.groupby(['time', 'motif']).size()
        .unstack('motif', fill_value=0)
        .reindex(np.unique(df['time']), fill_value=0)
        .rename_axis('time')
        .astype(np.uint16)
    )
    info = bound.drop_duplicates('motif').set_index('motif')[MOTIF_INFO]
    info['mg_indices'] = bound.groupby('motif')['index'].agg(lambda x: sorted(set(x)))
    return occupancy, info


def summarize_motifs(blocks, weights=None, rng=0, **bootstrap_kwargs):
    """
    Mean occupancy `p` of each motif over `blocks` [(occupancy, info), …] with a bayesian block
    bootstrap interval (`p_low`, `p_high`), joined with the motif table. `weights`: optional
    frame weights for each block (biased runs). Blocks without frames are ignored.
    """
    occupancies, infos = zip(*blocks)
    if weights is None:
        weights = [np.ones(len(occupancy)) for occupancy in occupancies]
    motifs = sorted(set().union(*(occupancy.columns for occupancy in occupancies)))
    block_means, block_weights = [], []
    for occupancy, w in zip(occupancies, weights):
        if np.sum(w) > 0:
            block_means.append(np.average(occupancy.reindex(columns=motifs, fill_value=0), axis=0, weights=w))
            block_weights.append(np.sum(w))
    res = block_bootstrap(
        np.array(block_means), lambda x, weights: np.average(x, axis=0, weights=weights),
        weights=block_weights, bayesian=True, rng=rng, **bootstrap_kwargs,
    )
    summary = pd.DataFrame(
        {'p': res.estimate, 'p_low': res.confidence_interval.low, 'p_high': res.confidence_interval.high},
        index=pd.Index(motifs, name='motif'),
    )
    info = pd.concat(infos)
    mg_indices = info.groupby(level='motif')['mg_indices'].agg(lambda x: sorted(set().union(*x)))
    info = info[~info.index.duplicated()].drop(columns='mg_indices').join(mg_indices)
    return summary.join(info)


def main():
    parser = argparse.ArgumentParser(
        prog="ion_motif_analysis",
        description="Anaylize ion-RNA contacts and classify them according to MgRNA nomenclature",
    )
    parser.add_argument('topfile', type=str, help="MDAnalysis-readable file for topology")
    parser.add_argument('trajfiles', nargs='*', type=str, help="MDAnalysis-readable file(s) for trajectory")
    parser.add_argument('--sel', type=str, default='name MG',
        help="('name MG') MDAnalysis selection string for ions")
    parser.add_argument('--ligsel', type=str, default='(nucleic and not type H) or name MG',
        help="('(nucleic and not type H) or name MG') MDAnalysis selection string for binding sites")
    parser.add_argument('--inner_cutoff', type=float, default=2.5, help="(2.5) Inner-sphere cutoff radius in Å")
    parser.add_argument('--outer_cutoff', type=float, default=4.8, help="(4.8) Outer-sphere cutoff radius in Å; set to zero or negative to skip outer-sphere analysis")
    parser.add_argument('--output-distances', action='store_true', help="Store distances to inner (and outer) ligand atoms")
    parser.add_argument('--drop-empty', action='store_true', help="Discard frame for ion with no ligand")
    parser.add_argument('--standardize', action='store_true', help="Standardize some atom names (O1P->OP1…)")
    parser.add_argument('-o', '--outfile', type=str, default='ion_analysis.csv', help="(ion_analysis.csv) Output file (.csv or .parquet)")
    parser.add_argument('--start', type=int, default=0, help="(0) Starting frame")
    parser.add_argument('--stop', type=int, default=None, help="(None) Ending frame")
    parser.add_argument('--step', type=int, default=1, help="(1) Subsample every N frame")
    parser.add_argument('-n', '--n_workers', type=int, default=None, help="(None) Number of parallel processes")

    args = parser.parse_args()
    if args.outer_cutoff <= 0:
        args.outer_cutoff = None

    u = mda.Universe(args.topfile, args.trajfiles) if args.trajfiles else mda.Universe(args.topfile)
    if args.standardize:
        standardize_names(u)
    anal = IonNucleicAnalysis(
        u,
        ionsel=args.sel,
        inner_cutoff=args.inner_cutoff,
        ligsel=args.ligsel,
        outer_cutoff=args.outer_cutoff,
        output_distances=args.output_distances,
        keep_empty=(not args.drop_empty),
    )
    print(f"Running analysis with {args.n_workers} processes...")
    if args.n_workers is not None:
        anal.run(start=args.start, stop=args.stop, step=args.step, backend='dask', n_workers=args.n_workers)
    else:
        anal.run(start=args.start, stop=args.stop, step=args.step, backend='serial')
    df = anal.to_pandas()

    if args.outfile.endswith(".parquet"):
        try:
            df.to_parquet(args.outfile)
        except Exception as e:
            print("Could not save parquet format:", e)
            print("Output as CSV instead:", args.outfile)
            df.to_csv(args.outfile)
    else:
        df.to_csv(args.outfile)

if __name__ == "__main__":
    main()
