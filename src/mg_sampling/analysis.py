from collections import defaultdict, Counter

from MDAnalysis import Universe
from MDAnalysis.core.groups import Atom, AtomGroup
from MDAnalysis.lib.distances import calc_angles, distance_array
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.results import ResultsGroup

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
            if calc_angles(b.position, a.position, c.position, box) <= 135:
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
            if angle_sum <= 202.5:
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
            if angle_sum <= 472.5:
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
            if ag_lig_inner or ag_lig_outer:
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
        import pandas as pd
        df = pd.DataFrame([data for frame in self.results.timeseries for data in frame])
        return df

