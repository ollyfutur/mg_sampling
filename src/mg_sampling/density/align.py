"""
Local or global alignment.
"""

import warnings
import argparse
from pathlib import Path
import numpy as np

import MDAnalysis as mda
from MDAnalysis import AtomGroup, Universe
from MDAnalysis.analysis import align
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.results import ResultsGroup
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.transformations.base import TransformationBase

from tqdm import tqdm

STANDARD_NAMES = {
    "O1P": "OP1",
    "O2P": "OP2",
}

def unwrap_sequential(ag: AtomGroup, inplace=True):
    pos = ag.positions
    Δpos = np.diff(pos, axis=0, prepend=pos[[0]])
    Δpos = mda.lib.distances.minimize_vectors(Δpos, ag.dimensions)
    new_pos = pos[0] + np.cumsum(Δpos, axis=0)
    if inplace:
        ag.positions = new_pos
    return new_pos


def standardize_names(ag: AtomGroup | Universe):
    """
    Standardize atom names according to PDB standard residues
    """
    ag = ag.atoms
    ag.names = [
        STANDARD_NAMES[name] if name in STANDARD_NAMES else name for name in ag.names
    ]


def reorder_residue_atoms(ag):
    """
    Reorder atoms in `ag` by resindex, then by name.
    """
    order = np.lexsort((ag.names, ag.resindices))
    return ag[order]


def check_alignable(mob_align, ref_align):
    if len(mob_align) != len(ref_align):
        raise ValueError(
            "Mobile and reference alignment selection should contain the same number of atoms. "
            f"({len(mob_align)} != {len(ref_align)})"
        )
    if not np.all(mob_align.names == ref_align.names):
        diff = mob_align.names != ref_align.names
        warnings.warn(
            "Mobile and reference alignment selections have the same number of atoms but different "
            "atom names, please make sure this is expected! "
            "Mob: " + " ".join(mob_align.names[diff][:10]) + "… "
            "Ref: " + " ".join(ref_align.names[diff][:10]) + "… "
        )


def wrap_around(ag1, position):
    """Wraps atoms in `ag1` in the unit cell centered at `position`"""
    position = np.asarray(position)
    box_center = ag1.ts.triclinic_dimensions.sum(axis=0) / 2
    ag1.translate(box_center - position)
    ag1.wrap()
    ag1.translate(position - box_center)


def wrap_around_atomgroup(ag1, ag2):
    """Wraps the atoms in `ag1` to the periodic image
    minimizing the distance to `ag2`."""
    d = distance_array(ag1, ag2, box=ag1.ts.dimensions)
    closest_j = np.argmin(d, axis=1)
    pos_j = ag2.positions[closest_j]
    wrap_around(ag1, pos_j)


def global_align(
    mobile: AtomGroup | Universe,
    reference: AtomGroup | Universe,
    align_sel: str = "nucleic",
    ref_align_sel: str | None = None,
    out_sel: str = "all",
    out_radius: float | None = None,
    reorder: bool = False,
    verbose: bool = False,
):
    """
    Aligns `mobile` to `reference` based on `align_sel` (and `ref_align_sel`),
    then returns `out_sel` in a radius of `out_radius` around the mobile's `align_sel`.
    (set `out_radius` to `None` to output the entire `out_sel`)
    Warning: make sure the ref and mobile are whole (not cut across boundaries). The
    PBC will be ignored altogether.
    """
    mobile = mobile.atoms
    reference = reference.atoms

    ref_align_sel = align_sel if ref_align_sel is None else ref_align_sel

    mob_align = mobile.select_atoms(align_sel, sorted=False)
    ref_align = reference.select_atoms(ref_align_sel, sorted=False)
    if reorder:
        mob_align = reorder_residue_atoms(mob_align)
        ref_align = reorder_residue_atoms(ref_align)

    check_alignable(mob_align, ref_align)

    mob_out = mobile.select_atoms(out_sel, sorted=False)
    if out_radius is not None:
        # NB: `around` is exclusive
        mob_out = mob_out.select_atoms(
            f"group mob_align or around {out_radius} group mob_align",
            mob_align=mob_align,
            sorted=False,
            periodic=False,
        )

    _ = align.alignto(mob_align, ref_align, subselection=mob_out)

    return mob_out


def local_align(
    mobile: AtomGroup | Universe,
    reference: AtomGroup | Universe,
    align_sel: str = "nucleic",
    ref_align_sel: str | None = None,
    align_radius: float | None = None,
    out_sel: str = "all",
    out_radius: float | None = None,
    reorder: bool = False,
    verbose: bool = False,
) -> mda.Universe:
    """
    Creates a new `MDAnalysis.Universe` by mapping local parts of a given `mobile` structure onto a `reference` structure.
    This is the algorithm:
    For each residue `i` of the `align_sel` selection:
        - The selection is expanded with `align_radius` (in `reference`) and used to align
          the corresponding `mobile` structure.
        - The mobile's `out_sel` around `out_radius` of the residue `i` is added
    If atoms are repeated (they are within the `out_radius` of multiple residues) their `occupancy`
    is divided accordingly.
    Warning: make sure the ref and mobile are whole (not cut across boundaries). The
    PBC will be ignored altogether.
    """

    # Non-unique attributes
    ATOM_ATTRNAMES = [
        "names",
        "ids",
        "types",
        "masses",
        "chainIDs",
        "elements",
        "tempfactors",
    ]

    RES_ATTRNAMES = [
        "resids",
        "resnums",
        "resnames",
    ]

    mobile = mobile.atoms
    reference = reference.atoms

    if reorder:
        mobile = reorder_residue_atoms(mobile)
        reference = reorder_residue_atoms(reference)

    ref_align_sel = align_sel if ref_align_sel is None else ref_align_sel

    mob_align = mobile.select_atoms(align_sel, sorted=False)
    ref_align = reference.select_atoms(ref_align_sel, sorted=False)
    check_alignable(mob_align, ref_align)

    # Topology attribute lists
    indices = []
    resindices = []
    positions = []
    attrs = {
        **{attrname: [] for attrname in ATOM_ATTRNAMES if hasattr(mobile, attrname)},
        **{attrname: [] for attrname in RES_ATTRNAMES if hasattr(mobile, attrname)},
    }

    # Could be extended to atom-wise alignment
    mob_out_or_align = mobile.select_atoms(
        f"({align_sel}) or ({out_sel})", sorted=False
    )
    iter = tqdm(ref_align.residues) if verbose else ref_align.residues
    for r in iter:
        # Select reference atoms
        _ref = ref_align.select_atoms(
            f"resindex {r.resindex} or around {align_radius} resindex {r.resindex}",
            periodic=False,
            sorted=False,
        )
        # Get corresponding mobile atoms
        # We match indexing inside `ref_align` and `mob_align`
        _, _, _inds = np.intersect1d(
            _ref.indices, ref_align.indices, assume_unique=True, return_indices=True
        )
        _mob = mob_align[_inds]

        # Select the atoms that need to be added
        # _mob_sel = f"index {' '.join(map(str, _mob.indices))}"
        # _toadd = mob_out_or_align.select_atoms(
        #     f"(around {out_radius} ({_mob_sel})) and ({out_sel})",
        #     sorted=False,
        # )
        _ref_r = r.atoms
        _, _, _inds = np.intersect1d(
            _ref_r.indices,
            ref_align.indices,
            assume_unique=True,
            return_indices=True,
        )
        _out_r = mob_align[_inds]
        _out_r_sel = f"index {' '.join(map(str, _out_r.indices))}"
        _toadd = mob_out_or_align.select_atoms(
            f"((around {out_radius} ({_out_r_sel})) or ({_out_r_sel})) and ({out_sel})",
            periodic=False,
            sorted=False,
        )
        # Trick: only apply alignment to _toadd speeds things up...
        # ... but we need to reverse it immediately after (*)
        _temp_pos = _toadd.positions.copy()
        _ = align.alignto(_mob, _ref, subselection=_toadd)

        indices.append(_toadd.indices)
        resindices.append(_toadd.resindices)
        positions.append(_toadd.positions)

        # (*)
        _toadd.positions = _temp_pos

        for attrname in attrs:
            attrs[attrname].append(getattr(_toadd, attrname))

    indices = np.concatenate(indices)
    n_atoms = len(indices)

    # Reorder based on original indices
    ordering = np.argsort(indices)
    indices = indices[ordering]
    resindices = np.concatenate(resindices)[ordering]
    positions = np.concatenate(positions)[ordering]
    for attrname in attrs:
        attrs[attrname] = np.concatenate(attrs[attrname])[ordering]

    # Account for repeated atoms using "occupancy"
    _, unique_inverse, unique_counts = np.unique(
        indices, return_inverse=True, return_counts=True
    )
    occupancies = 1 / unique_counts[unique_inverse]

    # Get unique resindices in new topology
    res_unique, res_unique_indices, res_unique_inverse = np.unique(
        resindices,
        return_index=True,
        return_inverse=True,
    )
    n_residues = len(res_unique)

    # Build topology
    unew = mda.Universe.empty(
        n_atoms,
        n_residues,
        atom_resindex=res_unique_inverse,
        trajectory=True,
    )
    unew.add_TopologyAttr("occupancy", occupancies)

    for attrname in attrs:
        # Add per-atom attributes
        if attrname in ATOM_ATTRNAMES:
            unew.add_TopologyAttr(attrname, attrs[attrname])
        # Add per-residue attributes
        if attrname in RES_ATTRNAMES:
            unew.add_TopologyAttr(attrname, attrs[attrname][res_unique_indices])

    # Add coordinates
    unew.atoms.positions = positions
    unew.dimensions = mobile.dimensions

    return unew


class Alignment(AnalysisBase):
    @classmethod
    def get_supported_backends(cls):
        return ("serial", "multiprocessing", "dask")

    _analysis_algorithm_is_parallelizable = True

    def __init__(
        self,
        mobile: AtomGroup | Universe,
        reference: AtomGroup | Universe,
        mode: str = "global",
        align_sel: str = "nucleic",
        ref_align_sel: str | None = None,
        align_radius: float = 12,
        out_sel: str = "all",
        out_radius: float = 5,
        reorder: bool = False,
        outpath: str = "aligned.pdb",
        verbose: bool = False,
    ):
        super(Alignment, self).__init__(mobile.universe.trajectory, verbose=verbose)

        if mode not in ("global", "residue"):
            raise ValueError("`mode` must be 'global' or 'residue'")

        self.mobile = mobile
        self.reference = reference
        self.mode = mode
        self.align_sel = align_sel
        self.ref_align_sel = ref_align_sel
        self.align_radius = align_radius
        self.out_sel = out_sel
        self.out_radius = out_radius
        self.reorder = reorder
        self.outpath = Path(outpath)
        self.verbose = verbose

    def _prepare(self):
        self.results.outfiles = []

    def _single_frame(self):
        if self.mode == "residue":
            self._single_frame_local()
        elif self.mode == "global":
            self._single_frame_global()

    def _single_frame_local(self):
        unew = local_align(
            self.mobile,
            self.reference,
            align_sel=self.align_sel,
            ref_align_sel=self.ref_align_sel,
            align_radius=self.align_radius,
            out_sel=self.out_sel,
            out_radius=self.out_radius,
            reorder=self.reorder,
            verbose=self.verbose,
        )
        outfile = self.outpath.with_stem(f"{self.outpath.stem}_{self._ts.frame}")
        with mda.Writer(outfile, multiframe=False) as W:
            W.write(unew.atoms)
            del W.obj  # This is needed to avoid a memory leak
        self.results.outfiles.append(outfile)

    def _single_frame_global(self):
        unew = global_align(
            self.mobile,
            self.reference,
            align_sel=self.align_sel,
            ref_align_sel=self.ref_align_sel,
            out_sel=self.out_sel,
            out_radius=self.out_radius,
            reorder=self.reorder,
            verbose=self.verbose,
        )
        outfile = self.outpath.with_stem(f"{self.outpath.stem}_{self._ts.frame}")
        unew.atoms.write(outfile)
        self.results.outfiles.append(outfile)

    def _get_aggregator(self):
        return ResultsGroup(lookup={"outfiles": ResultsGroup.flatten_sequence})


def main():
    parser = argparse.ArgumentParser(
        description="Align a structure onto another by superposing globally or "
        "aligned residue-wise local environments. "
        "The alignment selections for reference and for input structures must match exactly. "
        "Repeated atoms are weighted accordingly in the occupancy field.",
    )

    parser.add_argument("input", type=str, help="Structure to align")
    parser.add_argument("--traj", type=str, nargs="*", help="Trajectory file(s)")
    parser.add_argument("--start", type=int, help="Trajectory starting frame")
    parser.add_argument("--stop", type=int, help="Trajectory end frame (exclusive)")
    parser.add_argument("--step", type=int, help="Trajectory frame stride")
    parser.add_argument("--ref", type=str, help="Reference structure")
    parser.add_argument(
        "--mode",
        type=str,
        choices=["residue", "global"],
        default="residue",
        help="Alignment mode",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="aligned.pdb",
        help="(aligned.pdb) Output file",
    )
    parser.add_argument(
        "--sel",
        type=str,
        default="not type H and nucleicbackbone",
        help="('not type H and nucleicbackbone') Alignment MDAnalysis selection string",
    )
    parser.add_argument(
        "--refsel", type=str, help="Alignment selection for reference (if different)"
    )
    parser.add_argument(
        "--outsel",
        type=str,
        default="not type H and not nucleic",
        help="('not type H and not nucleic') Selection for output atoms after local alignment",
    )
    parser.add_argument(
        "--ralign",
        type=float,
        default=12,
        help="(12.0) Inclusion radius around residue for alignment selection [Å]",
    )
    parser.add_argument(
        "--rout",
        type=float,
        default=5,
        help="(5.0) Inclusion radius around residue for output selection [Å]",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        help="Number of processes for parallelizing trajectory alignment",
    )
    parser.add_argument(
        "--reorder",
        action="store_true",
        help="Reorder atoms within residues before alignment, usefull if ordering differs in reference and input files.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Show progress")
    args = parser.parse_args()

    u = (
        mda.Universe(args.input, *args.traj)
        if args.traj is not None
        else mda.Universe(args.input)
    )
    ref = mda.Universe(args.ref)

    analysis = Alignment(
        u,
        ref,
        mode=args.mode,
        outpath=args.output,
        align_sel=args.sel,
        ref_align_sel=args.refsel,
        align_radius=args.ralign,
        out_sel=args.outsel,
        out_radius=args.rout,
        reorder=args.reorder,
        verbose=args.verbose,
    )

    if args.processes is not None:
        analysis.run(
            backend="multiprocessing",
            n_workers=args.processes,
            start=args.start,
            step=args.step,
            stop=args.stop,
        )
    else:
        analysis.run(backend="serial", start=args.start, step=args.step, stop=args.stop)


if __name__ == "__main__":
    main()
