"""Atom naming, ordering and periodic wrapping helpers for MDAnalysis atom groups."""

import warnings

import numpy as np

import MDAnalysis as mda
from MDAnalysis import AtomGroup, Universe
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.lib.distances import minimize_vectors

STANDARD_NAMES = {
    "O1P": "OP1",
    "O2P": "OP2",
}


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


def unwrap_sequential(ag: AtomGroup, inplace=True):
    pos = ag.positions
    Δpos = np.diff(pos, axis=0, prepend=pos[[0]])
    Δpos = mda.lib.distances.minimize_vectors(Δpos, ag.dimensions)
    new_pos = pos[0] + np.cumsum(Δpos, axis=0)
    if inplace:
        ag.positions = new_pos
    return new_pos


def wrap_around(ag1, position):
    """Wraps atoms in `ag1` in the unit cell centered at `position`"""
    position = np.asarray(position)
    vectors = ag1.positions - position
    min_vectors = minimize_vectors(vectors, ag1.ts.dimensions)
    ag1.positions = position + min_vectors


def wrap_around_atomgroup(ag1, ag2):
    """Wraps the atoms in `ag1` to the periodic image
    minimizing the distance to `ag2`."""
    d = distance_array(ag1, ag2, box=ag1.ts.dimensions)
    closest_j = np.argmin(d, axis=1)
    pos_j = ag2.positions[closest_j]
    wrap_around(ag1, pos_j)
