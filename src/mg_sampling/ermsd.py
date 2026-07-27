import numpy as np
from scipy.spatial import distance

import warnings

from MDAnalysis import AtomGroup, Universe
from MDAnalysis.analysis.base import AnalysisBase, ResultsGroup
from MDAnalysis.exceptions import SelectionError, SelectionWarning

# Geometric definitions
F_FACTORS = np.array([5, 5, 3])  # Å


class eRMSD(AnalysisBase):
    @classmethod
    def get_supported_backends(cls):
        return (
            "serial",
            "multiprocessing",
            "dask",
        )

    _analysis_algorithm_is_parallelizable = True

    def _get_aggregator(self):
        return ResultsGroup(lookup={"ermsd": ResultsGroup.ndarray_hstack})

    def __init__(
        self,
        atomgroup: AtomGroup | Universe,
        reference: AtomGroup | Universe | None = None,
        cutoff: float = 2.4,
        ref_frame: int = 0,
        **kwargs,
    ):
        """
        Compute the eRMSD [1], a contact map-based distance suitable for the
        comparison of nucleic acid structures.

        [1] Bottaro S, Di Palma F, Bussi G. The role of nucleobase interactions in RNA structure and
        dynamics. Nucleic acids research. 2014; 42(21):13306–13314 (https://doi.org/10.1093/nar/gku972)

        Parameters
        ----------
        atomgroup : AtomGroup or Universe
            Group of atoms for which the εRMSD is calculated. It must
            be composed of nucleic residues.
        reference : AtomGroup or Universe (optional)
            Group of reference atoms; if ``None`` then the first frame
            of `atomgroup` is used (unless `ref_frame` is set).
        ref_frame : int (optional)
            frame index to select frame from `reference`
        """
        super().__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup.atoms
        self.reference = reference.atoms if reference is not None else self.atomgroup
        self.cutoff = cutoff
        self.ref_frame = ref_frame

        self.mobile_atoms = _get_lcs_atomgroup(atomgroup)
        self.ref_atoms = (
            _get_lcs_atomgroup(reference)
            if reference is not None
            else self.mobile_atoms
        )

        if len(self.ref_atoms) != len(self.mobile_atoms):
            msg = (
                "Reference and trajectory atom selections for local coordinate system "
                "do not contain the same number of atoms: "
                "N_ref={0:d}, N_traj={1:d}".format(
                    self.ref_atoms.n_atoms, self.mobile_atoms.n_atoms
                )
            )
            raise SelectionError(msg)

        self._n_residues = len(self.ref_atoms.residues)

    def _prepare(self):
        current_frame = self.reference.universe.trajectory.ts.frame
        try:
            # Go to reference frame
            self.reference.universe.trajectory[self.ref_frame]
            # Store reference local coordinate system
            self._ref_gmat = calc_gmat(
                self.ref_atoms.positions.reshape(3, -1, 3), self.cutoff
            )
        finally:
            self.reference.universe.trajectory[current_frame]

        self.results.ermsd = np.zeros(self.n_frames)

    def _single_frame(self):
        gmat = calc_gmat(self.mobile_atoms.positions.reshape(3, -1, 3), self.cutoff)

        self.results.ermsd[self._frame_index] = distance.euclidean(
            gmat.flatten(), self._ref_gmat.flatten()
        ) / np.sqrt(self._n_residues)


LCS_ATOMS = {
    # Purines
    "A": ("C2", "C6", "C4"),
    "G": ("C2", "C6", "C4"),
    # Pyrimidines
    "C": ("C2", "C4", "C6"),
    "T": ("C2", "C4", "C6"),
    "U": ("C2", "C4", "C6"),
}


def _get_lcs_atomgroup(atomgroup, lcs_atoms=None):
    if lcs_atoms is not None:
        lcs_atoms = dict(LCS_ATOMS).update(lcs_atoms)
    else:
        lcs_atoms = LCS_ATOMS

    unknown_residues = ~np.isin(atomgroup.residues.resnames, list(lcs_atoms))
    if np.any(unknown_residues):
        warnings.warn(
            "Ignoring unknown nucleic residues: %s" % atomgroup.residues.resnames[unknown_residues],
            SelectionWarning,
        ) 

    sel_strs = (" or ".join(
        f"(resname {resname} and name {lcs_atoms[resname][i]})"
        for resname in lcs_atoms
    ) for i in range(3))

    ag0, ag1, ag2 = (atomgroup.select_atoms(sel_str) for sel_str in sel_strs)

    if len(ag0) != len(ag1) or len(ag0) != len(ag2):
        msg = f"Mismatch in atom selections for local coordinate system. {len(ag0) = }, {len(ag1) = }, {len(ag2) = }."
        raise SelectionError(msg)

    return AtomGroup(ag0 + ag1 + ag2)

def calc_lcs(coords):
    """
    Calculate local coordinates system

    Calculate origin position and xyz orthonormal vectors in the six-member rings

    Parameters
    ----------
    coords : (3,n,3) numpy array
        (3,n,3) numpy array with positions of C2,C4 and C6 atoms for pyrimidines (C,U,T) and C2,C6,C4 for purines (A,G) (axis 0) relative to n nucleobases (axis 1). xyz coordinates in axis 2.

    Returns
    -------
    lcs : (n,3,3) numpy array
        x y z vectors defining a local coordinate system in the geometrical center of the nucleobase
    origo : (n,3) numpy array
        origin  of local coordinate systems
    """

    # calculate center of mass
    origo = np.mean(coords, axis=0)

    # CoM-C2 (x axis)
    x = coords[0] - origo
    x_norm = np.sqrt(np.sum(x * x, axis=1))
    x = x / x_norm[:, np.newaxis]
    # CoM-C4/C6
    c = coords[1] - origo
    # z/y axis
    z = np.cross(x, c, axis=1)
    z_norm = np.sqrt(np.sum(z * z, axis=1))
    z = z / z_norm[:, np.newaxis]

    y = np.cross(z, x, axis=1)
    lcs = np.array([x.T, y.T, z.T]).T
    return lcs, origo


def calc_3dmat(coords, cutoff):
    """
    Calculate relative position of nucleobases within an ellipsoidal cutoff

    Parameters
    ----------
    coords : (3,n,3) numpy array
        (3,n,3) numpy array with positions of C2,C4 and C6 atoms for pyrimidines (C,U,T) and C2,C6,C4 for purines (A,G) (axis 0) relative to n nucleobases (axis 1). xyz coordinates in axis 2.

    cutoff : float
       ellipsoidal cutoff

    Returns
    -------
    dotp : (x,3) numpy array
       xyz coordinates for each pair

    m_idx : (x,2) numpy array
       indeces of the pair
    """

    # prune search first
    lcs, origo = calc_lcs(coords)
    max_r = np.max(F_FACTORS) * cutoff
    dmat = distance.squareform(distance.pdist(origo))
    m_idx = np.array(np.where((dmat < max_r) & (dmat > 0.01))).T

    # calculate scaled distances
    # diff = [origo[y]-origo[x] for x,y in m_idx]
    diff = origo[m_idx[:, 1]] - origo[m_idx[:, 0]]

    dotp = np.array(
        [np.dot(diff[i], lcs[j]) for i, j in zip(range(len(diff)), m_idx[:, 0])]
    )
    return dotp, m_idx


def calc_gmat(coords, cutoff):
    """
    Calculate G-vectors for each pair of bases within ellipsoidal cutoff distance

    Parameters
    ----------
    coords : (3,n,3) numpy array
        (3,n,3) numpy array with positions of C2,C4 and C6 atoms for pyrimidines (C,U,T) and C2,C6,C4 for purines (A,G) (axis 0) relative to n nucleobases (axis 1). xyz coordinates in axis 2.

    cutoff : float
        ellipsoidal cutoff

    Returns
    -------
    dotp : (n,n,4) numpy array
        G coordinates for each pair. For pairs outside the cutoff the coordinates are (0,0,0,0)
    """

    ll = coords.shape[1]

    mat = np.zeros((ll, ll, 4))

    dotp, m_idx = calc_3dmat(coords, cutoff)

    # return zero matrix when there are no contacts
    if dotp.shape[0] == 0:
        return mat

    dotp /= F_FACTORS[np.newaxis, :]
    dotp_norm = np.sqrt(np.sum(dotp**2, axis=1))

    # calculate 4D g-vector
    ff = (np.pi * dotp_norm) / cutoff
    factor13 = np.sin(ff) / ff
    factor4 = ((1.0 + np.cos(ff)) * cutoff) / np.pi
    gmat = dotp * factor13[:, np.newaxis]
    gmat = np.concatenate((gmat, factor4[:, np.newaxis]), axis=1)

    # set to zero when norm is larger than cutoff
    gmat[dotp_norm > cutoff] = 0.0

    mat[m_idx[:, 0], m_idx[:, 1]] = gmat

    return mat
