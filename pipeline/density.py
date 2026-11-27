"""Forward density maps of the simulations and their comparison with the cryo-EM maps."""

import functools
import warnings

import MDAnalysis as mda
import numpy as np
import pandas as pd
import torch
from cryomap import Grid, read_ensemble, read_mrc
from cryomap.masks import within
from scipy.interpolate import interp1d

from mg_sampling.stats import BootstrapResult, ConfidenceInterval
from pipeline import PROJECT_DIR

# Map families: rigid fit of the RNA, and local fit (see density/density_maps.py)
FAMILIES = {"global": "global", "local": "local/align10_out10"}

REF_PDB = PROJECT_DIR / "structures" / "9CBU_fixed.pdb"
RNA_SEL = "nucleic and not element H and not (resid 22 and name P OP1 OP2)"
CORE_SEL = "not resid 63:88 229:245 284:295 332:335 364:402"
MASK_RADIUS = 6.0  # Å around the RNA atoms
SHELL_EXCLUDE = 1.8  # Å: voxels this close to an RNA atom are not "shell"
TEMPFACTORS = np.linspace(0, 100, 11)  # Å², B-factor scan
MOTIF_RADIUS = 6.0  # Å: voxels compared around the atoms of a motif
BLUR_MARGIN = 5.0  # Å: > 4σ of the largest B-factor of the scan (σ² = B / 8π²)
HALF_SEED = 0  # seed of the half-map weights, shared by all simulations (paired draws)


def family_dir(family):
    return PROJECT_DIR / "density" / FAMILIES[family]


def blocks_path(family, system):
    """Block means of `system` (written by density/block_maps.py)."""
    return PROJECT_DIR / "density" / "blocks" / FAMILIES[family] / system.replace("/", "_")


def load_blocks(family, system):
    """`GridEnsemble` with one frame per continuous trajectory, weighted by its number of frames."""
    return read_ensemble(blocks_path(family, system))


def motif_maps_path(family, system):
    """Per-motif maps of `system` (written by density/motif_maps.py)."""
    return PROJECT_DIR / "density" / "motifs" / FAMILIES[family] / system.replace("/", "_")


def load_motif_maps(family, system, motifs=None):
    """{motif: `GridEnsemble` of its block maps where it is occupied} for the `motifs`
    (default: all) found in `motifs.parquet`."""
    path = motif_maps_path(family, system)
    found = pd.read_parquet(path / "motifs.parquet").index
    motifs = found if motifs is None else found.intersection(motifs)
    return {motif: read_ensemble(path / motif) for motif in motifs}


def group_means(ensemble, labels, chunk_frames=64):
    """Mean map of the frames of each label (sorted unique labels), and their frame counts.

    Reads the (possibly memory-mapped) frames in chunks; frame weights are ignored.
    """
    groups, inverse = np.unique(np.asarray(labels), return_inverse=True)
    data = ensemble.data
    sums = np.zeros((len(groups), data.shape[1]))
    for start in range(0, len(data), chunk_frames):
        chunk = np.asarray(data[start : start + chunk_frames])
        # Sum of each group's frames in the chunk, as a product with its one-hot matrix
        onehot = (inverse[start : start + chunk_frames, None] == np.arange(len(groups))).astype(
            chunk.dtype
        )
        sums += onehot.T @ chunk
    counts = np.bincount(inverse, minlength=len(groups))
    return groups, sums / counts[:, None], counts


def half_maps():
    """The two Con2 half maps of EMD-42499."""
    return [read_mrc(PROJECT_DIR / "density" / f"emd_42499_Con2-2.2A_half_{h}.mrc") for h in "AB"]


def reference_maps():
    """Cryo-EM maps of EMD-42499: mean of the Con2 half maps ("half"), and the sharpened map."""
    half = half_maps()
    return {
        "half": 0.5 * (half[0] + half[1]),
        "sharpened": read_mrc(PROJECT_DIR / "density" / "emd_42499.mrc"),
    }


@functools.cache
def reference_rna():
    """RNA atoms of 9CBU (`RNA_SEL`)."""
    with warnings.catch_warnings():  # placeholder unit cell in the PDB file, not used here
        warnings.filterwarnings("ignore", message=".*CRYST1 record")
        u = mda.Universe(REF_PDB, to_guess=("elements",))
    return u.select_atoms(RNA_SEL)


def masks(params):
    """Dense boolean masks around the 9CBU RNA: all, core, core phosphates, and their shell only."""
    rna = reference_rna()
    core = rna.select_atoms(CORE_SEL)
    params = params.dense()
    result = {
        "All": within(params, rna.positions, MASK_RADIUS),
        "Core": within(params, core.positions, MASK_RADIUS),
        "Core Phosphates": within(params, core.select_atoms("name P").positions, MASK_RADIUS),
    }
    shell = result["All"] & ~within(params, rna.positions, SHELL_EXCLUDE)
    for name in list(result):
        result[f"{name} (Shell only)"] = result[name] & shell
    return result


def motif_coords(motif):
    """Positions in 9CBU of the inner-sphere ligands of `motif` (a row of a motif table),
    or None if one of them is not among the reference atoms."""
    rna = reference_rna()
    coords = []
    for resname, resid, name in zip(
        motif["inner_resnames"], motif["inner_resids"], motif["inner_names"]
    ):
        atom = rna.select_atoms(f"resname {resname} and resid {resid} and name {name}")
        if len(atom) != 1:
            return None
        coords.append(atom.positions[0])
    return np.array(coords)


def motif_voxels(params, coords):
    """Flat indices of the voxels of sparse `params` within `MOTIF_RADIUS + BLUR_MARGIN` of
    the motif atoms `coords`."""
    mask = within(params.dense(), coords, MOTIF_RADIUS + BLUR_MARGIN).ravel()
    return params.voxel_indices[mask[params.voxel_indices]]


def motif_values(grid, coords, tempfactor=0.0):
    """Grid or ensemble stored on the `motif_voxels` of the motif atoms `coords`, blurred,
    at its stored voxels within `MOTIF_RADIUS` of them (sparse, in the box around them).
    Unstored voxels count as zero."""
    if tempfactor and 4 * np.sqrt(tempfactor / (8 * np.pi**2)) > BLUR_MARGIN:
        raise ValueError(f"B = {tempfactor} Å² reaches beyond BLUR_MARGIN = {BLUR_MARGIN} Å")
    trimmed = grid.trim(BLUR_MARGIN)  # zero padding: the FFT is periodic
    box = trimmed.to_dense()
    if tempfactor:
        box = box.gaussian_filter(tempfactor=tempfactor)
    target = within(box.params, coords, MOTIF_RADIUS).ravel()
    target &= np.isin(np.arange(box.params.size), trimmed.params.voxel_indices)
    return box.to_sparse(mask=target)


def blur(grid, tempfactor):
    """Dense numpy copy of `grid` blurred with a B-factor (Å²), computed on the GPU if available."""
    grid = grid.to_dense()
    if not tempfactor:
        return grid
    if torch.cuda.is_available():
        return (
            grid.to(backend="torch", device="cuda")
            .gaussian_filter(tempfactor=tempfactor)
            .to(backend="numpy")
        )
    return grid.gaussian_filter(tempfactor=tempfactor)


def blurred_blocks(blocks, tempfactor, mask):
    """Values (n_blocks, n_mask_voxels) of each block map, blurred, within `mask`."""
    return np.stack(
        [
            np.asarray(blur(Grid(blocks.data[i], blocks.params), tempfactor).array)[mask]
            for i in range(blocks.n_frames)
        ]
    )


def correlation(x, y):
    return np.corrcoef(x, y)[0, 1]


def tempfactor_scan(grid, refs, masks, tempfactors=TEMPFACTORS):
    """Correlation of `grid` blurred with each B-factor, with each reference map within each
    mask: {(reference name, mask name): array over `tempfactors`}."""
    refs = {name: np.asarray(ref.to_dense().array) for name, ref in refs.items()}
    cc = {(r, m): [] for r in refs for m in masks}
    for tempfactor in tempfactors:
        values = np.asarray(blur(grid, tempfactor).array)
        for r, ref in refs.items():
            for m, mask in masks.items():
                cc[r, m].append(correlation(values[mask], ref[mask]))
    return {key: np.array(v) for key, v in cc.items()}


def best_tempfactor(tempfactors, cc, resolution=0.1):
    """B-factor maximizing a cubic interpolation of the correlation scan."""
    fine = np.arange(tempfactors[0], tempfactors[-1] + resolution / 2, resolution)
    return fine[np.argmax(interp1d(tempfactors, cc, kind="cubic")(fine))]


# Correlations with the experiment, bootstrapped over the blocks and the two half maps.
# With the centred rows S (block maps, then the two half maps) and their Gram matrix
# G = S Sᵀ, a map Σ c_i S_i with Σ c_i = 1 is centred, and its correlation with the
# reference Σ d_i S_i is cᵀGd / √(cᵀGc · dᵀGd): exact, without forming the maps.


def _gram(*rows):
    s = np.concatenate([np.asarray(r, dtype=np.float64) for r in rows])
    s -= s.mean(axis=1, keepdims=True)
    return s @ s.T


def _cc(gram, c, d):
    # Correlations of the maps of coefficients `c` (..., n) and `d` (..., n)
    gc, gd = c @ gram, d @ gram
    return (gc * d).sum(-1) / np.sqrt((gc * c).sum(-1) * (gd * d).sum(-1))


def _draws(n_blocks, n_resamples, rng):
    # Block weights (Dirichlet(1), first row: none) and half-map weights (first row: ½, ½)
    blocks = np.vstack(
        [np.ones(n_blocks), np.random.default_rng(rng).dirichlet(np.ones(n_blocks), n_resamples)]
    )
    halves = np.vstack(
        [[0.5, 0.5], np.random.default_rng(HALF_SEED).dirichlet(np.ones(2), n_resamples)]
    )
    return blocks, halves


def _result(values, confidence_level):
    # BootstrapResult from the statistic on the data (first row) and on the resamples
    q = (1 - confidence_level) / 2
    low, high = np.quantile(values[1:], (q, 1 - q), axis=0)
    return BootstrapResult(
        ConfidenceInterval(low, high), values[1:], np.std(values[1:], ddof=1, axis=0), values[0]
    )


def cc_bootstrap(values, weights, half_values, n_resamples=200, confidence_level=0.95, rng=None):
    """Correlation of the mean of block `values` (n_blocks, n_voxels), weighted by `weights`,
    with the mean of the two `half_values` (2, n_voxels). Bayesian bootstrap: Dirichlet(1)
    weights on the blocks (times `weights`) and on the half maps."""
    weights = np.asarray(weights, dtype=np.float64)
    gram = _gram(values, half_values)
    w, h = _draws(len(weights), n_resamples, rng)
    c = w * weights
    c /= c.sum(axis=1, keepdims=True)
    d = np.zeros((len(c), len(gram)))
    c = np.hstack([c, np.zeros((len(c), 2))])
    d[:, -2:] = h
    return _result(_cc(gram, c, d), confidence_level)


def mixture_cc_bootstrap(
    present,
    n_present,
    absent,
    n_absent,
    half_values,
    fractions,
    n_resamples=200,
    confidence_level=0.95,
    rng=None,
):
    """Correlation of f·P + (1 - f)·A with the mean of the two `half_values` (2, n_voxels),
    for each f of `fractions`. P and A are the means of the block maps `present` and `absent`
    (n_blocks, n_voxels) weighted by their frame counts; a block missing one of them has a
    count of 0. Bayesian bootstrap as `cc_bootstrap`, with the same weight for the present
    and absent maps of a block. Returns the `BootstrapResult` of the correlations (per f)
    and of the best f."""
    n_present, n_absent = (np.asarray(n, dtype=np.float64) for n in (n_present, n_absent))
    gram = _gram(present, absent, half_values)
    w, h = _draws(len(n_present), n_resamples, rng)
    p, a = w * n_present, w * n_absent
    p /= p.sum(axis=1, keepdims=True)
    a /= a.sum(axis=1, keepdims=True)
    f = np.asarray(fractions)[:, None, None]
    c = np.concatenate(
        [f * p, (1 - f) * a, np.zeros((*np.broadcast_shapes(f.shape, p.shape)[:2], 2))], axis=-1
    )
    d = np.zeros(c.shape[1:])
    d[:, -2:] = h
    cc = _cc(gram, c, d).T  # (1 + n_resamples, n_fractions)
    best = np.asarray(fractions)[np.argmax(cc, axis=1)]
    return _result(cc, confidence_level), _result(best, confidence_level)
