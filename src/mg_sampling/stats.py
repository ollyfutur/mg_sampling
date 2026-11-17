"""Histogram helpers and bootstrap error estimates."""

from collections import namedtuple
from dataclasses import dataclass
import inspect

import numpy as np


def finitemax(x):
    x = np.atleast_1d(x)
    return np.max(x[np.isfinite(x)])


def finitemin(x):
    x = np.atleast_1d(x)
    return np.min(x[np.isfinite(x)])


def edges_to_range(edges):
    return 0.5 * (edges[1:] + edges[:-1])


ConfidenceInterval = namedtuple("ConfidenceInterval", ["low", "high"])


@dataclass
class BootstrapResult:
    """Result object returned by `scipy.stats.bootstrap`.

    Attributes
    ----------
    confidence_interval : ConfidenceInterval
        The bootstrap confidence interval as an instance of
        `collections.namedtuple` with attributes `low` and `high`.
    bootstrap_distribution : ndarray
        The bootstrap distribution, that is, the value of `statistic` for
        each resample. The last dimension corresponds with the resamples
        (e.g. ``res.bootstrap_distribution.shape[-1] == n_resamples``).
    standard_error : float or ndarray
        The bootstrap standard error, that is, the sample standard
        deviation of the bootstrap distribution.
    estimate : float or ndarray
        The statistic on the original data (`block_bootstrap` only).

    """

    confidence_interval: ConfidenceInterval
    bootstrap_distribution: np.ndarray
    standard_error: float | np.ndarray
    estimate: float | np.ndarray | None = None


def bayesian_bootstrap(
    data,
    statistic,
    *,
    alpha=1,
    n_resamples=200,
    rng=None,
    confidence_level=0.95,
    weights=None,
    **kwargs,
):
    """
    Perform bayesian bootstrap on a sequence of observations `data` for a given `statistic`.
    A weight vector for each observation of `data` is drawn from a Dirichlet distribution
    of dimension `len(data)` with parameter `alpha`, and combined with the provided `weights`
    if any. If `statistic` does not accept the `weights` keyword argument, two-level sampling is
    performed.
    """
    rng = np.random.default_rng(rng)
    accepts_weights = "weights" in inspect.signature(statistic).parameters
    q_l, q_u = (1 - confidence_level) / 2, 1 - (1 - confidence_level) / 2

    n_samples = len(data)

    # If no weights passed, use uniform Dirichlet
    if weights is None:
        weights = np.ones(n_samples)
    else:
        weights = np.asarray(weights)

    # Normalize weights to mean concentration
    weights = weights / weights.mean() * alpha

    weights = rng.dirichlet(weights, n_resamples)

    if accepts_weights:
        res = [statistic(data, weights=weights[i], **kwargs) for i in range(n_resamples)]
    else:
        res = []
        for i in range(n_resamples):
            _data = [
                data[j]
                for j in np.random.choice(n_samples, size=n_samples, p=weights[i])
            ]
            res.append(statistic(_data, **kwargs))

    res = np.asarray(res)

    ci_l, ci_u = np.quantile(res, (q_l, q_u))

    return BootstrapResult(
        confidence_interval=ConfidenceInterval(ci_l, ci_u),
        bootstrap_distribution=res,
        standard_error=np.std(res, ddof=1, axis=-1),
    )


def bootstrap(data, statistic, weights=None, n_resamples=200, rng=None, confidence_level=0.95):
    """
    Classical bootstrap with weighted samples.
    """
    rng = np.random.default_rng(rng)
    q_l, q_u = (1 - confidence_level) / 2, 1 - (1 - confidence_level) / 2
    n_samples = len(data)

    if weights is not None:
        weights = np.asarray(weights)
        weights = weights / weights.sum()

    res = []
    for _ in range(n_resamples):
        inds = rng.choice(n_samples, size=n_samples, replace=True, p=weights)
        _data = data[inds]
        res.append(statistic(_data))

    ci_l, ci_u = np.quantile(res, (q_l, q_u))

    return BootstrapResult(
        confidence_interval=ConfidenceInterval(ci_l, ci_u),
        bootstrap_distribution=res,
        standard_error=np.std(res, ddof=1, axis=-1),
    )


def block_bootstrap(data, statistic, *, block_id=None, n_blocks=None, weights=None,
                    bayesian=False, n_resamples=200, rng=None, confidence_level=0.95):
    """
    Bootstrap `statistic` over blocks of rows of `data` (for correlated samples).

    Blocks are the unique values of `block_id` (one label per row), or `n_blocks` contiguous
    blocks; by default each row is a block. Classical: blocks are drawn with replacement.
    Bayesian: blocks get Dirichlet(1) weights, multiplied with the row `weights` and passed as
    `statistic(data, weights=...)`. Vector statistics get one interval per element (axis 0).
    """
    rng = np.random.default_rng(rng)
    q_l, q_u = (1 - confidence_level) / 2, 1 - (1 - confidence_level) / 2
    data = np.asarray(data)
    n_samples = len(data)
    if weights is not None:
        weights = np.asarray(weights, dtype=float)

    if block_id is not None:
        _, row_block = np.unique(np.asarray(block_id), return_inverse=True)
        n_blocks = row_block.max() + 1
        block_inds = [np.flatnonzero(row_block == b) for b in range(n_blocks)]
    else:
        block_inds = np.array_split(np.arange(n_samples), n_blocks or n_samples)
        n_blocks = len(block_inds)
        row_block = np.repeat(np.arange(n_blocks), [len(b) for b in block_inds])

    def call(x, w):
        return statistic(x) if w is None else statistic(x, weights=w)

    res = []
    for _ in range(n_resamples):
        if bayesian:
            w = rng.dirichlet(np.ones(n_blocks))[row_block]
            res.append(call(data, w if weights is None else w * weights))
        else:
            inds = np.concatenate([block_inds[i] for i in rng.choice(n_blocks, size=n_blocks, replace=True)])
            res.append(call(data[inds], None if weights is None else weights[inds]))
    res = np.asarray(res)

    ci_l, ci_u = np.quantile(res, (q_l, q_u), axis=0)

    return BootstrapResult(
        confidence_interval=ConfidenceInterval(ci_l, ci_u),
        bootstrap_distribution=res,
        standard_error=np.std(res, ddof=1, axis=0),
        estimate=call(data, weights),
    )
