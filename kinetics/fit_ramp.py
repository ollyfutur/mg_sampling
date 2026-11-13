"""Fit the bias ramp α(λ) and the replica ladder λ_k from constant-scaling kinetics runs."""

import argparse
from pathlib import Path

import numpy as np
from bussilab.wham import wham
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.special import logsumexp

from mg_sampling.constants import RT
from mg_sampling.tools import read_plumed

COEF = 1.039661  # rescales coordination numbers to integer values


def load(path, scalings, stride=10):
    """Bias energies (kT, unscaled), and n_OP, n_OW for all runs."""
    u, nop, now = [], [], []
    for scaling in scalings:
        colvar = read_plumed(path / f"{scaling:.1f}" / "COLVAR").iloc[::stride]
        u.append(np.asarray(colvar["bias.bias"]) / scaling)
        nop.append(np.asarray(colvar["cv.nb-0"] * COEF))
        now.append(np.asarray(colvar["cv.nc-0"] * COEF))
    return np.concatenate(u) / RT, np.concatenate(nop), np.concatenate(now)


def fit_alpha(logw_n, u_n, nop_n, n_lambda=50):
    """α(λ) on a λ grid such that the binding free energy does not depend on λ."""
    is_bound = nop_n[::10] > 0
    u_b = u_n[::10][is_bound]
    nop_b = nop_n[::10][is_bound]
    logw_b = logw_n[::10][is_bound]
    w_b = np.exp(logw_b - logsumexp(logw_b))

    def loss(α_k):
        Δg_k = [-logsumexp(-(u_b * λ + nop_b * α), b=w_b) for λ, α in zip(λ_k, α_k)]
        return np.square(Δg_k).sum()

    λ_k = np.linspace(0, 1, n_lambda)
    res = minimize(loss, x0=np.zeros_like(λ_k))
    return interp1d(λ_k, res.x, kind="cubic", fill_value="extrapolate")


def get_w_n(logw_n, u_n, λ, nop_n, α):
    _logw_n = logw_n - λ * u_n - α * nop_n
    return np.exp(_logw_n - logsumexp(_logw_n))


def get_p_acc(w1, w2):
    p_acc = 0.0
    for i in range(len(w1)):
        p_acc += 2 * np.sum(np.minimum(w1[i] * w2[i + 1:], w1[i + 1:] * w2[i])) + w1[i] * w2[i]
    return p_acc


def get_acceptances(logw_n, u_n, λ_k, nop_n, α_k):
    w_k = [get_w_n(logw_n, u_n, λ, nop_n, α) for λ, α in zip(λ_k, α_k)]
    return np.array([get_p_acc(w1, w2) for w1, w2 in zip(w_k[:-1], w_k[1:])])


def resample(logw_n, u_n, nop_n, α, λ_k, M, rng, size=200):
    """Samples of M-fold systems at each λ, reweighted with WHAM."""
    uM_n, nopM_n = [], []
    for λ in λ_k:
        p = get_w_n(logw_n, u_n, λ, nop_n, α(λ))
        idx_m = [rng.choice(len(u_n), size=size, p=p) for _ in range(M)]
        uM_n.append(np.sum([u_n[idx] for idx in idx_m], axis=0))
        nopM_n.append(np.sum([nop_n[idx] for idx in idx_m], axis=0))
    uM_n, nopM_n = np.concatenate(uM_n), np.concatenate(nopM_n)
    logwM_n = wham(uM_n[:, None] * λ_k + nopM_n[:, None] * α(λ_k), T=1.0).logW
    return logwM_n, uM_n, nopM_n


def fit_ladder(logwM_n, uM_n, nopM_n, α, nrep):
    """λ_k (λ_0 = 0, λ_{nrep-1} = 1) with uniform exchange acceptances."""

    def loss(x):
        λ_k = np.concatenate([[0], x, [1]])
        acceptances = get_acceptances(logwM_n, uM_n, λ_k, nopM_n, α(λ_k))
        reg = -np.diff(λ_k)
        return np.var(acceptances) + np.sum(reg * (reg > 0))

    res = minimize(loss, np.linspace(0, 1, nrep)[1:-1], tol=1e-5)
    return np.concatenate([[0], res.x, [1]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="Kinetics folder (e.g. op_mMg_tip3p)")
    parser.add_argument("--scalings", type=float, nargs="+", default=[0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
    parser.add_argument("--nrep", type=int, default=12, help="Number of replicas")
    parser.add_argument("-M", type=int, default=150, help="Number of independent systems for resampling")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-o", "--output-dir", type=Path, help="Output folder (default: path)")
    args = parser.parse_args()

    output_dir = args.output_dir or args.path
    out_scalings = output_dir / "opt_scalings_ramp.txt"
    out_alphas = output_dir / "opt_alpha_ramp.txt"
    for f in (out_scalings, out_alphas):
        if f.exists():
            raise SystemExit(f"{f} exists, not overwriting.")

    scalings = np.array(args.scalings)
    u_n, nop_n, _ = load(args.path, scalings)
    logw_n = wham(u_n[:, None] * scalings, T=1.0).logW

    α = fit_alpha(logw_n, u_n, nop_n)
    rng = np.random.default_rng(args.seed)
    logwM_n, uM_n, nopM_n = resample(logw_n, u_n, nop_n, α, np.linspace(0, 1, args.nrep), args.M, rng)
    λ_k = fit_ladder(logwM_n, uM_n, nopM_n, α, args.nrep)
    print("Acceptances:", np.round(get_acceptances(logwM_n, uM_n, λ_k, nopM_n, α(λ_k)), 3))

    np.savetxt(out_scalings, λ_k, fmt="%.6f")
    np.savetxt(out_alphas, RT * np.abs(α(λ_k)), fmt="%.6f")  # kJ/mol


if __name__ == "__main__":
    main()
