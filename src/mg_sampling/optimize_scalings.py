import argparse

import numpy as np
from numpy.typing import ArrayLike, NDArray

from scipy.constants import R
from scipy.special import logsumexp
from scipy.optimize import minimize

from bussilab.wham import wham
from mg_sampling.tools import read_plumed

def _get_w_n(logw_n: NDArray, u_n: NDArray, λ: float) -> NDArray:
        _logw_n = logw_n - λ * u_n
        _logw_n -= logsumexp(_logw_n)
        return np.exp(_logw_n)

# Broadcasting is slower than loop, probably because of memory access
# def get_p_acc(w1, w2):
#     ww = w1[:, None] * w2[None, :]
#     return np.sum(np.minimum(ww, ww.T))
def _get_p_acc(w1: NDArray, w2:NDArray) -> NDArray:
    p_acc = 0.0
    for i in range(len(w1)):
        p_acc += 2 * np.sum(np.minimum(w1[i] * w2[i+1:], w1[i+1:] * w2[i])) + w1[i] * w2[i]
    return p_acc

def _get_acceptances(logw_n: NDArray, u_n: NDArray, λ_k: NDArray) -> NDArray:
    w_k = [_get_w_n(logw_n, u_n, λ) for λ in λ_k]
    acceptances = [_get_p_acc(w1, w2) for w1, w2 in zip(w_k[:-1], w_k[1:])]
    return np.array(acceptances)

def _uniform_acc_loss(acceptances):
    return np.var(acceptances)

def _max_acc_loss(acceptances):
    return -np.sum(np.log(acceptances))


_losses = {
    'uniform': _uniform_acc_loss,
    'max': _max_acc_loss,
}

def optimize_scalings(
        logw_n: ArrayLike,
        u_n: ArrayLike,
        nrep: int,
        guess: ArrayLike | None = None,
        lmin: float = 0.0,
        lmax: float = 1.0,
        target: str = 'uniform',
        verbose: bool = False,
        tol: float = 1e-3,
        **kwargs,
    ) -> np.ndarray:
    """
    Finds scalings for uniform acceptance given a series of configurational log-weights `logw_n` in the reference (λ=0) state,
    the corresponding unscaled reduced bias (at λ=1), and the number of desired replicas `nrep` (from λ `lmin` to `lmax`).
    """
    logw_n = np.asarray(logw_n) - logsumexp(logw_n)
    u_n = np.asarray(u_n)

    if guess is None:
        guess = np.linspace(lmin, lmax, nrep)
    elif len(guess) != nrep:
        raise ValueError(f"`guess` should have {nrep=} elements, found: {len(guess)}")
        guess = np.asarray(guess)
    else:
        guess = np.asarray(guess)

    x0 = guess[1:-1]

    if target not in _losses:
        raise ValueError("Target must be one of:", ', '.join(_losses.keys()))
    loss_fn = _losses[target]

    def loss(x):
        λ_k = np.concatenate([[lmin], x, [lmax]])
        acceptances = _get_acceptances(logw_n, u_n, λ_k)
        # Regularize: strictly increasing λ
        reg = -np.diff(λ_k)
        reg = np.sum(reg * (reg > 0))
        out = loss_fn(acceptances) + reg
        if verbose: print("Loss:", f"{out:8.6f}", end='\r')
        return out

    res = minimize(loss, x0, tol=tol, **kwargs)
    if verbose: print("\nFinished minimization")

    if res.success:
        λ_k_opt = np.concatenate([[lmin], res.x, [lmax]])
        return λ_k_opt
    else:
        raise RuntimeError(f"Failed optimization: {res}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bias_files', type=str, nargs='+',
                        help="Plumed output files for each replica containing the instantaneous values of the *unscaled* (λ=1) bias.")
    parser.add_argument('--scalings', type=str, help="Text file with initial scalings (column format)")
    parser.add_argument('--nrep', type=int, help="(none) Number of replicas desired if different than `scalings`")
    parser.add_argument('--lmin', type=float, default=0, help="(0) λmin")
    parser.add_argument('--lmax', type=float, default=1, help="(1) λmax")
    parser.add_argument('--bias_col', type=str, nargs='+', default=['bias.bias'],
                        help="(bias.bias) Column(s) in the Plumed files corresponding to the *unscaled* (λ=1) bias(es) for each replica.")
    parser.add_argument('-o', '--outfile', type=str, default='opt_scalings.txt')
    parser.add_argument('--temp', type=float, default=300, help="(300) Temperature [K]")
    parser.add_argument('--target', type=str, choices=_losses.keys(), default='uniform', help="(uniform) Target scaling distribution")
    parser.add_argument('--start', type=float, default=0.0, help="(0) Starting time considered for optimization [ps]")
    parser.add_argument('--stop', type=float, default=float('inf'), help="(inf) Ending time [ps]")
    parser.add_argument('--stride', type=int, default=1, help="(1) Subsample every N frame")
    parser.add_argument('-v', '--verbose', action='store_true', help="Show loss")

    args = parser.parse_args()

    RT = R * args.temp * 1e-3

    λ_init = np.loadtxt(args.scalings)
    nrep = len(λ_init) if args.nrep is None else args.nrep
    guess = λ_init if len(λ_init) == nrep else np.linspace(args.lmin, args.lmax, nrep)

    print(f"Reading unscaled biases:")
    print(f"  columns:", *args.bias_col)
    print(f"  files:", *args.bias_files)
    print(f"  temperature:", args.temp, "K")
    print(f"  start time:", args.start, "ps")
    print(f"  stop time:", args.stop, "ps")
    print(f"  stride:", args.stride, "frames")

    U_kN = []
    for bias_file in args.bias_files:
        df = read_plumed(bias_file)
        df = df.loc[(df['time'] >= args.start) & (df['time'] <= args.stop)]
        df = df.iloc[::args.stride]
        U_N = df[args.bias_col].sum(axis=1).to_numpy()            
        U_kN.append(U_N)

    U_n = np.concatenate(U_kN)
    u_n = U_n / RT
    u_nk = λ_init * u_n[:, None]

    print(f"Analyzing {u_nk.shape[0]} frames × {u_nk.shape[1]} with WHAM")
    print(f"  reduced energies:")
    print("  [0]", *(f"{u:>6.4f}" for u in u_nk[0]))
    print("  …")
    print(f"  [{len(u_nk)}]", *(f"{u:>6.4f}" for u in u_nk[-1]))

    logw_n = wham(u_nk, T=1.0).logW

    print(f"Optimizing bias scalings:")
    print(f"  λmin = {args.lmin}")
    print(f"  λmax = {args.lmax}")
    print(f"  nrep = {nrep}")
    print(f"  Target distribution: {args.target}")
    print(f"  initial guess:", *(f"{λ:>6.2f}" for λ in guess))
    print("Guess predicted acceptances:")
    print("  ", *(f"{acc:>6.2f}" for acc in _get_acceptances(logw_n, u_n, guess)))

    λ_opt = optimize_scalings(
        logw_n,
        u_n,
        nrep,
        guess=guess,
        lmin=args.lmin,
        lmax=args.lmax,
        target=args.target,
        verbose=args.verbose
    )

    np.savetxt(args.outfile, λ_opt)

    print("Successfully optimized scalings.")
    print(f"Initial scalings ({args.scalings}):")
    print(*(f"{λ:>6.2f}" for λ in λ_init))
    print(f"Optimized scalings ({args.outfile}):")
    print(*(f"{λ:>6.2f}" for λ in λ_opt))
    print("Initial predicted acceptances:")
    print("  ", *(f"{acc:>6.2f}" for acc in _get_acceptances(logw_n, u_n, λ_init)))
    print("Optimized predicted acceptances:")
    print("  ", *(f"{acc:>6.2f}" for acc in _get_acceptances(logw_n, u_n, λ_opt)))


if __name__ == '__main__':
    main()

