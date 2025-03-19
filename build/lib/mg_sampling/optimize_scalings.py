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
    w_nk = [_get_w_n(logw_n, u_n, λ) for λ in λ_k]
    acceptances = [_get_p_acc(w1, w2) for w1, w2 in zip(w_nk[:-1], w_nk[1:])]
    return np.array(acceptances)

def _equal_acc_loss(acceptances):
    return np.var(acceptances)

_losses = {
    'equal': _equal_acc_loss,
}

def optimize_scalings(logw_n: ArrayLike, u_n: ArrayLike, nrep: int, guess: ArrayLike = None, target: str = 'equal', **kwargs) -> np.ndarray:
    """Finds scalings for uniform acceptance given a series of configurational log-weights `logw_n` in the reference (λ=0) state,
    the corresponding reduced bias (at λ=1), and the number of desired replicas `nrep` (including λ=0 and λ=1).
    If a `guess`is given `nrep` is ignored and the upper and lower scalings are kept.
    """
    logw_n = np.asarray(logw_n) - logsumexp(logw_n)
    u_n = np.asarray(u_n)

    if guess is None:
        guess = np.linspace(0, 1, nrep)
    else:
        guess = np.asarray(guess)

    λ0, x0, λ1 = guess[0], guess[1:-1], guess[-1]

    if target not in _losses:
        raise ValueError("Target must be one of:", ', '.join(_losses.keys()))
    _loss = _losses[target]

    def loss(x):
        λ_k = np.concatenate([[λ0], x, [λ1]])
        acceptances = _get_acceptances(logw_n, u_n, λ_k)
        return _loss(acceptances)

    res = minimize(loss, x0, **kwargs)

    if res.success:
        λ_k_opt = np.concatenate([[λ0], res.x, [λ1]])
        return λ_k_opt
    else:
        raise RuntimeError(f"Failed optimization: {res}")

def main(args):
    RT = R * args.temp * 1e-3

    λ_init = np.loadtxt(args.scalings)

    nrep = args.nrep

    U_kN = []
    for bias_file in args.bias_files:
        df = read_plumed(bias_file)
        df = df.loc[(df['time'] >= args.start) & (df['time'] <= args.stop)]
        U_kN.append(df[args.bias_col].to_numpy())

    U_n = np.concatenate(U_kN)
    u_n = U_n / RT
    u_nk = λ_init * u_n[:, None]
    logw_n = wham(u_nk, T=1.0).logW

    print("Initial scalings:", ' '.join(λ_init))
    print("Initial predicted acceptances:", ' '.join(_get_acceptances(logw_n, u_n, λ_init)))

    guess = λ_init[1:-1]
    λ_opt = optimize_scalings(logw_n, u_n, nrep, guess)

    print("Successfully optimized scalings.")
    print("Optimized scalings:", ' '.join(λ_opt))
    print("Predicted acceptances:", ' '.join(_get_acceptances(logw_n, u_n, λ_opt)))

    np.savetxt(args.outfile, λ_opt)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bias_files', type=str, nargs='+',
                        help="Plumed output files for each replica containing the instantaneous values of the bias.")
    parser.add_argument('--scalings_file', type=str, help="Text file with initial scalings (column format)")
    parser.add_argument('--bias_col', type=str, default='bias.bias',
                        help="(bias.bias) Column in the Plumed files corresponding to the applied bias for each replica.")
    parser.add_argument('-o', '--outfile', type=str, default='opt_scalings.txt')
    parser.add_argument('--temp', type=float, default=300, help="(300) Temperature [K]")
    parser.add_argument('--start', type=float, default=0.0, help="(0) Starting time considered for optimization [ps]")
    parser.add_argument('--stop', type=float, default=float('inf'), help="(inf) Ending time [ps]")

    args = parser.parse_args()

    main(args)

