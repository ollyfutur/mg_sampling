import numpy as np
from scipy.constants import R
from scipy.special import logsumexp
from scipy.optimize import minimize

import warnings

import plumed
from bussilab.wham import wham

def w_nk(λ_k):
    logw_nk = W.logW[:, None] - λ_k * U_n[:, None] / RT
    logw_nk -= logsumexp(logw_nk, axis=0)
    return np.exp(logw_nk)

def p_acc(w_nk, i, j, stride=1):
    w_i = w_nk[::stride, i]
    w_j = w_nk[::stride, j]
    w_i /= w_i.sum()
    w_j /= w_j.sum()
    ww = w_i[:, None] * w_j[None, :]
    return np.sum(np.minimum(ww, ww.T))

def read_p_acc(path):
    with open(path, 'r') as f:
        for line in f:
            if line.startswith("Repl  average probabilities:"):
                next(f)
                p_acc = np.array(next(f).split()[1:], dtype=float)
                return p_acc
    raise ValueError("Average probabilities not found in log file.")

def loss(scalings):
    _scalings = np.concatenate([[0], scalings, [1]])
    _w_nk = w_nk(_scalings)
    _p_acc = [p_acc(_w_nk, k, k+1, stride=2) for k in range(len(_scalings)-1)]
    return -np.prod(_p_acc)

def loss(scalings):
    _scalings = np.concatenate([[0], scalings, [1]])
    _w_nk = w_nk(_scalings)
    _p_acc = np.array([p_acc(_w_nk, k, k+1, stride=2) for k in range(len(_scalings)-1)])
    return -np.log(_p_acc).sum()

RT = R * 300 * 1e-3
scalings = np.loadtxt(f'scalings.txt')
empirical_p_acc = read_p_acc('rep00/run.log')

nrep = len(scalings)

with warnings.catch_warnings(action='ignore'):
    U_kN = [
        plumed.read_as_pandas(f'rep{k:02d}/bias.COLVAR.{k}')['cv.tot']
        for k in range(nrep)
    ]

U_n = np.concatenate(U_kN)
U_nk = scalings * U_n[:, None]
W = wham(U_nk, T=RT)

bounds = (nrep - 2) * [(0, 1)]
res = minimize(loss, scalings[1:-1], bounds=bounds, method='Nelder-Mead')
opt_scalings = np.concatenate([[0], res.x, [1]])
opt_p_acc = np.array([p_acc(w_nk(opt_scalings), k, k+1) for k in range(nrep-1)])
np.savetxt('opt_scalings.txt', opt_scalings)

print("Scalings:", scalings)
print("Acceptances:", empirical_p_acc)
if res.success:
    print("Successfully optimized scalings.")
    print("Optimized scalings:", opt_scalings)
    print("Predicted acceptances:", opt_p_acc)
