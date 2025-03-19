from functools import partial

import numpy as np

from scipy.stats import entropy
from scipy.optimize import minimize

from mg_sampling.tools import edges_to_range
from mg_sampling.classes import MetadSim
from mg_sampling.constants import RT, π

PARAMS_NAMES = ('E', 'sX', 'sY')
ASSOCIATIVE_P_0 = (1, 0.1, 0.2)
ASSOCIATIVE_BOUNDS = ((0, np.inf), (-1/π, 1/π), (0, 1/2))
DISSOCIATIVE_P_0 = (1, 0.1, -0.2)
DISSOCIATIVE_BOUNDS = ((0, np.inf), (-1/π, 1/π), (-1/π, 0))

def bias_fn(x, y, E, sX, sY):
    """
    Bias function based on a cosine envelope for the associative pathway.
    """
    X = x - y + 6
    Y = x + y - 6
    _X = X + sX * (1 - np.cos(X*π))
    _Y = Y - sY * (1 - np.cos(X*π))

    return -E * (_X > 0) * (_Y > 0) * (_Y < 1/2) * (1 - np.cos(_X*π)) * (1 - np.cos(_Y*π*4)) / 4

# def associative_loss_max_min(sim: MetadSim, bias_fn, params):
#     bias_weights = np.exp(-bias_fn(sim.X, sim.Y, params)/RT) * sim.mask
#     ρX, _ = np.histogram(sim.X, bins=X_EDGES, density=True, weights=sim.weights*bias_weights)
#     FX = -RT*np.log(ρX)
#     # Remove trend 
#     FX -= X_RANGE * (FX[-1] - FX[0])/(X_RANGE[-1] - X_RANGE[0]) + FX[0]
#     return (np.max(FX) - np.min(FX))/RT

# def associative_loss_max_min_kl_p(sim: MetadSim, bias_fn, params):
#     bias_weights = np.exp(-bias_fn(sim.X, sim.Y, params)/RT) * sim.mask
#     ρX, _ = np.histogram(sim.X, bins=X_EDGES, density=True, weights=sim.weights*bias_weights)
#     FX = -RT*np.log(ρX)
#     # Remove trend 
#     FX -= X_RANGE * (FX[-1] - FX[0])/(X_RANGE[-1] - X_RANGE[0]) + FX[0]
#     new_weights = sim.weights * bias_weights
#     rel_entropy = entropy(sim.weights * sim.mask, new_weights)
#     return (np.max(FX) - np.min(FX))/RT + rel_entropy

def associative_mask(x, y):
    X = x - y + 6
    Y = x + y - 6
    return  (X > 0) * (X < 2) * (Y > 0) * (Y < 2)

def dissociative_mask(x, y):
    X = x - y + 6
    Y = x + y - 6
    return  (X > 0) * (X < 2) * (Y < 2) * (x > 0) * (y > 5)

def associative_loss_fn(sim: MetadSim, bias_fn, params):
    """
    Loss for flattening the associative pathway.
    - FES restricted to the (0 < X < 2, 0 < Y < 2) square is projected on X
    - Trend between X=0.1 and X=1.9 is substracted
    - Loss targets barrier height (max - min)
    - Regularized with relative entropy DKL(new_weights, old_weights)
    """
    X_edges = np.linspace(0.1, 1.9, 2*20+1)
    X_range = edges_to_range(X_edges)
    mask = associative_mask(sim.x, sim.y)
    bias_weights = np.exp(-bias_fn(sim.x, sim.y, *params)/RT) * mask
    ρX, _ = np.histogram(sim.X, bins=X_edges, density=True, weights=sim.weights*bias_weights)
    FX = -RT*np.log(ρX)
    # Remove trend 
    FX -= X_range * (FX[-1] - FX[0])/(X_range[-1] - X_range[0]) + FX[0]
    new_weights = sim.weights * bias_weights
    rel_entropy = entropy(new_weights, sim.weights * mask)
    return (np.max(FX) - np.min(FX))/RT + rel_entropy

def dissociative_loss_fn(sim: MetadSim, bias_fn, params):
    """
    Loss for flattening the dissociative pathway.
    - FES restricted to the (0 < X < 2) and (x > 0) and (y > 5) is projected on X
    - Trend between X=0.1 and X=1.9 is substracted
    - Loss targets barrier height (max - min)
    - Regularized with relative entropy DKL(new_weights, old_weights)
    """
    X_edges = np.linspace(0.1, 1.9, 2*20+1)
    X_range = edges_to_range(X_edges)
    mask = dissociative_mask(sim.x, sim.y)
    bias_weights = np.exp(-bias_fn(sim.x, sim.y, *params)/RT) * mask
    ρX, _ = np.histogram(sim.X, bins=X_edges, density=True, weights=sim.weights*bias_weights)
    FX = -RT*np.log(ρX)
    # Remove trend 
    FX -= X_range * (FX[-1] - FX[0])/(X_range[-1] - X_range[0]) + FX[0]
    new_weights = sim.weights * bias_weights
    rel_entropy = entropy(new_weights, sim.weights * mask)
    return (np.max(FX) - np.min(FX))/RT


def fit_bias(sim: MetadSim, bias_fn, loss_fn, p0, bounds=None):
    _loss = partial(loss_fn, sim, bias_fn)
    res = minimize(_loss, p0, method='L-BFGS-B', bounds=bounds)
    print(res)
    return res.x
