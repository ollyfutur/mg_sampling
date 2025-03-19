import numpy as np
import pandas as pd

def finitemax(x):
    x = np.atleast_1d(x)
    return np.max(x[np.isfinite(x)])

def finitemin(x):
    x = np.atleast_1d(x)
    return np.min(x[np.isfinite(x)])

def edges_to_range(edges):
    return 0.5*(edges[1:]+edges[:-1])

def range_to_edges(arange):
    dx = arange[1] - arange[0]
    return np.linspace(arange[0] - dx/2, arange[-1] + dx/2, len(arange) + 1)

def block_analysis(*inputs, func, n_blocks=5, skip=0):
    """
    Perform a block analysis of the statistic `func` that takes as arguments one or more `inputs`
    timeseries, by splitting each series in `n_blocks` of equal length. Returns the mean and
    standard error of the mean on the `n_blocks`.
    """
    n_frames = len(inputs[0])
    if isinstance(skip, float):
        skip = int(n_frames * skip)
    block_inds = np.array_split(np.arange(skip, n_frames), n_blocks)

    outputs = []
    for inds in block_inds:
        _inputs = (inp[inds] for inp in inputs)
        outputs.append(func(*_inputs))
    
    mean = np.mean(outputs, axis=0)
    var = np.var(outputs, axis=0)
    stderr = np.sqrt(var/n_blocks)

    return mean, stderr

def to_wsv(df: pd.DataFrame, path):
    with open(path, 'w') as f:
        f.write(df.to_string(index=False))

def read_wsv(path):
    df = pd.read_csv(path, delimiter='\s+')
    return df

def read_plumed(path):
    with open(path, 'r') as f:
        header = next(f)
        columns = header.split()[2:]
    df = pd.read_csv(path, delimiter='\s+', comment='#', header=None)
    df.columns = columns
    return df

def read_plumed_fes(path):
    """
    Reads output of `sum_hills` (fes.dat).
    Returns
    -------
    cv_names : tuple of str
    cv_ranges : tuple of 1d ndarrays
    fes : ndarray
    """
    df = read_plumed(path)
    d = (df.shape[1] - 1) // 2
    df = df.set_index(list(df.columns[:d])).sort_index()
    fes_shape = tuple([len(x) for x in df.index.levels])
    fes = np.reshape(df['file.free'], fes_shape)
    cv_names = tuple([idx.name for idx in df.index.levels])
    cv_ranges = tuple([np.array(idx) for idx in df.index.levels])
    return cv_names, cv_ranges, fes


def count_round_trips(x, low=0, high=1):
    """
    Count the number of "round trips" in a series `x`.
    A round trip is defined as going from `x` < `low` to `x` >= `high` back to `low`, or vice versa.
    Note that `low` is excluded while `high` is included.
    """
    state = np.digitize(x, (low, high))
    # Only keep end states (0 and 2)
    state = state[state != 1]
    # Count transitions and divide by 2
    count = np.sum(state[1:] != state[:-1]) // 2
    return count

def detect_round_trips(x, low=0, high=1):
    x[np.isnan(x)] = (low + high) / 2
    state = np.digitize(x, (low, high))
    inds = []
    s0 = None
    for i, s in enumerate(state):
        if s == 0 and s0 == 2:
            inds.append(i)
            s0 = 0
        elif s == 2 and s0 == 0:
            inds.append(i)
            s0 = 2
        elif s in (0, 2) and s0 is None:
            s0 = s
    return inds

def kish_size(weights):
    weights = np.asarray(weights)
    return weights.sum()**2 / (weights**2).sum()

def kish_ratio(weights):
    return kish_size(weights) / len(weights)
