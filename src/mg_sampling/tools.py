import logging

import numpy as np
import pandas as pd
import re
from copy import copy

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
    standard error of the mean on the `n_blocks`, ignoring the first `skip` blocks.
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

def bootstrap(*inputs, func, n_resample=20, weights=None):
    n_frames = len(inputs[0])

    outputs = []
    for _ in range(n_resample):
        inds = np.random.choice(n_frames, size=n_frames, p=weights)
        _inputs = (inp[inds] for inp in inputs)
        outputs.append(func(*_inputs))
    
    mean = np.mean(outputs, axis=0)
    stderr = np.std(outputs, axis=0)

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

def read_xvg(path):
    """
    Reads .xvg as pandas dataframe.
    """
    columns = []
    with open(path, 'r') as f:
        for line in f:
            x = re.search('^@    xaxis  label \"(.*)\"', line)
            if x is not None:
                columns.append(x.group(1))
                break
        else:
            columns.append("")
        
        for line in f:
            x = re.search('^@ s\d+ legend \"(.*)\"', line)
            if x is not None:
                columns.append(x.group(1))
    
    data = np.loadtxt(path, comments=['#', '@'])
    df = pd.DataFrame(data, columns=columns).set_index(columns[0])
    return df

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

def demux(data_list, data_times=None, repindex_file='replica_index.xvg', rtol=1e-8, atol=1e-6):
    """
    From an array-like (n_rep, n_frames, n_features)
    outputs demuxed trajectories (n_rep, n_frames_processed, n_features)
    using the replice index file `repindex_file`.

    If `data_times` is `None` the first feature of the first replica is used.

    The parameter `rtol` and `atol` are the relative and absolute tolerances to consider `data_times`
    and times in `repindex_file` equal.
    """
    repindex = np.loadtxt(repindex_file)
    repindex_times, repindex = repindex[:, 0], repindex[:, 1:]
    unique, unique_indices, counts = np.unique(repindex_times, return_index=True, return_counts=True)
    if np.any(counts > 1):
        logging.warning(f"Found repeated times in {repindex_file}. Duplicates will be dropped, make sure they are consistent!")
    repindex_times = unique
    repindex = repindex[unique_indices]
    n_reps = repindex.shape[1]

    data_arr = np.array(data_list)

    if data_times is None:
        data_times = data_arr[0, :, 0]

    dt = min(np.diff(repindex_times).mean(), np.diff(data_times).mean())

    if data_arr.shape[0] != n_reps:
        raise ValueError(f"Number of replicas in {repindex_file} doesn't match `data_list`")

    # Gromacs saves then does replica exchange, hence at t=1*dt
    # the replicas have the configuration strictly before t in the
    # replica_index file, except at t=0 where they are exactly at t.

    max_time = np.max(repindex_times)

    demuxed_data = []
    for itime, time in enumerate(data_times):
        if time > max_time:
            logging.warning(f"Data times go above replica_index times. Data will be truncated after {max_time}.")
            break
        idx = np.argmin(repindex_times - time)
        if np.abs(repindex_times[idx] - time) < rtol * time + atol:  # Case 1: matching `time`
            idx = max(idx - 1, 0)  # Take one step before `time` in `replica_index` (or 0)
        else:  # Case 2: not matching `time`
            idx = np.argwhere(repindex_times < time)
            idx = idx[-1][0]  # Take first preceding `time`
        reptemp = np.argsort(repindex[idx]).astype(int)
        demuxed_data.append(data_arr[reptemp, itime, :])

    return np.array(demuxed_data).swapaxes(0, 1)

def read_performance(logfile):
    """
    Gets performance (hour/ns) from Gromacs log file.
    """
    with open(logfile) as f:
        for line in f:
            if line.startswith('Performance:'):
                return float(line.split()[2])
        else:
            raise RuntimeError(f"Could not find performance in {logfile}")