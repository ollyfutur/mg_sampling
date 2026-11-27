"""Replica exchange: multi-replica trajectory reader and demultiplexing of per-replica data."""

import numbers
import warnings

import numpy as np
from MDAnalysis.analysis.base import AnalysisBase, Results
from MDAnalysis.auxiliary.XVG import XVGReader
from MDAnalysis.coordinates import base, core
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.lib import util
from MDAnalysis.lib.distances import calc_bonds
from MDAnalysis.lib.util import asiterable, store_init_arguments
from tqdm import tqdm

TIME_TOL = 1e-3  # ps, tolerance on the times of the replica index file


def regular_times(times, name="times"):
    """`times` (ps) of a regularly spaced series, without the drift of times computed as
    step × float32(dt) (relative error up to 6e-8, i.e. 0.5 ps at 10 μs): rescaled so
    that the spacing is a multiple of 1e-6 ps. The drift is measured over the whole
    series, so it must have no gaps (`ValueError` otherwise)."""
    times = np.asarray(times, dtype=np.float64)
    n = len(times)
    if n < 2:
        return times
    dt = (times[-1] - times[0]) / (n - 1)
    if np.abs(times - np.linspace(times[0], times[-1], n)).max() > dt / 4:
        raise ValueError(
            f"{name} could not be interpreted as a linearly spaced timeseries "
            f"({n} frames, {times[0]} ps to {times[-1]} ps)"
        )
    # The drift is proportional to the time since step 0: rescale the times themselves
    return times * (round(dt, 6) / dt)


def _replica_index_table(data, name):
    # Times (without drift) and permutations of a replica index file (time, then one
    # column per replica). After a restart, rows from the checkpoint time on are written
    # again: keep the rows earlier than all the following ones.
    times = data[:, 0]
    later = np.minimum.accumulate(times[::-1])[::-1]
    keep = np.append(times[:-1] < later[1:] - TIME_TOL, True)
    if not np.all(keep):
        warnings.warn(f"Found repeated times in {name}. Rows superseded by a restart are dropped.")
    return regular_times(times[keep], name), data[keep, 1:].astype(int)


def _permutation_rows(index_times, times):
    # Frames are written before the exchange of their step: at time t the permutation
    # is the last one strictly before t (the first one at t = 0)
    idx = np.searchsorted(index_times, np.asarray(times) - TIME_TOL, side="left") - 1
    return np.maximum(idx, 0)


class ParallelAnalysis:
    """Performs an analysis on each replica of a ParallelReader trajectory.
    Individual results are accessible with:
        parallel_analysis.results[i]
    The base Analysis class attributes and methods can be accessed:
        parallel_analysis.compute_something()
    returns
        [analysis.compute_something() for analysis in parallel_analysis._analyses]
    """

    def get_supported_backends():
        return ("serial",)

    def __init__(self, analysis_cls, *args, **kwargs):
        if not issubclass(analysis_cls, AnalysisBase):
            raise TypeError("analysis_cls must be a AnalysisBase object.")

        # Let the Analysis class determine the trajectory object
        _analysis = analysis_cls(*args, **kwargs)
        trajectory = _analysis._trajectory
        if not isinstance(trajectory, ParallelReader):
            raise TypeError("ParallelReader trajectory not found.")

        self._analyses = []
        for repid in range(trajectory.n_reps):
            trajectory.set_replica(repid)  # In case __init__ depends on replica
            self._analyses.append(analysis_cls(*args, **kwargs))

        self.analysis_cls = analysis_cls
        self.args = args
        self.kwargs = kwargs
        self._verbose = kwargs.get("verbose", False)

    def _prepare(self):
        self.results = Results()

    def run(self, **kwargs):
        self._prepare()
        for repid, analysis in enumerate(self._analyses):
            self._trajectory.set_replica(repid)
            if self._verbose:
                print(f"Running analysis for replica {repid} ({repid + 1}/{self.n_reps})")
            analysis.run(**kwargs)
            self.results[repid] = analysis.results
        return self

    @property
    def n_reps(self):
        return self._trajectory.n_reps

    def __getattr__(self, attr):
        if hasattr(self._analyses[0], attr):
            if callable(getattr(self._analyses[0], attr)):

                def _func(*args, **kwargs):
                    return [getattr(analysis, attr)(*args, **kwargs) for analysis in self._analyses]

                return _func
            else:
                attrs = [getattr(analysis, attr) for analysis in self._analyses]
                if all(at is attrs[0] for at in attrs):
                    return attrs[0]
                else:
                    return attrs

    def __repr__(self):
        return f"{type(self).__name__}({self.analysis_cls.__name__}) for {self._trajectory}"


class ParallelFrameIterator(base.FrameIteratorBase):
    """Iterable over the frames and replicas of a `ParallelReader` trajectory.
    # Iterate frames from one replica
    u.trajectory[:, 2]

    # Iterate replicas from one frame
    u.trajectory[100, :]

    # Slice and array-like
    u.trajectory[:, [3, 2, 1]]
    """

    def __init__(self, trajectory, frames, repids):
        super().__init__(trajectory)
        self._frame_indices = self._to_indices(frames, len(trajectory))
        self._rep_indices = self._to_indices(repids, trajectory.n_reps)

    def __len__(self):
        if isinstance(self._frame_indices, numbers.Integral):
            return len(self._rep_indices)
        else:
            return len(self._frame_indices)

    def _to_indices(self, inds, length):
        if isinstance(inds, numbers.Integral):
            inds = self._enforce_bounds(inds, length)
            return inds
        elif isinstance(inds, slice):
            return tuple(range(*inds.indices(length)))
        else:
            _inds = []
            for idx in inds:
                if not isinstance(idx, numbers.Integral):
                    raise TypeError("Frame and replica indices must be integers or slices.")
                idx = self._enforce_bounds(idx, length)
                _inds.append(idx)
            return tuple(_inds)

    def _enforce_bounds(self, idx, length):
        if not (-length <= idx < length):
            raise IndexError(f"Index {idx} out of bounds for indexable of length {length}")
        return idx % length

    def __iter__(self):
        # If single frame iterate on replicas and yield frames
        if isinstance(self._frame_indices, numbers.Integral):
            for repid in asiterable(self._rep_indices):
                yield self.trajectory[self._frame_indices, repid]
        else:  # iterate on frames and yield single-frame ParallelFrameIterator (or frame)
            for frame in self._frame_indices:
                yield self.trajectory[frame, self._rep_indices]

    def __getitem__(self, inds):
        if isinstance(inds, tuple) and len(inds) == 2:
            return self.trajectory[self._frame_indices[inds[0]], self._rep_indices[inds[1]]]
        elif isinstance(self._frame_indices, numbers.Integral):
            return self.trajectory[self._frame_indices, self._rep_indices[inds]]
        else:
            return self.trajectory[self._frame_indices[inds], self._rep_indices]


class ParallelReader(base.ReaderBase):
    format = "PARALLEL"

    @store_init_arguments
    def __init__(
        self,
        filenames,
        skip=1,
        convert_units=True,
        fileformat=None,
        replica_index=None,
        **kwargs,
    ):
        super().__init__(filename="PARALLEL", skip=skip, convert_units=convert_units, **kwargs)
        filenames = asiterable(filenames)

        kwargs["format"] = fileformat
        self.readers = [
            core.reader(filename, convert_units=convert_units, **kwargs) for filename in filenames
        ]
        self.filenames = np.array([fn[0] if isinstance(fn, tuple) else fn for fn in filenames])
        self.skip = skip
        self.n_atoms = self._get_same("n_atoms")
        self.n_frames = self._get_same("n_frames")

        self.ts = None
        self._active_replica = 0

        self._is_demuxed = False
        self._replica_index = None
        self._replica_temp = None

        if replica_index is not None:
            self.demux(replica_index)

    @staticmethod
    def _format_hint(thing):
        return (
            not isinstance(thing, np.ndarray) and util.iterable(thing) and not util.isstream(thing)
        )

    def __getstate__(self):
        state = self.__dict__.copy()
        state["ts"] = self.ts.__deepcopy__()
        for reader in state["readers"]:
            reader.rewind()
        self.ts = state["ts"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def _apply(self, method, **kwargs):
        """Execute `method` with `kwargs` for all readers."""
        return [reader.__getattribute__(method)(**kwargs) for reader in self.readers]

    def _get(self, attr):
        return [reader.__getattribute__(attr) for reader in self.readers]

    def _get_same(self, attr):
        """Verify that `attr` has the same value for all readers and return value."""
        values = np.array(self._get(attr))
        value = values[0]
        if not np.allclose(values, value):
            bad_traj = np.array(self.filenames)[values != value]
            raise ValueError(
                f"The following trajectories do not have the correct {attr}  ({value}):\n{bad_traj}"
            )
        return value

    def _read_frame(self, frame, i):
        """Get frame from replica `i`.
        If the trajectory is not demuxed, `i` corresponds to the "temperature".
        If the trajectory is demuxed, `i` corresponds to the "replica" (continuous trajectory).
        """
        if not -self.n_reps <= i < self.n_reps:
            raise IndexError(f"No replica for index {i}")
        i = i % self.n_reps
        reader_idx = self.replica_temp[frame, i] if self._is_demuxed else i
        self.ts = self.readers[reader_idx][frame]
        self.ts = self._apply_transformations(self.ts)
        self.ts.data["replica_index"] = i if self._is_demuxed else None
        self.ts.data["replica_temp"] = reader_idx
        self.ts.data["demuxed"] = self._is_demuxed
        return self.ts

    def __getitem__(self, inds):
        if isinstance(inds, tuple) and len(inds) == 2:
            frames, repids = inds
            if isinstance(frames, numbers.Integral) and isinstance(repids, numbers.Integral):
                return self._read_frame(frames, repids)
            else:
                return ParallelFrameIterator(self.trajectory, frames, repids)
        else:
            return self[inds, self._active_replica]

    def _read_next_timestep(self, ts=None):
        """This implementation is consistent with indexing at fixed replica index
            u.trajectory[:, i]
        If trajectory is demuxed, keeps constant replica identity, otherwise keeps constant
        "temperature"."""

        if ts is None:
            ts = self.ts
        if ts is None:
            # By default `for ts in u.trajectory` will iterate on replica 0
            # similarly to `for ts in u.trajectory[:, 0]`,
            # unless _active_replica is set.
            frame = 0
            i = self._active_replica if self._active_replica is not None else 0
        else:
            frame = ts.frame + 1
            i = ts.data["replica_index"] if self._is_demuxed else ts.data["replica_temp"]
        if frame >= self.n_frames:
            raise StopIteration
        self.ts = self._read_frame(frame, i)
        return self.ts

    def set_replica(self, i):
        if not -self.n_reps <= i < self.n_reps:
            raise IndexError(f"No replica for index {i}")
        i = i % self.n_reps
        self._active_replica = i

    def add_auxiliary(self, *args, **kwargs):
        raise NotImplementedError(f"Auxiliary API not supported yet for {self.__class__.__name__}")

    def demux(self, xvgfile=None):
        """Demultiplexing of the replicas using the replica indices as provided by Gromac's
        `demux.pl` script, i.e. a XVG file with time (in ps) as first column and replica indices
        as the remaining columns."""

        aux = XVGReader(xvgfile)

        if len(aux[0].data) != self.n_reps + 1:
            raise ValueError(
                f"The provided replica index file implies {len(aux[0].data) - 1}"
                f"replicas but the current trajectory has {self.n_reps}."
            )

        aux_data = np.array([aux_ts.data for aux_ts in aux])
        aux_times, indices = _replica_index_table(aux_data, xvgfile)
        traj_times = np.array([ts.time for ts in self.trajectory])
        if traj_times.min() < aux_times[0] - TIME_TOL:
            raise IndexError(f"No replica indices at time {traj_times.min()}")
        self._times = traj_times
        self._replica_index = indices[_permutation_rows(aux_times, traj_times)]
        self._replica_temp = np.argsort(self._replica_index, axis=1)
        self._is_demuxed = True

    @property
    def n_reps(self):
        return len(self.readers)

    @property
    def replica_index(self):
        return self._replica_index

    @property
    def replica_temp(self):
        return self._replica_temp

    @property
    def is_demuxed(self):
        return self._is_demuxed

    def check_demuxing(self, indices=None, skip=0, verbose=False):
        """Checks continuity using a heuristic based on atom displacement.
        If the trajectory is continuous the instantaneous MSD will never be
        much higher than MSD between any two replicas at t and t + 1."""

        if indices is None:
            indices = np.arange(self.n_atoms)[:10]

        old_pos = [ts.positions[indices] for ts in self.trajectory[skip, :]]
        old_dims = [ts.dimensions for ts in self.trajectory[skip, :]]
        for timesteps in tqdm(self.trajectory[skip:, :]):
            new_pos = [ts.positions[indices] for ts in timesteps]
            new_dims = [ts.dimensions for ts in timesteps]

            dist_matrix = np.array(
                [
                    [np.sum(calc_bonds(new, old, box=dims) ** 2) for new in new_pos]
                    for old, dims in zip(old_pos, old_dims)
                ]
            )
            if verbose:
                print(f"At time {timesteps[0].time}")
                print(dist_matrix)

            # Assumption: total displacement between t and t+1
            # is always the smallest for diagonal elements if the trajectory is
            # correctly demuxed. This might not be true in the beginning if all replicas
            # start from the same configuration, in which case we should skip equilibration
            # with `skip` argument.
            if not np.all(np.argmin(dist_matrix, axis=0) == np.arange(self.n_reps)):
                print(f"Trajectories seem discontinuous at step {timesteps[0].frame}")
                print("Cumulated distance matrix:")
                print(dist_matrix)
                return False

            old_pos = new_pos
            old_dims = new_dims
        return True

    def round_trip_time(self):
        """Estimates the round trip time, that is the time needed for a replica to go from the
        lowest thermodynamic state to the highest, and back."""

        if self.replica_index is None:
            raise NoDataError("Cannot compute round trip time without replica information.")

        round_trips = 0
        times = []
        i0 = self.replica_index[0, 0]
        t0 = 0
        up = True
        for frame, repids in enumerate(self.replica_index):
            if up and repids[-1] == i0:
                up = False
            elif not up and repids[0] == i0:
                round_trips += 1
                times.append(frame * self.dt - t0)
                t0 = frame * self.dt
                up = True

        if round_trips == 0:
            warnings.warn("Cannot estimate round trip time because no round trip was completed.")
            return None, None
        else:
            if round_trips < 5:
                warnings.warn(
                    "Unreliable estimate of the round trip mean time and error (n_trips < 5)"
                )
            avg = np.mean(times)
            err = np.std(times) / np.sqrt(len(times))
            print(
                f"Round-trip time estimate over {round_trips} round trips: {avg:.2f} ± {err:.2f} ps."
            )
            return avg, err

    def __iter__(self):
        self._reopen()
        return self

    def _reopen(self):
        """Internal method: Rewind trajectories themselves and trj pointer."""
        self._apply("rewind")
        self.ts = None

    def close(self):
        self._apply("close")

    def __repr__(self):
        return f"<{self.__class__.__name__} with {self.n_reps} trajectories of {self.n_frames} frames of {self.n_atoms} atoms>"


def demux(data_list, data_times=None, repindex_file="replica_index.xvg"):
    """
    From an array-like (n_rep, n_frames, n_features)
    outputs demuxed trajectories (n_rep, n_frames_processed, n_features)
    using the replica index file `repindex_file`.

    If `data_times` is `None` the first feature of the first replica is used. Both series
    of times must be regularly spaced; they are freed from float32 drift (`regular_times`).
    """
    repindex_times, repindex = _replica_index_table(np.loadtxt(repindex_file), repindex_file)
    n_reps = repindex.shape[1]

    data_arr = np.array(data_list)

    if data_times is None:
        data_times = data_arr[0, :, 0]
    data_times = regular_times(data_times, "data times")

    if data_arr.shape[0] != n_reps:
        raise ValueError(f"Number of replicas in {repindex_file} doesn't match `data_list`")

    n_frames = np.searchsorted(data_times, repindex_times[-1] + TIME_TOL, side="right")
    if n_frames < len(data_times):
        warnings.warn(
            f"Data times go above replica_index times. Data will be truncated after {repindex_times[-1]}."
        )
    reptemp = np.argsort(repindex[_permutation_rows(repindex_times, data_times[:n_frames])], axis=1)
    return data_arr[reptemp, np.arange(n_frames)[:, None]].swapaxes(0, 1)
