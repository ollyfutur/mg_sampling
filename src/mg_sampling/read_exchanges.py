#!/usr/bin/env python3

import argparse
import numpy as np

class GMXLogError(Exception):
    pass

class ReplexLog:
    def mean_acceptance(self, start=0, stop=-1):
        return np.sum(self.exchanges[start:stop+1], axis=0) / (~np.isnan(self.probas[start:stop+1])).sum(axis=0)

    def count_half_trips(self, start=0, stop=-1):
        is_low = self.repindices[start:stop+1, 0, None] == np.arange(self.nrep)
        is_high = self.repindices[start:stop+1, -1, None] == np.arange(self.nrep)
        states = (-1 * is_low + 1 * is_high)
        states = [s[s != 0] for s in states.T]
        n_half_trips = [(s[1:] != s[:-1]).sum() for s in states]
        return np.array(n_half_trips)

    def count_round_trips(self, start=0, stop=-1):
        return self.count_half_trips(start, stop) // 2

    def occupancy_matrix(self, start=0, stop=-1):
        return np.array([
            (self.repindices[start:stop+1, i, None] == np.arange(self.nrep)).mean(axis=0)
            for i in range(self.nrep)
        ])

    def __len__(self):
        return len(self.steps)

def parse_logfile(path: str) -> ReplexLog:
    with open(path, 'r') as f:

        for line in f:
            if line.startswith("Repl  There are"):
                nrep = int(line.split()[3])
                break
        else:
            raise GMXLogError("Number of replicas not found")

        for line in f:
            if line.startswith("Replica exchange in temperature"):
                line = next(f)
                temps = list(map(float, line.split()))
                break
        else:
            raise GMXLogError("Temperatures not found")

        for line in f:
            if line.startswith("Replica exchange interval"):
                replex = int(line.split()[3])
                break
        else:
            raise GMXLogError("Temperatures not found")

        steps = []
        times = []
        parities = []
        exchanges = []
        probas = []
        repindices = [np.arange(nrep)]
        for line in f:
            if line.startswith("Replica exchange at step"):
                split = line.split()
                step, time = int(split[4]), float(split[6])
                parity = (step // replex) % 2  # 1 odd, 0 even
                steps.append(step)
                times.append(time)
                parities.append(parity)
                for line in f:
                    if line.startswith("Repl ex"):
                        repindex = repindices[-1].copy()
                        split = line.removeprefix("Repl ex").split("x")
                        exchange = (nrep - 1) * [False]
                        i = -1
                        for s in split[:-1]:
                            i += len(s.split())
                            exchange[i] = True
                            repindex[i], repindex[i+1] = repindex[i+1], repindex[i]
                        exchanges.append(exchange)
                        repindices.append(repindex)
                        break
                for line in f:
                    if line.startswith("Repl pr"):
                        split = list(map(float, line.removeprefix("Repl pr").split()))
                        proba = (nrep -1) * [None]
                        proba[(1-parity)::2] = split
                        probas.append(proba)
                        break

    rexlog = ReplexLog()
    rexlog.nrep = nrep
    rexlog.replex = replex
    rexlog.steps = np.array(steps, dtype=int)
    rexlog.times = np.array(times, dtype=float)
    rexlog.parities = np.array(parities, dtype=bool)
    rexlog.exchanges = np.array(exchanges, dtype=bool)
    rexlog.probas = np.array(probas, dtype=float)
    rexlog.repindices = np.array(repindices, dtype=int)

    return rexlog

def parse_logfile_biasvalue(path: str) -> ReplexLog:
    """
    Parse a GROMACS .log file where PLUMED's REPLEX_BIASVALUE was used.
    """
    with open(path, 'r') as f:

        for line in f:
            if line.startswith("PLUMED: Action REPLEX_BIASVALUE"):
                break
        else:
            raise GMXLogError("PLUMED Action REPLEX_BIASVALUE not found.")

        for line in f:
            if line.startswith("PLUMED:   Number of replicas"):
                nrep = int(line.split()[4])
                break
        else:
            raise GMXLogError("Number of replicas not found")

        for line in f:
            if line.startswith("PLUMED:   Exchange frequency"):
                replex = int(line.split()[3])
                break
        else:
            raise GMXLogError("Exchange frequency not found")

        for line in f:
            if line.startswith("PLUMED:   Number of exchange"):
                nex = int(line.split()[-1])
                break
        else:
            raise GMXLogError("Number of attempts per exchange step not found")

        for line in f:
            if line.startswith("PLUMED: Timestep"):
                dt = float(line.split()[2])
                break
        else:
            raise GMXLogError("Timestep not found")


        steps = []
        times = []
        parities = []
        exchanges = []
        probas = []
        repindices = [np.arange(nrep)]
        for line in f:
            if line.startswith("PLUMED: Exchange step"):
                split = line.split()
                step = int(split[3].removesuffix(','))
                time = float(split[5])
                parity = 1 - ((step // replex) - 1) * nex % 2  # 1 odd, 0 even
                for attempt in range(nex):
                    try:
                        for line in f:
                            if line.startswith("PLUMED: pacc:"):
                                split = list(map(float, line.removeprefix("PLUMED: pacc:").split()))
                                proba = (nrep-1) * [None]
                                proba[(1-parity)::2] = split
                                break
                        else:
                            raise GMXLogError("Expected line with pacc, not found")

                        line = next(f)
                        if line.startswith("PLUMED: acc:"):
                            exchange = (nrep - 1) * [False]
                            split = list(map(lambda x: True if x == '+' else False, line.removeprefix("PLUMED: acc:").split()))
                            exchange[(1-parity)::2] = split
                            repindex = repindices[-1].copy()
                            for i in range(nrep-1):
                                if exchange[i]:
                                    repindex[i], repindex[i+1] = repindex[i+1], repindex[i]
                        else:
                                raise GMXLogError(f"Expected acc after pacc, instead found: {line}")

                        steps.append(step)
                        times.append(time)
                        parities.append(parity)
                        probas.append(proba)
                        exchanges.append(exchange)
                        repindices.append(repindex)

                        parity = not parity
                    except (GMXLogError, ValueError) as e:
                        print(f"Stopped at step {step}: {e}")
                        break
                else:
                    for line in f:
                        if line.startswith("PLUMED: New replica indices:"):
                            _repindex = np.array(list(map(int, line.removeprefix("PLUMED: New replica indices:").split())))
                            if np.any(_repindex != repindex):
                                raise GMXLogError(f"Discrepancy between inferred repindex ({repindex}) and logged repindex ({_repindex})")
                            break
                    else:
                        raise GMXLogError("Did not find new replica indices")


    rexlog = ReplexLog()
    rexlog.nrep = nrep
    rexlog.replex = replex
    rexlog.nex = nex
    rexlog.steps = np.array(steps, dtype=int)
    rexlog.times = np.array(times, dtype=float)
    rexlog.parities = np.array(parities, dtype=bool)
    rexlog.exchanges = np.array(exchanges, dtype=bool)
    rexlog.probas = np.array(probas, dtype=float)
    rexlog.repindices = np.array(repindices, dtype=int)

    return rexlog


def print_summary(rexlog: ReplexLog, start:float=0, stop:float=np.inf):
    i = np.searchsorted(rexlog.times, start)
    j = min(np.searchsorted(rexlog.times, stop), len(rexlog.times)-1)

    print("Number of replicas:", rexlog.nrep)
    print("Exchange interval:", rexlog.replex)
    print("Between step", rexlog.steps[i],
          f"({rexlog.times[i]} ps) and step", rexlog.steps[j],
          f"({rexlog.times[j]} ps)")
    print("Total attempts:", len(rexlog.steps[i:j+1]))
    print("  odd:", rexlog.parities[i:j+1].sum())
    print("  even:", (~rexlog.parities[i:j+1]).sum())
    print("Half trips (number of times replica goes from cold to hot or vice-versa):")
    print("  ".join(f"{i:>2d}" for i in range(rexlog.nrep)))
    print("  ".join(f"{k:>2d}" for k in rexlog.count_half_trips(i, j)))
    print("Occupancy matrix:")
    print("  Replica:", "  ".join(f"{i:>4d}" for i in range(rexlog.nrep)))
    for temp, occupancy in enumerate(rexlog.occupancy_matrix(i, j)):
        print(f"  T{temp:<6d}:", "  ".join(f"{occ:.2f}" for occ in occupancy))
    print("Final replica indices:")
    print("  ".join(f"{i:>2d}" for i in rexlog.repindices[-1]))
    print("Average probabilities:")
    print("  ".join(f"{a:>2d} ↔ {a+1:<2d}" for a in range(rexlog.nrep - 1)))
    print("  ".join(f"{p:^7.2f}" for p in np.nanmean(rexlog.probas[i:j+1], axis=0)))
    print("Total exchanges:")
    print("  ".join(f"{a:>2d} ↔ {a+1:<2d}" for a in range(rexlog.nrep - 1)))
    print("  ".join(f"{n:^7d}" for n in np.sum(rexlog.exchanges[i:j+1], axis=0)))
    print("Average acceptance:")
    print("  ".join(f"{a:>2d} ↔ {a+1:<2d}" for a in range(rexlog.nrep - 1)))
    avg_acc = np.sum(rexlog.exchanges[i:j+1], axis=0) / (~np.isnan(rexlog.probas[i:j+1])).sum(axis=0)
    print("  ".join(f"{p:^7.2f}" for p in avg_acc))


def main():
    parser = argparse.ArgumentParser(
        prog="read_exchanges",
        description="Extract replica exchange statistics from a Gromacs log file.",
    )
    parser.add_argument('filepath', help="Gromacs .log file of a replica exchange simulation")
    parser.add_argument('--start', type=float, default=0.0, help="(0) Starting time [ps]")
    parser.add_argument('--stop', type=float, default=float('inf'), help="(inf) Ending time [ps]")
    parser.add_argument('--biasvalue', action='store_true', help="Reads REPLEX_BIASVALUE output.")

    args = parser.parse_args()

    if args.biasvalue:
        rexlog = parse_logfile_biasvalue(args.filepath)
    else:
        rexlog = parse_logfile(args.filepath)
    print_summary(rexlog, start=args.start, stop=args.stop)

if __name__ == "__main__":
    main()

