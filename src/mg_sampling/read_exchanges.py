#!/usr/bin/env python3

import argparse
import numpy as np

class GMXLogError(Exception):
    pass

class ReplexLog:
    def mean_acceptance(self, start=0, stop=-1):
        return np.sum(self.exchanges[start:stop+1], axis=0) / (~np.isnan(self.probas[start:stop+1])).sum(axis=0)

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
        for line in f:
            if line.startswith("Replica exchange at step"):
                split = line.split()
                step, time = int(split[4]), float(split[6])
                parity = (step // replex) % 2
                steps.append(step)
                times.append(time)
                parities.append(parity)
                for line in f:
                    if line.startswith("Repl ex"):
                        split = line.removeprefix("Repl ex").split("x")
                        exchange = (nrep - 1) * [False]
                        i = -1
                        for s in split[:-1]:
                            i += len(s.split())
                            exchange[i] = True
                        exchanges.append(exchange)
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

    args = parser.parse_args()

    rexlog = parse_logfile(args.filepath)
    print_summary(rexlog, start=args.start, stop=args.stop)

if __name__ == "__main__":
    main()

