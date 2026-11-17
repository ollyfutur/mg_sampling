"""Simulated systems: production window and replica folders."""

from dataclasses import dataclass
from pathlib import Path

from pipeline import PROJECT_DIR

N_REPLICAS = 8  # replicas per replex replicate


@dataclass(frozen=True)
class System:
    start: float  # ns, end of equilibration
    end: float  # ns, end of the run
    n_replicates: (
        int  # replex: replicates of N_REPLICAS replicas; plain MD: independent runs
    )


SYSTEMS = {
    "replex/7ez0/mMg_tip3p": System(75, 300, 2),
    "replex/7ez0_noion/mMg_tip3p": System(75, 300, 2),
    "plain_md/7ez0/mMg_tip3p": System(150, 600, 8),
    "plain_md/7ez0/nMg_tip3p": System(150, 600, 8),
    "plain_md/7ez0/nMg_tip3p_soft": System(150, 600, 8),
    "plain_md/7ez0_noion/mMg_tip3p": System(150, 600, 8),
    "plain_md/7ez0_noion/nMg_tip3p": System(150, 600, 8),
    "plain_md/7ez0_nomg/mMg_tip3p": System(100, 400, 8),
}


def is_replex(system: str) -> bool:
    return system.startswith("replex")


def replica_dirs(system: str, k: int) -> list[Path]:
    """Folders of the continuous trajectories of replicate k (replex: N_REPLICAS replicas; plain MD: one)."""
    path = PROJECT_DIR / system
    if is_replex(system):
        return [
            path / f"replicate{k}" / f"rep{i:02d}"
            for i in range(N_REPLICAS * k, N_REPLICAS * (k + 1))
        ]
    return [path / "replicate0" / f"rep{k:02d}"]
