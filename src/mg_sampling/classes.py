from pathlib import Path
import glob

import numpy as np
import pandas as pd

from mg_sampling.tools import read_plumed, read_plumed_fes
from mg_sampling.constants import RT

class Sim:
    """
    A base class wrapping a Gromacs MD simulation.
    """
    def __init__(self, path: Path | str, name: str = "Simulation", colvar_file: str | list[str] = 'COLVAR'):
        self.path = Path(path)
        self.name = name
        self.colvar_file = colvar_file

        if isinstance(colvar_file, list):
            df = read_plumed(self.path / colvar_file[0]).drop_duplicates(subset='time')
            for file in colvar_file[1:]:
                df = pd.merge(df, read_plumed(self.path / file).drop_duplicates(subset='time'))
        else:
            df = read_plumed(self.path / self.colvar_file).drop_duplicates(subset='time')
        self.colvar = df

class MetadSim(Sim):
    def __init__(self, path: Path | str, name: str = "Metadynamics", colvar_file: str = 'static.COLVAR', fes_file: str = 'fes.dat'):
        super().__init__(path, name, colvar_file)

        self.fes_file = fes_file
        (self.x_name, self.y_name), (self.x_range, self.y_range), self.fes = read_plumed_fes(self.path / fes_file)

        self.x = self.colvar[self.x_name].to_numpy()
        self.y = self.colvar[self.y_name].to_numpy()
        self.X = self.x - self.y + 6
        self.Y = self.x + self.y - 6
        self.bias = self.colvar['metad.bias'].to_numpy()
        self.weights = np.exp(self.bias / RT)
    
class KineticSim(Sim):
    def __init__(self, path: Path | str, name: str = "Kinetics", colvar_file: str = 'COLVAR'):
        super().__init__(path, name, colvar_file)
        self.scaling = self.path.name

class KineticSims:
    def __init__(self, path: Path | str, name: str = "Kinetics", colvar_file: str = 'COLVAR'):
        self.path = Path(path)
        self.scalings = glob.glob('[0-9].[0-9]', root_dir=self.path)
        self.scalings.sort(key=float)
        self.sims = [
            KineticSim(self.path / scaling, colvar_file=colvar_file)
            for scaling in self.scalings
        ]
    
    def __getitem__(self, i: int):
        return self.sims[i]

class ReplexSims:
    def __init__(self, path: Path | str, name: str = "Replex"):
        self.path = Path(path)
        self.repdirs = glob.glob('rep*', root_dir=self.path).sort()
        self.nrep = len(self.repdirs)
        self.sims = [
            Sim(
                self.path / repdir,
                f"{name} (Replica {i})",
                colvar_file=[f'bias.COLVAR.{i}', f'rmsd.COLVAR.{i}'],
            )
            for i, repdir in enumerate(self.repdirs)
        ]
