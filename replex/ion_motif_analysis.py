from mg_sampling import PROJECT_DIR
from mg_sampling.analysis import IonNucleicAnalysis
import MDAnalysis as mda

if __name__ == '__main__':

    #repids = [0, 11]
    #systems = ['7ez0', '7ez0_noion', '7ez0_2', '7ez0_noion_2']
    #forcefields = ['mMg_tip3p', 'mMg_opc', 'nMg_tip3p', 'ecc_spce']
    repids = [11]
    systems = ['7ez0_2', '7ez0_noion_2']
    forcefields = ['ecc_spce']

    n_workers = 32

    for repid in repids:
        for forcefield in forcefields:
            for system in systems:
                path = PROJECT_DIR / 'replex' / system / forcefield / f'rep{repid:02d}'
                u = mda.Universe(path / 'input.gro', path / 'run.xtc')
                anal = IonNucleicAnalysis(u, 'name MG NA CL', 2.5, ligsel='(nucleic and not type H) or name MG NA CL', outer_cutoff=4.8)
                anal.run(backend='dask', n_workers=n_workers)
                df = anal.to_pandas()
                df.to_csv(path / 'ion_analysis.csv')
