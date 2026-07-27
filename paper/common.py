import itertools
from ast import literal_eval
import os
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
from scipy.stats import linregress, binned_statistic
from scipy.constants import Avogadro
from scipy.special import logsumexp
from scipy.optimize import minimize_scalar

from tqdm import tqdm

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import matplotlib.patheffects as PathEffects

import seaborn as sns

from bussilab import wham
import MDAnalysis as mda
from MDAnalysis.analysis.align import alignto

from mg_sampling import PROJECT_DIR
from mg_sampling.constants import RT
from mg_sampling.tools import read_plumed, read_plumed_fes, range_to_edges, edges_to_range, finitemin, finitemax, read_wsv, demux, block_bootstrap, bayesian_bootstrap, bootstrap
from mg_sampling.bias import bias_fn
from mg_sampling.density.align import reorder_residue_atoms, unwrap_sequential, wrap_around, wrap_around_atomgroup, TransformationFromFunction, standardize_names, check_alignable
from mg_sampling.density.grid import Grid, SparseGrid
from mg_sampling.density.forward_density import atomgroup_to_map
from mg_sampling.density.utils import indices_to_mask, voxels_around
from mg_sampling.density.plot import plot_map_slice

os.makedirs('figures', exist_ok=True)

sns.set_style("whitegrid")


plt.rcParams.update({
    'font.size': 8.0, 
    'figure.constrained_layout.use': True,
})

cm = 1 / 2.54  # inches

simlabels = {
    'plain_md_7ez0_mMg_tip3p': "μMg-MD",
    'plain_md_7ez0_mMg_tip3p_cryo': "μMg-MD + Cryo",
    'plain_md_7ez0_nMg_tip3p': "nMg-MD",
    'plain_md_7ez0_nMg_tip3p_cryo': "nMg-MD + Cryo",
    'plain_md_7ez0_nMg_tip3p_soft': "nMg-MD-soft",
    'plain_md_7ez0_noion_nMg_tip3p': "nMg-novo-MD",
    'plain_md_7ez0_noion_mMg_tip3p': "μMg-novo-MD",
    'replex_7ez0_mMg_tip3p': "μMg-REx",
    'replex_7ez0_mMg_tip3p_cryo': "μMg-REx + Cryo",
    'replex_7ez0_noion_mMg_tip3p': "μMg-novo-REx",
    'plain_md_7ez0_nomg_mMg_tip3p': "no-Mg-MD",
    'equil_7ez0_mMg_tip3p': "μMg-Eq",
    'equil_7ez0_nMg_tip3p': "nMg-Eq",
    'equil_7ez0_noion_mMg_tip3p': "μMg-novo-Eq",
    'equil_7ez0_noion_nMg_tip3p': "nMg-novo-Eq",
    'replex_7ez0_mMg_tip3p_0': "μMg-REx [1]",
    'replex_7ez0_noion_mMg_tip3p_0': "μMg-novo-REx [1]",
    'replex_7ez0_mMg_tip3p_1': "μMg-REx [2]",
    'replex_7ez0_noion_mMg_tip3p_1': "μMg-novo-REx [2]",
}

simcolors = {
    'plain_md_7ez0_mMg_tip3p': "dimgray",
    'plain_md_7ez0_mMg_tip3p_cryo': "dimgray",
    'plain_md_7ez0_noion_mMg_tip3p': "lightgray",
    'plain_md_7ez0_nMg_tip3p': "mediumpurple",
    'plain_md_7ez0_nMg_tip3p_cryo': "mediumpurple",
    'plain_md_7ez0_nMg_tip3p_soft': "mediumseagreen",
    'plain_md_7ez0_noion_nMg_tip3p': "royalblue",
    'plain_md_7ez0_nomg_mMg_tip3p': "skyblue",
    'replex_7ez0_mMg_tip3p': "firebrick",
    'replex_7ez0_mMg_tip3p_cryo': "firebrick",
    'replex_7ez0_noion_mMg_tip3p': "tomato",
}

labelcolors = {
    v: simcolors[k] for k, v in simlabels.items() if k in simcolors
}

ioncolors = {
    'NA': 'darkorchid',
    'MG': 'darkgreen',
    'CL': 'greenyellow',
}

ionlabels = {
    'NA': r'$\mathrm{Na^{+}}$',
    'MG': r'$\mathrm{Mg^{2\!+}}$',
    'CL': r'$\mathrm{Cl^{-}}$',
    'ztot': r'$Z_\mathrm{tot}$',
}

# Discard equilibration [ns]
simstart = {
    'plain_md_7ez0_mMg_tip3p': 150,
    'plain_md_7ez0_nMg_tip3p': 150,
    'plain_md_7ez0_nMg_tip3p_soft': 150,
    'plain_md_7ez0_noion_mMg_tip3p': 150,
    'plain_md_7ez0_noion_nMg_tip3p': 150,
    'plain_md_7ez0_nomg_mMg_tip3p': 100,
    'replex_7ez0_mMg_tip3p': 50,
    'replex_7ez0_noion_mMg_tip3p': 50,
}

# Set this to True to (re)compute expensive stuff
compute = False

# simselection = [
#     'plain_md_7ez0_mMg_tip3p',
#     'plain_md_7ez0_nMg_tip3p',
#     'plain_md_7ez0_nMg_tip3p_soft',
#     'plain_md_7ez0_nomg_mMg_tip3p',
#     'replex_7ez0_mMg_tip3p',  
# ]
# suffix = '_talk'
simselection = [
    'plain_md_7ez0_mMg_tip3p',
    'plain_md_7ez0_noion_mMg_tip3p',
    'plain_md_7ez0_nMg_tip3p',
    'plain_md_7ez0_noion_nMg_tip3p',
    'replex_7ez0_mMg_tip3p',
    'replex_7ez0_noion_mMg_tip3p',
    'plain_md_7ez0_nMg_tip3p_soft',
    'plain_md_7ez0_nomg_mMg_tip3p',
]
labelselection = [simlabels[sim] for sim in simselection]
suffix = ''


roman = "i ii iii iv v vi vii viii ix x xi xii xiii xiv xv xvi xvii xviii xix xx".split()