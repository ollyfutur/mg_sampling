import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import scienceplots
plt.style.use(['science', 'notebook'])
plt.rcParams.update({
    "text.usetex": False,
})

from mg_sampling.classes import Sim, MetadSim, KineticSims
from mg_sampling.tools import range_to_edges, edges_to_range, finitemin, finitemax, detect_round_trips, count_round_trips, kish_ratio
from mg_sampling.constants import RT, SWITCH_COEFF

from scipy.signal import medfilt

LABELS = {
    'nop': r'$n_\mathrm{O_P}$',
    'nob': r'$n_\mathrm{O_B}$',
    'now': r'$n_\mathrm{O_W}$',
    'op': r'$\mathrm{O_P}$',
    'O_Phosphate': r'$\mathrm{O_P}$',
    'ob': r'$\mathrm{O_B}$',
    'O_Base': r'$\mathrm{O_B}$',
    'amber_tip3p': "Amber/TIP3P",
    'mMg_tip3p': "microMg/TIP3P",
    'nMg_tip3p': "nanoMg/TIP3P", 
    'mMg_opc': "microMg/OPC",
    'ecc_spce': "AmberECC/SPCE",
    'dop1': r'$d(\mathrm{Mg},\mathrm{O_P})$ [nm]',
    'bias.bias': r'Bias [kJ/mol]',
    'X': r'$X$',
    'Y': r'$Y$',
    'MG': 'Mg²⁺',
    'NA': 'Na⁺',
    '7ez0': '7ez0 (1)',
    '7ez0_2': '7ez0 (2)',
    '7ez0_noion': '7ez0 (1) (w/o PDB ions)',
    '7ez0_noion_2': '7ez0 (2) (w/o PDB ions)',
}

def plot_metad_fes(sim: MetadSim, method='hills', ax=None, bins=10, **kwargs):

    if ax is None:
        _, ax = plt.subplots()

    if method == 'hills':
        fes = sim.fes - sim.fes.min()
        x_edges, y_edges = range_to_edges(sim.x_range), range_to_edges(sim.y_range)
    elif method == 'reweight':
        x, y = sim.colvar[sim.x_name], sim.colvar[sim.y_name]
        bias = sim.colvar['metad.bias']
        weights = np.exp(bias/RT)
        weights /= weights.mean()
        ρ_avg, x_edges, y_edges = np.histogram2d(x, y, bins=bins, density=True, weights=weights)
        fes = -RT * np.log(ρ_avg)
        fes -= finitemin(fes)
    
    extent = [*x_edges[[0, -1]], *y_edges[[0, -1]]]

    im = ax.imshow(fes.T, extent=extent, origin='lower', aspect=1, cmap='RdPu_r')
    levels = RT * np.arange(np.floor(finitemin(fes)/RT), np.ceil(finitemax(fes)/RT))
    cf = ax.contour(fes.T, colors='black', linewidths=[1, 0.5], extent=extent, levels=levels)
    cbar = plt.colorbar(im)
    cbar.add_lines(cf)
    cbar.ax.set_ylabel(r'FES [kJ/mol]')
    ax.set_title(sim.name)
    ax.set_xlabel(LABELS[sim.x_name])
    ax.set_ylabel(LABELS[sim.y_name])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(color='black', linestyle='solid', linewidth=0.5)
    ax.xaxis.set_zorder(10)
    ax.yaxis.set_zorder(10)

    return ax

def plot_metad_fes_marginals(sim: MetadSim, which='x', method='hills', bins=10, ax=None, **kwargs):
    
    if ax is None:
        _, ax = plt.subplots()

    if method == 'hills':
        fes = sim.fes - sim.fes.min()
        x_range, y_range = sim.x_range, sim.y_range
    elif method == 'reweight':
        x, y = sim.colvar[sim.x_name], sim.colvar[sim.y_name]
        bias = sim.colvar['metad.bias']
        weights = np.exp(bias/RT)
        weights /= weights.mean()
        ρ_avg, x_edges, y_edges = np.histogram2d(x, y, bins=bins, density=True, weights=weights)
        x_range, y_range = edges_to_range(x_edges), edges_to_range(y_edges)
        fes = -RT * np.log(ρ_avg)
        fes -= finitemin(fes)
    
    ρ = np.exp(-fes/RT)

    match which:
        case 'x':
            i = 1
            i_range = x_range
            i_name = sim.x_name
        case 'y':
            i = 0
            i_range = y_range
            i_name = sim.y_name

    ρ_i = ρ.sum(axis=i)
    fes_i = -RT*np.log(ρ_i)
    fes_i -= fes_i.min()

    ax.plot(i_range, fes_i, **kwargs)
    ax.set_ylabel(r'FES [kJ/mol]')
    ax.set_xlabel(LABELS[i_name])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    return ax

def plot_metad_timeseries(sim: MetadSim, axarr=None):
    x_name, y_name = sim.x_name, sim.y_name
    x, y = sim.colvar[x_name], sim.colvar[y_name]
    bias = sim.colvar['metad.bias']
    t_range = sim.colvar['time'] / 1000

    data = [
        {'y': x, 'label': LABELS[x_name]},
        {'y': y, 'label': LABELS[y_name]},
        {'y': bias, 'label': r'Bias [kJ/mol]'},
    ]

    if axarr is None:
        _, axarr = plt.subplots(len(data), figsize=(6, 10), sharex=True)

    for ax, d in zip(axarr, data):
        ax.plot(t_range, d['y'], lw=1, alpha=0.5, rasterized=True)
        ax.scatter(t_range, d['y'], s=0.1, alpha=0.5, c='k', rasterized=True)
        ax.set_ylabel(d['label'])
        ax.grid(axis='y')

    axarr[-1].set_xlabel(r'Time [ns]')
    return axarr


def plot_bias(sim: MetadSim, bias_fn, params, bins=10, ax=None):

    if ax is None:
        ax = plt.subplot()

    bias_weights = np.exp(-bias_fn(sim.x, sim.y, *params)/RT)
    ρxy, x_edges, y_edges = np.histogram2d(sim.x, sim.y, bins=bins, density=True, weights=sim.weights*bias_weights)
    x_range, y_range = edges_to_range(x_edges), edges_to_range(y_edges)
    extent = [*x_edges[[0, -1]], *y_edges[[0, -1]]]
    F = -RT*np.log(ρxy)
    xx, yy = np.meshgrid(x_range, y_range, indexing='ij')
    u_bias = bias_fn(xx, yy, *params)

    im = ax.imshow(F.T, extent=extent, origin='lower', aspect=1, cmap='RdPu_r')
    levels = RT * np.arange(np.floor(finitemin(F)/RT), np.ceil(finitemax(F)/RT))
    colors = len(levels) * ['black']
    cf = ax.contour(F.T, colors=colors, linewidths=[1, 0.5], extent=extent, levels=levels)
    bias_levels = RT * np.arange(np.floor(finitemin(u_bias)/RT), np.ceil(finitemax(u_bias)/RT))
    ax.contour(u_bias.T, origin='lower', extent=extent, colors='white', linewidths=1, levels=bias_levels)
    ax.contour(u_bias.T, origin='lower', extent=extent, colors='red', linewidths=1, levels=[-0.00000000000001])
    cbar = plt.colorbar(im)
    cbar.add_lines(cf)
    cbar.ax.set_ylabel(r'FES [kJ/mol]')
    ax.set_title("Reweighting")
    ax.set_xlabel(LABELS[sim.x_name])
    ax.set_ylabel(LABELS[sim.y_name])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(color='black', linestyle='solid', linewidth=0.5)
    ax.xaxis.set_zorder(10)
    ax.yaxis.set_zorder(10)

    return ax

def plot_bias_XY_marginals(sim: MetadSim, bias_fn, params, X_bins=10, Y_bins=10, mask=None, axes=None):
    if axes is None:
        fig, axes = plt.subplots(ncols=2, sharey=True)

    if mask is None:
        mask = (sim.X > 0) * (sim.X < 2) * (sim.Y > 0) * (sim.Y < 2)
    elif callable(mask):
        mask = mask(sim.x, sim.y)

    E, *params = params
    for scale in np.linspace(0, 1, 5):
        bias_weights = np.exp(-bias_fn(sim.x, sim.y, scale*E, *params)/RT) * mask

        ρX, x_edges = np.histogram(sim.X, bins=X_bins, density=True, weights=sim.weights*bias_weights)
        x_range = edges_to_range(x_edges)
        FX = -RT*np.log(ρX)
        FX -= FX.min()
        ρY, y_edges = np.histogram(sim.Y, bins=Y_bins, density=True, weights=sim.weights*bias_weights)
        y_range = edges_to_range(y_edges)
        FY = -RT*np.log(ρY)
        FY -= FY.min()

        ax = axes[0]
        ax.plot(x_range, FX, label=f"{scale}")
        ax.set_xlabel("X")
        ax.set_ylabel("FES [kJ/mol]")

        ax = axes[1]
        ax.plot(y_range, FY)
        ax.set_xlabel("Y")
        ax.set_ylabel("FES [kJ/mol]")
    
    return axes

def plot_kinetics_timeseries(sims: KineticSims):

    fig, axarr = plt.subplots(nrows=len(sims.scalings), ncols=6, figsize=(15, 10), sharex=True, sharey='col')

    for i, (scaling, _axarr) in enumerate(zip(sims.scalings, axarr)):

        dop1 = sims[i].colvar['dop1']
        nop = sims[i].colvar['cv.nb-0']
        now = sims[i].colvar['cv.nc-0']
        now[now == 0] = 6
        X = SWITCH_COEFF * nop - SWITCH_COEFF * now + 6
        Y = SWITCH_COEFF * nop + SWITCH_COEFF * now - 6
        bias = sims[i].colvar['bias.bias']
        t_range = sims[i].colvar['time'] / 1000

        ax = _axarr[0]
        ax.plot(t_range, dop1, lw=1, alpha=0.5, rasterized=True)
        ax.scatter(t_range, dop1, s=0.1, alpha=0.1, c='k', rasterized=True)
        ax.grid(axis='y')
        ax.set_ylim(0, 1.5)
        

        ax = _axarr[1]
        ax.plot(t_range, nop, lw=1, alpha=0.5, rasterized=True)
        ax.scatter(t_range, nop, s=0.1, alpha=0.1, c='k', rasterized=True)
        ax.grid(axis='y')

        ax = _axarr[2]
        ax.plot(t_range, now, lw=1, alpha=0.5, rasterized=True)
        ax.scatter(t_range, now, s=0.1, alpha=0.1, c='k', rasterized=True)
        ax.grid(axis='y')
            
        ax = _axarr[3]
        ax.plot(t_range, X, lw=1, alpha=0.5, rasterized=True)
        ax.scatter(t_range, X, s=0.1, alpha=0.1, c='k', rasterized=True)
        Xfilt = medfilt(X, 9)
        ax.plot(t_range, Xfilt, lw=1, alpha=1, color='red', zorder=10, rasterized=True)
        inds = detect_round_trips(Xfilt, 0, 2)
        ax.scatter(t_range[inds], Xfilt[inds], s=0.5, c='red', marker='x', rasterized=True)
        ax.grid(axis='y')        
        
        ax = _axarr[4]
        ax.plot(t_range, Y, lw=1, alpha=0.5, rasterized=True)
        ax.scatter(t_range, Y, s=0.1, alpha=0.1, c='k', rasterized=True)
        ax.grid(axis='y')

        ax = _axarr[5]
        ax.plot(t_range, bias, lw=1, alpha=0.5, color='darkred', rasterized=True)
        ax.scatter(t_range, bias, s=0.1, alpha=0.1, c='k', rasterized=True)
        ax.grid(axis='y')
    
        _axarr[0].set_ylabel(f'λ = {scaling}')

    axarr[0, 0].set_title(LABELS['dop1'])
    axarr[0, 1].set_title(LABELS['nop'])
    axarr[0, 2].set_title(LABELS['now'])
    axarr[0, 3].set_title(LABELS['X'])
    axarr[0, 4].set_title(LABELS['Y'])
    axarr[0, 5].set_title("Bias [kJ/mol]")

    return fig, axarr

def analyse_kinetics(sims: KineticSims):
    data = []
    for i, scaling in enumerate(sims.scalings):
        dop1 = sims[i].colvar['dop1']
        nop = sims[i].colvar['cv.nb-0']
        now = sims[i].colvar['cv.nc-0']
        now[now == 0] = 6
        X = SWITCH_COEFF * nop - SWITCH_COEFF * now + 6
        Y = SWITCH_COEFF * nop + SWITCH_COEFF * now - 6
        bias = sims[i].colvar['bias.bias']
        weights = np.exp(bias/RT)
        t_range = sims[i].colvar['time'] / 1000
        
        Xfilt = medfilt(X, 9)
        n_trips = count_round_trips(X, 0, 2)
        trip_rate = n_trips / (t_range.iloc[-1] - t_range.iloc[0])
        kr = kish_ratio(weights)
        data.append({
            'scaling': scaling,
            'round_trips': n_trips,
            'trip_rate': trip_rate,
            'kish_ratio': kr,
        })
    df = pd.DataFrame(data)

    return df

def plot_kinetics(sims: KineticSims, ax=None):
    df = analyse_kinetics(sims)

    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 4))

    ax2 = ax.twinx()

    ax.plot(df['scaling'], df['kish_ratio'], linestyle='dotted', color='darkred', marker='.')
    ax.set_xlabel("Bias scaling")
    ax.set_ylabel("Kish ratio")
    ax.set_ylim(0, 1)
    ax.tick_params(axis='y', colors='darkred')

    ax2.plot(df['scaling'], df['trip_rate'], color='black', marker='.')
    ax2.set_ylabel(r"Binding event rate [ns$^{-1}$]")
    ax2.set_ylim(0, None)

    return ax

def plot_site_heatmap(sim: Sim, site='O_Phosphate', ion='MG', ax=None):
    data = np.loadtxt(f'{sim.path}/dist_{site}_{ion}.xvg', comments=['@', '#'])
    t_range = data[:, 0] / 1000
    dist_op_mg = data[:, 1:]
    state = np.digitize(dist_op_mg, (0.25, 0.5)) # 0 = inner shell, 1 = outer shell, 2 = unbound
    nres = dist_op_mg.shape[1] // 2
    extent = [t_range[0], t_range[-1], 1, nres]

    if ax is None:
        _, ax = plt.subplots()

    ax.imshow(state.T, origin='lower', aspect='auto', cmap='afmhot', extent=extent, vmin=0, vmax=2, interpolation='nearest')
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel("Site")
    ax.set_title(f"{site}—{ion}")

    return ax

def plot_site_occupancy(sim: Sim, site='O_Phosphate', ion='MG', ax=None):
    data = np.loadtxt(f'{sim.path}/dist_{site}_{ion}.xvg', comments=['@', '#'])
    t_range = data[:, 0] / 1000
    dist_op_mg = data[:, 1:]
    state = np.digitize(dist_op_mg, (0.25, 0.5)) # 0 = inner shell, 1 = outer shell, 2 = unbound
    nres = dist_op_mg.shape[1] // 2
    resids = range(1, nres+1)

    if ax is None:
        _, ax = plt.subplots()
    
    ref = (state[0] == 0)
    avg = (state == 0).mean(axis=0)
    for i, (a, b) in zip(resids, zip(avg, ref)):
        ax.plot([i, i], [a, b], color='gray')
        scref = ax.scatter(i, b, marker='x', color='gray', label="PDB")
        scavg = ax.scatter(i, a, marker='.', c='k', label="Average", zorder=10)
        if (a > 0.5 or b > 0.5): ax.text(i, a + 0.05*(i%2), f'{i}', fontsize=3, fontweight='bold', rotation=90, color='darkred', zorder=20)
    ax.set_xlabel("Site")
    ax.set_ylabel(r"$p_\mathrm{inner}$")
    ax.set_title(f"{site}—{ion}")

    return ax


def plot_site_global_occupancy(sim: Sim, site='O_Phosphate', ion='MG', ax=None, **kwargs):
    data = np.loadtxt(f'{sim.path}/dist_{site}_{ion}.xvg', comments=['@', '#'])
    t_range = data[:, 0] / 1000
    dist_op_mg = data[:, 1:]
    state = np.digitize(dist_op_mg, (0.25, 0.5)) # 0 = inner shell, 1 = outer shell, 2 = unbound
    nres = dist_op_mg.shape[1] // 2

    if ax is None:
        _, ax = plt.subplots()
    
    ax.plot(t_range, (state == 0).mean(axis=1), **kwargs)
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel(r"$p_\mathrm{inner}$")
    ax.set_title(f"{site}—{ion}")
    ax.set_ylim(0, 0.3)

    return ax

def plot_bound_ion(sim: Sim, site='O_Phosphate', ion='MG', window=100, ax=None):
    data = np.loadtxt(f'{sim.path}/dist_{ion}_{site}.xvg', comments=['@', '#'])
    t_range = data[:, 0] / 1000
    dist_ion_site = data[:, 1:]
    state = np.digitize(dist_ion_site, (0.25, 0.5))
    
    if ax is None:
        _, ax = plt.subplots()
    
    inner = (state == 0).sum(axis=1)
    outer = (state <= 1).sum(axis=1)

    ax.scatter(t_range, inner, s=0.1, alpha=0.2)
    ax.scatter(t_range, outer, s=0.1, alpha=0.2)
    ax.plot(*running_avg(t_range, inner, window=window), label='Inner')
    ax.plot(*running_avg(t_range, outer, window=window), label='Inner + Outer')
    ax.set_xlabel("Time [ns]")
    ax.set_ylabel(f"#{LABELS[ion]} ({LABELS[site]})")

    return ax

def running_avg(x, y, window=10):
    y = np.convolve(y, np.ones(window)/window, mode='valid')
    x = x[window//2:window//2+len(y)]
    return x, y