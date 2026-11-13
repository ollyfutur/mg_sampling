"""Labels, colors and settings shared by the paper notebooks."""

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

labelcolors = {v: simcolors[k] for k, v in simlabels.items() if k in simcolors}

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

# Equilibration discarded at the start of each simulation [ns]
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
