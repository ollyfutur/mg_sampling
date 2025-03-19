from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent.parent

FORCEFIELDS = [
    'amber_tip3p',
    'mMg_tip3p',
    'nMg_tip3p',
    'mMg_opc',
]

LABELS = {
    'nop': r'$n_\mathrm{O_P}$',
    'nob': r'$n_\mathrm{O_B}$',
    'now': r'$n_\mathrm{O_W}$',
    'op': r'$\mathrm{O_P}$',
    'ob': r'$\mathrm{O_B}$',
    'amber_tip3p': "Amber + TIP3P",
    'mMg_tip3p': "microMg + TIP3P",
    'nMg_tip3p': "nanoMg + TIP3P", 
    'mMg_opc': "microMg + OPC",
    'dop1': r'$d(\mathrm{Mg},\mathrm{O_P})$ [nm]',
    'bias.bias': r'Bias [kJ/mol]',
    'X': r'$X$',
    'Y': r'$Y$',
}
