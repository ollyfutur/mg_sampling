"""Physical constants at the simulation conditions."""

import numpy as np
from scipy.constants import R

π = np.pi

RT = R * 300 * 1e-3  # kJ/mol at 300 K

"""
Coefficient to multiply the coordination number to get approx. integer values.
"""
SWITCH_COEFF = 1.039661