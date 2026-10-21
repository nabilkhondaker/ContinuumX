"""Mass and momentum diagnostics."""
import numpy as np

def fluid_mass(density_field):
    return float(np.sum(density_field))

def relative_mass_drift(m0, m):
    return abs(m - m0) / max(abs(m0), 1e-12)
