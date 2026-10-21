"""Grid and time-step convergence utilities."""
import numpy as np

def relative_error(coarse, fine):
    return np.linalg.norm(fine - coarse) / max(np.linalg.norm(fine), 1e-12)
