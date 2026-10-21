"""Simple field export helpers."""
import numpy as np

def save_velocity_csv(ux, uy, path):
    np.savetxt(path, np.column_stack([ux.ravel(), uy.ravel()]), header="ux uy")
