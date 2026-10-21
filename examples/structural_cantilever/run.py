"""Structural-only cantilever under tip load."""
import numpy as np
from continuumx.structure.beam import FlexibleBeam
b = FlexibleBeam(1.0, 21, EI=0.05, mass_per_length=0.5, damping=0.05)
f = np.zeros(21); f[-1] = 0.005
for _ in range(200):
    b.step(0.01, f)
print("Tip displacement:", b.tip_displacement())
