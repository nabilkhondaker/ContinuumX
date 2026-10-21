"""Minimal fluid-only example."""
from continuumx.fluid.lbm import LatticeBoltzmann2D
lbm = LatticeBoltzmann2D(40, 20, tau=0.7)
lbm.initialize(1.0, 0.02, 0.0)
for _ in range(50):
    lbm.step()
print("Final mass:", lbm.density().sum())
