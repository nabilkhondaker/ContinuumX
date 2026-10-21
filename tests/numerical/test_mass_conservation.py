"""Mass conservation under periodic conditions."""
import numpy as np
from continuumx.fluid.lbm import LatticeBoltzmann2D

def test_periodic_mass_long():
    lbm = LatticeBoltzmann2D(32, 32, tau=0.9)
    lbm.initialize(1.0, 0.03, 0.01)
    m0 = np.sum(lbm.density())
    for _ in range(20):
        lbm.step()
    m1 = np.sum(lbm.density())
    assert abs(m1 - m0) / m0 < 1e-5
