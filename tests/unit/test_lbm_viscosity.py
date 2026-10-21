"""Viscosity-related checks for LBM."""
import numpy as np
from continuumx.fluid.lbm import LatticeBoltzmann2D

def test_tau_stability():
    import pytest
    with pytest.raises(ValueError):
        LatticeBoltzmann2D(10, 10, tau=0.4)

def test_density_recovery():
    lbm = LatticeBoltzmann2D(8, 8, tau=0.8)
    lbm.initialize(rho0=1.2, ux=0.0, uy=0.0)
    rho = lbm.density()
    assert np.allclose(rho, 1.2, atol=1e-10)
