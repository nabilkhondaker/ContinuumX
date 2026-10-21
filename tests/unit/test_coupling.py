"""Interface coupling tests."""
import numpy as np
from continuumx.fluid.lbm import LatticeBoltzmann2D
from continuumx.structure.beam import FlexibleBeam
from continuumx.coupling.interface import ImmersedBoundaryCoupler

def test_interpolate_constant_field():
    fluid = LatticeBoltzmann2D(20, 10, tau=0.7)
    fluid.initialize(rho0=1.0, ux=0.1, uy=0.0)
    beam = FlexibleBeam(0.5, 5, EI=1.0, mass_per_length=1.0)
    # Place beam in domain
    beam.x = np.linspace(0.1, 0.4, 5)
    coupler = ImmersedBoundaryCoupler(fluid, beam, 2.0, 1.0)
    Ux, Uy = coupler.interpolate_velocity()
    # Should recover approximately the uniform velocity
    assert np.mean(np.abs(Ux - 0.1)) < 0.05

def test_residual_zero_when_matching():
    fluid = LatticeBoltzmann2D(10, 5)
    beam = FlexibleBeam(0.3, 4, 1.0, 1.0)
    coupler = ImmersedBoundaryCoupler(fluid, beam, 1.0, 0.5)
    u = np.array([0.0, 0.1, 0.1, 0.0])
    r = coupler.compute_coupling_residual(u, u)
    assert r < 1e-12
