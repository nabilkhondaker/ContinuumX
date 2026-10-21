"""Unit tests for ContinuumX core components."""

import pytest
import numpy as np
from continuumx.fluid.lbm import LatticeBoltzmann2D
from continuumx.structure.beam import FlexibleBeam


def test_lbm_equilibrium_preservation():
    lbm = LatticeBoltzmann2D(10, 10, tau=0.6)
    lbm.initialize(rho0=1.0, ux=0.0, uy=0.0)
    rho_before = lbm.density()
    lbm.step()
    rho_after = lbm.density()
    # Equilibrium should be approximately preserved for zero velocity
    assert np.allclose(rho_before, rho_after, rtol=1e-5)


def test_lbm_mass_conservation_periodic():
    lbm = LatticeBoltzmann2D(20, 20, tau=0.8)
    lbm.initialize(rho0=1.0, ux=0.02, uy=0.0)
    mass0 = np.sum(lbm.density())
    for _ in range(10):
        lbm.step()
    mass1 = np.sum(lbm.density())
    assert abs(mass1 - mass0) / mass0 < 1e-6


def test_beam_clamped_equilibrium():
    beam = FlexibleBeam(length=1.0, n_nodes=11, EI=1.0, mass_per_length=1.0)
    beam.apply_clamped_bc()
    assert beam.y[0] == 0.0
    assert beam.v[0] == 0.0


def test_beam_tip_response_to_load():
    beam = FlexibleBeam(length=1.0, n_nodes=21, EI=0.1, mass_per_length=0.5, damping=0.1)
    force = np.zeros(21)
    force[-1] = 0.01
    for _ in range(100):
        beam.step(0.01, force)
    assert beam.tip_displacement() > 0.0
