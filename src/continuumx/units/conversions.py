"""Lattice to physical unit conversions (documented placeholders)."""
def viscosity_from_tau(tau, dx=1.0, dt=1.0):
    cs2 = 1.0 / 3.0
    return cs2 * (tau - 0.5) * dt
