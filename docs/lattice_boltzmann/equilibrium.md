# Lattice Boltzmann Mathematics

## D2Q9 Equilibrium

The equilibrium distribution for the isothermal BGK model is

\[
f_q^{\mathrm{eq}}(\rho, \mathbf{u}) = w_q \rho \left[
1 + \frac{\mathbf{c}_q \cdot \mathbf{u}}{c_s^2}
+ \frac{(\mathbf{c}_q \cdot \mathbf{u})^2}{2 c_s^4}
- \frac{\mathbf{u}\cdot\mathbf{u}}{2 c_s^2}
\right]
\]

with \(c_s^2 = 1/3\).

## Viscosity

\[
\nu = c_s^2 (\tau - \tfrac12) \Delta t
\]

Stability requires \(\tau > 1/2\).

## Forcing

A simple Guo-type forcing is used: the macroscopic velocity is corrected by \(\frac12 \mathbf{F}/\rho\) before equilibrium evaluation, and the force is added to the momentum.
