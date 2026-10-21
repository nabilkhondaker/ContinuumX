# ContinuumX

**A coupled solid–fluid simulation engine for fluid–structure interaction**

📢 **Release Notice:** This repository contains the complete ContinuumX codebase, engineered between *July 20, 2026* and *October 21, 2026*. It was developed intermittently next to other research software, then assembled and published in full so the partitioned FSI loop, kernels, and diagnostics are available for public viewing and use.

**Author:** Nabil Khondaker

## Overview

ContinuumX is a research-oriented computational framework for simulating the interaction between incompressible fluid flow and deformable structures. The central problem is fluid–structure interaction (FSI): hydrodynamic forces deform a flexible body, the deformed geometry alters the flow, and the modified flow produces new forces.

The initial demonstration targets two-dimensional channel flow past a flexible cantilever. The fluid is modeled with a lattice Boltzmann method (D2Q9 BGK). The structure is a reduced-order cantilever beam with geometric nonlinearity. Coupling is performed through an immersed-boundary formulation that interpolates fluid velocity onto the structure and spreads structural forces onto the fluid grid.

This repository implements the numerical kernels, coupling algorithms, diagnostics, experiment configuration system, and verification tests. It is not a production multiphysics package. Its purpose is to demonstrate a correctly formulated partitioned FSI loop with conservation diagnostics and reproducible experiments.

## Why fluid–structure interaction is difficult

Independently correct fluid and structural solvers can produce an incorrect coupled result. Common failure modes include:

- Inconsistent force and velocity transfer operators that violate discrete power identity.
- Temporal mismatch between fluid and structural time levels.
- Added-mass instability when the structure is light relative to the displaced fluid.
- Uncontrolled energy growth from inconsistent interface work.
- Boundary-condition treatment that does not respect the moving interface.

ContinuumX addresses these by making the interface operators explicit, recording coupling residuals, and tracking mass, momentum, and energy diagnostics.

## Research questions

The implemented experiments are designed to investigate:

1. How does the choice of explicit versus strongly coupled iteration affect stability for a given density ratio?
2. Does the discrete delta function used for interpolation and spreading satisfy approximate power consistency?
3. Under what parameter regimes does explicit staggered coupling become unstable?
4. How sensitive is tip displacement to fluid grid resolution and structural node count?
5. Can energy diagnostics detect a deliberately incorrect force-transfer sign?

## Scope and non-goals

**Implemented:**

- D2Q9 BGK lattice Boltzmann fluid solver (pure Python reference implementation; C++ source provided).
- Reduced-order flexible cantilever with explicit time integration and simple geometric stiffening.
- Immersed-boundary velocity interpolation and force spreading using a 4-point Peskin delta.
- Explicit staggered coupling loop.
- Basic mass conservation and tip-displacement diagnostics.
- YAML experiment configuration and CLI runner.
- Unit tests for equilibrium preservation, mass conservation, and beam boundary conditions.

**Not implemented (explicit non-goals for this version):**

- Full geometrically exact finite-element beam or shell formulation.
- Multiple-relaxation-time or entropic collision operators.
- Strongly coupled fixed-point iteration with Aitken acceleration (interface is prepared; full residual loop is partial).
- Three-dimensional domains.
- Adaptive mesh refinement.
- GPU acceleration.
- Production-grade checkpoint/restart with HDF5.
- Quantitative comparison against published experimental or high-fidelity reference data with uncertainty quantification.

## Architecture

```
Python orchestration (continuumx package)
    ├── Experiment YAML → validated configuration
    ├── LatticeBoltzmann2D (fluid)
    ├── FlexibleBeam (structure)
    └── ImmersedBoundaryCoupler
            ├── interpolate_velocity
            └── spread_force
```

The C++ sources under `cpp/` provide a parallel D2Q9 implementation that can be compiled with pybind11 when a full build environment is available. The pure-Python path is the primary tested pathway in this distribution.

## Mathematical formulation

### Fluid

The isothermal incompressible Navier–Stokes equations are approximated by the lattice Boltzmann equation on the D2Q9 lattice. The BGK collision operator relaxes distributions toward the local equilibrium

\[
f_q^{\mathrm{eq}} = w_q \rho \left(1 + \frac{\mathbf{c}_q\cdot\mathbf{u}}{c_s^2} + \frac{(\mathbf{c}_q\cdot\mathbf{u})^2}{2c_s^4} - \frac{\mathbf{u}\cdot\mathbf{u}}{2c_s^2}\right)
\]

with \(c_s^2 = 1/3\). Kinematic viscosity is related to the relaxation time by \(\nu = c_s^2(\tau - 1/2)\Delta t\).

### Structure

The cantilever is represented by a chain of nodes with transverse displacement \(y_i\). Elastic restoring force is obtained from a finite-difference approximation of bending moment divergence. A simplified geometric stiffening term accounts for moderate rotation. Time integration uses explicit central difference with optional linear damping.

### Interface

Velocity at structural node \(\mathbf{X}_k\) is interpolated as

\[
\mathbf{U}_k = \sum_j \mathbf{u}(\mathbf{x}_j) \delta_h(\mathbf{x}_j - \mathbf{X}_k) \Delta V_j
\]

Force is spread by the adjoint relation using the same discrete delta. The 4-point Peskin kernel is used.

## Explicit versus strong coupling

The current runner implements an explicit staggered scheme: fluid is advanced, interface velocity is interpolated, a force estimate is formed, the force is spread, and the structure is advanced. A residual between fluid-interpolated and structural velocity is computable but is not yet used to drive a fixed-point iteration inside each time step. Strong coupling can be added by wrapping the structural update and force evaluation in an inner loop until the residual falls below a tolerance.

## Conservation and energy

Fluid mass is recovered as the sum of distribution functions. On a periodic domain the implementation conserves mass to machine precision (verified by unit test). Structural kinetic energy is available from the nodal velocities and lumped masses. Interface work can be estimated from the inner product of interpolated velocity and structural force; a full discrete power-identity test is present in the coupler but requires consistent scaling of the force density.

## Verification

Unit tests cover:

- Preservation of uniform equilibrium under zero velocity.
- Mass conservation under periodic advection.
- Enforcement of clamped boundary condition.
- Positive tip deflection under tip loading.

These tests establish that the kernels behave as intended. They do not constitute physical validation against experimental data.

## Main numerical case study

The configuration `experiments/cantilever/vortex_shedding.yaml` defines a 200×50 lattice, a 21-node cantilever of length 0.4, and an inlet velocity of 0.05. Running

```bash
python -m continuumx.cli simulate experiments/cantilever/vortex_shedding.yaml
```

produces a time history of tip displacement. The magnitude of the response depends on the simplified force model currently used (relative-velocity loading). A production implementation would replace this with the momentum-exchange force obtained from the lattice populations at the interface.

## Limitations

- The force applied to the structure is a simplified relative-velocity model rather than the full momentum-exchange force from the LBM populations.
- Strong coupling iteration is not yet active in the time loop.
- No quantitative comparison to a published FSI benchmark has been performed.
- The beam model is a reduced-order finite-difference representation, not a full nonlinear finite-element formulation.
- Lattice Boltzmann stability requires \(\tau > 0.5\) and low Mach number; these constraints are enforced at construction but not dynamically monitored in every experiment.
- OpenMP parallelization and C++ bindings require a build environment that is not guaranteed in all deployment contexts.

## Installation

```bash
cd ContinuumX
pip install -e .
```

Requires Python ≥ 3.10, NumPy, PyYAML, and pytest (for tests).

## Quick start

```bash
python -m continuumx.cli simulate experiments/cantilever/vortex_shedding.yaml
pytest tests/unit/test_core.py -q
```

## Citation

```
@software{khondaker2026continuumx,
  author = {Khondaker, Nabil},
  title = {ContinuumX: A Coupled Solid–Fluid Simulation Engine},
  year = {2026},
  url = {https://github.com/nabilkhondaker/ContinuumX}
}
```

## Final assessment

ContinuumX demonstrates a correctly structured partitioned FSI loop with an immersed-boundary interface, mass-conserving lattice Boltzmann fluid, and a flexible cantilever that responds to hydrodynamic loading. The implementation is sufficient to explore coupling strategy, grid sensitivity, and diagnostic detection of interface errors. It is not yet a validated predictive tool for quantitative engineering design. Further work is required on the momentum-exchange force, strong coupling convergence, and comparison against independent reference solutions.
