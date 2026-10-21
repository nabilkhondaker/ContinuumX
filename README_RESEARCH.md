# ContinuumX Research Notes

## Motivation
Partitioned fluid-structure interaction remains sensitive to interface transfer operators and coupling strategy. ContinuumX provides a transparent implementation of a D2Q9 lattice Boltzmann fluid coupled to a flexible cantilever via an immersed-boundary method.

## Research questions addressed
- Effect of explicit versus strongly coupled iteration on tip response.
- Ability of residual monitoring to detect interface mismatch.
- Mass conservation properties of the lattice Boltzmann core.

## Limitations
See docs/limitations.md. The force model is a reduced relative-velocity approximation rather than a full momentum-exchange reconstruction from distribution functions.
