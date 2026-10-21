# Known Limitations

## Addressed in this version
- Explicit and strong coupling strategies are both implemented in the time loop.
- Coupling residual is computed and available for convergence checks.
- Additional experiments and tests cover mass conservation, interpolation, and parameter sweeps.

## Remaining limitations
- Hydrodynamic force is a reduced relative-velocity model scaled for stability. A full momentum-exchange force reconstructed from lattice populations is not yet implemented.
- The structural model remains a reduced-order finite-difference cantilever, not a geometrically exact finite-element beam.
- No quantitative comparison against a published experimental or high-fidelity FSI benchmark has been executed.
- C++ bindings require a manual compilation step; the pure-Python path is the tested default.
- Energy diagnostics track structural kinetic energy but do not yet close a full discrete power balance across the interface.
