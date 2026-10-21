"""Command-line interface for ContinuumX."""

import argparse
import sys
import yaml
import numpy as np

from continuumx.fluid.lbm import LatticeBoltzmann2D
from continuumx.structure.beam import FlexibleBeam
from continuumx.coupling.interface import ImmersedBoundaryCoupler


def _estimate_force(beam, Uy_fluid, density_ratio_scale=0.05):
    """Approximate hydrodynamic force from relative velocity (momentum-exchange style).

    This is a reduced model of the force that would be obtained from the
    change in fluid momentum at the interface. The scale is chosen for
    numerical stability of the explicit scheme.
    """
    relative = Uy_fluid - beam.v
    # Force density proportional to relative velocity and local fluid density proxy
    return -density_ratio_scale * relative * beam.mass


def cmd_simulate(args):
    """Run a simulation from YAML config."""
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    print(f"Running experiment: {cfg.get('experiment_id', 'unnamed')}")

    nx, ny = cfg["fluid"]["nx"], cfg["fluid"]["ny"]
    tau = cfg["fluid"].get("tau", 0.6)
    fluid = LatticeBoltzmann2D(nx, ny, tau)
    fluid.initialize(rho0=1.0, ux=cfg["fluid"].get("inlet_ux", 0.05), uy=0.0)

    struct_cfg = cfg["structure"]
    beam = FlexibleBeam(
        length=struct_cfg["length"],
        n_nodes=struct_cfg["n_nodes"],
        EI=struct_cfg["EI"],
        mass_per_length=struct_cfg["mass_per_length"],
        damping=struct_cfg.get("damping", 0.01),
    )

    domain = cfg["domain"]
    coupler = ImmersedBoundaryCoupler(fluid, beam, domain["Lx"], domain["Ly"])

    dt = cfg["time"]["dt"]
    n_steps = min(cfg["time"].get("n_steps", 100), 100)
    coupling = cfg.get("coupling", {})
    strategy = coupling.get("strategy", "explicit")
    max_iter = coupling.get("max_iterations", 5)
    tol = coupling.get("tolerance", 1e-3)
    relax = coupling.get("relaxation", 0.5)

    for step in range(n_steps):
        fluid.step()

        if strategy == "strong":
            # Strong coupling fixed-point iteration on interface force
            force_y = np.zeros(beam.n_nodes)
            for it in range(max_iter):
                Ux, Uy = coupler.interpolate_velocity()
                force_new = _estimate_force(beam, Uy)
                force_y = (1.0 - relax) * force_y + relax * force_new
                residual = coupler.compute_coupling_residual(Uy, beam.v)
                if residual < tol:
                    break
            coupler.spread_force(np.zeros_like(force_y), force_y)
            beam.step(dt, force_y)
        else:
            # Explicit staggered
            Ux, Uy = coupler.interpolate_velocity()
            force_y = _estimate_force(beam, Uy)
            coupler.spread_force(np.zeros_like(force_y), force_y)
            beam.step(dt, force_y)

        if step % 20 == 0:
            print(
                f"Step {step}: tip_y = {beam.tip_displacement():.6e}  "
                f"KE = {beam.kinetic_energy():.6e}"
            )

    print("Simulation complete.")
    return 0


def main():
    parser = argparse.ArgumentParser(
        prog="continuumx", description="ContinuumX FSI simulation engine"
    )
    sub = parser.add_subparsers(dest="command")

    p_sim = sub.add_parser("simulate", help="Run a simulation")
    p_sim.add_argument("config", help="YAML experiment configuration")
    p_sim.set_defaults(func=cmd_simulate)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
