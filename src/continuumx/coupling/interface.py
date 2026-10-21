"""Immersed boundary coupling between LBM fluid and beam structure."""

import numpy as np
from typing import Tuple
from ..fluid.lbm import LatticeBoltzmann2D
from ..structure.beam import FlexibleBeam

class ImmersedBoundaryCoupler:
    """Force spreading and velocity interpolation for FSI."""
    
    def __init__(self, fluid: LatticeBoltzmann2D, structure: FlexibleBeam,
                 domain_length: float, domain_height: float):
        self.fluid = fluid
        self.structure = structure
        self.Lx = domain_length
        self.Ly = domain_height
        self.dx = domain_length / fluid.nx
        self.dy = domain_height / fluid.ny
    
    def _discrete_delta(self, r: float) -> float:
        """Peskin 4-point discrete delta function."""
        r = abs(r)
        if r < 1.0:
            return (3.0 - 2.0*r + np.sqrt(1.0 + 4.0*r - 4.0*r*r)) / 8.0
        elif r < 2.0:
            return (5.0 - 2.0*r - np.sqrt(-7.0 + 12.0*r - 4.0*r*r)) / 8.0
        return 0.0
    
    def interpolate_velocity(self) -> Tuple[np.ndarray, np.ndarray]:
        """Interpolate fluid velocity to structural nodes."""
        ux, uy = self.fluid.velocity()
        n = self.structure.n_nodes
        Ux = np.zeros(n)
        Uy = np.zeros(n)
        
        positions = self.structure.positions()
        for k in range(n):
            X, Y = positions[k]
            i_c = int(X / self.dx)
            j_c = int(Y / self.dy)
            
            for di in range(-2, 3):
                for dj in range(-2, 3):
                    i = i_c + di
                    j = j_c + dj
                    if 0 <= i < self.fluid.nx and 0 <= j < self.fluid.ny:
                        r_x = (X - i * self.dx) / self.dx
                        r_y = (Y - j * self.dy) / self.dy
                        phi = self._discrete_delta(r_x) * self._discrete_delta(r_y)
                        Ux[k] += ux[j, i] * phi
                        Uy[k] += uy[j, i] * phi
        
        return Ux, Uy
    
    def spread_force(self, force_x: np.ndarray, force_y: np.ndarray) -> None:
        """Spread structural forces onto fluid grid."""
        fx_grid = np.zeros((self.fluid.ny, self.fluid.nx))
        fy_grid = np.zeros((self.fluid.ny, self.fluid.nx))
        
        positions = self.structure.positions()
        ds = self.structure.dx
        
        for k in range(self.structure.n_nodes):
            X, Y = positions[k]
            i_c = int(X / self.dx)
            j_c = int(Y / self.dy)
            Fk_x = force_x[k]
            Fk_y = force_y[k]
            
            for di in range(-2, 3):
                for dj in range(-2, 3):
                    i = i_c + di
                    j = j_c + dj
                    if 0 <= i < self.fluid.nx and 0 <= j < self.fluid.ny:
                        r_x = (X - i * self.dx) / self.dx
                        r_y = (Y - j * self.dy) / self.dy
                        phi = self._discrete_delta(r_x) * self._discrete_delta(r_y)
                        fx_grid[j, i] += Fk_x * phi / (self.dx * self.dy)
                        fy_grid[j, i] += Fk_y * phi / (self.dx * self.dy)
        
        self.fluid.add_force(fx_grid, fy_grid)
    
    def compute_coupling_residual(self, U_fluid: np.ndarray, U_struct: np.ndarray,
                                  scale: float = 1e-6) -> float:
        """Relative residual between interface velocities."""
        diff = np.linalg.norm(U_fluid - U_struct)
        denom = max(np.linalg.norm(U_struct), scale)
        return diff / denom
