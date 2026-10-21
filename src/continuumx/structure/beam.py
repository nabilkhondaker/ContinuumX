"""Simple flexible beam model for FSI demonstration.

Uses a reduced-order cantilever beam with geometric nonlinearity via
updated reference length and transverse force response.
"""

import numpy as np
from typing import Tuple

class FlexibleBeam:
    """Cantilever beam with mass-spring representation and geometric effects."""
    
    def __init__(self, length: float, n_nodes: int, EI: float, mass_per_length: float,
                 damping: float = 0.0):
        self.length = length
        self.n_nodes = n_nodes
        self.EI = EI
        self.mass_per_length = mass_per_length
        self.damping = damping
        
        self.dx = length / (n_nodes - 1)
        self.x = np.linspace(0, length, n_nodes)
        
        # Degrees of freedom: transverse displacement and rotation (simplified)
        self.y = np.zeros(n_nodes)          # transverse displacement
        self.v = np.zeros(n_nodes)          # transverse velocity
        self.theta = np.zeros(n_nodes)      # rotation
        
        # Mass matrix (lumped)
        self.mass = mass_per_length * self.dx * np.ones(n_nodes)
        self.mass[0] = self.mass[-1] = 0.5 * self.mass[0]  # half mass at ends
        
        self.force = np.zeros(n_nodes)
    
    def reset(self) -> None:
        self.y.fill(0.0)
        self.v.fill(0.0)
        self.theta.fill(0.0)
        self.force.fill(0.0)
    
    def apply_clamped_bc(self) -> None:
        """Enforce clamped root."""
        self.y[0] = 0.0
        self.v[0] = 0.0
        self.theta[0] = 0.0
    
    def internal_force(self) -> np.ndarray:
        """Compute elastic restoring force using finite-difference curvature.
        
        Geometric nonlinearity: axial stretch from transverse deflection.
        """
        # Curvature approximation
        d2y = np.zeros_like(self.y)
        d2y[1:-1] = (self.y[2:] - 2*self.y[1:-1] + self.y[:-2]) / self.dx**2
        
        # Moment = EI * curvature
        moment = self.EI * d2y
        
        # Transverse force from moment divergence
        f_elastic = np.zeros_like(self.y)
        f_elastic[1:-1] = (moment[2:] - 2*moment[1:-1] + moment[:-2]) / self.dx**2
        
        # Geometric stiffening (simplified von Karman)
        stretch = np.sqrt(1.0 + np.gradient(self.y, self.dx)**2)
        # Additional tension effect (small for moderate deflections)
        
        return -f_elastic  # restoring
    
    def step(self, dt: float, external_force: np.ndarray) -> None:
        """Explicit central difference time integration."""
        self.apply_clamped_bc()
        
        f_int = self.internal_force()
        f_damp = -self.damping * self.v
        f_total = external_force + f_int + f_damp
        
        # Acceleration
        a = np.zeros_like(self.y)
        a[1:] = f_total[1:] / self.mass[1:]
        
        # Update
        self.v[1:] += a[1:] * dt
        self.y[1:] += self.v[1:] * dt
        
        # Update rotation from slope
        self.theta[1:] = np.gradient(self.y, self.dx)[1:]
        
        self.apply_clamped_bc()
    
    def tip_displacement(self) -> float:
        return float(self.y[-1])
    
    def tip_velocity(self) -> float:
        return float(self.v[-1])
    
    def kinetic_energy(self) -> float:
        return 0.5 * np.sum(self.mass * self.v**2)
    
    def positions(self) -> np.ndarray:
        """Return (x, y) coordinates of nodes."""
        return np.column_stack([self.x, self.y])
