"""Pure-Python D2Q9 Lattice Boltzmann Method implementation.

This provides a reference and fallback implementation that matches the C++ core semantics.
"""

import numpy as np
from typing import Tuple, Optional

# Lattice parameters
Q = 9
CX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=np.int32)
CY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=np.int32)
W = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36], dtype=np.float64)
CS2 = 1.0 / 3.0

class LatticeBoltzmann2D:
    """D2Q9 BGK Lattice Boltzmann solver."""
    
    def __init__(self, nx: int, ny: int, tau: float = 0.6):
        if nx <= 0 or ny <= 0:
            raise ValueError("Grid dimensions must be positive")
        if tau <= 0.5:
            raise ValueError("tau must be > 0.5 for stability")
        self.nx = nx
        self.ny = ny
        self.tau = tau
        self.omega = 1.0 / tau
        self.f = np.zeros((ny, nx, Q), dtype=np.float64)
        self.force_x = np.zeros((ny, nx), dtype=np.float64)
        self.force_y = np.zeros((ny, nx), dtype=np.float64)
    
    def initialize(self, rho0: float = 1.0, ux: float = 0.0, uy: float = 0.0) -> None:
        """Initialize distributions to equilibrium."""
        for q in range(Q):
            cu = CX[q] * ux + CY[q] * uy
            uu = ux*ux + uy*uy
            feq = W[q] * rho0 * (1.0 + cu/CS2 + cu*cu/(2*CS2*CS2) - uu/(2*CS2))
            self.f[:, :, q] = feq
    
    def density(self) -> np.ndarray:
        return np.sum(self.f, axis=2)
    
    def velocity(self) -> Tuple[np.ndarray, np.ndarray]:
        rho = self.density()
        rho_safe = np.maximum(rho, 1e-10)
        ux = np.sum(self.f * CX[None, None, :], axis=2) / rho_safe
        uy = np.sum(self.f * CY[None, None, :], axis=2) / rho_safe
        return ux, uy
    
    def collide(self) -> None:
        rho = self.density()
        ux, uy = self.velocity()
        # Guo forcing correction
        rho_safe = np.maximum(rho, 1e-10)
        ux = ux + 0.5 * self.force_x / rho_safe
        uy = uy + 0.5 * self.force_y / rho_safe
        
        for q in range(Q):
            cu = CX[q] * ux + CY[q] * uy
            uu = ux*ux + uy*uy
            feq = W[q] * rho * (1.0 + cu/CS2 + cu*cu/(2*CS2*CS2) - uu/(2*CS2))
            self.f[:, :, q] -= self.omega * (self.f[:, :, q] - feq)
        
        self.force_x.fill(0.0)
        self.force_y.fill(0.0)
    
    def stream(self) -> None:
        f_new = np.zeros_like(self.f)
        for q in range(Q):
            f_new[:, :, q] = np.roll(np.roll(self.f[:, :, q], CX[q], axis=1), CY[q], axis=0)
        self.f = f_new
    
    def step(self) -> None:
        self.collide()
        self.stream()
    
    def add_force(self, fx: np.ndarray, fy: np.ndarray) -> None:
        self.force_x += fx
        self.force_y += fy
    
    def vorticity(self) -> np.ndarray:
        ux, uy = self.velocity()
        # Central difference approximation
        dux_dy = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) * 0.5
        duy_dx = (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1)) * 0.5
        return duy_dx - dux_dy
