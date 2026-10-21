"""ContinuumX: Coupled solid-fluid simulation engine.

Author: Nabil Khondaker
"""

__version__ = "0.1.0"
__author__ = "Nabil Khondaker"

from .fluid.lbm import LatticeBoltzmann2D
from .structure.beam import FlexibleBeam
from .coupling.interface import ImmersedBoundaryCoupler

__all__ = ["LatticeBoltzmann2D", "FlexibleBeam", "ImmersedBoundaryCoupler"]
