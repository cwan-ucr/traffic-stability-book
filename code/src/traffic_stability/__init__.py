"""Reproducible teaching experiments for traffic-flow stability."""

from .idm import IDMParameters, equilibrium_gap, equilibrium_speed, idm_acceleration
from .linear import idm_derivatives, ring_eigenvalues, transfer_gain
from .simulation import simulate_ring

__all__ = [
    "IDMParameters",
    "equilibrium_gap",
    "equilibrium_speed",
    "idm_acceleration",
    "idm_derivatives",
    "ring_eigenvalues",
    "transfer_gain",
    "simulate_ring",
]

__version__ = "0.1.0"
