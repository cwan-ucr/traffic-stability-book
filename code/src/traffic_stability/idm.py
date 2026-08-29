"""Intelligent Driver Model equations used by the book and examples."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class IDMParameters:
    """Baseline IDM parameters in SI units."""

    desired_speed: float = 33.3
    time_headway: float = 1.5
    minimum_gap: float = 2.0
    maximum_acceleration: float = 1.0
    comfortable_deceleration: float = 1.5
    acceleration_exponent: float = 4.0
    vehicle_length: float = 5.0


def equilibrium_gap(speed: float | np.ndarray, parameters: IDMParameters = IDMParameters()):
    """Return the net equilibrium gap for a uniform IDM flow."""

    value = np.asarray(speed, dtype=float)
    ratio = np.clip(value / parameters.desired_speed, 0.0, 1.0 - 1.0e-12)
    gap = (
        parameters.minimum_gap + parameters.time_headway * value
    ) / np.sqrt(1.0 - ratio ** parameters.acceleration_exponent)
    return float(gap) if np.ndim(value) == 0 else gap


def equilibrium_speed(net_gap: float, parameters: IDMParameters = IDMParameters()) -> float:
    """Invert the uniform-flow equilibrium relation by bisection."""

    if net_gap < parameters.minimum_gap:
        return 0.0
    lo, hi = 0.0, (1.0 - 1.0e-10) * parameters.desired_speed
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if equilibrium_gap(mid, parameters) < net_gap:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def idm_acceleration(
    speed: np.ndarray,
    gap: np.ndarray,
    leader_speed: np.ndarray,
    parameters: IDMParameters = IDMParameters(),
) -> np.ndarray:
    """Evaluate IDM acceleration for vectorized follower states."""

    p = parameters
    speed = np.asarray(speed, dtype=float)
    gap = np.asarray(gap, dtype=float)
    leader_speed = np.asarray(leader_speed, dtype=float)
    closing_speed = speed - leader_speed
    desired_gap = (
        p.minimum_gap
        + p.time_headway * np.maximum(speed, 0.0)
        + speed * closing_speed
        / (2.0 * np.sqrt(p.maximum_acceleration * p.comfortable_deceleration))
    )
    return p.maximum_acceleration * (
        1.0
        - (np.maximum(speed, 0.0) / p.desired_speed) ** p.acceleration_exponent
        - (desired_gap / np.maximum(gap, 0.1)) ** 2
    )
