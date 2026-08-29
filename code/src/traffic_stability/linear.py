"""Linear stability quantities for the baseline IDM."""

from __future__ import annotations

import numpy as np

from .idm import IDMParameters, equilibrium_gap


def idm_derivatives(speed: float, parameters: IDMParameters = IDMParameters()) -> dict[str, float]:
    """Return equilibrium gap, partial derivatives, and long-wave margin.

    The sign convention is ``Delta v = v_n - v_{n+1}``.  Consequently,
    ``f_v`` is the derivative with respect to the follower speed after the
    acceleration law is written as ``F(s_n, v_n, v_{n+1})`` and ``f_vl`` is
    the leader-speed derivative.
    """

    p = parameters
    gap = equilibrium_gap(speed, p)
    z = np.sqrt(1.0 - (speed / p.desired_speed) ** p.acceleration_exponent)
    f_s = 2.0 * p.maximum_acceleration * z**2 / gap
    free_speed_term = (
        p.maximum_acceleration
        * p.acceleration_exponent
        / p.desired_speed
        * (speed / p.desired_speed) ** (p.acceleration_exponent - 1.0)
    )
    headway_term = 2.0 * p.maximum_acceleration * z * p.time_headway / gap
    approach_term = (
        np.sqrt(p.maximum_acceleration / p.comfortable_deceleration)
        * z
        * speed
        / gap
    )
    f_v = -(free_speed_term + headway_term + approach_term)
    f_vl = approach_term
    margin = 0.5 * (f_v**2 - f_vl**2) - f_s
    equilibrium_slope = -f_s / (f_v + f_vl)
    return {
        "speed": float(speed),
        "gap": float(gap),
        "f_s": float(f_s),
        "f_v": float(f_v),
        "f_vl": float(f_vl),
        "margin": float(margin),
        "equilibrium_slope": float(equilibrium_slope),
    }


def ring_eigenvalues(
    speed: float,
    vehicle_count: int,
    mode: int,
    parameters: IDMParameters = IDMParameters(),
    leader_acceleration_gain: float = 0.0,
) -> np.ndarray:
    """Return the two temporal roots of one discrete ring mode.

    The optional gain corresponds to ``a_n = F_n + kappa a_{n+1}``.
    """

    d = idm_derivatives(speed, parameters)
    phase = np.exp(1j * 2.0 * np.pi * mode / vehicle_count)
    coefficients = [
        1.0 - leader_acceleration_gain * phase,
        -(d["f_v"] + d["f_vl"] * phase),
        -d["f_s"] * (phase - 1.0),
    ]
    return np.roots(coefficients)


def dominant_ring_mode(
    speed: float,
    vehicle_count: int,
    parameters: IDMParameters = IDMParameters(),
    leader_acceleration_gain: float = 0.0,
) -> dict[str, float]:
    """Find the non-trivial mode with the largest real eigenvalue."""

    best_mode, best_root = 1, None
    for mode in range(1, vehicle_count):
        roots = ring_eigenvalues(
            speed, vehicle_count, mode, parameters, leader_acceleration_gain
        )
        root = roots[np.argmax(roots.real)]
        if best_root is None or root.real > best_root.real:
            best_mode, best_root = mode, root
    assert best_root is not None
    represented_mode = min(best_mode, vehicle_count - best_mode)
    return {
        "mode": int(represented_mode),
        "growth_rate_per_s": float(best_root.real),
        "angular_frequency_rad_per_s": float(abs(best_root.imag)),
        "wavelength_vehicles": float(vehicle_count / represented_mode),
    }


def transfer_gain(
    speed: float,
    angular_frequency: float,
    parameters: IDMParameters = IDMParameters(),
) -> float:
    """Magnitude of the baseline vehicle-to-vehicle transfer function."""

    d = idm_derivatives(speed, parameters)
    s = 1j * angular_frequency
    gain = (d["f_vl"] * s + d["f_s"]) / (
        s**2 - d["f_v"] * s + d["f_s"]
    )
    return float(abs(gain))
