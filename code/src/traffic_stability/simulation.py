"""Small, explicit numerical solvers used by the reproducibility examples."""

from __future__ import annotations

import numpy as np

from .idm import IDMParameters, equilibrium_gap, idm_acceleration


def _cyclic_front_acceleration(base: np.ndarray, gain: float) -> np.ndarray:
    """Solve ``a_n - gain*a_(n+1) = base_n`` exactly on a ring."""

    if abs(gain) < 1.0e-15:
        return base
    count = len(base)
    powers = gain ** np.arange(count)
    result = np.empty_like(base)
    denominator = 1.0 - gain**count
    if abs(denominator) < 1.0e-10:
        raise ValueError("cyclic acceleration coupling is singular")
    for n in range(count):
        result[n] = np.dot(powers, base[(n + np.arange(count)) % count]) / denominator
    return result


def _rhs(
    position: np.ndarray,
    speed: np.ndarray,
    ring_length: float,
    parameters: IDMParameters,
    leader_acceleration_gain: float,
) -> tuple[np.ndarray, np.ndarray]:
    leader_position = np.roll(position, -1)
    gap = leader_position - position - parameters.vehicle_length
    gap[-1] += ring_length
    leader_speed = np.roll(speed, -1)
    base = idm_acceleration(speed, gap, leader_speed, parameters)
    acceleration = _cyclic_front_acceleration(base, leader_acceleration_gain)
    return speed, acceleration


def _rk4_step(
    position: np.ndarray,
    speed: np.ndarray,
    ring_length: float,
    dt: float,
    parameters: IDMParameters,
    leader_acceleration_gain: float,
) -> tuple[np.ndarray, np.ndarray]:
    args = (ring_length, parameters, leader_acceleration_gain)
    k1x, k1v = _rhs(position, speed, *args)
    k2x, k2v = _rhs(position + 0.5 * dt * k1x, speed + 0.5 * dt * k1v, *args)
    k3x, k3v = _rhs(position + 0.5 * dt * k2x, speed + 0.5 * dt * k2v, *args)
    k4x, k4v = _rhs(position + dt * k3x, speed + dt * k3v, *args)
    next_position = position + dt * (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0
    next_speed = speed + dt * (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0
    return next_position, np.maximum(next_speed, 0.0)


def simulate_ring(
    equilibrium_speed_value: float,
    vehicle_count: int = 30,
    mode: int = 1,
    amplitude: float = 0.002,
    duration: float = 600.0,
    dt: float = 0.05,
    sample_dt: float = 0.5,
    parameters: IDMParameters = IDMParameters(),
    gain_schedule=None,
) -> dict[str, np.ndarray | float]:
    """Integrate a periodic IDM ring and record a Fourier-mode amplitude.

    ``gain_schedule(t)`` may change the front-acceleration gain without
    resetting the state, which is the basis of mini-example M04.
    """

    gap = equilibrium_gap(equilibrium_speed_value, parameters)
    headway = gap + parameters.vehicle_length
    ring_length = vehicle_count * headway
    vehicle_index = np.arange(vehicle_count)
    position = headway * vehicle_index.astype(float)
    phase = 2.0 * np.pi * mode * vehicle_index / vehicle_count
    speed = equilibrium_speed_value + amplitude * np.sin(phase)
    schedule = gain_schedule or (lambda _t: 0.0)

    sample_every = max(1, round(sample_dt / dt))
    steps = round(duration / dt)
    samples = steps // sample_every + 1
    time = np.empty(samples)
    sigma = np.empty(samples)
    mode_amplitude = np.empty(samples, dtype=complex)
    gains = np.empty(samples)
    probe_speed = np.empty(samples)
    probe_gap = np.empty(samples)

    out = 0
    for step in range(steps + 1):
        current_time = step * dt
        if step % sample_every == 0:
            velocity_deviation = speed - np.mean(speed)
            time[out] = current_time
            sigma[out] = np.std(speed)
            mode_amplitude[out] = (
                2.0
                / vehicle_count
                * np.sum(velocity_deviation * np.exp(-1j * phase))
            )
            gains[out] = schedule(current_time)
            probe_speed[out] = speed[vehicle_count // 4]
            gaps = np.roll(position, -1) - position - parameters.vehicle_length
            gaps[-1] += ring_length
            probe_gap[out] = gaps[vehicle_count // 4]
            out += 1
        if step == steps:
            break
        position, speed = _rk4_step(
            position,
            speed,
            ring_length,
            dt,
            parameters,
            float(schedule(current_time)),
        )
        if not np.all(np.isfinite(speed)):
            raise RuntimeError(f"non-finite state at t={current_time:.3f} s")
        if np.min(speed) < -1.0e-9:
            raise RuntimeError("negative speed encountered")

    return {
        "time": time,
        "sigma_speed": sigma,
        "mode_amplitude": mode_amplitude,
        "gain": gains,
        "probe_speed": probe_speed,
        "probe_gap": probe_gap,
        "ring_length": float(ring_length),
        "equilibrium_gap": float(gap),
    }


def fit_growth_rate(time: np.ndarray, amplitude: np.ndarray, start: float, end: float) -> float:
    """Fit ``log(amplitude)`` over a declared linear time window."""

    mask = (time >= start) & (time <= end) & (amplitude > 1.0e-12)
    if np.count_nonzero(mask) < 4:
        raise ValueError("growth-rate window has too few valid samples")
    return float(np.polyfit(time[mask], np.log(amplitude[mask]), 1)[0])
