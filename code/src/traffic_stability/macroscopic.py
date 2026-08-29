"""A compact IDM-consistent LWR solver for teaching example M06."""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from .idm import IDMParameters, equilibrium_gap


def idm_flow(density_veh_per_m: np.ndarray, parameters: IDMParameters = IDMParameters()) -> np.ndarray:
    """IDM equilibrium flow ``q=rho*V(rho)`` in vehicles per second."""

    density = np.asarray(density_veh_per_m, dtype=float)
    output = np.zeros_like(density)
    jam_density = 1.0 / (parameters.minimum_gap + parameters.vehicle_length)
    valid = (density > 0.0) & (density < jam_density)
    if np.any(valid):
        target_gap = 1.0 / density[valid] - parameters.vehicle_length
        lo = np.zeros_like(target_gap)
        hi = np.full_like(target_gap, (1.0 - 1.0e-10) * parameters.desired_speed)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            below = equilibrium_gap(mid, parameters) < target_gap
            lo = np.where(below, mid, lo)
            hi = np.where(below, hi, mid)
        output[valid] = density[valid] * 0.5 * (lo + hi)
    return output


@lru_cache(maxsize=8)
def capacity_state(parameters: IDMParameters = IDMParameters()) -> tuple[float, float]:
    """Return capacity density and flow from a dense numerical scan."""

    jam_density = 1.0 / (parameters.minimum_gap + parameters.vehicle_length)
    grid = np.linspace(1.0e-7, jam_density, 5000)
    flow = idm_flow(grid, parameters)
    peak = int(np.argmax(flow))
    return float(grid[peak]), float(flow[peak])


def godunov_flux(
    left_density: np.ndarray,
    right_density: np.ndarray,
    parameters: IDMParameters = IDMParameters(),
) -> np.ndarray:
    """Supply-demand Godunov flux for the unimodal IDM fundamental diagram."""

    critical_density, capacity = capacity_state(parameters)
    left_flow = idm_flow(left_density, parameters)
    right_flow = idm_flow(right_density, parameters)
    demand = np.where(left_density <= critical_density, left_flow, capacity)
    supply = np.where(right_density <= critical_density, capacity, right_flow)
    return np.minimum(demand, supply)


def solve_lwr_riemann(
    left_density_veh_per_km: float = 20.0,
    right_density_veh_per_km: float = 60.0,
    length_m: float = 4000.0,
    cells: int = 160,
    duration_s: float = 60.0,
    dt_s: float = 0.5,
    parameters: IDMParameters = IDMParameters(),
) -> dict[str, np.ndarray | float]:
    """Solve a fixed-boundary Riemann problem with finite volumes."""

    dx = length_m / cells
    x = (np.arange(cells) + 0.5) * dx
    midpoint = 0.5 * length_m
    left = left_density_veh_per_km / 1000.0
    right = right_density_veh_per_km / 1000.0
    density = np.where(x < midpoint, left, right)
    initial = density.copy()
    steps = round(duration_s / dt_s)
    for _ in range(steps):
        extended = np.concatenate(([left], density, [right]))
        flux = godunov_flux(extended[:-1], extended[1:], parameters)
        density = density - dt_s / dx * (flux[1:] - flux[:-1])
        if np.min(density) < -1.0e-12:
            raise RuntimeError("negative LWR density; check the CFL number")

    q_left, q_right = idm_flow(np.array([left, right]), parameters)
    exact_speed = float((q_right - q_left) / (right - left))
    middle_density = 0.5 * (left + right)
    crossing = int(np.argmin(np.abs(density - middle_density)))
    measured_position = float(x[crossing])
    measured_speed = (measured_position - midpoint) / duration_s
    return {
        "x_m": x,
        "initial_density_veh_per_km": 1000.0 * initial,
        "final_density_veh_per_km": 1000.0 * density,
        "exact_shock_speed_mps": exact_speed,
        "measured_shock_speed_mps": measured_speed,
        "shock_position_m": measured_position,
        "dx_m": float(dx),
        "dt_s": float(dt_s),
    }
