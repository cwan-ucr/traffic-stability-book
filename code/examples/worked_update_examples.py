#!/usr/bin/env python3
"""W01--W03: hand-checkable microscopic and macroscopic update examples.

The examples deliberately use only three vehicles or three/four cells.  They
print every state needed to reproduce the first update and save the first two
steps as JSON for the worked examples in the book.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from traffic_stability.idm import IDMParameters, idm_acceleration


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def ring_rhs(position: np.ndarray, speed: np.ndarray, length: float, p: IDMParameters):
    leader_position = np.roll(position, -1)
    gap = leader_position - position - p.vehicle_length
    gap[-1] += length
    leader_speed = np.roll(speed, -1)
    acceleration = idm_acceleration(speed, gap, leader_speed, p)
    return speed.copy(), acceleration, gap


def rk4_with_stages(
    position: np.ndarray,
    speed: np.ndarray,
    length: float,
    dt: float,
    p: IDMParameters,
):
    k1x, k1v, gap1 = ring_rhs(position, speed, length, p)
    x2, v2 = position + 0.5 * dt * k1x, speed + 0.5 * dt * k1v
    k2x, k2v, gap2 = ring_rhs(x2, v2, length, p)
    x3, v3 = position + 0.5 * dt * k2x, speed + 0.5 * dt * k2v
    k3x, k3v, gap3 = ring_rhs(x3, v3, length, p)
    x4, v4 = position + dt * k3x, speed + dt * k3v
    k4x, k4v, gap4 = ring_rhs(x4, v4, length, p)
    next_position = position + dt * (k1x + 2 * k2x + 2 * k3x + k4x) / 6
    next_speed = speed + dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6
    return next_position, next_speed, {
        "K1": {"gap_m": gap1, "dxdt_mps": k1x, "dvdt_mps2": k1v},
        "K2": {"gap_m": gap2, "dxdt_mps": k2x, "dvdt_mps2": k2v},
        "K3": {"gap_m": gap3, "dxdt_mps": k3x, "dvdt_mps2": k3v},
        "K4": {"gap_m": gap4, "dxdt_mps": k4x, "dvdt_mps2": k4v},
    }


def microscopic_example():
    p = IDMParameters()
    length, dt = 75.0, 0.2
    position = np.array([0.0, 25.0, 50.0])
    speed = np.array([10.0, 9.5, 10.5])
    states = [{"step": 0, "time_s": 0.0, "position_m": position, "speed_mps": speed}]
    stages_first = None
    for step in range(2):
        position, speed, stages = rk4_with_stages(position, speed, length, dt, p)
        if step == 0:
            stages_first = stages
        states.append(
            {
                "step": step + 1,
                "time_s": (step + 1) * dt,
                "position_m": position,
                "speed_mps": speed,
            }
        )
    return {
        "vehicle_count": 3,
        "ring_length_m": length,
        "vehicle_length_m": p.vehicle_length,
        "dt_s": dt,
        "initial_leader_map": {"0": 1, "1": 2, "2": 0},
        "first_step_stages": stages_first,
        "states": states,
    }


def triangular_flow(density: np.ndarray):
    free_speed, congestion_speed, jam_density = 30.0, 5.0, 0.15
    critical = congestion_speed * jam_density / (free_speed + congestion_speed)
    capacity = free_speed * critical
    density = np.asarray(density)
    flow = np.minimum(free_speed * density, congestion_speed * (jam_density - density))
    return flow, critical, capacity


def triangular_godunov(left: np.ndarray, right: np.ndarray):
    left_flow, critical, capacity = triangular_flow(left)
    right_flow, _, _ = triangular_flow(right)
    demand = np.where(left <= critical, left_flow, capacity)
    supply = np.where(right <= critical, capacity, right_flow)
    return np.minimum(demand, supply)


def lwr_example():
    dx, dt = 500.0, 5.0
    left_ghost, right_ghost = 0.010, 0.050
    density = np.array([0.010, 0.020, 0.050, 0.050])
    records = [{"step": 0, "time_s": 0.0, "density_veh_per_km": 1000 * density}]
    first_flux = None
    for step in range(2):
        extended = np.r_[left_ghost, density, right_ghost]
        flux = triangular_godunov(extended[:-1], extended[1:])
        if step == 0:
            first_flux = flux.copy()
        density = density - dt / dx * (flux[1:] - flux[:-1])
        records.append(
            {
                "step": step + 1,
                "time_s": (step + 1) * dt,
                "density_veh_per_km": 1000 * density,
            }
        )
    _, critical, capacity = triangular_flow(np.array([0.0]))
    return {
        "cells": 4,
        "dx_m": dx,
        "dt_s": dt,
        "free_speed_mps": 30.0,
        "congestion_wave_speed_mps": 5.0,
        "jam_density_veh_per_km": 150.0,
        "critical_density_veh_per_km": 1000 * critical,
        "capacity_veh_per_h": 3600 * capacity,
        "boundary_density_veh_per_km": [10.0, 50.0],
        "first_interface_flux_veh_per_s": first_flux,
        "states": records,
    }


def arz_primitive_to_conserved(density: np.ndarray, speed: np.ndarray, pressure_slope: float):
    generalized_speed = speed + pressure_slope * density
    return np.column_stack([density, density * generalized_speed])


def arz_conserved_to_primitive(state: np.ndarray, pressure_slope: float):
    density = state[:, 0]
    generalized_speed = state[:, 1] / density
    speed = generalized_speed - pressure_slope * density
    return density, speed, generalized_speed


def arz_flux(state: np.ndarray, pressure_slope: float):
    density, speed, generalized_speed = arz_conserved_to_primitive(state, pressure_slope)
    return np.column_stack([density * speed, density * speed * generalized_speed])


def arz_rusanov_interfaces(state: np.ndarray, pressure_slope: float):
    flux = arz_flux(state, pressure_slope)
    interfaces = []
    for left_index in range(len(state)):
        right_index = (left_index + 1) % len(state)
        rho_l, v_l, _ = arz_conserved_to_primitive(state[[left_index]], pressure_slope)
        rho_r, v_r, _ = arz_conserved_to_primitive(state[[right_index]], pressure_slope)
        eigenvalues = [
            v_l[0],
            v_l[0] - pressure_slope * rho_l[0],
            v_r[0],
            v_r[0] - pressure_slope * rho_r[0],
        ]
        alpha = max(abs(item) for item in eigenvalues)
        numerical_flux = (
            0.5 * (flux[left_index] + flux[right_index])
            - 0.5 * alpha * (state[right_index] - state[left_index])
        )
        interfaces.append(numerical_flux)
    return np.asarray(interfaces)


def arz_equilibrium_speed(density: np.ndarray):
    return 30.0 * (1.0 - density / 0.15)


def arz_example():
    dx, dt, tau, pressure_slope = 50.0, 0.5, 10.0, 200.0
    density = np.array([0.020, 0.030, 0.040])
    speed = np.array([20.0, 16.0, 12.0])
    state = arz_primitive_to_conserved(density, speed, pressure_slope)
    records = []
    first_flux = None
    for step in range(3):
        rho, vel, generalized = arz_conserved_to_primitive(state, pressure_slope)
        records.append(
            {
                "step": step,
                "time_s": step * dt,
                "density_veh_per_km": 1000 * rho,
                "speed_mps": vel,
                "generalized_speed_mps": generalized,
                "conserved_rho_w": state[:, 1],
                "vehicles_in_periodic_domain": float(np.sum(rho) * dx),
            }
        )
        if step == 2:
            break
        interface_flux = arz_rusanov_interfaces(state, pressure_slope)
        if step == 0:
            first_flux = interface_flux.copy()
        source = np.zeros_like(state)
        source[:, 1] = rho * (arz_equilibrium_speed(rho) - vel) / tau
        state = (
            state
            - dt / dx * (interface_flux - np.roll(interface_flux, 1, axis=0))
            + dt * source
        )
    return {
        "cells": 3,
        "periodic_length_m": 3 * dx,
        "dx_m": dx,
        "dt_s": dt,
        "relaxation_time_s": tau,
        "pressure": "p(rho)=200 rho",
        "equilibrium_speed": "V(rho)=30(1-rho/0.15)",
        "first_interface_flux": first_flux,
        "states": records,
    }


def serializable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [serializable(item) for item in value]
    return value


def run():
    result = serializable(
        {
            "W01_three_vehicle_idm_rk4": microscopic_example(),
            "W02_four_cell_lwr": lwr_example(),
            "W03_three_cell_arz": arz_example(),
        }
    )
    RESULTS.mkdir(parents=True, exist_ok=True)
    destination = RESULTS / "worked_step_examples.json"
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    print(f"results: {destination}")


if __name__ == "__main__":
    run()
