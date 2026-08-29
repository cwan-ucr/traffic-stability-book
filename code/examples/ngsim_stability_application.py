#!/usr/bin/env python3
"""D02: stability-oriented calibration, validation, and controller design.

The example uses the eight-vehicle US-101 platoon extracted by D01.  It
compares a conventional one-follower calibration with a platoon calibration
that also matches the disturbance-amplitude profile, validates both on a
held-out time interval, and then designs the smallest front-acceleration
feedforward gain that meets a declared ring-stability margin.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution
from scipy.signal import detrend, savgol_filter

from traffic_stability.idm import IDMParameters, idm_acceleration
from traffic_stability.linear import dominant_ring_mode, idm_derivatives
from traffic_stability.simulation import simulate_ring


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "derived" / "us101_selected_platoon.csv"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
VEHICLE_LENGTH_M = 4.8
TRAIN_END_S = 12.0
SEED = 20260828

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Hiragino Sans GB", "Songti SC", "Arial Unicode MS", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 9,
    }
)


def prepare_platoon(path: Path):
    frame = pd.read_csv(path)
    times = np.sort(frame["episode_time_s"].unique())
    orders = np.sort(frame["vehicle_order"].unique())
    vehicle_ids = [
        int(frame.loc[frame["vehicle_order"] == order, "vehicle_id"].iloc[0])
        for order in orders
    ]
    speed = np.column_stack(
        [
            frame[frame["vehicle_order"] == order]
            .sort_values("episode_time_s")["speed_mps"]
            .to_numpy()
            for order in orders
        ]
    )
    headway = np.column_stack(
        [
            frame[frame["vehicle_order"] == order]
            .sort_values("episode_time_s")["space_headway_m"]
            .to_numpy()
            for order in orders
        ]
    )
    speed = savgol_filter(speed, 21, 3, axis=0)
    headway = savgol_filter(headway, 21, 3, axis=0)
    net_gap = np.maximum(headway[:, 1:] - VEHICLE_LENGTH_M, 0.2)
    return times, vehicle_ids, speed, net_gap


def parameters_from_vector(vector):
    time_headway, minimum_gap, acceleration, deceleration = vector
    return IDMParameters(
        desired_speed=33.3,
        time_headway=float(time_headway),
        minimum_gap=float(minimum_gap),
        maximum_acceleration=float(acceleration),
        comfortable_deceleration=float(deceleration),
        acceleration_exponent=4.0,
        vehicle_length=VEHICLE_LENGTH_M,
    )


def interpolate_series(values: np.ndarray, index: int, fraction: float):
    if index >= len(values) - 1:
        return values[-1]
    return (1.0 - fraction) * values[index] + fraction * values[index + 1]


def chain_rhs(gap, follower_speed, leader_speed, leader_acceleration, parameters, gain):
    all_leader_speed = np.r_[leader_speed, follower_speed[:-1]]
    base = idm_acceleration(follower_speed, gap, all_leader_speed, parameters)
    acceleration = np.empty_like(base)
    upstream_acceleration = leader_acceleration
    for index in range(len(base)):
        acceleration[index] = base[index] + gain * upstream_acceleration
        upstream_acceleration = acceleration[index]
    return all_leader_speed - follower_speed, acceleration


def simulate_chain(
    time: np.ndarray,
    observed_speed: np.ndarray,
    observed_gap: np.ndarray,
    start_index: int,
    end_index: int,
    parameters: IDMParameters,
    gain: float = 0.0,
):
    dt = float(time[1] - time[0])
    follower_speed = observed_speed[start_index, 1:].copy()
    gap = observed_gap[start_index].copy()
    samples = end_index - start_index + 1
    simulated_speed = np.empty((samples, follower_speed.size))
    simulated_gap = np.empty_like(simulated_speed)
    simulated_speed[0], simulated_gap[0] = follower_speed, gap
    leader_acceleration_series = np.gradient(observed_speed[:, 0], dt)

    def rhs(local_gap, local_speed, absolute_index, fraction):
        leader_speed = float(interpolate_series(observed_speed[:, 0], absolute_index, fraction))
        leader_acceleration = float(
            interpolate_series(leader_acceleration_series, absolute_index, fraction)
        )
        return chain_rhs(
            local_gap,
            local_speed,
            leader_speed,
            leader_acceleration,
            parameters,
            gain,
        )

    for output_index, absolute_index in enumerate(range(start_index, end_index), start=1):
        k1s, k1v = rhs(gap, follower_speed, absolute_index, 0.0)
        k2s, k2v = rhs(
            gap + 0.5 * dt * k1s,
            follower_speed + 0.5 * dt * k1v,
            absolute_index,
            0.5,
        )
        k3s, k3v = rhs(
            gap + 0.5 * dt * k2s,
            follower_speed + 0.5 * dt * k2v,
            absolute_index,
            0.5,
        )
        k4s, k4v = rhs(
            gap + dt * k3s,
            follower_speed + dt * k3v,
            absolute_index,
            1.0,
        )
        gap = np.maximum(gap + dt * (k1s + 2 * k2s + 2 * k3s + k4s) / 6.0, 0.2)
        follower_speed = np.maximum(
            follower_speed + dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0,
            0.0,
        )
        simulated_speed[output_index], simulated_gap[output_index] = follower_speed, gap
    return simulated_speed, simulated_gap


def amplitude_profile(values: np.ndarray):
    return np.std(detrend(values, axis=0), axis=0)


def calibration_objective(vector, time, speed, gap, end_index, mode):
    parameters = parameters_from_vector(vector)
    simulated_speed, simulated_gap = simulate_chain(
        time, speed, gap, 0, end_index, parameters
    )
    observed_speed = speed[: end_index + 1, 1:]
    observed_gap = gap[: end_index + 1]
    if not np.all(np.isfinite(simulated_speed)) or np.min(simulated_gap) <= 0.2:
        return 1.0e6

    if mode == "single":
        speed_rmse = np.sqrt(np.mean((simulated_speed[:, 0] - observed_speed[:, 0]) ** 2))
        gap_rmse = np.sqrt(np.mean((simulated_gap[:, 0] - observed_gap[:, 0]) ** 2))
        return speed_rmse / 1.0 + gap_rmse / 5.0

    speed_rmse = np.sqrt(np.mean((simulated_speed - observed_speed) ** 2))
    gap_rmse = np.sqrt(np.mean((simulated_gap - observed_gap) ** 2))
    observed_amplitude = amplitude_profile(np.column_stack([speed[: end_index + 1, 0], observed_speed]))
    simulated_amplitude = amplitude_profile(
        np.column_stack([speed[: end_index + 1, 0], simulated_speed])
    )
    log_profile_rmse = np.sqrt(
        np.mean(
            (
                np.log(np.maximum(simulated_amplitude, 0.05))
                - np.log(np.maximum(observed_amplitude, 0.05))
            )
            ** 2
        )
    )
    return speed_rmse / 1.0 + gap_rmse / 5.0 + 1.25 * log_profile_rmse


def fit_parameters(time, speed, gap, mode):
    end_index = int(np.searchsorted(time, TRAIN_END_S))
    bounds = [(0.4, 2.8), (0.5, 6.0), (0.3, 3.0), (0.5, 4.5)]
    result = differential_evolution(
        calibration_objective,
        bounds,
        args=(time, speed, gap, end_index, mode),
        seed=SEED,
        maxiter=45,
        popsize=10,
        tol=2.0e-5,
        polish=True,
        workers=1,
        updating="immediate",
    )
    return parameters_from_vector(result.x), float(result.fun), end_index


def validation_metrics(
    time, speed, gap, start_index, parameters, gain=0.0, end_index=None
):
    if end_index is None:
        end_index = len(time) - 1
    simulated_speed, simulated_gap = simulate_chain(
        time, speed, gap, start_index, end_index, parameters, gain=gain
    )
    observed_speed = speed[start_index : end_index + 1, 1:]
    observed_gap = gap[start_index : end_index + 1]
    external_leader = speed[start_index : end_index + 1, 0]
    observed_profile = amplitude_profile(np.column_stack([external_leader, observed_speed]))
    simulated_profile = amplitude_profile(np.column_stack([external_leader, simulated_speed]))
    return {
        "speed_rmse_mps": float(np.sqrt(np.mean((simulated_speed - observed_speed) ** 2))),
        "gap_rmse_m": float(np.sqrt(np.mean((simulated_gap - observed_gap) ** 2))),
        "amplitude_profile_rmse_mps": float(
            np.sqrt(np.mean((simulated_profile - observed_profile) ** 2))
        ),
        "last_vehicle_gain_observed": float(observed_profile[-1] / observed_profile[0]),
        "last_vehicle_gain_simulated": float(simulated_profile[-1] / simulated_profile[0]),
        "speed": simulated_speed,
        "gap": simulated_gap,
        "observed_profile": observed_profile,
        "simulated_profile": simulated_profile,
    }


def choose_feedforward_gain(parameters, representative_speed):
    candidates = np.linspace(0.0, 0.50, 101)
    rates = np.array(
        [
            dominant_ring_mode(
                representative_speed,
                vehicle_count=60,
                parameters=parameters,
                leader_acceleration_gain=float(gain),
            )["growth_rate_per_s"]
            for gain in candidates
        ]
    )
    target = -0.002
    feasible = np.flatnonzero(rates <= target)
    selected_index = int(feasible[0]) if len(feasible) else int(np.argmin(rates))
    return float(candidates[selected_index]), float(rates[0]), float(rates[selected_index]), target


def json_parameters(parameters):
    return {
        "desired_speed_mps": parameters.desired_speed,
        "time_headway_s": parameters.time_headway,
        "minimum_gap_m": parameters.minimum_gap,
        "maximum_acceleration_mps2": parameters.maximum_acceleration,
        "comfortable_deceleration_mps2": parameters.comfortable_deceleration,
        "acceleration_exponent": parameters.acceleration_exponent,
        "vehicle_length_m": parameters.vehicle_length,
    }


def make_figure(
    time,
    vehicle_ids,
    observed_speed,
    train_end_index,
    single_parameters,
    fleet_parameters,
    single_validation,
    fleet_validation,
    optimized_validation,
    selected_gain,
    ring_baseline,
    ring_controlled,
    path,
):
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.4), constrained_layout=True)
    validation_time = time[train_end_index:] - time[train_end_index]

    axis = axes[0, 0]
    names = ["T [s]", "s0 [m]", "a [m/s^2]", "b [m/s^2]"]
    x = np.arange(len(names))
    axis.bar(x - 0.18, [
        single_parameters.time_headway,
        single_parameters.minimum_gap,
        single_parameters.maximum_acceleration,
        single_parameters.comfortable_deceleration,
    ], width=0.36, color="#c04a3a", label="单车标定")
    axis.bar(x + 0.18, [
        fleet_parameters.time_headway,
        fleet_parameters.minimum_gap,
        fleet_parameters.maximum_acceleration,
        fleet_parameters.comfortable_deceleration,
    ], width=0.36, color="#274c77", label="稳定性车队标定")
    axis.set_xticks(x, names)
    axis.set(title="D02(a) 两种标定得到的 IDM 参数", ylabel="参数值")
    axis.legend(frameon=False)

    axis = axes[0, 1]
    index = np.arange(len(vehicle_ids))
    axis.plot(index, fleet_validation["observed_profile"], "o-", color="0.25", label="留出数据")
    axis.plot(index, single_validation["simulated_profile"], "s--", color="#c04a3a", label="单车标定")
    axis.plot(index, fleet_validation["simulated_profile"], "d-", color="#274c77", label="车队标定")
    axis.set_xticks(index, [str(item) for item in vehicle_ids], rotation=35)
    axis.set(title="D02(b) 留出时段的逐车扰动幅值", xlabel="从领车到上游跟车（车辆 ID）", ylabel="速度标准差 [m/s]")
    axis.legend(frameon=False)

    axis = axes[1, 0]
    probe = -1
    axis.plot(validation_time, observed_speed[train_end_index:, -1], color="0.25", lw=2.0, label="实测")
    axis.plot(validation_time, single_validation["speed"][:, probe], color="#c04a3a", ls="--", label="单车标定模拟")
    axis.plot(validation_time, fleet_validation["speed"][:, probe], color="#274c77", label="车队标定模拟")
    axis.set(title=f"D02(c) 留出车辆 {vehicle_ids[-1]} 的速度验证", xlabel="留出时段 [s]", ylabel="速度 [m/s]")
    axis.legend(frameon=False)

    axis = axes[1, 1]
    axis.semilogy(ring_baseline["time"], ring_baseline["sigma_speed"], color="#c04a3a", label="标定基线")
    axis.semilogy(ring_controlled["time"], ring_controlled["sigma_speed"], color="#2f7d67", label=f"κ={selected_gain:.2f}")
    axis.set(title="D02(d) 稳定性优化后的环道验证", xlabel="时间 [s]", ylabel="速度标准差 [m/s]")
    axis.legend(frameon=False)

    for label, axis in zip(["(a)", "(b)", "(c)", "(d)"], axes.flat):
        axis.text(0.01, 0.98, label, transform=axis.transAxes, ha="left", va="top", fontweight="bold")
        axis.grid(alpha=0.18)
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def run():
    time, vehicle_ids, speed, gap = prepare_platoon(INPUT)
    single_parameters, single_objective, train_end_index = fit_parameters(
        time, speed, gap, "single"
    )
    fleet_parameters, fleet_objective, _ = fit_parameters(time, speed, gap, "fleet")

    single_validation = validation_metrics(
        time, speed, gap, train_end_index, single_parameters
    )
    fleet_validation = validation_metrics(
        time, speed, gap, train_end_index, fleet_parameters
    )
    single_training = validation_metrics(
        time, speed, gap, 0, single_parameters, end_index=train_end_index
    )
    fleet_training = validation_metrics(
        time, speed, gap, 0, fleet_parameters, end_index=train_end_index
    )
    observed_speed_interval = np.quantile(speed, [0.05, 0.50, 0.95])
    representative_speed = float(observed_speed_interval[1])
    design_speed = float(observed_speed_interval[0])
    selected_gain, baseline_rate, controlled_rate, target_rate = choose_feedforward_gain(
        fleet_parameters, design_speed
    )
    optimized_validation = validation_metrics(
        time, speed, gap, train_end_index, fleet_parameters, gain=selected_gain
    )
    ring_baseline = simulate_ring(
        design_speed,
        vehicle_count=60,
        amplitude=0.02,
        duration=600.0,
        dt=0.05,
        parameters=fleet_parameters,
    )
    ring_controlled = simulate_ring(
        design_speed,
        vehicle_count=60,
        amplitude=0.02,
        duration=600.0,
        dt=0.05,
        parameters=fleet_parameters,
        gain_schedule=lambda _time: selected_gain,
    )

    def public_metrics(metrics):
        return {
            key: value
            for key, value in metrics.items()
            if key not in {"speed", "gap", "observed_profile", "simulated_profile"}
        }

    margins = {
        "single_at_median": idm_derivatives(representative_speed, single_parameters)["margin"],
        "fleet_at_median": idm_derivatives(representative_speed, fleet_parameters)["margin"],
        "fleet_at_fifth_percentile_speed": idm_derivatives(design_speed, fleet_parameters)["margin"],
    }
    result = {
        "dataset": "NGSIM US-101 D01 selected eight-vehicle platoon",
        "vehicle_ids_leader_to_follower": vehicle_ids,
        "training_interval_s": [0.0, float(time[train_end_index])],
        "validation_interval_s": [float(time[train_end_index]), float(time[-1])],
        "observed_speed_quantiles_mps": observed_speed_interval.tolist(),
        "single_vehicle_calibration": {
            "parameters": json_parameters(single_parameters),
            "objective": single_objective,
            "training": public_metrics(single_training),
            "validation": public_metrics(single_validation),
        },
        "stability_aware_platoon_calibration": {
            "parameters": json_parameters(fleet_parameters),
            "objective": fleet_objective,
            "training": public_metrics(fleet_training),
            "validation": public_metrics(fleet_validation),
        },
        "stability_optimization": {
            "decision": "front-acceleration feedforward gain kappa",
            "search_interval": [0.0, 0.5],
            "target_growth_rate_per_s": target_rate,
            "design_speed_mps": design_speed,
            "selected_gain": selected_gain,
            "baseline_growth_rate_per_s": baseline_rate,
            "controlled_growth_rate_per_s": controlled_rate,
            "held_out_validation": public_metrics(optimized_validation),
            "ring_sigma_ratio_baseline": float(
                ring_baseline["sigma_speed"][-1] / ring_baseline["sigma_speed"][0]
            ),
            "ring_sigma_ratio_controlled": float(
                ring_controlled["sigma_speed"][-1] / ring_controlled["sigma_speed"][0]
            ),
        },
        "equilibrium_long_wave_margin": margins,
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    result_path = RESULTS / "ngsim_stability_application_d02.json"
    figure_path = FIGURES / "ngsim_stability_application_d02.png"
    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    make_figure(
        time,
        vehicle_ids,
        speed,
        train_end_index,
        single_parameters,
        fleet_parameters,
        single_validation,
        fleet_validation,
        optimized_validation,
        selected_gain,
        ring_baseline,
        ring_controlled,
        figure_path,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    print(f"figure: {figure_path}")
    print(f"results: {result_path}")


if __name__ == "__main__":
    run()
