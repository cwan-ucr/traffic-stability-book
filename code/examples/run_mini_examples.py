#!/usr/bin/env python3
"""Run the six compact numerical examples M01--M06 and make one figure."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np

from traffic_stability.linear import dominant_ring_mode, ring_eigenvalues
from traffic_stability.macroscopic import solve_lwr_riemann
from traffic_stability.simulation import fit_growth_rate, simulate_ring


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Hiragino Sans GB", "Songti SC", "Arial Unicode MS", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 9,
    }
)


def rk4_local(dt: float = 0.05, duration: float = 40.0):
    """M01: integrate y''+0.2y'+0.25y=0 with y(0)=1,y'(0)=0."""

    fs, fv = 0.25, -0.2
    time = np.arange(0.0, duration + 0.5 * dt, dt)
    state = np.array([1.0, 0.0])
    numerical = np.empty_like(time)

    def rhs(value):
        return np.array([value[1], fv * value[1] - fs * value[0]])

    for index, _ in enumerate(time):
        numerical[index] = state[0]
        if index == len(time) - 1:
            break
        k1 = rhs(state)
        k2 = rhs(state + 0.5 * dt * k1)
        k3 = rhs(state + 0.5 * dt * k2)
        k4 = rhs(state + dt * k3)
        state += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    omega = np.sqrt(fs - 0.25 * fv**2)
    analytic = np.exp(0.5 * fv * time) * (
        np.cos(omega * time) - 0.5 * fv / omega * np.sin(omega * time)
    )
    roots = np.roots([1.0, -fv, fs])
    return time, numerical, analytic, roots


def run():
    metrics: dict[str, dict] = {}
    fig, axes = plt.subplots(3, 2, figsize=(12.0, 11.2), constrained_layout=True)

    # M01: local underdamped recovery.
    time, numerical, analytic, roots = rk4_local()
    axes[0, 0].plot(time, analytic, color="#c04a3a", lw=2.0, label="解析解")
    axes[0, 0].plot(time[::20], numerical[::20], "o", ms=2.7, color="#274c77", label="RK4")
    axes[0, 0].axhline(0.0, color="0.75", lw=0.8)
    axes[0, 0].set(title="M01  局部欠阻尼恢复", xlabel="时间 t [s]", ylabel="间距扰动 δs [m]")
    axes[0, 0].legend(frameon=False)
    metrics["M01"] = {
        "equation": "delta_s_ddot + 0.2 delta_s_dot + 0.25 delta_s = 0",
        "roots": [[float(root.real), float(root.imag)] for root in roots],
        "period_s": float(2.0 * np.pi / abs(roots[0].imag)),
        "decay_time_s": float(-1.0 / roots[0].real),
        "first_undershoot_m": float(np.min(numerical[(time > 0.0) & (time < 10.0)])),
        "rk4_max_abs_error_m": float(np.max(np.abs(numerical - analytic))),
    }

    # M02: a tiny ring demonstrates discrete wavelength sampling.
    counts = np.array([12, 20, 30])
    finite = [dominant_ring_mode(8.0, int(count)) for count in counts]
    growth = np.array([item["growth_rate_per_s"] for item in finite])
    colors = np.where(growth >= 0.0, "#c04a3a", "#2f7d67")
    axes[0, 1].bar(counts.astype(str), growth, color=colors)
    axes[0, 1].axhline(0.0, color="0.2", lw=0.8)
    axes[0, 1].set(title="M02  有限环道只允许离散波长", xlabel="车辆数 N", ylabel="最右增长率 [1/s]")
    metrics["M02"] = {f"N{count}": result for count, result in zip(counts, finite)}

    # M03: nonlinear IDM, same perturbation, stable and unstable speeds.
    stable = simulate_ring(25.0, duration=600.0, dt=0.05)
    unstable = simulate_ring(8.0, duration=600.0, dt=0.05)
    stable_amp = np.abs(stable["mode_amplitude"])
    unstable_amp = np.abs(unstable["mode_amplitude"])
    axes[1, 0].semilogy(stable["time"], stable_amp, color="#2f7d67", label="25 m/s：稳定")
    axes[1, 0].semilogy(unstable["time"], unstable_amp, color="#c04a3a", label="8 m/s：失稳")
    axes[1, 0].set(title="M03  原始 IDM 稳定/失稳对照", xlabel="时间 t [s]", ylabel="模态幅值 [m/s]")
    axes[1, 0].legend(frameon=False)
    metrics["M03"] = {}
    for name, speed, output, amp in [
        ("stable", 25.0, stable, stable_amp),
        ("unstable", 8.0, unstable, unstable_amp),
    ]:
        theory = dominant_ring_mode(speed, 30)
        fit_end = 250.0 if speed == 25.0 else 500.0
        metrics["M03"][name] = {
            "speed_mps": speed,
            "theoretical_growth_per_s": theory["growth_rate_per_s"],
            "fitted_growth_per_s": fit_growth_rate(output["time"], amp, 30.0, fit_end),
            "amplitude_ratio": float(amp[-1] / amp[0]),
        }

    # M04: preserve the state and ramp up front-acceleration feedback.
    def gain_schedule(t):
        if t < 300.0:
            return 0.0
        if t < 360.0:
            return 0.25 * (t - 300.0) / 60.0
        return 0.25

    intervention = simulate_ring(
        8.0,
        amplitude=0.2,
        duration=900.0,
        dt=0.05,
        gain_schedule=gain_schedule,
    )
    axis = axes[1, 1]
    axis.plot(intervention["time"], intervention["sigma_speed"], color="#274c77", label="速度标准差")
    axis.axvline(300.0, color="#c04a3a", ls="--", lw=1.0)
    axis.set(title="M04  仿真中途开启前车加速度反馈", xlabel="时间 t [s]", ylabel=r"$\sigma_v$ [m/s]")
    second = axis.twinx()
    second.plot(intervention["time"], intervention["gain"], color="#d19a31", lw=1.4, label="κ")
    second.set_ylabel("前馈增益 κ")
    before = (intervention["time"] >= 240.0) & (intervention["time"] <= 300.0)
    after = (intervention["time"] >= 780.0) & (intervention["time"] <= 900.0)
    metrics["M04"] = {
        "activation_time_s": 300.0,
        "ramp_duration_s": 60.0,
        "final_gain": 0.25,
        "uncontrolled_growth_per_s": dominant_ring_mode(8.0, 30, leader_acceleration_gain=0.0)["growth_rate_per_s"],
        "controlled_growth_per_s": dominant_ring_mode(8.0, 30, leader_acceleration_gain=0.25)["growth_rate_per_s"],
        "mean_sigma_before_mps": float(np.mean(intervention["sigma_speed"][before])),
        "mean_sigma_after_mps": float(np.mean(intervention["sigma_speed"][after])),
        "peak_sigma_mps": float(np.max(intervention["sigma_speed"])),
    }

    # M05: time-step convergence of a fitted physical growth rate.
    time_steps = np.array([0.40, 0.20, 0.10])
    fitted = []
    for dt in time_steps:
        output = simulate_ring(8.0, amplitude=2.0e-4, duration=600.0, dt=float(dt), sample_dt=0.5)
        fitted.append(fit_growth_rate(output["time"], np.abs(output["mode_amplitude"]), 150.0, 550.0))
    fitted = np.asarray(fitted)
    theory = dominant_ring_mode(8.0, 30)["growth_rate_per_s"]
    relative_error = 100.0 * np.abs(fitted - theory) / abs(theory)
    axes[2, 0].semilogy(time_steps, relative_error, "o-", color="#274c77", label="RK4 拟合误差")
    axes[2, 0].invert_xaxis()
    axes[2, 0].set(title="M05  时间步长收敛检查", xlabel="时间步长 Δt [s]", ylabel="增长率相对误差 [%]")
    axes[2, 0].legend(frameon=False)
    metrics["M05"] = {
        "theoretical_growth_per_s": theory,
        "runs": [
            {
                "dt_s": float(dt),
                "fitted_growth_per_s": float(value),
                "relative_error_percent": float(100.0 * abs(value - theory) / abs(theory)),
            }
            for dt, value in zip(time_steps, fitted)
        ],
    }

    # M06: a conservative LWR shock with an IDM equilibrium diagram.
    lwr = solve_lwr_riemann()
    axes[2, 1].plot(lwr["x_m"] / 1000.0, lwr["initial_density_veh_per_km"], color="0.6", ls="--", label="t=0 s")
    axes[2, 1].plot(lwr["x_m"] / 1000.0, lwr["final_density_veh_per_km"], color="#274c77", label="t=60 s")
    exact_position = 2.0 + lwr["exact_shock_speed_mps"] * 60.0 / 1000.0
    axes[2, 1].axvline(exact_position, color="#c04a3a", ls=":", label="Rankine–Hugoniot")
    axes[2, 1].set(title="M06  IDM 基本图下的 LWR 激波", xlabel="道路位置 x [km]", ylabel="密度 [veh/km]")
    axes[2, 1].legend(frameon=False)
    metrics["M06"] = {
        "left_density_veh_per_km": 20.0,
        "right_density_veh_per_km": 60.0,
        "exact_shock_speed_mps": lwr["exact_shock_speed_mps"],
        "measured_shock_speed_mps": lwr["measured_shock_speed_mps"],
        "absolute_speed_error_mps": abs(lwr["measured_shock_speed_mps"] - lwr["exact_shock_speed_mps"]),
        "dx_m": lwr["dx_m"],
        "dt_s": lwr["dt_s"],
    }

    for label, axis in zip(["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"], axes.flat):
        axis.text(0.01, 0.98, label, transform=axis.transAxes, ha="left", va="top", fontweight="bold")
        axis.grid(alpha=0.18)

    figure_path = FIGURES / "mini_examples_m01_m06.png"
    fig.savefig(figure_path, dpi=240, bbox_inches="tight")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    panel_axes = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1], axes[2, 0], axes[2, 1]]
    for number, panel_axis in enumerate(panel_axes, start=1):
        boxes = [panel_axis.get_tightbbox(renderer)]
        if number == 4:
            boxes.append(second.get_tightbbox(renderer))
        panel_box = Bbox.union(boxes)
        panel_box_inches = panel_box.transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(
            FIGURES / f"mini_example_m{number:02d}.png",
            dpi=240,
            bbox_inches=panel_box_inches,
            pad_inches=0.04,
        )
    plt.close(fig)
    metrics_path = RESULTS / "mini_examples_metrics.json"
    metrics_path.write_text(json.dumps(metrics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))
    print(f"figure: {figure_path}")
    print(f"metrics: {metrics_path}")


if __name__ == "__main__":
    run()
