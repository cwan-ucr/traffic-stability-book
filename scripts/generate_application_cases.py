#!/usr/bin/env python3
"""Generate two concrete stability-management cases for the lab and manuscript."""

from __future__ import annotations

import base64
import json
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ASSET_OPEN = ROOT / "book" / "source" / "assets" / "application_open_platoon.png"
ASSET_VSL = ROOT / "book" / "source" / "assets" / "application_moving_bottleneck.png"
METRICS = ROOT / "experiments" / "metrics" / "application-cases-metrics.json"
COMPACT = ROOT / "lab" / "application-cases-compact.json"
LAB_DATA = ROOT / "lab" / "application-cases-data.js"


@dataclass(frozen=True)
class IDM:
    v0: float = 33.3
    T: float = 1.5
    s0: float = 2.0
    a: float = 1.0
    b: float = 1.5
    delta: float = 4.0


P = IDM()
# High-demand motorway case: the statutory limit remains 120 km/h, but the
# short desired headway intentionally places the dense equilibrium near the
# long-wave stability boundary.
HIGHWAY_P = replace(P, v0=120.0 / 3.6, T=0.7)


@dataclass(frozen=True)
class WuEV:
    """Parameters reported in Wu et al. (TRD 34, 2015), Table 1."""

    mass: float = 1266.0
    rolling: float = 0.006
    aero_k: float = 1.30
    motor_k: float = 10.08
    motor_r: float = 0.11
    wheel_r: float = 0.50
    grade_rad: float = 0.0


EV = WuEV()


def equilibrium_gap(v: float, p: IDM = P) -> float:
    return (p.s0 + p.T * v) / np.sqrt(1.0 - (v / p.v0) ** p.delta)


def idm_acc(v: float, gap: float, leader_v: float, p: IDM = P) -> float:
    dv = v - leader_v
    desired = p.s0 + max(0.0, v * p.T + v * dv / (2.0 * np.sqrt(p.a * p.b)))
    return p.a * (1.0 - (max(v, 0.0) / p.v0) ** p.delta - (desired / max(gap, 0.1)) ** 2)


def wu_ev_power(v: np.ndarray, a: np.ndarray, ev: WuEV = EV) -> np.ndarray:
    """Wu et al. Eq. (11), battery-side instantaneous power in watts.

    Negative values represent regenerative charging.  The road is level in the
    two application cases, so the grade term is zero but retained explicitly.
    """
    g = 9.81
    resistance = ev.aero_k * v**2 + ev.rolling * ev.mass * g + ev.mass * g * np.sin(ev.grade_rad)
    force = ev.mass * a + resistance
    copper_loss = ev.motor_r * ev.wheel_r**2 / ev.motor_k**2 * force**2
    return copper_loss + resistance * v + ev.mass * a * v


def ev_energy_metrics(time: np.ndarray, v: np.ndarray, a: np.ndarray, reference_speed: float) -> dict:
    power = wu_ev_power(v, a)
    energy_each = np.trapezoid(power, time, axis=0) / 3.6e6
    distance_each = np.trapezoid(v, time, axis=0) / 1000.0
    ref_power = float(wu_ev_power(np.array(reference_speed), np.array(0.0)))
    ref_energy = ref_power * (time[-1] - time[0]) / 3.6e6
    return {
        "mean_ev_energy_kwh": float(np.mean(energy_each)),
        "ev_energy_intensity_kwh_per_100km": float(100.0 * np.sum(energy_each) / np.sum(distance_each)),
        "mean_excess_ev_energy_wh": float(1000.0 * np.mean(energy_each - ref_energy)),
        "min_regen_power_kw": float(np.min(power) / 1000.0),
    }


def reconstruct_positions(time: np.ndarray, speed: np.ndarray, gaps: np.ndarray, vehicle_length: float = 5.0) -> np.ndarray:
    """Reconstruct absolute open-road trajectories from leader motion and gaps."""
    x = np.zeros_like(speed)
    dt = np.diff(time)
    x[1:, 0] = np.cumsum(0.5 * (speed[:-1, 0] + speed[1:, 0]) * dt)
    for i in range(1, speed.shape[1]):
        x[:, i] = x[:, i - 1] - gaps[:, i - 1] - vehicle_length
    x -= x[0, -1]
    return x


def equilibrium_speed_for_gap(gap: float, p: IDM) -> float:
    lo, hi = 0.0, 0.999 * p.v0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if idm_acc(mid, gap, mid, p) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def idm_long_wave_margin(gap: float, p: IDM) -> tuple[float, float]:
    """Return equilibrium speed and Phi=1/2(f_v^2-f_vl^2)-f_s."""
    v = equilibrium_speed_for_gap(gap, p)
    hs, hv = 1.0e-3, 1.0e-4
    fs = (idm_acc(v, gap + hs, v, p) - idm_acc(v, gap - hs, v, p)) / (2.0 * hs)
    fv = (idm_acc(v + hv, gap, v, p) - idm_acc(v - hv, gap, v, p)) / (2.0 * hv)
    fvl = (idm_acc(v, gap, v + hv, p) - idm_acc(v, gap, v - hv, p)) / (2.0 * hv)
    return v, 0.5 * (fv**2 - fvl**2) - fs


def desired_speed_parameter_for_equilibrium(gap: float, speed: float, p: IDM = P) -> float:
    """Invert the uniform IDM equilibrium for the internal desired-speed parameter."""
    remaining = 1.0 - ((p.s0 + p.T * speed) / gap) ** 2
    if remaining <= 0.0:
        raise ValueError("requested equilibrium speed is infeasible at this gap")
    return speed / remaining ** (1.0 / p.delta)


def recovery_time(time: np.ndarray, signal: np.ndarray, start: float, threshold: float, hold: float) -> float | None:
    dt = float(np.median(np.diff(time)))
    k = max(1, round(hold / dt))
    good = (signal < threshold) & (time >= start)
    run = np.convolve(good.astype(int), np.ones(k, dtype=int), mode="valid")
    ids = np.flatnonzero(run == k)
    return float(time[ids[0]]) if ids.size else None


def simulate_open_platoon(controlled: bool, *, duration: float = 1200.0, dt: float = 0.1) -> dict:
    n = 60
    ve = 8.0
    gap0 = equilibrium_gap(ve)
    force_start, period, cycles, amp = 100.0, 48.0, 5, 0.6
    force_end = force_start + period * cycles

    def leader(t: float) -> tuple[float, float]:
        if force_start <= t <= force_end:
            phase = 2.0 * np.pi * (t - force_start) / period
            return ve + amp * np.sin(phase), amp * (2.0 * np.pi / period) * np.cos(phase)
        return ve, 0.0

    speed = np.full(n, ve, dtype=float)
    gaps = np.full(n - 1, gap0, dtype=float)

    def rhs(t: float, g: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        lv, la = leader(t)
        vv = v.copy()
        vv[0] = lv
        base = np.zeros(n)
        base[0] = la
        for i in range(1, n):
            base[i] = idm_acc(vv[i], g[i - 1], vv[i - 1])
        acc = base.copy()
        if controlled:
            # Cooperative AV: three-predecessor acceleration feed-forward plus
            # rear-vehicle speed feedback.  All communicated terms are bounded.
            kappa, alpha, gamma, eta = 0.34, 0.52, 0.045, 0.055
            for i in range(1, n):
                ff = 0.0
                align = 0.0
                for j in range(1, min(3, i) + 1):
                    w = alpha ** (j - 1)
                    ff += kappa * w * acc[i - j]
                    align += gamma * w * (vv[i - j] - vv[i])
                rear = eta * (vv[i + 1] - vv[i]) if i + 1 < n else 0.0
                acc[i] = base[i] + ff + align + rear
        dg = vv[:-1] - vv[1:]
        return dg, np.clip(acc, -5.0, 3.0)

    sample_every = round(2.0 / dt)
    times, speeds, accels, gap_hist = [], [], [], []
    nsteps = round(duration / dt)
    for step in range(nsteps + 1):
        t = step * dt
        lv, _ = leader(t)
        speed[0] = lv
        if step % sample_every == 0:
            _, a_now = rhs(t, gaps, speed)
            times.append(t); speeds.append(speed.copy()); accels.append(a_now.copy()); gap_hist.append(gaps.copy())
        if step == nsteps:
            break
        k1g, k1v = rhs(t, gaps, speed)
        k2g, k2v = rhs(t + dt / 2, gaps + dt * k1g / 2, speed + dt * k1v / 2)
        k3g, k3v = rhs(t + dt / 2, gaps + dt * k2g / 2, speed + dt * k2v / 2)
        k4g, k4v = rhs(t + dt, gaps + dt * k3g, speed + dt * k3v)
        gaps += dt * (k1g + 2*k2g + 2*k3g + k4g) / 6
        speed += dt * (k1v + 2*k2v + 2*k3v + k4v) / 6
        speed = np.maximum(speed, 0.0)
        speed[0] = leader(t + dt)[0]
        if np.min(gaps) < 0.2:
            raise RuntimeError("open platoon collision")

    time = np.asarray(times); V = np.asarray(speeds); A = np.asarray(accels); G = np.asarray(gap_hist)
    X = reconstruct_positions(time, V, G)
    sigma = np.std(V[:, 1:], axis=1)
    window = (time >= force_start) & (time <= 650)
    amplitudes = 0.5 * (V[window].max(axis=0) - V[window].min(axis=0))
    closing = V[:, 1:] - V[:, :-1]
    ttc = np.full_like(G, np.inf)
    np.divide(G, closing, out=ttc, where=closing > 0.05)
    rec = recovery_time(time, sigma, force_end, 0.05, 60.0)
    metrics = {
        "leader_amplitude": amp,
        "last_vehicle_amplitude": float(amplitudes[-1]),
        "amplification": float(amplitudes[-1] / amp),
        "peak_speed_std": float(np.max(sigma)),
        "recovery_time": rec,
        "min_gap": float(np.min(G)),
        "min_ttc": float(np.min(ttc)),
        "harsh_braking_vehicle_seconds": float(np.sum(A[:, 1:] < -2.0) * 2.0),
        "mean_speed": float(np.mean(V[:, 1:])),
    }
    metrics.update(ev_energy_metrics(time, V[:, 1:], A[:, 1:], ve))
    return {"time": time, "speed": V, "acc": A, "gap": G, "position": X, "sigma": sigma, "amplitudes": amplitudes, "metrics": metrics,
            "events": {"start": force_start, "end": force_end}, "probe_ids": [0, 15, 35, 59]}


def simulate_moving_bottleneck(controlled: bool, *, vsl_limit: float = 100.0 / 3.6,
                               incident_low_speed: float = 18.0, incident_duration: float = 120.0,
                               duration: float = 1200.0, dt: float = 0.1) -> dict:
    """Continuous open traffic with a fixed 1-km VSL zone upstream.

    Vehicles are initialized on an extended upstream reservoir at the spacing
    implied by a short-headway high-demand equilibrium.  They cross x=0 continuously throughout the
    1200-s horizon and leave at x=12 km.  Only vehicles physically on this road
    are included in performance measures and trajectory plots.
    """
    # The road is posted at 120 km/h.  With T=0.7 s and a net gap of 18 m, the
    # deliberately dense pre-incident equilibrium travels at about 75 km/h.
    hp = HIGHWAY_P
    ve = 75.0 / 3.6
    demand_vph = 3600.0 * ve / (18.0 + 5.0)
    road_length, vehicle_length = 12_000.0, 5.0
    gross_spacing = ve / (demand_vph / 3600.0)
    net_gap = gross_spacing - vehicle_length
    incident_start, incident_end = 300.0, 300.0 + incident_duration
    low_speed = incident_low_speed
    # Event-triggered control: the VSL is unavailable before the slow vehicle
    # starts the disturbance.  At t=incident_start the sign activates and the
    # displayed target ramps down over 30 s.  This removes the anticipatory
    # control that was present in the earlier data set.
    control_start, control_end, zone_length = incident_start, 900.0, 1000.0
    ramp_in, ramp_out = 30.0, 90.0
    slow_nominal_position = 6000.0
    zone_start, zone_end = slow_nominal_position - zone_length, slow_nominal_position
    slow_initial_position = slow_nominal_position - ve * incident_start
    slow_id = int(round((road_length - slow_initial_position) / gross_spacing))
    n = int(np.ceil((road_length + ve * duration + 2 * gross_spacing) / gross_spacing))
    position = road_length - np.arange(n, dtype=float) * gross_spacing
    speed = np.full(n, ve, dtype=float)

    def bottleneck_target(t: float) -> float:
        return low_speed if incident_start <= t <= incident_end else ve

    def speed_limit(t: float) -> float:
        if not controlled or t < control_start or t > control_end:
            return hp.v0
        if t < control_start + ramp_in:
            q = (t - control_start) / ramp_in
            return hp.v0 + q * (vsl_limit - hp.v0)
        if t > control_end - ramp_out:
            q = (t - (control_end - ramp_out)) / ramp_out
            return vsl_limit + q * (hp.v0 - vsl_limit)
        return vsl_limit

    def rhs(t: float, x: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        acc = np.zeros_like(v)
        acc[0] = (ve - v[0]) / 1.5
        gap = np.maximum(x[:-1] - x[1:] - vehicle_length, 0.1)
        follower_v = np.maximum(v[1:], 0.0)
        dv = follower_v - v[:-1]
        desired = hp.s0 + np.maximum(0.0, follower_v * hp.T + follower_v * dv / (2.0 * np.sqrt(hp.a * hp.b)))
        local_v0 = np.full(n - 1, hp.v0)
        if controlled and control_start <= t <= control_end:
            in_zone = (x[1:] >= zone_start) & (x[1:] < zone_end)
            # The fixed one-kilometre approach uses a smooth speed profile, so
            # the VSL does not create an artificial compression shock at entry.
            spatial_blend = np.clip((x[1:] - zone_start) / zone_length, 0.0, 1.0)
            local_v0[in_zone] = hp.v0 + spatial_blend[in_zone] * (speed_limit(t) - hp.v0)
        acc[1:] = hp.a * (1.0 - (follower_v / local_v0) ** hp.delta - (desired / gap) ** 2)
        target = bottleneck_target(t)
        acc[slow_id] = np.clip((target - v[slow_id]) / 2.2, -3.5, 2.0)
        return v, np.clip(acc, -5.0, 3.0)

    sample_every = round(2.0 / dt)
    times, positions, speeds, accels, limits = [], [], [], [], []
    nsteps = round(duration / dt)
    for step in range(nsteps + 1):
        t = step * dt
        if step % sample_every == 0:
            _, a_now = rhs(t, position, speed)
            times.append(t); positions.append(position.copy()); speeds.append(speed.copy()); accels.append(a_now.copy()); limits.append(speed_limit(t))
        if step == nsteps:
            break
        k1x, k1v = rhs(t, position, speed)
        k2x, k2v = rhs(t + dt/2, position + dt*k1x/2, speed + dt*k1v/2)
        k3x, k3v = rhs(t + dt/2, position + dt*k2x/2, speed + dt*k2v/2)
        k4x, k4v = rhs(t + dt, position + dt*k3x, speed + dt*k3v)
        position += dt * (k1x + 2*k2x + 2*k3x + k4x) / 6
        speed += dt * (k1v + 2*k2v + 2*k3v + k4v) / 6
        speed = np.maximum(speed, 0.0)
        min_gap_step = float(np.min(position[:-1] - position[1:] - vehicle_length))
        if min_gap_step < 0.2:
            raise RuntimeError(f"continuous moving-bottleneck collision at t={t:.1f}s, gap={min_gap_step:.3f}m")

    time = np.asarray(times); X = np.asarray(positions); V = np.asarray(speeds); A = np.asarray(accels); limits = np.asarray(limits)
    active = (X >= 0.0) & (X <= road_length)
    sigma = np.array([np.std(V[i, active[i]]) if np.any(active[i]) else 0.0 for i in range(len(time))])
    # The incident target is 64.8 km/h; a slightly higher threshold identifies
    # the connected low-speed packet and its upstream-moving tail.
    queue_threshold = 18.2
    queue = np.sum(active & (V < queue_threshold), axis=1)
    queue_tail = np.full(len(time), np.nan)
    for i in range(len(time)):
        queued = active[i] & (V[i] < queue_threshold)
        if np.any(queued):
            queue_tail[i] = np.min(X[i, queued])
    wave_fit = (time >= incident_start + 60.0) & (time <= incident_end) & np.isfinite(queue_tail)
    wave_speed = float(np.polyfit(time[wave_fit], queue_tail[wave_fit], 1)[0]) if np.sum(wave_fit) >= 5 else None
    rec = recovery_time(time, sigma, incident_end, 0.3, 90.0)
    gaps = X[:, :-1] - X[:, 1:] - vehicle_length
    pair_active = active[:, :-1] & active[:, 1:]
    closing = V[:, 1:] - V[:, :-1]
    ttc = np.full_like(gaps, np.inf)
    np.divide(gaps, closing, out=ttc, where=pair_active & (closing > 0.05))
    power = wu_ev_power(V, A)
    total_energy = float(np.trapezoid(np.sum(np.where(active, power, 0.0), axis=1), time) / 3.6e6)
    total_distance = float(np.trapezoid(np.sum(np.where(active, V, 0.0), axis=1), time) / 1000.0)
    arrivals = int(np.sum((X[0] < 0.0) & (X[-1] >= 0.0)))
    departures = int(np.sum((X[0] <= road_length) & (X[-1] > road_length)))
    metrics = {
        "arrival_vehicles": arrivals,
        "departure_vehicles": departures,
        "peak_vehicles_on_road": int(np.max(np.sum(active, axis=1))),
        "peak_queue_vehicles": int(np.max(queue)),
        "queue_speed_threshold": queue_threshold,
        "queue_tail_wave_speed_mps": wave_speed,
        "peak_speed_std": float(np.max(sigma)),
        "recovery_time": rec,
        "min_speed": float(np.min(V[active])),
        "min_gap": float(np.min(gaps[pair_active])),
        "min_ttc": float(np.min(ttc)),
        "delay_vehicle_hours": float(np.trapezoid(np.sum(np.where(active, np.maximum(ve - V, 0.0), 0.0), axis=1) / ve, time) / 3600.0),
        "harsh_braking_vehicle_seconds": float(np.sum(active & (A < -2.0)) * 2.0),
        "total_ev_energy_kwh": total_energy,
        "ev_energy_intensity_kwh_per_100km": float(100.0 * total_energy / total_distance),
        "mean_speed": float(np.sum(V[active]) / np.sum(active)),
    }
    active_with_slow = (X >= 0.0) & (X <= road_length)
    return {"time": time, "speed": V, "acc": A, "gap": gaps, "position": X, "active": active_with_slow, "sigma": sigma,
            "queue": queue, "queue_tail": queue_tail, "limit": limits, "metrics": metrics,
            "events": {"incident_start": incident_start, "incident_end": incident_end, "control_start": control_start, "control_end": control_end,
                       "zone_start": zone_start, "zone_end": zone_end, "zone_length": zone_length,
                       "road_length": road_length, "slow_id": slow_id},
            "probe_ids": [slow_id, slow_id + 20, slow_id + 60, slow_id + 120]}


def q8(values: np.ndarray, lo: float, hi: float) -> str:
    arr = np.clip(np.rint((values - lo) / (hi - lo) * 255), 0, 255).astype(np.uint8)
    return base64.b64encode(arr.tobytes()).decode("ascii")


def q16(values: np.ndarray, lo: float, hi: float) -> str:
    arr = np.clip(np.rint((values - lo) / (hi - lo) * 65535), 0, 65535).astype("<u2")
    return base64.b64encode(arr.tobytes()).decode("ascii")


def compact_case(result: dict, *, vmin: float, vmax: float) -> dict:
    # 4 s trajectory samples; 2 s probe samples.
    heat_idx = np.arange(0, len(result["time"]), 2)
    probes = result["probe_ids"]
    positions = result["position"][heat_idx]
    active = result.get("active", np.ones_like(result["position"], dtype=bool))[heat_idx]
    position_min = float(np.min(positions[active]))
    position_max = float(np.max(positions[active]))
    compact = {
        "time": np.round(result["time"], 2).tolist(),
        "heatTime": np.round(result["time"][heat_idx], 2).tolist(),
        "n": int(result["speed"].shape[1]),
        "vmin": vmin,
        "vmax": vmax,
        "speedQ8": q8(result["speed"][heat_idx], vmin, vmax),
        "activeQ8": q8(active.astype(float), 0.0, 1.0),
        "positionMin": position_min,
        "positionMax": position_max,
        "positionQ16": q16(positions, position_min, position_max),
        "probeIds": probes,
        "probes": np.round(result["speed"][:, probes], 4).T.tolist(),
        "sigma": np.round(result["sigma"], 5).tolist(),
        "metrics": result["metrics"],
        "events": result["events"],
    }
    if "queue_tail" in result:
        compact["queueTail"] = [None if not np.isfinite(x) else round(float(x), 2) for x in result["queue_tail"]]
    return compact


def pct(before: float, after: float) -> float:
    return 100.0 * (before - after) / before


def main() -> None:
    open_hdv = simulate_open_platoon(False)
    open_av = simulate_open_platoon(True)
    bottleneck_none = simulate_moving_bottleneck(False)

    # Online design uses stability only. Safety, queueing, delay and energy
    # are deliberately excluded because they are ex-post simulation outputs.
    # The statutory motorway limit is 120 km/h.  Candidate VSL commands span
    # 114--120 km/h, and the online objective remains stability-only.  The
    # lower bound is an exogenous operational constraint on a one-kilometre
    # control zone, not a simulation-derived safety objective.
    free_speed = 75.0 / 3.6
    demand_vph = 3600.0 * free_speed / (18.0 + 5.0)
    upstream_gap = free_speed / (demand_vph / 3600.0) - 5.0
    candidate_kph = np.arange(114.0, 120.1, 1.0)
    candidates = candidate_kph / 3.6
    scan = []
    for displayed_limit in candidates:
        # A posted VSL is used directly as the IDM desired-speed parameter.
        internal_v0 = float(displayed_limit)
        eq_speed, phi = idm_long_wave_margin(upstream_gap, replace(HIGHWAY_P, v0=internal_v0))
        scan.append({
            "limit": float(displayed_limit), "limit_kph": float(displayed_limit * 3.6), "idm_v0": float(internal_v0),
            "equilibrium_speed": float(eq_speed), "phi": float(phi), "feasible": bool(phi > 0.0),
            "objective": float(phi),
        })
    selected = max(scan, key=lambda row: row["phi"])
    selected_limit = selected["limit"]
    bottleneck_vsl = simulate_moving_bottleneck(True, vsl_limit=selected["idm_v0"])
    baseline_eq, baseline_phi = idm_long_wave_margin(upstream_gap, HIGHWAY_P)

    metrics = {
        "case1": {
        "setup": {"N": 60, "ve": 8.0, "leader_amplitude": 0.6, "period": 48.0, "cycles": 5},
            "HDV": open_hdv["metrics"], "AV": open_av["metrics"],
            "benefit": {
                "amplitude_reduction_percent": pct(open_hdv["metrics"]["last_vehicle_amplitude"], open_av["metrics"]["last_vehicle_amplitude"]),
                "peak_std_reduction_percent": pct(open_hdv["metrics"]["peak_speed_std"], open_av["metrics"]["peak_speed_std"]),
                "harsh_braking_reduction_percent": pct(open_hdv["metrics"]["harsh_braking_vehicle_seconds"], open_av["metrics"]["harsh_braking_vehicle_seconds"]) if open_hdv["metrics"]["harsh_braking_vehicle_seconds"] else None,
                "excess_ev_energy_reduction_percent": pct(open_hdv["metrics"]["mean_excess_ev_energy_wh"], open_av["metrics"]["mean_excess_ev_energy_wh"]),
                "ev_energy_intensity_reduction_percent": pct(open_hdv["metrics"]["ev_energy_intensity_kwh_per_100km"], open_av["metrics"]["ev_energy_intensity_kwh_per_100km"]),
            },
        },
        "case2": {
            "setup": {"duration": 1200.0, "demand_vph": demand_vph, "road_length_m": 12000.0,
                      "statutory_limit_kph": 120.0, "T": HIGHWAY_P.T, "ve": free_speed,
                      "slow_vehicle": bottleneck_none["events"]["slow_id"], "low_speed": 18.0, "incident": [300.0, 420.0],
                      "vsl": [bottleneck_vsl["events"]["control_start"], bottleneck_vsl["events"]["control_end"]],
                      "zone_start_m": 5000.0, "zone_end_m": 6000.0,
                      "zone_length_m": 1000.0, "limit": selected_limit, "limit_kph": selected["limit_kph"],
                      "idm_v0": selected["idm_v0"]},
            "vsl_design": {
                "baseline": {"equilibrium_speed": baseline_eq, "phi": baseline_phi, "idm_v0": HIGHWAY_P.v0},
                "admissible_range": [float(candidate_kph[0]), float(candidate_kph[-1])],
                "objective": "maximize Phi(u) over the admissible displayed VSL range; no safety or simulation metric in online optimization",
                "selected": selected,
                "scan": scan,
            },
            "uncontrolled": bottleneck_none["metrics"], "VSL": bottleneck_vsl["metrics"],
            "benefit": {
                "queue_reduction_percent": pct(bottleneck_none["metrics"]["peak_queue_vehicles"], bottleneck_vsl["metrics"]["peak_queue_vehicles"]),
                "peak_std_reduction_percent": pct(bottleneck_none["metrics"]["peak_speed_std"], bottleneck_vsl["metrics"]["peak_speed_std"]),
                "delay_reduction_percent": pct(bottleneck_none["metrics"]["delay_vehicle_hours"], bottleneck_vsl["metrics"]["delay_vehicle_hours"]),
                "harsh_braking_reduction_percent": pct(bottleneck_none["metrics"]["harsh_braking_vehicle_seconds"], bottleneck_vsl["metrics"]["harsh_braking_vehicle_seconds"]) if bottleneck_none["metrics"]["harsh_braking_vehicle_seconds"] else None,
                "ev_energy_intensity_reduction_percent": pct(bottleneck_none["metrics"]["ev_energy_intensity_kwh_per_100km"], bottleneck_vsl["metrics"]["ev_energy_intensity_kwh_per_100km"]),
            },
        },
    }
    METRICS.write_text(json.dumps(metrics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    compact = {
        "open": {"baseline": compact_case(open_hdv, vmin=0, vmax=14), "control": compact_case(open_av, vmin=0, vmax=14)},
        "bottleneck": {"baseline": compact_case(bottleneck_none, vmin=0, vmax=35), "control": compact_case(bottleneck_vsl, vmin=0, vmax=35)},
        "metrics": metrics,
    }
    COMPACT.write_text(json.dumps(compact, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    LAB_DATA.write_text(
        "window.APPLICATION_CASES_DATA="
        + json.dumps(compact, ensure_ascii=False, separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )

    cjk = "/System/Library/Fonts/Hiragino Sans GB.ttc"
    font_manager.fontManager.addfont(cjk)
    family = font_manager.FontProperties(fname=cjk).get_name()
    plt.rcParams.update({"font.family": family, "font.sans-serif": [family, "DejaVu Sans"], "axes.unicode_minus": False, "font.size": 8.2})
    blue, red, teal, amber = "#315b8a", "#b65043", "#1f7564", "#b3832d"

    def trajectory(a, r, title, vmin, vmax):
        t = r["time"]
        x = r["position"] / 1000.0
        active = r.get("active", np.ones_like(r["position"], dtype=bool))
        segments, colors = [], []
        for i in range(x.shape[1]):
            pts = np.column_stack([t, x[:, i]])
            keep = active[:-1, i] & active[1:, i]
            if np.any(keep):
                segments.extend(np.stack([pts[:-1][keep], pts[1:][keep]], axis=1))
                colors.extend(0.5 * (r["speed"][:-1, i][keep] + r["speed"][1:, i][keep]))
        lc = LineCollection(segments, cmap="turbo", norm=Normalize(vmin=vmin, vmax=vmax), linewidths=0.55, alpha=0.90)
        lc.set_array(np.asarray(colors))
        a.add_collection(lc); a.autoscale(); a.margins(x=0.01, y=0.02)
        a.set_title(title); a.set_xlabel("time [s]"); a.set_ylabel("road position x [km]")
        a.grid(alpha=.12)
        return lc

    ASSET_OPEN.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(12.4, 3.55), constrained_layout=True)
    lc1 = trajectory(ax[0], open_hdv, "(a) HDV/IDM 位置--时间轨迹", 0, 14)
    trajectory(ax[1], open_av, "(b) 协同 AV 位置--时间轨迹", 0, 14)
    for i, pid in enumerate(open_hdv["probe_ids"]):
        ax[2].plot(open_hdv["time"], open_hdv["speed"][:,pid], color=red, lw=1.0, alpha=.45, ls="--", label=f"HDV n={pid}" if i in (0,3) else None)
        ax[2].plot(open_av["time"], open_av["speed"][:,pid], color=blue, lw=1.0, alpha=.75, label=f"AV n={pid}" if i in (0,3) else None)
    ax[2].axvspan(100,340,color=amber,alpha=.12); ax[2].set_title("(c) 指定车辆速度响应"); ax[2].set_xlabel("time [s]"); ax[2].set_ylabel("speed [m/s]"); ax[2].legend(frameon=False,fontsize=7,ncol=2); ax[2].grid(alpha=.18)
    fig.colorbar(lc1, ax=ax[:2], shrink=.76, pad=.015, label="trajectory speed [m/s]")
    fig.savefig(ASSET_OPEN, dpi=220, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(2, 3, figsize=(12.4, 7.0), constrained_layout=True)
    selected_kph = selected["limit_kph"]
    lc2 = trajectory(ax[0,0], bottleneck_none, "(a) 120 km/h 法定上限，无控制", 0, 35)
    trajectory(ax[0,1], bottleneck_vsl, f"(b) 稳定性约束 VSL={selected_kph:.0f} km/h", 0, 35)
    for axis, result in zip(ax[0,:2], [bottleneck_none, bottleneck_vsl]):
        good_tail = np.isfinite(result["queue_tail"])
        axis.plot(result["time"][good_tail], result["queue_tail"][good_tail] / 1000.0,
                  color="black", ls="--", lw=1.0,
                  label=fr"queue tail; formation-stage $c_q={result['metrics']['queue_tail_wave_speed_mps']:.2f}$ m/s")
        axis.legend(frameon=False, fontsize=7, loc="lower left")
    for i, pid in enumerate(bottleneck_none["probe_ids"]):
        ax[0,2].plot(bottleneck_none["time"], bottleneck_none["speed"][:,pid], color=red, lw=1.0, alpha=.45, ls="--", label=f"none n={pid}" if i in (0,3) else None)
        ax[0,2].plot(bottleneck_vsl["time"], bottleneck_vsl["speed"][:,pid], color=teal, lw=1.0, alpha=.75, label=f"VSL n={pid}" if i in (0,3) else None)
    ax[0,2].axvspan(300,420,color=red,alpha=.09); ax[0,2].axvspan(300,900,color=teal,alpha=.07)
    ax[0,2].set_title("(c) 指定车辆速度响应"); ax[0,2].set_xlabel("time [s]"); ax[0,2].set_ylabel("speed [m/s]"); ax[0,2].legend(frameon=False,fontsize=7,ncol=2); ax[0,2].grid(alpha=.18)
    limits = np.array([row["limit_kph"] for row in scan]); phis = np.array([row["phi"] for row in scan])
    ax[1,0].plot(limits, phis, color=blue, marker="o", ms=3); ax[1,0].axhline(0.0,color=red,ls="--",lw=1,label=r"neutral boundary $\Phi=0$")
    ax[1,0].axvline(selected_kph,color=teal,ls=":",lw=1.4,label=f"max-margin VSL {selected_kph:.0f} km/h")
    ax[1,0].set_title("(d) 120 km/h 上限内最大化长波稳定裕度"); ax[1,0].set_xlabel("displayed VSL target speed [km/h]"); ax[1,0].set_ylabel(r"long-wave margin $\Phi$ [s$^{-2}$]"); ax[1,0].grid(alpha=.18); ax[1,0].legend(frameon=False,fontsize=7)
    bar_values = [baseline_phi, selected["phi"]]
    bars = ax[1,1].bar(["120 km/h limit", f"VSL {selected_kph:.0f} km/h"], bar_values, color=[red, teal], alpha=.8)
    ax[1,1].axhline(0.0,color="black",lw=.8)
    for bar, val in zip(bars, bar_values):
        ax[1,1].text(bar.get_x()+bar.get_width()/2, val + (0.00045 if val >= 0 else -0.00065), f"{val:+.4f}", ha="center", va="bottom" if val >= 0 else "top", fontsize=8)
    ax[1,1].set_title("(e) 小 T 失稳工作点被 VSL 推回稳定侧")
    ax[1,1].set_ylabel(r"long-wave margin $\Phi$ [s$^{-2}$]"); ax[1,1].grid(axis="y",alpha=.18)
    slow_id = bottleneck_vsl["events"]["slow_id"]
    slow_x = bottleneck_vsl["position"][:, slow_id] / 1000.0
    zone_on = ((bottleneck_vsl["time"] >= bottleneck_vsl["events"]["control_start"])
               & (bottleneck_vsl["time"] <= bottleneck_vsl["events"]["control_end"]))
    ax[1,2].plot(bottleneck_vsl["time"], slow_x,color=red,lw=1.5,label="slow vehicle")
    ax[1,2].fill_between(bottleneck_vsl["time"], 5.0, 6.0, where=zone_on,color=teal,alpha=.25,label="fixed VSL zone: x=5--6 km")
    ax[1,2].set_ylim(0,12); ax[1,2].set_title("(f) VSL 固定在事件点上游 1 km")
    ax[1,2].set_xlabel("time [s]"); ax[1,2].set_ylabel("road position x [km]"); ax[1,2].grid(alpha=.18); ax[1,2].legend(frameon=False,fontsize=7)
    fig.colorbar(lc2, ax=ax[0,:2], shrink=.72, pad=.015, label="trajectory speed [m/s]")
    fig.savefig(ASSET_VSL, dpi=220, facecolor="white")
    plt.close(fig)
    print(json.dumps(metrics, ensure_ascii=False, indent=2))
    print(COMPACT)
    print(ASSET_OPEN)
    print(ASSET_VSL)


if __name__ == "__main__":
    main()
