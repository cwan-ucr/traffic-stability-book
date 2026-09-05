#!/usr/bin/env python3
"""D01: a reproducible stability-oriented analysis of an official NGSIM slice.

The example uses three minutes from US-101 lane 2.  It reports descriptive
statistics, extracts one eight-vehicle platoon, and estimates natural
leader-to-follower disturbance gains from continuous episodes of at least
20 seconds.  These empirical gains are evidence, not a proof of global string
stability: lane changes, unobserved reactions, measurement noise, and
non-stationarity remain in the data.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlretrieve

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
from scipy.signal import detrend, savgol_filter


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "raw" / "us101_lane2_180s.csv"
DERIVED = ROOT / "data" / "derived"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

DATASET_ENDPOINT = "https://data.transportation.gov/resource/8ect-6jqj.csv"
DATASET_DOI = "https://doi.org/10.21949/1504477"
SELECT = ",".join(
    [
        "vehicle_id",
        "frame_id",
        "total_frames",
        "global_time",
        "local_y",
        "v_length",
        "v_class",
        "v_vel",
        "v_acc",
        "lane_id",
        "preceding",
        "following",
        "space_headway",
        "time_headway",
        "location",
    ]
)
WHERE = (
    'location="us-101" and lane_id="2" and '
    "global_time between 1118846979700 and 1118847159700"
)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Hiragino Sans GB", "Songti SC", "Arial Unicode MS", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 9,
    }
)


def download_official_slice(destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    query = urlencode(
        {
            "$select": SELECT,
            "$where": WHERE,
            "$order": "global_time,vehicle_id",
            "$limit": "50000",
        }
    )
    urlretrieve(f"{DATASET_ENDPOINT}?{query}", destination)


def load_and_convert(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    expected = set(SELECT.split(","))
    missing = expected.difference(data.columns)
    if missing:
        raise ValueError(f"missing NGSIM fields: {sorted(missing)}")
    data = data.sort_values(["global_time", "vehicle_id"]).reset_index(drop=True)
    data["time_s"] = (data["global_time"] - data["global_time"].min()) / 1000.0
    data["position_m"] = data["local_y"] * 0.3048
    data["speed_mps"] = data["v_vel"] * 0.3048
    data["acceleration_mps2"] = data["v_acc"] * 0.3048
    data["space_headway_m"] = data["space_headway"] * 0.3048
    return data


def longest_chain_at_time(data: pd.DataFrame, target_time_s: float) -> list[int]:
    available = np.sort(data["time_s"].unique())
    instant = float(available[np.argmin(abs(available - target_time_s))])
    snapshot = data[data["time_s"] == instant].set_index("vehicle_id")
    best: list[int] = []
    for follower in snapshot.index:
        current = int(follower)
        chain, seen = [], set()
        while current in snapshot.index and current not in seen:
            chain.append(current)
            seen.add(current)
            leader = int(snapshot.loc[current, "preceding"])
            if leader <= 0 or leader not in snapshot.index:
                break
            current = leader
        if len(chain) > len(best):
            best = chain
    return best


def select_platoon(data: pd.DataFrame, vehicle_count: int = 8) -> pd.DataFrame:
    """Select a continuous chain and return it in leader-to-follower order."""

    upstream_to_downstream = longest_chain_at_time(data, target_time_s=100.0)
    if len(upstream_to_downstream) < vehicle_count:
        raise RuntimeError("the selected snapshot contains no sufficiently long chain")
    chosen = upstream_to_downstream[:vehicle_count]
    intervals = []
    for follower, leader in zip(chosen[:-1], chosen[1:]):
        pair = data[(data["vehicle_id"] == follower) & (data["preceding"] == leader)]
        intervals.append((pair["global_time"].min(), pair["global_time"].max()))
    start = max(item[0] for item in intervals)
    end = min(item[1] for item in intervals)
    if (end - start) / 1000.0 < 20.0 - 1.0e-9:
        raise RuntimeError("the automatically selected chain is shorter than 20 s")

    leader_to_follower = list(reversed(chosen))
    selected = data[
        data["vehicle_id"].isin(leader_to_follower)
        & data["global_time"].between(start, end)
    ].copy()
    order = {vehicle: index for index, vehicle in enumerate(leader_to_follower)}
    selected["vehicle_order"] = selected["vehicle_id"].map(order)
    selected["episode_time_s"] = (selected["global_time"] - start) / 1000.0
    return selected.sort_values(["vehicle_order", "global_time"])


def _continuous_runs(times_ms: np.ndarray) -> list[tuple[int, int]]:
    breaks = np.flatnonzero(np.diff(times_ms) > 100)
    starts = np.r_[0, breaks + 1]
    ends = np.r_[breaks + 1, len(times_ms)]
    return [(int(start), int(end)) for start, end in zip(starts, ends)]


def pair_episode_metrics(data: pd.DataFrame) -> pd.DataFrame:
    """Measure gain and response lag for continuous natural following episodes."""

    records = []
    leader_rows = data[["vehicle_id", "global_time", "speed_mps"]].rename(
        columns={"vehicle_id": "preceding", "speed_mps": "leader_speed_mps"}
    )
    for (follower, leader), group in data[data["preceding"] > 0].groupby(
        ["vehicle_id", "preceding"], sort=False
    ):
        group = group.sort_values("global_time")
        times = group["global_time"].to_numpy()
        for start, end in _continuous_runs(times):
            run = group.iloc[start:end]
            if len(run) < 201:
                continue
            merged = run.merge(
                leader_rows[leader_rows["preceding"] == leader],
                on=["preceding", "global_time"],
                how="inner",
            ).sort_values("global_time")
            if len(merged) < 201:
                continue
            follower_speed = detrend(savgol_filter(merged["speed_mps"], 21, 3))
            leader_speed = detrend(savgol_filter(merged["leader_speed_mps"], 21, 3))
            leader_std = float(np.std(leader_speed))
            follower_std = float(np.std(follower_speed))
            if leader_std < 0.25:
                continue

            correlations = []
            maximum_lag = min(50, len(merged) // 3)
            for lag in range(maximum_lag + 1):
                leader_part = leader_speed[: len(leader_speed) - lag or None]
                follower_part = follower_speed[lag:]
                correlations.append(float(np.corrcoef(leader_part, follower_part)[0, 1]))
            best_lag = int(np.nanargmax(correlations))
            records.append(
                {
                    "follower_id": int(follower),
                    "leader_id": int(leader),
                    "samples": int(len(merged)),
                    "duration_s": float((merged["global_time"].iloc[-1] - merged["global_time"].iloc[0]) / 1000.0),
                    "leader_std_mps": leader_std,
                    "follower_std_mps": follower_std,
                    "gain": follower_std / leader_std,
                    "lag_s": 0.1 * best_lag,
                    "maximum_correlation": correlations[best_lag],
                }
            )
    return pd.DataFrame.from_records(records)


def make_figure(data: pd.DataFrame, platoon: pd.DataFrame, episodes: pd.DataFrame, path: Path, *, english: bool = False) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.6), constrained_layout=True)

    axis = axes[0, 0]
    collection = None
    for _, trajectory in data.groupby("vehicle_id"):
        trajectory = trajectory.sort_values("time_s")
        if len(trajectory) < 2:
            continue
        points = np.column_stack([trajectory["time_s"], trajectory["position_m"]]).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        collection = LineCollection(segments, cmap="turbo", norm=plt.Normalize(4.0, 24.0), linewidth=0.75)
        collection.set_array(trajectory["speed_mps"].to_numpy()[:-1])
        axis.add_collection(collection)
    axis.autoscale()
    if collection is not None:
        colorbar = fig.colorbar(collection, ax=axis, pad=0.02)
        colorbar.set_label("speed [m/s]" if english else "速度 [m/s]")
    axis.set(title="D01(a) US-101 lane 2: position--time trajectories" if english else "D01(a) US-101 车道 2：位置-时间轨迹", xlabel="time [s]" if english else "时间 [s]", ylabel="road position [m]" if english else "道路纵向位置 [m]")

    axis = axes[0, 1]
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, platoon["vehicle_order"].nunique()))
    amplitudes = []
    vehicle_ids = []
    for order, group in platoon.groupby("vehicle_order"):
        group = group.sort_values("episode_time_s")
        smooth = savgol_filter(group["speed_mps"].to_numpy(), 21, 3)
        axis.plot(group["episode_time_s"], smooth, color=colors[int(order)], label=f"{int(order)}: {int(group['vehicle_id'].iloc[0])}")
        amplitudes.append(float(np.std(detrend(smooth))))
        vehicle_ids.append(int(group["vehicle_id"].iloc[0]))
    axis.set(title="D01(b) continuous eight-vehicle speed response" if english else "D01(b) 连续八车队速度响应", xlabel="episode time [s]" if english else "片段时间 [s]", ylabel="speed [m/s]" if english else "速度 [m/s]")
    axis.legend(title="order: vehicle ID" if english else "次序: 车辆 ID", ncol=2, fontsize=7, frameon=False)

    axis = axes[1, 0]
    axis.bar(np.arange(len(amplitudes)), amplitudes, color=colors)
    axis.set_xticks(np.arange(len(amplitudes)), [str(item) for item in vehicle_ids], rotation=35)
    axis.set(title="D01(c) detrended fluctuation amplitude within one platoon" if english else "D01(c) 同一车队的去趋势波动幅值", xlabel="leader to upstream follower (vehicle ID)" if english else "从领车到上游跟车（车辆 ID）", ylabel="speed standard deviation [m/s]" if english else "速度标准差 [m/s]")

    axis = axes[1, 1]
    stable = episodes["gain"] <= 1.0
    axis.scatter(episodes.loc[stable, "lag_s"], episodes.loc[stable, "gain"], s=25, alpha=0.75, color="#2f7d67", label="gain <= 1")
    axis.scatter(episodes.loc[~stable, "lag_s"], episodes.loc[~stable, "gain"], s=25, alpha=0.75, color="#c04a3a", label="gain > 1")
    axis.axhline(1.0, color="0.2", lw=1.0, ls="--")
    axis.set(title="D01(d) empirical gain in natural following episodes" if english else "D01(d) 自然跟驰片段的逐车扰动增益", xlabel="maximum-correlation response lag [s]" if english else "最大相关响应滞后 [s]", ylabel="follower/leader speed standard deviation" if english else "跟车/前车速度标准差")
    axis.legend(frameon=False)

    for label, axis in zip(["(a)", "(b)", "(c)", "(d)"], axes.flat):
        axis.text(0.01, 0.98, label, transform=axis.transAxes, ha="left", va="top", fontweight="bold")
        axis.grid(alpha=0.16)
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def run(input_path: Path, *, english_figure: bool = False) -> None:
    data = load_and_convert(input_path)
    platoon = select_platoon(data)
    episodes = pair_episode_metrics(data)
    if episodes.empty:
        raise RuntimeError("no valid leader-follower episode was extracted")

    DERIVED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    platoon_columns = [
        "episode_time_s",
        "vehicle_order",
        "vehicle_id",
        "position_m",
        "speed_mps",
        "space_headway_m",
        "time_headway",
        "preceding",
    ]
    platoon[platoon_columns].to_csv(DERIVED / "us101_selected_platoon.csv", index=False)
    episodes.to_csv(DERIVED / "us101_pair_episode_metrics.csv", index=False)

    valid_headway = data[(data["preceding"] > 0) & data["time_headway"].between(0.1, 10.0)]["time_headway"]
    random = np.random.default_rng(20260828)
    gains = episodes["gain"].to_numpy()
    bootstrap = gains[random.integers(0, len(gains), size=(20000, len(gains)))]
    median_gain_ci = np.quantile(np.median(bootstrap, axis=1), [0.025, 0.975])
    amplified_fraction_ci = np.quantile(np.mean(bootstrap > 1.0, axis=1), [0.025, 0.975])
    metrics = {
        "source": DATASET_DOI,
        "location": "US-101 southbound, lane 2",
        "window_s": float(data["time_s"].max() - data["time_s"].min()),
        "rows": int(len(data)),
        "vehicles": int(data["vehicle_id"].nunique()),
        "median_speed_mps": float(data["speed_mps"].median()),
        "median_time_headway_s": float(valid_headway.median()),
        "time_headway_iqr_s": [float(valid_headway.quantile(0.25)), float(valid_headway.quantile(0.75))],
        "selected_platoon_vehicle_ids_leader_to_follower": [
            int(item) for item in platoon.sort_values("vehicle_order")["vehicle_id"].drop_duplicates()
        ],
        "selected_platoon_duration_s": float(platoon["episode_time_s"].max()),
        "valid_pair_episodes": int(len(episodes)),
        "median_pair_gain": float(episodes["gain"].median()),
        "pair_gain_iqr": [float(episodes["gain"].quantile(0.25)), float(episodes["gain"].quantile(0.75))],
        "median_pair_gain_bootstrap_95pct_ci": [float(value) for value in median_gain_ci],
        "fraction_pair_gain_above_one": float(np.mean(episodes["gain"] > 1.0)),
        "fraction_pair_gain_above_one_bootstrap_95pct_ci": [float(value) for value in amplified_fraction_ci],
        "median_response_lag_s": float(episodes["lag_s"].median()),
        "median_maximum_correlation": float(episodes["maximum_correlation"].median()),
        "method_note": "21-frame Savitzky-Golay smoothing, linear detrending, episodes >=20 s, leader std >=0.25 m/s",
        "interpretation_limit": "Natural episodes are observational and do not by themselves prove model-level string stability or causality.",
    }
    metrics_path = RESULTS / "ngsim_us101_d01_metrics.json"
    metrics_path.write_text(json.dumps(metrics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    figure_path = FIGURES / ("ngsim_us101_d01_en.png" if english_figure else "ngsim_us101_d01.png")
    make_figure(data, platoon, episodes, figure_path, english=english_figure)
    print(json.dumps(metrics, ensure_ascii=False, indent=2))
    print(f"figure: {figure_path}")
    print(f"metrics: {metrics_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--download", action="store_true", help="download the official three-minute slice first")
    parser.add_argument("--english-figure", action="store_true", help="write an English-labelled figure")
    args = parser.parse_args()
    if args.download or not args.input.exists():
        download_official_slice(args.input)
    run(args.input, english_figure=args.english_figure)


if __name__ == "__main__":
    main()
