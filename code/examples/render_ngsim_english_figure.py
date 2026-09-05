#!/usr/bin/env python3
"""Render the D01 derived NGSIM evidence using English-only figure labels."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
from scipy.signal import detrend, savgol_filter


ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "data" / "derived"
OUT = ROOT / "figures" / "ngsim_us101_d01_en.png"


def main() -> None:
    platoon = pd.read_csv(DERIVED / "us101_selected_platoon.csv")
    episodes = pd.read_csv(DERIVED / "us101_pair_episode_metrics.csv")
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.6), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, platoon.vehicle_order.nunique()))

    ax = axes[0, 0]
    collection = None
    for order, group in platoon.groupby("vehicle_order"):
        group = group.sort_values("episode_time_s")
        points = np.column_stack([group.episode_time_s, group.position_m]).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        collection = LineCollection(segments, cmap="turbo", norm=plt.Normalize(4, 16), linewidth=1.5)
        collection.set_array(group.speed_mps.to_numpy()[:-1])
        ax.add_collection(collection)
    ax.autoscale()
    fig.colorbar(collection, ax=ax, pad=.02, label="speed [m/s]")
    ax.set(title="D01(a) selected eight-vehicle position--time trajectories", xlabel="episode time [s]", ylabel="road position [m]")

    ax = axes[0, 1]
    amplitudes, ids = [], []
    for order, group in platoon.groupby("vehicle_order"):
        group = group.sort_values("episode_time_s")
        smooth = savgol_filter(group.speed_mps.to_numpy(), 21, 3)
        ax.plot(group.episode_time_s, smooth, color=colors[int(order)], label=f"{int(order)}: {int(group.vehicle_id.iloc[0])}")
        amplitudes.append(float(np.std(detrend(smooth))))
        ids.append(int(group.vehicle_id.iloc[0]))
    ax.set(title="D01(b) continuous eight-vehicle speed response", xlabel="episode time [s]", ylabel="speed [m/s]")
    ax.legend(title="order: vehicle ID", ncol=2, fontsize=7, frameon=False)

    ax = axes[1, 0]
    ax.bar(np.arange(len(amplitudes)), amplitudes, color=colors)
    ax.set_xticks(np.arange(len(ids)), [str(i) for i in ids], rotation=35)
    ax.set(title="D01(c) detrended fluctuation amplitude within the platoon", xlabel="leader to upstream follower (vehicle ID)", ylabel="speed standard deviation [m/s]")

    ax = axes[1, 1]
    stable = episodes.gain <= 1
    ax.scatter(episodes.loc[stable, "lag_s"], episodes.loc[stable, "gain"], s=25, alpha=.75, color="#2f7d67", label="gain ≤ 1")
    ax.scatter(episodes.loc[~stable, "lag_s"], episodes.loc[~stable, "gain"], s=25, alpha=.75, color="#c04a3a", label="gain > 1")
    ax.axhline(1, color="0.2", lw=1, ls="--")
    ax.set(title="D01(d) empirical gain in natural following episodes", xlabel="maximum-correlation response lag [s]", ylabel="follower/leader speed standard deviation")
    ax.legend(frameon=False)
    for panel, ax in zip("abcd", axes.flat):
        ax.text(.01, .98, f"({panel})", transform=ax.transAxes, ha="left", va="top", fontweight="bold")
        ax.grid(alpha=.16)
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=240, bbox_inches="tight")


if __name__ == "__main__":
    main()
