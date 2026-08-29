from __future__ import annotations

import json
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class ReferenceOutputTests(unittest.TestCase):
    def test_manifest_has_all_published_examples(self):
        manifest = json.loads((ROOT / "experiments" / "manifest.json").read_text(encoding="utf-8"))
        identifiers = {item["id"] for item in manifest["experiments"]}
        self.assertEqual(
            identifiers,
            {
                "M01", "M02", "M03", "M04", "M05", "M06",
                "W01", "W02", "W03", "D01", "D02",
            },
        )

    def test_ngsim_derived_table_matches_published_metrics(self):
        episodes = pd.read_csv(ROOT / "data" / "derived" / "us101_pair_episode_metrics.csv")
        metrics = json.loads(
            (ROOT / "results" / "ngsim_us101_d01_metrics.json").read_text(encoding="utf-8")
        )
        self.assertEqual(len(episodes), metrics["valid_pair_episodes"])
        self.assertAlmostEqual(float(episodes["gain"].median()), metrics["median_pair_gain"], places=12)
        self.assertAlmostEqual(
            float(np.mean(episodes["gain"] > 1.0)),
            metrics["fraction_pair_gain_above_one"],
            places=12,
        )

    def test_worked_updates_preserve_published_steps(self):
        metrics = json.loads(
            (ROOT / "results" / "worked_step_examples.json").read_text(encoding="utf-8")
        )
        micro = metrics["W01_three_vehicle_idm_rk4"]
        self.assertAlmostEqual(micro["states"][1]["speed_mps"][0], 10.019862100874903)
        lwr = metrics["W02_four_cell_lwr"]
        self.assertAlmostEqual(lwr["states"][1]["density_veh_per_km"][1], 18.0)
        arz = metrics["W03_three_cell_arz"]
        for state in arz["states"]:
            self.assertAlmostEqual(state["vehicles_in_periodic_domain"], 4.5, places=12)

    def test_ngsim_d02_keeps_calibration_validation_and_optimization_separate(self):
        metrics = json.loads(
            (ROOT / "results" / "ngsim_stability_application_d02.json").read_text(
                encoding="utf-8"
            )
        )
        single = metrics["single_vehicle_calibration"]
        fleet = metrics["stability_aware_platoon_calibration"]
        design = metrics["stability_optimization"]
        self.assertLess(
            fleet["training"]["speed_rmse_mps"],
            single["training"]["speed_rmse_mps"],
        )
        self.assertLess(
            fleet["validation"]["speed_rmse_mps"],
            single["validation"]["speed_rmse_mps"],
        )
        self.assertAlmostEqual(design["selected_gain"], 0.04, places=12)
        self.assertLess(
            design["ring_sigma_ratio_controlled"],
            design["ring_sigma_ratio_baseline"],
        )


if __name__ == "__main__":
    unittest.main()
