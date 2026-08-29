from __future__ import annotations

import unittest

import numpy as np

from traffic_stability.idm import IDMParameters, equilibrium_gap, idm_acceleration
from traffic_stability.linear import dominant_ring_mode, idm_derivatives
from traffic_stability.macroscopic import solve_lwr_riemann


class TrafficStabilityCoreTests(unittest.TestCase):
    def test_idm_equilibrium_residual(self):
        for speed in (3.0, 8.0, 25.0):
            gap = equilibrium_gap(speed)
            acceleration = idm_acceleration(
                np.array([speed]), np.array([gap]), np.array([speed])
            )[0]
            self.assertAlmostEqual(acceleration, 0.0, places=11)

    def test_long_wave_margin_changes_sign(self):
        self.assertLess(idm_derivatives(8.0)["margin"], 0.0)
        self.assertGreater(idm_derivatives(25.0)["margin"], 0.0)

    def test_finite_ring_can_miss_the_dangerous_long_wave(self):
        self.assertLess(dominant_ring_mode(8.0, 12)["growth_rate_per_s"], 0.0)
        self.assertGreater(dominant_ring_mode(8.0, 20)["growth_rate_per_s"], 0.0)

    def test_front_acceleration_feedback_stabilizes_example(self):
        uncontrolled = dominant_ring_mode(8.0, 30, leader_acceleration_gain=0.0)
        controlled = dominant_ring_mode(8.0, 30, leader_acceleration_gain=0.25)
        self.assertGreater(uncontrolled["growth_rate_per_s"], 0.0)
        self.assertLess(controlled["growth_rate_per_s"], 0.0)

    def test_lwr_rankine_hugoniot_speed(self):
        result = solve_lwr_riemann()
        self.assertAlmostEqual(result["exact_shock_speed_mps"], -2.42536, places=4)
        self.assertLess(
            abs(result["measured_shock_speed_mps"] - result["exact_shock_speed_mps"]),
            result["dx_m"] / 60.0,
        )


if __name__ == "__main__":
    unittest.main()
