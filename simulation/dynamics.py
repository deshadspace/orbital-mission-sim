"""
Asteroid Mission Simulator using Poliastro

NOTE: Poliastro 0.17.0 requires Python >=3.8,<3.11
      For Python 3.11+, use poliastro 0.18+ (currently in development)
      or downgrade Python to 3.10

Installation:
    pip install numpy matplotlib astropy poliastro

For Python 3.12+:
    Consider using a virtual environment with Python 3.10:
    conda create -n asteroid python=3.10
    conda activate asteroid
    pip install numpy matplotlib astropy poliastro
"""

import numpy as np
from typing import Dict, Callable

from astropy import units as u
from astropy.time import Time

from poliastro.bodies import Sun, Earth
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.iod import izzo
from poliastro.ephem import Ephem

# ===============================
# Monte Carlo + Polynomial Chaos
# ===============================

class UncertaintyPropagator:
    """
    Generic uncertainty propagator using Monte Carlo and Hermite PCE.
    """

    def __init__(
        self,
        param_distributions: Dict[str, Dict],
        model_fn: Callable[[Dict[str, float]], float],
        pce_order: int = 2,
        num_samples: int = 1000,
    ):
        self.param_distributions = param_distributions
        self.model_fn = model_fn
        self.pce_order = pce_order
        self.num_samples = num_samples

        self.sampled_inputs = []
        self.mc_outputs = []
        self.pce_coeffs = {}

    def sample_inputs(self):
        samples = []
        for _ in range(self.num_samples):
            s = {}
            for name, dist in self.param_distributions.items():
                if dist["type"] == "normal":
                    s[name] = np.random.normal(dist["mean"], dist["std"])
                else:
                    raise NotImplementedError
            samples.append(s)
        self.sampled_inputs = samples

    def run_monte_carlo(self):
        self.mc_outputs = [self.model_fn(s) for s in self.sampled_inputs]

    def fit_pce(self):
        for name, dist in self.param_distributions.items():
            xi = np.array([(s[name] - dist["mean"]) / dist["std"] for s in self.sampled_inputs])
            basis = np.vstack([xi**i for i in range(self.pce_order + 1)]).T
            y = np.array(self.mc_outputs)
            coeffs, *_ = np.linalg.lstsq(basis, y, rcond=None)
            self.pce_coeffs[name] = coeffs

    def summary(self):
        return {
            "mean": float(np.mean(self.mc_outputs)),
            "std": float(np.std(self.mc_outputs)),
            "samples": self.num_samples,
            "pce_order": self.pce_order,
        }

# =========================================
# Heliocentric Earth → Asteroid → Earth
# =========================================

class AsteroidMissionSimulator:
    """
    High-fidelity heliocentric mission simulator with fuel tracking.
    """

    def __init__(
        self,
        epoch_launch="2028-01-01",
        dry_mass_kg=900,
        fuel_mass_kg=600,
        isp_s=320,
        asteroid_name="99942 Apophis",
    ):
        self.epoch = Time(epoch_launch, scale="tdb")

        # Spacecraft properties
        self.dry_mass = dry_mass_kg * u.kg
        self.fuel_mass = fuel_mass_kg * u.kg
        self.isp = isp_s * u.s
        self.g0 = 9.80665 * u.m / u.s**2

        self.total_mass = self.dry_mass + self.fuel_mass

        # Initial heliocentric orbit = Earth state
        earth_ephem = Ephem.from_body(Earth, self.epoch)
        self.orbit = Orbit.from_ephem(Sun, earth_ephem, self.epoch)

        # Asteroid ephemeris
        self.asteroid_ephem = Ephem.from_horizons(
            asteroid_name,
            epochs=self.epoch,
            attractor=Sun,
        )

        self.dv_used = 0 * u.m / u.s

        print(f"[Init] Mission initialized at {self.epoch.iso}")
        print(f"[Init] Initial mass: {self.total_mass:.1f}")

    # ------------------------
    # Rocket equation
    # ------------------------

    def _consume_fuel(self, delta_v):
        m0 = self.total_mass
        mf = m0 / np.exp(delta_v / (self.isp * self.g0))
        fuel_used = m0 - mf

        if fuel_used > self.fuel_mass:
            raise RuntimeError("Out of fuel!")

        self.fuel_mass -= fuel_used
        self.total_mass = self.dry_mass + self.fuel_mass
        self.dv_used += delta_v

    # ------------------------
    # Lambert transfer
    # ------------------------

    def lambert_transfer(self, target_orbit, tof_days):
        r1, v1 = self.orbit.r, self.orbit.v
        r2 = target_orbit.r
        tof = tof_days * u.day

        v_depart, v_arrive = izzo.lambert(Sun.k, r1, r2, tof)

        dv_vec = v_depart - v1
        dv = np.linalg.norm(dv_vec)

        self._consume_fuel(dv)
        self.orbit = self.orbit.apply_maneuver(Maneuver.impulse(dv_vec))

        self.orbit = self.orbit.propagate(tof)

        return dv.to(u.m / u.s)

    # ------------------------
    # Mission phases
    # ------------------------

    def go_to_asteroid(self, tof_days=300):
        asteroid_orbit = Orbit.from_ephem(
            Sun,
            self.asteroid_ephem,
            epoch=self.epoch + tof_days * u.day,
        )
        dv = self.lambert_transfer(asteroid_orbit, tof_days)
        print(f"[Burn] Earth → Asteroid Δv = {dv:.1f}")
        return dv

    def return_to_earth(self, tof_days=300):
        earth_ephem = Ephem.from_body(Earth, self.epoch + tof_days * u.day)
        earth_orbit = Orbit.from_ephem(
            Sun,
            earth_ephem,
            self.epoch + tof_days * u.day,
        )
        dv = self.lambert_transfer(earth_orbit, tof_days)
        print(f"[Burn] Asteroid → Earth Δv = {dv:.1f}")
        return dv

    # ------------------------
    # Telemetry
    # ------------------------

    def state(self):
        return {
            "epoch": self.orbit.epoch.iso,
            "position_au": self.orbit.r.to(u.au).value.tolist(),
            "velocity_kms": self.orbit.v.to(u.km / u.s).value.tolist(),
            "fuel_remaining_kg": self.fuel_mass.to(u.kg).value,
            "delta_v_used_ms": self.dv_used.to(u.m / u.s).value,
        }

# ===============================
# Example Monte Carlo usage
# ===============================

def mission_delta_v_model(params):
    sim = AsteroidMissionSimulator(
        epoch_launch="2028-01-01",
        fuel_mass_kg=600 * (1 + params["fuel_bias"]),
        isp_s=320 * (1 + params["isp_bias"]),
    )

    sim.go_to_asteroid(300)
    sim.return_to_earth(300)

    return sim.dv_used.to(u.m / u.s).value