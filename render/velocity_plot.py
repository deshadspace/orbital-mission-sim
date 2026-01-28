# render/velocity_plot.py

import matplotlib.pyplot as plt
from astropy import units as u
import numpy as np

from simulation.dynamics import AsteroidMissionSimulator


def render_velocity():
    sim = AsteroidMissionSimulator()
    velocities = []

    sim.go_to_asteroid(300)
    velocities.append(np.linalg.norm(sim.orbit.v).to(u.km / u.s).value)

    sim.return_to_earth(300)
    velocities.append(np.linalg.norm(sim.orbit.v).to(u.km / u.s).value)

    plt.figure()
    plt.plot(["Outbound", "Inbound"], velocities, marker="o")
    plt.ylabel("Velocity (km/s)")
    plt.title("Heliocentric Velocity Profile")

    plt.savefig("assets/velocity.svg", format="svg")
    plt.close()


if __name__ == "__main__":
    render_velocity()