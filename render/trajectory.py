# render/trajectory.py

import matplotlib.pyplot as plt
from astropy import units as u

from simulation.dynamics import AsteroidMissionSimulator


def render_trajectory():
    sim = AsteroidMissionSimulator()
    positions = []

    # Earth → Asteroid
    sim.go_to_asteroid(300)
    positions.append(sim.orbit.r.to(u.au).value)

    # Asteroid → Earth
    sim.return_to_earth(300)
    positions.append(sim.orbit.r.to(u.au).value)

    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]

    plt.figure(figsize=(6, 6))
    plt.plot(xs, ys, "-o", label="Spacecraft")
    plt.scatter(0, 0, c="orange", label="Sun")
    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.legend()
    plt.title("Earth → Asteroid → Earth Trajectory")
    plt.axis("equal")

    plt.savefig("assets/trajectory.svg", format="svg")
    plt.close()


if __name__ == "__main__":
    render_trajectory()