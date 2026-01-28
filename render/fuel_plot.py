# render/fuel_plot.py

import matplotlib.pyplot as plt
from simulation.dynamics import AsteroidMissionSimulator


def render_fuel():
    sim = AsteroidMissionSimulator()
    fuel = []

    sim.go_to_asteroid(300)
    fuel.append(sim.fuel_mass.value)

    sim.return_to_earth(300)
    fuel.append(sim.fuel_mass.value)

    plt.figure()
    plt.plot(["After outbound", "After return"], fuel, marker="o")
    plt.ylabel("Fuel remaining (kg)")
    plt.title("Fuel Consumption")

    plt.savefig("assets/fuel.svg", format="svg")
    plt.close()


if __name__ == "__main__":
    render_fuel()