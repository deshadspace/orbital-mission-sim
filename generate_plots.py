
"""
Generate all mission visualization plots.
Run this script from the project root directory.
"""

import sys
import os

# Add project root to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Now import and run the render scripts
print(" Generating mission visualizations...\n")

print(" Rendering trajectory plot...")
from render.trajectory import render_trajectory
render_trajectory()
print("✓ Trajectory plot saved to assets/trajectory.svg\n")

print("Rendering fuel consumption plot...")
from render.fuel_plot import render_fuel
render_fuel()
print(" Fuel plot saved to assets/fuel.svg\n")

print("3/3 Rendering velocity profile plot...")
from render.velocity_plot import render_velocity
render_velocity()
print(" Velocity plot saved to assets/velocity.svg\n")

print(" All visualizations generated successfully!")
print("Check the assets/ folder for SVG files.")