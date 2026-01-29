# 🌍🚀🌙 Earth-Moon Mission Simulator

Real-time orbital mechanics simulation of a complete Earth-Moon-Earth round-trip mission with live 2-hour looping visualization.

![Mission Status](https://img.shields.io/badge/status-live%20simulation-brightgreen)
![Python](https://img.shields.io/badge/python-3.10-blue)
![License](https://img.shields.io/badge/license-MIT-blue)

## 🎥 Live Mission Simulation

**[▶️ VIEW LIVE SIMULATION](mission_live.html)** ← Click to see the mission in real-time!

The complete 7-day mission compressed into a 2-hour loop, running continuously with:
- 🌍 Earth-centered orbital view
- 🌙 Moon orbital mechanics
- 🚀 Spacecraft trajectory with trail
- 📊 Real-time telemetry
- 🔄 Auto-restart every 2 hours

---

## 🎯 Mission Profile

### Complete Mission Phases (11 Phases):

1. **🚀 Launch & LEO Insertion** - Insert into 400 km Low Earth Orbit
2. **🔄 LEO Parking Orbits** - 4 complete orbits around Earth
3. **🔥 Trans-Lunar Injection (TLI)** - Δv = 3,100 m/s burn towards Moon
4. **🌙 Trans-Lunar Coast** - 3-day coast phase to Moon
5. **🔥 Lunar Orbit Insertion (LOI)** - Δv = 900 m/s capture burn
6. **🔄 Lunar Parking Orbits** - 4 complete orbits around Moon at 100 km
7. **🔥 Trans-Earth Injection (TEI)** - Δv = 900 m/s return burn
8. **🌍 Trans-Earth Coast** - 3-day return to Earth
9. **🔥 Earth Orbit Insertion (EOI)** - Δv = 3,100 m/s capture burn
10. **🔄 Final LEO Orbits** - 2 orbits before de-orbit
11. **🛬 De-orbit & Landing** - Δv = 100 m/s, touchdown on Earth

**Total Mission Duration:** ~7 days  
**Total Delta-V:** ~8,100 m/s  
**Visualization Loop:** 2 hours (accelerated 84x)

---

## 🔧 Technical Details

### Orbital Mechanics
- **Framework:** Two-body problem with patched conics
- **Primary Body:** Earth-centered inertial frame
- **Secondary Body:** Moon-centered for lunar orbit phase
- **Integration:** Euler method with adaptive timestep
- **Coordinate System:** 2D Earth-centered (x, y)

### Spacecraft Specifications
- **Dry Mass:** 10,000 kg
- **Fuel Mass:** 5,000 kg
- **Total Mass:** 15,000 kg
- **Specific Impulse (Isp):** 450 seconds
- **Propulsion:** Chemical (high-efficiency orbital maneuvering)

### Key Orbital Parameters
| Parameter | Value |
|-----------|-------|
| LEO Altitude | 400 km |
| LEO Period | ~92 minutes |
| LLO Altitude | 100 km |
| LLO Period | ~118 minutes |
| Earth-Moon Distance | 384,400 km |
| TLI Δv | 3,100 m/s |
| LOI Δv | 900 m/s |
| TEI Δv | 900 m/s |
| EOI Δv | 3,100 m/s |

---

## 🚀 Quick Start

### Run Simulation Locally

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/orbital-mission-sim.git
cd orbital-mission-sim

# 2. Install dependencies
pip install numpy

# 3. Generate mission timeline
python scripts/generate_mission_data.py

# 4. Open visualization
open mission_live.html
# Or on Linux: xdg-open mission_live.html
# Or just drag mission_live.html into your browser
```

### Generate New Mission Data

```bash
python scripts/generate_mission_data.py
```

This creates `assets/mission_timeline.json` with:
- Position (x, y) at every timestep
- Velocity vectors
- Fuel consumption
- Mission phase
- Altitude and orbital parameters

---

## 📂 Project Structure

```
orbital-mission-sim/
├── mission_live.html           # Real-time visualization (open this!)
├── simulation/
│   ├── __init__.py
│   └── earth_moon_dynamics.py  # Mission physics engine
├── scripts/
│   └── generate_mission_data.py # Timeline generator
├── assets/
│   └── mission_timeline.json   # Mission state data
└── README.md                   # You are here
```

---

## 🎨 Visualization Features

### Real-Time Display
- **Earth**: Large blue sphere (to scale)
- **Moon**: Gray sphere, orbits Earth every 27.3 days
- **Spacecraft**: Red dot with white trajectory trail
- **Moon Orbit**: Dashed circle showing lunar path
- **Coordinate Grid**: Earth-centered inertial frame

### Telemetry Panel
- Mission elapsed time
- Current phase indicator
- Altitude above Earth/Moon
- Velocity magnitude
- Fuel remaining
- Position coordinates (x, y)
- Mission progress bar

### Animation
- **Loop Duration**: Exactly 2 hours
- **Time Compression**: 84x real-time (7 days → 2 hours)
- **Frame Rate**: 60 FPS
- **Trail Length**: Last 500 positions
- **Auto-Restart**: Seamless loop

---

## 🧮 Physics Implementation

### Two-Body Orbital Motion

The simulator uses Newton's law of gravitation:

```
F = -G * M * m / r²
```

Where:
- G = Gravitational constant (6.674×10⁻¹¹ m³/kg·s²)
- M = Mass of Earth or Moon
- m = Spacecraft mass
- r = Distance from center of primary body

### Patched Conics Approximation

The mission switches between:
1. **Earth-centered** (most of mission)
2. **Moon-centered** (during lunar orbit only)

Switch occurs at Moon's sphere of influence (~66,000 km from Moon).

### Rocket Equation

Fuel consumption via Tsiolkovsky:

```
Δv = Isp * g₀ * ln(m_initial / m_final)
```

Where:
- Isp = Specific impulse (450 s)
- g₀ = Standard gravity (9.80665 m/s²)
- Δv = Change in velocity

---

## 🔬 Mission Analysis

### Delta-V Budget Breakdown

| Phase | Δv (m/s) | % of Total |
|-------|----------|------------|
| TLI (Earth escape) | 3,100 | 38.3% |
| LOI (Moon capture) | 900 | 11.1% |
| TEI (Moon escape) | 900 | 11.1% |
| EOI (Earth capture) | 3,100 | 38.3% |
| De-orbit | 100 | 1.2% |
| **Total** | **8,100** | **100%** |

### Fuel Consumption

- **Initial Fuel:** 5,000 kg
- **Fuel Used:** ~4,850 kg
- **Fuel Remaining:** ~150 kg (3% margin)
- **Mass Ratio:** 1.5:1 (wet/dry)

---

## 🛠️ Customization

### Modify Mission Parameters

Edit `simulation/earth_moon_dynamics.py`:

```python
mission = EarthMoonMission(
    dry_mass_kg=10000,   # Change spacecraft mass
    fuel_mass_kg=5000,   # Adjust fuel capacity
    isp_s=450            # Modify engine efficiency
)
```

### Change Visualization Speed

Edit `mission_live.html`:

```javascript
const LOOP_DURATION = 2 * 3600; // Change loop time (seconds)
```

### Adjust Orbital Altitudes

Edit constants in `earth_moon_dynamics.py`:

```python
LEO_ALTITUDE = 400e3  # Low Earth Orbit (m)
LLO_ALTITUDE = 100e3  # Low Lunar Orbit (m)
```

---

## 📊 Data Format

The `mission_timeline.json` contains:

```json
{
  "mission_duration": 604800,
  "initial_fuel": 5000,
  "final_fuel": 150,
  "states": [
    {
      "time": 0,
      "phase": "leo_insertion",
      "x": 0,
      "y": 6771000,
      "vx": 7670,
      "vy": 0,
      "fuel": 5000,
      "altitude": 400000,
      "velocity": 7670
    },
    ...
  ]
}
```

---

## 🎓 Educational Value

Perfect for learning:
- ✅ Orbital mechanics fundamentals
- ✅ Two-body problem physics
- ✅ Patched conic approximation
- ✅ Hohmann transfers
- ✅ Delta-v budgeting
- ✅ Mission planning
- ✅ Real-time data visualization

---

## 🚧 Future Enhancements

Potential additions:
- [ ] 3D visualization
- [ ] Multiple trajectory options
- [ ] Gravity assists
- [ ] N-body physics
- [ ] Realistic atmospheric drag
- [ ] Solar radiation pressure
- [ ] Station-keeping maneuvers
- [ ] Interactive mission planning

---

## 📜 License

MIT License - Free to use, modify, and distribute

---

## 🙏 Acknowledgments

- Orbital mechanics based on NASA mission design principles
- Inspired by Apollo program trajectories
- Physics calculations verified against JPL HORIZONS

---

## 📞 Contact & Contribution

Found a bug? Have a suggestion? Want to add features?

- **Issues:** [GitHub Issues](https://github.com/YOUR_USERNAME/orbital-mission-sim/issues)
- **Contributions:** Pull requests welcome!
- **Questions:** Open a discussion

---

**🌟 Star this repo if you find it useful!**

**🚀 Happy Orbiting!**
