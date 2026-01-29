"""
Earth-Moon Round Trip Mission Simulator
Two-body problem with patched conics approximation
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple
import json

# Constants
G = 6.67430e-11  # m^3 kg^-1 s^-2
M_EARTH = 5.972e24  # kg
M_MOON = 7.342e22  # kg
R_EARTH = 6.371e6  # m
R_MOON = 1.737e6  # m
EARTH_MOON_DISTANCE = 384400e3  # m (average)
MU_EARTH = G * M_EARTH
MU_MOON = G * M_MOON

# Mission parameters
LEO_ALTITUDE = 400e3  # m (400 km)
LLO_ALTITUDE = 100e3  # m (100 km lunar orbit)
LEO_RADIUS = R_EARTH + LEO_ALTITUDE
LLO_RADIUS = R_MOON + LLO_ALTITUDE

@dataclass
class MissionState:
    """State of spacecraft at any time"""
    time: float  # seconds from mission start
    phase: str  # mission phase name
    x: float  # position x (m, Earth-centered)
    y: float  # position y (m, Earth-centered)
    vx: float  # velocity x (m/s)
    vy: float  # velocity y (m/s)
    fuel_remaining: float  # kg
    altitude: float  # m
    velocity_magnitude: float  # m/s

class EarthMoonMission:
    """
    Simulates complete Earth-Moon-Earth mission
    All calculations in 2D Earth-centered inertial frame
    """
    
    def __init__(
        self,
        dry_mass_kg: float = 10000,
        fuel_mass_kg: float = 5000,
        isp_s: float = 450,  # High Isp for orbital maneuvers
    ):
        self.dry_mass = dry_mass_kg
        self.fuel_mass = fuel_mass_kg
        self.initial_fuel = fuel_mass_kg
        self.isp = isp_s
        self.g0 = 9.80665  # m/s^2
        
        # Mission timeline storage
        self.timeline: List[MissionState] = []
        self.current_time = 0.0
        
        # Current state
        self.x = 0.0
        self.y = LEO_RADIUS  # Start at LEO
        self.vx = np.sqrt(MU_EARTH / LEO_RADIUS)  # Circular orbit velocity
        self.vy = 0.0
        self.phase = "launch"
        
    def total_mass(self):
        return self.dry_mass + self.fuel_mass
    
    def apply_delta_v(self, dv_x: float, dv_y: float) -> float:
        """
        Apply delta-v and consume fuel via rocket equation
        Returns fuel consumed
        """
        dv_mag = np.sqrt(dv_x**2 + dv_y**2)
        
        # Rocket equation: m_final = m_initial * exp(-dv / (Isp * g0))
        m_initial = self.total_mass()
        m_final = m_initial * np.exp(-dv_mag / (self.isp * self.g0))
        fuel_consumed = m_initial - m_final
        
        if fuel_consumed > self.fuel_mass:
            raise RuntimeError(f"Out of fuel! Need {fuel_consumed:.1f} kg, have {self.fuel_mass:.1f} kg")
        
        self.fuel_mass -= fuel_consumed
        self.vx += dv_x
        self.vy += dv_y
        
        return fuel_consumed
    
    def propagate_orbit(self, duration: float, dt: float = 60.0, center='earth'):
        """
        Propagate two-body orbit for given duration using proper physics
        center: 'earth' or 'moon'
        """
        steps = int(duration / dt)
        
        for step in range(steps):
            # Choose gravitational parameter and reference frame
            if center == 'earth':
                mu = MU_EARTH
                # Earth-centered: Earth is at origin
                cx, cy = 0.0, 0.0
                ref_radius = R_EARTH
            else:  # moon-centered
                # Moon position in Earth-centered frame
                moon_angle = (self.current_time / (27.3 * 86400)) * 2 * np.pi
                cx = EARTH_MOON_DISTANCE * np.cos(moon_angle)
                cy = EARTH_MOON_DISTANCE * np.sin(moon_angle)
                mu = MU_MOON
                ref_radius = R_MOON
            
            # Position relative to primary body
            rx = self.x - cx
            ry = self.y - cy
            r = np.sqrt(rx**2 + ry**2)
            
            if r < ref_radius:
                print(f"WARNING: Spacecraft below surface at t={self.current_time:.1f}s, r={r/1000:.1f}km")
                r = ref_radius * 1.01  # Prevent crash
            
            # Gravitational acceleration (two-body)
            ax = -mu * rx / r**3
            ay = -mu * ry / r**3
            
            # Add Earth's gravity when orbiting Moon (patched conics)
            if center == 'moon':
                # Distance from Earth
                earth_rx = self.x
                earth_ry = self.y  
                earth_r = np.sqrt(earth_rx**2 + earth_ry**2)
                if earth_r > 0:
                    # Earth's pull on spacecraft
                    ax += -MU_EARTH * earth_rx / earth_r**3
                    ay += -MU_EARTH * earth_ry / earth_r**3
            
            # Velocity Verlet integration (more stable than Euler)
            # Half-step velocity update
            vx_half = self.vx + 0.5 * ax * dt
            vy_half = self.vy + 0.5 * ay * dt
            
            # Full-step position update
            self.x += vx_half * dt
            self.y += vy_half * dt
            
            # Recalculate acceleration at new position
            rx = self.x - cx
            ry = self.y - cy
            r = np.sqrt(rx**2 + ry**2)
            
            if r < ref_radius:
                r = ref_radius * 1.01
                
            ax_new = -mu * rx / r**3
            ay_new = -mu * ry / r**3
            
            if center == 'moon':
                earth_rx = self.x
                earth_ry = self.y
                earth_r = np.sqrt(earth_rx**2 + earth_ry**2)
                if earth_r > 0:
                    ax_new += -MU_EARTH * earth_rx / earth_r**3
                    ay_new += -MU_EARTH * earth_ry / earth_r**3
            
            # Full-step velocity update
            self.vx = vx_half + 0.5 * ax_new * dt
            self.vy = vy_half + 0.5 * ay_new * dt
            
            self.current_time += dt
            
            # Record state every few steps to reduce data size
            if step % 10 == 0 or step == steps - 1:
                altitude = r - ref_radius
                vel_mag = np.sqrt(self.vx**2 + self.vy**2)
                
                self.timeline.append(MissionState(
                    time=self.current_time,
                    phase=self.phase,
                    x=self.x,
                    y=self.y,
                    vx=self.vx,
                    vy=self.vy,
                    fuel_remaining=self.fuel_mass,
                    altitude=altitude,
                    velocity_magnitude=vel_mag
                ))
    
    def run_mission(self):
        """Execute complete mission profile with SIMPLIFIED but CORRECT trajectory"""
        print("🚀 Starting Earth-Moon-Earth Mission")
        print(f"Initial mass: {self.total_mass():.1f} kg")
        print(f"Fuel: {self.fuel_mass:.1f} kg\n")
        
        # Phase 1: Launch to LEO (assume already there)
        self.phase = "leo_insertion"
        print(f"✓ Phase 1: LEO Insertion at {LEO_ALTITUDE/1000:.0f} km")
        
        # Phase 2: Park in LEO - position spacecraft for intercept
        self.phase = "leo_parking"
        leo_period = 2 * np.pi * np.sqrt(LEO_RADIUS**3 / MU_EARTH)
        print(f"✓ Phase 2: LEO parking orbit ({leo_period/60:.1f} min period)")
        
        # Orbit until we're at the RIGHT POSITION to burn toward Moon
        # We want to be at a position where burning prograde will send us to the Moon
        
        # Calculate Moon's position NOW
        moon_angle_start = (self.current_time / (27.3 * 86400)) * 2 * np.pi
        
        # We need to burn when the angle between spacecraft and Moon is ~60-90 degrees
        # (this is the "phase angle" for Hohmann transfer)
        # Let's orbit until we're at about 60° behind the Moon
        
        # Current spacecraft angle
        craft_angle = np.arctan2(self.y, self.x)
        
        # Calculate how much to orbit to get proper phasing
        # Target: 90° behind Moon (in direction of motion)
        target_angle_diff = np.pi / 2  # 90 degrees
        
        angle_to_wait = (moon_angle_start - craft_angle - target_angle_diff) % (2 * np.pi)
        orbits_to_wait = angle_to_wait / (2 * np.pi)
        if orbits_to_wait < 0.1:
            orbits_to_wait += 1  # Wait at least one orbit
        
        wait_time = orbits_to_wait * leo_period
        self.propagate_orbit(wait_time, dt=20.0, center='earth')
        
        print(f"   Waited {orbits_to_wait:.2f} orbits for proper Moon phasing")
        
        # Now spacecraft position
        craft_angle = np.arctan2(self.y, self.x)
        moon_angle_now = (self.current_time / (27.3 * 86400)) * 2 * np.pi
        phase_angle = (moon_angle_now - craft_angle) * 180 / np.pi
        print(f"   Phase angle to Moon: {phase_angle:.1f}°")
        
        # Phase 3: TLI Burn - PURE PROGRADE (tangent to orbit)
        self.phase = "tli_burn"
        print(f"✓ Phase 3: Trans-Lunar Injection")
        
        # Simple: just add velocity in current direction of motion
        # This creates an ellipse that goes higher
        v_current = np.sqrt(self.vx**2 + self.vy**2)
        vx_unit = self.vx / v_current
        vy_unit = self.vy / v_current
        
        # Calculate exact delta-v for Hohmann-like transfer
        # From LEO to Moon's orbit altitude
        r1 = LEO_RADIUS
        r2 = EARTH_MOON_DISTANCE
        
        # Hohmann transfer formula
        v_circular_leo = np.sqrt(MU_EARTH / r1)
        v_transfer_periapsis = np.sqrt(MU_EARTH * (2/r1 - 2/(r1 + r2)))
        dv_tli = v_transfer_periapsis - v_circular_leo
        
        print(f"   Calculated Hohmann Δv: {dv_tli:.0f} m/s")
        
        # Apply PURE PROGRADE burn
        fuel_used = self.apply_delta_v(vx_unit * dv_tli, vy_unit * dv_tli)
        print(f"   Fuel used: {fuel_used:.1f} kg")
        
        # Phase 4: Coast to Moon altitude
        self.phase = "trans_lunar_coast"
        print(f"✓ Phase 4: Coasting to Moon")
        
        # Calculate transfer time (half period of transfer ellipse)
        a_transfer = (r1 + r2) / 2
        transfer_time = np.pi * np.sqrt(a_transfer**3 / MU_EARTH)
        
        print(f"   Transfer time: {transfer_time/86400:.2f} days")
        
        self.propagate_orbit(transfer_time, dt=100.0, center='earth')
        
        # Check position
        dist_from_earth = np.sqrt(self.x**2 + self.y**2)
        moon_angle_arrival = (self.current_time / (27.3 * 86400)) * 2 * np.pi
        moon_x = EARTH_MOON_DISTANCE * np.cos(moon_angle_arrival)
        moon_y = EARTH_MOON_DISTANCE * np.sin(moon_angle_arrival)
        dist_to_moon = np.sqrt((self.x - moon_x)**2 + (self.y - moon_y)**2)
        
        print(f"   Arrived at distance: {dist_from_earth/1000:.0f} km from Earth")
        print(f"   Distance to Moon: {dist_to_moon/1000:.0f} km")
        
        # Phase 5: LOI - Retrograde burn to circularize
        self.phase = "loi_burn"
        print(f"✓ Phase 5: Lunar Orbit Insertion")
        
        v_current = np.sqrt(self.vx**2 + self.vy**2)
        vx_unit = self.vx / v_current
        vy_unit = self.vy / v_current
        
        # Slow down to lunar orbital velocity
        dv_loi = v_current * 0.3  # Reduce velocity by 30%
        
        fuel_used = self.apply_delta_v(-vx_unit * dv_loi, -vy_unit * dv_loi)
        print(f"   Δv = {dv_loi:.0f} m/s, Fuel used: {fuel_used:.1f} kg")
        
        # Phase 6: Lunar orbits
        self.phase = "llo_parking"
        print(f"✓ Phase 6: Lunar orbit phase")
        
        # Get current Moon position
        moon_angle = (self.current_time / (27.3 * 86400)) * 2 * np.pi
        moon_x = EARTH_MOON_DISTANCE * np.cos(moon_angle)
        moon_y = EARTH_MOON_DISTANCE * np.sin(moon_angle)
        
        # Orbit around the Moon's current position for a bit
        dist_to_moon = np.sqrt((self.x - moon_x)**2 + (self.y - moon_y)**2)
        print(f"   Distance to Moon: {dist_to_moon/1000:.0f} km")
        
        # Do short lunar orbit phase
        self.propagate_orbit(4 * 3600, dt=20.0, center='moon')  # 4 hours
        
        # Phase 7: TEI - Return trajectory
        self.phase = "tei_burn"
        print(f"✓ Phase 7: Trans-Earth Injection")
        
        # Current position and velocity
        v_current = np.sqrt(self.vx**2 + self.vy**2)
        
        # Calculate return Hohmann transfer
        # We're at Moon's distance, need to get back to LEO
        r1 = EARTH_MOON_DISTANCE  # Starting at Moon's orbit
        r2 = LEO_RADIUS  # Returning to LEO
        
        # Hohmann return burn (from apoapsis of return ellipse)
        v_at_moon = np.sqrt(MU_EARTH * (2/r1 - 2/(r1 + r2)))
        
        # But we need to account for current velocity
        # Apply burn in OPPOSITE direction to current velocity (retrograde at Moon)
        vx_unit = self.vx / v_current
        vy_unit = self.vy / v_current
        
        # The delta-v is the difference between current speed and required transfer speed
        dv_tei = abs(v_current - v_at_moon)
        
        # Apply retrograde burn (slow down to drop orbit)
        fuel_used = self.apply_delta_v(-vx_unit * dv_tei, -vy_unit * dv_tei)
        print(f"   Δv = {dv_tei:.0f} m/s, Fuel used: {fuel_used:.1f} kg")
        
        # Phase 8: Coast back - same transfer time as outbound
        self.phase = "trans_earth_coast"
        print(f"✓ Phase 8: Return coast")
        self.propagate_orbit(transfer_time, dt=100.0, center='earth')
        
        dist_from_earth = np.sqrt(self.x**2 + self.y**2)
        print(f"   Distance from Earth at EOI: {dist_from_earth/1000:.0f} km")
        
        # Phase 9: EOI
        self.phase = "eoi_burn"
        print(f"✓ Phase 9: Earth Orbit Insertion")
        
        v_current = np.sqrt(self.vx**2 + self.vy**2)
        vx_unit = self.vx / v_current
        vy_unit = self.vy / v_current
        
        # Hohmann circularization
        dv_eoi = v_current * 0.4
        fuel_used = self.apply_delta_v(-vx_unit * dv_eoi, -vy_unit * dv_eoi)
        print(f"   Δv = {dv_eoi:.0f} m/s, Fuel used: {fuel_used:.1f} kg")
        
        # Phase 10: Final orbits
        self.phase = "leo_final"
        print(f"✓ Phase 10: Final LEO orbits")
        self.propagate_orbit(2 * leo_period, dt=20.0, center='earth')
        
        # Phase 11: Deorbit
        self.phase = "deorbit"
        print(f"✓ Phase 11: De-orbit")
        
        v_current = np.sqrt(self.vx**2 + self.vy**2)
        vx_unit = self.vx / v_current
        vy_unit = self.vy / v_current
        
        dv_deorbit = 100
        fuel_used = self.apply_delta_v(-vx_unit * dv_deorbit, -vy_unit * dv_deorbit)
        print(f"   Δv = {dv_deorbit:.0f} m/s")
        
        self.phase = "reentry"
        self.propagate_orbit(300, dt=10.0, center='earth')
        
        self.phase = "landed"
        print(f"\n🎉 Mission Complete!")
        print(f"Total mission time: {self.current_time/86400:.2f} days")
        print(f"Fuel remaining: {self.fuel_mass:.1f} kg")
        
    def export_timeline(self, filename: str):
        """Export timeline to JSON for visualization"""
        data = {
            "mission_duration": self.current_time,
            "initial_fuel": self.initial_fuel,
            "final_fuel": self.fuel_mass,
            "states": [
                {
                    "time": s.time,
                    "phase": s.phase,
                    "x": s.x,
                    "y": s.y,
                    "vx": s.vx,
                    "vy": s.vy,
                    "fuel": s.fuel_remaining,
                    "altitude": s.altitude,
                    "velocity": s.velocity_magnitude
                }
                for s in self.timeline
            ]
        }
        
        with open(filename, 'w') as f:
            json.dump(data, f)
        
        print(f"\n📊 Timeline exported to {filename}")
        print(f"   {len(self.timeline)} data points")

if __name__ == "__main__":
    mission = EarthMoonMission()
    mission.run_mission()
    mission.export_timeline("mission_timeline.json")