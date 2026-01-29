#!/usr/bin/env python3
"""
Generate mission timeline data for visualization
Run this to create mission_timeline.json
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from simulation.earth_moon_dynamics import EarthMoonMission

def main():
    print("="*60)
    print("GENERATING EARTH-MOON MISSION TIMELINE")
    print("="*60)
    
    # Create mission
    mission = EarthMoonMission(
        dry_mass_kg=50,
        fuel_mass_kg=656.0,
        isp_s=450
    )
    
    # Run complete mission
    mission.run_mission()
    
    # Export timeline
    output_path = "assets/mission_timeline.json"
    os.makedirs("assets", exist_ok=True)
    mission.export_timeline(output_path)
    
    print(f"\n✅ Mission timeline generated successfully!")
    print(f"📁 Saved to: {output_path}")
    print(f"📊 Data points: {len(mission.timeline)}")
    print(f"⏱️  Mission duration: {mission.current_time/3600:.1f} hours")

if __name__ == "__main__":
    main()