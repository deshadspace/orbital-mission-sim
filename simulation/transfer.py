import numpy as np
from astropy import units as u
from poliastro.iod import izzo
from poliastro.maneuver import Maneuver
from poliastro.bodies import Sun

def lambert_transfer(current_orbit, target_orbit, tof_days):
    r1, v1 = current_orbit.r, current_orbit.v
    r2 = target_orbit.r
    tof = tof_days * u.day

    v_depart, v_arrive = izzo.lambert(Sun.k, r1, r2, tof)

    dv_vec = v_depart - v1
    dv = np.linalg.norm(dv_vec)

    new_orbit = current_orbit.apply_maneuver(
        Maneuver.impulse(dv_vec)
    ).propagate(tof)

    return new_orbit, dv.to(u.m / u.s)
