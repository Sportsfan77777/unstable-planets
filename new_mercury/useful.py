"""
Useful functions that will be useful in plotting and in general

(1) Convert single-star orbital elements to xyz coordinates
(2) Convert xyz coordinates to double-star orbital elements
(3) By transitive property, convert single-star orbital element to double-star orbital elements
(4) Convert single-star orbital elements to companion-star orbital elements

"""

import numpy as np

import avatar
import amon

from structures import *


""" Save initial conditions of planets in pickle files! (not here...) """


def convert_to_elements(position, velocity, mu):
    """
    Given coordinates, convert to the orbital elements 
  
    Input
    Position: (x,y,z), Velocity: (vx,vy,vz), mu = GM
  
    Returns
    OrbitalElements: (a, e, i, argw, node, M)
    """

    orbital = OrbitalElements() # Return orbital
    h = position.cross(velocity, ret_type = AngularMomentum)

    # (1) semi-major axis
    sm_axis = avatar.calculate_a(position, velocity, mu)
    orbital.set_a(sm_axis)

    # (2) eccentricity 
    ecc = avatar.calculate_ecc(position, velocity, mu, h)
    orbital.set_ecc(ecc)

    # (3) inclination
    inc = avatar.calculate_inc(h)
    orbital.set_inc(inc)

    # (5) node
    node = avatar.calculate_node(elements, h)
    orbital.set_node(node)

    # (4) argument of pericenter
    argw = avatar.calculate_argw(elements, position, velocity, h)
    orbital.set_argw(argw)
    
    # (6) Mean Anomaly
    mean_anom = avatar.calculate_mean_anom(elements, position, velocity, h)
    orbital.set_mean_anom(mean_anom)
    
    return orbital

