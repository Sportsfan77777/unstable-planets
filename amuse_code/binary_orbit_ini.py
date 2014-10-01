"""
binary on orbit with given parameters
"""

import numpy as np
import os
import math

from amuse.units.optparse import OptionParser
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.community.huayno.interface import Huayno
from amuse.community.kepler.interface import Kepler as Kepler_twobody
from amuse.io import write_set_to_file

from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import newton
from amuse.ext.orbital_elements import true_anomaly_from_eccentric_anomaly as true_from_ecc

pi_180 = np.pi/180.0

def orbital_parameters(rel_position, rel_velocity, total_mass):
    separation = rel_position.lengths()
    speed_squared = rel_velocity.lengths_squared()

    semimajor_axis = (constants.G * total_mass * separation / (2 * constants.G * total_mass - separation * speed_squared)).as_quantity_in(units.AU)
    eccentricity = numpy.sqrt( 1.0 - (rel_position.cross(rel_velocity)**2).sum(axis=0) / (constants.G * total_mass * semimajor_axis))
    period = (2 * numpy.pi * semimajor_axis**1.5 / (numpy.sqrt(constants.G * total_mass))).as_quantity_in(units.yr)
    
    return semimajor_axis, eccentricity, period
    
def true_anomaly_from_mean_anomaly(M,e):
    """ 
    converts the mean anomaly M to the true anomaly f 
    via the intermediate eccentric anomaly E that is calculated from eccentricity e
    """
    def k(E,e,M):
       """Kepler's Equation: solve for E"""
       return E - e*math.sin(E) - M 
    
    def dk(E,e,M):
	   """Derivative of Kepler's Equation: solve for E"""
	   return 1 - e*math.cos(E)
	
    guess = 3.0  # can be anything really (as far as I know)
    (E, flag) = newton(k, guess, fprime = dk, args = (e,M))
    
    print "Mean Anomaly:", M
    print "Eccentric Anomaly:", E, "(Flag:)", flag
	
    f = true_from_ecc(E, e)
    print "True Anomaly:", f
    return f
    
    
def kepler_orbital_elements_from_binary(two_bodies, angle = True):
    """ 
    Kepler orbital elements from two bodies 
    If angle, return angle. else return cos(angle) for angle orbital elements
    """
    if not angle:
       return (orbital_elements_from_binary(two_bodies))
    else:
       mass1, mass2, semimajor_axis, eccentricity, \
       true_anomaly, inclination, long_asc_node, arg_per \
          = orbital_elements_from_binary(two_bodies)
          
       return mass1, mass2, semimajor_axis, eccentricity, \
              math.acos(true_anomaly), math.acos(inclination), \
              math.acos(long_asc_node), math.acos(arg_per)
              
def kepler_binary_from_orbital_elements(m1, m2, sma, 
                        ecc, true_anom, inc, long_asc_node, arg_peri):
    """
    returns two particles (w/ coordinates in 'Center of Mass' frame)
    from given orbital elements
    """
    return new_binary_from_orbital_elements(m1, m2, sma, 
                        ecc, true_anom, inc, long_asc_node, arg_peri,
                        G = constants.G)
                        
def planet_from_orbital_elements(mass, sma, ecc, true_anom, inc, long_asc_node, arg_peri):
    """
    returns one particle (a planet) w/ coordinates in 'Center of Mass' frame of the binary
    This is not what you would expect from the typical binary_from_... function.
    """
    
    two_bodies = kepler_binary_from_orbital_elements(mass, 0.0 * mass, 
                        sma, ecc, true_anom, inc, long_asc_node, arg_peri)
                       
    # Move to 'Binary' Frame 
    two_bodies.position -= two_bodies[0].position
    two_bodies.velocity -= two_bodies[0].velocity
    
    planet = Particles(0)
    planet.add_particle(two_bodies[1])
    
    
    return planet


          
