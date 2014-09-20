"""
'struct' for storing orbital elements
"""

class Orb_Kepler:
    """ stores 6 Keplerian orbital elements (things that aren't cartesian coordinates) """
    def __init__(self, a, e, i, M, argw, node):
       self.a = a  # semi-major axis
       self.e = e  # eccentricity
       self.i = i  # inclination
       self.M = M  # mean anomaly
       self.argw = self.argw  # argument of pericenter
       self.node = self.node  # longitude of the ascending node
    
class Orb_Cartesian:
    """ 
    stores 6 orbital elements (things that actually are cartesian coordinates) 
    Note: with AMUSE, this should rarely be necessary
    """
    def __init__(self, x, y, z, vx, argw, node):
       self.x = x  # semi-major axis
       self.y = y  # eccentricity
       self.z = z  # inclination
       self.vx = vx  # mean anomaly
       self.vy = self.vy  # argument of pericenter
       self.vz = self.vz  # longitude of the ascending node