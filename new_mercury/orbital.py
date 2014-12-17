import numpy as np
#from pylab import *   ??????

from structures import *
import avatar
import amon

#set unit system
UnitMass_in_g            = 1.98892e33 # solar mass
UnitVelocity_in_cm_per_s = 4.74057581e5 # AU per year
UnitLength_in_cm         = 1.49598e13  # AU
G                        = 39.4751488


"""###################################################################"""

#functions
fac1 = 1.0

"""###############################################################################"""

def elements(position, velocity, mu):
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
    orbital.set_a()

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
    argw = avatar.calculate_argw()
    orbital.set_argw(argw)
    
    # (6) Mean Anomaly
    mean_anom = avatar.calculate_mean_anom()
    orbital.set_mean_anom(mean_anom)
    
    return orbital

"""#############################################################################"""

def convert_to_xyz(elements, mu):
    """
    Given the orbital elements, convert to coordinates
    
    Input
    (1) peri = distance to peri-center
    (2) e = eccentricity
    (3) I = inclination
    (4) omega = argument of peri-center
    (5) Omega = longitude of the ascending node
    (6) M = mean anomaly
    (7) mu = GM
    
    Returns
    (1,2,3) Position: (X,Y,Z)
    (4,5,6) Velocity: (VX,VY,VZ)
    """
 
    if (e < 1): # Elliptic
        ecc_anom = KeplerEquation(M,e)
        x = a * (np.cos(ecc_anom) - e)
        y = a * np.sqrt(1 - e * e) * np.sin(ecc_anom)
        xdot = -np.sqrt(mu/a) / (1.0 - e*np.cos(ecc_anom)) * np.sin(ecc_anom)
        ydot = np.sqrt(mu/a) / (1.0 - e*np.cos(ecc_anom)) * np.cos(ecc_anom) * np.sqrt(1.0 - e * e)

        flat_position =
        flat_velocity = 
    
    elif (e > 1): # Hyperbolic
        hyp_anom = HyperbolicKeplerEquation(M,e)
        x = a * (np.cosh(hyp_anom) - e)
        y = -a * np.sqrt(e * e - 1.0) * np.sinh(hyp_anom)
        xdot = -np.sqrt(mu/-a) / (e*np.cosh(hyp_anom)- 1.0) * np.sinh(hyp_anom)
        ydot = np.sqrt(mu/-a) / (e*np.cosh(hyp_anom) - 1.0) * np.cosh(hyp_anom) * np.sqrt(e * e - 1.0)

        flat_position = 
        flat_velocity = 
    
    else: # Parabolic (e = 1)
        tan_trueanom_2 = ParabolicKeplerEquation(M)
        x = peri * (1.0 - tan_trueanom_2**2)
        y = 2.0* peri * tan_trueanom_2
        xdot = -np.sqrt(2.0 * mu / peri) / (1.0 + tan_trueanom_2**2) * tan_trueanom_2
        ydot =  np.sqrt(2.0 * mu / peri) / (1.0 + tan_trueanom_2**2)

    flat_position = calculate_flat_position(elements, mu)
    flat_velocity = calculate_flat_velocity(elements, mu)

    position = rotate_xy_into_xyz(elements, flat_position, ret_type = Position)
    velocity = rotate_xy_into_xyz(elements, flat_velocity, ret_type = Velocity)
        
    return position, velocity


"""#############################################################################"""

def KeplerEquation(mean_anom, e):
    """
    Finds the eccentric anomaly 'ecc_anom' 
    given the mean anomaly 'mean_anom' and eccentricity 'e'
    """
    while (np.abs(mean_anom) > 2.0*np.pi):
    mean_anom-=2.0*np.pi*np.sign(mean_anom)
    if (mean_anom < 0.0): mean_anom = 2*np.pi + mean_anom

    k = 0.85
    ecc_anom = mean_anom + np.sign(np.sin(mean_anom))* k * e
    #ecc_anom = mean_anom
    if (e > 0.8):
    ecc_anom = np.pi

    abstol,reltol = 1.0e-8, 1.0e-8
    iter = 0
    while(True):
    f = ecc_anom -e * np.sin(ecc_anom) - mean_anom
    fprime = 1.0 - e * np.cos(ecc_anom)
    fprime2 = e * np.sin(ecc_anom)
    fprime3 = e * np.cos(ecc_anom)
    delta1 = - f / fprime
    delta2 = - f /(fprime + 0.5 * delta1 * fprime2)
    delta3 = - f /(fprime + 0.5 * delta2 * fprime2 + 0.16666666666 * delta2**2 * fprime3) 

    if (delta3 == 0): break

    if (np.abs(ecc_anom) > 0.0):
      abserr,relerr = np.abs(delta3),np.abs(delta3)/np.abs(ecc_anom)
    else:
      abserr,relerr = np.abs(delta3),1.0e40

    ecc_anom+=delta3
    #print iter,ecc_anom,e,delta3

    if (np.abs(ecc_anom) > abstol/reltol):
      if (abserr < abstol): break
    else:
      if (relerr < reltol): break
    iter+=1      

    return ecc_anom % (2*np.pi)

"""##############################################################################"""

def HyperbolicKeplerEquation(mean_anom,e):
    """
    Finds the hyperbolic anomaly 'hyp_anom' 
    given the mean anomaly 'mean_anom' and eccentricity 'e'
    """
    #while (np.abs(mean_anom) > 2.0*np.pi):
    #  mean_anom-=2.0*np.pi*np.sign(mean_anom)
    #if (np.abs(mean_anom) > np.pi): mean_anom-=2.0*np.pi*np.sign(mean_anom)

    if (mean_anom > 0):
    hyp_anom = np.log(mean_anom/np.e + 1.8)
    else:
    hyp_anom = -np.log(np.abs(mean_anom)/np.e + 1.8)
    abstol,reltol = 1.0e-10, 1.0e-10
    iter = 0
    while(True):
    f = e * np.sinh(hyp_anom) - hyp_anom - mean_anom
    fprime =  e * np.cosh(hyp_anom) - 1.0
    fprime2 = e * np.sinh(hyp_anom) 
    fprime3 = e * np.cosh(hyp_anom)
    delta1 = - f / fprime
    delta2 = - f /(fprime + 0.5 * delta1 * fprime2)
    delta3 = - f /(fprime + 0.5 * delta2 * fprime2 + 0.16666666666 * delta2**2 * fprime3) 

    if (delta3 == 0): break

    if (np.abs(hyp_anom) > 0.0):
      abserr,relerr = np.abs(delta3),np.abs(delta3)/np.abs(hyp_anom)
    else:
      abserr,relerr = np.abs(delta3),1.0e40
      
    hyp_anom+=delta3
    #print iter,hyp_anom,e,delta3

    if (np.abs(hyp_anom) > abstol/reltol):
      if (abserr < abstol): break
    else:
      if (relerr < reltol): break
    iter+=1
      
    return hyp_anom


"""##################################################################"""

def ParabolicKeplerEquation(mean_anom): #jerome cardan's method
    """
    Finds the eccentric anomaly 'tan_trueanom2' 
    given the mean anomaly 'mean_anom' (eccentricity 'e' = 1)
    """
    #while (np.abs(mean_anom) > 2.0*np.pi):
    #  mean_anom-=2.0*np.pi*np.sign(mean_anom)
    #if (np.abs(mean_anom) > np.pi): mean_anom-=2.0*np.pi*np.sign(mean_anom)

    B = 1.5 * mean_anom
    tan_trueanom_2 = (B + np.sqrt(1 + B *  B))**(1.0/3) - 1.0/(B + np.sqrt(1 + B *  B))**(1.0/3)
    return tan_trueanom_2

"""##################################################################"""

