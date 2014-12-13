import numpy as np
#from pylab import *   ??????

from structures import *
import avatar

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
    if ((incl > 1.e-3) & (incl < np.pi-1.0e-3)):
        sinw_plus_f = position.get('z') / position.length() / np.sin(orbital.get('i'))
        cosw_plus_f = 1.0 / np.cos(orbital.get('node')) * 
                      ( position.get('x') / position.length() + 
                        np.sin(orbital.get('node')) * sinw_plus_f * np.cos(orbital.get('i') )
                                      
        argw_plus_true_anom = np.arctan2(sinw_plus_f, cosw_plus_f)
    else:
        argw_plus_true_anom = np.cos(orbital.get('i')) * np.arctan2(position.get('y'), position.get('x'))

    true_anom = np.arctan2(orbital.get('a') * (1.0 - orbital.get('ecc') **2) / h.length() / orbital.get('ecc') * r_dot,
                             1.0/ecc*(orbital.get('a') / R * (1.0 - orbital.get('ecc')**2) - 1.0))
    argw = argw_plus_true_anom - true_anom

    true_anom %= 2 * np.pi
    if (true_anom < 0.0): 
        true_anom += 2.0 * np.pi
    
    argw %= 2 * np.pi
    if (argw < 0.0): 
        argw += 2.0 * np.pi

    orbital.set_w(argw)
    
    # (6) Mean Anomaly
    argw = avatar.calculate_mean_anom()
    orbital.set_argw(argw)
    

    return orbital

"""#############################################################################"""

def orbit(peri,e,I,omega,Omega,M,mu):
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
 
    if (e < 1): #Elliptic
        a = peri/(1.0 - e)
        ecc_anom = KeplerEquation(M,e)
        x = a * (np.cos(ecc_anom) - e)
        y = a * np.sqrt(1 - e * e) * np.sin(ecc_anom)
        xdot = -np.sqrt(mu/a) / (1.0 - e*np.cos(ecc_anom))* np.sin(ecc_anom)
        ydot = np.sqrt(mu/a) / (1.0 - e*np.cos(ecc_anom))* np.cos(ecc_anom) * np.sqrt(1.0 - e * e)
        f = np.arctan2(np.sqrt(1.0 - e * e) * np.sin(ecc_anom),np.cos(ecc_anom) - e)
        radius = a*(1.0 - e*np.cos(ecc_anom))
    
    elif (e > 1): #Hyperbolic
        a = peri/(1.0 - e)
        hyp_anom = HyperbolicKeplerEquation(M,e)
        x = a * (np.cosh(hyp_anom) - e)
        y = -a * np.sqrt(e * e - 1.0) * np.sinh(hyp_anom)
        xdot = -np.sqrt(mu/-a) / (e*np.cosh(hyp_anom)- 1.0)* np.sinh(hyp_anom)
        ydot = np.sqrt(mu/-a) / (e*np.cosh(hyp_anom) - 1.0)* np.cosh(hyp_anom) * np.sqrt(e * e - 1.0)
        f = np.arctan2(np.sqrt(e * e - 1.0) * np.sinh(hyp_anom),e - np.cosh(hyp_anom))
        radius = a*(1.0 - e * np.cosh(hyp_anom))
    
    elif (e == 1): #Parabolic
        tan_trueanom_2 = ParabolicKeplerEquation(M)
        x = peri * (1.0 - tan_trueanom_2**2)
        y = 2.0* peri * tan_trueanom_2
        if (peri == 0):
            print peri
        xdot = -np.sqrt(2.0 * mu / peri) / (1.0 + tan_trueanom_2**2) * tan_trueanom_2
        ydot =  np.sqrt(2.0 * mu / peri) / (1.0 + tan_trueanom_2**2)
        f = np.arctan2(2.0 * tan_trueanom_2/(1.0 + tan_trueanom_2**2),(1.0 - tan_trueanom_2**2)/(1.0 + tan_trueanom_2**2))
        radius = peri * (1.0 + tan_trueanom_2**2)
        
    else: # ?????
        return np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
   
    """  
    #fac = np.sqrt(mu/a/(1.0-e*e))
    
    
    #X = radius*(np.cos(Omega)*cos(omega+f) - \
    #                 np.sin(Omega)*sin(omega+f)*np.cos(I))
    #Y = radius*(np.sin(Omega)*cos(omega+f) + \
    #                 np.cos(Omega)*sin(omega+f)*np.cos(I))           
    #Z = radius*np.sin(omega+f)*np.sin(I)
    
    #VX = fac*(- np.sin(omega+f)- e*np.sin(omega))
    
    #VY = fac*np.cos(I)*(np.cos(omega+f) + e*np.cos(omega))
    
    #VZ = fac*np.sin(I)*(-np.cos(omega+f) - e*np.cos(omega))
    """
  
    #rotation matrix
    d11 =  np.cos(omega) * np.cos(Omega) - np.sin(omega) * np.sin(Omega) * np.cos(I)
    d12 =  np.cos(omega) * np.sin(Omega) + np.sin(omega) * np.cos(Omega) * np.cos(I)
    d13 =  np.sin(omega) * np.sin(I)
    d21 = -np.sin(omega) * np.cos(Omega) - np.cos(omega) * np.sin(Omega) * np.cos(I)
    d22 = -np.sin(omega) * np.sin(Omega) + np.cos(omega) * np.cos(Omega) * np.cos(I)
    d23 =  np.cos(omega) * np.sin(I)

    X = d11 * x + d21 * y
    Y = d12 * x + d22 * y
    Z = d13 * x + d23 * y
    VX = d11 * xdot + d21 * ydot
    VY = d12 * xdot + d22 * ydot
    VZ = d13 * xdot + d23 * ydot
        
    return X, Y, Z, VX, VY, VZ


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

