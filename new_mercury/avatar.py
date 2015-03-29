"""
masters (calculates) all 4 (really 6) elements
"""

import numpy as np
import kepler
from structures import *

def calculate_a(position, velocity, mu):
    """ 
    1st, a
    """
    return 1.0 / (2.0 / position.length() - velocity.norm_sq() / mu)

def calculate_sma(position, velocity, mu):
    return calculate_a(position, velocity, mu)

def calculate_e(position, velocity, mu, angular_momentum):
    """ 
    2nd, e
    """
    h = angular_momentum
    ecc_sq = 1.0 + h.norm_sq() / mu * (velocity.norm_sq() / mu - 2.0 / position.length())
    if (ecc_sq < 0): 
        return 0
    else: 
        return np.sqrt(ecc_sq)

def calculate_ecc(position, velocity, mu, angular_momentum):
    return calculate_e(position, velocity, mu, angular_momentum)

def calculate_i(angular_momentum):
    """ 
    3rd, i 
    """
    h = angular_momentum
    return np.arccos(h.get('hz') / h.length())

def calculate_inc(angular_momentum):
    return calculate_i(angular_momentum)

def calculate_W(elements, angular_momentum):
    """ 
    4th, W (node)
    Requires elements: 'inc'
    """
    h = angular_momentum
    if (elements.get('i') == 0.0):
        return 0.0
    else:
        if (h.get('hz') > 0):
            negate = -1
        else:
            negate = 1
        return np.arctan2(h.get('hx') / h.length() / np.sin(elements.get('i')), 
                          negate * h.get('hy') / h.length() / np.sin(elements.get('i')))

def calculate_node(elements, angular_momentum):
    return calculate_W(elements, angular_momentum)

def calculate_w(elements, position, velocity, angular_momentum, positive = True):
    """ 
    5th, w (argument of pericenter)
    Requires elements: 'a', 'inc', 'node'
    Optional parameter: 'positive' ensures 0 < w < 2*pi
    """
    
    if ((elements.get('i') > 1.e-3) and (elements.get('i') < np.pi-1.0e-3)):
        sinw_plus_f = position.get('z') / position.length() / np.sin(elements.get('i'))
        cosw_plus_f = 1.0 / np.cos(elements.get('node')) * ( position.get('x') / position.length() + np.sin(elements.get('node')) * sinw_plus_f * np.cos(elements.get('i')) )
        argw_plus_true_anom = np.arctan2(sinw_plus_f, cosw_plus_f)
    else:
        argw_plus_true_anom = np.cos(elements.get('i')) * np.arctan2(position.get('y'), position.get('x'))

    true_anom = calculate_true_anom(elements, position, velocity, angular_momentum)
    argw = argw_plus_true_anom - true_anom

    argw %= 2 * np.pi
    if positive and argw < 0.0:
        argw += 2.0 * np.pi

    return argw

def calculate_argw(elements, position, velocity, angular_momentum):
    return calculate_w(elements, position, velocity, angular_momentum)

def calculate_M(elements, position, velocity, angular_momentum):
    """ 
    6th, M (mean anomaly)
    Requires elements: 'a', 'ecc'
    """
    
    true_anom = calculate_true_anom(elements, position, velocity, angular_momentum)

    if (elements.get('ecc') < 1.0): # Elliptic
        tan_ecc_anom_2 = np.tan(0.5 * true_anom) * np.sqrt((1.0 - elements.get('ecc')) / (1.0 + elements.get('ecc')))
        tan_ecc_anom = 2.0 * tan_ecc_anom_2 / (1.0 - tan_ecc_anom_2**2)
        cos_ecc_anom = (1.0 - position.length() / elements.get('a') ) / elements.get('ecc')

        ecc_anom = np.arctan2(tan_ecc_anom * cos_ecc_anom, cos_ecc_anom)

        if (ecc_anom < 0): 
            ecc_anom = 2.0 * np.pi + ecc_anom
        
        mn_anom = ecc_anom - elements.get('ecc') * np.sin(ecc_anom)
        return mn_anom
    else: # Hyperbolic
        tanhhyp_anom_2 = np.tan(0.5 * true_anom) * np.sqrt((elements.get('ecc') - 1.0) / (1.0 + elements.get('ecc')))
        tanhhyp_anom = 2.0 * tanhhyp_anom_2 / (1.0 + tanhhyp_anom_2**2)
        hyp_anom = np.arctanh(tanhhyp_anom)
            
        mn_anom = 0.0 - hyp_anom + elements.get('ecc') * np.sinh(hyp_anom)
        return mn_anom

def calculate_mean_anom():
    return calculate_M(elements, position, velocity, angular_momentum)

####### SUPPLEMENTAL #######

def calculate_true_anom(elements, position, velocity, angular_momentum, positive = False):
    h = angular_momentum

    r_dot_v = position.dot(velocity)
    if (r_dot_v > 0 ):
        negate = 1
    else:
        negate = -1
    r_dot = negate * np.sqrt(abs(velocity.norm_sq() - h.norm_sq() / position.norm_sq()))

    true_anom = np.arctan2(elements.get('a') * (1.0 - elements.get('ecc') **2) / h.length() / elements.get('ecc') * r_dot,
                             1.0 / elements.get('ecc') * (elements.get('a') / position.length() * (1.0 - elements.get('ecc')**2) - 1.0))

    true_anom %= 2 * np.pi
    if positive and true_anom < 0.0: 
        true_anom += 2.0 * np.pi

    return true_anom

####### Kepler's Equation ########

def calculate_geometric_anomaly(elements):
    """
    Requires elements: 'e', 'M'
    if 'e' < 1, calculate eccentric anomaly
    if 'e' = 1, 
    if 'e' > 1, calculate Hyperbolic anomaly
    """

    if elements.get('ecc') < 1:
        kepler.calculate_eccentric_anomaly(elements)
    elif elements.get('ecc') == 1:
        kepler.calculate_parabolic_anomaly(elements)
    elif elements.get('ecc') > 1:
        kepler.calculate_hyperbolic_anomaly(elements)


