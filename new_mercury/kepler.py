"""
Kepler Equations for all geometries
"""

from structures import *

def newton_method(function, function_prime, guess = 1.234, num_iter = 8):
	"""
	finds the root of a function given its derivative function_prime, a guess, and a number of iterations
	"""
	xi = guess
	for i in range(num_iter):
		f_xi = function(xi)
		fp_xi = function_prime(xi)

		xi -= f_xi / fp_xi

	return x

def calculate_eccentric_anomaly(elements):
	"""
	Requires elements: 'e', 'M'
	"""

	def test_function(elements, ecc_anom):
		return ecc_anom - elements.get('ecc') * np.sin(ecc_anom)

	def test_function_prime(elements, ecc_anom):
		return 1.0 - elements.get('ecc') * np.cos(ecc_anom)

	function = lambda x : test_function(elements, x)
	function_prime = lambda x : test_function_prime(elements, x)

	ecc_anom = newton_method(function, function_prime)
	return ecc_anom

def calculate_hyperbolic_anomaly(elements):
	"""
	Requires elements: 'e', 'M'
	"""

	def test_function(elements, ecc_anom):
		return ecc_anom - elements.get('ecc') * np.sin(ecc_anom)

	def test_function_prime(elements, ecc_anom):
		return 1.0 - elements.get('ecc') * np.cos(ecc_anom)

	function = lambda x : test_function(elements, x)  ## ### FIX THIS  !!!!!
	function_prime = lambda x : test_function_prime(elements, x) ##### FIX THIS !!!!!!

	hyp_anom = newton_method(function, function_prime)
	return hyp_anom


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

def HyperbolicKeplerEquation(mean_anom,e):
    """
    Finds the hyperbolic anomaly 'hyp_anom' 
    given the mean anomaly 'mean_anom' and eccentricity 'e'
    """
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

