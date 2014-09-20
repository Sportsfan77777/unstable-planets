import math
import matplotlib.pyplot as mpy
from functools import partial

def newtonMethod(f, df, guess = 3.0, iter = 1):
	""" 
	Newton's Method for iteratively solving an equation
	Calculates the root of 'f' using its derivative 'df', and an initial guess 'guess'
	Returns the root of 'f' and 'iter', the number of iterations
	"""
	y = f(guess) / df(guess)
	next = guess - y  # next guess
	#print next
	if (abs(y) < pow(10,-8)):
		return next, iter
	else:
		return newtonMethod(f, df, guess = next, iter = iter + 1)
		
""""""

def n_f(mp, ms, a):
	""" n = 2pi/P, the mean motion """
	z = 4 * (math.pi)**2 * (1 + mp/ms)
	return math.sqrt(z)
		
def t_f(M,T,n):
    """ t = time """
    return M/n + T
    
def M_f(t,T,n):
	""" M = mean anomaly"""
	return n*(t-T)
	
def r_f(a,e,E):
	""" r = radius """
	return a*(1 - e*math.cos(E))
	
def v_f(e,E):
	""" v = true anomaly """
	return 2 * math.atan(math.sqrt((1+e)/(1-e))*math.tan(E/2))
	
""""""

def k(e,M,E):
    """Kepler's Equation: solve for E"""
    return E - e*math.sin(E) - M 
    
def dk(e,M,E):
	"""Derivative of Kepler's Equation: solve for E"""
	return 1 - e*math.cos(E)
	
""""""

def advanceOrbit(e,M):
	"""Solve for E"""
	k_p = partial(k, e, M)
	dk_p = partial(dk, e, M)
	
	next_E, iter = newtonMethod(k_p,dk_p)
	if (next_E >= math.pi):
		next_E = next_E - 2*math.pi 
	return next_E
	
def determineOrbit(a, e, w, ml, mp = 0.00001, ms = 1, T = 0, steps = 20000, yrs = 1.0):
    """
    Plot orbit from orbital elements
    
    a = semi-major axis
    e = eccentricity
    w = longitude of peri-center (in radians)
    ml = mean longitude at a reference epoch (in radians)
    
    mp = mass of planet (test mass)
    ms = mass of star # the way it is currently set up, this must be 1
    
    T = time of peri-center passage
    
    steps = number of time-steps
    yrs = number of years # must be a float
    """
    d  = steps / yrs
    
    ec = e
    sm_axis = a
    
    t = [j/d for j in xrange(steps)] # Initialize arrays
    
    r = [j for j in xrange(steps)]
    v = [j for j in xrange(steps)]
    E = [j for j in xrange(steps)]
    
    x = [j for j in xrange(steps)]
    y = [j for j in xrange(steps)]
    
    vx = [j for j in xrange(steps)]
    vy = [j for j in xrange(steps)]
    
    for i in xrange(len(t)):
    	#a = sm_axis # (1) Update a
    	#e = ec # (2) Update e
    	
    	n = n_f(mp, ms, a) # (3) Find n
    	M = M_f(t[i], T, n) # (4) Find M
    	
    	E[i] = advanceOrbit(e,M) # (5) Calculate new E
    	r[i] = r_f(a,e,E[i]) # (6) Calculate new r
    	v[i] = v_f(e,E[i]) # (7) Calculate true anomaly
    	
    	x_tmp = a*e + r[i] * math.cos(v[i])
    	y_tmp = r[i] * math.sin(v[i])
    	
    	x[i] = math.cos(w) * x_tmp - math.sin(w) * y_tmp # rotate by w
    	y[i] = math.sin(w) * x_tmp + math.cos(w) * y_tmp
    	    	    	    	    	
    	if (i != 0):
    		vx[i] = (x[i] - x[i-1])/(t[i] - t[i-1])
    		vy[i] = (y[i] - y[i-1])/(t[i] - t[i-1])
	    	     	
    return t, r, v, E, x, y, vx, vy
    	
tpp = 0
t, r, v, E, x, y, vx, vy = determineOrbit(1.0, 0, 0, 0, T = tpp)
t2, r2, v2, E2, x2, y2, vx2, vy2 = determineOrbit(1.0, 0.5, 0, 0, T = tpp)
t3, r3, v3, E3, x3, y3, vx3, vy3 = determineOrbit(1.0, 0.98, 0, 0, T = tpp)

"""
print r
print v
print E

print x
print y

print vx
print vy
"""

"""
mpy.title('r vs t')
mpy.plot(t, r)
mpy.show()

mpy.title('v vs t')
mpy.plot(t, v)
mpy.show()

mpy.title('x vs t')
mpy.plot(t, x)
mpy.show()

mpy.title('y vs t')
mpy.plot(t, y)
mpy.show()
"""

mpy.title('y vs x')
mpy.plot(x, y)
mpy.plot(x2, y2)
mpy.plot(x3, y3)
mpy.show()

"""
mpy.title('vx vs t')
mpy.plot(t, vx)
mpy.show()

mpy.title('vy vs t')
mpy.plot(t, vy)
mpy.show()
"""

mpy.title('E vs v')
mpy.plot(v, E, color = 'b')
mpy.plot(v2, E2, color = 'g')
mpy.plot(v3, E3, color = 'r')
mpy.show()