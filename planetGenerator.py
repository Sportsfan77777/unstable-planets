import numpy as np
import mercury
import orbital as orb
import sys
import numpy.random as random
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D

"""
This program places planets on orbits around a binary star system.
There are 1 to 5 sensible inputs.
The inputs are not necessary. There are default values in the code.

(1)
This program takes an input integer 'N' and produces N evenly spaced orbits
around a binary star system (e.g. M = 0, 2pi/N, 2*2pi/N, ..., (N-1)*2pi/N)
at each given semi-major axis, each given eccentricity, AND each given inclination. 

(2)
This program takes an input integer 'periapse'
Periapse = 0; Apoapse = 1

(3)
This program takes an input integer 'num_a' and produces num_a evenly spaced 'a'
around a binary star system. For each 'a', 'N' * 'num_e' * 'num_i' planets are produced
(e.g. a = a_min, a_min + step_size, ..., a_max = a_min + (N-1)*step_size)

(4)
This program takes an input integer 'num_i' and produces num_e evenly spaced 'e'
around a binary star system. For each 'i', 'N' * 'num_a' * 'num_e' planets are produced
(e.g. i = i_min, i_min + step_size, ..., i_max = i_min + (N-1)*step_size)

(5)
This program takes an input integer 'num_e' and produces num_e evenly spaced 'e'
around a binary star system. For each 'e', 'N' * 'num_a' * 'num_i' planets are produced
(e.g. e = e_min, e_min + step_size, ..., e_max = e_min + (N-1)*step_size)

"""

""" CONSTANTS AND UNITS """
G = mercury.G
random.seed(42)
""" CONSTANTS AND UNITS """

""" Imported Functions """
# These functions take the conversions from Cartesian to Orbital Elements
# and applies them to arrays of input at a time instead of just one input at a time
makeVectorOrbit = np.vectorize(orb.orbit)
makeVectorElements = np.vectorize(orb.elements)
""" Imported Functions """


"""Star-Related Parameters"""
periapse = 0 # Location of Companion Star: Periapse = 0, Apoapse = 1
"""Star-Related Parameters"""

# Setup entire thing as a function called by a main method.

# Re-write this as an AMUSE option parser

# Also, install AMUSE on this computer...

"""Planet-Related Parameters"""
N = 8       # Number of Planets (Planet = 'Test Mass' <--- No Interaction)

num_a = 10   # Number of Different Semi-Major Axes
			# There is a default range of semi-major axes (a_min, a_max)
			# 'num_a' evenly spaced semi-major axes are produced
			
num_i = 1   # Number of Different Inclinations
			# There is a default range of inclinations (i_min, i_max)
			# 'num_a' evenly spaced inclinations are produced
						
num_e = 1   # Number of Different Eccentricities
			# There is a default range of eccentricities (e_min, e_max)
			# 'num_a' evenly spaced eccentricities are produced
"""Planet-Related Parameters"""


""" Input Parameters """
args = sys.argv[1:]  # ignore 'planetGenerator.py' as input
num_args = len(args)
			
if num_args >= 1:
	N = int(args[0])        # First Parameter (N)
if num_args >= 2:
	periapse = int(args[2]) # Second (periapse = 0, apoapse = 1)
if num_args >= 3:
	num_a = int(args[3])    # Third (a)
if num_args >= 4:
	num_i = int(args[4])    # Fourth (i)
if num_args >= 5:
	num_e = int(args[5])    # Fifth (e)
""" Input Parameters """


""" Binary Orbital Parameters """
a_b = 1.00  # Semi-Major Axis
e_b = 0.40  # Eccentricity
M_b = 1.00  # Total Mass (M1 + M2)
u_b = 0.30  # Mass Ratio  M2/(M1 + M2) = M2/M_b

omega_b = 0 	  		  # Argument of Pericenter
OMEGA_b = 0 	 		  # Longitude of the Ascending Node
mean_b = periapse * 180.0 # Mean Anomaly  // THERE IS A BUG HERE!!!!
i_b = 0					  # Inclination (By definition, zero)
""" Binary Orbital Parameters """

	
""" Range of Planet Parameters """
M_min = 0
M_max = 2 * np.pi

a_min = 3.1 * a_b  # Note: should be greater than a_b !!!! (P-type)
a_max = 4.0 * a_b  

i_min = 0
i_max = 0 

e_min = 0.0001
e_max = 0.0001
""" Range of Planet Parameters """


""" (1) Below: Initialize Planet Parameters """

""" M """
mean_anomalies = np.ones(N) * 2 * np.pi / N
for i in xrange(N):
	mean_anomalies[i] *= i  # i = 0 to N - 1
	
""" a """
semi_major_axes = np.ones(num_a) * a_min
if num_a > 1:
	a_range = a_max - a_min
	step_size = a_range / (num_a - 1) 
	for i in xrange(num_a):
		semi_major_axes[i] += (i * step_size)
		
""" i """ 
inclinations = np.ones(num_i) * i_min
if num_i > 1:
	i_range = i_max - i_min
	step_size = i_range / (num_i - 1) 
	for i in xrange(num_i):
		inclinations[i] += (i * step_size)
		
""" e """
eccentricities = np.ones(num_e) * e_min
if num_e > 1:
	e_range = e_max - e_min
	step_size = e_range / (num_e - 1) 
	for i in xrange(num_e):
		eccentricities[i] += (i * step_size)
		
""" (1) Above: Initialize Planet Parameters """

""" (2) Below: Convert Orbital Elements to Cartesian Coordinates (Stars) """

# Two Stars
X0, Y0, Z0, VX0, VY0, VZ0 = \
	makeVectorOrbit(a_b * (1-e_b), e_b, i_b * np.pi/180.0, \
                    np.pi + omega_b * np.pi/180.0, OMEGA_b * np.pi/180.0, \
                    mean_b * np.pi/180.0, G * M_b)
X1, Y1, Z1, VX1, VY1, VZ1 = \
	makeVectorOrbit(a_b * (1-e_b), e_b, i_b * np.pi/180.0, \
                    omega_b * np.pi/180.0, OMEGA_b * np.pi/180.0, \
                    mean_b * np.pi/180.0, G * M_b)
                    
# Scale Individual Orbits by Mass Ratio u_b & (1 - u_b)
scale0 = u_b
X0, Y0, Z0, VX0, VY0, VZ0 = \
	X0*scale0, Y0*scale0, Z0*scale0, VX0*scale0, VY0*scale0, VZ0*scale0

scale1 = 1.0 - u_b
X1, Y1, Z1, VX1, VY1, VZ1 = \
	X1*scale1, Y1*scale1, Z1*scale1, VX1*scale1, VY1*scale1, VZ1*scale1

# Check that the barycenter is (still) at the origin
baryX = ((1 - u_b) * X0 + u_b * X1)
baryY = ((1 - u_b) * Y0 + u_b * Y1)
baryZ = ((1 - u_b) * Z0 + u_b * Z1)

baryVX = ((1 - u_b) * VX0 + u_b * VX1)
baryVY = ((1 - u_b) * VY0 + u_b * VY1)
baryVZ = ((1 - u_b) * VZ0 + u_b * VZ1)

print "Binary Barycenter:",baryX,baryY,baryZ
print "                  ",baryVX,baryVY,baryVZ

#Also check that the relative motion satisfies the given orbital parameters
deltaX = X1 - X0
deltaY = Y1 - Y0
deltaZ = Z1 - Z0
deltaVX = VX1 - VX0
deltaVY = VY1 - VY0
deltaVZ = VZ1 - VZ0
peri_b, ecc_b, inc_b, arg_b, node_b, true_b, mean_b, eAnom_b = \
	makeVectorElements(deltaX, deltaY, deltaZ, deltaVX, deltaVY, deltaVZ, G * M_b)
print "Central Binary Orbital Elements:", peri_b / (1 - ecc_b), ecc_b, inc_b

# Convert to coordinates with respect to the 1st body
# Big Bodies
x0, y0, z0, vx0, vy0, vz0 = X0-X0, Y0-Y0, Z0-Z0, VX0-VX0, VY0-VY0, VZ0-VZ0
x1, y1, z1, vx1, vy1, vz1 = X1-X0, Y1-Y0, Z1-Z0, VX1-VX0, VY1-VY0, VZ1-VZ0

# Data: Big Bodies (ignore the first one)
# If two (2) big bodies, list has one element (for 2nd big body)
x_b  = np.array([x1])
y_b  = np.array([y1])
z_b  = np.array([z1])
vx_b = np.array([vx1])
vy_b = np.array([vy1])
vz_b = np.array([vz1])
mass_b = np.array([u_b * M_b])
name_b = np.array(['S2']) # Central Object = S1

""" (2) Above: Convert Orbital Elements to Cartesian Coordinates (Stars) """

""" (3) Below: Set Up Arrays (Planets) """

total_planets = N * num_a * num_i * num_e

gather = [] # Cartesian Product of All Orbital Element Combinations
for i in itertools.product(mean_anomalies, semi_major_axes, \
							inclinations, eccentricities):
	gather.append(i)
	
# Mercury Input Arrays
mean_array = [m for (m,s,i,e) in gather]
sm_axis_array = [s for (m,s,i,e) in gather]
inc_array = [i  for (m,s,i,e) in gather]
ecc_array = [e for (m,s,i,e) in gather]
omega_array = np.ones(total_planets) * omega_b
OMEGA_array = np.ones(total_planets) * OMEGA_b


# For Naming the Planets (by m, s, i, e) (or really 100*m, 100*s, 100*i, 100*e)
m_deg_array = [int(round(x * 180.0 / np.pi)) for x in mean_array]
i_deg_array = [int(round(x * 180.0 / np.pi)) for x in inc_array]

planet_names = []
for i in xrange(total_planets):
	a = m_deg_array[i] 
	b = sm_axis_array[i] * 10
	c = i_deg_array[i] * 100 # to show two (2) decimal places
	d = ecc_array[i] * 100
	#planet_names.append(' Planet_M%1.0f_S%1.0f_I%1.0f_E%1.0f'%(a, b, c, d))
	planet_names.append(' M%1.0f_S%1.0f'%(a, b))
	#planet_names.append(' M%d'%(i))
	
#print planet_names

    
# Convert to Cartesian Coordinates
Xp, Yp, Zp, VXp, VYp, VZp = \
	makeVectorOrbit(np.multiply(sm_axis_array, (np.ones(total_planets) - ecc_array)), \
					ecc_array, inc_array, \
					omega_array, OMEGA_array, \
					mean_array, G * M_b)

# Convert to coordinates with respect to the 1st body
x_p, y_p, z_p, vx_p, vy_p, vz_p = Xp-X0, Yp-Y0, Zp-Z0, VXp-VX0, VYp-VY0, VZp-VZ0

print "error below?"

# Convert back to asteroidal coordinates
peri_p, ecc_p, inc_p, arg_p, node_p, true_p, mean_p, eAnom_p = \
	makeVectorElements(x_p, y_p, z_p, vx_p, vy_p, vz_p, G * M_b)
	
print "error above?"
	
a_p = peri_p / (1 - ecc_p) # semi-major axis

""" (3) Above: Set Up Arrays (Planets) """

""" (4) Write Mercury Files """

mercury.write_mercury_files(x_b, y_b, z_b, vx_b, vy_b, vz_b, mass_b,
                            np.zeros(2),np.zeros(2),np.zeros(2),
                            name_b, 'big.in')

spin_x = np.zeros(total_planets)
spin_y = np.zeros(total_planets)
spin_z = np.zeros(total_planets)
mass_p = np.zeros(total_planets)
"""
# Asteroidal Option
mercury.write_mercury_files(a_p, ecc_p, inc_p * 180.0/np.pi, \
							arg_p, node_p, mean_p , mass_p, \
                            spin_x, spin_y, spin_z, \
                            planet_names, 'small.in', \
                            body_type = 'small', coord_type = 'ast')

# Cartesian Option
"""
mercury.write_mercury_files(x_p, y_p, z_p, vx_p, vy_p, vz_p, mass_p, \
                            spin_x, spin_y, spin_z, \
                            planet_names, 'small.in', \
							body_type = 'small', coord_type = 'cart')

""" (4) Write Mercury Files """

""" (5) Plot Orbits (Optional) """

# Options
plot_orbits = 1
plot_barycentric = 1

print "before plot"

if plot_orbits:
    if plot_barycentric:
        plt.plot([0]+X0,[0]+Y0,'bo',mew=0,ms=7.0)
        plt.plot(x_b+X0,y_b+Y0,'bo',mew=0,ms=7.0)
        plt.plot(x_p+X0,y_p+Y0,'go',mew=0,ms=4.0)
    else:
        plt.plot([0],[0],'bo',mew=0,ms=7.0)
        plt.plot(x_b,y_b,'bo',mew=0,ms=7.0)
        plt.plot(x_p,y_p,'go',mew=0,ms=4.0)
    for k in range(0,x_p.shape[0]):
        #anomalies=np.arange(MeanAnom[k],MeanAnom[k]+2*np.pi,0.01)
        #xorb, yorb, zorb, vxorb, vyorb, vorb = vecOrbit(ap[k]*(1-ep[k]),ep[k],Ip[k]*np.pi/180.0,omegap[k],
        #                                                Omegap[k],anomalies,G*MB*(1-muB))
        #xorb, yorb, zorb, vxorb, vyorb, vorb = vecOrbit(a[k]*(1-e[k]),e[k],I[k]*np.pi/180.0,omega[k]*np.pi/180.0,
        #                                                Omega[k]*np.pi/180.0,anomalies,G*MB)
        
        if plot_barycentric:
            peridistp,ep,Ip,omegap,Omegap,TrueAnomp, MeanAnomp,EccAnomp = \
            		makeVectorElements(x_p[k] + X0, y_p[k] + Y0, z_p[k] + Z0, \
            						vx_p[k] + VX0, vy_p[k] + VY0, vz_p[k] + VZ0, G * M_b)
        else:
            peridistp,ep,Ip,omegap,Omegap,TrueAnomp, MeanAnomp,EccAnomp = \
            		makeVectorElements(x_p[k], y_p[k], z_p[k], \
            						vx_p[k], vy_p[k], vz_p[k], G * M_b)
                                        
        anomalies = np.arange(MeanAnomp,MeanAnomp+2*np.pi,0.01)
        
        xorb, yorb, zorb, vxorb, vyorb, vorb = \
        			makeVectorOrbit(peridistp, ep, Ip, omegap, Omegap, anomalies, G * M_b)
        
        plt.plot(xorb,yorb,'k-')

    plt.axis([-5,5,-5,5])
    plt.show()

    if plot_barycentric:
        plt.plot([0]+Y0,[0]+Z0,'bo',mew=0,ms=7.0)
        plt.plot(y_b+Y0,z_b+Z0,'bo',mew=0,ms=7.0)
        plt.plot(y_p+Y0,z_p+Z0,'go',mew=0,ms=4.0)
    else:
        plt.plot([0],[0],'bo',mew=0,ms=7.0)
        plt.plot(y_b,z_b,'bo',mew=0,ms=7.0)
        plt.plot(y_p,z_p,'go',mew=0,ms=4.0)
    for k in range(0,y_p.shape[0]):
        plt.plot(yorb,zorb,'k-')

    plt.axis([-5,5,-5,5])
    plt.show()

    if plot_barycentric:
        fig=plt.figure()
        ax=Axes3D(fig)
        ax.plot([0]+X0,[0]+Y0,[0]+Z0,'bo')
        ax.plot(x_b+X0,y_b+Y0,z_b+Z0,'bo')
        ax.plot(x_p+X0,y_p+Y0,z_p+Z0,'go')
        for k in range(0,y_p.shape[0]):
            ax.plot(xorb,yorb,zorb,'k-')
        ax.axis('equal')
        ax.auto_scale_xyz([-1,1],[-1,1],[-1,1])
        plt.show()

""" (5) Plot Orbits (Optional) """
