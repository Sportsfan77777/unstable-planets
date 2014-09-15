import numpy as np
import mercury
import orbital as orb
import sys
import numpy.random as random
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D

from optparse import OptionParser

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


def grid_planets(num_M = 8, num_a = 10, num_i = 1, num_e = 1,
                 min_M = 0, max_M = 360, min_sma = 3.1, max_sma = 4.0,
                 min_ecc = 0, max_ecc = 0, min_inc = 0, max_inc = 0,
                 mass_bin = 1.0, u_bin = 0.5, a_bin = 1.0, e_bin = 1.0,
                 argw_bin = 0.0, node_bin = 0.0, mean_anom_bin = 1.0, i_bin = 0,
                 plot_orbits = 0, plot_barycentric = 0):
    """
    Explain parameters here
    
    ####*** Star-Related Parameters ***####
    
    periapse = 0 # Location of Companion Star: Periapse = 0, Apoapse = 1
    
    ####*** Planet-Related Parameters ***####
    
    num_M = 8       # Number of Planets (Planet = 'Test Mass' <--- No Interaction)

	num_a = 10   # Number of Different Semi-Major Axes
				# There is a default range of semi-major axes (a_min, a_max)
				# 'num_a' evenly spaced semi-major axes are produced
			
	num_i = 1   # Number of Different Inclinations
				# There is a default range of inclinations (i_min, i_max)
				# 'num_a' evenly spaced inclinations are produced
						
	num_e = 1   # Number of Different Eccentricities
				# There is a default range of eccentricities (e_min, e_max)
				# 'num_a' evenly spaced eccentricities are produced
				
    ####*** Binary Orbital Parameters ***####
    
	a_b = 1.00  # Semi-Major Axis
	e_b = 0.40  # Eccentricity
	M_b = 1.00  # Total Mass (M1 + M2)
	u_b = 0.30  # Mass Ratio  M2/(M1 + M2) = M2/M_b

	omega_b = 0 	  		  # Argument of Pericenter
	OMEGA_b = 0 	 		  # Longitude of the Ascending Node
	mean_b = periapse * 180.0 # Mean Anomaly  // THERE IS A BUG HERE!!!!
	i_b = 0					  # Inclination (By definition, zero)
	
	####*** Range of Planet Parameters ***####
	
	M_min = 0
	M_max = 2 * np.pi

	a_min = 3.1 * a_b  # Note: should be greater than a_b !!!! (P-type)
	a_max = 4.0 * a_b  

	i_min = 0
	i_max = 0 

	e_min = 0.0001
	e_max = 0.0001
	
    """

    line_break = 0

    """ (1) Below: Initialize Planet Parameters """

    # M = mean anomaly
    mean_anomalies = np.ones(num_M) * 2 * np.pi / num_M
    for i in xrange(num_M):
    	mean_anomalies[i] *= i  # i = 0 to num_M - 1
	
	# a = semimajor axis
    semi_major_axes = np.ones(num_a) * min_sma
    if num_a > 1:
    	sma_range = max_sma - min_sma
    	step_size = sma_range / (num_a - 1) 
    	for i in xrange(num_a):
    		semi_major_axes[i] += (i * step_size)
		
    # i = inclination
    inclinations = np.ones(num_i) * min_inc
    if num_i > 1:
    	inc_range = max_inc - min_inc
    	step_size = inc_range / (num_i - 1) 
    	for i in xrange(num_i):
    		inclinations[i] += (i * step_size)
    	
    # e = eccentricity
    eccentricities = np.ones(num_e) * min_ecc
    if num_e > 1:
    	ecc_range = max_ecc - min_ecc
    	step_size = ecc_range / (num_e - 1) 
    	for i in xrange(num_e):
    		eccentricities[i] += (i * step_size)
    	
    """ (1) Above: Initialize Planet Parameters """

    """ (2) Below: Convert Orbital Elements to Cartesian Coordinates (Stars) """

    # Two Stars
    X0, Y0, Z0, VX0, VY0, VZ0 = \
    	makeVectorOrbit(a_bin * (1-e_bin), e_bin, i_bin * np.pi/180.0, \
    					np.pi + argw_bin * np.pi/180.0, node_bin * np.pi/180.0, \
    					mean_anom_bin * np.pi/180.0, G * mass_bin)
    X1, Y1, Z1, VX1, VY1, VZ1 = \
    	makeVectorOrbit(a_bin * (1-e_bin), e_bin, i_bin * np.pi/180.0, \
    					argw_bin * np.pi/180.0, node_bin * np.pi/180.0, \
    					mean_anom_bin * np.pi/180.0, G * mass_bin)
    				
    # Scale Individual Orbits by Mass Ratio u_bin & (1 - u_bin)
    scale0 = u_bin
    X0, Y0, Z0, VX0, VY0, VZ0 = \
    	X0*scale0, Y0*scale0, Z0*scale0, VX0*scale0, VY0*scale0, VZ0*scale0

    scale1 = 1.0 - u_bin
    X1, Y1, Z1, VX1, VY1, VZ1 = \
    	X1*scale1, Y1*scale1, Z1*scale1, VX1*scale1, VY1*scale1, VZ1*scale1

    # Check that the barycenter is (still) at the origin
    baryX = ((1 - u_bin) * X0 + u_bin * X1)
    baryY = ((1 - u_bin) * Y0 + u_bin * Y1)
    baryZ = ((1 - u_bin) * Z0 + u_bin * Z1)

    baryVX = ((1 - u_bin) * VX0 + u_bin * VX1)
    baryVY = ((1 - u_bin) * VY0 + u_bin * VY1)
    baryVZ = ((1 - u_bin) * VZ0 + u_bin * VZ1)

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
    	makeVectorElements(deltaX, deltaY, deltaZ, deltaVX, deltaVY, deltaVZ, G * mass_bin)
    print "Central Binary Orbital Elements:", peri_b / (1 - ecc_b), u_bin, ecc_b, inc_b

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
    mass_b = np.array([u_bin * mass_bin])
    name_b = np.array(['S2']) # Central Object = S1

    """ (2) Above: Convert Orbital Elements to Cartesian Coordinates (Stars) """

    """ (3) Below: Set Up Arrays (Planets) """

    total_planets = num_M * num_a * num_i * num_e

    gather = [] # Cartesian Product of All Orbital Element Combinations
    for i in itertools.product(mean_anomalies, semi_major_axes, \
    							inclinations, eccentricities):
    	gather.append(i)

    # Mercury Input Arrays
    mean_array = [m for (m,s,i,e) in gather]
    sm_axis_array = [s for (m,s,i,e) in gather]
    inc_array = [i  for (m,s,i,e) in gather]
    ecc_array = [e for (m,s,i,e) in gather]
    omega_array = np.ones(total_planets) * argw_bin
    OMEGA_array = np.ones(total_planets) * node_bin


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
    	planet_names.append(' M%1.0f_S%1.0f'%(a, b)) # (sm-axis, mean anomaly)
    	#planet_names.append(' M%d'%(i))
	
    #print planet_names
    
    # Convert to Cartesian Coordinates
    Xp, Yp, Zp, VXp, VYp, VZp = \
    	makeVectorOrbit(np.multiply(sm_axis_array, (np.ones(total_planets) - ecc_array)), \
    					ecc_array, inc_array, \
    					omega_array, OMEGA_array, \
						mean_array, G * mass_bin)

    # Convert to coordinates with respect to the 1st body
    x_p, y_p, z_p, vx_p, vy_p, vz_p = Xp-X0, Yp-Y0, Zp-Z0, VXp-VX0, VYp-VY0, VZp-VZ0

    print "error below?"

    # Convert back to asteroidal coordinates
    peri_p, ecc_p, inc_p, arg_p, node_p, true_p, mean_p, eAnom_p = \
    	makeVectorElements(x_p, y_p, z_p, vx_p, vy_p, vz_p, G * mass_bin)
	
    print "error above?"
	
    a_p = peri_p / (1 - ecc_p) # semi-major axis

    """ (3) Above: Set Up Arrays (Planets) """
    
    integration_dir = "simtest"

    """ (4) Write Mercury Files """

    mercury.write_mercury_files(x_b, y_b, z_b, vx_b, vy_b, vz_b, mass_b,
    							np.zeros(2),np.zeros(2),np.zeros(2),
    							name_b, 'big.in',
    							directory = integration_dir)

    spin_x = np.zeros(total_planets)
    spin_y = np.zeros(total_planets)
    spin_z = np.zeros(total_planets)
    mass_p = np.zeros(total_planets)
    
    
    central_mass = mass_bin - mass_b
    mercury.write_param_file(central_mass, dest_dir = integration_dir)
    
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
    							planet_names, 'small.in', directory = integration_dir, \
    							body_type = 'small', coord_type = 'cart')

    """ (4) Write Mercury Files """

    """ (5) Plot Orbits (Optional) """

    # Options
    #plot_orbits = 0
    #plot_barycentric = 0

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
    									vx_p[k] + VX0, vy_p[k] + VY0, vz_p[k] + VZ0, G * mass_bin)
    		else:
    			peridistp,ep,Ip,omegap,Omegap,TrueAnomp, MeanAnomp,EccAnomp = \
    					makeVectorElements(x_p[k], y_p[k], z_p[k], \
    									vx_p[k], vy_p[k], vz_p[k], G * mass_bin)
    									
    		anomalies = np.arange(MeanAnomp,MeanAnomp+2*np.pi,0.01)
    	
    		xorb, yorb, zorb, vxorb, vyorb, vorb = \
    					makeVectorOrbit(peridistp, ep, Ip, omegap, Omegap, anomalies, G * mass_bin)
    	
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
	
def new_option_parser():
  """ parameters for simulation of 1 to N independent planets around two binary stars """
  result = OptionParser()
  result.add_option("--dir", 
                    dest="dir", default="sim",
                    help="directory base name [%default]")
  result.add_option("--num_M", 
                    dest="num_M", type="int", default = 8,
                    help="number of different planet mean anomalies [%default]")
  result.add_option("--num_a", 
                    dest="num_a", type="int", default = 10,
                    help="number of different planet semi-major axes [%default]")
  result.add_option("--num_i", 
                    dest="num_i", type="int", default = 1,
                    help="number of different planet inclinations [%default]")
  result.add_option("--num_e", 
                    dest="num_e", type="int", default = 1,
                    help="number of different planet eccentricities [%default]")
  result.add_option("--min_M", 
                    dest="min_M", type="int", default = 0,
                    help="planet mean anomaly minimum [%default]")
  result.add_option("--max_M", 
                    dest="max_M", type="int", default = 360,
                    help="planet mean anomaly maximum [%default]")
  result.add_option("--min_sma", 
                    dest="min_sma", type="float", default = 3.1,
                    help="planet semi-major axis minimum [%default]")
  result.add_option("--max_sma", 
                    dest="max_sma", type="float", default = 4.0,
                    help="planet semi-major axis maximum [%default]")
  result.add_option("--min_ecc", 
                    dest="min_ecc", type="float", default = 0.0000001,
                    help="planet eccentricity minimum [%default]")
  result.add_option("--max_ecc", 
                    dest="max_ecc", type="float", default = 0.0000001,
                    help="planet eccentricity maximum [%default]")
  result.add_option("--min_inc", 
                    dest="min_inc", type="float", default = 0.0,
                    help="planet inclination minimum in deg [%default]")
  result.add_option("--max_inc", 
                    dest="max_inc", type="float", default = 180.0,
                    help="planet inclination maximum in deg [%default]")
  result.add_option("--mass_bin",
                    dest="mass_bin", type="float", default=1.0,
                    help="total mass of two binary stars [%default]")
  result.add_option("--u_bin",
                    dest="u_bin", type="float", default = 0.3,
                    help="binary star mass ratio M2/mass_b [%default]")
  result.add_option("--e_bin",
                    dest="e_bin", type="float", default = 0.4,
                    help="eccentricity of the binary stars [%default]")
  result.add_option("--a_bin",
                    dest="a_bin", type="float", default = 1.0,
                    help="semi-major axis of binary [%default]")
  result.add_option("--argw_bin",
                    dest="argw_bin", type="float", default = 0.0,
                    help="argument of pericenter of the binary stars in deg [%default]")
  result.add_option("--node_bin",
                    dest="node_bin", type="float", default = 0.0,
                    help="longitude of the ascending node of the binary stars in deg [%default]")
  result.add_option("--mean_anom_bin",
                    dest="mean_anom_bin", type="float", default = 0.0,
                    help="mean anomaly of the binary stars in deg, \
                          periapse = 0.0, apoapse = 180.0 [%default]")
  result.add_option("--i_bin",
                    dest="i_bin", type="float", default = 0.0,
                    help="inclination of the binary stars in deg [%default]")
  result.add_option("--plot_orb",
                    dest="plot_orb", type="int", default = 1,
                    help="plot orbits (0 or 1) [%default]")
  result.add_option("--plot_bary",
                    dest="plot_bary", type="int", default = 1,
                    help="plot barycentric orbits (0 or 1) (not totally sure what this does) [%default]")
                    
  return result
                    
                    
def execute_main(o, arguments):
   """ main method for command line call or pseudo_main call """
   
   print "Executing"
   grid_planets(num_M = o.num_M, num_a = o.num_a, num_i = o.num_i, num_e = o.num_e,
                min_M = o.min_M, max_M = o.max_M, min_sma = o.min_sma, max_sma = o.max_sma,
                min_ecc = o.min_ecc, max_ecc = o.max_ecc, 
                min_inc = o.min_inc, max_inc = o.max_inc,
                mass_bin = o.mass_bin, u_bin = o.u_bin, a_bin = o.a_bin, e_bin = o.e_bin,
                argw_bin = o.argw_bin, node_bin = o.node_bin, 
                mean_anom_bin = o.mean_anom_bin, i_bin = o.i_bin,
                plot_orbits = o.plot_orb, plot_barycentric = o.plot_bary)
                
   print "Done"

def psuedo_main(args):
   """ option parser input not from the command line """
   o, arguments  = new_option_parser().parse_args(args)
   execute_main(o, arguments)
                    
if __name__ in ('__main__'):
   o, arguments  = new_option_parser().parse_args()
   execute_main(o, arguments)
   
   
