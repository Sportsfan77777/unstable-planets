"""
Integrator of Circumbinary Planetary Systems (w/ One Planet)
"""

from amuse.units.optparse import OptionParser
from amuse.units import units
from amuse.units import constants
from amuse.units import quantities
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSuperset
from amuse.community.huayno.interface import Huayno
#from kepler6.interface import Kepler
from amuse.io import write_set_to_file

# Select 'nbody' file
import circumbinary_nbody_01 as nbody

from orbital_elements import Orb_Kepler

from misc import mkdir

def run_simulation(planet_elements, binary_elements, star_masses,
                   huayno_eta, total_sim_time, n_steps,
                   directory, outfile, file_redir = None):
  """ 
  runs simulation out of separate 'nbody' file 
  (which must have the methods (i) initialize and (ii) integrate)
  """
  
  system = nbody.initialize_circumbinary_system(
                           binary_elements, planet_elements, star_masses)
                           
  # n_steps? file_redir?
  nbody.integrate_circumbinary_system(system, total_sim_time, n_steps,
                         directory, outfile, file_redir, huayno_eta)
  

def new_option_parser():
  """ parameters for simulation of 1 to N independent planets around two binary stars """
  result = OptionParser()
  result.add_option("--dir", 
                    dest="dir", default="sim",
                    help="directory base name [%default]")
  result.add_option("-n", 
                    dest="n_steps", type="int", default = 10,
                    help="number of steps [%default]")
  result.add_option("--fout", 
                    dest="fout", default="system_output.hdf5",
                    help="output file [%default]")
  result.add_option("--file_redir", 
                    dest="file_redir", default=None,
                    help="output file [%default]")
  result.add_option("--a_pl", unit=units.AU,
                    dest="a_pl", type="float", default = 3.1|units.AU,
                    help="planet semi-major axis [%default]")
  result.add_option("--e_pl", 
                    dest="e_pl", type="float", default = 0,
                    help="planet eccentricity [%default]")
  result.add_option("--i_pl", 
                    dest="i_pl", type="float", default = 0.0,
                    help="planet inclination [%default]")
  result.add_option("--M_pl", 
                    dest="M_pl", type="float", default = 0,
                    help="planet mean anomaly (in deg) [%default]")
  result.add_option("--argw_pl", 
                    dest="argw_pl", type="float", default = 0.0,
                    help="planet argument of pericenter (in deg) [%default]")
  result.add_option("--node_pl", 
                    dest="node_pl", type="float", default = 0.0,
                    help="planet longitude of the ascending node (in deg) [%default]")
  result.add_option("--mass_bin", unit=units.MSun,
                    dest="mass_bin", type="float", default=1.0|units.MSun,
                    help="total mass of two binary stars [%default]")
  result.add_option("--u_bin",
                    dest="u_bin", type="float", default = 0.3,
                    help="binary star mass ratio M2/mass_b [%default]")
  result.add_option("--a_bin", unit=units.AU,
                    dest="a_bin", type="float", default = 1.0|units.AU,
                    help="binary semi-major axis [%default]")
  result.add_option("--e_bin",
                    dest="e_bin", type="float", default = 0.4,
                    help="binary eccentricity [%default]")
  result.add_option("--i_bin",
                    dest="i_bin", type="float", default = 0.1,
                    help="binary inclination in deg [%default]")
  result.add_option("--argw_bin",
                    dest="argw_bin", type="float", default = 1.0,
                    help="binary argument of pericenter in deg [%default]")
  result.add_option("--node_bin",
                    dest="node_bin", type="float", default = 0.0,
                    help="longitude of the ascending node of the binary stars in deg [%default]")
  result.add_option("--M_bin",
                    dest="M_bin", type="float", default = 0.0,
                    help="binary mean anomaly in deg [%default]")
  result.add_option("--eta",
                    dest="eta", type="float", default=0.001,
                    help="Huayno eta parameter (~timestep) [%default]") # Huayno timestep here
  result.add_option("--sim_time", unit=units.yr,
                    dest="sim_time", type="int", default=30000|units.yr,
                    help="total simulation time [%default]")

  return result
  

def execute_main(o, arguments):
   """ the main method (callable internally (pseudo_main(args)) or from command line) """
   
   # Orbital Elements
   planet_elements = Orb_Kepler(o.a_pl, o.e_pl, o.i_pl, o.M_pl, o.argw_pl, o.node_pl)
   binary_elements = Orb_Kepler(o.a_bin, o.e_bin, o.i_bin, o.M_bin, o.argw_bin, o.node_bin)
   
   # Star Masses
   mass_one = (o.mass_bin) * (1 - o.u_bin)
   mass_two = (o.mass_bin) * (o.u_bin)
      
   star_masses = quantities.AdaptingVectorQuantity()
   star_masses.append(mass_one)
   star_masses.append(mass_two)
   
   print star_masses, o.a_pl
   
   mkdir(o.dir)
   
   run_simulation(planet_elements, binary_elements, star_masses,
                    o.eta, o.sim_time, 
                    o.n_steps, o.dir, o.fout, file_redir = o.file_redir)
  
def pseudo_main(args):
   """ option parser input not from the command line """
   o, arguments  = new_option_parser().parse_args(args)
   execute_main(o, arguments)
                    
if __name__ in ('__main__'):
   o, arguments  = new_option_parser().parse_args()
   execute_main(o, arguments)
  
