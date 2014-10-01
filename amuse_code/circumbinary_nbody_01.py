import numpy 
import os
import time

from amuse.units.optparse import OptionParser
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSuperset
from amuse.community.huayno.interface import Huayno
#from kepler6.interface import Kepler
from amuse.io import write_set_to_file

import binary_orbit_ini as orbit

from orbital_elements import Orb_Kepler

def rm_file(file_name):
  """
  delete the file 'file_name', if it exists
  """  
  if(os.path.isfile(file_name)):
    os.system('rm %s' % (file_name))

def orbital_period(a, Mtot) :
    """ return Keplerian orbital period from orbital elements """
    return 2*numpy.pi*(a**3/(constants.G*Mtot)).sqrt()
    
def initialize_huayno(bodies, converter, huayno_eta):
  """ initialize huayno code """
  gravity = Huayno(converter,channel_type="sockets")
  gravity.particles.add_particles(bodies) # either bodies or just stars
  gravity.commit_particles()
  #gravity.set_timestep_parameter(huayno_eta)
  gravity.parameters.timestep_parameter = huayno_eta
  gravity.set_inttype_parameter(8) # CC_KEPLER
  
  """
  Huayno Options:
  (pass / hold / bridge / maybe shared? are described in the Huayno Paper)
        SHARED2=1
        EXTRAPOLATE=5
        PASS_KDK=2
        PASS_DKD=7
        HOLD_KDK=3
        HOLD_DKD=8 <<<<-------
        PPASS_DKD=9
        BRIDGE_KDK=4
        BRIDGE_DKD=10
        CC=11
        CC_KEPLER=12 <<<<-------
        OK=13
        KEPLER=14
        SHARED4=15
        SHARED6=18
        SHARED8=19
        SHARED10=20
        SHAREDBS=21    **BS**
        CCC=22
        CCC_KEPLER=23
        CC_BS=24    **BS**
        CCC_BS=25    **BS**
        BS_CC_KEPLER=26    **BS**
        CC_BSA=27    **BS** ???
        CCC_BSA=28    **BS** ???
        SHARED2_COLLISIONS=29
        SHARED4_COLLISIONS=30
        SHARED6_COLLISIONS=31
        SHARED8_COLLISIONS=32
        SHARED10_COLLISIONS=33
  """

  return gravity
  
def initialize_circumbinary_system(star_elements, planet_elements, star_masses):
  """
  set up two_stars and one planet orbiting both
  from given orbital parameters
  """
  
  se = star_elements
  
  pe = planet_elements
  
  se_true_anom = orbit.true_anomaly_from_mean_anomaly(se.M, se.e)
  
  stars = orbit.kepler_binary_from_orbital_elements(star_masses[0], star_masses[1], 
                    se.a, se.e, se_true_anom, se.i, se.node, se.argw)
                    
  pe_true_anom = orbit.true_anomaly_from_mean_anomaly(pe.M, pe.e)
                        
  binary_mass = star_masses[0] + star_masses[1] # Is this incorrect??????
  
  planet = orbit.planet_from_orbital_elements(binary_mass,
                       pe.a, pe.e, pe_true_anom, pe.i, pe.node, pe.argw)
                       
  system = Particles(0)
  system.add_particle(stars[0])
  system.add_particle(stars[1])
  system.add_particle(planet)
                        
  return system
    
def integrate_circumbinary_system(system, t_end, n_steps,
                         snap_dir, file_out, file_redir, huayno_eta):
  """
  integrate circumbinary system
  """
  #print "Stars"
  #print stars
  #print "planet"
  #print planet
  
  #bodies = ParticlesSuperset([stars, planet])
  bodies = system
  converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  
  # Initialize Code
  gravity = initialize_huayno(bodies, converter, huayno_eta)
  
  print ' ** timestep: ', gravity.get_timestep_parameter()
  print ' ** inttype: ', gravity.get_inttype_parameter()
  print ' ** eps2: ', gravity.get_eps2_parameter(), numpy.sqrt(gravity.get_eps2_parameter())
  
  t0 = time.time()
  evolve_circumbinary_system(bodies, gravity, 
                    t_end, n_steps, converter, snap_dir, file_out)
  t1= time.time()
  dt = t1-t0
  print "Performace data: N =", len(bodies), "dt=", dt, "s =", dt/60.0, "min"
  
  return

def evolve_circumbinary_system(bodies, gravity, 
                            t_end, n_steps, converter, snap_dir, file_out):
  """
  Evolve bodies (two stars and one circum-binary planet) using 'gravity' code
  """
  #bodies = ParticlesSuperset([stars, planet])
  
  channel_from_gr_to_framework = gravity.particles.new_channel_to(bodies)
  
  rm_file(file_out)
    
  Etot_init = gravity.kinetic_energy + gravity.potential_energy
  Etot = Etot_init
  
  dt = t_end / float(n_steps)
  time = 0.0 | units.yr

  stdout = (file_out.split('.'))[0]
  stdout += '.txt'
  f = open(snap_dir + "/" + stdout, 'a')
  
  print " ** evolving: t_end = ", t_end, ", dt = ", dt.in_(units.yr)
  print " \t", "time", "\t\t\t", "E" "\t\t\t", "dE"
  while time<=t_end:
    gravity.evolve_model(time)
    channel_from_gr_to_framework.copy()
    
    bodies.collection_attributes.timestamp = time

    Ekin = gravity.kinetic_energy 
    Epot = gravity.potential_energy
    Etot = Ekin + Epot
    dE = Etot_init-Etot

    nb_E = converter.to_nbody(Etot)
    #nb_J = converter.nbody_length ** 2 * units.nbody_mass * units.nbody_time ** -2 # not supposed to work

    # A formatted string would work better than tabs. (Tabs never work)
    line = " \t" + str(time.value_in(units.yr)) + "\t" + str(nb_E) + "\t" + str(dE/Etot_init)
    print line
    
    f.write(line + "\n")

    # Write coordinates in Center of Mass frame

    # Move Stars to CoM (initially in CoM)
    #### The stars are already in CoM coordinates ####

    # Move planets to CoM (initially in CoM)
    #### The stars are already in CoM coordinates #### (Does this mean the other method is wrong?)
    
    write_set_to_file(bodies, snap_dir + "/" + file_out, "hdf5")
    
    time += dt
  
  gravity.stop()
  
  return
  
#if __name__ in ('__main__', '__plot__'):
#  o, arguments  = new_option_parser().parse_args()
  