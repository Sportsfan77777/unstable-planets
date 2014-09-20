"""
Integrator of Circumbinary Planetary Systems
"""




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
                    dest="num_a", type="int", default = 8,
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
                    dest="min_ecc", type="int", default = 0,
                    help="planet eccentricity minimum [%default]")
  result.add_option("--max_ecc", 
                    dest="max_ecc", type="int", default = 1,
                    help="planet eccentricity maximum [%default]")
  result.add_option("--min_inc", 
                    dest="min_inc", type="int", default = 0,
                    help="planet inclination minimum in deg [%default]")
  result.add_option("--max_inc", 
                    dest="max_inc", type="int", default = 360,
                    help="planet inclination maximum in deg [%default]")
  result.add_option("--mass_b", units=units.MSun,
                    dest="mass_b", type="float", default=1.0|units.MSun,
                    help="total mass of two binary stars [%default]")
  result.add_option("--u_b",
                    dest="u_b", type="float", default = 0.3,
                    help="binary star mass ratio M2/mass_b [%default]")
  result.add_option("--a_b", unit=units.AU,
                    dest="a_b", type="float", default = 1.0|units.MSun,
                    help="mass of the passing star (with no disk) in MSun [%default]")
  result.add_option("--e_b", unit=units.AU,
                    dest="e_b", type="float", default = 0.4,
                    help="eccentricity of the binary stars [%default]")
  result.add_option("--argw_b",
                    dest="argw_b", type="float", default = 1.0,
                    help="argument of pericenter of the binary stars in deg [%default]")
  result.add_option("--node_b",
                    dest="node_b", type="float", default = 0.0,
                    help="longitude of the ascending node of the binary stars in deg [%default]")
  result.add_option("--mean_anom_b",
                    dest="mean_anom_b", type="float", default = 0.0,
                    help="mean anomaly of the binary stars in deg [%default]") # was 90.0 before
  result.add_option("--i_b",
                    dest="i_b", type="float", default = 0.1,
                    help="inclination of the binary stars in deg [%default]")
  result.add_option("--eta",
                    dest="eta", type="float", default=0.001,
                    help="Huayno eta parameter (~timestep) [%default]") # Huayno timestep here
  ###### EXPERIMENT #####   not used as far as I know
  result.add_option("--br_dt",
                    dest="br_dt", type="float", default=0.01,
                    help="bridge timestep -- fraction of inner disk period [%default]") # Bridge timestep here
  ###### EXPERIMENT #####   not used as far as I know
  result.add_option("--num_orbits", units=units.yr
                    dest="num_orbits", type="int", default=30000|units.yr,
                    help="total simulation time (in 'binary star orbit' units) [%default]")

  return result
  
