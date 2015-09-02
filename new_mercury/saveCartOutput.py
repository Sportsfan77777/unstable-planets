"""
This program converts the cartesian .aei files
into binary files that actually work and store all of the cartesian elements over time.

A different file is saved for each planet.
"""

import numpy as np
import string
import math
import sys
import os
import subprocess

import glob
import pickle

from id import ID_Manager
from structures import *
from avatar import *

from mercury import G as BigG


def isFloat(value):
    try:
      float(value)
      return True
    except ValueError:
      return False

""" Read from info.p file (parameters) """

pickle_fn = "info.p"
pickle_f = open(pickle_fn, "rb")
o = pickle.load(pickle_f)
pickle_f.close()

phase_str = ['Peri', 'Apo']
if o.mean_anom_bin == 0:
   phase = 0
else:
   phase = 1   

N = o.num_M
num_a = o.num_a

a_b = o.a_bin

a_min = o.min_sma
a_max = o.max_sma

mu = BigG * o.mass_bin # Used in calculate orbital elements from cartesian output files

""" Read IDs from ids.p file (ids map to the old format of M#a_S_ """

pickle_fn = "ids.p"
pickle_f = open(pickle_fn, "rb")
id_dict = pickle.load(pickle_f)
pickle_f.close()

ID_manager = ID_Manager()
ID_manager.read()

""" Collect *.aei files (output data from e.exe = mercury output parser) """

aei_path = "ID*.aei"
aei_files = sorted(glob.glob(aei_path))

if len(aei_files) == 0:
    # ./e.exe has not been run yet, so run it
    element_output = "./e.exe"
    subprocess.call(element_output)

    # Re-try globglob
    aei_files = sorted(glob.glob(aei_path))

""" For each .aei file, get all the xyz,uvw and convert to aeiWwM """

############# BARYCENTRIC CARTESIAN FORMAT!!!!!!!! (because Jacobi thing doesn't work...) ##########

# store each a_over_time in dictionary corresponding to IDs (or orbital parameters??????)
for aei_fn in aei_files:
    id_name = aei_fn[:aei_fn.rfind(".")] # for 'ID_0000.aei', return 'ID_0000'
    save_fn = "%s_cart.npy" % id_name

    positions_over_time = []
    velocities_over_time = []
    time_array = []

    if os.path.exists(save_fn):
        # This has already been done!
        pass
    else:
        elements_over_time = []
        time_array = []

        f = open(aei_fn)
        lines = f.readlines()

        for line in lines:
            split_line = line.split()

            if len(split_line) >= 7:
                t = split_line[0]

                x = split_line[1]
                y = split_line[2]
                z = split_line[3]

                u = split_line[4]
                v = split_line[5]
                w = split_line[6]
                

                if (isFloat(x) and isFloat(y) and isFloat(z) and isFloat(u) and isFloat(v) and isFloat(w)):
                   position = Position(float(x), float(y), float(z))
                   velocity = Velocity(float(u), float(v), float(w))

                   positions_over_time.append(position)
                   velocities_over_time.append(velocity)
                   time_array.append(t)

        timesteps = len(time_array)
        final_array = np.zeros((timesteps, 7)) # t, a, e, i, W, w, M   ### or do t, elements????

        # ---->>>> save time array and element array separately for each planet???? <<<<<----

        for i, (position, velocity, time) in enumerate(zip(positions_over_time, velocities_over_time, time_array)):
        	final_array[i, 0] = time
        	final_array[i, 1] = position.get('x')
        	final_array[i, 2] = position.get('y')
        	final_array[i, 3] = position.get('z')
        	final_array[i, 4] = velocity.get('vx')
        	final_array[i, 5] = velocity.get('vy')
        	final_array[i, 6] = velocity.get('vz')

        np.save(save_fn, final_array)

# Delete .aei files
for aei_file in aei_files:
	rm_command = ['rm', aei_file]
	subprocess.call(rm_command)
