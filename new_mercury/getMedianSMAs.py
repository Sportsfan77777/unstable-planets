import numpy as np
import string
import math
import sys
import os
import subprocess

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plot

import glob
import pickle

from id import ID_Manager
from structures import *
from avatar import *

from mercury import G as BigG
from getEjectionData import STABLE_VALUE

"""
This program makes a table of the median semi_major_axes for each planet.
There are 3 tables produced.
One contains all planets

This program converts the collision + ejection data from info.out
into a table to match the format from the Holman/Wiegert paper
showing the mean semimajor axes for each test mass.

The orbital elements are not real.
Thus, 'a' should oscillate over time as it goes through its orbit, 
even if its trajectory stays the same.
This program averages 'a' over time, neglecting the first term and final term for safety.

Certain parameters must be specified near the top of the file.
UPDATE: These are called from an "info.p" file to make things easier.
"""

def isFloat(value):
    try:
      float(value)
      return True
    except ValueError:
      return False

""" Check args for make_plots """
make_plots = True
if len(sys.argv) > 1:
    make_plots = False

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

""" Collect all ID names using .npy files"""

npy_path = "ID*.npy"
npy_files = sorted(glob.glob(npy_path))

if len(npy_files) == 0:
    # ./e.exe has not been run yet, so run it
    command = ['python', 'saveAvatarOutput.py']
    subprocess.call(command)

    # Re-try globglob
    npy_files = sorted(glob.glob(npy_path))

id_names = []
for npy_fn in npy_files:
    id_name = npy_fn[:npy_fn.rfind("_")]
    id_names.append(id_name)

""" M """
mean_anomalies = np.ones(N) * 2 * np.pi / N
for i in xrange(N):
    mean_anomalies[i] *= i  # i = 0 to N - 1
    
m_deg_array = [int(round(x * 180.0 / np.pi)) for x in mean_anomalies]
#m_deg_array = [("03d" % x) for x in mean_anomalies_array]  # <<<<<------- FIX THIS!!!!!!! (for id.py, not here...)

#print m_deg_array # M Names

""" a """
semi_major_axes = np.ones(num_a) * a_min
if num_a > 1:
    a_range = a_max - a_min
    step_size = a_range / (num_a - 1) 
    for i in xrange(num_a):
        semi_major_axes[i] += (i * step_size)
        semi_major_axes[i] = round(semi_major_axes[i], o.sep_sma) # <<--- Note the change from x.x to x.xx!!

sm_array = [round(x, o.sep_sma) for x in semi_major_axes] # <<<<---- ARCHAIC use!!!!!!
#print sm_array # S Names

sm_axis_table = np.zeros((num_a, N)) + STABLE_VALUE

# (1) Calculate Mean SMAs, (2) Plot Them, and (3) Parse Them Into Table

for id_name in id_names:
    save_fn = "%s_elements.npy" % id_name
    planet_name = ID_manager.get_name(id_name)

    id_split = planet_name.split('_')
    A_str = id_split[0]
    M_str = id_split[1]

    Ai = sm_array.index(float(A_str[1:]))
    Mi = m_deg_array.index(int(M_str[1:]))

    if os.path.exists(save_fn):
        final_array = np.load(save_fn)

        this_time = final_array[:, 0]
        this_a_over_time = final_array[:, 1]
        this_e_over_time = final_array[:, 2]

        if make_plots:
            # de-triggered by supplying an argument to the function call
            plot_fn = "%s_sm-axis_evolution.png" % id_name # Note: ID_name = e.g. A2.1
            plot_title = "%s, Planet: %s" % (id_name, planet_name)

            plot.plot(this_time, this_a_over_time)
            plot.plot(this_time, this_a_over_time, 'ro', alpha = 0.15)

            plot.title(plot_title)
            plot.xlabel("time (in years)")
            plot.ylabel("semimajor axis (in scaled AU)")

            plot.savefig(plot_fn)

        # Filter bad values too far from initial_sma
        initial_sma = float(A_str[1:])

        # If it's out of this range, something is horribly wrong... (hopefully, it is ejecting)
        min_sma = initial_sma - 0.5
        max_sma = initial_sma + 0.5

        # Filter according to 0 < 'e' < 1 and min_sma < 'a' < max_sma
        filtered_a_over_time = this_a_over_time[(this_e_over_time >= 0.0) & (this_e_over_time < 1.0)]
        filtered_a_over_time = filtered_a_over_time[(filtered_a_over_time >= min_sma) & (filtered_a_over_time <= max_sma)]

        median_sma = np.median(filtered_a_over_time)
        sm_axis_table[Ai][Mi] = median_sma   ### <<<<<----- Currently working on this (but it is better now)

# Get Ejection Table
pickle_f = open("table_of_ejections.p", "rb")
ejectionTable = pickle.load(pickle_f)

# Pretty Print
width = 7
rowzero = phase_str[phase].center(width) + '|'
dash = '-' * width
dash_row = dash
for x in m_deg_array:
    rowzero += (str(x)).center(width)
    dash_row += dash
    
rows = []
filtered_stable_rows = []
filtered_unstable_rows = []
count = 0
str_y_base = "%.0" + ("%d" % o.sep_sma) + "f"
for i,y in enumerate(sm_array):
    str_y = str_y_base % y
    row = str_y.center(width) + '|'
    filtered_stable_row = str_y.center(width) + '|'
    filtered_unstable_row = str_y.center(width) + '|'

    all_stable = True
    all_unstable = True

    for ej_t, mean_sma in zip(ejectionTable[i], sm_axis_table[i]):
        s = ""
        f_stable_s = ""
        f_unstable_s = ""

        # (1) set s
        if mean_sma == 0.0:
            s = "QUICK".center(width)
        else:
            s = str("%.2f" % mean_sma).center(width)

        # (2, 3) set f_stable_s and f_unstable_s
        if ej_t == STABLE_VALUE:
            # Stable
            all_unstable = False

            f_stable_s = s
            f_unstable_s = "+".center(width)
        else:
            # Unstable
            all_stable = False

            f_stable_s = "-".center(width)
            f_unstable_s = s

        row += s
        filtered_stable_row += f_stable_s
        filtered_unstable_row += f_unstable_s
    rows.append(row)
    if not all_stable:
        filtered_unstable_rows.append(filtered_unstable_row)
    if not all_unstable:
        filtered_stable_rows.append(filtered_stable_row)
    
# Write Pretty Print to Files
f = open("sm-axes.t", 'w')
g = open("sm-axes_stable.t", 'w')
h = open("sm-axes_unstable.t", 'w')

# Unfiltered    
header = "Directory: %s" % o.integration_dir
print header
f.write(header + "\n")
print rowzero
f.write(rowzero + "\n")
print dash_row
f.write(dash_row + "\n")
for r_str in rows:
    print r_str
    f.write(r_str + "\n")
f.close()

# Break
print

# Filtered: Only show planets that survive
header = "Directory: %s" % o.integration_dir
print header
g.write(header + "\n")
print rowzero
g.write(rowzero + "\n")
print dash_row
g.write(dash_row + "\n")
for r_str in filtered_stable_rows:
    print r_str
    g.write(r_str + "\n")
g.close()

# Break
print

# Filtered: Only show planets that eject
header = "Directory: %s" % o.integration_dir
print header
h.write(header + "\n")
print rowzero
h.write(rowzero + "\n")
print dash_row
h.write(dash_row + "\n")
for r_str in filtered_unstable_rows:
    print r_str
    h.write(r_str + "\n")
h.close()

# Write SM-Axes Array

pickle_f = open("sm_axes.p", "wb")
pickle.dump(semi_major_axes, pickle_f)
pickle_f.close()

# Write Mean Anomaly Array

pickle_f = open("mean_anomalies.p", "wb")
pickle.dump(m_deg_array, pickle_f)
pickle_f.close()

# Write Ejection Table!
# (i) as table
pickle_f = open("table_of_final_sm-axes.p", "wb")
pickle.dump(sm_axis_table, pickle_f)
pickle_f.close()

# (ii) as array of dictionaries

