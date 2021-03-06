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

from structures import *
from avatar import *

from mercury import G as BigG

# Make plotting optional in the future????

"""
This program converts the collision + ejection data from info.out
into a table to match the format from the Holman/Wiegert paper
showing the mean semimajor axes for each test mass.

The orbital elements are not real.
Thus, 'a' should oscillate over time as it goes through its orbit, 
even if its trajectory stays the same.
This program averages 'a' over time, neglecting the first term and final term for safety.

Certain parameters must be specified near the top of the file.
UPDATE: These are called from an "info.p" file to make things easier.

NOTE: This version uses the old format without IDs!!!!!
"""

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

""" Collect *.aei files (output data from e.exe = mercury output parser) """

aei_path = "*_*.aei"  # because the formt is M#_S#.aei
aei_files = sorted(glob.glob(aei_path))

if len(aei_files) == 0:
    # ./e.exe has not been run yet, so run it
    element_output = "./e.exe"
    subprocess.call(element_output)

    # Re-try globglob
    aei_files = sorted(glob.glob(aei_path))


""" For each .aei file, get only the value for 'a' """

# store each a_over_time in dictionary corresponding to IDs (or orbital parameters??????)
names = []
a_dict = {}
t_dict = {}
for aei_fn in aei_files:
    name = aei_fn[:aei_fn.rfind(".")] # for 'ID_0000.aei', return 'ID_0000'
    names.append(name)

    a_over_time = []
    time_array = []

    f = open(aei_fn)
    lines = f.readlines()

    for line in lines:
        split_line = line.split()

        if len(split_line) > 2:
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
               
               a = calculate_a(position, velocity, mu)

               a_over_time.append(a)
               time_array.append(t)

    if len(a_over_time) >= 5:
       time_array = time_array[1:-1]
       a_over_time = a_over_time[1:-1] # get rid of first and last entry IFF eject after 5000 years

    a_dict[name] = a_over_time
    t_dict[name] = time_array

print "Names:", names
        
""" 
Format: ID_#%04d collided with the central body at #c.#d years
Format: ID_#%04d ejected at #c.#d years

Read ID_# and look up corresponding mean_anomaly and semimajor_axis, #c as ejection time in 1000s of years
"""


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
        semi_major_axes[i] = round(semi_major_axes[i], 1) # <<--- Note the change from x.x to x.xx!!

sm_array = [int(10 * x) for x in semi_major_axes] # <<<<---- ARCHAIC use!!!!!!
#print sm_array # S Names

sm_axis_table = np.zeros((num_a, N)) + 99.9

# (1) Calculate Mean SMAs, (2) Plot Them, and (3) Parse Them Into Table

for name in names:
    name_split = name.split('_')

    M_str = name_split[0]
    S_str = name_split[1]

    Mi = m_deg_array.index(int(M_str[1:]))
    Si = sm_array.index(int(S_str[1:]))

    this_a_over_time = a_dict[name] # retrieve
    this_time = t_dict[name] # retrieve

    #print this_a_over_time
    #print this_time

    if len(this_a_over_time) > 0:
        plot_fn = "%s_sm-axis_evolution.png" % name # Note: ID_name = e.g. A2.1
        plot_title = "Planet: %s" % (name)

        plot.plot(this_time, this_a_over_time)
        #plot.plot(this_time, this_a_over_time, 'ro')

        plot.title(plot_title)
        plot.xlabel("time (in years)")
        plot.ylabel("semimajor axis (in scaled AU)")

        plot.ylim(0.0, 8.0) # in AU

        plot.savefig(plot_fn)

        plot.cla() # CLEAR FIGURE!

        median_sma = np.median(this_a_over_time)
        three_sigma = 3.0 * np.std(this_a_over_time)

        outliers = True
        tmp_a_over_time = this_a_over_time[:] 
        while outliers:
            #print median_sma, three_sigma
            #print

            prev_tmp = tmp_a_over_time[:]
            tmp_a_over_time = []
            for a in prev_tmp:
                if abs(a - median_sma) > three_sigma:
                    pass
                else:
                    tmp_a_over_time.append(a)

            median_sma = np.median(tmp_a_over_time)
            three_sigma = 3.0 * np.std(tmp_a_over_time)

            if len(tmp_a_over_time) == len(prev_tmp):
                outliers = False

        sm_axis_table[Si][Mi] = median_sma   ### <<<<<----- Currently working on this

# Simple Print  
#print ejectionTable, '\n'

# Pretty Print
width = 7
rowzero = phase_str[phase].center(width) + '|'
dash = '-' * width
dash_row = dash
for x in m_deg_array:
    rowzero += (str(x)).center(width)
    dash_row += dash
    
rows = []
count = 0
str_y_base = "%.1f"
for i,y in enumerate(sm_array):
    str_y = str_y_base % (y / 10.0)
    row = str_y.center(width) + '|'
    for median_sma in sm_axis_table[i]:
        s = ""
        if median_sma == 0.0:
            s = "QUICK".center(width)
        else:
            s = str("%.2f" % median_sma).center(width)
        row += s
    rows.append(row)
    
# Write Pretty Print to File
if len(sys.argv) > 1:
   fn = sys.argv[1]  # write to a file (if provided)
   f = open(fn, 'a')
else:
   f = open("sm-axes.t", 'w')
    
#header = "Directory: %s\n" % o.integration_dir
#f.write(header)
print rowzero
f.write(rowzero + "\n")
print dash_row
f.write(dash_row + "\n")
for r_str in rows:
    print r_str
    f.write(r_str + "\n")

f.close()

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

