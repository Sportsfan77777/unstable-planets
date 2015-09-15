"""
Creates a scatter plot of the minimum ejection times needed 
to find the critical semimajor critical semimajor axis
across the entire parameter space.

This is useful for knowing how long is necessary to run simulations.
"""

import numpy as np
import string
import math
import sys
import os
import subprocess

import matplotlib
#matplotlib.use('Agg') # for ssh-ed computers
from matplotlib import pyplot as plot

import glob
import pickle

from mercury import G as BigG

# Get List of Critical Semimajor Axes
pickle_f = open("sim_safe_critical_sma_inc_list.p", "rb")
crit_sm_axis_list = pickle.load(pickle_f)
pickle_f.close()

# Get List of Ejection Times
pickle_f = open("sim_safe_critical_sma_inc_eject_list.p", "rb")
ejectionTime_list = pickle.load(pickle_f)
pickle_f.close()

# Reshape to 1-D
ejectionTime_list = ejectionTime_list.reshape(np.size(ejectionTime_list))
crit_sm_axis_list = crit_sm_axis_list.reshape(np.size(crit_sm_axis_list))

x = ejectionTime_list[:]
y = crit_sm_axis_list[:]

# Set Up Plots
fontsize = 16
labelsize = 14
rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

min_y = 1.8
max_y = 7.4

color = "#B22222"
line_color = "#0000CD"
alpha_val = 0.28

size = 30
width = 3

### SUB PLOTS ###
f, (ax1, ax2) = plot.subplots(2, sharex=True)
f.subplots_adjust(hspace=0)

######## SCATTER ########

ax1.plot([300, 300], [min_y, max_y], color = line_color, linewidth = width)
ax1.scatter(x, y, color = color, alpha = alpha_val, marker = "o", s = size)

ax1.set_xlim(-20, 400)
ax1.set_ylim(min_y, max_y)

#ax1.set_xlabel("Sufficient Simulation Time [$10^3 T_b$] at $a_{st}$", fontsize = fontsize)
ax1.set_ylabel("$a_{st}$ [$a_b$]", fontsize = fontsize)
#plot.title("Simulation Times Sufficient for Finding $a_{st}$")

##### ADD HISTOGRAM #####
#plot.title("Distribution of Sufficient Simulation Times")
ax2.set_xlabel("Required Simulation Time [$10^3 T_b$]", fontsize = fontsize)
ax2.set_ylabel("Fraction of Simulations", fontsize = fontsize)

ax2.set_xlim(-5, 300)
ax2.set_ylim(0.00, 1.02)

ax2.hist(x, bins = range(0, 301, 25), normed = True, cumulative = True, histtype = "bar", color = color)

plot.show()
plot.savefig("combinedNescessity.png")
plot.savefig("combinedNescessity.pdf", format = 'pdf', dpi = 1000)

###### Quantities ######
total = len(y)

less_than_10 = len(x[x < 10])
less_than_50 = len(x[x < 50])
less_than_100 = len(x[x < 100])
less_than_200 = len(x[x < 200])

fraction_10 = str(1.0 * less_than_10 / total)
fraction_50 = str(1.0 * less_than_50 / total)
fraction_100 = str(1.0 * less_than_100 / total)
fraction_200 = str(1.0 * less_than_200 / total)

print less_than_10, total, fraction_10
print less_than_50, total, fraction_50
print less_than_100, total, fraction_100
print less_than_200, total, fraction_200

duration_f = open("duration_values.txt", "w")

duration_f.write("Total Simulations: " + str(total) + "\n")
duration_f.write("Fraction of Simulations that need  10 kyr: " + fraction_10 + "\n")
duration_f.write("Fraction of Simulations that need  50 kyr: " + fraction_50 + "\n")
duration_f.write("Fraction of Simulations that need 100 kyr: " + fraction_100 + "\n")
duration_f.write("Fraction of Simulations that need 200 kyr: " + fraction_200 + "\n")

duration_f.close()



