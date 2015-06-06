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
crit_sm_axis_list = crit_sm_axis_list.reshape(np.size(crit_sm_axis_list))
ejectionTime_list = ejectionTime_list.reshape(np.size(ejectionTime_list))

x = crit_sm_axis_list[:]
y = ejectionTime_list[:]

# Set Up Plot
min_x = 1.8
max_x = 7.4

color = "#B22222"
line_color = "#0000CD"
alpha_val = 0.35

size = 30
width = 3

plot.plot([min_x, max_x], [300, 300], color = line_color, linewidth = width)
plot.scatter(x, y, color = color, alpha = alpha_val, marker = "o", s = size)

plot.xlim(min_x, max_x)
plot.ylim(-20, 500)

plot.xlabel("$a_{st}$ (in units of binary separation)")
plot.ylabel("Minimum $t_{eject}$ near $a_{st}$ (in units of thousands of binary periods)")
plot.title("Simulation Times Necessary for Finding $a_{st}$")

plot.savefig("necessarySimulationTimes.png")
plot.show()

plot.cla()

# Set Up Log Plot
plot.plot([min_x, max_x], [300, 300], color = line_color, linewidth = width)
plot.scatter(x, y, color = color, alpha = alpha_val, marker = "o", s = size)

plot.xlim(min_x, max_x)
#plot.ylim(0.01, 500)

plot.xlabel("$a_{st}$ (in units of binary separation)")
plot.ylabel("Minimum $t_{eject}$ near $a_{st}$ (in units of thousands of binary periods)")
plot.title("Simulation Times Necessary for Finding $a_{st}$")

plot.yscale('log')
plot.ylim(0.01, 500)
plot.savefig("necessarySimulationTimes_log.png")
plot.show()


# Quantities
total = len(x)

less_than_50 = len(y[y < 50])
less_than_100 = len(y[y < 100])
less_than_200 = len(y[y < 200])

fraction_50 = str(1.0 * less_than_50 / total)
fraction_100 = str(1.0 * less_than_100 / total)
fraction_200 = str(1.0 * less_than_200 / total)

print less_than_50, total, fraction_50
print less_than_100, total, fraction_100
print less_than_200, total, fraction_200

duration_f = open("duration_values.txt", "w")

duration_f.write("Total Simulations: " + str(total) + "\n")
duration_f.write("Fraction of Simulations that need  50 kyr: " + fraction_50 + "\n")
duration_f.write("Fraction of Simulations that need 100 kyr: " + fraction_100 + "\n")
duration_f.write("Fraction of Simulations that need 200 kyr: " + fraction_200 + "\n")

duration_f.close()
