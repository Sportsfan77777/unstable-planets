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

# Get Table of Final Semimajor Axes
pickle_f = open("sim_safe_critical_sma_inc_list.p", "rb")
crit_sm_axis_list = pickle.load(pickle_f)
pickle_f.close()

# Get Table of Ejection Times
pickle_f = open("sim_safe_critical_sma_inc_eject_list.p", "rb")
ejectionTime_list = pickle.load(pickle_f)
pickle_f.close()

x = crit_sm_axis_list[:]
y = ejectionTime_list[:]

# Set Up Plot
plot.scatter(x, y, color = 'r', alpha = 0.7)

plot.xlim(1.8, 7.5)
plot.ylim(-20, 500)

plot.xlabel("$a_{st}$ (in units of binary separation)")
plot.ylabel("$t_{eject}$ near $a_{st}$ (in units of thousands of binary periods)")
plot.title("Simulation Times Necessary for Finding $a_{st}$")

plot.savefig("necessarySimulationTimes.png")
plot.show()

plot.cla()

# Set Up Log Plot

plot.scatter(x, y, color = 'r', alpha = 0.7)

plot.xlim(1.8, 7.5)
#plot.ylim(-20, 500)

plot.xlabel("$a_{st}$ (in units of binary separation)")
plot.ylabel("$t_{eject}$ near $a_{st}$ (in units of thousands of binary periods)")
plot.title("Simulation Times Necessary for Finding $a_{st}$")

plot.yscale('log')
plot.savefig("necessarySimulationTimes_log.png")
plot.show()

