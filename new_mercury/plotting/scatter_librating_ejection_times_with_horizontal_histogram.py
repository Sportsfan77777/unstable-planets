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
from matplotlib import gridspec
from matplotlib import rcParams as rc
from matplotlib.ticker import AutoMinorLocator

import glob
import pickle

from mercury import G as BigG

# Get List of Critical Semimajor Axes
pickle_f = open("sim_safe_librating_critical_sma_inc_list.p", "rb")
crit_sm_axis_list = pickle.load(pickle_f)
pickle_f.close()

# Get List of Ejection Times
pickle_f = open("sim_safe_librating_critical_sma_inc_eject_list.p", "rb")
ejectionTime_list = pickle.load(pickle_f)
pickle_f.close()

# Reshape to 1-D
crit_sm_axis_list = crit_sm_axis_list.reshape(np.size(crit_sm_axis_list))
ejectionTime_list = ejectionTime_list.reshape(np.size(ejectionTime_list))

x = crit_sm_axis_list[:]
y = ejectionTime_list[:]

# Filter non-librating data
FILTER_VALUE = -1
x = x[x != FILTER_VALUE]
y = y[y != FILTER_VALUE]

# Set Up Plots
fontsize = 17
labelsize = 15
rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

min_x = 1.8
max_x = 7.4

min_y = -6
max_y = 400

color = "#B22222"
line_color = "#0000CD"
subline_color = "#6495ED"
alpha_val = 0.28

dashes = [10, 4, 4, 4]

size = 30
width = 3

### SUB PLOTS ###
#f, (ax1, ax2) = plot.subplots(1, 2, sharey=True)
#f.subplots_adjust(wspace=0.2)
f = plot.figure(figsize =(12, 6))
gs = gridspec.GridSpec(1,2, width_ratios=[7, 4])
ax1 = f.add_subplot(gs[0])
ax2 = f.add_subplot(gs[1], sharey=ax1)

######## SCATTER ########

ax1.plot([min_x, max_x], [50, 50], color = subline_color, linewidth = width, linestyle = "-.", dashes = dashes)
ax1.plot([min_x, max_x], [200, 200], color = subline_color, linewidth = width, linestyle = "-.", dashes = dashes)
ax1.plot([min_x, max_x], [300, 300], color = line_color, linewidth = width)
ax1.scatter(x, y, color = color, alpha = alpha_val, marker = "o", s = size)

ax1.set_xlim(min_x, max_x)
ax1.set_ylim(min_y, max_y)
ax1.yaxis.set_ticks(np.arange(0, 401, 25))

ax1.set_xlabel(r"$a_{\rm st}$ [$a_{\rm b}$]", fontsize = fontsize)
ax1.set_ylabel(r"Required Integration Time [$10^3 T_{\rm b}$]", fontsize = fontsize)
#plot.title("Simulation Times Sufficient for Finding $a_{st}$")

# Set up Orbital Period on y-axis of ax1
ax3 = ax1.twiny()
min_r = int(np.ceil((min_x - 0.1)**1.5))
max_r = int(np.ceil((max_x - 0.1)**1.5))
new_tick_locations = [q**(0.66667) for q in np.array(range(min_r, max_r))]
def resonant_axis_labels(a_s):
    """ set up axis label numbers """
    r_s = [a**1.5 for a in a_s]
    r_labels = ["%d" % int(round(r,0)) for r in r_s]
    return r_labels

ax3.set_xlim(ax1.get_xlim())
ax3.set_xticks(new_tick_locations)
ax3.set_xticklabels(resonant_axis_labels(new_tick_locations))
ax3.set_xlabel(r"Orbital Period [$T_{\rm b}$]", fontsize = fontsize, labelpad = 7)

##### ADD HISTOGRAM #####
#plot.title("Distribution of Sufficient Simulation Times")
ax2.set_xlabel("Fraction of Integrations", fontsize = fontsize)
#ax2.set_ylabel("Required Simulation Time [$10^3 T_b$]", fontsize = fontsize)

ax2.set_xlim(0.00, 1.02)
ax2.xaxis.set_minor_locator(AutoMinorLocator())

ax2.plot([-1, 2], [50, 50], color = subline_color, linewidth = width, linestyle = "-.", dashes = dashes)
ax2.plot([-1, 2], [200, 200], color = subline_color, linewidth = width, linestyle = "-.", dashes = dashes)
ax2.plot([-1, 2], [300, 300], color = line_color, linewidth = width)
ax2.hist(y, bins = range(0, 401, 5), normed = True, cumulative = True, histtype = "bar", color = color, orientation = "horizontal")

plot.savefig("librating_combinedNecessity.png")
plot.savefig("librating_combinedNecessity.pdf", format = 'pdf', dpi = 1000)
plot.show()

###### Quantities ######
total = len(y)

less_than_10 = len(y[y < 10])
less_than_50 = len(y[y < 50])
less_than_100 = len(y[y < 100])
less_than_200 = len(y[y < 200])

fraction_10 = str(1.0 * less_than_10 / total)
fraction_50 = str(1.0 * less_than_50 / total)
fraction_100 = str(1.0 * less_than_100 / total)
fraction_200 = str(1.0 * less_than_200 / total)

print less_than_10, total, fraction_10
print less_than_50, total, fraction_50
print less_than_100, total, fraction_100
print less_than_200, total, fraction_200

duration_f = open("librating_duration_values.txt", "w")

duration_f.write("Total Simulations: " + str(total) + "\n")
duration_f.write("Fraction of Simulations that need  10 kyr: " + fraction_10 + "\n")
duration_f.write("Fraction of Simulations that need  50 kyr: " + fraction_50 + "\n")
duration_f.write("Fraction of Simulations that need 100 kyr: " + fraction_100 + "\n")
duration_f.write("Fraction of Simulations that need 200 kyr: " + fraction_200 + "\n")

duration_f.close()



