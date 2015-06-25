"""
Creates a scatter plot of ejection times against final sm-axes
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
from matplotlib import cm as cmap
from matplotlib.font_manager import FontProperties

import glob
import pickle

from structures import *
from avatar import *

from mercury import G as BigG

degree_sign= u'\N{DEGREE SIGN}'
SIM_TIME = 299999 # Simulation Time in Years (really T_b)

# Used to Annotate Plot with a_crit if 'annotate' is selected
def holman_fit(mass_ratio, ecc):
    """HW99 fitting formula for the LCO"""
    a_crit = 1.60 + 5.10 * ecc - 2.22 * ecc**2 \
             + 4.12 * mass_ratio - 4.27 * mass_ratio * ecc \
             -5.09 * mass_ratio**2 + 4.61 * mass_ratio**2 * ecc**2

    return a_crit

# Used to plot vertical line at a_st
def crit_sma(mass_ratio, ecc, inc, crit_type = "safe"):
    """ dictionary-like calls for the tables of critical sm-axes """
    mass_ratios = np.array([abs(0.05 * x - mass_ratio) for x in range(2,11)]) # [0.1 - mass_ratio, 0.15 - mass_ratio, ..., 0.5 - mass_ratio]
    mu_i = np.argmin(mass_ratios) # The selected mass ratio will be '0' after mass_ratio is subtracted from the array

    ecc_s = np.array([abs(0.05 * x - ecc) for x in range(15)])
    ecc_i = np.argmin(ecc_s)

    inc_s = np.array([abs(10.0 * x - inc) for x in range(10)])
    inc_i = np.argmin(inc_s)

    # 3 Possibilities (safe = UCO from initial, moving = UCO from median, else = LCO)
    if crit_type == "safe":
        crit_type = "safe_"
    elif crit_type == "moving":
        crit_type = "moving_safe_"
    else:
        crit_type = ""

    safe_crit_sma_table = pickle.load(open("/Users/Sportsfan77777/planets/critical/sim_%scritical_sma_inc_tables.p" % crit_type, "rb"))
    return safe_crit_sma_table[inc_i, ecc_i, mu_i] # stored as [i = inc, j = ecc, k = mass_ratio]


#### BEGIN HERE ####

# Initialize Arrays
sm_axis_table = []
ejectionTable = []

# Combine Zones
zone_directories = glob.glob("zone*/")

for zone in zone_directories:
    # Read info.p
    pickle_fn = "info.p"
    pickle_f = open(pickle_fn, "rb")
    o = pickle.load(pickle_f)
    pickle_f.close()

    # Get Parameters from 'o'
    u_bin = o.u_bin
    e_bin = o.e_bin
    i_bin = o.min_inc

    # Get Table of Final Semimajor Axes
    pickle_f = open("table_of_final_sm-axes.p", "rb")
    sm_axis_table = pickle.load(pickle_f)
    pickle_f.close()

    # Get Table of Ejection Times
    pickle_f = open("table_of_ejections.p", "rb")
    ejectionTable = pickle.load(pickle_f)
    pickle_f.close()

    # Reshape Arrays into 1-D (this is only useful for plotting purposes)
    sm_axis_table_i = np.array(sm_axis_table.reshape(sm_axis_table.size))
    ejectionTable_i = np.array(ejectionTable.reshape(ejectionTable.size))

    # Concatenate Into Master Zone Arrays
    sm_axis_table.concatenate(sm_axis_table, sm_axis_table_i)
    ejectionTable.concatenate(ejectionTable, ejectionTable_i)
    

# Re-label Surviving Planets with Ejection Time = SIM_TIME 
# Also, create ejectee mask and survivor mask
x = sm_axis_table[:]
mask = np.zeros(len(ejectionTable))
crit_sma_i = 0
for i,z in enumerate(ejectionTable):
    if z == 999888.9:
        ejectionTable[i] = SIM_TIME / 1000 # Set non-ejectees apart (to ~300,000 years)
        mask[i] = 1 # Create mask to color differently in the plot
    else:
        crit_sma_i = i # Save final 'i' to track critical sm-axis

reverse_mask = np.ones(len(ejectionTable))
reverse_mask -= mask # inverted mask

# Format Ejection Times for Plot
y = [1000*i for i in ejectionTable]
#colors = np.linspace(0,1,len(x)) # To color by initial sm-axis to mark incorrect median sm-axis

# Set Up Plot
fontsize = 14

#plot.scatter(x, y, c = colors, cmap = cmap.cool)
fig = plot.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

# Ax1 is the Main Plot
# Unstable Below
ax1.scatter(np.ma.masked_array(x, mask), np.ma.masked_array(y, mask), color = 'red', edgecolor = "black")
# Stable Below
ax1.scatter(np.ma.masked_array(x, reverse_mask), np.ma.masked_array(y, reverse_mask), color = 'blue', edgecolor = "black")

ax1.set_xlabel("Median Semimajor Axis $a$ [$a_b$]", fontsize = fontsize)
ax1.set_ylabel("Ejection Time [$T_b$]", fontsize = fontsize)
plot.title("System Parameters: ( $\mu$ = $%.02f$, $e$ = $%.02f$, $i$ = $%2d^{\circ}$)" % (u_bin, e_bin, i_bin), y = 1.1)

ax1.set_xlim(o.min_sma - 0.1, o.max_sma + 0.1)
ax2.set_ylim(10, 10**(np.ceil(np.log(SIM_TIME)/np.log(10))))
ax1.set_yscale('log')

# Ax2 is the Resonant Axis
min_r = int(np.ceil((o.min_sma - 0.1)**1.5))
max_r = int(np.ceil((o.max_sma + 0.1)**1.5))
new_tick_locations = [q**(0.66667) for q in np.array(range(min_r, max_r))]
def resonant_axis_labels(a_s):
    """ set up axis label numbers """
    r_s = [a**1.5 for a in a_s]
    r_labels = ["%d" % int(round(r,0)) for r in r_s]
    return r_labels

ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(resonant_axis_labels(new_tick_locations))
ax2.set_xlabel("Orbital Period $T$ [$T_b$]", fontsize = fontsize)

## Annotate if any argument is supplied
if (len(sys.argv) > 1):
    print " *** Annotating ***"

    ### STABLE LABEL ###
    #x_st = x[crit_sma_i] # Fake Method
    x_st = crit_sma(u_bin, e_bin, i_bin) # Real Method
    print "Critical Sm-Axis:",  x_st

    # Set 'Stable' Text Height
    annotate_y = SIM_TIME / 3

    # Set 'Stable' Text Font
    font0 = FontProperties()
    font = font0.copy()
    font.set_weight("bold")
    font.set_family("Arial Narrow")
    ax1.annotate("    ALL  \nSTABLE", xy = (o.max_sma + 0.075, annotate_y), xytext = (x_st + 0.05, annotate_y),
                 arrowprops = dict(facecolor = "green", edgecolor = "green",
                                   arrowstyle = "simple, tail_width = 0.7, head_width = 1.6"),
                 verticalalignment = "center",
                 fontproperties = font, fontsize = 18, color = "green")

    #### A_ST LINE LABEL ####
    # Plot 'A_ST' Line
    x_st_line = x_st + 0.03
    ax1.plot([x_st_line, x_st_line], [10, 2.4 * SIM_TIME], color = "green", linewidth = 4, linestyle = "-")

    # Set 'A_ST' Text Height
    annotate_y_top = 1000
    annotate_y_bottom = 250

    # Set 'A_ST' Text Font
    font2 = font0.copy()
    font2.set_weight("bold")
    ax1.annotate("$a_{st}$", xy = (x_st + 0.06, annotate_y_top), xytext = (x_st + 0.25, annotate_y_bottom),
                 arrowprops = dict(facecolor = "green", edgecolor = "green",
                                   arrowstyle = "simple, tail_width = 0.1, head_width = 0.7"),
                 verticalalignment = "center", 
                 fontproperties = font2, size = 25, color = "green")

    #### A_CRIT LINE LABEL ####
    # Plot 'A_CRIT' Line
    x_crit = holman_fit(u_bin, e_bin)
    ax1.plot([x_crit, x_crit], [10, 2.4 * SIM_TIME], color = "green", linewidth = 3, linestyle = "dashed")

    # Set 'A_CRIT' Text Height
    annotate_y_top = 100
    annotate_y_bottom = 25

    # Set 'A_CRIT' Text Font
    font3 = font0.copy()
    font3.set_weight("bold")
    ax1.annotate("$a_{crit,HW99}$", xy = (x_crit + 0.03, annotate_y_top), xytext = (x_st + 0.22, annotate_y_bottom),
                 arrowprops = dict(facecolor = "green", edgecolor = "green",
                                   arrowstyle = "simple, tail_width = 0.08, head_width = 0.6"),
                 verticalalignment = "center", 
                 fontproperties = font3, size = 20, color = "green")



name = "t_vs_a_and_T_for_%s_u%02d_e%02d_i%02d" % (o.dir, 100 * u_bin, 100 * e_bin, i_bin)

plot.savefig(name + ".png", bbox_inches='tight')
plot.savefig(name + ".pdf", format = "pdf", dpi = 1000, bbox_inches='tight')
plot.show()

plot.cla()

# Save Arrays
"""
# Save Arrays for Re-Use in Combined Plots
pickle.dump(x_r, open("xr" + name, 'wb'))
pickle.dump(y, open("y" + name, 'wb'))
"""
