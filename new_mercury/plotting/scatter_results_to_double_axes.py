"""
Creates a scatter plot of ejection times against final sm-axes

Arguments:
a_crit (optional) [or use non-float argument to annotate plot]
T_st (optional)
"""

import numpy as np
import string
import math
import sys
import os
import subprocess

import matplotlib
#matplotlib.use('Agg') # for ssh-ed computers
from matplotlib import rcParams as rc
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

# Used to check if input is float (for setting a_crit and a_st in the plot)
def isFloat(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

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

# Read info.p
pickle_fn = "info.p"
pickle_f = open(pickle_fn, "rb")
o = pickle.load(pickle_f)
pickle_f.close()

# Get Parameters from 'o'
u_bin = o.u_bin
e_bin = o.e_bin
i_bin = o.min_inc

##### BASE #####

# Get Table of Final Semimajor Axes
pickle_f = open("table_of_final_sm-axes.p", "rb")
sm_axis_table = pickle.load(pickle_f)
pickle_f.close()

# Get Table of Ejection Times
pickle_f = open("table_of_ejections.p", "rb")
ejectionTable = pickle.load(pickle_f)
pickle_f.close()

# Reshape Arrays into 1-D (this is only useful for plotting purposes)
sm_axis_table = sm_axis_table.reshape(sm_axis_table.size)
ejectionTable = ejectionTable.reshape(ejectionTable.size)

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
y = np.array([1000*i for i in ejectionTable])
#colors = np.linspace(0,1,len(x)) # To color by initial sm-axis to mark incorrect median sm-axis

##### FINE #####

# Get Table of Final Semimajor Axes
pickle_f = open("table_of_final_sm-axes_fine.p", "rb")
sm_axis_table_fine = pickle.load(pickle_f)
pickle_f.close()

# Get Table of Ejection Times
pickle_f = open("table_of_ejections_fine.p", "rb")
ejectionTable_fine = pickle.load(pickle_f)
pickle_f.close()

# Reshape Arrays into 1-D (this is only useful for plotting purposes)
sm_axis_table_fine = sm_axis_table_fine.reshape(sm_axis_table_fine.size)
ejectionTable_fine = ejectionTable_fine.reshape(ejectionTable_fine.size)

# Re-label Surviving Planets with Ejection Time = SIM_TIME 
# Also, create ejectee mask and survivor mask
x_fine = sm_axis_table_fine[:]
mask_fine = np.zeros(len(ejectionTable_fine))
crit_sma_i = 0
for i,z in enumerate(ejectionTable_fine):
    if z == 999888.9:
        ejectionTable_fine[i] = SIM_TIME / 1000 # Set non-ejectees apart (to ~300,000 years)
        mask_fine[i] = 1 # Create mask to color differently in the plot
    else:
        crit_sma_i = i # Save final 'i' to track critical sm-axis

reverse_mask_fine = np.ones(len(ejectionTable_fine))
reverse_mask_fine -= mask_fine # inverted mask

# Format Ejection Times for Plot
y_fine = np.array([1000*i for i in ejectionTable_fine])
#colors = np.linspace(0,1,len(x)) # To color by initial sm-axis to mark incorrect median sm-axis

# Set Up Plot
fontsize = 17
labelsize = 15
rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

fig = plot.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

# Ax1 is the Main Plot
alpha_val = 0.17
# Unstable Below
ax1.scatter(np.ma.masked_array(x, mask), np.ma.masked_array(y, mask), color = 'red', edgecolor = "black", zorder = 10)
ax1.scatter(np.ma.masked_array(x_fine, mask_fine), np.ma.masked_array(y_fine, mask_fine), color = 'red', edgecolor = "black", alpha = alpha_val, zorder = 2)
# Stable Below
ax1.scatter(np.ma.masked_array(x, reverse_mask), np.ma.masked_array(y, reverse_mask), color = 'blue', edgecolor = "black", zorder = 10)
ax1.scatter(np.ma.masked_array(x_fine, reverse_mask_fine), np.ma.masked_array(y_fine, reverse_mask_fine), color = 'blue', edgecolor = "black", alpha = alpha_val, zorder = 2)

ax1.set_xlabel(r"Median Semimajor Axis [$a_{\rm b}$]", fontsize = fontsize)
ax1.set_ylabel(r"Ejection Time [$T_{\rm b}$]", fontsize = fontsize)
plot.title(" $\mu$ = $%.02f$, $e$ = $%.02f$, $i$ = $%2d^{\circ}$" % (u_bin, e_bin, i_bin), y = 1.13, bbox=dict(facecolor = 'none', edgecolor = 'black', linewidth = 1.5, pad=7.0), fontsize = fontsize + 2)

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
ax2.set_xlabel(r"Orbital Period [$T_{\rm b}$]", fontsize = fontsize, labelpad = 7)

## Annotate if any argument is supplied
if (len(sys.argv) > 1):
    print " *** Annotating ***"

    ### STABLE LABEL ###
    if len(sys.argv) == 3:
        # Arg 1
        x_crit = sys.argv[1]
        if isFloat(x_crit):
            x_crit = float(x_crit)
        else:
            x_crit = crit_sma(u_bin, e_bin, i_bin, crit_type = "hw99")
        # Arg 2
        x_st = sys.argv[2]
        if isFloat(x_st):
            x_st = float(x_st)
            #x_st = (float(T_st))**(2.0/3) # convert from T to a
        else:
            x_st = crit_sma(u_bin, e_bin, i_bin)

    else:
        x_crit = crit_sma(u_bin, e_bin, i_bin, crit_type = "hw99")
        x_st = crit_sma(u_bin, e_bin, i_bin)
    print "Lower Critical Sm-Axis:",  x_crit
    print "Outer Critical Sm-Axis:",  x_st

    # Set 'Stable' Text Height
    annotate_y = SIM_TIME / 3

    # Set 'Stable' Text Font
    font0 = FontProperties()
    font = font0.copy()
    font.set_weight("bold")
    font.set_family("Arial Narrow")
    ax1.annotate("    ALL  \nSTABLE", xy = (o.max_sma + 0.075, annotate_y), xytext = (x_st + 0.35, annotate_y),
                 arrowprops = dict(facecolor = "green", edgecolor = "green",
                                   arrowstyle = "simple, tail_width = 0.7, head_width = 1.6"),
                 verticalalignment = "center",
                 fontproperties = font, fontsize = 18, color = "green")

    #### A_ST LINE LABEL ####
    color_st = "green"
    # Plot 'A_ST' Line
    num_points = 50
    x_st_line = x_st + 0.00
    ax1.scatter(np.linspace(x_st_line, x_st_line, num_points), np.logspace(1, np.log10(2.4 * SIM_TIME), num_points), color = color_st, marker = 'o', s = 6, zorder = 5)

    # Set 'A_ST' Text Height
    annotate_y_top = 3000
    annotate_y_bottom = 1500

    # Set 'A_ST' Text Font
    font2 = font0.copy()
    font2.set_weight("bold")
    ax1.annotate(r"$a_{\rm crit} = a_{\rm st}$", xy = (x_st + 0.04, annotate_y_top), xytext = (x_st + 0.4, annotate_y_bottom),
                 arrowprops = dict(facecolor = color_st, edgecolor = color_st,
                                   arrowstyle = "simple, tail_width = 0.1, head_width = 0.7"),
                 verticalalignment = "center", 
                 fontproperties = font2, size = 25, color = color_st)

    #### A_CRIT LINE LABEL ####
    color_crit = "green"
    # Plot 'A_CRIT' Line
    x_crit_line = x_crit + 0.00
    ax1.plot([x_crit_line, x_crit_line], [10, 2.4 * SIM_TIME], color = color_crit, linewidth = 3, linestyle = "dashed", dashes = [4, 3], zorder = 5)

    # Set 'A_CRIT' Text Height
    annotate_y_top = 100
    annotate_y_bottom = 50

    # Set 'A_CRIT_HW' Text Font
    font3 = font0.copy()
    font3.set_weight("bold")
    """
    ax1.annotate("$a_{crit}$", xy = (x_crit_line + 0.03, annotate_y_top), xytext = (x_crit + 0.42, annotate_y_bottom),
                 arrowprops = dict(facecolor = color_crit, edgecolor = color_crit,
                                   arrowstyle = "simple, tail_width = 0.08, head_width = 0.6"),
                 verticalalignment = "center", 
                 fontproperties = font3, size = 20, color = color_crit)
    """

    #### A_CRIT_HW LINE LABEL ####
    color_hw = "black"
    # Plot 'A_CRIT_HW' Line
    x_crit_hw = holman_fit(u_bin, e_bin)
    ax1.plot([x_crit_hw, x_crit_hw], [10, 2.4 * SIM_TIME], color = color_hw, linewidth = 3, linestyle = "-.", dashes = [12, 6, 4, 6], zorder = 5)

    # Set 'A_CRIT_HW' Text Height
    annotate_y_top = 45
    annotate_y_bottom = 23

    # Set 'A_CRIT_HW' Text Font
    font3 = font0.copy()
    font3.set_weight("bold")
    ax1.annotate(r"$a_{\rm crit,HW99}$", xy = (x_crit_hw + 0.014, annotate_y_top), xytext = (x_crit_hw + 0.50, annotate_y_bottom), zorder = 6,
                 arrowprops = dict(facecolor = color_hw, edgecolor = color_hw,
                                   arrowstyle = "simple, tail_width = 0.08, head_width = 0.6"),
                 verticalalignment = "center", 
                 fontproperties = font3, size = 20, color = color_hw)



name = "t_vs_a_and_T_for_%s_u%02d_e%02d_i%02d" % (o.dir, 100 * u_bin, 100 * e_bin, i_bin)

plot.savefig(name + ".png", bbox_inches='tight')
plot.savefig(name + ".eps", format = "eps", dpi = 1000, bbox_inches='tight')
plot.savefig(name + ".pdf", format = "pdf", dpi = 1000, bbox_inches='tight')
plot.show()

plot.cla()

# Save Arrays
"""
# Save Arrays for Re-Use in Combined Plots
pickle.dump(x_r, open("xr" + name, 'wb'))
pickle.dump(y, open("y" + name, 'wb'))
"""
