"""
plot just data points (because the fits don't add much)
a_st + a_crit + a_crit,HW99
"""

import sys
import numpy as np

from matplotlib import gridspec
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
from matplotlib.ticker import AutoMinorLocator

import pickle

ecc = [0, 0.15, 0.6]

# Real Data
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


# Plot Parameters
fontsize = 16
labelsize = 14
rc['xtick.labelsize'] = labelsize
rc['ytick.labelsize'] = labelsize

min_y = 1.8
max_y = 4.5

### SUB PLOTS ###
f, (ax1, ax2, ax3) = plot.subplots(1, 3, sharey=True)
f.subplots_adjust(wspace = 0)
f.set_size_inches(15, 6, forward = True)
#f = plot.figure(figsize =(12, 5))
#gs = gridspec.GridSpec(1,3)
#ax1 = f.add_subplot(gs[0])
#ax2 = f.add_subplot(gs[1], sharey=ax1)
#ax3 = f.add_subplot(gs[2], sharey=ax1)
ax4 = ax1.twinx()
ax5 = ax2.twinx()
ax6 = ax3.twinx()

# Set up each figure
def plot_ax(ax, ax_b, mu):
    """
    Note: Plots are on twin axes (ax_b) because that is where the grid lines are
    """

    # Annotate
    ax.set_xlim(-5, 75)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylim(min_y, max_y)
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax2.set_xlabel("Inclination", fontsize = fontsize)
    ax1.set_ylabel("Semimajor Axis", fontsize = fontsize)
    ax.set_title("$\mu = %.2f$" % mu, fontsize = fontsize)

    ### Plot Real Data Points ###
    x_real = np.array([10 * x for x in range(8)])

    y_crit_real_2 = np.array([crit_sma(mu, ecc[2], inc_i, crit_type = "hw99") for inc_i in x_real])
    y_st_real_2 = np.array([crit_sma(mu, ecc[2], inc_i, crit_type = "safe") for inc_i in x_real])

    y_crit_real_1 = np.array([crit_sma(mu, ecc[1], inc_i, crit_type = "hw99") for inc_i in x_real])
    y_st_real_1 = np.array([crit_sma(mu, ecc[1], inc_i, crit_type = "safe") for inc_i in x_real])

    y_crit_real_0 = np.array([crit_sma(mu, ecc[0], inc_i, crit_type = "hw99") for inc_i in x_real])
    y_st_real_0 = np.array([crit_sma(mu, ecc[0], inc_i, crit_type = "safe") for inc_i in x_real])

    # Real Data
    zorder_real = 10
    size_real = 8
    linewidth = 3
    ax_b.plot(x_real, y_st_real_2, color = "#660000", linewidth = linewidth, marker = "^", ms = size_real, zorder = zorder_real)
    ax_b.plot(x_real, y_crit_real_2, color = "#00004C", linewidth = linewidth, marker = "v", ms = size_real, zorder = zorder_real)

    ax_b.plot(x_real, y_st_real_1, color = "#FF0000", linewidth = linewidth, marker = "^", ms = size_real, zorder = zorder_real)
    ax_b.plot(x_real, y_crit_real_1, color = "#0000FF", linewidth = linewidth, marker = "v", ms = size_real, zorder = zorder_real)

    ax_b.plot(x_real, y_st_real_0, color = "#FF5050", linewidth = linewidth, marker = "^", ms = size_real, zorder = zorder_real)
    ax_b.plot(x_real, y_crit_real_0, color = "#33CCFF", linewidth = linewidth, marker = "v", ms = size_real, zorder = zorder_real)


    # Legend
    #l = ax4.legend(bbox_to_anchor = (0.35, 1))
    #l.get_frame().set_alpha(1.0)
    #l.set_zorder(20)
    #ax6.legend(bbox_to_anchor = (1, 0.40))

# Select 3 different values of mu
plot_ax(ax1, ax4, 0.1)
plot_ax(ax2, ax5, 0.3)
plot_ax(ax3, ax6, 0.5)

# Set up Orbital Period on y-axis of ax3
min_r = int(np.ceil((min_y - 0.1)**1.5))
max_r = int(np.ceil((max_y + 0.1)**1.5))
new_tick_locations = [q**(0.66667) for q in np.array(range(min_r, max_r))]
def resonant_axis_labels(a_s):
    """ set up axis label numbers """
    r_s = [a**1.5 for a in a_s]
    r_labels = ["%d" % int(round(r,0)) for r in r_s]
    return r_labels

ax4.set_ylim(ax1.get_ylim())
ax4.set_yticks(new_tick_locations)
ax4.set_yticklabels(resonant_axis_labels(new_tick_locations), visible = False)
ax4.grid(b = True, which = "major", color = "black", linestyle = "--")

ax5.set_ylim(ax2.get_ylim())
ax5.set_yticks(new_tick_locations)
ax5.set_yticklabels(resonant_axis_labels(new_tick_locations), visible = False)
ax5.grid(b = True, which = "major", color = "black", linestyle = "--")

ax6.set_axisbelow(True)
ax6.set_ylim(ax3.get_ylim())
ax6.set_yticks(new_tick_locations)
ax6.set_yticklabels(resonant_axis_labels(new_tick_locations))
ax6.set_ylabel("Orbital Period [$T_b$]   ", fontsize = fontsize, rotation = 270)
ax6.grid(b = True, which = "major", color = "black", linestyle = "--", zorder = 0)

# Annotate
ax1.text(0, 2.3, "$e=0.0$", fontsize = fontsize)
ax1.text(30, 3.1, "$e=0.15$", fontsize = fontsize)
ax1.text(0, 4.02, "$e=0.6$", fontsize = fontsize)

ax2.text(0, 2.38, "$e=0.0$", fontsize = fontsize)
ax2.text(30, 3.1, "$e=0.15$", fontsize = fontsize)
ax2.text(0, 4.02, "$e=0.6$", fontsize = fontsize)

ax3.text(0, 2.35, "$e=0.0$", fontsize = fontsize)
ax3.text(40, 3.1, "$e=0.15$", fontsize = fontsize)
ax3.text(0, 4.02, "$e=0.6$", fontsize = fontsize)

# Save figure
save_fn = "e_plots/fixed_e_comparison.%s"
plot.savefig(save_fn % ("png"), bbox_inches='tight')
plot.savefig(save_fn % ("pdf"), format = "pdf", dpi = 1000, bbox_inches='tight')
plot.show()

