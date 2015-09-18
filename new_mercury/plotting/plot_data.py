"""
plot fits with data points
a_st + a_crit + a_crit,HW99
"""

import sys
import numpy as np

from matplotlib import gridspec
from matplotlib import rcParams as rc
from matplotlib import pyplot as plot
from matplotlib.ticker import AutoMinorLocator

import pickle

inc = int(sys.argv[1])

info_fit_crit = pickle.load(open("fit_crit_inc%02d.p" % inc, "rb"))
info_fit_st = pickle.load(open("fit_st_inc%02d.p" % inc, "rb"))

popt_crit = info_fit_crit["popt"]
popt_st = info_fit_st["popt"]


# Fit Functions
def fit(x, const, e_c, esq_c, m_c, me_c, msq_c, esqmsq_c):
    m = x[0]
    e = x[1]
    return const + e_c * (e) + esq_c * (e*e) + m_c * (m) + me_c * (m*e) + msq_c * (m*m) + esqmsq_c * (m*m*e*e)

def fit_crit(mu, e):
    return fit([mu, e], popt_crit[0], popt_crit[1], popt_crit[2], popt_crit[3], popt_crit[4], popt_crit[5], popt_crit[6])

def fit_st(mu, e):
    return fit([mu, e], popt_st[0], popt_st[1], popt_st[2], popt_st[3], popt_st[4], popt_st[5], popt_st[6])

def holman_fit(mass_ratio, ecc):
    """HW99 fitting formula for the LCO"""
    a_crit = 1.60 + 5.10 * ecc - 2.22 * ecc**2 \
             + 4.12 * mass_ratio - 4.27 * mass_ratio * ecc \
             -5.09 * mass_ratio**2 + 4.61 * mass_ratio**2 * ecc**2

    return a_crit


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

min_x = -0.05
max_x = 0.75

min_y = 1.8
max_y = 7.2 # 7.2, 6.7

### SUB PLOTS ###
#f, ([ax1, ax2, ax3], [ax7, ax8, ax9]) = plot.subplots(nrows = 2, ncols = 3, sharey="row", gridspec_kw = dict(height_ratios=[0.35, 0.05, 1, 1, 1, 1]))
#f.subplots_adjust(wspace = 0)
#f.set_size_inches(14, 4, forward = True)
f = plot.figure(figsize =(14, 5.5))
gs = gridspec.GridSpec(2,3, height_ratios = [3,1], wspace = 0, hspace = 0.08)

ax1 = f.add_subplot(gs[0])
ax2 = f.add_subplot(gs[1], sharey=ax1)
ax3 = f.add_subplot(gs[2], sharey=ax1)

ax7 = f.add_subplot(gs[3])
ax8 = f.add_subplot(gs[4], sharey=ax7)
ax9 = f.add_subplot(gs[5], sharey=ax7)

ax4 = ax1.twinx()
ax5 = ax2.twinx()
ax6 = ax3.twinx()

# Set up each figure (with data)
def plot_ax(ax, ax_b, mu):
    """
    Note: Plots are on twin axes (ax_b) because that is where the grid lines are
    """

    # Annotate
    ax.set_xlim(-0.05, 0.75)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xticklabels("", visible = False)
    ax.set_ylim(min_y, max_y)
    min_y_tick = round(2 * min_y, 0) / 2.0 # Round to nearest half-integer
    max_y_tick = round(2 * max_y, 0) / 2.0
    if ax is ax1:
        ax1.set_yticks(np.linspace(min_y_tick, max_y_tick, 2.0 * (max_y_tick - min_y_tick) + 1)) # e.g. 2.0, 2.5, ..., 6.5, 7.0
    else:
        ax.set_yticklabels(np.linspace(min_y_tick, max_y_tick, 2.0 * (max_y_tick - min_y_tick) + 1), visible = False) # e.g. 2.0, 2.5, ..., 6.5, 7.0
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    #ax2.set_xlabel("Eccentricity", fontsize = fontsize)
    ax1.set_ylabel("Semimajor Axis", fontsize = fontsize)
    ax.set_title("( $\mu = %.2f$, $i = %d^{\circ}$)" % (mu, inc), fontsize = fontsize)

    label_hw = ""
    label_crit = ""
    label_st = ""
    label_crit_dots = ""
    label_st_dots = ""
    if ax_b == ax4:
        label_crit_dots = "$a_{crit}$"
        label_st_dots = "$a_{st}$"
    elif ax_b == ax6:
        label_hw = "$a_{crit,HW99}$"
        label_crit = "$a_{crit}$"
        label_st = "$a_{st}$"

    ### Plot Fits ###
    x_fit = np.linspace(0, 0.7, 100) # eccentricities
    y_crit_hw = np.array([holman_fit(mu, e_i) for e_i in x_fit])
    y_crit_fit = np.array([fit_crit(mu, e_i) for e_i in x_fit])
    y_st_fit = np.array([fit_st(mu, e_i) for e_i in x_fit])

    # Fitting Functions
    zorder_fit = 1
    linewidth = 4
    #ax_b.plot(x_fit, y_st_fit, linewidth = linewidth, color = "orange", zorder = zorder_fit, label = label_st)
    #ax_b.plot(x_fit, y_crit_fit, linewidth = linewidth, color = "#64C5ED", zorder = zorder_fit, label = label_crit)
    #ax_b.plot(x_fit, y_crit_hw, linewidth = linewidth - 1, color = "black", linestyle = "-.", dashes = [12, 6, 4, 6], zorder = zorder_fit, label = label_hw)
    
    ### Plot Real Data Points ###
    x_real = np.array([0.05 * x for x in range(15)])
    y_crit_real = np.array([crit_sma(mu, e_i, inc, crit_type = "hw99") for e_i in x_real])
    y_st_real = np.array([crit_sma(mu, e_i, inc, crit_type = "safe") for e_i in x_real])

    # Real Data
    zorder_real = 10
    size_real = 8
    linewidth = 3
    ax_b.plot(x_real, y_st_real, color = "red", linewidth = linewidth, marker = "^", ms = size_real, zorder = zorder_real, label = label_st_dots)
    ax_b.plot(x_real, y_crit_real, color = "blue", linewidth = linewidth, marker = "v", ms = size_real, zorder = zorder_real, label = label_crit_dots)

    # Legend
    l = ax4.legend(bbox_to_anchor = (0.35, 1))
    #l.get_frame().set_alpha(1.0)
    #l.set_zorder(20)
    ax6.legend(bbox_to_anchor = (1, 0.40))

# Select 3 different values of mu
mu_1 = 0.1; mu_2 = 0.3; mu_3 = 0.5
plot_ax(ax1, ax4, mu_1)
plot_ax(ax2, ax5, mu_2)
plot_ax(ax3, ax6, mu_3)

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
ax4.grid(b = True, which = "major", color = "black", linestyle = "--", dashes = [1, 2])

ax5.set_ylim(ax2.get_ylim())
ax5.set_yticks(new_tick_locations)
ax5.set_yticklabels(resonant_axis_labels(new_tick_locations), visible = False)
ax5.grid(b = True, which = "major", color = "black", linestyle = "--", dashes = [1, 2])

ax6.set_axisbelow(True)
ax6.set_ylim(ax3.get_ylim())
ax6.set_yticks(new_tick_locations)
ax6.set_yticklabels(resonant_axis_labels(new_tick_locations))
ax6.set_ylabel("Orbital Period [$T_b$]   ", fontsize = fontsize, rotation = 270)
ax6.grid(b = True, which = "major", color = "black", linestyle = "--", dashes = [1, 2], zorder = 0)


###### Set up each figure (with delta data) ######
min_y = -0.6
max_y = 1.24

def plot_ax2(ax, mu):
    ax.set_xlim(min_x, max_x)
    #ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.set_ylim(min_y, max_y)
    ax.yaxis.grid(b = True, which = "major", color = "black", linestyle = "--")
    min_y_tick = round(2 * min_y, 0) / 2.0 # Round to nearest half-integer
    max_y_tick = round(2 * max_y, 0) / 2.0
    if ax is ax7:
        ax7.set_yticks(np.linspace(min_y_tick, max_y_tick, 2.0 * (max_y_tick - min_y_tick) + 1)) # e.g. 0.0, 0.5, 1.0
    else:
        ax.set_yticklabels(np.linspace(min_y_tick, max_y_tick, 2.0 * (max_y_tick - min_y_tick) + 1), visible = False) # e.g. 2.0, 2.5, ..., 6.5, 7.0

    ax8.set_xlabel("Eccentricity", fontsize = fontsize)
    ax7.set_ylabel("$a_{%d^{\circ}} -$ $a_{0^{\circ}}$" % inc, fontsize = fontsize)

    ### Plot Real Data Points ###
    x_real = np.array([0.05 * x for x in range(15)])
    dy_crit_real = np.array([crit_sma(mu, e_i, inc, crit_type = "hw99") - crit_sma(mu, e_i, 0, crit_type = "hw99") for e_i in x_real])
    dy_st_real = np.array([crit_sma(mu, e_i, inc, crit_type = "safe") - crit_sma(mu, e_i, 0, crit_type = "safe") for e_i in x_real])

    zorder_real = 3
    size_real = 6
    linewidth = 2
    ax.plot(x_real, dy_st_real, color = "red", linewidth = linewidth, marker = "^", ms = size_real - 2, label = "", zorder = zorder_real)
    ax.plot(x_real, dy_crit_real, color = "blue", linewidth = linewidth, marker = "v", ms = size_real - 2, label = "", zorder = zorder_real)

    # Plot Horizontal Line
    ax.plot([min_x - 1, max_x + 1], [0, 0], color = "black", linewidth = 1, zorder = 1)

# Plot each delta_a
plot_ax2(ax7, mu_1)
plot_ax2(ax8, mu_2)
plot_ax2(ax9, mu_3)

# Save figure
save_fn = "data_plots/data_comparison_inc%02d.%s"
plot.savefig(save_fn % (inc, "png"), bbox_inches='tight')
plot.savefig(save_fn % (inc, "pdf"), format = "pdf", dpi = 1000, bbox_inches='tight')
plot.show()

