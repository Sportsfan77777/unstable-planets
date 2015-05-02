"""
Make plots of critical sm-axes for various conditions
"""

from matplotlib import pyplot as plot
from matplotlib import ticker
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import pickle

max_mass_i = 11
num_ecc = 5

mass_ratios = [0.05 * x for x in range(2, max_mass_i)]
ecc_s = [0.05 * x for x in range(num_ecc)]
inc_s = [10.0 * x for x in range(1,8)]
Ms = [0, 180]

degree_sign= u'\N{DEGREE SIGN}'

def contour_plot(cr_sma_table, inc, mass_ratios = mass_ratios, ecc_s = ecc_s):
    """ make a contour plot of table (use globals for axes) """

    levels = range(200, 461, 30)
    levels = [0.01 * l for l in levels]

    fig, ax = plot.subplots()

    #cmap = plot.get_cmap('YlGnBu_r')
    cmap = plot.get_cmap('RdYlBu')

    #contour_cmap = plot.get_cmap('YlGnBu')
    #contour_cmap = plot.get_cmap('winter_r')
    contour_cmap = plot.get_cmap('hot_r')

    e_bins = np.array([0.05 * x - 0.025 for x in range(num_ecc + 1)])
    m_bins = np.array([0.05 * x + 0.025 for x in range(max_mass_i)])

    p_map = ax.pcolormesh(m_bins, e_bins, cr_sma_table, cmap = cmap)
    p_con = ax.contour(mass_ratios, ecc_s, cr_sma_table, levels = levels, linewidths = 2, cmap = contour_cmap)

    ax.clabel(p_con, fmt = "%.1f")

    plot.gca().invert_yaxis()

    #ax.set_yscale('log')

    #plot.title("Transferred Particle Counts (by %)")
    plot.xlabel("Mass Ratio u", fontsize = 15)
    plot.ylabel("Eccentricity e", fontsize = 15)

    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.get_yaxis().get_major_formatter().labelOnlyBase = False

    ax.set_xticks(mass_ratios)
    ax.set_yticks(ecc_s)

    p_map.set_clim([levels[0], levels[-1]]) # <<---- USE THIS TO STANDARDIZE COLORBAR
    cbar = plot.colorbar(p_map)
    #cbar = plot.colorbar(p_map, ticks = levels, norm = plot.Normalize(vmin=0, vmax=1))
    cbar.set_label("Critical SMA (in AU) at Inclination %02d%s" % (inc, degree_sign), fontsize = 15)

    plot.tight_layout()

    #plot.show()
    plot.savefig("critical_sma_cmap_inc%02d.png" % inc)

    plot.clf()

    ########## MAKE A 3-D PLOT!!!!!! ##########

    fig = plot.figure()
    ax = Axes3D(fig)

    ax.set_xticks(mass_ratios)
    ax.set_yticks(ecc_s)

    X, Y = np.meshgrid(mass_ratios, ecc_s)
    Z = cr_sma_table

    surface = ax.plot_surface(X, Y, Z, rstride=1, cstride = 1, cmap = cmap, linewidth = 0, antialiased = False)

    #ax.set_yscale('log')

    #plot.title("Transferred Particle Counts (by %)")
    plot.xlabel("Mass Ratio u", fontsize = 15)
    plot.ylabel("Eccentricity e", fontsize = 15)

    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.get_yaxis().get_major_formatter().labelOnlyBase = False

    ax.set_xticks(mass_ratios)
    ax.set_yticks(ecc_s)

    cbar = fig.colorbar(surface)
    cbar.set_label("Critical SMA (in AU) at Inclination %02d%s" % (inc, degree_sign), fontsize = 15)

    #plot.colorbar(ticks = levels, norm = plot.Normalize(vmin=0, vmax=1))

    #plot.show()

    plot.clf()


def line_graph():
    pass