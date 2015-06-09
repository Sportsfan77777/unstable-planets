"""
Make plots of critical sm-axes for various conditions
"""

from matplotlib import pyplot as plot
from matplotlib import ticker
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import pickle

max_mass_i = 6
num_ecc = 15

mass_ratios = [round(0.10 * x, 1) for x in range(1, max_mass_i)]
ecc_s = [0.05 * x for x in range(num_ecc)]
inc_s = [10.0 * x for x in range(8)]
Ms = [0, 180]

degree_sign= u'\N{DEGREE SIGN}'

## LINE GRAPH FEATURES
colors = pickle.load(open("/Users/Sportsfan77777/planets/plots/features/colors.p", "rb"))
markers = pickle.load(open("/Users/Sportsfan77777/planets/plots/features/markers.p", "rb"))


def line_graph(cr_sma_table, inc, mass_ratios = mass_ratios, inc_s = inc_s, apo = False):
     """ make a set of line graphs from the table (use globals for axes) """

     # Something to do with legend placement

     title = "Delta Critical Semimajor Axes at Inclination $i$ = %2d%s" % (inc, degree_sign)
     plot.title(title)

     plot.xlabel("Eccentricity $e$")
     plot.ylabel("Delta Critical Semimajor Axis $a_{st}$ / $a_b$")

     x = ecc_s[:]
     plot.xlim(ecc_s[0] - 0.025, ecc_s[-1] + 0.025)

     # Divide by Apocenter
     if apo:
        plot.ylabel("Critical Semimajor Axis a_crit / (1 + e)\n(in units of binary separation) ")
        Q_s = 1 + np.array(ecc_s)
        cr_sma_table = np.transpose(np.transpose(cr_sma_table) / Q_s)

     cr_sma_list = cr_sma_table.reshape(np.size(cr_sma_table))

     plot.ylim(min(cr_sma_list) - 0.1, max(cr_sma_list) + 0.1)

     plots = []
     labels = []

     min_y = 0.0
     max_y = 0.2

     for k, mass_ratio in enumerate(mass_ratios):
        label = "$\mu$ = %0.02f" % mass_ratio

        y0 = pickle.load(open("y_i00_u%02d.p" % (100 * mass_ratio), "rb")) 
        y = np.array([yk + 0.002 * (k - len(mass_ratios) / 2) for yk in cr_sma_table[:,k]])

        # Update Bounds
        if min(y - y0) - 0.1 < min_y:
            min_y = min(y - y0) - 0.1
        if max(y - y0) + 0.1 > max_y:
            max_y = max(y - y0) + 0.1

        plot.ylim(min_y, max_y)

        this_plot, = plot.plot(x, y - y0, c = colors[mass_ratio], marker = markers[mass_ratio], linestyle = "dashed", linewidth = 1, markersize = 7)

        plots.append(this_plot)
        labels.append(label)

     plot.legend(plots, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

     #plot.show()
     fn = "line_graph_delta_inc%02d.png" % (inc) 
     if apo:
        fn = "line_graph_div_by_Q_delta_inc%02d.png" % (inc) 
     plot.savefig(fn, bbox_inches='tight')

     # Clear Figure
     plot.clf()



