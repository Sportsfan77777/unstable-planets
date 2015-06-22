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
inc_s = [10.0 * x for x in range(1)]
Ms = [0, 180]

degree_sign= u'\N{DEGREE SIGN}'

## LINE GRAPH FEATURES
colors = pickle.load(open("/Users/Sportsfan77777/planets/plots/features/colors.p", "rb"))
markers = pickle.load(open("/Users/Sportsfan77777/planets/plots/features/markers.p", "rb"))


def fit(mass_ratio, ecc_s):
    def holman_fit(mass_ratio, ecc):
        """HW99 fitting formula for the LCO"""
        a_crit = 1.60 + 5.10 * ecc - 2.22 * ecc**2 \
                 + 4.12 * mass_ratio - 4.27 * mass_ratio * ecc \
                 -5.09 * mass_ratio**2 + 4.61 * mass_ratio**2 * ecc**2

        return a_crit

    return np.array([holman_fit(mass_ratio, e) for e in ecc_s])


def line_graph(cr_sma_table, inc, mass_ratios = mass_ratios, inc_s = inc_s, apo = False):
     """ make a set of line graphs from the table (use globals for axes) """

     # Something to do with legend placement

     fontsize = 17

     #title = "Inclination $i$ = %2d%s" % (inc, degree_sign)
     title = "Comparison of UCO ($a_{st}$) to HW99 Fit of LCO ($a_{crit}$)"
     #title = "Comparison of our LCO to HW99 Fit of LCO ($a_{crit}$)"
     plot.title(title, fontsize = fontsize)

     plot.xlabel("Eccentricity $e$", fontsize = fontsize)
     plot.ylabel("$a_{st}$ [$a_b$]", fontsize = fontsize)

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

     for k, mass_ratio in enumerate(mass_ratios):
        label = "$\mu$ = %0.02f" % mass_ratio

        # Scatter points (do not connect the dots)
        y = [yk + 0.002 * (k - len(mass_ratios) / 2) for yk in cr_sma_table[:,k]]
        this_scatter, = plot.plot(x, y, c = colors[mass_ratio], marker = markers[mass_ratio], linewidth = 0, markersize = 7)

        # Draw Holman Fit
        y_hol = fit(mass_ratio, x)
        print y_hol
        this_plot, = plot.plot(x, y_hol, c = colors[mass_ratio], marker = markers[mass_ratio], linewidth = 2, markersize = 0)

        #pickle.dump(y, open("y_i%02d_u%02d.p" % (inc, 100 * mass_ratio), "wb")) ## <<<==== Only need to do this once

        plots.append(this_plot)
        labels.append(label)

     plot.legend(plots, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

     #plot.show()
     fn = "line_graph_holman_safe_inc%02d" % (inc) 
     if apo:
        fn = "line_graph_holman_div_by_Q_inc%02d" % (inc) 
     plot.savefig(fn + ".png", bbox_inches='tight')
     plot.savefig(fn + ".pdf", format = "pdf", dpi = 1000, bbox_inches='tight')

     # Clear Figure
     plot.clf()



