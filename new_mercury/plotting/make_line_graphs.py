"""
parse giant tables and lists
"""

# SPLIT INTO (i) CONTOUR PLOT AND (ii) FIT --- due to lack of scipy

import sys
import pickle
import math
import numpy as np

from make_plots_ecc import line_graph

dir_base = "sim"
if len(sys.argv) > 1:
    dir_base = sys.argv[0]

max_mass_i = 11
num_ecc = 5

mass_ratios = [0.05 * x for x in range(2, max_mass_i)]
ecc_s = [0.05 * x for x in range(num_ecc)]
inc_s = [10.0 * x for x in range(1,8)]
Ms = [0, 180]

tables_fn = "%s_safe_critical_sma_inc_tables.p" % dir_base
list_fn = "%s_safe_critical_sma_inc_list.p" % dir_base

tables_f = open(tables_fn, "rb")
critical_sma_tables = pickle.load(tables_f)
tables_f.close()

list_f = open(list_fn, "rb")
critical_sma_list = pickle.load(list_f)
list_f.close()

#print "TABLES"
#print critical_sma_tables # Note: Table is more useful (i, e, u)

#print "LIST"
#print critical_sma_list

print np.shape(critical_sma_tables)
print np.shape(critical_sma_list)

master_list = np.zeros(len(inc_s) * len(ecc_s) * len(mass_ratios)) # y_data
parameter_list = np.zeros((3, len(inc_s) * len(ecc_s) * len(mass_ratios))) # x_data
# Make (10 * 8 * 5)-sized master_list of data
# Make (3 * (10 * 8 * 5))-sized list of parameters

print np.shape(parameter_list)

count = 0
for i, inc in enumerate(inc_s):
    for j, ecc in enumerate(ecc_s):
        for k, mass_ratio in enumerate(mass_ratios):
            parameter_list[0, count] = math.cos(inc)
            parameter_list[1, count] = ecc
            parameter_list[2, count] = mass_ratio

            master_list[count] = critical_sma_tables[i, j, k]

            count += 1

for i in range(count):
    inc = parameter_list[0, i]
    ecc = parameter_list[1, i]
    mass_ratio = parameter_list[2, i]
    crit_sma = master_list[i]
    #print "Ecc:, %.1f, Mass Ratio: %.1f, Critical SM-Axis: %1.1f" % (ecc, mass_ratio, crit_sma)


####### MAKE A CONTOUR PLOT AND 3-D PLOT ########

for j, ecc in enumerate(ecc_s):
    line_graph(critical_sma_tables[:,j,:], ecc)


