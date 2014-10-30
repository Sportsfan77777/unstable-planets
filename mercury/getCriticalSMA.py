"""
gathers the critical sma over a given range of ((u, e), i)
and makes two tables representing the data
"""

import numpy as np
import pickle
import sys

dir_base = "sim"
mass_ratios = [0.1 * x for x in range(1,6)]
ecc_s = [0.1 * x for x in range(8)]
inc_s = [10.0 * x for x in range(1,10)]
Ms = [0, 180]

# Make 2-D table at each i
# Make 1-D table for each (u,e)

critical_sma_table = np.zeros((len(ecc_s), len(mass_ratios)), dtype = float)
critical_sma = np.zeros(len(ecc_s) * len(mass_ratios), dtype = float)

for inc in inc_s:
    for j,ecc in enumerate(ecc_s):
        for k,u in enumerate(mass_ratios):
            sm_axes = []
            stable_arrays = [[],[]]
            for i,M in enumerate(Ms):
                u_bin_str = int(round(u * 10, 0))
                e_bin_str = int(round(ecc * 10, 0))
                i_bin_str = int(round(inc, 0))
                M_bin_str = int(round(M, 0))
                stable_file = "%s_u%d_e%d_i%03d_M%03d/stability.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str)
                stable_file = "%s_u%d_e%d_i%03d_M%03d/sm_axes.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str) 
                
                stability_f = open(stable_file, "rb")
                stable_arrays[i] = stability_f.load()
                stability_f.close()
                
            sm_axes_f = open(stable_file, "rb")
            sm_axes = sm_axes_f.load()
            sm_axes_f.close()
            
            crit_sma = sm_axes[-1] + 900.0
            for sma, s0, s1 in zip(sm_axes, stable_arrays[0], stable_arrays[1]):
                 if s0 and s1:
                     crit_sma = sma
                     break
            critical_sma_table[j,k] = crit_sma
            critical_sma[j* len(mass_ratios) + k] = crit_sma
            
# Pickle Files!

pickle_f = open(("%s_critical_sma_table.p", "wb") % (dir_base))
pickle.dump(critical_sma_table, pickle_f)
pickle_f.close()

pickle_f = open(("%s_critical_sma_list.p", "wb") % (dir_base))
pickle.dump(critical_sma_list, pickle_f)
pickle_f.close()
            
# Pretty Print Table and List

table_fn = "%s_critical_sma_table.txt" % (dir_base)
list_fn = "%s_critical_sma_list.txt" % (dir_base)

table_f = open(table_fn, "w")
list_f = open(list_fn, "w")

# (1) Table

# (2) List



table_f.close()
list_f.close()
