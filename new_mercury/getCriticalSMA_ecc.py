"""
gathers the critical sma over a given range of ((u, e), i)
and makes two tables representing the data
"""

import numpy as np
import pickle
import sys

dir_base = "sim"
mass_ratios = [0.1 * x for x in range(1,6)]
ecc_s = [0.05 * x for x in range(15)]
inc_s = [10.0 * x for x in range(1,8)]
Ms = [0, 0]

# Make 2-D table at each i
# Make 1-D table for each (u,e)

critical_sma_tables = np.zeros((len(inc_s), len(ecc_s), len(mass_ratios)), dtype = float)
critical_sma_list = np.zeros((len(inc_s), (len(ecc_s) * len(mass_ratios))), dtype = float)

for ith, inc in enumerate(inc_s):
    for j,ecc in enumerate(ecc_s):
        for k,u in enumerate(mass_ratios):
            sm_axes = []
            stable_arrays = [[],[]]
            for i,M in enumerate(Ms):
                u_bin_str = int(round(u * 100, 0))
                e_bin_str = int(round(ecc * 100, 0))
                i_bin_str = int(round(inc, 0))
                M_bin_str = int(round(M, 0))
                stable_file = "storage/%s_u%02d_e%02d_i%03d_M%03d/stability.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str)
                sm_axes_fn = "storage/%s_u%02d_e%02d_i%03d_M%03d/sm_axes.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str) 
                
                stability_f = open(stable_file, "rb")
                stable_arrays[i] = pickle.load(stability_f)
                stability_f.close()
                
            sm_axes_f = open(sm_axes_fn, "rb")
            sm_axes = pickle.load(sm_axes_f)
            sm_axes_f.close()
            
            crit_sma = sm_axes[-1] + 900.0
            for sma, s0, s1 in zip(sm_axes, stable_arrays[0], stable_arrays[1]):
                 if s0 and s1:
                     crit_sma = sma
                     break
            critical_sma_tables[ith,j,k] = crit_sma
            critical_sma_list[ith, j* len(mass_ratios) + k] = crit_sma
            
# Pickle Files!

pickle_f = open(("%s_critical_sma_ecc_tables.p" % (dir_base)), "wb")
pickle.dump(critical_sma_tables, pickle_f)
pickle_f.close()

pickle_f = open(("%s_critical_sma_ecc_list.p" % (dir_base)), "wb")
pickle.dump(critical_sma_list, pickle_f)
pickle_f.close()
            
# Pretty Print Table and List

tables_fn = "%s_critical_sma_ecc_tables.txt" % (dir_base)
list_fn = "%s_critical_sma_ecc_list.txt" % (dir_base)

tables_f = open(tables_fn, "w")
list_f = open(list_fn, "w")

############### SORT BY ECCENTRICITY #################

# (1) Tables

width = 8
dash = '-' * width

for j, ecc in enumerate(ecc_s):
    row_zero = ("Ecc=%0.2f" % ecc).center(width) + '|'
    dash_row = dash + '-'
    for u in mass_ratios:
        row_zero += (str(u)).center(width)
        dash_row += dash
        
    rows = []
    for ith, inc in enumerate(inc_s):
        row = (str(inc)).center(width) + '|'
        for k, u in enumerate(mass_ratios):
            row += (str(critical_sma_tables[ith,j,k])).center(width)
        rows.append(row)
        
    tables_f.write(row_zero + "\n")
    tables_f.write(dash_row + "\n")
    for r_str in rows:
       tables_f.write(r_str + "\n")
    tables_f.write("\n")
    
# (2) List

row_zero = "(u=u0, e=e0)".center(2*width) + '|'
dash_row = dash + dash + '-'
for inc in inc_s:
    row_zero += (str(inc)).center(width)
    dash_row += dash
    
rows = []
for k, u in enumerate(mass_ratios):
  for ith, inc in enumerate(inc_s):
    row = ("(u=%.1f, i=%02d)" % (u, inc)).center(2*width) + '|'
    for j, ecc in enumerate(ecc_s):
        row += str(critical_sma_list[ith, j*len(mass_ratios) + k]).center(width)
    rows.append(row)
  rows.append(dash_row)
        
list_f.write(row_zero + "\n")
list_f.write(dash_row + "\n")
for r_str in rows:
   list_f.write(r_str + "\n")
list_f.write("\n")

tables_f.close()
list_f.close()

