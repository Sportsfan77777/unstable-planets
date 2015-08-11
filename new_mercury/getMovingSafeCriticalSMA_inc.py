"""
gathers the critical sma (a_st) over a given range of ((u, e), i)
and makes two tables representing the data

Definiton: UCO, with median 'a' instead of initial 'a'
"""

import numpy as np
import pickle
import sys

STABLE_VALUE = 999888.9

dir_base = "sim"
mass_ratios = [0.05 * x for x in range(2,11)]
ecc_s = [0.05 * x for x in range(15)]
inc_s = [10.0 * x for x in range(8)]
Ms = [0, 180]

# Make 2-D table at each i
# Make 1-D table for each (u,e)

critical_sma_tables = np.zeros((len(inc_s), len(ecc_s), len(mass_ratios)), dtype = float)
critical_sma_list = np.zeros((len(inc_s), (len(ecc_s) * len(mass_ratios))), dtype = float)

for ith, inc in enumerate(inc_s):
    for j,ecc in enumerate(ecc_s):
        for k,u in enumerate(mass_ratios):
            eject_table_arrays = [[],[]]
            sma_table_arrays = [[],[]]
            for i,M in enumerate(Ms):
                u_bin_str = int(round(u * 100, 0))
                e_bin_str = int(round(ecc * 100, 0))
                i_bin_str = int(round(inc, 0))
                M_bin_str = int(round(M, 0))
                eject_table_fn = "storage/%s_u%02d_e%02d_i%03d_M%03d/table_of_ejections.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str)
                sma_table_fn = "storage/%s_u%02d_e%02d_i%03d_M%03d/table_of_final_sm-axes.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str)
                
                eject_table_f = open(eject_table_fn, "rb")
                eject_table_arrays[i] = pickle.load(eject_table_f)
                eject_table_f.close()

                sma_table_f = open(sma_table_fn, "rb")
                sma_table_arrays[i] = pickle.load(sma_table_f)
                sma_table_f.close()
            
            crit_sma = 0.0
            for (e0, e0_a) in zip(eject_table_arrays[0], sma_table_arrays[0]):
                for (eject_time, median_sma) in zip(e0, e0_a):
                    if eject_time == STABLE_VALUE:
                        pass
                    else:
                        # Find the ejected planet with the largest median sm-axis (search through M = 0)
                        if median_sma > crit_sma:
                            crit_sma = round(median_sma, 3)

            for (e1, e1_a) in zip(eject_table_arrays[1], sma_table_arrays[1]):
                for (eject_time, median_sma) in zip(e0, e0_a):
                    if eject_time == STABLE_VALUE:
                        pass
                    else:
                        # Find the ejected planet with the largest median sm-axis (search through M = 180)
                        if median_sma > crit_sma:
                            crit_sma = round(median_sma, 3)

            critical_sma_tables[ith,j,k] = crit_sma
            critical_sma_list[ith, j* len(mass_ratios) + k] = crit_sma
            
# Pickle Files!

pickle_f = open(("%s_moving_safe_critical_sma_inc_tables.p" % (dir_base)), "wb")
pickle.dump(critical_sma_tables, pickle_f)
pickle_f.close()

pickle_f = open(("%s_moving_safe_critical_sma_inc_list.p" % (dir_base)), "wb")
pickle.dump(critical_sma_list, pickle_f)
pickle_f.close()
            
# Pretty Print Table and List

tables_fn = "%s_moving_safe_critical_sma_inc_tables.txt" % (dir_base)
list_fn = "%s_moving_safe_critical_sma_inc_list.txt" % (dir_base)

tables_f = open(tables_fn, "w")
list_f = open(list_fn, "w")


############### SORT BY INCLINATION #################

# (1) Tables

width = 8
dash = '-' * width

for ith, inc in enumerate(inc_s):
    row_zero = ("Inc%02d" % inc).center(width) + '|'
    dash_row = dash + '-'
    for u in mass_ratios:
        row_zero += (str(u)).center(width)
        dash_row += dash
        
    rows = []
    for j, ecc in enumerate(ecc_s):
        row = (str(ecc)).center(width) + '|'
        row_eject = (str(ecc)).center(width) + '|'
        for k, u in enumerate(mass_ratios):
            row += ("%.03f" % (critical_sma_tables[ith,j,k])).center(width)
        rows.append(row)

    # File 1: the critical 'a'
        
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
  for j, ecc in enumerate(ecc_s):
    row = ("(u=%.1f, e=%.1f)" % (u, ecc)).center(2*width) + '|'
    for ith, inc in enumerate(inc_s):
        row += str(critical_sma_list[ith, j*len(mass_ratios) + k]).center(width)
    rows.append(row)
  rows.append(dash_row)

# File 1: the critical 'a'
        
list_f.write(row_zero + "\n")
list_f.write(dash_row + "\n")
for r_str in rows:
   list_f.write(r_str + "\n")
list_f.write("\n")

# Close both files

tables_f.close()
list_f.close()


