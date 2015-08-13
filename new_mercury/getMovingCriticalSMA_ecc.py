"""
gathers the critical sma (a_st) over a given range of ((u, e), i)
and makes two tables representing the data

Definition: Refined LCO (better than HW99)
"""

import numpy as np
import pickle
import sys

STABLE_VALUE = 999888.9

dir_base = "sim"
mass_ratios = [0.1 * x for x in range(1,6)]
ecc_s = [0.1 * x for x in range(8)]
inc_s = [10.0 * x for x in range(6)]
Ms = [0, 180]

# Make 2-D table at each i
# Make 1-D table for each (u,e)

critical_sma_tables = np.zeros((len(inc_s), len(ecc_s), len(mass_ratios)), dtype = float)
critical_sma_list = np.zeros((len(inc_s), (len(ecc_s) * len(mass_ratios))), dtype = float)

#critical_sma_eject_time_tables = np.zeros((len(inc_s), len(ecc_s), len(mass_ratios)), dtype = float)
#critical_sma_eject_time_list = np.zeros((len(inc_s), (len(ecc_s) * len(mass_ratios))), dtype = float)

for ith, inc in enumerate(inc_s):
    for j,ecc in enumerate(ecc_s):
        for k,u in enumerate(mass_ratios):
            final_sm_axes = [[],[]]
            minimum_stable_sm_axes = [[],[]]
            for i,M in enumerate(Ms):
                u_bin_str = int(round(u * 100, 0))
                e_bin_str = int(round(ecc * 100, 0))
                i_bin_str = int(round(inc, 0))
                M_bin_str = int(round(M, 0))
                
                mask_ejection_fn = "storage/%s_u%02d_e%02d_i%03d_M%03d/mask_ejections.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str)
                final_sm_axes_fn = "storage/%s_u%02d_e%02d_i%03d_M%03d/table_of_final_sm-axes.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str)
                
                mask_ejection_f = open(mask_ejection_fn, "rb")
                ejection_mask = pickle.load(mask_ejection_f)
                mask_ejection_f.close()

                final_sm_axes_f = open(final_sm_axes_fn, "rb")
                table_of_final_sm_axes = pickle.load(final_sm_axes_f)
                final_sm_axes_f.close()

                #print table_of_final_sm_axes
                #print ejection_mask

                # Apply Mask (Only show ejection times of stable planets)
                final_sm_axes[i] = np.ma.array(table_of_final_sm_axes, mask = ejection_mask)

            # Find minimum in each column (for each mean anomaly -- 16 in total)
            minima = [[],[]]
            for i, final_sma in enumerate(final_sm_axes):
                minima[i] = (np.amin(final_sma, axis = 0)).tolist()
                minima[i] = np.round(minima[i], 3)  # round to 3 decimal places

            crit_sma = max([max(minima[0]), max(minima[1])]) # Take largest minimum out of all 16 columns

            print ecc, u, inc
            print minima[0]
            print minima[1]
            print crit_sma
            print

            # Store in tables
            critical_sma_tables[ith,j,k] = crit_sma
            critical_sma_list[ith, j* len(mass_ratios) + k] = crit_sma

            #print
            
# Pickle Files!

pickle_f = open(("%s_moving_critical_sma_ecc_tables.p" % (dir_base)), "wb")
pickle.dump(critical_sma_tables, pickle_f)
pickle_f.close()

pickle_f = open(("%s_moving_critical_sma_ecc_list.p" % (dir_base)), "wb")
pickle.dump(critical_sma_list, pickle_f)
pickle_f.close()

#pickle_f = open(("%s_moving_critical_sma_ecc_eject_tables.p" % (dir_base)), "wb")
#pickle.dump(critical_sma_eject_time_tables, pickle_f)
#pickle_f.close()

#pickle_f = open(("%s_moving_critical_sma_ecc_eject_list.p" % (dir_base)), "wb")
#pickle.dump(critical_sma_eject_time_list, pickle_f)
#pickle_f.close()
            
# Pretty Print Table and List

tables_fn = "%s_moving_critical_sma_ecc_tables.txt" % (dir_base)
list_fn = "%s_moving_critical_sma_ecc_list.txt" % (dir_base)

#eject_tables_fn = "%s_moving_critical_sma_inc_eject_tables.txt" % (dir_base)
#eject_list_fn = "%s_moving_critical_sma_inc_eject_list.txt" % (dir_base)

tables_f = open(tables_fn, "w")
list_f = open(list_fn, "w")

#eject_tables_f = open(eject_tables_fn, "w")
#eject_list_f = open(eject_list_fn, "w")


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