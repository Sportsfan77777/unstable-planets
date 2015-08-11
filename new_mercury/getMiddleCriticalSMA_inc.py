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
inc_s = [10.0 * x for x in range(10)]
Ms = [0, 180]

# Make 2-D table at each i
# Make 1-D table for each (u,e)

critical_sma_tables = np.zeros((len(inc_s), len(ecc_s), len(mass_ratios)), dtype = float)
critical_sma_list = np.zeros((len(inc_s), (len(ecc_s) * len(mass_ratios))), dtype = float)

critical_sma_eject_time_tables = np.zeros((len(inc_s), len(ecc_s), len(mass_ratios)), dtype = float)
critical_sma_eject_time_list = np.zeros((len(inc_s), (len(ecc_s) * len(mass_ratios))), dtype = float)

for ith, inc in enumerate(inc_s):
    for j,ecc in enumerate(ecc_s):
        for k,u in enumerate(mass_ratios):
            sm_axes = []
            stable_arrays = [[],[]]
            min_eject = [[],[]]
            for i,M in enumerate(Ms):
                u_bin_str = int(round(u * 100, 0))
                e_bin_str = int(round(ecc * 100, 0))
                i_bin_str = int(round(inc, 0))
                M_bin_str = int(round(M, 0))
                stable_file = "storage/%s_u%02d_e%02d_i%03d_M%03d/stability.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str)
                sm_axes_fn = "storage/%s_u%02d_e%02d_i%03d_M%03d/sm_axes.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str) 
                min_eject_fn = "storage/%s_u%02d_e%02d_i%03d_M%03d/min_ejection_times.p" % (dir_base, u_bin_str, e_bin_str, i_bin_str, M_bin_str) 
                
                stability_f = open(stable_file, "rb")
                stable_arrays[i] = pickle.load(stability_f)
                stability_f.close()

                min_eject_f = open(min_eject_fn, "rb")
                min_eject[i] = pickle.load(min_eject_f)
                min_eject_f.close()
                
            sm_axes_f = open(sm_axes_fn, "rb")
            sm_axes = pickle.load(sm_axes_f)
            sm_axes_f.close()
            
            crit_sma = sm_axes[-1] + 900.0
            crit_eject_t = STABLE_VALUE

            done = False

            length_stable = 0
            length_unstable = 0

            found_stable = False
            found_unstable = False
            for (sma, s0, s1, min_eject_t0, min_eject_t1) in zip(sm_axes, stable_arrays[0], stable_arrays[1], min_eject[0], min_eject[1]):
                if not done:
                    # Verify that stable island is critical sma
                    if found_stable:
                        if s0 and s1:
                            length_stable += 1
                        else:
                            found_stable = False
                            found_unstable = True
                            length_unstable = 1
                    # Check if there is a larger island of instability outside of stable island
                    elif found_unstable:
                        if s0 and s1:
                            # End of Island
                            if length_stable > length_unstable:
                                # Found Critical SMA!
                                done = True
                                #print "Done!"
                            else:
                                # Keep Looking! (This is the next new island!)
                                length_stable = 1; length_unstable = 0 
                                found_stable = True; found_unstable = False
                        else:
                            # Grow Island
                            length_unstable += 1
                    # Search for critical sma
                    else:
                        if s0 and s1:
                            crit_sma = sma
                            crit_eject_t = min(min_eject_t0, min_eject_t1)
                            found_stable = True
                            length_stable = 1

            # Store in tables
            critical_sma_tables[ith,j,k] = crit_sma
            critical_sma_list[ith, j* len(mass_ratios) + k] = crit_sma

            critical_sma_eject_time_tables[ith,j,k] = crit_eject_t
            critical_sma_eject_time_list[ith, j* len(mass_ratios) + k] = crit_eject_t

            #print
            
# Pickle Files!

pickle_f = open(("%s_middle_critical_sma_inc_tables.p" % (dir_base)), "wb")
pickle.dump(critical_sma_tables, pickle_f)
pickle_f.close()

pickle_f = open(("%s_middle_critical_sma_inc_list.p" % (dir_base)), "wb")
pickle.dump(critical_sma_list, pickle_f)
pickle_f.close()

pickle_f = open(("%s_middle_critical_sma_inc_eject_tables.p" % (dir_base)), "wb")
pickle.dump(critical_sma_eject_time_tables, pickle_f)
pickle_f.close()

pickle_f = open(("%s_middle_critical_sma_inc_eject_list.p" % (dir_base)), "wb")
pickle.dump(critical_sma_eject_time_list, pickle_f)
pickle_f.close()
            
# Pretty Print Table and List

tables_fn = "%s_middle_critical_sma_inc_tables.txt" % (dir_base)
list_fn = "%s_middle_critical_sma_inc_list.txt" % (dir_base)

eject_tables_fn = "%s_middle_critical_sma_inc_eject_tables.txt" % (dir_base)
eject_list_fn = "%s_middle_critical_sma_inc_eject_list.txt" % (dir_base)

tables_f = open(tables_fn, "w")
list_f = open(list_fn, "w")

eject_tables_f = open(eject_tables_fn, "w")
eject_list_f = open(eject_list_fn, "w")


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
    rows_eject = []
    for j, ecc in enumerate(ecc_s):
        row = (str(ecc)).center(width) + '|'
        row_eject = (str(ecc)).center(width) + '|'
        for k, u in enumerate(mass_ratios):
            row += (str(critical_sma_tables[ith,j,k])).center(width)
            row_eject += (str(critical_sma_eject_time_tables[ith,j,k])).center(width)
        rows.append(row)
        rows_eject.append(row_eject)

    # File 1: the critical 'a'
        
    tables_f.write(row_zero + "\n")
    tables_f.write(dash_row + "\n")
    for r_str in rows:
       tables_f.write(r_str + "\n")
    tables_f.write("\n")

    # File 2: the ejection time at that 'a'

    eject_tables_f.write(row_zero + "\n")
    eject_tables_f.write(dash_row + "\n")
    for r_str in rows_eject:
       eject_tables_f.write(r_str + "\n")
    eject_tables_f.write("\n")

    
# (2) List

row_zero = "(u=u0, e=e0)".center(2*width) + '|'
dash_row = dash + dash + '-'
for inc in inc_s:
    row_zero += (str(inc)).center(width)
    dash_row += dash
    
rows = []
rows_eject = []
for k, u in enumerate(mass_ratios):
  for j, ecc in enumerate(ecc_s):
    row = ("(u=%.1f, e=%.1f)" % (u, ecc)).center(2*width) + '|'
    row_eject = ("(u=%.1f, e=%.1f)" % (u, ecc)).center(2*width) + '|'
    for ith, inc in enumerate(inc_s):
        row += str(critical_sma_list[ith, j*len(mass_ratios) + k]).center(width)
        row_eject += str(critical_sma_eject_time_list[ith, j*len(mass_ratios) + k]).center(width)
    rows.append(row)
    rows_eject.append(row_eject)
  rows.append(dash_row)
  rows_eject.append(dash_row)

# File 1: the critical 'a'
        
list_f.write(row_zero + "\n")
list_f.write(dash_row + "\n")
for r_str in rows:
   list_f.write(r_str + "\n")
list_f.write("\n")

# File 2: the ejection time at that 'a'

eject_list_f.write(row_zero + "\n")
eject_list_f.write(dash_row + "\n")
for r_str in rows_eject:
   eject_list_f.write(r_str + "\n")
eject_list_f.write("\n")

# Close all four files

tables_f.close()
eject_tables_f.close()
list_f.close()
eject_list_f.close()


