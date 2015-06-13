import numpy as np
import string
import math
import sys

import pickle
from id import ID_Manager

"""
This program converts the collision + ejection data from info.out
into a table to match the format from the Holman/Wiegert paper
showing the ejection times for each test mass.

Certain parameters must be specified near the top of the file.
Outputs showing a '+' indicate the test mass survived.

Currently, it is only used to produce "looseStability.p"
that is used in plot_stability_regions for 'loose' = True
For loose stability, each 'a'-ring is recorded by the fraction of stable planets
instead of just whether there is a stable planet
"""

STABLE_VALUE = 999888.9

""" Read from info.p file (parameters) """

pickle_fn = "info.p"
pickle_f = open(pickle_fn, "rb")
o = pickle.load(pickle_f)
pickle_f.close()

phase_str = ['Peri', 'Apo']
if o.mean_anom_bin == 0:
   phase = 0
else:
   phase = 1   

N = o.num_M
num_a = o.num_a

a_b = o.a_bin

a_min = o.min_sma
a_max = o.max_sma


""" Read IDs from ids.p file (ids map to the old format of M#a_S_ """

pickle_fn = "ids.p"
pickle_f = open(pickle_fn, "rb")
id_dict = pickle.load(pickle_f)
pickle_f.close()

ID_manager = ID_Manager()
ID_manager.read()

""" Read from into.out file (output data from m.exe = mercury) """

fn = 'info.out'
f = open(fn)
lines = f.readlines()
f.close()

""" Get only the lines with the collision + ejection data """

data = []
for line in lines:
    a = line.find('ejected') > -1
    b = line.find('collided') > -1
    c = line.find('hit') > -1
    if a or b or c:
        data.append(line)
        #print line
        
""" 
Format: ID_#%04d collided with the central body at #c.#d years
Format: ID_#%04d ejected at #c.#d years

Read ID_# and look up corresponding mean_anomaly and semimajor_axis, #c as ejection time in 1000s of years
"""


""" M """
mean_anomalies = np.ones(N) * 2 * np.pi / N
for i in xrange(N):
    mean_anomalies[i] *= i  # i = 0 to N - 1
    
m_deg_array = [int(round(x * 180.0 / np.pi)) for x in mean_anomalies]
#m_deg_array = [("03d" % x) for x in mean_anomalies_array]  # <<<<<------- FIX THIS!!!!!!! (for id.py, not here...)

#print m_deg_array # M Names

""" a """
semi_major_axes = np.ones(num_a) * a_min
if num_a > 1:
    a_range = a_max - a_min
    step_size = a_range / (num_a - 1) 
    for i in xrange(num_a):
        semi_major_axes[i] += (i * step_size)
        semi_major_axes[i] = round(semi_major_axes[i], o.sep_sma) # <<--- Note the change from x.x to x.xx!!

sm_array = [round(x, o.sep_sma) for x in semi_major_axes] # <<<<---- ARCHAIC use!!!!!!
#print sm_array # S Names

ejectionTable = np.zeros((num_a, N)) + STABLE_VALUE

# Parse Ejection Data Into Table

for line in data:
    split_str = line.split('_')
    #ID_label = split_str[0] # This is 'ID' <<<---- Assert this?
    
    nxt_split = line.split(' ')
    
    ID_str = nxt_split[1] # This is the ID number
    ID_name = id_dict[ID_str]
    id_split = ID_name.split('_')

    A_str = id_split[0]
    M_str = id_split[1]

    eject_val = int(float(nxt_split[-2])) / 1000.0 # This is #c (in kyr)
    eject = round(eject_val, 2)
    
    #print M_str[1:], A_str[1:], Eject
    
    Ai = sm_array.index(float(A_str[1:]))
    Mi = m_deg_array.index(int(M_str[1:]))
    
    #print Mi, Si
    ejectionTable[Ai][Mi] = eject  ### <<<<------ Save this table somehow (as dictionary?)
    
    # parse into M and a (directly or indirectly?)

# Simple Print  
#print ejectionTable, '\n'

# Pretty Print
width = 9
rowzero = phase_str[phase].center(width) + '|'
dash = '-' * width
dash_row = dash
for x in m_deg_array:
    rowzero += (str(x)).center(width)
    dash_row += dash
    
# Stable Array
stable_array = np.zeros(len(sm_array))
# Minimum Ejection Time Array
min_eject_array = np.zeros(len(sm_array)) + STABLE_VALUE
    
rows = []
count = 0
str_y_base = "%.0" + ("%d" % o.sep_sma) + "f"
for i,y in enumerate(sm_array):
    str_y = str_y_base % y
    row = str_y.center(width) + '|'
    stable_count = 0
    for ej in ejectionTable[i]:
        s = ""
        if ej == 999888.9:
            s = "+".center(width)
            stable_count += 1
        else:
            s = ("%.2f" % ej).center(width)
        row += s
    rows.append(row)
    # Mark Stability
    stable_array[i] = 0.125 * stable_count
    # Find Min Eject Time (for one row back)
    if (i + 1) < len(sm_array):
        min_eject_array[i + 1] = min(ejectionTable[i])

    
# Write Pretty Print to File
"""
if len(sys.argv) > 1:
   fn = sys.argv[1]  # write to a file (if provided)
   f = open(fn, 'a')
else:
   f = open("eject.t", 'w')
    
header = "Directory: %s" % o.integration_dir
print header
f.write(header + "\n")
print rowzero
f.write(rowzero + "\n")
print dash_row
f.write(dash_row + "\n")
for r_str in rows:
    print r_str
    f.write(r_str + "\n")

f.close()
"""

# Write Minimum Ejection Time Array
"""
pickle_f = open("min_ejection_times.p", "wb")
pickle.dump(min_eject_array, pickle_f)
pickle_f.close()
"""
# Write Stability Array

pickle_f = open("looseStability.p", "wb")
pickle.dump(stable_array, pickle_f)
pickle_f.close()

# Write SM-Axes Array
"""
pickle_f = open("sm_axes.p", "wb")
pickle.dump(semi_major_axes, pickle_f)
pickle_f.close()
"""
# Write Mean Anomaly Array
"""
pickle_f = open("mean_anomalies.p", "wb")
pickle.dump(m_deg_array, pickle_f)
pickle_f.close()
"""
# Write Ejection Table!
# (i) as table
"""
pickle_f = open("table_of_ejections.p", "wb")
pickle.dump(ejectionTable, pickle_f)
pickle_f.close()
"""
# (ii) as array of dictionaries



