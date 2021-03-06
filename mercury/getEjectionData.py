import numpy as np
import string
import sys

import pickle

"""
This program converts the collision + ejection data from info.out
into a table to match the format from the Holman/Wiegert paper
showing the ejection times for each test mass.

Certain parameters must be specified near the top of the file.
Outputs showing a '+' indicate the test mass survived.
"""

"""
phase_str = ['Peri', 'Apo']
phase = 0  # Periapse = 0; Apoapse = 1

N = 8         # Number of Mean Anomalies
num_a = 10    

a_b = 1.00  # Semi-Major Axis of Binary

a_min = 3.1 * a_b  # Note: should be greater than a_b !!!! (P-type)
a_max = 4.0 * a_b  
"""

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
Format: M#a_S#b collided with the central body at #c.#d years
Format: M#a_S#b ejected at #c.#d years

Read #a as semi-major axis, #b as mean anomaly, #c as ejection time in 1000s of years
"""


""" M """
mean_anomalies = np.ones(N) * 2 * np.pi / N
for i in xrange(N):
	mean_anomalies[i] *= i  # i = 0 to N - 1
	
m_deg_array = [int(round(x * 180.0 / np.pi)) for x in mean_anomalies]
#print m_deg_array # M Names

""" a """
semi_major_axes = np.ones(num_a) * a_min
if num_a > 1:
	a_range = a_max - a_min
	step_size = a_range / (num_a - 1) 
	for i in xrange(num_a):
		semi_major_axes[i] += (i * step_size)
		semi_major_axes[i] = round(semi_major_axes[i], 1) # This is specific for a = %.1f
		
sm_array = [int(10 * x) for x in semi_major_axes]
#print sm_array # S Names

ejectionTable = np.zeros((num_a, N)) + 99.9

# Parse Ejection Data Into Table

for line in data:
	split_str = line.split('_')
	M_str = split_str[0] # This is M#a
	
	nxt_str = split_str[1]
	nxt_split = nxt_str.split(' ')
	
	S_str = nxt_split[0] # This is S#b
	Eject_val = int(float(nxt_split[-2])) / 1000.0 # This is #c (in kyr)
	Eject = round(Eject_val, 1)
	
	#print M_str[2:], S_str[1:], Eject
	
	Mi = m_deg_array.index(int(M_str[2:]))
	Si = sm_array.index(int(S_str[1:]))
	
	#print Mi, Si
	ejectionTable[Si][Mi] = Eject
	
	# parse into M and a (directly or indirectly?)

# Simple Print	
#print ejectionTable, '\n'

# Pretty Print
width = 7
rowzero = phase_str[phase].center(width) + '|'
dash = '-' * width
dash_row = dash
for x in m_deg_array:
	rowzero += (str(x)).center(width)
	dash_row += dash
	
stable_array = np.zeros(len(sm_array))
	
rows = []
count = 0
for i,y in enumerate(sm_array):
	row = str(y / 10.).center(width) + '|'
	stable = True
	for ej in ejectionTable[i]:
		s = ""
		if ej == 99.9:
			s = "+".center(width)
		else:
			s = str(ej).center(width)
			stable = False
		row += s
	rows.append(row)
	if stable == True:
	   stable_array[i] = 1
	
# Write Pretty Print to File
if len(sys.argv) > 1:
   fn = sys.argv[1]  # write to a file (if provided)
   f = open(fn, 'a')
else:
   f = open("eject.t", 'w')
	
print rowzero
f.write(rowzero + "\n")
print dash_row
f.write(dash_row + "\n")
for r_str in rows:
    print r_str
    f.write(r_str + "\n")

f.close()
   
# Write Stability Array

pickle_f = open("stability.p", "wb")
pickle.dump(stable_array, pickle_f)
pickle_f.close()

# Write SM-Axes Array

pickle_f = open("sm_axes.p", "wb")
pickle.dump(semi_major_axes, pickle_f)
pickle_f.close()

# Write Mean Anomaly Array

pickle_f = open("mean_anomalies.p", "wb")
pickle.dump(m_deg_array, pickle_f)
pickle_f.close()



