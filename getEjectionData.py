import numpy as np
import string

"""
This program converts the collision + ejection data from info.out
into a table to match the format from the Holman/Wiegert paper
showing the ejection times for each test mass.

Certain parameters must be specified near the top of the file.
Outputs showing a '+' indicate the test mass survived.
"""

phase_str = ['Peri', 'Apo']
phase = 0  # Periapse = 0; Apoapse = 1

N = 8         # Number of Mean Anomalies
num_a = 10    

a_b = 1.00  # Semi-Major Axis of Binary

a_min = 3.1 * a_b  # Note: should be greater than a_b !!!! (P-type)
a_max = 4.0 * a_b  


""" Read from file """

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
		
sm_array = [int(10 * x) for x in semi_major_axes]
#print sm_array # S Names

ejectionTable = np.zeros((num_a, N)) + 99.9

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
	
rows = []
count = 0
for y in sm_array:
	row = str(y / 10.).center(width) + '|'
	for ej in ejectionTable[count]:
		s = ""
		if ej == 99.9:
			s = "+".center(width)
		else:
			s = str(ej).center(width)
		row += s
	rows.append(row)
	count += 1
	
print rowzero
print dash_row
for r_str in rows:
	print r_str
