"""
Creates a scatter plot of ejection times against final sm-axes
"""

import numpy as np
import string
import math
import sys
import os
import subprocess

import matplotlib
matplotlib.use('Agg') # for ssh-ed computers
from matplotlib import pyplot as plot

import glob
import pickle

from structures import *
from avatar import *

from mercury import G as BigG

# Read info.p
pickle_fn = "info.p"
pickle_f = open(pickle_fn, "rb")
o = pickle.load(pickle_f)
pickle_f.close()

# Get Table of Final Semimajor Axes
pickle_f = open("table_of_final_sm-axes.p", "rb")
sm_axis_table = pickle.load(pickle_f)
pickle_f.close()

# Get Table of Ejection Times
pickle_f = open("table_of_ejections.p", "rb")
ejectionTable = pickle.load(pickle_f)
pickle_f.close()

# Reshape Arrays into 1-D (this is only useful for plotting purposes)
sm_axis_table = sm_axis_table.reshape(sm_axis_table.size)
ejectionTable = ejectionTable.reshape(ejectionTable.size)

x = sm_axis_table[:]
for i,z in enumerate(ejectionTable):
	if z == 99.9:
		ejectionTable[i] = 10 ** 5.95 # Set non-ejectees apart

y = [1000*i for i in ejectionTable]

# Set Up Plot
plot.scatter(x, y)

plot.xlim(o.min_sma - 0.1, o.max_sma + 0.1)
plot.yscale('log')

plot.savefig("t_vs_a.png")
plot.show()

plot.cla()

# Set Up Resonant Plot
x_r = [i**1.5 for i in sm_axis_table] # resonant thing

plot.xlim((o.min_sma - 0.1)**(1.5), (o.max_sma + 0.1)**(1.5))
plot.yscale('log')
plot.scatter(x_r, y)

plot.savefig("t_vs_T.png")
plot.show()


