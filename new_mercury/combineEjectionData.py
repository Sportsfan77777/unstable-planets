import numpy as np
import string
import math
import sys

import pickle

"""
This program combines 2 raw collision tables, stability, sm_axes pickle files
from an "inner" run and an "outer" run
"""

""" Read from info.p file (parameters) <<<<------ need to make new info.p file first """

inner = "inner"
outer = "outer"

info = "info.p"
info_out = "info.out"
stability = "stability.p"
sm_axes = "sm_axes.p"
table_of_ejections = "table_of_ejections.p"

directory_M000 = "%s/sim_u1_e4_i060_M000/%s"
directory_M180 = "%s/sim_u1_e4_i060_M180/%s"

up_M000 = directory_M000[3:]
up_M180 = directory_M180[3:]

###################################################################################

info_M000_inner = pickle.load(open(directory_M000 % (inner, info), "rb"))
info_M000_outer = pickle.load(open(directory_M000 % (outer, info), "rb"))

info_M000 = info_M000_inner
info_M000.num_a = info_M000_inner.num_a + info_M000_outer.num_a
info_M000.max_a = info_M000_outer.max_sma
pickle.dump(info_M000, open(up_M000 % (info), "wb"))

info_M180_inner = pickle.load(open(directory_M180 % (inner, info), "rb"))
info_M180_outer = pickle.load(open(directory_M180 % (outer, info), "rb"))

info_M180 = info_M180_inner
info_M180.num_a = info_M180_inner.num_a + info_M180_outer.num_a
info_M180.max_a = info_M180_outer.max_sma
pickle.dump(info_M180, open(up_M180 % (stability), "wb"))

###################################################################################

info_out_M000_inner = open(directory_M000 % (inner, info_out), "r")
info_out_M000_outer = open(directory_M000 % (outer, info_out), "r")
info_out_M000 = open(up_M000 % (info_out), "w")

for line in info_out_M000_inner:
	info_out_M000.write(line)

for line in info_out_M000_outer:
	info_out_M000.write(line)

info_out_M000_inner.close()
info_out_M000_outer.close()
info_out_M000.close()

info_out_M180_inner = open(directory_M180 % (inner, info_out), "r")
info_out_M180_outer = open(directory_M180 % (outer, info_out), "r")
info_out_M180 = open(up_M180 % (info_out), "w")

for line in info_out_M180_inner:
	info_out_M180.write(line)

for line in info_out_M180_outer:
	info_out_M180.write(line)

info_out_M180_inner.close()
info_out_M180_outer.close()
info_out_M180.close()


###################################################################################

stability_M000_inner = pickle.load(open(directory_M000 % (inner, stability), "rb"))
stability_M000_outer = pickle.load(open(directory_M000 % (outer, stability), "rb"))

stability_M000 = np.concatenate([stability_M000_inner, stability_M000_outer])
pickle.dump(stability_M000, open(up_M000 % (stability), "wb"))

stability_M180_inner = pickle.load(open(directory_M180 % (inner, stability), "rb"))
stability_M180_outer = pickle.load(open(directory_M180 % (outer, stability), "rb"))

stability_M180 = np.concatenate([stability_M180_inner, stability_M180_outer])
pickle.dump(stability_M180, open(up_M180 % (stability), "wb"))

####################################################################################

sm_axes_M000_inner = pickle.load(open(directory_M000 % (inner, sm_axes), "rb"))
sm_axes_M000_outer = pickle.load(open(directory_M000 % (outer, sm_axes), "rb"))

sm_axes_M000 = np.concatenate([sm_axes_M000_inner, sm_axes_M000_outer])
pickle.dump(sm_axes_M000, open(up_M000 % (sm_axes), "wb"))

sm_axes_M180_inner = pickle.load(open(directory_M180 % (inner, sm_axes), "rb"))
sm_axes_M180_outer = pickle.load(open(directory_M180 % (outer, sm_axes), "rb"))

sm_axes_M180 = np.concatenate([sm_axes_M180_inner, sm_axes_M180_outer])
pickle.dump(sm_axes_M180, open(up_M180 % (sm_axes), "wb"))

####################################################################################

num_a = len(sm_axes_M000)
num_a_inner = len(sm_axes_M000_inner)
num_a_outer = len(sm_axes_M000_outer)


table_of_ejections_M000_inner = pickle.load(open(directory_M000 % (inner, table_of_ejections), "rb"))
table_of_ejections_M000_outer = pickle.load(open(directory_M000 % (outer, table_of_ejections), "rb"))

table_of_ejections_M000 = np.zeros((num_a, 8))
table_of_ejections_M000[:num_a_inner,:] = table_of_ejections_M000_inner
table_of_ejections_M000[num_a_inner:,:] = table_of_ejections_M000_outer
pickle.dump(table_of_ejections_M000, open(up_M000 % (table_of_ejections), "wb"))


table_of_ejections_M180_inner = pickle.load(open(directory_M180 % (inner, table_of_ejections), "rb"))
table_of_ejections_M180_outer = pickle.load(open(directory_M180 % (outer, table_of_ejections), "rb"))

table_of_ejections_M180 = np.zeros((num_a, 8))
table_of_ejections_M180[:num_a_inner,:] = table_of_ejections_M180_inner
table_of_ejections_M180[num_a_inner:,:] = table_of_ejections_M180_outer
pickle.dump(table_of_ejections_M180, open(up_M180 % (table_of_ejections), "wb"))

####################################################################################



