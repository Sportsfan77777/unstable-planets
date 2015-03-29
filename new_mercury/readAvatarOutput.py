import numpy as np
import string
import math
import sys
import os
import subprocess

import glob
import pickle

from id import ID_Manager
from structures import *
from avatar import *

from mercury import G as BigG

def convert_to_deg(angle):
    return 180.0 * angle / np.pi

""" Read from info.p file (parameters) """

pickle_fn = "info.p"
pickle_f = open(pickle_fn, "rb")
o = pickle.load(pickle_f)
pickle_f.close()

""" Read IDs from ids.p file (ids map to the old format of M#a_S_ """

pickle_fn = "ids.p"
pickle_f = open(pickle_fn, "rb")
id_dict = pickle.load(pickle_f)
pickle_f.close()

ID_manager = ID_Manager()
ID_manager.read()

""" Handle Input Planet """

if (len(sys.argv) is not 3) and (len(sys.argv) is not 2):
    print "Proper Use: python readAvatarOutput.py id_num"
    print "Proper Use: python readAvatarOutput.py a M"
else:
    if len(sys.argv) == 2:
        id_num = int(sys.argv[1])
        id_name = "ID_%04d" % id_num
    elif len(sys.argv) == 3:
        sma = float(sys.argv[1])
        M = int(sys.argv[2])

        planet_name = ("A%1." + ("%df" % o.sep_sma) + "_M%1.0f") % (sma, M)
        id_name = ID_manager.get_id(planet_name)

    save_fn = "%s_elements.npy" % id_name

    # Final Array: a, e, i, W, w, M
    final_array = np.load(save_fn)

    # Format
    width = 10
    time_width = 12

    decimal_places_t = 1
    decimal_places_a = 3
    decimal_places_e = 4
    decimal_places_angle = 2

    # Header
    rowzero = "time".center(time_width)
    rowzero += "a".center(width) + "e".center(width) + "i".center(width)
    rowzero += "node".center(width) + "argw".center(width) + "M".center(width)

    print rowzero

    for i in range(len(final_array[:,0])):
        time = final_array[i, 0]

        elements = OrbitalElements()
        elements.set_a(final_array[i, 1])
        elements.set_e(final_array[i, 2])
        elements.set_i(final_array[i, 3])
        elements.set_W(final_array[i, 4])
        elements.set_w(final_array[i, 5])
        elements.set_M(final_array[i, 6])

        row_i = (str(round(time, decimal_places_t))).center(time_width)
        row_i += (str(round(elements.get('a'), decimal_places_a))).center(width)
        row_i += (str(round(elements.get('e'), decimal_places_e))).center(width)
        row_i += (str(round(convert_to_deg(elements.get('inc')), decimal_places_angle))).center(width)
        row_i += (str(round(convert_to_deg(elements.get('node')), decimal_places_angle))).center(width)
        row_i += (str(round(convert_to_deg(elements.get('argw')), decimal_places_angle))).center(width)
        row_i += (str(round(convert_to_deg(elements.get('M')), decimal_places_angle))).center(width)

        print row_i

    print





