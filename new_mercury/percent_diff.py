"""
Calculate the percent differences between a_st and a_crit
"""

import numpy as np
import pickle

import sys

# Parameters
mass_ratios = [0.05 * x for x in range(2, 11)]
ecc_s = [0.05 * x for x in range(15)]
inc_s = [10.0 * x for x in range(10)]

### Real Data ###
def crit_sma(mass_ratio, ecc, inc, crit_type = "safe"):
    """ dictionary-like calls for the tables of critical sm-axes """
    mass_ratios = np.array([abs(0.05 * x - mass_ratio) for x in range(2,11)]) # [0.1 - mass_ratio, 0.15 - mass_ratio, ..., 0.5 - mass_ratio]
    mu_i = np.argmin(mass_ratios) # The selected mass ratio will be '0' after mass_ratio is subtracted from the array

    ecc_s = np.array([abs(0.05 * x - ecc) for x in range(15)])
    ecc_i = np.argmin(ecc_s)

    inc_s = np.array([abs(10.0 * x - inc) for x in range(10)])
    inc_i = np.argmin(inc_s)

    # 3 Possibilities (safe = UCO from initial, moving = UCO from median, else = LCO)
    if crit_type == "safe":
        crit_type = "safe_"
    elif crit_type == "moving":
        crit_type = "moving_safe_"
    else:
        crit_type = ""

    safe_crit_sma_table = pickle.load(open("/Users/Sportsfan77777/planets/critical/sim_%scritical_sma_inc_tables.p" % crit_type, "rb"))
    return safe_crit_sma_table[inc_i, ecc_i, mu_i] # stored as [i = inc, j = ecc, k = mass_ratio]


#### Differences ####
data_points = len(inc_s) * len(ecc_s) * len(mass_ratios)

differences = np.zeros(data_points)


count = 0
for mass_ratio in mass_ratios:
	for ecc in ecc_s:
		for inc in inc_s:
			# Real Data
			real_crit = crit_sma(mass_ratio, ecc, inc, crit_type = "")
			real_st = crit_sma(mass_ratio, ecc, inc, crit_type = "safe")

			differences[count] = real_st / real_crit

			count += 1

# Sort
differences = np.sort(differences)

# Analysis
def analyze(array):
	print "Median: %.3f" % (np.median(array))
	print "Mean (w/ Std Dev): %.3f +/- %.3f" % (np.mean(array), np.std(array))
	print "Five Maximum Values: %.3f, %.3f, %.3f, %.3f, %.3f" % (array[-1], array[-2], array[-3], array[-4], array[-5])


print "Percent Difference"
analyze(differences)
print





