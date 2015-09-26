"""
Calculate the residuals between the co-planar fit and the actual values at any inclination

Parameter: inclination
"""

import numpy as np
import pickle

import sys

# Parameters
mass_ratios = [0.05 * x for x in range(2, 11)]
ecc_s = [0.05 * x for x in range(15)]
inc = int(sys.argv[1])


### Fit Functions ###
info_fit_crit = pickle.load(open("fit_crit_inc00.p", "rb"))
info_fit_st = pickle.load(open("fit_st_inc00.p", "rb"))

popt_crit = info_fit_crit["popt"]
popt_st = info_fit_st["popt"]

def fit(x, const, e_c, esq_c, m_c, me_c, msq_c, esqmsq_c):
    m = x[0]
    e = x[1]
    return const + e_c * (e) + esq_c * (e*e) + m_c * (m) + me_c * (m*e) + msq_c * (m*m) + esqmsq_c * (m*m*e*e)

def fit_crit(mu, e):
    return fit([mu, e], popt_crit[0], popt_crit[1], popt_crit[2], popt_crit[3], popt_crit[4], popt_crit[5], popt_crit[6])

def fit_st(mu, e):
    return fit([mu, e], popt_st[0], popt_st[1], popt_st[2], popt_st[3], popt_st[4], popt_st[5], popt_st[6])

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


#### Residuals ####
data_points = len(ecc_s) * len(mass_ratios)

residuals_crit = np.zeros(data_points)
residuals_st = np.zeros(data_points)

count = 0
for mass_ratio in mass_ratios:
	for ecc in ecc_s:
		# Real Data
		real_crit = crit_sma(mass_ratio, ecc, inc, crit_type = "")
		real_st = crit_sma(mass_ratio, ecc, inc, crit_type = "safe")

		# Fits
		fit_crit_value = fit_crit(mass_ratio, ecc)
		fit_st_value = fit_st(mass_ratio, ecc)

		residuals_crit[count] = abs((fit_crit_value - real_crit) / real_crit)
		residuals_st[count] = abs((fit_st_value - real_st) / real_st)

		count += 1

# Sort
residuals_crit = np.sort(residuals_crit)
residuals_st = np.sort(residuals_st)

# Analysis
def analyze(array):
	print "Median: %.3f" % (np.median(array))
	print "Mean (w/ Std Dev): %.3f +/- %.3f" % (np.mean(array), np.std(array))
	print "Five Maximum Values: %.3f, %.3f, %.3f, %.3f, %.3f" % (array[-1], array[-2], array[-3], array[-4], array[-5])

print
print " *** Inclination %d ***" % (inc)

print "a_crit:"
analyze(residuals_crit)
print

print "a_st"
analyze(residuals_st)
print





