"""
make fits for a_crit and a_st
"""

import pickle
import numpy as np
from scipy.optimize import curve_fit

# Parameters
mass_ratios = [0.05 * x for x in range(2, 11)]
ecc_s = [0.05 * x for x in range(15)]
inc_s = [10.0 * x for x in range(10)]

# Read in a_crit
tables_f = open("sim_critical_sma_inc_tables.p", "rb")
critical_sma_tables = pickle.load(tables_f)
tables_f.close()
# Read in a_st
safe_tables_f = open("sim_safe_critical_sma_inc_tables.p", "rb")
safe_critical_sma_tables = pickle.load(safe_tables_f)
safe_tables_f.close()

# Format xdata and ydata
data_points = len(mass_ratios) * len(ecc_s)

xdatas = np.zeros((len(inc_s), 2, data_points))
ydatas = np.zeros((len(inc_s), 2, data_points))

for i, inc in enumerate(inc_s):
    count = 0
    for j, ecc in enumerate(ecc_s):
        for k, mass_ratio in enumerate(mass_ratios):
            xdatas[i, 0, count] = mass_ratio
            xdatas[i, 1, count] = ecc

            ydatas[i, 0, count] = critical_sma_tables[i, j, k]
            ydatas[i, 1, count] = safe_critical_sma_tables[i, j, k]

            count += 1


# Fit Function
def fit(x, const, e_c, esq_c, m_c, me_c, msq_c, esqmsq_c):
    m = x[0]
    e = x[1]
    return const + e_c * (e) + esq_c * (e*e) + m_c * (m) + me_c * (m*e) + msq_c * (m*m) + esqmsq_c * (m*m*e*e)

def save_fit(popt, pcov):
    d = {}
    coeffs = {}

    # Format
    perr = np.sqrt(np.diag(pcov))

    popt_t = np.round(popt, 5)
    popt_r = np.round(popt, 2)

    pcov_t = np.round(pcov, 5)
    pcov_r = np.round(pcov, 2)

    perr_t = np.round(perr, 5)
    perr_r = np.round(perr, 2)

    pcombo = np.vstack((popt_r, perr_r)).reshape((-1,), order='F') # Interleaves popt with perr

    # By Coefficient
    coeffs['const'] = popt_r[0]
    coeffs['e_c'] = popt_r[1]
    coeffs['esq_c'] = popt_r[2]
    coeffs['m_c'] = popt_r[3]
    coeffs['me_c'] = popt_r[4]
    coeffs['msq_c'] = popt_r[5]
    coeffs['esqmsq_c'] = popt_r[6]

    # General
    d['popt'] = popt_t
    d['pcov'] = pcov_t
    d['perr'] = perr_t
    d['coeffs'] = coeffs
    d['print'] = "%.2f + %.2fe + %.2fe^2 + %.2fm + %.2fme + %.2fm^2 + %.2f(me)^2" % tuple(popt_r)
    d['printE'] = "(%.2f +/- %.2f) + (%.2f +/- %.2f)e + (%.2f +/- %.2f)e^2 + (%.2f +/- %.2f)m + (%.2f +/- %.2f)me + (%.2f +/- %.2f)m^2 + (%.2f +/- %.2f)(me)^2" % tuple(pcombo)
    return d

# Calculate Fit
save_fn = "fit_%s_inc%02d.p"
for i, inc in enumerate(inc_s):
    # a_crit = ydata[0]
    popt, pcov = curve_fit(fit, xdatas[i], ydatas[i,0])
    d = save_fit(popt, pcov)
    pickle.dump(d, open(save_fn % ("crit", inc), "wb"))

    # a_st = ydata[1]
    popt, pcov = curve_fit(fit, xdatas[i], ydatas[i,1])
    d = save_fit(popt, pcov)
    pickle.dump(d, open(save_fn % ("st", inc), "wb"))



