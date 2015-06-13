"""
Makes grid plots of stability regions
Stable regions indicated by Green
Half-Stable regions indicated by Yellow
Unstable regions indicated by Red
"""
import matplotlib
matplotlib.use('Agg') # for ssh-ed computers

from matplotlib import pyplot as plot
import numpy as np
import pickle

# Internal Parameters
this_mu = 0.1
this_ecc = 0.1
this_inc_s = [10.0 * x for x in range(8)]

def plot_stability_regions(mu = this_mu, ecc = this_ecc, inc_s = this_inc_s, dir = "sim"):
    # Assume this range of sm-axes contains actual range
    min_a = 1.5
    max_a = 5.0
    num_a = (max_a - min_a) * 30.0 + 1
    sm_axes_base = np.array([np.round(x, 3) for x in np.linspace(min_a, max_a, num_a)])

    # Final 2-D Stability Array (0 = Unstable, 1 = Semi-Stable, 2 = Stable) -- The sum of stability arrays
    # Individual Arrays have (0 = Unstable, 1 = Stable)
    full_stability_array = np.zeros((len(inc_s), num_a))

    for inc_i, inc in enumerate(inc_s):
        dir_base = "storage/%s_u%02d_e%02d_i%03d" % (dir, 100 * mu, 100 * ecc, inc)
        directories = [dir_base, dir_base]
        directories[0] += "_M000"
        directories[1] += "_M180"

        sm_axes_fns = []
        stable_fns = []
        info_fns = []

        for directory in directories:
            sm_axes_fn = directory + "/sm_axes.p"
            stable_fn = directory + "/stability.p"
            info_fn = directory + "/info.p"

            sm_axes_fns.append(sm_axes_fn)
            stable_fns.append(stable_fn)
            info_fns.append(info_fn)

        # Load Stability into Staggered Stability array than spans sm_axes_base, not sm_axes
        sm_axes = pickle.load(open(sm_axes_fns[0], "rb"))
        offset = (sm_axes[0] - min_a) * 30.0

        # Two Stability Arrays
        stability_one = pickle.load(open(stable_fns[0], "rb"))
        stability_two = pickle.load(open(stable_fns[1], "rb"))

        stability = np.zeros(len(sm_axes_base)) # Initialize to unstable (arbitrary)

        for i, (stable_one, stable_two) in enumerate(zip(stability_one, stability_two)):
            # Limit by max_a (Future: change max_a from the start if this happens)
            if (i + offset) < len(stability):
                stability[i + offset] = stable_one + stable_two
                final_i = int(round(i + offset, 0))

        # Fill in remainder of array
        for j in range(final_i, len(stability)):
            stability[j] = 2 # STABLE

        # Store in 2-D Stability Array
        full_stability_array[inc_i, :] = stability

    # Plot 2-D Stability Array
    cmap = plot.get_cmap('RdYlGn')

    ax = plot.gca()

    # Set up bins
    delta_a_bin = 0.033 / 2.0
    a_bins = np.linspace(min_a - delta_a_bin, max_a + delta_a_bin, num_a + 1)
    i_bins = np.array([10 * x - 5 for x in range(len(inc_s) + 2)])
    
    p_map = ax.pcolormesh(a_bins, i_bins, full_stability_array, cmap = cmap)

    plot.title("Stability Regions")
    plot.xlabel("$a_p$ / $a_b$", fontsize = 15)
    plot.ylabel("Inclination $i$", fontsize = 15)

    plot.xlim(min_a - delta_a_bin, max_a + delta_a_bin)
    plot.ylim(inc_s[0] - 5, inc_s[-1] + 5)

    save_fn = "%s_u%02d_e%02d_stabilityMap" % (dir, 100 * mu, 100 * ecc)
    save_fn += ".%s"
    plot.savefig(save_fn % "png")
    plot.savefig(save_fn % "eps", format = "eps", dpi = 1000)


if __name__ == "__main__":
    mass_ratios = [0.05 * x for x in range(2, 11)]
    ecc_s = [0.05 * x for x in range(15)]
    for ecc in ecc_s:
        for mass_ratio in mass_ratios:
            plot_stability_regions(mu = mass_ratio, ecc = ecc, inc_s = this_inc_s)
