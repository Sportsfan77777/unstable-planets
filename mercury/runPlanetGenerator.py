"""
runs planetGenerator.py

That is, supplies arguments for option parser.
Can do successive runs (in different integration directories)
"""

import planetGenerator as PG


if __name__ in ('__main__'):
    mass_ratios = [0.1 * x for x in range(1,6)]
    #mass_ratios = [0.1]
    eccentricities = [0.1 * x for x in range(8)]
    eccentricities = [0.7]
    
    min_a = 3.7
    
    for e in eccentricities:
        max_a = min_a + 0.9
        for m in mass_ratios:
            args = ["--u_bin", m, "--e_bin", e, "--min_sma", min_a, "--max_sma", max_a]
            PG.pseudo_main(args)
        min_a += 0.3