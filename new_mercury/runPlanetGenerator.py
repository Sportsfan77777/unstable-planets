"""
runs planetGenerator.py

That is, supplies arguments for option parser.
Can do successive runs (in different integration directories)
"""

import planetGenerator as PG

import sys

if __name__ in ('__main__'):
    mass_ratios = [0.1 * x for x in range(1,6)]
    #mass_ratios = [0.5]
    eccentricities = [0.1 * x for x in range(8)]
    #eccentricities = [0.6]
    
    # Mean Anomaly of Binary
    zero = round(float(sys.argv[1]), 0)

    # Inclination of System
    inc = round(float(sys.argv[2]), 1)
    inc = 60
    
    min_a = 1.5
    max_a = min_a + 3.9
    for e in eccentricities:
        for m in mass_ratios:
            args = ["--u_bin", m, "--e_bin", e, "--min_sma", min_a, "--max_sma", max_a, "--num_a", "40",
                        "--min_inc", inc, "--time", 50000]
            if zero == 0:
                print "Periapse"
                PG.pseudo_main(args)
            elif zero == 180:
                print "Apoapse"
                args.append("--mean_anom_bin")
                args.append(180.0)
                PG.pseudo_main(args)
            else:
                print "Improper Argument" + zero
          
        #min_a += 0.3
        
        