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
    
    zero = round(float(sys.argv[1]), 0)
    inc = sys.argv[2]
    
    ecc_pl = round(float(sys.argv[3]), 2)
    ecc_sc = round(100.0 * ecc_pl, 0)
    dir_name = "ecc%02d" % ecc_sc
    
    inclinations = [round(float(inc), 1)]
    
    min_a = 1.5
    
    for incl in inclinations:
      for e in eccentricities:
        #max_a = 2.0
        max_a = min_a + 3.9
        for m in mass_ratios:
          if zero == 0:
            print "Periapse"
            args = ["--u_bin", m, "--e_bin", e, "--min_sma", min_a, "--max_sma", max_a, "--num_a", "40",
                    "--min_inc", incl,
                    "--dir", dir_name, "--min_ecc", ecc_pl]
            PG.pseudo_main(args)
          elif zero == 180:
            print "Apoapse"
            args = ["--u_bin", m, "--e_bin", e, "--min_sma", min_a, "--max_sma", max_a, "--num_a", "40", 
                    "--min_inc", incl, "--mean_anom_bin", 180.0,
                    "--dir", dir_name, "--min_ecc", ecc_pl]
            PG.pseudo_main(args)
          else:
            print "Improper Argument" + zero
          
        #min_a += 0.3
        
        