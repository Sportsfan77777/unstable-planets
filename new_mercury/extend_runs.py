"""
Take all directories (in the future, a subset of directories)
and extend the duration of each run.

Then, execute ./m.exe in the appropriate directory to continue each run
"""

import os
import glob


# Get all directories
dir_path = "sim_u*_e*_i*_M*"
directories = sorted(glob.glob(dir_path))

print directories