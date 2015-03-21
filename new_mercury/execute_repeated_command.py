"""
Run a command to analyze the output for each and every directory (or a subset of directories)

The command must be specified.
Whether this is a background process is indicated with the inclusion of an argument
"""

import sys
import os
import subprocess
import shutil
import glob

# Background?
background = False
if len(sys.argv) > 1:
	background = True

# Get all directories (or a subset of directories)
dir_path = "sim_u*_e*_i*_M*"
directories = sorted(glob.glob(dir_path))

# Execute command in each directory
for directory in directories:
	os.chdir(directory)

	command = [] # <<<<<<<<<<<<<<<<<<========================== Insert command here ============
	if background:
		subprocess.Popen(command, stdout=subprocess.PIPE)
	else:
		subprocess.call(command)

	os.chdir("../")