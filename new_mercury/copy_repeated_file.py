"""
Copies a file to analyze the output for each and every directory (or a subset of directories)

The file must be specified.
It is assumed that the file is in the code directory.
"""

import sys
import os
import subprocess
import shutil
import glob

code_path = "data/mhammer/code"

# Get all directories (or a subset of directories)
dir_path = "sim_u*_e*_i*_M*"
directories = sorted(glob.glob(dir_path))

# Copy file into each directory
for directory in directories:
	filename = "copy_repeated_file.py" # <<<<<<<<<<<<<<======================== Insert filename here ============]]]]]]
	full_filepath = "%s/%s" % (code_path, filename)
	command = ["cp", full_filepath, directory] 
	subprocess.call(command)