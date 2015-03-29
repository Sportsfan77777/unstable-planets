"""
Take all directories (or, a subset of directories)
and extend the duration of each run.

Then, execute ./m.exe in the appropriate directory to continue each run
"""

import sys
import os
import subprocess
import shutil
import glob
import pickle

new_num_years = 500000
if len(sys.argv) > 1:
	new_num_years = sys.argv[1]

new_num_days = str(int(round(365 * new_num_years, 0)))

# Get all directories (or a subset of directories)
dir_path = "sim_u*_e00_i*_M*"
directories = sorted(glob.glob(dir_path))

# Re-write time in param.dmp
for directory in directories:
	os.chdir(directory)

	fn = "param.dmp"
	fn_tmp = "param-tmp.dmp"

	wr_lines = []

	# Open File
	with open(fn, 'r') as f:
		for line in f:
			# If line contains "stop time"
			if "stop time" in line:
				split_dot = line.split(".")[0]
				days_str = split_dot.split(" ")[-1] # Old Number of Days

				replacement = line.replace(days_str, new_num_days)
				wr_lines.append(replacement)
			else:
				if len(line) > 1:
				    wr_lines.append(line)

	# Open TMP File
	with open(fn_tmp, 'w') as f:
		for line in wr_lines:
			f.write(line)


	# Modify time in info.p also
	info = pickle.load(open('info.p', 'rb'))
	info.time = new_num_years
	pickle.dump(info, open('info.p', 'wb'))

	shutil.move(fn_tmp, fn)
	os.chdir("../")

# Run ./m.exe in background
for directory in directories:
	os.chdir(directory)
	command = ["./m.exe"]
	subprocess.call(command, stdout = subprocess.PIPE)
	os.chdir("../")

# Keep it going (not sure if necessary)
while True:
	pass
	


