import glob
import subprocess

# Get all the files
existing_files_path = "../../j*0*/sim_u*_e*_i*/"

existing_files = glob.glob(existing_files_path)

print existing_files

print len(existing_files)

for file in existing_files:
    print file
    split = file.split("/")
    link = split[-2] # just the sim_u*e*M* part
    print link
    print
    command = ["ln", "-s", file, link]
    subprocess.call(command)
