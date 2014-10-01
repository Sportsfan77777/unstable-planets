"""
General Useful Functions
"""

import os

def mkdir(directory, safety = False, multiple = False):
  """ make a directory if it doesn't exist, with options for protection and sub-directories """
  if not(multiple):
     try:
        os.mkdir(directory)
     except:
        print "\t(" , directory, "already exists)"
        if (safety):
           raise Exception("This directory already exists. Do not try to delete it!")
  else:
     # Parse "/": Try to create all directories
     count = directory.count("/")
     if count == 0:
         mkdir(directory, safety = safety) # only called with something like 'dir/sub_dir/' <--- Don't do this (end in '/')
     else:
         slash_indices = [i for i,x in enumerate(directory) if x == '/']
         for i in slash_indices:
             sub_dir = directory[:i] # technically, not a sub-directory
             mkdir(sub_dir, safety = safety, multiple = True) # tree directories
         mkdir(directory, safety = safety) # final directory