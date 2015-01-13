import numpy
import string

import os
import shutil

import sys


G=2.959122082855911E-4

template_path = "/Users/Sportsfan77777/Documents/Fall 2014/ASTRO 4940/planets/tests/template/"
#template_path = "/home/astrosun/mhammer/tests/template/"

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
         
def mkdir_integration(directory, template_dir = None):
    """ make a directory loaded with all of the files ready for an integration """
    if template_dir is None:
       src_path = template_path
    else:
       src_path = template_dir
       
    dst_path = directory
    
    print " *** Copying Directory *** "
    shutil.copytree(src_path, dst_path)
    
    #mkdir(directory, safety = True) # Note: multiple = False (only one directory created)
    
    
def write_param_file(central_mass, duration = 10000000, src_dir = None, dest_dir = None,
                     fn_param_base = "param_base.in"):
                     
    if src_dir is None:
       fn_param_base = template_path + fn_param_base
    else:
       if directory[:0] == "/":
         fn_param_base = src_dir + fn_param_base
       else:
         fn_param_base = src_dir + "/" + fn_param_base
         
    fn_param = "param.in"
       
    if dest_dir is None:
       fn_param = template_path + fn_param
    else:
       if dest_dir[:0] == "/":
         fn_param = dest_dir + fn_param
       else:
         fn_param = dest_dir + "/" + fn_param
    
    f1 = open(fn_param_base, 'r')
    f2 = open(fn_param, 'w')
    
    needs_to_be_filled = f1.read()

    # central mass string
    central_mass_string = "%0.1f" % central_mass

    # duration string
    sci_notation = "{:.2E}".format(duration)
    sci_decimal = sci_notation[:4]
    sci_power = str(int(sci_notation[6:]))
    duration_string = "%sd%s" % (sci_decimal, sci_power)

    complete = needs_to_be_filled % (duration_string, central_mass_string) # Insert (1) simulation time and (2) central mass
    
    print " *** Writing param.in file *** "
        
    f2.write(complete)
    
    f1.close()
    f2.close()
          

def write_mercury_files(A, B, C, D, E, F, mass,
                        spinx, spiny, spinz, 
                        name, filename, directory = None,
                        body_type='big',coord_type='cart'):

    #if coord_type = 'cart'
    # A,B,C,D,E,F = x, y, z, vx, vy ,vz
    #if coord_type = 'ast'
    # A,B,C,D,E,F = a, e, I, argw, Node, M
    
    if not os.path.exists(directory):
       mkdir_integration(directory) # Note: multiple = False (only one directory created)
    # else: it might have been created by write_param_files
    
    if directory[:0] == "/":
       filename = directory + filename
    else:
       filename = directory + "/" + filename

    f = open(filename,'w')
    
    if body_type == 'big':
        f.write(")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n")
        f.write(") Lines beginning with `)' are ignored.\n")
        f.write(")---------------------------------------------------------------------\n")
        f.write(" style (Cartesian, Asteroidal, Cometary) = Cartesian\n")
        f.write(" epoch (in days) = 0.0\n")
        f.write(")---------------------------------------------------------------------\n")
    elif body_type == 'small':
        f.write(")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n")
        f.write(") Lines beginning with `)' are ignored.\n")
        f.write(")---------------------------------------------------------------------\n")
        if coord_type == 'cart':
            f.write(" style (Cartesian, Asteroidal, Cometary) = Cart\n")
        elif coord_type == 'ast':
            f.write(" style (Cartesian, Asteroidal, Cometary) = Ast\n")
        else:
            f.close()
            print 'ERROR: Unrecognized coord_type input', coord_type
            return 0
        f.write(")---------------------------------------------------------------------\n")
    else:
        f.close()
        print 'ERROR: Unrecognized body_type input', body_type
        return 0
    
    Nbody = A.shape[0]
 
    if (coord_type == 'cart'):
        for i in range(0, Nbody):
            #bodylabel = ' M%d'%(i)
            bodylabel = name[i]
            bodylabel+='     m=%23.17e r=0.d0 d=1.1'%(mass[i])
            f.write(bodylabel+"\n")
            
            f.write(' %24.17e %24.17e %24.17e\n'% (A[i],B[i],C[i]))
            f.write(' %24.17e %24.17e %24.17e\n'% (D[i],E[i],F[i]))
            f.write(' %8.5f %8.5f %8.5f\n'% (spinx[i],spiny[i],spinz[i]))
    elif (coord_type == 'ast'): 
        for i in range(0, Nbody):
            #bodylabel = ' M%d'%(i)
            bodylabel = name[i]
            bodylabel+='     m=%23.17e r=0.d0 d=1.1'%(mass[i])
            f.write(bodylabel+"\n")
            
            f.write(' %14.7e %10.8f %9.5f %9.5f %9.5f %9.5f 0 0 0\n'% (A[i],B[i],C[i],E[i],E[i],F[i]))
    else:
        f.close()
        print 'ERROR: Unrecognized coord_type input',coord_type
        return 0
            

    f.close()














