import numpy
import string
import os
import sys


G=2.959122082855911E-4


def write_mercury_files(A,B,C,D,E,F,mass,spinx,spiny,spinz,name,filename,body_type='big',coord_type='cart'):

    #if coord_type = 'cart'
    # A,B,C,D,E,F = x, y, z, vx, vy ,vz
    #if coord_type = 'ast'
    # A,B,C,D,E,F = a, e, I, omega, Omega ,M
    

    f = open(filename,'w')
    
    if (body_type == 'big'):
        f.write(")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n")
        f.write(") Lines beginning with `)' are ignored.\n")
        f.write(")---------------------------------------------------------------------\n")
        f.write(" style (Cartesian, Asteroidal, Cometary) = Cartesian\n")
        f.write(" epoch (in days) = 0.0\n")
        f.write(")---------------------------------------------------------------------\n")
    elif (body_type == 'small'):
        f.write(")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n")
        f.write(") Lines beginning with `)' are ignored.\n")
        f.write(")---------------------------------------------------------------------\n")
        if (coord_type == 'cart'):
            f.write(" style (Cartesian, Asteroidal, Cometary) = Cart\n")
        elif (coord_type == 'ast'):
            f.write(" style (Cartesian, Asteroidal, Cometary) = Ast\n")
        else:
            f.close()
            print 'ERROR: Unrecognized coord_type input',coord_type
            return 0
        f.write(")---------------------------------------------------------------------\n")
    else:
        f.close()
        print 'ERROR: Unrecognized body_type input',body_type
        return 0
    
    Nbody = A.shape[0]
 
    if (coord_type == 'cart'):
        for i in range(0,Nbody):
            #bodylabel = ' M%d'%(i)
            bodylabel = name[i]
            bodylabel+='     m=%23.17e r=0.d0 d=1.1'%(mass[i])
            f.write(bodylabel+"\n")
            
            f.write(' %24.17e %24.17e %24.17e\n'% (A[i],B[i],C[i]))
            f.write(' %24.17e %24.17e %24.17e\n'% (D[i],E[i],F[i]))
            f.write(' %8.5f %8.5f %8.5f\n'% (spinx[i],spiny[i],spinz[i]))
    elif (coord_type == 'ast'): 
        for i in range(0,Nbody):
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














