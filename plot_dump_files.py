import numpy as np
import matplotlib.pyplot as plt
import glob
#import mercury
import orbital as orb
###############################################################
#CONSTANTS AND UNITS
G = 2.959122082855911e-4
#G=mercury.G
###############################################################
vecOrbit = np.vectorize(orb.orbit)
vecElements = np.vectorize(orb.elements)

small_file = 'small.in'
big_file = 'big.in'

small_file = 'small.dmp'
big_file = 'big.dmp'

bodytag='M'
delta_time= 10000 # in years
ind_count = 0
max_entries = 0
skip_big=6
skip_small=5


MB = 1.0
CentralMass = 0.7

while(True):
    #big bodies file
    f = open(big_file, 'r')
    line_count=0
    body_count=0
    xbig,ybig,zbig = [0.0], [0.0], [0.0]
    vxbig,vybig,vzbig = [0.0], [0.0], [0.0]
    massbig = [CentralMass]
    namebig = ['S1']
    for line in f:
        if (line_count >= skip_big):
            if ((line_count-skip_big)%4 == 0):
                name=line.split()[0]
                mass=np.array(line.split(' m=')[1].split()[0]).astype(np.float)
            if ((line_count-skip_big)%4 == 1):
                x, y, z = np.array(line.split()).astype(np.float)
            if ((line_count-skip_big)%4 == 2):
                vx, vy, vz = np.array(line.split()).astype(np.float)
            if ((line_count-skip_big)%4 == 3):
                body_count+=1
        
                xbig,ybig,zbig  = np.append(xbig,x),np.append(ybig,y),np.append(zbig,z)
                vxbig,vybig,vzbig  = np.append(vxbig,vx),np.append(vybig,vy),np.append(vzbig,vz)
                massbig = np.append(massbig,mass)
                namebig = np.append(namebig,name)

        line_count+=1


    #small bodies file
    f = open(small_file, 'r')
    line_count=0
    body_count=0
    xsmall,ysmall,zsmall = [], [], []
    vxsmall,vysmall,vzsmall = [], [], []
    masssmall = []
    namesmall = []
    for line in f:
        if (line_count >= skip_small):
            if ((line_count-skip_small)%4 == 0):
                name=line.split()[0]
                if ('m=' in line):
                    mass=np.array(line.split(' m=')[1].split()[0]).astype(np.float)
                else:
                    mass=0.0
            if ((line_count-skip_small)%4 == 1):
                x, y, z = np.array(line.split()).astype(np.float)
            if ((line_count-skip_small)%4 == 2):
                vx, vy, vz = np.array(line.split()).astype(np.float)
            if ((line_count-skip_small)%4 == 3):
                body_count+=1
                
                xsmall,ysmall,zsmall  = np.append(xsmall,x),np.append(ysmall,y),np.append(zsmall,z)
                vxsmall,vysmall,vzsmall  = np.append(vxsmall,vx),np.append(vysmall,vy),np.append(vzsmall,vz)
                masssmall = np.append(masssmall,mass)
                namesmall = np.append(namesmall,name)
            
        line_count+=1
    #print xbig,ybig,zbig
    #print vxbig,vybig,vzbig
    #print massbig
    print xsmall,ysmall,zsmall
    print vxsmall,vysmall,vzsmall
    print masssmall

    #calculate barycenter of the inner binary
    print xbig, xbig[0:2]
    X = ((massbig[0:2] * xbig[0:2]).sum() + (masssmall * xsmall).sum()) / (massbig[0:2].sum()  + masssmall.sum())
    Y = ((massbig[0:2] * ybig[0:2]).sum() + (masssmall * ysmall).sum()) / (massbig[0:2].sum()  + masssmall.sum())
    Z = ((massbig[0:2] * zbig[0:2]).sum() + (masssmall * zsmall).sum()) / (massbig[0:2].sum()  + masssmall.sum())
    VX = ((massbig[0:2] * vxbig[0:2]).sum() + (masssmall * vxsmall).sum()) / (massbig[0:2].sum()  + masssmall.sum())
    VY = ((massbig[0:2] * vybig[0:2]).sum() + (masssmall * vysmall).sum()) / (massbig[0:2].sum()  + masssmall.sum())
    VZ = ((massbig[0:2] * vzbig[0:2]).sum() + (masssmall * vzsmall).sum()) / (massbig[0:2].sum()  + masssmall.sum())

   
    p,e,I,g,h,f,l,E = vecElements(xbig-X,ybig-Y,zbig-Z,(vxbig-VX),(vybig-VY),(vzbig-VZ),G*MB)
    plt.plot(xbig-X,ybig-Y,'bo',mew=0,ms=8.0)
    for k in range(0,xbig.shape[0]):
        anomalies=np.arange(l[k],l[k]+2*np.pi,0.01)
        xorb, yorb, zorb, vxorb, vyorb, vzorb = vecOrbit(p[k],e[k],I[k],g[k],h[k],anomalies,G*MB)
        #plt.plot(xorb,yorb,'k-')

    p,e,I,g,h,f,l,E = vecElements(xsmall-X,ysmall-Y,zsmall-Z,vxsmall-VX,vysmall-VY,vzsmall-VZ,G*MB)
    print p/(1-e),e
    plt.plot(xsmall-X,ysmall-Y,'go',mew=0,ms=4.0)
    for k in range(0,xsmall.shape[0]):
        anomalies=np.arange(l[k],l[k]+2*np.pi,0.01)
        xorb, yorb, zorb, vxorb, vyorb, vzorb = vecOrbit(p[k],e[k],I[k],g[k],h[k],anomalies,G*MB)
        plt.plot(xorb,yorb,'k-')

    plt.plot([0],[0],'k+',ms=8.0)
    plt.axis([-4,4,-4,4])
    plt.show()
