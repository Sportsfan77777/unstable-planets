import numpy as np
import mercury
import orbital as orb
import numpy.random as ran
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

###############################################################
#CONSTANTS AND UNITS
G=mercury.G
ran.seed(42)
##############################################################
vecOrbit = np.vectorize(orb.orbit)
vecElements = np.vectorize(orb.elements)

#SYSTEM PARAMETERS


#central binary orbital parameters
aB=0.25 # semi-major axis of binary
eB=0.01 # eccentricity of binary
MB = 2.0 # total mass of binary
muB=0.5  # mass ratio of binary mu=m2/(m1+m2)


N_samples = 1 #number of planets/test particles

#If the planet semimajor axis is not random
ap=4.0

"""###############################################"""

#find the critical semimajor axis for stable outer orbits
#using the results of Holman & Weigert 1999
acrit = aB * (1.6 + 5.1 * eB - 2.22 * eB**2 + 4.12 * muB \
                  - 4.27 * eB * muB - 5.09 * muB**2 + 4.61 * eB**2 * muB**2)
print "\nCritical planetary (outer) orbit: acrit=",acrit,"\n"
a_min = acrit
a_max = 5*acrit

"""##################################################################################"""
#Generate cartesian coordinates for the stellar-mass bodies 
#(in a reference frame centered around the most massive of the inner stars


# (1) First, obtain motion of the inner binary respect to the binary center-of-mass

omegaB= ran.random()*180.0
OmegaB = ran.random()*180.0
MeanAnomB = ran.random()*180.0
#If only two stars, we define the reference plane as the binary plane
IB=0

X0, Y0, Z0, VX0, VY0, VZ0 = vecOrbit(aB*(1-eB),eB,IB*np.pi/180.0,
                                     np.pi+omegaB*np.pi/180.0,OmegaB*np.pi/180.0,
                                     MeanAnomB*np.pi/180.0,G*MB)
X1, Y1, Z1, VX1, VY1, VZ1 = vecOrbit(aB*(1-eB),eB,IB*np.pi/180.0,
                                     omegaB*np.pi/180.0,OmegaB*np.pi/180.0,
                                     MeanAnomB*np.pi/180.0,G*MB)
#scale individual orbits
scale0 = muB
X0, Y0, Z0, VX0, VY0, VZ0 = X0*scale0, Y0*scale0, Z0*scale0, VX0*scale0, VY0*scale0, VZ0*scale0
scale1 = 1.0 - muB
X1, Y1, Z1, VX1, VY1, VZ1 = X1*scale1, Y1*scale1, Z1*scale1, VX1*scale1, VY1*scale1, VZ1*scale1

#check that the barycenter is at the origin
baryX = ((1 - muB) * X0 + muB * X1)
baryY = ((1 - muB) * Y0 + muB * Y1)
baryZ = ((1 - muB) * Z0 + muB * Z1)
baryVX = ((1 - muB) * VX0 + muB * VX1)
baryVY = ((1 - muB) * VY0 + muB * VY1)
baryVZ = ((1 - muB) * VZ0 + muB * VZ1)

print "Binary barycenter:",baryX,baryY,baryZ
print "                  ",baryVX,baryVY,baryVZ

#Also check that the relative motion satisfies the given orbital parameters
deltaX = X1 - X0
deltaY = Y1 - Y0
deltaZ = Z1 - Z0
deltaVX = VX1 - VX0
deltaVY = VY1 - VY0
deltaVZ = VZ1 - VZ0
pB,eccB,iB,gB,nB,fB,lB,EB = vecElements(deltaX, deltaY,deltaZ, deltaVX, deltaVY, deltaVZ,G*MB)
print "Central binary orbital elements:",pB/(1-eccB),eccB,iB


#convert to coordinates respect to the first body
#first big bodies
x0, y0, z0, vx0, vy0, vz0 = X0-X0, Y0-Y0, Z0-Z0, VX0-VX0, VY0-VY0, VZ0-VZ0
x1, y1, z1, vx1, vy1, vz1 = X1-X0, Y1-Y0, Z1-Z0, VX1-VX0, VY1-VY0, VZ1-VZ0

#data big bodies (ignore the first one)
# if only two massive bodies, then this list only contains data for one object
x  = np.array([x1])
y  = np.array([y1])
z  = np.array([z1])
vx = np.array([vx1])
vy = np.array([vy1])
vz = np.array([vz1])
mass=np.array([muB*MB])
name=np.array(['S2']) #the central object's name is S1

#######################################################################################################################
# (2) Second stage
# Generate a population or orbiting test particles
# (2)(a) 1st, generate population of orbital elements

#eccentricities
e=ran.rayleigh(0.01,N_samples)
#inclinations
I=ran.rayleigh(0.01*180.0/np.pi,N_samples)
#masses (test particles)
m=np.zeros(N_samples)
#semimajor axes
#a=10.0**(ran.random_sample(N_samples)*(np.log10(a_max)-np.log10(a_min)) + np.log10(a_min))
a=10.0**(np.arange(0,1.0,1.0/N_samples)*(np.log10(a_max)-np.log10(a_min)) + np.log10(a_min))
#longitudes or pericenter
omega=ran.random_sample(N_samples)*180.0
#longitudes of ascending nodes
Omega=ran.random_sample(N_samples)*180.0
#mean anomalies
MeanAnom = ran.random_sample(N_samples)*180.0

# (2)(b) 2nd, generate a list of names for the bodies
label=[]
for i in range(0,N_samples):
    label.append(' M%d'%(i))

a[:]=ap
I[:]=IB
e[:]=0.0001
Omega[:]=OmegaB

# (2)(c) 3rd, coonvert to cartesian coordinates
Xp, Yp, Zp, VXp, VYp, VZp = vecOrbit(a*(1-e),e,I*np.pi/180.0,omega*np.pi/180.0,
                                     Omega*np.pi/180.0,MeanAnom*np.pi/180.0,G*MB)

#now small bodies
xp, yp, zp, vxp, vyp, vzp = Xp-X0, Yp-Y0, Zp-Z0, VXp-VX0, VYp-VY0, VZp-VZ0
#convert back to asteroidal coordinates
peridistp,ep,Ip,omegap,Omegap,TrueAnomp, MeanAnomp,EccAnomp = vecElements(xp, yp, zp, vxp, vyp, vzp,G*MB)
ap = peridistp/(1 - ep)


#Finally, generate a MERCURY-readable file for small and big objects
mercury.write_mercury_files(x, y, z, vx, vy, vz, mass,
                            np.zeros(2),np.zeros(2),np.zeros(2),
                            name,'big.in')

#mercury.write_mercury_files(ap,ep,Ip*180.0/np.pi,omegap,Omegap,MeanAnomp,m,
#                            np.zeros(N_samples),np.zeros(N_samples),np.zeros(N_samples),
#                            label,'small.in',body_type='small',coord_type='ast')

mercury.write_mercury_files(xp,yp,zp,vxp,vyp,vzp,m,
                            np.zeros(N_samples),np.zeros(N_samples),np.zeros(N_samples),
                            label,'small.in',body_type='small',coord_type='cart')
                            
"""####################################################"""

plot_orbits = 1
plot_barycentric = 1

if plot_orbits:
    if plot_barycentric:
        plt.plot([0]+X0,[0]+Y0,'bo',mew=0,ms=7.0)
        plt.plot(x+X0,y+Y0,'bo',mew=0,ms=7.0)
        plt.plot(xp+X0,yp+Y0,'go',mew=0,ms=4.0)
    else:
        plt.plot([0],[0],'bo',mew=0,ms=7.0)
        plt.plot(x,y,'bo',mew=0,ms=7.0)
        plt.plot(xp,yp,'go',mew=0,ms=4.0)
    for k in range(0,xp.shape[0]):
        #anomalies=np.arange(MeanAnom[k],MeanAnom[k]+2*np.pi,0.01)
        #xorb, yorb, zorb, vxorb, vyorb, vorb = vecOrbit(ap[k]*(1-ep[k]),ep[k],Ip[k]*np.pi/180.0,omegap[k],
        #                                                Omegap[k],anomalies,G*MB*(1-muB))
        #xorb, yorb, zorb, vxorb, vyorb, vorb = vecOrbit(a[k]*(1-e[k]),e[k],I[k]*np.pi/180.0,omega[k]*np.pi/180.0,
        #                                                Omega[k]*np.pi/180.0,anomalies,G*MB)
        
        if plot_barycentric:
            peridistp,ep,Ip,omegap,Omegap,TrueAnomp, MeanAnomp,EccAnomp = vecElements(xp[k]+X0, yp[k]+Y0, zp[k]+Z0, vxp[k]+VX0, 
                                                                                      vyp[k]+VY0, vzp[k]+VZ0,G*MB)
        else:
            peridistp,ep,Ip,omegap,Omegap,TrueAnomp, MeanAnomp,EccAnomp = vecElements(xp[k], yp[k], zp[k], vxp[k], 
                                                                                      vyp[k], vzp[k],G*MB)
        anomalies=np.arange(MeanAnomp,MeanAnomp+2*np.pi,0.01)
        xorb, yorb, zorb, vxorb, vyorb, vorb = vecOrbit(peridistp,ep,Ip,omegap,Omegap,anomalies,G*MB)
        
        plt.plot(xorb,yorb,'k-')

    plt.axis([-3,3,-3,3])
    plt.show()

    if plot_barycentric:
        plt.plot([0]+Y0,[0]+Z0,'bo',mew=0,ms=7.0)
        plt.plot(y+Y0,z+Z0,'bo',mew=0,ms=7.0)
        plt.plot(yp+Y0,zp+Z0,'go',mew=0,ms=4.0)
    else:
        plt.plot([0],[0],'bo',mew=0,ms=7.0)
        plt.plot(y,z,'bo',mew=0,ms=7.0)
        plt.plot(yp,zp,'go',mew=0,ms=4.0)
    for k in range(0,yp.shape[0]):
        plt.plot(yorb,zorb,'k-')

    plt.axis([-3,3,-3,3])
    plt.show()

    if plot_barycentric:
        fig=plt.figure()
        ax=Axes3D(fig)
        ax.plot([0]+X0,[0]+Y0,[0]+Z0,'bo')
        ax.plot(x+X0,y+Y0,z+Z0,'bo')
        ax.plot(xp+X0,yp+Y0,zp+Z0,'go')
        for k in range(0,yp.shape[0]):
            ax.plot(xorb,yorb,zorb,'k-')
        ax.axis('equal')
        ax.auto_scale_xyz([-1,1],[-1,1],[-1,1])
        plt.show()

