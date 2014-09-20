import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset']="stix"

bodyname='M0_S40'
filename=bodyname+'.aei'
elements = np.loadtxt(filename,skiprows=4)
time=elements[:,0]
I=elements[:,3]
e=elements[:,2]

bodyname='S2'
filename=bodyname+'.aei'
elements = np.loadtxt(filename,skiprows=4)
timeB=elements[:,0]
IB=elements[:,3]
eB=elements[:,2]

T_init = time[0]
T_final = time[-1]

n_frames=2
fig = plt.figure(1, figsize=(7,8))
bottom,left,top = 0.1, 0.15, 0.03
width= 0.82
height = (1.0 - bottom -top)/n_frames
axes1 = [left,height+bottom,width,height]
axes2 = [left,bottom,width,height]

ax= fig.add_axes(axes1)
ax.plot(time,I,'g-')
ax.plot(timeB,IB,'b-')
ax.axis([T_init,T_final,0.0,180])
ax.set_xticklabels([])
y_ticklabels = ax.get_yticks().tolist()
y_ticklabels[0]=' '
ax.set_yticklabels(y_ticklabels)
ax.set_ylabel(r"$i_1\;\mathrm{[deg]}$",size=20)


ax= fig.add_axes(axes2)
ax.semilogy(time,1.0-e,'g-')
ax.semilogy(timeB,1.0-eB,'b-')
ax.axis([T_init,T_final,1e-4,1.5])
y_ticks = ax.get_yticks().tolist()
y_ticks[-1]=' '
ax.set_yticklabels(y_ticks)
ax.set_ylabel(r"$1-e_1$",size=20)
ax.set_xlabel(r"$t\;\mathrm{[yr]}$",size=20)
ax.ticklabel_format(style='sci',axis='x',scilimits=(-10,0),useOffset=False)
ax.xaxis.major.formatter._useMathText = True

#plt.plot(time,I)
#plt.plot(timeB,IB)
#plt.axis([0,time.max(),0,180])
#plt.show()
#plt.semilogy(time,1-e)
#plt.semilogy(timeB,1-eB)
#plt.axis([0,time.max(),0.001,1.5])
plt.show()
