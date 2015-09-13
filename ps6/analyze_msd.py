import numpy as np
import matplotlib.pyplot as plt

MSDvt = np.genfromtxt("msd_rho_0.80.dat")
t = MSDvt[:,0]
MSD = MSDvt[:,1]
#Dest = 1./6 *  np.mean(np.diff(MSDvt[MSDvt.shape[0]/2:-1,1]))
#print Dest
npts = MSD.size
start = 9*npts/10
slope = np.polyfit(t[start:], MSD[start:],1)[0]
print slope / 6.
plt.plot(t,MSD)
plt.ylim([-0.1,max(MSDvt[:,1])])
plt.xlabel('t')
plt.ylabel('MSD')
plt.title(r'$\rho^{*} = 0.8$')
plt.savefig('1v_rho_0.80.png')
plt.show()
