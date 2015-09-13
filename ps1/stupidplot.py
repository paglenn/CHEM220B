import matplotlib.pyplot as plt
from scipy.interpolate import spline
import numpy as np

D = [2,3,4]
Nv = [2.09, 7.07, 16.76]
pv = list(Nv)
dNv2 = [0.83,2.11,4.30]
Dnew = np.linspace(min(D), max(D), 300 )
pv_smooth = spline(D,pv,Dnew)
#plt.plot(D,Nv,'bs',label=r'$\langle N_{v} \rangle $')
#plt.plot(Dnew,pv_smooth,'k-',label=r'$\bar{\rho}v $')
plt.plot(D,dNv2,'rs',label=r'$\langle \delta N_{v} ^{2} \rangle $')
plt.xlim([0.8*min(D),1.2*max(D)])
plt.xlabel('D/d')
plt.ylabel(r'$\langle \delta N_{v} ^{2} \rangle $')
#plt.legend()
plt.savefig('nv_var.png')
