
import numpy as np
theory  = np.genfromtxt('theory.dat')
data = np.genfromtxt('hist.dat')
#T = 0.2
dt = 0.01

import matplotlib.pyplot as plt


plt.plot(theory[:,0],theory[:,1], label="Boltzmann")
plt.plot(data[:,0],data[:,1],'ks', ms =5, label="data")
plt.legend()
plt.xlabel("q")
plt.ylabel("p(q)")
plt.xlim([-1.5,1.5])
plt.ylim([0.9*np.amin(data[:,1]),1.1*np.amax(data[:,1])])
plt.title("dt = %.2f"%dt)
plt.yscale('log')
plt.savefig("hist_dt_%.2f.png"%dt)
plt.show()

