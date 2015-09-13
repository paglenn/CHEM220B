#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt


thefile = sys.argv[1]
skip = int(sys.argv[2])
R = np.genfromtxt(thefile,skip_header=skip)
X = R[:,0]
Y = R[:,1]
G = R[:,2]

plt.plot(X,Y,'s', label="data")
plt.plot(X,G,'r-',label="gaussian")
plt.title('D=%sd'%thefile[8])
plt.xlabel(r'$N_{v}$')
plt.ylabel(r'P($N_{v}$)')
plt.legend()
plt.savefig('dist-%s'%thefile[8])
#plt.show()

