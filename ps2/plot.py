#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt


thefile = sys.argv[1]
if len(sys.argv) > 2: skip = int(sys.argv[2])
R = np.genfromtxt(thefile,skip_header=skip)
X = R[:,0]
Y1 = R[:,1]
Y2 = R[:,2]

plt.plot(X,Y2,'o-', label="data")
plt.xlabel(r'$k$')
plt.ylabel(r'$\langle \mathfrak{R}_{k}^{2} \rangle$')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.title('2 (ii)')
plt.legend(loc=3,ncol=1, mode="expand", borderaxespad=0.)
plt.savefig('2ii.png')
#plt.show()

