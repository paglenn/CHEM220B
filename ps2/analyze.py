#!/usr/bin/env python
import math
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('density_hist.dat')
Mu = [np.mean(data[:,i]) for i in range(3)]
Var = [np.var(data[:,i]) for i in range(3)]
print Var
X = [np.histogram(data[:,i],bins=20, density=True) for i in range(3)]
B = [x[1] for x in X  ]
P = [x[0] for x in X  ]
logP = [[] for x in X  ]
logPG = list()

for i in range(3):

    B[i] = [0.5*(B[i][j]+B[i][j+1]) for j in range(len(B[i])-1)]
    logPG.append([])
    logP[i] = [math.log(p) for p in P[i] if p != 0 ]

    for binCenter in B[i] :
        logPGx = - ( binCenter - Mu[i]) **2. / (2* Var[i])  \
                - 0.5*math.log(2*math.pi* Var[i])
        logPG[i].append(logPGx)

cm = ['r','b','g']
for i in range(3):
    plt.plot(B[i],logP[i],cm[i]+'s',label=r'k=%s$\pi$/L'%(i+1))
    plt.plot(B[i],logPG[i],cm[i],label='Gaussian' )

plt.xlabel(r'$\mathfrak{R}_{k}$')
plt.ylabel(r'$P(\mathfrak{R}_{k})$')
plt.title('2(iv)')
plt.legend(loc=3,ncol=3, mode="expand", borderaxespad=0.)
plt.savefig('hist.png')
plt.show()




