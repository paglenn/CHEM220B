#!/usr/bin/env python
import matplotlib.pyplot as plt
import os
import numpy as np

files = [f for f in os.listdir('.') if 'density_hist' in f ]
for fn in files:
    A = np.genfromtxt(fn,delimiter=' ')
    X = A[:,0]
    Y = A[:,1]
    tag = fn[13:-4]
    plt.ylim([-0.1,1.1*max(Y)])
    plt.plot(X,Y,'.',label=r'$\bar{\rho}$=%s'%tag)
    plt.xlabel('r/d')
    plt.ylabel(r'$g(r)$')
    plt.title(r'$\rho^{*}$=%s RDF'%tag)
    plt.savefig('plot_rho=%s.png'%tag)
    plt.clf()
#plt.legend()
#plt.show()



