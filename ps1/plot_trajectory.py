#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys

fn = sys.argv[1]
data = np.genfromtxt(fn)
T = data[:,0]
X = data[:,1]
Y = data[:,2]
Z = data[:,3]

plt.plot(T,X,label='x')
plt.plot(T,Y,label='y')
plt.plot(T,Z,label='z')
plt.show()

