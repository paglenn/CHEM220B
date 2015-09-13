import numpy as np
import matplotlib.pyplot as plt
#######################################
#Analysis of xvals
#######################################
exforce=5.0
a = np.genfromtxt('xvals.out')
X = a[:,0]
Y = a[:,1]

plt.xlabel('t')
plt.ylabel(r'$\bar{X}$')
plt.title("f=%.1f"%exforce)
plt.plot(X,Y)
plt.savefig("1ii_f_%.1f.png"%exforce)
