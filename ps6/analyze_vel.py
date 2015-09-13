import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps

data  = np.genfromtxt("vcf_rho_0.80.dat")
t = data[:,0]
V = data[:,1]
D = 1./3 * simps(V,t)
print D
plt.plot(t,V)
plt.xlabel('t')
plt.ylabel('VCF')
plt.title(r'$\rho^{*} = 0.8$')
plt.savefig('1v_rho_0.80.png')
plt.show()



