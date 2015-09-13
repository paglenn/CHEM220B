import numpy as np
import matplotlib.pyplot as plt

nbins = 50
T = 1.5
V = np.genfromtxt("vel.out")
binContents, binEdges = np.histogram(V, bins=nbins, density=True)
binCenters = np.zeros(nbins)
for i in range(nbins):
    binCenters[i] = 0.5* (binEdges[i] + binEdges[i+1])
pred = 1./np.sqrt(2*np.pi*T) * np.exp(- binCenters**2. / (2*T) )
plt.plot(binCenters, binContents,'s', label='data')
plt.plot(binCenters, pred , label='theory')
plt.legend()
plt.ylabel(r'P($v_{x}$)')
plt.xlabel(r'$v_{x}$')
plt.title(r'$\rho^{*}$ = 0.8')
plt.savefig('1ii_rho_p8.png')
plt.show()

