import numpy as np
import matplotlib.pyplot as plt

nbins = 501
N = 3375
rho = 0.3
R = np.genfromtxt("rdf.dat")
counts = 1.*R.size ;
binContents, binEdges = np.histogram(R, bins =nbins)
binContents = binContents / counts;
binCenters = np.zeros(nbins)
for i in range(nbins):
    binCenters[i] = 0.5* (binEdges[i] + binEdges[i+1])
dR = np.diff(binCenters)[0]; # should be a constant array
pf = (N-1.)/(4*np.pi*rho * dR )
g_R = pf * binContents / (binCenters*binCenters)
g_R /= 2
plt.plot(binCenters, g_R)
plt.xlabel('R')
plt.ylabel('g(R)')
plt.title(r'$\rho^{*}$ = 0.3')
plt.savefig('1iii_rho_p3.png')
plt.show()
