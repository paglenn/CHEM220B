import numpy as np
import matplotlib.pyplot as plt

dt = 0.1
dt2 = dt * dt
T = 2 * np.pi
nT = 10000
tf = nT * T
nsteps = int(tf/dt )
r = np.zeros(nsteps+1)
v = np.zeros(nsteps+1)
r[0] = 1.

for i in range(nsteps):
    r[i+1] = r[i] + dt * v[i] - 0.5 * dt2 * r[i]
    v[i+1] = v[i] - 0.5 * dt * (r[i+1] + r[i])

# construct histograms
binContents, bins = np.histogram(r,bins=100, density=True)
bins = [0.5 * ( bins[i-1] + bins[i]) for i in range(1,bins.size)]
bins = np.array(bins)
plt.plot(bins,binContents)
plt.plot(bins, (1./np.pi) * 1./(np.sqrt(1.-bins*bins)))
plt.xlabel('r')
plt.ylabel('P(r)')
plt.show()
