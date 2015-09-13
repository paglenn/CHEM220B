import numpy as np
import matplotlib.pyplot as plt

dt = 0.1
dt2 = dt * dt
T = 2 * np.pi
nT = 10
tf = nT * T
nsteps = int(tf/dt )
r = np.zeros(nsteps+1)
v = np.zeros(nsteps+1)
r[0] = 1.

for i in range(nsteps):
    r[i+1] = r[i] + dt * v[i] - 0.5 * dt2 * r[i]
    v[i+1] = v[i] - 0.5 * dt * (r[i+1] + r[i])

t = np.linspace(0.0,tf,nsteps+1)
K = 0.5 * v * v
U = 0.5 * r * r
E = K + U

plt.plot(r,v)
plt.axis('equal')
plt.title('Phase Portrait')
plt.xlabel(r'r [$\ell$]')
plt.ylabel(r'v [$\ell/\tau$]')
plt.savefig('2iii_portrait.png')
plt.clf()

