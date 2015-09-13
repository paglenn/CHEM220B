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

plt.plot(t/T,r,label='r(t)')
plt.plot(t/T,v,label='v(t)')
plt.title('Position/Velocity')
plt.xlabel(r'$t/2\pi [\sqrt{m/k}]$')
plt.xticks(np.arange(0.0,nT,1.0))
plt.legend()
plt.savefig('2ii_posvel.png')
plt.clf()

plt.subplot(211)
plt.plot(t/T,K,label='K')
plt.plot(t/T,U,label='U')
plt.xticks(np.arange(0.0,nT,1.0))
plt.legend()
#plt.xlabel('t')
plt.xlabel(r'$t/2\pi [\sqrt{m/k}]$')
plt.subplot(212)
plt.plot(t/T,E,label=r'$E_{T}$')
plt.xticks(np.arange(0.0,nT,1.0))
#plt.xlabel('t')
plt.xlabel(r'$t/2\pi [\sqrt{m/k}]$')
plt.savefig('2ii_energy.png')

