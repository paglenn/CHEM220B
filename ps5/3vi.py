import numpy as np
import matplotlib.pyplot as plt


dt = 3.0
for dt in [2.0, 2.5, 3.0, 4.0]:

    dt2 = dt * dt
    T = 2 * np.pi
    nT = 12
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

    plt.plot(t/T,r,label='dt = %.1f'%dt)
    #cosndt = [1]+ [np.cosh(n*dt ) for n in range(nsteps)]
    #plt.plot(t/T,cosndt,label='cos(n*dt)')
    plt.title('Analysis, dt > 2')
    plt.xlabel(r'$t/2\pi [\sqrt{m/k}]$')

plt.xticks(np.arange(0.0,nT,1.0))
plt.legend(loc=2)
plt.savefig('3vi.png')

