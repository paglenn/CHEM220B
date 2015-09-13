import numpy as np
import matplotlib.pyplot as plt
T = 0.4

a = np.genfromtxt('3a.dat')
plt.plot(a[:,0],a[:,1],'k.')
plt.xlabel('q_0')
plt.ylabel(r'$\phi_{B}(q_{0})$')
plt.title('T = %.2f'%T)
plt.ylim([-0.1,1.1])
plt.savefig('3iii_T_%.2f.png'%T)
plt.show()
