import numpy as np
import matplotlib.pyplot as plt
T = 0.2
sk = 10 # skip
a = np.genfromtxt('3b.dat')
plt.plot(a[:,0],a[:,1],'ro', label = 'approx')
plt.xlabel('q_0')
plt.ylabel(r'$\phi_{B}(q_{0})$')
plt.title('T = %.2f'%T)
#plt.legend(loc='best')
plt.ylim([-0.1,1.1])
plt.savefig('3iii_T_%.2f_b.png'%T)
plt.show()
