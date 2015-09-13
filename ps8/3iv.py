import numpy as np
import matplotlib.pyplot as plt
T = 0.4
sk = 10 # skip
num = np.genfromtxt('3a.dat')
a = np.genfromtxt('3iv.dat')
plt.plot(num[::2,0],num[::2,1],'ko', label = 'exact potential')
plt.plot(a[::sk,0],a[::sk,1],'ro', label = 'approx')
plt.xlabel('q_0')
plt.ylabel(r'$\phi_{B}(q_{0})$')
plt.title('T = %.2f'%T)
plt.legend(loc='best')
plt.savefig('3iv_T_%.2f.png'%T)
plt.ylim([-0.1,1.1])
plt.show()
