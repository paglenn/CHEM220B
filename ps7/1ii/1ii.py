import numpy as np
import matplotlib.pyplot as plt

f = 10.0
a = np.genfromtxt("xvals.out")

t = a[:,0]
Xk1 = a[:,1]
Xk2 = a[:,2]
Xk3 = a[:,3]

plt.plot(t,Xk1,'.', label=r"k = 12$\pi/L$ ")
plt.plot(t,Xk2, '.',label=r"k = 24$\pi/L$ ")
plt.plot(t,Xk3, '.',label=r"k = 36$\pi/L$ ")
plt.legend()
plt.title('f = %.1f'%f)
plt.xlabel("t")
plt.ylabel(r'$\bar{X}(t)/(\beta f_{0})$')
plt.savefig('1ii_f_%.1f.png'%f)
plt.show()
