import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt('data.txt')
plt.plot(a[:,0],a[:,1])
plt.xlabel('time')
plt.ylabel('q')
plt.savefig('2iv.png')
#plt.show()
