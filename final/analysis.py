import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt("data.txt")

ntraj = 12
F = a[:,0]
Z = 1. - a[:,1] /6. + a[:,2] / 120.

plt.plot(F,Z)
plt.show()
