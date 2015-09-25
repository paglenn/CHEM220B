
import numpy as np
import matplotlib.pyplot as plt

f = np.linspace(3.0,10,100)
T = 2.0
Z = 1 + 10*(3/f - 1) + 2/f

plt.plot(f,Z)
plt.show()
