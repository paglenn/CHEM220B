import numpy as np
import matplotlib.pyplot as plt

f = np.linspace(3.0,100.0,100)
L = 1.
a = 1.
K = np.sqrt(10.*(f * L * L - 3*a) / (L**4. * f ) )


plt.plot(f,K)
plt.show()
