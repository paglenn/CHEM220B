import scipy.special
import matplotlib.pyplot as plt
import numpy as np
import math

X = np.linspace(-2, 2.0, 100)
erfX = scipy.special.erf(X)
phi = 0.5 + 0.5 * erfX

vals = [x for x in phi if (x >= 0.1 and x <=0.9) ]
set1 = set( np.where(phi >= 0.05)[0]  )
set2 = set( np.where(phi <= 0.95)[0]  )
indices = list( set1 & set2 )
print X[indices]
print phi[indices]
#print vals


plt.plot(X,phi)
plt.xlabel(r'$q_{0} \sqrt{ \beta m \omega^{2}} $')
plt.ylabel(r'$\phi_{B}(q_{0})$')
plt.savefig('erf.png')


