import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import sys
integrate = scipy.integrate.simps

P = np.array([0.02 ,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
lenP = P.size
d = 1. ; d3= d * d * d
kd = d* np.linspace(0.01, 100, 1000 )
_kd = 1./kd
_KD = [_kd ** x for x in range(7)]
coskd = np.cos(kd)
sinkd = np.sin(kd)
phi_1 = -coskd * _KD[2] + sinkd * _KD[3]
phi_2 = -coskd * _KD[2] + 2 * sinkd * _KD[3] + 2 * (coskd-1)  * _KD[4]
phi_3 = -coskd * _KD[2] + 4 * sinkd * _KD[3] + 12 * coskd * _KD[4] - 24 * sinkd * _KD[5] - 24 * (coskd -1 ) * _KD[6]

print phi_1[:10]
print phi_2[:10]
print phi_3[:10]

npoints_r = 100
R = np.linspace(0.01, 5, npoints_r)
x = kd
#f, AX = plt.subplots(10)
for i in range(lenP):
    p = P[i]
    eta = np.pi * p /6.
    L_1 = - (1 + 2* eta) ** 2. / (1 - eta)**4.
    L_2 = 6 * eta * (1 + 0.5 * eta) **2. / (1 - eta )**4.
    L_3 = 0.5 * eta * L_1
    c_x = L_1 * phi_1 + L_2 * phi_2 + L_3 * phi_3
    fourpi = 4 * np.pi
    I_R = np.zeros(R.shape)
    c_R = np.copy(I_R)
    g_R = np.copy(I_R)

    I_x = fourpi ** 2. * p * d3 * c_x * c_x  / (1 - fourpi * p * c_x )

    for j in range(npoints_r) :
        r = R[j]
        sin_x = np.sin( x * r / d )
        f_x  = x * sin_x * I_x
        pre = 1./(2 * np.pi**2.  * r * d*d)
        I_R[j] = pre * integrate(f_x , x )

        if r < d : c_R[j] = L_1 + L_2 * (r/d) + L_3 * (r/d)**3.
        else: c_R[j] = 0.0

    #AX[i].plot(R,I_R)
    #AX[i].set_title(r'$\rho^{*}$ = %f'%p)
    h_R = I_R + c_R
    plt.plot(R,h_R)
    plt.title(r'$\rho^{*}$ = %f'%p)
    plt.savefig('1iii_rho_%.2f.png'%p)
    plt.clf()

    # write h_R to file
    g_R = h_R + 1.
    ofn = 'g_R_rho_%.2f.dat'%p
    fout = open(ofn, 'w')
    for j in range(npoints_r):
        r = R[j]
        g_r = g_R[j]
        fout.write('%f\t%f\n'%(r,g_r))
    fout.close()

    #if p == 0.5 :
    #    plt.plot(R,g_R)
    #    plt.show()
    #    sys.exit()



