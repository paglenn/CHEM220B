import numpy as np
import sys
import scipy.interpolate
import scipy.integrate
import matplotlib.pyplot as plt

integrator = scipy.integrate.simps
sin = np.sin
cos = np.cos
pi = np.pi

nmax = 3
dw = 2.7
rho = 0.033
sigma_a = np.array( [0.1*i for i in range(27,74) ] )

# read in structure data
data = np.genfromtxt('narten_levy.txt',skiprows=2)
k_old = data[:,0]
k_old[0] = 0.001
hww_k_old = data[:,1]

func = scipy.interpolate.interp1d(k_old,hww_k_old,kind='cubic')
k = np.linspace(min(k_old),max(k_old),1000)
hww_k = func(k)
print np.linalg.norm(hww_k_old - func(k_old))
#sys.exit()

#plt.plot(k,hww_k)
#plt.show()
#sys.exit()

# solve matrix eqn. for each radius
Y = []
for a in sigma_a :
    R = 0.5* (dw + a )
    R2 = R * R ; R3 = R2*R
    phi_int = lambda m : 8*np.pi * R3 * ((-1)**m) / ((m+1)*(m+2.)*(m+3))
    # analytic integrals
    I = map(phi_int, range(7))

    phi_k =[]
    sinkR = sin(k*R)
    coskR = cos(k*R)
    k2 = k* k ; k3 = k2 * k ; k4 = k2*k2;k5 = k4 * k ; k6 = k5* k

    p0 = sinkR / k3 - R * coskR / k2
    p0[0] = R3 / 3.

    p1 = sinkR / k3 + 2 * ( coskR - 1. ) / (R*k4)
    p1[0] = -p0[0]/2.

    p2 = -6*sinkR / (R2 * k5) + (2*coskR+4) / (R*k4)
    p2[0] = -p1[0]/5.

    p3 = -6 *sinkR / (R2 * k5) + 24 * (1-coskR)/(R3 * k6 ) - 6 / (R*k4)
    p3[0] = -p2[0]/2.

    phi_k = [p0, p1 , p2, p3 ]

    b = np.zeros(nmax+1)
    A = np.zeros([nmax+1,nmax+1])
    for n in range(nmax+1):
        b[n] = - I[n]
        for m in range(n,nmax+1):
            a_nm = I[n+m]
            #print phi_k[m]
            #plt.plot(k, 1./(4*pi)*phi_k[m])
            k_int = 8*integrator(k*k*phi_k[n]*phi_k[m]*hww_k,k)
            A[m,n] = a_nm + k_int
            A[n,m] = A[m, n]
            #print c
            #if a == 4.0 :
            #plt.plot(k, hww_k)
                #plt.plot(k,2*k*k*phi_k[n]*phi_k[m]*hww_k)
                #plt.show()
            #sys.exit()

    coeffs = np.linalg.solve(A,b)

    #---------------------------------
    # obtain g_aw(r) (first g_aw(k) )
    caw_k = 4*pi* (coeffs[0]*phi_k[0] + coeffs[1]*phi_k[1] +coeffs[2]*phi_k[2] +coeffs[3]*phi_k[3] )
    haw_k = caw_k*( 1. + hww_k)
    npts = 1000
    R = 0.5* (dw + a)
    r = np.linspace(1e-6+R , 6 + R ,npts)
    haw_r =np.zeros(npts)
    for i in range(npts) :
        haw_r[i] = 1./(2.*pi*pi*r[i])*integrator(k*sin(k*r[i])*haw_k,k)
    gaw_r = 1 + haw_r
    fn = "percus_yevick_%.2f.dat"%a
    fout = open(fn,'w')
    for i in range(gaw_r.size):
        fout.write("%f,%f\n"%(r[i],gaw_r[i]))

    fout.close()
    '''
    haa_k = 0.033 * caw_k * haw_k
    r = np.linspace(a , 9,npts)
    haa_r = np.zeros(npts)
    for i in range(npts) :
        haa_r[i] = 1./(2.*pi*pi*r[i])*integrator(k*sin(k*r[i])*haa_k,k)
    gaa_r = 1 + haa_r
    if a == 3.0 or a == 4.0 or a == 5.0 :
        plt.plot(r, gaa_r, label='%s'%a)
    '''
#plt.legend()
#plt.xlabel('r [Angstrom]')
#plt.ylabel(r'$g_{AA}(r)$')
#plt.savefig('1iv.png')
#plt.show()
