//detailed balance algorithm 
#include<iostream>
#include<cmath> 
#include<gsl/gsl_rng.h> 
using std::cout; 
using std::endl; 

int main() { 

	//system and simulation variables 
	const int np = 1000 ; 
	double L ,d ,d2,  Nside;  // volume is L^3 
	double x[np], y[np], z[np] ;
	double dx, dy, dz, dr2 ; 
	double xtrial, ytrial, ztrial ; 
	int i, j, k, nsweeps, ip, jp;
	int xcell, ycell, zcell ; // for PBCs
	double scale , acc, pacc ; 

	nsweeps = 1000; 
	d = 1 ; d2 = d*d ; 
	L = pow(2*np,1./3) * d ;  
	scale = 0.25 ; 
	Nside = pow(np,1./3) ; 

	// for gsl rng 
	double u ; 
	const gsl_rng_type * ranT;  
	gsl_rng * ranr; 
	gsl_rng_env_setup(); 
	ranT = gsl_rng_default;
	ranr = gsl_rng_alloc (ranT);

	// initialize box 
	ip = 0 ; 
	for(int i = 0 ; i < Nside ; i++) { 
		for(int j = 0 ; j < Nside; j++) { 
			for(int k = 0 ; k < Nside ; k++) {
			ip++; 	
			x[ip] = i*L/Nside;
			y[ip] = j*L/Nside;
			z[ip] = k*L/Nside;
			}	
		}
	} 

	//return 0; 
	// dynamics 
	for ( i = 1 ; i <= nsweeps; i++) { 

		for (j = 0 ; j < np ; j++) {
			
			// choose particle to move 
			u = gsl_rng_uniform(ranr); 
			ip = u*np ; 

			// trial move 
			u = gsl_rng_uniform(ranr);
			xtrial = x[ip] + scale * ( 2*u-1 )  ; 
			u = gsl_rng_uniform(ranr);
			ytrial = y[ip] + scale * ( 2*u-1 )  ; 
			u = gsl_rng_uniform(ranr);
			ztrial = z[ip] + scale * ( 2*u-1 )  ; 


			for (int jp = 0 ; jp < np ; jp++) { 
				if (jp == ip) continue; 
				dx = x[jp] - xtrial; 
				dx -= L * round(dx/L); 
				dy = y[jp] - ytrial; 
				dy -= L * round(dy/L); 
				dz = z[jp] - ztrial; 
				dz -= L * round(dz/L); 
				dr2 = dx*dx + dy*dy + dz*dz;  
				acc = (dr2 < d*d) ? 0 : 1 ; 
				if (acc == 0 ) break; 

			}	
			pacc += acc ; 

			if (acc) {
				x[ip] = xtrial; 
				y[ip] = ytrial; 
				z[ip] = ztrial; 
			}
			

		}
	printf("%d %f %f %f\n",i, x[np/2],y[np/2],z[np/2]) ; 
	
	} // end loop over sweeps 
	pacc /= ( nsweeps*np) ; 
	cout << pacc << endl ;  

	return 0; 
}




	






