//detailed balance algorithm 
#include<iostream>
#include<cmath> 
#include<gsl/gsl_rng.h> 
#include<ctime>
#include<fstream>
using std::cout; 
using std::ofstream; 
using std::endl; 

#define PI 3.14159265358979323

int main() { 
	int start_time = time(0); 
	//system and simulation variables 
	const int np = 1000 ; 
	double L ,d ,d2,  Nside ;  // volume is L^3 
	double x[np], y[np], z[np] ;
	double dx, dy, dz, dr2 ; 
	double xtrial, ytrial, ztrial ; 
	int i, j, k, ii, jj, kk, ip, jp;
	int xcell, ycell, zcell ; // for PBCs
	int nsweeps, nskip; 
	double scale , acc, pacc ; 
	int iwatch; // tagged particle  
	double freePath, Ncoll, mfp  ; 

	nsweeps = 10000; 
	nskip = 500 ; 
	d = 1. ; d2 = d*d ; 
	L = pow(2*np,1./3) * d ;  
	scale = 0.44 ; 
	Nside = pow(np,1./3) ; 
	iwatch = np/2 ; 
	
	// for gsl rng 
	double u ; 
	const gsl_rng_type * ranT;  
	gsl_rng * ranr; 
	gsl_rng_env_setup(); 
	ranT = gsl_rng_default;
	ranr = gsl_rng_alloc (ranT);
	gsl_rng_set(ranr, time(0)); 

	/* ofstream outfile ; 
	char fn[50];
	sprintf(fn,"dists_D=%.1f.dat",m) ; 
	outfile.open(fn) ; 
	*/

	// Mersenne Twister 
    // MTRand rng(90210);

	// initialize box 
	ip = 0 ; 
	for(i = 0 ; i < Nside ; i++) { 
		for( j = 0 ; j < Nside; j++) { 
			for( k = 0 ; k < Nside ; k++) {
			x[ip] = i*L/Nside;
			y[ip] = j*L/Nside;
			z[ip] = k*L/Nside;
			ip++; 	
			}	
		}
	} 

	// init histogram 
	Ncoll = 0 ; 
	freePath = 0; 

	//return 0; 
	// dynamics 
	for ( i = 1 ; i <= nsweeps; i++) { 

		for (j = 0 ; j < np ; j++) {
			
			// choose particle to move 
			u = gsl_rng_uniform(ranr); 
			ip = u*np ; // auto int cast 

			// trial move 
			u = gsl_rng_uniform(ranr);
			xtrial = x[ip] + scale * ( u- 0.5 )  ; 
			u = gsl_rng_uniform(ranr);
			ytrial = y[ip] + scale * ( u-0.5 )  ; 
			u = gsl_rng_uniform(ranr);
			ztrial = z[ip] + scale * ( u- 0.5 )  ; 

			// check for collisions 
			jp = 0 ; 
			acc = 1 ; 
			while ( jp < np && acc == 1 ) { 

				if ( jp != ip ) {

					dx = xtrial - x[jp]; 
					dy = ytrial - y[jp]; 
					dz = ztrial - z[jp]; 
					
					// PBCs minimum image convention 
					dx -= L * round(dx/L); 
					dy -= L * round(dy/L); 
					dz -= L * round(dz/L); 

					dr2 = dx*dx + dy*dy + dz*dz;  
					acc = (dr2 >= d2) ? 1 : 0 ; 
				}
				jp++; 
			
			}
			if (acc == 1 ) pacc += 1 ; 

			if (acc) {
				if (ip == iwatch) { 
					dx = xtrial - x[ip] ; 
					dy = ytrial - y[ip] ; 
					dz = ztrial - z[ip] ; 
					freePath += sqrt(dx*dx + dy*dy + dz*dz) ; 
				}
				x[ip] = xtrial; 
				y[ip] = ytrial; 
				z[ip] = ztrial; 
			} else { 
				if (ip == iwatch) {
					Ncoll += 1; 
				}
			}
		
		
		}

		// collect statistics 
		if ( i > nskip) {    
	
			
		}
	//outfile << "sweep: " << i << endl; 
	} // end loop over sweeps 
	mfp = freePath / Ncoll ; 
	cout << " mean free path : " << mfp << endl; 
	
	int tdiff = time(0) - start_time ; 
	cout << "time: " << tdiff/60. << " minutes" << endl; 



	//outfile.close();
	return 0; 
}




	






