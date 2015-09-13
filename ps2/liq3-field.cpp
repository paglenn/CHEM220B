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
	const int np = 125 ; 
	double L ,d ,d2,  Nside ;  // volume is L^3 
	double x[np], y[np], z[np] ;
	double dx, dy, dz, dr2 ; 
	double xtrial, ytrial, ztrial ; 
	int i, j, k, ii, jj, kk, ip, jp;
	int xcell, ycell, zcell ; // for PBCs
	int nsweeps, nskip; 
	double scale , acc, pacc ; 
	int iwatch; // tagged particle  

	nsweeps = 1000000; 
	nskip = 500 ; 
	d = 1. ; d2 = d*d ; 
	Nside = pow(np,1./3) ; 
	L = pow(2,1./3) * d * Nside;  
	scale = 0.44 ; 
	iwatch = np/2 ; 
	

	const int nk = 5*L / d + 2 ; 
	double pk[nk], Rk ; 
	double pk2[nk], Rk2, Sk; 
	double ki[nk] ; 
	double totalCounts = 0 ; 
	double cx, cy, cz, sx, sy, sz ; 
	double csum, ssum; 
	
	// for gsl rng 
	double u ; 
	const gsl_rng_type * ranT;  
	gsl_rng * ranr; 
	gsl_rng_env_setup(); 
	ranT = gsl_rng_default;
	ranr = gsl_rng_alloc (ranT);
	gsl_rng_set(ranr, time(0)); 

	ofstream outfile ; 
	outfile.open("density_field.dat") ; 
	
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
	for(i = 0; i < nk; i++) {
		pk[i] = 0 ; 
		pk2[i] = 0;
		ki[i] = 2*i * PI / L ; 

	}

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
				x[ip] = xtrial; 
				y[ip] = ytrial; 
				z[ip] = ztrial; 
			}
		
		
		}

		// collect statistics 
		if ( i > nskip) {    
			// check calculate FT of density field 
			for( ii = 0 ; ii < nk ; ii++) {
				csum = ssum = 0 ; 
				for(ip = 0 ; ip < np ; ip++ ) {
					csum += cos(x[ip]*ki[ii]);
					ssum += sin(x[ip]*ki[ii]);
					csum += cos(y[ip]*ki[ii]);
					ssum += sin(y[ip]*ki[ii]);
					csum += cos(z[ip]*ki[ii]);
					ssum += sin(z[ip]*ki[ii]);
				}
				pk[ii] += csum + ssum ; 
				pk2[ii] += csum *csum + ssum * ssum ; 
				//cout << csum * csum / np << endl ; 
			}
			//cout << csum << '\t\t' << ssum << endl; 

			/*
			for( jj = 0 ; jj < nk ; jj++) {
				csum = ssum = 0 ; 
				for(ip = 0 ;ip < np ; ip++ ) {
					csum += cos(y[ip]*ki[jj]);
					ssum += sin(y[ip]*ki[jj]);
				}
				pk[jj] += csum + ssum ; 
				pk2[jj] += csum*csum + ssum * ssum ; 
			}
			//cout << csum << '\t' << ssum << endl; 

			for(kk = 0; kk < nk ; kk++) { 
				csum = ssum = 0 ; 
				for(ip = 0 ;ip < np ; ip++ ) {
					csum += cos(z[ip]*ki[kk]);
					ssum += sin(z[ip]*ki[kk]);
				}
				pk[kk] += csum + ssum ; 
				pk2[kk] += csum*csum + ssum * ssum ; 
			}
			*/
			totalCounts += 6; 
			//cout << csum << '\t' << ssum << endl; 
	
		}
	} // end loop over sweeps 

	for( kk = 1 ; kk < nk ; kk++) { 
		Rk = pk[kk] / (totalCounts * sqrt(np) )  ; 
		Rk2 = pk2[kk] / (np*totalCounts) ; 
		Sk = 2 * Rk2 ; 
		outfile << ki[kk] << " " << Rk <<" "<< Sk <<  endl; 
	}

	
	int tdiff = time(0) - start_time ; 
	cout << "total counts: " << totalCounts << endl; 
	cout << "time: " << tdiff/60. << " minutes" << endl; 
	outfile.close(); 

	return 0; 
}




	






