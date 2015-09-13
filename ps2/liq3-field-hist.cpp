//detailed balance algorithm 
#include<iostream>
#include<iomanip>
#include<cmath> 
#include<gsl/gsl_rng.h> 
#include<ctime>
#include<fstream>
#include<vector>
using std::cout; 
using std::ofstream; 
using std::endl; 
using std::vector; 
using std::setprecision ; 

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

	nsweeps = 5000000; 
	nskip = 1000 ; 
	d = 1. ; d2 = d*d ; 
	Nside = pow(np,1./3) ; 
	L = pow(2,1./3) * d * Nside;  
	scale = 0.44 ; 
	iwatch = np/2 ; 
	

	//const int nk = 5 *L / d + 1 ; 
	const int nk = 4 ; 
	double pk[nk], Rk ; 
	double ki[nk] ; 
	double totalCounts = 0 ; 
	double cx, cy, cz, sx, sy, sz ; 
	double csum, ssum; 
	vector<double> pkHist[nk]; 
	
	// for gsl rng 
	double u ; 
	const gsl_rng_type * ranT;  
	gsl_rng * ranr; 
	gsl_rng_env_setup(); 
	ranT = gsl_rng_default;
	ranr = gsl_rng_alloc (ranT);
	gsl_rng_set(ranr, time(0)); 

	ofstream outfile ; 
	outfile.open("density_hist.dat") ; 

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
	for(i = 0; i < nk; i++) {
		pk[i] = 0 ; 
		ki[i] = i * PI / d ; 

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
			for( ii = 1 ; ii < nk ; ii++) {
				jj = kk = ii ; 
				csum = ssum = 0 ; 
				for(ip = 0 ;ip < np ; ip++ ) {
					csum += cos(x[ip]*ki[ii]);
					ssum += sin(x[ip]*ki[ii]);
				}
				pk[ii] += csum + ssum ; 
				pkHist[ii].push_back(csum) ; 
				pkHist[ii].push_back(ssum) ; 
				pkHist[jj].push_back(csum) ; 
				pkHist[jj].push_back(ssum) ; 
				pkHist[kk].push_back(csum) ; 
				pkHist[kk].push_back(ssum) ; 

			}
			totalCounts += 6; 
		}
	} // end loop over sweeps 

	for( kk = 1 ; kk < nk ; kk++) { 
		pk[kk] /= totalCounts ; 
		Rk = pk[kk] / sqrt(np) ; 

	}

	double sqrtN = sqrt(np) ; 
	for(i = 0 ; i < pkHist[1].size() ; i++) {
		outfile << setprecision(5) ; 
		outfile << pkHist[1][i]/sqrtN << " " ; 
		outfile << setprecision(5) ; 
		outfile << pkHist[2][i]/sqrtN << " " ; 
		outfile << setprecision(5) ; 
		outfile << pkHist[3][i]/sqrtN << " " ; 
		outfile << endl; 

	}

	
	int tdiff = time(0) - start_time ; 
	cout << "time: " << tdiff/60. << " minutes" << endl; 
	outfile.close(); 

	return 0; 
}




	






