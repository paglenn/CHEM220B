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
using std::cin ; 

//#define PI 3.14159265358979323

int main(int argc, char* argv[]) { 
	// density as input 
	if (argc < 2) return 0 ; 
	
	double _rho, rho; 
	_rho = atof(argv[1]) ; 
	//cout << " reduced density: " ; 
	//cin >> _rho ; 
	cout << " density " << _rho << endl; 

	int start_time = time(0); 
	//system and simulation variables 
	const int np = 343 ; 
	double L ,d ,d2,  Nside ;  // volume is L^3 
	double x[np], y[np], z[np] ;
	double dx, dy, dz, dr2 ; 
	double xtrial, ytrial, ztrial ; 
	int i, j, k, ii, jj, kk, ip, jp;
	int xcell, ycell, zcell ; // for PBCs
	int nsweeps, nskip; 
	double scale , acc, pacc ; 
	int iwatch; // tagged particle  

	nsweeps = 20000; 
	nskip = 500 ; 
	d = 1. ; d2 = d*d ; 
	Nside = pow(np,1./3) ; 
	L = pow(1./_rho,1./3) * d * Nside;  
	rho =   _rho/ (d*d*d) ; 
	scale = 0.44 * pow(0.5/ _rho,3); 
	iwatch = np/2 ; 
	

	//const int nk = 5 *L / d + 1 ; 
	const int nshells = 100 ; 
	double totalCounts = 0 ; 
	double Rmax = L /2.  ; 
	double dr; 
	int index; 
	int NRHist[nshells]; 
	
	// for gsl rng 
	double u ; 
	const gsl_rng_type * ranT;  
	gsl_rng * ranr; 
	gsl_rng_env_setup(); 
	ranT = gsl_rng_default;
	ranr = gsl_rng_alloc (ranT);
	gsl_rng_set(ranr, time(0)); 

	ofstream outfile ; 
	char fname[50]; 
	sprintf(fname, "density_hist-%.1f.dat",_rho) ; 
	outfile.open(fname) ; 

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

	// init hist 
	for(i = 0 ; i < nshells ; i++) {
		NRHist[i] = 0. ; 
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
			if (i == nskip +1) cout << "Hello" << endl; 
			for(ip = 0 ; ip < np; ip++ ) {
				for( jp = ip+1 ; jp < np; jp++) {
					dx = x[ip] - x[jp]; 
					dy = y[ip] - y[jp]; 
					dz = z[ip] - z[jp]; 
					
					// PBCs minimum image convention 
					dx -= L * round(dx/L); 
					dy -= L * round(dy/L); 
					dz -= L * round(dz/L); 

					dr2 = dx*dx + dy*dy + dz*dz;  
					dr = sqrt(dr2) ; 

					index = nshells * dr / Rmax ; 
					NRHist[index]++ ; 
					totalCounts++ ; 
				}
			}
		}

	} // end loop over sweeps 

	/* checksum 
	int s = 0 ; 
	for(i= 0 ; i < nshells ; i++ ) s += NRHist[i] ; 
	cout << "sum : "<< s ; 
	cout <<" expected: " << totalCounts << endl; 
	*/

	double R ,dR, P_R, g_R; 
	dR = Rmax / nshells; 
	for(i = 0 ; i < nshells ; i++) {
		//cout << dR << endl; 
		R = (i+0.5) * dR ;   
		if (NRHist[i] > 100 ){
			P_R =  NRHist[i]/ (double) totalCounts; 
			g_R = (np-1) * P_R / (4*M_PI*rho*R*R*dR)  ; 
		

			outfile << setprecision(5) << R  << " " ; 
			outfile << setprecision(5) << g_R << endl; 
			//outfile << NRHist[i] << endl; 
		}
	}
	cout << _rho << endl;  
	cout << rho << endl ; 

	
	int tdiff = time(0) - start_time ; 
	cout << "time: " << tdiff/60. << " minutes" << endl; 
	cout << " acceptance: "<< pacc / (np * nsweeps ) << endl ;  
	outfile.close(); 

	return 0; 
}




	






