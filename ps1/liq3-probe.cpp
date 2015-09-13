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
	double prx,pry,prz,m,pl,pl2; 
	double pv0_th; 
	int nv ; 
	double Nv, Nv2, Nvbar, Nv2bar, dNv2bar; 
	double NvHist[np+1],  totalCounts; 

	nsweeps = 10000; 
	nskip = 500 ; 
	d = 1. ; d2 = d*d ; 
	L = pow(2*np,1./3) * d ;  
	scale = 0.44 ; 
	Nside = pow(np,1./3) ; 
	prx = pry = prz = L/2.; 
	m = 4; pl = m * d ; 
	pl2 = pl*pl; 
	Nv = Nv2 = 0 ; 
	iwatch = np/2 ; 
	
	// for gsl rng 
	double u ; 
	const gsl_rng_type * ranT;  
	gsl_rng * ranr; 
	gsl_rng_env_setup(); 
	ranT = gsl_rng_default;
	ranr = gsl_rng_alloc (ranT);
	gsl_rng_set(ranr, time(0)); 

	ofstream outfile ; 
	char fn[50];
	sprintf(fn,"dists_D=%.1f.dat",m) ; 
	outfile.open(fn) ; 

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
	for(i = 0; i <= np; i++) NvHist[i] = 0 ; 
	totalCounts = 0 ; 

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
				x[ip] = xtrial; 
				y[ip] = ytrial; 
				z[ip] = ztrial; 
			}
		
		
		}

		// collect statistics 
		if ( i > nskip) {    
			// particle trajectory 
			//if (i%10 == 0 ) printf("%d %f %f %f\n",i, x[iwatch],y[iwatch],z[iwatch]) ; 

			// check probe volume 
			for( ii = 1; ii <= 35 ; ii++) { 
				for (jj = 1; jj <= 35 ; jj++ ) { 
					for ( kk = 1; kk <= 35 ; kk ++) { 

						prx = ii*L/Nside; 
						pry = jj*L/Nside; 
						prz = kk*L/Nside; 

						ip = 0;
						nv = 0 ; 
						while ( ip < np) { 

								dx = x[ip]  - prx ;
								dy = y[ip]  - pry ;
								dz = z[ip]  - prz ;
								dx -= L * round(dx/L); 
								dy -= L * round(dy/L); 
								dz -= L * round(dz/L); 
								dr2 = dx * dx + dy * dy + dz * dz ; 

								nv += (dr2 < pl2/4.) ? 1 : 0 ;  
								ip++ ; 
						}	
						Nv += nv; 
						Nv2 += nv*nv; 
						NvHist[nv]++; 
						totalCounts++; 
					}
				}
			}
			
	
		}
	//outfile << "sweep: " << i << endl; 
	} // end loop over sweeps 
	pacc /= ( (double) nsweeps*np) ; 
	Nvbar  = Nv / (double) totalCounts ;
	Nv2bar = Nv2 / (double) totalCounts ; 
	dNv2bar = Nv2bar - Nvbar * Nvbar ; 
	outfile << "<Nv> :\t" << Nvbar << endl; 
	outfile << "<dNv2> :\t" << dNv2bar << endl ; 
	//outfile << nsweeps << " sweeps "<< endl ; 
	//outfile << "acceptance ratio: " << pacc << endl ;  

	// 
	double dNv,dNv2, logP, logGauss ; 
	double logTC = log(totalCounts) ;
	for (i = 0; i <= np; i++) { 
		char data[80]; 
		if(NvHist[i] > 100 ) {
			dNv = i - Nvbar ; dNv2 = dNv * dNv ; 	
			logGauss = - dNv2 / (2.*dNv2bar) - 0.5 * log(2.*PI*dNv2bar) ;
			logP = log(NvHist[i]) - logTC; 
			sprintf(data, "%d\t%f\t%f\n", i, logP , logGauss ) ; 
			outfile << data ; 
		
		}
	
	}  

	int tdiff = time(0) - start_time ; 
	cout << "time: " << tdiff/60. << " minutes" << endl; 



	outfile.close();
	return 0; 
}




	






