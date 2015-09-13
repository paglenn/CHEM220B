//detailed balance algorithm 
// global event chain MC 
#include<iostream>
#include<cmath> 
#include<gsl/gsl_rng.h> 
#include<ctime>
#include<fstream>
#include<cstdlib>
using std::cout; 
using std::ofstream; 
using std::endl; 
using std::cin; 

#define PI 3.14159265358979323

int main() { 
	int start_time = time(0); 
	//system and simulation variables 
	const int np = 100; 
	double L ,d ,d2,  Nside ;  // volume is L^3 
	double x[np], y[np], z[np] ;
	double dx, dy, dz, dr , lmax; 
	double delx , dely, delz, delr2; 
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
	int coll ; // collision indicator  
	FILE * trajFile ; 

	nsweeps = 300; 
	nskip = 10 ; 
	d = 1. ; d2 = d*d ; 
	L = pow(2*np,1./3) * d ;  
	cout << "Input length scale on which to move:" ; 
	cin >>scale; 
	Nside = pow(np,1./3) ; 
	prx = pry = prz = L/2.; 
	m = 4; pl = m * d ; 
	pl2 = pl*pl; 
	Nv = Nv2 = 0 ; 
	iwatch = np/2 ; 
	lmax = 5*d; 
	
	// for gsl rng 
	double u ; 
	const gsl_rng_type * ranT;  
	gsl_rng * ranr; 
	gsl_rng_env_setup(); 
	ranT = gsl_rng_default;
	ranr = gsl_rng_alloc (ranT);
	gsl_rng_set(ranr, time(0)); 


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
	trajFile = fopen("ecmc-traj.dat","w") ; 
	fprintf(trajFile,"%d %f %f %f\n",0, x[iwatch],y[iwatch],z[iwatch]) ;

	for ( i = 1 ; i <= nsweeps; i++) { 

		for (j = 0 ; j < np ; j++) {
			
			// choose particle to move 
			u = gsl_rng_uniform(ranr); 
			ip = ceil(u*np) ;  

			// first trial move 
			u = gsl_rng_uniform(ranr);
			dx = scale * ( u- 0.5 ) ;
			u = gsl_rng_uniform(ranr);
			dy = scale * ( u- 0.5 ) ; 
			u = gsl_rng_uniform(ranr);
			dz = scale * ( u- 0.5 )  ;
			dr = 0 ; 
			
			while ( dr < lmax || coll > 0 ) { 
				x[ip] += dx ; 
				y[ip] += dy; 
				z[ip] += dz ; 
				//cout<< "dr " << dr  << " lmax " << lmax << endl; 
				//cout <<"ip: " << ip << " dr: " << dr << endl; 
				dr += sqrt(dx * dx + dy * dy + dz * dz) ; 
				coll = 0 ; 

				// check for collisions 
				jp = 0 ; 
				while ( jp < np && coll == 0 ) { 

					if ( jp != ip ) {

						delx = x[ip] - x[jp]; 
						dely = y[ip] - y[jp]; 
						delz = z[ip] - z[jp]; 
						
						// PBCs minimum image convention 
						delx -= L * round(delx/L); 
						dely -= L * round(dely/L); 
						delz -= L * round(delz/L); 

						delr2 = delx*delx + dely*dely + delz*delz;  
						if (delr2 <= d2 ) {
							coll = 1 ; 
							ip = jp 	; 
						}
					}
					jp++; 
				
				}
			}
			//cout << j << endl; 
		
		}
		//cout << i << endl; 

		// collect statistics 
		if ( i > nskip) {    
			// particle trajectory 
			if (i%1 == 0 ) fprintf(trajFile,"%d %f %f %f\n",i, x[iwatch],y[iwatch],z[iwatch]) ; 

			// check probe volume 
			for( ii = 1; ii <= 1 ; ii++) { 
				for (jj = 1; jj <= 1 ; jj++ ) { 
					for ( kk = 1; kk <= 2 ; kk ++) { 

						prx = ii*L/Nside; 
						pry = jj*L/Nside; 
						prz = kk*L/Nside; 

						ip = 0;
						nv = 0 ; 
						while ( ip < np) { 

								delx = x[ip]  - prx ;
								dely = y[ip]  - pry ;
								delz = z[ip]  - prz ;
								delx -= L * round(delx/L); 
								dely -= L * round(dely/L); 
								delz -= L * round(delz/L); 
								delr2 = delx * delx + dely * dely + delz * delz ; 

								nv += (delr2 < pl2/4.) ? 1 : 0 ;  
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
	fclose(trajFile); 

	ofstream outfile ; 
	char fn[100];
	sprintf(fn,"dists_D=%.1f.dat",m) ; 
	outfile.open(fn) ; 
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
			dNv = i - Nvbar ; 
			dNv2 = dNv * dNv ; 	
			logGauss = - dNv2 / (2.*dNv2bar) - 0.5 * log(2.*PI*dNv2bar) ;
			logP = log(NvHist[i]) - logTC; 
			sprintf(data, "%d\t%f\t%f", i, logP , logGauss ) ; 
			outfile << data << endl; 
		
		}
	
	}  
	outfile.close();
	int tdiff = time(0) - start_time ; 
	cout << "time: " << tdiff/60. << " minutes" << endl; 

	return 0; 
}




	






