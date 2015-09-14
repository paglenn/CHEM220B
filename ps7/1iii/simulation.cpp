//detailed balance algorithm 
#include<iostream>
#include<cmath> 
#include<fstream>
#include "engine.h"
#include<iomanip>
using std::cout; 
using std::endl; 
using std::ofstream; 
using std::setw; 

int main() { 

	//system and simulation variables 
	int np; 
	double L ,rho, d ,d2,  Nside;  // box Length (volume is L^3) , density, particle diameter
	//double x[np], y[np], z[np] ;
	//double vx[np], vy[np], vz[np] ; 
	double dx, dy, dz, dr2 ; 
	int i, j, k, step, real_step, nsteps, ip, jp;
	int xcell, ycell, zcell ; // for PBCs
	double rc, uij, uLJ, uLJ_ref ;
	double En, Pot, Kin, Sk ; 
	double dt, temp , ti, tfinal; // timestep, temperature, final time (determines nsteps)  
	int samplingFreq ; 
	int endRescale ; 
	int iwatch ; 
	double f0, K, t; 
	int nconfig ; // num configurations to try to get a converged answer for Xbar 

	ofstream fout("rdf.dat") ; 
	ofstream coorFile("traj.xyz") ; // output file for coordinates  
	ofstream energyFile("ener.out") ; // output file for energies 
	ofstream exFile("xvals.out") ;
	ofstream velFile("vel.out") ; // output file for velocities 


	// Read in Parameters from param.txt 
	FILE* inFile = fopen("param.txt","r") ; 
	while(!feof(inFile)) { 
		char str[100] ; 
		fscanf(inFile,"%s",str); 
		if (*str == '#') fgets(str, 100, inFile) ; 
		else if (!strcmp(str,"numParticles")) fscanf(inFile,"%d", &np) ; 
		else if (!strcmp(str,"tinit")) fscanf(inFile,"%lf", &ti) ; 
		else if (!strcmp(str,"tfinal")) fscanf(inFile,"%lf", &tfinal) ; 
		else if (!strcmp(str,"dt")) fscanf(inFile,"%lf", &dt) ; 
		else if (!strcmp(str,"endRescale")) fscanf(inFile,"%d", &endRescale) ; 
		else if (!strcmp(str,"density")) fscanf(inFile,"%lf", &rho) ; 
		else if (!strcmp(str,"diameter")) fscanf(inFile,"%lf", &d) ; 
		else if (!strcmp(str,"cutoff")) fscanf(inFile,"%lf", &rc) ; // cutoff in units of diameters  
		else if (!strcmp(str,"temperature")) fscanf(inFile,"%lf", &temp) ; 
		else if (!strcmp(str,"exforce")) fscanf(inFile,"%lf", &f0) ; 
		else if (!strcmp(str,"nconfig")) fscanf(inFile,"%d", &nconfig) ; 
		
	}
	nsteps = ceil((tfinal - ti)/dt) + endRescale; 
	d2 = d*d ; 
	L = pow(np/rho,1./3) * d ;  
	double Kvals[3] ; 
	double PI = acos(-1) ; 
	vector<double> Skbar[3] ;
	const int tmax = 1./dt ; 
	double tau = dt * samplingFreq ;  
	vector<int> counts[3];
	for(int i = 0 ; i < 3;  i++) {
		Kvals[i] = 12.*(i+1)*PI/L;
		counts[i] = vector<int>(tmax,0);
	}
	samplingFreq = 20; 
	const int numSamples = ( nsteps - endRescale) / samplingFreq ;  
	const int N = np ; 
	vector< vector<double> > Vx(np, vector<double>());
	vector< vector<double> > Vy(np, vector<double>() ); 
	vector< vector<double> > Vz(np, vector<double>()) ; 
	vector< vector<double> > X(np, vector<double>()) ; 
	vector< vector<double> > Y(np, vector<double>()) ; 
	vector< vector<double> > Z(np, vector<double>() ); 
	if (numSamples < tmax ) {
		cout << "Too few samples! " << endl ; 
		return 0 ; 
	}

	ofstream SkFile ; 
	SkFile.open("Sk.dat") ; 

	// begin simulations 
	double x, y, z; 
	for(int kk = 0 ; kk < 3 ; kk++ ) {
	K = Kvals[kk];	
	cout << "k = " << K << endl ; 

	for( int cc = 0 ; cc < nconfig ; cc++ ) {
	Particles system(np, L, rc, temp, f0, K, ti, cc) ; 
	//cout << "start " << ti << endl ; 

	for ( i = 1 ; i < nsteps ; i++) { 
		step = i - endRescale; 
		t = ti + step *dt; 
		//cout << t << endl ; 
		system.Advance(t, dt) ; 
		if (i%samplingFreq == 0) {
			if(i < endRescale ) {
				system.RescaleV(temp) ; 

			} else {
				double msd = 0 ; 
				coorFile << np << endl << endl ; 
				for(iwatch = 0 ; iwatch < np; iwatch ++ ) {
					
					vector<double> ipos = system.GetPosition(iwatch) ; 
					x = ipos[0] - L*round(ipos[0]/L) ;
					y = ipos[1] - L*round(ipos[1]/L) ;
					z = ipos[2] - L*round(ipos[2]/L); 
					/*
					if (kk == 0 && cc == 0 ){
					 coorFile << 1 << setw(15) << x << '\t' << y << '\t' << z << endl; 
					}
					*/
					X[iwatch].push_back(ipos[0]) ; 
					Y[iwatch].push_back(ipos[1]) ; 
					Z[iwatch].push_back(ipos[2]) ; 

					/*
					vector<double> ivel = system.GetVelocity(iwatch) ; 
					if ( kk == 0 && cc == 0) { 
					velFile << ivel[0] << endl; 
					velFile << ivel[1] << endl ; 
					velFile << ivel[2] << endl ; 
					}

					//Vx[iwatch].push_back( ivel[0]) ; 
					//Vy[iwatch].push_back( ivel[1]) ; 
					//Vz[iwatch].push_back( ivel[2]) ; 
					*/

				}
			}
		}
	} // end loop over sweeps 
	double coskx,sinkx; 
	double cosky,sinky; 
	double coskz,sinkz; 
	cout << " min time separation: " << tau << endl ; 
	Skbar[kk] = vector<double>(tmax,0.) ; 
	for(int s = 0 ; s < numSamples - tmax ; s++  ) { 
		for(int t = 0 ; t < tmax ; t++ ) {
			for(int ip = 0 ; ip < np ; ip++) {
				for(int jp = 0 ; jp < np ; jp++) {
				// averaging starts here 


				dx = X[ip][s] - X[jp][s+t] ; 
				dy = Y[ip][s] - Y[jp][s+t] ; 
				dz = Z[ip][s] - Z[jp][s+t] ; 
					
				coskx = cos(K *dx) ; 
				//sinkx = sin(K* dx)  ; 
				cosky = cos(K * dy) ; 
				//sinky = sin(K* dy)  ; 
				coskz = cos(K * dz) ; 
				//sinkz = sin(K* dz)  ; 

				Skbar[kk][t] += coskx ; 
				Skbar[kk][t] += cosky ; 
				Skbar[kk][t] += coskz ; 
				//Skbar[kk][t] += sinkx ; 
				//Skbar[kk][t] += sinky ; 
				//Skbar[kk][t] += sinkz ; 

				//if (kk == 0 ) counts[t] += 6 ; 
				counts[kk][t] += 3 ; 
				//counts[kk][t] += 1 ; 
			}
			}
			//cout << counts[t] << endl ; 
		}
	}
	
	} // end loop over configs 

	//return 0 ; 
	} // end loop over k's 

	// velocity correlation function 

	

	for(int t = 0 ; t < tmax; t++) { 
		double tau = dt * samplingFreq ;  
		Skbar[0][t] /= (double )(2 * counts[0][t])  ; 
		Skbar[1][t] /= (double )(2 * counts[1][t]) ; 
		Skbar[2][t] /= (double )(2 * counts[2][t]) ; 
		SkFile << t*tau ;  
		SkFile << setw(15) << Skbar[0][t] ;  
		SkFile << setw(15) << Skbar[1][t] ;  
		SkFile << setw(15) << Skbar[2][t] ;  
		SkFile << endl; 
	}

	


	// radial distribution function 
	/*
	double dr; 
	for(int ip = 0 ; ip < np ; ip++) {
		for(int jp = ip+1 ; jp < np; jp++) {
			vector<double> ipos = system.GetPosition(ip) ;
			vector<double> jpos = system.GetPosition(jp) ;
			dx = ipos[0] - jpos[0] ; 
			dx -= L * round(dx/L) ; 
			dy = ipos[1] - jpos[1] ; 
			dy -= L * round(dy/L) ; 
			dz = ipos[2] - jpos[2] ; 
			dz -= L * round(dz/L) ; 
			dr = sqrt(dx*dx+dy*dy+dz*dz) ; 
			if (dr < L/2.) {
				fout << dr << endl ; 
			}
		}
	}
	*/

	return 0; 
}




	






