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
	double En, Pot, Kin, Xk ; 
	double dt, temp , ti, tfinal; // timestep, temperature, final time (determines nsteps)  
	int samplingFreq ; 
	int endRescale ; 
	int iwatch ; 
	double f0, K, t; 
	int nconfig ; // num configurations to try to get a converged answer for Xbar 
	double Xb0, Xb1, Xb2; 

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
	vector<double> Xbar[3] ;
	for(int i = 0 ; i < 3;  i++) {
		Kvals[i] = 12.*(i+1)*PI/L;
		Xbar[i] = vector<double>(nsteps,0.0);
	}
	samplingFreq = 10; 
	const int numSamples = ( nsteps - endRescale) / samplingFreq ;  
	const int N = np ; 
	vector< vector<double> > Vx(np, vector<double>());
	vector< vector<double> > Vy(np, vector<double>() ); 
	vector< vector<double> > Vz(np, vector<double>()) ; 
	vector< vector<double> > X(np, vector<double>()) ; 
	vector< vector<double> > Y(np, vector<double>()) ; 
	vector< vector<double> > Z(np, vector<double>() ); 


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
				/*
				double msd = 0 ; 
				coorFile << np << endl << endl ; 
				for(iwatch = 0 ; iwatch < np; iwatch ++ ) {
					
					vector<double> ipos = system.GetPosition(iwatch) ; 
					x = ipos[0] - L*round(ipos[0]/L) ;
					y = ipos[1] - L*round(ipos[1]/L) ;
					z = ipos[2] - L*round(ipos[2]/L); 
					if (kk == 0 && cc == 0 ){
					 coorFile << 1 << setw(15) << x << '\t' << y << '\t' << z << endl; 
					}

					//X[iwatch].push_back(ipos[0]) ; 
					//Y[iwatch].push_back(ipos[1]) ; 
					//Z[iwatch].push_back(ipos[2]) ; 

					vector<double> ivel = system.GetVelocity(iwatch) ; 
					if ( kk == 0 && cc == 0) { 
					velFile << ivel[0] << endl; 
					velFile << ivel[1] << endl ; 
					velFile << ivel[2] << endl ; 
					}

					//Vx[iwatch].push_back( ivel[0]) ; 
					//Vy[iwatch].push_back( ivel[1]) ; 
					//Vz[iwatch].push_back( ivel[2]) ; 

				}
				// note: changed from t > 0 
				*/
				if ( t >= 0) {
					//cout << t << endl ; 
				
					if (kk == 0 && cc == 0 ) {
						En = system.GetEnergy() ; 
						Pot = system.GetPotential(); 
						Kin = system.GetKinetic(); 
						energyFile << t <<setw(15) << Pot/np << '\t' << Kin/np << '\t' << En/np << endl ; 
					}
					Xk = system.GetX(); 	
					Xbar[kk][step] += Xk ; 
				}
			}
		}
	} // end loop over sweeps 
	} // end loop over configs 
	//return 0 ; 
	} // end loop over k's 

	// print the X_bar's 
	for (int i = 0 ; i < nsteps ; i++) {
		step = i - endRescale; 
		t = ti + step* dt; 
		
		if (t >= 0 && i%samplingFreq == 0  ) { 
			Xb0 = Xbar[0][step] /  (double) nconfig; 
			Xb1 = Xbar[1][step] /  (double) nconfig; 
			Xb2 = Xbar[2][step] /  (double) nconfig; 
			exFile << t ; 
			exFile << setw(15) << Xb0 ; 
			exFile << setw(15) << Xb1 ; 
			exFile << setw(15) << Xb2 ; 
			exFile << endl ; 
		}
	}


	// velocity correlation function 

	/*
	char fn[100]; 
	sprintf(fn, "vcf_rho_%.2f.dat", rho ) ; 
	ofstream vcfFile(fn) ; // output file for velocities 
	
	sprintf(fn, "msd_rho_%.2f.dat", rho) ;
	ofstream msdFile(fn) ;  

	double tau = dt * samplingFreq ;  
	cout << " min time separation: " << tau << endl ; 
	const int tmax = numSamples / 5 ; 
	vector<double> MSD(tmax, 0.) ; 
	vector<double> Cvv(tmax, 0.) ; 
	vector<double> counts(tmax, 0.) ; 
	for(int s = 0 ; s < numSamples - tmax ; s++  ) { 
		for(int t = 0 ; t < tmax ; t++ ) {
			for(int ip = 0 ; ip < np ; ip++) {
				dx = X[ip][s] - X[ip][s+t] ; 
				dy = Y[ip][s] - Y[ip][s+t] ; 
				dz = Z[ip][s] - Z[ip][s+t] ; 


				MSD[t] += dx*dx + dy * dy + dz * dz ; 
				Cvv[t] += Vx[ip][s]* Vx[ip][s+t] +  Vy[ip][s]* Vy[ip][s+t] +  Vz[ip][s]* Vz[ip][s+t] ; 
				counts[t] += 1 ; 
			}
		}
	}

	for(int t = 0 ; t < tmax; t++) { 
		MSD[t] /= counts[t] ; 
		Cvv[t] /= counts[t] ; 
		vcfFile << t*tau <<setw(15)<< Cvv[t] << endl ; 
		msdFile << t*tau <<setw(15) <<MSD[t] << endl ; 
	}
	*/

	


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




	






