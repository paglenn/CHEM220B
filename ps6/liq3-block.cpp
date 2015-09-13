//detailed balance algorithm 
#include<iostream>
#include<cmath> 
#include<fstream>
#include "structure.h"
using std::cout; 
using std::endl; 
using std::ofstream; 

int main() { 

	//system and simulation variables 
	int np; 
	double L ,rho, d ,d2,  Nside;  // box Length (volume is L^3) , density, particle diameter
	//double x[np], y[np], z[np] ;
	//double vx[np], vy[np], vz[np] ; 
	double dx, dy, dz, dr2 ; 
	int i, j, k, nsteps, ip, jp;
	int xcell, ycell, zcell ; // for PBCs
	double rc, uij, uLJ, uLJ_ref ;
	double dt, temp , tfinal; // timestep, temperature, final time (determines nsteps)  
	int samplingFreq ; 
	int iwatch ; 
	ofstream fout("rdf.dat") ; 
	ofstream coorFile("traj.xyz") ; // output file for coordinates  
	ofstream energyFile("ener.out") ; // output file for energies 
	ofstream velFile("vel.out") ; // output file for velocities 


	// Read in Parameters from param.txt 
	FILE* inFile = fopen("param.txt","r") ; 
	while(!feof(inFile)) { 
		char str[100] ; 
		fscanf(inFile,"%s",str); 
		if (*str == '#') fgets(str, 100, inFile) ; 
		else if (!strcmp(str,"numParticles")) fscanf(inFile,"%d", &np) ; 
		else if (!strcmp(str,"tfinal")) fscanf(inFile,"%lf", &tfinal) ; 
		else if (!strcmp(str,"dt")) fscanf(inFile,"%lf", &dt) ; 
		else if (!strcmp(str,"density")) fscanf(inFile,"%lf", &rho) ; 
		else if (!strcmp(str,"diameter")) fscanf(inFile,"%lf", &d) ; 
		else if (!strcmp(str,"cutoff")) fscanf(inFile,"%lf", &rc) ; // cutoff in units of diameters  
		else if (!strcmp(str,"temperature")) fscanf(inFile,"%lf", &temp) ; 

	}
	//cout << "tfinal" << tfinal << endl ; 
	//cout << "rho: " << rho << endl ; 
	//cout <<"dt" << dt << endl ; 
	//cout <<"np" << np << endl ; 
	//cout <<"cutoff" << rc << endl ; 
	//cout << " diameter " << d << endl ; 
	//cout << "temp: " << temp << endl ; 

	// system params 
	//tfinal = 5; 
	//dt = 0.001; 
	//rho = 0.8 ; 
	nsteps = ceil(tfinal/dt) ; 
	d2 = d*d ; 
	L = pow(np/rho,1./3) * d ;  
	samplingFreq = 100; 
	double endRescale = 1000 ; 
	const int numSamples = nsteps/ samplingFreq ;  
	const int N = np ; 
	vector< vector<double> > Vx(np, vector<double>());
	vector< vector<double> > Vy(np, vector<double>() ); 
	vector< vector<double> > Vz(np, vector<double>()) ; 
	vector< vector<double> > X(np, vector<double>()) ; 
	vector< vector<double> > Y(np, vector<double>()) ; 
	vector< vector<double> > Z(np, vector<double>() ); 

	Particles system(np, L, rc, temp) ; 

	for ( i = 1 ; i <= nsteps; i++) { 
		system.Advance(dt) ; 
		if (i%samplingFreq == 0 ) {
			if(i < endRescale ) {
				system.RescaleV(temp) ; 
			} else {
				coorFile << np << endl << endl ; 
				double msd = 0 ; 
				double x, y, z; 
				for(iwatch = 0 ; iwatch < np; iwatch ++ ) {
					
					vector<double> ipos = system.GetPosition(iwatch) ; 
					x =  ipos[0] - L*round(ipos[0]/L) ;
					y = ipos[1] - L*round(ipos[1]/L) ;
					z = ipos[2] - L*round(ipos[2]/L); 
					coorFile << 1 << '\t' << x << '\t' << y << '\t' << z << endl; 

					X[iwatch].push_back(ipos[0]) ; 
					Y[iwatch].push_back(ipos[1]) ; 
					Z[iwatch].push_back(ipos[2]) ; 

					vector<double> ivel = system.GetVelocity(iwatch) ; 
					velFile << ivel[0] << endl; 
					velFile << ivel[1] << endl ; 
					velFile << ivel[2] << endl ; 

					Vx[iwatch].push_back( ivel[0]) ; 
					Vy[iwatch].push_back( ivel[1]) ; 
					Vz[iwatch].push_back( ivel[2]) ; 

				}
				double E = system.GetEnergy() ; 
				double U = system.GetPotential(); 
				double K = system.GetKinetic(); 

				energyFile << (i-endRescale) *dt <<'\t' << U/np << '\t' << K/np << '\t' << E/np << endl ; 
			}
		}

	

	

	} // end loop over sweeps 

	// velocity correlation function 

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
		vcfFile << t*tau <<'\t'<< Cvv[t] << endl ; 
		msdFile << t*tau <<'\t' <<MSD[t] << endl ; 
	}

	


	// radial distribution function 
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

	return 0; 
}




	






