#include<cmath> 
#include<iostream>
#include<vector>
#include<ctime>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
using std::vector; 
using std::ofstream; 

class Particles { 

	public: 
		vector<double> R ;  // positions 
		vector<double> V ; // velocities 
		vector<double> Fi ; 
		size_t N ;  // num particles 
		double U ; // potential energy 
		double L ; // length of box dimension  
		double rc2 ; 
		double meanVel2 ; 
		
		void ComputeForces() ; 
		void InitVelocities(double temp); 
		void InitPositions(); 
		double GetPotential(); 
		double GetKinetic(); 
		double GetEnergy() ; 
		void UpdatePositions(double dt) ; 
		void UpdateVelocities(double dt) ; 
		void RescaleV(double temp)  ; 
		void Advance(double dt) ; 
		vector<double> GetPosition(int i ) ; 
		vector<double> GetVelocity(int i ) ; 
		double Temperature() ; 

	Particles(size_t np, double boxSize, double cutoff, double temp ) : N(np), L(boxSize), rc2 (cutoff*cutoff)
	{ 
		std::cout << "rc: " << cutoff << " > " << rc2 << std::endl ; 
		R = vector<double>(3*N,0.0) ; 
		V = vector<double>(3*N,0.0) ; 
		Fi = vector<double>(3*N,0.0);
		InitPositions(); 
		ComputeForces(); 
		InitVelocities(temp); 
	} 

} ; 

double Particles::GetKinetic() { 
	return 0.5 * meanVel2 *N ; 
}

double Particles::GetEnergy() {
	return GetKinetic() + U ; 
}
double Particles::Temperature() { 
	return meanVel2/3.; 
}

vector<double> Particles::GetPosition(int i ) { 
	vector<double> myvec; 
	myvec.push_back(R[3*i]) ; 
	myvec.push_back(R[3*i+1]) ; 
	myvec.push_back(R[3*i+2]) ; 
	return myvec; 
}

vector<double> Particles::GetVelocity(int i ) { 
	vector<double> myvec; 
	myvec.push_back(V[3*i]) ; 
	myvec.push_back(V[3*i+1]) ; 
	myvec.push_back(V[3*i+2]) ; 
	return myvec; 
}

void Particles::RescaleV(double temp) { 
	double alpha = sqrt(3.*temp / meanVel2) ;  
	meanVel2 = 0.0 ; 
	for (int i = 0 ; i <N ; i++) { 
		int islot = 3*i ; 
		V[islot]*= alpha ; 
		V[islot+1]*= alpha ; 
		V[islot+2]*= alpha ; 
		meanVel2 += V[islot]*V[islot]; 
		meanVel2 += V[islot+1]*V[islot+1]; 
		meanVel2 += V[islot+2]*V[islot+2]; 
	}
	meanVel2 = meanVel2 / N ; 
}

void Particles::InitPositions() { 
	double Nside = pow(N, 1./3) ; 
	int pslot ;
	int ip = 0 ; 
	for(int i = 0 ; i < Nside ; i++) { 
		for(int j = 0 ; j < Nside; j++) { 
			for(int k = 0 ; k < Nside ; k++) {
			pslot = 3 *ip ; 
			R[pslot] = i*L/Nside;
			R[pslot+1] = j*L/Nside;
			R[pslot+2] = k*L/Nside;
			ip++; 	
			}	
		}
	} 

}

void Particles::ComputeForces() { 
	int islot, jslot; 
	double fij, fij_y, fij_x, fij_z ; 
	double  dx, dy, dz, r2, _r2, _r6, _r12; 
	double uLJ, uLJ_ref, uij ; 
	double _L = 1./L ; 
	double _rc2 = 1./rc2 ; 
	double _rc6 = _rc2*_rc2*_rc2 ; 
	U = 0. ;

	uLJ_ref = 4 * _rc6*(_rc6 - 1 ) ; 

	for(int i = 0 ; i < 3*N ; i++) Fi[i] = 0.0; 

	for (int i = 0 ; i < N ; i++ ) { 
		islot = 3* i ; 
		for(int j = i+1 ; j <N ; j++ ) {
			jslot = 3*j ; 
			
			dx = R[jslot] - R[islot] ;  
			dx -= L * round(dx*_L);
			dy = R[jslot+1] - R[islot+1] ;	
			dy -= L * round(dy*_L);
			dz = R[jslot+2] - R[islot+2] ;	
			dz -= L * round(dz*_L);
			r2 = dx*dx + dy*dy + dz*dz;	
			
			if ( r2 < rc2) { 
				_r2 = 1./r2 ; 
				_r6 = _r2 * _r2 * _r2 ; 
				uLJ = 4*_r6*(1. - _r6) ; 
				uij = uLJ - uLJ_ref ;  
				fij = 24*_r2*_r6*(1.-2*_r6) ;
				// force on i from j

				// projections 
				fij_x =  fij * dx ; 
				fij_y = fij * dy  ; 
				fij_z = fij * dz  ; 

				U += uij; 
				Fi[islot] += fij_x ; 
				Fi[islot+1 ] += fij_y ; 
				Fi[islot+2] += fij_z ; 
				Fi[jslot] -= fij_x ; 
				Fi[jslot+1 ] -= fij_y ; 
				Fi[jslot+2] -= fij_z ; 

			}
		}
	}
	
}

double Particles::GetPotential() { 
	return U ; 
}

void Particles::InitVelocities(double temp) { 
	double u , vx, vy, vz ; 
	const gsl_rng_type * ranT;  
	meanVel2 = 0; 
	gsl_rng * ranr; 
	gsl_rng_env_setup(); 
	ranT = gsl_rng_default;
	ranr = gsl_rng_alloc (ranT);
	gsl_rng_set(ranr, time(0)); 

	double sigma = sqrt(temp) ; 
	double _N = 1./N ; 

	for(int i = 0 ; i < N; i++ ) { 
		vx = gsl_ran_gaussian(ranr, sigma) ; 
		vy = gsl_ran_gaussian(ranr, sigma) ; 
		vz = gsl_ran_gaussian(ranr, sigma) ; 
		int islot = 3 * i ;
		V[islot] =  vx; 
		V[islot+1] = vy; 
		V[islot+2] = vz;
		meanVel2 += (vx*vx + vy*vy + vz*vz)*_N ; 
		
	}
	//std::cout << "Start" << meanVel2 << std::endl ; 
}

void Particles::UpdatePositions(double dt) { 
	for(int i = 0 ; i < N ; i++) { 
		int islot = 3*i ; 
		R[islot] += dt * V[islot] ; 
		R[islot+1] += dt * V[islot+1] ; 
		R[islot+2] += dt * V[islot+2] ; 
	}
}

void Particles::UpdateVelocities(double dt) { 
	meanVel2 = 0.0 ; 
	//std::cout << "mv2_before: " << meanVel2 << std::endl; 
	double _N = 1./N ; 
	double Fx, Fy, Fz ;  

	//std::cout << " velocities" << std::endl ; 
	for(int i = 0 ; i < N ; i++) {
		int islot = 3*i ; 
		Fx = Fi[islot] ; 
		Fy = Fi[islot+1] ; 
		Fz = Fi[islot+2] ; 
		//std::cout << Fx << Fy << Fz << std::endl ; 

		V[islot] += dt * Fx; 
		V[islot+1] += dt * Fy; 
		V[islot+2] += dt * Fz; 
		//std::cout <<V[islot]*V[islot] << ":"; 
		//std::cout <<V[islot+1]*V[islot+1] << ":"; 
		//std::cout <<V[islot+2]*V[islot+2] ; 
		//std::cout << std::endl; 
		meanVel2 += V[islot]*V[islot]; 
		meanVel2 += V[islot+1]*V[islot+1]; 
		meanVel2 += V[islot+2]*V[islot+2]; 
		/*if(i+1 ==N ) *///std::cout << " " << meanVel2 << std::endl ; 

	}
	meanVel2 = meanVel2 * _N ; 
	//std::cout << "mv2: " << meanVel2 << std::endl; 
}

void Particles::Advance(double dt ) { 
	// update velocities 
	UpdateVelocities(dt/2.) ; 
	UpdatePositions(dt ) ; 
	ComputeForces() ; 
	UpdateVelocities(dt/2.) ; 
	//std::cout << R[3*20] << ',' << R[61] << ',' << R[62] << std::endl; 
	
}
