#include<iostream>
#include<cmath> 
#include<fstream>
#include<iomanip>
#include<random>
#include<ctime> 

using namespace std; 

mt19937 generator; 

double w(double q ) {
	return (q*q - 1)*(q*q-1) ; 
}

double f(double q) {
	return 4*(q - q*q*q) ; 
}
double R(double sigma) {
	normal_distribution<double> gaus(0.0,sigma );
	return gaus(generator) ; 

}

int main()  {
	
	const int nbins = 20 ; 
	double _gamma, D, beta, _beta ;  
	double sigma , var ; 
	double q, t, dt, nsteps ; 
	double qa = -1 , qb = 1 ; 
	srand(time(0)) ; 

	ofstream fout ; 
	fout.open("3b.dat") ; 

	dt = 0.001 ; 
	D = 1 ; 
	_beta = 0.2 ; 
	beta = 1./_beta; 
	_gamma = beta * D ; 
	cout << _gamma << endl ; 
	var = 2* D * dt  ; 
	sigma = sqrt(var) ;
	nsteps = 1e7 ; 

	double phi ; 
	int ntraj = 1000 ;
	
	//int samplingFreq = 1000 ; 

	double q0, dq0 ; 
	double tol = 5e-4; 
	int stop ; 
	dq0  = 2.0/nbins ; 
	for(int n = 0 ; n < nbins; n++ ) { 
		phi  = 0.0 ; 
		q0 = qa + n*dq0 ; 
		int TC = 0 ; 
		cout << n << endl ; 
		


		for(int m = 0 ; m < ntraj; m++ ) {
			q = q0 ;
			t = 0 ; 
			stop = 0 ; 

			while ( stop == 0 ) { 
				t += dt ; 
				q = q + f(q) * _gamma*dt + R(sigma)  ; 
				if ( fabs(q - qa) < tol) { 
					stop = 1 ; 
					TC += 1 ;
				} else if (fabs(q-qb) < tol) {
					stop = 1 ; 
					TC+=1; 
					phi += 1 ; 
				}

			}	
		}
		phi /= (double ) TC ; 

		fout << q0 ;  
		fout << setw(15) << phi ;
		fout << endl ; 

	}
	fout.close(); 

	return 0 ; 

}
