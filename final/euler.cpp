#include<iostream>
#include<cmath> 
#include<fstream>
#include<iomanip>
#include<random>
#include<ctime> 

using namespace std; 

mt19937 generator; 

double F(double f0, double kappa  ) {
	double ans = - (1 - f0/3. ) * kappa - f0*pow(kappa, 3) / 30. ;
	return ans ; 
}

double R(double sigma) {
	normal_distribution<double> gaus(0.0,sigma );
	return gaus(generator) ; 
}

int main()  {
	
	double _gamma, D, beta, _beta ;  
	double sigma , var ; 
	double kappa, t, dt, nsteps ; 
	const int nbins = 80 ; 
	double qHist[nbins] ; 
	int bin, TC ; 

	ofstream fout ; 
	fout.open("data.txt") ; 

	dt = 0.001 ; 
	D = 1 ; 
	_beta = 0.4 ; 
	beta = 1./_beta; 
	_gamma = beta * D ; 
	var = 2* D * dt  ; 
	sigma = sqrt(var) ;
	kappa = 0.0 ; 
	nsteps = 5e6 ; 
	
	int samplingFreq = 100 ; 
	
	for(int k = 0 ; k < 12 ; k++) { 
		double f0 = k/2.0 ; 
		srand(time(0)) ; 
		vector<double> avg(2,0.0) ; 

		for(int i = 0 ; i < nsteps ; i++ ) { 

			for(int j = 0 ; j < samplingFreq; j++ ) {
				t = i * dt ; 
				kappa = kappa + F(f0, kappa)* _gamma *dt + R(sigma)  ; 
			}

			double kappa2 = kappa * kappa; 
			avg[0] += kappa2; 
			avg[1] += kappa2*kappa2 ; 
			TC += 1 ; 
		}	

		fout << setw(7) << f0 ; 
		fout << setw(15) << avg[0] / (double) TC ; 
		fout << setw(15) << avg[1] / (double) TC ; 
		fout << endl ; 
	}

	fout.close(); 
	

	return 0 ; 
	


}
