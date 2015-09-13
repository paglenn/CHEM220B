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
	
	double _gamma, D, beta, _beta ;  
	double sigma , var ; 
	double q, t, dt, nsteps ; 
	srand(time(0)) ; 

	ofstream fout ; 
	fout.open("data.txt") ; 

	dt = 0.001 ; 
	D = 1 ; 
	_beta = 0.4 ; 
	beta = 1./_beta; 
	_gamma = beta * D ; 
	cout << _gamma << endl ; 
	var = 2* D * dt  ; 
	sigma = sqrt(var) ;
	q = 0 ; 
	nsteps = 1e7 ; 
	
	int samplingFreq = nsteps / 100 ; 

	for(int i = 0 ; i < nsteps ; i++ ) { 
		t = i * dt ; 
		q = q + f(q) * _gamma*dt + R(sigma)  ; 

		if ( i % samplingFreq == 0 ) {  
			fout << setw(7) << t ; 
			fout << setw(15) << q ; 
			fout << endl ; 
		}
	}	

	fout.close(); 

	return 0 ; 
	


}
