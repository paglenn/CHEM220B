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
	const int nbins = 80 ; 
	double qHist[nbins] ; 
	int bin, TC ; 
	srand(time(0)) ; 

	ofstream fout ; 
	fout.open("data.txt") ; 

	dt = 0.03 ; 
	D = 1 ; 
	_beta = 0.4 ; 
	beta = 1./_beta; 
	_gamma = beta * D ; 
	cout << _gamma << endl ; 
	var = 2* D * dt  ; 
	sigma = sqrt(var) ;
	q = 0 ; 
	nsteps = 1e7 ; 

	for(int i = 0 ; i < nbins ; i++ ) { 
		qHist[i] = 0.0 ;  
	}
	
	int samplingFreq = 100 ; 

	for(int i = 0 ; i < nsteps ; i++ ) { 
		t = i * dt ; 
		q = q + f(q) * _gamma*dt + R(sigma)  ; 

		if ( i % samplingFreq == 0 ) {  
			//fout << setw(7) << t ; 
			//fout << setw(15) << q ; 
			//fout << endl ; 

			bin = floor(nbins / 4. * (q + 2)) ; 
			qHist[bin] += 1 ; 
		}
	
	}	

	fout.close(); 


	// write histogram 

	fout.open("hist.dat");
	double dq = 4./nbins ; 
	double histIntegral = 0.0; 
	for(int b = 0 ; b < nbins ; b++ ) {

		histIntegral += qHist[b] ; 
	
	}
	histIntegral *= dq ;

	double _I = 1./histIntegral ; 

	for(int b = 0 ; b < nbins ; b++) { 
	
		if ( qHist[b] > 100) {
			double q_center = 4* (b + 0.5) / nbins - 2.0 ;  
			fout << q_center; 
			fout << setw(15) << qHist[b] * _I  ; 
			fout << endl ; 
		}

	}
	fout.close() ; 


	// calculate and write theoretical curve 

	fout.open("theory.dat") ; 
	double Z = 0 ; 
	int npoints = 2*nbins; 
	dq = 4./npoints; 
	for(int i = 0; i < npoints; i++ ) { 
		double qi = -2 + (i+0.5) * dq ; 
		Z += exp( - beta * w(qi) )*dq ; 
	}

	double _Z = 1./Z ; 
	for(int b = 0 ; b < nbins ; b++ ) { 
		double q_center = 4*(b + 0.5) / nbins - 2.0 ;  
		double p_q = exp(-beta * w(q_center)) * _Z; 
		fout << q_center; 
		fout << setw(15) << p_q  ; 
		fout << endl ; 
	}
	fout.close() ;


	return 0 ; 
	

}
