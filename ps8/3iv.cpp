
#include<iostream>
#include<cmath> 
#include<fstream>
#include<iomanip>
#include<random>
#include<ctime> 

using namespace std; 

double w(double q ) {
	return -2*q*q ; 
}

int main() { 

	const int npoints = 2000 ; 

	double a = -1.0 ; 
	double b = 1.0 ; 
	double dq = (b - a) / (float) npoints; 
	double _beta = 0.4; 
	double beta = 1./_beta  ; 

	double I_n[npoints], I_d ;  
	I_d = 0.0 ;   

	double q, qp, dqp ; 

	for( int i = 0 ; i < npoints ; i++)  {
		q = -1 + (i+0.5) * dq ; 
		I_d += exp(beta * w(q)) * dq ;  
		I_n[i] = 0.0 ; 

		for(int j = 0 ; j < npoints ; j++ ) {
			double bp = q ; 
			dqp = (bp -a ) / (float) npoints; 
			qp = -1 + (j+0.5) * dqp ; 
			I_n[i] += exp(beta * w(qp)) * dqp ; 
		}

	}

	ofstream fout ; 
	fout.open("3iv.dat") ; 
	for(int i = 0 ; i < npoints; i++ ) {
		
		q = -1 + (i+0.5) * dq ; 
		double phi = I_n[i] / I_d ; 
		fout << q ; 
		fout << setw(15) << phi ; 
		fout << endl ; 

	}

	fout.close(); 
}
