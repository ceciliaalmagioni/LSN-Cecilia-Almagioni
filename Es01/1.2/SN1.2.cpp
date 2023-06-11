#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"


using namespace std;

//fun errore non dovrebbe servire in 1.2
double errore(double a, double a2, int n) {

    if (n==0) {
        return 0.;
    } else {
        return sqrt((a2-a*a)/n);
    }

};  
 
int main (int argc, char *argv[]){


	int M=10000;

	//CARICO M NUMERI CASUALI GENERATI UNIFORMEMENTE TRA 0 E 1 SU random

	Random rnd; //definisco un oggetto della classe random
	rnd.Initialize(1, "Primes", "seed.in");

	vector<double> S1uni(M,0.);
	vector<double> S2uni(M,0.);
	vector<double> S10uni(M,0.);
	vector<double> S100uni(M,0.);

	//riempio vector unif
	for(int i=0; i<M; i++){
	        S1uni[i]=rnd.Rannyu();
	}

	for(int i=0; i<M; i++){
	        S2uni[i]=0.5*(rnd.Rannyu()+rnd.Rannyu());
	}

	for(int i=0; i<M; i++){
	        double sum=0.;
	        for(int j=0; j<10; j++) {
			sum+=rnd.Rannyu();
        	}
        	S10uni[i]=sum/10.;
	}

	for(int i=0; i<M; i++){
	        double sum=0.;
	        for(int j=0; j<100; j++) {
			sum+=rnd.Rannyu();
        	}
        	S100uni[i]=sum/100.;
	}

	ofstream unif("unif.txt");

	for (int i=0; i<M; i++) {
	        unif << S1uni[i] << " " << S2uni[i] << " " << S10uni[i] << " " << S100uni[i] << "\n"; 
	}
	
	unif.close();

	vector<double> S1exp(M,0.);
	vector<double> S2exp(M,0.);
	vector<double> S10exp(M,0.);
	vector<double> S100exp(M,0.);

	//riempio vector unif
	for(int i=0; i<M; i++){
	        S1exp[i]=rnd.Exp(1.);
	}

	for(int i=0; i<M; i++){
	        S2exp[i]=0.5*(rnd.Exp(1.)+rnd.Exp(1.));
	}

	for(int i=0; i<M; i++){
	        double sum=0.;
	        for(int j=0; j<10; j++) {
			sum+=rnd.Exp(1.);
        	}
		S10exp[i]=sum/10.;
	}

	for(int i=0; i<M; i++){
	        double sum=0.;
	        for(int j=0; j<100; j++) {
			sum+=rnd.Exp(1.);
        	}
        	S100exp[i]=sum/100.;
	}

	ofstream exp("exp.txt");
	for (int i=0; i<M; i++) {
	        exp << S1exp[i] << " " << S2exp[i] << " " << S10exp[i] << " " << S100exp[i] << "\n";
	}
	exp.close();

	vector<double> S1lor(M,0.);
	vector<double> S2lor(M,0.);
	vector<double> S10lor(M,0.);
	vector<double> S100lor(M,0.);

	//riempio vector unif
	for(int i=0; i<M; i++){
	        S1lor[i]=rnd.Lorentz(1., 0.);
	}

	for(int i=0; i<M; i++){
	        S2lor[i]=0.5*(rnd.Lorentz(1., 0.)+rnd.Lorentz(1., 0.));
	}

	for(int i=0; i<M; i++){
		double sum=0.;
	        for(int j=0; j<10; j++) {
	            sum+=rnd.Lorentz(1., 0.);
	        }
        S10lor[i]=sum/10.;
	}

	for(int i=0; i<M; i++){
	        double sum=0.;
	        for(int j=0; j<100; j++) {
			sum+=rnd.Lorentz(1., 0.);
        	}
        	S100lor[i]=sum/100.;
	}

	ofstream lor("lor.txt");
	for (int i=0; i<M; i++) {
	        lor << S1lor[i] << " " << S2lor[i] << " " << S10lor[i] << " " << S100lor[i] << "\n";
	}
	lor.close();

	return 0;
}






