
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
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

//descrivi

	int M=10000;
//	int n[M];
//n array indice [0, 1..., N-1] NON serve?

//	double uni[M]={0.}; //array di num gen unif
//	double exp[M]={0.}; //array di num gen exp
//	double lor[M]={0.}; //array di num gen lor
//	int n[M]={0}; //array di indice

//	for (int i=0; i<M; i++) {

//		n[i]=i;

//	}
	

//CARICO M NUMERI CASUALI GENERATI UNIFORMEMENTE TRA 0 E 1 SU random

	Random rnd; //definisco un oggetto della classe random
	int seed[4]; // definisco array di 4 interi che faranno da seme; 
	int p1, p2; //due numeri primi che servono al generatore

	ifstream Primes("Primes"); //apro file di numeri primi, nomino il file Primes
	if (Primes.is_open()){
		Primes >> p1 >> p2 ; //inizializza p1 p2 con i primi due numeri del file primes
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in"); //apro file di input (RS+4numeri), nomino il file input
	string property; //def stringa chiamata property
	if (input.is_open()){
		while ( !input.eof() ){ //fino a che non finisce il file di input
			input >> property; //in property metto la stringa che leggo da input (RANDOMSEED)
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		} //leggo i numeri di seed.in eli metto nell'array di interi seed, do il seme e i numeri P1 P2 al generatore
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	rnd.SaveSeed();

//riempio array unif
//	for(int i=0; i<M; i++){

//		uni[i]=rnd.Rannyu();
	
//	}

//riempio array exp
//	for(int i=0; i<M; i++){

//		exp[i]=rnd.Exp(1.);
	
//	}

//riempio array lorentz
//	for(int i=0; i<M; i++){

//		lor[i]=rnd.Lorentz(1.);
	
//	}


	double S1uni[M]={0.};
	double S2uni[M]={0.};
	double S10uni[M]={0.};
	double S100uni[M]={0.};

//riempio array unif
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

	double S1exp[M]={0.};
	double S2exp[M]={0.};
	double S10exp[M]={0.};
	double S100exp[M]={0.};


//riempio array unif
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


	double S1lor[M]={0.};
	double S2lor[M]={0.};
	double S10lor[M]={0.};
	double S100lor[M]={0.};


//riempio array unif
	for(int i=0; i<M; i++){

		S1lor[i]=rnd.Lorentz(1., 0.);
cout << S1lor[i] << endl;
	
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
