
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



	Random rnd; //definisco un oggetto della classe random
	rnd.Initialize(1, "Primes", "seed.in");
//struttura a blocchi come 1.1


	int M=100000;
	int N=100;
	int L=M/N;

	vector<double> posizioni(M,0.);
	vector<double> angoli(M,0.);
	vector<double> x(M, 0.);

	double d=2.;
	double l=1.;

//	int Nhits[N];
	vector<double> pi(N, 0.);
	vector<double> pi2(N, 0.);
	vector<double> sum_prog(N,0.);
	vector<double> su2_prog(N, 0.);
	vector<double> err_prog(N, 0.);
	vector<double> n(N, 0.); 

	for (int i=0; i<N; i++) {

		n[i]=i;

	}
	

//genero posizione
	for(int i=0; i<M; i++){

		posizioni[i]=rnd.Rannyu(-1., 1.)*d/2.;
	
	}

//genero array di angoli:metodo A-R per generare angolo unif in (0,pi) senza utilizzare il valore di pi

	int i=0;

	while(i<M) {
	
		double Xp=rnd.Rannyu(-1., 1.);
		double Yp=rnd.Rannyu();
		
		if (Xp*Xp+Yp*Yp<1) { 
			angoli[i]=acos(Xp/sqrt(Xp*Xp+Yp*Yp));
			i++;
		}
			
	}

//calcolo ascissa 
	for(int i=0; i<M; i++){

		x[i]=abs((l/2.)*cos(angoli[i]));
	
	}

//contatore Nhits: stima pi sui blocchi

	for(int i=0; i<N; i++){

		double Nhits=0.;

		for (int j=0; j<L; j++) {

			int k=j+i*L;

			if (posizioni[k]+x[k]>(d/2.) || posizioni[k]-x[k]<(-d/2.)) {

				Nhits=Nhits+1.;

			}
		}

		pi[i]=(2.*l*L)/(d*Nhits);
		pi2[i]=pi[i]*pi[i];
	
	}


//errori

	for (int i=0; i<N; i++) {

		for (int j=0; j<i+1; j++) {

			
			sum_prog[i]+=pi[j]; // somma delle medie
			su2_prog[i]+=pi2[j]; //somma delle medie quadratiche
		}
		sum_prog[i]/=(i+1);
		su2_prog[i]/=(i+1);
		err_prog[i]=errore(sum_prog[i], su2_prog[i],i);
	}


	ofstream pidata("pi.txt");

	for (int i=0; i<N; i++) {
		
		pidata << n[i]*L << " " << sum_prog[i] << " " << err_prog[i] << "\n"; 
	}
	
	pidata.close();




	return 0;
}



