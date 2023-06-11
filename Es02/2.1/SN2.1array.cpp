#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"


using namespace std;


double errore(double a, double a2, int n) {

	if (n==0) {
		return 0.;
	} else {
		return sqrt((a2-a*a)/n);
	}

};  
 
int main (int argc, char *argv[]){

//descrivi

	int M=10000000;
	int N=100;
	int L=M/N;

	double random[M];
//x array indice [0, 1..., N-1]

	double ave[N]={0.}; //array delle medie Ai
	double av2[N]={0.}; //array delle medie al quadrato Ai2
	double sum_prog[N]={0.};
	double su2_prog[N]={0.};
	double err_prog[N]={0.};
	int x[N]={0};

	for (int i=0; i<N; i++) {

		x[i]=i;

	}
	

//CARICO M NUMERI CASUALI GENERATI UNIFORMEMENTE TRA 0 E 1 SU random

	Random rnd; //definisco un oggetto della classe random
	int seed[4]; // definisco array di 4 interi che faranno da seme; 
	int p1, p2; //due numeri primi che servono al generatore

	ifstream Primes("Primes"); //apro file di numeri primi, nomino il file Primes
	if (Primes.is_open()){
		Primes >> p1 >> p2 ; //inizializza p1 p2 con i primi due numeri del file primes
		Primes >> p1 >> p2 ;
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

//riempio array random
	for(int i=0; i<M; i++){

		random[i]=rnd.Rannyu();
	
	}


//medie dei blocchi 

	for (int i=0; i<N; i++) {

		double sum1=0.;
		double f=0.;
		double X=0.;

		for (int j=0; j<L; j++) {

			int k=j+i*L;
			X=random[k];
			f=M_PI*0.5*cos(0.5*M_PI*X);
			sum1+=f;
		}
		ave[i]=sum1/L;
		av2[i]=ave[i]*ave[i];

	}

//errori: perché la media sui blocchi e l'errore di un blocco progressivo da 0?

	for (int i=0; i<N; i++) {

		for (int j=0; j<i+1; j++) {

			
			sum_prog[i]+=ave[j]; // somma delle medie
			su2_prog[i]+=av2[j]; //somma delle medie quadratiche
		}
		sum_prog[i]/=(i+1);
		su2_prog[i]/=(i+1);
		err_prog[i]=errore(sum_prog[i], su2_prog[i],i);
	}

//fare file di OUTPUT con dati grezzi

	ofstream dati("dati2.txt");

	for (int i=0; i<N; i++) {
		
		dati << x[i] << " " << sum_prog[i] << " " << err_prog[i] << "\n"; 

		cout << i << " " << sum_prog[i] << endl;
	}
	
	dati.close();

	double aveIS[N]={0.}; //array delle medie Ai
	double av2IS[N]={0.}; //array delle medie al quadrato Ai2
	double sum_progIS[N]={0.};
	double su2_progIS[N]={0.};
	double err_progIS[N]={0.};
	double randomIS[M]={0.};
	int xIS[N]={0};

	for (int i=0; i<N; i++) {

		xIS[i]=i;

	}

	for(int i=0; i<M; i++){
//		cout << "riempio array" << endl;
		randomIS[i]=rnd.d();
	
	}

	for (int i=0; i<N; i++) {

		double sum2=0.;
		double f2=0.;
		double X2=0.;

		for (int j=0; j<L; j++) {

			int k=j+i*L;
			X2=randomIS[k];
			f2=M_PI*0.5*cos(0.5*M_PI*X2)/(2.*(1.-X2));
			sum2+=f2;
		}
		aveIS[i]=sum2/L;
		av2IS[i]=aveIS[i]*aveIS[i];

	}

//errori: perché la media sui blocchi e l'errore di un blocco progressivo da 0?

	for (int i=0; i<N; i++) {

		for (int j=0; j<i+1; j++) {

			
			sum_progIS[i]+=aveIS[j]; // somma delle medie
			su2_progIS[i]+=av2IS[j]; //somma delle medie quadratiche
		}
		sum_progIS[i]/=(i+1);
		su2_progIS[i]/=(i+1);
		err_progIS[i]=errore(sum_progIS[i], su2_progIS[i],i);
	}

//fare file di OUTPUT con dati grezzi

	ofstream datiIS("datiIS.txt");

	for (int i=0; i<N; i++) {
		
		datiIS << xIS[i] << " " << sum_progIS[i] << " " << err_progIS[i] << "\n"; 
	}
	
	datiIS.close();
//devo rimettere gli array a zero?



	return 0;
}
