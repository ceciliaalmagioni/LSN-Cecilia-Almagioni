
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



//genero cammino random discreto: credo funzioni 

	int nsteps=100; //passi di un singolo RW
	int M=1000; //num di RW
	int N=100; //blocchi
	int L=M/N; //RW in ciascun blocco



	double x=0., y=0., z=0.;
	double dx=0., dy=0., dz=0.;
	double r[M][nsteps]={0.}; //distanza dall'origine al passo i-esimo (nsteps) dell'n-esimo RW (M)




	for (int n=0; n<M; n++) { //M RW: n indicizza RW

		x=0.;	//ogni RW parte nell'origine
		y=0.;
		z=0.;	

		for (int i=0; i<nsteps; i++) { //genero l'ennesimo RW: i indicizza i passi di un RW


	
			double ran=6.*rnd.Rannyu(); 			
			int dir=static_cast<int>(floor(ran));
			switch (dir)  {
	
				case 0:
					dx=-1.; dy=0.; dz=0.;
						break;
				case 1:
					dx=1.; dy=0.; dz=0.;
						break;
				case 2:
					dx=0.; dy=-1.; dz=0.;
						break;
				case 3:
					dx=0.; dy=1.; dz=0.;
						break;
				case 4:
					dx=0.; dy=0.; dz=-1.;
						break;
				case 5:
					dx=0.; dy=0.; dz=1.;
					break;
			}
	
		x+=dx;
		y+=dy;
		z+=dz; // aggiorno la posizione a ogni passo

		r[n][i]=x*x+y*y+z*z; //distanza al quadrato dall'origine al passo iesimo dell'ennesimo RW

		}

	
	}

//N o nteps??


	for (int i=0; i<nsteps; i++) {

		ind[i]=i;

	}


//calcolo rms

	double rms[nsteps]={0.}; //rms in funzione del numero di passi rms[i] distanza media dall'origine all'iesimo passo (mediata sui blocchi)
	double errori[nsteps]={0.};

	for (int i=0; i<nsteps; i++) { // ciclo sui 100 passi: fissato il passo i calcolo media ed errore

		double ave[nsteps]={0.}; //array delle medie Ai: rms per ogni numero di passi i
		double av2[nsteps]={0.}; //array delle medie al quadrato Ai2
		double sum_prog[nsteps]={0.};
		double su2_prog[nsteps]={0.};
		double err_prog[nsteps]={0.};


		for (int j=0; j<N; j++) { //ciclo sui 100 blocchi

			double sum1=0.;

			for (int k=0; k<L; k++) { //ciclo sul j-esimo blocco

				int m=k+j*L;

				sum1+=r[m][i];

			}
		 
			ave[j]=sum1/L; //valor medio r2 del passo iesimo per il blocco j-esimo
			av2[j]=ave[j]*ave[j];
			
		}

//rms[i]=ave[j]; //? qui ok?

		for (int j=0; j<N; j++) {
	
			for (int k=0; k<j+1; k++) {
	
				
				sum_prog[j]+=ave[k]; // somma delle medie
				su2_prog[j]+=av2[k]; //somma delle medie quadratiche
			}
			sum_prog[j]/=(j+1);
			su2_prog[j]/=(j+1);
			err_prog[j]=errore(sum_prog[j], su2_prog[j],j);
		}

	rms[i]=sqrt(sum_prog[N-1]);
	errori[i]=err_prog[N-1]/(2.*sqrt(N-1)); //capire perché devo normalizzare: GPT "Inoltre, l'errore viene normalizzato dividendo per la radice quadrata del numero di blocchi N moltiplicato per un fattore 2. Questa normalizzazione viene fatta perché l'errore calcolato è l'errore sulla media di un valore stimato su N blocchi. La normalizzazione è necessaria per rendere l'errore confrontabile con il valore stimato su un singolo blocco." --> fare anche negli altri codici?
		

	}


//fare file di OUTPUT con dati grezzi

	ofstream dati("RWdisc.txt");

	for (int i=0; i<nsteps; i++) {
		
		dati << ind[i] << " " << rms[i] << " " << errori[i] << "\n"; 
	}
	
	dati.close();




	double x2=0., y2=0., z2=0.;
	double dx2=0., dy2=0., dz2=0.;
	double r2[M][nsteps]={0.}; //distanza dall'origine al passo i-esimo (nsteps) dell'n-esimo RW (M)


	for (int n=0; n<M; n++)	{

		x2=0.; //ogni RW incomincia dall'origine
		y2=0.;
		z2=0.;

		for (int i=0; i<nsteps; i++){
			
			dx2=0.;
			dy2=0.;
			dz2=0.;
	
			double theta=rnd.Rannyu(0., M_PI);
			double phi=rnd.Rannyu(0., 2.*M_PI);
	
			dx2=sin(phi)*cos(theta);
			dy2=sin(phi)*sin(theta);
			dz2=cos(phi);
	
			x2+=dx2;
			y2+=dy2;
			z2+=dz2;

			r2[n][i]=x2*x2+y2*y2+z2*z2;

	
		}
	


	}




//calcolo rms

	double rms2[nsteps]={0.}; //rms in funzione del numero di passi rms[i] distanza media dall'origine all'iesimo passo (mediata sui blocchi)
	double errori2[nsteps]={0.};

	for (int i=0; i<nsteps; i++) { // ciclo sui 100 passi: fissato il passo i calcolo media ed errore

		double ave2[nsteps]={0.}; //array delle medie Ai: rms per ogni numero di passi i
		double av22[nsteps]={0.}; //array delle medie al quadrato Ai2
		double sum_prog2[nsteps]={0.};
		double su2_prog2[nsteps]={0.};
		double err_prog2[nsteps]={0.};


		for (int j=0; j<N; j++) { //ciclo sui 100 blocchi

			double sum2=0.;

			for (int k=0; k<L; k++) { //ciclo sul j-esimo blocco

				int m=k+j*L;

				sum2+=r2[m][i];

			}
		 
			ave2[j]=sqrt(sum2/L); //valor medio r2 del passo iesimo per il blocco j-esimo
			av22[j]=ave2[j]*ave2[j];
			
		}

//rms[i]=ave[j]; //? qui ok?

		for (int j=0; j<N; j++) {
	
			for (int k=0; k<j+1; k++) {
	
				
				sum_prog2[j]+=ave2[k]; // somma delle medie
				su2_prog2[j]+=av22[k]; //somma delle medie quadratiche
			}
			sum_prog2[j]/=(j+1);
			su2_prog2[j]/=(j+1);
			err_prog2[j]=errore(sum_prog2[j], su2_prog2[j],j);
		}

	rms2[i]=sum_prog2[N-1];
	errori2[i]=err_prog2[N-1];
		

	}


	


//fare file di OUTPUT con dati grezzi

	ofstream dati2("RWcont.txt");

	for (int i=0; i<nsteps; i++) {
		
		dati2 << ind[i] << " " << rms2[i] << " " << errori2[i] << "\n"; 
	}
	
	dati2.close(); 



	return 0;
}





