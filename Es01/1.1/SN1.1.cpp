#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

double errore(double a, double a2, int n) {
    if (n == 0) {
        return 0.;
    } else {
        return sqrt((a2 - a*a) / n);
    }
}

int main (int argc, char *argv[]) {

    int M = 1000000;
    int N = 100;
    int L = M / N;

    vector<double> random(M, 0.);
    vector<double> ave(N, 0.);
    vector<double> av2(N, 0.);
    vector<double> sum_prog(N, 0.);
    vector<double> su2_prog(N, 0.);
    vector<double> err_prog(N, 0.);
    vector<int> x(N, 0);

    for (int i = 0; i < N; i++) {
        x[i] = i;
    }

    //CARICO M NUMERI CASUALI GENERATI UNIFORMEMENTE TRA 0 E 1 SU random
    Random rnd;

    rnd.Initialize(1, "Primes", "seed.in");

//-------------------------calcolo  <r>------------------------------------

    for (int i = 0; i < M; i++) {
        random[i] = rnd.Rannyu();
    }

    //medie dei blocchi 
    for (int i = 0; i < N; i++) {
        double sum1 = 0;
        for (int j = 0; j < L; j++) {
            int k = j + i*L;
            sum1 += random[k];
        }
        ave[i] = sum1 / L;
        av2[i] = ave[i] * ave[i];
    }

    //errori
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i+1; j++) {
            sum_prog[i] += ave[j]; // somma delle medie
            su2_prog[i] += av2[j]; //somma delle medie quadratiche
        }
        sum_prog[i] /= (i+1);
        su2_prog[i] /= (i+1);
        err_prog[i] = errore(sum_prog[i], su2_prog[i], i);
    }

    ofstream dati("dati1.txt");
    for (int i = 0; i < N; i++) {
        dati << x[i]  << " " << sum_prog[i] - 0.5 << " " << err_prog[i] << "\n"; 
    }
    dati.close();
    
//......................................................................................

//----------------calcolo sigma2---------------------------------------------------------

	vector<double> ave2(N, 0.);
	vector<double> av22(N, 0.);
	vector<double> sum_prog2(N, 0.);
	vector<double> su2_prog2(N, 0.);
	vector<double> err_prog2(N, 0.);

	for (int i=0; i<N; i++) {

		double sum2=0.;

		for (int j=0; j<L; j++) {

			int k=j+i*L;
			sum2+=(random[k]-0.5)*(random[k]-0.5);

		}
		
		ave2[i]=sum2/L;
		av22[i]=(ave2[i])*(ave2[i]);

	}

	for (int i=0; i<N; i++) {

		for (int j=0; j<i+1; j++) {

			sum_prog2[i]+=ave2[j];
			su2_prog2[i]+=av22[j];
		}
		
		sum_prog2[i]/=(i+1);
		su2_prog2[i]/=(i+1);
		err_prog2[i]=errore(sum_prog2[i], su2_prog2[i],i);

	}


	ofstream dati2("dati2.txt");
	for (int i=0; i<N; i++) {
		dati2 << x[i] << " " << sum_prog2[i]-1./12. << " " << err_prog2[i] << "\n"; 
	}
	dati2.close();

//.......................................................................................

//------------------------------chi-quadro----------------------------------------------


	vector<double> chi2(N, 0.); //array che contiene i 100valori del chi2

	for (int m=0; m<N; m++) { //ciclo sui 100 blocchi per calcolare i 100 chiq

		 vector<int> cont(100,0); //array contatore: immagazzino frequenze ni di chiascun intervallo

		for(int i=0; i<100; i++) { //ciclo che divide [0,1) in 100 int

			for (int j=0; j<L; j++) { //ciclo a blocchi di 1000 valori di random: incrementa contatore dell'intervallo iesimo

				int k=j+m*L;

				if (random[k]>=(double)i/100. && random[k]<(double)(i+1)/100.) {
				
					cont[i]=cont[i]+1;
					
				}
			}	
		}

		double chi_squared = 0.0;
    		for (int i = 0; i < 100; i++) {
        		double expected = (double)L / 100.0; // Valore atteso nel caso di distribuzione uniforme
       		double observed = (double)cont[i]; // Frequenza osservata nell'intervallo i-esimo
        		chi_squared += (observed - expected) * (observed - expected) / expected;
		}

    		chi2[m] = chi_squared;

	}

	ofstream chi2dat("chi2.txt");
	for (int i=0; i<N; i++) {
		chi2dat << chi2[i] << "\n"; 
	}
	chi2dat.close();

	return 0;
}





