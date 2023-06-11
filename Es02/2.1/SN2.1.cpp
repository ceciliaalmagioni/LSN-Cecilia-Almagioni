#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include <vector>

using namespace std;

double errore(double a, double a2, int n) {
    if (n==0) {
        return 0.;
    } else {
        return sqrt((a2-a*a)/n);
    }
};

int main (int argc, char *argv[]){

    int M=1000000;
    int N=100;
    int L=M/N;

    vector<double> random(M);
    vector<double> ave(N,0.); //array delle medie Ai
    vector<double> av2(N,0.); //array delle medie al quadrato Ai2
    vector<double> sum_prog(N,0.);
    vector<double> su2_prog(N,0.);
    vector<double> err_prog(N,0.);
    vector<int> x(N,0);

    for (int i=0; i<N; i++) {
        x[i]=i;
    }

    Random rnd;
    rnd.Initialize(3, "Primes", "seed.in");


    for(int i=0; i<M; i++){
        random[i]=rnd.Rannyu();
    }

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

    for (int i=0; i<N; i++) {
        for (int j=0; j<i+1; j++) {
            sum_prog[i]+=ave[j]; // somma delle medie
            su2_prog[i]+=av2[j]; //somma delle medie quadratiche
        }
        sum_prog[i]/=(i+1);
        su2_prog[i]/=(i+1);
        err_prog[i]=errore(sum_prog[i], su2_prog[i],i);
    }

    ofstream dati("dati2.txt");

    for (int i=0; i<N; i++) {
        dati << x[i] << " " << sum_prog[i] << " " << err_prog[i] << "\n";
 //       cout << i << " " << sum_prog[i] << endl;
    }
    dati.close();

    vector<double> randomIS(M);
    vector<double> aveIS(N,0.); //array delle medie Ai
    vector<double> av2IS(N,0.); //array delle medie
    vector<double> sum_progIS(N, 0.);
    vector<double> su2_progIS(N, 0.);
    vector<double> err_progIS(N, 0.);
  
    vector<int> xIS(N, 0);
;

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

//errori

	for (int i=0; i<N; i++) {

		for (int j=0; j<i+1; j++) {

			
			sum_progIS[i]+=aveIS[j]; // somma delle medie
			su2_progIS[i]+=av2IS[j]; //somma delle medie quadratiche
		}
		sum_progIS[i]/=(i+1);
		su2_progIS[i]/=(i+1);
		err_progIS[i]=errore(sum_progIS[i], su2_progIS[i],i);
	}


	ofstream datiIS("datiIS.txt");

	for (int i=0; i<N; i++) {
		
		datiIS << xIS[i] << " " << sum_progIS[i] << " " << err_progIS[i] << "\n"; 
	}
	
	datiIS.close();



	return 0;
}






