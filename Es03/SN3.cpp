#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

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

	rnd.Initialize(1, "Primes", "seed.in");

	int M=100000; 
	int N=100; //blocchi
	int L=M/N; 

	vector<int> x(N, 0);

	for (int i=0; i<N; i++) {

		x[i]=i+1;

	}

//parametri del problema
	const double T=1.;
	const double K=100.;
	const double r=0.1;
	const double sigma=0.25;


//=========================EVOLUZIONE SINGOLA============================================
	
	vector<double> call(M, 0.);
	vector<double> put(M, 0.);



	for (int i=0; i<M; i++) {

		double Z=rnd.Gauss(0., 1.);
	
		double S=100.*exp((r-sigma*sigma*0.5)*T+sigma*Z*sqrt(T));
	
		call[i]=exp(-r*T)*max(0., S-K);
		put[i]=exp(-r*T)*max(0., K-S);	

//		cout << call[i] << endl;	

	}


//----------media a blocchi per call e put----------------

        vector<double> aveC(N,0.); //vector delle medie Ai: rms per ogni numero di passi i
        vector<double> av2C(N,0.); //vector delle medie al quadrato Ai2
        vector<double> sum_progC(N,0.);
        vector<double> su2_progC(N,0.);
        vector<double> err_progC(N,0.);

        vector<double> aveP(N,0.); //vector delle medie Ai: rms per ogni numero di passi i
        vector<double> av2P(N,0.); //vector delle medie al quadrato Ai2
        vector<double> sum_progP(N,0.);
        vector<double> su2_progP(N,0.);
        vector<double> err_progP(N,0.);


	for (int i=0; i<N; i++) {

		double sum1=0., sum2 =0.;

		for (int j=0; j<L; j++) {

			int k=j+i*L;
			sum1+=call[k];
			sum2+=put[k];

		}
		
		aveC[i]=sum1/L;
		av2C[i]=(aveC[i])*(aveC[i]);
		
		aveP[i]=sum2/L;
		av2P[i]=(aveP[i])*(aveP[i]);


	}

	for (int i=0; i<N; i++) {

		for (int j=0; j<i+1; j++) {

			sum_progC[i]+=aveC[j];
			su2_progC[i]+=av2C[j];
			
			sum_progP[i]+=aveP[j];
			su2_progP[i]+=av2P[j];
		}
		
		sum_progC[i]/=(i+1);
		su2_progC[i]/=(i+1);
		err_progC[i]=errore(sum_progC[i], su2_progC[i],i);
		
		sum_progP[i]/=(i+1);
		su2_progP[i]/=(i+1);
		err_progP[i]=errore(sum_progP[i], su2_progP[i],i);

	}


	ofstream call1step("call1step.txt");
	ofstream put1step("put1step.txt");
	
	for (int i=0; i<N; i++) {

		
		call1step << x[i] << " " << sum_progC[i] << " " << err_progC[i] << "\n"; 
		put1step << x[i] << " " << sum_progP[i] << " " << err_progP[i] << "\n"; 
	}
	
	call1step.close();
	put1step.close();


//=============================EVOLUZIONE N=100===========================================

	vector<double> call2(M, 0.);
	vector<double> put2(M, 0.);



	for (int i=0; i<M; i++) {

		double S=0.;

		for (int j=0; j<100; j++) {

			double Z=rnd.Gauss(0., 1.);
			S=100.*exp((r-sigma*sigma*0.5)*((double)j/100.)+sigma*Z*sqrt(T));

		}

		
	
	
		call2[i]=exp(-r*T)*max(0., S-K);
		put2[i]=exp(-r*T)*max(0., K-S);		
		

	}


//-----------------media a blocchi per call e put--------------------------------------------

        vector<double> aveC2(N,0.); //vector delle medie Ai: rms per ogni numero di passi i
        vector<double> av2C2(N,0.); //vector delle medie al quadrato Ai2
        vector<double> sum_progC2(N,0.);
        vector<double> su2_progC2(N,0.);
        vector<double> err_progC2(N,0.);

        vector<double> aveP2(N,0.); //vector delle medie Ai: rms per ogni numero di passi i
        vector<double> av2P2(N,0.); //vector delle medie al quadrato Ai2
        vector<double> sum_progP2(N,0.);
        vector<double> su2_progP2(N,0.);
        vector<double> err_progP2(N,0.);


	for (int i=0; i<N; i++) {

		double sum1=0., sum2=0.;

		for (int j=0; j<L; j++) {

			int k=j+i*L;
			sum1+=call2[k];
			sum2+=put2[k];

		}
		
		aveC2[i]=sum1/L;
		av2C2[i]=(aveC2[i])*(aveC2[i]);
		
		aveP2[i]=sum2/L;
		av2P2[i]=(aveP2[i])*(aveP2[i]);

	}

	for (int i=0; i<N; i++) {

		for (int j=0; j<i+1; j++) {

			sum_progC2[i]+=aveC2[j];
			su2_progC2[i]+=av2C2[j];
			
			sum_progP2[i]+=aveP2[j];
			su2_progP2[i]+=av2P2[j];
		}
		
		sum_progC2[i]/=(i+1);
		su2_progC2[i]/=(i+1);
		err_progC2[i]=errore(sum_progC2[i], su2_progC2[i],i);
		
		sum_progP2[i]/=(i+1);
		su2_progP2[i]/=(i+1);
		err_progP2[i]=errore(sum_progP2[i], su2_progP2[i],i);

	}


	ofstream callmulti("callmultistep.txt");
	ofstream putmulti("putmultistep.txt");

	for (int i=0; i<N; i++) {

		
		callmulti << x[i] << " " << sum_progC2[i] << " " << err_progC2[i] << "\n"; 
		putmulti << x[i] << " " << sum_progP2[i] << " " << err_progP2[i] << "\n"; 
	}
	
	callmulti.close();
	putmulti.close();


	return 0;
}





