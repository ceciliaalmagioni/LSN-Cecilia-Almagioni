#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#include "random.h"


//r1 e r2G funzionano gli altri no

//fare accettanza e equilibrazione


using namespace std;


double errore(double a, double a2, int n) {

	if (n==0) {
		return 0.;
	} else {
		return sqrt((a2-a*a)/n);
	}

}; 

double psi1(vector<double> rp) {

	return pow(M_PI, -1)*exp(-2.*rp[0]);

};

double psi2(vector<double> rp) {

	return (2.0/M_PI)*(1./8.)*(1./8.)*exp(-rp[0])*cos(rp[1])*cos(rp[1])*rp[0]*rp[0];

}; 

double A1(vector<double> x, vector<double> y) {

	return min(1., psi1(x)/psi1(y));

};

double A2(vector<double> x, vector<double> y) {

	return min(1., psi2(x)/psi2(y));

};



vector<double> CarToPol(vector<double> rc) {

	vector<double> PolarCoord(3);
	
	PolarCoord[0]=sqrt(rc[0]*rc[0]+rc[1]*rc[1]+rc[2]*rc[2]);
	PolarCoord[1]=acos(rc[2]/sqrt(rc[0]*rc[0]+rc[1]*rc[1]+rc[2]*rc[2]));
	PolarCoord[2]=(signbit(rc[1]) ? -1 : 1) * acos(rc[0]/sqrt(rc[0]*rc[0] + rc[1]*rc[1])); //il primo termine restituisce il segno di y[1]

	return PolarCoord;;

}


 
int main (int argc, char *argv[]){




	Random rnd; //definisco un oggetto della classe random

	rnd.Initialize(1, "Primes", "seed.in");


	int M=1000000; //num di RW
	int N=100; //blocchi
	int L=M/N; //RW in ciascun blocco
	
	
	//Metropolis con generatore uniforme per ground state

	vector<int> x(N, 0);

	vector<vector<double>> r(M, vector<double>(3,0.));
	

	
	
	double acc=0.;
	
	r[0][0]=0.5;
	r[0][1]=1.0;
	r[0][2]=-0.75;

	for (int i=0; i<N; i++) {

		x[i]=i+1;

	}


	for (int i=0; i<1000; i++) { //equilibrazione Metropolis per 1000 step
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Rannyu(r[0][dim]-1.75, r[0][dim]+1.75);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A1(xppolar, CarToPol(r[0]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				r[0]=xp;

			
			}
	
		
	}



	for (int i=0; i<M-1; i++) {
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Rannyu(r[i][dim]-1.75, r[i][dim]+1.75);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A1(xppolar, CarToPol(r[i]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				r[i+1][0]=xp[0];
				r[i+1][1]=xp[1];
				r[i+1][2]=xp[2];
				
				acc+=+1.;
			
			} else {
			
				r[i+1][0]=r[i][0];
				r[i+1][1]=r[i][1];
				r[i+1][2]=r[i][2];			
				
				}
	
		
	}
	
	cout << "Acceptance 1 unif is " << acc/(double)M << endl;
	
//	vector<double> random(M, 0.);
	vector<double> ave(N, 0.);
	vector<double> av2(N, 0.);
	vector<double> sum_prog(N, 0.);
	vector<double> su2_prog(N, 0.);
	vector<double> err_prog(N, 0.);
	
	vector<double> r_m(M, 0.); //vettore per il valor medio del raggio
	
	for (int i=0; i<M; i++) {
	
		vector<double> pol(3, 0.);
		pol=CarToPol(r[i]);
		r_m[i] = pol[0]; 
	}


	
	for (int i = 0; i < N; i++) {
        	double sum1 = 0.;
        	for (int j = 0; j < L; j++) {
        	    int k = j + i*L;
        	    sum1 += r_m[k];
        	}
        	ave[i] = sum1 / L;
        	av2[i] = ave[i] * ave[i];
	}

    //errori: perché la media sui blocchi e l'errore di un blocco progressivo da 0?
	for (int i = 0; i < N; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_prog[i] += ave[j]; // somma delle medie
        	    su2_prog[i] += av2[j]; //somma delle medie quadratiche
        	}
        	sum_prog[i] /= (i+1);
        	su2_prog[i] /= (i+1);
        	err_prog[i] = errore(sum_prog[i], su2_prog[i], i);
    	}

    //fare file di OUTPUT con dati grezzi
	ofstream r1("r1.txt");
	
	for (int i = 0; i < N; i++) {
	        r1 << x[i]  << " " << sum_prog[i] << " " << err_prog[i] << "\n"; 
	}
	
	r1.close();
	
	ofstream PSI1("psi1.txt");
	
	for (int i = 0; i < M; i++) {
	    
	        PSI1 << r[i][0]  << " " << r[i][1] << " " << r[i][2] << "\n"; 
	}
	
	PSI1.close();
	
		
	//Metropolis con generatore gaussiano per ground state

	vector<int> xG(N, 0);

	vector<vector<double>> rG(M, vector<double>(3,0.));
	

	
	
	double accG=0.;
	
	rG[0][0]=0.5;
	rG[0][1]=1.2;
	rG[0][2]=-1.5;

	for (int i=0; i<N; i++) {

		xG[i]=i+1;

	}


	for (int i=0; i<1000; i++) { //equilibrazione Metropolis per 1000 step
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Gauss(rG[0][dim], 2.);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A1(xppolar, CarToPol(rG[0]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				rG[0][0]=xp[0];
				rG[0][1]=xp[1];
				rG[0][2]=xp[2];
				
			
			}
	
		
	}



	for (int i=0; i<M-1; i++) {
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Gauss(rG[i][dim], 1.75);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A1(xppolar, CarToPol(r[i]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				rG[i+1][0]=xp[0];
				rG[i+1][1]=xp[1];
				rG[i+1][2]=xp[2];
				
				accG+=+1.;
			
			} else {
			
				rG[i+1][0]=rG[i][0];
				rG[i+1][1]=rG[i][1];
				rG[i+1][2]=rG[i][2];			
				
				}
	
		
	}
	
	cout << "Acceptance 1 gauss is " << accG/(double)M << endl;
	
//	vector<double> random(M, 0.);
	vector<double> aveG(N, 0.);
	vector<double> av2G(N, 0.);
	vector<double> sum_progG(N, 0.);
	vector<double> su2_progG(N, 0.);
	vector<double> err_progG(N, 0.);
	
	vector<double> r_mG(M, 0.); //vettore per il valor medio del raggio
	
	for (int i=0; i<M; i++) {
	
		vector<double> pol(3, 0.);
		pol=CarToPol(rG[i]);
		r_mG[i] = pol[0]; 
	}


	
	for (int i = 0; i < N; i++) {
        	double sum1 = 0.;
        	for (int j = 0; j < L; j++) {
        	    int k = j + i*L;
        	    sum1 += r_mG[k];
        	}
        	aveG[i] = sum1 / L;
        	av2G[i] = aveG[i] * aveG[i];
	}

    //errori: perché la media sui blocchi e l'errore di un blocco progressivo da 0?
	for (int i = 0; i < N; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_progG[i] += aveG[j]; // somma delle medie
        	    su2_progG[i] += av2G[j]; //somma delle medie quadratiche
        	}
        	sum_progG[i] /= (i+1);
        	su2_progG[i] /= (i+1);
        	err_progG[i] = errore(sum_progG[i], su2_progG[i], i);
    	}

    //fare file di OUTPUT con dati grezzi
	ofstream r1G("r1G.txt");
	
	for (int i = 0; i < N; i++) {
	        r1G << xG[i]  << " " << sum_progG[i] << " " << err_progG[i] << "\n"; 
	}
	
	r1.close();
	
	ofstream PSI1G("psi1G.txt");
	
	for (int i = 0; i < M; i++) {
	    
	        PSI1G << rG[i][0]  << " " << rG[i][1] << " " << rG[i][2] << "\n"; 
	}
	
	PSI1G.close();
	
	
// Metropolis uniforme con stato eccitato



	vector<int> x2(N, 0);
	
//	vector<double> d2(M, 0.);
	vector<vector<double>> r2(M, vector<double>(3,0.));
	
	double acc2=0.;
	
	r2[0][0]=2.1;
	r2[0][1]=1.8;
	r2[0][2]=3.2;

	for (int i=0; i<N; i++) {

		x2[i]=i+1;

	}
	

	for (int i=0; i<1000; i++) { //equilibrazione Metropolis per 1000 step
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Rannyu(r2[0][dim]-1.75, r2[0][dim]+1.75);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A2(xppolar, CarToPol(r2[0]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				r2[0]=xp;

			
			}
	
		
	}

//step iesimo

	for (int i=0; i<M-1; i++) {
	
		
			vector<double> xp2(3, 0.);
			xp2[0]=rnd.Rannyu(4.0, 6.0);
			xp2[1]=rnd.Rannyu(4.0, 6.0);
			xp2[2]=rnd.Rannyu(4.0, 6.0);
			
			vector<double> xppolar2=CarToPol(xp2);
			
			double alpha2=0.;
			alpha2=A2(xppolar2, CarToPol(r2[i]));
			
			double random2=rnd.Rannyu();
			
			if (random2<=alpha2) {
			
				r2[i+1][0]=xp2[0];
				r2[i+1][1]=xp2[1];
				r2[i+1][2]=xp2[2];
				
				acc2+=1.;
			
			} else {
			
				r2[i+1][0]=r2[i][0];
				r2[i+1][1]=r2[i][1];
				r2[i+1][2]=r2[i][2];			
				
				}
	
		
	}
	cout << "Acceptance 2 unif is " << acc2/(double)M << endl;

//DATA BLOCKING
	
//	vector<double> random2(M, 0.);
	vector<double> ave2(N, 0.);
	vector<double> av22(N, 0.);
	vector<double> sum_prog2(N, 0.);
	vector<double> su2_prog2(N, 0.);
	vector<double> err_prog2(N, 0.);
	
	vector<double> r_m2(M, 0.); //vettore per il valor medio del raggio
	
	for (int i=0; i<M; i++) {
	
		vector<double> pol2(3, 0.);
		pol2=CarToPol(r2[i]);
		r_m2[i] = pol2[0]; 
	}



	
	for (int i = 0; i < N; i++) {
        	double sum2 = 0;
        	for (int j = 0; j < L; j++) {
        	    int k = j + i*L;
        	    sum2 += r_m2[k];
        	}
        	ave2[i] = sum2 / L;
        	av22[i] = ave2[i] * ave2[i];
	}

    //errori: perché la media sui blocchi e l'errore di un blocco progressivo da 0?
	for (int i = 0; i < N; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_prog2[i] += ave2[j]; // somma delle medie
        	    su2_prog2[i] += av22[j]; //somma delle medie quadratiche
        	}
        	sum_prog2[i] /= (i+1);
        	su2_prog2[i] /= (i+1);
        	err_prog2[i] = errore(sum_prog2[i], su2_prog2[i], i);
    	}

    //fare file di OUTPUT con dati grezzi
	ofstream R2("r2.txt");
	
	for (int i = 0; i < N; i++) {
	        R2 << x[i]  << " " << sum_prog2[i] << " " << err_prog2[i] << "\n"; 
	}
	
	R2.close();
	
	ofstream PSI2("psi2.txt");
	
	for (int i = 0; i < M; i++) {
	    
	        PSI2 << r2[i][0]  << " " << r2[i][1] << " " << r2[i][2] << "\n"; 
	}
	
	PSI2.close();
	
		
		
	//Metropolis con generatore gaussiano per stato eccitato

	vector<int> x2G(N, 0);

	vector<vector<double>> r2G(M, vector<double>(3,0.));
	

	
	
	double acc2G=0.;
	
	r2G[0][0]=0.5;
	r2G[0][1]=1.2;
	r2G[0][2]=-1.5;

	for (int i=0; i<N; i++) {

		x2G[i]=i+1;

	}


	for (int i=0; i<1000; i++) { //equilibrazione Metropolis per 1000 step
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Gauss(r2G[0][dim], 2.);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A2(xppolar, CarToPol(rG[0]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				r2G[0][0]=xp[0];
				r2G[0][1]=xp[1];
				r2G[0][2]=xp[2];
				
			
			}
	
		
	}



	for (int i=0; i<M-1; i++) {
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Gauss(r2G[i][dim], 1.75);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A2(xppolar, CarToPol(r2G[i]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				r2G[i+1][0]=xp[0];
				r2G[i+1][1]=xp[1];
				r2G[i+1][2]=xp[2];
				
				acc2G+=+1.;
			
			} else {
			
				r2G[i+1][0]=r2G[i][0];
				r2G[i+1][1]=r2G[i][1];
				r2G[i+1][2]=r2G[i][2];			
				
				}
	
		
	}
	
	cout << "Acceptance 2 Gauss is " << acc2G/(double)M << endl;
	
//	vector<double> random(M, 0.);
	vector<double> ave2G(N, 0.);
	vector<double> av22G(N, 0.);
	vector<double> sum_prog2G(N, 0.);
	vector<double> su2_prog2G(N, 0.);
	vector<double> err_prog2G(N, 0.);
	
	vector<double> r_m2G(M, 0.); //vettore per il valor medio del raggio
	
	for (int i=0; i<M; i++) {
	
		vector<double> pol(3, 0.);
		pol=CarToPol(r2G[i]);
		r_m2G[i] = pol[0]; 
	}


	
	for (int i = 0; i < N; i++) {
        	double sum1 = 0.;
        	for (int j = 0; j < L; j++) {
        	    int k = j + i*L;
        	    sum1 += r_m2G[k];
        	}
        	ave2G[i] = sum1 / L;
        	av22G[i] = ave2G[i] * ave2G[i];
	}

    //errori: perché la media sui blocchi e l'errore di un blocco progressivo da 0?
	for (int i = 0; i < N; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_prog2G[i] += ave2G[j]; // somma delle medie
        	    su2_prog2G[i] += av22G[j]; //somma delle medie quadratiche
        	}
        	sum_prog2G[i] /= (i+1);
        	su2_prog2G[i] /= (i+1);
        	err_prog2G[i] = errore(sum_prog2G[i], su2_prog2G[i], i);
    	}

    //fare file di OUTPUT con dati grezzi
	ofstream R2G("r2G.txt");
	
	for (int i = 0; i < N; i++) {
	        R2G << x2G[i]  << " " << sum_prog2G[i] << " " << err_prog2G[i] << "\n"; 
	}
	
	R2G.close();
	
	ofstream PSI2G("psi2G.txt");
	
	for (int i = 0; i < M; i++) {
	    
	        PSI2G << r2G[i][0]  << " " << r2G[i][1] << " " << r2G[i][2] << "\n"; 
	}
	
	PSI2G.close();



	return 0;
}


