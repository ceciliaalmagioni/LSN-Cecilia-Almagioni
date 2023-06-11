#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#include "random.h"

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


	int M=10e6; //num di passi
	int N=100; //blocchi
	int L=M/N; //passi in ciascun blocco
	
//========================GROUND STATE=======================================	
	
//-------------------------generatore uniforme-------------------------------

	vector<int> x(N, 0);
	vector<vector<double>> r(M, vector<double>(3,0.));
	
	
	double acc=0.;
//scelgo una posizione casuale di partenza	
	r[0][0]=0.5; 
	r[0][1]=0.0; 
	r[0][2]=-0.75; 

	for (int i=0; i<N; i++) {

		x[i]=i+1;

	}
	
	

	
	ofstream eq1("eq1.txt"); 


	for (int i=0; i<1000; i++) { //equilibrazione Metropolis per 1000 step
			
			
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) { //genero posizione casuale di partnza
			
				xp[dim]=rnd.Rannyu(r[0][dim]-1.25, r[0][dim]+1.25);
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A1(xppolar, CarToPol(r[0]));
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				r[0]=xp;
			}
			
			eq1 << i+1 << " " << CarToPol(r[0])[0] << "\n";
	}
	
	eq1.close();



	for (int i=0; i<M-1; i++) {//Genero M posizioni con Metropolis
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Rannyu(r[i][dim]-1.25, r[i][dim]+1.25);			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A1(xppolar, CarToPol(r[i]));
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				r[i+1]=xp;
				
				acc+=+1.;
			
			} else {
			
				r[i+1]=r[i];
				
				}
		
	}
	
	cout << "Acceptance 1 unif is " << acc/(double)M << endl;
	

//DATA BLOCKING PER VALOR MEDIO DEL RAGGIO

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


	for (int i = 0; i < N; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_prog[i] += ave[j]; 
        	    su2_prog[i] += av2[j]; 
        	}
        	sum_prog[i] /= (i+1);
        	su2_prog[i] /= (i+1);
        	err_prog[i] = errore(sum_prog[i], su2_prog[i], i);
    	}

    //valor medio del raggio
	ofstream r1("r1.txt");
	
	for (int i = 0; i < N; i++) {
	        r1 << x[i]  << " " << sum_prog[i] << " " << err_prog[i] << "\n"; 
	}
	
	r1.close();
	
	//posizioni 3d campionate 
	ofstream PSI1("psi1FAR.txt");
	
	for (int i = 0; i < M; i++) {
	
		if(i%1000==0){
		        PSI1 << r[i][0]  << " " << r[i][1] << " " << r[i][2] << "\n";
	   /    } 
	}
	
	PSI1.close();
	
	
	
	
//---------------------generatore gaussiano------------------------------------
	
	vector<vector<double>> rG(M, vector<double>(3,0.));
	
	double accG=0.;
	
	rG[0][0]=0.5; 
	rG[0][1]=1.0;
	rG[0][2]=-0.75;

	ofstream eq1G("eq1G.txt");


	for (int i=0; i<1000; i++) { //equilibrazione Metropolis per 1000 step
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Gauss(rG[0][dim], 0.75);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.; 
			alpha=A1(xppolar, CarToPol(rG[0]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				rG[0]=xp;

			
			}
	
		eq1G << i+1 << " " << CarToPol(rG[0])[0] << "\n";
		
	}
	
	eq1G.close();



	for (int i=0; i<M-1; i++) {
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Gauss(rG[i][dim], 0.75);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A1(xppolar, CarToPol(rG[i]));
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				rG[i+1]=xp;
				
				accG+=+1.;
			
			} else {
			
				rG[i+1]=rG[i];
				
				}
	
	}
	
	cout << "Acceptance 1 gauss is " << accG/(double)M << endl;
	

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


	for (int i = 0; i < N; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_progG[i] += aveG[j]; 
        	    su2_progG[i] += av2G[j]; 
        	}
        	sum_progG[i] /= (i+1);
        	su2_progG[i] /= (i+1);
        	err_progG[i] = errore(sum_progG[i], su2_progG[i], i);
    	}


	ofstream r1G("r1G.txt");
	
	for (int i = 0; i < N; i++) {
	        r1G << x[i]  << " " << sum_progG[i] << " " << err_progG[i] << "\n"; 
	}
	
	r1G.close();
	
	ofstream PSI1G("psi1G.txt");
	
	for (int i = 0; i < M; i++) {
		if (i%1000==0){  
	  	      PSI1G << rG[i][0]  << " " << rG[i][1] << " " << rG[i][2] << "\n"; 
	        }
	}
	
	PSI1G.close();

//========================STATO ECCITATO=======================================	
	
//-------------------------generatore uniforme-------------------------------	
		



	vector<vector<double>> r2(M, vector<double>(3,0.));
	
	double acc2=0.;
	
	r2[0][0]=0.5; 
	r2[0][1]=2.0;
	r2[0][2]=-0.75;


	ofstream eq2("eq2.txt");

	for (int i=0; i<1000; i++) { //equilibrazione Metropolis per 1000 step
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Rannyu(r2[0][dim]-3.0, r2[0][dim]+3.0);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A2(xppolar, CarToPol(r2[0]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				r2[0]=xp;

			
			}
	
		eq2 << i+1 << " " << CarToPol(r2[0])[0] << "\n";
		
	}
	
	eq2.close();



	for (int i=0; i<M-1; i++) {
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Rannyu(r2[i][dim]-3.0, r2[i][dim]+3.0);

			
			}
			
			vector<double> xppolar=CarToPol(xp);
			
			double alpha=0.;
			alpha=A2(xppolar, CarToPol(r2[i]));
			
			double random=rnd.Rannyu();
			
			if (random<=alpha) {
			
				r2[i+1]=xp;
				
				acc2+=+1.;
			
			} else {
			
				r2[i+1]=r2[i];
				
				}
	
		
	}
	
	cout << "Acceptance 2 unif is " << acc2/(double)M << endl;
	

	vector<double> ave2(N, 0.);
	vector<double> av22(N, 0.);
	vector<double> sum_prog2(N, 0.);
	vector<double> su2_prog2(N, 0.);
	vector<double> err_prog2(N, 0.);
	
	vector<double> r_m2(M, 0.); //vettore per il valor medio del raggio
	
	for (int i=0; i<M; i++) {
	
		vector<double> pol(3, 0.);
		pol=CarToPol(r2[i]);
		r_m2[i] = pol[0]; 
	}


	
	for (int i = 0; i < N; i++) {
        	double sum1 = 0.;
        	for (int j = 0; j < L; j++) {
        	    int k = j + i*L;
        	    sum1 += r_m2[k];
        	}
        	ave2[i] = sum1 / L;
        	av22[i] = ave2[i] * ave2[i];
	}


	for (int i = 0; i < N; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_prog2[i] += ave2[j]; 
        	    su2_prog2[i] += av22[j]; 
        	}
        	sum_prog2[i] /= (i+1);
        	su2_prog2[i] /= (i+1);
        	err_prog2[i] = errore(sum_prog2[i], su2_prog2[i], i);
    	}


	ofstream R2("r2.txt");
	
	for (int i = 0; i < N; i++) {
	        R2 << x[i]  << " " << sum_prog2[i] << " " << err_prog2[i] << "\n"; 
	}
	
	R2.close();
	
	ofstream PSI2("psi2.txt");
	
	for (int i = 0; i < M; i++) {
		if(i%1000==0){
		        PSI2 << r2[i][0]  << " " << r2[i][1] << " " << r2[i][2] << "\n"; 
	        }
	}
	
	PSI2.close();

		
//----------------------------------generatore gaussiano------------------------------


	vector<vector<double>> r2G(M, vector<double>(3,0.));
	
	
	double acc2G=0.;
	
	r2G[0][0]=0.5;  
	r2G[0][1]=1.2;
	r2G[0][2]=-1.5;


	ofstream eq2G("eq2G.txt");

	for (int i=0; i<1000; i++) { //equilibrazione Metropolis per 1000 step
	
		
			vector<double> xp(3, 0.);
			
			for (int dim=0; dim<3; dim++) {
			
				xp[dim]=rnd.Gauss(r2G[0][dim], 1.75);

			
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
			
		eq2G << i+1 << " " << CarToPol(r2G[0])[0] << "\n";
	
		
	}
	
	eq2G.close();



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


	for (int i = 0; i < N; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_prog2G[i] += ave2G[j]; 
        	    su2_prog2G[i] += av22G[j]; 
        	}
        	sum_prog2G[i] /= (i+1);
        	su2_prog2G[i] /= (i+1);
        	err_prog2G[i] = errore(sum_prog2G[i], su2_prog2G[i], i);
    	}


	ofstream R2G("r2G.txt");
	
	for (int i = 0; i < N; i++) {
	        R2G << x[i]  << " " << sum_prog2G[i] << " " << err_prog2G[i] << "\n"; 
	}
	
	R2G.close();
	
	ofstream PSI2G("psi2G.txt");
	
	for (int i = 0; i < M; i++) {
		if(i%1000==0){
		        PSI2G << r2G[i][0]  << " " << r2G[i][1] << " " << r2G[i][2] << "\n"; 
	        }
	}
	
	PSI2G.close();



	return 0;
}


