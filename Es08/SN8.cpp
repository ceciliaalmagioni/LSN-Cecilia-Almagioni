#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#include "random.h"
#include "SN8.h"

//fare accettanza e equilibrazione


using namespace std;


double errore(double a, double a2, int n) {

	if (n==0) {
		return 0.;
	} else {
		return sqrt((a2-a*a)/n);
	}

}; 

double psi(double x, double mi_psi, double sigma_psi) {

	double alpha = pow((x-mi_psi)/sigma_psi, 2)/2.;
	double beta  = pow((x+mi_psi)/sigma_psi, 2)/2.;

	return exp(-alpha) + exp(-beta);

};

//modulo di psi al quadrato 
double psi2(double x, double mi_psi2, double sigma_psi2) {

	return pow(psi(x, mi_psi2, sigma_psi2), 2);

};


//accettazione Metropolis
double A(double x, double y, double mi_A, double sigma_A) {

	return min(1., psi2(x, mi_A, sigma_A)/psi2(y, mi_A, sigma_A));

};

double Hpsi_over_psi(double x, double mi_H, double sigma_H) {

	double DDPsi_over_Psi = ( mi_H*mi_H - 2*mi_H*x*tanh(mi_H*x/(sigma_H*sigma_H)) - sigma_H*sigma_H + x*x)/(2*pow(sigma_H,4));
	return -DDPsi_over_Psi + pow(x, 4)- 5.*x*x/2.;

};


//calcolo energia con Metropolis
void E(double mi_E, double sigma_E, int index) {


	Random rand; 

	rand.Initialize(2, "Primes", "seed.in");
	
	int m=1e5; 
	int n=100; //blocchi ??
	int l=m/n; 
	
	vector<double> y(m, 0.);
	vector<double> res_integral(m, 0.);
	y[0]=0.;
	
	for (int i=0; i<m-1; i++) { //passi metropolis
	
		
		double yp=0.;
					
		yp=rand.Rannyu(y[i]-1.5, y[i]+1.5);
			
		double alfa=0.;
			
		alfa=A(yp, y[i], mi_E, sigma_E);
			
		double ran=rand.Rannyu();
			
		if (ran<=alfa) {
			
			y[i+1]=yp;
				
			//acc+=+1.;
			
		} else {
			
			y[i+1]=y[i];
				
			}
	
		res_integral[i] = Hpsi_over_psi(y[i], mi_E, sigma_E);
	
	}
	
//----------data blocking---------------------------------------
	
	vector<double> ave_E(n, 0.);
	vector<double> av2_E(n, 0.);
	vector<double> sum_prog_E(n, 0.);
	vector<double> su2_prog_E(n, 0.);
	vector<double> err_prog_E(n, 0.);
	

	
	for (int i = 0; i < n; i++) {
        	double sum = 0.;
        	for (int j = 0; j < l; j++) {
        	    int k = j + i*l;
        	    sum += res_integral[k];
        	}
        	ave_E[i] = sum / l;
        	av2_E[i] = ave_E[i] * ave_E[i];
	}

    //errori
	for (int i = 0; i < n; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_prog_E[i] += ave_E[j]; 
        	    su2_prog_E[i] += av2_E[j]; 
        	}
        	sum_prog_E[i] /= (i+1);
        	su2_prog_E[i] /= (i+1);
        	err_prog_E[i] = errore(sum_prog_E[i], su2_prog_E[i], i);
    	}

	energy_appo[index] = sum_prog_E[n-1];
	err_energy_appo[index] = err_prog_E[n-1];
	

	return;

};


//accettazione metropolis simulated annealing
double A_SA(double e_old, double e_new, double b) {

	return min(1., exp(-b*(e_new-e_old)));

};



 
int main (int argc, char *argv[]){

	Random rnd; 

	rnd.Initialize(2, "Primes", "seed.in");

//parametri per calibrazione iniziale
	/*double mi = 1.;
	double sigma = 0.5;*/

//parametri per stampare valori finali	
	double mi = 0.8;
	double sigma = -0.63;


	vector<double> x_psi(M, 0.); //vettore di posizioni 1D
	double acc=0.;
	
	x_psi[0]=0.; //dove iniziare
	
	vector<int> x(N, 0); //vettore indice
	
	for (int i=0; i<N; i++) {
	
		x[i]= i+1;
	}

//-----------Metropolis campionamento psi^2 ------------------------

	for (int i=0; i<M-1; i++) { //perché M-1?
	
		
		double xp=0.;
						
		xp=rnd.Rannyu(x_psi[i]-1.5, x_psi[i]+1.5);
			
		double alpha=0.;
			
		alpha=A(xp, x_psi[i], mi, sigma);
			
		double random=rnd.Rannyu();
			
		if (random<=alpha) {
			
			x_psi[i+1]=xp;
				
			acc+=+1.;
			
		} else {
			
			x_psi[i+1]=x_psi[i];
				
			}
	
		
	}
//----------------------------------------------------------------

//-------------Stampa valori per istogramma------------------------	

	//scommentare se mi=sigma=1
	/*ofstream pos("pos.txt");
	
	for (int i = 0; i < M; i++) {
	        pos << x_psi[i] << "\n"; 
	}
	
	pos.close();*/
	
	
	//scommentare se mi e sigma sono quelli finali
	ofstream pos2("pos2.txt");	
	for (int i = 0; i < M; i++) {
	        pos2 << x_psi[i] << "\n"; 
	}
	pos2.close();

	cout << "Acceptance 1 is " << acc/(double)M << endl;
	

	
//---------------------Calcolo dell'energia-------------------------- 	
	
//energia in ciascuna posizione
	vector<double> integral(M, 0.);
	
	for (int i=0; i<M-1; i++) {
	
		integral[i] = Hpsi_over_psi(x_psi[i], mi, sigma);
	
	}
	
//data blocking 
	
	vector<double> ave(N, 0.);
	vector<double> av2(N, 0.);
	vector<double> sum_prog(N, 0.);
	vector<double> su2_prog(N, 0.);
	vector<double> err_prog(N, 0.);
	

	
	for (int i = 0; i < N; i++) {
        	double sum1 = 0.;
        	for (int j = 0; j < L; j++) {
        	    int k = j + i*L;
        	    sum1 += integral[k];
        	}
        	ave[i] = sum1 / L;
        	av2[i] = ave[i] * ave[i];
	}

	for (int i = 0; i < N; i++) {
	       for (int j = 0; j < i+1; j++) {
        	    sum_prog[i] += ave[j]; // somma delle medie
        	    su2_prog[i] += av2[j]; //somma delle medie quadratiche
        	}
        	sum_prog[i] /= (i+1);
        	su2_prog[i] /= (i+1);
        	err_prog[i] = errore(sum_prog[i], su2_prog[i], i);
    	}
//-----------------------------------------------------------------------------------

//------------------Stampa media progressiva energia--------------------------------


	//scommenta se mi=sigma=1
/*	ofstream prova("prova.txt");	
	for (int i = 0; i < N; i++) {
	        prova << x[i]  << " " << sum_prog[i] << " " << err_prog[i] << "\n"; 
	}	
	prova.close();*/


	//scommenta se mi e sigma sono quelli finali
	ofstream fin_H("final.txt");
		for (int i = 0; i < N; i++) {
	        fin_H << x[i]  << " " << sum_prog[i] << " " << err_prog[i] << "\n"; 
	}	
	fin_H.close();
	
	
//==============================SIMULATED ANNEALING==========================================
	
	double beta_appo = 0.2; //beta di partenza

	
	//punto di partenza
	mi_sa[0] = 1.0;
	sigma_sa[0] = 1.0;
	
	//Valuto energia nel punto di partenza (la funzione E immagazzina il valore dell'energia in energy_appo)	
	E(mi_sa[0], sigma_sa[0], 0);
	energy[0] = energy_appo[0];

	
	double acc2=0.;
	
	for (int i=0; i<M_SA-1; i++) { 
	
		beta_appo+= 0.2;
		
		for(int j=0; j<N_beta; j++) {	//passi a temp fissata
		
			double mp=0., sp=0., half_step=2.5;
						
			mp=rnd.Rannyu(mi_sa[i]-half_step/(beta_appo), mi_sa[i]+half_step/(beta_appo)); // 
			sp=rnd.Rannyu(sigma_sa[i]-half_step/(beta_appo), sigma_sa[i]+half_step/(beta_appo));

			
			E(mi_sa[i], sigma_sa[i], i); // energia vecchia con i valori di mi e sigma vecchi
			double energy_old = energy_appo[i];
			double err_energy_old = err_energy_appo[i];	
			
			E(mp, sp, i);	//energia nuova con i mi e sigma estratti		
			double energy_new = energy_appo[i];
			double err_energy_new = err_energy_appo[i];			


			double alpha = A_SA(energy_old, energy_new, beta_appo); //Accettazione SA			
				
			beta[i] = (beta_appo); 
									
			if (rnd.Rannyu()<=alpha) { //accetto i valori nuovi
			
				mi_sa[i+1] = mp;
				sigma_sa[i+1] = sp;
				
				energy[i+1] = energy_new;
				err_energy[i+1] = err_energy_new;
							
				acc2+=1.;
			
			} else { //rifiuto la mossa
				
				mi_sa[i+1]=mi_sa[i];
				sigma_sa[i+1]=sigma_sa[i];
				energy[i+1] = energy_old;
				err_energy[i+1] = err_energy_old;
				}
	
		
			cout << "step " << i << " of SA completed" << endl;
		
		}
			
		//aggiorno la temperatura
		beta[i] = (beta_appo); 			
			
	
	}
	
	beta[M_SA-1] = beta[M_SA-2] +0.2; 
	beta[M_SA] = beta[M_SA-1] +0.2; 

	cout << "Acceptance SA is " << acc2/double(M_SA) << endl;

//----------------stampa i valori finali--------------------------------------------		

	ofstream mi_sigma("misigma.txt");
	for(int i=0; i<M_SA; i++) {
	
		mi_sigma << beta [i] << " "<< mi_sa[i] << " " << sigma_sa[i] << "\n";
	
	}
	mi_sigma.close();
	
	
	ofstream energySA("energy.txt");
	for(int i=0; i<M_SA; i++) {
	
		energySA << beta [i] << " " << energy[i] << " " << err_energy[i] << "\n";
	
	}
	energySA.close();
	
	

	return 0;
}


