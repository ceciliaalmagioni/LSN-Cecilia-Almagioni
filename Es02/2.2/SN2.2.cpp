#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

double errore(double a, double a2, int n) {
    if (n==0) {
        return 0.;
    } else {
        return sqrt((a2-a*a)/n);
    }
};  

int main (int argc, char *argv[]) {

    Random rnd; 

    rnd.Initialize(1, "Primes", "seed.in");

 

    int nsteps=100; 
    int M=10000;
    int N=100; 
    int L=M/N;

    double x=0., y=0., z=0.;
    double dx=0., dy=0., dz=0.;

    vector<vector<double>> r(M, vector<double>(nsteps+1,0.)); 


    for (int n=0; n<M; n++) { //parto da 1 così che il passo zero coincida con l'origine
        x=0.;	
        y=0.;
        z=0.;	
        for (int i=1; i<nsteps+1; i++) {
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

	vector<int> ind(nsteps+1, 0);

	for (int i=0; i<nsteps+1; i++) {

		ind[i]=i;

	}

//calcolo rms

    vector<double> rms(nsteps+1,0.); 
    vector<double> errori(nsteps+1,0.);

    for (int i=0; i<nsteps+1; i++) { // ciclo sui 100 passi: fissato il passo i calcolo media ed errore

        vector<double> ave(nsteps+1,0.); //vector delle medie Ai: rms per ogni numero di passi i
        vector<double> av2(nsteps+1,0.); //vector delle medie al quadrato Ai2
        vector<double> sum_prog(nsteps+1,0.);
        vector<double> su2_prog(nsteps+1,0.);
        vector<double> err_prog(nsteps+1,0.);


        for (int j=0; j<N; j++) { //ciclo sui 100 blocchi

            double sum1=0.;

            for (int k=0; k<L; k++) { //ciclo sul j-esimo blocco

                int m=k+j*L;

                sum1+=r[m][i];

            }
         
            ave[j]=sqrt(sum1/L); //valor medio r2 del passo iesimo per il blocco j-esimo
            av2[j]=ave[j]*ave[j];
            
        }

        for (int j=0; j<N; j++) {
    
            for (int k=0; k<j+1; k++) {
    
                
                sum_prog[j]+=ave[k]; // somma delle medie
                su2_prog[j]+=av2[k]; //somma delle medie quadratiche
            }
            sum_prog[j]/=(j+1);
            su2_prog[j]/=(j+1);
            err_prog[j]=errore(sum_prog[j], su2_prog[j],j);
        }

        rms[i]=sum_prog[N-1];
        errori[i]=err_prog[N-1]; 
		

	}  


	ofstream dati("RWdisc.txt");

	for (int i=0; i<nsteps+1; i++) {
		
		dati << ind[i] << " " << rms[i] << " " << errori[i] << "\n"; 
	}
	
	dati.close();




	double x2=0., y2=0., z2=0.;
	double dx2=0., dy2=0., dz2=0.;
	vector<vector<double>> r2(M, vector<double>(nsteps+1,0.)); //distanza dall'origine al passo i-esimo (nsteps) dell'n-esimo RW (M)


	for (int n=0; n<M; n++)	{

		x2=0.; //ogni RW incomincia dall'origine
		y2=0.;
		z2=0.;

		for (int i=1; i<nsteps+1; i++){
			
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

	vector<double> rms2(nsteps+1, 0.); 
	vector<double> errori2(nsteps, 0.);

	for (int i=0; i<nsteps+1; i++) { // ciclo sui 100 passi: fissato il passo i calcolo media ed errore

		vector<double> ave2(nsteps+1, 0.); //array delle medie Ai: rms per ogni numero di passi i
		vector<double> av22(nsteps+1, 0.); //array delle medie al quadrato Ai2
		vector<double> sum_prog2(nsteps+1, 0.);
		vector<double> su2_prog2(nsteps+1, 0.);
		vector<double> err_prog2(nsteps+1, 0.);


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



	ofstream dati2("RWcont.txt");

	for (int i=0; i<nsteps+1; i++) {
		
		dati2 << ind[i] << " " << rms2[i] << " " << errori2[i] << "\n"; 
	}
	
	dati2.close(); 



	return 0;
}





