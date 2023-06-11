#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <bits/stdc++.h>
#include "random.h"
#include "mpi.h"
//#include "SN9.h"

using namespace std;

int N_cities = 34;
int dim_pop = 1000;
int N_gen = 500;
int Nmigr = 10;

double pm1 = 0.3, pm2 = 0.1, pm3 = 0.1, pc = 0.55;

Random rnd0, rnd;

vector<vector<double>> xy_city(N_cities, vector<double>(2));
vector<vector<int>> gen0(dim_pop, vector<int>(N_cities));
vector<vector<int>> gen1(dim_pop, vector<int>(N_cities));

vector<double> best_L(N_gen);
vector<double> besthalf_L(N_gen);



int pbc(int a){
	if(a<N_cities){
	return a;
	}
	else
	 return a%N_cities + 1;
}




//funzione che genera le posizioni delle città
vector<vector<double>> Pos_circle(int N_cities) {

	vector<vector<double>> pos_city(34, vector<double>(2));
	
	//genero posizioni casuali su una circonferenza: può essere una funzione
	for (int i=0; i<N_cities; i++) {
	
		double x = rnd0.Rannyu(-1., 1.);
		double y =0.;
	
		double sgn_y =	rnd0.Rannyu(-1., 1.);
	
		if (sgn_y>0.) y = sqrt(1. - x*x);
		else  y = -sqrt(1. - x*x);
	
		pos_city[i][0] = x;
		pos_city[i][1] = y;
	
	}
	
	return pos_city;

}



//funzione che genera le posizioni delle città su un quadrato
vector<vector<double>> Pos_square(int N_cities) {

	vector<vector<double>> pos_city(N_cities, vector<double>(2));
	
	//genero posizioni casuali su una circonferenza: può essere una funzione
	for (int i=0; i<N_cities; i++) {
		
		double x = rnd0.Rannyu(-1., 1.);
		double y = rnd0.Rannyu(-1., 1.);
		
		pos_city[i][0] = x;
		pos_city[i][1] = y;		
		
		}
	
	
	return pos_city;

}




//funzione di scambio
void exchange(vector<int>&v, int i, int j) {

	int dim = v.size();
	if(i < dim and j < dim) {
	
		int appo;
		appo = v[i];
		v[i] = v[j];
		v[j] = appo;
		
	}
}



//funzione di check
bool check(vector<int> chromosome) {

	int truth_value = 1;
	if (chromosome[0] != 0) {
	truth_value = 0;
	}
	else {
		//controlla che ci siano TUTTE le città, una volta sola
		for(int i=0; i<N_cities; i++) {
			for(int j=0; j<N_cities; j++) {
				if(j!=i and chromosome[i]==chromosome[j]) {
				truth_value = 0;
				break;
				}
			}
		if (truth_value == 0) break;
		}
		
	}
	return truth_value;
}



//funzione di fitness
double fitness(vector<int> sol) {
	vector<int> sol_circ(sol.size()+1);
	double L = 0.;
	for(int i=0; i<N_cities-1; i++) {
		int n = sol[i];
		int m = sol[i+1];
		double xn = xy_city[n][0];
		double xm = xy_city[m][0];
		double yn = xy_city[n][1];
		double ym = xy_city[m][1];

		double sum = (xn-xm)*(xn-xm)+(yn-ym)*(yn-ym);
		L+= sum;
	}
	int n = sol[N_cities-1];
	int m = sol[0] ;
	double xn = xy_city[n][0];
	double xm = xy_city[m][0];
	double yn = xy_city[n][1];
	double ym = xy_city[m][1];
	
	L+=(xn-xm)*(xn-xm)+(yn-ym)*(yn-ym); //calcolo anche il tragitto di ritorno alla città di partenza
	return L;
}



//funzione che crea una popolazione
vector<vector<int>> Generate_pop(int dim_pop) {

	vector<int> sol(N_cities, 0);
	vector<vector<int>> pop (dim_pop, vector<int>(N_cities));
	

	for (int i=0; i<N_cities; i++) {
	
		sol[i] = i; //le città sono numerate da 0 a 33
	}
		
	pop[0] = sol;

	
	for (int j=1; j<dim_pop; j++) {
		
		int a = rnd.Rannyu(1, N_cities-1);
		int b = rnd.Rannyu(1, N_cities-1);
	
		exchange(sol, a, b);
	
		if(check(sol)) pop[j] = sol;
		else j = j-1;
	}
	
	return pop;

}



//funzione che ordina la pop in base al fitness
void order(vector<vector<int>> &v){
	for(int i=0; i<dim_pop; i++) {
		for(int j=i; j<dim_pop; j++) {
			if(fitness(v[j]) < fitness(v[i]))  v[i].swap(v[j]);
		}
	}
}



//funzione che stampa una popolazione

void print_pop(vector<vector<int>> gen) {
	ofstream prova("prova.txt");
	for (int i=0; i<dim_pop; i++) {
		for(int j=0; j<N_cities; j++) {
			prova << gen[i][j] << " ";
		}	
	}
	prova.close();
}



double besthalf_ave(vector<vector<int>> v) {


	int N = dim_pop/2;
	double sum = 0.;
	vector<vector<int>> v_half(N, vector<int>(N_cities)); 
	
	for(int i=0; i<N; i++) {
	
		v_half[i] = v[i];
	}
	
	for(int i=0; i<N; i++) {
	
		sum+= fitness(v_half[i]);
	}
	
	return sum/N;
}



//-------------->OPERATORI GENETICI<--------------


//operatore di selezione
int selector(void) {

	int j;
	double r = rnd.Rannyu();
	double p = 1.7;
	return j = int(dim_pop * pow(r, p));
	
}

void pair_permut(vector<int>&v) {

	int N = v.size();
	int a = rnd.Rannyu(1, N);
	int b;
	do {
	b = rnd.Rannyu(1, N);
	} while (b == a); //così sono sicura di fare lo scambio
	exchange(v, a, b);
}



void m_permut(vector<int>&v) {

	int N = v.size();
	int m = rnd.Rannyu(2, N/2-1); 
	
	int a = rnd.Rannyu(1, N);
	int b = rnd.Rannyu(1, N); 
//	cout << "a " << a<< "b " << b << "m " << m << endl;	
	for(int i=0; i<m; i++) {
		exchange(v, pbc(a+i), pbc(b+i));
	}

}//M fino a N/2-1 perchè 17 non va bene dato che se divido esattamente a metà scambio anche la prima posizione



void invert(vector<int>&v) {

	int N = v.size();
	int m = rnd.Rannyu(1, N);
	int a = rnd.Rannyu(1, N_cities);
	
//	cout << "a " << a << "m " << m << endl;
	for(int i=0; i<m/2; i++){
		exchange(v, pbc(a+i), pbc(a+m-1-i));
	}	
}



vector<int> crossover (vector<int> parent1, vector<int> parent2) {

	vector<int> son1(N_cities);
	int i_start = floor(rnd.Rannyu(1, N_cities));
	int appo_index = 0;
	
	
	for(int j=1; j<N_cities; j++) {                        
		son1[j] = parent1[j];
        }


	for(int j=1; j<N_cities; j++) {
			
		if( count(parent1.begin()+1, parent1.begin()+i_start, parent2[j]) == 0 ) { //controllo che l'elemento j di parent2 non sia nella parte fissata di parent1
			son1[i_start+appo_index] = parent2[j];
			appo_index ++;				
		}
	}
	return son1;

}



vector<vector<int>> new_gen (vector<vector<int>>old_gen) {

	vector<vector<int>> gen(dim_pop, vector<int>(N_cities));
	
	int check1 = 1, check2 = 1;
	


	for(int i=0; i<dim_pop/2; i++) {
	
		vector<int> parent1(N_cities);
		vector<int> parent2(N_cities);
		vector<int> son1(N_cities);
		vector<int> son2(N_cities);

		//scelgo il primo genitore e applico mutazioni
		int l = selector();
		parent1 = old_gen[l];
		
		
		//scelgo il secondo genitore e applico mutazioni
		int k;
		do {
		k = selector();
		} while(k == l);
		
		parent2 = old_gen[k];
		
		
		for (int m=0; m<N_cities; m++) {
			son1[m] = parent1[m];
			son2[m] = parent2[m];
		}
	
		
		
		//crossover
		
		double p = rnd.Rannyu();
		if (p < pc) {
		
			son1 = crossover(parent1, parent2);
			son2 = crossover(parent2, parent1);	
		
		}
		
	
		//mutazioni della prole
		
		double p1 = rnd.Rannyu();
		if (p1 <= pm1) pair_permut(son1);
		
		double p2 = rnd.Rannyu();
		if (p2 <= pm2) m_permut(son1);
	
		double p3 = rnd.Rannyu();
		if (p3 <= pm2) invert(son1);
		
		double p12 = rnd.Rannyu();
		if (p12 <= pm1) pair_permut(son2);
		
		double p22 = rnd.Rannyu();
		if (p22 <= pm2) m_permut(son2);
	
		double p32 = rnd.Rannyu();
		if (p32 <= pm2) invert(son2);
					
		check1 = check(son1);
		check2 = check(son2);
	
		if (check1 == 0 or check2 == 0) i = i-1;

		else{
			gen[i] = son1;
			gen[i+dim_pop/2] = son2;
		}
	}
	
		return gen;

}

void print_path(vector<int> v, string text) {

	ofstream path(text + "_path.txt");
	for(int i=0; i<N_cities; i++){
		int n = v[i];
		path << xy_city[n][0] << " " << xy_city[n][1] << "\n";
	}
//	path << xy_city[0][0] << " " << xy_city[0][1] << "\n";
	path.close();
}





int main (int argc, char *argv[]) {

	int size, rank;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);//ottengo num tot di processi
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//ogni processo ottiene il proprio rank
	MPI_Status stat1, stat2;

	rnd0.Initialize(1,"Primes", "seed.in"); //stesso seme per la generazione delle città
	rnd.Initialize(rank, "Primes", "seed.in"); //così ogni processo utilizza un seme diverso


	//genero la posizione delle città sulla circonferenza
	xy_city = Pos_square(N_cities);
			
	//creo la prima generazione
	gen0 = Generate_pop(dim_pop);
 

	//ordino la gen0 in base al fitness
	order(gen0);
	
	best_L[0] = fitness(gen0[0]);
	besthalf_L[0] = besthalf_ave(gen0);
	
	//Stampo il percorso iniziale
	string nome = "initial_square";
	print_path(gen0[dim_pop - 1], nome);
	
	vector<int>imesg (N_cities);
	vector<int>imesg2 (N_cities); 
	int itag=1; int itag2=2;
	vector <int> which_swap = {0, 1, 2, 3};
	//vector <bool> preso(4); //valore booleano per ciascun rank

	
	//creo nuove generazioni con genetica
	for(int i=1; i<N_gen; i++) { //parto da 1 perchè ho già creato la generazione 0
		vector<vector<int>> gen_old(dim_pop, vector<int>(N_cities));
		if(i == 1) gen_old = gen0;
		else gen_old = gen1;
		
		vector<vector<int>> gen_new(dim_pop, vector<int>(N_cities)); 
		gen_new = new_gen(gen_old);
		gen1 = gen_new;
		vector<vector<int>> gen_new2(dim_pop, vector<int>(N_cities));

		// Scambio dei migliori individui tra i nodi ogni Nmigr generazioni
		if (i% Nmigr == 0) {
	
			random_shuffle(which_swap.begin(), which_swap.end());
			 //MPI_Bcast(which_swap.data(), which_swap.size(), MPI_INT, 0, MPI_COMM_WORLD);
 

			//scambio i migliori delle prime due popolazioni: invece che fare due scambi non si può fare Comm_split?
			for(int j=0; j<N_cities; j++){	
			 	imesg[j] = gen1[0][j];
			 	imesg2[j] = gen1[0][j];
			}
			
			if(rank==which_swap[1]){
				MPI_Send(&imesg[0],N_cities,MPI_INTEGER,which_swap[0],itag,MPI_COMM_WORLD);
				MPI_Recv(&imesg2[0],N_cities,MPI_INTEGER,which_swap[0],itag2, MPI_COMM_WORLD,&stat2);
				//cout<<"messaggio1 = "<<imesg2[0]<<endl;
			}
			else if(rank==which_swap[0]){
				MPI_Send(&imesg2[0],N_cities,MPI_INTEGER,which_swap[1],itag2, MPI_COMM_WORLD);
				MPI_Recv(&imesg[0],N_cities,MPI_INTEGER,which_swap[1],itag, MPI_COMM_WORLD,&stat1);
				//cout<<"messaggio = "<<imesg[0]<<endl;
			}
			
			for(int j=0; j<N_cities; j++){
				if(rank==which_swap[1]){
					gen1[0][j] = imesg2[j];
				}
				else if(rank==which_swap[0]){
					gen1[0][j] = imesg[j];
				}
			}
			
			//scambio i migliori delle seconde due popolazioni
			for(int j=0; j<N_cities; j++){	
			 	imesg[j] = gen1[0][j];
			 	imesg2[j] = gen1[0][j];
			}
			 
			if(rank==which_swap[3]){
				MPI_Send(&imesg[0],N_cities,MPI_INTEGER,which_swap[2],itag,MPI_COMM_WORLD);
				MPI_Recv(&imesg2[0],N_cities,MPI_INTEGER,which_swap[2],itag2, MPI_COMM_WORLD,&stat2);
				//cout<<"messaggio1 = "<<imesg2[0]<<endl;
			}
			else if(rank==which_swap[2]){
				MPI_Send(&imesg2[0],N_cities,MPI_INTEGER,which_swap[3],itag2, MPI_COMM_WORLD);
				MPI_Recv(&imesg[0],N_cities,MPI_INTEGER,which_swap[3],itag, MPI_COMM_WORLD,&stat1);
				//cout<<"messaggio = "<<imesg[0]<<endl;
			}
			
			for(int j=0; j<N_cities; j++){
				if(rank==which_swap[3]){
					gen1[0][j] = imesg2[j];
				}
				else if(rank==which_swap[2]){
					gen1[0][j] = imesg[j];
				}
			}
		cout << "Scambiato " << which_swap[0] << " con " << which_swap[1] << ", " << which_swap[2] << " con " << which_swap[3] << endl;
        	}		
		
		order(gen1);
		
		best_L[i] = fitness(gen1[0]);
		besthalf_L[i] = besthalf_ave(gen1);
		
		cout << "numero di generazione: " << i << " per rank "<< rank << endl;

		
	}
	

	//stampo miglior L2
	ofstream out("L_rk"+to_string(rank)+".dat");
	for (int i=0; i<N_gen; i++) {
		out << i+1 << " " << best_L[i] << "\n";
	}
	out.close();

	//stampo
	ofstream out2("besthalfL_square"+to_string(rank)+".txt");	
	for (int i=0; i<N_gen; i++) {
		out2 << i+1 << " " << besthalf_L[i] << "\n";
	}
	out2.close();	
	
	
	string nome2 = "final_square"+to_string(rank);
	print_path(gen1[0], nome2);
	
	
	MPI_Finalize();

	

	return 0;
}


