
//cose

int M=1e5; 
int N=100; //blocchi ??
int N_beta = 10;
int L=M/N; 
int M_SA=300; 

vector<double> beta(M_SA);
vector<double> mi_sa(M_SA);
vector<double> sigma_sa(M_SA);
vector<double> energy(M_SA);
vector<double> err_energy(M_SA);
vector<double> energy_appo(M_SA);
vector<double> err_energy_appo(M_SA);

//funzioni

double errore(double, double, int);
double psi(double, double, double);
double psi2(double, double, double);
double A(double, double, double, double);
double Hpsi_over_psi(double, double, double);
void E(double, double, int);
double A_SA(double, double, double);
