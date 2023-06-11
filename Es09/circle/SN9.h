
//dimensioni

int N_cities = 34;
int dim_pop = 1000;
int N_gen = 500;

//probabilità
double pm1 = 0.3, pm2 = 0.1, pm3 = 0.1, pc = 0.55;
//double pm1 = 1., pm2 = 1., pm3 = 1., pc = 0.63;

vector<vector<double>> xy_city(N_cities, vector<double>(2));
vector<vector<int>> gen0(dim_pop, vector<int>(N_cities));
vector<vector<int>> gen1(dim_pop, vector<int>(N_cities));

vector<double> best_L(N_gen);
vector<double> besthalf_L(N_gen);


//funzioni
int pbc(int );
vector<vector<double>> Pos_circle(int);
void exchange(vector<int>&, int, int);
bool check(vector<int>);
double fitness(vector<int>);
vector<vector<int>> Generate_pop(int);
void order(vector<vector<int>>&);
void print_pop(vector<vector<int>>);
double besthalf_ave(vector<vector<int>>);
int selector(void);
void pair_permut(vector<int>&);
void m_permut(vector<int>&);
void invert(vector<int>&);
vector<int> crossover(vector<int>, vector<int>);
vector<vector<int>> new_gen (vector<vector<int>>);
void print_path(vector<int>, string);
