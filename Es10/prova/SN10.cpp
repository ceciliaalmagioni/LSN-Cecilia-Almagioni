#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
//include "random.h"
#include "mpi.h"

//#include "SN9.h"

using namespace std;

int main(int argc, char* argv[]) {

	int size, rank;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);//ottengo num tot di processi
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//ogni processo ottiene il proprio rank

/*//---------------------BROADCAST---------------------------------	
	int my_values[3];
	for(int i=0;i<3;i++) {
		if(rank==0) my_values[i]=i+1;
		else my_values[i]=0;
	} //a seconda del processo, my_values viene inizializzato in maniera differente
	
	cout<< "Prima: "<< my_values[0]<< " "<< my_values[1]<<
" "<< my_values[2]<< " per il processo "<< rank<< endl;

	MPI_Bcast(my_values,3,MPI_INTEGER,1, MPI_COMM_WORLD);
//opo la chiamata a MPI_Bcast, ogni processo avrà ricevuto i valori trasmessi e il contenuto dell'array my_values sarà lo stesso per tutti i processi: se metto 0 uguali a quelli del rank 0, se metto 1, 2, 3 quello dei rank 1 2 3
	
	cout<< "Dopo: "<< my_values[0]<< " "<< my_values[1]<< " " << my_values[2]<< " per il processo "<< rank<< endl;*/
	
//-------------------------------------------------------------------*/


/*//------------------------GATHER-------------------------------------

	if(size>3) {
		cout<<"Hai scelto troppi processi"<<endl;
	return 1;
	}

	int irecv[3]; //creo array di destinazione
	for(int i=0;i<3;i++) irecv[i]=0;
	int isend = rank + 1; //valore di invio; rank è il valore rank del processo corrente(processo 0 valore 1, processo 1 valore 2, processo 2 valore 3)

	MPI_Gather(&isend,1,MPI_INTEGER,irecv,1,MPI_INTEGER,0,MPI_COMM_WORLD); //chiamata per raccogliere i valori inviati da ogni processo nel processo root(cioè quello principale). Ogni processo invia il proprio valore (isend) e il processo principale li raccoglie nell'array irecv.(isend invia un elemento di tipo intero a irecv che attende un elemento da ogni processo di tipo intero, il processo root (riceve) è lo 0).

	if(rank==0) cout<< "irecv: " <<irecv[0] <<" "<<irecv[1] <<" " <<irecv[2] <<endl; //Infine, il processo principale (con rank 0) mostra il contenuto dell'array irecv, che conterrà tutti i valori raccolti da MPI_Gather.
	
//------------------------------------------------------------------*/


/*//---------------------------REDUCE----------------------------------

	int isend[2],irecv[2];
	for(int i=0;i<2;i++) isend[i]=rank+i+1; //{1, 2} e {2, 3} ho un isend per processo? sì
	
	cout << isend[0] << endl;

	MPI_Reduce(&isend[1],&irecv[0],1,MPI_INTEGER, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&isend[0],&irecv[1],1,MPI_INTEGER, MPI_PROD,0,MPI_COMM_WORLD);
	
	if(rank==0)cout<< "irecv somma: "<<irecv[0]<<" prodotto "<<irecv[1]<<endl;
	
//------------------------------------------------------------------*/


//--------------------COMM_SPLIT------------------------------------

	if(size!=4){cout<<"Servono 4 processi, non "<<size<<"!!"<<endl; return 1;}
	int icolor, ikey;
	if(rank==0){icolor=1;ikey=2;}
	if(rank==1){icolor=1;ikey=1;}
	if(rank==2){icolor=2;ikey=1;}
	if(rank==3){icolor=2;ikey=2;}
	//ho ottenuto due sottoprocessi (1, 0) e (2, 3)
	
	MPI_Comm nuovocom;
	MPI_Comm_split(MPI_COMM_WORLD,icolor,ikey,&nuovocom); //creazione di un nuovo comunicatore che contiene solo un sottoinsieme di processi di quello originale
	int newsize,newrank;
	MPI_Comm_size(nuovocom, &newsize);
	MPI_Comm_rank(nuovocom, &newrank);
	cout<<"Ero: "<<rank<<" di "<<size<<" ... e adesso sono:	"<< newrank<<" di "<<newsize<<endl;
	
	MPI_Finalize();
	
	return 0;
}
