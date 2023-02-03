#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAXSIZE 1000

int	myid, numprocs;

#define TRACE(X) fprintf(stderr,"\n[%d/%d]#%d: %s",myid,numprocs,__LINE__,X)
#define TRACEd(X,V) fprintf(stderr,"\n[%d/%d]#%d: %s%d",myid,numprocs,__LINE__,X,V)


int add(int *A, int low, int high)
{
	int res =0;
	for(int i=low; i<high; i++)
	{
		res += A[i];
	}
	
	return(res);
}



int main(int argc,char *argv[])
{
	int done = 0, n ;
	int data[MAXSIZE], i, low, high, myres, res;
	char fn[255];
	FILE *fp;
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);		
    
if(myid==0){
TRACE("starting") ;
}
	n = 0;
	while (!done)
	{
		if (myid == 0)
		{
			if (n==0) n=100; else n=0;

			/*
			strcpy(fn,getenv("HOME"));
			strcat(fn,"/lab2/rand_data.txt");
			*/
			strcpy(fn,"rand_data.txt");

TRACE("opening data file \n") ;
			/* Open Input File and Initialize Data */
			if ((fp = fopen(fn,"r")) == NULL)
			{
				printf("Can't open the input file: %s\n\n", fn);
				exit(1);
			}
			

			for(i=0; i<MAXSIZE; i++)
				fscanf(fp,"%d", &data[i]);
			
TRACE("data is read\n") ;

		} 
if(myid==0){
TRACE("calling MPI_Bcast\n") ;
}

		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(data, MAXSIZE, MPI_INT, 0, MPI_COMM_WORLD);                            /* Bcast data[] */
        
if(myid==0){
TRACE("MPI_Bcast done \n") ;
}
		
MPI_Barrier(MPI_COMM_WORLD);
		if (n == 0)
			done = 1;
		else
		{
			low = myid * 100;
			high = low + 100;
                           /* Calculate the low and high index for each processor */
if(myid==0){
TRACE("adding...\n") ;
}
			myres = add(data, low, high);                   /* Local addition for all processes */
			printf("I got %d from %d\n", myres, myid);

MPI_Barrier(MPI_COMM_WORLD);
if(myid==0){
TRACE("calling reduce...\n") ;
}
             MPI_Reduce(&myres, &res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);                   /* Global reduce */
			if (myid == 0)	
				printf("The sum is %d.\n", res);
		}
	}
	MPI_Finalize();                            /* MPI finalize */
}





