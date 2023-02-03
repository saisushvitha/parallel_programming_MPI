 

#include <stdio.h>

#include <stdlib.h>

#include <math.h>

 

/* copied from not-strictly-standard part of math.h */

#define M_PI               3.14159265358979323846

#include <mpi.h>                    /* MPI header file */

 

#define NUM_STEPS 800000000

 

/* main program */

int main(int argc, char *argv[]) {

 

    int nprocs;

    int myid;

	double start_time, end_time;
	double mypi;

	double start, end;

    double x, pi;

    double sum = 0.0;

    double step = 1.0/(double) NUM_STEPS;

 

    /* initialize for MPI */

	MPI_Init(&argc, &argv);

    /* get number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   

    /* get this process's number (ranges from 0 to nprocs - 1) */
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   

    /* record start time */
	start_time = MPI_Wtime();
       

    /* do computation */
	for(int i=myid; i < NUM_STEPS; i += nprocs) {
		x =(i+0.5)*step;
		sum += 4.0/(1.0 + x*x);
	}
 
	sum *= step;
	MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);	

    /* record end time */
	end_time = MPI_Wtime(); 
   

    /* print results */

    if (myid == 0) {

        printf("parallel program results with %d processes and %d steps:\n",

                nprocs, NUM_STEPS);

        printf("computed pi = %g  (%17.15f)\n",pi, pi);

        printf("difference between computed pi and math.h M_PI = %17.15f\n",

                fabs(pi - M_PI));

        printf("time to compute = %g seconds\n", end_time - start_time);

    }
	

	MPI_Finalize();  

    /* clean up for MPI */

 

    return 0;

}