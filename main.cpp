// This code generates homogeneously sheared turbulence
// Original code by J. Schumacher
// Refactored code by P.E. Hamlington
// September, 2016



#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>


main(int argc, char **argv) 
{
	long int sum, partial_sum;
	
	MPI_Status status;
	
	int rank, root_process, ierr, num_rows, num_procs,
	 an_id, num_rows_to_receive, avg_rows_per_process, 
	 sender, num_rows_received, start_row, end_row, num_rows_to_send;

	// Initialize MPI.
	ierr = MPI_Init(&argc, &argv);

	
	/* find out MY process ID, and how many processes were started. */
	
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); 			// Determine this process's rank.
	
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  		// Determine the number of available processes.

	if ( rank == 0 ) 
		printf("Hello! My rank (my_id) is %d of %d\n", rank, num_procs);
		
	// Terminate MPI.
	ierr = MPI_Finalize();

	if ( rank == 0 ) 
        printf ( "\nNormal end of execution.\n" );

  	return 0;
}