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
	int my_id, root_process, ierr, i, num_rows, num_procs,
	 an_id, num_rows_to_receive, avg_rows_per_process, 
	 sender, num_rows_received, start_row, end_row, num_rows_to_send;

	/* Now replicte this process to create parallel processes.
	* From this point on, every process executes a seperate copy
	* of this program */

	ierr = MPI_Init(&argc, &argv);

	root_process = 0;

	/* find out MY process ID, and how many processes were started. */
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	if(my_id == root_process) {
		printf("hello! My id is %d of %d\n", my_id, num_procs);
	}	
	ierr = MPI_Finalize();
}