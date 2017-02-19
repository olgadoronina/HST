#include "all.h"
#include <mpi.h>

//======================================================================================================================
void InitMPI(int *argc, char ***argv){ // Initialization of MPI
//======================================================================================================================
	MyID=0; num_procs=1;

	if(MPI_Init(argc,argv) != MPI_SUCCESS) {
        fprintf(stderr,"MPI_Init failed \n");
        exit(0);
    }
    if(MPI_Comm_rank(MPI_COMM_WORLD,&MyID) != MPI_SUCCESS) { fprintf(stderr,"MPI_Comm_rank failed \n"); exit(0);}
    if(MPI_Comm_size(MPI_COMM_WORLD,&num_procs) != MPI_SUCCESS) { fprintf(stderr,"MPI_Comm_size failed \n"); exit(0);}
    //int ll;
	// if(MPI_Get_processor_name(proc_name, &ll) != MPI_SUCCESS)
	// { fprintf(stderr,"MPI_Get_processor_name failed \n"); exit(0);}
    //MPIinitialized = 1;
	printf("Hello! My rank (my_id) is %d of %d\n", MyID, num_procs);
}
//======================================================================================================================

//======================================================================================================================
void ReadMPI(const char* fname, double* Array, int NumOfVar){
//======================================================================================================================
	MPI_File file;
	MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	MPI_File_set_view(file, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_read(file, Array, NumOfVar, MPI_DOUBLE, &status);
	MPI_File_close(&file);
}
//======================================================================================================================