#ifndef INIT_HEADER
#define INIT_HEADER

#ifdef GLOBAL_DATA
   #define EXTERNAL
   #define DEFAULT(X) = X
#else
#define EXTERNAL extern
   #define DEFAULT(X)
#endif

#include "global.h"

EXTERNAL MPI_Status status;
EXTERNAL int MyID, ierr, num_procs;
EXTERNAL int ipic;

EXTERNAL int NnLocalx, NnLocaly, NnLocalz, Nx;                 // Number of grid nodes for each processor
EXTERNAL int NnLocal;
EXTERNAL int Nn;                 // Number of grid nodes

void InitMeanShear();
void InitMPI(int *argc,char ***argv);
void InitCaseDim();
void InitCaseParams();
void CoorInit();
void WaveNumSetup();
void TurbFieldInit();

EXTERNAL double* ux_eq DEFAULT(NULL); 

#endif