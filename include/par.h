#ifndef PAR_HEADER
#define PAR_HEADER

#ifdef GLOBAL_DATA
   #define EXTERNAL
   #define DEFAULT(X) = X
#else
#define EXTERNAL extern
   #define DEFAULT(X)
#endif



void InitMPI(int *argc,char ***argv);
void ReadMPI(const char* fname, double* Array, int NumOfVar);

#endif