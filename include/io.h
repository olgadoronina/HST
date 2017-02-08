#ifndef IO_HEADER
#define IO_HEADER

#ifdef GLOBAL_DATA
   #define EXTERNAL
   #define DEFAULT(X) = X
#else
#define EXTERNAL extern
   #define DEFAULT(X)
#endif

#include "global.h"
// for system (mkdir)
#include <sys/stat.h>
#include <sys/types.h>


void MakeDir(string path);
void IO_MakeOutDirs();
void IO_WriteCaseDim(); 		//Write case dimention in the file ./dim.dat
void IO_ArrayToFile(const char* fname, const double* array, int size);
void IO_DivToFile();

#endif