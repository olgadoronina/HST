#ifndef UTILS_HEADER
#define UTILS_HEADER

#ifdef GLOBAL_DATA
   #define EXTERNAL
   #define DEFAULT(X) = X
#else
#define EXTERNAL extern
   #define DEFAULT(X)
#endif

#include "global.h"

EXTERNAL double* Div_Im DEFAULT(NULL);


int FillMask(int localNn);
double EnergyCalc(double k,double kmax);
void ResetUA();
void Symmetrize();
void CheckDivergence();
double MaxAbsValue(const double* array, int N);

#endif