#ifndef IO_HEADER
#define IO_HEADER

#ifdef GLOBAL_DATA
   #define EXTERNAL
   #define DEFAULT(X) = X
#else
#define EXTERNAL extern
   #define DEFAULT(X)
#endif

// Default system includes
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
using namespace std;

#include <mpi.h>

// for system (mkdir)
#include <sys/stat.h>
#include <sys/types.h>


// Constants
#define PiNumber  3.141592653589793238462643383279
#define Pi2       6.283185307179586476925286766558
#define kTiny     1e-16
#define kTinyFlt  1e-8
#define kHuge     1e+16
#define kHugeFlt  1e+8

#define TwoNinth 2/9

//----------------------------------------------------------------------------------------------------------------------
// Position of coordinates and variables
//----------------------------------------------------------------------------------------------------------------------
#define Coor_X 0                  //позиция X соординаты
#define Coor_Y 1                  //позиция Y координаты
#define Coor_Z 2                  //позиция Z координаты

//#define Var_R  0                  //позиция плотности
#define Var_U  0                  //позиция Х-скорости
#define Var_V  1                  //позиция Y-скорости
#define Var_W  2                  //позиция Z-скорости
// #define Mask   3                   //позиция флага Mask

//#define Var_P  4                  //позиция давления
//#define Var_N  5                  //число базовых переменных без учета турбуля
// #define Var_NN 25                 //квадрат предыдушего числа

// #define Var_Nu 5                  //позиция Ню       
// #define Var_K  5                  //позиция К
// #define Var_D  6                  //позиция Эпсилон

// #define Var_Reserv 7              //стратегический резерв
// #define Var_NumMax 8              //пусть уж кратно 64 байтам
// //----------------------------------------------------------------------------------------------------------------------



    
EXTERNAL int NnLocalx, NnLocaly, NnLocalz, Nx;                 // Number of grid nodes for each processor
EXTERNAL int NnLocal;
EXTERNAL int Nn;                 // Number of grid nodes
//----------------------------------------------------------------------------------------------------------------------

EXTERNAL int KT DEFAULT(0);                     // Step number
EXTERNAL double DT DEFAULT(0);                  // Time step						
EXTERNAL double TIME DEFAULT(0);					
EXTERNAL bool noinitialrandom;

#define NumCoords 3
#define crash(...) exit(fprintf(stderr, __VA_ARGS__))
//----------------------------------------------------------------------------------------------------------------------



void MakeDir(string path);
void IO_MakeOutDirs();
void IO_WriteCaseDim(); 		//Write case dimention in the file ./dim.dat
void IO_ArrayToFile(const char* fname, const double* array, int size);



EXTERNAL MPI_Status status;
EXTERNAL int MyID, ierr, num_procs;
EXTERNAL int ipic;

#define MPI_SUCCESS 0


void InitMeanShear();
void InitMPI(int *argc,char ***argv);
void InitCaseDim();
void InitCaseParams();
void CoorInit();
void WaveNumSetup();


double FillMask(int localNn);
double EnergyCalc(double k,double kmax);
// Constatnts
//===Grid dimensions=====================================================
//int nyg,nzg,my,mz;
//int nxh,nxhp,nyh,nyhp,nzh,nzhp;
//nyg = nyzg;
//nzg=nyzg;
//my=ny/nyg;
//mz=nz/nzg;
//nxh=nx/2;
//nxhp=nxh+1;
//nyh=ny/2;
//nyhp=nyh+1;
//nzh=nz/2;
//nzhp=nzh+1;



//Should be in params
#define Ny 64                 // Grid dimension y (ny must equal nz)
#define Nz 64                 // Grid dimension z (ny must equal nz)
#define Nx_init 64                 // Grid dimension x
//#define nyzg 2                 // Number of processors (must not be greater than ny/2)
#define R_inv 1./300.          // Viscosity (~inverse Reynolds number) 
#define shear 0.5                  // Shear rate for the flow = 2*S/pi
#define epsin 0.01             // Energy dissipation rate
#define Ly Pi2                 // Physical dimension y
#define ipics 0                // Initial snapshot (=0 for new)
#define dts 0.0001             // Initial time step
#define times 0.               // Initial time 
#define kshift 3               // Type of shift (see shifts.f)

#define Num_step 100              // Number of time steps
#define iout 10                // Time steps between outputs


#define SYMMETRIC 1
\

//======================================================================================================================
inline int IfInt(double number) { return fmod(number,1) == 0; } // Check if number is integer
//======================================================================================================================
inline int SQR(int x){ return x*x;}
//======================================================================================================================
// Memory interface
//======================================================================================================================
template <typename ValueType>
inline ValueType* GimmeMem(size_t N, const char* label = NULL){
    if(N==0) return NULL;
    ValueType* p=NULL;
    try{p = new ValueType[N];}
    catch (bad_alloc &ba){
        if(p!=NULL) delete[] p;
        //crash("GimmeMem: memory fuckup! N = %llu, sizeof = %d, name %s, text: %s\n", 
              //(unsigned long long) N,sizeof(ValueType), label ? label:"noname",ba.what());
    }
    return p;
}

template <typename ValueType> 
inline void FreeMem(ValueType* &Array) { 
    delete[] Array; 
    Array=NULL; 
}
//======================================================================================================================

//======================================================================================================================
// Default data type for block 2D arrays NxM on a base of flat array
//======================================================================================================================
template <typename ValueType> class tBlockArray {
public:
    typedef tBlockArray<ValueType> BaseType;
    int N; // number of blocks
    int M; // size of block
    ValueType * V; // data array

    tBlockArray(){M=0;N=0;V=NULL;}
    tBlockArray(int n, int m, const char* label = NULL){
        if((n<=0)||(m<=0))crash("tBlockArray constructor: wrong input %d %d", n,m);
        V = GimmeMem<ValueType>(n*m, label); M=m; N=n;
    }
    ~tBlockArray(){ 
        if(M>0 && V) FreeMem(V); 
    }
    
    inline bool Allocated()const{return (M>0); }
    
    inline int Alloc(int n, int m, const char* label = NULL){
        if((n<=0)||(m<=0))crash("tBlockArray Alloc: wrong input %d %d", n,m);
        if(M>0)crash("tBlockArray Alloc: already allocated!\n");
        unsigned long long Size = (unsigned long long)n * (unsigned long long)m;
        if(Size > 0x7FFFFFFF) crash("tBlockArray Alloc: size of array %llu exceeds maximal value %i", (unsigned long long) Size, 0x7FFFFFFF);
        V = GimmeMem<ValueType>(size_t(Size), label); M=m; N=n;
        return 0;
    }
    inline void Dealloc(){ if((M>0)&& V) FreeMem(V); N=0; M=0; V=NULL; }

    inline BaseType& operator=(const BaseType& object){
        if(this != &object){
            if((M==0)&&(object.M>0)) Alloc(object.N,object.M,"tBlockArray="); //if not yet allocated
            if(object.M!=M) crash("tBlockArray: assignment of arrays with different size of block %d %d", object.M,M);
            if(object.N!=N) crash("tBlockArray: assignment of arrays with different size", object.N, N);
            for(int i=0; i<N*M; i++) V[i]=object.V[i];
        }
        return *this;
    }
    inline BaseType& operator+=(BaseType& object){
        if(this != &object){
            if(object.M!=M) crash("tBlockArray: assignment of arrays with different size of block %d %d", object.M,M);
            if(object.N!=N) crash("tBlockArray: assignment of arrays with different size", object.N, N);
            for(int i=0; i<N*M; i++) V[i]+=object.V[i];
        }
        return *this;
    }
    inline BaseType& operator=(ValueType x){
        if(M==0) crash("tBlockArray =: not allocated!");
        for(int i=0; i<N*M; i++) V[i]=x;
        return *this;
    }
    inline BaseType& operator*=(ValueType x){
        if(M==0) crash("tBlockArray *=: not allocated!");
        for(int i=0; i<N*M; i++) V[i]*=x;
        return *this;
    }
    inline ValueType *operator[](int i) const {return V+i*M;}

    // prototypes for additional debug functionality for tBlockArray (in utils.cpp)
    //void PrintToFile(const char* fname);
 private:
    tBlockArray(const BaseType &b); // force compiler not to create default copy constructor    
};
//======================================================================================================================

EXTERNAL tBlockArray<double> UA;       // velocity field
EXTERNAL tBlockArray<double> UA_Im;     // Imaginary part of solution
EXTERNAL tBlockArray<double> RHS_UA;   // velocity equation rhs
EXTERNAL tBlockArray<int> Coor;     // Coordinates


EXTERNAL double* ux_eq DEFAULT(NULL); 

// EXTERNAL double* Ux_k DEFAULT(NULL);
// EXTERNAL double* Uy_k DEFAULT(NULL);
// EXTERNAL double* Uz_k DEFAULT(NULL);


EXTERNAL double* Wave_num_x DEFAULT(NULL);
EXTERNAL double* Wave_num_y DEFAULT(NULL);
EXTERNAL double* Wave_num_z DEFAULT(NULL);


//===Wavenumbers=========================================================


#endif