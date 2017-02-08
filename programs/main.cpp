#include "all.h"

int main(int argc, char *argv[]){

	IO_MakeOutDirs();
	printf ( "Step 1\n" );
	InitMPI(&argc,&argv);
	printf ( "Step 2\n" );
	InitCaseDim();
	printf ( "Step 3\n" );
	InitCaseParams();
	printf ( "Step 4\n" );
// //----------------------------------------------------------------------------------------------------------------------
// //  Initialize arrays
// //----------------------------------------------------------------------------------------------------------------------
// Coordinates
	CoorInit();
	printf ( "Step 5\n" );
// Zero fields
	UA.Alloc(Nn, NumCoords, "UA"); UA = 0.0;   					// real part of velocity field
	UA_Im.Alloc(Nn, NumCoords, "UA"); UA_Im = 0.0;   			// imaginary part of velocity field
	//RHS_UA.Alloc(Nn, NumCoords, "RHS_UA"); RHS_UA = 0.0;   		// velocity equation rhs
	printf ( "Step 6\n" );
// Set mean shear profile
	InitMeanShear();
	printf ( "Step 7\n" );
// Setup wavenumbers
	WaveNumSetup();
	printf ( "Step 8\n" );
// Initialize Fourier fields 
	TurbFieldInit();
	printf ( "Step 9\n" );
	ResetUA();
	printf ( "Step 10\n" );
	Symmetrize();
	printf ( "Step 11\n" );


// Terminate MPI.
	MPI_Finalize();

// Array deallocation	
	Coor.Dealloc();
	UA.Dealloc();
	//RHS_UA.Dealloc();
	FreeMem(ux_eq);

	if ( MyID == 0 ) 
        printf ( "Normal end of execution.\n" );
  	return 0;
}