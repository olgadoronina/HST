#include "io.h"

int main(int argc, char *argv[]){

	IO_MakeOutDirs();

	InitMPI(&argc,&argv);

	InitCaseDim();

	InitCaseParams();

// //----------------------------------------------------------------------------------------------------------------------
// //  Initialize arrays
// //----------------------------------------------------------------------------------------------------------------------
// 	complex ux_k(1:nxhp,1:ny,1:mz),ux1_k(1:nxhp,1:nz,1:my)
// 	complex uy_k(1:nxhp,1:ny,1:mz),uy1_k(1:nxhp,1:nz,1:my)
// 	complex uz_k(1:nxhp,1:ny,1:mz),uz1_k(1:nxhp,1:nz,1:my)
 //  	Ux_k = GimmeMem<double>(Nx, "Ux_k");
	// Uy_k = GimmeMem<double>(Ny, "Uz_k");
	// Uz_k = GimmeMem<double>(Nz, "Uz_k");

// Coordinates
	CoorInit();

// Zero fields
	UA.Alloc(Nn, NumCoords, "UA"); UA = 0.0;   					// real part of velocity field
	UA_Im.Alloc(Nn, NumCoords, "UA"); UA_Im = 0.0;   			// imaginary part of velocity field
	RHS_UA.Alloc(Nn, NumCoords, "RHS_UA"); RHS_UA = 0.0;   		// velocity equation rhs

// Set mean shear profile
	InitMeanShear();
// Setup wavenumbers
	WaveNumSetup();

// Initialize or load Fourier fields 




// Terminate MPI.
	MPI_Finalize();

// Array deallocation	
	Coor.Dealloc();
	UA.Dealloc();
	RHS_UA.Dealloc();
	FreeMem(ux_eq);

	if ( MyID == 0 ) 
        printf ( "Normal end of execution.\n" );
  	return 0;
}