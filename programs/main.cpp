#include "all.h"

int main(int argc, char *argv[]){

	InitMPI(&argc,&argv);
	//printf ( "MyID = %d\tStep 1\n", MyID );
	IO_MakeOutDirs();
	//printf ( "MyID = %d\tStep 2\n", MyID );
	InitCaseDim();
	//printf ( "MyID = %d\tStep 3\n", MyID );
	InitCaseParams();
	//printf ( "MyID = %d\tStep 4\n", MyID );
	VisuaInit();
	// printf ( "MyID = %d\tStep 5\n", MyID );
	CoorInit();
	// printf ( "MyID = %d\tStep 6\n", MyID );
// Zero fields
	UA.Alloc(Nn, NumCoords, "UA"); UA = 2.0;   					// real part of velocity field
	UA_Im.Alloc(Nn, NumCoords, "UA"); UA_Im = 2.0;   			// imaginary part of velocity field
	//RHS_UA.Alloc(Nn, NumCoords, "RHS_UA"); RHS_UA = 0.0;   		// velocity equation rhs
	// printf ( "MyID = %d\tStep 7\n", MyID );
// Set mean shear profile
	InitMeanShear(); // matches!
	// printf ( "MyID = %d\tStep 8\n", MyID );
// Setup wavenumbers
	WaveNumSetup(); // matches!
	// printf ( "MyID = %d\tStep 9\n", MyID );
// Initialize Fourier fields 
	TurbFieldInit();
	// printf ( "MyID = %d\tStep 10\n", MyID );
	if ( MyID ==0 ) UA.PrintToFile("./OUTPUT/initUA.dat");
	if ( MyID ==0 ) UA_Im.PrintToFile("./OUTPUT/initUA_Im.dat");
	if ( MyID ==0 ) {
		tVisuaRecord record;
		record.WriteStructuredGrid("./OUTPUT/VISUA/init.vts");
	}



	ResetUA();
	// printf ( "MyID = %d\tStep 8\n", MyID );
	Symmetrize();
	// printf ( "MyID = %d\tStep 9\n", MyID );

	//UA.PrintToFile("./OUTPUT/initUA.dat");
	// if ( MyID ==0 ) {
	// 	tVisuaRecord record;
	// 	record.WriteStructuredGrid("./OUTPUT/VISUA/init.vts");
	// }




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