#include "io.h"

//======================================================================================================================
void MakeDir(string path) { mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
//======================================================================================================================
//======================================================================================================================
void IO_MakeOutDirs(){  //creates output dirs - history, visua, recovery
//======================================================================================================================
    if( MyID == 0 ){
        MakeDir("OUTPUT/");
        MakeDir("OUTPUT/HISTORY/");
        MakeDir("OUTPUT/VISUA/");
        //MakeDir("OUTPUT/RECOVERY/");
        printf("IO_MakeOutDirs(): Output directories created\n");
    }
    //Barrier();
}
//======================================================================================================================

//======================================================================================================================
void IO_WriteCaseDim(){ // Write case dimensions to file
//======================================================================================================================
    if( MyID == 0 ){
		FILE *fp = fopen("./OUTPUT/dims.dat", "wt");
    	fprintf(fp, "%d\t%d\t%d\t%d\n", Nx, Ny, Nz, num_procs);
    	//fputs("This is testing for fputs...\n", fp);
    	fclose(fp);
    	printf("IO_WriteCaseDim(): Write case dimentions\n");
    }
    //Barrier();
}

//======================================================================================================================
void InitCaseDim(){ // Define case dimention according possible symmetry
//======================================================================================================================
	IO_WriteCaseDim();				// Write case dimensions to file

	if (SYMMETRIC == 1)
		if (Nx/2. == 0) {
			Nx = Nx/2.+1;
		    if( MyID == 0 ) printf("InitCaseDim(): Due simmetry new Nx = %d\n",Nx);
		} 
		else {
			Nx = Nx/2;
			if( MyID == 0 ) printf("InitCaseDim(): Due simmetry new Nx = %d\n",Nx);
		}

	if ( IfInt(Nx/(double)num_procs) && IfInt(Nx/(double)num_procs) && IfInt(Nx/(double)num_procs) ) {
		NnLocalx = Nx/num_procs;
		NnLocaly = Ny/num_procs;
		NnLocalz = Nz/num_procs;
		if( MyID == 0 ) printf("InitCaseDim():  NnLocalx = %d\tNnLocaly = %d\tNnLocalz = %d\n",NnLocalx,NnLocaly,NnLocalz);
	}
	else 
		crash("NnLocal: Nx or Ny or Nz should be divisible by num_procs\n");
	
	Nn = Nx*Ny*Nz; 
}

//======================================================================================================================
void IO_ArrayToFile(const char* fname, const double* array, int size){ // Write array in file
//======================================================================================================================
    if( MyID == 0 ){      
		FILE *fp = fopen(fname, "wt");
    	for(int i=0; i<size; i++) {
            fprintf(fp, "%25.15e     ", array[i]);
        	fprintf(fp, "\n");
    	}
    	fclose(fp);
    	printf("IO_ArrayToFile(): Write array in file\n");
    }
    //Barrier();
}

//======================================================================================================================

main(int argc, char **argv){

	//y_grid_half = Ny/2;
	//y_grid_half_plusone = nyh+1;

//----------------------------------------------------------------------------------------------------------------------
//  Initialization of MPI
//----------------------------------------------------------------------------------------------------------------------
	MPI_Status status;
	
	ierr = MPI_Init(&argc, &argv); 							// Initialize MPI.
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &MyID); 			// Determine this process's rank.
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  		// Determine the number of available processes.

	if ( MyID == 0 ) 
		printf("Hello! My rank (my_id) is %d of %d\n", MyID, num_procs);

	InitCaseDim();
//----------------------------------------------------------------------------------------------------------------------
//  Set case parameters
//----------------------------------------------------------------------------------------------------------------------
	ipic = ipics;  					// Initial snapshot (=0 for new)
	DT = dts; 						// Initial time step
	TIME = times;  					// Initial time
	noinitialrandom = false;

//----------------------------------------------------------------------------------------------------------------------
//  Set mean shear profile
//----------------------------------------------------------------------------------------------------------------------
	static double* ux_eq = NULL; 
	ux_eq = GimmeMem<double>(Ny, "ux_eq");
	//Set mean shear in lower half of domain
	for (int i=0; i<Ny/2+1; i++) {
  		ux_eq[i]=-8.*shear/PiNumber/PiNumber*
		  		(cos(1. *Pi2*i/Ny)     + 
			     cos(3. *Pi2*i/Ny)/9.  + 
			     cos(5. *Pi2*i/Ny)/25. + 
			     cos(7. *Pi2*i/Ny)/49. + 
			     cos(9. *Pi2*i/Ny)/81. + 
			     cos(11.*Pi2*i/Ny)/121.); 
  	}
  	//Set mean shear in upper half of domain
	for (int i=Ny/2+1;i<Ny; i++)
   		ux_eq[i] = ux_eq[Ny+1-i];                  

  	IO_ArrayToFile("./OUTPUT/ux_eq_cpp.txt", ux_eq, Ny);

// //----------------------------------------------------------------------------------------------------------------------
// //  Initialize arrays
// //----------------------------------------------------------------------------------------------------------------------
// 	complex ux_k(1:nxhp,1:ny,1:mz),ux1_k(1:nxhp,1:nz,1:my)
// 	complex uy_k(1:nxhp,1:ny,1:mz),uy1_k(1:nxhp,1:nz,1:my)
// 	complex uz_k(1:nxhp,1:ny,1:mz),uz1_k(1:nxhp,1:nz,1:my)


// Coordinates
	Coor.Alloc(Nn, 3, "Coor"); 									// Coordinates
	for (int ix=0; ix<Nx; ix++)	
		for (int iy=0; iy<Ny; iy++)
			for (int iz=0; iz<Nz; iz++) {
				Coor[ix+iy+iz][Coor_X] = (double)ix;
				Coor[ix+iy+iz][Coor_Y] = (double)iy;
				Coor[ix+iy+iz][Coor_Z] = (double)iz;
			}

// Zero fields
	UA.Alloc(Nn, NumCoords, "UA"); UA = 0.0;   					// velocity field
	RHS_UA.Alloc(Nn, NumCoords, "RHS_UA"); RHS_UA = 0.0;   		// velocity equation rhs

// Setup wavenumbers
	double* Wave_num_x = GimmeMem<double>(Nx, "Wave_num_x");
	double* Wave_num_y = GimmeMem<double>(Ny, "Wave_num_y");
	double* Wave_num_z = GimmeMem<double>(Nz, "Wave_num_z");
	for( int i=0; i<Nx; i++) Wave_num_x[i]=(double)i;							// Setup horizontal wavenumbers
	for( int i=0; i<Ny; i++) Wave_num_y[i]=(double)(i-(i+1)*Ny/(Ny/2+2));		//Setup vertical wavenumbers
	for( int i=0; i<Ny; i++) Wave_num_y[i]=(double)(i-(i+1)*Nz/(Nz/2+2));			//Setup spanwise wavenumbers

	// Terminate MPI.
	ierr = MPI_Finalize();
	Coor.Dealloc();
	UA.Dealloc();
	RHS_UA.Dealloc();
	FreeMem(ux_eq);

	if ( MyID == 0 ) 
        printf ( "Normal end of execution.\n" );
  	return 0;
}