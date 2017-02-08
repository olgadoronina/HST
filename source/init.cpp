#define GLOBAL_DATA
#include "all.h"


//======================================================================================================================
void InitCaseDim(){ // Define case dimention according possible symmetry
//======================================================================================================================
	IO_WriteCaseDim();				// Write case dimensions to file
	
	if (SYMMETRIC == 1)
		if (Nx_init/2. == 0) {
			Nx = Nx_init/2.+1;
		    if( MyID == 0 ) printf("InitCaseDim(): Due simmetry new Nx = %d\n",Nx);
		} 
		else {
			Nx = Nx_init/2;
			if( MyID == 0 ) printf("InitCaseDim(): Due simmetry new Nx = %d\n",Nx);
		}

	if ( IfInt(Nx/(double)num_procs) && IfInt(Nx/(double)num_procs) && IfInt(Nx/(double)num_procs) ) {
		NnLocalx = Nx/num_procs;
		NnLocaly = Ny/num_procs;
		NnLocalz = Nz/num_procs;
		if( MyID == 0 ) 
			printf("InitCaseDim():  NnLocalx = %d\tNnLocaly = %d\tNnLocalz = %d\n",NnLocalx,NnLocaly,NnLocalz);
	}
	else 
		crash("InitCaseDim(): NnLocal: Nx or Ny or Nz should be divisible by num_procs, DO IT MANUALLY\n");

	Nn = Nx*Ny*Nz;
	printf("InitCaseDim():  Nn = %d\n",Nn);
	
	if (IfInt(Nx/(double)num_procs)) {
		NnLocal = Nn/num_procs;
		printf("InitCaseDim():  NnLocal = %d\n",NnLocal);
	}
	else crash("InitCaseDim(): NnLocal: Nx*Ny*Nz should be divisible by num_procs, DO IT MANUALLY\n");
 
}
//======================================================================================================================

//======================================================================================================================
void InitCaseParams(){ // Set case parameters
//======================================================================================================================
	ipic = ipics;  					// Initial snapshot (=0 for new)
	DT = dts; 						// Initial time step
	TIME = times;  					// Initial time
	noinitialrandom = false;	
}
//======================================================================================================================

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
    //int ll; if(MPI_Get_processor_name(proc_name, &ll) != MPI_SUCCESS) { fprintf(stderr,"MPI_Get_processor_name failed \n"); exit(0);}
    //MPIinitialized = 1;
	// ierr = MPI_Init(&argc, &argv); 							// Initialize MPI.
	// ierr = MPI_Comm_rank(MPI_COMM_WORLD, &MyID); 			// Determine this process's rank.
	// ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  		// Determine the number of available processes.

	if ( MyID == 0 ) 
		printf("Hello! My rank (my_id) is %d of %d\n", MyID, num_procs);
}
//======================================================================================================================

//======================================================================================================================
void InitMeanShear(){ // Set mean shear profile
//======================================================================================================================

	if(ux_eq != NULL) crash("InitMeanShear(): array ux_eq already not NULL");

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
}

//======================================================================================================================
void CoorInit(){ // Fullfill the Coor array (array of coordinates)
//======================================================================================================================
	Coor.Alloc(Nn, 3, "Coor"); 									// Coordinates
	for (int ix=0; ix<Nx; ix++)	
		for (int iy=0; iy<Ny; iy++)
			for (int iz=0; iz<Nz; iz++) {
				Coor[ ix*Ny*Nz + iy*Nz + iz ][Coor_X] = ix;
				Coor[ ix*Ny*Nz + iy*Nz + iz ][Coor_Y] = iy;
				Coor[ ix*Ny*Nz + iy*Nz + iz ][Coor_Z] = iz;
			}
}
//======================================================================================================================


//======================================================================================================================
void WaveNumSetup(){ // Set case parameters
//======================================================================================================================
	Wave_num_x = GimmeMem<double>(Nx, "Wave_num_x");
	Wave_num_y = GimmeMem<double>(Ny, "Wave_num_y");
	Wave_num_z = GimmeMem<double>(Nz, "Wave_num_z");
	for( int ix=0; ix<Nx; ix++) Wave_num_x[ix]=(double)ix;							// Setup horizontal wavenumbers
	for( int iy=0; iy<Ny; iy++) Wave_num_y[iy]=(double)(iy-(iy+1)*Ny/(Ny/2+2));		// Setup vertical wavenumbers
	for( int iz=0; iz<Nz; iz++) Wave_num_z[iz]=(double)(iz-(iz+1)*Nz/(Nz/2+2));		// Setup spanwise wavenumbers
}
//======================================================================================================================


//======================================================================================================================
void TurbFieldInit(){ // //Create initial turbulent fields
//======================================================================================================================
	double kmax = sqrt(2)*Nx_init/3;

	for (int n=0; n<NnLocal; n++) {
		int in = MyID*NnLocal+n; 		//определяем номера для узла
	 	
	 	double mask = FillMask(in); 	// посчитали загадочный флаг
	 	//printf ( "Turb 1\n" );
		int X = Coor[in][Coor_X];		
		int Y = Coor[in][Coor_Y];
		int Z = Coor[in][Coor_Z];

		//if (Z == Nz/2) continue;  // если попали в середину? (проверить нужени ли +1) по z, то почкму-то ничего не делаем !!! Проверить

  		double Wave_num = sqrt(SQR(Wave_num_x[X]) + SQR(Wave_num_y[Y]) + SQR(Wave_num_z[Z]));
		
		//  Spectral function call	   
		double energy_k = EnergyCalc(Wave_num,kmax); 		// energy calculated using wave number
		double ef = mask*sqrt(energy_k)/sqrt(Pi2)/Wave_num;
		//printf ( "Turb 2 %d\n", in );
		if ( ef <= 0 || mask == 0. || in == 0) {   //also set the k=0 mode to (0.0.0)
			UA[in][Coor_X] = 0.0; 	UA_Im[in][Coor_X] = 0.0; 
			UA[in][Coor_Y] = 0.0; 	UA_Im[in][Coor_Y] = 0.0; 
			UA[in][Coor_Z] = 0.0; 	UA_Im[in][Coor_Z] = 0.0; 
			continue;
		}

		//----------Generate random phases
		double theta1 = (double) rand();
		double theta2 = (double) rand();
		double theta3 = (double) rand();   
		double phi    = (double) rand();  

		// original formula : argc=ima*Pi2*theta1; alpha =ef*cexp(argc)*cos(Pi2*phi) , where complex ima = (0,1). Use that e^(i* 2*pi*theta) = i*sin(2*pi*theta)
		double alpha_Im = ef*sin(theta1)*cos(phi);
		double beta_Im  = ef*sin(theta2)*sin(phi);
		double delta_Im = ef*sin(theta3);

		if ( Wave_num_x[X] ==0 && Wave_num_y[Y]) {
			UA_Im[in][Coor_X] = alpha_Im;    
			UA_Im[in][Coor_Y] = beta_Im;
			UA_Im[in][Coor_Z] = 0;   
		} 
		else {
			double coef1 = alpha_Im*Wave_num;
			double coef2 = beta_Im*Wave_num_z[Z];
			double coef3 = Wave_num*sqrt(SQR(Wave_num_x[X]) + SQR(Wave_num_y[Y]));

			UA_Im[in][Coor_X] = (coef1*Wave_num_y[Y] + coef2*Wave_num_x[X])/coef3 ;    // странная формула
			UA_Im[in][Coor_Y] = (coef2*Wave_num_y[Y] + coef1*Wave_num_x[X])/coef3;
			UA_Im[in][Coor_Z] = -beta_Im*sqrt(SQR(Wave_num_x[X]) + SQR(Wave_num_y[Y]))/Wave_num; 
	   	}
	   	//printf ( "Turb 3\n" );
	   	// Real part equal zero in any case
	   	// !!! Почему-то в оригинале заполняют только массив локального размера по z, то есть NnLocal_z (для комплексно часть тоже)
		UA[in][Coor_X] = 0.0; 
		UA[in][Coor_Y] = 0.0; 
		UA[in][Coor_Z] = 0.0; 
	}
	printf ( "Turb Finished\n" );

}

	

//======================================================================================================================




//======================================================================================================================
// //======================================================================================================================
// void ResetModes(){ // Restore first six modes
// //======================================================================================================================
//     if( MyID == 0 ){      
// 		ux_k(1,2,1)= CMPLX(-4.*S/pi/pi     ,0.)
// 	   	ux_k(1,4,1)= CMPLX(-4.*S/pi/pi/  9.,0.)
// 	   	ux_k(1,6,1)= CMPLX(-4.*S/pi/pi/ 25.,0.)
// 	   	ux_k(1,8,1)= CMPLX(-4.*S/pi/pi/ 49.,0.)
// 	   	ux_k(1,10,1)=CMPLX(-4.*S/pi/pi/ 81.,0.)
// 	   	ux_k(1,12,1)=CMPLX(-4.*S/pi/pi/121.,0.)
//     	printf("ResetModes(): Restore first six modes\n");
//     }
//     //Barrier();
// }