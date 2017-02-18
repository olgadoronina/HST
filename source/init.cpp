#define GLOBAL_DATA
#include "all.h"


//======================================================================================================================
void InitCaseDim(){ // Define case dimension according possible symmetry
//======================================================================================================================
	if (SYMMETRIC == 1) {
		if (Nx_init%2 != 0) {
			Nx = Nx_init/2+1;
		    if( MyID == 0 ) printf("InitCaseDim():\tDue symmetry new Nx = %d\n",Nx);
		} 
		else {
			Nx = Nx_init/2;
			if( MyID == 0 ) printf("InitCaseDim():\tDue symmetry new Nx = %d\n",Nx);
		}
	}
	Dims = GimmeMem<int>(NumCoords, "Dims");
	Dims[Coor_X] = Nx;
	Dims[Coor_Y] = Ny;
	Dims[Coor_Z] = Nz;

	if ( IfInt(Nx/(double)num_procs) && IfInt(Nx/(double)num_procs) && IfInt(Nx/(double)num_procs) ) {
		Dims_Local = GimmeMem<int>(NumCoords, "Dims_Local");
		Dims_Local[Coor_X] = Nx/num_procs;
		Dims_Local[Coor_Y] = Ny/num_procs;
		Dims_Local[Coor_Z] = Nz/num_procs;
		if( MyID == 0 ) 
			printf("InitCaseDim():\tNnLocalx = %d\tNnLocaly = %d\tNnLocalz = %d\n", 
						Dims[Coor_X], Dims[Coor_Y], Dims[Coor_Z]);
	}
	else 
		crash("InitCaseDim():\tNnLocal: Nx or Ny or Nz should be divisible by num_procs, DO IT MANUALLY\n");

	Nn = Nx*Ny*Nz;
	printf("InitCaseDim():\tNn = %d\n",Nn);
	
	if (IfInt(Nx/(double)num_procs)) {
		NnLocal = Nn/num_procs;
		printf("InitCaseDim():\tNnLocal = %d\n",NnLocal);
	}
	else crash("InitCaseDim():\tNnLocal: Nx*Ny*Nz should be divisible by num_procs, DO IT MANUALLY\n");

	IO_WriteCaseDim();				// Write case dimensions to file
 
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
void CoorInit(){ // Fulfill the Coor array (array of coordinates)
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
	for( int iy=0; iy<Ny; iy++) Wave_num_y[iy]=(double)(iy-(iy+1)/Ny/(Ny/2+2));		// Setup vertical wavenumbers
	for( int iz=0; iz<Nz; iz++) Wave_num_z[iz]=(double)(iz-(iz+1)/Nz/(Nz/2+2));		// Setup spanwise wavenumbers
	// IO_ArrayToFile("./OUTPUT/Wave_num_x_cpp.txt", Wave_num_x, Nx);
	// IO_ArrayToFile("./OUTPUT/Wave_num_y_cpp.txt", Wave_num_y, Ny);
	// IO_ArrayToFile("./OUTPUT/Wave_num_z_cpp.txt", Wave_num_z, Nz);

}
//======================================================================================================================


//======================================================================================================================
void TurbFieldInit(){ // //Create initial turbulent fields
//======================================================================================================================
	double kmax = sqrt(2)*Nx_init/3; // Почему определено здесь а не в EnergyCalc()?

	for (int n=0; n<NnLocal; n++) {
		int in = MyID*NnLocal+n; 		//определяем номер для узла

	 	int mask = FillMask(in); 	// посчитали загадочный флаг

		int X = Coor[in][Coor_X];		
		int Y = Coor[in][Coor_Y];
		int Z = Coor[in][Coor_Z];
		printf("mask = %d\t%d\t%d\t%d\n",mask, X, Y, Z);
		//if (Z == Nz/2) continue;  // кейс включен в if ef <= 0 
		if (X >= Nx) continue;  // убрать нафиг
		
  		double Wave_num = sqrt(SQR(Wave_num_x[X]) + SQR(Wave_num_y[Y]) + SQR(Wave_num_z[Z]));
		
		//  Spectral function call	   
		double energy_k = EnergyCalc(Wave_num,kmax); 		// energy calculated using wave number
		double ef = mask*sqrt(energy_k)/sqrt(Pi2)/Wave_num;

		if ( ef <= 0 || in == 0) {   //also set the k=0 mode to (0.0.0)
			UA[in][Coor_X] = 0.0; 	UA_Im[in][Coor_X] = 0.0; 
			UA[in][Coor_Y] = 0.0; 	UA_Im[in][Coor_Y] = 0.0; 
			UA[in][Coor_Z] = 0.0; 	UA_Im[in][Coor_Z] = 0.0; 
			continue;
		}
		//----------Generate random phases
		double theta1 = 2*PiNumber*0.1;//(double) rand();
		double theta2 = 2*PiNumber*0.2;//(double) rand();
		// double theta3 = (double) rand();   // нигде не используется?
		double phi    = 2*PiNumber*0.3;//(double) rand();  

		// original formula : argc=ima*Pi2*theta1; alpha =ef*cexp(argc)*cos(Pi2*phi) , where complex ima = (0,1). 
		// Use that e^(i* 2*pi*theta) = i*sin(2*pi*theta)
		double alpha = ef*cos(theta1)*cos(phi);
		double beta  = ef*cos(theta2)*sin(phi);
		//double delta = ef*cos(theta3);  // нигде не используется?

		double alpha_Im = ef*sin(theta1)*cos(phi);
		double beta_Im  = ef*sin(theta2)*sin(phi);
		//double delta_Im = ef*sin(theta3);  // нигде не используется?
		printf("%.10e\t%.10e\t%d\t%d\t%d\n", alpha, alpha_Im, X, Y, Z);
		printf("%.10e\t%.10e\t%d\t%d\t%d\n", beta, beta_Im, X, Y, Z);
		
		if ( Wave_num_x[X]==0 && Wave_num_y[Y]==0) {
			UA[in][Coor_X] = alpha;		UA_Im[in][Coor_X] = alpha_Im;    
			UA[in][Coor_Y] = beta;		UA_Im[in][Coor_Y] = beta_Im;
			UA[in][Coor_Z] = 0.0;		UA_Im[in][Coor_Z] = 0.0;   
		} 
		else {
			printf("else\n");
			double coef1 = alpha*Wave_num;
			double coef2 = beta*Wave_num_z[Z];
			double denominator = Wave_num*sqrt(SQR(Wave_num_x[X]) + SQR(Wave_num_y[Y]));
			
			UA[in][Coor_X] = (coef1*Wave_num_y[Y] + coef2*Wave_num_x[X])/denominator ;    // странная формула
			UA[in][Coor_Y] = (coef2*Wave_num_y[Y] + coef1*Wave_num_x[X])/denominator;
			UA[in][Coor_Z] = -beta*sqrt(SQR(Wave_num_x[X]) + SQR(Wave_num_y[Y]))/Wave_num; 
			
			coef1 = alpha_Im*Wave_num;
			coef2 = beta_Im*Wave_num_z[Z];
			printf("%.10e\t%.10e\t%d\t%d\t%d\n", coef1, coef2, X, Y, Z);
			UA_Im[in][Coor_X] = (coef1*Wave_num_y[Y] + coef2*Wave_num_x[X])/denominator ;    // странная формула
			UA_Im[in][Coor_Y] = (coef2*Wave_num_y[Y] + coef1*Wave_num_x[X])/denominator;
			UA_Im[in][Coor_Z] = -beta_Im*sqrt(SQR(Wave_num_x[X]) + SQR(Wave_num_y[Y]))/Wave_num; 
	   	}
	   	//printf ( "Turb 3\n" );
	   	// Real part equal zero in any case
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