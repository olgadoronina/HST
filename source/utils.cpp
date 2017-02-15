#include "all.h"
//======================================================================================================================
double FillMask(int local_n){
//======================================================================================================================
	//printf ( "MyID = %d\tFillMask 0\t%d\n", MyID,local_n);

	int X = Coor[local_n][Coor_X];		
	int Y = Coor[local_n][Coor_Y];
	int Z = Coor[local_n][Coor_Z];
	if (Coor[local_n][Coor_Z] == Nz/2) return 0.0;          // Wave_num_z=N/2 always gets mask=0. 
	if (Coor[local_n][Coor_Y] == Ny/2) return 0.0; 			// Wave_num_y=N/2 always gets mask=0.
	//printf ( "MyID = %d\tFillMask 1\t%d\n", MyID,local_n);

	double  Wave_num_SQR = SQR(Wave_num_x[X]/Nx_init) + SQR(Wave_num_y[Y]/Ny) + SQR(Wave_num_z[Z]/Nz);
	//printf ( "MyID = %d\tFillMask 2\t%d\n", MyID,local_n);
	if(Wave_num_SQR > TwoNinth) 
		return 0.0;
	else
	 	return 1.0;
}
//======================================================================================================================

//======================================================================================================================
double EnergyCalc(double k,double kmax) { // specification of initial 3-d energy spectrum in waveno space
//======================================================================================================================
	double k0 = 1.;
	//-------128/256/512 cube	 
	double kp = 6.87;
	double gamma = 7.5e-5;
	//-------64 cube	
	// double kp = 8.08;
	// double gamma = 7.888e-4;

	if (k >= k0 && k <= kp) 
	 	return gamma*k*k; 
	else if (k > kp && k <= kmax)
		return pow(gamma*kp,11./3.)*pow(k,-5./3.);
	else
		return 0.0;
}
//======================================================================================================================

//======================================================================================================================
void ResetUA() { // Restore lowest modes to mean shear
//======================================================================================================================
//Restore first six modes // kz=0 
	if ( MyID == 0 ) {
		int in = GetNodeNum(0,1,0);
		UA[in][Coor_X] = -4*shear/PiNumber/PiNumber;
		in = GetNodeNum(0,3,0);
		UA[in][Coor_X] = -4*shear/PiNumber/PiNumber/9.; 
		in = GetNodeNum(0,5,0);
		UA[in][Coor_X] = -4*shear/PiNumber/PiNumber/25.; 
		in = GetNodeNum(0,7,0);
		UA[in][Coor_X] = -4*shear/PiNumber/PiNumber/49.; 
		in = GetNodeNum(0,9,0);
		UA[in][Coor_X] = -4*shear/PiNumber/PiNumber/81.; 
		in = GetNodeNum(0,11,0);
		UA[in][Coor_X] = -4*shear/PiNumber/PiNumber/121.;
	}
}
//======================================================================================================================

//======================================================================================================================
void Symmetrize() { // Apply symmetry constraint
//======================================================================================================================
// Symmetrize fields
	int Nn_Ymiddle = Nn/2; // just because Coor array was fulfilled this way
	for (int in=0; in<Nn_Ymiddle; in++) {
		// Zero-out bottom plane       
	    if ( Coor[in][Coor_Y] == 0 ) {
	     	UA[in][Coor_Y] = 0.0; 
	     	UA_Im[in][Coor_Y] = 0.0;
	    }
	    else {
		int symm_nodeNum = GetNodeNum(Coor[in][Coor_X], Ny - Coor[in][Coor_Y], Coor[in][Coor_Z]);
		UA[symm_nodeNum][Coor_X] = UA[in][Coor_X];	UA_Im[symm_nodeNum][Coor_X] = UA_Im[in][Coor_X];
		UA[symm_nodeNum][Coor_Y] = UA[in][Coor_Y];	UA_Im[symm_nodeNum][Coor_Y] = UA_Im[in][Coor_Y]; 
		UA[symm_nodeNum][Coor_Z] = UA[in][Coor_Z];	UA_Im[symm_nodeNum][Coor_Z] = UA_Im[in][Coor_Z];
		}
	}

}
//======================================================================================================================
//======================================================================================================================
void CheckDivergence() { // Check divergence
//======================================================================================================================
	if ( Div_Im = NULL ) Div_Im = GimmeMem<double>(Nn, "Div_Im"); 
	// Calculate divergence
	for (int in; in<Nn; in++)
    	Div_Im[in]= Wave_num_x[Coor[in][Coor_X]] * UA[in][Coor_X] + 
    			    Wave_num_y[Coor[in][Coor_Y]] * UA[in][Coor_Y] +
    			    Wave_num_z[Coor[in][Coor_Z]] * UA[in][Coor_Z];

	IO_DivToFile();
}
//======================================================================================================================
	

//======================================================================================================================
double MaxAbsValue(const double* array, int N) { // Check divergence
//======================================================================================================================
	double max = 0.0;
	for (int i; i<N; i++) 
		if ( fabs(array[i]) > max ) max = fabs(array[i]);
	return max;
}
//======================================================================================================================

