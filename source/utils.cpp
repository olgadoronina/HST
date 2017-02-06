/*
* @Author: Olga Doronina
* @Date:   2017-02-03 14:17:03
* @Last Modified by:   Olga Doronina
* @Last Modified time: 2017-02-06 14:15:27
*/

#include "io.h"
//======================================================================================================================
double FillMask(int local_n){
//======================================================================================================================
	if (Coor[local_n][Coor_Z] == Nz/2) return 0.0;          // Wave_num_z=N/2 always gets mask=0. 
	if (Coor[local_n][Coor_Y] == Ny/2) return 0.0; 			// Wave_num_y=N/2 always gets mask=0.
	
	double  Wave_num_SQR = SQR(Wave_num_x[local_n]/Nx_init) + SQR(Wave_num_y[local_n]/Ny) + SQR(Wave_num_z[local_n]/Nz);

	if(Wave_num_SQR > TwoNinth) 
		return 0.0;
	else
	 	return 1.0;
}
//======================================================================================================================

//======================================================================================================================
double EnergyCalc(double k,double kmax) {//specification of initial 3-d energy spectrum in waveno space
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