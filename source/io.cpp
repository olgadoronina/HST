#include "all.h"


//======================================================================================================================
void MakeDir(string path) { mkdir(path.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); }
//======================================================================================================================
//======================================================================================================================
void IO_MakeOutDirs(){  //creates output dirs - history, visua, recovery
//======================================================================================================================
    if( MyID == 0 ){
        MakeDir("./OUTPUT/");
        MakeDir("./OUTPUT/HISTORY/");
        MakeDir("./OUTPUT/VISUA/");
        //MakeDir("OUTPUT/RECOVERY/");
        printf("MyID = %d\tIO_MakeOutDirs(): Output directories created\n", MyID);
    }
    //Barrier();
}
//======================================================================================================================

//======================================================================================================================
void IO_WriteCaseDim(){ // Write case dimensions to file
//======================================================================================================================
    if( MyID == 0 ){
		FILE *fp = fopen("./OUTPUT/dims.dat", "wt");
		if (fp == NULL) crash("IO_WriteCaseDim(): Can't open file (probably there is no OUTPUT folder)\n");
    	fprintf(fp, "%d\t%d\t%d\t%d\n", Dims[Coor_X] , Dims[Coor_Y] , Dims[Coor_Z] , num_procs);
    	//fputs("This is testing for fputs...\n", fp);
    	fclose(fp);
    	printf("MyID = %d\tIO_WriteCaseDim(): Write case dimensions\n", MyID);
    }
    //Barrier();
}

//======================================================================================================================
void IO_ArrayToFile(const char* fname, const double* array, int size){ // Write array in file
//======================================================================================================================
    if( MyID == 0 ){      
		FILE *fp = fopen(fname, "wt");
		if (fp == NULL) crash("IO_ArrayToFile: Can't open file");
    	for(int i=0; i<size; i++) {
            fprintf(fp, "%25.15e     ", array[i]);
        	fprintf(fp, "\n");
    	}
    	fclose(fp);
    	printf("MyID = %d\tIO_ArrayToFile(): Write array in file\n", MyID);
    }
    //Barrier();
}

//======================================================================================================================
void IO_DivToFile() { // Write divergence to file  // надо распаллалелить !!!!
//======================================================================================================================
    if( MyID == 0 ){      
        double div_max = MaxAbsValue(Div_Im, Nn);
        FILE *fp = fopen("./OUTPUT/div.dat", "at");
        if (fp == NULL) crash("IO_DivToFile(): Can't open file (probably there is no OUTPUT folder)\n");
        fprintf(fp, "Maximum divergence (time, max div):\n%f\t%f\n", TIME, div_max);
        fclose(fp);
        printf("MyID = %d\tMaximum divergence (time, max div):\n%f\t%f\n", MyID, TIME, div_max);
    }
    //Barrier();
}

  
//======================================================================================================================
//--- B l o c k A r r a y   m e t h o d s -------------------------------------
//======================================================================================================================

//-----------------------------------------------------------------------------
template<> void tBlockArray<int>::PrintToFile(const char* fname) {
    FILE *in = fopen(fname,"wt");
    int i, j;
    for(i=0; i<this->N; i++) {
        for(j=0; j<this->M; j++)
            fprintf(in, "%9i ", (*this)[i][j]);
        fprintf(in, "\n");
    }
    fclose(in);
}

//-----------------------------------------------------------------------------
template<> void tBlockArray<double>::PrintToFile(const char* fname) {
    FILE *in = fopen(fname,"wt");
    int i, j;
    for(i=0; i<this->N; i++) {
        for(j=0; j<this->M; j++)
            fprintf(in, "%25.15e     ", (*this)[i][j]);
        fprintf(in, "\n");
    }
    fclose(in);
}