#include "io.h"


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
		if (fp == NULL) crash("IO_ArrayToFile: Can't open file (probably there is no OUTPUT folder)\n");
    	fprintf(fp, "%d\t%d\t%d\t%d\n", Nx, Ny, Nz, num_procs);
    	//fputs("This is testing for fputs...\n", fp);
    	fclose(fp);
    	printf("IO_WriteCaseDim(): Write case dimentions\n");
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
    	printf("IO_ArrayToFile(): Write array in file\n");
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
        printf("Maximum divergence (time, max div):\n%f\t%f\n", TIME, div_max);
    }
    //Barrier();
}


   
