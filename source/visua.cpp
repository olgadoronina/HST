#include "all.h"
#include <string>




//======================================================================================================================
void VisuaInit() {
//======================================================================================================================
	Num_Var = NUMVAR_VISUA;
	Var_Names = GimmeMem<string>(Num_Var, "VarNames");
	Var_Names[0] = "V_x";
	Var_Names[1] = "V_y";
	Var_Names[2] = "V_z";
	printf("VisuaInit()\n");

}
//======================================================================================================================
//======================================================================================================================
void tVisuaRecord::WriteHeader(){
//======================================================================================================================

	string str_dim = "0 "+ToString(Dims[Coor_Z]-1)+" 0 "+ToString(Dims[Coor_Y]-1)+" 0 "+ToString(Dims[Coor_X]-1);

	StartOpenTag("VTKFile");
	WriteAttribute("type = \"StructuredGrid\"");
	WriteAttribute("version = \"0.1\"");
	WriteAttribute("byte_order = \"LittleEndian\"");
	EndStartedTag();

	StartOpenTag("StructuredGrid");
	WriteAttribute("WholeExtent = \"" + str_dim + "\"");
	EndStartedTag();

	StartOpenTag("Piece");
	WriteAttribute("Extent = \"" + str_dim + "\"");
	EndStartedTag();

	printf("tVisuaRecord::WriteHeader()\n");
}
//======================================================================================================================

//======================================================================================================================
void tVisuaRecord::WritePointData(){
//======================================================================================================================
	StartOpenTag("PointData");
	// for (int var=0;var<Num_Var; var++)
	WriteAttribute("Scalars=\""+Var_Names[Coor_X]+"\"");
	EndStartedTag();
	for (int var=0;var<Num_Var; var++){
		StartOpenTag("DataArray");
		WriteAttribute("type=\"Float32\"");
		WriteAttribute("Name=\""+Var_Names[var]+"\"");
		WriteAttribute("format = \" ascii\"");
		WriteAttribute("NumberOfComponents = \"1\"");
		EndStartedTag();
		string str = "";
		for (int in=0; in<Nn; in++) {
			str += ToString(UA[in][var]) + "\t";
		}
		WriteData(str+"\n");
		CloseElement(); 
	}
	
	CloseElement();
	printf("tVisuaRecord::WritePointData()\n");
}
//====================================================================================================================== 


//======================================================================================================================
void tVisuaRecord::WriteCoord(){
//======================================================================================================================

	OpenElement("Points");

	StartOpenTag("DataArray");
	WriteAttribute("type=\"Int32\"");
	WriteAttribute("Name=\"Points\"");
	WriteAttribute("format = \" ascii\"");
	WriteAttribute("NumberOfComponents = \"3\"");
	EndStartedTag();
	for (int in=0; in<Nn; in++) {
		string str = ToString(Coor[in][Coor_X]) + "\t" + 
			  ToString(Coor[in][Coor_Y]) + "\t" +
			  ToString(Coor[in][Coor_Z]) + "\n";
		WriteData(str);
	}

	CloseElement();
	printf("tVisuaRecord::WriteCoord()\n");
}
//====================================================================================================================== 

 
//======================================================================================================================
void tVisuaRecord::WriteStructuredGrid(const char *filename) {
//======================================================================================================================
      
    OpenFile(filename); 

    WriteHeader();
    WritePointData();
    WriteCoord();
    FinishFile();
    CloseFile();

    printf("tVisuaRecord::WriteStructuredGridSuccess!\n");
}
//======================================================================================================================


