#ifndef VISUA_HEADER
#define VISUA_HEADER


#ifdef GLOBAL_DATA
   #define EXTERNAL
   #define DEFAULT(X) = X
#else
#define EXTERNAL extern
   #define DEFAULT(X)
#endif


#define NUMVAR_VISUA 3

EXTERNAL int Num_Var DEFAULT(0);
EXTERNAL string* Var_Names DEFAULT(NULL);

void VisuaInit();

class tVisuaRecord:private XmlWriter{
private:
	void WriteHeader();
	void WriteCoord();
	void WritePointData();
public:
	void WriteStructuredGrid(const char *filename);
};






#endif

