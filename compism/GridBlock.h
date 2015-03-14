class GridBlock {
private:
	FloatType Porosity;
	FloatType Permeability[3];		//Absolute permeability in X, Y and Z directions in miliDarcy
	FloatType RelPerm[3];			//Relative permeability for water, oil (liquid hydrocarbon) and gas (gaseous hydrocarbon)
	FloatType Saturation[3];		//water, oil (liquid hydrocarbon) and gas (gaseous hydrocarbon) phase saturations
	FloatType Pressure;				//Block pressure in Pascal
	FloatType Dimension[3];			//Block length in X, Y and Z direction in meter
	FloatType *Componenet;			//Dynamic array for components present in the block
	int Index;						//Grid Block Index in Cartesian coordinate system

public:
	GridBlock();
	~GridBlock();
	void SetIndex(int, int, int);
	int ReadGridProperties(ifstream);
};

GridBlock::GridBlock(void) {
	Componenet=new FloatType[Nc];
}

GridBlock::~GridBlock(void){
	delete[] Componenet;
	Componenet=NULL;
}

void GridBlock::SetIndex(int Ix, int Iy, int Iz) {
	Index=Iz*(Ny*Nx)+Iy*Nx+Ix;
}

int GridBlock::ReadGridProperties(ifstream InputFile) {
	char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
	register int i;
	

	if (!File_Search(InputFile, "DI")) TerM("No DI keyword in the input file!");
	if (!Read_Word(InputFile, str)) TerM("Incorrect DI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=0; i<Index; i++) { 
		if (!Read_Word(InputFile, str1)) TerM("Incorrect DI keyword format in the input file!");		
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(InputFile, str1)) TerM("Incorrect DI keyword format in the input file!");		
	}
	else {
		TerM("Incorrect DI keyword format in the input file!");
	}
	Dimension[0]=atof(str1);

	if (!File_Search(InputFile, "DJ")) TerM("No DJ keyword in the input file!");
	if (!Read_Word(InputFile, str)) TerM("Incorrect DJ keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=0; i<Index; i++) { 
		if (!Read_Word(InputFile, str1)) TerM("Incorrect DJ keyword format in the input file!");		
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(InputFile, str1)) TerM("Incorrect DJ keyword format in the input file!");		
	}
	else {
		TerM("Incorrect DJ keyword format in the input file!");
	}
	Dimension[1]=atof(str1);

	if (!File_Search(InputFile, "DK")) TerM("No DK keyword in the input file!");
	if (!Read_Word(InputFile, str)) TerM("Incorrect DK keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=0; i<Index; i++) { 
		if (!Read_Word(InputFile, str1)) TerM("Incorrect DK keyword format in the input file!");		
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(InputFile, str1)) TerM("Incorrect DK keyword format in the input file!");		
	}
	else {
		TerM("Incorrect DK keyword format in the input file!");
	}
	Dimension[2]=atof(str1);
	
	
	//////porosity
	//////////////
if (!File_Search(InputFile, "POR")) TerM("No POR keyword in the input file!");
if (!Read_Word(InputFile, str1)) TerM("Incorrect POR keyword format in the input file!");
if (!strcmp(str1, "VAR")) for (i = 0; i < Index; i++) {
if (!Read_Word(InputFile, str1)) TerM("Incorrect POR keyword format in the input file!");
}
else if (!strcmp(str1, "CON")){
if (!Read_Word(InputFile, str1)) TerM("Incorrect POR keyword format in the input file!");
}
					// i < Ix
else if (!strcmp(str1, "IVAR")) for (i = 0; i < (Index%Nx); i++){
if (!Read_Word(InputFile, str1)) TerM("Incorrect POR keyword format in the input file!");
}
					// i < Iy
else if (!strcmp(str1, "JVAR")) for (i = 0; i < (((Index - (Index%Nx)) / Nx) % Ny); i++){
if (!Read_Word(InputFile, str1)) TerM("Incorrect POR keyword format in the input file!");
}
					// i << Iz
else if (!strcmp(str1, "KVAR")) for (i = 0; i < (Index - ((Index%Nx)*Nx) - (((((Index - (Index%Nx)) / Nx) % Ny)*Ny) / (Nx*Ny)); i++){
if (!Read_Word(InputFile, str1)) TerM("Incorrect POR keyword format in the input file!");
}
else {
TerM("Incorrect POR keyword format in the input file!");
}
Porosity = atof(str1);

	
	return 0;
}
