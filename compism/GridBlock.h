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


void GridBlock::SetDimI(){
	if (!InputFile.FileSearch("DI")) TerM("No DI keyword in the input file!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect DI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=0; i<Index+1; i++) { 
			if (!InputFile.ReadWord(str1)) TerM("Incorrect DI keyword format in the input file!");	
	}
	else if (!strcmp(str, "CON")){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect DI keyword format in the input file!");		
	}
	else {
		TerM("Incorrect DI keyword format in the input file!");
	}
	Dimension[0]=atof(str1);
}


int GridBlock::ReadGridProperties(ifstream InputFile) {
	char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
	register int i;
	MIfstream InputFile;
	

	if (!InputFile.FileSearch("DI")) TerM("No DI keyword in the input file!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect DI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=0; i<Index+1; i++) { 
			if (!InputFile.ReadWord(str1)) TerM("Incorrect DI keyword format in the input file!");	
	}
	else if (!strcmp(str, "CON")){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect DI keyword format in the input file!");		
	}
	else {
		TerM("Incorrect DI keyword format in the input file!");
	}
	Dimension[0]=atof(str1);

	if (!InputFile.FileSearch("DJ")) TerM("No DJ keyword in the input file!");
	if (!InputFile.ReadWord(str))  TerM("Incorrect DJ keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=0; i<Index+1; i++) { 
		if (!InputFile.ReadWord(str1))  TerM("Incorrect DJ keyword format in the input file!");		
	}
	else if (!strcmp(str, "CON")){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect DJ keyword format in the input file!");		
	}
	else {
		TerM("Incorrect DJ keyword format in the input file!");
	}
	Dimension[1]=atof(str1);

	if (!InputFile.FileSearch"DK") TerM("No DK keyword in the input file!");
	if (!InputFile.ReadWord(str))  TerM("Incorrect DK keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=0; i<Index+1; i++) { 
	if (!InputFile.ReadWord(str1))  TerM("Incorrect DK keyword format in the input file!");
	}
	else if (!strcmp(str, "CON")){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect DK keyword format in the input file!");		
	}
	else {
		TerM("Incorrect DK keyword format in the input file!");
	}
	Dimension[2]=atof(str1);
	
	//////porosity
	//////////////
	if (!InputFile.FileSearch("POR")) TerM("No POR keyword in the input file!");
	if (!InputFile.ReadWord(str))  TerM("Incorrect POR keyword format in the input file!");
	if (!strcmp(str1, "VAR")) for (i = 0; i < Index+1; i++) {
		if (!InputFile.ReadWord(str1))  TerM("Incorrect POR keyword format in the input file!");
	}
	else if (!strcmp(str1, "CON")){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect POR keyword format in the input file!");
	}
						// i < Ix
	else if (!strcmp(str1, "IVAR")) for (i = 0; i < (Index%Nx)+1; i++){
		if (!InputFile.ReadWord(str1)) TerM("Incorrect POR keyword format in the input file!");
	}
						// i < Ij
	else if (!strcmp(str1, "JVAR")) for (i = 0; i < (((Index - (Index%Nx)) / Nx) % Ny)+1; i++){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect POR keyword format in the input file!");
	}
						// i << Ik
	else if (!strcmp(str1, "KVAR")) for (i = 0; i < (Index - ((Index%Nx)*Nx) - ((((Index - (Index%Nx)) / Nx) % Ny)*Ny) / (Nx*Ny))+1; i++){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect POR keyword format in the input file!");
	}
	else {
	TerM("Incorrect POR keyword format in the input file!");
	}
	Porosity = atof(str1);
	
	//PERMEABILITY
	if (!InputFile.FileSearch("PERMI")) TerM("No PERMI keyword in the input file!");
	if (!InputFile.ReadWord(str))  TerM("Incorrect PERMI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i = 0; i<Index+1; i++) {
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMI keyword format in the input file!");
	}
	else if (!strcmp(str, "CON")) {
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMI keyword format in the input file!");
	}
	else if (!strcmp(str1, "IVAR")) for (i = 0; i < (Index%Nx)+1; i++){
		iif (!InputFile.Read_Word(str1))  TerM("Incorrect PERMI keyword format in the input file!");
	}
	else if (!strcmp(str1, "JVAR")) for (i = 0; i < (((Index - (Index%Nx)) / Nx) % Ny)+1; i++){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMI keyword format in the input file!");
	}
	else if (!strcmp(str1, "KVAR")) for (i = 0; i < (Index - ((Index%Nx)*Nx) - ((((Index - (Index%Nx)) / Nx) % Ny)*Ny) / (Nx*Ny))+1; i++){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMI keyword format in the input file!");
	}
	else {
		TerM("Incorrect PERMI keyword format in the input file!");
	}
	Permeability[0] = atof(str1);

	if (!InputFile.FileSearch("PERMJ")) TerM("No PERMJ keyword in the input file!");
	if (!InputFile.ReadWord(str))  TerM("Incorrect PERMJ keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i = 0; i<Index+1; i++) {
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMJ keyword format in the input file!");
	}
	else if (!strcmp(str, "CON")) {
		if (!ReadWord(InputFile, str1)) TerM("Incorrect PERMJ keyword format in the input file!");
	}
	else if (!strcmp(str1, "IVAR")) for (i = 0; i < (Index%Nx)+1; i++){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMJ keyword format in the input file!");
	}
	else if (!strcmp(str1, "JVAR")) for (i = 0; i < (((Index - (Index%Nx)) / Nx) % Ny)+1; i++){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMJ keyword format in the input file!");
	}
	else if (!strcmp(str1, "KVAR")) for (i = 0; i < (Index - ((Index%Nx)*Nx) - ((((Index - (Index%Nx)) / Nx) % Ny)*Ny) / (Nx*Ny))+1; i++){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMJ keyword format in the input file!");
	}
	else {
		TerM("Incorrect PERMJ keyword format in the input file!");
	}
	Permeability[1] = atof(str1);

	if (!InputFile.FileSearch("PERMK")) TerM("No PERMK keyword in the input file!");
	if (!InputFile.ReadWord(str))  TerM("Incorrect PERMK keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i = 0; i<Index+1; i++) {
	if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMK keyword format in the input file!");
	}
	else if (!strcmp(str, "CON")) {
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMK keyword format in the input file!");
	}
	else if (!strcmp(str1, "IVAR")) for (i = 0; i < (Index%Nx)+1; i++){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMK keyword format in the input file!");
	}
	else if (!strcmp(str1, "JVAR")) for (i = 0; i < (((Index - (Index%Nx)) / Nx) % Ny)+1; i++){
		if (!InputFile.ReadWord(str1)) TerM("Incorrect PERMK keyword format in the input file!");
	}
	else if (!strcmp(str1, "KVAR")) for (i = 0; i < (Index - ((Index%Nx)*Nx) - ((((Index - (Index%Nx)) / Nx) % Ny)*Ny) / (Nx*Ny))+1; i++){
		if (!InputFile.ReadWord(str1))  TerM("Incorrect PERMK keyword format in the input file!");
	}
	else {
		TerM("Incorrect PERMK keyword format in the input file!");
	}
	Permeability[2] = atof(str1);

	//////initial condition(pressure,water saturation,component present in the block)
	//case 2 needs to be completed
	if (!InputFile.FileSearch("INITCOND")) TerM("No INITCOND keyword in the input file!");
	if (!InputFile.ReadWord(str))  TerM("Incorrect INITCOND keyword format in the input file!");
	initCond = atoi(str);
	switch (initCond) {
		case 0:				//all
			if (!InputFile.FileSearch("IPRESS")) TerM("No IPRESS keyword in the input file!"); //initial pressure
			eclint = 0;
			for (i = 0; i < Index+1; i++) {
				if (!eclint) {
					if (!InputFile.ReadWord(str))  TerM("Incorrect IPRESS keyword format in the input file!");
					ECLStar(str, &eclint, &tempL);
					Pressure = tempL;//because p is the same for the first eclint gridblock 
							//so it just needs to evaluate one time for them 
				}
					eclint--;
			}
			
			
			if (!InputFile.FileSearch("IWS")) TerM("No IWS keyword in the input file!"); //Initial water saturation
			eclint = 0;
			for (i = 0; i < Index + 1; i++) {
				if (!eclint) {
					if (!InputFile.ReadWord(str)) TerM("Incorrect IWS keyword format in the input file!");
					ECLStar(str, &eclint, &tempL);
					Saturation[0] = tempL;
				}
				eclint--;
			}
			
			
			
			if (!InputFile.FileSearch("IGC")) TerM("No IGC keyword in the input file!"); //Initial global composition
			eclint = 0;
			for (i = 0; i < Index + 1; i++) for (int n = 0; n<Nc; n++)  {
				if (!eclint) {
					if (!InputFile.ReadWord(str)) TerM("Incorrect IGC keyword format in the input file!");
					ECLStar(str, &eclint, &tempL);
					Componenet[n] = templ;
				}
				eclint--;
			}
			AllFlash();  //we didn't define it yet
			break;
		case 1:				//same as 0 but depth variation only
			if (!FileSearch(InputFile, "IPRESS")) TerM("No IPRESS keyword in the input file!"); //initial pressure
			for (i = 0; i < (Index - ((Index%Nx)*Nx) - ((((Index - (Index%Nx)) / Nx) % Ny)*Ny) / (Nx*Ny))+1; i++){
				if (!InputFile.ReadWord(str)) TerM("Incorrect IPRESS keyword format in the input file!");
			}
			Pressure = atof(str);
			
			if (!InputFile.FileSearch("IWS")) TerM("No IWS keyword in the input file!"); //Initial water saturation
			for (i = 0; i < (Index - ((Index%Nx)*Nx) - ((((Index - (Index%Nx)) / Nx) % Ny)*Ny) / (Nx*Ny)) + 1; i++){
				if (!InputFile.ReadWord(str)) TerM("Incorrect IWS keyword format in the input file!");
			}
			Saturation[0] = atof(str);

			if (!InputFile.FileSearch("IGC")) TerM("No IGC keyword in the input file!"); //Initial global composition
			for (i = 0; i < (Index - ((Index%Nx)*Nx) - ((((Index - (Index%Nx)) / Nx) % Ny)*Ny) / (Nx*Ny)) + 1; i++)
				for (int n = 0; n<Nc; n++)  {
				if (!InputFile.ReadWord(str)) TerM("Incorrect IGC keyword format in the input file!");
				Componenet[n] = templ;
				}
			AllFlash();
			break;
				
		case 2: //compositional grading needs to be completed
	}
	return 0;
}
