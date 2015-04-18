#include"flash.h"

class GridBlock {
private:
	FloatType Porosity;
	FloatType Permeability[3];		//Absolute permeability in X, Y and Z directions in miliDarcy
	FloatType RelPerm[3];			//Relative permeability for water, oil (liquid hydrocarbon) and gas (gaseous hydrocarbon)
	FloatType Saturation[3];		//water, oil (liquid hydrocarbon) and gas (gaseous hydrocarbon) phase saturations
	FloatType Pressure;				//Block pressure in Pascal
	FloatType Dimension[3];			//Block length in X, Y and Z direction in meter
	FloatType WOCHeight;
	FloatType *MW;
	FloatType *TCRIT;
	FloatType *PCRIT;
	FloatType *VCRIT;
	FloatType *AC;
	FloatType *PARACHOR;
	FloatType *Componenet;			//Dynamic array for components present in the block
	FloatType *P;
	double resTemp;
	int refL;
	int Index;						//Grid Block Index in Cartesian coordinate system

public:
	GridBlock();
	~GridBlock();
	void SetIndex(int, int, int);
	int ReadGridProperties(ifstream);
	void SetDimX(FloatType);
};

GridBlock::GridBlock(void) {
	Componenet=new FloatType[Nc+1];
	MW = new FloatType[Nc+1];
	TCRIT = new FloatType[Nc+1];
	PCRIT = new FloatType[Nc+1];
	VCRIT = new FloatType[Nc+1];
	AC = new FloatType[Nc+1];
	PARACHOR = new FloatType[Nc+1];
}

GridBlock::~GridBlock(void){
	delete[] Componenet;
	Componenet=NULL;
	delete[] TCRIT;
	TCRIT = NULL;
	delete[] VCRIT;
	VCRIT = NULL;
	delete[] AC;
	AC = NULL;
	delete[] PARACHOR;
	PARACHOR = NULL;
	delete[] MW;
	MW = NULL;
	delete[] P;
	P = NULL;
}

void GridBlock::SetIndex(int Ix, int Iy, int Iz) {
	Index=Iz*(Ny*Nx)+Iy*Nx+Ix;
}

void GridBlock::SetDimX(FloatType Dx) {
	Dimension[0]=Dx;
}
void GridBlock::SetDimY(FloatType Dy) {
	Dimension[1]=Dy;
}
void GridBlock::SetDimZ(FloatType Dz) {
	Dimension[2]=Dz;
}

void GridBlock::SetPorosity(FloatType Por) {
	Porosity=Por;
}

/*int GridBlock::ReadGridProperties(ifstream InputFile) {
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
	
	//Reservoir temp.
	if (!InputFile.FileSearch("RESTEMP")) TerM("No RESTEMP keyword in the input file!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect RESTEMP keyword format in the input file!"); {
		resTemp = atof(str) + 273.15;
	}
	
	/////////////////////////////////
	if (!InputFile.FileSearch("COMPNAME")) TerM("No COMPNAME keyword in the input file!");
	for (i = 0; i < PNc; i++){
		if (!InputFile.ReadWord(str)) TerM("Incorrect COMPNAME  keyword format in the input file!");
		if (!strcmp(str, "C1")) {
			MW[i] = 16.043;
			TCRIT[i] = 190.56;
			PCRIT[i] = 4599000;
			VCRIT[i] = 98.6;
			AC[i] = 0.008; //edited
			PARACHOR[i] = 74.05;
		}
		else if (!strcmp(str, "C2")) {
			MW[i] = 30.070;
			TCRIT[i] = 305.32;
			PCRIT[i] = 4872000;
			VCRIT[i] = 145.5;
			AC[i] = 0.0995;
			PARACHOR[i] = 112.91;
		}
		else if (!strcmp(str, "C3")) {
			MW[i] = 44.096;
			TCRIT[i] = 369.83;
			PCRIT[i] = 4248000;
			VCRIT[i] = 200;
			AC[i] = 0.1523;
			PARACHOR[i] = 154.03;
		}
		else if (!strcmp(str, "iC4")) {
			MW[i] = 58.123;
			TCRIT[i] = 408.14;
			PCRIT[i] = 3648000;
			VCRIT[i] = 262.7;
			AC[i] = 0.1770;
			PARACHOR[i] = 185.32;
		}
		else if (!strcmp(str, "nC4")) {
			MW[i] = 58.123;
			TCRIT[i] = 425.12;
			PCRIT[i] = 3796000;
			VCRIT[i] = 255;
			AC[i] = 0.2002;
			PARACHOR[i] = 193.90;
		}
		else if (!strcmp(str, "iC5")) {
			MW[i] = 72.150;
			TCRIT[i] = 460.43;
			PCRIT[i] = 3381000;
			VCRIT[i] = 305.8;
			AC[i] = 0.2275;
			PARACHOR[i] = 229.37;
		}
		else if (!strcmp(str, "nC5")) {
			MW[i] = 72.150;
			TCRIT[i] = 469.7;
			PCRIT[i] = 3370000;
			VCRIT[i] = 313;
			AC[i] = 0.2515;
			PARACHOR[i] = 236.00;
		}
		else if (!strcmp(str, "nC6")) {
			MW[i] = 86.177;
			TCRIT[i] = 507.6;
			PCRIT[i] = 3025000;
			VCRIT[i] = 371;
			AC[i] = 0.3013;
			PARACHOR[i] = 276.71;
		}
		else if (!strcmp(str, "nC7")) {
			MW[i] = 100.204;
			TCRIT[i] = 540.2;
			PCRIT[i] = 2740000;
			VCRIT[i] = 428;
			AC[i] = 0.3495;
			PARACHOR[i] = 318.44;
		}
		else if (!strcmp(str, "nC8")) {
			MW[i] = 114.231;
			TCRIT[i] = 568.7;
			PCRIT[i] = 2490000;
			VCRIT[i] = 486;
			AC[i] = 0.3996;
			PARACHOR[i] = 359.33;
		}
		else if (!strcmp(str, "nC9")) {
			MW[i] = 128.258;
			TCRIT[i] = 594.6;
			PCRIT[i] = 2290000;
			VCRIT[i] = 544;
			AC[i] = 0.4435;
			PARACHOR[i] = 399.57;
		}
		else if (!strcmp(str, "nC10")) {
			MW[i] = 142.285;
			TCRIT[i] = 617.7;
			PCRIT[i] = 2110000;
			VCRIT[i] = 600;
			AC[i] = 0.4923;
			PARACHOR[i] = 440.69;
		}
		else if (!strcmp(str, "nC11")) {
			MW[i] = 156.312;
			TCRIT[i] = 639;
			PCRIT[i] = 1949000;
			VCRIT[i] = 659;
			AC[i] = 0.5303;
			PARACHOR[i] = 482.00;
		}
		else if (!strcmp(str, "nC12")) {
			MW[i] = 170.338;
			TCRIT[i] = 658;
			PCRIT[i] = 1820000;
			VCRIT[i] = 716;
			AC[i] = 0.5764;
			PARACHOR[i] = 522.26;
		}
		else if (!strcmp(str, "nC13")) {
			MW[i] = 184.365;
			TCRIT[i] = 675;
			PCRIT[i] = 1680000;
			VCRIT[i] = 775;
			AC[i] = 0.6174;
			PARACHOR[i] = 536.77;
		}
		else if (!strcmp(str, "nC14")) {
			MW[i] = 198.392;
			TCRIT[i] = 693;
			PCRIT[i] = 1570000;
			VCRIT[i] = 830;
			AC[i] = 0.6430;
			PARACHOR[i] = 606.05;
		}
		else if (!strcmp(str, "nC15")) {
			MW[i] = 212.419;
			TCRIT[i] = 708;
			PCRIT[i] = 1480000;
			VCRIT[i] = 889;
			AC[i] = 0.6863;
			PARACHOR[i] = 647.43;
		}
		else if (!strcmp(str, "nC16")) {
			MW[i] = 226.446;
			TCRIT[i] = 723;
			PCRIT[i] = 1400000;
			VCRIT[i] = 944;
			AC[i] = 0.7174;
			PARACHOR[i] = 688.50;
		}
		else if (!strcmp(str, "nC17")) {
			MW[i] = 240.473;
			TCRIT[i] = 736;
			PCRIT[i] = 1340000;
			VCRIT[i] = 1000;
			AC[i] = 0.7697;
			PARACHOR[i] = 730.05;
		}
		else if (!strcmp(str, "nC18")) {
			MW[i] = 254.5;
			TCRIT[i] = 747;
			PCRIT[i] = 1270000;
			VCRIT[i] = 1060;
			AC[i] = 0.8114;
			PARACHOR[i] = 771.95;
		}
		else if (!strcmp(str, "nC19")) {
			MW[i] = 268.527;
			TCRIT[i] = 758;
			PCRIT[i] = 1210000;
			VCRIT[i] = 1120;
			AC[i] = 0.8522;
			PARACHOR[i] = 813.85;
		}
		else if (!strcmp(str, "nC20")) {
			MW[i] = 282.553;
			TCRIT[i] = 768;
			PCRIT[i] = 1160000;
			VCRIT[i] = 1170;
			AC[i] = 0.9069;
			PARACHOR[i] = 853.67;
		}
		else if (!strcmp(str, "nC21")) {
			MW[i] = 296.580;
			TCRIT[i] = 781.7;
			PCRIT[i] = 1147000;
			VCRIT[i] = 1198;
			AC[i] = 0.9220;
			PARACHOR[i] = 897.64;
		}
		else if (!strcmp(str, "nC22")) {
			MW[i] = 310.610;
			TCRIT[i] = 791.8;
			PCRIT[i] = 1101000;
			VCRIT[i] = 1253;
			AC[i] = 0.9550;
			PARACHOR[i] = 939.55;
		}
		else if (!strcmp(str, "nC23")) {
			MW[i] = 324.630;
			TCRIT[i] = 801.3;
			PCRIT[i] = 1059000;
			VCRIT[i] = 1307;
			AC[i] = 0.9890;
			PARACHOR[i] = 981.43;
		}
		else if (!strcmp(str, "nC24")) {
			MW[i] = 338.680;
			TCRIT[i] = 810.4;
			PCRIT[i] = 1019000;
			VCRIT[i] = 1362;
			AC[i] = 1.0190;
			PARACHOR[i] = 1023.40;
		}
		else if (!strcmp(str, "CO2")) {
			MW[i] = 44.010;
			TCRIT[i] = 304.19;
			PCRIT[i] = 7382000;
			VCRIT[i] = 94;
			AC[i] = 0.2276;
			PARACHOR[i] = 82.00;
		}
		else if (!strcmp(str, "O2")) {
			MW[i] = 31.999;
			TCRIT[i] = 154.58;
			PCRIT[i] = 5043000;
			VCRIT[i] = 73.4;
			AC[i] = 0.0218;
		}
		else if (!strcmp(str, "N2")) {
			MW[i] = 28.014;
			TCRIT[i] = 126.1;
			PCRIT[i] = 3394000;
			VCRIT[i] = 90.1;
			AC[i] = 0.0403;
			PARACHOR[i] = 61.12;
		}
		else if (!strcmp(str, "H2S")) {
			MW[i] = 34.082;
			TCRIT[i] = 373.53;
			PCRIT[i] = 8963000;
			VCRIT[i] = 98.5;
			AC[i] = 0.0827;
			PARACHOR[i] = 85.50;
		}
		else if (!strcmp(str, "SO2")) {
			MW[i] = 64.065;
			TCRIT[i] = 430.75;
			PCRIT[i] = 7884000;
			VCRIT[i] = 122;
			AC[i] = 0.2451;
		}
		else if (!strcmp(str, "H2")) {
			MW[i] = 2.016;
			TCRIT[i] = 33.18;
			PCRIT[i] = 1313000;
			VCRIT[i] = 64.2;
			AC[i] = 0.2150;
		}
		else if (!strcmp(str, "H2O")) {
			MW[i] = 18.015;
			TCRIT[i] = 647.13;
			PCRIT[i] = 22055000;
			VCRIT[i] = 56;
			AC[i] = 0.3449;
		}
		else {
			TerM("Unknown Component!");
		}

	}

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
		if (!InputFile.FileSearch("REFLAYER")) TerM("No REFLAYER keyword in the input file!");
		if (!InputFile.ReadWord(str)) TerM("Incorrect REFLAYER keyword format in the input file!");
		refL = (atoi(str) - 1);

		if (!InputFile.FileSearch("REFPRES")) TerM("No REFPRES keyword in the input file!");
		if (!InputFile.ReadWord(str)) TerM("Incorrect REFPRES keyword format in the input file!");
		P[refL] = atof(str);

		if (!InputFile.FileSearch("WOCHEIGHT")) TerM("No WOCHEIGHT keyword in the input file!");
		if (!InputFile.ReadWord(str)) TerM("Incorrect WOCHEIGHT keyword format in the input file!");
		WOCHeight = atof(str);

		if (!InputFile.FileSearch("REFCOMP")) TerM("No REFCOMP keyword in the input file!");
		for (i = 0; i < Nc; i++){
		if (!InputFile.ReadWord(str)) TerM("Incorrect REFCOMP keyword format in the input file!");
		composition[i] = atof(str);
	}
	firstflash(0.5, composition, AC, resTemp, TCRIT, PCRIT, P[refL], true);
	return 0;
}*/
