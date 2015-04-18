#include <iostream>
#include <string.h>
#include "MIfstream.h"
extern int Nx, Ny, Nz, Nc;

int InitialRead(MIfstream &InputFile) {
	char str[MAX_STRING_LENGTH], str1[MAX_STRING_LENGTH];
	register int i, j, k;
	GridBlock ***block;
	FloatType TempL;

	if (!InputFile.FileSearch("GRID")) TerM("No GRID keyword in the input file!");
	if (!InputFile.ReadWord(str)) TerM("InputFile.ReadWord GRID keyword format in the input file!");
	Nx = atoi(str);
	cout << "Nx=" << Nx << endl;
	if (!InputFile.ReadWord(str)) TerM("InputFile.ReadWord GRID keyword format in the input file!");
	Ny = atoi(str);
	if (!InputFile.ReadWord(str)) TerM("InputFile.ReadWord GRID keyword format in the input file!");
	Nz = atoi(str);

	if (!InputFile.FileSearch("NC")) TerM("No NC keyword in the input file!");
	if (!InputFile.ReadWord(str)) TerM("InputFile.ReadWord NC keyword format in the input file!");
	PNc = atoi(str);
	if (!InputFile.ReadWord(str)) TerM("InputFile.ReadWord NC keyword format in the input file!");
	UNc = atoi(str);
	Nc = PNc + UNc;


	//Allocate GridBlock Objects
	block=new GridBlock**[Nx];
	for (i=0;i<Nx;i++) {
		block[i]=new GridBlock*[Ny];
		for (j=0;j<Ny;j++) block[i][j]=new GridBlock[Nz];
	}

	//Tell blocks to read their own properties from input
	for (k=0; k<Nz; k++) {
		for (j=0; j<Ny; j++) {
			for (i=0; i<Nx; i++) {
				block[i][j][k].SetIndex(i, j, k);
				//block[i][j][k].ReadGridProperties(InputFile);
			}
		}
	}

	if (!InputFile.FileSearch("DI")) TerM("No DI keyword in the input file!");
	if (!Read_Word(str)) TerM("InputFile.ReadWord DI keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=0; i<Nx; i++) { 
		if (!Read_Word(str1)) TerM("InputFile.ReadWord DI keyword format in the input file!");
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				GridBlock[i][j][k].SetDimX(atof(str1));
			}
		}
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(str1)) TerM("InputFile.ReadWord DI keyword format in the input file!");
		TempL=atof(str1);
		for (i=0; i<Nx; i++) {
			for (j=0; j<Ny; j++) {
				for (k=0; k<Nz; k++) {
					GridBlock[i][j][k].SetDimX(TempL);
				}
			}
		}
	else {
		TerM("InputFile.ReadWord DI keyword format in the input file!");
	}

	if (!InputFile.FileSearch(fp, "DJ")) TerM("No DJ keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("InputFile.ReadWord DJ keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=Nx; i<(Nx+Ny); i++) { 
		if (!Read_Word(fp, str1)) TerM("InputFile.ReadWord DJ keyword format in the input file!");
		gridDim[i]=atof(str1);
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("InputFile.ReadWord DJ keyword format in the input file!");
		tempL=atof(str1);
		for (i=Nx; i<(Nx+Ny); i++) gridDim[i]=tempL;
	}
	else {
		TerM("InputFile.ReadWord DJ keyword format in the input file!");
	}

	if (!InputFile.FileSearch(fp, "DK")) TerM("No DK keyword in the input file!");
	if (!Read_Word(fp, str)) TerM("InputFile.ReadWord DK keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (i=(Nx+Ny); i<(Nx+Ny+Nz); i++) { 
		if (!Read_Word(fp, str1)) TerM("InputFile.ReadWord DK keyword format in the input file!");
		gridDim[i]=atof(str1);
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("InputFile.ReadWord DK keyword format in the input file!");
		TempL=atof(str1);
		for (i=(Nx+Ny); i<(Nx+Ny+Nz); i++) gridDim[i]=tempL;
	}
	else {
		TerM("InputFile.ReadWord DK keyword format in the input file!");
	}

	//Porosity
	if (!InputFile.FileSearch("POR")) TerM("No POR keyword in the input file!");
	if (!Read_Word(str)) TerM("InputFile.ReadWord POR keyword format in the input file!");
	if (!strcmp(str, "VAR")) for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) { 
		if (!Read_Word(str1)) TerM("Invalid POR keyword format in the input file!");
		GridBlock[i][j][k].SetPorosity(atof(str1));
		//if (!porosity[i][j][k]) porosity[i][j][k]=1e-5;
	}
	else if (!strcmp(str, "CON")){
		if (!Read_Word(fp, str1)) TerM("InputFile.ReadWord POR keyword format in the input file!");
		tempL=atof(str1);
		for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) porosity[i][j][k]=tempL;
	}
	else if (!strcmp(str, "IVAR")){
		for (i=0; i<Nx; i++) {
			if (!Read_Word(fp, str1)) TerM("InputFile.ReadWord POR keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (j=0; j<Ny; j++) porosity[i][j][k]=tempL;
		}
	}
	else if (!strcmp(str, "JVAR")){
		for (j=0; i<Ny; j++) {
			if (!Read_Word(fp, str1)) TerM("InputFile.ReadWord POR keyword format in the input file!");
			tempL=atof(str1);
			for (k=0; k<Nz; k++) for (i=0; i<Nx; i++) porosity[i][j][k]=tempL;
		}
	}
	else if (!strcmp(str, "KVAR")){
		for (k=0; k<Nz; k++) {
			if (!Read_Word(fp, str1)) TerM("InputFile.ReadWord POR keyword format in the input file!");
			tempL=atof(str1);
			for (j=0; j<Ny; j++) for (i=0; i<Nx; i++) porosity[i][j][k]=tempL;
		}
	}
	else {
		TerM("InputFile.ReadWord POR keyword format in the input file!");
	}


	

	return 0;
}
