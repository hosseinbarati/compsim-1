#include <iostream>
#include <string.h>
extern int Nx, Ny, Nz, Nc;

int InitialRead(ifstream InputFile) {
	char str[MAX_STRING_LENGTH];
	register int i, j, k;
	GridBlock ***block;

	if (!InputFile.FileSearch("GRID")) TerM("No GRID keyword in the input file!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect GRID keyword format in the input file!");
	Nx = atoi(str);
	cout << "Nx=" << Nx << endl;
	if (!InputFile.ReadWord(str)) TerM("Incorrect GRID keyword format in the input file!");
	Ny = atoi(str);
	if (!InputFile.ReadWord(str)) TerM("Incorrect GRID keyword format in the input file!");
	Nz = atoi(str);

	if (!InputFile.FileSearch("NC")) TerM("No NC keyword in the input file!");
	if (!InputFile.ReadWord(str)) TerM("Incorrect NC keyword format in the input file!");
	PNc = atoi(str);
	if (!InputFile.ReadWord(str)) TerM("Incorrect NC keyword format in the input file!");
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
				block[i][j][k].ReadGridProperties(InputFile);
		}
	}
	

	return 0;
}
int File_Search(ifstream InputFile, char *rSeek) {
	register int i;
	char str[MAX_STRING_LENGTH];

	InputFile.clear();                 // clear fail and eof bits
	InputFile.seekg(0, std::ios::beg); // back to the start!
	*str='\0';
	do {
		i=Read_Word(InputFile, str);
		if (!strcmp(str, rSeek)) return -1;
	} while (i);

	return 0;		//Nothing found
}

int Read_Word(ifstream InputFile, char *rWord) {
	register unsigned char ch;
	register int i=0;

	*rWord='\0';
	do {
		InputFile.get(ch);
		if (InputFile.eof()) {
			return 0;		//Nothing has been read
		}
	} while ((ch<33) || (ch>126));

	while ((ch>32) && (ch<127)) {
		*(rWord+i)=ch;
		i++;
		InputFile.get(ch)
		if (InputFile.eof()) {
			*(rWord+i)='\0';
			return 1;		//Read, but end of file also encountered
		}		
	}

	*(rWord+i)='\0';
	return -1;		//Correct execution
}

