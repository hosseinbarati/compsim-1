#include <iostream>
#include <stdlib.h>

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#define MAX_STRING_LENGTH 256

typedef double FloatType;

int Nx, Ny, Nz;		//Number of grid blocks in X, Y, and Z directions
int Nc;				//Number of components
int PNc, UNc;		//Predefined Components, User Components

void TerM(char *ErrorMessage) {
	std::cout<<ErrorMessage<<std::endl;
	exit(-1):
}