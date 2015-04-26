#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef _WIN32
#include <CL/cl.h>
#else
#include <CL\cl.h>
#endif
#include <time.h>

//#define USE_FLOATS


#define MAX_NR_CYCLE 10
#define MAX_STRING_LENGTH 80
#define WELL_I 5		//well integer properties
//#define WELL_F 3		//well FType properties

#define COMP_PROPS 14
#define SG 0
#define TB 1
#define MW 2
#define AC 3
#define PCRIT 4
#define TCRIT 5
#define VCRIT 6
#define TR 7
#define MU 8
#define EOS_A 9
#define EOS_B 10
#define COLDIA 11
#define EPDIA 12
#define PARACHOR 13

#define BLOCK_F_PROPS 4
#define RO 0		//molar density
#define BLOCK_PC 1
#define BMU 2
#define BGAM 3



#define SAT_TABLE 4
#define SAT 0
#define KR 1
#define KRO 2
#define PC 3
#define MAX_FLASH_N_ITR 100

#ifdef USE_FLOATS

#define FType float

#define G_ACC 9.80665f

#define WAT_M_RO 55.5084f
#define SI_C 9.869232667160128e-13f		// (1e-7/101325)
//Volume=Cubic meter
//Area=Square meter
//Permeability=miliDarcy
//Length=Meter
//Viscosity=centipoise
//pressure=Pascal
//Rate=Cubic meter per second
//Time=Second
#define RGAS 8.3144621f		//Universal gas constant
#define PI 3.1415926535897932384626433f
#define SQRT2 1.4142135623730950488016887f
#define EPSILON 0.001f

#define BTOL 1e-5f

#define XMTOL 0.0005f

#define EPS_SAT 0.001f		//Saturation @ bubble and dew
/*#define DKRGDSG 0;
#define DKRODSO 1;
#define DKRODSG 2;
#define DKRWDS 3;
#define DPCGDSG 4;
#define DPCWDS 5;*/

#define LN10 2.302585092994046f

#define NEWTON_TOL 1e-8f
#define FLASH_TOL 1e-8f

#define DZERO 1e-20f

#define TOLMAIN 1e-5f

#define TIMECHFACT 3600.0f		//Hour

//#define NO_ROOT_L 2000.0f

#define PZERO 1000000.0f

#define PRESSTOL 1000.0f

#define SXMTOL 1e-10f

#define MAXDP 0.01f
#define MAXPP 1e6f
#define MINDT 1000.0f

#define TIME_STEP_INC 1.05f
#define TIME_STEP_DEC 1.7f


#else


#define FType double

#define G_ACC 9.80665
#define G_ACCG 9.80665e-3

#define WAT_M_RO 55.5084
#define SI_C 9.869232667160128e-13		// (1e-7/101325)
//Volume=Cubic meter
//Area=Square meter
//Permeability=miliDarcy
//Length=Meter
//Viscosity=centipoise
//pressure=Pascal
//Rate=Cubic meter per second
//Time=Second
#define RGAS 8.3144621		//Universal gas constant
//#define BOLTZ 1.380648813e-16		//Boltzmann's Constant (erg/K)
#define PI 3.1415926535897932384626433
#define SQRT2 1.4142135623730950488016887
#define EPSILON 0.001

#define BTOL 1e-100

#define XMTOL 0.0005

#define EPS_SAT 1e-7		//Saturation @ bubble and dew
/*#define DKRGDSG 0;
#define DKRODSO 1;
#define DKRODSG 2;
#define DKRWDS 3;
#define DPCGDSG 4;
#define DPCWDS 5;*/

#define LN10 2.302585092994046

#define NEWTON_TOL 1e-8
#define FLASH_TOL 1e-8

#define DZERO 1e-20

#define TOLMAIN 1e-2

#define TIMECHFACT 3600		//Hour

//#define NO_ROOT_L 10

#define PZERO 1000000

#define PRESSTOL 1000

#define SXMTOL 1e-14

#define MAXDP 1e6
#define MAXPP 5e6
#define MINDT 1e-13
#define MAXDT 1e10

#define TIME_STEP_INC 1.05
#define TIME_STEP_DEC 1.5
#define TSTEPTOL 1
#define TSTEPIF 1e-5
#define RELPERM0 0	//1e-11
#define MAXINITER 5000
#define NTIMESTERINC 10
#define TRIV_K 1e-5
#define MAX_MNR_ITER 100

#define MRINERITER 3
#define MRMAXNONZERO 128
#define BICGSTABTOL	1e-30

#define PSGROUP_NO 56


#define IFTBASECASE 50		//mN/m
#define IFTPOWER 7
#define MISCIBLEIFT 1		//mN/m

#define GPU_TEMP_THRESHOLD 80
#define CPU_TEMP_THRESHOLD 80

//#define pCoeff 1e7

#endif





///////////////////GLOBALS//////////////////////////////////
int Nx, Ny, Nz;		//Reservoir Dimension
FType *gridDim;	//Grid size
FType ***porosity;		//Porosity Distribution
FType ****perm;		//Permeability tensor distribution
FType refP, cpor, dcpor;		//Compressibility
int PNc, UNc, Nc;		//Number of Components
unsigned char PR, SRK;		//EOS Type, boolean
FType **fluidProp;			//Fluids components properties
int Nswt, Nsgt;		//Number of saturation tables inputs
FType **swt, **sgt;			//Saturation tables
int initCond;		//Initial Condition State
FType ****sat, ****bsat;		//Saturation and its backup
FType ****P;		//Pressure
FType *****comp;		//Phases Composition
int refL;		//Reference layer in initial condition
int wellNO;		//Number of wells
int **welli;	//integer well data: X-Well, Y_Well, Well_Type
FType **wellf;	//Double well data: Pwf or Q
FType ****IFTran;		//Permeability average on blocks boundary, kA/h
FType resTemp;		//Reservoir Temperature
FType *****blockFProps;		//Calculated fluids properties in each gridblock
FType *****trans;		//Transmissibility
FType ****Wtran;		//Water Transmissibility
FType ****relPerm;		//Relative Permeability
FType watRo, watMu;	//Water properties
//FType **jac;		//Jacobian matrix
FType *ans;		//Answer matrix
FType ****preProp;		//Multiplication of properties in previous time step
FType Dt;		//Timestep size, seconds
FType *****dE;		//Phase density derivative
char *****transS;		//upstream weighting condition
FType **satJac, *satAns, *Xm, *Xms;	//Minor NR Matrixes
char ***phaseStat;		//Number and type of the phases in GBs, 1->L, -1->G, 0->2Ph
//int curSize;		//Jacobian Matrix current size
//FType *Unk;		//Unknowns matrix
FType WOCHeight;	//WOC
FType totalTime;	//total simulation time, days input, seconds change
FType *Ki;
FType QoT, QgT, QwT;		//Well flow rates
FType cumQo, cumQg, cumQw;		//cummulative rates
FType ***blockH;		//Block Heights
int incCount=0;		//Time Step Control Counter
FType *****bcomp;		//Composition Backup
FType ****bP;		//Pressure Backup
FType Eps;		//Machine Epsilon
FType *CSRjac;	//Compressed Sparse Row Jacobian
int *CSRrow, *CSRcol;		//Compressed sparse row location holders
int bCSRSize;		//Compressed Sparse Row Size for a block
//int jacIndex;		//Compressed Sparse Row Jacobian Index Holder
//int CSRrowIndex;		//Compressed Sparse Row Index Holder
int ***pJHolder;		//Jacobian Position Holder: Parallel design
FType pCoeff;		//Mean P for Preconditioning
FType *preCon;		//Preconditioner
int *preConRow;		//Preconditioner row holder
int *preConIndex;		//Preconditioner index holder
cl_platform_id *platformIds;
cl_device_id deviceId;
cl_context context;
cl_program program;
cl_kernel xZero, xZeroF, init, f1, MatMul, f2, f3, clBuildPrecon, dot_persist_kernel, MatMulT, Pf1, Pf2, Pf3, Pf4, dot_persist_kernelF, MatMulTF, MatMulF, ClScaling, clBinPartialSort;
cl_command_queue commandQueue;
//char chWellStat=0;
int totalJac, totalRow;
FType **bic;
FType *TStepMarker;
int TSMSize;
char repStat=0;
FType **STcomp;
FType wellrate[50];
char ***bphaseStat;
int TSMCtrl=0;
char maxNR=0;
FType ****dRelPerm;
FType ****dWtran;
char gasInjStat;
//FType **fullPre;
int *preConCSRrow;
int *preConCSRcol;
double solverTimeAcc=0;
double solverTime=0;
clock_t BICGStart, BICGEnd;
FType ***tor, *****diffusion;
FType Qw, Qo, Qg;
FType wGLR;
char buildPreconFlag=-1;
cl_mem clM, clMCSRrow, clMCSRcol;
FType ***Bift;
FType QgCumInj=0, QgInj;
FType QoCumProd=0, QoProd;
FType QoProduced, QgInjected, SumQoProduced=0, SumQgInjected=0;
FType nesbat=0;

long int nprocs;
int *CPUTempArray;



void TerM(char *str);
FType dSgn(FType);
void CLFinish(void);

FType Solve_Z(FType a1, FType a2, FType a3, char state) {
	FType Q, J, D, x1, x2, x3, t, t1, Z;
	FType a13, qsqrtemp;

	Q=(3*a2-a1*a1)/9;
	J=(9*a1*a2-27*a3-2*a1*a1*a1)/54;
	D=Q*Q*Q+J*J;

	if (D<0) {
		qsqrtemp=2*sqrt(-Q);
		a13=a1/3;
		t=acos(J/sqrt(-Q*Q*Q))/3;
		x1=qsqrtemp*cos(t)-a13;
		x2=qsqrtemp*cos(t+2*PI/3)-a13;
		x3=qsqrtemp*cos(t+4*PI/3)-a13;
		if ((state=='g') || (state=='G')){
			if ((x1>=x2) && (x1>=x3)) Z=x1;
			else if ((x2>=x1) && (x2>=x3)) Z=x2;
			else Z=x3;
		}
		else {		//Liquid state
			if ((x1<=x2) && (x1<=x3) && (x1>0)) Z=x1;
			else if ((x2<=x1) && (x2<=x3) && (x2>0)) Z=x2;
			else Z=x3;
		}
	}
	else if (D==0) {
		if (J>0) t=pow(J, 1.0/3);
		else t=-pow(-J,1.0/3);
		x1=2*t-a1/3;
		x2=-t-a1/3;
		if ((state=='g') || (state=='G')) {
			if (x2>=x1) Z=x2;
			else Z=x1;
		}
		else {
			if ((x2<=x1) && (x2>0)) Z=x2;
			else Z=x1;
		}		
	}
	else {
		t1=sqrt(D);
		t=J+t1;
		if (t>0) Z=pow(t, 1.0/3);
		else Z=-pow(-t,1.0/3);
		t=J-t1;
		if (t>0) Z+=pow(t, 1.0/3);
		else Z-=pow(-t,1.0/3);
		Z-=a1/3;
	}

	return Z;
}

void TerM(char *str) {
	register int i, j, k, n;

	puts(str);

	free(gridDim);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) free(porosity[i][j]);
		free(porosity[i]);
	}
	free(porosity);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) free(Bift[i][j]);
		free(Bift[i]);
	}
	free(Bift);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) free(tor[i][j]);
		free(tor[i]);
	}
	free(tor);

	for (i=0; i<(Nx+2); i++) {
		for (j=0; j<(Ny+2); j++) {
			for(k=0; k<(Nz+2); k++) free(P[i][j][k]);
			free(P[i][j]);
		}
		free(P[i]);
	}
	free(P);


	for (i=0; i<(Nx+2); i++) {
		for (j=0; j<(Ny+2); j++) {
			for(k=0; k<(Nz+2); k++) free(bP[i][j][k]);
			free(bP[i][j]);
		}
		free(bP[i]);
	}
	free(bP);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) free(perm[i][j][k]);
			free(perm[i][j]);
		}
		free(perm[i]);
	}
	free(perm);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) free(sat[i][j][k]);
			free(sat[i][j]);
		}
		free(sat[i]);
	}
	free(sat);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) free(bsat[i][j][k]);
			free(bsat[i][j]);
		}
		free(bsat[i]);
	}
	free(bsat);

	for (i=0; i<Nc; i++) {
		free(fluidProp[i]);
	}
	free(fluidProp);

	for (i=0; i<Nswt; i++) {
		free(swt[i]);
	}
	free(swt);

	for (i=0; i<Nsgt; i++) {
		free(sgt[i]);
	}
	free(sgt);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) {
				for(n=0; n<Nc; n++) free(comp[i][j][k][n]);
				free(comp[i][j][k]);
			}
			free(comp[i][j]);
		}
		free(comp[i]);
	}
	free(comp);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) {
				for(n=0; n<Nc; n++) free(diffusion[i][j][k][n]);
				free(diffusion[i][j][k]);
			}
			free(diffusion[i][j]);
		}
		free(diffusion[i]);
	}
	free(diffusion);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) {
				for(n=0; n<Nc; n++) free(bcomp[i][j][k][n]);
				free(bcomp[i][j][k]);
			}
			free(bcomp[i][j]);
		}
		free(bcomp[i]);
	}
	free(bcomp);

	for (i=0; i<wellNO; i++) {
		free(welli[i]);
		free(wellf[i]);
	}
	free(welli);
	free(wellf);
	
	
	//////////////////////////////////////////////////
	//////////////NON_INPUT///////////////////////////
	/////////////////////////////////////////////////
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) free(IFTran[i][j][k]);
			free(IFTran[i][j]);
		}
		free(IFTran[i]);
	}
	free(IFTran);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) {
				for(n=0; n<BLOCK_F_PROPS; n++) free(blockFProps[i][j][k][n]);
				free(blockFProps[i][j][k]);
			}
			free(blockFProps[i][j]);
		}
		free(blockFProps[i]);
	}
	free(blockFProps);

	for (i=0; i<(Nx+1); i++) {
		for (j=0; j<(Ny+1); j++) {
			for(k=0; k<(Nz+1); k++) {
				for(n=0; n<12; n++) {
					free(trans[i][j][k][n]);
				}
				free(trans[i][j][k]);
			}
			free(trans[i][j]);
		}
		free(trans[i]);
	}
	free(trans);

	for (i=0; i<(Nx+1); i++) {
		for (j=0; j<(Ny+1); j++) {
			for(k=0; k<(Nz+1); k++) {
				for(n=0; n<3; n++) {
					free(transS[i][j][k][n]);
				}
				free(transS[i][j][k]);
			}
			free(transS[i][j]);
		}
		free(transS[i]);
	}
	free(transS);

	for (i=0; i<(Nx+1); i++) {
		for (j=0; j<(Ny+1); j++) {
			for(k=0; k<(Nz+1); k++) free(Wtran[i][j][k]);
			free(Wtran[i][j]);
		}
		free(Wtran[i]);
	}
	free(Wtran);

	for (i=0; i<(Nx+1); i++) {
		for (j=0; j<(Ny+1); j++) {
			for(k=0; k<(Nz+1); k++) free(dWtran[i][j][k]);
			free(dWtran[i][j]);
		}
		free(dWtran[i]);
	}
	free(dWtran);
	
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) free(relPerm[i][j][k]);
			free(relPerm[i][j]);
		}
		free(relPerm[i]);
	}
	free(relPerm);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) free(dRelPerm[i][j][k]);
			free(dRelPerm[i][j]);
		}
		free(dRelPerm[i]);
	}
	free(dRelPerm);


	/*for (i=0; i<(Nx*Ny*Nz*(2*Nc+4)); i++) {
		free(preCon[i]);
	}
	free(preCon);*/

	free(CSRjac);
	free(CSRrow);
	free(CSRcol);

	free(preCon);

	free(ans);
	free(preConRow);
	free(preConIndex);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) free(pJHolder[i][j]);
		free(pJHolder[i]);
	}
	free(pJHolder);

	//free(Unk);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) free(preProp[i][j][k]);
			free(preProp[i][j]);
		}
		free(preProp[i]);
	}
	free(preProp);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for(k=0; k<Nz; k++) {
				for(n=0; n<(Nc+1); n++) free(dE[i][j][k][n]);
				free(dE[i][j][k]);
			}
			free(dE[i][j]);
		}
		free(dE[i]);
	}
	free(dE);

	////////////////////////////////////////////////////

	for (i=0; i<(Nc+1); i++) {
		free(satJac[i]);
	}
	free(satJac);

	free(satAns);
	free(Xm);
	free(Xms);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) free(phaseStat[i][j]);
		free(phaseStat[i]);
	}
	free(phaseStat);

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) free(bphaseStat[i][j]);
		free(bphaseStat[i]);
	}
	free(bphaseStat);

	free(Ki);

	for (i=0; i<(Nx+2); i++) {
		for (j=0; j<(Ny+2); j++) free(blockH[i][j]);
		free(blockH[i]);
	}
	free(blockH);

	for (i=0; i<Nc; i++) {
		free(bic[i]);
	}
	free(bic);

	free(TStepMarker);

	for (i=0; i<Nc; i++) {
		free(STcomp[i]);
	}
	free(STcomp);

	free(preConCSRrow);
	free(preConCSRcol);

	free(CPUTempArray);

	//CLFinish();

	//system("shutdown -s -f");

	exit(0);
}

FType dSgn(FType x) {
	if (x) return x/fabs(x);
	else return 0;
}

void CLFinish(void) {
	free(platformIds);
	if (xZero) clReleaseKernel(xZero);
	if (xZeroF) clReleaseKernel(xZeroF);
	if (init) clReleaseKernel(init);
	if (f1) clReleaseKernel(f1);
	if (MatMul) clReleaseKernel(MatMul);
	if (MatMulT) clReleaseKernel(MatMulT);
	if (MatMulTF) clReleaseKernel(MatMulTF);
	if (MatMulF) clReleaseKernel(MatMulF);
	if (f2) clReleaseKernel(f2);
	if (f3) clReleaseKernel(f3);
	if (clBuildPrecon) clReleaseKernel(clBuildPrecon);
	if (dot_persist_kernel) clReleaseKernel(dot_persist_kernel);
	if (dot_persist_kernelF) clReleaseKernel(dot_persist_kernelF);
	if (Pf1) clReleaseKernel(Pf1);
	if (Pf2) clReleaseKernel(Pf2);
	if (Pf3) clReleaseKernel(Pf3);
	if (Pf4) clReleaseKernel(Pf4);
	if (ClScaling) clReleaseKernel(ClScaling);
	if (clBinPartialSort) clReleaseKernel(clBinPartialSort);	


	if (program) clReleaseProgram(program);

	if (context) clReleaseContext(context);
	if (commandQueue) clReleaseCommandQueue(commandQueue);

	if (clMCSRrow) clReleaseMemObject(clMCSRrow);
	if (clMCSRcol) clReleaseMemObject(clMCSRcol);
	if (clM) clReleaseMemObject(clM);	
}


#endif