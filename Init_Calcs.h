#ifndef INIT_CALC_H
#define INIT_CALC_H

#include <math.h>
#include <stdlib.h>
#include <malloc.h>


#include "Globals.h"
#include "MKL_Solver.h"
#include "Do_Time_Step.h"


extern FType ****IFTran;
extern FType ****perm;
extern int Nx, Ny, Nz;
extern FType **fluidProp;
extern FType resTemp;
extern FType *gridDim;
extern FType **jac;
extern FType *Unk;
extern FType *****trans;
extern FType ****Wtran;
extern int PNc, UNc, Nc;
extern FType WOCHeight;
extern int refL;
extern char ***phaseStat;
extern FType *****blockFProps;
extern FType *Ki;
extern char *****transS;	
extern FType cumQo, cumQg, cumQw;
extern FType ***blockH;
extern int incCount;
//extern int jacIndex;
//extern int CSRrowIndex;
extern int ***pJHolder;
extern FType pCoeff;
extern int totalJac, totalRow;
extern int *preConRow;
//extern int *preConIndex;
extern int wellNO;
extern FType *TStepMarker;
extern int TSMSize;
extern FType **bic;
extern FType wellrate[50];
extern FType Dt;
extern int *preConCSRrow;


void Tr_Mu(void);
void EOS_Init(void);
//void Jac_Init(void);
void Edge_Trans(void);
void CPlus_Props(void);
void WaterSat(void);
void SuccFalsh(int Ix, int Iy, int Iz);
//void GCE(void);
//FType SuccFalshFug(int Ix, int Iy, int Iz, FType *Fug);
//char SuccFalshFugDout(int Ix, int Iy, int Iz, FType PressIn, FType *Fug, FType *dFugdP);
FType Pcgo(int Ix, int Iy, int Iz);
FType Pcow(int Ix, int Iy, int Iz);
void AllFlash(void);
void QInitValue(void);
void CalcBlockHeight(void);
//void GCE_NR(void);
//FType GCE_NR_eqns(FType *, FType *, FType, char);
//void GOC(FType *, FType *, FType, char);
void PJac(void);
void CalcPCoeff(void);
void ManageTSMarker(void);



void IFK_Calc(void) {	//calculate kA/h
	register int i, j, k;
	
	for (i=0;i<Nx;i++) 
		for (j=0;j<Ny;j++) 
			for (k=0;k<Nz;k++) {
				if ((Nx!=1) && (i!=(Nx-1))) IFTran[i][j][k][0]=SI_C*2*gridDim[Nx+j]*gridDim[Nx+Ny+k]/(gridDim[i]/perm[i][j][k][0]+gridDim[i+1]/perm[i+1][j][k][0]);
				else IFTran[i][j][k][0]=0;
				if ((Ny!=1) && (j!=(Ny-1))) IFTran[i][j][k][1]=SI_C*2*gridDim[i]*gridDim[Nx+Ny+k]/(gridDim[Nx+j]/perm[i][j][k][1]+gridDim[Nx+j+1]/perm[i][j+1][k][1]);
				else IFTran[i][j][k][1]=0;
				if ((Nz!=1) && (k!=(Nz-1))) IFTran[i][j][k][2]=SI_C*2*gridDim[i]*gridDim[Nx+j]/(gridDim[Nx+Ny+k]/perm[i][j][k][2]+gridDim[Nx+Ny+k+1]/perm[i][j][k+1][2]);
				else IFTran[i][j][k][2]=0;
			}

	for (i=0; i<50; i++) {
		wellrate[i]=((FType) (i+1))/50.0;
	}
}

void Tr_Mu(void) {
	register int i;
	FType gamma;
	for (i=0; i<Nc; i++) {
		//fluidProp[i][TR]=resTemp/fluidProp[i][TCRIT];
		gamma=pow(fluidProp[i][TCRIT],(1.0/6))/(sqrt(fluidProp[i][MW])*pow((fluidProp[i][PCRIT]/101325), (2.0/3)));
		/*if (fluidProp[i][TR]>1.5) {
			fluidProp[i][MU]=(1.778e-4)*pow((4.58*fluidProp[i][TR]-1.67), 0.625)/gamma;
		}
		else {
			fluidProp[i][MU]=(3.4e-4)*pow(fluidProp[i][TR], 0.94)/gamma;
		}*/

		fluidProp[i][MU]=(4.61*pow(fluidProp[i][TR], 0.618)-2.04*exp(-0.449*fluidProp[i][TR])+1.94*exp(-4.058*fluidProp[i][TR])+0.1)*1e-4/gamma;

	}
}

void EOS_Init(void) {
	register int i;
	FType m, ac;

	for (i=0; i<Nc; i++) {
		fluidProp[i][TR]=resTemp/fluidProp[i][TCRIT];
		fluidProp[i][VCRIT]/=1000000.0;
		//fluidProp[i][SE]=1-2.258/pow(fluidProp[i][MW], 0.1823);
		fluidProp[i][COLDIA]=(2.3551-0.087*fluidProp[i][AC])*pow((fluidProp[i][TCRIT]*101325/fluidProp[i][PCRIT]), (1.0/3.0));
		fluidProp[i][EPDIA]=(0.7915+0.1963*fluidProp[i][AC])*fluidProp[i][TCRIT];
	}

	if (SRK) {
		for (i=0; i<Nc; i++) {
			ac=0.42747*RGAS*RGAS*fluidProp[i][TCRIT]*fluidProp[i][TCRIT]/fluidProp[i][PCRIT];
			m=0.48508+1.55171*fluidProp[i][AC]-0.15613*fluidProp[i][AC]*fluidProp[i][AC];
			fluidProp[i][EOS_A]=ac*(1+m*(1-sqrt(fluidProp[i][TR])))*(1+m*(1-sqrt(fluidProp[i][TR])));
			fluidProp[i][EOS_B]=0.08664*RGAS*fluidProp[i][TCRIT]/fluidProp[i][PCRIT];
		}
	}
	else if (PR) {
		for (i=0; i<Nc; i++) {
			ac=0.457235*RGAS*RGAS*fluidProp[i][TCRIT]*fluidProp[i][TCRIT]/fluidProp[i][PCRIT];
			m=0.3796+1.485*fluidProp[i][AC]-0.1644*fluidProp[i][AC]*fluidProp[i][AC]+0.01667*fluidProp[i][AC]*fluidProp[i][AC]*fluidProp[i][AC];
			fluidProp[i][EOS_A]=ac*(1+m*(1-sqrt(fluidProp[i][TR])))*(1+m*(1-sqrt(fluidProp[i][TR])));
			fluidProp[i][EOS_B]=0.077796*RGAS*fluidProp[i][TCRIT]/fluidProp[i][PCRIT];
		}
	}
}

void Edge_Trans(void) {
	register int Ix, Iy, Iz;

	for (Iy=0; Iy<(Ny+1); Iy++)
		for (Iz=0; Iz<(Nz+1); Iz++) {
			trans[0][Iy][Iz][0][0]=0;
			trans[0][Iy][Iz][1][0]=0;

			trans[Nx][Iy][Iz][0][0]=0;
			trans[Nx][Iy][Iz][1][0]=0;

			Wtran[0][Iy][Iz][0]=0;
			Wtran[Nx][Iy][Iz][0]=0;

			dWtran[0][Iy][Iz][0]=0;
			dWtran[Nx][Iy][Iz][0]=0;

			transS[0][Iy][Iz][0][0]=1;
			transS[0][Iy][Iz][1][0]=1;
			transS[0][Iy][Iz][2][0]=1;


			transS[Nx][Iy][Iz][0][0]=0;
			transS[Nx][Iy][Iz][1][0]=0;
			transS[Nx][Iy][Iz][2][0]=0;

			trans[0][Iy][Iz][2][0]=0;
			trans[0][Iy][Iz][3][0]=0;

			trans[Nx][Iy][Iz][2][0]=0;
			trans[Nx][Iy][Iz][3][0]=0;
			//////////////////////////
			trans[0][Iy][Iz][4][0]=0;
			trans[0][Iy][Iz][5][0]=0;

			trans[Nx][Iy][Iz][4][0]=0;
			trans[Nx][Iy][Iz][5][0]=0;

			trans[0][Iy][Iz][6][0]=0;
			trans[0][Iy][Iz][7][0]=0;

			trans[Nx][Iy][Iz][6][0]=0;
			trans[Nx][Iy][Iz][7][0]=0;

			trans[0][Iy][Iz][8][0]=0;
			trans[0][Iy][Iz][9][0]=0;

			trans[Nx][Iy][Iz][8][0]=0;
			trans[Nx][Iy][Iz][9][0]=0;

			trans[0][Iy][Iz][10][0]=0;
			trans[0][Iy][Iz][11][0]=0;

			trans[Nx][Iy][Iz][10][0]=0;
			trans[Nx][Iy][Iz][11][0]=0;
		}

		for (Ix=0; Ix<(Nx+1); Ix++)
			for (Iz=0; Iz<(Nz+1); Iz++) {
				trans[Ix][0][Iz][0][1]=0;			
				trans[Ix][0][Iz][1][1]=0;

				trans[Ix][Ny][Iz][0][1]=0;
				trans[Ix][Ny][Iz][1][1]=0;

				Wtran[Ix][0][Iz][1]=0;
				Wtran[Ix][Ny][Iz][1]=0;

				dWtran[Ix][0][Iz][1]=0;
				dWtran[Ix][Ny][Iz][1]=0;


				transS[Ix][0][Iz][0][1]=1;
				transS[Ix][0][Iz][1][1]=1;
				transS[Ix][0][Iz][2][1]=1;

				transS[Ix][Ny][Iz][0][1]=0;
				transS[Ix][Ny][Iz][1][1]=0;
				transS[Ix][Ny][Iz][2][1]=0;

				trans[Ix][0][Iz][2][1]=0;			
				trans[Ix][0][Iz][3][1]=0;

				trans[Ix][Ny][Iz][2][1]=0;
				trans[Ix][Ny][Iz][3][1]=0;
				/////////////////////////////////
				trans[Ix][0][Iz][4][1]=0;			
				trans[Ix][0][Iz][5][1]=0;

				trans[Ix][Ny][Iz][4][1]=0;
				trans[Ix][Ny][Iz][5][1]=0;

				trans[Ix][0][Iz][6][1]=0;			
				trans[Ix][0][Iz][7][1]=0;

				trans[Ix][Ny][Iz][6][1]=0;
				trans[Ix][Ny][Iz][7][1]=0;

				trans[Ix][0][Iz][8][1]=0;			
				trans[Ix][0][Iz][9][1]=0;

				trans[Ix][Ny][Iz][8][1]=0;
				trans[Ix][Ny][Iz][9][1]=0;

				trans[Ix][0][Iz][10][1]=0;			
				trans[Ix][0][Iz][11][1]=0;

				trans[Ix][Ny][Iz][10][1]=0;
				trans[Ix][Ny][Iz][11][1]=0;
			}

			for (Ix=0; Ix<(Nx+1); Ix++)
				for (Iy=0; Iy<(Ny+1); Iy++) {
					trans[Ix][Iy][0][0][2]=0;
					trans[Ix][Iy][0][1][2]=0;

					trans[Ix][Iy][Nz][0][2]=0;
					trans[Ix][Iy][Nz][1][2]=0;

					Wtran[Ix][Iy][0][2]=0;
					Wtran[Ix][Iy][Nz][2]=0;

					dWtran[Ix][Iy][0][2]=0;
					dWtran[Ix][Iy][Nz][2]=0;


					transS[Ix][Iy][0][0][2]=1;
					transS[Ix][Iy][0][1][2]=1;
					transS[Ix][Iy][0][2][2]=1;

					transS[Ix][Iy][Nz][0][2]=0;
					transS[Ix][Iy][Nz][1][2]=0;
					transS[Ix][Iy][Nz][2][2]=0;

					trans[Ix][Iy][0][2][2]=0;
					trans[Ix][Iy][0][3][2]=0;

					trans[Ix][Iy][Nz][2][2]=0;
					trans[Ix][Iy][Nz][3][2]=0;
					///////////////////////////////
					trans[Ix][Iy][0][4][2]=0;
					trans[Ix][Iy][0][5][2]=0;

					trans[Ix][Iy][Nz][4][2]=0;
					trans[Ix][Iy][Nz][5][2]=0;

					trans[Ix][Iy][0][6][2]=0;
					trans[Ix][Iy][0][7][2]=0;

					trans[Ix][Iy][Nz][6][2]=0;
					trans[Ix][Iy][Nz][7][2]=0;

					trans[Ix][Iy][0][8][2]=0;
					trans[Ix][Iy][0][9][2]=0;

					trans[Ix][Iy][Nz][8][2]=0;
					trans[Ix][Iy][Nz][9][2]=0;

					trans[Ix][Iy][0][10][2]=0;
					trans[Ix][Iy][0][11][2]=0;

					trans[Ix][Iy][Nz][10][2]=0;
					trans[Ix][Iy][Nz][11][2]=0;
				}
				////////////////////////////////////////////////
				for (Iy=1; Iy<(Ny+1); Iy++)
					for (Iz=1; Iz<(Nz+1); Iz++) {
						P[0][Iy][Iz][0]=0;
						P[0][Iy][Iz][1]=0;
						P[0][Iy][Iz][2]=0;

						P[Nx+1][Iy][Iz][0]=0;
						P[Nx+1][Iy][Iz][1]=0;
						P[Nx+1][Iy][Iz][2]=0;

						blockH[0][Iy][Iz]=0;			
						blockH[Nx+1][Iy][Iz]=0;
					}

					for (Ix=1; Ix<(Nx+1); Ix++)
						for (Iz=1; Iz<(Nz+1); Iz++) {
							P[Ix][0][Iz][0]=0;
							P[Ix][0][Iz][1]=0;
							P[Ix][0][Iz][2]=0;

							P[Ix][Ny+1][Iz][0]=0;
							P[Ix][Ny+1][Iz][1]=0;
							P[Ix][Ny+1][Iz][2]=0;


							blockH[Ix][0][Iz]=0;
							blockH[Ix][Ny+1][Iz]=0;
						}

						for (Ix=1; Ix<(Nx+1); Ix++)
							for (Iy=1; Iy<(Ny+1); Iy++) {
								P[Ix][Iy][0][0]=0;
								P[Ix][Iy][0][1]=0;
								P[Ix][Iy][0][2]=0;

								P[Ix][Iy][Nz+1][0]=0;
								P[Ix][Iy][Nz+1][1]=0;
								P[Ix][Iy][Nz+1][2]=0;


								blockH[Ix][Iy][0]=0;
								blockH[Ix][Iy][Nz+1]=0;
							}
}
void CPlus_Props(void) {
	int i;
	FType Tb, API, Tc, Pc, Zc;

	for(i=PNc; i<Nc; i++) {		
		Tb=fluidProp[i][TB]*5/9-459.67;
		API=141.5/fluidProp[i][SG]-131.5;

		//Cavett Correlation
		Tc=768.071+1.7134*Tb-1.0834e-3*Tb*Tb+3.889e-7*Tb*Tb*Tb-8.9213e-3*Tb*API+5.3095e-6*Tb*Tb*API+3.2712e-8*Tb*Tb*API*API;
		fluidProp[i][TCRIT]=(9/5)*Tc;
	
		Pc=2.829+9.412e-4*Tb-3.0475e-6*Tb*Tb+1.5184e-9*Tb*Tb*Tb-2.0876e-5*Tb*API+1.1048e-8*Tb*Tb*API-4.827e-8*Tb*API*API+1.395e-10*Tb*Tb*API*API;
		fluidProp[i][PCRIT]=6897.549353301566*exp(Pc*LN10);

		//Edmister Correlation
		fluidProp[i][AC]=(3.0/7.0)*log(fluidProp[i][PCRIT]/101325)/(LN10*(fluidProp[i][TCRIT]/fluidProp[i][TB]-1))-1;

		//Pitzer correlation
		Zc=0.2901-0.0879*fluidProp[i][AC];
		

		fluidProp[i][VCRIT]=Zc*RGAS*fluidProp[i][TCRIT]/fluidProp[i][PCRIT];
	}
	
}
void WaterSat(void){
	register int i, j, k, n;
	FType watP;

	for (i=0; i<Nx; i++)
		for (j=0; j<Ny; j++) 
			for (k=0; k<Nz; k++) {
				watP=G_ACC*(blockH[i+1][j+1][k+1]-WOCHeight)*watRo;


				for (n=0; n<Nswt; n++) {
					if (watP>swt[n][PC]) break;
				}
				if (n==0) {
					sat[i][j][k][0]=swt[0][SAT];
				}
				else if (n==Nswt) {
					sat[i][j][k][0]=1;
				}
				else {
					sat[i][j][k][0]=(swt[n][SAT]-swt[n-1][SAT])/(swt[n][PC]-swt[n-1][PC])*(watP-swt[n][PC])+swt[n][SAT];
				}
			}
}
void SuccFalsh(int Ix, int Iy, int Iz) {
	register int i, j, k;
	FType L, F, dF, tempD;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG;
	FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType ESum;
	FType tempSQR;
	FType newD;
	FType Lb=10;
	
	for (i=0; i<Nc; i++) {
		Ki[i]=(fluidProp[i][PCRIT]/P[Ix+1][Iy+1][Iz+1][1])*exp(5.37*(1+fluidProp[i][AC])*(1-fluidProp[i][TCRIT]/resTemp));
	}

	k=0;
	do {
		L=0.6;
		j=0;
		do{
			F=0;
			dF=0;
			for (i=0; i<Nc; i++) {
				if (Ki[i]>=0) tempD=(1-Ki[i])/(L+(1-L)*Ki[i]);
				else tempD=1/(L-1);
				F+=comp[Ix][Iy][Iz][i][2]*tempD;
				dF-=comp[Ix][Iy][Iz][i][2]*tempD*tempD;
			}
			tempD=F/dF;
			L-=tempD;
			j++;		
		} while (((tempD*tempD)>NEWTON_TOL) && (j<MAX_FLASH_N_ITR));

		if (L>1) {
			L=1;
			//break;
		}
		else if (L<0) {
			L=0;
			//break;
		}

		ESum=fabs(Lb-L);
		Lb=L;

		for (i=0; i<Nc; i++) {
			if (Ki[i]>=0){
				comp[Ix][Iy][Iz][i][0]=comp[Ix][Iy][Iz][i][2]/(L+Ki[i]*(1-L));
				comp[Ix][Iy][Iz][i][1]=Ki[i]*comp[Ix][Iy][Iz][i][0];
			}
			else{
				comp[Ix][Iy][Iz][i][0]=0;
				comp[Ix][Iy][Iz][i][1]=comp[Ix][Iy][Iz][i][2]/(1-L);
			}
		}

		/////////////////////////////////////////////////////////////

		Al=0;
		Bl=0;
		Ag=0;
		Bg=0;
		for (i=0; i<Nc; i++) {
			Bl+=comp[Ix][Iy][Iz][i][0]*fluidProp[i][EOS_B];
			Bg+=comp[Ix][Iy][Iz][i][1]*fluidProp[i][EOS_B];
			for (j=0; j<Nc; j++) {
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A])*bic[i][j];
				Al+=comp[Ix][Iy][Iz][i][0]*comp[Ix][Iy][Iz][j][0]*tempSQR;
				Ag+=comp[Ix][Iy][Iz][i][1]*comp[Ix][Iy][Iz][j][1]*tempSQR;
			}
		}

		Cal=Al;
		Cbl=Bl;
		Cag=Ag;
		Cbg=Bg;
		Al*=P[Ix+1][Iy+1][Iz+1][1]/(RGAS*RGAS*resTemp*resTemp);
		Bl*=P[Ix+1][Iy+1][Iz+1][1]/(RGAS*resTemp);
		Ag*=P[Ix+1][Iy+1][Iz+1][1]/(RGAS*RGAS*resTemp*resTemp);		//gas-oil
		Bg*=P[Ix+1][Iy+1][Iz+1][1]/(RGAS*resTemp);


		if (SRK) {
			Ul=Bl;
			Wl=0;
			Ug=Bg;
			Wg=0;
			d1=1;
			d2=0;
		}
		else if (PR) {
			Ul=2*Bl;
			Wl=Bl;
			Ug=2*Bg;
			Wg=Bg;
			d1=1+SQRT2;
			d2=1-SQRT2;
		}

		Zl=Solve_Z(-(1+Bl-Ul), Al-Bl*Ul-Ul-Wl*Wl, -(Al*Bl-Bl*Wl*Wl-Wl*Wl), 'l');
		Zg=Solve_Z(-(1+Bg-Ug), Ag-Bg*Ug-Ug-Wg*Wg, -(Ag*Bg-Bg*Wg*Wg-Wg*Wg), 'g');


		for (i=0; i<Nc; i++) {
			tl=0;
			tg=0;
			for (j=0; j<Nc; j++) {	//FFSS
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A])*bic[i][j];
				tl+=comp[Ix][Iy][Iz][j][0]*tempSQR;
				tg+=comp[Ix][Iy][Iz][j][1]*tempSQR;
			}
			t1l=fluidProp[i][EOS_B]/Cbl;
			t1g=fluidProp[i][EOS_B]/Cbg;

			EEl=2*tl/Cal-t1l;
			EEg=2*tg/Cag-t1g;

			DDl=Al/(Bl*(d1-d2));
			DDg=Ag/(Bg*(d1-d2));

			FFl=log((Zl+d2*Bl)/(Zl+d1*Bl));
			FFg=log((Zg+d2*Bg)/(Zg+d1*Bg));

			phiL=exp(t1l*(Zl-1)-log(Zl-Bl)+DDl*EEl*FFl);
			phiG=exp(t1g*(Zg-1)-log(Zg-Bg)+DDg*EEg*FFg);

			fL=comp[Ix][Iy][Iz][i][0]*phiL;		//eliminated P
			fG=comp[Ix][Iy][Iz][i][1]*phiG;	//gas-oil

			if (comp[Ix][Iy][Iz][i][2]) {
				if (fG) {
					Ki[i]*=fL/fG;
				}
				else Ki[i]=-1;
				if (fL!=fL) Ki[i]=-1;
				if (fG!=fG) Ki[i]=1;
			}
		}
		k++;
	} while ((ESum>FLASH_TOL) && (k<500));
	//////////////////////////////////////////////////////////////

	if (fabs(L)<DZERO) {
		phaseStat[Ix][Iy][Iz]=-1;
		//L=0;
	}
	else if (fabs(1-L)<DZERO){
		phaseStat[Ix][Iy][Iz]=1;
		//L=1;
	}
	else {
		phaseStat[Ix][Iy][Iz]=0;
	}


	newD=Zl*L/(Zl*L+Zg*(1-L));

	sat[Ix][Iy][Iz][1]=newD*(1-sat[Ix][Iy][Iz][0]);
	sat[Ix][Iy][Iz][2]=(1-newD)*(1-sat[Ix][Iy][Iz][0]);

	blockFProps[Ix][Iy][Iz][RO][0]=P[Ix+1][Iy+1][Iz+1][1]/(Zl*RGAS*resTemp);
	blockFProps[Ix][Iy][Iz][RO][1]=P[Ix+1][Iy+1][Iz+1][1]/(Zg*RGAS*resTemp);	//gas-oil
	
}

/*void GCE(void) {
	register int Ix, Iy, Iz, n;
	FType height, refHeight, Ppress, SYi, _SYi, dQdP, _Tol;
	FType *Fug0, *Yi, *Zi, *Fug1, *dFug, *r;
	char PhStat;


	if ((Fug0=(FType *) malloc(Nc*sizeof(FType)))==NULL) {
		free(Fug0);
		TerM("Can not allocate memory for fugacity matrix");
	}
	if ((Yi=(FType *) malloc(Nc*sizeof(FType)))==NULL) {
		free(Fug0);
		free(Yi);
		TerM("Can not allocate memory for Yi matrix");
	}
	if ((Zi=(FType *) malloc(Nc*sizeof(FType)))==NULL) {
		free(Fug0);
		free(Yi);
		free(Zi);
		TerM("Can not allocate memory for Zi matrix");
	}
	if ((Fug1=(FType *) malloc(Nc*sizeof(FType)))==NULL) {
		free(Fug1);
		free(Fug0);
		free(Yi);
		free(Zi);
		TerM("Can not allocate memory for fugacity1 matrix");
	}
	if ((dFug=(FType *) malloc(Nc*sizeof(FType)))==NULL) {
		free(Fug1);
		free(Fug0);
		free(Yi);
		free(Zi);
		free(dFug);
		TerM("Can not allocate memory for dFug matrix");
	}
	if ((r=(FType *) malloc(Nc*sizeof(FType)))==NULL) {
		free(Fug1);
		free(Fug0);
		free(Yi);
		free(Zi);
		free(dFug);
		free(r);
		TerM("Can not allocate memory for r matrix");
	}



	WaterSat();
	SuccFalshFug(0, 0, refL, Fug0);

	for (Iz=(Nz-1); Iz>=0; Iz--) {
		if (Iz==(Nz-1)) {
			height=0.5*gridDim[Nx+Ny+Nz-1];
			refHeight=height;
		}
		else {
			height+=0.5*(gridDim[Nx+Ny+Iz]+gridDim[Nx+Ny+Iz+1]);
			if (Iz>=refL) refHeight+=0.5*(gridDim[Nx+Ny+Iz]+gridDim[Nx+Ny+Iz+1]);
		}
		
		
		Ppress=P[1][1][refL+1][1];
		
		do {
			PhStat=SuccFalshFugDout(0, 0, refL, Ppress, Fug1, dFug);

			SYi=0;
			for (n=0; n<Nc; n++) {
				Yi[n]=Zi[n]*Fug0[n]*exp(-fluidProp[n][MW]*G_ACC*(height-refHeight)/(RGAS*resTemp))/Fug1[n];
				SYi+=Yi[n];
			}

			_SYi=0;
			dQdP=0;
			for (n=0; n<Nc; n++) {
				r[n]=Fug0[n]*exp(-fluidProp[n][MW]*G_ACC*(height-refHeight)/(RGAS*resTemp))/(Fug1[n]*SYi);
				Yi[n]*=r[n];		//Without Lambda, No Acceleration
				_SYi+=Yi[n];
				dQdP+=Yi[n]*r[n]*dFug[n]/Fug1[n];
			}

			_Tol=0;
			for (n=0; n<Nc; n++) {
				Zi[n]=Yi[n]/_SYi;
				_Tol+=log(r[n])/log(Yi[n]/Zi[n]);
			}
			_SYi=1-_SYi;
			Ppress-=_SYi/dQdP;

		} while ((_SYi>1e-13) || (_Tol>1e-8));

		if (PhStat==-1) {
			//P[0][0][Iz][2]=Ppress;
			P[1][1][Iz+1][1]=Ppress-Pcgo(0, 0, refL);
			//P[0][0][Iz][0]=P[0][0][Iz][1]-Pcow(0, 0, refL);
		}
		else if (PhStat==1) {
			P[1][1][Iz+1][1]=Ppress;
			//P[0][0][Iz][2]=Ppress+Pcgo(0, 0, refL);
			//P[0][0][Iz][0]=P[0][0][Iz][1]-Pcow(0, 0, refL);
		}
		else {
			P[1][1][Iz+1][1]=Ppress;
			//P[0][0][Iz][2]=Ppress+Pcgo(0, 0, refL);
			//P[0][0][Iz][0]=P[0][0][Iz][1]-Pcow(0, 0, refL);
		}
		

		for(Ix=1; Ix<Nx; Ix++)
			for(Iy=1; Iy<Ny; Iy++) {
				P[Ix+1][Iy+1][Iz+1][0]=P[1][1][Iz+1][0];
				//P[Ix][Iy][Iz][1]=P[0][0][Iz][1];
				//P[Ix][Iz][Iz][2]=P[0][0][Iz][2];

				sat[Ix][Iy][Iz][1]=sat[0][0][Iz][1];

				phaseStat[Ix][Iy][Iz]=phaseStat[0][0][Iz];

				for (n=0; n<Nc; n++) {
					comp[Ix][Iy][Iz][n][0]=comp[0][0][Iz][n][0];
					comp[Ix][Iy][Iz][n][1]=comp[0][0][Iz][n][1];
				}
			}
	}

	free(Fug0);
	free(Fug1);
	free(Yi);
	free(Zi);
	free(dFug);
	free(r);
}

FType SuccFalshFug(int Ix, int Iy, int Iz, FType *Fug) {
	register int i, j;
	FType L, F, dF, tempD;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG;
	FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType ESum;
	FType tempSQR;
	FType newD;

	for (i=0; i<Nc; i++) {
		Ki[i]=(fluidProp[i][PCRIT]/P[Ix+1][Iy+1][Iz+1][1])*exp(5.37*(1+fluidProp[i][AC])*(1-fluidProp[i][TCRIT]/resTemp));
	}

	do {
		L=0.5;
		do{
			F=0;
			dF=0;
			for (i=0; i<Nc; i++) {
				tempD=(1-Ki[i])/(L+(1-L)*Ki[i]);
				F+=comp[Ix][Iy][Iz][i][2]*tempD;
				dF-=comp[Ix][Iy][Iz][i][2]*tempD*tempD;
			}
			tempD=F/dF;
			L-=tempD;
			if (L>NO_ROOT_L) {
				L=1;
				break;
			}
			else if (L<(-NO_ROOT_L)) {
				L=0;
				break;
			}
		} while ((tempD*tempD)>NEWTON_TOL);

		for (i=0; i<Nc; i++) {
			comp[Ix][Iy][Iz][i][0]=comp[Ix][Iy][Iz][i][2]/(L+Ki[i]*(1-L));
			comp[Ix][Iy][Iz][i][1]=Ki[i]*comp[Ix][Iy][Iz][i][0];
		}

		/////////////////////////////////////////////////////////////

		Al=0;
		Bl=0;
		Ag=0;
		Bg=0;
		for (i=0; i<Nc; i++) {
			Bl+=comp[Ix][Iy][Iz][i][0]*fluidProp[i][EOS_B];
			Bg+=comp[Ix][Iy][Iz][i][1]*fluidProp[i][EOS_B];
			for (j=0; j<Nc; j++) {
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A]);
				Al+=comp[Ix][Iy][Iz][i][0]*comp[Ix][Iy][Iz][j][0]*tempSQR;
				Ag+=comp[Ix][Iy][Iz][i][1]*comp[Ix][Iy][Iz][j][1]*tempSQR;
			}
		}

		Cal=Al;
		Cbl=Bl;
		Cag=Ag;
		Cbg=Bg;
		Al*=P[Ix+1][Iy+1][Iz+1][1]/(RGAS*RGAS*resTemp*resTemp);
		Bl*=P[Ix+1][Iy+1][Iz+1][1]/(RGAS*resTemp);
		Ag*=P[Ix+1][Iy+1][Iz+1][1]/(RGAS*RGAS*resTemp*resTemp);	//Oil press note
		Bg*=P[Ix+1][Iy+1][Iz+1][1]/(RGAS*resTemp);		//oil press note


		if (SRK) {
			Ul=Bl;
			Wl=0;
			Ug=Bg;
			Wg=0;
			d1=1;
			d2=0;
		}
		else if (PR) {
			Ul=2*Bl;
			Wl=Bl;
			Ug=2*Bg;
			Wg=Bg;
			d1=1+SQRT2;
			d2=1-SQRT2;
		}

		Zl=Solve_Z(-(1+Bl-Ul), Al-Bl*Ul-Ul-Wl*Wl, -(Al*Bl-Bl*Wl*Wl-Wl*Wl), 'l');
		Zg=Solve_Z(-(1+Bg-Ug), Ag-Bg*Ug-Ug-Wg*Wg, -(Ag*Bg-Bg*Wg*Wg-Wg*Wg), 'g');



		ESum=0;
		for (i=0; i<Nc; i++) {
			tl=0;
			tg=0;
			for (j=0; j<Nc; j++) {	//FFSS
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A]);
				tl+=comp[Ix][Iy][Iz][j][0]*tempSQR;
				tg+=comp[Ix][Iy][Iz][j][1]*tempSQR;
			}
			t1l=fluidProp[i][EOS_B]/Cbl;
			t1g=fluidProp[i][EOS_B]/Cbg;

			EEl=2*tl/Cal-t1l;
			EEg=2*tg/Cag-t1g;

			DDl=Al/(Bl*(d1-d2));
			DDg=Ag/(Bg*(d1-d2));

			FFl=log((Zl+d2*Bl)/(Zl+d1*Bl));
			FFg=log((Zg+d2*Bg)/(Zg+d1*Bg));

			phiL=exp(t1l*(Zl-1)-log(Zl-Bl)+DDl*EEl*FFl);
			phiG=exp(t1g*(Zg-1)-log(Zg-Bg)+DDg*EEg*FFg);

			fL=comp[Ix][Iy][Iz][i][0]*P[Ix+1][Iy+1][Iz+1][1]*phiL;
			fG=comp[Ix][Iy][Iz][i][1]*P[Ix+1][Iy+1][Iz+1][1]*phiG;

			if (L<DZERO) {
				Fug[i]=fG;
				phaseStat[Ix][Iy][Iz]=-1;
			}
			else if ((1-L)<DZERO){
				Fug[i]=fL;
				phaseStat[Ix][Iy][Iz]=1;
			}
			else {
				Fug[i]=fL;
				phaseStat[Ix][Iy][Iz]=0;
			}

			Ki[i]*=fL/fG;

			ESum+=(1-fL/fG)*(1-fL/fG);
		}
	} while (ESum>FLASH_TOL);
	//////////////////////////////////////////////////////////////

	newD=Zl*L/(Zl*L+Zg*(1-L));

	sat[Ix][Iy][Iz][1]=newD*(1-sat[Ix][Iy][Iz][0]);
	sat[Ix][Iy][Iz][2]=(1-newD)*(1-sat[Ix][Iy][Iz][0]);

	blockFProps[Ix][Iy][Iz][RO][0]=P[Ix+1][Iy+1][Iz+1][1]/(Zl*RGAS*resTemp);
	blockFProps[Ix][Iy][Iz][RO][1]=P[Ix+1][Iy+1][Iz+1][1]/(Zg*RGAS*resTemp);

	return phaseStat[Ix][Iy][Iz]=0;	
}

char SuccFalshFugDout(int Ix, int Iy, int Iz, FType PressIn, FType *Fug, FType *dFugdP) {
	register int i, j;
	char rCh;
	FType L, F, dF, tempD;
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG;
	FType DDl, DDg, EEl, EEg, FFl, FFg;
	FType fL, fG;
	FType ESum;
	FType dAdpl, dBdpl, dAdpg, dBdpg,  dZdpl, dZdpg;
	FType dAAdpl, dAAdpg, dBBdpl, dBBdpg, dDDdpl, dDDdpg, dFFdpl, dFFdpg, dPhidpl, dPhidpg;
	FType dUdpg, dWdpg, dUdpl, dWdpl;
	FType tempSQR;
	FType newD;

	for (i=0; i<Nc; i++) {
		Ki[i]=(fluidProp[i][PCRIT]/PressIn)*exp(5.37*(1+fluidProp[i][AC])*(1-fluidProp[i][TCRIT]/resTemp));
	}

	do {
		L=0.5;
		do{
			F=0;
			dF=0;
			for (i=0; i<Nc; i++) {
				tempD=(1-Ki[i])/(L+(1-L)*Ki[i]);
				F+=comp[Ix][Iy][Iz][i][2]*tempD;
				dF-=comp[Ix][Iy][Iz][i][2]*tempD*tempD;
			}
			tempD=F/dF;
			L-=tempD;
			if (L>NO_ROOT_L) {
				L=1;
				break;
			}
			else if (L<(-NO_ROOT_L)) {
				L=0;
				break;
			}
		} while ((tempD*tempD)>NEWTON_TOL);

		for (i=0; i<Nc; i++) {
			comp[Ix][Iy][Iz][i][0]=comp[Ix][Iy][Iz][i][2]/(L+Ki[i]*(1-L));
			comp[Ix][Iy][Iz][i][1]=Ki[i]*comp[Ix][Iy][Iz][i][2];
		}

		/////////////////////////////////////////////////////////////

		Al=0;
		Bl=0;
		Ag=0;
		Bg=0;
		for (i=0; i<Nc; i++) {
			Bl+=comp[Ix][Iy][Iz][i][0]*fluidProp[i][EOS_B];
			Bg+=comp[Ix][Iy][Iz][i][1]*fluidProp[i][EOS_B];
			for (j=0; j<Nc; j++) {
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A]);
				Al+=comp[Ix][Iy][Iz][i][0]*comp[Ix][Iy][Iz][j][0]*tempSQR;
				Ag+=comp[Ix][Iy][Iz][i][1]*comp[Ix][Iy][Iz][j][1]*tempSQR;
			}
		}

		Cal=Al;
		Cbl=Bl;
		Cag=Ag;
		Cbg=Bg;
		Al*=PressIn/(RGAS*RGAS*resTemp*resTemp);
		Bl*=PressIn/(RGAS*resTemp);
		Ag*=PressIn/(RGAS*RGAS*resTemp*resTemp);
		Bg*=PressIn/(RGAS*resTemp);

		dAdpl=Cal/(RGAS*RGAS*resTemp*resTemp);
		dAdpg=Cag/(RGAS*RGAS*resTemp*resTemp);
		dBdpl=Cbl/(RGAS*resTemp);
		dBdpg=Cbg/(RGAS*resTemp);



		if (SRK) {
			Ul=Bl;
			Wl=0;
			Ug=Bg;
			Wg=0;
			d1=1;
			d2=0;

			dUdpl=dBdpl;
			dWdpl=0;
			dUdpg=dBdpg;
			dWdpg=0;
		}
		else if (PR) {
			Ul=2*Bl;
			Wl=Bl;
			Ug=2*Bg;
			Wg=Bg;
			d1=1+SQRT2;
			d2=1-SQRT2;

			dUdpg=2*dBdpg;
			dWdpg=dBdpg;
			dUdpl=2*dBdpl;
			dWdpl=dBdpl;
		}

		Zl=Solve_Z(-(1+Bl-Ul), Al-Bl*Ul-Ul-Wl*Wl, -(Al*Bl-Bl*Wl*Wl-Wl*Wl), 'l');
		Zg=Solve_Z(-(1+Bg-Ug), Ag-Bg*Ug-Ug-Wg*Wg, -(Ag*Bg-Bg*Wg*Wg-Wg*Wg), 'g');


		dZdpg=(dAdpg*(Zg-Bg)+dBdpg*(-Zg*Zg-Ug*Zg-Ag+Wg*Wg)+dUdpg*(Zg*Zg-Bg*Zg-Zg)+dWdpg*(-2*Wg*Zg+2*Bg*Wg+2*Wg))/(3*Zg*Zg-2*Zg*(1+Bg-Ug)+Ag+Bg*Ug-Ug-Wg*Wg);
		dZdpl=(dAdpl*(Zl-Bl)+dBdpl*(-Zl*Zl-Ul*Zl-Al+Wl*Wl)+dUdpl*(Zl*Zl-Bl*Zl-Zl)+dWdpl*(-2*Wl*Zl+2*Bl*Wl+2*Wl))/(3*Zl*Zl-2*Zl*(1+Bl-Ul)+Al+Bl*Ul-Ul-Wl*Wl);



		ESum=0;
		for (i=0; i<Nc; i++) {
			tl=0;
			tg=0;
			for (j=0; j<Nc; j++) {	//FFSS
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A]);
				tl+=comp[Ix][Iy][Iz][j][0]*tempSQR;
				tg+=comp[Ix][Iy][Iz][j][1]*tempSQR;
			}
			t1l=fluidProp[i][EOS_B]/Cbl;
			t1g=fluidProp[i][EOS_B]/Cbg;

			EEl=2*tl/Cal-t1l;
			EEg=2*tg/Cag-t1g;

			DDl=Al/(Bl*(d1-d2));
			DDg=Ag/(Bg*(d1-d2));

			FFl=log((Zl+d2*Bl)/(Zl+d1*Bl));
			FFg=log((Zg+d2*Bg)/(Zg+d1*Bg));

			phiL=exp(t1l*(Zl-1)-log(Zl-Bl)+DDl*EEl*FFl);
			phiG=exp(t1g*(Zg-1)-log(Zg-Bg)+DDg*EEg*FFg);

			fL=comp[Ix][Iy][Iz][i][0]*PressIn*phiL;
			fG=comp[Ix][Iy][Iz][i][1]*PressIn*phiG;



			dAAdpl=fluidProp[i][EOS_B]*dZdpl/Cbl;
			dAAdpg=fluidProp[i][EOS_B]*dZdpg/Cbg;

			dBBdpl=(dBdpl-dZdpl)/(Zl-Bl);
			dBBdpg=(dBdpg-dZdpg)/(Zg-Bg);

			dDDdpl=(Bl*dAdpl-Al*dBdpl)/(Bl*Bl*(d1-d2));
			dDDdpg=(Bg*dAdpg-Ag*dBdpg)/(Bg*Bg*(d1-d2));

			dFFdpl=(dZdpl+d2*dBdpl)/(Zl+d2*Bl)-(dZdpl+d1*dBdpl)/(Zl+d1*Bl);
			dFFdpg=(dZdpg+d2*dBdpg)/(Zg+d2*Bg)-(dZdpg+d1*dBdpg)/(Zg+d1*Bg);

			dPhidpl=phiL*(dAAdpl+dBBdpl+EEl*(dDDdpl*FFl+DDl*dFFdpl));
			dPhidpg=phiG*(dAAdpg+dBBdpg+EEg*(dDDdpg*FFg+DDg*dFFdpg));




			if (L<DZERO) {
				Fug[i]=fG;
				dFugdP[i]=comp[Ix][Iy][Iz][i][1]*(phiG+PressIn*dPhidpg);
				rCh=-1;
				phaseStat[Ix][Iy][Iz]=-1;
			}
			else if ((1-L)<DZERO) {
				Fug[i]=fL;
				dFugdP[i]=comp[Ix][Iy][Iz][i][0]*(phiL+PressIn*dPhidpl);
				rCh=1;
				phaseStat[Ix][Iy][Iz]=1;
			}
			else {
				Fug[i]=fL;
				dFugdP[i]=comp[Ix][Iy][Iz][i][0]*(phiL+PressIn*dPhidpl);
				rCh=0;
				phaseStat[Ix][Iy][Iz]=0;
			}

			Ki[i]*=fL/fG;

			ESum+=(1-fL/fG)*(1-fL/fG);
		}
	} while (ESum>FLASH_TOL);
	//////////////////////////////////////////////////////////////

	newD=Zl*L/(Zl*L+Zg*(1-L));

	sat[Ix][Iy][Iz][1]=newD*(1-sat[Ix][Iy][Iz][0]);
	sat[Ix][Iy][Iz][2]=(1-newD)*(1-sat[Ix][Iy][Iz][0]);

	blockFProps[Ix][Iy][Iz][RO][0]=PressIn/(Zl*RGAS*resTemp);
	blockFProps[Ix][Iy][Iz][RO][1]=PressIn/(Zg*RGAS*resTemp);
	
	return rCh;
}*/


FType Pcgo(int Ix, int Iy, int Iz) {
	register int n;
	FType rd;

	for (n=0; n<Nsgt; n++) {
		if (sat[Ix][Iy][Iz][2]<sgt[n][SAT]) break;
	}
	if (n==0) {
		rd=0;
	}
	else if (n==Nsgt) {
		rd=sgt[Nsgt-1][PC];
	}
	else {
		rd=(sgt[n][PC]-sgt[n-1][PC])/(sgt[n][SAT]-sgt[n-1][SAT])*(sat[Ix][Iy][Iz][2]-sgt[n][SAT])+sgt[n][PC];
	}
	
	return rd;
}

FType Pcow(int Ix, int Iy, int Iz) {
	register int n;
	FType rd;

	for (n=0; n<Nswt; n++) {
		if (sat[Ix][Iy][Iz][0]<swt[n][SAT]) break;
	}
	if (n==0) {
		rd=swt[0][PC];
	}
	else if (n==Nswt) {
		rd=0;
	}
	else {
		rd=(swt[n][PC]-swt[n-1][PC])/(swt[n][SAT]-swt[n-1][SAT])*(sat[Ix][Iy][Iz][0]-swt[n][SAT])+swt[n][PC];
	}

	return rd;
}
void AllFlash(void) {
	register int Ix, Iy, Iz;
	for (Ix=0; Ix<Nx; Ix++)
		for (Iy=0; Iy<Ny; Iy++)
			for (Iz=0; Iz<Nz; Iz++) {
				SuccFalsh(Ix, Iy, Iz);
				blockFProps[Ix][Iy][Iz][BLOCK_PC][1]=Pcgo(Ix, Iy, Iz);
				blockFProps[Ix][Iy][Iz][BLOCK_PC][0]=Pcow(Ix, Iy, Iz);
			}

}

void QInitValue(void) {
	cumQo=0;
	cumQg=0;
	cumQw=0;

	incCount=0;
}

void CalcBlockHeight(void) {
	register int i, j, k;
	FType height;

	for (i=0; i<Nx; i++)
		for (j=0; j<Ny; j++) 
			for (k=(Nz-1); k>=0; k--) {
				if (k==(Nz-1)) {
					height=0.5*gridDim[Nx+Ny+Nz-1];
				}
				else {
					height+=0.5*(gridDim[Nx+Ny+k]+gridDim[Nx+Ny+k+1]);
				}

				blockH[i+1][j+1][k+1]=height;				
			}				
}

/*void GCE_NR(void) {
	register int Ix, Iy, Iz, i;
	FType *Fug0, *Ci0, *iComp;
	int L2Phase;
	char ChS, ChChS;
	FType Ip, PTol;

	if ((Fug0=(FType *) malloc(Nc*sizeof(FType)))==NULL) {
		free(Fug0);
		TerM("Can not allocate memory for fugacity matrix");
	}
	if ((Ci0=(FType *) malloc(Nc*sizeof(FType)))==NULL) {
		free(Fug0);
		free(Ci0);
		TerM("Can not allocate memory for Ci0 matrix");
	}
	if ((iComp=(FType *) malloc(Nc*sizeof(FType)))==NULL) {
		free(iComp);
		free(Fug0);
		free(Ci0);
		TerM("Can not allocate memory for iComp matrix");
	}
	
	WaterSat();
	SuccFalshFug(0, 0, refL, Fug0);

	for (i=0; i<Nc; i++) {		
		iComp[i]=comp[0][0][refL][i][0];
	}
	ChS=0;
	GOC(iComp, Fug0, P[1][1][refL+1][1], ChS);

	/*if (!phaseStat[0][0][refL]) {
		for (Iz=0; Iz<Nz; Iz++) {
			if (Iz<refL) ChS=2;		//gas
			else if (Iz>refL) ChS=1;	//liquid
			else continue;

			for (i=0; i<Nc; i++) {
				Ci0[i]=Fug0[i]*exp(fluidProp[i][MW]*G_ACC*(blockH[1][1][refL+1]-blockH[1][1][Iz+1])/(RGAS*resTemp));
				iComp[i]=comp[0][0][refL][i][2];
			}			
			//P[1][1][Iz+1][ChS]=GCE_NR_eqns(iComp, Ci0, P[1][1][refL], ChS);
			if (ChS==2) P[1][1][Iz+1][1]=P[1][1][Iz+1][2]-Pcgo(0, 0, Iz);
		}
	}

	else if (phaseStat[0][0][refL]==1) {
		for (Iz=(Nz-1); Iz>0; Iz--) {
			if (Iz<refL) ChS=1;		//gas
			else if (Iz>refL) ChS=1;	//liquid
			else continue;

			for (i=0; i<Nc; i++) {
				Ci0[i]=Fug0[i]*exp(fluidProp[i][MW]*G_ACC*(blockH[1][1][refL+1]-blockH[1][1][Iz+1])/(RGAS*resTemp));
				iComp[i]=comp[0][0][refL][i][2];
			}			
			//P[1][1][Iz+1][ChS]=GCE_NR_eqns(iComp, Ci0, P[1][1][refL+1], ChS);
			if (Iz<refL){
				SuccFalsh(0, 0, Iz);
				if (phaseStat[0][0][Iz]==0) {
					do {
						for (i=0; i<Nc; i++) {
							iComp[i]=comp[0][0][Iz][i][0];
						}
						PTol=P[1][1][Iz+1][1];
						//P[1][1][Iz+1][1]=GCE_NR_eqns(iComp, Ci0, P[1][1][Iz+1], 1);
						PTol=fabs(P[1][1][Iz+1][1]-PTol);
					} while (PTol>PRESSTOL);
				}

				
			}
		}
	}*/

	/* ezafe ast


	free(Fug0);
	free(Ci0);
	free(iComp);
}

FType GCE_NR_eqns(FType *initComp, FType *Ci, FType Ipress, char PhStat) {
	FType Zg, Ag, Bg;
	FType Ug, Wg;
	FType tg, d1, d2, t1g;
	FType Cag, Cbg;
	FType phiG, dAdxg, dBdxg, dZdxg;
	FType dAAdxg, dBBdxg, DDg, EEg, FFg;
	FType dFFdxg, dEEdxg, dDDdxg, dPhidxg;
	FType dAdpg, dBdpg,  dZdpg;
	FType dAAdpg, dBBdpg, dDDdpg, dFFdpg, dPhidpg;
	FType dUdpg, dWdpg;
	FType dUdxg, dWdxg;
	register int i, j, k, n;

	FType SXm;
	FType XmTol;
	FType tempSQR;

	char IPhStat;

	if (PhStat==2) IPhStat='g';
	else if (PhStat==1) IPhStat='l';


	Xm[Nc]=Ipress;
	for (i=0; i<Nc; i++){
		Xm[i]=initComp[i];
	}

	do {
		SXm=0;
		for (i=0; i<Nc; i++){
			if (Xm[i]<0) Xm[i]=0;
			else if (Xm[i]>1) Xm[i]=1;
			SXm+=Xm[i];
			satJac[Nc][i]=1;
		}

		if (Xm[Nc]<0) Xm[Nc]=PZERO;

		satAns[Nc]=-SXm+1;
		satJac[Nc][Nc]=0;


		Ag=0;
		Bg=0;
		for (i=0; i<Nc; i++) {
			Bg+=Xm[i]*fluidProp[i][EOS_B];
			for (j=0; j<Nc; j++) {
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A]);
				Ag+=Xm[i]*Xm[j]*tempSQR;
			}
		}

		Cag=Ag;
		Cbg=Bg;
		Ag*=Xm[Nc]/(RGAS*RGAS*resTemp*resTemp);
		Bg*=Xm[Nc]/(RGAS*resTemp);


		dAdpg=Cag/(RGAS*RGAS*resTemp*resTemp);
		dBdpg=Cbg/(RGAS*resTemp);
		
		if (SRK) {
			Ug=Bg;
			Wg=0;

			dUdpg=dBdpg;
			dWdpg=0;
			
			d1=1;
			d2=0;
		}
		else if (PR) {
			Ug=2*Bg;
			Wg=Bg;

			dUdpg=2*dBdpg;
			dWdpg=dBdpg;
			
			d1=1+SQRT2;
			d2=1-SQRT2;
		}

		Zg=Solve_Z(-(1+Bg-Ug), Ag-Bg*Ug-Ug-Wg*Wg, -(Ag*Bg-Bg*Wg*Wg-Wg*Wg), IPhStat);


		dZdpg=-(dAdpg*(Zg-Bg)+dBdpg*(-Zg*Zg-Ug*Zg-Ag+Wg*Wg)+dUdpg*(Zg*Zg-Bg*Zg-Zg)+dWdpg*(-2*Wg*Zg+2*Bg*Wg+2*Wg))/(3*Zg*Zg-2*Zg*(1+Bg-Ug)+Ag-Bg*Ug-Ug-Wg*Wg);		//Modified
		
		for (i=0; i<Nc; i++) {
			tg=0;
			for (j=0; j<Nc; j++) {	//FFSS
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A]);
				tg+=Xm[j]*tempSQR;
			}
			t1g=fluidProp[i][EOS_B]/Cbg;

			EEg=2*tg/Cag-t1g;

			DDg=Ag/(Bg*(d1-d2));

			FFg=log((Zg+d2*Bg)/(Zg+d1*Bg));

			phiG=exp(t1g*(Zg-1)-log(Zg-Bg)+DDg*EEg*FFg);

			satAns[i]=Ci[i]-Xm[i]*Xm[Nc]*phiG;

			dAAdpg=fluidProp[i][EOS_B]*dZdpg/Cbg;
			
			dBBdpg=(dBdpg-dZdpg)/(Zg-Bg);
			
			dDDdpg=(Bg*dAdpg-Ag*dBdpg)/(Bg*Bg*(d1-d2));
			
			dFFdpg=(dZdpg+d2*dBdpg)/(Zg+d2*Bg)-(dZdpg+d1*dBdpg)/(Zg+d1*Bg);
			
			dPhidpg=phiG*(dAAdpg+dBBdpg+EEg*(dDDdpg*FFg+DDg*dFFdpg));
			
			satJac[i][Nc]=Xm[i]*(phiG+Xm[Nc]*dPhidpg);		//Modified
			
			//i: Row of Jac
			//k: Col of Jac
			for (k=0; k<Nc; k++) {
				dAdxg=0;
				for (n=0; n<Nc; n++) {		//SSFF
					dAdxg+=Xm[n]*sqrt(fluidProp[n][EOS_A]*fluidProp[k][EOS_A]);
					//This loop can be reduced
				}
				dAdxg*=2*Xm[Nc]/(RGAS*RGAS*resTemp*resTemp);
				dBdxg=fluidProp[k][EOS_B]*Xm[Nc]/(RGAS*resTemp);

				if (SRK) {
					dUdxg=dBdxg;
					dWdxg=0;
				}
				else if (PR) {
					dUdxg=2*dBdxg;
					dWdxg=dBdxg;
				}


				dZdxg=-(dAdxg*(Zg-Bg)+dBdxg*(-Zg*Zg-Ug*Zg-Ag+Wg*Wg)+dUdxg*(Zg*Zg-Bg*Zg-Zg)+dWdxg*(-2*Wg*Zg+2*Bg*Wg+2*Wg))/(3*Zg*Zg-2*Zg*(1+Bg-Ug)+Ag-Bg*Ug-Ug-Wg*Wg);		//Modified
				
				dAAdxg=-fluidProp[k][EOS_B]*fluidProp[i][EOS_B]*(Zg-1)/(Cbg*Cbg)+fluidProp[i][EOS_B]*dZdxg/Cbg;

				dBBdxg=(dBdxg-dZdxg)/(Zg-Bg);

				dDDdxg=(dAdxg*Bg-dBdxg*Ag)/(Bg*Bg*(d1-d2));

				dFFdxg=(dZdxg+d2*dBdxg)/(Zg+d2*Bg)-(dZdxg+d1*dBdxg)/(Zg+d1*Bg);

				dEEdxg=2*(sqrt(fluidProp[i][EOS_A]*fluidProp[k][EOS_A])/Cag-dAdxg*RGAS*RGAS*resTemp*resTemp*tg/(Xm[Nc]*Cag*Cag))+fluidProp[i][EOS_B]*fluidProp[k][EOS_B]/(Cbg*Cbg);		//modified
				
				dPhidxg=phiG*(dAAdxg+dBBdxg+dDDdxg*EEg*FFg+DDg*dEEdxg*FFg+DDg*EEg*dFFdxg);

				if (i==k) {
					satJac[i][k]=Xm[Nc]*(phiG+Xm[i]*dPhidxg);
				}
				else {
					satJac[i][k]=Xm[Nc]*Xm[i]*dPhidxg;
				}


			}


		}



		MKLS(satJac, satAns, Xms, Nc+1);

		XmTol=0;
		for (i=0; i<(Nc+1); i++) {
			Xm[i]+= Xms[i];
			XmTol+=fabs(Xms[i]/Xm[i]);
		}
		XmTol/=(Nc+1);

	} while ((XmTol>XMTOL) && (fabs(Xms[Nc])>PRESSTOL));

	return Xm[Nc];
}

void GOC(FType *initComp, FType *Fug0, FType Ipress, char PhStat) {
	FType Ul, Wl, Al, Bl, Zl, Zg, Ag, Bg;
	FType Ug, Wg;
	FType tl, tg, d1, d2, t1l, t1g;
	FType Cal, Cbl, Cag, Cbg;
	FType phiL, phiG, dAdxg, dBdxg, dZdxg, dZdxl;
	FType dAAdxg, dBBdxg, DDl, DDg, EEl, EEg, FFl, FFg, dAAdxl, dBBdxl;
	FType dFFdxg, dEEdxg, dDDdxg, dPhidxg, dFFdxl, dEEdxl, dDDdxl, dPhidxl;
	FType dAdpg, dBdpg,  dZdpg;
	FType dAAdpg, dBBdpg, dDDdpg, dFFdpg, dPhidpg;
	FType dUdpg, dWdpg;
	FType dUdxg, dWdxg, dUdxl, dWdxl;
	register int i, j, k, n;

	FType NRKi;
	FType SXm, SYm;
	FType XmTol;

	FType dAdpl, dBdpl, dZdpl;
	FType dAAdpl, dBBdpl, dDDdpl, dFFdpl, dPhidpl;
	FType dUdpl, dWdpl;
	FType tempSQR;

	FType dAdxl, dBdxl;

	FType **JacWOC, *AnsWOC, *XmWOC, *XmsWOC;

	FType Ci;

	if ((JacWOC=(FType **) malloc((2*Nc+2)*sizeof(FType *)))==NULL) TerM("Can not allocate memory for JacWOC matrix");
	for (i=0;i<((2*Nc+2));i++) {
		if ((JacWOC[i]=(FType *) malloc((2*Nc+4)*sizeof(FType)))==NULL) TerM("Can not allocate memory for JacWOC matrix");
	}

	if ((AnsWOC=(FType *) malloc((2*Nc+2)*sizeof(FType)))==NULL) {
		for (i=0; i<(2*Nc+2); i++) {
			free(JacWOC[i]);
		}
		free(JacWOC);
		TerM("Can not allocate memory for AnsWOC matrix");
	}

	if ((XmWOC=(FType *) malloc((2*Nc+2)*sizeof(FType)))==NULL) {
		for (i=0; i<(2*Nc+2); i++) {
			free(JacWOC[i]);
		}
		free(JacWOC);
		free(AnsWOC);
		TerM("Can not allocate memory for XmWOC matrix");
	}

	if ((XmsWOC=(FType *) malloc((2*Nc+2)*sizeof(FType)))==NULL) {
		for (i=0; i<(2*Nc+2); i++) {
			free(JacWOC[i]);
		}
		free(JacWOC);
		free(AnsWOC);
		free(XmWOC);
		TerM("Can not allocate memory for XmsWOC matrix");
	}


	XmWOC[2*Nc+1]=blockH[1][1][refL+1];
	XmWOC[2*Nc]=Ipress;
	for (i=0; i<Nc; i++){
		NRKi=(fluidProp[i][PCRIT]/XmWOC[2*Nc])*exp(5.37*(1+fluidProp[i][AC])*(1-fluidProp[i][TCRIT]/resTemp));		//Wilson Equation
		XmWOC[Nc+i]=NRKi*initComp[i];
		XmWOC[i]=initComp[i];
	}

	do {
		for (i=0; i<(2*Nc+2); i++)
			for (j=0; j<(2*Nc+2); j++) JacWOC[i][j]=0;

		SXm=0;
		SYm=0;
		for (i=0; i<Nc; i++){
			if (XmWOC[i]<0) XmWOC[i]=0;
			else if (XmWOC[i]>1) XmWOC[i]=1;
			if (XmWOC[Nc+i]<0) XmWOC[Nc+i]=0;
			else if (XmWOC[Nc+i]>1) XmWOC[Nc+i]=1;
			SXm+=XmWOC[i];
			SYm+=XmWOC[Nc+i];
			JacWOC[2*Nc][i]=1;
			JacWOC[2*Nc+1][Nc+i]=1;
		}

		if (XmWOC[2*Nc]<0) XmWOC[2*Nc]=PZERO;

		AnsWOC[2*Nc]=-SXm+1;
		AnsWOC[2*Nc+1]=-SYm+1;


		Al=0;
		Bl=0;
		Ag=0;
		Bg=0;
		for (i=0; i<Nc; i++) {
			Bl+=XmWOC[i]*fluidProp[i][EOS_B];
			Bg+=XmWOC[Nc+i]*fluidProp[i][EOS_B];
			for (j=0; j<Nc; j++) {
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A]);
				Al+=XmWOC[i]*XmWOC[j]*tempSQR;
				Ag+=XmWOC[Nc+i]*XmWOC[Nc+j]*tempSQR;
			}
		}

		Cal=Al;
		Cbl=Bl;
		Cag=Ag;
		Cbg=Bg;
		Al*=XmWOC[2*Nc]/(RGAS*RGAS*resTemp*resTemp);
		Bl*=XmWOC[2*Nc]/(RGAS*resTemp);
		Ag*=XmWOC[2*Nc]/(RGAS*RGAS*resTemp*resTemp);
		Bg*=XmWOC[2*Nc]/(RGAS*resTemp);


		dAdpg=Cag/(RGAS*RGAS*resTemp*resTemp);
		dBdpg=Cbg/(RGAS*resTemp);
		dAdpl=Cal/(RGAS*RGAS*resTemp*resTemp);
		dBdpl=Cbl/(RGAS*resTemp);



		if (SRK) {
			Ul=Bl;
			Wl=0;
			Ug=Bg;
			Wg=0;

			dUdpg=dBdpg;
			dWdpg=0;
			dUdpl=dBdpl;
			dWdpl=0;

			d1=1;
			d2=0;
		}
		else if (PR) {
			Ul=2*Bl;
			Wl=Bl;
			Ug=2*Bg;
			Wg=Bg;

			dUdpg=2*dBdpg;
			dWdpg=dBdpg;
			dUdpl=2*dBdpl;
			dWdpl=dBdpl;

			d1=1+SQRT2;
			d2=1-SQRT2;
		}

		Zl=Solve_Z(-(1+Bl-Ul), Al-Bl*Ul-Ul-Wl*Wl, -(Al*Bl-Bl*Wl*Wl-Wl*Wl), 'l');
		Zg=Solve_Z(-(1+Bg-Ug), Ag-Bg*Ug-Ug-Wg*Wg, -(Ag*Bg-Bg*Wg*Wg-Wg*Wg), 'g');


		dZdpg=-(dAdpg*(Zg-Bg)+dBdpg*(-Zg*Zg-Ug*Zg-Ag+Wg*Wg)+dUdpg*(Zg*Zg-Bg*Zg-Zg)+dWdpg*(-2*Wg*Zg+2*Bg*Wg+2*Wg))/(3*Zg*Zg-2*Zg*(1+Bg-Ug)+Ag-Bg*Ug-Ug-Wg*Wg);		//Modified
		dZdpl=-(dAdpl*(Zl-Bl)+dBdpl*(-Zl*Zl-Ul*Zl-Al+Wl*Wl)+dUdpl*(Zl*Zl-Bl*Zl-Zl)+dWdpl*(-2*Wl*Zl+2*Bl*Wl+2*Wl))/(3*Zl*Zl-2*Zl*(1+Bl-Ul)+Al-Bl*Ul-Ul-Wl*Wl);		//Modified


		for (i=0; i<Nc; i++) {
			tl=0;
			tg=0;
			for (j=0; j<Nc; j++) {	//FFSS
				tempSQR=sqrt(fluidProp[i][EOS_A]*fluidProp[j][EOS_A]);
				tl+=XmWOC[j]*tempSQR;
				tg+=XmWOC[Nc+j]*tempSQR;
			}
			t1l=fluidProp[i][EOS_B]/Cbl;
			t1g=fluidProp[i][EOS_B]/Cbg;

			EEl=2*tl/Cal-t1l;
			EEg=2*tg/Cag-t1g;

			DDl=Al/(Bl*(d1-d2));
			DDg=Ag/(Bg*(d1-d2));

			FFl=log((Zl+d2*Bl)/(Zl+d1*Bl));
			FFg=log((Zg+d2*Bg)/(Zg+d1*Bg));

			phiL=exp(t1l*(Zl-1)-log(Zl-Bl)+DDl*EEl*FFl);
			phiG=exp(t1g*(Zg-1)-log(Zg-Bg)+DDg*EEg*FFg);

			Ci=Fug0[i]*exp(fluidProp[i][MW]*G_ACC*(blockH[1][1][refL+1]-XmWOC[2*Nc+1])/(RGAS*resTemp));

			AnsWOC[i]=Ci-XmWOC[i]*XmWOC[2*Nc]*phiL;
			AnsWOC[Nc+i]=-XmWOC[i]*phiL+XmWOC[Nc+i]*phiG;

			dAAdpg=fluidProp[i][EOS_B]*dZdpg/Cbg;
			dAAdpl=fluidProp[i][EOS_B]*dZdpl/Cbl;

			dBBdpg=(dBdpg-dZdpg)/(Zg-Bg);
			dBBdpl=(dBdpl-dZdpl)/(Zl-Bl);

			dDDdpg=(Bg*dAdpg-Ag*dBdpg)/(Bg*Bg*(d1-d2));
			dDDdpl=(Bl*dAdpl-Al*dBdpl)/(Bl*Bl*(d1-d2));

			dFFdpg=(dZdpg+d2*dBdpg)/(Zg+d2*Bg)-(dZdpg+d1*dBdpg)/(Zg+d1*Bg);
			dFFdpl=(dZdpl+d2*dBdpl)/(Zl+d2*Bl)-(dZdpl+d1*dBdpl)/(Zl+d1*Bl);

			dPhidpg=phiG*(dAAdpg+dBBdpg+EEg*(dDDdpg*FFg+DDg*dFFdpg));
			dPhidpl=phiL*(dAAdpl+dBBdpl+EEl*(dDDdpl*FFl+DDl*dFFdpl));

			
			JacWOC[i][2*Nc+1]=Fug0[i]*fluidProp[i][MW]*G_ACC/(RGAS*resTemp)*exp(fluidProp[i][MW]*G_ACC*(blockH[1][1][refL+1]-XmWOC[2*Nc+1])/(RGAS*resTemp));
			JacWOC[i][2*Nc]=XmWOC[i]*(phiL+XmWOC[2*Nc]*dPhidpl);
			JacWOC[Nc+i][2*Nc]=XmWOC[i]*dPhidpl-XmWOC[Nc+i]*dPhidpg;
				
			

			//i: Row of Jac
			//k: Col of Jac
			for (k=0; k<Nc; k++) {
				dAdxl=0;
				dAdxg=0;
				for (n=0; n<Nc; n++) {		//SSFF
					tempSQR=sqrt(fluidProp[n][EOS_A]*fluidProp[k][EOS_A]);
					dAdxl+=XmWOC[n]*tempSQR;
					dAdxg+=XmWOC[Nc+n]*tempSQR;
					//This loop can be reduced
				}
				dAdxg*=2*XmWOC[2*Nc]/(RGAS*RGAS*resTemp*resTemp);
				dAdxl*=2*XmWOC[2*Nc]/(RGAS*RGAS*resTemp*resTemp);
				dBdxg=fluidProp[k][EOS_B]*XmWOC[2*Nc]/(RGAS*resTemp);				
				dBdxl=fluidProp[k][EOS_B]*XmWOC[2*Nc]/(RGAS*resTemp);

				if (SRK) {
					dUdxg=dBdxg;
					dWdxg=0;

					dUdxl=dBdxl;
					dWdxl=0;
				}
				else if (PR) {
					dUdxg=2*dBdxg;
					dWdxg=dBdxg;

					dUdxl=2*dBdxl;
					dWdxl=dBdxl;
				}


				dZdxl=-(dAdxl*(Zl-Bl)+dBdxl*(-Zl*Zl-Ul*Zl-Al+Wl*Wl)+dUdxl*(Zl*Zl-Bl*Zl-Zl)+dWdxl*(-2*Wl*Zl+2*Bl*Wl+2*Wl))/(3*Zl*Zl-2*Zl*(1+Bl-Ul)+Al-Bl*Ul-Ul-Wl*Wl);		//Modified
				dZdxg=-(dAdxg*(Zg-Bg)+dBdxg*(-Zg*Zg-Ug*Zg-Ag+Wg*Wg)+dUdxg*(Zg*Zg-Bg*Zg-Zg)+dWdxg*(-2*Wg*Zg+2*Bg*Wg+2*Wg))/(3*Zg*Zg-2*Zg*(1+Bg-Ug)+Ag-Bg*Ug-Ug-Wg*Wg);		//Modified
				
				dAAdxl=-fluidProp[k][EOS_B]*fluidProp[i][EOS_B]*(Zl-1)/(Cbl*Cbl)+fluidProp[i][EOS_B]*dZdxl/Cbl;
				dAAdxg=-fluidProp[k][EOS_B]*fluidProp[i][EOS_B]*(Zg-1)/(Cbg*Cbg)+fluidProp[i][EOS_B]*dZdxg/Cbg;

				dBBdxl=(dBdxl-dZdxl)/(Zl-Bl);
				dBBdxg=(dBdxg-dZdxg)/(Zg-Bg);

				dDDdxl=(dAdxl*Bl-dBdxl*Al)/(Bl*Bl*(d1-d2));
				dDDdxg=(dAdxg*Bg-dBdxg*Ag)/(Bg*Bg*(d1-d2));

				dFFdxl=(dZdxl+d2*dBdxl)/(Zl+d2*Bl)-(dZdxl+d1*dBdxl)/(Zl+d1*Bl);
				dFFdxg=(dZdxg+d2*dBdxg)/(Zg+d2*Bg)-(dZdxg+d1*dBdxg)/(Zg+d1*Bg);

				dEEdxl=2*(sqrt(fluidProp[i][EOS_A]*fluidProp[k][EOS_A])/Cal-dAdxl*RGAS*RGAS*resTemp*resTemp*tl/(XmWOC[2*Nc]*Cal*Cal))+fluidProp[i][EOS_B]*fluidProp[k][EOS_B]/(Cbl*Cbl);
				dEEdxg=2*(sqrt(fluidProp[i][EOS_A]*fluidProp[k][EOS_A])/Cag-dAdxg*RGAS*RGAS*resTemp*resTemp*tg/(XmWOC[2*Nc]*Cag*Cag))+fluidProp[i][EOS_B]*fluidProp[k][EOS_B]/(Cbg*Cbg);		//modified
				

				dPhidxg=phiG*(dAAdxg+dBBdxg+dDDdxg*EEg*FFg+DDg*dEEdxg*FFg+DDg*EEg*dFFdxg);
				dPhidxl=phiL*(dAAdxl+dBBdxl+dDDdxl*EEl*FFl+DDl*dEEdxl*FFl+DDl*EEl*dFFdxl);

				if (i==k) {
					JacWOC[i][k]=XmWOC[2*Nc]*(phiL+XmWOC[i]*dPhidxl);
					JacWOC[Nc+i][k]=phiL+XmWOC[i]*dPhidxl;
					JacWOC[Nc+i][Nc+k]=-phiG-XmWOC[Nc+i]*dPhidxg;
				}
				else {
					JacWOC[i][k]=XmWOC[2*Nc]*XmWOC[i]*dPhidxl;
					JacWOC[Nc+i][k]=XmWOC[i]*dPhidxl;
					JacWOC[Nc+i][Nc+k]=-XmWOC[Nc+i]*dPhidxg;
				}
				

			}


		}



		MKLS(JacWOC, AnsWOC, XmsWOC, 2*Nc+2);

		XmTol=0;
		for (i=0; i<(2*Nc+2); i++) {
			if (XmWOC[i]) XmTol+=fabs(XmsWOC[i]/XmWOC[i]);
			XmWOC[i]+= XmsWOC[i];			
		}
		XmTol/=(2*Nc+2);		

	} while ((XmTol>1e-5) || (fabs(XmsWOC[2*Nc])>1));



	for (i=0; i<(2*Nc+2); i++) {
		free(JacWOC[i]);
	}
	free(JacWOC);
	free(AnsWOC);
	free(XmWOC);
	free(XmsWOC);


	//return Xm[Nc];
}*/


void PJac(void) {
	register int Ix, Iy, Iz;
	int blockN, bSize;
	preConRow[0]=0;
	preConIndex[0]=0;

	blockN=16*Nc*Nc+31*Nc+17;
	totalRow=Nx*Ny*Nz*(2*Nc+4);


	totalJac=0;
	for (Iz=0; Iz<Nz; Iz++)
		for (Iy=0; Iy<Ny; Iy++)
			for (Ix=0; Ix<Nx; Ix++) {
				pJHolder[Ix][Iy][Iz]=totalJac;
				totalJac+=blockN;

				if (Ix==0) totalJac-=2*Nc*Nc+4*Nc+2;
				if (Ix==(Nx-1)) totalJac-=2*Nc*Nc+4*Nc+2;
				if (Iy==0) totalJac-=2*Nc*Nc+4*Nc+2;
				if (Iy==(Ny-1)) totalJac-=2*Nc*Nc+4*Nc+2;
				if (Iz==0) totalJac-=2*Nc*Nc+4*Nc+2;
				if (Iz==(Nz-1)) totalJac-=2*Nc*Nc+4*Nc+2;
				if (phaseStat[Ix][Iy][Iz]==1) {
					totalJac-=3*Nc*Nc+3*Nc+1;
					totalRow-=Nc+1;
					bSize=Nc+3;
				}
				else if (phaseStat[Ix][Iy][Iz]==-1) {
					totalJac-=3*Nc*Nc+4*Nc+1;
					totalRow-=Nc+1;
					bSize=Nc+3;
				}
				else bSize=2*Nc+4;

				preConRow[Iz*(Nx*Ny)+Iy*Nx+Ix+1]=preConRow[Iz*(Nx*Ny)+Iy*Nx+Ix]+bSize;
				preConIndex[Iz*(Nx*Ny)+Iy*Nx+Ix+1]=preConIndex[Iz*(Nx*Ny)+Iy*Nx+Ix]+bSize*(bSize-1)+MRMAXNONZERO;

			}
			preConCSRrow[preConRow[Nx*Ny*Nz]]=preConIndex[Nx*Ny*Nz];
			//totalJac=jacIndex;
}

void CalcPCoeff(void) {
	/*register int Ix, Iy, Iz;
	
	pCoeff=0;
	for (Iz=1; Iz<=Nz; Iz++)
		for (Iy=1; Iy<=Ny; Iy++)
			for (Ix=1; Ix<=Nx; Ix++) {
				if (P[Ix][Iy][Iz][1]>pCoeff) pCoeff=P[Ix][Iy][Iz][1];
			}*/
	pCoeff=AvgPReport();

	pCoeff*=1.5;
}


void ManageTSMarker(void) {
	register int i, j, k, n;
	FType temp;

	n=2*wellNO;
	//Sort
	for (i=0; i<n; i++)
		for (j=(i+1); j<=n; j++) if (TStepMarker[i]>TStepMarker[j]) {
				temp=TStepMarker[i];
				TStepMarker[i]=TStepMarker[j];
				TStepMarker[j]=temp;
			}

	j=0;
	for (i=0; i<n; i++) {
		if ((fabs(TStepMarker[i-j]-TStepMarker[i-j+1])<TSTEPTOL) || (!TStepMarker[i-j])) {
			for (k=(i-j); k<(n-j); k++) TStepMarker[k]=TStepMarker[k+1];
			j++;			
		}
	}

	TSMSize=n-j+1;

	if (Dt>=TStepMarker[0]) Dt=TStepMarker[0]/2.1;
}

#endif