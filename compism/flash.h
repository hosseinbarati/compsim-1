
#ifndef Globals_h
#define Globals_h
#include"Globals.h"
#endif

#define radical2 1.414213562

using namespace std;

double solve(double a1, double a2, double a3, bool Gas){
	double z[3];
	double q = (3 * a2 - pow(a1, 2)) / 9;
	double j = (9 * a1*a2 - 27 * a3 - 2 * pow(a1, 3)) / 54;
	double D = pow(q, 3) + pow(j, 2);
	if (D > 0){
		z[0] = pow(j + pow(D, 0.5), 1 / 3) + pow(j - pow(D, 0.5), 1 / 3) - a1 / 3;
		cout << "Z is equal to" << z[0] << endl;
	}
	else if (D < 0){
		double theta = acos(j / pow(-q, 1.5));
		z[0] = 2 * pow(-q, 0.5)*cos((theta / 3)) - a1 / 3;
		z[1] = 2 * pow(-q, 0.5)*cos(theta / 3 + (120)* 3.141592654 / 180) - a1 / 3;
		z[2] = 2 * pow(-q, 0.5)*cos(theta / 3 + (240)* 3.141592654 / 180) - a1 / 3;

		if (Gas){
			if ((z[0]>z[1]) && (z[0] > z[2]))
				z[0] = z[0];
			else if ((z[1]) > (z[2]))
				z[0] = z[1];
			else
				z[0] = z[2];
			cout << "Z gas is equal to = " << z[0] << endl;
		}

		else
		{
			if (((z[0]) > 0) && ((z[0]) < (z[1])) && ((z[0]) < (z[2])))
				z[0];
			else if (((z[1])>0) && ((z[1]) < (z[2])))
				z[0] = z[1];
			else
				z[0] = z[3];
			cout << "Z liquid is equal to = " << z[0] << endl;
		}
	}
	else{
		z[0] = 2 * pow(j, 1 / 3) - a1 / 3;
		z[1] = (-1)*pow(j, 1 / 3) - a1 / 3;
		if (Gas){
			if ((z[0]) > (z[1]))
				z[0] = z[0];
			else
				z[0] = z[1];
			cout << "Z gas is equal to = " << z[0] << endl;
		}
		else{
			if (((z[0]) > 0) && ((z[0]) < (z[1])))
				z[0] = z[0];
			else
				z[0] = z[1];
			cout << "Z liquid is equal to = " << z[0] << endl;
		}
	}
	return z[0];
}

double f(double l,double* compo,double* Km){
	double y = 0;
	for (int i = 0; i < Nc; i++){
		y += (((compo[i])*(1 - (Km[i]))) / (l + (1 - l)*Km[i]));
	}
	return y;
}
double df(double l,double* compo,double* Km){
	double y = 0;
	for (int i = 0; i < Nc; i++){
		y += (((compo[i])*(1 - Km[i])*(1 - Km[i])) / (pow((l + (1 - l)*Km[i]), 2)));
	}
	return y;
}

double* Approx(double* pc,double* tc,double* w,double p,double t){
	double* Kvalue = new double[Nc];
	for (int i = 0; i < Nc; i++){
		Kvalue[i] = ((pc[i]) / (p))*exp(5.37*(1 + w[i])*(1 - ((tc[i]) / t)));
	}
	return Kvalue;
}

bool condition(double L,double l){
	bool doit;
	for (int i = 0; i < Nc; i++)
	if (abs(L-l)>pow(10, -3)){
		doit = true;
	}
	return doit;
}

void firstflash(double L,double* compo,double* w,double t,double* tc,double* pc,double p,bool doit) {

	double* Xmo = new double[Nc];
	double* Xmg = new double[Nc];
	double* lambdam = new double[Nc];
	double* alpham = new double[Nc];
	double* am = new double[Nc];
	double* bm = new double[Nc];
	double bAlphaO = 0;
	double bAlphaG = 0;
	double aAlphaG = 0;
	double aAlphaO = 0;
	double AalphaO = 0;
	double AalphaG = 0;
	double BalphaO = 0;
	double BalphaG = 0;
	double sigmaG = 0;
	double l = 0;
	
	double* PhiG = new double[Nc];
	double* PhiO = new double[Nc];

	double* km = new double[Nc];
	double* Km = new double[Nc];

	bool gas = true;
	bool liquid = true;

	
	
	Km = Approx(pc, tc, w, p, t);
	do{
		l = L;
		if (L == 0){
			exit(-1);
		}
		
		if (!df(L, compo, Km)){
			cout << "Choose another L please.(because derivative is 0) this code must be modified later" << endl;
			exit(-1);
		}
		else{
			double L2;
			do {
				L2 = L;
				L += (-f(L, compo, Km)) / df(L, compo, Km);
			} while (abs(L - L2) > pow(10, -10));
			cout << "L=" << L << endl;
			if (L < 0){
				L = 0;
			}	
		}	
		if ((L>0)&(L < 1)){
			for (int i = 0; i < Nc; i++){
				Xmo[i] = compo[i] / (L + (1 - L)*Km[i]);
				Xmg[i] = Km[i] * Xmo[i];
				cout << "we have two phase" << endl;
			}
		}
		else if(!L){
			for (int i = 0; i < Nc; i++){
				Xmg[i] = compo[i];
				Xmo[i] = 0;
				liquid = false;
			}
			cout << "we have only one vapor phase" << endl;
		}
		
		for (int i = 0; i < Nc; i++)
		{
			lambdam[i] = 0.37464 + 1.5422*w[i] - 0.26992*pow(w[i], 2);
			alpham[i] = pow((1 + lambdam[i] * (1 - pow((t / tc[i]), 0.5))), 2);
			am[i] = ((0.457235*alpham[i] * pow((Rconst*tc[i]), 2)) / pc[i]);
			bm[i] = ((0.077796*Rconst*tc[i]) / (pc[i]));
		}
		for (int i = 0; i < Nc; i++){
			if (liquid){
				bAlphaO += Xmo[i] * bm[i];
			}
			if (gas){
				bAlphaG += Xmg[i] * bm[i];
			}
		}

		for (int i = 0; i < Nc; i++){
			for (int j = 0; j < Nc; j++){
				if (liquid){
					aAlphaO += Xmo[i] * Xmo[j] * pow((am[i] * am[j]), 0.5);
				}
				if (gas){
					aAlphaG += Xmg[i] * Xmg[j] * pow((am[i] * am[j]), 0.5);
				}
				
			}
		}

		if (liquid){
			double AalphaO = (aAlphaO*p) / pow((Rconst*t), 2);
			double BalphaO = (bAlphaO*p) / (Rconst*t);
		}
		if (gas){
			AalphaG = (aAlphaG*p) / pow((Rconst*t), 2);
			BalphaG = (bAlphaG*p) / (Rconst*t);
		}
		

		if (gas){
			double a1G = (BalphaG - 1);
			double a2G = (AalphaG - 2 * BalphaG - 3 * pow(BalphaG, 2));
			double a3G = (-(AalphaG*BalphaG) + pow(BalphaG, 2) + pow(BalphaG, 3));
			double Zg = solve(a1G, a2G, a3G, true);
			
			for (int i = 0; i < Nc; i++){
				for (int j = 0; j < Nc; j++){
					sigmaG += Xmg[j] * pow((am[i] * am[j]), 0.5);
				}
				PhiG[i] = exp((bm[i] / bAlphaG)*(Zg - 1) - log(Zg - BalphaG) - (AalphaG / (2 * radical2*BalphaG))*((2 / aAlphaG)*sigmaG - (bm[i] / bAlphaG))*log((Zg + (1 + radical2)*BalphaG) / (Zg - (1 - radical2)*BalphaG)));
			}

		}
		
		if (liquid){
			double a1O = (BalphaO - 1);
			double a2O = (AalphaO - 2 * BalphaO - 3 * pow(BalphaO, 2));
			double a3O = (-AalphaO*BalphaO + pow(BalphaO, 2) + pow(BalphaO, 3));
			double Zo = solve(a1O, a2O, a3O, false);
		}
		
		if (liquid){
			for (int i = 0; i < Nc; i++){
				cout << "phiO[" << i << "]= " << PhiO[i] << endl;
			}
		}

		if (gas){
			for (int i = 0; i < Nc; i++){
				cout << "phig[" << i << "]= " << PhiG[i] << endl;
			}
		}
		if (liquid&gas)
		for (int i = 0; i < Nc; i++){
			Km[i] = (PhiO[i]) / (PhiG[i]);
			cout << "km[" << i << "]= " << Km[i] << endl;
		}
		doit = condition(L, l);
	} while (doit);
	for(int i = 0; i < Nc; i++){
		cout << "Xmo[" << i << "] =" << Xmo[i] << endl;
		cout << "Xmg[" << i << "] =" << Xmg[i] << endl;
	}

}


