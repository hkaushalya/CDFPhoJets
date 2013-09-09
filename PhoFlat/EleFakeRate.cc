#include "PhoFlat/EleFakeRate.hh"
#include <iostream>
#include "TMath.h"

EleFakeRate::EleFakeRate()
{
}


//_________________________________ returns fake rate & systematics
//_________________________________ Main function _______________________
//__________ input:  eleEt= Pho->Etc(); (where Pho must be a TStnPhoton object)
//		     goodrun_bit=0 --- for goodrun_v13_pho_00.txt good run list
//		     goodrun_bit=1 --- for goodrun_v13_pho_lepphx_00.txt good run list		     
//		     phnx_bit=0 --- if you DON'T require photon-to-phoenix rejection
//		     phnx_bit=1 --- if you DO require photon-to-phoenix rejection
//____________________________________________________________________________________________________
//__________ output: wght_d=fake rate; wght_m=(fake rate) - (1 sigma); wght_p=(fake rate) + (1 sigma) 

void EleFakeRate::MyFakeRate(double eleEt, int goodrun_bit, int phnx_bit,
											double &wght_d, double &wght_m, double &wght_p)
{
	wght_d=0.0;
	wght_m=0.0;
	wght_p=0.0;

	if (goodrun_bit==0 && phnx_bit==0) { // with phnx, all runs
		phnx_par[0]=0.0; phnx_par[1]=0.0; phnx_par[2]=0.0; 
		phnx_parerr[0]=0.0; phnx_parerr[1]=0.0; phnx_parerr[2]=0.0; 
		phnx_corr_coeff[0][1]=phnx_corr_coeff[1][0]=0.0;
		phnx_corr_coeff[0][2]=phnx_corr_coeff[2][0]=0.0;
		phnx_corr_coeff[1][2]=phnx_corr_coeff[2][1]=0.0;
		phnx_corr_coeff[0][0]=phnx_corr_coeff[0][0]=1.0;
		phnx_corr_coeff[1][1]=phnx_corr_coeff[1][1]=1.0;
		phnx_corr_coeff[2][2]=phnx_corr_coeff[2][2]=1.0;
		fakerate_par[0]=-3.089; fakerate_par[1]=-0.0304; fakerate_par[2]=0.0059; 
		fakerate_parerr[0]=0.067; fakerate_parerr[1]=0.0055; fakerate_parerr[2]=0.0025;
		fakerate_par[3]=1.09;      // data/MC scale factor (using 40<Et<50)
		fakerate_parerr[3]=0.07; 
		fakerate_corr_coeff[0][0]=1.0; // diagonal covariance matrix correlation coefficients
		fakerate_corr_coeff[1][1]=1.0; 
		fakerate_corr_coeff[2][2]=1.0; 
		fakerate_corr_coeff[3][3]=1.0; 
		fakerate_corr_coeff[0][1]=fakerate_corr_coeff[1][0]=-0.514;
		fakerate_corr_coeff[0][2]=fakerate_corr_coeff[2][0]=0.237;
		fakerate_corr_coeff[1][2]=fakerate_corr_coeff[2][1]=-0.951;
		fakerate_corr_coeff[0][3]=fakerate_corr_coeff[3][0]=0.0;
		fakerate_corr_coeff[1][3]=fakerate_corr_coeff[3][1]=0.0;
		fakerate_corr_coeff[2][3]=fakerate_corr_coeff[3][2]=0.0;
	}
	
	if (goodrun_bit==1 && phnx_bit==0) { // with phnx, good Silicon //sam: isnt this and the following case same
		phnx_par[0]=0.0; phnx_par[1]=0.0; phnx_par[2]=0.0; 
		phnx_parerr[0]=0.0; phnx_parerr[1]=0.0; phnx_parerr[2]=0.0; 
		phnx_corr_coeff[0][1]=phnx_corr_coeff[1][0]=0.0;
		phnx_corr_coeff[0][2]=phnx_corr_coeff[2][0]=0.0;
		phnx_corr_coeff[1][2]=phnx_corr_coeff[2][1]=0.0;
		phnx_corr_coeff[0][0]=phnx_corr_coeff[0][0]=1.0;
		phnx_corr_coeff[1][1]=phnx_corr_coeff[1][1]=1.0;
		phnx_corr_coeff[2][2]=phnx_corr_coeff[2][2]=1.0;
		fakerate_par[0]=-3.086; fakerate_par[1]=-0.0270; fakerate_par[2]=0.0040; // if v11_lepphx
		fakerate_parerr[0]=0.063; fakerate_parerr[1]=0.0054; fakerate_parerr[2]=0.0031;
		fakerate_par[3]=1.03; 
		fakerate_parerr[3]=0.08; 
		fakerate_corr_coeff[0][0]=1.0; // diagonal covariance matrix correlation coefficients
		fakerate_corr_coeff[1][1]=1.0; 
		fakerate_corr_coeff[2][2]=1.0; 
		fakerate_corr_coeff[3][3]=1.0; 
		fakerate_corr_coeff[0][1]=fakerate_corr_coeff[1][0]=-0.286;
		fakerate_corr_coeff[0][2]=fakerate_corr_coeff[2][0]=0.0;
		fakerate_corr_coeff[1][2]=fakerate_corr_coeff[2][1]=-0.955;
		fakerate_corr_coeff[0][3]=fakerate_corr_coeff[3][0]=0.0;
		fakerate_corr_coeff[1][3]=fakerate_corr_coeff[3][1]=0.0;
		fakerate_corr_coeff[2][3]=fakerate_corr_coeff[3][2]=0.0;
	}
	
	//if (goodrun_bit==1 && phnx_bit==0) { // no phnx, all runs //sam::this shouldn't this be goodrun=0 phx=1 
	if (goodrun_bit==0 && phnx_bit==1) { // no phnx, all runs //sam::this shouldn't this be goodrun=0 phx=1 
		phnx_par[0]=0.394; phnx_par[1]=0.106; phnx_par[2]=24.54; 
		phnx_parerr[0]=0.004; phnx_parerr[1]=0.006; phnx_parerr[2]=0.32; 
		phnx_corr_coeff[0][1]=phnx_corr_coeff[1][0]=-0.482;
		phnx_corr_coeff[0][2]=phnx_corr_coeff[2][0]=0.450;
		phnx_corr_coeff[1][2]=phnx_corr_coeff[2][1]=-0.265;
		phnx_corr_coeff[0][0]=phnx_corr_coeff[0][0]=1.0;
		phnx_corr_coeff[1][1]=phnx_corr_coeff[1][1]=1.0;
		phnx_corr_coeff[2][2]=phnx_corr_coeff[2][2]=1.0;
		fakerate_par[0]=-3.089; fakerate_par[1]=-0.0304; fakerate_par[2]=0.0059; 
		fakerate_parerr[0]=0.067; fakerate_parerr[1]=0.0055; fakerate_parerr[2]=0.0025;
		fakerate_par[3]=1.42;      // data/MC scale factor (using 40<Et<50)
		fakerate_parerr[3]=0.15; 
		fakerate_corr_coeff[0][0]=1.0; // diagonal covariance matrix correlation coefficients
		fakerate_corr_coeff[1][1]=1.0; 
		fakerate_corr_coeff[2][2]=1.0; 
		fakerate_corr_coeff[3][3]=1.0; 
		fakerate_corr_coeff[0][1]=fakerate_corr_coeff[1][0]=-0.514;
		fakerate_corr_coeff[0][2]=fakerate_corr_coeff[2][0]=0.237;
		fakerate_corr_coeff[1][2]=fakerate_corr_coeff[2][1]=-0.951;
		fakerate_corr_coeff[0][3]=fakerate_corr_coeff[3][0]=0.0;
		fakerate_corr_coeff[1][3]=fakerate_corr_coeff[3][1]=0.0;
		fakerate_corr_coeff[2][3]=fakerate_corr_coeff[3][2]=0.0;
	}
	
	if (goodrun_bit==1 && phnx_bit==1) { // no phnx, good Silicon
		phnx_par[0]=0.421; phnx_par[1]=0.104; phnx_par[2]=24.61; 
		phnx_parerr[0]=0.004; phnx_parerr[1]=0.006; phnx_parerr[2]=0.32; 
		phnx_corr_coeff[0][1]=phnx_corr_coeff[1][0]=-0.487;
		phnx_corr_coeff[0][2]=phnx_corr_coeff[2][0]=0.426;
		phnx_corr_coeff[1][2]=phnx_corr_coeff[2][1]=-0.185;
		phnx_corr_coeff[0][0]=phnx_corr_coeff[0][0]=1.0;
		phnx_corr_coeff[1][1]=phnx_corr_coeff[1][1]=1.0;
		phnx_corr_coeff[2][2]=phnx_corr_coeff[2][2]=1.0;
		fakerate_par[0]=-3.086; fakerate_par[1]=-0.0270; fakerate_par[2]=0.0040; // if v11_lepphx
		fakerate_parerr[0]=0.063; fakerate_parerr[1]=0.0054; fakerate_parerr[2]=0.0031;
		fakerate_par[3]=1.19; 
		fakerate_parerr[3]=0.18; 
		fakerate_corr_coeff[0][0]=1.0; // diagonal covariance matrix correlation coefficients
		fakerate_corr_coeff[1][1]=1.0; 
		fakerate_corr_coeff[2][2]=1.0; 
		fakerate_corr_coeff[3][3]=1.0; 
		fakerate_corr_coeff[0][1]=fakerate_corr_coeff[1][0]=-0.286;
		fakerate_corr_coeff[0][2]=fakerate_corr_coeff[2][0]=0.0;
		fakerate_corr_coeff[1][2]=fakerate_corr_coeff[2][1]=-0.955;
		fakerate_corr_coeff[0][3]=fakerate_corr_coeff[3][0]=0.0;
		fakerate_corr_coeff[1][3]=fakerate_corr_coeff[3][1]=0.0;
		fakerate_corr_coeff[2][3]=fakerate_corr_coeff[3][2]=0.0;
	}
	
	if (goodrun_bit<0 || goodrun_bit>1 || phnx_bit<0 || phnx_bit>1) return;

	if (eleEt>=13.0) {
		double term1 = MyFakeRateFuncErr(eleEt);
		double multp = (1.0 - MyPhnxEffFunc(eleEt));
		double term2 = MyPhnxEffFuncErr(eleEt) * MyFakeRateFunc(eleEt);
		double my_err = MyShiftCorrFunc(eleEt) * sqrt(term1 * term1 * multp * multp + term2 * term2);
		wght_d = MyShiftCorrFunc(eleEt) * MyFakeRateFunc(eleEt) * (1.0 - MyPhnxEffFunc(eleEt));
		wght_m = wght_d - my_err;
		wght_p = wght_d + my_err;
		if (wght_d < 0.0) wght_d = 0.0;
		if (wght_m < 0.0) wght_m = 0.0;
		if (wght_p < 0.0) wght_p = 0.0;
	}

	return;
}


//____________ correction for Et(det)/Et(genp) difference
double EleFakeRate::MyShiftCorrFunc(double eleEt) 
{
	double arg1 = 0.0048 + exp(-3.031-0.027 * eleEt);
	double arg2 = 0.0056 + exp(-3.015-0.030 * eleEt);
	double value = arg1 / arg2;
	return value;
}

//__________________________________ phoenix efficiency function
double EleFakeRate::MyPhnxEffFunc(double eleEt)
{
	double arg = phnx_par[1] * (phnx_par[2] - eleEt);
	double value = phnx_par[0] * TMath::Erfc(arg);
	if(value<0.0 || value>1.0) value = 0.0;
	return value;
}

//__________________________________ uncertainty on phoenix efficiency
double EleFakeRate::MyPhnxEffFuncErr(double eleEt)
{
	double arg = phnx_par[1] * (phnx_par[2]-eleEt);
	double erfc_drv = 2.0 * exp(-1.0 * arg * arg) / TMath::Pi();
	double term1 = phnx_parerr[0] * TMath::Erfc(arg);
	double term2 = phnx_parerr[1] * phnx_par[0] * (phnx_par[2]-eleEt) * erfc_drv;
	double term3 = phnx_parerr[2] * phnx_par[0] * phnx_par[1] * erfc_drv;
	double term4 = 2.0 * phnx_corr_coeff[0][1] * phnx_parerr[0] * phnx_parerr[1] * TMath::Erfc(arg)
						* erfc_drv * phnx_par[0] * (phnx_par[2]-eleEt);
	double term5 = 2.0 * phnx_corr_coeff[0][2] * phnx_parerr[0] * phnx_parerr[2] * TMath::Erfc(arg)
						* erfc_drv * phnx_par[0] * phnx_par[1];
	double term6 = 2.0 * phnx_corr_coeff[1][2] * phnx_parerr[1] * phnx_parerr[2]
						* arg * erfc_drv * erfc_drv * phnx_par[0] * phnx_par[0];
	double term = term1 * term1 + term2 * term2 + term3 * term3 + term4 + term5 + term6;
	double err = 0.0;
	if(term>0.0) err = sqrt(term);
	return err;
}


//__________________ returns fake rate
double EleFakeRate::MyFakeRateFunc(double eleEt)
{
	double fkrt = 0.0;
	if (eleEt >= 13.0) {
		double arg1 = fakerate_par[0] + fakerate_par[1] * eleEt;
		double lin1 = fakerate_par[2];
		double scale1 = fakerate_par[3];
		fkrt = scale1 * (exp(arg1) + lin1);
	}  
	return fkrt;
}

//________________________________________ returns fake rates
double EleFakeRate::MyFakeRateFuncErr(double eleEt)
{
	double err=0.0;
	double _expo=0.0;
	
	if (fakerate_par[3] > 0.0) {
		_expo=MyFakeRateFunc(eleEt)/fakerate_par[3];
	}
	double mult1 = _expo - fakerate_par[2];
	double term1 = _expo * _expo * fakerate_parerr[3] * fakerate_parerr[3];
	double sub_term1 = mult1 * mult1 * (fakerate_parerr[0] * fakerate_parerr[0]
							+ 2.0 * fakerate_corr_coeff[1][0] * fakerate_parerr[0] * fakerate_parerr[1] * eleEt
							+ fakerate_parerr[1] * fakerate_parerr[1] * eleEt * eleEt);
	double sub_term2 = 2.0 * fakerate_parerr[2] * mult1 * (fakerate_corr_coeff[2][0] * fakerate_parerr[0]
							+ fakerate_corr_coeff[2][1] * fakerate_parerr[1] * eleEt);
	double sub_term3 = fakerate_parerr[2] * fakerate_parerr[2];
	double term2 = fakerate_par[3] * fakerate_par[3] * (sub_term1 + sub_term2 + sub_term3);
	if ((term1 + term2)>0.0) err = sqrt(term1 + term2);
	return err;
}



