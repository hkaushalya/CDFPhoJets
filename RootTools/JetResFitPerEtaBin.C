/////////////////////////////////////////////////////////
// This is Sasha's JET energy resolution fit method.   //
// You can derive JER parametirization from the fits   //
// to be put in the JetFilter module to model the JER. //
// See elog#1227.                                      //
// I renamed this macro from myJetResFit_fnl_response.C//
// This fits for one eta bin at a time.                //
// The 2nd script JetResFitResponse.C (renamed from    //
// myJetResFit_fnl_response_v4.C) does the fitting for //
// a particular eta bin.                               //
// The ROOT file JER_hadE_sample1_dR01_072507.root is  //
// required to make the fits.                          //
/////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TFolder.h"
#include "TMinuit.h"
#include "TGaxis.h"
#include "TStyle.h"
#include <sstream>
#include <vector>
#include "TLegend.h"
#include <algorithm>
#include <cmath>
#include "TFrame.h"

//extern void JetResFitPerEPerEtaBin(int e_bin_id=0, int eta_bin_id=0,
//		int folder_id=0, TVirtualPad* pad =0);

using namespace std;

//----- Another new fit function: sqrt(A0/Pt+A2) * step function combined
double fitFunc3_new(double *x, double *par)
{
	double fitval = 0.0;
	if (x[0] > 0.0)
	{
      double p1 = par[0] / x[0] + par[1];
      double p2 = par[2] / x[0] + par[3];
		
		double ws = 1 / (1 + exp((x[0]-par[4])/25.));
		fitval = p2 * (1 - ws) + p1 * ws;
	}

	return fitval;
}

// //-------------- fit function: poly3 * step func + poly3 * (step func)
// this is the new fit function used to describe the final
// Gausian Mean, Landau MPV, and Gaus Normalization
// See elog#
double dJayFunc1(double *x, double *par)
{

	double p1 = par[0]  + par[1] * x[0] + par[2] * x[0] * x[0];
	double p2 = par[3]  + par[4] * x[0] + par[5] * x[0] * x[0];
	//double p2 = par[3]  + par[4] * x[0];

	double ws = 1 / (1 + exp((x[0]-par[6])/10.5));
	double fitval = p1 * ws + p2 * (1 - ws);
	return fitval;
}

// //-------------- fit function: [Const*Gaus(-x/(1+x))+Landau(-x/(1+x))]/(1+Const)
double myFunc1(double *x, double *par)
{
  double fitval 	= 0.0;
  double arg 		= -x[0]/(x[0]+1.0);
  double arg_L 	= -x[0]/(x[0]+1.0);
  double mean 		= par[1];
  double sigmaG 	= par[2];
  double mpv 		= par[3];
  double sigmaL 	= par[4];
  double normL 	= par[5];
  if (normL < 0.0) return -1.0E6;
  double f1 = normL / (1.0 + normL) * TMath::Gaus(arg,mean,sigmaG);
  double f2 = TMath::Landau(arg_L,mpv,sigmaL) / (1.0+normL);
  fitval = par[0] * (f1 + f2);
  return fitval;
}

//------ Gauss mean: A0+A1*Pt+A3/sqrt(Pt) ---- default
double fitFunc0(double *x, double *par)
{
  double fitval = 0.0;
  if (x[0] > 0.0) fitval = par[0] + par[1] * x[0] + par[2] / x[0];
  else fitval = -1.0E6;
  return fitval;
}
//------ Gauss mean: A0+A1*Pt+A3/sqrt(Pt) ---- default
double fitFunc0_syst(double *x, double *par)
{
  double fitval=0.0;
  if (x[0] > 0.0) fitval = par[0] + par[1] * x[0] + par[2] / sqrt(x[0]);
  else fitval = -1.0E6;
  return fitval;
}
//------ Gauss norm: A0+A1*Pt+A3/ln(Pt) ---- default
double fitFunc5(double *x, double *par)
{
  double fitval = 0.0;
  if (x[0] > 0.0) fitval = (par[0] + par[1] * sqrt(x[0])) / x[0] + par[2];
  else fitval = -1.0E6;
  return fitval;
}

//----- Gauss Sigma fit function: sqrt(A0*A0/Pt+A1*A1/(Pt*Pt)+A2*A2)
double fitFunc1(double *x, double *par)
{
  double fitval = 0.0;
  if (x[0] > 0.0)
    {
      double arg = 0.0;
      arg = (par[0] / x[0] + par[1] / (x[0] * x[0]) + par[2]);
      fitval = sqrt(arg);
    }
  else fitval = -1.0E6;
  return fitval;
}
// //----- Landau Sigma fit function: sqrt(A0/Pt+A2) ---- try as default
double fitFunc3(double *x, double *par)
{
  double fitval = 0.0;
  if (x[0] > 0.0)
    {
      double arg = 0.0;
      arg = par[0] / x[0] + par[1];
      fitval = sqrt(arg);
    }
  else fitval = -1.0E6;
  return fitval;
}
//----- Landau Sigma fit function: A0/sqrt(Pt)+A2 ---- try as systematics 
double fitFunc3_syst(double *x, double *par)
{
  double fitval = 0.0;
  if (x[0] > 0.0)
    {
      double arg = 0.0;
      arg = par[0] / sqrt(x[0]) + par[1];
      fitval = arg;
    }
  else fitval = -1.0E6;
  return fitval;
}


//=====================================================================
//-------------- Driving script as of 07/24/07
// --- fits to all fJetRes_TrueBal for particular eta_det bin 
//
// Fit function is [Const*Gaus(-x/(1+x))+Landau(-x/(1+x))]/(1+Const)
//
// 6 plots: 5 fit param +fitMPV
//_____________________________________________________________________
std::vector<TF1*> JetResFitPerEtaBin(const int eta_bin_id=0, const bool bNewfits=1, const int folder_id=0)
		
{
	if (eta_bin_id <0 || eta_bin_id>14)
	{
		std::cout << "Invalid Eta bin. must be between 0 and 14. Returning"<< std::endl;
		exit (0);
	}

	if (! (folder_id == 0 || folder_id == 3 || folder_id == 4))
	{
		std::cout << "Invalid folder selection. must be 1,3 or 4. Returning" << std::endl;
		exit(0);
	}


	// to store list of functions to be trturned
	std::vector<TF1*> funcList;

	gROOT->Reset();

	char foldername[200];
	char fObject[200];

	if (folder_id==0) sprintf(foldername,"JetRes_def");
	//if (folder_id==1) sprintf(foldername,"JetRes_syst1");
	//if (folder_id==2) sprintf(foldername,"JetRes_syst2");
	if (folder_id==3) sprintf(foldername,"JetRes_def_nvx1");
	if (folder_id==4) sprintf(foldername,"JetRes_def_nvx2");

	double gauss_mean[100];
	double gauss_mean_er[100];
	double gauss_sigma[100];
	double gauss_sigma_er[100];
	double landau_mean[100];
	double landau_mean_er[100];
	double landau_sigma[100];
	double landau_sigma_er[100];
	double landau_norm[100];
	double landau_norm_er[100];
	double fit_chi2[100];
	double fit_chi2_er[100];
	double avePt[100];
	double avePt_er[100];
	double histo_Nentry[100];
	int fit_status[100];

	double maxFit_X[100]; // MPV for fit function
	double histoRMS[100]; // RMS of histogram 

	int h_bin_first=0;
	int h_bin_last=100;


	//________________________________________________________________________________
	//------------------- stopped here -----------------------------------------------
	//________________________________________________________________________________

	TH1F *h1[100];

	gDirectory->pwd();

	//-------------- matching in dR<0.1
	TFile myFile1("JER_hadE_sample1_dR01_072507.root","r");

	TFolder* myFolder1 = (TFolder* ) myFile1.FindObjectAny("Ana");

	//how/what is this hist??? there is one for each Eta bin
	sprintf(fObject,"MyJetFilter/Hist/%s/fJetRes_AveJetE[%i]",foldername,eta_bin_id); // get object name
	TProfile *_avePt_histo= (TProfile*) myFolder1->FindObjectAny(fObject);

	// get all the hists of different E bins for this eta bin.
	// initialize the fit var arrays
	for (int i = 0; i < 100; i++)
	{
		sprintf(fObject,"MyJetFilter/Hist/%s/fJetRes_TrueBal[%i][%i]",foldername,eta_bin_id,i); // get object name
		h1[i]           = (TH1F*) myFolder1->FindObjectAny(fObject);
		avePt[i]        = _avePt_histo->GetBinContent(i+1);
		avePt_er[i]     = _avePt_histo->GetBinError(i+1);
		histo_Nentry[i] = h1[i]->GetEntries();
		//if (i == 2)	std::cout << "i,avgPt, entries = " << i << "\t" << avePt[i] <<"\t" << histo_Nentry[i] << std::endl;
		
		if (histo_Nentry[i] > 99) h_bin_last = i;
		
		gauss_mean[i]    = 0.0;
		gauss_mean_er[i] = 0.0;
		gauss_sigma[i]   = 0.0;
		gauss_sigma_er[i]= 0.0;
		landau_mean[i]   = 0.0;
		landau_mean_er[i]= 0.0;
		landau_sigma[i]  = 0.0;
		landau_sigma_er[i] = 0.0;
		landau_norm[i]   = 0.0;
		landau_norm_er[i]= 0.0;
		fit_chi2[i]      = 0.0;
		fit_chi2_er[i]   = 0.0; 
		fit_status[i]    = 0;

		maxFit_X[i] = 0.0;
		histoRMS[i] = 0.01;
	}

	gROOT->cd();
	gDirectory->pwd();
	gStyle->SetAxisColor(kBlue);
	gStyle->SetPadTickY(1);
	gStyle->SetPadTickX(1);
	gStyle->SetOptFit(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetTitleFontSize(0.035);
	gStyle->SetStatFontSize(0.03);
	gStyle->SetTitleFont(2);


	
	//here we fit each energy bin separately (we fit twice)
	//and get the fit parameters needed to make the final 5 hists.
	//why is necessary to build new hist out of the original?
	//why not format(rebin) it and fit it straight?
	//we wont change it(saved hist) by doing do!
	for (int i=0; i<100; i++)
	{
		/*
		std::cout<<"---------------------------------------------"<<std::endl;
		std::cout<<"     Fitting bin="<<i<<std::endl;
		std::cout<<"_____________________________________________"<<std::endl;
		*/

		if (histo_Nentry[i] >= 10000.0) h1[i]->Rebin(2);
		if (histo_Nentry[i] >= 1000.0 && histo_Nentry[i] < 10000.0) h1[i]->Rebin(3);
		if (histo_Nentry[i] >= 200.0 && histo_Nentry[i] < 1000.0) h1[i]->Rebin(3);
		if (histo_Nentry[i] >= 100.0 && histo_Nentry[i] < 200.0) h1[i]->Rebin(4);
		if (histo_Nentry[i] < 100.0) h1[i]->Rebin(4);
		
		if (eta_bin_id < 14 && histo_Nentry[i] < 100.0)
		{
			continue; // skipping bins with Nevnt<100

		}
		if (eta_bin_id == 14 && histo_Nentry[i] < 80.0)
		{
			std::cout << "ETABIN=14 HIST_ENTRY = " << histo_Nentry[i] << std::endl;
			continue; // skipping bins with Nevnt<100
		}
		

		int   nbins  	= h1[i]->GetNbinsX();
		float binwidth = h1[i]->GetBinWidth(1);
		float lowedge  = (h1[i]->GetBinCenter(1))-binwidth/2.0;
		float highedge = lowedge + nbins * binwidth;

		char title[200];
		char name[200];

		sprintf(name,"bal");
		if (eta_bin_id < 14) sprintf(title,"E^{had}/E^{det}-1.0, %2.1f<|#eta_{det}^{jet}|<%2.1f and %i<E^{jet}<%i",0.2*eta_bin_id,0.2*(eta_bin_id+1),5*i,5*(i+1));
		else sprintf(title,"E^{had}/E^{det}-1.0, |#eta_{det}^{jet}|>%2.1f and %i<E^{jet}<%i",0.2*eta_bin_id,5*i,5*(i+1));

		TH1F* _h = new TH1F(name,title,nbins,lowedge,highedge);
		_h->Sumw2();

		float cut1 = -1.0;
		float cut2 = 2.0;
		float val  = 0.0;
		float err  = 0.0;
		int bin_first = 0;
		int bin_last  = 0;

		for (int j=0; j<nbins; j++)
		{
			val = h1[i]->GetBinContent(j+1);
			err = h1[i]->GetBinError(j+1);
			_h->SetBinContent(j+1,val);
			_h->SetBinError(j+1,err);
			if (bin_first == 0 && val>0.0) bin_first = j + 1;
			//       if(val>0.0) bin_last=i+1; 
			if (val > 0.0 && h1[i]->GetBinCenter(j+1) < 3.5) bin_last = j + 1; 
		}
		
		//these are the x axis limits
		cut1 = h1[i]->GetBinCenter(bin_first);
		cut2 = h1[i]->GetBinCenter(bin_last);

		int maxbin   = h1[i]->GetMaximumBin();
		double hight = h1[i]->GetBinContent(maxbin);
		double mean_histo = h1[i]->GetMean();
		double rms_histo  = h1[i]->GetRMS();

		TF1 *fitFun1 = new TF1("fit1",myFunc1,cut1-binwidth,cut2+binwidth,6);

		if (eta_bin_id < 14)
		{
			fitFun1->SetParameter(0,0.9*hight);			// norm 
			fitFun1->SetParameter(1,mean_histo);		// mean of Gauss
			fitFun1->SetParameter(2,0.7*rms_histo);	// Sigma of Gauss
			fitFun1->SetParameter(3,mean_histo);		// MPV of Landau
			fitFun1->SetParameter(4,0.5*rms_histo);	// Sigma of Landau
			fitFun1->SetParameter(5,0.2);					// Landau norm
		} else {
			fitFun1->SetParameter(0,0.9*hight);			// norm 
			fitFun1->FixParameter(1,0.0);					// mean of Gauss
			fitFun1->FixParameter(2,0.0);					// Sigma of Gauss
			fitFun1->SetParameter(3,mean_histo);		// MPV of Landau
			fitFun1->SetParameter(4,0.5*rms_histo);	// Sigma of Landau
			fitFun1->FixParameter(5,0.0);					// Landau norm
		}
		
		fitFun1->SetParName(0,"Gauss: const");
		fitFun1->SetParName(1,"Gauss: mean");
		fitFun1->SetParName(2,"Gauss: #sigma");
		fitFun1->SetParName(3,"Landau: MPV");
		fitFun1->SetParName(4,"Landau: #sigma");
		fitFun1->SetParName(5,"Gaus: const");

		_h->Fit(fitFun1,"BQ","",cut1,cut2);

		//---------- re-fitting second time --------------------------
		if (eta_bin_id<14)
		{
			fitFun1->SetParameter(0,fitFun1->GetParameter(0)); // norm 
			fitFun1->SetParameter(1,fitFun1->GetParameter(1)); // mean of Gauss
			fitFun1->SetParameter(2,fitFun1->GetParameter(2)); // Sigma of Gauss
			fitFun1->SetParameter(3,fitFun1->GetParameter(3)); // MPV of Landau
			fitFun1->SetParameter(4,fitFun1->GetParameter(4)); // Sigma of Landau
			fitFun1->SetParameter(5,fitFun1->GetParameter(5)); // Landau norm
		} else {
			fitFun1->SetParameter(0,fitFun1->GetParameter(0)); // norm 
			fitFun1->FixParameter(1,0.0); // mean of Gauss
			fitFun1->FixParameter(2,0.0); // Sigma of Gauss
			fitFun1->SetParameter(3,fitFun1->GetParameter(3)); // MPV of Landau
			fitFun1->SetParameter(4,fitFun1->GetParameter(4)); // Sigma of Landau
			fitFun1->FixParameter(5,0.0); // Landau norm	  
		}

		_h->Fit(fitFun1,"EBQ","",cut1,cut2);
		
		//---------- end of re-fitting -------------------------------

		if ((gMinuit->fCstatu).Contains("SUCCESSFUL")==true) fit_status[i]=3;
		if ((gMinuit->fCstatu).Contains("PROBLEMS")==true) fit_status[i]=2;
		if ((gMinuit->fCstatu).Contains("FAILURE")==true) fit_status[i]=0;

		//std::cout<<" Minuit Error Matrix Status ="<<gMinuit->fCstatu<<std::endl;
		//std::cout<<" Minuit Error Matrix Status ="<<fit_status[i]<<std::endl;
		//       gMinuit->mnmatu(1); 
		
		// after the fit, get all the info about the fit parameters 
		maxFit_X[i] 	    = fitFun1->GetMaximumX(-1.0,1.0);
		histoRMS[i] 	    = rms_histo;
		gauss_mean[i] 	    = fitFun1->GetParameter(1);
		gauss_mean_er[i]   = fitFun1->GetParError(1);
		gauss_sigma[i]     = fitFun1->GetParameter(2);
		gauss_sigma_er[i]  = fitFun1->GetParError(2);
		landau_mean[i]     = fitFun1->GetParameter(3);
		landau_mean_er[i]  = fitFun1->GetParError(3);
		landau_sigma[i] 	 = fitFun1->GetParameter(4);
		landau_sigma_er[i] = fitFun1->GetParError(4);
		landau_norm[i] 	 = fitFun1->GetParameter(5);
		landau_norm_er[i]  = fitFun1->GetParError(5);

		
		

		if (fitFun1->GetNDF() > 0)
			fit_chi2[i] = (1.0 * fitFun1->GetChisquare())/(1.0 * fitFun1->GetNDF());
		else fit_chi2[i] = -1.0;

		delete fitFun1;
		delete _h;
	}


	// NOW WE ARE GENERATING THE FINAL 5 HISTS (GAUS MEAN, GAUS SIGMA, etc)
	
	//TCanvas* c1 = new TCanvas("c1","Jet energy resolution study",200,10,1200,900);
	TCanvas* c1 = new TCanvas("c1","Jet energy resolution study",200,10,1100,850);
	c1->Divide(2,3);

	char title[200];
	char name[200]; 

	//--------------------- Gauss mean ---------------------------------------------
	sprintf(name,"h1");
	if (eta_bin_id<14) sprintf(title,"Energy dependence of Gauss mean for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
	else sprintf(title,"Energy dependence of Gauss mean for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
	TH1F* h1_fnl= new TH1F(name,title,500,0.0,500.0);
	h1_fnl->Sumw2();
	//--------------------- Gauss sigma ---------------------------------------------
	sprintf(name,"h2");
	if (eta_bin_id<14) sprintf(title,"Energy dependence of Gauss #sigma for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
	else sprintf(title,"Energy dependence of Gauss #sigma for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
	TH1F* h2_fnl= new TH1F(name,title,500,0.0,500.0);
	h2_fnl->Sumw2();
	//--------------------- Landau MPV ---------------------------------------------
	sprintf(name,"h3");
	if (eta_bin_id<14) sprintf(title,"Energy dependence of Landau MPV for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
	else sprintf(title,"Energy dependence of Landau MPV for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
	TH1F* h3_fnl= new TH1F(name,title,500,0.0,500.0);
	h3_fnl->Sumw2();
	//--------------------- Landau sigma ---------------------------------------------
	sprintf(name,"h4");
	if (eta_bin_id<14) sprintf(title,"Energy dependence of Landau #sigma for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
	else sprintf(title,"Energy dependence of Landau #sigma for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
	TH1F* h4_fnl= new TH1F(name,title,500,0.0,500.0);
	h4_fnl->Sumw2();
	//--------------------- Gauss normalization ---------------------------------------------
	sprintf(name,"h5");
	if (eta_bin_id<14) sprintf(title,"Energy dependence of Gauss normalization for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
	else sprintf(title,"Energy dependence of Gauss normalization for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
	TH1F* h5_fnl= new TH1F(name,title,500,0.0,500.0);
	h5_fnl->Sumw2();
	//--------------------- Fit MPV ---------------------------------------------
	sprintf(name,"h6");
	if (eta_bin_id<14) sprintf(title,"Energy dependence of MPV of Fit function for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
	else sprintf(title,"Energy dependence of MPV of Fit function for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
	TH1F* h6_fnl= new TH1F(name,title,500,0.0,500.0);
	h6_fnl->Sumw2();
	//--------------------- Filling up histograms ---------------------------------------------
	int last_E_bin=0;
	int first_E_bin=100;


	double fit_cut1 = 3.0;
	double fit_cut2 = 500.0;
	if (eta_bin_id < 5) fit_cut1 = 12.0;
	if (eta_bin_id > 4 && eta_bin_id < 8) fit_cut1 = 17.0;
	if (eta_bin_id == 8) fit_cut1 = 22.0;
	if (eta_bin_id == 9) fit_cut1 = 27.0;
	if (eta_bin_id == 10) fit_cut1 = 32.0;
	if (eta_bin_id == 11) fit_cut1 = 37.0;
	if (eta_bin_id == 12) fit_cut1 = 47.0;
	if (eta_bin_id == 13) fit_cut1 = 57.0;
	//if(eta_bin_id == 14) fit_cut1 = 112.0;
	if (eta_bin_id == 14) fit_cut1 = 102.0;


	//here we fill the final 5 hists (GAUS MEAN, GAUS SIGMA etc)
	
	std::cout << "\n=========== Apply quality cuts to Eta Bin = " << eta_bin_id << std::endl;
	
	for (int i = 0; i < 100; i++)  //100 energy bins from 0-500GeV
	{
		//here we exclude some E bins that disrupts final fits
		//I assume Sasha had closer look at these plots and varibles.
		
		bool good_fit = true;
		if (eta_bin_id == 2 && gauss_mean[i] > 0.1) good_fit = false;
		if (eta_bin_id < 10 && gauss_sigma[i] < 0.05) good_fit = false;
		if (eta_bin_id == 7 && (landau_norm[i] > 0.4 || gauss_mean[i] < -0.02)) good_fit = false; 
		if (eta_bin_id > 8 && eta_bin_id<13 && landau_mean[i]>0.0) good_fit = false; 
		
		if (eta_bin_id < 13)
		{	  
			if (fabs(landau_mean[i] / histoRMS[i]) > 1.1) good_fit = false;
			if (gauss_sigma[i] > 1.0 * histoRMS[i] || gauss_sigma[i] < 0.2*histoRMS[i]) good_fit = false;
			if (landau_sigma[i] > 0.9 * histoRMS[i] || landau_sigma[i] < 0.1*histoRMS[i]) good_fit = false;
			if (gauss_sigma_er[i] < 0.01 * gauss_sigma[i] || gauss_sigma_er[i] > 0.3*gauss_sigma[i]) good_fit = false;
			if (landau_sigma_er[i] < 0.01 * landau_sigma[i] || landau_sigma_er[i] > 0.3*landau_sigma[i]) good_fit = false;
			if (landau_mean_er[i] > 0.06) good_fit = false;
			if (gauss_mean_er[i] > 0.15) good_fit = false;
			if (landau_norm_er[i] < 0.05 * landau_norm[i]) good_fit = false;
		}
		
		if (eta_bin_id == 13)
		{
			if (landau_norm[i] > 0.2) good_fit = false;
			if (gauss_sigma[i] < 0.06) good_fit = false;
			if (gauss_mean[i] < -0.1) good_fit = false;
		}
		
		if (eta_bin_id < 14 && histo_Nentry[i] < 100) good_fit = false;
		if (eta_bin_id == 14 && (histo_Nentry[i] < 80 || landau_mean[i] < -0.17)) good_fit = false;
		if (fit_status[i] != 3) good_fit = false;		//fit was not SUCCESSFULL!
		if (avePt[i] < fit_cut1) good_fit = false;
		
		// if this enrgy bin is a good one to make into the final plot ...
	
		//this is a condition for my last test to use exact bin values for <60GeV
		//and use a fit to >60 region.
		if (i<=12) good_fit = true;  //include all points below 60GeV temporily
		if (good_fit == true)
		{ 
			// Gauss mean
			//need this NAN check to avoid hist not being drawn.
			//some of the fit values returns NAN and it causes trouble 
			//when hist is drawn. (it does not draw).
			//this happens when I include all the data points <60GeV. 
			//
			if (! std::isnan(gauss_mean[i]))
			{
				h1_fnl->Fill(avePt[i],gauss_mean[i]);
				h1_fnl->SetBinError(h1_fnl->FindBin(avePt[i]),gauss_mean_er[i]);
			}

			// Gauss sigma	  
			if (! std::isnan(gauss_sigma[i]))
			{
				h2_fnl->Fill(avePt[i],gauss_sigma[i]);
				h2_fnl->SetBinError(h2_fnl->FindBin(avePt[i]),gauss_sigma_er[i]);
			}
			//
			// Landau mpv
			if (! std::isnan(landau_mean[i]))
			{
				h3_fnl->Fill(avePt[i],landau_mean[i]);
				h3_fnl->SetBinError(h3_fnl->FindBin(avePt[i]),landau_mean_er[i]);
			}
			
			// Landau sigma
			if (! std::isnan(landau_sigma[i]))
			{
				h4_fnl->Fill(avePt[i],landau_sigma[i]);
				h4_fnl->SetBinError(h4_fnl->FindBin(avePt[i]),landau_sigma_er[i]);
			}
			
			// Gauss norm
			if (! std::isnan(landau_norm[i]))
			{
				h5_fnl->Fill(avePt[i],landau_norm[i]);
				h5_fnl->SetBinError(h5_fnl->FindBin(avePt[i]),landau_norm_er[i]);
			}
			
			// Fit MPV (comment: uncertainties need to be correctly estimated)
			h6_fnl->Fill(avePt[i],maxFit_X[i]);
			h6_fnl->SetBinError(h6_fnl->FindBin(avePt[i]),landau_mean_er[i]);
			last_E_bin = i;
			
			fit_cut2 = avePt[last_E_bin] + 1.0;
			
			if (i < first_E_bin) 
			{
				first_E_bin = i;
				if (fit_cut1 < avePt[first_E_bin] - 2.0) fit_cut1 = avePt[first_E_bin] - 1.0;
			}

		} else { // if good_fit
			if (i<6)
			{
				std::cout << __LINE__ << "E :: EXCLUDED THIS POINT: E bin = " << i <<std::endl;
				//new TCanvas();
				//JetResFitPerEPerEtaBin(i, eta_bin_id,0, gPad);
			}
		}
		
	}

	//temporary
	fit_cut1 = 60;
	if (eta_bin_id == 13) fit_cut1 = 57.0;
	if (eta_bin_id == 14) fit_cut1 = 102.0;

	//--------------------- Gauss mean ---------------------------------------------
	c1->cd(1);
	//new TCanvas();
	gPad->SetTicky();
	gPad->SetTickx();
	h1_fnl->SetLineColor(2);
	h1_fnl->SetLineWidth(1);
	h1_fnl->SetMarkerColor(2);
	h1_fnl->SetMarkerStyle(20);
	h1_fnl->SetMarkerSize(0.6);
	h1_fnl->SetMinimum(-0.2);
	h1_fnl->SetMaximum(0.2);

	if (bNewfits)
	{
		TF1 *JayFitFunc1 = new TF1("JayFitFunc1",dJayFunc1,2,500,7);
		JayFitFunc1->SetParameter(0,0.0005);
		JayFitFunc1->SetParameter(1,0.0045);
		JayFitFunc1->SetParameter(2,-0.00005);
		JayFitFunc1->SetParameter(3,0.09); 
		JayFitFunc1->SetParameter(4,-0.0004);
		JayFitFunc1->SetParameter(5,0.003);
		JayFitFunc1->SetParameter(6,60);
		JayFitFunc1->SetParLimits(6,50,70);

		h1_fnl->Fit(JayFitFunc1,"EBQ","",fit_cut1,fit_cut2);

	} else {

		TF1 *param_fitFun1=new TF1("fit_1",fitFunc0,2.0,500.0,3);
		param_fitFun1->SetParameter(0,0.09); // A0 param: term=A0
		param_fitFun1->SetParameter(1,-0.0002); // param A1: term=A1*Pt
		param_fitFun1->SetParameter(2,-1.0); // param A2: term=A2/Pt
		h1_fnl->Fit(param_fitFun1,"EBQ","",fit_cut1,fit_cut2);

			for (int i=1; i<=h1_fnl->GetNbinsX(); ++i)
			{
				if (h1_fnl->GetBinContent(i))
					std::cout << "bin [" << i << ", " 
						<< h1_fnl->GetBinContent(i) << "]" 
					<< "     --- Binedges = " 
					<<  h1_fnl->GetBinLowEdge(i) << ", " 
					<< h1_fnl->GetXaxis()->GetBinUpEdge(i) 
						<< std::endl;
			}
	}
	

	c1->cd(1);
	h1_fnl->Draw("AP");
	h1_fnl->GetXaxis()->CenterTitle(kTRUE);
	h1_fnl->GetYaxis()->CenterTitle(kTRUE);
	h1_fnl->GetXaxis()->SetTitle("average E^{jet}");
	h1_fnl->GetYaxis()->SetTitle("Gauss mean");

	//return funcList;
	
	//--------------------- Gauss sigma ---------------------------------------------
	c1->cd(2);
	//new TCanvas();
	gPad->SetTicky();
	gPad->SetTickx();
	h2_fnl->SetLineColor(2);
	h2_fnl->SetLineWidth(1);
	h2_fnl->SetMarkerColor(2);
	h2_fnl->SetMarkerStyle(20);
	h2_fnl->SetMarkerSize(0.6);
	h2_fnl->SetMinimum(0.0);
	h2_fnl->SetMaximum(0.3);

	if (bNewfits)
	{
		std::cout << "Using fitFunc3_new for Gsigma" << std::endl;
		TF1 *stitchTest2 = new TF1("stitchTest2",fitFunc3_new,2,500,5);
		stitchTest2->SetParameter(0,1.6);
		stitchTest2->SetParameter(1,0.00005);
		stitchTest2->SetParameter(2,2); 
		stitchTest2->SetParameter(3,0.004);
		stitchTest2->SetParameter(4,30);
		stitchTest2->SetParLimits(4,20,40.0);

		h2_fnl->Fit(stitchTest2,"EBQ","",fit_cut1,fit_cut2);

	} else {
		TF1 *param_fitFun2=new TF1("fit_2",fitFunc3,2.0,500.0,2);
		param_fitFun2->SetParameter(0,1.0); // A0 param: term=A0/Pt
		param_fitFun2->SetParameter(1,0.0025); // param A1: term=A1
		h2_fnl->Fit(param_fitFun2,"EBQ","",fit_cut1,fit_cut2);
	}

	h2_fnl->Draw("AP");
	h2_fnl->GetXaxis()->CenterTitle(kTRUE);
	h2_fnl->GetYaxis()->CenterTitle(kTRUE);
	h2_fnl->GetXaxis()->SetTitle("average E^{jet}");
	h2_fnl->GetYaxis()->SetTitle("Gauss #sigma");

	
	//--------------------- Landau MPV ---------------------------------------------
	c1->cd(3);
	//new TCanvas();
	gPad->SetTicky();
	gPad->SetTickx();
	h3_fnl->SetLineColor(2);
	h3_fnl->SetLineWidth(1);
	h3_fnl->SetMarkerColor(2);
	h3_fnl->SetMarkerStyle(20);
	h3_fnl->SetMarkerSize(0.6);
	h3_fnl->SetMinimum(-0.2);
	h3_fnl->SetMaximum(0.2);

	if (bNewfits)
	{
		TF1 *JayFitFunc3 = new TF1("JayFitFunc3",dJayFunc1,2,500,7);
		JayFitFunc3->SetParameter(0,0.0005);
		JayFitFunc3->SetParameter(1,0.0045);
		JayFitFunc3->SetParameter(2,-0.00005);
		JayFitFunc3->SetParameter(3,0.09); 
		JayFitFunc3->SetParameter(4,-0.0004);
		JayFitFunc3->SetParameter(5,0.003);
		JayFitFunc3->SetParameter(6,60);
		JayFitFunc3->SetParLimits(6,50,70);
		h3_fnl->Fit(JayFitFunc3,"EBQ","",fit_cut1,fit_cut2);

	} else {
		
		TF1 *param_fitFun3=new TF1("fit_3",fitFunc0,2.0,500.0,3);
		param_fitFun3->SetParameter(0,-0.01); // A0 param: term=A0
		param_fitFun3->SetParameter(1,0.0002); // param A1: term=A1*Pt
		param_fitFun3->SetParameter(2,-1.0); // param A2: term=A2/sqrt(Pt)
		h3_fnl->Fit(param_fitFun3,"EBQ","",fit_cut1,fit_cut2);
	}
	
	h3_fnl->Draw("AP");
	h3_fnl->GetXaxis()->CenterTitle(kTRUE);
	h3_fnl->GetYaxis()->CenterTitle(kTRUE);
	h3_fnl->GetXaxis()->SetTitle("average E^{jet}");
	h3_fnl->GetYaxis()->SetTitle("Landau mpv");
	
	//--------------------- Landau sigma ---------------------------------------------
	c1->cd(4);
	//new TCanvas();
	gPad->SetTicky();
	gPad->SetTickx();
	h4_fnl->SetLineColor(2);
	h4_fnl->SetLineWidth(1);
	h4_fnl->SetMarkerColor(2);
	h4_fnl->SetMarkerStyle(20);
	h4_fnl->SetMarkerSize(0.6);
	h4_fnl->SetMinimum(0.0);
	h4_fnl->SetMaximum(0.2);
	
	if (bNewfits)
	{
		std::cout << "Using fitFunc3_new for LandauSigma" << std::endl;
		TF1 *stitchTest4 = new TF1("stitchTest4",fitFunc3_new,2,500,5);
		stitchTest4->SetParameter(0,1.6);
		//stitchTest4->SetParameter(1,0.0045);
		stitchTest4->SetParameter(1,0.00005);
		stitchTest4->SetParameter(2,2); 
		stitchTest4->SetParameter(3,0.004);
		stitchTest4->SetParameter(4,30.0);
		stitchTest4->SetParLimits(4,20,40.0);
		h4_fnl->Fit(stitchTest4,"EBQ","",fit_cut1,fit_cut2);

	} else {

		TF1 *param_fitFun4=new TF1("fit_4",fitFunc3,2.0,500.0,2);
		if(eta_bin_id<8)
		{
			param_fitFun4->SetParameter(0,0.2); // A0 param: term=A0/sqrt(Pt)
			param_fitFun4->SetParameter(1,0.003); // param A1: term=A1/Pt
		} else {
			param_fitFun4->SetParameter(0,0.1); // A0 param: term=A0/sqrt(Pt)
			param_fitFun4->SetParameter(1,0.0005); // param A1: term=A1/Pt
		}

		h4_fnl->Fit(param_fitFun4,"EBQ","",fit_cut1,fit_cut2);
	}

	h4_fnl->Draw("AP");
	h4_fnl->GetXaxis()->CenterTitle(kTRUE);
	h4_fnl->GetYaxis()->CenterTitle(kTRUE);
	h4_fnl->GetXaxis()->SetTitle("average E^{jet}");
	h4_fnl->GetYaxis()->SetTitle("Landau #sigma");

	//--------------------- Gauss normalization ---------------------------------------------
	//new TCanvas();
	c1->cd(5);
	gPad->SetTicky();
	gPad->SetTickx();
	h5_fnl->SetLineColor(2);
	h5_fnl->SetLineWidth(1);
	h5_fnl->SetMarkerColor(2);
	h5_fnl->SetMarkerStyle(20);
	h5_fnl->SetMarkerSize(0.6);
	h5_fnl->SetMinimum(0.0);
	h5_fnl->SetMaximum(1.0);
	
	if (bNewfits)
	{
		TF1 *JayFitFunc5 = new TF1("JayFitFunc5",dJayFunc1,2,500,7);
		JayFitFunc5->SetParameter(0,0.0005);
		JayFitFunc5->SetParameter(1,0.0045);
		JayFitFunc5->SetParameter(2,0.00005);
		JayFitFunc5->SetParameter(3,0.09); 
		JayFitFunc5->SetParameter(4,-0.0004);
		JayFitFunc5->SetParameter(5,0.003);
		JayFitFunc5->SetParameter(6,60);
		JayFitFunc5->SetParLimits(6,50,70);
		h5_fnl->Fit(JayFitFunc5,"EBQ","",fit_cut1,fit_cut2);

	} else {

		TF1 *param_fitFun5=new TF1("fit_5",fitFunc5,2.0,500.0,3);
		param_fitFun5->SetParameter(0,0.22); // A0 param: term=A0
		param_fitFun5->SetParameter(1,-0.00043); // param A1: term=A1*Pt
		param_fitFun5->SetParameter(2,0.22); // param A2: term=A2/sqrt(Pt)
		h5_fnl->Fit(param_fitFun5,"EB","",fit_cut1,fit_cut2);

	}

	h5_fnl->Draw("AP");
	h5_fnl->GetXaxis()->CenterTitle(kTRUE);
	h5_fnl->GetYaxis()->CenterTitle(kTRUE);
	h5_fnl->GetXaxis()->SetTitle("average E^{jet}");
	h5_fnl->GetYaxis()->SetTitle("Gauss norm");

	//--------------------- Fit MPV ---------------------------------------------
	c1->cd(6);
	gPad->SetTicky();
	gPad->SetTickx();
	h6_fnl->SetLineColor(2);
	h6_fnl->SetLineWidth(1);
	h6_fnl->SetMarkerColor(2);
	h6_fnl->SetMarkerStyle(20);
	h6_fnl->SetMarkerSize(0.6);
	h6_fnl->SetMinimum(-0.2);
	h6_fnl->SetMaximum(0.2);
	h6_fnl->Draw("AP");
	h6_fnl->GetXaxis()->CenterTitle(kTRUE);
	h6_fnl->GetYaxis()->CenterTitle(kTRUE);
	h6_fnl->GetXaxis()->SetTitle("average E^{jet}");
	h6_fnl->GetYaxis()->SetTitle("Fit mpv");


	//funcList.push_back(JayFitFunc1);
	//funcList.push_back(param_fitFun1);
	//funcList.push_back(param_fitFun2);
	//funcList.push_back(JayFitFunc3);
	//funcList.push_back(param_fitFun3);
	//funcList.push_back(param_fitFun4);
	//funcList.push_back(JayFitFunc5);
	//funcList.push_back(param_fitFun5);
	
	//bring the focus out of the last pad to the whole canvas and print
	c1->cd();
	std::stringstream giffile, epsfile;
	giffile << "~/TMP/JERfit_EtaBin"<< eta_bin_id <<  ".gif";
	epsfile << "~/TMP/JERfit_EtaBin"<< eta_bin_id <<  ".eps";
	//c1->Print(giffile.str().c_str());
	c1->Print(epsfile.str().c_str(),"Preview");
	
	//c1->Close();
	
	std::cout<<"Closing Files"<<std::endl;

	myFile1.cd();
	gDirectory->pwd();
	myFile1.Close();

	return funcList;
}


void JetResFitAllEtaBins(const bool bNewfits=1)
{
	const int folder_id=0;
	for (int etabin=0; etabin<=14;++etabin)
	{
		JetResFitPerEtaBin(etabin, bNewfits, folder_id);
	}
}

//------------------------------------------------------------------------------------
// ONLY FOR TEST MET MODEL IMPLEMENTATION
//------------------------------------------------------------------------------------
//this is to test what I copied into METmodel is correct and
//is same as I have it in this code. I am going over lay the
//two fit fuctions, one in this final fit and once I can get 
//MEtmodel code.

//-------- Jet Energy Resolution: centered at zero
// //-------------- fit function: [Const*Gaus(-x/(1+x))+Landau(-x/(1+x))]/(1+Const)
double MyJER(double* x, double* par) {
	double value = 0.0;
	double arg    = -x[0]/(x[0] + 1.0);
	double arg_L  = -x[0]/(x[0] + 1.0);
	double mean   = par[0];
	double sigmaG = par[1];
	double mpv    = par[2];
	double sigmaL = par[3];
	double normL  = par[4];
	if (normL < 0.0) normL = 0.0;
	double f1 = normL / (1.0 + normL) * TMath::Gaus(arg, mean, sigmaG);
	double f2 = TMath::Landau(arg_L, mpv, sigmaL) / (1.0 + normL);
	value = f1 + f2;
	if (value<0.0) return 0.0;
	return value;
}

//________________ returns JER param or its uncertainty
double JERparam(int code, int i, int j)
{
	double jer_param[15][19]; // [eta][param]
	double jer_param_er[15][19]; // [eta][param]
	double value=0.0;

	// NEW parameters 07-14-2009
	jer_param[0][0]=0.119343;  jer_param_er[0][0]=0.00314468;
	jer_param[0][1]=-0.00105953;  jer_param_er[0][1]=6.95625e-05;
	jer_param[0][2]=3.45175e-05;  jer_param_er[0][2]=1.36251e-06;
	jer_param[0][3]=0.0306632;  jer_param_er[0][3]=0.00256322;
	jer_param[0][4]=-0.000131089;  jer_param_er[0][4]=8.70229e-06;
	jer_param[0][5]=1.03604;  jer_param_er[0][5]=0.0162871;
	jer_param[0][6]=0.000363984;  jer_param_er[0][6]=9.99763e-05;
	jer_param[0][7]=-0.0361686;  jer_param_er[0][7]=0.0185456;
	jer_param[0][8]=0.00128243;  jer_param_er[0][8]=4.24584e-05;
	jer_param[0][9]=-1.52357e-06;  jer_param_er[0][9]=1.00586e-05;
	jer_param[0][10]=0.0215701;  jer_param_er[0][10]=0.00392601;
	jer_param[0][11]=0.000136695;  jer_param_er[0][11]=2.02636e-05;
	jer_param[0][12]=0.157408;  jer_param_er[0][12]=0.00948233;
	jer_param[0][13]=0.00745406;  jer_param_er[0][13]=0.000132492;
	jer_param[0][14]=-0.0272588;  jer_param_er[0][14]=0.0342609;
	jer_param[0][15]=0.00449702;  jer_param_er[0][15]=0.000111858;
	jer_param[0][16]=-8.08661e-06;  jer_param_er[0][16]=2.2715e-05;
	jer_param[0][17]=0.0368499;  jer_param_er[0][17]=0.0109694;
	jer_param[0][18]=0.000367415;  jer_param_er[0][18]=5.70766e-05;
	jer_param[1][0]=0.0743446;  jer_param_er[1][0]=0.0183135;
	jer_param[1][1]=0.000247436;  jer_param_er[1][1]=4.80135e-05;
	jer_param[1][2]=2.87667e-06;  jer_param_er[1][2]=8.85104e-06;
	jer_param[1][3]=0.0904459;  jer_param_er[1][3]=0.00177181;
	jer_param[1][4]=-0.000188996;  jer_param_er[1][4]=8.124e-06;
	jer_param[1][5]=0.670063;  jer_param_er[1][5]=0.0119732;
	jer_param[1][6]=0.00308811;  jer_param_er[1][6]=9.34083e-05;
	jer_param[1][7]=-0.147575;  jer_param_er[1][7]=0.0106425;
	jer_param[1][8]=0.00551904;  jer_param_er[1][8]=0.000567946;
	jer_param[1][9]=-6.02653e-05;  jer_param_er[1][9]=6.89239e-06;
	jer_param[1][10]=-0.00799666;  jer_param_er[1][10]=0.0012872;
	jer_param[1][11]=1.88714e-05;  jer_param_er[1][11]=5.32563e-06;
	jer_param[1][12]=0.137285;  jer_param_er[1][12]=0.00491652;
	jer_param[1][13]=0.0020658;  jer_param_er[1][13]=3.46802e-05;
	jer_param[1][14]=0.0494562;  jer_param_er[1][14]=0.0278234;
	jer_param[1][15]=0.00238938;  jer_param_er[1][15]=0.00206893;
	jer_param[1][16]=3.70568e-05;  jer_param_er[1][16]=2.8803e-05;
	jer_param[1][17]=0.235407;  jer_param_er[1][17]=0.00688744;
	jer_param[1][18]=-0.00045021;  jer_param_er[1][18]=2.72024e-05;
	jer_param[2][0]=-0.0275628;  jer_param_er[2][0]=0.0193417;
	jer_param[2][1]=0.00348816;  jer_param_er[2][1]=4.0658e-05;
	jer_param[2][2]=-2.58846e-05;  jer_param_er[2][2]=9.08771e-06;
	jer_param[2][3]=0.088725;  jer_param_er[2][3]=0.00178948;
	jer_param[2][4]=-0.000219068;  jer_param_er[2][4]=7.83811e-06;
	jer_param[2][5]=0.655545;  jer_param_er[2][5]=0.0103474;
	jer_param[2][6]=0.00300742;  jer_param_er[2][6]=8.73698e-05;
	jer_param[2][7]=-0.175947;  jer_param_er[2][7]=0.0147411;
	jer_param[2][8]=0.00545286;  jer_param_er[2][8]=4.14128e-05;
	jer_param[2][9]=-5.21046e-05;  jer_param_er[2][9]=7.58811e-06;
	jer_param[2][10]=-0.00656557;  jer_param_er[2][10]=0.00133665;
	jer_param[2][11]=-2.77298e-06;  jer_param_er[2][11]=5.36955e-06;
	jer_param[2][12]=0.139571;  jer_param_er[2][12]=0.00510727;
	jer_param[2][13]=0.00227237;  jer_param_er[2][13]=3.75247e-05;
	jer_param[2][14]=0.105951;  jer_param_er[2][14]=0.0507535;
	jer_param[2][15]=0.00108184;  jer_param_er[2][15]=0.000209741;
	jer_param[2][16]=4.83183e-05;  jer_param_er[2][16]=3.42105e-05;
	jer_param[2][17]=0.222155;  jer_param_er[2][17]=0.00708712;
	jer_param[2][18]=-0.000384935;  jer_param_er[2][18]=2.76368e-05;
	jer_param[3][0]=0.170694;  jer_param_er[3][0]=0.0184655;
	jer_param[3][1]=-0.00400434;  jer_param_er[3][1]=5.12281e-05;
	jer_param[3][2]=4.51275e-05;  jer_param_er[3][2]=9.53991e-06;
	jer_param[3][3]=0.094768;  jer_param_er[3][3]=0.00222356;
	jer_param[3][4]=-0.000184981;  jer_param_er[3][4]=9.72646e-06;
	jer_param[3][5]=0.731283;  jer_param_er[3][5]=0.0128143;
	jer_param[3][6]=0.00341218;  jer_param_er[3][6]=0.00011154;
	jer_param[3][7]=-0.203812;  jer_param_er[3][7]=0.010309;
	jer_param[3][8]=0.00677219;  jer_param_er[3][8]=0.000524366;
	jer_param[3][9]=-6.90694e-05;  jer_param_er[3][9]=6.44028e-06;
	jer_param[3][10]=0.00968912;  jer_param_er[3][10]=0.00153634;
	jer_param[3][11]=-2.1932e-05;  jer_param_er[3][11]=6.07106e-06;
	jer_param[3][12]=0.179229;  jer_param_er[3][12]=0.00515981;
	jer_param[3][13]=0.0024204;  jer_param_er[3][13]=4.113e-05;
	jer_param[3][14]=0.0861196;  jer_param_er[3][14]=0.0238194;
	jer_param[3][15]=2.85701e-05;  jer_param_er[3][15]=0.00158696;
	jer_param[3][16]=5.64572e-05;  jer_param_er[3][16]=2.28046e-05;
	jer_param[3][17]=0.200478;  jer_param_er[3][17]=0.00701648;
	jer_param[3][18]=-0.000353516;  jer_param_er[3][18]=2.73206e-05;
	jer_param[4][0]=0.11804;  jer_param_er[4][0]=0.0184238;
	jer_param[4][1]=-0.00232059;  jer_param_er[4][1]=5.23777e-05;
	jer_param[4][2]=4.23751e-05;  jer_param_er[4][2]=1.05065e-05;
	jer_param[4][3]=0.104218;  jer_param_er[4][3]=0.00235565;
	jer_param[4][4]=-0.000118992;  jer_param_er[4][4]=9.97604e-06;
	jer_param[4][5]=0.795373;  jer_param_er[4][5]=0.016243;
	jer_param[4][6]=0.00635758;  jer_param_er[4][6]=0.000146956;
	jer_param[4][7]=-0.192637;  jer_param_er[4][7]=0.013221;
	jer_param[4][8]=0.00701698;  jer_param_er[4][8]=4.63014e-05;
	jer_param[4][9]=-7.554e-05;  jer_param_er[4][9]=8.76423e-06;
	jer_param[4][10]=0.0207201;  jer_param_er[4][10]=0.00200845;
	jer_param[4][11]=-3.80897e-05;  jer_param_er[4][11]=7.4085e-06;
	jer_param[4][12]=0.212192;  jer_param_er[4][12]=0.00725209;
	jer_param[4][13]=0.00349351;  jer_param_er[4][13]=6.06418e-05;
	jer_param[4][14]=0.138705;  jer_param_er[4][14]=0.03147;
	jer_param[4][15]=-0.00334811;  jer_param_er[4][15]=0.0020746;
	jer_param[4][16]=9.2358e-05;  jer_param_er[4][16]=2.96098e-05;
	jer_param[4][17]=0.255835;  jer_param_er[4][17]=0.00874165;
	jer_param[4][18]=-0.000434759;  jer_param_er[4][18]=3.08506e-05;
	jer_param[5][0]=0.0744578;  jer_param_er[5][0]=0.0246458;
	jer_param[5][1]=0.000218351;  jer_param_er[5][1]=4.40605e-05;
	jer_param[5][2]=-6.74408e-06;  jer_param_er[5][2]=7.54392e-07;
	jer_param[5][3]=0.0866891;  jer_param_er[5][3]=0.00163009;
	jer_param[5][4]=-7.76345e-05;  jer_param_er[5][4]=6.14488e-06;
	jer_param[5][5]=0.706037;  jer_param_er[5][5]=0.0137402;
	jer_param[5][6]=0.0066424;  jer_param_er[5][6]=0.000110101;
	jer_param[5][7]=-0.230767;  jer_param_er[5][7]=0.0214145;
	jer_param[5][8]=0.00789996;  jer_param_er[5][8]=4.30248e-05;
	jer_param[5][9]=-9.88319e-05;  jer_param_er[5][9]=7.74325e-07;
	jer_param[5][10]=0.00343542;  jer_param_er[5][10]=0.00250523;
	jer_param[5][11]=-3.3107e-05;  jer_param_er[5][11]=8.44783e-06;
	jer_param[5][12]=0.163174;  jer_param_er[5][12]=0.00755835;
	jer_param[5][13]=0.00291679;  jer_param_er[5][13]=6.87186e-05;
	jer_param[5][14]=0.176869;  jer_param_er[5][14]=0.0678367;
	jer_param[5][15]=-0.00259215;  jer_param_er[5][15]=0.000196127;
	jer_param[5][16]=2.45926e-05;  jer_param_er[5][16]=3.83147e-06;
	jer_param[5][17]=0.479562;  jer_param_er[5][17]=0.0175298;
	jer_param[5][18]=-0.000678321;  jer_param_er[5][18]=5.60546e-05;
	jer_param[6][0]=0.00938799;  jer_param_er[6][0]=0.0214178;
	jer_param[6][1]=0.00286763;  jer_param_er[6][1]=5.49692e-05;
	jer_param[6][2]=-4.64671e-05;  jer_param_er[6][2]=9.12845e-07;
	jer_param[6][3]=0.110373;  jer_param_er[6][3]=0.00285928;
	jer_param[6][4]=-0.000322779;  jer_param_er[6][4]=1.1545e-05;
	jer_param[6][5]=0.84678;  jer_param_er[6][5]=0.0199186;
	jer_param[6][6]=0.00418531;  jer_param_er[6][6]=0.00022283;
	jer_param[6][7]=-0.261284;  jer_param_er[6][7]=0.0164498;
	jer_param[6][8]=0.00726371;  jer_param_er[6][8]=4.72184e-05;
	jer_param[6][9]=-7.91647e-05;  jer_param_er[6][9]=1.01505e-05;
	jer_param[6][10]=-0.0219306;  jer_param_er[6][10]=0.00243926;
	jer_param[6][11]=9.96534e-05;  jer_param_er[6][11]=1.15067e-05;
	jer_param[6][12]=0.0758951;  jer_param_er[6][12]=0.00708642;
	jer_param[6][13]=0.00415232;  jer_param_er[6][13]=8.11982e-05;
	jer_param[6][14]=0.332473;  jer_param_er[6][14]=0.0524506;
	jer_param[6][15]=-0.0139085;  jer_param_er[6][15]=0.000198343;
	jer_param[6][16]=0.000222125;  jer_param_er[6][16]=3.63124e-05;
	jer_param[6][17]=0.168506;  jer_param_er[6][17]=0.00873082;
	jer_param[6][18]=-0.00037344;  jer_param_er[6][18]=3.87933e-05;
	jer_param[7][0]=-0.021659;  jer_param_er[7][0]=0.0183932;
	jer_param[7][1]=0.00176846;  jer_param_er[7][1]=6.69075e-05;
	jer_param[7][2]=-1.70584e-05;  jer_param_er[7][2]=1.08547e-05;
	jer_param[7][3]=0.0550046;  jer_param_er[7][3]=0.00234551;
	jer_param[7][4]=-0.000172672;  jer_param_er[7][4]=1.18965e-05;
	jer_param[7][5]=0.889453;  jer_param_er[7][5]=0.015712;
	jer_param[7][6]=0.00152891;  jer_param_er[7][6]=0.000140114;
	jer_param[7][7]=-0.355076;  jer_param_er[7][7]=0.0175825;
	jer_param[7][8]=0.0126345;  jer_param_er[7][8]=6.849e-05;
	jer_param[7][9]=-0.000141083;  jer_param_er[7][9]=1.11936e-05;
	jer_param[7][10]=-0.0399519;  jer_param_er[7][10]=0.00299535;
	jer_param[7][11]=0.000117128;  jer_param_er[7][11]=1.67435e-05;
	jer_param[7][12]=0.170338;  jer_param_er[7][12]=0.00800203;
	jer_param[7][13]=0.00119121;  jer_param_er[7][13]=7.89355e-05;
	jer_param[7][14]=0.362403;  jer_param_er[7][14]=0.067433;
	jer_param[7][15]=-0.0112584;  jer_param_er[7][15]=0.00400909;
	jer_param[7][16]=0.000155376;  jer_param_er[7][16]=5.09827e-05;
	jer_param[7][17]=0.27005;  jer_param_er[7][17]=0.0196289;
	jer_param[7][18]=-0.000376317;  jer_param_er[7][18]=0.00011391;
	jer_param[8][0]=-0.147359;  jer_param_er[8][0]=0.030064;
	jer_param[8][1]=0.00860543;  jer_param_er[8][1]=9.03519e-05;
	jer_param[8][2]=-0.00012039;  jer_param_er[8][2]=1.40453e-05;
	jer_param[8][3]=0.0405336;  jer_param_er[8][3]=0.00134233;
	jer_param[8][4]=-6.27006e-05;  jer_param_er[8][4]=5.70287e-06;
	jer_param[8][5]=0.87612;  jer_param_er[8][5]=0.0135608;
	jer_param[8][6]=0.000862571;  jer_param_er[8][6]=8.72723e-05;
	jer_param[8][7]=-0.340401;  jer_param_er[8][7]=0.0332533;
	jer_param[8][8]=0.0106851;  jer_param_er[8][8]=0.000115048;
	jer_param[8][9]=-0.000118701;  jer_param_er[8][9]=1.85101e-05;
	jer_param[8][10]=-0.0537098;  jer_param_er[8][10]=0.00254475;
	jer_param[8][11]=0.000105504;  jer_param_er[8][11]=1.11912e-05;
	jer_param[8][12]=0.218132;  jer_param_er[8][12]=0.00804389;
	jer_param[8][13]=-0.000115495;  jer_param_er[8][13]=5.26047e-05;
	jer_param[8][14]=0.291315;  jer_param_er[8][14]=0.145622;
	jer_param[8][15]=-0.00487384;  jer_param_er[8][15]=0.000420369;
	jer_param[8][16]=1.11175e-05;  jer_param_er[8][16]=0.0001095;
	jer_param[8][17]=0.452202;  jer_param_er[8][17]=0.0306732;
	jer_param[8][18]=-0.00016021;  jer_param_er[8][18]=0.000159679;
	jer_param[9][0]=-0.299768;  jer_param_er[9][0]=0.0536478;
	jer_param[9][1]=0.0158593;  jer_param_er[9][1]=0.000136572;
	jer_param[9][2]=-0.00021132;  jer_param_er[9][2]=2.47158e-05;
	jer_param[9][3]=0.0301886;  jer_param_er[9][3]=0.00206139;
	jer_param[9][4]=-1.59033e-05;  jer_param_er[9][4]=9.69668e-06;
	jer_param[9][5]=0.90659;  jer_param_er[9][5]=0.0207463;
	jer_param[9][6]=0.00138775;  jer_param_er[9][6]=0.000128634;
	jer_param[9][7]=-0.476516;  jer_param_er[9][7]=0.0607404;
	jer_param[9][8]=0.0153612;  jer_param_er[9][8]=0.000125725;
	jer_param[9][9]=-0.000159452;  jer_param_er[9][9]=3.31762e-05;
	jer_param[9][10]=-0.0590107;  jer_param_er[9][10]=0.0034089;
	jer_param[9][11]=0.000149436;  jer_param_er[9][11]=1.73205e-05;
	jer_param[9][12]=0.235095;  jer_param_er[9][12]=0.0112121;
	jer_param[9][13]=-0.000176939;  jer_param_er[9][13]=7.11021e-05;
	jer_param[9][14]=0.819459;  jer_param_er[9][14]=0.0137594;
	jer_param[9][15]=-0.0314006;  jer_param_er[9][15]=0.000374014;
	jer_param[9][16]=0.000308394;  jer_param_er[9][16]=9.04288e-06;
	jer_param[9][17]=0.325153;  jer_param_er[9][17]=0.0360778;
	jer_param[9][18]=0.000131873;  jer_param_er[9][18]=0.000206635;
	jer_param[10][0]=-0.953732;  jer_param_er[10][0]=0.00869531;
	jer_param[10][1]=0.0410471;  jer_param_er[10][1]=0.000196105;
	jer_param[10][2]=-0.000444301;  jer_param_er[10][2]=3.91545e-06;
	jer_param[10][3]=0.0307947;  jer_param_er[10][3]=0.00219417;
	jer_param[10][4]=-5.3727e-05;  jer_param_er[10][4]=9.27327e-06;
	jer_param[10][5]=1.06671;  jer_param_er[10][5]=0.0316977;
	jer_param[10][6]=0.00105544;  jer_param_er[10][6]=0.000168547;
	jer_param[10][7]=-0.268869;  jer_param_er[10][7]=0.00658315;
	jer_param[10][8]=0.00743125;  jer_param_er[10][8]=0.000153269;
	jer_param[10][9]=-9.31827e-05;  jer_param_er[10][9]=3.2684e-06;
	jer_param[10][10]=-0.0614147;  jer_param_er[10][10]=0.00327184;
	jer_param[10][11]=0.000104405;  jer_param_er[10][11]=1.53567e-05;
	jer_param[10][12]=0.268162;  jer_param_er[10][12]=0.0148949;
	jer_param[10][13]=-6.26511e-05;  jer_param_er[10][13]=8.37528e-05;
	jer_param[10][14]=-0.403697;  jer_param_er[10][14]=0.0180041;
	jer_param[10][15]=0.0298899;  jer_param_er[10][15]=0.000429542;
	jer_param[10][16]=-0.000414844;  jer_param_er[10][16]=9.54328e-06;
	jer_param[10][17]=0.28128;  jer_param_er[10][17]=0.0329415;
	jer_param[10][18]=0.000116748;  jer_param_er[10][18]=0.000171607;
	jer_param[11][0]=-0.721504;  jer_param_er[11][0]=0.0123063;
	jer_param[11][1]=0.0285797;  jer_param_er[11][1]=0.000254757;
	jer_param[11][2]=-0.000301013;  jer_param_er[11][2]=4.87484e-06;
	jer_param[11][3]=-0.0105836;  jer_param_er[11][3]=0.00269672;
	jer_param[11][4]=6.65018e-05;  jer_param_er[11][4]=1.06714e-05;
	jer_param[11][5]=1.15381;  jer_param_er[11][5]=0.0466366;
	jer_param[11][6]=0.000660203;  jer_param_er[11][6]=0.000219101;
	jer_param[11][7]=-0.610956;  jer_param_er[11][7]=0.00768133;
	jer_param[11][8]=0.0191861;  jer_param_er[11][8]=0.000172439;
	jer_param[11][9]=-0.000196723;  jer_param_er[11][9]=3.61353e-06;
	jer_param[11][10]=-0.0809709;  jer_param_er[11][10]=0.00406365;
	jer_param[11][11]=0.00014617;  jer_param_er[11][11]=1.73409e-05;
	jer_param[11][12]=0.28416;  jer_param_er[11][12]=0.0195007;
	jer_param[11][13]=-4.88533e-07;  jer_param_er[11][13]=0.000101491;
	jer_param[11][14]=-0.882516;  jer_param_er[11][14]=0.0179587;
	jer_param[11][15]=0.0451445;  jer_param_er[11][15]=0.000415322;
	jer_param[11][16]=-0.000517426;  jer_param_er[11][16]=9.05856e-06;
	jer_param[11][17]=0.196388;  jer_param_er[11][17]=0.0366708;
	jer_param[11][18]=0.00033582;  jer_param_er[11][18]=0.0001734;
	jer_param[12][0]=-0.0417805;  jer_param_er[12][0]=0.010778;
	jer_param[12][1]=0.00107649;  jer_param_er[12][1]=0.00021101;
	jer_param[12][2]=-6.06845e-05;  jer_param_er[12][2]=3.99126e-06;
	jer_param[12][3]=-0.0556253;  jer_param_er[12][3]=0.00408475;
	jer_param[12][4]=0.000160627;  jer_param_er[12][4]=1.62227e-05;
	jer_param[12][5]=1.01168;  jer_param_er[12][5]=0.0580327;
	jer_param[12][6]=0.00124107;  jer_param_er[12][6]=0.000266261;
	jer_param[12][7]=-1.81198;  jer_param_er[12][7]=0.00882027;
	jer_param[12][8]=0.0585738;  jer_param_er[12][8]=0.000165491;
	jer_param[12][9]=-0.000541543;  jer_param_er[12][9]=2.99558e-06;
	jer_param[12][10]=-0.124696;  jer_param_er[12][10]=0.006;
	jer_param[12][11]=0.000249786;  jer_param_er[12][11]=2.38289e-05;
	jer_param[12][12]=0.236166;  jer_param_er[12][12]=0.0233694;
	jer_param[12][13]=0.000272798;  jer_param_er[12][13]=0.000115376;
	jer_param[12][14]=-1.49496;  jer_param_er[12][14]=0.0268586;
	jer_param[12][15]=0.0732949;  jer_param_er[12][15]=0.000502876;
	jer_param[12][16]=-0.000831473;  jer_param_er[12][16]=9.13321e-06;
	jer_param[12][17]=0.289883;  jer_param_er[12][17]=0.0557916;
	jer_param[12][18]=-0.000132964;  jer_param_er[12][18]=0.000224823;
	jer_param[13][0]=-2.94925;  jer_param_er[13][0]=0.0262199;
	jer_param[13][1]=0.085559;  jer_param_er[13][1]=0.000430568;
	jer_param[13][2]=-0.000628295;  jer_param_er[13][2]=6.9836e-06;
	jer_param[13][3]=-0.0540865;  jer_param_er[13][3]=0.0136904;
	jer_param[13][4]=0.000139313;  jer_param_er[13][4]=4.936e-05;
	jer_param[13][5]=1.16241;  jer_param_er[13][5]=0.154151;
	jer_param[13][6]=0.00449765;  jer_param_er[13][6]=0.000680983;
	jer_param[13][7]=-0.909461;  jer_param_er[13][7]=0.0136892;
	jer_param[13][8]=0.0269072;  jer_param_er[13][8]=0.00022884;
	jer_param[13][9]=-0.000263558;  jer_param_er[13][9]=3.78923e-06;
	jer_param[13][10]=-0.120874;  jer_param_er[13][10]=0.00759237;
	jer_param[13][11]=0.000149562;  jer_param_er[13][11]=2.68079e-05;
	jer_param[13][12]=0.401773;  jer_param_er[13][12]=0.0401819;
	jer_param[13][13]=0.00048766;  jer_param_er[13][13]=0.000171866;
	jer_param[13][14]=-2.45903;  jer_param_er[13][14]=0.0337559;
	jer_param[13][15]=0.080548;  jer_param_er[13][15]=0.000562604;
	jer_param[13][16]=-0.000635591;  jer_param_er[13][16]=9.28515e-06;
	jer_param[13][17]=0.123775;  jer_param_er[13][17]=0.0329547;
	jer_param[13][18]=-6.06731e-05;  jer_param_er[13][18]=0.00011711;
	jer_param[14][0]=0;  jer_param_er[14][0]=0;
	jer_param[14][1]=0;  jer_param_er[14][1]=0;
	jer_param[14][2]=0;  jer_param_er[14][2]=0;
	jer_param[14][3]=0;  jer_param_er[14][3]=0;
	jer_param[14][4]=0;  jer_param_er[14][4]=0;
	jer_param[14][5]=0;  jer_param_er[14][5]=0;
	jer_param[14][6]=0;  jer_param_er[14][6]=0;
	jer_param[14][7]=-430.789;  jer_param_er[14][7]=0.377943;
	jer_param[14][8]=8.10679;  jer_param_er[14][8]=0.00358588;
	jer_param[14][9]=-0.0382566;  jer_param_er[14][9]=3.38665e-05;
	jer_param[14][10]=-0.143086;  jer_param_er[14][10]=0.00730127;
	jer_param[14][11]=6.88638e-05;  jer_param_er[14][11]=1.81481e-05;
	jer_param[14][12]=0.9097;  jer_param_er[14][12]=0.0515384;
	jer_param[14][13]=-0.000149165;  jer_param_er[14][13]=0.000158924;
	jer_param[14][14]=0;  jer_param_er[14][14]=0;
	jer_param[14][15]=0;  jer_param_er[14][15]=0;
	jer_param[14][16]=0;  jer_param_er[14][16]=0;
	jer_param[14][17]=0;  jer_param_er[14][17]=0;
	jer_param[14][18]=0;  jer_param_er[14][18]=0;

	if(i<15 && j<19)
	{
		if(code==0) value=jer_param[i][j];
		if(code==1) value=jer_param_er[i][j];
	}
	return value;
}
// //-------------- fit function: poly3 * step func + poly2 * (step func)
// this is the new fit function used for 
// Gaus Mean, Landau MPV, and Gaus Norm
double dJayFunc1_metcode(double *x, double *p)
{
		double f1 = p[0]  + p[1] * x[0] + p[2] * x[0] * x[0];
		double f2 = p[3]  + p[4] * x[0];

		double ws = 1 / (1 + exp((x[0]-60.0)/10.5));
		double fitval = f1 * ws + f2 * (1 - ws);
		return fitval;
}
void METModelCodeTest()
{
	gStyle->SetAxisColor(kBlue);
	gStyle->SetPadTickY(1);
	gStyle->SetPadTickX(1);
	gStyle->SetOptFit(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetTitleFontSize(0.035);
	gStyle->SetStatFontSize(0.03);
	gStyle->SetTitleFont(2);

	for (int eta_bin = 0; eta_bin <=14 ; ++eta_bin)
	{
		TF1 *fit;
		std::cout << "fit = " << fit << std::endl;
		std::vector<TF1*> funcList;
		funcList = JetResFitPerEtaBin(eta_bin);
		assert(funcList.size() == 5 && "Did not get 5 functions from original fits");
		std::cout << "func size = " << funcList.size() << std::endl;

		TCanvas* c1 = new TCanvas("refitcanvas","Jet energy resolution study: Testing;Jet E(GeV);",200,10,1100,850);
		c1->Divide(2,3);

		TLegend *leg = new TLegend (0.7,0.7,0.90,0.90);
		//leg->SetTextFont(42);
		//leg->SetTextSize(0.025);
		leg->SetBorderSize (1);
		leg->SetFillColor (10);

		float etamin = eta_bin * 0.2;
		float etamax = etamin + 0.2;
		std::stringstream rattitle;
		if (eta_bin<14) rattitle << "Compare new and old fit functions for " << etamin << "<|#Eta|<" <<etamax << ";Jet E(GeV); New Fit/Old Fit - 1";
		else rattitle << "Compare new and old fit functions for | #Eta |>" <<etamin << ";Jet E(GeV); New Fit/Old Fit - 1";
		TH1F *ratiohist;
		//TH1F *ratiohist = new TH1F("overlay",rattitle.str().c_str(),1000,0,500);
		TCanvas *c2 = new TCanvas("oldnewcanvas","Jet energy resolution study: NewOld fit functions",200,10,1100,850);
		c2->Divide(2,3);

		std::string title;
	
		for (int iFunc=0; iFunc<=4; ++iFunc)
		{
			double p0=0,p1=0,p2=0,p3=0,p4=0;
			std::cout << "etabin, iFunc = " << eta_bin << ", " << iFunc << std::endl;

			if (iFunc == 0)  //pick up Gaus Mean parameters
			{
				p0=JERparam(0,eta_bin,0);
				p1=JERparam(0,eta_bin,1);
				p2=JERparam(0,eta_bin,2);
				p3=JERparam(0,eta_bin,3);
				p4=JERparam(0,eta_bin,4);
			} else if (iFunc == 1)  //pick up Gaus Sigma parameterization
			{
				p0=JERparam(0,eta_bin,5);
				p1=JERparam(0,eta_bin,6);
			} else if (iFunc == 2)  //pick up Landau MPV parameterization
			{
				p0=JERparam(0,eta_bin,7);
				p1=JERparam(0,eta_bin,8);
				p2=JERparam(0,eta_bin,9);
				p3=JERparam(0,eta_bin,10);
				p4=JERparam(0,eta_bin,11);
			} else if (iFunc == 3)  //pick up Landau Sigma parameterization
			{
				p0=JERparam(0,eta_bin,12);
				p1=JERparam(0,eta_bin,13);
			} else if (iFunc == 4)	//pick up Gaus Norm Parameterization
			{
				p0=JERparam(0,eta_bin,14);
				p1=JERparam(0,eta_bin,15);
				p2=JERparam(0,eta_bin,16);
				p3=JERparam(0,eta_bin,17);
				p4=JERparam(0,eta_bin,18);
			}

			if (iFunc==0) title = "Gaus Mean";
			else if (iFunc==1) title = "Gaus Sigma";
			else if (iFunc==2) title = "Landau MPV";
			else if (iFunc==3) title = "LandauSigma";
			else if (iFunc==4) title = "Gaus Norm";

			std::stringstream funcname;
			if (iFunc == 0 || iFunc == 2 || iFunc == 4)
			{
				funcname << "JayFitFunc3_" << iFunc << std::endl;
				fit = new TF1(funcname.str().c_str(),dJayFunc1_metcode,5,500,5);
				fit->SetParameter(0,p0);
				fit->SetParameter(1,p1);
				fit->SetParameter(2,p2);
				fit->SetParameter(3,p3);
				fit->SetParameter(4,p4);
			} else {
				funcname << "fit_4_" << iFunc << std::endl;
				fit = new TF1(funcname.str().c_str() ,fitFunc3,2.0,500.0,2);
				fit->SetParameter(0,p0);
				fit->SetParameter(1,p1);
			}

			fit->SetTitle(title.c_str());
			fit->SetLineColor(kRed);
			fit->SetLineWidth(3);

			funcList.at(iFunc)->SetTitle(title.c_str());
			funcList.at(iFunc)->SetLineWidth(5);

			//float ymax = max(fit->GetMaximum(),funcList.at(iFunc)->GetMaximum());
			//float ymin = min(fit->GetMinimum(),funcList.at(iFunc)->GetMinimum());

			c1->cd(iFunc+1);
			funcList.at(iFunc)->Draw();
			fit->SetLineStyle(2);
			fit->Draw("same");
			if (iFunc == 0)
			{
				leg->AddEntry(fit, "METModel code");
				leg->AddEntry(funcList.at(iFunc),"JER Fit Code");
			}
			leg->Draw();

			//to make the ratio plots of the Sasha's fits and mine
			//funcList 0,2,4 has hist f1s
			
			if (iFunc == 0 || iFunc == 2 || iFunc == 4)
			{

				ratiohist = new TH1F("overlay",rattitle.str().c_str(),1000,0,500);
				ratiohist->SetMinimum(-1);
				ratiohist->SetMaximum(1);
				ratiohist->SetLineColor(20);

				TLegend *leg2 = new TLegend (0.7,0.7,0.90,0.90);
				leg2->SetBorderSize (1);
				leg2->SetFillColor (10);

				//gStyle->UseCurrentStyle();
				gStyle->SetOptStat(0);
				c2->cd(iFunc+1);
				funcList.at(iFunc)->Draw();
				fit->Draw("same");

				leg2->AddEntry(fit, "NEW Fit");
				leg2->AddEntry(funcList.at(iFunc),"OLD fit");
				leg2->Draw();

				for (int i=1; i<= ratiohist->GetNbinsX(); ++i)
				{
					if (i%2!=0) continue;
					double oldval = funcList.at(iFunc)->Eval(ratiohist->GetBinCenter(i));
					double newval = fit->Eval(ratiohist->GetBinCenter(i));
					float rat = newval/oldval -1;
					ratiohist->SetBinContent(i,rat);
				}
				c2->cd(iFunc+2);
				ratiohist->Draw();

			}
			
		} // for each function

		c2->cd();
		std::stringstream giffile2, epsfile2;
		giffile2 << "~/TMP/JERoldnewfit_EtaBin"<< eta_bin << "_" << title <<  ".gif";
		epsfile2 << "~/TMP/JERoldnewfit_EtaBin"<< eta_bin << "_" << title <<  ".eps";
		//c2->Print(giffile2.str().c_str());
		c2->Print(epsfile2.str().c_str(),"Preview");

		
		c1->cd();
		
		std::stringstream giffile, epsfile;
		giffile << "~/TMP/JERfit_EtaBin"<< eta_bin <<  ".gif";
		epsfile << "~/TMP/JERfit_EtaBin"<< eta_bin <<  ".eps";
		//c1->Print(giffile.str().c_str());
		//c1->Print(epsfile.str().c_str(),"Preview");
		c1->Close();

	} //for each eta bin

}
