/////////////////////////////////////////////////////////
// This is Sasha's JET energy resolution fit method.   //
// You can derive JER parametirization from the fits   //
// to be put in the JetFilter module to model the JER. //
// See elog#1227.                                      //
// I renamed this macro from myJetResFit_fnl_default.C //
// This fits for all eta bins.                         //
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
#include <cmath>

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
  double fitval=0.0;
  double arg=-x[0]/(x[0]+1.0);
  double arg_L=-x[0]/(x[0]+1.0);
  double mean=par[1];
  double sigmaG=par[2];
  double mpv=par[3];
  double sigmaL=par[4];
  double normL=par[5];
  if(normL<0.0) return -1.0E6;
  double f1=normL/(1.0+normL)*TMath::Gaus(arg,mean,sigmaG);
  double f2=TMath::Landau(arg_L,mpv,sigmaL)/(1.0+normL);
  fitval=par[0]*(f1+f2);
  return fitval;
}

//------ Gauss mean and Landau MPV: A0+A1*Pt+A3/Pt ---- default
double fitFunc0(double *x, double *par)
{
  double fitval=0.0;
  if(x[0]>0.0) fitval=par[0]+par[1]*x[0]+par[2]/x[0];
  else fitval=-1.0E6;
  return fitval;
}
//------ Gauss norm: [A0+A1*sqrt(Pt)]/Pt+A3 ---- default
double fitFunc5(double *x, double *par)
{
  double fitval=0.0;
  if(x[0]>0.0) fitval=(par[0]+par[1]*sqrt(x[0]))/x[0]+par[2];
  else fitval=-1.0E6;
  return fitval;
}
// //----- Landau Sigma fit function: sqrt(A0/Pt+A2) ---- default
double fitFunc3(double *x, double *par)
{
  double fitval=0.0;
  if(x[0]>0.0)
    {
      double arg=0.0;
      arg=par[0]/x[0]+par[1];
      fitval=sqrt(arg);
    }
  else fitval=-1.0E6;
  return fitval;
}

//=====================================================================
//-------------- Driving script as of 08/02/07
// --- fits to all fJetRes_TrueBal for all eta_det bins 
//
// Fit function is [Const*Gaus(-x/(1+x))+Landau(-x/(1+x))]/(1+Const)
//_____________________________________________________________________
void JetResFitAllEtaBins(const bool bNewfits = 1, int folder_id=0)
{

	if (! (folder_id == 0 || folder_id == 3 || folder_id == 4))
	{
		std::cout << "Invalid folder selection. must be 1,3 or 4. Returning" << std::endl;
		return;
	}


	gROOT->Reset();

	char foldername[200];
	char fObject[200];

	if(folder_id==0) sprintf(foldername,"JetRes_def");
	//these folders no longer exists
	//if(folder_id==1) sprintf(foldername,"JetRes_syst1");
	//if(folder_id==2) sprintf(foldername,"JetRes_syst2");
	if(folder_id==3) sprintf(foldername,"JetRes_def_nvx1");
	if(folder_id==4) sprintf(foldername,"JetRes_def_nvx2");

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
	
	//old fit parameters
	double fnl_fit_param[15][13];
	double fnl_fit_paramER[15][13];
	double corr_param[15][11];//???????
	double binvals[15][5][60];		//eta bin, 5 parameters/60 bin values (1GeV bins)


	
	//new fit parameters
	double fnl_newfit_param[15][31];
	double fnl_newfit_paramER[15][31];



	if (bNewfits ==1 )
	{
		for(int i=0; i<15; i++)
		{
			for(int j=0; j<31; j++)
			{
				fnl_newfit_param[i][j]=0.0;
				fnl_newfit_paramER[i][j]=0.0;
			}
		}
	} else {

		for(int i=0; i<15; i++)
		{
			for(int j=0; j<13; j++)
			{
				fnl_fit_param[i][j]=0.0;
				fnl_fit_paramER[i][j]=0.0;
				if(j<11) corr_param[i][j]=0.0;
			}

				
			for (int par=0; par <5; ++par)
			{
				for (int bin=0;bin<60;++bin)
				{
					binvals[i][par][bin] = 0;
				}
			}

			
		}

		
	}
	//________________________________________________________________________________
	//------------------- stopped here -----------------------------------------------
	//________________________________________________________________________________

	TH1F *h1[100];

	gDirectory->pwd();

	//-------------- matching in dR<0.1
	TFile myFile1("JER_hadE_sample1_dR01_072507.root","r");


	TFolder* myFolder1 = (TFolder* ) myFile1.FindObjectAny("Ana");

	for (int eta_bin_id = 0; eta_bin_id < 15; eta_bin_id++)
	{
		sprintf(fObject,"MyJetFilter/Hist/%s/fJetRes_AveJetE[%i]",foldername,eta_bin_id); // get object name
		TProfile *_avePt_histo= (TProfile*) myFolder1->FindObjectAny(fObject);

		for(int i=0; i<100; i++)
		{
			sprintf(fObject,"MyJetFilter/Hist/%s/fJetRes_TrueBal[%i][%i]",foldername,eta_bin_id,i); // get object name
			h1[i]= (TH1F*) myFolder1->FindObjectAny(fObject);
			avePt[i]=_avePt_histo->GetBinContent(i+1);
			avePt_er[i]=_avePt_histo->GetBinError(i+1);
			histo_Nentry[i]=h1[i]->GetEntries();
			if(histo_Nentry[i]>99) h_bin_last=i;
			gauss_mean[i]=0.0;
			gauss_mean_er[i]=0.0;
			gauss_sigma[i]=0.0;
			gauss_sigma_er[i]=0.0;
			landau_mean[i]=0.0;
			landau_mean_er[i]=0.0;
			landau_sigma[i]=0.0;
			landau_sigma_er[i]=0.0;
			landau_norm[i]=0.0;
			landau_norm_er[i]=0.0;
			fit_chi2[i]=0.0;
			fit_chi2_er[i]=0.0; 
			fit_status[i]=0;

			maxFit_X[i]=0.0;
			histoRMS[i]=0.01;
		}

		for(int i=0; i<100; i++)
		{
			std::cout<<"---------------------------------------------"<<std::endl;
			std::cout<<"     Fitting bin="<<i<<std::endl;
			std::cout<<"_____________________________________________"<<std::endl;

			if(histo_Nentry[i]>=10000.0) h1[i]->Rebin(2);
			if(histo_Nentry[i]>=1000.0 && histo_Nentry[i]<10000.0) h1[i]->Rebin(3);
			if(histo_Nentry[i]>=200.0 && histo_Nentry[i]<1000.0) h1[i]->Rebin(3);
			if(histo_Nentry[i]>=100.0 && histo_Nentry[i]<200.0) h1[i]->Rebin(4);
			if(eta_bin_id<14 && histo_Nentry[i]<100.0) continue; // skipping bins with Nevnt<100
			if(eta_bin_id==14 && histo_Nentry[i]<80.0) continue; // skipping bins with Nevnt<100
			if(histo_Nentry[i]<100.0) h1[i]->Rebin(4);

			int nbins=h1[i]->GetNbinsX();
			float binwidth=h1[i]->GetBinWidth(1);
			float lowedge=(h1[i]->GetBinCenter(1))-binwidth/2.0;
			float highedge=lowedge+nbins*binwidth;

			char title[200];
			char name[200];

			sprintf(name,"bal");
			if(eta_bin_id<14) sprintf(title,"E^{had}/E^{det}-1.0, %2.1f<|#eta_{det}^{jet}|<%2.1f and %i<E^{jet}<%i",0.2*eta_bin_id,0.2*(eta_bin_id+1),5*i,5*(i+1));
			else sprintf(title,"E^{had}/E^{det}-1.0, |#eta_{det}^{jet}|>%2.1f and %i<E^{jet}<%i",0.2*eta_bin_id,5*i,5*(i+1));
			TH1F* _h= new TH1F(name,title,nbins,lowedge,highedge);
			_h->Sumw2();

			float cut1=-1.0;
			float cut2=2.0;
			float val=0.0;
			float err=0.0;
			int bin_first=0;
			int bin_last=0;

			for(int j=0; j<nbins; j++)
			{
				val=h1[i]->GetBinContent(j+1);
				err=h1[i]->GetBinError(j+1);
				_h->SetBinContent(j+1,val);
				_h->SetBinError(j+1,err);
				if(bin_first==0 && val>0.0) bin_first=j+1;
				if(val>0.0 && h1[i]->GetBinCenter(j+1)<3.5) bin_last=j+1; 
			}
			cut1=h1[i]->GetBinCenter(bin_first);
			cut2=h1[i]->GetBinCenter(bin_last);

			int maxbin=h1[i]->GetMaximumBin();
			double hight=h1[i]->GetBinContent(maxbin);
			double mean_histo=h1[i]->GetMean();
			double rms_histo=h1[i]->GetRMS();

			TF1 *fitFun1=new TF1("fit1",myFunc1,cut1-binwidth,cut2+binwidth,6);
			if(eta_bin_id<14)
			{
				fitFun1->SetParameter(0,0.9*hight); // norm 
				fitFun1->SetParameter(1,mean_histo); // mean of Gauss
				fitFun1->SetParameter(2,0.7*rms_histo); // Sigma of Gauss
				fitFun1->SetParameter(3,mean_histo); // MPV of Landau
				fitFun1->SetParameter(4,0.5*rms_histo); // Sigma of Landau
				fitFun1->SetParameter(5,0.2); // Landau norm
			}
			else
			{
				fitFun1->SetParameter(0,0.9*hight); // norm 
				fitFun1->FixParameter(1,0.0); // mean of Gauss
				fitFun1->FixParameter(2,0.0); // Sigma of Gauss
				fitFun1->SetParameter(3,mean_histo); // MPV of Landau
				fitFun1->SetParameter(4,0.5*rms_histo); // Sigma of Landau
				fitFun1->FixParameter(5,0.0); // Landau norm
			}
			//       fitFun1->SetParLimits(0,0.0,100.0*hight);      
			//       fitFun1->SetParLimits(1,mean_histo-2.0*rms_histo,mean_histo+2.0*rms_histo);
			//       fitFun1->SetParLimits(2,0.0,3.0*rms_histo);
			//       fitFun1->SetParLimits(3,mean_histo-2.0*rms_histo,mean_histo+2.0*rms_histo);      
			//       fitFun1->SetParLimits(4,0.0,3.0*rms_histo);      
			//       fitFun1->SetParLimits(5,0.0,1000.0*hight);     
			fitFun1->SetParName(0,"Gauss: const");
			fitFun1->SetParName(1,"Gauss: mean");
			fitFun1->SetParName(2,"Gauss: #sigma");
			fitFun1->SetParName(3,"Landau: MPV");
			fitFun1->SetParName(4,"Landau: #sigma");
			fitFun1->SetParName(5,"Gaus: const");

			_h->Fit(fitFun1,"BQ","",cut1,cut2);

			//---------- re-fitting second time --------------------------
			if(eta_bin_id<14)
			{
				fitFun1->SetParameter(0,fitFun1->GetParameter(0)); // norm 
				fitFun1->SetParameter(1,fitFun1->GetParameter(1)); // mean of Gauss
				fitFun1->SetParameter(2,fitFun1->GetParameter(2)); // Sigma of Gauss
				fitFun1->SetParameter(3,fitFun1->GetParameter(3)); // MPV of Landau
				fitFun1->SetParameter(4,fitFun1->GetParameter(4)); // Sigma of Landau
				fitFun1->SetParameter(5,fitFun1->GetParameter(5)); // Landau norm
			}
			else
			{
				fitFun1->SetParameter(0,fitFun1->GetParameter(0)); // norm 
				fitFun1->FixParameter(1,0.0); // mean of Gauss
				fitFun1->FixParameter(2,0.0); // Sigma of Gauss
				fitFun1->SetParameter(3,fitFun1->GetParameter(3)); // MPV of Landau
				fitFun1->SetParameter(4,fitFun1->GetParameter(4)); // Sigma of Landau
				fitFun1->FixParameter(5,0.0); // Landau norm	  
			}

			//       _h->Fit(fitFun1,"LEB","",cut1,cut2);
			_h->Fit(fitFun1,"EBQ","",cut1,cut2);
			//---------- end of re-fitting -------------------------------

			if((gMinuit->fCstatu).Contains("SUCCESSFUL")==true) fit_status[i]=3;
			if((gMinuit->fCstatu).Contains("PROBLEMS")==true) fit_status[i]=2;
			if((gMinuit->fCstatu).Contains("FAILURE")==true) fit_status[i]=0;
			std::cout<<" Minuit Error Matrix Status ="<<gMinuit->fCstatu<<std::endl;
			std::cout<<" Minuit Error Matrix Status ="<<fit_status[i]<<std::endl;
			//       gMinuit->mnmatu(1); 
			maxFit_X[i]=fitFun1->GetMaximumX(-1.0,1.0);
			histoRMS[i]=rms_histo;
			gauss_mean[i]=fitFun1->GetParameter(1);
			gauss_mean_er[i]=fitFun1->GetParError(1);
			gauss_sigma[i]=fitFun1->GetParameter(2);
			gauss_sigma_er[i]=fitFun1->GetParError(2);
			landau_mean[i]=fitFun1->GetParameter(3);
			landau_mean_er[i]=fitFun1->GetParError(3);
			landau_sigma[i]=fitFun1->GetParameter(4);
			landau_sigma_er[i]=fitFun1->GetParError(4);
			landau_norm[i]=fitFun1->GetParameter(5);
			landau_norm_er[i]=fitFun1->GetParError(5);

			if(fitFun1->GetNDF()>0)
				fit_chi2[i]=(1.0*fitFun1->GetChisquare())/(1.0*fitFun1->GetNDF());
			else fit_chi2[i]=-1.0;

			delete fitFun1;
			delete _h;
		}

		float cut1=avePt[h_bin_first]-2.0;
		float cut2=avePt[h_bin_last]+2.0;

		char title[200];
		char name[200]; 

		//--------------------- Gauss mean ---------------------------------------------
		sprintf(name,"h1");
		if(eta_bin_id<14) sprintf(title,"Energy dependence of Gauss mean for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
		else sprintf(title,"Energy dependence of Gauss mean for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
		TH1F* h1_fnl= new TH1F(name,title,500,0.0,500.0);
		h1_fnl->Sumw2();
		//--------------------- Gauss sigma ---------------------------------------------
		sprintf(name,"h2");
		if(eta_bin_id<14) sprintf(title,"Energy dependence of Gauss #sigma for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
		else sprintf(title,"Energy dependence of Gauss #sigma for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
		TH1F* h2_fnl= new TH1F(name,title,500,0.0,500.0);
		h2_fnl->Sumw2();
		//--------------------- Landau MPV ---------------------------------------------
		sprintf(name,"h3");
		if(eta_bin_id<14) sprintf(title,"Energy dependence of Landau MPV for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
		else sprintf(title,"Energy dependence of Landau MPV for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
		TH1F* h3_fnl= new TH1F(name,title,500,0.0,500.0);
		h3_fnl->Sumw2();
		//--------------------- Landau sigma ---------------------------------------------
		sprintf(name,"h4");
		if(eta_bin_id<14) sprintf(title,"Energy dependence of Landau #sigma for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
		else sprintf(title,"Energy dependence of Landau #sigma for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
		TH1F* h4_fnl= new TH1F(name,title,500,0.0,500.0);
		h4_fnl->Sumw2();
		//--------------------- Gauss normalization ---------------------------------------------
		sprintf(name,"h5");
		if(eta_bin_id<14) sprintf(title,"Energy dependence of Gauss normalization for %2.1f<|#eta_{det}^{jet}|<%2.1f",0.2*eta_bin_id,0.2*(eta_bin_id+1));
		else sprintf(title,"Energy dependence of Gauss normalization for |#eta_{det}^{jet}|>%2.1f",0.2*eta_bin_id);
		TH1F* h5_fnl= new TH1F(name,title,500,0.0,500.0);
		h5_fnl->Sumw2();
		//--------------------- Filling up histograms ---------------------------------------------
		int last_E_bin=0;
		int first_E_bin=100;


		double fit_cut1=3.0;
		double fit_cut2=500.0;
		if(eta_bin_id<5) fit_cut1=12.0;
		if(eta_bin_id>4 && eta_bin_id<8) fit_cut1=17.0;
		if(eta_bin_id==8) fit_cut1=22.0;
		if(eta_bin_id==9) fit_cut1=27.0;
		if(eta_bin_id==10) fit_cut1=32.0;
		if(eta_bin_id==11) fit_cut1=37.0;
		if(eta_bin_id==12) fit_cut1=47.0;
		if(eta_bin_id==13) fit_cut1=57.0;
		if(eta_bin_id==14) fit_cut1=102.0;

		for(int i=0; i<100; i++)
		{
			bool good_fit=true;

			if(eta_bin_id==2 && gauss_mean[i]>0.1) good_fit=false;
			if(eta_bin_id<10 && gauss_sigma[i]<0.05) good_fit=false;
			if(eta_bin_id==7 && (landau_norm[i]>0.4 || gauss_mean[i]<-0.02)) good_fit=false; 
			if(eta_bin_id>8 && eta_bin_id<13 && landau_mean[i]>0.0) good_fit=false; 
			if(eta_bin_id<13)
			{	  
				if(fabs(landau_mean[i]/histoRMS[i])>1.1) good_fit=false;
				if(gauss_sigma[i]>1.0*histoRMS[i] || gauss_sigma[i]<0.2*histoRMS[i]) good_fit=false;
				if(landau_sigma[i]>0.9*histoRMS[i] || landau_sigma[i]<0.1*histoRMS[i]) good_fit=false;
				if(gauss_sigma_er[i]<0.01*gauss_sigma[i] || gauss_sigma_er[i]>0.3*gauss_sigma[i]) good_fit=false;
				if(landau_sigma_er[i]<0.01*landau_sigma[i] || landau_sigma_er[i]>0.3*landau_sigma[i]) good_fit=false;
				if(landau_mean_er[i]>0.06) good_fit=false;
				if(gauss_mean_er[i]>0.15) good_fit=false;
				if(landau_norm_er[i]<0.05*landau_norm[i]) good_fit=false;
			}
			if(eta_bin_id==13)
			{
				if(landau_norm[i]>0.2) good_fit=false;
				if(gauss_sigma[i]<0.06) good_fit=false;
				if(gauss_mean[i]<-0.1) good_fit=false;
			}
			if(eta_bin_id<14 && histo_Nentry[i]<100) good_fit=false;
			if(eta_bin_id==14 && (histo_Nentry[i]<80 || landau_mean[i]<-0.17)) good_fit=false;
			if(fit_status[i]!=3) good_fit=false;
			if(avePt[i]<fit_cut1) good_fit=false;

		//this is a condition for my last test to use exact bin values for <60GeV
		//and use a fit to >60 region.
			if (i<=12) good_fit = true;  //include all points below 60GeV temporily

			if(good_fit==true)
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

				last_E_bin=i;
				fit_cut2=avePt[last_E_bin]+1.0;
				if(i<first_E_bin) 
				{
					first_E_bin=i;
					if(fit_cut1<avePt[first_E_bin]-2.0) fit_cut1=avePt[first_E_bin]-1.0;
				}
			}
		}


	//temporary
	fit_cut1 = 60;

		
		//________ parameterization for Gauss mean __________________
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

			if(eta_bin_id<14)
			{
				fnl_newfit_param[eta_bin_id][0]=JayFitFunc1->GetParameter(0);
				fnl_newfit_paramER[eta_bin_id][0]=JayFitFunc1->GetParError(0);
				fnl_newfit_param[eta_bin_id][1]=JayFitFunc1->GetParameter(1);
				fnl_newfit_paramER[eta_bin_id][1]=JayFitFunc1->GetParError(1);
				fnl_newfit_param[eta_bin_id][2]=JayFitFunc1->GetParameter(2);
				fnl_newfit_paramER[eta_bin_id][2]=JayFitFunc1->GetParError(2);
				fnl_newfit_param[eta_bin_id][3]=JayFitFunc1->GetParameter(3);
				fnl_newfit_paramER[eta_bin_id][3]=JayFitFunc1->GetParError(3);
				fnl_newfit_param[eta_bin_id][4]=JayFitFunc1->GetParameter(4);
				fnl_newfit_paramER[eta_bin_id][4]=JayFitFunc1->GetParError(4);
				fnl_newfit_param[eta_bin_id][5]=JayFitFunc1->GetParameter(5);
				fnl_newfit_paramER[eta_bin_id][5]=JayFitFunc1->GetParError(5);
				fnl_newfit_param[eta_bin_id][6]=JayFitFunc1->GetParameter(6);
				fnl_newfit_paramER[eta_bin_id][6]=JayFitFunc1->GetParError(6);
			} else {
				fnl_newfit_param[eta_bin_id][0]=0.0;
				fnl_newfit_paramER[eta_bin_id][0]=0.0;
				fnl_newfit_param[eta_bin_id][1]=0.0;
				fnl_newfit_paramER[eta_bin_id][1]=0.0;
				fnl_newfit_param[eta_bin_id][2]=0.0;
				fnl_newfit_paramER[eta_bin_id][2]=0.0;	  
				fnl_newfit_param[eta_bin_id][3]=0.0;
				fnl_newfit_paramER[eta_bin_id][3]=0.0;	  
				fnl_newfit_param[eta_bin_id][4]=0.0;
				fnl_newfit_paramER[eta_bin_id][4]=0.0;	  
				fnl_newfit_param[eta_bin_id][5]=0.0;
				fnl_newfit_paramER[eta_bin_id][5]=0.0;	  
				fnl_newfit_param[eta_bin_id][6]=0.0;
				fnl_newfit_paramER[eta_bin_id][6]=0.0;	  
			}

			delete JayFitFunc1;

		} else {

			TF1 *param_fitFun1=new TF1("fit_1",fitFunc0,2.0,500.0,3);
			param_fitFun1->SetParameter(0,0.09); // A0 param: term=A0
			param_fitFun1->SetParameter(1,-0.0002); // param A1: term=A1*Pt
			param_fitFun1->SetParameter(2,-1.0); // param A2: term=A2/Pt
			h1_fnl->Fit(param_fitFun1,"EB","",fit_cut1,fit_cut2);

			gMinuit->mnmatu(1);
			if(TMath::Sqrt(TMath::Abs((gMinuit->fVhmat[2])*((gMinuit->fVhmat[0]))))>0.0) 
				corr_param[eta_bin_id][0]=(gMinuit->fVhmat[1])/TMath::Sqrt(TMath::Abs((gMinuit->fVhmat[2])*((gMinuit->fVhmat[0]))));
			corr_param[eta_bin_id][1]=gMinuit->fMATUvline[0];
			corr_param[eta_bin_id][2]=gMinuit->fMATUvline[1];

			if(eta_bin_id<14)
			{
				fnl_fit_param[eta_bin_id][0]=param_fitFun1->GetParameter(0);
				fnl_fit_paramER[eta_bin_id][0]=param_fitFun1->GetParError(0);
				fnl_fit_param[eta_bin_id][1]=param_fitFun1->GetParameter(1);
				fnl_fit_paramER[eta_bin_id][1]=param_fitFun1->GetParError(1);
				fnl_fit_param[eta_bin_id][2]=param_fitFun1->GetParameter(2);
				fnl_fit_paramER[eta_bin_id][2]=param_fitFun1->GetParError(2);
			} else {
				fnl_fit_param[eta_bin_id][0]=0.0;
				fnl_fit_paramER[eta_bin_id][0]=0.0;
				fnl_fit_param[eta_bin_id][1]=0.0;
				fnl_fit_paramER[eta_bin_id][1]=0.0;
				fnl_fit_param[eta_bin_id][2]=0.0;
				fnl_fit_paramER[eta_bin_id][2]=0.0;	  
			}

			for (int i=1; i<=h1_fnl->GetNbinsX(); ++i)
			{
				if (h1_fnl->GetBinContent(i))
					std::cout << "bin [" << i << ", " << h1_fnl->GetBinContent(i) << "]" << std::endl;
			}
				
			for (int bin=1;bin <=60; ++bin)
			{
				binvals[eta_bin_id][0][bin-1] = h1_fnl->GetBinContent(bin);
				std::cout << "binvals[" << eta_bin_id << "][0]["<< bin-1 
					<< "] = " << binvals[eta_bin_id][0][bin-1] 
					<< "     --- Binedges = " 
					<<  h1_fnl->GetBinLowEdge(bin) << ", " 
					<< h1_fnl->GetXaxis()->GetBinUpEdge(bin) 
					<< std::endl;
			}

			delete param_fitFun1;
		}

		//________ parameterization for Gauss sigma __________________

		if (bNewfits)
		{
			TF1 *stitchTest2 = new TF1("stitchTest2",fitFunc3_new,2,500,5);
			stitchTest2->SetParameter(0,1.6);
			stitchTest2->SetParameter(1,0.00005);
			stitchTest2->SetParameter(2,2); 
			stitchTest2->SetParameter(3,0.004);
			stitchTest2->SetParameter(4,30);
			stitchTest2->SetParLimits(4,20,40.0);
			h2_fnl->Fit(stitchTest2,"EBQ","",fit_cut1,fit_cut2);

			if(eta_bin_id<14)
			{
				fnl_newfit_param[eta_bin_id][7]=stitchTest2->GetParameter(0);
				fnl_newfit_paramER[eta_bin_id][7]=stitchTest2->GetParError(0);
				fnl_newfit_param[eta_bin_id][8]=stitchTest2->GetParameter(1);
				fnl_newfit_paramER[eta_bin_id][8]=stitchTest2->GetParError(1);
				fnl_newfit_param[eta_bin_id][9]=stitchTest2->GetParameter(2);
				fnl_newfit_paramER[eta_bin_id][9]=stitchTest2->GetParError(2);
				fnl_newfit_param[eta_bin_id][10]=stitchTest2->GetParameter(3);
				fnl_newfit_paramER[eta_bin_id][10]=stitchTest2->GetParError(3);
				fnl_newfit_param[eta_bin_id][11]=stitchTest2->GetParameter(4);
				fnl_newfit_paramER[eta_bin_id][11]=stitchTest2->GetParError(4);
			} else {
				fnl_newfit_param[eta_bin_id][7]=0.0;
				fnl_newfit_paramER[eta_bin_id][7]=0.0;
				fnl_newfit_param[eta_bin_id][8]=0.0;
				fnl_newfit_paramER[eta_bin_id][8]=0.0;
				fnl_newfit_param[eta_bin_id][9]=0.0;
				fnl_newfit_paramER[eta_bin_id][9]=0.0;	  
				fnl_newfit_param[eta_bin_id][10]=0.0;
				fnl_newfit_paramER[eta_bin_id][10]=0.0;	  
				fnl_newfit_param[eta_bin_id][11]=0.0;
				fnl_newfit_paramER[eta_bin_id][11]=0.0;	  
			}

			delete stitchTest2;

		} else {

			TF1 *param_fitFun2=new TF1("fit_2",fitFunc3,2.0,500.0,2);
			param_fitFun2->SetParameter(0,1.0); // A0 param: term=A0/Pt
			param_fitFun2->SetParameter(1,0.0025); // param A1: term=A1
			h2_fnl->Fit(param_fitFun2,"EB","",fit_cut1,fit_cut2);
			gMinuit->mnmatu(1);
			corr_param[eta_bin_id][3]=gMinuit->fMATUvline[0];
			if(eta_bin_id<14)
			{
				fnl_fit_param[eta_bin_id][3]=param_fitFun2->GetParameter(0);
				fnl_fit_paramER[eta_bin_id][3]=param_fitFun2->GetParError(0);
				fnl_fit_param[eta_bin_id][4]=param_fitFun2->GetParameter(1);
				fnl_fit_paramER[eta_bin_id][4]=param_fitFun2->GetParError(1);
			} else {
				fnl_fit_param[eta_bin_id][3]=0.0;
				fnl_fit_paramER[eta_bin_id][3]=0.0;
				fnl_fit_param[eta_bin_id][4]=0.0;
				fnl_fit_paramER[eta_bin_id][4]=0.0;
			}


			for (int bin=1;bin <=60; ++bin)
			{
				binvals[eta_bin_id][1][bin-1] = h2_fnl->GetBinContent(bin);
			}

			
			delete param_fitFun2;
		}

		//________ parameterization for Landau mpv __________________
	
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

			gMinuit->mnmatu(1);
			fnl_newfit_param[eta_bin_id][12]=JayFitFunc3->GetParameter(0);
			fnl_newfit_paramER[eta_bin_id][12]=JayFitFunc3->GetParError(0);
			fnl_newfit_param[eta_bin_id][13]=JayFitFunc3->GetParameter(1);
			fnl_newfit_paramER[eta_bin_id][13]=JayFitFunc3->GetParError(1);
			fnl_newfit_param[eta_bin_id][14]=JayFitFunc3->GetParameter(2);
			fnl_newfit_paramER[eta_bin_id][14]=JayFitFunc3->GetParError(2);
			fnl_newfit_param[eta_bin_id][15]=JayFitFunc3->GetParameter(3);
			fnl_newfit_paramER[eta_bin_id][15]=JayFitFunc3->GetParError(3);
			fnl_newfit_param[eta_bin_id][16]=JayFitFunc3->GetParameter(4);
			fnl_newfit_paramER[eta_bin_id][16]=JayFitFunc3->GetParError(4);
			fnl_newfit_param[eta_bin_id][17]=JayFitFunc3->GetParameter(5);
			fnl_newfit_paramER[eta_bin_id][17]=JayFitFunc3->GetParError(5);
			fnl_newfit_param[eta_bin_id][18]=JayFitFunc3->GetParameter(6);
			fnl_newfit_paramER[eta_bin_id][18]=JayFitFunc3->GetParError(6);

			delete JayFitFunc3;


		} else {
			TF1 *param_fitFun3=new TF1("fit_3",fitFunc0,2.0,500.0,3);
			param_fitFun3->SetParameter(0,-0.01); // A0 param: term=A0
			param_fitFun3->SetParameter(1,0.0002); // param A1: term=A1*Pt
			param_fitFun3->SetParameter(2,-1.0); // param A2: term=A2/sqrt(Pt)
			h3_fnl->Fit(param_fitFun3,"EB","",fit_cut1,fit_cut2);

			gMinuit->mnmatu(1);
			if(TMath::Sqrt(TMath::Abs((gMinuit->fVhmat[2])*((gMinuit->fVhmat[0]))))>0.0) 
				corr_param[eta_bin_id][4]=(gMinuit->fVhmat[1])/TMath::Sqrt(TMath::Abs((gMinuit->fVhmat[2])*((gMinuit->fVhmat[0]))));
			corr_param[eta_bin_id][5]=gMinuit->fMATUvline[0];
			corr_param[eta_bin_id][6]=gMinuit->fMATUvline[1];
			fnl_fit_param[eta_bin_id][5]=param_fitFun3->GetParameter(0);
			fnl_fit_paramER[eta_bin_id][5]=param_fitFun3->GetParError(0);
			fnl_fit_param[eta_bin_id][6]=param_fitFun3->GetParameter(1);
			fnl_fit_paramER[eta_bin_id][6]=param_fitFun3->GetParError(1);
			fnl_fit_param[eta_bin_id][7]=param_fitFun3->GetParameter(2);
			fnl_fit_paramER[eta_bin_id][7]=param_fitFun3->GetParError(2);

			delete param_fitFun3;

			for (int bin=1;bin <=60; ++bin)
			{
				binvals[eta_bin_id][2][bin-1] = h3_fnl->GetBinContent(bin);
			}


		}

		
		//________ parameterization for Landau sigma __________________      

		if (bNewfits)
		{
			TF1 *stitchTest4 = new TF1("stitchTest4",fitFunc3_new,2,500,5);
			stitchTest4->SetParameter(0,1.6);
			stitchTest4->SetParameter(1,0.00005);
			stitchTest4->SetParameter(2,2); 
			stitchTest4->SetParameter(3,0.004);
			stitchTest4->SetParameter(4,30.0);
			stitchTest4->SetParLimits(4,20,40.0);
			h4_fnl->Fit(stitchTest4,"EBQ","",fit_cut1,fit_cut2);

			gMinuit->mnmatu(1);
			fnl_newfit_param[eta_bin_id][19]=stitchTest4->GetParameter(0);
			fnl_newfit_paramER[eta_bin_id][19]=stitchTest4->GetParError(0);
			fnl_newfit_param[eta_bin_id][20]=stitchTest4->GetParameter(1);
			fnl_newfit_paramER[eta_bin_id][20]=stitchTest4->GetParError(1);
			fnl_newfit_param[eta_bin_id][21]=stitchTest4->GetParameter(2);
			fnl_newfit_paramER[eta_bin_id][21]=stitchTest4->GetParError(2);
			fnl_newfit_param[eta_bin_id][22]=stitchTest4->GetParameter(3);
			fnl_newfit_paramER[eta_bin_id][22]=stitchTest4->GetParError(3);
			fnl_newfit_param[eta_bin_id][23]=stitchTest4->GetParameter(4);
			fnl_newfit_paramER[eta_bin_id][23]=stitchTest4->GetParError(4);

			delete stitchTest4;

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
			
			h4_fnl->Fit(param_fitFun4,"EB","",fit_cut1,fit_cut2);
			gMinuit->mnmatu(1);
			corr_param[eta_bin_id][7]=gMinuit->fMATUvline[0];
			fnl_fit_param[eta_bin_id][8]=param_fitFun4->GetParameter(0);
			fnl_fit_paramER[eta_bin_id][8]=param_fitFun4->GetParError(0);
			fnl_fit_param[eta_bin_id][9]=param_fitFun4->GetParameter(1);
			fnl_fit_paramER[eta_bin_id][9]=param_fitFun4->GetParError(1);

			delete param_fitFun4;
			
			for (int bin=1;bin <=60; ++bin)
			{
				binvals[eta_bin_id][3][bin-1] = h4_fnl->GetBinContent(bin);
			}

		}

		//________ parameterization for Gauss norm ____________________      
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

			gMinuit->mnmatu(1);
			if(eta_bin_id<14)
			{
				fnl_newfit_param[eta_bin_id][24]=JayFitFunc5->GetParameter(0);
				fnl_newfit_paramER[eta_bin_id][24]=JayFitFunc5->GetParError(0);
				fnl_newfit_param[eta_bin_id][25]=JayFitFunc5->GetParameter(1);
				fnl_newfit_paramER[eta_bin_id][25]=JayFitFunc5->GetParError(1);
				fnl_newfit_param[eta_bin_id][26]=JayFitFunc5->GetParameter(2);
				fnl_newfit_paramER[eta_bin_id][26]=JayFitFunc5->GetParError(2);
				fnl_newfit_param[eta_bin_id][27]=JayFitFunc5->GetParameter(3);
				fnl_newfit_paramER[eta_bin_id][27]=JayFitFunc5->GetParError(3);
				fnl_newfit_param[eta_bin_id][28]=JayFitFunc5->GetParameter(4);
				fnl_newfit_paramER[eta_bin_id][28]=JayFitFunc5->GetParError(4);
				fnl_newfit_param[eta_bin_id][29]=JayFitFunc5->GetParameter(5);
				fnl_newfit_paramER[eta_bin_id][29]=JayFitFunc5->GetParError(5);
				fnl_newfit_param[eta_bin_id][30]=JayFitFunc5->GetParameter(6);
				fnl_newfit_paramER[eta_bin_id][30]=JayFitFunc5->GetParError(6);
			} else {
				fnl_newfit_param[eta_bin_id][24]=0.0;
				fnl_newfit_paramER[eta_bin_id][24]=0.0;
				fnl_newfit_param[eta_bin_id][25]=0.0;
				fnl_newfit_paramER[eta_bin_id][25]=0.0;
				fnl_newfit_param[eta_bin_id][26]=0.0;
				fnl_newfit_paramER[eta_bin_id][26]=0.0;
				fnl_newfit_param[eta_bin_id][27]=0.0;
				fnl_newfit_paramER[eta_bin_id][27]=0.0;
				fnl_newfit_param[eta_bin_id][28]=0.0;
				fnl_newfit_paramER[eta_bin_id][28]=0.0;
				fnl_newfit_param[eta_bin_id][29]=0.0;
				fnl_newfit_paramER[eta_bin_id][29]=0.0;
				fnl_newfit_param[eta_bin_id][30]=0.0;
				fnl_newfit_paramER[eta_bin_id][30]=0.0;
			}

			delete JayFitFunc5;

		} else {

			TF1 *param_fitFun5=new TF1("fit_5",fitFunc5,2.0,500.0,3);
			param_fitFun5->SetParameter(0,0.22); // A0 param: term=A0
			param_fitFun5->SetParameter(1,-0.00043); // param A1: term=A1*Pt
			param_fitFun5->SetParameter(2,0.22); // param A2: term=A2/sqrt(Pt)
			h5_fnl->Fit(param_fitFun5,"EB","",fit_cut1,fit_cut2);

			gMinuit->mnmatu(1);
			if(TMath::Sqrt(TMath::Abs((gMinuit->fVhmat[2])*((gMinuit->fVhmat[0]))))>0.0) 
				corr_param[eta_bin_id][8]=(gMinuit->fVhmat[1])/TMath::Sqrt(TMath::Abs((gMinuit->fVhmat[2])*((gMinuit->fVhmat[0]))));
			corr_param[eta_bin_id][9]=gMinuit->fMATUvline[0];
			corr_param[eta_bin_id][10]=gMinuit->fMATUvline[1];
			if(eta_bin_id<14)
			{
				fnl_fit_param[eta_bin_id][10]=param_fitFun5->GetParameter(0);
				fnl_fit_paramER[eta_bin_id][10]=param_fitFun5->GetParError(0);
				fnl_fit_param[eta_bin_id][11]=param_fitFun5->GetParameter(1);
				fnl_fit_paramER[eta_bin_id][11]=param_fitFun5->GetParError(1);
				fnl_fit_param[eta_bin_id][12]=param_fitFun5->GetParameter(2);
				fnl_fit_paramER[eta_bin_id][12]=param_fitFun5->GetParError(2);
			} else {
				fnl_fit_param[eta_bin_id][10]=0.0;
				fnl_fit_paramER[eta_bin_id][10]=0.0;
				fnl_fit_param[eta_bin_id][11]=0.0;
				fnl_fit_paramER[eta_bin_id][11]=0.0;
				fnl_fit_param[eta_bin_id][12]=0.0;
				fnl_fit_paramER[eta_bin_id][12]=0.0;
			}

			delete param_fitFun5;

			for (int bin=1;bin <=60; ++bin)
			{
				binvals[eta_bin_id][4][bin-1] = h5_fnl->GetBinContent(bin);
			}

		}

		delete h1_fnl;
		delete h2_fnl;
		delete h3_fnl;
		delete h4_fnl;
		delete h5_fnl;

	}
	std::cout<<"Closing Files"<<std::endl;


	gROOT->cd();
	gDirectory->pwd();
	myFile1.cd();
	gDirectory->pwd();
	myFile1.Close();

	if (bNewfits)
	{
			std::cout << " ========================== NEW FIT FUNCTION PARAMETERS =========" << std::endl;
			for(int i=0; i<15; i++)  //loop over eta bin
			{
				for(int j=0; j<31; j++)  //loop over fit parameters
				{
					std::cout << "  jer_param[" << i << "][" << j << "]=" 
						<< fnl_newfit_param[i][j] << ";  jer_param_er[" << i << "]["
						<< j << "]=" << fnl_newfit_paramER[i][j] << ";" << std::endl;
				}
			}

	} else {

		std::cout << " ========================== OLD FIT FUNCTION PARAMETERS =========" << std::endl;
		for(int i=0; i<15; i++)  //loop over eta bin
		{
			for(int j=0; j<13; j++)  //loop over fit parameters
			{
				std::cout<<"  jer_param["<<i<<"]["<<j<<"]="<<fnl_fit_param[i][j]<<";  jer_param_er["<<i<<"]["<<j<<"]="<<fnl_fit_paramER[i][j]<<";"<<std::endl;
			}
		}

		for(int i=0; i<15; i++)
		{
			for(int j=0; j<11; j++)
			{
				if(j==0) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // Gauss mean: c01"<<std::endl;
				if(j==1) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // Gauss mean: c02"<<std::endl;
				if(j==2) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // Gauss mean: c12"<<std::endl;
				if(j==3) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // Gauss sigma: c01"<<std::endl;
				if(j==4) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // Landau mpv: c01"<<std::endl;
				if(j==5) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // Landau mpv: c02"<<std::endl;
				if(j==6) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // Landau mpv: c12"<<std::endl;
				if(j==7) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // Landau sigma: c01"<<std::endl;
				if(j==8) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // norm: c01"<<std::endl;
				if(j==9) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // norm: c02"<<std::endl;
				if(j==10) std::cout<<"  corr_param["<<i<<"]["<<j<<"]="<<corr_param[i][j]<<"; // norm: c12"<<std::endl;
			}
		}
	}

	// exact bin values for Ebin<12
	// i'll use exact values for E<6GeV
	
	if (! bNewfits )
	{
		std::cout << "\n" << " ====== Exact bin values for E<60GeV" << std::endl;
		std::cout << " ======= par0-4 == Gaus Mean/Sigma, Landau MPV/Sigma and Relative fraction " << std::endl;

		for (int etabin = 0; etabin<15; ++etabin)
		{
			for (int par=0; par <5; ++par)
			{
				for (int bin=0; bin<60;++bin)
				{
					if (binvals[etabin][par][bin] ==0) // find the closest upper/lower non-zero bin and use its value
					{
						//std::cout << "bin " << bin << " is zero!" << std::endl;
						int upjumps=0, downjumps=0;
						double upval = 0, downval = 0;
						bool foundupval = false, founddownval = false;
						
						for (int upbin=bin+1; upbin<60; ++upbin)
						{
							if (binvals[etabin][par][upbin] ==0) { upjumps++; continue; }
							else 
							{ 
								foundupval = true; 
								upval=  binvals[etabin][par][upbin] ; 
								break; 
							}
						}
						for (int downbin=bin-1; downbin>=0; --downbin)
						{
							if (binvals[etabin][par][downbin] ==0) { downjumps++; continue; }
							else 
							{ 
								founddownval = true;
								downval=  binvals[etabin][par][downbin]; 
								break; 
							}
						}

						//std::cout << "   upfound | val | jumps = " << foundupval 
						//			<< " | " << upval << " | " << upjumps << std::endl;
						//std::cout << "down found | val | jumps = " << founddownval 
						//			<< " | " << downval << " | " << downjumps << std::endl;
						
							
						if (foundupval && founddownval)
						{
							//std::cout << "case 1" << std::endl;
							if (upjumps<downjumps) binvals[etabin][par][bin] = upval;
							else binvals[etabin][par][bin] = downval;
						} else if (foundupval && !founddownval)
						{
							//std::cout << "case 2" << std::endl;
							binvals[etabin][par][bin] = upval;
						} else if (!foundupval && founddownval)
						{
							//std::cout << "case 3" << std::endl;
							binvals[etabin][par][bin] = downval;
						}

						//std::cout << "assigned value for bin [" << bin << "] =" << binvals[etabin][par][bin] << std::endl;

					}
					std::cout << "vBinVals["<< etabin << "][" << par << "][" << bin << "] = " << binvals[etabin][par][bin] << ";" << std::endl;
				}
			}

		}	
	}
	return;
}
