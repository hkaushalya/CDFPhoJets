#include <iomanip>
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TFolder.h"
#include "TMinuit.h"
#include "TStyle.h"
#include <sstream>
#include "TVirtualPad.h"
#include <TPaveText.h>



// Jay's two 1/2 gausian fit function
/*
 * 
 * 
 */
/*
	double JayFuntion()
	{
//slp1 is positive, slp2 is negative
// Evaluate first smeared exponential (upwards):
double smrpb1 = smrpb(dPt,sigma,slp1,slp2);

// Evaluate second smeared exponential (downward):
double smrpb2=smrpb(dPt,sigma,slp2,slp1);

return (smrpb1+smrpb2)/2.0
}
*/


//    function smrpb(dPt,sigma,a,b)
//    how do I derive
/*double smrpb(const double dPt, const double sigma,
  const double a, const double b)
  {

  double precision  smrpb,dPt,sigma,a,b,sigmaTmp
  double precision  dfreq,arg,x_log,p_log,errfunc
  double precision  ROOT2PI
  double ROOT2PI= 2*TMath::Pi();

  if(a == 0)
  {
  arg   =   dPt/sigma;
  arg   = - (arg**2.)/2;
//smrpb = 1.0d+00/(sigma*ROOT2PI) * dexp(arg)
smrpd = Exp(arg)/(sigma*ROOT2PI);
} else
{
arg      = (dPt - (a+b)/2)/sigma + (sigma/a);
//errfunc  = dfreq(arg)
errfunc  = Erf(arg);
if (a>0.0) errfunc = 1-errfunc

if(errfunc>1E-37)
{
//x_log =-0.5+(0.5*(sigma/a)**2.)+(dPt-b/2)/a
x_log = -0.5 + (0.5 * pow(sigma/a,2)) + (dPt-b/2)/a;
//p_log = x_log + dlog(errfunc)  ??is this log=log10
p_log = x_log + Log(errfunc);
//smrpb = (1./abs(a)) * dexp(p_log)
smrpb = Exp(p_log)/fabs(a); 
} else
{
sigmaTmp = sqrt( pow(sigma,2) + pow(a,2));
arg      = dPt/sigmaTmp;
arg      = -1*(pow(arg,2))/2.;
smrpb    = Exp(arg)/(sigma*ROOT2PI);
}
}
return smrpb;
}
*/

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

//------ Gauss mean: A0+A1*Pt+A3/sqrt(Pt) ---- default
double fitFunc0(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0) fitval=par[0]+par[1]*x[0]+par[2]/x[0];
	else fitval=-1.0E6;
	return fitval;
}
//------ Gauss mean: A0+A1*Pt+A3/sqrt(Pt) ---- default
double fitFunc0_syst(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0) fitval=par[0]+par[1]*x[0]+par[2]/sqrt(x[0]);
	else fitval=-1.0E6;
	return fitval;
}
//------ Gauss norm: A0+A1*Pt+A3/ln(Pt) ---- default
double fitFunc5(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0) fitval=(par[0]+par[1]*sqrt(x[0]))/x[0]+par[2];
	else fitval=-1.0E6;
	return fitval;
}

//----- Gauss Sigma fit function: sqrt(A0*A0/Pt+A1*A1/(Pt*Pt)+A2*A2)
double fitFunc1(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0)
	{
		double arg=0.0;
		arg=(par[0]/x[0]+par[1]/(x[0]*x[0])+par[2]);
		fitval=sqrt(arg);
	}
	else fitval=-1.0E6;
	return fitval;
}
// //----- Landau Sigma fit function: sqrt(A0/Pt+A2) ---- try as default
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
//----- Landau Sigma fit function: A0/sqrt(Pt)+A2 ---- try as systematics 
double fitFunc3_syst(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0)
	{
		double arg=0.0;
		arg=par[0]/sqrt(x[0])+par[1];
		fitval=arg;
	}
	else fitval=-1.0E6;
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
void JetResFitPerEPerEtaBin(int e_bin_id=0, int eta_bin_id=0,
		int folder_id=0, TVirtualPad* pad =0)
{
	if (e_bin_id <0 || e_bin_id>99)
	{
		std::cout << "Invalid E bin. must be between 0 and 99. Returning"<< std::endl;
		return;
	}
	if (eta_bin_id <0 || eta_bin_id>14)
	{
		std::cout << "Invalid Eta bin. must be between 0 and 14. Returning"<< std::endl;
		return;
	}

	if (! (folder_id == 0 || folder_id == 3 || folder_id == 4))
	{
		std::cout << "Invalid folder selection. must be 1,3 or 4. Returning" << std::endl;
		return;
	}


	gROOT->Reset();

	char foldername[200];
	char fObject[200];

	// folders JetRes_syst1, JetRes_syst2 does not exits in the
	// root file. 
	// JetRes_def is same as JetRes_def_nvx1. So these were 
	// derived from events with a good vertex.
	// These JERs were derived from diJet MC.
	// 
	if(folder_id==0) sprintf(foldername,"JetRes_def");
	//if(folder_id==1) sprintf(foldername,"JetRes_syst1");
	//if(folder_id==2) sprintf(foldername,"JetRes_syst2");
	if(folder_id==3) sprintf(foldername,"JetRes_def_nvx1");
	if(folder_id==4) sprintf(foldername,"JetRes_def_nvx2");

	double gauss_mean;
	double gauss_mean_er;
	double gauss_sigma;
	double gauss_sigma_er;
	double landau_mean;
	double landau_mean_er;
	double landau_sigma;
	double landau_sigma_er;
	double landau_norm;
	double landau_norm_er;
	double fit_chi2;
	double fit_chi2_er;
	double avePt;
	double avePt_er;
	double histo_Nentry;
	int fit_status;

	double maxFit_X; // MPV for fit function
	double histoRMS; // RMS of histogram 

	//int h_bin_first=0;
	int h_bin_last=100;

	//________________________________________________________________________________
	//------------------- stopped here -----------------------------------------------
	//________________________________________________________________________________

	TH1F *h1=0;

	gDirectory->pwd();

	//-------------- matching in dR<0.3
	//   TFile myFile1("results/JER_hadE_sample1_071307.root","r");
	//   TFile myFile1("results/JER_hadE_sample2_071307.root","r");
	//   TFile myFile1("results/JER_hadE_sample3_071307.root","r");
	//   TFile myFile1("results/JER_hadE_sample4_071307.root","r");

	//-------------- matching in dR<0.1
	TFile myFile1("JER_hadE_sample1_dR01_072507.root","r");
	//   TFile myFile1("results/JER_hadE_sample2_dR01_072507.root","r");
	//   TFile myFile1("results/JER_hadE_sample3_dR01_072507.root","r");
	//   TFile myFile1("results/JER_hadE_sample4_dR01_072507.root","r");


	TFolder* myFolder1 = (TFolder* ) myFile1.FindObjectAny("Ana");

	sprintf(fObject,"MyJetFilter/Hist/%s/fJetRes_AveJetE[%i]",foldername,eta_bin_id); // get object name
	TProfile *_avePt_histo= dynamic_cast<TProfile*> (myFolder1->FindObjectAny(fObject));
	assert (_avePt_histo != NULL && "avgPt profile hist is null!");

	sprintf(fObject,"MyJetFilter/Hist/%s/fJetRes_TrueBal[%i][%i]",foldername,eta_bin_id,e_bin_id); // get object name
	h1= dynamic_cast<TH1F*> (myFolder1->FindObjectAny(fObject));
	assert (h1 != NULL && "h1 is null");

	avePt 		= _avePt_histo->GetBinContent(e_bin_id+1);
	avePt_er 	= _avePt_histo->GetBinError(e_bin_id+1);
	histo_Nentry = h1->GetEntries();

	if(histo_Nentry>99) h_bin_last = e_bin_id;

	gauss_mean     = 0.0;
	gauss_mean_er  = 0.0;
	gauss_sigma    = 0.0;
	gauss_sigma_er = 0.0;
	landau_mean    = 0.0;
	landau_mean_er = 0.0;
	landau_sigma   = 0.0;
	landau_sigma_er= 0.0;
	landau_norm 	= 0.0;
	landau_norm_er = 0.0;
	fit_chi2 	   = 0.0;
	fit_chi2_er    = 0.0; 
	fit_status     = 0;

	maxFit_X = 0.0;
	histoRMS = 0.01;

	gROOT->cd();
	gDirectory->pwd();

	std::cout<<"---------------------------------------------"<<std::endl;
	std::cout<<"     Fitting bin="<<e_bin_id<<std::endl;
	std::cout<<"_____________________________________________"<<std::endl;

	if (histo_Nentry >= 10000.0) h1->Rebin(2);
	if (histo_Nentry >= 1000.0 && histo_Nentry < 10000.0) h1->Rebin(3);
	if (histo_Nentry >= 200.0 && histo_Nentry < 1000.0) h1->Rebin(3);
	if (histo_Nentry >= 100.0 && histo_Nentry < 200.0) h1->Rebin(4);
	if (eta_bin_id < 14 && histo_Nentry < 100.0) // skipping bins with Nevnt<100
	{
		std::cout << "Returning because I failed the 'eta_bin_id<14 && histo_Nentry<100.0' cut" << std::endl; 
		std::cout << "Eta bin = " << eta_bin_id << "\thist Nentries = " << histo_Nentry << std::endl; 
		//return;
	}
	if (eta_bin_id == 14 && histo_Nentry < 80.0) // skipping bins with Nevnt<100
	{
		std::cout << "Returning because I failed the 'eta_bin_id<14 && histo_Nentry<80.0' cut" << std::endl; 
		std::cout << "Eta bin = " << eta_bin_id << "\thist Nentries = " << histo_Nentry << std::endl; 
		//return;
	}

	if (histo_Nentry < 100.0) h1->Rebin(4);

	int nbins   	= h1->GetNbinsX();
	float binwidth = h1->GetBinWidth(1);
	float lowedge  = (h1->GetBinCenter(1)) - binwidth/2.0;
	float highedge = lowedge + nbins * binwidth;

	char title[200];
	char name[200];

	sprintf(name,"bal");
	if (eta_bin_id < 14) 
	{
		sprintf(title,"E^{det}/E^{had}-1.0, %2.1f<|#eta_{det}^{jet}|<%2.1f and %i<E^{jet}<%i",0.2*eta_bin_id,0.2*(eta_bin_id+1),5*e_bin_id,5*(e_bin_id+1));
	} else
	{
		sprintf(title,"E^{det}/E^{had}-1.0, |#eta_{det}^{jet}|>%2.1f and %i<E^{jet}<%i",0.2*eta_bin_id,5*e_bin_id,5*(e_bin_id+1));
	}
	TH1F* _h= new TH1F(name,title,nbins,lowedge,highedge);
	_h->Sumw2();

	float cut1 = -1.0;
	float cut2 = 2.0;
	float val  = 0.0;
	float err  = 0.0;
	int bin_first = 0;
	int bin_last  = 0;

	for (int j=0; j<nbins; j++)
	{
		val = h1->GetBinContent(j+1);
		err = h1->GetBinError(j+1);
		_h->SetBinContent(j+1,val);
		_h->SetBinError(j+1,err);
		if (bin_first == 0 && val > 0.0) bin_first = j + 1;
		//       if(val>0.0) bin_last=i+1; 
		if (val > 0.0 && h1->GetBinCenter(j+1) < 3.5) bin_last = j + 1; 
	}

	cut1 = h1->GetBinCenter(bin_first);
	cut2 = h1->GetBinCenter(bin_last);

	int maxbin = h1->GetMaximumBin();
	double hight = h1->GetBinContent(maxbin);
	double mean_histo = h1->GetMean();
	double rms_histo = h1->GetRMS();



	TF1 *fitFun1 = new TF1("fit1", myFunc1, cut1 - binwidth, cut2 + binwidth, 6);
	fitFun1->SetLineWidth(2);
	if (eta_bin_id < 14)
	{
		fitFun1->SetParameter(0, 0.9 * hight); 	// norm 
		fitFun1->SetParameter(1, mean_histo); 		// mean of Gauss
		fitFun1->SetParameter(2, 0.7 * rms_histo);// Sigma of Gauss
		fitFun1->SetParameter(3, mean_histo); 		// MPV of Landau
		fitFun1->SetParameter(4, 0.5 * rms_histo); // Sigma of Landau
		fitFun1->SetParameter(5, 0.2); 				// Landau norm
	} else
	{
		fitFun1->SetParameter(0, 0.9 * hight); 	// norm 
		fitFun1->FixParameter(1, 0.0); 				// mean of Gauss
		fitFun1->FixParameter(2, 0.0); 				// Sigma of Gauss
		fitFun1->SetParameter(3, mean_histo); 		// MPV of Landau
		fitFun1->SetParameter(4, 0.5 * rms_histo); // Sigma of Landau
		fitFun1->FixParameter(5, 0.0); 				// Landau norm
	}

	fitFun1->SetParName(0,"Gauss: const");
	fitFun1->SetParName(1,"Gauss: mean");
	fitFun1->SetParName(2,"Gauss: #sigma");
	fitFun1->SetParName(3,"Landau: MPV");
	fitFun1->SetParName(4,"Landau: #sigma");
	fitFun1->SetParName(5,"Gaus: const");


	// if no canvas is given, use these settings.
	// the Fit method will autmatically
	// create a hist pad with these settings
	if (pad == 0)
	{
		//default canvas settings
		gStyle->SetCanvasDefH(600);
		gStyle->SetCanvasDefW(900);
		gStyle->SetAxisColor(kBlue);
		gStyle->SetPadTickY(1);
		gStyle->SetPadTickX(1);
		gStyle->SetOptLogy(1);
		gStyle->SetOptFit(1);
		gStyle->SetPadGridX(1);
		gStyle->SetPadGridY(1);
		gStyle->SetTitleFontSize(0.035);
		gStyle->SetStatFontSize(0.03);
		gStyle->SetTitleFont(2);
	}
	_h->Fit(fitFun1,"B","",cut1,cut2);

	//---------- re-fitting second time --------------------------
	if (eta_bin_id < 14)
	{
		fitFun1->SetParameter(0,fitFun1->GetParameter(0)); // norm 
		fitFun1->SetParameter(1,fitFun1->GetParameter(1)); // mean of Gauss
		fitFun1->SetParameter(2,fitFun1->GetParameter(2)); // Sigma of Gauss
		fitFun1->SetParameter(3,fitFun1->GetParameter(3)); // MPV of Landau
		fitFun1->SetParameter(4,fitFun1->GetParameter(4)); // Sigma of Landau
		fitFun1->SetParameter(5,fitFun1->GetParameter(5)); // Landau norm
	} else
	{
		fitFun1->SetParameter(0, fitFun1->GetParameter(0)); // norm 
		fitFun1->FixParameter(1, 0.0); // mean of Gauss
		fitFun1->FixParameter(2, 0.0); // Sigma of Gauss
		fitFun1->SetParameter(3, fitFun1->GetParameter(3)); // MPV of Landau
		fitFun1->SetParameter(4, fitFun1->GetParameter(4)); // Sigma of Landau
		fitFun1->FixParameter(5, 0.0); // Landau norm	  
	}

	//	std::stringstream canvasname;
	//	canvasname << "HadJER_" << e_bin_id << "_" << eta_bin_id;
	//	TCanvas* c1 = new TCanvas(canvasname.str().c_str(),"Jet energy resolution study",200,10,1000,700);
	//	gPad->SetTicky();
	//	gPad->SetTickx();
	//	gPad->SetLogy();
	//	gStyle->SetOptFit(1);
	//	gPad->SetGridx();
	//	gPad->SetGridy();
	//	gStyle->SetTitleFontSize(0.035);
	//	gStyle->SetStatFontSize(0.03);
	//	gStyle->SetTitleFont(2);

	_h->GetXaxis()->SetAxisColor(kBlue);
	_h->GetYaxis()->SetAxisColor(kBlue);

	_h->Fit(fitFun1,"EB","",cut1,cut2);

	_h->SetLineColor(2);
	_h->SetLineWidth(1);
	_h->SetMarkerColor(2);
	_h->SetMarkerStyle(20);
	_h->SetMarkerSize(0.6);
	_h->SetMinimum(0.1);
	//   _h->SetMaximum(0.2);
	//_h->Draw("AP");

	//---------- end of re-fitting -------------------------------

	if((gMinuit->fCstatu).Contains("SUCCESSFUL")==true) fit_status=3;
	if((gMinuit->fCstatu).Contains("PROBLEMS")==true) fit_status=2;
	if((gMinuit->fCstatu).Contains("FAILURE")==true) fit_status=0;
	std::cout<<" Minuit Error Matrix Status ="<<gMinuit->fCstatu<<std::endl;
	std::cout<<" Minuit Error Matrix Status ="<<fit_status<<std::endl;


	TPaveText *tp = new TPaveText(0.6,0.4,0.9,0.5,"NDC");
	tp->SetLineColor(kRed);
	tp->AddText(gMinuit->fCstatu);
	tp->Draw();
	if (pad == 0)
	{
		std::stringstream filename;
		filename << "HadJER_" << e_bin_id << "_" << eta_bin_id << ".pdf";
		//gPad->Print(filename.str().c_str(),"pdf");
	}




	//       gMinuit->mnmatu(1); 
	maxFit_X=fitFun1->GetMaximumX(-1.0,1.0);
	histoRMS=rms_histo;
	gauss_mean=fitFun1->GetParameter(1);
	gauss_mean_er=fitFun1->GetParError(1);
	gauss_sigma=fitFun1->GetParameter(2);
	gauss_sigma_er=fitFun1->GetParError(2);
	landau_mean=fitFun1->GetParameter(3);
	landau_mean_er=fitFun1->GetParError(3);
	landau_sigma=fitFun1->GetParameter(4);
	landau_sigma_er=fitFun1->GetParError(4);
	landau_norm=fitFun1->GetParameter(5);
	landau_norm_er=fitFun1->GetParError(5);

	if(fitFun1->GetNDF()>0)
		fit_chi2=(1.0*fitFun1->GetChisquare())/(1.0*fitFun1->GetNDF());
	else fit_chi2=-1.0;

	//   delete fitFun1;
	//   delete _h;



	//   TCanvas* c1 = new TCanvas("c1","Jet energy resolution study",200,10,700,500);
	// //   c1->Divide(2,3);
	// //   c1->cd(1);
	//   gPad->SetTicky();
	//   gPad->SetTickx();
	//   h1_fnl->SetLineColor(2);
	//   h1_fnl->SetLineWidth(1);
	//   h1_fnl->SetMarkerColor(2);
	//   h1_fnl->SetMarkerStyle(20);
	//   h1_fnl->SetMarkerSize(0.6);
	//   h1_fnl->SetMinimum(-0.2);
	//   h1_fnl->SetMaximum(0.2);
	//   TF1 *param_fitFun1=new TF1("fit_1",fitFunc0,2.0,500.0,3);
	//   param_fitFun1->SetParameter(0,0.09); // A0 param: term=A0
	//   param_fitFun1->SetParameter(1,-0.0002); // param A1: term=A1*Pt
	//   param_fitFun1->SetParameter(2,-1.0); // param A2: term=A2/Pt
	//   h1_fnl->Fit(param_fitFun1,"EB","",fit_cut1,fit_cut2);
	//   h1_fnl->Draw("AP");
	//   h1_fnl->GetXaxis()->SetTitle("average E_{T}^{jet}");
	//   h1_fnl->GetYaxis()->SetTitle("Gauss mean");

	std::stringstream filenamegif;
	float minE = e_bin_id * 5;
	float maxE = minE + 5.;
	filenamegif << "HadJes_ExcludedEta" << eta_bin_id <<"_"<< minE << "E" << maxE << ".gif";
	if (pad != NULL)	pad->Print(filenamegif.str().c_str());


	
	std::cout<<"Closing Files"<<std::endl;
	gDirectory->pwd();

	myFile1.cd();
	gDirectory->pwd();
	myFile1.Close();
}



void AllEtaInPad(const int iEbin=0, const int iFolder=0)
{
	// This is a hack I did to get all the Eta bins for a
	// given E bin into one canvas. The canvas will have all
	// 15 Eta bins.
	//
	if (iEbin <0 || iEbin>99)
	{
		std::cout << "Invalid E bin. must be between 0 and 99. Returning"<< std::endl;
		return;
	}

	gStyle->SetCanvasDefH(1000);
	gStyle->SetCanvasDefW(1200);
	gStyle->SetAxisColor(kBlue);
	gStyle->SetPadTickY(1);
	gStyle->SetPadTickX(1);
	gStyle->SetOptLogy(1);
	gStyle->SetOptFit(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	//gStyle->SetTitleFontSize(0.035);
	//gStyle->SetStatFontSize(0.03);
	//gStyle->SetTitleFont(2);

	TCanvas *c1 = new TCanvas("C1","Jet Energy Resolution");
	c1->Divide(3,5);

	for (int eta=0; eta<=14; ++eta)
	{
		c1->cd(eta+1);
		TVirtualPad* pad = c1->GetSelectedPad();
		JetResFitPerEPerEtaBin(iEbin,eta,iFolder,pad);
	}
	c1->cd();
	std::stringstream filenamepdf, filenamegif;
	float minE = iEbin*5;
	float maxE = minE + 5.;
	filenamepdf << "HadJes_AllEtaBinsForEbin_" << minE << "_" << maxE << ".pdf";
	filenamegif << "HadJes_AllEtaBinsForEbin_" << minE << "_" << maxE << ".gif";
	c1->Print(filenamepdf.str().c_str());
	c1->Print(filenamegif.str().c_str());
}

void AllEbinsByEta(const int e_from=0, const int e_to=99,const int iFolder=0)
{
	for (int ebin=e_from;ebin<=e_to; ++ebin)
	{
		AllEtaInPad(ebin,iFolder);
	}
}
