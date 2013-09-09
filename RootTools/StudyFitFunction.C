#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TPad.h>
#include <sstream>
#include <string>

// created this to study the fitting a custom function
// to a hist. This requires one of the calibration hist
// generated from MEtModel (JetFilterV2.cc/hh)


double FitFunc1(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0)
	{
		fitval=par[0]*exp(-1.0*par[1]*x[0]);
		//fitval=par[0]*exp(-1.0*par[1]*x[0]+par[2]*x[0]*x[0]) + par[3]*exp(-1.0*par[4]*x[0]);
		//fitval=par[0]*exp(-1.0*par[1]*x[0]+par[2]*x[0]*x[0]) + par[3]*exp(-1.0*par[4]*x[0]) + 1/(par[5]*x[0]*x[0]);
		//fitval=par[0]*exp(-1.0*par[1]*x[0]+par[2]*x[0]*x[0]);
	}
	else fitval=-1.0E6;
	
	return fitval;
}


void StudyFitFunction(const int njet=1, const int rebin=20, const float lox=0, const float upx=14)
{

	assert( (njet>=1 && njet<=4) && "invalid njets!");
	assert( (rebin > 0) && "Rebin must be > 0!");
	assert( (lox < upx) && "xlo must be smaller than xup limit!");

	//TFile *f = new TFile("MetTest_PhotonMC.root_1");
	TFile *f = new TFile("Merged.root");
	assert ( f != NULL && "File did not open !");
	
	std::stringstream name;
	name << "fMetSigCalib_estimate_njet" << njet;
	TH1* hist = (TH1*) f->FindObjectAny(name.str().c_str());
	hist->SetDirectory(0);
	assert (hist != NULL && "hist not found!");


	gStyle->SetOptStat("nemourm");
	gStyle->SetOptFit(11111);

	new TCanvas("MyCanvas");
	//gPad->SetLogy();
	hist->Rebin(rebin);
	hist->SetMarkerSize(1);
	hist->SetMarkerStyle(20);
	hist->SetMarkerColor(4);
	//hist->Fit("gaus");
	//TF1 *f1 = new TF1("F1","[0]+[1]*x+[2]x*x+[3]x*x*x+[4]x*x*x*x4+[5]*exp([6]*x)+",0,16);
	//TF1 *f1 = new TF1("F1",FitFunc1,0,20.0,5);
	TF1 *f1 = new TF1("F1",FitFunc1,0,20.0,2);
	f1->Draw();
		//f1->SetParameters(0.05,1.5,3,-1,3);
	//	f1->SetParameter(0,0.03 * rebin);
	//	f1->SetParameter(1,0.1);
	//	f1->SetParameter(2,-0.2);
	//	f1->SetParameter(3,0.01);
	//	f1->SetParameter(4,0.01);
		//f1->SetParameter(5,0.003);
	
	//hist->Fit("F1");
	//hist->Fit(f1,"+E","",lox,upx);
/*
		f1->SetParameter(0,0.03);
		f1->SetParameter(1,0.1);
		f1->SetParameter(2,-0.2);
		f1->SetParameter(3,0.01);
		f1->SetParameter(4,0.01);
*/	
	//hist->Fit("F1");
	//hist->Fit(f1,"E","",lox,upx);
	//hist->Draw("ep");
	delete f;
}

