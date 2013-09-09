#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TPad.h"

/*
 * Run this on caf to get limits. This generates a text file which
 * is needed by getlimits.C to derive the observed and expected limits.
 * 04-30-2010 - Sam 
 */


using namespace std;

//const double pi = 3.14159265;
const double pi = TMath::Pi();

double fun(double *x, double *par){
	//x[0]: background
	//par[0]: mean background
	//par[1]: sys error
	//par[2]: num of events
	//double f = pow(x[0],par[2])*exp(-x[0])/TMath::Factorial(int(par[2]+0.001))/sqrt(2*pi)/par[1]*exp(-pow(x[0]-par[0],2)/2.0/par[1]/par[1]);

	double xx = x[0];
	double gmean = par[0];
	double gsigma = par[1];
	double pmean = par[2];
	
	
	// Double_t Gaus(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
   double pois = TMath::Poisson(pmean, xx);
	double gaus = TMath::Gaus(xx, gmean, gsigma, 1);
	double ff = pois * gaus;
	
	//return f;
	return ff;
}

void scan_caf()
{

	/*	TFile f1("hist.root");
		TFile f2("SysTot.root");

		TH1D *datahist = (TH1D*)f1.Get("datahist");
		TH1D *fithist = (TH1D*)f1.Get("fithist");

		TH1D *totlow = (TH1D*)f2.Get("totlow");
		TH1D *tothigh = (TH1D*)f2.Get("tothigh");
		*/

	TFile f("Output.root");
	TH1F *datahist = dynamic_cast<TH1F*> (f.Get("InvMass"));
	TH1F *fithist = dynamic_cast<TH1F*> (f.Get("hist_err"));
	TH1F *tothigh = dynamic_cast<TH1F*> (f.Get("SystPlus"));
	TH1F *totlow = dynamic_cast<TH1F*> (f.Get("SystMinus"));

	TF1 *fun1 = new TF1("fun",fun,0,1000,3);

	TH1D *p = new TH1D("p","Probability",100,150,650);

	/*
		for (int i = 1; i<=fithist->GetNbinsX(); i++){
		double mass = fithist->GetBinCenter(i);
		if (mass>150&&mass<650){
		double sigma = sqrt(2.0)*sqrt(pow(0.135*sqrt(mass/2.0),2)+
		pow(0.02*(mass/2.0),2));
	//cout<<mass<<" "<<sigma<<endl;
	int bin1 = fithist->FindBin(mass-sigma/2);
	int bin2 = fithist->FindBin(mass+sigma/2);
	//cout<<mass<<" "<<bin1<<" "<<bin2<<endl;
	double data = 0;
	double bg = 0;
	double err = 0;
	for (int j = bin1; j<=bin2; j++){
	data+=datahist->GetBinContent(j);
	bg+=fithist->GetBinContent(j);
	double err1 = -totlow->GetBinContent(j);
	double err2 = tothigh->GetBinContent(j);
	err+=TMath::Max(err1,err2)*bg;
	}
	//cout<<mass<<" "<<data<<" "<<bg<<" "<<err<<endl;
	double prob = 0;
	fun1->SetParameter(0,bg);
	fun1->SetParameter(1,err);
	for (int j = int(data+0.001); j<100; j++){
	fun1->SetParameter(2,j);
	prob += fun1->Integral(TMath::Max(0.0,bg-10*err),bg+10*err);
	}
	//cout<<mass<<" "<<prob<<endl;
	p->SetBinContent(p->FindBin(mass),prob);
	}
	}
	*/

	delete gRandom;
	gRandom = (TRandom*) new TRandom3;
	cout<<"test"<<endl;
	int ind=0;
	if (gSystem->Getenv("CAF_SECTION")) sscanf(gSystem->Getenv("CAF_SECTION"),"%d",&ind);
	cout<< "caf_section = " << ind <<endl;

	gRandom->SetSeed(ind);

	TH1D *minp = new TH1D("minp","Minimum Probability of Each PseudoExpt",100,0,0.2);

	//int nexp = 50000;
	int nexp = 100;

	TH1D *htemp = (TH1D*)datahist->Clone("htemp");

	ofstream out;
	out.open(Form("minprob_%d.txt",ind));

	for (int iexp = 0; iexp<nexp; iexp++){
		//if (iexp%10==0) cout<<iexp<<endl;
		//generate pseudo-experiments
		htemp->Reset();
		for (int i = 1; i<=htemp->GetNbinsX(); i++){
			double mass = htemp->GetBinCenter(i);
			if (mass>150&&mass<650){
				double bg = fithist->GetBinContent(i);
				double err1 = -totlow->GetBinContent(i);
				double err2 = tothigh->GetBinContent(i);
				double err = TMath::Max(err1,err2)*bg;
				double mean = gRandom->Gaus(bg,err);
				if (mean<0) mean = 0;
				htemp->SetBinContent(i,gRandom->Poisson(mean));
			}
		}
		double minprob = 2.;
		for (int i = 1; i<=fithist->GetNbinsX(); i++){
			double mass = fithist->GetBinCenter(i);
			if (mass>150&&mass<650){
				double sigma = sqrt(2.0)*sqrt(pow(0.135*sqrt(mass/2.0),2)+
						pow(0.02*(mass/2.0),2));
				//cout<<mass<<" "<<sigma<<endl;
				int bin1 = fithist->FindBin(mass-sigma/2);
				int bin2 = fithist->FindBin(mass+sigma/2);
				//cout<<mass<<" "<<bin1<<" "<<bin2<<endl;
				double data = 0;
				double bg = 0;
				double err = 0;
				for (int j = bin1; j<=bin2; j++){
					data+=htemp->GetBinContent(j);
					bg+=fithist->GetBinContent(j);
					double err1 = -totlow->GetBinContent(j);
					double err2 = tothigh->GetBinContent(j);
					//err+=TMath::Max(err1,err2)*bg;
					err+=TMath::Max(err1,err2)*fithist->GetBinContent(j); //above line has a bug, bg should be this bin
				}
				//cout<<mass<<" "<<data<<" "<<bg<<" "<<err<<endl;
				double prob = 0;
				fun1->SetParameter(0,bg);
				fun1->SetParameter(1,err);
				for (int j = int(data+0.001); j<100; j++){
					fun1->SetParameter(2,j);
					prob += fun1->Integral(TMath::Max(0.0,bg-10*err),bg+10*err);
				}
				if (prob<minprob) minprob=prob;
			}
		}
		minp->Fill(minprob);
		out<<minprob<<endl;
	}

	out.close();
	//  TCanvas *c1 = new TCanvas("c1","c1");
	//  TH2D *fr = new TH2D("fr","",100,150,650,100,1e-6,2);
	//  fr->DrawCopy();
	//  p->DrawCopy("same");
	//  gPad->SetLogy();
	//
	TCanvas *c2 = new TCanvas("c2","c2");
	minp->DrawCopy();
	cout<<minp->GetMean()<<endl;
}
