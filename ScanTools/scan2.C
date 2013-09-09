//CVS Log 
/*{{{*/
/************************************************************************
 * $Id: scan2.C,v 1.1 2013/02/28 03:43:46 samantha Exp $
 * $Log: scan2.C,v $
 * Revision 1.1  2013/02/28 03:43:46  samantha
 * Final commit. no checks.! these were never commited to cvs before!
 *
 ************************************************************************/
/*}}}*/

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLatex.h"
#include <iostream>
#include <math.h>

using namespace std;

const double pi = 3.14159265;
int c=0;

double fun(double *x, double *par){
  //x[0]: background
  //par[0]: mean background
  //par[1]: sys error
  //par[2]: num of events
  //double f = pow(x[0],par[2])*exp(-x[0])/TMath::Factorial(int(par[2]+0.001))/sqrt(2*pi)/par[1]*exp(-pow(x[0]-par[0],2)/2.0/par[1]/par[1]);
  double mean = par[0];
  double sigma = par[1];
  double N = par[2];
  const bool normalize = true;		//normalize the gausian
  
  double numer1  = pow(x[0],par[2]);
  double f = 0;
  if (isinf(numer1)) f=0;
  else f = pow(x[0],par[2])*exp(-x[0])/TMath::Factorial(int(par[2]+0.001))/sqrt(2*pi)/par[1]*exp(-pow(x[0]-par[0],2)/2.0/par[1]/par[1]);
 

	double xx = x[0];
	double pmean = par[2];
	double gmean = par[0];
	double gsigma = par[1];
	
	
	// Double_t Gaus(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
   double pois = TMath::Poisson(pmean, xx);
	double gaus = TMath::Gaus(xx, gmean, gsigma, 1);

	double ff = pois * gaus;
	cout << "f/ff= " << f << "/ " << ff  << endl;
	//if (ff/f==0) cout << "MATCH [old/new = " << f << "/ " << ff << "]" << endl;
	//else cout << " NO NO MATCH [old/new = " << f << "/ " << ff << "]" << endl;
	 
/*  if (isnan(f))
	{
		if (c<4)
		{

			
  double numer  = pow(x[0],par[2])*exp(-x[0]);
  double numer1  = pow(x[0],par[2]);
	double numer2	= exp(-x[0]);
  double denom1 = TMath::Factorial(int(par[2]+0.001));
  double denom2 = sqrt(2*pi);
  double denom3 = par[1]*exp(-pow(x[0]-par[0],2)/2.0/par[1]/par[1]);
	std::cout << "x, parp2 = " << x[0] << ", " <<  par[2] << std::endl;
  cout << "numer = " << numer << endl;
  cout << "numer1 = " << numer1 << endl;
  cout << "numer2 = " << numer2 << endl;
  cout << "denom1 = " << denom1 << endl;
  cout << "denom2 = " << denom2 << endl;
  cout << "denom3 = " << denom3 << endl;
  //double f = numer/denom1/denom2/denom3;
  cout << "f = " << f << endl;
  ++c;
		}
	}
  */
  return f;
}

void scan2(){

/*  TFile f1("hist.root");
  TFile f2("SysTot.root");

  TH1D *datahist = (TH1D*)f1.Get("datahist");
  TH1D *fithist = (TH1D*)f1.Get("fithist");
  
  TH1D *totlow = (TH1D*)f2.Get("totlow");
  TH1D *tothigh = (TH1D*)f2.Get("tothigh");
*/
  
	TF1 *fun1 = new TF1("fun",fun,0,1000,3);

	TFile f("Output.root");
	TH1F *datahist = dynamic_cast<TH1F*> (f.Get("InvMass"));
	TH1F *fithist = dynamic_cast<TH1F*> (f.Get("hist_err"));
	TH1F *tothigh = dynamic_cast<TH1F*> (f.Get("SystPlus"));
	TH1F *totlow = dynamic_cast<TH1F*> (f.Get("SystMinus"));
	new TCanvas();
	datahist->DrawCopy();
	fithist->DrawCopy("same");
	gPad->SetEditable(0);

/*	if (datahist->IsZombie()) { std::cout << "Err! datahist not found " << std::endl; return; }
	if (fithist->IsZombie()) { std::cout << "Err! err_hist not found " << std::endl; return; }
	if (tothigh->IsZombie()) { std::cout << "Err! syst_plus hist not found " << std::endl; return; }
	if (totlow->IsZombie()) { std::cout << "Err! syst_minus hist not found " << std::endl; return; }
*/

  TH1D *p = new TH1D("p","Probability",100,150,900);

  for (int i = 1; i<=fithist->GetNbinsX(); i++)
  {
    double mass = fithist->GetBinCenter(i);
  
	 if (mass>440 && mass<450)
	 {
      double sigma = sqrt(2.0)*sqrt(pow(0.135*sqrt(mass/2.0),2)+
				    pow(0.02*(mass/2.0),2));
      cout<< " ===================== mass +/- sigma = " << mass<<"+/-"<<sigma<<endl;
      int bin1 = fithist->FindBin(mass-sigma/2);
      int bin2 = fithist->FindBin(mass+sigma/2);
      cout<<mass<<" "<<bin1<<" "<<bin2<<endl;
      double data = 0;
      double bg = 0;
      double err = 0;

      for (int j = bin1; j<=bin2; j++)
		{
			data+=datahist->GetBinContent(j);
			bg+=fithist->GetBinContent(j);
			double err1 = -totlow->GetBinContent(j);
			double err2 = tothigh->GetBinContent(j);
			err+=TMath::Max(err1,err2)*bg; //why multiply by bg???
		}
      cout << "Total Data/Bg+/-err in mass window[" << mass << "] = "<< data <<"/ "<< bg << "+/-" << err <<endl;
      double prob = 0;
      fun1->SetParameter(0,bg);
      fun1->SetParameter(1,err);
      
		for (int j = int(data+0.001); j<100; j++){
			fun1->SetParameter(2,j);
			//fun1->Print();
			//cout << "Evaluating Intrgral for j = " << j << " from x0= " << TMath::Max(0.0,bg-10*err) << " to x1 = " << bg+10*err << endl;
			double val = fun1->Integral(TMath::Max(0.0,bg-10*err),bg+10*err);
			//double val = fun1->Integral(TMath::Max(0.0,bg-2*err),bg+2*err);
			/*for (int z=TMath::Max(0.0,bg-2*err); z < bg+2*err; ++z)
			{
				if (c<4)
				{
				std::cout << "func at [" << z << "]=" << 	fun1->Eval(z) << std::endl;
				}
			}
			*/
			prob += val;
      }
      cout<< "Prob for mass[" << mass<<"]="<< prob <<endl;
      p->SetBinContent(p->FindBin(mass),prob);
    }
  }
  /*
  delete gRandom;
  gRandom = (TRandom*) new TRandom3;
  gRandom->SetSeed(3);

  TH1D *minp = new TH1D("minp","Minimum Probability of Each PseudoExpt",100,0,0.2);

  //int nexp = 50000;
  int nexp = 10;

  TH1D *htemp = (TH1D*)datahist->Clone("htemp");

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
	if (prob<minprob) minprob=prob;
      }
    }
    minp->Fill(minprob);
  }
  */
  TCanvas *c1 = new TCanvas("c1","c1");
  TH2D *fr = new TH2D("fr","",100,150,900,100,1e-5,2);
  fr->SetStats(0);
  fr->SetXTitle("M(#gamma,lead jet)(GeV)");
  fr->SetYTitle("Prob of fluctuation #geq N_{obs}");
  fr->GetXaxis()->CenterTitle();
  fr->GetYaxis()->CenterTitle();
  fr->DrawCopy();
  p->SetLineWidth(2);
  p->DrawCopy("same");
  double minp=0;
  double mgg=0;
  double minc = 10;
  for (int i = 1; i<=p->GetNbinsX(); i++){
    double bin = p->GetBinCenter(i);
    double binc = p->GetBinContent(i);
    if (binc<minc){
      minp = binc;
      mgg = bin;
      minc = binc;
    }
  }
  cout<<mgg<<" "<<minp<<endl;
  gPad->SetLogy();
  double p1s = 0.00458319;
  double m1s = 0.0435982;
  double s3s = 0.000100319;
  TLine *l1 = new TLine(150,p1s,900,p1s);
  TLine *l2 = new TLine(150,m1s,900,m1s);
  TLine *l3 = new TLine(150,s3s,900,s3s);
  l1->SetLineColor(4);
  l2->SetLineColor(4);
  l3->SetLineColor(2);
  l1->SetLineWidth(2);
  l2->SetLineWidth(2);
  l3->SetLineWidth(2);
  l1->Draw();
  l2->Draw();
  l3->Draw();
  TLatex *t1 = new TLatex(250,m1s/4,"Expected Range for Min. Obs. Prob.");
  t1->SetTextColor(4);
  t1->SetTextSize(0.05);
  t1->Draw();
  TLatex *t2 = new TLatex(350,s3s*1.5,"3 #sigma evidence level");
  t2->SetTextColor(2);
  t2->SetTextSize(0.05);
  t2->Draw();

  TLatex *t3 = new TLatex(0.3,0.93,"CDF Run II Preliminary, 2.0 fb^{-1}");
  t3->SetNDC(true);
  t3->SetTextSize(0.06);
  t3->Draw();

//
//  TCanvas *c2 = new TCanvas("c2","c2");
//  minp->DrawCopy();
  //cout<<minp->GetMean()<<endl;
}
