//CVS Log 
/*{{{*/
/************************************************************************
 * $Id: scan.C,v 1.3 2011/05/26 17:40:20 samantha Exp $
 * $Log: scan.C,v $
 * Revision 1.3  2011/05/26 17:40:20  samantha
 * Revamped to fit my needs.
 *
 * Revision 1.2  2010/04/30 23:32:25  samantha
 * MODIFIED: For my work. I have added many comments about the code and the
 * calculations. See elog#1625 for details.
 * CHANGED: 1. The probability function is modified to use TMath functions which
 * can handle large statistics.
 * BUGFIX: 1 Minor bug is fixed replacing 'err+=TMath::Max(err1,err2)*bg' with
 * 'err+=TMath::Max(err1,err2)*fithist->GetBinContent(j)'
 *
 * This error calculation is correct only if the systematic error hists contains
 * fractional errors. Since I am saving absolute values I need to either modify the
 * input or the code to account for this.
 *
 ************************************************************************/
/*}}}*/

/* This is for scanning mass distributions.
 * Require: 1. data hist
 * 			2. background hist
 * 			3. +sigma systematic hist 
 *  			4. -sigma systematic hist
 *		The systematic hists should contain fractional error.
 *		Not absolute.
 */
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
#include "TMath.h"
#include "TPad.h"

using namespace std;

const double pi = 3.14159265;

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
	//cout << "f/ff = " << f  << "/" << ff << endl;


 // return f;
  return ff;
}

void scan(){

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
	TCanvas *c = new TCanvas();
	datahist->DrawCopy();
	c->SetEditable(0);

/*	if (datahist->IsZombie()) { std::cout << "Err! datahist not found " << std::endl; return; }
	if (fithist->IsZombie()) { std::cout << "Err! err_hist not found " << std::endl; return; }
	if (tothigh->IsZombie()) { std::cout << "Err! syst_plus hist not found " << std::endl; return; }
	if (totlow->IsZombie()) { std::cout << "Err! syst_minus hist not found " << std::endl; return; }
*/

  TH1D *p = new TH1D("p","Probability",100,150,1000);

  for (int i = 1; i<=fithist->GetNbinsX(); i++)
  {
    double mass = fithist->GetBinCenter(i);
  
	 if (mass>100 && mass<1000)
	 {
      double sigma = sqrt(2.0)*sqrt(pow(0.135*sqrt(mass/2.0),2)+
				    pow(0.02*(mass/2.0),2));	//this was taken from the previous RS graviton analysis note for electrons/photons. which should work for pho+jets too.
														//this is about ~2GeV. I am not sure about the M/2 in the constant term. (see elog#1625)
										
		//const double mysigma = 0.135 * sqrt(mass) + 0.014 * mass;
      //cout<< "Mass +- sigma = " << mass<<"+-"<<sigma << " (mysigma = " << mysigma << ")" <<endl;
		
      cout<< "Mass +- sigma = " << mass<<"+-"<<sigma  <<endl;
      int bin1 = fithist->FindBin(mass-sigma/2);
      int bin2 = fithist->FindBin(mass+sigma/2);
      cout<<"\t bin [" << bin1 << "-" << bin2 << "] low , hi edges ["<< fithist->GetBinLowEdge(bin1) <<", "<< fithist->GetXaxis()->GetBinUpEdge(bin2) << "]" <<endl;
      double data = 0;
      double bg = 0;
      double err = 0;

      for (int j = bin1; j<=bin2; j++)
		{
			data+=datahist->GetBinContent(j);
			bg+=fithist->GetBinContent(j);
			double err1 = -totlow->GetBinContent(j);
			double err2 = tothigh->GetBinContent(j);
			//err+=TMath::Max(err1,err2)*fithist->GetBinContent(j);
			err = TMath::Max(err1,err2); //as I am passing in absolute errors - sam - 05-03-2010
		}
      cout<<"\t data / bg+-err = "<< data<<" / "<<bg<<" +- "<<err<<endl;
      double prob = 0;
      fun1->SetParameter(0,bg);
      fun1->SetParameter(1,err);
      
		//for (int j = int(data+0.001); j<100; j++){	//this 100 should be changed to fit my stat region. if my data events =1000 and sigma ~30. So I should go up 5*sigma
		//for (int j = int(data+0.001); j< data + sqrt(data) * 5 ; j++){	//this 100 should be changed to fit my stat region. if my data events =1000 and sigma ~30. So I should go up 5*sigma
		for (int j = int(data+0.001); j< data + sqrt(data) * 100 ; j++){	//this 100 should be changed to fit my stat region. if my data events =1000 and sigma ~30. So I should go up 5*sigma
			fun1->SetParameter(2,j);
			prob += fun1->Integral(TMath::Max(0.0,bg-10*err),bg+10*err);	//the factor of 10 is 10*sigma, the integration region. ideally we want to go from -inf to + inf. but this should be enough.
			//std::cout << "\t\t inst prob = " << prob << endl;
      }
      cout<<"\t Prob = "<< prob <<endl;
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
  //TH2D *fr = new TH2D("fr","",100,150,650,100,1e-13,2);
  TH2D *fr = new TH2D("fr","",100,150,650,100,1e-5,2);
  fr->SetStats(0);
  fr->SetXTitle("M(#gamma,lead jet)(GeV)");
  fr->SetYTitle("Prob of fluctuation #geq N_{obs}");
  fr->GetXaxis()->CenterTitle();
  fr->GetYaxis()->CenterTitle();
  fr->DrawCopy();
  p->SetLineWidth(2);
  p->DrawCopy("same");
  double minp;
  double mgg;
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
  TLine *l1 = new TLine(150,p1s,650,p1s);
  TLine *l2 = new TLine(150,m1s,650,m1s);
  TLine *l3 = new TLine(150,s3s,650,s3s);
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
