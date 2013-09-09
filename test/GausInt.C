#include<iostream>
#include"TH1F.h"
#include "TF1.h"
#include"TMath.h"
#include"TLegend.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TRandom1.h"

void GausInt(const unsigned int pts =1000)
{
	TH1F* hist = new TH1F("HIST","GAUS RANDOM",100,-5,5);

	hist->FillRandom("landau",pts);
	hist->Draw();
	std::cout << "Integral = " << hist->Integral() << std::endl;
	
	TF1 *f1 = new TF1("fit","landau",-5,5);
	hist->Fit(f1);
	std::cout << f1->GetParameter(0) << std::endl;
	std::cout << f1->GetParameter(1) << std::endl;
	std::cout << f1->GetParameter(2) << std::endl;

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
std::cout << "x, val = " << x[0] << ", " << fitval << std::endl;
  return fitval;
}

//gaus part
double myFunc1_f1(double *x, double *par)
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
  //double f1 = normL / (1.0 + normL) * TMath::Gaus(arg,mean,sigmaG);
  double f1 = TMath::Gaus(x[0],mean,sigmaG);
  //fitval = par[0] * f1;
  fitval = f1;
  return fitval;
}

//landau part
double myFunc1_f2(double *x, double *par)
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
  double f2 = TMath::Landau(arg_L,mpv,sigmaL) / (1.0+normL);
  //fitval = par[0] * f2;
  fitval = f2;
  return fitval;
}




void DrawmyFunc1(const float norm=0.9, const float Gmean=1.5, const float Gsigma=0.7,
						const float mpv=1.5, const float Lsigma=0.5, const float Lnorm=0.5)
{
	TLegend *leg = new TLegend (0.5,0.8,0.90,0.90);
	leg->SetTextFont(42);
	leg->SetTextSize(0.025);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);

	gStyle->SetAxisColor(kBlue);
	gStyle->SetPadTickY(1);
	gStyle->SetPadTickX(1);
	gStyle->SetOptFit(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetTitleFontSize(0.035);
	gStyle->SetStatFontSize(0.03);
	gStyle->SetTitleFont(2);

	
	TF1 *fitFun1 = new TF1("fit1",myFunc1,-2,5,6);
	fitFun1->SetLineWidth(2);
	fitFun1->SetParameter(0,norm);			// norm 
	fitFun1->SetParameter(1,Gmean);		// mean of Gauss
	fitFun1->SetParameter(2,Gsigma);	// Sigma of Gauss
	fitFun1->SetParameter(3,mpv);		// MPV of Landau
	fitFun1->SetParameter(4,Lsigma);	// Sigma of Landau
	fitFun1->SetParameter(5,Lnorm);					// Landau norm

	fitFun1->Draw();

	TF1 *fitFun1_f1 = new TF1("fit1_f1",myFunc1_f1,-2,5,6);
	fitFun1_f1->SetLineColor(kRed);
	fitFun1_f1->SetLineWidth(2);
	fitFun1_f1->SetParameter(0,norm);			// norm 
	fitFun1_f1->SetParameter(1,Gmean);		// mean of Gauss
	fitFun1_f1->SetParameter(2,Gsigma);	// Sigma of Gauss
	fitFun1_f1->SetParameter(3,mpv);		// MPV of Landau
	fitFun1_f1->SetParameter(4,Lsigma);	// Sigma of Landau
	fitFun1_f1->SetParameter(5,Lnorm);					// Landau norm
	fitFun1_f1->Draw("same");

	TF1 *fitFun1_f2 = new TF1("fit1_f2",myFunc1_f2,-2,5,6);
	fitFun1_f2->SetLineColor(kBlue);
	fitFun1_f2->SetLineWidth(2);
	fitFun1_f2->SetParameter(0,norm);			// norm 
	fitFun1_f2->SetParameter(1,Gmean);		// mean of Gauss
	fitFun1_f2->SetParameter(2,Gsigma);	// Sigma of Gauss
	fitFun1_f2->SetParameter(3,mpv);		// MPV of Landau
	fitFun1_f2->SetParameter(4,Lsigma);	// Sigma of Landau
	fitFun1_f2->SetParameter(5,Lnorm);					// Landau norm
	fitFun1_f2->Draw("same");


	leg->AddEntry(fitFun1, "Gaus+Landau");
	leg->AddEntry(fitFun1_f1, "Gaus");
	leg->AddEntry(fitFun1_f2, "Landau");
	leg->Draw();
}
	

int main(int argc, char* argv[]) {


	TApplication* rootapp = new TApplication("example",&argc, argv);

	TRandom1* myrand = new TRandom1();

	TH1F* myhist = new TH1F("stats","",100,0,10);

	for(int i=0;i<10000;++i) {

		myhist->Fill(myrand->Gaus(5,1));

	}

	myhist->Draw();

	rootapp->Run();

	return 0;

}
