#include <iostream>
#include "TF1.h"
#include "TCanvas.h"

double MyJER(double* x, double* par)
{
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




void FitFunc1(const int badone = 1)
{
	TF1* JerFun=new TF1("jer",MyJER,-1.0,5.,5);

	double p0=0,p1=0,p2=0,p3=0,p4=0;
	if (badone)
	{
		std::cout << "Drawing the BAD one!" << std::endl;
		//breaks for this
		p0 = -0.0191989;
		p1 = 0.40405;
		p2 = -0.239216;
		p3 = 0.000724924;
		p4 = 0.28215;
	} else {
		std::cout << "Drawing the GOOD one!" << std::endl;
		p0 = -0.0613734;
		p1 = 0.282392;
		p2 = -0.220161;
		p3 = 0.0816814;
		p4 = 0.37754;

	}

	JerFun->SetParameter(0,p0);
	JerFun->SetParameter(1,p1);
	JerFun->SetParameter(2,p2);
	JerFun->SetParameter(3,p3);
	JerFun->SetParameter(4,p4);

	new TCanvas();
	//TF1* Landau=new TF1("Landau_func","TMath::Landau(x,[0],[1],0)",-1.0,5.);
	TF1* Landau=new TF1("Landau_func","landau",-1.0,5.);
	//Landau->SetParameters(p2,p3);
	Landau->SetParameter(0,0.5);
	Landau->SetParameter(1,0);
	Landau->SetParameter(2,0.5);
	Landau->Print();
	Landau->SetLineColor(kRed);
	Landau->Draw();

//	TF1* Gauss=new TF1("Gauss","TMath::Gaus(x,[0],[1],0)",-1.0,5.);
//	Gauss->SetParameters(p0,p1);
//	Gauss->SetLineColor(kBlue);
	//Gauss->Draw("same");
//	new TCanvas();
//	Gauss->Draw();

	JerFun->SetTitle("JER before jet before the problem");
	//JerFun->Draw("SAME");


	TF1* myf = new TF1("myf","1/x",10,20);
	
	for (float i=19900; i < 100000; ++i) std::cout << i << " -> " << myf->Eval(i) << std::endl;
	

}
