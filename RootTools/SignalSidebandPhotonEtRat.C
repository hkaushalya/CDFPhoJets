/********************************************************************
 * This is to compare the sideband photon Et to signal photon
 * Et spectrum. This ratio will be used to reweight the sideband
 * events and try to see if it can model the Ht plot
 * better. PYTHIA Photon MC does not describe second leading
 * jet well enough and it also causes the Ht plot to disagree.
 * This will pick up the sideband and signal photon Et plots
 * and normalize sideband to signal and take the ratio.
 * see elog#https://hep.baylor.edu/hep/samantha/1500 for the 
 * final results of this. 
 * Author: Sam Hewamanage <samantha@fnal.gov> 12-07-2009
 ********************************************************************
 *
 * $Id: SignalSidebandPhotonEtRat.C,v 1.1 2013/02/28 03:27:20 samantha Exp $
 * $Log: SignalSidebandPhotonEtRat.C,v $
 * Revision 1.1  2013/02/28 03:27:20  samantha
 * Final commit. no checks. These were never commited to cvs.
 *
 *
 *******************************************************************/
#include <iostream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLegend.h"
#include <sstream>
#include <iomanip>
#include <cmath>
#include "TF1.h"
#include "TFile.h"
#include "CommonTools.hh"
#include "TMath.h"


Double_t fitFunction(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
		double pol1 = par[0]  + par[1] * x[0] * x[0]; 
		double pol2 = par[2]  + par[3] * x[0];

		double wgt = 1 / (1 + exp((x[0]-par[4])/10.5));
		val = pol1 * wgt + pol2 * (1 - wgt);
	}
	return val;
}
void SignalSidebandPhotonEtRat(
					const std::string title,
				 	const float xmin, const float xpoint1, 
					const float xpoint2, const float xpoint3, 
					const float xpoint4, const float width1, 
					const float width2, const float width3, 
				const float width4, 
					const std::string brname
		)
{



	
	//TFile *file = new TFile("PhoJets_data.root");
	//TFile *file = new TFile("SidebandPhoReweightForHtPlot.root");
	//TFile *file = new TFile("../PhoFlat/SidebandPhoReweightForHtPlot.root");
	//TFile *file = new TFile("../PhoFlat/PhoJets_phodata.root");
	TFile *file = new TFile("../PhoFlat/SidebandPhoReweightForHtPlot.root");
	assert (file != NULL && "file not found!");

	TH1* hsig = dynamic_cast<TH1F*> (file->Get("Hist/SIGNAL/1Jet/Photon/EtCorr"));	
	TH1* hside = dynamic_cast<TH1F*> (file->Get("Hist/SIDEBAND/1Jet/Photon/EtCorr"));	
	assert (hsig != NULL && "signal hist is null");
	assert (hside != NULL && "sideband hist is null");


  	hsig = (TH1F*) MakeVariableBinHist (hsig, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
  	hside = (TH1F*) MakeVariableBinHist (hside, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
	hsig->SetLineColor(kRed);
	hside->SetLineColor(kBlue);
	
	hsig->Print();
	hside->Print();
	//sumw2 is already established in the above hist rebinning
	hside->Scale(hsig->Integral()/(double)hside->Integral());
	hside->Print();
	hsig->Divide(hside);
	new TCanvas();
	gPad->SetGridx();
	gPad->SetGridy();
	hsig->SetTitle(title.c_str());
	//gPad->SetLogy();
	hsig->DrawClone("PE");
	//hside->DrawClone("SAME PE");

	return;
	
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(0);
	hsig->SetMarkerColor(kRed);
	hsig->SetLineColor(kBlue);
	hsig->SetLineWidth(2);
	hsig->SetMarkerStyle(22);
	hsig->GetXaxis()->CenterTitle(kTRUE);
	hsig->GetYaxis()->CenterTitle(kTRUE);
	hsig->SetTitle(title.c_str());
	TCanvas *c1 = new TCanvas();
	c1->SetName("Fit");
	gPad->SetGridx();
	gPad->SetGridy();
	//find limits for fitting range to improve fit
	int xminbin=0,xmaxbin=0;
	FindNonZeroXrange(hsig, xminbin,  xmaxbin);
	
	
	hsig->SetMaximum(20);
	hsig->SetMinimum(0);
	//hsig->Fit("pol3","","",hsig->GetBinCenter(xminbin), hsig->GetBinCenter(xmaxbin));
	//TF1 *f1 = new TF1("fit_func","[0]*x*x*x",hsig->GetBinCenter(xminbin),hsig->GetXaxis()->GetBinUpEdge(xmaxbin-3));
	//TF1 *f2 = new TF1("fit_func2","pol1",hsig->GetBinCenter(xminbin+3),hsig->GetXaxis()->GetBinUpEdge(xmaxbin));
	//TF1 *f3 = new TF1("fit_func3","pol2",hsig->GetBinCenter(xminbin),hsig->GetXaxis()->GetBinUpEdge(xminbin+4));
	TF1 *f4 = new TF1("fit_func4",fitFunction,hsig->GetBinCenter(xminbin),hsig->GetXaxis()->GetBinUpEdge(xmaxbin),5);
	//TF1 *f4 = new TF1("fit_func4",fitFunction,0,500,5);
	f4->SetParameter(0,0.5);
	f4->SetParameter(1,0.5);
	f4->SetParameter(2,2);
	f4->SetParameter(3,1.5);
	f4->SetParameter(4,50);
	f4->SetParLimits(4,30,35);
	f4->SetLineColor(8);
	
	//hsig->Fit(f1,"+R");
	//hsig->Fit(f2,"+R");
	//hsig->Fit(f3,"+R");
	hsig->Fit(f4,"+R");
	//hsig->Draw("AP");
	hsig->Draw("P");

	TList *funcList = hsig->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1* fitfunc = (TF1*) it->Next();
	std::string func_name("tf1_Pho_signal_side_ratio");
	fitfunc->SetName(func_name.c_str());

	TFile f("PhoEt_SigSideRatio_FitFunc.root","UPDATE");
	if ( f.Get(func_name.c_str()) != NULL)
	{
		std::stringstream old_c1, old_func;
		old_c1 << c1->GetName() << ";1";
		f.Delete(old_c1.str().c_str());
		old_func << fitfunc->GetName() << ";1";
		f.Delete(old_func.str().c_str());
	}
	
	fitfunc->Write();
	c1->Write();
	f.ls();
	f.Close();
	

}

void SignalSidebandPhotonEtRat()
{
	SignalSidebandPhotonEtRat ("#gamma+ #geq 1Jet;E_{T}^{#gamma} (GeV); E_{T}^{#gamma-signal}/E_{T}^{#gamma-sideband} (GeV);", 
			//0,100,180,300,400, 10,20,120,100, //these are the bin sized I use in the final plots.
			0,100,180,300,400, 10,10,120,100, //these are the bin sized I use in the final plots.
			//0,100,180,500,650, 20,20,320,150, //these are the bin sized I use in the final plots.
			"pj1_Et_pho");
}
