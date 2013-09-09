#include<iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "../RootTools/CommonTools.hh"
#include "TSystem.h"
#include "TStyle.h"

//testing how the TH1::Integral() values change when a hist is scaled
void scale()
{
	gSystem->Load("../RootTools/CommonTools.C.so");
	gStyle->SetOptStat("nemouri");
	
	TH1F *h = new TH1F("H1","Test Integral",100,0,10);
	h->Fill (1);
	h->Fill (1);
	h->Fill (1);
	h->Fill (1);
	//TH1F *h2 = new TH1F("H2","Test Integral copy",100,0,1);
	//h->FillRandom("gaus",10000);
	//h2->FillRandom("gaus",10000);

//	std::cout << "nbins, Entries, Integral = " << h->GetNbinsX() << ", " 
//				<< h->GetEntries() << ", " << h->Integral("width") << std::endl;
	//h = (TH1F*) MakeVariableBinHist(h, 0, 0.2, 0.6, 0.8, 1, 0.1, 0.2, 0.2, 0.2);
	//std::cout << "nbins, Entries, Integral = " << h->GetNbinsX() << ", " 
	//			<< h->GetEntries() << ", " << h->Integral("width") << std::endl;

	new TCanvas();
	h->Draw();
	//new TCanvas();
	//h2->Draw();

	//h->Scale(1,"width"); //this replaces the loop below
	/*
	for (int bin = 1; bin <= h->GetNbinsX(); ++bin)
	{
		const double val = h->GetBinContent(bin);
		const double err = h->GetBinError(bin);
		const double width = h->GetBinWidth(bin);
		h->SetBinContent(bin,val/width);
		h->SetBinError(bin,err/width);
	}
	*/
	//h->DrawClone("E");
	h->Print();
	
	std::cout << "nbins, Entries, Integral = " << h->GetNbinsX() << ", " 
				<< h->GetEntries() << ", " << h->Integral("width") << std::endl;


	//h->Scale(10);
	h->Scale(10/(double)h->Integral());
	//h_copy->Scale(h2_copy->Integral()/(double) h2->Integral());
	//h_copy->Draw();

	std::cout << "after scaling to " << 10 << std::endl;
	h->Print();
	std::cout << "nbins, Entries, Integral = " << h->GetNbinsX() << ", " 
				<< h->GetEntries() << ", " << h->Integral("width") << std::endl;
	
};

