#include<iostream>
#include "TF1.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TList.h"
#include "TIterator.h"

void Weight()
{
		TF1* line = new TF1("LINE","0.02 * x", 0,50);
		TF1* line2 = new TF1("LINE2","100 * x", 0,50);
		TF1* arc = new TF1("ARC","x*x", 0,50);

		TH1F* hline = new TH1F("hline","line hist",50,0,100);
		TH1F* hline2 = new TH1F("hline2","line2 hist",50,0,100);
		TH1F* harc = new TH1F("harc","arc hist",50,0,100);
		hline->Sumw2();
		hline2->Sumw2();
		harc->Sumw2();

		hline->SetLineColor(kRed);
		hline2->SetLineColor(7);
		harc->SetLineColor(kBlue);

		hline->FillRandom("LINE",5000);
		hline2->FillRandom("LINE2",5000);
		harc->FillRandom("ARC",5000);

		new TCanvas();
		gPad->SetGridx();
		gPad->SetGridy();
		hline->DrawClone();
		harc->DrawClone("same");
		//hline2->Draw("same");
		
		/*new TCanvas();
		gPad->SetGridx();
		gPad->SetGridy();
		harc->Draw();
		*/

		new TCanvas();
		gPad->SetGridx();
		gPad->SetGridy();
		TH1F* hline_copy = (TH1F*) hline->Clone("hline_copy");
		hline->Scale(harc->Integral()/(double)hline->Integral());
		hline->Divide(harc);
		hline->Fit("pol1");
		hline->Draw();

		TList *funcList = hline->GetListOfFunctions();
		TIterator *it = funcList->MakeIterator();
		TF1* fitfunc = (TF1*) it->Next();


		TH1F* warc = (TH1F*) harc->Clone("arc_copy");
		warc->Print();

		std::cout << "nbins = " << warc->GetNbinsX() << std::endl; 
		for (int bin = 1; bin <= warc->GetNbinsX(); ++bin)
		{
			std::cout << "bin/val/time/final= " << bin  << ", " 
				<< harc->GetBinContent(bin) 
				<< ", " << hline->GetBinContent(bin) << ", "
				<< harc->GetBinContent(bin) * hline->GetBinContent(bin)
				<< std::endl;
				//warc->SetBinContent(bin,harc->GetBinContent(bin) * hline->GetBinContent(bin));
				warc->SetBinContent(bin,harc->GetBinContent(bin) * fitfunc->Eval(warc->GetBinCenter(bin)));
				warc->SetBinError(bin,0);
		}

		new TCanvas();
		gPad->SetGridx();
		gPad->SetGridy();
		warc->SetMarkerStyle(23);
		warc->DrawClone("P");
		hline_copy->DrawClone("P same");


		warc->Divide(hline_copy);
		new TCanvas();
		gPad->SetGridx();
		gPad->SetGridy();
		warc->DrawClone("P");
		
};
