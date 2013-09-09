#include<iostream>
#include "TH1.h"
#include "TCanvas.h"

//testing if the TH1::Integral() include sthe overflow bin or not
void integral()
{

	//TH1F *h = new TH1F("H1","Test Integral",5,0,5);
	TH1F *h = new TH1F("H1","Test Integral",10,-1,1);
	h->Sumw2();
	//h->Sumw2();

/*	h->Fill(-2);
	h->Fill(-2);
	h->Fill(3);
	h->Fill(1);
	h->Fill(2);
	h->Fill(2);
	*/
	h->FillRandom("gaus",10);
	
	h->DrawClone("E");

	h->Print();
	std::cout << "nbins, Entries, Integral = " << h->GetNbinsX() << ", " 
				<< h->GetEntries() << ", " << h->Integral() << std::endl;

	for (int i = 1; i <= h->GetNbinsX()+1; ++i)
	{
		std::cout << "bin, content, err, err^2 = " << i  << "\t" 
			<< h->GetBinContent(i) << "\t" 
			<< h->GetBinError(i) << "\t" 
			<< pow(h->GetBinError(i),2) <<std::endl;
	}



	//TH1F *h2 = new TH1F("H2","Test Integral2",5,0,5);
	TH1F *h2 = new TH1F("H2","Test Integral2",10,-1,1);
	h2->Sumw2();
/*	h2->Fill(2);
	h2->Fill(1);
	h2->Fill(3);
	h2->Fill(2);
	*/
	h2->FillRandom("pol2",50);
	new TCanvas();
	h2->DrawClone("E");
	
	h2->Print();
	std::cout << "nbins, Entries, Integral = " << h2->GetNbinsX() << ", " 
				<< h2->GetEntries() << ", " << h2->Integral() << std::endl;

	for (int i = 1; i <= h2->GetNbinsX()+1; ++i)
	{
		std::cout << "bin,content, err, err^2, rel.err = " << i  << "\t" 
			<< h2->GetBinContent(i) << "\t" 
			<< h2->GetBinError(i) << "\t" 
			<< pow(h2->GetBinError(i),2) << "\t" 
			<< h2->GetBinError(i)/h2->GetBinContent(i)
			<<std::endl;
	}

	//h2->Add(h);
	TH1* h2_copy = dynamic_cast<TH1*> (h2->Clone("h2copy"));
	TH1* h_copy = dynamic_cast<TH1*> (h->Clone("hcopy"));

	h2->Print();
	h2->Scale(100.0/(1.0 * h2->Integral()));
	std::cout << "nbins, Entries, Integral = " << h2->GetNbinsX() << ", " 
				<< h2->GetEntries() << ", " << h2->Integral() << std::endl;

	for (int i = 0; i <= h2->GetNbinsX()+1; ++i)
	{
		std::cout << "bin, content, err, err^2, rel.err = " << i  << "\t" 
			<< h2->GetBinContent(i) << "\t" 
			<< h2->GetBinError(i) << "\t" 
			<< pow(h2->GetBinError(i),2) << "\t"
			<< h2->GetBinError(i)/h2->GetBinContent(i)
			<<std::endl;
	}


	h2_copy->Scale(10.0/(1.0 * h2_copy->Integral()));


	h->Scale(h2->Integral()/(double) h->Integral());
	h_copy->Scale(h2_copy->Integral()/(double) h2->Integral());

	h->Divide(h2);
	h_copy->Divide(h2_copy);

	new TCanvas();
	h->Draw();
	new TCanvas();
	h_copy->Draw();

	
};

