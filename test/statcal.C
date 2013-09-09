#include<iostream>
#include "TH1.h"
#include "TCanvas.h"

//testing if the TH1::Integral() include sthe overflow bin or not
void statcal()
{

	TH1F *h = new TH1F("H1","Test Integral",5,0,5);
	h->Sumw2();

	h->Fill(-2);
	h->Fill(-2);
	h->Fill(3);
	h->Fill(1);
	h->Fill(2);
	h->Fill(2);
	
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



	TH1F *h2 = new TH1F("H2","Test Integral2",5,0,5);
	h2->Sumw2();
	h2->Fill(2);
	h2->Fill(1);
	h2->Fill(3);
	h2->Fill(2);
	new TCanvas();
	h2->DrawClone("E");
	
	h2->Print();
	std::cout << "nbins, Entries, Integral = " << h2->GetNbinsX() << ", " 
				<< h2->GetEntries() << ", " << h2->Integral() << std::endl;

	for (int i = 1; i <= h2->GetNbinsX()+1; ++i)
	{
		if (h2->GetBinContent(i)) 
		std::cout << "bin,content, err, err^2, rel.err = " << i  << "\t" 
			<< h2->GetBinContent(i) << "\t" 
			<< h2->GetBinError(i) << "\t" 
			<< pow(h2->GetBinError(i),2) << "\t" 
			<< h2->GetBinError(i)/h2->GetBinContent(i)
			<<std::endl;
	}

	h2->Print();
	h2->Scale(100/(1.0 * h2->Integral()));
	new TCanvas();
	h2->DrawClone("E");
	std::cout << "nbins, Entries, Integral = " << h2->GetNbinsX() << ", " 
				<< h2->GetEntries() << ", " << h2->Integral() << std::endl;

	for (int i = 0; i <= h2->GetNbinsX()+1; ++i)
	{
		if (h2->GetBinContent(i)) 
		std::cout << "bin, content, err, err^2, rel.err = " << i  << "\t" 
			<< h2->GetBinContent(i) << "\t" 
			<< h2->GetBinError(i) << "\t" 
			<< pow(h2->GetBinError(i),2) << "\t"
			<< h2->GetBinError(i)/(double)h2->GetBinContent(i)
			<<std::endl;
	}

	
};

