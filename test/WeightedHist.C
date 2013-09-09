#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"

//study the weighing of hists.


void Dump(TH1* hist)
{
	assert (hist != NULL && "Dump:: hist is null");
	hist->Print();
	std::cout << "bin / edges / val / err " << std::endl;
	for (int i=1;i<=hist->GetNbinsX(); ++i)
	{
		std::cout << i << " / " << hist->GetBinLowEdge(i) << ", " 
			<< hist->GetXaxis()->GetBinUpEdge(i) << " / " << hist->GetBinContent(i) 
			<< ", " << hist->GetBinError(i) << std::endl;
	}

}

void WeightedHist()
{

	float wgt[3]={0.3,0.5,0.2};
	TH1F *h1 = new TH1F("hist1","fixed bin size hist",20,0,20);
	TH1F *h2 = new TH1F("hist2","fixed bin size hist",20,0,20);
		h1->Fill(1);
		h1->Fill(1);
		h1->Fill(1);
		h1->Fill(2,wgt[0]);
		h1->Fill(2,wgt[1]);
		h1->Fill(2,wgt[2]);
/*		h1->Fill(1);
		h1->Fill(18);
		h1->Fill(18);
		h1->Fill(18);
		h1->Fill(19);
		h1->Fill(20);
		h1->Fill(20);
		h1->Fill(20);
		h1->Fill(20);
		h1->Fill(20);
*/
	
		
	for (int bin=1;bin<=h1->GetNbinsX(); ++bin)
	{
		for (int cont = 0; cont< h1->GetBinContent(bin); ++cont)
		{	
			float val = h1->GetBinLowEdge(bin);
			h2->Fill(val,wgt[cont]);
		}
		
	}

	new TCanvas();
	h1->Draw("E");
	new TCanvas();
	h2->Draw("E");

	for (int bin=1;bin<h1->GetNbinsX(); ++bin)
		std::cout << "bin, loedge = " << bin << "\t" << h1->GetBinLowEdge(bin) << "\t" <<  h1->GetBinContent(bin) << std::endl;

	//const Float_t bin[]={0,1,3,6,10};
	//TH1F *h2 = new TH1F("varbin","variable bin size hist",4,bin);
	//h2->Fill(3);
	//h2->Fill(4);
	//h2->Fill(4);

}
