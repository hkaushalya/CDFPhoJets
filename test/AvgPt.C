#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"

//study the normalizing of hist with var bin sizes.


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

void AvgPt()
{

	TH1F *h1 = new TH1F("fixedbin","fixed bin size hist",20,0,20);
		h1->Fill(1);
		h1->Fill(1);
		h1->Fill(1);
		h1->Fill(18);
		h1->Fill(18);
		h1->Fill(18);
		h1->Fill(19);
		h1->Fill(20);
		h1->Fill(20);
		h1->Fill(20);
		h1->Fill(20);
		h1->Fill(20);

		

	new TCanvas();
	h1->Draw("E");

	for (int bin=1;bin<h1->GetNbinsX(); ++bin)
		std::cout << "bin, loedge = " << bin << "\t" << h1->GetBinLowEdge(bin) << "\t" <<  h1->GetBinContent(bin) << std::endl;

	//const Float_t bin[]={0,1,3,6,10};
	//TH1F *h2 = new TH1F("varbin","variable bin size hist",4,bin);
	//h2->Fill(3);
	//h2->Fill(4);
	//h2->Fill(4);

}
