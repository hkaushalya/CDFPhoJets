#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"

//study the normalizing of hist with var bin sizes.


void NormH2(TH1* hist)
{

	assert (hist != NULL && "Dump:: hist is null");
	hist->Print();
	for (int i=1;i<=hist->GetNbinsX(); ++i)
	{
			double val = hist->GetBinContent(i)/hist->GetBinWidth(i);
			hist->SetBinContent(i, val);
	}
}


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

void VarBinHist()
{

	TH1F *h1 = new TH1F("fixedbin","fixed bin size hist",20,0,20);

	for (int i=0;i<20;++i){
		h1->Fill(i);
		h1->Fill(i);
	}

	Dump(h1);
	new TCanvas();
	h1->Draw("E");


	Float_t bin[]={0,1,3,6,10};
	TH1F *h2 = new TH1F("varbin","variable bin size hist",4,bin);

	for (int newbin = 1; newbin <=h2->GetNbinsX(); ++newbin)
	{
		for (int oldbin = 1; oldbin <=h1->GetNbinsX(); ++oldbin)
		{
			std::cout << "newbin " << newbin << " [" << h2->GetXaxis()->GetBinLowEdge(newbin) << ", " << 
					h2->GetXaxis()->GetBinUpEdge(newbin) << "] oldbin " << oldbin <<" center[" <<
					h1->GetBinCenter(oldbin) << "]" << std::endl;
			if ( (h1->GetBinCenter(oldbin) >= h2->GetXaxis()->GetBinLowEdge(newbin))
				 && (h1->GetBinCenter(oldbin) < h2->GetXaxis()->GetBinUpEdge(newbin)) )
			{
				double val = h1->GetBinContent(oldbin) + h2->GetBinContent(newbin);
				h2->SetBinContent(newbin, val);

			}

		}
	}


	NormH2(h2);
	
	Dump(h2);
	new TCanvas();
	h2->Draw();
}
