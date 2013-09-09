#include<iostream>
#include<string>
#include<TFile.h>
#include<TH1.h>
#include<TCanvas.h>
#include<TLegend.h>

void PhxPlots(const std::string name, const std::string path, const int rebin=1)
{
	TFile *pb4 = new TFile("../PhoFlat/phx_b4.root");
	TFile *pa4 = new TFile("../PhoFlat/phx_a4.root");

	if (pb4->IsZombie() || pa4->IsZombie())
	{
		std::cout <<  "One of the files not found" << std::endl;
	}

	pb4->cd(path.c_str());
	TH1* tpho = (TH1*) gDirectory->FindObjectAny(name.c_str());
	assert (tpho != NULL && "tpho is NULL");

	pa4->cd(path.c_str());
	TH1* spho = (TH1*) gDirectory->FindObjectAny(name.c_str());
	assert (spho != NULL && "tpho is NULL");

	//tpho->Scale(1.0/tpho->Integral());
	//spho->Scale(1.0/spho->Integral());

	tpho->Rebin(rebin);
	spho->Rebin(rebin);


	new TCanvas();
	gPad->SetLogy();
	tpho->Draw();
	new TCanvas();
	gPad->SetLogy();
	spho->Draw();
	new TCanvas();
	TH1* ntpho = (TH1*) tpho->Clone("tpho_copy");
	TH1* ratio = (TH1*) tpho->Clone("tpho_copy2");

	ntpho->SetMaximum(1);
	ntpho->SetMinimum(-1);
	ntpho->SetMarkerStyle(kMultiply);
	gPad->SetGridx();
	gPad->SetGridy();
	//ntpho->Draw("P");


	new TCanvas();
	//ntpho->Divide(spho);
	ntpho->Draw("P");
	for (unsigned i=1;i <=tpho->GetNbinsX(); ++i)
	{
		int val = spho->GetBinContent(i);
		std::cout << "bin, val=" << i << "\t" << val << std::endl;
		ratio->SetBinContent(i, val ? val : 0);
	}

	new TCanvas();
	ratio->Draw("P");
}


void PhxPlots()
{
	PhxPlots("Ht", "/Hist/SIGNAL/1Jet/Event/", 5);
}
