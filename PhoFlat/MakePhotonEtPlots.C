#include<iostream>
#include<TFile.h>
#include<TLegend.h>
#include<TH1F.h>
#include<TCanvas.h>
#include<TStyle.h>


void MakePhotonEtPlots()
{

	TFile *f = new TFile("PhoJets_data_TrigPhoEt.root");
	if (f->IsZombie())
	{
		std::cout << "File not found!" << std::endl;
		exit (1);
	}

	TH1F* Et25;
	Et25 = (TH1F*) f->FindObjectAny("lphoEt_Iso25");
	if (!Et25) std::cout << "Et25 hist not found!" << std::endl;
	TH1F* Et50 = (TH1F*) f->FindObjectAny("lphoEt_50");
	if (!Et50) std::cout << "Et50 hist not found!" << std::endl;
	TH1F* Et70 = (TH1F*) f->FindObjectAny("lphoEt_70");
	if (!Et70) std::cout << "Et70 hist not found!" << std::endl;

	Et25->SetLineColor(kRed);
	Et50->SetLineColor(kBlue);

	int rebin=4;
	Et25->Rebin(rebin);
	Et50->Rebin(rebin);
	Et70->Rebin(rebin);

	TLegend *leg = new TLegend (0.6,0.7,0.9,0.9);
	leg->AddEntry(Et25,"pho_iso25");
	leg->AddEntry(Et50,"pho_50");
	leg->AddEntry(Et70,"pho_70");

	Et25->SetTitle("Loose Photon Et by trigger");

	new TCanvas();
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	
	Et25->Draw();
	Et50->Draw("SAME");
	Et70->Draw("SAME");
	leg->Draw();
}
