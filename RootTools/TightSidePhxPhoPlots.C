#include<iostream>
#include<string>
#include<TFile.h>
#include<TH1.h>
#include<TCanvas.h>
#include<TLegend.h>
#include<sstream>
#include<TStyle.h>
#include<TF1.h>
#include<TPaveText.h>

void TightSidePhxPhoPlots(const std::string name, const std::string xtitle, const std::string path1,
						const float xmin, const float xmax, const int rebin=1)
{
	gStyle->SetOptStat("");
	gStyle->SetCanvasColor(10);
	gStyle->SetCanvasBorderMode(0);
		
	std::stringstream sigPath, sidePath, phxPath;
	sigPath << "Hist/SIGNAL/" << path1;
	sidePath << "Hist/SIDEBAND/" << path1;
	phxPath << "Hist/PHOENIX/" << path1;
	TFile *file;
	//file = new TFile("../PhoFlat/TightSidePhx.root");
	file = new TFile("../PhoFlat/Temp.root");

	if (file->IsZombie())
	{
		std::cout <<  "file not found" << std::endl;
		exit (1);
	}

	file->cd(sigPath.str().c_str());
	TH1* sigPho = (TH1*) gDirectory->FindObjectAny(name.c_str());
	assert (sigPho != NULL && "sigPho is NULL");

	file->cd();
	file->cd(sidePath.str().c_str());
	TH1* sidePho = (TH1*) gDirectory->FindObjectAny(name.c_str());
	assert (sidePho != NULL && "sidePho is NULL");

	file->cd();
	file->cd(phxPath.str().c_str());
	TH1* phxPho = (TH1*) gDirectory->FindObjectAny(name.c_str());
	assert (sidePho != NULL && "phxPho is NULL");



	// this method hinders info. need to normalized to lum. 
	// so need to run on the full dataset.
	sigPho->Scale(1.0/sigPho->Integral());
	sidePho->Scale(1.0/sidePho->Integral());
	phxPho->Scale(1.0/phxPho->Integral());

	sigPho->Rebin(rebin);
	sidePho->Rebin(rebin);
	phxPho->Rebin(rebin);



	//new TCanvas();
	//gPad->SetLogy();
	//sigPho->Draw();
	new TCanvas();
	//gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();



	sigPho->SetMarkerStyle(8);
	sigPho->SetMarkerColor(kBlue);
	sidePho->SetMarkerStyle(22);
	sidePho->SetMarkerColor(kRed);
	sidePho->SetLineColor(kRed);
	phxPho->SetMarkerStyle(19);
	phxPho->SetMarkerColor(kGreen);
	phxPho->SetLineColor(kGreen);


	std::stringstream ytitle;

	ytitle << "BinSize = " << sidePho->GetXaxis()->GetBinWidth(1);
	
	//sidePho->SetTitle(title.str().c_str());
	sigPho->SetTitle("");
	sigPho->GetYaxis()->SetTitle(ytitle.str().c_str());
	sigPho->GetXaxis()->SetTitle(xtitle.c_str());
	sigPho->GetXaxis()->CenterTitle(true);
	sigPho->GetYaxis()->CenterTitle(true);
	int labelfont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	int titlefont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	sigPho->GetXaxis()->SetLabelFont(labelfont);
	sigPho->GetYaxis()->SetLabelFont(labelfont);
	sigPho->GetXaxis()->SetTitleFont(titlefont);
	sigPho->GetYaxis()->SetTitleFont(titlefont);
 	sigPho->GetYaxis()->SetLabelSize(0.05);
 	sigPho->GetXaxis()->SetLabelSize(0.05);
	sigPho->GetXaxis()->SetTitleSize(0.05);
	sigPho->GetYaxis()->SetTitleSize(0.05);
	sigPho->GetXaxis()->SetTitleOffset(0.9);
	sigPho->GetYaxis()->SetTitleOffset(0.9);


	sigPho->GetXaxis()->SetRangeUser(xmin, xmax);
	sidePho->GetXaxis()->SetRangeUser(xmin, xmax);
	phxPho->GetXaxis()->SetRangeUser(xmin, xmax);
	sigPho->Draw("P");
	sidePho->Draw("SAME P");
	phxPho->Draw("SAME P");


	gStyle->SetTitleBorderSize(0);

	TLegend *leg = new TLegend (0.6,0.75,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);
 	leg->AddEntry(sigPho,"#gamma_{tight}");
 	leg->AddEntry(sidePho,"#gamma_{sideband}");
 	leg->AddEntry(phxPho,"#gamma_{phoenix}");
	leg->SetFillColor(10);
	leg->Draw();
 	


  	TPaveText *tp = new TPaveText(0.05,0.91,0.98,0.99,"NDC");
	tp->SetLineColor(10);
	tp->SetFillColor(10);
	tp->SetTextFont(41);
	tp->SetTextSize(0.035);
	tp->AddText("Both Signal and Sideband normalized to 1. After Phoenix Rejection (from 1.2M event sample)");
	tp->Draw();


	if (gPad)
	{
		std::stringstream cname, cnamepdf;
		cname << "Phos_TightSidebandPhx_"<< name  << ".gif";
		gPad->Print(cname.str().c_str());
		cnamepdf << "Phos_TightSidebandPhx_"<< name  << ".pdf";
		gPad->Print(cnamepdf.str().c_str(),"pdf");
	}

}


void TightSidePhxPhoPlots()
{
	//TightSidePhxPhoPlots("Ht", "H_{T} [GeV]", "1Jet/Event/",0,1000, 10);
	//TightSidePhxPhoPlots("Sumet", "#Sigma E_{T} [GeV]", "1Jet/Event/"0,600, 10);
	//TightSidePhxPhoPlots("EtCorr","E_{T}^{#gamma} [GeV]", "1Jet/Photon/",0,350, 10);
	//TightSidePhxPhoPlots("NJet15", "NJet", "1Jet/Event/", 0,10,1, 1);
	TightSidePhxPhoPlots("Met", "MET [GeV]", "1Jet/Event/",0,400, 5);
	//TightSidePhxPhoPlots("N12Vertices", "Class 12 Vertices", "1Jet/Event/" 0,15,1, 1);
	//TightSidePhxPhoPlots("Ces2Strip","#gamma 2^{nd} CES STRIP Energy [GeV]", "1Jet/Photon/",0,10, 1);
	//TightSidePhxPhoPlots("Ces2Wire","#gamma 2^{nd} CES WIRE Energy [GeV]", "1Jet/Photon/",0,20, 1);
	//TightSidePhxPhoPlots("Chi2Mean","#gamma Chi^{2} Mean", "1Jet/Photon/",0,25, 1);
	//TightSidePhxPhoPlots("EmTime","#gamma EM Timing [ns]", "1Jet/Photon/",-100,150, 5);
	//TightSidePhxPhoPlots("HadEm","#gamma Had/Em", "1Jet/Photon/",0.0,1.0, 1);
	//TightSidePhxPhoPlots("IsoEtCorr","#gamma IsoEtCorr [GeV]", "1Jet/Photon/",-5,8, 1);
	//TightSidePhxPhoPlots("N3d","#gamma N3d Tracks", "1Jet/Photon/",0,5, 1);
	//TightSidePhxPhoPlots("PhiWedge","#gamma #Phi wedge", "1Jet/Photon/",0,25, 1);
	//TightSidePhxPhoPlots("TrkIso","#gamma Track Isolation [GeV]", "1Jet/Photon/",0,8, 4);
	//TightSidePhxPhoPlots("Trkpt","#gamma Track P_{T} [GeV]", "1Jet/Photon/",0,7, 1);
	//TightSidePhxPhoPlots("XCes","#gamma X-CES [cm]", "1Jet/Photon/",-30,30, 10);
	//TightSidePhxPhoPlots("ZCes","#gamma Z-CES [cm]", "1Jet/Photon/",-250,250, 10);
}
