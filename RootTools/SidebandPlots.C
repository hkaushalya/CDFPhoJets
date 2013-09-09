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

void SidebandPlots(const std::string name, const std::string xtitle, const std::string path1, const std::string path2,
						const float xmin, const float xmax, const int phxreject = 1, const int rebin=1)
{
	gStyle->SetOptStat("");
	gStyle->SetCanvasColor(10);
	gStyle->SetCanvasBorderMode(0);
		
	TFile *file;
	//if (phxreject) file = new TFile("../PhoFlat/phx_a4.root");
	//else file = new TFile("../PhoFlat/phx_b4.root");
	file = new TFile("../PhoFlat/TightSidePhx5M.root");

	if (file->IsZombie())
	{
		std::cout <<  "file not found" << std::endl;
		exit (1);
	}

	file->cd(path1.c_str());
	TH1* sigPho = (TH1*) gDirectory->FindObjectAny(name.c_str());
	assert (sigPho != NULL && "sigPho is NULL");

	file->cd();
	file->cd(path2.c_str());
	TH1* sidePho = (TH1*) gDirectory->FindObjectAny(name.c_str());
	assert (sidePho != NULL && "sidePho is NULL");

	sigPho->Scale(1.0/sigPho->Integral());
	sidePho->Scale(1.0/sidePho->Integral());

	sigPho->Rebin(rebin);
	sidePho->Rebin(rebin);



	//new TCanvas();
	//gPad->SetLogy();
	//sigPho->Draw();
	new TCanvas();
	//gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	sigPho->Divide(sidePho);
	sigPho->SetMarkerStyle(8);
	sigPho->SetMarkerColor(kBlue);



	std::stringstream title, ytitle;

	title << "Both Signal and Sideband normalized to 1. (from 1.2M event sample)";
	ytitle << "#gamma^{Tight} / #gamma^{Sideband}     BinSize = " << sidePho->GetXaxis()->GetBinWidth(1);
	
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


	sigPho->Draw("E1P");
	gStyle->SetTitleBorderSize(0);
	sigPho->GetXaxis()->SetRangeUser(xmin, xmax);

	TF1* line = new TF1("line","1",-500,1200);
	line->SetLineWidth(1);
	line->SetLineColor(kRed);
	line->Draw("SAME");


	TLegend *leg = new TLegend (0.6,0.75,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);
 	leg->AddEntry(sigPho,"DATA (w. stat error)");
 	leg->AddEntry(line,"y = 1");
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
		cname << "TightSidebandRatio_5M"<< name  << ".gif";
		gPad->Print(cname.str().c_str());
		cnamepdf << "TightSidebandRatio_5M"<< name  << ".pdf";
		gPad->Print(cnamepdf.str().c_str(),"pdf");
	}

}


void SidebandPlots()
{
	SidebandPlots("Ht", "H_{T} [GeV]", "/Hist/SIGNAL/1Jet/Event/", "/Hist/SIDEBAND/1Jet/Event/",0,1000, 1, 10);
	SidebandPlots("Sumet", "#Sigma E_{T} [GeV]", "/Hist/SIGNAL/1Jet/Event/", "/Hist/SIDEBAND/1Jet/Event/",0,600, 1, 10);
	SidebandPlots("EtCorr","E_{T}^{#gamma} [GeV]", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",0,350, 1, 10);
	SidebandPlots("NJet15", "NJet", "/Hist/SIGNAL/1Jet/Event/", "/Hist/SIDEBAND/1Jet/Event/", 0,10,1, 1);
	SidebandPlots("Met", "MET [GeV]", "/Hist/SIGNAL/1Jet/Event/", "/Hist/SIDEBAND/1Jet/Event/",0,400, 1, 5);
	SidebandPlots("N12Vertices", "Class 12 Vertices", "/Hist/SIGNAL/1Jet/Event/", "/Hist/SIDEBAND/1Jet/Event/", 0,15,1, 1);
	SidebandPlots("Ces2Strip","#gamma 2^{nd} CES STRIP Energy [GeV]", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",0,10, 1, 1);
	SidebandPlots("Ces2Wire","#gamma 2^{nd} CES WIRE Energy [GeV]", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",0,20, 1, 1);
	SidebandPlots("Chi2Mean","#gamma Chi^{2} Mean", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",0,25, 1, 1);
	SidebandPlots("EmTime","#gamma EM Timing [ns]", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",-100,150, 1, 5);
	SidebandPlots("HadEm","#gamma Had/Em", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",0.0,1.0, 1, 1);
	SidebandPlots("IsoEtCorr","#gamma IsoEtCorr [GeV]", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",-5,8, 1, 1);
	SidebandPlots("N3d","#gamma N3d Tracks", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",0,5, 1, 1);
	SidebandPlots("PhiWedge","#gamma #Phi wedge", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",0,25, 1, 1);
	SidebandPlots("TrkIso","#gamma Track Isolation [GeV]", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",0,8, 1, 4);
	SidebandPlots("Trkpt","#gamma Track P_{T} [GeV]", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",0,7, 1, 1);
	SidebandPlots("XCes","#gamma X-CES [cm]", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",-30,30, 1, 10);
	SidebandPlots("ZCes","#gamma Z-CES [cm]", "/Hist/SIGNAL/1Jet/Photon/", "/Hist/SIDEBAND/1Jet/Photon/",-250,250, 1, 10);
}
