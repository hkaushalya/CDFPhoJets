#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1.h"
#include "Rtypes.h"
#include <sstream>
#include "TLatex.h"

// this is to run on the final flat ntuples job output and to make plots
// of all photon  varaibles for signal/qcd/halo/cosmic/ewk
// like the ones I have on the walls!
// To Run: compile and load within the directory you have all the 
// required ROOT files and then execute.

/*************************** CVS LOG *******************************
 * $Id: DrawPhoJetVars.C,v 1.1 2009/12/23 23:11:45 samantha Exp $
 * $Log: DrawPhoJetVars.C,v $
 * Revision 1.1  2009/12/23 23:11:45  samantha
 * This is to create distributions of quantities (photon Et, Jet Et etc) that
 * goes in to making my final plots.
 *
 *
 *******************************************************************/

double GetMax(const std::vector<TH1*> hist)
{
	double dMax = 0;

		for (unsigned int i =0 ;i< (hist.size() - 1); i++)
		{
			int bin = hist.at(i)->GetMaximumBin();
			if (hist.at(i)->GetBinContent(bin) > dMax)
			{
				dMax = hist.at(i)->GetBinContent(bin);
			}
		}
	
	return dMax;

}


void DrawPhoJetVars(const int njet, const std::string which, const std::string path,
						const std::string title, const int rebin)
{
	std::string phofile("PhoJets_data.root");
	std::string mcphofile("PhoJets_phomc.root");
	std::string zeefile("PhoJets_zeemc.root");
	std::string zmmfile("PhoJets_zmmmc.root");
	std::string zttfile("PhoJets_zttmc.root");
	std::string wenfile("PhoJets_wenmc.root");
	std::string wmnfile("PhoJets_wmnmc.root");
	std::string wtnfile("PhoJets_wtnmc.root");


	TFile* fpho = new TFile(phofile.c_str());
	TFile* fmcpho = new TFile(mcphofile.c_str());
	TFile* fzee = new TFile(zeefile.c_str());
	TFile* fzmm = new TFile(zmmfile.c_str());
	TFile* fztt = new TFile(zttfile.c_str());
	TFile* fwen = new TFile(wenfile.c_str());
	TFile* fwmn = new TFile(wmnfile.c_str());
	TFile* fwtn = new TFile(wtnfile.c_str());


	if (fpho->IsZombie() ||fmcpho->IsZombie() ||
			fzee->IsZombie() ||fzmm->IsZombie() ||
			fztt->IsZombie() ||fwen->IsZombie() ||
			fwmn->IsZombie() ||fwtn->IsZombie() ) 
	{
		std::cout  << "a file not found. pl check" <<std::endl; return;
	}

	// ok path :1Jet/Photon/
	std::string sHaloDir     ("Hist/HALO/"     + path);
	std::string sCosmicDir   ("Hist/COSMIC/"   + path);
	std::string sQcdDir      ("Hist/SIDEBAND/" + path);
	std::string sSigDir      ("Hist/SIGNAL/"   + path);
	std::string sMcCentralDir("Hist/CENTRAL/"  + path);
	std::string sMcUpDir     ("Hist/EMJESUP/"  + path);
	std::string sMcDownDir   ("Hist/EMJESDOWN/"+ path);


	std::string name;
	if (njet == 1)
	{
		if (which == "PhotonEt") name = "EtCorr";
		else if (which == "PhotonChi2Mean") name = "Chi2Mean";
		else if (which == "PhotonDetEta") name = "DetEta";
		else if (which == "PhotonDetPhi") name = "DetPhi";
		else if (which == "PhotonEmTime") name = "EmTime";
		else if (which == "PhotonHadEm") name = "HadEm";
		else if (which == "PhotonIso") name = "IsoEtCorr";
		else if (which == "PhotonN3d") name = "N3d";
		else if (which == "PhotonPhiWedge") name = "PhiWedge";
		else if (which == "PhotonTrkPt") name = "Trkpt";
		else if (which == "PhotonTrkIso") name = "TrkIso";
		else if (which == "PhotonXCes") name = "XCes";
		else if (which == "PhotonZCes") name = "ZCes";
		else if (which == "LeadJetEt") name = "EtCorr";
		else if (which == "LeadJetEmfr") name = "Emfr";
		else if (which == "LeadJetNtwrs") name = "NTowers";
		else if (which == "LeadJetNtrks") name = "NTracks";
		else
		{
			std::cout << __FUNCTION__ << "::" << __LINE__ << "::unkown hist name " << std::endl;
			return;
		}

	} else if (njet == 2)
	{
		if (which == "PhotonEt") name = "EtCorr";
		else if (which == "PhotonChi2Mean") name = "Chi2Mean";
		else if (which == "PhotonDetEta") name = "DetEta";
		else if (which == "PhotonDetPhi") name = "DetPhi";
		else if (which == "PhotonEmTime") name = "EmTime";
		else if (which == "PhotonHadEm") name = "HadEm";
		else if (which == "PhotonIso") name = "IsoEtCorr";
		else if (which == "PhotonN3d") name = "N3d";
		else if (which == "PhotonPhiWedge") name = "PhiWedge";
		else if (which == "PhotonTrkPt") name = "Trkpt";
		else if (which == "PhotonTrkIso") name = "TrkIso";
		else if (which == "PhotonXCes") name = "XCes";
		else if (which == "PhotonZCes") name = "ZCes";
		else if (which == "LeadJetEt") name = "EtCorr";
		else if (which == "LeadJetEmfr") name = "Emfr";
		else if (which == "LeadJetNtwrs") name = "NTowers";
		else if (which == "LeadJetNtrks") name = "NTracks";
		else if (which == "SubLeadJetEt") name = "EtCorr";
		else if (which == "SubLeadJetEmfr") name = "Emfr";
		else if (which == "SubLeadJetNtwrs") name = "NTowers";
		else if (which == "SubLeadJetNtrks") name = "NTracks";
		else
		{
			std::cout << __FUNCTION__ << "::" << __LINE__ << "::unkown hist name " << std::endl;
			return;
		}

		
	} else
	{

		std::cout << __FUNCTION__ << "::" << __LINE__ << "::unkown njet " << std::endl;
		return;
	}

	fpho->cd();
	std::cout << "path="<<path<< std::endl;
	std::cout << "dir="<<sSigDir<< std::endl;

	if (! gDirectory->cd(sSigDir.c_str())) {
		std::cout << "path not found "<< std::endl;
		return;
	}


	TH1F* phojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (! phojet){
		std::cout << "hist not found in the dir" <<std::endl;
		return;
	}

	gDirectory->pwd();
	fpho->cd();
	gDirectory->cd(sHaloDir.c_str());
	gDirectory->pwd();
	TH1F* halojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	fpho->cd();
	gDirectory->cd(sCosmicDir.c_str());
	gDirectory->pwd();
	TH1F* cosmicjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	fpho->cd();
	gDirectory->cd(sQcdDir.c_str());
	gDirectory->pwd();
	TH1F* qcdjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());			//for QCD+MC combined method


	//MC HISTS: jesup= jesup && emup : jesdown = jesdown && emdown

	fmcpho->cd();
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* mcphojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fmcpho->cd();
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* mcphojetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fmcpho->cd();
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* mcphojetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());


	fzee->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* zeejet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fzee->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* zeejetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fzee->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* zeejetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	fzmm->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* zmmjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fzmm->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* zmmjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fzmm->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* zmmjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	fztt->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* zttjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fztt->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* zttjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fztt->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* zttjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	fwen->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* wenjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwen->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* wenjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwen->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* wenjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	fwmn->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* wmnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwmn->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* wmnjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwmn->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* wmnjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	fwtn->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* wtnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwtn->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* wtnjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwtn->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* wtnjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());




  phojet->Rebin(rebin);
  halojet->Rebin(rebin);
  cosmicjet->Rebin(rebin);
  qcdjet->Rebin(rebin);
  mcphojet->Rebin(rebin);

	
	gStyle->SetOptStat(0);
	gStyle->SetMarkerSize(2);


	TCanvas *c = new TCanvas();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();

  phojet->SetMarkerColor(kBlue);
  halojet->SetMarkerColor(kGreen);
  cosmicjet->SetMarkerColor(kCyan);
  qcdjet->SetMarkerColor(kMagenta);
  mcphojet->SetMarkerColor(kRed);

	phojet->SetMarkerStyle (20);
	halojet->SetMarkerStyle (5);
	cosmicjet->SetMarkerStyle (2);
	qcdjet->SetMarkerStyle (26);
	mcphojet->SetMarkerStyle (22);

	
//  zeejet->SetMarkerColor();
  //zmmjet->SetMarkerColor();
//  zttjet->SetMarkerColor();
//  wenjet->SetMarkerColor();
//  wmnjet->SetMarkerColor();
  //wtnjet->SetMarkerColor();

	
	std::stringstream ytitle;
	ytitle << "Bin Width = " << phojet->GetBinWidth(1);

	
	std::vector<TH1*> vHist;
	vHist.push_back(phojet);
	vHist.push_back(halojet);
	vHist.push_back(cosmicjet);
	vHist.push_back(qcdjet);
	vHist.push_back(mcphojet);
	double ymax = GetMax(vHist);
	
	phojet->SetMaximum(1);
	//phojet->GetYaxis()->SetMinimum(0);
	phojet->GetYaxis()->CenterTitle(kTRUE);
	phojet->GetYaxis()->SetTitle(ytitle.str().c_str());
	phojet->SetTitle(title.c_str());
	
	
	phojet->DrawNormalized("P");
	halojet->DrawNormalized("SAME P");
	cosmicjet->DrawNormalized("SAME P");
	qcdjet->DrawNormalized("SAME P");
	mcphojet->DrawNormalized("SAME P");
	//zeejet->DrawNormalized("SAME");
	//zmmjet->DrawNormalized("SAME");
	//zttjet->DrawNormalized("SAME");
	//wenjet->DrawNormalized("SAME");
	//wmnjet->DrawNormalized("SAME");
	//wtnjet->DrawNormalized("SAME");
	//
	//
	//
	TLegend *leg = new TLegend (0.7,0.72,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);

    leg->AddEntry(phojet,"#gamma DATA");
    leg->AddEntry(mcphojet,"PYTHIA #gamma MC");
    leg->AddEntry(qcdjet, "QCD (#gamma sideband)");
    leg->AddEntry(cosmicjet,"Cosmic #gamma");
    leg->AddEntry(halojet,"Beamhalo #gamma");
	 leg->Draw();

}

void DrawPhoJetVars()
{
/*
	DrawPhoJetVars(1, "PhotonEt", "1Jet/Photon", "E_{T}^{#gamma} : #gamma+#geq 1Jet (Normalized to unity)",5);
	DrawPhoJetVars(1, "PhotonChi2Mean", "1Jet/Photon", "#gamma Mean Chi2 : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "PhotonDetEta", "1Jet/Photon", "#eta^{#gamma} : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "PhotonDetPhi", "1Jet/Photon", "DetPhi^{#gamma} : #gamma+#geq 1Jet (Normalized to unity)",3);
	DrawPhoJetVars(1, "PhotonEmTime", "1Jet/Photon", "EM Timing^{#gamma} : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "PhotonHadEm", "1Jet/Photon", "Had/Em ratio {#gamma} : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "PhotonIso", "1Jet/Photon", " #gamma Isolation Energy : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "PhotonN3d", "1Jet/Photon", "Number of Tracks, N3d^{#gamma} : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "PhotonPhiWedge", "1Jet/Photon", "Azimuthal distribution of #gamma : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "PhotonTrkPt", "1Jet/Photon", "Track P_{T}^{#gamma} : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "PhotonTrkIso", "1Jet/Photon", "Track Isolation^{#gamma} : #gamma+#geq 1Jet (Normalized to unity)",2);
	DrawPhoJetVars(1, "PhotonXCes", "1Jet/Photon", "CES X^{#gamma} : #gamma+#geq 1Jet (Normalized to unity)",5);
	DrawPhoJetVars(1, "PhotonZCes", "1Jet/Photon", "CES Z^{#gamma} : #gamma+#geq 1Jet (Normalized to unity)",5);
*/
/*	
	DrawPhoJetVars(1, "LeadJetEt", "1Jet/LeadJet", "Lead Jet E_{T} : #gamma+#geq 1Jet (Normalized to unity)",5);
	DrawPhoJetVars(1, "LeadJetEmfr", "1Jet/LeadJet", "Lead Jet EM fraction : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "LeadJetNtwrs", "1Jet/LeadJet", "Lead Jet number of towers : #gamma+#geq 1Jet (Normalized to unity)",1);
	DrawPhoJetVars(1, "LeadJetNtrks", "1Jet/LeadJet", "Lead Jet number of tracks : #gamma+#geq 1Jet (Normalized to unity)",1);
*/
	DrawPhoJetVars(2, "LeadJetEt", "2Jet/LeadJet", "Lead Jet E_{T} : #gamma+#geq 2Jet (Normalized to unity)",5);
	DrawPhoJetVars(2, "LeadJetEmfr", "2Jet/LeadJet", "Lead Jet EM fraction : #gamma+#geq 2Jet (Normalized to unity)",1);
	DrawPhoJetVars(2, "LeadJetNtwrs", "2Jet/LeadJet", "Lead Jet number of towers : #gamma+#geq 2Jet (Normalized to unity)",1);
	DrawPhoJetVars(2, "LeadJetNtrks", "2Jet/LeadJet", "Lead Jet number of tracks : #gamma+#geq 2Jet (Normalized to unity)",1);
	DrawPhoJetVars(2, "SubLeadJetEt", "2Jet/SecondLeadJet", "Sub-Lead Jet E_{T} : #gamma+#geq 2Jet (Normalized to unity)",5);
	DrawPhoJetVars(2, "SubLeadJetEmfr", "2Jet/SecondLeadJet", "Sub-Lead Jet EM fraction : #gamma+#geq 2Jet (Normalized to unity)",1);
	DrawPhoJetVars(2, "SubLeadJetNtwrs", "2Jet/SecondLeadJet", "Sub-Lead Jet number of towers : #gamma+#geq 2Jet (Normalized to unity)",1);
	DrawPhoJetVars(2, "SubLeadJetNtrks", "2Jet/SecondLeadJet", "Sub-Lead Jet number of tracks : #gamma+#geq 2Jet (Normalized to unity)",1);

	
}
