#include<iostream>
#include"../RootTools/CommonTools.hh"
#include "TH1F.h"
#include "TH1.h"
#include<string>
#include "TF1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TList.h"
#include <sstream>


//this is to match the reweihged iso added sideband photon et 
//to total background prediction.
//I am comparing sideband + pho mc (split according to  to data in this test to
//see if the reweighing is done correctly.

void MatchReweighedSidebandToDataPhoEt()
{
	gSystem->Load("../RootTools/CommonTools_cc.so");
const float xmin=0,xpoint1=100,xpoint2=150,xpoint3=200,xpoint4=300, width1=20,width2=20,width3=50,width4=50;
float fake_pho_frac_d = 0.21; //fo pho et>40GeV

std::string sidebandpath("/Hist/SIDEBAND/1Jet/Photon/EtCorr");
std::string rewgtsidebandpath("/Hist/SIDEBAND/1Jet/Photon/EtCorr");
std::string mcpath("/Hist/CENTRAL/1Jet/Photon/EtCorr");
std::string title("#gamma^{central,E_{T}>40GeV}+=1Jet; E_{T}^{#gamma}; #epsilon . MC + (1-#epsilon) SB/E_{T}^{#gamma data}");
//std::string xtitle("E_{T}^{#gamma}");
//std::string ytitle("E_{T}^{#gamma MC}/E_{T}^{#gamma sideband};");
//std::string ytitle("E_{T}^{#gamma MC}/E_{T}^{#gamma iso addeds ideband};");

TFile *g = new TFile("PhoJets_data_gammaEt40_RewgtIsoAddedSidebandByhNominalwgts.root");
TH1* hrwgtsideband = (TH1*) g->Get(rewgtsidebandpath.c_str());
hrwgtsideband->SetTitle(title.c_str());
hrwgtsideband  = (TH1F*) MakeVariableBinHist (hrwgtsideband, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4); 
hrwgtsideband->SetLineColor(kBlack);
hrwgtsideband->SetLineWidth(2);
hrwgtsideband->Sumw2();


TH1* hsideband = (TH1*) g->Get(sidebandpath.c_str());
hsideband  = (TH1F*) MakeVariableBinHist (hsideband, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4); 
hsideband->Sumw2();


TFile *f = new TFile("PhoJets_phomc_gammaEt40_1538.root");
TH1* hmc = (TH1*) f->Get(mcpath.c_str());
hmc->SetTitle(title.c_str());
hmc  = (TH1F*) MakeVariableBinHist (hmc, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4); 
hmc->SetLineColor(kCyan);
hmc->SetLineWidth(2);



hsideband->Scale((fake_pho_frac_d)/(double) hsideband->Integral());
hmc->Scale((1-fake_pho_frac_d)/(double) hmc->Integral());
hsideband->Add(hmc);

hsideband->Scale(hrwgtsideband->Integral()/(double) hsideband->Integral());
hrwgtsideband->Divide(hsideband);

TCanvas *c = new TCanvas();
gStyle->SetOptStat(0);
gPad->SetGridx();
gPad->SetGridy();
hrwgtsideband->Draw();
hsideband->Print();
hrwgtsideband->Print();
}
