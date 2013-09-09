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


//this is to match the sideband photon et spectrum to 
//photon mc. the hists were generated with
//delphi(p,j)>2.7 and veto jets with et>8GeV.
//see elog#1533
//I'll use this weights to reweight the same sideband
//events and see if the p=j blance peak shifts closer to
//the truth seen in MC.

//editing histoty
//initial version used to sideband (no iso add) to match photon MC
//2nd edit for iso added sideband to match photon MC
void MatchSidebandToMCPhoEt()
{
	gSystem->Load("../RootTools/CommonTools_cc.so");
const float xmin=0,xpoint1=150,xpoint2=250,xpoint3=300,xpoint4=650, width1=10,width2=20,width3=50,width4=250;

std::string datapath("/Hist/SIDEBAND/1Jet/Photon/EtCorr");
std::string mcpath("/Hist/CENTRAL/1Jet/Photon/EtCorr");
//std::string title("#gamma^{central,E_{T}>40GeV}+=1Jet:#Delta#Phi(#gamma,jet)>2.7 && veto extra jets(lev6) P_{T}>8GeV; E_{T}^{#gamma}; E_{T}^{#gamma MC}/E_{T}^{#gamma sideband};");
std::string title("#gamma^{central,E_{T}>40GeV}+=1Jet:#Delta#Phi(#gamma,jet)>2.7 && veto extra jets(lev6) P_{T}>8GeV; E_{T}^{#gamma}; E_{T}^{#gamma MC}/E_{T}^{#gamma ISO added sideband};");
std::string xtitle("E_{T}^{#gamma}");
//std::string ytitle("E_{T}^{#gamma MC}/E_{T}^{#gamma sideband};");
std::string ytitle("E_{T}^{#gamma MC}/E_{T}^{#gamma iso addeds ideband};");

TFile *g = new TFile("PhoJets_phomc_noiso_pjbalance_1533.root");
TH1* hnoiso_mc = (TH1*) g->Get(mcpath.c_str());
hnoiso_mc->SetTitle(title.c_str());
hnoiso_mc->GetXaxis()->SetTitle(xtitle.c_str());
hnoiso_mc->GetYaxis()->SetTitle(ytitle.c_str());

hnoiso_mc  = (TH1F*) MakeVariableBinHist (hnoiso_mc, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4); 
hnoiso_mc->SetLineColor(kBlack);
hnoiso_mc->SetLineWidth(2);


//TFile *f = new TFile("PhoJets_phodata_noiso_pjbalance_1533.root");
TFile *f = new TFile("PhoJets_phodata_sidebandisoadded_1536.root");
TH1* hnoiso_data = (TH1*) f->Get(datapath.c_str());
hnoiso_data  = (TH1F*) MakeVariableBinHist (hnoiso_data, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4); 
hnoiso_data->Sumw2();
hnoiso_data->Scale(hnoiso_mc->Integral()/(double) hnoiso_data->Integral());
hnoiso_mc->Divide(hnoiso_data);

TCanvas *c = new TCanvas();
hnoiso_mc->Draw();
gStyle->SetOptStat(0);
gStyle->SetOptFit(1);
gPad->SetGridx();
gPad->SetGridy();

int minbin=0,maxbin=0;
FindNonZeroXrange(hnoiso_mc,minbin,maxbin);
hnoiso_mc->GetXaxis()->SetRangeUser(0,hnoiso_mc->GetBinCenter(maxbin+2));

TF1 *func = new TF1("tf1_matchSidebandToPhoMC","pol2",hnoiso_mc->GetBinCenter(minbin),hnoiso_mc->GetXaxis()->GetBinCenter(maxbin));
func->SetLineColor(8);
hnoiso_mc->Fit(func,"R+");
TPaveText *tp = new TPaveText(0.2,0.75,0.3,0.8,"NDC");
tp->AddText("Pol2 fit");
tp->Draw("same");

	TList *funcList = hnoiso_mc->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1* fitfunc = (TF1*) it->Next();
	
	//TFile file("MatchSidebandToPhoMC_withDelPhiJetVeto_1533.root","UPDATE");
	TFile file("MatchIsoAddedSidebandToPhoMC_withDelPhiJetVeto_1533.root","UPDATE");
	if ( file.Get(func->GetName()) != NULL)
	{
		std::stringstream old_tf;
		old_tf << func->GetName()<< ";1";
		file.Delete(old_tf.str().c_str());
	}
	fitfunc->Write();
	c->Write();
	hnoiso_mc->Write();
	file.ls();

	file.Close();

}
