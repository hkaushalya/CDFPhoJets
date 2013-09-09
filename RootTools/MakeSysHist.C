#include <iostream>
#include <sstream>
#include <string>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TF1.h>

using namespace std;


/*
I used this module to get the systematic error plots where I used 100% QCD sideband.
I tighten up the commong cuts in loose and tight pho id cuts one at a time and make a new QCD plot
for each change. new plots are normalized to default QCD number of events and devided by it.
the max in each bin of the plot as taken as the systematic error for that bin.

*/




TH1 *move_overflow (TH1 *hist, float lolimit, float hilimit, bool use_errors)
{
	gStyle->SetCanvasColor (10);
	gStyle->SetCanvasBorderSize (0);
	gStyle->SetCanvasBorderMode (0);
				 
	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (1);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);


	TH1 *result = new TH1F ((std::string (hist->GetName())+"_tmp").c_str(),hist->GetTitle(),(int)((hilimit - lolimit)/hist->GetBinWidth(1)),lolimit, hilimit);
	result->SetDirectory (NULL);

	if (use_errors)
		result->Sumw2();
	for (unsigned bin = 1; bin <= unsigned (hist->GetNbinsX()); ++ bin)
	{
		unsigned target = bin;
		if (hist->GetXaxis()->GetBinCenter (bin) > hilimit)
			target = result->GetNbinsX();
		
		float val0 = result->GetBinContent (target);
		float val1 = hist->GetBinContent (bin);
		
		result->SetBinContent (target, val0 + val1);
		
		if (use_errors)
		{
			float err0 = result->GetBinError (target);
			float err1 = hist->GetBinError (bin);
			result->SetBinError (target, sqrt (err0*err0 + err1*err1));
		};
	};

	return result;
};

void MakeSysHist (int jets, const std::string& name, const std::string& title,
				 float lolimit, float hilimit, int rebin, const std::string& path, const int qcd_opt)
{

	//std::string qcdfile("PhotonJets.root");
	std::string qcdfile("QCDJets.root");
	std::string hademfile("Systematics_TightHadEm.root");
	std::string isofile("Systematics_TightIso.root");
	std::string trkptfile("Systematics_TightTrkPt.root");
	std::string trkisofile("Systematics_TightTrkIso.root");
	
	TFile* fqcd = new TFile(qcdfile.c_str());
	TFile* fhadem = new TFile(hademfile.c_str());
	TFile* fiso = new TFile(isofile.c_str());
	TFile* ftrkpt = new TFile(trkptfile.c_str());
	TFile* ftrkiso = new TFile(trkisofile.c_str());


	fqcd->cd();
			std::ostringstream newpath;
			newpath << "/Hist/" << path << "/";
			std::cout << newpath.str() << std::endl;
			
			if (! gDirectory->cd(newpath.str().c_str())) {
				std::cout << "path not found "<<std::endl;
				return;
			}
	
			TH1F* qcdjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			if (! qcdjet){
				std::cout << "hist not found in the dir" <<std::endl;
				return;
			}

	fhadem->cd();		
			std::ostringstream newpath1;
			newpath1 << "/Hist/" << path << "/";
			gDirectory->cd(newpath1.str().c_str());
			TH1F* hademjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fiso->cd();		
			std::ostringstream newpath2;
			newpath2 << "/Hist/" << path << "/";
			gDirectory->cd(newpath2.str().c_str());
			TH1F* isojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	ftrkpt->cd();		
			std::ostringstream newpath3;
			newpath3 << "/Hist/" << path << "/";
			gDirectory->cd(newpath3.str().c_str());
			TH1F* trkptjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	ftrkiso->cd();		
			std::ostringstream newpath4;
			newpath4 << "/Hist/" << path << "/";
			gDirectory->cd(newpath4.str().c_str());
			TH1F* trkisojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());



gStyle->SetOptStat(0);


//REBIN IF NEEDED ------------------------
qcdjet->Rebin(rebin);
hademjet->Rebin(rebin);
isojet->Rebin(rebin);
trkptjet->Rebin(rebin);
trkisojet->Rebin(rebin);

std::cout << "rebin  " << rebin << std::endl;



// NORMALIZING ---------------------------

hademjet->Scale (qcdjet->Integral() * 1.0 /hademjet->Integral());
isojet->Scale (qcdjet->Integral() * 1.0 /isojet->Integral());
trkptjet->Scale (qcdjet->Integral() * 1.0 /trkptjet->Integral());
trkisojet->Scale (qcdjet->Integral() * 1.0 /trkisojet->Integral());

//adjust X scale and move overflow to last bin
qcdjet    = (TH1F*) move_overflow (qcdjet, lolimit, hilimit, true);
hademjet   = (TH1F*) move_overflow (hademjet, lolimit,  hilimit,false);
isojet   = (TH1F*) move_overflow (isojet, lolimit,  hilimit,false);
trkptjet   = (TH1F*) move_overflow (trkptjet, lolimit,  hilimit,false);
trkisojet   = (TH1F*) move_overflow (trkisojet, lolimit,  hilimit,false);
					 
//qcdjet->SetLineColor(kBlack);
hademjet->SetLineColor(kGreen);
isojet->SetLineColor(kRed);
trkptjet->SetLineColor(kBlue);
trkisojet->SetLineColor(kBlack);
hademjet->SetMarkerStyle(22);
isojet->SetMarkerStyle(22);
trkptjet->SetMarkerStyle(22);
trkisojet->SetMarkerStyle(22);
hademjet->SetMarkerColor(kGreen);
isojet->SetMarkerColor(kRed);
trkptjet->SetMarkerColor(kBlue);
trkisojet->SetMarkerColor(kBlack);


TLegend *leg = new TLegend (0.25,0.7,0.45,0.9);
std::string str_pho,str_hadem, str_iso, str_trkpt, str_trkiso;
if (jets==1) {
	str_pho = "#gamma + >=1 Jet";	
	str_hadem = "HadEm -> Tight";	
	str_iso = "Iso -> Tight";	
	str_trkpt = "TrkPt -> Tight";	
	str_trkiso = "TrkIso -> Tight";	
}
if (jets==2) {
	str_pho = "#gamma + >=2 Jet";	
	str_hadem = "HadEm -> Tight";	
	str_iso = "Iso -> Tight";	
	str_trkpt = "TrkPt -> Tight";	
	str_trkiso = "TrkIso -> Tight";	
}
	
	leg->AddEntry(qcdjet,str_pho.c_str());
	leg->AddEntry(hademjet,str_hadem.c_str());
	leg->AddEntry(isojet,str_iso.c_str());
	leg->AddEntry(trkptjet,str_trkpt.c_str());
	leg->AddEntry(trkisojet,str_trkiso.c_str());
	
leg->SetBorderSize (1);



hademjet->Divide(qcdjet);
isojet->Divide(qcdjet);
trkptjet->Divide(qcdjet);
trkisojet->Divide(qcdjet);
qcdjet->Divide(qcdjet);


//hademjet->GetYaxis()->SetRange(-1,4);
hademjet->SetMaximum(6);
hademjet->SetTitle(qcdjet->GetTitle());
hademjet->Draw();
isojet->Draw("same");
trkptjet->Draw("same");
trkisojet->Draw("same");
leg->Draw();



};



void MakeSysHist (int jets,std::string name,  std::string opt = "")
{
  const bool opt_qcdup = (opt.find ("QCDUP") != std::string::npos);
  const bool opt_qcddown = (opt.find ("QCDDOWN") != std::string::npos);


	int qcd_scale = 0;
	TVirtualPad *pad_old = gPad;	

	
	if (opt_qcdup) qcd_scale = 1;
	if (opt_qcddown) qcd_scale = -1;
	
	if (jets == 1) {
	  if (name == "InvMass") 	MakeSysHist (jets, "InvMass","Invariant Mass(#gamma,Lead Jet)",0,1000,5,"1Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "pet") 	MakeSysHist (jets, "EtCorr","E_{T}^{#gamma}",0,400,5,"1Jet/Photon", qcd_scale);
	  else if (name == "met")	MakeSysHist (jets, "Met","MET",0,200,1,"1Jet/Event", qcd_scale);
	  else if (name == "njet") 	MakeSysHist (jets, "NJet15","NJet15",0,15,1,"1Jet/Event", qcd_scale);
	  else if (name == "ht")  	MakeSysHist (jets, "Ht","H_{T}",0,800,2,"1Jet/Event", qcd_scale);
	  else if (name == "jet") 	MakeSysHist (jets, "EtCorr","E_{T}^{lead Jet}",0,400,5,"1Jet/LeadJet", qcd_scale);
	  else if (name == "etratio") MakeSysHist (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,2,"1Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "delphi") 	MakeSysHist (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,5,1,"1Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "delr") MakeSysHist (jets, "DelR", "#Delta R^{#gamma, lead Jet}",0,7,1,"1Jet/PhotonLeadJet", qcd_scale);
	}
	if (jets == 2) {
	  if (name == "InvMass") 	MakeSysHist (jets, "InvMass","Invariant Mass(#gamma,Two Lead Jets)",0,1000,5,"2Jet/Photon2Jets", qcd_scale);
	  else if (name == "Met")	MakeSysHist (jets, "Met","MET",0,200,1,"2Jet/Event", qcd_scale);
	  if (name == "jetsInvMass") 	MakeSysHist (jets, "InvMass","Invariant Mass(Lead two Lead Jets)",0,1000,5,"2Jet/2Jets", qcd_scale);
	  else if (name == "NJet15") 	MakeSysHist (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale);
	  else if (name == "PhoEt") 	MakeSysHist (jets, "EtCorr","E_{T}^{#gamma}",0,400,5,"2Jet/Photon", qcd_scale);
	  else if (name == "pj1DelPhi") 	MakeSysHist (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,4,1,"2Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "pj2DelPhi") 	MakeSysHist (jets, "DelPhi","#Delta #phi^{#gamma,2nd lead jet}",0,4,1,"2Jet/Photon2ndLeadJet", qcd_scale);
	  else if (name == "j1j2DelPhi") 	MakeSysHist (jets, "DelPhi","#Delta #phi^{lead jet,2nd lead jet}",0,4,1,"2Jet/2Jets", qcd_scale);
	  else if (name == "Ht")  	MakeSysHist (jets, "Ht","H_{T}",0,800,2,"2Jet/Event", qcd_scale);
	  else if (name == "j1Et") 	MakeSysHist (jets, "EtCorr","E_{T}^{lead Jet}",0,400,5,"2Jet/LeadJet", qcd_scale);
	  else if (name == "j2Et") 	MakeSysHist (jets, "EtCorr","E_{T}^{2nd lead Jet}",0,400,5,"2Jet/SecondLeadJet", qcd_scale);
	  else if (name == "pj1Etratio") MakeSysHist (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,1,"2Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "pj2Etratio") MakeSysHist (jets, "EtRatio","E_{T} ratio (2nd lead Jet/ #gamma)",0,10,1,"2Jet/Photon2ndLeadJet", qcd_scale);
	  else if (name == "j1j2Etratio") MakeSysHist (jets, "EtRatio","E_{T} ratio (2nd lead Jet/lead Jet)",0,2.5,3,"2Jet/2Jets", qcd_scale);
	}
	

if (gPad != pad_old)
{
	TCanvas *c = dynamic_cast<TCanvas*>(gPad);
	std::ostringstream str;
	str << "plot" << jets << "_" << name << ".gif";
	c->Print (str.str().c_str());
};

	
};




void MakeSysHist ()
{
	MakeSysHist (1, "InvMass");
	MakeSysHist (1, "pet");
	MakeSysHist (1, "met");
	MakeSysHist (1, "njet");
	MakeSysHist (1, "ht");
	MakeSysHist (1, "jet");
	MakeSysHist (1, "etratio");
	//MakeSysHist (1, "delphi");
	MakeSysHist (1, "delr");
	MakeSysHist (2, "InvMass");
	MakeSysHist (2, "Met");
	MakeSysHist (2, "jetsInvMass");
	MakeSysHist (2, "NJet15");
	MakeSysHist (2, "PhoEt");
	//MakeSysHist (2, "pj1DelPhi");
	//MakeSysHist (2, "pj2DelPhi");
	//MakeSysHist (2, "j1j2DelPhi");
	MakeSysHist (2, "Ht");
	MakeSysHist (2, "j1Et");
	MakeSysHist (2, "j2Et");
	MakeSysHist (2, "pj1Etratio");
	MakeSysHist (2, "pj2Etratio");
	MakeSysHist (2, "j1j2Etratio");
};
