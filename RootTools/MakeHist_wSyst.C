#include <iostream>
#include <sstream>
#include <string>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>

using namespace std;


/*
FOR 30/70 VARIATION PLOTS
This version of MakeHist is used to make the plots for QCD/MC mixtures.
Manually varied the two fractions and this will plots the default and the 
changed side-by-side for comparison. 03-23-08


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

TVirtualPad* MakeHist_wSyst (int jets, const std::string& name, const std::string& title,
				 float lolimit, float hilimit, int rebin, const std::string& path, const int qcd_opt)
{
	std::string phofile("PhotonJets.root");
	std::string mcphofile("MCPhotonJets.root");
	std::string zeefile("zeeMCPhotonJets.root");
	std::string zmmfile("zmmMCPhotonJets.root");
	std::string zttfile("zttMCPhotonJets.root");
	std::string wenfile("wenMCPhotonJets.root");
	std::string wmnfile("wmnMCPhotonJets.root");
	std::string wtnfile("wtnMCPhotonJets.root");
	std::string halofile("HaloJets.root");
	std::string cosmicfile("CosmicJets.root");
	std::string qcdfile("QCDJets.root");
	
	TFile* fpho = new TFile(phofile.c_str());
	TFile* fmcpho = new TFile(mcphofile.c_str());
	TFile* fzee = new TFile(zeefile.c_str());
	TFile* fzmm = new TFile(zmmfile.c_str());
	TFile* fztt = new TFile(zttfile.c_str());
	TFile* fwen = new TFile(wenfile.c_str());
	TFile* fwmn = new TFile(wmnfile.c_str());
	TFile* fwtn = new TFile(wtnfile.c_str());
	TFile* fhalo = new TFile(halofile.c_str());
	TFile* fcosmic = new TFile(cosmicfile.c_str());
	TFile* fqcd = new TFile(qcdfile.c_str());

	if (fwen->IsZombie()) {std::cout  << "wen file not found" <<std::endl; return NULL;}


	fpho->cd();
			std::ostringstream newpath;
			newpath << "/Hist/" << path << "/";
			std::cout << newpath.str() << std::endl;
			
			if (! gDirectory->cd(newpath.str().c_str())) {
				std::cout << "path not found "<<std::endl;
				return NULL;
			}
	
			TH1F* phojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			if (! phojet){
				std::cout << "hist not found in the dir" <<std::endl;
				return NULL;
			}

			
	fhalo->cd();
			std::ostringstream newpath1;
			newpath1 << "/Hist/" << path << "/";
			gDirectory->cd(newpath1.str().c_str());
			TH1F* halojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
//	fele->cd();		
//			std::ostringstream newpath2;
//			newpath2 << "/Hist/" << path << "/";
//			gDirectory->cd(newpath2.str().c_str());
//			TH1F* elejet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			

	fzee->cd();		
			std::ostringstream newpath10;
			newpath10 << "/Hist/" << path << "/";
			gDirectory->cd(newpath10.str().c_str());
			TH1F* zeejet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fzmm->cd();		
			std::ostringstream newpath11;
			newpath11 << "/Hist/" << path << "/";
			gDirectory->cd(newpath11.str().c_str());
			TH1F* zmmjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fztt->cd();		
			std::ostringstream newpath12;
			newpath12 << "/Hist/" << path << "/";
			gDirectory->cd(newpath12.str().c_str());
			TH1F* zttjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fwen->cd();		
			std::ostringstream newpath13;
			newpath13 << "/Hist/" << path << "/";
			gDirectory->cd(newpath13.str().c_str());
			TH1F* wenjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fwmn->cd();		
			std::ostringstream newpath14;
			newpath14 << "/Hist/" << path << "/";
			gDirectory->cd(newpath14.str().c_str());
			TH1F* wmnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fwtn->cd();		
			std::ostringstream newpath15;
			newpath15 << "/Hist/" << path << "/";
			gDirectory->cd(newpath15.str().c_str());
			TH1F* wtnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			


			
	fcosmic->cd();		
			std::ostringstream newpath3;
			newpath3 << "/Hist/" << path << "/";
			gDirectory->cd(newpath3.str().c_str());
			TH1F* cosmicjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fqcd->cd();		
			std::ostringstream newpath4;
			newpath4 << "/Hist/" << path << "/";
			gDirectory->cd(newpath4.str().c_str());
			TH1F* qcdjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			TH1F* qcdjet_100 = (TH1F*) qcdjet->Clone("qcdjet_100");
			TH1F* qcdchange = (TH1F*) qcdjet->Clone("qcdchange");

	fmcpho->cd();
			std::ostringstream newpath5;
			newpath5 << "/Hist/" << path << "/";
			gDirectory->cd(newpath5.str().c_str());
			TH1F* mcphojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			TH1F* mcphochange = (TH1F*) mcphojet->Clone("mcphochange");
			
			
gStyle->SetOptStat(0);


//REBIN IF NEEDED ------------------------
phojet->Rebin(rebin);
halojet->Rebin(rebin);
//elejet->Rebin(rebin);
zeejet->Rebin(rebin);
zmmjet->Rebin(rebin);
zttjet->Rebin(rebin);
wenjet->Rebin(rebin);
wmnjet->Rebin(rebin);
wtnjet->Rebin(rebin);
cosmicjet->Rebin(rebin);
qcdjet->Rebin(rebin);
mcphojet->Rebin(rebin);
qcdjet_100->Rebin(rebin);
mcphochange->Rebin(rebin);
qcdchange->Rebin(rebin);
std::cout << "rebin  " << rebin << std::endl;


//QCD Scaling options ---------------------
//fake photon fraction = 0.319+/-0.068(syst)
float qcd_d = 0.319;
float qcd_m = 0.251;
float qcd_p = 0.387;
float mc_d = 1 - qcd_d;	//.681	//real pho + qcd fake pho = 1
float mc_p = 1 - qcd_p; //.749
float mc_m = 1 - qcd_m; //.613

std::string qcd_str, qcd_title, mc_str;

float qcd_scale = qcd_d;
float mc_scale  = mc_d;

qcd_str = "QCD (#gamma sideband, 0.319(def) of signal)";
//qcd_str = "QCD (#gamma sideband)";
mc_str = "MC (0.681(def) of signal)";
if (qcd_opt == 1) {
	qcd_scale = qcd_p;		//one goes up and other goes down
	mc_scale  = mc_m;
	qcd_str = "QCD (#gamma sideband, 0.387(def+) of signal)";
	mc_str = "#gamma MC (0.613 (def-) of signal)";
}
if (qcd_opt == -1) {
	qcd_scale = qcd_m;
	mc_scale  = mc_p;
	qcd_str = "QCD (#gamma sideband, 0.251(def-) of signal)";
	mc_str = "#gamma MC (0.749 (def+) of signal)";
}



//tune the QCD AND PHO MC ratio
float qcdchange_scale = 0.9;
float mcphochange_scale  = 0.1;
std::ostringstream qcdnum, mcnum;
qcdnum << "QCD (#gamma sideband, " << qcdchange_scale<< " of signal)";
mcnum << "#gamma MC (" << mcphochange_scale << " of signal)";
std::string qcdchange_str = qcdnum.str();
std::string mcphochange_str = mcnum.str();




// NORMALIZING ---------------------------
float haloEst 		= 0;
float cosmicEst 	= 0;

float haloNorm    = 0;
float cosmicNorm  = 0;
float qcdNorm     = 0; 
float mcphoNorm 	= 0;

float qcd100Norm     = 1;


if (jets ==1 ) {
	haloEst = 9;
	cosmicEst = 110;
}
if (jets ==2 ) {
	haloEst = 1 ;
	cosmicEst = 7 ;
}


	haloNorm    = haloEst / halojet->Integral();
	cosmicNorm  = cosmicEst / cosmicjet->Integral();
	qcdNorm     = ((phojet->Integral()) * qcd_scale)/qcdjet->Integral();  //i am looking at only 1/10 of the signal
	std::cout << "qcd scale= " << qcd_scale << "\tqcdNorm="<< qcdNorm <<std::endl;
	qcd100Norm     = ((phojet->Integral()))/qcdjet_100->Integral();  //i am looking at only 1/10 of the signal
	mcphoNorm = ((phojet->Integral()) * mc_scale)/mcphojet->Integral();

float kFac = 1.4;
float zeenorm = (2043./47605)*kFac;                   // for EWK mc see. log book#2 pp.72
float zmmnorm = (2043./49577)*kFac;                   // for EWK mc see. log book#2 pp.72
float zttnorm = (2043./28991)*kFac;                   // for EWK mc see. log book#2 pp.72
float wennorm = (2043./15408)*kFac;                   // for EWK mc see. log book#2 pp.72
float wmnnorm = (2043./7704)*kFac;                   // for EWK mc see. log book#2 pp.72
float wtnnorm = (2043./7704)*kFac;                   // for EWK mc see. log book#2 pp.72



halojet->Scale (haloNorm);
cosmicjet->Scale (cosmicNorm);
qcdjet->Scale (qcdNorm);
qcdjet_100->Scale (qcd100Norm);
//qcdjet_100ErrBand->Scale (qcd100Norm);
mcphojet->Scale (mcphoNorm);
//elejet->Scale (elenorm);
zeejet->Scale (zeenorm);
zmmjet->Scale (zmmnorm);
zttjet->Scale (zttnorm);
wenjet->Scale (wennorm);
wmnjet->Scale (wmnnorm);
wtnjet->Scale (wtnnorm);

qcdchange->Scale (phojet->Integral() * qcdchange_scale / qcdchange->Integral());
mcphochange->Scale (phojet->Integral() * mcphochange_scale / mcphochange->Integral());

//adjust X scale and move overflow to last bin
phojet    = (TH1F*) move_overflow (phojet, lolimit, hilimit, true);
halojet   = (TH1F*) move_overflow (halojet, lolimit,  hilimit,false);
//elejet    = (TH1F*) move_overflow (elejet, lolimit, hilimit, false);
zeejet    = (TH1F*) move_overflow (zeejet, lolimit, hilimit, false);
zmmjet    = (TH1F*) move_overflow (zmmjet, lolimit, hilimit, false);
zttjet    = (TH1F*) move_overflow (zttjet, lolimit, hilimit, false);
wenjet    = (TH1F*) move_overflow (wenjet, lolimit, hilimit, false);
wmnjet    = (TH1F*) move_overflow (wmnjet, lolimit, hilimit, false);
wtnjet    = (TH1F*) move_overflow (wtnjet, lolimit, hilimit, false);
cosmicjet = (TH1F*) move_overflow (cosmicjet, lolimit, hilimit, false);
qcdjet    = (TH1F*) move_overflow (qcdjet, lolimit, hilimit, false);
mcphojet  = (TH1F*) move_overflow (mcphojet, lolimit, hilimit, false);
qcdjet_100= (TH1F*) move_overflow (qcdjet_100, lolimit, hilimit, false);
mcphochange  = (TH1F*) move_overflow (mcphochange, lolimit, hilimit, false);
qcdchange= (TH1F*) move_overflow (qcdchange, lolimit, hilimit, false);
					 

mcphojet->SetLineColor(1);
mcphojet->SetFillColor(5);
qcdjet->SetLineColor(29);
qcdjet->SetFillColor(6);
cosmicjet->SetLineColor(kBlue);
cosmicjet->SetFillColor(9);
halojet->SetLineColor(kGreen);
halojet->SetFillColor(8);

zeejet->SetLineColor(kRed);
zeejet->SetFillColor(50);
zmmjet->SetLineColor(1);
zmmjet->SetFillColor(16);
zttjet->SetLineColor(1);
zttjet->SetFillColor(19);
wenjet->SetLineColor(1);
wenjet->SetFillColor(41);
wmnjet->SetLineColor(1);
wmnjet->SetFillColor(42);
wtnjet->SetLineColor(1);
wtnjet->SetFillColor(38);

mcphochange->SetLineColor(1);
mcphochange->SetFillColor(5);
qcdchange->SetLineColor(29);
qcdchange->SetFillColor(6);

phojet->SetLineColor(kBlack);
phojet->SetMarkerStyle (8);

//this is the line histo showing 100% of QCD 0% MC pho
qcdjet_100->SetLineColor(kRed);
qcdjet_100->SetFillColor(10);
TH1F* halojet_4qcd100 = (TH1F*) halojet->Clone("halojet_4qcd100");
halojet_4qcd100->SetFillColor(0);
halojet_4qcd100->SetLineColor(kRed);
halojet_4qcd100->Add(wtnjet);
halojet_4qcd100->Add(zttjet);
halojet_4qcd100->Add(wenjet);
halojet_4qcd100->Add(zeejet);
halojet_4qcd100->Add(cosmicjet);
halojet_4qcd100->Add(qcdjet_100);



std::cout << "Ewk: " << (wtnjet->Integral()+zttjet->Integral()+wenjet->Integral()+zeejet->Integral()) << std::endl;





THStack *hs = new THStack ("hs", NULL);
cosmicjet->SetXTitle (title.c_str());
//elejet->SetXTitle (title.c_str());
halojet->SetXTitle (title.c_str());
//hs->Add(elejet);
//hs->Add(zmmjet);
hs->Add(halojet);
hs->Add(wtnjet);
hs->Add(zttjet);
//hs->Add(wmnjet);
hs->Add(wenjet);
hs->Add(zeejet);
hs->Add(cosmicjet);
hs->Add(qcdjet);
hs->Add(mcphojet);

gStyle->SetLabelOffset(0.01,"XY");
gStyle->SetLabelSize(0.04,"XY");
//gStyle->SetLabelColor(kRed,"XY");
gStyle->SetLabelFont(42,"XY");

gStyle->SetTitleOffset(1.2,"X");
gStyle->SetTitleOffset(1.8,"Y");
gStyle->SetTitleSize(0.04,"X");
gStyle->SetTitleSize(0.03,"Y");
//gStyle->SetTitleColor(kBlue,"XY");
gStyle->SetTitleFont(42,"XY");

TLegend *leg = new TLegend (0.54,0.65,0.95,0.95);
std::string str_pho,str_cosmic,str_halo,str_zee,str_zmm,str_ztt,str_wen,str_wmn,str_wtn;
if (jets==1) {
	str_pho = "#gamma + >=1 Jet (1/10)";	
	str_cosmic = "Cosmic";
	str_zee = "Z->ee MC";
	str_wen = "W->e#nu MC";
	str_ztt = "Z->#tau#tau MC";
	str_wtn = "W->#tau#nu MC";
	str_wmn = "W->#mu#nu MC";
	str_zmm = "Z->#mu#mu MC";
	str_halo = "Beam Halo)";
}
if (jets==2) {
	str_pho = "#gamma + >=2 Jet (1/10)";	
	str_cosmic = "Cosmic";
	str_zee = "Z->ee MC";
	str_wen = "W->e#nu MC";
	str_ztt = "Z->#tau#tau MC";
	str_wtn = "W->#tau#nu MC";
	str_wmn = "W->#mu#nu MC";
	str_zmm = "Z->#mu#mu MC";
	str_halo = "Beam Halo";
}
	
	leg->AddEntry(phojet,str_pho.c_str());
	leg->AddEntry(mcphojet,mc_str.c_str());
	leg->AddEntry(qcdjet,qcd_str.c_str());
	leg->AddEntry(cosmicjet,str_cosmic.c_str());
	leg->AddEntry(zeejet,str_zee.c_str());
	leg->AddEntry(wenjet,str_wen.c_str());
	leg->AddEntry(zttjet,str_ztt.c_str());
	leg->AddEntry(wtnjet,str_wtn.c_str());
	//leg->AddEntry(wmnjet,str_wmn.c_str());
	//leg->AddEntry(zmmjet,str_zmm.c_str()");
	leg->AddEntry(halojet,str_halo.c_str());
	//leg->AddEntry(elejet,"e + >=1Jet");
	//leg->AddEntry(elejet,"EWK MC (e + >=1Jet)");
	//leg->AddEntry(mcphojet,"MC (normalized by Lum)");
	leg->AddEntry(halojet_4qcd100, "QCD (100% #gamma sideband)");

leg->SetBorderSize (1);



THStack *hs2 = new THStack ("hs2", NULL);
cosmicjet->SetXTitle (title.c_str());
halojet->SetXTitle (title.c_str());
hs2->Add(halojet);
hs2->Add(wtnjet);
hs2->Add(zttjet);
hs2->Add(wenjet);
hs2->Add(zeejet);
hs2->Add(cosmicjet);
hs2->Add(qcdchange);
hs2->Add(mcphochange);

TLegend *leg2 = new TLegend (0.54,0.65,0.95,0.95);
	
	leg2->AddEntry(phojet,str_pho.c_str());
	leg2->AddEntry(mcphochange,mcphochange_str.c_str());
	leg2->AddEntry(qcdchange,qcdchange_str.c_str());
	leg2->AddEntry(cosmicjet,str_cosmic.c_str());
	leg2->AddEntry(zeejet,str_zee.c_str());
	leg2->AddEntry(wenjet,str_wen.c_str());
	leg2->AddEntry(zttjet,str_ztt.c_str());
	leg2->AddEntry(wtnjet,str_wtn.c_str());
	leg2->AddEntry(halojet,str_halo.c_str());
	leg2->AddEntry(halojet_4qcd100, "QCD (100% #gamma sideband)");

leg2->SetBorderSize (1);



//TCanvas *c1 = new TCanvas("mycanvas","",1200,600);	//if i use this the obj is destroyed after each call to this function
//TCanvas *c1;
gStyle->SetCanvasDefW(1200);
gStyle->SetCanvasDefH(600);
new TCanvas();
TVirtualPad *c1 = gPad;



gPad->Divide(2,1);
gPad->cd(1);
hs->SetMinimum(0.5);
hs->SetMaximum(0.5e6);
gPad->SetBorderMode(0);
gPad->SetRightMargin(0.03);
gPad->SetTopMargin(0.03);
gPad->SetLogy();
gPad->SetTickx();
gPad->SetTicky();
hs->Draw();

hs->GetXaxis()->SetTitle (title.c_str());
std::ostringstream ytitle;
ytitle << "Bin size = " << phojet->GetBinWidth(1);
hs->GetYaxis()->SetTitle(ytitle.str().c_str());
//phojet->Scale(10);
//qcdjet_100ErrBand->Draw("SAMEE3");
halojet_4qcd100->Draw("same");
phojet->Draw("same");
leg->Draw ();


//gPad->cd();
//gPad->cd(2);
c1->cd(2);
hs2->SetMinimum(0.5);
hs2->SetMaximum(0.5e6);
gPad->SetBorderMode(0);
gPad->SetRightMargin(0.03);
gPad->SetTopMargin(0.03);
gPad->SetLogy();
gPad->SetTickx();
gPad->SetTicky();
hs2->Draw();

hs2->GetXaxis()->SetTitle (title.c_str());
hs2->GetYaxis()->SetTitle(ytitle.str().c_str());
halojet_4qcd100->Draw("same");
phojet->Draw("same");
leg2->Draw ();

//std::ostringstream str;
//str << "plot" << jets << "_" << name << ".gif";
//std::cout << str.str() << std::endl;
//c1->Print (str.str().c_str());
return c1;
};





void MakeHist_wSyst (int jets,std::string name,  std::string opt = "")
{
  const bool opt_qcdup = (opt.find ("QCDUP") != std::string::npos);
  const bool opt_qcddown = (opt.find ("QCDDOWN") != std::string::npos);


	int qcd_scale = 0;
	//TVirtualPad *pad_old = gPad;	
	TVirtualPad *c1;	

	
	if (opt_qcdup) qcd_scale = 1;
	if (opt_qcddown) qcd_scale = -1;
	
	if (jets == 1) {
	  if (name == "InvMass") 	c1 = MakeHist_wSyst (jets, "InvMass","Invariant Mass( #gamma ,Lead Jet)",0,1000,5,"1Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "pet") 	c1 = MakeHist_wSyst (jets, "EtCorr","E_{T}^{ #gamma}",0,400,5,"1Jet/Photon", qcd_scale);
	  else if (name == "met")	c1 = MakeHist_wSyst (jets, "Met","MET",0,200,1,"1Jet/Event", qcd_scale);
	  else if (name == "njet") c1 = 	MakeHist_wSyst (jets, "NJet15","NJet15",0,15,1,"1Jet/Event", qcd_scale);
	  else if (name == "ht")  	c1 = MakeHist_wSyst (jets, "Ht","H_{T}",0,800,2,"1Jet/Event", qcd_scale);
	  else if (name == "jet") 	c1 = MakeHist_wSyst (jets, "EtCorr","E_{T}^{lead Jet}",0,400,5,"1Jet/LeadJet", qcd_scale);
	  else if (name == "etratio") c1 = MakeHist_wSyst (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,2,"1Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "delphi") 	c1 = MakeHist_wSyst (jets, "DelPhi","#Delta #phi^{ #gamma,lead jet}",0,5,1,"1Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "delr") c1 = MakeHist_wSyst (jets, "DelR", "#Delta R^{ #gamma, lead Jet}",0,7,1,"1Jet/PhotonLeadJet", qcd_scale);
	}
	if (jets == 2) {
	  if (name == "InvMass") 	c1 = MakeHist_wSyst (jets, "InvMass","Invariant Mass( #gamma,Two Lead Jets)",0,1000,5,"2Jet/Photon2Jets", qcd_scale);
	  else if (name == "Met")	c1 = MakeHist_wSyst (jets, "Met","MET",0,200,1,"2Jet/Event", qcd_scale);
	  if (name == "jetsInvMass") 	c1 = MakeHist_wSyst (jets, "InvMass","Invariant Mass(Lead two Lead Jets)",0,1000,5,"2Jet/2Jets", qcd_scale);
	  else if (name == "NJet15") 	c1 = MakeHist_wSyst (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale);
	  else if (name == "PhoEt") 	c1 = MakeHist_wSyst (jets, "EtCorr","E_{T}^{ #gamma}",0,400,5,"2Jet/Photon", qcd_scale);
	  else if (name == "pj1DelPhi") 	c1 = MakeHist_wSyst (jets, "DelPhi","#Delta #phi^{ #gamma,lead jet}",0,4,1,"2Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "pj2DelPhi") 	c1 = MakeHist_wSyst (jets, "DelPhi","#Delta #phi^{ #gamma,2nd lead jet}",0,4,1,"2Jet/Photon2ndLeadJet", qcd_scale);
	  else if (name == "j1j2DelPhi") 	c1 = MakeHist_wSyst (jets, "DelPhi","#Delta #phi^{lead jet,2nd lead jet}",0,4,1,"2Jet/2Jets", qcd_scale);
	  else if (name == "Ht")  	c1 = MakeHist_wSyst (jets, "Ht","H_{T}",0,800,2,"2Jet/Event", qcd_scale);
	  else if (name == "j1Et") 	c1 = MakeHist_wSyst (jets, "EtCorr","E_{T}^{lead Jet}",0,400,5,"2Jet/LeadJet", qcd_scale);
	  else if (name == "j2Et") 	c1 = MakeHist_wSyst (jets, "EtCorr","E_{T}^{2nd lead Jet}",0,400,5,"2Jet/SecondLeadJet", qcd_scale);
	  else if (name == "pj1Etratio") c1 = MakeHist_wSyst (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,1,"2Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "pj2Etratio") c1 = MakeHist_wSyst (jets, "EtRatio","E_{T} ratio (2nd lead Jet/ #gamma)",0,10,1,"2Jet/Photon2ndLeadJet", qcd_scale);
	  else if (name == "j1j2Etratio") c1 = MakeHist_wSyst (jets, "EtRatio","E_{T} ratio (2nd lead Jet/lead Jet)",0,2.5,3,"2Jet/2Jets", qcd_scale);
	}
	

//if (gPad != pad_old)
//{
	//TCanvas *c = dynamic_cast<TCanvas*>(gPad);
	if (c1)
	{
		std::ostringstream str;
		str << "plot" << jets << "_" << name << ".gif";
		c1->Print (str.str().c_str());
	};
//};

	
};




void MakeHist_wSyst ()
{
	MakeHist_wSyst (1, "InvMass");
	MakeHist_wSyst (1, "pet");
	MakeHist_wSyst (1, "met");
	MakeHist_wSyst (1, "njet");
	MakeHist_wSyst (1, "ht");
	MakeHist_wSyst (1, "jet");
	MakeHist_wSyst (1, "etratio");
	//MakeHist_wSyst (1, "delphi");
	MakeHist_wSyst (1, "delr");
	MakeHist_wSyst (2, "InvMass");
	MakeHist_wSyst (2, "Met");
	MakeHist_wSyst (2, "jetsInvMass");
	MakeHist_wSyst (2, "NJet15");
	MakeHist_wSyst (2, "PhoEt");
	//MakeHist_wSyst (2, "pj1DelPhi");
	//MakeHist_wSyst (2, "pj2DelPhi");
	//MakeHist_wSyst (2, "j1j2DelPhi");
	MakeHist_wSyst (2, "Ht");
	MakeHist_wSyst (2, "j1Et");
	MakeHist_wSyst (2, "j2Et");
	MakeHist_wSyst (2, "pj1Etratio");
	MakeHist_wSyst (2, "pj2Etratio");
	MakeHist_wSyst (2, "j1j2Etratio");

};
