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

void MakeHist (int jets, const std::string& name, const std::string& title,
				 float lolimit, float hilimit, int rebin, const std::string& path, const int qcd_opt)
{
	std::string phofile("PhotonJets.root");
	std::string mcphofile("MCPhotonJets.root");
	//std::string elefile("EleJets.root");
	//std::string elefile("EWKEleJets.root");
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
	//TFile* fele = new TFile(elefile.c_str());
	TFile* fzee = new TFile(zeefile.c_str());
	TFile* fzmm = new TFile(zmmfile.c_str());
	TFile* fztt = new TFile(zttfile.c_str());
	TFile* fwen = new TFile(wenfile.c_str());
	TFile* fwmn = new TFile(wmnfile.c_str());
	TFile* fwtn = new TFile(wtnfile.c_str());
	TFile* fhalo = new TFile(halofile.c_str());
	TFile* fcosmic = new TFile(cosmicfile.c_str());
	TFile* fqcd = new TFile(qcdfile.c_str());

	if (fwen->IsZombie()) {std::cout  << "wen file not found" <<std::endl; return;}


	fpho->cd();
			std::ostringstream newpath;
			newpath << "/Hist/" << path << "/";
			std::cout << newpath.str() << std::endl;
			
			if (! gDirectory->cd(newpath.str().c_str())) {
				std::cout << "path not found "<<std::endl;
				return;
			}
	
			TH1F* phojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			if (! phojet){
				std::cout << "hist not found in the dir" <<std::endl;
				return;
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
			//gDirectory->cd(newpath10.str().c_str());
			gDirectory->cd(newpath1.str().c_str());
			TH1F* zeejet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fzmm->cd();		
			std::ostringstream newpath11;
			newpath11 << "/Hist/" << path << "/";
			//gDirectory->cd(newpath11.str().c_str());
			gDirectory->cd(newpath1.str().c_str());
			TH1F* zmmjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fztt->cd();		
			std::ostringstream newpath12;
			newpath12 << "/Hist/" << path << "/";
			//gDirectory->cd(newpath12.str().c_str());
			gDirectory->cd(newpath1.str().c_str());
			TH1F* zttjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fwen->cd();		
			std::ostringstream newpath13;
			newpath13 << "/Hist/" << path << "/";
			//gDirectory->cd(newpath13.str().c_str());
			gDirectory->cd(newpath1.str().c_str());
			TH1F* wenjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fwmn->cd();		
			std::ostringstream newpath14;
			newpath14 << "/Hist/" << path << "/";
			//gDirectory->cd(newpath14.str().c_str());
			gDirectory->cd(newpath1.str().c_str());
			TH1F* wmnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fwtn->cd();		
			std::ostringstream newpath15;
			newpath15 << "/Hist/" << path << "/";
			//gDirectory->cd(newpath15.str().c_str());
			gDirectory->cd(newpath1.str().c_str());
			TH1F* wtnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			


			
	fcosmic->cd();		
			std::ostringstream newpath3;
			newpath3 << "/Hist/" << path << "/";
			//gDirectory->cd(newpath3.str().c_str());
			gDirectory->cd(newpath1.str().c_str());
			TH1F* cosmicjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fqcd->cd();		
			std::ostringstream newpath4;
			newpath4 << "/Hist/" << path << "/";
			//gDirectory->cd(newpath4.str().c_str());
			gDirectory->cd(newpath1.str().c_str());
			TH1F* qcdjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			TH1F* qcdjet_100 = (TH1F*) qcdjet->Clone("qcdjet_100");
			//TH1F* qcdjet_100ErrBand = (TH1F*) qcdjet->Clone("qcdjet_100ErrBand");

	fmcpho->cd();
			std::ostringstream newpath5;
			newpath5 << "/Hist/" << path << "/";
			//gDirectory->cd(newpath5.str().c_str());
			gDirectory->cd(newpath1.str().c_str());
			TH1F* mcphojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
			
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
//qcdjet_100ErrBand->Rebin(rebin);
std::cout << "rebin  " << rebin << std::endl;
std::cout << "bin size= " << mcphojet->GetBinWidth(1) << std::endl;


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
	mc_str = "MC (0.749 (def+) of signal)";
}

qcd_scale = 0.4;
mc_scale  = 0.6;
std::ostringstream qcdnum, mcnum;
qcdnum << "QCD (#gamma sideband, " << qcd_scale<< " of signal)";
mcnum << "#gamma MC (" << mc_scale << " of signal)";
	qcd_str = qcdnum.str();
	mc_str = mcnum.str();




// NORMALIZING ---------------------------
float haloEst 		= 0;
float cosmicEst 	= 0;

float haloNorm    = 0;
float cosmicNorm  = 0;
float qcdNorm     = 0; 
float mcphoNorm 	= 0;
float elenorm 		= 0;

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
	//float qcdNorm     = (10*phojet->Integral())/qcdjet->Integral();  //i am looking at only 1/10 of the signal
	std::cout << "qcd scale= " << qcd_scale << "\tqcdNorm="<< qcdNorm <<std::endl;
	qcd100Norm     = ((phojet->Integral()))/qcdjet_100->Integral();  //i am looking at only 1/10 of the signal
	mcphoNorm = ((phojet->Integral()) * mc_scale)/mcphojet->Integral();
	//float mcphoNorm   = (2043.32)/104.3296;		//total lum of signal(p1-p13(exclude 1st 400pb-1) = 2.0433E6, LUM(mc) = 104.3296818
	//float mcphoNorm   = (2043.32 * 0.7)/104.3296;		//total lum of signal(p1-p13(exclude 1st 400pb-1) = 2.0433E6, LUM(mc) = 104.3296818

	//float elenorm = 2963/elejet->Integral();                   // data est=2963+/-343
	//float elenorm = 2043/183467;                   // for EWK mc see. log book#2 pp.72


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
//qcdjet_100ErrBand= (TH1F*) move_overflow (qcdjet_100ErrBand, lolimit, hilimit, false);
					 

//set up err band 
//for (Int_t i = 1; i < qcdjet_100->GetNbinsX(); i++) {
//	qcdjet_100ErrBand->SetBinError(i, qcdjet_100->GetBinContent(i) * 0.1);
//}

					 

mcphojet->SetLineColor(1);
mcphojet->SetFillColor(5);
qcdjet->SetLineColor(29);
qcdjet->SetFillColor(6);
cosmicjet->SetLineColor(kBlue);
cosmicjet->SetFillColor(9);
halojet->SetLineColor(kGreen);
halojet->SetFillColor(8);
//elejet->SetLineColor(kRed);
//elejet->SetFillColor(50);

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

TLegend *leg = new TLegend (0.5,0.6,0.9,0.9);
std::string str_pho,str_cosmic,str_halo,str_zee,str_zmm,str_ztt,str_wen,str_wmn,str_wtn;
if (jets==1) {
	str_pho = "#gamma + >=1 Jet (1/10)";	
	str_cosmic = "#gamma^{cosmic} + >=1 Jet";
	str_zee = "Z->ee MC (e + >=1 Jet)";
	str_wen = "W->e#nu MC (e + >=1 Jet)";
	str_ztt = "Z->#tau#tau MC (e + >=1 Jet)";
	str_wtn = "W->#tau#nu MC (e + >=1 Jet)";
	str_wmn = "W->#mu#nu MC (e + >=1 Jet)";
	str_zmm = "Z->#mu#mu MC (e + >=1 Jet)";
	str_halo = "#gamma^{halo} + >=1 Jet";
}
if (jets==2) {
	str_pho = "#gamma + >=2 Jet (1/10)";	
	str_cosmic = "#gamma^{cosmic} + >=2 Jet";
	str_zee = "Z->ee MC (e + >=2 Jet)";
	str_wen = "W->e#nu MC (e + >=2 Jet)";
	str_ztt = "Z->#tau#tau MC (e + >=2 Jet)";
	str_wtn = "W->#tau#nu MC (e + >=2 Jet)";
	str_wmn = "W->#mu#nu MC (e + >=2 Jet)";
	str_zmm = "Z->#mu#mu MC (e + >=2 Jet)";
	str_halo = "#gamma^{halo} + >=2 Jet";
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

/*
new TCanvas;
gPad->SetLogy();
phojet->Draw();

new TCanvas;
gPad->SetLogy();
halojet->Draw();
*/
/*new TCanvas;
gPad->SetLogy();
elejet->Draw();
*/
/*:
new TCanvas;
gPad->SetLogy();
cosmicjet->Draw();

new TCanvas;
gPad->SetLogy();
qcdjet->Draw();

new TCanvas;
gPad->SetLogy();
mcphojet->Draw();
*/

//TCanvas *c1= new TCanvas("mycanvas","",1000,500);
//c1->Divide(2,1);
//c1->cd(1);
new TCanvas();
hs->SetMinimum(0.5);
hs->SetMaximum(0.5e6);
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
/*
new TCanvas();
gPad->SetLogy();
hs->Draw("nostack");
phojet->Draw("same");
leg->Draw();
*/
};

void MakeHist (int jets,std::string name,  std::string opt = "")
{
  const bool opt_qcdup = (opt.find ("QCDUP") != std::string::npos);
  const bool opt_qcddown = (opt.find ("QCDDOWN") != std::string::npos);


	int qcd_scale = 0;
	TVirtualPad *pad_old = gPad;	

	
	if (opt_qcdup) qcd_scale = 1;
	if (opt_qcddown) qcd_scale = -1;
	
	if (jets == 1) {
	  if (name == "InvMass") 	MakeHist (jets, "InvMass","Invariant Mass(#gamma,Lead Jet)",0,1000,5,"1Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "pet") 	MakeHist (jets, "EtCorr","E_{T}^{#gamma}",0,400,5,"1Jet/Photon", qcd_scale);
	  else if (name == "met")	MakeHist (jets, "Met","MET",0,200,1,"1Jet/Event", qcd_scale);
	  else if (name == "njet") 	MakeHist (jets, "NJet15","NJet15",0,15,1,"1Jet/Event", qcd_scale);
	  else if (name == "ht")  	MakeHist (jets, "Ht","H_{T}",0,800,2,"1Jet/Event", qcd_scale);
	  else if (name == "jet") 	MakeHist (jets, "EtCorr","E_{T}^{lead Jet}",0,400,5,"1Jet/LeadJet", qcd_scale);
	  else if (name == "etratio") MakeHist (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,2,"1Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "delphi") 	MakeHist (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,5,1,"1Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "delr") MakeHist (jets, "DelR", "#Delta R^{#gamma, lead Jet}",0,7,1,"1Jet/PhotonLeadJet", qcd_scale);
	}
	if (jets == 2) {
	  if (name == "InvMass") 	MakeHist (jets, "InvMass","Invariant Mass(#gamma,Two Lead Jets)",0,1000,5,"2Jet/Photon2Jets", qcd_scale);
	  else if (name == "Met")	MakeHist (jets, "Met","MET",0,200,1,"2Jet/Event", qcd_scale);
	  if (name == "jetsInvMass") 	MakeHist (jets, "InvMass","Invariant Mass(Lead two Lead Jets)",0,1000,5,"2Jet/2Jets", qcd_scale);
	  else if (name == "NJet15") 	MakeHist (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale);
	  else if (name == "PhoEt") 	MakeHist (jets, "EtCorr","E_{T}^{#gamma}",0,400,5,"2Jet/Photon", qcd_scale);
	  else if (name == "pj1DelPhi") 	MakeHist (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,4,1,"2Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "pj2DelPhi") 	MakeHist (jets, "DelPhi","#Delta #phi^{#gamma,2nd lead jet}",0,4,1,"2Jet/Photon2ndLeadJet", qcd_scale);
	  else if (name == "j1j2DelPhi") 	MakeHist (jets, "DelPhi","#Delta #phi^{lead jet,2nd lead jet}",0,4,1,"2Jet/2Jets", qcd_scale);
	  else if (name == "Ht")  	MakeHist (jets, "Ht","H_{T}",0,800,2,"2Jet/Event", qcd_scale);
	  else if (name == "j1Et") 	MakeHist (jets, "EtCorr","E_{T}^{lead Jet}",0,400,5,"2Jet/LeadJet", qcd_scale);
	  else if (name == "j2Et") 	MakeHist (jets, "EtCorr","E_{T}^{2nd lead Jet}",0,400,5,"2Jet/SecondLeadJet", qcd_scale);
	  else if (name == "pj1Etratio") MakeHist (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,1,"2Jet/PhotonLeadJet", qcd_scale);
	  else if (name == "pj2Etratio") MakeHist (jets, "EtRatio","E_{T} ratio (2nd lead Jet/ #gamma)",0,10,1,"2Jet/Photon2ndLeadJet", qcd_scale);
	  else if (name == "j1j2Etratio") MakeHist (jets, "EtRatio","E_{T} ratio (2nd lead Jet/lead Jet)",0,2.5,3,"2Jet/2Jets", qcd_scale);
	}
	

if (gPad != pad_old)
{
	TCanvas *c = dynamic_cast<TCanvas*>(gPad);
	if (c)
	{
		std::ostringstream str;
		str << "plot" << jets << "_" << name << ".gif";
		c->Print (str.str().c_str());
	};
};

	
};




void MakeHist ()
{
	MakeHist (1, "InvMass");
	MakeHist (1, "pet");
	MakeHist (1, "met");
	MakeHist (1, "njet");
	MakeHist (1, "ht");
	MakeHist (1, "jet");
	MakeHist (1, "etratio");
	//MakeHist (1, "delphi");
	MakeHist (1, "delr");
	MakeHist (2, "InvMass");
	MakeHist (2, "Met");
	MakeHist (2, "jetsInvMass");
	MakeHist (2, "NJet15");
	MakeHist (2, "PhoEt");
	//MakeHist (2, "pj1DelPhi");
	//MakeHist (2, "pj2DelPhi");
	//MakeHist (2, "j1j2DelPhi");
	MakeHist (2, "Ht");
	MakeHist (2, "j1Et");
	MakeHist (2, "j2Et");
	MakeHist (2, "pj1Etratio");
	MakeHist (2, "pj2Etratio");
	MakeHist (2, "j1j2Etratio");
};
