/**********************************************************
 * This is a strip down version of JetFilterModuleV2. I
 * removed all the hists and stuff to run faster in Flat ntuple
 * making.
 *********************************************************/
/*{{{*/
/*
 * $Id: JetFilterModuleV3.cc,v 1.1 2011/05/25 21:39:12 samantha Exp $
 * $Log: JetFilterModuleV3.cc,v $
 * Revision 1.1  2011/05/25 21:39:12  samantha
 * This is a strip down version of JetFilterModuleV2. I removed all the
 * hists and stuff to run faster in Flat ntuple making.
 *
 */
/*}}}*/

#include <iomanip>
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Stntuple/loop/TStnAna.hh"
#include "samantha/Pho/JetFilterModuleV3.hh"
#include "samantha/utils/FreeFunctions.hh"
#include <set>
#include <cmath>
#include "TRandom3.h"
#include "../../RootTools/IOColors.hh"
#include <assert.h>

ClassImp(JetFilterModuleV3)
//_____________________________________________________________________________
JetFilterModuleV3::JetFilterModuleV3(const char* name, const char* title):
		TStnModule(name,title),
		bRunPermit(true),
		debug(false),
	   bNoSummary(false),
		iPrintLvl(100),
		fMinisculeEta(-1e11),
		bRem_TightPho(0), bRem_LoosePho(0), bRem_SidebandPho(0),
		bRem_TightPhoLikeEle(0), bRem_LoosePhoLikeEle(0), bRem_StdLooseEle(0)
{

	fMyMetScenario=3; // by default met is corrected for jets with Et>15
	//______________________________________ setting default cut values
	MinEtThr=15.0;
	MaxEtThr=1000.0;
	MinNjetThr=0;
	MaxNjetThr=100;
	MinJetEta=0.0;
	MaxJetEta=3.0;
	MinEt1st=0.0;
	MaxEt1st=1000.0;
	MinEt2nd=0.0;
	MaxEt2nd=1000.0;
	MinDeltaRgjMin=0.0;
	MaxDeltaRgjMin=1.0E6;
	MinDeltaEtagjMin=0.0;
	MaxDeltaEtagjMin=1.0E6;
	MinDeltaPhigjMin=0.0;
	MaxDeltaPhigjMin=1.0E6;
	MinDeltaRgj=0.0;
	MaxDeltaRgj=1.0E6;
	MinDeltaEtagj=0.0;
	MaxDeltaEtagj=1.0E6;
	MinDeltaPhigj=0.0;
	MaxDeltaPhigj=1.0E6;
	//   MinDeltaPhiJMet=0.056*2.0*asin(1.0); // ~10 degrees
	//   MaxDeltaPhiJMet=0.944*2.0*asin(1.0); // ~170 degrees
	MinDeltaPhiJMet=-10.0;
	MaxDeltaPhiJMet=0.3;
	MinMet2Et=0.0;
	MaxMet2Et=1.0E6;
	MinMet2Etlev6=0.0;
	MaxMet2Etlev6=1.0E6;
	MindZ=0.0;
	MaxdZ=1.0E6;
	MinMjj=0.0;
	MaxMjj=1.0E6;
	MinMj1g1=0.0;
	MaxMj1g1=1.0E6;
	MinMj1g2=0.0;
	MaxMj1g2=1.0E6;
	MinMj2g1=0.0;
	MaxMj2g1=1.0E6;
	MinMj2g2=0.0;
	MaxMj2g2=1.0E6;
	MinMjg=0.0;
	MaxMjg=1.0E6;
	MinQtJet=0.0;
	MaxQtJet=1.0E6;
	MinHt=0.0;
	MaxHt=1.0E6;
	MinNjet=0;
	MaxNjet=100;
	MinNjet5=0;
	MaxNjet5=100;
	MinNjet10=0;
	MaxNjet10=100;
	MinNjet15=0;
	MaxNjet15=100;
	MinNjet20=0;
	MaxNjet20=100;
	MinNjet25=0;
	MaxNjet25=100;
	MinNjet30=0;
	MaxNjet30=100;
	MinNjet35=0;
	MaxNjet35=100;
	//________________________
	fJTC_coneSize=1;  // JetClue R=0.7
	fJTC_version=5;   // version of correction, suggested by Anwar for 5.3.1, may change later
	fJTC_level=7;     // full correction
	fJTC_systcode=0;  // default correction
	//fJTC_imode=1;     // DATA is default  //this will be automatically set at run time - sam 12-09-2008
	fJTC_imode=-1;     // DATA is default  //this will be automatically set at run time - sam 12-09-2008
	std::cout<<"Hello I am JetFilterModuleV3"<<std::endl;
}

//_____________________________________________________________________________
JetFilterModuleV3::~JetFilterModuleV3() {
}

//_____________________________________________________________________________
void JetFilterModuleV3::BookHistograms() {
	char folder_name[200];
	TFolder* fol;
	TFolder* hist_folder;

	//-----------------------------------------------------------------------------
	//  clear the histogram list
	//-----------------------------------------------------------------------------
	DeleteHistograms();
	hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

	//------ booking generic histograms


	sprintf(folder_name,"JetClu-0.4"); //----- Cone 0.4
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	BookJetHistograms(fHistJet04,Form("Hist/%s",folder_name),folder_name);

	//------ booking matching histograms
	sprintf(folder_name,"MatchPhoJet0.4"); //----- Cone 0.4
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	BookMatchStudyHistograms(fMatchPhoJet04,Form("Hist/%s",folder_name),folder_name);

	return;
}

//______________________ Booking final analysis histograms

//______________________ Booking general jet histograms
void JetFilterModuleV3::BookJetHistograms(JetGeneral_t& Hist, const char* Folder, const char* algoname) {
	char name [200];
	char title[200];
	// book histograms

	sprintf(name,"toyMET_all");
	sprintf(title,"#slash{E}_{T}^{toy}, all events");
	Hist.fEvnt_toyMET_all=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fEvnt_toyMET_all,Folder);
	sprintf(name,"toyMET_cut");
	sprintf(title,"#slash{E}_{T}^{toy}, events after Cleanup cuts");
	Hist.fEvnt_toyMET_cut=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fEvnt_toyMET_cut,Folder);
	sprintf(name,"corMET_all");
	sprintf(title,"#slash{E}_{T}^{cor}, all events before Cleanup cuts");
	Hist.fEvnt_corMET_all=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fEvnt_corMET_all,Folder);

	//__________________________________ before cuts

	sprintf(title,"%s%s",algoname,": Number of Jets");
	Hist.fEvntNjet_b= new TH1F("Njet_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>5 GeV");
	Hist.fEvntNjet5_b= new TH1F("Njet5_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet5_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>10 GeV");
	Hist.fEvntNjet10_b= new TH1F("Njet10_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet10_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>15 GeV");
	Hist.fEvntNjet15_b= new TH1F("Njet15_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet15_b,Folder);
	Hist.fEvntNjet20_b= new TH1F("Njet20_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet20_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>25 GeV");
	Hist.fEvntNjet25_b= new TH1F("Njet25_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet25_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>30 GeV");
	Hist.fEvntNjet30_b= new TH1F("Njet30_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet30_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>35 GeV");
	Hist.fEvntNjet35_b= new TH1F("Njet35_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet35_b,Folder);

	sprintf(title," %s%s",algoname,": dZ, jet & best vertex");
	Hist.fEvntdZ_b= new TH1F("dZ_b",title,400,-200.0,200.0);
	AddHistogram(Hist.fEvntdZ_b,Folder);
	sprintf(title," %s%s",algoname,": raw Et, 1st jet");
	Hist.fEvntEt0_b[0]= new TH1F("Et0_b[0]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_b[0],Folder);
	sprintf(title," %s%s",algoname,": raw Et, 2nd jet");
	Hist.fEvntEt0_b[1]= new TH1F("Et0_b[1]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_b[1],Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), 1st jet");
	Hist.fEvntEt_b[0]= new TH1F("Et_b[0]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_b[0],Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), 2nd jet");
	Hist.fEvntEt_b[1]= new TH1F("Et_b[1]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_b[1],Folder);
	sprintf(title," %s%s",algoname,": detector Eta, 1st jet");
	Hist.fEvntEtaDet_b[0]= new TH1F("EtaDet_b[0]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDet_b[0],Folder);
	sprintf(title," %s%s",algoname,": detector Eta, 2nd jet");
	Hist.fEvntEtaDet_b[1]= new TH1F("EtaDet_b[1]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDet_b[1],Folder);
	sprintf(title," %s%s",algoname,": Eta, 1st jet");
	Hist.fEvntEta_b[0]= new TH1F("Eta_b[0]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEta_b[0],Folder);
	sprintf(title," %s%s",algoname,": Eta, 2nd jet");
	Hist.fEvntEta_b[1]= new TH1F("Eta_b[1]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEta_b[1],Folder);
	sprintf(title," %s%s",algoname,": Phi, 1st jet");
	Hist.fEvntPhi_b[0]= new TH1F("Phi_b[0]",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhi_b[0],Folder);
	sprintf(title," %s%s",algoname,": Phi, 2nd jet");
	Hist.fEvntPhi_b[1]= new TH1F("Phi_b[1]",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhi_b[1],Folder);
	sprintf(title," %s%s",algoname,": Theta, 1st jet");
	Hist.fEvntTheta_b[0]= new TH1F("Theta_b[0]",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntTheta_b[0],Folder);
	sprintf(title," %s%s",algoname,": Theta, 2nd jet");
	Hist.fEvntTheta_b[1]= new TH1F("Theta_b[1]",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntTheta_b[1],Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, 1st jet");
	Hist.fEvntEmFr_b[0]= new TH1F("EmFr_b[0]",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFr_b[0],Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, 2nd jet");
	Hist.fEvntEmFr_b[1]= new TH1F("EmFr_b[1]",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFr_b[1],Folder);
	sprintf(title," %s%s",algoname,": Number of towers, 1st jet");
	Hist.fEvntNTowers_b[0]= new TH1F("NTowers_b[0]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowers_b[0],Folder);
	sprintf(title," %s%s",algoname,": Number of towers, 2nd jet");
	Hist.fEvntNTowers_b[1]= new TH1F("NTowers_b[1]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowers_b[1],Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, 1st jet");
	Hist.fEvntNTracks_b[0]= new TH1F("NTracks_b[0]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracks_b[0],Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, 2nd jet");
	Hist.fEvntNTracks_b[1]= new TH1F("NTracks_b[1]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracks_b[1],Folder);
	sprintf(title," %s%s",algoname,": raw Et, extra jets");
	Hist.fEvntEt0X_b= new TH1F("Et0X_b",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0X_b,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), extra jets");
	Hist.fEvntEtX_b= new TH1F("EtX_b",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEtX_b,Folder);
	sprintf(title," %s%s",algoname,": detector Eta, extra jets");
	Hist.fEvntEtaDetX_b= new TH1F("EtaDetX_b",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDetX_b,Folder);
	sprintf(title," %s%s",algoname,": Eta, extra jets");
	Hist.fEvntEtaX_b= new TH1F("EtaX_b",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaX_b,Folder);
	sprintf(title," %s%s",algoname,": Phi, extra jets");
	Hist.fEvntPhiX_b= new TH1F("PhiX_b",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhiX_b,Folder);
	sprintf(title," %s%s",algoname,": Theta, extra jets");
	Hist.fEvntThetaX_b= new TH1F("ThetaX_b",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntThetaX_b,Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, extra jets");
	Hist.fEvntEmFrX_b= new TH1F("EmFrX_b",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFrX_b,Folder);
	sprintf(title," %s%s",algoname,": Number of towers, extra jets");
	Hist.fEvntNTowersX_b= new TH1F("NTowersX_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowersX_b,Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, extra jets");
	Hist.fEvntNTracksX_b= new TH1F("NTracksX_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracksX_b,Folder);
	sprintf(title," %s%s",algoname,": My Invariant Mass, 2 leading jets");
	Hist.fEvntMjj_b= new TH1F("Mjj_b",title,2000,0.0,2000.0);
	AddHistogram(Hist.fEvntMjj_b,Folder);
	sprintf(title," %s%s",algoname,": Theta* of 2 leading jets");
	Hist.fEvntThetaStar_b= new TH1F("ThetaStar_b",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntThetaStar_b,Folder);
	sprintf(title," %s%s",algoname,": Delta Phi of 2 leading jets");
	Hist.fEvntDeltaPhi_b= new TH1F("DeltaPhi_b",title,350,0.0,3.5);
	AddHistogram(Hist.fEvntDeltaPhi_b,Folder);
	sprintf(title," %s%s",algoname,": Delta Eta of 2 leading jets");
	Hist.fEvntDeltaEta_b= new TH1F("DeltaEta_b",title,800,-8.0,8.0);
	AddHistogram(Hist.fEvntDeltaEta_b,Folder);
	sprintf(title," %s%s",algoname,": Delta R of 2 leading jets");
	Hist.fEvntDeltaR_b= new TH1F("DeltaR_b",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR_b,Folder);
	sprintf(title," %s%s",algoname,": dR, 1st jet & any extra jet");
	Hist.fEvntDeltaR1x_b= new TH1F("DeltaR1x_b",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR1x_b,Folder);
	sprintf(title," %s%s",algoname,": dR, 2nd jet & any extra jet");
	Hist.fEvntDeltaR2x_b= new TH1F("DeltaR2x_b",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR2x_b,Folder);
	sprintf(title," %s%s",algoname,": Kt of two best jets");
	Hist.fEvntKt2jet_b= new TH1F("Kt2jet_b",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntKt2jet_b,Folder);
	sprintf(title," %s%s",algoname,": Kt of all jets");
	Hist.fEvntKtAll_b= new TH1F("KtAll_b",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntKtAll_b,Folder);

	sprintf(title," %s%s",algoname,": Njet vs. Nvx(class>=12)");
	Hist.fEvntNJet_Nvx_b= new TProfile("NJet_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet5 vs. Nvx(class>=12)");
	Hist.fEvntNJet5_Nvx_b= new TProfile("NJet5_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet5_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet10 vs. Nvx(class>=12)");
	Hist.fEvntNJet10_Nvx_b= new TProfile("NJet10_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet10_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet15 vs. Nvx(class>=12)");
	Hist.fEvntNJet15_Nvx_b= new TProfile("NJet15_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet15_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet20 vs. Nvx(class>=12)");
	Hist.fEvntNJet20_Nvx_b= new TProfile("NJet20_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet20_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet25 vs. Nvx(class>=12)");
	Hist.fEvntNJet25_Nvx_b= new TProfile("NJet25_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet25_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet30 vs. Nvx(class>=12)");
	Hist.fEvntNJet30_Nvx_b= new TProfile("NJet30_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet30_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet35 vs. Nvx(class>=12)");
	Hist.fEvntNJet35_Nvx_b= new TProfile("NJet35_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet35_Nvx_b,Folder);

	sprintf(title," %s%s",algoname,": 2 leading jets, Kt vs. Mjj");
	Hist.fEvntKt_Mjj_b= new TProfile("Kt_Mjj_b",title,400,0.0,1000.0,0.0,1000.0);
	AddHistogram(Hist.fEvntKt_Mjj_b,Folder);
	sprintf(title," %s%s",algoname,": KtAll vs. Mjj");
	Hist.fEvntKtAll_Mjj_b= new TProfile("KtAll_Mjj_b",title,400,0.0,1000.0,0.0,1000.0);
	AddHistogram(Hist.fEvntKtAll_Mjj_b,Folder);
	sprintf(title," %s%s",algoname,": raw Et vs. Njet");
	Hist.fEvntEt0_Njet_b= new TProfile("Et0_Njet_b",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_Njet_b,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6) vs. Njet");
	Hist.fEvntEt_Njet_b= new TProfile("Et_Njet_b",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_Njet_b,Folder);
	sprintf(title," %s%s",algoname,": raw Et vs. Nvx(class>=12)");
	Hist.fEvntEt0_Nvx12_b= new TProfile("Et0_Nvx12_b",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_Nvx12_b,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6) vs. Nvx(class>=12)");
	Hist.fEvntEt_Nvx12_b= new TProfile("Et_Nvx12_b",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_Nvx12_b,Folder);
	sprintf(title," %s%s",algoname,": Njet vs. Lum for events with Nvx12=1");
	Hist.fEvntNjet_Lum_b= new TProfile("Njet_Lum_b",title,50,0.0,100.0,0.0,1000.0);
	AddHistogram(Hist.fEvntNjet_Lum_b,Folder);


	//________________________________________________________________ after cuts
	sprintf(title,"%s%s",algoname,": Number of Jets");
	Hist.fEvntNjet_a= new TH1F("Njet_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>5 GeV");
	Hist.fEvntNjet5_a= new TH1F("Njet5_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet5_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>10 GeV");
	Hist.fEvntNjet10_a= new TH1F("Njet10_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet10_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>15 GeV");
	Hist.fEvntNjet15_a= new TH1F("Njet15_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet15_a,Folder);
	Hist.fEvntNjet20_a= new TH1F("Njet20_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet20_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>25 GeV");
	Hist.fEvntNjet25_a= new TH1F("Njet25_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet25_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>30 GeV");
	Hist.fEvntNjet30_a= new TH1F("Njet30_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet30_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>35 GeV");
	Hist.fEvntNjet35_a= new TH1F("Njet35_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet35_a,Folder);

	sprintf(title," %s%s",algoname,": dZ, jet & best vertex");
	Hist.fEvntdZ_a= new TH1F("dZ_a",title,400,-200.0,200.0);
	AddHistogram(Hist.fEvntdZ_a,Folder);
	sprintf(title," %s%s",algoname,": raw Et, 1st jet");
	Hist.fEvntEt0_a[0]= new TH1F("Et0_a[0]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_a[0],Folder);
	sprintf(title," %s%s",algoname,": raw Et, 2nd jet");
	Hist.fEvntEt0_a[1]= new TH1F("Et0_a[1]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_a[1],Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), 1st jet");
	Hist.fEvntEt_a[0]= new TH1F("Et_a[0]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_a[0],Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), 2nd jet");
	Hist.fEvntEt_a[1]= new TH1F("Et_a[1]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_a[1],Folder);
	sprintf(title," %s%s",algoname,": detector Eta, 1st jet");
	Hist.fEvntEtaDet_a[0]= new TH1F("EtaDet_a[0]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDet_a[0],Folder);
	sprintf(title," %s%s",algoname,": detector Eta, 2nd jet");
	Hist.fEvntEtaDet_a[1]= new TH1F("EtaDet_a[1]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDet_a[1],Folder);
	sprintf(title," %s%s",algoname,": Eta, 1st jet");
	Hist.fEvntEta_a[0]= new TH1F("Eta_a[0]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEta_a[0],Folder);
	sprintf(title," %s%s",algoname,": Eta, 2nd jet");
	Hist.fEvntEta_a[1]= new TH1F("Eta_a[1]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEta_a[1],Folder);
	sprintf(title," %s%s",algoname,": Phi, 1st jet");
	Hist.fEvntPhi_a[0]= new TH1F("Phi_a[0]",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhi_a[0],Folder);
	sprintf(title," %s%s",algoname,": Phi, 2nd jet");
	Hist.fEvntPhi_a[1]= new TH1F("Phi_a[1]",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhi_a[1],Folder);
	sprintf(title," %s%s",algoname,": Theta, 1st jet");
	Hist.fEvntTheta_a[0]= new TH1F("Theta_a[0]",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntTheta_a[0],Folder);
	sprintf(title," %s%s",algoname,": Theta, 2nd jet");
	Hist.fEvntTheta_a[1]= new TH1F("Theta_a[1]",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntTheta_a[1],Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, 1st jet");
	Hist.fEvntEmFr_a[0]= new TH1F("EmFr_a[0]",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFr_a[0],Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, 2nd jet");
	Hist.fEvntEmFr_a[1]= new TH1F("EmFr_a[1]",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFr_a[1],Folder);
	sprintf(title," %s%s",algoname,": Number of towers, 1st jet");
	Hist.fEvntNTowers_a[0]= new TH1F("NTowers_a[0]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowers_a[0],Folder);
	sprintf(title," %s%s",algoname,": Number of towers, 2nd jet");
	Hist.fEvntNTowers_a[1]= new TH1F("NTowers_a[1]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowers_a[1],Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, 1st jet");
	Hist.fEvntNTracks_a[0]= new TH1F("NTracks_a[0]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracks_a[0],Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, 2nd jet");
	Hist.fEvntNTracks_a[1]= new TH1F("NTracks_a[1]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracks_a[1],Folder);
	sprintf(title," %s%s",algoname,": raw Et, extra jets");
	Hist.fEvntEt0X_a= new TH1F("Et0X_a",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0X_a,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), extra jets");
	Hist.fEvntEtX_a= new TH1F("EtX_a",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEtX_a,Folder);
	sprintf(title," %s%s",algoname,": detector Eta, extra jets");
	Hist.fEvntEtaDetX_a= new TH1F("EtaDetX_a",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDetX_a,Folder);
	sprintf(title," %s%s",algoname,": Eta, extra jets");
	Hist.fEvntEtaX_a= new TH1F("EtaX_a",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaX_a,Folder);
	sprintf(title," %s%s",algoname,": Phi, extra jets");
	Hist.fEvntPhiX_a= new TH1F("PhiX_a",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhiX_a,Folder);
	sprintf(title," %s%s",algoname,": Theta, extra jets");
	Hist.fEvntThetaX_a= new TH1F("ThetaX_a",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntThetaX_a,Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, extra jets");
	Hist.fEvntEmFrX_a= new TH1F("EmFrX_a",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFrX_a,Folder);
	sprintf(title," %s%s",algoname,": Number of towers, extra jets");
	Hist.fEvntNTowersX_a= new TH1F("NTowersX_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowersX_a,Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, extra jets");
	Hist.fEvntNTracksX_a= new TH1F("NTracksX_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracksX_a,Folder);
	sprintf(title," %s%s",algoname,": My Invariant Mass, 2 leading jets");
	Hist.fEvntMjj_a= new TH1F("Mjj_a",title,2000,0.0,2000.0);
	AddHistogram(Hist.fEvntMjj_a,Folder);
	sprintf(title," %s%s",algoname,": Theta* of 2 leading jets");
	Hist.fEvntThetaStar_a= new TH1F("ThetaStar_a",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntThetaStar_a,Folder);
	sprintf(title," %s%s",algoname,": Delta Phi of 2 leading jets");
	Hist.fEvntDeltaPhi_a= new TH1F("DeltaPhi_a",title,350,0.0,3.5);
	AddHistogram(Hist.fEvntDeltaPhi_a,Folder);
	sprintf(title," %s%s",algoname,": Delta Eta of 2 leading jets");
	Hist.fEvntDeltaEta_a= new TH1F("DeltaEta_a",title,800,-8.0,8.0);
	AddHistogram(Hist.fEvntDeltaEta_a,Folder);
	sprintf(title," %s%s",algoname,": Delta R of 2 leading jets");
	Hist.fEvntDeltaR_a= new TH1F("DeltaR_a",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR_a,Folder);
	sprintf(title," %s%s",algoname,": dR, 1st jet & any extra jet");
	Hist.fEvntDeltaR1x_a= new TH1F("DeltaR1x_a",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR1x_a,Folder);
	sprintf(title," %s%s",algoname,": dR, 2nd jet & any extra jet");
	Hist.fEvntDeltaR2x_a= new TH1F("DeltaR2x_a",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR2x_a,Folder);
	sprintf(title," %s%s",algoname,": Kt of two best jets");
	Hist.fEvntKt2jet_a= new TH1F("Kt2jet_a",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntKt2jet_a,Folder);
	sprintf(title," %s%s",algoname,": Kt of all jets");
	Hist.fEvntKtAll_a= new TH1F("KtAll_a",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntKtAll_a,Folder);

	sprintf(title," %s%s",algoname,": Njet vs. Nvx(class>=12)");
	Hist.fEvntNJet_Nvx_a= new TProfile("NJet_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet5 vs. Nvx(class>=12)");
	Hist.fEvntNJet5_Nvx_a= new TProfile("NJet5_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet5_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet10 vs. Nvx(class>=12)");
	Hist.fEvntNJet10_Nvx_a= new TProfile("NJet10_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet10_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet15 vs. Nvx(class>=12)");
	Hist.fEvntNJet15_Nvx_a= new TProfile("NJet15_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet15_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet20 vs. Nvx(class>=12)");
	Hist.fEvntNJet20_Nvx_a= new TProfile("NJet20_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet20_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet25 vs. Nvx(class>=12)");
	Hist.fEvntNJet25_Nvx_a= new TProfile("NJet25_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet25_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet30 vs. Nvx(class>=12)");
	Hist.fEvntNJet30_Nvx_a= new TProfile("NJet30_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet30_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet35 vs. Nvx(class>=12)");
	Hist.fEvntNJet35_Nvx_a= new TProfile("NJet35_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet35_Nvx_a,Folder);

	sprintf(title," %s%s",algoname,": 2 leading jets, Kt vs. Mjj");
	Hist.fEvntKt_Mjj_a= new TProfile("Kt_Mjj_a",title,400,0.0,1000.0,0.0,1000.0);
	AddHistogram(Hist.fEvntKt_Mjj_a,Folder);
	sprintf(title," %s%s",algoname,": KtAll vs. Mjj");
	Hist.fEvntKtAll_Mjj_a= new TProfile("KtAll_Mjj_a",title,400,0.0,1000.0,0.0,1000.0);
	AddHistogram(Hist.fEvntKtAll_Mjj_a,Folder);
	sprintf(title," %s%s",algoname,": raw Et vs. Njet");
	Hist.fEvntEt0_Njet_a= new TProfile("Et0_Njet_a",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_Njet_a,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6) vs. Njet");
	Hist.fEvntEt_Njet_a= new TProfile("Et_Njet_a",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_Njet_a,Folder);
	sprintf(title," %s%s",algoname,": raw Et vs. Nvx(class>=12)");
	Hist.fEvntEt0_Nvx12_a= new TProfile("Et0_Nvx12_a",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_Nvx12_a,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6) vs. Nvx(class>=12)");
	Hist.fEvntEt_Nvx12_a= new TProfile("Et_Nvx12_a",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_Nvx12_a,Folder);
	sprintf(title," %s%s",algoname,": Njet vs. Lum for events with Nvx12=1");
	Hist.fEvntNjet_Lum_a= new TProfile("Njet_Lum_a",title,50,0.0,100.0,0.0,1000.0);
	AddHistogram(Hist.fEvntNjet_Lum_a,Folder);

	return;
}


//_________________________________________ filling general histo for particular Jet Cone
void JetFilterModuleV3::FillJetHistogramsB(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12) {


	//   double toyMet=MetToy_min+(MetToy_max-MetToy_min)*(gRandom->Rndm());
	//   double toyMetPhi=TMath::TwoPi()*(gRandom->Rndm());
	//   double toyMet_x=toyMet*cos(toyMetPhi);
	//   double toyMet_y=toyMet*sin(toyMetPhi);
	//   TVector2 MetToy(toyMet_x,toyMet_y);
	//   TVector2 MetToyEvnt(stuff.myMETcorr_th15.Px()+toyMet_x,stuff.myMETcorr_th15.Py()+toyMet_y);
	//   if(MyMetCleanUpCut(stuff,allstuff,stuff.Jet_lev6_noEMobj,MetToyEvnt)==0) Hist.fEvnt_toyMET_cut->Fill(MetToy.Mod());
	//   Hist.fEvnt_toyMET_all->Fill(MetToy.Mod());
	Hist.fEvnt_corMET_all->Fill(stuff.myMETcorr_th15.Mod());

	TLorentzVector jetsum(0.0,0.0,0.0,0.0);
	double ave_et=0.0;
	double ave_et0=0.0;
	Hist.fEvntNjet_b->Fill(stuff.myNjet);
	Hist.fEvntNjet5_b->Fill(stuff.myNjet_th5);
	Hist.fEvntNjet10_b->Fill(stuff.myNjet_th10);
	Hist.fEvntNjet15_b->Fill(stuff.myNjet_th15);
	Hist.fEvntNjet20_b->Fill(stuff.myNjet_th20);
	Hist.fEvntNjet25_b->Fill(stuff.myNjet_th25);
	Hist.fEvntNjet30_b->Fill(stuff.myNjet_th30);
	Hist.fEvntNjet35_b->Fill(stuff.myNjet_th35);
	Hist.fEvntdZ_b->Fill(dz);

	for(unsigned int i=0; i<stuff.Jet_raw_noEMobj.size(); i++)
	{
		if((stuff.Npho_match.at(i)+stuff.Nele_match.at(i))==0)
		{
			jetsum=jetsum+stuff.Jet_lev6_noEMobj.at(i);
			ave_et=ave_et+stuff.Jet_lev6_noEMobj.at(i).Pt();
			ave_et0=ave_et0+stuff.Jet_raw_noEMobj.at(i).Pt();
			Hist.fEvntEt0_Njet_b->Fill(stuff.myNjet,stuff.Jet_raw_noEMobj.at(i).Pt());
			Hist.fEvntEt_Njet_b->Fill(stuff.myNjet,stuff.Jet_lev6_noEMobj.at(i).Pt());
			if(i<2)
			{
				Hist.fEvntEt0_b[i]->Fill(stuff.Jet_raw_noEMobj.at(i).Pt());
				Hist.fEvntEt_b[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Pt());
				Hist.fEvntEtaDet_b[i]->Fill(stuff.EtaDet.at(i));
				Hist.fEvntEta_b[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Eta());
				Hist.fEvntPhi_b[i]->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj.at(i).Phi()));
				Hist.fEvntTheta_b[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Theta());
				Hist.fEvntEmFr_b[i]->Fill(stuff.EmFrRaw.at(i));
				Hist.fEvntNTowers_b[i]->Fill(stuff.JetNtwr.at(i));
				Hist.fEvntNTracks_b[i]->Fill(stuff.JetNtrk.at(i));
			}
			else
			{
				Hist.fEvntEt0X_b->Fill(stuff.Jet_raw_noEMobj.at(i).Pt());
				Hist.fEvntEtX_b->Fill(stuff.Jet_lev6_noEMobj.at(i).Pt());
				Hist.fEvntEtaDetX_b->Fill(stuff.EtaDet.at(i));
				Hist.fEvntEtaX_b->Fill(stuff.Jet_lev6_noEMobj.at(i).Eta());
				Hist.fEvntPhiX_b->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj.at(i).Phi()));
				Hist.fEvntThetaX_b->Fill(stuff.Jet_lev6_noEMobj.at(i).Theta());
				Hist.fEvntEmFrX_b->Fill(stuff.EmFrRaw.at(i));
				Hist.fEvntNTowersX_b->Fill(stuff.JetNtwr.at(i));
				Hist.fEvntNTracksX_b->Fill(stuff.JetNtrk.at(i));
				Hist.fEvntDeltaR1x_b->Fill(stuff.Jet_lev6_noEMobj.at(i).DeltaR(stuff.Jet_lev6[0]));
				Hist.fEvntDeltaR2x_b->Fill(stuff.Jet_lev6_noEMobj.at(i).DeltaR(stuff.Jet_lev6[1]));
			}
		}
	}
	if(stuff.Jet_raw_noEMobj.size()>1
			&& (stuff.Npho_match[0]+stuff.Nele_match[0])==0
			&& (stuff.Npho_match[1]+stuff.Nele_match[1])==0)
	{
		Hist.fEvntDeltaPhi_b->Fill(fabs(TVector2::Phi_mpi_pi(stuff.Jet_lev6_noEMobj[0].DeltaPhi(stuff.Jet_lev6_noEMobj[1]))));
		Hist.fEvntDeltaEta_b->Fill(stuff.Jet_lev6_noEMobj[0].Eta()-stuff.Jet_lev6_noEMobj[1].Eta());
		Hist.fEvntDeltaR_b->Fill(stuff.Jet_lev6_noEMobj[1].DeltaR(stuff.Jet_lev6_noEMobj[0]));
		Hist.fEvntMjj_b->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M());
		Hist.fEvntKtAll_b->Fill(jetsum.Pt());
		Hist.fEvntKtAll_Mjj_b->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),jetsum.Pt());
	}
	Hist.fEvntNJet_Nvx_b->Fill(nvx12,stuff.myNjet);
	Hist.fEvntNJet5_Nvx_b->Fill(nvx12,stuff.myNjet_th5);
	Hist.fEvntNJet10_Nvx_b->Fill(nvx12,stuff.myNjet_th10);
	Hist.fEvntNJet15_Nvx_b->Fill(nvx12,stuff.myNjet_th15);
	Hist.fEvntNJet20_Nvx_b->Fill(nvx12,stuff.myNjet_th20);
	Hist.fEvntNJet25_Nvx_b->Fill(nvx12,stuff.myNjet_th25);
	Hist.fEvntNJet30_Nvx_b->Fill(nvx12,stuff.myNjet_th30);
	Hist.fEvntNJet35_Nvx_b->Fill(nvx12,stuff.myNjet_th35);
	if(stuff.myNjet>0)
	{
		Hist.fEvntEt0_Nvx12_b->Fill(nvx12,ave_et0/(stuff.myNjet));
		Hist.fEvntEt_Nvx12_b->Fill(nvx12,ave_et/(stuff.myNjet));
	}

	if(nvx12==1) Hist.fEvntNjet_Lum_b->Fill((GetHeaderBlock()->InstLum())*1.0E-30,stuff.myNjet);
	return;
}
//_________________________________________ filling general histo for particular Jet Cone
void JetFilterModuleV3::FillJetHistogramsA(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12) {
	TLorentzVector jetsum(0.0,0.0,0.0,0.0);
	double ave_et=0.0;
	double ave_et0=0.0;
	Hist.fEvntNjet_a->Fill(stuff.myNjet);
	Hist.fEvntNjet5_a->Fill(stuff.myNjet_th5);
	Hist.fEvntNjet10_a->Fill(stuff.myNjet_th10);
	Hist.fEvntNjet15_a->Fill(stuff.myNjet_th15);
	Hist.fEvntNjet20_a->Fill(stuff.myNjet_th20);
	Hist.fEvntNjet25_a->Fill(stuff.myNjet_th25);
	Hist.fEvntNjet30_a->Fill(stuff.myNjet_th30);
	Hist.fEvntNjet35_a->Fill(stuff.myNjet_th35);
	Hist.fEvntdZ_a->Fill(dz);

	for(unsigned int i=0; i<stuff.Jet_raw_noEMobj.size(); i++)
	{
		if((stuff.Npho_match.at(i)+stuff.Nele_match.at(i))==0)
		{
			jetsum=jetsum+stuff.Jet_lev6_noEMobj.at(i);
			ave_et=ave_et+stuff.Jet_lev6_noEMobj.at(i).Pt();
			ave_et0=ave_et0+stuff.Jet_raw_noEMobj.at(i).Pt();
			Hist.fEvntEt0_Njet_a->Fill(stuff.myNjet,stuff.Jet_raw_noEMobj.at(i).Pt());
			Hist.fEvntEt_Njet_a->Fill(stuff.myNjet,stuff.Jet_lev6_noEMobj.at(i).Pt());
			if(i<2)
			{
				Hist.fEvntEt0_a[i]->Fill(stuff.Jet_raw_noEMobj.at(i).Pt());
				Hist.fEvntEt_a[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Pt());
				Hist.fEvntEtaDet_a[i]->Fill(stuff.EtaDet.at(i));
				Hist.fEvntEta_a[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Eta());
				Hist.fEvntPhi_a[i]->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj.at(i).Phi()));
				Hist.fEvntTheta_a[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Theta());
				Hist.fEvntEmFr_a[i]->Fill(stuff.EmFrRaw.at(i));
				Hist.fEvntNTowers_a[i]->Fill(stuff.JetNtwr.at(i));
				Hist.fEvntNTracks_a[i]->Fill(stuff.JetNtrk.at(i));
			}
			else
			{
				Hist.fEvntEt0X_a->Fill(stuff.Jet_raw_noEMobj.at(i).Pt());
				Hist.fEvntEtX_a->Fill(stuff.Jet_lev6_noEMobj.at(i).Pt());
				Hist.fEvntEtaDetX_a->Fill(stuff.EtaDet.at(i));
				Hist.fEvntEtaX_a->Fill(stuff.Jet_lev6_noEMobj.at(i).Eta());
				Hist.fEvntPhiX_a->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj.at(i).Phi()));
				Hist.fEvntThetaX_a->Fill(stuff.Jet_lev6_noEMobj.at(i).Theta());
				Hist.fEvntEmFrX_a->Fill(stuff.EmFrRaw.at(i));
				Hist.fEvntNTowersX_a->Fill(stuff.JetNtwr.at(i));
				Hist.fEvntNTracksX_a->Fill(stuff.JetNtrk.at(i));
				Hist.fEvntDeltaR1x_a->Fill(stuff.Jet_lev6_noEMobj.at(i).DeltaR(stuff.Jet_lev6_noEMobj[0]));
				Hist.fEvntDeltaR2x_a->Fill(stuff.Jet_lev6_noEMobj.at(i).DeltaR(stuff.Jet_lev6_noEMobj[1]));
			}
		}
	}
	if(stuff.Jet_raw_noEMobj.size()>1
			&& (stuff.Npho_match[0]+stuff.Nele_match[0])==0
			&& (stuff.Npho_match[1]+stuff.Nele_match[1])==0)
	{
		Hist.fEvntDeltaPhi_a->Fill(fabs(TVector2::Phi_mpi_pi(stuff.Jet_lev6_noEMobj[0].DeltaPhi(stuff.Jet_lev6_noEMobj[1]))));
		Hist.fEvntDeltaEta_a->Fill(stuff.Jet_lev6_noEMobj[0].Eta()-stuff.Jet_lev6_noEMobj[1].Eta());
		Hist.fEvntDeltaR_a->Fill(stuff.Jet_lev6_noEMobj[1].DeltaR(stuff.Jet_lev6_noEMobj[0]));
		Hist.fEvntMjj_a->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M());
		Hist.fEvntKtAll_a->Fill(jetsum.Pt());
		Hist.fEvntKtAll_Mjj_a->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),jetsum.Pt());
	}
	Hist.fEvntNJet_Nvx_a->Fill(nvx12,stuff.myNjet);
	Hist.fEvntNJet5_Nvx_a->Fill(nvx12,stuff.myNjet_th5);
	Hist.fEvntNJet10_Nvx_a->Fill(nvx12,stuff.myNjet_th10);
	Hist.fEvntNJet15_Nvx_a->Fill(nvx12,stuff.myNjet_th15);
	Hist.fEvntNJet20_Nvx_a->Fill(nvx12,stuff.myNjet_th20);
	Hist.fEvntNJet25_Nvx_a->Fill(nvx12,stuff.myNjet_th25);
	Hist.fEvntNJet30_Nvx_a->Fill(nvx12,stuff.myNjet_th30);
	Hist.fEvntNJet35_Nvx_a->Fill(nvx12,stuff.myNjet_th35);

	if(stuff.myNjet>0)
	{
		Hist.fEvntEt0_Nvx12_a->Fill(nvx12,ave_et0/(stuff.myNjet));
		Hist.fEvntEt_Nvx12_a->Fill(nvx12,ave_et/(stuff.myNjet));
	}
	if(nvx12==1) Hist.fEvntNjet_Lum_a->Fill((GetHeaderBlock()->InstLum())*1.0E-30,stuff.myNjet);
	return;
}

//______________________________________________________ filling match histo
void JetFilterModuleV3::FillMatchingHistograms(MatchStudyHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, MatchStuff match)
{
	//hist of EM objects matched and removed from jet list
	
	Hist.fMatchNtwr->Fill(jetstuff.JetNtwr[match.JetInd_match]);
	
	if (match.EmObjType_match==0)
	{
		Hist.fMatchDelR->Fill(jetstuff.Jet_raw[match.JetInd_match].DeltaR(miscstuff.myRawPhoton[match.EmInd_match]));
		Hist.fMatchDelPhi->Fill(fabs(TVector2::Phi_mpi_pi(jetstuff.Jet_raw[match.JetInd_match].Phi()-
						miscstuff.myRawPhoton[match.EmInd_match].Phi())));
		Hist.fMatchDelEta->Fill(fabs(jetstuff.Jet_raw[match.JetInd_match].Eta()-
					miscstuff.myRawPhoton[match.EmInd_match].Eta()));
		Hist.fMatchDelEtaDet->Fill(fabs(jetstuff.EtaDet[match.JetInd_match]-
					miscstuff.myPhoEtaDet[match.EmInd_match]));

		Hist.fMatchEtJet2EtPho_b->Fill(jetstuff.Jet_raw[match.JetInd_match].Pt()/miscstuff.myRawPhoton[match.EmInd_match].Pt());
		Hist.fMatchEtJet2EtPho_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].Pt()/miscstuff.myRawPhoton[match.EmInd_match].Pt());
		//       Hist.fMatchEtJet2EtPho_b->Fill(jetstuff.Jet_raw[match.JetInd_match].E()/miscstuff.myRawPhoton[match.EmInd_match].E());
		//       Hist.fMatchEtJet2EtPho_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].E()/miscstuff.myRawPhoton[match.EmInd_match].E());
	}
	if(match.EmObjType_match==1)
	{
		Hist.fMatchDelR->Fill(jetstuff.Jet_raw[match.JetInd_match].DeltaR(miscstuff.myRawElectron[match.EmInd_match]));
		Hist.fMatchDelPhi->Fill(fabs(TVector2::Phi_mpi_pi(jetstuff.Jet_raw[match.JetInd_match].Phi()-
						miscstuff.myRawElectron[match.EmInd_match].Phi())));
		Hist.fMatchDelEta->Fill(fabs(jetstuff.Jet_raw[match.JetInd_match].Eta()-
					miscstuff.myRawElectron[match.EmInd_match].Eta()));
		Hist.fMatchDelEtaDet->Fill(fabs(jetstuff.EtaDet[match.JetInd_match]-
					miscstuff.myEleEtaDet[match.EmInd_match]));

		Hist.fMatchEtJet2EtPho_b->Fill(jetstuff.Jet_raw[match.JetInd_match].Pt()/miscstuff.myRawElectron[match.EmInd_match].Pt());
		Hist.fMatchEtJet2EtPho_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].Pt()/miscstuff.myRawElectron[match.EmInd_match].Pt());
		//       Hist.fMatchEtJet2EtPho_b->Fill(jetstuff.Jet_raw[match.JetInd_match].E()/miscstuff.myRawElectron[match.EmInd_match].E());
		//       Hist.fMatchEtJet2EtPho_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].E()/miscstuff.myRawElectron[match.EmInd_match].E());
	}
	Hist.fMatchNmatch->Fill(jetstuff.Npho_match[match.JetInd_match]+jetstuff.Nele_match[match.JetInd_match]);
	Hist.fMatchEt_raw_b->Fill(jetstuff.Jet_raw[match.JetInd_match].Pt());
	Hist.fMatchEt_raw_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].Pt());
	Hist.fMatchEt_lev6_b->Fill(jetstuff.Jet_lev6[match.JetInd_match].Pt());
	Hist.fMatchEt_lev6_a->Fill(jetstuff.Jet_lev6_noEMobj[match.JetInd_match].Pt());
	if(jetstuff.Jet_raw_noEMobj[match.JetInd_match].Pt()>0.0)
	{
		Hist.fMatchDelRoldnew->Fill(jetstuff.Jet_raw[match.JetInd_match].DeltaR(jetstuff.Jet_raw_noEMobj[match.JetInd_match]));
		Hist.fMatchDelEtaDetoldnew->Fill(jetstuff.EtaDet[match.JetInd_match]-jetstuff.EtaDetCorr[match.JetInd_match]);
	}
	return;
}


//______________________________________ reads raw Met, pho, ele info and fills CommonStuff
void JetFilterModuleV3::DoCommonStuff(CommonStuff &miscstuff) {

	//_______________________________________talk to InitSuperPhotons and get Sam's photon info

	int _myNpho=initSpMod->GetSuperPhoSize();

	for(int i=0; i<_myNpho; i++)
	{

		//SetRemoveSidebandPhotons(1) or RemoveSidebandPhoton()
		// this does not work like this.
		// I need to move this selection to the point when the objects are
		//actually removed. not when they are stored in JetFilterV2.

		if ( (RemoveTightPho() && initSpMod->GetSuperPhoton(i)->IsTightPhoton() )
			 	|| (RemoveLoosePho() && initSpMod->GetSuperPhoton(i)->IsLoosePhoton()) 
			 //|| (RemoveSidebandPho() && ! initSpMod->GetSuperPhoton(i)->IsTightPhoton()  -->not complete yet
			 //    && initSpMod->GetSuperPhoton(i)->IsLoosePhoton())
			 )
		{
			if (PrintLevel()>10)
			{
				std::cout <<cyan << __FUNCTION__ << "::phos tight, loose = " 
							<< "\t" << initSpMod->GetSuperPhoton(i)->IsTightPhoton()
							<< "\t"<< initSpMod->GetSuperPhoton(i)->IsLoosePhoton() 
							<< clearatt << std::endl;
			}

			TLorentzVector _pho_raw = initSpMod->GetSuperPhoton(i)->GetRawVec();
			TLorentzVector _pho_cor = initSpMod->GetSuperPhoton(i)->GetCorVec();
			miscstuff.myPhoInd.push_back(initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex());
			miscstuff.myRawPhoton.push_back(_pho_raw);
			miscstuff.myCorrPhoton.push_back(_pho_cor);
			double hadem=initSpMod->GetSuperPhoton(i)->GetHadEm();
			miscstuff.myPhoEmFr.push_back(1.0/(hadem+1.0));
			miscstuff.myPhoEtaDet.push_back(initSpMod->GetSuperPhoton(i)->GetDetEta());
			miscstuff.myPhoXces.push_back(initSpMod->GetSuperPhoton(i)->GetXCes());
			miscstuff.myPhoZces.push_back(initSpMod->GetSuperPhoton(i)->GetZCes());
			miscstuff.tightId.push_back(initSpMod->GetSuperPhoton(i)->GetTightPhotonId());
		}	

	}


	int _myNele=initSpMod->GetSuperPhoSize();
	for(int i=0; i<_myNele; i++) {
		if ( (RemoveTightPhoLikeEle() && initSpMod->GetSuperPhoton(i)->IsTightElectron())
				|| (RemoveLoosePhoLikeEle() && initSpMod->GetSuperPhoton(i)->IsLooseElectron())
				// || (RemoveStdTightEle() && initSpMod->GetSuperPhoton(i)->IsStdTightElectron())		//not imp. yet
				|| (RemoveStdLooseEle() && initSpMod->GetSuperPhoton(i)->IsStdLooseElectron()) 
				 )
		{
			if (PrintLevel()>10)
			{
			std::cout << "index - eles tight, loose - Std.Tight, Std. Loose = " 
						<< i << " -\t" 
						<< initSpMod->GetSuperPhoton(i)->IsTightElectron()
						<< "\t"<< initSpMod->GetSuperPhoton(i)->IsLooseElectron() 
			//			<< "\t"<< initSpMod->GetSuperPhoton(i)->IsStdTightElectron() 
						<< "\t"<< initSpMod->GetSuperPhoton(i)->IsStdLooseElectron() 
						<< std::endl;
			}

			TLorentzVector _ele_raw = initSpMod->GetSuperPhoton(i)->GetRawVec();
			TLorentzVector _ele_cor = initSpMod->GetSuperPhoton(i)->GetCorVec();
			miscstuff.myEleInd.push_back(initSpMod->GetSuperPhoton(i)->GetElectronBlockIndex());
			miscstuff.myRawElectron.push_back(_ele_raw);
			miscstuff.myCorrElectron.push_back(_ele_cor);
			double hadem=initSpMod->GetSuperPhoton(i)->GetHadEm();
			miscstuff.myEleEmFr.push_back(1.0/(hadem+1.0));
			miscstuff.myEleEtaDet.push_back(initSpMod->GetSuperPhoton(i)->GetDetEta());
			miscstuff.myElePhiDet.push_back(initSpMod->GetSuperPhoton(i)->GetDetPhi());
			miscstuff.myEleXces.push_back(initSpMod->GetSuperPhoton(i)->GetXCes());
			miscstuff.myEleZces.push_back(initSpMod->GetSuperPhoton(i)->GetZCes());

			double factor=1.0; // added on 11/07/05
			double extra_PemEleEt=0.0; // extra PEM ele Et due to difference in energy scale, added on 11/07/04
			if (fabs(initSpMod->GetSuperPhoton(i)->GetDetEta())>1.1) // added on 11/07/05
			{
				if (fabs(initSpMod->GetSuperPhoton(i)->GetDetEta())<1.78) factor *= (initSpMod->GetSuperPhoton(i)->GetDetEta()>0.0? 1.020 : 1.015);
				else factor *= (initSpMod->GetSuperPhoton(i)->GetDetEta()>0.0? 1.010 : 1.007);
				extra_PemEleEt=(factor-1.0)*(_ele_raw.Pt());
			}
			else factor=0.0;
			// need to verify with sasha: as I use photon variables and he used electron variables.
			//double pprEt=extra_PemEleEt+factor*(MyZee->GetmyEleCprPpr(i))*fabs((MyZee->GetmyEleFid(i))*sin(_ele_raw->Theta())); // added on 11/04/05
			float Ele_Det = initSpMod->GetSuperPhoton(i)->GetDetector();
			double pprEt= -1.0;
			if (Ele_Det == 0)
			{
				extra_PemEleEt+factor*(fabs(initSpMod->GetSuperPhoton(i)->GetPhoton()->CprEnergy())*sin(_ele_raw.Theta())); // added on 11/04/05
			} else if (Ele_Det == 1)
			{
				extra_PemEleEt+factor*(fabs(initSpMod->GetSuperPhoton(i)->GetPhoton()->PprE())*sin(_ele_raw.Theta())); // added on 11/04/05

			}
			if(pprEt<0.0) pprEt=0.0; // added on 11/04/05
			miscstuff.myElePprEt.push_back(pprEt); // added on 11/04/05
		}
	}



	//_______________________________________ get my raw MET & SumET info

	miscstuff.mySumEt_raw=fMetBlock->Sumet(2);
	miscstuff.myMET_raw.Set(fMetBlock->MetX(4),fMetBlock->MetY(4)); //met calculated at highest sumPt vtx - sam
	miscstuff.myMET0_raw.Set(fMetBlock->MetX(0),fMetBlock->MetY(0)); //raw met calculated at z=0? - sam

	//____________________________ these are to be added later ______________________________
	//
	//     std::vector<TLorentzVector> myRawMuon;      // raw muons
	//     std::vector<TLorentzVector> myCorrMuon;     // corrected muons (corrections???)
	//     std::vector<TLorentzVector> myRawTau;       // raw taus
	//     std::vector<TLorentzVector> myCorrTau;      // corrected taus  (corrections???, for consistency)
	//     std::vector<TLorentzVector> myRawBjet;      // raw b-jets
	//     std::vector<TLorentzVector> myCorrBjet;     // corrected b-jets  (b-specific corrections???, for consistency)

	return;
}
//___________________________________ corrects Met & SumEt for jets
void JetFilterModuleV3::DoMyMet(CommonStuff miscstuff, JetStuff &jetstuff) {
	//-----correcting SumET
	jetstuff.mySumEtCorr_th5=miscstuff.mySumEt_raw;  // make it the same for a moment
	jetstuff.mySumEtCorr_th10=miscstuff.mySumEt_raw;
	jetstuff.mySumEtCorr_th15=miscstuff.mySumEt_raw;
	jetstuff.mySumEtCorr_th20=miscstuff.mySumEt_raw;
	jetstuff.mySumEtJet_th5=0.0;  // make it the same for a moment
	jetstuff.mySumEtJet_th10=0.0;
	jetstuff.mySumEtJet_th15=0.0;
	jetstuff.mySumEtJet_th20=0.0;
	//_________ correcting all SumEt for photons
	for(unsigned int i=0; i<miscstuff.myRawPhoton.size(); i++)
	{
		jetstuff.mySumEtCorr_th5=jetstuff.mySumEtCorr_th5-miscstuff.myRawPhoton.at(i).Pt();
		jetstuff.mySumEtCorr_th10=jetstuff.mySumEtCorr_th10-miscstuff.myRawPhoton.at(i).Pt();
		jetstuff.mySumEtCorr_th15=jetstuff.mySumEtCorr_th15-miscstuff.myRawPhoton.at(i).Pt();
		jetstuff.mySumEtCorr_th20=jetstuff.mySumEtCorr_th20-miscstuff.myRawPhoton.at(i).Pt();
	}
	//_________ correcting all SumEt for electrons
	for(unsigned int i=0; i<miscstuff.myRawElectron.size(); i++)
	{
		jetstuff.mySumEtCorr_th5=jetstuff.mySumEtCorr_th5-miscstuff.myRawElectron.at(i).Pt();
		jetstuff.mySumEtCorr_th10=jetstuff.mySumEtCorr_th10-miscstuff.myRawElectron.at(i).Pt();
		jetstuff.mySumEtCorr_th15=jetstuff.mySumEtCorr_th15-miscstuff.myRawElectron.at(i).Pt();
		jetstuff.mySumEtCorr_th20=jetstuff.mySumEtCorr_th20-miscstuff.myRawElectron.at(i).Pt();
	}
	//_________ correcting Met for photons and electrons, new as of 09/20/06
	TLorentzVector MetPho(0.0,0.0,0.0,0.0);
	TLorentzVector MetEle(0.0,0.0,0.0,0.0);
	for(unsigned int i=0; i<miscstuff.myRawPhoton.size(); i++)
	{
		MetPho=MetPho+miscstuff.myCorrPhoton.at(i)-miscstuff.myRawPhoton.at(i);
	}
	for(unsigned int i=0; i<miscstuff.myRawElectron.size(); i++)
	{
		MetEle=MetEle+miscstuff.myCorrElectron.at(i)-miscstuff.myRawElectron.at(i);
	}
	if(MetPho.E()<MetPho.P() || MetPho.E()<0.0) MetPho.SetPxPyPzE(0.0,0.0,0.0,0.0);
	if(MetEle.E()<MetEle.P() || MetEle.E()<0.0) MetEle.SetPxPyPzE(0.0,0.0,0.0,0.0);

	//_________ correcting all SumEt and calculating Met correction for jets
	TLorentzVector MetJet5(0.0,0.0,0.0,0.0);
	TLorentzVector MetJet10(0.0,0.0,0.0,0.0);
	TLorentzVector MetJet15(0.0,0.0,0.0,0.0);
	TLorentzVector MetJet20(0.0,0.0,0.0,0.0);
	for(unsigned int i=0; i<jetstuff.Jet_raw_noEMobj.size(); i++)
	{
		//--------- may need an eta cut here // this is ok as we treat everything (or leftover) as unclustered energy- sam
		if(fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())>=MinJetEta && fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<=MaxJetEta)
		{
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>5.0)
			{
				// corr=U.E.+M.I-raw=lev5-lev6-raw+lev1-lev4,
				jetstuff.mySumEtCorr_th5=jetstuff.mySumEtCorr_th5-jetstuff.Jet_raw_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				jetstuff.mySumEtJet_th5=jetstuff.mySumEtJet_th5+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt();
				MetJet5=MetJet5+jetstuff.Jet_lev5_noEMobj.at(i)
					+jetstuff.Jet_lev1_noEMobj.at(i)
					-jetstuff.Jet_lev4_noEMobj.at(i)
					-jetstuff.Jet_raw_noEMobj.at(i);
			}
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>10.0)
			{
				// corr=lev5-raw+M.I.=lev5-raw+lev1-lev4,
				jetstuff.mySumEtCorr_th10=jetstuff.mySumEtCorr_th10-jetstuff.Jet_raw_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				jetstuff.mySumEtJet_th10=jetstuff.mySumEtJet_th10+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt();
				MetJet10=MetJet10+jetstuff.Jet_lev5_noEMobj.at(i)
					+jetstuff.Jet_lev1_noEMobj.at(i)
					-jetstuff.Jet_lev4_noEMobj.at(i)
					-jetstuff.Jet_raw_noEMobj.at(i);
			}
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>15.0)
			{
				// corr=lev5-raw+M.I.=lev5-raw+lev1-lev4,
				jetstuff.mySumEtCorr_th15=jetstuff.mySumEtCorr_th15-jetstuff.Jet_raw_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				jetstuff.mySumEtJet_th15=jetstuff.mySumEtJet_th15+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt();
				MetJet15=MetJet15+jetstuff.Jet_lev5_noEMobj.at(i)
					+jetstuff.Jet_lev1_noEMobj.at(i)
					-jetstuff.Jet_lev4_noEMobj.at(i)
					-jetstuff.Jet_raw_noEMobj.at(i);
			}
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>20.0)
			{
				// corr=lev5-raw+M.I.=lev5-raw+lev1-lev4,
				jetstuff.mySumEtCorr_th20=jetstuff.mySumEtCorr_th20-jetstuff.Jet_raw_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				jetstuff.mySumEtJet_th20=jetstuff.mySumEtJet_th20+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt();
				MetJet20=MetJet20+jetstuff.Jet_lev5_noEMobj.at(i)
					+jetstuff.Jet_lev1_noEMobj.at(i)
					-jetstuff.Jet_lev4_noEMobj.at(i)
					-jetstuff.Jet_raw_noEMobj.at(i);
			}
		}
	}
	//---- making sure all MetJet makes sense:
	if(MetJet5.E()<MetJet5.P()) MetJet5.SetPxPyPzE(0.0,0.0,0.0,0.0);
	if(MetJet10.E()<MetJet10.P()) MetJet10.SetPxPyPzE(0.0,0.0,0.0,0.0);
	if(MetJet15.E()<MetJet15.P()) MetJet15.SetPxPyPzE(0.0,0.0,0.0,0.0);
	if(MetJet20.E()<MetJet20.P()) MetJet20.SetPxPyPzE(0.0,0.0,0.0,0.0);

	//---- this part was modified on 09/20/06
	MetJet5=MetJet5+MetEle+MetPho; // modified on 09/20/06
	MetJet10=MetJet10+MetEle+MetPho; // modified on 09/20/06
	MetJet15=MetJet15+MetEle+MetPho; // modified on 09/20/06
	MetJet20=MetJet20+MetEle+MetPho; // modified on 09/20/06

	jetstuff.myMETcorr_th5.Set(MetJet5.Px(),MetJet5.Py());
	jetstuff.myMETcorr_th10.Set(MetJet10.Px(),MetJet10.Py());
	jetstuff.myMETcorr_th15.Set(MetJet15.Px(),MetJet15.Py());
	jetstuff.myMETcorr_th20.Set(MetJet20.Px(),MetJet20.Py());
	//---- finally, calculating corrected MET
	jetstuff.myMETcorr_th5=miscstuff.myMET_raw-jetstuff.myMETcorr_th5;
	jetstuff.myMETcorr_th10=miscstuff.myMET_raw-jetstuff.myMETcorr_th10;
	jetstuff.myMETcorr_th15=miscstuff.myMET_raw-jetstuff.myMETcorr_th15;
	jetstuff.myMETcorr_th20=miscstuff.myMET_raw-jetstuff.myMETcorr_th20;

	return;
}
//__________________________________ does my jets; main routine
void JetFilterModuleV3::DoMyJet(TStnJetBlock* fJetBlock, CommonStuff miscstuff, JetStuff &jetstuff, MatchStudyHisto_t& Hist) {
	DoMyJetNoMatch(fJetBlock,jetstuff);
	DoMyJetWithMatch(fJetBlock,miscstuff,jetstuff,Hist);
	ReorderMyJets(jetstuff);
	//DumpJets(__FUNCTION__,__LINE__,jetstuff, 6);
	return;
}



//--the following three functions could have been replaced by one which works with generic data type (to be done)
void JetFilterModuleV3::myExchange_tlv(TLorentzVector& val1, TLorentzVector& val2) { //exchanges val1 and val2
	TLorentzVector dummy=val1;
	val1=val2;
	val2=dummy;
	return;
}
void JetFilterModuleV3::myExchange_dbl(double& val1, double& val2) { //exchanges val1 and val2
	double dummy=val1;
	val1=val2;
	val2=dummy;
	return;
}
void JetFilterModuleV3::myExchange_int(int& val1, int& val2) { //exchanges val1 and val2
	int dummy=val1;
	val1=val2;
	val2=dummy;
	return;
}
//_________________________________________ reorders jets after removing EM objects
void JetFilterModuleV3::ReorderMyJets(JetStuff &jetstuff) {
	int reorder_word=0;
	int reorder_word_last=0;
	for(int j=0; j<jetstuff.myNjet-1; j++)
	{
		reorder_word_last=reorder_word;
		for(int i=1; i<jetstuff.myNjet; i++)
		{
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>jetstuff.Jet_lev6_noEMobj.at(i-1).Pt())
			{
				myExchange_tlv(jetstuff.Jet_raw.at(i-1),jetstuff.Jet_raw.at(i));
				myExchange_tlv(jetstuff.Jet_lev1.at(i-1),jetstuff.Jet_lev1.at(i));
				myExchange_tlv(jetstuff.Jet_lev4.at(i-1),jetstuff.Jet_lev4.at(i));
				myExchange_tlv(jetstuff.Jet_lev5.at(i-1),jetstuff.Jet_lev5.at(i));
				myExchange_tlv(jetstuff.Jet_lev6.at(i-1),jetstuff.Jet_lev6.at(i));
				myExchange_tlv(jetstuff.Jet_lev7.at(i-1),jetstuff.Jet_lev7.at(i));
				myExchange_tlv(jetstuff.Jet_raw_noEMobj.at(i-1),jetstuff.Jet_raw_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev1_noEMobj.at(i-1),jetstuff.Jet_lev1_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev4_noEMobj.at(i-1),jetstuff.Jet_lev4_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev5_noEMobj.at(i-1),jetstuff.Jet_lev5_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev6_noEMobj.at(i-1),jetstuff.Jet_lev6_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev7_noEMobj.at(i-1),jetstuff.Jet_lev7_noEMobj.at(i));

				myExchange_dbl(jetstuff.EtaDet.at(i-1),jetstuff.EtaDet.at(i));
				myExchange_dbl(jetstuff.EtaDetCorr.at(i-1),jetstuff.EtaDetCorr.at(i));
				myExchange_dbl(jetstuff.EmFrRaw.at(i-1),jetstuff.EmFrRaw.at(i));
				myExchange_dbl(jetstuff.EmFrCorr.at(i-1),jetstuff.EmFrCorr.at(i));

				myExchange_int(jetstuff.JetNtrk.at(i-1),jetstuff.JetNtrk.at(i));
				myExchange_int(jetstuff.JetNtwr.at(i-1),jetstuff.JetNtwr.at(i));
				myExchange_int(jetstuff.Nobj_match.at(i-1),jetstuff.Nobj_match.at(i));
				myExchange_int(jetstuff.Npho_match.at(i-1),jetstuff.Npho_match.at(i));
				myExchange_int(jetstuff.Nele_match.at(i-1),jetstuff.Nele_match.at(i));
				myExchange_int(jetstuff.Nmu_match.at(i-1),jetstuff.Nmu_match.at(i));
				myExchange_int(jetstuff.Ntau_match.at(i-1),jetstuff.Ntau_match.at(i));
				myExchange_int(jetstuff.Nbtag_match.at(i-1),jetstuff.Nbtag_match.at(i));
				myExchange_int(jetstuff.JetBlockInd.at(i-1),jetstuff.JetBlockInd.at(i));

				reorder_word++;
			}
		}
		if(reorder_word_last==reorder_word) break;
	}

	//dump jets if the sorting failed (check only leading two jets for now)
	
	if (jetstuff.Jet_lev6_noEMobj.size()>1) 
	{

			//std::cout << "===================";
			//GetHeaderBlock()->Print();
			//for (int i=1; i<jetstuff.myNjet; i++)
			for (int i=1; i<jetstuff.Jet_lev6_noEMobj.size(); i++)
			{
				if (jetstuff.Jet_lev6_noEMobj.at(i-1).Pt() > jetstuff.Jet_lev6_noEMobj.at(i).Pt()) continue;
				
				
				if (PrintLevel()>10)
				{
				std::cout << __FUNCTION__ << ": Sorting of Jets has failed!  ";
				std::cout << "i[" << i << "] jet Et (i-1, i) = " << jetstuff.Jet_lev6_noEMobj.at(i).Pt(i-1) << ", " << jetstuff.Jet_lev6_noEMobj.at(i).Pt(i)<< std::endl;
				}
			}
			//std::cout << "===================" << std::endl;
			//assert (jetstuff.Jet_lev6_noEMobj.at(0).Pt() < jetstuff.Jet_lev6_noEMobj.at(1).Pt() &&
			//		"JetFilterModule::ReorderMyJets::Sorting of Jets has failed!");
	}
	return;
}

//____________________________________ re-calculates jet EmFr after removing EM object
double JetFilterModuleV3::CorrectEmFr(double emfr_old, double emfr_obj, double e_old, double e_obj)
{
	double emfr_new = emfr_old;
	
	if(e_old > e_obj)
	{
		double term1 = emfr_obj * e_obj / e_old;
		double denom = 1.0-e_obj / e_old;
		emfr_new = (emfr_old - term1) / denom;
	}
	return emfr_new;
}
//____________________________________ re-calculates jet EtaDet after removing EM object
//_____________ This function has been updated on 06/27/08.
//_____________ Implemented fix for large unphysical DelEta=etaOld-etaNew
double JetFilterModuleV3::CorrectEtaDet(TLorentzVector *vec_pho, TLorentzVector *vec_old,
										double jetetadet_old, double pho_etadet) 
{
	double E_jet 	= vec_old->E();
	double Pt_jet 	= vec_old->Pt();
	double Phi_jet = vec_old->Phi();
	
	TLorentzVector dummy_jet;
	dummy_jet.SetPtEtaPhiE(Pt_jet,jetetadet_old,Phi_jet,E_jet);
	
	double E_pho 	= vec_pho->E();
	double Pt_pho 	= vec_pho->Pt();
	double Phi_pho = vec_pho->Phi();
	
	TLorentzVector dummy_pho;
	dummy_pho.SetPtEtaPhiE(Pt_pho,pho_etadet,Phi_pho,E_pho);
	dummy_jet = dummy_jet - dummy_pho;
	
	double etadet_new = jetetadet_old;
	
	if (dummy_jet.E() > 0.0 && dummy_jet.Pt() > 0.0 && Pt_jet > Pt_pho) etadet_new = dummy_jet.Eta();
	if (fabs(etadet_new-jetetadet_old) > 0.8) etadet_new = jetetadet_old;

	return etadet_new;
}

//__________________________________ does my jets after removing EM objetcs, to be called after DoMyJetNoMatch
void JetFilterModuleV3::DoMyJetWithMatch(TStnJetBlock* fJetBlock, CommonStuff miscstuff, JetStuff &jetstuff, MatchStudyHisto_t& Hist) {

	int _jetcone=0;

	if(fabs(fJetBlock->ConeSize()-0.4)<0.01) _jetcone=0;
	if(fabs(fJetBlock->ConeSize()-0.7)<0.01) _jetcone=1;
	if(fabs(fJetBlock->ConeSize()-1.0)<0.01) _jetcone=2;
	fJTC_coneSize=_jetcone;

	//_________________________ Initializing Jet corrections
	CorInit lev1;
	CorInit lev4;
	CorInit lev5;
	CorInit lev6;
	CorInit lev7;
	lev1.level=1;
	lev1.nvx=myNvx_class12;
	lev1.cone=fJTC_coneSize;
	lev1.version=fJTC_version;
	lev1.sys=fJTC_systcode;
	lev1.imode=fJTC_imode;
	lev1.Nrun=myJetRun;
	lev4.level=4;
	lev4.nvx=myNvx_class12;
	lev4.cone=fJTC_coneSize;
	lev4.version=fJTC_version;
	lev4.sys=fJTC_systcode;
	lev4.imode=fJTC_imode;
	lev4.Nrun=myJetRun;
	lev5.level=5;
	lev5.nvx=myNvx_class12;
	lev5.cone=fJTC_coneSize;
	lev5.version=fJTC_version;
	lev5.sys=fJTC_systcode;
	lev5.imode=fJTC_imode;
	lev5.Nrun=myJetRun;
	lev6.level=6;
	lev6.nvx=myNvx_class12;
	lev6.cone=fJTC_coneSize;
	lev6.version=fJTC_version;
	lev6.sys=fJTC_systcode;
	lev6.imode=fJTC_imode;
	lev6.Nrun=myJetRun;
	lev7.level=7;
	lev7.nvx=myNvx_class12;
	lev7.cone=fJTC_coneSize;
	lev7.version=fJTC_version;
	lev7.sys=fJTC_systcode;
	lev7.imode=fJTC_imode;
	lev7.Nrun=myJetRun;

	MatchStuff dummymatch, sam_match;
	ClearMatchStuff();

	int _Npho_match=0;
	int _Nele_match=0;
	double epsilon=1.0E-10;

	if (PrintLevel()>10)
	{
		std::cout << __FUNCTION__ << "::" << __LINE__ << ":: njets=" << jetstuff.myNjet << std::endl;
	}

	for (int i=0; i<jetstuff.myNjet; i++)
	{
		float _etadet_tmp=jetstuff.EtaDet.at(i); // new line, 10/20/05
		TLorentzVector *_jet=&jetstuff.Jet_raw.at(i);
		//_____ removing photons
		for (unsigned int j=0; j<miscstuff.myRawPhoton.size() && _Npho_match<(int)miscstuff.myRawPhoton.size(); j++)
		{
			TLorentzVector *_pho=&miscstuff.myRawPhoton[j];
			double _match_dR;
			double _match_dEta;
			double _match_dPhi;
			int match_stat=MyMatchPhoJet(jetstuff.JetBlockInd.at(i),miscstuff.myPhoInd[j],-1,
					fJetBlock,_pho,_jet,_match_dR,_match_dPhi,_match_dEta);

			if (PrintLevel()>10)
			{
				std::cout << __FUNCTION__ << "::" << __LINE__ 
								<< ":: match_stat = " << match_stat << "\tjet i, pho j=" 
								<< i << ", " << j << std::endl;
			}

			if (match_stat==1)
			{
				_Npho_match++;
				jetstuff.Nobj_match.at(i)=jetstuff.Nobj_match.at(i)+1;
				jetstuff.Npho_match.at(i)=jetstuff.Npho_match.at(i)+1;

				//i added this to get matching info correctly, above stores a dummy value for EmInd_match
				// sam -- 03/15/2008
				sam_match.JetInd_match=jetstuff.JetBlockInd.at(i) ;
				sam_match.EmInd_match= miscstuff.myPhoInd[j];
				sam_match.EmObjType_match=0; // photons
				sam_matchstuff.push_back(sam_match);

				if ((jetstuff.Npho_match.at(i)+jetstuff.Nele_match.at(i))==1)
				{ // filling only first match
					dummymatch.JetInd_match=i;
					dummymatch.EmInd_match=j;
					dummymatch.EmObjType_match=0; // photons
					matchstuff.push_back(dummymatch);
				}

				// here, I re-calculate jet EmFr
				jetstuff.EmFrCorr.at(i) = CorrectEmFr(jetstuff.EmFrCorr.at(i),miscstuff.myPhoEmFr[j],
						jetstuff.Jet_raw_noEMobj.at(i).E(),miscstuff.myRawPhoton[j].E());
				
				// corrected EtaDet is to be calculated before removing EM object
				_etadet_tmp = CorrectEtaDet(_pho,&jetstuff.Jet_raw_noEMobj.at(i),
						jetstuff.EtaDetCorr.at(i),miscstuff.myPhoEtaDet[j]); // new line, 10/20/05

				jetstuff.EtaDetCorr.at(i) = _etadet_tmp;  // new line, 10/20/05
				double _jet_phi = jetstuff.Jet_raw.at(i).Phi();
				double _jet_eta = jetstuff.Jet_raw.at(i).Eta();
				if (jetstuff.Jet_raw_noEMobj.at(i).E()>miscstuff.myRawPhoton[j].E())
					jetstuff.Jet_raw_noEMobj.at(i)=jetstuff.Jet_raw_noEMobj.at(i)-miscstuff.myRawPhoton[j]; // removing raw photon from raw jet
				else jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);

				// updated from V3
				/*	else jetstuff.Jet_raw_noEMobj.at(i).SetPxPyPzE(0.0,0.0,0.0,0.0);

					if (jetstuff.Jet_raw_noEMobj.at(i).E()<jetstuff.Jet_raw_noEMobj.at(i).P()
					|| jetstuff.Jet_raw_noEMobj.at(i).E()<=0.0
					|| (jetstuff.Jet_raw.at(i).DeltaR(jetstuff.Jet_raw_noEMobj.at(i)))>2.0*(fJetBlock->ConeSize())) {

					jetstuff.Jet_raw_noEMobj.at(i).SetPxPyPzE(0.0,0.0,0.0,0.0); // in case removing EM object leads to M()<0.0
					}
				 */
				if (jetstuff.Jet_raw_noEMobj.at(i).Pt()<epsilon) jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon); // new, 05/28/08

				if (jetstuff.Jet_raw_noEMobj.at(i).E()<jetstuff.Jet_raw_noEMobj.at(i).P()
						|| jetstuff.Jet_raw_noEMobj.at(i).E()<epsilon)
					jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);

				if (jetstuff.Jet_raw.at(i).DeltaR(jetstuff.Jet_raw_noEMobj.at(i))>2.0*(fJetBlock->ConeSize()))
				{
					jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);
				}

			} // matchstat

		}	//photon loop

		//_____ removing electrons
		for (int j=0; j< (int) miscstuff.myRawElectron.size() && _Nele_match<(int) miscstuff.myRawElectron.size(); j++) {
			TLorentzVector *_ele=&miscstuff.myRawElectron[j];
			double _match_dR;
			double _match_dEta;
			double _match_dPhi;

			//std::cout << " processing ele i, ind" << i << "\t" << miscstuff.myEleInd[j] << std::endl;
			int match_stat=MyMatchPhoJet(jetstuff.JetBlockInd.at(i),-1,miscstuff.myEleInd[j],
					fJetBlock,_ele,_jet,_match_dR,_match_dPhi,_match_dEta);

			//std::cout << "jet i, delR, dEta, dPhi= " << i <<  "\t" << _match_dR << "\t" << _match_dEta << "," << _match_dPhi << std::endl;

			if (match_stat==1)
			{
				_Nele_match++;
				jetstuff.Nobj_match.at(i)=jetstuff.Nobj_match.at(i)+1;
				jetstuff.Nele_match.at(i)=jetstuff.Nele_match.at(i)+1;

				//i added this to get matching info correctly, above stores a dummy value for EmInd_match
				// sam -- 03/02/2008
				sam_match.JetInd_match=jetstuff.JetBlockInd.at(i);
				sam_match.EmInd_match= miscstuff.myEleInd[j];
				sam_match.EmObjType_match=1; // electrons
				sam_matchstuff.push_back(sam_match);

				if ((jetstuff.Nele_match.at(i)+jetstuff.Npho_match.at(i))==1)
				{ // filling only first match
					dummymatch.JetInd_match=i;
					dummymatch.EmInd_match=j;
					dummymatch.EmObjType_match=1; // electrons
					matchstuff.push_back(dummymatch);
				}

				// here, I re-calculate jet EmFr
				jetstuff.EmFrCorr.at(i)=CorrectEmFr(jetstuff.EmFrCorr.at(i),miscstuff.myEleEmFr[j],
						jetstuff.Jet_raw_noEMobj.at(i).E(),miscstuff.myRawElectron[j].E());
				// corrected EtaDet is to be calculated before removing EM object
				_etadet_tmp=CorrectEtaDet(_ele,&jetstuff.Jet_raw_noEMobj.at(i),
						jetstuff.EtaDetCorr.at(i),miscstuff.myEleEtaDet[j]); // new line, 10/20/05
				jetstuff.EtaDetCorr.at(i)=_etadet_tmp;  // new line, 10/20/05
				double _jet_phi=jetstuff.Jet_raw.at(i).Phi();
				double _jet_eta=jetstuff.Jet_raw.at(i).Eta();

				if (jetstuff.Jet_raw_noEMobj.at(i).E()>miscstuff.myRawElectron[j].E())
					jetstuff.Jet_raw_noEMobj.at(i)=jetstuff.Jet_raw_noEMobj.at(i)-miscstuff.myRawElectron[j]; // removing raw electron from raw jet
				else jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);

				if (jetstuff.Jet_raw_noEMobj.at(i).Pt()<epsilon) jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon); // new, 05/28/08
				if (jetstuff.Jet_raw_noEMobj.at(i).E()<jetstuff.Jet_raw_noEMobj.at(i).P()
						|| jetstuff.Jet_raw_noEMobj.at(i).E()<epsilon) jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);
				if(jetstuff.Jet_raw.at(i).DeltaR(jetstuff.Jet_raw_noEMobj.at(i))>2.0*(fJetBlock->ConeSize()))
				{
					jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);
				}

			}	//match stat
		} // electron loop


		if ((jetstuff.Nele_match.at(i)+jetstuff.Npho_match.at(i))>0) // correcting jet energy after removing EM objetcs
		{
			TLorentzVector _rawJet=jetstuff.Jet_raw_noEMobj.at(i);
			TLorentzVector _rawJet_tmp;
			_rawJet_tmp=_rawJet;
			float _myEmFr_tmp=jetstuff.EmFrCorr.at(i);
			double corr7=GetCorrection(lev7,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			_myEmFr_tmp=jetstuff.EmFrCorr.at(i); // do this again because EmFr can be potentially changed in GetCorrection
			_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
			double corr6=GetCorrection(lev6,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			_myEmFr_tmp=jetstuff.EmFrCorr.at(i); // do this again because EmFr can be potentially changed in GetCorrection
			_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
			double corr5=GetCorrection(lev5,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			_myEmFr_tmp=jetstuff.EmFrCorr.at(i); // do this again because EmFr can be potentially changed in GetCorrection
			_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
			double corr4=GetCorrection(lev4,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			_myEmFr_tmp=jetstuff.EmFrCorr.at(i); // do this again because EmFr can be potentially changed in GetCorrection
			_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
			double corr1=GetCorrection(lev1,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			if(corr1<0.0) corr1=epsilon;
			if(corr4<0.0) corr4=epsilon;
			if(corr5<0.0) corr5=epsilon;
			if(corr6<0.0) corr6=epsilon;
			if(corr7<0.0) corr7=epsilon;
			jetstuff.Jet_lev1_noEMobj.at(i)=_rawJet*corr1;
			jetstuff.Jet_lev4_noEMobj.at(i)=_rawJet*corr4;
			jetstuff.Jet_lev5_noEMobj.at(i)=_rawJet*corr5;
			jetstuff.Jet_lev6_noEMobj.at(i)=_rawJet*corr6;
			jetstuff.Jet_lev7_noEMobj.at(i)=_rawJet*corr7;
		}

		if(fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())>=MinJetEta && fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<=MaxJetEta)
		{
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>5.0) jetstuff.myNjet_th5++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>10.0) jetstuff.myNjet_th10++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>15.0) jetstuff.myNjet_th15++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>20.0) jetstuff.myNjet_th20++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>25.0) jetstuff.myNjet_th25++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>30.0) jetstuff.myNjet_th30++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>35.0) jetstuff.myNjet_th35++;
		}
	} // jet loop

	if(_Nele_match!=(int)miscstuff.myRawElectron.size() || _Npho_match!=(int)miscstuff.myRawPhoton.size()) bad_EMjet_match_flag=0;
	for(unsigned int i=0; i<matchstuff.size(); i++)
	{
		FillMatchingHistograms(Hist,jetstuff,miscstuff,matchstuff[i]);
	}

	if (_Npho_match != (int)miscstuff.myRawPhoton.size())
	{
		std::cout << __FILE__ << ":" << __LINE__ << "::JET:: PHO NO MATCH: Npho, Nmatched ="
			<< miscstuff.myRawPhoton.size() << ", " << _Npho_match << "\t- ";
		GetHeaderBlock()->Print();
	}
	if (_Nele_match != (int) miscstuff.myRawElectron.size())
	{
		std::cout << __FILE__ << ":" << __LINE__ << "::JET:: ELE NO MATCH: Nele, Nmatched ="
			<< miscstuff.myRawElectron.size() << ", " << _Nele_match << "\t- ";
		GetHeaderBlock()->Print();
	}


	return;
}
//__________________________________ does my jets, but doesn't remove EM objetcs
void JetFilterModuleV3::DoMyJetNoMatch(TStnJetBlock* fJetBlock, JetStuff &jetstuff) {

	int _jetcone=0;
	if(fabs(fJetBlock->ConeSize()-0.4)<0.01) _jetcone=0;
	if(fabs(fJetBlock->ConeSize()-0.7)<0.01) _jetcone=1;
	if(fabs(fJetBlock->ConeSize()-1.0)<0.01) _jetcone=2;

	fJTC_coneSize=_jetcone;
	//_________________________ Initializing Jet corrections
	CorInit lev1;
	CorInit lev4;
	CorInit lev5;
	CorInit lev6;
	CorInit lev7;
	lev1.level=1;
	lev1.nvx=myNvx_class12;
	lev1.cone=fJTC_coneSize;
	lev1.version=fJTC_version;
	lev1.sys=fJTC_systcode;
	lev1.imode=fJTC_imode;
	lev1.Nrun=myJetRun;
	lev4.level=4;
	lev4.nvx=myNvx_class12;
	lev4.cone=fJTC_coneSize;
	lev4.version=fJTC_version;
	lev4.sys=fJTC_systcode;
	lev4.imode=fJTC_imode;
	lev4.Nrun=myJetRun;
	lev5.level=5;
	lev5.nvx=myNvx_class12;
	lev5.cone=fJTC_coneSize;
	lev5.version=fJTC_version;
	lev5.sys=fJTC_systcode;
	lev5.imode=fJTC_imode;
	lev5.Nrun=myJetRun;
	lev6.level=6;
	lev6.nvx=myNvx_class12;
	lev6.cone=fJTC_coneSize;
	lev6.version=fJTC_version;
	lev6.sys=fJTC_systcode;
	lev6.imode=fJTC_imode;
	lev6.Nrun=myJetRun;
	lev7.level=7;
	lev7.nvx=myNvx_class12;
	lev7.cone=fJTC_coneSize;
	lev7.version=fJTC_version;
	lev7.sys=fJTC_systcode;
	lev7.imode=fJTC_imode;
	lev7.Nrun=myJetRun;

	jetstuff.myNjet=fJetBlock->NJets();
	if (PrintLevel()>10)	
	{
		std::cout << __FUNCTION__ << ":: jets found = " << jetstuff.myNjet << std::endl;
		std::cout << std::setw(3) << "i" << std::setw(10) <<  "Et" << std::setw(10) << "E" << std::setw(10) << "Eta" << std::setw(10) << "JetBlkInd" << std::endl;
	}
	for(int i=0; i<jetstuff.myNjet; i++)
	{
		TStnJet* jet = fJetBlock->Jet(i);
		jetstuff.JetBlockInd.push_back(i);
		jetstuff.JetNtrk.push_back(jet->NTracks());
		jetstuff.JetNtwr.push_back(jet->NTowers());
		jetstuff.EtaDet.push_back(jet->DetEta());
		jetstuff.EtaDetCorr.push_back(jet->DetEta()); // for a moment, make it the same as "EtaDet"
		jetstuff.EmFrRaw.push_back(jet->Emfr());
		jetstuff.EmFrCorr.push_back(jet->Emfr()); // for a moment, make it the same as "EmFrRaw"

		jetstuff.Nobj_match.push_back(0); // these are packed with zero's for a moment
		jetstuff.Npho_match.push_back(0);
		jetstuff.Nele_match.push_back(0);
		jetstuff.Nmu_match.push_back(0);
		jetstuff.Ntau_match.push_back(0);
		jetstuff.Nbtag_match.push_back(0);


		//__________________ getting raw jets
		TLorentzVector _rawJet;
		TLorentzVector _rawJet_tmp;
		_rawJet.SetPx(jet->Momentum()->Px());
		_rawJet.SetPy(jet->Momentum()->Py());
		_rawJet.SetPz(jet->Momentum()->Pz());
		_rawJet.SetE(jet->Momentum()->E());
		jetstuff.Jet_raw.push_back(_rawJet);
		jetstuff.Jet_raw_noEMobj.push_back(_rawJet); // for a moment, make it the same as "Jet_raw"

		if (PrintLevel()>10)	std::cout << std::setw(3) << i << std::setw(10) <<  _rawJet.Pt() << std::setw(10) << _rawJet.E() << std::setw(10) << jet->DetEta() << std::setw(10) << i << std::endl;
		
		float _myEmFr_tmp=jet->Emfr();
		float _etadet_tmp=jet->DetEta();
		_rawJet_tmp=_rawJet;
		double corr7=GetCorrection(lev7,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
		_myEmFr_tmp=jet->Emfr(); // do this again because EmFr can be potentially changed in GetCorrection
		_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
		double corr6=GetCorrection(lev6,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
		_myEmFr_tmp=jet->Emfr(); // do this again because EmFr can be potentially changed in GetCorrection
		_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
		double corr5=GetCorrection(lev5,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
		_myEmFr_tmp=jet->Emfr(); // do this again because EmFr can be potentially changed in GetCorrection
		_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
		double corr4=GetCorrection(lev4,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
		_myEmFr_tmp=jet->Emfr(); // do this again because EmFr can be potentially changed in GetCorrection
		_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
		double corr1=GetCorrection(lev1,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);

		if(corr1<0.0) corr1=0.0;
		if(corr4<0.0) corr4=0.0;
		if(corr5<0.0) corr5=0.0;
		if(corr6<0.0) corr6=0.0;
		if(corr7<0.0) corr7=0.0;
		jetstuff.Jet_lev1.push_back(_rawJet*corr1);
		jetstuff.Jet_lev1_noEMobj.push_back(_rawJet*corr1); // for a moment, make it the same as "Jet_lev1"
		jetstuff.Jet_lev4.push_back(_rawJet*corr4);
		jetstuff.Jet_lev4_noEMobj.push_back(_rawJet*corr4); // for a moment, make it the same as "Jet_lev4"
		jetstuff.Jet_lev5.push_back(_rawJet*corr5);
		jetstuff.Jet_lev5_noEMobj.push_back(_rawJet*corr5); // for a moment, make it the same as "Jet_lev5"
		jetstuff.Jet_lev6.push_back(_rawJet*corr6);
		jetstuff.Jet_lev6_noEMobj.push_back(_rawJet*corr6); // for a moment, make it the same as "Jet_lev6"
		jetstuff.Jet_lev7.push_back(_rawJet*corr7);
		jetstuff.Jet_lev7_noEMobj.push_back(_rawJet*corr7); // for a moment, make it the same as "Jet_lev7"

	}
	return;
}


//___ returns relative position of EM shower inside CEM tower
double JetFilterModuleV3::MyInCemTowerEta(double eta_det) {
	double rel_eta=-10.0;
	double theta_boundary[12]={90.0,82.526,75.297,68.516,62.310,56.735,51.790,47.436,43.614,40.261,36.822,33.524};
	double eta_l, eta_r, d_eta;
	double _pi=TMath::Pi();
	for(int i=1; i<12; i++)
	{
		eta_l=fabs(log(tan(_pi*theta_boundary[i-1]/(2.0*180.0))));
		eta_r=fabs(log(tan(_pi*theta_boundary[i]/(2.0*180.0))));
		d_eta=eta_r-eta_l;
		if(fabs(eta_det)>eta_l && fabs(eta_det)<=eta_r)
		{
			rel_eta=fabs(fabs(eta_det)-eta_l)/d_eta;
			break;
		}
	}
	return rel_eta;
}



//________ creates a list of towers for all jets
void JetFilterModuleV3::MatchCalorTowers(int jet_ind,
		TStnJetBlock* fJetBlock,
		TCalDataBlock *fCalDataBlock,
		CalDataArray* cdholder) {
	cdholder->clear();
	TStnLinkBlock* links = fJetBlock->TowerLinkList();
	int nptow = links->NLinks(jet_ind);

	TCalTower* ctower = NULL;
	for(int j=0; j<nptow; j++)
	{
		int iptow = links->Index(jet_ind,j);
		int iphi = iptow & 0x3F;
		int ieta = (iptow>>8) & 0x3F;

		ctower = fCalDataBlock->Tower(ieta,iphi);
		cdholder->push_back(ctower);
	}

	std::sort(cdholder->begin(), cdholder->end(),SortTowersByEnergy);

	return;
}

//________ creates a list of towers for all EM objects
void JetFilterModuleV3::MatchCalorTowers(int pho_ind,
		TStnPhotonBlock* fPhotonBlock,
		TCalDataBlock *fCalDataBlock,
		CalDataArray* cdholder) {
	cdholder->clear();
	TStnLinkBlock* links = fPhotonBlock->TowerLinkList();
	int nptow = links->NLinks(pho_ind);
	TCalTower* ctower = NULL;
	for(int j=0; j<nptow; j++)
	{
		int iptow = links->Index(pho_ind,j);
		int iphi = iptow & 0x3F;
		int ieta = (iptow>>8) & 0x3F;

		ctower = fCalDataBlock->Tower(ieta,iphi);
		cdholder->push_back(ctower);
	}

	std::sort(cdholder->begin(),cdholder->end(),SortTowersByEnergy);

	return;
}

//________ creates a list of towers for electrons
void JetFilterModuleV3::MatchCalorTowers(int ele_ind,
		TStnElectronBlock* fElectronBlock,
		TCalDataBlock *fCalDataBlock,
		CalDataArray* cdholder) {
	cdholder->clear();
	TStnLinkBlock* links = fElectronBlock->TowerLinkList();
	int nptow = links->NLinks(ele_ind);
	TCalTower* ctower = NULL;
	for(int j=0; j<nptow; j++)
	{
		int iptow = links->Index(ele_ind,j);
		int iphi = iptow & 0x3F;
		int ieta = (iptow>>8) & 0x3F;

		ctower = fCalDataBlock->Tower(ieta,iphi);
		cdholder->push_back(ctower);
	}

	std::sort(cdholder->begin(),cdholder->end(),SortTowersByEnergy);

	return;
}

//_____________________________________________________________________________
int JetFilterModuleV3::BeginJob() {


	//_____________________________________________________ register the data block
	RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlockClu04);
	RegisterDataBlock("PROD@JetCluModule-had-cone0.4","TStnJetBlock",&fHadJetBlockClu04);

	//---------- need this for jet-EM matching
	RegisterDataBlock("CalDataBlock","TCalDataBlock",&fCalData);
	RegisterDataBlock("PhotonBlock","TStnPhotonBlock",&fPhotonBlock);
	RegisterDataBlock("ElectronBlock","TStnElectronBlock",&fElectronBlock);
	RegisterDataBlock("MetBlock","TStnMetBlock",&fMetBlock);		//added to replace Sasha's EventFiltermod 03-14-2008
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");

	//_____________________________________________________ book histograms
	BookHistograms();

	EventCount_b=0;
	EventCount_a=0;

	// sam's stuff -- 03-14-2008
	mycounter.evtsRunOver = 0;
	mycounter.evtsPassModule = 0;
	mycounter.evtsMoreHad = 0;

	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initSpMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
		bRunPermit = false;
	}
	tightMod = (TagTightPhotons*) ((TStnAna*) GetAna()->GetModule("TagTightPhotons"));
	if (tightMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TagTightPhotons!");
		bRunPermit = false;
	}
	looseMod = (TagLoosePhotons*) ((TStnAna*) GetAna()->GetModule("TagLoosePhotons"));
	if (looseMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TagLoosePhotons!");
		bRunPermit = false;
	}
	tightEleMod = (TagElectrons*) ((TStnAna*) GetAna()->GetModule("TagElectrons"));
	if (tightEleMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TagElectrons!");
		bRunPermit = false;
	}
	looseEleMod = (TagLooseElectrons*) ((TStnAna*) GetAna()->GetModule("TagLooseElectrons"));
	if (looseEleMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TagElectrons!");
		bRunPermit = false;
	}
	trigMod = (TriggerModule*) ((TStnAna*) GetAna()->GetModule("TriggerModule"));
	if (trigMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TriggerModule!");
		bRunPermit = false;
	}

	//removed these from Event loop
	if (!fPhotonBlock) {
		StdOut(__FILE__,__LINE__,3,"PhotonBlock not found!");
		bRunPermit = false;
		return 0;
	}
	if (!fElectronBlock) {
		StdOut(__FILE__,__LINE__,3,"ElectronBlock not found!");
		bRunPermit = false;
		return 0;
	}

	if (!fCalData) {
		StdOut(__FILE__,__LINE__,3,"CalDataBlock Block not found!");
		bRunPermit = false;
		return 0;
	}
	if (!fMetBlock)	{
		StdOut(__FILE__,__LINE__,3,"MetBlock not found!");
		bRunPermit = false;
		return 0;
	}
	if (!fJetBlockClu04)	{
		StdOut(__FILE__,__LINE__,3,"JetBlock not found!");
		bRunPermit = false;
		return 0;
	}

	return 0;
}


//_____________________________________________________________________________
int JetFilterModuleV3::BeginRun() {

	fJTC_imode = ! (fHeaderBlock->McFlag());		// fJTC_imode = 1(DATA), 0 (MC)
	std::cout << __FILE__ << "::" << __LINE__ << ":AUTOMATIC SETTING OF MC_FLAG = " << fHeaderBlock->McFlag() << " FOR JTCmodei = " << fJTC_imode << "(1=DATA, 0=MC)"  << std::endl;
	return 0;
}

//_____________________________________________________________________________
int JetFilterModuleV3::Event(int ientry)
{
	SetPassed(0); 
	mycounter.evtsRunOver++;
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"One or more dependencies not found. pl check. Run Permit not cleared!");
		SetPassed(0);
		exit (1);
		return 0;
	}


	bSaveThisEvent = false;
	EventCount_b++;
	ClearModuleOutput();

	myJetRun=GetHeaderBlock()->RunNumber(); // obtaining run number, myJetRun is globaly defined
	if (debug)
	{
		std::cout << green << "================" << myJetRun << ", " << GetHeaderBlock()->EventNumber() << clearatt << std::endl; 
	}

	//_________________________________________________________________________
	//---------  Accessing Jet Blocks
	fJetBlockClu04->GetEntry(ientry);
	fHadJetBlockClu04->GetEntry(ientry);
	fPhotonBlock->GetEntry(ientry);
	fCalData->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);	// sam 03-14-2008
	fGenpBlock->GetEntry(ientry);

	myNvx_class12     = trigMod->GetNClass12Vtx();

	DoCommonStuff(allstuff); // reads raw Met, pho, ele info and fills CommonStuff
	DoMyJet(fJetBlockClu04,allstuff,jet04stuff,fMatchPhoJet04); // does my JetClu04 jets
	DoMyMet(allstuff,jet04stuff); // corrects Met & SumEt for JetClu04 jets
	//temporary to study jetfilerv3
	float dz=0;
	FillJetHistogramsB(fHistJet04,jet04stuff,dz ,myNvx_class12); // filling general histo for JetClu04

//	const float jf3_Ht = GetMyHtCorr(0,3,6);
//	const float jf3_Met = GetMyMetCorr(0,3);
/*	if (fabs(jf3_Ht - jf2_Ht) > 0.1 || fabs(jf3_Ht - jf2_Ht) > 0.1)
	{
	std::cout << green << __FILE__ << ":";
	fHeaderBlock->Print();
		std::cout << "Ht2/Met2 = " << jf2_Ht << ", " << jf2_Met << std::endl;
		std::cout << "Ht3/Met3 = " << jf3_Ht << ", " << jf3_Met << std::endl;
		DumpEMobjects(__FUNCTION__, __LINE__, allstuff, 0);
		DumpJets(__FUNCTION__,__LINE__,jet04stuff,16);
		std::cout << clearatt << std::endl;
	}
	
*/		// I want to pass many events as possible when making flat ntuples
		// I'll pass events only if there is a jet with
		// Et>MinEt and Eta<MaxEta -- Sam, Sep 7,2009
		
		if (MyNjetCut(jet04stuff)==1)
		{

			for (unsigned int i=0; i<jet04stuff.Jet_lev6_noEMobj.size(); i++)
			{
				if (fabs(jet04stuff.Jet_lev6_noEMobj.at(i).Eta())<=MaxJetEta)
				{
					FillJetHistogramsA(fHistJet04,jet04stuff,dz ,myNvx_class12); // filling general histo for JetClu04
					SetPassed(1);
					break;
				}
			}
		}


	if (GetPassed()) 	mycounter.evtsPassModule++;


	//std::cout << "======= END EVENT - ";
	return 0;
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int JetFilterModuleV3::EndJob() {

	if (GetSummaryStat()) return 0;

	std::string sMsg;
	if (GetJTC_imode() == 0) sMsg = "MC";
	else if (GetJTC_imode() == 1) sMsg = "DATA";
	else sMsg = "UNKNOWN";


	std::string sJTCsys;
	if (GetJTC_systcode() == 0) sJTCsys = "(0=default)";
	else if (GetJTC_systcode() == 1) sJTCsys = "(+1 sigma(JES))";
	else if (GetJTC_systcode() == -1) sJTCsys = "(-1 sigma(JES))";
	else sJTCsys = "(UNKNOWN)";
	
	
	printf("[TMJ:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[TMJ:01:] Events Processed ----------- = " << mycounter.evtsRunOver << std::endl;
	std::cout << "[TMJ:02:] Events Passed -------------- = " << mycounter.evtsPassModule << std::endl;
	
	std::string njetcut("");
	if (GetMinNjet5()>0) njetcut += "[5 GeV, >=" + ToStr(GetMinNjet5()) + " & =<" +  ToStr(GetMaxNjet5()) + "] ";
	if (GetMinNjet10()>0) njetcut += "[10 GeV, >=" + ToStr(GetMinNjet10()) + " & =<" +  ToStr(GetMaxNjet10()) + "] ";
	if (GetMinNjet15()>0) njetcut += "[15 GeV, >=" + ToStr(GetMinNjet15()) + " & =<" +  ToStr(GetMaxNjet15()) + "] ";
	if (GetMinNjet20()>0) njetcut += "[20 GeV, >=" + ToStr(GetMinNjet20()) + " & =<" +  ToStr(GetMaxNjet20()) + "] ";
	if (GetMinNjet25()>0) njetcut += "[25 GeV, >=" + ToStr(GetMinNjet25()) + " & =<" +  ToStr(GetMaxNjet25()) + "] ";
	if (GetMinNjet30()>0) njetcut += "[30 GeV, >=" + ToStr(GetMinNjet30()) + " & =<" +  ToStr(GetMaxNjet30()) + "] ";
	if (GetMinNjet35()>0) njetcut += "[35 GeV, >=" + ToStr(GetMinNjet35()) + " & =<" +  ToStr(GetMaxNjet35()) + "] ";
	if (GetMinNjet5()==0 && GetMinNjet10()==0 && GetMinNjet15()==0 && GetMinNjet20()==0 && GetMinNjet25()==0
			&& GetMinNjet30()==0 && GetMinNjet35() ==0) njetcut += "No Jet Cut";
	
	std::cout << "[TMJ:03:] (Et> - Min/Max) NJets ------ = " << njetcut << std::endl;
	std::cout << "[TMJ:04:] Max Jet Detector Eta         = " << GetMaxJetEta() <<std::endl;
	std::cout << "[TMJ:08:] JTC_imode / Syst Code ------ = " << GetJTC_imode() << " (" << sMsg 
				<< ") / " << GetJTC_systcode() << sJTCsys << std::endl;


		std::string sEmObj;
		if (RemoveTightPho()) sEmObj += "[Tight Pho]";
		if (RemoveLoosePho()) sEmObj += "[Loose Pho]";
		//if (RemoveSidebandPho()) sEmObj += "[Sideband Pho]";
		if (RemoveTightPhoLikeEle()) sEmObj += "[Tight PhoLikeEle]";
		if (RemoveLoosePhoLikeEle()) sEmObj += "[Loose PhoLikeEle]";
		if (RemoveStdLooseEle()) sEmObj += "[Std. Loose Ele]";
		if (! (RemoveTightPho() || RemoveLoosePho() || RemoveTightPhoLikeEle()
				|| RemoveLoosePhoLikeEle() || RemoveStdLooseEle()) ) sEmObj += "None";
		std::cout << "[TMJ:30:] Types of EM objs removed     = " << sEmObj << std::endl;


	std::cout << "[TMJ:32:] MinEt1st                     = " << GetMinEt1st() << std::endl;
	std::cout << "[TMJ:33:] MinEt2nd                     = " << GetMinEt2nd() << std::endl;
	std::cout << "[TMJ:34:] Min/Max EtThr for Ht/MEt calc= " << GetMinEtThr() << ", " << GetMaxEtThr()  << " (this is a cur on jet when calcuating Ht/MET, def=15Gev)"<< std::endl; 

	
	printf("---------------------------------------------------\n");

	return 0;
}

double JetFilterModuleV3::GetCorrection(CorInit settings, TLorentzVector vec,
		float& emf, float etad) {
	int nrun = settings.Nrun;
	int nVertex = settings.nvx;
	int coneSize=settings.cone;
	int version=settings.version;
	int syscode=settings.sys;
	int level=settings.level;
	int imode=settings.imode;

	JetEnergyCorrections myJetEnergyCorrections=
		JetEnergyCorrections("JetCorrections","JetCorrections",
				level,nVertex,coneSize,version,syscode,nrun,imode);

	HepLorentzVector P4Jet;
	P4Jet.setPx(vec.Px());
	P4Jet.setPy(vec.Py());
	P4Jet.setPz(vec.Pz());
	P4Jet.setE(vec.E());
	if(abs(syscode)>0)
	{
		int n_sigma= (syscode>0) ? 1 : -1 ;
		myJetEnergyCorrections.setTotalSysUncertainties(n_sigma);
	}
	double scaleFactor = myJetEnergyCorrections.doEnergyCorrections(P4Jet,emf,etad);
	return scaleFactor;
}

double JetFilterModuleV3::GetCorrection(TLorentzVector vec, float& emf, float etad) {

	int nrun = myJetRun;
	int nVertex = myNvx_class12;
	int coneSize=fJTC_coneSize;
	int version=fJTC_version;
	int syscode=fJTC_systcode;
	int level=fJTC_level;
	int imode=fJTC_imode;

	JetEnergyCorrections myJetEnergyCorrections=
		JetEnergyCorrections("JetCorrections","JetCorrections",
				level,nVertex,coneSize,version,syscode,nrun,imode);

	HepLorentzVector P4Jet;
	P4Jet.setPx(vec.Px());
	P4Jet.setPy(vec.Py());
	P4Jet.setPz(vec.Pz());
	P4Jet.setE(vec.E());
	if(abs(syscode)>0)
	{
		int n_sigma= (syscode>0) ? 1 : -1 ;
		myJetEnergyCorrections.setTotalSysUncertainties(n_sigma);
	}
	double scaleFactor = myJetEnergyCorrections.doEnergyCorrections(P4Jet,emf,etad);
	return scaleFactor;
}

//__________________________________________________ returns 1 if jet is matched
//                                                   calculates dR, dPhi, dEta between pho & jet
int JetFilterModuleV3::MyMatchPhoJet(int jet_ind, int pho_ind, int ele_ind,TStnJetBlock* fJetBlock,
		TLorentzVector *mypho,TLorentzVector *myjet,
		double &dR_fj, double &dPhi_fj, double &dEta_fj)
{
	int match_code=0;

	if(jet_ind<0) return match_code;
	if(pho_ind<0 && ele_ind<0) return match_code;

	if(pho_ind>=0)
	{
		CalDataArray jet_towers;
		CalDataArray pho_towers;
		MatchCalorTowers(jet_ind,fJetBlock,fCalData,&jet_towers);
		MatchCalorTowers(pho_ind,fPhotonBlock,fCalData,&pho_towers);
		CalDataArrayI tt_em;
		CalDataArrayI tt_jet;
		for(tt_em = pho_towers.begin(); tt_em != pho_towers.end(); tt_em++)
		{
			TCalTower* calt_em = *tt_em;
			int iEta_em = calt_em->IEta();
			int iPhi_em = calt_em->IPhi();
			for(tt_jet = jet_towers.begin(); tt_jet != jet_towers.end(); tt_jet++)
			{
				TCalTower* calt_jet = *tt_jet;
				int iEta_jet = calt_jet->IEta();
				int iPhi_jet = calt_jet->IPhi();
				if(iEta_em==iEta_jet && iPhi_em==iPhi_jet) match_code=1;
				if(match_code==1) break;
			}
			if(match_code==1) break;
		}
	}

	if(ele_ind>=0)
	{
		CalDataArray jet_towers;
		CalDataArray ele_towers;
		MatchCalorTowers(jet_ind,fJetBlock,fCalData,&jet_towers);
		MatchCalorTowers(ele_ind,fElectronBlock,fCalData,&ele_towers);
		CalDataArrayI tt_em;
		CalDataArrayI tt_jet;
		for(tt_em = ele_towers.begin(); tt_em != ele_towers.end(); tt_em++)
		{
			TCalTower* calt_em = *tt_em;
			int iEta_em = calt_em->IEta();
			int iPhi_em = calt_em->IPhi();
			for(tt_jet = jet_towers.begin(); tt_jet != jet_towers.end(); tt_jet++)
			{
				TCalTower* calt_jet = *tt_jet;
				int iEta_jet = calt_jet->IEta();
				int iPhi_jet = calt_jet->IPhi();
				if(iEta_em==iEta_jet && iPhi_em==iPhi_jet) match_code=1;
				if(match_code==1) break;
			}
			if(match_code==1) break;
		}
	}

	TLorentzVector vec=*myjet;
	double dR=mypho->DeltaR(vec);
	dR_fj=dR;
	dPhi_fj=mypho->DeltaPhi(vec);
	dEta_fj=mypho->Eta()-myjet->Eta();
	return match_code;
}

void JetFilterModuleV3::ClearModuleOutput()       // clears the Module output from
{                              // previous event.

	//-------- my new global output parameters
	ClearJetStuff(jet04stuff); // clears JetClu-0.4 stuff
	ClearCommonStuff(allstuff); // clears common stuff

	return;
}

void JetFilterModuleV3::ClearJetStuff(JetStuff &stuff)
{ 
	// clears the JetStuff structures
	stuff.myMETcorr_th5.Set(-1.0E6,-1.0E6);
	stuff.mySumEtCorr_th5=-1.0E6;
	stuff.mySumEtJet_th5=-1.0E6;
	stuff.myMETcorr_th10.Set(-1.0E6,-1.0E6);
	stuff.mySumEtCorr_th10=-1.0E6;
	stuff.mySumEtJet_th10=-1.0E6;
	stuff.myMETcorr_th15.Set(-1.0E6,-1.0E6);
	stuff.mySumEtCorr_th15=-1.0E6;
	stuff.mySumEtJet_th15=-1.0E6;
	stuff.myMETcorr_th20.Set(-1.0E6,-1.0E6);
	stuff.mySumEtCorr_th20=-1.0E6;
	stuff.mySumEtJet_th20=-1.0E6;
	stuff.myNjet=0;
	stuff.myNjet_th5=0;
	stuff.myNjet_th10=0;
	stuff.myNjet_th15=0;
	stuff.myNjet_th20=0;
	stuff.myNjet_th25=0;
	stuff.myNjet_th30=0;
	stuff.myNjet_th35=0;
	stuff.smeared_Njet_th0=0;
	stuff.smeared_Njet_th5=0;
	stuff.smeared_Njet_th10=0;
	stuff.smeared_Njet_th15=0;
	stuff.smeared_Njet_th20=0;
	stuff.smeared_Njet_th25=0;
	stuff.smeared_Njet_th30=0;
	stuff.smeared_Njet_th35=0;
	stuff.Jet_raw.clear();
	stuff.Jet_lev1.clear();
	stuff.Jet_lev4.clear();
	stuff.Jet_lev5.clear();
	stuff.Jet_lev6.clear();
	stuff.Jet_lev7.clear();
	stuff.Jet_raw_noEMobj.clear();
	stuff.Jet_lev1_noEMobj.clear();
	stuff.Jet_lev4_noEMobj.clear();
	stuff.Jet_lev5_noEMobj.clear();
	stuff.Jet_lev6_noEMobj.clear();
	stuff.Jet_lev7_noEMobj.clear();
	stuff.newJetLev6.clear();
	stuff.smeared_newJetLev6_noEMObj.clear();
	stuff.bvMetAdded.clear();
	stuff.iL6MatchedHadJetIndex.clear();
	stuff.fL6NoEm_MisMeasureProb.clear();
	stuff.JetNtrk.clear();
	stuff.JetNtwr.clear();
	stuff.EtaDet.clear();
	stuff.EtaDetCorr.clear();
	stuff.EmFrRaw.clear();
	stuff.EmFrCorr.clear();
	stuff.Nobj_match.clear();
	stuff.Npho_match.clear();
	stuff.Nele_match.clear();
	stuff.Nmu_match.clear();
	stuff.Ntau_match.clear();
	stuff.Nbtag_match.clear();
	stuff.JetBlockInd.clear();
	stuff.smear_factor.clear();
	return;
}
//_________________________ clears parameters which do not depend on Jet Cone size
void JetFilterModuleV3::ClearCommonStuff(CommonStuff &stuff)
{ 
	// clears CommonStuff structures
	stuff.myRawPhoton.clear();
	stuff.myCorrPhoton.clear();
	stuff.myPhoInd.clear();
	stuff.tightId.clear();
	stuff.myPhoEmFr.clear();
	stuff.myPhoEtaDet.clear();
	stuff.myPhoXces.clear();
	stuff.myPhoZces.clear();
	stuff.myRawElectron.clear();
	stuff.myEleInd.clear();
	stuff.myCorrElectron.clear();
	stuff.myEleEmFr.clear();
	stuff.myEleEtaDet.clear();
	stuff.myElePhiDet.clear();
	stuff.myElePprEt.clear();
	stuff.myEleXces.clear();
	stuff.myEleZces.clear();
	stuff.myRawMuon.clear();
	stuff.myCorrMuon.clear();
	stuff.myRawTau.clear();
	stuff.myCorrTau.clear();
	stuff.myRawBjet.clear();
	stuff.myCorrBjet.clear();
	stuff.myMET_raw.Set(-1.0E6,-1.0E6);
	stuff.myMET0_raw.Set(-1.0E6,-1.0E6);
	stuff.mySumEt_raw=-1.0E6;
	stuff.newVertexZ=-1000.0;

	return;
}

void JetFilterModuleV3::ClearMatchStuff() {
	matchstuff.clear();
	return;
}


//__________________________________________________________ accessors to the
//__________________________________________________________ output params
//_________ need to specify which jets to return (variable "cone"): 0--cone 0.4; 1--cone 0.7; 2--cone 1.0

TLorentzVector* JetFilterModuleV3::GetMyJet_raw(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_raw.size()) return &_stuff.Jet_raw.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJet_lev1(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev1.size()) return &_stuff.Jet_lev1.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJet_lev4(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev4.size()) return &_stuff.Jet_lev4.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJet_lev5(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev5.size()) return &_stuff.Jet_lev5.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJet_lev6(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev6.size()) return &_stuff.Jet_lev6.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJetNoEMobj_raw(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_raw_noEMobj.size()) return &_stuff.Jet_raw_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJetNoEMobj_lev1(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev1_noEMobj.size()) return &_stuff.Jet_lev1_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJetNoEMobj_lev4(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev4_noEMobj.size()) return &_stuff.Jet_lev4_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJetNoEMobj_lev5(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev5_noEMobj.size()) return &_stuff.Jet_lev5_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJetNoEMobj_lev6(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev6_noEMobj.size()) return &_stuff.Jet_lev6_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV3::GetMyJetNoEMobj_lev7(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev7_noEMobj.size()) return &_stuff.Jet_lev7_noEMobj.at(i);
	else return NULL;
}

double JetFilterModuleV3::GetMyJetEtaDet(int cone, int i) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EtaDet.size()) return _stuff.EtaDet.at(i);
	else return -1.0E6;
}
double JetFilterModuleV3::GetMyJetEtaDetCorr(int cone, int i) { // after removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EtaDetCorr.size()) return _stuff.EtaDetCorr.at(i);
	else return -1.0E6;
}
double JetFilterModuleV3::GetMyJetEmFrRaw(int cone, int i) {  // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EmFrRaw.size()) return _stuff.EmFrRaw.at(i);
	else return -1.0E6;
}
double JetFilterModuleV3::GetMyJetEmFrCorr(int cone, int i) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EmFrCorr.size()) return _stuff.EmFrCorr.at(i);
	else return -1.0E6;
}
int JetFilterModuleV3::GetMyJetNtrk(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetNtrk.size()) return _stuff.JetNtrk.at(i);
	else return -1;
}
int JetFilterModuleV3::GetMyJetNtwr(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetNtwr.size()) return _stuff.JetNtwr.at(i);
	else return -1;
}
int JetFilterModuleV3::GetMyJetNobjMatch(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nobj_match.size()) return _stuff.Nobj_match.at(i);
	else return -1;
}
int JetFilterModuleV3::GetMyJetNphoMatch(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Npho_match.size()) return _stuff.Npho_match.at(i);
	else return -1;
}
int JetFilterModuleV3::GetMyJetNeleMatch(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nele_match.size()) return _stuff.Nele_match.at(i);
	else return -1;
}
int JetFilterModuleV3::GetMyJetBlockInd(int cone, int i)
{ 
	// original jet index in JetBlock, need this after jet reordering
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetBlockInd.size()) return _stuff.JetBlockInd.at(i);
	else return -1;
}
int JetFilterModuleV3::GetMyNjet(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25,
	// 5--et>30, 5--et>35
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return _stuff.myNjet;
	if(threshold==1) return _stuff.myNjet_th5;
	if(threshold==2) return _stuff.myNjet_th10;
	if(threshold==3) return _stuff.myNjet_th15;
	if(threshold==4) return _stuff.myNjet_th20;
	if(threshold==5) return _stuff.myNjet_th25;
	if(threshold==6) return _stuff.myNjet_th30;
	if(threshold==7) return _stuff.myNjet_th35;
	//if(threshold<0 || threshold>7) return -1;
	return -1;
}
double JetFilterModuleV3::GetMySumEtCorr(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.mySumEt_raw; // for consistency
	if(threshold==1) return _stuff.mySumEtCorr_th5;
	if(threshold==2) return _stuff.mySumEtCorr_th10;
	if(threshold==3) return _stuff.mySumEtCorr_th15;
	if(threshold==4) return _stuff.mySumEtCorr_th20;
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}
double JetFilterModuleV3::GetMyHtCorr(int cone, int threshold, const int iJetCorrLev)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	// CorrLev decides the lev of corrected jet to be used to calcualte Ht.
	// default would be Lev6 as used through out MEtModel.
	// I may need use lev7 corrected jets(BE CAREFULL! need to change how SumEt 
	// is calculated too). Besides I had to do this chage since
	// I am using the MyHtAll() for both DATA/BG calcualtions and in each time 
	// a different set of jets will be passed in. Sep7,2009-Sam
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;

	double Ht= -1.0E6;
	if (iJetCorrLev==6) Ht = MyHtAll(_stuff.Jet_lev6_noEMobj, _stuff,allstuff,threshold);
	else if (iJetCorrLev==7) 
	{
		//Ht = MyHtAll(_stuff.Jet_lev7_noEMobj, _stuff,allstuff,threshold);
		std::cout << __FILE__ << ":" << __LINE__ << 
			":: Not implemented yet. You need chage MEt/SumEt calculations to use Lev7 jets too!" 
			<< std::endl;
		exit (1);
	} else
	{
		StdOut(__FILE__,__LINE__,3,
		" Not yet defined level of corrected jets requested in GetMyHtCorr()!");
		exit (1);
	}

	return Ht;
}

double JetFilterModuleV3::GetMyMetCorr(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.myMET_raw.Mod(); // for consistency
	if(threshold==1) return _stuff.myMETcorr_th5.Mod();
	if(threshold==2) return _stuff.myMETcorr_th10.Mod();
	if(threshold==3) return _stuff.myMETcorr_th15.Mod();
	if(threshold==4) return _stuff.myMETcorr_th20.Mod();
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}
double JetFilterModuleV3::GetMyMetXCorr(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.myMET_raw.Px(); // for consistency
	if(threshold==1) return _stuff.myMETcorr_th5.Px();
	if(threshold==2) return _stuff.myMETcorr_th10.Px();
	if(threshold==3) return _stuff.myMETcorr_th15.Px();
	if(threshold==4) return _stuff.myMETcorr_th20.Px();
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}
//___________________________________________________________________
double JetFilterModuleV3::GetMyMetYCorr(int cone, int threshold) 
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.myMET_raw.Py(); // for consistency
	if(threshold==1) return _stuff.myMETcorr_th5.Py();
	if(threshold==2) return _stuff.myMETcorr_th10.Py();
	if(threshold==3) return _stuff.myMETcorr_th15.Py();
	if(threshold==4) return _stuff.myMETcorr_th20.Py();
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}
//___________________________________________________________________
double JetFilterModuleV3::GetMyMetPhiCorr(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.myMET_raw.Phi(); // for consistency
	if(threshold==1) return _stuff.myMETcorr_th5.Phi();
	if(threshold==2) return _stuff.myMETcorr_th10.Phi();
	if(threshold==3) return _stuff.myMETcorr_th15.Phi();
	if(threshold==4) return _stuff.myMETcorr_th20.Phi();
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}


//__________________________________________________________ setters of the
//__________________________________________________________ output params

//_________ need to specify which jets to return (variable "cone"): 
// 0--cone 0.4; 1--cone 0.7; 2--cone 1.0
void JetFilterModuleV3::SetMyJet_raw(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_raw.size()) _stuff.Jet_raw.at(i)=vec;
	else _stuff.Jet_raw.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJet_lev1(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev1.size()) _stuff.Jet_lev1.at(i)=vec;
	else _stuff.Jet_lev1.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJet_lev4(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev4.size()) _stuff.Jet_lev4.at(i)=vec;
	else _stuff.Jet_lev4.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJet_lev5(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev5.size()) _stuff.Jet_lev5.at(i)=vec;
	else _stuff.Jet_lev5.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJet_lev6(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev6.size()) _stuff.Jet_lev6.at(i)=vec;
	else _stuff.Jet_lev6.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNoEMobj_raw(int cone, int i, TLorentzVector vec)  {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_raw_noEMobj.size()) _stuff.Jet_raw_noEMobj.at(i)=vec;
	else _stuff.Jet_raw_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNoEMobj_lev1(int cone, int i, TLorentzVector vec)  {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev1_noEMobj.size()) _stuff.Jet_lev1_noEMobj.at(i)=vec;
	else _stuff.Jet_lev1_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNoEMobj_lev4(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev4_noEMobj.size()) _stuff.Jet_lev4_noEMobj.at(i)=vec;
	else _stuff.Jet_lev4_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNoEMobj_lev5(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev5_noEMobj.size()) _stuff.Jet_lev5_noEMobj.at(i)=vec;
	else _stuff.Jet_lev5_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNoEMobj_lev6(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev6_noEMobj.size()) _stuff.Jet_lev6_noEMobj.at(i)=vec;
	else _stuff.Jet_lev6_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetEtaDet(int cone, int i, double param) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EtaDet.size()) _stuff.EtaDet.at(i)=param;
	else _stuff.EtaDet.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetEtaDetCorr(int cone, int i, double param) { // after removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EtaDetCorr.size()) _stuff.EtaDetCorr.at(i)=param;
	else _stuff.EtaDetCorr.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetEmFrRaw(int cone, int i, double param) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EmFrRaw.size()) _stuff.EmFrRaw.at(i)=param;
	else _stuff.EmFrRaw.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetEmFrCorr(int cone, int i, double param) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EmFrCorr.size()) _stuff.EmFrCorr.at(i)=param;
	else _stuff.EmFrCorr.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNtrk(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetNtrk.size()) _stuff.JetNtrk.at(i)=param;
	else _stuff.JetNtrk.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNtwr(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetNtwr.size()) _stuff.JetNtwr.at(i)=param;
	else _stuff.JetNtwr.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNobjMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nobj_match.size()) _stuff.Nobj_match.at(i)=param;
	else _stuff.Nobj_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNphoMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Npho_match.size()) _stuff.Npho_match.at(i)=param;
	else _stuff.Npho_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNeleMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nele_match.size()) _stuff.Nele_match.at(i)=param;
	else _stuff.Nele_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNmuMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nmu_match.size()) _stuff.Nmu_match.at(i)=param;
	else _stuff.Nmu_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNtauMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Ntau_match.size()) _stuff.Ntau_match.at(i)=param;
	else _stuff.Ntau_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetNbtagMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nbtag_match.size()) _stuff.Nbtag_match.at(i)=param;
	else _stuff.Nbtag_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyJetBlockInd(int cone, int i, int param) {
	// original jet index in JetBlock, need this after jet reordering
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetBlockInd.size()) _stuff.JetBlockInd.at(i)=param;
	else _stuff.JetBlockInd.push_back(param);
	return;
}

//___________________________________________________________________
//this is not used???28-12-2008,sam
void JetFilterModuleV3::SetMyNjet(int cone, int threshold, int param) {
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25, 6--et>30, 7--et>35
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) _stuff.myNjet=param;
	if(threshold==1) _stuff.myNjet_th5=param;
	if(threshold==2) _stuff.myNjet_th10=param;
	if(threshold==3) _stuff.myNjet_th15=param;
	if(threshold==4) _stuff.myNjet_th20=param;
	if(threshold==5) _stuff.myNjet_th25=param;
	if(threshold==6) _stuff.myNjet_th30=param;
	if(threshold==7) _stuff.myNjet_th35=param;
	if(threshold<0 || threshold>7) return;
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMySumEtCorr(int cone, int threshold, double param) {
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) allstuff.mySumEt_raw=param; // for consistency
	if(threshold==1) _stuff.mySumEtCorr_th5=param;
	if(threshold==2) _stuff.mySumEtCorr_th10=param;
	if(threshold==3) _stuff.mySumEtCorr_th15=param;
	if(threshold==4) _stuff.mySumEtCorr_th20=param;
	if(threshold<0 || threshold>4) return;
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyMetCorr(int cone, int threshold, double px, double py) {
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) allstuff.myMET_raw.Set(px,py); // for consistency
	if(threshold==1) _stuff.myMETcorr_th5.Set(px,py);
	if(threshold==2) _stuff.myMETcorr_th10.Set(px,py);
	if(threshold==3) _stuff.myMETcorr_th15.Set(px,py);
	if(threshold==4) _stuff.myMETcorr_th20.Set(px,py);
	if(threshold<0 || threshold>4) return;
	return;
}
//___________________________________________________________________
void JetFilterModuleV3::SetMyMetCorr(int cone, int threshold, const TVector2& vec) {
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) allstuff.myMET_raw.Set(vec); // for consistency
	if(threshold==1) _stuff.myMETcorr_th5.Set(vec);
	if(threshold==2) _stuff.myMETcorr_th10.Set(vec);
	if(threshold==3) _stuff.myMETcorr_th15.Set(vec);
	if(threshold==4) _stuff.myMETcorr_th20.Set(vec);
	if(threshold<0 || threshold>4) return;
	return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   My Cuts are defined here
//___________________________________________________________________
int JetFilterModuleV3::MyNjetCut(JetStuff jetstuff)
{
/*uncomment for debugging - sam 03-26-2010
 * { std::cout << __FUNCTION__ << "\n"
			<< " jetstuff.myNjet      = " << jetstuff.myNjet << "\n"
			<< " jetstuff.myNjet_th5  = " << jetstuff.myNjet_th5 << "\n"
			<< " jetstuff.myNjet_th10 = " << jetstuff.myNjet_th10 << "\n"
			<< " jetstuff.myNjet_th15 = " << jetstuff.myNjet_th15 << "\n"
			<< " jetstuff.myNjet_th20 = " << jetstuff.myNjet_th20 << std::endl;
	}
	if(jetstuff.myNjet<MinNjet || jetstuff.myNjet>MaxNjet) { std::cout << " failed MinNjet " <<std::endl;  return 0;}
	if(jetstuff.myNjet_th5<MinNjet5 || jetstuff.myNjet_th5>MaxNjet5) { std::cout << " failed MinNjet " <<std::endl; return 0;}
	if(jetstuff.myNjet_th10<MinNjet10 || jetstuff.myNjet_th10>MaxNjet10) { std::cout << " failed MinNjet " <<std::endl; return 0;}
	if(jetstuff.myNjet_th15<MinNjet15 || jetstuff.myNjet_th15>MaxNjet15) { std::cout << " failed MinNjet " <<std::endl; return 0;}
	if(jetstuff.myNjet_th20<MinNjet20 || jetstuff.myNjet_th20>MaxNjet20) { std::cout << " failed MinNjet " <<std::endl; return 0;}
	if(jetstuff.myNjet_th25<MinNjet25 || jetstuff.myNjet_th25>MaxNjet25) { std::cout << " failed MinNjet " <<std::endl; return 0;}
	if(jetstuff.myNjet_th30<MinNjet30 || jetstuff.myNjet_th30>MaxNjet30) { std::cout << " failed MinNjet " <<std::endl; return 0;}
	if(jetstuff.myNjet_th35<MinNjet35 || jetstuff.myNjet_th35>MaxNjet35) { std::cout << " failed MinNjet " <<std::endl; return 0;}
*/
	if(jetstuff.myNjet<MinNjet || jetstuff.myNjet>MaxNjet) return 0;
	if(jetstuff.myNjet_th5<MinNjet5 || jetstuff.myNjet_th5>MaxNjet5) return 0;
	if(jetstuff.myNjet_th10<MinNjet10 || jetstuff.myNjet_th10>MaxNjet10) return 0;
	if(jetstuff.myNjet_th15<MinNjet15 || jetstuff.myNjet_th15>MaxNjet15) return 0;
	if(jetstuff.myNjet_th20<MinNjet20 || jetstuff.myNjet_th20>MaxNjet20) return 0;
	if(jetstuff.myNjet_th25<MinNjet25 || jetstuff.myNjet_th25>MaxNjet25) return 0;
	if(jetstuff.myNjet_th30<MinNjet30 || jetstuff.myNjet_th30>MaxNjet30) return 0;
	if(jetstuff.myNjet_th35<MinNjet35 || jetstuff.myNjet_th35>MaxNjet35) return 0;

	
	return 1;

}

//___________________________________________________________________
// sets the number of jets. We must cut on jets
// after we smear them. 12-28-2008, sam
//___________________________________________________________________
void JetFilterModuleV3::GetNjets(
			const std::vector<TLorentzVector>& vJet,	std::vector<int>& Njets)
{

	int Njet_th0=0, Njet_th5=0, Njet_th10=0, Njet_th15=0, Njet_th20=0, 
		 Njet_th25=0, Njet_th30=0, Njet_th35=0;

	for (unsigned int i=0; i < vJet.size(); ++i)
	{
		float fJetPt = vJet.at(i).Pt();
		if (fJetPt>0.0)  Njet_th0++;
		if (fJetPt>5.0)  Njet_th5++;
		if (fJetPt>10.0) Njet_th10++;
		if (fJetPt>15.0) Njet_th15++;
		if (fJetPt>20.0) Njet_th20++;
		if (fJetPt>25.0) Njet_th25++;
		if (fJetPt>30.0) Njet_th30++;
		if (fJetPt>35.0) Njet_th35++;
	}
	Njets.push_back(Njet_th0);
	Njets.push_back(Njet_th5);
	Njets.push_back(Njet_th10);
	Njets.push_back(Njet_th15);
	Njets.push_back(Njet_th20);
	Njets.push_back(Njet_th25);
	Njets.push_back(Njet_th30);
	Njets.push_back(Njet_th35);

	if (debug)
	{
		std::cout << __FUNCTION__ << "::" << __LINE__ 
			<< "::GetNjets njets 0,5,10,15,20,25,30,35=" << Njet_th0 << "," 
			<< Njet_th5 << "," << Njet_th10 << "," << Njet_th15 << "," << Njet_th20 << "," 
			<< Njet_th25 << "," << Njet_th30 << ", " << Njet_th35 << std::endl;
	}
}
//___________________________________________________________________
int JetFilterModuleV3::GetNjets(
			const std::vector<TLorentzVector>& vJet,	const float fEt)
{
// for a given set of jets and Et threshold
	int iNjets = 0; 

	for (unsigned int i=0; i < vJet.size(); ++i)
	{
		if (vJet.at(i).Pt()>=fEt) ++iNjets;
	}

	return iNjets;
}

//--------------------------------------------------------------------------
int JetFilterModuleV3::MyGenericCut(JetStuff jetstuff, double dz) {
	//__________________ cut on dZ=Zvx-Zjet
	if(fabs(dz)<MindZ || fabs(dz)>MaxdZ) return 0;
	//__________________ cut on 1st jet Et
	if(jetstuff.myNjet>0)
	{
		if(fabs(jetstuff.Jet_lev6_noEMobj[0].Eta())>MinJetEta
			&& fabs(jetstuff.Jet_lev6_noEMobj[0].Eta())<MaxJetEta
			&& (jetstuff.Jet_lev6_noEMobj[0].Pt()<MinEt1st 
			    || jetstuff.Jet_lev6_noEMobj[0].Pt()>MaxEt1st))
		{
			return 0;
		}
	}
	//__________________ cut on 2nd jet Et
	if(jetstuff.myNjet>1)
	{
		if(fabs(jetstuff.Jet_lev6_noEMobj[1].Eta())>MinJetEta
			&& fabs(jetstuff.Jet_lev6_noEMobj[1].Eta())<MaxJetEta
			&& (jetstuff.Jet_lev6_noEMobj[1].Pt()<MinEt2nd 
			    || jetstuff.Jet_lev6_noEMobj[1].Pt()>MaxEt2nd))
		{
			return 0;
		}
	}
	return 1;
}

//check to make sure that there is no overlap of photons and electrons in the event
//--------------------------------------------------------------------------
int JetFilterModuleV3::MyDuplicateCut(CommonStuff miscstuff) {

	double dRMin=0.2;
	int passcode=0;
	if(fRemoveDuplicate!=0)
	{
		if(bad_EMjet_match_flag==0) return 1;
		for(unsigned int j=0; j<miscstuff.myCorrPhoton.size(); j++)
		{
			for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
			{
				if(fabs(miscstuff.myEleEtaDet.at(i))<1.2)
				{
					double dPhi=fabs(TVector2::Phi_mpi_pi(miscstuff.myElePhiDet.at(i)-miscstuff.myCorrPhoton[j].Phi()));
					double dEta=fabs(miscstuff.myEleEtaDet.at(i)-miscstuff.myPhoEtaDet[j]);
					double dR=sqrt(dPhi*dPhi+dEta*dEta);
					if(dR<dRMin) return 1;
				}
			}
		}
	}
	return passcode;
}



//----- Test version as of 04/04/07

//________________________________ returns SumEt=jets+photons+leptons+taus+met (jets above Eta & Et thresholds)
//                                 metscenario= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
double JetFilterModuleV3::MyHtAll(const std::vector<TLorentzVector> vJets, 
									JetStuff jetstuff, CommonStuff miscstuff, 
									int metscenario)
{
	double ht=0.0;
	//_________________ contribution from jets
	for(int i=0; i<jetstuff.myNjet; i++)
	{
		if(fabs(vJets.at(i).Eta())>MinJetEta
			&& fabs(vJets.at(i).Eta())<MaxJetEta
			&& fabs(vJets.at(i).Pt())>MinEtThr
			&& fabs(vJets.at(i).Pt())<MaxEtThr)
		{
			ht=ht+vJets.at(i).Pt();
		}
	}
	//_________________ contribution from photons
	for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
	{
		ht=ht+miscstuff.myCorrPhoton.at(i).Pt();
	}
	//_________________ contribution from electrons
	for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
	{
		ht=ht+miscstuff.myCorrElectron.at(i).Pt();
	}
	//_________________ contribution from muons
	for(unsigned int i=0; i<miscstuff.myCorrMuon.size(); i++)
	{
		ht=ht+miscstuff.myCorrMuon.at(i).Pt();
	}
	//_________________ contribution from taus
	for(unsigned int i=0; i<miscstuff.myCorrTau.size(); i++)
	{
		ht=ht+miscstuff.myCorrTau.at(i).Pt();
	}
	//_________________ contribution from MET
	if(metscenario==0) ht=ht+miscstuff.myMET_raw.Mod();
	if(metscenario==1) ht=ht+jetstuff.myMETcorr_th5.Mod();
	if(metscenario==2) ht=ht+jetstuff.myMETcorr_th10.Mod();
	if(metscenario==3) ht=ht+jetstuff.myMETcorr_th15.Mod();
	if(metscenario==4) ht=ht+jetstuff.myMETcorr_th20.Mod();
	if(metscenario<0 || metscenario>4) ht=ht+miscstuff.myMET_raw.Mod();
	return ht;
}

//_________________________________________ dump jet info for debugging
void JetFilterModuleV3::DumpJets(const std::string func, const unsigned int line,
								const JetStuff &jetstuff, const int lev)
{
	// lev: 0=raw, 1=lev1, 4=lev4, 5=lev5, 6=lev6, 7=lev7
	// lev: 10=raw noEM, 11=lev1 noEM, 14=lev4 noEM, 15=lev5 noEM, 16=lev6 noEM, 17=lev7 noEM
	// lev 20 = newJetLev6, 21 = smeared_newJetLev6_noEMObj
	
	if (jetstuff.myNjet>0)
	{
		std::string sLev;
		switch (lev)
		{
			case 0: { sLev = "RAW"; break;}
			case 1: { sLev = "LEV1"; break;}
			case 4: { sLev = "LEV4"; break;}
			case 5: { sLev = "LEV5"; break;}
			case 6: { sLev = "LEV6"; break;}
			case 10: { sLev = "RAW NO EM"; break;}
			case 11: { sLev = "LEV1 NO EM"; break;}
			case 14: { sLev = "LEV4 NO EM"; break;}
			case 15: { sLev = "LEV5 NO EM"; break;}
			case 16: { sLev = "LEV6 NO EM"; break;}
			case 20: { sLev = "newJetL6 NO EM"; break;}
			case 21: { sLev = "smeared_newJetL6 NO EM"; break;}
			default: { sLev = "UNKOWN!"; }
		};
		
		std::cout << "=== Caller,line::" <<func << ":" << line <<": Jets For of " << lev  << " - " << sLev << std::endl;
		std::cout << std::setw(3) << "ind" <<  std::setw(10) << "jetPt" <<  std::setw(10) << "eta"  <<  std::setw(6) << "JBind" << std::endl;

		for (int i=0; i<jetstuff.myNjet; i++)
		{
			float jetEt = 0.0;
			float eta = -99.99;
			int index = jetstuff.JetBlockInd.at(i);
			switch (lev)
			{
				case 6:
					{
						jetEt = jetstuff.Jet_lev6.at(i).Pt();
						eta   = jetstuff.Jet_lev6.at(i).Eta();

						break;
					}
				
				case 16:
					{
						jetEt = jetstuff.Jet_lev6_noEMobj.at(i).Pt();
						eta   = jetstuff.Jet_lev6_noEMobj.at(i).Eta();
						break;
					}
				case 20:
					{
						jetEt = jetstuff.newJetLev6.at(i).Pt();
						eta   = jetstuff.newJetLev6.at(i).Eta();
						break;
					}
				case 21:
					{
						jetEt = jetstuff.smeared_newJetLev6_noEMObj.at(i).Pt();
						eta   = jetstuff.smeared_newJetLev6_noEMObj.at(i).Eta();
						break;
					}
				default:
					{
						std::cout << "function not defined for lev= " << lev << " jets!\n";
						return;
					}

			}
			std::cout << std::setw(3) << i <<  std::setw(10) << jetEt <<  std::setw(10) << eta  <<  std::setw(6) << index << std::endl;
		}
		
	} else {
		std::cout << "=== Caller,line::" <<func << ":" << line <<": NO Jets For lev " << lev << std::endl;
	}

}
//________________________________________________booking pho-jet match histo
void JetFilterModuleV3::BookMatchStudyHistograms(MatchStudyHisto_t& Hist, const char* Folder, const char* algoname) {

	char name [200];
	char title[200];

	sprintf(name,"MatchNtwr");
	sprintf(title,"%s: Number of towers in matched jet",algoname);
	Hist.fMatchNtwr=new TH1F(name,title,50,-0.5,49.5);
	AddHistogram(Hist.fMatchNtwr,Folder);
	sprintf(name,"MatchDelR");
	sprintf(title,"%s: #Delta R between jet and matched EM object",algoname);
	Hist.fMatchDelR=new TH1F(name,title,500,0.0,10.0);
	AddHistogram(Hist.fMatchDelR,Folder);
	sprintf(name,"MatchDelPhi");
	sprintf(title,"%s: #Delta#phi between jet and matched EM object",algoname);
	Hist.fMatchDelPhi=new TH1F(name,title,200,0.0,4.0);
	AddHistogram(Hist.fMatchDelPhi,Folder);
	sprintf(name,"MatchDelEta");
	sprintf(title,"%s: #Delta#eta between jet and matched EM object",algoname);
	Hist.fMatchDelEta=new TH1F(name,title,200,0.0,4.0);
	AddHistogram(Hist.fMatchDelEta,Folder);
	sprintf(name,"MatchDelEtaDet");
	sprintf(title,"%s: #Delta#eta_{det} between jet and matched EM object",algoname);
	Hist.fMatchDelEtaDet=new TH1F(name,title,200,0.0,4.0);
	AddHistogram(Hist.fMatchDelEtaDet,Folder);
	sprintf(name,"MatchNmatch");
	sprintf(title," %s%s",algoname,": Number of matched EM objects");
	Hist.fMatchNmatch=new TH1F(name,title,10,-0.5,9.5);
	AddHistogram(Hist.fMatchNmatch,Folder);
	sprintf(name,"MatchEt_raw_b");
	sprintf(title," %s%s",algoname,": raw E_{T} of jet before matching");
	Hist.fMatchEt_raw_b=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fMatchEt_raw_b,Folder);
	sprintf(name,"MatchEt_raw_a");
	sprintf(title," %s%s",algoname,": raw E_{T} of jet after matching");
	Hist.fMatchEt_raw_a=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fMatchEt_raw_a,Folder);
	sprintf(name,"MatchEt_lev6_b");
	sprintf(title," %s%s",algoname,": lev6 E_{T} of jet before matching");
	Hist.fMatchEt_lev6_b=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fMatchEt_lev6_b,Folder);
	sprintf(name,"MatchEt_lev6_a");
	sprintf(title," %s%s",algoname,": lev6 E_{T} of jet after matching");
	Hist.fMatchEt_lev6_a=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fMatchEt_lev6_a,Folder);
	sprintf(name,"MatchEtJet2EtPho_b");
	sprintf(title," %s%s",algoname,": raw E^{jet}_{T}/E^{EMobj}_{T} before removing EM object");
	Hist.fMatchEtJet2EtPho_b=new TH1F(name,title,1000,0.0,10.0);
	AddHistogram(Hist.fMatchEtJet2EtPho_b,Folder);
	sprintf(name,"MatchEtJet2EtPho_a");
	sprintf(title," %s%s",algoname,": raw E^{jet}_{T}/E^{EMobj}_{T} after removing EM object");
	Hist.fMatchEtJet2EtPho_a=new TH1F(name,title,1000,0.0,10.0);
	AddHistogram(Hist.fMatchEtJet2EtPho_a,Folder);
	sprintf(name,"MatchDelRoldnew");
	sprintf(title,"%s: #Delta R between new and old jet",algoname);
	Hist.fMatchDelRoldnew=new TH1F(name,title,500,0.0,10.0);
	AddHistogram(Hist.fMatchDelRoldnew,Folder);
	sprintf(name,"MatchDelEtaDetoldnew");
	sprintf(title,"%s: #Delta#eta_{det} between new and old jet",algoname);
	Hist.fMatchDelEtaDetoldnew=new TH1F(name,title,200,-4.0,4.0);
	AddHistogram(Hist.fMatchDelEtaDetoldnew,Folder);

	return;
}

//_________________________________________ dump EM objects info for debugging
// prints all EM Objects by default
//------------------------------------------------------------------------------------
void JetFilterModuleV3::DumpEMobjects(const std::string func, const unsigned int line,
								const CommonStuff& commonStuff, const int type)
{
	//type 0=all (pho+ele+++), 1 = pho only, 2= ele only
	if (type == 0 || type == 1)
	{
		if (commonStuff.myCorrPhoton.size())
		{
			std::cout << __FUNCTION__ << " Caller,line::" <<func << ":" << line <<": Photons" << std::endl;
			std::cout << setw(10) << "Index" << setw(10) << "Pt" << setw(10) << "Eta" << setw(10) << "Phi" << std::endl;
			for (unsigned int i=0; i<commonStuff.myCorrPhoton.size(); ++i)
			{
				std::cout << setw(10) << i 
							<< setw(10) << commonStuff.myCorrPhoton.at(i).Pt() 
							<< setw(10) << commonStuff.myPhoEtaDet.at(i) 
							<< setw(10) << commonStuff.myCorrPhoton.at(i).Phi() 
							<< std::endl;
			}
		}
	}
	
	if (type == 0 || type == 2)
	{
		if (commonStuff.myCorrElectron.size())
		{
			std::cout << "=== Caller,line::" <<func << ":" << line <<": Electrons" << std::endl;
			std::cout << setw(10) << "Index" << setw(10) << "Pt" << setw(10) << "Eta" << setw(10) << "Phi" << std::endl;
			for (unsigned int i=0; i<commonStuff.myCorrElectron.size(); i++)
			{
				std::cout << setw(10) << i 
					<< setw(10) << commonStuff.myCorrElectron.at(i).Pt() 
					<< setw(10) << commonStuff.myEleEtaDet.at(i) 
					<< setw(10) << commonStuff.myCorrElectron.at(i).Phi() 
					<< std::endl;
			}
		}
	}
	
}

//---------------------------------------------------------------------------
// Returns original jet block index of the jet matched to an EM object
//---------------------------------------------------------------------------
int JetFilterModuleV3::GetMatch_JetIndex(int ObjType, int EmInd) const
{
	if (sam_matchstuff.size() >0) {
		for (unsigned int j=0; j < sam_matchstuff.size(); j++) {
			if (sam_matchstuff[j].EmObjType_match == ObjType) {
				if (sam_matchstuff[j].EmInd_match == EmInd ) {
					return sam_matchstuff[j].JetInd_match;
				}
			}
		}
		std::cerr << __FILE__ << "::" << __LINE__ << "::" << __FUNCTION__
			<< ":: ERROR::: EM obj did not match" <<std::endl;
		std::cout << "\tERROR::Run,Evt:: " << fHeaderBlock->RunNumber() << ", " << fHeaderBlock->EventNumber() << std::endl;
		//	exit (1);
	}

	std::cerr << __FILE__ << "::" << __LINE__ << "::" << __FUNCTION__
		<< ":: ERROR:: NO MATCHING INFO FOUND. Returining -1 ! " <<std::endl;
	std::cout << "\tERROR::Run,Evt:: " << fHeaderBlock->RunNumber() << ", " << fHeaderBlock->EventNumber() << std::endl;
	return -1;
}

