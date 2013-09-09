/*

This is a modified version of TMyJetFilterModule. It is to be used in
g+2j+X analysis.
  
*/
//_____________________________________________________________________________
#include <iomanip>
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <TRandom.h>
#include "Stntuple/loop/TStnAna.hh"

#include "samantha/PhoJetsMet/TMyJetFilterModule.hh"
#include "samantha/PhoJetsMet/TInit.hh" // Sam's module

//------- need this for jet-EM matching

/*bool SortTowersByEnergy(TCalTower* c1, TCalTower* c2){
  if(c1->Energy() > c2->Energy()) return true;
  else return false;
}
*/
//ClassImp(TMyJetFilterModule)

//_____________________________________________________________________________
// HERE WE DEFINE ALL THE CUTS
//_____________________________________________________________________________
TMyJetFilterModule::TMyJetFilterModule(const char* name, const char* title):
  TStnModule(name,title)
/*{{{*/
{
  fUnclParamSwitch=0; // 0=sideband (default); 1=Z->ee unclustered energy parametrizations

  Npoints=1; // number of random points to generate

//______________________________________ setting default cut values
  fUseVerbose=0;
  fDumpEvent=0;
  fJetAlgo=0;      // by default it's JetClu(0.4)
  fMyMetScenario=3; // by default met is corrected for jets with Et>15
  fUseMetPDF=0; // MetPDF scenarios: 0=no MetPDF; 
                // 1= 3 sigma (0.27%) & Npoints=370; 
                // 2= 3.29 sigma (0.1%) & Npoints=1000;
                // 3= 3.89 sigma (0.01%) & Npoints=10000;
                // 4= 4 sigma (0.0063%) & Npoints=15873;
                // 5= 2.58 sigma (1%), to calculate Mean & Sigma, Npoints=100.
  MetToy_max=200.0; // max value of toy Met 
  MetToy_min=0.0; // min value of toy Met
  fSelectSigMetEvent=0; // 0= do nothing; -1=reject event; 1=select event

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
  MinDeltaPhiJMet=-10.0; // every event should pass this cut
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
  fSelectMetEvent=0; // 0=no selection; 1=select events
  fRemoveDuplicate=0; // 0=no selection; 1=remove events; -1=select event
  fRemoveBadMet=0; // 0=do nothing; 1=remove event; -1=select event  
  //________________________
  //________________________
  fDatFileName="results/MyJetFilter_output.dat";
  fLargeMetEventFileName="results/largeMetEventList.dat";
  fDumpEventFileName="results/dumpEvent.dat";
  //________________________
  fJTC_coneSize=0;  // JetClue R=0.4
  fJTC_version=5;   // version of correction, suggested by Anwar for 5.3.1, may change later
  fJTC_level=7;     // full correction
  fJTC_systcode=0;  // default correction
  fJTC_imode=1;     // DATA is default
  std::cout<<"Hi, entering TMyJetFilterModule"<<std::endl;
}
/*}}}*/


//_____________________________________________________________________________
TMyJetFilterModule::~TMyJetFilterModule() {
}


//_____________________________________________________________________________
// BOOK ALL YOUR HISTOGRAMS HERE
//_____________________________________________________________________________
void TMyJetFilterModule::BookHistograms()
/*{{{*/
{
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

  //_______________________ booking met study histograms
  //
  //----- Cone 0.4, Et>15
  sprintf(folder_name,"METJetClu0.4_scen3"); 
  fol = (TFolder*) hist_folder->FindObject(folder_name);
  if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
  BookMetStudyHistograms(fMetStudyJet04sc3,Form("Hist/%s",folder_name),folder_name);

  //------ booking matching histograms
  sprintf(folder_name,"MatchPhoJet0.4"); //----- Cone 0.4
  fol = (TFolder*) hist_folder->FindObject(folder_name);
  if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
  BookMatchStudyHistograms(fMatchPhoJet04,Form("Hist/%s",folder_name),folder_name);

//------- booking met cleanup histograms
  sprintf(folder_name,"MetCleanup0.4"); //----- Cone 0.4
  fol = (TFolder*) hist_folder->FindObject(folder_name);
  if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
  BookMetCleanupHistograms(fCleanup,Form("Hist/%s",folder_name),folder_name);

//------- booking generated met cleanup histograms
  sprintf(folder_name,"GenMetCleanup0.4"); //----- Cone 0.4, default MetModel
  fol = (TFolder*) hist_folder->FindObject(folder_name);
  if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
  BookMetCleanupHistograms(fGenCleanup,Form("Hist/%s",folder_name),folder_name);

  sprintf(folder_name,"GenMetDevmCleanup0.4"); //----- Cone 0.4, MetModel-SigmaJet
  fol = (TFolder*) hist_folder->FindObject(folder_name);
  if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
  BookMetCleanupHistograms(fGenDevmCleanup,Form("Hist/%s",folder_name),folder_name);

  sprintf(folder_name,"GenMetDevpCleanup0.4"); //----- Cone 0.4, MetModel+SigmaJet
  fol = (TFolder*) hist_folder->FindObject(folder_name);
  if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
  BookMetCleanupHistograms(fGenDevpCleanup,Form("Hist/%s",folder_name),folder_name);

  sprintf(folder_name,"GenMetDevmUnclCleanup0.4"); //----- Cone 0.4, MetModel-SigmaUncl
  fol = (TFolder*) hist_folder->FindObject(folder_name);
  if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
  BookMetCleanupHistograms(fGenDevmUnclCleanup,Form("Hist/%s",folder_name),folder_name);

  sprintf(folder_name,"GenMetDevpUnclCleanup0.4"); //----- Cone 0.4, MetModel+SigmaUncl
  fol = (TFolder*) hist_folder->FindObject(folder_name);
  if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
  BookMetCleanupHistograms(fGenDevpUnclCleanup,Form("Hist/%s",folder_name),folder_name);

  sprintf(folder_name,"PhoMetStudy0.4"); //----- Cone 0.4, correlation between Photon and MET
  fol = (TFolder*) hist_folder->FindObject(folder_name);
  if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
  BookPhoMetStudyHistograms(fPhoMetStudy,Form("Hist/%s",folder_name));

//------- booking Met probability histograms
  if(fUseMetPDF>0)
    {
      sprintf(folder_name,"MetProbability0.4"); //----- Cone 0.4
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      BookMetProbHistograms(fMetProbability,Form("Hist/%s",folder_name),folder_name);
    }

  return;
}
/*}}}*/

//_____________________________________________________________________________
//  CREATE  GENERAL JET HISTOGRAMS
//_____________________________________________________________________________
void TMyJetFilterModule::BookJetHistograms(JetGeneral_t& Hist,
										const char* Folder, const char* algoname)
/*{{{*/
{
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

  //------------------------------ Jet Resolution studies --------------------------------
  sprintf(name,"JetRes_MetSig_EtaDet_j5_b");
  sprintf(title,"#slash{E}_{T} significance, #slash{E}_{T}/#sigma_{E_{T}}, vs. #eta_{det} for E_{T}^{lev6}>5 GeV");
  Hist.fJetRes_MetSig_EtaDet_j5_b= new TProfile(name,title,104,0.0,2.6,0.0,20.0);  
  AddHistogram(Hist.fJetRes_MetSig_EtaDet_j5_b,Folder);
  sprintf(name,"JetRes_MetSig_EtaDet_j10_b");
  sprintf(title,"#slash{E}_{T} significance, #slash{E}_{T}/#sigma_{E_{T}}, vs. #eta_{det} for E_{T}^{lev6}>10 GeV");
  Hist.fJetRes_MetSig_EtaDet_j10_b= new TProfile(name,title,104,0.0,2.6,0.0,20.0);  
  AddHistogram(Hist.fJetRes_MetSig_EtaDet_j10_b,Folder);
  sprintf(name,"JetRes_MetSig_EtaDet_j15_b");
  sprintf(title,"#slash{E}_{T} significance, #slash{E}_{T}/#sigma_{E_{T}}, vs. #eta_{det} for E_{T}^{lev6}>15 GeV");
  Hist.fJetRes_MetSig_EtaDet_j15_b= new TProfile(name,title,104,0.0,2.6,0.0,20.0);  
  AddHistogram(Hist.fJetRes_MetSig_EtaDet_j15_b,Folder);
  sprintf(name,"JetRes_DetEta_j5_b");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>5 GeV");
  Hist.fJetRes_DetEta_j5_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_b,Folder);
  sprintf(name,"JetRes_DetEta_j10_b");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>10 GeV");
  Hist.fJetRes_DetEta_j10_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_b,Folder);
  sprintf(name,"JetRes_DetEta_j15_b");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>15 GeV");
  Hist.fJetRes_DetEta_j15_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_b,Folder);
  sprintf(name,"JetRes_DetEta_j5_Sig3dPhi04_b");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>5 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j5_Sig3dPhi04_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_Sig3dPhi04_b,Folder);
  sprintf(name,"JetRes_DetEta_j10_Sig3dPhi04_b");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>10 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j10_Sig3dPhi04_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_Sig3dPhi04_b,Folder);
  sprintf(name,"JetRes_DetEta_j15_Sig3dPhi04_b");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>15 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j15_Sig3dPhi04_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_Sig3dPhi04_b,Folder);
  sprintf(name,"JetRes_DetEta_j5_Sig3dPhi01_b");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>5 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j5_Sig3dPhi01_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_Sig3dPhi01_b,Folder);
  sprintf(name,"JetRes_DetEta_j10_Sig3dPhi01_b");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>10 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j10_Sig3dPhi01_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_Sig3dPhi01_b,Folder);
  sprintf(name,"JetRes_DetEta_j15_Sig3dPhi01_b");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>15 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j15_Sig3dPhi01_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_Sig3dPhi01_b,Folder);
  sprintf(name,"JetRes_DetEta_j5_FrSig3dPhi04_b");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>5 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j5_FrSig3dPhi04_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_FrSig3dPhi04_b,Folder);
  sprintf(name,"JetRes_DetEta_j10_FrSig3dPhi04_b");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>10 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j10_FrSig3dPhi04_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_FrSig3dPhi04_b,Folder);
  sprintf(name,"JetRes_DetEta_j15_FrSig3dPhi04_b");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>15 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j15_FrSig3dPhi04_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_FrSig3dPhi04_b,Folder);
  sprintf(name,"JetRes_DetEta_j5_FrSig3dPhi01_b");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>5 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j5_FrSig3dPhi01_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_FrSig3dPhi01_b,Folder);
  sprintf(name,"JetRes_DetEta_j10_FrSig3dPhi01_b");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>10 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j10_FrSig3dPhi01_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_FrSig3dPhi01_b,Folder);
  sprintf(name,"JetRes_DetEta_j15_FrSig3dPhi01_b");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>15 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j15_FrSig3dPhi01_b=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_FrSig3dPhi01_b,Folder);


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

  //------------------------------ Jet Resolution studies --------------------------------
  sprintf(name,"JetRes_MetSig_EtaDet_j5_a");
  sprintf(title,"#slash{E}_{T} significance, #slash{E}_{T}/#sigma_{E_{T}}, vs. #eta_{det} for E_{T}^{lev6}>5 GeV");
  Hist.fJetRes_MetSig_EtaDet_j5_a= new TProfile(name,title,104,0.0,2.6,0.0,20.0);  
  AddHistogram(Hist.fJetRes_MetSig_EtaDet_j5_a,Folder);
  sprintf(name,"JetRes_MetSig_EtaDet_j10_a");
  sprintf(title,"#slash{E}_{T} significance, #slash{E}_{T}/#sigma_{E_{T}}, vs. #eta_{det} for E_{T}^{lev6}>10 GeV");
  Hist.fJetRes_MetSig_EtaDet_j10_a= new TProfile(name,title,104,0.0,2.6,0.0,20.0);  
  AddHistogram(Hist.fJetRes_MetSig_EtaDet_j10_a,Folder);
  sprintf(name,"JetRes_MetSig_EtaDet_j15_a");
  sprintf(title,"#slash{E}_{T} significance, #slash{E}_{T}/#sigma_{E_{T}}, vs. #eta_{det} for E_{T}^{lev6}>15 GeV");
  Hist.fJetRes_MetSig_EtaDet_j15_a= new TProfile(name,title,104,0.0,2.6,0.0,20.0);  
  AddHistogram(Hist.fJetRes_MetSig_EtaDet_j15_a,Folder);
  sprintf(name,"JetRes_DetEta_j5_a");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>5 GeV");
  Hist.fJetRes_DetEta_j5_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_a,Folder);
  sprintf(name,"JetRes_DetEta_j10_a");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>10 GeV");
  Hist.fJetRes_DetEta_j10_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_a,Folder);
  sprintf(name,"JetRes_DetEta_j15_a");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>15 GeV");
  Hist.fJetRes_DetEta_j15_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_a,Folder);
  sprintf(name,"JetRes_DetEta_j5_Sig3dPhi04_a");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>5 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j5_Sig3dPhi04_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_Sig3dPhi04_a,Folder);
  sprintf(name,"JetRes_DetEta_j10_Sig3dPhi04_a");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>10 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j10_Sig3dPhi04_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_Sig3dPhi04_a,Folder);
  sprintf(name,"JetRes_DetEta_j15_Sig3dPhi04_a");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>15 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j15_Sig3dPhi04_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_Sig3dPhi04_a,Folder);
  sprintf(name,"JetRes_DetEta_j5_Sig3dPhi01_a");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>5 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j5_Sig3dPhi01_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_Sig3dPhi01_a,Folder);
  sprintf(name,"JetRes_DetEta_j10_Sig3dPhi01_a");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>10 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j10_Sig3dPhi01_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_Sig3dPhi01_a,Folder);
  sprintf(name,"JetRes_DetEta_j15_Sig3dPhi01_a");
  sprintf(title,"#eta_{det} of leading jet with E_{T}^{lev6}>15 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j15_Sig3dPhi01_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_Sig3dPhi01_a,Folder);
  sprintf(name,"JetRes_DetEta_j5_FrSig3dPhi04_a");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>5 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j5_FrSig3dPhi04_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_FrSig3dPhi04_a,Folder);
  sprintf(name,"JetRes_DetEta_j10_FrSig3dPhi04_a");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>10 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j10_FrSig3dPhi04_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_FrSig3dPhi04_a,Folder);
  sprintf(name,"JetRes_DetEta_j15_FrSig3dPhi04_a");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>15 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.4");
  Hist.fJetRes_DetEta_j15_FrSig3dPhi04_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_FrSig3dPhi04_a,Folder);
  sprintf(name,"JetRes_DetEta_j5_FrSig3dPhi01_a");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>5 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j5_FrSig3dPhi01_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j5_FrSig3dPhi01_a,Folder);
  sprintf(name,"JetRes_DetEta_j10_FrSig3dPhi01_a");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>10 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j10_FrSig3dPhi01_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j10_FrSig3dPhi01_a,Folder);
  sprintf(name,"JetRes_DetEta_j15_FrSig3dPhi01_a");
  sprintf(title,"Fraction of leading jets vs #eta_{det}: E_{T}^{lev6}>15 GeV, MetSig>3 and |#Delta#phi_{#slash{E}_{T}-jet}|<0.1");
  Hist.fJetRes_DetEta_j15_FrSig3dPhi01_a=new TH1F(name,title,104,0.0,2.6);
  AddHistogram(Hist.fJetRes_DetEta_j15_FrSig3dPhi01_a,Folder);


  return;
}
/*}}}*/

//_____________________________________________________________________________
//  CREATE MET HISTOGRAMS       
//_____________________________________________________________________________
void TMyJetFilterModule::BookMetStudyHistograms(MetStudyHisto_t& Hist,
										const char* Folder, const char* algoname) 
/*{{{*/
{

  char name [200];
  char title[200];

  sprintf(name,"SumEtCorr");
  sprintf(title," %s%s",algoname,": corrected #sum E_{T}");
  Hist.fSumEtCorr=new TH1F(name,title,2000,0.0,2000.0);
  AddHistogram(Hist.fSumEtCorr,Folder);
  sprintf(name,"SqrtSumEtCorr");
  sprintf(title," %s%s",algoname,": corrected #sqrt{#sum E_{T}}");
  Hist.fSqrtSumEtCorr=new TH1F(name,title,180,0.0,45.0);
  AddHistogram(Hist.fSqrtSumEtCorr,Folder);

  sprintf(name,"SumEtCorr_noJet_vx1");
  sprintf(title," %s%s",algoname,": corrected #sum E_{T}, N_{jet15}=0, N_{vx}=1");
  Hist.fSumEtCorr_noJet_vx1=new TH1F(name,title,2000,0.0,2000.0);
  AddHistogram(Hist.fSumEtCorr_noJet_vx1,Folder);
  sprintf(name,"SumEtCorr_noJet_vx2");
  sprintf(title," %s%s",algoname,": corrected #sum E_{T}, N_{jet15}=0, N_{vx}=2");
  Hist.fSumEtCorr_noJet_vx2=new TH1F(name,title,2000,0.0,2000.0);
  AddHistogram(Hist.fSumEtCorr_noJet_vx2,Folder);
  sprintf(name,"SumEtCorr_noJet_vx3");
  sprintf(title," %s%s",algoname,": corrected #sum E_{T}, N_{jet15}=0, N_{vx}=3");
  Hist.fSumEtCorr_noJet_vx3=new TH1F(name,title,2000,0.0,2000.0);
  AddHistogram(Hist.fSumEtCorr_noJet_vx3,Folder);
  sprintf(name,"SumEtCorrNoJet_vs_Nvx");
  sprintf(title," %s%s",algoname,": corrected #sum E_{T} .vs. N_{vx} for N_{jet15}=0 events");
  Hist.fSumEtCorrNoJet_vs_Nvx=new TProfile(name,title,11,-0.5,10.5,0.0,2000.0);
  AddHistogram(Hist.fSumEtCorrNoJet_vs_Nvx,Folder);


  sprintf(name,"SumEtCorr_withJet");
  sprintf(title," %s%s",algoname,": corrected #sum E_{T} for events with N_{jet}(E_{T}>Threshold)>0");
  Hist.fSumEtCorr_withJet=new TH1F(name,title,2000,0.0,2000.0);
  AddHistogram(Hist.fSumEtCorr_withJet,Folder);
  sprintf(name,"SqrtSumEtCorr_withJet");
  sprintf(title," %s%s",algoname,": corrected #sqrt{#sum E_{T}} for events with N_{jet}(E_{T}>Threshold)>0");
  Hist.fSqrtSumEtCorr_withJet=new TH1F(name,title,180,0.0,45.0);
  AddHistogram(Hist.fSqrtSumEtCorr_withJet,Folder);
  sprintf(name,"SumEtCorr_noJet");
  sprintf(title," %s%s",algoname,": corrected #sum E_{T} for events with N_{jet}(E_{T}>Threshold)=0");
  Hist.fSumEtCorr_noJet=new TH1F(name,title,2000,0.0,2000.0);
  AddHistogram(Hist.fSumEtCorr_noJet,Folder);
  sprintf(name,"SqrtSumEtCorr_noJet");
  sprintf(title," %s%s",algoname,": corrected #sqrt{#sum E_{T}} for events with N_{jet}(E_{T}>Threshold)=0");
  Hist.fSqrtSumEtCorr_noJet=new TH1F(name,title,180,0.0,45.0);
  AddHistogram(Hist.fSqrtSumEtCorr_noJet,Folder);

  sprintf(name,"MetRaw");
  sprintf(title," %s%s",algoname,": MEt_4 before corr. for jets (corrected events)");
  Hist.fMetRaw=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fMetRaw,Folder);
  sprintf(name,"MetCorr");
  sprintf(title," %s%s",algoname,": MEt_4 after corr. for jets (corrected events)");
  Hist.fMetCorr=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fMetCorr,Folder);
  sprintf(name,"MetPhiRaw");
  sprintf(title," %s%s",algoname,": #phi of MEt_4 before corr. for jets (corrected events)");
  Hist.fMetPhiRaw=new TH1F(name,title,128,0.0,6.4);
  AddHistogram(Hist.fMetPhiRaw,Folder);
  sprintf(name,"MetPhiCorr");
  sprintf(title," %s%s",algoname,": #phi of MEt_4 after corr. for jets (corrected events)");
  Hist.fMetPhiCorr=new TH1F(name,title,128,0.0,6.4);
  AddHistogram(Hist.fMetPhiCorr,Folder);
  sprintf(name,"MetRawAll");
  sprintf(title," %s%s",algoname,": MEt_4 before corr. for jets (all events)" );
  Hist.fMetRawAll=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fMetRawAll,Folder);
  sprintf(name,"MetCorrAll");
  sprintf(title," %s%s",algoname,": MEt_4 after corr. for jets (all events)");
  Hist.fMetCorrAll=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fMetCorrAll,Folder);
  sprintf(name,"MetPhiRawAll");
  sprintf(title," %s%s",algoname,": #phi of MEt_4 before corr. for jets (all events)");
  Hist.fMetPhiRawAll=new TH1F(name,title,128,0.0,6.4);
  AddHistogram(Hist.fMetPhiRawAll,Folder);
  sprintf(name,"MetPhiCorrAll");
  sprintf(title," %s%s",algoname,": #phi of MEt_4 after corr. for jets (all events)");
  Hist.fMetPhiCorrAll=new TH1F(name,title,128,0.0,6.4);
  AddHistogram(Hist.fMetPhiCorrAll,Folder);
  sprintf(name,"MetRawAll_X");
  sprintf(title," %s%s",algoname,": MEtX_4 before corr. for jets (all events)" );
  Hist.fMetRawAll_X=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fMetRawAll_X,Folder);
  sprintf(name,"MetCorrAll_X");
  sprintf(title," %s%s",algoname,": MEtX_4 after corr. for jets (all events)");
  Hist.fMetCorrAll_X=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fMetCorrAll_X,Folder);
  sprintf(name,"MetRawAll_Y");
  sprintf(title," %s%s",algoname,": MEtY_4 before corr. for jets (all events)" );
  Hist.fMetRawAll_Y=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fMetRawAll_Y,Folder);
  sprintf(name,"MetCorrAll_Y");
  sprintf(title," %s%s",algoname,": MEtY_4 after corr. for jets (all events)");
  Hist.fMetCorrAll_Y=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fMetCorrAll_Y,Folder);
  sprintf(name,"Met4VsSqrtCorrSumEt");
  sprintf(title," %s%s",algoname,": MEt_4 vs #sqrt{#sum E_{T}^{corr}} for all events");
  Hist.fMet4VsSqrtCorrSumEt=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsSqrtCorrSumEt,Folder);
  sprintf(name,"Met4VsCorrSumEt");
  sprintf(title," %s%s",algoname,": MEt_4 vs #sum E_{T}^{corr} for all events");
  Hist.fMet4VsCorrSumEt=new TProfile(name,title,40,0.0,400.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsCorrSumEt,Folder);
  sprintf(name,"Met4VsSqrtCorrSumEt_nvx1");
  sprintf(title," %s%s",algoname,": MEt_4 vs #sqrt{#sum E_{T}^{corr}} for N_{vx12}=1");
  Hist.fMet4VsSqrtCorrSumEt_nvx1=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsSqrtCorrSumEt_nvx1,Folder);
  sprintf(name,"Met4VsSqrtCorrSumEt_nvx2");
  sprintf(title," %s%s",algoname,": MEt_4 vs #sqrt{#sum E_{T}^{corr}} for N_{vx12}=2");
  Hist.fMet4VsSqrtCorrSumEt_nvx2=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsSqrtCorrSumEt_nvx2,Folder);
  sprintf(name,"Met4VsSqrtCorrSumEt_nvx3");
  sprintf(title," %s%s",algoname,": MEt_4 vs #sqrt{#sum E_{T}^{corr}} for N_{vx12}=3");
  Hist.fMet4VsSqrtCorrSumEt_nvx3=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsSqrtCorrSumEt_nvx3,Folder);
  sprintf(name,"Met4VsSqrtCorrSumEt_nvx4");
  sprintf(title," %s%s",algoname,": MEt_4 vs #sqrt{#sum E_{T}^{corr}} for N_{vx12}=4");
  Hist.fMet4VsSqrtCorrSumEt_nvx4=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsSqrtCorrSumEt_nvx4,Folder);
  sprintf(name,"Met4VsSqrtCorrSumEt_nvx5");
  sprintf(title," %s%s",algoname,": MEt_4 vs #sqrt{#sum E_{T}^{corr}} for N_{vx12}=5");
  Hist.fMet4VsSqrtCorrSumEt_nvx5=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsSqrtCorrSumEt_nvx5,Folder);

  sprintf(name,"Met4VsNjet");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} vs N_{jet}");
  Hist.fMet4VsNjet=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsNjet,Folder);
  sprintf(name,"Met4VsNjet_th5");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} vs N_{jet} with E_{T}>5 GeV");
  Hist.fMet4VsNjet_th5=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsNjet_th5,Folder);
  sprintf(name,"Met4VsNjet_th10");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} vs N_{jet} with E_{T}>10 GeV");
  Hist.fMet4VsNjet_th10=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsNjet_th10,Folder);
  sprintf(name,"Met4VsNjet_th15");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} vs N_{jet} with E_{T}>15 GeV");
  Hist.fMet4VsNjet_th15=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsNjet_th15,Folder);
  sprintf(name,"Met4VsNjet_th20");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} vs N_{jet} with E_{T}>20 GeV");
  Hist.fMet4VsNjet_th20=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsNjet_th20,Folder);
  sprintf(name,"Met4VsNjet_th25");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} vs N_{jet} with E_{T}>25 GeV");
  Hist.fMet4VsNjet_th25=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsNjet_th25,Folder);
  sprintf(name,"MetGenVsNjet");
  sprintf(title," %s%s",algoname,": generated (def parametrization) ME_{T} vs N_{jet}");
 
  sprintf(name,"Met4SigVsNjet");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet}");
  Hist.fMet4SigVsNjet=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsNjet,Folder);
  sprintf(name,"Met4SigVsNjet_th5");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>5 GeV");
  Hist.fMet4SigVsNjet_th5=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsNjet_th5,Folder);
  sprintf(name,"Met4SigVsNjet_th10");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>10 GeV");
  Hist.fMet4SigVsNjet_th10=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsNjet_th10,Folder);
  sprintf(name,"Met4SigVsNjet_th15");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>15 GeV");
  Hist.fMet4SigVsNjet_th15=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsNjet_th15,Folder);
  sprintf(name,"Met4SigVsNjet_th20");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>20 GeV");
  Hist.fMet4SigVsNjet_th20=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsNjet_th20,Folder);
  sprintf(name,"Met4SigVsNjet_th25");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>25 GeV");
  Hist.fMet4SigVsNjet_th25=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsNjet_th25,Folder);

  sprintf(name,"Met4VsRun");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} Vs. Run for N_{vx12}=1");
  Hist.fMet4VsRun=new TProfile(name,title,100,130000.0,230000.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsRun,Folder);
  sprintf(name,"Met4SigVsRun");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} Vs. Run for N_{vx12}=1");
  Hist.fMet4SigVsRun=new TProfile(name,title,100,130000.0,230000.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsRun,Folder);
  sprintf(name,"Nvx12VsRun");
  sprintf(title,"N_{vx12} Vs. Run");
  Hist.fNvx12VsRun=new TProfile(name,title,200,130000.0,230000.0,-1.0,100.0);
  AddHistogram(Hist.fNvx12VsRun,Folder);

  sprintf(name,"Met4VsNvx12");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} vs N_{vx12}");
  Hist.fMet4VsNvx12=new TProfile(name,title,10,-0.5,9.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsNvx12,Folder);
  sprintf(name,"Met4SigVsNvx12");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} vs N_{vx12}");
  Hist.fMet4SigVsNvx12=new TProfile(name,title,10,-0.5,9.5,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsNvx12,Folder);

  sprintf(name,"SumEtJetFrac");
  sprintf(title," %s%s",algoname,": #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T}");
  Hist.fSumEtJetFrac=new TH1F(name,title,100,0.0,1.0);
  AddHistogram(Hist.fSumEtJetFrac,Folder);
  sprintf(name,"SumEtJet");
  sprintf(title," %s%s",algoname,": corrected #sum E^{jets}_{T}");
  Hist.fSumEtJet=new TH1F(name,title,2000,0.0,2000.0);
  AddHistogram(Hist.fSumEtJet,Folder);
  sprintf(name,"SqrtSumEtJet");
  sprintf(title," %s%s",algoname,": corrected #sqrt{#sum E^{jets}_{T}}");
  Hist.fSqrtSumEtJet=new TH1F(name,title,180,0.0,45.0);
  AddHistogram(Hist.fSqrtSumEtJet,Folder);
  sprintf(name,"Met4VsJetFr");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} vs #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T}");
  Hist.fMet4VsJetFr=new TProfile(name,title,20,0.0,1.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsJetFr,Folder);
  sprintf(name,"Met4SigVsJetFr");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} vs #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T}");
  Hist.fMet4SigVsJetFr=new TProfile(name,title,20,0.0,1.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsJetFr,Folder);


  sprintf(name,"Met4VsSqrtSumEtJet");
  sprintf(title," %s%s",algoname,": meth4 ME_{T} vs #sqrt{#Sigma E^{jets}_{T}}");
  Hist.fMet4VsSqrtSumEtJet=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4VsSqrtSumEtJet,Folder);
  sprintf(name,"Met4SigVsSqrtSumEtJet");
  sprintf(title," %s%s",algoname,": meth4 ME_{T}/#sqrt{#Sigma E_{T}} vs #sqrt{#Sigma E^{jets}_{T}}");
  Hist.fMet4SigVsSqrtSumEtJet=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMet4SigVsSqrtSumEtJet,Folder);
  sprintf(name,"JetFrVsSumEt");
  sprintf(title," %s%s",algoname,": #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T} vs #Sigma E^{tot}_{T}}");
  Hist.fJetFrVsSumEt=new TProfile(name,title,40,0.0,400.0,-1.0,10.0);
  AddHistogram(Hist.fJetFrVsSumEt,Folder);

  sprintf(name,"Met4X_vs_Met4Y");
  sprintf(title," %s%s",algoname,": meth4 ME_{T,X} vs. ME_{T,Y}");
  Hist.fMet4X_vs_Met4Y=new TH2F(name,title,40,-40.0,40.0,40,-40.0,40.0);
  AddHistogram(Hist.fMet4X_vs_Met4Y,Folder);

  for(int i=0; i<10; i++)
    {

      sprintf(name,"Met4[%i]",i);
      sprintf(title,"%s: meth4 ME_{T}, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4[i]=new TH1F(name,title,500,0.0,500.0);
      AddHistogram(Hist.fMet4[i],Folder);
      sprintf(name,"Met4Phi[%i]",i);
      sprintf(title,"%s: #phi of meth4 ME_{T}, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Phi[i]=new TH1F(name,title,128,0.0,6.4);
      AddHistogram(Hist.fMet4Phi[i],Folder);

      //---- any number of vertices
      sprintf(name,"Met4X[%i]",i);
      sprintf(title,"%s: meth4 MEtX, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4X[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4X[i],Folder);
      sprintf(name,"Met4Y[%i]",i);
      sprintf(title,"%s: meth4 MEtY, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Y[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4Y[i],Folder);
      //---- nvx12=1
      sprintf(name,"Met4X_nvx1[%i]",i);
      sprintf(title,"%s: meth4 MEtX for N_{vx12}=1, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4X_nvx1[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4X_nvx1[i],Folder);
      sprintf(name,"Met4Y_nvx1[%i]",i);
      sprintf(title,"%s: meth4 MEtY for N_{vx12}=1, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Y_nvx1[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4Y_nvx1[i],Folder);
      //---- nvx12=2
      sprintf(name,"Met4X_nvx2[%i]",i);
      sprintf(title,"%s: meth4 MEtX for N_{vx12}=2, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4X_nvx2[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4X_nvx2[i],Folder);
      sprintf(name,"Met4Y_nvx2[%i]",i);
      sprintf(title,"%s: meth4 MEtY for N_{vx12}=2, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Y_nvx2[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4Y_nvx2[i],Folder);
      //---- nvx12=3
      sprintf(name,"Met4X_nvx3[%i]",i);
      sprintf(title,"%s: meth4 MEtX for N_{vx12}=3, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4X_nvx3[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4X_nvx3[i],Folder);
      sprintf(name,"Met4Y_nvx3[%i]",i);
      sprintf(title,"%s: meth4 MEtY for N_{vx12}=3, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Y_nvx3[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4Y_nvx3[i],Folder);
      //---- nvx12=4
      sprintf(name,"Met4X_nvx4[%i]",i);
      sprintf(title,"%s: meth4 MEtX for N_{vx12}=4, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4X_nvx4[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4X_nvx4[i],Folder);
      sprintf(name,"Met4Y_nvx4[%i]",i);
      sprintf(title,"%s: meth4 MEtY for N_{vx12}=4, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Y_nvx4[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4Y_nvx4[i],Folder);
      //---- nvx12=5
      sprintf(name,"Met4X_nvx5[%i]",i);
      sprintf(title,"%s: meth4 MEtX for N_{vx12}=5, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4X_nvx5[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4X_nvx5[i],Folder);
      sprintf(name,"Met4Y_nvx5[%i]",i);
      sprintf(title,"%s: meth4 MEtY for N_{vx12}=5, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Y_nvx5[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4Y_nvx5[i],Folder);
      //---- nvx12>1
      sprintf(name,"Met4X_nvxM1[%i]",i);
      sprintf(title,"%s: meth4 MEtX for N_{vx12}>1, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4X_nvxM1[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4X_nvxM1[i],Folder);
      sprintf(name,"Met4Y_nvxM1[%i]",i);
      sprintf(title,"%s: meth4 MEtY for N_{vx12}>1, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Y_nvxM1[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4Y_nvxM1[i],Folder);

      //---- Njet_thr*=0
      sprintf(name,"Met4X_noJet[%i]",i);
      sprintf(title,"%s: meth4 MEtX for N_{jet}(E_{T}>Threshold)=0, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4X_noJet[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4X_noJet[i],Folder);
      sprintf(name,"Met4Y_noJet[%i]",i);
      sprintf(title,"%s: meth4 MEtY for N_{jet}(E_{T}>Threshold)=0, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Y_noJet[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4Y_noJet[i],Folder);
      //---- Njet_thr*>0
      sprintf(name,"Met4X_withJet[%i]",i);
      sprintf(title,"%s: meth4 MEtX for N_{jet}(E_{T}>Threshold)>0, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4X_withJet[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4X_withJet[i],Folder);
      sprintf(name,"Met4Y_withJet[%i]",i);
      sprintf(title,"%s: meth4 MEtY for N_{jet}(E_{T}>Threshold)>0, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
      Hist.fMet4Y_withJet[i]=new TH1F(name,title,400,-200.0,200.0);
      AddHistogram(Hist.fMet4Y_withJet[i],Folder);
    }

  //______________________________________________________________________ histograms for MetModel-2
  sprintf(name,"GenV2Met_def");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T} (default parametrization)");
  Hist.fGenV2Met_def=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenV2Met_def,Folder);
  sprintf(name,"GenV2MetPhi_def");
  sprintf(title," %s%s",algoname,": #phi of generated Model-2 ME_{T} (default parametrization)");
  Hist.fGenV2MetPhi_def=new TH1F(name,title,128,0.0,6.4);
  AddHistogram(Hist.fGenV2MetPhi_def,Folder);
  sprintf(name,"GenV2MetX_def");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,X} (default parametrization)");
  Hist.fGenV2MetX_def=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetX_def,Folder);
  sprintf(name,"GenV2MetY_def");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,Y} (default parametrization)");
  Hist.fGenV2MetY_def=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetY_def,Folder);

  sprintf(name,"GenV2Met_devm");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T} (default-sigmaJet parametrization)");
  Hist.fGenV2Met_devm=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenV2Met_devm,Folder);
  sprintf(name,"GenV2MetPhi_devm");
  sprintf(title," %s%s",algoname,": #phi of generated Model-2 ME_{T} (default-sigmaJet parametrization)");
  Hist.fGenV2MetPhi_devm=new TH1F(name,title,128,0.0,6.4);
  AddHistogram(Hist.fGenV2MetPhi_devm,Folder);
  sprintf(name,"GenV2MetX_devm");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,X} (default-sigmaJet parametrization)");
  Hist.fGenV2MetX_devm=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetX_devm,Folder);
  sprintf(name,"GenV2MetY_devm");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,Y} (default-sigmaJet parametrization)");
  Hist.fGenV2MetY_devm=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetY_devm,Folder);

  sprintf(name,"GenV2Met_devp");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T} (default+sigmaJet parametrization)");
  Hist.fGenV2Met_devp=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenV2Met_devp,Folder);
  sprintf(name,"GenV2MetPhi_devp");
  sprintf(title," %s%s",algoname,": #phi of generated Model-2 ME_{T} (default+sigmaJet parametrization)");
  Hist.fGenV2MetPhi_devp=new TH1F(name,title,128,0.0,6.4);
  AddHistogram(Hist.fGenV2MetPhi_devp,Folder);
  sprintf(name,"GenV2MetX_devp");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,X} (default+sigmaJet parametrization)");
  Hist.fGenV2MetX_devp=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetX_devp,Folder);
  sprintf(name,"GenV2MetY_devp");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,Y} (default+sigmaJet parametrization)");
  Hist.fGenV2MetY_devp=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetY_devp,Folder);

  sprintf(name,"GenV2Met_devmUn");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T} (default-sigmaUncl parametrization)");
  Hist.fGenV2Met_devmUn=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenV2Met_devmUn,Folder);
  sprintf(name,"GenV2MetPhi_devmUn");
  sprintf(title," %s%s",algoname,": #phi of generated Model-2 ME_{T} (default-sigmaUncl parametrization)");
  Hist.fGenV2MetPhi_devmUn=new TH1F(name,title,128,0.0,6.4);
  AddHistogram(Hist.fGenV2MetPhi_devmUn,Folder);
  sprintf(name,"GenV2MetX_devmUn");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,X} (default-sigmaUncl parametrization)");
  Hist.fGenV2MetX_devmUn=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetX_devmUn,Folder);
  sprintf(name,"GenV2MetY_devmUn");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,Y} (default-sigmaUncl parametrization)");
  Hist.fGenV2MetY_devmUn=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetY_devmUn,Folder);

  sprintf(name,"GenV2Met_devpUn");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T} (default+sigmaUncl parametrization)");
  Hist.fGenV2Met_devpUn=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenV2Met_devpUn,Folder);
  sprintf(name,"GenV2MetPhi_devpUn");
  sprintf(title," %s%s",algoname,": #phi of generated Model-2 ME_{T} (default+sigmaUncl parametrization)");
  Hist.fGenV2MetPhi_devpUn=new TH1F(name,title,128,0.0,6.4);
  AddHistogram(Hist.fGenV2MetPhi_devpUn,Folder);
  sprintf(name,"GenV2MetX_devpUn");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,X} (default+sigmaUncl parametrization)");
  Hist.fGenV2MetX_devpUn=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetX_devpUn,Folder);
  sprintf(name,"GenV2MetY_devpUn");
  sprintf(title," %s%s",algoname,": generated Model-2 ME_{T,Y} (default+sigmaUncl parametrization)");
  Hist.fGenV2MetY_devpUn=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fGenV2MetY_devpUn,Folder);


  sprintf(name,"MetGenV2VsNjet");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} vs N_{jet}");
  Hist.fMetGenV2VsNjet=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsNjet,Folder);
  sprintf(name,"MetGenV2VsNjet_th5");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} vs N_{jet} with E_{T}>5 GeV");
  Hist.fMetGenV2VsNjet_th5=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsNjet_th5,Folder);
  sprintf(name,"MetGenV2VsNjet_th10");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} vs N_{jet} with E_{T}>10 GeV");
  Hist.fMetGenV2VsNjet_th10=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsNjet_th10,Folder);
  sprintf(name,"MetGenV2VsNjet_th15");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} vs N_{jet} with E_{T}>15 GeV");
  Hist.fMetGenV2VsNjet_th15=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsNjet_th15,Folder);
  sprintf(name,"MetGenV2VsNjet_th20");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} vs N_{jet} with E_{T}>20 GeV");
  Hist.fMetGenV2VsNjet_th20=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsNjet_th20,Folder);
  sprintf(name,"MetGenV2VsNjet_th25");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} vs N_{jet} with E_{T}>25 GeV");
  Hist.fMetGenV2VsNjet_th25=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsNjet_th25,Folder);

  sprintf(name,"MetGenV2SigVsNjet");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet}");
  Hist.fMetGenV2SigVsNjet=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsNjet,Folder);
  sprintf(name,"MetGenV2SigVsNjet_th5");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>5 GeV");
  Hist.fMetGenV2SigVsNjet_th5=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsNjet_th5,Folder);
  sprintf(name,"MetGenV2SigVsNjet_th10");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>10 GeV");
  Hist.fMetGenV2SigVsNjet_th10=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsNjet_th10,Folder);
  sprintf(name,"MetGenV2SigVsNjet_th15");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>15 GeV");
  Hist.fMetGenV2SigVsNjet_th15=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsNjet_th15,Folder);
  sprintf(name,"MetGenV2SigVsNjet_th20");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>20 GeV");
  Hist.fMetGenV2SigVsNjet_th20=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsNjet_th20,Folder);
  sprintf(name,"MetGenV2SigVsNjet_th25");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} vs N_{jet} with E_{T}>25 GeV");
  Hist.fMetGenV2SigVsNjet_th25=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsNjet_th25,Folder);

  sprintf(name,"MetGenV2VsRun");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} Vs. Run for N_{vx12}=1");
  Hist.fMetGenV2VsRun=new TProfile(name,title,100,130000.0,230000.0,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsRun,Folder);
  sprintf(name,"MetGenV2SigVsRun");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} Vs. Run for N_{vx12}=1");
  Hist.fMetGenV2SigVsRun=new TProfile(name,title,100,130000.0,230000.0,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsRun,Folder);
  sprintf(name,"MetGenV2VsNvx12");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} vs N_{vx12}");
  Hist.fMetGenV2VsNvx12=new TProfile(name,title,10,-0.5,9.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsNvx12,Folder);
  sprintf(name,"MetGenV2SigVsNvx12");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} vs N_{vx12}");
  Hist.fMetGenV2SigVsNvx12=new TProfile(name,title,10,-0.5,9.5,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsNvx12,Folder);
    
  sprintf(name,"MetGenV2VsJetFr");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} vs #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T}");
  Hist.fMetGenV2VsJetFr=new TProfile(name,title,20,0.0,1.0,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsJetFr,Folder);
  sprintf(name,"MetGenV2SigVsJetFr");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} vs #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T}");
  Hist.fMetGenV2SigVsJetFr=new TProfile(name,title,20,0.0,1.0,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsJetFr,Folder);
  sprintf(name,"MetGenV2VsSqrtSumEtJet");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T} vs #sqrt{#Sigma E^{jets}_{T}}");
  Hist.fMetGenV2VsSqrtSumEtJet=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2VsSqrtSumEtJet,Folder);
  sprintf(name,"MetGenV2SigVsSqrtSumEtJet");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T}/#sqrt{#Sigma E_{T}} vs #sqrt{#Sigma E^{jets}_{T}}");
  Hist.fMetGenV2SigVsSqrtSumEtJet=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
  AddHistogram(Hist.fMetGenV2SigVsSqrtSumEtJet,Folder);
  sprintf(name,"GenV2MetX_vs_MetY");
  sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) ME_{T,X} vs. ME_{T,Y}");
  Hist.fGenV2MetX_vs_MetY=new TH2F(name,title,40,-40.0,40.0,40,-40.0,40.0);
  AddHistogram(Hist.fGenV2MetX_vs_MetY,Folder);
  
  sprintf(name,"Met4GenV2Met");
  sprintf(title," %s%s",algoname,": #slash{E}_{T}-#slash{E}_{T}^{gen}");
  Hist.fMet4GenV2Met=new TH1F(name,title,400,-200.0,200.0);
  AddHistogram(Hist.fMet4GenV2Met,Folder);
  sprintf(name,"Met4_vs_Met4GenV2Met");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} .vs. #slash{E}_{T}-#slash{E}_{T}^{gen}");
  Hist.fMet4_vs_Met4GenV2Met=new TH2F(name,title,60,-150.0,150.0,30,0.0,150.0);
  AddHistogram(Hist.fMet4_vs_Met4GenV2Met,Folder);
  sprintf(name,"Met4_vs_GenV2Met");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} .vs. #slash{E}_{T}^{gen}");
  Hist.fMet4_vs_GenV2Met=new TH2F(name,title,50,0.0,250.0,50,0.0,250.0);
  AddHistogram(Hist.fMet4_vs_GenV2Met,Folder);



  return;
}
/*}}}*/


//_____________________________________________________________________________
//  CREATE PHOTON-JET MATCH HISTOGRAMS
//_____________________________________________________________________________
void TMyJetFilterModule::BookMatchStudyHistograms(MatchStudyHisto_t& Hist,
											const char* Folder, const char* algoname)
/*{{{*/
{

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

/*}}}*/


//_____________________________________________________________________________
//   CREATE MET CLEANUP HISTOGRAMS (FOR STUDIES)
//_____________________________________________________________________________
void TMyJetFilterModule::BookMetCleanupHistograms(MetCleanupHisto_t& Hist,
											const char* Folder, const char* algoname)
/*{{{*/
{

  char name [200];
  char title[200];

  sprintf(name,"CleanupDelPhiVsEtaDet_jet");
  sprintf(title," %s%s",algoname,": d#phi=#phi_{jet}-#phi_{#slash{E}_{T}} vs. #eta^{det}_{jet}");          
  Hist.fCleanupDelPhiVsEtaDet_jet= new TH2F(name,title,37,-3.7,3.7,16,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhiVsEtaDet_jet,Folder);
  sprintf(name,"CleanupDelPhiVsEtaDet_em");
  sprintf(title," %s%s",algoname,": d#phi=#phi_{e,#gamma}-#phi_{#slash{E}_{T}} vs. #eta^{det}_{e,#gamma}");          
  Hist.fCleanupDelPhiVsEtaDet_em= new TH2F(name,title,37,-3.7,3.7,16,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhiVsEtaDet_em,Folder);
  
  sprintf(name,"CleanupMetEtemVsDelPhi");
  sprintf(title," %s%s",algoname,": #slash{E}_{T}/E_{T}^{e,#gamma} vs. d#phi=#phi_{e,#gamma}-#phi_{#slash{E}_{T}}");          
  Hist.fCleanupMetEtemVsDelPhi= new TH2F(name,title,16,0.0,TMath::Pi(),50,0.0,10.0);          
  AddHistogram(Hist.fCleanupMetEtemVsDelPhi,Folder);
  sprintf(name,"CleanupDelPhi_metem");
  sprintf(title," %s%s",algoname,": d#phi=#phi_{e,#gamma}-#phi_{#slash{E}_{T}}");          
  Hist.fCleanupDelPhi_metem= new TH1F(name,title,128,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_metem,Folder);
  sprintf(name,"CleanupMetEtem");
  sprintf(title," %s%s",algoname,": #slash{E}_{T}/E_{T}^{e,#gamma}");          
  Hist.fCleanupMetEtem= new TH1F(name,title,800,0.0,20.0);          
  AddHistogram(Hist.fCleanupMetEtem,Folder);
  
  for(int i=0; i<3; i++)
    {
      char range[200];
      if(i==0) sprintf(range,", jet in crack");
      if(i==1) sprintf(range,", jet away from crack");
      if(i==2) sprintf(range,", jet in |#eta_{det}|>2.6");
      
      sprintf(name,"CleanupMetEtjetVsDelPhi[%i]",i);
      sprintf(title," %s: #slash{E}_{T}/E_{T}^{jet} vs. d#phi=#phi_{jet}-#phi_{#slash{E}_{T}}%s",algoname,range);          
      Hist.fCleanupMetEtjetVsDelPhi[i]= new TH2F(name,title,16,0.0,TMath::Pi(),50,0.0,10.0);          
      AddHistogram(Hist.fCleanupMetEtjetVsDelPhi[i],Folder);
      sprintf(name,"CleanupDelPhi_metjet[%i]",i);
      sprintf(title," %s: d#phi=#phi_{jet}-#phi_{#slash{E}_{T}}%s",algoname,range);          
      Hist.fCleanupDelPhi_metjet[i]= new TH1F(name,title,128,0.0,TMath::Pi());          
      AddHistogram(Hist.fCleanupDelPhi_metjet[i],Folder);
      sprintf(name,"CleanupMetEtjet[%i]",i);
      sprintf(title," %s: #slash{E}_{T}/E_{T}^{jet}%s",algoname,range);          
      Hist.fCleanupMetEtjet[i]= new TH1F(name,title,800,0.0,20.0);          
      AddHistogram(Hist.fCleanupMetEtjet[i],Folder);
    }
  
  sprintf(name,"CleanupDelPhi_RawMetRawJet3");
  sprintf(title," %s%s",algoname,": #Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}^{raw}}, closest jet with E_{T}^{raw}>3 GeV");          
  Hist.fCleanupDelPhi_RawMetRawJet3= new TH1F(name,title,126,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_RawMetRawJet3,Folder);
  sprintf(name,"CleanupDelPhi_RawMetRawJet5");
  sprintf(title," %s%s",algoname,": #Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}^{raw}}, closest jet with E_{T}^{raw}>5 GeV");          
  Hist.fCleanupDelPhi_RawMetRawJet5= new TH1F(name,title,126,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_RawMetRawJet5,Folder);
  sprintf(name,"CleanupDelPhi_RawMetRawJet10");
  sprintf(title," %s%s",algoname,": #Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}^{raw}}, closest jet with E_{T}^{raw}>10 GeV");          
  Hist.fCleanupDelPhi_RawMetRawJet10= new TH1F(name,title,126,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_RawMetRawJet10,Folder);
  sprintf(name,"CleanupDelPhi_RawMet1stRawJet");
  sprintf(title," %s%s",algoname,": #Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}^{raw}}, 1^{st} jet with E_{T}^{raw}>3 GeV");          
  Hist.fCleanupDelPhi_RawMet1stRawJet= new TH1F(name,title,126,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_RawMet1stRawJet,Folder);
  
  //----------------------------------------- histograms added on 04/03/07
  sprintf(name,"CleanupDelPhi_metjet5");
  sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>5 GeV");          
  Hist.fCleanupDelPhi_metjet5= new TH1F(name,title,128,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_metjet5,Folder);
  sprintf(name,"CleanupDelPhi_metjet15");
  sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>15 GeV");          
  Hist.fCleanupDelPhi_metjet15= new TH1F(name,title,128,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_metjet15,Folder);
  sprintf(name,"CleanupDelPhi_metjet25");
  sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>25 GeV");          
  Hist.fCleanupDelPhi_metjet25= new TH1F(name,title,128,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_metjet25,Folder);
  sprintf(name,"CleanupDelPhi_met10jet15");
  sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>15 GeV and #slash{E}_{T}>10 GeV");          
  Hist.fCleanupDelPhi_met10jet15= new TH1F(name,title,128,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_met10jet15,Folder);
  sprintf(name,"CleanupDelPhi_metjet1st15");
  sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, 1^{st} jet with E_{T}^{lev6}>15 GeV");          
  Hist.fCleanupDelPhi_metjet1st15= new TH1F(name,title,128,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_metjet1st15,Folder);
  sprintf(name,"CleanupDelPhi_met10jet1st15");
  sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, 1^{st} jet with E_{T}^{lev6}>15 GeV and #slash{E}_{T}>10 GeV");          
  Hist.fCleanupDelPhi_met10jet1st15= new TH1F(name,title,128,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_met10jet1st15,Folder);
  sprintf(name,"CleanupDelPhi_metjet15_dPhiMin");
  sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>15 GeV and min(#Delta#phi)");          
  Hist.fCleanupDelPhi_metjet15_dPhiMin= new TH1F(name,title,128,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_metjet15_dPhiMin,Folder);
  sprintf(name,"CleanupDelPhi_met10jet15_dPhiMin");
  sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>15 GeV, min(#Delta#phi), and #slash{E}_{T}>10 GeV");          
  Hist.fCleanupDelPhi_met10jet15_dPhiMin= new TH1F(name,title,128,0.0,TMath::Pi());          
  AddHistogram(Hist.fCleanupDelPhi_met10jet15_dPhiMin,Folder);
  
  return;
}

/*}}}*/


//_____________________________________________________________________________
//  CREATE MET ???????? HISTOGRAMS
//_____________________________________________________________________________
void TMyJetFilterModule::BookMetProbHistograms(MetProbHisto_t& Hist,
											const char* Folder, const char* algoname)
/*{{{*/
{
  char name [200];
  char title[200];

  sprintf(name,"DataMet_vs_GenSigmaMet");
  sprintf(title," %s%s",algoname,": #slash{E}_{T}.vs.#sigma_{#slash{E}_{T}}, data events");
  Hist.fDataMet_vs_GenSigmaMet=new TH2F(name,title,50,0.0,100.0,60,0.0,300.0);
  AddHistogram(Hist.fDataMet_vs_GenSigmaMet,Folder);
  sprintf(name,"ToyEvntMet_vs_GenSigmaMet");
  sprintf(title," %s%s",algoname,": #slash{E}^{meas.}_{T}.vs.#sigma_{#slash{E}_{T}}, Toy MC events");
  Hist.fToyEvntMet_vs_GenSigmaMet=new TH2F(name,title,50,0.0,100.0,60,0.0,300.0);
  AddHistogram(Hist.fToyEvntMet_vs_GenSigmaMet,Folder);

  sprintf(name,"DataMet_highLnProb");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in data events with -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}})>2.56");
  Hist.fDataMet_highLnProb=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fDataMet_highLnProb,Folder);
  sprintf(name,"GenMet_def_highLnProb");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (def.) events with -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}})>2.56");
  Hist.fGenMet_def_highLnProb=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_def_highLnProb,Folder);
  sprintf(name,"GenMet_dvm_highLnProb");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (-#sigma) events with -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}})>2.56");
  Hist.fGenMet_dvm_highLnProb=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_dvm_highLnProb,Folder);
  sprintf(name,"GenMet_dvp_highLnProb");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (+#sigma) events with -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}})>2.56");
  Hist.fGenMet_dvp_highLnProb=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_dvp_highLnProb,Folder);

  sprintf(name,"DataMet_lowLnProb");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in data events with -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}})<2.56");
  Hist.fDataMet_lowLnProb=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fDataMet_lowLnProb,Folder);
  sprintf(name,"GenMet_def_lowLnProb");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (def.) events with -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}})<2.56");
  Hist.fGenMet_def_lowLnProb=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_def_lowLnProb,Folder);
  sprintf(name,"GenMet_dvm_lowLnProb");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (-#sigma) events with -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}})<2.56");
  Hist.fGenMet_dvm_lowLnProb=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_dvm_lowLnProb,Folder);
  sprintf(name,"GenMet_dvp_lowLnProb");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (+#sigma) events with -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}})<2.56");
  Hist.fGenMet_dvp_lowLnProb=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_dvp_lowLnProb,Folder);

  sprintf(name,"DataMet_highSigma");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in data events with (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}>3.0");
  Hist.fDataMet_highSigma=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fDataMet_highSigma,Folder);
  sprintf(name,"GenMet_def_highSigma");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (def.) events with (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}>3.0");
  Hist.fGenMet_def_highSigma=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_def_highSigma,Folder);
  sprintf(name,"GenMet_dvm_highSigma");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (-#sigma) events with (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}>3.0");
  Hist.fGenMet_dvm_highSigma=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_dvm_highSigma,Folder);
  sprintf(name,"GenMet_dvp_highSigma");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (+#sigma) events with (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}>3.0");
  Hist.fGenMet_dvp_highSigma=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_dvp_highSigma,Folder);

  sprintf(name,"DataMet_lowSigma");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in data events with (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}<3.0");
  Hist.fDataMet_lowSigma=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fDataMet_lowSigma,Folder);
  sprintf(name,"GenMet_def_lowSigma");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (def.) events with (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}<3.0");
  Hist.fGenMet_def_lowSigma=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_def_lowSigma,Folder);
  sprintf(name,"GenMet_dvm_lowSigma");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (-#sigma) events with (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}<3.0");
  Hist.fGenMet_dvm_lowSigma=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_dvm_lowSigma,Folder);
  sprintf(name,"GenMet_dvp_lowSigma");
  sprintf(title," %s%s",algoname,": #slash{E}_{T} in generated (+#sigma) events with (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}<3.0");
  Hist.fGenMet_dvp_lowSigma=new TH1F(name,title,500,0.0,500.0);
  AddHistogram(Hist.fGenMet_dvp_lowSigma,Folder);

  sprintf(name,"DataMetProb_def");
  sprintf(title," %s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}}) in data, default parametrization");          
  Hist.fDataMetProb_def= new TH1F(name,title,130,0.0,6.5);          
  AddHistogram(Hist.fDataMetProb_def,Folder);
  sprintf(name,"DataMetProb_dvm");
  sprintf(title," %s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}}) in data, -sigma parametrization");          
  Hist.fDataMetProb_dvm= new TH1F(name,title,130,0.0,6.5);          
  AddHistogram(Hist.fDataMetProb_dvm,Folder);
  sprintf(name,"DataMetProb_dvp");
  sprintf(title," %s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}}) in data, +sigma parametrization");          
  Hist.fDataMetProb_dvp= new TH1F(name,title,130,0.0,6.5);          
  AddHistogram(Hist.fDataMetProb_dvp,Folder);

  sprintf(name,"DataMetSigma_def");
  sprintf(title," %s%s",algoname,": (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}} in data, default parametrization");          
  Hist.fDataMetSigma_def= new TH1F(name,title,250,-3.0,22.0);          
  AddHistogram(Hist.fDataMetSigma_def,Folder);
  sprintf(name,"DataMetSigma_dvm");
  sprintf(title," %s%s",algoname,": (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}} in data, -sigma parametrization");          
  Hist.fDataMetSigma_dvm= new TH1F(name,title,250,-3.0,22.0);          
  AddHistogram(Hist.fDataMetSigma_dvm,Folder);
  sprintf(name,"DataMetSigma_dvp");
  sprintf(title," %s%s",algoname,": (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}} in data, +sigma parametrization");          
  Hist.fDataMetSigma_dvp= new TH1F(name,title,250,-3.0,22.0);          
  AddHistogram(Hist.fDataMetSigma_dvp,Folder);

  for(int i=0; i<20; i++)
    {
      char range[200];
      sprintf(range,"%4.1f<#slash{E}_{T}^{toy}<%4.1f",MetToy_min+(MetToy_max-MetToy_min)*i/20.0,MetToy_min+(MetToy_max-MetToy_min)*(i+1)/20.0);
      sprintf(name,"ToyMetProb_def[%i]",i);
      sprintf(title," %s%s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{meas.}>#slash{E}_{T}^{gen}}), default parametrization, ",range);          
      Hist.fToyMetProb_def[i]= new TH1F(name,title,130,0.0,6.5);          
      AddHistogram(Hist.fToyMetProb_def[i],Folder);
      sprintf(name,"ToyMetProb_dvm[%i]",i);
      sprintf(title," %s%s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{meas.}>#slash{E}_{T}^{gen}}), -sigma parametrization, ",range);          
      Hist.fToyMetProb_dvm[i]= new TH1F(name,title,130,0.0,6.5);          
      AddHistogram(Hist.fToyMetProb_dvm[i],Folder);
      sprintf(name,"ToyMetProb_dvp[%i]",i);
      sprintf(title," %s%s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{meas.}>#slash{E}_{T}^{gen}}), +sigma parametrization, ",range);          
      Hist.fToyMetProb_dvp[i]= new TH1F(name,title,130,0.0,6.5);          
      AddHistogram(Hist.fToyMetProb_dvp[i],Folder);

      sprintf(name,"ToyMetSigma_def[%i]",i);
      sprintf(title," %s%s%s",algoname,": (#slash{E}_{T}^{meas.}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}, default parametrization, ",range);          
      Hist.fToyMetSigma_def[i]= new TH1F(name,title,250,-3.0,22.0);          
      AddHistogram(Hist.fToyMetSigma_def[i],Folder);
      sprintf(name,"ToyMetSigma_dvm[%i]",i);
      sprintf(title," %s%s%s",algoname,": (#slash{E}_{T}^{meas.}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}, -sigma parametrization, ",range);          
      Hist.fToyMetSigma_dvm[i]= new TH1F(name,title,250,-3.0,22.0);          
      AddHistogram(Hist.fToyMetSigma_dvm[i],Folder);
      sprintf(name,"ToyMetSigma_dvp[%i]",i);
      sprintf(title," %s%s%s",algoname,": (#slash{E}_{T}^{meas.}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}, +sigma parametrization, ",range);          
      Hist.fToyMetSigma_dvp[i]= new TH1F(name,title,250,-3.0,22.0);          
      AddHistogram(Hist.fToyMetSigma_dvp[i],Folder);
    }

  sprintf(name,"DataMetProb_Met_def");
  sprintf(title," %s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}}).vs.#slash{E}_{T}^{data}, default parametrization");          
  Hist.fDataMetProb_Met_def= new TProfile(name,title,40,0.0,200.0,-1.0,10.0);          
  AddHistogram(Hist.fDataMetProb_Met_def,Folder);
  sprintf(name,"DataMetProb_Met_dvm");
  sprintf(title," %s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}}).vs.#slash{E}_{T}^{data}, -sigma parametrization");          
  Hist.fDataMetProb_Met_dvm= new TProfile(name,title,40,0.0,200.0,-1.0,10.0);          
  AddHistogram(Hist.fDataMetProb_Met_dvm,Folder);
  sprintf(name,"DataMetProb_Met_dvp");
  sprintf(title," %s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{data}>#slash{E}_{T}^{gen}}).vs.#slash{E}_{T}^{data}, +sigma parametrization");          
  Hist.fDataMetProb_Met_dvp= new TProfile(name,title,40,0.0,200.0,-1.0,10.0);          
  AddHistogram(Hist.fDataMetProb_Met_dvp,Folder);
  sprintf(name,"ToyMetProb_toyMet_def");
  sprintf(title," %s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{meas.}>#slash{E}_{T}^{gen}}).vs.#slash{E}_{T}^{toy}, default parametrization");          
  Hist.fToyMetProb_toyMet_def= new TProfile(name,title,40,0.0,200.0,-1.0,10.0);          
  AddHistogram(Hist.fToyMetProb_toyMet_def,Folder);
  sprintf(name,"ToyMetProb_toyMet_dvm");
  sprintf(title," %s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{meas.}>#slash{E}_{T}^{gen}}).vs.#slash{E}_{T}^{toy}, -sigma parametrization");          
  Hist.fToyMetProb_toyMet_dvm= new TProfile(name,title,40,0.0,200.0,-1.0,10.0);          
  AddHistogram(Hist.fToyMetProb_toyMet_dvm,Folder);
  sprintf(name,"ToyMetProb_toyMet_dvp");
  sprintf(title," %s%s",algoname,": -log_{10}(1-P_{#slash{E}_{T}^{meas.}>#slash{E}_{T}^{gen}}).vs.#slash{E}_{T}^{toy}, +sigma parametrization");          
  Hist.fToyMetProb_toyMet_dvp= new TProfile(name,title,40,0.0,200.0,-1.0,10.0);          
  AddHistogram(Hist.fToyMetProb_toyMet_dvp,Folder);

  sprintf(name,"DataMetSigma_Met_def");
  sprintf(title," %s%s",algoname,": (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}.vs.#slash{E}_{T}^{data}, default parametrization");        
  Hist.fDataMetSigma_Met_def= new TProfile(name,title,40,0.0,200.0,-10.0,30.0);          
  AddHistogram(Hist.fDataMetSigma_Met_def,Folder);
  sprintf(name,"DataMetSigma_Met_dvm");
  sprintf(title," %s%s",algoname,": (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}.vs.#slash{E}_{T}^{data}, -sigma parametrization");          
  Hist.fDataMetSigma_Met_dvm= new TProfile(name,title,40,0.0,200.0,-10.0,30.0);          
  AddHistogram(Hist.fDataMetSigma_Met_dvm,Folder);
  sprintf(name,"DataMetSigma_Met_dvp");
  sprintf(title," %s%s",algoname,": (#slash{E}_{T}^{data}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}.vs.#slash{E}_{T}^{data}, +sigma parametrization");          
  Hist.fDataMetSigma_Met_dvp= new TProfile(name,title,40,0.0,200.0,-10.0,30.0);          
  AddHistogram(Hist.fDataMetSigma_Met_dvp,Folder);
  sprintf(name,"ToyMetSigma_toyMet_def");
  sprintf(title," %s%s",algoname,": (#slash{E}_{T}^{meas.}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}.vs.#slash{E}_{T}^{toy}, default parametrization");        
  Hist.fToyMetSigma_toyMet_def= new TProfile(name,title,40,0.0,200.0,-10.0,30.0);          
  AddHistogram(Hist.fToyMetSigma_toyMet_def,Folder);
  sprintf(name,"ToyMetSigma_toyMet_dvm");
  sprintf(title," %s%s",algoname,": (#slash{E}_{T}^{meas.}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}.vs.#slash{E}_{T}^{toy}, -sigma parametrization");          
  Hist.fToyMetSigma_toyMet_dvm= new TProfile(name,title,40,0.0,200.0,-10.0,30.0);          
  AddHistogram(Hist.fToyMetSigma_toyMet_dvm,Folder);
  sprintf(name,"ToyMetSigma_toyMet_dvp");
  sprintf(title," %s%s",algoname,": (#slash{E}_{T}^{meas.}-|#slash{E}_{T}^{gen}|)/#sigma_{#slash{E}_{T}^{gen}}.vs.#slash{E}_{T}^{toy}, +sigma parametrization");          
  Hist.fToyMetSigma_toyMet_dvp= new TProfile(name,title,40,0.0,200.0,-10.0,30.0);          
  AddHistogram(Hist.fToyMetSigma_toyMet_dvp,Folder);

  return;
}
/*}}}*/

//_____________________________________________________________________________
//  CREATE PHOTON-MET STUDY HISTOGRAMS
//_____________________________________________________________________________
void TMyJetFilterModule::BookPhoMetStudyHistograms(PhoMetStudyHisto_t& Hist,
																		const char* Folder)
/*{{{*/
{
  char name [200];
  char title[200];

  sprintf(name,"eta");
  sprintf(title,"|#eta_{det}^{#gamma}|");
  Hist.fPhoMet_eta=new TH1F(name,title,120,0.0,1.2);
  AddHistogram(Hist.fPhoMet_eta,Folder);
  sprintf(name,"xces");
  sprintf(title,"|X_{CES}|");
  Hist.fPhoMet_xces=new TH1F(name,title,86,0.0,23.0);
  AddHistogram(Hist.fPhoMet_xces,Folder);
  sprintf(name,"zces");
  sprintf(title,"|Z_{CES}|");
  Hist.fPhoMet_zces=new TH1F(name,title,240,0.0,240.0);
  AddHistogram(Hist.fPhoMet_zces,Folder);

  sprintf(name,"signeta");
  sprintf(title,"sign(Z_{vx})*#eta_{det}^{#gamma}");
  Hist.fPhoMet_signeta=new TH1F(name,title,240,-1.2,1.2);
  AddHistogram(Hist.fPhoMet_signeta,Folder);
  sprintf(name,"signzces");
  sprintf(title,"sign(Z_{vx})*Z_{CES}");
  Hist.fPhoMet_signzces=new TH1F(name,title,480,-240.0,240.0);
  AddHistogram(Hist.fPhoMet_signzces,Folder);
  sprintf(name,"MetSig5dPhi02_signeta");
  sprintf(title,"sign(Z_{vx})*#eta_{det}^{#gamma} if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_MetSig5dPhi02_signeta=new TH1F(name,title,240,-1.2,1.2);
  AddHistogram(Hist.fPhoMet_MetSig5dPhi02_signeta,Folder);
  sprintf(name,"MetSig5dPhi02_signzces");
  sprintf(title,"sign(Z_{vx})*Z_{CES} if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_MetSig5dPhi02_signzces=new TH1F(name,title,480,-240.0,240.0);
  AddHistogram(Hist.fPhoMet_MetSig5dPhi02_signzces,Folder);
  sprintf(name,"Fr_MetSig5dPhi02_signeta");
  sprintf(title,"Fraction of photons as a function of sign(Z_{vx})*#eta_{det}^{#gamma} if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_Fr_MetSig5dPhi02_signeta=new TH1F(name,title,240,-1.2,1.2);
  AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_signeta,Folder);
  sprintf(name,"Fr_MetSig5dPhi02_signzces");
  sprintf(title,"Fraction of photons as a function of sign(Z_{vx})*Z_{CES} if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_Fr_MetSig5dPhi02_signzces=new TH1F(name,title,480,-240.0,240.0);
  AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_signzces,Folder);

  sprintf(name,"dPhi02_eta");
  sprintf(title,"|#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
  Hist.fPhoMet_dPhi02_eta=new TH1F(name,title,120,0.0,1.2);
  AddHistogram(Hist.fPhoMet_dPhi02_eta,Folder);
  sprintf(name,"dPhi02_xces");
  sprintf(title,"|X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
  Hist.fPhoMet_dPhi02_xces=new TH1F(name,title,86,0.0,23);
  AddHistogram(Hist.fPhoMet_dPhi02_xces,Folder);
  sprintf(name,"dPhi02_zces");
  sprintf(title,"|Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
  Hist.fPhoMet_dPhi02_zces=new TH1F(name,title,240,0.0,240);
  AddHistogram(Hist.fPhoMet_dPhi02_zces,Folder);
  sprintf(name,"dPhi01_eta");
  sprintf(title,"|#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
  Hist.fPhoMet_dPhi01_eta=new TH1F(name,title,120,0.0,1.2);
  AddHistogram(Hist.fPhoMet_dPhi01_eta,Folder);
  sprintf(name,"dPhi01_xces");
  sprintf(title,"|X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
  Hist.fPhoMet_dPhi01_xces=new TH1F(name,title,86,0.0,23);
  AddHistogram(Hist.fPhoMet_dPhi01_xces,Folder);
  sprintf(name,"dPhi01_zces");
  sprintf(title,"|Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
  Hist.fPhoMet_dPhi01_zces=new TH1F(name,title,240,0.0,240);
  AddHistogram(Hist.fPhoMet_dPhi01_zces,Folder);

  sprintf(name,"Fr_dPhi02_eta");
  sprintf(title,"Fraction of photons as a function of |#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
  Hist.fPhoMet_Fr_dPhi02_eta=new TH1F(name,title,120,0.0,1.2);
  AddHistogram(Hist.fPhoMet_Fr_dPhi02_eta,Folder);
  sprintf(name,"Fr_dPhi02_xces");
  sprintf(title,"Fraction of photons as a function of |X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
  Hist.fPhoMet_Fr_dPhi02_xces=new TH1F(name,title,86,0.0,23);
  AddHistogram(Hist.fPhoMet_Fr_dPhi02_xces,Folder);
  sprintf(name,"Fr_dPhi02_zces");
  sprintf(title,"Fraction of photons as a function of |Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
  Hist.fPhoMet_Fr_dPhi02_zces=new TH1F(name,title,240,0.0,240);
  AddHistogram(Hist.fPhoMet_Fr_dPhi02_zces,Folder);
  sprintf(name,"Fr_dPhi01_eta");
  sprintf(title,"Fraction of photons as a function of |#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
  Hist.fPhoMet_Fr_dPhi01_eta=new TH1F(name,title,120,0.0,1.2);
  AddHistogram(Hist.fPhoMet_Fr_dPhi01_eta,Folder);
  sprintf(name,"Fr_dPhi01_xces");
  sprintf(title,"Fraction of photons as a function of |X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
  Hist.fPhoMet_Fr_dPhi01_xces=new TH1F(name,title,86,0.0,23);
  AddHistogram(Hist.fPhoMet_Fr_dPhi01_xces,Folder);
  sprintf(name,"Fr_dPhi01_zces");
  sprintf(title,"Fraction of photons as a function of |Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
  Hist.fPhoMet_Fr_dPhi01_zces=new TH1F(name,title,240,0.0,240);
  AddHistogram(Hist.fPhoMet_Fr_dPhi01_zces,Folder);

  sprintf(name,"MetSig5dPhi02_eta");
  sprintf(title,"|#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_MetSig5dPhi02_eta=new TH1F(name,title,120,0.0,1.2);
  AddHistogram(Hist.fPhoMet_MetSig5dPhi02_eta,Folder);
  sprintf(name,"MetSig5dPhi02_xces");
  sprintf(title,"|X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_MetSig5dPhi02_xces=new TH1F(name,title,86,0.0,23);
  AddHistogram(Hist.fPhoMet_MetSig5dPhi02_xces,Folder);
  sprintf(name,"MetSig5dPhi02_zces");
  sprintf(title,"|Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_MetSig5dPhi02_zces=new TH1F(name,title,240,0.0,240);
  AddHistogram(Hist.fPhoMet_MetSig5dPhi02_zces,Folder);
  sprintf(name,"MetSig5dPhi01_eta");
  sprintf(title,"|#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_MetSig5dPhi01_eta=new TH1F(name,title,120,0.0,1.2);
  AddHistogram(Hist.fPhoMet_MetSig5dPhi01_eta,Folder);
  sprintf(name,"MetSig5dPhi01_xces");
  sprintf(title,"|X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_MetSig5dPhi01_xces=new TH1F(name,title,86,0.0,23);
  AddHistogram(Hist.fPhoMet_MetSig5dPhi01_xces,Folder);
  sprintf(name,"MetSig5dPhi01_zces");
  sprintf(title,"|Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_MetSig5dPhi01_zces=new TH1F(name,title,240,0.0,240);
  AddHistogram(Hist.fPhoMet_MetSig5dPhi01_zces,Folder);

  sprintf(name,"Fr_MetSig5dPhi02_eta");
  sprintf(title,"Fraction of photons as a function of |#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_Fr_MetSig5dPhi02_eta=new TH1F(name,title,120,0.0,1.2);
  AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_eta,Folder);
  sprintf(name,"Fr_MetSig5dPhi02_xces");
  sprintf(title,"Fraction of photons as a function of |X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_Fr_MetSig5dPhi02_xces=new TH1F(name,title,86,0.0,23);
  AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_xces,Folder);
  sprintf(name,"Fr_MetSig5dPhi02_zces");
  sprintf(title,"Fraction of photons as a function of |Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_Fr_MetSig5dPhi02_zces=new TH1F(name,title,240,0.0,240);
  AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_zces,Folder);
  sprintf(name,"Fr_MetSig5dPhi01_eta");
  sprintf(title,"Fraction of photons as a function of |#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_Fr_MetSig5dPhi01_eta=new TH1F(name,title,120,0.0,1.2);
  AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi01_eta,Folder);
  sprintf(name,"Fr_MetSig5dPhi01_xces");
  sprintf(title,"Fraction of photons as a function of |X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_Fr_MetSig5dPhi01_xces=new TH1F(name,title,86,0.0,23);
  AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi01_xces,Folder);
  sprintf(name,"Fr_MetSig5dPhi01_zces");
  sprintf(title,"Fraction of photons as a function of |Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
  Hist.fPhoMet_Fr_MetSig5dPhi01_zces=new TH1F(name,title,240,0.0,240);
  AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi01_zces,Folder);


  //--------------------------- correlation between electrons and met
  sprintf(name,"Ele_eta");
  sprintf(title,"|#eta_{det}^{ele}|");
  Hist.fEleMet_eta=new TH1F(name,title,300,0.0,3.0);
  AddHistogram(Hist.fEleMet_eta,Folder);
  sprintf(name,"Ele_phi");
  sprintf(title,"#phi_{det}");
  Hist.fEleMet_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
  AddHistogram(Hist.fEleMet_phi,Folder);

  sprintf(name,"Ele_dPhi02_eta");
  sprintf(title,"|#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.2");
  Hist.fEleMet_dPhi02_eta=new TH1F(name,title,300,0.0,3.0);
  AddHistogram(Hist.fEleMet_dPhi02_eta,Folder);
  sprintf(name,"Ele_dPhi02_phi");
  sprintf(title,"#phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.2");
  Hist.fEleMet_dPhi02_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
  AddHistogram(Hist.fEleMet_dPhi02_phi,Folder);
  sprintf(name,"Ele_dPhi01_eta");
  sprintf(title,"|#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.1");
  Hist.fEleMet_dPhi01_eta=new TH1F(name,title,300,0.0,3.0);
  AddHistogram(Hist.fEleMet_dPhi01_eta,Folder);
  sprintf(name,"Ele_dPhi01_phi");
  sprintf(title,"#phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.1");
  Hist.fEleMet_dPhi01_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
  AddHistogram(Hist.fEleMet_dPhi01_phi,Folder);

  sprintf(name,"Ele_Fr_dPhi02_eta");
  sprintf(title,"Fraction of electrons as a function of |#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.2");
  Hist.fEleMet_Fr_dPhi02_eta=new TH1F(name,title,300,0.0,3.0);
  AddHistogram(Hist.fEleMet_Fr_dPhi02_eta,Folder);
  sprintf(name,"Ele_Fr_dPhi02_phi");
  sprintf(title,"Fraction of electrons as a function of #phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.2");
  Hist.fEleMet_Fr_dPhi02_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
  AddHistogram(Hist.fEleMet_Fr_dPhi02_phi,Folder);
  sprintf(name,"Ele_Fr_dPhi01_eta");
  sprintf(title,"Fraction of electrons as a function of |#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.1");
  Hist.fEleMet_Fr_dPhi01_eta=new TH1F(name,title,300,0.0,3.0);
  AddHistogram(Hist.fEleMet_Fr_dPhi01_eta,Folder);
  sprintf(name,"Ele_Fr_dPhi01_phi");
  sprintf(title,"Fraction of electrons as a function of #phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.1");
  Hist.fEleMet_Fr_dPhi01_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
  AddHistogram(Hist.fEleMet_Fr_dPhi01_phi,Folder);

  sprintf(name,"Ele_MetSig5dPhi02_eta");
  sprintf(title,"|#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
  Hist.fEleMet_MetSig5dPhi02_eta=new TH1F(name,title,300,0.0,3.0);
  AddHistogram(Hist.fEleMet_MetSig5dPhi02_eta,Folder);
  sprintf(name,"Ele_MetSig5dPhi02_phi");
  sprintf(title,"#phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
  Hist.fEleMet_MetSig5dPhi02_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
  AddHistogram(Hist.fEleMet_MetSig5dPhi02_phi,Folder);
  sprintf(name,"Ele_MetSig5dPhi01_eta");
  sprintf(title,"|#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
  Hist.fEleMet_MetSig5dPhi01_eta=new TH1F(name,title,300,0.0,3.0);
  AddHistogram(Hist.fEleMet_MetSig5dPhi01_eta,Folder);
  sprintf(name,"Ele_MetSig5dPhi01_phi");
  sprintf(title,"#phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
  Hist.fEleMet_MetSig5dPhi01_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
  AddHistogram(Hist.fEleMet_MetSig5dPhi01_phi,Folder);

  sprintf(name,"Ele_Fr_MetSig5dPhi02_eta");
  sprintf(title,"Fraction of electrons as a function of |#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
  Hist.fEleMet_Fr_MetSig5dPhi02_eta=new TH1F(name,title,300,0.0,3.0);
  AddHistogram(Hist.fEleMet_Fr_MetSig5dPhi02_eta,Folder);
  sprintf(name,"Ele_Fr_MetSig5dPhi02_phi");
  sprintf(title,"Fraction of electrons as a function of #phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
  Hist.fEleMet_Fr_MetSig5dPhi02_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
  AddHistogram(Hist.fEleMet_Fr_MetSig5dPhi02_phi,Folder);
  sprintf(name,"Ele_Fr_MetSig5dPhi01_eta");
  sprintf(title,"Fraction of electrons as a function of |#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
  Hist.fEleMet_Fr_MetSig5dPhi01_eta=new TH1F(name,title,300,0.0,3.0);
  AddHistogram(Hist.fEleMet_Fr_MetSig5dPhi01_eta,Folder);
  sprintf(name,"Ele_Fr_MetSig5dPhi01_phi");
  sprintf(title,"Fraction of electrons as a function of #phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
  Hist.fEleMet_Fr_MetSig5dPhi01_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
  AddHistogram(Hist.fEleMet_Fr_MetSig5dPhi01_phi,Folder);

  return;
}
/*}}}*/

//_____________________________________________________________________________
//  FILLING PHOTON-MET STUDY HISTOS      
//_____________________________________________________________________________
void TMyJetFilterModule::FillPhoMetStudyHistograms(PhoMetStudyHisto_t& Hist,
										TVector2 Met, CommonStuff miscstuff, double zvx)
/*{{{*/
{

	for(int i=0; i<miscstuff.myCorrPhoton.size(); i++)
  	{
      double _dphi=fabs(TVector2::Phi_mpi_pi(Met.Phi()-miscstuff.myCorrPhoton[i].Phi()));
      double p0, p1;
      if(fabs(miscstuff.myPhoEtaDet[i])<1.1)
	{
	  p0=0.135*0.135;
	  p1=0.015*0.015;
	}
      else
	{
	  p0=0.16*0.16;
	  p1=0.01*0.01;	  
	}
      double metsig=Met.Mod()*cos(_dphi)/sqrt(miscstuff.myCorrPhoton[i].Pt()*p0+miscstuff.myCorrPhoton[i].Pt()*miscstuff.myCorrPhoton[i].Pt()*p1);
      Hist.fPhoMet_eta->Fill(fabs(miscstuff.myPhoEtaDet[i]));
      Hist.fPhoMet_xces->Fill(fabs(miscstuff.myPhoXces[i]));
      Hist.fPhoMet_zces->Fill(fabs(miscstuff.myPhoZces[i])); 

      if(zvx*miscstuff.myPhoEtaDet[i]>0.0) Hist.fPhoMet_signeta->Fill(fabs(miscstuff.myPhoEtaDet[i]));
      else Hist.fPhoMet_signeta->Fill(-1.0*fabs(miscstuff.myPhoEtaDet[i]));
      if(zvx*miscstuff.myPhoZces[i]>0.0) Hist.fPhoMet_signzces->Fill(fabs(miscstuff.myPhoZces[i]));
      else Hist.fPhoMet_signzces->Fill(-1.0*fabs(miscstuff.myPhoZces[i]));

      if(Met.Mod()>0.0 && _dphi<0.2)
	{
	  Hist.fPhoMet_dPhi02_eta->Fill(fabs(miscstuff.myPhoEtaDet[i]));
	  Hist.fPhoMet_dPhi02_xces->Fill(fabs(miscstuff.myPhoXces[i]));
	  Hist.fPhoMet_dPhi02_zces->Fill(fabs(miscstuff.myPhoZces[i]));
	  if(metsig>5.0)
	    {
	      Hist.fPhoMet_MetSig5dPhi02_eta->Fill(fabs(miscstuff.myPhoEtaDet[i]));
	      Hist.fPhoMet_MetSig5dPhi02_xces->Fill(fabs(miscstuff.myPhoXces[i]));
	      Hist.fPhoMet_MetSig5dPhi02_zces->Fill(fabs(miscstuff.myPhoZces[i]));	      
	      if(zvx*miscstuff.myPhoEtaDet[i]>0.0) Hist.fPhoMet_MetSig5dPhi02_signeta->Fill(fabs(miscstuff.myPhoEtaDet[i]));
	      else Hist.fPhoMet_MetSig5dPhi02_signeta->Fill(-1.0*fabs(miscstuff.myPhoEtaDet[i]));
	      if(zvx*miscstuff.myPhoZces[i]>0.0) Hist.fPhoMet_MetSig5dPhi02_signzces->Fill(fabs(miscstuff.myPhoZces[i]));
	      else Hist.fPhoMet_MetSig5dPhi02_signzces->Fill(-1.0*fabs(miscstuff.myPhoZces[i]));
	    }
	  if(_dphi<0.1)
	    {
	      Hist.fPhoMet_dPhi01_eta->Fill(fabs(miscstuff.myPhoEtaDet[i]));
	      Hist.fPhoMet_dPhi01_xces->Fill(fabs(miscstuff.myPhoXces[i]));
	      Hist.fPhoMet_dPhi01_zces->Fill(fabs(miscstuff.myPhoZces[i]));
	      if(metsig>5.0)
		{
		  Hist.fPhoMet_MetSig5dPhi01_eta->Fill(fabs(miscstuff.myPhoEtaDet[i]));
		  Hist.fPhoMet_MetSig5dPhi01_xces->Fill(fabs(miscstuff.myPhoXces[i]));
		  Hist.fPhoMet_MetSig5dPhi01_zces->Fill(fabs(miscstuff.myPhoZces[i]));
		}
	    }
	}
    }

  for(int i=0; i<miscstuff.myCorrElectron.size(); i++)
    {
      double _dphi=fabs(TVector2::Phi_mpi_pi(Met.Phi()-miscstuff.myCorrElectron[i].Phi()));
      double p0, p1;
      if(fabs(miscstuff.myEleEtaDet[i])<1.1)
	{
	  p0=0.135*0.135;
	  p1=0.015*0.015;
	}
      else
	{
	  p0=0.16*0.16;
	  p1=0.01*0.01;	  
	}
      double metsig=Met.Mod()*cos(_dphi)/sqrt(miscstuff.myCorrElectron[i].Pt()*p0+miscstuff.myCorrElectron[i].Pt()*miscstuff.myCorrElectron[i].Pt()*p1);
      Hist.fEleMet_eta->Fill(fabs(miscstuff.myEleEtaDet[i]));
      Hist.fEleMet_phi->Fill(fabs(miscstuff.myElePhiDet[i]));

      if(Met.Mod()>0.0 && _dphi<0.2)
	{
	  Hist.fEleMet_dPhi02_eta->Fill(fabs(miscstuff.myEleEtaDet[i]));
	  Hist.fEleMet_dPhi02_phi->Fill(fabs(miscstuff.myElePhiDet[i]));
	  if(metsig>5.0)
	    {
	      Hist.fEleMet_MetSig5dPhi02_eta->Fill(fabs(miscstuff.myEleEtaDet[i]));
	      Hist.fEleMet_MetSig5dPhi02_phi->Fill(fabs(miscstuff.myElePhiDet[i]));
	    }
	  if(_dphi<0.1)
	    {
	      Hist.fEleMet_dPhi01_eta->Fill(fabs(miscstuff.myEleEtaDet[i]));
	      Hist.fEleMet_dPhi01_phi->Fill(fabs(miscstuff.myElePhiDet[i]));
	      if(metsig>5.0)
		{
		  Hist.fEleMet_MetSig5dPhi01_eta->Fill(fabs(miscstuff.myEleEtaDet[i]));
		  Hist.fEleMet_MetSig5dPhi01_phi->Fill(fabs(miscstuff.myElePhiDet[i]));
		}
	    }
	}
    }

  return;
}
/*}}}*/


//_____________________________________________________________________________
//  FINALIZE HISTOS FOR PHOTON-MET STUDY  
//_____________________________________________________________________________
void TMyJetFilterModule::DoFinalPhoMetHisto(PhoMetStudyHisto_t& Hist)
/*{{{*/
{

  Hist.fPhoMet_Fr_dPhi02_eta->Sumw2(); 
  Hist.fPhoMet_Fr_dPhi02_xces->Sumw2();
  Hist.fPhoMet_Fr_dPhi02_zces->Sumw2();
  Hist.fPhoMet_Fr_dPhi01_eta->Sumw2(); 
  Hist.fPhoMet_Fr_dPhi01_xces->Sumw2();
  Hist.fPhoMet_Fr_dPhi01_zces->Sumw2();

  Hist.fPhoMet_Fr_MetSig5dPhi02_eta->Sumw2(); 
  Hist.fPhoMet_Fr_MetSig5dPhi02_xces->Sumw2();
  Hist.fPhoMet_Fr_MetSig5dPhi02_zces->Sumw2();
  Hist.fPhoMet_Fr_MetSig5dPhi01_eta->Sumw2(); 
  Hist.fPhoMet_Fr_MetSig5dPhi01_xces->Sumw2();
  Hist.fPhoMet_Fr_MetSig5dPhi01_zces->Sumw2();

  Hist.fPhoMet_Fr_MetSig5dPhi02_signeta->Sumw2(); 
  Hist.fPhoMet_Fr_MetSig5dPhi02_signzces->Sumw2();

  Hist.fPhoMet_Fr_dPhi02_eta->Divide(Hist.fPhoMet_dPhi02_eta,Hist.fPhoMet_eta); 
  Hist.fPhoMet_Fr_dPhi02_xces->Divide(Hist.fPhoMet_dPhi02_xces,Hist.fPhoMet_xces);
  Hist.fPhoMet_Fr_dPhi02_zces->Divide(Hist.fPhoMet_dPhi02_zces,Hist.fPhoMet_zces);
  Hist.fPhoMet_Fr_dPhi01_eta->Divide(Hist.fPhoMet_dPhi01_eta,Hist.fPhoMet_eta); 
  Hist.fPhoMet_Fr_dPhi01_xces->Divide(Hist.fPhoMet_dPhi01_xces,Hist.fPhoMet_xces);
  Hist.fPhoMet_Fr_dPhi01_zces->Divide(Hist.fPhoMet_dPhi01_zces,Hist.fPhoMet_zces);

  Hist.fPhoMet_Fr_MetSig5dPhi02_eta->Divide(Hist.fPhoMet_MetSig5dPhi02_eta,Hist.fPhoMet_eta); 
  Hist.fPhoMet_Fr_MetSig5dPhi02_xces->Divide(Hist.fPhoMet_MetSig5dPhi02_xces,Hist.fPhoMet_xces);
  Hist.fPhoMet_Fr_MetSig5dPhi02_zces->Divide(Hist.fPhoMet_MetSig5dPhi02_zces,Hist.fPhoMet_zces);
  Hist.fPhoMet_Fr_MetSig5dPhi01_eta->Divide(Hist.fPhoMet_MetSig5dPhi01_eta,Hist.fPhoMet_eta); 
  Hist.fPhoMet_Fr_MetSig5dPhi01_xces->Divide(Hist.fPhoMet_MetSig5dPhi01_xces,Hist.fPhoMet_xces);
  Hist.fPhoMet_Fr_MetSig5dPhi01_zces->Divide(Hist.fPhoMet_MetSig5dPhi01_zces,Hist.fPhoMet_zces);

  Hist.fPhoMet_Fr_MetSig5dPhi02_signeta->Divide(Hist.fPhoMet_MetSig5dPhi02_signeta,Hist.fPhoMet_signeta); 
  Hist.fPhoMet_Fr_MetSig5dPhi02_signzces->Divide(Hist.fPhoMet_MetSig5dPhi02_signzces,Hist.fPhoMet_signzces);


  //-------------------- correlation between electron and met ---------------------

  Hist.fEleMet_Fr_dPhi02_eta->Sumw2(); 
  Hist.fEleMet_Fr_dPhi02_phi->Sumw2();
  Hist.fEleMet_Fr_dPhi01_eta->Sumw2(); 
  Hist.fEleMet_Fr_dPhi01_phi->Sumw2();

  Hist.fEleMet_Fr_MetSig5dPhi02_eta->Sumw2(); 
  Hist.fEleMet_Fr_MetSig5dPhi02_phi->Sumw2();
  Hist.fEleMet_Fr_MetSig5dPhi01_eta->Sumw2(); 
  Hist.fEleMet_Fr_MetSig5dPhi01_phi->Sumw2();

  Hist.fEleMet_Fr_dPhi02_eta->Divide(Hist.fEleMet_dPhi02_eta,Hist.fEleMet_eta); 
  Hist.fEleMet_Fr_dPhi02_phi->Divide(Hist.fEleMet_dPhi02_phi,Hist.fEleMet_phi);
  Hist.fEleMet_Fr_dPhi01_eta->Divide(Hist.fEleMet_dPhi01_eta,Hist.fEleMet_eta); 
  Hist.fEleMet_Fr_dPhi01_phi->Divide(Hist.fEleMet_dPhi01_phi,Hist.fEleMet_phi);

  Hist.fEleMet_Fr_MetSig5dPhi02_eta->Divide(Hist.fEleMet_MetSig5dPhi02_eta,Hist.fEleMet_eta); 
  Hist.fEleMet_Fr_MetSig5dPhi02_phi->Divide(Hist.fEleMet_MetSig5dPhi02_phi,Hist.fEleMet_phi);
  Hist.fEleMet_Fr_MetSig5dPhi01_eta->Divide(Hist.fEleMet_MetSig5dPhi01_eta,Hist.fEleMet_eta); 
  Hist.fEleMet_Fr_MetSig5dPhi01_phi->Divide(Hist.fEleMet_MetSig5dPhi01_phi,Hist.fEleMet_phi);

  return;
}
/*}}}*/


//_____________________________________________________________________________
//---------- Filling histograms for Met cleanup studies -----------------------
//_____________________________________________________________________________
void TMyJetFilterModule::FillMetCleanupHistograms(MetCleanupHisto_t& Hist,
						JetStuff jetstuff, CommonStuff miscstuff, int metcode)
/*{{{*/
{
  
  double dphi;
  double ratio;
  double etadet;
  double ptcut=0.0;
  TVector2 metvec(0.0,0.0);
  TVector2 metvec_raw(0.0,0.0);
  
  double dphi_3=10.0;
  double dphi_5=10.0;
  double dphi_10=10.0;
  double dphi_1st=10.0;
  
  if(metcode==1)
    {
      metvec.Set(jetstuff.myMETcorr_th5.Px(),jetstuff.myMETcorr_th5.Py());
      ptcut=5.0;
    }
  if(metcode==2)
    {
      metvec.Set(jetstuff.myMETcorr_th10.Px(),jetstuff.myMETcorr_th10.Py());
      ptcut=10.0;
    }
  if(metcode==3)
    {
      metvec.Set(jetstuff.myMETcorr_th15.Px(),jetstuff.myMETcorr_th15.Py());
      ptcut=15.0;
    }
  if(metcode==4)
    {
      metvec.Set(jetstuff.myMETcorr_th20.Px(),jetstuff.myMETcorr_th20.Py());
      ptcut=20.0;
    }
  if(metcode<1 || metcode>4)
    {
      metvec.Set(miscstuff.myMET_raw.Px(),miscstuff.myMET_raw.Py());
      ptcut=0.0;
    }
  metvec_raw.Set(miscstuff.myMET_raw.Px(),miscstuff.myMET_raw.Py());
  
  if(metvec.Mod()>0.0)
    {
      for(int i=0; i<miscstuff.myRawPhoton.size(); i++)
	{
	  dphi=fabs(TVector2::Phi_mpi_pi(metvec.Phi()-miscstuff.myRawPhoton[i].Phi()));
	  ratio=metvec.Mod()/miscstuff.myRawPhoton[i].Pt();
	  etadet=miscstuff.myPhoEtaDet[i];
	  Hist.fCleanupDelPhiVsEtaDet_em->Fill(etadet,dphi);
	  Hist.fCleanupDelPhi_metem->Fill(dphi);
	  Hist.fCleanupMetEtem->Fill(ratio);
	  Hist.fCleanupMetEtemVsDelPhi->Fill(dphi,ratio);
	}
      for(int i=0; i<miscstuff.myRawElectron.size(); i++)
	{
	  dphi=fabs(TVector2::Phi_mpi_pi(metvec.Phi()-miscstuff.myRawElectron[i].Phi()));
	  ratio=metvec.Mod()/miscstuff.myRawElectron[i].Pt();
	  etadet=miscstuff.myEleEtaDet[i];
	  Hist.fCleanupDelPhiVsEtaDet_em->Fill(etadet,dphi);
	  Hist.fCleanupDelPhi_metem->Fill(dphi);
	  Hist.fCleanupMetEtem->Fill(ratio);
	  Hist.fCleanupMetEtemVsDelPhi->Fill(dphi,ratio);
	}
      
      double maxEt=15.0; // new part as of 04/03/07
      double minDelPhi=10.0;  // new part as of 04/03/07
      double minDelPhi_1st=10.0;  // new part as of 04/03/07
      for(int i=0; i<jetstuff.Jet_lev6_noEMobj.size(); i++)
 	{
	  dphi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(metvec.Phi())-TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj[i].Phi()))); 
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>maxEt && fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<3.0)  // new part as of 04/03/07
	    {
	      maxEt=jetstuff.Jet_lev6_noEMobj[i].Pt();
	      minDelPhi_1st=dphi;
	    }
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>15.0 
	     && fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<3.0 
	     && dphi<minDelPhi) minDelPhi=dphi;  // new part as of 04/03/07
	  
 	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>ptcut)
	    {
	      ratio=metvec.Mod()/jetstuff.Jet_lev6_noEMobj[i].Pt();
	      etadet=jetstuff.EtaDetCorr[i];
	      Hist.fCleanupDelPhiVsEtaDet_jet->Fill(etadet,dphi);	     
	      int j=1;
	      if(fabs(etadet)>2.6) j=2;
	      if(fabs(etadet)<0.2 || (fabs(etadet)>0.9 && fabs(etadet)<1.3)) j=0;	      
	      Hist.fCleanupDelPhi_metjet[j]->Fill(dphi);
	      Hist.fCleanupMetEtjet[j]->Fill(ratio);
	      Hist.fCleanupMetEtjetVsDelPhi[j]->Fill(dphi,ratio);
	    }
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>5.0 && fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<3.0) Hist.fCleanupDelPhi_metjet5->Fill(dphi);  // new part as of 04/03/07
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>15.0 && fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<3.0)   // new part as of 04/03/07
	    {
	      Hist.fCleanupDelPhi_metjet15->Fill(dphi);
	      if(metvec.Mod()>10.0) Hist.fCleanupDelPhi_met10jet15->Fill(dphi);
	    }
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>25 && fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<3.0) Hist.fCleanupDelPhi_metjet25->Fill(dphi);  // new part as of 04/03/07
	}  
      if(minDelPhi_1st<=TMath::Pi())  // new part as of 04/03/07
	{
	  Hist.fCleanupDelPhi_metjet1st15->Fill(minDelPhi_1st);
	  if(metvec.Mod()>10.0) Hist.fCleanupDelPhi_met10jet1st15->Fill(minDelPhi_1st);
	}
      if(minDelPhi<=TMath::Pi())  // new part as of 04/03/07
	{
	  Hist.fCleanupDelPhi_metjet15_dPhiMin->Fill(minDelPhi);
	  if(metvec.Mod()>10.0) Hist.fCleanupDelPhi_met10jet15_dPhiMin->Fill(minDelPhi);
	}
    }  
  
  if(metvec_raw.Mod()>0.0)
    {
      double rawEtmax=0.0;
      for(int i=0; i<jetstuff.Jet_raw_noEMobj.size(); i++)
 	{
	  double _dphi=15.0; // dummy value 
	  if(jetstuff.Jet_raw_noEMobj[i].Pt()>3.0)
	    {
	      _dphi=fabs(TVector2::Phi_mpi_pi(metvec_raw.Phi()-jetstuff.Jet_raw_noEMobj[i].Phi()));
	      if(jetstuff.Jet_raw_noEMobj[i].Pt()>rawEtmax)
		{
		  rawEtmax=jetstuff.Jet_raw_noEMobj[i].Pt();
		  dphi_1st=_dphi;
		}
	      if(_dphi<dphi_3) dphi_3=_dphi;
	      if(jetstuff.Jet_raw_noEMobj[i].Pt()>5.0 && _dphi<dphi_5) dphi_5=_dphi;
	      if(jetstuff.Jet_raw_noEMobj[i].Pt()>10.0 && _dphi<dphi_10) dphi_10=_dphi;
	    }
	}
      if(dphi_3<3.2) Hist.fCleanupDelPhi_RawMetRawJet3->Fill(dphi_3);
      if(dphi_5<3.2) Hist.fCleanupDelPhi_RawMetRawJet5->Fill(dphi_5);
      if(dphi_10<3.2) Hist.fCleanupDelPhi_RawMetRawJet10->Fill(dphi_10);
      if(dphi_1st<3.2) Hist.fCleanupDelPhi_RawMet1stRawJet->Fill(dphi_1st);
    } 
  
  return;
}
/*}}}*/


//_____________________________________________________________________________
//__________ Filling histograms for GenMet cleanup studies --------------------
//_____________________________________________________________________________
void TMyJetFilterModule::FillGenMetCleanupHistograms(MetCleanupHisto_t& Hist,
						JetStuff jetstuff, CommonStuff miscstuff, TVector2 genMet,
						std::vector<TLorentzVector> vec)
/*{{{*/
{
  double dphi;
  double ratio;
  double etadet;
  double ptcut=15.0;
  if(genMet.Mod()>0.0 && genMet.Mod()<2000.0)
    {
      for(int i=0; i<miscstuff.myRawPhoton.size(); i++)
	{
	  dphi=fabs(TVector2::Phi_mpi_pi(genMet.Phi()-miscstuff.myRawPhoton[i].Phi()));
	  ratio=genMet.Mod()/miscstuff.myRawPhoton[i].Pt();
	  etadet=miscstuff.myPhoEtaDet[i];
	  Hist.fCleanupDelPhiVsEtaDet_em->Fill(etadet,dphi);
	  Hist.fCleanupDelPhi_metem->Fill(dphi);
	  Hist.fCleanupMetEtem->Fill(ratio);
	  Hist.fCleanupMetEtemVsDelPhi->Fill(dphi,ratio);
	}
      for(int i=0; i<miscstuff.myRawElectron.size(); i++)
	{
	  dphi=fabs(TVector2::Phi_mpi_pi(genMet.Phi()-miscstuff.myRawElectron[i].Phi()));
	  ratio=genMet.Mod()/miscstuff.myRawElectron[i].Pt();
	  etadet=miscstuff.myEleEtaDet[i];
	  Hist.fCleanupDelPhiVsEtaDet_em->Fill(etadet,dphi);
	  Hist.fCleanupDelPhi_metem->Fill(dphi);
	  Hist.fCleanupMetEtem->Fill(ratio);
	  Hist.fCleanupMetEtemVsDelPhi->Fill(dphi,ratio);
	}

      double maxEt=15.0; // new part as of 04/03/07
      double minDelPhi=10.0;  // new part as of 04/03/07
      double minDelPhi_1st=10.0;  // new part as of 04/03/07
      for(int i=0; i<vec.size(); i++)
 	{
	  dphi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(genMet.Phi())-TVector2::Phi_0_2pi(vec[i].Phi())));
	  if(vec[i].Pt()>0.0) ratio=genMet.Mod()/vec[i].Pt();
	  else ratio=9999.0; // dummy value
	  etadet=jetstuff.EtaDetCorr[i];
	  if(vec[i].Pt()>maxEt && fabs(vec[i].Eta())<3.0)  // new part as of 04/03/07
	    {
	      maxEt=vec[i].Pt();
	      minDelPhi_1st=dphi;
	    }
	  if(vec[i].Pt()>15.0 
	     && fabs(vec[i].Eta())<3.0 
	     && dphi<minDelPhi) minDelPhi=dphi;  // new part as of 04/03/07

	  if(vec[i].Pt()>ptcut)
	    {
	      Hist.fCleanupDelPhiVsEtaDet_jet->Fill(etadet,dphi);	     
	      int j=1;
	      if(fabs(etadet)>2.6) j=2;
	      if(fabs(etadet)<0.2 || (fabs(etadet)>0.9 && fabs(etadet)<1.3)) j=0;	      
	      Hist.fCleanupDelPhi_metjet[j]->Fill(dphi);
	      Hist.fCleanupMetEtjet[j]->Fill(ratio);
	      Hist.fCleanupMetEtjetVsDelPhi[j]->Fill(dphi,ratio);
	    }
	  if(vec[i].Pt()>5.0 && fabs(vec[i].Eta())<3.0) Hist.fCleanupDelPhi_metjet5->Fill(dphi);  // new part as of 04/03/07
	  if(vec[i].Pt()>15.0 && fabs(vec[i].Eta())<3.0)   // new part as of 04/03/07
	    {
	      Hist.fCleanupDelPhi_metjet15->Fill(dphi);
	      if(genMet.Mod()>10.0) Hist.fCleanupDelPhi_met10jet15->Fill(dphi);
	    }
	  if(vec[i].Pt()>25 && fabs(vec[i].Eta())<3.0) Hist.fCleanupDelPhi_metjet25->Fill(dphi);  // new part as of 04/03/07
	}
      if(minDelPhi_1st<=TMath::Pi())  // new part as of 04/03/07
	{
	  Hist.fCleanupDelPhi_metjet1st15->Fill(minDelPhi_1st);
	  if(genMet.Mod()>10.0) Hist.fCleanupDelPhi_met10jet1st15->Fill(minDelPhi_1st);
	}
      if(minDelPhi<=TMath::Pi())  // new part as of 04/03/07
	{
	  Hist.fCleanupDelPhi_metjet15_dPhiMin->Fill(minDelPhi);
	  if(genMet.Mod()>10.0) Hist.fCleanupDelPhi_met10jet15_dPhiMin->Fill(minDelPhi);
	}
    }  
  return;
}
/*}}}*/


//_____________________________________________________________________________
//__________ filling met histo _________________________________________________
//_____________________________________________________________________________
void TMyJetFilterModule::FillMetStudyHistograms(MetStudyHisto_t& Hist,
						JetStuff jetstuff, CommonStuff miscstuff, int metcode,
						int nvx)
/*{{{*/
{
  
  double mgg=0.0;
  if(miscstuff.myCorrPhoton.size()>1) mgg=(miscstuff.myCorrPhoton[0]+miscstuff.myCorrPhoton[1]).M();
  double ht=MyHtAll(jetstuff,miscstuff,metcode);

  double sqrtSumEt;
  double SumEt;
  double sqrtSumEt_jet;
  double SumEt_jet;
  double JetFr=0.0;
  double met;
  double metR=miscstuff.myMET_raw.Mod();
  double metR0=miscstuff.myMET0_raw.Mod();
  double metx;
  double mety;
  double metRx=miscstuff.myMET_raw.Px();
  double metRy=miscstuff.myMET_raw.Py();
  double metphiR=miscstuff.myMET_raw.Phi();
  double metphiC;
  int njet=0;

  if(metcode==1)
    {
      SumEt=jetstuff.mySumEtCorr_th5;
      sqrtSumEt=sqrt(jetstuff.mySumEtCorr_th5);
      SumEt_jet=jetstuff.mySumEtJet_th5;
      sqrtSumEt_jet=sqrt(jetstuff.mySumEtJet_th5);
      met=jetstuff.myMETcorr_th5.Mod();
      metx=jetstuff.myMETcorr_th5.Px();
      mety=jetstuff.myMETcorr_th5.Py();
      metphiC=jetstuff.myMETcorr_th5.Phi();
      njet=jetstuff.myNjet_th5;
    }
  if(metcode==2)
    {
      SumEt=jetstuff.mySumEtCorr_th10;
      sqrtSumEt=sqrt(jetstuff.mySumEtCorr_th10);
      SumEt_jet=jetstuff.mySumEtJet_th10;
      sqrtSumEt_jet=sqrt(jetstuff.mySumEtJet_th10);
      met=jetstuff.myMETcorr_th10.Mod();
      metx=jetstuff.myMETcorr_th10.Px();
      mety=jetstuff.myMETcorr_th10.Py();
      metphiC=jetstuff.myMETcorr_th10.Phi();
      njet=jetstuff.myNjet_th10;
    }
  if(metcode==3)
    {
      SumEt=jetstuff.mySumEtCorr_th15;
      sqrtSumEt=sqrt(jetstuff.mySumEtCorr_th15);
      SumEt_jet=jetstuff.mySumEtJet_th15;
      sqrtSumEt_jet=sqrt(jetstuff.mySumEtJet_th15);
      met=jetstuff.myMETcorr_th15.Mod();
      metx=jetstuff.myMETcorr_th15.Px();
      mety=jetstuff.myMETcorr_th15.Py();
      metphiC=jetstuff.myMETcorr_th15.Phi();
      njet=jetstuff.myNjet_th15;
    }
  if(metcode==4)
    {
      SumEt=jetstuff.mySumEtCorr_th20;
      sqrtSumEt=sqrt(jetstuff.mySumEtCorr_th20);
      SumEt_jet=jetstuff.mySumEtJet_th20;
      sqrtSumEt_jet=sqrt(jetstuff.mySumEtJet_th20);
      met=jetstuff.myMETcorr_th20.Mod();
      metx=jetstuff.myMETcorr_th20.Px();
      mety=jetstuff.myMETcorr_th20.Py();
      metphiC=jetstuff.myMETcorr_th20.Phi();
      njet=jetstuff.myNjet_th20;
    }

  if(metcode<1 || metcode>4)
    {
      SumEt=miscstuff.mySumEt_raw;
      sqrtSumEt=sqrt(miscstuff.mySumEt_raw);
      SumEt_jet=0.0; // not defined for this case
      sqrtSumEt_jet=0.0; // not defined for this case
      met=miscstuff.myMET_raw.Mod();
      metx=miscstuff.myMET_raw.Px();
      mety=miscstuff.myMET_raw.Py();
      metphiC=miscstuff.myMET_raw.Phi();
      njet=jetstuff.myNjet;
    }
  if(SumEt>0.0) JetFr=SumEt_jet/SumEt;

  if(fabs(met-metR)>0.0)
    {
      Hist.fMetRaw->Fill(metR); 
      Hist.fMetCorr->Fill(met);
      Hist.fMetPhiRaw->Fill(metphiR); 
      Hist.fMetPhiCorr->Fill(metphiC);
    }

  Hist.fSumEtCorr->Fill(SumEt);    
  Hist.fSqrtSumEtCorr->Fill(sqrtSumEt);

  if(njet>0)
    {
      Hist.fSumEtCorr_withJet->Fill(SumEt);    
      Hist.fSqrtSumEtCorr_withJet->Fill(sqrtSumEt);
    }
  else
    {
      Hist.fSumEtCorr_noJet->Fill(SumEt);    
      Hist.fSqrtSumEtCorr_noJet->Fill(sqrtSumEt);   
      if(nvx==1) Hist.fSumEtCorr_noJet_vx1->Fill(SumEt);     
      if(nvx==2) Hist.fSumEtCorr_noJet_vx2->Fill(SumEt);     
      if(nvx==3) Hist.fSumEtCorr_noJet_vx3->Fill(SumEt);
      Hist.fSumEtCorrNoJet_vs_Nvx->Fill(nvx,SumEt);
    }
 
  Hist.fMetRawAll->Fill(metR); 
  Hist.fMetCorrAll->Fill(met);
  Hist.fMetRawAll_X->Fill(metRx); 
  Hist.fMetCorrAll_X->Fill(metx);
  Hist.fMetRawAll_Y->Fill(metRy); 
  Hist.fMetCorrAll_Y->Fill(mety);
  Hist.fMetPhiRawAll->Fill(metphiR); 
  Hist.fMetPhiCorrAll->Fill(metphiC);
    
  Hist.fMet4VsSqrtCorrSumEt->Fill(sqrtSumEt,met);
  Hist.fMet4VsCorrSumEt->Fill(SumEt,met);

  if(nvx==1) Hist.fMet4VsSqrtCorrSumEt_nvx1->Fill(sqrtSumEt,met);  
  if(nvx==2) Hist.fMet4VsSqrtCorrSumEt_nvx2->Fill(sqrtSumEt,met);  
  if(nvx==3) Hist.fMet4VsSqrtCorrSumEt_nvx3->Fill(sqrtSumEt,met);  
  if(nvx==4) Hist.fMet4VsSqrtCorrSumEt_nvx4->Fill(sqrtSumEt,met);  
  if(nvx==5) Hist.fMet4VsSqrtCorrSumEt_nvx5->Fill(sqrtSumEt,met);  

  //---- Below this point, all histograms for met model will be filled using the last generated value of met 

  Hist.fMet4VsNjet->Fill(jetstuff.myNjet,met);
  Hist.fMet4VsNjet_th5->Fill(jetstuff.myNjet_th5,met);
  Hist.fMet4VsNjet_th10->Fill(jetstuff.myNjet_th10,met);
  Hist.fMet4VsNjet_th15->Fill(jetstuff.myNjet_th15,met);
  Hist.fMet4VsNjet_th20->Fill(jetstuff.myNjet_th20,met);
  Hist.fMet4VsNjet_th25->Fill(jetstuff.myNjet_th25,met);

  Hist.fMet4VsNvx12->Fill(nvx,met);
  Hist.fMet4X_vs_Met4Y->Fill(metx,mety);

  if(sqrtSumEt>0.0)
    {
      Hist.fMet4SigVsNjet->Fill(jetstuff.myNjet,met/sqrtSumEt);
      Hist.fMet4SigVsNjet_th5->Fill(jetstuff.myNjet_th5,met/sqrtSumEt);
      Hist.fMet4SigVsNjet_th10->Fill(jetstuff.myNjet_th10,met/sqrtSumEt);
      Hist.fMet4SigVsNjet_th15->Fill(jetstuff.myNjet_th15,met/sqrtSumEt);
      Hist.fMet4SigVsNjet_th20->Fill(jetstuff.myNjet_th20,met/sqrtSumEt);
      Hist.fMet4SigVsNjet_th25->Fill(jetstuff.myNjet_th25,met/sqrtSumEt);
      Hist.fMet4SigVsNvx12->Fill(nvx,met/sqrtSumEt);
    }

  if(nvx==1)
    {
      Hist.fMet4VsRun->Fill(GetHeaderBlock()->RunNumber(),met);
      if(sqrtSumEt>0.0) Hist.fMet4SigVsRun->Fill(GetHeaderBlock()->RunNumber(),met/sqrtSumEt); 
    }
  Hist.fNvx12VsRun->Fill(GetHeaderBlock()->RunNumber(),nvx);

  if(SumEt_jet>0.0)
    {
      Hist.fSumEtJetFrac->Fill(JetFr);
      Hist.fSumEtJet->Fill(SumEt_jet);
      Hist.fSqrtSumEtJet->Fill(sqrtSumEt_jet);
    }
  Hist.fMet4VsJetFr->Fill(JetFr,met);
  Hist.fMet4VsSqrtSumEtJet->Fill(sqrtSumEt_jet,met);
  if(sqrtSumEt>0.0)
    {
      Hist.fMet4SigVsJetFr->Fill(JetFr,met/sqrtSumEt); 
      Hist.fMet4SigVsSqrtSumEtJet->Fill(sqrtSumEt_jet,met/sqrtSumEt);
    }
  Hist.fJetFrVsSumEt->Fill(SumEt,JetFr); 


  for(int i=0; i<10; i++)
    {
      if(sqrtSumEt>i*2.0 && sqrtSumEt<=(i+1)*2.0)
	{

	  Hist.fMet4[i]->Fill(met);
	  Hist.fMet4Phi[i]->Fill(metphiC);

	  if(njet==0)
	    {
	      Hist.fMet4X_noJet[i]->Fill(metx);
	      Hist.fMet4Y_noJet[i]->Fill(mety);
	    }
	  else
	    {
	      Hist.fMet4X_withJet[i]->Fill(metx);
	      Hist.fMet4Y_withJet[i]->Fill(mety);
	    }
	  
	  Hist.fMet4X[i]->Fill(metx);
	  Hist.fMet4Y[i]->Fill(mety);
	  if(nvx>1)
	    {
	      Hist.fMet4X_nvxM1[i]->Fill(metx);
	      Hist.fMet4Y_nvxM1[i]->Fill(mety);
	    }
	  if(nvx==1)
	    {
	      Hist.fMet4X_nvx1[i]->Fill(metx);
	      Hist.fMet4Y_nvx1[i]->Fill(mety);
	      break;
	    }
	  if(nvx==2)
	    {
	      Hist.fMet4X_nvx2[i]->Fill(metx);
	      Hist.fMet4Y_nvx2[i]->Fill(mety);
	      break;
	    }
	  if(nvx==3)
	    {
	      Hist.fMet4X_nvx3[i]->Fill(metx);
	      Hist.fMet4Y_nvx3[i]->Fill(mety);
	      break;
	    }
	  if(nvx==4)
	    {
	      Hist.fMet4X_nvx4[i]->Fill(metx);
	      Hist.fMet4Y_nvx4[i]->Fill(mety);
	      break;
	    }
	  if(nvx==5)
	    {
	      Hist.fMet4X_nvx5[i]->Fill(metx);
	      Hist.fMet4Y_nvx5[i]->Fill(mety);
	      break;
	    }
	  break; // just in case if nvx>5
	}
    }

  //__________________________ filling histograms for Met Model-2
  TVector2 myGenV2Met_def(0.0,0.0);
  TVector2 myGenV2Met_devm(0.0,0.0);
  TVector2 myGenV2Met_devp(0.0,0.0);
  TVector2 myGenV2Met_devmUn(0.0,0.0);
  TVector2 myGenV2Met_devpUn(0.0,0.0);

  for(int k=0; k<Npoints; k++)
    {
      //------------------------ default parametrization

      //______________ generating no-clustered component

      //______ this part is new; added on 11/17/05 to improve dPhi=Phi(met)-Phi(jet) modeling
      GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,0,0,0,metcode,myGenV2Met_def);
//       FillGenMetCleanupHistograms(fGenCleanup,jet04stuff,allstuff,myGenV2Met_def); // histo for generated met cleanup studies		 
      //__________________ end of new part _______________________________________________________
      if(myGenV2Met_def.Mod()<2000.0) // new as of 04/04/07
	{
	  
	  Hist.fMet4GenV2Met->Fill(met-myGenV2Met_def.Mod()); // Met4-MetGen
	  Hist.fMet4_vs_Met4GenV2Met->Fill(met-myGenV2Met_def.Mod(),met); // Met4 .vs. Met4-MetGen
	  Hist.fMet4_vs_GenV2Met->Fill(myGenV2Met_def.Mod(),met); // Met4 .vs. MetGen
	  
	  Hist.fGenV2Met_def->Fill(myGenV2Met_def.Mod());    
	  Hist.fGenV2MetX_def->Fill(myGenV2Met_def.Px());   
	  Hist.fGenV2MetY_def->Fill(myGenV2Met_def.Py());   
	  Hist.fGenV2MetPhi_def->Fill(myGenV2Met_def.Phi()); 
	  
	  //---- profile histograms
	  Hist.fMetGenV2VsNjet->Fill(jetstuff.myNjet,myGenV2Met_def.Mod()); 
	  Hist.fMetGenV2VsNjet_th5->Fill(jetstuff.myNjet_th5,myGenV2Met_def.Mod());
	  Hist.fMetGenV2VsNjet_th10->Fill(jetstuff.myNjet_th10,myGenV2Met_def.Mod());
	  Hist.fMetGenV2VsNjet_th15->Fill(jetstuff.myNjet_th15,myGenV2Met_def.Mod());
	  Hist.fMetGenV2VsNjet_th20->Fill(jetstuff.myNjet_th20,myGenV2Met_def.Mod());
	  Hist.fMetGenV2VsNjet_th25->Fill(jetstuff.myNjet_th25,myGenV2Met_def.Mod());
	  Hist.fMetGenV2VsNvx12->Fill(nvx,myGenV2Met_def.Mod());
	  Hist.fGenV2MetX_vs_MetY->Fill(myGenV2Met_def.Px(),myGenV2Met_def.Py());
	  if(sqrtSumEt>0.0)
	    {
	      Hist.fMetGenV2SigVsNjet->Fill(jetstuff.myNjet,myGenV2Met_def.Mod()/sqrtSumEt); 
	      Hist.fMetGenV2SigVsNjet_th5->Fill(jetstuff.myNjet_th5,myGenV2Met_def.Mod()/sqrtSumEt);
	      Hist.fMetGenV2SigVsNjet_th10->Fill(jetstuff.myNjet_th10,myGenV2Met_def.Mod()/sqrtSumEt);
	      Hist.fMetGenV2SigVsNjet_th15->Fill(jetstuff.myNjet_th15,myGenV2Met_def.Mod()/sqrtSumEt);
	      Hist.fMetGenV2SigVsNjet_th20->Fill(jetstuff.myNjet_th20,myGenV2Met_def.Mod()/sqrtSumEt);
	      Hist.fMetGenV2SigVsNjet_th25->Fill(jetstuff.myNjet_th25,myGenV2Met_def.Mod()/sqrtSumEt);
	      Hist.fMetGenV2SigVsNvx12->Fill(nvx,myGenV2Met_def.Mod()/sqrtSumEt);
	    }
	  if(nvx==1)
	    {
	      Hist.fMetGenV2VsRun->Fill(GetHeaderBlock()->RunNumber(),myGenV2Met_def.Mod());
	      if(sqrtSumEt>0.0) Hist.fMetGenV2SigVsRun->Fill(GetHeaderBlock()->RunNumber(),myGenV2Met_def.Mod()/sqrtSumEt); 
	    }
	  Hist.fMetGenV2VsJetFr->Fill(JetFr,myGenV2Met_def.Mod());
	  Hist.fMetGenV2VsSqrtSumEtJet->Fill(sqrtSumEt_jet,myGenV2Met_def.Mod());
	  if(sqrtSumEt>0.0)
	    {
	      Hist.fMetGenV2SigVsJetFr->Fill(JetFr,myGenV2Met_def.Mod()/sqrtSumEt);
	      Hist.fMetGenV2SigVsSqrtSumEtJet->Fill(sqrtSumEt_jet,myGenV2Met_def.Mod()/sqrtSumEt);
	    }
	}
      //------------------------ "-sigma" deviated parametrization
      //______ this part is new; added on 11/17/05 to improve dPhi=Phi(met)-Phi(jet) modeling
      GenerateMyTotalMet(jetstuff,miscstuff,fGenDevmCleanup,0,-1,0,metcode,myGenV2Met_devm);
      //__________________ end of new part _______________________________________________________
      if(myGenV2Met_devm.Mod()<2000.0) // new as of 04/04/07
	{	  
	  Hist.fGenV2Met_devm->Fill(myGenV2Met_devm.Mod());    
	  Hist.fGenV2MetX_devm->Fill(myGenV2Met_devm.Px());   
	  Hist.fGenV2MetY_devm->Fill(myGenV2Met_devm.Py());   
	  Hist.fGenV2MetPhi_devm->Fill(myGenV2Met_devm.Phi()); 
	}
      //------------------------ "+sigma" deviated parametrization
      //______ this part is new; added on 11/17/05 to improve dPhi=Phi(met)-Phi(jet) modeling
      GenerateMyTotalMet(jetstuff,miscstuff,fGenDevpCleanup,0,1,0,metcode,myGenV2Met_devp);
      //__________________ end of new part _______________________________________________________
      if(myGenV2Met_devp.Mod()<2000.0) // new as of 04/04/07
	{
	  Hist.fGenV2Met_devp->Fill(myGenV2Met_devp.Mod());    
	  Hist.fGenV2MetX_devp->Fill(myGenV2Met_devp.Px());   
	  Hist.fGenV2MetY_devp->Fill(myGenV2Met_devp.Py());   
	  Hist.fGenV2MetPhi_devp->Fill(myGenV2Met_devp.Phi()); 
	}

      //------------------------ "-sigmaUncl" deviated parametrization
      //______ this part is new; added on 11/17/05 to improve dPhi=Phi(met)-Phi(jet) modeling
      GenerateMyTotalMet(jetstuff,miscstuff,fGenDevmUnclCleanup,0,0,-1,metcode,myGenV2Met_devmUn);
      //__________________ end of new part _______________________________________________________
      if(myGenV2Met_devmUn.Mod()<2000.0)  // new as of 04/04/07
	{
	  Hist.fGenV2Met_devmUn->Fill(myGenV2Met_devmUn.Mod());    
	  Hist.fGenV2MetX_devmUn->Fill(myGenV2Met_devmUn.Px());   
	  Hist.fGenV2MetY_devmUn->Fill(myGenV2Met_devmUn.Py());   
	  Hist.fGenV2MetPhi_devmUn->Fill(myGenV2Met_devmUn.Phi()); 
	}
      //------------------------ "+sigmaUncl" deviated parametrization
      //______ this part is new; added on 11/17/05 to improve dPhi=Phi(met)-Phi(jet) modeling
      GenerateMyTotalMet(jetstuff,miscstuff,fGenDevpUnclCleanup,0,0,1,metcode,myGenV2Met_devpUn);
      //__________________ end of new part _______________________________________________________
      if(myGenV2Met_devpUn.Mod()<2000.0)  // new as of 04/04/07
	{
	  Hist.fGenV2Met_devpUn->Fill(myGenV2Met_devpUn.Mod());    
	  Hist.fGenV2MetX_devpUn->Fill(myGenV2Met_devpUn.Px());   
	  Hist.fGenV2MetY_devpUn->Fill(myGenV2Met_devpUn.Py());   
	  Hist.fGenV2MetPhi_devpUn->Fill(myGenV2Met_devpUn.Phi()); 
	}

      //______ Warning!!! this has to be changed in future to be independent of "if(metcode==3)" condition
      if(metcode==3) 
	{
	  MyMetModelEventCount(jetstuff,myGenV2Met_def,myGenV2Met_devm,myGenV2Met_devp,myGenV2Met_devmUn,myGenV2Met_devpUn,met_results); // count generated events with Met>cut
	  MyMetPDF(k,jetstuff.myMETcorr_th15,myGenV2Met_def,myGenV2Met_devm,myGenV2Met_devp,fMetProbability,SigMetEvent_status); // met probability stuff
	}
    }
  return;
}

/*}}}*/


//_____________________________________________________________________________
//_________ filling general histo for particular Jet Cone _______________________
//_____________________________________________________________________________
void TMyJetFilterModule::FillJetHistogramsB(JetGeneral_t& Hist, JetStuff stuff,
						double dz, int nvx12)
/*{{{*/
{ 


  double toyMet=MetToy_min+(MetToy_max-MetToy_min)*(gRandom->Rndm());
  double toyMetPhi=TMath::TwoPi()*(gRandom->Rndm());
  double toyMet_x=toyMet*cos(toyMetPhi);
  double toyMet_y=toyMet*sin(toyMetPhi);
  TVector2 MetToy(toyMet_x,toyMet_y);
  TVector2 MetToyEvnt(stuff.myMETcorr_th15.Px()+toyMet_x,stuff.myMETcorr_th15.Py()+toyMet_y);
  if(MyMetCleanUpCut(stuff,allstuff,stuff.Jet_lev6_noEMobj,MetToyEvnt)==0) Hist.fEvnt_toyMET_cut->Fill(MetToy.Mod());
  Hist.fEvnt_toyMET_all->Fill(MetToy.Mod());
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
  Hist.fEvntdZ_b->Fill(dz);

//________________ Jet Resolution study part ___________________________ 
  if(stuff.Jet_lev6_noEMobj.size()>0)
    {
      if((stuff.Npho_match[0]+stuff.Nele_match[0])==0) 
	{
	  if(stuff.Jet_lev6_noEMobj[0].Pt()>5.0)  Hist.fJetRes_DetEta_j5_b->Fill(fabs(stuff.EtaDet[0]));
	  if(stuff.Jet_lev6_noEMobj[0].Pt()>10.0) Hist.fJetRes_DetEta_j10_b->Fill(fabs(stuff.EtaDet[0]));
	  if(stuff.Jet_lev6_noEMobj[0].Pt()>15.0) Hist.fJetRes_DetEta_j15_b->Fill(fabs(stuff.EtaDet[0]));

	  float _dphi1=fabs(TVector2::Phi_mpi_pi(stuff.myMETcorr_th15.Phi()-stuff.Jet_lev6_noEMobj[0].Phi()));
	  float _dphi2;
	  float dPhi1=100.0;
	  float dPhi2=100.0;
	  if(_dphi1<TMath::Pi()/2.0) dPhi1=_dphi1;
	  else dPhi1=TMath::Pi()-_dphi1;
	  if(stuff.Jet_lev6_noEMobj.size()>1)
	    {
	      if((stuff.Npho_match[1]+stuff.Nele_match[1])==0 && stuff.Jet_lev6_noEMobj[1].Pt()>5.0) 
		{
		  _dphi2=fabs(TVector2::Phi_mpi_pi(stuff.myMETcorr_th15.Phi()-stuff.Jet_lev6_noEMobj[1].Phi()));
		  if(_dphi1<TMath::Pi()/2.0) dPhi1=_dphi1;
		  else dPhi2=TMath::Pi()-_dphi2;
		}
	    }
	  if(dPhi1<0.4 && dPhi2>dPhi1)
	    {
	      double sigma=stuff.Jet_lev6_noEMobj[0].Pt()*GetMyJetResolution(stuff.Jet_lev6_noEMobj[0].Pt(),stuff.EtaDet[0],0,0);
	      double MetSig=-1.0;
	      if(sigma>0.0) MetSig=stuff.myMETcorr_th15.Mod()*cos(dPhi1)/sigma;
	      if(MetSig>=0.0)
		{
		  if(stuff.Jet_lev6_noEMobj[0].Pt()>5.0)  Hist.fJetRes_MetSig_EtaDet_j5_b->Fill(fabs(stuff.EtaDet[0]),MetSig);
		  if(stuff.Jet_lev6_noEMobj[0].Pt()>10.0) Hist.fJetRes_MetSig_EtaDet_j10_b->Fill(fabs(stuff.EtaDet[0]),MetSig);
		  if(stuff.Jet_lev6_noEMobj[0].Pt()>15.0) Hist.fJetRes_MetSig_EtaDet_j15_b->Fill(fabs(stuff.EtaDet[0]),MetSig);
		  if(MetSig>=3.0)
		    {
		      if(stuff.Jet_lev6_noEMobj[0].Pt()>5.0)  Hist.fJetRes_DetEta_j5_Sig3dPhi04_b->Fill(fabs(stuff.EtaDet[0]));
		      if(stuff.Jet_lev6_noEMobj[0].Pt()>10.0) Hist.fJetRes_DetEta_j10_Sig3dPhi04_b->Fill(fabs(stuff.EtaDet[0]));
		      if(stuff.Jet_lev6_noEMobj[0].Pt()>15.0) Hist.fJetRes_DetEta_j15_Sig3dPhi04_b->Fill(fabs(stuff.EtaDet[0]));
		      if(dPhi1<0.1)
			{
			  if(stuff.Jet_lev6_noEMobj[0].Pt()>5.0)  Hist.fJetRes_DetEta_j5_Sig3dPhi01_b->Fill(fabs(stuff.EtaDet[0]));
			  if(stuff.Jet_lev6_noEMobj[0].Pt()>10.0) Hist.fJetRes_DetEta_j10_Sig3dPhi01_b->Fill(fabs(stuff.EtaDet[0]));
			  if(stuff.Jet_lev6_noEMobj[0].Pt()>15.0) Hist.fJetRes_DetEta_j15_Sig3dPhi01_b->Fill(fabs(stuff.EtaDet[0])); 
			}
		    }
		}
	    }
	}
    }

  for(int i=0; i<stuff.Jet_raw_noEMobj.size(); i++)
    { 
      if((stuff.Npho_match[i]+stuff.Nele_match[i])==0) 
	{
	  jetsum=jetsum+stuff.Jet_lev6_noEMobj[i];
	  ave_et=ave_et+stuff.Jet_lev6_noEMobj[i].Pt();
	  ave_et0=ave_et0+stuff.Jet_raw_noEMobj[i].Pt();
	  Hist.fEvntEt0_Njet_b->Fill(stuff.myNjet,stuff.Jet_raw_noEMobj[i].Pt());
	  Hist.fEvntEt_Njet_b->Fill(stuff.myNjet,stuff.Jet_lev6_noEMobj[i].Pt());
	  if(i<2)
	    {
	      Hist.fEvntEt0_b[i]->Fill(stuff.Jet_raw_noEMobj[i].Pt()); 
	      Hist.fEvntEt_b[i]->Fill(stuff.Jet_lev6_noEMobj[i].Pt()); 
	      Hist.fEvntEtaDet_b[i]->Fill(stuff.EtaDet[i]);
	      Hist.fEvntEta_b[i]->Fill(stuff.Jet_lev6_noEMobj[i].Eta());
	      Hist.fEvntPhi_b[i]->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj[i].Phi()));
	      Hist.fEvntTheta_b[i]->Fill(stuff.Jet_lev6_noEMobj[i].Theta());
	      Hist.fEvntEmFr_b[i]->Fill(stuff.EmFrRaw[i]);
	      Hist.fEvntNTowers_b[i]->Fill(stuff.JetNtwr[i]);
	      Hist.fEvntNTracks_b[i]->Fill(stuff.JetNtrk[i]); 
	    }
	  else
	    {
	      Hist.fEvntEt0X_b->Fill(stuff.Jet_raw_noEMobj[i].Pt());
	      Hist.fEvntEtX_b->Fill(stuff.Jet_lev6_noEMobj[i].Pt());        
	      Hist.fEvntEtaDetX_b->Fill(stuff.EtaDet[i]);    
	      Hist.fEvntEtaX_b->Fill(stuff.Jet_lev6_noEMobj[i].Eta());       
	      Hist.fEvntPhiX_b->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj[i].Phi()));       
	      Hist.fEvntThetaX_b->Fill(stuff.Jet_lev6_noEMobj[i].Theta());     
	      Hist.fEvntEmFrX_b->Fill(stuff.EmFrRaw[i]);      
	      Hist.fEvntNTowersX_b->Fill(stuff.JetNtwr[i]);   
	      Hist.fEvntNTracksX_b->Fill(stuff.JetNtrk[i]);
	      Hist.fEvntDeltaR1x_b->Fill(stuff.Jet_lev6_noEMobj[i].DeltaR(stuff.Jet_lev6[0]));
	      Hist.fEvntDeltaR2x_b->Fill(stuff.Jet_lev6_noEMobj[i].DeltaR(stuff.Jet_lev6[1]));
	    }   
	}
    }
  if(stuff.Jet_raw_noEMobj.size()>1 
     && (stuff.Npho_match[0]+stuff.Nele_match[0])==0 
     && (stuff.Npho_match[1]+stuff.Nele_match[1])==0)
    {
      Hist.fEvntThetaStar_b->Fill(GetThetaStar(stuff.Jet_lev6_noEMobj));
      Hist.fEvntDeltaPhi_b->Fill(fabs(TVector2::Phi_mpi_pi(stuff.Jet_lev6_noEMobj[0].DeltaPhi(stuff.Jet_lev6_noEMobj[1]))));
      Hist.fEvntDeltaEta_b->Fill(stuff.Jet_lev6_noEMobj[0].Eta()-stuff.Jet_lev6_noEMobj[1].Eta());
      Hist.fEvntDeltaR_b->Fill(stuff.Jet_lev6_noEMobj[1].DeltaR(stuff.Jet_lev6_noEMobj[0])); 
      Hist.fEvntMjj_b->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M()); 
      Hist.fEvntKt2jet_b->Fill(GetMyKtKick(stuff.Jet_lev6_noEMobj));  
      Hist.fEvntKtAll_b->Fill(jetsum.Pt());
      Hist.fEvntKt_Mjj_b->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),GetMyKtKick(stuff.Jet_lev6_noEMobj));
      Hist.fEvntKtAll_Mjj_b->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),jetsum.Pt());      
    }
  Hist.fEvntNJet_Nvx_b->Fill(nvx12,stuff.myNjet);        
  Hist.fEvntNJet5_Nvx_b->Fill(nvx12,stuff.myNjet_th5);        
  Hist.fEvntNJet10_Nvx_b->Fill(nvx12,stuff.myNjet_th10);        
  Hist.fEvntNJet15_Nvx_b->Fill(nvx12,stuff.myNjet_th15);        
  Hist.fEvntNJet20_Nvx_b->Fill(nvx12,stuff.myNjet_th20);        
  Hist.fEvntNJet25_Nvx_b->Fill(nvx12,stuff.myNjet_th25);        
  if(stuff.myNjet>0)
    {
      Hist.fEvntEt0_Nvx12_b->Fill(nvx12,ave_et0/(stuff.myNjet)); 
      Hist.fEvntEt_Nvx12_b->Fill(nvx12,ave_et/(stuff.myNjet));
    }
  
  if(nvx12==1) Hist.fEvntNjet_Lum_b->Fill((GetHeaderBlock()->InstLum())*1.0E-30,stuff.myNjet); 
  return;
}
/*}}}*/


//_____________________________________________________________________________
//______ filling general histo for particular Jet Cone
//_____________________________________________________________________________
void TMyJetFilterModule::FillJetHistogramsA(JetGeneral_t& Hist,
						JetStuff stuff, double dz, int nvx12)
/*{{{*/
{
  TLorentzVector jetsum(0.0,0.0,0.0,0.0);
  double ave_et=0.0;
  double ave_et0=0.0;
  Hist.fEvntNjet_a->Fill(stuff.myNjet);
  Hist.fEvntNjet5_a->Fill(stuff.myNjet_th5);
  Hist.fEvntNjet10_a->Fill(stuff.myNjet_th10);
  Hist.fEvntNjet15_a->Fill(stuff.myNjet_th15);
  Hist.fEvntNjet20_a->Fill(stuff.myNjet_th20);
  Hist.fEvntNjet25_a->Fill(stuff.myNjet_th25);
  Hist.fEvntdZ_a->Fill(dz);

//________________ Jet Resolution study part ___________________________ 
  if(stuff.Jet_lev6_noEMobj.size()>0)
    {
      if((stuff.Npho_match[0]+stuff.Nele_match[0])==0) 
	{
	  if(stuff.Jet_lev6_noEMobj[0].Pt()>5.0)  Hist.fJetRes_DetEta_j5_a->Fill(fabs(stuff.EtaDet[0]));
	  if(stuff.Jet_lev6_noEMobj[0].Pt()>10.0) Hist.fJetRes_DetEta_j10_a->Fill(fabs(stuff.EtaDet[0]));
	  if(stuff.Jet_lev6_noEMobj[0].Pt()>15.0) Hist.fJetRes_DetEta_j15_a->Fill(fabs(stuff.EtaDet[0]));

	  float _dphi1=fabs(TVector2::Phi_mpi_pi(stuff.myMETcorr_th15.Phi()-stuff.Jet_lev6_noEMobj[0].Phi()));
	  float _dphi2;
	  float dPhi1=100.0;
	  float dPhi2=100.0;
	  if(_dphi1<TMath::Pi()/2.0) dPhi1=_dphi1;
	  else dPhi1=TMath::Pi()-_dphi1;
	  if(stuff.Jet_lev6_noEMobj.size()>1)
	    {
	      if((stuff.Npho_match[1]+stuff.Nele_match[1])==0 && stuff.Jet_lev6_noEMobj[1].Pt()>5.0) 
		{
		  _dphi2=fabs(TVector2::Phi_mpi_pi(stuff.myMETcorr_th15.Phi()-stuff.Jet_lev6_noEMobj[1].Phi()));
		  if(_dphi1<TMath::Pi()/2.0) dPhi1=_dphi1;
		  else dPhi2=TMath::Pi()-_dphi2;
		}
	    }
	  if(dPhi1<0.4 && dPhi2>dPhi1)
	    {
	      double sigma=stuff.Jet_lev6_noEMobj[0].Pt()*GetMyJetResolution(stuff.Jet_lev6_noEMobj[0].Pt(),stuff.EtaDet[0],0,0);
	      double MetSig=-1.0;
	      if(sigma>0.0) MetSig=stuff.myMETcorr_th15.Mod()*cos(dPhi1)/sigma;
	      if(MetSig>=0.0)
		{
		  if(stuff.Jet_lev6_noEMobj[0].Pt()>5.0)  Hist.fJetRes_MetSig_EtaDet_j5_a->Fill(fabs(stuff.EtaDet[0]),MetSig);
		  if(stuff.Jet_lev6_noEMobj[0].Pt()>10.0) Hist.fJetRes_MetSig_EtaDet_j10_a->Fill(fabs(stuff.EtaDet[0]),MetSig);
		  if(stuff.Jet_lev6_noEMobj[0].Pt()>15.0) Hist.fJetRes_MetSig_EtaDet_j15_a->Fill(fabs(stuff.EtaDet[0]),MetSig);
		  if(MetSig>=3.0)
		    {
		      if(stuff.Jet_lev6_noEMobj[0].Pt()>5.0)  Hist.fJetRes_DetEta_j5_Sig3dPhi04_a->Fill(fabs(stuff.EtaDet[0]));
		      if(stuff.Jet_lev6_noEMobj[0].Pt()>10.0) Hist.fJetRes_DetEta_j10_Sig3dPhi04_a->Fill(fabs(stuff.EtaDet[0]));
		      if(stuff.Jet_lev6_noEMobj[0].Pt()>15.0) Hist.fJetRes_DetEta_j15_Sig3dPhi04_a->Fill(fabs(stuff.EtaDet[0]));
		      if(dPhi1<0.1)
			{
			  if(stuff.Jet_lev6_noEMobj[0].Pt()>5.0)  Hist.fJetRes_DetEta_j5_Sig3dPhi01_a->Fill(fabs(stuff.EtaDet[0]));
			  if(stuff.Jet_lev6_noEMobj[0].Pt()>10.0) Hist.fJetRes_DetEta_j10_Sig3dPhi01_a->Fill(fabs(stuff.EtaDet[0]));
			  if(stuff.Jet_lev6_noEMobj[0].Pt()>15.0) Hist.fJetRes_DetEta_j15_Sig3dPhi01_a->Fill(fabs(stuff.EtaDet[0])); 
			}
		    }
		}
	    }
	}
    }


  for(int i=0; i<stuff.Jet_raw_noEMobj.size(); i++)
    { 
      if((stuff.Npho_match[i]+stuff.Nele_match[i])==0) 
	{
	  jetsum=jetsum+stuff.Jet_lev6_noEMobj[i];
	  ave_et=ave_et+stuff.Jet_lev6_noEMobj[i].Pt();
	  ave_et0=ave_et0+stuff.Jet_raw_noEMobj[i].Pt();
	  Hist.fEvntEt0_Njet_a->Fill(stuff.myNjet,stuff.Jet_raw_noEMobj[i].Pt());
	  Hist.fEvntEt_Njet_a->Fill(stuff.myNjet,stuff.Jet_lev6_noEMobj[i].Pt());
	  if(i<2)
	    {
	      Hist.fEvntEt0_a[i]->Fill(stuff.Jet_raw_noEMobj[i].Pt()); 
	      Hist.fEvntEt_a[i]->Fill(stuff.Jet_lev6_noEMobj[i].Pt()); 
	      Hist.fEvntEtaDet_a[i]->Fill(stuff.EtaDet[i]);
	      Hist.fEvntEta_a[i]->Fill(stuff.Jet_lev6_noEMobj[i].Eta());
	      Hist.fEvntPhi_a[i]->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj[i].Phi()));
	      Hist.fEvntTheta_a[i]->Fill(stuff.Jet_lev6_noEMobj[i].Theta());
	      Hist.fEvntEmFr_a[i]->Fill(stuff.EmFrRaw[i]);
	      Hist.fEvntNTowers_a[i]->Fill(stuff.JetNtwr[i]);
	      Hist.fEvntNTracks_a[i]->Fill(stuff.JetNtrk[i]); 
	    }
	  else
	    {
	      Hist.fEvntEt0X_a->Fill(stuff.Jet_raw_noEMobj[i].Pt());
	      Hist.fEvntEtX_a->Fill(stuff.Jet_lev6_noEMobj[i].Pt());        
	      Hist.fEvntEtaDetX_a->Fill(stuff.EtaDet[i]);    
	      Hist.fEvntEtaX_a->Fill(stuff.Jet_lev6_noEMobj[i].Eta());       
	      Hist.fEvntPhiX_a->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj[i].Phi()));       
	      Hist.fEvntThetaX_a->Fill(stuff.Jet_lev6_noEMobj[i].Theta());     
	      Hist.fEvntEmFrX_a->Fill(stuff.EmFrRaw[i]);      
	      Hist.fEvntNTowersX_a->Fill(stuff.JetNtwr[i]);   
	      Hist.fEvntNTracksX_a->Fill(stuff.JetNtrk[i]);
	      Hist.fEvntDeltaR1x_a->Fill(stuff.Jet_lev6_noEMobj[i].DeltaR(stuff.Jet_lev6_noEMobj[0]));
	      Hist.fEvntDeltaR2x_a->Fill(stuff.Jet_lev6_noEMobj[i].DeltaR(stuff.Jet_lev6_noEMobj[1]));
	    }   
	}
    }
  if(stuff.Jet_raw_noEMobj.size()>1 
     && (stuff.Npho_match[0]+stuff.Nele_match[0])==0 
     && (stuff.Npho_match[1]+stuff.Nele_match[1])==0) 
    {
      Hist.fEvntThetaStar_a->Fill(GetThetaStar(stuff.Jet_lev6_noEMobj));
      Hist.fEvntDeltaPhi_a->Fill(fabs(TVector2::Phi_mpi_pi(stuff.Jet_lev6_noEMobj[0].DeltaPhi(stuff.Jet_lev6_noEMobj[1]))));
      Hist.fEvntDeltaEta_a->Fill(stuff.Jet_lev6_noEMobj[0].Eta()-stuff.Jet_lev6_noEMobj[1].Eta());
      Hist.fEvntDeltaR_a->Fill(stuff.Jet_lev6_noEMobj[1].DeltaR(stuff.Jet_lev6_noEMobj[0])); 
      Hist.fEvntMjj_a->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M()); 
      Hist.fEvntKt2jet_a->Fill(GetMyKtKick(stuff.Jet_lev6_noEMobj));  
      Hist.fEvntKtAll_a->Fill(jetsum.Pt());
      Hist.fEvntKt_Mjj_a->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),GetMyKtKick(stuff.Jet_lev6_noEMobj));
      Hist.fEvntKtAll_Mjj_a->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),jetsum.Pt());      
    }
  Hist.fEvntNJet_Nvx_a->Fill(nvx12,stuff.myNjet);        
  Hist.fEvntNJet5_Nvx_a->Fill(nvx12,stuff.myNjet_th5);        
  Hist.fEvntNJet10_Nvx_a->Fill(nvx12,stuff.myNjet_th10);        
  Hist.fEvntNJet15_Nvx_a->Fill(nvx12,stuff.myNjet_th15);        
  Hist.fEvntNJet20_Nvx_a->Fill(nvx12,stuff.myNjet_th20);        
  Hist.fEvntNJet25_Nvx_a->Fill(nvx12,stuff.myNjet_th25);        

  if(stuff.myNjet>0)
    {
      Hist.fEvntEt0_Nvx12_a->Fill(nvx12,ave_et0/(stuff.myNjet)); 
      Hist.fEvntEt_Nvx12_a->Fill(nvx12,ave_et/(stuff.myNjet));
    }
  if(nvx12==1) Hist.fEvntNjet_Lum_a->Fill((GetHeaderBlock()->InstLum())*1.0E-30,stuff.myNjet); 
  return;
}

/*}}}*/

//_____________________________________________________________________________
//______ filling match histo 
//_____________________________________________________________________________
void TMyJetFilterModule::FillMatchingHistograms(MatchStudyHisto_t& Hist,
						JetStuff jetstuff, CommonStuff miscstuff, MatchStuff match)
/*{{{*/
{
  Hist.fMatchNtwr->Fill(jetstuff.JetNtwr[match.JetInd_match]);
  if(match.EmObjType_match==0) 
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

/*}}}*/

//_____________________________________________________________________________
//_____ returns normalized kT_perp. for resolution studies 
//_____________________________________________________________________________
float TMyJetFilterModule::GetMyKtPerp(std::vector<TLorentzVector> vec)
/*{{{*/
{
  float kT_perp=0.0;
  if(vec.size()>=2)
    {
      float phi12=vec[0].DeltaPhi(vec[1]);
      kT_perp=2*cos(phi12/2.0);
    }  
  return kT_perp;
}  
/*}}}*/

//_____________________________________________________________________________
//_____ returns normalized kT_parl. for resolution studies 
//_____________________________________________________________________________
float TMyJetFilterModule::GetMyKtParl(std::vector<TLorentzVector> vec)
/*{{{*/
{
  float kT_parl=0.0;
    if(vec.size()>=2)
    {
      float pt_1=vec[0].Pt();
      float pt_2=vec[1].Pt();
      float phi12=vec[0].DeltaPhi(vec[1]);
      kT_parl=2.0*(pt_1-pt_2)*sin(phi12/2.0)/(pt_1+pt_2);
    }  
  return kT_parl;
}  
/*}}}*/

//_____________________________________________________________________________
//___ reads raw Met, pho, ele info and fills CommonStuff
//_____________________________________________________________________________
void TMyJetFilterModule::DoCommonStuff(CommonStuff &miscstuff)
{

//------------------------------------- get vertex info
  float sumpt=0.0;
	for (int j=0; j<fZVertexBlock->NVertices(); j++) 
	{
  		TStnVertex *vertex= fZVertexBlock->Vertex(j);
      if(vertex->VClass()>=12) myNvx_class12++;
      if((vertex->SumPt())>=sumpt)
		{
			sumpt=vertex->SumPt();
			zvx_best=vertex->Z();
		}
	}


//_______________________________________get my raw MET & SumET info

	// this mat be where the LORENTZVECTOR WARNING IS COMING FROM!
  miscstuff.mySumEt_raw=fMetBlock->Sumet(2);  //sam? raw met is 0, so sumet_raw should be from 0?
  miscstuff.myMET_raw.SetMagPhi(fMetBlock->Met(4),TVector2::Phi_0_2pi(fMetBlock->MetPhi(4)));
  miscstuff.myMET0_raw.SetMagPhi(fMetBlock->Met(0),TVector2::Phi_0_2pi(fMetBlock->MetPhi(0)));


//_______________________________________talk to TInit and get my photon info
  TInit* MyPhoton=
    (TInit*) ((TStnAna*) GetAna()->GetModule("GammaJetsInit")); // "GammaJetsInit" is a default name
	if (MyPhoton == NULL) std::cout << "MyPhoton is NULL" <<std::endl;
	else std::cout << "MyPhoton is GOOD" <<std::endl;
  	int _myNpho=MyPhoton->GetmyNpho();
	std::cout << "_myNpho=" <<  _myNpho << std::endl;

  for(int i=0; i<_myNpho; i++)
    {
      TLorentzVector* _pho_raw = MyPhoton->GetmyUncorrPho(i);
      TLorentzVector* _pho_cor = MyPhoton->GetmyCorrPho(i);
      miscstuff.myPhoInd.push_back(MyPhoton->GetmyPhoInd(i));
      miscstuff.myRawPhoton.push_back(*_pho_raw);
      miscstuff.myCorrPhoton.push_back(*_pho_cor);
      double hadem=MyPhoton->GetphoHadEm3(i);

      miscstuff.myPhoEmFr.push_back(1.0/(hadem+1.0));
      miscstuff.myPhoEtaDet.push_back(MyPhoton->GetphoEtaDet(i));
      miscstuff.myPhoXces.push_back(MyPhoton->GetphoCesX(i));
      miscstuff.myPhoZces.push_back(MyPhoton->GetphoCesZ(i));
    }

  // This part is commented out because accessor to them are not provided yet
// //_______________________________________talk to TMyZeeFilterModule and get my electron info
//   TMyZeeFilterModule* MyZee=
//     (TMyZeeFilterModule*) ((TStnAna*) GetAna()->GetModule("MyZee")); // "MyZee" is a default name
//   int _myNele=MyZee->GetmyNele();
//   for(int i=0; i<_myNele; i++)
//     {
//       TLorentzVector* _ele_raw=MyZee->GetmyEleRaw(i);
//       TLorentzVector* _ele_cor=MyZee->GetmyEle(i);
//       miscstuff.myEleInd.push_back(MyZee->GetmyZeeEleInd(i));
//       miscstuff.myRawElectron.push_back(*_ele_raw); 
//       miscstuff.myCorrElectron.push_back(*_ele_cor);
//       double hadem=MyZee->GetmyEleHadEm(i);
//       miscstuff.myEleEmFr.push_back(1.0/(hadem+1.0));   
//       miscstuff.myEleEtaDet.push_back(MyZee->GetmyEleDetEta(i));
//       miscstuff.myElePhiDet.push_back(MyZee->GetmyEleClusPhi(i));
//       double factor=1.0; // added on 11/07/05
//       double extra_PemEleEt=0.0; // extra PEM ele Et due to difference in energy scale, added on 11/07/04 
//       if(fabs(MyZee->GetmyEleDetEta(i))>1.1) // added on 11/07/05
// 	{
// 	  if(fabs(MyZee->GetmyEleDetEta(i))<1.78) factor *= (MyZee->GetmyEleDetEta(i)>0.0? 1.020 : 1.015);
// 	  else factor *= (MyZee->GetmyEleDetEta(i)>0.0? 1.010 : 1.007);
// 	  extra_PemEleEt=(factor-1.0)*(_ele_raw->Pt());
// 	}
//       else factor=0.0;
//       double pprEt=extra_PemEleEt+factor*(MyZee->GetmyEleCprPpr(i))*
// 	fabs((MyZee->GetmyEleFid(i))*sin(_ele_raw->Theta())); // added on 11/04/05
//       if(pprEt<0.0) pprEt=0.0; // added on 11/04/05
//       miscstuff.myElePprEt.push_back(pprEt); // added on 11/04/05
//     }


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

//_____________________________________________________________________________
//___ corrects Met & SumEt for jets
//_____________________________________________________________________________
void TMyJetFilterModule::DoMyMet(CommonStuff miscstuff, JetStuff &jetstuff)
/*{{{*/
{
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
  for(int i=0; i<miscstuff.myRawPhoton.size(); i++)
    {
      jetstuff.mySumEtCorr_th5=jetstuff.mySumEtCorr_th5-miscstuff.myRawPhoton[i].Pt();  
      jetstuff.mySumEtCorr_th10=jetstuff.mySumEtCorr_th10-miscstuff.myRawPhoton[i].Pt(); 
      jetstuff.mySumEtCorr_th15=jetstuff.mySumEtCorr_th15-miscstuff.myRawPhoton[i].Pt(); 
      jetstuff.mySumEtCorr_th20=jetstuff.mySumEtCorr_th20-miscstuff.myRawPhoton[i].Pt();       
    }
  //_________ correcting all SumEt for electrons 
  for(int i=0; i<miscstuff.myRawElectron.size(); i++)
    {
      jetstuff.mySumEtCorr_th5=jetstuff.mySumEtCorr_th5-miscstuff.myRawElectron[i].Pt();  
      jetstuff.mySumEtCorr_th10=jetstuff.mySumEtCorr_th10-miscstuff.myRawElectron[i].Pt(); 
      jetstuff.mySumEtCorr_th15=jetstuff.mySumEtCorr_th15-miscstuff.myRawElectron[i].Pt(); 
      jetstuff.mySumEtCorr_th20=jetstuff.mySumEtCorr_th20-miscstuff.myRawElectron[i].Pt();  
    }
  //_________ correcting Met for photons and electrons, new as of 09/20/06
  TLorentzVector MetPho(0.0,0.0,0.0,0.0);
  TLorentzVector MetEle(0.0,0.0,0.0,0.0);
  for(int i=0; i<miscstuff.myRawPhoton.size(); i++)
    {
      MetPho=MetPho+miscstuff.myCorrPhoton[i]-miscstuff.myRawPhoton[i];
    }
  for(int i=0; i<miscstuff.myRawElectron.size(); i++)
    {
      MetEle=MetEle+miscstuff.myCorrElectron[i]-miscstuff.myRawElectron[i];
    }
  if(MetPho.E()<MetPho.P() || MetPho.E()<0.0) MetPho.SetPxPyPzE(0.0,0.0,0.0,0.0);
  if(MetEle.E()<MetEle.P() || MetEle.E()<0.0) MetEle.SetPxPyPzE(0.0,0.0,0.0,0.0);

  //_________ correcting all SumEt and calculating Met correction for jets
  TLorentzVector MetJet5(0.0,0.0,0.0,0.0);
  TLorentzVector MetJet10(0.0,0.0,0.0,0.0);
  TLorentzVector MetJet15(0.0,0.0,0.0,0.0);
  TLorentzVector MetJet20(0.0,0.0,0.0,0.0);
  for(int i=0; i<jetstuff.Jet_raw_noEMobj.size(); i++)
    {
      //--------- may need an eta cut here
      if(fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())>=MinJetEta && fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<=MaxJetEta)
	{
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>5.0) 
	    {
	      // corr=U.E.+M.I-raw=lev5-lev6-raw+lev1-lev4, 
	      jetstuff.mySumEtCorr_th5=jetstuff.mySumEtCorr_th5-jetstuff.Jet_raw_noEMobj[i].Pt()
		+jetstuff.Jet_lev5_noEMobj[i].Pt()
		+jetstuff.Jet_lev1_noEMobj[i].Pt()
		-jetstuff.Jet_lev4_noEMobj[i].Pt()
		-jetstuff.Jet_lev6_noEMobj[i].Pt();
	      jetstuff.mySumEtJet_th5=jetstuff.mySumEtJet_th5+jetstuff.Jet_lev5_noEMobj[i].Pt()
		+jetstuff.Jet_lev1_noEMobj[i].Pt()
		-jetstuff.Jet_lev4_noEMobj[i].Pt();
	      MetJet5=MetJet5+jetstuff.Jet_lev5_noEMobj[i]
		+jetstuff.Jet_lev1_noEMobj[i]
		-jetstuff.Jet_lev4_noEMobj[i]
		-jetstuff.Jet_raw_noEMobj[i];
	    } 
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>10.0) 
	    {
	      // corr=lev5-raw+M.I.=lev5-raw+lev1-lev4, 
	      jetstuff.mySumEtCorr_th10=jetstuff.mySumEtCorr_th10-jetstuff.Jet_raw_noEMobj[i].Pt()
		+jetstuff.Jet_lev5_noEMobj[i].Pt()
		+jetstuff.Jet_lev1_noEMobj[i].Pt()
		-jetstuff.Jet_lev4_noEMobj[i].Pt()
		-jetstuff.Jet_lev6_noEMobj[i].Pt();
	      jetstuff.mySumEtJet_th10=jetstuff.mySumEtJet_th10+jetstuff.Jet_lev5_noEMobj[i].Pt()
		+jetstuff.Jet_lev1_noEMobj[i].Pt()
		-jetstuff.Jet_lev4_noEMobj[i].Pt();
	      MetJet10=MetJet10+jetstuff.Jet_lev5_noEMobj[i]
		+jetstuff.Jet_lev1_noEMobj[i]
		-jetstuff.Jet_lev4_noEMobj[i]
		-jetstuff.Jet_raw_noEMobj[i];
	    } 
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>15.0) 
	    {
	      // corr=lev5-raw+M.I.=lev5-raw+lev1-lev4, 
	      jetstuff.mySumEtCorr_th15=jetstuff.mySumEtCorr_th15-jetstuff.Jet_raw_noEMobj[i].Pt()
		+jetstuff.Jet_lev5_noEMobj[i].Pt()
		+jetstuff.Jet_lev1_noEMobj[i].Pt()
		-jetstuff.Jet_lev4_noEMobj[i].Pt()
		-jetstuff.Jet_lev6_noEMobj[i].Pt();
	      jetstuff.mySumEtJet_th15=jetstuff.mySumEtJet_th15+jetstuff.Jet_lev5_noEMobj[i].Pt()
		+jetstuff.Jet_lev1_noEMobj[i].Pt()
		-jetstuff.Jet_lev4_noEMobj[i].Pt();
	      MetJet15=MetJet15+jetstuff.Jet_lev5_noEMobj[i]
		+jetstuff.Jet_lev1_noEMobj[i]
		-jetstuff.Jet_lev4_noEMobj[i]
		-jetstuff.Jet_raw_noEMobj[i];
	    } 
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>20.0) 
	    {
	      // corr=lev5-raw+M.I.=lev5-raw+lev1-lev4, 
	      jetstuff.mySumEtCorr_th20=jetstuff.mySumEtCorr_th20-jetstuff.Jet_raw_noEMobj[i].Pt()
		+jetstuff.Jet_lev5_noEMobj[i].Pt()
		+jetstuff.Jet_lev1_noEMobj[i].Pt()
		-jetstuff.Jet_lev4_noEMobj[i].Pt()
		-jetstuff.Jet_lev6_noEMobj[i].Pt();
	      jetstuff.mySumEtJet_th20=jetstuff.mySumEtJet_th20+jetstuff.Jet_lev5_noEMobj[i].Pt()
		+jetstuff.Jet_lev1_noEMobj[i].Pt()
		-jetstuff.Jet_lev4_noEMobj[i].Pt();
	      MetJet20=MetJet20+jetstuff.Jet_lev5_noEMobj[i]
		+jetstuff.Jet_lev1_noEMobj[i]
		-jetstuff.Jet_lev4_noEMobj[i]
		-jetstuff.Jet_raw_noEMobj[i];
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
/*}}}*/

//_____________________________________________________________________________
//___ does my jets; main routine  
//_____________________________________________________________________________
void TMyJetFilterModule::DoMyJet(TStnJetBlock* fJetBlock, CommonStuff miscstuff,
							JetStuff &jetstuff, MatchStudyHisto_t& Hist)
{
  DoMyJetNoMatch(fJetBlock,jetstuff);
  DoMyJetWithMatch(fJetBlock,miscstuff,jetstuff,Hist);
  ReorderMyJets(jetstuff);
  return;
}

//_____________________________________________________________________________
//--the following three functions could have been replaced by one which works
//--with generic data type (to be done)
//____ SWAPPING FUNCTIONS _________________________________________________
//_____________________________________________________________________________
void TMyJetFilterModule::myExchange_tlv(TLorentzVector& val1, TLorentzVector& val2) { //exchanges val1 and val2
  TLorentzVector dummy=val1;
  val1=val2;
  val2=dummy;
  return;
}
//_____________________________________________________________________________
//_____________________________________________________________________________
void TMyJetFilterModule::myExchange_dbl(double& val1, double& val2) { //exchanges val1 and val2
  double dummy=val1;
  val1=val2;
  val2=dummy;
  return;
}
//_____________________________________________________________________________
//_____________________________________________________________________________
void TMyJetFilterModule::myExchange_int(int& val1, int& val2) { //exchanges val1 and val2
  int dummy=val1;
  val1=val2;
  val2=dummy;
  return;
}


//_____________________________________________________________________________
//______ reorders jets after removing EM objects
//_____________________________________________________________________________
void TMyJetFilterModule::ReorderMyJets(JetStuff &jetstuff)
/*{{{*/
{
  int reorder_word=0;
  int reorder_word_last=0;
  for(int j=0; j<jetstuff.myNjet-1; j++)
    {
      reorder_word_last=reorder_word;
      for(int i=1; i<jetstuff.myNjet; i++)
	{
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>jetstuff.Jet_lev6_noEMobj[i-1].Pt()) 
	    {
	      myExchange_tlv(jetstuff.Jet_raw[i-1],jetstuff.Jet_raw[i]);
	      myExchange_tlv(jetstuff.Jet_lev1[i-1],jetstuff.Jet_lev1[i]);
	      myExchange_tlv(jetstuff.Jet_lev4[i-1],jetstuff.Jet_lev4[i]);
	      myExchange_tlv(jetstuff.Jet_lev5[i-1],jetstuff.Jet_lev5[i]);
	      myExchange_tlv(jetstuff.Jet_lev6[i-1],jetstuff.Jet_lev6[i]);
	      myExchange_tlv(jetstuff.Jet_lev7[i-1],jetstuff.Jet_lev7[i]);
	      myExchange_tlv(jetstuff.Jet_raw_noEMobj[i-1],jetstuff.Jet_raw_noEMobj[i]);
	      myExchange_tlv(jetstuff.Jet_lev1_noEMobj[i-1],jetstuff.Jet_lev1_noEMobj[i]);
	      myExchange_tlv(jetstuff.Jet_lev4_noEMobj[i-1],jetstuff.Jet_lev4_noEMobj[i]);
	      myExchange_tlv(jetstuff.Jet_lev5_noEMobj[i-1],jetstuff.Jet_lev5_noEMobj[i]);
	      myExchange_tlv(jetstuff.Jet_lev6_noEMobj[i-1],jetstuff.Jet_lev6_noEMobj[i]);
	      myExchange_tlv(jetstuff.Jet_lev7_noEMobj[i-1],jetstuff.Jet_lev7_noEMobj[i]);
	      
	      myExchange_dbl(jetstuff.EtaDet[i-1],jetstuff.EtaDet[i]);
	      myExchange_dbl(jetstuff.EtaDetCorr[i-1],jetstuff.EtaDetCorr[i]);
	      myExchange_dbl(jetstuff.EmFrRaw[i-1],jetstuff.EmFrRaw[i]);
	      myExchange_dbl(jetstuff.EmFrCorr[i-1],jetstuff.EmFrCorr[i]);

	      myExchange_int(jetstuff.JetNtrk[i-1],jetstuff.JetNtrk[i]);
	      myExchange_int(jetstuff.JetNtwr[i-1],jetstuff.JetNtwr[i]);
	      myExchange_int(jetstuff.Nobj_match[i-1],jetstuff.Nobj_match[i]); 
	      myExchange_int(jetstuff.Npho_match[i-1],jetstuff.Npho_match[i]); 
	      myExchange_int(jetstuff.Nele_match[i-1],jetstuff.Nele_match[i]); 
	      myExchange_int(jetstuff.Nmu_match[i-1],jetstuff.Nmu_match[i]);
	      myExchange_int(jetstuff.Ntau_match[i-1],jetstuff.Ntau_match[i]); 
	      myExchange_int(jetstuff.Nbtag_match[i-1],jetstuff.Nbtag_match[i]);
	      myExchange_int(jetstuff.JetBlockInd[i-1],jetstuff.JetBlockInd[i]);

	      reorder_word++;
	    }
	}
      if(reorder_word_last==reorder_word) break;
    }  
  return;
}
/*}}}*/


//_____________________________________________________________________________
//______ re-calculates jet EmFr after removing EM object
//_____________________________________________________________________________
double TMyJetFilterModule::CorrectEmFr(double emfr_old, double emfr_obj,
						double e_old, double e_obj)
/*{{{*/
{
  double emfr_new=emfr_old;
  if(e_old>e_obj)
    {
      double term1=emfr_obj*e_obj/e_old;
      double denom=1.0-e_obj/e_old;
      emfr_new=(emfr_old-term1)/denom;
    }
  return emfr_new;
}
/*}}}*/

//_____________________________________________________________________________
//____ re-calculates jet EtaDet after removing EM object
//  this function uses Snowmass convention to calculate jet centroid (true for JetClu only)
//_____________________________________________________________________________
double TMyJetFilterModule::CorrectEtaDet(TLorentzVector *vec_pho,
						TLorentzVector *vec_old, double jetetadet_old,
						double pho_etadet)
/*{{{*/
{
  double E_old=vec_old->E();
  double E_pho=vec_pho->E();
  double theta_det_old=2.0*atan(exp(-jetetadet_old));
  double theta_det_pho=2.0*atan(exp(-pho_etadet));
  double etadet_new=jetetadet_old;
  double det=E_old*fabs(sin(theta_det_old))-E_pho*fabs(sin(theta_det_pho));
  double num=jetetadet_old*E_old*fabs(sin(theta_det_old))-pho_etadet*E_pho*fabs(sin(theta_det_pho));
  if(det>0.1) etadet_new=num/det;
  return etadet_new;
}
/*}}}*/


//_____________________________________________________________________________
//_____ does my jets after removing EM objetcs, to be called after DoMyJetNoMatch
//_____________________________________________________________________________
void TMyJetFilterModule::DoMyJetWithMatch(TStnJetBlock* fJetBlock,
				CommonStuff miscstuff, JetStuff &jetstuff, MatchStudyHisto_t& Hist)
/*{{{*/
{

  int _jetcone=0;

  if(fabs(fJetBlock->ConeSize()-0.4)<0.01) _jetcone=0;
  if(fabs(fJetBlock->ConeSize()-0.7)<0.01) _jetcone=1;
  if(fabs(fJetBlock->ConeSize()-1.0)<0.01) _jetcone=2;  
  fJTC_coneSize=_jetcone;
  //_________________________ Initializing Jet corrections
  CorInit lev1;
  CorInit lev4;		//sam? why no lev 2,3???
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

  MatchStuff dummymatch;
  ClearMatchStuff();

  int _Npho_match=0;
  int _Nele_match=0;
  for(int i=0; i<jetstuff.myNjet; i++)
    {
      float _etadet_tmp=jetstuff.EtaDet[i]; // new line, 10/20/05 
      TLorentzVector *_jet=&jetstuff.Jet_raw[i];
      //_____ removing photons
      for(int j=0; j<miscstuff.myRawPhoton.size() && _Npho_match<miscstuff.myRawPhoton.size(); j++)
	{
	  TLorentzVector *_pho=&miscstuff.myRawPhoton[j];
	  double _match_dR;
	  double _match_dEta;
	  double _match_dPhi;
	  int match_stat=MyMatchPhoJet(jetstuff.JetBlockInd[i],miscstuff.myPhoInd[j],-1,
				       fJetBlock,_pho,_jet,fJetBlock->ConeSize(),
				       _match_dR,_match_dPhi,_match_dEta);
	  if(match_stat==1)
	    {
	      _Npho_match++;
	      jetstuff.Nobj_match[i]=jetstuff.Nobj_match[i]+1;
	      jetstuff.Npho_match[i]=jetstuff.Npho_match[i]+1;
	      if((jetstuff.Npho_match[i]+jetstuff.Nele_match[i])==1) // filling only first match
		{
		  dummymatch.JetInd_match=i;
		  dummymatch.EmInd_match=j;
		  dummymatch.EmObjType_match=0; // photons
		  matchstuff.push_back(dummymatch);
		}
	      // here, I re-calculate jet EmFr 
	      jetstuff.EmFrCorr[i]=CorrectEmFr(jetstuff.EmFrCorr[i],miscstuff.myPhoEmFr[j],
					       jetstuff.Jet_raw_noEMobj[i].E(),miscstuff.myRawPhoton[j].E());
	      // corrected EtaDet is to be calculated before removing EM object 
	      _etadet_tmp=CorrectEtaDet(_pho,&jetstuff.Jet_raw_noEMobj[i], 
					jetstuff.EtaDetCorr[i],miscstuff.myPhoEtaDet[j]); // new line, 10/20/05
	      jetstuff.EtaDetCorr[i]=_etadet_tmp;  // new line, 10/20/05
	      if(jetstuff.Jet_raw_noEMobj[i].E()>miscstuff.myRawPhoton[j].E()) 
		jetstuff.Jet_raw_noEMobj[i]=jetstuff.Jet_raw_noEMobj[i]-miscstuff.myRawPhoton[j]; // removing raw photon from raw jet
	      else jetstuff.Jet_raw_noEMobj[i].SetPxPyPzE(0.0,0.0,0.0,0.0);
	      if(jetstuff.Jet_raw_noEMobj[i].E()<jetstuff.Jet_raw_noEMobj[i].P()
		 || jetstuff.Jet_raw_noEMobj[i].E()<=0.0
		 || (jetstuff.Jet_raw[i].DeltaR(jetstuff.Jet_raw_noEMobj[i]))>2.0*(fJetBlock->ConeSize())) 
		{
		  jetstuff.Jet_raw_noEMobj[i].SetPxPyPzE(0.0,0.0,0.0,0.0); // in case removing EM object leads to M()<0.0
		}
	    }
	}
      //_____ removing electrons
      for(int j=0; j<miscstuff.myRawElectron.size() && _Nele_match<miscstuff.myRawElectron.size(); j++)
	{
	  TLorentzVector *_ele=&miscstuff.myRawElectron[j];
	  double _match_dR;
	  double _match_dEta;
	  double _match_dPhi;
	  int match_stat=MyMatchPhoJet(jetstuff.JetBlockInd[i],-1,miscstuff.myEleInd[j],
				       fJetBlock,_ele,_jet,fJetBlock->ConeSize(),
				       _match_dR,_match_dPhi,_match_dEta);
	  if(match_stat==1)
	    {
	      _Nele_match++;
	      jetstuff.Nobj_match[i]=jetstuff.Nobj_match[i]+1;
	      jetstuff.Nele_match[i]=jetstuff.Nele_match[i]+1;
	      if((jetstuff.Nele_match[i]+jetstuff.Npho_match[i])==1) // filling only first match
		{
		  dummymatch.JetInd_match=i;
		  dummymatch.EmInd_match=j;
		  dummymatch.EmObjType_match=1; // electrons
		  matchstuff.push_back(dummymatch);
		}
	      // here, I re-calculate jet EmFr 
	      jetstuff.EmFrCorr[i]=CorrectEmFr(jetstuff.EmFrCorr[i],miscstuff.myEleEmFr[j],
					       jetstuff.Jet_raw_noEMobj[i].E(),miscstuff.myRawElectron[j].E());
	      // corrected EtaDet is to be calculated before removing EM object 
	      _etadet_tmp=CorrectEtaDet(_ele,&jetstuff.Jet_raw_noEMobj[i], 
					jetstuff.EtaDetCorr[i],miscstuff.myEleEtaDet[j]); // new line, 10/20/05
	      jetstuff.EtaDetCorr[i]=_etadet_tmp;  // new line, 10/20/05
	      if(jetstuff.Jet_raw_noEMobj[i].E()>miscstuff.myRawElectron[j].E()) 
		jetstuff.Jet_raw_noEMobj[i]=jetstuff.Jet_raw_noEMobj[i]-miscstuff.myRawElectron[j]; // removing raw electron from raw jet
	      else jetstuff.Jet_raw_noEMobj[i].SetPxPyPzE(0.0,0.0,0.0,0.0);
	      if(jetstuff.Jet_raw_noEMobj[i].E()<jetstuff.Jet_raw_noEMobj[i].P() 
		 || jetstuff.Jet_raw_noEMobj[i].E()<=0.0
		 || (jetstuff.Jet_raw[i].DeltaR(jetstuff.Jet_raw_noEMobj[i]))>2.0*(fJetBlock->ConeSize())) 
		{
		  jetstuff.Jet_raw_noEMobj[i].SetPxPyPzE(0.0,0.0,0.0,0.0); // in case removing EM object leads to M()<0.0
		}
	    }
	}
      if((jetstuff.Nele_match[i]+jetstuff.Npho_match[i])>0) // correcting jet energy after removing EM objetcs
	{
	  TLorentzVector _rawJet=jetstuff.Jet_raw_noEMobj[i];
	  TLorentzVector _rawJet_tmp;
	  _rawJet_tmp=_rawJet;
	  float _myEmFr_tmp=jetstuff.EmFrCorr[i];
	  double corr7=GetCorrection(lev7,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
	  _myEmFr_tmp=jetstuff.EmFrCorr[i]; // do this again because EmFr can be potentially changed in GetCorrection
	  _rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
	  double corr6=GetCorrection(lev6,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
	  _myEmFr_tmp=jetstuff.EmFrCorr[i]; // do this again because EmFr can be potentially changed in GetCorrection
	  _rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
	  double corr5=GetCorrection(lev5,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
	  _myEmFr_tmp=jetstuff.EmFrCorr[i]; // do this again because EmFr can be potentially changed in GetCorrection
	  _rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
	  double corr4=GetCorrection(lev4,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
	  _myEmFr_tmp=jetstuff.EmFrCorr[i]; // do this again because EmFr can be potentially changed in GetCorrection
	  _rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
	  double corr1=GetCorrection(lev1,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);	  
	  jetstuff.Jet_lev1_noEMobj[i]=_rawJet*corr1; 
	  jetstuff.Jet_lev4_noEMobj[i]=_rawJet*corr4; 
	  jetstuff.Jet_lev5_noEMobj[i]=_rawJet*corr5; 
	  jetstuff.Jet_lev6_noEMobj[i]=_rawJet*corr6; 
	  jetstuff.Jet_lev7_noEMobj[i]=_rawJet*corr7; 	  
	}
      if(fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())>=MinJetEta && fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<=MaxJetEta)
	{
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>5.0) jetstuff.myNjet_th5++;
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>10.0) jetstuff.myNjet_th10++;
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>15.0) jetstuff.myNjet_th15++;
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>20.0) jetstuff.myNjet_th20++;
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()>25.0) jetstuff.myNjet_th25++;
	}
    }
  if(_Nele_match!=miscstuff.myRawElectron.size() || _Npho_match!=miscstuff.myRawPhoton.size()) bad_EMjet_match_flag=0;
  for(int i=0; i<matchstuff.size(); i++)
    {
      FillMatchingHistograms(Hist,jetstuff,miscstuff,matchstuff[i]);
    }
  return;
}

/*}}}*/


//_____________________________________________________________________________
//______ does my jets, but doesn't remove EM objetcs
//_____________________________________________________________________________
void TMyJetFilterModule::DoMyJetNoMatch(TStnJetBlock* fJetBlock,
				JetStuff &jetstuff)
/*{{{*/
{

  int _jetcone=0;
  if(fabs(fJetBlock->ConeSize()-0.4)<0.01) _jetcone=0;  // why reluctant to use equality check??
  if(fabs(fJetBlock->ConeSize()-0.7)<0.01) _jetcone=1;
  if(fabs(fJetBlock->ConeSize()-1.0)<0.01) _jetcone=2;

  fJTC_coneSize=_jetcone;
  //_________________________ Initializing Jet corrections
  CorInit lev1;	//CorInit-- param for jet correction inits
  CorInit lev4;
  CorInit lev5;
  CorInit lev6;
  CorInit lev7;
  lev1.level	= 1;
  lev1.nvx		= myNvx_class12;
  lev1.cone		= fJTC_coneSize;
  lev1.version	= fJTC_version;	//sam: version of jet corrections: 0--MC, 1..5(?)--Data
  lev1.sys		= fJTC_systcode; 	//sam: code foe systematic corrections
  lev1.imode	= fJTC_imode;  	//sam: correction mode in V5* data and MC: 0=MC, 1=data
  lev1.Nrun		= myJetRun;			//sam: to be used in corrections and for crosschecks
  lev4.level	= 4;
  lev4.nvx		= myNvx_class12;
  lev4.cone		= fJTC_coneSize;
  lev4.version	= fJTC_version;
  lev4.sys		= fJTC_systcode;  
  lev4.imode	= fJTC_imode;  
  lev4.Nrun		= myJetRun;
  lev5.level	= 5;
  lev5.nvx		= myNvx_class12;
  lev5.cone		= fJTC_coneSize;
  lev5.version	= fJTC_version;
  lev5.sys		= fJTC_systcode;
  lev5.imode	= fJTC_imode;    
  lev5.Nrun		= myJetRun; 
  lev6.level	= 6;
  lev6.nvx		= myNvx_class12;
  lev6.cone		= fJTC_coneSize;
  lev6.version	= fJTC_version;
  lev6.sys		= fJTC_systcode;
  lev6.imode	= fJTC_imode;  
  lev6.Nrun		= myJetRun;
  lev7.level	= 7;
  lev7.nvx		= myNvx_class12;
  lev7.cone		= fJTC_coneSize;
  lev7.version	= fJTC_version;
  lev7.sys		= fJTC_systcode;
  lev7.imode	= fJTC_imode;  
  lev7.Nrun		= myJetRun;

  jetstuff.myNjet = fJetBlock->NJets();
  for(int i=0; i<jetstuff.myNjet; i++)
    {
      TStnJet* jet = fJetBlock->Jet(i);
      jetstuff.JetBlockInd.push_back(i);
      jetstuff.JetNtrk.push_back(jet->NTracks());
      jetstuff.JetNtwr.push_back(jet->NTowers());
      jetstuff.EtaDet.push_back(jet->DetEta());
      jetstuff.EtaDetCorr.push_back(jet->DetEta()); // for a moment, make it the same as "EtaDet"
      jetstuff.EmFrRaw.push_back(jet->Emfr());   //sam: check if these 2 r used differently
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
/*}}}*/

//____________________ Jet Resolution: default function _________________
//__________ Function: sqrt(p0/Et+p1/(Et*Et)+p2), using 0d,h,i datasets of dijet and Z+jet events
//__________ def cuts: dphi>2.7, Njet(extra)<=2
//____________________ input params: Pt(lev6), jet cone size (0=0.4, 1=0.7, 2=1.0)
double TMyJetFilterModule::MyJetRes_def(double jetPtLev6, int jetcone)
/*{{{*/
{
  double def=0.0;
  double par[3];
  double pg[3];
  if(jetcone==0) // JetClu-0.4 at this moment
    {
      par[0]=0.020;
      par[1]=21.9;
      par[2]=0.0092;      
      //_______________ params for Gauss at low Pt
      pg[0]=0.655;
      pg[1]=1.49;
      pg[2]=11.02;
    }
  else
    {
      par[0]=0.0;
      par[1]=0.0;
      par[2]=0.0;      
      //_______________ params for Gauss at low Pt
      pg[0]=0.0;
      pg[1]=0.0;
      pg[2]=1.0E6;
    }
  if(jetPtLev6>0.0)
    {
      double arg=0.0;
      arg=(par[0]/jetPtLev6+par[1]/(jetPtLev6*jetPtLev6)+par[2]);
      double val1=sqrt(arg);
      arg=0.5*(jetPtLev6-pg[1])*(jetPtLev6-pg[1])/(pg[2]*pg[2]);
      double val2=pg[0]*exp(-arg);
      if(val2<val1) def=val2; // lowPt fix
      else def=val1;
      if(jetPtLev6>12.5) def=val1; 
    }
  return def;
}
/*}}}*/

//____________________ Jet Resolution: step function _________________
//_________ returns 1 in the region where fit function is used; 
//_________ returns 0 in the region where a fix to low Pt is applied 
double TMyJetFilterModule::MyJetRes_step(double jetPtLev6, int jetcone)
/*{{{*/
{
  double def=1.0;
  double par[3];
  double pg[3];
  if(jetcone==0) // JetClu-0.4 at this moment
    {
      par[0]=0.020;
      par[1]=21.9;
      par[2]=0.0092;      
      //_______________ params for Gauss at low Pt
      pg[0]=0.655;
      pg[1]=1.49;
      pg[2]=11.02;
    }
  else
    {
      par[0]=0.0;
      par[1]=0.0;
      par[2]=0.0;      
      //_______________ params for Gauss at low Pt
      pg[0]=0.0;
      pg[1]=0.0;
      pg[2]=1.0E6;
    }
  if(jetPtLev6>0.0)
    {
      double arg=0.0;
      arg=(par[0]/jetPtLev6+par[1]/(jetPtLev6*jetPtLev6)+par[2]);
      double val1=sqrt(arg);
      arg=0.5*(jetPtLev6-pg[1])*(jetPtLev6-pg[1])/(pg[2]*pg[2]);
      double val2=pg[0]*exp(-arg);
      if(val2<val1) def=0.0; // lowPt fix
      else def=1.0;
      if(jetPtLev6>12.5) def=1.0; 
    }
  return def;
}
/*}}}*/

//_________________________ Jet Resolution: systematics-2 function ________
//__________Syst2: "syst1(0d,h,i)"-def(0d,h,i)
//__________ syst1 cuts: dphi>2.9, Njet(extra)<=2
double TMyJetFilterModule::MyJetRes_syst2(double jetPtLev6, int jetcone)
/*{{{*/
{
  double syst=0.0;
  double def=MyJetRes_def(jetPtLev6,jetcone);
  double par[3];
  double pg[3];
  if(jetcone==0) // JetClu-0.4 at this moment
    {
      par[0]=0.279;
      par[1]=22.0;
      par[2]=0.00812;
      //_______________ params for Gauss at low Pt
      pg[0]=0.689;
      pg[1]=1.63;
      pg[2]=11.3;      
    }
  else
    {
      par[0]=0.0;
      par[1]=0.0;
      par[2]=0.0;
      //_______________ params for Gauss at low Pt
      pg[0]=0.0;
      pg[1]=0.0;
      pg[2]=1.0E6;      
    }
  if(jetPtLev6>0.0)
    {
      double arg=0.0;
      arg=(par[0]/jetPtLev6+par[1]/(jetPtLev6*jetPtLev6)+par[2]);
      double val1=sqrt(arg);
      arg=0.5*(jetPtLev6-pg[1])*(jetPtLev6-pg[1])/(pg[2]*pg[2]);
      double val2=pg[0]*exp(-arg);
      double dev=def;
      if(val2<val1) dev=val2; // lowPt fix
      else dev=val1;
      if(jetPtLev6>12.5) dev=val1; 
      syst=fabs(dev-def);
    }
  return syst;
}
/*}}}*/

//_________________________ Jet Resolution: systematics-3 function ________
//__________Syst3: "syst2(0d,h,i)"-def(0d,h,i)
//__________ syst3 cuts: dphi>2.7, Njet(extra)<=1
double TMyJetFilterModule::MyJetRes_syst3(double jetPtLev6, int jetcone)
/*{{{*/
{
  double syst=0.0;
  double def=MyJetRes_def(jetPtLev6,jetcone);
  double par[3];
  double pg[3];
  if(jetcone==0) // JetClu-0.4 at this moment
    {
      par[0]=-0.124;
      par[1]=20.8;
      par[2]=0.00991;
      //_______________ params for Gauss at low Pt
      pg[0]=0.655;
      pg[1]=2.06;
      pg[2]=9.09;      
    }
  else
    {
      par[0]=0.0;
      par[1]=0.0;
      par[2]=0.0;
      //_______________ params for Gauss at low Pt
      pg[0]=0.0;
      pg[1]=0.0;
      pg[2]=1.0E6;            
    }
  if(jetPtLev6>0.0)
    {
      double arg=0.0;
      arg=(par[0]/jetPtLev6+par[1]/(jetPtLev6*jetPtLev6)+par[2]);
      double val1=sqrt(arg);
      arg=0.5*(jetPtLev6-pg[1])*(jetPtLev6-pg[1])/(pg[2]*pg[2]);
      double val2=pg[0]*exp(-arg);
      double dev=def;
      if(val2<val1) dev=val2; // lowPt fix
      else dev=val1;
      if(jetPtLev6>12.5) dev=val1; 
      syst=fabs(dev-def);
    }
  return syst;
}
/*}}}*/

//_____________________ Jet Resolution: systematics-4 function ____________
// Syst4="p0/Et+p1"-def(0d,h,i)
//__________ def cuts: dphi>2.7, Njet(extra)<=2
double TMyJetFilterModule::MyJetRes_syst4(double jetPtLev6, int jetcone)
/*{{{*/
{
  double syst=0.0;
  double def=MyJetRes_def(jetPtLev6,jetcone);
  double def_switch=MyJetRes_step(jetPtLev6,jetcone);
  double par[2];
  if(jetcone==0) // JetClu-0.4 at this moment
    {
      par[0]=3.942;
      par[1]=0.08323;
    }
  else
    {
      par[0]=0.0;
      par[1]=0.0;
    }
  if(jetPtLev6>0.0)
    {
      double arg=0.0;
      arg=(par[0]/jetPtLev6+par[1]);
      syst=fabs(arg-def)*def_switch;
    }
  return syst;
}
/*}}}*/

//_____________________ Jet Resolution: stat. uncertainty function __________
double TMyJetFilterModule::MyJetRes_stat(double jetPtLev6, int jetcone)
/*{{{*/
{
  double stat=0.0;
  double def=MyJetRes_def(jetPtLev6,jetcone);
  double par_er[3];
  double c[3][3];
  c[0][0]=c[1][1]=c[2][2]=1.0;
  if(jetcone==0) // JetClu-0.4 at this moment
    {
      if(MyJetRes_step(jetPtLev6,jetcone)>0.5)
	{
	  par_er[0]=0.014;
	  par_er[1]=0.3;
	  par_er[2]=0.0001;
	  c[0][0]=c[1][1]=c[2][2]=1.0;
	  c[0][1]=c[1][0]=-0.958;
	  c[0][2]=c[2][0]=-0.924;
	  c[1][2]=c[2][1]=0.825;
	}
      else // low Et fix
	{
	  par_er[0]=0.014;
	  par_er[1]=0.99;
	  par_er[2]=0.94;
	  c[0][0]=c[1][1]=c[2][2]=1.0;
	  c[0][1]=c[1][0]=-0.834;
	  c[0][2]=c[2][0]=0.657;
	  c[1][2]=c[2][1]=-0.945;
	}
    }
  else
    {
      par_er[0]=0.0;
      par_er[1]=0.0;
      par_er[2]=0.0;
      c[0][1]=c[1][0]=0.0;
      c[0][2]=c[2][0]=0.0;
      c[1][2]=c[2][1]=0.0;      
    }
  if(jetPtLev6>0.0 && def>0.0)
    {
      if(MyJetRes_step(jetPtLev6,jetcone)>0.5)
	{
	  double arg=0.0;
	  double x=jetPtLev6;
	  double xx=jetPtLev6*jetPtLev6;
	  double xxx=jetPtLev6*jetPtLev6*jetPtLev6;
	  double xxxx=jetPtLev6*jetPtLev6*jetPtLev6*jetPtLev6;
	  arg=(par_er[0]*par_er[0]/xx+par_er[1]*par_er[1]/xxxx+par_er[2]*par_er[2]
	       +2.0*c[0][1]*par_er[0]*par_er[1]/xxx
	       +2.0*c[0][2]*par_er[0]*par_er[2]/x
	       +2.0*c[1][2]*par_er[1]*par_er[2]/xx);
	  if(arg>0.0) stat=sqrt(arg)/(2.0*def);
	  else stat=0.0;
	}
      else // low Et fix
	{
	  double x=jetPtLev6;
	  double t1=(x-1.49)/(11.02*11.02);
	  double t2=(x-1.49)*(x-1.49)/(11.02*11.02*11.02);
	  double arg=(par_er[0]*par_er[0]+par_er[1]*par_er[1]*t1*t1+par_er[2]*par_er[2]*t2*t2
		      +2.0*c[0][1]*par_er[0]*par_er[1]*t1
		      +2.0*c[0][2]*par_er[0]*par_er[2]*t2
		      +2.0*c[1][2]*par_er[1]*par_er[2]*t1*t2);
	  if(arg>0.0) stat=sqrt(arg)*exp(-0.5*(x-1.49)*(x-1.49)/(11.02*11.02));
	  else stat=0.0;
	}
    }
  return stat;
}
/*}}}*/

//_______________________ Jet Resolution: total uncertainty function _______
double TMyJetFilterModule::MyJetRes_total(double jetPtLev6, int jetcone)
/*{{{*/
{
//   double syst1=MyJetRes_syst1(jetPtLev6,jetcone);
  double syst1=0.0; // there is no "syst1" in this version. "syst1"=def(0i)-def(0d,h)
  double syst2=MyJetRes_syst2(jetPtLev6,jetcone);
  double syst3=MyJetRes_syst3(jetPtLev6,jetcone);
  double stat=MyJetRes_stat(jetPtLev6,jetcone);
  double tot=sqrt(syst1*syst1+syst2*syst2+syst3*syst3+stat*stat);
  return tot;
}
/*}}}*/


//______ Jet Resolution: correction for eta dependence of resolution ______
double TMyJetFilterModule::MyJetRes_EtaCor(double deteta, int jetcone)
{
  double corr=1.0; // no correction for now 
  return corr;
}
//______ Jet Resolution: scale factor Had/Det _______________________________
double TMyJetFilterModule::MyJetRes_Had2Det(double jetPtLev6, int jetcone)
{
  // p0=0.513+-0.016; p1=0.0431+-0.0036
//   double scale_Had2Det=1.0-0.513*exp(-0.0431*jetPtLev6); // based on Zee MC: Resolution(response)/Resolution(Bisector) 
  double scale_Had2Det=1.0; // no scalling 
  return scale_Had2Det;
}

//____________________ this function returns my jet energy resolution
//____________________ input params: Pt(lev6), jet cone size (0=0.4, 1=0.7, 2=1.0), 
//____________________               status code (0=default, -1=-sigma, 1=+sigma)  
double TMyJetFilterModule::GetMyJetResolution(double jetPtLev6, double deteta,
						int jetcone, int statcode)
{
  double jetSigma=MyJetRes_def(jetPtLev6,jetcone);
  double jetSigma_err=statcode*MyJetRes_total(jetPtLev6,jetcone);
  if((jetSigma+jetSigma_err)>0.0) jetSigma=(jetSigma+jetSigma_err)*MyJetRes_EtaCor(deteta,jetcone)*MyJetRes_Had2Det(jetPtLev6,jetcone);
  else jetSigma=jetSigma*MyJetRes_EtaCor(deteta,jetcone)*MyJetRes_Had2Det(jetPtLev6,jetcone);
  if(jetSigma<0.0) jetSigma=0.15; 
  return jetSigma;
}


//___ generates Sigma and Mean for "unclustered" Met
//___ takes into account correlations between params
void TMyJetFilterModule::GetMyUnclMetResolution(double sumEt, int systcode,
						double &sigmaMx, double &sigmaMy, double &meanMx,
						double &meanMy)
/*{{{*/
{
  double sigmaX_err=0.0;
  double meanX_err=0.0;
  double sigmaY_err=0.0;
  double meanY_err=0.0;
  double dummypar=0.0;

  double systX_sigma=metmodel_sigmaX[fUnclParamSwitch][0]+metmodel_sigmaX[fUnclParamSwitch][1]*sqrt(sumEt)-
    metmodel_sigmaX[1-fUnclParamSwitch][0]-metmodel_sigmaX[1-fUnclParamSwitch][1]*sqrt(sumEt);
  double systY_sigma=metmodel_sigmaY[fUnclParamSwitch][0]+metmodel_sigmaY[fUnclParamSwitch][1]*sqrt(sumEt)-
    metmodel_sigmaY[1-fUnclParamSwitch][0]-metmodel_sigmaY[1-fUnclParamSwitch][1]*sqrt(sumEt);
  double systX_mean=metmodel_meanX[fUnclParamSwitch][0]+metmodel_meanX[fUnclParamSwitch][1]*sumEt-
    metmodel_meanX[1-fUnclParamSwitch][0]-metmodel_meanX[1-fUnclParamSwitch][1]*sumEt;
  double systY_mean=metmodel_meanY[fUnclParamSwitch][0]+metmodel_meanY[fUnclParamSwitch][1]*sumEt-
    metmodel_meanY[1-fUnclParamSwitch][0]-metmodel_meanY[1-fUnclParamSwitch][1]*sumEt;

  if(sumEt<0.0) sumEt=0.0;
  //___________ X-axis
  dummypar=metmodel_sigmaX_er[fUnclParamSwitch][0]*metmodel_sigmaX_er[fUnclParamSwitch][0]
    + metmodel_sigmaX_er[fUnclParamSwitch][1]*metmodel_sigmaX_er[fUnclParamSwitch][1]*sumEt
    + 2.0*metmodel_sigmaXcorr[fUnclParamSwitch]*metmodel_sigmaX_er[fUnclParamSwitch][0]*metmodel_sigmaX_er[fUnclParamSwitch][1]*sqrt(sumEt);
  if(dummypar<0.0) dummypar=0.0;
  sigmaX_err=sqrt(dummypar+systX_sigma*systX_sigma);
  sigmaMx=metmodel_sigmaX[fUnclParamSwitch][0]+metmodel_sigmaX[fUnclParamSwitch][1]*sqrt(sumEt)+systcode*sigmaX_err;
  if(sigmaMx<0.0) sigmaMx=0.0;

  dummypar=metmodel_meanX_er[fUnclParamSwitch][0]*metmodel_meanX_er[fUnclParamSwitch][0]
    + metmodel_meanX_er[fUnclParamSwitch][1]*metmodel_meanX_er[fUnclParamSwitch][1]*sumEt*sumEt
    + 2.0*metmodel_meanXcorr[fUnclParamSwitch]*metmodel_meanX_er[fUnclParamSwitch][0]*metmodel_meanX_er[fUnclParamSwitch][1]*sumEt;
  if(dummypar<0.0) dummypar=0.0;
  meanX_err=sqrt(dummypar+systX_mean*systX_mean);
  meanMx=metmodel_meanX[fUnclParamSwitch][0]+metmodel_meanX[fUnclParamSwitch][1]*sumEt+systcode*meanX_err;

  //___________ Y-axis
  dummypar=metmodel_sigmaY_er[fUnclParamSwitch][0]*metmodel_sigmaY_er[fUnclParamSwitch][0]
    + metmodel_sigmaY_er[fUnclParamSwitch][1]*metmodel_sigmaY_er[fUnclParamSwitch][1]*sumEt
    + 2.0*metmodel_sigmaYcorr[fUnclParamSwitch]*metmodel_sigmaY_er[fUnclParamSwitch][0]*metmodel_sigmaY_er[fUnclParamSwitch][1]*sqrt(sumEt);
  if(dummypar<0.0) dummypar=0.0;
  sigmaY_err=sqrt(dummypar+systY_sigma*systY_sigma);
  sigmaMy=metmodel_sigmaY[fUnclParamSwitch][0]+metmodel_sigmaY[fUnclParamSwitch][1]*sqrt(sumEt)+systcode*sigmaY_err;
  if(sigmaMy<0.0) sigmaMy=0.0;

  dummypar=metmodel_meanY_er[fUnclParamSwitch][0]*metmodel_meanY_er[fUnclParamSwitch][0]
    + metmodel_meanY_er[fUnclParamSwitch][1]*metmodel_meanY_er[fUnclParamSwitch][1]*sumEt*sumEt
    + 2.0*metmodel_meanYcorr[fUnclParamSwitch]*metmodel_meanY_er[fUnclParamSwitch][0]*metmodel_meanY_er[fUnclParamSwitch][1]*sumEt;
  if(dummypar<0.0) dummypar=0.0;
  meanY_err=sqrt(dummypar+systY_mean*systY_mean);  
  meanMy=metmodel_meanY[fUnclParamSwitch][0]+metmodel_meanY[fUnclParamSwitch][1]*sumEt+systcode*meanY_err;

  return;
}
/*}}}*/

//--------------------------- 09/20/06 ----------------------------------
//____ This is new function to generate Met due to EM object energy resolution
void TMyJetFilterModule::GenerateMyEMobjMet(CommonStuff miscstuff, int systcode,
						TVector2 &myEMGenMet)
/*{{{*/
{
  // Resolution: G/Et=sqrt(p0*p0/Et+p1*p1)
  double p0[2]; // resolution params: [0]=CEM; [1]=PEM
  double p1[2];
  double p0_er[2]; 
  double p1_er[2];
  //______________ CEM parametrization (temporary)
  p0[0]=0.135; // according to hep-ex/0510047 (jes nim)
  p1[0]=0.015;  // according to hep-ex/0510047 (jes nim)
  p0_er[0]=0.01*p0[0]; // ? guesstimate 
  p1_er[0]=0.01*p1[0]; // ? guesstimate 
  //______________ PEM parametrization (temporary)
  p0[1]=0.16; // according to hep-ex/0510047 (jes nim)
  p1[1]=0.01;  // according to hep-ex/0510047 (jes nim)
  p0_er[1]=0.01*p0[1]; // ? guesstimate 
  p1_er[1]=0.01*p1[1]; // ? guesstimate 
  //--------------------------------------
  int det_code=0;
  TLorentzVector emmetsum(0.0,0.0,0.0,0.0);
  TLorentzVector emmetsum_devm(0.0,0.0,0.0,0.0);
  TLorentzVector emmetsum_devp(0.0,0.0,0.0,0.0);
  for(int i=0; i<miscstuff.myRawPhoton.size(); i++)
    {
      if(fabs(miscstuff.myPhoEtaDet[i])<=1.1) det_code=0;
      else det_code=1;
      double obj_Pt=miscstuff.myCorrPhoton[i].Pt();
      double P_0=(p0[det_code]+systcode*p0_er[det_code])*(p0[det_code]+systcode*p0_er[det_code]);
      double P_1=(p1[det_code]+systcode*p1_er[det_code])*(p1[det_code]+systcode*p1_er[det_code]);
      double sigmaEM=sqrt(P_0/obj_Pt+P_1);
      double scale_factor=gRandom->Gaus(1.0,sigmaEM); // get jet energy scale factor
      if(miscstuff.myRawElectron.size()<2)
	{
	  if(i<2) while(scale_factor*miscstuff.myCorrPhoton[i].Pt()<13.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
	  if(i>=2) while(scale_factor*miscstuff.myCorrPhoton[i].Pt()<7.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
        }
      if(miscstuff.myRawElectron.size()>=2) while(scale_factor*miscstuff.myCorrPhoton[i].Pt()<7.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
      emmetsum=emmetsum+(1.0-scale_factor)*miscstuff.myCorrPhoton[i];
    }
  for(int i=0; i<miscstuff.myRawElectron.size(); i++)
    {
      if(fabs(miscstuff.myEleEtaDet[i])<=1.1) det_code=0;
      else det_code=1;
      double obj_Pt=miscstuff.myCorrElectron[i].Pt();
      double P_0=(p0[det_code]+systcode*p0_er[det_code])*(p0[det_code]+systcode*p0_er[det_code]);
      double P_1=(p1[det_code]+systcode*p1_er[det_code])*(p1[det_code]+systcode*p1_er[det_code]);
      double sigmaEM=sqrt(P_0/obj_Pt+P_1);
      double scale_factor=gRandom->Gaus(1.0,sigmaEM); // get jet energy scale factor
      if(miscstuff.myRawElectron.size()<2) while(scale_factor*miscstuff.myCorrElectron[i].Pt()<13.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
      if(miscstuff.myRawElectron.size()>=2)
	{
	  if(i==0) while(scale_factor*miscstuff.myCorrElectron[i].Pt()<20.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
	  if(i>0) while(scale_factor*miscstuff.myCorrElectron[i].Pt()<10.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
	}
      emmetsum=emmetsum+(1.0-scale_factor)*miscstuff.myCorrElectron[i];
    }
  myEMGenMet.Set(emmetsum.Px(),emmetsum.Py());

  return;
}
/*}}}*/


//_________ This function adds Met to closest jet, the output is "new" JetLev6
void TMyJetFilterModule::MyNewJetWithMet(JetStuff jetstuff, CommonStuff miscstuff,
						int jetcone, int jes_systcode,
						std::vector<TLorentzVector>& newJet)
/*{{{*/
{
  std::vector<double> dphi_jetmet;
  dphi_jetmet.clear();
  double phi_met=jetstuff.myMETcorr_th15.Phi();
  int case_switchJ=0; // for jets not matched to EM objects
  int closejet_indJ=-1; // for jets not matched to EM objects
  int NrawJet=0;
  for(int i=0; i<jetstuff.Jet_lev6_noEMobj.size(); i++) // packing "newJet" with "old" jets 
    {
      newJet.push_back(jetstuff.Jet_lev6_noEMobj[i]);
//       double jes=GetMyJetResolution(jetstuff.Jet_lev6_noEMobj[i].Pt(),jetstuff.EtaDetCorr[i],jetcone,jes_systcode); // get energy resolution
      NrawJet++;
      if(jetstuff.myMETcorr_th15.Mod()>0.0) 
	{
	  double phi_jet=jetstuff.Jet_lev6_noEMobj[i].Phi();
	  double dphi=fabs(TVector2::Phi_mpi_pi(phi_jet-phi_met));
	  dphi_jetmet.push_back(dphi);
	  double metjetRatio=jetstuff.myMETcorr_th15.Mod()*cos(dphi)/jetstuff.Jet_lev6_noEMobj[i].Pt();
	  double jetRatio=jetstuff.Jet_raw_noEMobj[i].Pt()/jetstuff.Jet_lev6_noEMobj[i].Pt();

	  double low_lim1=-0.8; // metjetRatio>-0.8 --- from Et(true)/Et(lev6) studies
	  double low_lim2;
	  double up_lim=3.0; // metjetRatio<2.5 --- from Et(true)/Et(lev6) studies
	  if(jetstuff.Jet_lev6_noEMobj[i].Pt()<15.0) low_lim2=-1.0*jetRatio;
	  else low_lim2=-1.0;
	  if((dphi<0.4 || (TMath::Pi()-dphi)<0.4) 
	     && (jetstuff.Npho_match[i]+jetstuff.Nele_match[i])==0 
	     && metjetRatio>low_lim1 && metjetRatio>low_lim2 && metjetRatio<up_lim) 
	    {
	      case_switchJ++;
	      closejet_indJ=i;
	    }
	}
    }
  if(dphi_jetmet.size()!=jetstuff.Jet_lev6_noEMobj.size() || NrawJet==0) return; 
  if(NrawJet==1 && case_switchJ==0) return;
  if(case_switchJ>1) return; // to be removed if case-2 is considered
  if(case_switchJ==1 && jetstuff.Jet_lev6_noEMobj[closejet_indJ].Pt()>15.0) // CASE-1 for unmatched jets only
    {
      double phi_jet=jetstuff.Jet_lev6_noEMobj[closejet_indJ].Phi();
      double eta_jet=jetstuff.Jet_lev6_noEMobj[closejet_indJ].Eta();
      double pt_jet=jetstuff.Jet_lev6_noEMobj[closejet_indJ].Pt();
      double e_jet=jetstuff.Jet_lev6_noEMobj[closejet_indJ].E();
      double dphi=fabs(TVector2::Phi_mpi_pi(phi_jet-phi_met));
      double pt_scale=(jetstuff.myMETcorr_th15.Mod()*cos(dphi)+pt_jet)/pt_jet;
      if(pt_scale>0.0) newJet[closejet_indJ].SetPtEtaPhiE(pt_jet*pt_scale,eta_jet,phi_jet,e_jet*pt_scale);
    }
  return;
}
/*}}}*/


//___ This is new function to generate Met; it takes care of unclustered & jet components of Met.
//___ cleanup histograms for generated met are also filled here.
void TMyJetFilterModule::GenerateMyTotalMet(JetStuff jetstuff,
						CommonStuff miscstuff, MetCleanupHisto_t& Hist, int jetcone,
						int systcode, int systcodeUncl, int metcode,
						TVector2 &myGenMet)
/*{{{*/
{
  int bad_met=1; // new as of 04/04/07
  std::vector<TLorentzVector> jet_smeared; // new as of 04/03/07
  jet_smeared.clear(); // new as of 04/03/07
  myGenMet.Set(0.0,0.0);
  int it_counter=0;
  it_counter++;
  //___________ generating met due to jets
  double SumEtJet_smear=0.0; // U.E.+M.I. for jets above threshold after smearing
  double SumEtJet=0.0;       // U.E.+M.I. for jets above threshold before smearing 
  double SumEtRawJet=0.0;
  double SumEtRawJet_smear=0.0;
  double Pt_cut[5]={2000.0,5.0,10.0,15.0,20.0};
  std::vector<int> jet_ind; // index of smeared jets fluctuated above Pt-cut
  double cut;
  double scale_factor=1.0;
  double sigma;
  TLorentzVector jetmetsum(0.0,0.0,0.0,0.0);
  jet_smeared.clear(); // new as of 04/03/07, need to clear again
  //_______________ new part as of 11/02/06
  std::vector<TLorentzVector> newJetLev6; // vector of new jets after adding Met and Jet
  newJetLev6.clear(); // making sure it's empty in the beginning
  MyNewJetWithMet(jetstuff,miscstuff,jetcone,systcode,newJetLev6); // adding jet and Met
  //_______________ end of new part 
  TVector2 myJetMet(0.0,0.0);
  if(metcode>0 && metcode<5) 
    {
      cut=Pt_cut[metcode];
      for(int i=0; i<newJetLev6.size(); i++) // "jetstuff.Jet_lev6_noEMobj" is replaced by "newJetLev6" (11/02/06)
	{
	  if(newJetLev6[i].Pt()>0.0 && fabs(jetstuff.EtaDetCorr[i])<=3.6) // "jetstuff.Jet_lev6_noEMobj" is replaced by "newJetLev6" (11/02/06) 
	    {
	      sigma=GetMyJetResolution(newJetLev6[i].Pt(),jetstuff.EtaDetCorr[i],jetcone,systcode); // get energy resolution
	      // "jetstuff.Jet_lev6_noEMobj" is replaced by "newJetLev6" on above line (11/02/06)
	      scale_factor=gRandom->Gaus(1.0,sigma); // get jet energy scale factor
	      while(scale_factor<0.2) // 0.2 is motivated by Et(had)/Et(det) studies (added on 04/30/07)
		{
		  scale_factor=gRandom->Gaus(1.0,sigma); 
		}
	      jet_smeared.push_back(scale_factor*newJetLev6[i]); // new as of 04/03/07
	      if((scale_factor*newJetLev6[i]).Pt()>cut) // "jetstuff.Jet_lev6_noEMobj" is replaced by "newJetLev6" (11/02/06)
		{
		  jet_ind.push_back(i); // temporary; commented out on 06/27/06 for studies
		  jetmetsum=jetmetsum+(1.0-scale_factor)*newJetLev6[i]; // "jetstuff.Jet_lev6_noEMobj" is replaced by "newJetLev6" (11/02/06)
		  // SumEtJet_smear= Sum(U.E.+M.I.)
		  SumEtJet_smear=SumEtJet_smear+jetstuff.Jet_lev5_noEMobj[i].Pt()
		    +jetstuff.Jet_lev1_noEMobj[i].Pt()
		    -jetstuff.Jet_lev4_noEMobj[i].Pt()
		    -jetstuff.Jet_lev6_noEMobj[i].Pt(); // Jet_lev6_noEMobj added on 06/26/06
		  SumEtRawJet_smear=SumEtRawJet_smear+jetstuff.Jet_raw_noEMobj[i].Pt();
		}
	      if(jetstuff.Jet_lev6_noEMobj[i].Pt()>cut) 
		{
		  SumEtJet=SumEtJet+jetstuff.Jet_lev5_noEMobj[i].Pt()
		    +jetstuff.Jet_lev1_noEMobj[i].Pt()
		    -jetstuff.Jet_lev4_noEMobj[i].Pt()
		    -jetstuff.Jet_lev6_noEMobj[i].Pt(); // Jet_lev6_noEMobj added on 06/26/06
		  SumEtRawJet=SumEtRawJet+jetstuff.Jet_raw_noEMobj[i].Pt();
		}
	    }
	}
      myJetMet.Set(jetmetsum.Px(),jetmetsum.Py());   
    }
  //_____________ generating unclustered component of Met
  double x;
  double y;
  double _meanX=0.0;
  double _sigmaX=0.0;
  double _meanY=0.0;
  double _sigmaY=0.0;
  double sumEt=0.0;
  double sumEt_raw=0.0;  
  //--------------------- this chunk of code is modified on 11/09/06
  if(metcode==1) sumEt_raw=jetstuff.mySumEtCorr_th5-SumEtJet+SumEtRawJet;
  if(metcode==2) sumEt_raw=jetstuff.mySumEtCorr_th10-SumEtJet+SumEtRawJet;
  if(metcode==3) sumEt_raw=jetstuff.mySumEtCorr_th15-SumEtJet+SumEtRawJet;
  if(metcode==4) sumEt_raw=jetstuff.mySumEtCorr_th20-SumEtJet+SumEtRawJet;
  if(metcode<1 || metcode>4) sumEt_raw=miscstuff.mySumEt_raw;
  sumEt=sumEt_raw+SumEtJet_smear-SumEtRawJet_smear;
  //------------------------------------------------ end of 11/09/06 part
  if(sumEt<0.0) sumEt=0.0;
  GetMyUnclMetResolution(sumEt,systcodeUncl,_sigmaX,_sigmaY,_meanX,_meanY);
  x=gRandom->Gaus(_meanX,_sigmaX);
  y=gRandom->Gaus(_meanY,_sigmaY);
  myGenMet.Set(x,y); 
  myGenMet=myGenMet+myJetMet;
  //_______________ new part, 09/20/06
  TVector2 myEMGenMet(0.0,0.0);
  GenerateMyEMobjMet(miscstuff,systcode,myEMGenMet);
  myGenMet=myGenMet+myEMGenMet;
  //------- end of new part, 09/20/06 
  bad_met=MyMetCleanUpCut(jetstuff,miscstuff,jet_smeared,myGenMet); // new as of 04/04/07  
  //___________________________________ filling GenMet cleanup histograms
  FillGenMetCleanupHistograms(Hist,jetstuff,miscstuff,myGenMet,jetstuff.Jet_lev6_noEMobj); // histo for generated met cleanup studies

  return;
}
/*}}}*/


//_______________ ZEROes met_results arrays; to be called in BeginJob
void TMyJetFilterModule::CleanMetResults(MetResults &metstuff)
/*{{{*/
{
  
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<10; j++)
	{
	  metstuff.Met20Njet_dt[i][j]=0; 
	  metstuff.Met25Njet_dt[i][j]=0; 
	  metstuff.Met30Njet_dt[i][j]=0; 
	  metstuff.Met35Njet_dt[i][j]=0; 
	  metstuff.Met40Njet_dt[i][j]=0; 
	  metstuff.Met45Njet_dt[i][j]=0; 
	  metstuff.Met50Njet_dt[i][j]=0; 
	  metstuff.Met75Njet_dt[i][j]=0; 
	  metstuff.Met100Njet_dt[i][j]=0; 
	  metstuff.Met150Njet_dt[i][j]=0; 
	  metstuff.Met20Njet_def[i][j]=0; 
	  metstuff.Met25Njet_def[i][j]=0; 
	  metstuff.Met30Njet_def[i][j]=0; 
	  metstuff.Met35Njet_def[i][j]=0; 
	  metstuff.Met40Njet_def[i][j]=0; 
	  metstuff.Met45Njet_def[i][j]=0; 
	  metstuff.Met50Njet_def[i][j]=0;
	  metstuff.Met75Njet_def[i][j]=0; 
	  metstuff.Met100Njet_def[i][j]=0; 
	  metstuff.Met150Njet_def[i][j]=0;

	  metstuff.Met20Njet_dvm[i][j]=0; 
	  metstuff.Met25Njet_dvm[i][j]=0; 
	  metstuff.Met30Njet_dvm[i][j]=0; 
	  metstuff.Met35Njet_dvm[i][j]=0; 
	  metstuff.Met40Njet_dvm[i][j]=0; 
	  metstuff.Met45Njet_dvm[i][j]=0; 
	  metstuff.Met50Njet_dvm[i][j]=0; 
	  metstuff.Met75Njet_dvm[i][j]=0; 
	  metstuff.Met100Njet_dvm[i][j]=0; 
	  metstuff.Met150Njet_dvm[i][j]=0; 
	  metstuff.Met20Njet_dvp[i][j]=0; 
	  metstuff.Met25Njet_dvp[i][j]=0; 
	  metstuff.Met30Njet_dvp[i][j]=0; 
	  metstuff.Met35Njet_dvp[i][j]=0; 
	  metstuff.Met40Njet_dvp[i][j]=0; 
	  metstuff.Met45Njet_dvp[i][j]=0; 
	  metstuff.Met50Njet_dvp[i][j]=0; 
	  metstuff.Met75Njet_dvp[i][j]=0; 
	  metstuff.Met100Njet_dvp[i][j]=0; 
	  metstuff.Met150Njet_dvp[i][j]=0; 

	  metstuff.Met20Njet_dvmUn[i][j]=0; 
	  metstuff.Met25Njet_dvmUn[i][j]=0; 
	  metstuff.Met30Njet_dvmUn[i][j]=0; 
	  metstuff.Met35Njet_dvmUn[i][j]=0; 
	  metstuff.Met40Njet_dvmUn[i][j]=0; 
	  metstuff.Met45Njet_dvmUn[i][j]=0; 
	  metstuff.Met50Njet_dvmUn[i][j]=0; 
	  metstuff.Met75Njet_dvmUn[i][j]=0; 
	  metstuff.Met100Njet_dvmUn[i][j]=0; 
	  metstuff.Met150Njet_dvmUn[i][j]=0; 
	  metstuff.Met20Njet_dvpUn[i][j]=0; 
	  metstuff.Met25Njet_dvpUn[i][j]=0; 
	  metstuff.Met30Njet_dvpUn[i][j]=0; 
	  metstuff.Met35Njet_dvpUn[i][j]=0; 
	  metstuff.Met40Njet_dvpUn[i][j]=0; 
	  metstuff.Met45Njet_dvpUn[i][j]=0; 
	  metstuff.Met50Njet_dvpUn[i][j]=0; 
	  metstuff.Met75Njet_dvpUn[i][j]=0; 
	  metstuff.Met100Njet_dvpUn[i][j]=0; 
	  metstuff.Met150Njet_dvpUn[i][j]=0; 

	  metstuff.Met20Njet_gen[i][j]=0.0; 
	  metstuff.Met25Njet_gen[i][j]=0.0; 
	  metstuff.Met30Njet_gen[i][j]=0.0; 
	  metstuff.Met35Njet_gen[i][j]=0.0; 
	  metstuff.Met40Njet_gen[i][j]=0.0; 
	  metstuff.Met45Njet_gen[i][j]=0.0; 
	  metstuff.Met50Njet_gen[i][j]=0.0; 
	  metstuff.Met75Njet_gen[i][j]=0.0; 
	  metstuff.Met100Njet_gen[i][j]=0.0; 
	  metstuff.Met150Njet_gen[i][j]=0.0; 
	  metstuff.Met20Njet_genstat[i][j]=0.0; 
	  metstuff.Met25Njet_genstat[i][j]=0.0; 
	  metstuff.Met30Njet_genstat[i][j]=0.0; 
	  metstuff.Met35Njet_genstat[i][j]=0.0; 
	  metstuff.Met40Njet_genstat[i][j]=0.0; 
	  metstuff.Met45Njet_genstat[i][j]=0.0; 
	  metstuff.Met50Njet_genstat[i][j]=0.0; 
	  metstuff.Met75Njet_genstat[i][j]=0.0; 
	  metstuff.Met100Njet_genstat[i][j]=0.0; 
	  metstuff.Met150Njet_genstat[i][j]=0.0; 
	  metstuff.Met20Njet_gensyst[i][j]=0.0; 
	  metstuff.Met25Njet_gensyst[i][j]=0.0; 
	  metstuff.Met30Njet_gensyst[i][j]=0.0; 
	  metstuff.Met35Njet_gensyst[i][j]=0.0; 
	  metstuff.Met40Njet_gensyst[i][j]=0.0; 
	  metstuff.Met45Njet_gensyst[i][j]=0.0; 
	  metstuff.Met50Njet_gensyst[i][j]=0.0;
	  metstuff.Met75Njet_gensyst[i][j]=0.0; 
	  metstuff.Met100Njet_gensyst[i][j]=0.0; 
	  metstuff.Met150Njet_gensyst[i][j]=0.0;
	} 
    }
  return;
}
/*}}}*/

//_____________________ calculates MetPDF and MetProb according to Clustered/Unclustered MetModel
void TMyJetFilterModule::MyMetPDF(int counter, TVector2 MetEvnt,
						TVector2 GenMet_def, TVector2 GenMet_dvm, TVector2 GenMet_dvp,
						MetProbHisto_t& Hist, int &event_code)
/*{{{*/
{
  if(fUseMetPDF>0 && Npoints>0)
    {
      if(counter==0) // clears global vectors at the beginning 
	{
	  MetGen_def.clear();
	  MetGen_dvm.clear();
	  MetGen_dvp.clear();  
	}
      //____fills global vectors with generated values
      MetGen_def.push_back(GenMet_def);
      MetGen_dvm.push_back(GenMet_dvm);
      MetGen_dvp.push_back(GenMet_dvp);
      if(counter==(Npoints-1)) // calculates MetPDFs and does all the job
	{
	  double toyMet=MetToy_min+(MetToy_max-MetToy_min)*(gRandom->Rndm());
	  double toyMetPhi=TMath::TwoPi()*(gRandom->Rndm());
	  double toyMet_x=toyMet*cos(toyMetPhi);
	  double toyMet_y=toyMet*sin(toyMetPhi);	  
	  TVector2 MetToy(toyMet_x,toyMet_y);
	  TVector2 MetToyEvnt(MetEvnt.Px()+toyMet_x,MetEvnt.Py()+toyMet_y);

	  int evnt_def=0;
	  int evnt_dvm=0;
	  int evnt_dvp=0;
	  int toy_def=0;
	  int toy_dvm=0;
	  int toy_dvp=0;
	  double ave2_def=0.0;
	  double ave2_dvm=0.0;
	  double ave2_dvp=0.0;
	  double ave_def=0.0;
	  double ave_dvm=0.0;
	  double ave_dvp=0.0;
	  double mean_def=0.0;
	  double mean_dvm=0.0;
	  double mean_dvp=0.0;
	  double sigma_def=0.0;
	  double sigma_dvm=0.0;
	  double sigma_dvp=0.0;
	  double lnProb_def=0.0;
	  double lnProb_dvm=0.0;
	  double lnProb_dvp=0.0;
	  double toy_lnProb_def=0.0;
	  double toy_lnProb_dvm=0.0;
	  double toy_lnProb_dvp=0.0;
	  for(int ind=0; ind<Npoints; ind++)
	    {
	      ave2_def+=MetGen_def[ind].Mod()*MetGen_def[ind].Mod();
	      ave2_dvm+=MetGen_dvm[ind].Mod()*MetGen_dvm[ind].Mod();
	      ave2_dvp+=MetGen_dvp[ind].Mod()*MetGen_dvp[ind].Mod();
	      ave_def+=MetGen_def[ind].Mod();
	      ave_dvm+=MetGen_dvm[ind].Mod();
	      ave_dvp+=MetGen_dvp[ind].Mod();
	      if(MetGen_def[ind].Mod()<MetEvnt.Mod()) evnt_def++;
	      if(MetGen_dvm[ind].Mod()<MetEvnt.Mod()) evnt_dvm++;
	      if(MetGen_dvp[ind].Mod()<MetEvnt.Mod()) evnt_dvp++;
	      if(MetGen_def[ind].Mod()<MetToyEvnt.Mod()) toy_def++;
	      if(MetGen_dvm[ind].Mod()<MetToyEvnt.Mod()) toy_dvm++;
	      if(MetGen_dvp[ind].Mod()<MetToyEvnt.Mod()) toy_dvp++;
	    }
	  mean_def=ave_def/Npoints;
	  mean_dvm=ave_dvm/Npoints;
	  mean_dvp=ave_dvp/Npoints;
	  sigma_def=sqrt(ave2_def/Npoints-mean_def*mean_def);
	  sigma_dvm=sqrt(ave2_dvm/Npoints-mean_dvm*mean_dvm);
	  sigma_dvp=sqrt(ave2_dvp/Npoints-mean_dvp*mean_dvp);
	  if(evnt_def<Npoints) lnProb_def=-1.0*log10(1.0-(1.0*evnt_def)/Npoints);
	  else lnProb_def=-1.0*log10(1.0/Npoints);
	  if(evnt_dvm<Npoints) lnProb_dvm=-1.0*log10(1.0-(1.0*evnt_dvm)/Npoints);
	  else lnProb_dvm=-1.0*log10(1.0/Npoints);
	  if(evnt_dvp<Npoints) lnProb_dvp=-1.0*log10(1.0-(1.0*evnt_dvp)/Npoints);
	  else lnProb_dvp=-1.0*log10(1.0/Npoints);

	  if(toy_def<Npoints) toy_lnProb_def=-1.0*log10(1.0-(1.0*toy_def)/Npoints);
	  else toy_lnProb_def=-1.0*log10(1.0/Npoints);
	  if(toy_dvm<Npoints) toy_lnProb_dvm=-1.0*log10(1.0-(1.0*toy_dvm)/Npoints);
	  else toy_lnProb_dvm=-1.0*log10(1.0/Npoints);
	  if(toy_dvp<Npoints) toy_lnProb_dvp=-1.0*log10(1.0-(1.0*toy_dvp)/Npoints);
	  else toy_lnProb_dvp=-1.0*log10(1.0/Npoints);

	  //______________ filling histograms
	  Hist.fDataMetProb_def->Fill(lnProb_def); 
	  Hist.fDataMetProb_dvm->Fill(lnProb_dvm);  
	  Hist.fDataMetProb_dvp->Fill(lnProb_dvp);  
	  Hist.fDataMetSigma_def->Fill((MetEvnt.Mod()-mean_def)/sigma_def);
	  Hist.fDataMetSigma_dvm->Fill((MetEvnt.Mod()-mean_dvm)/sigma_dvm);
	  Hist.fDataMetSigma_dvp->Fill((MetEvnt.Mod()-mean_dvp)/sigma_dvp);
	  Hist.fDataMetProb_Met_def->Fill(MetEvnt.Mod(),lnProb_def); 
	  Hist.fDataMetProb_Met_dvm->Fill(MetEvnt.Mod(),lnProb_dvm);  
	  Hist.fDataMetProb_Met_dvp->Fill(MetEvnt.Mod(),lnProb_dvp);  
	  Hist.fDataMetSigma_Met_def->Fill(MetEvnt.Mod(),(MetEvnt.Mod()-mean_def)/sigma_def);
	  Hist.fDataMetSigma_Met_dvm->Fill(MetEvnt.Mod(),(MetEvnt.Mod()-mean_dvm)/sigma_dvm);
	  Hist.fDataMetSigma_Met_dvp->Fill(MetEvnt.Mod(),(MetEvnt.Mod()-mean_dvp)/sigma_dvp);
	  Hist.fDataMet_vs_GenSigmaMet->Fill(sigma_def,MetEvnt.Mod());
	  
	  //______________ filling histograms for toy MC
	  Hist.fToyEvntMet_vs_GenSigmaMet->Fill(sigma_def,MetToyEvnt.Mod());
	  Hist.fToyMetProb_toyMet_def->Fill(MetToy.Mod(),toy_lnProb_def); 
	  Hist.fToyMetProb_toyMet_dvm->Fill(MetToy.Mod(),toy_lnProb_dvm); 
	  Hist.fToyMetProb_toyMet_dvp->Fill(MetToy.Mod(),toy_lnProb_dvp); 
	  Hist.fToyMetSigma_toyMet_def->Fill(MetToy.Mod(),(MetToyEvnt.Mod()-mean_def)/sigma_def);
	  Hist.fToyMetSigma_toyMet_dvm->Fill(MetToy.Mod(),(MetToyEvnt.Mod()-mean_dvm)/sigma_dvm);
	  Hist.fToyMetSigma_toyMet_dvp->Fill(MetToy.Mod(),(MetToyEvnt.Mod()-mean_dvp)/sigma_dvp);	  
	  for(int ic=0; ic<20; ic++)
	    {
	      if(MetToy.Mod()>(MetToy_min+(MetToy_max-MetToy_min)*ic/20.0) && MetToy.Mod()<=(MetToy_min+(MetToy_max-MetToy_min)*(ic+1)/20.0))
		{
		  Hist.fToyMetProb_def[ic]->Fill(toy_lnProb_def); 
		  Hist.fToyMetProb_dvm[ic]->Fill(toy_lnProb_dvm); 
		  Hist.fToyMetProb_dvp[ic]->Fill(toy_lnProb_dvp); 
		  Hist.fToyMetSigma_def[ic]->Fill((MetToyEvnt.Mod()-mean_def)/sigma_def);
		  Hist.fToyMetSigma_dvm[ic]->Fill((MetToyEvnt.Mod()-mean_dvm)/sigma_dvm);
		  Hist.fToyMetSigma_dvp[ic]->Fill((MetToyEvnt.Mod()-mean_dvp)/sigma_dvp);
		  break;
		}
	    }
	  //_________ setting event status code
	  if(fUseMetPDF==1 && lnProb_def>2.56) event_code=1; // 3*sigma significance
	  if(fUseMetPDF==2 && lnProb_def>3.0) event_code=1; // 3.29*sigma significance
	  if(fUseMetPDF==3 && lnProb_def>4.0) event_code=1; // 3.89*sigma significance
	  if(fUseMetPDF==4 && lnProb_def>4.2) event_code=1; // 4.0*sigma significance
	  if(fUseMetPDF==5 && lnProb_def>2.0) event_code=1; // 2.58*sigma significance

	  //____________________ selection according to probability
	  if(lnProb_def>2.56) Hist.fDataMet_highLnProb->Fill(MetEvnt.Mod()); // filling histo for low probability events (def parametrization) 
	  else Hist.fDataMet_lowLnProb->Fill(MetEvnt.Mod()); // filling histo for high probability events (def parametrization)
	  //____________________ selection according to sigma
	  if((MetEvnt.Mod()-mean_def)/sigma_def>3.0) Hist.fDataMet_highSigma->Fill(MetEvnt.Mod()); // filling histo for events fluctuated above 3sigma (def parametrization) 
	  else Hist.fDataMet_lowSigma->Fill(MetEvnt.Mod()); // filling histo for events fluctuated below 3sigma (def parametrization)

	  for(int ind=0; ind<Npoints; ind++)
	    {
	      //____________________ selection according to probability
	      if(lnProb_def>2.56) // filling histo for low probability events (def parametrization) 
		{
		  Hist.fGenMet_def_highLnProb->Fill(MetGen_def[ind].Mod()); 
		  Hist.fGenMet_dvm_highLnProb->Fill(MetGen_dvm[ind].Mod()); 
		  Hist.fGenMet_dvp_highLnProb->Fill(MetGen_dvp[ind].Mod()); 
		}
	      else  // filling histo for high probability events (def parametrization)
		{
		  Hist.fGenMet_def_lowLnProb->Fill(MetGen_def[ind].Mod()); 
		  Hist.fGenMet_dvm_lowLnProb->Fill(MetGen_dvm[ind].Mod()); 
		  Hist.fGenMet_dvp_lowLnProb->Fill(MetGen_dvp[ind].Mod()); 
		}
	      //____________________ selection according to sigma
	      if((MetEvnt.Mod()-mean_def)/sigma_def>3.0) // filling histo for events fluctuated above 3sigma (def parametrization) 
		{
		  Hist.fGenMet_def_highSigma->Fill(MetGen_def[ind].Mod()); 
		  Hist.fGenMet_dvm_highSigma->Fill(MetGen_dvm[ind].Mod()); 
		  Hist.fGenMet_dvp_highSigma->Fill(MetGen_dvp[ind].Mod()); 
		}
	      else  // filling histo for events fluctuated below 3sigma (def parametrization)
		{
		  Hist.fGenMet_def_lowSigma->Fill(MetGen_def[ind].Mod()); 
		  Hist.fGenMet_dvm_lowSigma->Fill(MetGen_dvm[ind].Mod()); 
		  Hist.fGenMet_dvp_lowSigma->Fill(MetGen_dvp[ind].Mod()); 
		} 
	    }
	}
    }
  return;
}
/*}}}*/


void TMyJetFilterModule::MyMetEventCount(JetStuff jetstuff,
						CommonStuff miscstuff, int metcode, MetResults &metstuff)
/*{{{*/
{
  double met=0.0;
  int njet15=0;
  int njet20=0;
  int njet25=0;
  
  if(metcode==1) met=jetstuff.myMETcorr_th5.Mod(); 
  if(metcode==2) met=jetstuff.myMETcorr_th10.Mod(); 
  if(metcode==3) met=jetstuff.myMETcorr_th15.Mod(); 
  if(metcode==4) met=jetstuff.myMETcorr_th20.Mod();
  if(metcode<1 || metcode>4) met=miscstuff.myMET_raw.Mod();
  njet15= (jetstuff.myNjet_th15 < 10) ? jetstuff.myNjet_th15 : 9; 
  njet20= (jetstuff.myNjet_th20 < 10) ? jetstuff.myNjet_th20 : 9; 
  njet25= (jetstuff.myNjet_th25 < 10) ? jetstuff.myNjet_th25 : 9; 
  
  if(met>20) 
    {
      metstuff.Met20Njet_dt[0][njet15]++;
      metstuff.Met20Njet_dt[1][njet20]++;
      metstuff.Met20Njet_dt[2][njet25]++;
    }
  if(met>25) 
    {
      metstuff.Met25Njet_dt[0][njet15]++;
      metstuff.Met25Njet_dt[1][njet20]++;
      metstuff.Met25Njet_dt[2][njet25]++;
    }
  if(met>30) 
    {
      metstuff.Met30Njet_dt[0][njet15]++;
      metstuff.Met30Njet_dt[1][njet20]++;
      metstuff.Met30Njet_dt[2][njet25]++;
    }
  if(met>35) 
    {
      metstuff.Met35Njet_dt[0][njet15]++;
      metstuff.Met35Njet_dt[1][njet20]++;
      metstuff.Met35Njet_dt[2][njet25]++;
    }
  if(met>40) 
    {
      metstuff.Met40Njet_dt[0][njet15]++;
      metstuff.Met40Njet_dt[1][njet20]++;
      metstuff.Met40Njet_dt[2][njet25]++;
    }
  if(met>45) 
    {
      metstuff.Met45Njet_dt[0][njet15]++;
      metstuff.Met45Njet_dt[1][njet20]++;
      metstuff.Met45Njet_dt[2][njet25]++;
    }
  if(met>50) 
    {
      metstuff.Met50Njet_dt[0][njet15]++;
      metstuff.Met50Njet_dt[1][njet20]++;
      metstuff.Met50Njet_dt[2][njet25]++;
    }
  if(met>75) 
    {
      metstuff.Met75Njet_dt[0][njet15]++;
      metstuff.Met75Njet_dt[1][njet20]++;
      metstuff.Met75Njet_dt[2][njet25]++;
    }
  if(met>100) 
    {
      metstuff.Met100Njet_dt[0][njet15]++;
      metstuff.Met100Njet_dt[1][njet20]++;
      metstuff.Met100Njet_dt[2][njet25]++;
    }
  if(met>150) 
    {
      metstuff.Met150Njet_dt[0][njet15]++;
      metstuff.Met150Njet_dt[1][njet20]++;
      metstuff.Met150Njet_dt[2][njet25]++;
    }
  return;
}

/*}}}*/


//______________ counts generated events above Met cut
void TMyJetFilterModule::MyMetModelEventCount(JetStuff jetstuff,
							TVector2 myGenMet_def, TVector2 myGenMet_dvm, 
					      TVector2 myGenMet_dvp, TVector2 myGenMet_dvmUn, 
					      TVector2 myGenMet_dvpUn, MetResults &metstuff)
/*{{{*/
{
  int njet15=0;
  int njet20=0;
  int njet25=0;
  
  double met_def=myGenMet_def.Mod();
  double met_dvm=myGenMet_dvm.Mod();
  double met_dvp=myGenMet_dvp.Mod();
  double met_dvmUn=myGenMet_dvmUn.Mod();
  double met_dvpUn=myGenMet_dvpUn.Mod();
  
  njet15= (jetstuff.myNjet_th15 < 10) ? jetstuff.myNjet_th15 : 9; 
  njet20= (jetstuff.myNjet_th20 < 10) ? jetstuff.myNjet_th20 : 9; 
  njet25= (jetstuff.myNjet_th20 < 10) ? jetstuff.myNjet_th25 : 9; 
  
  //_____________ default parametrization
  if(met_def>20.0 && met_def<2000.0) 
    {
      metstuff.Met20Njet_def[0][njet15]++;
      metstuff.Met20Njet_def[1][njet20]++;
      metstuff.Met20Njet_def[2][njet25]++;
    }
  if(met_def>25.0 && met_def<2000.0) 
    {
      metstuff.Met25Njet_def[0][njet15]++;
      metstuff.Met25Njet_def[1][njet20]++;
      metstuff.Met25Njet_def[2][njet25]++;
    }
  if(met_def>30.0 && met_def<2000.0) 
    {
      metstuff.Met30Njet_def[0][njet15]++;
      metstuff.Met30Njet_def[1][njet20]++;
      metstuff.Met30Njet_def[2][njet25]++;
    }
  if(met_def>35.0 && met_def<2000.0) 
    {
      metstuff.Met35Njet_def[0][njet15]++;
      metstuff.Met35Njet_def[1][njet20]++;
      metstuff.Met35Njet_def[2][njet25]++;
    }
  if(met_def>40.0 && met_def<2000.0) 
    {
      metstuff.Met40Njet_def[0][njet15]++;
      metstuff.Met40Njet_def[1][njet20]++;
      metstuff.Met40Njet_def[2][njet25]++;
    }
  if(met_def>45.0 && met_def<2000.0) 
    {
      metstuff.Met45Njet_def[0][njet15]++;
      metstuff.Met45Njet_def[1][njet20]++;
      metstuff.Met45Njet_def[2][njet25]++;
    }
  if(met_def>50.0 && met_def<2000.0) 
    {
      metstuff.Met50Njet_def[0][njet15]++;
      metstuff.Met50Njet_def[1][njet20]++;
      metstuff.Met50Njet_def[2][njet25]++;
    }
  if(met_def>75.0 && met_def<2000.0) 
    {
      metstuff.Met75Njet_def[0][njet15]++;
      metstuff.Met75Njet_def[1][njet20]++;
      metstuff.Met75Njet_def[2][njet25]++;
    }
  if(met_def>100.0 && met_def<2000.0) 
    {
      metstuff.Met100Njet_def[0][njet15]++;
      metstuff.Met100Njet_def[1][njet20]++;
      metstuff.Met100Njet_def[2][njet25]++;
    }
  if(met_def>150.0 && met_def<2000.0) 
    {
      metstuff.Met150Njet_def[0][njet15]++;
      metstuff.Met150Njet_def[1][njet20]++;
      metstuff.Met150Njet_def[2][njet25]++;
    }
  //_____________ -SigmaJet parametrization
  if(met_dvm>20.0 && met_dvm<2000.0) 
    {
      metstuff.Met20Njet_dvm[0][njet15]++;
      metstuff.Met20Njet_dvm[1][njet20]++;
      metstuff.Met20Njet_dvm[2][njet25]++;
    }
  if(met_dvm>25.0 && met_dvm<2000.0) 
    {
      metstuff.Met25Njet_dvm[0][njet15]++;
      metstuff.Met25Njet_dvm[1][njet20]++;
      metstuff.Met25Njet_dvm[2][njet25]++;
    }
  if(met_dvm>30.0 && met_dvm<2000.0) 
    {
      metstuff.Met30Njet_dvm[0][njet15]++;
      metstuff.Met30Njet_dvm[1][njet20]++;
      metstuff.Met30Njet_dvm[2][njet25]++;
    }
  if(met_dvm>35.0 && met_dvm<2000.0) 
    {
      metstuff.Met35Njet_dvm[0][njet15]++;
      metstuff.Met35Njet_dvm[1][njet20]++;
      metstuff.Met35Njet_dvm[2][njet25]++;
    }
  if(met_dvm>40.0 && met_dvm<2000.0) 
    {
      metstuff.Met40Njet_dvm[0][njet15]++;
      metstuff.Met40Njet_dvm[1][njet20]++;
      metstuff.Met40Njet_dvm[2][njet25]++;
    }
  if(met_dvm>45.0 && met_dvm<2000.0) 
    {
      metstuff.Met45Njet_dvm[0][njet15]++;
      metstuff.Met45Njet_dvm[1][njet20]++;
      metstuff.Met45Njet_dvm[2][njet25]++;
    }
  if(met_dvm>50.0 && met_dvm<2000.0) 
    {
      metstuff.Met50Njet_dvm[0][njet15]++;
      metstuff.Met50Njet_dvm[1][njet20]++;
      metstuff.Met50Njet_dvm[2][njet25]++;
    }
  if(met_dvm>75.0 && met_dvm<2000.0) 
    {
      metstuff.Met75Njet_dvm[0][njet15]++;
      metstuff.Met75Njet_dvm[1][njet20]++;
      metstuff.Met75Njet_dvm[2][njet25]++;
    }
  if(met_dvm>100.0 && met_dvm<2000.0) 
    {
      metstuff.Met100Njet_dvm[0][njet15]++;
      metstuff.Met100Njet_dvm[1][njet20]++;
      metstuff.Met100Njet_dvm[2][njet25]++;
    }
  if(met_dvm>150.0 && met_dvm<2000.0) 
    {
      metstuff.Met150Njet_dvm[0][njet15]++;
      metstuff.Met150Njet_dvm[1][njet20]++;
      metstuff.Met150Njet_dvm[2][njet25]++;
    }
  //_____________ +SigmaJet parametrization
  if(met_dvp>20.0 && met_dvp<2000.0) 
    {
      metstuff.Met20Njet_dvp[0][njet15]++;
      metstuff.Met20Njet_dvp[1][njet20]++;
      metstuff.Met20Njet_dvp[2][njet25]++;
    }
  if(met_dvp>25.0 && met_dvp<2000.0) 
    {
      metstuff.Met25Njet_dvp[0][njet15]++;
      metstuff.Met25Njet_dvp[1][njet20]++;
      metstuff.Met25Njet_dvp[2][njet25]++;
    }
  if(met_dvp>30.0 && met_dvp<2000.0) 
    {
      metstuff.Met30Njet_dvp[0][njet15]++;
      metstuff.Met30Njet_dvp[1][njet20]++;
      metstuff.Met30Njet_dvp[2][njet25]++;
    }
  if(met_dvp>35.0 && met_dvp<2000.0) 
    {
      metstuff.Met35Njet_dvp[0][njet15]++;
      metstuff.Met35Njet_dvp[1][njet20]++;
      metstuff.Met35Njet_dvp[2][njet25]++;
    }
  if(met_dvp>40.0 && met_dvp<2000.0) 
    {
      metstuff.Met40Njet_dvp[0][njet15]++;
      metstuff.Met40Njet_dvp[1][njet20]++;
      metstuff.Met40Njet_dvp[2][njet25]++;
    }
  if(met_dvp>45.0 && met_dvp<2000.0) 
    {
      metstuff.Met45Njet_dvp[0][njet15]++;
      metstuff.Met45Njet_dvp[1][njet20]++;
      metstuff.Met45Njet_dvp[2][njet25]++;
    }
  if(met_dvp>50.0 && met_dvp<2000.0) 
    {
      metstuff.Met50Njet_dvp[0][njet15]++;
      metstuff.Met50Njet_dvp[1][njet20]++;
      metstuff.Met50Njet_dvp[2][njet25]++;
    }
  if(met_dvp>75.0 && met_dvp<2000.0) 
    {
      metstuff.Met75Njet_dvp[0][njet15]++;
      metstuff.Met75Njet_dvp[1][njet20]++;
      metstuff.Met75Njet_dvp[2][njet25]++;
    }
  if(met_dvp>100.0 && met_dvp<2000.0) 
    {
      metstuff.Met100Njet_dvp[0][njet15]++;
      metstuff.Met100Njet_dvp[1][njet20]++;
      metstuff.Met100Njet_dvp[2][njet25]++;
    }
  if(met_dvp>150.0 && met_dvp<2000.0) 
    {
      metstuff.Met150Njet_dvp[0][njet15]++;
      metstuff.Met150Njet_dvp[1][njet20]++;
      metstuff.Met150Njet_dvp[2][njet25]++;
    }
  //_____________ -SigmaUncl parametrization
  if(met_dvmUn>20.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met20Njet_dvmUn[0][njet15]++;
      metstuff.Met20Njet_dvmUn[1][njet20]++;
      metstuff.Met20Njet_dvmUn[2][njet25]++;
    }
  if(met_dvmUn>25.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met25Njet_dvmUn[0][njet15]++;
      metstuff.Met25Njet_dvmUn[1][njet20]++;
      metstuff.Met25Njet_dvmUn[2][njet25]++;
    }
  if(met_dvmUn>30.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met30Njet_dvmUn[0][njet15]++;
      metstuff.Met30Njet_dvmUn[1][njet20]++;
      metstuff.Met30Njet_dvmUn[2][njet25]++;
    }
  if(met_dvmUn>35.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met35Njet_dvmUn[0][njet15]++;
      metstuff.Met35Njet_dvmUn[1][njet20]++;
      metstuff.Met35Njet_dvmUn[2][njet25]++;
    }
  if(met_dvmUn>40.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met40Njet_dvmUn[0][njet15]++;
      metstuff.Met40Njet_dvmUn[1][njet20]++;
      metstuff.Met40Njet_dvmUn[2][njet25]++;
    }
  if(met_dvmUn>45.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met45Njet_dvmUn[0][njet15]++;
      metstuff.Met45Njet_dvmUn[1][njet20]++;
      metstuff.Met45Njet_dvmUn[2][njet25]++;
    }
  if(met_dvmUn>50.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met50Njet_dvmUn[0][njet15]++;
      metstuff.Met50Njet_dvmUn[1][njet20]++;
      metstuff.Met50Njet_dvmUn[2][njet25]++;
    }
  if(met_dvmUn>75.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met75Njet_dvmUn[0][njet15]++;
      metstuff.Met75Njet_dvmUn[1][njet20]++;
      metstuff.Met75Njet_dvmUn[2][njet25]++;
    }
  if(met_dvmUn>100.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met100Njet_dvmUn[0][njet15]++;
      metstuff.Met100Njet_dvmUn[1][njet20]++;
      metstuff.Met100Njet_dvmUn[2][njet25]++;
    }
  if(met_dvmUn>150.0 && met_dvmUn<2000.0) 
    {
      metstuff.Met150Njet_dvmUn[0][njet15]++;
      metstuff.Met150Njet_dvmUn[1][njet20]++;
      metstuff.Met150Njet_dvmUn[2][njet25]++;
    }
  //_____________ +SigmaUncl parametrization
  if(met_dvpUn>20.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met20Njet_dvpUn[0][njet15]++;
      metstuff.Met20Njet_dvpUn[1][njet20]++;
      metstuff.Met20Njet_dvpUn[2][njet25]++;
    }
  if(met_dvpUn>25.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met25Njet_dvpUn[0][njet15]++;
      metstuff.Met25Njet_dvpUn[1][njet20]++;
      metstuff.Met25Njet_dvpUn[2][njet25]++;
    }
  if(met_dvpUn>30.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met30Njet_dvpUn[0][njet15]++;
      metstuff.Met30Njet_dvpUn[1][njet20]++;
      metstuff.Met30Njet_dvpUn[2][njet25]++;
    }
  if(met_dvpUn>35.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met35Njet_dvpUn[0][njet15]++;
      metstuff.Met35Njet_dvpUn[1][njet20]++;
      metstuff.Met35Njet_dvpUn[2][njet25]++;
    }
  if(met_dvpUn>40.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met40Njet_dvpUn[0][njet15]++;
      metstuff.Met40Njet_dvpUn[1][njet20]++;
      metstuff.Met40Njet_dvpUn[2][njet25]++;
    }
  if(met_dvpUn>45.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met45Njet_dvpUn[0][njet15]++;
      metstuff.Met45Njet_dvpUn[1][njet20]++;
      metstuff.Met45Njet_dvpUn[2][njet25]++;
    }
  if(met_dvpUn>50.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met50Njet_dvpUn[0][njet15]++;
      metstuff.Met50Njet_dvpUn[1][njet20]++;
      metstuff.Met50Njet_dvpUn[2][njet25]++;
    }
  if(met_dvpUn>75.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met75Njet_dvpUn[0][njet15]++;
      metstuff.Met75Njet_dvpUn[1][njet20]++;
      metstuff.Met75Njet_dvpUn[2][njet25]++;
    }
  if(met_dvpUn>100.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met100Njet_dvpUn[0][njet15]++;
      metstuff.Met100Njet_dvpUn[1][njet20]++;
      metstuff.Met100Njet_dvpUn[2][njet25]++;
    }
  if(met_dvpUn>150.0 && met_dvpUn<2000.0) 
    {
      metstuff.Met150Njet_dvpUn[0][njet15]++;
      metstuff.Met150Njet_dvpUn[1][njet20]++;
      metstuff.Met150Njet_dvpUn[2][njet25]++;
    }
  return;
}

/*}}}*/


//______________ prints out a comparison of Met in Data to Met Model predictions
void TMyJetFilterModule::MyMetResults(MetResults &metstuff)
/*{{{*/
{

  if(Npoints>0)
    {
      double dummy1=0.0;
      double dummy2=0.0;
      double syst_J1=0.0;
      double syst_J2=0.0;
      double syst_U1=0.0;
      double syst_U2=0.0;
      for(int i=0; i<3; i++)
	{
	  for(int j=0; j<10; j++)
	    {
	      metstuff.Met20Njet_genstat[i][j]=sqrt(1.0*metstuff.Met20Njet_def[i][j])/Npoints; 
	      metstuff.Met25Njet_genstat[i][j]=sqrt(1.0*metstuff.Met25Njet_def[i][j])/Npoints; 
	      metstuff.Met30Njet_genstat[i][j]=sqrt(1.0*metstuff.Met30Njet_def[i][j])/Npoints; 
	      metstuff.Met35Njet_genstat[i][j]=sqrt(1.0*metstuff.Met35Njet_def[i][j])/Npoints; 
	      metstuff.Met40Njet_genstat[i][j]=sqrt(1.0*metstuff.Met40Njet_def[i][j])/Npoints; 
	      metstuff.Met45Njet_genstat[i][j]=sqrt(1.0*metstuff.Met45Njet_def[i][j])/Npoints; 
	      metstuff.Met50Njet_genstat[i][j]=sqrt(1.0*metstuff.Met50Njet_def[i][j])/Npoints;
	      metstuff.Met75Njet_genstat[i][j]=sqrt(1.0*metstuff.Met75Njet_def[i][j])/Npoints;
	      metstuff.Met100Njet_genstat[i][j]=sqrt(1.0*metstuff.Met100Njet_def[i][j])/Npoints;
	      metstuff.Met150Njet_genstat[i][j]=sqrt(1.0*metstuff.Met150Njet_def[i][j])/Npoints;

	      //_____Met=20
	      syst_J1=fabs(1.0*(metstuff.Met20Njet_def[i][j]-metstuff.Met20Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met20Njet_def[i][j]-metstuff.Met20Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met20Njet_def[i][j]-metstuff.Met20Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met20Njet_def[i][j]-metstuff.Met20Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met20Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;
	      //_____Met=25
	      syst_J1=fabs(1.0*(metstuff.Met25Njet_def[i][j]-metstuff.Met25Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met25Njet_def[i][j]-metstuff.Met25Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met25Njet_def[i][j]-metstuff.Met25Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met25Njet_def[i][j]-metstuff.Met25Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met25Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;
	      //_____Met=30
	      syst_J1=fabs(1.0*(metstuff.Met30Njet_def[i][j]-metstuff.Met30Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met30Njet_def[i][j]-metstuff.Met30Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met30Njet_def[i][j]-metstuff.Met30Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met30Njet_def[i][j]-metstuff.Met30Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met30Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;
	      //_____Met=35
	      syst_J1=fabs(1.0*(metstuff.Met35Njet_def[i][j]-metstuff.Met35Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met35Njet_def[i][j]-metstuff.Met35Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met35Njet_def[i][j]-metstuff.Met35Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met35Njet_def[i][j]-metstuff.Met35Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met35Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;
	      //_____Met=40
	      syst_J1=fabs(1.0*(metstuff.Met40Njet_def[i][j]-metstuff.Met40Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met40Njet_def[i][j]-metstuff.Met40Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met40Njet_def[i][j]-metstuff.Met40Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met40Njet_def[i][j]-metstuff.Met40Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met40Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;
	      //_____Met=45
	      syst_J1=fabs(1.0*(metstuff.Met45Njet_def[i][j]-metstuff.Met45Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met45Njet_def[i][j]-metstuff.Met45Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met45Njet_def[i][j]-metstuff.Met45Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met45Njet_def[i][j]-metstuff.Met45Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met45Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;
	      //_____Met=50
	      syst_J1=fabs(1.0*(metstuff.Met50Njet_def[i][j]-metstuff.Met50Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met50Njet_def[i][j]-metstuff.Met50Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met50Njet_def[i][j]-metstuff.Met50Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met50Njet_def[i][j]-metstuff.Met50Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met50Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;
	      //_____Met=75
	      syst_J1=fabs(1.0*(metstuff.Met75Njet_def[i][j]-metstuff.Met75Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met75Njet_def[i][j]-metstuff.Met75Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met75Njet_def[i][j]-metstuff.Met75Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met75Njet_def[i][j]-metstuff.Met75Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met75Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;
	      //_____Met=100
	      syst_J1=fabs(1.0*(metstuff.Met100Njet_def[i][j]-metstuff.Met100Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met100Njet_def[i][j]-metstuff.Met100Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met100Njet_def[i][j]-metstuff.Met100Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met100Njet_def[i][j]-metstuff.Met100Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met100Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;
	      //_____Met=150
	      syst_J1=fabs(1.0*(metstuff.Met150Njet_def[i][j]-metstuff.Met150Njet_dvm[i][j])/Npoints);
	      syst_J2=fabs(1.0*(metstuff.Met150Njet_def[i][j]-metstuff.Met150Njet_dvp[i][j])/Npoints);
	      syst_U1=fabs(1.0*(metstuff.Met150Njet_def[i][j]-metstuff.Met150Njet_dvmUn[i][j])/Npoints);
	      syst_U2=fabs(1.0*(metstuff.Met150Njet_def[i][j]-metstuff.Met150Njet_dvpUn[i][j])/Npoints);
	      dummy1=sqrt(syst_J1*syst_J1+syst_U1*syst_U1);
	      dummy2=sqrt(syst_J2*syst_J2+syst_U2*syst_U2);
	      metstuff.Met150Njet_gensyst[i][j]= (dummy1 > dummy2) ? dummy1 : dummy2;

	      metstuff.Met20Njet_gen[i][j]=1.0*metstuff.Met20Njet_def[i][j]/Npoints; 
	      metstuff.Met25Njet_gen[i][j]=1.0*metstuff.Met25Njet_def[i][j]/Npoints; 
	      metstuff.Met30Njet_gen[i][j]=1.0*metstuff.Met30Njet_def[i][j]/Npoints; 
	      metstuff.Met35Njet_gen[i][j]=1.0*metstuff.Met35Njet_def[i][j]/Npoints; 
	      metstuff.Met40Njet_gen[i][j]=1.0*metstuff.Met40Njet_def[i][j]/Npoints; 
	      metstuff.Met45Njet_gen[i][j]=1.0*metstuff.Met45Njet_def[i][j]/Npoints; 
	      metstuff.Met50Njet_gen[i][j]=1.0*metstuff.Met50Njet_def[i][j]/Npoints;
	      metstuff.Met75Njet_gen[i][j]=1.0*metstuff.Met75Njet_def[i][j]/Npoints;
	      metstuff.Met100Njet_gen[i][j]=1.0*metstuff.Met100Njet_def[i][j]/Npoints;
	      metstuff.Met150Njet_gen[i][j]=1.0*metstuff.Met150Njet_def[i][j]/Npoints;
	    }
	}
    }

  ofstream outfile(fDatFileName, std::ios::out);
  if (! outfile) 
    {
      std::cerr<<"Error in openning of .DAT file\n";
    }
  if(outfile)
    {
      outfile<<"---------------------------- Results ---------------------------"<<"\n";
      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 20 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met20Njet_dt[0][i]<<"            "
			 <<metstuff.Met20Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met20Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met20Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met20Njet_dt[0][i]<<"            "
		      <<metstuff.Met20Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met20Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met20Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met20Njet_dt[1][i]<<"            "
			 <<metstuff.Met20Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met20Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met20Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met20Njet_dt[1][i]<<"            "
		      <<metstuff.Met20Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met20Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met20Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met20Njet_dt[2][i]<<"            "
			 <<metstuff.Met20Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met20Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met20Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met20Njet_dt[2][i]<<"            "
		      <<metstuff.Met20Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met20Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met20Njet_gensyst[2][i]<<std::endl; 
	}

      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 25 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met25Njet_dt[0][i]<<"            "
			 <<metstuff.Met25Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met25Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met25Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met25Njet_dt[0][i]<<"            "
		      <<metstuff.Met25Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met25Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met25Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met25Njet_dt[1][i]<<"            "
			 <<metstuff.Met25Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met25Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met25Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met25Njet_dt[1][i]<<"            "
		      <<metstuff.Met25Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met25Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met25Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met25Njet_dt[2][i]<<"            "
			 <<metstuff.Met25Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met25Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met25Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met25Njet_dt[2][i]<<"            "
		      <<metstuff.Met25Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met25Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met25Njet_gensyst[2][i]<<std::endl; 
	}

      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 30 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met30Njet_dt[0][i]<<"            "
			 <<metstuff.Met30Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met30Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met30Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met30Njet_dt[0][i]<<"            "
		      <<metstuff.Met30Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met30Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met30Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met30Njet_dt[1][i]<<"            "
			 <<metstuff.Met30Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met30Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met30Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met30Njet_dt[1][i]<<"            "
		      <<metstuff.Met30Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met30Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met30Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met30Njet_dt[2][i]<<"            "
			 <<metstuff.Met30Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met30Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met30Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met30Njet_dt[2][i]<<"            "
		      <<metstuff.Met30Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met30Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met30Njet_gensyst[2][i]<<std::endl; 
	}

      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 35 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met35Njet_dt[0][i]<<"            "
			 <<metstuff.Met35Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met35Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met35Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met35Njet_dt[0][i]<<"            "
		      <<metstuff.Met35Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met35Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met35Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met35Njet_dt[1][i]<<"            "
			 <<metstuff.Met35Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met35Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met35Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met35Njet_dt[1][i]<<"            "
		      <<metstuff.Met35Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met35Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met35Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met35Njet_dt[2][i]<<"            "
			 <<metstuff.Met35Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met35Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met35Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met35Njet_dt[2][i]<<"            "
		      <<metstuff.Met35Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met35Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met35Njet_gensyst[2][i]<<std::endl; 
	}

      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 40 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met40Njet_dt[0][i]<<"            "
			 <<metstuff.Met40Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met40Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met40Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met40Njet_dt[0][i]<<"            "
		      <<metstuff.Met40Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met40Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met40Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met40Njet_dt[1][i]<<"            "
			 <<metstuff.Met40Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met40Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met40Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met40Njet_dt[1][i]<<"            "
		      <<metstuff.Met40Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met40Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met40Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met40Njet_dt[2][i]<<"            "
			 <<metstuff.Met40Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met40Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met40Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met40Njet_dt[2][i]<<"            "
		      <<metstuff.Met40Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met40Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met40Njet_gensyst[2][i]<<std::endl; 
	}

      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 45 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met45Njet_dt[0][i]<<"            "
			 <<metstuff.Met45Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met45Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met45Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met45Njet_dt[0][i]<<"            "
		      <<metstuff.Met45Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met45Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met45Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met45Njet_dt[1][i]<<"            "
			 <<metstuff.Met45Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met45Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met45Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met45Njet_dt[1][i]<<"            "
		      <<metstuff.Met45Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met45Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met45Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met45Njet_dt[2][i]<<"            "
			 <<metstuff.Met45Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met45Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met45Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met45Njet_dt[2][i]<<"            "
		      <<metstuff.Met45Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met45Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met45Njet_gensyst[2][i]<<std::endl; 
	}

      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 50 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met50Njet_dt[0][i]<<"            "
			 <<metstuff.Met50Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met50Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met50Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met50Njet_dt[0][i]<<"            "
		      <<metstuff.Met50Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met50Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met50Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met50Njet_dt[1][i]<<"            "
			 <<metstuff.Met50Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met50Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met50Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met50Njet_dt[1][i]<<"            "
		      <<metstuff.Met50Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met50Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met50Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met50Njet_dt[2][i]<<"            "
			 <<metstuff.Met50Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met50Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met50Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met50Njet_dt[2][i]<<"            "
		      <<metstuff.Met50Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met50Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met50Njet_gensyst[2][i]<<std::endl; 
	}

      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 75 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met75Njet_dt[0][i]<<"            "
			 <<metstuff.Met75Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met75Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met75Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met75Njet_dt[0][i]<<"            "
		      <<metstuff.Met75Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met75Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met75Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met75Njet_dt[1][i]<<"            "
			 <<metstuff.Met75Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met75Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met75Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met75Njet_dt[1][i]<<"            "
		      <<metstuff.Met75Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met75Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met75Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met75Njet_dt[2][i]<<"            "
			 <<metstuff.Met75Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met75Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met75Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met75Njet_dt[2][i]<<"            "
		      <<metstuff.Met75Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met75Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met75Njet_gensyst[2][i]<<std::endl; 
	}
      
      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 100 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met100Njet_dt[0][i]<<"            "
			 <<metstuff.Met100Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met100Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met100Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met100Njet_dt[0][i]<<"            "
		      <<metstuff.Met100Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met100Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met100Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met100Njet_dt[1][i]<<"            "
			 <<metstuff.Met100Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met100Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met100Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met100Njet_dt[1][i]<<"            "
		      <<metstuff.Met100Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met100Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met100Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met100Njet_dt[2][i]<<"            "
			 <<metstuff.Met100Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met100Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met100Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met100Njet_dt[2][i]<<"            "
		      <<metstuff.Met100Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met100Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met100Njet_gensyst[2][i]<<std::endl; 
	}
      
      outfile<<"________________________________________________________________"<<"\n";
      outfile<<"                                                                "<<"\n";
      outfile<<" Met > 150 GeV         observed         predicted"<<"\n";
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet15="<<i<<"            "
			 <<metstuff.Met150Njet_dt[0][i]<<"            "
			 <<metstuff.Met150Njet_gen[0][i]<<" +/- "
			 <<metstuff.Met150Njet_genstat[0][i]<<" +/- "
			 <<metstuff.Met150Njet_gensyst[0][i]<<"\n";
	  else outfile<<" Njet15>="<<i<<"            "
		      <<metstuff.Met150Njet_dt[0][i]<<"            "
		      <<metstuff.Met150Njet_gen[0][i]<<" +/- "
		      <<metstuff.Met150Njet_genstat[0][i]<<" +/- "
		      <<metstuff.Met150Njet_gensyst[0][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet20="<<i<<"            "
			 <<metstuff.Met150Njet_dt[1][i]<<"            "
			 <<metstuff.Met150Njet_gen[1][i]<<" +/- "
			 <<metstuff.Met150Njet_genstat[1][i]<<" +/- "
			 <<metstuff.Met150Njet_gensyst[1][i]<<"\n";
	  else outfile<<" Njet20>="<<i<<"            "
		      <<metstuff.Met150Njet_dt[1][i]<<"            "
		      <<metstuff.Met150Njet_gen[1][i]<<" +/- "
		      <<metstuff.Met150Njet_genstat[1][i]<<" +/- "
		      <<metstuff.Met150Njet_gensyst[1][i]<<std::endl; 
	}
      outfile<<"                                                                "<<"\n";
      for(int i=0; i<10; i++)
	{
	  if(i<9) outfile<<"  Njet25="<<i<<"            "
			 <<metstuff.Met150Njet_dt[2][i]<<"            "
			 <<metstuff.Met150Njet_gen[2][i]<<" +/- "
			 <<metstuff.Met150Njet_genstat[2][i]<<" +/- "
			 <<metstuff.Met150Njet_gensyst[2][i]<<"\n";
	  else outfile<<" Njet25>="<<i<<"            "
		      <<metstuff.Met150Njet_dt[2][i]<<"            "
		      <<metstuff.Met150Njet_gen[2][i]<<" +/- "
		      <<metstuff.Met150Njet_genstat[2][i]<<" +/- "
		      <<metstuff.Met150Njet_gensyst[2][i]<<std::endl; 
	}

      outfile<<"---- The End!------"<<std::endl;
    }
  return;
}
/*}}}*/

//_______________ selects large MET events; writes them to file
int TMyJetFilterModule::MyLargeMetEvent(JetStuff jetstuff, CommonStuff miscstuff,
					int metscenario, int Nvx12, int runN, int eventN)
/*{{{*/
{
  int largemet_code=0;
  int Nexo_obj=0;
  //__________________cuts
  double cut_MEt=50.0;
  
  //_________________ contribution from MET
  if(metscenario==0 && miscstuff.myMET_raw.Mod() > cut_MEt) Nexo_obj++;
  if(metscenario==1 && jetstuff.myMETcorr_th5.Mod() > cut_MEt) Nexo_obj++;
  if(metscenario==2 && jetstuff.myMETcorr_th10.Mod() > cut_MEt) Nexo_obj++;
  if(metscenario==3 && jetstuff.myMETcorr_th15.Mod() > cut_MEt) Nexo_obj++;
  if(metscenario==4 && jetstuff.myMETcorr_th20.Mod() > cut_MEt) Nexo_obj++;
  if((metscenario<0 || metscenario>4) && miscstuff.myMET_raw.Mod() > cut_MEt) Nexo_obj++;

  //__________________________ printing out interesting event
  if(Nexo_obj>0)
    {
      largemet_code=1;
      //____________________ opening existing file
      ofstream outfile(fLargeMetEventFileName, std::ios::app);
      if (! outfile) 
	{
	  std::cerr<<"Error in openning of .DAT file\n";
	}
      if(outfile)
	{  
	  char outputstring[200];
	  outfile<<"_____________________________________________________________________________"<<"\n";
	  outfile<<"Run "<<runN<<"  Event "<<eventN<<"\n";
	  outfile<<"Nvx12 "<<Nvx12<<"\n";	  
	  outfile<<"_______ Npho="<<miscstuff.myCorrPhoton.size()<<"\n";
	  //_________________ printing photons
	  for(int i=0; i<miscstuff.myCorrPhoton.size(); i++)
	    {
	      outfile<<"...... pho-"<<i+1
		     <<" Et="<<miscstuff.myCorrPhoton[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrPhoton[i].Phi())
		     <<" eta="<<miscstuff.myCorrPhoton[i].Eta()
		     <<" eta_det="<<miscstuff.myPhoEtaDet[i]<<"\n";
	    }
	  outfile<<"_______ Nele="<<miscstuff.myCorrElectron.size()<<"\n";
	  //_________________ printing electrons
	  for(int i=0; i<miscstuff.myCorrElectron.size(); i++)
	    {
	      outfile<<"...... ele-"<<i+1
		     <<" Et="<<miscstuff.myCorrElectron[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrElectron[i].Phi())
		     <<" eta="<<miscstuff.myCorrElectron[i].Eta()
		     <<" eta_det="<<miscstuff.myEleEtaDet[i]<<"\n";
	    }
	  outfile<<"_______ Njet5="<<jetstuff.myNjet_th5<<"\n";
	  //_________________ printing jets
	  for(int i=0; i<jetstuff.myNjet; i++)
	    {
	      if(jetstuff.Jet_lev6_noEMobj[i].Pt() >=5.0) outfile<<"...... jet-"<<i+1
								  <<" Et="<<jetstuff.Jet_lev6_noEMobj[i].Pt()
								  <<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj[i].Phi())
								  <<" eta="<<jetstuff.Jet_lev6_noEMobj[i].Eta()
 								  <<" eta_det="<<jetstuff.EtaDet[i]
 								  <<" raw_em_fr="<<jetstuff.EmFrRaw[i]<<"\n";

	    }    
	  //_________________ printing MET
	  if(metscenario==0) outfile<<"_______ MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";
	  if(metscenario==1) outfile<<"_______ MET="<<jetstuff.myMETcorr_th5.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th5.Phi())<<"\n";
	  if(metscenario==2) outfile<<"_______ MET="<<jetstuff.myMETcorr_th10.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th10.Phi())<<"\n";
	  if(metscenario==3) outfile<<"_______ MET="<<jetstuff.myMETcorr_th15.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th15.Phi())<<"\n";
	  if(metscenario==4) outfile<<"_______ MET="<<jetstuff.myMETcorr_th20.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th20.Phi())<<"\n";
	  if(metscenario<0 || metscenario>4) outfile<<"_______ MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";
	  outfile<<"............................................................................."<<std::endl;	  
	}
    }
  return largemet_code;
}
/*}}}*/

//_______________ event print-out
void TMyJetFilterModule::MyDumpEvent(JetStuff jetstuff, CommonStuff miscstuff, 
						int metscenario, int Nvx12, int runN, int eventN)
/*{{{*/
{
  
  //__________________________ printing out interesting event
  if(fDumpEvent==1)
    {
      //____________________ opening existing file
      ofstream outfile(fDumpEventFileName, std::ios::app);
      if (! outfile) 
	{
	  std::cerr<<"Error in openning of .DAT file\n";
	}
      if(outfile)
	{  
	  char outputstring[200];
	  outfile<<"_____________________________________________________________________________"<<"\n";
	  outfile<<"Run "<<runN<<"  Event "<<eventN<<"\n";
	  outfile<<"Nvx12 "<<Nvx12<<"\n";	  
	  outfile<<"_______ Npho="<<miscstuff.myCorrPhoton.size()<<"\n";
	  //_________________ printing photons
	  outfile<<"___ printing corrected photons ____"<<"\n"; 
	  for(int i=0; i<miscstuff.myCorrPhoton.size(); i++)
	    {
	      outfile<<"...... pho-"<<i+1
		     <<" Et="<<miscstuff.myCorrPhoton[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrPhoton[i].Phi())
		     <<" eta="<<miscstuff.myCorrPhoton[i].Eta()
		     <<" eta_det="<<miscstuff.myPhoEtaDet[i]<<"\n";
	    }
	  outfile<<"___ printing raw photons ____"<<"\n";
	  for(int i=0; i<miscstuff.myCorrPhoton.size(); i++)
	    {
	      outfile<<"...... pho-"<<i+1
		     <<" raw Et="<<miscstuff.myRawPhoton[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myRawPhoton[i].Phi())
		     <<" eta="<<miscstuff.myRawPhoton[i].Eta()
		     <<" index="<<miscstuff.myPhoInd[i]<<"\n";
	    }
	  
	  outfile<<"_______ Nele="<<miscstuff.myCorrElectron.size()<<"\n";
	  //_________________ printing electrons
	  outfile<<"___ printing corrected electrons ____"<<"\n";
	  for(int i=0; i<miscstuff.myCorrElectron.size(); i++)
	    {
	      outfile<<"...... ele-"<<i+1
		     <<" Et="<<miscstuff.myCorrElectron[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrElectron[i].Phi())
		     <<" eta="<<miscstuff.myCorrElectron[i].Eta()
		     <<" eta_det="<<miscstuff.myEleEtaDet[i]<<"\n";
	    }
	  outfile<<"___ printing raw electrons ____"<<"\n";
	  for(int i=0; i<miscstuff.myCorrElectron.size(); i++)
	    {
	      outfile<<"...... ele-"<<i+1
		     <<" raw Et="<<miscstuff.myRawElectron[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myRawElectron[i].Phi())
		     <<" eta="<<miscstuff.myRawElectron[i].Eta()
		     <<" index="<<miscstuff.myEleInd[i]<<"\n";
	    }
	  
	  outfile<<"_______ Njet="<<jetstuff.myNjet<<"\n";
	  //_________________ printing jets
	  outfile<<"___ printing corrected jets after removing EM object ____"<<"\n";
	  for(int i=0; i<jetstuff.myNjet; i++)
	    {
	      outfile<<"...... jet-"<<i+1
		     <<" Et="<<jetstuff.Jet_lev6_noEMobj[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj[i].Phi())
		     <<" eta="<<jetstuff.Jet_lev6_noEMobj[i].Eta()
		     <<" cor_eta_det="<<jetstuff.EtaDetCorr[i]
		     <<" cor_em_fr="<<jetstuff.EmFrCorr[i]<<"\n";
	    }
	  outfile<<"___ printing raw jets after removing EM object ____"<<"\n";
	  for(int i=0; i<jetstuff.myNjet; i++)
	    {
	      outfile<<"...... jet-"<<i+1
		     <<" Et="<<jetstuff.Jet_raw_noEMobj[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_raw_noEMobj[i].Phi())
		     <<" eta="<<jetstuff.Jet_raw_noEMobj[i].Eta()
		     <<" cor_eta_det="<<jetstuff.EtaDetCorr[i]
		     <<" cor_em_fr="<<jetstuff.EmFrCorr[i]<<"\n";
	    }

	  outfile<<"___ printing corrected jets before removing EM object ____"<<"\n";
	  for(int i=0; i<jetstuff.myNjet; i++)
	    {
	      outfile<<"...... jet-"<<i+1
		     <<" Et="<<jetstuff.Jet_lev6[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_lev6[i].Phi())
		     <<" eta="<<jetstuff.Jet_lev6[i].Eta()
		     <<" raw_eta_det="<<jetstuff.EtaDet[i]
		     <<" raw_em_fr="<<jetstuff.EmFrRaw[i]<<"\n";
	    }
	  outfile<<"___ printing raw jets before removing EM object ____"<<"\n";
	  for(int i=0; i<jetstuff.myNjet; i++)
	    {
	      outfile<<"...... jet-"<<i+1
		     <<" Et="<<jetstuff.Jet_raw[i].Pt()
		     <<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_raw[i].Phi())
		     <<" eta="<<jetstuff.Jet_raw[i].Eta()
		     <<" raw_eta_det="<<jetstuff.EtaDet[i]
		     <<" raw_em_fr="<<jetstuff.EmFrRaw[i]<<"\n";
	    }

    
	  //_________________ printing MET
	  outfile<<"_______ raw MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";
	  outfile<<"_______ corr MET="<<jetstuff.myMETcorr_th15.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th15.Phi())<<"\n";
	  outfile<<"............................................................................."<<std::endl;	  
	}
    }
  return;
}
/*}}}*/


//_______  JetResolution histograms
void TMyJetFilterModule::DoJetResHisto1(JetGeneral_t& Hist) 
/*{{{*/
{ 
  
  Hist.fJetRes_DetEta_j5_FrSig3dPhi04_b->Sumw2();
  Hist.fJetRes_DetEta_j10_FrSig3dPhi04_b->Sumw2();
  Hist.fJetRes_DetEta_j15_FrSig3dPhi04_b->Sumw2();
  Hist.fJetRes_DetEta_j5_FrSig3dPhi01_b->Sumw2();
  Hist.fJetRes_DetEta_j10_FrSig3dPhi01_b->Sumw2();
  Hist.fJetRes_DetEta_j15_FrSig3dPhi01_b->Sumw2();

  Hist.fJetRes_DetEta_j5_FrSig3dPhi04_b->Divide(Hist.fJetRes_DetEta_j5_Sig3dPhi04_b,Hist.fJetRes_DetEta_j5_b,1.0,1.0);
  Hist.fJetRes_DetEta_j10_FrSig3dPhi04_b->Divide(Hist.fJetRes_DetEta_j10_Sig3dPhi04_b,Hist.fJetRes_DetEta_j10_b,1.0,1.0);
  Hist.fJetRes_DetEta_j15_FrSig3dPhi04_b->Divide(Hist.fJetRes_DetEta_j15_Sig3dPhi04_b,Hist.fJetRes_DetEta_j15_b,1.0,1.0);
  Hist.fJetRes_DetEta_j5_FrSig3dPhi01_b->Divide(Hist.fJetRes_DetEta_j5_Sig3dPhi01_b,Hist.fJetRes_DetEta_j5_b,1.0,1.0);
  Hist.fJetRes_DetEta_j10_FrSig3dPhi01_b->Divide(Hist.fJetRes_DetEta_j10_Sig3dPhi01_b,Hist.fJetRes_DetEta_j10_b,1.0,1.0);
  Hist.fJetRes_DetEta_j15_FrSig3dPhi01_b->Divide(Hist.fJetRes_DetEta_j15_Sig3dPhi01_b,Hist.fJetRes_DetEta_j15_b,1.0,1.0);
  
  Hist.fJetRes_DetEta_j5_FrSig3dPhi04_a->Sumw2();
  Hist.fJetRes_DetEta_j10_FrSig3dPhi04_a->Sumw2();
  Hist.fJetRes_DetEta_j15_FrSig3dPhi04_a->Sumw2();
  Hist.fJetRes_DetEta_j5_FrSig3dPhi01_a->Sumw2();
  Hist.fJetRes_DetEta_j10_FrSig3dPhi01_a->Sumw2();
  Hist.fJetRes_DetEta_j15_FrSig3dPhi01_a->Sumw2();

  Hist.fJetRes_DetEta_j5_FrSig3dPhi04_a->Divide(Hist.fJetRes_DetEta_j5_Sig3dPhi04_a,Hist.fJetRes_DetEta_j5_a,1.0,1.0);
  Hist.fJetRes_DetEta_j10_FrSig3dPhi04_a->Divide(Hist.fJetRes_DetEta_j10_Sig3dPhi04_a,Hist.fJetRes_DetEta_j10_a,1.0,1.0);
  Hist.fJetRes_DetEta_j15_FrSig3dPhi04_a->Divide(Hist.fJetRes_DetEta_j15_Sig3dPhi04_a,Hist.fJetRes_DetEta_j15_a,1.0,1.0);
  Hist.fJetRes_DetEta_j5_FrSig3dPhi01_a->Divide(Hist.fJetRes_DetEta_j5_Sig3dPhi01_a,Hist.fJetRes_DetEta_j5_a,1.0,1.0);
  Hist.fJetRes_DetEta_j10_FrSig3dPhi01_a->Divide(Hist.fJetRes_DetEta_j10_Sig3dPhi01_a,Hist.fJetRes_DetEta_j10_a,1.0,1.0);
  Hist.fJetRes_DetEta_j15_FrSig3dPhi01_a->Divide(Hist.fJetRes_DetEta_j15_Sig3dPhi01_a,Hist.fJetRes_DetEta_j15_a,1.0,1.0);

  return;
}
/*}}}*/

//________ creates a list of towers for all jets
void TMyJetFilterModule::MatchCalorTowers(int jet_ind,
					    TStnJetBlock* fJetBlock,
					    TCalDataBlock *fCalDataBlock,
					    CalDataArray* cdholder)
/*{{{*/
{
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

  std::sort(cdholder->begin(), cdholder->end(), ::SortTowersByEnergy);
  
  return;
}
/*}}}*/

//________ creates a list of towers for all EM objects
void TMyJetFilterModule::MatchCalorTowers(int pho_ind,
					    TStnPhotonBlock* fPhotonBlock,
					    TCalDataBlock *fCalDataBlock,
					    CalDataArray* cdholder)
/*{{{*/
{
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
/*}}}*/

//________ creates a list of towers for electrons
void TMyJetFilterModule::MatchCalorTowers(int ele_ind,
					  TStnElectronBlock* fElectronBlock,
					  TCalDataBlock *fCalDataBlock,
					  CalDataArray* cdholder)
/*{{{*/
{
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
/*}}}*/

//_____________________________________________________________________________
int TMyJetFilterModule::BeginJob() {


//--------------------------------------------------------------------------------
//_______________ ***** met model parametrization.
//                for MetX & MetY: sigma=p0+p1*sqrt(SumEt)
//                for MetX & MetY: mean=p0+p1*SumEt  

	if(fJTC_imode==1) //data
	{
      //____________________________ di-pho sideband parametrization for events with Njet15 = 0
      //------------- sideband parametrization after 01/11/07
      metmodel_sigmaX[0][0]=1.02;
      metmodel_sigmaX[0][1]=0.397;
      metmodel_sigmaX_er[0][0]=0.26;
      metmodel_sigmaX_er[0][1]=0.032;
      metmodel_sigmaY[0][0]=0.96;
      metmodel_sigmaY[0][1]=0.394;
      metmodel_sigmaY_er[0][0]=0.25;
      metmodel_sigmaY_er[0][1]=0.03;
      metmodel_meanX[0][0]=0.050;
      metmodel_meanX[0][1]=0.0056;
      metmodel_meanX_er[0][0]=0.066;
      metmodel_meanX_er[0][1]=0.0010;
      metmodel_meanY[0][0]=-0.02;
      metmodel_meanY[0][1]=-0.0041;
      metmodel_meanY_er[0][0]=0.06;
      metmodel_meanY_er[0][1]=0.0009;
      
      metmodel_sigmaXcorr[0]=-0.927;
      metmodel_sigmaYcorr[0]=-0.923;
      metmodel_meanXcorr[0]=-0.850;
      metmodel_meanYcorr[0]=-0.855;
      
      //------------- Z->ee parametrization after 01/11/07
      metmodel_sigmaX[1][0]=1.78;
      metmodel_sigmaX[1][1]=0.356;
      metmodel_sigmaX_er[1][0]=0.23;
      metmodel_sigmaX_er[1][1]=0.026;
      metmodel_sigmaY[1][0]=1.57;
      metmodel_sigmaY[1][1]=0.391;
      metmodel_sigmaY_er[1][0]=0.26;
      metmodel_sigmaY_er[1][1]=0.03;
      metmodel_meanX[1][0]=-0.087;
      metmodel_meanX[1][1]=0.0071;
      metmodel_meanX_er[1][0]=0.068;
      metmodel_meanX_er[1][1]=0.001;
      metmodel_meanY[1][0]=-0.092;
      metmodel_meanY[1][1]=-0.0038;
      metmodel_meanY_er[1][0]=0.064;
      metmodel_meanY_er[1][1]=0.0009;
      
      metmodel_sigmaXcorr[1]=-0.917;
      metmodel_sigmaYcorr[1]=-0.921;
      metmodel_meanXcorr[1]=-0.798;
      metmodel_meanYcorr[1]=-0.839;
   } 

	if(fJTC_imode==0)  // MC
   {
   	if(fUnclParamSwitch==0) // 0=dipho parametrization
		{
		  //____________________________ Pythia di-pho signal parametrization for events with Njet15 = 0
		  //------------- parametrization after 01/31/07
		  metmodel_sigmaX[0][0]=0.45;
		  metmodel_sigmaX[0][1]=0.512;
		  metmodel_sigmaX_er[0][0]=0.42;
		  metmodel_sigmaX_er[0][1]=0.061;
		  metmodel_sigmaY[0][0]=0.67;
		  metmodel_sigmaY[0][1]=0.470;
		  metmodel_sigmaY_er[0][0]=0.36;
		  metmodel_sigmaY_er[0][1]=0.051;
		  metmodel_meanX[0][0]=0.066;
		  metmodel_meanX[0][1]=0.0065;
		  metmodel_meanX_er[0][0]=0.047;
		  metmodel_meanX_er[0][1]=0.0013;
		  metmodel_meanY[0][0]=-0.150;
		  metmodel_meanY[0][1]=-0.0116;
		  metmodel_meanY_er[0][0]=0.058;
		  metmodel_meanY_er[0][1]=0.0015;
		  
		  metmodel_sigmaXcorr[0]=-0.942;
		  metmodel_sigmaYcorr[0]=-0.934;
		  metmodel_meanXcorr[0]=-0.811;
		  metmodel_meanYcorr[0]=-0.762;
		  
		  //------------- Pythia Z->ee (CC) parametrization is the same as Pythia di-pho (for now)
		  metmodel_sigmaX[1][0]=0.45;
		  metmodel_sigmaX[1][1]=0.512;
		  metmodel_sigmaX_er[1][0]=0.42;
		  metmodel_sigmaX_er[1][1]=0.061;
		  metmodel_sigmaY[1][0]=0.67;
		  metmodel_sigmaY[1][1]=0.470;
		  metmodel_sigmaY_er[1][0]=0.36;
		  metmodel_sigmaY_er[1][1]=0.051;
		  metmodel_meanX[1][0]=0.066;
		  metmodel_meanX[1][1]=0.0065;
		  metmodel_meanX_er[1][0]=0.047;
		  metmodel_meanX_er[1][1]=0.0013;
		  metmodel_meanY[1][0]=-0.150;
		  metmodel_meanY[1][1]=-0.0116;
		  metmodel_meanY_er[1][0]=0.058;
		  metmodel_meanY_er[1][1]=0.0015;
		  
		  metmodel_sigmaXcorr[1]=-0.942;
		  metmodel_sigmaYcorr[1]=-0.934;
		  metmodel_meanXcorr[1]=-0.811;
		  metmodel_meanYcorr[1]=-0.762;      
		}
      else // 1= cem-cem Z->ee parametrization
		{
		  //------------- Pythia Z->ee (CC) parametrization (as of 04/25/07)
		  metmodel_sigmaX[1][0]=1.57;
		  metmodel_sigmaX[1][1]=0.401;
		  metmodel_sigmaX_er[1][0]=0.21;
		  metmodel_sigmaX_er[1][1]=0.025;
		  metmodel_sigmaY[1][0]=1.66;
		  metmodel_sigmaY[1][1]=0.411;
		  metmodel_sigmaY_er[1][0]=0.19;
		  metmodel_sigmaY_er[1][1]=0.025;
		  metmodel_meanX[1][0]=0.281;
		  metmodel_meanX[1][1]=-0.002;
		  metmodel_meanX_er[1][0]=0.048;
		  metmodel_meanX_er[1][1]=0.0012;
		  metmodel_meanY[1][0]=-0.150;
		  metmodel_meanY[1][1]=-0.00422;
		  metmodel_meanY_er[1][0]=0.028;
		  metmodel_meanY_er[1][1]=0.00055;
		  
		  metmodel_sigmaXcorr[1]=-0.901;
		  metmodel_sigmaYcorr[1]=-0.879;
		  metmodel_meanXcorr[1]=-0.964;
		  metmodel_meanYcorr[1]=-0.787;
		  
		  //----------------------- di-pho papamertization for studies with Pythia Zee events
		  metmodel_sigmaX[0][0]=1.57;
		  metmodel_sigmaX[0][1]=0.401;
		  metmodel_sigmaX_er[0][0]=0.21;
		  metmodel_sigmaX_er[0][1]=0.025;
		  metmodel_sigmaY[0][0]=1.66;
		  metmodel_sigmaY[0][1]=0.411;
		  metmodel_sigmaY_er[0][0]=0.19;
		  metmodel_sigmaY_er[0][1]=0.025;
		  metmodel_meanX[0][0]=0.281;
		  metmodel_meanX[0][1]=-0.002;
		  metmodel_meanX_er[0][0]=0.048;
		  metmodel_meanX_er[0][1]=0.0012;
		  metmodel_meanY[0][0]=-0.150;
		  metmodel_meanY[0][1]=-0.00422;
		  metmodel_meanY_er[0][0]=0.028;
		  metmodel_meanY_er[0][1]=0.00055;
		  
		  metmodel_sigmaXcorr[0]=-0.901;
		  metmodel_sigmaYcorr[0]=-0.879;
		  metmodel_meanXcorr[0]=-0.964;
		  metmodel_meanYcorr[0]=-0.787;
		}
	}

//_____________________________________________________ register the data block
	RegisterDataBlock("ZVertexBlock" ,"TStnVertexBlock"  ,&fZVertexBlock); // ZVertex Collection 
	RegisterDataBlock("MetBlock","TStnMetBlock",&fMetBlock);
	RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlockClu04);

//---------- need this for jet-EM matching
	RegisterDataBlock("CalDataBlock","TCalDataBlock",&fCalData);
	RegisterDataBlock("PhotonBlock","TStnPhotonBlock",&fPhotonBlock);
	RegisterDataBlock("ElectronBlock","TStnElectronBlock",&fElectronBlock);

//_____________________________________________________ book histograms

	BookHistograms();
	CleanMetResults(met_results); // init all arrays to zeros

	if(fUseMetPDF>0) // MetPDF is to be used and Npoints is over-written
	{
		switch (fUseMetPDF)
		{
			case 1: 
			  Npoints=370;
			  break;
			case 2:
			  Npoints=1000;
			  break;
			case 3:
			  Npoints=10000;
			  break;
			case 4:
			  Npoints=15873;
			  break;
			case 5:
			  Npoints=100;
			  break;
			default:
			  Npoints=100;
		}
	}
	else fSelectSigMetEvent=0; // to select events fUseMetPDF must be >0 

  ofstream outfile1(fLargeMetEventFileName, std::ios::out); // re-creating file
  if (! outfile1) 
    {
      std::cerr<<"Error in openning of .DAT file\n";
    }
  ofstream outfile0(fDumpEventFileName, std::ios::out); // re-creating file
  if (! outfile0) 
    {
      std::cerr<<"Error in openning of .DAT file\n";
    }
  if(outfile1)
    {
      outfile1<<"--------------- My large MET event list ------------------------"<<"\n";
      outfile1<<"________________________________________________________________"<<std::endl;
    }
 
  return 0;
}


//_____________________________________________________________________________
int TMyJetFilterModule::BeginRun() {
  return 0;
}

//_____________________________________________________________________________
int TMyJetFilterModule::Event(int ientry) {

	cout << "\t -----  IN TMyJetFilterModule EVENT LOOP ------\n";

	zvx_best = 0.0;
	myNvx_class12 = 0;

	ClearModuleOutput();
	SigMetEvent_status = 0; // re-setting status code for Significan Met event
	bad_EMjet_match_flag = 1; // by default no events flagged "bad" in the beginning
	bool pass   = true;

	myJetRun = GetHeaderBlock()->RunNumber(); // obtaining run number, myJetRun is globaly defined 
	int myJetEvent = GetHeaderBlock()->EventNumber(); // obtaining event number
	double fLum = (GetHeaderBlock()->InstLum())*1.0E-30;

	fZVertexBlock->GetEntry(ientry); // Accessing Veterx Block
	fMetBlock->GetEntry(ientry); //  Accessing MET Block 

	//_________________________________________________________________________
	//---------  Accessing Jet Blocks 
	fJetBlockClu04->GetEntry(ientry);

	if (!fPhotonBlock) {
		pass=false;
		SetPassed(0);
		return -1;
	}
  fPhotonBlock->GetEntry(ientry);

  if (!fCalData) {
      pass=false;
      SetPassed(0);
      return -1;
	}
	fCalData->GetEntry(ientry);

	if (!fElectronBlock) {
		pass=false;
		SetPassed(0);
		return -1;      
	}
	fElectronBlock->GetEntry(ientry);

	DoCommonStuff(allstuff);  // _____________________________ reads raw Met, pho, ele info and fills CommonStuff
	DoMyJet(fJetBlockClu04,allstuff,jet04stuff,fMatchPhoJet04); // does my JetClu04 jets
	DoMyMet(allstuff,jet04stuff); // _________________________ corrects Met & SumEt for JetClu04 jets

	//-------------------------------------------------------------------------------
	// Filling histograms before event cuts >>>>>>>>
	//_______________________________________________________________________________
	FillJetHistogramsB(fHistJet04, jet04stuff, zvx_best - fJetBlockClu04->ZVertex(), myNvx_class12); // filling general histo for JetClu04

  //-------------------------------------------------------------------------------
  // Applying event Cuts
  //_______________________________________________________________________________

  //....... temporary cuts
  pass=false;

  //if(allstuff.myCorrPhoton.size()>=2)		// for sasha's di pho ana
  /*	sam: i already check for one good photon in Init mod
  	if (allstuff.myCorrPhoton.size()>=1) {	//for me just one
      //if(allstuff.myCorrPhoton[0].Pt()>=13.0 && allstuff.myCorrPhoton[1].Pt()>=13.0) pass=true;
      if(allstuff.myCorrPhoton[0].Pt()>=13.0 && allstuff.myCorrPhoton[1].Pt()>=13.0) {
			pass=true;
		} else {
	  		SetPassed(0);
			return 0;
		}
 	}
*/
/*	sam; does not do anything
	if(MyAnalysisCut(jet04stuff,allstuff,fMyMetScenario)==1) pass=true;
  else 
    {
      SetPassed(0);
      return 0;      
    }
*/

	//no harm. just limits the number of min and max jets per event
	//min=0 max=100
	
	if (MyNjetCut(jet04stuff)==1) pass=true; // added on 11/07/05
  	else {
      SetPassed(0);
      return 0;      
	} 

	if (MyAngularCut(jet04stuff,allstuff)==1) pass=true; // added on 03/08/06 (temporary)
	else{
      SetPassed(0);
      return 0;      
	}

  int Nduplicate=MyDuplicateCut(jet04stuff,allstuff);
  if(fRemoveDuplicate==1 && Nduplicate>0) // removing duplicates
    {
      SetPassed(0);
      return 0;      
    }
  else pass=true;
  if((fRemoveDuplicate+Nduplicate)<0) // selecting duplicates
    {
      SetPassed(0);
      return 0;      
    }
  else pass=true;

  int Nbadmet=MyMetCleanUpCut(jet04stuff,allstuff,jet04stuff.Jet_lev6_noEMobj,jet04stuff.myMETcorr_th15); // new as of 04/04/07
  if(fRemoveBadMet==1 && Nbadmet>0) // removing events with bad MET
    {
      SetPassed(0);
      return 0;      
    }
  else pass=true;
  if((fRemoveBadMet+Nbadmet)<0) // selecting events with bad MET
    {
      SetPassed(0);
      return 0;      
    }
  else pass=true;
    
  if(pass)
    {
      SetPassed(1);
      //-------------------------------------------------------------------------------
      // Filling histos for events which pass the cuts
      //_______________________________________________________________________________
      FillMetCleanupHistograms(fCleanup,jet04stuff,allstuff,3); // histo for met cleanup studies
      FillJetHistogramsA(fHistJet04,jet04stuff,zvx_best-fJetBlockClu04->ZVertex(),myNvx_class12); // filling general histo for JetClu04
      FillMetStudyHistograms(fMetStudyJet04sc3,jet04stuff,allstuff,3,myNvx_class12); // met study histo for JetClu 0.4 and Et>15 GeV
      MyMetEventCount(jet04stuff,allstuff,3,met_results); // counting number of data events with Met>cut
      FillPhoMetStudyHistograms(fPhoMetStudy,jet04stuff.myMETcorr_th15,allstuff,zvx_best);
      int met_event=MyLargeMetEvent(jet04stuff,allstuff,3,myNvx_class12,GetHeaderBlock()->RunNumber(),GetHeaderBlock()->EventNumber());
      MyDumpEvent(jet04stuff,allstuff,3,myNvx_class12,GetHeaderBlock()->RunNumber(),GetHeaderBlock()->EventNumber());
      if(fSelectMetEvent==1 && met_event==0) SetPassed(0); 
      
      //_____ selecting events according to MetSignificance
      if(fSelectSigMetEvent==1 && SigMetEvent_status==0) SetPassed(0); // selecting event with Significant Met
      if(fSelectSigMetEvent==-1 && SigMetEvent_status==1) SetPassed(0); // rejecting event with Significant Met
    }   
  else SetPassed(0);
  return 0;
}


//_____________________________________________________________________________
void TMyJetFilterModule::Display() {

  return;
}

int TMyJetFilterModule::EndJob() {

  DoJetResHisto1(fHistJet04);
  MyMetResults(met_results);
  DoFinalPhoMetHisto(fPhoMetStudy);
  return 0;
}

double TMyJetFilterModule::GetCorrection(CorInit settings, TLorentzVector vec, 
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
  double scaleFactor = myJetEnergyCorrections.doEnergyCorrections(P4Jet,emf,etad); 
  return scaleFactor; 
}

double TMyJetFilterModule::GetCorrection(TLorentzVector vec, float& emf, float etad) {

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
  double scaleFactor = myJetEnergyCorrections.doEnergyCorrections(P4Jet,emf,etad); 
  return scaleFactor; 
}

//_____________________________________________________________ Balance of two leading jets
// !!!!!!!!      OLD stuff, but I still need it  !!!!!!!!! 
double TMyJetFilterModule::GetBalance2(double pt1, double pt2, double phi1, double phi2) {
  double balX;
  double balY;
  double bal=9999.0;
  balX=pt1*TMath::Cos(phi1)+pt2*TMath::Cos(phi2);
  balY=pt1*TMath::Sin(phi1)+pt2*TMath::Sin(phi2);
  if(fabs(pt1+pt2)>0.0) bal=TMath::Sqrt(balX*balX+balY*balY)/(pt1+pt2);
  return bal;
}
double TMyJetFilterModule::GetBalance2(std::vector<TLorentzVector> vec){
  double bal=9999.0;
  if(vec.size()>=2)
    {
      double pt_1=vec[0].Pt();
      double pt_2=vec[1].Pt();
      double Phi_1=vec[0].Phi();
      double Phi_2=vec[1].Phi();
      bal=GetBalance2(pt_1,pt_2,Phi_1,Phi_2);
    }
  return bal;
}

//________________________________________________ Theta* or dijet scattering angle
//------------------- Input params are eta's or true rapidities
double TMyJetFilterModule::GetThetaStar(double y1, double y2) {
  double tcos;
  double y=(y1-y2)/2.0;
  double t=0.0;
  tcos=(TMath::Exp(y)-TMath::Exp(-y))/(TMath::Exp(y)+TMath::Exp(-y));
  if(fabs(tcos)<=1.0) t=TMath::ACos(tcos);
  return t;
}
double TMyJetFilterModule::GetThetaStar(std::vector<TLorentzVector> vec) {
  double t=0.0;
  if(vec.size()>=2)
    {
      double Y1=vec[0].Rapidity();
      double Y2=vec[1].Rapidity();
      t=GetThetaStar(Y1, Y2);
    }
  return t;
}


//__________________________________________________ returns Kt of two leading jets
double TMyJetFilterModule::GetMyKtKick(std::vector<TLorentzVector> vec) {
  double kt_kick=0.0;
  if(vec.size()>=2)
    {
      double pt_1=vec[0].Pt();
      double pt_2=vec[1].Pt();
      kt_kick=(pt_1+pt_2)*GetBalance2(vec);
    }  
  return kt_kick;
}
//__________________________________________________ returns Phi of two leading jets
double TMyJetFilterModule::GetMyPhiKt2L(std::vector<TLorentzVector> vec) {
  double phi_kick=-1.0E6;
  if(vec.size()>=2)
    {
      double pt=(vec[0]+vec[1]).Pt();
      if(pt>0.0) phi_kick=TVector2::Phi_0_2pi((vec[0]+vec[1]).Phi());
    }  
  return phi_kick;
}
//__________________________________________________ returns Phi of sumKt of all jets
double TMyJetFilterModule::GetMyPhiKtAll(std::vector<TLorentzVector> vec) {
  double phi_extra=-1.0E6;
  if(vec.size()>0)
    {
      TLorentzVector extra_jets(0.0,0.0,0.0,0.0);
      for(int i=0; i<vec.size(); i++)
	{ 
	  extra_jets=extra_jets+vec[i]; 
	}
      double pt=extra_jets.Pt();
      if(pt>0.0) phi_extra=TVector2::Phi_0_2pi(extra_jets.Phi());
    }  
  return phi_extra;
}

//__________________________________________________ returns 1 if jet is matched
//                                                   calculates dR, dPhi, dEta between pho & jet   
int TMyJetFilterModule::MyMatchPhoJet(int jet_ind, int pho_ind, int ele_ind,TStnJetBlock* fJetBlock,
				      TLorentzVector *mypho,TLorentzVector *myjet, double algo_cone, 
				      double &dR_fj, double &dPhi_fj, double &dEta_fj)
/*{{{*/
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
//   if(fabs(dR)<algo_cone) match_code=1;
  dR_fj=dR;
  dPhi_fj=mypho->DeltaPhi(vec);
  dEta_fj=mypho->Eta()-myjet->Eta();
  return match_code;
}
/*}}}*/

void TMyJetFilterModule::ClearModuleOutput()       // clears the Module output from 		 
{                              // previous event.				        	 

  //-------- my new global output parameters 
  ClearJetStuff(jet04stuff); // clears JetClu-0.4 stuff 
//   ClearJetStuff(jet07stuff); // clears JetClu-0.7 stuff
//   ClearJetStuff(jet10stuff); // clears JetClu-1.0 stuff
  ClearCommonStuff(allstuff); // clears common stuff
  //----- end of my new output params 
  //.............................................................
  //  ClearMatchStuff(matchstuff);

  return;
}

void TMyJetFilterModule::ClearJetStuff(JetStuff &stuff)
{ // clears the JetStuff structures
/*{{{*/  
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
  return; 		 
}
/*}}}*/

//_________________________ clears parameters which do not depend on Jet Cone size
void TMyJetFilterModule::ClearCommonStuff(CommonStuff &stuff)
// clears CommonStuff structures 
/*{{{*/
{  
  stuff.myRawPhoton.clear();    
  stuff.myCorrPhoton.clear();
  stuff.myPhoInd.clear();
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
  stuff.myRawMuon.clear();      
  stuff.myCorrMuon.clear();     
  stuff.myRawTau.clear();       
  stuff.myCorrTau.clear();      
  stuff.myRawBjet.clear();      
  stuff.myCorrBjet.clear();     
  stuff.myMET_raw.Set(-1.0E6,-1.0E6); 
  stuff.myMET0_raw.Set(-1.0E6,-1.0E6); 
  stuff.mySumEt_raw=-1.0E6;

  return;
}
/*}}}*/

void TMyJetFilterModule::ClearMatchStuff() {
  matchstuff.clear();
  return;
}
  

//__________________________________________________________ accessors to the
//__________________________________________________________ output params
//_________ need to specify which jets to return (variable "cone"): 0--cone 0.4; 1--cone 0.7; 2--cone 1.0

TLorentzVector* TMyJetFilterModule::GetMyJet_raw(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_raw.size()) return &_stuff.Jet_raw[i];
  else return NULL;
}
TLorentzVector* TMyJetFilterModule::GetMyJet_lev1(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_lev1.size()) return &_stuff.Jet_lev1[i];
  else return NULL;
} 
TLorentzVector* TMyJetFilterModule::GetMyJet_lev4(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_lev4.size()) return &_stuff.Jet_lev4[i];
  else return NULL;
} 
TLorentzVector* TMyJetFilterModule::GetMyJet_lev5(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_lev5.size()) return &_stuff.Jet_lev5[i];
  else return NULL;
} 
TLorentzVector* TMyJetFilterModule::GetMyJet_lev6(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_lev6.size()) return &_stuff.Jet_lev6[i];
  else return NULL;
} 
TLorentzVector* TMyJetFilterModule::GetMyJetNoEMobj_raw(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_raw_noEMobj.size()) return &_stuff.Jet_raw_noEMobj[i];
  else return NULL;
} 
TLorentzVector* TMyJetFilterModule::GetMyJetNoEMobj_lev1(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_lev1_noEMobj.size()) return &_stuff.Jet_lev1_noEMobj[i];
  else return NULL;
} 
TLorentzVector* TMyJetFilterModule::GetMyJetNoEMobj_lev4(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_lev4_noEMobj.size()) return &_stuff.Jet_lev4_noEMobj[i];
  else return NULL;
}  
TLorentzVector* TMyJetFilterModule::GetMyJetNoEMobj_lev5(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_lev5_noEMobj.size()) return &_stuff.Jet_lev5_noEMobj[i];
  else return NULL;
}  
TLorentzVector* TMyJetFilterModule::GetMyJetNoEMobj_lev6(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_lev6_noEMobj.size()) return &_stuff.Jet_lev6_noEMobj[i];
  else return NULL;
} 
TLorentzVector* TMyJetFilterModule::GetMyJetNoEMobj_lev7(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return NULL;
  if(i>=0 && i<_stuff.Jet_lev7_noEMobj.size()) return &_stuff.Jet_lev7_noEMobj[i];
  else return NULL;
} 
double TMyJetFilterModule::GetMyJetEtaDet(int cone, int i) { // before removing EM object 
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  if(i>=0 && i<_stuff.EtaDet.size()) return _stuff.EtaDet[i];
  else return -1.0E6;
}
double TMyJetFilterModule::GetMyJetEtaDetCorr(int cone, int i) { // after removing EM object
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  if(i>=0 && i<_stuff.EtaDetCorr.size()) return _stuff.EtaDetCorr[i];
  else return -1.0E6;
}
double TMyJetFilterModule::GetMyJetEmFrRaw(int cone, int i) {  // before removing EM object
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  if(i>=0 && i<_stuff.EmFrRaw.size()) return _stuff.EmFrRaw[i];
  else return -1.0E6;
}
double TMyJetFilterModule::GetMyJetEmFrCorr(int cone, int i) { // before removing EM object
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  if(i>=0 && i<_stuff.EmFrCorr.size()) return _stuff.EmFrCorr[i];
  else return -1.0E6;
}
int TMyJetFilterModule::GetMyJetNtrk(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(i>=0 && i<_stuff.JetNtrk.size()) return _stuff.JetNtrk[i];
  else return -1;
}
int TMyJetFilterModule::GetMyJetNtwr(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(i>=0 && i<_stuff.JetNtwr.size()) return _stuff.JetNtwr[i];
  else return -1;
}
int TMyJetFilterModule::GetMyJetNobjMatch(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(i>=0 && i<_stuff.Nobj_match.size()) return _stuff.Nobj_match[i];
  else return -1;
}
int TMyJetFilterModule::GetMyJetNphoMatch(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(i>=0 && i<_stuff.Npho_match.size()) return _stuff.Npho_match[i];
  else return -1;
}
int TMyJetFilterModule::GetMyJetNeleMatch(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(i>=0 && i<_stuff.Nele_match.size()) return _stuff.Nele_match[i];
  else return -1;
}
int TMyJetFilterModule::GetMyJetNmuMatch(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(i>=0 && i<_stuff.Nmu_match.size()) return _stuff.Nmu_match[i];
  else return -1;
}
int TMyJetFilterModule::GetMyJetNtauMatch(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(i>=0 && i<_stuff.Ntau_match.size()) return _stuff.Ntau_match[i];
  else return -1;
}
int TMyJetFilterModule::GetMyJetNbtagMatch(int cone, int i) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(i>=0 && i<_stuff.Nbtag_match.size()) return _stuff.Nbtag_match[i];
  else return -1;
}
int TMyJetFilterModule::GetMyJetBlockInd(int cone, int i) { // original jet index in JetBlock, need this after jet reordering
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(i>=0 && i<_stuff.JetBlockInd.size()) return _stuff.JetBlockInd[i];
  else return -1;
}
int TMyJetFilterModule::GetMyNjet(int cone, int threshold) {// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25  
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  if(cone<0 || cone>2) return -1;
  if(threshold==0) return _stuff.myNjet;
  if(threshold==1) return _stuff.myNjet_th5;
  if(threshold==2) return _stuff.myNjet_th10;
  if(threshold==3) return _stuff.myNjet_th15;
  if(threshold==4) return _stuff.myNjet_th20;
  if(threshold==5) return _stuff.myNjet_th25;
  if(threshold<0 || threshold>5) return -1;
}
double TMyJetFilterModule::GetMySumEtCorr(int cone, int threshold) {// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20  
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  if(threshold==0) return allstuff.mySumEt_raw; // for consistency
  if(threshold==1) return _stuff.mySumEtCorr_th5;
  if(threshold==2) return _stuff.mySumEtCorr_th10;
  if(threshold==3) return _stuff.mySumEtCorr_th15;
  if(threshold==4) return _stuff.mySumEtCorr_th20;
  if(threshold<0 || threshold>4) return -1.0E6;
}
double TMyJetFilterModule::GetMyHtCorr(int cone, int threshold) {// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20 
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  return MyHtAll(_stuff,allstuff,threshold); 
}
double TMyJetFilterModule::GetMyMetCorr(int cone, int threshold) {// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20 
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  if(threshold==0) return allstuff.myMET_raw.Mod(); // for consistency
  if(threshold==1) return _stuff.myMETcorr_th5.Mod();
  if(threshold==2) return _stuff.myMETcorr_th10.Mod();
  if(threshold==3) return _stuff.myMETcorr_th15.Mod();
  if(threshold==4) return _stuff.myMETcorr_th20.Mod();
  if(threshold<0 || threshold>4) return -1.0E6;
}
double TMyJetFilterModule::GetMyMetXCorr(int cone, int threshold) {// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  if(threshold==0) return allstuff.myMET_raw.Px(); // for consistency
  if(threshold==1) return _stuff.myMETcorr_th5.Px();
  if(threshold==2) return _stuff.myMETcorr_th10.Px();
  if(threshold==3) return _stuff.myMETcorr_th15.Px();
  if(threshold==4) return _stuff.myMETcorr_th20.Px();
  if(threshold<0 || threshold>4) return -1.0E6;
}
double TMyJetFilterModule::GetMyMetYCorr(int cone, int threshold) {// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  if(threshold==0) return allstuff.myMET_raw.Py(); // for consistency
  if(threshold==1) return _stuff.myMETcorr_th5.Py();
  if(threshold==2) return _stuff.myMETcorr_th10.Py();
  if(threshold==3) return _stuff.myMETcorr_th15.Py();
  if(threshold==4) return _stuff.myMETcorr_th20.Py();
  if(threshold<0 || threshold>4) return -1.0E6;
}
double TMyJetFilterModule::GetMyMetPhiCorr(int cone, int threshold) {// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return -1.0E6;
  if(threshold==0) return allstuff.myMET_raw.Phi(); // for consistency
  if(threshold==1) return _stuff.myMETcorr_th5.Phi();
  if(threshold==2) return _stuff.myMETcorr_th10.Phi();
  if(threshold==3) return _stuff.myMETcorr_th15.Phi();
  if(threshold==4) return _stuff.myMETcorr_th20.Phi();
  if(threshold<0 || threshold>4) return -1.0E6;
}


//__________________________________________________________ setters of the
//__________________________________________________________ output params

//_________ need to specify which jets to return (variable "cone"): 0--cone 0.4; 1--cone 0.7; 2--cone 1.0
void TMyJetFilterModule::SetMyJet_raw(int cone, int i, TLorentzVector vec) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_raw.size()) _stuff.Jet_raw[i]=vec;
  else _stuff.Jet_raw.push_back(vec);
  return;
}
void TMyJetFilterModule::SetMyJet_lev1(int cone, int i, TLorentzVector vec) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_lev1.size()) _stuff.Jet_lev1[i]=vec;
  else _stuff.Jet_lev1.push_back(vec);
  return;
}
void TMyJetFilterModule::SetMyJet_lev4(int cone, int i, TLorentzVector vec) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_lev4.size()) _stuff.Jet_lev4[i]=vec;
  else _stuff.Jet_lev4.push_back(vec);
  return;
} 
void TMyJetFilterModule::SetMyJet_lev5(int cone, int i, TLorentzVector vec) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_lev5.size()) _stuff.Jet_lev5[i]=vec;
  else _stuff.Jet_lev5.push_back(vec);
  return;
} 
void TMyJetFilterModule::SetMyJet_lev6(int cone, int i, TLorentzVector vec) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_lev6.size()) _stuff.Jet_lev6[i]=vec;
  else _stuff.Jet_lev6.push_back(vec);
  return;
} 
void TMyJetFilterModule::SetMyJetNoEMobj_raw(int cone, int i, TLorentzVector vec)  {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_raw_noEMobj.size()) _stuff.Jet_raw_noEMobj[i]=vec;
  else _stuff.Jet_raw_noEMobj.push_back(vec);
  return;
} 
void TMyJetFilterModule::SetMyJetNoEMobj_lev1(int cone, int i, TLorentzVector vec)  {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_lev1_noEMobj.size()) _stuff.Jet_lev1_noEMobj[i]=vec;
  else _stuff.Jet_lev1_noEMobj.push_back(vec);
  return;
}
void TMyJetFilterModule::SetMyJetNoEMobj_lev4(int cone, int i, TLorentzVector vec) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_lev4_noEMobj.size()) _stuff.Jet_lev4_noEMobj[i]=vec;
  else _stuff.Jet_lev4_noEMobj.push_back(vec);
  return;
}  
void TMyJetFilterModule::SetMyJetNoEMobj_lev5(int cone, int i, TLorentzVector vec) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_lev5_noEMobj.size()) _stuff.Jet_lev5_noEMobj[i]=vec;
  else _stuff.Jet_lev5_noEMobj.push_back(vec);
  return;
}  
void TMyJetFilterModule::SetMyJetNoEMobj_lev6(int cone, int i, TLorentzVector vec) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Jet_lev6_noEMobj.size()) _stuff.Jet_lev6_noEMobj[i]=vec;
  else _stuff.Jet_lev6_noEMobj.push_back(vec);
  return;
} 
void TMyJetFilterModule::SetMyJetEtaDet(int cone, int i, double param) { // before removing EM object
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.EtaDet.size()) _stuff.EtaDet[i]=param;
  else _stuff.EtaDet.push_back(param);
  return;
} 
void TMyJetFilterModule::SetMyJetEtaDetCorr(int cone, int i, double param) { // after removing EM object
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.EtaDetCorr.size()) _stuff.EtaDetCorr[i]=param;
  else _stuff.EtaDetCorr.push_back(param);
  return;
} 
void TMyJetFilterModule::SetMyJetEmFrRaw(int cone, int i, double param) { // before removing EM object
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.EmFrRaw.size()) _stuff.EmFrRaw[i]=param;
  else _stuff.EmFrRaw.push_back(param);
  return;
} 
void TMyJetFilterModule::SetMyJetEmFrCorr(int cone, int i, double param) { // before removing EM object
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.EmFrCorr.size()) _stuff.EmFrCorr[i]=param;
  else _stuff.EmFrCorr.push_back(param);
  return;
}
void TMyJetFilterModule::SetMyJetNtrk(int cone, int i, int param) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.JetNtrk.size()) _stuff.JetNtrk[i]=param;
  else _stuff.JetNtrk.push_back(param);
  return;
}
void TMyJetFilterModule::SetMyJetNtwr(int cone, int i, int param) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.JetNtwr.size()) _stuff.JetNtwr[i]=param;
  else _stuff.JetNtwr.push_back(param);
  return;
}
void TMyJetFilterModule::SetMyJetNobjMatch(int cone, int i, int param) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Nobj_match.size()) _stuff.Nobj_match[i]=param;
  else _stuff.Nobj_match.push_back(param);
  return;
 }
void TMyJetFilterModule::SetMyJetNphoMatch(int cone, int i, int param) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Npho_match.size()) _stuff.Npho_match[i]=param;
  else _stuff.Npho_match.push_back(param);
  return;
}
void TMyJetFilterModule::SetMyJetNeleMatch(int cone, int i, int param) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Nele_match.size()) _stuff.Nele_match[i]=param;
  else _stuff.Nele_match.push_back(param);
  return;
}
void TMyJetFilterModule::SetMyJetNmuMatch(int cone, int i, int param) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Nmu_match.size()) _stuff.Nmu_match[i]=param;
  else _stuff.Nmu_match.push_back(param);
  return;
}
void TMyJetFilterModule::SetMyJetNtauMatch(int cone, int i, int param) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Ntau_match.size()) _stuff.Ntau_match[i]=param;
  else _stuff.Ntau_match.push_back(param);
  return;
}
void TMyJetFilterModule::SetMyJetNbtagMatch(int cone, int i, int param) {
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.Nbtag_match.size()) _stuff.Nbtag_match[i]=param;
  else _stuff.Nbtag_match.push_back(param);
  return;
}
void TMyJetFilterModule::SetMyJetBlockInd(int cone, int i, int param) { 
// original jet index in JetBlock, need this after jet reordering
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(i>=0 && i<_stuff.JetBlockInd.size()) _stuff.JetBlockInd[i]=param;
  else _stuff.JetBlockInd.push_back(param);
  return;
}
void TMyJetFilterModule::SetMyNjet(int cone, int threshold, int param) {
  // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(threshold==0) _stuff.myNjet=param;
  if(threshold==1) _stuff.myNjet_th5=param;
  if(threshold==2) _stuff.myNjet_th10=param;
  if(threshold==3) _stuff.myNjet_th15=param;
  if(threshold==4) _stuff.myNjet_th20=param;
  if(threshold==5) _stuff.myNjet_th25=param;
  if(threshold<0 || threshold>5) return;
  return;
}
void TMyJetFilterModule::SetMySumEtCorr(int cone, int threshold, double param) {
// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(threshold==0) allstuff.mySumEt_raw=param; // for consistency
  if(threshold==1) _stuff.mySumEtCorr_th5=param;
  if(threshold==2) _stuff.mySumEtCorr_th10=param;
  if(threshold==3) _stuff.mySumEtCorr_th15=param;
  if(threshold==4) _stuff.mySumEtCorr_th20=param;
  if(threshold<0 || threshold>4) return;
  return;
}
void TMyJetFilterModule::SetMyMetCorr(int cone, int threshold, double px, double py) { 
// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(threshold==0) allstuff.myMET_raw.Set(px,py); // for consistency
  if(threshold==1) _stuff.myMETcorr_th5.Set(px,py);
  if(threshold==2) _stuff.myMETcorr_th10.Set(px,py);
  if(threshold==3) _stuff.myMETcorr_th15.Set(px,py);
  if(threshold==4) _stuff.myMETcorr_th20.Set(px,py);
  if(threshold<0 || threshold>4) return;
  return;
}
void TMyJetFilterModule::SetMyMetCorr(int cone, int threshold, const TVector2& vec) {
 // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  JetStuff _stuff;
  if(cone==0) _stuff=jet04stuff;
  else return;
  if(threshold==0) allstuff.myMET_raw.Set(vec); // for consistency
  if(threshold==1) _stuff.myMETcorr_th5.Set(vec);
  if(threshold==2) _stuff.myMETcorr_th10.Set(vec);
  if(threshold==3) _stuff.myMETcorr_th15.Set(vec);
  if(threshold==4) _stuff.myMETcorr_th20.Set(vec);
  if(threshold<0 || threshold>4) return;
  return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//   My Cuts are defined here
//___________________________________________________________________
int TMyJetFilterModule::MyNjetCut(JetStuff jetstuff) {
  if(jetstuff.myNjet<MinNjet || jetstuff.myNjet>MaxNjet) return 0;
  if(jetstuff.myNjet_th5<MinNjet5 || jetstuff.myNjet_th5>MaxNjet5) return 0;
  if(jetstuff.myNjet_th10<MinNjet10 || jetstuff.myNjet_th10>MaxNjet10) return 0;
  if(jetstuff.myNjet_th15<MinNjet15 || jetstuff.myNjet_th15>MaxNjet15) return 0;
  if(jetstuff.myNjet_th20<MinNjet20 || jetstuff.myNjet_th20>MaxNjet20) return 0;
  if(jetstuff.myNjet_th25<MinNjet25 || jetstuff.myNjet_th25>MaxNjet25) return 0;
  return 1;
}

//--------------------------------------------------------------------------
int TMyJetFilterModule::MyGenericCut(JetStuff jetstuff, double dz)
/*{{{*/
{
  //__________________ cut on dZ=Zvx-Zjet
  if(fabs(dz)<MindZ || fabs(dz)>MaxdZ) return 0;
  //__________________ cut on 1st jet Et
  if(jetstuff.myNjet>0) 
    if(fabs(jetstuff.Jet_lev6_noEMobj[0].Eta())>MinJetEta && 
       fabs(jetstuff.Jet_lev6_noEMobj[0].Eta())<MaxJetEta &&
       (jetstuff.Jet_lev6_noEMobj[0].Pt()<MinEt1st || jetstuff.Jet_lev6_noEMobj[0].Pt()>MaxEt1st)) return 0;
  //__________________ cut on 2nd jet Et
  if(jetstuff.myNjet>1) 
    if(fabs(jetstuff.Jet_lev6_noEMobj[1].Eta())>MinJetEta && 
       fabs(jetstuff.Jet_lev6_noEMobj[1].Eta())<MaxJetEta &&
       (jetstuff.Jet_lev6_noEMobj[1].Pt()<MinEt2nd || jetstuff.Jet_lev6_noEMobj[1].Pt()>MaxEt2nd)) return 0;
  return 1;
}
/*}}}*/

//--------------------------------------------------------------------------
int TMyJetFilterModule::MyAngularCut(JetStuff jetstuff, CommonStuff miscstuff)
/*{{{*/
{
	double dPhiMin=100.0;
	double dEtaMin=100.0;
	double dRMin=100.0;

	for (int i=0; i<jetstuff.myNjet; i++) {
      if (fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())>MinJetEta && 
	 		fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<MaxJetEta &&
	 		fabs(jetstuff.Jet_lev6_noEMobj[i].Pt())>MinEtThr && 
	 		fabs(jetstuff.Jet_lev6_noEMobj[i].Pt())<MaxEtThr) {
 	  		// looping over all photons
 	  		for (int j=0; j<miscstuff.myCorrPhoton.size(); j++) {
				double dPhi=fabs(TVector2::Phi_mpi_pi(jetstuff.Jet_lev6_noEMobj[i].Phi()-miscstuff.myCorrPhoton[j].Phi()));
				double dEta=fabs(jetstuff.Jet_lev6_noEMobj[i].Eta()-miscstuff.myCorrPhoton[j].Eta());
				double dR=sqrt(dPhi*dPhi+dEta*dEta);
				if (dPhi < dPhiMin) dPhiMin=dPhi;
				if (dEta < dEtaMin) dEtaMin=dEta;
				if (dR < dRMin) dRMin=dR;
				if (dPhi < MinDeltaPhigj || dPhi > MaxDeltaPhigj) return 0;
				if (dEta < MinDeltaEtagj || dEta > MaxDeltaEtagj) return 0;
				if (dR < MinDeltaRgj || dR > MaxDeltaRgj) return 0;
	    	} //for
		} //if 
	} //for
	
  if (dPhiMin < MinDeltaPhigjMin || dPhiMin > MaxDeltaPhigjMin) return 0;
  if (dEtaMin < MinDeltaEtagjMin || dEtaMin > MaxDeltaEtagjMin) return 0;
  if (dRMin < MinDeltaRgjMin || dRMin > MaxDeltaRgjMin) return 0;
  return 1;
}
/*}}}*/
 
//--------------------------------------------------------------------------
int TMyJetFilterModule::MyDuplicateCut(JetStuff jetstuff, CommonStuff miscstuff)
/*{{{*/
{
  double dRMin=0.2;
  int passcode=0;
  if(fRemoveDuplicate!=0)
    {
      if(bad_EMjet_match_flag==0) return 1;
      for(int j=0; j<miscstuff.myCorrPhoton.size(); j++)
 	{
 	  for(int i=0; i<miscstuff.myCorrElectron.size(); i++)
 	    {
 	      if(fabs(miscstuff.myEleEtaDet[i])<1.2)
 		{
 		  double dPhi=fabs(TVector2::Phi_mpi_pi(miscstuff.myElePhiDet[i]-miscstuff.myCorrPhoton[j].Phi()));
 		  double dEta=fabs(miscstuff.myEleEtaDet[i]-miscstuff.myPhoEtaDet[j]);
 		  double dR=sqrt(dPhi*dPhi+dEta*dEta);
 		  if(dR<dRMin) return 1;
 		}
 	    }
 	}
    }
  return passcode;
}
/*}}}*/ 

//--------------------------------------------------------------------------
int TMyJetFilterModule::MyMassCut(JetStuff jetstuff, CommonStuff miscstuff)
/*{{{*/
{
  double Mj1g1=0.0;
  double Mj1g2=0.0;
  double Mj2g1=0.0;
  double Mj2g2=0.0;
  double Mjj=0.0;
  int jetid=-1;
  int ind1=-1;
  int ind2=-1;
  for(int i=0; i<jetstuff.myNjet; i++)
    {
      if(fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())>MinJetEta && 
	 fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<MaxJetEta &&
	 fabs(jetstuff.Jet_lev6_noEMobj[i].Pt())>MinEtThr && 
	 fabs(jetstuff.Jet_lev6_noEMobj[i].Pt())<MaxEtThr)
	{
	  jetid++;
	  if(jetid==0) ind1=i;
	  if(jetid==1) ind2=i;
	  // assuming that there are at least two photons
	  for(int j=0; j<2; j++)
	    {
	      double M=(jetstuff.Jet_lev6_noEMobj[i]+miscstuff.myCorrPhoton[j]).M();
	      if(jetid==0 && j==0) Mj1g1=M;
	      if(jetid==0 && j==1) Mj1g2=M;
	      if(jetid==1 && j==0) Mj2g1=M;
	      if(jetid==1 && j==1) Mj2g2=M;
	      if(M<MinMjg || M>MaxMjg) return 0;
	    }
	} 
    }
  if(jetid>-1)
    {
      if(Mj1g1<MinMj1g1 || Mj1g1>MaxMj1g1) return 0;
      if(Mj1g2<MinMj1g2 || Mj1g2>MaxMj1g2) return 0;
      if(jetid>0) 
	{
	  Mjj=(jetstuff.Jet_lev6_noEMobj[ind1]+jetstuff.Jet_lev6_noEMobj[ind2]).M();
	  if(Mj2g1<MinMj2g1 || Mj2g1>MaxMj2g1) return 0;
	  if(Mj2g2<MinMj2g2 || Mj2g2>MaxMj2g2) return 0;
	  if(Mjj<MinMjj || Mjj>MaxMjj) return 0;
	}
    }
  return 1;
}
/*}}}*/


//--------------------------------------------------------------------------
int TMyJetFilterModule::MyAnalysisCut(JetStuff jetstuff, CommonStuff miscstuff, int metcode) {
  //--- temporarely commented out
//   if(MyQtJet(jetstuff)<MinQtJet || MyQtJet(jetstuff)>MaxQtJet) return 0;
//   if(MyHtAll(jetstuff,miscstuff,metcode)<MinHt || MyHtAll(jetstuff,miscstuff,metcode)>MaxHt) return 0;
  return 1;
}


//----- Test version as of 04/04/07
//--------------------------------------------------------------------------
// ATTN!!! corrected jet & met are used in this cut
//
// This version of cuts can be used for generated MET in MetModel studies.
//
int TMyJetFilterModule::MyMetCleanUpCut(JetStuff jetstuff, CommonStuff miscstuff,
					std::vector<TLorentzVector> jetvec, TVector2 metvec)
/*{{{*/
{
  int Nbadmetjet=0;
  if(abs(fRemoveBadMet)>0) // new as of 04/04/07
    {
      std::vector<double> phi_ObjInCrack;
      phi_ObjInCrack.clear();
      double MetEM_dPhi_max=0.3;
//       //______________ Bad MET due to mis-measured photons (new part)
//       for(int i=0; i<miscstuff.myCorrPhoton.size(); i++)
// 	{
// 	  int pho_in_twr9=0;
// 	  int pho_in_chimney=0;
// 	  if(fabs(miscstuff.myPhoZces[i])<9.0 || fabs(miscstuff.myPhoZces[i])>215.0) 
// 	    {
// 	      pho_in_twr9=1;
// // 	      phi_ObjInCrack.push_back(TVector2::Phi_0_2pi(miscstuff.myCorrPhoton[i].Phi()));
// 	    }
// // 	  if(miscstuff.myPhoEtaDet[i]>0.77 
// // 	     && miscstuff.myPhoEtaDet[i]<1.0
// // 	     && TVector2::Phi_0_2pi(miscstuff.myCorrPhoton[i].Phi())>1.30
// // 	     && TVector2::Phi_0_2pi(miscstuff.myCorrPhoton[i].Phi())>1.571) 
// // 	    {
// // 	      pho_in_chimney=1;
// // 	      phi_ObjInCrack.push_back(TVector2::Phi_0_2pi(miscstuff.myCorrPhoton[i].Phi()));
// // 	    }
// 	  if((pho_in_twr9==1 || pho_in_chimney==1) && metvec.Mod()>0.0) 
// 	    {
// 	      double dPhi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(miscstuff.myCorrPhoton[i].Phi())-
// 						    TVector2::Phi_0_2pi(metvec.Phi())));
// 	      if(dPhi<MetEM_dPhi_max) Nbadmetjet++;
// 	    }
// 	}
//       //______________ Bad MET due to mis-measured electrons (new 03/07/07)
//       for(int i=0; i<miscstuff.myCorrElectron.size(); i++)
// 	{
// 	  int ele_in_twr9=0;
// 	  if(fabs(miscstuff.myEleEtaDet[i])>0.95 
// 	     && fabs(miscstuff.myEleEtaDet[i])<1.2) 
// 	    {
// 	      ele_in_twr9=1;
// // 	      phi_ObjInCrack.push_back(TVector2::Phi_0_2pi(miscstuff.myElePhiDet[i]));
// 	    }
// 	  if(ele_in_twr9==1 && metvec.Mod()>0.0) // added on 03/08/07
// 	    {
// 	      double dPhi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(miscstuff.myElePhiDet[i])-
// 						    TVector2::Phi_0_2pi(metvec.Phi())));
// 	      if(dPhi<MetEM_dPhi_max) Nbadmetjet++;
// 	    }
// 	}
      
      //______________ Bad MET due to mis-measured jets
      for(int i=0; i<jetvec.size(); i++)
	{
	  //       if(jetvec[i].Pt()>=15.0) // removing bad MET due to jet in chimney
	  // 	{
	  // 	  int jet_in_chimney=0;
	  // 	  if(jetstuff.EtaDetCorr[i]>0.5 
	  // 	     && jetstuff.EtaDetCorr[i]<1.0 
	  // 	     && TVector2::Phi_0_2pi(jetvec[i].Phi())>1.05 
	  // 	     && TVector2::Phi_0_2pi(jetvec[i].Phi())<1.66) jet_in_chimney=1; // jet in chimney (same as in LED jet+MET search)
	  // 	  if(metvec.Mod()>0.0 && jet_in_chimney==1) // chimney
	  // 	    {
	  // 	      double dPhi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetvec[i].Phi())-
	  // 						    TVector2::Phi_0_2pi(metvec.Phi())));
	  // 	      if(dPhi>MinDeltaPhiJMet && dPhi<MaxDeltaPhiJMet) Nbadmetjet++; // jet is chimney and corrMET is along its direction 
	  // 	    }
	  // 	}
	  if(jetvec[i].Pt()>=3.0 && (jetstuff.Npho_match[i]+jetstuff.Nele_match[i])==0) // new cut (03/09/07)
	    {
	      int jet_in_crack=0;
	      int jet_in_lastTwr=0;
	      int jet_in_largeEmFr=0;
	      if(jetstuff.EmFrRaw[i]>0.875) jet_in_largeEmFr=1;
	      if(fabs(jetstuff.EtaDetCorr[i])>2.6) 
		{
		  jet_in_lastTwr=1;
		  phi_ObjInCrack.push_back(TVector2::Phi_0_2pi(jetvec[i].Phi()));
		}
	      if(fabs(jetstuff.EtaDetCorr[i])<0.2) 
		{
		  jet_in_crack=1;
		  phi_ObjInCrack.push_back(TVector2::Phi_0_2pi(jetvec[i].Phi()));
		}
	      else if(fabs(jetstuff.EtaDetCorr[i])>0.9 && fabs(jetstuff.EtaDetCorr[i])<1.3) 
		{
		  jet_in_crack=1;
		  phi_ObjInCrack.push_back(TVector2::Phi_0_2pi(jetvec[i].Phi()));
		}
	      if(metvec.Mod()>0.0) 
		{
		  double dPhi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetvec[i].Phi())-
							TVector2::Phi_0_2pi(metvec.Phi())));
		  if(jet_in_crack==1 && dPhi<MaxDeltaPhiJMet) Nbadmetjet++; // jet in crack and corrMET is along its direction
		  if(jet_in_lastTwr==1 && dPhi<0.4) Nbadmetjet++; // jet is in last tower
		  if(jet_in_largeEmFr==1 && dPhi<0.1) Nbadmetjet++; // mismeasured EM object
		}
	    }
	}
      double min_phi=100.0;
      double max_phi=-1.0;
      for(int i=0; i<phi_ObjInCrack.size(); i++)
	{
	  if(phi_ObjInCrack[i]>max_phi) max_phi=phi_ObjInCrack[i];
	  if(phi_ObjInCrack[i]<min_phi) min_phi=phi_ObjInCrack[i];      
	}
      if(phi_ObjInCrack.size()>1)
	{
	  double dphi=fabs(TVector2::Phi_mpi_pi(max_phi-min_phi));
	  double dphi1=fabs(TVector2::Phi_mpi_pi(max_phi-TVector2::Phi_0_2pi(metvec.Phi())));
	  double dphi2=fabs(TVector2::Phi_mpi_pi(min_phi-TVector2::Phi_0_2pi(metvec.Phi())));      
	  if(dphi1<=dphi && dphi2<dphi) Nbadmetjet++; // at least two mismeasured objects in uninstrumented regions
	}
    }
  return Nbadmetjet;
}
/*}}}*/

//________________________________ returns Qt of all jets above Eta & Et thresholds
double TMyJetFilterModule::MyQtJet(JetStuff jetstuff)
/*{{{*/
{
  double qt=-1.0E6;
  TLorentzVector vec(0.0,0.0,0.0,0.0);
  for(int i=0; i<jetstuff.myNjet; i++)
    {
      if(fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())>MinJetEta && 
	 fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<MaxJetEta &&
	 fabs(jetstuff.Jet_lev6_noEMobj[i].Pt())>MinEtThr && 
	 fabs(jetstuff.Jet_lev6_noEMobj[i].Pt())<MaxEtThr)
	{
	  vec=vec+jetstuff.Jet_lev6_noEMobj[i];
	}
    }  
  if(vec.E()>0.0) qt=vec.Pt();
  return qt;
}
/*}}}*/

//________________________________ returns SumEt=jets+photons+leptons+taus+met (jets above Eta & Et thresholds)
//                                 metscenario= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
double TMyJetFilterModule::MyHtAll(JetStuff jetstuff, CommonStuff miscstuff,
						int metscenario)
/*{{{*/
{
  double ht=0.0;
  //_________________ contribution from jets
  for(int i=0; i<jetstuff.myNjet; i++)
    {
      if(fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())>MinJetEta && 
	 fabs(jetstuff.Jet_lev6_noEMobj[i].Eta())<MaxJetEta &&
	 fabs(jetstuff.Jet_lev6_noEMobj[i].Pt())>MinEtThr && 
	 fabs(jetstuff.Jet_lev6_noEMobj[i].Pt())<MaxEtThr) ht=ht+jetstuff.Jet_lev6_noEMobj[i].Pt();
    }  
  //_________________ contribution from photons
  for(int i=0; i<miscstuff.myCorrPhoton.size(); i++)
    {
      ht=ht+miscstuff.myCorrPhoton[i].Pt();
    }
  //_________________ contribution from electrons
  for(int i=0; i<miscstuff.myCorrElectron.size(); i++)
    {
      ht=ht+miscstuff.myCorrElectron[i].Pt();
    }
  //_________________ contribution from muons
  for(int i=0; i<miscstuff.myCorrMuon.size(); i++)
    {
      ht=ht+miscstuff.myCorrMuon[i].Pt();
    }
  //_________________ contribution from taus
  for(int i=0; i<miscstuff.myCorrTau.size(); i++)
    {
      ht=ht+miscstuff.myCorrTau[i].Pt();
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
/*}}}*/
