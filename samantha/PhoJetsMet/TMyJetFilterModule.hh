/* this is SAM's stuff */

#ifndef TMYJET_HH
#define TMYJET_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include <vector>
//#include "CLHEP/Vector/LorentzVector.h" 
#include "TLorentzVector.h"
#include "TMath.h"

#include <Stntuple/loop/TStnModule.hh>
//------declaration of blocks to be read
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include <Stntuple/obj/TStnJetBlock.hh>
#include "JetUser/JetEnergyCorrections.hh"
//----- need this for jet-EM matching
#include <Stntuple/obj/TCalDataBlock.hh>
#include <CalorGeometry/CalConstants.hh>
#include <Stntuple/obj/TStnPhotonBlock.hh>
#include <Stntuple/obj/TStnElectronBlock.hh>
#include <Stntuple/obj/TStnMetBlock.hh>
#include <Stntuple/obj/TStnVertexBlock.hh>

#endif

class TGammaJets; // Sam's module

class TMyJetFilterModule: public TStnModule {
public:
 
  //------------------- parameters for Jet Correction initialization
  struct CorInit {
    int level;					
    int nvx;
    int cone;
    int version;
    int sys;
    int Nrun;
    int imode;
  };

  //------------------- this is a general structure for jet params which depend on Jet Cone size
  struct JetStuff {
    TVector2 myMETcorr_th5; // my MEt(4) after corrections for jets (muons?), using jets with Et>5 GeV
    double mySumEtCorr_th5; // SumEt after corrections for photons, electrons, JetClu-0.4 jets (muons?), using jets with Et>5 GeV  
    double mySumEtJet_th5;  // SumEt of JetClu-0.4 jets (lev6), using jets with Et>5 GeV  
    TVector2 myMETcorr_th10; // my MEt(4) after corrections for jets (muons?), using jets with Et>10 GeV
    double mySumEtCorr_th10; // SumEt after corrections for photons, electrons, JetClu-0.4 jets (muons?), using jets with Et>10 GeV  
    double mySumEtJet_th10;  // SumEt of JetClu-0.4 jets (lev6), using jets with Et>10 GeV   
    TVector2 myMETcorr_th15; // my MEt(4) after corrections for jets (muons?), using jets with Et>15 GeV
    double mySumEtCorr_th15; // SumEt after corrections for photons, electrons, JetClu-0.4 jets (muons?), using jets with Et>15 GeV  
    double mySumEtJet_th15;  // SumEt of JetClu-0.4 jets (lev6), using jets with Et>15 GeV   
    TVector2 myMETcorr_th20; // my MEt(4) after corrections for jets (muons?), using jets with Et>20 GeV
    double mySumEtCorr_th20; // SumEt after corrections for photons, electrons, JetClu-0.4 jets (muons?), using jets with Et>20 GeV  
    double mySumEtJet_th20;  // SumEt of JetClu-0.4 jets (lev6), using jets with Et>20 GeV   
    //_____________________________________________________________
    //----- my new output params (old one's to be removed later)
    
    int myNjet; // number of clusters ("raw jets")
    int myNjet_th5; // number of level-6 jets with Et>5, after matching with EM object
    int myNjet_th10; // number of level-6 jets with Et>10, after matching with EM object
    int myNjet_th15; // number of level-6 jets with Et>15, after matching with EM object
    int myNjet_th20; // number of level-6 jets with Et>20, after matching with EM object
    int myNjet_th25; // number of level-6 jets with Et>25, after matching with EM object
    
    std::vector<TLorentzVector> Jet_raw;      // uncorrected jets 
    std::vector<TLorentzVector> Jet_lev1;     // jets corrected to level-1
    std::vector<TLorentzVector> Jet_lev4;     // jets corrected to level-4
    std::vector<TLorentzVector> Jet_lev5;     // jets corrected to level-5
    std::vector<TLorentzVector> Jet_lev6;     // jets corrected to level-6
    std::vector<TLorentzVector> Jet_lev7;     // jets corrected to level-6 ?sam? should be 7?
    std::vector<TLorentzVector> Jet_raw_noEMobj;      // uncorrected jets, matched EM object removed 
    std::vector<TLorentzVector> Jet_lev1_noEMobj;     // jets corrected to level-1, matched EM object removed
    std::vector<TLorentzVector> Jet_lev4_noEMobj;     // jets corrected to level-4, matched EM object removed
    std::vector<TLorentzVector> Jet_lev5_noEMobj;     // jets corrected to level-5, matched EM object removed
    std::vector<TLorentzVector> Jet_lev6_noEMobj;     // jets corrected to level-6, matched EM object removed
    std::vector<TLorentzVector> Jet_lev7_noEMobj;     // jets corrected to level-7, matched EM object removed
    std::vector<int> JetNtrk;  // number of tracks in jet (Ntrk as defined in Stntuple)
    std::vector<int> JetNtwr;  // number of towers in jet (Ntwr as defined in Stntuple)
    std::vector<double> EtaDet;   // jet detector eta (before removing matched EM objects)
    std::vector<double> EtaDetCorr;   // jet detector eta (after removing matched EM objects)
    std::vector<double> EmFrRaw;  // EM fraction (before removing matched EM objects)
    std::vector<double> EmFrCorr; // EM fraction (after removing matched EM objects)
    std::vector<int> Nobj_match; // total number of matched objects (pho,ele,mu,tau,btag)
    std::vector<int> Npho_match; // number of matched photons
    std::vector<int> Nele_match; // number of matched electrons
    std::vector<int> Nmu_match; // number of matched muons
    std::vector<int> Ntau_match; // number of matched tau's (hadronic???)
    std::vector<int> Nbtag_match; // number of matched btags    
    std::vector<int> JetBlockInd; // original jet index in JetBlock, need this after jet reordering    
  };

  struct CommonStuff { // params which do not depend on Jet Cone size
    std::vector<TLorentzVector> myRawPhoton;    // raw photons
    std::vector<TLorentzVector> myCorrPhoton;   // corrected photons (leakage, em scale,etc.)
    std::vector<int> myPhoInd; // original photon index
    std::vector<double> myPhoEmFr;      // EM energy fraction for photons
    std::vector<double> myPhoEtaDet;    // detector eta for photons
    std::vector<double> myPhoXces;    // X_ces for photons
    std::vector<double> myPhoZces;    // Z_ces for photons
    std::vector<TLorentzVector> myRawElectron;  // raw electrons
    std::vector<TLorentzVector> myCorrElectron; // corrected electrons (leakage, em scale,etc.)
    std::vector<int> myEleInd; // original electron index
    std::vector<double> myEleEmFr;      // EM energy fraction for electrons
    std::vector<double> myEleEtaDet;    // detector eta for electrons
    std::vector<double> myElePhiDet;    // detector eta for electrons
    std::vector<double> myElePprEt;    // Ppr Et for electrons
    std::vector<TLorentzVector> myRawMuon;      // raw muons
    std::vector<TLorentzVector> myCorrMuon;     // corrected muons (corrections???)
    std::vector<TLorentzVector> myRawTau;       // raw taus
    std::vector<TLorentzVector> myCorrTau;      // corrected taus  (corrections???, for consistency)
    std::vector<TLorentzVector> myRawBjet;      // raw b-jets
    std::vector<TLorentzVector> myCorrBjet;     // corrected b-jets  (b-specific corrections???, for consistency)
    TVector2 myMET_raw; // MEt(4) before any corrections
    TVector2 myMET0_raw; // MEt(0) before any corrections
    double mySumEt_raw; // SumEt before any corrections    
  };

  struct MatchStuff {
    int JetInd_match;
    int EmInd_match;
    int EmObjType_match; // pho=0, ele=1
  };
  
  struct MetResults {
    //_______ data
    int Met20Njet_dt[3][10]; // DATA: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
    int Met25Njet_dt[3][10]; // DATA: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
    int Met30Njet_dt[3][10]; // DATA: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
    int Met35Njet_dt[3][10]; // DATA: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
    int Met40Njet_dt[3][10]; // DATA: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
    int Met45Njet_dt[3][10]; // DATA: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
    int Met50Njet_dt[3][10]; // DATA: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
    int Met75Njet_dt[3][10]; // DATA: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
    int Met100Njet_dt[3][10]; // DATA: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
    int Met150Njet_dt[3][10]; // DATA: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25
    
    //__________ met model default
    int Met20Njet_def[3][10]; // Met Model Default: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
    int Met25Njet_def[3][10]; // Met Model Default: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
    int Met30Njet_def[3][10]; // Met Model Default: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
    int Met35Njet_def[3][10]; // Met Model Default: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
    int Met40Njet_def[3][10]; // Met Model Default: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
    int Met45Njet_def[3][10]; // Met Model Default: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
    int Met50Njet_def[3][10]; // Met Model Default: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
    int Met75Njet_def[3][10]; // Met Model Default: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
    int Met100Njet_def[3][10]; // Met Model Default: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
    int Met150Njet_def[3][10]; // Met Model Default: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25

    //__________ met model -Sigma
    int Met20Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
    int Met25Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
    int Met30Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
    int Met35Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
    int Met40Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
    int Met45Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
    int Met50Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
    int Met75Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
    int Met100Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
    int Met150Njet_dvm[3][10]; // Met Model -SigmaJet: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25
    
    //__________ met model +Sigma
    int Met20Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
    int Met25Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
    int Met30Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
    int Met35Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
    int Met40Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
    int Met45Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
    int Met50Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
    int Met75Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
    int Met100Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
    int Met150Njet_dvp[3][10]; // Met Model +SigmaJet: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25

    //__________ met model -Sigma
    int Met20Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
    int Met25Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
    int Met30Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
    int Met35Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
    int Met40Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
    int Met45Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
    int Met50Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
    int Met75Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
    int Met100Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
    int Met150Njet_dvmUn[3][10]; // Met Model -SigmaUncl: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25
    
    //__________ met model +Sigma
    int Met20Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
    int Met25Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
    int Met30Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
    int Met35Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
    int Met40Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
    int Met45Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
    int Met50Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
    int Met75Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
    int Met100Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
    int Met150Njet_dvpUn[3][10]; // Met Model +SigmaUncl: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25
    
    //_______ normalized predictions
    double Met20Njet_gen[3][10]; // Met Model prediction: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
    double Met25Njet_gen[3][10]; // Met Model prediction: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
    double Met30Njet_gen[3][10]; // Met Model prediction: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
    double Met35Njet_gen[3][10]; // Met Model prediction: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
    double Met40Njet_gen[3][10]; // Met Model prediction: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
    double Met45Njet_gen[3][10]; // Met Model prediction: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
    double Met50Njet_gen[3][10]; // Met Model prediction: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
    double Met75Njet_gen[3][10]; // Met Model prediction: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
    double Met100Njet_gen[3][10]; // Met Model prediction: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
    double Met150Njet_gen[3][10]; // Met Model prediction: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25

    double Met20Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
    double Met25Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
    double Met30Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
    double Met35Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
    double Met40Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
    double Met45Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
    double Met50Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
    double Met75Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
    double Met100Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
    double Met150Njet_genstat[3][10]; // stat. uncertainty of Met Model prediction: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25

    double Met20Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
    double Met25Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
    double Met30Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
    double Met35Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
    double Met40Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
    double Met45Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
    double Met50Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
    double Met75Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
    double Met100Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
    double Met150Njet_gensyst[3][10]; // syst. uncertainty of Met Model prediction: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25
  };

  struct JetGeneral_t {
    //___________________________ before cuts
    // jet parameters are not corrected for matched EM objects
    // "raw" jet is a raw jet
    // corrected jet is level6 (not level7!!!) 
    TH1F*  fEvntNjet_b;
    TH1F*  fEvntNjet5_b;
    TH1F*  fEvntNjet10_b;
    TH1F*  fEvntNjet15_b;
    TH1F*  fEvntNjet20_b;
    TH1F*  fEvntNjet25_b;

    TH1F*  fEvntdZ_b;          // dZ=Z_jet-Z_best (Z_best is defined in MyUtil)
    TH1F*  fEvntEt0_b[2];      // raw Et of two leading jets, 1st jet is the highest Et0 jet
    TH1F*  fEvntEt_b[2];       // corrected Et of two leading jets, 1st jet is the highest Et0 jet
    TH1F*  fEvntEtaDet_b[2];   // EtaDet of two leading jets
    TH1F*  fEvntEta_b[2];      // Eta or true rapidity (if MidPoint is used) of two leading jets
    TH1F*  fEvntPhi_b[2];      // Phi of two leading jets
    TH1F*  fEvntTheta_b[2];    // theta of two leading jets
    TH1F*  fEvntEmFr_b[2];     // EmFraction of two leading jets
    TH1F*  fEvntNTowers_b[2];  // NTowers in two leading jets
    TH1F*  fEvntNTracks_b[2];  // NTracks in two leading jets
    TH1F*  fEvntEt0X_b;        // This histo is for extra jets
    TH1F*  fEvntEtX_b;        // This histo is for extra jets
    TH1F*  fEvntEtaDetX_b;    
    TH1F*  fEvntEtaX_b;       
    TH1F*  fEvntPhiX_b;       
    TH1F*  fEvntThetaX_b;     
    TH1F*  fEvntEmFrX_b;      
    TH1F*  fEvntNTowersX_b;   
    TH1F*  fEvntNTracksX_b;   
    TH1F*  fEvntThetaStar_b;    // ThetaStar of two leading jets
    TH1F*  fEvntDeltaPhi_b;     // dPhi of two leading jets
    TH1F*  fEvntDeltaEta_b;     // dEta of two leading jets
    TH1F*  fEvntDeltaR_b;       // dR of two leading jets
    TH1F*  fEvntMjj_b;          // Mjj of two leading jets
    TH1F*  fEvntDeltaR1x_b;     // dR between 1st jet and any extra jet
    TH1F*  fEvntDeltaR2x_b;     // dR between 2nd jet and any extra jet
    TH1F*  fEvntKt2jet_b;       // Kt of 2 leading jets
    TH1F*  fEvntKtAll_b;        // Kt of all jets

    TProfile*  fEvntNJet_Nvx_b;          // Njet vs Nvx(class>12) 
    TProfile*  fEvntNJet5_Nvx_b;         // Njet5 vs Nvx(class>12) 
    TProfile*  fEvntNJet10_Nvx_b;        // Njet10 vs Nvx(class>12) 
    TProfile*  fEvntNJet15_Nvx_b;        // Njet15 vs Nvx(class>12) 
    TProfile*  fEvntNJet20_Nvx_b;        // Njet20 vs Nvx(class>12) 
    TProfile*  fEvntNJet25_Nvx_b;        // Njet25 vs Nvx(class>12) 

    TProfile*  fEvntKt_Mjj_b;            // "Kt-kick" of dijet system vs. Mjj
    TProfile*  fEvntKtAll_Mjj_b;         // "Kt" of extra vs. Mjj
    TProfile*  fEvntEt0_Njet_b;          // Average raw Et of jets vs. Njet 
    TProfile*  fEvntEt_Njet_b;           // Average corrected Et of jets vs. Njet 
    TProfile*  fEvntEt0_Nvx12_b;         // Average raw Et of jets vs. Nvx(class>=12) 
    TProfile*  fEvntEt_Nvx12_b;          // Average corrected Et of jets vs. Nvx(class>=12)                                          
    TProfile*  fEvntNjet_Lum_b;          // number of jets vs. Inst. lum for Nvx12=1

    TProfile*  fJetRes_MetSig_EtaDet_j5_b;  // resolution study:  corrMET(|dPhi|<0.4)/Sigma(jet) vs. EtaDet for Et(lev6)>5 
    TProfile*  fJetRes_MetSig_EtaDet_j10_b;  // resolution study:  corrMET(|dPhi|<0.4)/Sigma(jet) vs. EtaDet for Et(lev6)>10 
    TProfile*  fJetRes_MetSig_EtaDet_j15_b;  // resolution study:  corrMET(|dPhi|<0.4)/Sigma(jet) vs. EtaDet for Et(lev6)>15 
    TH1F* fJetRes_DetEta_j5_b; // resolution study: det eta of jets with Et(lev6)>5
    TH1F* fJetRes_DetEta_j10_b; // resolution study: det eta of jets with Et(lev6)>10
    TH1F* fJetRes_DetEta_j15_b; // resolution study: det eta of jets with Et(lev6)>15
    TH1F* fJetRes_DetEta_j5_Sig3dPhi04_b; // resolution study: det eta of jets with Et(lev6)>5 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j10_Sig3dPhi04_b; // resolution study: det eta of jets with Et(lev6)>10 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j15_Sig3dPhi04_b; // resolution study: det eta of jets with Et(lev6)>15 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j5_Sig3dPhi01_b; // resolution study: det eta of jets with Et(lev6)>5 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j10_Sig3dPhi01_b; // resolution study: det eta of jets with Et(lev6)>10 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j15_Sig3dPhi01_b; // resolution study: det eta of jets with Et(lev6)>15 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j5_FrSig3dPhi04_b; // resolution study: det eta of jets with Et(lev6)>5 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j10_FrSig3dPhi04_b; // resolution study: det eta of jets with Et(lev6)>10 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j15_FrSig3dPhi04_b; // resolution study: det eta of jets with Et(lev6)>15 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j5_FrSig3dPhi01_b; // resolution study: det eta of jets with Et(lev6)>5 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j10_FrSig3dPhi01_b; // resolution study: det eta of jets with Et(lev6)>10 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j15_FrSig3dPhi01_b; // resolution study: det eta of jets with Et(lev6)>15 and MetSig>3 if |dPhi|<0.1

    TH1F* fEvnt_toyMET_all; // "toy" MET for all events before cuts. 
    TH1F* fEvnt_toyMET_cut; // "toy" MET for events passing MetCleanup(using toyMET) cuts. 
    TH1F* fEvnt_corMET_all; // corrected MET before any MetCleanup cuts

    //______________________________ after cuts 
    TH1F*  fEvntNjet_a;
    TH1F*  fEvntNjet5_a;
    TH1F*  fEvntNjet10_a;
    TH1F*  fEvntNjet15_a;
    TH1F*  fEvntNjet20_a;
    TH1F*  fEvntNjet25_a;

    TH1F*  fEvntdZ_a;        
    TH1F*  fEvntEt0_a[2];    
    TH1F*  fEvntEt_a[2];     
    TH1F*  fEvntEtaDet_a[2]; 
    TH1F*  fEvntEta_a[2];    
    TH1F*  fEvntPhi_a[2];    
    TH1F*  fEvntTheta_a[2];  
    TH1F*  fEvntEmFr_a[2];   
    TH1F*  fEvntNTowers_a[2];
    TH1F*  fEvntNTracks_a[2];
    TH1F*  fEvntEt0X_a;      
    TH1F*  fEvntEtX_a;      
    TH1F*  fEvntEtaDetX_a;    
    TH1F*  fEvntEtaX_a;       
    TH1F*  fEvntPhiX_a;       
    TH1F*  fEvntThetaX_a;     
    TH1F*  fEvntEmFrX_a;      
    TH1F*  fEvntNTowersX_a;   
    TH1F*  fEvntNTracksX_a;   
    TH1F*  fEvntThetaStar_a; 
    TH1F*  fEvntDeltaPhi_a;  
    TH1F*  fEvntDeltaEta_a;  
    TH1F*  fEvntDeltaR_a;   
    TH1F*  fEvntMjj_a;  
    TH1F*  fEvntDeltaR1x_a; 
    TH1F*  fEvntDeltaR2x_a; 
    TH1F*  fEvntKt2jet_a;   
    TH1F*  fEvntKtAll_a;    

    TProfile*  fEvntNJet_Nvx_a;  
    TProfile*  fEvntNJet5_Nvx_a;         
    TProfile*  fEvntNJet10_Nvx_a;         
    TProfile*  fEvntNJet15_Nvx_a;         
    TProfile*  fEvntNJet20_Nvx_a;         
    TProfile*  fEvntNJet25_Nvx_a;         

    TProfile*  fEvntKt_Mjj_a;    
    TProfile*  fEvntKtAll_Mjj_a; 
    TProfile*  fEvntEt0_Njet_a;  
    TProfile*  fEvntEt_Njet_a;   
    TProfile*  fEvntEt0_Nvx12_a; 
    TProfile*  fEvntEt_Nvx12_a;  
    TProfile*  fEvntNjet_Lum_a;  

    TProfile*  fJetRes_MetSig_EtaDet_j5_a;  // resolution study:  corrMET(|dPhi|<0.4)/Sigma(jet) vs. EtaDet for Et(lev6)>5 
    TProfile*  fJetRes_MetSig_EtaDet_j10_a;  // resolution study:  corrMET(|dPhi|<0.4)/Sigma(jet) vs. EtaDet for Et(lev6)>10 
    TProfile*  fJetRes_MetSig_EtaDet_j15_a;  // resolution study:  corrMET(|dPhi|<0.4)/Sigma(jet) vs. EtaDet for Et(lev6)>15 
    TH1F* fJetRes_DetEta_j5_a; // resolution study: det eta of jets with Et(lev6)>5
    TH1F* fJetRes_DetEta_j10_a; // resolution study: det eta of jets with Et(lev6)>10
    TH1F* fJetRes_DetEta_j15_a; // resolution study: det eta of jets with Et(lev6)>15
    TH1F* fJetRes_DetEta_j5_Sig3dPhi04_a; // resolution study: det eta of jets with Et(lev6)>5 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j10_Sig3dPhi04_a; // resolution study: det eta of jets with Et(lev6)>10 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j15_Sig3dPhi04_a; // resolution study: det eta of jets with Et(lev6)>15 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j5_Sig3dPhi01_a; // resolution study: det eta of jets with Et(lev6)>5 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j10_Sig3dPhi01_a; // resolution study: det eta of jets with Et(lev6)>10 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j15_Sig3dPhi01_a; // resolution study: det eta of jets with Et(lev6)>15 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j5_FrSig3dPhi04_a; // resolution study: det eta of jets with Et(lev6)>5 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j10_FrSig3dPhi04_a; // resolution study: det eta of jets with Et(lev6)>10 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j15_FrSig3dPhi04_a; // resolution study: det eta of jets with Et(lev6)>15 and MetSig>3 if |dPhi|<0.4
    TH1F* fJetRes_DetEta_j5_FrSig3dPhi01_a; // resolution study: det eta of jets with Et(lev6)>5 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j10_FrSig3dPhi01_a; // resolution study: det eta of jets with Et(lev6)>10 and MetSig>3 if |dPhi|<0.1
    TH1F* fJetRes_DetEta_j15_FrSig3dPhi01_a; // resolution study: det eta of jets with Et(lev6)>15 and MetSig>3 if |dPhi|<0.1

  };

  struct MetStudyHisto_t {  //------------ Met study histograms 

    TH1F* fSumEtCorr;     // corrected SumEt
    TH1F* fSqrtSumEtCorr; // corrected sqrt(SumEt)

    TH1F* fSumEtCorr_withJet;     // corrected SumEt for events with jets above threshold
    TH1F* fSqrtSumEtCorr_withJet; // corrected sqrt(SumEt) for events with jets above threshold
    TH1F* fSumEtCorr_noJet;     // corrected SumEt for events without jets above threshold
    TH1F* fSqrtSumEtCorr_noJet; // corrected sqrt(SumEt) for events without jets above threshold

    TH1F* fSumEtCorr_noJet_vx1;     // corrected SumEt for events without jets above threshold and Nvx=1
    TH1F* fSumEtCorr_noJet_vx2;     // corrected SumEt for events without jets above threshold and Nvx=2
    TH1F* fSumEtCorr_noJet_vx3;     // corrected SumEt for events without jets above threshold and Nvx=3
    TProfile* fSumEtCorrNoJet_vs_Nvx; // corrected SumEt for events without jets vs. Nvx


    TH1F* fMetRaw;  // for events where met was corrected
    TH1F* fMetCorr; // for events where met was corrected
    TH1F* fMetPhiRaw;  // for events where met was corrected
    TH1F* fMetPhiCorr; // for events where met was corrected

    TH1F* fMetRawAll;  // for all events (whether or not met was corrected)
    TH1F* fMetCorrAll; // for all events (whether or not met was corrected)
    TH1F* fMetRawAll_X;  // for all events (whether or not met was corrected)
    TH1F* fMetCorrAll_X; // for all events (whether or not met was corrected)
    TH1F* fMetRawAll_Y;  // for all events (whether or not met was corrected)
    TH1F* fMetCorrAll_Y; // for all events (whether or not met was corrected)
    TH1F* fMetPhiRawAll;  // for all events (whether or not met was corrected)
    TH1F* fMetPhiCorrAll; // for all events (whether or not met was corrected)
    
    TProfile* fMet4VsSqrtCorrSumEt; // Met4 vs. Sqrt(corr SumEt) for all events  
    TProfile* fMet4VsCorrSumEt;     // Met4 vs. (corr SumEt) for all events  
    
    TProfile* fMet4VsSqrtCorrSumEt_nvx1; // Met4 vs. Sqrt(corr SumEt) for Nvx12=1  
    TProfile* fMet4VsSqrtCorrSumEt_nvx2; // Met4 vs. Sqrt(corr SumEt) for Nvx12=2  
    TProfile* fMet4VsSqrtCorrSumEt_nvx3; // Met4 vs. Sqrt(corr SumEt) for Nvx12=3  
    TProfile* fMet4VsSqrtCorrSumEt_nvx4; // Met4 vs. Sqrt(corr SumEt) for Nvx12=4  
    TProfile* fMet4VsSqrtCorrSumEt_nvx5; // Met4 vs. Sqrt(corr SumEt) for Nvx12=5  

    TProfile* fMet4VsNjet; // Met4 vs. Njet (number of clusters)
    TProfile* fMet4VsNjet_th5; // Met4 vs. Njet(lev6>5GeV)
    TProfile* fMet4VsNjet_th10; // Met4 vs. Njet(lev6>10GeV)
    TProfile* fMet4VsNjet_th15; // Met4 vs. Njet(lev6>15GeV)
    TProfile* fMet4VsNjet_th20; // Met4 vs. Njet(lev6>20GeV)
    TProfile* fMet4VsNjet_th25; // Met4 vs. Njet(lev6>25GeV)

    TProfile* fMet4SigVsNjet; // Met4/sqrt(SumEt) vs. Njet (number of clusters)
    TProfile* fMet4SigVsNjet_th5; // Met4/sqrt(SumEt) vs. Njet(lev6>5GeV)
    TProfile* fMet4SigVsNjet_th10; // Met4/sqrt(SumEt) vs. Njet(lev6>10GeV)
    TProfile* fMet4SigVsNjet_th15; // Met4/sqrt(SumEt) vs. Njet(lev6>15GeV)
    TProfile* fMet4SigVsNjet_th20; // Met4/sqrt(SumEt) vs. Njet(lev6>20GeV)
    TProfile* fMet4SigVsNjet_th25; // Met4/sqrt(SumEt) vs. Njet(lev6>25GeV)

    TProfile* fMet4VsRun; // Met4 vs. Run number, for nvx12=1
    TProfile* fMet4SigVsRun; // Met4/sqrt(SumEt) vs. Run number, for nvx12=1
    TProfile* fMet4VsNvx12; // Met4 vs. Nvx12
    TProfile* fMet4SigVsNvx12; // Met4/sqrt(SumEt) vs. Nvx12
    TProfile* fNvx12VsRun; // Nvx12 vs. Run number

    TH1F* fMet4[10]; // Met4 for events with any Nvx12 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
    TH1F* fMet4Phi[10]; // Phi of Met4 for events with any Nvx12 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
    
    TH1F* fMet4X[10]; // Met4X for events with any Nvx12 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
    TH1F* fMet4X_nvx1[10]; // Met4X for Nvx12=1 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4X_nvx2[10]; // Met4X for Nvx12=2 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4X_nvx3[10]; // Met4X for Nvx12=3 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4X_nvx4[10]; // Met4X for Nvx12=4 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4X_nvx5[10]; // Met4X for Nvx12=5 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4X_nvxM1[10]; // Met4X for Nvx12>1 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 

    TH1F* fMet4Y[10]; // Met4Y for events with any Nvx12 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
    TH1F* fMet4Y_nvx1[10]; // Met4X for Nvx12=1 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4Y_nvx2[10]; // Met4X for Nvx12=2 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4Y_nvx3[10]; // Met4X for Nvx12=3 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4Y_nvx4[10]; // Met4X for Nvx12=4 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4Y_nvx5[10]; // Met4X for Nvx12=5 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 
    TH1F* fMet4Y_nvxM1[10]; // Met4X for Nvx12>1 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0 

    TH1F* fMet4X_noJet[10]; // Met4X for events with Njet_thr*=0 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
    TH1F* fMet4X_withJet[10]; // Met4X for events with Njet_thr*>0 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
    TH1F* fMet4Y_noJet[10]; // Met4Y for events with Njet_thr*=0 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
    TH1F* fMet4Y_withJet[10]; // Met4Y for events with Njet_thr*>0 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0

    TH1F* fSumEtJetFrac; // fraction of SumEt carried by jets with Et>threshold
    TH1F* fSumEtJet; // SumEt carried by jets with Et>threshold
    TH1F* fSqrtSumEtJet; // sqrt(SumEt) carried by jets with Et>threshold
    TProfile* fJetFrVsSumEt; // fraction of SumEt carried by jets with Et>threshold vs. total SumEt
    TProfile* fMet4VsJetFr; // meth4 Met .vs. fraction of SumEt carried by jets with Et>threshold
    TProfile* fMet4SigVsJetFr; // meth4 Met/sqrt(SumEt) .vs. fraction of SumEt carried by jets with Et>threshold
    TProfile* fMet4VsSqrtSumEtJet; // meth4 Met .vs. sqrt(SumEtJet) (SumEt carried by jets with Et>threshold)
    TProfile* fMet4SigVsSqrtSumEtJet; // meth4 Met/sqrt(SumEt) .vs. sqrt(SumEtJet) (SumEt carried by jets with Et>threshold)
    TH2F* fMet4X_vs_Met4Y; // meth4 Met_X vs. Met_Y for all events
    TH2F* fGenMetX_vs_MetY; // generated (def. parametrization) Met_X vs. Met_Y for all events

    //__________________ histograms for MetModel2
    TProfile* fMetGenV2VsNjet; // generated (def parametrization) Met vs. Njet (number of clusters)
    TProfile* fMetGenV2VsNjet_th5; // generated (def parametrization) Met vs. Njet(lev6>5GeV)
    TProfile* fMetGenV2VsNjet_th10; // generated (def parametrization) Met vs. Njet(lev6>10GeV)
    TProfile* fMetGenV2VsNjet_th15; // generated (def parametrization) Met vs. Njet(lev6>15GeV)
    TProfile* fMetGenV2VsNjet_th20; // generated (def parametrization) Met vs. Njet(lev6>20GeV)
    TProfile* fMetGenV2VsNjet_th25; // generated (def parametrization) Met vs. Njet(lev6>25GeV)

    TProfile* fMetGenV2SigVsNjet; // generated (def parametrization) Met/sqrt(SumEt) vs. Njet (number of clusters)
    TProfile* fMetGenV2SigVsNjet_th5; // generated (def parametrization) Met/sqrt(SumEt) vs. Njet(lev6>5GeV)
    TProfile* fMetGenV2SigVsNjet_th10; // generated (def parametrization) Met/sqrt(SumEt) vs. Njet(lev6>10GeV)
    TProfile* fMetGenV2SigVsNjet_th15; // generated (def parametrization) Met/sqrt(SumEt) vs. Njet(lev6>15GeV)
    TProfile* fMetGenV2SigVsNjet_th20; // generated (def parametrization) Met/sqrt(SumEt) vs. Njet(lev6>20GeV)
    TProfile* fMetGenV2SigVsNjet_th25; // generated (def parametrization) Met/sqrt(SumEt) vs. Njet(lev6>25GeV)

    TProfile* fMetGenV2VsRun; // generated (def parametrization) Met vs. Run number, for nvx12=1
    TProfile* fMetGenV2SigVsRun; // generated (def parametrization) Met/sqrt(SumEt) vs. Run number, for nvx12=1
    TProfile* fMetGenV2VsNvx12; // generated (def parametrization) Met vs. Nvx12
    TProfile* fMetGenV2SigVsNvx12; // generated (def parametrization) Met/sqrt(SumEt) vs. Nvx12
    
    TH1F* fGenV2Met_def;    // generated MET using default parametrization
    TH1F* fGenV2MetX_def;    // generated MET_X using default parametrization
    TH1F* fGenV2MetY_def;    // generated MET_X using default parametrization
    TH1F* fGenV2MetPhi_def; // PhiMET calculated using generated MET and default parametrization

    TH1F* fGenV2Met_devm;    // generated MET using deviated (-sigma) jet parametrization
    TH1F* fGenV2MetX_devm;    // generated MET_X using deviated (-sigma) jet parametrization
    TH1F* fGenV2MetY_devm;    // generated MET_X using deviated (-sigma) jet parametrization
    TH1F* fGenV2MetPhi_devm; // PhiMET calculated using generated MET and deviated (-sigma) jet parametrization

    TH1F* fGenV2Met_devp;    // generated MET using deviated (+sigma) jet parametrization
    TH1F* fGenV2MetX_devp;    // generated MET_X using deviated (+sigma) jet parametrization
    TH1F* fGenV2MetY_devp;    // generated MET_X using deviated (+sigma) jet parametrization
    TH1F* fGenV2MetPhi_devp; // PhiMET calculated using generated MET and deviated (+sigma) jet parametrization

    TH1F* fGenV2Met_devmUn;    // generated MET using deviated (-sigma) unclustered parametrization
    TH1F* fGenV2MetX_devmUn;    // generated MET_X using deviated (-sigma) unclustered parametrization
    TH1F* fGenV2MetY_devmUn;    // generated MET_X using deviated (-sigma) unclustered  parametrization
    TH1F* fGenV2MetPhi_devmUn; // PhiMET calculated using generated MET and deviated (-sigma) unclustered parametrization

    TH1F* fGenV2Met_devpUn;    // generated MET using deviated (+sigma) unclustered parametrization
    TH1F* fGenV2MetX_devpUn;    // generated MET_X using deviated (+sigma) unclustered parametrization
    TH1F* fGenV2MetY_devpUn;    // generated MET_X using deviated (+sigma) unclustered parametrization
    TH1F* fGenV2MetPhi_devpUn; // PhiMET calculated using generated MET and deviated (+sigma) unclustered parametrization

    TProfile* fMetGenV2VsJetFr; // generated (def. parametrization) Met .vs. fraction of SumEt carried by jets with Et>threshold
    TProfile* fMetGenV2SigVsJetFr; // generated (def. parametrization) Met/sqrt(SumEt) .vs. fraction of SumEt carried by jets with Et>threshold
    TProfile* fMetGenV2VsSqrtSumEtJet; // generated (def. parametrization) Met .vs. sqrt(SumEtJet) (SumEt carried by jets with Et>threshold)
    TProfile* fMetGenV2SigVsSqrtSumEtJet; // generated (def. parametrization) Met/sqrt(SumEt) .vs. sqrt(SumEtJet) (SumEt carried by jets with Et>threshold)
    TH2F* fGenV2MetX_vs_MetY; // generated (def. parametrization) Met_X vs. Met_Y for all events

    TH1F* fMet4GenV2Met; // Met4-MetGen
    TH2F* fMet4_vs_GenV2Met; // Met4 .vs. MetGen
    TH2F* fMet4_vs_Met4GenV2Met; // Met4 .vs. Met4-MetGen

  };

  struct MatchStudyHisto_t {  //------------ histograms for study of Pho-Jet  
    TH1F* fMatchNtwr; // number of towers in matched jet (1st match)
    TH1F* fMatchDelR; // angular distance between jet and matched photon (1st match)
    TH1F* fMatchDelPhi; // dPhi between jet and matched photon (1st match)
    TH1F* fMatchDelEta; // dEta between jet and matched photon (1st match)
    TH1F* fMatchDelEtaDet; // dEtaDet between jet and matched photon (1st match)
    TH1F* fMatchNmatch; // number of matched objects
    TH1F* fMatchEt_raw_b; // raw Et of jet before matching
    TH1F* fMatchEt_raw_a; // raw Et of jet after matching
    TH1F* fMatchEt_lev6_b; // lev6 Et of jet before matching
    TH1F* fMatchEt_lev6_a; // lev6 Et of jet after matching
    TH1F* fMatchEtJet2EtPho_b; // ratio Et(raw jet)/Et(raw pho) of jet before matching
    TH1F* fMatchEtJet2EtPho_a; // ratio Et(raw jet)/Et(raw pho) of jet after matching
    TH1F* fMatchDelRoldnew; // angular distance between "old" and "new" jet
    TH1F* fMatchDelEtaDetoldnew; // dEtaDet between "old" and "new" jet
  };


  struct MetCleanupHisto_t {  //------------ histograms for met cleanup studies

    TH1F* fCleanupDelPhi_metjet[3]; // dPhi=Phi(met)-Phi(jet) in 3 regions: 
                                    //      [0]-- if jet is close to cracks (|dEta|<0.2); [1]-- if jet is away from cracks
                                    //      [2]-- if jet is in last towers (|eta|>2.6)  
    TH1F* fCleanupMetEtjet[3]; // MET/Et(jet) in 3 regions: 
                                    //      [0]-- if jet is close to cracks (|dEta|<0.2); [1]-- if jet is away from cracks
                                    //      [2]-- if jet is in last towers (|eta|>2.6)  
    TH2F* fCleanupMetEtjetVsDelPhi[3]; // MET/Et(jet).vs.DelPhi(met-jet) in 3 regions: 
                                    //      [0]-- if jet is close to cracks (|dEta|<0.2); [1]-- if jet is away from cracks
                                    //      [2]-- if jet is in last towers (|eta|>2.6)  
    TH2F* fCleanupDelPhiVsEtaDet_jet; // DelPhi(met-jet).vs.dPhi=Phi(met)-Phi(jet)

    TH2F* fCleanupDelPhiVsEtaDet_em; // DelPhi(met-em).vs.dPhi=Phi(met)-Phi(em)
    TH1F* fCleanupDelPhi_metem; // dPhi=Phi(met)-Phi(em)
    TH1F* fCleanupMetEtem; // MET/Et(met) 
    TH2F* fCleanupMetEtemVsDelPhi; // MET/Et(em).vs.DelPhi(met-em)

    TH1F* fCleanupDelPhi_RawMetRawJet3; // dPhi=Phi(rawMet)-Phi(rawJet), raw Et>3
    TH1F* fCleanupDelPhi_RawMetRawJet5; // dPhi=Phi(rawMet)-Phi(rawJet), raw Et>5
    TH1F* fCleanupDelPhi_RawMetRawJet10; // dPhi=Phi(rawMet)-Phi(rawJet), raw Et>10
    TH1F* fCleanupDelPhi_RawMet1stRawJet; // dPhi=Phi(rawMet)-Phi(1st rawJet) raw Et>3

    TH1F* fCleanupDelPhi_metjet5; // dPhi=Phi(met)-Phi(jet), lev6 Et(jet)>5, all jets after MetCleanup cuts 
    TH1F* fCleanupDelPhi_metjet15; // dPhi=Phi(met)-Phi(jet), lev6 Et(jet)>15, all jets after MetCleanup cuts 
    TH1F* fCleanupDelPhi_metjet25; // dPhi=Phi(met)-Phi(jet), lev6 Et(jet)>25, all jets after MetCleanup cuts 
    TH1F* fCleanupDelPhi_met10jet15; // dPhi=Phi(met)-Phi(jet), lev6 Et(jet)>15, corMET>10, all jets after MetCleanup cuts 
    TH1F* fCleanupDelPhi_metjet1st15; // dPhi=Phi(met)-Phi(jet), highest lev6 Et(jet)>15 after MetCleanup cuts 
    TH1F* fCleanupDelPhi_met10jet1st15; // dPhi=Phi(met)-Phi(jet), highest lev6 Et(jet)>15 after MetCleanup cuts, corMET>10 
    TH1F* fCleanupDelPhi_metjet15_dPhiMin; // dPhi=Phi(met)-Phi(jet), closest lev6 Et(jet)>15 after MetCleanup cuts 
    TH1F* fCleanupDelPhi_met10jet15_dPhiMin; // dPhi=Phi(met)-Phi(jet), closest lev6 Et(jet)>15 after MetCleanup cuts, corMET>10 

  };


  struct MetProbHisto_t { //----- histograms for Met probability study;

    TH2F* fDataMet_vs_GenSigmaMet; // Met(data) vs. SigmaMet(generated) 
    TH2F* fToyEvntMet_vs_GenSigmaMet; // measured Met(Toy) vs. SigmaMet(generated) 

    TH1F* fDataMet_highLnProb; // met in data events with -log10(1-P)>2.56 (>~3sigma fluctuation)
    TH1F* fDataMet_highSigma; // met in data events with (Met-meanGenMet)/sigmaGenMet>3.0 (>~3sigma fluctuation)
    TH1F* fDataMet_lowLnProb; // met in data events with -log10(1-P)<2.56 (<~3sigma fluctuation)
    TH1F* fDataMet_lowSigma; // met in data events with (Met-meanGenMet)/sigmaGenMet<3.0 (<~3sigma fluctuation)
    TH1F* fGenMet_def_highLnProb; // Generated def. met for data events where met -log10(1-P)>2.56 (>~3sigma fluctuation)
    TH1F* fGenMet_def_highSigma; // Generated def. met for data events where met (Met-meanGenMet)/sigmaGenMet>3.0 (>~3sigma fluctuation)
    TH1F* fGenMet_def_lowLnProb; // Generated def. met for data events where met -log10(1-P)<2.56 (<~3sigma fluctuation)
    TH1F* fGenMet_def_lowSigma; // Generated def. met for data events where met (Met-meanGenMet)/sigmaGenMet<3.0 (<~3sigma fluctuation)
    TH1F* fGenMet_dvm_highLnProb; // Generated "-sigma" met for data events where met -log10(1-P)>2.56 (>~3sigma fluctuation)
    TH1F* fGenMet_dvm_highSigma; // Generated "-sigma" met for data events where met (Met-meanGenMet)/sigmaGenMet>3.0 (>~3sigma fluctuation)
    TH1F* fGenMet_dvm_lowLnProb; // Generated "-sigma" met for data events where met -log10(1-P)<2.56 (<~3sigma fluctuation)
    TH1F* fGenMet_dvm_lowSigma; // Generated "-sigma" met for data events where met (Met-meanGenMet)/sigmaGenMet<3.0 (<~3sigma fluctuation)
    TH1F* fGenMet_dvp_highLnProb; // Generated "+sigma" met for data events where met -log10(1-P)>2.56 (>~3sigma fluctuation)
    TH1F* fGenMet_dvp_highSigma; // Generated "+sigma" met for data events where met (Met-meanGenMet)/sigmaGenMet>3.0 (>~3sigma fluctuation)
    TH1F* fGenMet_dvp_lowLnProb; // Generated "+sigma" met for data events where met -log10(1-P)<2.56 (<~3sigma fluctuation)
    TH1F* fGenMet_dvp_lowSigma; // Generated "+sigma" met for data events where met (Met-meanGenMet)/sigmaGenMet<3.0 (<~3sigma fluctuation)

    TH1F* fDataMetProb_def;  // Data: -log10(1-P) for default parametrization 
    TH1F* fDataMetProb_dvm;  // Data: -log10(1-P) for -sigma parametrization
    TH1F* fDataMetProb_dvp;  // Data: -log10(1-P) for +sigma parametrization
    TH1F* fToyMetProb_def[20];  // Toy MC in bins of toyMet: -log10(1-P) for default parametrization 
    TH1F* fToyMetProb_dvm[20];  // Toy MC in bins of toyMet: -log10(1-P) for -sigma parametrization
    TH1F* fToyMetProb_dvp[20];  // Toy MC in bins of toyMet: -log10(1-P) for +sigma parametrization
    TProfile* fDataMetProb_Met_def; // Data: -log10(1-P).vs.Met for default parametrization
    TProfile* fDataMetProb_Met_dvm; // Data: -log10(1-P).vs.Met for -sigma parametrization
    TProfile* fDataMetProb_Met_dvp; // Data: -log10(1-P).vs.Met for +sigma parametrization
    TProfile* fToyMetProb_toyMet_def; // Toy MC: -log10(1-P).vs.toyMet for default parametrization
    TProfile* fToyMetProb_toyMet_dvm; // Toy MC: -log10(1-P).vs.toyMet for -sigma parametrization
    TProfile* fToyMetProb_toyMet_dvp; // Toy MC: -log10(1-P).vs.toyMet for +sigma parametrization
    TH1F* fDataMetSigma_def; // Data: (Met-meanGenMet)/sigmaGenMet for default parametrization
    TH1F* fDataMetSigma_dvm; // Data: (Met-meanGenMet)/sigmaGenMet for -sigma parametrization
    TH1F* fDataMetSigma_dvp; // Data: (Met-meanGenMet)/sigmaGenMet for -sigma parametrization
    TH1F* fToyMetSigma_def[20]; // Toy MC in bins of toyMet: (toyMet-meanGenMet)/sigmaGenMet for default parametrization
    TH1F* fToyMetSigma_dvm[20]; // Toy MC in bins of toyMet: (toyMet-meanGenMet)/sigmaGenMet for -sigma parametrization
    TH1F* fToyMetSigma_dvp[20]; // Toy MC in bins of toyMet: (toyMet-meanGenMet)/sigmaGenMet for -sigma parametrization
    TProfile* fDataMetSigma_Met_def; // Data: (Met-meanGenMet)/sigmaGenMet.vs.Met for default parametrization
    TProfile* fDataMetSigma_Met_dvm; // Data: (Met-meanGenMet)/sigmaGenMet.vs.Met for -sigma parametrization
    TProfile* fDataMetSigma_Met_dvp; // Data: (Met-meanGenMet)/sigmaGenMet.vs.Met for +sigma parametrization
    TProfile* fToyMetSigma_toyMet_def; // Toy MC: (toyMet-meanGenMet)/sigmaGenMet.vs.toyMet for default parametrization
    TProfile* fToyMetSigma_toyMet_dvm; // Toy MC: (toyMet-meanGenMet)/sigmaGenMet.vs.toyMet for -sigma parametrization
    TProfile* fToyMetSigma_toyMet_dvp; // Toy MC: (toyMet-meanGenMet)/sigmaGenMet.vs.toyMet for +sigma parametrization
  };


  struct PhoMetStudyHisto_t {  //------------ histograms for study of Pho-Met correlations  

    TH1F* fPhoMet_eta; // EtaDet of photons
    TH1F* fPhoMet_xces; // Xces of photons
    TH1F* fPhoMet_zces; // Zces of photons
    TH1F* fPhoMet_signeta; // sign(Zvx)*EtaDet of photons
    TH1F* fPhoMet_signzces; // sign(Zvx)*Zces of photons

    TH1F* fPhoMet_MetSig5dPhi02_signeta; // sign(Zvx)*EtaDet of photons if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_MetSig5dPhi02_signzces; // sign(Zvx)*Zces of photons if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_Fr_MetSig5dPhi02_signeta; // Fraction of photons as a function of sign(Zvx)*EtaDet if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_Fr_MetSig5dPhi02_signzces; // Fraction of photons as a function of sign(Zvx)*Zces if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5

    TH1F* fPhoMet_dPhi02_eta; // EtaDet of photons if dPhi=Phi(pho)-Phi(met)<0.2
    TH1F* fPhoMet_dPhi02_xces; // Xces of photons if dPhi=Phi(pho)-Phi(met)<0.2
    TH1F* fPhoMet_dPhi02_zces; // Zces of photons if dPhi=Phi(pho)-Phi(met)<0.2
    TH1F* fPhoMet_dPhi01_eta; // EtaDet of photons if dPhi=Phi(pho)-Phi(met)<0.1
    TH1F* fPhoMet_dPhi01_xces; // Xces of photons if dPhi=Phi(pho)-Phi(met)<0.1
    TH1F* fPhoMet_dPhi01_zces; // Zces of photons if dPhi=Phi(pho)-Phi(met)<0.1

    TH1F* fPhoMet_Fr_dPhi02_eta; // Fraction of photons as a function of EtaDet if dPhi=Phi(pho)-Phi(met)<0.2
    TH1F* fPhoMet_Fr_dPhi02_xces; // Fraction of photons as a function of Xces if dPhi=Phi(pho)-Phi(met)<0.2
    TH1F* fPhoMet_Fr_dPhi02_zces; // Fraction of photons as a function of Zces if dPhi=Phi(pho)-Phi(met)<0.2
    TH1F* fPhoMet_Fr_dPhi01_eta; // Fraction of photons as a function of EtaDet if dPhi=Phi(pho)-Phi(met)<0.1
    TH1F* fPhoMet_Fr_dPhi01_xces; // Fraction of photons as a function of Xces if dPhi=Phi(pho)-Phi(met)<0.1
    TH1F* fPhoMet_Fr_dPhi01_zces; // Fraction of photons as a function of Zces if dPhi=Phi(pho)-Phi(met)<0.1

    TH1F* fPhoMet_MetSig5dPhi02_eta; // EtaDet of photons if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_MetSig5dPhi02_xces; // Xces of photons if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_MetSig5dPhi02_zces; // Zces of photons if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_MetSig5dPhi01_eta; // EtaDet of photons if dPhi=Phi(pho)-Phi(met)<0.1 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_MetSig5dPhi01_xces; // Xces of photons if dPhi=Phi(pho)-Phi(met)<0.1 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_MetSig5dPhi01_zces; // Zces of photons if dPhi=Phi(pho)-Phi(met)<0.1 and Met/Sigma(Et_pho)>5

    TH1F* fPhoMet_Fr_MetSig5dPhi02_eta; // Fraction of photons as a function of EtaDet if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_Fr_MetSig5dPhi02_xces; // Fraction of photons as a function of Xces if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_Fr_MetSig5dPhi02_zces; // Fraction of photons as a function of Zces if dPhi=Phi(pho)-Phi(met)<0.2 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_Fr_MetSig5dPhi01_eta; // Fraction of photons as a function of EtaDet if dPhi=Phi(pho)-Phi(met)<0.1 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_Fr_MetSig5dPhi01_xces; // Fraction of photons as a function of Xces if dPhi=Phi(pho)-Phi(met)<0.1 and Met/Sigma(Et_pho)>5
    TH1F* fPhoMet_Fr_MetSig5dPhi01_zces; // Fraction of photons as a function of Zces if dPhi=Phi(pho)-Phi(met)<0.1 and Met/Sigma(Et_pho)>5


    //----------------------------------------- Correlation between electrons and met
    TH1F* fEleMet_eta; // EtaDet of electrons
    TH1F* fEleMet_phi; // PhiDet of electrons

    TH1F* fEleMet_dPhi02_eta; // EtaDet of electrons if dPhi=PhiDet(ele)-Phi(met)<0.2
    TH1F* fEleMet_dPhi02_phi; // PhiDet of electrons if dPhi=PhiDet(ele)-Phi(met)<0.2
    TH1F* fEleMet_dPhi01_eta; // EtaDet of electrons if dPhi=PhiDet(ele)-Phi(met)<0.1
    TH1F* fEleMet_dPhi01_phi; // PhiDet of electrons if dPhi=PhiDet(ele)-Phi(met)<0.1

    TH1F* fEleMet_Fr_dPhi02_eta; // Fraction of electrons as a function of EtaDet if dPhi=PhiDet(ele)-Phi(met)<0.2
    TH1F* fEleMet_Fr_dPhi02_phi; // Fraction of electrons as a function of Xces if dPhi=PhiDet(ele)-Phi(met)<0.2
    TH1F* fEleMet_Fr_dPhi01_eta; // Fraction of electrons as a function of EtaDet if dPhi=PhiDet(ele)-Phi(met)<0.1
    TH1F* fEleMet_Fr_dPhi01_phi; // Fraction of electrons as a function of Xces if dPhi=PhiDet(ele)-Phi(met)<0.1

    TH1F* fEleMet_MetSig5dPhi02_eta; // EtaDet of electrons if dPhi=PhiDet(ele)-Phi(met)<0.2 and Met/Sigma(Et_ele)>5
    TH1F* fEleMet_MetSig5dPhi02_phi; // PhiDet of electrons if dPhi=PhiDet(ele)-Phi(met)<0.2 and Met/Sigma(Et_ele)>5
    TH1F* fEleMet_MetSig5dPhi01_eta; // EtaDet of electrons if dPhi=PhiDet(ele)-Phi(met)<0.1 and Met/Sigma(Et_ele)>5
    TH1F* fEleMet_MetSig5dPhi01_phi; // PhiDet of electrons if dPhi=PhiDet(ele)-Phi(met)<0.1 and Met/Sigma(Et_ele)>5

    TH1F* fEleMet_Fr_MetSig5dPhi02_eta; // Fraction of electrons as a function of EtaDet if dPhi=PhiDet(ele)-Phi(met)<0.2 and Met/Sigma(Et_ele)>5
    TH1F* fEleMet_Fr_MetSig5dPhi02_phi; // Fraction of electrons as a function of PhiDet if dPhi=PhiDet(ele)-Phi(met)<0.2 and Met/Sigma(Et_ele)>5
    TH1F* fEleMet_Fr_MetSig5dPhi01_eta; // Fraction of electrons as a function of EtaDet if dPhi=PhiDet(ele)-Phi(met)<0.1 and Met/Sigma(Et_ele)>5
    TH1F* fEleMet_Fr_MetSig5dPhi01_phi; // Fraction of electrons as a function of PhiDet if dPhi=PhiDet(ele)-Phi(met)<0.1 and Met/Sigma(Et_ele)>5
 
  };

  enum { myNjetArray = 100 };
protected:

  char histo_outputstring[200];
  //____________________________________   pointers to the data blocks used,
					// header block is always available via
                                        // TStnModule::GetHeaderBlock()
  TStnVertexBlock*   fZVertexBlock;
  TStnMetBlock*      fMetBlock;
  TStnJetBlock*      fJetBlockClu04;

//----------------------- need this for jet-EM matching
  TCalDataBlock* fCalData;
  TStnPhotonBlock* fPhotonBlock;
  TStnElectronBlock* fElectronBlock;
  typedef std::vector<TCalTower*> CalDataArray; // need for beam halo
  typedef std::vector<TCalTower*>::iterator CalDataArrayI; // need for beam halo
  

  //______________________________________  ***** histograms filled
  JetGeneral_t       fHistJet04; // general histograms for JetClu 0.4

  //______________________________________ Met study histograms
  MetStudyHisto_t  fMetStudyJet04sc3; // met study histograms for JetClu 0.4 and Et>15 GeV

  //_____________________________________ histograms for study of Pho-Jet matching
  MatchStudyHisto_t fMatchPhoJet04; 

//_____________________________________ histograms for met cleanup studies
  MetCleanupHisto_t fCleanup;
  MetCleanupHisto_t fGenCleanup;
  MetCleanupHisto_t fGenDevmCleanup;
  MetCleanupHisto_t fGenDevpCleanup;
  MetCleanupHisto_t fGenDevmUnclCleanup;
  MetCleanupHisto_t fGenDevpUnclCleanup;
//___________________________________ histograms for MetProbability
  MetProbHisto_t fMetProbability;
  
  //_________________________________ histograms to study correlation between pho and met
  PhoMetStudyHisto_t fPhoMetStudy;

  //______________________________________ ***** jet params for different jet cones
  JetStuff jet04stuff;
  CommonStuff allstuff;
  std::vector<MatchStuff> matchstuff;
  //_____________________________________ ***** met results
  MetResults met_results;

//--------------------------------------------------------------------------------------------
//_________________________________________ ***** service parameters

  int fUseVerbose;   // status code to invoke Debugging mode (print out)
  int fDumpEvent;   // status code to invoke Debugging mode (print out)
  
  int myJetRun;         // to be used in corrections and for crosschecks
  int myNvx_class12;     // Nvx(class12) to be used in corrections
  double zvx_best;
  //_____________________________________________________________________________________________________
  int fJetAlgo;      // identifier of jet algorithm: 0-JetClu(0.4), 1-JetClu(0.7), 2-JetClu(1.0)
                       //                              3-MidPoint(0.4), 4-MidPoint(0.7), 5-MidPoint(1.0)
                       //                              6-KtCl(0.4), 7-KtCl(0.7), 8-KtCl(1.0)
  //_____________________________________________________________________________________________________
  int fMyMetScenario; // 0--raw Met; corrected for jets with:
                      // 1-- et>5, 2-- et>10, 3-- et>15, 4-- et>20
                        
  //_____________________________________________________________________________________________________
  int fJTC_coneSize; // parameters for Jet Corrections

  int fJTC_version;  // version of jet corrections: 0--Monte Carlo
                       //                             1..5(?)--Data
                       // for 5.3.1 it's recommended to use version 5  

  int fJTC_level;    // Level of corrections: 
                       //       Level=1 Relative Energy Corrections 
                       //       Level=2 Relative Energy Corrections and Time Dependence 
                       //       Level=3 Raw Energy Scale corrections, Relative Energy 
                       //               Corrections and Time Dependence 
                       //       Level=4 Multiple Interaction corrections, and all of above 
                       //       Level=5 Absolute Energy Scale and all of above 
                       //       Level=6 Underlying Event corrections and all of above. 
                       //       Level=7 Out-of-cone and all of above.

  int fJTC_systcode; // code for systematics in corrections:
                       //     0--default corrections                       
                       //     1--Relative Correction
                       //     2--Central Cal Stability
                       //     3--Scale Correction
                       //     4--Multiple Interaction Correction
                       //     5--Absolute Correction
                       //     6--Underlying Event Correction
                       //     7--Out-Of-Cone Correction
                       //     8--Splash Out Correction
                       //     100--Total Uncertainty
                       //--- if systcode>0 -- plus sigma deviation(?)
                       //--- if systcode<0 -- minos sigma deviation(?)
                       //-------- All these codes have to be double-checked!!!!
  int fJTC_imode;      // correction mode in V5* data and MC: 0=MC, 1=data


//--------------------------------------------------------------------------------
//_______________ ***** met model parametrization.
//                for MetX & MetY: sigma=p0+p1*sqrt(SumEt)
//                for MetX & MetY: mean=p0+p1*SumEt
//                first index is for di-pho sideband ([0]) & central Z->ee ([1]) parametrizations   
  double metmodel_sigmaX[2][2];    // sigmaX=p[0]+p[1]*sqrt(corrSumEt)
  double metmodel_sigmaX_er[2][2];
  double metmodel_sigmaY[2][2];    // sigmaY=p[0]+p[1]*sqrt(corrSumEt) 
  double metmodel_sigmaY_er[2][2];
  double metmodel_meanX[2][2];     // meanX=p[0]+p[1]*corrSumEt
  double metmodel_meanX_er[2][2];
  double metmodel_meanY[2][2];     // meanY=p[0]+p[1]*corrSumEt
  double metmodel_meanY_er[2][2];

  //-------------- For unclustered MetModel with Njet(cut)=0
  double metmodel_sigmaXcorr[2]; // correlation coefficient for sigmaX
  double metmodel_sigmaYcorr[2]; // correlation coefficient for sigmaY
  double metmodel_meanXcorr[2]; // correlation coefficient for meanX
  double metmodel_meanYcorr[2]; // correlation coefficient for meanY

  int fUnclParamSwitch; // switch to set sideband or Z parametrization as default
                        // 0=sideband; 1=Z

  int Npoints; // number of points to be generated per value of corrSumEt  
  int fUseMetPDF; // scenario how to calculate MetPDF
  std::vector<TVector2> MetGen_def; // global vectors to be used in MetPDF
  double MetToy_max; // max value of toy Met
  double MetToy_min; // min value of toy Met
  std::vector<TVector2> MetGen_dvm;
  std::vector<TVector2> MetGen_dvp;
  int fSelectSigMetEvent; // to select events with Significant MET according to MyMetPDF
  int SigMetEvent_status; // global status code for Significant MET according to MyMetPDF

  int bad_EMjet_match_flag; // global flag to mark events where not all EM objects are matched to jets
//--------------------------------------------------------------------------------
//______________________________________ ***** cuts

  //________________________
  double MinEtThr; // Et(lev6) threshold on jets to be considered in the analysis, after eta-cut & removing EM objects
  double MaxEtThr;
  int MinNjetThr; // number of jets above the threshold
  int MaxNjetThr; // number of jets above the threshold
  double MinJetEta; // cut of the eta-range for jets to be counted and used in MET corrections
  double MaxJetEta;
  double MinEt1st; // cut on the leading jet Et(lev6) after removing EM objects, must also pass eta-cut
  double MaxEt1st;
  double MinEt2nd; // cut on the 2nd jet Et(lev6) after removing EM objects, must also pass eta-cut
  double MaxEt2nd;
  double MinDeltaRgjMin; // cut on min dR between photons and jets after removing EM objects, must also pass eta-cut
  double MaxDeltaRgjMin; 
  double MinDeltaEtagjMin; // cut on min dEta between photons and jets after removing EM objects, must also pass eta-cut
  double MaxDeltaEtagjMin; 
  double MinDeltaPhigjMin; // cut on min dPhi between photons and jets after removing EM objects, must also pass eta-cut
  double MaxDeltaPhigjMin; 
  double MinDeltaRgj; // cut on dR between any photons and jets after removing EM objects, must also pass eta-cut
  double MaxDeltaRgj; 
  double MinDeltaEtagj; // cut on dEta between any photons and jets after removing EM objects, must also pass eta-cut
  double MaxDeltaEtagj; 
  double MinDeltaPhigj; // cut on dPhi between any photons and jets after removing EM objects, must also pass eta-cut
  double MaxDeltaPhigj; 
  double MinDeltaPhiJMet; // cut on min dPhi between raw jet and raw Met
  double MaxDeltaPhiJMet; // cut on max dPhi between raw jet and raw Met
  double MinMet2Et; // cut on MET/Et(raw jet)
  double MaxMet2Et;
  double MinMet2Etlev6; // cut on MET/Et(lev6 jet)
  double MaxMet2Etlev6;
  double MindZ; // cut on dZ=Zvx-Zjet
  double MaxdZ;
  double MinMjj; // cut on Mjj of two leading jets
  double MaxMjj;
  double MinMj1g1; // cut on M(jet1-pho1)
  double MaxMj1g1;
  double MinMj1g2; // cut on M(jet1-pho2)
  double MaxMj1g2;
  double MinMj2g1; // cut on M(jet2-pho1)
  double MaxMj2g1;
  double MinMj2g2; // cut on M(jet2-pho2)
  double MaxMj2g2;
  double MinMjg; // cut on Mass of any pho-jet system
  double MaxMjg; 
  double MinQtJet; // cut on Qt of jets
  double MaxQtJet;
  double MinHt; // cut on Ht of jets, photons, Met, leptons
  double MaxHt;
  int MinNjet; // cut on the number of jets(clusters) before removing EM objects
  int MaxNjet;
  int MinNjet5; // cut on the number of jets with Et(lev6)>5 after removing EM objects
  int MaxNjet5;
  int MinNjet10; // cut on the number of jets with Et(lev6)>10 after removing EM objects
  int MaxNjet10;
  int MinNjet15; // cut on the number of jets with Et(lev6)>15 after removing EM objects
  int MaxNjet15;
  int MinNjet20; // cut on the number of jets with Et(lev6)>20 after removing EM objects
  int MaxNjet20;
  int MinNjet25; // cut on the number of jets with Et(lev6)>25 after removing EM objects
  int MaxNjet25;

  int fSelectMetEvent; // switch to select large met events
  int fRemoveDuplicate; // switch to remove duplicate objects
  int fRemoveBadMet; // switch to remove events with bad met (met along the jet which is close to a crack)  
                              
  char* fDatFileName;
  char* fLargeMetEventFileName;
  char* fDumpEventFileName;

public:
  TMyJetFilterModule(const char* name ="MyJetFilter", 
		     const char* title="MyJetFilter");
  ~TMyJetFilterModule();
  //_______________________________________________________ ****** accessors

  JetGeneral_t*       GetHistJet04()          { return &fHistJet04; }
  TStnJetBlock*       GetJetBlockClu04() { return fJetBlockClu04; }

  //_______________________________________________________ ****** cut values
  int    GetUseVerbose()           { return fUseVerbose; }
  int    GetDumpEvent()           { return fDumpEvent; }
  int    GetmyJetRun()           { return myJetRun; } 
  int    GetJetAlgo()      	   { return fJetAlgo; }      	    
  int    GetMyMetScenario() { return fMyMetScenario; }
  int    GetUseMetPDFscenario() { return fUseMetPDF; }
  int    GetSelectSigMetEvent() { return fSelectSigMetEvent; }

  double GetMinEtThr() { return MinEtThr; }
  double GetMaxEtThr() { return MaxEtThr; }
  int GetMinNjetThr() { return MinNjetThr; }
  int GetMaxNjetThr() { return MaxNjetThr; }
  double GetMinJetEta() { return MinJetEta; }
  double GetMaxJetEta() { return MaxJetEta; }
  double GetMinEt1st() { return MinEt1st; }
  double GetMaxEt1st() { return MaxEt1st; }
  double GetMinEt2nd() { return MinEt2nd; }
  double GetMaxEt2nd() { return MaxEt2nd; }
  double GetMinDeltaRgjMin() { return MinDeltaRgjMin; }
  double GetMaxDeltaRgjMin() { return MaxDeltaRgjMin; } 
  double GetMinDeltaEtagjMin() { return MinDeltaEtagjMin; }
  double GetMaxDeltaEtagjMin() { return MaxDeltaEtagjMin; } 
  double GetMinDeltaPhigjMin() { return MinDeltaPhigjMin; }
  double GetMaxDeltaPhigjMin() { return MaxDeltaPhigjMin; } 
  double GetMinDeltaRgj() { return MinDeltaRgj; }
  double GetMaxDeltaRgj() { return MaxDeltaRgj; } 
  double GetMinDeltaEtagj() { return MinDeltaEtagj; }
  double GetMaxDeltaEtagj() { return MaxDeltaEtagj; } 
  double GetMinDeltaPhigj() { return MinDeltaPhigj; }
  double GetMaxDeltaPhigj() { return MaxDeltaPhigj; } 
  double GetMinDeltaPhiJMet() { return MinDeltaPhiJMet; }
  double GetMaxDeltaPhiJMet() { return MaxDeltaPhiJMet; }
  double GetMinMet2Et() { return MinMet2Et; }
  double GetMaxMet2Et() { return MaxMet2Et; }
  double GetMinMet2Etlev6() { return MinMet2Etlev6; }
  double GetMaxMet2Etlev6() { return MaxMet2Etlev6; }
  double GetMindZ() { return MindZ; }
  double GetMaxdZ() { return MaxdZ; }
  double GetMinMjj() { return MinMjj; }
  double GetMaxMjj() { return MaxMjj; }
  double GetMinMj1g1() { return MinMj1g1; }
  double GetMaxMj1g1() { return MaxMj1g1; }
  double GetMinMj1g2() { return MinMj1g2; }
  double GetMaxMj1g2() { return MaxMj1g2; }
  double GetMinMj2g1() { return MinMj2g1; }
  double GetMaxMj2g1() { return MaxMj2g1; }
  double GetMinMj2g2() { return MinMj2g2; }
  double GetMaxMj2g2() { return MaxMj2g2; }
  double GetMinMjg() { return MinMjg; }
  double GetMaxMjg() { return MaxMjg; } 
  double GetMinQtJet() { return MinQtJet; }
  double GetMaxQtJet() { return MaxQtJet; }
  double GetMinHt() { return MinHt; }
  double GetMaxHt() { return MaxHt; }
  int GetMinNjet() { return MinNjet; }
  int GetMaxNjet() { return MaxNjet; }
  int GetMinNjet5() { return MinNjet5; }
  int GetMaxNjet5() { return MaxNjet5; }
  int GetMinNjet10() { return MinNjet10; }
  int GetMaxNjet10() { return MaxNjet10; }
  int GetMinNjet15() { return MinNjet15; }
  int GetMaxNjet15() { return MaxNjet15; }
  int GetMinNjet20() { return MinNjet20; }
  int GetMaxNjet20() { return MaxNjet20; }
  int GetMinNjet25() { return MinNjet25; }
  int GetMaxNjet25() { return MaxNjet25; }
  double GetMetToy_max() { return MetToy_max; }
  double GetMetToy_min() { return MetToy_min; }
  int GetSelectMetEvent() { return fSelectMetEvent; }
  int GetRemoveDuplicate() { return fRemoveDuplicate; }
  int GetRemoveBadMet() { return fRemoveBadMet; }

  //__________________________________________________________ accessors to jet 
  //__________________________________________________________ correction params
  int    GetJTC_coneSize()       { return fJTC_coneSize; } 
  int    GetJTC_version()  	 { return fJTC_version; }
  int    GetJTC_level()    	 { return fJTC_level; }  
  int    GetJTC_systcode()	 { return fJTC_systcode; }
  int    GetJTC_imode()	         { return fJTC_imode; }
  
  char* GetDatFileName()       { return fDatFileName; }
  char* GetLargeMetEventFileName() { return fLargeMetEventFileName; }
  char* GetDumpEventFileName() { return fDumpEventFileName; }
  //__________________________________________________________ accessors to the
  //__________________________________________________________ output params

  //_________ need to specify which jets to return (variable "cone"): 0--cone 0.4; 1--cone 0.7; 2--cone 1.0
  TLorentzVector* GetMyJet_raw(int cone, int i); 
  TLorentzVector* GetMyJet_lev1(int cone, int i); 
  TLorentzVector* GetMyJet_lev4(int cone, int i); 
  TLorentzVector* GetMyJet_lev5(int cone, int i); 
  TLorentzVector* GetMyJet_lev6(int cone, int i); 
  TLorentzVector* GetMyJetNoEMobj_raw(int cone, int i); 
  TLorentzVector* GetMyJetNoEMobj_lev1(int cone, int i); 
  TLorentzVector* GetMyJetNoEMobj_lev4(int cone, int i); 
  TLorentzVector* GetMyJetNoEMobj_lev5(int cone, int i); 
  TLorentzVector* GetMyJetNoEMobj_lev6(int cone, int i); 
  int GetMyJetNoEMobj_lev7size() { return (jet04stuff.Jet_lev7_noEMobj).size();}
  TLorentzVector* GetMyJetNoEMobj_lev7(int cone, int i); 
  double GetMyJetEtaDet(int cone, int i); // before removing EM object 
  double GetMyJetEtaDetCorr(int cone, int i); // after removing EM object
  double GetMyJetEmFrRaw(int cone, int i); // before removing EM object
  double GetMyJetEmFrCorr(int cone, int i); // before removing EM object
  int GetMyJetNtrk(int cone, int i);
  int GetMyJetNtwr(int cone, int i);
  int GetMyJetNobjMatch(int cone, int i);
  int GetMyJetNphoMatch(int cone, int i);
  int GetMyJetNeleMatch(int cone, int i);
  int GetMyJetNmuMatch(int cone, int i);
  int GetMyJetNtauMatch(int cone, int i);
  int GetMyJetNbtagMatch(int cone, int i);
  int GetMyJetBlockInd(int cone, int i); // original jet index in JetBlock, need this after jet reordering
  int GetMyNjet(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25
  double GetMySumEtCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  double GetMyHtCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  double GetMyMetCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  double GetMyMetXCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  double GetMyMetYCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  double GetMyMetPhiCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20

  //________________________________________________________ ****** setters
  void SetUseVerbose(int cut) { fUseVerbose = cut; }
  void SetDumpEvent(int cut) { fDumpEvent = cut; }
  void SetmyJetRun(int run) { myJetRun = run; } 
  void SetJetAlgo(int cut) { fJetAlgo = cut; }      	    
  void SetMyMetScenario(int cut) { fMyMetScenario = cut; }
  void SetUseMetPDFscenario(int cut) { fUseMetPDF = cut; }
  void SetSelectSigMetEvent(int cut) { fSelectSigMetEvent = cut; }


  void SetJTC_coneSize    (int param)           { fJTC_coneSize = param; } 
  void SetJTC_version     (int param)           { fJTC_version = param; }
  void SetJTC_level       (int param)           { fJTC_level = param; }  
  void SetJTC_systcode    (int param)	        { fJTC_systcode = param; }
  void SetJTC_imode       (int param)	        { fJTC_imode = param; }

  //------------------------ Setting Met Model parameters

  void SetMetModelSigmaCorr(int Axis, int modelInd, double par)
  {
    if(Axis==0) metmodel_sigmaXcorr[modelInd]=par; 
    if(Axis==1) metmodel_sigmaYcorr[modelInd]=par; 
    if(Axis>1 || Axis<0) std::cout<<"Warning in SetMetModelSigmaCorr: wrong axis code!!! Axis: 0-x or 1-y"<<std::endl;
    return;
  }
  void SetMetModelMeanCorr(int Axis, int modelInd, double par)
  {
    if(Axis==0) metmodel_meanXcorr[modelInd]=par; 
    if(Axis==1) metmodel_meanYcorr[modelInd]=par; 
    if(Axis>1 || Axis<0) std::cout<<"Warning in SetMetModelMeanCorr: wrong axis code!!! Axis: 0-x or 1-y"<<std::endl;
    return;
  }

  void SetMetModelSigma(int modelInd, int Axis, double par0, double par1, double er_par0, double er_par1)
  {
    // modelInd: 0=noJet parametrization; 1=withJet parametrization 
    if(Axis==0 && modelInd>-1 && modelInd<2)
      {
	metmodel_sigmaX[modelInd][0]=par0;
	metmodel_sigmaX[modelInd][1]=par1;
	metmodel_sigmaX_er[modelInd][0]=er_par0;
	metmodel_sigmaX_er[modelInd][1]=er_par1;
      }
    if(Axis==1 && modelInd>-1 && modelInd<2)
      {
	metmodel_sigmaY[modelInd][0]=par0;
	metmodel_sigmaY[modelInd][1]=par1;
	metmodel_sigmaY_er[modelInd][0]=er_par0;
	metmodel_sigmaY_er[modelInd][1]=er_par1;
      }
    if(Axis>1 || Axis<0) std::cout<<"Warning in SetMetModelSigma: wrong axis code!!! Axis: 0-x or 1-y"<<std::endl;
    if(modelInd>1 || modelInd<0) 
      std::cout<<"Warning in SetMetModelSigma: wrong Met Model code!!! modelInd: 0-sideband or 1-Zee"<<std::endl;
    return;
  }
  void SetMetModelMean(int modelInd, int Axis, double par0, double par1, double er_par0, double er_par1)
  {
    if(Axis==0 && modelInd>-1 && modelInd<2)
      {
	metmodel_meanX[modelInd][0]=par0;
	metmodel_meanX[modelInd][1]=par1;
	metmodel_meanX_er[modelInd][0]=er_par0;
	metmodel_meanX_er[modelInd][1]=er_par1;
      }
    if(Axis==1 && modelInd>-1 && modelInd<2)
      {
	metmodel_meanY[modelInd][0]=par0;
	metmodel_meanY[modelInd][1]=par1;
	metmodel_meanY_er[modelInd][0]=er_par0;
	metmodel_meanY_er[modelInd][1]=er_par1;
      }
    if(Axis>1 || Axis<0) std::cout<<"Warning in SetMetModelMean: wrong axis code!!! Axis: 0-x or 1-y"<<std::endl;
    if(modelInd>1 || modelInd<0) 
      std::cout<<"Warning in SetMetModelMean: wrong Met Model code!!! modelInd: 0-sideband or 1-Zee"<<std::endl;
    return;
  }
  void SetNpointsToGenerate(int n) { Npoints=n; }
  void SetUnclParamSwitch(int n) { if(n>-1 && n<2) fUnclParamSwitch=n; }

  //_________________________________________________ end of setting met model params

  void  SetDatFileName(char* fname)       { fDatFileName = fname; }    
  void  SetLargeMetEventFileName( char* fname) { fLargeMetEventFileName = fname; }
  void  SetDumpEventFileName( char* fname) { fDumpEventFileName = fname; }

  //_________________________ setters of cuts
  void SetMinEtThr(double param) { MinEtThr=param; }
  void SetMaxEtThr(double param) { MaxEtThr=param; }
  void SetMinNjetThr(int param) { MinNjetThr=param; }
  void SetMaxNjetThr(int param) { MaxNjetThr=param; }
  void SetMinJetEta(double param) { MinJetEta=param; }
  void SetMaxJetEta(double param) { MaxJetEta=param; }
  void SetMinEt1st(double param) { MinEt1st=param; }
  void SetMaxEt1st(double param) { MaxEt1st=param; }
  void SetMinEt2nd(double param) { MinEt2nd=param; }
  void SetMaxEt2nd(double param) { MaxEt2nd=param; }
  void SetMinDeltaRgjMin(double param) { MinDeltaRgjMin=param; }
  void SetMaxDeltaRgjMin(double param) { MaxDeltaRgjMin=param; } 
  void SetMinDeltaEtagjMin(double param) { MinDeltaEtagjMin=param; }
  void SetMaxDeltaEtagjMin(double param) { MaxDeltaEtagjMin=param; } 
  void SetMinDeltaPhigjMin(double param) { MinDeltaPhigjMin=param; }
  void SetMaxDeltaPhigjMin(double param) { MaxDeltaPhigjMin=param; } 
  void SetMinDeltaRgj(double param) { MinDeltaRgj=param; }
  void SetMaxDeltaRgj(double param) { MaxDeltaRgj=param; } 
  void SetMinDeltaEtagj(double param) { MinDeltaEtagj=param; }
  void SetMaxDeltaEtagj(double param) { MaxDeltaEtagj=param; } 
  void SetMinDeltaPhigj(double param) { MinDeltaPhigj=param; }
  void SetMaxDeltaPhigj(double param) { MaxDeltaPhigj=param; } 
  void SetMinDeltaPhiJMet(double param) { MinDeltaPhiJMet=param; }
  void SetMaxDeltaPhiJMet(double param) { MaxDeltaPhiJMet=param; }
  void SetMinMet2Et(double param) { MinMet2Et=param; }
  void SetMaxMet2Et(double param) { MaxMet2Et=param; }
  void SetMinMet2Etlev6(double param) { MinMet2Etlev6=param; }
  void SetMaxMet2Etlev6(double param) { MaxMet2Etlev6=param; }
  void SetMindZ(double param) { MindZ=param; }
  void SetMaxdZ(double param) { MaxdZ=param; }
  void SetMinMjj(double param) { MinMjj=param; }
  void SetMaxMjj(double param) { MaxMjj=param; }
  void SetMinMj1g1(double param) { MinMj1g1=param; }
  void SetMaxMj1g1(double param) { MaxMj1g1=param; }
  void SetMinMj1g2(double param) { MinMj1g2=param; }
  void SetMaxMj1g2(double param) { MaxMj1g2=param; }
  void SetMinMj2g1(double param) { MinMj2g1=param; }
  void SetMaxMj2g1(double param) { MaxMj2g1=param; }
  void SetMinMj2g2(double param) { MinMj2g2=param; }
  void SetMaxMj2g2(double param) { MaxMj2g2=param; }
  void SetMinMjg(double param) { MinMjg=param; }
  void SetMaxMjg(double param) { MaxMjg=param; } 
  void SetMinQtJet(double param) { MinQtJet=param; }
  void SetMaxQtJet(double param) { MaxQtJet=param; }
  void SetMinHt(double param) { MinHt=param; }
  void SetMaxHt(double param) { MaxHt=param; }
  void SetMinNjet(int param) { MinNjet=param; }
  void SetMaxNjet(int param) { MaxNjet=param; }
  void SetMinNjet5(int param) { MinNjet5=param; }
  void SetMaxNjet5(int param) { MaxNjet5=param; }
  void SetMinNjet10(int param) { MinNjet10=param; }
  void SetMaxNjet10(int param) { MaxNjet10=param; }
  void SetMinNjet15(int param) { MinNjet15=param; }
  void SetMaxNjet15(int param) { MaxNjet15=param; }
  void SetMinNjet20(int param) { MinNjet20=param; }
  void SetMaxNjet20(int param) { MaxNjet20=param; }
  void SetMinNjet25(int param) { MinNjet25=param; }
  void SetMaxNjet25(int param) { MaxNjet25=param; }
  void SetMetToy_max(double param) { MetToy_max=param; }
  void SetMetToy_min(double param) { MetToy_min=param; }
  void SetSelectMetEvent(int param) { fSelectMetEvent=param; }
  void SetRemoveDuplicate(int param) { fRemoveDuplicate=param; }
  void SetRemoveBadMet(int param) { fRemoveBadMet=param; }

  //__________________________________________________________ setters of the
  //__________________________________________________________ output params

  //_________ need to specify which jets to return (variable "cone"): 0--cone 0.4; 1--cone 0.7; 2--cone 1.0
  void SetMyJet_raw(int cone, int i, TLorentzVector vec); 
  void SetMyJet_lev1(int cone, int i, TLorentzVector vec); 
  void SetMyJet_lev4(int cone, int i, TLorentzVector vec); 
  void SetMyJet_lev5(int cone, int i, TLorentzVector vec); 
  void SetMyJet_lev6(int cone, int i, TLorentzVector vec); 
  void SetMyJetNoEMobj_raw(int cone, int i, TLorentzVector vec); 
  void SetMyJetNoEMobj_lev1(int cone, int i, TLorentzVector vec); 
  void SetMyJetNoEMobj_lev4(int cone, int i, TLorentzVector vec); 
  void SetMyJetNoEMobj_lev5(int cone, int i, TLorentzVector vec); 
  void SetMyJetNoEMobj_lev6(int cone, int i, TLorentzVector vec); 
  void SetMyJetEtaDet(int cone, int i, double param); // before removing EM object 
  void SetMyJetEtaDetCorr(int cone, int i, double param); // after removing EM object
  void SetMyJetEmFrRaw(int cone, int i, double param); // before removing EM object
  void SetMyJetEmFrCorr(int cone, int i, double param); // before removing EM object
  void SetMyJetNtrk(int cone, int i, int param);
  void SetMyJetNtwr(int cone, int i, int param);
  void SetMyJetNobjMatch(int cone, int i, int param);
  void SetMyJetNphoMatch(int cone, int i, int param);
  void SetMyJetNeleMatch(int cone, int i, int param);
  void SetMyJetNmuMatch(int cone, int i, int param);
  void SetMyJetNtauMatch(int cone, int i, int param);
  void SetMyJetNbtagMatch(int cone, int i, int param);
  void SetMyJetBlockInd(int cone, int i, int param); // original jet index in JetBlock, need this after jet reordering
  void SetMyNjet(int cone, int threshold, int param); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25
  void SetMySumEtCorr(int cone, int threshold, double param); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  void SetMyMetCorr(int cone, int threshold, double px, double py); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  void SetMyMetCorr(int cone, int threshold, const TVector2& vec); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20

  //____________________________________________ Service methods


  void ClearModuleOutput(); // clears the Module output  		 
  void ClearJetStuff(JetStuff &stuff);     // clears the JetStuff structures 		 
  void ClearCommonStuff(CommonStuff &stuff); // clears CommonStuff structures 
  void ClearMatchStuff(); // clears information about pho-jet matching

  char*   HistoSummary(TH1F* myhisto)
  {
    //    char* outputstring="";
    sprintf(histo_outputstring,"...%s   %f   %f",
	    myhisto->GetName(),
	    myhisto->GetMean(),
	    myhisto->GetRMS());
    return histo_outputstring;
  }
  
  void    Display();
 //___________________________________________ ****** overloaded methods of 
 //___________________________________________        TStnModule
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
 //___________________________________________ ****** other methods
  void BookHistograms();
  void BookJetHistograms(JetGeneral_t& Hist, const char* Folder, const char* algoname); // booking general jet histo
  void BookMetStudyHistograms(MetStudyHisto_t& Hist, const char* Folder, const char* algoname); // booking met histo
  void BookMatchStudyHistograms(MatchStudyHisto_t& Hist, const char* Folder, const char* algoname); // booking match histo
  void BookMetCleanupHistograms(MetCleanupHisto_t& Hist, const char* Folder, const char* algoname); // met cleanup studies
  void BookMetProbHistograms(MetProbHisto_t& Hist, const char* Folder, const char* algoname); // MetProbability histograms
  void BookPhoMetStudyHistograms(PhoMetStudyHisto_t& Hist, const char* Folder); // booking Pho-Met study histograms

  void FillJetHistogramsB(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12); // filling general histo for particular Jet Cone
  void FillJetHistogramsA(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12); // filling general histo for particular Jet Cone
  void FillMetStudyHistograms(MetStudyHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, int metcode, int nvx);// filling met histo
  void FillMatchingHistograms(MatchStudyHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, MatchStuff match);// filling match histo
  void FillMetCleanupHistograms(MetCleanupHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, int metcode); // met cleanup studies
  void FillGenMetCleanupHistograms(MetCleanupHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, 
				   TVector2 genMet, std::vector<TLorentzVector> vec); // cleanup studies using generated met
  void FillPhoMetStudyHistograms(PhoMetStudyHisto_t& Hist, TVector2 Met, CommonStuff miscstuff, double zvx);  // filling histograms for pho-met study
  void DoFinalPhoMetHisto(PhoMetStudyHisto_t& Hist);  // finalize histograms for pho-met study

  void DoJetResHisto1(JetGeneral_t& Hist); // finalize JetResolution histograms
  int MyDuplicateCut(JetStuff jetstuff, CommonStuff miscstuff);
  int MyNjetCut(JetStuff jetstuff);
  int MyGenericCut(JetStuff jetstuff, double dz);
  int MyAngularCut(JetStuff jetstuff, CommonStuff miscstuff);
  int MyMassCut(JetStuff jetstuff, CommonStuff miscstuff);
  int MyAnalysisCut(JetStuff jetstuff, CommonStuff miscstuff, int metcode); // metcode is Met scenario
//   int MyMetCleanUpCut(JetStuff jetstuff, CommonStuff miscstuff); 
  int MyMetCleanUpCut(JetStuff jetstuff, CommonStuff miscstuff, std::vector<TLorentzVector> jetvec, TVector2 metvec); // Test version as of 04/04/07
  int MyLargeMetEvent(JetStuff jetstuff, CommonStuff miscstuff, int metscenario, int Nvx12, int runN, int eventN); // selects my large MET events.
  void MyDumpEvent(JetStuff jetstuff, CommonStuff miscstuff, int metscenario, int Nvx12, int runN, int eventN); // event dump.

  void DoCommonStuff(CommonStuff &miscstuff); // reads raw Met, pho, ele info and fills CommonStuff  
  void DoMyMet(CommonStuff miscstuff, JetStuff &jetstuff); // corrects Met & SumEt for jets 
  void DoMyJet(TStnJetBlock* fJetBlock, CommonStuff miscstuff, JetStuff &jetstuff, MatchStudyHisto_t& Hist); // does my jets
  void DoMyJetNoMatch(TStnJetBlock* fJetBlock, JetStuff &jetstuff); // does my jets, but doesn't remove EM objetcs
  void DoMyJetWithMatch(TStnJetBlock* fJetBlock, CommonStuff miscstuff, 
			JetStuff &jetstuff, MatchStudyHisto_t& Hist); // does my jets after removing EM objects
  void ReorderMyJets(JetStuff &jetstuff); // reorders jets after removing EM objects

  double CorrectEmFr(double emfr_old, double emfr_obj, double e_old, double e_obj); // re-calculates jet EmFr after removing EM object
  double CorrectEtaDet(TLorentzVector *vec_pho, TLorentzVector *vec_old, 
		       double jetetadet_old, double pho_etadet); // re-calculates jet EtaDet after removing EM object

  //--the following three functions could have been replaced by one which works with generic data type (to be done)
  void myExchange_tlv(TLorentzVector& val1, TLorentzVector& val2); //exchanges val1 and val2
  void myExchange_dbl(double& val1, double& val2); //exchanges val1 and val2
  void myExchange_int(int& val1, int& val2); //exchanges val1 and val2

  //------------------ need this for jet-EM matching
  void MatchCalorTowers(int jet_ind,TStnJetBlock* fJetBlock,TCalDataBlock *fCalDataBlock,CalDataArray* cdholder); // for jets
  void MatchCalorTowers(int pho_ind,TStnPhotonBlock* fPhotonBlock,TCalDataBlock *fCalDataBlock,CalDataArray* cdholder); // for photons
  void MatchCalorTowers(int ele_ind,TStnElectronBlock* fElectronBlock,TCalDataBlock *fCalDataBlock,CalDataArray* cdholder); // for electrons

  double GetMyJetResolution(double jetPtLev6, double deteta, int jetcone, int statcode); /* this
function returns my jet energy resolution */ 
  double MyJetRes_def(double jetPtLev6, int jetcone); // Jet Resolution: default function
  double MyJetRes_step(double jetPtLev6, int jetcone); // Jet Resolution: step function
//   double MyJetRes_syst1(double jetPtLev6, int jetcone); // Jet Resolution: systematics-1 function
  double MyJetRes_syst2(double jetPtLev6, int jetcone); // Jet Resolution: systematics-2 function
  double MyJetRes_syst3(double jetPtLev6, int jetcone); // Jet Resolution: systematics-3 function
  double MyJetRes_syst4(double jetPtLev6, int jetcone); // Jet Resolution: systematics-4 function
  double MyJetRes_stat(double jetPtLev6, int jetcone); // Jet Resolution: stat. uncertainty function
  double MyJetRes_total(double jetPtLev6, int jetcone); // Jet Resolution: total uncertainty function
  double MyJetRes_EtaCor(double deteta, int jetcone); // Jet Resolution: correction for eta dependence of resolution
  double MyJetRes_Had2Det(double jetPtLev6, int jetcone); // Jet Resolution: scale factor Had/Det

  void GenerateMyTotalMet(JetStuff jetstuff, CommonStuff miscstuff, MetCleanupHisto_t& Hist, 
			  int jetcone, int systcode, int systcodeUncl, int metcode, TVector2 &myGenMet); 
  void GenerateMyEMobjMet(CommonStuff miscstuff, 
			  int systcode, TVector2 &myEMGenMet); // new function (09/20/06), MET due to EM obj. resolution
  void GetMyUnclMetResolution(double sumEt, int systcode, 
			      double &sigmaMx, double &sigmaMy, 
			      double &meanMx, double &meanMy); // generates Sigma and Mean for "unclustered" Met
  void MyNewJetWithMet(JetStuff jetstuff, CommonStuff miscstuff, int jetcone, int jes_systcode, std::vector<TLorentzVector>& newJet); // adding jet and Met
  //______________ counts generated events above Met cut
  void MyMetModelEventCount(JetStuff jetstuff, TVector2 myGenMet_def, TVector2 myGenMet_dvm, 
			    TVector2 myGenMet_dvp, TVector2 myGenMet_dvmUn, 
			    TVector2 myGenMet_dvpUn, MetResults &metstuff);
  //______________ counts data events above Met cut
  void MyMetEventCount(JetStuff jetstuff, CommonStuff miscstuff, int metcode, MetResults &metstuff);
  //______________ prints out a comparison of Met in Data to Met Model predictions
  void MyMetResults(MetResults &metstuff);
  //______________ calculates probability of Met according to MetModel for clustered and unclustered energy
  void MyMetPDF(int counter, TVector2 MetEvnt, TVector2 GenMet_def, 
		TVector2 GenMet_dvm, TVector2 GenMet_dvp, MetProbHisto_t& Hist, int &event_code);

  void CleanMetResults(MetResults &metstuff); // cleans met_results vectors
  
//   double GetCorrection(CorInit settings, TLorentzVector vec, double& emf, double etad); /* this
// function returns correction based on parameters localy set in the Module */ 
//   double GetCorrection(TLorentzVector vec, double& emf, double etad); /* this function returns 
// correction based on globaly initialized parameters */ 
  double GetCorrection(CorInit settings, TLorentzVector vec, float& emf, float etad); /* this
function returns correction based on parameters localy set in the Module */ 
  double GetCorrection(TLorentzVector vec, float& emf, float etad); /* this function returns 
correction based on globaly initialized parameters */ 

  //_____________________________________________________________ Balance of two leading jets
  //------------- !!!!!!!!      OLD stuff, but I still need it  !!!!!!!!!
  double GetBalance2(double pt1, double pt2, double phi1, double phi2);
  double GetBalance2(std::vector<TLorentzVector> vec);

  //_____________________________________________________________ theta* -- dijet scattering angle
  double GetThetaStar(double y1, double y2);
  double GetThetaStar(std::vector<TLorentzVector> vec);

  //__________________________________________________ returns Kt of two leading jets
  double GetMyKtKick(std::vector<TLorentzVector> vec);
  
  //__________________________________________________ returns Phi of two leading jets
  double GetMyPhiKt2L(std::vector<TLorentzVector> vec);
  //__________________________________________________ returns Phi of sumKt of all jets 
  double GetMyPhiKtAll(std::vector<TLorentzVector> vec);  
  //__________________________________________________ returns 1 if jet is matched
  int MyMatchPhoJet(int jet_ind, int pho_ind, int ele_ind, TStnJetBlock* fJetBlock,
		    TLorentzVector *mypho, TLorentzVector *myjet, double algo_cone, 
		    double &dR_fj, double &dPhi_fj, double &dEta_fj);
  //________________________________ returns Qt of all jets above Eta & Et thresholds
  double MyQtJet(JetStuff jetstuff);
  //________________________________ returns SumEt=jets+photons+leptons+taus+met (jets above Eta & Et thresholds)
  //                                 metscenario= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
  double MyHtAll(JetStuff jetstuff, CommonStuff miscstuff, int metscenario);
  //_________________________________ returns normalized kT_perp. for resolution studies 
  float GetMyKtPerp(std::vector<TLorentzVector> vec);  
  //_________________________________ returns normalized kT_parl. for resolution studies 
  float GetMyKtParl(std::vector<TLorentzVector> vec);  

  ClassDef(TMyJetFilterModule,0)
};

#endif
