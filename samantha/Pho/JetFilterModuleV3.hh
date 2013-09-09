#ifndef JETFILTERV3_HH
#define JETFILTERV3_HH

/**********************************************************
 * This is a strip down version of JetFilterModuleV2. I
 * removed all the hists and stuff to run faster in Flat ntuple
 * making.
 *********************************************************/
/*{{{*/
/*
 * $Id: JetFilterModuleV3.hh,v 1.1 2011/05/25 21:38:36 samantha Exp $
 * $Log: JetFilterModuleV3.hh,v $
 * Revision 1.1  2011/05/25 21:38:36  samantha
 * This is a strip down version of JetFilterModuleV2. I removed all the
 * hists and stuff to run faster in Flat ntuple making.
 *
 */
/*}}}*/


#if !defined (__CINT__) || defined (__MAKECINT__)

#include "TF1.h"
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include <vector>
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

#endif

//sam's stuff

#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/Pho/TagLoosePhotons.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "samantha/Pho/TagLooseElectrons.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/TriggerModule.hh"
#include <Stntuple/obj/TStnMetBlock.hh>
#include "samantha/obj/Histograms.hh"
#include "Stntuple/obj/TGenpBlock.hh"

class JetFilterModuleV3: public TStnModule {
	public:

		//to control the print statements
		enum {
			PRN_INFO = 30,
			PRN_WARN = 20,
			PRN_FATAL = 10
		};

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
			int myNjet_th30; // number of level-6 jets with Et>30, after matching with EM object
			int myNjet_th35; // number of level-6 jets with Et>35, after matching with EM object
			//below are the number of smeared jets for given threshold for pseudo expts only. this is reused. so do not use this for any real purpose , 12-28-2008, sam
			// smeared jets are lev6 no em obj jets with added MET and then smeared by the smear factors - 1-14-2008,sam 
			int smeared_Njet_th0;
			int smeared_Njet_th5;
			int smeared_Njet_th10;
			int smeared_Njet_th15;
			int smeared_Njet_th20;
			int smeared_Njet_th25;
			int smeared_Njet_th30;
			int smeared_Njet_th35;

			std::vector<TLorentzVector> Jet_raw;      // uncorrected jets 
			std::vector<TLorentzVector> Jet_lev1;     // jets corrected to level-1
			std::vector<TLorentzVector> Jet_lev4;     // jets corrected to level-4
			std::vector<TLorentzVector> Jet_lev5;     // jets corrected to level-5
			std::vector<TLorentzVector> Jet_lev6;     // jets corrected to level-6
			std::vector<TLorentzVector> Jet_lev7;     // jets corrected to level-7
			std::vector<TLorentzVector> Jet_raw_noEMobj;      // uncorrected jets, matched EM object removed 
			std::vector<TLorentzVector> Jet_lev1_noEMobj;     // jets corrected to level-1, matched EM object removed
			std::vector<TLorentzVector> Jet_lev4_noEMobj;     // jets corrected to level-4, matched EM object removed
			std::vector<TLorentzVector> Jet_lev5_noEMobj;     // jets corrected to level-5, matched EM object removed
			std::vector<TLorentzVector> Jet_lev6_noEMobj;     // jets corrected to level-6, matched EM object removed
			std::vector<TLorentzVector> Jet_lev7_noEMobj;     // jets corrected to level-7, matched EM object removed
			std::vector<TLorentzVector> newJetLev6;           // pseudo-jets for MetModel: noEM_lev6+MET*cos(dPhi)
			std::vector<TLorentzVector> smeared_newJetLev6_noEMObj;   // smeared pseudo-jets for MetModel: noEM_lev6 //sam added on 1-1-2009
			std::vector<bool> bvMetAdded; // keep track of the Jet/s that got Met Added 0=no met added 1=met is added
			std::vector<int> iL6MatchedHadJetIndex; // Index of the HAD jet matched to a Lev6NoEmObj jet
			std::vector<double> smear_factor; // smear factor for each jet (for pseudo-experiments) 
			std::vector<float> fL6NoEm_MisMeasureProb;          //for MC, calculated prob. for each jet to be mismeasured
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
			std::vector<int> tightId;		//tight photon ID (=0 if all cuts were passed)
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
			std::vector<double> myEleXces;    // X_ces for electrons
			std::vector<double> myEleZces;    // Z_ces for electrons
			std::vector<TLorentzVector> myRawMuon;      // raw muons
			std::vector<TLorentzVector> myCorrMuon;     // corrected muons (corrections???)
			std::vector<TLorentzVector> myRawTau;       // raw taus
			std::vector<TLorentzVector> myCorrTau;      // corrected taus  (corrections???, for consistency)
			std::vector<TLorentzVector> myRawBjet;      // raw b-jets
			std::vector<TLorentzVector> myCorrBjet;     // corrected b-jets  (b-specific corrections???, for consistency)
			TVector2 myMET_raw; // MEt(4) before any corrections
			TVector2 myMET0_raw; // MEt(0) before any corrections
			double mySumEt_raw; // SumEt before any corrections
			double newVertexZ; // new vertex position after swap
		};

		struct MatchStuff {
			int JetInd_match;
			int EmInd_match;
			int EmObjType_match; // pho=0, ele=1
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
			TH1F*  fEvntNjet30_b;
			TH1F*  fEvntNjet35_b;

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
			TProfile*  fEvntNJet30_Nvx_b;        // Njet30 vs Nvx(class>12) 
			TProfile*  fEvntNJet35_Nvx_b;        // Njet35 vs Nvx(class>12) 

			TProfile*  fEvntKt_Mjj_b;            // "Kt-kick" of dijet system vs. Mjj
			TProfile*  fEvntKtAll_Mjj_b;         // "Kt" of extra vs. Mjj
			TProfile*  fEvntEt0_Njet_b;          // Average raw Et of jets vs. Njet 
			TProfile*  fEvntEt_Njet_b;           // Average corrected Et of jets vs. Njet 
			TProfile*  fEvntEt0_Nvx12_b;         // Average raw Et of jets vs. Nvx(class>=12) 
			TProfile*  fEvntEt_Nvx12_b;          // Average corrected Et of jets vs. Nvx(class>=12)                                          
			TProfile*  fEvntNjet_Lum_b;          // number of jets vs. Inst. lum for Nvx12=1

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
			TH1F*  fEvntNjet30_a;
			TH1F*  fEvntNjet35_a;

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
			TProfile*  fEvntNJet30_Nvx_a;         
			TProfile*  fEvntNJet35_Nvx_a;         

			TProfile*  fEvntKt_Mjj_a;    
			TProfile*  fEvntKtAll_Mjj_a; 
			TProfile*  fEvntEt0_Njet_a;  
			TProfile*  fEvntEt_Njet_a;   
			TProfile*  fEvntEt0_Nvx12_a; 
			TProfile*  fEvntEt_Nvx12_a;  
			TProfile*  fEvntNjet_Lum_a;  

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

		enum { myNjetArray = 100 };
	protected:

		char histo_outputstring[200];

		int EventCount_b; // global event count before cuts
		int EventCount_a; // global event count after cuts
		//____________________________________   pointers to the data blocks used,
		// header block is always available via
		// TStnModule::GetHeaderBlock()
		TStnJetBlock*      fJetBlockClu04;
		TStnJetBlock*      fHadJetBlockClu04;
		TStnJetBlock*      fJetBlockClu07;
		TStnJetBlock*      fJetBlockClu10;
		TStnMetBlock*      fMetBlock;		// sam added
		TStnHeaderBlock*    fHeaderBlock;	// sam added
		TGenpBlock* 		  fGenpBlock;

		//----------------------- need this for jet-EM matching
		TCalDataBlock* fCalData;
		TStnPhotonBlock* fPhotonBlock;
		TStnElectronBlock* fElectronBlock;
		typedef std::vector<TCalTower*> CalDataArray; // need for beam halo
		typedef std::vector<TCalTower*>::iterator CalDataArrayI; // need for beam halo


		//______________________________________  ***** histograms filled
		JetGeneral_t       fHistJet04; // general histograms for JetClu 0.4

		//_____________________________________ histograms for study of Pho-Jet matching
		MatchStudyHisto_t fMatchPhoJet04; 

		//______________________________________ ***** jet params for different jet cones
		JetStuff jet04stuff;
		CommonStuff allstuff;
		std::vector<MatchStuff> matchstuff;

		int myJetRun;         // to be used in corrections and for crosschecks
		int myNvx_class12;     // Nvx(class12) to be used in corrections

		//_____________________________________________________________________________________________________
		int fJetAlgo;      // identifier of jet algorithm: 0-JetClu(0.4), 1-JetClu(0.7), 2-JetClu(1.0)
		//                              3-MidPoint(0.4), 4-MidPoint(0.7), 5-MidPoint(1.0)
		//                              6-KtCl(0.4), 7-KtCl(0.7), 8-KtCl(1.0)
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

		int bad_EMjet_match_flag; // global flag to mark events where not all EM objects are matched to jets

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
		int MinNjet30; // cut on the number of jets with Et(lev6)>30 after removing EM objects
		int MaxNjet30;
		int MinNjet35; // cut on the number of jets with Et(lev6)>30 after removing EM objects
		int MaxNjet35;
		int fRemoveDuplicate; // switch to remove duplicate objects

	public:
		JetFilterModuleV3(const char* name ="JetFilterV3", 
				const char* title="JetFilterV3");
		~JetFilterModuleV3();

		//_______________________________________________________ ****** accessors

		JetGeneral_t*       GetHistJet04()          { return &fHistJet04; }

		TStnJetBlock*       GetJetBlockClu04() { return fJetBlockClu04; }
		TStnJetBlock*       GetJetBlockClu07() { return fJetBlockClu07; }
		TStnJetBlock*       GetJetBlockClu10() { return fJetBlockClu10; }
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }

		int    GetMyMetScenario() { return fMyMetScenario; }
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
		int GetMinNjet30() { return MinNjet30; }
		int GetMaxNjet30() { return MaxNjet30; }
		int GetMinNjet35() { return MinNjet35; }
		int GetMaxNjet35() { return MaxNjet35; }

		int GetRemoveDuplicate() { return fRemoveDuplicate; }
		//__________________________________________________________ accessors to jet 
		//__________________________________________________________ correction params
		int    GetJTC_coneSize()       { return fJTC_coneSize; } 
		int    GetJTC_version()  	 { return fJTC_version; }
		int    GetJTC_level()    	 { return fJTC_level; }  
		int    GetJTC_systcode()	 { return fJTC_systcode; }
		int    GetJTC_imode()	         { return fJTC_imode; }


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
		int GetMyJetBlockInd(int cone, int i); // original jet index in JetBlock, need this after jet reordering
		int GetMyNjet(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25, 5--et>30, 5--et>35
		double GetMySumEtCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double GetMyHtCorr(int cone, int threshold, const int iJetCorrLev=6); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20  // iJetCorrLev=forces to used jets with different corrections when calculating Ht. Default is lev6NOEM.
		double GetMyMetCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double GetMyMetXCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double GetMyMetYCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double GetMyMetPhiCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		//________________________________________________________ ****** setters

		void SetMyMetScenario(int cut) { fMyMetScenario = cut; }
		void SetJTC_coneSize    (int param)           { fJTC_coneSize = param; } 
		void SetJTC_version     (int param)           { fJTC_version = param; }
		void SetJTC_level       (int param)           { fJTC_level = param; }  
		void SetJTC_systcode    (int param)	        { fJTC_systcode = param; }
		void SetJTC_imode       (int param)	        { fJTC_imode = param; }


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
		void SetMinNjet30(int param) { MinNjet30=param; }
		void SetMaxNjet30(int param) { MaxNjet30=param; }
		void SetMinNjet35(int param) { MinNjet35=param; }
		void SetMaxNjet35(int param) { MaxNjet35=param; }

		void SetRemoveDuplicate(int param) { fRemoveDuplicate=param; }

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
		void SetMyNjet(int cone, int threshold, int param); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25, 5--et>30, 5--et>35
		void SetMySumEtCorr(int cone, int threshold, double param); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		void SetMyMetCorr(int cone, int threshold, double px, double py); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		void SetMyMetCorr(int cone, int threshold, const TVector2& vec); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20

		//____________________________________________ Service methods


		void ClearModuleOutput(); // clears the Module output  		 
		void ClearJetStuff(JetStuff &stuff);     // clears the JetStuff structures 		 
		void ClearCommonStuff(CommonStuff &stuff); // clears CommonStuff structures 
		void ClearMatchStuff(); // clears information about pho-jet matching

		//___________________________________________        TStnModule
		int     BeginJob();
		int     BeginRun();
		int     Event   (int ientry);
		int     EndJob  ();
		//___________________________________________ ****** other methods
		void BookHistograms();
		void BookJetHistograms(JetGeneral_t& Hist, const char* Folder, const char* algoname); // booking general jet histo
		void BookMatchStudyHistograms(MatchStudyHisto_t& Hist, const char* Folder, const char* algoname); // booking match histo

		void FillJetHistogramsB(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12); // filling general histo for particular Jet Cone
		void FillJetHistogramsA(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12); // filling general histo for particular Jet Cone
		void FillMatchingHistograms(MatchStudyHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, MatchStuff match);// filling match histo

		int MyDuplicateCut(CommonStuff miscstuff);
		int MyNjetCut(JetStuff jetstuff);
		void GetNjets(const std::vector<TLorentzVector>& JetVec, std::vector<int>& njets);
		int GetNjets(const std::vector<TLorentzVector>& JetVec, const float fJetEt);
		int MyGenericCut(JetStuff jetstuff, double dz);
		int MyAngularCut(JetStuff jetstuff, CommonStuff miscstuff);

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

		double MyInCemTowerEta(double eta_det); // returns relative position of EM shower inside CEM tower

		double GetCorrection(CorInit settings, TLorentzVector vec, float& emf, float etad); /* this
																															function returns correction based on parameters localy set in the Module */ 
		double GetCorrection(TLorentzVector vec, float& emf, float etad); /* this function returns 
																									correction based on globaly initialized parameters */ 

		int MyMatchPhoJet(int jet_ind, int pho_ind, int ele_ind, TStnJetBlock* fJetBlock,
				TLorentzVector *mypho, TLorentzVector *myjet, 
				double &dR_fj, double &dPhi_fj, double &dEta_fj);
		//________________________________ returns SumEt=jets+photons+leptons+taus+met (jets above Eta & Et thresholds)
		//                                 metscenario= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double MyHtAll(std::vector<TLorentzVector> vJets, JetStuff jetstuff, CommonStuff miscstuff, int metscenario);

		/// sam's stuff
		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int evtsMoreHad;
		};

		int GetNMatchedObj() const { return sam_matchstuff.size(); }
		int GetMatch_JetIndex(int ObjType, int EmInd) const;//Returns original jet block index of the jet matched to an EM object
		void SetDebug(const int num) { debug = num; }
		int Debug() const { return debug; }
		void DebugInfo(const int iAnaMode, const std::string sFunction, const int printlvl, const std::string message);

		int GetMyNCorrPho() const {return allstuff.myCorrPhoton.size(); }
		int GetMyNCorrEle() const {return allstuff.myCorrElectron.size(); }
		
		TLorentzVector GetMyCorrPhoton(const int i) const
		{
			if (i>=0 && i< allstuff.myCorrPhoton.size()) return allstuff.myCorrPhoton.at(i);
		}
		int GetMyCorrPhoInd(const int i) const 
		{
			if (i>=0 && i< allstuff.myPhoInd.size()) return allstuff.myPhoInd.at(i);
		}
		double GetMyCorrPhoDetEta(const int i) const 
		{
			if (i>=0 && i< allstuff.myPhoEtaDet.size()) return allstuff.myPhoEtaDet.at(i);
		}
		TLorentzVector GetMyCorrElectron(const int i) const
		{
			if (i>=0 && i< allstuff.myCorrElectron.size()) return allstuff.myCorrElectron.at(i);
		}
		int GetMyCorrEleInd(const int i) const 
		{
			if (i>=0 && i< allstuff.myEleInd.size()) return allstuff.myEleInd.at(i);
		}
		double GetMyCorrEleDetEta(const int i) const 
		{
			if (i>=0 && i< allstuff.myEleEtaDet.size()) return allstuff.myEleEtaDet.at(i);
		}


		void DumpJets(const std::string func_name, const unsigned int line_num,
							const JetStuff&, const int jetLev);
		void DumpEMobjects(const std::string func, const unsigned int line, const CommonStuff&, const int type=0);
		
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
		void DumpTLvec(const std::string func, const unsigned int line,
								const std::vector<TLorentzVector> vObj , const std::string label="");
		void SetPrintLevel(const int p) { iPrintLvl = p; }
		int PrintLevel() const { return iPrintLvl; }

		void SetRemoveTightPhotons(const bool b) { bRem_TightPho = b; }
		void SetRemoveLoosePhotons(const bool b) { bRem_LoosePho = b; }
		//void SetRemoveSidebandPhotons(const bool b) { bRem_SidebandPho = b; }
		void SetRemoveTightPhoLikeElectrons(const bool b) { bRem_TightPhoLikeEle = b; }
		void SetRemoveLoosePhoLikeElectrons(const bool b) { bRem_LoosePhoLikeEle = b; }
		void SetRemoveStdLooseElectrons(const bool b) { bRem_StdLooseEle = b; }

		bool RemoveTightPho() const { return bRem_TightPho; }
		bool RemoveLoosePho() const { return bRem_LoosePho; }
		//bool RemoveSidebandPho() const { return bRem_SidebandPho; }
		bool RemoveTightPhoLikeEle() const { return bRem_TightPhoLikeEle; }
		bool RemoveLoosePhoLikeEle() const { return bRem_LoosePhoLikeEle; }
		bool RemoveStdLooseEle() const { return bRem_StdLooseEle; }

		
	private:
		Global_Counters_t mycounter;
		bool bRunPermit;		//__ make sure we have all dependencies met before running
		InitSuperPhotons *initSpMod;
		TriggerModule *trigMod;
		TagTightPhotons *tightMod;
		TagLoosePhotons *looseMod;
		TagElectrons *tightEleMod;
		TagLooseElectrons *looseEleMod;

		std::vector<MatchStuff> sam_matchstuff;  // to get the correct match jet ind. 

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		int debug;								// control the print statements
		bool bSaveThisEvent;			// to skim few stntuples for debugging
		bool bUseHadJets;				// 1=replaces Detector jets with matching HAD jets. 0=does nothing
		bool bNoSummary;				// control the end job summary print
		int iPrintLvl;					//controls the print statements	

		//a very small eta value to remove the warning from Root when Pt=0 and called for Eta()
		float fMinisculeEta;
		
		// decides what types of EM objects are removed from the jet list.
		bool bRem_TightPho, bRem_LoosePho, bRem_SidebandPho, bRem_TightPhoLikeEle;
		bool bRem_LoosePhoLikeEle, bRem_StdLooseEle;


		ClassDef(JetFilterModuleV3,0)
};

#endif
