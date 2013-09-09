#ifndef JETFILTERV2_HH
#define JETFILTERV2_HH


/*
 * $Id: JetFilterModuleV2.hh,v 1.43 2009/11/10 19:28:41 samantha Exp $
 * $Log: JetFilterModuleV2.hh,v $
 * Revision 1.43  2009/11/10 19:28:41  samantha
 * Uncommented the old GetMyHtCorr() as I am not using the modified MEtModel stuff
 * anymore.
 * ADDED: GetMyJetNoEMobj_lev7() to get lev7 corrected jets.
 * MOVED: GetMatch_JetIndex() implementation to cc file
 * DELETED: unused 'unsigned iSavedEvents'
 *
 * Revision 1.42  2009/08/25 16:17:22  samantha
 * ADDED:1:Few multi-dimensional std::vectors to hold the JER parameters.
 * The 2D vectors seems to work fine. But when I tried to create 3D vector
 * the same, and the compiler says it can't find it and gives linking error.
 * So I have commented it out. Other vectors are not used in the cc yet.
 *       2:There is a 3D hist I added to debug list which I don't remember why?
 *
 * Revision 1.41  2009/08/07 20:50:35  samantha
 * ADDED: A Photon energy resolution plot to the standard list of analysis hists.
 * NOT COMPLETE: Created two 2-D vectors to hold the JER parameters. This way,accessing
 * JER parameters for JER reconstruction would be much faster(I think). Right now,
 * we are creating and destroying this JER parameter array for every iteration.
 * I have to figure out a way to read them in once within the class and keep it
 * memory. Right now they are not part of the class. If I take them in
 * the JER recontruction runs into trouble as they are not included in class either.
 *
 * Revision 1.40  2009/07/24 23:07:28  samantha
 * ADDED: 1: Choose between old and new parameterization for individual final 5 fit
 * parameters (Gauss Mean, Gauss Sigma, Landau MPV, Landau Sigma, and Gauss Norm).
 * Many enums to make things streamline and simple.
 *
 * Revision 1.39  2009/07/21 16:50:59  samantha
 * MODIFIED: Name of the inclusion gaurd, TMYJETFILTERV2_HH -> JETFILTERV2_HH.
 * DELETED: 1: Xsection hists/methods and  and the Vtx swapping stuff.
 *
 * Revision 1.38  2009/06/26 23:01:34  samantha
 * ADDED:  1: Option to remove only sideband photons from jet list. Not complete
 * yet.
 *
 * Revision 1.37  2009/06/18 22:25:43  samantha
 * MAJOR CHANGES: Controls the EM objects that get removed from the jet list.
 *
 * Revision 1.36  2009/06/16 04:28:50  samantha
 * ADDED: 	1: enum to control print statements.
 * 	2: enum to decide how I 'll add MEt to a jet.
 * 	3. Set/Get to choose how we add MEt to a jet.
 *
 * Revision 1.35  2009/06/03 01:32:35  samantha
 * ADDED: 	1: A dalitz plot to the std. set of hists for g+>=2 jets case
 * 	2: New set of debug hists to look at Photon energy resolution effect
 * 	on the DelPhi(Met,Lead Jet) plot. See Elog#1203, 1204, 1207, 1208.
 * 	3: Method to find a matching HEPG particle for a given detector particle
 * 	(4-vec).
 * 	4. Method to find the matching HAD jets for given det had jet(4-vec)
 * 	5. iPrintLvl to control the amount of debugging info printed out. Higher
 * 	number means more print statements.
 * DELETED: 1: All the hDbl_xxxx debugging hists.
 *
 * Revision 1.34  2009/05/25 18:32:38  samantha
 * ADDED: A hard MetCut option.
 *
 * Revision 1.33  2009/05/21 18:53:47  samantha
 * MODIFIED: LeastMismeasuredJet() is introduced to find the least mismeasured jet
 * index from a given set of jets.
 *
 * Revision 1.32  2009/05/13 22:32:37  samantha
 * MAJOR CHNAGES: 1. Added stuff to smear the jets for pseudo experiments twice.
 * 2. Modified things to use the same hist fill routine FillAnalysisHistograms()
 * for both DATA and pseudo experiments. 3. MyHtAll is modified to be used in both
 * DATA and BG. 4. JER off_set is setup as an option now.
 * ADDED:  1. Switch to enable/disbale additional smearing and function Smear_newJetL6FirstTime().
 * 2. Added 2nd Jet Et hist and Eta hists for two lead jets.  3. GetNjets() returns number of
 * jet above a given Pt threshold for given set of jets (a vector with 4-momenta) 4. Expanded the
 * DumpJets() to dump more info about other jets.
 * DELETED: 1. MyMaddCut() which is not used.  2. FillBackgroundAnalysisHistograms()
 * which is now reaplced by FillAnalysisHistograms(). 3. Some of the debug hists and stuff in
 * FillAnalysisHistograms() and DoMyAnalysis()  4. SmearedHtAll() is deleted and MyHtAll
 * is modified to use in all occasions.
 *
 * Revision 1.31  2009/05/07 02:20:32  samantha
 * ADDED: photon tight IDs so I can make sure I use only leading tight photon when I fill final hists. Temporary counter to count the events with more HAD jets that Detector jets. In such cases I use only the leading HAD jets which poses a problem for MEt and Ht calculation etc.
 *
 * Revision 1.30  2009/04/17 18:03:20  samantha
 * DELETED:  ReplaceLev6JetWithHadJets(). CHANGED: MyNewJetWithMet() to pass in EM object info. ADDED: Dump_newJetLev6() to dump MetAdded jets.
 *
 * Revision 1.29  2009/04/14 20:40:13  samantha
 * CHNAGED: redeclare a loop iterator from int to 'unsigned int' to remove the warninig from compiler. DELETED: some hists that are not used like the generated JerRandom point distribution (complete list - hDbl_JerRandom[6], hDbl_JerRandom_JetE[6], hDbl_JerRandom_JetDetEta[6], smearfacVsGenJer,smearfacVsJetE,smearfacVsGenJerRatio,smearfac_withOffset,offsets,JerRandom.)
 *
 * Revision 1.28  2009/04/10 20:58:15  samantha
 * ADDED: hDbl_MetAddedJetMisMeasureProb_ByJetEt to plot the mismeasure probability of the MEt added jet by Jet Et bins. See Elog#1146. RENAMED: MostMismeasuredProb() to JetMismeasureProb().
 *
 * Revision 1.27  2009/04/07 17:54:02  samantha
 * MAJOR CHANGES: Now I am adding MEt to the least mismeasured jet. Results are in elog#1125.  ADDED: :MostMisMeasuredJet() to measure the mismeasure probability of a jet. Information is stored for later used under JetStuff.
 *
 *
 */

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "TF1.h"
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include <vector>
//#include "CLHEP/Vector/LorentzVector.h" //??????????????????? sam
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

class JetFilterModuleV2: public TStnModule {
	public:

		//to control the print statements
		enum {
			PRN_INFO = 30,
			PRN_WARN = 20,
			PRN_FATAL = 10
		};

		// method to add MEt to a jet
		enum {
			DO_NOT_ADD_MET = 0,
			TOMOST_MISMEASURED_JET = 1,
			TO_JET_CLOSEST_TO_MET = 2
		};
	
		// to select different JER functions
		enum JERFuncTypes_t
		{
			JER_OLDFITS     = 0,
			JER_GAUSSMEAN   = 1,
			JER_GAUSSSIGMA  = 2,
			JER_LANDAUMPV   = 3,
			JER_LANDAUSIGMA = 4,
			JER_GAUSSNORM   = 5
		};

		//set up JER to use old/new fit functions
		struct JERFuncOldNew_t
		{
			unsigned int GaussMean_New;
			unsigned int GaussSigma_New;
			unsigned int LandauMPV_New;
			unsigned int LandauSigma_New;
			unsigned int GaussNorm_New;
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

		//------------------- parameters for metsig studies
		struct MetSigDelPhiStuff {
			double dPhi_det[5]; // i=0- closest obj; i=1- closest EM obj; i=2- closest jet; i=3- 1st EM; i=4- 1st jet.
			double eta_det_em; // detEta of closest EM object
			double eta_det_jet; // detEta of closest jet
			double xces_em; // X_ces of closest EM    
		};

		//------------------- parameters for Met Probability calculation
		struct MetProbStuff {
			double max_MetProbInteg; // upper limit on MetProbability to have MET(fluct)<MET(measured)
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

		struct MetResults {

			//_______ data: MetSig>cut
			//__ [jetThr][Njet]; jetThr=Njet15,Njet20,Njet25; Njet=0,1,...9
			int ana20Njet_dt[3][10]; // DATA: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
			int ana25Njet_dt[3][10]; // DATA: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
			int ana30Njet_dt[3][10]; // DATA: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
			int ana35Njet_dt[3][10]; // DATA: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
			int ana40Njet_dt[3][10]; // DATA: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
			int ana45Njet_dt[3][10]; // DATA: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
			int ana50Njet_dt[3][10]; // DATA: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
			int ana75Njet_dt[3][10]; // DATA: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
			int ana100Njet_dt[3][10]; // DATA: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
			int ana150Njet_dt[3][10]; // DATA: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25

			//______ background: MetSig(pseudo-exp)>cut
			//__ [20]=0-def; 1=z->ee vs. dipho sideband; 2=UnclEn:mean-G; 3=UnclEn:mean+G; 4=UnclEn:sigma-G; 5=UnclEn:sigma+G;
			//__  6=UnclEn:scale-G; 7=UnclEn:scale+G; 8=UnclEn:norm-G; 9=UnclEn:norm+G;           
			//__  10=JER:meanG-G; 11=JER:meanG+G; 12=JER:sigmaG-G; 13=JER:sigmaG+G;
			//__  14=JER:mpvL-G; 15=JER:mpvL+G; 16=JER:sigmaL-G; 17=JER:sigmaL+G;
			//__  18=JER:norm-G; 19=JER:norm+G
			int ana20Njet_bg[20][3][10]; // Met Model: number of events with Met>20 GeV in bins of Njet15, Njet20, Njet25
			int ana25Njet_bg[20][3][10]; // Met Model: number of events with Met>25 GeV in bins of Njet15, Njet20, Njet25
			int ana30Njet_bg[20][3][10]; // Met Model: number of events with Met>30 GeV in bins of Njet15, Njet20, Njet25
			int ana35Njet_bg[20][3][10]; // Met Model: number of events with Met>35 GeV in bins of Njet15, Njet20, Njet25
			int ana40Njet_bg[20][3][10]; // Met Model: number of events with Met>40 GeV in bins of Njet15, Njet20, Njet25
			int ana45Njet_bg[20][3][10]; // Met Model: number of events with Met>45 GeV in bins of Njet15, Njet20, Njet25
			int ana50Njet_bg[20][3][10]; // Met Model: number of events with Met>50 GeV in bins of Njet15, Njet20, Njet25
			int ana75Njet_bg[20][3][10]; // Met Model: number of events with Met>75 GeV in bins of Njet15, Njet20, Njet25
			int ana100Njet_bg[20][3][10]; // Met Model: number of events with Met>100 GeV in bins of Njet15, Njet20, Njet25
			int ana150Njet_bg[20][3][10]; // Met Model: number of events with Met>150 GeV in bins of Njet15, Njet20, Njet25

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

		struct MetStudyHisto_t {  //------------ Met study histograms 

			//------------ Met significance histograms
			TH2F* fMetSig_vs_sqrtHt; // metSig vs. sqrt(Ht)
			TH2F* fMetSig_vs_sqrtMet; // metSig vs. sqrt(Met)
			TH2F* fMetSig_vs_sqrtSumEt; // metSig vs. sqrt(corSumEt)

			TH2F* fMetSig_vs_dPhi[5]; // metSig vs. dPhi[i]
			TH2F* fMetSig_vs_Xces; // metSig vs. Xces for closest photon/electron if dPhi<0.4
			TH2F* fMetSig_vs_etaEM; // metSig vs. eta_det for closest photon/electron if dPhi<0.4
			TH2F* fMetSig_Xces_vs_etaEM; // Xces vs. relative_eta_det for closest photon/electron if dPhi<0.4
			TH2F* fMetSig_vs_etaJET; // metSig vs. eta_det for closest jet if dphi<0.4
			TH2F* fMetSigCalib_vs_dPhi[5]; // pseudo-experiments: metSig vs. dPhi[i] for maximum genMET
			// i=0- closest obj; i=1- closest EM obj; i=2- closest jet; i=3- 1st EM; i=4- 1st jet.
			TH1F* fObjType; // type of closest to MET object in events with significant MET 
			TH1F* fMetSig_corrFactor; // MetSig correction factor; based on pseudo-experiments

			TH1F* fMetSigCalib_estimate; // for MET from pseudo-experiments: estimated upper limit on MetSig
			TH1F* fMetSigCalib_estimate_njet0; // for MET from pseudo-experiments: estimated upper limit on MetSig if Njet15=0
			TH1F* fMetSigCalib_estimate_njet1; // for MET from pseudo-experiments: estimated upper limit on MetSig if Njet15=1
			TH1F* fMetSigCalib_estimate_njet2; // for MET from pseudo-experiments: estimated upper limit on MetSig if Njet15=2
			TH1F* fMetSigCalib_estimate_njet3; // for MET from pseudo-experiments: estimated upper limit on MetSig if Njet15=3
			TH1F* fMetSigCalib_estimate_njet4; // for MET from pseudo-experiments: estimated upper limit on MetSig if Njet15>=4
			TH1F* fMetSigCalib_Met[5]; // for MET from pseudo-experiments: Met for 0==all events; i== -log10(1-MetProb)>i

			TH1F* fMetSigToy_Met[5]; // for toy MET: ToyMet for 0==all events; i== -log10(1-MetProb)>i, i=3,4,5,6
			TH1F* fMetSigToy_eff[4]; // for toy MET: efficiency of MetSig cut: i== -log10(1-MetProb)>i, i=3,4,5,6
			TH1F* fMetSigToy_estimate; // for toy MET: estimated upper limit on MetSig

			TH1F* fMetSig_estimate; // estimated upper limit on MetSig
			TH1F* fMetSig_estimate_njet0; // estimated upper limit on MetSig if Njet15=0
			TH1F* fMetSig_estimate_njet1; // estimated upper limit on MetSig if Njet15=1
			TH1F* fMetSig_estimate_njet2; // estimated upper limit on MetSig if Njet15=2
			TH1F* fMetSig_estimate_njet3; // estimated upper limit on MetSig if Njet15=3
			TH1F* fMetSig_estimate_njet4; // estimated upper limit on MetSig if Njet15>=4
			TH1F* fMetSig_Met[5]; // Met for 0==all events; i== -log10(1-MetProb)>i
			//________________________________________

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

			TProfile* fMet4VsNjet; // Met4 vs. Njet (number of clusters)
			TProfile* fMet4VsNjet_th5; // Met4 vs. Njet(lev6>5GeV)
			TProfile* fMet4VsNjet_th10; // Met4 vs. Njet(lev6>10GeV)
			TProfile* fMet4VsNjet_th15; // Met4 vs. Njet(lev6>15GeV)
			TProfile* fMet4VsNjet_th20; // Met4 vs. Njet(lev6>20GeV)
			TProfile* fMet4VsNjet_th25; // Met4 vs. Njet(lev6>25GeV)
			TProfile* fMet4VsNjet_th30; // Met4 vs. Njet(lev6>30GeV)
			TProfile* fMet4VsNjet_th35; // Met4 vs. Njet(lev6>35GeV)

			TProfile* fMet4VsRun; // Met4 vs. Run number, for nvx12=1
			TProfile* fMet4VsNvx12; // Met4 vs. Nvx12
			TProfile* fNvx12VsRun; // Nvx12 vs. Run number

			TH1F* fMet4[10]; // Met4 for events with any Nvx12 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
			TH1F* fMet4Phi[10]; // Phi of Met4 for events with any Nvx12 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
			TH1F* fMet4X[10]; // Met4X for events with any Nvx12 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
			TH1F* fMet4Y[10]; // Met4Y for events with any Nvx12 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
			TH1F* fMet4X_noJet[10]; // Met4X for events with Njet_thr*=0 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
			TH1F* fMet4X_withJet[10]; // Met4X for events with Njet_thr*>0 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
			TH1F* fMet4Y_noJet[10]; // Met4Y for events with Njet_thr*=0 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0
			TH1F* fMet4Y_withJet[10]; // Met4Y for events with Njet_thr*>0 in 10 bins of Sqrt(CorrSumEt); bins size is 2.0

			TH1F* fSumEtJetFrac; // fraction of SumEt carried by jets with Et>threshold
			TH1F* fSumEtJet; // SumEt carried by jets with Et>threshold
			TH1F* fSqrtSumEtJet; // sqrt(SumEt) carried by jets with Et>threshold
			TProfile* fJetFrVsSumEt; // fraction of SumEt carried by jets with Et>threshold vs. total SumEt
			TProfile* fMet4VsJetFr; // meth4 Met .vs. fraction of SumEt carried by jets with Et>threshold
			TProfile* fMet4VsSqrtSumEtJet; // meth4 Met .vs. sqrt(SumEtJet) (SumEt carried by jets with Et>threshold)
			TH2F* fMet4X_vs_Met4Y; // meth4 Met_X vs. Met_Y for all events
			TH2F* fGenMetX_vs_MetY; // generated (def. parametrization) Met_X vs. Met_Y for all events

			//__________________ histograms for MetModel2
			TH1F* fGenV2Met_def_proj; // projection of GenMet on direction of detMet
			TH2F* fGenV2MetXY_def_proj; // GenMet is recalculated respect to detMet
			TH1F* fGenV2Met_def_projNjet15; // projection of GenMet on direction of detMet for events with Njet15>0
			TH2F* fGenV2MetXY_def_projNjet15; // GenMet is recalculated respect to detMet for events with Njet15>0
			TH1F* fGenV2_dPhiMet_def; // dPhi(genMet-detMet)
			TH1F* fGenV2_dPhiMet_def_Njet15; // dPhi(genMet-detMet) for events with Njet15>0

			TH1F* fGenV2Met_max_def_proj; // projection of GenMet on direction of detMet, using maximum generated MET
			TH2F* fGenV2MetXY_max_def_proj; // GenMet is recalculated respect to detMet, using maximum generated MET
			TH1F* fGenV2Met_max_def_projNjet15; // projection of GenMet on direction of detMet for events with Njet15>0, using maximum generated MET
			TH2F* fGenV2MetXY_max_def_projNjet15; // GenMet is recalculated respect to detMet for events with Njet15>0, using maximum generated MET
			TH1F* fGenV2_dPhiMet_max_def; // dPhi(genMet-detMet), using maximum generated MET
			TH1F* fGenV2_dPhiMet_max_def_Njet15; // dPhi(genMet-detMet) for events with Njet15>0, using maximum generated MET

			TH1F* fGenV2Met_max_def_projMet10; // projection of GenMet on direction of detMet>10, using maximum generated MET
			TH2F* fGenV2MetXY_max_def_projMet10; // GenMet is recalculated respect to detMet>10, using maximum generated MET
			TH1F* fGenV2Met_max_def_projMet10Njet15; // projection of GenMet on direction of detMet>10 for events with Njet15>0, using maximum generated MET
			TH2F* fGenV2MetXY_max_def_projMet10Njet15; // GenMet is recalculated respect to detMet>10 for events with Njet15>0, using maximum generated MET
			TH1F* fGenV2_dPhiMet_max_def_Met10; // dPhi(genMet-detMet), if detMet>10 and using maximum generated MET
			TH1F* fGenV2_dPhiMet_max_def_Met10Njet15; // dPhi(genMet-detMet) for events with Njet15>0 and detMet>10, using maximum generated MET

			TProfile* fMetGenV2VsNjet; // generated (def parametrization) Met vs. Njet (number of clusters)
			TProfile* fMetGenV2VsNjet_th5; // generated (def parametrization) Met vs. Njet(lev6>5GeV)
			TProfile* fMetGenV2VsNjet_th10; // generated (def parametrization) Met vs. Njet(lev6>10GeV)
			TProfile* fMetGenV2VsNjet_th15; // generated (def parametrization) Met vs. Njet(lev6>15GeV)
			TProfile* fMetGenV2VsNjet_th20; // generated (def parametrization) Met vs. Njet(lev6>20GeV)
			TProfile* fMetGenV2VsNjet_th25; // generated (def parametrization) Met vs. Njet(lev6>25GeV)
			TProfile* fMetGenV2VsNjet_th30; // generated (def parametrization) Met vs. Njet(lev6>30GeV)
			TProfile* fMetGenV2VsNjet_th35; // generated (def parametrization) Met vs. Njet(lev6>35GeV)

			TProfile* fMetGenV2VsRun; // generated (def parametrization) Met vs. Run number, for nvx12=1
			TProfile* fMetGenV2VsNvx12; // generated (def parametrization) Met vs. Nvx12

			TH1F* fGenV2Met_def;    // generated MET using default parametrization
			TH1F* fGenV2MetX_def;    // generated MET_X using default parametrization
			TH1F* fGenV2MetY_def;    // generated MET_X using default parametrization
			TH1F* fGenV2MetPhi_def; // PhiMET calculated using generated MET and default parametrization

			TProfile* fMetGenV2VsJetFr; // generated (def. parametrization) Met .vs. fraction of SumEt carried by jets with Et>threshold
			TProfile* fMetGenV2VsSqrtSumEtJet; // generated (def. parametrization) Met .vs. sqrt(SumEtJet) (SumEt carried by jets with Et>threshold)
			TH2F* fGenV2MetX_vs_MetY; // generated (def. parametrization) Met_X vs. Met_Y for all events

			TH1F* fMet4GenV2Met; // Met4-MetGen
			TH2F* fMet4_vs_GenV2Met; // Met4 .vs. MetGen
			TH2F* fMet4_vs_Met4GenV2Met; // Met4 .vs. Met4-MetGen

			TH2F* fMetCorr_vs_dZvx; // CorrMet vs. dZ=Zvx1-Zvx2 for events with Nvx12>1
			TH2F* fMetCorr_vs_dZvxWorse; // CorrMet vs. worse dZ for events with Nvx>1
			TH2F* fMetCorr_vs_ZvxBest; // CorrMet vs. Zvx of best vertex
			TH2F* fMetCorr_vs_Zvx2ndBest; // CorrMet vs. Zvx of 2nd best vertex, if Nvx>1
			TH2F* fMetCorr_vs_ZvxWorse; // CorrMet vs. largest Zvx, if Nvx>1 
			TH2F* fMet0Met4_vs_Zvx; // rawMet0-rawMet4 vs. Zvx, if Nvx12=1 
			TH2F* fMet0_vs_Zvx; // rawMet0 vs. Zvx, if Nvx12=1 
			TProfile* fMetCorrVsZvxBest; // CorrMet vs. Zvx of best vertex
			TProfile* fMetCorrVsZvx2ndBest; // CorrMet vs. Zvx of 2nd best vertex, if Nvx>1
			TProfile* fMetCorrVsZvxWorse; // CorrMet vs. largest Zvx, if Nvx>1 

		};

		struct AnalysisHisto_t { // final analysis histograms
			TH1F* fAna_MetAll; // MET for all events
			TH1F* fAna_MetSig; // MET-significance for all events
			TH2F* fAna_MetAllVsMetSig; // MET VS MET-sig  for all events //sam added 11/16/2008


			// the following histrograms are only for events that pass MetSig cut
			TH1F* fAna_Met;		// MET 
			TH1F* fAna_M; 			// invariant mass of di-EM system 
			TH1F* fAna_dPhi; 		// dPhi between two EM objects
			TH1F* fAna_Njet15; 	// Njet(Et>15)
			TH1F* fAna_Njet20; 	// Njet(Et>20)
			TH1F* fAna_Njet25; 	// Njet(Et>25)
			TH1F* fAna_Njet30; 	// Njet(Et>30)
			TH1F* fAna_Njet35; 	// Njet(Et>35)
			TH1F* fAna_Qt; 		// qT of di-EM system
			TH1F* fAna_Et1; 		// Et of first EM object
			TH1F* fAna_Et2; 		// Et of second EM object
			TH1F* fAna_Etjet; 	// Et of 1st jet with Et>15
			TH1F* fAna_Etjet2; 	// Et of 2nd jet with Et>0
			TH1F* fAna_Ht; 		// Ht--sum Et of all objects
			TH1F* fAna_Mjj; 		// Mjj (of two highest Et jets) if Njet15>=2
			TH1F* fAna_Nem; 		// number of EM objects: ele+pho
			TH1F* fAna_Mej; 		// M(em-jet1) if Njet15>=1
			TH1F* fAna_Mextra; 	// M(em1-emN) & M(em2-emN) if Nem>=3
			TH1F* fAna_Etem; 		// Et of extra EM objetcs
			TH1F* fAna_dPhi1; 	// dPhi(met-em1)
			TH1F* fAna_dPhi2; 	// dPhi(met-em2)
			TH1F* fAna_dPhi3; 	// dPhi(met-jet1), Et>15
			TH1F* fAna_dPhi4; 	// dPhi(met-jet1), Et>20
			TH1F* fAna_dPhi5; 	// dPhi(met-jet1), Et>25
			TH1F* fAna_dPhi6; 	// dPhi(met-jet1), Et>30
			TH1F* fAna_dPhi7; 	// dPhi(met-jet1), Et>35
			TH1F* fAna_Jet1Eta;  // lead jet Eta
			TH1F* fAna_Jet2Eta;  // 2nd lead jet Eta
			TH1F* fAna_Pho1Eta;  // lead photon Eta
			TH2F* fAna_Dalitz;  	// 2d plot of InvMass(pho,jet1) and InvMass(pho,jet2)

			//for debugging only

			TH1F* hDebug_PtRatio_JMetdelphiLow;				// Det.Pho/Hepg Pho ratio for events in DelPhi(Jet,MEt)<0.4
			TH1F* hDebug_PhoJetPtRatio_JMetdelphiLow;		// Lead Jet/Hepg Pho ratio for events in DelPhi(Jet,MEt)<0.4
			TH1F* hDebug_PtRatio_JMetdelphiMid;				// Det.Pho/Hepg Pho ratio for events in 1.3<DelPhi(Jet,MEt)<1.7
			TH1F* hDebug_PhoJetPtRatio_JMetdelphiMid;		// Lead Jet/Hepg Pho ratio for events in 1.3<DelPhi(Jet,MEt)<1.7
			TH1F* hDebug_PtRatio_JMetdelphiTop;				// Det.Pho/Hepg Pho ratio for events in DelPhi(Jet,MEt)>2.8
			TH1F* hDebug_PhoJetPtRatio_JMetdelphiTop;		// Lead Jet/Hepg Pho ratio for events in DelPhi(Jet,MEt)>2.8

			TH1F* hDebug_Jet1PtRatio_dphi4;			//lead jet/had jet ratio in DelPhi(Lead Jet,MEt)>2.8 & DelPhi(Jet1,Jet2)<0.4
			TH1F* hDebug_Jet2PtRatio_dphi4;			//sublead jet/had jet ratio in DelPhi(Lead Jet,MEt)>2.8 & DelPhi(Jet1,Jet2)<0.4
			TH1F* hDebug_LeadJet;						// lead jet hist after removing only the lead photon from the jet list

			TH1F* hDebug_PhoRes;							// to study photon energu resolution effect after MEtsig cut

			//2-D plots with MEtSig
			TH2F* fAna_Met_MetSig;			// MET 
			TH2F* fAna_Njet15_MetSig; 		// Njet(Et>15)
			TH2F* fAna_Et1_MetSig; 			// Et of first EM object
			TH2F* fAna_Etjet_MetSig; 		// Et of 1st jet with Et>15
			TH2F* fAna_Etjet2_MetSig; 		// Et of 2nd jet with Et>0
			TH2F* fAna_Ht_MetSig; 			// Ht--sum Et of all objects
			TH2F* fAna_dPhi1_MetSig; 		// dPhi(met-em1)
			TH2F* fAna_dPhi3_MetSig; 		// dPhi(met-jet1), Et>15
			TH2F* fAna_Jet1Eta_MetSig;  	// lead jet Eta
			TH2F* fAna_Jet2Eta_MetSig;  	// 2nd lead jet Eta
			TH2F* fAna_Pho1Eta_MetSig;  	// lead photon Eta
			
			TH3F* hDebug_offsetVsJetE;		// cases when off_set ==0
			//end debugging hists
			
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

		//______________________________________ Met study histograms
		MetStudyHisto_t  fMetStudyJet04sc3; // met study histograms for JetClu 0.4 and Et>15 GeV

		//______________________________________ Final analysis histograms
		AnalysisHisto_t fAna_data; // data
		AnalysisHisto_t fAna_bckg; // Met Model total prediction (stat+syst)
		AnalysisHisto_t fAna_def; // Met Model default prediction
		AnalysisHisto_t fAna_ue1; // Met Model systematics: z->ee vs. dipho sideband parameterization
		AnalysisHisto_t fAna_ue2; // Met Model Uncl. En. systematics: mean-G
		AnalysisHisto_t fAna_ue3; // Met Model Uncl. En. systematics: mean+G
		AnalysisHisto_t fAna_ue4; // Met Model Uncl. En. systematics: sigma-G
		AnalysisHisto_t fAna_ue5; // Met Model Uncl. En. systematics: sigma+G
		AnalysisHisto_t fAna_ue6; // Met Model Uncl. En. systematics: scale-G
		AnalysisHisto_t fAna_ue7; // Met Model Uncl. En. systematics: scale+G
		AnalysisHisto_t fAna_ue8; // Met Model Uncl. En. systematics: norm-G
		AnalysisHisto_t fAna_ue9; // Met Model Uncl. En. systematics: norm+G
		AnalysisHisto_t fAna_jer1; // Met Model JER systematics: meanG-G
		AnalysisHisto_t fAna_jer2; // Met Model JER systematics: meanG+G
		AnalysisHisto_t fAna_jer3; // Met Model JER systematics: sigmaG-G
		AnalysisHisto_t fAna_jer4; // Met Model JER systematics: sigmaG+G
		AnalysisHisto_t fAna_jer5; // Met Model JER systematics: mpvL-G
		AnalysisHisto_t fAna_jer6; // Met Model JER systematics: mpvL+G
		AnalysisHisto_t fAna_jer7; // Met Model JER systematics: sigmaL-G
		AnalysisHisto_t fAna_jer8; // Met Model JER systematics: sigmaL+G
		AnalysisHisto_t fAna_jer9; // Met Model JER systematics: norm-G
		AnalysisHisto_t fAna_jer10; // Met Model JER systematics: norm+G

		//_____________________________________ histograms for study of Pho-Jet matching
		MatchStudyHisto_t fMatchPhoJet04; 

		//_____________________________________ histograms for met cleanup studies
		MetCleanupHisto_t fCleanup;
		MetCleanupHisto_t fGenCleanup;

		//_________________________________ histograms to study correlation between pho and met
		PhoMetStudyHisto_t fPhoMetStudy;

		//______________________________________ ***** jet params for different jet cones
		JetStuff jet04stuff;
		CommonStuff allstuff;
		std::vector<MatchStuff> matchstuff;
		//_____________________________________ ***** met results
		MetResults met_results;
		//_____________________________________ ***** met significance
		MetProbStuff metsigstuff;
		double myMetSig; // significance of MET
		double myMetSig_gen; // significance of MET from pseudo-experimenet
		TVector2 myGenMetVec; // last generated MET from  pseudo-experiment

		//--------------------------------------------------------------------------------------------
		//_________________________________________ ***** service parameters

		int fAnalysisMode; // status code to turn ON/OFF analysis mode 
		int fUseVerbose;   // status code to invoke Debugging mode (print out)
		int fDumpEvent;   // status code to invoke Debugging mode (print out)

		int myJetRun;         // to be used in corrections and for crosschecks
		int myNvx_class12;     // Nvx(class12) to be used in corrections

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


		int fUnclParamSwitch; // switch to set sideband or Z parametrization as default
		// 0=sideband; 1=Z
		int fSampleID; // sample ID: 0=ggX signal; 1=ggX sideband; 2=Z->ee

		int Npoints; // number of points to be generated per value of corrSumEt  
		std::vector<TVector2> MetGen_def; // global vectors to be used in MetPDF
		double MetToy_max; // max value of toy Met
		double MetToy_min; // min value of toy Met
		int fSelectSigMetEvent; // to select events with Significant MET according to MyMetPDF
		int SigMetEvent_status; // global status code for Significant MET according to MyMetPDF
		double fMetSig_cut; // cut on MetSignificance
		double fEventWeight; // event weight for background samples
		int fUseEventWeight; // code to turn ON/OFF event weight

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
		int MinNjet30; // cut on the number of jets with Et(lev6)>30 after removing EM objects
		int MaxNjet30;
		int MinNjet35; // cut on the number of jets with Et(lev6)>30 after removing EM objects
		int MaxNjet35;

		int fSelectMetEvent; // switch to select large met events
		int fSelectExoticEvent; // switch to select exotic events 
		int fRemoveDuplicate; // switch to remove duplicate objects
		int fRemoveBadMet; // switch to remove events with bad met (met along the jet which is close to a crack)  
		int fDoVxSwap; // swap vertices if it minimizes MET  

		char* fDatFileName;
		char* fRearEventFileName;
		char* fExoticEventFileName;
		char* fLargeMetEventFileName;
		char* fDumpEventFileName;

	public:
		JetFilterModuleV2(const char* name ="JetFilterV2", 
				const char* title="JetFilterV2");
		~JetFilterModuleV2();

		//_______________________________________________________ ****** accessors

		JetGeneral_t*       GetHistJet04()          { return &fHistJet04; }
		//   JetGeneral_t*       GetHistJet07()          { return &fHistJet07; }
		//   JetGeneral_t*       GetHistJet10()          { return &fHistJet10; }

		TStnJetBlock*       GetJetBlockClu04() { return fJetBlockClu04; }
		TStnJetBlock*       GetJetBlockClu07() { return fJetBlockClu07; }
		TStnJetBlock*       GetJetBlockClu10() { return fJetBlockClu10; }
		//   TStnJetBlock*       GetJetBlockMid04() { return fJetBlockMid04; }
		//   TStnJetBlock*       GetJetBlockMid07() { return fJetBlockMid07; }
		//   TStnJetBlock*       GetJetBlockMid10() { return fJetBlockMid10; }
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }

		//_______________________________________________________ ****** cut values
		int    GetUseVerbose()           { return fUseVerbose; }
		int    GetDumpEvent()           { return fDumpEvent; }
		int    GetmyJetRun()           { return myJetRun; } 
		int    GetJetAlgo()      	   { return fJetAlgo; }      	    
		int    GetMyMetScenario() { return fMyMetScenario; }
		int    GetSelectSigMetEvent() { return fSelectSigMetEvent; }

		double GetMetSigCut() { return fMetSig_cut; }
		double GetEventWeight() { return fEventWeight; }
		int GetUseEventWeight() { return fUseEventWeight; }
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
		double GetMetToy_max() { return MetToy_max; }
		double GetMetToy_min() { return MetToy_min; }
		int GetSelectMetEvent() { return fSelectMetEvent; }
		int GetSelectExoticEvent() { return fSelectExoticEvent; }
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
		char* GetRearEventFileName() { return fRearEventFileName; }
		char* GetVeryExoticEventFileName() { return fExoticEventFileName; }
		char* GetLargeMetEventFileName() { return fLargeMetEventFileName; }
		char* GetDumpEventFileName() { return fDumpEventFileName; }
		//sam added - 10-08-2008

		int GetNpointsToGenerate() const { return Npoints; }
		int GetUnclParamSwitch() const { return fUnclParamSwitch; }
		int GetAnalysisMode() const { return fAnalysisMode; }
		int GetSampleID() const { return fSampleID; }
		void GenerateSmearedJets(JetStuff& jetstuff);
		void Smear_newJetL6FirstTime(JetStuff& jetstuff, const int systcode, const int metcode,
							const int rnd_seed);

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
		int GetMyNjet(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25, 5--et>30, 5--et>35
		double GetMySumEtCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double GetMyHtCorr(int cone, int threshold, const int iJetCorrLev=6); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20  // iJetCorrLev=forces to used jets with different corrections when calculating Ht. Default is lev6NOEM.
		double GetMyMetCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double GetMyMetXCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double GetMyMetYCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double GetMyMetPhiCorr(int cone, int threshold); // thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double GetMyMetSig() { return myMetSig; } // significance of MET
		double GetMyMetSigGen() { return myMetSig_gen; } // significance of MET from the last pseudo-experiment (default)
		TVector2 GetMyGenMetVec() { return myGenMetVec; } // returns last generated MET from pseudo-experiment
		//________________________________________________________ ****** setters
		void SetEventWeight(double weight) { fEventWeight = weight; }
		void SetUseEventWeight(int cut) { fUseEventWeight = cut; }
		void SetAnalysisMode(int cut) { fAnalysisMode = cut; }
		void SetUseVerbose(int cut) { fUseVerbose = cut; }
		void SetDumpEvent(int cut) { fDumpEvent = cut; }
		void SetmyJetRun(int run) { myJetRun = run; } 
		void SetJetAlgo(int cut) { fJetAlgo = cut; }      	    
		void SetMyMetScenario(int cut) { fMyMetScenario = cut; }
		void SetSelectSigMetEvent(int cut) { fSelectSigMetEvent = cut; }
		void SetMetSigCut(double cut) { fMetSig_cut = cut; }

		void SetJTC_coneSize    (int param)           { fJTC_coneSize = param; } 
		void SetJTC_version     (int param)           { fJTC_version = param; }
		void SetJTC_level       (int param)           { fJTC_level = param; }  
		void SetJTC_systcode    (int param)	        { fJTC_systcode = param; }
		void SetJTC_imode       (int param)	        { fJTC_imode = param; }

		//------------------------ Setting Met Model parameters
		void SetNpointsToGenerate(int n) { Npoints=n; }
		void SetUnclParamSwitch(int n) { if(n>-1 && n<2) fUnclParamSwitch=n; }
		void SetSampleID(int n) { fSampleID=n; }
		//_________________________________________________ end of setting met model params

		void  SetDatFileName(char* fname)       { fDatFileName = fname; }    
		void  SetRearEventFileName( char* fname) { fRearEventFileName = fname; }
		void  SetVeryExoticEventFileName( char* fname) { fExoticEventFileName = fname; }
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
		void SetMinNjet30(int param) { MinNjet30=param; }
		void SetMaxNjet30(int param) { MaxNjet30=param; }
		void SetMinNjet35(int param) { MinNjet35=param; }
		void SetMaxNjet35(int param) { MaxNjet35=param; }
		void SetMetToy_max(double param) { MetToy_max=param; }
		void SetMetToy_min(double param) { MetToy_min=param; }
		void SetSelectMetEvent(int param) { fSelectMetEvent=param; }
		void SetSelectExoticEvent(int param) { fSelectExoticEvent=param; }
		void SetRemoveDuplicate(int param) { fRemoveDuplicate=param; }
		void SetRemoveBadMet(int param) { fRemoveBadMet=param; }
		void SetDoVxSwap(int param) { fDoVxSwap=param; }


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
		void BookPhoMetStudyHistograms(PhoMetStudyHisto_t& Hist, const char* Folder); // booking Pho-Met study histograms
		void BookAnalysisHistograms(AnalysisHisto_t& Hist, const char* Folder); // booking analysis histograms 

		void FillJetHistogramsB(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12); // filling general histo for particular Jet Cone
		void FillJetHistogramsA(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12); // filling general histo for particular Jet Cone
		void FillMetStudyHistograms(MetStudyHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, 
				int metcode, int nvx, double dzvx, double dzvx_worse, double zvx[3]);// filling met histo
		void FillMatchingHistograms(MatchStudyHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, MatchStuff match);// filling match histo
		void FillMetCleanupHistograms(MetCleanupHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, int metcode); // met cleanup studies
		void FillGenMetCleanupHistograms(MetCleanupHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, 
				TVector2 genMet, std::vector<TLorentzVector> vec); // cleanup studies using generated met
		void FillPhoMetStudyHistograms(PhoMetStudyHisto_t& Hist, TVector2 Met, CommonStuff miscstuff, double zvx);  // filling histograms for pho-met study
		void FillAnalysisHistograms(AnalysisHisto_t& Hist,const std::vector<TLorentzVector>, JetStuff, const CommonStuff miscstuff, const TVector2 Met, const float metsig);// filling analysis histograms 
		//void FillBackgroundAnalysisHistograms(AnalysisHisto_t& Hist,const JetStuff& jetstuff,const CommonStuff& miscstuff, const TVector2& Met,const double metsig); //fill background hists with smeared jets in pseudo experiments // sam-1/1/2009

		void DoMyAnalysis(JetStuff jetstuff,CommonStuff miscstuff); // generates Met predictions and fills out analysis histograms
		void DoFinalAnalysisHisto(AnalysisHisto_t& HistF,AnalysisHisto_t Hist1,AnalysisHisto_t Hist2); // finalize analysis histograms;
		void FinalAnalysisHistoStep1(AnalysisHisto_t& Hist); // calls Sumw2 for final analysis histograms;
		void FinalAnalysisHistoStep2(AnalysisHisto_t& bckg, const AnalysisHisto_t& data);// divides the final analysis histograms by Npoints;
		void FinalAnalysisHistoStep3(AnalysisHisto_t& Hist); // divides the final analysis histograms by Npoints;
		double HistoBinDiff(TH1F* h,TH1F* h1,TH1F* h2,int bin_ind); // returns sqrt(stat^2 + max_diff(h-h1,h-h2)^2) in bin=bin_ind
		void DoFinalPhoMetHisto(PhoMetStudyHisto_t& Hist);  // finalize histograms for pho-met study
		void MyEventCount(double met,double metsig,int systcode,MetResults &metstuff); // counts data & background(metmodel) events with MetSig>cut 

		int MyDuplicateCut(CommonStuff miscstuff);
		int MyNjetCut(JetStuff jetstuff);
		int MySmeared_NjetCut(JetStuff jetstuff);    //this is the same as MyNjetCut but for pseudo expts aftert he jets are smeared 12-28-2008, sam
		void GetNjets(const std::vector<TLorentzVector>& JetVec, std::vector<int>& njets);
		int GetNjets(const std::vector<TLorentzVector>& JetVec, const float fJetEt);
		int MyGenericCut(JetStuff jetstuff, double dz);
		int MyAngularCut(JetStuff jetstuff, CommonStuff miscstuff);
		int MyMetCleanUpCut(JetStuff jetstuff, CommonStuff miscstuff, std::vector<TLorentzVector> jetvec, TVector2 metvec); // Test version as of 04/04/07
		int MyExoticEvent(JetStuff jetstuff, CommonStuff miscstuff, int metscenario, int runN, int eventN); // selects my exotic events.
		int MyLargeMetEvent(JetStuff jetstuff, CommonStuff miscstuff, int metscenario, int Nvx12, int runN, int eventN); // selects my large MET events.
		int MyVeryExoticEvent(JetStuff jetstuff, CommonStuff miscstuff, int metscenario, 
				int runN, int eventN, int nvx12, std::vector<double> zvx_vec); // selects my very exotic events.
		void MyDumpEvent(JetStuff jetstuff, CommonStuff miscstuff, int Nvx12, int runN, int eventN); // event dump.

		void DoCommonStuff(CommonStuff &miscstuff); // reads raw Met, pho, ele info and fills CommonStuff  
		void DoMyMet(CommonStuff miscstuff, JetStuff &jetstuff); // corrects Met & SumEt for jets 
		void DoMyJet(TStnJetBlock* fJetBlock, CommonStuff miscstuff, JetStuff &jetstuff, MatchStudyHisto_t& Hist); // does my jets
		void DoMyJetNoMatch(TStnJetBlock* fJetBlock, JetStuff &jetstuff); // does my jets, but doesn't remove EM objetcs
		void DoMyJetWithMatch(TStnJetBlock* fJetBlock, CommonStuff miscstuff, 
				JetStuff &jetstuff, MatchStudyHisto_t& Hist); // does my jets after removing EM objects
		void ReorderMyJets(JetStuff &jetstuff); // reorders jets after removing EM objects
		void ReorderSmearedJets(std::vector<TLorentzVector>&); // reorders smeared jets (newLev6 smeared jets) in pseudo experiments // sam 1/2/2009

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

		double MyTotalMetSignificance(JetStuff jetstuff, TVector2 myMetVec, int systcode, int systcodeUncl); // My Met significance for all events
		double MyInCemTowerEta(double eta_det); // returns relative position of EM shower inside CEM tower
		double MyIntegralUpLimit(double jet_E, double eta_det, int stat_code, int syst_code);
		void ClearMetProbStuff(MetProbStuff& stuff); // clears params for MetProb
		void InitMetProbStuff(MetProbStuff& stuff,JetStuff jetstuff, TVector2 myMetVec, int systcode, int systcodeUncl); // inits params for MetProb;
		double MetSigCorrPAR(int ind, int njet15, int isMC, int sample_id); // returns input parameters for MetSig correction
		double MetSigCorrection(double rawSig); // returns corrected MetSig

		void CalculateMetDelPhi(MetSigDelPhiStuff &msdp_stuff,JetStuff jetstuff,CommonStuff miscstuff,TVector2 MetVec); // calculates dPhi for metsig.vs.dPhi studies

		void FinalMetSigCorrHisto(MetStudyHisto_t& Hist); // filling final histogram for MetSig correction 

		void GenerateMyTotalMet(JetStuff& jetstuff, CommonStuff miscstuff, MetCleanupHisto_t& Hist, 
				int systcode, int systcodeUncl, int metcode, int rnd_seed, TVector2 &myGenMet); 
		void GenerateMyEMobjMet(CommonStuff miscstuff, 
				int systcode, TVector2 &myEMGenMet); // new function (09/20/06), MET due to EM obj. resolution
		void MyNewJetWithMet(JetStuff &jetstuff, const CommonStuff&); // adding jet and Met
		//______________ counts generated events above Met cut
		void MyMetModelEventCount(JetStuff jetstuff, TVector2 myGenMet_def, MetResults &metstuff);
		//______________ counts data events above Met cut
		void MyMetEventCount(JetStuff jetstuff, CommonStuff miscstuff, int metcode, MetResults &metstuff);
		//______________ prints out a comparison of Met in Data to Met Model predictions
		void MyMetResults(MetResults &metstuff);

		void CleanMetResults(MetResults &metstuff); // cleans met_results vectors

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
				TLorentzVector *mypho, TLorentzVector *myjet, 
				double &dR_fj, double &dPhi_fj, double &dEta_fj);
		//________________________________ returns Qt of all jets above Eta & Et thresholds
		double MyQtJet(JetStuff jetstuff);
		//________________________________ returns SumEt=jets+photons+leptons+taus+met (jets above Eta & Et thresholds)
		//                                 metscenario= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
		double MyHtAll(std::vector<TLorentzVector> vJets, JetStuff jetstuff, CommonStuff miscstuff, int metscenario);
		double SmearedHtAll(const JetStuff& jetstuff, const CommonStuff& miscstuff, const int metscenario);

		//_________________________________ returns normalized kT_perp. for resolution studies 
		float GetMyKtPerp(std::vector<TLorentzVector> vec);  
		//_________________________________ returns normalized kT_parl. for resolution studies 
		float GetMyKtParl(std::vector<TLorentzVector> vec);  

		/// sam's stuff
		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int evtsMoreHad;
		};

		int GetNMatchedObj() const { return sam_matchstuff.size(); }
		int GetMatch_JetIndex(int ObjType, int EmInd) const;//Returns original jet block index of the jet matched to an EM object

		TVector2 GetGenMet_Def() const { return samGenV2Met_def; }
		TVector2 GetGenMet_devm() const { return samGenV2Met_devm; }
		TVector2 GetGenMet_devp() const { return samGenV2Met_devp; }
		TVector2 GetGenMet_devmUn() const { return samGenV2Met_devmUn; }
		TVector2 GetGenMet_devpUn() const { return samGenV2Met_devpUn; }

		void SetDebug(const int num) { debug = num; }
		int Debug() const { return debug; }
		void DebugInfo(const int iAnaMode, const std::string sFunction, const int printlvl, const std::string message);



		void GetGenLevelInfo(); // to study generator level met stuff	
		void SetUseHadJets(const bool b) { bUseHadJets = b; }
		bool GetUseHadJets() const { return bUseHadJets; }
		void DumpJets(const std::string func_name, const unsigned int line_num,
							const JetStuff&, const int jetLev);
		void DumpSmearedJets(const std::string func_name, const unsigned int line_num,
							const JetStuff&, const int jetLev);
		void DumpEMobjects(const std::string func, const unsigned int line, const CommonStuff&, const int type=0);
		
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
		void JetMismeasureProb(JetStuff&, const TVector2& MetVec, const int systcode=0);
		void Dump_newJetLev6(const std::string funcname, const int line, const JetStuff& jetstuff); // Dumps newJetLev6 jet info
		void SetDoubleSmearJets(const bool b) { bDoubleSmearJets = b; }
		bool DoubleSmearJets() const { return bDoubleSmearJets; }
		void SetZeroJERoffset(const bool b) { bOffsetJER = b; }
		bool ZeroJERoffset() const { return bOffsetJER; }
		int LeastMismeasuredJet(const std::vector<TLorentzVector> vJets, const TVector2 MetVec);
		void SetMetCut(float fM) { fMetCut = fM; }
		float MetCut() const { return fMetCut; }
		TLorentzVector FindMatchingHEPGPar(const TLorentzVector tlObj,
									const int iPDGcode, const int iStatus, const float fDelR);
		TLorentzVector FindMatchingHADJet(const TLorentzVector tlJetVec,
				const float fDelR);
		void DumpTLvec(const std::string func, const unsigned int line,
								const std::vector<TLorentzVector> vObj , const std::string label="");
		void SetPrintLevel(const int p) { iPrintLvl = p; }
		int PrintLevel() const { return iPrintLvl; }
		int JetClosestToMet(std::vector<TLorentzVector> vJet, const TVector2 Met);
		void SetMEtAddMethod(const int i) { iMetAddMethod = i; }
		int MEtAddMethod() const { return iMetAddMethod; }

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

		void UseOldJERForGaussMean(const bool t)
		{ 
			if (t) abUseJERnew.GaussMean_New = JER_OLDFITS;
			else abUseJERnew.GaussMean_New = JER_GAUSSMEAN;
		}
		void UseOldJERForGaussSigma(const bool t)
		{ 
			if (t) abUseJERnew.GaussSigma_New = JER_OLDFITS;
			else abUseJERnew.GaussSigma_New = JER_GAUSSSIGMA;
		}
		void UseOldJERForLandauMPV(const bool t) 
		{ 
			if (t) abUseJERnew.LandauMPV_New = JER_OLDFITS;
			else abUseJERnew.LandauMPV_New = JER_LANDAUMPV;
		}
		void UseOldJERForLandauSigma(const bool t)
		{ 
			if (t) abUseJERnew.LandauSigma_New = JER_OLDFITS; 
			else abUseJERnew.LandauSigma_New = JER_LANDAUSIGMA; 
		}
		void UseOldJERForGaussNorm(const bool t) 
		{ 
			if (t) abUseJERnew.GaussNorm_New = JER_OLDFITS; 
			else abUseJERnew.GaussNorm_New = JER_GAUSSNORM; 
		}
		unsigned int JERForGaussMean() const { return abUseJERnew.GaussMean_New; }
		unsigned int JERForGaussSigma() const { return abUseJERnew.GaussSigma_New; }
		unsigned int JERForLandauMPV() const { return abUseJERnew.LandauMPV_New; }
		unsigned int JERForLandauSigma() const { return abUseJERnew.LandauSigma_New; }
		unsigned int JERForGaussNorm() const { return abUseJERnew.GaussNorm_New; }

		
		//thes will hold the JET fit paramaters during the job
		//vJerParam, vJerParamErr are default JER parameters and corresponding errors
		//vJerParamSysy and vJerParamSystErr are JER parameters for systematic
		//calculation for the JER fits shapes
		static std::vector< std::vector<double> > vJerParam, vJerParamErr, vJerSystParam, vJerSystParamErr;
		//3D vector for exact bin values <60GeV for old fit
//		static std::vector<std::vector<std::vector<double> > > vJerBinValue;


		
		
	private:
		Global_Counters_t mycounter;
		bool bRunPermit;		//__ make sure we have all dependencies met before running
		InitSuperPhotons *initSpMod;
		TriggerModule *trigMod;
		TagTightPhotons *tightMod;
		TagLoosePhotons *looseMod;
		TagElectrons *tightEleMod;
		TagLooseElectrons *looseEleMod;

		TVector2 samGenV2Met_def;
		TVector2 samGenV2Met_devm;
		TVector2 samGenV2Met_devp;
		TVector2 samGenV2Met_devmUn;
		TVector2 samGenV2Met_devpUn;

		std::vector<MatchStuff> sam_matchstuff;  // to get the correct match jet ind. 

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		int debug;								// control the print statements
		bool bSaveThisEvent;			// to skim few stntuples for debugging
		bool bUseHadJets;				// 1=replaces Detector jets with matching HAD jets. 0=does nothing
		bool bNoSummary;				// control the end job summary print
		bool bDoubleSmearJets;		// smear jets an additional time for the pseudo experiments
		bool bOffsetJER;				// shifts the gausian mean position of the JER
		float fMetCut;					//Required minimum hard cut on Met
		int iPrintLvl;					//controls the print statements	
		int iMetAddMethod;			// decide how we add MEt to a jet

		//a very small eta value to remove the warning from Root when Pt=0 and called for Eta()
		float fMinisculeEta;
		
		// decides what types of EM objects are removed from the jet list.
		bool bRem_TightPho, bRem_LoosePho, bRem_SidebandPho, bRem_TightPhoLikeEle;
		bool bRem_LoosePhoLikeEle, bRem_StdLooseEle;

		JERFuncOldNew_t abUseJERnew; 	//switch JER to use old/new fit fucntions
												//default will be old

		ClassDef(JetFilterModuleV2,0)
};

#endif
