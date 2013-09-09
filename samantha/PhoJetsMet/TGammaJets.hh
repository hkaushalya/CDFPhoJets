#ifndef TGAMMAJET_HH
#define TGAMMAJET_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

//#include <iostream>
//#include <fstream>

#include <vector>
#include "TTree.h"
#include <Stntuple/loop/TStnModule.hh>

#include <Stntuple/obj/TStnHeaderBlock.hh>
#include <Stntuple/obj/TCprDataBlock.hh>
#include <Stntuple/obj/TStnJetBlock.hh>
#include "Stntuple/obj/TStnTrackBlock.hh"
#include <Stntuple/obj/TStnPhotonBlock.hh>
#include <Stntuple/obj/TStnTriggerBlock.hh>
#include "Stntuple/obj/TStnElectronBlock.hh"
#include "Stntuple/obj/TStnVertexBlock.hh"
#include "Stntuple/obj/TStnMetBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"

#include <Stntuple/obj/TStnDBManager.hh>
#include "Stntuple/obj/TL3SummaryBlock.hh"

#include "Stntuple/alg/TStntuple.hh"
#include <Stntuple/alg/TStnElectronID.hh>
#include <Stntuple/alg/TStnPhotonID.hh>

#include "Stntuple/obj/TStnEvent.hh"

//for jay
#include "Stntuple/obj/TCalDataBlock.hh"

class TGammaJetsInit;
#endif


class TGammaJets: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TCprDataBlock*      fCprDataBlock;
		TStnJetBlock*       fJetBlock;
		TStnTrackBlock*     fTrackBlock;
		TStnPhotonBlock*    fPhotonBlock;
		TStnElectronBlock*  fElectronBlock;
		TStnTriggerBlock*   fTriggerBlock;
		TStnVertexBlock*    fVertexBlock;
		TStnMetBlock*       fMetBlock;
		TL3SummaryBlock*    fL3SummaryBlock;
		TGenpBlock*			  fGenpBlock;
		TCalDataBlock*      fCalDataBlock;
		
	public:
		TGammaJets(const char* name="GammaJets", const char* title="GammaJets");
		~TGammaJets();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TCprDataBlock*     GetCprDataBlock() { return fCprDataBlock; }
		TStnJetBlock*      GetJetBlock()     { return fJetBlock;     }
		TStnTrackBlock*    GetTrackBlock()   { return fTrackBlock;   }
		TStnPhotonBlock*   GetPhotonBlock()  { return fPhotonBlock;  }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnTriggerBlock*  GetTriggerBlock() { return fTriggerBlock; }
		TStnVertexBlock*   GetVertexBlock()  { return fVertexBlock;  }
		TL3SummaryBlock* GetL3SummaryBlock() { return fL3SummaryBlock;}

		// ****** setters
		void SetPhotonList(int ientry);
		void SetJetList(int ientry);
		
		// ****** getters
		unsigned int GetLoosePhotonIDWord(int i) {
			if ((i>=0) && (i<superpho_vec.size())) return superpho_vec[i].LooseID;
		}
		unsigned int GetTightPhotonIDWord(int i) {
			if ((i>=0) && (i<superpho_vec.size())) return superpho_vec[i].TightID;
		}

		std::vector<TStnPhoton*> GetSortedPhotonList();		//the sorted photon list from superpho_vec
		std::vector<TStnJet*> GetSortedJetList();				//the sorted jet list from the superjet_vec

		TStnPhoton* GetSuperPhoton(int i) {	//______________ get the i th super photon			
			if (i>=0 && i < superpho_vec.size() ) return superpho_vec[i].pho;
		}
		TStnJet* GetSuperJet(int i) {		//_________________ get the i th super jet
			if (i>=0 && i < superjet_vec.size() ) return superjet_vec[i].jet;
		}
		// --------- end getters -----------------------

		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
		void FillDataBlocks(int);  // put all GetEntry statements here
		void GammaJetsSelection(int);
		unsigned int LoosePhotonCuts(TStnPhoton*);
		unsigned int TightPhotonCuts(TStnPhoton*);
		double BestVertex_Z();
		int    Class12Vertices();
		void   JetMassPlots(int ientry);
		void   TriggerEfficiency(int ientry);
		void RemovePhotonsFromJetBlock();

		void  L3Summary(int);
		double Get3BodyMass(TLorentzVector,TLorentzVector,TLorentzVector);
		void 	RemoveOverlap(int ientry);
		void  JetPlots(TStnJet*, TStnJet*);
		void  PlotJetPhotonSeperation(int);
		void  PlotPhotonJetsPhi(int);
		void  MissingEtPlots();
		void  GenpBlockStudy(int);
		void  PhotonIDPlots(int);
		void  PhotonCutsCumulativePlot();
		void  PhotonIDcutsAna(int);
		double GetInvariantMass(TLorentzVector*,TLorentzVector*);
		double GetInvariantMass(TLorentzVector*,TLorentzVector*,TLorentzVector*);
		unsigned int MatchPho2GenLevel(TStnPhoton*);
		void  DoCommonThings(int ientry);

			// background study methods
		void  Background_W();
		void  Background_Z();
		void  Background_ElectronFakeRate(int ientry);

			//fake rate functions
		int   CEMPhoEleIDcut(TStnPhoton*, TStnEvent*, double);
		enum { FkRtParArray = 4 };
		enum { PhnxParArray = 3 };

		double fakerate_par[FkRtParArray]; // fake rate parametrization
		double fakerate_parerr[FkRtParArray]; // fake rate parametrization uncertainty
		double fakerate_corr_coeff[FkRtParArray][FkRtParArray]; // covariance matrix coefficients
		double phnx_par[PhnxParArray]; // phnx eff. parametrization
		double phnx_parerr[PhnxParArray]; // phnx eff. parametrization uncertainty
		double phnx_corr_coeff[PhnxParArray][PhnxParArray]; // covariance matrix coefficients

		void FakeRate(double eleEt, int goodrun_bit, int phnx_bit, 
			      double &wght_d, double &wght_m, double &wght_p); // returns fake rate & syst
		double ShiftCorrFunc(double eleEt); // correction for Et(det)/Et(genp) difference
		double PhnxEffFunc(double eleEt);  // phoenix efficiency function
		double PhnxEffFuncErr(double eleEt);  // uncertainty on phoenix efficiency 
		double FakeRateFunc(double eleEt); // returns fake rate
		double FakeRateFuncErr(double eleEt); // returns uncertainty on fake rate

		void FindConversion();// get mathing electron to check if its a conversion
		bool IsConversionElectron(TStnElectron* ele, double& minsep, double& mindt, double& radius);
		double ConversionSep(TStnTrack* t0, TStnTrack* t1, double& rconv);

		void JayTower(int ientry);
	
	
		//methods for Sasha's code
		bool  InitForSashaCode(int);
		void  Cleanup();
 		bool IsWhatIWant();
 
		// for making strip ntuples
		void StripCuts(int ientry);
 
				/****************** struct's FOR HISTROGRAMS **************/
		// structs for ALL THE HISTOS
		// convention: j1 = Leading Jet, j2 = 2nd Leading Jet
		// photon, is always the leading photon.

/*{{{*/		

		struct MassPlots_t {
			TH1F* photon_j1_j2;
			TH1F* photon_j1_j2_fake;
			TH1F* photon_j1;
			TH1F* photon_j2;
			TH1F* j1_j2;
		}; //MassPlots_t

		struct GeneralPlots_t {
			TH1F* photon_loose_and_tight_cuts_cumm;
			TH1F* photon_loose_and_tight_cuts_cumm_temp;
			TH1F* j1_to_photon_energy_ratio;
			TH1F* j2_to_photon_energy_ratio;
			TH1F* j2_to_j1_energy_ratio;
			TH1F* j1_momentum_to_energy_ratio;
			TH1F* j2_momentum_to_energy_ratio;
			TH1F* j1_to_photon_momentum_ratio;
			TH1F* j2_to_photon_momentum_ratio;
			TH1F* j2_to_j1_momentum_ratio;
			TH1F* j1_photon_phi_seperation;
			TH1F* j2_photon_phi_seperation;
			TH1F* j2_j1_phi_seperation;
			TH1F* j2_and_j1_plus_photon_phi_seperation;  // j2 w.r.t to sum vector of j1 and photon
			TH1F* Nvertices;
			TH1F* vertexZ;
			TH1F* pho_IDbins_pass;
			TH1F* pho_IDbins_fail;
		}; //GeneralPlots_t

		struct PhotonPlots_t {
			TH1F* Detector_b;
			TH1F* EtCorr_b;
			TH1F* XCes_b;
			TH1F* ZCes_b;
			TH1F* HadEm_b;
			TH1F* IsoEtCorr_b;
			TH1F* Chi2Mean_b;
			TH1F* N3d_b;
			TH1F* TrkPt_b;
			TH1F* TrkIso_b;
			TH1F* Ces2Wire_b;
			TH1F* Ces2Strip_b;

			TH1F* Detector_aL;
			TH1F* EtCorr_aL;		// AFTER LOOSE CUTS
			TH1F* XCes_aL;
			TH1F* ZCes_aL;
			TH1F* HadEm_aL;
			TH1F* IsoEtCorr_aL;
			TH1F* Chi2Mean_aL;
			TH1F* N3d_aL;
			TH1F* TrkPt_aL;
			TH1F* TrkIso_aL;
			TH1F* Ces2Wire_aL;
			TH1F* Ces2Strip_aL;

			TH1F* Detector_aT;
			TH1F* EtCorr_aT;		//AFTER TIGHT CUTS
			TH1F* XCes_aT;
			TH1F* ZCes_aT;
			TH1F* HadEm_aT;
			TH1F* IsoEtCorr_aT;
			TH1F* Chi2Mean_aT;
			TH1F* N3d_aT;
			TH1F* TrkPt_aT;
			TH1F* TrkIso_aT;
			TH1F* Ces2Wire_aT;
			TH1F* Ces2Strip_aT;
		}; //PhotonPlots_t

		struct JetPlots_t {
			TH1F* EtRatio;
			TH1F* DelPhi;
			TH1F* DelEta;
			TH1F* DelR;
		};
		struct Jet1Plots_t {
			TH1F* Etraw;
			TH1F* Phi;
			TH1F* DetEta;
			TH1F* Emfr;
			TH1F* EOverP;
		};
		struct Jet2Plots_t {
			TH1F* Etraw;
			TH1F* Phi;
			TH1F* DetEta;
			TH1F* Emfr;
			TH1F* EOverP;
		};

		// each sub mod should have their own set of plots
		// so things will be easy to swicth on and off
		struct L3SummaryPlots_t {
			TH1F* NEm;
			TH1F* NMuon;
			TH1F* NTau;
			TH1F* NJet4;
			TH1F* NJet7;
			TH1F* NMet;
			TH1F* NCotTrack;
			TH1F* NSvtTrack;
			TH1F* NTau_to_NEm;
			TH1F* NMuon_to_NEm;
			TH1F* NJet4_to_NJet7;
			TH1F* NJet7_to_NJet4;
			TH1F* NSvtTrack_to_NCotTrack;
			TH1F* NCotTrack_to_NSvtTrack;
			TH2F* NCotTrack_vs_NEm;
			TH2F* NSvtTrack_vs_NEm;

		}; //L3SummaryPlots_t

		struct GenpStudyPlots_t {
			TH1F* Photon_Et;	
			TH1F* Photon_Eta;	
			TH1F* Photon_Phi;	
			TH1F* JetsLostInTheDetector_Eta;
			TH1F* JetsLostInTheDetector_Phi;
			TH1F* PhotonsLostInTheDetector_Eta;
			TH1F* PhotonsLostInTheDetector_Phi;
		}; //GenpStudyPlots_t


		struct MetPlots_t {
			TH1F* Met_0;
			TH1F* MetPhi_0;
			TH1F* MetX_0;
			TH1F* MetY_0;
			TH1F* Sumet_0; // Z=0
			TH1F* Sumet_1;	// high pt lepton/ 1st vertex
			TH1F* Sumet_2;	//high pt vertex
			TH1F* Metsig;		//what????
			TH1F* Z0;
			
		}; //MetPlots_t

		struct PhoIDcutsPlots_t {			// for n-1 ID analyses
			TH1F* Detector_b;		//pass the trigger but before any other cut
			TH1F* EtCorr_b;
			TH1F* XCes_b;
			TH1F* ZCes_b;
			TH1F* HadEm_b;
			TH1F* IsoEtCorr_b;
			TH1F* TrkPt_b;
			TH1F* TrkIso_b;

			TH1F* Detector_a_n_1;		// after applying n-1 cuts
			TH1F* EtCorr_a_n_1;
			TH1F* XCes_a_n_1;
			TH1F* ZCes_a_n_1;
			TH1F* HadEm_a_n_1;
			TH1F* IsoEtCorr_a_n_1;
			TH1F* TrkPt_a_n_1;
			TH1F* TrkIso_a_n_1;

			TH1F* Detector_a;		// after only this cut is applied
			TH1F* EtCorr_a;
			TH1F* XCes_a;
			TH1F* ZCes_a;
			TH1F* HadEm_a;
			TH1F* IsoEtCorr_a;
			TH1F* TrkPt_a;
			TH1F* TrkIso_a;
		};

		struct HEPGPhoMatchPlots_t {	// when a offline photon is matched ti HEPG level
			TH1F* DelEta;
			TH1F* DelPhi;
			TH1F* DelR;
			TH1F* Nmatches;			// number of matches found per search
			TH1F* matched_Object;	//if matched photon mom is a Pi0/Eta/Omega...
			TH1F* HEPGObj_Et;				// et of all matched photons
			TH1F* HEPGObj_Eta;
			TH1F* Et_ratio;				// et of all matched photons
			
		};

		struct 	Em2JetMatchPlots_t {	// for removing EM OBJECT FROM jet block
			TH1F* DelEta;
			TH1F* DelPhi;
			TH1F* DelR;
			TH1F* NOmatches;			// number em objects/per event, that did not find match
			TH1F* Et_ratio;				// et ratio of matched em and jet
		};
		struct BackgroundPlots_t {	//__ all background plots
			// objects passing photon-like electron id cuts
			TH1F* Wele_EtCorr;
			TH1F* Wele_HadEm;
			TH1F* Wele_IsoEtCorr;
			TH1F* Wele_TrkPt;
			TH1F* W_Met;
			TH1F* W_eleEt_to_Met_ratio;
			TH1F* W_Sumet_0;
			TH1F* W_Sumet_1_Htc;
			TH1F* W_Sumet_2_Sumetjet;
			TH1F* W_Mass;
			TH1F* W_Pt;
			
			TH1F* Zele1_EtCorr;
			TH1F* Zele2_EtCorr;
			TH1F* Zele1_HadEm;
			TH1F* Zele1_IsoEtCorr;
			TH1F* Zele1_TrkPt;
			TH1F* Zele2_TrkPt;
			TH1F* Z_Met;
			TH1F* Z_Sumet_0;
			TH1F* Z_Sumet_1_Htc;
			TH1F* Z_Sumet_2_Sumetjet;
			TH1F* Z_Mass;
			TH1F* Z_Pt;
			TH1F* Z_eleDelPhi;
			TH1F* Z_eleDelEta;
			TH1F* Z_eleEtratio;

		};
/*}}}*/
				/**************END  struct's FOR HISTROGRAMS **************/

				/************** abstract objects *************************/
/*{{{*/
		struct SuperPhoton_t {		//adds more info to the photon
			unsigned int index;     //orginal order in the photon Block
			TStnPhoton*  pho;
			TLorentzVector corvec;
			TLorentzVector rawvec;
			int Detector;
			double DetEta;
			double Phi;
			double Etraw;
			double ECorr;
			double EtCorr;
			double XCes;
			double ZCes;
			double HadEm;
			double Chi2Mean;
			int N3d;
			double IsoEtCorr;
			double TrkPt;
			double TrkIso;
			double CesWireE2;
			double CesStripE2;
			unsigned int LooseID;
			unsigned int TightID;
			unsigned int PhotonLikeID;
			bool BeamHalo;
			bool Cosmic;
		};

		struct SuperJet_t {
			unsigned int index;
			TStnJet* jet;
			double Et;

		};

		struct HEPGmatchedObject_t {
			int index;	//_____________________ index of the l3 object that is macthed to hepg
			TObject* l3obj;	//________________ lets not use these for now. just worry about photons only
			TObject* hepgobj;
			double DelEta;
			double DelPhi;
			double DelR;
			double l3E_to_hepgE_ratio;		//energy ratio
		};
		

		struct LorentzVectors_t {
			TLorentzVector photon;
			TLorentzVector j1;
			TLorentzVector j2;
			TLorentzVector photon_j1_sum;
			TLorentzVector photon_j2_sum;
			TLorentzVector photon_j1_j2sum;
			TLorentzVector j1_j2_sum;
		};


		enum PhotonIDBits_t {  // bit i will be set when a cut is failed
			kCentral_L  	= 0x1 << 0,
			kEtCorr30_L		= 0x1 << 1,
			kXCes_L			= 0x1 << 2,
			kZCes_L			= 0x1 << 3,
			kHadEm_L 		= 0x1 << 4,
			kTrkPt_L			= 0x1 << 5,
			kTrkIso_L		= 0x1 << 6,
			kIsoEtCorr_L	= 0x1 << 7,
			kEtCorr7_T		= 0x1 << 8,
			kXCes_T			= 0x1 << 9,
			kZCes_T			= 0x1 <<10,
			kHadEm_T 		= 0x1 <<11,
			kChi2Mean_T		= 0x1 <<12,
			kN3d_T 			= 0x1 <<13,
			kTrkPt_T			= 0x1 <<14,
			kTrkIso_T		= 0x1 <<15,
			kIsoEtCorr_T	= 0x1 <<16,
			kCesWireE2_T 	= 0x1 <<17,
			kCesStripE2_T 	= 0x1 <<18,
			kCentral_T  	= 0x1 <<19
		};
		
		struct Global_Counters_t {
			unsigned int evtsWithoutAphoton;
			unsigned int evtsWithoutAjet;
			unsigned int evtsPassTrigger;
			unsigned int evtsPass25Trigger;
			unsigned int evtsPass50Trigger;
			unsigned int evtsPass70Trigger;
			unsigned int evtsWithPhotonAnd2Jets;
			unsigned int evtsPassLooseCuts;
			unsigned int evtsPassTightCuts;
			unsigned int evtsPass2JetsCut;
			unsigned int sidebandevts;			//number of side band events
			unsigned int total_evts_b4_strip;
			unsigned int total_evts_a4_strip;
			unsigned int phos_match2_genp;
			unsigned int phos_nomatch2_genp;
			unsigned int HEPGmomPi0;
			unsigned int HEPGmomElectron;
			unsigned int HEPGmomMuon;
			unsigned int HEPGmomTau;
			unsigned int HEPGmomChargedPi;
			unsigned int HEPGmomEta;
			unsigned int HEPGmomOmega;
			unsigned int HEPGmomOther;
			unsigned int Wevents;
			unsigned int Zevents;
			unsigned int pho_2j_events;
			unsigned int goodrunevts;
			
		};
/*}}}*/
				/**********end abstract objects *************************/

		// create folders for histograms
		TFolder* GetHistoFolder(char* name, char* title);

		//book my histograms
		void  BookInvariantMassHistograms();
		void  BookGeneralHistograms();
		void  BookL3SummaryHistograms();
		void  BookGenpStudyHistograms();
		void  BookPhotonHistograms();
		void  BookMetHistograms();
		void  BookBackgroundHistograms();
		void  BookPhotonIDcutsHistograms();
		void  BookHEPGPhotonMatchingHistograms();
		void  BookEmObjectRemovedFromJetsHistograms();
		void  BookJetHistograms();

		//_________________ accesors for sasha's code
	  	int GetmyNpho() { return superpho_vec.size(); }
		
		TLorentzVector* GetmyUncorrPho(int i) {     
		  if(i>=0 && i<superpho_vec.size()) return &(superpho_vec[i].rawvec);
		  else return NULL;	  
		}
		
		TLorentzVector* GetmyCorrPho(int i) { 
		  if(i>=0 && i<superpho_vec.size()) return &(superpho_vec[i].corvec);
		  else return NULL;	  
		}

		int GetmyPhoInd(int i) { return superpho_vec[i].index; } 
		double GetphoHadEm3(int i) { return superpho_vec[i].HadEm; }  // photon had em 3 TOWER
		double GetphoEtaDet(int i) { return superpho_vec[i].DetEta; }  // photon DetEta
		double GetphoCesX(int i) { return superpho_vec[i].XCes; }
		double GetphoCesZ(int i) { return superpho_vec[i].ZCes; }
		//----- end sasha's stuff ------
		
		TStnPhoton* GetLeadingPhoton() { if (superpho_vec.size() >0) return superpho_vec[0].pho; else NULL;}
		TStnPhoton* GetPhoton(int i) { if ( (i>=0) && (i<superpho_vec.size())) return superpho_vec[i].pho; else NULL;}
		TStnJet* GetLeadingJet1() { if (superjet_vec.size() >0) return superjet_vec[0].jet; else NULL;}
		TStnJet* GetLeadingJet2() { if (superjet_vec.size() >1) return superjet_vec[1].jet; else NULL;}
		TStnJet* GetJet(int i) { if ( (i>=0) && (i<superjet_vec.size()) ) return superjet_vec[i].jet; else NULL;}



	private:
		int 	run, evt, sec;
		bool 	qMc;  				// true if MonteCarlo data
		bool 	trig25, trig50, trig70;
		bool 	ignore_pho_zero_evts; // remind me that events without a photon/jet are removed or not
		Global_Counters_t counter;
		TStnPhoton* leading_pho;
		TStnJet *leading_jet1,*leading_jet2;
		TLorentzVector leadPhoVec, leadJet1Vec, leadJet2Vec;			// get these from TGammaJetsInit and TMyJetFilterModule
		//std::vector<TStnPhoton*> photon;										// store Et sorted photons
		//std::vector<TStnPhoton*> photon;										// store Et sorted photons
		//std::vector<TStnJet*> jet;
		TLorentzVector phovec,j1vec,j2vec, pj1sum, pj2sum, pj1j2sum, j1j2sum; //no need. i put the vectos in the super object

		MassPlots_t MassPlot;
		GeneralPlots_t GeneralPlot;
		PhotonPlots_t PhotonPlot;
		PhoIDcutsPlots_t PhoIDcutPlot;
		HEPGPhoMatchPlots_t HEPGPhoMatchPlot;	// when a offline photon is matched ti HEPG level
		Em2JetMatchPlots_t EmJetMatchPlot;
		JetPlots_t JetPlot;
		Jet1Plots_t Jet1Plot;
		Jet2Plots_t Jet2Plot;
		L3SummaryPlots_t L3SummaryPlot;
		GenpStudyPlots_t GenpPlot;
		LorentzVectors_t LorentzVector;
		MetPlots_t MetPlot;
		BackgroundPlots_t BackgroundPlot;


		//___________ for sasha's code
	  	//int Nphotons;
		//std::vector<TLorentzVector> phoUncorrected, phoCorrected;
		//std::vector<float> phoHadEm;
		//std::vector<float> phoEtaDet, phoCesX, phoCesZ;
		//std::vector<int> phoIndex;
		
		//my super photons vector
		std::vector<SuperPhoton_t> superpho_vec;
		std::vector<SuperJet_t> superjet_vec;


		//photon candidate index, so i don't have to call Init mod to get this again
		int leadPhoIndex;

		double fakeRateSum;			//sum of the fake rates

		std::vector<TLorentzVector> PhoVec;
		std::vector<int> PhoTightID, PhoEleID;
	ClassDef(TGammaJets,1)
};

#endif
