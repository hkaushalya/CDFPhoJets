#ifndef TGAMMAJETINIT_HH
#define TGAMMAJETINIT_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

//#include <iostream>
//#include <fstream>

#include <vector>
#include "TTree.h"
#include "TGraph.h"
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

#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/ana/TSuperPhoton.hh"
#include "rlc/Pho/TGoodRun.hh"
#include "Stntuple/photon/TPhoTrigBits.hh"

#include "Stntuple/obj/TStnEvent.hh"

//for PMT spike removal
#include "Stntuple/obj/TCalDataBlock.hh"
#include "CalorGeometry/CalConstants.hh"

#endif

class TGammaJetsInit: public TStnModule {
//--------------------------------------------------------------

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
		TGenpBlock* 		  fGenpBlock;
		TCalDataBlock*      fCalData;
		typedef std::vector<TCalTower*> CalDataArray;  // need for beam halo
		typedef std::vector<TCalTower*>::iterator CalDataArrayI;   // need for BH

	public:
		TGammaJetsInit(const char* name="GammaJetsInit", const char* title="GammaJetsInit");
		~TGammaJetsInit();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TCprDataBlock*     GetCprDataBlock() { return fCprDataBlock; }
		TStnJetBlock*      GetJetBlock()     { return fJetBlock;     }
		TStnTrackBlock*    GetTrackBlock()   { return fTrackBlock;   }
		TStnPhotonBlock*   GetPhotonBlock()  { return fPhotonBlock;  }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnTriggerBlock*  GetTriggerBlock() { return fTriggerBlock; }
		TStnVertexBlock*   GetVertexBlock()  { return fVertexBlock;  }
		TGenpBlock* 		 GetGenpBlock() 	 { return fGenpBlock;    }


				/************** abstract objects *************************/
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
			bool BeamHalo;				// 0 = no 1=yes
			bool Cosmic;				// 0 = no 1= yes
			bool Conversion;			// 0= no 1= yes
			float HadTDCtime;				// timinf of Had TDCs
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
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassGoodrun;
			unsigned int evtsPassTrigger;
			unsigned int evtsPass25Trigger;
			unsigned int evtsPass50Trigger;
			unsigned int evtsPass70Trigger;
			unsigned int evtsPassVertexCut;
			unsigned int evtsPassZvCut;
			unsigned int evtsPassModule;	//___number of evts pass this module

			unsigned int evtsPassLooseIDCuts;
			unsigned int evtsPassTightIDCuts;
         unsigned int evtsPassPhoLikeIDCuts;

			unsigned int phos_match2_genp;
			unsigned int phos_nomatch2_genp;

			unsigned int found_hepg_Ws;

			unsigned int photons;		//__ all photons passing the photon id cuts
			unsigned int electrons;		//__ all photons passing electron id cuts(i.e. photon like electron id cuts)
			unsigned int none;			//__ all photons that did not pass either of the above (what are these)
			unsigned int HepgMatchFailFiducial;	//__

			//W ele id cut eff
			unsigned int ideff_hepgWs;				//__ hepg W events
			unsigned int ideff_hepgWenus;			//__ W that decayed to electron+ nutrino
			unsigned int ideff_hepgWs_passFid;	//__ Wenu, electron found in central region
			unsigned int ideff_hepgWs_passEt30;	//__ that electron having Et > 30
			unsigned int ideff_passGoodrun;		//__ from the events with a W, number of events pass good run
			unsigned int ideff_passVertex;		//__ then pass 1 good vertex
			unsigned int ideff_passZv;				//__ then pass Zv cut
			unsigned int ideff_passMetcut;		//___ need met>20
			unsigned int ideff_match2OfflineObject;	//__ then findind a matching offline photon in delR space
			unsigned int ideff_detWs;						//__ then that matching photon passing electron id cuts

			//Z ele id cut eff
			unsigned int zideff_hepgZs;				//__ hepg Z events
			unsigned int zideff_hepgZees;				//__ Z that decayed to electron+ nutrino
			unsigned int zideff_hepgZs_passFid;		//__ ZeE, electrons found in central region
			unsigned int zideff_hepgZs_passEt30;	//__ that electrons having Et > 30
			unsigned int zideff_passGoodrun;			//__ from the events with a Z, number of events pass good run
			unsigned int zideff_passVertex;			//__ then pass 1 good vertex
			unsigned int zideff_preZs;					// pre tag zee
			unsigned int zideff_passZv;					//__ then pass Zv cut
			unsigned int zideff_match2OfflineObject;	//__ then findind a matching offline photon in delR space
			unsigned int zideff_passConv;						//__matching photon are not conversion
			unsigned int zideff_detZs;						//__ then that matching photon passing electron id cuts
			//data W/Z estimate
			unsigned int dWe_detWs;
			unsigned int dZe_detZs;

		};


		// all the histograms

		struct GeneralPlots_t {
			TH1F* photon_loose_and_tight_cuts_cumm;
			TH1F* Nvertices;
			TH1F* vertexZ;
			TH1F* pho_passIDs;
		}; //GeneralPlots_t


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

		struct PhotonPlots_t {
			TH1F* Detector;
			TH1F* EtCorr;
			TH1F* XCes;
			TH1F* ZCes;
			TH1F* HadEm;
			TH1F* IsoEtCorr;
			TH1F* Chi2Mean;
			TH1F* N3d;
			TH1F* TrkPt;
			TH1F* TrkIso;
			TH1F* Ces2Wire;
			TH1F* Ces2Strip;
			TH1F* HadTDCtime;
		};

		struct EleIDeffPlots_t {	// W electron id eff plots
			TH1F* MinDelRMatch;
			TH1F* Et_ratio;
			TH1F* HepgEleEt;
			TH1F* OfflineEleEt;
			TH1F* OfflineEleTrkpt;
			TH1F* Wele_IDeff;
			TH1F* CEMPHOEleIDeff;
		};

		struct ZEleIDeffPlots_t {	// Z electrons id eff plots
			TH1F* ele1Et_ratio;
			TH1F* ele2Et_ratio;
			TH1F* hepgElesPhi;
			TH1F* detElesDelPhi;
			TH1F* Zeles_IDeff;
			TH1F* hepgZpt;
			TH1F* hepgZmass;
		};

		struct HEPGMatchPlots_t {	// offline objects  matched to HEPG level
			TH1F* Et_ratio;				// et of  matched objects
		};

		struct BackgroundPlots_t {	//__ all background plots

			// objects passing photon-like electron id cuts
			TH1F* Wele_EtCorr;
			TH1F* Wele_HadEm;
			TH1F* Wele_IsoEtCorr;
			TH1F* Wele_TrkPt;
			TH1F* W_eleEtMetratio;
			TH1F* W_Sumet0;
			TH1F* W_Mass;
			TH1F* W_Pt;
			TH1F* W_Met;

			TH1F* Zele1_EtCorr;
			TH1F* Zele1_HadEm;
			TH1F* Zele1_IsoEtCorr;
			TH1F* Zele1_TrkPt;
			TH1F* Zele2_EtCorr;
			TH1F* Zele2_HadEm;
			TH1F* Zele2_IsoEtCorr;
			TH1F* Zele2_TrkPt;
			TH1F* Z_Mass;
			TH1F* Z_Pt;
			TH1F* Z_eleEtratio;
			TH1F* Z_Sumet0;
			TH1F* Z_Met;

			TGraph* RunAvgWMet;
			TGraph* RunAvgZMet;
			TGraph* RunEvents;
			TGraph* RunAvgWs;
			TGraph* RunAvgZs;
		};


		struct RunAvg_t {
			int runNumber;
			unsigned int events;		// total number of events in the run
			long double sumWMet;				// sum of Met of W events
			long double sumZMet;				// sum of Met of Z events
			unsigned int Wcount;			//sum of all W events upto and including the current run
			unsigned int Zcount;			//sum of all Z events upto and including the current run
		};

		// for PMT spike removal
		struct CalHist_t {
			//all histograms
		};
		struct BadTower {};

		// end PMT spike removal



				/**********end abstract objects *************************/



		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void 		BookHistograms();
		void 		SaveHistograms();
		void  	Cleanup();
		bool 		PrePass();
		bool 		PassGoodRun();
		bool 		PassTrigger();
		bool 		PassZvCut();
		int 		Class12Vertices();
		double 	BestVertex_Z();

		void 		TagConversions();
		double 	ConversionSep(TStnTrack* t0, TStnTrack* t1, double & rconv);
		bool 		IsConversionElectron(TStnElectron* ele, 
				     double& minsep, double& mindt, double& radius);
		int  		MatchPhotonToHEPGLevel(TStnPhoton*); // return mom_pdg code
		bool 		FindZ(int&, int&);			// look for Z and return the vectors of two electrons
		bool 		FindW(int&);			// look for a W and return the vector of the electron
		void 		WEleIDcutEff();
		void 		ZEleIDcutEff();
		void 		W_mcEstimate();
		void 		W_dataEstimate();
		void     FillWplots(int);
		void     FillZplots(int,int);

		int 		FiducialLepton(TLorentzVector, double Zv);
		void 		SetPhotonList();

		unsigned int CEMPhoEleIDcut(TStnPhoton*, TStnEvent*, double Zvtx);
		unsigned int TightPhotonIDcut(TStnPhoton*);
		unsigned int LoosePhotonIDcut(TStnPhoton*);

		TFolder* GetHistoFolder(char* name, char* title);
		void 		BookBackgroundHistograms();
		void 		BookGeneralHistograms();
		void 		BookMetHistograms();
		void 		BookPhotonHistograms();
		void 		BookElectronHistograms();
		void 		BookHEPGMatchingHistograms();
		void     BookWeleIDeffHistograms();
		void     BookZeleIDeffHistograms();

		void 		FillDataBlocks(int ientry);  // put all GetEntry statements here
		void 		FillWmassPlots(SuperPhoton_t&);
		void 		FillPhotonPlots(SuperPhoton_t&);
		void 		FillElectronPlots(SuperPhoton_t&);
		void 		FillMetPlots();
		void	   ElectronPhotonPlots(SuperPhoton_t & , int);
	

		//methods for Sasha's code
		void  	InitForSashaCode(int);




		//_________________ accesors for sasha's code
	  	int GetmyNpho() { 
			std::cout << "TGI: GetmyNpho()" << std::endl;
			return Nphotons;
		}
		
		TLorentzVector* GetmyUncorrPho(int i) {     
		  if(i>=0 && i<phoUncorrected.size()) return &(phoUncorrected[i]);
		  else return NULL;	  
		}
		
		TLorentzVector* GetmyCorrPho(int i) { 
		  if(i>=0 && i<phoCorrected.size()) return &(phoCorrected[i]);
		  else return NULL;	  
		}

		int GetmyPhoInd(int i) { return phoIndex[i]; } 
		double GetphoHadEm3(int i) { return phoHadEm[i]; }  // photon had em 3 TOWER
		double GetphoEtaDet(int i) { return phoEtaDet[i]; }  // photon DetEta
		double GetphoCesX(int i) { return phoCesX[i]; }
		double GetphoCesZ(int i) { return phoCesZ[i]; }
		//----- end sasha's stuff ------
	
		// accessors for my 2nd mod
		int GetSuperPhotonSize() { return superpho.size(); }
		int GetSuperPhotonTightID(int i) {
			if ( (i >=0) && (i<superpho.size())) { return superpho[i].TightID; }
		}
		int GetSuperPhotonPhoLikeID(int i) {
			if ( (i >=0) && (i<superpho.size())) { return superpho[i].PhotonLikeID; }
		}
		TLorentzVector GetSuperPhotonCorVec(int i) {
			if ( (i >=0) && (i<superpho.size())) { return superpho[i].corvec; }
		}
		TLorentzVector GetSuperPhotonRawVec(int i) {
			if ( (i >=0) && (i<superpho.size())) { return superpho[i].rawvec; }
		}

		//stuff for PMT SPIKE REMOVAL CODE

	  TCalDataBlock*      GetCalDataBlock()  { return fCalData; }
	  void SetUseFindSpike  (int cut)  { fUseFindSpike = cut; }
	  void SetUseSpikeFilter(int cut)  { fUseSpikeFilter = cut; }
	  void SetMinEmPMTasym(double cut)  { fMinEmPMTasym = cut; }
	  void SetMinHadPMTasym(double cut) { fMinHadPMTasym = cut; }
	  void SetMaxEmPMTasym(double cut)  { fMaxEmPMTasym = cut; }
	  void SetMaxHadPMTasym(double cut) { fMaxHadPMTasym = cut; }
	  void SetMinEmPmtE(double cut)     { fMinEmPmtE = cut; }
	  void SetMinHadPmtE(double cut)    { fMinHadPmtE = cut; }

	  void MatchCalorTowers(TStnJet* jet, TStnJetBlock* fJetBlock,
				              TCalDataBlock *fCalDataBlock,CalDataArray* cdholder);
	  void MatchCalorTowers(int pho_ind, TStnPhotonBlock* fPhotonBlock,
									TCalDataBlock *fCalDataBlock,CalDataArray* cdholder);

	  int FindPmtSpikes(CalDataArray* towers, double& EmESpike, double& HadESpike, double& spikeHadTdc, CalHist_t& Hist);
	  void DoMyCalStuff(TStnJetBlock* fJetBlock,TStnPhotonBlock* fPhotonBlock,TCalDataBlock *fCalDataBlock,CalHist_t& Hist);
	  static int towerInd(int iEta) { return  iEta >= TENETA/2 ? iEta - TENETA/2 : -(TENETA/2 - 1 - iEta); }

		// END STUFF FOR PMT SPIKE REMOVAL CODE



	private:
		bool 	qMc;  				// true if MonteCarlo data
		bool 	trig25, trig50, trig70;
		Global_Counters_t counter;
		
		TGoodRun 		goodrun;
		TStnDBManager 	DBman;
		TString 			goodRunListFile;
		TPhoTrigBits 	trigbits;

		//___________ for sasha's code
	  	int Nphotons;
		std::vector<TLorentzVector> phoUncorrected, phoCorrected;
		std::vector<float> phoHadEm;
		std::vector<float> phoEtaDet, phoCesX, phoCesZ;
		std::vector<int> phoIndex;
		
		//my super photons vector
		std::vector<SuperPhoton_t> superpho;

		// histograms
		GeneralPlots_t GeneralPlot;
		MetPlots_t MetPlot;
		PhotonPlots_t PhotonPlot,ElectronPlot;
		HEPGMatchPlots_t HEPGMatchPlot;
		EleIDeffPlots_t EleIDeffPlot;
		ZEleIDeffPlots_t ZEleIDeffPlot;
		BackgroundPlots_t BackgroundPlot;

		//for run dependent plots
		double sumWeventMet, sumZeventMet;
		int currRun, prevRun;
		std::vector<RunAvg_t> RunAvg;		
		
		bool FoundGoodPhoton;		// if there is a good photon, pass this module
	

		// stuff for PMT spike removal
		int fUseFindSpike;		// status code for PMT spike studies
		int fUseSpikeFilter; 	// status code for PMT spike filter

		double fMinEmPMTasym;  	//minimum EM PMT asym
		double fMinHadPMTasym;	// minimum HAD PMT asym
		double fMaxEmPMTasym;	// maximum EM PMT asym
		double fMaxHadPMTasym;	// maximum HAD PMT asym

		double fMinEmPmtE;		//minimum EM PMT energy
		double fMinHadPmtE;		//minimum HAD PMT energy

	  std::vector<BadTower> badEmTwr;  // list of bad EM towers
	  std::vector<BadTower> badHadTwr; // list of bad HAD towers

	  int Nspike_global;   // global number of spikes; to be used in cut


	
	ClassDef(TGammaJetsInit,1)
};

#endif
