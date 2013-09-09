#ifndef TINIT_HH
#define TINIT_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <vector>
#include "TLorentzVector.h"

#include <Stntuple/loop/TStnModule.hh>

#include <Stntuple/obj/TStnHeaderBlock.hh>
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

//for BH/CR/PMT spike removal
#include "Stntuple/obj/TCalDataBlock.hh"
#include "CalorGeometry/CalConstants.hh"

// for CR
#include <Stntuple/obj/TStnEmTimingBlock.hh>
#include "Stntuple/ana/TEmTimingModule.hh"

//for Phoenix-to-photon rejection
#include <Stntuple/obj/TPhoenixElectronBlock.hh>

#endif

class TInit: public TStnModule {
//--------------------------------------------------------------

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnJetBlock*       fJetBlock;
		TStnTrackBlock*     fTrackBlock;
		TStnPhotonBlock*    fPhotonBlock;
		TStnElectronBlock*  fElectronBlock;
		TStnTriggerBlock*   fTriggerBlock;
		TStnVertexBlock*    fVertexBlock;
		TStnMetBlock*       fMetBlock;
		TCalDataBlock*      fCalBlock;
		TStnEmTimingBlock* fTimingBlock;
		TPhoenixElectronBlock* fPhxEleBlock;
		TStnTrackBlock* fPhxSiTrackBlock;
		TGenpBlock* 		  fGenpBlock;

		typedef std::vector<TCalTower*> CalDataArray;  // need for beam halo
		typedef std::vector<TCalTower*>::iterator CalDataArrayI;   // need for BH

	public:
		TInit(const char* name="Init", const char* title="Init");
		~TInit();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnJetBlock*      GetJetBlock()     { return fJetBlock;     }
		TStnTrackBlock*    GetTrackBlock()   { return fTrackBlock;   }
		TStnPhotonBlock*   GetPhotonBlock()  { return fPhotonBlock;  }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnTriggerBlock*  GetTriggerBlock() { return fTriggerBlock; }
		TStnVertexBlock*   GetVertexBlock()  { return fVertexBlock;  }


				/************** abstract objects *************************/


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
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int evtsPassLooseIDCuts;
			unsigned int evtsPassTightIDCuts;
			unsigned int evtsPassPhoLikeIDCuts;
			unsigned int haloPassPhoId;
			unsigned int haloCandidates;
			unsigned int cosmicCandidates;
			unsigned int cosmicPassPhoId;
			unsigned int phoenixPhotons;
			unsigned int phoenixPhoPassPID;
		};


		// all the histograms

		struct GeneralPlots_t {
			TH1F* photon_loose_and_tight_cuts_cumm;
			TH1F* Nvertices;
			TH1F* vertexZ;
			TH1F* pho_passIDs;
		}; //GeneralPlots_t

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

		// for PMT spike removal
	  struct CalHist_t {

		 TH2F*  fCalEMtwr;  // spike tower map for EM 
		 TH2F*  fCalCHAtwr;  // spike tower map for CHA 
		 TH2F*  fCalWHAtwr;  // spike tower map for WHA 
    		TH1F*  fCalHadTDCspike; // timing of HAD spikes;

		 TH1F*  fCalPmtEmAssm_bad;  // pmt assymetry for bad EM towers
		 TH1F*  fCalPmtHadAssm_bad;  // pmt assymetry for bad HAD towers
		 TProfile* fCalPmtAssmVsRun[19]; // pmt assymetry vs. Run num ber for bad towers
		 TH1F*  fCalBad_dPhi[19]; // dPhi=Phi(met)-Phi(bad tower) if |asymmetry|>0.55
		 TH1F*  fCalBad_MetRatio[19]; // MET/Et(bad tower) if |asymmetry|>0.55

		 TH1F*  fCalEm_dPhi; // dPhi=Phi(met)-Phi(em spike) if EmEnergy>10
		 TH1F*  fCalEm_MetRatio; // MET/Et(em spike) if EmEnergy>10 
		 TH1F*  fCalHad_dPhi; // dPhi=Phi(met)-Phi(had spike) if HadEnergy>10
		 TH1F*  fCalHad_MetRatio; // MET/Et(had spike) if HadEnergy>10 

		 TH1F*  fCalPmtEmAssm_le10;  // pmt assymetry for towers with EmEnergy<10
		 TH1F*  fCalPmtHadAssm_le10; // pmt assymetry for towers with HadEnergy<10
		 TH1F*  fCalPmtEmAssm_ge10;  // pmt assymetry for towers with EmEnergy>10
		 TH1F*  fCalPmtHadAssm_ge10; // pmt assymetry for towers with HadEnergy>10

		 TH1F*  fCalPmtCemAssm_le10;  //CEM pmt assymetry for towers with EmEnergy<10
		 TH1F*  fCalPmtChaAssm_le10; // CHA pmt assymetry for towers with HadEnergy<10, outside overlap
		 TH1F*  fCalPmtWhaAssm_le10; // WHA pmt assymetry for towers with HadEnergy<10, outside overlap
		 TH1F*  fCalPmtCemAssm_ge10;  // CEM pmt assymetry for towers with EmEnergy>10
		 TH1F*  fCalPmtChaAssm_ge10; // CHA pmt assymetry for towers with HadEnergy>10, outside overlap
		 TH1F*  fCalPmtWhaAssm_ge10; // WHA pmt assymetry for towers with HadEnergy>10, outside overlap

		 TH1F*  fCalPmtChaAssm_overlap_le10; // CHA pmt assymetry for towers with HadEnergy<10, inside overlap
		 TH1F*  fCalPmtWhaAssm_overlap_le10; // WHA pmt assymetry for towers with HadEnergy<10, inside overlap
		 TH1F*  fCalPmtChaAssm_overlap_ge10; // CHA pmt assymetry for towers with HadEnergy>10, inside overlap
		 TH1F*  fCalPmtWhaAssm_overlap_ge10; // WHA pmt assymetry for towers with HadEnergy>10, inside overlap

		 TH1F*  fCalPmtHadAssm_ge10_intime;  // pmt assymetry for towers with HadEnergy>10 and |HadTdc|<40.0 ns
		 TH1F*  fCalPmtHadAssm_ge10_outtime; // pmt assymetry for towers with HadEnergy>10 and |HadTdc|>40.0 ns
		 TH2F*  fCalEmPmtAssm_vs_PmtSum; // (pmt1-pmt2)/(pmt1+pmt2) for em energy
		 TH2F*  fCalHadPmtAssm_vs_PmtSum; // (pmt1-pmt2)/(pmt1+pmt2) for had energy
		 
		};

		struct BadTower {
			 int phiI;
			 int etaI;
			 int histoI;
			 int run1;
			 int run2;
		};

		// end PMT spike removal

		/*--------------------------- BEAM HALO STUFF -------------------*/
	  //-------------------------- parameters for Beam Halo
		struct HaloStuff {
			int seedwedge;
			int sidewedge;
			int seedwedgeH;
			int sidewedgeH;
			int eastNhad;
			int westNhad;
			double emseedE;
			double emsideE;
			double hadseedE;
			double seedwedgeHadE;
		};

		//______________________________ Beam Halo study histograms
		struct BeamHaloStudyHisto_t {

			//_____________________ My params
			TH1F* fHaloSeedWedgeH;   // number of seed wedge towers with Had>thresh 
			TH1F* fHaloSideWedgeH;   // number of side wedge towers with Em>thresh
			//______________________ Max's params
			TH1F* fHaloSeedWedge;   // number of seed wedge towers with Em>thresh 
			TH1F* fHaloSideWedge;   // number of side wedge towers with Em>thresh
			TH1F* fHaloEastNHad;    // number of had towers on East side
			TH1F* fHaloWestNHad;    // number of had towers on West side
			TH1F* fHaloEmSeedE;     // EM energy of seed wedge towers
			TH1F* fHaloEmSideE;     // EM energy of side wedge towers
			TH1F* fHaloHadSeedE;    // Plug HAD energy of seed wegde towers
			TH1F* fHaloSeedWedgeHadE; // HAD energy of seed wegde towers
			TH1F* fHaloNHad;        // number of towers on East & West: EastNHad+WestNHad

			TH2F* fHaloSeedWedgeEM_SeedWedgeHAD; // number of EM seed wedge towers vs. number of HAD towers 

			//______________________ Ray's params
			TH1F* fHaloEast;        // output of Pho->HaloEast()
			TH1F* fHaloWest;        // output of Pho->HaloWest()
			TH1F* fHaloEastWest;    // Pho->HaloEast()+Pho->HaloWest()
			TH1F* fHaloSideWedgeR;  // number of towers in side wedges according to Ray's code  
			TH1F* fHaloSeedWedgeR;  // number of towers in seed wedge according to Ray's code  
			TH1F* fHaloTwrInRawR;   // number of continues towers in seed wedge according to Ray's code  

			//______________________ my params
			TH1F* fHaloCesStripEoverE; // CesStripE/E(pho)
			TH1F* fHaloCesWireEoverE;  // CesWireE/E(pho)
			TH1F* fHaloHadTDC;         // Had TDC timing   

			TH1F* fHaloPhiWedge;  // wedge number of halo candidate
			TH1F* fHaloEta;       // eta of halo candidate
			TH1F* fHaloIso;       // iso of halo candidate
			TH1F* fHaloHadEm;     // Had/Em of halo candidate
			TH1F* fHaloCesChi2;   // CES Chi2 of halo candidate
			TH1F* fHaloEt;        // Et of halo candidate

			//my plots
			TH1F*  fHaloEt_to_Met;		//Beam halo candidate Et to Event Met ratio
			TH1F*  fHalo_Met_delPhi;		//Phi seperation between BH photon and MEt
			TH1F*  fHalo_Met_InvMass;		//Inv. mass of BH and MEt
			TH1F*  fnHalo_to_nPho;			// # of BH to photon ratio
			TH1F*	 fHalo_nPho_to_nJets;	// number of Photons to Jets in BH events
			TH2F*  fSeedWedge_HadTowers; //seed towers vs N Had towers

		};

		/*--------------------------- END BH     STUFF -------------------*/

		/*--------------------------- COSMIC RAY STUFF -------------------*/

		struct EmTimingHisto_t {
			 TH1F* fEmtTimePho; // timing of highest energy tower in photon
		};

		void SetTimeWindow(int i, double cut)          { if(i<2) timeWindow[i] = cut;}
		void SetPhoTimeThr(double cut)                 { phoTimeThr = cut; }
		void SetDelTimeThr(double cut)                 { dTimeThr = cut; }
		void SetUseCosmicsCut(int cut)                 { fUseCosmicsCut = cut; }
		void SetUseTimingCut(int cut)                  { fUseTimingCut = cut; }

		double  GetEmTiming(const int phoInd);
		int   MyCosmicsCut(TStnPhoton* Pho);
		void  BookTimingHistograms(EmTimingHisto_t& Hist, char* Folder); 
															

		/*------------------------END COSMIC RAY STUFF -------------------*/


		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void 	BookHistograms();
		void 	SaveHistograms();
		void 	Cleanup();
		void 	TagConversions();
		double 	ConversionSep(TStnTrack* t0, TStnTrack* t1, double & rconv);
		bool 	IsConversionElectron(TStnElectron* ele,
				     double& minsep, double& mindt, double& radius);

		int 	FiducialLepton(TLorentzVector, double Zv);
		void 	SetPhotonList();

		unsigned int CEMPhoEleIDcut(TStnPhoton*, TStnEvent*, double Zvtx);
		unsigned int TightPhotonIDcut(TStnPhoton*);
		unsigned int LoosePhotonIDcut(TStnPhoton*);

		TFolder* GetHistoFolder(char* name, char* title);
		void 	BookGeneralHistograms();
		void 	BookPhotonHistograms(PhotonPlots_t&, char* folder_name);
		void 	BookElectronHistograms(PhotonPlots_t&, char* folder_name);
		void 	BookCalHistograms(CalHist_t& Hist);

		void 	FillDataBlocks(int ientry);  // put all GetEntry statements here
		void 	FillPhotonPlots(TSuperPhoton&, PhotonPlots_t&);

		void SpikeStudy();
		void HaloStudy();
		void CosmicRayStudy();
		//methods for Sasha's code
		void  	InitForSashaCode(int);
		void SetMcDataType(int mc_type = 0) { fUseMcDataType = mc_type; } // 0=photon,1=Zee, 2=Wenu




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
		//need function to return an entire superphoton object
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

		TCalDataBlock*      GetCalDataBlock()  { return fCalBlock; }
		void SetUseFindSpike  (int cut)   { fUseFindSpike = cut; }
		void SetUseSpikeFilter(int cut)   { fUseSpikeFilter = cut; }
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

		void CreateBadTowerList();
		int FindPmtSpikes(CalDataArray* towers, double& EmESpike, double& HadESpike, double& spikeHadTdc, CalHist_t& Hist);
		void DoMyCalStuff(TStnJetBlock* fJetBlock,TStnPhotonBlock* fPhotonBlock,TCalDataBlock *fCalDataBlock,CalHist_t& Hist);
		static int towerInd(int iEta) { return  iEta >= TENETA/2 ? iEta - TENETA/2 : -(TENETA/2 - 1 - iEta); }

		// END STUFF FOR PMT SPIKE REMOVAL CODE

		// BEAM HALO STUFF
		
		void SetUseBeamHaloCut(int cut)	{ fUseBeamHaloCut=cut; }
		void SetHaloCutScenario(int cut)	{ HaloCutScenario=cut; }
		
		void  DoMyBeamHalo(int pho_ind);
		void  getBeamHaloInfo(int pho_ind, int& seedWedge, int& seedWedgeH, 
							int& sideWedge, int& sideWedgeH, int& eastNHads,
							int& westNHads, double& emSeedE, double& emSideE, 
							double& hadSeedE, double& seedWedgeHadE);
		int   BeamHaloCut(HaloStuff *beamhalo, int nvx12, int scenario=0);
		void  MatchCalorTowers(int pho_ind, CalDataArray* cdholder);
		int   Class12Vertices();
		void  FillBeamHaloStudyHistograms(BeamHaloStudyHisto_t& Hist, 
							 HaloStuff *beamhalo, TStnPhoton* Pho);

		void  BookBeamHaloStudyHistograms(BeamHaloStudyHisto_t& Hist, char* Folder);


		// END BEAM HALO STUFF



		// PHOENIX REJECTION
		int MyPhoFakeIDcut(TStnPhoton* Pho, TStnEvent* event);
		void PhoenixRejectionStudy();

		//END PHOENIX REJECTION




	private:
		bool 	qMc;  				// true if MonteCarlo data
		bool fUseMcDataType;	// tell which MC is used, Z/W...
		Global_Counters_t counter;

		TGoodRun 	goodrun;
		TStnDBManager 	DBman;
		TString 	goodRunListFile;
		TPhoTrigBits 	trigbits;
		

		//___________ for sasha's code
	  	int Nphotons;
		std::vector<TLorentzVector> phoUncorrected, phoCorrected;
		std::vector<float> phoHadEm;
		std::vector<float> phoEtaDet, phoCesX, phoCesZ;
		std::vector<int> phoIndex;

		//my super photons vector
		std::vector<TSuperPhoton> superpho;

		// histograms
		GeneralPlots_t GeneralPlot;
		PhotonPlots_t fPhotonHist, fPhoenixPhoHist, fPhoenixPhoPassPIDHist,fElectronHist;

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
		CalHist_t CalHist;

		// BH STUFF
		HaloStuff myHaloStuff;
		int fUseBeamHaloCut;      // parameter to invoke beam halo cut
		int HaloCutScenario;      // beam halo cut scenario
		BeamHaloStudyHisto_t fNoBeamHaloHist;   	// histograms for photons without beam halo candidate
		BeamHaloStudyHisto_t fWithBeamHaloHist; 	// histograms for photons with beam halo candidate
		BeamHaloStudyHisto_t fBH_passPIDHist; 		// histograms for beam halo candidate passing photon id
		BeamHaloStudyHisto_t fBH_failPIDHist; 		// histograms for beam halo candidate failing photon id
		BeamHaloStudyHisto_t fnotBH_passPIDHist; 	// histograms for non-beam halo candidate passing photon id
		BeamHaloStudyHisto_t fnotBH_failPIDHist; 	// histograms for non-beam halo candidate failing photon id

		//CR STUFF
		BeamHaloStudyHisto_t fCosmicHist; 	// histograms for cosmic ray candidates
		int fUseCosmicsCut;       // parameter to invoke cosmics cut
		int fUseTimingCut;        // parameter to invoke EM timing cut
		double timeWindow[2];     // time window for EM timing
		double phoTimeThr;        // time threshold for photons
		double dTimeThr;          // dT=T1-T2 time threshold
											 
		EmTimingHisto_t fTimingHisto; // EM timing histo
		int fUseVerbose;				// control print statements

	ClassDef(TInit,1)
};

bool SortTowersByEnergy(TCalTower* c1, TCalTower* c2);

bool FindInHepg(TGenpBlock* _fGenpBlock, std::vector<int> _parPdgStat,
					std::vector<TLorentzVector>& parvec, bool _Zmasscut = false);

#endif
