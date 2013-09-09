#ifndef TWZCROSSSECTION_HH
#define TWZCROSSSECTION_HH

#if !defined (__CINT__) || defined (__MAKECINT__)


#include <vector>
//#include "TTree.h"
//#include "TGraph.h"
#include <Stntuple/loop/TStnModule.hh>

#include <Stntuple/obj/TStnHeaderBlock.hh>
#include "Stntuple/obj/TStnTrackBlock.hh"
#include <Stntuple/obj/TStnPhotonBlock.hh>
#include <Stntuple/obj/TStnTriggerBlock.hh>
#include "Stntuple/obj/TStnElectronBlock.hh"
#include "Stntuple/obj/TStnVertexBlock.hh"
#include "Stntuple/obj/TStnMetBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"

#include <Stntuple/obj/TStnDBManager.hh>

#include "Stntuple/alg/TStntuple.hh"
#include "rlc/Pho/TGoodRun.hh"
#include "Stntuple/photon/TPhoTrigBits.hh"

#include "Stntuple/obj/TStnEvent.hh"

#endif

class TWZcrossSection: public TStnModule {
//--------------------------------------------------------------

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnTrackBlock*     fTrackBlock;
		TStnPhotonBlock*    fPhotonBlock;
		TStnElectronBlock*  fElectronBlock;
		TStnTriggerBlock*   fTriggerBlock;
		TStnVertexBlock*    fVertexBlock;
		TStnMetBlock*       fMetBlock;
		TGenpBlock* 		  fGenpBlock;

	public:
		TWZcrossSection(const char* name="WZcrossSection", const char* title="WZ cross section");
		~TWZcrossSection();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
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

			//W ele id cut eff
			unsigned int ideff_hepgWs;				//__ hepg W events
			unsigned int ideff_hepgWenus;			//__ W that decayed to electron+ nutrino
			unsigned int ideff_passGoodrun;		//__ from the events with a W, number of events pass good run
			unsigned int ideff_passVertex;		//__ then pass 1 good vertex
			unsigned int ideff_passZv;				//__ then pass Zv cut
			unsigned int ideff_passMetcut;		//___ need met>20
			unsigned int ideff_detWs;						//__ then that matching photon passing electron id cuts

			//Z ele id cut eff
			unsigned int zideff_hepgZs;				//__ hepg Z events
			unsigned int zideff_hepgZees;				//__ Z that decayed to electron+ nutrino
			unsigned int zideff_passGoodrun;			//__ from the events with a Z, number of events pass good run
			unsigned int zideff_passVertex;			//__ then pass 1 good vertex
			unsigned int zideff_preZs;					// pre tag zee
			unsigned int zideff_passZv;					//__ then pass Zv cut
			unsigned int zideff_detZs;						//__ then that matching photon passing electron id cuts
			// for data W/Z estimate
			unsigned int dWe_detWs;
			unsigned int dZe_detZs;

		};


		// all the histograms

		struct GeneralPlots_t {
			TH1F* Nvertices;
			TH1F* vertexZ;
			TH1F* CEMPHOEleIDeff;
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


		struct WEleIDeffPlots_t {	// W electron id eff plots
			TH1F* MinDelRMatch;
			TH1F* Et_ratio;
			TH1F* HepgEleEt;
			TH1F* Wele_IDeff;
		};

		struct ZEleIDeffPlots_t {	// Z electrons id eff plots
			TH1F* hepgZpt;
			TH1F* hepgZmass;
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

		};


		struct RunAvg_t {
			int runNumber;
			unsigned int events;		// total number of events in the run
			long double sumWMet;				// sum of Met of W events
			long double sumZMet;				// sum of Met of Z events
			unsigned int Wcount;			//sum of all W events upto and including the current run
			unsigned int Zcount;			//sum of all Z events upto and including the current run
		};

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
		void 		WEleIDcutEff();
		void 		ZEleIDcutEff();
		void 		WZ_dataEstimate();
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
		void     BookWeleIDeffHistograms();
		void     BookZeleIDeffHistograms();

		void 		FillDataBlocks(int ientry);  // put all GetEntry statements here
		void 		FillMetPlots();
	

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


	private:
		bool 	qMc;  				// true if MonteCarlo data
		bool 	trig25, trig50, trig70;
		Global_Counters_t counter;
		
		TGoodRun 		goodrun;
		TStnDBManager 	DBman;
		TString 			goodRunListFile;
		TPhoTrigBits 	trigbits;

		//my super photons vector
		std::vector<SuperPhoton_t> superpho;

		// histograms
		GeneralPlots_t GeneralPlot;
		MetPlots_t MetPlot;
		WEleIDeffPlots_t WEleIDeffPlot;
		ZEleIDeffPlots_t ZEleIDeffPlot;
		BackgroundPlots_t BackgroundPlot;

		//for collecting run avg values
		double sumWeventMet, sumZeventMet;
		int currRun, prevRun;
		std::vector<RunAvg_t> RunAvg;		
		
	ClassDef(TWZcrossSection,1)
};

#endif
