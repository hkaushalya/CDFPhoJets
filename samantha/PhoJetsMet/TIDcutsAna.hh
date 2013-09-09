#ifndef TIDCUTSANA_HH
#define TIDCUTSANA_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <vector>
#include "TTree.h"
#include "TGraph.h"
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
#include "Stntuple/ana/TSuperPhoton.hh"
#include "rlc/Pho/TGoodRun.hh"
#include "Stntuple/photon/TPhoTrigBits.hh"

#include "Stntuple/obj/TStnEvent.hh"


#endif

class TIDcutsAna: public TStnModule {
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
		TIDcutsAna(const char* name="GammaJetsInit", const char* title="GammaJetsInit");
		~TIDcutsAna();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnTrackBlock*    GetTrackBlock()   { return fTrackBlock;   }
		TStnPhotonBlock*   GetPhotonBlock()  { return fPhotonBlock;  }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnTriggerBlock*  GetTriggerBlock() { return fTriggerBlock; }
		TStnVertexBlock*   GetVertexBlock()  { return fVertexBlock;  }
		TGenpBlock* 		 GetGenpBlock() 	 { return fGenpBlock;    }


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

		enum PhoLikeEleIDBits_t{
			kCentralPho_E 	= 0x1 << 0,
			kEtCorr30Pho_E = 0x1 << 1,
			kHadEmPho_E 	= 0x1 << 2,
			kChi2Pho_E 		= 0x1 << 3,
			kN3d0_E 			= 0x1 << 4,
			kN3d2_E 			= 0x1 << 5,
			kPt2_E 			= 0x1 << 6,
			kFoundMatch_E  = 0x1 << 7,
			kEleTrkNum_E 	= 0x1 << 8,
			kEleTrkZ0_E 	= 0x1 << 9,
			kEleTrkPt_E 	= 0x1 << 10,
			kEleEoverP_E 	= 0x1 << 11,
			kPhoCalIso4_E 	= 0x1 << 12,
			kPhoTrkIso_E 	= 0x1 << 13,
			kPhoCes2Et_E 	= 0x1 << 14,
			kPhoXCes_E 		= 0x1 << 15,
			kPhoZCes_E 		= 0x1 << 16
		};

		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over 
			unsigned int evtsPassGoodrun;
			unsigned int evtsPassTrigger;
			unsigned int evtsPass25Trigger;
			unsigned int evtsPass50Trigger;
			unsigned int evtsPass70Trigger;
			unsigned int evtsPassAllTriggers;
			unsigned int evtsPass2550Triggers;
			unsigned int evtsPass2570Triggers;
			unsigned int evtsPass5070Triggers;
			unsigned int evtsPassVertexCut;
			unsigned int evtsPassZvCut;		
			unsigned int evtsPassModule;	//___number of evts pass this module

			unsigned int evtsPassLooseIDCuts;
			unsigned int evtsPassTightIDCuts;
         unsigned int evtsPassPhoLikeIDCuts;
		};

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

		struct GeneralPlots_t {
			TH1F* IdCummPlot;
			TH1F* Nvertices;
			TH1F* vertexZ;
			TH1F* pho_passIDs;
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

		unsigned int CEMPhoEleIDcut(TStnPhoton*, TStnEvent*, double Zvtx);

		TFolder* GetHistoFolder(char* name, char* title);
		void 		BookPhotonHistograms();
		void 		BookElectronHistograms();
		void 		BookGeneralHistograms();

		void 		FillDataBlocks(int ientry);  // put all GetEntry statements here
		void 		FillPhotonPlots(const TSuperPhoton&);
		void 		FillElectronPlots(const TSuperPhoton&);
		void		PhoLikeEleIDeff();
		void 		FillEleIDcummPlot(unsigned int);



	private:
		bool 	qMc;  				// true if MonteCarlo data
		bool 	trig25, trig50, trig70;
		Global_Counters_t counter;
		
		TGoodRun 		goodrun;
		TStnDBManager 	DBman;
		TString 			goodRunListFile;
		TPhoTrigBits 	trigbits;


		// histograms
		PhotonPlots_t PhotonPlot,ElectronPlot;
		GeneralPlots_t GeneralPlot;

	ClassDef(TIDcutsAna,1)
};

#endif
