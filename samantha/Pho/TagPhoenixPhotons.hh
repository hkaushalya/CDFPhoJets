#ifndef TAGPHOENIXPHOTONS_HH
#define TAGPHOENIXPHOTONS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TPhoenixElectronBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnMetBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "TFolder.h"
#include <string>
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#endif

class TagPhoenixPhotons: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TPhoenixElectronBlock* fPhxEleBlock;
		TStnTrackBlock* fPhxSiTrackBlock;
		TStnMetBlock* fMetBlock;


	public:
		TagPhoenixPhotons(const char* name="TagPhoenixPhotons", const char* title="TagPhoenixPhotons");
		~TagPhoenixPhotons();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TPhoenixElectronBlock* GetPhoenixBlock() { return fPhxEleBlock; }
		TStnTrackBlock* GetPhxSiTrkBlock() { return fPhxSiTrackBlock; }
		TStnMetBlock* GetMetBlock() { return fMetBlock; }


				/************** abstract objects *************************/


		struct GlobalCounters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int phoenixCandidatesFound;
		};

		// all the histograms
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
			TH1F* EvtMet;
		};

		struct EventHist_t {
			TH1F* Met;
			TH1F* Sumet;
		};

		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
		void Cleanup();
		void FillDataBlocks(int ientry);
		int PhoenixCut(SuperPhoton*, TStnEvent* event);
		void BookPhotonHistograms(PhotonPlots_t&, std::string sFoldName, std::string sFoldTitle);
		void FillPhotonPlots(SuperPhoton*, PhotonPlots_t&);
		void FillEventHist(EventHist_t& hist);


		//setters
		int GetMode() { return iMode; }

		//getters
		void SetMode(int mode=0) {
			if (mode >=0 && mode <=3) iMode = mode;
			else {
				StdOut("Invalid mode!, valid 0-3. Using default(0) tagging only mode.",GetName());
			}
		}
		
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }

	private:
		bool 	qMc;  				// true if MonteCarlo data
		GlobalCounters_t counter;
		int iMode;
		PhotonPlots_t fPhoenixHist, fNoPhoenixHist;
		EventHist_t  fPhxHist, fNoPhxHist;	
		InitSuperPhotons *initSpMod;
		bool bRunPermit;		//to make sure we have all dependecies
		bool bNoSummary;			//control the end job summary print
	ClassDef(TagPhoenixPhotons,1)
};

#endif
