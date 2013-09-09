#ifndef WSEL_HH
#define WSEL_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TStnElectronBlock.hh"
#include "Stntuple/obj/TStnMetBlock.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "samantha/Pho/TagConvElectrons.hh"
#endif

class WSel: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnElectronBlock*  fElectronBlock;
		TStnMetBlock*       fMetBlock;


	public:
		WSel(const char* name="WSel", const char* title="WSel");
		~WSel();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnMetBlock*      GetMetBlock()     { return fMetBlock; }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int Wevents;
		};

		struct RunDep_Counters_t {		//run dependent counters
			unsigned int run;				// run number
			unsigned int evts;			// events in the run
			unsigned int Ws;				//w candidates 
		};

		// all the histograms
		struct WPlots_t {	//__  plots
			TH1F* WMass;
			TH1F* WPt;
			TH1F* WEleEtMetRatio;
			TH1F* WSumet0;
			TH1F* WHt;
			TH1F* WMet;
			TH1F* WeleMetDelPhi;
		};
		struct ElePlots_t { //__ vars from the matching electron
			TH1F* EtCorr;
			TH1F* HadEm;
			TH1F* IsoEtCorr;
			TH1F* TrkPt;
			TH1F* EoverP;
			TH1F* Chi2Mean;
			TH1F* Ces2Wire;
			TH1F* Ces2Strip;
		};
		
		struct AcceptPlots_t {  // acceptance plots
			TH1F* Cumm;
			TH1F* All;
			TH1F* EventSum;		//events summary after each cut
		};

		struct RunDep_Hists_t {			//run dependent histograms counters
			TH1F* Nevts;			// events in the run
			TH1F* NWs;				//w candidates with run number 
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
		void BookWHist(WPlots_t& hist, std::string sFoldName, std::string sFoldTitle);
		void BookEleHist(ElePlots_t& hist, std::string sFoldName, std::string sFoldTitle, std::string sText);
		void BookAcceptanceHist(AcceptPlots_t& hist, std::string sFoldName, std::string sFoldTitle);
		void FillWPlots(WPlots_t& hist, SuperPhoton* sp);
		void FillElePlots(ElePlots_t& hist, SuperPhoton* sp);
		void FillAcceptancePlots(AcceptPlots_t& hist, SuperPhoton* sp);
		void BookRunDependentHistograms(RunDep_Hists_t& hist,
				std::string sFoldName, std::string sFoldTitle);

	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		WPlots_t hWPlot;
		ElePlots_t hWElePlot, hAllElePlot;
		InitSuperPhotons* initSpMod;
		TagConvElectrons* tagConvEleMod;
		TagElectrons* tagEleMod;
		TH1F *metb4cut,*meta4cut;
		AcceptPlots_t  hAccept;
		RunDep_Counters_t aRDCount_prev, aRDCount_curr;		//run dependent counters
		RunDep_Hists_t hRDHist;

		
	ClassDef(WSel,1)
};

#endif
