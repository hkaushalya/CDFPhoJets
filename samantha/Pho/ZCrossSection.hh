#ifndef ZCROSSSECTION_HH
#define ZCROSSSECTION_HH

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

class ZCrossSection: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnElectronBlock*  fElectronBlock;
		TStnMetBlock*       fMetBlock;


	public:
		ZCrossSection(const char* name="ZCrossSection", const char* title="ZCrossSection");
		~ZCrossSection();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnMetBlock*      GetMetBlock()     { return fMetBlock; }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int Zevents;
		};


		// all the histograms
		struct ZPlots_t {	//__  plots
			TH1F* ZMass;
			TH1F* ZPt;
			TH1F* ZEleEtratio;
			TH1F* ZSumet0;
			TH1F* ZMet;
			TH1F* ZeleDelPhi;
		};
		struct ElePlots_t { //__ vars from the matching electron
			TH1F* EtCorr;
			TH1F* HadEm;
			TH1F* IsoEtCorr;
			TH1F* TrkPt;
			TH1F* EoverP;
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
		void BookZHist(ZPlots_t& hist, std::string sFoldName, std::string sFoldTitle);
		void BookEleHist(ElePlots_t& hist, std::string sFoldName, std::string sFoldTitle);
		void FillZPlots(ZPlots_t& hist, SuperPhoton* sp1, SuperPhoton* sp2);
		void FillElePlots(ElePlots_t& hist, SuperPhoton* sp);

	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		ZPlots_t hZPlot;
		ElePlots_t hEle1Plot, hEle2Plot;
		InitSuperPhotons* initSpMod;
		TagConvElectrons* tagConvEleMod;
		TagElectrons* tagEleMod;
	
	ClassDef(ZCrossSection,1)
};

#endif
