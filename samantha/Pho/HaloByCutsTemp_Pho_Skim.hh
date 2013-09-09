#ifndef HALOBYCUTSTEMP_PHO_SKIM_HH
#define HALOBYCUTSTEMP_PHO_SKIM_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/Pho/TagBeamHalo.hh"
#include "samantha/obj/Histograms.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/SetEmTimes.hh"
#endif

class HaloByCutsTemp_Pho_Skim: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		HaloByCutsTemp_Pho_Skim(const char* name="HaloByCutsTemp_Pho_Skim", const char* title="HaloByCutsTemp_Pho_Skim");
		~HaloByCutsTemp_Pho_Skim();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int haloEvtsFound;
		};


		// all the histograms


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


	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TagTightPhotons* tightMod;
		TagBeamHalo *haloMod;
		SetEmTimes *setTimeMod;

		bool bRunPermit;			//make sure we have all dependencies

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;
		

	ClassDef(HaloByCutsTemp_Pho_Skim,1)
};

#endif
