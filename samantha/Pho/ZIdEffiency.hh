#ifndef ZIDEFFIENCY_HH
#define ZIDEFFIENCY_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/utils/FreeFunctions.hh"
#endif

class ZIdEffiency: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TGenpBlock* 		  fGenpBlock;


	public:
		ZIdEffiency(const char* name="ZIdEffiency", const char* title="ZIdEffiency");
		~ZIdEffiency();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TGenpBlock* 		 GetGenpBlock() 	 { return fGenpBlock;    }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			//Z ele id cut eff
			unsigned int zideff_hepgZs;				//__ hepg Z events
			unsigned int zideff_hepgZees;				//__ Z that decayed to electron+ nutrino
			unsigned int zideff_passGoodrun;			//__ from the events with a Z, number of events pass good run
			unsigned int zideff_passVertex;			//__ then pass 1 good vertex
			unsigned int zideff_preZs;					// pre tag zee
			unsigned int zideff_passZv;					//__ then pass Zv cut
			unsigned int zideff_detZs;						//__ then that matching photon passing electron id cuts
		};


		// all the histograms
		struct ZPlots_t {	// Z electrons id eff plots
			TH1F* hepgZpt;
			TH1F* hepgZmass_b4masscut;
			TH1F* hepgZmass_a4masscut;
			TH1F* detZmass;
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
		bool FoundHepgZee();
		void BookZHistograms();

	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		ZPlots_t ZPlot;

	
	ClassDef(ZIdEffiency,1)
};

#endif
