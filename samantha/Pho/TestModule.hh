///////////////////////////////////////////////////////////
//  //
//  //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#ifndef TESTMODULE_HH
#define TESTMODULE_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#endif

class TestModule: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		TestModule(const char* name="TestModule", const char* title="TestModule");
		~TestModule();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsProcessed;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
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
		bool bRunPermit;			// make sure all dependencies are met before running
	ClassDef(TestModule,1)
};

#endif
