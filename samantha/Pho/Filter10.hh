#ifndef FILTER10_HH
#define FILTER10_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/obj/Histograms.hh"
#endif

class Filter10: public TStnModule {

	protected:

		TStnHeaderBlock*    fHeaderBlock;

	public:
		Filter10(const char* name="Filter10", const char* title="Filter10");
		~Filter10();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
		};



		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void SaveHistograms();

	private:
	
	Global_Counters_t counter;

	ClassDef(Filter10,1)
};

#endif
