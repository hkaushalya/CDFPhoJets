///////////////////////////////////////////////////////////
// 1 in 10E8 events in STNTUPLES are known to return a   //
// burst of bad pointer.                                 //
// see, BU HEP Elog: STNtuple : entry #4                 //
// This check all poniters returned and skip that event  //
// if any bad pointer is found.                          //
// Got this from Ray Culbertson on 01-17-2008            //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#ifndef BADPOINTERCHECKER_HH
#define BADPOINTERCHECKER_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include <Stntuple/obj/TStnJetBlock.hh>
#include <Stntuple/obj/TCalDataBlock.hh>
#include "Stntuple/obj/TStnTrackBlock.hh"
#endif

class BadPointerChecker: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnTrackBlock*     fTrackBlock;
		TStnJetBlock*		  fJetBlockClu04;
		TCalDataBlock*      fCalDataBlock;


	public:
		BadPointerChecker(const char* name="BadPointerChecker", const char* title="BadPointerChecker");
		~BadPointerChecker();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock() { return fHeaderBlock;  }
		TStnTrackBlock*    GetTrackBlock()	{ return fTrackBlock; 	}
		TStnJetBlock*		 GetJetBlock()    { return fJetBlockClu04;}
		TCalDataBlock*     GetCalDataBlock() { return fCalDataBlock;}


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
		void FillDataBlocks(int ientry);

	private:
		Global_Counters_t counter;

	ClassDef(BadPointerChecker,1)
};

#endif
