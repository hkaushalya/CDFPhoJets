#ifndef EVENTFILTER_HH
#define EVENTFILTER_HH

/*
 * $Id: EventFilter.hh,v 1.1 2009/05/22 17:08:09 samantha Exp $ 
 * $Log: EventFilter.hh,v $
 * Revision 1.1  2009/05/22 17:08:09  samantha
 * To filter an specific run/event or combinations thereof.
 *
 *
 */

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include <string>
#include <iostream>
#include <vector>
#include <set>

class EventFilter: public TStnModule {

	protected:

		TStnHeaderBlock*    fHeaderBlock;

	public:
		EventFilter(const char* name="EventFilter", const char* title="EventFilter");
		~EventFilter();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsProcessed;		//___number of evts we ran over
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
		void AddRun(const unsigned int run);
		void AddRunRange(const unsigned int run1, const unsigned int run2);
		void AddEvent(const unsigned int run, const unsigned int evt);
		void AddEventRange(const unsigned int run1, const unsigned int evt1, const unsigned int evt2);
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
		void SetDebug(const int de) { iDebug = de; }

	private:
		bool bRunPermit;			//make sure all dependencies are met
		bool bNoSummary;			//control the end job summary print
		Global_Counters_t counter;
		std::set<unsigned int> vRun;		//list of runs to be processed
											//list of runs and event number to be processed
		std::set<std::pair<unsigned int, unsigned int> > vRunEvent;		
		int iDebug;

	ClassDef(EventFilter,1)
};

#endif
