///////////////////////////////////////////////////////////
// I careted this to run over few selected events when   //
// debugging etc. I can use the TStnAna->ProcessEvent()  //
// for this. But if I want to save the event in Stntuple //
// format, it does not work with ProcessEvent().         //
// Provide all the run/event numbers to be processed.    //
// Everything else will be skipped. 04-28-2009           //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/*
 * $Id: EventFilter.cc,v 1.2 2009/06/26 20:35:13 samantha Exp $ 
 * $Log: EventFilter.cc,v $
 * Revision 1.2  2009/06/26 20:35:13  samantha
 * MINOR CHANGES: Header info in Event() is printed only in the debug mode.
 *
 * Revision 1.1  2009/05/22 17:08:26  samantha
 * To filter an specific run/event or combinations thereof.
 *
 *
 */


#include "samantha/Pho/EventFilter.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "Stntuple/loop/TStnAna.hh"


ClassImp(EventFilter)

//_____________________________________________________________________________
EventFilter::EventFilter(const char* name, const char* title):
  TStnModule(name,title),
  bRunPermit(true),
  iDebug(0)
{
	std::cout << "Hello I am EventFilter module" << std::endl;
}

//_____________________________________________________________________________
EventFilter::~EventFilter() {
}

//_____________________________________________________________________________
void EventFilter::SaveHistograms() {
}

//_____________________________________________________________________________
void EventFilter::BookHistograms()
{
	DeleteHistograms();
}


//_____________________________________________________________________________
int EventFilter::BeginJob()
{
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	if (!fHeaderBlock) {
		StdOut(__FILE__,__LINE__,3,"Required Header Block not found. pl. check!.");
		bRunPermit = false;
	}

	BookHistograms();

	// initialize vars
	counter.evtsProcessed 		= 0;
	counter.evtsPassModule 		= 0;

	return 0;
}

//_____________________________________________________________________________
int EventFilter::BeginRun()
{
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int EventFilter::Event(int ientry)
{
	SetPassed(0);
	counter.evtsProcessed++;

	if (iDebug)
	{
		if (vRun.size())
		{
			std::cout << "\tRun " << std::endl;
			for (std::set<unsigned>::iterator it= vRun.begin(); it != vRun.end(); it++)
			{
				std::cout << "\t" << *it << std::endl;
			}
		}
		if (vRunEvent.size())
		{
			std::cout << "\tRun\tEvent" << std::endl;
			for (std::set<std::pair<unsigned, unsigned> >::iterator it= vRunEvent.begin(); it != vRunEvent.end(); it++)
			{
				std::cout << "\t" << it->first << "\t" << it->second << std::endl;
			}
		}
	}

	std::pair<unsigned, unsigned> evt(fHeaderBlock->RunNumber(),fHeaderBlock->EventNumber());
	if ( (vRun.find(fHeaderBlock->RunNumber()) == vRun.end())
			&& (vRunEvent.find(evt) == vRunEvent.end())) return 0;
 
	if (iDebug) fHeaderBlock->Print();
 
 	SetPassed(1);
	if (GetPassed()) counter.evtsPassModule++;
	
	return 0;
} // Event


//-------------------------------------------------------------------
void EventFilter::AddRun(const unsigned int run)
{
	// Insert a run number to the list of runs to be processed.
	
	if (run >0)
	{
		if (vRun.find(run) == vRun.end()) vRun.insert(run);
	} else
	{
		std::string msg = "Run number " + ToStr(run) + " is <0! Not added!";
		StdOut(__FILE__,__LINE__,3,msg);
	}
}

//-------------------------------------------------------------------
void EventFilter::AddRunRange(const unsigned int lorun, const unsigned int hirun)
{
	// Insert a run range to the list of runs to be processed.

	if (lorun>0 && hirun>0 && (hirun>=lorun))
	{
		for (unsigned int run = lorun; run <= hirun; ++run)
		{
			AddRun(run);
		}
	} else {
		std::string msg("Invalid Run range. Not added to list. ");
			msg += "Please check. Runs must be >0 and lorun>=hirun.";

		StdOut(__FILE__,__LINE__,3,msg);
	}
}

//-------------------------------------------------------------------
void EventFilter::AddEvent(const unsigned int run, const unsigned int evt)
{
	// Insert an event to the list of runs to be processed.

	if (run>0 && evt>0)
	{
		std::pair<unsigned int, unsigned int> event(run,evt);
		vRunEvent.insert(event);
		//AddRun(run);
	} else {
		std::string msg("Invalid Run range. Not added to list. ");
						msg+=	"Please check. Runs must be >0 and lorun>=hirun.";
		StdOut(__FILE__,__LINE__,3,msg);
	}
}

//-------------------------------------------------------------------
void EventFilter::AddEventRange(const unsigned int run, const unsigned int evt1,
										const unsigned int evt2)
{
	// Insert a range of events to the list of runs to be processed.

	if (run>0 && evt1>0 &&  evt2>0 && (evt1>=evt2))
	{
		//AddRun(run);
		for (unsigned int evt = evt1; evt<=evt2; ++evt)
		{
			AddEvent(run,evt);
		}
	} else {
		std::string msg("Invalid Event range. Not added to list. ");
						msg+=	"Please check. Runs must be >0 and lorun>=hirun.";
		StdOut(__FILE__,__LINE__,3,msg);
	}
}


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int EventFilter::EndJob() {

	if (GetSummaryStat()) return 0;
	printf("[EVF:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[EVF:01:] Events Porcessed ----------- = " << counter.evtsProcessed << std::endl;
	std::cout << "[EVF:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;
	std::cout << "[EVF:03:] Processed Events ----------- = " << std::endl;

	if (vRun.size())
	{
		std::cout << "\tRun " << std::endl;
		for (std::set<unsigned>::iterator it= vRun.begin(); it != vRun.end(); it++)
		{
			std::cout << "\t" << *it << std::endl;
		}
	}
	if (vRunEvent.size())
	{
		std::cout << "\tRun\tEvent" << std::endl;
		for (std::set<std::pair<unsigned, unsigned> >::iterator it= vRunEvent.begin(); it != vRunEvent.end(); it++)
		{
			std::cout << "\t" << it->first << "\t" << it->second << std::endl;
		}
	}

	printf("---------------------------------------------------\n");
	return 0;
}
