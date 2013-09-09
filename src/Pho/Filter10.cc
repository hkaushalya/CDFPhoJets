#include "samantha/Pho/Filter10.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>

ClassImp(Filter10)

//_____________________________________________________________________________
Filter10::Filter10(const char* name, const char* title):
  TStnModule(name,title)
{
	std::cout << "Hello I am Filter10 module" << std::endl;
}

//_____________________________________________________________________________
Filter10::~Filter10() {
}

//_____________________________________________________________________________
void Filter10::SaveHistograms() {
}

//_____________________________________________________________________________
void Filter10::BookHistograms()
{
	DeleteHistograms();
}


//_____________________________________________________________________________
int Filter10::BeginJob()
{
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;

	return 0;
}

//_____________________________________________________________________________
int Filter10::BeginRun()
{
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int Filter10::Event(int ientry)
{
	SetPassed(0);
	counter.evtsRunOver++;
 
 	if (fHeaderBlock->EventNumber()%10 == 0) SetPassed(1);
	if (GetPassed()) counter.evtsPassModule++;
	
	return 0;
} // Event


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int Filter10::EndJob() {

	printf("[F10:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[F10:01:] Events Porcessed ----------- = " << counter.evtsRunOver << std::endl;
	std::cout << "[F10:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
