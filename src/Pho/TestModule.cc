///////////////////////////////////////////////////////////
//  //
//  //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "samantha/Pho/TestModule.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>

ClassImp(TestModule)

//_____________________________________________________________________________
TestModule::TestModule(const char* name, const char* title):
  TStnModule(name,title), qMc(false), bRunPermit(true)
{
	std::cout << "Hello I am TestModule module" << std::endl;
}

//_____________________________________________________________________________
TestModule::~TestModule() {
}

//_____________________________________________________________________________
void TestModule::SaveHistograms() {
}

//_____________________________________________________________________________
void TestModule::BookHistograms()
{
	DeleteHistograms();
  	//char name [200];
  	//char title[200];
	TFolder* new_folder = GetHistFolder(this, "Hist","Histos");
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}

}


//_____________________________________________________________________________
int TestModule::BeginJob()
{
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();

	// initialize vars
	counter.evtsProcessed		= 0;
	counter.evtsPassModule 		= 0;

	return 0;
}

//_____________________________________________________________________________
int TestModule::BeginRun()
{

	//int currRun =  fHeaderBlock->RunNumber();
	//std::cout << " BEGINING RUN# " << currRun << std::endl;
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TestModule::Event(int ientry)
{
	SetPassed(1);
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"One or more dependencies not found. pl check. Run Permit not cleared!");
		return 0;
	}
	counter.evtsProcessed++;
  
  	int NsuperPho = initSpMod->GetSuperPhoSize();
	
	for (int i=0; i < NsuperPho; i++) {
	}

	if (GetPassed()) counter.evtsPassModule++;
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TestModule::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void TestModule::FillDataBlocks(int ientry)
{
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TestModule::EndJob() {

	printf("[TST:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[TST:01:] Events Processed ------------ = " << counter.evtsProcessed << std::endl;
	std::cout << "[TST:02:] Events Pass this module ---- = " << counter.evtsPassModule << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
