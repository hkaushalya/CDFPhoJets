////////////////////////////////////////////////////////////
// This will create a root file with a tree. this tree   //
// will have Run# Event# and a branch for each module in //
// the TStnAna loop. initially all mods status will be   //
// set to 0(zero=fail). one they are passed, they will   //
// 1 reflect.                                            //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "samantha/utils/PassFailTree.hh"
#include "Stntuple/loop/TStnAna.hh"
#include <iostream>

ClassImp(PassFailTree)

//_____________________________________________________________________________
PassFailTree::PassFailTree(const char* name, const char* title):
  TStnModule(name,title),
  bMcFlag(false),
  sFileName("PassFailTree.root")
{
	std::cout << "Hello I am PassFailTree module" << std::endl;
}

//_____________________________________________________________________________
PassFailTree::~PassFailTree() {
}

//_____________________________________________________________________________
void PassFailTree::SaveHistograms() {
}

//_____________________________________________________________________________
void PassFailTree::BookHistograms()
{
	DeleteHistograms();
}


//_____________________________________________________________________________
int PassFailTree::BeginJob()
{
	rootFile = new TFile(sFileName.c_str(),"RECREATE");
	tree = new TTree("PassFail","Pass/fail status of the modules.");

  	modList = GetAna()->GetModuleList();
	
	if (modList->GetSize() > 0) {
		it = modList->MakeIterator();
		TStnModule* tstnMod;
		while (tstnMod = (TStnModule*) it->Next()) {
			vPassFail.push_back(0);
			vTstnMods.push_back(tstnMod);
			//tstnMod->SetPassed(0);
		}
		tstnMod =0;
	
		tree->Branch("mcflag",&bMcFlag,"mcflag/B");
		tree->Branch("run",&iPrevRun,"run/I");
		tree->Branch("event",&iPrevEvent,"event/I");
		// now create a brach for each module
		for (int i=0; i < vTstnMods.size(); i++) {
			tree->Branch(vTstnMods[i]->GetName(), &vPassFail[i],
							 (std::string(vTstnMods[i]->GetName())+"/I").c_str());
		}
	}

	
				// register the data blocks
	//fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();

	counter.evtsProcessed		= 0;
	counter.evtsPassModule 		= 0;
	bSkipThis = true;   // skip saving the first event. start from the following.
	iCurrRun = 0;
	iCurrEvent = 0;
	iPrevRun = 0;
	iPrevEvent = 0;
	
	return 0;
}

//_____________________________________________________________________________
int PassFailTree::BeginRun()
{
	
	if (GetHeaderBlock()->McFlag()) 
   	bMcFlag = true;
	else
   	bMcFlag = false;
 
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int PassFailTree::Event(int ientry)
{

	SetPassed(1);
	counter.evtsProcessed++;

	// postpone recording very first event. fill it in the next loop.
	// there is no other solution, i can think of, for now.
	if (bSkipThis) {
		iCurrRun = GetHeaderBlock()->RunNumber();
		iCurrEvent = GetHeaderBlock()->EventNumber();
		bSkipThis = false;		//skip only very first event
	} else {
		CheckAndFill();
	}
 
	if (GetPassed()) counter.evtsPassModule++;
	return 0;

} // Event


/*-------------------------------------------------------------------*/
void PassFailTree::CheckAndFill()
{
	for (int i=0; i <vTstnMods.size(); i++) {
		if (vTstnMods[i]->GetPassed()) vPassFail[i] = 1;		//mods passed
		else {
			vPassFail[i] = 0;		// first mod to fail
			for (int j = i+1; j <vTstnMods.size(); j++) {
				vPassFail[j] = -2;		// others did not get to process that event
												// why -2, bcos it is easy to see than -1 
			}
			break;
		}
		 //if (i!=0) vTstnMods[i]->SetPassed(0);		//reset all mods status to fail, except PassFailMod
	}
	iPrevRun = iCurrRun;
	iPrevEvent = iCurrEvent;
	
	tree->Fill();
	//tree->AutoSave();		//do not use unless debugging. just too slow!
	
	iCurrRun = GetHeaderBlock()->RunNumber();
	iCurrEvent = GetHeaderBlock()->EventNumber();
 
}

/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void PassFailTree::Cleanup()
{
} //Cleanup

//____________________________________________________________________
//  END JOB SUMMARY
//____________________________________________________________________
int PassFailTree::EndJob() {

	CheckAndFill();		//do this for one last time, to record the last event processed
	if (rootFile->IsOpen()) {
		rootFile->Write();
		rootFile->Close();
	}
	printf("[TST:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[TST:01:] Events Processed ----------- = " << counter.evtsProcessed << std::endl;
	std::cout << "[TST:02:] Events Pass this module ---- = " << counter.evtsPassModule << std::endl;
	std::cout << "[TST:03:] Pass/Fail Root file -------- = " << GetFileName() << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
