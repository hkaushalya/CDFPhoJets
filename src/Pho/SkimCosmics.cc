///////////////////////////////////////////////////////////
// This is Skim cosmics events to study.                 //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <numeric>
#include <vector>
#include <algorithm>
#include "samantha/Pho/SkimCosmics.hh"
#include "samantha/utils/FreeFunctions.hh"

ClassImp(SkimCosmics)

//_____________________________________________________________________________
SkimCosmics::SkimCosmics(const char* name, const char* title):
  TStnModule(name,title),
  bPermitRun(true),
  fCosmicTimeMin(30),
  fCosmicTimeMax(90)
{
	std::cout << "Hello I am SkimCosmics module" << std::endl;
}

//_____________________________________________________________________________
SkimCosmics::~SkimCosmics() {
}

//_____________________________________________________________________________
void SkimCosmics::SaveHistograms() {
}

//_____________________________________________________________________________
void SkimCosmics::BookHistograms()
{
	DeleteHistograms();
}

//_____________________________________________________________________________
int SkimCosmics::BeginJob()
{
	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bPermitRun = false;
	}

  	jetMod = (JetFilterModuleV2*) ((TStnAna*) GetAna()->GetModule("JetFilterV2"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		bPermitRun = false;
	}
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();

	// initialize vars
	iEvtsProcessed = 0;
	iEvtsPassed = 0;
	iCosmicEvts = 0;
	
	
	return 0;
}

//_____________________________________________________________________________
int SkimCosmics::BeginRun()
{
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int SkimCosmics::Event(int ientry)
{
	SetPassed(0);

	// this makes sure that all required mods exists
	if (!bPermitRun) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
		return 0;
	}
	
	iEvtsProcessed++;

		int Nsp = initSpMod->GetSuperPhoSize();
		for (int i = 0; i < Nsp; i++) {
			if (! initSpMod->GetSuperPhoton(i)->IsTightPhoton()) continue;
			if (initSpMod->GetSuperPhoton(i)->IsPhoenix()) continue;
			float fTime = initSpMod->GetSuperPhoton(i)->GetEmTime();

			if ((fTime > GetCosmicTimeMin() && fTime < GetCosmicTimeMax())) {
				SetPassed(1);
				iCosmicEvts++;
				break;
			}
		}

	if (GetPassed()) {
		iEvtsPassed++;
	}
	return 0;

} // Event


/*-------------------------------------------------------------------*/
void SkimCosmics::FillDataBlocks(int ientry)
{
}
//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int SkimCosmics::EndJob() {

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[CSD:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[CSD:01:] Events Processed ----------- = " + ToStr(iEvtsProcessed) + "\n";
	sMsg += "[CSD:02:] Events Passed -------------- = " + ToStr(iEvtsPassed) + "\n";
	sMsg += "[CSD:04:] Cosmic Event --------------- = " + ToStr(iCosmicEvts) + "\n";
	sMsg += "[CSD:08:] EM time window used  ------- = " + ToStr(GetCosmicTimeMin()) + "," +
																			ToStr(GetCosmicTimeMax()) + "\n";
	
	sMsg += "---------------------------------------------------\n";
	std::cout << sMsg;
	return 0;
}
