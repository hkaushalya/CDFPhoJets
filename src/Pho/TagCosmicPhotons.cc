///////////////////////////////////////////////////////////
// this module will tag cosmic photons in the super      //
// photon list. use em timing for this. for the period   //
// without em timing, look for trackless muon stub       //
// around a 30 degree cone of the photon. -01-24-08      //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/* $Id: TagCosmicPhotons.cc,v 1.11 2009/11/10 19:43:09 samantha Exp $
 * $Log: TagCosmicPhotons.cc,v $
 * Revision 1.11  2009/11/10 19:43:09  samantha
 * ADDED: Auto CVS Log info.
 *
 *
 */

#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <numeric>
#include <vector>
#include <algorithm>
#include "Stntuple/loop/TStnAna.hh"
#include "samantha/Pho/TagCosmicPhotons.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/utils/FreeFunctions.hh"

ClassImp(TagCosmicPhotons)

//_____________________________________________________________________________
TagCosmicPhotons::TagCosmicPhotons(const char* name, const char* title):
	TStnModule(name,title),
	fCosmicTimeMin(30),
	fCosmicTimeMax(90),
	bRunPermit(true)
{
	std::cout << "Hello I am TagCosmicPhotons module" << std::endl;
}

//_____________________________________________________________________________
TagCosmicPhotons::~TagCosmicPhotons() {
}

//_____________________________________________________________________________
void TagCosmicPhotons::SaveHistograms() {
}

//_____________________________________________________________________________
void TagCosmicPhotons::BookHistograms()
{
	DeleteHistograms();
	
  	//char name [200];
  	//char title[200];
	TFolder* new_folder = GetHistFolder(this, "CountingRoom","Couting Histograms");
	
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		bRunPermit = false;
		return;
	} else {
		//list of counters
		Histograms HistoMan;
		vCounterHistLabels.push_back("Events Processed");							//=0
		vCounterHistLabels.push_back("Events Passed");								//=1
		vCounterHistLabels.push_back("Cosmics w.stubs(1^{st}400pb^{-1}");		//=2
		vCounterHistLabels.push_back("Cosmics in EM Timing Win.");				//=3
		vCounterHistLabels.push_back("Cosmics in EM Timing Win. w. stubs");	//=4
		vCounterHistLabels.push_back("Cosmics out of EM Timing Win. w. stubs");	//=5
		HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	}

}


//_____________________________________________________________________________
int TagCosmicPhotons::BeginJob()
{
	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initSpMod == NULL) {
		StdOut("Err ! Could not find module InitSuperPhotons!",0);
		bRunPermit = false;
		return 0;
	}
	
	setEmTimeMod = (SetEmTimes*) ((TStnAna*) GetAna()->GetModule("SetEmTimes"));
	if (setEmTimeMod == NULL) {
		StdOut(__FILE__,__LINE__,3,"Cannot find required SetEmTime module. pl. check.");
		bRunPermit = false;
		return 0;
	}
				// register the data blocks
				
	BookHistograms();

	// initialize vars
	
	
	return 0;
}

//_____________________________________________________________________________
int TagCosmicPhotons::BeginRun()
{
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TagCosmicPhotons::Event(int ientry)
{
	SetPassed(1);

	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"RunPermit failed!. pl. make sure all depndencies are met");
		SetPassed(1);
		return 0;
	}
	hCount.iCounter->Fill(0);
	
	int Nsp = initSpMod->GetSuperPhoSize();
	int Run = GetHeaderBlock()->RunNumber();	
	
	// how i reject possible cosmic events:
	// in 1st 400pb-1 with no em timing (RunNumber < 190851)
	// use if there is a trackless muon stub within the cone of 30 degrees of the photon, remove it.
	// in rest of the data use EM timing info
	// -6 ns to 6ns for prompt photon selection, may be i should make this narrower, from -4 to 4, to reduce BH.

	int iNcosmics = 0;

	if (Run<190851) { // first 400pb-1 with out em timing
		for (int i = 0; i < Nsp; i++) {
			int iNCosStubs = initSpMod->GetSuperPhoton(i)->GetPhoton()->NCosStubPho();  // look for any muon stubs with the 30 degree cone of the photon
			if (iNCosStubs>0) {
				initSpMod->GetSuperPhoton(i)->SetCosmicId(1);  //1=true=this is cosmic
				hCount.iCounter->Fill(2);
				iNcosmics++;
			} else {
				initSpMod->GetSuperPhoton(i)->SetCosmicId(0);  //0=false =this is not cosmic
			}
		}
	} else {

		for (int i = 0; i < Nsp; i++) {
			float fCTime = initSpMod->GetSuperPhoton(i)->GetEmTime();
			int iNCosStubs = initSpMod->GetSuperPhoton(i)->GetPhoton()->NCosStubPho();  // look for any muon stubs with the 30 degree cone of the photon
			//use en time and cosmic stub here too.
			if (fCTime > GetCosmicTimeMin() && fCTime < GetCosmicTimeMax()) {
				hCount.iCounter->Fill(3);
				initSpMod->GetSuperPhoton(i)->SetCosmicId(1);  //1=true=this is cosmic
				if (iNCosStubs>0) {
					hCount.iCounter->Fill(4);
				}
				iNcosmics++;
			} else {
				if (iNCosStubs>0) {
					hCount.iCounter->Fill(5);
					initSpMod->GetSuperPhoton(i)->SetCosmicId(1);  //1=true=this is cosmic
					iNcosmics++;
				} else {
					initSpMod->GetSuperPhoton(i)->SetCosmicId(0);  //0=false =this is not cosmic
				}
			}
		}
	}


	if (GetPassed()) {
		hCount.iCounter->Fill(1);
	}

	return 0;

} // Event

/*-------------------------------------------------------------------*/
void TagCosmicPhotons::FillDataBlocks(int ientry)
{
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TagCosmicPhotons::EndJob() {

	std::string sMsg;
	std::string sModName(GetName());

	sMsg += "[TCP:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[TCP:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n"; 
	sMsg += "[TCP:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n"; 
	sMsg += "[TCP:03:] 1st 400pb: w. stubs -------- = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n"; 
	sMsg += "[TCP:04:] in EM Time Window ---------- = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n"; 
	sMsg += "[TCP:05:] out of EM Time Win(w. stubs) = " + ToStr(hCount.iCounter->GetBinContent(5)) + "\n"; 
	printf("---------------------------------------------------\n");
	return 0;
}
