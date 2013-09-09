///////////////////////////////////////////////////////////
// Stntuples some time have bad pointers to data blocks. //
// This will skip over those events.                     //
// Instead of running on every single event, since I use //
// fix set of datasets I'll create bad run,sec,event list//
// and use it at run time. The unpacking of cal and track//
// blocks are very time expensive. As I add more datasets//
// I'll have to run this over them and update this bad   //
// list. 12-27-2008                                      //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "samantha/utils/BadPointerChecker.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnTrack.hh"
#include "Stntuple/obj/TStnJet.hh"
#include "Stntuple/data/TCalTower.hh"
#include <iostream>

ClassImp(BadPointerChecker)

//_____________________________________________________________________________
BadPointerChecker::BadPointerChecker(const char* name, const char* title):
  TStnModule(name,title)
{
	std::cout << "Hello I am BadPointerChecker module" << std::endl;
}

//_____________________________________________________________________________
BadPointerChecker::~BadPointerChecker() {
}

//_____________________________________________________________________________
int BadPointerChecker::BeginJob()
{
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*) RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
  	RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlockClu04);
	RegisterDataBlock("CalDataBlock","TCalDataBlock",&fCalDataBlock);
	fTrackBlock 	= (TStnTrackBlock*) RegisterDataBlock("TrackBlock","TStnTrackBlock");

	// initialize vars
	counter.evtsProcessed		= 0;
	counter.evtsPassModule 		= 0;

	return 0;
}

//_____________________________________________________________________________
int BadPointerChecker::BeginRun()
{

	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int BadPointerChecker::Event(int ientry)
{
	SetPassed(0);
	counter.evtsProcessed++;
	
	int rc = 0;
	bool flagJet = false;
	bool flagTrk = false;
	bool flagCal = false;
	
	fJetBlockClu04->GetEntry(ientry);
	for(int ij=0; ij<fJetBlockClu04->NJets(); ij++) {
		TStnJet* jet = fJetBlockClu04->Jet(ij);
		if(int(jet)<1000) {
			flagJet=true;
		}
	}

	fTrackBlock->GetEntry(ientry);
	for(int it=0; it<fTrackBlock->NTracks(); it++) {
		TStnTrack* trk = fTrackBlock->Track(it);
		if(int(trk)<1000) {
			flagTrk=true;
		}
	}
	
	fCalDataBlock->GetEntry(ientry);
	for(int ic=0; ic<fCalDataBlock->NTowers(); ic++) {
		TCalTower* tow = fCalDataBlock->Tower(ic);
		if(int(tow)<1000) {
			flagCal=true;
		}
	}

	if (flagTrk || flagJet || flagCal) {
		std::cout<< __FILE__ << "::" << __LINE__ << "::" << GetHeaderBlock()->RunNumber() << "," << GetHeaderBlock()->SectionNumber() << "," << GetHeaderBlock()->EventNumber()<< ", BadPtCheck::Jet/Trk/Cal="<< flagJet << ","<< flagTrk << "," << flagCal << std::cout << std::endl;
		rc = 1;
	}
// if you want to know what is happening in this burst of bad pointer,
// uncomment following.
/* 
	if(flagTrk) {
		printf("N of trks:%d  %d\n",fTrackBlock->NTracks(),
		fTrackBlock->TrackList()->GetEntries());
		for(int it=0; it<fTrackBlock->NTracks(); it++) {
			TStnTrack* trk = fTrackBlock->Track(it);
			printf("%8x\n",int(trk));
		}
	} else {
		for(int it=0; it<fTrackBlock->NTracks(); it++) {
			TStnTrack* trk = fTrackBlock->Track(it);
			float x = trk->Pt();
		}
	}
	
	if(flagJet) {
		printf("N of jets:%d  %d\n",fJetBlockClu04->NJets(),
		fJetBlockClu04->ListOfJets()->GetEntries());
		for(int ij=0; ij<fJetBlockClu04->NJets(); ij++) {
			TStnJet* jet = fJetBlockClu04->Jet(ij);
			printf("%8x\n",int(jet));
		}
	} else {
		for(int ij=0; ij<fJetBlockClu04->NJets(); ij++) {
			TStnJet* jet = fJetBlockClu04->Jet(ij);
			float x = jet->Et();
		}
	}

	if(flagCal) {
		printf("N of towerss:%d  %d\n",fCalDataBlock->NTowers(),
		fCalDataBlock->fTowerList->GetEntries());
		for(int ic=0; ic<fCalDataBlock->NTowers(); ic++) {
			TCalTower* tow = fCalDataBlock->Tower(ic);
			printf("%8x\n",int(tow));
		}
	} else {
		for(int ic=0; ic<fCalDataBlock->NTowers(); ic++) {
			TCalTower* tow = fCalDataBlock->Tower(ic);
			float x = tow->EmEnergy();
		}
	}
*/																																																																																																				 

	if (!rc) SetPassed(1);
	if (GetPassed()) counter.evtsPassModule++;
	
	return rc;

} // Event


/*-------------------------------------------------------------------*/
void BadPointerChecker::FillDataBlocks(int ientry)
{
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int BadPointerChecker::EndJob() {

	printf("[BPC:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[BPC:01:] Events Processed ----------- = " << counter.evtsProcessed << std::endl;
	std::cout << "[BPC:02:] Events Pass this module ---- = " << counter.evtsPassModule << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
