////////////////////////////////////////////////////////////
// This will create a flat stuples for me.                //
////////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "Stntuple/loop/TStnAna.hh"
#include "samantha/Pho/PhoEleFilter.hh"
#include <iostream>
#include "samantha/utils/FreeFunctions.hh"
#include <cmath>

ClassImp(PhoEleFilter)

//_____________________________________________________________________________
PhoEleFilter::PhoEleFilter(const char* name, const char* title):
  TStnModule(name,title),
  bRunPermit(true),
  bNoSummary(false),
  fMinEmEt(7),
  fMaxEmDetEta(1.1)
{
	std::cout << "Hello I am PhoEleFilter module" << std::endl;
}

//_____________________________________________________________________________
PhoEleFilter::~PhoEleFilter() {
}

//_____________________________________________________________________________
void PhoEleFilter::SaveHistograms() {
}

//_____________________________________________________________________________
void PhoEleFilter::BookHistograms()
{
}


//_____________________________________________________________________________
int PhoEleFilter::BeginJob()
{


  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}
	BookHistograms();

	counter.evtsProcessed		= 0;
	counter.evtsPassModule 		= 0;
	

	phoOnlyEvts = 0;
	phoEleEvts = 0;
	EleOnlyEvts = 0;

	cemEle =0;
	pemEle=0;
	cemPho=0; 
	extraPho=0; 
	
	return 0;
}

//_____________________________________________________________________________
int PhoEleFilter::BeginRun()
{
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int PhoEleFilter::Event(int ientry)
{
	if (! bRunPermit) exit (1);
	SetPassed(0);
	counter.evtsProcessed++;

	int Nele=0, Npho=0;
	for (int i=0; i < initSpMod->GetSuperPhoSize(); i++) {
		if (initSpMod->GetSuperPhoton(i)->IsTightPhoton() 
			|| initSpMod->GetSuperPhoton(i)->IsLoosePhoton()
			|| initSpMod->GetSuperPhoton(i)->IsTightElectron()			//photon-like electron
			|| initSpMod->GetSuperPhoton(i)->IsLooseElectron()			//photon-like electron, loose e/p
		//	|| initSpMod->GetSuperPhoton(i)->IsStdTightElectron()		// std. elecron cuts //not implemented yet
			|| initSpMod->GetSuperPhoton(i)->IsStdLooseElectron() )	// std. electron cuts
		{
			if (initSpMod->GetSuperPhoton(i)->GetEtCorr() > GetMinEmEt()
					&& fabs(initSpMod->GetSuperPhoton(i)->GetDetEta()) < GetMaxEmDetEta() )
			{
				SetPassed(1);
				break; //for this case, as I need to save any event with an EM obj
			}
		}
	}


	if (Npho >0 && Nele==0) ++phoOnlyEvts;
	if (Npho >0 && Nele > 0) ++phoEleEvts;
	if (Npho ==0 && Nele > 0) ++EleOnlyEvts;

	bool foundPho = false;
	int phoInd = -1;
	for (int i=0; i < initSpMod->GetSuperPhoSize(); i++)
	{
		if (initSpMod->GetSuperPhoton(i)->IsTightPhoton() 
		    && (initSpMod->GetSuperPhoton(i)->GetDetector() == 0)
			 && (initSpMod->GetSuperPhoton(i)->GetEtCorr() > GetMinEmEt())
			 && (!initSpMod->GetSuperPhoton(i)->IsPhoenix()) )
		{
			 foundPho = true;
			 phoInd = i;
			 cemPho++;
			 break;
		}
	}

	if (foundPho)
	{
		for (int i=0; i < initSpMod->GetSuperPhoSize(); i++)
		{
			if (i == phoInd) continue;
			if (initSpMod->GetSuperPhoton(i)->IsLoosePhoton()
				&& (initSpMod->GetSuperPhoton(i)->GetDetector() == 0) )
			{
				extraPho++;
			}

			if ( initSpMod->GetSuperPhoton(i)->IsStdLooseElectron()
				|| initSpMod->GetSuperPhoton(i)->IsLooseElectron()
				|| initSpMod->GetSuperPhoton(i)->IsPhoenix())
			{
				if (initSpMod->GetSuperPhoton(i)->GetDetector() == 0) 
				{
					cemEle++;
				} else if (initSpMod->GetSuperPhoton(i)->GetDetector() == 1) 
				{
					pemEle++;
				}
			}
		}
	}


	if (GetPassed()) counter.evtsPassModule++;
	return 0;

} // Event


//____________________________________________________________________
//  END JOB SUMMARY
//____________________________________________________________________
int PhoEleFilter::EndJob() {

	if (GetSummaryStat()) return 0;

	printf("[PEF:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[PEF:01:] Events Processed ----------- = " << counter.evtsProcessed << std::endl;
	std::cout << "[PEF:02:] Events Pass this module ---- = " << counter.evtsPassModule << std::endl;
	std::cout << "[PEF:03:] Min EM obj Et -------------- = " << GetMinEmEt() << " GeV" << std::endl;
	std::cout << "[PEF:04:] Max EM obj DetEta ---------- = " << GetMaxEmDetEta() << std::endl;
	std::cout << "[PEF:05:] Photons only evts   (7GeV) - = " << phoOnlyEvts << std::endl;
	std::cout << "[PEF:06:] Photons & Ele evts  (7GeV) - = " << phoEleEvts << std::endl;
	std::cout << "[PEF:07:] Electrons only evts (7GeV) - = " << EleOnlyEvts << std::endl;
	std::cout << "[PEF:08:] cemPho, extraPho    (7GeV) - = " << cemPho << ", " << extraPho << std::endl;
	std::cout << "[PEF:09:] cemEle, pemEle      (7GeV) - = " << cemEle<< ", " << pemEle << std::endl;
	printf("---------------------------------------------------\n");
	return 0;
}
