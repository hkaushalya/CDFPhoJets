///////////////////////////////////////////////////////////
//	This will tag the cosmic photons in the SuperPhoton   //
// list. for the early 400pb-1, where there is no EMtime //
// info, use trkless muon stubs to id cosmics. else use  //
// specified EM time window to identify the cosmics.     //
// - 01-24-2008                                          //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <numeric>
#include <vector>
#include <algorithm>
#include "Stntuple/loop/TStnAna.hh"
#include "samantha/Pho/SetEmTimes.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "../../RootTools/IOColors.hh"

ClassImp(SetEmTimes)

//_____________________________________________________________________________
SetEmTimes::SetEmTimes(const char* name, const char* title):
  TStnModule(name,title),
  bRunPermit(true),
  bNoSummary(false)
{
	std::cout << "Hello I am SetEmTimes module" << std::endl;
}

//_____________________________________________________________________________
SetEmTimes::~SetEmTimes() {
}

//_____________________________________________________________________________
void SetEmTimes::SaveHistograms() {
}

//_____________________________________________________________________________
void SetEmTimes::BookHistograms()
{
	DeleteHistograms();
}


//_____________________________________________________________________________
int SetEmTimes::BeginJob()
{
	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initSpMod == NULL) {
		StdOut(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
		bRunPermit = false;
	}
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fCalBlock 		= (TCalDataBlock*)     RegisterDataBlock("CalDataBlock","TCalDataBlock");
	fPhotonBlock 	= (TStnPhotonBlock*)   RegisterDataBlock("PhotonBlock","TStnPhotonBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	
	EmTimeCorrection::OperMode opmode = EmTimeCorrection::READM;
	if (sTimeCorrectionFile.length() >0) {
		cemT0Runs= new EmTimeCorrection(opmode,sTimeCorrectionFile);
	} else {
		StdOut(__FILE__,__LINE__,3,"File with timing corrections not specified!");
		//bRunPermit = false;
	}
	
	return 0;
}

//_____________________________________________________________________________
int SetEmTimes::BeginRun()
{
	fCorrection = 0;
	if (sTimeCorrectionFile.length() >0) 
	{
		if( cemT0Runs->hasRunCorrection(fHeaderBlock->RunNumber())){
			fCorrection = cemT0Runs->getRunCorrection(fHeaderBlock->RunNumber());
		}
	}

	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int SetEmTimes::Event(int ientry)
{
	SetPassed(1);
	counter.evtsRunOver++;
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"One or more requires mods not found!");
		SetPassed(0);
		exit (1);
		return 0;
	}

	FillDataBlocks(ientry);
	
	//std::cout << __FILE__ << ":" << __FUNCTION__ << ": ";
	//fHeaderBlock->Print();
	int Nsp = initSpMod->GetSuperPhoSize();
	for (int i = 0; i < Nsp; i++) {
		double time = GetEmTiming(initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex(), 
						initSpMod->GetSuperPhoton(i)->GetDetEta());
		double fCTime = time - fCorrection; 
		initSpMod->GetSuperPhoton(i)->SetEmTime(fCTime);

		//if (time>-150) std::cout << cyan << "pho " << i << " " << time << clearatt << std::endl;
	}

	if (GetPassed()) counter.evtsPassModule++;

	return 0;

} // Event


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void SetEmTimes::MatchCalorTowers(int pho_ind, CalDataArray* cdholder)
{
  cdholder->clear();
  TStnLinkBlock* links = fPhotonBlock->TowerLinkList();
  int nptow = links->NLinks(pho_ind);

  TCalTower* ctower = NULL;
  for(int j=0; j<nptow; j++) 
    {
      int iptow = links->Index(pho_ind,j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;

      ctower = fCalBlock->Tower(ieta,iphi);
      cdholder->push_back(ctower);
    }

  std::sort(cdholder->begin(), cdholder->end(), SortTowersByEnergy);
  //temporary : use only highest et tower
  //if (nptow>1)
  //{
  //for (int i=nptow; i>1 ;--i) cdholder->pop_back();
	//}
  //std::cout << "energy = " << cdholder->at(0)->Et() << std::endl;
  return;
}


/*-------------------------------------------------------------------*/
// runN>=190851  first run for EM timing system
/*-------------------------------------------------------------------*/
double SetEmTimes::GetEmTiming(const int phoInd, const float phoDetEta)
{
	double time = -999999.99;
	if (phoInd < 0) {
		StdOut(__FILE__,__LINE__,3,"Invalid photon index!.");
		return time;
	}
	
	int runN = GetHeaderBlock()->RunNumber();
	if (runN >= 190851 && GetHeaderBlock()->McFlag() == false) { // first run for EM timing system
		TEmTimingModule* myEmTiming = (TEmTimingModule*) ((TStnAna*) GetAna()->GetModule("EmTimingAna"));
		
      if (myEmTiming == NULL) {
			StdOut(__FILE__,__LINE__,3," EmTiminigModule not found.");
			return time;
		} else {		
      	//__________________ timing for photons
			int Nhits=-1;
			int NGoodHits=-1;
			EmTimingTower* _emTimeTower = NULL;
			CalDataArray calholder;
			
			MatchCalorTowers(phoInd, &calholder);
				
		//std::cout << __FUNCTION__ << ": photon towers" << std::endl;
		//std::cout << " i \t eta \t phi " << std::endl;

			int iCEMPEM = 0; // 0=CEM; 1=PEM;
			if (fabs(phoDetEta) > 1.1) iCEMPEM = 1;
			
	      for (unsigned int icount = 0 ; icount < calholder.size(); icount++) {
				TCalTower* ctower = calholder[icount];
				int iEta = ctower->IEta();
				int iPhi = ctower->IPhi();
				//std::cout << icount  << "\t" << iEta << "\t" << iPhi << std::endl;
		  
				_emTimeTower = myEmTiming->TEmTimingModule::getTimingTower(iCEMPEM,iEta,iPhi); // 0=CEM; 1=PEM; 2=CHA; 3=WHA; 4=PHA 
		
				if (TETTYP[iEta] == 5 && iPhi%2 != 0) _emTimeTower = myEmTiming->TEmTimingModule::getTimingTower(0,iEta,iPhi-1);
					
				if (_emTimeTower != NULL) {
					Nhits=_emTimeTower->nHits();
						
					for (int i=0; i< Nhits; i++) {
						int status = _emTimeTower->getStatus(i);
						  
						if (TStnEmTimingBlock::EnergyTooLow(status)) continue;
						if (TStnEmTimingBlock::SigmasFromThreshold(status) < 3) continue;
						if (_emTimeTower->getT0(i)<-80.0 || _emTimeTower->getT0(i)>160.0) continue;
						
						NGoodHits++;
						time=_emTimeTower->getT0(i);
						break;
		   	 	}
				}
			}
		}

	}
  return time;
}

/*-------------------------------------------------------------------*/
void SetEmTimes::FillDataBlocks(int ientry)
{
	fCalBlock->GetEntry(ientry);
	fPhotonBlock->GetEntry(ientry);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int SetEmTimes::EndJob() {

	if (GetSummaryStat()) return 0;
		

	if (sTimeCorrectionFile.length() >0) cemT0Runs->endJob();  // stops the EM time correction module
	std::string sMsg;


	printf("[SET:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[SET:01:] Events Processed ----------- = " << counter.evtsRunOver << std::endl;
	std::cout << "[SET:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;
	if (sTimeCorrectionFile.length() >0)
	std::cout << "[SET:03:] Correction File Used ------- = " << "NOT GIVEN. NO CORRECTIONS DONE." << std::endl;
	else 
	std::cout << "[SET:03:] Correction File Used ------- = " << GetEmTimeCorrectionFileName() << std::endl;
	printf("---------------------------------------------------\n");
	return 0;
}
