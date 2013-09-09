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
#include "samantha/Pho/PlotEmTimes.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/utils/FreeFunctions.hh"

ClassImp(PlotEmTimes)

//_____________________________________________________________________________
PlotEmTimes::PlotEmTimes(const char* name, const char* title):
  TStnModule(name,title),
  bRunPermit(true)
{
	std::cout << "Hello I am PlotEmTimes module" << std::endl;
}

//_____________________________________________________________________________
PlotEmTimes::~PlotEmTimes() {
}

//_____________________________________________________________________________
void PlotEmTimes::SaveHistograms() {
}

//_____________________________________________________________________________
void PlotEmTimes::BookHistograms()
{
	DeleteHistograms();
	TFolder *folder = GetHistFolder(this, "TimingHist","AllTiming plots");
hist.PhoEmTimeVsRun_b4      = new TProfile("PhoEmTimeVsRun_b4","All Photons EM Time (before corrections)",270000,130000,400000,-200,200);
hist.TPhoEmTimeVsRun_b4     = new TProfile("TPhoEmTimeVsRun_b4","Tight Photons EM Time (before corrections)",270000,130000,400000,-200,200);
hist.TPho1JetEmTimeVsRun_b4 = new TProfile("TPho1JetEmTimeVsRun_b4","Tight Photons (+>=1Jet) EM Time (before corrections)",270000,130000,400000,-200,200);
hist.TPho2JetEmTimeVsRun_b4 = new TProfile("TPho2JetEmTimeVsRun_b4","Tight Photons (+2=1Jet) EM Time (before corrections)",270000,130000,400000,-200,200);
hist.PhoEmTimeVsRun_a4      = new TProfile("PhoEmTimeVsRun_a4","All Photons EM Time (after corrections)",270000,130000,400000,-200,200);
hist.TPhoEmTimeVsRun_a4     = new TProfile("TPhoEmTimeVsRun_a4","Tight Photons EM Time (after corrections)",270000,130000,400000,-200,200);
hist.TPho1JetEmTimeVsRun_a4 = new TProfile("TPho1JetEmTimeVsRun_a4","Tight Photons (+>=1Jet) EM Time (after corrections)",270000,130000,400000,-200,200);
hist.TPho2JetEmTimeVsRun_a4 = new TProfile("TPho2JetEmTimeVsRun_a4","Tight Photons (+2=1Jet) EM Time (after corrections)",270000,130000,400000,-200,200);
hist.PhoEmTime_b4  = new TH1F("PhoEmTime_b4","All L3 photons EM Times (before corrections)",300,-100,200);
hist.TPhoEmTime_b4  = new TH1F("TPhoEmTime_b4","All Tight photons EM Times (before corrections)",300,-100,200);
hist.TPho1JetEmTime_b4  = new TH1F("Pho1JetEmTime_b4","All tight photons (+>=1Jet) EM Times (before corrections)",300,-100,200);
hist.TPho2JetEmTime_b4  = new TH1F("Pho2JetEmTime_b4","All tight photons (+>=2Jets) EM Times (before corrections)",300,-100,200);
hist.PhoEmTime_a4  = new TH1F("PhoEmTime_a4","All L3 photons EM Times (after corrections)",300,-100,200);
hist.TPhoEmTime_a4  = new TH1F("TPhoEmTime_a4","All Tight photons EM Times (after corrections)",300,-100,200);
hist.TPho1JetEmTime_a4  = new TH1F("Pho1JetEmTime_a4","All tight photons (+>=1Jet) EM Times (after corrections)",300,-100,200);
hist.TPho2JetEmTime_a4  = new TH1F("Pho2JetEmTime_a4","All tight photons (+>=2Jets) EM Times (after corrections)",300,-100,200);

folder->Add(hist.PhoEmTimeVsRun_b4);
folder->Add(hist.TPhoEmTimeVsRun_b4);
folder->Add(hist.TPho1JetEmTimeVsRun_b4);
folder->Add(hist.TPho2JetEmTimeVsRun_b4);
folder->Add(hist.PhoEmTimeVsRun_a4);
folder->Add(hist.TPhoEmTimeVsRun_a4);
folder->Add(hist.TPho1JetEmTimeVsRun_a4);
folder->Add(hist.TPho2JetEmTimeVsRun_a4);
folder->Add(hist.PhoEmTime_b4);
folder->Add(hist.TPhoEmTime_b4);
folder->Add(hist.TPho1JetEmTime_b4);
folder->Add(hist.TPho2JetEmTime_b4);
folder->Add(hist.PhoEmTime_a4);
folder->Add(hist.TPhoEmTime_a4);
folder->Add(hist.TPho1JetEmTime_a4);
folder->Add(hist.TPho2JetEmTime_a4);
	
}


//_____________________________________________________________________________
int PlotEmTimes::BeginJob()
{
	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initSpMod == NULL) {
		StdOut(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
		bRunPermit = false;
	}
  	jetMod = (JetFilterModuleV2*) ((TStnAna*) GetAna()->GetModule("JetFilterV2"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
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
		bRunPermit = false;
	}
	
	return 0;
}

//_____________________________________________________________________________
int PlotEmTimes::BeginRun()
{
		if( cemT0Runs->hasRunCorrection(fHeaderBlock->RunNumber())){
			fCorrection = cemT0Runs->getRunCorrection(fHeaderBlock->RunNumber());
		} else {
			fCorrection = 0;
		}

	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int PlotEmTimes::Event(int ientry)
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
	std::vector<float> rawtime, cortime;
	
	int Nsp = initSpMod->GetSuperPhoSize();

	for (int i = 0; i < Nsp; i++) {
		double time = GetEmTiming(initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex());
		double fCTime = time - fCorrection; 

		hist.PhoEmTime_b4->Fill(time);
		hist.PhoEmTime_a4->Fill(fCTime);
		if (time >-300)  {  // some number >-999999, which indicate there is no timinig ingo
			hist.PhoEmTimeVsRun_b4->Fill(GetHeaderBlock()->RunNumber(),time);
			hist.PhoEmTimeVsRun_a4->Fill(GetHeaderBlock()->RunNumber(),fCTime);
		}
		
		if (initSpMod->GetSuperPhoton(i)->IsTightPhoton()) {
			hist.TPhoEmTime_b4->Fill(time);
			hist.TPhoEmTime_a4->Fill(fCTime);
			if (time >-300)  {  // some number >-999999, which indicate there is no timinig ingo
				hist.TPhoEmTimeVsRun_b4->Fill(GetHeaderBlock()->RunNumber(),time);
				hist.TPhoEmTimeVsRun_a4->Fill(GetHeaderBlock()->RunNumber(),fCTime);
			}
			rawtime.push_back(time);
			cortime.push_back(fCTime);
		}
	}

	if (rawtime.size()>0) {
		if (jetMod->GetMyNjet(0,3) >= 1) {
			hist.TPho1JetEmTime_b4->Fill(rawtime[0]);
			hist.TPho1JetEmTime_a4->Fill(cortime[0]);
			if (rawtime[0] >-300)  {  // some number >-999999, which indicate there is no timinig ingo
				hist.TPho1JetEmTimeVsRun_b4->Fill(GetHeaderBlock()->RunNumber(),rawtime[0]);
				hist.TPho1JetEmTimeVsRun_a4->Fill(GetHeaderBlock()->RunNumber(),cortime[0]);
			}
		}
		if (jetMod->GetMyNjet(0,3) >= 2) {
			hist.TPho2JetEmTime_b4->Fill(rawtime[0]);
			hist.TPho2JetEmTime_a4->Fill(cortime[0]);
			if (rawtime[0] >-300)  {  // some number >-999999, which indicate there is no timinig ingo
				hist.TPho2JetEmTimeVsRun_b4->Fill(GetHeaderBlock()->RunNumber(),rawtime[0]);
				hist.TPho2JetEmTimeVsRun_a4->Fill(GetHeaderBlock()->RunNumber(),cortime[0]);
			}
		}
	}
		
	if (GetPassed()) counter.evtsPassModule++;

	return 0;

} // Event


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void PlotEmTimes::MatchCalorTowers(int pho_ind, CalDataArray* cdholder)
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
  return;
}


/*-------------------------------------------------------------------*/
// runN>=190851  first run for EM timing system
/*-------------------------------------------------------------------*/
double PlotEmTimes::GetEmTiming(const int phoInd)
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
				
	      for (unsigned int icount = 0 ; icount < calholder.size(); icount++) {
				TCalTower* ctower = calholder[icount];
				int iEta = ctower->IEta();
				int iPhi = ctower->IPhi();
		  
				_emTimeTower = myEmTiming->TEmTimingModule::getTimingTower(0,iEta,iPhi); // 0=CEM; 1=PEM; 2=CHA; 3=WHA; 4=PHA 
		
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
void PlotEmTimes::FillDataBlocks(int ientry)
{
	fCalBlock->GetEntry(ientry);
	fPhotonBlock->GetEntry(ientry);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int PlotEmTimes::EndJob() {

	cemT0Runs->endJob();  // stops the EM time correction module
	std::string sMsg;


	printf("[PET:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[PET:01:] Events Processed ----------- = " << counter.evtsRunOver << std::endl;
	std::cout << "[PET:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;
	std::cout << "[PET:03:] Correction File Used ------- = " << GetEmTimeCorrectionFileName() << std::endl;
	printf("---------------------------------------------------\n");
	return 0;
}
