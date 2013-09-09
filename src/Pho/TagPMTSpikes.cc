/////////////////////////////////////////////////////////////////////
// Searched for fake photons from PMT spikes. Talks to             //
// InitSuperPhotons classc and tag the photons accordingly.        // 
/////////////////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/*
 *	$Log: TagPMTSpikes.cc,v $
 *	Revision 1.9  2009/06/08 03:25:07  samantha
 *	MODIFIED: No tagging is done for MC samples.
 *	
 *
 */

#include "samantha/Pho/TagPMTSpikes.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include <iostream>

ClassImp(TagPMTSpikes)

//_____________________________________________________________________________
TagPMTSpikes::TagPMTSpikes(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  fUseFindSpike(0), 		// 0--don't look for spikes; 1--look for spikes
  fUseSpikeFilter(0), 	// 0==do nothing; 1==remove spikes; -1==select spikes
  fMinEmPMTasym(0.65),
  fMinHadPMTasym(0.85), 
  fMaxEmPMTasym(1.0),
  fMaxHadPMTasym(1.0),
  fMinEmPmtE(10.0), 
  fMinHadPmtE(10.0),
  iMode(0),
  bNoSummary(false)
{
	std::cout << "Hello I am TagPMTSpikes module" << std::endl;
}

//_____________________________________________________________________________
TagPMTSpikes::~TagPMTSpikes() {
}

//_____________________________________________________________________________
void TagPMTSpikes::SaveHistograms() {
}

//_____________________________________________________________________________
void TagPMTSpikes::BookHistograms()
{
	DeleteHistograms();
	BookCalHistograms(CalHist,"CalHist","CalHists");
}


//_____________________________________________________________________________
int TagPMTSpikes::BeginJob()
{
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fJetBlock 		= (TStnJetBlock*)      RegisterDataBlock("JetBlock","TStnJetBlock");
	fCalBlock 		= (TCalDataBlock*)     RegisterDataBlock("CalDataBlock","TCalDataBlock");
	fPhotonBlock 	= (TStnPhotonBlock*)   RegisterDataBlock("PhotonBlock","TStnPhotonBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 	  = 0;
	counter.evtsPassModule = 0;
	counter.spikeEvts      = 0;

	return 0;
}

//_____________________________________________________________________________
int TagPMTSpikes::BeginRun()
{

	//int currRun =  fHeaderBlock->RunNumber();
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TagPMTSpikes::Event(int ientry)
{
	SetPassed(0);
	counter.evtsRunOver++;
  

	if (fHeaderBlock->McFlag())
	{
		SetPassed(1);
	} else {

		InitSuperPhotons* initphomod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
		if (initphomod == NULL) std::cout << "InitPhotons is NULL" <<std::endl;

		int NsuperPho = initphomod->GetSuperPhoSize();

		FillDataBlocks(ientry);
		bool bSpikeFound = false;

		for (int i=0; i < NsuperPho; i++) {
			int iPhoIndex = initphomod->GetSuperPhoton(i)->GetPhotonBlockIndex();
			double EmESpike=0,HadESpike=0,spikeHadTdc=0;
			CalDataArray cdarray;
			MatchCalorTowers(iPhoIndex, fPhotonBlock, fCalBlock, &cdarray);
			int spikecode = FindPmtSpikes(&cdarray, EmESpike, HadESpike, spikeHadTdc, CalHist);
			if (spikecode > 0) bSpikeFound = true;
			initphomod->GetSuperPhoton(i)->SetPMTSpikeId(spikecode);
		}

		if (bSpikeFound) counter.spikeEvts++;

		if (GetMode() == 0) {		// by-pass mode
			SetPassed(1);
		} else if (GetMode() == 1 && !bSpikeFound) {     // remove PMT spike events
			SetPassed(1);
		} else if (GetMode() == 2 && bSpikeFound) {     // select PMT spike events
			SetPassed(1);
		}
	}
	
	if (GetPassed()) counter.evtsPassModule++;
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TagPMTSpikes::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
//________ creates a list of towers for all jets
/*-------------------------------------------------------------------*/
void TagPMTSpikes::MatchCalorTowers(TStnJet* jet,  TStnJetBlock* fJetBlock,
				TCalDataBlock *fCalDataBlock, CalDataArray* cdholder)
{
  cdholder->clear();
  TStnLinkBlock* links = fJetBlock->TowerLinkList();
  int nptow = links->NLinks(jet->Number());

  TCalTower* ctower = NULL;
  for(int j=0; j<nptow; j++)
    {
      int iptow = links->Index(jet->Number(),j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;

      ctower = fCalDataBlock->Tower(ieta,iphi);
      cdholder->push_back(ctower);
    }

  std::sort(cdholder->begin(), cdholder->end(),SortTowersByEnergy);

  return;
}

/*-------------------------------------------------------------------*/
//________ creates a list of towers for all EM objects
/*-------------------------------------------------------------------*/
void TagPMTSpikes::MatchCalorTowers(int pho_ind, TStnPhotonBlock* fPhotonBlock,
				TCalDataBlock *fCalDataBlock, CalDataArray* cdholder)
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

      ctower = fCalDataBlock->Tower(ieta,iphi);
      cdholder->push_back(ctower);
    }

  std::sort(cdholder->begin(),cdholder->end(),SortTowersByEnergy);

  return;
}


/*-------------------------------------------------------------------*/
// ___________ create a list of bad EM/HAD towers
/*-------------------------------------------------------------------*/
void TagPMTSpikes::CreateBadTowerList() 
{
	//_____________________________ creating a list of BAD EM towers
	BadTower _twr;
	_twr.phiI=0;
	_twr.etaI=2;
	_twr.histoI=0;
	_twr.run1=-1;
	_twr.run2=-1;
	badEmTwr.push_back(_twr);
	_twr.phiI=4;
	_twr.etaI=8;
	_twr.histoI=1;
	_twr.run1=-1;
	_twr.run2=-1;
	badEmTwr.push_back(_twr);
	_twr.phiI=4;
	_twr.etaI=7;
	_twr.histoI=2;
	_twr.run1=-1;
	_twr.run2=-1;
	badEmTwr.push_back(_twr);
	_twr.phiI=23;
	_twr.etaI=3;
	_twr.histoI=3;
	_twr.run1=-1;
	_twr.run2=-1;
	badEmTwr.push_back(_twr);
	//_____________________________ creating a list of BAD HAD towers
	_twr.phiI=2;
	_twr.etaI=-11;
	_twr.histoI=4;
	_twr.run1=182800;
	_twr.run2=186600;
	badHadTwr.push_back(_twr);
	_twr.phiI=0;
	_twr.etaI=-10;
	_twr.histoI=5;
	_twr.run1=164600;
	_twr.run2=169000;
	badHadTwr.push_back(_twr);
	_twr.phiI=23;
	_twr.etaI=-1;
	_twr.histoI=6;
	_twr.run1=138500;
	_twr.run2=156500;
	badHadTwr.push_back(_twr);
	_twr.phiI=5;
	_twr.etaI=7;
	_twr.histoI=7;
	_twr.run1=138000;
	_twr.run2=188000;
	badHadTwr.push_back(_twr);
	_twr.phiI=5;
	_twr.etaI=6;
	_twr.histoI=8;
	_twr.run1=138000;
	_twr.run2=188000;
	badHadTwr.push_back(_twr);
	_twr.phiI=22;
	_twr.etaI=8;
	_twr.histoI=9;
	_twr.run1=160500;
	_twr.run2=164500;
	badHadTwr.push_back(_twr);
	_twr.phiI=1;
	_twr.etaI=9;
	_twr.histoI=10;
	_twr.run1=152400;
	_twr.run2=154400;
	badHadTwr.push_back(_twr);
	_twr.phiI=12;
	_twr.etaI=10;
	_twr.histoI=11;
	_twr.run1=151200;
	_twr.run2=156600;
	badHadTwr.push_back(_twr);
	
}

/*-------------------------------------------------------------------*/
//___ finds sipkes, fills histo, returns spike Et
// loop over all towers occupies by the object,
// calculate pmt asym for each tower.
/*-------------------------------------------------------------------*/
int TagPMTSpikes::FindPmtSpikes(CalDataArray* towers, double& EmESpike,
								 double& HadESpike, double& spikeHadTdc,
								 CalHist_t& Hist) 
{
	int spike_code = 0;

	CalDataArrayI tti;
	EmESpike 	= 0.0;
	HadESpike 	= 0.0;
	spikeHadTdc = -1.0E6;
  
	for (tti = towers->begin(); tti != towers->end(); tti++) {
		TCalTower* calt = *tti;
      int iEta = calt->IEta();
      int iPhi = calt->IPhi();
      
		if (fabs(1.0 * towerInd(iEta)) > 11) continue; // cover all 2-pmt towers  //sam: pick only central towers

      int good_CEMtower	= 0; // should be ZERO for MY good towers  //sam? what is a good/bad tower
      int good_CHAtower	= 0; // should be ZERO for MY good towers
      int good_WHAtower	= 0; // should be ZERO for MY good towers
      float asymCEM		= -10.0;
      float asymCHAD		= -10.0;
      float asymWHAD		= -10.0;
      float pmt0			= 0.0;
      float pmt1			= 0.0;
      float cha_TDC		= -1.0E6;
      float wha_TDC		= -1.0E6;

      if (abs(towerInd(iEta)) <= 5) {	//_ central rowers
	  		pmt0 = calt->EmPmt(0);
	  		pmt1 = calt->EmPmt(1);
	  		
			if (pmt0+pmt1 > 200) asymCEM = (pmt0-pmt1)/(pmt0+pmt1);
	  
	  		pmt0 = calt->HadPmt(0);
	  		pmt1 = calt->HadPmt(1);
	  		cha_TDC = calt->HadTdc(0);
	  
	  		if (pmt0+pmt1 > 200) asymCHAD = (pmt0-pmt1)/(pmt0+pmt1);
		} //if
      
		if (abs(towerInd(iEta)) > 5 && abs(towerInd(iEta)) <= 7) {	//___ overlap region
	  		pmt0 = calt->EmPmt(0);
	  		pmt1 = calt->EmPmt(1);
	  		
			if (pmt0+pmt1 > 200) asymCEM = (pmt0-pmt1)/(pmt0+pmt1);
	  
	  		pmt0 = calt->HadPmt(0);
	  		pmt1 = calt->HadPmt(1);
	  		cha_TDC = calt->HadTdc(0);
	  
	  		if (pmt0+pmt1 > 200) asymCHAD = (pmt0-pmt1)/(pmt0+pmt1);
	  
	  		pmt0 = calt->HadPmt(2);
	  		pmt1 = calt->HadPmt(3);
	  		wha_TDC = calt->HadTdc(1);
	  		
			if (pmt0+pmt1 > 200) asymWHAD = (pmt0-pmt1)/(pmt0+pmt1);
		}
      
		if (abs(towerInd(iEta)) > 7 && abs(towerInd(iEta)) <=9) {	//___ overlap
	  		pmt0 = calt->EmPmt(0);
	  		pmt1 = calt->EmPmt(1);
	  
	  		if (pmt0+pmt1 > 200) asymCEM = (pmt0-pmt1)/(pmt0+pmt1);
	  		
			pmt0 = calt->HadPmt(0);
	  		pmt1 = calt->HadPmt(1);
	  		wha_TDC = calt->HadTdc(1);
	  		
			if (pmt0+pmt1 > 200) asymWHAD = (pmt0-pmt1)/(pmt0+pmt1);
		}
      
		if (abs(towerInd(iEta)) == 10) {
	  		pmt0 = calt->HadPmt(0);
	  		pmt1 = calt->HadPmt(1);
	  		wha_TDC = calt->HadTdc(1);
	  		
			if (pmt0+pmt1 > 200) asymWHAD = (pmt0-pmt1)/(pmt0+pmt1);
		}
      
		if (abs(towerInd(iEta)) == 11) {
			pmt0 = calt->HadPmt(2);
	  		pmt1 = calt->HadPmt(3);
	  		wha_TDC = calt->HadTdc(1);
	  
	  		if (pmt0+pmt1 > 200) asymWHAD = (pmt0-pmt1)/(pmt0+pmt1);
		}
		
      float EmEng		= calt->EmEnergy();
      float HadEng	= calt->HadEnergy();

      float twrEmEt 		= EmEng*fabs(sin(calt->Theta()*TMath::Pi()/180.0));
      float twrHadEt 	= HadEng*fabs(sin(calt->Theta()*TMath::Pi()/180.0));
      float twrEmPhi 	= calt->Phi();
      float twrHadPhi 	= calt->Phi();
      float dphimetEM 	= -1.0*TMath::Pi()/2.0;
      float dphimetHAD 	= -1.0*TMath::Pi()/2.0;
      float metEM			= -1.0;
      float metHAD		= -1.0;
      
		if (twrEmEt > 0.0) {
	  		metEM = fMetBlock->Met(0)/twrEmEt;			// MEt to Emt ratio and delPhi
	  		dphimetEM = fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(twrEmPhi) - TVector2::Phi_0_2pi(fMetBlock->MetPhi(0))));
		}
		
      if (twrHadEt>0.0) {
	  		metHAD=fMetBlock->Met(0)/twrHadEt;			// MEt to HadEt ratio and delPhi
	  		dphimetHAD=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(twrHadPhi)-TVector2::Phi_0_2pi(fMetBlock->MetPhi(0))));
		}

      if (fabs(asymCEM) <= 1.0) {	// why do we need this , asym is always <= 1 (=1 if one PMT is ZERO)
	  
	 		for (unsigned int ik = 0; ik < badEmTwr.size(); ik++) {
	      	if (towerInd(iEta) == badEmTwr[ik].etaI && iPhi == badEmTwr[ik].phiI) {
		  			good_CEMtower++;		//___ sam" this should be zero if all towers are good towers
		  			
					//sam: this is all bad tower stuff
					//if (EmEng>10.0) Hist.fCalPmtEmAssm_bad->Fill(asymCEM);
					
		  			//Hist.fCalPmtAssmVsRun[badEmTwr[ik].histoI]->Fill(GetHeaderBlock()->RunNumber(),asymCEM);
		  			
					if (fabs(asymCEM)>0.55 && EmEng>fMinEmPmtE) {
		      		//Hist.fCalBad_dPhi[badEmTwr[ik].histoI]->Fill(dphimetEM);
		      		//Hist.fCalBad_MetRatio[badEmTwr[ik].histoI]->Fill(metEM);
		    		}
		 			break;
				}
			} // for

			if (good_CEMtower == 0) {	
			
				if (EmEng<10.0) {
					Hist.fCalPmtEmAssm_le10->Fill(asymCEM); //sam: keep this, if object is central, this will be for central an so on
  					Hist.fCalPmtCemAssm_le10->Fill(asymCEM);
				} else {
					Hist.fCalPmtEmAssm_ge10->Fill(asymCEM);
		  			Hist.fCalPmtCemAssm_ge10->Fill(asymCEM);
				}
		
				Hist.fCalEmPmtAssm_vs_PmtSum->Fill(asymCEM,EmEng);
			} 
 	 
			if (fabs(asymCEM) > fMinEmPMTasym && fabs(asymCEM) <= fMaxEmPMTasym
					&& EmEng>fMinEmPmtE)
// 	  if(fabs(asymCEM)>0.5 && fabs(asymCEM)<=1.0 && EmEng>10.0)
//  	  if(fabs(asymCEM)>0.55 && fabs(asymCEM)<0.85 && EmEng>10.0)
//  	  if(asymCEM>0.55 && asymCEM<0.85 && EmEng>10.0)
//  	  if(asymCEM>-0.85 && asymCEM<-0.65 && EmEng>10.0)
			{
				spike_code++;
	      	EmESpike+=EmEng; // does it make sense to sum up energy for spikes?
	      	Hist.fCalEMtwr->Fill(towerInd(iEta),iPhi);
	      	Hist.fCalEm_dPhi->Fill(dphimetEM);
	      	Hist.fCalEm_MetRatio->Fill(metEM);
	    	}
		
		} // if asymCEM
      
		if (fabs(asymCHAD) <= 1.0) {
			
			for (unsigned int ik = 0; ik < badHadTwr.size(); ik++) {
	      	
				if (towerInd(iEta) == badHadTwr[ik].etaI && iPhi == badHadTwr[ik].phiI) {
		  			good_CHAtower++;
		  			//Hist.fCalPmtAssmVsRun[badHadTwr[ik].histoI]->Fill(GetHeaderBlock()->RunNumber(),asymCHAD);
		  			
					if ( (GetHeaderBlock()->RunNumber()) > badHadTwr[ik].run1 &&
						  (GetHeaderBlock()->RunNumber()) < badHadTwr[ik].run2 ) {
		     			
						//if (HadEng>10.0) Hist.fCalPmtHadAssm_bad->Fill(asymCHAD);
		      		
						if (fabs(asymCHAD)>0.55 && HadEng>fMinHadPmtE) {
			  				//Hist.fCalBad_dPhi[badHadTwr[ik].histoI]->Fill(dphimetHAD);
			  				//Hist.fCalBad_MetRatio[badHadTwr[ik].histoI]->Fill(metHAD);
						}
		    		}
					break;
				} 
	    	} // for 

			if (good_CHAtower == 0) {		
	      	if (HadEng<10.0) {		//___since we are looking at photons (electrons)
		  			
					Hist.fCalPmtHadAssm_le10->Fill(asymCHAD);
		  			if (abs(towerInd(iEta)) < 6) Hist.fCalPmtChaAssm_le10->Fill(asymCHAD);
		  			else Hist.fCalPmtChaAssm_overlap_le10->Fill(asymCHAD);
				
				} else {
		  		
					Hist.fCalPmtHadAssm_ge10->Fill(asymCHAD);
		  			
					if (abs(towerInd(iEta))<6) Hist.fCalPmtChaAssm_ge10->Fill(asymCHAD);
		  			else Hist.fCalPmtChaAssm_overlap_ge10->Fill(asymCHAD);
		  			
					if (fabs(cha_TDC)<40.0) Hist.fCalPmtHadAssm_ge10_intime->Fill(asymCHAD);
		  			if (fabs(cha_TDC)>40.0 && fabs(cha_TDC)<80.0) Hist.fCalPmtHadAssm_ge10_outtime->Fill(asymCHAD);
				}
	      	Hist.fCalHadPmtAssm_vs_PmtSum->Fill(asymCHAD,HadEng);
			} // if
 	  
			if (fabs(asymCHAD) > fMinHadPMTasym &&
				 fabs(asymCHAD) <= fMaxHadPMTasym && HadEng > fMinHadPmtE)
				// if(fabs(asymCHAD)>0.5 && fabs(asymCHAD)<=1.0 && HadEng>10.0)
	    		{
	      	spike_code++;
	      	HadESpike += HadEng; // does it make sense to sum up energy for spikes?
	      	spikeHadTdc=cha_TDC;
	      	Hist.fCalHadTDCspike->Fill(cha_TDC); // timing of HAD spikes;
	      	Hist.fCalCHAtwr->Fill(towerInd(iEta),iPhi);
	      	Hist.fCalHad_dPhi->Fill(dphimetHAD);
	      	Hist.fCalHad_MetRatio->Fill(metHAD);
	    	}
		} // if asymCHAD
      
		
		if (fabs(asymWHAD) <= 1.0) {
		
			for (unsigned int ik=0; ik<badHadTwr.size(); ik++) {
	      
				if (towerInd(iEta) == badHadTwr[ik].etaI && iPhi == badHadTwr[ik].phiI) {
		  			good_WHAtower++;
		  			//Hist.fCalPmtAssmVsRun[badHadTwr[ik].histoI]->Fill(GetHeaderBlock()->RunNumber(),asymWHAD);
		  
		  			if ( (GetHeaderBlock()->RunNumber()) > badHadTwr[ik].run1 && 
						  (GetHeaderBlock()->RunNumber()) < badHadTwr[ik].run2 ) {
		      
						//if (HadEng>10.0) Hist.fCalPmtHadAssm_bad->Fill(asymWHAD);
		      
						if (fabs(asymWHAD)>0.55 && HadEng>fMinHadPmtE) {
			  				//Hist.fCalBad_dPhi[badHadTwr[ik].histoI]->Fill(dphimetHAD);
			  				//Hist.fCalBad_MetRatio[badHadTwr[ik].histoI]->Fill(metHAD);
						}
		    		}
		  			break;
				}
			} // for

			if (good_WHAtower == 0) {
	      	if (HadEng<10.0) {
		 			Hist.fCalPmtHadAssm_le10->Fill(asymWHAD);
		  			
					if (abs(towerInd(iEta)) > 7 && abs(towerInd(iEta)) < 11) {
						Hist.fCalPmtWhaAssm_le10->Fill(asymWHAD);
		  			} else {
						Hist.fCalPmtWhaAssm_overlap_le10->Fill(asymWHAD);
					}
			
				} else {
		  		
					Hist.fCalPmtHadAssm_ge10->Fill(asymWHAD);
		  			if (abs(towerInd(iEta)) > 7 && abs(towerInd(iEta)) < 11) {
						Hist.fCalPmtWhaAssm_ge10->Fill(asymWHAD);
		  			} else {
						Hist.fCalPmtWhaAssm_overlap_ge10->Fill(asymWHAD);
					}
		  
		  			if (fabs(wha_TDC) < 40.0) {
						Hist.fCalPmtHadAssm_ge10_intime->Fill(asymWHAD);
		  			}
				
					if (fabs(wha_TDC) > 40.0 && fabs(wha_TDC) < 80.0) {
						Hist.fCalPmtHadAssm_ge10_outtime->Fill(asymWHAD);
					}
				}
	      	Hist.fCalHadPmtAssm_vs_PmtSum->Fill(asymWHAD,HadEng);
			}
	  	
			if (fabs(asymWHAD) > fMinHadPMTasym &&
				 fabs(asymWHAD) <= fMaxHadPMTasym && 
				 HadEng > fMinHadPmtE ) {

				spike_code++;
				HadESpike += HadEng; // does it make sense to sum up energy for spikes?
				spikeHadTdc = wha_TDC;
				Hist.fCalHadTDCspike->Fill(wha_TDC); // timing of HAD spikes;
				Hist.fCalWHAtwr->Fill(towerInd(iEta),iPhi);
				Hist.fCalHad_dPhi->Fill(dphimetHAD);
				Hist.fCalHad_MetRatio->Fill(metHAD);
			}
		}
   }

  return spike_code;
}


/*-------------------------------------------------------------------*/
void TagPMTSpikes::FillDataBlocks(int ientry)
{
	fJetBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);
	fCalBlock->GetEntry(ientry);
}

/*-------------------------------------------------------------------*/
void TagPMTSpikes::BookCalHistograms(CalHist_t& Hist, std::string sFoldName,
								std::string sFoldTitle)
{
	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}

  char name [200];
  char title[200];

  sprintf(name,"CalEMtwr");
  sprintf(title," EM towers with spikes: iPhi vs Tower Index");
  Hist.fCalEMtwr=new TH2F(name,title,23,-11.5,11.5,24,-0.5,23.5);
  Hist.fCalEMtwr->SetMarkerStyle(kPlus);
  new_folder->Add(Hist.fCalEMtwr);
 
  sprintf(name,"CalCHAtwr");
  sprintf(title," CHA towers with spikes: iPhi vs Tower Index");
  Hist.fCalCHAtwr=new TH2F(name,title,23,-11.5,11.5,24,-0.5,23.5);
  Hist.fCalCHAtwr->SetMarkerStyle(kPlus);
  new_folder->Add(Hist.fCalCHAtwr);
  
  sprintf(name,"CalWHAtwr");
  sprintf(title," WHA towers with spikes: iPhi vs Tower Index");
  Hist.fCalWHAtwr=new TH2F(name,title,23,-11.5,11.5,24,-0.5,23.5);
  Hist.fCalWHAtwr->SetMarkerStyle(kPlus);
  new_folder->Add(Hist.fCalWHAtwr);

  sprintf(name,"CalHadTDCspike");
  sprintf(title,"Timing of HAD spikes (CHA & WHA)");
  Hist.fCalHadTDCspike=new TH1F(name,title,1000,-200.0,200.0);
  new_folder->Add(Hist.fCalHadTDCspike);

	//____________________________________ bad tower studies
/*{{{*/
/*
  sprintf(name,"CalPmtAssmVsRun[0]");
  sprintf(title,"CEM tower (0,2): (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[0]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[0]);
  sprintf(name,"CalPmtAssmVsRun[1]");
  sprintf(title,"CEM tower (4,8): (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[1]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[1]);
  sprintf(name,"CalPmtAssmVsRun[2]");
  sprintf(title,"CEM tower (4,7): (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[2]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[2]);
  sprintf(name,"CalPmtAssmVsRun[3]");
  sprintf(title,"CEM tower (23,3): (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[3]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[3]);
  sprintf(name,"CalPmtAssmVsRun[4]");
  sprintf(title,"WHA tower (2,-11): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[4]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[4]);
  sprintf(name,"CalPmtAssmVsRun[5]");
  sprintf(title,"WHA tower (0,-10): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[5]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[5]);
  sprintf(name,"CalPmtAssmVsRun[6]");
  sprintf(title,"CHA tower (23,-1): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[6]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[6]);
  sprintf(name,"CalPmtAssmVsRun[7]");
  sprintf(title,"CHA tower (5,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[7]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[7]);
  sprintf(name,"CalPmtAssmVsRun[8]");
  sprintf(title,"CHA tower (5,6): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[8]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[8]);
  sprintf(name,"CalPmtAssmVsRun[9]");
  sprintf(title,"WHA tower (22,8): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[9]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[9]);
  sprintf(name,"CalPmtAssmVsRun[10]");
  sprintf(title,"WHA tower (1,9): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[10]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[10]);
  sprintf(name,"CalPmtAssmVsRun[11]");
  sprintf(title,"WHA tower (12,10): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[11]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[11]);
  sprintf(name,"CalPmtAssmVsRun[12]");
  sprintf(title,"CHA tower (11,-7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[12]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[12]);
  sprintf(name,"CalPmtAssmVsRun[13]");
  sprintf(title,"CHA tower (12,-7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[13]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[13]);
  sprintf(name,"CalPmtAssmVsRun[14]");
  sprintf(title,"CHA tower (7,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[14]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[14]);
  sprintf(name,"CalPmtAssmVsRun[15]");
  sprintf(title,"CHA tower (9,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[15]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[15]);
  sprintf(name,"CalPmtAssmVsRun[16]");
  sprintf(title,"CHA tower (10,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[16]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[16]);
  sprintf(name,"CalPmtAssmVsRun[17]");
  sprintf(title,"CHA tower (14,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[17]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[17]);
  sprintf(name,"CalPmtAssmVsRun[18]");
  sprintf(title,"CHA tower (19,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[18]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[18]);
*/
/*}}}*/

/*{{{*/
/*
  sprintf(name,"CalBad_dPhi[0]");
  sprintf(title,"CEM tower (0,2): d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[0]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[0]);
  sprintf(name,"CalBad_MetRatio[0]");
  sprintf(title,"CEM tower (0,2): #slash{E}_{T}/E_{T}^{EM spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[0]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[0]);
  sprintf(name,"CalBad_dPhi[1]");
  sprintf(title,"CEM tower (4,8): d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[1]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[1]);
  sprintf(name,"CalBad_MetRatio[1]");
  sprintf(title,"CEM tower (4,8): #slash{E}_{T}/E_{T}^{EM spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[1]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[1]);
  sprintf(name,"CalBad_dPhi[2]");
  sprintf(title,"CEM tower (4,7): d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[2]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[2]);
  sprintf(name,"CalBad_MetRatio[2]");
  sprintf(title,"CEM tower (4,7): #slash{E}_{T}/E_{T}^{EM spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[2]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[2]);
  sprintf(name,"CalBad_dPhi[3]");
  sprintf(title,"CEM tower (23,3): d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[3]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[3]);
  sprintf(name,"CalBad_MetRatio[3]");
  sprintf(title,"CEM tower (23,3): #slash{E}_{T}/E_{T}^{EM spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[3]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[3]);
  sprintf(name,"CalBad_dPhi[4]");
  sprintf(title,"WHA tower (2,-11): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[4]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[4]);
  sprintf(name,"CalBad_MetRatio[4]");
  sprintf(title,"WHA tower (2,-11): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[4]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[4]);
  sprintf(name,"CalBad_dPhi[5]");
  sprintf(title,"WHA tower (0,-10): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[5]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[5]);
  sprintf(name,"CalBad_MetRatio[5]");
  sprintf(title,"WHA tower (0,-10): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[5]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[5]);
  sprintf(name,"CalBad_dPhi[6]");
  sprintf(title,"CHA tower (23,-1): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[6]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[6]);
  sprintf(name,"CalBad_MetRatio[6]");
  sprintf(title,"CHA tower (23,-1): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[6]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[6]);
  sprintf(name,"CalBad_dPhi[7]");
  sprintf(title,"CHA tower (5,7): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[7]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[7]);
  sprintf(name,"CalBad_MetRatio[7]");
  sprintf(title,"CHA tower (5,7): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[7]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[7]);
  sprintf(name,"CalBad_dPhi[8]");
  sprintf(title,"CHA tower (5,6): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[8]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[8]);
  sprintf(name,"CalBad_MetRatio[8]");
  sprintf(title,"CHA tower (5,6): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[8]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[8]);
  sprintf(name,"CalBad_dPhi[9]");
  sprintf(title,"WHA tower (22,8): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[9]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[9]);
  sprintf(name,"CalBad_MetRatio[9]");
  sprintf(title,"WHA tower (22,8): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[9]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[9]);
  sprintf(name,"CalBad_dPhi[10]");
  sprintf(title,"WHA tower (1,9): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[10]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[10]);
  sprintf(name,"CalBad_MetRatio[10]");
  sprintf(title,"WHA tower (1,9): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[10]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[10]);
  sprintf(name,"CalBad_dPhi[11]");
  sprintf(title,"WHA tower (12,10): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[11]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[11]);
  sprintf(name,"CalBad_MetRatio[11]");
  sprintf(title,"WHA tower (12,10): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[11]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[11]);
*/
/*}}}*/

   //___________________________________________________________________________________________________________

  sprintf(name,"CalEm_dPhi");
  sprintf(title,"d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if E_{em}^{twr}>10 GeV");
  Hist.fCalEm_dPhi=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalEm_dPhi);

  sprintf(name,"CalEm_MetRatio");
  sprintf(title,"#slash{E}_{T}/E_{T}^{EM spike} if E_{em}^{twr}>10 GeV");
  Hist.fCalEm_MetRatio=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalEm_MetRatio);
  
  sprintf(name,"CalHad_dPhi");
  sprintf(title,"d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if E_{had}^{twr}>10 GeV");
  Hist.fCalHad_dPhi=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalHad_dPhi);
  
  sprintf(name,"CalHad_MetRatio");
  sprintf(title,"#slash{E}_{T}/E_{T}^{HAD spike} if E_{had}^{twr}>10 GeV");
  Hist.fCalHad_MetRatio=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalHad_MetRatio);


/*  sprintf(name,"CalPmtEmAssm_bad");
  sprintf(title,"pmt asymmetry for bad EM towers: (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2})");
  Hist.fCalPmtEmAssm_bad=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtEmAssm_bad);
  sprintf(name,"CalPmtHadAssm_bad");
  sprintf(title,"pmt asymmetry for bad HAD towers: (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2})");
  Hist.fCalPmtHadAssm_bad=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_bad);
*/

  sprintf(name,"CalPmtEmAssm_le10");
  sprintf(title,"pmt asymmetry: (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) if E_{em}^{twr}<10 GeV");
  Hist.fCalPmtEmAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtEmAssm_le10);
  sprintf(name,"CalPmtHadAssm_le10");
  sprintf(title,"pmt asymmetry: (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) if E_{had}^{twr}<10 GeV");
  Hist.fCalPmtHadAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_le10);
  sprintf(name,"CalPmtEmAssm_ge10");
  sprintf(title,"EM pmt asymmetry: (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) if E_{em}^{twr}>10 GeV");
  Hist.fCalPmtEmAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtEmAssm_ge10);
  sprintf(name,"CalPmtHadAssm_ge10");
  sprintf(title,"HAD pmt asymmetry: (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) if E_{had}^{twr}>10 GeV");
  Hist.fCalPmtHadAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_ge10);

  sprintf(name,"CalPmtCemAssm_le10");
  sprintf(title,"CEM pmt asymmetry: (E_{cem}^{pmt1}-E_{cem}^{pmt2})/(E_{cem}^{pmt1}+E_{cem}^{pmt2}) if E_{cem}^{twr}<10 GeV");
  Hist.fCalPmtCemAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtCemAssm_le10);
  sprintf(name,"CalPmtChaAssm_le10");
  sprintf(title,"CHA pmt asymmetry outside overlap: (E_{cha}^{pmt1}-E_{cha}^{pmt2})/(E_{cha}^{pmt1}+E_{cha}^{pmt2}) if E_{cha}^{twr}<10 GeV");
  Hist.fCalPmtChaAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtChaAssm_le10);
  sprintf(name,"CalPmtWhaAssm_le10");
  sprintf(title,"WHA pmt asymmetry outside overlap: (E_{wha}^{pmt1}-E_{wha}^{pmt2})/(E_{wha}^{pmt1}+E_{wha}^{pmt2}) if E_{wha}^{twr}<10 GeV");
  Hist.fCalPmtWhaAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtWhaAssm_le10);

  sprintf(name,"CalPmtCemAssm_ge10");
  sprintf(title,"CEM pmt asymmetry: (E_{cem}^{pmt1}-E_{cem}^{pmt2})/(E_{cem}^{pmt1}+E_{cem}^{pmt2}) if E_{cem}^{twr}>10 GeV");
  Hist.fCalPmtCemAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtCemAssm_ge10);
  sprintf(name,"CalPmtChaAssm_ge10");
  sprintf(title,"CHA pmt asymmetry outside overlap: (E_{cha}^{pmt1}-E_{cha}^{pmt2})/(E_{cha}^{pmt1}+E_{cha}^{pmt2}) if E_{cha}^{twr}>10 GeV");
  Hist.fCalPmtChaAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtChaAssm_ge10);
  sprintf(name,"CalPmtWhaAssm_ge10");
  sprintf(title,"WHA pmt asymmetry outside overlap: (E_{wha}^{pmt1}-E_{wha}^{pmt2})/(E_{wha}^{pmt1}+E_{wha}^{pmt2}) if E_{wha}^{twr}>10 GeV");
  Hist.fCalPmtWhaAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtWhaAssm_ge10);

  sprintf(name,"CalPmtChaAssm_overlap_le10");
  sprintf(title,"CHA pmt asymmetry in overlap region: (E_{cha}^{pmt1}-E_{cha}^{pmt2})/(E_{cha}^{pmt1}+E_{cha}^{pmt2}) if E_{had}^{twr}<10 GeV");
  Hist.fCalPmtChaAssm_overlap_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtChaAssm_overlap_le10);
  sprintf(name,"CalPmtWhaAssm_overlap_le10");
  sprintf(title,"WHA pmt asymmetry in overlap region: (E_{wha}^{pmt1}-E_{wha}^{pmt2})/(E_{wha}^{pmt1}+E_{wha}^{pmt2}) if E_{had}^{twr}<10 GeV");
  Hist.fCalPmtWhaAssm_overlap_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtWhaAssm_overlap_le10);
  sprintf(name,"CalPmtChaAssm_overlap_ge10");
  sprintf(title,"CHA pmt asymmetry in overlap region: (E_{cha}^{pmt1}-E_{cha}^{pmt2})/(E_{cha}^{pmt1}+E_{cha}^{pmt2}) if E_{had}^{twr}>10 GeV");
  Hist.fCalPmtChaAssm_overlap_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtChaAssm_overlap_ge10);
  sprintf(name,"CalPmtWhaAssm_overlap_ge10");
  sprintf(title,"WHA pmt asymmetry in overlap region: (E_{wha}^{pmt1}-E_{wha}^{pmt2})/(E_{wha}^{pmt1}+E_{wha}^{pmt2}) if E_{had}^{twr}>10 GeV");
  Hist.fCalPmtWhaAssm_overlap_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtWhaAssm_overlap_ge10);

  sprintf(name,"CalPmtHadAssm_ge10_intime");
  sprintf(title,"HAD pmt asymmetry: (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) if E_{had}^{twr}>10 GeV and |HadTdc|<40 ns");
  Hist.fCalPmtHadAssm_ge10_intime=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_ge10_intime);
  sprintf(name,"CalPmtHadAssm_ge10_outtime");
  sprintf(title,"HAD pmt asymmetry: (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) if E_{had}^{twr}>10 GeV and 40<|HadTdc|<80 ns");
  Hist.fCalPmtHadAssm_ge10_outtime=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_ge10_outtime);

  sprintf(name,"CalEmPmtAssm_vs_PmtSum");
  sprintf(title,"CEM Tower Energy vs CEM Tower PMT Asymmetry");
  Hist.fCalEmPmtAssm_vs_PmtSum=new TH2F(name,title,100,-1.0,1.0,300,0.0,300.0);
  new_folder->Add(Hist.fCalEmPmtAssm_vs_PmtSum);
  
  sprintf(name,"CalHadPmtAssm_vs_PmtSum");
  sprintf(title,"Had (CHA & WHA) Tower Energy vs Had Tower PMT Asymmetry");
  Hist.fCalHadPmtAssm_vs_PmtSum=new TH2F(name,title,100,-1.0,1.0,300,0.0,300.0);
  new_folder->Add(Hist.fCalHadPmtAssm_vs_PmtSum);

  return;
}


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TagPMTSpikes::EndJob() {
	
	if (GetSummaryStat()) return 0;
	
	std::string sMsg;
	if (GetMode() == 0) sMsg = "Tag Only and Tag all photons.";
	if (GetMode() == 1) sMsg = "(Reject PMT spike events)";
	if (GetMode() == 2) sMsg = "(Tag all and pass only if EVENT has PMT spikes)";
	
	printf("[PMT:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[PMT:01:] Events Processed------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "[PMT:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;
	if (fHeaderBlock->McFlag())
	{
		std::cout << "[PMT:03:] NO TAGGING DONE FOR MC ****" << std::endl;
	} else {
		std::cout << "[PMT:03:] Mode used ------------------ = " << GetMode() << " (" << sMsg << ")" << std::endl;
		std::cout << "[PMT:04:] Spike Events Found --------- = " << counter.spikeEvts << std::endl;
	}

	printf("---------------------------------------------------\n");
	return 0;
}
