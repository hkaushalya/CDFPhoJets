///////////////////////////////////////////////////////////
// This is to study the Jet Energy Corrections.          //
// I want to look at how much of correction is applied   //
// by looking in different Jet Et bins.                  //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/* 
 *$Log: JERStudy.cc,v $
 *Revision 1.6  2009/05/21 22:17:38  samantha
 *Commiting all the minor changes in one.
 *
 *Revision 1.5  2009/05/07 17:10:30  samantha
 *testing cvs editor.
 *
 *Revision 1.4  2009/05/07 16:57:02  samantha
 *testing cvs comment editor.
 *
 *Revision 1.3  2009/05/07 16:55:50  samantha
 *Testing CVS editor for entering log info.
 *
 * Revision 1.2  2009/05/07 16:53:45  samantha
 * *** empty log message ***
 *
 * Revision 1.1  2009/03/28 19:01:47  samantha
 * This is to study the JER/JES by different Jet Et bins. Not complete yet.
 *
 */    

#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <vector>
#include <algorithm>
#include "samantha/Pho/JERStudy.hh"
#include "samantha/utils/FreeFunctions.hh"

ClassImp(JERStudy)

//_____________________________________________________________________________
JERStudy::JERStudy(const char* name, const char* title):
  TStnModule(name,title),
  bPermitRun(true)
{
	std::cout << "Hello I am JERStudy module" << std::endl;
}

//_____________________________________________________________________________
JERStudy::~JERStudy() {
}

//_____________________________________________________________________________
void JERStudy::SaveHistograms() {
}

//_____________________________________________________________________________
void JERStudy::BookHistograms()
{
	DeleteHistograms();

	hJetEmTime = new TH1F("Jet1EmTime","Jet  EM Time",800,-50,150);
	hLev6RawJetEtRatio = new TH1F("L6RawJetEtRatio","L6 Corrected Jet Pt/Raw Jet Pt(>3GeV) - EM objects removed",2000,-100,100);
	hL6RawJetEtRatioVsRawJetEt = new TH2F("L6RawJetEtRatioVsRawJetEt","Jet^{Level 6} Pt/ Jet^{Raw} Pt(>3GeV) VS Jet^{Raw} Pt(3GeV)- EM objects removed",0,300,300,2000,-100,100);

  	//char name [200];
  	//char title[200];
	TFolder *new_folder, *sub_folder;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");
	vCounterHistLabels.push_back("Events Passed");
	HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);

	new_folder = GetHistFolder(this, "1Jet","#gamma + >=1 Jet case");
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt,"#gamma + >=1 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet,"#gamma + >=1 Jets(15GeV)");
	sub_folder->Add(hJetEmTime);
	sub_folder->Add(hLev6RawJetEtRatio);
	sub_folder->Add(hL6RawJetEtRatioVsRawJetEt);
	
}


//_____________________________________________________________________________
int JERStudy::BeginJob()
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
  	evtPropMod = (EventProperties*) ((TStnAna*) GetAna()->GetModule("EventProperties"));
	if (!evtPropMod) {
		StdOut(__FILE__,__LINE__,2,"EventProperties tag module not found.");
		bPermitRun = false;
	}
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fJetBlock 		= (TStnJetBlock*)      RegisterDataBlock("JetBlock","TStnJetBlock");
	RegisterDataBlock("CalDataBlock","TCalDataBlock",&fCalBlock);

	BookHistograms();

	// initialize vars
	
	
	return 0;
}

//_____________________________________________________________________________
int JERStudy::BeginRun()
{

	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int JERStudy::Event(int ientry)
{
	SetPassed(1);

	// this makes sure that all required mods exists
	if (!bPermitRun) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
	}
	
	hCount.iCounter->Fill(0);

	FillDataBlocks(ientry);

	//int Run = fHeaderBlock->RunNumber();
	
	GetHeaderBlock()->Print();
	if (jetMod->GetMyNjet(0,3) >= 1) hCount.iCounter->Fill(4);
	
	for (int i=0; i < jetMod->GetMyNjet(0,0); ++i)
	{
		std::cout << i << "\t" << jetMod->GetMyJet_raw(0,i)->Pt() << "\t" << jetMod->GetMyJet_lev6(0,i)->Pt() << "\t" << jetMod->GetMyJetNoEMobj_raw(0,i)->Pt() << "\t" << jetMod->GetMyJetNoEMobj_lev6(0,i)->Pt() << std::endl;

		float pt_l6 = jetMod->GetMyJetNoEMobj_lev6(0,i)->Pt();
		float  pt_raw = jetMod->GetMyJetNoEMobj_raw(0,i)->Pt(); 
		float pt_ratio = pt_l6/pt_raw;
		//float pt_residual = (pt_l6 - pt_raw)/pt_raw; 
		hLev6RawJetEtRatio->Fill(pt_ratio);
		hL6RawJetEtRatioVsRawJetEt->Fill(pt_raw,pt_ratio);

	}
	

	if (GetPassed()) {
		hCount.iCounter->Fill(1);
	}
	return 0;

} // Event


/*-------------------------------------------------------------------*/
void JERStudy::DoJetSelection(int iSpIndex)
{
	//float fWht = 1;
	//int iLeadJetIndex = 0;

	//float fPhoTime = initSpMod->GetSuperPhoton(iSpIndex)->GetEmTime();

	if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
		hCount.iCounter->Fill(2);
		evtPropMod->FillEventHists(h1_Evt);
		evtPropMod->FillJetHists(h1_Jet,0);		//lead jet
		float jet1_emtime = GetEmTiming(0);
		hJetEmTime->Fill(jet1_emtime);
		
			
	}
}

/*-------------------------------------------------------------------*/
// runN>=190851  first run for EM timing system
/*-------------------------------------------------------------------*/
double JERStudy::GetEmTiming(const int jetInd)
{
	double time = -999999.99;
	if (jetInd < 0) {
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
			float sumtime = 0;
			EmTimingTower* _emTimeTower = NULL;
			CalDataArray calholder;
			
			MatchCalorTowers(jetInd, &calholder);
				
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
						float t=_emTimeTower->getT0(i);
						sumtime+=t;
						//break;
		   	 	}
				}
			}
			
			time = sumtime/(float)NGoodHits;
			
		}

	}
  return time;
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void JERStudy::MatchCalorTowers(int jet_ind, CalDataArray* cdholder)
{
  cdholder->clear();
  TStnLinkBlock* links = fJetBlock->TowerLinkList();
  int nptow = links->NLinks(jet_ind);

  TCalTower* ctower = NULL;
  for(int j=0; j<nptow; j++) 
    {
      int iptow = links->Index(jet_ind,j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;

      ctower = fCalBlock->Tower(ieta,iphi);
      cdholder->push_back(ctower);
    }

  std::sort(cdholder->begin(), cdholder->end(), SortTowersByEnergy);
  return;
}



/*-------------------------------------------------------------------*/
void JERStudy::FillDataBlocks(int ientry)
{
	fCalBlock->GetEntry(ientry);
	fJetBlock->GetEntry(ientry);
}
//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int JERStudy::EndJob() {

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[JER:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[JER:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[JER:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	
	sMsg += "---------------------------------------------------\n";
	std::cout << sMsg;
	return 0;
}
