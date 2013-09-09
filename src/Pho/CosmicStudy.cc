///////////////////////////////////////////////////////////
// This is to study the Cosmic events. I need to reduce  //
// effect of cosmics in the high MET tail.               //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/*
 *	$Log: CosmicStudy.cc,v $
 *	Revision 1.3  2009/06/08 03:29:55  samantha
 *	ADDED: Auto CVS log in file.
 *	
 *
 */

#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <numeric>
#include <vector>
#include <algorithm>
#include "samantha/Pho/CosmicStudy.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "Stntuple/obj/TStnVertex.hh"

ClassImp(CosmicStudy)

//_____________________________________________________________________________
CosmicStudy::CosmicStudy(const char* name, const char* title):
  TStnModule(name,title),
  bPermitRun(true),
  fCosmicTimeMin(30),
  fCosmicTimeMax(90)
{
	std::cout << "Hello I am CosmicStudy module" << std::endl;
}

//_____________________________________________________________________________
CosmicStudy::~CosmicStudy() {
}

//_____________________________________________________________________________
void CosmicStudy::SaveHistograms() {
}

//_____________________________________________________________________________
void CosmicStudy::BookHistograms()
{
	DeleteHistograms();


	hPho1Jet_Jet1EmTime = new TH1F("Jet1EmTime","#gamma+>=1 Jet : Lead Jet EmTime",800,-50,150);
	hPho1Jet_PhoJet1DelEmTime = new TH1F("PhoJet1DelEmTime","#gamma+>=1 Jet : #Delta t (#gamma,Lead Jet)",150,0,300);
	hPho2Jet_Jet1EmTime = new TH1F("Jet1EmTime","#gamma+>=2 Jets : Lead Jet EmTime",800,-50,150);
	hPho2Jet_Jet2EmTime = new TH1F("Jet2EmTime","#gamma+>2 Jets: 2nd Lead Jet EmTime",800,-50,150);
	hPho2Jet_PhoJet1DelEmTime = new TH1F("PhoJet1DelEmTime","#gamma+>=2 Jet : #Delta t (#gamma,Lead Jet)",150,0,300);
	hPho2Jet_PhoJet2DelEmTime = new TH1F("PhoJet2DelEmTime","#gamma+>=2 Jet : #Delta t (#gamma,Second Lead Jet)",150,0,300);

  	//char name [200];
  	//char title[200];
	TFolder *new_folder, *sub_folder;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");
	vCounterHistLabels.push_back("Events Passed");
	vCounterHistLabels.push_back("Cosmic #gamma +>=1 Jets Events");
	vCounterHistLabels.push_back("Cosmic #gamma +>=2 Jets Events");
	vCounterHistLabels.push_back("EM Cos. w. stubs (1>=Jets)");
	vCounterHistLabels.push_back("EM Cos. w. stubs (2>=Jets)");
	HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);

	new_folder = GetHistFolder(this, "1Jet","#gamma + >=1 Jet case");
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt,"#gamma + >=1 Jets(15GeV)");
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h1_Pho,"#gamma + >=1 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet,"#gamma + >=1 Jets(15GeV)");
	sub_folder->Add(hPho1Jet_Jet1EmTime);
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	sub_folder->Add(hPho1Jet_PhoJet1DelEmTime);
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet,"#gamma + >=1 Jets(15GeV)");



	new_folder = GetHistFolder(this, "2Jet","#gamma + >=2 Jets case");
	
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h2_Evt,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h2_Pho,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet1,"#gamma + >=2 Jets(15GeV)");
	sub_folder->Add(hPho2Jet_Jet1EmTime);
	
	sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet2,"#gamma + >=2 Jets(15GeV)");
	sub_folder->Add(hPho2Jet_Jet2EmTime);
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	sub_folder->Add(hPho2Jet_PhoJet1DelEmTime);
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet1,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
	sub_folder->Add(hPho2Jet_PhoJet2DelEmTime);
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet2,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
	HistoMan.GetPhoton2JetsHistograms(sub_folder,h2_PhoJets,"#gamma + >=2 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
	HistoMan.GetTwoJetsHistograms(sub_folder,h2_TwoJets,"#gamma + >=2 Jets(15GeV)");
	
}


//_____________________________________________________________________________
int CosmicStudy::BeginJob()
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
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");

	BookHistograms();

	// initialize vars
	
	
	return 0;
}

//_____________________________________________________________________________
int CosmicStudy::BeginRun()
{
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int CosmicStudy::Event(int ientry)
{
	SetPassed(1);

	// this makes sure that all required mods exists
	if (!bPermitRun) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
	}
	
	hCount.iCounter->Fill(0);

	FillDataBlocks(ientry);


//	if (NClass12Vertices() < 1) return 1;		//don't need this. can control this from TriggerModule
//	if (BestVertexZ() > 60) return 1;

	int Run = fHeaderBlock->RunNumber();
	
	if (Run>=190851) { // exclude first 400pb-1 with out em timing
		
		int Nsp = initSpMod->GetSuperPhoSize();
		std::cout << "Nsp= " << Nsp << std::endl;
		
		for (int i = 0; i < Nsp; i++) {
			if (! initSpMod->GetSuperPhoton(i)->IsTightPhoton()) continue;
			int iHaloId = initSpMod->GetSuperPhoton(i)->GetBeamHaloId();
			//remove halo for a given type. must be same type as used for 
			// halo template!!
			if ( ( (iHaloId >> 0) & 0x1) ) continue;
			
			int iNCosStubs = initSpMod->GetSuperPhoton(i)->GetPhoton()->NCosStubPho();  // look for any muon stubs with the 30 degree cone of the photon
			float fTime = initSpMod->GetSuperPhoton(i)->GetEmTime();

			std::cout << "ftime = " << fTime << std::endl;
			if ((fTime > GetCosmicTimeMin() && fTime < GetCosmicTimeMax())) {
				DoJetSelection(i);
				//count how many cosmics from em time has a trkless muon stub
				if (iNCosStubs>0) {
					if (jetMod->GetMyNjet(0,3) >= 1) hCount.iCounter->Fill(4);
					if (jetMod->GetMyNjet(0,3) >= 2)	hCount.iCounter->Fill(5);
				}
			}
			break; // don't care about rest of the photons
					//if not i'll be putting multiple entries into met.sumet plots
		}
	}

	if (GetPassed()) {
		hCount.iCounter->Fill(1);
	}
	return 0;

} // Event


/*-------------------------------------------------------------------*/
void CosmicStudy::DoJetSelection(int iSpIndex)
{
	//float fWht = 1;
	//int iLeadJetIndex = 0;

	float fPhoTime = initSpMod->GetSuperPhoton(iSpIndex)->GetEmTime();

	if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
		hCount.iCounter->Fill(2);
		evtPropMod->FillPhotonHists(h1_Pho,iSpIndex);
		evtPropMod->FillEventHists(h1_Evt);
		evtPropMod->FillJetHists(h1_Jet,0);		//lead jet
		evtPropMod->FillPhoton1JetHists(h1_PhoJet,iSpIndex,0);
		float jet1_emtime = GetEmTiming(0);
		hPho1Jet_Jet1EmTime->Fill(jet1_emtime);
		hPho1Jet_PhoJet1DelEmTime->Fill(fabs(jet1_emtime - fPhoTime));
		
			
	}
	
	if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
		hCount.iCounter->Fill(3);
		evtPropMod->FillPhotonHists(h2_Pho,iSpIndex);
		evtPropMod->FillEventHists(h2_Evt);
		evtPropMod->FillJetHists(h2_Jet1,0);		//lead jet
		evtPropMod->FillJetHists(h2_Jet2,1);		//2nd lead jet
		evtPropMod->FillPhoton1JetHists(h2_PhoJet1,iSpIndex,0);
		evtPropMod->FillPhoton1JetHists(h2_PhoJet2,iSpIndex,1);
		evtPropMod->FillPhoton2JetsHists(h2_PhoJets,iSpIndex,0,1);
		evtPropMod->FillTwoJetsHists(h2_TwoJets,0,1);
		float jet1_emtime = GetEmTiming(0);
		float jet2_emtime = GetEmTiming(1);
		hPho2Jet_Jet1EmTime->Fill(jet1_emtime);
		hPho2Jet_Jet2EmTime->Fill(jet2_emtime);
		hPho2Jet_PhoJet1DelEmTime->Fill(fabs(jet1_emtime - fPhoTime));
		hPho2Jet_PhoJet2DelEmTime->Fill(fabs(jet2_emtime - fPhoTime));
	}

}


/*-------------------------------------------------------------------*/
// runN>=190851  first run for EM timing system
/*-------------------------------------------------------------------*/
double CosmicStudy::GetEmTiming(const int jetInd)
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
void CosmicStudy::MatchCalorTowers(int jet_ind, CalDataArray* cdholder)
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
// z-position of BestClass12 vertex
/*-------------------------------------------------------------------*/
float CosmicStudy::BestVertexZ()
{
  float zvx =0;
  if (fVertexBlock->GetBestVertex(12,1) != NULL) {
    zvx = fVertexBlock->GetBestVertex(12,1)->Z();
  }
  
  return zvx;
} //BestVertex_Z

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int CosmicStudy::NClass12Vertices()
{
 	int iN12Vertices = 0;
	int iNVertices = fVertexBlock->NVertices();

	for (int iVtx = 0; iVtx < iNVertices; ++iVtx) {
		TStnVertex* vert = fVertexBlock->Vertex(iVtx);
		if (vert->VClass() >= 12) ++iN12Vertices;
	}
  
  return iN12Vertices;

}  // Class12Vertices





/*-------------------------------------------------------------------*/
void CosmicStudy::FillDataBlocks(int ientry)
{
	fCalBlock->GetEntry(ientry);
	fJetBlock->GetEntry(ientry);
	fVertexBlock->GetEntry(ientry);
}
//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int CosmicStudy::EndJob() {

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[CSD:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[CSD:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[CSD:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	sMsg += "[CSD:04:] Cosmic Pho + >=1 Jets Events = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
	sMsg += "[CSD:05:] Cosmic Pho + >=2 Jets Events = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n";
	sMsg += "[CSD:06:] Em Cos. w. stubs (+>=1 Jets) = " + ToStr(hCount.iCounter->GetBinContent(5)) + "\n";
	sMsg += "[CSD:07:] Em Cos. w. stubs (+>=2 Jets) = " + ToStr(hCount.iCounter->GetBinContent(6)) + "\n";
	sMsg += "[CSD:08:] EM time window used  ------- = " + ToStr(GetCosmicTimeMin()) + "," +
																			ToStr(GetCosmicTimeMax()) + "\n";
	
	sMsg += "---------------------------------------------------\n";
	std::cout << sMsg;
	return 0;
}
