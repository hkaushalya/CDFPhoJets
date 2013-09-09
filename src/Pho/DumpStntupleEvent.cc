///////////////////////////////////////////////////////////
// I am writing this to dump em obj and jets info        //
// (towers etc)from Stntuples. There are some events     //
// where an em obj is not matched to a jet.              //
// see elog: https://hep.baylor.edu/elog/samantha/535    //
// 05-07-2008                                            //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#include "samantha/Pho/DumpStntupleEvent.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnTrack.hh"
#include "Stntuple/obj/TStnJet.hh"
#include "Stntuple/data/TCalTower.hh"
#include <iostream>

ClassImp(DumpStntupleEvent)

//_____________________________________________________________________________
DumpStntupleEvent::DumpStntupleEvent(const char* name, const char* title):
TStnModule(name,title),
bRunPermit(true)
{
	std::cout << "Hello I am DumpStntupleEvent module" << std::endl;
}

//_____________________________________________________________________________
DumpStntupleEvent::~DumpStntupleEvent() {
}

//_____________________________________________________________________________
int DumpStntupleEvent::BeginJob()
{
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*) RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
  	RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlockClu04);
	RegisterDataBlock("CalDataBlock","TCalDataBlock",&fCalDataBlock);
	fTrackBlock 	= (TStnTrackBlock*) RegisterDataBlock("TrackBlock","TStnTrackBlock");
	RegisterDataBlock("PhotonBlock","TStnPhotonBlock",&fPhotonBlock);
	RegisterDataBlock("ElectronBlock","TStnElectronBlock",&fElectronBlock);


  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}

  	jetMod = (JetFilterModuleV2*) ((TStnAna*) GetAna()->GetModule("JetFilterV2"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		bRunPermit = false;
	}


	// initialize vars
	counter.evtsProcessed		= 0;
	counter.evtsPassModule 		= 0;

	return 0;
}

//_____________________________________________________________________________
int DumpStntupleEvent::BeginRun()
{

	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int DumpStntupleEvent::Event(int ientry)
{
	SetPassed(0);
	counter.evtsProcessed++;
	FillDataBlocks(ientry);
	
	std::cout << " ========================================================= " << std::endl;
	std::cout << __FILE__ << std::endl;
	std::cout << "\t" <<  fHeaderBlock->RunNumber() << ", "
				 << fHeaderBlock->EventNumber() << std::endl;

/*	if (fPhotonBlock->NPhotons() > 0) {
		std::cout << " ====== Photon Block Info: NPho = " << fPhotonBlock->NPhotons() << std::endl;
		std::cout << "\t\tEta\tphi\tEt" << std::endl;

		for (int i=0; i < fPhotonBlock->NPhotons(); ++i) DumpPhotonTowers(i);
	}
*/

	std::cout << "\t\t Photon Block Info: NPho = " << fPhotonBlock->NPhotons() << std::endl;

	for (int i=0; i < initSpMod->GetSuperPhoSize(); i++) {
		if (initSpMod->GetSuperPhoton(i)->IsTightPhoton()) 
				std::cout << "\t\t\tphoton " << initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex() << " is tight photon." << std::endl;
		if (initSpMod->GetSuperPhoton(i)->IsLoosePhoton()) 
				std::cout << "\t\t\tphoton " << initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex() << " is loose photon." << std::endl;
		if (initSpMod->GetSuperPhoton(i)->IsTightElectron()) 
				std::cout << "\t\t\tphoton " << initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex() << " is tight electron." << std::endl;
		if (initSpMod->GetSuperPhoton(i)->IsLooseElectron()) 
				std::cout << "\t\t\tphoton " << initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex() << " is loose electron.\n" << std::endl;
	
		if (initSpMod->GetSuperPhoton(i)->IsTightPhoton() || 
			 initSpMod->GetSuperPhoton(i)->IsLoosePhoton() ||
			 initSpMod->GetSuperPhoton(i)->IsTightElectron() || 
			 initSpMod->GetSuperPhoton(i)->IsLooseElectron() ) 
			DumpPhotonTowers(initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex());
	}	


	DumpJetMatchStuff();

//std::cout << " ====== Electron Info ==================================== " << std::endl;
	//for (int i=0; i < fElectronBlock->NElectrons(); ++i) DumpElectronTowers(i);
	
	if (fJetBlockClu04->NJets() > 0) {
		std::cout << " ====== Jet Info : NJet = " << fJetBlockClu04->NJets() << std::endl;
		std::cout << "\t\tEta\tphi\tEt" << std::endl;
		for (int i=0; i < fJetBlockClu04->NJets(); ++i) DumpJetTowers(i);
	}

	std::cout << " ========================================================= \n\n" << std::endl;

	
	return 0;
} // Event




void DumpStntupleEvent::DumpJetMatchStuff()
{
	//jetMod->DumpMatchStuff();
	//jetMod->DumpSamsMatchStuff();
}




//________ creates a list of towers for all jets
void DumpStntupleEvent::DumpJetTowers(int jet_ind)
{
  //CalDataArray* cdholder;
  TStnLinkBlock* links = fJetBlockClu04->TowerLinkList();
  int nptow = links->NLinks(jet_ind);

  TCalTower* ctower = NULL;

	std::cout << "\tInd="<< jet_ind << std::endl;
  
  for(int j=0; j<nptow; j++) 
    {
      int iptow = links->Index(jet_ind,j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;
      
      ctower = fCalDataBlock->Tower(ieta,iphi);
      //cdholder->push_back(ctower);

		//std::cout << "\t\t" << ieta << "\t" << iphi  << "\t" << ctower->EmEnergy() << "\t" << ctower->HadEnergy() << std::endl;
		std::cout << "\t\t" << ctower->Eta() << "\t" << ctower->Phi()  << "\t" << ctower->Et() << std::endl;
    }

  //std::sort(cdholder->begin(), cdholder->end(),SortTowersByEnergy);
  
  return;
}

//________ creates a list of towers for all EM objects
void DumpStntupleEvent::DumpPhotonTowers(int pho_ind)
{
  //CalDataArray* cdholder;
  TStnLinkBlock* links = fPhotonBlock->TowerLinkList();
  int nptow = links->NLinks(pho_ind);
  TCalTower* ctower = NULL;

	std::cout << "\tInd="<< pho_ind << std::endl;
  
  for(int j=0; j<nptow; j++) 
    {
      int iptow = links->Index(pho_ind,j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;
      
      ctower = fCalDataBlock->Tower(ieta,iphi);
      //cdholder->push_back(ctower);

		//std::cout << "\t\t" << ieta << "\t" << iphi  << "\t" << ctower->EmEnergy() << "\t" << ctower->HadEnergy() << std::endl;
		std::cout << "\t\t" << ctower->Eta() << "\t" << ctower->Phi()  << "\t" << ctower->Et() << std::endl;
		
    }
  
 // std::sort(cdholder->begin(),cdholder->end(),SortTowersByEnergy);
  
  
  return;
}

//________ creates a list of towers for electrons
void DumpStntupleEvent::DumpElectronTowers(int ele_ind)
{
  //CalDataArray* cdholder;
  TStnLinkBlock* links = fElectronBlock->TowerLinkList();
  int nptow = links->NLinks(ele_ind);
  //TCalTower* ctower = NULL;

  
  for(int j=0; j<nptow; j++) 
    {
      int iptow = links->Index(ele_ind,j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;
      
      //ctower = fCalDataBlock->Tower(ieta,iphi);
      //cdholder->push_back(ctower);

		std::cout << "\t\t" << ieta << "\t" << iphi << std::endl;
    }
  
 // std::sort(cdholder->begin(),cdholder->end(),SortTowersByEnergy);
  
  return;
}


/*-------------------------------------------------------------------*/
void DumpStntupleEvent::FillDataBlocks(int ientry)
{
  fJetBlockClu04->GetEntry(ientry);
  fPhotonBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  fElectronBlock->GetEntry(ientry);
  fTrackBlock->GetEntry(ientry);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int DumpStntupleEvent::EndJob() {

	printf("[DSE:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[DSE:01:] Events Processed ----------- = " << counter.evtsProcessed << std::endl;
	std::cout << "[DSE:02:] Events Pass this module ---- = " << counter.evtsPassModule << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
