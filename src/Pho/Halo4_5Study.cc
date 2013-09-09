///////////////////////////////////////////////////////////
// this is to study halos from scenario 4 and 5 with     //
// Met/PhoEt > 1 and <1 to figure out if a halo can      //
// produce a jet not only a photons.                     //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "samantha/Pho/Halo4_5Study.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>

ClassImp(Halo4_5Study)

//_____________________________________________________________________________
Halo4_5Study::Halo4_5Study(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bRunPermit(true),
  iHaloScenario(0),
  fEmTimeMin(-4.8),
  fEmTimeMax(4.8)
{
	// Default constructor
	std::cout << "Hello I am Halo4_5Study module" << std::endl;
}

//_____________________________________________________________________________
Halo4_5Study::~Halo4_5Study()
{
	// Default destructor
}

//_____________________________________________________________________________
void Halo4_5Study::SaveHistograms()
{
}

//_____________________________________________________________________________
void Halo4_5Study::BookHistograms()
{
	// Book all the histograms
	
	DeleteHistograms();
	
  //	char name [200];
  	//char title[200];
	TFolder *folder, *new_folder, *sub_folder;
	Histograms HistoMan;

	folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");
	vCounterHistLabels.push_back("Events Passed");
	vCounterHistLabels.push_back("Met/Pho Et>1");
	vCounterHistLabels.push_back("Met/Pho Et<1");
	HistoMan.GetCounterHistograms(folder,hCount,this->GetName(),vCounterHistLabels);

	folder = GetHistFolder(this, "MetEtRatio>1","Met/Photon Et Ratio > 1");
	
	new_folder = folder->AddFolder("1Jet","#gamma + >=1 Jet case");

	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt,"#gamma + >=1 Jets(15GeV)");
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h1_Pho,"#gamma + >=1 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet,"#gamma + >=1 Jets(15GeV)");
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet,"#gamma + >=1 Jets(15GeV)");


	new_folder = folder->AddFolder("2Jet","#gamma + >=2 Jet case");
	
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h2_Evt,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h2_Pho,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet1,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet2,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet1,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet2,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
	HistoMan.GetPhoton2JetsHistograms(sub_folder,h2_PhoJets,"#gamma + >=2 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
	HistoMan.GetTwoJetsHistograms(sub_folder,h2_TwoJets,"#gamma + >=2 Jets(15GeV)");


	folder = GetHistFolder(this, "MetEtRatio<1","Met/Photon Et Ratio < 1");
	
	new_folder = folder->AddFolder("1Jet","#gamma + >=1 Jet case");

	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt_,"#gamma + >=1 Jets(15GeV)");
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h1_Pho_,"#gamma + >=1 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet_,"#gamma + >=1 Jets(15GeV)");
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet_,"#gamma + >=1 Jets(15GeV)");


	new_folder = folder->AddFolder("2Jet","#gamma + >=2 Jet case");
	
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h2_Evt_,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h2_Pho_,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet1_,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet2_,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet1_,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet2_,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
	HistoMan.GetPhoton2JetsHistograms(sub_folder,h2_PhoJets_,"#gamma + >=2 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
	HistoMan.GetTwoJetsHistograms(sub_folder,h2_TwoJets_,"#gamma + >=2 Jets(15GeV)");



}


//_____________________________________________________________________________
int Halo4_5Study::BeginJob()
{
	// Initialize and check for all dependecies
	
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

  	tightMod = (TagTightPhotons*) ((TStnAna*) GetAna()->GetModule("TagTightPhotons"));
	if (!tightMod) {
		StdOut(__FILE__,__LINE__,2,"Tight photon tag module not found.");
		bRunPermit = false;
	}

  	tagHaloMod = (TagBeamHalo*) ((TStnAna*) GetAna()->GetModule("TagBeamHalo"));
	if (!tagHaloMod) {
		StdOut(__FILE__,__LINE__,2,"Halo tag module not found.");
		bRunPermit = false;
	}
	
  	evtPropMod = (EventProperties*) ((TStnAna*) GetAna()->GetModule("EventProperties"));
	if (!evtPropMod) {
		StdOut(__FILE__,__LINE__,2,"EventProperties tag module not found.");
		bRunPermit = false;
	}

  	setTimeMod = (SetEmTimes*) ((TStnAna*) GetAna()->GetModule("SetEmTimes"));
	if (!setTimeMod) {
		StdOut(__FILE__,__LINE__,2,"SetEmTimes tag module not found.");
		bRunPermit = false;
	}


				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.haloEventFound = 0;

	return 0;
}

//_____________________________________________________________________________
int Halo4_5Study::BeginRun()
{

	//int currRun =  fHeaderBlock->RunNumber();
	//std::cout << " BEGINING RUN# " << currRun << std::endl;
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int Halo4_5Study::Event(int ientry)
{
	// Event loop
	
	SetPassed(1);
	
	// this makes sure that all required mods exists
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
		return 0;
	}
		
	counter.evtsRunOver++;
	hCount.iCounter->Fill(0);
  	int NsuperPho = initSpMod->GetSuperPhoSize();
	
	for (int i=0; i < NsuperPho; i++) {
		if (! initSpMod->GetSuperPhoton(i)->IsTightPhoton()) continue;
		if (! initSpMod->GetSuperPhoton(i)->IsBeamHalo()) continue;		// select BH
			
		float fEmTime = initSpMod->GetSuperPhoton(i)->GetEmTime();
		if (fEmTime > fEmTimeMax || fEmTime < fEmTimeMin) continue;		//cut on em time, same as real photon template

		int iHaloId = initSpMod->GetSuperPhoton(i)->GetBeamHaloId();
		if ( ! ((iHaloId >> GetHaloScenario()) & 0x1) ) continue;   //ignore if it is not the scenario picked

		float MetEtRatio =  jetMod->GetMyMetCorr(0,3)/ initSpMod->GetSuperPhoton(i)->GetEtCorr();

		if (MetEtRatio >= 1.0) {
				hCount.iCounter->Fill(2);

			if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
				evtPropMod->FillPhotonHists(h1_Pho,i,1);
				evtPropMod->FillEventHists(h1_Evt,1);
				evtPropMod->FillJetHists(h1_Jet,0,1);		//lead jet
				evtPropMod->FillPhoton1JetHists(h1_PhoJet,i,0,1);
			}
			
			if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
				evtPropMod->FillPhotonHists(h2_Pho,i,1);
				evtPropMod->FillEventHists(h2_Evt,1);
				evtPropMod->FillJetHists(h2_Jet1,0,1);		//lead jet
				evtPropMod->FillJetHists(h2_Jet2,1,1);		//2nd lead jet
				evtPropMod->FillPhoton1JetHists(h2_PhoJet1,i,0,1);
				evtPropMod->FillPhoton1JetHists(h2_PhoJet2,i,1,1);
				evtPropMod->FillPhoton2JetsHists(h2_PhoJets,i,0,1,1);
				evtPropMod->FillTwoJetsHists(h2_TwoJets,0,1,1);
			}

		} else {
				
				hCount.iCounter->Fill(3);

			if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
				evtPropMod->FillPhotonHists(h1_Pho_,i,1);
				evtPropMod->FillEventHists(h1_Evt_,1);
				evtPropMod->FillJetHists(h1_Jet_,0,1);		//lead jet
				evtPropMod->FillPhoton1JetHists(h1_PhoJet_,i,0,1);
			}
			
			if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
				evtPropMod->FillPhotonHists(h2_Pho_,i,1);
				evtPropMod->FillEventHists(h2_Evt_,1);
				evtPropMod->FillJetHists(h2_Jet1_,0,1);		//lead jet
				evtPropMod->FillJetHists(h2_Jet2_,1,1);		//2nd lead jet
				evtPropMod->FillPhoton1JetHists(h2_PhoJet1_,i,0,1);
				evtPropMod->FillPhoton1JetHists(h2_PhoJet2_,i,1,1);
				evtPropMod->FillPhoton2JetsHists(h2_PhoJets_,i,0,1,1);
				evtPropMod->FillTwoJetsHists(h2_TwoJets_,0,1,1);
			}




		}

			break; // don't care about rest of the photons
					//if not i'll be putting multiple entries into met.sumet plots
	}

	if (GetPassed()) {
		counter.evtsPassModule++;
		hCount.iCounter->Fill(1);
	}
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void Halo4_5Study::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void Halo4_5Study::FillDataBlocks(int ientry)
{
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int Halo4_5Study::EndJob() {

	std::string sModName(GetName());
	std::string sMsg;
	
	sMsg  = "[H45:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[H45:01:] Events Processed ----------- = " + ToStr(counter.evtsRunOver) + "\n";
	sMsg += "[H45:02:] Events Pass this module ---- = " + ToStr(counter.evtsPassModule) + "\n";
	sMsg += "[H45:03:] Scenario Used -------------- = " + ToStr(GetHaloScenario()) + "\n";
	sMsg += "[H45:04:] MEt/PhoEt > 1 -------------- = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
	sMsg += "[H45:05:] MEt/PhoEt < 1 -------------- = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n";

	sMsg += "---------------------------------------------------\n";
	std::cout << sMsg;
	return 0;
}
