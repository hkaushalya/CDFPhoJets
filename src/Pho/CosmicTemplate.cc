///////////////////////////////////////////////////////////
// This is the Cosmic template.                          //
// here is the plan for this.                            //
// let cosmic tag module tag them.                       //
// then in here, find tags and count how many cosmic     //
// photons identified by em time has a trkless muon stub //
// to estimate the remaider of cosmic in the signal.!    //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <numeric>
#include <vector>
#include <algorithm>
#include "Stntuple/loop/TStnAna.hh"
#include "samantha/Pho/CosmicTemplate.hh"
#include "samantha/utils/FreeFunctions.hh"

ClassImp(CosmicTemplate)

//_____________________________________________________________________________
CosmicTemplate::CosmicTemplate(const char* name, const char* title):
  TStnModule(name,title),
  bPermitRun(true),
  fCosmicTimeMin(30),
  fCosmicTimeMax(90)
{
	std::cout << "Hello I am CosmicTemplate module" << std::endl;
}

//_____________________________________________________________________________
CosmicTemplate::~CosmicTemplate() {
}

//_____________________________________________________________________________
void CosmicTemplate::SaveHistograms() {
}

//_____________________________________________________________________________
void CosmicTemplate::BookHistograms()
{
	DeleteHistograms();

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
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet,"#gamma + >=1 Jets(15GeV)");



	new_folder = GetHistFolder(this, "2Jet","#gamma + >=2 Jets case");
	
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
	
}


//_____________________________________________________________________________
int CosmicTemplate::BeginJob()
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

	BookHistograms();

	// initialize vars
	
	
	return 0;
}

//_____________________________________________________________________________
int CosmicTemplate::BeginRun()
{

	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int CosmicTemplate::Event(int ientry)
{
	SetPassed(1);

	// this makes sure that all required mods exists
	if (!bPermitRun) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
		return 0;
	}
	
	hCount.iCounter->Fill(0);

	int Run = fHeaderBlock->RunNumber();
	
	if (Run>=190851) { // exclude first 400pb-1 with out em timing
		
		int Nsp = initSpMod->GetSuperPhoSize();
		
		for (int i = 0; i < Nsp; i++) {
			if (! initSpMod->GetSuperPhoton(i)->IsTightPhoton()) continue;
			
			int iNCosStubs = initSpMod->GetSuperPhoton(i)->GetPhoton()->NCosStubPho();  // look for any muon stubs with the 30 degree cone of the photon
			float fTime = initSpMod->GetSuperPhoton(i)->GetEmTime();

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
void CosmicTemplate::DoJetSelection(int iSpIndex)
{
	//float fWht = 1;
	//int iLeadJetIndex = 0;

	if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
		SetPassed(0);			//making exclusive templates
		hCount.iCounter->Fill(2);
		evtPropMod->FillPhotonHists(h1_Pho,iSpIndex);
		evtPropMod->FillEventHists(h1_Evt);
		evtPropMod->FillJetHists(h1_Jet,0);		//lead jet
		evtPropMod->FillPhoton1JetHists(h1_PhoJet,iSpIndex,0);
	}
	
	if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
		SetPassed(0);			//making exclusive templates
		hCount.iCounter->Fill(3);
		evtPropMod->FillPhotonHists(h2_Pho,iSpIndex);
		evtPropMod->FillEventHists(h2_Evt);
		evtPropMod->FillJetHists(h2_Jet1,0);		//lead jet
		evtPropMod->FillJetHists(h2_Jet2,1);		//2nd lead jet
		evtPropMod->FillPhoton1JetHists(h2_PhoJet1,iSpIndex,0);
		evtPropMod->FillPhoton1JetHists(h2_PhoJet2,iSpIndex,1);
		evtPropMod->FillPhoton2JetsHists(h2_PhoJets,iSpIndex,0,1);
		evtPropMod->FillTwoJetsHists(h2_TwoJets,0,1);
	}

}



//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int CosmicTemplate::EndJob() {

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[COT:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[COT:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[COT:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	sMsg += "[COT:04:] Cosmic Pho + >=1 Jets Events = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
	sMsg += "[COT:05:] Cosmic Pho + >=2 Jets Events = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n";
	sMsg += "[COT:06:] Em Cos. w. stubs (+>=1 Jets) = " + ToStr(hCount.iCounter->GetBinContent(5)) + "\n";
	sMsg += "[COT:07:] Em Cos. w. stubs (+>=2 Jets) = " + ToStr(hCount.iCounter->GetBinContent(6)) + "\n";
	sMsg += "[COT:08:] EM time window used  ------- = " + ToStr(GetCosmicTimeMin()) + "," +
																			ToStr(GetCosmicTimeMax()) + "\n";
	
	sMsg += "---------------------------------------------------\n";
	std::cout << sMsg;
	return 0;
}
