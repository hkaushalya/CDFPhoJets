#include "samantha/Pho/PhoJets.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>
#include <iomanip>
#include <Stntuple/obj/TStnJet.hh>
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include <sstream>
#include "TAxis.h"
#include "samantha/utils/FreeFunctions.hh"
#include <assert.h>


ClassImp(PhoJets)


//_____________________________________________________________________________
PhoJets::PhoJets(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bPermitRun(true),
  fEmTimeMin(-4.8),
  fEmTimeMax(4.8),
  iPhotonType(0),
  fPhoMaxDetEta(1.1),
  fPhoMinDetEta(0),
  fPhoMinEt(30),
  bRejPhxPhos(true),
  iHaloType(5),
  fJetMinEt(15),
  fJetMaxDetEta(3.0),
  fMinMEt(0.0),
  bNoSummary(false)
{
	std::cout << "Hello I am PhoJets module" << std::endl;
}

//_____________________________________________________________________________
PhoJets::~PhoJets() {
}

//_____________________________________________________________________________
void PhoJets::SaveHistograms() {
}

//_____________________________________________________________________________
void PhoJets::BookHistograms()
{
	DeleteHistograms();

	TFolder *new_folder, *sub_folder;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");				// 0
	vCounterHistLabels.push_back("Events Passed");					// 1
	vCounterHistLabels.push_back("All Tight #gamma Events");		// 2
	vCounterHistLabels.push_back("Tight #gamma+#geq 1 Events");		// 3
	vCounterHistLabels.push_back("Tight #gamma+#geq 2 Events");		// 4
	vCounterHistLabels.push_back("One tight #gamma only Events");	// 5
	vCounterHistLabels.push_back("Two tight #gamma only events");	// 6
	vCounterHistLabels.push_back("One tight #gamma and any Nloose pho events");	// 7
	HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	

	new_folder = GetHistFolder(this, "1Jet","#gamma + >=1 Jet case");
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt,"#gamma + >=1 Jets(15GeV)");

	metb4 = dynamic_cast<TH1F*> (h1_Evt.Met->Clone("Met_b4cut"));
	meta4 = dynamic_cast<TH1F*> (h1_Evt.Met->Clone("Met_a4cut"));
	sub_folder->Add(metb4);
	sub_folder->Add(meta4);
	
	
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h1_Pho,"#gamma + >=1 Jets(15GeV)");
	
	phometdelphi_a4 = dynamic_cast<TH1F*> (h1_Pho.PhoMetDelPhi->Clone("PhoMetDelPhi_a4"));
	sub_folder->Add(phometdelphi_a4);

	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet,"#gamma + >=1 Jets(15GeV)");
	

	jetmetdelphi_a4 = dynamic_cast<TH1F*> (h1_Jet.JetMetDelPhi->Clone("JetMetDelPhi_a4"));
	sub_folder->Add(jetmetdelphi_a4);
	
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
int PhoJets::BeginJob()
{
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bPermitRun = false;
	}

	jetMod = (JetFilterModuleV3*) ((TStnAna*) GetAna()->GetModule("JetFilterV3"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		bPermitRun = false;
	}

  	tightMod = (TagTightPhotons*) ((TStnAna*) GetAna()->GetModule("TagTightPhotons"));
	if (!tightMod) {
		StdOut(__FILE__,__LINE__,2,"Tight photon tag module not found.");
		bPermitRun = false;
	}

  	evtPropMod = (EventProperties*) ((TStnAna*) GetAna()->GetModule("EventProperties"));
	if (!evtPropMod) {
		StdOut(__FILE__,__LINE__,2,"EventProperties tag module not found.");
		bPermitRun = false;
	}

  	pmtMod = (TagPMTSpikes*) ((TStnAna*) GetAna()->GetModule("TagPMTSpikes"));
	if (!pmtMod) {
		StdOut(__FILE__,__LINE__,2,"PMT Spike module not found.");
		bPermitRun = false;
	}

  	eleMod = (TagElectrons*) ((TStnAna*) GetAna()->GetModule("TagElectrons"));
	if (!eleMod) {
		StdOut(__FILE__,__LINE__,2,"Electron tagging module not found.");
		bPermitRun = false;
	}

  	haloMod = (TagBeamHalo*) ((TStnAna*) GetAna()->GetModule("TagBeamHalo"));
	if (!haloMod) {
		StdOut(__FILE__,__LINE__,2,"Halo tagging module not found.");
		bPermitRun = false;
	}

  	setEMTMod = (SetEmTimes*) ((TStnAna*) GetAna()->GetModule("SetEmTimes"));
	if (!setEMTMod) {
		StdOut(__FILE__,__LINE__,2,"SetEMtimes module not found.");
		bPermitRun = false;
	}

  	phoenixMod = (TagPhoenixPhotons*) ((TStnAna*) GetAna()->GetModule("TagPhoenixPhotons"));
	if (!phoenixMod) {
		StdOut(__FILE__,__LINE__,2,"PhoenixPhotons Tagging module not found.");
		bPermitRun = false;
	}

	
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	//fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");
	//RegisterDataBlock("PROD@JetCluModule-had-cone0.4","TStnJetBlock",&fHadJetBlockClu04);
//	RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlockClu04);

	BookHistograms();

	// initialize vars
  

	return 0;
}

//_____________________________________________________________________________
int PhoJets::BeginRun()
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
int PhoJets::Event(int ientry)
{
  	SetPassed(0);

	// this makes sure that all required mods exists
	if (!bPermitRun) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
	}
	

	Cleanup();		//reinit the tree branch values
	hCount.iCounter->Fill(EVENTS_PROCESSED);
	FillDataBlocks(ientry);

	//event count 
	int iNtightPho = 0;
	int iNloosePho = 0;
	for (int i=0; i < initSpMod->GetSuperPhoSize(); i++)
	{
		if (initSpMod->GetSuperPhoton(i)->GetEtCorr() < GetPhoMinEt()) continue;
		if (initSpMod->GetSuperPhoton(i)->IsLoosePhoton())
		{
			if (initSpMod->GetSuperPhoton(i)->IsSidebandPhoton()) ++iNloosePho;
			if (initSpMod->GetSuperPhoton(i)->IsTightPhoton()) ++iNtightPho;
		}
	}

	if (iNloosePho>0 || iNtightPho>0) 
	{
		if (iNtightPho>0) hCount.iCounter->Fill(TIGHT_PHO_EVENTS);
		if (iNtightPho==2 && iNloosePho == 0) hCount.iCounter->Fill(TWO_TIGHT_PHO_ONLY);
		if (iNtightPho==1 && iNloosePho == 0) hCount.iCounter->Fill(ONE_TIGHT_PHO_ONLY);
		if (iNtightPho==1 && iNloosePho > 0) hCount.iCounter->Fill(ONE_TIGHT_ANY_LOOSE_PHO);
	}

	
	thePhoton = 0;			//reset the Photon
	std::vector<int> vPhoIndex, vEleIndex;
	bool bSkipThisEvent = false;

  	int NsuperPho = initSpMod->GetSuperPhoSize();
	for (int i=0; i < NsuperPho; i++) {
		
		//Det Eta cut
		if (fabs(initSpMod->GetSuperPhoton(i)->GetDetEta()) > GetPhoMaxDetEta()
			|| fabs(initSpMod->GetSuperPhoton(i)->GetDetEta()) < GetPhoMinDetEta())
		{
			continue;
		}

		 // reject phoenix photons
		if (GetRejPhxPhos() && initSpMod->GetSuperPhoton(i)->IsPhoenix()) continue;
		

		//remove halo for a given type. must be same type as used for 
		// halo template!!
		if (! fHeaderBlock->McFlag())
		{
			int iHaloId = initSpMod->GetSuperPhoton(i)->GetBeamHaloId();
			if ( ((iHaloId >> GetHaloType()) & 0x1) ) {
				bSkipThisEvent = true;
				break;
			}
		}

		// I should check Et>30 after BH cut. beam halo photons below 30 GeV
		// can be the extra photons. Then I need to throw that out. --12-11-2208
		if (initSpMod->GetSuperPhoton(i)->GetEtCorr() < GetPhoMinEt()) continue;

		bool bIsPho = 0;
		if (GetPhotonType() == iSignalPhoton)
		{
			bIsPho = initSpMod->GetSuperPhoton(i)->IsTightPhoton();
		} else if (GetPhotonType() == iSidebandPhoton)
		{
			bIsPho =  (! initSpMod->GetSuperPhoton(i)->IsTightPhoton()) & initSpMod->GetSuperPhoton(i)->IsLoosePhoton();
		} else {
			std::cout << __FILE__ << "::" << __LINE__ << "::Invalid Photon type requested!" << std::endl;
			exit (1);
		}

		//now if this is a good photon upto now, check if it is cosmic
		//photon
		if (bIsPho)
		{
				// to veto cosmics: if MC do nothing. if data: use muon stubs if run before EM timing system
				// else using EM timing cut
				if (GetHeaderBlock()->McFlag() == false)  
				{
					if (GetHeaderBlock()->RunNumber() <190851)
					{
						int iNCosStubs = initSpMod->GetSuperPhoton(i)->GetPhoton()->NCosStubPho();  // look for any muon stubs with the 30 degree cone of the photon
						if (iNCosStubs == 0) {
							vPhoIndex.push_back(i);
						} else 
						{
							bSkipThisEvent = true;
							break;
						}
					
					} else
					{
						float fEmTime = initSpMod->GetSuperPhoton(i)->GetEmTime();
						if (fEmTime > fEmTimeMin && fEmTime < fEmTimeMax)
						{
							vPhoIndex.push_back(i);
							//break; // NO! keep going so I can check if there is a BH photon. if so throw away the event!
						}
					}
			
				} else
				{
					vPhoIndex.push_back(i);
					//break; //after finding highest et tight photon, must not let override!
					// NO! keep going so I can check if there is a BH photon. if so throw away the event!
				}
				
			}
		
	}


	if (!bSkipThisEvent) {
		//exclude di-photon and photon-electron events only for the MEt plot
		//if (vPhoIndex.size() > 0 && vEleIndex.size() == 0) {
		if (vPhoIndex.size() > 0) {
			//make sure the photons are Et sorted
			if (vPhoIndex.size() > 1)
			{
				assert (initSpMod->GetSuperPhoton(vPhoIndex.at(0))->GetEtCorr() >
						initSpMod->GetSuperPhoton(vPhoIndex.at(1))->GetEtCorr()
						&& "PhoJets::Photons are not sorted!");
			}
			thePhoton = initSpMod->GetSuperPhoton(vPhoIndex.at(0));
			FillPhotonIDHist(vPhoIndex.at(0));		//Fill all the ID information of my photon for debugging
																// this way I can know if my selection is correct.
			hCount.iCounter->Fill(TIGHT_PHO_EVENTS);
			const bool bPassJetSelection = DoJetSelection(vPhoIndex.at(0));
			if (bPassJetSelection) SetPassed(1);
		}
	}

	if (GetPassed()) {
		hCount.iCounter->Fill(EVENTS_PASSED);
	}
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// Returns a pointer to the photon selected for this event. The photon
// that passed all selections criteria. Client should check for pointer
// status.
/*-------------------------------------------------------------------*/
SuperPhoton* PhoJets::GetPhoton() const
{
	if (thePhoton != NULL)
	{
		return thePhoton;
	} else
	{
		StdOut(__FILE__,__LINE__,2,"No qualified photon in this event! returning zero.");
		return 0;
	}
}

/*-------------------------------------------------------------------*/
void PhoJets::FillPhotonIDHist(const int iSpInd)
{
	if (iSpInd < 0 || iSpInd > initSpMod->GetSuperPhoSize())
	{
		StdOut(__FILE__,__LINE__,3,"Invalid SuperPhoton index!");
		return;
	}
	
	int TightPhoId 	= initSpMod->GetSuperPhoton(iSpInd)->GetTightPhotonId();
	int LoosePhoId 	= initSpMod->GetSuperPhoton(iSpInd)->GetLoosePhotonId();
	int PhoLikeEleId 	= initSpMod->GetSuperPhoton(iSpInd)->GetTightElectronId();
	int StdLooseEleId = initSpMod->GetSuperPhoton(iSpInd)->GetStdLooseElectronId();
	int StdTightEleId = initSpMod->GetSuperPhoton(iSpInd)->GetStdTightElectronId();
	int iHaloId 			= initSpMod->GetSuperPhoton(iSpInd)->GetBeamHaloId();
	int PhoenixId 		= initSpMod->GetSuperPhoton(iSpInd)->GetPhoenixId();
	int ConversionId 	= initSpMod->GetSuperPhoton(iSpInd)->GetConversionId();

	if (TightPhoId == 0) h1_Pho.PhoIDs->Fill(0);
	if (LoosePhoId == 0) h1_Pho.PhoIDs->Fill(1);
	if (TightPhoId>0 && LoosePhoId == 0) h1_Pho.PhoIDs->Fill(2);
	if (PhoLikeEleId == 0)	h1_Pho.PhoIDs->Fill(3);
	if (StdTightEleId == 0) h1_Pho.PhoIDs->Fill(4);
	if (StdLooseEleId == 0) h1_Pho.PhoIDs->Fill(5);
	if (PhoenixId > 0) h1_Pho.PhoIDs->Fill(6);
	if (ConversionId > 0) h1_Pho.PhoIDs->Fill(7);
	if ( ((iHaloId >> GetHaloType()) & 0x1) )
	{
		if (iHaloId > 0) h1_Pho.PhoIDs->Fill(8);
	}
	
}


			
/*-------------------------------------------------------------------*/
bool PhoJets::DoJetSelection(int iSpIndex)
{
	
	//try met cut here

	bool bPassJetSelection = false;
	//find the leading detector jet after removing the leading photon
	TLorentzVector tlPhoVec = initSpMod->GetSuperPhoton(iSpIndex)->GetCorVec();
	TVector2 tv2MetVec(jetMod->GetMyMetXCorr(0,3), jetMod->GetMyMetYCorr(0,3));

	if (tv2MetVec.Mod() < GetMinMEt()) return 0;

	if (jetMod->GetMyNjet(0,3) >= 1)
	{	// mono-jet case
		//hPhoLeadJetPtRatio->Fill(tlLeadJetVec.Pt()/tlPhoVec.Pt() - 1);
		//std::cout << "========= "; fHeaderBlock->Print();
		h1_Pho.EtCorr->Fill(tlPhoVec.Pt());
		const float dphi_pm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlPhoVec.Phi())));
		h1_Pho.PhoMetDelPhi->Fill(dphi_pm);
		metb4->Fill(tv2MetVec.Mod());
			
		//std::cout << "pho pt, met, delphi = " << tlPhoVec.Pt() << ", " << tv2MetVec.Mod() << ", " <<  dphi_pm << std::endl;

		float dphi_jm_closest = 10.0;
		TLorentzVector tlVecClosestJet(0,0,0,0);
		for (int ijet=0; ijet<jetMod->GetMyNjet(0,0); ++ijet)
		{
			TLorentzVector tlLeadJetVec = *(jetMod->GetMyJetNoEMobj_lev6(0,ijet));
			if (tlLeadJetVec.Pt()>15.)
			{
				//hLeadJet_Et->Fill(tlLeadJetVec.Pt());
				//h1_Jet.EtCorr->Fill(tlLeadJetVec.Pt());
				//TLorentzVector sum = tlPhoVec + tlLeadJetVec;
				//h1_PhoJet.InvMass->Fill(sum.M());
				const float dphi_jm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlLeadJetVec.Phi())));
				if (dphi_jm < dphi_jm_closest) 
				{
					dphi_jm_closest = dphi_jm;
					tlVecClosestJet = tlLeadJetVec;
					//std::cout <<  ijet << "\t"<< tlLeadJetVec.Pt() << "\t" << dphi_jm << " -> " <<  dphi_jm_closest << std::endl;
				} else 
				{
					//std::cout <<  ijet << "\t"<< tlLeadJetVec.Pt() << "\t" << dphi_jm  << std::endl;
				}
			} else 
			{
				//std::cout <<  ijet << "\t"<< tlLeadJetVec.Pt() << "\t low pt, ignored!" << std::endl;
			}
		}

		if (dphi_jm_closest<10.0) 
		{
			h1_Jet.JetMetDelPhi->Fill(dphi_jm_closest);
			h1_Jet.EtCorr->Fill(tlVecClosestJet.Pt());
			hCount.iCounter->Fill(2);

		}

		//if (dphi_jm_closest>0.4 && dphi_pm>0.4)
		if (dphi_jm_closest>0.4) //not bad met event
		{
			hCount.iCounter->Fill(TIGHT_PHO1JET_EVENTS);
			bPassJetSelection = true;
			meta4->Fill(tv2MetVec.Mod());
			jetmetdelphi_a4->Fill(dphi_jm_closest);
			phometdelphi_a4->Fill(dphi_pm);

			evtPropMod->FillPhotonHists(h1_Pho,iSpIndex);
			evtPropMod->FillEventHists(h1_Evt);
			evtPropMod->FillJetHists(h1_Jet,0);		//lead jet
			evtPropMod->FillPhoton1JetHists(h1_PhoJet,iSpIndex,0);

			if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
				hCount.iCounter->Fill(TIGHT_PHO2JET_EVENTS);
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

		
	}

	return bPassJetSelection;
	
		//evtPropMod->FillPhotonHists(h1_Pho,iSpIndex);
		//evtPropMod->FillEventHists(h1_Evt);
/*	if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
		write = true;
		SetPassed(1);			//making exclusive templates
		hCount.iCounter->Fill(2);
		evtPropMod->FillPhotonHists(h1_Pho,iSpIndex);
		evtPropMod->FillEventHists(h1_Evt);
		evtPropMod->FillJetHists(h1_Jet,0);		//lead jet
		evtPropMod->FillPhoton1JetHists(h1_PhoJet,iSpIndex,0);
	}
	
	if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
		write = true;
		SetPassed(1);			//making exclusive templates
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
	
*/	
}


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void PhoJets::Cleanup()
{
} //Cleanup


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void PhoJets::SetHaloType(const int iT)
{ 
	if (iT <0 || iT >5) {
		StdOut(__FILE__,__LINE__,3,"No such halo type. pl. pick from 0-5.");
		exit (1);
	} else {
		iHaloType = iT;
	}
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void PhoJets::SetPhotonType(const int type)
{
	if (type  < 0 || type > 1)
	{
		StdOut(__FILE__,__LINE__, 3, "Invalid photon type requested! must be 0(signal) or 1 (sideband)!. exiting!");
		exit (1);
	} else
	{
		iPhotonType = type;
	}
}



/*-------------------------------------------------------------------*/
void PhoJets::FillDataBlocks(int ientry)
{
	//fGenpBlock->GetEntry(ientry);
//	fJetBlockClu04->GetEntry(ientry);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int PhoJets::EndJob() {

	if (GetSummaryStat()) return 0;
		

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[PJS:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[PJS:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(EVENTS_PROCESSED+1)) + "\n";
	sMsg += "[PJS:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(EVENTS_PASSED+1)) + "\n";
	sMsg += "[PJS:03:] Photon selection ----------- = " + ToStr(GetPhotonType()) + " (" + GetPhotonTypeLabel() + ")"+ "\n";
	sMsg += "[PJS:04:] Pho/Jet Min Et ------------- = " + ToStr(GetPhoMinEt()) + " / " + ToStr(GetJetMinEt()) + " GeV" + "\n";
	sMsg += "[PJS:05:] Pho/Jet |DetEta|<Max ------- = " + ToStr(GetPhoMaxDetEta()) + " / " + ToStr(GetJetMaxDetEta()) +  "\n";
	if (GetRejPhxPhos())
	sMsg += "[PJS:06:] Phoenix Photons ------------ = " + ToStr(GetRejPhxPhos()) + " (Rejected)" + "\n";
	else 
	sMsg += "[PJS:06:] Phoenix Photons ------------ = " + ToStr(GetRejPhxPhos()) + " (NOT Rejected)" + "\n";
	sMsg += "[PJS:07:] Halo Rejection Type -------- = " + ToStr(GetHaloType()) + "\n";
	if (GetHeaderBlock()->McFlag())
	sMsg += "[PJS:08:] EM Time window used -------- = NO CUT (MC)\n";
	else
	sMsg += "[PJS:08:] EM Time window used -------- = " + ToStr(GetEmTimeWindowMin()) + "," + ToStr(GetEmTimeWindowMax())+ "\n";
	sMsg += "[PJS:09:] Incl 1 Tight Photon Events - = " + ToStr(hCount.iCounter->GetBinContent(TIGHT_PHO_EVENTS+1)) + "\n";
	sMsg += "[PJS:10:] Incl 1g >=1jet event ------- = " + ToStr(hCount.iCounter->GetBinContent(TIGHT_PHO1JET_EVENTS+1)) + "\n";
	sMsg += "[PJS:11:] Incl 1g >=2jet event ------- = " + ToStr(hCount.iCounter->GetBinContent(TIGHT_PHO2JET_EVENTS+1)) + "\n";
	sMsg += "[PJS:12:] 1 Tight Photon only Events - = " + ToStr(hCount.iCounter->GetBinContent(ONE_TIGHT_PHO_ONLY+1)) + "\n";
	sMsg += "[PJS:13:] 2 Tight Photons Only Events  = " + ToStr(hCount.iCounter->GetBinContent(TWO_TIGHT_PHO_ONLY+1)) + "\n";
	sMsg += "[PJS:14:] 1 T+>0 Loose Pho Only Events = " + ToStr(hCount.iCounter->GetBinContent(ONE_TIGHT_ANY_LOOSE_PHO+1)) + "\n";
	sMsg += "---------------------------------------------------\n";

	std::cout << sMsg;

	return 0;
}



void PhoJets::SetPhoMaxDetEta(const float eta)
{
	std::cout << "eta = " << fabs(eta) << std::endl;
	if (fabs(eta) > 1.101)		//1.101 is a hack. for some reason when I require
	{									// eta=1.1, and compare (>) to 1.1, this condition becomes true!
										// sam -11-12-2009
		StdOut(__FILE__,__LINE__,3,"Photon Detector eta you required is > than what is returned in TStnPhoton::Detector() == 0??");
		exit (1);
	} else {
 		fPhoMaxDetEta = fabs(eta);
	}
}

void PhoJets::SetPhoMinDetEta(const float eta)
{
	if (fabs(eta) > 1.1)
	{
		StdOut(__FILE__,__LINE__,3,"Photon Detector eta require is > than what is returned in TStnPhoton::Detector() == 0??");
		exit (1);
	} else {
 		fPhoMinDetEta = fabs(eta);
	}
}

void PhoJets::SetPhoMinEt(const float et)
{ 
	if (et < 0)
	{
		StdOut(__FILE__,__LINE__,3,"Photon Et must be > 0. Using default value.");
		exit(-1);
	} else
	{
		fPhoMinEt = et;
	}
}
