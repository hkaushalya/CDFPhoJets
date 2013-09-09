////////////////////////////////////////////////////////////
// This is to do the event count as I apply the cuts. I   //
// this is get the final event counts in my dissertation. //
////////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/*{{{*/
/* $Id: PhoJetsCount.cc,v 1.1 2011/05/25 21:26:32 samantha Exp $
 * $Log: PhoJetsCount.cc,v $
 * Revision 1.1  2011/05/25 21:26:32  samantha
 * This is to do the event count as I apply the cuts. I
 * this is get the final event counts in my dissertation.
 *
 */ 
/*}}}*/

#include "samantha/Pho/PhoJetsCount.hh"
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


ClassImp(PhoJetsCount)


//_____________________________________________________________________________
PhoJetsCount::PhoJetsCount(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bPermitRun(true),
  fEmTimeMin(-4.8),
  fEmTimeMax(4.8),
  iPhotonType(0),
  fPhoMaxDetEta(1.1),
  fPhoMinDetEta(0),
  fPhoMinEt(30.),
  fJetMinEt(15.),
  fJetMaxDetEta(3.0),
  bRejPhxPhos(true),
  iHaloType(5),
  bNoSummary(false)
{
	std::cout << "Hello I am PhoJetsCount module" << std::endl;
}

//_____________________________________________________________________________
PhoJetsCount::~PhoJetsCount() {
}

//_____________________________________________________________________________
void PhoJetsCount::SaveHistograms() {
}

//_____________________________________________________________________________
void PhoJetsCount::BookHistograms()
{
	DeleteHistograms();

	TFolder *new_folder, *sub_folder, *new_folder2;;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");				// 0
	vCounterHistLabels.push_back("Pass Goodrun");
	vCounterHistLabels.push_back("Pass Trigger");
	vCounterHistLabels.push_back("Pass Nvtx");
	vCounterHistLabels.push_back("Pass Vtxz");
	vCounterHistLabels.push_back("Pass Pho Et");
	vCounterHistLabels.push_back("Pass Tight Pho Id");
	vCounterHistLabels.push_back("Pass Not PMT");
	vCounterHistLabels.push_back("Pass EM Time");
	vCounterHistLabels.push_back("Pass No Muon Stubs");
	vCounterHistLabels.push_back("Pass Not Phoenix");
	vCounterHistLabels.push_back("Pass Not BH5");
	vCounterHistLabels.push_back("Pass 1Njet15");
	//vCounterHistLabels.push_back("Pass Jet Et");
	//vCounterHistLabels.push_back("Pass Jet Eta");
	vCounterHistLabels.push_back("Pass Met 20");
	vCounterHistLabels.push_back("Pass Met 50");
	vCounterHistLabels.push_back("Events Passed");					// 1
	vCounterHistLabels.push_back("All Tight #gamma Events");		// 2
	vCounterHistLabels.push_back("One tight #gamma only Events");	// 3
	vCounterHistLabels.push_back("Two tight #gamma only events");	// 4
	vCounterHistLabels.push_back("One tight #gamma and any Nloose pho events");	// 5
	HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	

	new_folder = GetHistFolder(this, "1Jet","#gamma + >=1 Jet case");
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt,"#gamma + >=1 Jets(15GeV)");

	//metb4 = dynamic_cast<TH1F*> (h1_Evt.Met->Clone("Met_b4cut"));
	//meta4 = dynamic_cast<TH1F*> (h1_Evt.Met->Clone("Met_a4cut"));
	//sub_folder->Add(metb4);
	//sub_folder->Add(meta4);
	
	
	//before cuts

	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h1_Pho,"#gamma + >=1 Jets(15GeV)");
	
	//phometdelphi_a4 = dynamic_cast<TH1F*> (h1_Pho.PhoMetDelPhi->Clone("PhoMetDelPhi_a4"));
	//sub_folder->Add(phometdelphi_a4);

	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet,"#gamma + >=1 Jets(15GeV)");
	

	//jetmetdelphi_a4 = dynamic_cast<TH1F*> (h1_Jet.JetMetDelPhi->Clone("JetMetDelPhi_a4"));
	//sub_folder->Add(jetmetdelphi_a4);
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet,"#gamma + >=1 Jets(15GeV)");

	//after n-1 cuts
	
	//new_folder2 = GetHistFolder(this, "AfterCuts","After N-1 cuts");
	new_folder = GetHistFolder(this, "1JetAfter","#gamma + >=1 Jet case after cuts");

	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt_a4,"#gamma + >=1 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h1_Pho_a4,"#gamma + >=1 Jets(15GeV)");
	
	//phometdelphi_a4 = dynamic_cast<TH1F*> (h1_Pho.PhoMetDelPhi->Clone("PhoMetDelPhi_a4"));
	//sub_folder->Add(phometdelphi_a4);

	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet_a4,"#gamma + >=1 Jets(15GeV)");
	

	//jetmetdelphi_a4 = dynamic_cast<TH1F*> (h1_Jet.JetMetDelPhi->Clone("JetMetDelPhi_a4"));
	//sub_folder->Add(jetmetdelphi_a4);
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet_a4,"#gamma + >=1 Jets(15GeV)");


	

/*
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
	HistoMan.GetPhoton2JetsHistograms(sub_folder,h2_PhoJetsCount,"#gamma + >=2 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
	HistoMan.GetTwoJetsHistograms(sub_folder,h2_TwoJets,"#gamma + >=2 Jets(15GeV)");
	*/

}


//_____________________________________________________________________________
int PhoJetsCount::BeginJob()
{

	//check if the photon Eta range is correct
	if (fPhoMaxDetEta < fPhoMinDetEta)
	{
		StdOut(__FILE__,__LINE__,3,"min photon DetEta larger than max DetEta. Correct this!.");
		assert (false);
	}

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

  	//evtPropMod = (EventProperties*) ((TStnAna*) GetAna()->GetModule("EventProperties"));
	//if (!evtPropMod) {
	//	StdOut(__FILE__,__LINE__,2,"EventProperties tag module not found.");
	//	bPermitRun = false;
	//}

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
  	trigMod = (TriggerModule*) ((TStnAna*) GetAna()->GetModule("TriggerModule"));
	if (!trigMod) {
		StdOut(__FILE__,__LINE__,2,"Trigger module not found.");
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
int PhoJetsCount::BeginRun()
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
int PhoJetsCount::Event(int ientry)
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
	
	if (trigMod->PassGoodRun())
	{ /*std::cout << "pass good run" << std::endl;*/  
		hCount.iCounter->Fill(PASS_GOORUN);
	} else return 0;

	if (trigMod->PassAnyPhoTriggers()) hCount.iCounter->Fill(PASS_TRIGGER); else return 0;
	if (trigMod->GetN12vertex()>0) hCount.iCounter->Fill(PASS_NVTX); else return 0;
	if (trigMod->GetBestVertexZ()<60.) hCount.iCounter->Fill(PASS_VTXZ); else return 0;

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
	std::vector<int> vPhoIndex;
	bool bSkipThisEvent = false;
	bool bPhoWithEt30 = false;

  	int NsuperPho = initSpMod->GetSuperPhoSize();
	//for each photon candidate make n-1 plots
	for (int i=0; i < NsuperPho; i++) 
	{
		SuperPhoton* sp = initSpMod->GetSuperPhoton(i);
		//if (sp->GetTightPhotonId() != 0) continue;
		//if (fabs(sp->GetDetEta()) >1.1) continue;
		//if (sp->GetEtCorr() < GetPhoMinEt()) continue;

		//try to find the leading photon for analysis
		//const float EtCorr = sp->GetEtCorr();
		//bool pass = false;
		//if (trigMod->GetTrig70Bit() == 1 && EtCorr < 70.0) pass = true;
		//else if (trigMod->GetTrig50Bit() == 1 && EtCorr < 50.0) pass = true;
		//else if (trigMod->GetTrig25IsoBit() == 1 && EtCorr < 25.0) pass = true;

		FillPhotonIDHist(h1_Pho, i);
		FillPhotonHists(h1_Pho, sp); //fill photon ID variables before any cut
		//n-1 plots
		//detector
		//std::cout << "====================== looping for photon " << i << std::endl;
		if (PassEtTrigMatched(sp) && PassXCes(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassN3d(sp) && PassTrkPt(sp) && PassTrkIso(sp)
			&& PassCes(sp)) 
		{
			std::cout << "pho["<< i << "] Eta = " << sp->GetDetEta() << std::endl;
			h1_Pho_a4.DetEta->Fill(sp->GetDetEta());
		}
			 //h1_Pho_a4.Detector->Print();

		if (PassDetEta(sp) && PassXCes(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassN3d(sp) && PassTrkPt(sp) && PassTrkIso(sp)
			&& PassCes(sp)) h1_Pho_a4.EtCorr->Fill(sp->GetEtCorr());
			//h1_Pho_a4.EtCorr->Print();
		
		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassN3d(sp) && PassTrkPt(sp) && PassTrkIso(sp)
			&& PassCes(sp)) h1_Pho_a4.XCes->Fill(sp->GetXCes());
		
		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassXCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassN3d(sp) && PassTrkPt(sp) && PassTrkIso(sp)
			&& PassCes(sp)) h1_Pho_a4.ZCes->Fill(sp->GetZCes());

		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassXCes(sp) && PassZCes(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassN3d(sp) && PassTrkPt(sp) && PassTrkIso(sp)
			&& PassCes(sp)) h1_Pho_a4.HadEm->Fill(sp->GetHadEm());
	
		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassXCes(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassChi2Mean(sp) && PassN3d(sp) && PassTrkPt(sp) && PassTrkIso(sp)
			&& PassCes(sp)) h1_Pho_a4.Iso4->Fill(sp->GetIso4());
	
		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassXCes(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassN3d(sp) && PassTrkPt(sp) && PassTrkIso(sp)
			&& PassCes(sp)) h1_Pho_a4.Chi2Mean->Fill(sp->GetChi2Mean());
	
		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassXCes(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassTrkPt(sp) && PassTrkIso(sp)
			&& PassCes(sp)) h1_Pho_a4.N3d->Fill(sp->GetN3d());
	
		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassXCes(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassN3d(sp) && PassTrkIso(sp)
			&& PassCes(sp)) h1_Pho_a4.TrkPt->Fill(sp->GetTrkPt());
	
		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassXCes(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassN3d(sp) && PassTrkPt(sp) 
			&& PassCes(sp)) h1_Pho_a4.TrkIso->Fill(sp->GetTrkIso());
	
		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassXCes(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassN3d(sp) && PassTrkPt(sp) && PassTrkIso(sp)) 
		{
			h1_Pho_a4.Ces2Wire->Fill(sp->GetCesWireE2());
			h1_Pho_a4.Ces2Strip->Fill(sp->GetCesStripE2());
		}

		//FILL ONCE IF A PHOTON WITH min et is found
		if (! bPhoWithEt30 && sp->GetEtCorr()> GetPhoMinEt()) 
		{ 
			hCount.iCounter->Fill(PASS_ET_TRIGMATCHED);
			bPhoWithEt30 = true;
		}
		//cross check to see if I have implemented the cut scorrectly here.
		if (PassDetEta(sp) && PassEtTrigMatched(sp) && PassXCes(sp) && PassZCes(sp) && PassHadEm(sp)
			&& PassIso(sp) && PassChi2Mean(sp) && PassN3d(sp) && PassTrkPt(sp) && PassTrkIso(sp)
		//try to find the leading photon for analysis
			&& PassCes(sp)) 
		{
			vPhoIndex.push_back(i); 
			assert(sp->GetTightPhotonId() == 0 && "Tight pho found but does not match the tag!");
		}
	}

	if (vPhoIndex.size() >0) hCount.iCounter->Fill(PASS_TIGHTPHOID); else return 0;  
	std::cout << "Tight pho found" << std::endl;
	SuperPhoton* sp = initSpMod->GetSuperPhoton(vPhoIndex.at(0));
	thePhoton = initSpMod->GetSuperPhoton(vPhoIndex.at(0));

	//pmt
	if (sp->GetPMTSpikeId() == 0) hCount.iCounter->Fill(PASS_NOT_PMT); else return 0;

	//cosmic
	if (fabs(sp->GetEmTime())<4.8) hCount.iCounter->Fill(PASS_EMTIME); else return 0;
	if (sp->GetNMuonStubs()<1) hCount.iCounter->Fill(PASS_NO_MUON_STUBS); else return 0;
	
	//phoenix track //leptons
	if (sp->GetPhoenixId() == 0) hCount.iCounter->Fill(PASS_NOT_PHOENIX); else return 0;
						
	
	//bh
	if (sp->GetBeamHaloId() == 0) hCount.iCounter->Fill(PASS_NOT_BH5); else return 0;
						
	
	///jet selections
	///first look for extra em objects that is removed from the jet list.
	int hasExtraEm = jetMod->GetMyNCorrPho() + jetMod->GetMyNCorrEle();
	if (vPhoIndex.size()>1) hasExtraEm = 1;
	for (int i=0; i < NsuperPho; i++) 
	{
		SuperPhoton* sp = initSpMod->GetSuperPhoton(i);
		if (sp->GetEtCorr() < GetJetMinEt()) continue;
		if (sp->GetStdLooseElectronId() == 0 || fabs(sp->GetDetEta()) >1.1 || sp->GetEtCorr() < 30.) continue;
	}

	DoJetSelection(vPhoIndex.at(0));
	

	//delphi(jet,met)
	

		
	SetPassed(1);

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
SuperPhoton* PhoJetsCount::GetPhoton() const
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
void PhoJetsCount::FillPhotonIDHist(Histograms::PhotonHists_t& hist, const int iSpInd)
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

	if (TightPhoId == 0) hist.PhoIDs->Fill(0);
	if (LoosePhoId == 0) hist.PhoIDs->Fill(1);
	if (TightPhoId>0 && LoosePhoId == 0) hist.PhoIDs->Fill(2);
	if (PhoLikeEleId == 0)	hist.PhoIDs->Fill(3);
	if (StdTightEleId == 0) hist.PhoIDs->Fill(4);
	if (StdLooseEleId == 0) hist.PhoIDs->Fill(5);
	if (PhoenixId > 0) hist.PhoIDs->Fill(6);
	if (ConversionId > 0) hist.PhoIDs->Fill(7);
	if ( ((iHaloId >> GetHaloType()) & 0x1) )
	{
		if (iHaloId > 0) hist.PhoIDs->Fill(8);
	}
	
}


			
/*-------------------------------------------------------------------*/
void PhoJetsCount::DoJetSelection(int iSpIndex)
{
	
	//find the leading detector jet after removing the leading photon
	TLorentzVector tlPhoVec = initSpMod->GetSuperPhoton(iSpIndex)->GetCorVec();
	TVector2 tv2MetVec(jetMod->GetMyMetXCorr(0,3), jetMod->GetMyMetYCorr(0,3));
	
	if (jetMod->GetMyNjet(0,3) >= 1)
	{	
		hCount.iCounter->Fill(PASS_1NJET15);
		// mono-jet case
		//hPhoLeadJetPtRatio->Fill(tlLeadJetVec.Pt()/tlPhoVec.Pt() - 1);
		//std::cout << "========= "; fHeaderBlock->Print();
		//h1_Pho.EtCorr->Fill(tlPhoVec.Pt());
		const float dphi_pm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlPhoVec.Phi())));
		//h1_Pho.PhoMetDelPhi->Fill(dphi_pm);
		//metb4->Fill(tv2MetVec.Mod());
			
		//std::cout << "pho pt, met, delphi = " << tlPhoVec.Pt() << ", " << tv2MetVec.Mod() << ", " <<  dphi_pm << std::endl;

		float dphi_jm_closest = 10.0;
		TLorentzVector tlVecClosestJet(0,0,0,0);
		for (int ijet=0; ijet<jetMod->GetMyNjet(0,0); ++ijet)
		{
			TLorentzVector tlLeadJetVec = *(jetMod->GetMyJetNoEMobj_lev6(0,ijet));
			if (tlLeadJetVec.Pt()>GetJetMinEt() && GetJetMaxDetEta())
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
			//h1_Jet.JetMetDelPhi->Fill(dphi_jm_closest);
			//h1_Jet.EtCorr->Fill(tlVecClosestJet.Pt());
			//hCount.iCounter->Fill(2);

		}

		//if (dphi_jm_closest>0.4 && dphi_pm>0.4)
		if (dphi_jm_closest>0.4)
		{
			
			hCount.iCounter->Fill(PASS_DELPHI_JETMET);
			//meta4->Fill(tv2MetVec.Mod());
			//jetmetdelphi_a4->Fill(dphi_jm_closest);
			//phometdelphi_a4->Fill(dphi_pm);
		}

		
	}

	if (tv2MetVec.Mod()>20.0) hCount.iCounter->Fill(PASS_MET_20);
	if (tv2MetVec.Mod()>50.0) hCount.iCounter->Fill(PASS_MET_50);
	
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
		evtPropMod->FillPhoton2JetsHists(h2_PhoJetsCount,iSpIndex,0,1);
		evtPropMod->FillTwoJetsHists(h2_TwoJets,0,1);
	}
	
*/	
}


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void PhoJetsCount::Cleanup()
{
} //Cleanup


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void PhoJetsCount::SetHaloType(const int iT)
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
void PhoJetsCount::SetPhotonType(const int type)
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
void PhoJetsCount::FillDataBlocks(int ientry)
{
	//fGenpBlock->GetEntry(ientry);
//	fJetBlockClu04->GetEntry(ientry);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int PhoJetsCount::EndJob() {

	if (GetSummaryStat()) return 0;
		

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[P2J:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[P2J:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[P2J:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	sMsg += "[P2J:03:] Photon selection ----------- = " + ToStr(GetPhotonType()) + " (" + GetPhotonTypeLabel() + ")"+ "\n";
	sMsg += "[P2J:04:] Pho/Jet Min Et ------------- = " + ToStr(GetPhoMinEt()) + " / " + ToStr(GetJetMinEt()) + " GeV" + "\n";
	sMsg += "[P2J:05:] Pho/Jet |DetEta|<Max ------- = " + ToStr(GetPhoMaxDetEta()) + " / " + ToStr(GetJetMaxDetEta()) +  "\n";
	if (GetRejPhxPhos())
	sMsg += "[P2J:06:] Phoenix Photons ------------ = " + ToStr(GetRejPhxPhos()) + " (Rejected)" + "\n";
	else 
	sMsg += "[P2J:06:] Phoenix Photons ------------ = " + ToStr(GetRejPhxPhos()) + " (NOT Rejected)" + "\n";
	sMsg += "[P2J:07:] Halo Rejection Type -------- = " + ToStr(GetHaloType()) + "\n";
	if (GetHeaderBlock()->McFlag())
	sMsg += "[P2J:07:] EM Time window used -------- = NO CUT (MC)\n";
	else
	sMsg += "[P2J:08:] EM Time window used -------- = " + ToStr(GetEmTimeWindowMin()) + "," + ToStr(GetEmTimeWindowMax())+ "\n";
	sMsg += "[P2J:09:] Incl 1 Tight Photon Events - = " + ToStr(hCount.iCounter->GetBinContent(TIGHT_PHO_EVENTS+1)) + "\n";
	sMsg += "[P2J:10:] 1 Tight Photon only Events - = " + ToStr(hCount.iCounter->GetBinContent(ONE_TIGHT_PHO_ONLY+1)) + "\n";
	sMsg += "[P2J:11:] 2 Tight Photons Only Events  = " + ToStr(hCount.iCounter->GetBinContent(TWO_TIGHT_PHO_ONLY+1)) + "\n";
	sMsg += "[P2J:12:] 1 T+>0 Loose Pho Only Events = " + ToStr(hCount.iCounter->GetBinContent(ONE_TIGHT_ANY_LOOSE_PHO+1)) + "\n";
	sMsg += "---------------------------------------------------\n";

	std::cout << sMsg;
	std::cout << "Processed          = " << hCount.iCounter->GetBinContent(EVENTS_PROCESSED+1) << std::endl;
	std::cout << "Goodrun            = " << hCount.iCounter->GetBinContent(PASS_GOORUN+1) << std::endl;
	std::cout << "Trigger            = " << hCount.iCounter->GetBinContent(PASS_TRIGGER+1) << std::endl;
	std::cout << "Nvtx               = " << hCount.iCounter->GetBinContent(PASS_NVTX+1) << std::endl;
	std::cout << "VtxZ               = " << hCount.iCounter->GetBinContent(PASS_VTXZ+1) << std::endl;
	std::cout << "Et (trig mathced)  = " << hCount.iCounter->GetBinContent(PASS_ET_TRIGMATCHED+1) << std::endl;
	std::cout << "tight pho        	= " << hCount.iCounter->GetBinContent(PASS_TIGHTPHOID+1) << std::endl;
	std::cout << "NOT PMT            = " << hCount.iCounter->GetBinContent(PASS_NOT_PMT+1) << std::endl;
	std::cout << "Pass EM TIME   		= " << hCount.iCounter->GetBinContent(PASS_EMTIME+1) << std::endl;
	std::cout << "NO MUON STUBS      = " << hCount.iCounter->GetBinContent(PASS_NO_MUON_STUBS+1) << std::endl;
	std::cout << "NO TRACK           = " << hCount.iCounter->GetBinContent(PASS_NOT_PHOENIX+1) << std::endl;
	std::cout << "NOT HALO           = " << hCount.iCounter->GetBinContent(PASS_NOT_BH5+1) << std::endl;
	std::cout << "njet15>=1          = " << hCount.iCounter->GetBinContent(PASS_1NJET15+1) << std::endl;
	//std::cout << "NO JET ET          = " << hCount.iCounter->GetBinContent(PASS_JET_ET+1) << std::endl;
	//std::cout << "NO JET ETA         = " << hCount.iCounter->GetBinContent(PASS_JET_ETA+1) << std::endl;
	std::cout << "PASS DelPhi(JetMet)= " << hCount.iCounter->GetBinContent(PASS_DELPHI_JETMET+1) << std::endl;
	std::cout << "PASS MET 20          = " << hCount.iCounter->GetBinContent(PASS_MET_20+1) << std::endl;
	std::cout << "PASS MET 50          = " << hCount.iCounter->GetBinContent(PASS_MET_50+1) << std::endl;

	return 0;
}



void PhoJetsCount::SetPhoMaxDetEta(const double eta)
{
	//std::cout << "eta = " << fabs(eta) << std::endl;
	fPhoMaxDetEta = fabs(eta);
}

void PhoJetsCount::SetPhoMinDetEta(const double eta)
{
	fPhoMinDetEta = fabs(eta);
}

void PhoJetsCount::SetPhoMinEt(const float et)
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


bool PhoJetsCount::PassDetEta(const SuperPhoton *sp)
{
	bool pass = false;
	if (fabs(sp->GetDetEta())< fabs(GetPhoMaxDetEta())) pass = true;
	//if (!pass) std::cout << "\tdetector failed "<< std::endl;
	return pass;
}

bool PhoJetsCount::PassEtTrigMatched(const SuperPhoton *sp)
{
	bool pass = false;
	const float EtCorr = sp->GetEtCorr();
	
	if (trigMod->GetTrig70Bit() == 1 && EtCorr > 70.0) pass = true;
	else if (trigMod->GetTrig50Bit() == 1 && EtCorr > 50.0) pass = true;
	else if (trigMod->GetTrig25IsoBit() == 1 && EtCorr > 25.0) pass = true;
	
	if (EtCorr<GetPhoMinEt()) pass = false;

	//if (!pass) std::cout << "\tEtc failed  "<< std::endl;
	return pass;
}
bool PhoJetsCount::PassXCes(const SuperPhoton* sp)
{
	bool pass = false;
	if (fabs(sp->GetXCes()) < 21.) pass = true;
	//if (!pass) std::cout << "\t xces failed!  xces =" << sp->GetXCes() << std::endl;
	return pass;
}
bool PhoJetsCount::PassZCes(const SuperPhoton* sp)
{
	bool pass = false;
	if (fabs(sp->GetZCes()) > 9 &&  fabs(sp->GetZCes()) < 230)  pass = true;
	//if (!pass) std::cout << "\t zces failed! zces =" << sp->GetZCes() << std::endl;
	return pass;
}
bool PhoJetsCount::PassHadEm(const SuperPhoton* sp)
{
	bool pass = false;
	const float HadEm = sp->GetHadEm();
	const float ECorr = sp->GetECorr();
	
	if ( (HadEm < 0.125) || (HadEm < (0.055 + 0.00045 * ECorr)) ) pass = true;
	//if (!pass) std::cout << "\thadem failed  "<< std::endl;
	return pass;
}

bool PhoJetsCount::PassIso(const SuperPhoton *sp)
{
	bool pass = false;
	const float EtCorr = sp->GetEtCorr();
	const float IsoEtCorr = sp->GetIso4();
	
	if (EtCorr < 20) 
	{
		if (IsoEtCorr < 0.1*EtCorr) pass = true;
		
	} else {
		if (IsoEtCorr < (2.0+0.02*(EtCorr-20.0)) ) pass = true;
	}
	//if (!pass) std::cout << "\tiso failed  "<< std::endl;
	return pass;
}

bool PhoJetsCount::PassChi2Mean(const SuperPhoton* sp)
{
	bool pass = false;
	if (sp->GetChi2Mean() < 20.) pass = true;
	//if (!pass) std::cout << "\tchi2 failed  "<< std::endl;
	return pass;
}

bool PhoJetsCount::PassN3d(const SuperPhoton* sp)
{
	bool pass = false;
	const int N3d = sp->GetN3d();
	if (N3d>=0 && N3d <= 1) pass = true;
	//if (!pass) std::cout << "\tn3d failed  "<< std::endl;
	return pass;
}

bool PhoJetsCount::PassTrkPt(const SuperPhoton* sp)
{
	bool pass = false;
	const int N3d = sp->GetN3d();
	const float TrkPt = sp->GetTrkPt();
	const float EtCorr = sp->GetEtCorr();
	if (N3d >= 1) //if there is a track
	{ 
		if (TrkPt < (1+0.005*EtCorr))  pass = true;
	} else pass = true;

	//if (!pass) std::cout << "\t trk pt failed  "<< std::endl;
	return pass;
}

bool PhoJetsCount::PassTrkIso(const SuperPhoton* sp)
{
	bool pass = false;
	const float TrkIso = sp->GetTrkIso();
	const float EtCorr = sp->GetEtCorr();
	if (TrkIso < (2.0+0.005*EtCorr)) pass = true;
	//if (!pass) std::cout << "\tTrkiso failed  "<< std::endl;
	return pass;
}

bool PhoJetsCount::PassCes(const SuperPhoton* sp)
{
	bool pass = false;
	const float CesWireE2 = sp->GetCesWireE2();
	const float CesStripE2 = sp->GetCesStripE2();
	const float SinTheta = sp->GetSinTheta();
	const float EtCorr = sp->GetEtCorr();
	if (EtCorr<18) 
	{
		if ( ((CesWireE2 * SinTheta) < (0.14*EtCorr)) &&
				((CesStripE2 * SinTheta) < (0.14*EtCorr)) ) pass = true;
	} else 
	{
		if ( ((CesWireE2 * SinTheta) < (2.4+0.01*EtCorr))
				&& ((CesStripE2 * SinTheta) < (2.4+0.01*EtCorr)) ) pass = true;
	}

	//if (!pass) std::cout << "\tCESE2 failed  "<< std::endl;
	return pass;
}

void PhoJetsCount::FillPhotonHists(Histograms::PhotonHists_t& hist, const SuperPhoton* sp, const float fWeight)
{
	hist.Detector->Fill(sp->GetDetector(), fWeight);
	hist.DetEta->Fill(sp->GetDetEta(), fWeight);
	hist.DetPhi->Fill(sp->GetDetPhi(), fWeight);
	hist.EtCorr->Fill(sp->GetEtCorr(), fWeight);
	hist.XCes->Fill(sp->GetXCes(), fWeight);
	hist.ZCes->Fill(sp->GetZCes(), fWeight);
	hist.HadEm->Fill(sp->GetHadEm(), fWeight);
	hist.Chi2Mean->Fill(sp->GetChi2Mean(), fWeight);
	hist.N3d->Fill(sp->GetN3d(), fWeight);
	hist.Iso4->Fill(sp->GetIso4(), fWeight);
	hist.TrkPt->Fill(sp->GetTrkPt(), fWeight);
	hist.TrkIso->Fill(sp->GetTrkIso(), fWeight);
	hist.Ces2Wire->Fill(sp->GetCesWireE2(), fWeight);
	hist.Ces2Strip->Fill(sp->GetCesStripE2(), fWeight);
	hist.EmTime->Fill(sp->GetEmTime(), fWeight);
	hist.PhiWedge->Fill(sp->GetPhoton()->PhiSeedIndex(), fWeight);
}

