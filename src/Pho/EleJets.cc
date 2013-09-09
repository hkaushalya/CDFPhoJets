//////////////////////////////////////////////////////////////////////
// This is identical to Pho2JetTemp. Here we look for tight electron//
// +jets.                                                           //
//////////////////////////////////////////////////////////////////////

/*
 *	$Log: EleJets.cc,v $
 *	Revision 1.3  2011/05/25 20:52:12  samantha
 *	ADDED: 1. To weight the EWK MC with nvtx weights, SetNvtxWeights().
 *	2. Included the optional weights when filling hists.
 *	
 *	Revision 1.2  2010/06/28 22:21:50  samantha
 *	MODIFIED: this to check how well the EWK MC represents e+jet events in data.
 *	We are ignoring m+jet and t+jet. I only try to ID electrons. See the results of
 *	this version in elog:1678. I have looked at e+jet(15GeV) events and rejected all
 *	events where the MET is aligned with a jet (Et>5GeV) or the electron.
 *	I have added many setters so this is similar to PhoJets and some histos.
 *	
 *	Revision 1.1  2009/05/22 17:06:37  samantha
 *	Identical to Pho2JetsTemp but looks for Tight Ele+Jets.
 *	
 *
 */

#include "samantha/Pho/EleJets.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>
#include <iomanip>
#include <assert.h>

ClassImp(EleJets)

//_____________________________________________________________________________
EleJets::EleJets(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bPermitRun(true),
  fEmTimeMin(-4.8),
  fEmTimeMax(4.8),
  fEleMaxDetEta(1.1),
  fEleMinDetEta(0.),
  fEleMinEt(30),
  fEleMaxEt(1500),
  bRemConvEles(0),
  bNoSummary(false)
{
	std::cout << "Hello I am EleJets module" << std::endl;
}

//_____________________________________________________________________________
EleJets::~EleJets() {
}

//_____________________________________________________________________________
void EleJets::SaveHistograms() {
}

//_____________________________________________________________________________
void EleJets::BookHistograms()
{
	DeleteHistograms();

  	//char name [200];
  	//char title[200];
	TFolder *new_folder, *sub_folder;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");				// 0
	vCounterHistLabels.push_back("Events Passed");					// 1
	vCounterHistLabels.push_back("e+>=1 Jets Events");		// 2
	vCounterHistLabels.push_back("e+>=2 Jets Events");		// 3
	vCounterHistLabels.push_back("e Events");					// 4
	HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	

	new_folder = GetHistFolder(this, "1Jet","e + >=1 Jet case");
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt,"e + >=1 Jets(15GeV)");

	
	metb4 = dynamic_cast<TH1F*> (h1_Evt.Met->Clone("Met_b4cut"));
	meta4 = dynamic_cast<TH1F*> (h1_Evt.Met->Clone("Met_a4cut"));
	sub_folder->Add(metb4);
	sub_folder->Add(meta4);
	Nvtx12_a4 = dynamic_cast<TH1F*> (h1_Evt.N12Vertices->Clone("Nvtx12_a4cut"));
	sub_folder->Add(Nvtx12_a4);
	
	sub_folder = new_folder->AddFolder("Photon","e Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h1_Pho,"e + >=1 Jets(15GeV)");


	phoEt_a4 = dynamic_cast<TH1F*> (h1_Pho.EtCorr->Clone("PhoEtCorr_a4"));
	sub_folder->Add(phoEt_a4);
	phometdelphi_a4 = dynamic_cast<TH1F*> (h1_Pho.PhoMetDelPhi->Clone("PhoMetDelPhi_a4"));
	sub_folder->Add(phometdelphi_a4);

	
	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet,"e + >=1 Jets(15GeV)");

	jetEt_a4 = dynamic_cast<TH1F*> (h1_Jet.EtCorr->Clone("JetEtCorr_a4"));
	sub_folder->Add(jetEt_a4);
	jetmetdelphi_a4 = dynamic_cast<TH1F*> (h1_Jet.JetMetDelPhi->Clone("JetMetDelPhi_a4"));
	sub_folder->Add(jetmetdelphi_a4);
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","e and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet,"e + >=1 Jets(15GeV)");
	phoLeadJetInvMass_a4 = dynamic_cast<TH1F*> (h1_PhoJet.InvMass->Clone("PhoJetMass_a4"));
	phoMetTrnMass_b4 = dynamic_cast<TH1F*> (h1_PhoJet.InvMass->Clone("PhoMetTrnMass_b4"));
	phoMetTrnMass_a4 = dynamic_cast<TH1F*> (h1_PhoJet.InvMass->Clone("PhoMetTrnMass_a4"));

	sub_folder->Add(phoLeadJetInvMass_a4);
	sub_folder->Add(phoMetTrnMass_b4);
	sub_folder->Add(phoMetTrnMass_a4);

/*
	new_folder = GetHistFolder(this, "2Jet","e + >=2 Jets case");
	
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h2_Evt,"e + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon","e Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h2_Pho,"e + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet1,"e + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet2,"e + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","e and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet1,"e + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2ndLeadJet","e and 2nd Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet2,"e + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2Jets","e + 2Jets Histograms");
	HistoMan.GetPhoton2JetsHistograms(sub_folder,h2_PhoJets,"e + >=2 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
	HistoMan.GetTwoJetsHistograms(sub_folder,h2_TwoJets,"e + >=2 Jets(15GeV)");
	*/
}


//_____________________________________________________________________________
int EleJets::BeginJob()
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

//  	evtPropMod = (EventProperties*) ((TStnAna*) GetAna()->GetModule("EventProperties"));
//	if (!evtPropMod) {
//		StdOut(__FILE__,__LINE__,2,"EventProperties tag module not found.");
//		bPermitRun = false;
//	}

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
		StdOut(__FILE__,__LINE__,3,"TriggerModule required!.");
		bPermitRun = false;
	}


				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();

	// initialize vars
	if (GetUseNvtxWeigts())
	{
		SetNvtxWeights();
	}
  

	return 0;
}

//_____________________________________________________________________________
int EleJets::BeginRun()
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
int EleJets::Event(int ientry)
{
  	SetPassed(0);

	// this makes sure that all required mods exists
	if (!bPermitRun) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
		return 0;
	}

	hCount.iCounter->Fill(0);
	std::vector<int> vEleIndex;
	bool bSkipThisEvent = false;

  	int NsuperPho = initSpMod->GetSuperPhoSize();
	for (int i=0; i < NsuperPho; i++) {
		if (fabs(initSpMod->GetSuperPhoton(i)->GetDetEta()) > GetEleMaxDetEta()) continue;
		if (initSpMod->GetSuperPhoton(i)->GetEtCorr() < GetEleMinEt() 
			 || initSpMod->GetSuperPhoton(i)->GetEtCorr() > GetEleMaxEt()) continue;
		if (GetRemoveConvEles() && initSpMod->GetSuperPhoton(i)->IsConversion()) continue;
		if (! initSpMod->GetSuperPhoton(i)->IsTightElectron()) continue;		//pho-like electrons
		if (GetHeaderBlock()->McFlag() == false)  
		{
			float fEmTime = initSpMod->GetSuperPhoton(i)->GetEmTime();
			if (fEmTime < fEmTimeMin || fEmTime > fEmTimeMax) continue;
		}
		vEleIndex.push_back(i);
		break;
	}


	if (!bSkipThisEvent) {
		if (vEleIndex.size() > 0) {
			FillPhotonIDHist(vEleIndex.at(0));		//Fill all the ID information of my photon for debugging
																// this way I can know if my selection is correct.
			hCount.iCounter->Fill(4);
			DoJetSelection(vEleIndex.at(0));
			SetPassed(1);
		}
	}

	if (GetPassed()) {
		hCount.iCounter->Fill(1);
	}
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
void EleJets::FillPhotonIDHist(const int iSpInd)
{
	if (iSpInd < 0 || iSpInd > initSpMod->GetSuperPhoSize())
	{
		StdOut(__FILE__,__LINE__,3,"Invalid SuperPhoton index!");
		return;
	}
	
	int TightPhoId = initSpMod->GetSuperPhoton(iSpInd)->GetTightPhotonId();
	int LoosePhoId = initSpMod->GetSuperPhoton(iSpInd)->GetLoosePhotonId();
	int PhoLikeEleId = initSpMod->GetSuperPhoton(iSpInd)->GetTightElectronId();
	int StdLooseEleId = initSpMod->GetSuperPhoton(iSpInd)->GetStdLooseElectronId();
	int StdTightEleId = initSpMod->GetSuperPhoton(iSpInd)->GetStdTightElectronId();
	int HaloId = initSpMod->GetSuperPhoton(iSpInd)->GetBeamHaloId();
	int PhoenixId = initSpMod->GetSuperPhoton(iSpInd)->GetPhoenixId();
	int ConversionId = initSpMod->GetSuperPhoton(iSpInd)->GetConversionId();

	if (TightPhoId == 0) h1_Pho.PhoIDs->Fill(0);
	if (LoosePhoId == 0) h1_Pho.PhoIDs->Fill(1);
	if (TightPhoId>0 && LoosePhoId == 0) h1_Pho.PhoIDs->Fill(2);
	if (PhoLikeEleId == 0)	h1_Pho.PhoIDs->Fill(3);
	if (StdTightEleId == 0) h1_Pho.PhoIDs->Fill(4);
	if (StdLooseEleId == 0) h1_Pho.PhoIDs->Fill(5);
	if (PhoenixId > 0) h1_Pho.PhoIDs->Fill(6);
	if (ConversionId > 0) h1_Pho.PhoIDs->Fill(7);
	if (HaloId > 0) h1_Pho.PhoIDs->Fill(8);
	
}


			
/*-------------------------------------------------------------------*/
void EleJets::DoJetSelection(int iSpIndex)
{
		
	//find the leading detector jet after removing the leading photon
	TLorentzVector tlPhoVec = initSpMod->GetSuperPhoton(iSpIndex)->GetCorVec();
	TVector2 tv2MetVec(jetMod->GetMyMetXCorr(0,3), jetMod->GetMyMetYCorr(0,3));
	
   TLorentzVector tlMetVec(tv2MetVec.Px(), tv2MetVec.Py(), 0, tv2MetVec.Mod());

	float fWgt = 1;
	if (GetUseNvtxWeigts())
	{
		const int iNvtx12 = trigMod->GetN12vertex();
		if (iNvtx12 <= vNvtxWgts.size()) fWgt = vNvtxWgts.at(iNvtx12 - 1);
	}
	

	
	
	if (jetMod->GetMyNjet(0,3) >= 1)
	{	// mono-jet case
		//hPhoLeadJetPtRatio->Fill(tlLeadJetVec.Pt()/tlPhoVec.Pt() - 1);
		//std::cout << "========= "; fHeaderBlock->Print();
		const float dphi_pm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlPhoVec.Phi())));
			
		//std::cout << "pho pt, met, delphi = " << tlPhoVec.Pt() << ", " << tv2MetVec.Mod() << ", " <<  dphi_pm << std::endl;

		float dphi_jm_closest = 10.0;
		float dphi_jm_furthest = 0.0;
		TLorentzVector tlVecClosestJet(0,0,0,0);
		for (int ijet=0; ijet<jetMod->GetMyNjet(0,0); ++ijet)
		{
			TLorentzVector tlJetVec = *(jetMod->GetMyJetNoEMobj_lev6(0,ijet));
			if (tlJetVec.Pt()<5.0)
			{
				//std::cout <<  ijet << "\t"<< tlJetVec.Pt() << "\t low pt, ignored!" << std::endl;
				continue;
			}
			const float dphi_jm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlJetVec.Phi())));
			//std::cout << "dphi_jm [" << dphi_jm << "]" << std::endl;

			if (dphi_jm < dphi_jm_closest)
			{
				dphi_jm_closest = dphi_jm;
			   tlVecClosestJet = tlJetVec;
			}

			/*if (dphi_jm <0.4 || dphi_jm>2.8) 
			  {
			  std::cout << " FOUND A CLOSE JET dphi_jm [" << dphi_jm << "]" << std::endl;
			  if (dphi_jm < 0.4 || dphi_jm > (TMath::Pi() - 0.4)) 
			  {
			  std::cout << "\t******************** DelPhi<0.4 || Delphi>2.8 (" << dphi_jm << ") " << std::endl;

			  if (dphi_jm < dphi_jm_closest || (TMath::Pi()-dphi_jm) < dphi_jm_closest) 
			  {
			  if (dphi_jm < dphi_jm_closest) dphi_jm_closest = dphi_jm;
			  else if ((TMath::Pi()-dphi_jm)< dphi_jm_closest) dphi_jm_closest = TMath::Pi() - dphi_jm;
			  tlVecClosestJet = tlJetVec;
			  std::cout << "\t******************** ijet/pt/dphi_jm/dphi_jm_closest = " 
			  << ijet << "\t"<< tlJetVec.Pt() << "\t" << dphi_jm << " -> " <<  dphi_jm_closest << std::endl;
			  } else 
			  {
			//	std::cout <<  ijet << "\t"<< tlJetVec.Pt() << "\t" << dphi_jm  << std::endl;
			}
			}
			}
			*/
		}


		//before delphi cuts
		TLorentzVector tlLeadJetVec = *(jetMod->GetMyJetNoEMobj_lev6(0,0));
		const float fEleMet_TrnMass = (tlPhoVec + tlMetVec).M();
		h1_Pho.EtCorr->Fill(tlPhoVec.Pt());
		h1_Pho.PhoMetDelPhi->Fill(dphi_pm);
		h1_Jet.JetMetDelPhi->Fill(dphi_jm_closest);
		h1_Jet.EtCorr->Fill(tlVecClosestJet.Pt());
		metb4->Fill(tv2MetVec.Mod());
		h1_PhoJet.InvMass->Fill((tlPhoVec + tlLeadJetVec).M());
		phoMetTrnMass_b4->Fill(fEleMet_TrnMass);
		hCount.iCounter->Fill(2);
		h1_Evt.N12Vertices->Fill(trigMod->GetN12vertex());

		//if there is a jet/photon aligned with MET reject those events
/*		if (dphi_jm_closest < 10 
				&& (dphi_jm_closest > 0.4 && dphi_jm_closest < (TMath::Pi() - 0.4))
				//&& (dphi_pm > 0.4) && (dphi_pm < (TMath::Pi() - 0.4))
				&& fEleMet_TrnMass > 50.0
				)
*/
		if (fEleMet_TrnMass > 50.0)
		{
			//std::cout << " $$$$$$$$$$$$$$$$$$$$$$ Delphi > 0.4 (" << dphi_jm_closest << ")" << std::endl;
			phoEt_a4->Fill(tlPhoVec.Pt(), fWgt);
			//jetEt_a4->Fill(tlVecClosestJet.Pt(), fWgt);
			meta4->Fill(tv2MetVec.Mod(), fWgt);
			jetmetdelphi_a4->Fill(dphi_jm_closest, fWgt);
			phometdelphi_a4->Fill(dphi_pm, fWgt);

			if (tlLeadJetVec.Pt() <15.) 
			{
				std::cout << " ############ JET WITH LOW pT " <<  tlLeadJetVec.Pt() << std::endl;
				GetHeaderBlock()->Print();
				assert( false);

			}
			jetEt_a4->Fill(tlVecClosestJet.Pt(), fWgt);
			phoLeadJetInvMass_a4->Fill((tlPhoVec + tlLeadJetVec).M(), fWgt);
			h1_PhoJet.DelPhi->Fill(fabs(tlPhoVec.DeltaPhi(tlLeadJetVec)), fWgt);
			phoMetTrnMass_a4->Fill(fEleMet_TrnMass, fWgt);
			Nvtx12_a4->Fill(trigMod->GetN12vertex(), fWgt);
			
		} else 
		{
			//std::cout << " xxxxxxxxxxxxxxxxxxxxx Delphi < 0.4 (" << dphi_jm_closest << ")" << std::endl;
		}

		
	}

	
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
void EleJets::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void EleJets::FillDataBlocks(int ientry)
{
}

/*-------------------------------------------------------------------*/
void EleJets::SetNvtxWeights()
{

	vNvtxWgts.clear();

	if (iDataSample == 20) // Zee MC
	{
		vNvtxWgts.push_back(0.542458);
		vNvtxWgts.push_back(1.00602);
		vNvtxWgts.push_back(1.68042);
		vNvtxWgts.push_back(2.5309);
		vNvtxWgts.push_back(3.46027);
		vNvtxWgts.push_back(4.77576);
		vNvtxWgts.push_back(12.6794);
		vNvtxWgts.push_back(16.9325);
		vNvtxWgts.push_back(13.886);

	} else if (iDataSample == 21) //Zmm MC
	{
		vNvtxWgts.push_back(0.527496);
		vNvtxWgts.push_back(1.15733);
		vNvtxWgts.push_back(1.56931);
		vNvtxWgts.push_back(1.60449);
		vNvtxWgts.push_back(5.4945);
	} else if (iDataSample == 22) //Ztt MC
	{
		vNvtxWgts.push_back(0.467678);
		vNvtxWgts.push_back(1.09859);
		vNvtxWgts.push_back(2.15647);
		vNvtxWgts.push_back(3.32276);
		vNvtxWgts.push_back(4.79241);
		vNvtxWgts.push_back(10.4261);
		vNvtxWgts.push_back(10.8621);

	} else if (iDataSample == 30) //Wen MC)
	{
		vNvtxWgts.push_back(0.496098);
		vNvtxWgts.push_back(1.03245);
		vNvtxWgts.push_back(1.96743);
		vNvtxWgts.push_back(3.23924);
		vNvtxWgts.push_back(4.76124);
		vNvtxWgts.push_back(8.08731);
		vNvtxWgts.push_back(16.4894);
		vNvtxWgts.push_back(60.4164);
		vNvtxWgts.push_back(176.51);
	} else if (iDataSample == 31) //Wmn MC
	{
		vNvtxWgts.push_back(0.569584);
		vNvtxWgts.push_back(0.823644);
		vNvtxWgts.push_back(2.2337);
		vNvtxWgts.push_back(2.53752);
	} else if (iDataSample == 32) //Wtn MC
	{
		vNvtxWgts.push_back(0.51911);
		vNvtxWgts.push_back(1.04513);
		vNvtxWgts.push_back(1.83233);
		vNvtxWgts.push_back(2.5778);
		vNvtxWgts.push_back(2.84779);
		vNvtxWgts.push_back(3.36053);
		vNvtxWgts.push_back(7.3911);
		vNvtxWgts.push_back(7.64155);
	} else {
		assert (false && "Unknown sample, No Nvtx weights found!");
	}
		

}



//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int EleJets::EndJob() {

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[ELJ:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[ELJ:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[ELJ:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	sMsg += "[ELJ:03:] Ele + >=1 Jets Events ------ = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
	sMsg += "[ELJ:04:] Ele + >=2 Jets Events ------ = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n";
	sMsg += "[ELJ:05:] Electron Events ------------ = " + ToStr(hCount.iCounter->GetBinContent(5)) + "\n";
	sMsg += "[ELJ:06:] Electron Min<|DetEta|<Max -- = " + ToStr(GetEleMinDetEta()) + " - "+ ToStr(GetEleMaxDetEta()) + "\n";
	sMsg += "[ELJ:07:] Electron Min/Max Et -------- = " + ToStr(GetEleMinEt())  + " - " + ToStr(GetEleMaxEt()) + " GeV" + "\n";
	sMsg += "[ELJ:08:] Remove Conversion Electrons? = ";
	if (GetRemoveConvEles()) sMsg += "YES\n";
	else sMsg += "NO\n";
	if (GetHeaderBlock()->McFlag())
	sMsg += "[ELJ:09:] EM Time window used -------- = NO CUT (MC)\n";
	else
	sMsg += "[ELJ:09:] EM Time window used -------- = " + ToStr(GetEmTimeWindowMin()) + "," + ToStr(GetEmTimeWindowMax())+ "\n";
	sMsg += "[ELJ:10:] Uset Nvtx Wegihts?  -------- = " + ToStr(GetUseNvtxWeigts()) + " (1==YES)\n";
	sMsg += "[ELJ:11:] iDataSample ---------------- = " + ToStr(GetDataSample()) + " (20=Zee,21=Zmm,22=Ztt,30=Wen,31=Wmn,32=Wtn)\n";
	sMsg += "---------------------------------------------------\n";

	std::cout << sMsg;
	return 0;
}


void EleJets::SetEleMaxDetEta(const float eta)
{
	std::cout << "eta = " << fabs(eta) << std::endl;
	if (fabs(eta) > 1.1001)		//1.101 is a hack. for some reason when I require
	{									// eta=1.1, and compare (>) to 1.1, this condition becomes true!
										// sam -11-12-2009
		StdOut(__FILE__,__LINE__,3,"Electron Detector eta you required is > than what is returned in TStnPhoton::Detector() == 0??");
		exit (1);
	} else {
 		fEleMaxDetEta = fabs(eta);
	}
}

void EleJets::SetEleMinDetEta(const float eta)
{
	if (fabs(eta) > 1.1)
	{
		StdOut(__FILE__,__LINE__,3,"Electron Detector eta require is > than what is returned in TStnPhoton::Detector() == 0??");
		exit (1);
	} else {
 		fEleMinDetEta = fabs(eta);
	}
}

void EleJets::SetEleMinEt(const float et)
{ 
	if (et < 0)
	{
		StdOut(__FILE__,__LINE__,3,"Electron Et must be > 0. Using default value.");
		exit(-1);
	} else
	{
		fEleMinEt = et;
	}
}
void EleJets::SetEleMaxEt(const float et)
{ 
	if (et < 0)
	{
		StdOut(__FILE__,__LINE__,3,"Electron Et must be > 0. Using default value.");
		exit(-1);
	} else
	{
		fEleMaxEt = et;
	}
}
