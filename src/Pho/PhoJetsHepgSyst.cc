/**********************************************************/
// Created this to do a study of the ISR/FSR and Q2
// systematics using only HEPG objetcs only . Sep 27,2009
/**********************************************************/

/*	
 *	$Id: PhoJetsHepgSyst.cc,v 1.3 2011/05/23 20:35:04 samantha Exp $
 *	$Log: PhoJetsHepgSyst.cc,v $
 *	Revision 1.3  2011/05/23 20:35:04  samantha
 *	ADDED:#include <assert.h>
 *	
 *	Revision 1.2  2009/11/19 20:48:16  samantha
 *	ADDED: A tree to hold all my final kinematic varibales. This is used
 *	to generate systemtic uncertainy for ISR/FSR/Q.
 *	
 *	Revision 1.1  2009/11/10 19:08:00  samantha
 *	First version used to study ISR/FSR/Q2 uncertainties. All results/plots were
 *	generated till NOv 10,2009 is from this version. This module is similar to
 *	Pho2JetsTemp but uses only HEPG objects. Looked at only the scale of
 *	systematic on photon Et plot only.
 *	
 */

#include "samantha/Pho/PhoJetsHepgSyst.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>
#include <iomanip>
#include <Stntuple/obj/TStnJet.hh>
#include <algorithm>
#include <assert.h>

ClassImp(PhoJetsHepgSyst)

//_____________________________________________________________________________
PhoJetsHepgSyst::PhoJetsHepgSyst(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  fPhoMaxEta(1.1),
  fPhoMinEta(0),
  fPhoMinEt(30.0),
  fJetMaxEta(3.0),
  fJetMinEt(15.0),
  bNoSummary(false),
  iJetCone(0)
{
	std::cout << "Hello I am PhoJetsHepgSyst module" << std::endl;
}

//_____________________________________________________________________________
PhoJetsHepgSyst::~PhoJetsHepgSyst() {
}

//_____________________________________________________________________________
void PhoJetsHepgSyst::SaveHistograms() {
}

//_____________________________________________________________________________
void PhoJetsHepgSyst::BookHistograms()
{
	DeleteHistograms();

	TFolder *new_folder, *sub_folder;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");				// 0
	vCounterHistLabels.push_back("Events Passed");					// 1
	vCounterHistLabels.push_back("All Tight #gamma Events");		// 2
	vCounterHistLabels.push_back("One tight #gamma only Events");	// 3
	vCounterHistLabels.push_back("Two tight #gamma only events");	// 4
	vCounterHistLabels.push_back("One tight #gamma and any Nloose pho events");	// 5
	HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	

	new_folder = GetHistFolder(this, "1Jet","#gamma + >=1 Jet case");
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt,"#gamma + >=1 Jets(15GeV)");
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h1_Pho,"#gamma + >=1 Jets(15GeV)");
	sub_folder = new_folder->AddFolder("LeadJet","#gamma Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet,"#gamma + >=1 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet,"#gamma + >=1 Jets(15GeV)");


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
	HistoMan.GetPhoton2JetsHistograms(sub_folder,h2_PhoJets,"#gamma + >=2 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
	HistoMan.GetTwoJetsHistograms(sub_folder,h2_TwoJets,"#gamma + >=2 Jets(15GeV)");
	*/
}


//_____________________________________________________________________________
int PhoJetsHepgSyst::BeginJob()
{
  	//evtPropMod = (EventProperties*) ((TStnAna*) GetAna()->GetModule("EventProperties"));
	//if (!evtPropMod) {
	//	StdOut(__FILE__,__LINE__,2,"EventProperties tag module not found.");
	//	bPermitRun = false;
	//}


				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");
	//RegisterDataBlock("PROD@JetCluModule-par-cone0.4","TStnJetBlock",&fJetBlockHadClu04);
	if (iJetCone == 0)
	{
		RegisterDataBlock("PROD@JetCluModule-had-cone0.4","TStnJetBlock",&fJetBlockHadClu04);
	} else if (iJetCone == 1)
	{
		RegisterDataBlock("PROD@JetCluModule-had-cone0.7","TStnJetBlock",&fJetBlockHadClu04);
	} else if (iJetCone == 2)
	{
		RegisterDataBlock("PROD@JetCluModule-had-cone1.0","TStnJetBlock",&fJetBlockHadClu04);
	}

	BookHistograms();

	// initialize vars
  
	tree = new TTree("PhoJetsHepgSyst_FlatTree","Holds final values of PhoJetsHepgSyst");
	TFolder *new_folder = GetHistFolder(this, "HistValuesTree","HistValuesTree");
	new_folder->Add(tree);

	tree->Branch("RunNumber",&aFinalVal.RunNumber, "RunNumber/I");
	tree->Branch("EvtNumber",&aFinalVal.EvtNumber, "EvtNumber/I");
	tree->Branch("pj1_Et_pho",&aFinalVal.pj1_Et_pho, "pj1_Et_pho/F");
	tree->Branch("pj1_Et_j1",&aFinalVal.pj1_Et_j1, "pj1_Et_j1/F");
	tree->Branch("pj1_InvM_pj1",&aFinalVal.pj1_InvM_pj1, "pj1_InvM_pj1/F");
	tree->Branch("pj1_Ht",&aFinalVal.pj1_Ht, "pj1_Ht/F");
	tree->Branch("pj1_Njet15",&aFinalVal.pj1_Njet15, "pj1_Njet15/I");
	tree->Branch("EtRatLeadJetLeadPho",&aFinalVal.EtRatLeadJetLeadPho, "EtRatLeadJetLeadPho/F");
	tree->Branch("DelRLeadPhoLeadJet",&aFinalVal.DelRLeadPhoLeadJet, "DelRLeadPhoLeadJet/F");
	tree->Branch("DelPhiLeadPhoLeadJet",&aFinalVal.DelPhiLeadPhoLeadJet, "DelPhiLeadPhoLeadJet/F");
	tree->Branch("pj2_Et_pho",&aFinalVal.pj2_Et_pho, "pj2_Et_pho/F");
	tree->Branch("pj2_Et_j1",&aFinalVal.pj2_Et_j1, "pj2_Et_j1/F");
	tree->Branch("pj2_Et_j2",&aFinalVal.pj2_Et_j2, "pj2_Et_j2/F");
	tree->Branch("pj2_InvM_pj1",&aFinalVal.pj2_InvM_pj1, "pj2_InvM_pj1/F");
	tree->Branch("pj2_InvM_pj2",&aFinalVal.pj2_InvM_pj2, "pj2_InvM_pj2/F");
	tree->Branch("pj2_InvM_pj1j2",&aFinalVal.pj2_InvM_pj1j2, "pj2_InvM_pj1j2/F");
	tree->Branch("pj2_InvM_j1j2",&aFinalVal.pj2_InvM_j1j2, "pj2_InvM_j1j2/F");
	tree->Branch("pj2_Ht",&aFinalVal.pj2_Ht, "pj2_Ht/F");
	

	return 0;
}

//_____________________________________________________________________________
int PhoJetsHepgSyst::BeginRun()
{

	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
	{
		StdOut(__FILE__,__LINE__,3,"Module is design to run only on MC data! exiting.");
	  	qMc = false;
		exit (1);
	}
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int PhoJetsHepgSyst::Event(int ientry)
{
  	SetPassed(0);
	hCount.iCounter->Fill(EVENTS_PROCESSED);
	//fHeaderBlock->Print();
	assert (fJetBlockHadClu04!=NULL && "jet block null");
	FillDataBlocks(ientry);
	//std::cout << "ientry, njet =" << ientry << ", " << fJetBlockHadClu04->NJets() << std::endl;
	
	TLorentzVector tlPhoVec = FindHEPGPar(fGenpBlock, 22,3,50,0);
	//std::cout << "Photon Pt = " << tlPhoVec.Pt() << std::endl;

	if (tlPhoVec.Pt() < GetPhoMinEt() || fabs(tlPhoVec.Eta()) > GetPhoMaxEta()) return 0;

	std::vector<TLorentzVector> vJets;  //store only jet with min pt and pass eta cut
	
	for (int ijet = 0; ijet < fJetBlockHadClu04->NJets(); ++ijet)
	{
		TStnJet *jet = fJetBlockHadClu04->Jet(ijet);
		//assert (jet != NULL && "jet is null");
		TLorentzVector jetvec = *(jet->Momentum());
		//std::cout << "i, Jet = " << ijet << ", " << jetvec.Pt() << std::endl;
		//remove the lead photon by mathching to GENP
		if (jetvec.DeltaR(tlPhoVec) < 0.1) continue;
		if (jetvec.Pt() > GetMinJetEt() && fabs(jetvec.Eta()) < GetMaxJetEta()) 
		{
			vJets.push_back(jetvec);
		}
	}
	//require at least 1 jet above threshold
	if (vJets.size() == 0) return 0;
	SetPassed(1);
	hCount.iCounter->Fill(TIGHT_PHO_EVENTS);
	

	//Sort jets in pt
	std::sort(vJets.begin(),vJets.end(),CompareTLorentzVectorPt);
	
	std::vector<TLorentzVector>::iterator it;
	/*fHeaderBlock->Print();
	for (it=vJets.begin(); it!=vJets.end(); ++it)
		std::cout << " " << (*it).Pt();
	std::cout << std::endl;
	*/

	
	Cleanup();		//reinit the tree branch values
	aFinalVal.RunNumber = fHeaderBlock->RunNumber();
	aFinalVal.EvtNumber = fHeaderBlock->EventNumber();
	aFinalVal.pj1_Et_pho = tlPhoVec.Pt();
	aFinalVal.pj1_Njet15    = vJets.size();
	
	if (vJets.size()>=1)
	{
		TLorentzVector tlLeadJetVec(0,0,0,0);
		tlLeadJetVec                   = vJets.at(0);
		aFinalVal.pj1_Et_j1            = tlLeadJetVec.Pt();
		TLorentzVector sum             = tlPhoVec + tlLeadJetVec;
		aFinalVal.pj1_InvM_pj1         = sum.M();
		aFinalVal.EtRatLeadJetLeadPho  = tlLeadJetVec.Pt()/tlPhoVec.Pt();
		aFinalVal.DelRLeadPhoLeadJet   = tlLeadJetVec.DeltaR(tlPhoVec);
		aFinalVal.DelPhiLeadPhoLeadJet = tlLeadJetVec.DeltaPhi(tlPhoVec);
		aFinalVal.pj1_Ht               = (GetHEPGHt(tlPhoVec,vJets));
	
		h1_Pho.EtCorr->Fill(tlPhoVec.Pt());
		h1_Pho.DetEta->Fill(tlPhoVec.Eta());
		h1_Jet.EtCorr->Fill(tlLeadJetVec.Pt());
		h1_Jet.DetEta->Fill(tlLeadJetVec.Eta());
		h1_PhoJet.InvMass->Fill(sum.M());
		h1_PhoJet.DelPhi->Fill(tlLeadJetVec.DeltaPhi(tlPhoVec));
		h1_PhoJet.DelR->Fill(tlLeadJetVec.DeltaR(tlPhoVec));
		h1_PhoJet.EtRatio->Fill(tlLeadJetVec.Pt()/tlPhoVec.Pt());

	}
	if (vJets.size()>=2)
	{
		TLorentzVector tlLeadJetVec(0,0,0,0),tlSubLeadJetVec(0,0,0,0);
		tlLeadJetVec = vJets.at(0);
		tlSubLeadJetVec = vJets.at(1);
		
		aFinalVal.pj2_Et_pho   = tlPhoVec.Pt();
		aFinalVal.pj2_Et_j1    = tlLeadJetVec.Pt();
		aFinalVal.pj2_Et_j2    = tlSubLeadJetVec.Pt();
		aFinalVal.pj2_InvM_pj1 = (tlPhoVec + tlLeadJetVec).M();
		aFinalVal.pj2_InvM_pj2 = (tlPhoVec + tlSubLeadJetVec).M();
		aFinalVal.pj2_InvM_pj1j2= (tlPhoVec + tlLeadJetVec + tlSubLeadJetVec).M();
		aFinalVal.pj2_InvM_j1j2 = (tlLeadJetVec + tlSubLeadJetVec).M();
		aFinalVal.pj2_Ht = (GetHEPGHt(tlPhoVec,vJets));

	}
	
	tree->Fill();
	
	//std::cout << "Photon Pt/Eta = " << tlPhoVec.Pt() << ", " << tlPhoVec.Eta() << std::endl;
	//std::cout << "Lead Jet      = " << tlLeadJetVec.Pt() << ", " << tlLeadJetVec.Eta() << std::endl;

	if (GetPassed()) {
		hCount.iCounter->Fill(EVENTS_PASSED);
	}
	
	return 0;

} // Event

/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void PhoJetsHepgSyst::Cleanup()
{
	aFinalVal.RunNumber = -999;
	aFinalVal.EvtNumber = -999;
	aFinalVal.pj1_Et_pho = -999;
	aFinalVal.pj1_Et_j1 = -999;
	aFinalVal.pj1_InvM_pj1 = -999;
	aFinalVal.pj1_Ht = -999;
	aFinalVal.pj1_Njet15 = -999;
	aFinalVal.EtRatLeadJetLeadPho = -999;
	aFinalVal.DelRLeadPhoLeadJet = -999;
	aFinalVal.DelPhiLeadPhoLeadJet = -999;
	aFinalVal.pj2_Et_pho = -999;
	aFinalVal.pj2_Et_j1 = -999;
	aFinalVal.pj2_Et_j2 = -999;
	aFinalVal.pj2_InvM_pj1 = -999;
	aFinalVal.pj2_InvM_pj2 = -999;
	aFinalVal.pj2_InvM_pj1j2 = -999;
	aFinalVal.pj2_InvM_j1j2 = -999;
	aFinalVal.pj2_Ht = -999;

} //Cleanup

/*-------------------------------------------------------------------*/
void PhoJetsHepgSyst::FillDataBlocks(int ientry)
{
	fGenpBlock->GetEntry(ientry);
	//fJetBlockClu04->GetEntry(ientry);
	fJetBlockHadClu04->GetEntry(ientry);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int PhoJetsHepgSyst::EndJob() {

	if (GetSummaryStat()) return 0;
		

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[PJS:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[PJS:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[PJS:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	sMsg += "[PJS:03:] Photon Min Et -------------- = " + ToStr(GetPhoMinEt()) + " GeV" + "\n";
	sMsg += "[PJS:04:] Photon Min<|DetEta|<Max ---- = " + ToStr(GetPhoMinEta()) + ", "+ ToStr(GetPhoMaxEta()) + "\n";
	sMsg += "[PJS:05:] Jet Min Et ----------------- = " + ToStr(GetMinJetEt()) + " GeV" + "\n";
	sMsg += "[PJS:06:] Jet |DetEta|<Max ----------- = " + ToStr(GetMaxJetEta()) + "\n";
	sMsg += "[PJS:10:] Photon+>=1 Jet Events ------ = " + ToStr(hCount.iCounter->GetBinContent(TIGHT_PHO_EVENTS+1)) + "\n";
	sMsg += "[PJS:11:] Jets' with cone size ------- = " + ToStr(GetJetCone()) + "\n";
	sMsg += "---------------------------------------------------\n";

	std::cout << sMsg;

	return 0;
}



void PhoJetsHepgSyst::SetPhoMaxEta(const float eta)
{
	if (fabs(eta) > 1.1)
	{
		StdOut(__FILE__,__LINE__,3,"Photon Detector eta require is > than what is returned in TStnPhoton::Detector() == 0??");
		exit (1);
	} else {
 		fPhoMaxEta = fabs(eta);
	}
}

void PhoJetsHepgSyst::SetPhoMinEta(const float eta)
{
	if (fabs(eta) > 1.1)
	{
		StdOut(__FILE__,__LINE__,3,"Photon Detector eta require is > than what is returned in TStnPhoton::Detector() == 0??");
		exit (1);
	} else {
 		fPhoMinEta = fabs(eta);
	}
}

void PhoJetsHepgSyst::SetPhoMinEt(const float et)
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

void PhoJetsHepgSyst::SetMaxJetEta(const float eta)
{
	fJetMaxEta = fabs(eta);
}

void PhoJetsHepgSyst::SetMinJetEt(const float et)
{ 
	if (et < 0)
	{
		StdOut(__FILE__,__LINE__,3,"Lead Jet Et must be > 0. Using default value.");
		exit(-1);
	} else
	{
		fJetMinEt = et;
	}
}

//return sumPt of lead photon and the jet with Et>15GeV
float PhoJetsHepgSyst::GetHEPGHt(const TLorentzVector tlPho, const std::vector<TLorentzVector>& vJet)
{
	float sumPt = 0;
	sumPt += tlPho.Pt();
	std::vector<TLorentzVector>::const_iterator it;
	
	for (it = vJet.begin(); it!= vJet.end(); ++it)
	{
		sumPt += (*it).Pt();
	}
	return sumPt;
	
}
	
