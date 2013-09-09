///////////////////////////////////////////////////////////
// This will hold common event properties like MEt, etc. //
// This will have function/s to fill a given Event Hists //
// derived from Histograms.cc/hh.                        //
// This talks to the JetFilterModule for correct MEt and //
// SumEt etc.                                            //
// All other mods can access these info and can get their//
// own histograms filled by functions in here.           //
//                                                       //
// now this can fill properties of a SuperPhoton too.    //
// histos must be derived from Histograms::EmObjHists_t  //
// and pass in the SuperPhoton for this.                 //
// 01-27-2008                                            //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "samantha/Pho/EventProperties.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnJet.hh"
#include "samantha/utils/FreeFunctions.hh"

ClassImp(EventProperties)

//_____________________________________________________________________________
EventProperties::EventProperties(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bRunPermit(true),
  bNoSummary(false)
{
	std::cout << "Hello I am EventProperties module" << std::endl;
}

//_____________________________________________________________________________
EventProperties::~EventProperties() {
}

//_____________________________________________________________________________
void EventProperties::SaveHistograms() {
}

//_____________________________________________________________________________
void EventProperties::BookHistograms()
{
	DeleteHistograms();
	TFolder *new_folder;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");				// 0
	vCounterHistLabels.push_back("Events Passed");					// 1
	HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	
}


//_____________________________________________________________________________
int EventProperties::BeginJob()
{
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}
  	jetMod = (JetFilterModuleV3*) ((TStnAna*) GetAna()->GetModule("JetFilterV3"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		bRunPermit = false;
	}
  	trigMod = (TriggerModule*) ((TStnAna*) GetAna()->GetModule("TriggerModule"));
	if (!trigMod) {
		StdOut(__FILE__,__LINE__,3,"Trigger module required!.");
		bRunPermit = false;
	}
	
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fJetBlock 		= (TStnJetBlock*)      RegisterDataBlock("JetBlock","TStnJetBlock");

	BookHistograms();

	// initialize vars
	counter.evtsProcessed 	= 0;
	counter.evtsPassed 		= 0;

	return 0;
}

//_____________________________________________________________________________
int EventProperties::BeginRun()
{
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int EventProperties::Event(int ientry)
{
	SetPassed(1);
	hCount.iCounter->Fill(0);
  	if (!bRunPermit) {
		SetPassed(0);
		StdOut(__FILE__,__LINE__,3,"RunPermit failed. One or more pre-req. mods not found!.");
		exit (1);
		return 0;
	}
	
	Cleanup();
	FillDataBlocks(ientry);
	
	aEvtProp.MetCorr 		= jetMod->GetMyMetCorr(0,3);			// cone 0.4 ,Et>15GeV
	aEvtProp.SumetCorr 	= jetMod->GetMySumEtCorr(0,3);
	aEvtProp.HtCorr 		= jetMod->GetMyHtCorr(0,3);
	aEvtProp.NVertices	= trigMod->GetNvertex();
	aEvtProp.N12Vertices	= trigMod->GetN12vertex();
	aEvtProp.NJet15		= jetMod->GetMyNjet(0,3);
	aEvtProp.BestVertexZ	= trigMod->GetBestVertexZ();
	aEvtProp.RunNumber	= fHeaderBlock->RunNumber();
	aEvtProp.EvtNumber	= fHeaderBlock->EventNumber();

	if (GetPassed()) {
		hCount.iCounter->Fill(1);
	}
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void EventProperties::Cleanup()
{
	// Clear or reinitialize data memebers
	aEvtProp.MetCorr		= 0;
	aEvtProp.SumetCorr	= 0;
	aEvtProp.NVertices	= 0;
	aEvtProp.N12Vertices	= 0;
	aEvtProp.NJet15		= 0;
	aEvtProp.RunNumber	= 0;
	aEvtProp.EvtNumber	= 0;
	aEvtProp.BestVertexZ	= -99999.99;
	
} //Cleanup


/*-------------------------------------------------------------------*/
void EventProperties::FillDataBlocks(int ientry)
{
	// Fill all data blocks needed.
	fJetBlock->GetEntry(ientry);		// need to get Jet tower IEta/IPhi
}

//____________________________________________________________________________________________________
void EventProperties::FillEventHists(Histograms::EventHists_t& hist, const float fWeight)
{ 
	// public method to be used by other mods to fill event properties histos
	hist.Met->Fill(aEvtProp.MetCorr, fWeight);
	hist.Sumet->Fill(aEvtProp.SumetCorr, fWeight);
	hist.Ht->Fill(aEvtProp.HtCorr, fWeight);
	hist.NVertices->Fill(aEvtProp.NVertices);
	hist.N12Vertices->Fill(aEvtProp.N12Vertices);
	hist.NJet15->Fill(aEvtProp.NJet15);
	hist.MetRun->Fill(aEvtProp.RunNumber, aEvtProp.MetCorr);
	hist.SumetRun->Fill(aEvtProp.RunNumber, aEvtProp.SumetCorr);
	hist.NVerticesRun->Fill(aEvtProp.RunNumber, aEvtProp.NVertices);
	hist.N12VerticesRun->Fill(aEvtProp.RunNumber, aEvtProp.N12Vertices);
	hist.NJet15Run->Fill(aEvtProp.RunNumber, aEvtProp.NJet15);
	hist.NJet15VsMet->Fill(aEvtProp.MetCorr, aEvtProp.NJet15);
	hist.NJet15VsSumet->Fill(aEvtProp.SumetCorr, aEvtProp.NJet15);
	hist.N12VerticesVsMet->Fill(aEvtProp.MetCorr, aEvtProp.N12Vertices);
	hist.N12VerticesVsSumet->Fill(aEvtProp.SumetCorr, aEvtProp.N12Vertices);
	hist.MetSumetRatio->Fill(aEvtProp.MetCorr/aEvtProp.SumetCorr, fWeight);
	hist.BestVertexZ->Fill(aEvtProp.BestVertexZ, fWeight);
}

//____________________________________________________________________________________________________
void EventProperties::FillPhotonHists(Histograms::PhotonHists_t& hist, int iSuperPhoIndex, float fWght)
{
	//templates pass in the index of the SuperPhoton and the hists to filled.
	// this is basically to look at the photon in each of my templates.

	hist.Detector->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetDetector(),fWght);
	hist.DetEta->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetDetEta(),fWght);
	hist.DetPhi->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetDetPhi(),fWght);
	hist.EtCorr->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetEtCorr(),fWght);
	hist.XCes->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetXCes(),fWght);
	hist.ZCes->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetZCes(),fWght);
	hist.HadEm->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetHadEm(),fWght);
	hist.Chi2Mean->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetChi2Mean(),fWght);
	hist.N3d->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetN3d(),fWght);
	hist.Iso4->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetIso4(),fWght);
	hist.TrkPt->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetTrkPt(),fWght);
	hist.TrkIso->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetTrkIso(),fWght);
	hist.Ces2Wire->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetCesWireE2(),fWght);
	hist.Ces2Strip->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetCesStripE2(),fWght);
	hist.EmTime->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetEmTime(),fWght);
	hist.EmTimeVsRun->Fill(GetHeaderBlock()->RunNumber(), initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetEmTime());
	hist.EtCorrVsRun->Fill(GetHeaderBlock()->RunNumber(), initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetEtCorr());
	hist.PhiWedge->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetPhoton()->PhiSeedIndex(),fWght);

	
	float fPhoPhi = initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetDetPhi();
	float fMetPhi = jetMod->GetMyMetPhiCorr(0,3);		//goes from 0-2pi
	float fPhoMetDelPhi = DelPhi(fPhoPhi, fMetPhi);				//this is a free function
	hist.PhoMetDelPhi->Fill(fabs(fPhoMetDelPhi), fWght);		// 0 to 2pi
	
}

//____________________________________________________________________________________________________
void EventProperties::FillJetHists(Histograms::JetHists_t& hist, int iJetIndex, float fWght)
{
	// iIndex is the position of the jet in the JetFilterV2 jet list.
	if (iJetIndex <0 || iJetIndex > jetMod->GetMyNjet(0,3)) {
		StdOut(__FILE__,__LINE__,3,"Inavalid Jet index. Cannot find Jet in JetFilterModuleV3!");
		return;
	}

	int iJetCone4 = 0;			//get cone 0.4 jets
	
	hist.Emfr->Fill(jetMod->GetMyJetEmFrCorr(iJetCone4,iJetIndex), fWght);
	hist.DetEta->Fill(jetMod->GetMyJetEtaDetCorr(iJetCone4,iJetIndex), fWght);
	hist.NTracks->Fill(jetMod->GetMyJetNtrk(iJetCone4,iJetIndex), fWght);
	hist.NTowers->Fill(jetMod->GetMyJetNtwr(iJetCone4,iJetIndex), fWght);
	TLorentzVector tlVec = *(jetMod->GetMyJetNoEMobj_lev6(iJetCone4,iJetIndex));
	hist.EtCorr->Fill(tlVec.Energy() * TMath::Sin(tlVec.Theta()), fWght);
	hist.EvtPhi->Fill(fabs(tlVec.Phi()), fWght);
	
	//get info not provided by JetFilterModule
	int iJetOrigIndex = jetMod->GetMyJetBlockInd(0,iJetIndex);
	TStnJet* jet = fJetBlock->Jet(iJetOrigIndex);
	hist.SeedIEta->Fill(jet->SeedIEta(), fWght);
	hist.SeedIPhi->Fill(jet->SeedIPhi(), fWght);

	float fMetPhi = jetMod->GetMyMetPhiCorr(0,3);		//goes from 0-2pi
	float fJetMetDelPhi = DelPhi(tlVec.Phi(), fMetPhi);				//this is a free function
	hist.JetMetDelPhi->Fill(fabs(fJetMetDelPhi), fWght);		// 0 to 2pi

	
}


//____________________________________________________________________________________________________
void EventProperties::FillPhoton1JetHists(Histograms::Photon1JetHists_t& hist,
								int iSuperPhoIndex, int iJetIndex, const float fWeight)
{
	int iJetCone4 = 0;
	TLorentzVector tlPho =  initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetCorVec();
	TLorentzVector tlJet =  *(jetMod->GetMyJetNoEMobj_lev6(0,iJetIndex));
	
	TLorentzVector tlSum(0,0,0,0);
	tlSum = tlPho + tlJet; 
	float fEtRatio = (tlJet.E() * TMath::Sin(tlJet.Theta())) / (tlPho.E() * TMath::Sin(tlPho.Theta()));
	
	hist.InvMass->Fill(tlSum.M(), fWeight);
	hist.DelPhi->Fill(fabs(tlPho.DeltaPhi(tlJet)), fWeight);
	hist.DelEta->Fill(fabs(tlPho.Eta() - tlJet.Eta()), fWeight);
	hist.DelEtaPlus->Fill(fabs(tlPho.Eta() + tlJet.Eta()), fWeight);
	hist.DelR->Fill(tlPho.DeltaR(tlJet), fWeight);
	hist.EtRatio->Fill(fEtRatio, fWeight);
	hist.JetEmfrVsDelPhi->Fill(fabs(tlPho.DeltaPhi(tlJet))
								, jetMod->GetMyJetEmFrCorr(iJetCone4,iJetIndex), fWeight);
	hist.JetEmfrVsPhoPhiWedge->Fill(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetPhoton()->PhiSeedIndex()
								, jetMod->GetMyJetEmFrCorr(iJetCone4,iJetIndex), fWeight);

	hist.DelEtaDet->Fill(fabs(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetDetEta() - jetMod->GetMyJetEtaDetCorr(0,iJetIndex)), fWeight);
	hist.DelEtaDetPlus->Fill(fabs(initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetDetEta() + jetMod->GetMyJetEtaDetCorr(0,iJetIndex)), fWeight);
}

//____________________________________________________________________________________________________
void EventProperties::FillPhoton2JetsHists(Histograms::Photon2JetsHists_t& hist,
								int iSuperPhoIndex, int iJet1Index, int iJet2Index, const float fWeight)
{
	TLorentzVector tlPho =  initSpMod->GetSuperPhoton(iSuperPhoIndex)->GetCorVec();
	TLorentzVector tlJet1 =  *(jetMod->GetMyJetNoEMobj_lev6(0,iJet1Index));
	TLorentzVector tlJet2 =  *(jetMod->GetMyJetNoEMobj_lev6(0,iJet2Index));
	
	TLorentzVector tlSum(0,0,0,0);
	tlSum = tlPho + tlJet1 + tlJet2; 
	
	hist.InvMass->Fill(tlSum.M(), fWeight);
}

//____________________________________________________________________________________________________
void EventProperties::FillTwoJetsHists(Histograms::TwoJetsHists_t& hist,
				int iJet1Index, int iJet2Index, const float fWeight)
{
	TLorentzVector tlJet1 =  *(jetMod->GetMyJetNoEMobj_lev6(0,iJet1Index));
	TLorentzVector tlJet2 =  *(jetMod->GetMyJetNoEMobj_lev6(0,iJet2Index));
	
	TLorentzVector tlSum(0,0,0,0);
	tlSum = tlJet1 + tlJet2; 
	float fEtRatio = (tlJet2.E() * TMath::Sin(tlJet2.Theta())) / (tlJet1.E() * TMath::Sin(tlJet1.Theta()));
	
	hist.InvMass->Fill(tlSum.M(), fWeight);
	hist.DelPhi->Fill(fabs(tlJet1.DeltaPhi(tlJet2)), fWeight);
	hist.DelEta->Fill(fabs(tlJet1.Eta() - tlJet2.Eta()), fWeight);
	hist.DelR->Fill(tlJet1.DeltaR(tlJet2), fWeight);
	hist.EtRatio->Fill(fEtRatio, fWeight);
}


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int EventProperties::EndJob() {

	if (GetSummaryStat()) return 0;
		

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[EVP:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[EVP:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[EVP:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	std::cout << sMsg;
	printf("---------------------------------------------------\n");
	return 0;
}
