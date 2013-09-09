#include "samantha/Pho/Pho2JetsTemp.hh"
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


ClassImp(Pho2JetsTemp)

const UInt_t Pho2JetsTemp::iNPDFs;			//maximum number of variations for any PDF set

//_____________________________________________________________________________
Pho2JetsTemp::Pho2JetsTemp(const char* name, const char* title):
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
  bNoSummary(false),
  bDoPDFSyst(false)
{
	std::cout << "Hello I am Pho2JetsTemp module" << std::endl;
}

//_____________________________________________________________________________
Pho2JetsTemp::~Pho2JetsTemp() {
}

//_____________________________________________________________________________
void Pho2JetsTemp::SaveHistograms() {
}

//_____________________________________________________________________________
void Pho2JetsTemp::BookHistograms()
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
	

	//new_folder = GetHistFolder(this, "TEMP","GEN 6 MC SAMPLE INEFF STUDY");
	//hPhoEt = new TH1F("phoEt","all pho Et",40,0,200);
	//hDetEta = new TH1F("phoEta","all pho Eta",100,-4,4);
	//hLoosePhoEt = new TH1F("loosephoEt","loose pho Et",40,0,200);
	//hLooseDetEta = new TH1F("loosephoEta","loose pho Eta",100,-4,4);
	//hTightPhoEt = new TH1F("tightphoEt","tight pho Et",40,0,200);
	//hHepgPhoEt = new TH1F("HepgphoEt","Hepg pho Et",40,0,200);
	//hHepgPhoEta = new TH1F("HepgphoEta","Hepg pho Eta",100,-4,4);
	
	//new_folder->Add(hPhoEt);
	//new_folder->Add(hLoosePhoEt);
	//new_folder->Add(hDetEta);
	//new_folder->Add(hLooseDetEta);
	//new_folder->Add(hTightPhoEt);
	//new_folder->Add(hHepgPhoEt);
	//new_folder->Add(hHepgPhoEta);

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
//	hLeadJet_Et = new TH1F("TempLeadJetEt","Lead jet Et after removing only the lead photon from jet list",1000,0,500);
//	hPhoLeadJetPtRatio = new TH1F("TempPhoLeadJetPtRatio","Lead jet Pt/Lead Pho Pt - 1 :  after removing only the lead photon from jet list",1000,-2,2);
//	sub_folder->Add(hLeadJet_Et);
//	sub_folder->Add(hPhoLeadJetPtRatio);

	jetmetdelphi_a4 = dynamic_cast<TH1F*> (h1_Jet.JetMetDelPhi->Clone("JetMetDelPhi_a4"));
	sub_folder->Add(jetmetdelphi_a4);
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet,"#gamma + >=1 Jets(15GeV)");

	if (pdfWmod)
	{
		iNaltenativePDFs = pdfWmod->GetNalternativePDFs();
		if (iNaltenativePDFs<1)
		{
			StdOut(__FILE__,__LINE__,3,"Cannot detemined the number of alternative PDFs to generate hitograms. Make sure PDFReweight is created before Pho2JetsTemp. exiting!");
			exit (1);
		}

		std::string sTreeName;
		sTreeName = "PhoJet Events with " + pdfWmod->GetReweightedwithPDF() + " PDF Weights";
		treePDFw = new TTree ("PhoJets",sTreeName.c_str());
		treePDFw->Branch("RunNumber", &aPDFsyst.RunNumber,"RunNumber/I");
		treePDFw->Branch("EventNumber", &aPDFsyst.EventNumber,"EventNumber/I");
		treePDFw->Branch("Npdfs", &aPDFsyst.iNpdfs,"Npdfs/i");
		treePDFw->Branch("PDFWeights", &aPDFsyst.fPdfWeight,"PDFWeights[Npdfs]/I");
		treePDFw->Branch("pj1_Et_pho", &aPDFsyst.pj1_Et_pho,"pj1_Et_pho/F");
		treePDFw->Branch("pj1_Et_j1", &aPDFsyst.pj1_Et_j1,"pj1_Et_j1/F");
		treePDFw->Branch("pj1_InvM_pj1",&aPDFsyst.pj1_InvM_pj1, "pj1_InvM_pj1/F");
		treePDFw->Branch("pj1_Met",&aPDFsyst.pj1_Met, "pj1_Met/F");
		treePDFw->Branch("pj1_Ht",&aPDFsyst.pj1_Ht, "pj1_Ht/F");
		treePDFw->Branch("pj1_Njet15",&aPDFsyst.pj1_Njet15, "pj1_Njet15/I");
		treePDFw->Branch("EtRatLeadJetLeadPho",&aPDFsyst.EtRatLeadJetLeadPho, "EtRatLeadJetLeadPho/F");
		treePDFw->Branch("DelRLeadPhoLeadJet",&aPDFsyst.DelRLeadPhoLeadJet, "DelRLeadPhoLeadJet/F");
		treePDFw->Branch("DelPhiLeadPhoLeadJet",&aPDFsyst.DelPhiLeadPhoLeadJet, "DelPhiLeadPhoLeadJet/F");
		treePDFw->Branch("pj2_Et_pho",&aPDFsyst.pj2_Et_pho, "pj2_Et_pho/F");
		treePDFw->Branch("pj2_Et_j1",&aPDFsyst.pj2_Et_j1, "pj2_Et_j1/F");
		treePDFw->Branch("pj2_Et_j2",&aPDFsyst.pj2_Et_j2, "pj2_Et_j2/F");
		treePDFw->Branch("pj2_InvM_pj1",&aPDFsyst.pj2_InvM_pj1, "pj2_InvM_pj1/F");
		treePDFw->Branch("pj2_InvM_pj2",&aPDFsyst.pj2_InvM_pj2, "pj2_InvM_pj2/F");
		treePDFw->Branch("pj2_InvM_pj1j2",&aPDFsyst.pj2_InvM_pj1j2, "pj2_InvM_pj1j2/F");
		treePDFw->Branch("pj2_InvM_j1j2",&aPDFsyst.pj2_InvM_j1j2, "pj2_InvM_j1j2/F");
		treePDFw->Branch("pj2_Met",&aPDFsyst.pj2_Met, "pj2_Met/F");
		treePDFw->Branch("pj2_Ht",&aPDFsyst.pj2_Ht, "pj2_Ht/F");

		new_folder = GetHistFolder(this, "PJevtsWithPDFw","#gamma + >=1 Jet case");
		new_folder->Add(treePDFw);

		//sub_folder = new_folder->AddFolder("PDFSyst","Sample PDFSystHists, phoEt and jetEt");
		sub_folder = new_folder->AddFolder("PDFSyst","Sample PDFSystHists, phoEt");
		for (int i =0; i < iNaltenativePDFs; ++i)
		{
			std::stringstream pname, ptitle, jname, jtitle;
			pname << "PDFSystPara"<< i << "_PhoEt";
			if (i==0)
			{
				ptitle << "PDF Syst. Parameter (" << i << ") BASE ;E_{T}^{#gamma};Weighted Events";
			} else 
			{
				ptitle << "PDF Syst. Parameter (" << i << ");E_{T}^{#gamma};Weighted Events";
			}

			TH1F* hist = dynamic_cast<TH1F*> (h1_Pho.EtCorr->Clone());
			hist->SetName(pname.str().c_str());
			hist->SetTitle(ptitle.str().c_str());
			hist->Sumw2();
			hist->SetLineColor(kBlue);
			hist->SetMarkerStyle(24);
			hist->SetMarkerColor(kBlue);
			vPDFhist.push_back(hist);
			
			sub_folder->Add(hist);
		}			

		hQvsE = new TH2F("QvsE","#gamma^{E_{T}>30GeV} + >= 1 Jet^{E_{T}>15GeV};Sum E(#gamma, lead jet);Q from initial partons",800,0,800,800,0,800);
		sub_folder->Add(hQvsE);
		hQERat = new TH1F("QErat","#gamma^{E_{T}>30GeV} + >= 1 Jet^{E_{T}>15GeV};Q/Sum E(#gamma,lead jet);Events",500,0,10);
		sub_folder->Add(hQERat);
		hQSumPtRat = new TH1F("QSumPtRat","#gamma^{E_{T}>30GeV} + >= 1 Jet^{E_{T}>15GeV};Q/Sum Pt(#gamma,lead jet);Events",500,0,10);
		sub_folder->Add(hQSumPtRat);
		hQPhoPtRat = new TH1F("QPhoPtRat","#gamma^{E_{T}>30GeV} + >= 1 Jet^{E_{T}>15GeV};Q/P_{T}^{#gamma};Events",500,0,10);
		sub_folder->Add(hQPhoPtRat);
		hQJetPtRat = new TH1F("QJetPtRat","#gamma^{E_{T}>30GeV} + >= 1 Jet^{E_{T}>15GeV};Q/P_{T}^{#jet};Events",500,0,10);
		sub_folder->Add(hQJetPtRat);
	}


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
int Pho2JetsTemp::BeginJob()
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

	pdfWmod = (TPDFReweight*) ((TStnAna*) GetAna()->GetModule("PDFReweight"));
	if (!pdfWmod) {
		StdOut(__FILE__,__LINE__,2,"TPDFReweight module not found.");
		bPermitRun = false;
	}

	
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");
	//RegisterDataBlock("PROD@JetCluModule-had-cone0.4","TStnJetBlock",&fHadJetBlockClu04);
	RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlockClu04);

	BookHistograms();

	// initialize vars
  

	return 0;
}

//_____________________________________________________________________________
int Pho2JetsTemp::BeginRun()
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
int Pho2JetsTemp::Event(int ientry)
{
  	SetPassed(0);

	// this makes sure that all required mods exists
	if (!bPermitRun) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
	}
	

	Cleanup();		//reinit the tree branch values
	vPDFweights.clear();
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
						&& "Pho2JetsTemp::Photons are not sorted!");
			}
			thePhoton = initSpMod->GetSuperPhoton(vPhoIndex.at(0));
			FillPhotonIDHist(vPhoIndex.at(0));		//Fill all the ID information of my photon for debugging
																// this way I can know if my selection is correct.
			hCount.iCounter->Fill(TIGHT_PHO_EVENTS);

			bool bPassJetSelection = false;
			// this is bcos I need to run this module before JetFilter to select the events with a photon.
			// thne do the jet stuff. And then fill these hists with those jets.
			if (DoPDFSyst()) bPassJetSelection = FillPhoPDFSystHists(vPhoIndex.at(0));
			else bPassJetSelection = DoJetSelection(vPhoIndex.at(0));
			
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
SuperPhoton* Pho2JetsTemp::GetPhoton() const
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
void Pho2JetsTemp::FillPhotonIDHist(const int iSpInd)
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
bool Pho2JetsTemp::DoJetSelection(int iSpIndex)
{
	
	bool bPassJetSelection = false;
	//find the leading detector jet after removing the leading photon
	TLorentzVector tlPhoVec = initSpMod->GetSuperPhoton(iSpIndex)->GetCorVec();
	TVector2 tv2MetVec(jetMod->GetMyMetXCorr(0,3), jetMod->GetMyMetYCorr(0,3));
	
	if (jetMod->GetMyNjet(0,3) >= 1)
	{	
		
		// mono-jet case
		//hPhoLeadJetPtRatio->Fill(tlLeadJetVec.Pt()/tlPhoVec.Pt() - 1);
		//std::cout << "========= "; fHeaderBlock->Print();
		h1_Pho.EtCorr->Fill(tlPhoVec.Pt());
		const float dphi_pm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlPhoVec.Phi())));
		h1_Pho.PhoMetDelPhi->Fill(dphi_pm);
			
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
					//std::cout <<  ijet << "\t"<< tlLeadJetVec.Pt() << "\t" <<  dphi_jm_closest << std::endl;
					tlVecClosestJet = tlLeadJetVec;
				}
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
			bPassJetSelection = true;
			meta4->Fill(tv2MetVec.Mod());
			jetmetdelphi_a4->Fill(dphi_jm_closest);
			phometdelphi_a4->Fill(dphi_pm);

			hCount.iCounter->Fill(2);
			evtPropMod->FillPhotonHists(h1_Pho,iSpIndex);
			evtPropMod->FillEventHists(h1_Evt);
			evtPropMod->FillJetHists(h1_Jet,0);		//lead jet
			evtPropMod->FillPhoton1JetHists(h1_PhoJet,iSpIndex,0);

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
			}

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

	bPassJetSelection;
}


/*-------------------------------------------------------------------*/
// Fill photon PDF systematic hists
/*-------------------------------------------------------------------*/
bool Pho2JetsTemp::FillPhoPDFSystHists(const int iPhoIndex)
{

	bool bPassJetSelection = false;
	//need a jet with 15GeV
	if (jetMod->GetMyNjet(0,3) >= 1)
	{
		bPassJetSelection = true;
		int iJetCone4 = 0;         //get cone 0.4 jets
		TLorentzVector tlLeadJetVec = *(jetMod->GetMyJetNoEMobj_lev6(iJetCone4,0));

		vPDFweights = pdfWmod->GetEventPDFWeights();
		assert ((int) vPDFweights.size() == iNaltenativePDFs && "Required number of PDF weights not found! pl. check");

		TLorentzVector tlPhoVec = initSpMod->GetSuperPhoton(iPhoIndex)->GetCorVec();
		h1_Pho.EtCorr->Fill(tlPhoVec.Pt());
		h1_Pho.DetEta->Fill(tlPhoVec.Eta());
		//fHeaderBlock->Print();
		for (int i=0; i<iNaltenativePDFs; ++i)
		{
			vPDFhist.at(i)->Fill(tlPhoVec.Pt(), vPDFweights.at(i));
		}

		float fQ = pdfWmod->GetEventQ();
		float fSumE = (tlPhoVec+tlLeadJetVec).Energy();
		hQvsE->Fill (fQ, fSumE);
		hQERat->Fill (fQ/fSumE);
		hQSumPtRat->Fill (fQ/(tlPhoVec+tlLeadJetVec).Pt());
		hQPhoPtRat->Fill (fQ/tlPhoVec.Pt());
		hQJetPtRat->Fill (fQ/tlLeadJetVec.Pt());


		aPDFsyst.RunNumber = fHeaderBlock->RunNumber();
		aPDFsyst.EventNumber = fHeaderBlock->EventNumber();
		aPDFsyst.iNpdfs = iNaltenativePDFs;
		for (int i=0; i < iNaltenativePDFs; ++i)
		{
			aPDFsyst.fPdfWeight[i]= vPDFweights.at(i);
		}
		aPDFsyst.pj1_Et_pho  = tlPhoVec.Pt();
		aPDFsyst.pj1_Et_j1   = tlLeadJetVec.Pt();
		aPDFsyst.pj1_InvM_pj1= (tlPhoVec + tlLeadJetVec).M();
		aPDFsyst.pj1_Njet15  = jetMod->GetMyNjet(0,3);
		aPDFsyst.pj1_Ht      = jetMod->GetMyHtCorr(0,3);
		aPDFsyst.pj1_Met     = jetMod->GetMyMetCorr(0,3);

		aPDFsyst.EtRatLeadJetLeadPho  = tlLeadJetVec.Pt()/tlPhoVec.Pt();
		aPDFsyst.DelRLeadPhoLeadJet   = tlLeadJetVec.DeltaR(tlPhoVec);
		aPDFsyst.DelPhiLeadPhoLeadJet = tlLeadJetVec.DeltaPhi(tlPhoVec);

		TLorentzVector tlSubLeadJetVec(0,0,0,0);
		if (jetMod->GetMyNjet(0,3)>=2)
		{
			tlSubLeadJetVec = *(jetMod->GetMyJetNoEMobj_lev6(iJetCone4,1));
			assert (tlLeadJetVec.Pt() >= tlSubLeadJetVec.Pt() && "Pho2Jets:: jets are not sorted!");

			aPDFsyst.pj2_Et_pho = tlPhoVec.Pt();
			aPDFsyst.pj2_Et_j1  = tlLeadJetVec.Pt();
			aPDFsyst.pj2_Et_j2  = tlSubLeadJetVec.Pt();
			aPDFsyst.pj2_InvM_pj1   = (tlPhoVec + tlLeadJetVec).M();
			aPDFsyst.pj2_InvM_pj2   = (tlPhoVec + tlSubLeadJetVec).M();
			aPDFsyst.pj2_InvM_pj1j2 = (tlPhoVec + tlLeadJetVec + tlSubLeadJetVec).M();
			aPDFsyst.pj2_InvM_j1j2  = (tlLeadJetVec + tlSubLeadJetVec).M();
			aPDFsyst.pj2_Met        = jetMod->GetMyMetCorr(0,3);
			aPDFsyst.pj2_Ht         = jetMod->GetMyHtCorr(0,3);
		}

		treePDFw->Fill(); //fill the tree only if the event has a pho+>=1jet
	}

	return bPassJetSelection;
	
}

/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void Pho2JetsTemp::Cleanup()
{
	aPDFsyst.RunNumber = -999;
	aPDFsyst.EventNumber = -999;
	aPDFsyst.pj1_Et_pho = -999;
	aPDFsyst.pj1_Et_j1 = -999;
	aPDFsyst.pj1_InvM_pj1 = -999;
	aPDFsyst.pj1_Met = -999;
	aPDFsyst.pj1_Ht = -999;
	aPDFsyst.pj1_Njet15 = -999;
	aPDFsyst.EtRatLeadJetLeadPho = -999;
	aPDFsyst.DelRLeadPhoLeadJet = -999;
	aPDFsyst.DelPhiLeadPhoLeadJet = -999;
	aPDFsyst.pj2_Et_pho = -999;
	aPDFsyst.pj2_Et_j1 = -999;
	aPDFsyst.pj2_Et_j2 = -999;
	aPDFsyst.pj2_InvM_pj1 = -999;
	aPDFsyst.pj2_InvM_pj2 = -999;
	aPDFsyst.pj2_InvM_pj1j2 = -999;
	aPDFsyst.pj2_InvM_j1j2 = -999;
	aPDFsyst.pj2_Met = -999;
	aPDFsyst.pj2_Ht = -999;


} //Cleanup


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void Pho2JetsTemp::SetHaloType(const int iT)
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
void Pho2JetsTemp::SetPhotonType(const int type)
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
void Pho2JetsTemp::FillDataBlocks(int ientry)
{
	fGenpBlock->GetEntry(ientry);
	fJetBlockClu04->GetEntry(ientry);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int Pho2JetsTemp::EndJob() {

	if (GetSummaryStat()) return 0;
		

	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[P2J:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[P2J:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[P2J:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	sMsg += "[P2J:03:] DoPDFSyst        ----------- = " + ToStr((int)DoPDFSyst()) + " (1 =YES, 0= NO)"+ "\n";
	sMsg += "[P2J:04:] Photon selection ----------- = " + ToStr(GetPhotonType()) + " (" + GetPhotonTypeLabel() + ")"+ "\n";
	sMsg += "[P2J:05:] Photon Min Et -------------- = " + ToStr(GetPhoMinEt()) + " GeV" + "\n";
	sMsg += "[P2J:06:] Photon Min<|DetEta|<Max ---- = " + ToStr(GetPhoMinDetEta()) + ", "+ ToStr(GetPhoMaxDetEta()) + "\n";
	if (GetRejPhxPhos())
	sMsg += "[P2J:07:] Phoenix Photons ------------ = " + ToStr(GetRejPhxPhos()) + " (Rejected)" + "\n";
	else 
	sMsg += "[P2J:07:] Phoenix Photons ------------ = " + ToStr(GetRejPhxPhos()) + " (NOT Rejected)" + "\n";
	sMsg += "[P2J:08:] Halo Rejection Type -------- = " + ToStr(GetHaloType()) + "\n";
	if (GetHeaderBlock()->McFlag())
	sMsg += "[P2J:09:] EM Time window used -------- = NO CUT (MC)\n";
	else
	sMsg += "[P2J:09:] EM Time window used -------- = " + ToStr(GetEmTimeWindowMin()) + "," + ToStr(GetEmTimeWindowMax())+ "\n";
	sMsg += "[P2J:10:] Incl 1 Tight Photon Events - = " + ToStr(hCount.iCounter->GetBinContent(TIGHT_PHO_EVENTS+1)) + "\n";
	sMsg += "[P2J:11:] 1 Tight Photon only Events - = " + ToStr(hCount.iCounter->GetBinContent(ONE_TIGHT_PHO_ONLY+1)) + "\n";
	sMsg += "[P2J:12:] 2 Tight Photons Only Events  = " + ToStr(hCount.iCounter->GetBinContent(TWO_TIGHT_PHO_ONLY+1)) + "\n";
	sMsg += "[P2J:13:] 1 T+>0 Loose Pho Only Events = " + ToStr(hCount.iCounter->GetBinContent(ONE_TIGHT_ANY_LOOSE_PHO+1)) + "\n";
	sMsg += "---------------------------------------------------\n";

	std::cout << sMsg;

	return 0;
}



void Pho2JetsTemp::SetPhoMaxDetEta(const float eta)
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

void Pho2JetsTemp::SetPhoMinDetEta(const float eta)
{
	if (fabs(eta) > 1.1)
	{
		StdOut(__FILE__,__LINE__,3,"Photon Detector eta require is > than what is returned in TStnPhoton::Detector() == 0??");
		exit (1);
	} else {
 		fPhoMinDetEta = fabs(eta);
	}
}

void Pho2JetsTemp::SetPhoMinEt(const float et)
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
