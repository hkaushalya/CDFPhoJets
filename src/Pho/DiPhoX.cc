/////////////////////////////////////////////////////////////////////
// This is designed to emulate Sasha' gg+X search. I needed to     //
// the curretn version of MEt Model I have (JetFilterV2) is not    //
// broken. I reproduced the MEtSig from this with current settings //
// I used in JetFilterV2. Here I choose two central photons with   //
// Et>13GeV. I do not remove tri-photon events and do not swap     //
// vertices to minimize the MET. I remove all the EM objects       //
// identified from the jet list. Sasha removed only th two tight   //
// photons. In this very first test MEtSig is better than for my   //
// pho+jets sample. But was not good enough according to Sasha.    //
// See elog: 1206                                                  //
/////////////////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/*
 * $Log: DiPhoX.cc,v $
 * Revision 1.3  2009/06/26 19:57:58  samantha
 * ADDED: Option to manually set the MC flag. This way I can save some CPU time
 * by not running the halo/PMT tagging module when I run over MC samples. I also
 * checks this manual setting at run time to make sure there is not conflict
 * of setting the MC flag manually.
 *
 * Revision 1.2  2009/06/09 21:49:31  samantha
 * MODIFIED: 1. I somehow had dropped the TightPhoton requirement from the photon
 * selection! All the tests I did was wrong. Now it is fixed. I also fixed the
 * event counting to see how many tri-pho events are in the sample. Now the beam
 * halo and PMT spike check will be done only for MC. MC does not simulate either
 * of the cases. Though we could use MC to figure out the misId rate of halo cuts.
 * But using electrons is better.
 *
 * Revision 1.1  2009/06/03 04:02:20  samantha
 * My own version of gg+X to check if my version of MetModel (JetFilterV2).
 * Require two photons (central and Et>13GeV). See elog:1206 for the first
 * result. For Sasha MetSig is not that good. But the difference is smaller
 * than for my pho+jets.
 *
 *
 */

#include "samantha/Pho/DiPhoX.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>
#include <iomanip>

ClassImp(DiPhoX)

//_____________________________________________________________________________
DiPhoX::DiPhoX(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bPermitRun(true),
  iPhotonType(0),
  fPhoMaxDetEta(1.1),
  fPhoMinDetEta(0),
  fPhoMinEt(30),
  bRejPhxPhos(true),
  bNoSummary(false),
  bManualMCflag(false)
{
	std::cout << "Hello I am DiPhoX module" << std::endl;
}

//_____________________________________________________________________________
DiPhoX::~DiPhoX() {
}

//_____________________________________________________________________________
void DiPhoX::SaveHistograms() {
}

//_____________________________________________________________________________
void DiPhoX::BookHistograms()
{
	DeleteHistograms();

	TFolder *new_folder;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");				// 0
	vCounterHistLabels.push_back("Events Passed");					// 1
	vCounterHistLabels.push_back("All Tight #gamma Events");		// 2
	vCounterHistLabels.push_back("One tight #gamma only Events");	// 3
	vCounterHistLabels.push_back("Two tight #gamma only events");	// 4
	vCounterHistLabels.push_back("Two or more tight #gamma events");	// 5
	HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	
	new_folder = GetHistFolder(this, "DiPho","DiPhoX hists");

	Hist.LeadPhoEt = new TH1F("LeadPhoEt","DiPhoX: Lead #gamma E_{T}",250,0,500);
	Hist.SubPhoEt = new TH1F("SubLeadPhoEt","DiPhoX: Subleading #gamma E_{T}",250,0,500);
	Hist.LeadPhoEta = new TH1F("LeadPhoEta","DiPhoX: Lead #gamma #Eta",320,4,-4);
	Hist.SubPhoEta = new TH1F("SubLeadPhoEta","DiPhoX: Subleading #gamma #Eta",320,4,-4);
	Hist.DiPhoInvMass = new TH1F("DiPhoInvMass","DiPhoX: DiPhoton Invariant Mass",500,0,1000);
	Hist.DelPhi = new TH1F("DiPhoDelPhi","DiPhoX: #Delta#Phi^{#gammma 1,#gamma 2}",320,0,3.2);
	new_folder->Add(Hist.LeadPhoEt);
	new_folder->Add(Hist.SubPhoEt);
	new_folder->Add(Hist.LeadPhoEta);
	new_folder->Add(Hist.SubPhoEta);
	new_folder->Add(Hist.DiPhoInvMass);
	new_folder->Add(Hist.DelPhi);
}

//_____________________________________________________________________________
int DiPhoX::BeginJob()
{
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bPermitRun = false;
	}

  	tightMod = (TagTightPhotons*) ((TStnAna*) GetAna()->GetModule("TagTightPhotons"));
	if (!tightMod) {
		StdOut(__FILE__,__LINE__,2,"Tight photon tag module not found.");
		bPermitRun = false;
	}

	if (!ManualMCflag())
	{
		pmtMod = (TagPMTSpikes*) ((TStnAna*) GetAna()->GetModule("TagPMTSpikes"));
		if (!pmtMod) {
			StdOut(__FILE__,__LINE__,2,"PMT Spike module not found.");
			bPermitRun = false;
		}

		haloMod = (TagBeamHalo*) ((TStnAna*) GetAna()->GetModule("TagBeamHalo"));
		if (!haloMod) {
			StdOut(__FILE__,__LINE__,2,"Halo tagging module not found.");
			bPermitRun = false;
		}
	}

  	eleMod = (TagElectrons*) ((TStnAna*) GetAna()->GetModule("TagElectrons"));
	if (!eleMod) {
		StdOut(__FILE__,__LINE__,2,"Electron tagging module not found.");
		bPermitRun = false;
	}

	phoenixMod = (TagPhoenixPhotons*) ((TStnAna*) GetAna()->GetModule("TagPhoenixPhotons"));
	if (!phoenixMod) {
		StdOut(__FILE__,__LINE__,2,"PhoenixPhotons Tagging module not found.");
		bPermitRun = false;
	}

				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();


	return 0;
}

//_____________________________________________________________________________
int DiPhoX::BeginRun()
{

	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
  	// checks if there is a conflict in settings
	if (qMc && !ManualMCflag())
	{
		StdOut(__FILE__,__LINE__,3,"Conflicting settings for MC flag! please check");
		exit (1);
	}
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int DiPhoX::Event(int ientry)
{
  	SetPassed(0);

	// this makes sure that all required mods exists
	if (!bPermitRun) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
		return 0;
	}

	hCount.iCounter->Fill(EVENTS_PROCESSED);

	std::vector<int> vPhoIndex;
	bool bSkipThisEvent = false;

  	int NsuperPho = initSpMod->GetSuperPhoSize();
	for (int i=0; i < NsuperPho; i++) {
		if (! (initSpMod->GetSuperPhoton(i)->IsTightPhoton())) continue;

		//Det Eta cut
		if (fabs(initSpMod->GetSuperPhoton(i)->GetDetEta()) > GetPhoMaxDetEta() ||
			fabs(initSpMod->GetSuperPhoton(i)->GetDetEta()) < GetPhoMinDetEta()) continue;


		// I should check Et>30 after BH cut. beam halo photons below 30 GeV
		// can be the extra photons. Then I need to throw that out. --12-11-2208
		if (initSpMod->GetSuperPhoton(i)->GetEtCorr() < GetPhoMinEt()) continue;

		if (GetRejPhxPhos())
		{
			if (initSpMod->GetSuperPhoton(i)->IsPhoenix()) continue;		// reject phoenix photons
		}

		// Beamhalo check
		if (! (fHeaderBlock->McFlag()))
		{
			int iHaloId = initSpMod->GetSuperPhoton(i)->GetBeamHaloId();
			//remove halo for a given type. must be same type as used for 
			// halo template!!
			if ( ((iHaloId >> GetHaloType()) & 0x1) ) {
				bSkipThisEvent = true;
				break;
			}
		}

		vPhoIndex.push_back(i);
	}


	if (!bSkipThisEvent) {
		//di-photon events

		if (vPhoIndex.size()>0) hCount.iCounter->Fill(TIGHT_PHO_EVENTS);
		if (vPhoIndex.size()==1) hCount.iCounter->Fill(ONE_TIGHT_PHO_ONLY);
		if (vPhoIndex.size()==2) hCount.iCounter->Fill(TWO_TIGHT_PHO_ONLY);
		if (vPhoIndex.size()>2) hCount.iCounter->Fill(MORETHAN2_TIGHT_PHO);
		
		if (vPhoIndex.size() >= 2) {
			FillPhotonHists(vPhoIndex.at(0), vPhoIndex.at(1));
			SetPassed(1);
		}
	}

	if (GetPassed()) {
		hCount.iCounter->Fill(EVENTS_PASSED);
	}
	
	return 0;

} // Event



/*-------------------------------------------------------------------*/
void DiPhoX::FillPhotonHists(const int iPho1, const int iPho2)
{
	if ( (iPho1 < 0 || iPho1 > initSpMod->GetSuperPhoSize())
		|| (iPho2 < 0 || iPho2 > initSpMod->GetSuperPhoSize()))
	{
		StdOut(__FILE__,__LINE__,3,"Invalid SuperPhoton index/indices!");
		return;
	}
	
	float fPho1Et = initSpMod->GetSuperPhoton(iPho1)->GetEtCorr();
	float fPho2Et = initSpMod->GetSuperPhoton(iPho2)->GetEtCorr();

	TLorentzVector leadPhoVec, subPhoVec;
	if (fPho1Et > fPho2Et)
	{
		leadPhoVec = initSpMod->GetSuperPhoton(iPho1)->GetCorVec();
		subPhoVec = initSpMod->GetSuperPhoton(iPho2)->GetCorVec();
	} else {
		leadPhoVec = initSpMod->GetSuperPhoton(iPho2)->GetCorVec();
		subPhoVec = initSpMod->GetSuperPhoton(iPho1)->GetCorVec();
	}


	Hist.LeadPhoEt->Fill(leadPhoVec.Pt());
	Hist.SubPhoEt->Fill(subPhoVec.Pt());
	Hist.LeadPhoEta->Fill(leadPhoVec.Eta());
	Hist.SubPhoEta->Fill(subPhoVec.Eta());
	Hist.DiPhoInvMass->Fill((leadPhoVec+subPhoVec).M());
	Hist.DelPhi->Fill (fabs(leadPhoVec.DeltaPhi(subPhoVec)));

}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void DiPhoX::SetHaloType(const int iT)
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
void DiPhoX::SetPhotonType(const int type)
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
void DiPhoX::FillDataBlocks(int ientry)
{
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int DiPhoX::EndJob() {

	if (GetSummaryStat()) return 0;
		
	std::string sModName(GetName());
	std::string sMsg;
	sMsg  = "[GGX:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[GGX:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[GGX:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	sMsg += "[GGX:03:] Photon selection ----------- = " + ToStr(GetPhotonType()) + " (" + GetPhotonTypeLabel() + ")"+ "\n";
	sMsg += "[GGX:04:] Photon Min Et -------------- = " + ToStr(GetPhoMinEt()) + " GeV" + "\n";
	sMsg += "[GGX:05:] Photon Min<|DetEta|<Max ---- = " + ToStr(GetPhoMinDetEta()) + ", "+ ToStr(GetPhoMaxDetEta()) + "\n";
	if (GetRejPhxPhos())
	sMsg += "[GGX:06:] Phoenix Photons ------------ = " + ToStr(GetRejPhxPhos()) + " (Rejected)" + "\n";
	else 
	sMsg += "[GGX:06:] Phoenix Photons ------------ = " + ToStr(GetRejPhxPhos()) + " (NOT Rejected)" + "\n";
	sMsg += "[GGX:07:] Halo Rejection Type -------- = " + ToStr(GetHaloType()) + "\n";
	sMsg += "[GGX:08:] Incl 1 Tight Photon Events - = " + ToStr(hCount.iCounter->GetBinContent(TIGHT_PHO_EVENTS+1)) + "\n";
	sMsg += "[GGX:09:] 1 Tight Photon only Events - = " + ToStr(hCount.iCounter->GetBinContent(ONE_TIGHT_PHO_ONLY+1)) + "\n";
	sMsg += "[GGX:10:] 2 Tight Photons Only Events  = " + ToStr(hCount.iCounter->GetBinContent(TWO_TIGHT_PHO_ONLY+1)) + "\n";
	sMsg += "[GGX:11:] >2 Tight Photons Events      = " + ToStr(hCount.iCounter->GetBinContent(MORETHAN2_TIGHT_PHO+1)) + "\n";
	sMsg += "---------------------------------------------------\n";

	std::cout << sMsg;
	return 0;
}
