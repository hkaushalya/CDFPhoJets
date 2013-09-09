#include <iostream>
#include <sstream>
#include "PhoFlat/PhotonTriggerStudy.hh"
#include "PhoFlat/ReadInAll.hh"
#include "TBenchmark.h"

//-------------------------------------------------------------------
PhotonTriggerStudy::PhotonTriggerStudy():
iProgressBy(50000),
sHistFileName("PhotonTriggerStudy.root"),
bMcFlag(true)
{
}

//-------------------------------------------------------------
void PhotonTriggerStudy::Init(TChain *tree)
{
	readIn = new ReadInAll(tree, &stuple);
	readIn->CentralPhotons(1);
	readIn->Init();

	//cannot add folder without creating a dir first!
	histFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = histFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

} // Init

//-------------------------------------------------------------------
void PhotonTriggerStudy::CleanUp()
{
	histFile->Write();
	histFile->Close();
} 


//------ main -------------------------------------------------------
void PhotonTriggerStudy::Main(TChain *ch, int iRunEvents)
{
	
	TBenchmark timer;
	timer.Start("phoana_time");
	
	if (ch) {
		myChain = ch;
	} else {
		std::cout << "NULL chain!" << std::endl;
		exit(1);
	}

	gROOT->Reset();
	Init(myChain);
	BookHistograms();
	
	unsigned int iENTRIES = (unsigned int) myChain->GetEntries();
	unsigned int iEvt2Process = 0; 	//_____________________________________________ events to be processed
	if (iRunEvents < 0 ) {
		iRunEvents = iENTRIES;			//______________________________________________ requested all entries
		iEvt2Process = iENTRIES;
	} else if (iRunEvents > 0 && iRunEvents <= iENTRIES) iEvt2Process = iRunEvents;
	else iEvt2Process = iENTRIES;		//______________________________________________ default

	std::cout << " ==> Entries found :: " << iENTRIES << std::endl;


	unsigned int iEvtProc = 0;			//______________________________________________ number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iFirst400 =0;
	std::cout << "Init done. Ready, Set , GO... " << std::endl; 
	
	timer.Start("looptimer");
	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) {
	  	readIn->GetEntry(iEvtProc);

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << "evt, run,Mc " << stuple.evt_RunNumber << "\t" << stuple.evt_EventNumber << "\t" << stuple.evt_McFlag <<  std::endl;
			std::cout << iEvtProc << "\t events processed.";
			timer.Show("looptimer");
			timer.Start("looptimer");
		}

		  //___________________________________________________________________________ drop the first 400pb-1
		if (stuple.evt_McFlag == 0) {  //______________________________________________ if data
			if (stuple.evt_RunNumber < 190851) {
				iFirst400++;
				continue;
			}
		}
		if (stuple.evt_RunNumber >246231) continue;  //drop the newly added data, 0k dataset in StupleV4
	
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;
	

		if (stuple.vtx_NClass12 < 1) continue;
		if (stuple.vtx_z > 60) continue;

		for (unsigned int i = 0; i < stuple.pho_num; ++i )
		{
			//________________________________________________________________________ this is a hack to find the trigger photon, as I don't have the L3 info to match reconstructed photon to trigger photon.
			if (stuple.tri_pho70 == 1 && stuple.pho_Etc[i] < 70.0) continue;
			if (stuple.tri_pho50 == 1 && stuple.pho_Etc[i] < 50.0) continue;
			if (stuple.tri_pho25iso == 1 && stuple.pho_Etc[i] < 25.0) continue;

			


			if (stuple.pho_TightId[i] == 0)
			{
				if (stuple.tri_pho70 == 1) phoEt_70->Fill(stuple.pho_Etc[i]);
				if (stuple.tri_pho50 == 1) phoEt_50->Fill(stuple.pho_Etc[i]);
				if (stuple.tri_pho25iso == 1) phoEt_iso25->Fill(stuple.pho_Etc[i]);
				cprWgtTight->Fill(stuple.pho_CprWgt[i]);
			} 
			if (stuple.pho_LooseId[i] == 0)
			{
				if (stuple.tri_pho70==1) lphoEt_70->Fill(stuple.pho_Etc[i]);
				if (stuple.tri_pho50==1) lphoEt_50->Fill(stuple.pho_Etc[i]);
				if (stuple.tri_pho25iso==1) lphoEt_iso25->Fill(stuple.pho_Etc[i]);
				cprWgtLoose->Fill(stuple.pho_CprWgt[i]);
			} 
			if (stuple.pho_TightId[i] > 0 && stuple.pho_LooseId[i] == 0)			// >0 means that this has failed tight cuts. must do this as I init them to large negative values and != could mean <0!!
			{
				if (stuple.tri_pho70==1) sphoEt_70->Fill(stuple.pho_Etc[i]);
				if (stuple.tri_pho50==1) sphoEt_50->Fill(stuple.pho_Etc[i]);
				if (stuple.tri_pho25iso==1) sphoEt_iso25->Fill(stuple.pho_Etc[i]);
				cprWgtSideband->Fill(stuple.pho_CprWgt[i]);
			}

			
			if (stuple.pho_PhoenixId[i] != 0) continue;
			if (stuple.pho_TightId[i] == 0)
			{
				if (stuple.tri_pho70 == 1) pphoEt_70->Fill(stuple.pho_Etc[i]);
				else if (stuple.tri_pho50 == 1) pphoEt_50->Fill(stuple.pho_Etc[i]);
				else if (stuple.tri_pho25iso == 1) pphoEt_iso25->Fill(stuple.pho_Etc[i]);
				cprWgtTight_a->Fill(stuple.pho_CprWgt[i]);
			} 
			if (stuple.pho_LooseId[i] == 0)
			{
				if (stuple.tri_pho70==1) plphoEt_70->Fill(stuple.pho_Etc[i]);
				else if (stuple.tri_pho50==1) plphoEt_50->Fill(stuple.pho_Etc[i]);
				else if (stuple.tri_pho25iso==1) plphoEt_iso25->Fill(stuple.pho_Etc[i]);
				cprWgtLoose_a->Fill(stuple.pho_CprWgt[i]);
			} 
			if (stuple.pho_TightId[i] > 0 && stuple.pho_LooseId[i] == 0)			// >0 means that this has failed tight cuts. must do this as I init them to large negative values and != could mean <0!!
			{
				if (stuple.tri_pho70==1) psphoEt_70->Fill(stuple.pho_Etc[i]);
				else if (stuple.tri_pho50==1) psphoEt_50->Fill(stuple.pho_Etc[i]);
				else if (stuple.tri_pho25iso==1) psphoEt_iso25->Fill(stuple.pho_Etc[i]);
				cprWgtSideband_a->Fill(stuple.pho_CprWgt[i]);
			}
	
			
		}	

	}
	
	std::cout << "======== PhotonTriggerStudy =======================" << std::endl;
	std::cout << "Entries found    = " << iENTRIES << std::endl;
	std::cout << "Entries request  = " << iRunEvents << std::endl;
	std::cout << "Events Processed = " << iEvtProc << std::endl;
	std::cout << "Ele./Pho Events  = " << iPhoEleEvts << std::endl;
	std::cout << "Ele Events       = " << iEleEvts << std::endl;
	std::cout << "Pho Events       = " << iPhoEvts << std::endl;
	std::cout << "Other Events     = " << iOther << std::endl;
	std::cout << "First400 Events  = " << iFirst400 << std::endl;
	std::cout << "sum Events       = " << iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400 << std::endl;
	std::cout << "__________________________________________" << std::endl;

	std::cout << "MC Flag setting        = " << bMcFlag << std::endl;
	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	std::cout << "======= END PhotonTriggerStudy ====================" << std::endl;

	CleanUp();
	timer.Show("phoana_time");
} // Main




//-------------------------------------------------------------------
void PhotonTriggerStudy::BookHistograms()
{
	phoEt_iso25 = new TH1F("phoEt_Iso25","E_{T}^{tight #gamma} (pho_iso25 trigger and before phoenix rejection)",300,0,1200);
	phoEt_50 = new TH1F("phoEt_50","E_{T}^{tight #gamma} (pho_50 trigger and before phoenix rejection)",300,0,1200);
	phoEt_70 = new TH1F("phoEt_70","E_{T}^{tight #gamma} (pho_70 trigger and before phoenix rejection)",300,0,1200);

	lphoEt_iso25 = new TH1F("lphoEt_Iso25","E_{T}^{loose #gamma} (pho_iso25 trigger and before phoenix rejection)",300,0,1200);
	lphoEt_50 = new TH1F("lphoEt_50","E_{T}^{loose #gamma} (pho_50 trigger and before phoenix rejection)",300,0,1200);
	lphoEt_70 = new TH1F("lphoEt_70","E_{T}^{loose #gamma} (pho_70 trigger and before phoenix rejection)",300,0,1200);

	sphoEt_iso25 = new TH1F("sphoEt_Iso25","E_{T}^{sideband #gamma} (pho_iso25 trigger and before phoenix rejection)",300,0,1200);
	sphoEt_50 = new TH1F("sphoEt_50","E_{T}^{sideband #gamma} (pho_50 trigger and before phoenix rejection)",300,0,1200);
	sphoEt_70 = new TH1F("sphoEt_70","E_{T}^{sideband #gamma} (pho_70 trigger and before phoenix rejection)",300,0,1200);
	

	pphoEt_iso25 = new TH1F("pphoEt_Iso25","E_{T}^{tight #gamma} (pho_iso25 trigger and after phoenix rejection)",300,0,1200);
	pphoEt_50 = new TH1F("pphoEt_50","E_{T}^{tight #gamma} (pho_50 trigger and after phoenix rejection)",300,0,1200);
	pphoEt_70 = new TH1F("pphoEt_70","E_{T}^{tight #gamma} (pho_70 trigger and after phoenix rejection)",300,0,1200);

	plphoEt_iso25 = new TH1F("plphoEt_Iso25","E_{T}^{loose #gamma} (pho_iso25 trigger and after phoenix rejection)",300,0,1200);
	plphoEt_50 = new TH1F("plphoEt_50","E_{T}^{loose #gamma} (pho_50 trigger and after phoenix rejection)",300,0,1200);
	plphoEt_70 = new TH1F("plphoEt_70","E_{T}^{loose #gamma} (pho_70 trigger and after phoenix rejection)",300,0,1200);

	psphoEt_iso25 = new TH1F("psphoEt_Iso25","E_{T}^{sideband #gamma} (pho_iso25 trigger and after phoenix rejection)",300,0,1200);
	psphoEt_50 = new TH1F("psphoEt_50","E_{T}^{sideband #gamma} (pho_50 trigger and after phoenix rejection)",300,0,1200);
	psphoEt_70 = new TH1F("psphoEt_70","E_{T}^{sideband #gamma} (pho_70 trigger and after phoenix rejection)",300,0,1200);
	

	cprWgtTight = new TH1F("cprWgtTight","CPR weight - #gamma^{Tight}",120,-10,20);
	cprWgtLoose = new TH1F("cprWgtLoose","CPR weight -  #gamma^{Loose}",120,-10,20);
	cprWgtSideband = new TH1F("cprWgtSideband","CPR weight -  #gamma^{Sideband}",120,-10,20);
	
	cprWgtTight_a = new TH1F("cprWgtTight_a","CPR weight - #gamma^{Tight} (after phoenix rejection)",120,-10,20);
	cprWgtLoose_a = new TH1F("cprWgtLoose_a","CPR weight -  #gamma^{Loose} (after phoenix rejection)",120,-10,20);
	cprWgtSideband_a = new TH1F("cprWgtSideband_a","CPR weight -  #gamma^{Sideband} (after phoenix rejection)",120,-10,20);
	
	
	
} //BookHistograms

