
/*  $Id: PhotonJets2.cc,v 1.1 2013/02/28 03:37:30 samantha Exp $
 *  $Log: PhotonJets2.cc,v $
 *  Revision 1.1  2013/02/28 03:37:30  samantha
 *  Final commit. no checks.! these were never commited to cvs before!
 *
 */


#include <iostream>
#include <sstream>
#include "PhoFlat/PhotonJets2.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
#include "TBenchmark.h"
#include "../RootTools/CommonTools.hh"

//-------------------------------------------------------------------
PhotonJets2::PhotonJets2():
iProgressBy(50000),
sHistFileName("PhotonJets2.root"),
bMcFlag(true)
{
}

//-------------------------------------------------------------
void PhotonJets2::Init(TChain *tree)
{
	readIn = new ReadInAll(tree, &stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->Met(1);

	if (bMcFlag) {
		readIn->EmUpPhotons(1);
		readIn->EmUpElectrons(1);
		readIn->EmDownPhotons(1);
		readIn->EmDownElectrons(1);
		readIn->JesUpJets(1);
		readIn->JesDownJets(1);
		readIn->GenLevel(0);
		
		for (int i=0; i<3; i++) {
			vMcCount.push_back(0);
			vMcUpCount.push_back(0);
			vMcDownCount.push_back(0);
		}
	} else {
		for (int i=0; i<3; i++) {
			vHaloCount.push_back(0);
			vCosmicCount.push_back(0);
			vQCDCount.push_back(0);
			vSignalCount.push_back(0);
		}

	}

	readIn->Init();

	//cannot add folder without creating a dir first!
	histFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = histFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

	//get the sideband systematic function
	//TFile f("../RootTools/DiJet_Sideband_FitFuncForPhoEt.root");
	//tf1_sideband = (TF1*) f.Get("tf1_sideband");
	
	//for Ht plot reweighing of photon Et
	TFile f("PhoEt_SigSideRatio_FitFunc.root");
	tf1_sideband = (TF1*) f.Get("tf1_Pho_signal_side_ratio");

	std::stringstream sEr;
	sEr << __FILE__ << ":" << __FUNCTION__ << ":" << __LINE__ << ":" 
			<< "weight function for sideband systematic is not found!";
	assert (tf1_sideband != NULL && sEr.str().c_str());
	
} // Init

//-------------------------------------------------------------------
void PhotonJets2::CleanUp()
{
	vMcCount.clear();
	vMcUpCount.clear();
	vMcDownCount.clear();
	vHaloCount.clear();
	vCosmicCount.clear();
	vSignalCount.clear();
	
	histFile->Write();
	histFile->Close();
} 


//------ main -------------------------------------------------------
void PhotonJets2::Main(TChain *ch, int iRunEvents)
{
	
	TBenchmark timer;
	timer.Start("phoana_time");
	
	if (ch) {
		myChain = ch;
		//ch->Print();	
	} else {
		std::cout << "NULL chain!" << std::endl;
		assert(false);
	}
	if (ch->GetEntries()<1)
	{
		std::cout << "No entries found!" << std::endl;
		return;
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

	//do a check for the MC. this is a common mistake I make
	
	  readIn->GetEntry(1);
	  if (stuple.evt_McFlag != GetMcFlag()) 
	  {
		  std::cout << red << "MC Flag setting does not match what is in Stuple. pleas check. returning!" << clearatt << std::endl;
		  //std::cout << "Changing the MCFlag to Stuple setting McFlag = " << stuple.evt_McFlag << clearatt << std::endl;
		  //SetMcFlag(stuple.evt_McFlag); //for some reason this seems to overwrite counting arrays? 12-07-2009
		  //return; // simple return causes trouble when rerun. so drastic measure needed, even exit(1) would not work
		  assert (false);
	  }

	unsigned int iEvtProc = 0;			//______________________________________________ number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iFirst400 =0;
	std::cout << "Init done. Ready, Set , GO... " << std::endl; 
	
	timer.Start("looptimer");
	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) {
	  	readIn->GetEntry(iEvtProc);
		//if (stuple.evt_RunNumber != 196983 || stuple.evt_EventNumber !=1533562) continue;
		if (stuple.evt_RunNumber != 196989 || stuple.evt_EventNumber != 8696978) continue;
		std::cout << cyan << "****************************************************************************" << clearatt << std::endl;
		std::cout << "evt, run,Mc " << stuple.evt_RunNumber << "\t" << stuple.evt_EventNumber << "\t" << stuple.evt_McFlag <<  std::endl;
		//stuple.DumpJetBlock();
		//stuple.DumpPhotonBlock();
		//stuple.DumpElectronBlock();
		stuple.Dump(1);

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << iEvtProc << "\t events processed [" << (int)( iEvtProc/(double)iEvt2Process * 100)<< "%] ";
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
	
		//std::cout << "pho_num=" << stuple.pho_num << std::endl;
		//std::cout << "pho_num status=" << myChain->GetBranchStatus("pho_num") << std::endl;
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;

		CommonVars commVars;
		PhotonList centralPhos, emupPhos, emdownPhos;
		ElectronList centralEles, emupEles, emdownEles;	//do not comment these out. I need these to recreate the jet list with all unused EM objects
		JetList centralJets, jesupJets, jesdownJets;

		GetCommonVariables(stuple, commVars);
		CreatePhotonLists(stuple, centralPhos,0);
		CreateElectronLists(stuple, centralEles,0);
		CreateJetLists(stuple, centralJets,0);
		//centralJets.Dump();
		if (bMcFlag)		//__________________________ need EM/JES UP/DOWN collections
		{
			CreatePhotonLists(stuple, emupPhos,1);
			CreatePhotonLists(stuple, emdownPhos,-1);
			CreateElectronLists(stuple, emupEles,1);
			CreateElectronLists(stuple, emdownEles,-1);
			CreateJetLists(stuple, jesupJets,1);
			CreateJetLists(stuple, jesdownJets,-1);
		}

		DoMyStuff(commVars, centralPhos, emupPhos, emdownPhos, centralEles, emupEles,
						emdownEles, centralJets, jesupJets, jesdownJets);
		//DoZJetsStuff(commVars, centralPhos, emupPhos, emdownPhos, centralEles, emupEles,
						//emdownEles, centralJets, jesupJets, jesdownJets);


	}
	
	std::cout << "======== PhotonJets2 =======================" << std::endl;
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
	if (bMcFlag) {
		std::cout << "-------------- CENTRAL --------------------- "<< std::endl;
		std::cout << "Pho Evts               = " << vMcCount[0] << std::endl;
		std::cout << "Pho+>=1Jet Evts        = " << vMcCount[1] << std::endl; 
		std::cout << "Pho+>=2Jet Evts        = " << vMcCount[2] << std::endl; 
		std::cout << "-------------- EM/JES UP ------------------- "<< std::endl;
		std::cout << "Pho Evts               = " << vMcUpCount[0] << std::endl;
		std::cout << "Pho+>=1Jet Evts        = " << vMcUpCount[1] << std::endl; 
		std::cout << "Pho+>=2Jet Evts        = " << vMcUpCount[2] << std::endl; 
		std::cout << "-------------- EM/JES DOWN ----------------- "<< std::endl;
		std::cout << "Pho Evts               = " << vMcDownCount[0] << std::endl;
		std::cout << "Pho+>=1Jet Evts        = " << vMcDownCount[1] << std::endl; 
		std::cout << "Pho+>=2Jet Evts        = " << vMcDownCount[2] << std::endl; 
	} else {
		std::cout << "-------------- BEAM HALO ------------------- "<< std::endl;
		std::cout << "Halo Evts              = " << vHaloCount[0] << std::endl;
		std::cout << "Halo+>=1Jet Evts       = " << vHaloCount[1] << std::endl; 
		std::cout << "Halo+>=2Jet Evts       = " << vHaloCount[2] << std::endl; 
		std::cout << "-------------- COSMIC ---------------------- "<< std::endl;
		std::cout << "Cosmic Evts            = " << vCosmicCount[0] << std::endl;
		std::cout << "Cosmic+>=1Jet Evts     = " << vCosmicCount[1] << std::endl; 
		std::cout << "Cosmic+>=2Jet Evts     = " << vCosmicCount[2] << std::endl; 
		std::cout << "-------------- QCD ------------------------- "<< std::endl;
		std::cout << "QCD Evts               = " << vQCDCount[0] << std::endl;
		std::cout << "QCD+>=1Jet Evts        = " << vQCDCount[1] << std::endl; 
		std::cout << "QCD+>=2Jet Evts        = " << vQCDCount[2] << std::endl; 
		std::cout << "-------------- SIGNAL----------------------- "<< std::endl;
		std::cout << "Pho Evts               = " << vSignalCount[0] << std::endl;
		std::cout << "Pho+>=1Jet Evts        = " << vSignalCount[1] << std::endl; 
		std::cout << "Pho+>=2Jet Evts        = " << vSignalCount[2] << std::endl; 
	}
	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	std::cout << "======= END PhotonJets2 ====================" << std::endl;

/*	std::cout << "Events in mutiple bins::" << std::endl;
	std::cout << "Run" << "\t" << "Evt" << "\t" << "Halo" << "\t" << "Cosmic" << "\t" << "Signal" << "\t" << "QCD" << std::endl; 
	for (int i=0; i < EVT.size(); i++)
	{
		if (EVT[i].Halo + EVT[i].Cosmic + EVT[i].Signal + EVT[i].QCD >= 2)
			std::cout << EVT[i].Run << "\t" << EVT[i].Evt << "\t" << EVT[i].Halo << "\t" << EVT[i].Cosmic << "\t" << EVT[i].Signal << "\t" << EVT[i].QCD << std::endl; 
	}
*/
	EVT.clear();
	
	CleanUp();
	timer.Show("phoana_time");
} // Main


//-------------------------------------------------------------------
void PhotonJets2::BookHistograms()
{
	//for halo estimate
	haloPhiWedge_1j = new TH1F("haloPhiWedge_1j","Halo(5)+>=1 Jet: #phi^{#gamma}",24,0,24);
	haloEmTime_1j = new TH1F("haloEmTime_1j","Halo(5)+>=1 Jet: #gamma EmTime",30,-10,10);
	haloPhiWedge_2j = new TH1F("haloPhiWedge_2j","Halo(5)+>=2 Jet: #phi^{#gamma}",24,0,24);
	haloEmTime_2j = new TH1F("haloEmTime_2j","Halo(5)+>=2 Jet: #gamma EmTime",30,-10,10);
	

	phoEt_iso25 = new TH1F("phoEt_Iso25","E_{T}^{#gamma} (pho_iso25 trigger)",300,0,1200);
	phoEt_50 = new TH1F("phoEt_50","E_{T}^{#gamma} (pho_50 trigger)",300,0,1200);
	phoEt_70 = new TH1F("phoEt_70","E_{T}^{#gamma} (pho_70 trigger)",300,0,1200);

	Hsignal.z1j_zpt	= new TH1F("z1j_zpt","#gamma+>=1Jet :- p_{T}^{Z} ",200,-500,500);
	Hsignal.z1j_zet	= new TH1F("z1j_zet","#gamma+>=1Jet :- E_{T}^{Z} ",250,0,500);
	Hsignal.z1j_zmass	= new TH1F("z1j_zmass","#gamma+>=1Jet :- Z mass ",200,0,1000);
	Hsignal.z1j_zeta	= new TH1F("z1j_zeta","#gamma+>=1Jet :- Event #Eta^{Z} ",50,-5,5);
	Hsignal.z1j_zj_mass	= new TH1F("z1j_zj1mass","#gamma+>=1Jet :- Invariant Mass (Z, Lead Jet)",200,0,1000);
	
	Hsignal.z2j_zpt	= new TH1F("z2j_zpt","#gamma+>=2Jets :- p_{T}^{Z} ",200,-500,500);
	Hsignal.z2j_zet	= new TH1F("z2j_zet","#gamma+>=2Jets :- E_{T}^{Z} ",250,0,500);
	Hsignal.z2j_zmass	= new TH1F("z2j_zmass","#gamma+>=2Jets :- Z mass ",200,0,1000);
	Hsignal.z2j_zeta	= new TH1F("z2j_zeta","#gamma+>=2Jets :- Event #Eta^{Z} ",50,-5,5);
	Hsignal.z2j_z1j_mass	= new TH1F("z2j_zj1mass","#gamma+>=2Jets :- Invariant Mass (Z, Lead Jet)",200,0,1000);
	Hsignal.z2j_z2j_mass	= new TH1F("z2j_zj2mass","#gamma+>=2Jets :- Invariant Mass (Z, Second Lead Jet)",200,0,1000);
	Hsignal.z2j_zj1j2_mass	= new TH1F("z2j_zj1j2mass","#gamma+>=2Jets :- Invariant Mass (Z, Two Lead Jets)",200,0,1000);

	
	Histograms HistoMan;
	std::string text2title1, text2title2;

	TDirectory *l2dir, *l3dir, *l4dir;

	if (bMcFlag) {			//plots for MC

		//CENTRAL
		text2title1 = "CENTRAL #gamma_{MC}+>=1Jet :";
		text2title2 = "CENTRAL #gamma_{MC}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("CENTRAL","Central values");
		l2dir->cd();

		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hmc.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hmc.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hmc.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hmc.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hmc.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hmc.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hmc.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hmc.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hmc.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hmc.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hmc.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hmc.p2j_TwoJets,text2title2.c_str());



		//EM AND JES UP
		text2title1 = "EM/JES UP #gamma_{MC}+>=1Jet :";
		text2title2 = "EM/JES UP #gamma_{MC}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("EMJESUP","EM and JES up values");
		l2dir->cd();

		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HmcUp.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HmcUp.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcUp.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(HmcUp.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HmcUp.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HmcUp.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcUp.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcUp.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HmcUp.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HmcUp.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(HmcUp.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(HmcUp.p2j_TwoJets,text2title2.c_str());



		//EM AND JES DOWN
		text2title1 = "EM/JES DOWN #gamma_{MC}+>=1Jet :";
		text2title2 = "EM/JES DOWN #gamma_{MC}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("EMJESDOWN","EM and JES down values");
		l2dir->cd();

		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HmcDown.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HmcDown.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcDown.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(HmcDown.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HmcDown.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HmcDown.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcDown.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcDown.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HmcDown.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HmcDown.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(HmcDown.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(HmcDown.p2j_TwoJets,text2title2.c_str());


	} else {

			//plots for signal
		text2title1 = "#gamma_{signal}+>=1Jet :";
		text2title2 = "#gamma_{signal}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("SIGNAL","Signal");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsignal.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsignal.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsignal.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hsignal.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsignal.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsignal.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsignal.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsignal.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsignal.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsignal.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hsignal.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hsignal.p2j_TwoJets,text2title2.c_str());



			//plots for halo
		text2title1 = "#gamma_{halo}+>=1Jet :";
		text2title2 = "#gamma_{halo}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("HALO","Beam Halo");
		l2dir->cd();
	
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hhalo.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hhalo.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hhalo.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hhalo.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hhalo.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hhalo.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hhalo.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hhalo.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hhalo.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hhalo.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hhalo.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hhalo.p2j_TwoJets,text2title2.c_str());



			//plots for cosmic
		text2title1 = "#gamma_{cosmic}+>=1Jet :";
		text2title2 = "#gamma_{cosmic}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("COSMIC","Cosmic");
		l2dir->cd();

		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hcosmic.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hcosmic.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hcosmic.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hcosmic.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hcosmic.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hcosmic.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hcosmic.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hcosmic.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hcosmic.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hcosmic.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hcosmic.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hcosmic.p2j_TwoJets,text2title2.c_str());



			//qcd sideband
		text2title1 = "#gamma_{sideband}+>=1Jet :";
		text2title2 = "#gamma_{sideband}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("SIDEBAND","QCD sideband");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hqcd.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hqcd.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hqcd.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hqcd.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hqcd.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hqcd.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hqcd.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hqcd.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hqcd.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hqcd.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hqcd.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hqcd.p2j_TwoJets,text2title2.c_str());

		
	}
	
} //BookHistograms



void PhotonJets2::CreatePhotonLists(const Stuple stuple, PhotonList& pholist, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
			pholist.pho_num 	 =  stuple.pho_num;
			pholist.pho_Ntight = stuple.pho_Ntight;
			pholist.pho_Nloose = stuple.pho_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_E[i]);
				pholist.pho_Px.push_back(stuple.pho_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_matchJetIndex[i]);
				
				pholist.pho_CprWgt.push_back(stuple.pho_CprWgt[i]);
				pholist.pho_CprSys1.push_back(stuple.pho_CprSys1[i]);
				pholist.pho_CprSys2.push_back(stuple.pho_CprSys2[i]);
				pholist.pho_CprSys3.push_back(stuple.pho_CprSys3[i]);
				pholist.pho_CprSys4.push_back(stuple.pho_CprSys4[i]);
				pholist.pho_CprSys5.push_back(stuple.pho_CprSys5[i]);
				pholist.pho_CprSys6.push_back(stuple.pho_CprSys6[i]);
				pholist.pho_CprSys7.push_back(stuple.pho_CprSys7[i]);
				pholist.pho_CprSys8.push_back(stuple.pho_CprSys8[i]);
				
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		case 1:
		{
			pholist.pho_num 	 =  stuple.pho_up_num;
			pholist.pho_Ntight = stuple.pho_up_Ntight;
			pholist.pho_Nloose = stuple.pho_up_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_up_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_up_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_up_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_up_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_up_E[i]);
				pholist.pho_Px.push_back(stuple.pho_up_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_up_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_up_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_up_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_up_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_up_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_up_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_up_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_up_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_up_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_up_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_up_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_up_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_up_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_up_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_up_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_up_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_up_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_up_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_up_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_up_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_up_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_up_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_up_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_up_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_up_matchJetIndex[i]);
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		case -1:
		{
			pholist.pho_num 	 =  stuple.pho_down_num;
			pholist.pho_Ntight = stuple.pho_down_Ntight;
			pholist.pho_Nloose = stuple.pho_down_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_down_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_down_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_down_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_down_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_down_E[i]);
				pholist.pho_Px.push_back(stuple.pho_down_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_down_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_down_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_down_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_down_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_down_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_down_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_down_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_down_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_down_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_down_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_down_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_down_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_down_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_down_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_down_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_down_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_down_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_down_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_down_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_down_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_down_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_down_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_down_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_down_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_down_matchJetIndex[i]);
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid photon collection requested! can make central/em up/em down photon lists only!");
			assert(false);

		}
		
	}	//switch

/*}}}*/
}	//CreatePhotonLists



void PhotonJets2::CreateElectronLists(const Stuple stuple, ElectronList& list, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
			list.ele_num    = stuple.ele_num;
			list.ele_Ntight = stuple.ele_Ntight;
			list.ele_Nloose = stuple.ele_Nloose;

			for (unsigned int i=0; i < stuple.ele_num; i++)
			{
				list.ele_Index.push_back(stuple.ele_Index[i]);
				list.ele_PhoBlockIndex.push_back(stuple.ele_PhoBlockIndex[i]);
				list.ele_EleBlockIndex.push_back(stuple.ele_EleBlockIndex[i]);
				list.ele_Etc.push_back(stuple.ele_Etc[i]);
				list.ele_E.push_back(stuple.ele_E[i]);
				list.ele_Px.push_back(stuple.ele_Px[i]);
				list.ele_Py.push_back(stuple.ele_Py[i]);
				list.ele_Pz.push_back(stuple.ele_Pz[i]);
				list.ele_Detector.push_back(stuple.ele_Detector[i]);
				list.ele_DetEta.push_back(stuple.ele_DetEta[i]);
				list.ele_DetPhi.push_back(stuple.ele_DetPhi[i]);
				list.ele_XCes.push_back(stuple.ele_XCes[i]);
				list.ele_ZCes.push_back(stuple.ele_ZCes[i]);
				list.ele_HadEm.push_back(stuple.ele_HadEm[i]);
				list.ele_Chi2Mean.push_back(stuple.ele_Chi2Mean[i]);
				list.ele_N3d.push_back(stuple.ele_N3d[i]);
				list.ele_Iso4.push_back(stuple.ele_Iso4[i]);
				list.ele_TrkIso.push_back(stuple.ele_TrkIso[i]);		
				list.ele_CesWireE2.push_back(stuple.ele_CesWireE2[i]);
				list.ele_CesStripE2.push_back(stuple.ele_CesStripE2[i]);
				list.ele_PhiWedge.push_back(stuple.ele_PhiWedge[i]);
				list.ele_NMuonStubs.push_back(stuple.ele_NMuonStubs[i]);
				list.ele_EmTime.push_back(stuple.ele_EmTime[i]);
				list.ele_PhoenixId.push_back(stuple.ele_PhoenixId[i]);
				list.ele_Halo_seedWedge.push_back(stuple.ele_Halo_seedWedge[i]);
				list.ele_Halo_eastNhad.push_back(stuple.ele_Halo_eastNhad[i]);
				list.ele_Halo_westNhad.push_back(stuple.ele_Halo_westNhad[i]);
				list.ele_matchJetIndex.push_back(stuple.ele_matchJetIndex[i]);
				list.ele_Ntracks.push_back(stuple.ele_Ntracks[i]);	
				list.ele_Emfr.push_back(stuple.ele_Emfr[i]);
				list.ele_EoverP.push_back(stuple.ele_EoverP[i]);
				list.ele_TrackPt.push_back(stuple.ele_TrackPt[i]);
				list.ele_TrackBcPt.push_back(stuple.ele_TrackBcPt[i]);
				list.ele_TrackPhi.push_back(stuple.ele_TrackPhi[i]);
				list.ele_Nssl.push_back(stuple.ele_Nssl[i]);
				list.ele_Nasl.push_back(stuple.ele_Nasl[i]);
				list.ele_TightId.push_back(stuple.ele_TightId[i]);
				list.ele_LooseId.push_back(stuple.ele_LooseId[i]);
				list.ele_ConversionId.push_back(stuple.ele_ConversionId[i]);
			}

			assert( (int)list.ele_num == (int)list.ele_Index.size());
			break;
		}

		case 1:
		{
			list.ele_num    = stuple.ele_up_num;
			list.ele_Ntight = stuple.ele_up_Ntight;
			list.ele_Nloose = stuple.ele_up_Nloose;

			for (unsigned int i=0; i < stuple.ele_up_num; i++)
			{
				list.ele_Index.push_back(stuple.ele_up_Index[i]);
				list.ele_PhoBlockIndex.push_back(stuple.ele_up_PhoBlockIndex[i]);
				list.ele_EleBlockIndex.push_back(stuple.ele_up_EleBlockIndex[i]);
				list.ele_Etc.push_back(stuple.ele_up_Etc[i]);
				list.ele_E.push_back(stuple.ele_up_E[i]);
				list.ele_Px.push_back(stuple.ele_up_Px[i]);
				list.ele_Py.push_back(stuple.ele_up_Py[i]);
				list.ele_Pz.push_back(stuple.ele_up_Pz[i]);
				list.ele_Detector.push_back(stuple.ele_up_Detector[i]);
				list.ele_DetEta.push_back(stuple.ele_up_DetEta[i]);
				list.ele_DetPhi.push_back(stuple.ele_up_DetPhi[i]);
				list.ele_XCes.push_back(stuple.ele_up_XCes[i]);
				list.ele_ZCes.push_back(stuple.ele_up_ZCes[i]);
				list.ele_HadEm.push_back(stuple.ele_up_HadEm[i]);
				list.ele_Chi2Mean.push_back(stuple.ele_up_Chi2Mean[i]);
				list.ele_N3d.push_back(stuple.ele_up_N3d[i]);
				list.ele_Iso4.push_back(stuple.ele_up_Iso4[i]);
				list.ele_TrkIso.push_back(stuple.ele_up_TrkIso[i]);		
				list.ele_CesWireE2.push_back(stuple.ele_up_CesWireE2[i]);
				list.ele_CesStripE2.push_back(stuple.ele_up_CesStripE2[i]);
				list.ele_PhiWedge.push_back(stuple.ele_up_PhiWedge[i]);
				list.ele_NMuonStubs.push_back(stuple.ele_up_NMuonStubs[i]);
				list.ele_EmTime.push_back(stuple.ele_up_EmTime[i]);
				list.ele_PhoenixId.push_back(stuple.ele_up_PhoenixId[i]);
				list.ele_Halo_seedWedge.push_back(stuple.ele_up_Halo_seedWedge[i]);
				list.ele_Halo_eastNhad.push_back(stuple.ele_up_Halo_eastNhad[i]);
				list.ele_Halo_westNhad.push_back(stuple.ele_up_Halo_westNhad[i]);
				list.ele_matchJetIndex.push_back(stuple.ele_up_matchJetIndex[i]);
				list.ele_Ntracks.push_back(stuple.ele_up_Ntracks[i]);	
				list.ele_Emfr.push_back(stuple.ele_up_Emfr[i]);
				list.ele_EoverP.push_back(stuple.ele_up_EoverP[i]);
				list.ele_TrackPt.push_back(stuple.ele_up_TrackPt[i]);
				list.ele_TrackBcPt.push_back(stuple.ele_up_TrackBcPt[i]);
				list.ele_TrackPhi.push_back(stuple.ele_up_TrackPhi[i]);
				list.ele_Nssl.push_back(stuple.ele_up_Nssl[i]);
				list.ele_Nasl.push_back(stuple.ele_up_Nasl[i]);
				list.ele_TightId.push_back(stuple.ele_up_TightId[i]);
				list.ele_LooseId.push_back(stuple.ele_up_LooseId[i]);
				list.ele_ConversionId.push_back(stuple.ele_up_ConversionId[i]);
			}

			assert(list.ele_num == list.ele_Index.size());
			break;

		}

		case -1:
		{
			list.ele_num    = stuple.ele_down_num;
			list.ele_Ntight = stuple.ele_down_Ntight;
			list.ele_Nloose = stuple.ele_down_Nloose;

			for (unsigned int i=0; i < stuple.ele_down_num; i++)
			{
				list.ele_Index.push_back(stuple.ele_down_Index[i]);
				list.ele_PhoBlockIndex.push_back(stuple.ele_down_PhoBlockIndex[i]);
				list.ele_EleBlockIndex.push_back(stuple.ele_down_EleBlockIndex[i]);
				list.ele_Etc.push_back(stuple.ele_down_Etc[i]);
				list.ele_E.push_back(stuple.ele_down_E[i]);
				list.ele_Px.push_back(stuple.ele_down_Px[i]);
				list.ele_Py.push_back(stuple.ele_down_Py[i]);
				list.ele_Pz.push_back(stuple.ele_down_Pz[i]);
				list.ele_Detector.push_back(stuple.ele_down_Detector[i]);
				list.ele_DetEta.push_back(stuple.ele_down_DetEta[i]);
				list.ele_DetPhi.push_back(stuple.ele_down_DetPhi[i]);
				list.ele_XCes.push_back(stuple.ele_down_XCes[i]);
				list.ele_ZCes.push_back(stuple.ele_down_ZCes[i]);
				list.ele_HadEm.push_back(stuple.ele_down_HadEm[i]);
				list.ele_Chi2Mean.push_back(stuple.ele_down_Chi2Mean[i]);
				list.ele_N3d.push_back(stuple.ele_down_N3d[i]);
				list.ele_Iso4.push_back(stuple.ele_down_Iso4[i]);
				list.ele_TrkIso.push_back(stuple.ele_down_TrkIso[i]);		
				list.ele_CesWireE2.push_back(stuple.ele_down_CesWireE2[i]);
				list.ele_CesStripE2.push_back(stuple.ele_down_CesStripE2[i]);
				list.ele_PhiWedge.push_back(stuple.ele_down_PhiWedge[i]);
				list.ele_NMuonStubs.push_back(stuple.ele_down_NMuonStubs[i]);
				list.ele_EmTime.push_back(stuple.ele_down_EmTime[i]);
				list.ele_PhoenixId.push_back(stuple.ele_down_PhoenixId[i]);
				list.ele_Halo_seedWedge.push_back(stuple.ele_down_Halo_seedWedge[i]);
				list.ele_Halo_eastNhad.push_back(stuple.ele_down_Halo_eastNhad[i]);
				list.ele_Halo_westNhad.push_back(stuple.ele_down_Halo_westNhad[i]);
				list.ele_matchJetIndex.push_back(stuple.ele_down_matchJetIndex[i]);
				list.ele_Ntracks.push_back(stuple.ele_down_Ntracks[i]);	
				list.ele_Emfr.push_back(stuple.ele_down_Emfr[i]);
				list.ele_EoverP.push_back(stuple.ele_down_EoverP[i]);
				list.ele_TrackPt.push_back(stuple.ele_down_TrackPt[i]);
				list.ele_TrackBcPt.push_back(stuple.ele_down_TrackBcPt[i]);
				list.ele_TrackPhi.push_back(stuple.ele_down_TrackPhi[i]);
				list.ele_Nssl.push_back(stuple.ele_down_Nssl[i]);
				list.ele_Nasl.push_back(stuple.ele_down_Nasl[i]);
				list.ele_TightId.push_back(stuple.ele_down_TightId[i]);
				list.ele_LooseId.push_back(stuple.ele_down_LooseId[i]);
				list.ele_ConversionId.push_back(stuple.ele_down_ConversionId[i]);
			}

			assert(list.ele_num == list.ele_Index.size());
			break;
		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid collection requested! can make central/em up/em down  lists only!");
			assert(false);
		}
	}
/*}}}*/
}	//CreateElectronLists

void PhotonJets2::CreateJetLists(const Stuple stuple, JetList& list, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
		  list.jet_num  = stuple.jet_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_NJet15;
			for (unsigned int i=0; i < stuple.jet_num; i++)
			{
				if ( fabs(stuple.jet_Px[i]) > 1500 || 
					  fabs(stuple.jet_Py[i]) > 1500 || 
					  fabs(stuple.jet_Pz[i]) > 1500 || 
					  fabs(stuple.jet_E[i]) > 1500)
				{
					stuple.DumpJetBlock();
					assert(false);
				}

				list.jet_Index.push_back(stuple.jet_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_E[i]);
				list.jet_Px.push_back(stuple.jet_Px[i]);
				list.jet_Py.push_back(stuple.jet_Py[i]);
				list.jet_Pz.push_back(stuple.jet_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_SeedIEta[i]);
			}

			assert(list.jet_num == list.jet_Index.size());
			break;
		}

		case 1:
		{
		  list.jet_num  = stuple.jet_up_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_up_NJet15;
			for (unsigned int i=0; i < stuple.jet_up_num; i++)
			{
				list.jet_Index.push_back(stuple.jet_up_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_up_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_up_E[i]);
				list.jet_Px.push_back(stuple.jet_up_Px[i]);
				list.jet_Py.push_back(stuple.jet_up_Py[i]);
				list.jet_Pz.push_back(stuple.jet_up_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_up_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_up_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_up_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_up_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_up_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_up_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_up_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_up_SeedIEta[i]);
				
			}

			assert(list.jet_num == list.jet_Index.size());
			break;

			
		}

		case -1:
		{
		  list.jet_num  = stuple.jet_down_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_down_NJet15;
			for (unsigned int i=0; i < stuple.jet_down_num; i++)
			{
				list.jet_Index.push_back(stuple.jet_down_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_down_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_down_E[i]);
				list.jet_Px.push_back(stuple.jet_down_Px[i]);
				list.jet_Py.push_back(stuple.jet_down_Py[i]);
				list.jet_Pz.push_back(stuple.jet_down_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_down_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_down_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_down_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_down_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_down_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_down_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_down_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_down_SeedIEta[i]);
				
			}

			assert(list.jet_num == list.jet_Index.size());
			break;

		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid collection requested! can make central and JES up/down lists only!");
			assert(false);

		}

	}
/*}}}*/
} // CreateJetLists

void PhotonJets2::GetCommonVariables(const Stuple stuple, CommonVars& cVars)
{
/*{{{*/
	cVars.evt_McFlag 		= stuple.evt_McFlag;
	cVars.evt_RunNumber 	= stuple.evt_RunNumber;
	cVars.evt_EventNumber= stuple.evt_EventNumber;
	cVars.tri_pho25iso 	= stuple.tri_pho25iso;
	cVars.tri_pho50 		= stuple.tri_pho50;
	cVars.tri_pho70 		= stuple.tri_pho70;

	cVars.vtx_N 			= stuple.vtx_N;
	cVars.vtx_NClass12 	= stuple.vtx_NClass12;
	cVars.vtx_z 			= stuple.vtx_z;
	cVars.vtx_Ntracks 	= stuple.vtx_Ntracks;
	cVars.vtx_SumPt 		= stuple.vtx_SumPt;

	cVars.met_Met 			= stuple.met_Met;	
	cVars.met_SumEt 		= stuple.met_SumEt;
	cVars.met_Ht 			= stuple.met_Ht;	
	cVars.met_MetPhi 		= stuple.met_MetPhi;
//	cVars.met_Gen_d 		= stuple.met_Gen_d;	
//	cVars.met_Gen_m 		= stuple.met_Gen_m;
//	cVars.met_Gen_p 		= stuple.met_Gen_p;
//	cVars.met_Gen_mUn 	= stuple.met_Gen_mUn;
//	cVars.met_Gen_pUn 	= stuple.met_Gen_pUn;
/*}}}*/
}


void PhotonJets2::DoMyStuff(CommonVars cVars, PhotonList cPhos, PhotonList uPhos, PhotonList dPhos,
									ElectronList cEles, ElectronList uEles, ElectronList dEles,
									JetList cJets, JetList uJets, JetList dJets)
{	


	std::cout << "___begin DOMyStuff " << std::endl;

		//cPhos.Dump();
		cJets.Dump();

	std::cout << "_____________________ " << std::endl;

	if (!bMcFlag)   		//_________________ 0=data, 1=mc
	{
		// DATA path
		
		//______________________________________________________________ trigger
		if (  (cVars.tri_pho25iso != 1) && (cVars.tri_pho50 != 1) 	//_______________ only data has L1-3 trig. info
				&& (cVars.tri_pho70 != 1) ) 										//_______________ mc has L1-2.
			return;												
			
		//require at least one tight photon
		//remove e+pho events
		//if (stuple.pho_Ntight < 1) return;			//this will screw up the qcd template. so drop this
		//if (cEles.ele_Ntight > 0) return;
		
		int iPhoIndex = -1; //leading photon candidate
		int iNele = 0;
		EvtTag_t evt;
		evt.Run 		= cVars.evt_RunNumber;
		evt.Evt 		= cVars.evt_EventNumber;
		evt.Halo 	= 0;
		evt.Cosmic 	= 0;
		evt.Signal 	= 0;
		evt.QCD 		= 0;
	
		//______________________________________________________________ select the photon
		for (unsigned int i = 0; i < cPhos.pho_num; ++i )
		{
			if (cPhos.pho_PhoenixId[i] != 0)
			{
				++iNele;
				continue;  //____________________________ reject phoenix photons
			}
			//________________________________________________________________________ this is a hack to find the trigger photon, as I don't have the L3 info to match reconstructed photon to trigger photon.
			if (cVars.tri_pho70 == 1 && cPhos.pho_Etc[i] < 70.0) continue;
			if (cVars.tri_pho50 == 1 && cPhos.pho_Etc[i] < 50.0) continue;
			if (cVars.tri_pho25iso == 1 && cPhos.pho_Etc[i] < 25.0) continue;

			if (cPhos.pho_Etc[i]<30) continue;			//need this as I am tagging Et>7GeV em objects. Need to update rest of this code
																	// when looking for electrons. 11-19-2009
			if (fabs(cPhos.DetEta) > 1.1) continue;	//central
			
			if (cPhos.pho_TightId[i] == 0)
			{ 

				//cosmic selection
				if ( cPhos.pho_EmTime[UsedPho[0]] > 30 
						&& cPhos.pho_EmTime[UsedPho[0]] < 90  //em time cut
						&& cVars.vtx_NClass12 >= 1 
						&& fabs(cVars.vtx_z) < 60.	//vtx cut
						)
				{
					JetSelAndHistFill(cVars, cPhos, cEles, cJets, &Hcosmic, i, 
							vCosmicCount, false);
				 	evt.Cosmic = 1;
					breakl
					
				} else if ( fabs(cPhos.pho_EmTime[UsedPho[0]]) < 4.8)
				{
					
					//_________________________________________ pick events from no vtx for halo 
					if (cVars.vtx_NClass12 == 0)
					{
						//_______________________________________________________________ select for Halo template
						HaloTypes ht(cPhos.pho_Halo_seedWedge[UsedPho[0]], 
								cPhos.pho_Halo_eastNhad[UsedPho[0]] + cPhos.pho_Halo_westNhad[UsedPho[0]]);

						if (ht.IsType(5))
						{
							if (cPhos.pho_PhiWedge[0] == 0 || cPhos.pho_PhiWedge[0] == 23)
							{
								used = JetSelAndHistFill(cVars, cPhos, cEles, cJets,
										&Hhalo, 1, vHaloCount, false); 
								evt.Halo = 1;
								break;
							}
						}
					
					} else if (cVars.vtx_NClass12 >= 1)
					{
						if (fabs(cVars.vtx_z) < 60.)
						{
							std::cout << "SIGNAL JETS" << std::endl;
							cJets.Dump();
							int tmp = JetSelAndHistFill(cVars, cPhos, cEles, cJets,
													&Hsignal, i, vSignalCount, false);

							if (cVars.tri_pho70 == 1) phoEt_70->Fill(cPhos.pho_Etc.at(UsedPho.at(0))); 
							else if (cVars.tri_pho50 == 1) phoEt_50->Fill(cPhos.pho_Etc.at(UsedPho.at(0))); 
							else if (cVars.tri_pho25iso == 1) phoEt_iso25->Fill(cPhos.pho_Etc.at(UsedPho.at(0))); 
							evt.Signal = 1;
							break;
						}
						
					}

				}
			} else if (cPhos.pho_TightId[i] > 0 && cPhos.pho_LooseId[i] == 0)			// >0 means that this has failed tight cuts. must do this as I init them to large negative values and != could mean <0!!
			{
				//sideband selection (reject bh) ewk?
				HaloTypes ht(cPhos.pho_Halo_seedWedge[UsedPho[0]], 
						cPhos.pho_Halo_eastNhad[UsedPho[0]] + cPhos.pho_Halo_westNhad[UsedPho[0]]);

				if ( fabs(cPhos.pho_EmTime[UsedPho[0]]) < 4.8   // same as signal
						&& cVars.vtx_NClass12 >= 1 && cVars.vtx_z < 60.   //same as signal
						&&  ! ht.IsType(5)	// reject BH
						&& iNele == 0 		//reject lepton+jets
					)
				{

						int tmp = JetSelAndHistFill(cVars, qcdcPhos, qcdcEles, qcdcJets,
													&Hqcd, i, vQCDCount, false);
						evt.QCD = 1;
						break;

				}
			}
		}	
	
		// now store the summary of this event
		EVT.push_back(evt);
		
	} else {
/*
		// I had this at beginnig of this method. this causes NO BH events as they are
		// choosen from nvtx==0 events! So had to apply these cuts here.
		//12-02-2009
		//_________________________________________________________________ require a one good vertex 		  
		if (cVars.vtx_NClass12 <= 0) return;
		if (cVars.vtx_z > 60.)  return;

		
		// MC path
		std::vector<int> UsedPhoc, UsedElec;
		std::vector<int> UsedPhou, UsedEleu;
		std::vector<int> UsedPhod, UsedEled;

		
		// CENTRAL
		//first select the photon, then add un-used em obj to jet list, then fill the hists
		SelectPhoton(cVars, cPhos, cEles, cJets, UsedPhoc, UsedElec);
		if (UsedPhoc.size() > 0)
		{
			//________________________________________________________________ now put back the un-used em objs to the jet list and reorder the jets
			NewJetList nj;
			nj.AddUnused(cPhos, cEles, cJets, UsedPhoc, UsedElec);
			//________________________________________________________________ Select Jets and fill hists
			JetSelAndHistFill(cVars, cPhos, cEles, cJets, &Hmc, UsedPhoc, UsedElec, vMcCount, false);
		}

		// JES AND EM UP
		//first select the photon, then add un-used em obj to jet list, then fill the hists
		SelectPhoton(cVars, uPhos, uEles, uJets, UsedPhou, UsedEleu);
		if (UsedPhou.size() > 0)
		{
			//________________________________________________________________ now put back the un-used em objs to the jet list and reorder the jets
			NewJetList nj;
			nj.AddUnused(uPhos, uEles, uJets, UsedPhou, UsedEleu);
			//________________________________________________________________ Select Jets and fill hists
			JetSelAndHistFill(cVars, uPhos, uEles, uJets, &HmcUp, UsedPhou, UsedEleu, vMcUpCount, false);
		}

		// JES AND EM DOWN
		//first select the photon, then add un-used em obj to jet list, then fill the hists
		SelectPhoton(cVars, dPhos, dEles, dJets, UsedPhod, UsedEled);
		if (UsedPhod.size() > 0)
		{
			//________________________________________________________________ now put back the un-used em objs to the jet list and reorder the jets
			NewJetList nj;
			nj.AddUnused(dPhos, dEles, dJets, UsedPhod, UsedEled);
			//________________________________________________________________ Select Jets and fill hists
			JetSelAndHistFill(cVars, dPhos, dEles, dJets, &HmcDown, UsedPhod, 
					UsedEled, vMcDownCount, false);
		}
*/
	}	// if data/mc
		

}


void PhotonJets2::SelectPhoton(const CommonVars Vars, const PhotonList phos,
									const ElectronList eles, const JetList jets, 
									std::vector<int>& UsedPho, std::vector<int>& UsedEle)
{

		//after general evt selection (run, trigger, vtx) everything else should be same as signal selection
		//require at least one tight photon
		//remove e+pho events
		if (phos.pho_Ntight < 1) return;
		
		//______________________________________________________________ select the photon
		for (unsigned int i = 0; i < phos.pho_num; ++i ) {
			if (phos.pho_TightId[i] != 0) continue;		// ___________________________ tight photons
			if (phos.pho_PhoenixId[i] != 0) continue;  //____________________________ reject phoenix photons
			UsedPho.push_back(i);	
			break;		//______________________________ one is enough. and this is the highest tight Et photon in the event as I have sorted the SuperPhoton list.
		}	
	
		if (UsedPho.size() > 0) {					//______ no photons were found
			//_______________________________________________________________ reject if for Halo template
			HaloTypes ht(phos.pho_Halo_seedWedge[UsedPho[0]], 
					phos.pho_Halo_eastNhad[UsedPho[0]] + phos.pho_Halo_westNhad[UsedPho[0]]);

			if (ht.IsType(5)) {
				UsedPho.clear();
			}
		}

}



int PhotonJets2::JetSelAndHistFill(const CommonVars Vars, const PhotonList phos, 
											const ElectronList eles, 
											const JetList jets, Hist_t* hist,
											const int iPhoIndex, 
											std::vector<unsigned int>& count, 
											const bool UseCprWgtForHt, 
											const float fWeight)
{
	int used = 0;
	
	if (vPhoIndex.size() < 1) {
		StdOut(__FILE__,__LINE__,3,"::No photons were selected to make plots!");
		assert(false);
	}
	if (count.size() != 3) {
		StdOut(__FILE__,__LINE__,3,"::counter vector size must be 3!");
		assert(false);
	}
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"::Hist pointer is NULL!");
		assert(false);

	}
	
	int iPhoIndex = vPhoIndex[0];	//assuming the lead photon index is stored first

	++count[0];
	
	//if (jets.jet_NJet15 >= 1) {	// mono-jet case
	
	//to study exclusive njet==1 case for Ht plot problem to see if 
	//the discrepency is due to the MC higher order modeling failure. - 11-27-2009
	std::cout << __LINE__ << ":NJET15= "<< jets.jet_NJet15 << std::endl;
	if (jets.jet_NJet15 >= 1) {	// mono-jet case
		++count[1];
		used = 1;
		FillPhotonHists(Vars, phos, &(hist->p1j_Pho), iPhoIndex, fWeight);
		FillEventHists(Vars, phos, eles, jets, &(hist->p1j_Evt), iPhoIndex, fWeight, UseCprWgtForHt);
		FillJetHists(Vars, jets, &(hist->p1j_Jet),0, fWeight);		//lead jet
		FillPhoton1JetHists(Vars, phos, jets, &(hist->p1j_PhoJet),iPhoIndex, 0, fWeight);
		
	}



/*	if (jets.jet_NJet15 >=2) {		// di-jet case
		++count[2];
		used = 1;
		FillPhotonHists(Vars, phos, &(hist->p2j_Pho), iPhoIndex, fWeight);
		FillEventHists(Vars, phos, eles, jets, &(hist->p2j_Evt), iPhoIndex, fWeight, UseCprWgtForHt);
		FillJetHists(Vars, jets, &(hist->p2j_Jet1),0, fWeight);		//lead jet
		FillJetHists(Vars, jets, &(hist->p2j_Jet2),1, fWeight);		//2nd lead jet
		FillPhoton1JetHists(Vars, phos, jets, &(hist->p2j_PhoJet1),iPhoIndex, 0, fWeight);
		FillPhoton1JetHists(Vars, phos, jets, &(hist->p2j_PhoJet2),iPhoIndex, 1, fWeight);
		FillPhoton2JetsHists(Vars, phos, jets, &(hist->p2j_PhoJets),iPhoIndex, 0,1, fWeight);
		FillTwoJetsHists(Vars, jets, &(hist->p2j_TwoJets), 0,1, fWeight);
	}
*/
	return used;		//1=this event has been used!
}


void PhotonJets2::FillEventHists(const CommonVars& vars, const PhotonList& phos,
								const ElectronList& eles, const JetList& jets, 
								Histograms::EventHists_t* hist, const int iPhoIndex,
								const float fWeight, const bool UseCprWgtForHt)
{ 
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	hist->Met->Fill(vars.met_Met, fWeight);
	hist->Met_Gen_d->Fill(vars.met_Gen_d, fWeight);
	hist->Sumet->Fill(vars.met_SumEt, fWeight);
	hist->NVertices->Fill(vars.vtx_N, fWeight);
	hist->N12Vertices->Fill(vars.vtx_NClass12, fWeight);
	hist->NJet15->Fill(jets.jet_NJet15, fWeight);
	hist->MetRun->Fill(vars.evt_RunNumber, vars.met_Met, fWeight);
	hist->SumetRun->Fill(vars.evt_RunNumber, vars.met_SumEt, fWeight);
	hist->NVerticesRun->Fill(vars.evt_RunNumber, vars.vtx_N, fWeight);
	hist->N12VerticesRun->Fill(vars.evt_RunNumber, vars.vtx_NClass12, fWeight);
	hist->NJet15Run->Fill(vars.evt_RunNumber, jets.jet_NJet15, fWeight);
	hist->NJet15VsMet->Fill(vars.met_Met, jets.jet_NJet15, fWeight);
	hist->NJet15VsSumet->Fill(vars.met_SumEt, jets.jet_NJet15, fWeight);
	hist->N12VerticesVsMet->Fill(vars.met_Met, vars.vtx_NClass12, fWeight);
	hist->N12VerticesVsSumet->Fill(vars.met_SumEt, vars.vtx_NClass12, fWeight);
	if (vars.tri_pho25iso == 1) hist->Triggers->Fill(0);
	if (vars.tri_pho50 == 1) hist->Triggers->Fill(1);
	if (vars.tri_pho70 == 1) hist->Triggers->Fill(2);

	hist->Ht->Fill(vars.met_Ht, fWeight);
	if (UseCprWgtForHt)			// these will be filled only if CPR weight is required
	{
		
//		hist->HtWgt->Fill(vars.met_Ht, phos.pho_CprWgt.at(iPhoIndex));	
//		hist->HtSys[0]->Fill(vars.met_Ht, phos.pho_CprSys1.at(iPhoIndex));
//		hist->HtSys[1]->Fill(vars.met_Ht, phos.pho_CprSys2.at(iPhoIndex));
//		hist->HtSys[2]->Fill(vars.met_Ht, phos.pho_CprSys3.at(iPhoIndex));
//		hist->HtSys[3]->Fill(vars.met_Ht, phos.pho_CprSys4.at(iPhoIndex));
//		hist->HtSys[4]->Fill(vars.met_Ht, phos.pho_CprSys5.at(iPhoIndex));
//		hist->HtSys[5]->Fill(vars.met_Ht, phos.pho_CprSys6.at(iPhoIndex));
//		hist->HtSys[6]->Fill(vars.met_Ht, phos.pho_CprSys7.at(iPhoIndex));
//		hist->HtSys[7]->Fill(vars.met_Ht, phos.pho_CprSys8.at(iPhoIndex));
	}




}

void PhotonJets2::FillPhotonHists(const CommonVars& cVars, const PhotonList& phos, Histograms::PhotonHists_t* hist,
								const int iPhoIndex, const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	hist->Detector->Fill(phos.pho_Detector[iPhoIndex], fWeight);
	hist->DetEta->Fill(phos.pho_DetEta[iPhoIndex], fWeight);
	hist->DetPhi->Fill(phos.pho_DetPhi[iPhoIndex], fWeight);
	hist->EtCorr->Fill(phos.pho_Etc[iPhoIndex], fWeight);
	hist->XCes->Fill(phos.pho_XCes[iPhoIndex], fWeight);
	hist->ZCes->Fill(phos.pho_ZCes[iPhoIndex], fWeight);
	hist->HadEm->Fill(phos.pho_HadEm[iPhoIndex], fWeight);
	hist->Chi2Mean->Fill(phos.pho_Chi2Mean[iPhoIndex], fWeight);
	hist->N3d->Fill(phos.pho_N3d[iPhoIndex], fWeight);
	hist->Iso4->Fill(phos.pho_Iso4[iPhoIndex], fWeight);
	hist->TrkPt->Fill(phos.pho_TrkPt[iPhoIndex], fWeight);
	hist->TrkIso->Fill(phos.pho_TrkIso[iPhoIndex], fWeight);
	hist->Ces2Wire->Fill(phos.pho_CesWireE2[iPhoIndex], fWeight);
	hist->Ces2Strip->Fill(phos.pho_CesStripE2[iPhoIndex], fWeight);
	hist->EmTime->Fill(phos.pho_EmTime[iPhoIndex], fWeight);
	hist->PhiWedge->Fill(phos.pho_PhiWedge[iPhoIndex], fWeight);
	hist->EmTimeVsRun->Fill(cVars.evt_RunNumber, phos.pho_EmTime[iPhoIndex], fWeight);
	hist->EtCorrVsRun->Fill(cVars.evt_RunNumber, phos.pho_Etc[iPhoIndex], fWeight);
	if (phos.pho_CprWgt[iPhoIndex] != 0) //0 weight means the photon failed the 
	{												// CES/CPR fiducial cuts. see elog#1491 and #1495
													// these photon need to be excluded to get the
													// correct photon fraction. 
		hist->CprWeight->Fill(phos.pho_CprWgt[iPhoIndex]);	
	}
}
	

void PhotonJets2::FillJetHists(const CommonVars& cVars, const JetList& jets, Histograms::JetHists_t* hist, 
								const int iJetIndex,	const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	hist->Emfr->Fill(jets.jet_Emfr[iJetIndex], fWeight);
	hist->DetEta->Fill(jets.jet_DetEta[iJetIndex], fWeight);
	hist->NTracks->Fill(jets.jet_Ntracks[iJetIndex], fWeight);
	hist->NTowers->Fill(jets.jet_Ntowers[iJetIndex], fWeight);
	hist->EtCorr->Fill(jets.jet_Pt[iJetIndex], fWeight);

	TLorentzVector tlVec(jets.jet_Px[iJetIndex], jets.jet_Py[iJetIndex], jets.jet_Pz[iJetIndex], jets.jet_E[iJetIndex]);
	hist->EvtPhi->Fill(tlVec.Phi(), fWeight);
	hist->EvtEta->Fill(tlVec.Eta(), fWeight);
	//hist->SeedIEta->Fill(jets.jet_SeedIEta[iJetIndex], fWeight);
	//hist->SeedIPhi->Fill(jets.jet_SeedIPhi[iJetIndex], fWeight);
}

void PhotonJets2::FillPhoton1JetHists(const CommonVars& cVars, const PhotonList& phos, const JetList& jets, 
							Histograms::Photon1JetHists_t* hist, const int iPhoIndex, const int iJetIndex,
							const float fWeight)
{
	//std::cout << "\tfillpho1jethist:: npho=" << phos.pho_num << std::endl;
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	if (iPhoIndex < 0 || iPhoIndex >= (int)phos.pho_num)
	{
		StdOut(__FILE__,__LINE__,3,"Invalid Photon index!");
		assert(false);
	}
	if (iJetIndex < 0 || iJetIndex >= (int)jets.jet_num)
	{
		StdOut(__FILE__,__LINE__,3,"Invalid Jet index!");
		assert(false);
	}
	TLorentzVector tlPho, tlJet, tlSum;
	float InvMass 	= 0;
	float DelPhi 	= 0;
	float DelEta 	= 0;
	float DelR 		= 0;
	float EtRatio 	= 0;
	float JetEmFr 	= 0;
	float PhoPhiWedge = 0; 

	std::cout << "-----------jet dump in hist fill " << std::endl;
	jets.Dump();

	tlJet.SetPxPyPzE(jets.jet_Px[iJetIndex], jets.jet_Py[iJetIndex], jets.jet_Pz[iJetIndex], jets.jet_E[iJetIndex]);
	tlPho.SetPxPyPzE(phos.pho_Px[iPhoIndex], phos.pho_Py[iPhoIndex], phos.pho_Pz[iPhoIndex], phos.pho_E[iPhoIndex]);

	tlSum = tlPho + tlJet; 
	InvMass 	= tlSum.M();
	if (InvMass< 0) 
	{
		std::cout << __FILE__ <<":" << __FUNCTION__ << ":" << __LINE__ <<	"::Event with InvMass = " << InvMass << " <0" << std::endl;
		std::cout << "Run,Event Number::" << cVars.evt_RunNumber << ", " << cVars.evt_EventNumber << std::endl;
		std::cout << "Jet Pt/E = " << tlJet.Pt() << ", " << tlJet.E() << std::endl;
		std::cout << "Pho Pt/E = " << tlPho.Pt() << ", " << tlPho.E() << std::endl;
		std::cout << "Sum Pt/E = " << tlSum.Pt() << ", " << tlSum.E() << std::endl;
		jets.Dump();
	}
	
	//DelPhi 	= fabs(tlPho.DeltaPhi(tlJet));
	DelEta 	= fabs(phos.pho_DetEta[iPhoIndex] - jets.jet_DetEta[iJetIndex]);
	DelR 		= tlPho.DeltaR(tlJet);
	EtRatio 	= (tlJet.Energy() * TMath::Sin(tlJet.Theta())) / (tlPho.Energy() * TMath::Sin(tlPho.Theta()));
	//JetEmFr 	= jets.jet_Emfr[iJetIndex];
	PhoPhiWedge = phos.pho_PhiWedge[iPhoIndex];
	
	hist->InvMass->Fill(InvMass, fWeight);
	//hist->DelPhi->Fill(DelPhi, fWeight);
	hist->DelEta->Fill(DelEta, fWeight);
	hist->DelR->Fill(DelR, fWeight);
	hist->EtRatio->Fill(EtRatio, fWeight);
	//hist->JetEmfrVsDelPhi->Fill(DelPhi, JetEmFr, fWeight);
	//hist->JetEmfrVsPhoPhiWedge->Fill(PhoPhiWedge, JetEmFr, fWeight);

}


void PhotonJets2::FillPhoton2JetsHists(const CommonVars& cVars, const PhotonList& phos,
						const JetList& jets, Histograms::Photon2JetsHists_t* hist,
						const int iPhoIndex, const int iJet1Index,
						const int iJet2Index, const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	TLorentzVector tlJet1, tlJet2, tlPho, tlSum;

	tlJet1.SetPxPyPzE(jets.jet_Px[iJet1Index], jets.jet_Py[iJet1Index], jets.jet_Pz[iJet1Index], jets.jet_E[iJet1Index]);
	tlJet2.SetPxPyPzE(jets.jet_Px[iJet2Index], jets.jet_Py[iJet2Index], jets.jet_Pz[iJet2Index], jets.jet_E[iJet2Index]);
	tlPho.SetPxPyPzE(phos.pho_Px[iPhoIndex], phos.pho_Py[iPhoIndex], phos.pho_Pz[iPhoIndex], phos.pho_E[iPhoIndex]);
	
	tlSum = tlPho + tlJet1 + tlJet2; 
	hist->InvMass->Fill(tlSum.M(), fWeight);

}

void PhotonJets2::FillTwoJetsHists(const CommonVars& cVars, const JetList& jets, Histograms::TwoJetsHists_t* hist,
						const int iJet1Index, const int iJet2Index,
						const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	TLorentzVector tlJet1, tlJet2, tlSum;
	float DelEta = 0;

	tlJet1.SetPxPyPzE(jets.jet_Px[iJet1Index], jets.jet_Py[iJet1Index], jets.jet_Pz[iJet1Index], jets.jet_E[iJet1Index]);
	tlJet2.SetPxPyPzE(jets.jet_Px[iJet2Index], jets.jet_Py[iJet2Index], jets.jet_Pz[iJet2Index], jets.jet_E[iJet2Index]);
	DelEta = fabs(jets.jet_DetEta[iJet1Index] - jets.jet_DetEta[iJet2Index]);
	
	tlSum = tlJet1 + tlJet2; 
	float fEtRatio = (tlJet2.Energy() * TMath::Sin(tlJet2.Theta())) / (tlJet1.Energy() * TMath::Sin(tlJet1.Theta()));
	float DelPhi = fabs(tlJet1.DeltaPhi(tlJet2));
	float DelR = tlJet1.DeltaR(tlJet2);
	
	hist->InvMass->Fill(tlSum.M(), fWeight);
	hist->DelPhi->Fill(DelPhi, fWeight);
	hist->DelEta->Fill(DelEta, fWeight);
	hist->DelR->Fill(DelR, fWeight);
	hist->EtRatio->Fill(fEtRatio, fWeight);

}

void PhotonJets2::PrintHeader(const CommonVars& cVars)
{
	std::cout << " ============== " << cVars.evt_RunNumber << ", " << cVars.evt_EventNumber << std::endl;
}
