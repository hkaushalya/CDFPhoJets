#include <iostream>
#include <sstream>
#include "PhoFlat/PhoenixPhoStudy.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
#include "TBenchmark.h"

//-------------------------------------------------------------------
PhoenixPhoStudy::PhoenixPhoStudy():
iProgressBy(50000),
sHistFileName("PhoenixPhoStudy.root"),
bMcFlag(true),
bRejectPhoenix(true),
sPhxStatus(" After Phoenix Rejection ")
{
}

//-------------------------------------------------------------
void PhoenixPhoStudy::Init(TChain *tree)
{
	
	readIn = new ReadInAll(tree, &stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->Met(1);

	for (int i=0; i<3; i++) {
		vQCDCount.push_back(0);
		vSignalCount.push_back(0);
		vPhxCount.push_back(0);
	}

	readIn->Init();

	//cannot add folder without creating a dir first!
	histFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = histFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

} // Init

//-------------------------------------------------------------------
void PhoenixPhoStudy::CleanUp()
{
	vSignalCount.clear();
	vQCDCount.clear();
	EVT.clear();
	
	histFile->Write();
	histFile->Close();
} 


//------ main -------------------------------------------------------
void PhoenixPhoStudy::Main(TChain *ch, int iRunEvents)
{
	
	TBenchmark timer;
	timer.Start("phoana_time");
	
	if (ch) {
		myChain = ch;
		//ch->Print();	
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
	
		//std::cout << "pho_num=" << stuple.pho_num << std::endl;
		//std::cout << "pho_num status=" << myChain->GetBranchStatus("pho_num") << std::endl;
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;
	
		//CheckOrder(stuple);
		DoMyStuff(stuple);

	}
	
	std::cout << "======== PhoenixPhoStudy =======================" << std::endl;
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
	//std::cout << "Phoenix reject status  = " << bRejectPhoenix << " (" << sPhxStatus <<")" << std::endl;
	std::cout << "-------------- PHOENIX ---------------------- "<< std::endl;
	std::cout << "Phx Evts               = " << vPhxCount[0] << std::endl;
	std::cout << "Phx+>=1Jet Evts        = " << vPhxCount[1] << std::endl; 
	std::cout << "Phx+>=2Jet Evts        = " << vPhxCount[2] << std::endl; 
	std::cout << "-------------- QCD ------------------------- "<< std::endl;
	std::cout << "QCD Evts               = " << vQCDCount[0] << std::endl;
	std::cout << "QCD+>=1Jet Evts        = " << vQCDCount[1] << std::endl; 
	std::cout << "QCD+>=2Jet Evts        = " << vQCDCount[2] << std::endl; 
	std::cout << "-------------- SIGNAL----------------------- "<< std::endl;
	std::cout << "Pho Evts               = " << vSignalCount[0] << std::endl;
	std::cout << "Pho+>=1Jet Evts        = " << vSignalCount[1] << std::endl; 
	std::cout << "Pho+>=2Jet Evts        = " << vSignalCount[2] << std::endl; 
	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	std::cout << "======= END PhoenixPhoStudy ====================" << std::endl;

	
	CleanUp();
	timer.Show("phoana_time");
} // Main


//-------------------------------------------------------------------
void PhoenixPhoStudy::BookHistograms()
{
/*{{{*/
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


		// PHOENIX PHOTON PLOTS
		text2title1 = "#gamma_{Phoenix} + >=1Jet :";
		text2title2 = "#gamma_{Phoenix} + >=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("PHOENIX","Phoenix Photons");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hphx.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hphx.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hphx.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hphx.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hphx.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hphx.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hphx.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hphx.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hphx.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hphx.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hphx.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hphx.p2j_TwoJets,text2title2.c_str());

		
	
/*}}}*/		
} //BookHistograms

void PhoenixPhoStudy::DoMyStuff(Stuple st)
{	
		EvtTag_t evt;
		evt.Run 		= st.evt_RunNumber;
		evt.Evt 		= st.evt_EventNumber;
		evt.Signal 	= 0;
		evt.QCD 		= 0;
		evt.Phx 		= 0;

	//_________________________________________________________________ require a one good vertex 		  
	if (st.vtx_NClass12 <= 0) return;
	if (st.vtx_z > 60.)  return;

		//______________________________________________________________ trigger
		if (!bMcFlag)
		{
			if (  (st.tri_pho25iso != 1) && (st.tri_pho50 != 1) 	//_______________ only data has L1-3 trig. info
					&& (st.tri_pho70 != 1) ) 										//_______________ mc has L1-2.
				return;												
		}
			
		std::vector<int> vTPhosInd, vSPhosInd, vPhxInd;
	
		//______________________________________________________________ select the photon
		for (unsigned int i = 0; i < st.pho_num; ++i )
		{
			//if (bRejectPhoenix)
			//{
			//	if (st.pho_PhoenixId[i] != 0) continue;  //____________________________ reject phoenix photons
			//}
			//________________________________________________________________________ this is a hack to find the trigger photon, as I don't have the L3 info to match reconstructed photon to trigger photon.
			//if (st.tri_pho70 == 1 && st.pho_Etc[i] < 70.0) continue;
			//if (st.tri_pho50 == 1 && st.pho_Etc[i] < 50.0) continue;
			//if (st.tri_pho25iso == 1 && st.pho_Etc[i] < 25.0) continue;
			if (st.pho_Etc[i] < 180 || st.pho_Etc[i]>200) continue;
			//std::cout << "run, evt=" << st.evt_RunNumber << ", " << st.evt_EventNumber << std::endl;

			if (st.pho_PhoenixId[i] != 0)
			{
				vPhxInd.push_back(i);
				continue;  //____________________________ reject phoenix photons
			}

			if (st.pho_TightId[i] == 0) vTPhosInd.push_back(i);
			else if (st.pho_LooseId[i] == 0) vSPhosInd.push_back(i);
		}	


		if (vTPhosInd.size()) { evt.Signal 	= 1; vSignalCount.at(0)++; }
		if (vSPhosInd.size()) { evt.QCD 		= 1; vQCDCount.at(0)++; }
		if (vPhxInd.size()) { evt.Phx 		= 1; vPhxCount.at(0)++; }

		for (unsigned i=0; i < vTPhosInd.size(); ++i)
		{
			FillPhotonHists(stuple, (&Hsignal.p1j_Pho), vTPhosInd.at(i));
			FillEventHists(stuple, &(Hsignal.p1j_Evt), 1.0/vTPhosInd.size());	//weight to avoid multiple entries
		}

		for (unsigned i=0; i < vSPhosInd.size(); ++i)
		{
			FillPhotonHists(stuple, (&Hqcd.p1j_Pho), vSPhosInd.at(i));
			FillEventHists(stuple, &(Hqcd.p1j_Evt), 1.0/vSPhosInd.size());
		}

		for (unsigned i=0; i < vPhxInd.size(); ++i)
		{
			FillPhotonHists(stuple, (&Hphx.p1j_Pho), vPhxInd.at(i));
			FillEventHists(stuple, &(Hphx.p1j_Evt), 1.0/vPhxInd.size());
		}


		// now store the summary of this event
		EVT.push_back(evt);
		
}



void PhoenixPhoStudy::FillEventHists(Stuple& st, Histograms::EventHists_t* hist,
									const float fWeight)
{ 
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		exit (1);
	}
	hist->Met->Fill(st.met_Met, fWeight);
	hist->Met_Gen_d->Fill(st.met_Gen_d, fWeight);
	hist->Sumet->Fill(st.met_SumEt, fWeight);
	hist->NVertices->Fill(st.vtx_N, fWeight);
	hist->N12Vertices->Fill(st.vtx_NClass12, fWeight);
	hist->NJet15->Fill(st.jet_NJet15, fWeight);
	hist->MetRun->Fill(st.evt_RunNumber, st.met_Met, fWeight);
	hist->SumetRun->Fill(st.evt_RunNumber, st.met_SumEt, fWeight);
	hist->NVerticesRun->Fill(st.evt_RunNumber, st.vtx_N, fWeight);
	hist->N12VerticesRun->Fill(st.evt_RunNumber, st.vtx_NClass12, fWeight);
	hist->NJet15Run->Fill(st.evt_RunNumber, st.jet_NJet15, fWeight);
	hist->NJet15VsMet->Fill(st.met_Met, st.jet_NJet15, fWeight);
	hist->NJet15VsSumet->Fill(st.met_SumEt, st.jet_NJet15, fWeight);
	hist->N12VerticesVsMet->Fill(st.met_Met, st.vtx_NClass12, fWeight);
	hist->N12VerticesVsSumet->Fill(st.met_SumEt, st.vtx_NClass12, fWeight);
	if (st.tri_pho25iso == 1) hist->Triggers->Fill(0);
	if (st.tri_pho50 == 1) hist->Triggers->Fill(1);
	if (st.tri_pho70 == 1) hist->Triggers->Fill(2);

	hist->Ht->Fill(st.met_Ht, fWeight);
}

void PhoenixPhoStudy::FillPhotonHists(Stuple& st, Histograms::PhotonHists_t* hist, 
							const int iPhoIndex, const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		exit (1);
	}
	hist->Detector->Fill(st.pho_Detector[iPhoIndex], fWeight);
	hist->DetEta->Fill(st.pho_DetEta[iPhoIndex], fWeight);
	hist->DetPhi->Fill(st.pho_DetPhi[iPhoIndex], fWeight);
	hist->EtCorr->Fill(st.pho_Etc[iPhoIndex], fWeight);
	hist->XCes->Fill(st.pho_XCes[iPhoIndex], fWeight);
	hist->ZCes->Fill(st.pho_ZCes[iPhoIndex], fWeight);
	hist->HadEm->Fill(st.pho_HadEm[iPhoIndex], fWeight);
	hist->Chi2Mean->Fill(st.pho_Chi2Mean[iPhoIndex], fWeight);
	hist->N3d->Fill(st.pho_N3d[iPhoIndex], fWeight);
	hist->Iso4->Fill(st.pho_Iso4[iPhoIndex], fWeight);
	hist->TrkPt->Fill(st.pho_TrkPt[iPhoIndex], fWeight);
	hist->TrkIso->Fill(st.pho_TrkIso[iPhoIndex], fWeight);
	hist->Ces2Wire->Fill(st.pho_CesWireE2[iPhoIndex], fWeight);
	hist->Ces2Strip->Fill(st.pho_CesStripE2[iPhoIndex], fWeight);
	hist->EmTime->Fill(st.pho_EmTime[iPhoIndex], fWeight);
	hist->PhiWedge->Fill(st.pho_PhiWedge[iPhoIndex], fWeight);
	hist->EmTimeVsRun->Fill(st.evt_RunNumber, st.pho_EmTime[iPhoIndex], fWeight);
	hist->EtCorrVsRun->Fill(st.evt_RunNumber, st.pho_Etc[iPhoIndex], fWeight);
}
	

bool PhoenixPhoStudy::CheckOrder(Stuple st)
{
	if (st.pho_num >=2)
	{
		for (unsigned i=0; i < (st.pho_num-1); ++i)
			assert(st.pho_Etc[i] >= st.pho_Etc[i+1]);
	}
	if (st.ele_num >=2)
	{
		for (unsigned i=0; i < (st.ele_num-1); ++i)
			assert(st.ele_Etc[i] >= st.ele_Etc[i+1]);
	}
	if (st.jet_num >=2)
	{
		for (unsigned i=0; i < (st.jet_num-1); ++i)
			assert(st.jet_Pt[i] >= st.jet_Pt[i+1]);
	}

	return true;
}
