//////////////////////////////////////////////////////////////
// This is to get the systematic error for plots with       //
// QCD background. I used the photon sideband to predict it //
// and hence I am probing the sideband and trying to make it//
// more like the tight photon sample by tightning the loose //
// photon cuts. -- 04-21-08                                 //
//////////////////////////////////////////////////////////////


#include <iostream>
#include <sstream>
#include "PhoFlat/QCDSystematics.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
#include "PhoFlat/LoosePhoCuts.hh"
#include <iomanip>
#include "../RootTools/IOColors.hh"
#include "TROOT.h"

//-------------------------------------------------------------------
QCDSystematics::QCDSystematics():
iProgressBy(500000),
sHistFileName("QCDSystematics.root")
{
	fMinPhoEt=30;
	fMaxPhoEt=1500;
	fMaxPhoEta=1.1;
	fMinJetEt=15;
	fMaxJetEta=3.2;
	fMinMet = 0.0;
}

//--- do all my stuff here ------------------------------------------
void QCDSystematics::DoMyStuff(Stuple& stuple)
{
	if ( (stuple.tri_pho25iso != 1)  && (stuple.tri_pho50 != 1) &&
		  (stuple.tri_pho70 != 1) ) return;
		  
	if (stuple.vtx_NClass12 <= 0) return;
	if (stuple.vtx_z > 60.)  return;
	//get rid of events with extra em objects.
	//To include these events I need to double check 
	//my jet list creation. So for now it is easier 
	//not to use those events -01-14-2010
	//
	//if (stuple.ele_num>1 | stuple.pho_num>1) return;
	
	//met cut
	if ( fabs(stuple.met_Met) < fMinMet)  return;
	
	
	hCounter->Fill(1);	
	
	std::vector<int> UsedPho, UsedEle;
	
	for (unsigned int i = 0; i < stuple.pho_num; ++i ) 
	{
		if (stuple.pho_EmTime[i] < -4.8 || stuple.pho_EmTime[i] > 4.8)  continue;
		if (stuple.pho_PhoenixId[i] != 0) continue;

		HaloTypes ht(stuple.pho_Halo_seedWedge[i], 
				stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);
		if (ht.IsType(5)) return;	//discard halo events 01-13-2010
		
		if (stuple.pho_LooseId[i] == 0) {
			hCounter->Fill(2);
			if (stuple.pho_TightId[i] == 0) {
				hCounter->Fill(4);
			} else {
				hCounter->Fill(5);
			}

		}
		if (stuple.pho_TightId[i] == 0) {
			hCounter->Fill(3);
			if (stuple.pho_LooseId[i] != 0) {
				hCounter->Fill(6);
			}
		}

	}
	
	
	int ntight =0, nloose=0;
	
	for (unsigned int i = 0; i < stuple.pho_num; ++i )
	{
		//reapply all cuts just to make sure
		if (stuple.pho_Etc[i] < fMinPhoEt || stuple.pho_Etc[i] > fMaxPhoEt ) continue;
		if (fabs(stuple.pho_DetEta[i]) > fMaxPhoEta) continue;
		if (stuple.pho_TightId[i] == 0) continue;		// selectin sideband photons
		//if (stuple.pho_TightId[i] == 0) return;		// removed events with a tight photon
		if (stuple.tri_pho70 == 1 && stuple.pho_Etc[i] < 70.0) continue;
		if (stuple.tri_pho50 == 1 && stuple.pho_Etc[i] < 50.0) continue;
		if (stuple.tri_pho25iso == 1 && stuple.pho_Etc[i] < 25.0) continue;
		if (stuple.pho_LooseId[i] != 0) continue;		// selectin sideband photons
		if (stuple.pho_EmTime[i] < -4.8 || stuple.pho_EmTime[i] > 4.8) continue;  // select in-time photons
		if (stuple.pho_PhoenixId[i] != 0) continue;

		HaloTypes ht(stuple.pho_Halo_seedWedge[i], 
				stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);
		if (ht.IsType(5)) {
			if (stuple.jet_NJet15 >= 1) {	// for halo estimate
				haloPhiWedge_1j->Fill(stuple.pho_PhiWedge[i]);
				haloEmTime_1j->Fill(stuple.pho_EmTime[i]);
			}
			if (stuple.jet_NJet15 >= 2) {	// for halo estimate
				haloPhiWedge_2j->Fill(stuple.pho_PhiWedge[i]);
				haloEmTime_2j->Fill(stuple.pho_EmTime[i]);
			}
			//continue;
			return;
		}

		
		Stuple st_temp = stuple;
		LoosePhoCuts lcut(st_temp);
		lcut.SetTightenUpCut(iCutUsed);
		
		if (lcut.IDword(i) == 0) {
			//std::cout << "pass cuts " << std::endl;	
		} else {
			//std::cout << "fail cuts " << std::endl;	
			continue;
		}
		UsedPho.push_back(i);	
		NewJetList nj;
		nj.AddUnused(stuple,UsedPho,UsedEle);

		iEventsSelected++;
		HistManager myHistMan;

		const int collection = 0; //photon/jet collection
											//now there are central and +/- sigma phos and jets
											//in the Stuples 
		if (stuple.jet_NJet15 >= 1) 
		{	// mono-jet case
			iPho1JetEvts++;
			myHistMan.FillPhotonHists(stuple, collection, GA_Pho, i);
			myHistMan.FillEventHists(stuple, GA_Evt);
			myHistMan.FillJetHists(stuple, collection, GA_Jet,0);		//lead jet
			myHistMan.FillPhoton1JetHists(stuple, collection, GA_PhoJet,i,0);
		}
		
		if (stuple.jet_NJet15 >= 2) 
		{	// di-jet case
			iPho2JetEvts++;
			myHistMan.FillPhotonHists(stuple, collection, GB_Pho,i);
			myHistMan.FillEventHists(stuple, GB_Evt);
			myHistMan.FillJetHists(stuple, collection, GB_Jet1,0);		//lead jet
			myHistMan.FillJetHists(stuple, collection, GB_Jet2,1);		//2nd lead jet
			myHistMan.FillPhoton1JetHists(stuple, collection, GB_PhoJet1,i,0);
			myHistMan.FillPhoton1JetHists(stuple, collection, GB_PhoJet2,i,1);
			myHistMan.FillPhoton2JetsHists(stuple, collection, GB_PhoJets,i,0,1);
			myHistMan.FillTwoJetsHists(stuple, collection, GB_TwoJets,0,1);
		}


		break;
	} //for 


}

//-------------------------------------------------------------
void QCDSystematics::Init(TChain *tree)
{
	PrintJobSettings();

	iEventsSelected = 0;
	iPho1JetEvts = 0;
	iPho2JetEvts = 0;
	readIn = new ReadInAll(tree,&stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->Met(1);
	readIn->Init();

	//cannot add folder without creating a dir first!
	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = rootFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

} // Init

//-------------------------------------------------------------------
void QCDSystematics::CleanUp()
{
	rootFile->Write();
	rootFile->Close();
} 

//-------------------------------------------------------------------
void QCDSystematics::BookHistograms()
{

	//for halo estimate
	haloPhiWedge_1j = new TH1F("haloPhiWedge_1j","Halo(5)+>=1 Jet: #phi^{#gamma}",24,0,24);
	haloEmTime_1j = new TH1F("haloEmTime_1j","Halo(5)+>=1 Jet: #gamma EmTime",30,-10,10);
	haloPhiWedge_2j = new TH1F("haloPhiWedge_2j","Halo(5)+>=2 Jet: #phi^{#gamma}",24,0,24);
	haloEmTime_2j = new TH1F("haloEmTime_2j","Halo(5)+>=2 Jet: #gamma EmTime",30,-10,10);
	
	hCounter = new TH1F("PhoCounter","Loose/Tight Photon counter",10,0,10);
	hCounter->GetXaxis()->SetBinLabel(1,"Pass Gen. cut");
	hCounter->GetXaxis()->SetBinLabel(2,"Pass Loose");
	hCounter->GetXaxis()->SetBinLabel(3,"Pass Tight");
	hCounter->GetXaxis()->SetBinLabel(4,"Pass Loose & Tight");
	hCounter->GetXaxis()->SetBinLabel(5,"Pass Loose & Fail Tight");
	hCounter->GetXaxis()->SetBinLabel(6,"Fail Loose & Pass Tight");

	
	Histograms HistoMan;
	std::ostringstream foldername, text2title1, text2title2;
	text2title1 << "Photon+>=1Jets :";
	text2title2 << "Photon+>=2Jets :";
	

	TDirectory *sub, *subsub;

		topDir->cd();
		sub = gDirectory->mkdir("1Jet","1 Jet case");
		sub->cd();
		
		subsub = gDirectory->mkdir("Event","Event Histograms");
		subsub->cd();
				HistoMan.GetEventHistograms(GA_Evt,text2title1.str().c_str());
				
		sub->cd();
		subsub = gDirectory->mkdir("Photon","Photon Histograms");
		subsub->cd();
				HistoMan.GetPhotonHistograms(GA_Pho,text2title1.str().c_str());

		sub->cd();
		subsub = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetJetHistograms(GA_Jet,text2title1.str().c_str());
				
		sub->cd();		
		subsub = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		subsub->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(GA_PhoJet,text2title1.str().c_str());


	
		topDir->cd();
		sub = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		sub->cd();

		subsub = gDirectory->mkdir("Event","Event Histograms");
		subsub->cd();
				HistoMan.GetEventHistograms(GB_Evt,text2title2.str().c_str());

		sub->cd();
		subsub = gDirectory->mkdir("Photon","Photon Histograms");
		subsub->cd();
				HistoMan.GetPhotonHistograms(GB_Pho,text2title2.str().c_str());
	
		sub->cd();
		subsub = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetJetHistograms(GB_Jet1,text2title2.str().c_str());

		sub->cd();
		subsub = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetJetHistograms(GB_Jet2,text2title2.str().c_str());
	
		sub->cd();
		subsub = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetPhoton1JetHistograms(GB_PhoJet1,text2title2.str().c_str());
				
		sub->cd();
		subsub = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetPhoton1JetHistograms(GB_PhoJet2,text2title2.str().c_str());
					
		sub->cd();
		subsub = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		subsub->cd();
				HistoMan.GetPhoton2JetsHistograms(GB_PhoJets,text2title2.str().c_str());
				
		sub->cd();
		subsub = gDirectory->mkdir("2Jets","2Jets Histograms");
		subsub->cd();
				HistoMan.GetTwoJetsHistograms(GB_TwoJets,text2title2.str().c_str());


		
} //BookHistograms


//------ main -------------------------------------------------------
void QCDSystematics::Loop(TChain *ch, int iRunEvents)
{
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
	unsigned int iEvt2Process = 0; 	 //_______________ events to be processed
	if (iRunEvents < 0 ) {
		iRunEvents = iENTRIES;	//____ requested all entries
		iEvt2Process = iENTRIES;
	} else if (iRunEvents > 0 && iRunEvents <= iENTRIES) iEvt2Process = iRunEvents;
	else iEvt2Process = iENTRIES;		//_____ default

	std::cout << " ==> Entries found :: " << iENTRIES << std::endl;


	unsigned int iEvtProc = 0;		// number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iFirst400 =0;

	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) {
	  	readIn->GetEntry(iEvtProc);

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << iEvtProc << "\t events processed." << std::endl;
		}
		
		if (stuple.evt_RunNumber < 190851) {
			iFirst400++;
			continue;  //drop the first 400pb-1
		}
	
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;


		DoMyStuff(stuple);
	}
	
	std::cout << "======== QCDSystematics =======================" << std::endl;
	std::cout << "Entries found          = " << std::setw(8) << iENTRIES << std::endl;
	std::cout << "Entries request        = " << std::setw(8) << iRunEvents << std::endl;
	std::cout << "Events Processed       = " << std::setw(8) << iEvtProc << std::endl;
	std::cout << "Ele./Pho Events        = " << std::setw(8) << iPhoEleEvts << std::endl;
	std::cout << "Ele Events             = " << std::setw(8) << iEleEvts << std::endl;
	std::cout << "Pho Events             = " << std::setw(8) << iPhoEvts << std::endl;
	std::cout << "Other Events           = " << std::setw(8) << iOther << std::endl;
	std::cout << "First400 Events        = " << std::setw(8) << iFirst400 << std::endl;
	std::cout << "sum Events             = " << std::setw(8) << iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400 << std::endl;
	std::cout << "__________________________________________" << std::endl;

	std::cout << "Pho Evts               = " << std::setw(8) << iEventsSelected << std::endl;
	std::cout << "Pho+>=1Jet Evts        = " << std::setw(8) << iPho1JetEvts << std::endl; 
	std::cout << "Pho+>=2Jet Evts        = " << std::setw(8) << iPho2JetEvts << std::endl; 

	PrintJobSettings();
	

	std::cout << "Evts Pass General cuts = " << hCounter->GetBinContent(1) << std::endl;
	std::cout << "Evts w. Loose phos     = " << hCounter->GetBinContent(2) << std::endl;
	std::cout << "Evts w. Tight phos     = " << hCounter->GetBinContent(3) << std::endl;
	std::cout << "Evts w. pass Loose & Tight phos            = " << hCounter->GetBinContent(4) << std::endl;
	std::cout << "Evts w. pass Loose & fail Tight phos  cuts = " << hCounter->GetBinContent(5) << std::endl;
	std::cout << "Evts w. fail Loose & pass Tight phos  cuts = " << hCounter->GetBinContent(6) << std::endl;
	
	std::cout << "======= END QCDSystematics ====================" << std::endl;
	CleanUp();
} // Main


void QCDSystematics::PrintJobSettings()
{
	std::cout << red << "___________ Job Settings _________________" << std::endl;
	std::cout << "Pho Min-Max Et, Max Eta= " << std::setw(4) << fMinPhoEt << "," << std::setw(4) << fMaxPhoEt << "," << std::setw(4) << fMaxPhoEta << std::endl;
	std::cout << "Jet Min Et, Max Eta    = " << std::setw(4) <<fMinJetEt  << "," << std::setw(4) << fMaxJetEta << std::endl;
	std::cout << "Min MEt                = " << std::setw(4) << fMinMet << std::endl;
	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	std::string str;
	if (iCutUsed == 1) str = "HADEM";
	if (iCutUsed == 2) str = "ISO";
	if (iCutUsed == 3) str = "TRKPT";
	if (iCutUsed == 4) str = "TRKISO";
	std::cout << "Tighten up cut         = " << str << std::endl;
	std::cout << "__________________________________________" << clearatt << std::endl;


}
