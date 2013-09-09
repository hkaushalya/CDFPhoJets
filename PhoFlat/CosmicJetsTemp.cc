#include <iostream>
#include <sstream>
#include "PhoFlat/CosmicJetsTemp.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
//-------------------------------------------------------------------
CosmicJetsTemp::CosmicJetsTemp():
iProgressBy(500000),
sHistFileName("CosmicJetsTemp.root")
{
}

//--- do all my stuff here ------------------------------------------
void CosmicJetsTemp::DoMyStuff(Stuple& stuple)
{
	if ( (stuple.tri_pho25iso != 1)  && (stuple.tri_pho50 != 1) &&
		  (stuple.tri_pho70 != 1) ) return;
		  
	if (stuple.vtx_NClass12 <= 0) return;
	if (stuple.vtx_z > 60.)  return;
	if (stuple.pho_Ntight < 1) return;
	
	std::vector<int> UsedPho, UsedEle;
	

	for (unsigned int i = 0; i < stuple.pho_num; ++i ) {
		if (stuple.pho_TightId[i] != 0) continue;		// tight photons
		if (stuple.pho_EmTime[i] < 30 || stuple.pho_EmTime[i] > 90) continue;  // select in-time photons
		if (stuple.pho_PhoenixId[i] != 0) continue;  //____________________________ reject phoenix photons
		HaloTypes ht(stuple.pho_Halo_seedWedge[i], 
				stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);

		if (ht.IsType(5))	continue;



		iEventsSelected++;
		UsedPho.push_back(i);	
		NewJetList nj;
		nj.AddUnused(stuple,UsedPho,UsedEle);

		// now look for jets>15 as i did not make a jet cut when skimming
		std::vector<int> vJet15Ind;
		TLorentzVector tlPho(stuple.pho_Px[i], stuple.pho_Py[i], stuple.pho_Pz[i], stuple.pho_E[i]);
		for (int j=0; j < stuple.jet_num; j++) {
			TLorentzVector tlJet(stuple.jet_Px[j], stuple.jet_Py[j], stuple.jet_Pz[j], stuple.jet_E[j]);
			float DelR = tlPho.DeltaR(tlJet);
			if (DelR<0.2) continue;			// _________________________this is a temp fix to em objects not being removed properly
													// _________________________from jet list. few events for photons. but a lot for eles  -- 04-04-08
			if (stuple.jet_Pt[j] < 15.0) continue;
			if (fabs(tlJet.Eta()) <3.0) vJet15Ind.push_back(j);		//Sasha cuts on eta only in MET clean up cuts 04-04-08
		}


		//now set the NJet15 correctly to fix the plots
		stuple.jet_NJet15 = (UInt_t) vJet15Ind.size();

		
	//if (stuple.jet_num != vJet15Ind.size()) stuple.Dump();  need to check this. I am saving all jets.
	// But I don't find any jets <15GeV in the skim??? -- 04-04-08

		const int collection = 0; //central photon/jet collections
		HistManager myHistMan;
		if (vJet15Ind.size() >= 1) {	// mono-jet case
			iCosmic1JetEvts++;
			myHistMan.FillPhotonHists(stuple, collection, GA_Pho,i);
			myHistMan.FillEventHists(stuple, GA_Evt);
			myHistMan.FillJetHists(stuple, collection, GA_Jet,vJet15Ind[0]);		//lead jet
			myHistMan.FillPhoton1JetHists(stuple, collection, GA_PhoJet,i,vJet15Ind[0]);
		}
		
		if (vJet15Ind.size() >= 2) {		// di-jet case
			iCosmic2JetEvts++;
			myHistMan.FillPhotonHists(stuple, collection, GB_Pho,i);
			myHistMan.FillEventHists(stuple, GB_Evt);
			myHistMan.FillJetHists(stuple, collection, GB_Jet1,vJet15Ind[0]);		//lead jet
			myHistMan.FillJetHists(stuple, collection, GB_Jet2,vJet15Ind[1]);		//2nd lead jet
			myHistMan.FillPhoton1JetHists(stuple, collection, GB_PhoJet1,i,vJet15Ind[0]);
			myHistMan.FillPhoton1JetHists(stuple, collection, GB_PhoJet2,i,vJet15Ind[0]);
			myHistMan.FillPhoton2JetsHists(stuple, collection, GB_PhoJets,i,vJet15Ind[0],vJet15Ind[1]);
			myHistMan.FillTwoJetsHists(stuple, collection, GB_TwoJets,vJet15Ind[0],vJet15Ind[1]);
		}



		break;
	} //for 

}

//-------------------------------------------------------------
void CosmicJetsTemp::Init(TChain *tree)
{
	iEventsSelected = 0;
	iCosmic1JetEvts = 0;
	iCosmic2JetEvts = 0;
	iNCosLow1 = 0;
	iNCosUp1 = 0;
	iNCosLow2 = 0;
	iNCosUp2 = 0;
	readIn = new ReadInAll(tree,&stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->Met(1);

	//cannot add folder without creating a dir first!
	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = rootFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

} // Init

//-------------------------------------------------------------------
void CosmicJetsTemp::CleanUp()
{
	rootFile->Write();
	rootFile->Close();
} 

//-------------------------------------------------------------------
void CosmicJetsTemp::BookHistograms()
{
	Histograms HistoMan;
	std::ostringstream foldername, text2title1, text2title2;
	text2title1 << "Cosmic Photon+>=1Jets :";
	text2title2 << "Cosmic Photon+>=2Jets :";
	

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
void CosmicJetsTemp::Main(TChain *ch, int iRunEvents)
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
	};
	
	std::cout << "======= CosmicJets ========================" << std::endl;
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
	
	std::cout << "Events Selected  = " << iEventsSelected << std::endl;
	std::cout << "Cosmic Pho+>=1Jet= " << iCosmic1JetEvts << std::endl; 
	std::cout << "Cosmic Pho+>=2Jet= " << iCosmic2JetEvts << std::endl; 
	std::cout << "Data Written to  = " << GetHistFileName() << std::endl;
	std::cout << "__________________________________________" << std::endl;
	std::cout << "Cosmic+>=1j >50ns & <70ns = " << iNCosLow1 << std::endl; 
	std::cout << "Cosmic+>=1j >70ns & <90ns = " << iNCosUp1 << std::endl; 
	std::cout << "Cosmic+>=2j >50ns & <70ns = " << iNCosLow2 << std::endl; 
	std::cout << "Cosmic+>=2j >70ns & <90ns = " << iNCosUp2 << std::endl; 
	std::cout << "======= END CosmicJets ====================" << std::endl;
	CleanUp();
} // Main

