#include <iostream>
#include <sstream>
#include "PhoFlat/SidebandHalos.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
//-------------------------------------------------------------------
SidebandHalos::SidebandHalos():
iProgressBy(1000),
sHistFileName("SidebandHalos.root")
{
}

//--- do all my stuff here ------------------------------------------
void SidebandHalos::DoMyStuff(Stuple& stuple)
{
	for (unsigned int i = 0; i < stuple.pho_num; ++i ) {
		//select sideband photons
		if (stuple.pho_LooseId[i]  != 0) continue;
		if (stuple.pho_TightId[i] == 0) continue;
		if (stuple.pho_EmTime[i] < -4.8 || stuple.pho_EmTime[i] > 4.8) continue;
		//if (stuple.ele_num > 0) NewJetList(stuple,1); // add electrons to jet list

		HaloTypes ht(stuple.pho_Halo_seedWedge[i], 
				stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);

		for (int j =0; j < ht.iHALOTYPES; ++j) {
			if (ht.IsType(j)) {
				vHalosFound[j]++;
				++iEventsSelected;
			
				if (stuple.jet_NJet15 >= 1) {	// mono-jet case
					myHistMan.FillPhotonHists(stuple, GA_Pho[j],i);
					myHistMan.FillEventHists(stuple, GA_Evt[j]);
					myHistMan.FillJetHists(stuple, GA_Jet[j],0);		//lead jet
					myHistMan.FillPhoton1JetHists(stuple, GA_PhoJet[j],i,0);
				}
				
				if (stuple.jet_NJet15 >= 2) {	// di-jet case
					myHistMan.FillPhotonHists(stuple, GB_Pho[j],i);
					myHistMan.FillEventHists(stuple, GB_Evt[j]);
					myHistMan.FillJetHists(stuple, GB_Jet1[j],0);		//lead jet
					myHistMan.FillJetHists(stuple, GB_Jet2[j],1);		//2nd lead jet
					myHistMan.FillPhoton1JetHists(stuple, GB_PhoJet1[j],i,0);
					myHistMan.FillPhoton1JetHists(stuple, GB_PhoJet2[j],i,1);
					myHistMan.FillPhoton2JetsHists(stuple, GB_PhoJets[j],i,0,1);
					myHistMan.FillTwoJetsHists(stuple, GB_TwoJets[j],0,1);
				}
			}
		} //for
		break;
	} //for 


}

//-------------------------------------------------------------
void SidebandHalos::Init(TChain *tree)
{
	HaloTypes ht;
	for (int i =0; i <ht.iHALOTYPES; i++) vHalosFound.push_back(0);
	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	//cannot add folder without creating a dir first!
	topDir = rootFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();

	iEventsSelected = 0;
	readIn = new ReadInAll(tree,&stuple);
} // Init

//-------------------------------------------------------------------
void SidebandHalos::CleanUp()
{
	//vHalosFound.clear();
	HaloTypes ht;
	for (int i =0; i <ht.iHALOTYPES; i++) vHalosFound.push_back(0);
	rootFile->Write();
	rootFile->Close();
} 

//-------------------------------------------------------------------
void SidebandHalos::BookHistograms()
{
	Histograms::EventHists_t hEvt;
	Histograms::PhotonHists_t hPho;
	Histograms::JetHists_t hJet;
	Histograms::Photon1JetHists_t hPhoJet;
	Histograms::Photon2JetsHists_t hPhoJets;
	Histograms::TwoJetsHists_t hTwoJets;

	Histograms HistoMan;

	HaloTypes ht;
	TDirectory *top, *sub, *subsub;

	for (int i =0; i <ht.iHALOTYPES; i++) {
		
		GA_Evt.push_back(hEvt);
		GB_Evt.push_back(hEvt);
		GA_Pho.push_back(hPho);
		GB_Pho.push_back(hPho);
		GA_Jet.push_back(hJet);
		GB_Jet1.push_back(hJet);
		GB_Jet2.push_back(hJet);
		GA_PhoJet.push_back(hPhoJet);
		GB_PhoJet1.push_back(hPhoJet);
		GB_PhoJet2.push_back(hPhoJet);
		GB_PhoJets.push_back(hPhoJets);
		GB_TwoJets.push_back(hTwoJets);

		std::ostringstream foldername, text2title1, text2title2;
		foldername << "Type"<< i;
		text2title1 << "Sideband Halo("<< i << ")+>=1Jets :";
		text2title2 << "Sideband Halo("<< i << ")+>=2Jets :";
			
		topDir->cd();
		top = topDir->mkdir(foldername.str().c_str(), foldername.str().c_str());
		top->cd();
	
		sub = gDirectory->mkdir("1Jet","1 Jet case");
		sub->cd();
		
		subsub = gDirectory->mkdir("Event","Event Histograms");
		subsub->cd();
				HistoMan.GetEventHistograms(GA_Evt[i],text2title1.str().c_str());
				
		sub->cd();
		subsub = gDirectory->mkdir("Photon","Photon Histograms");
		subsub->cd();
				HistoMan.GetPhotonHistograms(GA_Pho[i],text2title1.str().c_str());

		sub->cd();
		subsub = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetJetHistograms(GA_Jet[i],text2title1.str().c_str());
				
		sub->cd();		
		subsub = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		subsub->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(GA_PhoJet[i],text2title1.str().c_str());


	
		top->cd();
		sub = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		sub->cd();

		subsub = gDirectory->mkdir("Event","Event Histograms");
		subsub->cd();
				HistoMan.GetEventHistograms(GB_Evt[i],text2title2.str().c_str());

		sub->cd();
		subsub = gDirectory->mkdir("Photon","Photon Histograms");
		subsub->cd();
				HistoMan.GetPhotonHistograms(GB_Pho[i],text2title2.str().c_str());
	
		sub->cd();
		subsub = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetJetHistograms(GB_Jet1[i],text2title2.str().c_str());

		sub->cd();
		subsub = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetJetHistograms(GB_Jet2[i],text2title2.str().c_str());
	
		sub->cd();
		subsub = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetPhoton1JetHistograms(GB_PhoJet1[i],text2title2.str().c_str());
				
		sub->cd();
		subsub = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetPhoton1JetHistograms(GB_PhoJet2[i],text2title2.str().c_str());
					
		sub->cd();
		subsub = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		subsub->cd();
				HistoMan.GetPhoton2JetsHistograms(GB_PhoJets[i],text2title2.str().c_str());
				
		sub->cd();
		subsub = gDirectory->mkdir("2Jets","2Jets Histograms");
		subsub->cd();
				HistoMan.GetTwoJetsHistograms(GB_TwoJets[i],text2title2.str().c_str());

	}


} //BookHistograms

//------ main -------------------------------------------------------
bool SidebandHalos::PrePass(const Stuple& stuple) {

	//make event preselections here
	bool bTrig = false, bVtx = false, bVtx_z = false;

	if ( (stuple.tri_pho25iso == 1) || (stuple.tri_pho50 == 1) || 
		  (stuple.tri_pho70 == 1) )
		bTrig = true;

	if (stuple.vtx_NClass12 >=1) bVtx = true;
	if (stuple.vtx_z < 60.)  bVtx_z = true;

	return (bTrig & bVtx & bVtx_z);

}


//------ main -------------------------------------------------------
void SidebandHalos::Main(TChain *ch, int iRunEvents)
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
	
	unsigned int entries = (unsigned int) myChain->GetEntries();

	if (iRunEvents > 0 && iRunEvents < entries) entries = iRunEvents;
	else std::cout << __FILE__ << "::" << __LINE__ << "::" << "Only " << entries << " found and those will be processed." << std::endl;

	if (iRunEvents < 0 ) iRunEvents = entries;	//requested all entries

	unsigned int iEvtProc = 0;		// number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iFirst400 =0;

	for ( ; iEvtProc < entries; ++iEvtProc) {
	  	readIn->GetEntry(iEvtProc);
		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) std::cout << iEvtProc << "\t events processed." << std::endl;
		if (stuple.evt_RunNumber < 190851) {
			iFirst400++;
			continue;  //drop the first 400pb-1
		}
			
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;
			
		if (! PrePass(stuple)) continue;
		DoMyStuff(stuple);
	};
	
	std::cout << "===========================================" << entries << std::endl;
	std::cout << "Entries found    = " << entries << std::endl;
	std::cout << "Entries request  = " << iRunEvents << std::endl;
	std::cout << "Events Processed = " << iEvtProc << std::endl;
	std::cout << "Ele./Pho Events  = " << iPhoEleEvts << std::endl;
	std::cout << "Ele Events       = " << iEleEvts << std::endl;
	std::cout << "Pho Events       = " << iPhoEvts << std::endl;
	std::cout << "Other Events     = " << iOther << std::endl;
	std::cout << "First400 Events  = " << iFirst400 << std::endl;
	std::cout << "sum Events       = " << iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400 << std::endl;
	
	std::cout << "Events Selected  = " << iEventsSelected << std::endl;
	std::cout << "Data Written to  = " << GetHistFileName() << std::endl;
	HaloTypes ht;
	for (int i =0; i < ht.iHALOTYPES; ++i) std::cout << "Halos Found: type " << i << "\t= "  << vHalosFound[i] << std::endl;
	std::cout << "===========================================" << entries << std::endl;
	CleanUp();
} // Main

