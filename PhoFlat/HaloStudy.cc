#include <iostream>
#include <sstream>
#include "PhoFlat/HaloStudy.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TDirectory.h"
//-------------------------------------------------------------------
HaloStudy::HaloStudy():
iProgressBy(1000),
sHistFileName("HaloStudy.root")
{
}

//--- do all my stuff here ------------------------------------------
void HaloStudy::DoMyStuff(Stuple& stuple)
{
	for (unsigned int i = 0; i < stuple.pho_num; ++i ) {
		if (stuple.pho_TightId[i] != 0) continue;		// tight photons
		if (stuple.pho_Etc[i] < 30.) continue;
		if (fabs(stuple.pho_DetEta[i]) >1.1) continue;
		if (stuple.pho_EmTime[i] < -4.8 || stuple.pho_EmTime[i] > 4.8) continue;  // select in-time photons

		HaloTypes ht(stuple.pho_Halo_seedWedge[i], 
				stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);

		for (int j =0; j < ht.iNHALOTYPES; ++j) {
			if (ht.IsType(j)) {
				vHalosFound[j]++;
				++iEventsSelected;
				float fMetEtRatio = stuple.met_Met / stuple.pho_Etc[i];
				//std::cout << "type ="<< j << "\tratio=" << fMetEtRatio << std::endl;
				//std::cout << "\tjets = " << stuple.jet_NJet15 << std::endl;
			
				TLorentzVector tlPho(stuple.pho_Px[i], stuple.pho_Py[i], stuple.pho_Pz[i], stuple.pho_E[i]);
				TLorentzVector tlJet(stuple.jet_Px[0], stuple.jet_Py[0], stuple.jet_Pz[0], stuple.jet_E[0]);
							
				float fDelR = tlJet.DeltaR(tlPho);
							
				//if (fMetEtRatio >= 1.0) {
				if (fDelR > 0.2) {
					if (stuple.jet_NJet15 >= 1) {	// mono-jet case
						myHistMan.FillPhotonHists(stuple,0, GA_Pho[j],i);
						myHistMan.FillEventHists(stuple, GA_Evt[j]);
						myHistMan.FillJetHists(stuple,0, GA_Jet[j],0);		//lead jet
						myHistMan.FillPhoton1JetHists(stuple,0, GA_PhoJet[j],i,0);
					}
					
					if (stuple.jet_NJet15 >= 2) {	// di-jet case
						myHistMan.FillPhotonHists(stuple, 0, GB_Pho[j],i);
						myHistMan.FillEventHists(stuple, GB_Evt[j]);
						myHistMan.FillJetHists(stuple, 0, GB_Jet1[j], 0);		//lead jet
						myHistMan.FillJetHists(stuple, 0, GB_Jet2[j], 1);		//2nd lead jet
						myHistMan.FillPhoton1JetHists(stuple,0, GB_PhoJet1[j],i, 0);
						myHistMan.FillPhoton1JetHists(stuple,0, GB_PhoJet2[j],i, 1);
						myHistMan.FillPhoton2JetsHists(stuple,0, GB_PhoJets[j],i, 0,1);
						myHistMan.FillTwoJetsHists(stuple,0, GB_TwoJets[j], 0, 1);
					}

				} else {

					if (stuple.jet_NJet15 >= 1) {	// mono-jet case
						myHistMan.FillPhotonHists(stuple,0, LA_Pho[j],1);
						myHistMan.FillEventHists(stuple, LA_Evt[j]);
						myHistMan.FillJetHists(stuple,0, LA_Jet[j],0);		//lead jet
						myHistMan.FillPhoton1JetHists(stuple,0, LA_PhoJet[j],i,0);
					}
					
					if (stuple.jet_NJet15 >= 2) {	// di-jet case
						myHistMan.FillPhotonHists(stuple,0, LB_Pho[j],i);
						myHistMan.FillEventHists(stuple, LB_Evt[j]);
						myHistMan.FillJetHists(stuple,0, LB_Jet1[j],0);		//lead jet
						myHistMan.FillJetHists(stuple,0, LB_Jet2[j],1);		//2nd lead jet
						myHistMan.FillPhoton1JetHists(stuple,0, LB_PhoJet1[j],i,0);
						myHistMan.FillPhoton1JetHists(stuple,0, LB_PhoJet2[j],i,1);
						myHistMan.FillPhoton2JetsHists(stuple,0, LB_PhoJets[j],i,0,1);
						myHistMan.FillTwoJetsHists(stuple,0, LB_TwoJets[j],0,1);
					}
				}

			}
		} //for
		break;
	} //for 


}

//-------------------------------------------------------------
void HaloStudy::Init(TChain *tree)
{
	HaloTypes ht;
	for (int i =0; i <ht.iNHALOTYPES; i++) vHalosFound.push_back(0);
	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	//cannot add folder without creating a dir first!
	topDir = gDirectory->mkdir("Hist","All Histograms");	//create top level dir.

	iEventsSelected = 0;
	readIn = new ReadInAll(tree,&stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->Met(1);
	readIn->Init();
} // Init

//-------------------------------------------------------------------
void HaloStudy::CleanUp()
{
	//vHalosFound.clear();
	HaloTypes ht;
	for (int i =0; i <ht.iNHALOTYPES; i++) vHalosFound.push_back(0);
	rootFile->Write();
	rootFile->Close();
} 

//-------------------------------------------------------------------
void HaloStudy::BookHistograms()
{
	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

	Histograms::EventHists_t hEvt;
	Histograms::PhotonHists_t hPho;
	Histograms::JetHists_t hJet;
	Histograms::Photon1JetHists_t hPhoJet;
	Histograms::Photon2JetsHists_t hPhoJets;
	Histograms::TwoJetsHists_t hTwoJets;

	Histograms HistoMan;
	// for MEtEtRatio < 1  : halo type 0-5
	//folder structure
	// HaloType/MetEtRatio_gt1/1Jet/...

	HaloTypes ht;
	for (int i =0; i <ht.iNHALOTYPES; i++) {
		
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

		LA_Evt.push_back(hEvt);
		LB_Evt.push_back(hEvt);
		LA_Pho.push_back(hPho);
		LB_Pho.push_back(hPho);
		LA_Jet.push_back(hJet);
		LB_Jet1.push_back(hJet);
		LB_Jet2.push_back(hJet);
		LA_PhoJet.push_back(hPhoJet);
		LB_PhoJet1.push_back(hPhoJet);
		LB_PhoJet2.push_back(hPhoJet);
		LB_PhoJets.push_back(hPhoJets);
		LB_TwoJets.push_back(hTwoJets);


		std::ostringstream foldername, text2title1G, text2title1L,text2title2G,text2title2L;
		foldername << "Type"<< i;
		text2title1G << "Halo("<< i << ")+>=1Jets : #Delta R^{#gamma,Lead jet}>0.2:";
		text2title1L << "Halo("<< i << ")+>=1Jets : #Delta R^{#gamma,Lead jet}<0.2:";
		text2title2G << "Halo("<< i << ")+>=2Jets : #Delta R^{#gamma,Lead jet}>0.2:";
		text2title2L << "Halo("<< i << ")+>=2Jets : #Delta R^{#gamma,Lead jet}<0.2:";
			
		halofolder.push_back(new TFolder(foldername.str().c_str(),foldername.str().c_str()));
		topDir->Add(halofolder[i]);
	
		MetEtRatioHist.push_back(new TH1F("MetEtRatio",foldername.str().c_str(),100,0,2));
			MetEtRatioHist[i]->GetXaxis()->SetTitle("ME_{T}/E_{T}^{#gamma}");
			halofolder[i]->Add(MetEtRatioHist[i]);

		TDirectory *folder = halofolder[i]->AddFolder("MetEtRatio_gt1","Met/Photon Et Ratio > 1");
		TDirectory *new_folder = folder->AddFolder("1Jet","1 Jet case");
		TDirectory *sub_folder = new_folder->AddFolder("Event","Event Histograms");
		
				HistoMan.GetEventHistograms(GA_Evt[i],text2title1G.str().c_str());
				sub_folder = new_folder->AddFolder("Photon",text2title1G.str().c_str());
				HistoMan.GetPhotonHistograms(GA_Pho[i],text2title1G.str().c_str());

				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(GA_Jet[i],text2title1G.str().c_str());
				sub_folder = new_folder->AddFolder("PhotonLeadJet",text2title1G.str().c_str());
				HistoMan.GetPhoton1JetHistograms(GA_PhoJet[i],text2title1G.str().c_str());


			new_folder = folder->AddFolder("2Jet","#gamma + >=2 Jet case");
	
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
				HistoMan.GetEventHistograms(GB_Evt[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
				HistoMan.GetPhotonHistograms(GB_Pho[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(GB_Jet1[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
				HistoMan.GetJetHistograms(GB_Jet2[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(GB_PhoJet1[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(GB_PhoJet2[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
				HistoMan.GetPhoton2JetsHistograms(GB_PhoJets[i],text2title2G.str().c_str());

				sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
				HistoMan.GetTwoJetsHistograms(GB_TwoJets[i],text2title2G.str().c_str());


		folder = halofolder[i]->AddFolder("MetEtRatio_lt1","Met/Photon Et Ratio < 1");
			new_folder = folder->AddFolder("1Jet","1 Jet case");
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
		
				HistoMan.GetEventHistograms(LA_Evt[i],text2title1L.str().c_str());
				sub_folder = new_folder->AddFolder("Photon",text2title1L.str().c_str());
				HistoMan.GetPhotonHistograms(LA_Pho[i],text2title1L.str().c_str());

				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(LA_Jet[i],text2title1L.str().c_str());
				sub_folder = new_folder->AddFolder("PhotonLeadJet",text2title1L.str().c_str());
				HistoMan.GetPhoton1JetHistograms(LA_PhoJet[i],text2title1L.str().c_str());


			new_folder = folder->AddFolder("2Jet","#gamma + >=2 Jet case");
	
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
				HistoMan.GetEventHistograms(LB_Evt[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
				HistoMan.GetPhotonHistograms(LB_Pho[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(LB_Jet1[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
				HistoMan.GetJetHistograms(LB_Jet2[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(LB_PhoJet1[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(LB_PhoJet2[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
				HistoMan.GetPhoton2JetsHistograms(LB_PhoJets[i],text2title2L.str().c_str());

				sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
				HistoMan.GetTwoJetsHistograms(LB_TwoJets[i],text2title2L.str().c_str());
		}


} //BookHistograms

//------ main -------------------------------------------------------
void HaloStudy::Main(TChain *ch, int iRunEvents)
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
			
	  	if (stuple.pho_num < 1) continue;
	  	if (stuple.ele_num > 0) continue;
		//want to see with and without electrons. remember to put back the electron to jet when needed!
		// lets see by removing e event as ele can become my jet!
		//if (stuple.ele_num > 0) {
			//std::cout << "found ele" << std::endl;
			//NewJetList(stuple,"ele");
		//}
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
	for (int i =0; i < ht.iNHALOTYPES; ++i) std::cout << "Halos Found: type " << i << "\t= "  << vHalosFound[i] << std::endl;
	std::cout << "===========================================" << entries << std::endl;
	CleanUp();
} // Main

