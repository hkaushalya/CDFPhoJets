#include <iostream>
#include <sstream>
#include "PhoFlat/HaloStudy_DelR.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
//-------------------------------------------------------------------
HaloStudy_DelR::HaloStudy_DelR():
iProgressBy(50000000),
sHistFileName("HaloStudy_DelR.root")
{
}

//--- do all my stuff here ------------------------------------------
void HaloStudy_DelR::DoMyStuff(Stuple& stuple)
{
	for (unsigned int i = 0; i < stuple.pho_num; ++i ) {
		if (stuple.pho_TightId[i] != 0) continue;		// tight photons
		if (stuple.pho_EmTime[i] < -4.8 || stuple.pho_EmTime[i] > 4.8) continue;  // select in-time photons
		TLorentzVector tlPho(stuple.pho_Px[i], stuple.pho_Py[i], stuple.pho_Pz[i], stuple.pho_E[i]);
		TLorentzVector tlJet(stuple.jet_Px[0], stuple.jet_Py[0], stuple.jet_Pz[0], stuple.jet_E[0]);
							
		float fDelR = tlJet.DeltaR(tlPho);
		if (fDelR < 0.4) {
			std::cout << stuple.evt_RunNumber << "\t" << stuple.evt_EventNumber <<"\t\t"
						 << stuple.pho_PhiWedge[i] << "\t" << stuple.met_Met << "\t"
						 << stuple.met_SumEt << "\t" << stuple.jet_num << std::endl;

		}
		if (fDelR > 0.4) {
			HaloTypes ht(stuple.pho_Halo_seedWedge[i], stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);
					
			for (int j =0; j < ht.iHALOTYPES; ++j) {
				if (ht.IsType(j)) {
					vHalosFound_GT[j]++;
					DelRHist_GT[j]->Fill(fDelR);
							
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

			}

		} else {
			HaloTypes ht(stuple.pho_Halo_seedWedge[i], stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);
					
			for (int j =0; j < ht.iHALOTYPES; ++j) {
				if (ht.IsType(j)) {
					vHalosFound_LT[j]++;
					DelRHist_LT[j]->Fill(fDelR);
	
					if (stuple.jet_NJet15 >= 1) {	// mono-jet case
						myHistMan.FillPhotonHists(stuple, LA_Pho[j],i);
						myHistMan.FillEventHists(stuple, LA_Evt[j]);
						myHistMan.FillJetHists(stuple, LA_Jet[j],0);		//lead jet
						myHistMan.FillPhoton1JetHists(stuple, LA_PhoJet[j],i,0);
					}
					
					if (stuple.jet_NJet15 >= 2) {	// di-jet case
						myHistMan.FillPhotonHists(stuple, LB_Pho[j],i);
						myHistMan.FillEventHists(stuple, LB_Evt[j]);
						myHistMan.FillJetHists(stuple, LB_Jet1[j],0);		//lead jet
						myHistMan.FillJetHists(stuple, LB_Jet2[j],1);		//2nd lead jet
						myHistMan.FillPhoton1JetHists(stuple, LB_PhoJet1[j],i,0);
						myHistMan.FillPhoton1JetHists(stuple, LB_PhoJet2[j],i,1);
						myHistMan.FillPhoton2JetsHists(stuple, LB_PhoJets[j],i,0,1);
						myHistMan.FillTwoJetsHists(stuple, LB_TwoJets[j],0,1);
					}

				}
			}
			
		} // if
		break;
	} //for 


}

//-------------------------------------------------------------
void HaloStudy_DelR::Init(TChain *tree)
{
	HaloTypes ht;
	for (int i =0; i <ht.iHALOTYPES; i++) {
		vHalosFound_GT.push_back(0);
		vHalosFound_LT.push_back(0);
	}
	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	//cannot add folder without creating a dir first!
	topDir = gDirectory->mkdir("Hist","All Histograms");	//create top level dir.

	iEventsSelected = 0;
	readIn = new ReadInAll(tree,&stuple);
} // Init

//-------------------------------------------------------------------
void HaloStudy_DelR::CleanUp()
{
	HaloTypes ht;
	for (int i =0; i <ht.iHALOTYPES; i++) {
		vHalosFound_GT.push_back(0);
		vHalosFound_LT.push_back(0);
	}
	rootFile->Write();
	rootFile->Close();
} 

//-------------------------------------------------------------------
void HaloStudy_DelR::BookHistograms()
{
	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

	Histograms::EventHists_t hEvt;
	Histograms::PhotonHists_t hPho;
	Histograms::JetHists_t hJet;
	Histograms::Photon1JetHists_t hPhoJet;
	Histograms::Photon2JetsHists_t hPhoJets;
	Histograms::TwoJetsHists_t hTwoJets;

	Histograms HistoMan;

	HaloTypes ht;

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
		text2title1G << "Halo("<< i << ")+>=1Jets : #Delta R^{#gamma,Lead jet}>0.4:";
		text2title1L << "Halo("<< i << ")+>=1Jets : #Delta R^{#gamma,Lead jet}<0.4:";
		text2title2G << "Halo("<< i << ")+>=2Jets : #Delta R^{#gamma,Lead jet}>0.4:";
		text2title2L << "Halo("<< i << ")+>=2Jets : #Delta R^{#gamma,Lead jet}<0.4:";
			
		halofolder.push_back(new TFolder(foldername.str().c_str(),foldername.str().c_str()));
		topDir->Add(halofolder[i]);

	
		TFolder *folder = halofolder[i]->AddFolder("DelR_gt04","DelR>0.4");

			DelRHist_GT.push_back(new TH1F("DelR",foldername.str().c_str(),100,0,2));
				DelRHist_GT[i]->GetXaxis()->SetTitle("#Delta R^{#gamma,Lead Jet}");
				folder->Add(DelRHist_GT[i]);

		TFolder *new_folder = folder->AddFolder("1Jet","1 Jet case");
		TFolder *sub_folder = new_folder->AddFolder("Event","Event Histograms");
		
				HistoMan.GetEventHistograms(sub_folder,GA_Evt[i],text2title1G.str().c_str());
				sub_folder = new_folder->AddFolder("Photon",text2title1G.str().c_str());
				HistoMan.GetPhotonHistograms(sub_folder,GA_Pho[i],text2title1G.str().c_str());

				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,GA_Jet[i],text2title1G.str().c_str());
				sub_folder = new_folder->AddFolder("PhotonLeadJet",text2title1G.str().c_str());
				HistoMan.GetPhoton1JetHistograms(sub_folder,GA_PhoJet[i],text2title1G.str().c_str());


			new_folder = folder->AddFolder("2Jet","#gamma + >=2 Jet case");
	
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
				HistoMan.GetEventHistograms(sub_folder,GB_Evt[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
				HistoMan.GetPhotonHistograms(sub_folder,GB_Pho[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,GB_Jet1[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,GB_Jet2[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(sub_folder,GB_PhoJet1[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(sub_folder,GB_PhoJet2[i],text2title2G.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
				HistoMan.GetPhoton2JetsHistograms(sub_folder,GB_PhoJets[i],text2title2G.str().c_str());

				sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
				HistoMan.GetTwoJetsHistograms(sub_folder,GB_TwoJets[i],text2title2G.str().c_str());


		folder = halofolder[i]->AddFolder("DelR_lt04","DelR<0.4");

			DelRHist_LT.push_back(new TH1F("DelR",foldername.str().c_str(),100,0,2));
				DelRHist_LT[i]->GetXaxis()->SetTitle("#Delta R^{#gamma,Lead Jet}");
				folder->Add(DelRHist_LT[i]);
		
			new_folder = folder->AddFolder("1Jet","1 Jet case");
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
		
		
				HistoMan.GetEventHistograms(sub_folder,LA_Evt[i],text2title1L.str().c_str());
				sub_folder = new_folder->AddFolder("Photon",text2title1L.str().c_str());
				HistoMan.GetPhotonHistograms(sub_folder,LA_Pho[i],text2title1L.str().c_str());

				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,LA_Jet[i],text2title1L.str().c_str());
				sub_folder = new_folder->AddFolder("PhotonLeadJet",text2title1L.str().c_str());
				HistoMan.GetPhoton1JetHistograms(sub_folder,LA_PhoJet[i],text2title1L.str().c_str());


			new_folder = folder->AddFolder("2Jet","#gamma + >=2 Jet case");
	
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
				HistoMan.GetEventHistograms(sub_folder,LB_Evt[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
				HistoMan.GetPhotonHistograms(sub_folder,LB_Pho[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,LB_Jet1[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,LB_Jet2[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(sub_folder,LB_PhoJet1[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(sub_folder,LB_PhoJet2[i],text2title2L.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
				HistoMan.GetPhoton2JetsHistograms(sub_folder,LB_PhoJets[i],text2title2L.str().c_str());

				sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
				HistoMan.GetTwoJetsHistograms(sub_folder,LB_TwoJets[i],text2title2L.str().c_str());
		}


} //BookHistograms

//------ main -------------------------------------------------------
void HaloStudy_DelR::Main(TChain *ch, int iRunEvents)
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
	unsigned int iEntries = iENTRIES;
	if (iRunEvents > 0 && iRunEvents <= iENTRIES) iEntries = iRunEvents;
	else std::cout << __FILE__ << "::" << __LINE__ << "::" << "Only " << iEntries << " found and those will be processed." << std::endl;

	if (iRunEvents < 0 ) iRunEvents = iENTRIES;	//requested all entries

	unsigned int iEvtProc = 0;		// number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iEle1JetEvts = 0, iPho1JetEvts = 0, iPhoEle1JetEvts =0;
	unsigned int iEle2JetEvts = 0, iPho2JetEvts = 0, iPhoEle2JetEvts =0;
	unsigned int iFirst400 =0;
	std::cout << "Run\tEvt\tPhiWedge\tMETcor\tNJet15"  << std::endl;
	

	for ( ; iEvtProc < iEntries; ++iEvtProc) {
	  	readIn->GetEntry(iEvtProc);
		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) std::cout << iEvtProc << "\t events processed." << std::endl;
		if (stuple.evt_RunNumber < 190851) {
			iFirst400++;
			continue;  //drop the first 400pb-1
		}
			
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) {
			iPhoEleEvts++;
			if (stuple.jet_num >=1) iPhoEle1JetEvts++;
			if (stuple.jet_num >=2) iPhoEle2JetEvts++;
	  	} else if (stuple.ele_num > 0 && stuple.pho_num < 1) {
			iEleEvts++;
			if (stuple.jet_num >=1) iEle1JetEvts++;
			if (stuple.jet_num >=2) iEle2JetEvts++;
	  	} else if (stuple.ele_num < 1 && stuple.pho_num > 0) {
			iPhoEvts++;
			if (stuple.jet_num >=1) iPho1JetEvts++;
			if (stuple.jet_num >=2) iPho2JetEvts++;
		}	else iOther++;
			
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
	
	std::cout << "===========================================" << std::endl;
	std::cout << "Entries found    = " << iENTRIES << std::endl;
	std::cout << "Entries request  = " << iRunEvents << std::endl;
	std::cout << "Events Processed = " << iEvtProc << std::endl;
	std::cout << "Ele./Pho Events, +1Jet, +2Jet = " << iPhoEleEvts << "\t" << iPhoEle1JetEvts<< "\t" << iPhoEle2JetEvts << std::endl;
	std::cout << "Ele Events,      +1Jet, +2Jet = " << iEleEvts << "\t" << iEle1JetEvts<< "\t" << iEle2JetEvts << std::endl;
	std::cout << "Pho Events,      +1Jet, +2Jet = " << iPhoEvts << "\t" << iPho1JetEvts<< "\t" << iPho2JetEvts << std::endl;
	std::cout << "Other Events     = " << iOther << std::endl;
	std::cout << "First400 Events  = " << iFirst400 << std::endl;
	std::cout << "sum Events       = " << iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400 << std::endl;
	
	std::cout << "Events Selected  = " << iEventsSelected << std::endl;
	std::cout << "Data Written to  = " << GetHistFileName() << std::endl;
	HaloTypes ht;
	std::cout << "====== DelR >0.2 ==================" << std::endl;
	for (int i =0; i < ht.iHALOTYPES; ++i) std::cout << "Halos Found: type " << i << "\t= "  << vHalosFound_GT[i] << std::endl;
	std::cout << "====== DelR <0.2 ==================" << std::endl;
	for (int i =0; i < ht.iHALOTYPES; ++i) std::cout << "Halos Found: type " << i << "\t= "  << vHalosFound_LT[i] << std::endl;
	std::cout << "===========================================" << std::endl;
	CleanUp();
} // Main

