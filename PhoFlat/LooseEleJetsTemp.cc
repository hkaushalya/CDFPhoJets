#include <iostream>
#include <sstream>
#include "PhoFlat/LooseEleJetsTemp.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
#include "PhoFlat/EleFakeRate.hh"
#include "TMath.h"

//-------------------------------------------------------------------
LooseEleJetsTemp::LooseEleJetsTemp():
iProgressBy(100000),
sHistFileName("LooseEleJetsTemp.root")
{
}

//--- do all my stuff here ------------------------------------------
void LooseEleJetsTemp::DoMyStuff(Stuple& stuple)
{

	if ( (stuple.tri_pho25iso != 1)  && (stuple.tri_pho50 != 1) &&
		  (stuple.tri_pho70 != 1) ) return;
		  
	if (stuple.vtx_NClass12 <= 0) return;
	if (stuple.vtx_z > 60.)  return;
	if (stuple.ele_Nloose < 1) return;
	
	std::vector<int> UsedPho, UsedEle;
	
	std::vector<float> vWght, vWght_m, vWght_p;

	//i want to get the sum of weights if there is more than one electron in the event.
	//but still add the electrons(loose) that are not used(if any).

	for (unsigned int i = 0; i < stuple.ele_num; ++i ) {
		if (stuple.ele_LooseId[i] != 0) continue;		// tight photons
		if (stuple.ele_EmTime[i] < -4.8 || stuple.pho_EmTime[i] > 4.8) continue;  // select in-time photons

		UsedEle.push_back(i);	

		double eleEt = stuple.ele_Etc[i];
		double wght_d = 0, wght_m =0, wght_p=0;
		int iGoodRunBit = 0;
		int iPhoenixBit = 1;
		EleFakeRate fr;
		fr.MyFakeRate(eleEt, iGoodRunBit, iPhoenixBit, wght_d, wght_m, wght_p);

		vWght.push_back(wght_d);
		vWght_m.push_back(wght_m);
		vWght_p.push_back(wght_p);
		iEventsSelected++;
	}
	
	NewJetList nj;
	nj.AddUnused(stuple, UsedPho, UsedEle);
	HistManager myHistMan;

	for (int i=0; i < UsedEle.size(); ++i) {
		if (stuple.jet_NJet15 >= 1) {	// mono-jet case
			iEle1JetEvts++;
			fSFR_1j += vWght[i];
			fSFR_m_1j += vWght_m[i];
			fSFR_p_1j += vWght_p[i];
			myHistMan.FillPhotonHists(stuple, GA_Pho,UsedEle[i], vWght[i]);
			myHistMan.FillEventHists(stuple, GA_Evt, vWght[i]);
			myHistMan.FillJetHists(stuple, GA_Jet,0, vWght[i]);		//lead jet
			myHistMan.FillPhoton1JetHists(stuple, GA_PhoJet,UsedEle[i],0, vWght[i]);
		}
		
		if (stuple.jet_NJet15 >= 2) {	// di-jet case
			iEle2JetEvts++;
			fSFR_2j += vWght[i];
			fSFR_m_2j += vWght_m[i];
			fSFR_p_2j += vWght_p[i];
			myHistMan.FillPhotonHists(stuple, GB_Pho,UsedEle[i], vWght[i]);
			myHistMan.FillEventHists(stuple, GB_Evt, vWght[i]);
			myHistMan.FillJetHists(stuple, GB_Jet1,0, vWght[i]);		//lead jet
			myHistMan.FillJetHists(stuple, GB_Jet2,1, vWght[i]);		//2nd lead jet
			myHistMan.FillPhoton1JetHists(stuple, GB_PhoJet1,UsedEle[i],0, vWght[i]);
			myHistMan.FillPhoton1JetHists(stuple, GB_PhoJet2,UsedEle[i],1, vWght[i]);
			myHistMan.FillPhoton2JetsHists(stuple, GB_PhoJets,UsedEle[i],0,1, vWght[i]);
			myHistMan.FillTwoJetsHists(stuple, GB_TwoJets,0,1, vWght[i]);
		}
	}

	}

//-------------------------------------------------------------
void LooseEleJetsTemp::Init(TChain *tree)
{
	iEventsSelected = 0;
	iEle1JetEvts = 0;
	iEle2JetEvts = 0;

	fSFR_1j =0; fSFR_m_1j=0; fSFR_p_1j=0;			// sum of fake rates
	fSFR_2j=0; fSFR_m_2j=0; fSFR_p_2j=0;			// sum of fake rates

	
	readIn = new ReadInAll(tree,&stuple);

	//cannot add folder without creating a dir first!
	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = rootFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

} // Init

//-------------------------------------------------------------------
void LooseEleJetsTemp::CleanUp()
{
	rootFile->Write();
	rootFile->Close();
} 

//-------------------------------------------------------------------
void LooseEleJetsTemp::BookHistograms()
{
	Histograms HistoMan;
	std::ostringstream foldername, text2title1, text2title2;
	text2title1 << "e+>=1Jets :";
	text2title2 << "e+>=2Jets :";
	

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
void LooseEleJetsTemp::Main(TChain *ch, int iRunEvents)
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
	
	std::cout << "======== LooseEleJetsTemp ======================" << std::endl;
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
	std::cout << "Ele Evts         = " << iEventsSelected << std::endl;
	std::cout << "Ele+>=1Jet Evts  = " << iEle1JetEvts << std::endl; 
	std::cout << "Ele+>=2Jet Evts  = " << iEle2JetEvts << std::endl; 
	std::cout << "Sum FR_d (1Jet)  = " << fSFR_1j << " +/- " << TMath::Max(fabs(fSFR_m_1j-fSFR_1j), fabs(fSFR_p_1j-fSFR_1j)) << std::endl;
	std::cout << "Sum FR_d (2Jet)  = " << fSFR_2j<< " +/- " << TMath::Max(fabs(fSFR_m_2j-fSFR_2j), fabs(fSFR_p_2j -fSFR_2j)) << std::endl;
	std::cout << "Data Written to  = " << GetHistFileName() << std::endl;
	std::cout << "======= END LooseEleJetsTemp ===================" << std::endl;
	CleanUp();
} // Main

