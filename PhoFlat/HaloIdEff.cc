#include <iostream>
#include <sstream>
#include "PhoFlat/HaloIdEff.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "TLorentzVector.h"
#include <iomanip>
#include "../RootTools/CommonTools.hh"
//-------------------------------------------------------------------
HaloIdEff::HaloIdEff():
iProgressBy(10000),
sHistFileName("HaloIdEff_novtx_nojet.root")
{
	fMinPhoEt = 30.0;
	fMaxPhoEta = 1.1;
	fMinMet = 0.0;
	iHaloType = 5;
	iNMinVtx = 0;
	iNMaxVtx = 0;
	iMinNjet15 = 0;
}

//--- do all my stuff here ------------------------------------------
void HaloIdEff::DoMyStuff(Stuple& stuple)
{
	if (stuple.vtx_NClass12 != 0)  return;		// select no vtx events
	if (stuple.pho_num < 1) return;
	if (stuple.ele_num > 0) return;  // no need. remove this and see. this should not matter. all we care about is a pure halo event. effect of this should be small as  i have few e+p evets
	if (stuple.met_Met < fMinMet) return; //met cut
	if (stuple.jet_NJet15 < iMinNjet15) return;

	const int iPhoCollection = 0; //central photon collection
	
	for (unsigned int i = 0; i < stuple.pho_num; ++i ) {
		if (stuple.pho_TightId[i] != 0) continue;		// tight photons
		if (fabs(stuple.pho_DetEta[i]) > fMaxPhoEta) continue;		// eta cut
		if (fabs(stuple.pho_Etc[i]) < fMinPhoEt) continue;		// et cut

		int j = -1;
		if (stuple.pho_EmTime[i] > -4.8 && stuple.pho_EmTime[i] < 4.8) {
			// this halo+cosmics with no vtx events b4/a4 cuts
			j = 0; //b4 halo cuts
					vHalosFound[j]++;
					myHistMan.FillPhotonHists(stuple, iPhoCollection, GA_Pho[j],i);
					myHistMan.FillEventHists(stuple, GA_Evt[j]);
				HaloTypes ht(stuple.pho_Halo_seedWedge[i], 
					stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);

				j = 2; //halo+cosmic a4 cuts
				if (ht.IsType(iHaloType)) {
						vHalosFound[j]++;
						myHistMan.FillPhotonHists(stuple, iPhoCollection, GA_Pho[j],i);
						myHistMan.FillEventHists(stuple, GA_Evt[j]);
				} //a4 cuts
				
			} // in-time halo+cosmic		



		if (stuple.pho_EmTime[i] > 30. && stuple.pho_EmTime[i] < 90.) {
			// cosmics with no vtx events b4/a4 cuts
			j = 1; //b4 halo cuts
					vHalosFound[j]++;
					myHistMan.FillPhotonHists(stuple, iPhoCollection, GA_Pho[j],i);
					myHistMan.FillEventHists(stuple, GA_Evt[j]);

				HaloTypes ht(stuple.pho_Halo_seedWedge[i], 
					stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);

				j = 3; //cosmic a4 cuts
				if (ht.IsType(iHaloType)) {
			
						vHalosFound[j]++;
						myHistMan.FillPhotonHists(stuple, iPhoCollection, GA_Pho[j],i);
						myHistMan.FillEventHists(stuple, GA_Evt[j]);
				} //a4 cuts
				
			} // cosmic		

		break;
	} //for 


}

//-------------------------------------------------------------
void HaloIdEff::Init(TChain *tree)
{
	HaloTypes ht;
	for (int i =0; i <HaloTypes::iNHALOTYPES; i++) vHalosFound.push_back(0);
	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	//cannot add folder without creating a dir first!
	topDir = gDirectory->mkdir("Hist","All Histograms");	//create top level dir.

	iEventsSelected = 0;
	readIn = new ReadInAll(tree, &stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->Met(1);
	readIn->Init();

} // Init

//-------------------------------------------------------------------
void HaloIdEff::CleanUp()
{
	//vHalosFound.clear();
	for (int i =0; i < HaloTypes::iNHALOTYPES; i++) vHalosFound.push_back(0);
	rootFile->Write();
	rootFile->Close();
} 

//-------------------------------------------------------------------
void HaloIdEff::BookHistograms()
{
	Histograms::EventHists_t hEvt;
	Histograms::PhotonHists_t hPho;
	Histograms::JetHists_t hJet;
	Histograms::Photon1JetHists_t hPhoJet;
	Histograms::Photon2JetsHists_t hPhoJets;
	Histograms::TwoJetsHists_t hTwoJets;

	Histograms HistoMan;
	//hist array
	// [0] = Halo+cosmic b4 halo cuts
	// [1] = cosmic only b4 halo cuts
	// [2] = Halo+cosmic a4 halo cuts
	// [3] = cosmic only a4 halo cuts

	TDirectory *top, *sub,*subsub;
	for (int i =0; i <4; i++) {
		
		GA_Evt.push_back(hEvt);
//		GB_Evt.push_back(hEvt);
		GA_Pho.push_back(hPho);
//		GB_Pho.push_back(hPho);
/*		GA_Jet.push_back(hJet);
		GB_Jet1.push_back(hJet);
		GB_Jet2.push_back(hJet);
		GA_PhoJet.push_back(hPhoJet);
		GB_PhoJet1.push_back(hPhoJet);
		GB_PhoJet2.push_back(hPhoJet);
		GB_PhoJets.push_back(hPhoJets);
		GB_TwoJets.push_back(hTwoJets);
*/
		std::ostringstream foldername, text2title1, text2title2;
		if (i == 0) {
			foldername << "HaloAndCosmic_b4";
			//text2title1 << "Nvtx=0 : EMtime < |4.8|ns : MET>25 : Halo+Cosmics (+>= 1 Jet) before Halo(5) cuts";
			//text2title2 << "Nvtx=0 : EMtime < |4.8|ns : MET>25 : Halo+Cosmics (+>= 2 Jets) before Halo(5) cuts";
			text2title1 << "EMtime < |4.8|ns : Jets>=0 : Halo+Cosmics before Halo(5) cuts";
		} else if (i == 1) {
			foldername << "Cosmic_b4";
			//text2title1 << "Nvtx=0 : 30ns<EMtime<90ns : MET>25 : Cosmics (+>= 1 Jet) before Halo(5) cuts";
			//text2title2 << "Nvtx=0 : 30ns<EMtime<90ns : MET>25 : Cosmics (+>= 2 Jets) before Halo(5) cuts";
			text2title1 << "Nvtx=0 : 30ns<EMtime<90ns :Jets>=0 : Cosmics before Halo(5) cuts";
		} else if (i == 2) {
			foldername << "HaloAndCosmic_a4";
			//text2title1 << "Nvtx=0 : EMtime < |4.8|ns : MET>25 : Halo+Cosmics (+>= 1 Jet) after Halo(5) cuts";
			//text2title2 << "Nvtx=0 : EMtime < |4.8|ns : MET>25 : Halo+Cosmics (+>= 2 Jets) after Halo(5) cuts";
			text2title1 << "Nvtx=0 : EMtime < |4.8|ns :Jets>=0 : Halo+Cosmics after Halo(5) cuts";
		} else if (i == 3) {
			foldername << "Cosmic_a4";
			//text2title1 << "Nvtx=0 : 30ns<EMtime<90ns : MET>25 : Cosmics (+>= 1 Jet) after Halo(5) cuts";
			//text2title2 << "Nvtx=0 : 30ns<EMtime<90ns : MET>25 : Cosmics (+>= 2 Jets) after Halo(5) cuts";
			text2title1 << "Nvtx=0 : 30ns<EMtime<90ns :Jets>=0 : Cosmics after Halo(5) cuts";
		}
			
		topDir->cd();
		top = topDir->mkdir(foldername.str().c_str(), foldername.str().c_str());
		top->cd();
	
//		sub = gDirectory->mkdir("1Jet","1 Jet case");
//		sub->cd();
		
		subsub = gDirectory->mkdir("Event","Event Histograms");
		subsub->cd();
				HistoMan.GetEventHistograms(GA_Evt[i],text2title1.str().c_str());
				
		//sub->cd();
		top->cd();
		subsub = gDirectory->mkdir("Photon","Photon Histograms");
		subsub->cd();
				HistoMan.GetPhotonHistograms(GA_Pho[i],text2title1.str().c_str());

/*		sub->cd();
		subsub = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		subsub->cd();
				HistoMan.GetJetHistograms(GA_Jet[i],text2title1.str().c_str());
				
		sub->cd();		
		subsub = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		subsub->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(GA_PhoJet[i],text2title1.str().c_str());
*/

	/*
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
				*/

	}


} //BookHistograms

//------ main -------------------------------------------------------
void HaloIdEff::Main(TChain *ch, int iRunEvents)
{
	if (ch) {
		myChain = ch;
		//ch->Print();	
	} else {
		std::cout << "NULL chain!" << std::endl;
		exit(1);
	}
	
	PrintJobSettings();
	
	gROOT->Reset();
	Init(myChain);
	BookHistograms();
	
	unsigned int entries = (unsigned int) myChain->GetEntries();
	unsigned int ENTRIES = entries;

	if (iRunEvents > 0 && iRunEvents < entries) entries = iRunEvents;
	else std::cout << __FILE__ << "::" << __LINE__ << "::" << "Only " << entries << " found and those will be processed." << std::endl;

	if (iRunEvents < 0 ) iRunEvents = entries;	//requested all entries


	  readIn->GetEntry(1);
	  if (stuple.evt_McFlag != 0) 
	  {
		  std::cout << red << "THIS NOT DATA. Read from stuple " <<
			  stuple.evt_McFlag << ". pleas check. returning!" << clearatt << std::endl;
		  assert (false);
	  }

	
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
			
		//std::cout << "phonum,nvtx " << stuple.pho_num << ", " << stuple.vtx_NClass12 << std::endl;
			
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;
			
		DoMyStuff(stuple);
	};
	
	std::cout << "===========================================" << entries << std::endl;
	std::cout << "Entries found    = " << ENTRIES<< std::endl;
	std::cout << "Entries request  = " << iRunEvents << std::endl;
	std::cout << "Events Processed = " << iEvtProc << std::endl;
	std::cout << "Ele./Pho Events  = " << iPhoEleEvts << std::endl;
	std::cout << "Ele Events       = " << iEleEvts << std::endl;
	std::cout << "Pho Events       = " << iPhoEvts << std::endl;
	std::cout << "Other Events     = " << iOther << std::endl;
	std::cout << "First400 Events  = " << iFirst400 << std::endl;
	std::cout << "sum Events       = " << iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400 << std::endl;
	
	PrintJobSettings();
	
	std::cout << "Events Selected  = " << iEventsSelected << std::endl;
	std::cout << "Data Written to  = " << GetHistFileName() << std::endl;
	HaloTypes ht;
	for (int i =0; i < 4; ++i) {
		if (i == 0) std::cout << "Halo+cosmic (+1Jet) b4 cuts =\t"  << vHalosFound[i] << std::endl;
		if (i == 1) std::cout << "cosmic (+1Jet)      b4 cuts =\t"  << vHalosFound[i] << std::endl;
		if (i == 2) std::cout << "Halo+cosmic (+1Jet) a4 cuts =\t"  << vHalosFound[i] << std::endl;
		if (i == 3) std::cout << "cosmic (+1Jet)      a4 cuts =\t"  << vHalosFound[i] << std::endl;
	}
	std::cout << "===========================================" << entries << std::endl;
	CleanUp();
} // Main

void HaloIdEff::PrintJobSettings()
{
	std::cout << "Min Pho Et, Eta  = " << std::setw(4) << fMinPhoEt << std::setw(4) << fMaxPhoEta << std::endl;
	std::cout << "Min/Max Vtx      = " << std::setw(4) << iNMinVtx << ", " << std::setw(3) << iNMaxVtx << std::endl;
	std::cout << "Halo Type        = " << std::setw(4) << iHaloType << std::endl;
	std::cout << "Min MEt          = " << std::setw(4) << fMinMet << std::endl;
	std::cout << "Min NJet15       = " << std::setw(4) << iMinNjet15 << std::endl;
}
