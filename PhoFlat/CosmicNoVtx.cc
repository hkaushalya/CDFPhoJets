////////////////////////////////////////////////////////////////
// To look at cosmics in no vtx sample. I think I was trying  //
// get a systematic error for cosmic estimate ann Max         //
// suggested that using cosmic in no vtx sample and fitting it//
// to cosmic in signal and then making another estimate from  //
// it. The difference can be taken as the systematic error.   //
////////////////////////////////////////////////////////////////


#include <iostream>
#include <sstream>
#include "PhoFlat/CosmicNoVtx.hh"
#include "PhoFlat/ReadInAll.hh"
#include "TLorentzVector.h"
#include "TROOT.h"
//-------------------------------------------------------------------
CosmicNoVtx::CosmicNoVtx():
iProgressBy(100000),
sHistFileName("CosmicNoVtx.root")
{
}

//--- do all my stuff here ------------------------------------------
void CosmicNoVtx::DoMyStuff(Stuple& stuple)
{
	if (stuple.vtx_NClass12 != 0)  return;		// select no vtx events
	if (stuple.pho_num < 1) return;
	if (stuple.ele_num > 0) return;  // no need. remove this and see. this should not matter. all we care about is a pure halo event. effect of this should be small as  i have few e+p evets
	for (unsigned int i = 0; i < stuple.pho_num; ++i ) {
		if (stuple.pho_TightId[i] != 0) continue;		// tight photons
		if (stuple.pho_Etc[i] < 30.) continue;
		if (fabs(stuple.pho_DetEta[i]) >1.1) continue;
		if (stuple.pho_NMuonStubs[i] < 1) continue;
		
		//if (stuple.pho_EmTime[i] > 30.0 && stuple.pho_EmTime[i] < 90.0) {
			myHistMan.FillPhotonHists(stuple,0, GA_Pho[0],i);
			myHistMan.FillEventHists(stuple, GA_Evt[0]);
			iEventsSelected++;
		//} // cosmic		
			
		//break;
	} //for 


}

//-------------------------------------------------------------
void CosmicNoVtx::Init(TChain *tree)
{
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
void CosmicNoVtx::CleanUp()
{
	//vHalosFound.clear();
	for (int i =0; i <4; i++) vHalosFound.push_back(0);
	rootFile->Write();
	rootFile->Close();
} 

//-------------------------------------------------------------------
void CosmicNoVtx::BookHistograms()
{
	Histograms::EventHists_t hEvt;
	Histograms::PhotonHists_t hPho;

	Histograms HistoMan;

	TDirectory *top, *sub;
	for (int i =0; i <1; i++) {
		
		GA_Evt.push_back(hEvt);
		GA_Pho.push_back(hPho);
		std::ostringstream foldername, text2title1, text2title2;
		if (i == 0) {
			foldername << "Cosmics";
			text2title1 << "EMtime >30ns & <90 ns : Jets>=0 : Cosmics from no vtx evts";
		}
			
		topDir->cd();
		top = topDir->mkdir(foldername.str().c_str(), foldername.str().c_str());
		top->cd();
	
		sub = gDirectory->mkdir("Event","Event Histograms");
		sub->cd();
				HistoMan.GetEventHistograms(GA_Evt[i],text2title1.str().c_str());
				
		top->cd();
		sub = gDirectory->mkdir("Photon","Photon Histograms");
		sub->cd();
				HistoMan.GetPhotonHistograms(GA_Pho[i],text2title1.str().c_str());

	}


} //BookHistograms

//------ main -------------------------------------------------------
void CosmicNoVtx::Main(TChain *ch, int iRunEvents)
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
	unsigned int ENTRIES = entries;

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
			
//		std::cout << "njet15=" << stuple.jet_NJet15 << std::endl;
//		std::cout << "jetnum=" << stuple.jet_num << std::endl;
			
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;
			
		DoMyStuff(stuple);
	};
	
	std::cout << "======= CosmicNoVtx =======================" << entries << std::endl;
	std::cout << "Entries found    = " << ENTRIES<< std::endl;
	std::cout << "Entries request  = " << iRunEvents << std::endl;
	std::cout << "Events Processed = " << iEvtProc << std::endl;
	std::cout << "Ele./Pho Events  = " << iPhoEleEvts << std::endl;
	std::cout << "Ele Events       = " << iEleEvts << std::endl;
	std::cout << "Pho Events       = " << iPhoEvts << std::endl;
	std::cout << "Other Events     = " << iOther << std::endl;
	std::cout << "First400 Events  = " << iFirst400 << std::endl;
	std::cout << "sum Events       = " << iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400 << std::endl;
	std::cout << "___________________________________________" << iEventsSelected << std::endl;
	std::cout << "Events Selected  = " << iEventsSelected << std::endl;
	std::cout << "Data Written to  = " << GetHistFileName() << std::endl;
	std::cout << "======= end CosmicNoVtx ===================" << entries << std::endl;
	CleanUp();
} // Main

