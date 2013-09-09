#include <iostream>
#include <sstream>
#include "PhoFlat/ewkMCPhotonJets.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
#include "PhoFlat/FreeFunctions.hh"
#include "PhoFlat/EleFakeRate.hh"

//-------------------------------------------------------------------
ewkMCPhotonJets::ewkMCPhotonJets():
iProgressBy(500000)
{
	iJES = 0;
	sHistFileName = "ewkMCPhotonJets.root";
}

//--- do all my stuff here ------------------------------------------
void ewkMCPhotonJets::DoMyStuff(Stuple& st)
{
	//if (st.pho_Ntight + st.pho_Nloose > 0 && st.ele_Ntight + st.ele_Nloose > 0)
	 // st.Dump ();
	
	Stuple stuple = st;

	//MC has no trig info
	//if ( (stuple.tri_pho25iso != 1)  && (stuple.tri_pho50 != 1) &&
	//	  (stuple.tri_pho70 != 1) ) return;
		  
	if (stuple.vtx_NClass12 <= 0) return;
	if (stuple.vtx_z > 60.)  return;
	if (stuple.pho_Ntight < 1) return;
	//if (stuple.ele_Ntight < 1) return;

	std::vector<int> UsedPho, UsedEle;

	if (iJES) CorrectPhotonEnergy(stuple, iJES);	
	
	std::vector<float> vWght, vWght_m, vWght_p;
	for (unsigned int i = 0; i < stuple.pho_num; ++i ) {
		//if (stuple.ele_TightId[i] != 0) continue;		// tight electrons
		//if (stuple.pho_TightId[i] != 0) continue;		// tight photons.
		if (!PassPhoTightIdCuts(stuple, i)) continue;	// this is needed with JES up and down 04-04-08
		
		//em time not simulated
		//if (stuple.pho_EmTime[i] < -4.8 || stuple.pho_EmTime[i] > 4.8) continue;  // select in-time photons
		//if (stuple.pho_PhoenixId[i] != 0) continue;

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
			continue;
		}

//		double eleEt = stuple.ele_Etc[i];
//		double wght_d = 0, wght_m =0, wght_p=0;
//		int iGoodRunBit = 0;
//		int iPhoenixBit = 1;
//		EleFakeRate fr;
//		fr.MyFakeRate(eleEt, iGoodRunBit, iPhoenixBit, wght_d, wght_m, wght_p);

//		vWght.push_back(wght_d);
//		vWght_m.push_back(wght_m);
//		vWght_p.push_back(wght_p);
	
		
		
		UsedPho.push_back(i);	
		NewJetList nj;
		nj.AddUnused(stuple,UsedPho,UsedEle);

		iEventsSelected++;


		// now look for jets>15 as i did not make a jet cut when skimming
		std::vector<int> vJet15Ind;
		//TLorentzVector tlPho(stuple.pho_Px[i], stuple.pho_Py[i], stuple.pho_Pz[i], stuple.pho_E[i]);
		TLorentzVector tlEle(stuple.ele_Px[i], stuple.ele_Py[i], stuple.ele_Pz[i], stuple.ele_E[i]);
		for (int j=0; j < stuple.jet_num; j++) {
			TLorentzVector tlJet(stuple.jet_Px[j], stuple.jet_Py[j], stuple.jet_Pz[j], stuple.jet_E[j]);
			float DelR = tlEle.DeltaR(tlJet);
			if (DelR<0.2) continue;			// _________________________this is a temp fix to em objects not being removed properly
													// _________________________from jet list. few events for photons. but a lot for eles  -- 04-04-08
			if (stuple.jet_Pt[j] < 15.0) continue;
			if (fabs(tlJet.Eta()) <3.0) vJet15Ind.push_back(j);		//Sasha cuts on eta only in MET clean up cuts 04-04-08
		}

		
		//now set the NJet15 correctly to fix the plots
		stuple.jet_NJet15 = (UInt_t) vJet15Ind.size();


		
	//if (stuple.jet_num != vJet15Ind.size()) stuple.Dump();  need to check this. I am saving all jets.
	// But I don't find any jets <15GeV in the skim??? -- 04-04-08

		HistManager myHistMan;
		if (vJet15Ind.size() >= 1) {	// mono-jet case
		//if (stuple.jet_num >= 1) {	// mono-jet case
			iPho1JetEvts++;
			myHistMan.FillPhotonHists(stuple, GA_Pho,i);
			myHistMan.FillEventHists(stuple, GA_Evt);
			myHistMan.FillJetHists(stuple, GA_Jet,vJet15Ind[0]);		//lead jet
			myHistMan.FillPhoton1JetHists(stuple, GA_PhoJet,i,vJet15Ind[0]);
//			myHistMan.FillElectronHists(stuple, GA_Pho,i, wght_d);
//			myHistMan.FillEventHists(stuple, GA_Evt, wght_d);
//			myHistMan.FillJetHists(stuple, GA_Jet,vJet15Ind[0], wght_d);		//lead jet
//			myHistMan.FillElectron1JetHists(stuple, GA_PhoJet,i,vJet15Ind[0], wght_d);
		}
		
		if (vJet15Ind.size() >= 2) {		// di-jet case
		//if (stuple.jet_num >= 2) {	// mono-jet case
			iPho2JetEvts++;
			myHistMan.FillPhotonHists(stuple, GB_Pho,i);
			myHistMan.FillEventHists(stuple, GB_Evt);
			myHistMan.FillJetHists(stuple, GB_Jet1,vJet15Ind[0]);		//lead jet
			myHistMan.FillJetHists(stuple, GB_Jet2,vJet15Ind[1]);		//2nd lead jet
			myHistMan.FillPhoton1JetHists(stuple, GB_PhoJet1,i,vJet15Ind[0]);
			myHistMan.FillPhoton1JetHists(stuple, GB_PhoJet2,i,vJet15Ind[0]);
			myHistMan.FillPhoton2JetsHists(stuple, GB_PhoJets,i,vJet15Ind[0],vJet15Ind[1]);
			myHistMan.FillTwoJetsHists(stuple, GB_TwoJets,vJet15Ind[0],vJet15Ind[1]);
			
			//myHistMan.FillElectronHists(stuple, GB_Pho,i, wght_d);
			//myHistMan.FillEventHists(stuple, GB_Evt, wght_d);
			//myHistMan.FillJetHists(stuple, GB_Jet1,vJet15Ind[0], wght_d);		//lead jet
			//myHistMan.FillJetHists(stuple, GB_Jet2,vJet15Ind[1], wght_d);		//2nd lead jet
			//myHistMan.FillElectron1JetHists(stuple, GB_PhoJet1,i,vJet15Ind[0], wght_d);
			//myHistMan.FillElectron1JetHists(stuple, GB_PhoJet2,i,vJet15Ind[0], wght_d);
			//myHistMan.FillElectron2JetsHists(stuple, GB_PhoJets,i,vJet15Ind[0],vJet15Ind[1], wght_d);
			//myHistMan.FillTwoJetsHists(stuple, GB_TwoJets,vJet15Ind[0],vJet15Ind[1], wght_d);
		}


		break;
	} //for 


}


//-------------------------------------------------------------
bool ewkMCPhotonJets::PassPhoTightIdCuts(Stuple& st, int i)
{
	if (st.pho_Detector[i] !=0 ) return false;
//	std::cout << "1 pass " << "det " << std::endl;
	if (st.pho_Etc[i] < 30.) return false; 
//	std::cout << "2 pass " << "et " << std::endl;
	if (fabs(st.pho_XCes[i]) > 21)  return false;
//	std::cout << "3 pass " << "xcez " << std::endl;
	if (fabs(st.pho_ZCes[i]) < 9 ||  fabs(st.pho_ZCes[i]) > 230)  return false;
//	std::cout << "4 pass " << "zces " << std::endl;
	if ( !((st.pho_HadEm[i] < 0.125) || (st.pho_HadEm[i] < (0.055 + 0.00045 * st.pho_E[i]))))  return false;
//	std::cout << "5 pass " << "hadem " << std::endl;
	if (st.pho_Etc[i] < 20) {
		if (st.pho_Iso4[i] > 0.1*st.pho_Etc[i]) return false;
	} else {
		if (st.pho_Iso4[i] > (2.0+0.02*(st.pho_Etc[i]-20.0)) )  return false;
	}
//	std::cout << "6 pass " << "iso " << std::endl;
	
	if (st.pho_Chi2Mean[i] > 20) return false;
//	std::cout << "7 pass " << "chi2mean " << std::endl;
	if (st.pho_N3d[i] > 1)  return false;
	//std::cout << "8 pass " << "n3d " << std::endl;
	if (st.pho_TrkPt[i] > (1+0.005*st.pho_Etc[i])) return false;
//	std::cout << "9 pass " << "trkpt " << std::endl;
	if (st.pho_TrkIso[i] > (2.0+0.005*st.pho_Etc[i])) return false;
//	std::cout << "10 pass " << "trkiso " << std::endl;

	TLorentzVector tlEle(st.pho_Px[i], st.pho_Py[i], st.pho_Pz[i], st.pho_E[i]);

	if (st.pho_Etc[i]<18) {
		if ( (st.pho_CesWireE2[i] * TMath::Sin(tlEle.Theta())) > (0.14*st.pho_Etc[i]))  return false;
		if ( (st.pho_CesStripE2[i] * TMath::Sin(tlEle.Theta())) > (0.14*st.pho_Etc[i]))  return false;
	} else {
		if ( (st.pho_CesWireE2[i] * TMath::Sin(tlEle.Theta())) > (2.4+0.01* st.pho_Etc[i]))  return false;
		if ( (st.pho_CesStripE2[i] * TMath::Sin(tlEle.Theta())) > (2.4+0.01* st.pho_Etc[i]))  return false;
	}
//	std::cout << "12 pass " << " all" << std::endl;
	
	return true;


}


//-------------------------------------------------------------
void ewkMCPhotonJets::Init(TChain *tree)
{
	iEventsSelected = 0;
	iPho1JetEvts = 0;
	iPho2JetEvts = 0;
	readIn = new ReadInAll(tree,&stuple);

	//cannot add folder without creating a dir first!
	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = rootFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

} // Init

//-------------------------------------------------------------------
void ewkMCPhotonJets::CleanUp()
{
	rootFile->Write();
	rootFile->Close();
} 

//-------------------------------------------------------------------
void ewkMCPhotonJets::BookHistograms()
{

	//for halo estimate
	haloPhiWedge_1j = new TH1F("haloPhiWedge_1j","Halo(5)+>=1 Jet: #phi^{#gamma}",24,0,24);
	haloEmTime_1j = new TH1F("haloEmTime_1j","Halo(5)+>=1 Jet: #gamma EmTime",30,-10,10);
	haloPhiWedge_2j = new TH1F("haloPhiWedge_2j","Halo(5)+>=2 Jet: #phi^{#gamma}",24,0,24);
	haloEmTime_2j = new TH1F("haloEmTime_2j","Halo(5)+>=2 Jet: #gamma EmTime",30,-10,10);
	
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
void ewkMCPhotonJets::Main(TChain *ch, int iRunEvents)
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
		
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;


		DoMyStuff(stuple);
	};
	
	std::cout << "======== ewkMCPhotonJets =======================" << std::endl;
	std::cout << "Entries found    = " << iENTRIES << std::endl;
	std::cout << "Entries request  = " << iRunEvents << std::endl;
	std::cout << "Events Processed = " << iEvtProc << std::endl;
	std::cout << "Ele./Pho Events  = " << iPhoEleEvts << std::endl;
	std::cout << "Ele Events       = " << iEleEvts << std::endl;
	std::cout << "Pho Events       = " << iPhoEvts << std::endl;
	std::cout << "Other Events     = " << iOther << std::endl;
	std::cout << "sum Events       = " << iPhoEvts + iEleEvts + iPhoEleEvts + iOther << std::endl;
	std::cout << "__________________________________________" << std::endl;
	std::cout << "\n!!!CHECK THAT YOU ARE RUNNING ON MC!!!\n" << std::endl;
	std::cout << "Pho Evts               = " << iEventsSelected << std::endl;
	std::cout << "Pho+>=1Jet Evts        = " << iPho1JetEvts << std::endl; 
	std::cout << "Pho+>=2Jet Evts        = " << iPho2JetEvts << std::endl; 
	std::cout << "JES(0=def,1=up,-1=down)= " << GetJES() << std::endl; 
	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	HaloTypes ht;
	std::cout << "======= END ewkMCPhotonJets ====================" << std::endl;
	CleanUp();
} // Main

