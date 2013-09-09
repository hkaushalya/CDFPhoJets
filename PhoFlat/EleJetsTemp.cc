#include <iostream>
#include <sstream>
#include "PhoFlat/EleJetsTemp.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
#include "PhoFlat/EleFakeRate.hh"
#include "TMath.h"
#include "TROOT.h"
#include <iomanip>
#include "TLatex.h"
#include "TCanvas.h"
#include <sstream>
#include "TDatime.h"

//-------------------------------------------------------------------
EleJetsTemp::EleJetsTemp():
iProgressBy(100000)
{
	sHistFileName = "EleJetsTemp.root";
	fMaxEleEta = 1.1;
	fMinEleEt = 30.;
	iMinClass12Vtx = 1;
	iMaxClass12Vtx = 100;
	fMinMet = 0.;
	fMinSignalEmTime = -4.8;
	fMaxSignalEmTime = 4.8;
	fMaxVtxz = 60.0;
	
}

//--- do all my stuff here ------------------------------------------
void EleJetsTemp::DoMyStuff(Stuple& stuple)
{
	if (! bMcFlag)
	{
		if ( (stuple.tri_pho25iso != 1)  && (stuple.tri_pho50 != 1) &&
				(stuple.tri_pho70 != 1) ) return;
	}

	if (stuple.vtx_NClass12 < iMinClass12Vtx 
			|| stuple.vtx_NClass12 > iMaxClass12Vtx
			|| fabs(stuple.vtx_z) > fMaxVtxz) return;
	if (stuple.met_Met < GetMinMet()) return;
	
	std::vector<int> UsedPho, UsedEle;
	
	std::vector<float> vWght, vWght_m, vWght_p;

	//i want to get the sum of weights if there is more than one electron in the event.
	//but still add the electrons(loose) that are not used(if any).

	for (unsigned int i = 0; i < stuple.ele_num; ++i ) {
		if (stuple.ele_TightId[i] != 0) continue;		// tight electrons
		if (fabs(stuple.ele_DetEta[i]) > fMaxEleEta) continue;
		if (fabs(stuple.ele_Etc[i]) < fMinEleEt) continue;
		if (fabs(stuple.ele_ConversionId[i]) != 0) continue;
		
		if (! bMcFlag) 
		{
			if (stuple.ele_EmTime[i] < -4.8 || stuple.pho_EmTime[i] > 4.8) continue;  // select in-time electrons
		}

		UsedEle.push_back(i);	

//I am commenting it out for now to study EWK component in data and MC - 03-16-2010
/*		double eleEt = stuple.ele_Etc[i];
		double wght_d = 0, wght_m =0, wght_p=0;
		int iGoodRunBit = 0;
		int iPhoenixBit = 1;
		EleFakeRate fr;
		fr.MyFakeRate(eleEt, iGoodRunBit, iPhoenixBit, wght_d, wght_m, wght_p);

		vWght.push_back(wght_d);
		vWght_m.push_back(wght_m);
		vWght_p.push_back(wght_p);
*/
		iEventsSelected++;
		break;  //look for only one electron
	}
	

	NewJetList nj;
	nj.AddUnused(stuple, UsedPho, UsedEle);

	HistManager myHistMan;

	//to study EWK component in DATA and MC
	if (UsedEle.size()>0) {
		const int obj_collection = 0; //pick central values of all objects
		if (stuple.jet_NJet15 >= 1) {	// mono-jet case
			iEle1JetEvts++;
			myHistMan.FillElectronHists(stuple,obj_collection, GA_Pho,UsedEle[0], 1);
			myHistMan.FillEventHists(stuple, GA_Evt, 1);
			myHistMan.FillJetHists(stuple,obj_collection, GA_Jet,0, 1);		//lead jet
			myHistMan.FillElectron1JetHists(stuple,obj_collection, GA_PhoJet,UsedEle[0],0, 1);

		}
		
		if (stuple.jet_NJet15 >= 2) {	// di-jet case
			iEle2JetEvts++;
			assert (stuple.jet_Pt[0]>=stuple.jet_Pt[1] && "EleJetsTemp::DoMyStuff()::Jets are not sorted!");
			myHistMan.FillElectronHists(stuple,obj_collection, GB_Pho,UsedEle[0], 1);
			myHistMan.FillEventHists(stuple, GB_Evt, 1);
			myHistMan.FillJetHists(stuple,obj_collection, GB_Jet1,0, 1);		//lead jet
			myHistMan.FillJetHists(stuple,obj_collection, GB_Jet2,1, 1);		//2nd lead jet
			myHistMan.FillElectron1JetHists(stuple,obj_collection, GB_PhoJet1,UsedEle[0],0, 1);
			myHistMan.FillElectron1JetHists(stuple,obj_collection, GB_PhoJet2,UsedEle[0],1, 1);
			myHistMan.FillElectron2JetsHists(stuple,obj_collection, GB_PhoJets,UsedEle[0],0,1, 1);
			myHistMan.FillTwoJetsHists(stuple,obj_collection, GB_TwoJets,0,1, 1);
		}
	}



	

// I am not sure what this is. This was already commented out. Now I need to EWK component in data and MC - 03-16-2010
	
/*
	for (int i=0; i < UsedEle.size(); ++i) {
		if (stuple.jet_NJet15 >= 1) {	// mono-jet case
			iEle1JetEvts++;
			fSFR_1j += vWght[i];
			fSFR_m_1j += vWght_m[i];
			fSFR_p_1j += vWght_p[i];
			myHistMan.FillElectronHists(stuple, GA_Pho,UsedEle[i], vWght[i]);
			myHistMan.FillEventHists(stuple, GA_Evt, vWght[i]);
			myHistMan.FillJetHists(stuple, GA_Jet,0, vWght[i]);		//lead jet
			myHistMan.FillElectron1JetHists(stuple, GA_PhoJet,UsedEle[i],0, vWght[i]);
		}
		
		if (stuple.jet_NJet15 >= 2) {	// di-jet case
			iEle2JetEvts++;
			fSFR_2j += vWght[i];
			fSFR_m_2j += vWght_m[i];
			fSFR_p_2j += vWght_p[i];
			myHistMan.FillElectronHists(stuple, GB_Pho,UsedEle[i], vWght[i]);
			myHistMan.FillEventHists(stuple, GB_Evt, vWght[i]);
			myHistMan.FillJetHists(stuple, GB_Jet1,0, vWght[i]);		//lead jet
			myHistMan.FillJetHists(stuple, GB_Jet2,1, vWght[i]);		//2nd lead jet
			myHistMan.FillElectron1JetHists(stuple, GB_PhoJet1,UsedEle[i],0, vWght[i]);
			myHistMan.FillElectron1JetHists(stuple, GB_PhoJet2,UsedEle[i],1, vWght[i]);
			myHistMan.FillElectron2JetsHists(stuple, GB_PhoJets,UsedEle[i],0,1, vWght[i]);
			myHistMan.FillTwoJetsHists(stuple, GB_TwoJets,0,1, vWght[i]);
		}
	}
*/

// I am not sure what this is. So I am commenting it out for now to study EWK component in data and MC - 03-16-2010
	//for EWK MC only. weghts == 1, then noramize to data
/*	int nJet15 = 0;
	std::vector<int> JetInd;
	for (int i=0; i < UsedEle.size(); ++i) {
		nJet15 = 0;
		JetInd.clear();
		for (int j=0; j < stuple.jet_num; j++) {
				
			int iJetIndex = j;
			int iEleIndex = UsedEle[i];
				
			TLorentzVector tlJet(stuple.jet_Px[iJetIndex], stuple.jet_Py[iJetIndex], stuple.jet_Pz[iJetIndex], stuple.jet_E[iJetIndex]);
			TLorentzVector tlEle(stuple.ele_Px[iEleIndex], stuple.ele_Py[iEleIndex], stuple.ele_Pz[iEleIndex], stuple.ele_E[iEleIndex]);
	
			TLorentzVector tlSum = tlEle + tlJet; 
			float DelR = tlEle.DeltaR(tlJet);
			float etratio = stuple.jet_Pt[j]/stuple.ele_Etc[i];
			delREt->Fill(DelR, etratio);
			if (DelR < 0.1)  continue;
			nJet15++;
			JetInd.push_back(j);
		}

		const int obj_collection = 0; //pick central values of all objects
		//if (stuple.jet_NJet15 >= 1) {	// mono-jet case
		if (nJet15 >= 1) {	// mono-jet case
			iEle1JetEvts++;
			fSFR_1j += 1;
			fSFR_m_1j += vWght_m[i];
			fSFR_p_1j += vWght_p[i];
			myHistMan.FillElectronHists(stuple,obj_collection, GA_Pho,UsedEle[i], 1);
			myHistMan.FillEventHists(stuple, GA_Evt, 1);
			myHistMan.FillJetHists(stuple,obj_collection, GA_Jet,JetInd[0], 1);		//lead jet
			myHistMan.FillElectron1JetHists(stuple,obj_collection, GA_PhoJet,UsedEle[i],JetInd[0], 1);
		}
		
		//if (stuple.jet_NJet15 >= 2) {	// di-jet case
		if (nJet15 >= 2) {	// di-jet case
			iEle2JetEvts++;
			fSFR_2j += 1;
			fSFR_m_2j += vWght_m[i];
			fSFR_p_2j += vWght_p[i];
			myHistMan.FillElectronHists(stuple,obj_collection, GB_Pho,UsedEle[i], 1);
			myHistMan.FillEventHists(stuple, GB_Evt, 1);
			myHistMan.FillJetHists(stuple,obj_collection, GB_Jet1,JetInd[0], 1);		//lead jet
			myHistMan.FillJetHists(stuple,obj_collection, GB_Jet2,JetInd[1], 1);		//2nd lead jet
			myHistMan.FillElectron1JetHists(stuple,obj_collection, GB_PhoJet1,UsedEle[i],JetInd[0], 1);
			myHistMan.FillElectron1JetHists(stuple,obj_collection, GB_PhoJet2,UsedEle[i],JetInd[1], 1);
			myHistMan.FillElectron2JetsHists(stuple,obj_collection, GB_PhoJets,UsedEle[i],JetInd[0],JetInd[1], 1);
			myHistMan.FillTwoJetsHists(stuple,obj_collection, GB_TwoJets,JetInd[0],JetInd[1], 1);
		}
	}
*/




	}

//-------------------------------------------------------------
void EleJetsTemp::Init(TChain *tree)
{
	iEventsSelected = 0;
	iEle1JetEvts = 0;
	iEle2JetEvts = 0;
	iEvtsDelR2 = 0;

	fSFR_1j =0; fSFR_m_1j=0; fSFR_p_1j=0;			// sum of fake rates
	fSFR_2j=0; fSFR_m_2j=0; fSFR_p_2j=0;			// sum of fake rates

	
	readIn = new ReadInAll(tree,&stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->Met(1);
/*	if (bMcFlag) {
		readIn->EmUpPhotons(1);
		readIn->EmUpElectrons(1);
		readIn->EmDownPhotons(1);
		readIn->EmDownElectrons(1);
		readIn->JesUpJets(1);
		readIn->JesDownJets(1);
		readIn->GenLevel(0);
	}
*/
	readIn->Init();

	//cannot add folder without creating a dir first!
	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = rootFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

	TDatime dt;
	std::string sTime(dt.AsString());
	sVec.push_back(" Job Summary for EleJetsTemp : " + sTime);
	sVec.push_back("Data Written to        = " + GetHistFileName());
	sVec.push_back("MC Flag setting        = " + ToStr(bMcFlag) + " (0==DATA, 1==MC)");
	std::stringstream m1, m2, m3, m4;
	m1 << "min,max Vtx for Signal = " << std::setw(4) << GetMinClass12Vtx() << ", " << std::setw(4) << GetMaxClass12Vtx();
	m2 << "Ele Et, Eta            = " << std::setw(4) << GetMinEleEt() << ", " << std::setw(4)<< GetMaxEleEta();
	m3 << "Signal Time Window     = " << std::setw(4) << GetSignalTimeLoLimit() << ", " << std::setw(4) << GetSignalTimeUpLimit();
	m4 << "Min. MEt               = " << std::setw(4) << GetMinMet();
	sVec.push_back(m1.str());
	sVec.push_back(m2.str());
	sVec.push_back(m3.str());
	sVec.push_back(m4.str());
	sVec.print();

	

} // Init

//-------------------------------------------------------------------
void EleJetsTemp::CleanUp()
{
	rootFile->Write();
	rootFile->Close();
	sVec.clear();
} 

//-------------------------------------------------------------------
void EleJetsTemp::BookHistograms()
{
	Histograms HistoMan;
	std::ostringstream foldername, text2title1, text2title2;
	text2title1 << "e+>=1Jets :";
	text2title2 << "e+>=2Jets :";
	

	TDirectory *sub, *subsub;

		topDir->cd();

	delREt = new TH2F("DELRET"," dELeT RATIO" , 50,0,5,50,0,5);
		
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
void EleJetsTemp::DoDump(Stuple& stuple)
{
	
	if (stuple.ele_num <1) return;
	if (stuple.pho_num <1) return;
	
	
	std::cout << "============= " << stuple.evt_RunNumber << ", " << stuple.evt_EventNumber << std::endl;
	std::cout << "elenum ,jetnum = " << stuple.ele_num << ", " << stuple.jet_num << std::endl;
	std::cout << "---- ELE INFO " << std::endl;
	if (stuple.ele_num > 0 && stuple.jet_num > 0) {
		for (int i =0 ; i < stuple.ele_num ; i++ ) {
			//if (stuple.ele_TightId[i] != 0) continue;
			std::cout << i << "\t EleEt=" << stuple.ele_Etc[i] << "\tDetEta=" << stuple.ele_DetEta[i] << "\tDetPhi=" << stuple.ele_DetPhi[i] << std::endl;
			for (int j=0; j < stuple.jet_num; j++) {
				
				int iJetIndex = j;
				int iEleIndex = i;
				
				TLorentzVector tlJet(stuple.jet_Px[iJetIndex], stuple.jet_Py[iJetIndex], stuple.jet_Pz[iJetIndex], stuple.jet_E[iJetIndex]);
				TLorentzVector tlEle(stuple.ele_Px[iEleIndex], stuple.ele_Py[iEleIndex], stuple.ele_Pz[iEleIndex], stuple.ele_E[iEleIndex]);
	
				TLorentzVector tlSum = tlEle + tlJet; 
				float DelR = tlEle.DeltaR(tlJet);
				float etratio = stuple.jet_Pt[j]/stuple.ele_Etc[i];
				delREt->Fill(DelR, etratio);
				
				if (DelR < 0.1) { iEvtsDelR2++;
				std::cout << j << "\t JetEt=" << stuple.jet_Pt[j] << "\tDetEta " << stuple.jet_DetEta[j] << "\tDetPhi=" << stuple.jet_DetPhi[i] << "\tDelR= " << DelR<< std::endl;
				
						std::cout << "---- PHO INFO " << std::endl;
						std::cout << "phonum ,jetnum = " << stuple.pho_num << ", " << stuple.jet_num << std::endl;
						for (int i =0 ; i < stuple.pho_num ; i++ ) {
							std::cout << i << "\t PhoEt=" << stuple.pho_Etc[i] << "\tDetEta=" << stuple.pho_DetEta[i] << "\tDetPhi=" << stuple.pho_DetPhi[i] << std::endl;
						}
				}
			}
		}

	}
	
}

//------ main -------------------------------------------------------
void EleJetsTemp::Main(TChain *ch, int iRunEvents)
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

	//do a check for the MC. this is a common mistake I make
	
	  readIn->GetEntry(1);
	  if (stuple.evt_McFlag != GetMcFlag()) 
	  {
		  std::cout << red << "MC Flag," << GetMcFlag()
			  << ", setting does not match what is in Stuple, " 
			  << stuple.evt_McFlag << ". pleas check." 
			  << " returning!" << clearatt << std::endl;
		  //std::cout << "Changing the MCFlag to Stuple setting McFlag = " << stuple.evt_McFlag << clearatt << std::endl;
		  //SetMcFlag(stuple.evt_McFlag); //for some reason this seems to overwrite counting arrays? 12-07-2009
		  //return; // simple return causes trouble when rerun. so drastic measure needed, even exit(1) would not work
		  assert (false);
	  }

	unsigned int iEvtProc = 0;		// number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iFirst400 =0;

	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) {
	  	readIn->GetEntry(iEvtProc);

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << iEvtProc << "\t events processed." << std::endl;
		}
		
		if (stuple.evt_McFlag == 0) {  //______________________________________________ if data
			if (stuple.evt_RunNumber < 190851) {
				iFirst400++;
				continue;
			}
		}
		if (stuple.evt_RunNumber >246231) continue;  //drop the newly added data, 0k dataset in StupleV4

	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;

		DoMyStuff(stuple);
		//DoDump(stuple);
		
	};
	
	sVec.push_back("__________________________________________");
	sVec.push_back("Entries found    = " + ToStr(iENTRIES));
	sVec.push_back("Entries request  = " + ToStr(iRunEvents));
	sVec.push_back("Events Processed = " + ToStr(iEvtProc));
	sVec.push_back("Ele./Pho Events  = " + ToStr(iPhoEleEvts));
	sVec.push_back("Ele Events       = " + ToStr(iEleEvts));
	sVec.push_back("Pho Events       = " + ToStr(iPhoEvts));
	sVec.push_back("Other Events     = " + ToStr(iOther));
	sVec.push_back("First400 Events  = " + ToStr(iFirst400));
	sVec.push_back("sum Events       = " + ToStr(iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400));
	sVec.push_back("__________________________________________");
	sVec.push_back("Ele Evts         = " + ToStr(iEventsSelected));
	sVec.push_back("Ele+>=1Jet Evts  = " + ToStr(iEle1JetEvts)); 
	sVec.push_back("Ele+>=2Jet Evts  = " + ToStr(iEle2JetEvts)); 
	sVec.push_back("Sum FR_d (1Jet)  = " + ToStr(fSFR_1j) + " +/- " + ToStr(TMath::Max(fabs(fSFR_m_1j-fSFR_1j), fabs(fSFR_p_1j-fSFR_1j))));
	sVec.push_back("Sum FR_d (2Jet)  = " + ToStr(fSFR_2j) + " +/- " + ToStr(TMath::Max(fabs(fSFR_m_2j-fSFR_2j), fabs(fSFR_p_2j -fSFR_2j))));
	sVec.push_back("Data Written to  = " + GetHistFileName());
	sVec.push_back("DelR<0.2 evts    = " + ToStr(iEvtsDelR2)); 
	sVec.push_back("======= END EleJetsTemp ===================");
	
	sVec.print();
	topDir->cd();
	(sVec.writeToCanvas())->Write();
	delete sVec.canvas();

	CleanUp();
} // Main

