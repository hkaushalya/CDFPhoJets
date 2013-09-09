/*{{{*/
/*  $Id: PhoEleStudy.cc,v 1.1 2013/02/28 03:37:30 samantha Exp $
 *  $Log: PhoEleStudy.cc,v $
 *  Revision 1.1  2013/02/28 03:37:30  samantha
 *  Final commit. no checks.! these were never commited to cvs before!
 *
 */

/*}}}*/

#include <iostream>
#include <sstream>
#include "PhoFlat/PhoEleStudy.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
#include "TBenchmark.h"
#include "TCanvas.h"
#include <iomanip>
#include "TROOT.h"
//for debugging
#include <exception>
#include <stdexcept> // out_of_range exception

void DEBUG(const std::string func, const double line, const std::string msg)
{
	std::cout << func << ":" << line << ": " << 	msg << std::endl;
}

//-------------------------------------------------------------------
PhoEleStudy::PhoEleStudy():
iProgressBy(50000),
sHistFileName("PhoEleStudy.root"),
bMcFlag(true)
{

}

//-------------------------------------------------------------
void PhoEleStudy::Init(TChain *tree)
{
	readIn = new ReadInAll(tree, &stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	//readIn->CentralJets(1);
	readIn->Met(1);

		for (int i=0; i<4; i++) {
			vSignalCount.push_back(0);
		}
	
	readIn->Init();

	//cannot add folder without creating a dir first!
	histFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = histFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

	
} // Init

//-------------------------------------------------------------------
void PhoEleStudy::CleanUp()
{
	vSignalCount.clear();

	histFile->Write();
	histFile->Close();
} 


//------ main -------------------------------------------------------
void PhoEleStudy::Main(TChain *ch, int iRunEvents)
{
	TBenchmark timer;
	timer.Start("phoana_time");
	
	
	if (ch) {
		myChain = ch;
		//ch->Print();	
	} else {
		std::cout << "NULL chain!" << std::endl;
		assert(false);
	}
	if (ch->GetEntries()<1)
	{
		std::cout << "No entries found!" << std::endl;
		return;
	}

	gROOT->Reset();
	Init(myChain);
	
	BookHistograms();
	
	unsigned int iENTRIES = (unsigned int) myChain->GetEntries();
	unsigned int iEvt2Process = 0; 	//_____________________________________________ events to be processed
	if (iRunEvents < 0 ) {
		iRunEvents = iENTRIES;			//______________________________________________ requested all entries
		iEvt2Process = iENTRIES;
	} else if (iRunEvents > 0 && iRunEvents <= iENTRIES) iEvt2Process = iRunEvents;
	else iEvt2Process = iENTRIES;		//______________________________________________ default

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

	unsigned int iEvtProc = 0;			//______________________________________________ number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iFirst400 =0;
	std::cout << "Init done. Ready, Set , GO... " << std::endl; 
	
	timer.Start("looptimer");
	
	unsigned tightPhoEle=0, sidePhoEle=0;
	iCurrTree = -1;
	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) 
	{
		if (iCurrTree != myChain->GetTreeNumber())
		{
			std::cout << green << "Opening file " << myChain->GetFile()->GetName() << clearatt << std::endl;
			iCurrTree = myChain->GetTreeNumber();
		}
		readIn->GetEntry(iEvtProc);

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << iEvtProc << "\t events processed [" << (int)( iEvtProc/(double)iEvt2Process * 100)<< "%] ";
			timer.Show("looptimer");
			timer.Start("looptimer");
		}

		//___________________________________________________________________________ drop the first 400pb-1
		if (stuple.evt_McFlag == 0) {  //______________________________________________ if data
			if (stuple.evt_RunNumber < 190851) {
				iFirst400++;
				continue;
			}
		}

		//for halo id eff only 3-13-2010
		//std::cout << "pho_num status=" << myChain->GetBranchStatus("pho_num") << std::endl;
		if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
		else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
		else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;

		CommonVars cVars;
		PhotonList cPhos;
		ElectronList cEles;
		//JetList centralJets, jesupJets, jesdownJets;

		GetCommonVariables(stuple, cVars);
		CreatePhotonLists(stuple, cPhos,0);
		CreateElectronLists(stuple, cEles,0);
		//RemoveDuplicateEMobjects(centralPhos,centralEles); 
		//CreateJetLists(stuple, centralJets,0);
		//centralJets.Dump();
		//centralPhos.Dump();
		//centralEles.Dump();

		//_________________________________________________________________ require a one good vertex 		  
		if (cVars.vtx_NClass12 <= 0) continue;
		if (cVars.vtx_z > 60.)  continue;


		//______________________________________________________________ trigger
		if (  (cVars.tri_pho25iso != 1) && (cVars.tri_pho50 != 1) 	//_______________ only data has L1-3 trig. info
				&& (cVars.tri_pho70 != 1) ) 										//_______________ mc has L1-2.
			continue;												

		//______________________________________________________________ select the photon
		for (unsigned int i = 0; i < cPhos.pho_num; ++i )
		{

			if (cPhos.pho_Etc[i] <10) continue;
			if (fabs(cPhos.pho_DetEta.at(i))>1.1) continue;
			const TLorentzVector tlPho(cPhos.pho_Px.at(i), cPhos.pho_Py.at(i), cPhos.pho_Pz.at(i), cPhos.pho_E.at(i));
			if (cPhos.pho_TightId[i] == 0)
			{
				for (unsigned int i = 0; i < cEles.ele_num; ++i )
				{
					const TLorentzVector tlEle(cEles.ele_Px.at(i), cEles.ele_Py.at(i), cEles.ele_Pz.at(i), cEles.ele_E.at(i));
					if (tlPho.DeltaR(tlEle)<0.3)
					{
						//if (cEles.ele_TightId[i] == 0) 
						++tightPhoEle;
						break;
					}
				}

			} 
			if (cPhos.pho_TightId[i] > 0 && cPhos.pho_LooseId[i] == 0)
			{
				for (unsigned int i = 0; i < cEles.ele_num; ++i )
				{
					const TLorentzVector tlEle(cEles.ele_Px.at(i), cEles.ele_Py.at(i), cEles.ele_Pz.at(i), cEles.ele_E.at(i));
					if (tlPho.DeltaR(tlEle)<0.3)
					{
						if (cEles.ele_TightId[i] == 0) ++sidePhoEle;
						break;
					}
				}
			}
		}	



	}
	
	//print out the job summary
	std::cout << magenta << "======== PhoEleStudy =======================" << std::endl;
	std::cout << "Entries found    = " << iENTRIES << std::endl;
	std::cout << "Entries request  = " << iRunEvents << " (" << iRunEvents/(double) iENTRIES * 100 << "%)"<< std::endl;
	std::cout << "Events Processed = " << iEvtProc << " (" << iEvtProc/(double) iENTRIES * 100 << "%)"<< std::endl;
	std::cout << "Ele./Pho Events  = " << iPhoEleEvts << std::endl;
	std::cout << "Ele Events       = " << iEleEvts << std::endl;
	std::cout << "Pho Events       = " << iPhoEvts << std::endl;
	std::cout << "Other Events     = " << iOther << std::endl;
	std::cout << "First400 Events  = " << iFirst400 << std::endl;
	std::cout << "sum Events       = " << iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400 << std::endl;
	std::cout << red << "__________________________________________" << std::endl;

		const int itextWidth = 20;
		std::cout << std::setw(itextWidth) << "tight pho and tight ele = " << tightPhoEle << std::endl;
		std::cout << std::setw(itextWidth) << "sidepho and tight ele   = " << sidePhoEle << std::endl; 
	std::cout << cyan << "Data Written to        = " << GetHistFileName() << clearatt << std::endl;
	std::cout << "======= END PhoEleStudy ====================" << std::endl;

	CleanUp();
	timer.Show("phoana_time");
} // Main


//-------------------------------------------------------------------
void PhoEleStudy::BookHistograms()
{
	
} //BookHistograms



void PhoEleStudy::CreatePhotonLists(const Stuple stuple, PhotonList& pholist, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
			pholist.pho_num 	 =  stuple.pho_num;
			pholist.pho_Ntight = stuple.pho_Ntight;
			pholist.pho_Nloose = stuple.pho_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_E[i]);
				pholist.pho_Px.push_back(stuple.pho_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_matchJetIndex[i]);
				
				pholist.pho_CprWgt.push_back(stuple.pho_CprWgt[i]);
				pholist.pho_CprSys1.push_back(stuple.pho_CprSys1[i]);
				pholist.pho_CprSys2.push_back(stuple.pho_CprSys2[i]);
				pholist.pho_CprSys3.push_back(stuple.pho_CprSys3[i]);
				pholist.pho_CprSys4.push_back(stuple.pho_CprSys4[i]);
				pholist.pho_CprSys5.push_back(stuple.pho_CprSys5[i]);
				pholist.pho_CprSys6.push_back(stuple.pho_CprSys6[i]);
				pholist.pho_CprSys7.push_back(stuple.pho_CprSys7[i]);
				pholist.pho_CprSys8.push_back(stuple.pho_CprSys8[i]);
				
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		case 1:
		{
			pholist.pho_num 	 =  stuple.pho_up_num;
			pholist.pho_Ntight = stuple.pho_up_Ntight;
			pholist.pho_Nloose = stuple.pho_up_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_up_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_up_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_up_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_up_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_up_E[i]);
				pholist.pho_Px.push_back(stuple.pho_up_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_up_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_up_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_up_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_up_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_up_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_up_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_up_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_up_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_up_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_up_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_up_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_up_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_up_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_up_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_up_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_up_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_up_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_up_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_up_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_up_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_up_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_up_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_up_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_up_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_up_matchJetIndex[i]);
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		case -1:
		{
			pholist.pho_num 	 =  stuple.pho_down_num;
			pholist.pho_Ntight = stuple.pho_down_Ntight;
			pholist.pho_Nloose = stuple.pho_down_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_down_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_down_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_down_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_down_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_down_E[i]);
				pholist.pho_Px.push_back(stuple.pho_down_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_down_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_down_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_down_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_down_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_down_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_down_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_down_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_down_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_down_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_down_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_down_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_down_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_down_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_down_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_down_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_down_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_down_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_down_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_down_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_down_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_down_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_down_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_down_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_down_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_down_matchJetIndex[i]);
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid photon collection requested! can make central/em up/em down photon lists only!");
			assert(false);

		}
		
	}	//switch

/*}}}*/
}	//CreatePhotonLists



void PhoEleStudy::CreateElectronLists(const Stuple stuple, ElectronList& list, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
			list.ele_num    = stuple.ele_num;
			list.ele_Ntight = stuple.ele_Ntight;
			list.ele_Nloose = stuple.ele_Nloose;

			for (unsigned int i=0; i < stuple.ele_num; i++)
			{
				list.ele_Index.push_back(stuple.ele_Index[i]);
				list.ele_PhoBlockIndex.push_back(stuple.ele_PhoBlockIndex[i]);
				list.ele_EleBlockIndex.push_back(stuple.ele_EleBlockIndex[i]);
				list.ele_Etc.push_back(stuple.ele_Etc[i]);
				list.ele_E.push_back(stuple.ele_E[i]);
				list.ele_Px.push_back(stuple.ele_Px[i]);
				list.ele_Py.push_back(stuple.ele_Py[i]);
				list.ele_Pz.push_back(stuple.ele_Pz[i]);
				list.ele_Detector.push_back(stuple.ele_Detector[i]);
				list.ele_DetEta.push_back(stuple.ele_DetEta[i]);
				list.ele_DetPhi.push_back(stuple.ele_DetPhi[i]);
				list.ele_XCes.push_back(stuple.ele_XCes[i]);
				list.ele_ZCes.push_back(stuple.ele_ZCes[i]);
				list.ele_HadEm.push_back(stuple.ele_HadEm[i]);
				list.ele_Chi2Mean.push_back(stuple.ele_Chi2Mean[i]);
				list.ele_N3d.push_back(stuple.ele_N3d[i]);
				list.ele_Iso4.push_back(stuple.ele_Iso4[i]);
				list.ele_TrkIso.push_back(stuple.ele_TrkIso[i]);		
				list.ele_CesWireE2.push_back(stuple.ele_CesWireE2[i]);
				list.ele_CesStripE2.push_back(stuple.ele_CesStripE2[i]);
				list.ele_PhiWedge.push_back(stuple.ele_PhiWedge[i]);
				list.ele_NMuonStubs.push_back(stuple.ele_NMuonStubs[i]);
				list.ele_EmTime.push_back(stuple.ele_EmTime[i]);
				list.ele_PhoenixId.push_back(stuple.ele_PhoenixId[i]);
				list.ele_Halo_seedWedge.push_back(stuple.ele_Halo_seedWedge[i]);
				list.ele_Halo_eastNhad.push_back(stuple.ele_Halo_eastNhad[i]);
				list.ele_Halo_westNhad.push_back(stuple.ele_Halo_westNhad[i]);
				list.ele_matchJetIndex.push_back(stuple.ele_matchJetIndex[i]);
				list.ele_Ntracks.push_back(stuple.ele_Ntracks[i]);	
				list.ele_Emfr.push_back(stuple.ele_Emfr[i]);
				list.ele_EoverP.push_back(stuple.ele_EoverP[i]);
				list.ele_TrackPt.push_back(stuple.ele_TrackPt[i]);
				list.ele_TrackBcPt.push_back(stuple.ele_TrackBcPt[i]);
				list.ele_TrackPhi.push_back(stuple.ele_TrackPhi[i]);
				list.ele_Nssl.push_back(stuple.ele_Nssl[i]);
				list.ele_Nasl.push_back(stuple.ele_Nasl[i]);
				list.ele_TightId.push_back(stuple.ele_TightId[i]);
				list.ele_LooseId.push_back(stuple.ele_LooseId[i]);
				list.ele_ConversionId.push_back(stuple.ele_ConversionId[i]);
			}

			assert( (int)list.ele_num == (int)list.ele_Index.size());
			break;
		}

		case 1:
		{
			list.ele_num    = stuple.ele_up_num;
			list.ele_Ntight = stuple.ele_up_Ntight;
			list.ele_Nloose = stuple.ele_up_Nloose;

			for (unsigned int i=0; i < stuple.ele_up_num; i++)
			{
				list.ele_Index.push_back(stuple.ele_up_Index[i]);
				list.ele_PhoBlockIndex.push_back(stuple.ele_up_PhoBlockIndex[i]);
				list.ele_EleBlockIndex.push_back(stuple.ele_up_EleBlockIndex[i]);
				list.ele_Etc.push_back(stuple.ele_up_Etc[i]);
				list.ele_E.push_back(stuple.ele_up_E[i]);
				list.ele_Px.push_back(stuple.ele_up_Px[i]);
				list.ele_Py.push_back(stuple.ele_up_Py[i]);
				list.ele_Pz.push_back(stuple.ele_up_Pz[i]);
				list.ele_Detector.push_back(stuple.ele_up_Detector[i]);
				list.ele_DetEta.push_back(stuple.ele_up_DetEta[i]);
				list.ele_DetPhi.push_back(stuple.ele_up_DetPhi[i]);
				list.ele_XCes.push_back(stuple.ele_up_XCes[i]);
				list.ele_ZCes.push_back(stuple.ele_up_ZCes[i]);
				list.ele_HadEm.push_back(stuple.ele_up_HadEm[i]);
				list.ele_Chi2Mean.push_back(stuple.ele_up_Chi2Mean[i]);
				list.ele_N3d.push_back(stuple.ele_up_N3d[i]);
				list.ele_Iso4.push_back(stuple.ele_up_Iso4[i]);
				list.ele_TrkIso.push_back(stuple.ele_up_TrkIso[i]);		
				list.ele_CesWireE2.push_back(stuple.ele_up_CesWireE2[i]);
				list.ele_CesStripE2.push_back(stuple.ele_up_CesStripE2[i]);
				list.ele_PhiWedge.push_back(stuple.ele_up_PhiWedge[i]);
				list.ele_NMuonStubs.push_back(stuple.ele_up_NMuonStubs[i]);
				list.ele_EmTime.push_back(stuple.ele_up_EmTime[i]);
				list.ele_PhoenixId.push_back(stuple.ele_up_PhoenixId[i]);
				list.ele_Halo_seedWedge.push_back(stuple.ele_up_Halo_seedWedge[i]);
				list.ele_Halo_eastNhad.push_back(stuple.ele_up_Halo_eastNhad[i]);
				list.ele_Halo_westNhad.push_back(stuple.ele_up_Halo_westNhad[i]);
				list.ele_matchJetIndex.push_back(stuple.ele_up_matchJetIndex[i]);
				list.ele_Ntracks.push_back(stuple.ele_up_Ntracks[i]);	
				list.ele_Emfr.push_back(stuple.ele_up_Emfr[i]);
				list.ele_EoverP.push_back(stuple.ele_up_EoverP[i]);
				list.ele_TrackPt.push_back(stuple.ele_up_TrackPt[i]);
				list.ele_TrackBcPt.push_back(stuple.ele_up_TrackBcPt[i]);
				list.ele_TrackPhi.push_back(stuple.ele_up_TrackPhi[i]);
				list.ele_Nssl.push_back(stuple.ele_up_Nssl[i]);
				list.ele_Nasl.push_back(stuple.ele_up_Nasl[i]);
				list.ele_TightId.push_back(stuple.ele_up_TightId[i]);
				list.ele_LooseId.push_back(stuple.ele_up_LooseId[i]);
				list.ele_ConversionId.push_back(stuple.ele_up_ConversionId[i]);
			}

			assert(list.ele_num == list.ele_Index.size());
			break;

		}

		case -1:
		{
			list.ele_num    = stuple.ele_down_num;
			list.ele_Ntight = stuple.ele_down_Ntight;
			list.ele_Nloose = stuple.ele_down_Nloose;

			for (unsigned int i=0; i < stuple.ele_down_num; i++)
			{
				list.ele_Index.push_back(stuple.ele_down_Index[i]);
				list.ele_PhoBlockIndex.push_back(stuple.ele_down_PhoBlockIndex[i]);
				list.ele_EleBlockIndex.push_back(stuple.ele_down_EleBlockIndex[i]);
				list.ele_Etc.push_back(stuple.ele_down_Etc[i]);
				list.ele_E.push_back(stuple.ele_down_E[i]);
				list.ele_Px.push_back(stuple.ele_down_Px[i]);
				list.ele_Py.push_back(stuple.ele_down_Py[i]);
				list.ele_Pz.push_back(stuple.ele_down_Pz[i]);
				list.ele_Detector.push_back(stuple.ele_down_Detector[i]);
				list.ele_DetEta.push_back(stuple.ele_down_DetEta[i]);
				list.ele_DetPhi.push_back(stuple.ele_down_DetPhi[i]);
				list.ele_XCes.push_back(stuple.ele_down_XCes[i]);
				list.ele_ZCes.push_back(stuple.ele_down_ZCes[i]);
				list.ele_HadEm.push_back(stuple.ele_down_HadEm[i]);
				list.ele_Chi2Mean.push_back(stuple.ele_down_Chi2Mean[i]);
				list.ele_N3d.push_back(stuple.ele_down_N3d[i]);
				list.ele_Iso4.push_back(stuple.ele_down_Iso4[i]);
				list.ele_TrkIso.push_back(stuple.ele_down_TrkIso[i]);		
				list.ele_CesWireE2.push_back(stuple.ele_down_CesWireE2[i]);
				list.ele_CesStripE2.push_back(stuple.ele_down_CesStripE2[i]);
				list.ele_PhiWedge.push_back(stuple.ele_down_PhiWedge[i]);
				list.ele_NMuonStubs.push_back(stuple.ele_down_NMuonStubs[i]);
				list.ele_EmTime.push_back(stuple.ele_down_EmTime[i]);
				list.ele_PhoenixId.push_back(stuple.ele_down_PhoenixId[i]);
				list.ele_Halo_seedWedge.push_back(stuple.ele_down_Halo_seedWedge[i]);
				list.ele_Halo_eastNhad.push_back(stuple.ele_down_Halo_eastNhad[i]);
				list.ele_Halo_westNhad.push_back(stuple.ele_down_Halo_westNhad[i]);
				list.ele_matchJetIndex.push_back(stuple.ele_down_matchJetIndex[i]);
				list.ele_Ntracks.push_back(stuple.ele_down_Ntracks[i]);	
				list.ele_Emfr.push_back(stuple.ele_down_Emfr[i]);
				list.ele_EoverP.push_back(stuple.ele_down_EoverP[i]);
				list.ele_TrackPt.push_back(stuple.ele_down_TrackPt[i]);
				list.ele_TrackBcPt.push_back(stuple.ele_down_TrackBcPt[i]);
				list.ele_TrackPhi.push_back(stuple.ele_down_TrackPhi[i]);
				list.ele_Nssl.push_back(stuple.ele_down_Nssl[i]);
				list.ele_Nasl.push_back(stuple.ele_down_Nasl[i]);
				list.ele_TightId.push_back(stuple.ele_down_TightId[i]);
				list.ele_LooseId.push_back(stuple.ele_down_LooseId[i]);
				list.ele_ConversionId.push_back(stuple.ele_down_ConversionId[i]);
			}

			assert(list.ele_num == list.ele_Index.size());
			break;
		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid collection requested! can make central/em up/em down  lists only!");
			assert(false);
		}
	}
/*}}}*/
}	//CreateElectronLists

void PhoEleStudy::CreateJetLists(const Stuple stuple, JetList& list, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
		  list.jet_num  = stuple.jet_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_NJet15;
			for (unsigned int i=0; i < stuple.jet_num; i++)
			{
				if ( fabs(stuple.jet_Px[i]) > 1500 || 
					  fabs(stuple.jet_Py[i]) > 1500 || 
					  fabs(stuple.jet_Pz[i]) > 1500 || 
					  fabs(stuple.jet_E[i]) > 1500)
				{
					stuple.DumpJetBlock();
					assert(false);
				}

				list.jet_Index.push_back(stuple.jet_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_E[i]);
				list.jet_Px.push_back(stuple.jet_Px[i]);
				list.jet_Py.push_back(stuple.jet_Py[i]);
				list.jet_Pz.push_back(stuple.jet_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_SeedIEta[i]);
			}

			assert(list.jet_num == list.jet_Index.size());

			list.jet_raw_num  = stuple.jet_raw_num;							// number of raw jets
			for (unsigned int i=0; i < stuple.jet_raw_num; i++)
			{
				list.jet_raw_Index.push_back(stuple.jet_raw_Index[i]);				//original index in jet block
				list.jet_raw_Pt.push_back(stuple.jet_raw_Pt[i]);
				list.jet_raw_E.push_back(stuple.jet_raw_E[i]);
				list.jet_raw_Px.push_back(stuple.jet_raw_Px[i]);
				list.jet_raw_Py.push_back(stuple.jet_raw_Py[i]);
				list.jet_raw_Pz.push_back(stuple.jet_raw_Pz[i]);
			}
			
			break;
		}

		case 1:
		{
		  list.jet_num  = stuple.jet_up_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_up_NJet15;
			for (unsigned int i=0; i < stuple.jet_up_num; i++)
			{
				list.jet_Index.push_back(stuple.jet_up_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_up_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_up_E[i]);
				list.jet_Px.push_back(stuple.jet_up_Px[i]);
				list.jet_Py.push_back(stuple.jet_up_Py[i]);
				list.jet_Pz.push_back(stuple.jet_up_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_up_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_up_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_up_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_up_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_up_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_up_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_up_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_up_SeedIEta[i]);
				
			}

			assert(list.jet_num == list.jet_Index.size());
			break;

			
		}

		case -1:
		{
		  list.jet_num  = stuple.jet_down_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_down_NJet15;
			for (unsigned int i=0; i < stuple.jet_down_num; i++)
			{
				list.jet_Index.push_back(stuple.jet_down_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_down_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_down_E[i]);
				list.jet_Px.push_back(stuple.jet_down_Px[i]);
				list.jet_Py.push_back(stuple.jet_down_Py[i]);
				list.jet_Pz.push_back(stuple.jet_down_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_down_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_down_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_down_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_down_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_down_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_down_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_down_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_down_SeedIEta[i]);
				
			}

			assert(list.jet_num == list.jet_Index.size());
			break;

		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid collection requested! can make central and JES up/down lists only!");
			assert(false);

		}

	}
/*}}}*/
} // CreateJetLists

void PhoEleStudy::GetCommonVariables(const Stuple stuple, CommonVars& cVars)
{
/*{{{*/
	cVars.evt_McFlag 		= stuple.evt_McFlag;
	cVars.evt_RunNumber 	= stuple.evt_RunNumber;
	cVars.evt_EventNumber= stuple.evt_EventNumber;
	cVars.tri_pho25iso 	= stuple.tri_pho25iso;
	cVars.tri_pho50 		= stuple.tri_pho50;
	cVars.tri_pho70 		= stuple.tri_pho70;

	cVars.vtx_N 			= stuple.vtx_N;
	cVars.vtx_NClass12 	= stuple.vtx_NClass12;
	cVars.vtx_z 			= stuple.vtx_z;
	cVars.vtx_Ntracks 	= stuple.vtx_Ntracks;
	cVars.vtx_SumPt 		= stuple.vtx_SumPt;

	cVars.met_Met 			= stuple.met_Met;	
	cVars.met_MetX 		= stuple.met_MetX;	
	cVars.met_MetY 		= stuple.met_MetY;	
	cVars.met_SumEt 		= stuple.met_SumEt;
	cVars.met_Ht 			= stuple.met_Ht;	
	cVars.met_MetPhi 		= stuple.met_MetPhi;
//	cVars.met_Gen_d 		= stuple.met_Gen_d;	
//	cVars.met_Gen_m 		= stuple.met_Gen_m;
//	cVars.met_Gen_p 		= stuple.met_Gen_p;
//	cVars.met_Gen_mUn 	= stuple.met_Gen_mUn;
//	cVars.met_Gen_pUn 	= stuple.met_Gen_pUn;
/*}}}*/
}

