/*{{{*/
/*  $Id: ZSel.cc,v 1.1 2011/05/31 17:57:49 samantha Exp $
 *  $Log: ZSel.cc,v $
 *  Revision 1.1  2011/05/31 17:57:49  samantha
 *  This is to get corrections to EWK MC by comparing Z data and Z MC.
 *
 */

/*}}}*/

#include <iostream>
#include <sstream>
#include "PhoFlat/ZSel.hh"
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
ZSel::ZSel():
iProgressBy(50000),
sHistFileName("ZSel.root"),
bMcFlag(true),
bUseNvtxWgts(false)
{

	fMinEleEt = 15.; // GeV  //do not lower this any more. See header file Set method for more details!
	fMaxEleEt = 1500.; // GeV
	iMinClass12Vtx = 1;
	iMaxClass12Vtx = 100;
	fMaxVtxz = 60.; //cm
	fMinMet = 0.0;	//minimum Missing transverse energy
	fMaxMet = 1500.0;	//maximum Missing transverse energy

}

//-------------------------------------------------------------
void ZSel::Init(TChain *tree)
{
	iExtraEle = 0;

	readIn = new ReadInAll(tree, &stuple);
	readIn->CentralElectrons(1);
	readIn->Met(1);

	for (int i=0; i<4; i++) {
		vMcCount.push_back(0);
	}

	readIn->Init();

	//cannot add folder without creating a dir first!
	histFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = histFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

	vNvtxWgts.clear();
	if (UseNvtxWeights())
	{
		vNvtxWgts.push_back(0.772324);
		vNvtxWgts.push_back(1.04891);
		vNvtxWgts.push_back(1.33888);
		vNvtxWgts.push_back(1.71057);
		vNvtxWgts.push_back(2.63263);
		vNvtxWgts.push_back(5.95169);
		vNvtxWgts.push_back(13.9581);
		vNvtxWgts.push_back(32.4681);
	}

} // Init

//-------------------------------------------------------------------
void ZSel::CleanUp()
{
	vMcCount.clear();

	histFile->Write();
	histFile->Close();
} 


//------ main -------------------------------------------------------
void ZSel::Main(TChain *ch, int iRunEvents)
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
	std::cout << yellow << "******* Job Settings ***********" << std::endl;
	PrintJobSettings();
	std::cout << "********************************" << clearatt << std::endl;

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

	iCurrTree = -1;
	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) {
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
		if (stuple.evt_RunNumber >246231) continue;  //drop the newly added data, 0k dataset in StupleV4

		//std::cout << "ele_num=" << stuple.ele_num << std::endl;
		//std::cout << "ele_num status=" << myChain->GetBranchStatus("ele_num") << std::endl;
		if (stuple.ele_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
		else if (stuple.ele_num > 0 && stuple.ele_num < 1) iEleEvts++;
		else if (stuple.ele_num < 1 && stuple.ele_num > 0) iPhoEvts++;
		else iOther++;

		//apply the met cut here to speed things up
		//met is not calculated again in this so this should be fine
		//2-9-2010
		if (stuple.met_Met < fMinMet || stuple.met_Met > fMaxMet) continue;

		CommonVars commVars;
		ElectronList centralEles;	//do not comment these out. I need these to recreate the jet list with all unused EM objects

		GetCommonVariables(stuple, commVars);
		CreateElectronLists(stuple, centralEles,0);

		DoMyStuff(commVars, centralEles);


	}

	//print out the job summary
	std::cout << magenta << "======== ZSel =======================" << std::endl;
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

	PrintJobSettings();

	const int itextWidth = 25;
	std::cout << "Pho Evts               = " << vMcCount[0] << std::endl;
	std::cout << "Pho+>=1Jet Evts        = " << vMcCount[1] << std::endl; 
	std::cout << "Event with > 2e        = " << iExtraEle << std::endl; 
	std::cout << cyan << "Data Written to        = " << GetHistFileName() << clearatt << std::endl;
	std::cout << "======= END ZSel ====================" << std::endl;

	EVT.clear();

	CleanUp();
	timer.Show("phoana_time");
} // Main


//-------------------------------------------------------------------
void ZSel::BookHistograms()
{

	/*Hsignal.z1j_zpt	= new TH1F("z1j_zpt","#gamma+>=1Jet :- p_{T}^{Z} ",200,-500,500);
	  Hsignal.z1j_zet	= new TH1F("z1j_zet","#gamma+>=1Jet :- E_{T}^{Z} ",250,0,500);
	  Hsignal.z1j_zmass	= new TH1F("z1j_zmass","#gamma+>=1Jet :- Z mass ",200,0,1000);
	  Hsignal.z1j_zeta	= new TH1F("z1j_zeta","#gamma+>=1Jet :- Event #Eta^{Z} ",50,-5,5);
	  Hsignal.z1j_zj_mass	= new TH1F("z1j_zj1mass","#gamma+>=1Jet :- Invariant Mass (Z, Lead Jet)",200,0,1000);

	  Hsignal.z2j_zpt	= new TH1F("z2j_zpt","#gamma+>=2Jets :- p_{T}^{Z} ",200,-500,500);
	  Hsignal.z2j_zet	= new TH1F("z2j_zet","#gamma+>=2Jets :- E_{T}^{Z} ",250,0,500);
	  Hsignal.z2j_zmass	= new TH1F("z2j_zmass","#gamma+>=2Jets :- Z mass ",200,0,1000);
	  Hsignal.z2j_zeta	= new TH1F("z2j_zeta","#gamma+>=2Jets :- Event #Eta^{Z} ",50,-5,5);
	  Hsignal.z2j_z1j_mass	= new TH1F("z2j_zj1mass","#gamma+>=2Jets :- Invariant Mass (Z, Lead Jet)",200,0,1000);
	  Hsignal.z2j_z2j_mass	= new TH1F("z2j_zj2mass","#gamma+>=2Jets :- Invariant Mass (Z, Second Lead Jet)",200,0,1000);
	  Hsignal.z2j_zj1j2_mass	= new TH1F("z2j_zj1j2mass","#gamma+>=2Jets :- Invariant Mass (Z, Two Lead Jets)",200,0,1000);
	  */

	Histograms HistoMan;
	std::string text2title1, text2title2;

	TDirectory *l2dir, *l3dir, *l4dir;

	//CENTRAL
	text2title1 = "Z :";
	text2title2 = "Z :";

	topDir->cd();

	l2dir = gDirectory->mkdir("Event","Event Histograms");
	l2dir->cd();
	HistoMan.GetEventHistograms(Hmc.Evt,text2title1.c_str());

	topDir->cd();
	l2dir = gDirectory->mkdir("Electron1","Electron-1 Histograms");
	l2dir->cd();
	HistoMan.GetPhotonHistograms(Hmc.Ele1,text2title1.c_str());


	topDir->cd();
	l2dir = gDirectory->mkdir("Electron2","Electron-2 Histograms");
	l2dir->cd();
	HistoMan.GetPhotonHistograms(Hmc.Ele2,text2title1.c_str());

	topDir->cd();
	l2dir = gDirectory->mkdir("2Jets","2Jets Histograms");
	l2dir->cd();
	HistoMan.GetTwoJetsHistograms(Hmc.TwoEles,text2title2.c_str());

} //BookHistograms

void ZSel::CreateElectronLists(const Stuple stuple, ElectronList& list, const int collection)
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

void ZSel::GetCommonVariables(const Stuple stuple, CommonVars& cVars)
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


void ZSel::DoMyStuff(CommonVars cVars, ElectronList cEles)
{	
	//_________________________________________________________________ require a one good vertex 		  
	if (cVars.vtx_NClass12 < 1) return;
	if (fabs(cVars.vtx_z) > 60.)  return;

	std::vector<int> UsedEle;

	//______________________________________________________________ select the photon
	for (unsigned int i = 0; i < cEles.ele_num; ++i )
	{
		if (cEles.ele_Etc[i]<20) continue;	//________ need only one. else problem
		if (fabs(cEles.ele_DetEta.at(i))>3.2) continue;
		if (cEles.ele_ConversionId[i] != 0) continue;
		if (cEles.ele_TightId[i] == 0) UsedEle.push_back(i);	//________ need only one. else problem
	}	

	int used = 0;				//____ to make sure evts are not fallen into multiple bins

	EvtTag_t evt;
	evt.Run 		= cVars.evt_RunNumber;
	evt.Evt 		= cVars.evt_EventNumber;

	if (UsedEle.size()>=2)
	{
		if (UsedEle.size()>=3) ++iExtraEle; 

		const int iEle1Index = UsedEle.at(0);
		const int iEle2Index = UsedEle.at(1);
		TLorentzVector tlEle1, tlEle2, tlSum;
		tlEle1.SetPxPyPzE(cEles.ele_Px[iEle1Index], cEles.ele_Py[iEle1Index], cEles.ele_Pz[iEle1Index], cEles.ele_E[iEle1Index]);
		tlEle2.SetPxPyPzE(cEles.ele_Px[iEle2Index], cEles.ele_Py[iEle2Index], cEles.ele_Pz[iEle2Index], cEles.ele_E[iEle2Index]);
		tlSum = tlEle1 + tlEle2; 
		if (tlSum.M() < 76. || tlSum.M() > 106.) return;
		used = JetSelAndHistFill(cVars, cEles, UsedEle, &Hmc, vMcCount);
	}

}

int ZSel::JetSelAndHistFill(const CommonVars Vars, 
		const ElectronList eles, 
		std::vector<int> vEleIndex, 
		Hist_t* hist,
		std::vector<unsigned int>& count, 
		const float fWgt,
		const bool debug
		)
{
	int iEle1Index = vEleIndex[0];	//assuming the electrons are sorted in Et
	int iEle2Index = vEleIndex[1];

	++count[0];
	float fWeight = 1;
	if (UseNvtxWeights())
	{
		if (Vars.vtx_NClass12 <= vNvtxWgts.size()) fWeight = vNvtxWgts.at(Vars.vtx_NClass12 - 1);
	}

	
	FillElectronHists(Vars, eles, &(hist->Ele1), iEle1Index, fWeight);
	FillElectronHists(Vars, eles, &(hist->Ele2), iEle2Index, fWeight);
	FillEventHists(Vars, eles, &(hist->Evt), fWeight);
	FillZHists(Vars, eles, &(hist->TwoEles), 0,1, fWeight);

	return 1;		//1=this event has been used!
}


void ZSel::FillEventHists(const CommonVars& vars,
								const ElectronList& eles, 
								Histograms::EventHists_t* hist,
								const float fWeight)
{ 
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}

	hist->Met->Fill(vars.met_Met, fWeight);
	hist->MetX->Fill(vars.met_MetX, fWeight);
	hist->MetY->Fill(vars.met_MetY, fWeight);
	//hist->Met_Gen_d->Fill(vars.met_Gen_d, fWeight);
	hist->Sumet->Fill(vars.met_SumEt, fWeight);
	hist->NVertices->Fill(vars.vtx_N, fWeight);
	hist->N12Vertices->Fill(vars.vtx_NClass12, fWeight);
	hist->VertexZ->Fill(vars.vtx_z, fWeight);
	hist->MetRun->Fill(vars.evt_RunNumber, vars.met_Met, fWeight);
	hist->SumetRun->Fill(vars.evt_RunNumber, vars.met_SumEt, fWeight);
	hist->NVerticesRun->Fill(vars.evt_RunNumber, vars.vtx_N, fWeight);
	hist->N12VerticesRun->Fill(vars.evt_RunNumber, vars.vtx_NClass12, fWeight);
	hist->N12VerticesVsMet->Fill(vars.met_Met, vars.vtx_NClass12, fWeight);
	hist->N12VerticesVsSumet->Fill(vars.met_SumEt, vars.vtx_NClass12, fWeight);
	hist->N12VerticesMetProf->Fill(vars.met_Met, vars.vtx_NClass12, fWeight);
	hist->N12VerticesSumetProf->Fill(vars.met_SumEt, vars.vtx_NClass12, fWeight);
	if (vars.tri_pho25iso == 1) hist->Triggers->Fill(0);
	if (vars.tri_pho50 == 1) hist->Triggers->Fill(1);
	if (vars.tri_pho70 == 1) hist->Triggers->Fill(2);

	hist->Ht->Fill(vars.met_Ht, fWeight);

}

void ZSel::FillElectronHists(const CommonVars& cVars, const ElectronList& eles, Histograms::PhotonHists_t* hist,
								const int iEleIndex, const float fWeight, const bool debug)
{
	hist->Detector->Fill(eles.ele_Detector[iEleIndex], fWeight);
	hist->DetEta->Fill(eles.ele_DetEta[iEleIndex], fWeight);
	hist->DetPhi->Fill(eles.ele_DetPhi[iEleIndex], fWeight);
	hist->EtCorr->Fill(eles.ele_Etc[iEleIndex], fWeight);
	hist->EventWeights->Fill(eles.ele_Etc[iEleIndex], fWeight);
	hist->XCes->Fill(eles.ele_XCes[iEleIndex], fWeight);
	hist->ZCes->Fill(eles.ele_ZCes[iEleIndex], fWeight);
	hist->HadEm->Fill(eles.ele_HadEm[iEleIndex], fWeight);
	hist->Chi2Mean->Fill(eles.ele_Chi2Mean[iEleIndex], fWeight);
	hist->N3d->Fill(eles.ele_N3d[iEleIndex], fWeight);
	hist->Iso4->Fill(eles.ele_Iso4[iEleIndex], fWeight);
	hist->EoverP->Fill(eles.ele_EoverP[iEleIndex], fWeight);
	hist->TrkPt->Fill(eles.ele_TrackPt[iEleIndex], fWeight);
	hist->TrkIso->Fill(eles.ele_TrkIso[iEleIndex], fWeight);
	hist->Ces2Wire->Fill(eles.ele_CesWireE2[iEleIndex], fWeight);
	hist->Ces2Strip->Fill(eles.ele_CesStripE2[iEleIndex], fWeight);
	hist->EmTime->Fill(eles.ele_EmTime[iEleIndex], fWeight);
	hist->PhiWedge->Fill(eles.ele_PhiWedge[iEleIndex], fWeight);
	hist->EmTimeVsRun->Fill(cVars.evt_RunNumber, eles.ele_EmTime[iEleIndex], fWeight);
	hist->EtCorrVsRun->Fill(cVars.evt_RunNumber, eles.ele_Etc[iEleIndex], fWeight);
	TLorentzVector tlEleVec(0,0,0,0);
	TVector2 tv2MetVec(cVars.met_MetX, cVars.met_MetY);
	tlEleVec.SetPxPyPzE(eles.ele_Px[iEleIndex], eles.ele_Py[iEleIndex], eles.ele_Pz[iEleIndex], eles.ele_E[iEleIndex]);

	//from V7 stuples and above has complete MET vector to do this
	const float dphi_pm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlEleVec.Phi())));
	hist->PhoMetDelPhi->Fill(dphi_pm);
	hist->ClosestPhoMetDelPhi_a4->Fill(dphi_pm);
	
	hist->PhoEtMetProf->Fill(cVars.met_Met, eles.ele_Etc[iEleIndex], fWeight);
	hist->PhoEMetProf->Fill(cVars.met_Met, tlEleVec.E(), fWeight);
	hist->PhoEtaMetProf->Fill(cVars.met_Met, tlEleVec.Eta(), fWeight);
}

void ZSel::PrintHeader(const CommonVars& cVars) const
{
	std::cout << " ============== " << cVars.evt_RunNumber << ", " << cVars.evt_EventNumber << std::endl;
}

void ZSel::PrintJobSettings()
{	

	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	std::cout << "MC Flag setting        = ";
	if (bMcFlag) std::cout << "MC" << std::endl;
	else std::cout << "DATA" << std::endl;
	std::cout << "min,max Vtx for Signal = " << std::setw(4) << GetMinClass12Vtx() << ", " << std::setw(4) << GetMaxClass12Vtx() << std::endl;
	std::cout << "Min/max MEt            = " << std::setw(4) << GetMinMet() << ", " << std::setw(4) << GetMaxMet()<< std::endl;
	std::cout << "Nvtx Weghting?         = " << std::setw(4) << UseNvtxWeights() << " (1=YES)"<< std::endl;
}
void ZSel::FillZHists(const CommonVars& cVars, const ElectronList& eles, 
						Histograms::TwoJetsHists_t* hist,
						const int iEle1Index, const int iEle2Index,
						const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	TLorentzVector tlEle1, tlEle2, tlSum;
	float DelEta = 0;

	tlEle1.SetPxPyPzE(eles.ele_Px[iEle1Index], eles.ele_Py[iEle1Index], eles.ele_Pz[iEle1Index], eles.ele_E[iEle1Index]);
	tlEle2.SetPxPyPzE(eles.ele_Px[iEle2Index], eles.ele_Py[iEle2Index], eles.ele_Pz[iEle2Index], eles.ele_E[iEle2Index]);
	DelEta = fabs(eles.ele_DetEta[iEle1Index] - eles.ele_DetEta[iEle2Index]);
	
	tlSum = tlEle1 + tlEle2; 
	float fEtRatio = (tlEle2.Energy() * TMath::Sin(tlEle2.Theta())) / (tlEle1.Energy() * TMath::Sin(tlEle1.Theta()));
	float DelPhi = fabs(tlEle1.DeltaPhi(tlEle2));
	float DelR = tlEle1.DeltaR(tlEle2);
	
	hist->InvMass->Fill(tlSum.M(), fWeight);
	hist->DelPhi->Fill(DelPhi, fWeight);
	hist->DelEta->Fill(DelEta, fWeight);
	hist->DelR->Fill(DelR, fWeight);
	hist->EtRatio->Fill(fEtRatio, fWeight);
}


