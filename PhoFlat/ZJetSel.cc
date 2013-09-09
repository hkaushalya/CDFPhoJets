/*{{{*/
/*  $Id: ZJetSel.cc,v 1.1 2011/05/31 18:00:20 samantha Exp $
 *  $Log: ZJetSel.cc,v $
 *  Revision 1.1  2011/05/31 18:00:20  samantha
 *  This is to get corrections to EWK MC by comparing Z+Jets data and W MC.
 *  I am simply treating Z as the photon and applying cuts similar to the photon.
 *  Also standard mass windows cuts on the di-electron mass to reduce the
 *  bacgrounds.
 *
 */

/*}}}*/

#include <iostream>
#include <sstream>
#include "PhoFlat/ZJetSel.hh"
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
ZJetSel::ZJetSel():
iProgressBy(50000),
sHistFileName("ZJetSel.root"),
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
	fMinZpt = 30.;
	fMaxZeta = 1.1;
	fMinJetEt = 15.;
	fMaxJetEta = 3.;

}

//-------------------------------------------------------------
void ZJetSel::Init(TChain *tree)
{
	iExtraEle = 0;

	readIn = new ReadInAll(tree, &stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
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
		//Z+jets without Z pt cut
/*		vNvtxWgts.push_back(0.758408);
		vNvtxWgts.push_back(1.04703);
		vNvtxWgts.push_back(1.27509);
		vNvtxWgts.push_back(1.62549);
		vNvtxWgts.push_back(2.86066);
		vNvtxWgts.push_back(6.33622);
		vNvtxWgts.push_back(14.9584);
		vNvtxWgts.push_back(29.0621);
*/

		vNvtxWgts.push_back(0.745858);
		vNvtxWgts.push_back(1.06972);
		vNvtxWgts.push_back(1.3202);
		vNvtxWgts.push_back(1.75381);
		vNvtxWgts.push_back(3.24774);
		vNvtxWgts.push_back(6.86493);
		vNvtxWgts.push_back(24.2292);
		vNvtxWgts.push_back(32.7094);


		
	}

} // Init

//-------------------------------------------------------------------
void ZJetSel::CleanUp()
{
	vMcCount.clear();

	histFile->Write();
	histFile->Close();
} 


//------ main -------------------------------------------------------
void ZJetSel::Main(TChain *ch, int iRunEvents)
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

		//std::cout << "ele_num=" << stuple.ele_num << std::endl;
		//std::cout << "ele_num status=" << myChain->GetBranchStatus("ele_num") << std::endl;
		if (stuple.ele_num > 0 && stuple.pho_num > 0) iPhoEleEvts++;
		else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
		else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;

		//apply the met cut here to speed things up
		//met is not calculated again in this so this should be fine
		//2-9-2010
		if (stuple.met_Met < fMinMet || stuple.met_Met > fMaxMet) continue;
		if (stuple.ele_num<1) continue;

		CommonVars cVars;
		ElectronList cEles;	//do not comment these out. I need these to recreate the jet list with all unused EM objects

		GetCommonVariables(stuple, cVars);
		CreateElectronLists(stuple, cEles,0);

		/*****************************************************************/
		//main routine that selections Z+JETS and fills hists
		
		//_________________________________________________________________ require a one good vertex 		  
		if (cVars.vtx_NClass12 < 1) continue;
		if (fabs(cVars.vtx_z) > 60.)  continue;

		std::vector<int> UsedEle, UsedPho;

		//______________________________________________________________ select the electrons
		for (unsigned int i = 0; i < cEles.ele_num; ++i )
		{
			if (cEles.ele_Etc[i]<20) continue;
			if (fabs(cEles.ele_DetEta.at(i))>3.2) continue;
			if (cEles.ele_ConversionId[i] != 0) continue;
			if (cEles.ele_TightId[i] == 0) UsedEle.push_back(i);
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
			TLorentzVector tlEle1, tlEle2, tlZVec;
			tlEle1.SetPxPyPzE(cEles.ele_Px[iEle1Index], cEles.ele_Py[iEle1Index], cEles.ele_Pz[iEle1Index], cEles.ele_E[iEle1Index]);
			tlEle2.SetPxPyPzE(cEles.ele_Px[iEle2Index], cEles.ele_Py[iEle2Index], cEles.ele_Pz[iEle2Index], cEles.ele_E[iEle2Index]);
			tlZVec = tlEle1 + tlEle2; 
			if (tlZVec.M() < 76. || tlZVec.M() > 106.) continue;  //mass cut to reduce background
			if (fabs(tlZVec.Eta())<GetMaxZEta()) continue;  //keep it central like the photon
			if (tlZVec.Pt()<GetMinZPt()) continue;


			++vMcCount[0]; //Z events;

			//to create new jet list need the list of photons and jets
			PhotonList cPhos;
			JetList cJets;
			CreatePhotonLists(stuple, cPhos,0);
			CreateJetLists(stuple, cJets,0);
			
			NewJetList nj;
			nj.AddUnused(cPhos, cEles, cJets, UsedPho, UsedEle, fMinJetEt, fMaxJetEta);

			//cut on njets
			if (cJets.jet_NJet15 < 1) continue;

			//MET clean up cut. require MET to be away from any jet above 15GeV.
			const float fDelPhiCut = 0.4;
			bool bad = false;
			float dphi_jm_closest = 10.0;
			float closest_jetPt = -99999.;
			const TLorentzVector tlMet(cVars.met_MetX,cVars.met_MetY,0,cVars.met_Met);
			for (int i=0; i <cJets.jet_num ; ++i)
			{
				TLorentzVector tlJetVec(cJets.jet_Px[i], cJets.jet_Py[i], cJets.jet_Pz[i], cJets.jet_E[i]);
				if (tlJetVec.Pt()< GetMinJetEt() || fabs(tlJetVec.Eta())> GetMaxJetEta()) continue;
				const float dphi_jm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tlMet.Phi())-TVector2::Phi_0_2pi(tlJetVec.Phi())));
				//std::cout << " ijet, dphi = " << i << ", " << dphi_jm << std::endl;
				if (dphi_jm <0.4 || dphi_jm>(TMath::Pi()-0.4))
				{
					//std::cout << red << "Failed MET CLEAN UP : jet Pt = " << tlJetVec.Pt() << clearatt << std::endl;
					bad = true;
				}
				if (dphi_jm < dphi_jm_closest) 
				{
					dphi_jm_closest = dphi_jm;
					closest_jetPt = tlJetVec.Pt();
				}
			}

			const float dphi_wm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tlMet.Phi())-TVector2::Phi_0_2pi(tlZVec.Phi())));
			hZMetDelPhi_b4->Fill(dphi_wm);
			hClosestJetMetDelPhi_b4->Fill(dphi_jm_closest);

			if (bad) continue;

			hZMetDelPhi_a4->Fill(dphi_wm);
			hClosestJetMetDelPhi_a4->Fill(dphi_jm_closest);
			hClosestJetPt_a4->Fill(closest_jetPt);

			++vMcCount[1];

			used = JetSelAndHistFill(cVars, cEles, cJets, UsedEle, &Hmc);
		}

		/*****************************************************************/

	}

	//print out the job summary
	std::cout << magenta << "======== ZJetSel =======================" << std::endl;
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
	std::cout << "Z Events               = " << vMcCount[0] << std::endl;
	std::cout << "Z+Jets Events          = " << vMcCount[1] << std::endl; 
	std::cout << "Event with > 2e        = " << iExtraEle << std::endl; 
	std::cout << cyan << "Data Written to        = " << GetHistFileName() << clearatt << std::endl;
	std::cout << "======= END ZJetSel ====================" << std::endl;

	EVT.clear();

	CleanUp();
	timer.Show("phoana_time");
} // Main


//-------------------------------------------------------------------
void ZJetSel::BookHistograms()
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

	TDirectory *l2dir;

	//CENTRAL
	text2title1 = "Z :";
	text2title2 = "Z :";

	topDir->cd();
		std::stringstream title1;
		title1 <<  text2title1 << " - Closest Jet Et>15GeV before cut: #delta#phi(jet,#slash{E}_{T})";
		hClosestJetMetDelPhi_a4	= new TH1F("ClosesetJetMetDelPhi_a4", title1.str().c_str(), 700, -3.5, 3.5);
		hClosestJetMetDelPhi_b4	= new TH1F("ClosesetJetMetDelPhi_b4", title1.str().c_str(), 700, -3.5, 3.5);
		hClosestJetPt_a4	= new TH1F("ClosesetJetPt_a4", title1.str().c_str(), 500, 0, 500);

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
	l2dir = gDirectory->mkdir("Z","Z Histograms");
	l2dir->cd();
	HistoMan.GetTwoJetsHistograms(Hmc.Z,text2title2.c_str());

		std::stringstream title2;
		title2 <<  text2title1 << "#delta#phi(Z,#slash{E}_{T})";
		hZMetDelPhi_b4	= new TH1F("ZMetDelPhi_b4",title2.str().c_str(), 700, -3.5, 3.5);
		hZMetDelPhi_a4	= new TH1F("ZMetDelPhi_a4",title2.str().c_str(), 700, -3.5, 3.5);

	topDir->cd();
	l2dir = gDirectory->mkdir("ZJet","ZJet Histograms");
	l2dir->cd();
	HistoMan.GetTwoJetsHistograms(Hmc.ZJet,text2title2.c_str());

	topDir->cd();
	l2dir = gDirectory->mkdir("Jet","Jet Histograms");
	l2dir->cd();
	HistoMan.GetJetHistograms(Hmc.Jet,text2title1.c_str());

	
} //BookHistograms
void ZJetSel::CreatePhotonLists(const Stuple stuple, PhotonList& pholist, const int collection)
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


void ZJetSel::CreateElectronLists(const Stuple stuple, ElectronList& list, const int collection)
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
void ZJetSel::CreateJetLists(const Stuple stuple, JetList& list, const int collection)
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
				if ( fabs(stuple.jet_Px[i]) > 5500 || 
					  fabs(stuple.jet_Py[i]) > 5500 || 
					  fabs(stuple.jet_Pz[i]) > 5500 || 
					  fabs(stuple.jet_E[i]) > 5500) //there is one event with E>2500?
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


void ZJetSel::GetCommonVariables(const Stuple stuple, CommonVars& cVars)
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


void ZJetSel::DoMyStuff(CommonVars cVars, ElectronList cEles, JetList cJets, PhotonList cPhos)
{	

}

int ZJetSel::JetSelAndHistFill(const CommonVars Vars, 
		const ElectronList eles, 
		const JetList jets, 
		std::vector<int> vEleIndex, 
		Hist_t* hist,
		const float fWgt,
		const bool debug
		)
{
	int iEle1Index = vEleIndex[0];	//assuming the electrons are sorted in Et
	int iEle2Index = vEleIndex[1];

	float fWeight = 1;
	if (UseNvtxWeights())
	{
		if (Vars.vtx_NClass12 <= vNvtxWgts.size()) fWeight = vNvtxWgts.at(Vars.vtx_NClass12 - 1);
	}

	
	FillElectronHists(Vars, eles, &(hist->Ele1), iEle1Index, fWeight);
	FillElectronHists(Vars, eles, &(hist->Ele2), iEle2Index, fWeight);
	FillEventHists(Vars, eles, &(hist->Evt), fWeight);
	FillZHists(Vars, eles, &(hist->Z), 0,1, fWeight);
	FillJetHists(Vars, jets, &(hist->Jet),0, fWeight);		//lead jet
	FillZJetHists(Vars, eles, jets, &(hist->ZJet), iEle1Index, iEle2Index, 0, fWeight);

	return 1;		//1=this event has been used!
}


void ZJetSel::FillEventHists(const CommonVars& vars,
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

void ZJetSel::FillElectronHists(const CommonVars& cVars, const ElectronList& eles, Histograms::PhotonHists_t* hist,
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

void ZJetSel::PrintHeader(const CommonVars& cVars) const
{
	std::cout << " ============== " << cVars.evt_RunNumber << ", " << cVars.evt_EventNumber << std::endl;
}

void ZJetSel::PrintJobSettings()
{	

	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	std::cout << "MC Flag setting        = ";
	if (bMcFlag) std::cout << "MC" << std::endl;
	else std::cout << "DATA" << std::endl;
	std::cout << "min,max Vtx for Signal = " << std::setw(4) << GetMinClass12Vtx() << ", " << std::setw(4) << GetMaxClass12Vtx() << std::endl;
	std::cout << "Min/max MEt            = " << std::setw(4) << GetMinMet() << ", " << std::setw(4) << GetMaxMet()<< std::endl;
	std::cout << "Min Z pt/ Max Eta      = " << std::setw(4) << GetMinZPt() << ", " << GetMaxZEta() << std::endl;
	std::cout << "Jet Min Et/ Max Eta    = " << std::setw(4) << GetMinJetEt() << ", " << std::setw(4) << GetMaxJetEta()<< std::endl;
	std::cout << "Nvtx Weghting?         = " << std::setw(4) << UseNvtxWeights() << " (1=YES)"<< std::endl;
}
void ZJetSel::FillZHists(const CommonVars& cVars, const ElectronList& eles, 
						Histograms::TwoJetsHists_t* hist,
						const int iEle1Index, const int iEle2Index,
						const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	TLorentzVector tlEle1, tlEle2, tlZVec;
	float DelEta = 0;

	tlEle1.SetPxPyPzE(eles.ele_Px[iEle1Index], eles.ele_Py[iEle1Index], eles.ele_Pz[iEle1Index], eles.ele_E[iEle1Index]);
	tlEle2.SetPxPyPzE(eles.ele_Px[iEle2Index], eles.ele_Py[iEle2Index], eles.ele_Pz[iEle2Index], eles.ele_E[iEle2Index]);
	DelEta = fabs(eles.ele_DetEta[iEle1Index] - eles.ele_DetEta[iEle2Index]);
	
	tlZVec = tlEle1 + tlEle2; 
	float fEtRatio = (tlEle2.Energy() * TMath::Sin(tlEle2.Theta())) / (tlEle1.Energy() * TMath::Sin(tlEle1.Theta()));
	float DelPhi = fabs(tlEle1.DeltaPhi(tlEle2));
	float DelR = tlEle1.DeltaR(tlEle2);
	
	hist->InvMass->Fill(tlZVec.M(), fWeight);
	hist->DelPhi->Fill(DelPhi, fWeight);
	hist->DelEta->Fill(DelEta, fWeight);
	hist->DelR->Fill(DelR, fWeight);
	hist->EtRatio->Fill(fEtRatio, fWeight);
	hist->Pt->Fill(tlZVec.Pt(), fWeight);
}

void ZJetSel::FillJetHists(const CommonVars& cVars, const JetList& jets, Histograms::JetHists_t* hist, 
								const int iJetIndex,	const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	hist->Emfr->Fill(jets.jet_Emfr[iJetIndex], fWeight);
	hist->DetEta->Fill(jets.jet_DetEta[iJetIndex], fWeight);
	hist->NTracks->Fill(jets.jet_Ntracks[iJetIndex], fWeight);
	if (jets.jet_Ntowers[iJetIndex]==0)
	{
		std::cout << __FUNCTION__ << "::" << __LINE__ << "A jet with ntower==0 found!!! ";
		//cVars.PrintHeader();
		//jets.DumpAll();
	}
	hist->NTowers->Fill(jets.jet_Ntowers[iJetIndex], fWeight);
	hist->EtCorr->Fill(jets.jet_Pt[iJetIndex], fWeight);

	TLorentzVector tlVec(jets.jet_Px[iJetIndex], jets.jet_Py[iJetIndex], jets.jet_Pz[iJetIndex], jets.jet_E[iJetIndex]);
	hist->EvtPhi->Fill(tlVec.Phi(), fWeight);
	hist->EvtEta->Fill(tlVec.Eta(), fWeight);
	//hist->SeedIEta->Fill(jets.jet_SeedIEta[iJetIndex], fWeight);
	//hist->SeedIPhi->Fill(jets.jet_SeedIPhi[iJetIndex], fWeight);
	
	TVector2 tv2MetVec(cVars.met_MetX, cVars.met_MetY);

	//from V7 stuples and above has complete MET vector to do this
	const float dphi_jm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlVec.Phi())));
	hist->JetMetDelPhi->Fill(dphi_jm);
	hist->ClosestJetMetDelPhi_a4->Fill(dphi_jm);
	
	hist->JetEtMetProf->Fill(cVars.met_Met, tlVec.Pt(), fWeight);
	hist->JetEMetProf->Fill(cVars.met_Met, tlVec.E(), fWeight);
	hist->JetEtaMetProf->Fill(cVars.met_Met, tlVec.Eta(), fWeight);
}

void ZJetSel::FillZJetHists(const CommonVars& cVars, const ElectronList& eles, const JetList& jets, 
					Histograms::TwoJetsHists_t* hist, const int iEle1Index, const int iEle2Index,
					const int iJetIndex, const float fWeight)
{
	//std::cout << "\tfillpho1jethist:: npho=" << phos.pho_num << std::endl;
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	if (iJetIndex < 0 || iJetIndex > (int)jets.jet_NJet15)
	{
		StdOut(__FILE__,__LINE__,3,"Invalid Jet index!");
		assert(false);
	}
	float InvMass 	= 0;
	float DelPhi 	= 0;
	float DelEta 	= 0;
	float DelR 		= 0;
	float EtRatio 	= 0;
	float JetEmFr 	= 0;

	TLorentzVector tlEle1, tlEle2, tlZVec, tlJet, tlZJetVec;

	tlEle1.SetPxPyPzE(eles.ele_Px[iEle1Index], eles.ele_Py[iEle1Index], eles.ele_Pz[iEle1Index], eles.ele_E[iEle1Index]);
	tlEle2.SetPxPyPzE(eles.ele_Px[iEle2Index], eles.ele_Py[iEle2Index], eles.ele_Pz[iEle2Index], eles.ele_E[iEle2Index]);
	tlZVec = tlEle1 + tlEle2; 

	tlJet.SetPxPyPzE(jets.jet_Px[iJetIndex], jets.jet_Py[iJetIndex], jets.jet_Pz[iJetIndex], jets.jet_E[iJetIndex]);

	tlZJetVec = tlJet + tlZVec; 
	InvMass 	= tlZJetVec.M();
	if (InvMass< 1) 
	{
		std::cout << __FILE__ <<":" << __FUNCTION__ << ":" << __LINE__ <<	"::Event with InvMass = " << InvMass << " <5" << std::endl;
		std::cout << "Run,Event Number::" << cVars.evt_RunNumber << ", " << cVars.evt_EventNumber << std::endl;
		std::cout << "Jet Pt/E = " << tlJet.Pt() << ", " << tlJet.E() << std::endl;
		std::cout << "Pho Pt/E = " << tlZVec.Pt() << ", " << tlZVec.E() << std::endl;
		std::cout << "Sum Pt/E = " << tlZJetVec.Pt() << ", " << tlZJetVec.E() << std::endl;
		PrintHeader(cVars);
		jets.Dump();
		stuple.Dump(0);
	}
	
	DelPhi 	= fabs(tlZVec.DeltaPhi(tlJet));
	DelEta 	= fabs(tlZVec.Eta() - jets.jet_DetEta[iJetIndex]);
	DelR 		= tlZVec.DeltaR(tlJet);
	EtRatio 	= tlJet.Pt()/ tlZVec.Pt();
	//JetEmFr 	= jets.jet_Emfr[iJetIndex];
	
	hist->InvMass->Fill(InvMass, fWeight);
	hist->DelPhi->Fill(DelPhi, fWeight);
	hist->DelEta->Fill(DelEta, fWeight);
	hist->DelR->Fill(DelR, fWeight);
	hist->EtRatio->Fill(EtRatio, fWeight);

}


