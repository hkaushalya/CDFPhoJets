#include <iostream>
#include "PhoFlat/ReadInAll.hh"
/*{{{*/
/*******************************************************************
 *
 * $Id: ReadInAll.cc,v 1.9 2011/05/02 17:29:16 samantha Exp $
 * $Log: ReadInAll.cc,v $
 * Revision 1.9  2011/05/02 17:29:16  samantha
 * ADDED: SecVtx b-tag of jets.
 *
 * Revision 1.8  2011/04/26 23:49:26  samantha
 * ADDED: Raw jet info for read-in.
 *
 * Revision 1.7  2010/04/13 16:16:09  samantha
 * ADDED: MetX, MetY and RawMet branches.
 *
 * Revision 1.6  2010/02/08 20:56:42  samantha
 * UNCOMMENTED: The cpr weights. I need those to calculate the fake fraction
 * in different samples. (like with high phot et)
 * COMMETED: The generated MET stuff met_gen_xxx as MET model is out of the way
 * now.
 *
 * Revision 1.5  2010/01/12 20:33:23  samantha
 * COMMENTED: CPR weights are not read as they are not used any more.
 *
 * $----------------------------$
 *******************************************************************
 * revision 1.4
 * date: 2009/12/29 04:15:50;  author: samantha;  state: Exp;  lines: +19 -5
 * ADDED: Reading-in PhoBlockIndex and EleBlockIndex.
 * MODIFIED: Commented out (or not reading) Generated MEt,'met_Gen_xxx'.
 * ----------------------------
 * revision 1.3
 * date: 2008/07/16 21:10:53;  author: samantha;  state: Exp;  lines: +59 -1
 * added new CES/CPR weight branches for read in.
 * ----------------------------
 * revision 1.2
 * date: 2008/06/28 01:05:11;  author: samantha;  state: Exp;  lines: +656 -212
 * updated with new central/em up/em down photon and electron branches 
 * and central/jes up/jes down jet branches. can enable/disable any 
 * collection (photon/elecron/jet).
 * ----------------------------
 * revision 1.1
 * date: 2008/03/18 23:32:23;  author: samantha;  state: Exp;
 * this handles my data read back from the Stuples. Need to init 
 * with the TChain of Stuples and call for each iteration to 
 * fill up the branches
 *******************************************************************/
/*}}}*/

//-------------------------------------------------------------------
void ReadInAll::GetEntry(int entry)
{
	myStuple->Init();   			// must reset before storing next event. 
	myTree->GetEntry(entry);
//	std::cout << "pho_num=" << myStuple->pho_num << std::endl;

}

//-------------------------------------------------------------------
ReadInAll::ReadInAll(TChain* tree, Stuple* st)
{
	if (!tree) {
		std::cout << __FILE__ << "::" << __LINE__ <<"::tree not found" << std::endl;
		return;
	}
	if (!st) {
		std::cout << __FILE__ << "::" << __LINE__ <<"::Stuple not found" << std::endl;
		return;
	}
	myTree = tree;
	myStuple = st;

	myTree->SetBranchStatus("*", 0);

	bAllPho     = false; 
	bCentralPho = false;
	bEmUpPho    = false;
	bEmDownPho  = false;
	bAllEle     = false;
	bCentralEle = false;
	bEmUpEle    = false;
	bEmDownEle  = false;
	bAllJet     = false;
	bCentralJet = false;
	bJesUpJet   = false;
	bJesDownJet = false;
	bMet = false;
	bGen = false;
}

//-------------------------------------------------------------------
void ReadInAll::Init()
{
	myTree->SetBranchStatus("evt_RunNumber", 1);
	myTree->SetBranchStatus("evt_EventNumber", 1);
	myTree->SetBranchStatus("evt_McFlag", 1);
	myTree->SetBranchStatus("tri_pho25iso", 1);
	myTree->SetBranchStatus("tri_pho50", 1);
	myTree->SetBranchStatus("tri_pho70", 1);

	myTree->SetBranchAddress("evt_RunNumber", &myStuple->evt_RunNumber);
	myTree->SetBranchAddress("evt_EventNumber", &myStuple->evt_EventNumber);
	myTree->SetBranchAddress("evt_McFlag", &myStuple->evt_McFlag);
	myTree->SetBranchAddress("tri_pho25iso", &myStuple->tri_pho25iso);
	myTree->SetBranchAddress("tri_pho50", &myStuple->tri_pho50);
	myTree->SetBranchAddress("tri_pho70", &myStuple->tri_pho70);
	
	myTree->SetBranchStatus("vtx_N", 1);
	myTree->SetBranchStatus("vtx_NClass12", 1);
	myTree->SetBranchStatus("vtx_z", 1);				// highest sum pt vertex Z
	myTree->SetBranchStatus("vtx_Ntracks", 1);		// in highest sum pt vertex
	myTree->SetBranchStatus("vtx_SumPt", 1);			// highest

	myTree->SetBranchAddress("vtx_N", &myStuple->vtx_N);
	myTree->SetBranchAddress("vtx_NClass12", &myStuple->vtx_NClass12);
	myTree->SetBranchAddress("vtx_z", &myStuple->vtx_z);				// highest sum pt vertex Z
	myTree->SetBranchAddress("vtx_Ntracks", &myStuple->vtx_Ntracks);		// in highest sum pt vertex
	myTree->SetBranchAddress("vtx_SumPt", &myStuple->vtx_SumPt);			// highest

	
	if (bAllPho || bCentralPho) {
		std::cout << "Enabling Central Photons Read-in." << std::endl;
		myTree->SetBranchStatus("pho_num", 1);
		myTree->SetBranchStatus("pho_Ntight", 1);
		myTree->SetBranchStatus("pho_Nloose", 1);
		myTree->SetBranchStatus("pho_Index", 1);
		myTree->SetBranchStatus("pho_PhoBlockIndex", 1);
		myTree->SetBranchStatus("pho_Etc", 1);
		myTree->SetBranchStatus("pho_E", 1);
		myTree->SetBranchStatus("pho_Px", 1);
		myTree->SetBranchStatus("pho_Py", 1);
		myTree->SetBranchStatus("pho_Pz", 1);
		myTree->SetBranchStatus("pho_Detector", 1);
		myTree->SetBranchStatus("pho_DetEta", 1);
		myTree->SetBranchStatus("pho_DetPhi", 1);
		myTree->SetBranchStatus("pho_XCes", 1);
		myTree->SetBranchStatus("pho_ZCes", 1);
		myTree->SetBranchStatus("pho_HadEm", 1);
		myTree->SetBranchStatus("pho_Chi2Mean", 1);
		myTree->SetBranchStatus("pho_N3d", 1);
		myTree->SetBranchStatus("pho_Iso4", 1);
		myTree->SetBranchStatus("pho_TrkPt", 1);		//from photon block
		myTree->SetBranchStatus("pho_TrkIso", 1);		//from photon block
		myTree->SetBranchStatus("pho_CesWireE2", 1);
		myTree->SetBranchStatus("pho_CesStripE2", 1);
		myTree->SetBranchStatus("pho_PhiWedge", 1);
		myTree->SetBranchStatus("pho_NMuonStubs", 1);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchStatus("pho_EmTime", 1);
		myTree->SetBranchStatus("pho_TightId", 1);
		myTree->SetBranchStatus("pho_LooseId", 1);
		myTree->SetBranchStatus("pho_PhoenixId", 1);
		myTree->SetBranchStatus("pho_Halo_seedWedge", 1);
		myTree->SetBranchStatus("pho_Halo_eastNhad", 1);
		myTree->SetBranchStatus("pho_Halo_westNhad", 1);
		myTree->SetBranchStatus("pho_matchJetIndex", 1); 	//if a matching jet is found and removed from jet list

		//do not disable COR wgts. I fill the weights
		//with each job so I'll have the true photon fraction 
		//for normalization. 01-14-2010
		myTree->SetBranchStatus("pho_CprWgt", 1);
		myTree->SetBranchStatus("pho_CprSys1", 1);
		myTree->SetBranchStatus("pho_CprSys2", 1);
		myTree->SetBranchStatus("pho_CprSys3", 1);
		myTree->SetBranchStatus("pho_CprSys4", 1);
		myTree->SetBranchStatus("pho_CprSys5", 1);
		myTree->SetBranchStatus("pho_CprSys6", 1);
		myTree->SetBranchStatus("pho_CprSys7", 1);
		myTree->SetBranchStatus("pho_CprSys8", 1);
	
		
		myTree->SetBranchAddress("pho_num", &myStuple->pho_num);
		myTree->SetBranchAddress("pho_Ntight", &myStuple->pho_Ntight);			// # of tight photons
		myTree->SetBranchAddress("pho_Nloose", &myStuple->pho_Nloose);			// # of loose photons
		myTree->SetBranchAddress("pho_Index", &myStuple->pho_Index);
		myTree->SetBranchAddress("pho_PhoBlockIndex", &myStuple->pho_PhoBlockIndex);
		myTree->SetBranchAddress("pho_Etc", &myStuple->pho_Etc);
		myTree->SetBranchAddress("pho_E", &myStuple->pho_E);
		myTree->SetBranchAddress("pho_Px", &myStuple->pho_Px);
		myTree->SetBranchAddress("pho_Py", &myStuple->pho_Py);
		myTree->SetBranchAddress("pho_Pz", &myStuple->pho_Pz);
		myTree->SetBranchAddress("pho_Detector", &myStuple->pho_Detector);
		myTree->SetBranchAddress("pho_DetEta", &myStuple->pho_DetEta);
		myTree->SetBranchAddress("pho_DetPhi", &myStuple->pho_DetPhi);
		myTree->SetBranchAddress("pho_XCes", &myStuple->pho_XCes);
		myTree->SetBranchAddress("pho_ZCes", &myStuple->pho_ZCes);
		myTree->SetBranchAddress("pho_HadEm", &myStuple->pho_HadEm);
		myTree->SetBranchAddress("pho_Chi2Mean", &myStuple->pho_Chi2Mean);
		myTree->SetBranchAddress("pho_N3d", &myStuple->pho_N3d);
		myTree->SetBranchAddress("pho_Iso4", &myStuple->pho_Iso4);
		myTree->SetBranchAddress("pho_TrkPt", &myStuple->pho_TrkPt);		//from photon block
		myTree->SetBranchAddress("pho_TrkIso", &myStuple->pho_TrkIso);		//from photon block
		myTree->SetBranchAddress("pho_CesWireE2", &myStuple->pho_CesWireE2);
		myTree->SetBranchAddress("pho_CesStripE2", &myStuple->pho_CesStripE2);
		myTree->SetBranchAddress("pho_PhiWedge", &myStuple->pho_PhiWedge);
		myTree->SetBranchAddress("pho_NMuonStubs", &myStuple->pho_NMuonStubs);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchAddress("pho_EmTime", &myStuple->pho_EmTime);
		myTree->SetBranchAddress("pho_TightId", &myStuple->pho_TightId);
		myTree->SetBranchAddress("pho_LooseId", &myStuple->pho_LooseId);
		myTree->SetBranchAddress("pho_PhoenixId", &myStuple->pho_PhoenixId);
		myTree->SetBranchAddress("pho_Halo_seedWedge", &myStuple->pho_Halo_seedWedge);
		myTree->SetBranchAddress("pho_Halo_eastNhad", &myStuple->pho_Halo_eastNhad);
		myTree->SetBranchAddress("pho_Halo_westNhad", &myStuple->pho_Halo_westNhad);
		myTree->SetBranchAddress("pho_matchJetIndex", &myStuple->pho_matchJetIndex); 	//if a matching jet is found and removed from jet list

		myTree->SetBranchAddress("pho_CprWgt", &myStuple->pho_CprWgt);
		myTree->SetBranchAddress("pho_CprSys1", &myStuple->pho_CprSys1);
		myTree->SetBranchAddress("pho_CprSys2", &myStuple->pho_CprSys2);
		myTree->SetBranchAddress("pho_CprSys3", &myStuple->pho_CprSys3);
		myTree->SetBranchAddress("pho_CprSys4", &myStuple->pho_CprSys4);
		myTree->SetBranchAddress("pho_CprSys5", &myStuple->pho_CprSys5);
		myTree->SetBranchAddress("pho_CprSys6", &myStuple->pho_CprSys6);
		myTree->SetBranchAddress("pho_CprSys7", &myStuple->pho_CprSys7);
		myTree->SetBranchAddress("pho_CprSys8", &myStuple->pho_CprSys8);
		
	}

	if (bAllPho || bEmUpPho) {
			// Pho EM E uncertainty up 1%
		std::cout << "Enabling EM UP Photons Read-in." << std::endl;
		myTree->SetBranchStatus("pho_up_num", 1);
		myTree->SetBranchStatus("pho_up_Ntight", 1);
		myTree->SetBranchStatus("pho_up_Nloose", 1);
		myTree->SetBranchStatus("pho_up_Index", 1);
		myTree->SetBranchStatus("pho_up_PhoBlockIndex", 1);
		myTree->SetBranchStatus("pho_up_Etc", 1);
		myTree->SetBranchStatus("pho_up_E", 1);
		myTree->SetBranchStatus("pho_up_Px", 1);
		myTree->SetBranchStatus("pho_up_Py", 1);
		myTree->SetBranchStatus("pho_up_Pz", 1);
		myTree->SetBranchStatus("pho_up_Detector", 1);
		myTree->SetBranchStatus("pho_up_DetEta", 1);
		myTree->SetBranchStatus("pho_up_DetPhi", 1);
		myTree->SetBranchStatus("pho_up_XCes", 1);
		myTree->SetBranchStatus("pho_up_ZCes", 1);
		myTree->SetBranchStatus("pho_up_HadEm", 1);
		myTree->SetBranchStatus("pho_up_Chi2Mean", 1);
		myTree->SetBranchStatus("pho_up_N3d", 1);
		myTree->SetBranchStatus("pho_up_Iso4", 1);
		myTree->SetBranchStatus("pho_up_TrkPt", 1);		//from photon block
		myTree->SetBranchStatus("pho_up_TrkIso", 1);		//from photon block
		myTree->SetBranchStatus("pho_up_CesWireE2", 1);
		myTree->SetBranchStatus("pho_up_CesStripE2", 1);
		myTree->SetBranchStatus("pho_up_PhiWedge", 1);
		myTree->SetBranchStatus("pho_up_NMuonStubs", 1);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchStatus("pho_up_EmTime", 1);
		myTree->SetBranchStatus("pho_up_TightId", 1);
		myTree->SetBranchStatus("pho_up_LooseId", 1);
		myTree->SetBranchStatus("pho_up_PhoenixId", 1);
		myTree->SetBranchStatus("pho_up_Halo_seedWedge", 1);
		myTree->SetBranchStatus("pho_up_Halo_eastNhad", 1);
		myTree->SetBranchStatus("pho_up_Halo_westNhad", 1);
		myTree->SetBranchStatus("pho_up_matchJetIndex", 1); 	//if a matching jet is found and removed from jet list

		myTree->SetBranchStatus("pho_up_CprWgt", 1);
		myTree->SetBranchStatus("pho_up_CprSys1", 1);
		myTree->SetBranchStatus("pho_up_CprSys2", 1);
		myTree->SetBranchStatus("pho_up_CprSys3", 1);
		myTree->SetBranchStatus("pho_up_CprSys4", 1);
		myTree->SetBranchStatus("pho_up_CprSys5", 1);
		myTree->SetBranchStatus("pho_up_CprSys6", 1);
		myTree->SetBranchStatus("pho_up_CprSys7", 1);
		myTree->SetBranchStatus("pho_up_CprSys8", 1);


		myTree->SetBranchAddress("pho_up_num", &myStuple->pho_up_num);
		myTree->SetBranchAddress("pho_up_Ntight", &myStuple->pho_up_Ntight);			// # of tight photons
		myTree->SetBranchAddress("pho_up_Nloose", &myStuple->pho_up_Nloose);			// # of loose photons
		myTree->SetBranchAddress("pho_up_Index", &myStuple->pho_up_Index);
		myTree->SetBranchAddress("pho_up_PhoBlockIndex", &myStuple->pho_up_PhoBlockIndex);
		myTree->SetBranchAddress("pho_up_Etc", &myStuple->pho_up_Etc);
		myTree->SetBranchAddress("pho_up_E", &myStuple->pho_up_E);
		myTree->SetBranchAddress("pho_up_Px", &myStuple->pho_up_Px);
		myTree->SetBranchAddress("pho_up_Py", &myStuple->pho_up_Py);
		myTree->SetBranchAddress("pho_up_Pz", &myStuple->pho_up_Pz);
		myTree->SetBranchAddress("pho_up_Detector", &myStuple->pho_up_Detector);
		myTree->SetBranchAddress("pho_up_DetEta", &myStuple->pho_up_DetEta);
		myTree->SetBranchAddress("pho_up_DetPhi", &myStuple->pho_up_DetPhi);
		myTree->SetBranchAddress("pho_up_XCes", &myStuple->pho_up_XCes);
		myTree->SetBranchAddress("pho_up_ZCes", &myStuple->pho_up_ZCes);
		myTree->SetBranchAddress("pho_up_HadEm", &myStuple->pho_up_HadEm);
		myTree->SetBranchAddress("pho_up_Chi2Mean", &myStuple->pho_up_Chi2Mean);
		myTree->SetBranchAddress("pho_up_N3d", &myStuple->pho_up_N3d);
		myTree->SetBranchAddress("pho_up_Iso4", &myStuple->pho_up_Iso4);
		myTree->SetBranchAddress("pho_up_TrkPt", &myStuple->pho_up_TrkPt);		//from photon block
		myTree->SetBranchAddress("pho_up_TrkIso", &myStuple->pho_up_TrkIso);		//from photon block
		myTree->SetBranchAddress("pho_up_CesWireE2", &myStuple->pho_up_CesWireE2);
		myTree->SetBranchAddress("pho_up_CesStripE2", &myStuple->pho_up_CesStripE2);
		myTree->SetBranchAddress("pho_up_PhiWedge", &myStuple->pho_up_PhiWedge);
		myTree->SetBranchAddress("pho_up_NMuonStubs", &myStuple->pho_up_NMuonStubs);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchAddress("pho_up_EmTime", &myStuple->pho_up_EmTime);
		myTree->SetBranchAddress("pho_up_TightId", &myStuple->pho_up_TightId);
		myTree->SetBranchAddress("pho_up_LooseId", &myStuple->pho_up_LooseId);
		myTree->SetBranchAddress("pho_up_PhoenixId", &myStuple->pho_up_PhoenixId);
		myTree->SetBranchAddress("pho_up_Halo_seedWedge", &myStuple->pho_up_Halo_seedWedge);
		myTree->SetBranchAddress("pho_up_Halo_eastNhad", &myStuple->pho_up_Halo_eastNhad);
		myTree->SetBranchAddress("pho_up_Halo_westNhad", &myStuple->pho_up_Halo_westNhad);
		myTree->SetBranchAddress("pho_up_matchJetIndex", &myStuple->pho_up_matchJetIndex); 	//if a matching jet is found and removed from jet list

		myTree->SetBranchAddress("pho_up_CprWgt", &myStuple->pho_up_CprWgt);
		myTree->SetBranchAddress("pho_up_CprSys1", &myStuple->pho_up_CprSys1);
		myTree->SetBranchAddress("pho_up_CprSys2", &myStuple->pho_up_CprSys2);
		myTree->SetBranchAddress("pho_up_CprSys3", &myStuple->pho_up_CprSys3);
		myTree->SetBranchAddress("pho_up_CprSys4", &myStuple->pho_up_CprSys4);
		myTree->SetBranchAddress("pho_up_CprSys5", &myStuple->pho_up_CprSys5);
		myTree->SetBranchAddress("pho_up_CprSys6", &myStuple->pho_up_CprSys6);
		myTree->SetBranchAddress("pho_up_CprSys7", &myStuple->pho_up_CprSys7);
		myTree->SetBranchAddress("pho_up_CprSys8", &myStuple->pho_up_CprSys8);
		
	}

	if (bAllPho || bEmDownPho) {
			// Pho EM E uncertainty down 1%
		std::cout << "Enabling EM DOWN Photons Read-in." << std::endl;
		myTree->SetBranchStatus("pho_down_num", 1);
		myTree->SetBranchStatus("pho_down_Ntight", 1);
		myTree->SetBranchStatus("pho_down_Nloose", 1);
		myTree->SetBranchStatus("pho_down_Index", 1);
		myTree->SetBranchStatus("pho_down_PhoBlockIndex", 1);
		myTree->SetBranchStatus("pho_down_Etc", 1);
		myTree->SetBranchStatus("pho_down_E", 1);
		myTree->SetBranchStatus("pho_down_Px", 1);
		myTree->SetBranchStatus("pho_down_Py", 1);
		myTree->SetBranchStatus("pho_down_Pz", 1);
		myTree->SetBranchStatus("pho_down_Detector", 1);
		myTree->SetBranchStatus("pho_down_DetEta", 1);
		myTree->SetBranchStatus("pho_down_DetPhi", 1);
		myTree->SetBranchStatus("pho_down_XCes", 1);
		myTree->SetBranchStatus("pho_down_ZCes", 1);
		myTree->SetBranchStatus("pho_down_HadEm", 1);
		myTree->SetBranchStatus("pho_down_Chi2Mean", 1);
		myTree->SetBranchStatus("pho_down_N3d", 1);
		myTree->SetBranchStatus("pho_down_Iso4", 1);
		myTree->SetBranchStatus("pho_down_TrkPt", 1);		//from photon block
		myTree->SetBranchStatus("pho_down_TrkIso", 1);		//from photon block
		myTree->SetBranchStatus("pho_down_CesWireE2", 1);
		myTree->SetBranchStatus("pho_down_CesStripE2", 1);
		myTree->SetBranchStatus("pho_down_PhiWedge", 1);
		myTree->SetBranchStatus("pho_down_NMuonStubs", 1);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchStatus("pho_down_EmTime", 1);
		myTree->SetBranchStatus("pho_down_TightId", 1);
		myTree->SetBranchStatus("pho_down_LooseId", 1);
		myTree->SetBranchStatus("pho_down_PhoenixId", 1);
		myTree->SetBranchStatus("pho_down_Halo_seedWedge", 1);
		myTree->SetBranchStatus("pho_down_Halo_eastNhad", 1);
		myTree->SetBranchStatus("pho_down_Halo_westNhad", 1);
		myTree->SetBranchStatus("pho_down_matchJetIndex", 1); 	//if a matching jet is found and removed from jet list

		myTree->SetBranchStatus("pho_down_CprWgt", 1);
		myTree->SetBranchStatus("pho_down_CprSys1", 1);
		myTree->SetBranchStatus("pho_down_CprSys2", 1);
		myTree->SetBranchStatus("pho_down_CprSys3", 1);
		myTree->SetBranchStatus("pho_down_CprSys4", 1);
		myTree->SetBranchStatus("pho_down_CprSys5", 1);
		myTree->SetBranchStatus("pho_down_CprSys6", 1);
		myTree->SetBranchStatus("pho_down_CprSys7", 1);
		myTree->SetBranchStatus("pho_down_CprSys8", 1);

		myTree->SetBranchAddress("pho_down_num", &myStuple->pho_down_num);
		myTree->SetBranchAddress("pho_down_Ntight", &myStuple->pho_down_Ntight);			// # of tight photons
		myTree->SetBranchAddress("pho_down_Nloose", &myStuple->pho_down_Nloose);			// # of loose photons
		myTree->SetBranchAddress("pho_down_Index", &myStuple->pho_down_Index);
		myTree->SetBranchAddress("pho_down_PhoBlockIndex", &myStuple->pho_down_PhoBlockIndex);
		myTree->SetBranchAddress("pho_down_Etc", &myStuple->pho_down_Etc);
		myTree->SetBranchAddress("pho_down_E", &myStuple->pho_down_E);
		myTree->SetBranchAddress("pho_down_Px", &myStuple->pho_down_Px);
		myTree->SetBranchAddress("pho_down_Py", &myStuple->pho_down_Py);
		myTree->SetBranchAddress("pho_down_Pz", &myStuple->pho_down_Pz);
		myTree->SetBranchAddress("pho_down_Detector", &myStuple->pho_down_Detector);
		myTree->SetBranchAddress("pho_down_DetEta", &myStuple->pho_down_DetEta);
		myTree->SetBranchAddress("pho_down_DetPhi", &myStuple->pho_down_DetPhi);
		myTree->SetBranchAddress("pho_down_XCes", &myStuple->pho_down_XCes);
		myTree->SetBranchAddress("pho_down_ZCes", &myStuple->pho_down_ZCes);
		myTree->SetBranchAddress("pho_down_HadEm", &myStuple->pho_down_HadEm);
		myTree->SetBranchAddress("pho_down_Chi2Mean", &myStuple->pho_down_Chi2Mean);
		myTree->SetBranchAddress("pho_down_N3d", &myStuple->pho_down_N3d);
		myTree->SetBranchAddress("pho_down_Iso4", &myStuple->pho_down_Iso4);
		myTree->SetBranchAddress("pho_down_TrkPt", &myStuple->pho_down_TrkPt);		//from photon block
		myTree->SetBranchAddress("pho_down_TrkIso", &myStuple->pho_down_TrkIso);		//from photon block
		myTree->SetBranchAddress("pho_down_CesWireE2", &myStuple->pho_down_CesWireE2);
		myTree->SetBranchAddress("pho_down_CesStripE2", &myStuple->pho_down_CesStripE2);
		myTree->SetBranchAddress("pho_down_PhiWedge", &myStuple->pho_down_PhiWedge);
		myTree->SetBranchAddress("pho_down_NMuonStubs", &myStuple->pho_down_NMuonStubs);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchAddress("pho_down_EmTime", &myStuple->pho_down_EmTime);
		myTree->SetBranchAddress("pho_down_TightId", &myStuple->pho_down_TightId);
		myTree->SetBranchAddress("pho_down_LooseId", &myStuple->pho_down_LooseId);
		myTree->SetBranchAddress("pho_down_PhoenixId", &myStuple->pho_down_PhoenixId);
		myTree->SetBranchAddress("pho_down_Halo_seedWedge", &myStuple->pho_down_Halo_seedWedge);
		myTree->SetBranchAddress("pho_down_Halo_eastNhad", &myStuple->pho_down_Halo_eastNhad);
		myTree->SetBranchAddress("pho_down_Halo_westNhad", &myStuple->pho_down_Halo_westNhad);
		myTree->SetBranchAddress("pho_down_matchJetIndex", &myStuple->pho_down_matchJetIndex); 	//if a matching jet is found and removed from jet list
		
		myTree->SetBranchAddress("pho_down_CprWgt", &myStuple->pho_down_CprWgt);
		myTree->SetBranchAddress("pho_down_CprSys1", &myStuple->pho_down_CprSys1);
		myTree->SetBranchAddress("pho_down_CprSys2", &myStuple->pho_down_CprSys2);
		myTree->SetBranchAddress("pho_down_CprSys3", &myStuple->pho_down_CprSys3);
		myTree->SetBranchAddress("pho_down_CprSys4", &myStuple->pho_down_CprSys4);
		myTree->SetBranchAddress("pho_down_CprSys5", &myStuple->pho_down_CprSys5);
		myTree->SetBranchAddress("pho_down_CprSys6", &myStuple->pho_down_CprSys6);
		myTree->SetBranchAddress("pho_down_CprSys7", &myStuple->pho_down_CprSys7);
		myTree->SetBranchAddress("pho_down_CprSys8", &myStuple->pho_down_CprSys8);
		
	}

	if (bAllEle || bCentralEle) {
		std::cout << "Enabling Central Electrons Read-in." << std::endl;
		myTree->SetBranchStatus("ele_num", 1);
		myTree->SetBranchStatus("ele_Index", 1);
		myTree->SetBranchStatus("ele_EleBlockIndex", 1);
		myTree->SetBranchStatus("ele_Etc", 1);
		myTree->SetBranchStatus("ele_E", 1);
		myTree->SetBranchStatus("ele_Px", 1);
		myTree->SetBranchStatus("ele_Py", 1);
		myTree->SetBranchStatus("ele_Pz", 1);
		myTree->SetBranchStatus("ele_Detector", 1);
		myTree->SetBranchStatus("ele_DetEta", 1);
		myTree->SetBranchStatus("ele_DetPhi", 1);
		myTree->SetBranchStatus("ele_XCes", 1);
		myTree->SetBranchStatus("ele_ZCes", 1);
		myTree->SetBranchStatus("ele_HadEm", 1);
		myTree->SetBranchStatus("ele_Chi2Mean", 1);
		myTree->SetBranchStatus("ele_N3d", 1);
		myTree->SetBranchStatus("ele_Iso4", 1);
		myTree->SetBranchStatus("ele_TrkIso", 1);		//from photon block
		myTree->SetBranchStatus("ele_CesWireE2", 1);
		myTree->SetBranchStatus("ele_CesStripE2", 1);
		myTree->SetBranchStatus("ele_PhiWedge", 1);
		myTree->SetBranchStatus("ele_NMuonStubs", 1);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchStatus("ele_EmTime", 1);
		myTree->SetBranchStatus("ele_PhoenixId", 1);
		myTree->SetBranchStatus("ele_Halo_seedWedge", 1);
		myTree->SetBranchStatus("ele_Halo_eastNhad", 1);
		myTree->SetBranchStatus("ele_Halo_westNhad", 1);
		myTree->SetBranchStatus("ele_Ntight", 1);
		myTree->SetBranchStatus("ele_Nloose", 1);
		myTree->SetBranchStatus("ele_Ntracks", 1);				//these are from TStnElectron
		myTree->SetBranchStatus("ele_Emfr", 1);
		myTree->SetBranchStatus("ele_EoverP", 1);
		myTree->SetBranchStatus("ele_TrackPt", 1);
		myTree->SetBranchStatus("ele_TrackBcPt", 1);
		myTree->SetBranchStatus("ele_TrackPhi", 1);
		myTree->SetBranchStatus("ele_Nssl", 1);
		myTree->SetBranchStatus("ele_Nasl", 1);
		myTree->SetBranchStatus("ele_TightId", 1);
		myTree->SetBranchStatus("ele_LooseId", 1);
		myTree->SetBranchStatus("ele_ConversionId", 1);
		myTree->SetBranchStatus("ele_matchJetIndex", 1);	//if a matching jet is found and removed from jet list


		myTree->SetBranchAddress("ele_num", &myStuple->ele_num);
		myTree->SetBranchAddress("ele_Ntight", &myStuple->ele_Ntight);
		myTree->SetBranchAddress("ele_Nloose", &myStuple->ele_Nloose);
		myTree->SetBranchAddress("ele_Ntracks", &myStuple->ele_Ntracks);				//these are from TStnElectron
		myTree->SetBranchAddress("ele_Index", &myStuple->ele_Index);
		myTree->SetBranchAddress("ele_EleBlockIndex", &myStuple->ele_EleBlockIndex);
		myTree->SetBranchAddress("ele_Etc", &myStuple->ele_Etc);
		myTree->SetBranchAddress("ele_E", &myStuple->ele_E);
		myTree->SetBranchAddress("ele_Px", &myStuple->ele_Px);
		myTree->SetBranchAddress("ele_Py", &myStuple->ele_Py);
		myTree->SetBranchAddress("ele_Pz", &myStuple->ele_Pz);
		myTree->SetBranchAddress("ele_Detector", &myStuple->ele_Detector);
		myTree->SetBranchAddress("ele_DetEta", &myStuple->ele_DetEta);
		myTree->SetBranchAddress("ele_DetPhi", &myStuple->ele_DetPhi);
		myTree->SetBranchAddress("ele_XCes", &myStuple->ele_XCes);
		myTree->SetBranchAddress("ele_ZCes", &myStuple->ele_ZCes);
		myTree->SetBranchAddress("ele_HadEm", &myStuple->ele_HadEm);
		myTree->SetBranchAddress("ele_Chi2Mean", &myStuple->ele_Chi2Mean);
		myTree->SetBranchAddress("ele_N3d", &myStuple->ele_N3d);
		myTree->SetBranchAddress("ele_Iso4", &myStuple->ele_Iso4);
		myTree->SetBranchAddress("ele_TrkIso", &myStuple->ele_TrkIso);		//from photon block
		myTree->SetBranchAddress("ele_CesWireE2", &myStuple->ele_CesWireE2);
		myTree->SetBranchAddress("ele_CesStripE2", &myStuple->ele_CesStripE2);
		myTree->SetBranchAddress("ele_PhiWedge", &myStuple->ele_PhiWedge);
		myTree->SetBranchAddress("ele_NMuonStubs", &myStuple->ele_NMuonStubs);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchAddress("ele_EmTime", &myStuple->ele_EmTime);
		myTree->SetBranchAddress("ele_PhoenixId", &myStuple->ele_PhoenixId);
		myTree->SetBranchAddress("ele_Halo_seedWedge", &myStuple->ele_Halo_seedWedge);
		myTree->SetBranchAddress("ele_Halo_eastNhad", &myStuple->ele_Halo_eastNhad);
		myTree->SetBranchAddress("ele_Halo_westNhad", &myStuple->ele_Halo_westNhad);
		myTree->SetBranchAddress("ele_Emfr", &myStuple->ele_Emfr);
		myTree->SetBranchAddress("ele_EoverP", &myStuple->ele_EoverP);
		myTree->SetBranchAddress("ele_TrackPt", &myStuple->ele_TrackPt);
		myTree->SetBranchAddress("ele_TrackBcPt", &myStuple->ele_TrackBcPt);
		myTree->SetBranchAddress("ele_TrackPhi", &myStuple->ele_TrackPhi);
		myTree->SetBranchAddress("ele_Nssl", &myStuple->ele_Nssl);
		myTree->SetBranchAddress("ele_Nasl", &myStuple->ele_Nasl);
		myTree->SetBranchAddress("ele_TightId", &myStuple->ele_TightId);
		myTree->SetBranchAddress("ele_LooseId", &myStuple->ele_LooseId);
		myTree->SetBranchAddress("ele_ConversionId", &myStuple->ele_ConversionId);
		myTree->SetBranchAddress("ele_matchJetIndex", &myStuple->ele_matchJetIndex);	//if a matching jet is found and removed from jet list

	}

	if (bAllEle || bEmUpEle) {
		std::cout << "Enabling EM UP Electrons Read-in." << std::endl;
		myTree->SetBranchStatus("ele_up_num", 1);
		myTree->SetBranchStatus("ele_up_Index", 1);
		myTree->SetBranchStatus("ele_up_EleBlockIndex", 1);
		myTree->SetBranchStatus("ele_up_Etc", 1);
		myTree->SetBranchStatus("ele_up_E", 1);
		myTree->SetBranchStatus("ele_up_Px", 1);
		myTree->SetBranchStatus("ele_up_Py", 1);
		myTree->SetBranchStatus("ele_up_Pz", 1);
		myTree->SetBranchStatus("ele_up_Detector", 1);
		myTree->SetBranchStatus("ele_up_DetEta", 1);
		myTree->SetBranchStatus("ele_up_DetPhi", 1);
		myTree->SetBranchStatus("ele_up_XCes", 1);
		myTree->SetBranchStatus("ele_up_ZCes", 1);
		myTree->SetBranchStatus("ele_up_HadEm", 1);
		myTree->SetBranchStatus("ele_up_Chi2Mean", 1);
		myTree->SetBranchStatus("ele_up_N3d", 1);
		myTree->SetBranchStatus("ele_up_Iso4", 1);
		myTree->SetBranchStatus("ele_up_TrkIso", 1);		//from photon block
		myTree->SetBranchStatus("ele_up_CesWireE2", 1);
		myTree->SetBranchStatus("ele_up_CesStripE2", 1);
		myTree->SetBranchStatus("ele_up_PhiWedge", 1);
		myTree->SetBranchStatus("ele_up_NMuonStubs", 1);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchStatus("ele_up_EmTime", 1);
		myTree->SetBranchStatus("ele_up_PhoenixId", 1);
		myTree->SetBranchStatus("ele_up_Halo_seedWedge", 1);
		myTree->SetBranchStatus("ele_up_Halo_eastNhad", 1);
		myTree->SetBranchStatus("ele_up_Halo_westNhad", 1);
		myTree->SetBranchStatus("ele_up_Ntight", 1);
		myTree->SetBranchStatus("ele_up_Nloose", 1);
		myTree->SetBranchStatus("ele_up_Ntracks", 1);				//these are from TStnElectron
		myTree->SetBranchStatus("ele_up_Emfr", 1);
		myTree->SetBranchStatus("ele_up_EoverP", 1);
		myTree->SetBranchStatus("ele_up_TrackPt", 1);
		myTree->SetBranchStatus("ele_up_TrackBcPt", 1);
		myTree->SetBranchStatus("ele_up_TrackPhi", 1);
		myTree->SetBranchStatus("ele_up_Nssl", 1);
		myTree->SetBranchStatus("ele_up_Nasl", 1);
		myTree->SetBranchStatus("ele_up_TightId", 1);
		myTree->SetBranchStatus("ele_up_LooseId", 1);
		myTree->SetBranchStatus("ele_up_ConversionId", 1);
		myTree->SetBranchStatus("ele_up_matchJetIndex", 1);	//if a matching jet is found and removed from jet list


		myTree->SetBranchAddress("ele_up_num", &myStuple->ele_up_num);
		myTree->SetBranchAddress("ele_up_Ntight", &myStuple->ele_up_Ntight);
		myTree->SetBranchAddress("ele_up_Nloose", &myStuple->ele_up_Nloose);
		myTree->SetBranchAddress("ele_up_Ntracks", &myStuple->ele_up_Ntracks);				//these are from TStnElectron
		myTree->SetBranchAddress("ele_up_Index", &myStuple->ele_up_Index);
		myTree->SetBranchAddress("ele_up_EleBlockIndex", &myStuple->ele_up_EleBlockIndex);
		myTree->SetBranchAddress("ele_up_Etc", &myStuple->ele_up_Etc);
		myTree->SetBranchAddress("ele_up_E", &myStuple->ele_up_E);
		myTree->SetBranchAddress("ele_up_Px", &myStuple->ele_up_Px);
		myTree->SetBranchAddress("ele_up_Py", &myStuple->ele_up_Py);
		myTree->SetBranchAddress("ele_up_Pz", &myStuple->ele_up_Pz);
		myTree->SetBranchAddress("ele_up_Detector", &myStuple->ele_up_Detector);
		myTree->SetBranchAddress("ele_up_DetEta", &myStuple->ele_up_DetEta);
		myTree->SetBranchAddress("ele_up_DetPhi", &myStuple->ele_up_DetPhi);
		myTree->SetBranchAddress("ele_up_XCes", &myStuple->ele_up_XCes);
		myTree->SetBranchAddress("ele_up_ZCes", &myStuple->ele_up_ZCes);
		myTree->SetBranchAddress("ele_up_HadEm", &myStuple->ele_up_HadEm);
		myTree->SetBranchAddress("ele_up_Chi2Mean", &myStuple->ele_up_Chi2Mean);
		myTree->SetBranchAddress("ele_up_N3d", &myStuple->ele_up_N3d);
		myTree->SetBranchAddress("ele_up_Iso4", &myStuple->ele_up_Iso4);
		myTree->SetBranchAddress("ele_up_TrkIso", &myStuple->ele_up_TrkIso);		//from photon block
		myTree->SetBranchAddress("ele_up_CesWireE2", &myStuple->ele_up_CesWireE2);
		myTree->SetBranchAddress("ele_up_CesStripE2", &myStuple->ele_up_CesStripE2);
		myTree->SetBranchAddress("ele_up_PhiWedge", &myStuple->ele_up_PhiWedge);
		myTree->SetBranchAddress("ele_up_NMuonStubs", &myStuple->ele_up_NMuonStubs);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchAddress("ele_up_EmTime", &myStuple->ele_up_EmTime);
		myTree->SetBranchAddress("ele_up_PhoenixId", &myStuple->ele_up_PhoenixId);
		myTree->SetBranchAddress("ele_up_Halo_seedWedge", &myStuple->ele_up_Halo_seedWedge);
		myTree->SetBranchAddress("ele_up_Halo_eastNhad", &myStuple->ele_up_Halo_eastNhad);
		myTree->SetBranchAddress("ele_up_Halo_westNhad", &myStuple->ele_up_Halo_westNhad);
		myTree->SetBranchAddress("ele_up_Emfr", &myStuple->ele_up_Emfr);
		myTree->SetBranchAddress("ele_up_EoverP", &myStuple->ele_up_EoverP);
		myTree->SetBranchAddress("ele_up_TrackPt", &myStuple->ele_up_TrackPt);
		myTree->SetBranchAddress("ele_up_TrackBcPt", &myStuple->ele_up_TrackBcPt);
		myTree->SetBranchAddress("ele_up_TrackPhi", &myStuple->ele_up_TrackPhi);
		myTree->SetBranchAddress("ele_up_Nssl", &myStuple->ele_up_Nssl);
		myTree->SetBranchAddress("ele_up_Nasl", &myStuple->ele_up_Nasl);
		myTree->SetBranchAddress("ele_up_TightId", &myStuple->ele_up_TightId);
		myTree->SetBranchAddress("ele_up_LooseId", &myStuple->ele_up_LooseId);
		myTree->SetBranchAddress("ele_up_ConversionId", &myStuple->ele_up_ConversionId);
		myTree->SetBranchAddress("ele_up_matchJetIndex", &myStuple->ele_up_matchJetIndex);	//if a matching jet is found and removed from jet list


	}

	if (bAllEle || bEmDownEle) {
		std::cout << "Enabling EM Down Electrons Read-in." << std::endl;
		myTree->SetBranchStatus("ele_down_num", 1);
		myTree->SetBranchStatus("ele_down_Index", 1);
		myTree->SetBranchStatus("ele_down_EleBlockIndex", 1);
		myTree->SetBranchStatus("ele_down_Etc", 1);
		myTree->SetBranchStatus("ele_down_E", 1);
		myTree->SetBranchStatus("ele_down_Px", 1);
		myTree->SetBranchStatus("ele_down_Py", 1);
		myTree->SetBranchStatus("ele_down_Pz", 1);
		myTree->SetBranchStatus("ele_down_Detector", 1);
		myTree->SetBranchStatus("ele_down_DetEta", 1);
		myTree->SetBranchStatus("ele_down_DetPhi", 1);
		myTree->SetBranchStatus("ele_down_XCes", 1);
		myTree->SetBranchStatus("ele_down_ZCes", 1);
		myTree->SetBranchStatus("ele_down_HadEm", 1);
		myTree->SetBranchStatus("ele_down_Chi2Mean", 1);
		myTree->SetBranchStatus("ele_down_N3d", 1);
		myTree->SetBranchStatus("ele_down_Iso4", 1);
		myTree->SetBranchStatus("ele_down_TrkIso", 1);		//from photon block
		myTree->SetBranchStatus("ele_down_CesWireE2", 1);
		myTree->SetBranchStatus("ele_down_CesStripE2", 1);
		myTree->SetBranchStatus("ele_down_PhiWedge", 1);
		myTree->SetBranchStatus("ele_down_NMuonStubs", 1);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchStatus("ele_down_EmTime", 1);
		myTree->SetBranchStatus("ele_down_PhoenixId", 1);
		myTree->SetBranchStatus("ele_down_Halo_seedWedge", 1);
		myTree->SetBranchStatus("ele_down_Halo_eastNhad", 1);
		myTree->SetBranchStatus("ele_down_Halo_westNhad", 1);
		myTree->SetBranchStatus("ele_down_Ntight", 1);
		myTree->SetBranchStatus("ele_down_Nloose", 1);
		myTree->SetBranchStatus("ele_down_Ntracks", 1);				//these are from TStnElectron
		myTree->SetBranchStatus("ele_down_Emfr", 1);
		myTree->SetBranchStatus("ele_down_EoverP", 1);
		myTree->SetBranchStatus("ele_down_TrackPt", 1);
		myTree->SetBranchStatus("ele_down_TrackBcPt", 1);
		myTree->SetBranchStatus("ele_down_TrackPhi", 1);
		myTree->SetBranchStatus("ele_down_Nssl", 1);
		myTree->SetBranchStatus("ele_down_Nasl", 1);
		myTree->SetBranchStatus("ele_down_TightId", 1);
		myTree->SetBranchStatus("ele_down_LooseId", 1);
		myTree->SetBranchStatus("ele_down_ConversionId", 1);
		myTree->SetBranchStatus("ele_down_matchJetIndex", 1);	//if a matching jet is found and removed from jet list


		myTree->SetBranchAddress("ele_down_num", &myStuple->ele_down_num);
		myTree->SetBranchAddress("ele_down_Ntight", &myStuple->ele_down_Ntight);
		myTree->SetBranchAddress("ele_down_Nloose", &myStuple->ele_down_Nloose);
		myTree->SetBranchAddress("ele_down_Ntracks", &myStuple->ele_down_Ntracks);				//these are from TStnElectron
		myTree->SetBranchAddress("ele_down_Index", &myStuple->ele_down_Index);
		myTree->SetBranchAddress("ele_down_EleBlockIndex", &myStuple->ele_down_EleBlockIndex);
		myTree->SetBranchAddress("ele_down_Etc", &myStuple->ele_down_Etc);
		myTree->SetBranchAddress("ele_down_E", &myStuple->ele_down_E);
		myTree->SetBranchAddress("ele_down_Px", &myStuple->ele_down_Px);
		myTree->SetBranchAddress("ele_down_Py", &myStuple->ele_down_Py);
		myTree->SetBranchAddress("ele_down_Pz", &myStuple->ele_down_Pz);
		myTree->SetBranchAddress("ele_down_Detector", &myStuple->ele_down_Detector);
		myTree->SetBranchAddress("ele_down_DetEta", &myStuple->ele_down_DetEta);
		myTree->SetBranchAddress("ele_down_DetPhi", &myStuple->ele_down_DetPhi);
		myTree->SetBranchAddress("ele_down_XCes", &myStuple->ele_down_XCes);
		myTree->SetBranchAddress("ele_down_ZCes", &myStuple->ele_down_ZCes);
		myTree->SetBranchAddress("ele_down_HadEm", &myStuple->ele_down_HadEm);
		myTree->SetBranchAddress("ele_down_Chi2Mean", &myStuple->ele_down_Chi2Mean);
		myTree->SetBranchAddress("ele_down_N3d", &myStuple->ele_down_N3d);
		myTree->SetBranchAddress("ele_down_Iso4", &myStuple->ele_down_Iso4);
		myTree->SetBranchAddress("ele_down_TrkIso", &myStuple->ele_down_TrkIso);		//from photon block
		myTree->SetBranchAddress("ele_down_CesWireE2", &myStuple->ele_down_CesWireE2);
		myTree->SetBranchAddress("ele_down_CesStripE2", &myStuple->ele_down_CesStripE2);
		myTree->SetBranchAddress("ele_down_PhiWedge", &myStuple->ele_down_PhiWedge);
		myTree->SetBranchAddress("ele_down_NMuonStubs", &myStuple->ele_down_NMuonStubs);	//trkless muon stubs around 30 degree cone of the photon
		myTree->SetBranchAddress("ele_down_EmTime", &myStuple->ele_down_EmTime);
		myTree->SetBranchAddress("ele_down_PhoenixId", &myStuple->ele_down_PhoenixId);
		myTree->SetBranchAddress("ele_down_Halo_seedWedge", &myStuple->ele_down_Halo_seedWedge);
		myTree->SetBranchAddress("ele_down_Halo_eastNhad", &myStuple->ele_down_Halo_eastNhad);
		myTree->SetBranchAddress("ele_down_Halo_westNhad", &myStuple->ele_down_Halo_westNhad);
		myTree->SetBranchAddress("ele_down_Emfr", &myStuple->ele_down_Emfr);
		myTree->SetBranchAddress("ele_down_EoverP", &myStuple->ele_down_EoverP);
		myTree->SetBranchAddress("ele_down_TrackPt", &myStuple->ele_down_TrackPt);
		myTree->SetBranchAddress("ele_down_TrackBcPt", &myStuple->ele_down_TrackBcPt);
		myTree->SetBranchAddress("ele_down_TrackPhi", &myStuple->ele_down_TrackPhi);
		myTree->SetBranchAddress("ele_down_Nssl", &myStuple->ele_down_Nssl);
		myTree->SetBranchAddress("ele_down_Nasl", &myStuple->ele_down_Nasl);
		myTree->SetBranchAddress("ele_down_TightId", &myStuple->ele_down_TightId);
		myTree->SetBranchAddress("ele_down_LooseId", &myStuple->ele_down_LooseId);
		myTree->SetBranchAddress("ele_down_ConversionId", &myStuple->ele_down_ConversionId);
		myTree->SetBranchAddress("ele_down_matchJetIndex", &myStuple->ele_down_matchJetIndex);	//if a matching jet is found and removed from jet list


	}


	if (bAllJet || bCentralJet) {
		std::cout << "Enabling Central Jets to Read-in." << std::endl;
		myTree->SetBranchStatus("jet_num", 1);
		myTree->SetBranchStatus("jet_NJet15", 1);
		myTree->SetBranchStatus("jet_Index", 1);				//original index in jet block
		myTree->SetBranchStatus("jet_Pt", 1);					//easy to get from JetFilter
		myTree->SetBranchStatus("jet_E", 1);
		myTree->SetBranchStatus("jet_Px", 1);
		myTree->SetBranchStatus("jet_Py", 1);
		myTree->SetBranchStatus("jet_Pz", 1);
		myTree->SetBranchStatus("jet_DetEta", 1);
		myTree->SetBranchStatus("jet_DetPhi", 1);
		myTree->SetBranchStatus("jet_HadEm", 1);
		myTree->SetBranchStatus("jet_Emfr", 1);
		myTree->SetBranchStatus("jet_Ntowers", 1);
		myTree->SetBranchStatus("jet_Ntracks", 1);
		myTree->SetBranchStatus("jet_SeedIPhi", 1);
		myTree->SetBranchStatus("jet_SeedIEta", 1);
		myTree->SetBranchStatus("jet_EmTime", 1);
		myTree->SetBranchStatus("jet_SecVtxTag", 1);
		myTree->SetBranchStatus("jet_SecVtxppb", 1);
		myTree->SetBranchStatus("jet_SecVtxnpb", 1);
		myTree->SetBranchStatus("jet_SecVtxTrkmass", 1);
		myTree->SetBranchStatus("jet_raw_num", 1);
		myTree->SetBranchStatus("jet_raw_Index", 1);				//original index in jet block
		myTree->SetBranchStatus("jet_raw_Pt", 1);
		myTree->SetBranchStatus("jet_raw_E", 1);
		myTree->SetBranchStatus("jet_raw_Px", 1);
		myTree->SetBranchStatus("jet_raw_Py", 1);
		myTree->SetBranchStatus("jet_raw_Pz", 1);

		
		myTree->SetBranchAddress("jet_num", &myStuple->jet_num);
		myTree->SetBranchAddress("jet_NJet15", &myStuple->jet_NJet15);
		myTree->SetBranchAddress("jet_Index", &myStuple->jet_Index);				//original index in jet block
		myTree->SetBranchAddress("jet_Pt", &myStuple->jet_Pt);					//easy to get from JetFilter
		myTree->SetBranchAddress("jet_E", &myStuple->jet_E);
		myTree->SetBranchAddress("jet_Px", &myStuple->jet_Px);
		myTree->SetBranchAddress("jet_Py", &myStuple->jet_Py);
		myTree->SetBranchAddress("jet_Pz", &myStuple->jet_Pz);
		myTree->SetBranchAddress("jet_DetEta", &myStuple->jet_DetEta);
		myTree->SetBranchAddress("jet_DetPhi", &myStuple->jet_DetPhi);
		myTree->SetBranchAddress("jet_HadEm", &myStuple->jet_HadEm);
		myTree->SetBranchAddress("jet_Emfr", &myStuple->jet_Emfr);
		myTree->SetBranchAddress("jet_Ntowers", &myStuple->jet_Ntowers);
		myTree->SetBranchAddress("jet_Ntracks", &myStuple->jet_Ntracks);
		myTree->SetBranchAddress("jet_SeedIPhi", &myStuple->jet_SeedIPhi);
		myTree->SetBranchAddress("jet_SeedIEta", &myStuple->jet_SeedIEta);
		myTree->SetBranchAddress("jet_EmTime", &myStuple->jet_EmTime);
		myTree->SetBranchAddress("jet_SecVtxTag", &myStuple->jet_SecVtxTag);
		myTree->SetBranchAddress("jet_SecVtxppb", &myStuple->jet_SecVtxppb);
		myTree->SetBranchAddress("jet_SecVtxnpb", &myStuple->jet_SecVtxnpb);
		myTree->SetBranchAddress("jet_SecVtxTrkmass", &myStuple->jet_SecVtxTrkmass);

		myTree->SetBranchAddress("jet_raw_num", &myStuple->jet_raw_num);
		myTree->SetBranchAddress("jet_raw_Index", &myStuple->jet_raw_Index);				//original index in jet block
		myTree->SetBranchAddress("jet_raw_Pt", &myStuple->jet_raw_Pt);
		myTree->SetBranchAddress("jet_raw_E", &myStuple->jet_raw_E);
		myTree->SetBranchAddress("jet_raw_Px", &myStuple->jet_raw_Px);
		myTree->SetBranchAddress("jet_raw_Py", &myStuple->jet_raw_Py);
		myTree->SetBranchAddress("jet_raw_Pz", &myStuple->jet_raw_Pz);

	}

	if (bAllJet || bJesUpJet) {
		std::cout << "Enabling JES Up Jets to Read-in." << std::endl;
		myTree->SetBranchStatus("jet_up_num", 1);
		myTree->SetBranchStatus("jet_up_NJet15", 1);
		myTree->SetBranchStatus("jet_up_Index", 1);				//original index in jet block
		myTree->SetBranchStatus("jet_up_Pt", 1);					//easy to get from JetFilter
		myTree->SetBranchStatus("jet_up_E", 1);
		myTree->SetBranchStatus("jet_up_Px", 1);
		myTree->SetBranchStatus("jet_up_Py", 1);
		myTree->SetBranchStatus("jet_up_Pz", 1);
		myTree->SetBranchStatus("jet_up_DetEta", 1);
		myTree->SetBranchStatus("jet_up_DetPhi", 1);
		myTree->SetBranchStatus("jet_up_HadEm", 1);
		myTree->SetBranchStatus("jet_up_Emfr", 1);
		myTree->SetBranchStatus("jet_up_Ntowers", 1);
		myTree->SetBranchStatus("jet_up_Ntracks", 1);
		myTree->SetBranchStatus("jet_up_SeedIPhi", 1);
		myTree->SetBranchStatus("jet_up_SeedIEta", 1);
		myTree->SetBranchStatus("jet_up_EmTime", 1);
		myTree->SetBranchStatus("jet_up_SecVtxTag", 1);
		myTree->SetBranchStatus("jet_up_SecVtxppb", 1);
		myTree->SetBranchStatus("jet_up_SecVtxnpb", 1);
		myTree->SetBranchStatus("jet_up_SecVtxTrkmass", 1);

		myTree->SetBranchAddress("jet_up_num", &myStuple->jet_up_num);
		myTree->SetBranchAddress("jet_up_NJet15", &myStuple->jet_up_NJet15);
		myTree->SetBranchAddress("jet_up_Index", &myStuple->jet_up_Index);				//original index in jet block
		myTree->SetBranchAddress("jet_up_Pt", &myStuple->jet_up_Pt);					//easy to get from JetFilter
		myTree->SetBranchAddress("jet_up_E", &myStuple->jet_up_E);
		myTree->SetBranchAddress("jet_up_Px", &myStuple->jet_up_Px);
		myTree->SetBranchAddress("jet_up_Py", &myStuple->jet_up_Py);
		myTree->SetBranchAddress("jet_up_Pz", &myStuple->jet_up_Pz);
		myTree->SetBranchAddress("jet_up_DetEta", &myStuple->jet_up_DetEta);
		myTree->SetBranchAddress("jet_up_DetPhi", &myStuple->jet_up_DetPhi);
		myTree->SetBranchAddress("jet_up_HadEm", &myStuple->jet_up_HadEm);
		myTree->SetBranchAddress("jet_up_Emfr", &myStuple->jet_up_Emfr);
		myTree->SetBranchAddress("jet_up_Ntowers", &myStuple->jet_up_Ntowers);
		myTree->SetBranchAddress("jet_up_Ntracks", &myStuple->jet_up_Ntracks);
		myTree->SetBranchAddress("jet_up_SeedIPhi", &myStuple->jet_up_SeedIPhi);
		myTree->SetBranchAddress("jet_up_SeedIEta", &myStuple->jet_up_SeedIEta);
		myTree->SetBranchAddress("jet_up_EmTime", &myStuple->jet_up_EmTime);
		myTree->SetBranchAddress("jet_up_SecVtxTag", &myStuple->jet_up_SecVtxTag);
		myTree->SetBranchAddress("jet_up_SecVtxppb", &myStuple->jet_up_SecVtxppb);
		myTree->SetBranchAddress("jet_up_SecVtxnpb", &myStuple->jet_up_SecVtxnpb);
		myTree->SetBranchAddress("jet_up_SecVtxTrkmass", &myStuple->jet_up_SecVtxTrkmass);
	}

	if (bAllJet || bJesDownJet) {
		std::cout << "Enabling JES Down Jets to Read-in." << std::endl;
		myTree->SetBranchStatus("jet_down_num", 1);
		myTree->SetBranchStatus("jet_down_NJet15", 1);
		myTree->SetBranchStatus("jet_down_Index", 1);				//original index in jet block
		myTree->SetBranchStatus("jet_down_Pt", 1);					//easy to get from JetFilter
		myTree->SetBranchStatus("jet_down_E", 1);
		myTree->SetBranchStatus("jet_down_Px", 1);
		myTree->SetBranchStatus("jet_down_Py", 1);
		myTree->SetBranchStatus("jet_down_Pz", 1);
		myTree->SetBranchStatus("jet_down_DetEta", 1);
		myTree->SetBranchStatus("jet_down_DetPhi", 1);
		myTree->SetBranchStatus("jet_down_HadEm", 1);
		myTree->SetBranchStatus("jet_down_Emfr", 1);
		myTree->SetBranchStatus("jet_down_Ntowers", 1);
		myTree->SetBranchStatus("jet_down_Ntracks", 1);
		myTree->SetBranchStatus("jet_down_SeedIPhi", 1);
		myTree->SetBranchStatus("jet_down_SeedIEta", 1);
		myTree->SetBranchStatus("jet_down_EmTime", 1);
		myTree->SetBranchStatus("jet_down_SecVtxTag", 1);
		myTree->SetBranchStatus("jet_down_SecVtxppb", 1);
		myTree->SetBranchStatus("jet_down_SecVtxnpb", 1);
		myTree->SetBranchStatus("jet_down_SecVtxTrkmass", 1);

		myTree->SetBranchAddress("jet_down_num", &myStuple->jet_down_num);
		myTree->SetBranchAddress("jet_down_NJet15", &myStuple->jet_down_NJet15);
		myTree->SetBranchAddress("jet_down_Index", &myStuple->jet_down_Index);				//original index in jet block
		myTree->SetBranchAddress("jet_down_Pt", &myStuple->jet_down_Pt);					//easy to get from JetFilter
		myTree->SetBranchAddress("jet_down_E", &myStuple->jet_down_E);
		myTree->SetBranchAddress("jet_down_Px", &myStuple->jet_down_Px);
		myTree->SetBranchAddress("jet_down_Py", &myStuple->jet_down_Py);
		myTree->SetBranchAddress("jet_down_Pz", &myStuple->jet_down_Pz);
		myTree->SetBranchAddress("jet_down_DetEta", &myStuple->jet_down_DetEta);
		myTree->SetBranchAddress("jet_down_DetPhi", &myStuple->jet_down_DetPhi);
		myTree->SetBranchAddress("jet_down_HadEm", &myStuple->jet_down_HadEm);
		myTree->SetBranchAddress("jet_down_Emfr", &myStuple->jet_down_Emfr);
		myTree->SetBranchAddress("jet_down_Ntowers", &myStuple->jet_down_Ntowers);
		myTree->SetBranchAddress("jet_down_Ntracks", &myStuple->jet_down_Ntracks);
		myTree->SetBranchAddress("jet_down_SeedIPhi", &myStuple->jet_down_SeedIPhi);
		myTree->SetBranchAddress("jet_down_SeedIEta", &myStuple->jet_down_SeedIEta);
		myTree->SetBranchAddress("jet_down_EmTime", &myStuple->jet_down_EmTime);
		myTree->SetBranchAddress("jet_down_SecVtxTag", &myStuple->jet_down_SecVtxTag);
		myTree->SetBranchAddress("jet_down_SecVtxppb", &myStuple->jet_down_SecVtxppb);
		myTree->SetBranchAddress("jet_down_SecVtxnpb", &myStuple->jet_down_SecVtxnpb);
		myTree->SetBranchAddress("jet_down_SecVtxTrkmass", &myStuple->jet_down_SecVtxTrkmass);
	}





	if (bMet) {
		std::cout << "Enabling MET to Read-in." << std::endl;
		
		myTree->SetBranchStatus("met_Met", 1);				// corrected MET
		myTree->SetBranchStatus("met_RawMet", 1);				// corrected MET
		myTree->SetBranchStatus("met_MetX", 1);				// corrected MET X 
		myTree->SetBranchStatus("met_MetY", 1);				// corrected MET Y 
		myTree->SetBranchStatus("met_SumEt", 1);			// corrected
		myTree->SetBranchStatus("met_Ht", 1);				// corrected
		myTree->SetBranchStatus("met_MetPhi", 1);			// corrected
//		myTree->SetBranchStatus("met_Gen_d", 1);
//		myTree->SetBranchStatus("met_Gen_m", 1);
//		myTree->SetBranchStatus("met_Gen_p", 1);
//		myTree->SetBranchStatus("met_Gen_mUn", 1);
//		myTree->SetBranchStatus("met_Gen_pUn", 1);

		myTree->SetBranchAddress("met_Met", &myStuple->met_Met);				// corrected MET
		myTree->SetBranchAddress("met_RawMet", &myStuple->met_RawMet);		// uncorrecte MET
		myTree->SetBranchAddress("met_MetX", &myStuple->met_MetX);				// corrected MET X
		myTree->SetBranchAddress("met_MetY", &myStuple->met_MetY);				// corrected MET Y
		myTree->SetBranchAddress("met_SumEt", &myStuple->met_SumEt);			// corrected
		myTree->SetBranchAddress("met_Ht", &myStuple->met_Ht);				// corrected
		myTree->SetBranchAddress("met_MetPhi", &myStuple->met_MetPhi);			// corrected
		
		//myTree->SetBranchAddress("met_Gen_d", &myStuple->met_Gen_d);
		//myTree->SetBranchAddress("met_Gen_m", &myStuple->met_Gen_m);
		//myTree->SetBranchAddress("met_Gen_p", &myStuple->met_Gen_p);
		//myTree->SetBranchAddress("met_Gen_mUn", &myStuple->met_Gen_mUn);
		//myTree->SetBranchAddress("met_Gen_pUn", &myStuple->met_Gen_pUn);
		
	}
	std::cout << "pho_num status=" << myTree->GetBranchStatus("pho_num") << std::endl;
	std::cout << "pho_up_num status=" << myTree->GetBranchStatus("pho_up_num") << std::endl;
	std::cout << "pho_down_num status=" << myTree->GetBranchStatus("pho_down_num") << std::endl;

}

