#include "PhoFlat/Stuple.hh"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "TChain.h"
//created this to make a simple ntuple from my BIG FAT ntuple 
//for an undergrad as requested by Jay-06-22-2009

void MakeASmallerStuple(const int nevts=100, const bool debug = false)
{
	Stuple *myStuple = new Stuple;

	//this is the new smaller ntuple. need to create the tree
	//from scratch for now
	//TFile *rootFile = new TFile("~/STUPLES/ForUndergrads/PHOJETS_DATA_NTUPLE.root","RECREATE");
	TFile *rootFile = new TFile("~/STUPLES/ForUndergrads/PHOJETS_MC_NTUPLE.root","RECREATE");
	TTree *newTree = new TTree("Stuple","Sam's flat Stuple.");

		newTree->Branch("evt_McFlag"   	 , &myStuple->evt_McFlag   	, "evt_McFlag/I");
		newTree->Branch("evt_RunNumber"   , &myStuple->evt_RunNumber   , "evt_RunNumber/i");
		newTree->Branch("evt_EventNumber" , &myStuple->evt_EventNumber , "ect_EventNumber/i");

		newTree->Branch("tri_pho25iso" , &myStuple->tri_pho25iso , "tri_pho25iso/I");
		newTree->Branch("tri_pho50"    , &myStuple->tri_pho50    , "tri_pho50/I");
		newTree->Branch("tri_pho70"    , &myStuple->tri_pho70    , "tri_pho70/I");
	
		newTree->Branch("vtx_N", &myStuple->vtx_N, "vtx_N/I");
		newTree->Branch("vtx_NClass12", &myStuple->vtx_NClass12, "vtx_NClass12/I");
		newTree->Branch("vtx_z", &myStuple->vtx_z, "vtx_z/F");				// highest sum pt vertex Z
		newTree->Branch("vtx_Ntracks", &myStuple->vtx_Ntracks, "vtx_Ntracks/I");		// in highest sum pt vertex
		newTree->Branch("vtx_SumPt", &myStuple->vtx_SumPt,"vtx_SumPt/F");			// highest

	
		newTree->Branch("pho_num"     , &myStuple->pho_num      , "pho_num/i");
		newTree->Branch("pho_Ntight"  , &myStuple->pho_Ntight   , "pho_Ntight/i");
		newTree->Branch("pho_Nloose"  , &myStuple->pho_Nloose   , "pho_Nloose/i");
		newTree->Branch("pho_Index"   , &myStuple->pho_Index    , "pho_Index[pho_num]/I");
		newTree->Branch("pho_PhoBlockIndex"   , &myStuple->pho_PhoBlockIndex    , "pho_PhoBlockIndex[pho_num]/I");
		newTree->Branch("pho_Etc"     , &myStuple->pho_Etc      , "pho_Etc[pho_num]/F");
		newTree->Branch("pho_E"       , &myStuple->pho_E        , "pho_E[pho_num]/F");
		newTree->Branch("pho_Px"      , &myStuple->pho_Px       , "pho_Px[pho_num]/F");
		newTree->Branch("pho_Py"      , &myStuple->pho_Py       , "pho_Py[pho_num]/F");
		newTree->Branch("pho_Pz"      , &myStuple->pho_Pz       , "pho_Pz[pho_num]/F");
		newTree->Branch("pho_Detector", &myStuple->pho_Detector , "pho_Detector[pho_num]/I");
		newTree->Branch("pho_DetEta"  , &myStuple->pho_DetEta   , "pho_DetEta[pho_num]/F");
		newTree->Branch("pho_DetPhi"  , &myStuple->pho_DetPhi   , "pho_DetPhi[pho_num]/F");
		newTree->Branch("pho_XCes"    , &myStuple->pho_XCes     , "pho_XCes[pho_num]/F");
		newTree->Branch("pho_ZCes"    , &myStuple->pho_ZCes     , "pho_ZCes[pho_num]/F");
		newTree->Branch("pho_HadEm"   , &myStuple->pho_HadEm    , "pho_HadEm[pho_num]/F");
		newTree->Branch("pho_Chi2Mean", &myStuple->pho_Chi2Mean , "pho_Chi2Mean[pho_num]/F");
		newTree->Branch("pho_N3d"     , &myStuple->pho_N3d      , "pho_N3d[pho_num]/I");
		newTree->Branch("pho_Iso4"    , &myStuple->pho_Iso4     , "pho_Iso4[pho_num]/F");
		newTree->Branch("pho_TrkPt"   , &myStuple->pho_TrkPt    , "pho_TrkPt[pho_num]/F");
		newTree->Branch("pho_TrkIso"  , &myStuple->pho_TrkIso   , "pho_TrkIso[pho_num]/F");
		newTree->Branch("pho_CesWireE2"      , &myStuple->pho_CesWireE2      , "pho_CesWireE2[pho_num]/F");
		newTree->Branch("pho_CesStripE2"     , &myStuple->pho_CesStripE2     , "pho_CesStripE2[pho_num]/F");
		newTree->Branch("pho_PhiWedge"       , &myStuple->pho_PhiWedge       , "pho_PhiWedge[pho_num]/I");
		newTree->Branch("pho_NMuonStubs"     , &myStuple->pho_NMuonStubs     , "pho_NMuonStubs[pho_num]/I");
		newTree->Branch("pho_EmTime"         , &myStuple->pho_EmTime         , "pho_EmTime[pho_num]/F");
		newTree->Branch("pho_TightId"        , &myStuple->pho_TightId        , "pho_TightId[pho_num]/I");
		newTree->Branch("pho_LooseId"        , &myStuple->pho_LooseId        , "pho_LooseId[pho_num]/I");

		newTree->Branch("jet_num"      , &myStuple->jet_num      , "jet_num/i");
		newTree->Branch("jet_NJet15"   , &myStuple->jet_NJet15   , "jet_Njet15/I");
		newTree->Branch("jet_Index"    , &myStuple->jet_Index    , "jet_Index[jet_num]/I");
		newTree->Branch("jet_Pt"       , &myStuple->jet_Pt       , "jet_Pt[jet_num]/F");
		newTree->Branch("jet_E"        , &myStuple->jet_E        , "jet_E[jet_num]/F");
		newTree->Branch("jet_Px"       , &myStuple->jet_Px       , "jet_Px[jet_num]/F");
		newTree->Branch("jet_Py"       , &myStuple->jet_Py       , "jet_Py[jet_num]/F");
		newTree->Branch("jet_Pz"       , &myStuple->jet_Pz       , "jet_Pz[jet_num]/F");
		newTree->Branch("jet_DetEta"   , &myStuple->jet_DetEta   , "jet_DetEta[jet_num]/F");
		newTree->Branch("jet_DetPhi"   , &myStuple->jet_DetPhi   , "jet_DetPhi[jet_num]/F");
		newTree->Branch("jet_HadEm"    , &myStuple->jet_HadEm    , "jet_HadEm[jet_num]/F");
		newTree->Branch("jet_Emfr"     , &myStuple->jet_Emfr     , "jet_Emfr[jet_num]/F");
		newTree->Branch("jet_Ntowers"  , &myStuple->jet_Ntowers  , "jet_Ntowers[jet_num]/I");
		newTree->Branch("jet_Ntracks"  , &myStuple->jet_Ntracks  , "jet_Ntracks[jet_num]/I");
		newTree->Branch("jet_SeedIPhi" , &myStuple->jet_SeedIPhi , "jet_SeedIPhi[jet_num]/I");
		newTree->Branch("jet_SeedIEta" , &myStuple->jet_SeedIEta , "jet_SeedIEta[jet_num]/I");



	//this is the BIG newTree (ntuple)
	//add all the files to a chain
	TChain *oldTree = new TChain("Stuple"); 
//	oldTree->Add("/data/nbay02/a/samantha/STUPLES/StupleV5_PhoData_1_1.root");
//	oldTree->Add("/data/nbay02/a/samantha/STUPLES/StupleV5_PhoData_2_0.root");
//	oldTree->Add("/data/nbay02/a/samantha/STUPLES/StupleV5_PhoData_2_1.root");
//	oldTree->Add("/data/nbay02/a/samantha/STUPLES/StupleV5_PhoData_3_0.root");
//	oldTree->Add("/data/nbay02/a/samantha/STUPLES/StupleV5_PhoData_3_1.root");
	oldTree->Add("/data/nbay02/a/samantha/STUPLES/StupleV4_PhoMC.root");
	
	std::cout << "Entries found = " << oldTree->GetEntries() << std::endl;
	oldTree->SetBranchStatus("*", 0);
	oldTree->SetBranchStatus("evt_RunNumber", 1);
	oldTree->SetBranchStatus("evt_EventNumber", 1);
	oldTree->SetBranchStatus("evt_McFlag", 1);
	oldTree->SetBranchStatus("tri_pho25iso", 1);
	oldTree->SetBranchStatus("tri_pho50", 1);
	oldTree->SetBranchStatus("tri_pho70", 1);

	oldTree->SetBranchAddress("evt_RunNumber", &myStuple->evt_RunNumber);
	oldTree->SetBranchAddress("evt_EventNumber", &myStuple->evt_EventNumber);
	oldTree->SetBranchAddress("evt_McFlag", &myStuple->evt_McFlag);
	oldTree->SetBranchAddress("tri_pho25iso", &myStuple->tri_pho25iso);
	oldTree->SetBranchAddress("tri_pho50", &myStuple->tri_pho50);
	oldTree->SetBranchAddress("tri_pho70", &myStuple->tri_pho70);
	
	oldTree->SetBranchStatus("vtx_N", 1);
	oldTree->SetBranchStatus("vtx_NClass12", 1);
	oldTree->SetBranchStatus("vtx_z", 1);				// highest sum pt vertex Z
	oldTree->SetBranchStatus("vtx_Ntracks", 1);		// in highest sum pt vertex
	oldTree->SetBranchStatus("vtx_SumPt", 1);			// highest

	oldTree->SetBranchAddress("vtx_N", &myStuple->vtx_N);
	oldTree->SetBranchAddress("vtx_NClass12", &myStuple->vtx_NClass12);
	oldTree->SetBranchAddress("vtx_z", &myStuple->vtx_z);				// highest sum pt vertex Z
	oldTree->SetBranchAddress("vtx_Ntracks", &myStuple->vtx_Ntracks);		// in highest sum pt vertex
	oldTree->SetBranchAddress("vtx_SumPt", &myStuple->vtx_SumPt);			// highest

	
	oldTree->SetBranchStatus("pho_num", 1);
	oldTree->SetBranchStatus("pho_Ntight", 1);
	oldTree->SetBranchStatus("pho_Nloose", 1);
	oldTree->SetBranchStatus("pho_Index", 1);
	oldTree->SetBranchStatus("pho_Etc", 1);
	oldTree->SetBranchStatus("pho_E", 1);
	oldTree->SetBranchStatus("pho_Px", 1);
	oldTree->SetBranchStatus("pho_Py", 1);
	oldTree->SetBranchStatus("pho_Pz", 1);
	oldTree->SetBranchStatus("pho_Detector", 1);
	oldTree->SetBranchStatus("pho_DetEta", 1);
	oldTree->SetBranchStatus("pho_DetPhi", 1);
	oldTree->SetBranchStatus("pho_XCes", 1);
	oldTree->SetBranchStatus("pho_ZCes", 1);
	oldTree->SetBranchStatus("pho_HadEm", 1);
	oldTree->SetBranchStatus("pho_Chi2Mean", 1);
	oldTree->SetBranchStatus("pho_N3d", 1);
	oldTree->SetBranchStatus("pho_Iso4", 1);
	oldTree->SetBranchStatus("pho_TrkPt", 1);		//from photon block
	oldTree->SetBranchStatus("pho_TrkIso", 1);		//from photon block
	oldTree->SetBranchStatus("pho_CesWireE2", 1);
	oldTree->SetBranchStatus("pho_CesStripE2", 1);
	oldTree->SetBranchStatus("pho_PhiWedge", 1);
	oldTree->SetBranchStatus("pho_NMuonStubs", 1);	//trkless muon stubs around 30 degree cone of the photon
	oldTree->SetBranchStatus("pho_EmTime", 1);
	oldTree->SetBranchStatus("pho_TightId", 1);
	oldTree->SetBranchStatus("pho_LooseId", 1);

	oldTree->SetBranchAddress("pho_num", &myStuple->pho_num);
	oldTree->SetBranchAddress("pho_Ntight", &myStuple->pho_Ntight);			// # of tight photons
	oldTree->SetBranchAddress("pho_Nloose", &myStuple->pho_Nloose);			// # of loose photons
	oldTree->SetBranchAddress("pho_Index", &myStuple->pho_Index);
	oldTree->SetBranchAddress("pho_Etc", &myStuple->pho_Etc);
	oldTree->SetBranchAddress("pho_E", &myStuple->pho_E);
	oldTree->SetBranchAddress("pho_Px", &myStuple->pho_Px);
	oldTree->SetBranchAddress("pho_Py", &myStuple->pho_Py);
	oldTree->SetBranchAddress("pho_Pz", &myStuple->pho_Pz);
	oldTree->SetBranchAddress("pho_Detector", &myStuple->pho_Detector);
	oldTree->SetBranchAddress("pho_DetEta", &myStuple->pho_DetEta);
	oldTree->SetBranchAddress("pho_DetPhi", &myStuple->pho_DetPhi);
	oldTree->SetBranchAddress("pho_XCes", &myStuple->pho_XCes);
	oldTree->SetBranchAddress("pho_ZCes", &myStuple->pho_ZCes);
	oldTree->SetBranchAddress("pho_HadEm", &myStuple->pho_HadEm);
	oldTree->SetBranchAddress("pho_Chi2Mean", &myStuple->pho_Chi2Mean);
	oldTree->SetBranchAddress("pho_N3d", &myStuple->pho_N3d);
	oldTree->SetBranchAddress("pho_Iso4", &myStuple->pho_Iso4);
	oldTree->SetBranchAddress("pho_TrkPt", &myStuple->pho_TrkPt);		//from photon block
	oldTree->SetBranchAddress("pho_TrkIso", &myStuple->pho_TrkIso);		//from photon block
	oldTree->SetBranchAddress("pho_CesWireE2", &myStuple->pho_CesWireE2);
	oldTree->SetBranchAddress("pho_CesStripE2", &myStuple->pho_CesStripE2);
	oldTree->SetBranchAddress("pho_PhiWedge", &myStuple->pho_PhiWedge);
	oldTree->SetBranchAddress("pho_NMuonStubs", &myStuple->pho_NMuonStubs);	//trkless muon stubs around 30 degree cone of the photon
	oldTree->SetBranchAddress("pho_EmTime", &myStuple->pho_EmTime);
	oldTree->SetBranchAddress("pho_TightId", &myStuple->pho_TightId);
	oldTree->SetBranchAddress("pho_LooseId", &myStuple->pho_LooseId);
	


	oldTree->SetBranchStatus("jet_num", 1);
	oldTree->SetBranchStatus("jet_NJet15", 1);
	oldTree->SetBranchStatus("jet_Index", 1);				//original index in jet block
	oldTree->SetBranchStatus("jet_Pt", 1);					//easy to get from JetFilter
	oldTree->SetBranchStatus("jet_E", 1);
	oldTree->SetBranchStatus("jet_Px", 1);
	oldTree->SetBranchStatus("jet_Py", 1);
	oldTree->SetBranchStatus("jet_Pz", 1);
	oldTree->SetBranchStatus("jet_DetEta", 1);
	oldTree->SetBranchStatus("jet_DetPhi", 1);
	oldTree->SetBranchStatus("jet_HadEm", 1);
	//		oldTree->SetBranchStatus("jet_Emfr", 1);
	oldTree->SetBranchStatus("jet_Ntowers", 1);
	oldTree->SetBranchStatus("jet_Ntracks", 1);
	oldTree->SetBranchStatus("jet_SeedIPhi", 1);
	oldTree->SetBranchStatus("jet_SeedIEta", 1);

	oldTree->SetBranchAddress("jet_num", &myStuple->jet_num);
	oldTree->SetBranchAddress("jet_NJet15", &myStuple->jet_NJet15);
	oldTree->SetBranchAddress("jet_Index", &myStuple->jet_Index);				//original index in jet block
	oldTree->SetBranchAddress("jet_Pt", &myStuple->jet_Pt);					//easy to get from JetFilter
	oldTree->SetBranchAddress("jet_E", &myStuple->jet_E);
	oldTree->SetBranchAddress("jet_Px", &myStuple->jet_Px);
	oldTree->SetBranchAddress("jet_Py", &myStuple->jet_Py);
	oldTree->SetBranchAddress("jet_Pz", &myStuple->jet_Pz);
	oldTree->SetBranchAddress("jet_DetEta", &myStuple->jet_DetEta);
	oldTree->SetBranchAddress("jet_DetPhi", &myStuple->jet_DetPhi);
	oldTree->SetBranchAddress("jet_HadEm", &myStuple->jet_HadEm);
	//		oldTree->SetBranchAddress("jet_Emfr", &myStuple->jet_Emfr);
	oldTree->SetBranchAddress("jet_Ntowers", &myStuple->jet_Ntowers);
	oldTree->SetBranchAddress("jet_Ntracks", &myStuple->jet_Ntracks);
	oldTree->SetBranchAddress("jet_SeedIPhi", &myStuple->jet_SeedIPhi);
	oldTree->SetBranchAddress("jet_SeedIEta", &myStuple->jet_SeedIEta);

	unsigned int entries = oldTree->GetEntries();
	if (nevts>0 && nevts<entries) entries = nevts;
	std::cout << "Entries to be processed = " << entries << std::endl;

	for (unsigned int i = 0; i < entries; ++i)
	{
		oldTree->GetEntry(i);

		if (debug)
		{
			int phonum = myStuple->pho_num;
			std::cout << myStuple->pho_num << std::endl;
			for (int j=0; j < phonum; j++)
			{
				std::cout << "PHOTON["<< j << "] HadEm = " << myStuple->pho_HadEm[j] << std::endl;
			}

			int jetnum = myStuple->jet_num;
			std::cout << myStuple->jet_num << std::endl;
			for (int j=0; j < jetnum; j++)
			{
				std::cout << "JET["<<j << "] HadEm = " << myStuple->jet_HadEm[j] << std::endl;
			}
		}

		newTree->Fill();
		if (i%1000 == 0)
		{
			int pct = ((i * 1.0)/(entries * 1.0)) * 100.0;
			std::cout << "Processed " << i  << " [" << pct << "%]" << std::endl;

		}
	}



	newTree->Write();
	std::cout << "***************** SUMMARY ****************** " << std::endl;
	std::cout << newTree->GetEntries() << " entries written to " << rootFile->GetName() << std::endl;
	rootFile->Close();

}	
