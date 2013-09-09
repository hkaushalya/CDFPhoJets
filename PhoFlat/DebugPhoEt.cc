#include <iostream>
#include "PhoFlat/DebugPhoEt.hh"
#include "TFile.h"


//-------------------------------------------------------------------
DebugPhoEt::DebugPhoEt(std::string stupleFile)
{

	if (stupleFile.length() <1) {
		std::cout << __FILE__ << "::" << __LINE__ <<":: no file specified.!" << std::endl;
		exit (1);
	}
	
	TFile *file = new TFile(stupleFile.c_str());
	if (file->IsZombie()) {
		std::cout << __FILE__ << "::" << __LINE__ <<":: File not found.!" << std::endl;
		exit (1);
	}
	
	
	myTree = (TChain*) file->Get("Stuple");
		std::cout << "tree not found" << std::endl;
		return;
	}
	myStuple = st;

	
	myTree->SetBranchStatus("evt_RunNumber", 1);
	myTree->SetBranchStatus("evt_EventNumber", 1);
	myTree->SetBranchStatus("tri_pho25iso", 1);
	myTree->SetBranchStatus("tri_pho50", 1);
	myTree->SetBranchStatus("tri_pho70", 1);

	myTree->SetBranchStatus("pho_num", 1);		// number of electrons so we know how much to read from arrays
	myTree->SetBranchStatus("ele_num", 1);
	myTree->SetBranchStatus("jet_num", 1);
	myTree->SetBranchStatus("pho_Ntight", 1);			// # of tight photons
	myTree->SetBranchStatus("pho_Nloose", 1);			// # of loose photons
	myTree->SetBranchStatus("pho_Index", 1);
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

	myTree->SetBranchStatus("ele_Index", 1);
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

	myTree->SetBranchStatus("vtx_N", 1);
	myTree->SetBranchStatus("vtx_NClass12", 1);
	myTree->SetBranchStatus("vtx_z", 1);				// highest sum pt vertex Z
	myTree->SetBranchStatus("vtx_Ntracks", 1);		// in highest sum pt vertex
	myTree->SetBranchStatus("vtx_SumPt", 1);			// highest

	myTree->SetBranchStatus("met_Met", 1);				// corrected MET
	myTree->SetBranchStatus("met_SumEt", 1);			// corrected
	myTree->SetBranchStatus("met_Ht", 1);				// corrected
	myTree->SetBranchStatus("met_MetPhi", 1);			// corrected
	myTree->SetBranchStatus("met_Gen_d", 1);
	myTree->SetBranchStatus("met_Gen_m", 1);
	myTree->SetBranchStatus("met_Gen_p", 1);
	myTree->SetBranchStatus("met_Gen_mUn", 1);
	myTree->SetBranchStatus("met_Gen_pUn", 1);
	

	myTree->SetBranchAddress("evt_RunNumber", &myStuple->evt_RunNumber);
	myTree->SetBranchAddress("evt_EventNumber", &myStuple->evt_EventNumber);
	myTree->SetBranchAddress("tri_pho25iso", &myStuple->tri_pho25iso);
	myTree->SetBranchAddress("tri_pho50", &myStuple->tri_pho50);
	myTree->SetBranchAddress("tri_pho70", &myStuple->tri_pho70);

	myTree->SetBranchAddress("pho_num", &myStuple->pho_num);		// number of electrons so we know how much to read from arrays
	myTree->SetBranchAddress("ele_num", &myStuple->ele_num);
	myTree->SetBranchAddress("jet_num", &myStuple->jet_num);
	myTree->SetBranchAddress("pho_Ntight", &myStuple->pho_Ntight);			// # of tight photons
	myTree->SetBranchAddress("pho_Nloose", &myStuple->pho_Nloose);			// # of loose photons
	myTree->SetBranchAddress("pho_Index", &myStuple->pho_Index);
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


	myTree->SetBranchAddress("ele_Ntight", &myStuple->ele_Ntight);
	myTree->SetBranchAddress("ele_Nloose", &myStuple->ele_Nloose);
	myTree->SetBranchAddress("ele_Ntracks", &myStuple->ele_Ntracks);				//these are from TStnElectron
	myTree->SetBranchAddress("ele_Index", &myStuple->ele_Index);
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

	myTree->SetBranchAddress("vtx_N", &myStuple->vtx_N);
	myTree->SetBranchAddress("vtx_NClass12", &myStuple->vtx_NClass12);
	myTree->SetBranchAddress("vtx_z", &myStuple->vtx_z);				// highest sum pt vertex Z
	myTree->SetBranchAddress("vtx_Ntracks", &myStuple->vtx_Ntracks);		// in highest sum pt vertex
	myTree->SetBranchAddress("vtx_SumPt", &myStuple->vtx_SumPt);			// highest

	myTree->SetBranchAddress("met_Met", &myStuple->met_Met);				// corrected MET
	myTree->SetBranchAddress("met_SumEt", &myStuple->met_SumEt);			// corrected
	myTree->SetBranchAddress("met_Ht", &myStuple->met_Ht);				// corrected
	myTree->SetBranchAddress("met_MetPhi", &myStuple->met_MetPhi);			// corrected
	
	myTree->SetBranchAddress("met_Gen_d", &myStuple->met_Gen_d);
	myTree->SetBranchAddress("met_Gen_m", &myStuple->met_Gen_m);
	myTree->SetBranchAddress("met_Gen_p", &myStuple->met_Gen_p);
	myTree->SetBranchAddress("met_Gen_mUn", &myStuple->met_Gen_mUn);
	myTree->SetBranchAddress("met_Gen_pUn", &myStuple->met_Gen_pUn);
}

