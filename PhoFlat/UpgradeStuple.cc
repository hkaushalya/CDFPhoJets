#include <iostream>
#include <sstream>
#include "PhoFlat/UpgradeStuple.hh"

/*
This will add new braches to v1 Stuples to make them v2 compatible. v2 has Generator level
info and MC flag. Otherwise there will be trouble Reading in uses the latest Stuple.cc/hh
available. This way I can easily upgrade old Stuples quickly.

In the TTree web page, there is a way to update the existing TTree by adding braches directly
to the same tree wihout creating a new file. But I prefer this as I can still keep the old
one if the job crashes or something weired happens during this process. You never know 
what happens with my bad luck!
So it is always better to be safe than sorry!. 03-13-2008

ALSO READ THE IMPORTANT NOTE AT THE END OF THIS FILE!

*/
//-------------------------------------------------------------------
UpgradeStuple::UpgradeStuple():
iProgressBy(500000),
sHistFileName("UpgradeStuple.root")
{
}

//-------------------------------------------------------------
void UpgradeStuple::Init(TChain *tree)
{
	iEventsSelected = 0;
	readIn = new ReadInAll_4StUpgrade(tree,&oldStuple);

	rootFile = new TFile(sHistFileName.c_str(),"RECREATE");
	if (rootFile->IsZombie()) {
		std::cout << __FILE__ << "::" << __LINE__ <<"::" << sHistFileName << "did not open.may be it already exists! remove it and try again!" << std::endl;
		exit(1);
	}

	//TTree::SetMaxTreeSize(1000000000);	//new file size ==1G just for fun!
	newTree = new TTree("Stuple","Sam's flat Stuple.");

		newTree->Branch("evt_McFlag"   , &newStuple.evt_McFlag   , "evt_McFlag/i");
		newTree->Branch("evt_RunNumber"   , &newStuple.evt_RunNumber   , "evt_RunNumber/i");
		newTree->Branch("evt_EventNumber" , &newStuple.evt_EventNumber , "ect_EventNumber/i");

		newTree->Branch("tri_pho25iso" , &newStuple.tri_pho25iso , "tri_pho25iso/I");
		newTree->Branch("tri_pho50"    , &newStuple.tri_pho50    , "tri_pho50/I");
		newTree->Branch("tri_pho70"    , &newStuple.tri_pho70    , "tri_pho70/I");
	
		newTree->Branch("pho_num"     , &newStuple.pho_num      , "pho_num/i");
		
		newTree->Branch("pho_Ntight"  , &newStuple.pho_Ntight   , "pho_Ntight/i");
		newTree->Branch("pho_Nloose"  , &newStuple.pho_Nloose   , "pho_Nloose/i");
		newTree->Branch("pho_Index"   , &newStuple.pho_Index    , "pho_Index[pho_num]/I");
		newTree->Branch("pho_Etc"     , &newStuple.pho_Etc      , "pho_Etc[pho_num]/F");
		newTree->Branch("pho_E"       , &newStuple.pho_E        , "pho_E[pho_num]/F");
		newTree->Branch("pho_Px"      , &newStuple.pho_Px       , "pho_Px[pho_num]/F");
		newTree->Branch("pho_Py"      , &newStuple.pho_Py       , "pho_Py[pho_num]/F");
		newTree->Branch("pho_Pz"      , &newStuple.pho_Pz       , "pho_Pz[pho_num]/F");
		newTree->Branch("pho_Detector", &newStuple.pho_Detector , "pho_Detector[pho_num]/I");
		newTree->Branch("pho_DetEta"  , &newStuple.pho_DetEta   , "pho_DetEta[pho_num]/F");
		newTree->Branch("pho_DetPhi"  , &newStuple.pho_DetPhi   , "pho_DetPhi[pho_num]/F");
		newTree->Branch("pho_XCes"    , &newStuple.pho_XCes     , "pho_XCes[pho_num]/F");
		newTree->Branch("pho_ZCes"    , &newStuple.pho_ZCes     , "pho_ZCes[pho_num]/F");
		newTree->Branch("pho_HadEm"   , &newStuple.pho_HadEm    , "pho_HadEm[pho_num]/F");
		newTree->Branch("pho_Chi2Mean", &newStuple.pho_Chi2Mean , "pho_Chi2Mean[pho_num]/F");
		newTree->Branch("pho_N3d"     , &newStuple.pho_N3d      , "pho_N3d[pho_num]/I");
		newTree->Branch("pho_Iso4"    , &newStuple.pho_Iso4     , "pho_Iso4[pho_num]/F");
		newTree->Branch("pho_TrkPt"   , &newStuple.pho_TrkPt    , "pho_TrkPt[pho_num]/F");
		newTree->Branch("pho_TrkIso"  , &newStuple.pho_TrkIso   , "pho_TrkIso[pho_num]/F");
		newTree->Branch("pho_CesWireE2"      , &newStuple.pho_CesWireE2      , "pho_CesWireE2[pho_num]/F");
		newTree->Branch("pho_CesStripE2"     , &newStuple.pho_CesStripE2     , "pho_CesStripE2[pho_num]/F");
		newTree->Branch("pho_PhiWedge"       , &newStuple.pho_PhiWedge       , "pho_PhiWedge[pho_num]/I");
		newTree->Branch("pho_NMuonStubs"     , &newStuple.pho_NMuonStubs     , "pho_NMuonStubs[pho_num]/I");
		newTree->Branch("pho_EmTime"         , &newStuple.pho_EmTime         , "pho_EmTime[pho_num]/F");
		newTree->Branch("pho_TightId"        , &newStuple.pho_TightId        , "pho_TightId[pho_num]/I");
		newTree->Branch("pho_LooseId"        , &newStuple.pho_LooseId        , "pho_LooseId[pho_num]/I");
		newTree->Branch("pho_PhoenixId"      , &newStuple.pho_PhoenixId      , "pho_PhoenixId[pho_num]/I");
		newTree->Branch("pho_Halo_seedWedge" , &newStuple.pho_Halo_seedWedge , "pho_Halo_seedWedge[pho_num]/I");
		newTree->Branch("pho_Halo_eastNhad"  , &newStuple.pho_Halo_eastNhad  , "pho_Halo_eastNhad[pho_num]/I");
		newTree->Branch("pho_Halo_westNhad"  , &newStuple.pho_Halo_westNhad  , "pho_Halo_westNhad[pho_num]/I");
		newTree->Branch("pho_matchJetIndex"  , &newStuple.pho_matchJetIndex  , "pho_matchJetIndex[pho_num]/I");
		

		newTree->Branch("ele_num"     , &newStuple.ele_num      , "ele_num/i");
		newTree->Branch("ele_Ntight"  , &newStuple.ele_Ntight   , "ele_Ntight/i");
		newTree->Branch("ele_Nloose"  , &newStuple.ele_Nloose   , "ele_Nloose/i");
		newTree->Branch("ele_Index"   , &newStuple.ele_Index    , "ele_Index[ele_num]/I");
		newTree->Branch("ele_Etc"     , &newStuple.ele_Etc      , "ele_Etc[ele_num]/F");
		newTree->Branch("ele_E"       , &newStuple.ele_E        , "ele_E[ele_num]/F");
		newTree->Branch("ele_Px"      , &newStuple.ele_Px       , "ele_Px[ele_num]/F");
		newTree->Branch("ele_Py"      , &newStuple.ele_Py       , "ele_Py[ele_num]/F");
		newTree->Branch("ele_Pz"      , &newStuple.ele_Pz       , "ele_Pz[ele_num]/F");
		newTree->Branch("ele_Detector", &newStuple.ele_Detector , "ele_Detector[ele_num]/I");
		newTree->Branch("ele_DetEta"  , &newStuple.ele_DetEta   , "ele_DetEta[ele_num]/F");
		newTree->Branch("ele_DetPhi"  , &newStuple.ele_DetPhi   , "ele_DetPhi[ele_num]/F");
		newTree->Branch("ele_XCes"    , &newStuple.ele_XCes     , "ele_XCes[ele_num]/F");
		newTree->Branch("ele_ZCes"    , &newStuple.ele_ZCes     , "ele_ZCes[ele_num]/F");
		newTree->Branch("ele_HadEm"   , &newStuple.ele_HadEm    , "ele_HadEm[ele_num]/F");
		newTree->Branch("ele_Chi2Mean", &newStuple.ele_Chi2Mean , "ele_Chi2Mean[ele_num]/F");
		newTree->Branch("ele_N3d"     , &newStuple.ele_N3d      , "ele_N3d[ele_num]/I");
		newTree->Branch("ele_Iso4"    , &newStuple.ele_Iso4     , "ele_Iso4[ele_num]/F");
		newTree->Branch("ele_TrkIso"  , &newStuple.ele_TrkIso   , "ele_TrkIso[ele_num]/F");
		newTree->Branch("ele_CesWireE2"      , &newStuple.ele_CesWireE2      , "ele_CesWireE2[ele_num]/F");
		newTree->Branch("ele_CesStripE2"     , &newStuple.ele_CesStripE2     , "ele_CesStripE2[ele_num]/F");
		newTree->Branch("ele_PhiWedge"       , &newStuple.ele_PhiWedge       , "ele_PhiWedge[ele_num]/I");
		newTree->Branch("ele_NMuonStubs"     , &newStuple.ele_NMuonStubs     , "ele_NMuonStubs[ele_num]/I");
		newTree->Branch("ele_EmTime"         , &newStuple.ele_EmTime         , "ele_EmTime[ele_num]/F");
		newTree->Branch("ele_PhoenixId"      , &newStuple.ele_PhoenixId      , "ele_PhoenixId[ele_num]/I");
		newTree->Branch("ele_Halo_seedWedge" , &newStuple.ele_Halo_seedWedge , "ele_Halo_seedWedge[ele_num]/I");
		newTree->Branch("ele_Halo_eastNhad"  , &newStuple.ele_Halo_eastNhad  , "ele_Halo_eastNhad[ele_num]/I");
		newTree->Branch("ele_Halo_westNhad"  , &newStuple.ele_Halo_westNhad  , "ele_Halo_westNhad[ele_num]/I");
		newTree->Branch("ele_matchJetIndex"  , &newStuple.ele_matchJetIndex  , "ele_matchJetIndex[ele_num]/I");
	
		newTree->Branch("ele_Ntracks"      , &newStuple.ele_Ntracks      , "ele_Ntracks[ele_num]/I");
		newTree->Branch("ele_Emfr"         , &newStuple.ele_Emfr         , "ele_Emfr[ele_num]/F");
		newTree->Branch("ele_EoverP"       , &newStuple.ele_EoverP       , "ele_EoverP[ele_num]/F");
		newTree->Branch("ele_TrackPt"      , &newStuple.ele_TrackPt      , "ele_TrackPt[ele_num]/F");
		newTree->Branch("ele_TrackBcPt"    , &newStuple.ele_TrackBcPt    , "ele_TrackBcPt[ele_num]/F");
		newTree->Branch("ele_TrackPhi"     , &newStuple.ele_TrackPhi     , "ele_TrackPhi[ele_num]/F");
		newTree->Branch("ele_Nssl"         , &newStuple.ele_Nssl         , "ele_Nssl[ele_num]/I");
		newTree->Branch("ele_Nasl"         , &newStuple.ele_Nasl         , "ele_Nasl[ele_num]/I");
		newTree->Branch("ele_TightId"      , &newStuple.ele_TightId      , "ele_TightId[ele_num]/I");
		newTree->Branch("ele_LooseId"      , &newStuple.ele_LooseId      , "ele_LooseId[ele_num]/I");
		newTree->Branch("ele_ConversionId" , &newStuple.ele_ConversionId , "ele_ConversionId[ele_num]/I");
		newTree->Branch("ele_matchJetIndex", &newStuple.ele_matchJetIndex, "ele_matchJetIndex[ele_num]/I");
		
		newTree->Branch("jet_num"      , &newStuple.jet_num      , "jet_num/i");
		newTree->Branch("jet_NJet15"   , &newStuple.jet_NJet15   , "jet_Njet15/I");
		newTree->Branch("jet_Index"    , &newStuple.jet_Index    , "jet_Index[jet_num]/I");
		newTree->Branch("jet_Pt"       , &newStuple.jet_Pt       , "jet_Pt[jet_num]/F");
		newTree->Branch("jet_E"        , &newStuple.jet_E        , "jet_E[jet_num]/F");
		newTree->Branch("jet_Px"       , &newStuple.jet_Px       , "jet_Px[jet_num]/F");
		newTree->Branch("jet_Py"       , &newStuple.jet_Py       , "jet_Py[jet_num]/F");
		newTree->Branch("jet_Pz"       , &newStuple.jet_Pz       , "jet_Pz[jet_num]/F");
		newTree->Branch("jet_DetEta"   , &newStuple.jet_DetEta   , "jet_DetEta[jet_num]/F");
		newTree->Branch("jet_DetPhi"   , &newStuple.jet_DetPhi   , "jet_DetPhi[jet_num]/F");
		newTree->Branch("jet_HadEm"    , &newStuple.jet_HadEm    , "jet_HadEm[jet_num]/F");
		newTree->Branch("jet_Emfr"     , &newStuple.jet_Emfr     , "jet_Emfr[jet_num]/F");
		newTree->Branch("jet_Ntowers"  , &newStuple.jet_Ntowers  , "jet_Ntowers[jet_num]/I");
		newTree->Branch("jet_Ntracks"  , &newStuple.jet_Ntracks  , "jet_Ntracks[jet_num]/I");
		newTree->Branch("jet_SeedIPhi" , &newStuple.jet_SeedIPhi , "jet_SeedIPhi[jet_num]/I");
		newTree->Branch("jet_SeedIEta" , &newStuple.jet_SeedIEta , "jet_SeedIEta[jet_num]/I");
		
		newTree->Branch("vtx_N"        , &newStuple.vtx_N        , "vtx_N/I");
		newTree->Branch("vtx_NClass12" , &newStuple.vtx_NClass12 , "vtx_NClass12/I");
		newTree->Branch("vtx_z"        , &newStuple.vtx_z        , "vtx_z/F");
		newTree->Branch("vtx_Ntracks"  , &newStuple.vtx_Ntracks  , "vtx_Ntracks/I");
		newTree->Branch("vtx_SumPt"    , &newStuple.vtx_SumPt    , "vtx_SumPt/F");

		newTree->Branch("met_Met"    , &newStuple.met_Met    , "met_Met/F");
		newTree->Branch("met_SumEt"  , &newStuple.met_SumEt  , "met_SumEt/F");
		newTree->Branch("met_Ht"     , &newStuple.met_Ht     , "met_Ht/F");
		newTree->Branch("met_MetPhi" , &newStuple.met_MetPhi , "met_MetPhi/F");

		newTree->Branch("met_Gen_d" , &newStuple.met_Gen_d , "met_Gen_d/F");
		newTree->Branch("met_Gen_m" , &newStuple.met_Gen_m , "met_Gen_m/F");
		newTree->Branch("met_Gen_p" , &newStuple.met_Gen_p , "met_Gen_p/F");
		newTree->Branch("met_Gen_mUn" , &newStuple.met_Gen_mUn , "met_Gen_mUn/F");
		newTree->Branch("met_Gen_pUn" , &newStuple.met_Gen_pUn , "met_Gen_pUn/F");

		

		newTree->Branch("gen_elenum" , &newStuple.gen_elenum , "gen_elenum/i");
		newTree->Branch("gen_phonum" , &newStuple.gen_phonum , "gen_phonum/i");
		
		newTree->Branch("gen_MomIndex" , &newStuple.gen_MomIndex , "gen_MomIndex/I");
		newTree->Branch("gen_MomPDG" , &newStuple.gen_MomPDG , "gen_MomPDG/I");
		newTree->Branch("gen_MomStaus" , &newStuple.gen_MomStatus , "gen_MomStatus/I");
		newTree->Branch("gen_MomEtc" , &newStuple.gen_MomEtc , "gen_MomEtc/F");
		newTree->Branch("gen_MomE" , &newStuple.gen_MomE , "gen_MomE/F");
		newTree->Branch("gen_MomPx" , &newStuple.gen_MomPx , "gen_MomPx/F");
		newTree->Branch("gen_MomPy" , &newStuple.gen_MomPy , "gen_MomPy/F");
		newTree->Branch("gen_MomPz" , &newStuple.gen_MomPz , "gen_MomPz/F");

		newTree->Branch("gen_ProdVtxX" , &newStuple.gen_ProdVtxX , "gen_ProdVtxX/F");
		newTree->Branch("gen_ProdVtxY" , &newStuple.gen_ProdVtxY , "gen_ProdVtxY/F");
		newTree->Branch("gen_ProdVtxZ" , &newStuple.gen_ProdVtxZ , "gen_ProdVtxZ/F");
		newTree->Branch("gen_ProdVtxT" , &newStuple.gen_ProdVtxT , "gen_ProdVtxT/F");

		newTree->Branch("gen_pho_Index" , &newStuple.gen_pho_Index , "gen_pho_Index[gen_phonum]/I");
		newTree->Branch("gen_pho_PDG" , &newStuple.gen_pho_PDG , "gen_pho_PDG[gen_phonum]/I");
		newTree->Branch("gen_pho_Status" , &newStuple.gen_pho_Status , "gen_pho_Status[gen_phonum]/I");
		newTree->Branch("gen_pho_Etc" , &newStuple.gen_pho_Etc , "gen_pho_Etc[gen_phonum]/F");
		newTree->Branch("gen_pho_E" , &newStuple.gen_pho_E , "gen_pho_E[gen_phonum]/F");
		newTree->Branch("gen_pho_Px" , &newStuple.gen_pho_Px , "gen_pho_Px[gen_phonum]/F");
		newTree->Branch("gen_pho_Py" , &newStuple.gen_pho_Py , "gen_pho_Py[gen_phonum]/F");
		newTree->Branch("gen_pho_Pz" , &newStuple.gen_pho_Pz , "gen_pho_Pz[gen_phonum]/F");

		newTree->Branch("gen_ele_Index" , &newStuple.gen_ele_Index , "gen_ele_Index[gen_elenum]/I");
		newTree->Branch("gen_ele_PDG" , &newStuple.gen_ele_PDG , "gen_ele_PDG[gen_elenum]/I");
		newTree->Branch("gen_ele_Status" , &newStuple.gen_ele_Status , "gen_ele_Status[gen_elenum]/I");
		newTree->Branch("gen_ele_Etc" , &newStuple.gen_ele_Etc , "gen_ele_Etc[gen_elenum]/F");
		newTree->Branch("gen_ele_E" , &newStuple.gen_ele_E , "gen_ele_E[gen_elenum]/F");
		newTree->Branch("gen_ele_Px" , &newStuple.gen_ele_Px , "gen_ele_Px[gen_elenum]/F");
		newTree->Branch("gen_ele_Py" , &newStuple.gen_ele_Py , "gen_ele_Py[gen_elenum]/F");
		newTree->Branch("gen_ele_Pz" , &newStuple.gen_ele_Pz , "gen_ele_Pz[gen_elenum]/F");
	

} // Init

//-------------------------------------------------------------------
void UpgradeStuple::CleanUp()
{
} 


//------ main -------------------------------------------------------
void UpgradeStuple::Main(TChain *ch, int iRunEvents)
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
	
	unsigned int iENTRIES = (unsigned int) myChain->GetEntries();
	unsigned int iEvt2Process = 0; 	 //_______________ events to be processed
	if (iRunEvents < 0 ) {
		iRunEvents = iENTRIES;	//____ requested all entries
		iEvt2Process = iENTRIES;
	} else if (iRunEvents > 0 && iRunEvents <= iENTRIES) iEvt2Process = iRunEvents;
	else iEvt2Process = iENTRIES;		//_____ default

	std::cout << " ==> Entries found :: " << iENTRIES << std::endl;


	unsigned int iEvtProc = 0;		// number of events processed

	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) {

	  	readIn->GetEntry(iEvtProc);

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << iEvtProc << "\t events processed." << std::endl;
		}
		newStuple.Init();
		newStuple = oldStuple;
		newTree->Fill();
	};
	
	// this is the way to save the tree. if you just call TFile->Write() it carshes when switching to to new file
	// when it reached max file size (1.9G);
	//see TTree web instructions and I have pasted part of it at the end of this file.
	rootFile = newTree->GetCurrentFile();
	rootFile->Write();
	rootFile->Close();

	std::cout << "======== UpgradeStuple =======================" << std::endl;
	std::cout << "Entries found    = " << iENTRIES << std::endl;
	std::cout << "Events Processed = " << iEvtProc << std::endl;
	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	std::cout << "======= END UpgradeStuple ====================" << std::endl;
} // Main

/*  

IMPORTANT NOTE:
Be careful when writing the final Tree header to the file!
Don't do:
TFile *file = new TFile("myfile.root","recreate");
TTree *T = new TTree("T","title");
T->Fill(); //loop
file->Write();
file->Close();
but do the following:
TFile *file = new TFile("myfile.root","recreate");
TTree *T = new TTree("T","title");
T->Fill(); //loop
file = T->GetCurrentFile(); //to get the pointer to the current file
file->Write();
file->Close();

*/
