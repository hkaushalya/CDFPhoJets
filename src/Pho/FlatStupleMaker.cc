////////////////////////////////////////////////////////////
// This will create flat stuples for me.                //
////////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/*{{{*/
/* $Id: FlatStupleMaker.cc,v 1.15 2011/05/24 00:46:46 samantha Exp $
 * $Log: FlatStupleMaker.cc,v $
 * Revision 1.15  2011/05/24 00:46:46  samantha
 * MODIFIED: To get the b-tag (SecVtx) jet info using TStnSecVtxTagBlock.
 * This required to call TUtil::FillJetTags() in rlc/Pho/TUtil.hh
 *
 * Revision 1.14  2011/04/13 16:31:18  samantha
 * ADDED: 2nd Track Pt poiting to photons or eletrons. For electrons, it seems little
 * harder to access track info witihin my framework. So I used the 2nd track info I
 * get for photons from photon block instead.
 *
 * Revision 1.13  2010/08/26 15:49:19  samantha
 * DROPPED: The min njet15 requirement as I need no nvtx/no jet events fo BH
 * template and rejection calcualtions. This also means I should NOT cut on nvtx.
 * ADDED: 1. Jet EM timing.
 *        2. EM obj threshold can be set.
 *        3. Raw jet info to the tree.
 * MODIFIED: 1. Replaced jetmod2 with jetmod3 which is stripdown version generating
 * only the jets/met etc I need.
 *         2. Halo info is saved only for data.
 *
 * Revision 1.12  2009/11/10 20:04:04  samantha
 * ADDED:  1. A message in case RunPermit fails.
 * 	2. Uncommented the TagForEmSysUp/Down
 * 	3. Auto CVS Log info.
 * MODIFIED: 1. Update the name changes in TagElectron and TagLooseElectron modules
 * function names, GetElectronIdWord() -> GetPhoLikeEleIdWord() and
 * GetElectronIdWord() -> GetPhoLikeEleIdWord()
 * 	2. Uncommented Ht accessor from jetfilter module.
 * 	3. Commented out the MEt model varibales.
 *
 *
 */
/*}}}*/

#include "samantha/Pho/FlatStupleMaker.hh"
#include "Stntuple/loop/TStnAna.hh"
#include <iostream>
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/obj/TStnElectron.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/Pho/TagLoosePhotons.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "samantha/Pho/TagLooseElectrons.hh"
#include "Stntuple/obj/TGenParticle.hh"
#include "TLorentzVector.h"
#include <iomanip>
#include <memory>
#include "TString.h"
#include "TObjString.h"
#include <assert.h>
#include "rlc/Pho/TUtil.hh" 

ClassImp(FlatStupleMaker)

//_____________________________________________________________________________
FlatStupleMaker::FlatStupleMaker(const char* name, const char* title):
  TStnModule(name,title),
  sFileName("Stuple.root"),
  bRunPermit(true),
  bNoSummary(false),
  fStupleVersion(1.0),
  fMinEtThr(15.)
{
	std::cout << "Hello I am FlatStupleMaker module" << std::endl;
}

//_____________________________________________________________________________
FlatStupleMaker::~FlatStupleMaker() {
}

//_____________________________________________________________________________
void FlatStupleMaker::SaveHistograms() {
}

//_____________________________________________________________________________
void FlatStupleMaker::BookHistograms()
{
	Hist.Npho = new TH1F("Npho","Photons per event",10,0,10);
	Hist.Nele = new TH1F("Nele","Electrons per event",10,0,10);
		
	TFolder* new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	new_folder->Add(Hist.Npho);
	new_folder->Add(Hist.Nele);

	
	Hist.Ncj = new TH1F("NJet","NJet>15GeV after removing all EM objs",20,0,20);
	Hist.Nupj = new TH1F("NJetup","NJet_{JES+}>15GeV after removing all EM objs",20,0,20);
	Hist.Ndownj = new TH1F("NJetdown","NJet_{JES-}>15GeV after removing all EM objs",20,0,20);
	Hist.NcjUpj = new TH2F("NJet_c_up","NJet_{JES+}>15GeV vs NJet>15GeV after removing all EM objs",200,0,20, 200, 0,20);
	Hist.NcjDownj = new TH2F("NJet_c_down","NJet_{JES-}>15GeV vs NJet>15GeV after removing all EM objs",200,0,20, 200, 0,20);
	Hist.cj_EperJet = new TH1F("EperJ_cen","#Sigma E_{T}^{Jets>15GeV} after removing all EM objs", 50, 0, 250);
	Hist.upj_EperJet = new TH1F("EperJ_up","#Sigma E_{T}^{Jets_{JES+}>15GeV} after removing all EM objs", 50, 0, 250);
	Hist.downj_EperJet = new TH1F("EperJ_down","#Sigma E_{T}^{Jets{JES-}>15GeV} after removing all EM objs", 50, 0, 250);


	new_folder->Add(Hist.Ncj);
	new_folder->Add(Hist.Nupj);
	new_folder->Add(Hist.Ndownj);
	new_folder->Add(Hist.NcjUpj);
	new_folder->Add(Hist.NcjDownj);
	new_folder->Add(Hist.cj_EperJet);
	new_folder->Add(Hist.upj_EperJet);
	new_folder->Add(Hist.downj_EperJet);

	Hist.w_met_c = new TH1F("wmet_c","std. loose cen. ele. selection : W MET",40,0,200);
	Hist.w_tmass_c = new TH1F("w_tmass_c","std. loose cen. ele. selection : W transverse mass",40,0,200);
	Hist.w_met_p = new TH1F("wmet_p","std. loose plug ele. selection : W MET",40,0,200);
	Hist.w_tmass_p = new TH1F("w_tmass_p","std. loose plug ele. selection : W transverse mass",40,0,200);
	Hist.z_met_cc = new TH1F("zmet_cc","std. loose cen. cen. ele. selection : Z MET",40,0,200);
	Hist.z_mass_cc = new TH1F("z_mass_cc","std. loose cen. cen. ele. selection : Z mass",40,0,200);
	Hist.z_met_cp = new TH1F("zmet_cp","std. loose cen. plug ele. selection : Z MET",40,0,200);
	Hist.z_mass_cp = new TH1F("z_mass_cp","std. loose cen. plug ele. selection : Z mass",40,0,200);
	Hist.z_met_pp = new TH1F("zmet_pp","std. loose plug plug ele. selection : Z MET",40,0,200);
	Hist.z_mass_pp = new TH1F("z_mass_pp","std. loose plug plug ele. selection : Z mass",40,0,200);
	Hist.ele_stdVsPholike = new TH1F("ele_stdVsPholike","std. loose plug ele. selection and photon-like electron selection",10,0,10);

	new_folder->Add(Hist.w_met_c);
	new_folder->Add(Hist.w_tmass_c);
	new_folder->Add(Hist.w_met_p);
	new_folder->Add(Hist.w_tmass_p);
	new_folder->Add(Hist.z_met_cc);
	new_folder->Add(Hist.z_mass_cc);
	new_folder->Add(Hist.z_met_cp);
	new_folder->Add(Hist.z_mass_cp);
	new_folder->Add(Hist.z_met_pp);
	new_folder->Add(Hist.z_mass_pp);
	new_folder->Add(Hist.ele_stdVsPholike);

}


//_____________________________________________________________________________
int FlatStupleMaker::BeginJob()
{

  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}

  	jetMod = (JetFilterModuleV3*) ((TStnAna*) GetAna()->GetModule("JetFilterV3"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		bRunPermit = false;
	}

  	jetModUp = (JetFilterModuleV3*) ((TStnAna*) GetAna()->GetModule("JetFilterV3JESup"));
	if (!jetModUp) {
		StdOut(__FILE__,__LINE__,3,"JetFilterJESup module required!.");
		//bRunPermit = false;
	}

  	jetModDown = (JetFilterModuleV3*) ((TStnAna*) GetAna()->GetModule("JetFilterV3JESdown"));
	if (!jetModDown) {
		StdOut(__FILE__,__LINE__,3,"JetFilterJESdown module required!.");
		//bRunPermit = false;
	}

  	trigMod = (TriggerModule*) ((TStnAna*) GetAna()->GetModule("TriggerModule"));
	if (!trigMod) {
		StdOut(__FILE__,__LINE__,3,"TriggerModule required!.");
		bRunPermit = false;
	}

	tagHaloMod = (TagBeamHalo*) ((TStnAna*) GetAna()->GetModule("TagBeamHalo"));
	if (!tagHaloMod) {
		StdOut(__FILE__,__LINE__,3,"TagBeamHalo module required!.");
		bRunPermit = false;
	}

		
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");
	fJetBlock 		= (TStnJetBlock*)      RegisterDataBlock("JetBlock","TStnJetBlock");
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");
	fTrackBlock 	= (TStnTrackBlock*)    RegisterDataBlock("TrackBlock","TStnTrackBlock");
	fSecVtxTagBlock = (TStnSecVtxTagBlock*)   RegisterDataBlock("SecVtxTagBlock","TStnSecVtxTagBlock");
		
	fCalBlock 		= (TCalDataBlock*)     RegisterDataBlock("CalDataBlock","TCalDataBlock");
//	if (GetHeaderBlock()->McFlag()) {
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");
//	}
	BookHistograms();

	counter.evtsProcessed	= 0;
	counter.evtsPassModule 	= 0;
	counter.evtsWritten 		= 0;
	counter.phoOnlyEvts		= 0;
	counter.eleOnlyEvts		= 0;
	counter.phoEleEvts		= 0;

	//to get b-tag info filled in TStnJet
	 TString mistag_file = "TightSECVTXparam_5600invpb.root"; 
	 matrix = new jetMistagJun2009(mistag_file,"tagrates");
	

	rootFile = new TFile(sFileName.c_str(),"RECREATE");
	std::stringstream treeTitle;
	treeTitle << "Sam's flat Stuple- version-" << fStupleVersion;
	tree = new TTree("Stuple", treeTitle.str().c_str());
//	TObjString str("sadaksjd asf; daksfj sd");
	//tree->GetUserInfo()->Add(&str);

		tree->Branch("evt_McFlag"   	 , &stuple.evt_McFlag   	, "evt_McFlag/I");
		tree->Branch("evt_RunNumber"   , &stuple.evt_RunNumber   , "evt_RunNumber/i");
		tree->Branch("evt_EventNumber" , &stuple.evt_EventNumber , "ect_EventNumber/i");

		tree->Branch("tri_pho25iso" , &stuple.tri_pho25iso , "tri_pho25iso/I");
		tree->Branch("tri_pho50"    , &stuple.tri_pho50    , "tri_pho50/I");
		tree->Branch("tri_pho70"    , &stuple.tri_pho70    , "tri_pho70/I");
	
		tree->Branch("pho_num"     , &stuple.pho_num      , "pho_num/i");
		tree->Branch("pho_Ntight"  , &stuple.pho_Ntight   , "pho_Ntight/i");
		tree->Branch("pho_Nloose"  , &stuple.pho_Nloose   , "pho_Nloose/i");
		tree->Branch("pho_Index"   , stuple.pho_Index    , "pho_Index[pho_num]/I");
		tree->Branch("pho_PhoBlockIndex"   , stuple.pho_PhoBlockIndex    , "pho_PhoBlockIndex[pho_num]/I");
		tree->Branch("pho_Etc"     , stuple.pho_Etc      , "pho_Etc[pho_num]/F");
		tree->Branch("pho_E"       , stuple.pho_E        , "pho_E[pho_num]/F");
		tree->Branch("pho_Px"      , stuple.pho_Px       , "pho_Px[pho_num]/F");
		tree->Branch("pho_Py"      , stuple.pho_Py       , "pho_Py[pho_num]/F");
		tree->Branch("pho_Pz"      , stuple.pho_Pz       , "pho_Pz[pho_num]/F");
		tree->Branch("pho_Detector", stuple.pho_Detector , "pho_Detector[pho_num]/I");
		tree->Branch("pho_DetEta"  , stuple.pho_DetEta   , "pho_DetEta[pho_num]/F");
		tree->Branch("pho_DetPhi"  , stuple.pho_DetPhi   , "pho_DetPhi[pho_num]/F");
		tree->Branch("pho_XCes"    , stuple.pho_XCes     , "pho_XCes[pho_num]/F");
		tree->Branch("pho_ZCes"    , stuple.pho_ZCes     , "pho_ZCes[pho_num]/F");
		tree->Branch("pho_HadEm"   , stuple.pho_HadEm    , "pho_HadEm[pho_num]/F");
		tree->Branch("pho_Chi2Mean", stuple.pho_Chi2Mean , "pho_Chi2Mean[pho_num]/F");
		tree->Branch("pho_N3d"     , stuple.pho_N3d      , "pho_N3d[pho_num]/I");
		tree->Branch("pho_Iso4"    , stuple.pho_Iso4     , "pho_Iso4[pho_num]/F");
		tree->Branch("pho_TrkPt"   , stuple.pho_TrkPt    , "pho_TrkPt[pho_num]/F");
		tree->Branch("pho_TrkPt2"  , stuple.pho_TrkPt2   , "pho_TrkPt2[pho_num]/F");
		tree->Branch("pho_TrkIso"  , stuple.pho_TrkIso   , "pho_TrkIso[pho_num]/F");
		tree->Branch("pho_CesWireE2"      , stuple.pho_CesWireE2      , "pho_CesWireE2[pho_num]/F");
		tree->Branch("pho_CesStripE2"     , stuple.pho_CesStripE2     , "pho_CesStripE2[pho_num]/F");
		tree->Branch("pho_PhiWedge"       , stuple.pho_PhiWedge       , "pho_PhiWedge[pho_num]/I");
		tree->Branch("pho_NMuonStubs"     , stuple.pho_NMuonStubs     , "pho_NMuonStubs[pho_num]/I");
		tree->Branch("pho_EmTime"         , stuple.pho_EmTime         , "pho_EmTime[pho_num]/F");
		tree->Branch("pho_TightId"        , stuple.pho_TightId        , "pho_TightId[pho_num]/I");
		tree->Branch("pho_LooseId"        , stuple.pho_LooseId        , "pho_LooseId[pho_num]/I");
		tree->Branch("pho_PhoenixId"      , stuple.pho_PhoenixId      , "pho_PhoenixId[pho_num]/I");
		tree->Branch("pho_Halo_seedWedge" , stuple.pho_Halo_seedWedge , "pho_Halo_seedWedge[pho_num]/I");
		tree->Branch("pho_Halo_eastNhad"  , stuple.pho_Halo_eastNhad  , "pho_Halo_eastNhad[pho_num]/I");
		tree->Branch("pho_Halo_westNhad"  , stuple.pho_Halo_westNhad  , "pho_Halo_westNhad[pho_num]/I");
		tree->Branch("pho_matchJetIndex"  , stuple.pho_matchJetIndex  , "pho_matchJetIndex[pho_num]/I");
		tree->Branch("pho_CprWgt"  , stuple.pho_CprWgt  , "pho_CprWgt[pho_num]/F");
		tree->Branch("pho_CprSys1"  , stuple.pho_CprSys1  , "pho_CprSys1[pho_num]/F");
		tree->Branch("pho_CprSys2"  , stuple.pho_CprSys2  , "pho_CprSys2[pho_num]/F");
		tree->Branch("pho_CprSys3"  , stuple.pho_CprSys3  , "pho_CprSys3[pho_num]/F");
		tree->Branch("pho_CprSys4"  , stuple.pho_CprSys4  , "pho_CprSys4[pho_num]/F");
		tree->Branch("pho_CprSys5"  , stuple.pho_CprSys5  , "pho_CprSys5[pho_num]/F");
		tree->Branch("pho_CprSys6"  , stuple.pho_CprSys6  , "pho_CprSys6[pho_num]/F");
		tree->Branch("pho_CprSys7"  , stuple.pho_CprSys7  , "pho_CprSys7[pho_num]/F");
		tree->Branch("pho_CprSys8"  , stuple.pho_CprSys8  , "pho_CprSys8[pho_num]/F");
		

			// Pho EM E uncertainty up 1%
		tree->Branch("pho_up_num"     , &stuple.pho_up_num      , "pho_up_num/i");
		tree->Branch("pho_up_Ntight"  , &stuple.pho_up_Ntight   , "pho_up_Ntight/i");
		tree->Branch("pho_up_Nloose"  , &stuple.pho_up_Nloose   , "pho_up_Nloose/i");
		tree->Branch("pho_up_Index"   , stuple.pho_up_Index    , "pho_up_Index[pho_up_num]/I");
		tree->Branch("pho_up_PhoBlockIndex"   , stuple.pho_up_PhoBlockIndex    , "pho_up_PhoBlockIndex[pho_up_num]/I");
		tree->Branch("pho_up_Etc"     , stuple.pho_up_Etc      , "pho_up_Etc[pho_up_num]/F");
		tree->Branch("pho_up_E"       , stuple.pho_up_E        , "pho_up_E[pho_up_num]/F");
		tree->Branch("pho_up_Px"      , stuple.pho_up_Px       , "pho_up_Px[pho_up_num]/F");
		tree->Branch("pho_up_Py"      , stuple.pho_up_Py       , "pho_up_Py[pho_up_num]/F");
		tree->Branch("pho_up_Pz"      , stuple.pho_up_Pz       , "pho_up_Pz[pho_up_num]/F");
		tree->Branch("pho_up_Detector", stuple.pho_up_Detector , "pho_up_Detector[pho_up_num]/I");
		tree->Branch("pho_up_DetEta"  , stuple.pho_up_DetEta   , "pho_up_DetEta[pho_up_num]/F");
		tree->Branch("pho_up_DetPhi"  , stuple.pho_up_DetPhi   , "pho_up_DetPhi[pho_up_num]/F");
		tree->Branch("pho_up_XCes"    , stuple.pho_up_XCes     , "pho_up_XCes[pho_up_num]/F");
		tree->Branch("pho_up_ZCes"    , stuple.pho_up_ZCes     , "pho_up_ZCes[pho_up_num]/F");
		tree->Branch("pho_up_HadEm"   , stuple.pho_up_HadEm    , "pho_up_HadEm[pho_up_num]/F");
		tree->Branch("pho_up_Chi2Mean", stuple.pho_up_Chi2Mean , "pho_up_Chi2Mean[pho_up_num]/F");
		tree->Branch("pho_up_N3d"     , stuple.pho_up_N3d      , "pho_up_N3d[pho_up_num]/I");
		tree->Branch("pho_up_Iso4"    , stuple.pho_up_Iso4     , "pho_up_Iso4[pho_up_num]/F");
		tree->Branch("pho_up_TrkPt"   , stuple.pho_up_TrkPt    , "pho_up_TrkPt[pho_up_num]/F");
		tree->Branch("pho_up_TrkPt2"  , stuple.pho_up_TrkPt2   , "pho_up_TrkPt2[pho_up_num]/F");
		tree->Branch("pho_up_TrkIso"  , stuple.pho_up_TrkIso   , "pho_up_TrkIso[pho_up_num]/F");
		tree->Branch("pho_up_CesWireE2"      , stuple.pho_up_CesWireE2      , "pho_up_CesWireE2[pho_up_num]/F");
		tree->Branch("pho_up_CesStripE2"     , stuple.pho_up_CesStripE2     , "pho_up_CesStripE2[pho_up_num]/F");
		tree->Branch("pho_up_PhiWedge"       , stuple.pho_up_PhiWedge       , "pho_up_PhiWedge[pho_up_num]/I");
		tree->Branch("pho_up_NMuonStubs"     , stuple.pho_up_NMuonStubs     , "pho_up_NMuonStubs[pho_up_num]/I");
		tree->Branch("pho_up_EmTime"         , stuple.pho_up_EmTime         , "pho_up_EmTime[pho_up_num]/F");
		tree->Branch("pho_up_TightId"        , stuple.pho_up_TightId        , "pho_up_TightId[pho_up_num]/I");
		tree->Branch("pho_up_LooseId"        , stuple.pho_up_LooseId        , "pho_up_LooseId[pho_up_num]/I");
		tree->Branch("pho_up_PhoenixId"      , stuple.pho_up_PhoenixId      , "pho_up_PhoenixId[pho_up_num]/I");
		tree->Branch("pho_up_Halo_seedWedge" , stuple.pho_up_Halo_seedWedge , "pho_up_Halo_seedWedge[pho_up_num]/I");
		tree->Branch("pho_up_Halo_eastNhad"  , stuple.pho_up_Halo_eastNhad  , "pho_up_Halo_eastNhad[pho_up_num]/I");
		tree->Branch("pho_up_Halo_westNhad"  , stuple.pho_up_Halo_westNhad  , "pho_up_Halo_westNhad[pho_up_num]/I");
		tree->Branch("pho_up_matchJetIndex"  , stuple.pho_up_matchJetIndex  , "pho_up_matchJetIndex[pho_up_num]/I");
		tree->Branch("pho_up_CprWgt"  , stuple.pho_up_CprWgt  , "pho_up_CprWgt[pho_num]/F");
		tree->Branch("pho_up_CprSys1"  , stuple.pho_up_CprSys1  , "pho_up_CprSys1[pho_num]/F");
		tree->Branch("pho_up_CprSys2"  , stuple.pho_up_CprSys2  , "pho_up_CprSys2[pho_num]/F");
		tree->Branch("pho_up_CprSys3"  , stuple.pho_up_CprSys3  , "pho_up_CprSys3[pho_num]/F");
		tree->Branch("pho_up_CprSys4"  , stuple.pho_up_CprSys4  , "pho_up_CprSys4[pho_num]/F");
		tree->Branch("pho_up_CprSys5"  , stuple.pho_up_CprSys5  , "pho_up_CprSys5[pho_num]/F");
		tree->Branch("pho_up_CprSys6"  , stuple.pho_up_CprSys6  , "pho_up_CprSys6[pho_num]/F");
		tree->Branch("pho_up_CprSys7"  , stuple.pho_up_CprSys7  , "pho_up_CprSys7[pho_num]/F");
		tree->Branch("pho_up_CprSys8"  , stuple.pho_up_CprSys8  , "pho_up_CprSys8[pho_num]/F");
	
			
			// Pho EM E uncertainty down 1%
		tree->Branch("pho_down_num"     , &stuple.pho_down_num      , "pho_down_num/i");
		tree->Branch("pho_down_Ntight"  , &stuple.pho_down_Ntight   , "pho_down_Ntight/i");
		tree->Branch("pho_down_Nloose"  , &stuple.pho_down_Nloose   , "pho_down_Nloose/i");
		tree->Branch("pho_down_Index"   , stuple.pho_down_Index    , "pho_down_Index[pho_down_num]/I");
		tree->Branch("pho_down_PhoBlockIndex"   , stuple.pho_down_PhoBlockIndex    , "pho_down_PhoBlockIndex[pho_down_num]/I");
		tree->Branch("pho_down_Etc"     , stuple.pho_down_Etc      , "pho_down_Etc[pho_down_num]/F");
		tree->Branch("pho_down_E"       , stuple.pho_down_E        , "pho_down_E[pho_down_num]/F");
		tree->Branch("pho_down_Px"      , stuple.pho_down_Px       , "pho_down_Px[pho_down_num]/F");
		tree->Branch("pho_down_Py"      , stuple.pho_down_Py       , "pho_down_Py[pho_down_num]/F");
		tree->Branch("pho_down_Pz"      , stuple.pho_down_Pz       , "pho_down_Pz[pho_down_num]/F");
		tree->Branch("pho_down_Detector", stuple.pho_down_Detector , "pho_down_Detector[pho_down_num]/I");
		tree->Branch("pho_down_DetEta"  , stuple.pho_down_DetEta   , "pho_down_DetEta[pho_down_num]/F");
		tree->Branch("pho_down_DetPhi"  , stuple.pho_down_DetPhi   , "pho_down_DetPhi[pho_down_num]/F");
		tree->Branch("pho_down_XCes"    , stuple.pho_down_XCes     , "pho_down_XCes[pho_down_num]/F");
		tree->Branch("pho_down_ZCes"    , stuple.pho_down_ZCes     , "pho_down_ZCes[pho_down_num]/F");
		tree->Branch("pho_down_HadEm"   , stuple.pho_down_HadEm    , "pho_down_HadEm[pho_down_num]/F");
		tree->Branch("pho_down_Chi2Mean", stuple.pho_down_Chi2Mean , "pho_down_Chi2Mean[pho_down_num]/F");
		tree->Branch("pho_down_N3d"     , stuple.pho_down_N3d      , "pho_down_N3d[pho_down_num]/I");
		tree->Branch("pho_down_Iso4"    , stuple.pho_down_Iso4     , "pho_down_Iso4[pho_down_num]/F");
		tree->Branch("pho_down_TrkPt"   , stuple.pho_down_TrkPt    , "pho_down_TrkPt[pho_down_num]/F");
		tree->Branch("pho_down_TrkPt2"  , stuple.pho_down_TrkPt2   , "pho_down_TrkPt2[pho_down_num]/F");
		tree->Branch("pho_down_TrkIso"  , stuple.pho_down_TrkIso   , "pho_down_TrkIso[pho_down_num]/F");
		tree->Branch("pho_down_CesWireE2"      , stuple.pho_down_CesWireE2      , "pho_down_CesWireE2[pho_down_num]/F");
		tree->Branch("pho_down_CesStripE2"     , stuple.pho_down_CesStripE2     , "pho_down_CesStripE2[pho_down_num]/F");
		tree->Branch("pho_down_PhiWedge"       , stuple.pho_down_PhiWedge       , "pho_down_PhiWedge[pho_down_num]/I");
		tree->Branch("pho_down_NMuonStubs"     , stuple.pho_down_NMuonStubs     , "pho_down_NMuonStubs[pho_down_num]/I");
		tree->Branch("pho_down_EmTime"         , stuple.pho_down_EmTime         , "pho_down_EmTime[pho_down_num]/F");
		tree->Branch("pho_down_TightId"        , stuple.pho_down_TightId        , "pho_down_TightId[pho_down_num]/I");
		tree->Branch("pho_down_LooseId"        , stuple.pho_down_LooseId        , "pho_down_LooseId[pho_down_num]/I");
		tree->Branch("pho_down_PhoenixId"      , stuple.pho_down_PhoenixId      , "pho_down_PhoenixId[pho_down_num]/I");
		tree->Branch("pho_down_Halo_seedWedge" , stuple.pho_down_Halo_seedWedge , "pho_down_Halo_seedWedge[pho_down_num]/I");
		tree->Branch("pho_down_Halo_eastNhad"  , stuple.pho_down_Halo_eastNhad  , "pho_down_Halo_eastNhad[pho_down_num]/I");
		tree->Branch("pho_down_Halo_westNhad"  , stuple.pho_down_Halo_westNhad  , "pho_down_Halo_westNhad[pho_down_num]/I");
		tree->Branch("pho_down_matchJetIndex"  , stuple.pho_down_matchJetIndex  , "pho_down_matchJetIndex[pho_down_num]/I");
		tree->Branch("pho_down_CprWgt"  , stuple.pho_down_CprWgt  , "pho_down_CprWgt[pho_num]/F");
		tree->Branch("pho_down_CprSys1"  , stuple.pho_down_CprSys1  , "pho_down_CprSys1[pho_num]/F");
		tree->Branch("pho_down_CprSys2"  , stuple.pho_down_CprSys2  , "pho_down_CprSys2[pho_num]/F");
		tree->Branch("pho_down_CprSys3"  , stuple.pho_down_CprSys3  , "pho_down_CprSys3[pho_num]/F");
		tree->Branch("pho_down_CprSys4"  , stuple.pho_down_CprSys4  , "pho_down_CprSys4[pho_num]/F");
		tree->Branch("pho_down_CprSys5"  , stuple.pho_down_CprSys5  , "pho_down_CprSys5[pho_num]/F");
		tree->Branch("pho_down_CprSys6"  , stuple.pho_down_CprSys6  , "pho_down_CprSys6[pho_num]/F");
		tree->Branch("pho_down_CprSys7"  , stuple.pho_down_CprSys7  , "pho_down_CprSys7[pho_num]/F");
		tree->Branch("pho_down_CprSys8"  , stuple.pho_down_CprSys8  , "pho_down_CprSys8[pho_num]/F");
	


		//electron collections	

		tree->Branch("ele_num"     , &stuple.ele_num      , "ele_num/i");
		tree->Branch("ele_Ntight"  , &stuple.ele_Ntight   , "ele_Ntight/i");
		tree->Branch("ele_Nloose"  , &stuple.ele_Nloose   , "ele_Nloose/i");
		tree->Branch("ele_Index"   , stuple.ele_Index    , "ele_Index[ele_num]/I");
		tree->Branch("ele_PhoBlockIndex", stuple.ele_PhoBlockIndex , "ele_PhoBlockIndex[ele_num]/I");
		tree->Branch("ele_EleBlockIndex", stuple.ele_EleBlockIndex, "ele_EleBlockIndex[ele_num]/I");
		tree->Branch("ele_Etc"     , stuple.ele_Etc      , "ele_Etc[ele_num]/F");
		tree->Branch("ele_E"       , stuple.ele_E        , "ele_E[ele_num]/F");
		tree->Branch("ele_Px"      , stuple.ele_Px       , "ele_Px[ele_num]/F");
		tree->Branch("ele_Py"      , stuple.ele_Py       , "ele_Py[ele_num]/F");
		tree->Branch("ele_Pz"      , stuple.ele_Pz       , "ele_Pz[ele_num]/F");
		tree->Branch("ele_Detector", stuple.ele_Detector , "ele_Detector[ele_num]/I");
		tree->Branch("ele_DetEta"  , stuple.ele_DetEta   , "ele_DetEta[ele_num]/F");
		tree->Branch("ele_DetPhi"  , stuple.ele_DetPhi   , "ele_DetPhi[ele_num]/F");
		tree->Branch("ele_XCes"    , stuple.ele_XCes     , "ele_XCes[ele_num]/F");
		tree->Branch("ele_ZCes"    , stuple.ele_ZCes     , "ele_ZCes[ele_num]/F");
		tree->Branch("ele_HadEm"   , stuple.ele_HadEm    , "ele_HadEm[ele_num]/F");
		tree->Branch("ele_Chi2Mean", stuple.ele_Chi2Mean , "ele_Chi2Mean[ele_num]/F");
		tree->Branch("ele_N3d"     , stuple.ele_N3d      , "ele_N3d[ele_num]/I");
		tree->Branch("ele_Iso4"    , stuple.ele_Iso4     , "ele_Iso4[ele_num]/F");
		tree->Branch("ele_TrkIso"  , stuple.ele_TrkIso   , "ele_TrkIso[ele_num]/F");
		tree->Branch("ele_CesWireE2"      , stuple.ele_CesWireE2      , "ele_CesWireE2[ele_num]/F");
		tree->Branch("ele_CesStripE2"     , stuple.ele_CesStripE2     , "ele_CesStripE2[ele_num]/F");
		tree->Branch("ele_PhiWedge"       , stuple.ele_PhiWedge       , "ele_PhiWedge[ele_num]/I");
		tree->Branch("ele_NMuonStubs"     , stuple.ele_NMuonStubs     , "ele_NMuonStubs[ele_num]/I");
		tree->Branch("ele_EmTime"         , stuple.ele_EmTime         , "ele_EmTime[ele_num]/F");
		tree->Branch("ele_PhoenixId"      , stuple.ele_PhoenixId      , "ele_PhoenixId[ele_num]/I");
		tree->Branch("ele_Halo_seedWedge" , stuple.ele_Halo_seedWedge , "ele_Halo_seedWedge[ele_num]/I");
		tree->Branch("ele_Halo_eastNhad"  , stuple.ele_Halo_eastNhad  , "ele_Halo_eastNhad[ele_num]/I");
		tree->Branch("ele_Halo_westNhad"  , stuple.ele_Halo_westNhad  , "ele_Halo_westNhad[ele_num]/I");
		tree->Branch("ele_matchJetIndex"  , stuple.ele_matchJetIndex  , "ele_matchJetIndex[ele_num]/I");
		tree->Branch("ele_Ntracks"      , stuple.ele_Ntracks      , "ele_Ntracks[ele_num]/I");
		tree->Branch("ele_Emfr"         , stuple.ele_Emfr         , "ele_Emfr[ele_num]/F");
		tree->Branch("ele_EoverP"       , stuple.ele_EoverP       , "ele_EoverP[ele_num]/F");
		tree->Branch("ele_TrackPt"      , stuple.ele_TrackPt      , "ele_TrackPt[ele_num]/F");
		tree->Branch("ele_TrackPt2"     , stuple.ele_TrackPt2     , "ele_TrackPt2[ele_num]/F");
		tree->Branch("ele_TrackBcPt"    , stuple.ele_TrackBcPt    , "ele_TrackBcPt[ele_num]/F");
		tree->Branch("ele_TrackPhi"     , stuple.ele_TrackPhi     , "ele_TrackPhi[ele_num]/F");
		tree->Branch("ele_Nssl"         , stuple.ele_Nssl         , "ele_Nssl[ele_num]/I");
		tree->Branch("ele_Nasl"         , stuple.ele_Nasl         , "ele_Nasl[ele_num]/I");
		tree->Branch("ele_TightId"      , stuple.ele_TightId      , "ele_TightId[ele_num]/I");
		tree->Branch("ele_LooseId"      , stuple.ele_LooseId      , "ele_LooseId[ele_num]/I");
		tree->Branch("ele_ConversionId" , stuple.ele_ConversionId , "ele_ConversionId[ele_num]/I");
		tree->Branch("ele_StdTightId"      , stuple.ele_StdTightId      , "ele_StdTightId[ele_num]/I");
		tree->Branch("ele_StdLooseId"      , stuple.ele_StdLooseId      , "ele_StdLooseId[ele_num]/I");


			// Ele EM E uncertainty up 1%
		tree->Branch("ele_up_num"     , &stuple.ele_up_num      , "ele_up_num/i");
		tree->Branch("ele_up_Ntight"  , &stuple.ele_up_Ntight   , "ele_up_Ntight/i");
		tree->Branch("ele_up_Nloose"  , &stuple.ele_up_Nloose   , "ele_up_Nloose/i");
		tree->Branch("ele_up_Index"   , stuple.ele_up_Index    , "ele_up_Index[ele_up_num]/I");
		tree->Branch("ele_up_PhoBlockIndex", stuple.ele_up_PhoBlockIndex , "ele_up_PhoBlockIndex[ele_up_num]/I");
		tree->Branch("ele_up_EleBlockIndex", stuple.ele_up_EleBlockIndex, "ele_up_EleBlockIndex[ele_up_num]/I");
		tree->Branch("ele_up_Etc"     , stuple.ele_up_Etc      , "ele_up_Etc[ele_up_num]/F");
		tree->Branch("ele_up_E"       , stuple.ele_up_E        , "ele_up_E[ele_up_num]/F");
		tree->Branch("ele_up_Px"      , stuple.ele_up_Px       , "ele_up_Px[ele_up_num]/F");
		tree->Branch("ele_up_Py"      , stuple.ele_up_Py       , "ele_up_Py[ele_up_num]/F");
		tree->Branch("ele_up_Pz"      , stuple.ele_up_Pz       , "ele_up_Pz[ele_up_num]/F");
		tree->Branch("ele_up_Detector", stuple.ele_up_Detector , "ele_up_Detector[ele_up_num]/I");
		tree->Branch("ele_up_DetEta"  , stuple.ele_up_DetEta   , "ele_up_DetEta[ele_up_num]/F");
		tree->Branch("ele_up_DetPhi"  , stuple.ele_up_DetPhi   , "ele_up_DetPhi[ele_up_num]/F");
		tree->Branch("ele_up_XCes"    , stuple.ele_up_XCes     , "ele_up_XCes[ele_up_num]/F");
		tree->Branch("ele_up_ZCes"    , stuple.ele_up_ZCes     , "ele_up_ZCes[ele_up_num]/F");
		tree->Branch("ele_up_HadEm"   , stuple.ele_up_HadEm    , "ele_up_HadEm[ele_up_num]/F");
		tree->Branch("ele_up_Chi2Mean", stuple.ele_up_Chi2Mean , "ele_up_Chi2Mean[ele_up_num]/F");
		tree->Branch("ele_up_N3d"     , stuple.ele_up_N3d      , "ele_up_N3d[ele_up_num]/I");
		tree->Branch("ele_up_Iso4"    , stuple.ele_up_Iso4     , "ele_up_Iso4[ele_up_num]/F");
		tree->Branch("ele_up_TrkIso"  , stuple.ele_up_TrkIso   , "ele_up_TrkIso[ele_up_num]/F");
		tree->Branch("ele_up_CesWireE2"      , stuple.ele_up_CesWireE2      , "ele_up_CesWireE2[ele_up_num]/F");
		tree->Branch("ele_up_CesStripE2"     , stuple.ele_up_CesStripE2     , "ele_up_CesStripE2[ele_up_num]/F");
		tree->Branch("ele_up_PhiWedge"       , stuple.ele_up_PhiWedge       , "ele_up_PhiWedge[ele_up_num]/I");
		tree->Branch("ele_up_NMuonStubs"     , stuple.ele_up_NMuonStubs     , "ele_up_NMuonStubs[ele_up_num]/I");
		tree->Branch("ele_up_EmTime"         , stuple.ele_up_EmTime         , "ele_up_EmTime[ele_up_num]/F");
		tree->Branch("ele_up_PhoenixId"      , stuple.ele_up_PhoenixId      , "ele_up_PhoenixId[ele_up_num]/I");
		tree->Branch("ele_up_Halo_seedWedge" , stuple.ele_up_Halo_seedWedge , "ele_up_Halo_seedWedge[ele_up_num]/I");
		tree->Branch("ele_up_Halo_eastNhad"  , stuple.ele_up_Halo_eastNhad  , "ele_up_Halo_eastNhad[ele_up_num]/I");
		tree->Branch("ele_up_Halo_westNhad"  , stuple.ele_up_Halo_westNhad  , "ele_up_Halo_westNhad[ele_up_num]/I");
		tree->Branch("ele_up_matchJetIndex"  , stuple.ele_up_matchJetIndex  , "ele_up_matchJetIndex[ele_up_num]/I");
		tree->Branch("ele_up_Ntracks"      , stuple.ele_up_Ntracks      , "ele_up_Ntracks[ele_up_num]/I");
		tree->Branch("ele_up_Emfr"         , stuple.ele_up_Emfr         , "ele_up_Emfr[ele_up_num]/F");
		tree->Branch("ele_up_EoverP"       , stuple.ele_up_EoverP       , "ele_up_EoverP[ele_up_num]/F");
		tree->Branch("ele_up_TrackPt"      , stuple.ele_up_TrackPt      , "ele_up_TrackPt[ele_up_num]/F");
		tree->Branch("ele_up_TrackPt2"     , stuple.ele_up_TrackPt2     , "ele_up_TrackPt2[ele_up_num]/F");
		tree->Branch("ele_up_TrackBcPt"    , stuple.ele_up_TrackBcPt    , "ele_up_TrackBcPt[ele_up_num]/F");
		tree->Branch("ele_up_TrackPhi"     , stuple.ele_up_TrackPhi     , "ele_up_TrackPhi[ele_up_num]/F");
		tree->Branch("ele_up_Nssl"         , stuple.ele_up_Nssl         , "ele_up_Nssl[ele_up_num]/I");
		tree->Branch("ele_up_Nasl"         , stuple.ele_up_Nasl         , "ele_up_Nasl[ele_up_num]/I");
		tree->Branch("ele_up_TightId"      , stuple.ele_up_TightId      , "ele_up_TightId[ele_up_num]/I");
		tree->Branch("ele_up_LooseId"      , stuple.ele_up_LooseId      , "ele_up_LooseId[ele_up_num]/I");
		tree->Branch("ele_up_ConversionId" , stuple.ele_up_ConversionId , "ele_up_ConversionId[ele_up_num]/I");
		tree->Branch("ele_up_StdTightId"      , stuple.ele_up_StdTightId      , "ele_up_StdTightId[ele_up_num]/I");
		tree->Branch("ele_up_StdLooseId"      , stuple.ele_up_StdLooseId      , "ele_up_StdLooseId[ele_up_num]/I");

			
			// Ele EM E uncertainty down 1%
		tree->Branch("ele_down_num"     , &stuple.ele_down_num      , "ele_down_num/i");
		tree->Branch("ele_down_Ntight"  , &stuple.ele_down_Ntight   , "ele_down_Ntight/i");
		tree->Branch("ele_down_Nloose"  , &stuple.ele_down_Nloose   , "ele_down_Nloose/i");
		tree->Branch("ele_down_Index"   , stuple.ele_down_Index    , "ele_down_Index[ele_down_num]/I");
		tree->Branch("ele_down_PhoBlockIndex", stuple.ele_down_PhoBlockIndex , "ele_down_PhoBlockIndex[ele_down_num]/I");
		tree->Branch("ele_down_EleBlockIndex", stuple.ele_down_EleBlockIndex, "ele_down_EleBlockIndex[ele_down_num]/I");
		tree->Branch("ele_down_Etc"     , stuple.ele_down_Etc      , "ele_down_Etc[ele_down_num]/F");
		tree->Branch("ele_down_E"       , stuple.ele_down_E        , "ele_down_E[ele_down_num]/F");
		tree->Branch("ele_down_Px"      , stuple.ele_down_Px       , "ele_down_Px[ele_down_num]/F");
		tree->Branch("ele_down_Py"      , stuple.ele_down_Py       , "ele_down_Py[ele_down_num]/F");
		tree->Branch("ele_down_Pz"      , stuple.ele_down_Pz       , "ele_down_Pz[ele_down_num]/F");
		tree->Branch("ele_down_Detector", stuple.ele_down_Detector , "ele_down_Detector[ele_down_num]/I");
		tree->Branch("ele_down_DetEta"  , stuple.ele_down_DetEta   , "ele_down_DetEta[ele_down_num]/F");
		tree->Branch("ele_down_DetPhi"  , stuple.ele_down_DetPhi   , "ele_down_DetPhi[ele_down_num]/F");
		tree->Branch("ele_down_XCes"    , stuple.ele_down_XCes     , "ele_down_XCes[ele_down_num]/F");
		tree->Branch("ele_down_ZCes"    , stuple.ele_down_ZCes     , "ele_down_ZCes[ele_down_num]/F");
		tree->Branch("ele_down_HadEm"   , stuple.ele_down_HadEm    , "ele_down_HadEm[ele_down_num]/F");
		tree->Branch("ele_down_Chi2Mean", stuple.ele_down_Chi2Mean , "ele_down_Chi2Mean[ele_down_num]/F");
		tree->Branch("ele_down_N3d"     , stuple.ele_down_N3d      , "ele_down_N3d[ele_down_num]/I");
		tree->Branch("ele_down_Iso4"    , stuple.ele_down_Iso4     , "ele_down_Iso4[ele_down_num]/F");
		tree->Branch("ele_down_TrkIso"  , stuple.ele_down_TrkIso   , "ele_down_TrkIso[ele_down_num]/F");
		tree->Branch("ele_down_CesWireE2"      , stuple.ele_down_CesWireE2      , "ele_down_CesWireE2[ele_down_num]/F");
		tree->Branch("ele_down_CesStripE2"     , stuple.ele_down_CesStripE2     , "ele_down_CesStripE2[ele_down_num]/F");
		tree->Branch("ele_down_PhiWedge"       , stuple.ele_down_PhiWedge       , "ele_down_PhiWedge[ele_down_num]/I");
		tree->Branch("ele_down_NMuonStubs"     , stuple.ele_down_NMuonStubs     , "ele_down_NMuonStubs[ele_down_num]/I");
		tree->Branch("ele_down_EmTime"         , stuple.ele_down_EmTime         , "ele_down_EmTime[ele_down_num]/F");
		tree->Branch("ele_down_PhoenixId"      , stuple.ele_down_PhoenixId      , "ele_down_PhoenixId[ele_down_num]/I");
		tree->Branch("ele_down_Halo_seedWedge" , stuple.ele_down_Halo_seedWedge , "ele_down_Halo_seedWedge[ele_down_num]/I");
		tree->Branch("ele_down_Halo_eastNhad"  , stuple.ele_down_Halo_eastNhad  , "ele_down_Halo_eastNhad[ele_down_num]/I");
		tree->Branch("ele_down_Halo_westNhad"  , stuple.ele_down_Halo_westNhad  , "ele_down_Halo_westNhad[ele_down_num]/I");
		tree->Branch("ele_down_matchJetIndex"  , stuple.ele_down_matchJetIndex  , "ele_down_matchJetIndex[ele_down_num]/I");
		tree->Branch("ele_down_Ntracks"      , stuple.ele_down_Ntracks      , "ele_down_Ntracks[ele_down_num]/I");
		tree->Branch("ele_down_Emfr"         , stuple.ele_down_Emfr         , "ele_down_Emfr[ele_down_num]/F");
		tree->Branch("ele_down_EoverP"       , stuple.ele_down_EoverP       , "ele_down_EoverP[ele_down_num]/F");
		tree->Branch("ele_down_TrackPt"      , stuple.ele_down_TrackPt      , "ele_down_TrackPt[ele_down_num]/F");
		tree->Branch("ele_down_TrackPt2"     , stuple.ele_down_TrackPt2     , "ele_down_TrackPt2[ele_down_num]/F");
		tree->Branch("ele_down_TrackBcPt"    , stuple.ele_down_TrackBcPt    , "ele_down_TrackBcPt[ele_down_num]/F");
		tree->Branch("ele_down_TrackPhi"     , stuple.ele_down_TrackPhi     , "ele_down_TrackPhi[ele_down_num]/F");
		tree->Branch("ele_down_Nssl"         , stuple.ele_down_Nssl         , "ele_down_Nssl[ele_down_num]/I");
		tree->Branch("ele_down_Nasl"         , stuple.ele_down_Nasl         , "ele_down_Nasl[ele_down_num]/I");
		tree->Branch("ele_down_TightId"      , stuple.ele_down_TightId      , "ele_down_TightId[ele_down_num]/I");
		tree->Branch("ele_down_LooseId"      , stuple.ele_down_LooseId      , "ele_down_LooseId[ele_down_num]/I");
		tree->Branch("ele_down_ConversionId" , stuple.ele_down_ConversionId , "ele_down_ConversionId[ele_down_num]/I");
		tree->Branch("ele_down_StdTightId"      , stuple.ele_down_StdTightId      , "ele_down_StdTightId[ele_down_num]/I");
		tree->Branch("ele_down_StdLooseId"      , stuple.ele_down_StdLooseId      , "ele_down_StdLooseId[ele_down_num]/I");


		// JET COLLECTIONS
		tree->Branch("jet_num"      , &stuple.jet_num      , "jet_num/i");
		tree->Branch("jet_NJet15"   , &stuple.jet_NJet15   , "jet_Njet15/I");
		tree->Branch("jet_Index"    , stuple.jet_Index    , "jet_Index[jet_num]/I");
		tree->Branch("jet_Pt"       , stuple.jet_Pt       , "jet_Pt[jet_num]/F");
		tree->Branch("jet_E"        , stuple.jet_E        , "jet_E[jet_num]/F");
		tree->Branch("jet_Px"       , stuple.jet_Px       , "jet_Px[jet_num]/F");
		tree->Branch("jet_Py"       , stuple.jet_Py       , "jet_Py[jet_num]/F");
		tree->Branch("jet_Pz"       , stuple.jet_Pz       , "jet_Pz[jet_num]/F");
		tree->Branch("jet_DetEta"   , stuple.jet_DetEta   , "jet_DetEta[jet_num]/F");
		tree->Branch("jet_DetPhi"   , stuple.jet_DetPhi   , "jet_DetPhi[jet_num]/F");
		tree->Branch("jet_HadEm"    , stuple.jet_HadEm    , "jet_HadEm[jet_num]/F");
		tree->Branch("jet_Emfr"     , stuple.jet_Emfr     , "jet_Emfr[jet_num]/F");
		tree->Branch("jet_Ntowers"  , stuple.jet_Ntowers  , "jet_Ntowers[jet_num]/I");
		tree->Branch("jet_Ntracks"  , stuple.jet_Ntracks  , "jet_Ntracks[jet_num]/I");
		tree->Branch("jet_SeedIPhi" , stuple.jet_SeedIPhi , "jet_SeedIPhi[jet_num]/I");
		tree->Branch("jet_SeedIEta" , stuple.jet_SeedIEta , "jet_SeedIEta[jet_num]/I");
		tree->Branch("jet_EmTime"   , stuple.jet_EmTime   , "jet_EmTime[jet_num]/F");
		tree->Branch("jet_SecVtxTag" , stuple.jet_SecVtxTag , "jet_SecVtxTag[jet_num]/I");
		tree->Branch("jet_SecVtxppb" , stuple.jet_SecVtxppb  , "jet_SecVtxppb[jet_num]/F");
		tree->Branch("jet_SecVtxnpb" , stuple.jet_SecVtxnpb  , "jet_SecVtxnpb[jet_num]/F");
		tree->Branch("jet_SecVtxTrkmass" , stuple.jet_SecVtxTrkmass , "jet_SecVtxTrkmass[jet_num]/F");

		// RAW JET COLLECTIONS
		tree->Branch("jet_raw_num"   , &stuple.jet_raw_num   , "jet_raw_num/i");
		tree->Branch("jet_raw_Index" , stuple.jet_raw_Index , "jet_raw_Index[jet_raw_num]/I");
		tree->Branch("jet_raw_Pt"    , stuple.jet_raw_Pt    , "jet_raw_Pt[jet_raw_num]/F");
		tree->Branch("jet_raw_E"     , stuple.jet_raw_E     , "jet_raw_E[jet_raw_num]/F");
		tree->Branch("jet_raw_Px"    , stuple.jet_raw_Px    , "jet_raw_Px[jet_raw_num]/F");
		tree->Branch("jet_raw_Py"    , stuple.jet_raw_Py    , "jet_raw_Py[jet_raw_num]/F");
		tree->Branch("jet_raw_Pz"    , stuple.jet_raw_Pz    , "jet_raw_Pz[jet_raw_num]/F");


		
			// JES UP JET COLLECTION
		tree->Branch("jet_up_num"     , &stuple.jet_up_num      , "jet_up_num/i");
		tree->Branch("jet_up_NJet15"  , &stuple.jet_up_NJet15   , "jet_up_NJet15/i");
		tree->Branch("jet_up_Index"   , stuple.jet_up_Index    , "jet_up_Index[jet_up_num]/I");
		tree->Branch("jet_up_Pt"      , stuple.jet_up_Pt       , "jet_up_Pt[jet_up_num]/F");
		tree->Branch("jet_up_E"       , stuple.jet_up_E        , "jet_up_E[jet_up_num]/F");
		tree->Branch("jet_up_Px"      , stuple.jet_up_Px       , "jet_up_Px[jet_up_num]/F");
		tree->Branch("jet_up_Py"      , stuple.jet_up_Py       , "jet_up_Py[jet_up_num]/F");
		tree->Branch("jet_up_Pz"      , stuple.jet_up_Pz       , "jet_up_Pz[jet_up_num]/F");
		tree->Branch("jet_up_DetEta"  , stuple.jet_up_DetEta   , "jet_up_DetEta[jet_up_num]/F");
		tree->Branch("jet_up_DetPhi"  , stuple.jet_up_DetPhi   , "jet_up_DetPhi[jet_up_num]/F");
		tree->Branch("jet_up_HadEm"   , stuple.jet_up_HadEm    , "jet_up_HadEm[jet_up_num]/F");
		tree->Branch("jet_up_Emfr"    , stuple.jet_up_Emfr     , "jet_up_Emfr[jet_up_num]/F");
		tree->Branch("jet_up_Ntowers" , stuple.jet_up_Ntowers  , "jet_up_Ntowers[jet_up_num]/I");
		tree->Branch("jet_up_Ntracks" , stuple.jet_up_Ntracks  , "jet_up_Ntracks[jet_up_num]/I");
		tree->Branch("jet_up_SeedIPhi", stuple.jet_up_SeedIPhi , "jet_up_SeedIPhi[jet_up_num]/I");
		tree->Branch("jet_up_SeedIEta", stuple.jet_up_SeedIEta , "jet_up_SeedIEta[jet_up_num]/I");
		tree->Branch("jet_up_EmTime"  , stuple.jet_up_EmTime   , "jet_up_EmTime[jet_up_num]/F");
		tree->Branch("jet_up_SecVtxTag" , stuple.jet_up_SecVtxTag , "jet_up_SecVtxTag[jet_num]/I");
		tree->Branch("jet_up_SecVtxppb" , stuple.jet_up_SecVtxppb  , "jet_up_SecVtxppb[jet_num]/F");
		tree->Branch("jet_up_SecVtxnpb" , stuple.jet_up_SecVtxnpb  , "jet_up_SecVtxnpb[jet_num]/F");
		tree->Branch("jet_up_SecVtxTrkmass" , stuple.jet_up_SecVtxTrkmass , "jet_up_SecVtxTrkmass[jet_num]/F");
			
			// JES DOWN JET COLLECTION
		tree->Branch("jet_down_num"     , &stuple.jet_down_num      , "jet_down_num/i");
		tree->Branch("jet_down_NJet15"  , &stuple.jet_down_NJet15   , "jet_down_NJet15/i");
		tree->Branch("jet_down_Index"   , stuple.jet_down_Index    , "jet_down_Index[jet_down_num]/I");
		tree->Branch("jet_down_Pt"      , stuple.jet_down_Pt       , "jet_down_Pt[jet_down_num]/F");
		tree->Branch("jet_down_E"       , stuple.jet_down_E        , "jet_down_E[jet_down_num]/F");
		tree->Branch("jet_down_Px"      , stuple.jet_down_Px       , "jet_down_Px[jet_down_num]/F");
		tree->Branch("jet_down_Py"      , stuple.jet_down_Py       , "jet_down_Py[jet_down_num]/F");
		tree->Branch("jet_down_Pz"      , stuple.jet_down_Pz       , "jet_down_Pz[jet_down_num]/F");
		tree->Branch("jet_down_DetEta"  , stuple.jet_down_DetEta   , "jet_down_DetEta[jet_down_num]/F");
		tree->Branch("jet_down_DetPhi"  , stuple.jet_down_DetPhi   , "jet_down_DetPhi[jet_down_num]/F");
		tree->Branch("jet_down_HadEm"   , stuple.jet_down_HadEm    , "jet_down_HadEm[jet_down_num]/F");
		tree->Branch("jet_down_Emfr"    , stuple.jet_down_Emfr     , "jet_down_Emfr[jet_down_num]/F");
		tree->Branch("jet_down_Ntowers" , stuple.jet_down_Ntowers  , "jet_down_Ntowers[jet_down_num]/I");
		tree->Branch("jet_down_Ntracks" , stuple.jet_down_Ntracks  , "jet_down_Ntracks[jet_down_num]/I");
		tree->Branch("jet_down_SeedIPhi", stuple.jet_down_SeedIPhi , "jet_down_SeedIPhi[jet_down_num]/I");
		tree->Branch("jet_down_SeedIEta", stuple.jet_down_SeedIEta , "jet_down_SeedIEta[jet_down_num]/I");
		tree->Branch("jet_down_EmTime"  , stuple.jet_down_EmTime   , "jet_down_EmTime[jet_down_num]/F");
		tree->Branch("jet_down_SecVtxTag" , stuple.jet_down_SecVtxTag , "jet_down_SecVtxTag[jet_num]/I");
		tree->Branch("jet_down_SecVtxppb" , stuple.jet_down_SecVtxppb  , "jet_down_SecVtxppb[jet_num]/F");
		tree->Branch("jet_down_SecVtxnpb" , stuple.jet_down_SecVtxnpb  , "jet_down_SecVtxnpb[jet_num]/F");
		tree->Branch("jet_down_SecVtxTrkmass" , stuple.jet_down_SecVtxTrkmass , "jet_down_SecVtxTrkmass[jet_num]/F");


			// VERTEX INFO
		tree->Branch("vtx_N"        , &stuple.vtx_N        , "vtx_N/I");
		tree->Branch("vtx_NClass12" , &stuple.vtx_NClass12 , "vtx_NClass12/I");
		tree->Branch("vtx_z"        , &stuple.vtx_z        , "vtx_z/F");
		tree->Branch("vtx_Ntracks"  , &stuple.vtx_Ntracks  , "vtx_Ntracks/I");
		tree->Branch("vtx_SumPt"    , &stuple.vtx_SumPt    , "vtx_SumPt/F");

			// MET INFO
		tree->Branch("met_Met"    , &stuple.met_Met    , "met_Met/F");
		tree->Branch("met_RawMet"    , &stuple.met_RawMet    , "met_RawMet/F");
		tree->Branch("met_MetX"   , &stuple.met_MetX   , "met_MetX/F");
		tree->Branch("met_MetY"   , &stuple.met_MetY   , "met_MetY/F");
		tree->Branch("met_SumEt"  , &stuple.met_SumEt  , "met_SumEt/F");
		tree->Branch("met_Ht"     , &stuple.met_Ht     , "met_Ht/F");
		tree->Branch("met_MetPhi" , &stuple.met_MetPhi , "met_MetPhi/F");

			// GENERATED MET FROM MET MODEL
//		tree->Branch("met_Gen_d"   , &stuple.met_Gen_d   , "met_Gen_d/F");
//		tree->Branch("met_Gen_m"   , &stuple.met_Gen_m   , "met_Gen_m/F");
//		tree->Branch("met_Gen_p"   , &stuple.met_Gen_p   , "met_Gen_p/F");
//		tree->Branch("met_Gen_mUn" , &stuple.met_Gen_mUn , "met_Gen_mUn/F");
//		tree->Branch("met_Gen_pUn" , &stuple.met_Gen_pUn , "met_Gen_pUn/F");

		
			// Generator level info
		tree->Branch("gen_elenum" , &stuple.gen_elenum , "gen_elenum/i");
		tree->Branch("gen_phonum" , &stuple.gen_phonum , "gen_phonum/i");
		
		tree->Branch("gen_MomIndex" , &stuple.gen_MomIndex  , "gen_MomIndex/I");
		tree->Branch("gen_MomPDG"   , &stuple.gen_MomPDG    , "gen_MomPDG/I");
		tree->Branch("gen_MomStaus" , &stuple.gen_MomStatus , "gen_MomStatus/I");
		tree->Branch("gen_MomEtc"   , &stuple.gen_MomEtc    , "gen_MomEtc/F");
		tree->Branch("gen_MomE"     , &stuple.gen_MomE      , "gen_MomE/F");
		tree->Branch("gen_MomPx"    , &stuple.gen_MomPx     , "gen_MomPx/F");
		tree->Branch("gen_MomPy"    , &stuple.gen_MomPy     , "gen_MomPy/F");
		tree->Branch("gen_MomPz"    , &stuple.gen_MomPz     , "gen_MomPz/F");

		tree->Branch("gen_ProdVtxX" , &stuple.gen_ProdVtxX , "gen_ProdVtxX/F");
		tree->Branch("gen_ProdVtxY" , &stuple.gen_ProdVtxY , "gen_ProdVtxY/F");
		tree->Branch("gen_ProdVtxZ" , &stuple.gen_ProdVtxZ , "gen_ProdVtxZ/F");
		tree->Branch("gen_ProdVtxT" , &stuple.gen_ProdVtxT , "gen_ProdVtxT/F");

		tree->Branch("gen_pho_Index"  , stuple.gen_pho_Index  , "gen_pho_Index[gen_phonum]/I");
		tree->Branch("gen_pho_PDG"    , stuple.gen_pho_PDG    , "gen_pho_PDG[gen_phonum]/I");
		tree->Branch("gen_pho_Status" , stuple.gen_pho_Status , "gen_pho_Status[gen_phonum]/I");
		tree->Branch("gen_pho_Etc"    , stuple.gen_pho_Etc    , "gen_pho_Etc[gen_phonum]/F");
		tree->Branch("gen_pho_E"      , stuple.gen_pho_E      , "gen_pho_E[gen_phonum]/F");
		tree->Branch("gen_pho_Px"     , stuple.gen_pho_Px     , "gen_pho_Px[gen_phonum]/F");
		tree->Branch("gen_pho_Py"     , stuple.gen_pho_Py     , "gen_pho_Py[gen_phonum]/F");
		tree->Branch("gen_pho_Pz"     , stuple.gen_pho_Pz     , "gen_pho_Pz[gen_phonum]/F");

		tree->Branch("gen_ele_Index"  , stuple.gen_ele_Index  , "gen_ele_Index[gen_elenum]/I");
		tree->Branch("gen_ele_PDG"    , stuple.gen_ele_PDG    , "gen_ele_PDG[gen_elenum]/I");
		tree->Branch("gen_ele_Status" , stuple.gen_ele_Status , "gen_ele_Status[gen_elenum]/I");
		tree->Branch("gen_ele_Etc"    , stuple.gen_ele_Etc    , "gen_ele_Etc[gen_elenum]/F");
		tree->Branch("gen_ele_E"      , stuple.gen_ele_E      , "gen_ele_E[gen_elenum]/F");
		tree->Branch("gen_ele_Px"     , stuple.gen_ele_Px     , "gen_ele_Px[gen_elenum]/F");
		tree->Branch("gen_ele_Py"     , stuple.gen_ele_Py     , "gen_ele_Py[gen_elenum]/F");
		tree->Branch("gen_ele_Pz"     , stuple.gen_ele_Pz     , "gen_ele_Pz[gen_elenum]/F");


	return 0;
}

//_____________________________________________________________________________
int FlatStupleMaker::BeginRun()
{
	bMc = GetHeaderBlock()->McFlag();

	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int FlatStupleMaker::Event(int ientry)
{
	SetPassed(1);
	counter.evtsProcessed++;
	if (! bRunPermit) 
	{
		StdOut(__FILE__,__LINE__,3,"Run premit failed. A required module is not present. Pl. check!.");
		exit (1);
	}

	//fHeaderBlock->Print();
	
	//to get b-tag info
	float zVertex = 0.0;
	int nzVertex = 0;
	TStnEvent* event = GetEvent();
   bool gVertex = TPhotonUtil::SelectVertex(event, zVertex, nzVertex); 
	TUtil::FillJetTags(event, matrix, zVertex, nzVertex);

	
	spUp.clear();
	spDown.clear();

	bCosmicEvent = false;
	bJetEmMismatch = false;
	
	// first copy central photons to EM up/down lists

	for (int i=0; i < initSpMod->GetSuperPhoSize(); i++) {
		spUp.push_back(*(initSpMod->GetSuperPhoton(i)));
		spDown.push_back(*(initSpMod->GetSuperPhoton(i)));
	}

	
	// now tag the new photon lists
	TagForEmSys(spUp,1);  		//shift energy up 1% and retag
	TagForEmSys(spDown,-1); 	// shift energy down 1% and retag

	// now do the same steps as the central photons to make the flat ntuples!

	stuple.Init();

	std::vector<int> vPho,vEle, vPhoUp, vEleUp, vPhoDown, vEleDown;
	
	// first save the central pho list
	// this is little different from up/down list info access
	// these pho/ele will get info from SuperPhoton list
	// in the InitSuperPhotons module directly.
	// up/dow list are local modified copies of that list.
	for (int i=0; i < initSpMod->GetSuperPhoSize(); i++) {
		if (initSpMod->GetSuperPhoton(i)->GetEtCorr() < fMinEtThr) continue;	//Et cut for em objects. tagging is done for Et>7 em objects
		if ( initSpMod->GetSuperPhoton(i)->IsTightPhoton() 
			|| initSpMod->GetSuperPhoton(i)->IsLoosePhoton())
		{
			stuple.pho_num++;
			if (initSpMod->GetSuperPhoton(i)->IsTightPhoton())	stuple.pho_Ntight++;
			if (initSpMod->GetSuperPhoton(i)->IsLoosePhoton()) stuple.pho_Nloose++;
			vPho.push_back(i);
		} else if (initSpMod->GetSuperPhoton(i)->IsTightElectron() //this way I'll remove the objects passing both loose pho and ele cuts 01-21-2010
			|| initSpMod->GetSuperPhoton(i)->IsLooseElectron()
			//|| initSpMod->GetSuperPhoton(i)->IsStdTightElectron()
			|| initSpMod->GetSuperPhoton(i)->IsStdLooseElectron()) // this is a temp solution instead of creating a whole new electron collection for std. electron cuts
		{

			stuple.ele_num++;
			if (initSpMod->GetSuperPhoton(i)->IsTightElectron()
				 //|| initSpMod->GetSuperPhoton(i)->IsStdTightElectron()
				 ) stuple.ele_Ntight++;
			if (initSpMod->GetSuperPhoton(i)->IsLooseElectron()
				|| initSpMod->GetSuperPhoton(i)->IsStdLooseElectron() ) stuple.ele_Nloose++;
			vEle.push_back(i);
		}
		//float fEmTime = initSpMod->GetSuperPhoton(i)->GetEmTime();
		//if (fEmTime >20 && fEmTime < 100) bCosmicEvent = true;
	}

	Hist.Npho->Fill(stuple.pho_num);
	Hist.Nele->Fill(stuple.ele_num);

	for (unsigned int i=0; i < spUp.size(); i++) {
		if (spUp[i].GetEtCorr() < fMinEtThr) continue;	//Et cut for em objects. tagging is done for Et>7 em objects
		if (spUp[i].IsTightPhoton() || spUp[i].IsLoosePhoton()) {
			stuple.pho_up_num++;
			if (spUp[i].IsTightPhoton()) stuple.pho_up_Ntight++;
			if (spUp[i].IsLoosePhoton()) stuple.pho_up_Nloose++;
			vPhoUp.push_back(i);
		}
		if (spUp[i].IsTightElectron() || spUp[i].IsLooseElectron()) {
			stuple.ele_up_num++;
			if (spUp[i].IsTightElectron()) { 
				//|| spUp[i].IsStdTightElectron()
				 stuple.ele_up_Ntight++;
			}
			if (spUp[i].IsLooseElectron() || spUp[i].IsStdLooseElectron()) stuple.ele_up_Nloose++;
			vEleUp.push_back(i);
		}
	}

	for (unsigned int i=0; i < spDown.size(); i++) {
		if (spDown[i].GetEtCorr() < fMinEtThr) continue;	//Et cut for em objects. tagging is done for Et>7 em objects
		if (spDown[i].IsTightPhoton() || spDown[i].IsLoosePhoton()) {
			stuple.pho_down_num++;
			if (spDown[i].IsTightPhoton()) stuple.pho_down_Ntight++;
			if (spDown[i].IsLoosePhoton()) stuple.pho_down_Nloose++;
			vPhoDown.push_back(i);
		}
		if (spDown[i].IsTightElectron() || spDown[i].IsLooseElectron()) {
			stuple.ele_down_num++;
			if (spDown[i].IsTightElectron()) { 
				//|| spDown[i].IsStdTightElectron()
				 stuple.ele_down_Ntight++;
			}
			if (spDown[i].IsLooseElectron() || spDown[i].IsStdLooseElectron()) stuple.ele_down_Nloose++;
			vEleDown.push_back(i);
		}
	}

	if (stuple.pho_num > 0 && stuple.ele_num == 0) counter.phoOnlyEvts++;
	if (stuple.pho_num == 0 && stuple.ele_num > 0) counter.eleOnlyEvts++;
	if (stuple.pho_num > 0 && stuple.ele_num > 0) counter.phoEleEvts++;

	assert (stuple.pho_num <= Stuple::Npho);
	assert (stuple.ele_num <= Stuple::Nele);
	assert (stuple.pho_up_num <= Stuple::Npho);
	assert (stuple.ele_up_num <= Stuple::Nele);

	if ( (stuple.pho_num + stuple.ele_num +
		 stuple.pho_up_num + stuple.ele_up_num +
		 stuple.pho_down_num + stuple.ele_down_num > 0)  //any em object
		//drop the jet requirement as for BH template/and rejections power we need no vtx /no jet events 02-12-2010
		// and a jet Et>15GeV 
		// I figured that no matter what the setting I used (njet15==0) the the final jet specturm starts from Et>15
		// an i cannot figure out where this cuts is imposed within the JetFilterModule. I changed the MET scenario,
		// MinEtThr=0 etc but it did not change the out come.
		//
		//&& ( jetMod->GetMyNjet(0,3)>=1     //drop the jet requirement as for BH template/and rejections
		//	|| jetModUp->GetMyNjet(0,3)>=1  //power we need no vtx /no jet events 02-12-2010
		//	|| jetModDown->GetMyNjet(0,3)>=1 )		//and a jet Et>15GeV 
		 )
	{
		FillDataBlocks(ientry);
		GetCommonValues();										// get the event header info, run,event,vtx etc	
		if (stuple.pho_num + stuple.ele_num > 0)
		{
			//GetHeaderBlock()->Print();
			GetCentralValues(vPho, vEle);						// get central pho/ele info
		}
		if (stuple.pho_up_num + stuple.ele_up_num > 0) GetEmUpValues(vPhoUp, vEleUp);				// get em up pho/ele info
		if (stuple.pho_down_num + stuple.ele_down_num > 0) GetEmDownValues(vPhoDown, vEleDown);// get em down pho/ele info
		
		GetRawJetInfo();
		
		//discard any event that had problems in removing EM Objects from the jet list
		if (! bJetEmMismatch) {
			tree->Fill();
			counter.evtsWritten++;
			//stuple.Dump();

		} else { 
			std::cout << "JETEMMISMATCH::" << GetHeaderBlock()->RunNumber() << ", " << GetHeaderBlock()->EventNumber() << std::endl;
			SetPassed(0);		// save this event Stntuple format
		}
		
	}

/*	std::cout << red << " jet em time " << std::endl;
	
	for (int i=0; i < fJetBlock->NJets(); ++i)
	{
		TStnJet *jet = fJetBlock->Jet(i);
		float time = GetEmTiming(i);
		std::cout << i << "\t" << jet->Et() << "/" << jet->Eta() << " --> " << time << std::endl;
	}
	std::cout << "\n" << clearatt << std::endl;
*/
	if (GetPassed()) counter.evtsPassModule++;
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// for debugging only
/*-------------------------------------------------------------------*/
void FlatStupleMaker::DumpSpAndJetBlock()
{
		std::cout << "=============== photon dump " << std::endl;
		std::cout << " i\t E \t  Et \t Eta \t Phi " << std::endl;
		for (int i=0; i < initSpMod->GetSuperPhoSize(); ++i)
		{
			SuperPhoton* sp = initSpMod->GetSuperPhoton(i);
			if (sp->IsTightPhoton()
				|| sp->IsLoosePhoton()
				|| sp->IsTightElectron()
				|| sp->IsLooseElectron()
				//|| sp->IsStdTightElectron()
				|| sp->IsStdLooseElectron() )
				{
					std::cout << i << "\t" << sp->GetECorr() << "\t" << sp->GetEtCorr() 
									<< "\t" << sp->GetDetEta() << "\t" << sp->GetDetPhi()
									<< std::endl;
				}

		}

		std::cout << "\t ===== jet dump " << std::endl;
		for (int i =0; i < fJetBlock->NJets(); ++i)
		{
			TStnJet* j = fJetBlock->Jet(i);
			//TLorentzVector* tlj = j->Momentum();
			std::cout << i << "\t" << j->Et() << "\t" << j->EtCorr() << "\t"
						<< "\t" << j->DetEta() << "\t" << j->Phi() << std::endl;


		}
		std::cout << " ==== after removing em objs " << std::endl;
		std::cout << "\t ===== jet dump after removing em objs" << std::endl;
		std::cout << "pho, ele, jet num = " << stuple.pho_num << "\t"
						<< stuple.ele_num << "\t" << stuple.jet_num << std::endl;
		for (unsigned int i =0; i < stuple.jet_num; ++i)
		{
			std::cout << i << "\t" << stuple.jet_E[i] << "\t" << stuple.jet_Pt[i] << "\t"
						<< "\t" << stuple.jet_DetEta[i]  << std::endl;


		}


		
	Hist.Ncj->Fill(stuple.jet_num);	
	Hist.Nupj->Fill(stuple.jet_up_num);	
	Hist.Ndownj->Fill(stuple.jet_down_num);	
	Hist.NcjUpj->Fill(stuple.jet_num, stuple.jet_up_num);	
	Hist.NcjDownj->Fill(stuple.jet_up_num, stuple.jet_down_num);	

	
	float sumPt = 0, sumPt_u=0, sumPt_d=0;
	for (unsigned int i = 0; i < stuple.jet_num; ++i) sumPt += stuple.jet_Pt[i];
	for (unsigned int i = 0; i < stuple.jet_up_num; ++i) sumPt_u += stuple.jet_up_Pt[i];
	for (unsigned int i = 0; i < stuple.jet_down_num; ++i) sumPt_d += stuple.jet_down_Pt[i];
	Hist.cj_EperJet->Fill(sumPt/(float)stuple.jet_num);	
	Hist.upj_EperJet->Fill(sumPt_u/(float)stuple.jet_up_num);	
	Hist.downj_EperJet->Fill(sumPt_d/(float)stuple.jet_down_num);	
		
}

/*-------------------------------------------------------------------*/
// for debugging only
/*-------------------------------------------------------------------*/
void FlatStupleMaker::CheckLooseElectronTag()
{
	// check if the electron tagged using the std. loose electron id
	// is acceptable by comparing it to photon-like electron id
	// and by looking at the W/Z peaks.

	std::vector<int> vInd;
	for (int i=0; i < initSpMod->GetSuperPhoSize(); ++i)
	{
		if (initSpMod->GetSuperPhoton(i)->IsConversion()) continue; 
		if (initSpMod->GetSuperPhoton(i)->IsStdLooseElectron()) vInd.push_back(i);
	}

	// select w events
	if (vInd.size() == 1)
	{
		if (initSpMod->GetSuperPhoton(vInd.at(0))->GetDetector() == 0)
			Hist.w_met_c->Fill(jetMod->GetMyMetCorr(0,3));
		if (initSpMod->GetSuperPhoton(vInd.at(0))->GetDetector() == 1)
			Hist.w_met_p->Fill(jetMod->GetMyMetCorr(0,3));
	}
	// select z events
	if (vInd.size() >= 2)
	{
		int iDet_e1 = initSpMod->GetSuperPhoton(vInd.at(0))->GetDetector();
		int iDet_e2 = initSpMod->GetSuperPhoton(vInd.at(1))->GetDetector();
		if (iDet_e1+ iDet_e2 == 0)
			Hist.z_met_cc->Fill(jetMod->GetMyMetCorr(0,3));
		if ((iDet_e1 + iDet_e2) == 1)
			Hist.z_met_cp->Fill(jetMod->GetMyMetCorr(0,3));
		if (iDet_e1 + iDet_e2 == 2)
			Hist.z_met_pp->Fill(jetMod->GetMyMetCorr(0,3));
	}


}

/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetCommonValues()
{
	stuple.evt_McFlag = GetHeaderBlock()->McFlag();
	stuple.evt_RunNumber = GetHeaderBlock()->RunNumber();
	stuple.evt_EventNumber = GetHeaderBlock()->EventNumber();
	
	if (!bMc) GetTriggerInfo();
	GetVertexInfo();
	GetMetInfo();
	if (bMc) GetGenLevelInfo();
	//if (stuple.jet_num + stuple.jet_up_num + stuple.jet_down_num)
	//DumpSpAndJetBlock();  //for debugging only
	CheckLooseElectronTag();
}


/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetCentralValues(std::vector<int> vPho, std::vector<int> vEle)
{
	for (unsigned int i=0; i< vPho.size(); i++) {
		GetPhotonInfo(i, *(initSpMod->GetSuperPhoton(vPho[i])));
		if (! GetHeaderBlock()->McFlag())	//NO BH tagging for MC samples -01-10-2010
		{
			TagBeamHalo::HaloStuff aHaloPara = tagHaloMod->GetHaloStuff(vPho[i]);
			stuple.pho_Halo_seedWedge[i] = aHaloPara.seedwedge;
			stuple.pho_Halo_eastNhad[i] = aHaloPara.eastNhad;
			stuple.pho_Halo_westNhad[i] = aHaloPara.westNhad;
		}
	}

	for (unsigned int i=0; i< vEle.size(); i++) {
		GetElectronInfo(i, *(initSpMod->GetSuperPhoton(vEle[i])));
		if (! GetHeaderBlock()->McFlag())
		{
			TagBeamHalo::HaloStuff aHaloPara = tagHaloMod->GetHaloStuff(vEle[i]);
			stuple.ele_Halo_seedWedge[i] = aHaloPara.seedwedge;
			stuple.ele_Halo_eastNhad[i] = aHaloPara.eastNhad;
			stuple.ele_Halo_westNhad[i] = aHaloPara.westNhad;
		}
	}

	GetJetInfo();

}

/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetEmUpValues(std::vector<int> vPhoUpIndex, std::vector<int> vEleUpIndex)
{
	// vPhoUpIndex and vEleUpIndex holds the indices of the SuperPhotons

	for (unsigned int i=0; i< vPhoUpIndex.size(); i++) {		// this is the index of the photon in the SuperPhoton list which is sorted in Et, highest to lowest
		GetEmUpPhotonInfo(i, spUp[vPhoUpIndex[i]]);
		if (! GetHeaderBlock()->McFlag())	//NO BH tagging for MC samples -01-10-2010
		{
			TagBeamHalo::HaloStuff aHaloPara = tagHaloMod->GetHaloStuff(vPhoUpIndex[i]);
			stuple.pho_up_Halo_seedWedge[i] = aHaloPara.seedwedge;
			stuple.pho_up_Halo_eastNhad[i] = aHaloPara.eastNhad;
			stuple.pho_up_Halo_westNhad[i] = aHaloPara.westNhad;
		}
	}

	for (unsigned int i=0; i< vEleUpIndex.size(); i++) {
		GetEmUpElectronInfo(i, spUp[vEleUpIndex[i]]);
		if (! GetHeaderBlock()->McFlag())
		{
			TagBeamHalo::HaloStuff aHaloPara = tagHaloMod->GetHaloStuff(vEleUpIndex[i]);
			stuple.ele_up_Halo_seedWedge[i] = aHaloPara.seedwedge;
			stuple.ele_up_Halo_eastNhad[i] = aHaloPara.eastNhad;
			stuple.ele_up_Halo_westNhad[i] = aHaloPara.westNhad;
		}
	}

	GetJesUpJetInfo();

}


/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetEmDownValues(std::vector<int> vPhoDownIndex, std::vector<int> vEleDownIndex)
{
	// vPhoDownIndex and vEleDownIndex holds the indices of the SuperPhotons

	for (unsigned int i=0; i< vPhoDownIndex.size(); i++) {		// this is the index of the photon in the SuperPhoton list which is sorted in Et, highest to lowest
		GetEmDownPhotonInfo(i, spUp[vPhoDownIndex[i]]);
		if (! GetHeaderBlock()->McFlag())	//NO BH tagging for MC samples -01-10-2010
		{
			TagBeamHalo::HaloStuff aHaloPara = tagHaloMod->GetHaloStuff(vPhoDownIndex[i]);
			stuple.pho_down_Halo_seedWedge[i] = aHaloPara.seedwedge;
			stuple.pho_down_Halo_eastNhad[i] = aHaloPara.eastNhad;
			stuple.pho_down_Halo_westNhad[i] = aHaloPara.westNhad;
		}
	}

	for (unsigned int i=0; i< vEleDownIndex.size(); i++) {
		GetEmDownElectronInfo(i, spUp[vEleDownIndex[i]]);
		if (! GetHeaderBlock()->McFlag())
		{
			TagBeamHalo::HaloStuff aHaloPara = tagHaloMod->GetHaloStuff(vEleDownIndex[i]);
			stuple.ele_down_Halo_seedWedge[i] = aHaloPara.seedwedge;
			stuple.ele_down_Halo_eastNhad[i] = aHaloPara.eastNhad;
			stuple.ele_down_Halo_westNhad[i] = aHaloPara.westNhad;
		}
	}

	GetJesDownJetInfo();
}


/*-------------------------------------------------------------------*/
// iSys=+1/-1 ->shift energy up/down
// Here I am only varying the energy bu 1% and checking if it 
// pass/fail the Et cut. But I should check all the other cuts that
// depends on energy (HadEm,Iso,TrkPt,TrkIso,2nd CES cluster).
// But this is very small effect. So lets not worry about it for now
// - Sep 7, 2009
/*-------------------------------------------------------------------*/
void FlatStupleMaker::TagForEmSys(std::vector<SuperPhoton>& vSp, int iSys)
{
	if (iSys != 1 && iSys != -1) {
		StdOut(__FILE__,__LINE__,2,"EM systematic must be either +1 or -1. pl. check!");
		return;
	}

	double zvx  = 0;
	if (fVertexBlock->GetBestVertex(12,1) != NULL) {
			zvx = fVertexBlock->GetBestVertex(12,1)->Z();  //highest sum pt vertex
	}

	TStnEvent* event = GetEvent();

	for (unsigned i=0; i< vSp.size(); ++i) {
		//first make local temp copy of the SuperPhoton
		SuperPhoton sp = vSp[i];
		TLorentzVector tl(sp.GetCorVec());

		TagTightPhotons tagTightPho("TagEmSysTightPhotons","TagEmSysTightPhotons");
		TagLoosePhotons tagLoosePho("TagEmSysLoosePhotons","TagEmSysLoosePhotons");

		if (iSys == 1) tl = tl + tl * 0.01;   			// 1% systematic in Energy measurement
		else if (iSys == -1) tl = tl - tl * 0.01;    // 1% systematic in Energy measurement
		
		// now set the update the SuperPhoton with new values
		sp.SetCorVec(tl);
		sp.SetEtCorr(tl.E() * TMath::Sin(tl.Theta()));
		sp.SetECorr(tl.E());
		
		int iPhoTightId = -1, iPhoLooseId = -1;
		if (sp.GetDetector() == 0)
		{
			iPhoTightId = tagTightPho.GetCenTightIdWord(&sp);
			iPhoLooseId = tagLoosePho.GetCenLooseIdWord(&sp);
			sp.SetTightPhotonId(iPhoTightId);
			sp.SetLoosePhotonId(iPhoLooseId);
		} else if (sp.GetDetector() == 1)
		{
			iPhoTightId = tagTightPho.GetPlugTightIdWord(&sp);
			iPhoLooseId = tagLoosePho.GetPlugLooseIdWord(&sp);
			sp.SetTightPhotonId(iPhoTightId);
			sp.SetLoosePhotonId(iPhoLooseId);
		}
		
		TagElectrons tagTightEle("TagEmSysTightElectrons","TagEmSysTightElectrons");
		TagLooseElectrons tagLooseEle("TagEmSysLooseElectrons","TagEmSysLooseElectrons");;

		int iEleTightId = tagTightEle.GetPhoLikeEleIdWord(sp.GetPhoton(), event, fTrackBlock, zvx);
		int iEleLooseId = tagLooseEle.GetPhoLikeEleIdWord(sp.GetPhoton(), event, fTrackBlock, zvx);
		sp.SetTightElectronId(iEleTightId);
		sp.SetLooseElectronId(iEleLooseId);


		//now get the Std. electron IDs
		int iStdEleLooseId = -1;
		TStnElectron* Ele = TPhotonUtil::MatchElectron(event,sp.GetPhoton()); // get matching electron

		if (sp.GetDetector() == 0)
		{
			iStdEleLooseId = tagLooseEle.GetStdCenLooseEleID(Ele, fTrackBlock,trigMod->GetN12vertex());	
			sp.SetStdLooseElectronId(iStdEleLooseId);
		} else if (sp.GetDetector() == 1)
		{
			iStdEleLooseId = tagLooseEle.GetStdPlugLooseEleID(Ele, event,trigMod->GetN12vertex());	
			sp.SetStdLooseElectronId(iStdEleLooseId);
		}

		// To see the effect of std. loose electron id cuts and the photon-like electron id cuts
		//bcos, it seems when the photon-like electron id finds an electron it is not found by std. electron
		// id cuts


		// now copy the temp to real
		vSp[i] = sp;	
	} // for
	
} 

/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetTriggerInfo()
{
	stuple.tri_pho25iso 	= trigMod->GetTrig25IsoBit();
	stuple.tri_pho50 		= trigMod->GetTrig50Bit();
	stuple.tri_pho70 		= trigMod->GetTrig70Bit();
}
	
/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetVertexInfo()
{
	stuple.vtx_N = trigMod->GetNvertex();
	stuple.vtx_NClass12 = trigMod->GetN12vertex();
	stuple.vtx_z = trigMod->GetBestVertexZ();  //hisghest sum pt vertex	
	stuple.vtx_Ntracks = trigMod->GetBestVertexNTracks();
	stuple.vtx_SumPt = trigMod->GetBestVertexSumPt();
}

/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetMetInfo()
{
	stuple.met_Met 	= jetMod->GetMyMetCorr(0,3);
	stuple.met_RawMet	= jetMod->GetMyMetCorr(0,0);
	stuple.met_MetX 	= jetMod->GetMyMetXCorr(0,3);
	stuple.met_MetY 	= jetMod->GetMyMetYCorr(0,3);
	stuple.met_SumEt 	= jetMod->GetMySumEtCorr(0,3);
	stuple.met_Ht 		= jetMod->GetMyHtCorr(0,3);
	stuple.met_MetPhi = jetMod->GetMyMetPhiCorr(0,3);
	
	//not using MET MODEL anymore. commenting out
	//Sep 07,2009
	
	//TVector2 GenMet_d(jetMod->GetGenMet_Def());
	//TVector2 GenMet_m(jetMod->GetGenMet_devm());
	//TVector2 GenMet_p(jetMod->GetGenMet_devp());
	//TVector2 GenMet_mUn(jetMod->GetGenMet_devmUn());
	//TVector2 GenMet_pUn(jetMod->GetGenMet_devpUn());
}
	

/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetJetInfo()
{
	stuple.jet_num = jetMod->GetMyNjet(0,3);
	stuple.jet_NJet15 = jetMod->GetMyNjet(0,3);
	//stuple.jet_NJet15 = jetMod->GetMyNjet(0,0);


/*	if (stuple.jet_num > 1)
	{
		std::cout << "\t====== CENTRAL JETS =====" << std::endl;
		GetHeaderBlock()->Print();
	}
*/

	for (unsigned int i=0; i < stuple.jet_NJet15; i++) {
		stuple.jet_Index[i] = jetMod->GetMyJetBlockInd(0,i);
		TLorentzVector tlJetVec = *(jetMod->GetMyJetNoEMobj_lev6(0,i));	
		stuple.jet_Pt[i] = tlJetVec.Pt();	
		stuple.jet_E[i]  = tlJetVec.Energy();
		stuple.jet_Px[i] = tlJetVec.Px();
		stuple.jet_Py[i] = tlJetVec.Py();
		stuple.jet_Pz[i] = tlJetVec.Pz();
		stuple.jet_DetEta[i] = jetMod->GetMyJetEtaDetCorr(0,i);
		//stuple.jet_DetPhi[i] = jetMod->GetMy(0,i);
		stuple.jet_Emfr[i] = jetMod->GetMyJetEmFrCorr(0,i);
		if (stuple.jet_Emfr[i] !=0) stuple.jet_HadEm[i] = 1/stuple.jet_Emfr[i] - 1;
		stuple.jet_Ntowers[i] = jetMod->GetMyJetNtwr(0,i);
		stuple.jet_Ntracks[i] = jetMod->GetMyJetNtrk(0,i);

		TStnJet* jet = fJetBlock->Jet(stuple.jet_Index[i]);
		stuple.jet_SeedIPhi[i] = jet->SeedIEta();
		stuple.jet_SeedIEta[i] = jet->SeedIPhi();
		stuple.jet_EmTime[i] = GetEmTiming(stuple.jet_Index[i], stuple.jet_DetEta[i]);

		const float fSecVtx_Tag = jet->Tag(); //tag
		const float fSecVtx_ptag = jet->Scbpb(); //+ve tag
		const float fSecVtx_ntag = jet->Scnpb(); //-ve tag
		const float fSecVtx_Trkmass = jet->TrackMass();

		stuple.jet_SecVtxTag[i] = fSecVtx_Tag;
		stuple.jet_SecVtxppb[i] = fSecVtx_ptag;
		stuple.jet_SecVtxnpb[i] = fSecVtx_ntag;
		stuple.jet_SecVtxTrkmass[i] = fSecVtx_Trkmass;

		/*if (fSecVtx_Tag>0)
		{
		std::cout << "tag/ptag/ntag/mass=" << fSecVtx_Tag << "/ " << fSecVtx_ptag << "/ "
				<< fSecVtx_ntag << "/ " << fSecVtx_Trkmass << std::endl;
		}*/

				
		/*if (stuple.jet_num > 1)
		{
			std::cout << "jet i, Pt, E, eta, hadem=" << stuple.jet_Index[i] 
					<< "\t" << stuple.jet_Pt[i] << "\t" 
					<< stuple.jet_E[i] << "\t" << stuple.jet_DetEta[i] 
					<< "\t" << stuple.jet_HadEm[i]<< std::endl;
		}
		*/
		
		
	}

}


/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetJesUpJetInfo()
{
	stuple.jet_up_num = jetModUp->GetMyNjet(0,3);
	stuple.jet_up_NJet15 = jetModUp->GetMyNjet(0,3);
	
/*	if (stuple.jet_up_num > 1)
	{
		std::cout << "\t====== JESUP =====" << std::endl;
		GetHeaderBlock()->Print();
	}
*/

	for (unsigned int i=0; i < stuple.jet_up_NJet15; i++) {
		stuple.jet_up_Index[i] = jetModUp->GetMyJetBlockInd(0,i);
		stuple.jet_up_Pt[i] = jetModUp->GetMyJetNoEMobj_lev6(0,i)->Pt();	
		stuple.jet_up_E[i] = jetModUp->GetMyJetNoEMobj_lev6(0,i)->E();
		stuple.jet_up_Px[i] = jetModUp->GetMyJetNoEMobj_lev6(0,i)->Px();
		stuple.jet_up_Py[i] = jetModUp->GetMyJetNoEMobj_lev6(0,i)->Py();
		stuple.jet_up_Pz[i] = jetModUp->GetMyJetNoEMobj_lev6(0,i)->Pz();
		stuple.jet_up_DetEta[i] = jetModUp->GetMyJetEtaDetCorr(0,i);
		//stuple.jet_up_DetPhi[i] = jetModUp->GetMy(0,i);
		stuple.jet_up_Emfr[i] = jetModUp->GetMyJetEmFrCorr(0,i);
		if (stuple.jet_up_Emfr[i] !=0) stuple.jet_up_HadEm[i] = 1/stuple.jet_up_Emfr[i] - 1;
		stuple.jet_up_Ntowers[i] = jetModUp->GetMyJetNtwr(0,i);
		stuple.jet_up_Ntracks[i] = jetModUp->GetMyJetNtrk(0,i);

		TStnJet* jet = fJetBlock->Jet(stuple.jet_up_Index[i]);
		stuple.jet_up_SeedIPhi[i] = jet->SeedIEta();
		stuple.jet_up_SeedIEta[i] = jet->SeedIPhi();
		stuple.jet_up_EmTime[i] = GetEmTiming(stuple.jet_up_Index[i], stuple.jet_up_DetEta[i]);

		stuple.jet_up_SecVtxTag[i] = jet->Tag();
		stuple.jet_up_SecVtxppb[i] = jet->Scbpb();
		stuple.jet_up_SecVtxnpb[i] = jet->Scnpb();
		stuple.jet_up_SecVtxTrkmass[i] = jet->TrackMass();

/*		if (stuple.jet_up_num > 1)
		{
			std::cout << "jet up i, Pt, E, Eta, HadEm = " << stuple.jet_up_Index[i] 
					<< "\t" << stuple.jet_up_Pt[i] << "\t" 
					<< stuple.jet_up_E[i] << "\t" << stuple.jet_up_DetEta[i]
					<< "\t" << stuple.jet_up_HadEm[i] << std::endl;
		}
*/
		
	}

}


/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetJesDownJetInfo()
{
	stuple.jet_down_num = jetModDown->GetMyNjet(0,3);
	stuple.jet_down_NJet15 = jetModDown->GetMyNjet(0,3);
	
/*	if (stuple.jet_down_num > 1)
	{
		std::cout << "\t====== JESDOWN =====" << std::endl;
		GetHeaderBlock()->Print();
	}
*/

	for (unsigned int i=0; i < stuple.jet_down_NJet15; i++) {
		stuple.jet_down_Index[i] = jetModDown->GetMyJetBlockInd(0,i);
		stuple.jet_down_Pt[i] = jetModDown->GetMyJetNoEMobj_lev6(0,i)->Pt();	
		stuple.jet_down_E[i] = jetModDown->GetMyJetNoEMobj_lev6(0,i)->E();
		stuple.jet_down_Px[i] = jetModDown->GetMyJetNoEMobj_lev6(0,i)->Px();
		stuple.jet_down_Py[i] = jetModDown->GetMyJetNoEMobj_lev6(0,i)->Py();
		stuple.jet_down_Pz[i] = jetModDown->GetMyJetNoEMobj_lev6(0,i)->Pz();
		stuple.jet_down_DetEta[i] = jetModDown->GetMyJetEtaDetCorr(0,i);
		//stuple.jet_down_DetPhi[i] = jetModDown->GetMy(0,i);
		stuple.jet_down_Emfr[i] = jetModDown->GetMyJetEmFrCorr(0,i);
		if (stuple.jet_down_Emfr[i] !=0) stuple.jet_down_HadEm[i] = 1/stuple.jet_down_Emfr[i] - 1;
		stuple.jet_down_Ntowers[i] = jetModDown->GetMyJetNtwr(0,i);
		stuple.jet_down_Ntracks[i] = jetModDown->GetMyJetNtrk(0,i);

		TStnJet* jet = fJetBlock->Jet(stuple.jet_down_Index[i]);
		stuple.jet_down_SeedIPhi[i] = jet->SeedIEta();
		stuple.jet_down_SeedIEta[i] = jet->SeedIPhi();

		stuple.jet_down_EmTime[i] = GetEmTiming(stuple.jet_down_Index[i], stuple.jet_down_DetEta[i]);
		stuple.jet_down_SecVtxTag[i] = jet->Tag();
		stuple.jet_down_SecVtxppb[i] = jet->Scbpb();
		stuple.jet_down_SecVtxnpb[i] = jet->Scnpb();
		stuple.jet_down_SecVtxTrkmass[i] = jet->TrackMass();

		/*if (stuple.jet_down_num > 0)
		{
			std::cout << "jet up i, Pt, E, Eta, HadEm = " << stuple.jet_down_Index[i] 
					<< "\t" << stuple.jet_down_Pt[i] << "\t" 
					<< stuple.jet_down_E[i] << "\t" << stuple.jet_down_DetEta[i]
					<< "\t" << stuple.jet_down_HadEm[i] << std::endl;
		}
		*/
	
	}

}

/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetRawJetInfo()
{


/*	if (stuple.jet_num > 1)
	{
		std::cout << "\t====== CENTRAL JETS =====" << std::endl;
		GetHeaderBlock()->Print();
	}
*/

	for (unsigned int i=0; i < jetMod->GetMyNjet(0,0); i++) {
		if (jetMod->GetMyJetNoEMobj_raw(0,i)->Pt()<5.) continue;
		stuple.jet_raw_Index[stuple.jet_raw_num] = jetMod->GetMyJetBlockInd(0,i);
		stuple.jet_raw_Pt[stuple.jet_raw_num] = jetMod->GetMyJetNoEMobj_raw(0,i)->Pt();	
		stuple.jet_raw_E[stuple.jet_raw_num]  = jetMod->GetMyJetNoEMobj_raw(0,i)->Energy();
		stuple.jet_raw_Px[stuple.jet_raw_num] = jetMod->GetMyJetNoEMobj_raw(0,i)->Px();
		stuple.jet_raw_Py[stuple.jet_raw_num] = jetMod->GetMyJetNoEMobj_raw(0,i)->Py();
		stuple.jet_raw_Pz[stuple.jet_raw_num] = jetMod->GetMyJetNoEMobj_raw(0,i)->Pz();
		stuple.jet_raw_num++ ;
	}

}



/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetPhotonInfo(int iArrInd, SuperPhoton sp)
{
	stuple.pho_PhoBlockIndex[iArrInd] = sp.GetPhotonBlockIndex();
	stuple.pho_Index[iArrInd] 		= sp.GetPhotonBlockIndex();
	stuple.pho_Etc[iArrInd]   		= sp.GetEtCorr();
	stuple.pho_E[iArrInd]  			= sp.GetCorVec().E();
	stuple.pho_Px[iArrInd] 			= sp.GetCorVec().Px();		//do we need to get his from electron Block for electrons!
	stuple.pho_Py[iArrInd] 			= sp.GetCorVec().Py();
	stuple.pho_Pz[iArrInd] 			= sp.GetCorVec().Pz();
	stuple.pho_Detector[iArrInd] 	= sp.GetDetector();
	stuple.pho_DetEta[iArrInd] 	= sp.GetDetEta();
	stuple.pho_DetPhi[iArrInd] 	= sp.GetDetPhi();
	stuple.pho_XCes[iArrInd] 		= sp.GetXCes();
	stuple.pho_ZCes[iArrInd] 		= sp.GetZCes();
	stuple.pho_HadEm[iArrInd] 		= sp.GetHadEm();
	stuple.pho_Chi2Mean[iArrInd] 	= sp.GetChi2Mean();
	stuple.pho_N3d[iArrInd] 		= sp.GetN3d();
	stuple.pho_Iso4[iArrInd] 		= sp.GetIso4();
	stuple.pho_TrkPt[iArrInd] 		= sp.GetTrkPt();
	stuple.pho_TrkPt2[iArrInd] 	= sp.GetTrkPt2();
	stuple.pho_TrkIso[iArrInd] 	= sp.GetTrkIso();
	stuple.pho_CesWireE2[iArrInd] 	= sp.GetCesWireE2();
	stuple.pho_CesStripE2[iArrInd] 	= sp.GetCesStripE2();
	stuple.pho_PhiWedge[iArrInd] 		= sp.GetPhoton()->PhiSeedIndex();
	stuple.pho_NMuonStubs[iArrInd] 	= sp.GetNMuonStubs();
	stuple.pho_EmTime[iArrInd] 		= sp.GetEmTime();
	stuple.pho_TightId[iArrInd] 		= sp.GetTightPhotonId();
	stuple.pho_LooseId[iArrInd] 		= sp.GetLoosePhotonId();
	stuple.pho_PhoenixId[iArrInd] 	= sp.GetPhoenixId();
	
	//std::cout << cyan << "pho " << iArrInd << " " << sp.GetEtCorr() << "/" << sp.GetDetEta() << "\t" << sp.GetEmTime() << clearatt << std::endl;
	
	stuple.pho_matchJetIndex[iArrInd] = jetMod->GetMatch_JetIndex(0,sp.GetPhotonBlockIndex());//0=ObjType(photon) 
	if (stuple.pho_matchJetIndex[iArrInd] < 0) {
		StdOut(__FILE__,__LINE__,3,"pho did not match to a jet.");
		std::cout << "ERROR_PHO::"<< GetHeaderBlock()->RunNumber() << "," << GetHeaderBlock()->EventNumber() << std::endl;
		bJetEmMismatch = true;
		//exit (1);
	} 
		
	//stuff to get fake photon fraction from CES/CPR method
	TStnEvent* event = GetEvent();
	stuple.pho_CprWgt[iArrInd]    = cesCpr.AddPhoton(event, sp.GetPhoton());

	stuple.pho_CprSys1[iArrInd]  = cesCpr.CesCprSys(0);
	stuple.pho_CprSys2[iArrInd]  = cesCpr.CesCprSys(1);
	stuple.pho_CprSys3[iArrInd]  = cesCpr.CesCprSys(2);
	stuple.pho_CprSys4[iArrInd]  = cesCpr.CesCprSys(3);
	stuple.pho_CprSys5[iArrInd]  = cesCpr.CesCprSys(4);
	stuple.pho_CprSys6[iArrInd]  = cesCpr.CesCprSys(5);
	stuple.pho_CprSys7[iArrInd]  = cesCpr.CesCprSys(6);
	stuple.pho_CprSys8[iArrInd]  = cesCpr.CesCprSys(7);
}


/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetElectronInfo(int iArrInd, SuperPhoton sp)
{
	stuple.ele_EleBlockIndex[iArrInd] = sp.GetElectronBlockIndex();
	stuple.ele_Index[iArrInd] = sp.GetPhotonBlockIndex();
	stuple.ele_Etc[iArrInd]   		= sp.GetEtCorr();
	stuple.ele_E[iArrInd]  			= sp.GetCorVec().E();
	stuple.ele_Px[iArrInd] 			= sp.GetCorVec().Px();
	stuple.ele_Py[iArrInd] 			= sp.GetCorVec().Py();
	stuple.ele_Pz[iArrInd] 			= sp.GetCorVec().Pz();
	stuple.ele_Detector[iArrInd] 	= sp.GetDetector();
	stuple.ele_DetEta[iArrInd] 	= sp.GetDetEta();
	stuple.ele_DetPhi[iArrInd] 	= sp.GetDetPhi();
	stuple.ele_XCes[iArrInd] 		= sp.GetXCes();
	stuple.ele_ZCes[iArrInd] 		= sp.GetZCes();
	stuple.ele_HadEm[iArrInd] 		= sp.GetHadEm();
	stuple.ele_Chi2Mean[iArrInd] 	= sp.GetChi2Mean();
	stuple.ele_N3d[iArrInd] 		= sp.GetN3d();
	stuple.ele_Iso4[iArrInd] 		= sp.GetIso4();
	stuple.ele_TrkIso[iArrInd] 	= sp.GetTrkIso();
	stuple.ele_CesWireE2[iArrInd] 	= sp.GetCesWireE2();
	stuple.ele_CesStripE2[iArrInd] 	= sp.GetCesStripE2();
	stuple.ele_PhiWedge[iArrInd] 		= sp.GetPhoton()->PhiSeedIndex();
	stuple.ele_NMuonStubs[iArrInd] 	= sp.GetNMuonStubs();
	stuple.ele_EmTime[iArrInd] 		= sp.GetEmTime();
	stuple.ele_PhoenixId[iArrInd] 	= sp.GetPhoenixId();
	stuple.ele_TightId[iArrInd]      = sp.GetTightElectronId();
	stuple.ele_LooseId[iArrInd]      = sp.GetLooseElectronId();
	stuple.ele_StdTightId[iArrInd]   = sp.GetStdTightElectronId();
	stuple.ele_StdLooseId[iArrInd]   = sp.GetStdLooseElectronId();
	stuple.ele_ConversionId[iArrInd] = sp.GetConversionId();
	//if there is no 2nd track TStnPhoton::TrkPt2() returns value of the 1st track
	stuple.ele_TrackPt2[iArrInd]     = sp.GetTrkPt2(); //there is no simple way to get the 2nd track from electron block.
																		//so I'll use what is given in photon block
	

	TStnEvent* event = GetEvent();
	TStnElectron* Ele=TPhotonUtil::MatchElectron(event,sp.GetPhoton()); // get matching electron
	if (Ele) {
		if (Ele->TrackNumber() >=0) {
			stuple.ele_Ntracks[iArrInd]   = Ele->NTracks();
			stuple.ele_Emfr[iArrInd]      = Ele->Emfr();
			stuple.ele_EoverP[iArrInd]    = Ele->EOverP();
			stuple.ele_TrackPt[iArrInd]   = Ele->TrackPt();
			stuple.ele_TrackBcPt[iArrInd] = Ele->TrackBcPt();
			stuple.ele_TrackPhi[iArrInd]  = Ele->TrackPhi();
			stuple.ele_Nssl[iArrInd]      = (Int_t) Ele->Nssl();  //TStnElectron returns Float_t??
			stuple.ele_Nasl[iArrInd]      = (Int_t) Ele->Nasl();
			/*std::cout << __LINE__ << "Ele trk pt from ele, from pho = " << Ele->TrackPt() << "/ " << sp.GetTrkPt() << " = " << Ele->TrackPt()/sp.GetTrkPt() 
				 << " (ntrks = " << Ele->NTracks() << ", sp.n3d = " << sp.GetN3d() << std::endl;
			if (Ele->NTracks()>1)
			{
				 std::cout << "\t" << __LINE__ << "Ele trk pt2 from pho block   = " << sp.GetTrkPt2() << std::endl;
			} else if (Ele->NTracks() <=1) 
			{
				 std::cout << __LINE__ << "Ele trk pt2 from pho block   = " << sp.GetTrkPt()  << " :::: this must <0 as no 2nd trk is found !" << std::endl;
			}
			*/
		}
	}

	stuple.ele_matchJetIndex[iArrInd] = jetMod->GetMatch_JetIndex(1,sp.GetElectronBlockIndex());//1=ObjType(electron) 
	if (stuple.ele_matchJetIndex[iArrInd] < 0) {
		StdOut(__FILE__,__LINE__,3," ele did not match to a jet.");
		std::cout << "ERROR_ELE::"<< GetHeaderBlock()->RunNumber() << "," << GetHeaderBlock()->RunNumber() << std::endl;
		bJetEmMismatch = true;
		//exit (1);
	} 
															 
}


/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetEmUpPhotonInfo(int iArrInd, SuperPhoton sp)
{
	stuple.pho_up_PhoBlockIndex[iArrInd] = sp.GetPhotonBlockIndex();
	stuple.pho_up_Index[iArrInd] 		= sp.GetPhotonBlockIndex();
	stuple.pho_up_Etc[iArrInd]   		= sp.GetEtCorr();
	stuple.pho_up_E[iArrInd]  			= sp.GetCorVec().E();
	stuple.pho_up_Px[iArrInd] 			= sp.GetCorVec().Px();		//do we need to get his from electron Block for electrons!
	stuple.pho_up_Py[iArrInd] 			= sp.GetCorVec().Py();
	stuple.pho_up_Pz[iArrInd] 			= sp.GetCorVec().Pz();
	stuple.pho_up_Detector[iArrInd] 	= sp.GetDetector();
	stuple.pho_up_DetEta[iArrInd] 	= sp.GetDetEta();
	stuple.pho_up_DetPhi[iArrInd] 	= sp.GetDetPhi();
	stuple.pho_up_XCes[iArrInd] 		= sp.GetXCes();
	stuple.pho_up_ZCes[iArrInd] 		= sp.GetZCes();
	stuple.pho_up_HadEm[iArrInd] 		= sp.GetHadEm();
	stuple.pho_up_Chi2Mean[iArrInd] 	= sp.GetChi2Mean();
	stuple.pho_up_N3d[iArrInd] 		= sp.GetN3d();
	stuple.pho_up_Iso4[iArrInd] 		= sp.GetIso4();
	stuple.pho_up_TrkPt[iArrInd] 		= sp.GetTrkPt();
	stuple.pho_up_TrkPt2[iArrInd] 	= sp.GetTrkPt2();
	stuple.pho_up_TrkIso[iArrInd] 	= sp.GetTrkIso();
	stuple.pho_up_CesWireE2[iArrInd] 	= sp.GetCesWireE2();
	stuple.pho_up_CesStripE2[iArrInd] 	= sp.GetCesStripE2();
	stuple.pho_up_PhiWedge[iArrInd] 		= sp.GetPhoton()->PhiSeedIndex();
	stuple.pho_up_NMuonStubs[iArrInd] 	= sp.GetNMuonStubs();
	stuple.pho_up_EmTime[iArrInd] 		= sp.GetEmTime();
	stuple.pho_up_TightId[iArrInd] 		= sp.GetTightPhotonId();
	stuple.pho_up_LooseId[iArrInd] 		= sp.GetLoosePhotonId();
	stuple.pho_up_PhoenixId[iArrInd] 	= sp.GetPhoenixId();

	stuple.pho_up_matchJetIndex[iArrInd] = jetMod->GetMatch_JetIndex(0,sp.GetPhotonBlockIndex());//0=ObjType(photon) 
	if (stuple.pho_up_matchJetIndex[iArrInd] < 0) {
		StdOut(__FILE__,__LINE__,3,"pho did not match to a jet.");
		std::cout << "ERROR_PHO::"<< GetHeaderBlock()->RunNumber() << "," << GetHeaderBlock()->EventNumber() << std::endl;
		bJetEmMismatch = true;
		//exit (1);
	} 

	TStnEvent* event = GetEvent();
	stuple.pho_up_CprWgt[iArrInd]    = cesCprUp.AddPhoton(event, sp.GetPhoton());
	stuple.pho_up_CprSys1[iArrInd]  = cesCprUp.CesCprSys(0);
	stuple.pho_up_CprSys2[iArrInd]  = cesCprUp.CesCprSys(1);
	stuple.pho_up_CprSys3[iArrInd]  = cesCprUp.CesCprSys(2);
	stuple.pho_up_CprSys4[iArrInd]  = cesCprUp.CesCprSys(3);
	stuple.pho_up_CprSys5[iArrInd]  = cesCprUp.CesCprSys(4);
	stuple.pho_up_CprSys6[iArrInd]  = cesCprUp.CesCprSys(5);
	stuple.pho_up_CprSys7[iArrInd]  = cesCprUp.CesCprSys(6);
	stuple.pho_up_CprSys8[iArrInd]  = cesCprUp.CesCprSys(7);
	
}


/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetEmUpElectronInfo(int iArrInd, SuperPhoton sp)
{
	stuple.ele_up_EleBlockIndex[iArrInd] = sp.GetElectronBlockIndex();
	stuple.ele_up_Index[iArrInd] 		= sp.GetPhotonBlockIndex();
	stuple.ele_up_Etc[iArrInd]   		= sp.GetEtCorr();
	stuple.ele_up_E[iArrInd]  			= sp.GetCorVec().E();
	stuple.ele_up_Px[iArrInd] 			= sp.GetCorVec().Px();
	stuple.ele_up_Py[iArrInd] 			= sp.GetCorVec().Py();
	stuple.ele_up_Pz[iArrInd] 			= sp.GetCorVec().Pz();
	stuple.ele_up_Detector[iArrInd] 	= sp.GetDetector();
	stuple.ele_up_DetEta[iArrInd] 	= sp.GetDetEta();
	stuple.ele_up_DetPhi[iArrInd] 	= sp.GetDetPhi();
	stuple.ele_up_XCes[iArrInd] 		= sp.GetXCes();
	stuple.ele_up_ZCes[iArrInd] 		= sp.GetZCes();
	stuple.ele_up_HadEm[iArrInd] 		= sp.GetHadEm();
	stuple.ele_up_Chi2Mean[iArrInd] 	= sp.GetChi2Mean();
	stuple.ele_up_N3d[iArrInd] 		= sp.GetN3d();
	stuple.ele_up_Iso4[iArrInd] 		= sp.GetIso4();
	stuple.ele_up_TrkIso[iArrInd] 	= sp.GetTrkIso();
	stuple.ele_up_CesWireE2[iArrInd] 	= sp.GetCesWireE2();
	stuple.ele_up_CesStripE2[iArrInd] 	= sp.GetCesStripE2();
	stuple.ele_up_PhiWedge[iArrInd] 		= sp.GetPhoton()->PhiSeedIndex();
	stuple.ele_up_NMuonStubs[iArrInd] 	= sp.GetNMuonStubs();
	stuple.ele_up_EmTime[iArrInd] 		= sp.GetEmTime();
	stuple.ele_up_PhoenixId[iArrInd] 	= sp.GetPhoenixId();
	stuple.ele_up_TightId[iArrInd]      = sp.GetTightElectronId();
	stuple.ele_up_LooseId[iArrInd]      = sp.GetLooseElectronId();
	stuple.ele_up_StdTightId[iArrInd]   = sp.GetStdTightElectronId();
	stuple.ele_up_StdLooseId[iArrInd]   = sp.GetStdLooseElectronId();
	stuple.ele_up_ConversionId[iArrInd] = sp.GetConversionId();
	stuple.ele_up_TrackPt2[iArrInd]     = sp.GetTrkPt2(); //there is no simple way to get the 2nd track from electron block.
																			//so I'll use what is given in photon block
	
	TStnEvent* event = GetEvent();
	TStnElectron* Ele=TPhotonUtil::MatchElectron(event,sp.GetPhoton()); // get matching electron
	if (Ele) {
		if (Ele->TrackNumber() >=0) {
			stuple.ele_up_Ntracks[iArrInd]   = Ele->NTracks();
			stuple.ele_up_Emfr[iArrInd]      = Ele->Emfr();
			stuple.ele_up_EoverP[iArrInd]    = Ele->EOverP();
			stuple.ele_up_TrackPt[iArrInd]   = Ele->TrackPt();
			stuple.ele_up_TrackBcPt[iArrInd] = Ele->TrackBcPt();
			stuple.ele_up_TrackPhi[iArrInd]  = Ele->TrackPhi();
			stuple.ele_up_Nssl[iArrInd]      = (Int_t) Ele->Nssl();  //TStnElectron returns Float_t??
			stuple.ele_up_Nasl[iArrInd]      = (Int_t) Ele->Nasl();
		}
	}

	stuple.ele_up_matchJetIndex[iArrInd] = jetMod->GetMatch_JetIndex(1,sp.GetElectronBlockIndex());//1=ObjType(electron) 
	if (stuple.ele_up_matchJetIndex[iArrInd] < 0) {
		StdOut(__FILE__,__LINE__,3," ele did not match to a jet.");
		std::cout << "ERROR_ELE::"<< GetHeaderBlock()->RunNumber() << "," << GetHeaderBlock()->RunNumber() << std::endl;
		bJetEmMismatch = true;
		//exit (1);
	} 
															 
}


/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetEmDownPhotonInfo(int iArrInd, SuperPhoton sp)
{
	stuple.pho_down_PhoBlockIndex[iArrInd] = sp.GetPhotonBlockIndex();
	stuple.pho_down_Index[iArrInd] 		= sp.GetPhotonBlockIndex();
	stuple.pho_down_Etc[iArrInd]   		= sp.GetEtCorr();
	stuple.pho_down_E[iArrInd]  			= sp.GetCorVec().E();
	stuple.pho_down_Px[iArrInd] 			= sp.GetCorVec().Px();		//do we need to get his from electron Block for electrons!
	stuple.pho_down_Py[iArrInd] 			= sp.GetCorVec().Py();
	stuple.pho_down_Pz[iArrInd] 			= sp.GetCorVec().Pz();
	stuple.pho_down_Detector[iArrInd] 	= sp.GetDetector();
	stuple.pho_down_DetEta[iArrInd] 	= sp.GetDetEta();
	stuple.pho_down_DetPhi[iArrInd] 	= sp.GetDetPhi();
	stuple.pho_down_XCes[iArrInd] 		= sp.GetXCes();
	stuple.pho_down_ZCes[iArrInd] 		= sp.GetZCes();
	stuple.pho_down_HadEm[iArrInd] 		= sp.GetHadEm();
	stuple.pho_down_Chi2Mean[iArrInd] 	= sp.GetChi2Mean();
	stuple.pho_down_N3d[iArrInd] 		= sp.GetN3d();
	stuple.pho_down_Iso4[iArrInd] 		= sp.GetIso4();
	stuple.pho_down_TrkPt[iArrInd] 		= sp.GetTrkPt();
	stuple.pho_down_TrkPt2[iArrInd] 		= sp.GetTrkPt2();
	stuple.pho_down_TrkIso[iArrInd] 	= sp.GetTrkIso();
	stuple.pho_down_CesWireE2[iArrInd] 	= sp.GetCesWireE2();
	stuple.pho_down_CesStripE2[iArrInd] 	= sp.GetCesStripE2();
	stuple.pho_down_PhiWedge[iArrInd] 		= sp.GetPhoton()->PhiSeedIndex();
	stuple.pho_down_NMuonStubs[iArrInd] 	= sp.GetNMuonStubs();
	stuple.pho_down_EmTime[iArrInd] 		= sp.GetEmTime();
	stuple.pho_down_TightId[iArrInd] 		= sp.GetTightPhotonId();
	stuple.pho_down_LooseId[iArrInd] 		= sp.GetLoosePhotonId();
	stuple.pho_down_PhoenixId[iArrInd] 	= sp.GetPhoenixId();

	stuple.pho_down_matchJetIndex[iArrInd] = jetMod->GetMatch_JetIndex(0,sp.GetPhotonBlockIndex());//0=ObjType(photon) 
	if (stuple.pho_down_matchJetIndex[iArrInd] < 0) {
		StdOut(__FILE__,__LINE__,3,"pho did not match to a jet.");
		std::cout << "ERROR_PHO::"<< GetHeaderBlock()->RunNumber() << "," << GetHeaderBlock()->EventNumber() << std::endl;
		bJetEmMismatch = true;
		//exit (1);
	} 

	TStnEvent* event = GetEvent();
	stuple.pho_down_CprWgt[iArrInd]   = cesCprDown.AddPhoton(event, sp.GetPhoton());
	stuple.pho_down_CprSys1[iArrInd]  = cesCprDown.CesCprSys(0);
	stuple.pho_down_CprSys2[iArrInd]  = cesCprDown.CesCprSys(1);
	stuple.pho_down_CprSys3[iArrInd]  = cesCprDown.CesCprSys(2);
	stuple.pho_down_CprSys4[iArrInd]  = cesCprDown.CesCprSys(3);
	stuple.pho_down_CprSys5[iArrInd]  = cesCprDown.CesCprSys(4);
	stuple.pho_down_CprSys6[iArrInd]  = cesCprDown.CesCprSys(5);
	stuple.pho_down_CprSys7[iArrInd]  = cesCprDown.CesCprSys(6);
	stuple.pho_down_CprSys8[iArrInd]  = cesCprDown.CesCprSys(7);
}


/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetEmDownElectronInfo(int iArrInd, SuperPhoton sp)
{
	stuple.ele_down_EleBlockIndex[iArrInd] = sp.GetElectronBlockIndex();
	stuple.ele_down_Index[iArrInd] 		= sp.GetPhotonBlockIndex();
	stuple.ele_down_Etc[iArrInd]   		= sp.GetEtCorr();
	stuple.ele_down_E[iArrInd]  			= sp.GetCorVec().E();
	stuple.ele_down_Px[iArrInd] 			= sp.GetCorVec().Px();
	stuple.ele_down_Py[iArrInd] 			= sp.GetCorVec().Py();
	stuple.ele_down_Pz[iArrInd] 			= sp.GetCorVec().Pz();
	stuple.ele_down_Detector[iArrInd] 	= sp.GetDetector();
	stuple.ele_down_DetEta[iArrInd] 	   = sp.GetDetEta();
	stuple.ele_down_DetPhi[iArrInd] 	   = sp.GetDetPhi();
	stuple.ele_down_XCes[iArrInd] 		= sp.GetXCes();
	stuple.ele_down_ZCes[iArrInd] 		= sp.GetZCes();
	stuple.ele_down_HadEm[iArrInd] 		= sp.GetHadEm();
	stuple.ele_down_Chi2Mean[iArrInd] 	= sp.GetChi2Mean();
	stuple.ele_down_N3d[iArrInd] 		   = sp.GetN3d();
	stuple.ele_down_Iso4[iArrInd] 		= sp.GetIso4();
	stuple.ele_down_TrkIso[iArrInd] 	   = sp.GetTrkIso();
	stuple.ele_down_CesWireE2[iArrInd] 	= sp.GetCesWireE2();
	stuple.ele_down_CesStripE2[iArrInd]   = sp.GetCesStripE2();
	stuple.ele_down_PhiWedge[iArrInd] 	  = sp.GetPhoton()->PhiSeedIndex();
	stuple.ele_down_NMuonStubs[iArrInd]   = sp.GetNMuonStubs();
	stuple.ele_down_EmTime[iArrInd] 		  = sp.GetEmTime();
	stuple.ele_down_PhoenixId[iArrInd] 	  = sp.GetPhoenixId();
	stuple.ele_down_TightId[iArrInd]      = sp.GetTightElectronId();
	stuple.ele_down_LooseId[iArrInd]      = sp.GetLooseElectronId();
	stuple.ele_down_StdTightId[iArrInd]   = sp.GetStdTightElectronId();
	stuple.ele_down_StdLooseId[iArrInd]   = sp.GetStdLooseElectronId();
	stuple.ele_down_ConversionId[iArrInd] = sp.GetConversionId();
	stuple.ele_down_TrackPt2[iArrInd]     = sp.GetTrkPt2(); //there is no simple way to get the 2nd track from electron block.
																				//so I'll use what is given in photon block

	TStnEvent* event = GetEvent();
	TStnElectron* Ele=TPhotonUtil::MatchElectron(event,sp.GetPhoton()); // get matching electron
	if (Ele) {
		if (Ele->TrackNumber() >=0) {
			stuple.ele_down_Ntracks[iArrInd]   = Ele->NTracks();
			stuple.ele_down_Emfr[iArrInd]      = Ele->Emfr();
			stuple.ele_down_EoverP[iArrInd]    = Ele->EOverP();
			stuple.ele_down_TrackPt[iArrInd]   = Ele->TrackPt();
			stuple.ele_down_TrackBcPt[iArrInd] = Ele->TrackBcPt();
			stuple.ele_down_TrackPhi[iArrInd]  = Ele->TrackPhi();
			stuple.ele_down_Nssl[iArrInd]      = (Int_t) Ele->Nssl();  //TStnElectron returns Float_t??
			stuple.ele_down_Nasl[iArrInd]      = (Int_t) Ele->Nasl();
		}
	}

	stuple.ele_down_matchJetIndex[iArrInd] = jetMod->GetMatch_JetIndex(1,sp.GetElectronBlockIndex());//1=ObjType(electron) 
	if (stuple.ele_down_matchJetIndex[iArrInd] < 0) {
		StdOut(__FILE__,__LINE__,3," ele did not match to a jet.");
		std::cout << "ERROR_ELE::"<< GetHeaderBlock()->RunNumber() << "," << GetHeaderBlock()->RunNumber() << std::endl;
		bJetEmMismatch = true;
		//exit (1);
	} 
															 
}



/*-------------------------------------------------------------------*/
void FlatStupleMaker::GetGenLevelInfo()
{
	int Nparticles = fGenpBlock->NParticles();
	//int mom_id = 0;
	TLorentzVector helevec(0,0,0,0), pv(0,0,0,0);	// temp ele vec and production vertex vector for fid cut
	std::vector<TLorentzVector> elevec;					// two electron 4-vectors
	std::vector<TLorentzVector> elepv;					// two electron productions vectors
	//int ndau =0;	// number of daughter from Zee = 2

	TGenParticle *par, *mom;
	
	for (int i = 0 ; i < Nparticles ; i++) {	//_____ not much diff from running over all. my match is mostly within the fist 10 objects
		par = fGenpBlock->Particle(i);
		int im = par->GetFirstMother();

		if (im < 0) continue; 	//_________________________________im=-1 means no mother for imcoming particles
		mom = fGenpBlock->Particle(im);

		if (mom == 0) continue;
		
		int par_id   = par->GetPdgCode();
		int par_stat = par->GetStatusCode();
				
		if ( (abs(par_id) == 23 || abs(par_id) == 24) && (par_stat == 2)) { 	//found a Z/W 
			TLorentzVector vec;
			par->Momentum(vec);
			stuple.gen_MomIndex 	= i;
			stuple.gen_MomPDG 	= par_id;
			stuple.gen_MomStatus	= par_stat;
			stuple.gen_MomEtc		= vec.E() * TMath::Sin(vec.Theta());
			stuple.gen_MomE		= vec.E();
			stuple.gen_MomPx		= vec.Px();
			stuple.gen_MomPy		= vec.Py();
			stuple.gen_MomPz		= vec.Pz();
			
			TLorentzVector pv(0,0,0,0);
			par->ProductionVertex(pv);
			stuple.gen_ProdVtxX	= pv.X();
			stuple.gen_ProdVtxY	= pv.Y();
			stuple.gen_ProdVtxZ	= pv.Z();
			stuple.gen_ProdVtxT	= pv.T();
		}
					
		if (par->IsPhoton()) { 	//a photon Et>10GeV
			TLorentzVector vec;
			par->Momentum(vec);
			float Etc		= vec.E() * TMath::Sin(vec.Theta());

			if (Etc>10) {
				stuple.gen_pho_Index[stuple.gen_phonum]	= i;
				stuple.gen_pho_PDG[stuple.gen_phonum] 		= par_id;
				stuple.gen_pho_Status[stuple.gen_phonum]	= par_stat;
				stuple.gen_pho_Etc[stuple.gen_phonum]		= vec.E() * TMath::Sin(vec.Theta());
				stuple.gen_pho_E[stuple.gen_phonum]			= vec.E();
				stuple.gen_pho_Px[stuple.gen_phonum]		= vec.Px();
				stuple.gen_pho_Py[stuple.gen_phonum]		= vec.Py();
				stuple.gen_pho_Pz[stuple.gen_phonum]		= vec.Pz();
				
				TLorentzVector pv(0,0,0,0);
				par->ProductionVertex(pv);
				stuple.gen_pho_ProdVtxX[stuple.gen_phonum]	= pv.X();
				stuple.gen_pho_ProdVtxY[stuple.gen_phonum]	= pv.Y();
				stuple.gen_pho_ProdVtxZ[stuple.gen_phonum]	= pv.Z();
				stuple.gen_pho_ProdVtxT[stuple.gen_phonum]	= pv.T();

				if (stuple.gen_phonum == (unsigned) stuple.Ngenpho) {
					StdOut(__FILE__,__LINE__,3,"Stuple Gen Pho array size exceeded!");
					exit (1);
				}
				//must increment this at the end. otherwise i wont have anything for the 1st element in the array!	
				stuple.gen_phonum++;
			}
		}

		if (fabs(par->GetPdgCode()) == 11
			 || fabs(par->GetPdgCode()) == 13
			 || fabs(par->GetPdgCode()) == 15) { 	//e,nu,tau
			
			TLorentzVector vec;
			par->Momentum(vec);

			stuple.gen_ele_Index[stuple.gen_elenum]	= i;
			stuple.gen_ele_PDG[stuple.gen_elenum] 		= par_id;
			stuple.gen_ele_Status[stuple.gen_elenum]	= par_stat;
			stuple.gen_ele_Etc[stuple.gen_elenum]		= vec.E() * TMath::Sin(vec.Theta());
			stuple.gen_ele_E[stuple.gen_elenum]			= vec.E();
			stuple.gen_ele_Px[stuple.gen_elenum]		= vec.Px();
			stuple.gen_ele_Py[stuple.gen_elenum]		= vec.Py();
			stuple.gen_ele_Pz[stuple.gen_elenum]		= vec.Pz();
			
			TLorentzVector pv(0,0,0,0);
			par->ProductionVertex(pv);
			stuple.gen_ele_ProdVtxX[stuple.gen_elenum]	= pv.X();
			stuple.gen_ele_ProdVtxY[stuple.gen_elenum]	= pv.Y();
			stuple.gen_ele_ProdVtxZ[stuple.gen_elenum]	= pv.Z();
			stuple.gen_ele_ProdVtxT[stuple.gen_elenum]	= pv.T();

			if (stuple.gen_elenum == (unsigned) stuple.Ngenele) {
				StdOut(__FILE__,__LINE__,3,"Stuple Gen Ele array size exceeded!");
				exit (1);
			}

			//must increment this at the end. otherwise i wont have anything for the 1st element in the array!	
			stuple.gen_elenum++;
			
		}
					
	} //for


}


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void FlatStupleMaker::Cleanup()
{

} //Cleanup

/*-------------------------------------------------------------------*/
void FlatStupleMaker::FillDataBlocks(int ientry)
{
	fElectronBlock->GetEntry(ientry);
	fJetBlock->GetEntry(ientry);		// need to get Jet tower IEta/IPhi
	fVertexBlock->GetEntry(ientry);
	fTrackBlock->GetEntry(ientry);
	fSecVtxTagBlock->GetEntry(ientry);

	if (GetHeaderBlock()->McFlag()) {
		fGenpBlock->GetEntry(ientry);
	}

	fCalBlock->GetEntry(ientry);
}
//____________________________________________________________________
//  END JOB SUMMARY
//____________________________________________________________________
int FlatStupleMaker::EndJob() {

	if (rootFile->IsOpen()) {
		rootFile->Write();
		rootFile->Close();
	}

	if (GetSummaryStat()) return 0;

	float nTrue = cesCpr.Signal(false); // number of true pho passing cleanup
	float nTrueTot = cesCpr.Signal();      // number of true pho in the whole sample

	float nstaterrPass = cesCpr.SignalErr(false); // statistical uncertainty for events passing cleanup
	float	nstaterrTot = cesCpr.SignalErr(); // statistical uncertainty for events in the whole sample
		  


	printf("[FST:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[FST:01:] Events Processed ----------- = " << counter.evtsProcessed << std::endl;
	std::cout << "[FST:02:] Events Pass this module ---- = " << counter.evtsPassModule << std::endl;
	std::cout << "[FST:03:] MC flag  ------------------- = " << GetHeaderBlock()->McFlag() << std::endl;
	std::cout << "[FST:04:] EM Et threshold ------------ = " << GetEmEtThr() << std::endl;
	std::cout << "[FST:05:] Photon only evts ----------- = " << counter.phoOnlyEvts << std::endl;
	std::cout << "[FST:06:] Electron only evts --------- = " << counter.eleOnlyEvts << std::endl;
	std::cout << "[FST:07:] Photon and Electron evts --- = " << counter.phoEleEvts << std::endl;
	std::cout << "[FST:08:] CES/CPR N true photons ----- = " << nTrue << "+/-" << nstaterrPass << std::endl;
	std::cout << "[FST:10:] CES/CPR N true photons Total = " << nTrueTot << "+/-" << nstaterrTot << std::endl;
	std::cout << "[FST:11:] Stuple written to ---------- = " << sFileName << std::endl;
	printf("---------------------------------------------------\n");
	return 0;
}



/*-------------------------------------------------------------------*/
// runN>=190851  first run for EM timing system
// I talked with Lee today and he confirmed that CHA/WHA/PHA has no meaning 
// at all. Only CEM/PEM has EM timinig info. We choose the timing of the
// highest Et tower.
// _emTimeTower = myEmTiming->TEmTimingModule::getTimingTower(0,iEta,iPhi); // 0=CEM; 1=PEM; 2=CHA; 3=WHA; 4=PHA 
// So for jets, I need to choose CEM/PEM option based on the detector
// eta of the jet! 03-25-2010
// PEM covers only up eta<3.4
/*-------------------------------------------------------------------*/
double FlatStupleMaker::GetEmTiming(const int jetInd, const float jetDetEta)
{
	double time = -999999.99;

	if (jetInd < 0) {
		StdOut(__FILE__,__LINE__,3,"Invalid jet index!.");
		return time;
	}
	
	int runN = GetHeaderBlock()->RunNumber();
	if (runN >= 190851 && GetHeaderBlock()->McFlag() == false) { // first run for EM timing system
		TEmTimingModule* myEmTiming = (TEmTimingModule*) ((TStnAna*) GetAna()->GetModule("EmTimingAna"));
		
      if (myEmTiming == NULL) {
			StdOut(__FILE__,__LINE__,3," EmTiminigModule not found.");
			return time;
		} else {		
      	//__________________ timing for photons
			int Nhits=-1;
			int NGoodHits=-1;
			EmTimingTower* _emTimeTower = NULL;
			CalDataArray calholder;
			
			MatchCalorTowers(jetInd, &calholder);
			//std::cout << __FUNCTION__ << ": jet towers (jet ind[" << jetInd << "] ntwrs[" << calholder.size() << "]"  << std::endl;
			//std::cout << " i \t eta \t phi " << std::endl;
				
			int iCEMPEM = 0; // 0=CEM; 1=PEM;
			if (fabs(jetDetEta) > 1.1) iCEMPEM = 1;
			
	      for (unsigned int icount = 0 ; icount < calholder.size(); icount++) {
				TCalTower* ctower = calholder[icount];
				int iEta = ctower->IEta();
				int iPhi = ctower->IPhi();
				//std::cout << icount  << "\teta[" << iEta << "]\tphi[" << iPhi << "]\tEt["<< ctower->Et() << "]" << std::endl;
		  
				_emTimeTower = myEmTiming->TEmTimingModule::getTimingTower(iCEMPEM,iEta,iPhi); // 0=CEM; 1=PEM; 2=CHA; 3=WHA; 4=PHA 
		
				if (TETTYP[iEta] == 5 && iPhi%2 != 0) _emTimeTower = myEmTiming->TEmTimingModule::getTimingTower(0,iEta,iPhi-1);
					
				if (_emTimeTower != NULL) {
					Nhits=_emTimeTower->nHits();
						
					for (int i=0; i< Nhits; i++) {
						int status = _emTimeTower->getStatus(i);
						  
						if (TStnEmTimingBlock::EnergyTooLow(status)) continue;
						if (TStnEmTimingBlock::SigmasFromThreshold(status) < 3) continue;
						if (_emTimeTower->getT0(i)<-80.0 || _emTimeTower->getT0(i)>160.0) continue;
						
						NGoodHits++;
						time=_emTimeTower->getT0(i);
						break;
		   	 	}
				}
			}
		}

	}
  return time;
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void FlatStupleMaker::MatchCalorTowers(int jet_ind, CalDataArray* cdholder)
{
  cdholder->clear();
  TStnLinkBlock* links = fJetBlock->TowerLinkList();
  int nptow = links->NLinks(jet_ind);

  TCalTower* ctower = NULL;
  for(int j=0; j<nptow; j++) 
    {
      int iptow = links->Index(jet_ind,j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;

      ctower = fCalBlock->Tower(ieta,iphi);
      cdholder->push_back(ctower);
    }

  std::sort(cdholder->begin(), cdholder->end(), SortTowersByEnergy);

  //temporary : use only highest et tower
 /* if (nptow>1)
  {
  for (int i=nptow; i>1 ;--i) cdholder->pop_back();
	}
  std::cout << "energy = " << cdholder->at(0)->Et() << std::endl;
  */
  return;
}
