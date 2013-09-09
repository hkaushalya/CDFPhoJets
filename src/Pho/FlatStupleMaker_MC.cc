////////////////////////////////////////////////////////////
// This will create a flat stuples for me.                //
// this is only good for Photon MC as I try to match the  //
// detector photon to HEPG photon.                        //
////////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "samantha/Pho/FlatStupleMaker_MC.hh"
#include "Stntuple/loop/TStnAna.hh"
#include <iostream>
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/obj/TStnElectron.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "Stntuple/obj/TGenParticle.hh"
#include "TMath.h"
#include <assert.h>

ClassImp(FlatStupleMaker_MC)

//_____________________________________________________________________________
FlatStupleMaker_MC::FlatStupleMaker_MC(const char* name, const char* title):
  TStnModule(name,title),
  sFileName("Stuple.root"),
  bRunPermit(true)
{
	std::cout << "Hello I am FlatStupleMaker_MC module" << std::endl;
}

//_____________________________________________________________________________
FlatStupleMaker_MC::~FlatStupleMaker_MC() {
}

//_____________________________________________________________________________
void FlatStupleMaker_MC::SaveHistograms() {
}

//_____________________________________________________________________________
void FlatStupleMaker_MC::BookHistograms()
{
	TFolder *new_folder;

	new_folder = GetHistFolder(this, "DelR","DelR");

	hDelR = new TH1F("DelR","DelR(detPho,genPho)",100000,0,2);
	new_folder->Add(hDelR);
}


//_____________________________________________________________________________
int FlatStupleMaker_MC::BeginJob()
{

	rootFile = new TFile(sFileName.c_str(),"RECREATE");
	tree = new TTree("Stuple","Sam's flat Stuple.");

		tree->Branch("evt_McFlag"   , &stuple.evt_McFlag   , "evt_McFlag/I");
		tree->Branch("evt_RunNumber"   , &stuple.evt_RunNumber   , "evt_RunNumber/i");
		tree->Branch("evt_EventNumber" , &stuple.evt_EventNumber , "ect_EventNumber/i");

		tree->Branch("tri_pho25iso" , &stuple.tri_pho25iso , "tri_pho25iso/I");
		tree->Branch("tri_pho50"    , &stuple.tri_pho50    , "tri_pho50/I");
		tree->Branch("tri_pho70"    , &stuple.tri_pho70    , "tri_pho70/I");
	
		tree->Branch("pho_num"     , &stuple.pho_num      , "pho_num/i");
		
		tree->Branch("pho_Ntight"  , &stuple.pho_Ntight   , "pho_Ntight/i");
		tree->Branch("pho_Nloose"  , &stuple.pho_Nloose   , "pho_Nloose/i");
		tree->Branch("pho_Index"   , &stuple.pho_Index    , "pho_Index[pho_num]/I");
		tree->Branch("pho_PhoBlockIndex"   , &stuple.pho_PhoBlockIndex    , "pho_PhoBlockIndex[pho_num]/I");
		tree->Branch("pho_Etc"     , &stuple.pho_Etc      , "pho_Etc[pho_num]/F");
		tree->Branch("pho_E"       , &stuple.pho_E        , "pho_E[pho_num]/F");
		tree->Branch("pho_Px"      , &stuple.pho_Px       , "pho_Px[pho_num]/F");
		tree->Branch("pho_Py"      , &stuple.pho_Py       , "pho_Py[pho_num]/F");
		tree->Branch("pho_Pz"      , &stuple.pho_Pz       , "pho_Pz[pho_num]/F");
		tree->Branch("pho_Detector", &stuple.pho_Detector , "pho_Detector[pho_num]/I");
		tree->Branch("pho_DetEta"  , &stuple.pho_DetEta   , "pho_DetEta[pho_num]/F");
		tree->Branch("pho_DetPhi"  , &stuple.pho_DetPhi   , "pho_DetPhi[pho_num]/F");
		tree->Branch("pho_XCes"    , &stuple.pho_XCes     , "pho_XCes[pho_num]/F");
		tree->Branch("pho_ZCes"    , &stuple.pho_ZCes     , "pho_ZCes[pho_num]/F");
		tree->Branch("pho_HadEm"   , &stuple.pho_HadEm    , "pho_HadEm[pho_num]/F");
		tree->Branch("pho_Chi2Mean", &stuple.pho_Chi2Mean , "pho_Chi2Mean[pho_num]/F");
		tree->Branch("pho_N3d"     , &stuple.pho_N3d      , "pho_N3d[pho_num]/I");
		tree->Branch("pho_Iso4"    , &stuple.pho_Iso4     , "pho_Iso4[pho_num]/F");
		tree->Branch("pho_TrkPt"   , &stuple.pho_TrkPt    , "pho_TrkPt[pho_num]/F");
		tree->Branch("pho_TrkIso"  , &stuple.pho_TrkIso   , "pho_TrkIso[pho_num]/F");
		tree->Branch("pho_CesWireE2"      , &stuple.pho_CesWireE2      , "pho_CesWireE2[pho_num]/F");
		tree->Branch("pho_CesStripE2"     , &stuple.pho_CesStripE2     , "pho_CesStripE2[pho_num]/F");
		tree->Branch("pho_PhiWedge"       , &stuple.pho_PhiWedge       , "pho_PhiWedge[pho_num]/I");
		tree->Branch("pho_NMuonStubs"     , &stuple.pho_NMuonStubs     , "pho_NMuonStubs[pho_num]/I");
		tree->Branch("pho_EmTime"         , &stuple.pho_EmTime         , "pho_EmTime[pho_num]/F");
		tree->Branch("pho_TightId"        , &stuple.pho_TightId        , "pho_TightId[pho_num]/I");
		tree->Branch("pho_LooseId"        , &stuple.pho_LooseId        , "pho_LooseId[pho_num]/I");
		tree->Branch("pho_PhoenixId"      , &stuple.pho_PhoenixId      , "pho_PhoenixId[pho_num]/I");
		tree->Branch("pho_Halo_seedWedge" , &stuple.pho_Halo_seedWedge , "pho_Halo_seedWedge[pho_num]/I");
		tree->Branch("pho_Halo_eastNhad"  , &stuple.pho_Halo_eastNhad  , "pho_Halo_eastNhad[pho_num]/I");
		tree->Branch("pho_Halo_westNhad"  , &stuple.pho_Halo_westNhad  , "pho_Halo_westNhad[pho_num]/I");
		tree->Branch("pho_matchJetIndex"  , &stuple.pho_matchJetIndex  , "pho_matchJetIndex[pho_num]/I");
		

		tree->Branch("ele_num"     , &stuple.ele_num      , "ele_num/i");
		tree->Branch("ele_Ntight"  , &stuple.ele_Ntight   , "ele_Ntight/i");
		tree->Branch("ele_Nloose"  , &stuple.ele_Nloose   , "ele_Nloose/i");
		tree->Branch("ele_Index"   , &stuple.ele_Index    , "ele_Index[ele_num]/I");
		tree->Branch("ele_EleBlockIndex", &stuple.ele_EleBlockIndex, "ele_EleBlockIndex[ele_num]/I");
		tree->Branch("ele_Etc"     , &stuple.ele_Etc      , "ele_Etc[ele_num]/F");
		tree->Branch("ele_E"       , &stuple.ele_E        , "ele_E[ele_num]/F");
		tree->Branch("ele_Px"      , &stuple.ele_Px       , "ele_Px[ele_num]/F");
		tree->Branch("ele_Py"      , &stuple.ele_Py       , "ele_Py[ele_num]/F");
		tree->Branch("ele_Pz"      , &stuple.ele_Pz       , "ele_Pz[ele_num]/F");
		tree->Branch("ele_Detector", &stuple.ele_Detector , "ele_Detector[ele_num]/I");
		tree->Branch("ele_DetEta"  , &stuple.ele_DetEta   , "ele_DetEta[ele_num]/F");
		tree->Branch("ele_DetPhi"  , &stuple.ele_DetPhi   , "ele_DetPhi[ele_num]/F");
		tree->Branch("ele_XCes"    , &stuple.ele_XCes     , "ele_XCes[ele_num]/F");
		tree->Branch("ele_ZCes"    , &stuple.ele_ZCes     , "ele_ZCes[ele_num]/F");
		tree->Branch("ele_HadEm"   , &stuple.ele_HadEm    , "ele_HadEm[ele_num]/F");
		tree->Branch("ele_Chi2Mean", &stuple.ele_Chi2Mean , "ele_Chi2Mean[ele_num]/F");
		tree->Branch("ele_N3d"     , &stuple.ele_N3d      , "ele_N3d[ele_num]/I");
		tree->Branch("ele_Iso4"    , &stuple.ele_Iso4     , "ele_Iso4[ele_num]/F");
		tree->Branch("ele_TrkIso"  , &stuple.ele_TrkIso   , "ele_TrkIso[ele_num]/F");
		tree->Branch("ele_CesWireE2"      , &stuple.ele_CesWireE2      , "ele_CesWireE2[ele_num]/F");
		tree->Branch("ele_CesStripE2"     , &stuple.ele_CesStripE2     , "ele_CesStripE2[ele_num]/F");
		tree->Branch("ele_PhiWedge"       , &stuple.ele_PhiWedge       , "ele_PhiWedge[ele_num]/I");
		tree->Branch("ele_NMuonStubs"     , &stuple.ele_NMuonStubs     , "ele_NMuonStubs[ele_num]/I");
		tree->Branch("ele_EmTime"         , &stuple.ele_EmTime         , "ele_EmTime[ele_num]/F");
		tree->Branch("ele_PhoenixId"      , &stuple.ele_PhoenixId      , "ele_PhoenixId[ele_num]/I");
		tree->Branch("ele_Halo_seedWedge" , &stuple.ele_Halo_seedWedge , "ele_Halo_seedWedge[ele_num]/I");
		tree->Branch("ele_Halo_eastNhad"  , &stuple.ele_Halo_eastNhad  , "ele_Halo_eastNhad[ele_num]/I");
		tree->Branch("ele_Halo_westNhad"  , &stuple.ele_Halo_westNhad  , "ele_Halo_westNhad[ele_num]/I");
		tree->Branch("ele_matchJetIndex"  , &stuple.ele_matchJetIndex  , "ele_matchJetIndex[ele_num]/I");
	
		tree->Branch("ele_Ntracks"      , &stuple.ele_Ntracks      , "ele_Ntracks[ele_num]/I");
		tree->Branch("ele_Emfr"         , &stuple.ele_Emfr         , "ele_Emfr[ele_num]/F");
		tree->Branch("ele_EoverP"       , &stuple.ele_EoverP       , "ele_EoverP[ele_num]/F");
		tree->Branch("ele_TrackPt"      , &stuple.ele_TrackPt      , "ele_TrackPt[ele_num]/F");
		tree->Branch("ele_TrackBcPt"    , &stuple.ele_TrackBcPt    , "ele_TrackBcPt[ele_num]/F");
		tree->Branch("ele_TrackPhi"     , &stuple.ele_TrackPhi     , "ele_TrackPhi[ele_num]/F");
		tree->Branch("ele_Nssl"         , &stuple.ele_Nssl         , "ele_Nssl[ele_num]/I");
		tree->Branch("ele_Nasl"         , &stuple.ele_Nasl         , "ele_Nasl[ele_num]/I");
		tree->Branch("ele_TightId"      , &stuple.ele_TightId      , "ele_TightId[ele_num]/I");
		tree->Branch("ele_LooseId"      , &stuple.ele_LooseId      , "ele_LooseId[ele_num]/I");
		tree->Branch("ele_ConversionId" , &stuple.ele_ConversionId , "ele_ConversionId[ele_num]/I");
		tree->Branch("ele_matchJetIndex", &stuple.ele_matchJetIndex, "ele_matchJetIndex[ele_num]/I");
		
		tree->Branch("jet_num"      , &stuple.jet_num      , "jet_num/i");
		tree->Branch("jet_NJet15"   , &stuple.jet_NJet15   , "jet_Njet15/I");
		tree->Branch("jet_Index"    , &stuple.jet_Index    , "jet_Index[jet_num]/I");
		tree->Branch("jet_Pt"       , &stuple.jet_Pt       , "jet_Pt[jet_num]/F");
		tree->Branch("jet_E"        , &stuple.jet_E        , "jet_E[jet_num]/F");
		tree->Branch("jet_Px"       , &stuple.jet_Px       , "jet_Px[jet_num]/F");
		tree->Branch("jet_Py"       , &stuple.jet_Py       , "jet_Py[jet_num]/F");
		tree->Branch("jet_Pz"       , &stuple.jet_Pz       , "jet_Pz[jet_num]/F");
		tree->Branch("jet_DetEta"   , &stuple.jet_DetEta   , "jet_DetEta[jet_num]/F");
		tree->Branch("jet_DetPhi"   , &stuple.jet_DetPhi   , "jet_DetPhi[jet_num]/F");
		tree->Branch("jet_HadEm"    , &stuple.jet_HadEm    , "jet_HadEm[jet_num]/F");
		tree->Branch("jet_Emfr"     , &stuple.jet_Emfr     , "jet_Emfr[jet_num]/F");
		tree->Branch("jet_Ntowers"  , &stuple.jet_Ntowers  , "jet_Ntowers[jet_num]/I");
		tree->Branch("jet_Ntracks"  , &stuple.jet_Ntracks  , "jet_Ntracks[jet_num]/I");
		tree->Branch("jet_SeedIPhi" , &stuple.jet_SeedIPhi , "jet_SeedIPhi[jet_num]/I");
		tree->Branch("jet_SeedIEta" , &stuple.jet_SeedIEta , "jet_SeedIEta[jet_num]/I");
		
		tree->Branch("vtx_N"        , &stuple.vtx_N        , "vtx_N/I");
		tree->Branch("vtx_NClass12" , &stuple.vtx_NClass12 , "vtx_NClass12/I");
		tree->Branch("vtx_z"        , &stuple.vtx_z        , "vtx_z/F");
		tree->Branch("vtx_Ntracks"  , &stuple.vtx_Ntracks  , "vtx_Ntracks/I");
		tree->Branch("vtx_SumPt"    , &stuple.vtx_SumPt    , "vtx_SumPt/F");

		tree->Branch("met_Met"    , &stuple.met_Met    , "met_Met/F");
		tree->Branch("met_SumEt"  , &stuple.met_SumEt  , "met_SumEt/F");
		tree->Branch("met_Ht"     , &stuple.met_Ht     , "met_Ht/F");
		tree->Branch("met_MetPhi" , &stuple.met_MetPhi , "met_MetPhi/F");

		//tree->Branch("met_Gen_d" , &stuple.met_Gen_d , "met_Gen_d/F");
		//tree->Branch("met_Gen_m" , &stuple.met_Gen_m , "met_Gen_m/F");
		//tree->Branch("met_Gen_p" , &stuple.met_Gen_p , "met_Gen_p/F");
		//tree->Branch("met_Gen_mUn" , &stuple.met_Gen_mUn , "met_Gen_mUn/F");
		//tree->Branch("met_Gen_pUn" , &stuple.met_Gen_pUn , "met_Gen_pUn/F");

		

  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}

  	jetMod = (JetFilterModuleV2*) ((TStnAna*) GetAna()->GetModule("JetFilterV2"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		bRunPermit = false;
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
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
		
	BookHistograms();

	counter.evtsProcessed		= 0;
	counter.evtsPassModule 		= 0;
	
	return 0;
}

//_____________________________________________________________________________
int FlatStupleMaker_MC::BeginRun()
{
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int FlatStupleMaker_MC::Event(int ientry)
{
	SetPassed(0);
	counter.evtsProcessed++;
	if (! bRunPermit) exit (1);

	stuple.Init();
	//std::cout << "=====================";
	//GetHeaderBlock()->Print();
	
	std::vector<int> vPho,vEle ;
	for (int i=0; i < initSpMod->GetSuperPhoSize(); i++) {
		
		if (initSpMod->GetSuperPhoton(i)->IsTightPhoton() || initSpMod->GetSuperPhoton(i)->IsLoosePhoton()) {
		//std::cout << "\t i, et, Eta=" << i << ", " << initSpMod->GetSuperPhoton(i)->GetEtCorr() << ", " <<initSpMod->GetSuperPhoton(i)->GetDetEta()<< std::endl;
			stuple.pho_num++;
			if (initSpMod->GetSuperPhoton(i)->IsTightPhoton())	stuple.pho_Ntight++;
			if (initSpMod->GetSuperPhoton(i)->IsLoosePhoton()) stuple.pho_Nloose++;
			vPho.push_back(i);
		}
		if (initSpMod->GetSuperPhoton(i)->IsTightElectron() || initSpMod->GetSuperPhoton(i)->IsLooseElectron()) {
			stuple.ele_num++;
			if (initSpMod->GetSuperPhoton(i)->IsTightElectron()) stuple.ele_Ntight++;
			if (initSpMod->GetSuperPhoton(i)->IsLooseElectron()) stuple.ele_Nloose++;
			vEle.push_back(i);
		}
	}


	assert (stuple.pho_num <= Stuple::Npho);
	assert (stuple.ele_num <= Stuple::Nele);
	if (stuple.pho_num + stuple.ele_num > 0) {
		//std::cout << "tights=" << stuple.pho_Ntight << std::endl; 
		//std::cout << "loose=" << stuple.pho_Nloose << std::endl; 
		//std::cout << "phosize, elesize=" << vPho.size() << ", " << vEle.size() << std::endl; 
		FillDataBlocks(ientry);
		// now I need to match the tight photon to gen level to make sure this is the photon. is so save it. else reject!
		if (MatchToGen(vPho)) {
			SkimEvent(vPho, vEle);
			counter.evtsWritten++;
		} else SetPassed(0);

		//this is to debug 2 things. ele removal problem and jet det eta >3
		for (unsigned int i=0; i < stuple.jet_num; ++i) {
			if (stuple.jet_DetEta[i] > 2.8) SetPassed(1);
			break;
		}
		for (unsigned int i=0; i < stuple.ele_num; ++i) {
			if (stuple.ele_matchJetIndex[i] <0) SetPassed(1);
			std::cout << __FILE__ << "::"<< __LINE__ <<":: No jet match event is " << std::endl;
			GetHeaderBlock()->Print();
			break;
		}
		
		
	} else SetPassed(0);



	if (GetPassed()) counter.evtsPassModule++;
	return 0;

} // Event



/*-------------------------------------------------------------------*/
bool FlatStupleMaker_MC::MatchToGen(std::vector<int> vPho)
{

	for (unsigned int i =0 ;i < vPho.size(); i++) {
		if (! initSpMod->GetSuperPhoton(vPho[i])->IsTightPhoton()) continue;
		TLorentzVector detvec(0,0,0,0), genvec(0,0,0,0);
		detvec = initSpMod->GetSuperPhoton(vPho[i])->GetCorVec();
		
		
		//int mom_id = 0;
		int Nparticles = fGenpBlock->NParticles();

		TGenParticle *par, *mom;
	
		for (int i = 0 ; i < Nparticles ; i++) {	//_____ not much diff from running over all. my match is mostly within the fist 10 objects
			par = fGenpBlock->Particle(i);
			int im = par->GetFirstMother();

			if (im >= 0) {	//_________________________________im=-1 means no mother for imcoming particles
				mom = fGenpBlock->Particle(im);

				if (mom != 0) {
					//int par_id   = par->GetPdgCode();
					int par_stat = par->GetStatusCode();
				
					//if (abs(par_id) == 22 && (par_stat == 2)) { 
					if (par->IsPhoton() && par_stat == 1) { 
						par->Momentum(genvec);
						//std::cout << "i=" << i << "\tid=" << par_id << "\tstat=" << par_stat  << "\teta=" << genvec.Eta() 
						//	<< "\tEt=" << genvec.E() * TMath::Sin(genvec.Theta()) << std::endl;
						//std::cout << "\tDELR=" << genvec.DeltaR(detvec);						
						hDelR->Fill(genvec.DeltaR(detvec));	
						if (genvec.DeltaR(detvec) < 0.1) return true;
					}
				} 
			} 
		} // for 
	} // for
	
	
	//std::cout << "ERROR::NO MATCHING PHOTON FOUND TO THE DET PHO." << std::endl;
	GetHeaderBlock()->Print();
	return false;

} // MatchToGen


/*-------------------------------------------------------------------*/
void FlatStupleMaker_MC::SkimEvent(std::vector<int> vPho, std::vector<int> vEle)
{

	stuple.evt_McFlag = GetHeaderBlock()->McFlag();
	stuple.evt_RunNumber = GetHeaderBlock()->RunNumber();
	stuple.evt_EventNumber = GetHeaderBlock()->EventNumber();
	
	for (unsigned int i=0; i< vPho.size(); i++) {
		GetPhotonInfo(i, *(initSpMod->GetSuperPhoton(vPho[i])));
		TagBeamHalo::HaloStuff aHaloPara = tagHaloMod->GetHaloStuff(vPho[i]);
		stuple.pho_Halo_seedWedge[i] = aHaloPara.seedwedge;
		stuple.pho_Halo_eastNhad[i] = aHaloPara.eastNhad;
		stuple.pho_Halo_westNhad[i] = aHaloPara.westNhad;
	}
	for (unsigned int i=0; i< vEle.size(); i++) {
		GetElectronInfo(i, *(initSpMod->GetSuperPhoton(vEle[i])));
		TagBeamHalo::HaloStuff aHaloPara = tagHaloMod->GetHaloStuff(vEle[i]);
		stuple.ele_Halo_seedWedge[i] = aHaloPara.seedwedge;
		stuple.ele_Halo_eastNhad[i] = aHaloPara.eastNhad;
		stuple.ele_Halo_westNhad[i] = aHaloPara.westNhad;
	
	}

	//GetTriggerInfo();  // Only L1/L2 triggers are simulated(if any). not L3. so drop it -- 03-27-08
	GetVertexInfo();
	GetJetInfo();
	GetMetInfo();
	tree->Fill();

}


/*-------------------------------------------------------------------*/
void FlatStupleMaker_MC::GetTriggerInfo()
{
	stuple.tri_pho25iso 	= trigMod->GetTrig25IsoBit();
	stuple.tri_pho50 		= trigMod->GetTrig50Bit();
	stuple.tri_pho70 		= trigMod->GetTrig70Bit();
}
	
/*-------------------------------------------------------------------*/
void FlatStupleMaker_MC::GetVertexInfo()
{
	stuple.vtx_N = trigMod->GetNvertex();
	stuple.vtx_NClass12 = trigMod->GetN12vertex();
	stuple.vtx_z = trigMod->GetBestVertexZ();	
	stuple.vtx_Ntracks = trigMod->GetBestVertexNTracks();
	stuple.vtx_SumPt = trigMod->GetBestVertexSumPt();
}

/*-------------------------------------------------------------------*/
void FlatStupleMaker_MC::GetMetInfo()
{
/*
	stuple.met_Met 	= jetMod->GetMyMetCorr(0,3);
	stuple.met_SumEt 	= jetMod->GetMySumEtCorr(0,3);
	stuple.met_Ht 		= jetMod->GetMyHtCorr(0,3);
	stuple.met_MetPhi = jetMod->GetMyMetPhiCorr(0,3);

	TVector2 GenMet_d(jetMod->GetGenMet_Def());
	TVector2 GenMet_m(jetMod->GetGenMet_devm());
	TVector2 GenMet_p(jetMod->GetGenMet_devp());
	TVector2 GenMet_mUn(jetMod->GetGenMet_devmUn());
	TVector2 GenMet_pUn(jetMod->GetGenMet_devpUn());
	
	stuple.met_Gen_d = GenMet_d.Mod();
	stuple.met_Gen_m = GenMet_m.Mod();
	stuple.met_Gen_p = GenMet_p.Mod();
	stuple.met_Gen_mUn = GenMet_mUn.Mod();
	stuple.met_Gen_pUn = GenMet_pUn.Mod();
*/
}
	

/*-------------------------------------------------------------------*/
void FlatStupleMaker_MC::GetJetInfo()
{
/*	stuple.jet_num = jetMod->GetMyNjet(0,3);
	stuple.jet_NJet15 = jetMod->GetMyNjet(0,3);


	for (int i=0; i < stuple.jet_NJet15; i++) {
		
		stuple.jet_Index[i] = jetMod->GetMyJetBlockInd(0,i);
		stuple.jet_Pt[i] = jetMod->GetMyJetNoEMobj_lev6(0,i)->Pt();	
		stuple.jet_E[i] = jetMod->GetMyJetNoEMobj_lev6(0,i)->E();
		stuple.jet_Px[i] = jetMod->GetMyJetNoEMobj_lev6(0,i)->Px();
		stuple.jet_Py[i] = jetMod->GetMyJetNoEMobj_lev6(0,i)->Py();
		stuple.jet_Pz[i] = jetMod->GetMyJetNoEMobj_lev6(0,i)->Pz();
		stuple.jet_DetEta[i] = jetMod->GetMyJetEtaDetCorr(0,i);
		//stuple.jet_DetPhi[i] = jetMod->GetMy(0,i);
		//stuple.jet_HadEm[i] = jetMod->GetMy(0,i);
		stuple.jet_Emfr[i] = jetMod->GetMyJetEmFrCorr(0,i);
		stuple.jet_Ntowers[i] = jetMod->GetMyJetNtwr(0,i);
		stuple.jet_Ntracks[i] = jetMod->GetMyJetNtrk(0,i);

		TStnJet* jet = fJetBlock->Jet(stuple.jet_Index[i]);
		stuple.jet_SeedIPhi[i] = jet->SeedIEta();
		stuple.jet_SeedIEta[i] = jet->SeedIPhi();
	}
*/
}

/*-------------------------------------------------------------------*/
void FlatStupleMaker_MC::GetPhotonInfo(int iArrInd, SuperPhoton sp)
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
	stuple.pho_TrkIso[iArrInd] 	= sp.GetTrkIso();
	stuple.pho_CesWireE2[iArrInd] 	= sp.GetCesWireE2();
	stuple.pho_CesStripE2[iArrInd] 	= sp.GetCesStripE2();
	stuple.pho_PhiWedge[iArrInd] 		= sp.GetPhoton()->PhiSeedIndex();
	stuple.pho_NMuonStubs[iArrInd] 	= sp.GetNMuonStubs();
	stuple.pho_EmTime[iArrInd] 		= sp.GetEmTime();
	stuple.pho_TightId[iArrInd] 		= sp.GetTightPhotonId();
	stuple.pho_LooseId[iArrInd] 		= sp.GetLoosePhotonId();
	stuple.pho_PhoenixId[iArrInd] 	= sp.GetPhoenixId();
	stuple.pho_matchJetIndex[iArrInd] = jetMod->GetMatch_JetIndex(0,sp.GetPhotonBlockIndex());//0=ObjType(photon) 
	if (stuple.pho_matchJetIndex[iArrInd] < 0) {
		StdOut(__FILE__,__LINE__,3,"pho did not match to a jet.");
		std::cout << "ERROR_PHO::"<< GetHeaderBlock()->RunNumber() << "," << GetHeaderBlock()->RunNumber() << std::endl;
		//exit (1);
	} 
	

}
/*-------------------------------------------------------------------*/
void FlatStupleMaker_MC::GetElectronInfo(int iArrInd, SuperPhoton sp)
{
	stuple.ele_EleBlockIndex[iArrInd] = sp.GetElectronBlockIndex();
	stuple.ele_Index[iArrInd] 		= sp.GetPhotonBlockIndex();
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

	stuple.ele_TightId[iArrInd] = sp.GetTightElectronId();
	stuple.ele_LooseId[iArrInd] = sp.GetLooseElectronId();
	stuple.ele_ConversionId[iArrInd] = sp.GetConversionId();

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
		}
	}
	stuple.ele_matchJetIndex[iArrInd] = jetMod->GetMatch_JetIndex(1,sp.GetElectronBlockIndex());//1=ObjType(electron) 
	if (stuple.ele_matchJetIndex[iArrInd] < 0) {
		StdOut(__FILE__,__LINE__,3," ele did not match to a jet.");
		std::cout << "ERROR_ELE::"<< GetHeaderBlock()->RunNumber() << "," << GetHeaderBlock()->RunNumber() << std::endl;
		//exit (1);
	} 
															 
}

/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void FlatStupleMaker_MC::Cleanup()
{
} //Cleanup

/*-------------------------------------------------------------------*/
void FlatStupleMaker_MC::FillDataBlocks(int ientry)
{
	fElectronBlock->GetEntry(ientry);
	fJetBlock->GetEntry(ientry);		// need to get Jet tower IEta/IPhi
	fGenpBlock->GetEntry(ientry);
}
//____________________________________________________________________
//  END JOB SUMMARY
//____________________________________________________________________
int FlatStupleMaker_MC::EndJob() {

	if (rootFile->IsOpen()) {
		rootFile->Write();
		rootFile->Close();
	}
	printf("[FSM:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[FSM:01:] Events Processed ----------- = " << counter.evtsProcessed << std::endl;
	std::cout << "[FSM:02:] Events Pass this module ---- = " << counter.evtsPassModule << std::endl;
	std::cout << "[FSM:03:] Events Written to Stuple --- = " << counter.evtsWritten << std::endl;
	std::cout << "[FSM:04:] MC flag  ------------------- = " << GetHeaderBlock()->McFlag() << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}