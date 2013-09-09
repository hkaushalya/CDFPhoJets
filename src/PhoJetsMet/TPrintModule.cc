/* HOW TO USE:
1. init with header Block
2. add objects (optional)
3. call Print with options
*/
#include "samantha/PhoJetsMet/TPrintModule.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/ana/TGammaJetsInit.hh"
#include <string>

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
TPrintModule::TPrintModule(const char* name, const char* title,TStnHeaderBlock* fHeaderBlock)
{
	Init(fHeaderBlock);
}

//------------------------------------------------------------------------------
// Overloaded Constructor - 1
//------------------------------------------------------------------------------
TPrintModule::TPrintModule(TStnHeaderBlock* fHeaderBlock)
{
	Init(fHeaderBlock);
}

//------------------------------------------------------------------------------
// Overloaded Constructor - 2  (for quick dump of one photon)
//------------------------------------------------------------------------------
TPrintModule::TPrintModule(TStnHeaderBlock* fHeaderBlock, TStnPhoton* pho)
{
	Init(fHeaderBlock);
	if (pho) photon_list.push_back(pho);
	else cout << "TPrintModule(TStnHeaderBlock*,TStnPhoton*): NULL Photon object.!\n"; 
}


//------------------------------------------------------------------------------
// Overloaded Constructor - 3 (for quick dump of one jet)
//------------------------------------------------------------------------------
TPrintModule::TPrintModule(TStnHeaderBlock* fHeaderBlock, TStnJet* jet)
{
	Init(fHeaderBlock);
	if (jet) jet_list.push_back(jet);
	else cout << "TPrintModule(TStnHeaderBlock*,TStnJet*): NULL Jet object!\n"; 
}

//------------------------------------------------------------------------------
// Overloaded Constructor - 4 (for quick dump of one hepg particle)
//------------------------------------------------------------------------------
TPrintModule::TPrintModule(TStnHeaderBlock* fHeaderBlock, TGenParticle* par)
{
	Init(fHeaderBlock);
	if (par) hepgpar_list.push_back(par);
	else cout << "TPrintModule(TStnHeaderBlock*,TGenParticle*): NULL HEPG particle!\n"; 
}

//------------------------------------------------------------------------------
// Overloaded Constructor - 5  (for quick dump of one electronon)
//------------------------------------------------------------------------------
TPrintModule::TPrintModule(TStnHeaderBlock* fHeaderBlock, TStnElectron* ele)
{
	Init(fHeaderBlock);
	if (ele) electron_list.push_back(ele);
	else cout << "TPrintModule(TStnHeaderBlock*,TStnElectron*): NULL Electron object.!\n"; 
}

//------------------------------------------------------------------------------
// Overloaded Constructor - 6
//------------------------------------------------------------------------------
TPrintModule::TPrintModule(TStnHeaderBlock* fHeaderBlock, TStnPhoton* pho, TStnJet* jet)
{
	Init(fHeaderBlock);
	if (pho)	photon_list.push_back(pho);
	else cout << "TPrintModule(TStnPhoton*,TStnJet*): No Photon to add!\n"; 

	if (jet) jet_list.push_back(jet);
	else cout << "TPrintModule(TStnPhoton*,TStnJet*): No Jet to add!\n"; 
	
}



//------------------------------------------------------------------------------
// Default destructor
//------------------------------------------------------------------------------
TPrintModule::~TPrintModule()
{
}


//------------------------------------------------------------------------------
// Initialize Class Variables
//------------------------------------------------------------------------------
void TPrintModule::Init(TStnHeaderBlock* fHeaderBlock)
{
	
	photonheader 	= true;
	jetheader 		= true;
	electronheader = true;
	hepgheader		= true;
	photoncounter 	= 0;
	jetcounter 		= 0;
	electroncounter= 0;
	hepgcounter    = 0;
	if (fHeaderBlock)	{
		eventNumber = fHeaderBlock->EventNumber();
		runNumber   = fHeaderBlock->RunNumber();
	} else {
		eventNumber = -0;
		runNumber = -0;
	}
}


//------------------------------------------------------------------------------
// Reset All Class Variables
// this will reset all to initial values
// we'll keep the fHeaderBlock.
//------------------------------------------------------------------------------
void TPrintModule::Reset()
{
	
	photonheader 	= true;
	jetheader 		= true;
	electronheader = true;
	hepgheader		= true;
	photoncounter 	= 0;
	jetcounter 		= 0;
	electroncounter= 0;
	hepgcounter		= 0;
	eventNumber = -99999999;
	runNumber = -99999999;
	jet_list.clear();
	photon_list.clear();
	hepgpar_list.clear();
}


//------------------------------------------------------------------------------
// Add a Photon to the dump 
//------------------------------------------------------------------------------
void TPrintModule::AddThisPhoton(TStnPhoton *pho)
{
	if (pho)	photon_list.push_back(pho);
	else cout << "TPrintModule::AddThisPhoton(TStnPhoton*): No Photon to add!\n"; 
}
	

//------------------------------------------------------------------------------
// Add a Jet to the dump 
//------------------------------------------------------------------------------
void TPrintModule::AddThisJet(TStnJet *jet)
{
	if (jet) jet_list.push_back(jet);
	else cout << "TPrintModule::AddThisJet(TStnJet*): No Jet to add!\n" ;
}


//------------------------------------------------------------------------------
// Add an Electron to the dump 
//------------------------------------------------------------------------------
void TPrintModule::AddThisElectron(TStnElectron *ele)
{
	if (ele) electron_list.push_back(ele);
	else cout << "TPrintModule::AddThisElectron(TStnElectron*): No Electron to add!\n" ;
}


//------------------------------------------------------------------------------
// Add a HEPG particle to the dump 
// warning: make sure that u pass in a valid stable particle.
// see the method used in HEPG Block dump method.
//------------------------------------------------------------------------------
void TPrintModule::AddThisHEPGparticle(TGenParticle *par)
{
	if (par) hepgpar_list.push_back(par);
	else cout << "TPrintModule::AddThisHEPGparticle(TGenParticle*): No HEPG particle to add!\n" ;
}

//------------------------------------------------------------------------------
//  print method: spit out things in the object vectors 
//	 Options : a=all  p=photons j=jets
//------------------------------------------------------------------------------
void TPrintModule::Print(Option_t *opt) 
{

	
	printf("\n<<<<<<<<<<<<<<<<<<<<< Event Summary for Run,Event");
	printf(": %8i,%8i >>>>>>>>>>>>>>>>>>>>>>>>>>>\n",runNumber,eventNumber);
	if ( ( strchr(opt,'a') || strchr(opt,'p') ) && (photon_list.size()>0) ) {
		for (int i=0; i < photon_list.size() ; i++ ) PhotonInfo(photon_list[i]);
	} 
	if ( ( strchr(opt,'a') || strchr(opt,'e') ) && (electron_list.size()>0) ) {
		for (int i=0; i < electron_list.size() ; i++ ) ElectronInfo(electron_list[i]);
	} 
	if ( ( strchr(opt,'a') || strchr(opt,'j') ) && (jet_list.size()>0) ) {
		for (int i=0; i < jet_list.size() ; i++ ) JetInfo(jet_list[i]);
	}
	if ( ( strchr(opt,'a') || strchr(opt,'h') ) && (hepgpar_list.size()>0) ) {
		for (int i=0; i < hepgpar_list.size() ; i++ ) HepgInfo(hepgpar_list[i]);
	}
}


//------------------------------------------------------------------------------
// 				DUMP WHOLE PHOTON BLOCK
//------------------------------------------------------------------------------
void TPrintModule::Print(TStnPhotonBlock* fPhotonBlock) 
{
	photon_list.clear();
	for (int i=0; i< fPhotonBlock->NPhotons() ;i++) {
		TStnPhoton* pho = fPhotonBlock->Photon(i);
		AddThisPhoton(pho);
	}
	Print("p");
}


//------------------------------------------------------------------------------
// 				DUMP WHOLE ELECTRON BLOCK
//------------------------------------------------------------------------------
void TPrintModule::Print(TStnElectronBlock* fElectronBlock) 
{
	electron_list.clear();
	for (int i=0; i< fElectronBlock->NElectrons() ;i++) {
		TStnElectron* ele = fElectronBlock->Electron(i);
		AddThisElectron(ele);
	}
	Print("e");
}


//------------------------------------------------------------------------------
// 						DUMP WHOLE JET BLOCK
//------------------------------------------------------------------------------
void TPrintModule::Print(TStnJetBlock* fJetBlock) 
{
	jet_list.clear();
	for (int i=0; i< fJetBlock->NJets() ;i++) {
		TStnJet* jet = fJetBlock->Jet(i);
		AddThisJet(jet);
	}
	Print("j");
}

//------------------------------------------------------------------------------
// 						DUMP WHOLE GENP BLOCK
//------------------------------------------------------------------------------
void TPrintModule::Print(TGenpBlock* fGenpBlock) 
{
	hepgpar_list.clear();
	for (int i=0; i< fGenpBlock->NParticles() ;i++) {
		TGenParticle* par = fGenpBlock->Particle(i);
		int im = par->GetFirstMother();
		if (im>=0) {
			TGenParticle* mom = fGenpBlock->Particle(im);
			if (mom != 0) {
				AddThisHEPGparticle(par);
			}
		}
	}
	Print("h");
}

//------------------------------------------------------------------------------
//             HEPG PARTICLE DUMP
//------------------------------------------------------------------------------
void TPrintModule::HepgInfo(TGenParticle* par)
{

	int pdgcode = par->GetPdgCode();
	double Et	= par->Energy() * TMath::Sin(par->Theta());
	int    status = par->GetStatusCode(); 
	double eta	= 0;
	double phi	= 0;
	TLorentzVector vec;
	par->Momentum(vec);
	//double charge = par->Charge();
	double charge = -9.9;
	int im = par->GetFirstMother();
/*	int mom_id = -99999;
	if (im >=0) {
		if (mom!=0) {
			mom_id = mom->GetPdgCode();
		}
	}
*/	

if (hepgheader) {
	printf("\n\t------------------ HEPG Block Summary ---------------------\n");
	
	printf("      PDGcode   Et    Status   Eta    Phi");
	printf("    E        Px       Py       Pz\n");
	hepgheader = false;
}
	if (vec.Pt() !=0) {
		eta = vec.Eta();	
		phi = vec.Phi();
	}
	printf(" %2i# %6i   %6.2f   %2i   %6.2f  %6.2f ",hepgcounter, pdgcode,Et, status, eta, phi);
	printf(" %6.2f   %6.2f  %6.2f   %6.2f\n", vec.Energy(), vec.Px(), vec.Py(), vec.Pz());
//	printf(" %5i\n", mom_id);

	hepgcounter++;
}
//------------------------------------------------------------------------------
//             PHOTON DUMP
//------------------------------------------------------------------------------
void TPrintModule::PhotonInfo(TStnPhoton *pho)
{
	int   Detector 	= -99999999;
	float VertexZ 	= -99999999;

	float Et 			= -99999999.99;
	float Etc 		= -99999999.99;
	float ECorr		= -99999999.99;
	float EtCorr		= -99999999.99;
	float DetEta 	= -99999999.99;
	float Eta 		= -99999999.99;
	float Phi 		= -99999999.99;
	float SinTheta 	= -99999999.99;

	float HadEm 		= -99999999.99;
	float HadEmT 	= -99999999.99;
	float LShr 		= -99999999.99;
	float Time 		= -99999999.99;

	int   PhiSeedIndex = -99999999;
	int   EtaSeedIndex = -99999999;

	float Iso 		= -99999999.99;
	float SumEt4		= -99999999.99;
	float IsoEtCorr	= -99999999.99;
											// corr for phi leakage and multiple vertices
	float IsoCorr 	= -99999999.99;
	float Tiso 		= -99999999.99;
	float SumPt4		= -99999999.99;
	float Pt 			= -99999999.99;
	float Pt2			= -99999999.99;
	int   N3d			= -99999999;
	float Iso7		= -99999999.99;
	float SumEt7		= -99999999.99;
	int   NIsoTrk 	= -99999999;
	float MaxIsoTrk = -99999999.99;

	float XCes 		= -99999999.99;
	float ZCes 		= -99999999.99;
	float Chi2		= -99999999.99;
	float Chi2Mean	= -99999999.99;
	float Chi2Strip	= -99999999.99;
	float Chi2Wire	= -99999999.99;
	float CesEnergy = -99999999.99;
	float CesWireE 	= -99999999.99;
	float CesStripE = -99999999.99;
	float ZCes2		= -99999999.99;
	float CesStripE2= -99999999.99;
	float XCes2		= -99999999.99;
	float CesWireE2	= -99999999.99;
	float CesWht		= -99999999.99;
	float CesSlide	= -99999999.99;

	Detector = pho->Detector();
	VertexZ  = pho->VertexZ();

	Et  		= pho->Et();
	Etc 		= pho->Etc();
	ECorr  	= pho->ECorr();
	EtCorr  	= pho->ECorr() * pho->SinTheta();
	DetEta	= pho->DetEta();
	Eta		= pho->Eta();
	Phi		= pho->Phi();
	SinTheta	= pho->SinTheta();

	HadEm		= pho->HadEm();
	HadEmT	= pho->HadEmT();
	LShr		= pho->LShr();
	Time		= pho->Time();

	PhiSeedIndex	= pho->PhiSeedIndex();
	EtaSeedIndex	= pho->EtaSeedIndex();

	Iso		= pho->Iso();
	SumEt4	= pho->SumEt4();
	IsoEtCorr = pho->EIso4(2);
											// corr for phi leakage and multiple vertices
	IsoCorr	= pho->IsoCorr();
	Tiso		= pho->Tiso();
	SumPt4	= pho->SumPt4();
	Pt			= pho->Pt();
	Pt2		= pho->Pt2();
	N3d		= pho->N3d();
	Iso7		= pho->Iso7();
	SumEt7	= pho->SumEt7();
	NIsoTrk	= pho->NIsoTrk();
	MaxIsoTrk= pho->MaxIsoTrk();

	XCes		= pho->XCes();
	ZCes		= pho->ZCes();
	Chi2		= pho->Chi2();
	Chi2Mean	= pho->Chi2Mean();
	Chi2Strip= pho->Chi2Strip();
	Chi2Wire	= pho->Chi2Wire();
	CesEnergy= pho->CesEnergy();
	CesWireE	= pho->CesWireE();
	CesStripE= pho->CesStripE();
	ZCes2		= pho->ZCes2();
	CesStripE2	= pho->CesStripE2();
	XCes2		= pho->XCes2();
	CesWireE2= pho->CesWireE2();
	CesWht	= pho->CesWht();
	CesSlide	= pho->CesSlide();

if (photonheader) {
	printf("\n\t----------------- Photon Block Summary --------------------\n");
	
	printf("      Det    Et   EtCorr     Zv   ");
	printf(" DetEta  Phi    EvtEta  HadEm   ");
	printf(" Pt   IsoEtCorr   LShr\n");
	photonheader = false;
}
	printf(" %2i#  %2i  %6.2f  %6.2f  %6.2f  ",photoncounter, Detector,Et,EtCorr,VertexZ);
	printf(" %5.2f   %4.2f   %5.2f   %4.3f ",DetEta, Phi,Eta,HadEm);
	printf(" %6.2f  %6.2f    %6.2f\n",Pt, IsoEtCorr,LShr);

	photoncounter++;
}




//------------------------------------------------------------------------------
// 						JET DUMP
//------------------------------------------------------------------------------
void TPrintModule::JetInfo(TStnJet *jet) 
{
	
	Double_t M 		= -99999999.99;
	Double_t Et 		= -99999999.99;
	float  EtCorr = -99999999.99;
	float  Vz 		= -99999999.99;
	float  Dz 		= -99999999.99;
	Double_t Eta 	= -99999999.99;
	Double_t Phi 	= -99999999.99;
	float  DetEta = -99999999.99;
	float  Emfr 	= -99999999.99;
	float  EOverP = -99999999.99;
	float  SeedEt = -99999999.99;
	float  TrackMass 	= -99999999.99;
	float  MaxTrackPt 	= -99999999.99;

	float Corfd 	= -99999999.99;
	float Corfm 	= -99999999.99;
	float Corfa 	= -99999999.99;
	float Tag 		= -99999999.99;
	float Scbpb 	= -99999999.99;


	M		= jet->M();
	Et		= jet->Et();
	EtCorr = jet->EtCorr();
	Vz     = jet->Vz();
	Dz     = jet->Dz();
	Eta		= jet->Eta();
	Phi		= jet->Phi();
	DetEta	= jet->DetEta();
	Emfr	= jet->Emfr();
	EOverP	= jet->EOverP();
	SeedEt	= jet->SeedEt();
	TrackMass	= jet->TrackMass();
	MaxTrackPt	= jet->MaxTrackPt();

	Corfd   = jet->Corfd();
	Corfm   = jet->Corfm();
	Corfa   = jet->Corfa();
	Tag     = jet->Tag();
	Scbpb   = jet->Scbpb();


	if (jetheader) {
		printf("\n\t------------------- Jet Block Summary ---------------------\n");
		
		printf("\t Et   EtCorr    Vz   ");
		printf(" DetEta   Phi   EvtEta   Emfr   P/E    ");
		printf(" MaxTrkPt   TrkMass\n");
		jetheader = false;
	}

	printf(" %2i#  %6.2f  %6.2f  %6.2f ",jetcounter,Et, EtCorr,Vz);
	printf(" %5.2f   %4.3f   %5.2f   %4.3f   %4.3f ",DetEta, Phi,Eta,Emfr,EOverP);
	printf(" %8.2f  %8.2f\n",MaxTrackPt,TrackMass);

	jetcounter++;

}

//------------------------------------------------------------------------------
// 						ELECTRON DUMP
//------------------------------------------------------------------------------
void TPrintModule::ElectronInfo(TStnElectron *ele) 
{
	int    DetectorCode = 0;
	int    TrackNumber  = 0;
	int    NTracks      = 0;
	float  DetEta      = .0;
	float  EvtEta      = .0;
	float  Phi       	= .0;
	float  Et          = .0;
	float  HadEm        = .0;
	float  Emfr         = .0;
	float  EOverP       = .0;
	float  Iso_4        = .0;
	float  IsoCorr      = .0;
	float  TrackPt      = .0;
	float  Z0           = .0;


	DetectorCode = ele->DetectorCode();
	TrackNumber  = ele->TrackNumber ();
	NTracks      = ele->NTracks();
	DetEta      = ele->DetEta();
	EvtEta      = ele->EmClusEvEta();
	Phi       	= ele->EmClusPhi();
	Et          = ele->Et();
	HadEm        = ele->HadEm();
	Emfr         = ele->Emfr();
	EOverP       = ele->EOverP();
	Iso_4        = ele->Iso();
	IsoCorr      = ele->IsoCorr();
	TrackPt      = ele->TrackPt();
	Z0           = ele->Z0();
/*
	int    DetectorCode = ele->DetectorCode();
	int    TrackNumber  = ele->TrackNumber ();
	int    NTracks      = ele->NTracks();
	float  DetEta      = ele->DetEta();
	float  EvtEta      = ele->EmClusEvEta();
	float  Phi       	= ele->EmClusPhi();
	float  Et          = ele->Et();
	float  HadEm        = ele->HadEm();
	float  Emfr         = ele->Emfr();
	float  EOverP       = ele->EOverP();
	float  Iso_4        = ele->Iso();
	float  IsoCorr      = ele->IsoCorr();
	float  TrackPt      = ele->TrackPt();
	float  Z0           = ele->Z0();
*/
	if (electronheader) {
		printf("\n\t---------------- Electron Block Summary -------------------\n");
		
		printf("      Det    Et   ");
		printf(" DetEta   Phi   EvtEta    Trk#   TrkPt ");
		printf(" HadEm  Iso(0.4)\n");
		electronheader = false;
	}

	printf(" %2i#   %1i   %6.2f ",electroncounter,DetectorCode, Et);
	printf(" %5.2f   %4.3f   %5.2f   %5.2f   %5.2f  ",DetEta, Phi,EvtEta, TrackNumber,TrackPt);
	printf(" %4.2f  %4.2f\n",HadEm,Iso_4);

	electroncounter++;

}
