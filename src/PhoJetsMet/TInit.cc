/*

READ THIS EVERYDAY!
keep it simple and readable. forget about performance!
*/
/*
this module will setup photons to be used later. all it will do is sort photons, if any found,
set up 'super photons'  and tag them for later use.
NO good run/trigger/veretx or any other selection should be done here. do it your module!

*/

#include <iostream>
#include <fstream>
#include <TMath.h>
#include <numeric>
#include <vector>
#include <algorithm>
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnTriggerTable.hh"
#include "samantha/PhoJetsMet/TInit.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/ana/TPrintModule.hh"
#include "Stntuple/ana/TMathModule.hh"
#include "Stntuple/ana/TTightPhotonID.hh"

double const kPI = TMath::Pi();
double const kTWOPI = 2.*kPI;
const static unsigned int kPHO_ID_BITS = 19;      //total number photon ID bits (loose + tight)
const static unsigned int kPHO_LOOSE_ID_BITS = 8; //number of Loose ID variables
const static double       kPHO_ET_CUT = 30; 	  //my leading photon. so i make the cut on 30GeV!


const static double 	kETCORR20= 20.0;	// Et SPLIT for ISO CUT
const static double 	kETCORR18= 18.0;	// Et SPLIT for 2nd CES CLUSTER CUT

//loose photon id cut values
const static int 		kL_CENTRAL 	= 0; 		// CENTRAL
const static double	kL_ETCORR	= 30.0;
const static double	kL_XCES		= 21.0;
const static double	kL_ZCES_MIN	= 9.0;
const static double	kL_ZCES_MAX	= 230.0;
const static double	kL_HADEM		= 0.125;
const static double	kL_ISOETCORR_MF_LTETC20 = 0.15;	//ISOETCORR Multiplicative factor when Etc<20
const static double	kL_ISOETCORR_MF_GTETC20 = 3.0;	//when Etc>20
const static double	kL_TRKPT_MF	= 0.25;		// TRKPT Multiplicative factor (*EtCorr)
const static double	kL_TRKISO	= 5.0;

//tight photon id cut values

const static int 		kT_CENTRAL 	= 0; 		// CENTRAL
const static double	kT_ETCORR	= 30.0;
const static double	kT_XCES		= 21.0;
const static double	kT_ZCES_MIN	= 9.0;
const static double	kT_ZCES_MAX	= 230.0;
const static double	kT_HADEM_CONST1	= 0.125;
const static double	kT_HADEM_CONST2	= 0.055;
const static double	kT_HADEM_MF	= 0.00045;
const static double	kT_ISOETCORR_MF_LTETC20 = 0.1;		//ISOETCORR Multiplicative factor when Etc<20
const static double	kT_ISOETCORR_CONST1_GTETC20 = 2.0;	//ADD when Etc>20
const static double	kT_ISOETCORR_CONST2_GTETC20 = 20.0;	//SUBSTRACT when Etc>20
const static double	kT_ISOETCORR_MF_GTETC20 = 0.02;		//mf when Etc>20
const static double	kT_CHI2MEAN	= 20.0;
const static int 		kT_N3DTRACKS	= 1;
const static double	kT_TRKPT_CONST	= 1.0;		// ADD
const static double	kT_TRKPT_MF	= 0.005;	//MF TRKPT (*EtCorr)
const static double	kT_TRKISO_CONST	= 2.0;		// ADD TO TRKISO
const static double	kT_TRKISO_MF	= 0.005;	// MF TRKISO
const static double	kT_CESWIREE2_MF_LTETC18	= 0.14;
const static double	kT_CESWIREE2_CONST_GTETC18= 2.4;
const static double	kT_CESWIREE2_MF_GTETC18	= 0.01;

/*
function list in order

void Cleanup()

void TagCoversions()
int FiducialLepton()
InitForSashaCode()
SetPhotonList

CEMPhoEleIDcut(pho,event, zvx)
TightPhotonIDcut(pho)
LoosePhotonIDcut(pho)

FillDataBlock(ientry)


FillPhotonPlots(TSuperPhoton)
FillMetPlots

TFolder* GetHistoFolder
BookGeneralHistograms
BookPhotonHistograms
BookElectronHistograms

*/

bool SortTowersByEnergy(TCalTower* c1, TCalTower* c2)
{
	return (c1->Energy() > c2->Energy()) ? true : false;
}

ClassImp(TInit)
/*
need: _parPdgStat{ momPDG, mom Status, dau1 PDG, dau1 status , dau2 PDG, dau2 status}
return: parvec{mom 4-mom, dau1 4-mom, dau2 4-mom}


*/

bool FindInHepg(TGenpBlock* _fGenpBlock, std::vector<int> _parPdgStat,
					std::vector<TLorentzVector>& parvec, bool _Zmasscut)
{
	if (_fGenpBlock == NULL) {
		std::cout << "Err TGenpBlock null" << std::endl;
		return false;
	}

	if (_parPdgStat.size() != 6) {
		std::cout << "Err Not enough info!" << std::endl;
	} else {
		if ( _parPdgStat[1] < 1 || _parPdgStat[1] > 3 ||
			  _parPdgStat[3] < 1 || _parPdgStat[3] > 3 ||
		     _parPdgStat[5] < 1 || _parPdgStat[5] > 3 ) {
			
			std::cout << "Err invalid particle status!" << std::endl;
			return false;
		}
	}

	int Nparticles = _fGenpBlock->NParticles();

	TGenParticle *par, *mom, *dau;
	
	for (int i = 0 ; i < Nparticles ; i++) {
		par = _fGenpBlock->Particle(i);
		int im = par->GetFirstMother();

		if (im >= 0) {	//_________________________________im=-1 means no mother for imcoming particles
			mom = _fGenpBlock->Particle(im);

			if (mom != 0) {
				int par_id   = par->GetPdgCode();
				int par_stat = par->GetStatusCode();
				
				if ( (abs(par_id) == _parPdgStat[0]) && (par_stat == _parPdgStat[1]) ) { 	//found mom 
					TLorentzVector momvec;
					par->Momentum(momvec);
					parvec.push_back(momvec);
				
					if (_Zmasscut) { // Z mass window cut
						if (momvec.M() < 66. || momvec.M() > 116.) return false;
					}
					
					int p1st = par->GetFirstDaughter();
					int plast = par->GetLastDaughter();
					
					for (int i= p1st ; i <= plast || parvec.size() == 3; i++) {
						if (parvec.size() != 3) { // redundant ??????
							dau = _fGenpBlock->Particle(i);
							int dau_id = dau->GetPdgCode();
							int dau_stat = dau->GetStatusCode();
							
							TLorentzVector dauvec;
							dau->Momentum(dauvec);
							if ( ( (abs(dau_id) == _parPdgStat[2]) && (dau_stat == _parPdgStat[3]) ) ||
								  ( (abs(dau_id) == _parPdgStat[4]) && (dau_stat == _parPdgStat[5]) ) ) {
								
								parvec.push_back(dauvec);	
							}
							
						} else {
							break;
						}// if (parvec.size)
						
					} //for
					
					break;	//mom found, but did not decay in to expected daughters or failed Fid cut, so quit
				} 
			} //if mom
		} // if
	} //for
	
	if (parvec.size() != 3) {
		parvec.clear();
		return false;
	} else return true;

}


//_____________________________________________________________________________
TInit::TInit(const char* name, const char* title):
  TStnModule(name,title), qMc(false)
{
	std::cout << "Hello I am TInit module" << std::endl;
  fUseFindSpike=0; // 0--don't look for spikes; 1--look for spikes
  fUseSpikeFilter=0; // 0==do nothing; 1==remove spikes; -1==select spikes
  fMinEmPMTasym=0.65;
  fMinHadPMTasym=0.85; 
  fMaxEmPMTasym=1.0; 
  fMaxHadPMTasym=1.0;
  fMinEmPmtE=10.0; 
  fMinHadPmtE=10.0;
}

//_____________________________________________________________________________
TInit::~TInit() {
}

//_____________________________________________________________________________
void TInit::SaveHistograms() {
}

//_____________________________________________________________________________
void TInit::BookHistograms()
{
	DeleteHistograms();
	BookGeneralHistograms();
	BookElectronHistograms(fElectronHist, "Electrons");
	BookPhotonHistograms(fPhotonHist,"TightPhotons");
	BookPhotonHistograms(fPhoenixPhoHist,"PhoenixPhotons");
	BookPhotonHistograms(fPhoenixPhoPassPIDHist,"PhoenixPhoPassPID");
	BookCalHistograms(CalHist);
	BookBeamHaloStudyHistograms(fNoBeamHaloHist, "NoBeamHalo");
	BookBeamHaloStudyHistograms(fWithBeamHaloHist, "WithBeamHalo");
	BookBeamHaloStudyHistograms(fBH_passPIDHist, "BHpassPID");
	BookBeamHaloStudyHistograms(fBH_failPIDHist, "BHfailPID");
	BookBeamHaloStudyHistograms(fnotBH_passPIDHist, "noBHpassPID");
	BookBeamHaloStudyHistograms(fnotBH_failPIDHist, "noBHfailPID");
	BookTimingHistograms(fTimingHisto,"Cosmic");
	BookBeamHaloStudyHistograms(fCosmicHist, "Cosmic");
}


//_____________________________________________________________________________
int TInit::BeginJob()
{
				// register the data blocks

	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fJetBlock 		= (TStnJetBlock*)      RegisterDataBlock("JetBlock","TStnJetBlock");
	fTrackBlock 	= (TStnTrackBlock*)    RegisterDataBlock("TrackBlock","TStnTrackBlock");
	fPhotonBlock 	= (TStnPhotonBlock*)   RegisterDataBlock("PhotonBlock","TStnPhotonBlock");
	fTriggerBlock 	= (TStnTriggerBlock*)  RegisterDataBlock("TriggerBlock","TStnTriggerBlock");
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");
	fCalBlock 		= (TCalDataBlock*)     RegisterDataBlock("CalDataBlock","TCalDataBlock");
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");
	RegisterDataBlock("Phoenix_Electrons","TPhoenixElectronBlock",&fPhxEleBlock);
   RegisterDataBlock("PROD@PhoenixSI_Tracking","TStnTrackBlock",&fPhxSiTrackBlock);
	 


	BookHistograms();
	CreateBadTowerList();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.haloPassPhoId 		= 0;
	counter.haloCandidates 		= 0;
	counter.cosmicCandidates 	= 0;
	counter.cosmicPassPhoId 	= 0;
	counter.phoenixPhotons 		= 0;
	counter.phoenixPhoPassPID 	= 0;

	if (goodRunListFile.Length()<=0) goodRunListFile = "goodrun_v17_pho_02.txt";
	goodrun.Read(goodRunListFile.Data());

	trigbits.SetTriggerBlock(fTriggerBlock);

	return 0;
}

//_____________________________________________________________________________
int TInit::BeginRun()
{

	int currRun =  fHeaderBlock->RunNumber();
	std::cout << " BEGINING RUN# " << currRun << std::endl;
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	trigbits.BeginRun();	
	TStntuple::Init(currRun);
	
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TInit::Event(int ientry)
{
	SetPassed(0);
	counter.evtsRunOver++;
	//_____________________________ REQUIRED Thingies 
  FillDataBlocks(ientry);
  
  //if (fPhotonBlock->NPhotons() >0 && fVertexBlock->NVertices() == 0) {
  if (fPhotonBlock->NPhotons() >0) {
		//trigbits.Event();
    	//bool trig25 = trigbits.Pho25Iso();
	   //bool trig50 = trigbits.Pho50();
   	//bool trig70 = trigbits.Pho70();

	    //if (!(trig25 || trig50 || trig70) ) {
   	//	return 0;
	    //}
		SetPassed(1);
		counter.evtsPassModule++;
		Cleanup();
		
		double zvx  = 0;
		if (fVertexBlock->GetBestVertex(12,1) != NULL) {
			zvx = fVertexBlock->GetBestVertex(12,1)->Z();
		}
		TStnEvent* event = GetEvent();
		TStntuple::CorrectPhotonEnergy(event,zvx);	//	sam: std minor corrections!
		TPhotonUtil::CorrectPhotonIsolation(event);	// sam,09/01/07: ray's isolation corrections
		
		SetPhotonList();
		//TagConversions();
		//SpikeStudy();
		//HaloStudy();
		CosmicRayStudy();
		//PhoenixRejectionStudy();
	}

	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TInit::Cleanup()
{
	Nphotons =0;
	phoUncorrected.clear();
	phoCorrected.clear();
	phoHadEm.clear();
	phoEtaDet.clear();
	phoCesX.clear();
	phoCesZ.clear();
	phoIndex.clear();

	superpho.clear();


} //Cleanup


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TInit::SpikeStudy()
{

	CalDataArray cdarray;
	std::vector<int> parPdgStat;
	std::vector<TLorentzVector> parvec;
	bool ZmassCut = false;
	int looptimes = 0;

	if (qMc) {
		switch (fUseMcDataType) {
			case 0: 	//______ photon sample
				{
				std::cout << "SpikeStudy::NOT YET IMPLIMENTED. RETURNING." << std::endl;
				return;
				}
			case 1:	//______ Zee sample
				{
					parPdgStat.push_back(23); 	//Z
					parPdgStat.push_back(2);
					parPdgStat.push_back(11);	//electron
					parPdgStat.push_back(1);
					parPdgStat.push_back(11);
					parPdgStat.push_back(1);
					ZmassCut = true;
					looptimes = 2;
					break;
				}
			case 2:	//______ Wenu sample
				{
					std::cout << "SpikeStudy:: W  MC SAMPLE." << std::endl;
					parPdgStat.push_back(24);	//W
					parPdgStat.push_back(2);
					parPdgStat.push_back(11);	// electron
					parPdgStat.push_back(1);
					parPdgStat.push_back(12);//electron nutrino
					parPdgStat.push_back(1);
					looptimes = 1;
					break;
				}
			default :
				{
					std::cout << "SpikeStudy::I DON'T KNOW WHAT THIS MC SAMPLE IS!!!." << std::endl;
					return;
				}	
			
		}
		

		if (!FindInHepg(fGenpBlock,parPdgStat,parvec,ZmassCut)) return;

		
		for (; looptimes !=0; looptimes--) {	// loop for two daughters
			for (int i=0; i < superpho.size(); i++) {
				if (parvec[looptimes].DeltaR(superpho[i].corvec) < 0.2) {	
					double EmESpike=0,HadESpike=0,spikeHadTdc=0;
					MatchCalorTowers(superpho[i].index, fPhotonBlock, fCalBlock, &cdarray);
					FindPmtSpikes(&cdarray, EmESpike, HadESpike, spikeHadTdc, CalHist);
					
					if ( (EmESpike + HadESpike) > 0 ) {
						std::cout <<i << "# Em,Had, Tdc=" << EmESpike << "," << HadESpike 
							<< "," << spikeHadTdc <<std::endl;
					}
					cdarray.clear();
				}
			}
		}
				
	} else {

		for (int i=0; i < superpho.size(); i++) {
			//if (superpho[i].PhotonLikeID != 0) continue;	//____ pick only the electrons to test Spike removal

			double EmESpike=0,HadESpike=0,spikeHadTdc=0;
			MatchCalorTowers(superpho[i].index, fPhotonBlock, fCalBlock, &cdarray);
			FindPmtSpikes(&cdarray, EmESpike, HadESpike, spikeHadTdc, CalHist);
			
			if ( (EmESpike + HadESpike) > 0 ) 
				std::cout <<i << "# Em,Had, Tdc=" << EmESpike << "," << HadESpike 
							<< "," << spikeHadTdc <<std::endl;
			
			cdarray.clear();
		}
	}

}

/*-------------------------------------------------------------------*/
//___ finds sipkes, fills histo, returns spike Et
// loop over all towers occupies by the object,
// calculate pmt asym for each tower.
/*-------------------------------------------------------------------*/
int TInit::FindPmtSpikes(CalDataArray* towers, double& EmESpike,
								 double& HadESpike, double& spikeHadTdc,
								 CalHist_t& Hist) 
{
	int spike_code = 0;

	CalDataArrayI tti;
	EmESpike 	= 0.0;
	HadESpike 	= 0.0;
	spikeHadTdc = -1.0E6;
  
	for (tti = towers->begin(); tti != towers->end(); tti++) {
		TCalTower* calt = *tti;
      int iEta = calt->IEta();
      int iPhi = calt->IPhi();
      
		if (fabs(1.0 * towerInd(iEta)) > 11) continue; // cover all 2-pmt towers  //sam: pick only central towers

      int good_CEMtower	= 0; // should be ZERO for MY good towers  //sam? what is a good/bad tower
      int good_CHAtower	= 0; // should be ZERO for MY good towers
      int good_WHAtower	= 0; // should be ZERO for MY good towers
      float asymCEM		= -10.0;
      float asymCHAD		= -10.0;
      float asymWHAD		= -10.0;
      float pmt0			= 0.0;
      float pmt1			= 0.0;
      float cha_TDC		= -1.0E6;
      float wha_TDC		= -1.0E6;

      if (abs(towerInd(iEta)) <= 5) {	//_ central rowers
	  		pmt0 = calt->EmPmt(0);
	  		pmt1 = calt->EmPmt(1);
	  		
			if (pmt0+pmt1 > 200) asymCEM = (pmt0-pmt1)/(pmt0+pmt1);
	  
	  		pmt0 = calt->HadPmt(0);
	  		pmt1 = calt->HadPmt(1);
	  		cha_TDC = calt->HadTdc(0);
	  
	  		if (pmt0+pmt1 > 200) asymCHAD = (pmt0-pmt1)/(pmt0+pmt1);
		} //if
      
		if (abs(towerInd(iEta)) > 5 && abs(towerInd(iEta)) <= 7) {	//___ overlap region
	  		pmt0 = calt->EmPmt(0);
	  		pmt1 = calt->EmPmt(1);
	  		
			if (pmt0+pmt1 > 200) asymCEM = (pmt0-pmt1)/(pmt0+pmt1);
	  
	  		pmt0 = calt->HadPmt(0);
	  		pmt1 = calt->HadPmt(1);
	  		cha_TDC = calt->HadTdc(0);
	  
	  		if (pmt0+pmt1 > 200) asymCHAD = (pmt0-pmt1)/(pmt0+pmt1);
	  
	  		pmt0 = calt->HadPmt(2);
	  		pmt1 = calt->HadPmt(3);
	  		wha_TDC = calt->HadTdc(1);
	  		
			if (pmt0+pmt1 > 200) asymWHAD = (pmt0-pmt1)/(pmt0+pmt1);
		}
      
		if (abs(towerInd(iEta)) > 7 && abs(towerInd(iEta)) <=9) {	//___ overlap
	  		pmt0 = calt->EmPmt(0);
	  		pmt1 = calt->EmPmt(1);
	  
	  		if (pmt0+pmt1 > 200) asymCEM = (pmt0-pmt1)/(pmt0+pmt1);
	  		
			pmt0 = calt->HadPmt(0);
	  		pmt1 = calt->HadPmt(1);
	  		wha_TDC = calt->HadTdc(1);
	  		
			if (pmt0+pmt1 > 200) asymWHAD = (pmt0-pmt1)/(pmt0+pmt1);
		}
      
		if (abs(towerInd(iEta)) == 10) {
	  		pmt0 = calt->HadPmt(0);
	  		pmt1 = calt->HadPmt(1);
	  		wha_TDC = calt->HadTdc(1);
	  		
			if (pmt0+pmt1 > 200) asymWHAD = (pmt0-pmt1)/(pmt0+pmt1);
		}
      
		if (abs(towerInd(iEta)) == 11) {
			pmt0 = calt->HadPmt(2);
	  		pmt1 = calt->HadPmt(3);
	  		wha_TDC = calt->HadTdc(1);
	  
	  		if (pmt0+pmt1 > 200) asymWHAD = (pmt0-pmt1)/(pmt0+pmt1);
		}
		
      float EmEng		= calt->EmEnergy();
      float HadEng	= calt->HadEnergy();

      float twrEmEt 		= EmEng*fabs(sin(calt->Theta()*TMath::Pi()/180.0));
      float twrHadEt 	= HadEng*fabs(sin(calt->Theta()*TMath::Pi()/180.0));
      float twrEmPhi 	= calt->Phi();
      float twrHadPhi 	= calt->Phi();
      float dphimetEM 	= -1.0*TMath::Pi()/2.0;
      float dphimetHAD 	= -1.0*TMath::Pi()/2.0;
      float metEM			= -1.0;
      float metHAD		= -1.0;
      
		if (twrEmEt > 0.0) {
	  		metEM = fMetBlock->Met(0)/twrEmEt;			// MEt to Emt ratio and delPhi
	  		dphimetEM = fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(twrEmPhi) - TVector2::Phi_0_2pi(fMetBlock->MetPhi(0))));
		}
		
      if (twrHadEt>0.0) {
	  		metHAD=fMetBlock->Met(0)/twrHadEt;			// MEt to HadEt ratio and delPhi
	  		dphimetHAD=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(twrHadPhi)-TVector2::Phi_0_2pi(fMetBlock->MetPhi(0))));
		}

      if (fabs(asymCEM) <= 1.0) {	// why do we need this , asym is always <= 1 (=1 if one PMT is ZERO)
	  
	 		for (int ik = 0; ik < badEmTwr.size(); ik++) {
	      	if (towerInd(iEta) == badEmTwr[ik].etaI && iPhi == badEmTwr[ik].phiI) {
		  			good_CEMtower++;		//___ sam" this should be zero if all towers are good towers
		  			
					//sam: this is all bad tower stuff
					//if (EmEng>10.0) Hist.fCalPmtEmAssm_bad->Fill(asymCEM);
					
		  			//Hist.fCalPmtAssmVsRun[badEmTwr[ik].histoI]->Fill(GetHeaderBlock()->RunNumber(),asymCEM);
		  			
					if (fabs(asymCEM)>0.55 && EmEng>fMinEmPmtE) {
		      		//Hist.fCalBad_dPhi[badEmTwr[ik].histoI]->Fill(dphimetEM);
		      		//Hist.fCalBad_MetRatio[badEmTwr[ik].histoI]->Fill(metEM);
		    		}
		 			break;
				}
			} // for

			if (good_CEMtower == 0) {	
			
				if (EmEng<10.0) {
					Hist.fCalPmtEmAssm_le10->Fill(asymCEM); //sam: keep this, if object is central, this will be for central an so on
  					Hist.fCalPmtCemAssm_le10->Fill(asymCEM);
				} else {
					Hist.fCalPmtEmAssm_ge10->Fill(asymCEM);
		  			Hist.fCalPmtCemAssm_ge10->Fill(asymCEM);
				}
		
				Hist.fCalEmPmtAssm_vs_PmtSum->Fill(asymCEM,EmEng);
			} 
 	 
			if (fabs(asymCEM) > fMinEmPMTasym && fabs(asymCEM) <= fMaxEmPMTasym
					&& EmEng>fMinEmPmtE)
// 	  if(fabs(asymCEM)>0.5 && fabs(asymCEM)<=1.0 && EmEng>10.0)
//  	  if(fabs(asymCEM)>0.55 && fabs(asymCEM)<0.85 && EmEng>10.0)
//  	  if(asymCEM>0.55 && asymCEM<0.85 && EmEng>10.0)
//  	  if(asymCEM>-0.85 && asymCEM<-0.65 && EmEng>10.0)
			{
				spike_code++;
	      	EmESpike+=EmEng; // does it make sense to sum up energy for spikes?
	      	Hist.fCalEMtwr->Fill(towerInd(iEta),iPhi);
	      	Hist.fCalEm_dPhi->Fill(dphimetEM);
	      	Hist.fCalEm_MetRatio->Fill(metEM);
	    	}
		
		} // if asymCEM
      
		if (fabs(asymCHAD) <= 1.0) {
			
			for (int ik = 0; ik < badHadTwr.size(); ik++) {
	      	
				if (towerInd(iEta) == badHadTwr[ik].etaI && iPhi == badHadTwr[ik].phiI) {
		  			good_CHAtower++;
		  			//Hist.fCalPmtAssmVsRun[badHadTwr[ik].histoI]->Fill(GetHeaderBlock()->RunNumber(),asymCHAD);
		  			
					if ( (GetHeaderBlock()->RunNumber()) > badHadTwr[ik].run1 &&
						  (GetHeaderBlock()->RunNumber()) < badHadTwr[ik].run2 ) {
		     			
						//if (HadEng>10.0) Hist.fCalPmtHadAssm_bad->Fill(asymCHAD);
		      		
						if (fabs(asymCHAD)>0.55 && HadEng>fMinHadPmtE) {
			  				//Hist.fCalBad_dPhi[badHadTwr[ik].histoI]->Fill(dphimetHAD);
			  				//Hist.fCalBad_MetRatio[badHadTwr[ik].histoI]->Fill(metHAD);
						}
		    		}
					break;
				} 
	    	} // for 

			if (good_CHAtower == 0) {		
	      	if (HadEng<10.0) {		//___since we are looking at photons (electrons)
		  			
					Hist.fCalPmtHadAssm_le10->Fill(asymCHAD);
		  			if (abs(towerInd(iEta)) < 6) Hist.fCalPmtChaAssm_le10->Fill(asymCHAD);
		  			else Hist.fCalPmtChaAssm_overlap_le10->Fill(asymCHAD);
				
				} else {
		  		
					Hist.fCalPmtHadAssm_ge10->Fill(asymCHAD);
		  			
					if (abs(towerInd(iEta))<6) Hist.fCalPmtChaAssm_ge10->Fill(asymCHAD);
		  			else Hist.fCalPmtChaAssm_overlap_ge10->Fill(asymCHAD);
		  			
					if (fabs(cha_TDC)<40.0) Hist.fCalPmtHadAssm_ge10_intime->Fill(asymCHAD);
		  			if (fabs(cha_TDC)>40.0 && fabs(cha_TDC)<80.0) Hist.fCalPmtHadAssm_ge10_outtime->Fill(asymCHAD);
				}
	      	Hist.fCalHadPmtAssm_vs_PmtSum->Fill(asymCHAD,HadEng);
			} // if
 	  
			if (fabs(asymCHAD) > fMinHadPMTasym &&
				 fabs(asymCHAD) <= fMaxHadPMTasym && HadEng > fMinHadPmtE)
				// if(fabs(asymCHAD)>0.5 && fabs(asymCHAD)<=1.0 && HadEng>10.0)
	    		{
	      	spike_code++;
	      	HadESpike += HadEng; // does it make sense to sum up energy for spikes?
	      	spikeHadTdc=cha_TDC;
	      	Hist.fCalHadTDCspike->Fill(cha_TDC); // timing of HAD spikes;
	      	Hist.fCalCHAtwr->Fill(towerInd(iEta),iPhi);
	      	Hist.fCalHad_dPhi->Fill(dphimetHAD);
	      	Hist.fCalHad_MetRatio->Fill(metHAD);
	    	}
		} // if asymCHAD
      
		
		if (fabs(asymWHAD) <= 1.0) {
		
			for (int ik=0; ik<badHadTwr.size(); ik++) {
	      
				if (towerInd(iEta) == badHadTwr[ik].etaI && iPhi == badHadTwr[ik].phiI) {
		  			good_WHAtower++;
		  			//Hist.fCalPmtAssmVsRun[badHadTwr[ik].histoI]->Fill(GetHeaderBlock()->RunNumber(),asymWHAD);
		  
		  			if ( (GetHeaderBlock()->RunNumber()) > badHadTwr[ik].run1 && 
						  (GetHeaderBlock()->RunNumber()) < badHadTwr[ik].run2 ) {
		      
						//if (HadEng>10.0) Hist.fCalPmtHadAssm_bad->Fill(asymWHAD);
		      
						if (fabs(asymWHAD)>0.55 && HadEng>fMinHadPmtE) {
			  				//Hist.fCalBad_dPhi[badHadTwr[ik].histoI]->Fill(dphimetHAD);
			  				//Hist.fCalBad_MetRatio[badHadTwr[ik].histoI]->Fill(metHAD);
						}
		    		}
		  			break;
				}
			} // for

			if (good_WHAtower == 0) {
	      	if (HadEng<10.0) {
		 			Hist.fCalPmtHadAssm_le10->Fill(asymWHAD);
		  			
					if (abs(towerInd(iEta)) > 7 && abs(towerInd(iEta)) < 11) {
						Hist.fCalPmtWhaAssm_le10->Fill(asymWHAD);
		  			} else {
						Hist.fCalPmtWhaAssm_overlap_le10->Fill(asymWHAD);
					}
			
				} else {
		  		
					Hist.fCalPmtHadAssm_ge10->Fill(asymWHAD);
		  			if (abs(towerInd(iEta)) > 7 && abs(towerInd(iEta)) < 11) {
						Hist.fCalPmtWhaAssm_ge10->Fill(asymWHAD);
		  			} else {
						Hist.fCalPmtWhaAssm_overlap_ge10->Fill(asymWHAD);
					}
		  
		  			if (fabs(wha_TDC) < 40.0) {
						Hist.fCalPmtHadAssm_ge10_intime->Fill(asymWHAD);
		  			}
				
					if (fabs(wha_TDC) > 40.0 && fabs(wha_TDC) < 80.0) {
						Hist.fCalPmtHadAssm_ge10_outtime->Fill(asymWHAD);
					}
				}
	      	Hist.fCalHadPmtAssm_vs_PmtSum->Fill(asymWHAD,HadEng);
			}
	  	
			if (fabs(asymWHAD) > fMinHadPMTasym &&
				 fabs(asymWHAD) <= fMaxHadPMTasym && 
				 HadEng > fMinHadPmtE ) {

				spike_code++;
				HadESpike += HadEng; // does it make sense to sum up energy for spikes?
				spikeHadTdc = wha_TDC;
				Hist.fCalHadTDCspike->Fill(wha_TDC); // timing of HAD spikes;
				Hist.fCalWHAtwr->Fill(towerInd(iEta),iPhi);
				Hist.fCalHad_dPhi->Fill(dphimetHAD);
				Hist.fCalHad_MetRatio->Fill(metHAD);
			}
		}
   }

  return spike_code;
}


/*-------------------------------------------------------------------*/
// ___________ create a list of bad EM/HAD towers
/*-------------------------------------------------------------------*/
void TInit::CreateBadTowerList() 
{
	//_____________________________ creating a list of BAD EM towers
	BadTower _twr;
	_twr.phiI=0;
	_twr.etaI=2;
	_twr.histoI=0;
	_twr.run1=-1;
	_twr.run2=-1;
	badEmTwr.push_back(_twr);
	_twr.phiI=4;
	_twr.etaI=8;
	_twr.histoI=1;
	_twr.run1=-1;
	_twr.run2=-1;
	badEmTwr.push_back(_twr);
	_twr.phiI=4;
	_twr.etaI=7;
	_twr.histoI=2;
	_twr.run1=-1;
	_twr.run2=-1;
	badEmTwr.push_back(_twr);
	_twr.phiI=23;
	_twr.etaI=3;
	_twr.histoI=3;
	_twr.run1=-1;
	_twr.run2=-1;
	badEmTwr.push_back(_twr);
	//_____________________________ creating a list of BAD HAD towers
	_twr.phiI=2;
	_twr.etaI=-11;
	_twr.histoI=4;
	_twr.run1=182800;
	_twr.run2=186600;
	badHadTwr.push_back(_twr);
	_twr.phiI=0;
	_twr.etaI=-10;
	_twr.histoI=5;
	_twr.run1=164600;
	_twr.run2=169000;
	badHadTwr.push_back(_twr);
	_twr.phiI=23;
	_twr.etaI=-1;
	_twr.histoI=6;
	_twr.run1=138500;
	_twr.run2=156500;
	badHadTwr.push_back(_twr);
	_twr.phiI=5;
	_twr.etaI=7;
	_twr.histoI=7;
	_twr.run1=138000;
	_twr.run2=188000;
	badHadTwr.push_back(_twr);
	_twr.phiI=5;
	_twr.etaI=6;
	_twr.histoI=8;
	_twr.run1=138000;
	_twr.run2=188000;
	badHadTwr.push_back(_twr);
	_twr.phiI=22;
	_twr.etaI=8;
	_twr.histoI=9;
	_twr.run1=160500;
	_twr.run2=164500;
	badHadTwr.push_back(_twr);
	_twr.phiI=1;
	_twr.etaI=9;
	_twr.histoI=10;
	_twr.run1=152400;
	_twr.run2=154400;
	badHadTwr.push_back(_twr);
	_twr.phiI=12;
	_twr.etaI=10;
	_twr.histoI=11;
	_twr.run1=151200;
	_twr.run2=156600;
	badHadTwr.push_back(_twr);
	
}


/*-------------------------------------------------------------------*/
//________ creates a list of towers for all jets
/*-------------------------------------------------------------------*/
void TInit::MatchCalorTowers(TStnJet* jet,  TStnJetBlock* fJetBlock,
				TCalDataBlock *fCalDataBlock, CalDataArray* cdholder) {
/*{{{*/
  cdholder->clear();
  TStnLinkBlock* links = fJetBlock->TowerLinkList();
  int nptow = links->NLinks(jet->Number());

  TCalTower* ctower = NULL;
  for(int j=0; j<nptow; j++)
    {
      int iptow = links->Index(jet->Number(),j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;

      ctower = fCalDataBlock->Tower(ieta,iphi);
      cdholder->push_back(ctower);
    }

  std::sort(cdholder->begin(), cdholder->end(),SortTowersByEnergy);

  return;
}
/*}}}*/

/*-------------------------------------------------------------------*/
//________ creates a list of towers for all EM objects
/*-------------------------------------------------------------------*/
void TInit::MatchCalorTowers(int pho_ind, TStnPhotonBlock* fPhotonBlock,
				TCalDataBlock *fCalDataBlock, CalDataArray* cdholder)
/*{{{*/
{
  cdholder->clear();
  TStnLinkBlock* links = fPhotonBlock->TowerLinkList();
  int nptow = links->NLinks(pho_ind);
  TCalTower* ctower = NULL;
  for(int j=0; j<nptow; j++)
    {
      int iptow = links->Index(pho_ind,j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;

      ctower = fCalDataBlock->Tower(ieta,iphi);
      cdholder->push_back(ctower);
    }

  std::sort(cdholder->begin(),cdholder->end(),SortTowersByEnergy);

  return;
}
/*}}}*/


void TInit::BookCalHistograms(CalHist_t& Hist) {
/*{{{*/
	
	TFolder *new_folder = GetHistoFolder("PMTSpikes", "PMTSpikes");

  char name [200];
  char title[200];

  sprintf(name,"CalEMtwr");
  sprintf(title," EM towers with spikes: iPhi vs Tower Index");
  Hist.fCalEMtwr=new TH2F(name,title,23,-11.5,11.5,24,-0.5,23.5);
  Hist.fCalEMtwr->SetMarkerStyle(kPlus);
  new_folder->Add(Hist.fCalEMtwr);
 
  sprintf(name,"CalCHAtwr");
  sprintf(title," CHA towers with spikes: iPhi vs Tower Index");
  Hist.fCalCHAtwr=new TH2F(name,title,23,-11.5,11.5,24,-0.5,23.5);
  Hist.fCalCHAtwr->SetMarkerStyle(kPlus);
  new_folder->Add(Hist.fCalCHAtwr);
  
  sprintf(name,"CalWHAtwr");
  sprintf(title," WHA towers with spikes: iPhi vs Tower Index");
  Hist.fCalWHAtwr=new TH2F(name,title,23,-11.5,11.5,24,-0.5,23.5);
  Hist.fCalWHAtwr->SetMarkerStyle(kPlus);
  new_folder->Add(Hist.fCalWHAtwr);

  sprintf(name,"CalHadTDCspike");
  sprintf(title,"Timing of HAD spikes (CHA & WHA)");
  Hist.fCalHadTDCspike=new TH1F(name,title,1000,-200.0,200.0);
  new_folder->Add(Hist.fCalHadTDCspike);

	//____________________________________ bad tower studies
/*{{{*/
/*
  sprintf(name,"CalPmtAssmVsRun[0]");
  sprintf(title,"CEM tower (0,2): (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[0]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[0]);
  sprintf(name,"CalPmtAssmVsRun[1]");
  sprintf(title,"CEM tower (4,8): (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[1]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[1]);
  sprintf(name,"CalPmtAssmVsRun[2]");
  sprintf(title,"CEM tower (4,7): (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[2]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[2]);
  sprintf(name,"CalPmtAssmVsRun[3]");
  sprintf(title,"CEM tower (23,3): (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[3]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[3]);
  sprintf(name,"CalPmtAssmVsRun[4]");
  sprintf(title,"WHA tower (2,-11): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[4]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[4]);
  sprintf(name,"CalPmtAssmVsRun[5]");
  sprintf(title,"WHA tower (0,-10): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[5]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[5]);
  sprintf(name,"CalPmtAssmVsRun[6]");
  sprintf(title,"CHA tower (23,-1): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[6]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[6]);
  sprintf(name,"CalPmtAssmVsRun[7]");
  sprintf(title,"CHA tower (5,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[7]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[7]);
  sprintf(name,"CalPmtAssmVsRun[8]");
  sprintf(title,"CHA tower (5,6): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[8]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[8]);
  sprintf(name,"CalPmtAssmVsRun[9]");
  sprintf(title,"WHA tower (22,8): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[9]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[9]);
  sprintf(name,"CalPmtAssmVsRun[10]");
  sprintf(title,"WHA tower (1,9): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[10]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[10]);
  sprintf(name,"CalPmtAssmVsRun[11]");
  sprintf(title,"WHA tower (12,10): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[11]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[11]);
  sprintf(name,"CalPmtAssmVsRun[12]");
  sprintf(title,"CHA tower (11,-7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[12]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[12]);
  sprintf(name,"CalPmtAssmVsRun[13]");
  sprintf(title,"CHA tower (12,-7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[13]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[13]);
  sprintf(name,"CalPmtAssmVsRun[14]");
  sprintf(title,"CHA tower (7,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[14]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[14]);
  sprintf(name,"CalPmtAssmVsRun[15]");
  sprintf(title,"CHA tower (9,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[15]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[15]);
  sprintf(name,"CalPmtAssmVsRun[16]");
  sprintf(title,"CHA tower (10,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[16]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[16]);
  sprintf(name,"CalPmtAssmVsRun[17]");
  sprintf(title,"CHA tower (14,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[17]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[17]);
  sprintf(name,"CalPmtAssmVsRun[18]");
  sprintf(title,"CHA tower (19,7): (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) vs. Run");
  Hist.fCalPmtAssmVsRun[18]=new TProfile(name,title,500,130000.0,230000.0,-1.01,1.01);
  new_folder->Add(Hist.fCalPmtAssmVsRun[18]);
*/
/*}}}*/

/*{{{*/
/*
  sprintf(name,"CalBad_dPhi[0]");
  sprintf(title,"CEM tower (0,2): d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[0]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[0]);
  sprintf(name,"CalBad_MetRatio[0]");
  sprintf(title,"CEM tower (0,2): #slash{E}_{T}/E_{T}^{EM spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[0]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[0]);
  sprintf(name,"CalBad_dPhi[1]");
  sprintf(title,"CEM tower (4,8): d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[1]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[1]);
  sprintf(name,"CalBad_MetRatio[1]");
  sprintf(title,"CEM tower (4,8): #slash{E}_{T}/E_{T}^{EM spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[1]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[1]);
  sprintf(name,"CalBad_dPhi[2]");
  sprintf(title,"CEM tower (4,7): d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[2]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[2]);
  sprintf(name,"CalBad_MetRatio[2]");
  sprintf(title,"CEM tower (4,7): #slash{E}_{T}/E_{T}^{EM spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[2]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[2]);
  sprintf(name,"CalBad_dPhi[3]");
  sprintf(title,"CEM tower (23,3): d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[3]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[3]);
  sprintf(name,"CalBad_MetRatio[3]");
  sprintf(title,"CEM tower (23,3): #slash{E}_{T}/E_{T}^{EM spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[3]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[3]);
  sprintf(name,"CalBad_dPhi[4]");
  sprintf(title,"WHA tower (2,-11): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[4]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[4]);
  sprintf(name,"CalBad_MetRatio[4]");
  sprintf(title,"WHA tower (2,-11): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[4]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[4]);
  sprintf(name,"CalBad_dPhi[5]");
  sprintf(title,"WHA tower (0,-10): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[5]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[5]);
  sprintf(name,"CalBad_MetRatio[5]");
  sprintf(title,"WHA tower (0,-10): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[5]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[5]);
  sprintf(name,"CalBad_dPhi[6]");
  sprintf(title,"CHA tower (23,-1): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[6]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[6]);
  sprintf(name,"CalBad_MetRatio[6]");
  sprintf(title,"CHA tower (23,-1): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[6]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[6]);
  sprintf(name,"CalBad_dPhi[7]");
  sprintf(title,"CHA tower (5,7): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[7]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[7]);
  sprintf(name,"CalBad_MetRatio[7]");
  sprintf(title,"CHA tower (5,7): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[7]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[7]);
  sprintf(name,"CalBad_dPhi[8]");
  sprintf(title,"CHA tower (5,6): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[8]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[8]);
  sprintf(name,"CalBad_MetRatio[8]");
  sprintf(title,"CHA tower (5,6): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[8]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[8]);
  sprintf(name,"CalBad_dPhi[9]");
  sprintf(title,"WHA tower (22,8): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[9]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[9]);
  sprintf(name,"CalBad_MetRatio[9]");
  sprintf(title,"WHA tower (22,8): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[9]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[9]);
  sprintf(name,"CalBad_dPhi[10]");
  sprintf(title,"WHA tower (1,9): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[10]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[10]);
  sprintf(name,"CalBad_MetRatio[10]");
  sprintf(title,"WHA tower (1,9): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[10]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[10]);
  sprintf(name,"CalBad_dPhi[11]");
  sprintf(title,"WHA tower (12,10): d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_dPhi[11]=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalBad_dPhi[11]);
  sprintf(name,"CalBad_MetRatio[11]");
  sprintf(title,"WHA tower (12,10): #slash{E}_{T}/E_{T}^{HAD spike} if (Q_{pmt1}-Q_{pmt2})/(Q_{pmt1}+Q_{pmt2})>0.55");
  Hist.fCalBad_MetRatio[11]=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalBad_MetRatio[11]);
*/
/*}}}*/

   //___________________________________________________________________________________________________________

  sprintf(name,"CalEm_dPhi");
  sprintf(title,"d#phi=#phi_{EM spike}-#phi_{#slash{E}_{T}} if E_{em}^{twr}>10 GeV");
  Hist.fCalEm_dPhi=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalEm_dPhi);

  sprintf(name,"CalEm_MetRatio");
  sprintf(title,"#slash{E}_{T}/E_{T}^{EM spike} if E_{em}^{twr}>10 GeV");
  Hist.fCalEm_MetRatio=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalEm_MetRatio);
  
  sprintf(name,"CalHad_dPhi");
  sprintf(title,"d#phi=#phi_{HAD spike}-#phi_{#slash{E}_{T}} if E_{had}^{twr}>10 GeV");
  Hist.fCalHad_dPhi=new TH1F(name,title,100,0.0,3.15);
  new_folder->Add(Hist.fCalHad_dPhi);
  
  sprintf(name,"CalHad_MetRatio");
  sprintf(title,"#slash{E}_{T}/E_{T}^{HAD spike} if E_{had}^{twr}>10 GeV");
  Hist.fCalHad_MetRatio=new TH1F(name,title,500,0.0,10.0);
  new_folder->Add(Hist.fCalHad_MetRatio);


/*  sprintf(name,"CalPmtEmAssm_bad");
  sprintf(title,"pmt asymmetry for bad EM towers: (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2})");
  Hist.fCalPmtEmAssm_bad=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtEmAssm_bad);
  sprintf(name,"CalPmtHadAssm_bad");
  sprintf(title,"pmt asymmetry for bad HAD towers: (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2})");
  Hist.fCalPmtHadAssm_bad=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_bad);
*/

  sprintf(name,"CalPmtEmAssm_le10");
  sprintf(title,"pmt asymmetry: (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) if E_{em}^{twr}<10 GeV");
  Hist.fCalPmtEmAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtEmAssm_le10);
  sprintf(name,"CalPmtHadAssm_le10");
  sprintf(title,"pmt asymmetry: (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) if E_{had}^{twr}<10 GeV");
  Hist.fCalPmtHadAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_le10);
  sprintf(name,"CalPmtEmAssm_ge10");
  sprintf(title,"EM pmt asymmetry: (E_{em}^{pmt1}-E_{em}^{pmt2})/(E_{em}^{pmt1}+E_{em}^{pmt2}) if E_{em}^{twr}>10 GeV");
  Hist.fCalPmtEmAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtEmAssm_ge10);
  sprintf(name,"CalPmtHadAssm_ge10");
  sprintf(title,"HAD pmt asymmetry: (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) if E_{had}^{twr}>10 GeV");
  Hist.fCalPmtHadAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_ge10);

  sprintf(name,"CalPmtCemAssm_le10");
  sprintf(title,"CEM pmt asymmetry: (E_{cem}^{pmt1}-E_{cem}^{pmt2})/(E_{cem}^{pmt1}+E_{cem}^{pmt2}) if E_{cem}^{twr}<10 GeV");
  Hist.fCalPmtCemAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtCemAssm_le10);
  sprintf(name,"CalPmtChaAssm_le10");
  sprintf(title,"CHA pmt asymmetry outside overlap: (E_{cha}^{pmt1}-E_{cha}^{pmt2})/(E_{cha}^{pmt1}+E_{cha}^{pmt2}) if E_{cha}^{twr}<10 GeV");
  Hist.fCalPmtChaAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtChaAssm_le10);
  sprintf(name,"CalPmtWhaAssm_le10");
  sprintf(title,"WHA pmt asymmetry outside overlap: (E_{wha}^{pmt1}-E_{wha}^{pmt2})/(E_{wha}^{pmt1}+E_{wha}^{pmt2}) if E_{wha}^{twr}<10 GeV");
  Hist.fCalPmtWhaAssm_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtWhaAssm_le10);

  sprintf(name,"CalPmtCemAssm_ge10");
  sprintf(title,"CEM pmt asymmetry: (E_{cem}^{pmt1}-E_{cem}^{pmt2})/(E_{cem}^{pmt1}+E_{cem}^{pmt2}) if E_{cem}^{twr}>10 GeV");
  Hist.fCalPmtCemAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtCemAssm_ge10);
  sprintf(name,"CalPmtChaAssm_ge10");
  sprintf(title,"CHA pmt asymmetry outside overlap: (E_{cha}^{pmt1}-E_{cha}^{pmt2})/(E_{cha}^{pmt1}+E_{cha}^{pmt2}) if E_{cha}^{twr}>10 GeV");
  Hist.fCalPmtChaAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtChaAssm_ge10);
  sprintf(name,"CalPmtWhaAssm_ge10");
  sprintf(title,"WHA pmt asymmetry outside overlap: (E_{wha}^{pmt1}-E_{wha}^{pmt2})/(E_{wha}^{pmt1}+E_{wha}^{pmt2}) if E_{wha}^{twr}>10 GeV");
  Hist.fCalPmtWhaAssm_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtWhaAssm_ge10);

  sprintf(name,"CalPmtChaAssm_overlap_le10");
  sprintf(title,"CHA pmt asymmetry in overlap region: (E_{cha}^{pmt1}-E_{cha}^{pmt2})/(E_{cha}^{pmt1}+E_{cha}^{pmt2}) if E_{had}^{twr}<10 GeV");
  Hist.fCalPmtChaAssm_overlap_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtChaAssm_overlap_le10);
  sprintf(name,"CalPmtWhaAssm_overlap_le10");
  sprintf(title,"WHA pmt asymmetry in overlap region: (E_{wha}^{pmt1}-E_{wha}^{pmt2})/(E_{wha}^{pmt1}+E_{wha}^{pmt2}) if E_{had}^{twr}<10 GeV");
  Hist.fCalPmtWhaAssm_overlap_le10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtWhaAssm_overlap_le10);
  sprintf(name,"CalPmtChaAssm_overlap_ge10");
  sprintf(title,"CHA pmt asymmetry in overlap region: (E_{cha}^{pmt1}-E_{cha}^{pmt2})/(E_{cha}^{pmt1}+E_{cha}^{pmt2}) if E_{had}^{twr}>10 GeV");
  Hist.fCalPmtChaAssm_overlap_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtChaAssm_overlap_ge10);
  sprintf(name,"CalPmtWhaAssm_overlap_ge10");
  sprintf(title,"WHA pmt asymmetry in overlap region: (E_{wha}^{pmt1}-E_{wha}^{pmt2})/(E_{wha}^{pmt1}+E_{wha}^{pmt2}) if E_{had}^{twr}>10 GeV");
  Hist.fCalPmtWhaAssm_overlap_ge10=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtWhaAssm_overlap_ge10);

  sprintf(name,"CalPmtHadAssm_ge10_intime");
  sprintf(title,"HAD pmt asymmetry: (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) if E_{had}^{twr}>10 GeV and |HadTdc|<40 ns");
  Hist.fCalPmtHadAssm_ge10_intime=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_ge10_intime);
  sprintf(name,"CalPmtHadAssm_ge10_outtime");
  sprintf(title,"HAD pmt asymmetry: (E_{had}^{pmt1}-E_{had}^{pmt2})/(E_{had}^{pmt1}+E_{had}^{pmt2}) if E_{had}^{twr}>10 GeV and 40<|HadTdc|<80 ns");
  Hist.fCalPmtHadAssm_ge10_outtime=new TH1F(name,title,200,-1.0,1.0);
  new_folder->Add(Hist.fCalPmtHadAssm_ge10_outtime);

  sprintf(name,"CalEmPmtAssm_vs_PmtSum");
  sprintf(title,"CEM Tower Energy vs CEM Tower PMT Asymmetry");
  Hist.fCalEmPmtAssm_vs_PmtSum=new TH2F(name,title,100,-1.0,1.0,300,0.0,300.0);
  new_folder->Add(Hist.fCalEmPmtAssm_vs_PmtSum);
  
  sprintf(name,"CalHadPmtAssm_vs_PmtSum");
  sprintf(title,"Had (CHA & WHA) Tower Energy vs Had Tower PMT Asymmetry");
  Hist.fCalHadPmtAssm_vs_PmtSum=new TH2F(name,title,100,-1.0,1.0,300,0.0,300.0);
  new_folder->Add(Hist.fCalHadPmtAssm_vs_PmtSum);

  return;
}
/*}}}*/



/*========================================== BH STUFF =============*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TInit::HaloStudy()
{
	int N_beamhalo = 0;
	float haloIso_cut = 0.;
	int my_nvx12 = Class12Vertices();

	for (int i=0; i< superpho.size(); i++) {
			
		if (superpho[i].EtCorr < 20.0) haloIso_cut = 0.15 * superpho[i].EtCorr;
      else haloIso_cut = 3.0 + 0.02 * (superpho[i].EtCorr - 20.0);
			
      float haloIso = superpho[i].IsoEtCorr;
		
      if ( (superpho[i].Detector == 0) && // minimum quality cuts for beam halo candidates
 	 		  (superpho[i].HadEm < 0.125) && 
	 		  (haloIso < haloIso_cut    ) ) {

			DoMyBeamHalo(superpho[i].index);
			
			if (BeamHaloCut(&myHaloStuff,my_nvx12,HaloCutScenario) == 0) {
				superpho[i].BeamHalo = true;
				
				FillBeamHaloStudyHistograms(fNoBeamHaloHist,&myHaloStuff,superpho[i].pho);
				if (superpho[i].TightID == 0 && superpho[i].PhotonLikeID != 0 ) {
					FillBeamHaloStudyHistograms(fnotBH_passPIDHist,&myHaloStuff,superpho[i].pho);
				} else {
					FillBeamHaloStudyHistograms(fnotBH_failPIDHist,&myHaloStuff,superpho[i].pho);
				}
				
	    	} else {
				
	      	N_beamhalo++;
				counter.haloCandidates++;
	      	FillBeamHaloStudyHistograms(fWithBeamHaloHist,&myHaloStuff,superpho[i].pho);
				
				if (superpho[i].TightID == 0 && superpho[i].PhotonLikeID != 0 ) {
					counter.haloPassPhoId++;
					FillBeamHaloStudyHistograms(fBH_passPIDHist,&myHaloStuff,superpho[i].pho);
				} else {
					FillBeamHaloStudyHistograms(fBH_failPIDHist,&myHaloStuff,superpho[i].pho);
				}

			}
		}
	}
	//if (N_beamhalo >0)
	//std::cout << "N_beamhalo =" << N_beamhalo << std::endl;

}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TInit::DoMyBeamHalo(int pho_ind)
{
  int _seedWedge	= 0;		//EM 
  int _sideWedge	= 0;
  int _seedWedgeH	= 0;		//HAD
  int _sideWedgeH	= 0;
  int _eastNHads	= 0;		
  int _westNHads	= 0;
  double _emSeedE	= 0.0; 
  double _emSideE	= 0.0; 
  double _hadSeedE= 0.0;
  double _hadSeedWedgeE = 0.0;

  getBeamHaloInfo(pho_ind,_seedWedge,_seedWedgeH,_sideWedge,_sideWedgeH,
		  _eastNHads,_westNHads,_emSeedE,_emSideE,_hadSeedE,_hadSeedWedgeE);
  
  myHaloStuff.seedwedgeH=_seedWedgeH;
  myHaloStuff.sidewedgeH=_sideWedgeH;
  myHaloStuff.seedwedge=_seedWedge;
  myHaloStuff.sidewedge=_sideWedge;
  myHaloStuff.eastNhad=_eastNHads;
  myHaloStuff.westNhad=_westNHads;
  myHaloStuff.emseedE=_emSeedE;
  myHaloStuff.emsideE=_emSideE;
  myHaloStuff.hadseedE=_hadSeedE;
  myHaloStuff.seedwedgeHadE=_hadSeedWedgeE;
 
  return;
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TInit::getBeamHaloInfo(int pho_ind, int& seedWedge, int& seedWedgeH,
					    int& sideWedge, int& sideWedgeH, int& eastNHads,
					    int& westNHads, double& emSeedE, double& emSideE, 
					    double& hadSeedE, double& seedWedgeHadE)
{
	CalDataArray calholder;
	MatchCalorTowers(pho_ind, &calholder);
	TCalTower* tower1 = calholder[0];
	int iPhiEmO = tower1->IPhi();

	float thresh 	= 0.1; //sam: EM energy threshold (0.1GeV)
	float threshH 	= 0.1; // added by Sasha : sam: HAD energy threshold (0.1 GeV)
	seedWedge 		= 0;
	sideWedge 		= 0;
	seedWedgeH 		= 0;  // added by Sasha
	sideWedgeH 		= 0;  // added by Sasha
	eastNHads 		= 0;
	westNHads 		= 0;
	emSeedE 			= 0.;
	emSideE 			= 0.;
	hadSeedE 		= 0.;
	seedWedgeHadE	= 0.0;
	
			//		iEta=12 to iEta< 40						
	for (int iEta = TENETA/2 - 14; iEta < TENETA/2 + 14; iEta++ ) {	//TENETA=52

      int iPhi = (TETTYP[iEta] == 5) ? iPhiEmO * 2 : iPhiEmO;
      TCalTower* c1 = fCalBlock->Tower(iEta,iPhi);

      if (iEta - TENETA/2 > -11 && iEta - TENETA/2 < 10) {  //_ sam: are we looking only at a specific region in the detecctor?
			if (c1 == NULL || c1->EmEnergy() < thresh) continue;
			seedWedge++;		//___________________________________total number of seed wedges
			emSeedE = emSeedE + c1->EmEnergy();	//________________ total ME energy in the seed wedges
			seedWedgeHadE = seedWedgeHadE + c1->HadEnergy(); //___ total HAD energy in the seed towers
			
			if (c1->HadEnergy() > threshH) seedWedgeH++;  //______ total number of had towers
			
		} else {
			//iEta < 12 and iEta >= 40
			TCalTower* c2 = (TETTYP[iEta] == 5) ? fCalBlock->Tower(iEta,iPhi+1) : NULL;
			bool pass = false;
			if (c1 != NULL && c1->HadEnergy() > thresh) pass = true;
			if (c2 != NULL && c2->HadEnergy() > thresh) pass = true;
			if (!pass) continue;

			if (c1) hadSeedE = hadSeedE + c1->HadEnergy();
			if (c2) hadSeedE = hadSeedE + c2->HadEnergy();

			if(iEta < TENETA/2) westNHads++;
			if(iEta > TENETA/2) eastNHads++;
		}
	}

  int iphi1 = iPhiEmO==0 ? 23 : iPhiEmO-1;
  int iphi2 = iPhiEmO==23 ? 0 : iPhiEmO+1;

	for (int iEta = TENETA/2-10; iEta < TENETA/2 + 10; iEta++ ) {
   	TCalTower* ctower1 = fCalBlock->Tower(iEta,iphi1);
		TCalTower* ctower2 = fCalBlock->Tower(iEta,iphi2);
		if (ctower1 != NULL && ctower1->EmEnergy() > thresh) {
			sideWedge++;
			emSideE = emSideE + ctower1->EmEnergy();
			if (ctower1->HadEnergy() > threshH) sideWedgeH++; // added by Sasha;
		}
      
		if(ctower2 != NULL && ctower2->Energy() > thresh) {
			sideWedge++;
			emSideE = emSideE + ctower2->EmEnergy();
			if (ctower2->HadEnergy() > threshH) sideWedgeH++; // added by Sasha;
		}
	}
  
  return;
}

/*-------------------------------------------------------------------*/
//default scenario=0
// we are not correcting for MI and UI. but we may want to later on.
// so leave the nvx12 in.
// returns 1 if 
/*-------------------------------------------------------------------*/

int TInit::BeamHaloCut(HaloStuff *beamhalo, int nvx12, int scenario)
{
  	//double ue=0.013; // underlying event EM+HAD energy per tower
	//double mi=0.019; // multiple interaction EM+HAD energy per tower
	  
	if (scenario == 0) { // rejection~100%; misID=1.6%
      if ( (beamhalo->seedwedge) > 8 ||
			  ((beamhalo->eastNhad + beamhalo->westNhad) > 2) ) return 1;
	}
  
	if (scenario == 1) { // rejection=99%; misID=0.3%
 		if ( (beamhalo->seedwedge) > 4 && 
			  ((beamhalo->eastNhad + beamhalo->westNhad) > 1) ) return 1; 
	}
  
	if (scenario == 2) { // rejection=98%; misID=0.05%
		if ( (beamhalo->seedwedge) > 4 && 
			  ((beamhalo->eastNhad + beamhalo->westNhad) > 2) ) return 1; 
    }
	 
	if (scenario == 3) { // rejection=95%; misID~0.006%
      if ( (beamhalo->seedwedge) > 7 && 
			  ((beamhalo->eastNhad + beamhalo->westNhad) > 2) ) return 1; 
	}
  
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

int TInit::Class12Vertices()
{
 	int nvtx = 0;
	int Nvtx = fVertexBlock->NVertices();

	for (int ivtx = 0; ivtx < Nvtx; ++ivtx) {
		TStnVertex* vert = fVertexBlock->Vertex(ivtx);
		if (vert->VClass() >= 12) ++nvtx;
	}
  
  return nvtx;

}  // Class12Vertices

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TInit::MatchCalorTowers(int pho_ind, CalDataArray* cdholder)
{
  cdholder->clear();
  TStnLinkBlock* links = fPhotonBlock->TowerLinkList();
  int nptow = links->NLinks(pho_ind);

  TCalTower* ctower = NULL;
  for(int j=0; j<nptow; j++) 
    {
      int iptow = links->Index(pho_ind,j);
      int iphi = iptow & 0x3F;
      int ieta = (iptow>>8) & 0x3F;

      ctower = fCalBlock->Tower(ieta,iphi);
      cdholder->push_back(ctower);
    }

  std::sort(cdholder->begin(), cdholder->end(), SortTowersByEnergy);
  return;
}




/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TInit::BookBeamHaloStudyHistograms(BeamHaloStudyHisto_t& Hist, char* Folder)
{

	TFolder* new_folder = GetHistoFolder(Folder,Folder);

  char name [200];
  char title[200];

  sprintf(name,"HaloSeedWedgeEM_SeedWedgeHAD");
  sprintf(title,"N_{EMtwr} vs. N_{HADtwr} in seed wedge");
  Hist.fHaloSeedWedgeEM_SeedWedgeHAD = new TH2F(name,title,25,-0.5,24.5,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloSeedWedgeEM_SeedWedgeHAD);

  sprintf(name,"HaloSeedWedge");
  sprintf(title,"number of EM towers in seed wedge, Max's version");
  Hist.fHaloSeedWedge=new TH1F(name,title,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloSeedWedge);
  sprintf(name,"HaloSideWedge");
  sprintf(title,"number of EM towers in side wedges, Max's version");
  Hist.fHaloSideWedge=new TH1F(name,title,45,-0.5,44.5);
  new_folder->Add(Hist.fHaloSideWedge);
  sprintf(name,"HaloSeedWedgeH");
  sprintf(title,"number of HAD towers in seed wedge, my version");
  Hist.fHaloSeedWedgeH=new TH1F(name,title,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloSeedWedgeH);
  sprintf(name,"HaloSideWedgeH");
  sprintf(title,"number of HAD towers in side wedges, my version");
  Hist.fHaloSideWedgeH=new TH1F(name,title,45,-0.5,44.5);
  new_folder->Add(Hist.fHaloSideWedgeH);
  sprintf(name,"HaloEastNHad");
  sprintf(title,"number of HAD towers on East side");
  Hist.fHaloEastNHad=new TH1F(name,title,35,-0.5,34.5);
  new_folder->Add(Hist.fHaloEastNHad);
  sprintf(name,"HaloWestNHad");
  sprintf(title,"number of HAD towers on West side");
  Hist.fHaloWestNHad=new TH1F(name,title,35,-0.5,34.5);
  new_folder->Add(Hist.fHaloWestNHad);
  sprintf(name,"HaloNHad");
  sprintf(title,"number of HAD towers on East+West");
  Hist.fHaloNHad=new TH1F(name,title,70,-0.5,69.5);
  new_folder->Add(Hist.fHaloNHad);
  sprintf(name,"HaloEmSeedE");
  sprintf(title,"EM energy of seed wedge towers");
  Hist.fHaloEmSeedE=new TH1F(name,title,1000,0.0,1000.0);
  new_folder->Add(Hist.fHaloEmSeedE);
  sprintf(name,"HaloEmSideE");
  sprintf(title,"EM energy of side wedge towers");
  Hist.fHaloEmSideE=new TH1F(name,title,1000,0.0,1000.0);
  new_folder->Add(Hist.fHaloEmSideE);
  sprintf(name,"HaloHadSeedE");
  sprintf(title,"Plug HAD energy of seed wedge towers");
  Hist.fHaloHadSeedE=new TH1F(name,title,1000,0.0,1000.0);
  new_folder->Add(Hist.fHaloHadSeedE);
  sprintf(name,"HaloSeedWedgeHadE");
  sprintf(title,"HAD energy of seed wedge towers");
  Hist.fHaloSeedWedgeHadE=new TH1F(name,title,5000,0.0,100.0);
  new_folder->Add(Hist.fHaloSeedWedgeHadE);
  sprintf(name,"HaloEast");
  sprintf(title,"HaloEast");
  Hist.fHaloEast=new TH1F(name,title,10,-0.5,9.5);
  new_folder->Add(Hist.fHaloEast);
  sprintf(name,"HaloWest");
  sprintf(title,"HaloWest");
  Hist.fHaloWest=new TH1F(name,title,10,-0.5,9.5);
  new_folder->Add(Hist.fHaloWest);
  sprintf(name,"HaloEastWest");
  sprintf(title,"HaloWest+HaloEast");
  Hist.fHaloEastWest=new TH1F(name,title,10,-0.5,9.5);
  new_folder->Add(Hist.fHaloEastWest);
  sprintf(name,"HaloSeedWedgeR");
  sprintf(title,"number of towers in seed wedge, Ray's version");
  Hist.fHaloSeedWedgeR=new TH1F(name,title,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloSeedWedgeR);
  sprintf(name,"HaloSideWedgeR");
  sprintf(title,"number of towers in side wedges, Ray's version");
  Hist.fHaloSideWedgeR=new TH1F(name,title,45,-0.5,44.5);
  new_folder->Add(Hist.fHaloSideWedgeR);
  sprintf(name,"HaloTwrInRawR");
  sprintf(title,"number of continues towers in seed wedge, Ray's version");
  Hist.fHaloTwrInRawR=new TH1F(name,title,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloTwrInRawR);

  sprintf(name,"HaloCesStripEoverE");
  sprintf(title,"strip energy E_{CES}/E_{#gamma}");
  Hist.fHaloCesStripEoverE=new TH1F(name,title,500,0.0,5.0);
  new_folder->Add(Hist.fHaloCesStripEoverE);
  sprintf(name,"HaloCesWireEoverE");
  sprintf(title,"wire energy E_{CES}/E_{#gamma}");
  Hist.fHaloCesWireEoverE=new TH1F(name,title,500,0.0,5.0);
  new_folder->Add(Hist.fHaloCesWireEoverE);
  sprintf(name,"HaloHadTDC");
  sprintf(title,"photons Had TDC");
  Hist.fHaloHadTDC=new TH1F(name,title,800,-100.0,100.0);
  new_folder->Add(Hist.fHaloHadTDC);

  sprintf(name,"HaloPhiWedge");
  sprintf(title,"wedge number of halo candidate");
  Hist.fHaloPhiWedge=new TH1F(name,title,24,-0.5,23.5);
  new_folder->Add(Hist.fHaloPhiWedge);
  sprintf(name,"HaloEta");
  sprintf(title,"#eta of halo candidate");
  Hist.fHaloEta=new TH1F(name,title,120,-3.0,3.0);
  new_folder->Add(Hist.fHaloEta);
  sprintf(name,"HaloIso");
  sprintf(title," CalIso of halo candidate");
  Hist.fHaloIso= new TH1F(name,title,140,-4.0,10.0);
  new_folder->Add(Hist.fHaloIso);
  sprintf(name,"HaloHadEm");
  sprintf(title,"Had/Em of halo candidate");
  Hist.fHaloHadEm= new TH1F(name,title,400,-0.5,1.5);
  new_folder->Add(Hist.fHaloHadEm);
  sprintf(name,"HaloCesChi2");
  sprintf(title," %s%s",Folder,": CES(strip+wire) #chi^{2} of halo candidate");
  Hist.fHaloCesChi2= new TH1F(name,title,500,0.0,100.0);
  new_folder->Add(Hist.fHaloCesChi2);
  sprintf(name,"HaloEt");
  sprintf(title," corrected E_{T} of halo candidate");
  Hist.fHaloEt= new TH1F(name,title,1000,0.0,1000.0);
  new_folder->Add(Hist.fHaloEt);


	// my stuff
  sprintf(name,"Halo_Et_to_MetRatio");
  sprintf(title,"halo candidate E_{T}^{corr} to ME_{T} ratio");
  Hist.fHaloEt_to_Met = new TH1F(name,title,500,0,50);
  new_folder->Add(Hist.fHaloEt_to_Met);

  sprintf(name,"Halo_Met_delPhi");
  sprintf(title,"halo candidate and ME_{T} #Delta#phi");
  Hist.fHalo_Met_delPhi = new TH1F(name,title,140,0,3.5);
  new_folder->Add(Hist.fHalo_Met_delPhi);

  sprintf(name,"Halo_Met_InvMass");
  sprintf(title,"halo candidate and ME_{T} invriant mass");
  Hist.fHalo_Met_InvMass = new TH1F(name,title,100,0,100);
  new_folder->Add(Hist.fHalo_Met_InvMass);

  sprintf(name,"NHalo_to_NPho");
  sprintf(title,"# of halo candidates to reconstructed #gamma");
  Hist.fnHalo_to_nPho = new TH1F(name,title,20,0,2);
  new_folder->Add(Hist.fnHalo_to_nPho);

  sprintf(name,"HaloEvt_nPho_to_nJets");
  sprintf(title,"# of reconstructed #gamma to Jets ratio");
  Hist.fHalo_nPho_to_nJets = new TH1F(name,title,20,0,2);
  new_folder->Add(Hist.fHalo_nPho_to_nJets);
  
  sprintf(name,"HaloEvt_SeedWedge_Nhad");
  sprintf(title,"Em towers in seed wedge and N Had towers (eat+west)");
  Hist.fSeedWedge_HadTowers = new TH2F(name,title,20,0,20,20,0,20);
  Hist.fSeedWedge_HadTowers->GetXaxis()->SetTitle("N Had");
  Hist.fSeedWedge_HadTowers->GetYaxis()->SetTitle("Seed Em Towers");
  new_folder->Add(Hist.fSeedWedge_HadTowers);
	
  return;
}












/*-------------------------------------------------------------------*/
//_____________________________ Filling beam halo histograms 
void  TInit::FillBeamHaloStudyHistograms(BeamHaloStudyHisto_t& Hist, 
							 HaloStuff *beamhalo, TStnPhoton* Pho)
{

  Hist.fHaloSeedWedgeEM_SeedWedgeHAD->Fill(beamhalo->seedwedgeH,beamhalo->seedwedge);
    //______________________ Max's params
  Hist.fHaloSeedWedgeH->Fill(beamhalo->seedwedgeH);   
  Hist.fHaloSideWedgeH->Fill(beamhalo->sidewedgeH); 
  Hist.fHaloSeedWedge->Fill(beamhalo->seedwedge);   
  Hist.fHaloSideWedge->Fill(beamhalo->sidewedge); 
  Hist.fHaloEastNHad->Fill(beamhalo->eastNhad);
  Hist.fHaloWestNHad->Fill(beamhalo->westNhad);
  Hist.fHaloEmSeedE->Fill(beamhalo->emseedE);  
  Hist.fHaloEmSideE->Fill(beamhalo->emsideE);   
  Hist.fHaloHadSeedE->Fill(beamhalo->hadseedE); 
  Hist.fHaloSeedWedgeHadE->Fill(beamhalo->seedwedgeHadE);
  Hist.fHaloNHad->Fill(beamhalo->eastNhad+beamhalo->westNhad);  
    //______________________ Ray's params
  Hist.fHaloEast->Fill(Pho->HaloEast()); 
  Hist.fHaloWest->Fill(Pho->HaloWest());
  Hist.fHaloEastWest->Fill(Pho->HaloWest()+Pho->HaloEast());
  Hist.fHaloCesStripEoverE->Fill(Pho->CesStripE()/Pho->Momentum()->E());
  Hist.fHaloCesWireEoverE->Fill(Pho->CesWireE()/Pho->Momentum()->E()); 
  Hist.fHaloHadTDC->Fill(Pho->Time());     
  Hist.fHaloSideWedgeR->Fill(Pho->Seedwedge()); 
  Hist.fHaloSeedWedgeR->Fill(Pho->Sidewedges()); 
  Hist.fHaloTwrInRawR->Fill(Pho->Ncontig());  

  Hist.fHaloPhiWedge->Fill(Pho->PhiSeedIndex()); 
  Hist.fHaloEta->Fill(Pho->Eta());    
  Hist.fHaloIso->Fill(Pho->EIso4(2));    
  Hist.fHaloHadEm->Fill(Pho->HadEm());
  Hist.fHaloCesChi2->Fill(Pho->Chi2()); 
  Hist.fHaloEt->Fill(Pho->Etc());
 
 	// my stuff
	Hist.fHaloEt_to_Met->Fill(Pho->Etc() / fMetBlock->Met(0));
 	float phi_1 = Pho->Phi();
	float phi_2 = fMetBlock->MetPhi(0);
 
	float DelPhi = phi_1 - phi_2;

	if (DelPhi >kPI) {
		DelPhi -= kTWOPI;
	}
	if (DelPhi < -kPI) {
		DelPhi += kTWOPI;
	}

	Hist.fHalo_Met_delPhi->Fill(fabs(DelPhi));
	Hist.fSeedWedge_HadTowers->Fill(beamhalo->eastNhad+beamhalo->westNhad,
			beamhalo->seedwedge);   
  return;
}





/*=====================================end== BH STUFF =============*/


/*===================================BEGIN CR STUFF =============*/
void TInit::CosmicRayStudy()
{
	for (int i=0; i< superpho.size(); i++) { // loop only over candidates which passed ID cuts
		double time = GetEmTiming(superpho[i].index);
		fTimingHisto.fEmtTimePho->Fill(time);
		if (MyCosmicsCut(superpho[i].pho) == 1 || (time > 30 && time < 90) ) { // cosmic photons
			counter.cosmicCandidates++;
			DoMyBeamHalo(superpho[i].index);
			FillBeamHaloStudyHistograms(fCosmicHist,&myHaloStuff,superpho[i].pho);
			
		}
	} 
}

/*-------------------------------------------------------------------*/
// runN>=190851  first run for EM timing system
/*-------------------------------------------------------------------*/
double TInit::GetEmTiming(const int phoInd)
{
		
		
	double time = -99999.99;
	if (phoInd < 0) {
		std::cout << "ERR: invalid photon index!." <<std::endl;
		return time;
	}
	
	int runN = GetHeaderBlock()->RunNumber();

	if (runN >= 190851 && GetHeaderBlock()->McFlag() == false) { // first run for EM timing system
		TEmTimingModule* myEmTiming = (TEmTimingModule*) ((TStnAna*) GetAna()->GetModule("EmTimingAna"));
		
      if (myEmTiming == NULL) {
			std::cout << "ERR: EmTiminigModule not found." <<std::endl;
		} else {		
      	//__________________ timing for photons
			int Nhits=-1;
			int NGoodHits=-1;
			EmTimingTower* _emTimeTower = NULL;
			CalDataArray calholder;
			
			MatchCalorTowers(phoInd, &calholder);
				
	      for (int icount = 0 ; icount < calholder.size(); icount++) {
				TCalTower* ctower = calholder[icount];
				int iEta = ctower->IEta();
				int iPhi = ctower->IPhi();
		  
				_emTimeTower = myEmTiming->TEmTimingModule::getTimingTower(0,iEta,iPhi); // 0=CEM; 1=PEM; 2=CHA; 3=WHA; 4=PHA 
		
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
int TInit::MyCosmicsCut(TStnPhoton* Pho)
{
	int stat = 0;
	int runN = GetHeaderBlock()->RunNumber();
	if (runN<190851) { // first run for EM timing system
	 
		if((Pho->NCosStubPho()) > 0) stat=1;
	}

return stat;
}



/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TInit::BookTimingHistograms(EmTimingHisto_t& Hist,
															 char* Folder) {
	
	TFolder* new_folder = GetHistoFolder(Folder,Folder);
  char name [200];
  char title[200];

  sprintf(name,"EmtTimePho");
  sprintf(title,"Times of all reconstructed photons");
  Hist.fEmtTimePho=new TH1F(name,title,480,-80.0,160.0);
  new_folder->Add(Hist.fEmtTimePho);

  return;
}




/*=================================== END  CR STUFF =============*/
/*=================================== PHOENIX REJECTION STUFF =============*/
void TInit::PhoenixRejectionStudy()
{
	TStnEvent *evt = GetEvent();
	
	for (int i=0; i < superpho.size(); i++) {
		if (superpho[i].Detector != 0) continue;
		if (MyPhoFakeIDcut(superpho[i].pho,evt) == 1) {
			counter.phoenixPhotons++;
			FillPhotonPlots(superpho[i], fPhoenixPhoHist);
		}
		
		if ( (superpho[i].TightID == 0) &&  (MyPhoFakeIDcut(superpho[i].pho,evt) == 1) ) { // found a phoenix photon
			//make phoenix plots
			counter.phoenixPhoPassPID++;
			FillPhotonPlots(superpho[i], fPhoenixPhoPassPIDHist);
		} if (superpho[i].TightID == 0 ) {
			// make non-phoenix plots
			FillPhotonPlots(superpho[i], fPhotonHist);
		} 
	}

}
/*-------------------------------------------------------------------*/
//__ returns 1 if photon has a matched phoenix track
/*-------------------------------------------------------------------*/
int TInit::MyPhoFakeIDcut(TStnPhoton* Pho, TStnEvent* event)
{
	int passcode=0;
	TStnElectron* phele = TPhotonUtil::MatchPhoenixElectron(event,Pho);
	TStnTrack*    phTrk = TPhotonUtil::PhoenixTrack(event,phele);
	if(phTrk!=NULL) passcode=1; 
	return passcode;
}

/*=================================== END  PHOENIX REJECTION =============*/
/*-------------------------------------------------------------------*/
//tag conversion electrons faking photons
/*-------------------------------------------------------------------*/
void TInit::TagConversions()
{
	TStnEvent* event = GetEvent();
	for (int i =0; i < superpho.size(); i++) {
		TStnElectron* Ele = TPhotonUtil::MatchElectron(event,superpho[i].pho); // get matching electroni
		double minsep =0, mindt=0, radius=0;
		if (IsConversionElectron(Ele, minsep,mindt,radius)) {
			superpho[i].Conversion = true;
		}
	}
}

/*-------------------------------------------------------------------*/
// Isconversion electrons
/*-------------------------------------------------------------------*/
bool TInit::IsConversionElectron(TStnElectron* ele, 
				     double& minsep, double& mindt, double& radius) 
/*{{{*/
{

	const double dlamCut = 0.04;
	const double sepCut  = 0.2;
	const double ptCut   = 0.0;
	const bool reqCot    = false;

	int nconv = 0, ntrident = 0;
	bool hasCot = false;

	minsep =  999;
	mindt  =  999;
	radius = -999;

	int ntrk = fTrackBlock->NTracks();

	// this should not be happening
	int iEleTrk = ele->TrackNumber();
	if(iEleTrk<0 || iEleTrk>=ntrk) return false;
	TStnTrack* eleTrk = fTrackBlock->Track(iEleTrk);

	double  sep=-999, dt=-999;
	float conv_sign=-999;
	bool cand = false;

	for(int j=0; j<ntrk; j++) {
		TStnTrack* oTrk = fTrackBlock->Track(j);
		
		// only ask for a track from the electron track
		if ( j==iEleTrk ) continue;

		// require number of COT axial and stereo segments > 1
		hasCot = ((oTrk->NAxSeg()> 1) && (oTrk->NStSeg() > 1));
		if (!hasCot && reqCot)continue;
		 
		dt        = eleTrk->Lam0() - oTrk->Lam0();
		conv_sign = eleTrk->Charge() + oTrk->Charge();

    	bool eleCand = false;
    	
		for (int i=0; i < fElectronBlock->NElectrons(); i++) {
      	if (j == fElectronBlock->Electron(i)->TrackNumber()) eleCand = true;
    	}

  
		// opposite charge
		if (conv_sign != 0)continue;

		// delta cottheta cut
		if (fabs(dt) > dlamCut)continue;

		double rconv;
		sep = ConversionSep(eleTrk,oTrk,rconv);

		//separation cuts
		if (fabs(sep)> sepCut) continue;
    	
		nconv++;
    
		if (fabs(sep) <= minsep && fabs(dt) <= mindt) {
				minsep = fabs(sep);
				mindt  = fabs(dt);
				radius = rconv;
		}

		 // now check for trident
		 // by looking for a partner of second conversion 
		for (int k=0; k<ntrk; k++) {
			if ( (k != iEleTrk) && (k != j) ) {
				TStnTrack* kTrk = fTrackBlock->Track(k);
				dt        = oTrk->Lam0() - kTrk->Lam0();
				conv_sign = oTrk->Charge() + kTrk->Charge();
				sep       = ConversionSep(oTrk,kTrk,rconv);
				hasCot 	 = ( (kTrk->NAxSeg() > 1) && (kTrk->NStSeg() > 1) );

				if ( (kTrk->Pt() > ptCut) && (hasCot || (!reqCot)) &&
			   	conv_sign == 0 && fabs(dt) < dlamCut && fabs(sep) < sepCut) {
			  		ntrident++;
				}
			}
		} // end trident loop
	}

  
  if (ntrident > 0) cand = false;
  else if (nconv > 0) cand = true;
  else cand = false;
  return cand;
}
/*}}}*/

/*-------------------------------------------------------------------*/
// calculate seperation between the twoo traks for the conversion electrons
/*-------------------------------------------------------------------*/
double TInit::ConversionSep(TStnTrack* t0, TStnTrack* t1, double & rconv)
/*{{{*/
{
  double  r1, r2, x1, x2, y1, y2, ff, dd, sep;

  r1 = 1.0/fabs(2*t0->C0());
  ff = t0->Charge()*r1+t0->D0();
  x1 = -ff*sin(t0->Phi0());
  y1 =  ff*cos(t0->Phi0());

  r2 = 1.0/fabs(2*t1->C0());
  ff = t1->Charge()*r2+t1->D0();
  x2 = -ff*sin(t1->Phi0());
  y2 =  ff*cos(t1->Phi0());

  dd  = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));

  if(t0->Charge()!=t1->Charge()) {
    sep = dd-r1-r2;
  } else {
    if (r1<r2) sep = r2 - dd - r1;
    else       sep = r1 - dd - r2;
  }

  double dx = x2 - x1;
  double dy = y2 - y1;
  double dr = sqrt( dx*dx + dy*dy );
    
  double xconv = x1 + dx/dr * (r1 + sep/2);
  double yconv = y1 + dy/dr * (r1 + sep/2);
  rconv = sqrt( xconv*xconv + yconv*yconv );	

  return sep;
}
/*}}}*/


/*-------------------------------------------------------------------*/
// check if the photon is fidual before matching to HEPG Level for W estimate.
/*-------------------------------------------------------------------*/
int TInit::FiducialLepton(TLorentzVector vec_old, double z_new)
/*{{{*/
{
  int fid=0;
  double cmx_Z=240.0;
  double cmx_R=183.9;

  double hemi = (vec_old.Eta()>0.0? 1.0 : -1.0);
  double sint = sin(vec_old.Theta());
  double cott = hemi*sqrt(1.0-sint*sint)/sint;

  double z = cmx_R*cott + z_new;
  if(fabs(z)<cmx_Z) fid=1;

  return fid;
} //FiducialLeopton
/*}}}*/

/*-------------------------------------------------------------------*/
// generate stuff needed by sasha's code
/*-------------------------------------------------------------------*/
void TInit::InitForSashaCode(int ientry)
/*{{{*/
{
 	//std::cout << " IN initfor sasha" << std::endl;
	int Npho = fPhotonBlock->NPhotons();
	//std::cout << "Npho=" << Npho << std::endl;

	for (int i =0 ; i < Npho; i++) {
		TStnPhoton* pho = fPhotonBlock->Photon(i);

		if	(TightPhotonIDcut(pho) ==0) {
			Nphotons++;  //done

			TLorentzVector corvec;
			corvec.SetPxPyPzE( pho->Momentum()->Px(),pho->Momentum()->Py(),
					   pho->Momentum()->Pz(),pho->Momentum()->E());
			phoCorrected.push_back(corvec);

			TLorentzVector rawvec;
			rawvec.SetPtEtaPhiM(pho->Et(),pho->Eta(),pho->Phi(),0.0);
			phoUncorrected.push_back(rawvec);

			phoEtaDet.push_back(pho->DetEta());
			phoHadEm.push_back(pho->HadEm());
			phoCesX.push_back(pho->XCes());
			phoCesZ.push_back(pho->ZCes());
			phoIndex.push_back(i);


		} // if
	} // for


/*	if (Nphotons > 0) {
	//std::cout << "Nphotons=" << Nphotons <<std::endl;
	for (int i =0 ; i < Nphotons ; i++) {
		std::cout << i << "# INDEX=" << GetmyPhoInd(i) << std::endl;
		std::cout << i << "# DetEta=" << GetphoEtaDet(i) << "\tHadEm=" << GetphoHadEm3(i) << std::endl;
		std::cout << i << "# CesX  =" << GetphoCesX(i) <<   "\tCesZ =" << GetphoCesZ(i) << std::endl;
		TLorentzVector corvec = *(GetmyCorrPho(i));
		TLorentzVector rawvec = *(GetmyUncorrPho(i));
		std::cout << "Corr.E,PX,PY,PZ=" << corvec.E()<<"," << corvec.Px()<< "," << corvec.Py() <<","<< corvec.Pz() << std::endl;
		std::cout << "RAW.E,PX,PY,PZ=" << rawvec.E() <<"," << rawvec.Px() <<","<< rawvec.Py() <<","<< rawvec.Pz() << std::endl;
	}
	}
*/
	//std::cout << "out of sasha init " << std::endl;

} //InitForSashaCode
/*}}}*/

/*-------------------------------------------------------------------*/
//	Sorts the Photon Block according to EtCorr, from highest to lowest
// & fill the photon info to my SuperPhoton list
/*-------------------------------------------------------------------*/
void TInit::SetPhotonList()
/*{{{*/
{
	int Npho = fPhotonBlock->NPhotons();
	//std::cout << " TGI: Npho=" << Npho << std::endl;
	std::vector<TStnPhoton*> photon;
	std::vector<int> index;

	for (int i = 0; i < Npho ; i++) {
		photon.push_back(fPhotonBlock->Photon(i));
		index.push_back(i);
	}

	if (Npho>1) {   //________________________ dont really need this since i check for this at the beginnig!
		int loopupto = Npho-1;
		while (loopupto !=0) {
			for (int i=0; i<loopupto;i++) {
				if ( (photon[i]->ECorr()*photon[i]->SinTheta()) <
						(photon[i+1]->ECorr()*photon[i+1]->SinTheta()) ) {

					TStnPhoton *temp = photon[i];
					photon[i] = photon[i+1];
					photon[i+1] = temp;

					int t = index[i];
					index[i] = index[i+1];
					index[i+1] = t;

				}
			} // for
			loopupto--;
		} // while
	} // if

	double zvx  = 0;
	if (fVertexBlock->GetBestVertex(12,1) != NULL) {
		zvx = fVertexBlock->GetBestVertex(12,1)->Z();
	}

	TStnEvent* event = GetEvent();
	TStntuple::CorrectPhotonEnergy(event,zvx);
	//std::cout << "pho size=" << photon.size() << std::endl;

	TTightPhotonID tid;

	for (int i=0; i < photon.size(); i++) {
		TSuperPhoton sp;
		sp.index   		= index[i];
		sp.pho 			= photon[i];
		sp.Detector 	= photon[i]->Detector();
		sp.DetEta 		= photon[i]->DetEta();
		sp.Phi 			= photon[i]->Phi();
		sp.Etraw 		= photon[i]->Et();
		sp.ECorr 		= photon[i]->ECorr();
		sp.EtCorr 		= photon[i]->ECorr() * photon[i]->SinTheta();
		sp.XCes 			= photon[i]->XCes();
		sp.ZCes 			= photon[i]->ZCes();
		sp.HadEm 		= photon[i]->HadEm();
		sp.Chi2Mean 	= photon[i]->Chi2Mean();
		sp.N3d 			= photon[i]->N3d();
		sp.IsoEtCorr	= photon[i]->EIso4(2);
		sp.TrkPt 		= photon[i]->Pt();
		sp.TrkIso 		= photon[i]->SumPt4();
		sp.CesWireE2	= photon[i]->CesWireE2();
		sp.CesStripE2	= photon[i]->CesStripE2();
		sp.LooseID		= LoosePhotonIDcut(photon[i]);
		sp.TightID		= TightPhotonIDcut(photon[i]);
		sp.PhotonLikeID= CEMPhoEleIDcut(photon[i], event, zvx);  // i am using the negative logic, 0=pass 1= fail
		sp.BeamHalo 	= false;				// for now???? call "TagBeamHalo" TODO
		sp.Cosmic   	= false;				// for now???	call "TagCosmics" TODO
		sp.Conversion 	= false;				// call "TagConversion" seperately to set this
		sp.HadTDCtime 	= photon[i]->Time();
		sp.corvec 		= *(photon[i]->Momentum());
		sp.rawvec.SetPtEtaPhiM(photon[i]->Et(),photon[i]->Eta(),photon[i]->Phi(),0.0);

/*		
		tm = new TMathModule();
		int id = tid.GetIdWord(photon[i]);
		int tmp = (sp.TightID >> 8) & 0x7ff;
		int tmp1 = (sp.TightID >> 19) & 0x1;
		int oid = (tmp << 1) | tmp1;
		std::cout << " sp.TightID =" << sp.TightID << std::endl;
		tm->binary(sp.TightID);
		std::cout << std::endl;
		std::cout << " oid        =" << oid << std::endl;
		tm->binary(oid);
		std::cout << std::endl;
		std::cout << " id         =" << id << std::endl;
		tm->binary(id);
		std::cout << std::endl;
	//	assert(id == oid);
		std::cout << "pushing back .." << std::endl;
*/		
		superpho.push_back(sp);
		//std::cout << "DONE." << std::endl;
		//delete tm;
	}

} //SetPhotonList
/*}}}*/


/*****************  PHOTON LIKE ID CUT FUNCTION ************************/
//_____ function returns 0 if electrons passes cuts
//      requires zvx of highest Pt class12 vertex
unsigned int TInit::CEMPhoEleIDcut(TStnPhoton* Pho, TStnEvent* event, double zvx)
/*{{{*/
{
	//______________________________________ only CEM photons with Et>7 GeV
	if (Pho->Detector()!=0) {
		return 1;
	}
	if (Pho->Etc()< 30.) {
		return 2;
	}

	//______________________________________ HADEM cut using 3 towers
	float cutMin_HADEM=0.0;
	float cutMax_HADEM=0.055+0.00045*(Pho->Momentum()->E());
	if ((Pho->HadEm())<cutMin_HADEM || (Pho->HadEm())>cutMax_HADEM){
		return 3;
	}

	//______________________________________ CES Chi^2 cut
	if ((Pho->Chi2())<0.0 || (Pho->Chi2())>20.0) {
		return 4;
	}

	//______________________________________ N3D cut
	if ((Pho->N3d())==0) {
		return 5; // it is not an electron if N3d=0
	}

	//______________________________________ second N3D cut
	if ((Pho->N3d())>2) {
		return 6;  // no more than 2 tracks (1st--ele,2nd--to match cuts for pho)
	}

	//______________________________________ cut on 2nd max Pt track in cluster if N3D=2
	if ((Pho->N3d())==2) {
		float trkPtcut_max=1.0+0.005*(Pho->Etc());

		if ((Pho->Pt2())>trkPtcut_max) {
			return 7;
		}
	}

	TStnElectron* Ele=TPhotonUtil::MatchElectron(event,Pho); // get matching electron

	//______________________________________ E/p cut for electron
	if (Ele==NULL) {
		return 8;
	}

	if (Ele->TrackNumber()<0) {
		return 9; // there should be a track
	}
	if (fabs(Ele->Z0()-zvx)>3.0) {
		return 10; // only events with ele from best vertex
	}

	int ele_trk = Ele->TrackNumber();
	TStnTrack* Trk = fTrackBlock->Track(ele_trk);

	if ((Trk->Algorithm())==2) {
		if ((Ele->TrackBcPt())<50.0 &&
	 			(Ele->EOverP()<0.8 || Ele->EOverP()>1.2)) return 11;
	} else {
		if ((Ele->TrackPt())<50.0 &&
	 			(Ele->EOverP()<0.8 || Ele->EOverP()>1.2)) return 11;
	}

	//______________________________________ CalIso4 cut
	float cutMax_CalIso=0.0;

	if ((Pho->Etc()) < 20.0) cutMax_CalIso = 0.1 * (Pho->Etc());
	else cutMax_CalIso = 2.0 + 0.02 * ( Pho->Etc() - 20.0 );

	if ((Pho->EIso4(2)) > cutMax_CalIso) return 12;

	//______________________________________ TrkIso4 cut
	float trkIso=Pho->SumPt4() - Pho->Pt();

	if (trkIso < 0.0) trkIso = Pho->SumPt4();

	if (trkIso < (trkIso > (2.0 + 0.005 * (Pho->Etc()) ) ) ) return 13;

	//______________________________________ Energy of 2nd CES cluster (Wire & Strip)
	float _Et2ndCES=0.0;
	if ((Pho->CesStripE2()) > (Pho->CesWireE2())) _Et2ndCES = (Pho->CesStripE2()) * (Pho->SinTheta());
	else _Et2ndCES = (Pho->CesWireE2()) * (Pho->SinTheta());

	float cutMax_2ndCes = 0.0;
	if ((Pho->Etc()) < 18.0) cutMax_2ndCes = 0.14 * (Pho->Etc());
	else cutMax_2ndCes = 2.4 + 0.01 * (Pho->Etc());

	if (_Et2ndCES > cutMax_2ndCes) return 14;

	//_______________________________________ Fiducial cuts
	if (fabs(Pho->XCes()) > 21.0) return 15;
	if (fabs(Pho->ZCes()) < 9.0 || fabs(Pho->ZCes()) > 230.0) return 16;

	return 0;

} // CEMPhoEleIDcut

/*}}}*/


/*-------------------------------------------------------------------*/
// Tight Photon ID CUTS
/*-------------------------------------------------------------------*/
unsigned int TInit::TightPhotonIDcut(TStnPhoton *pho)
/*{{{*/
{
	unsigned int IDWord = 0x0;
	unsigned int IDWord1 = 0x0;

	int   detector 	= pho->Detector();
	float ECorr 		= pho->ECorr();
	float EtCorr 		= ECorr * pho->SinTheta(); //total corrected energy // this is waht i want to cut on
	float XCes			= pho->XCes();
	float ZCes 			= pho->ZCes();
	float HadEm 		= pho->HadEm();
	float Chi2Mean 	= pho->Chi2Mean();  			// (strip+wire)/2
	int   N3d			= pho->N3d();
	float IsoEtCorr	= pho->EIso4(2);				//Corrected for leakage and MI
	float TrkPt  		= pho->Pt();     				//using max pt track in the cluster
	float TrkIso 		= pho->SumPt4();
	float CesWireE2	= pho->CesWireE2();
	float CesStripE2	= pho->CesStripE2();

	if (detector != 0						) IDWord |= kCentral_T;
	if (EtCorr < 30. 						) IDWord |= kEtCorr7_T; 		// has no effect! looser cut
	if (fabs(XCes) > 21					) IDWord |= kXCes_T;
	if (fabs(ZCes) < 9 ||  fabs(ZCes) > 230) IDWord |= kZCes_T;  	// same as loose cuts
	if ( !((HadEm < 0.125) || (HadEm < (0.055 + 0.00045 * ECorr)))) IDWord |= kHadEm_T;
	if (Chi2Mean > 20						) IDWord |= kChi2Mean_T;
	if (N3d > 1								) IDWord |= kN3d_T;
	if (TrkPt > (1+0.005*EtCorr)		) IDWord |= kTrkPt_T;
	if (TrkIso > (2.0+0.005*EtCorr)	) IDWord |= kTrkIso_T;

	if (EtCorr < 20) {
		if (IsoEtCorr > 0.1*EtCorr) IDWord |= kIsoEtCorr_T;
	} else {
		if (IsoEtCorr > (2.0+0.02*(EtCorr-20.0)) ) IDWord |= kIsoEtCorr_T;
	}

	if (EtCorr<18) {
		if ( (CesWireE2 * pho->SinTheta()) > (0.14*EtCorr)) IDWord |= kCesWireE2_T;
		if ( (CesStripE2 * pho->SinTheta()) > (0.14*EtCorr)) IDWord |= kCesStripE2_T;
	} else {
		if ( (CesWireE2 * pho->SinTheta()) > (2.4+0.01*EtCorr)) IDWord |= kCesWireE2_T;
		if ( (CesStripE2 * pho->SinTheta()) > (2.4+0.01*EtCorr)) IDWord |= kCesStripE2_T;
	}


	if (detector != kT_CENTRAL					) IDWord1 |= kCentral_T;
	if (EtCorr < kT_ETCORR 						) IDWord1 |= kEtCorr7_T; // has no effect! looser cut
	if (fabs(XCes) > kT_XCES					) IDWord1 |= kXCes_T;
	if (fabs(ZCes) < kT_ZCES_MIN ||  fabs(ZCes) > kT_ZCES_MAX	) IDWord1 |= kZCes_T;  // same as loose cuts
	if ( !((HadEm < kT_HADEM_CONST1) || (HadEm < (kT_HADEM_CONST2 + kT_HADEM_MF * ECorr))) 	) IDWord1 |= kHadEm_T;
	if (Chi2Mean > kT_CHI2MEAN					) IDWord1 |= kChi2Mean_T;
	if (N3d > kT_N3DTRACKS						) IDWord1 |= kN3d_T;
	if (TrkPt > ( kT_TRKPT_CONST+ kT_TRKPT_MF * EtCorr)		) IDWord1 |= kTrkPt_T;
	if (TrkIso > (kT_TRKISO_CONST + kT_TRKISO_MF * EtCorr)		) IDWord1 |= kTrkIso_T;

	if (EtCorr < kETCORR20) {
		if (IsoEtCorr > kT_ISOETCORR_MF_LTETC20 * EtCorr) IDWord1 |= kIsoEtCorr_T;
	} else {
	  if (IsoEtCorr > (kT_ISOETCORR_CONST1_GTETC20 + kT_ISOETCORR_MF_GTETC20 * (EtCorr- kT_ISOETCORR_CONST2_GTETC20)) )
		  	IDWord1 |= kIsoEtCorr_T;
	}

	if (EtCorr< kETCORR18) {
		if ( (CesWireE2 * pho->SinTheta()) > (kT_CESWIREE2_MF_LTETC18 * EtCorr)) IDWord1 |= kCesWireE2_T;
		if ( (CesStripE2 * pho->SinTheta()) > (kT_CESWIREE2_MF_LTETC18 * EtCorr)) IDWord1 |= kCesStripE2_T;
	} else {
		if ( (CesWireE2 * pho->SinTheta()) > (kT_CESWIREE2_CONST_GTETC18 + kT_CESWIREE2_MF_GTETC18 * EtCorr)) IDWord1 |= kCesWireE2_T;
		if ( (CesStripE2 * pho->SinTheta()) > (kT_CESWIREE2_CONST_GTETC18 + kT_CESWIREE2_MF_GTETC18 * EtCorr)) IDWord1 |= kCesStripE2_T;
	}

	//assert(IDWord == IDWord1);
	return IDWord;

} // TightPhotonCuts
/*}}}*/

/*-------------------------------------------------------------------*/
//	Loose Photon ID CUTS
/*-------------------------------------------------------------------*/
unsigned int TInit::LoosePhotonIDcut(TStnPhoton *pho)
/*{{{*/
{
	unsigned int IDWord = 0x0;
	unsigned int IDWord1 = 0x0;

	int detector	= pho->Detector();
	float EtCorr 	= pho->ECorr() * pho->SinTheta();	// corrected energy this is what i want to cut on
	float XCes 		= pho->XCes();
	float ZCes 		= pho->ZCes();
	float HadEm 	= pho->HadEm();
	float IsoEtCorr= pho->EIso4(2);	//__________________ corrected for leakage and MI
	float TrkPt  	= pho->Pt();    	//__________________ using max pt track in the cluster
	float TrkIso 	= pho->SumPt4();

	if (detector != 0					) IDWord |= kCentral_L;
	if (EtCorr < 30					) IDWord |= kEtCorr30_L;
	if (fabs(XCes) > 21				) IDWord |= kXCes_L;
	if ((fabs(ZCes) < 9) || (fabs(ZCes) > 230)	) IDWord |= kZCes_L;
	if (HadEm > 0.125					) IDWord |= kHadEm_L;
	if (TrkPt > 0.25*EtCorr				) IDWord |= kTrkPt_L;
	if (TrkIso > 5.0					) IDWord |= kTrkIso_L;

	if (EtCorr < 20) {
		if (IsoEtCorr > 0.15*EtCorr) IDWord |= kIsoEtCorr_L;
	} else {
		if (IsoEtCorr > 3.0*EtCorr	) IDWord |= kIsoEtCorr_L;
	}

	if (detector != kL_CENTRAL	) IDWord1 |= kCentral_L;
	if (EtCorr < kL_ETCORR		) IDWord1 |= kEtCorr30_L;
	if (fabs(XCes) > kL_XCES	) IDWord1 |= kXCes_L;
	if ((fabs(ZCes) < kL_ZCES_MIN) || (fabs(ZCes) > kL_ZCES_MAX)) IDWord1 |= kZCes_L;
	if (HadEm > kL_HADEM		) IDWord1 |= kHadEm_L;

	if (EtCorr < kETCORR20) {
		if (IsoEtCorr > kL_ISOETCORR_MF_LTETC20 * EtCorr) IDWord1 |= kIsoEtCorr_L;
	} else {
		if (IsoEtCorr > kL_ISOETCORR_MF_GTETC20 * EtCorr) IDWord1 |= kIsoEtCorr_L;
	}

	if (TrkPt > kL_TRKPT_MF * EtCorr) IDWord1 |= kTrkPt_L;
	if (TrkIso > kL_TRKISO		) IDWord1 |= kTrkIso_L;


	assert(IDWord == IDWord1);

	return IDWord;

} //LoosePhotonCuts
/*}}}*/

/*-------------------------------------------------------------------*/
void TInit::FillDataBlocks(int ientry)
{
	fTriggerBlock->GetEntry(ientry);
	fPhotonBlock->GetEntry(ientry);
	fJetBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
	fVertexBlock->GetEntry(ientry);
	fTrackBlock->GetEntry(ientry);
	fCalBlock->GetEntry(ientry);
	fGenpBlock->GetEntry(ientry);
	fPhxEleBlock->GetEntry(ientry);
	fPhxSiTrackBlock->GetEntry(ientry);
}

/*-------------------------------------------------------------------*/
// fill photons plots
/*-------------------------------------------------------------------*/
void TInit::FillPhotonPlots(TSuperPhoton& sp, PhotonPlots_t& PhotonPlot)
{
	PhotonPlot.Detector->Fill(sp.Detector);
	PhotonPlot.EtCorr->Fill(sp.EtCorr);
	PhotonPlot.XCes->Fill(sp.XCes);
	PhotonPlot.ZCes->Fill(sp.ZCes);
	PhotonPlot.HadEm->Fill(sp.HadEm);
	PhotonPlot.IsoEtCorr->Fill(sp.IsoEtCorr);
	PhotonPlot.Chi2Mean->Fill(sp.Chi2Mean);
	PhotonPlot.N3d->Fill(sp.N3d);
	PhotonPlot.TrkPt->Fill(sp.TrkPt);
	PhotonPlot.TrkIso->Fill(sp.TrkIso);
	PhotonPlot.Ces2Wire->Fill(sp.CesWireE2);
	PhotonPlot.HadTDCtime->Fill(sp.HadTDCtime); // this qty in TStnPhoton is useless
	PhotonPlot.Ces2Strip->Fill(sp.CesStripE2);
}

/*===================================================================*/
//	Create a folder in"Hist" 
/*===================================================================*/
TFolder* TInit::GetHistoFolder(char *name, char* title)
{
	char folder_name[200];
	char folder_title[200];
	TFolder* hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
	sprintf(folder_name,name);
	sprintf(folder_title,title);
	TFolder* new_folder = (TFolder*) hist_folder->FindObject(folder_name);
	if (! new_folder) new_folder = hist_folder->AddFolder(folder_name,folder_title,NULL);
	return new_folder;
}
//_____________________________________________________________________________
void TInit::BookGeneralHistograms()
/*{{{*/
{
	std::cout << "Booking General Plots...";
	
	GeneralPlot.photon_loose_and_tight_cuts_cumm = new TH1F("loose_tight_cuts_cummulative","Cummulative Photon Cuts (Loose and Tight)",22,0,22);
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(1,"Trigger");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(2,"L:Central");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(3,"L:Etc>30");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(4,"L:XCes");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(5,"L:ZCes");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(6,"L:HadEm");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(7,"L:TrkPt");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(8,"L:TrkIso");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(9,"L:IsoEtCorr");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(10,"T:EtCorr");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(11,"T:XCes");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(12,"T:ZCes");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(13,"T:HadEm");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(14,"T:Chi2Mean");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(15,"T:N3d");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(16,"T:TrkPt");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(17,"T:TrkIso");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(18,"T:IsoEt4Corr");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(19,"T:2ndCesWire");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(20,"T:2ndCesStrip");
	//GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(21,"2 Jets>15GeV");

	GeneralPlot.Nvertices = new TH1F("Nvertices","Number of Class12 vertices per event",10,0,10);
	GeneralPlot.vertexZ = new TH1F("BestVertex_Z","Best Class12 vertex Z position",100,-100,100);

	TFolder* new_folder = GetHistoFolder("General_Plots","General Plots");
	new_folder->Add(GeneralPlot.photon_loose_and_tight_cuts_cumm);	
	new_folder->Add(GeneralPlot.pho_passIDs);
	new_folder->Add(GeneralPlot.Nvertices);
	new_folder->Add(GeneralPlot.vertexZ);
	std::cout << "DONE\n";

} //BookGeneralHistograms
/*}}}*/

//_____________________________________________________________________________
void TInit::BookPhotonHistograms(PhotonPlots_t& PhotonPlot, char* Folder)
/*{{{*/
{
	std::cout << "Booking Photon Plots...";

	PhotonPlot.Detector 	= new TH1F("Detector","Photon Detector Region",3,0,3);
	PhotonPlot.EtCorr 	= new TH1F("EtCorr","Photon Corrected Et",600,0,120);
	PhotonPlot.XCes 		= new TH1F("XCes","Photon XCes",640,-32,32);
	PhotonPlot.ZCes 		= new TH1F("ZCes","Photon ZCes",2000,-250,250);
	PhotonPlot.HadEm 		= new TH1F("HadEm","Photon HadEm",50,0,0.5);
	PhotonPlot.IsoEtCorr = new TH1F("IsoEtCorr","Photon IsoCorr",1500,-5,10);
	PhotonPlot.Chi2Mean 	= new TH1F("Chi2Mean","Photon Chi2Mean (Wire+Strip/2)",800,0,80);
	PhotonPlot.N3d 		= new TH1F("N3d","Photon N3d",10,0,10);
	PhotonPlot.TrkPt 		= new TH1F("Trkpt","Photon TrkPt",1000,0,100);
	PhotonPlot.TrkIso 	= new TH1F("TrkIso","Photon TrkIso",150,0,15);
	PhotonPlot.Ces2Wire 	= new TH1F("Ces2Wire","Photon CES(2nd) Wire",400,0,40);
	PhotonPlot.Ces2Strip = new TH1F("Ces2Strip","Photon CES(2nd) Strip",400,0,40);
	PhotonPlot.HadTDCtime = new TH1F("HadTDCtime","Photon Had TDC time",200,-10,10);

	//use GetBinWidth() in here!!!! u moran
	PhotonPlot.EtCorr->SetYTitle("Events/0.2");
	PhotonPlot.XCes->SetYTitle("Events/0.1");
	PhotonPlot.ZCes->SetYTitle("Events/0.25");
	PhotonPlot.HadEm->SetYTitle("Events/0.01");
	PhotonPlot.IsoEtCorr->SetYTitle("Events/0.01");
	PhotonPlot.Chi2Mean->SetYTitle("Events/0.1");
	PhotonPlot.N3d->SetYTitle("Events/1.0");
	PhotonPlot.TrkPt->SetYTitle("Events/0.1");
	PhotonPlot.TrkIso->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Wire->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Strip->SetYTitle("Events/0.1");
	PhotonPlot.HadTDCtime->SetYTitle("Events/0.1");

	TFolder* new_folder = GetHistoFolder(Folder,Folder);
	new_folder->Add(PhotonPlot.Detector);
	new_folder->Add(PhotonPlot.EtCorr);
	new_folder->Add(PhotonPlot.XCes);
	new_folder->Add(PhotonPlot.ZCes);
	new_folder->Add(PhotonPlot.HadEm);
	new_folder->Add(PhotonPlot.IsoEtCorr);
	new_folder->Add(PhotonPlot.Chi2Mean);
	new_folder->Add(PhotonPlot.N3d);
	new_folder->Add(PhotonPlot.TrkPt);
	new_folder->Add(PhotonPlot.TrkIso);
	new_folder->Add(PhotonPlot.Ces2Wire);
	new_folder->Add(PhotonPlot.Ces2Strip);
	new_folder->Add(PhotonPlot.HadTDCtime);

	std::cout << "DONE.\n";

} //BookPhotonHistograms
/*}}}*/

//_____________________________________________________________________________
void TInit::BookElectronHistograms(PhotonPlots_t& ElectronPlot, char* Folder)
/*{{{*/
{
	std::cout << "Booking Electron Plots...";

	ElectronPlot.Detector 	= new TH1F("Detector","Wenu MC:  Electron Detector Region",3,0,3);
	ElectronPlot.EtCorr 	= new TH1F("EtCorr","Wenu MC:  Electron Corrected Et",600,0,120);
	ElectronPlot.XCes 		= new TH1F("XCes","Wenu MC:  Electron XCes",640,-32,32);
	ElectronPlot.ZCes 		= new TH1F("ZCes","Wenu MC:  Electron ZCes",2000,-250,250);
	ElectronPlot.HadEm 		= new TH1F("HadEm","Wenu MC:  Electron HadEm",50,0,0.5);
	ElectronPlot.IsoEtCorr = new TH1F("IsoEtCorr","Wenu MC:  Electron IsoCorr",1500,-5,10);
	ElectronPlot.Chi2Mean 	= new TH1F("Chi2Mean","Wenu MC:  Electron Chi2Mean (Wire+Strip/2)",800,0,80);
	ElectronPlot.N3d 		= new TH1F("N3d","Wenu MC:  Electron N3d",10,0,10);
	ElectronPlot.TrkPt 		= new TH1F("Trkpt","Wenu MC:  Electron TrkPt",1000,0,100);
	ElectronPlot.TrkIso 	= new TH1F("TrkIso","Wenu MC:  Electron TrkIso",150,0,150);
	ElectronPlot.Ces2Wire 	= new TH1F("Ces2Wire","Wenu MC:  Electron CES(2nd) Wire",400,0,40);
	ElectronPlot.Ces2Strip = new TH1F("Ces2Strip","Wenu MC:  Electron CES(2nd) Strip",400,0,40);
	ElectronPlot.HadTDCtime = new TH1F("HadTDCtime","Wenu MC:  Electron Had TDC time",200,-10,10);

	ElectronPlot.EtCorr->SetYTitle("Events/0.2");
	ElectronPlot.XCes->SetYTitle("Events/0.1");
	ElectronPlot.ZCes->SetYTitle("Events/0.25");
	ElectronPlot.HadEm->SetYTitle("Events/0.01");
	ElectronPlot.IsoEtCorr->SetYTitle("Events/0.01");
	ElectronPlot.Chi2Mean->SetYTitle("Events/0.1");
	ElectronPlot.N3d->SetYTitle("Events/1.0");
	ElectronPlot.TrkPt->SetYTitle("Events/0.1");
	ElectronPlot.TrkIso->SetYTitle("Events/0.1");
	ElectronPlot.Ces2Wire->SetYTitle("Events/0.1");
	ElectronPlot.Ces2Strip->SetYTitle("Events/0.1");
	ElectronPlot.HadTDCtime->SetYTitle("Events/0.1");

	TFolder* new_folder = GetHistoFolder(Folder,Folder);
	new_folder->Add(ElectronPlot.Detector);
	new_folder->Add(ElectronPlot.EtCorr);
	new_folder->Add(ElectronPlot.XCes);
	new_folder->Add(ElectronPlot.ZCes);
	new_folder->Add(ElectronPlot.HadEm);
	new_folder->Add(ElectronPlot.IsoEtCorr);
	new_folder->Add(ElectronPlot.Chi2Mean);
	new_folder->Add(ElectronPlot.N3d);
	new_folder->Add(ElectronPlot.TrkPt);
	new_folder->Add(ElectronPlot.TrkIso);
	new_folder->Add(ElectronPlot.Ces2Wire);
	new_folder->Add(ElectronPlot.Ces2Strip);
	new_folder->Add(ElectronPlot.HadTDCtime);

	std::cout << "DONE.\n";

} //BookElectronHistograms
/*}}}*/


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TInit::EndJob() {

	printf("----- end job: ---- %s\n",GetName());
	if (qMc)	std::cout << "RUN IS ON MC SAMPLE" << std::endl;
	else	std::cout << "RUN IS ON DATA SAMPLE" << std::endl;
	std::cout << "EVENTS RUN OVER ------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "EVENTS PASS THIS MODULE ---- = " << counter.evtsPassModule << std::endl;
	std::cout << "Halo candidates found ------ = " << counter.haloCandidates << std::endl;
	std::cout << "Halo cand. pass tight ID CUTS= " << counter.haloPassPhoId << std::endl;
	std::cout << "Cosmic candidates -----------= " << counter.cosmicCandidates << std::endl;
	std::cout << "Cosmic cand. pass tight ID   = " << counter.cosmicPassPhoId << std::endl;
	std::cout << "Phoenix candidates --------- = " << counter.phoenixPhotons << std::endl;
	std::cout << "Phoenix cand. pass tight ID  = " << counter.phoenixPhoPassPID << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
