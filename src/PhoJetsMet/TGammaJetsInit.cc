/*
when I say "pho_", generally it means the leading photon (most of the time)

READ THIS EVERYDAY!
keep it simple and readable. forget about performance!
*/

#include <iostream>
#include <fstream>
#include <TMath.h>
#include <numeric>
#include <algorithm>
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnTriggerTable.hh"
#include "samantha/PhoJetsMet/TGammaJetsInit.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/ana/TPrintModule.hh"
#include "Stntuple/ana/TMathModule.hh"


double const kPI = TMath::Pi();
double const kTWOPI = 2.*kPI;
const static unsigned int kPHO_ID_BITS = 19;      //total number photon ID bits (loose + tight)
const static unsigned int kPHO_LOOSE_ID_BITS = 8; //number of Loose ID variables
const static double       kPHO_ET_CUT = 30; 	  //my leading photon. so i make the cut on 30GeV!


const static double 	kETCORR20= 20.0;	// Et SPLIT for ISO CUT
const static double 	kETCORR18= 18.0;	// Et SPLIT for 2nd CES CLUSTER CUT

//loose photon id cut values
const static int 	kL_CENTRAL 	= 0; 		// CENTRAL
const static double	kL_ETCORR	= 30.0;
const static double	kL_XCES		= 21.0;
const static double	kL_ZCES_MIN	= 9.0;
const static double	kL_ZCES_MAX	= 230.0;
const static double	kL_HADEM	= 0.125;
const static double	kL_ISOETCORR_MF_LTETC20 = 0.15;	//ISOETCORR Multiplicative factor when Etc<20
const static double	kL_ISOETCORR_MF_GTETC20 = 3.0;	//when Etc>20
const static double	kL_TRKPT_MF	= 0.25;		// TRKPT Multiplicative factor (*EtCorr)
const static double	kL_TRKISO	= 5.0;

//tight photon id cut values

const static int 	kT_CENTRAL 	= 0; 		// CENTRAL
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
const static int 	kT_N3DTRACKS	= 1;
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
bool PrePass()
bool PassGoodRun()
bool PassTrigger()
bool PassZvCut()
int Class12Vertices()
double BestVetex_Z()
bool FindZ(TL&, TL&)
bool FindW(TL&)

void W_dataEstimate()
void WEleIDcutEff();
void ZEleIDcutEff();
void W_mcEstimate()

void TagCoversions()
int FiducialLepton()
MatchPhotonToHEPGLevel(pho)
InitForSashaCode()
SetPhotonList

CEMPhoEleIDcut(pho,event, zvx)
TightPhotonIDcut(pho)
LoosePhotonIDcut(pho)

FillDataBlock(ientry)


FillPhotonPlots(SuperPhoton_t)
FillElectronPlots(SuperPhoton_t)
FillMetPlots

TFolder* GetHistoFolder
BookGeneralHistograms
BookMetHistograms
BookPhotonHistograms
BookElectronHistograms
BookHEPGMatchHistograms

*/



ClassImp(TGammaJetsInit)

//_____________________________________________________________________________
TGammaJetsInit::TGammaJetsInit(const char* name, const char* title):
  TStnModule(name,title), qMc(false)
{
	std::cout << "Hello I am TGammaJetsInit module" << std::endl;
}

//_____________________________________________________________________________
TGammaJetsInit::~TGammaJetsInit() {
}

//_____________________________________________________________________________
void TGammaJetsInit::SaveHistograms() {
}

//_____________________________________________________________________________
void TGammaJetsInit::BookHistograms() {

	DeleteHistograms();
	BookGeneralHistograms();
	BookMetHistograms();
	BookElectronHistograms();
	//BookPhotonHistograms();
	//BookHEPGMatchingHistograms();
	BookBackgroundHistograms();
	BookWeleIDeffHistograms();
	BookZeleIDeffHistograms();
}


//_____________________________________________________________________________
int TGammaJetsInit::BeginJob()
/*{{{*/
{
	std::cout << " BEGIN JOB " << std::endl;
				// register the data block, why?  to unpack them!

	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fCprDataBlock 	= (TCprDataBlock*)     RegisterDataBlock("CprDataBlock","TCprDataBlock");
	fJetBlock 		= (TStnJetBlock*)      RegisterDataBlock("JetBlock","TStnJetBlock");
	fTrackBlock 	= (TStnTrackBlock*)    RegisterDataBlock("TrackBlock","TStnTrackBlock");
	fPhotonBlock 	= (TStnPhotonBlock*)   RegisterDataBlock("PhotonBlock","TStnPhotonBlock");
	fTriggerBlock 	= (TStnTriggerBlock*)  RegisterDataBlock("TriggerBlock","TStnTriggerBlock");
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");

	BookHistograms();

	// initialize vars
	counter.evtsPassTrigger = 0;
	counter.evtsPass25Trigger = 0;
	counter.evtsPass50Trigger = 0;
	counter.evtsPass70Trigger = 0;
	counter.evtsPassGoodrun = 0;
	counter.dWe_detWs = 0;
	counter.dZe_detZs = 0;
	counter.found_hepg_Ws = 0;
	counter.evtsPassVertexCut = 0;
	counter.evtsPassZvCut = 0;

	counter.photons = 0;
	counter.electrons = 0;
	counter.none = 0;
	counter.HepgMatchFailFiducial = 0;
	counter.evtsPassModule = 0;
	counter.evtsRunOver = 0;

	counter.ideff_hepgWs = 0;				//__ hepg W events
	counter.ideff_hepgWenus = 0;			//__ W that decayed to electron+ nutrino
	counter.ideff_hepgWs_passFid = 0;	//__ Wenu, electron found in central region
	counter.ideff_hepgWs_passEt30 = 0;	//__  that electron having Et > 30
	counter.ideff_passGoodrun = 0;		//__ from the events with a W, number of events pass good run
	counter.ideff_passVertex = 0;		//__ then pass 1 good vertex
	counter.ideff_passZv = 0;				//__ then pass Zv cut
	counter.ideff_passMetcut = 0;
	counter.ideff_match2OfflineObject = 0;
	counter.ideff_detWs = 0;


	counter.zideff_hepgZs = 0;				//__ hepg Z events
	counter.zideff_hepgZees = 0;			//__ Z that decayed to electron+ nutrino
	counter.zideff_hepgZs_passFid = 0;	//__ Zee, electrons found in central region
	counter.zideff_hepgZs_passEt30 = 0;	//__  that electron having Et > 30
	counter.zideff_passGoodrun = 0;		//__ from the events with a Z, number of events pass good run
	counter.zideff_passVertex = 0;		//__ then pass 1 good vertex
	counter.zideff_passZv = 0;
	counter.zideff_passZv = 0;				//__ then pass Zv cut
	counter.zideff_preZs = 0;				//__ then pass Zv cut
	counter.zideff_match2OfflineObject = 0;
	counter.zideff_detZs = 0;
	counter.zideff_passConv = 0;



	currRun = 0;
	prevRun = 0;			//cannot init with 1st run number here. fHeaderBlock returns -1
	
	sumWeventMet = 0.;
	sumZeventMet = 0.;

	if (goodRunListFile.Length()<=0)
		goodRunListFile = "goodrun_v17_pho_02.txt";

	goodrun.Read(goodRunListFile.Data());

	trigbits.SetTriggerBlock(fTriggerBlock);

  return 0;
}
/*}}}*/

//_____________________________________________________________________________
int TGammaJetsInit::BeginRun()
{

	currRun = fHeaderBlock->RunNumber();
	std::cout << " BEGIN RUN " << currRun << std::endl;
	
	RunAvg_t ra;
	ra.runNumber = prevRun;
	ra.events = counter.evtsRunOver;
	ra.sumWMet = sumWeventMet; 
	ra.sumZMet = sumZeventMet; 
	ra.Wcount = counter.dWe_detWs;
	ra.Zcount = counter.dZe_detZs;
	RunAvg.push_back(ra); 		//___ first element is all zeros!
	sumWeventMet = 0.;
	sumZeventMet = 0.;
	prevRun = currRun;

  
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
int TGammaJetsInit::Event(int ientry)
{
	SetPassed(0);
	counter.evtsRunOver++;
  //std::cout << "\n\n>>>>>>>>>>>>> entry:"<< ientry << " IN TGI EVENT LOOP\n";

	//_____________________________ REQUIRED Thingies 
  FillDataBlocks(ientry);
  
	if (PrePass()) {
		SetPassed(1);
		counter.evtsPassModule++;
		Cleanup();
		SetPhotonList();
		TagConversions();
		//InitForSashaCode(ientry);
		//if (qMc) W_mcEstimate();
		//W_dataEstimate();
		//if (qMc) WEleIDcutEff();
		if (qMc) ZEleIDcutEff();
		//InstLumMinBiasRelation();
	}

	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TGammaJetsInit::Cleanup()
{
  	trig25 = false;
	trig50 = false;
	trig70 = false;

	Nphotons =0;
	phoUncorrected.clear();
	phoCorrected.clear();
	phoHadEm.clear();
	phoEtaDet.clear();
	phoCesX.clear();
	phoCesZ.clear();
	phoIndex.clear();

	superpho.clear();
	FoundGoodPhoton = false;


} //Cleanup

/*-------------------------------------------------------------------*/
// Event pre-selection cuts
// REMOVE EVENTS WITHOUT PHOTONS/JETS
// do not cut on number of jets yet. do it in the TGammaJets mod
/*-------------------------------------------------------------------*/
bool TGammaJetsInit::PrePass()
{
	if (!qMc) {
		if (fPhotonBlock->NPhotons() < 1) return false;
		else return (PassGoodRun() && PassTrigger() && PassZvCut());
	} else {
		//return (PassGoodRun() && PassZvCut());
		return true;
	}
} // PrePass

/*-------------------------------------------------------------------*/
// good run preselection
// remove the bad runs/sections
/*-------------------------------------------------------------------*/
bool TGammaJetsInit::PassGoodRun()
{
  int run = fHeaderBlock->RunNumber();
  int sec = fHeaderBlock->SectionNumber();
	
  bool good = goodrun.Good(run,sec);
	
  if (!good) {
    return false;
  } else {
    counter.evtsPassGoodrun++;
    return true;
  }
} //PassGoodRun

/*-------------------------------------------------------------------*/
// trigger preselection
/*-------------------------------------------------------------------*/
bool TGammaJetsInit::PassTrigger()
{
	if (!qMc) { //___________________ pypho22_dfc MC sample does not have trigger info
		trigbits.Event();
    	trig25 = trigbits.Pho25Iso();
	   trig50 = trigbits.Pho50();
   	trig70 = trigbits.Pho70();

	    if (!(trig25 || trig50 || trig70) ) {
   		return false; //____ use all 3 triggers. may see something at hight pt.
	    } else {
      	counter.evtsPassTrigger++;
   	   if (trig25) counter.evtsPass25Trigger++;
	      if (trig50) counter.evtsPass50Trigger++;
      	if (trig70) counter.evtsPass70Trigger++;
      	return true;
   	}
	}
}// PassTrigger


/*-------------------------------------------------------------------*/
// Z vertex preselection
// cut on Z position on the vertex,
// so we are within the acceptance of the silicon.
/*-------------------------------------------------------------------*/
bool TGammaJetsInit::PassZvCut()
{
	int nvtx = Class12Vertices();
	//assert(fVertexBlock->NVertices() ==  nvtx );
	
/*	if (fVertexBlock->NVertices() < 1) {
		std::cout << "Nvertices = " << fVertexBlock->NVertices() << "   Class12=" << nvtx << std::endl;
		TPrintModule* tp = new TPrintModule(fHeaderBlock);
		tp->Print(fGenpBlock);
		delete tp;
	}
*/

	if (nvtx < 1) {
		return false;
	} else {
		counter.evtsPassVertexCut++;
		double best_z = BestVertex_Z(); 
    	GeneralPlot.vertexZ->Fill(best_z);
		if (fabs(best_z) < 60.) {
	  		counter.evtsPassZvCut++;
			return true;
  		} else {
			return false;
		}
	}
} //PassZvCut

/*-------------------------------------------------------------------*/
// for cut on vetices: require 1 vertex now. can remove some
// background by doing so.
/*-------------------------------------------------------------------*/
int TGammaJetsInit::Class12Vertices()
{
 	int nvtx = 0;
	double zvx  = 0;

	int Nvtx = fVertexBlock->NVertices();
	GeneralPlot.Nvertices->Fill(Nvtx);

	for (int ivtx = 0; ivtx < Nvtx; ++ivtx) {
		TStnVertex* vert = fVertexBlock->Vertex(ivtx);
		if (vert->VClass() >= 12) ++nvtx;
	}
  
  return nvtx;
} //Class12Vertices

/*-------------------------------------------------------------------*/
// z-position cut on vertex
/*-------------------------------------------------------------------*/
double TGammaJetsInit::BestVertex_Z()
{
  //double sumpt = 0;
  //TStnVertex* bestvert;
  double zvx =0;
  if (fVertexBlock->GetBestVertex(12,1) != NULL) {
    zvx = fVertexBlock->GetBestVertex(12,1)->Z();
  }


  /*for ( int ivtx = 0; ivtx < fVertexBlock->NVertices(); ++ivtx) {
    TStnVertex* ver = fVertexBlock->Vertex(ivtx);
    if (ver->SumPt() > sumpt) {
    bestvert = ver;
    }
    }
    if (bestvert) { std::cout << "found vertex\n"; }
    else { std::cout <<"no vertex\n"; }
    if (bestvert != NULL) return bestvert->Z();
    else return 0.0;
  */
  
  return zvx;
} //BestVertex_Z


/*-------------------------------------------------------------------*/
// look for a Z
// 1. find 2 photons passing Photon like Electron ID cuts
/*-------------------------------------------------------------------*/
bool TGammaJetsInit::FindZ(int& index1, int& index2)
{
	std::vector<int> index;

	if (superpho.size() < 2) return false;

	for (int i = 0; i < superpho.size(); i++) {
		if (index.size() < 2) {
	   	if ((superpho[i].Conversion == 0) && (superpho[i].PhotonLikeID == 0)) {
				index.push_back(i);		
			}
		} else break;
	}

	if (index.size() == 2) {
		index1 = index[0];
		index2 = index[1];
		return true;
	} 	else return false;
}

/*-------------------------------------------------------------------*/
// look for a W
// 1.look for a photon passing Photon like Electron ID cuts
/*-------------------------------------------------------------------*/
bool TGammaJetsInit::FindW(int& index)
{
	bool foundW = false;
	for (int i = 0; i < superpho.size(); i++) {
   	if ((superpho[i].Conversion == 0) && (superpho[i].PhotonLikeID == 0)) {
			std::cout << " FindW index = " << i << std::endl;
			index = i;
			foundW = true;
      	break;
		}
	}
	return foundW;
}

/*-------------------------------------------------------------------*/
// W background estimate from data
// 1. remover conversions
// 2. find photons passing EleLikeIdCUTS
// 3. calculate W/Zs
/*-------------------------------------------------------------------*/
void TGammaJetsInit::W_dataEstimate()
{
	std::vector<int> index; 
	
	for (int i=0; i < superpho.size(); i++) {
   	if ((superpho[i].PhotonLikeID == 0) && (superpho[i].Conversion == 0)) {
			index.push_back(i);
		}
	}

	
	switch (index.size()) {
		case 2:
		{	//_________________ Z event, do not count these as W events
			counter.dZe_detZs++;
			sumZeventMet += fMetBlock->Met(0);
			FillZplots(index[0], index[1]);
			break;
		}
		
		case 1: 
		{	//_________________ W candidate
			if (fMetBlock->Met(0) >20.) {
				counter.dWe_detWs++;
				sumWeventMet += fMetBlock->Met(0);
				FillWplots(index[0]);
			}
			break;
		}


	} //switch

}

void TGammaJetsInit::FillWplots(int index)
{
	SuperPhoton_t sp = superpho[index];
	
	BackgroundPlot.Wele_EtCorr->Fill(sp.EtCorr);
	BackgroundPlot.Wele_HadEm->Fill(sp.HadEm);
	BackgroundPlot.Wele_IsoEtCorr->Fill(sp.IsoEtCorr);
	BackgroundPlot.Wele_TrkPt->Fill(sp.TrkPt);

	double Met_0 = fMetBlock->Met(0);
	BackgroundPlot.W_Met->Fill(Met_0);
	if (Met_0 != 0) BackgroundPlot.W_eleEtMetratio->Fill(Met_0/sp.EtCorr);	//ele et is measurement is accurate than Met

	TLorentzVector ele_vec = *(sp.pho->Momentum());
	TLorentzVector eleT_vec(0,0,0,0);	//transverse component
	eleT_vec.SetPxPyPzE(ele_vec.Px(),ele_vec.Py(),0,ele_vec.E() * sp.pho->SinTheta());
	TLorentzVector met_vec(0,0,0,0);
	met_vec.SetPxPyPzE(Met_0 * TMath::Cos(fMetBlock->MetPhi(0)),Met_0 * TMath::Sin(fMetBlock->MetPhi(0)),0,Met_0);
	
	TLorentzVector sum = eleT_vec + met_vec;
	BackgroundPlot.W_Mass->Fill(sum.M());
	BackgroundPlot.W_Pt->Fill(sum.Perp());
	
	BackgroundPlot.W_Sumet0->Fill(fMetBlock->Sumet(0));
	
}


void TGammaJetsInit::FillZplots(int index1, int index2)
{
	SuperPhoton_t sp1 = superpho[index1];
	SuperPhoton_t sp2 = superpho[index2];
	
	BackgroundPlot.Zele1_EtCorr->Fill(sp1.EtCorr);
	BackgroundPlot.Zele1_HadEm->Fill(sp1.HadEm);
	BackgroundPlot.Zele1_IsoEtCorr->Fill(sp1.IsoEtCorr);
	BackgroundPlot.Zele1_TrkPt->Fill(sp1.TrkPt);
	BackgroundPlot.Zele2_EtCorr->Fill(sp2.EtCorr);
	BackgroundPlot.Zele2_HadEm->Fill(sp2.HadEm);
	BackgroundPlot.Zele2_IsoEtCorr->Fill(sp2.IsoEtCorr);
	BackgroundPlot.Zele2_TrkPt->Fill(sp2.TrkPt);
	BackgroundPlot.Z_eleEtratio->Fill(sp2.EtCorr/sp1.EtCorr);

	BackgroundPlot.Z_Met->Fill(fMetBlock->Met(0));

	TLorentzVector ele1_vec = *(sp1.pho->Momentum());
	TLorentzVector ele2_vec = *(sp2.pho->Momentum());
	TLorentzVector sum = ele1_vec + ele2_vec;
	BackgroundPlot.Z_Mass->Fill(sum.M());
	BackgroundPlot.Z_Pt->Fill(sum.Perp());

	BackgroundPlot.Z_Sumet0->Fill(fMetBlock->Sumet(0));
}





//--------------------------------------------------------------	
// measure eff of CEM PHOTON LIKE ELE ID CUTS for Zee 
/*
HEPG LEVEL
pick events passing good run.
1. find hepg Z
2. find its decay e's 
3. check if both are fiducial
4. check Et>30
4. if so keep their vectors
DETECTOR LEVEL
5. apply Zv cut
5. find matching clusters to both hepg electrons
6. check if both clusters pass ID CUTs

*/
//--------------------------------------------------------------	
void TGammaJetsInit::ZEleIDcutEff()
{
	int Nparticles = fGenpBlock->NParticles();
	TPrintModule p(fHeaderBlock);
	int mom_id = 0;
	TLorentzVector helevec(0,0,0,0), pv(0,0,0,0);	// temp ele vec and production vertex vector for fid cut
	std::vector<TLorentzVector> elevec;					// two electron 4-vectors
	std::vector<TLorentzVector> elepv;					// two electron productions vectors
	int ndau =0;	// number of daughter from Zee = 2

	TGenParticle *par, *mom, *dau;
	
	for (int i = 0 ; i < Nparticles ; i++) {	//_____ not much diff from running over all. my match is mostly within the fist 10 objects
		par = fGenpBlock->Particle(i);
		int im = par->GetFirstMother();

		if (im >= 0) {	//_________________________________im=-1 means no mother for imcoming particles
			mom = fGenpBlock->Particle(im);

			if (mom != 0) {
				int par_id   = par->GetPdgCode();
				int par_stat = par->GetStatusCode();
				if (abs(par_id) == 23 && (par_stat == 2)) { 	//found a Z 
					ZEleIDeffPlot.Zeles_IDeff->Fill(0);
					TLorentzVector zvec;
					par->Momentum(zvec);
					ZEleIDeffPlot.hepgZpt->Fill(zvec.Pt());
					
					ZEleIDeffPlot.hepgZmass->Fill(zvec.M());
					if (zvec.M() < 66. || zvec.M() > 116.) return; 
					std::cout << "Z mass=" << zvec.M() << std::endl;
					counter.zideff_hepgZs++;
					int p1 = par->GetFirstDaughter();
					int pl = par->GetLastDaughter();
					std::cout << "p1,pl=" << p1 << "," << pl << std::endl;
					for (int i= p1 ; i <= pl || ndau == 2; i++) {
						if (ndau != 2) {
							dau = fGenpBlock->Particle(i);
							int dau_id = dau->GetPdgCode();
							int dau_stat = dau->GetStatusCode();
							
							TLorentzVector temp;
							dau->Momentum(temp);
							if ( (abs(dau_id) == 11) && (dau_stat == 1)) {
								ndau++;
								dau->Momentum(helevec);
								dau->ProductionVertex(pv);
								elevec.push_back(helevec);	
								elepv.push_back(pv);	
							} //if (electron)
						} else {
							break;
						}// if (ndau)
						
					} //for
					
					break;	//W is found, but did not decay in to enu or failed Fid cut, so quit
				} 
			} //if W
		} // if
	} //for
	
	if (ndau == 2) {
		counter.zideff_hepgZees++;
		ZEleIDeffPlot.Zeles_IDeff->Fill(1);
	} else return;
	

	if (PassGoodRun()) {
		counter.zideff_passGoodrun++;
		ZEleIDeffPlot.Zeles_IDeff->Fill(4);
	} else return;
	
	if (Class12Vertices() >= 1) {
		counter.zideff_passVertex++;
		ZEleIDeffPlot.Zeles_IDeff->Fill(5);
	} else return;

	if (fabs(BestVertex_Z()) < 60) {
		counter.zideff_passZv++;
		ZEleIDeffPlot.Zeles_IDeff->Fill(6);
	} else return;

/*	int np = 0;
	for (int i = 0; i < fPhotonBlock->NPhotons(); i++) {
		TStnPhoton *pho = fPhotonBlock->Photon(i);
		if (pho->Etc() >20) {
			np++;
		}
	}
	if (np ==2 ) counter.zideff_preZs++;
*/

	int ne = 0;
	for (int i=0; i < superpho.size(); i++ ) {
  		if ( (superpho[i].Conversion == 0) && (superpho[i].PhotonLikeID == 0) ) {
			ne++;
		}
	}

	if (ne ==2) counter.zideff_detZs++;


} //ZEleIDcutEff













//--------------------------------------------------------------	
// measure eff of CEM PHOTON LIKE ELE ID CUTS for Wenu
//--------------------------------------------------------------	
void TGammaJetsInit::WEleIDcutEff()

{
	int Nparticles = fGenpBlock->NParticles();

	int mom_id = 0;
	TLorentzVector helevec(0,0,0,0);		// hepg stable electron
	double heleEt = 0;						// hepg electron et
	TLorentzVector pv(0,0,0,0);			// production vertex vector for fid cut
	bool foundWelectron = false;

	TGenParticle *par, *mom, *dau;
	
	for (int i = 0 ; i < Nparticles ; i++) {	//_____ not much diff from running over all. my match is mostly within the fist 10 objects
		par = fGenpBlock->Particle(i);
		int im = par->GetFirstMother();

		if (im >= 0) {	//_________________________________im=-1 means no mother for imcoming particles
			mom = fGenpBlock->Particle(im);

			if (mom != 0) {
				int par_id   = par->GetPdgCode();
				int par_stat = par->GetStatusCode();
				if (abs(par_id) == 24 && (par_stat == 2)) { 	//found a W 
					EleIDeffPlot.Wele_IDeff->Fill(0);
					counter.ideff_hepgWs++;
					int p1 = par->GetFirstDaughter();
					int pl = par->GetLastDaughter();

					for (int i= p1 ; i <= pl; i++) {
						dau = fGenpBlock->Particle(i);
						int dau_id = dau->GetPdgCode();
						int dau_stat = dau->GetStatusCode();
						if ( (abs(dau_id) == 11) && (dau_stat == 1)) {
							EleIDeffPlot.Wele_IDeff->Fill(1);
							counter.ideff_hepgWenus++;
							dau->Momentum(helevec);
							dau->ProductionVertex(pv);
							heleEt = helevec.Energy() * TMath::Sin(helevec.Theta());
							foundWelectron = true;
							break;
						} //if electron
					} //for
					
					break;	//W is found, but did not decay in to enu or failed Fid cut, so quit
				} 
			} //if W
		} // if
	} //for
	
	if (!foundWelectron) return;

	if (PassGoodRun()) {
		counter.ideff_passGoodrun++;
		EleIDeffPlot.Wele_IDeff->Fill(3);
	} else return;
	
	if (Class12Vertices() >= 1) {
		counter.ideff_passVertex++;
		EleIDeffPlot.Wele_IDeff->Fill(4);
	} else return;

	if (fabs(BestVertex_Z()) < 60) {
		counter.ideff_passZv++;
		EleIDeffPlot.Wele_IDeff->Fill(5);
	} else return;

	if (fMetBlock->Met(0) > 20) {
		counter.ideff_passMetcut++;
	} else return;
	
	int ne = 0;	
	for (int i=0; i< superpho.size(); i++) {	
  		if ((superpho[i].Conversion == 0) && (superpho[i].PhotonLikeID == 0)) {
			ne++;
			EleIDeffPlot.Wele_IDeff->Fill(7);
		} //for
	} //if

	if (ne == 1) counter.ideff_detWs++;

} //WEleIDcutEff



	
/*-------------------------------------------------------------------*/
//  W background estimate from W mc
//	1. remove conversion 
// 2. find a photon passing EleLikeIDcuts
// 3. match to hepg
// 4. see if mathces to an electron
// 5. check if its fiducial
// 6. Et cuts
// 7. see if it is from a W

//for denominator: just try to find a match
// numerator: see if that match is from a W

//i want to OR for muons(13) and tau(15) for decays of W
/*-------------------------------------------------------------------*/
void TGammaJetsInit::W_mcEstimate()
{
	int eleIndex = -99;

	if (FindW(eleIndex)) {
		std::cout << " W estimate index = " << eleIndex << std::endl;

		FillMetPlots();
		FillElectronPlots(superpho[eleIndex]);
		FillWmassPlots(superpho[eleIndex]);
		
		double eleEt = superpho[eleIndex].EtCorr;
		TLorentzVector elevec(superpho[eleIndex].corvec);
	
		//___________________________________________________ now find a match in hepg
		int Nparticles = fGenpBlock->NParticles();
	
		int mom_id = 0;
		int np = Nparticles;
	
		//if (np > 100) np =100; 		//________________________ i can find what i am looking for withing the 20 iterations

		for (int i = 0 ; i < np ; i++) {		//_______________ not much diff from running over all. my match is mostly within the fist 10 objects
			TGenParticle* par = fGenpBlock->Particle(i);
			int im = par->GetFirstMother();

			if (im >= 0) {	//_________________________________im=-1 means no mother for imcoming particles
				TGenParticle* mom = fGenpBlock->Particle(im);

				if (mom != 0) {
					mom_id = mom->GetPdgCode();
					int mom_stable = mom->GetStatusCode();
					int par_id 		= par->GetPdgCode();
					int par_stable = par->GetStatusCode();

					TLorentzVector parvec(0,0,0,0);
					par->Momentum(parvec);

					double DelR = parvec.DeltaR(elevec);

					TLorentzVector pv(0,0,0,0);
					par->ProductionVertex(pv);

					if (DelR < 0.1) {
						if (FiducialLepton(parvec,pv.Z())) { 
							if ( (fabs(par_id) == 11 ) && (par_stable == 1) ) { //2= when W is about to decay 3= when w is produced
								counter.found_hepg_Ws++;
								HEPGMatchPlot.Et_ratio->Fill(eleEt / (par->Energy() * TMath::Sin(par->Theta())));
								break;
							}
						} else counter.HepgMatchFailFiducial++;
					} 
				} //if
				delete mom;
			} // if
			delete par;
		} //for
	} // if

} // W_mcEstimate


/*-------------------------------------------------------------------*/
//tag conversion electrons for W calculations
/*-------------------------------------------------------------------*/
void TGammaJetsInit::TagConversions()
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
bool TGammaJetsInit::IsConversionElectron(TStnElectron* ele, 
				     double& minsep, double& mindt, double& radius) 
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


/*-------------------------------------------------------------------*/
// calculate seperation between the twoo traks for the conversion electrons
/*-------------------------------------------------------------------*/
double TGammaJetsInit::ConversionSep(TStnTrack* t0, TStnTrack* t1, double & rconv) {

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



/*-------------------------------------------------------------------*/
// check if the photon is fidual before matching to HEPG Level for W estimate.
/*-------------------------------------------------------------------*/
int TGammaJetsInit::FiducialLepton(TLorentzVector vec_old, double z_new) {

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

/*-------------------------------------------------------------------*/
// Match a photon to HEPG Level and return mom's PDG code
/*-------------------------------------------------------------------*/
int TGammaJetsInit::MatchPhotonToHEPGLevel(TStnPhoton* pho)
{

  // need to change this to use super photon!
	unsigned int found_matches = 0;
	bool found_match = false;
	double pho_eta = pho->Eta();		//_ use event eta. not detEta
	double pho_phi = pho->Phi();
	double pho_EtCorr = pho->ECorr() * pho->SinTheta();
	
	int Nparticles = fGenpBlock->NParticles();
	if (Nparticles < 1) {
	  return false;
	}
	
	
	int mom_id;
	for (int i = 0 ; i < Nparticles; i++) {
		TGenParticle* par = fGenpBlock->Particle(i);
		int im = par->GetFirstMother();

		if(im>=0) {	//_________________________________im=-1 means, no mother for imcoming particles
			TGenParticle* mom = fGenpBlock->Particle(im);
			if (mom != 0) {
				mom_id = mom->GetPdgCode();
				int mom_stable = mom->GetStatusCode();
				int mom_number = mom->Number();
				//double mom_charge = mom->Charge();  // give an error in TGenParticle.hh line 56??
				int par_id = par->GetPdgCode();
				int par_stable = par->GetStatusCode();
				int par_number = par->Number();
				//double par_charge = par->Charge();
				
				if ( par->GetStatusCode()==1 ) {
					
						float gen_eta = par->Eta();
						float gen_phi = par->Phi();
						TMathModule math;
						found_match = math.MatchEtaPhi(pho_phi,pho_eta,gen_phi,gen_eta,0.1); 		 // match with this precision
						
						if (found_match) {
							double delEta = math.GetDelEta(pho_eta,gen_eta);
							double delPhi = math.GetDelPhi(pho_phi,gen_phi);
							double delR   = math.GetDelR(pho_phi,pho_eta,gen_phi,gen_eta);
							double par_Et = par->Energy() * TMath::Sin(par->Theta());
							//HEPGPhoMatchPlot.DelEta->Fill(delEta);	
							//HEPGPhoMatchPlot.DelPhi->Fill(delPhi);	
							//HEPGPhoMatchPlot.DelR->Fill(delR);
							//HEPGPhoMatchPlot.HEPGObj_Et->Fill(par->Energy() * TMath::Sin(par->Theta()));
							//HEPGPhoMatchPlot.Et_ratio->Fill(pho_EtCorr / (par->Energy() * TMath::Sin(par->Theta())) );
							//HEPGPhoMatchPlot.HEPGObj_Eta->Fill(par->Eta());
							double match_et_ratio = par_Et/pho_EtCorr;
						      
							break;							

						}
				}
			}
		}
	}

	return mom_id;

} //MatchPhotonToHepg



/*-------------------------------------------------------------------*/
// generate stuff needed by sasha's code
/*-------------------------------------------------------------------*/
void TGammaJetsInit::InitForSashaCode(int ientry)
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




/*-------------------------------------------------------------------*/
//	Sorts the Photon Block according to EtCorr, from highest to lowest
// & fill the photon info to my SuperPhoton list
/*-------------------------------------------------------------------*/
void TGammaJetsInit::SetPhotonList()
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
					//delete temp;

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
	bool FoundElectron = false;
	bool FoundPhoton = false;
	
	for (int i=0; i < photon.size(); i++) {
		SuperPhoton_t sp;
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

		int fillto = sp.PhotonLikeID;
		if (sp.PhotonLikeID ==0) fillto = 17;
		for (int j=0; j < fillto; j++) EleIDeffPlot.CEMPHOEleIDeff->Fill(j);

		superpho.push_back(sp);

		if (!FoundGoodPhoton) {
			if ( sp.TightID ==0 ) {
				counter.photons++;
				FoundPhoton = true;
			}
		}
		if (!FoundElectron) {	
		 	if (sp.PhotonLikeID == 0) {
				counter.electrons++;
				FoundElectron = true;
			}
		}
		//else if (sp.PhotonLikeID == 0) counter.electrons++;
		//else counter.none++;		//__ what ever other stuff is?
	
	}

} //SetPhotonList



/*****************  PHOTON LIKE ID CUT FUNCTION ************************/
//_____ function returns 0 if electrons passes cuts
//      requires zvx of highest Pt class12 vertex
unsigned int TGammaJetsInit::CEMPhoEleIDcut(TStnPhoton* Pho, TStnEvent* event, double zvx)
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




/*-------------------------------------------------------------------*/
// Tight Photon ID CUTS
/*-------------------------------------------------------------------*/
unsigned int TGammaJetsInit::TightPhotonIDcut(TStnPhoton *pho)
{
	unsigned int IDWord = 0x0;
	unsigned int IDWord1 = 0x0;

	int   detector 	= pho->Detector();
	float ECorr 		= pho->ECorr(); //total corrected energy // this is what i want to cut on
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
	if (EtCorr < 30 						) IDWord |= kEtCorr7_T; // has no effect! looser cut
	if (fabs(XCes) > 21					) IDWord |= kXCes_T;
	if (fabs(ZCes) < 9 ||  fabs(ZCes) > 230) IDWord |= kZCes_T;  // same as loose cuts
	if ( !((HadEm < 0.125) || (HadEm < (0.055 + 0.00045 * ECorr))) 	) IDWord |= kHadEm_T;
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

	assert(IDWord == IDWord1);
	return IDWord;

} // TightPhotonCuts


/*-------------------------------------------------------------------*/
//	Loose Photon ID CUTS
/*-------------------------------------------------------------------*/
unsigned int TGammaJetsInit::LoosePhotonIDcut(TStnPhoton *pho)
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


/*-------------------------------------------------------------------*/
void TGammaJetsInit::FillDataBlocks(int ientry)
{
	fTriggerBlock->GetEntry(ientry);
	fPhotonBlock->GetEntry(ientry);
	fJetBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
	fVertexBlock->GetEntry(ientry);
	fTrackBlock->GetEntry(ientry);
	fGenpBlock->GetEntry(ientry);
}

/*-------------------------------------------------------------------*/
// fill W mass plots
/*-------------------------------------------------------------------*/
void TGammaJetsInit::FillWmassPlots(SuperPhoton_t& sp)
{
}
/*-------------------------------------------------------------------*/
// fill photons plots
/*-------------------------------------------------------------------*/
void TGammaJetsInit::FillPhotonPlots(SuperPhoton_t& sp)
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
	PhotonPlot.HadTDCtime->Fill(sp.HadTDCtime);
	PhotonPlot.Ces2Strip->Fill(sp.CesStripE2);
}
	
/*-------------------------------------------------------------------*/
// fill electron plots
/*-------------------------------------------------------------------*/
void TGammaJetsInit::FillElectronPlots(SuperPhoton_t& sp)
{
	ElectronPlot.Detector->Fill(sp.Detector);
	ElectronPlot.EtCorr->Fill(sp.EtCorr);
	ElectronPlot.XCes->Fill(sp.XCes);
	ElectronPlot.ZCes->Fill(sp.ZCes);
	ElectronPlot.HadEm->Fill(sp.HadEm);
	ElectronPlot.IsoEtCorr->Fill(sp.IsoEtCorr);
	ElectronPlot.Chi2Mean->Fill(sp.Chi2Mean);
	ElectronPlot.N3d->Fill(sp.N3d);
	ElectronPlot.TrkPt->Fill(sp.TrkPt);
	ElectronPlot.TrkIso->Fill(sp.TrkIso);
	ElectronPlot.Ces2Wire->Fill(sp.CesWireE2);
	ElectronPlot.HadTDCtime->Fill(sp.HadTDCtime);
	ElectronPlot.Ces2Strip->Fill(sp.CesStripE2);
}
/*-------------------------------------------------------------------*/
// fill  Missing Et plots
/*-------------------------------------------------------------------*/
void TGammaJetsInit::FillMetPlots()
{
	MetPlot.Met_0->Fill(fMetBlock->Met(0));
	MetPlot.MetPhi_0->Fill(fMetBlock->MetPhi(0));
	MetPlot.MetX_0->Fill(fMetBlock->MetX(0));
	MetPlot.MetY_0->Fill(fMetBlock->MetY(0));
	MetPlot.Sumet_0->Fill(fMetBlock->Sumet(0));
	MetPlot.Sumet_1->Fill(fMetBlock->Sumet(1));
	MetPlot.Sumet_2->Fill(fMetBlock->Sumet(2));
	MetPlot.Metsig->Fill(fMetBlock->MetSig());
	MetPlot.Z0->Fill(fMetBlock->Z0());
}

/*===================================================================*/
//	Create a folder in"Hist" 
/*===================================================================*/
TFolder* TGammaJetsInit::GetHistoFolder(char *name, char* title)
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
void TGammaJetsInit::BookGeneralHistograms()
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
	
/*	
	GeneralPlot.pho_IDbins_pass = new TH1F("pho_IDbins_pass","The Leading Photon Pass Binning",8,0,8);
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(1,"Fail all");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(2,"Pass Loose Pho");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(3,"Pass Tight Pho");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(4,"Pass Pho-Like Ele");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(5,"Pass Loose && Tight");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(6,"Pass Loose && Pho-Like");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(7,"Pass Tight && Pho-Like");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(8,"Pass all");
*/
	GeneralPlot.Nvertices = new TH1F("Nvertices","Number of Class12 vertices per event",10,0,10);
	GeneralPlot.vertexZ = new TH1F("BestVertex_Z","Best Class12 vertex Z position",100,-100,100);

	TFolder* new_folder = GetHistoFolder("General_Plots","General Plots");
	new_folder->Add(GeneralPlot.photon_loose_and_tight_cuts_cumm);	
	new_folder->Add(GeneralPlot.pho_passIDs);
	new_folder->Add(GeneralPlot.Nvertices);
	new_folder->Add(GeneralPlot.vertexZ);
	std::cout << "DONE\n";

} //BookGeneralHistograms

//_____________________________________________________________________________
void TGammaJetsInit::BookMetHistograms()
{
	std::cout << "BOOKING Met Plots...";
	
	MetPlot.Met_0   	= new TH1F("Met_0","Met(0) of W events",100,0,100);
	MetPlot.MetPhi_0  = new TH1F("MetPhi_0","Met(0) Phi of W events",70,0,7);
	MetPlot.MetX_0   	= new TH1F("MetX_0","MetX(0) of W events",50,0,50);
	MetPlot.MetY_0   	= new TH1F("MetY_0","MetY(0) of W events",50,0,50);
	MetPlot.Sumet_0   = new TH1F("Sumet_0","Sumet_0 of W events",100,0,500);
	MetPlot.Sumet_1   = new TH1F("Sumet_1","Sumet_1:Htc of W events",100,0,500);
	MetPlot.Sumet_2   = new TH1F("Sumet_2","Sumet_2:Sumetjet of W events",100,0,500);
	MetPlot.Metsig   	= new TH1F("Metsig","MetSig of W events",50,0,50);
	MetPlot.Z0   		= new TH1F("Z0","Z0 position of W events",50,0,50);

	TFolder* new_folder = GetHistoFolder("Met_Plots"," Missing Et Plots");
	new_folder->Add(MetPlot.Met_0);	
	new_folder->Add(MetPlot.MetPhi_0);	
	new_folder->Add(MetPlot.MetX_0);	
	new_folder->Add(MetPlot.MetY_0);	
	new_folder->Add(MetPlot.Sumet_0);	
	new_folder->Add(MetPlot.Sumet_1);	
	new_folder->Add(MetPlot.Sumet_2);	
	new_folder->Add(MetPlot.Metsig);	
	new_folder->Add(MetPlot.Z0);	
	std::cout << "DONE.\n";
}
//_____________________________________________________________________________
void TGammaJetsInit::BookBackgroundHistograms()
{
	std::cout << "BOOKING Background Plots...";

	BackgroundPlot.Wele_EtCorr   	= new TH1F("Wele_EtCorr","W Electron Et",100,0,100);
	BackgroundPlot.Wele_HadEm		= new TH1F("Wele_HadEm","W Electron HadEm",50,0,1);
	BackgroundPlot.Wele_IsoEtCorr	= new TH1F("Wele_IsoEtCorr","W Electron IsoEtCorr",140,-4,10);
	BackgroundPlot.Wele_TrkPt 		= new TH1F("Wele_Trkpt","W Electron TrkPt",200,0,100);
	BackgroundPlot.W_eleEtMetratio= new TH1F("W_eleEtMetratio","W: Met to Electron Et Ratio",40,0,2);
	BackgroundPlot.W_Mass 			= new TH1F("W_mass","W Transverse Invriant Mass",120,0,120);
	BackgroundPlot.W_Pt 				= new TH1F("W_Pt","W Pt",100,0,100);
	BackgroundPlot.W_Sumet0 		= new TH1F("W_Sumet0","W Events Sumet(0)",100,0,500);
	BackgroundPlot.W_Met 			= new TH1F("W_Met","Met for W events",140,0,70);

	BackgroundPlot.Zele1_EtCorr   	= new TH1F("Zele1_EtCorr","Z Leading Electron Et",100,0,100);
	BackgroundPlot.Zele1_HadEm			= new TH1F("Zele1_HadEm","Z Leading Electron HadEm",50,0,1);
	BackgroundPlot.Zele1_IsoEtCorr	= new TH1F("Zele1_IsoEtCorr","Z Leading Electron IsoEtCorr",140,-4,10);
	BackgroundPlot.Zele1_TrkPt 		= new TH1F("Zele1_Trkpt","Z Leading Electron TrkPt",1000,0,100);
	BackgroundPlot.Zele2_EtCorr   	= new TH1F("Zele2_EtCorr","Z 2nd Leading Electron Et",100,0,100);
	BackgroundPlot.Zele2_HadEm			= new TH1F("Zele2_HadEm","Z 2nd Leading Electron HadEm",50,0,1);
	BackgroundPlot.Zele2_IsoEtCorr	= new TH1F("Zele2_IsoEtCorr","Z 2nd Leading Electron IsoEtCorr",140,-4,10);
	BackgroundPlot.Zele2_TrkPt 		= new TH1F("Zele2_Trkpt","Z 2nd Leading Electron TrkPt",200,0,100);
	BackgroundPlot.Z_Mass 				= new TH1F("Z_mass","Z Invriant Mass",120,0,120);
	BackgroundPlot.Z_Pt 					= new TH1F("Z_Pt","Z Pt",100,0,200);
	BackgroundPlot.Z_eleEtratio		= new TH1F("Z_eleEtratio","Z electrons Et ratio (2nd leading to leading) ",100,0,2);
	BackgroundPlot.Z_Sumet0 			= new TH1F("Z_Sumet_0","Z Events Sumet(0)",100,0,500);
	BackgroundPlot.Z_Met 				= new TH1F("Z_Met","Met for Z events",100,0,50);


	TFolder* new_folder = GetHistoFolder("Background_Plots","Background Plots");
	new_folder->Add(BackgroundPlot.Wele_EtCorr);	
	new_folder->Add(BackgroundPlot.Wele_HadEm);	
	new_folder->Add(BackgroundPlot.Wele_IsoEtCorr);	
	new_folder->Add(BackgroundPlot.Wele_TrkPt);	
	new_folder->Add(BackgroundPlot.W_eleEtMetratio);	
	new_folder->Add(BackgroundPlot.W_Mass);
	new_folder->Add(BackgroundPlot.W_Pt);
	new_folder->Add(BackgroundPlot.W_Sumet0);
	new_folder->Add(BackgroundPlot.W_Met);

	new_folder->Add(BackgroundPlot.Zele1_EtCorr);	
	new_folder->Add(BackgroundPlot.Zele2_EtCorr);	
	new_folder->Add(BackgroundPlot.Zele1_HadEm);	
	new_folder->Add(BackgroundPlot.Zele2_HadEm);	
	new_folder->Add(BackgroundPlot.Zele1_IsoEtCorr);	
	new_folder->Add(BackgroundPlot.Zele2_IsoEtCorr);	
	new_folder->Add(BackgroundPlot.Zele1_TrkPt);	
	new_folder->Add(BackgroundPlot.Zele2_TrkPt);	
	new_folder->Add(BackgroundPlot.Z_Mass);	
	new_folder->Add(BackgroundPlot.Z_Pt);
	new_folder->Add(BackgroundPlot.Z_eleEtratio);	
	new_folder->Add(BackgroundPlot.Z_Sumet0);
	new_folder->Add(BackgroundPlot.Z_Met);
	
	std::cout << "DONE.\n";
}
//_____________________________________________________________________________
void TGammaJetsInit::BookPhotonHistograms()
{
	std::cout << "Booking Photon Plots...";

	PhotonPlot.Detector 	= new TH1F("Detector","Wenu MC: Tight Photon Detector Region",3,0,3);
	PhotonPlot.EtCorr 	= new TH1F("EtCorr","Wenu MC: Tight Photon Corrected Et",600,0,120);
	PhotonPlot.XCes 		= new TH1F("XCes","Wenu MC: Tight Photon XCes",640,-32,32);
	PhotonPlot.ZCes 		= new TH1F("ZCes","Wenu MC: Tight Photon ZCes",2000,-250,250);
	PhotonPlot.HadEm 		= new TH1F("HadEm","Wenu MC: Tight Photon HadEm",50,0,0.5);
	PhotonPlot.IsoEtCorr = new TH1F("IsoEtCorr","Wenu MC: Tight Photon IsoCorr",1500,-5,10);
	PhotonPlot.Chi2Mean 	= new TH1F("Chi2Mean","Wenu MC: Tight Photon Chi2Mean (Wire+Strip/2)",800,0,80);
	PhotonPlot.N3d 		= new TH1F("N3d","Wenu MC: Tight Photon N3d",10,0,10);
	PhotonPlot.TrkPt 		= new TH1F("Trkpt","Wenu MC: Tight Photon TrkPt",1000,0,100);
	PhotonPlot.TrkIso 	= new TH1F("TrkIso","Wenu MC: Tight Photon TrkIso",150,0,15);
	PhotonPlot.Ces2Wire 	= new TH1F("Ces2Wire","Wenu MC: Tight Photon CES(2nd) Wire",400,0,40);
	PhotonPlot.Ces2Strip = new TH1F("Ces2Strip","Wenu MC: Tight Photon CES(2nd) Strip",400,0,40);
	PhotonPlot.HadTDCtime = new TH1F("HadTDCtime","Wenu MC: Tight Photon Had TDC time",200,-10,10);

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

	TFolder* new_folder = GetHistoFolder("Photon_Plots","Photon Plots");
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


//_____________________________________________________________________________
void TGammaJetsInit::BookElectronHistograms()
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

	TFolder* new_folder = GetHistoFolder("Electron_Plots","Electron Plots");
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

//_____________________________________________________________________________
void TGammaJetsInit::BookWeleIDeffHistograms()
{
	std::cout << "Booking W ele ID Eff Plots...";
	EleIDeffPlot.MinDelRMatch	= new TH1F("Wele2CloseOffPhoDelR","DelR between W Electron and Closest offline Photon",500,0,5);
	EleIDeffPlot.HepgEleEt		= new TH1F("HepgEleEt","W electron Et",120,0,60);
	EleIDeffPlot.OfflineEleEt	= new TH1F("OfflineEleEt","Matched offline electron Et",120,0,60);
	EleIDeffPlot.OfflineEleTrkpt	= new TH1F("OfflineEleTrkpt","Matched offline electron Track Pt",120,0,120);
	EleIDeffPlot.Et_ratio		= new TH1F("Et_Offele2Hepgele","Matched Offline Electron to HEPG Object Et ratio",40,0,2);
	EleIDeffPlot.Et_ratio->SetYTitle("Events/0.05");
	EleIDeffPlot.Wele_IDeff  	= new TH1F("Wid_eff","W ID eff",15,0,15);
	TAxis* x= EleIDeffPlot.Wele_IDeff->GetXaxis();
	x->SetBinLabel(1,"Ws found");
	x->SetBinLabel(2,"Wenu");
	x->SetBinLabel(3,"Fid");
	x->SetBinLabel(4,"goodrun");
	x->SetBinLabel(5,"vertex");
	x->SetBinLabel(6,"Zv");
	x->SetBinLabel(7,"DelR<0.2");
	x->SetBinLabel(8,"passEleID");
	EleIDeffPlot.CEMPHOEleIDeff 	= new TH1F("EleIDeff","Electron ID eff for W",18,0,18);
	TAxis* x2= EleIDeffPlot.CEMPHOEleIDeff->GetXaxis();
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(1,"PhoCentral");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(2,"PhoEtc30");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(3,"PhoHadEm");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(4,"PhoChi2");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(5,"PhoN3d=0");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(6,"PhoN3d>2");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(7,"PhoPt2");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(8,"FindMatchEle");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(9,"EleHasTrk");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(10,"EleTrkZ");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(11,"EleEoveP");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(12,"PhoCalIso4");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(13,"PhoTrkIso");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(14,"PhoCes2Et");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(15,"PhoXCes");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(16,"PhoZCes");
	EleIDeffPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(17,"Electrons");

	TFolder* new_folder = GetHistoFolder("WeleIDeff_Plots","Electron ID eff Plots");
	new_folder->Add(EleIDeffPlot.HepgEleEt);
	new_folder->Add(EleIDeffPlot.OfflineEleEt);
	new_folder->Add(EleIDeffPlot.OfflineEleTrkpt);
	new_folder->Add(EleIDeffPlot.MinDelRMatch);
	new_folder->Add(EleIDeffPlot.Et_ratio);
	new_folder->Add(EleIDeffPlot.Wele_IDeff);
	new_folder->Add(EleIDeffPlot.CEMPHOEleIDeff);
	std::cout << "DONE.\n";
}


//_____________________________________________________________________________
void TGammaJetsInit::BookZeleIDeffHistograms()
{
	std::cout << "Booking Z ele ID Eff Plots...";
	ZEleIDeffPlot.ele1Et_ratio	= new TH1F("ele1Etratio","Z hepg electron(1) to matched offline photon's Et",50,0,2);
	ZEleIDeffPlot.ele2Et_ratio	= new TH1F("ele2Etratio","Z hepg electron(2) to matched offline photon's Et",50,0,2);
	ZEleIDeffPlot.Zeles_IDeff  	= new TH1F("Zid_eff","Z ID eff",15,0,15);
	ZEleIDeffPlot.hepgZpt  	= new TH1F("hepgZpt","HEPG Z pt",15,0,15);
	ZEleIDeffPlot.hepgZmass  	= new TH1F("hepgZmass","HEPG Z mass",60,0,120);
	TAxis* x= ZEleIDeffPlot.Zeles_IDeff->GetXaxis();
	x->SetBinLabel(1,"Zs found");
	x->SetBinLabel(2,"Z->ee");
	x->SetBinLabel(3,"Fid(both)");
	x->SetBinLabel(4,"eles Et>30");
	x->SetBinLabel(5,"goodrun");
	x->SetBinLabel(6,"vertex");
	x->SetBinLabel(7,"Zv");
	x->SetBinLabel(8,"DelR<0.2");
	x->SetBinLabel(9,"passEleID &Conv");

	TFolder* new_folder = GetHistoFolder("ZeleIDeff_Plots","Z Electrons ID eff Plots");
	new_folder->Add(ZEleIDeffPlot.ele1Et_ratio);
	new_folder->Add(ZEleIDeffPlot.ele2Et_ratio);
	new_folder->Add(ZEleIDeffPlot.Zeles_IDeff);
	new_folder->Add(ZEleIDeffPlot.hepgZpt);
	new_folder->Add(ZEleIDeffPlot.hepgZmass);
	std::cout << "DONE.\n";
}
//_____________________________________________________________________________
void TGammaJetsInit::BookHEPGMatchingHistograms()
{
	std::cout << "Booking HEPG Photon Matching Plots...";
	HEPGMatchPlot.Et_ratio		= new TH1F("Et_OffObj_to_HepgObj","Offline Electron to Matched HEPG Object Et",40,0,2);
	HEPGMatchPlot.Et_ratio->SetYTitle("Events/0.05");

	TFolder* new_folder = GetHistoFolder("HEPGMatching_Plots","Photon HEPG Matching Plots");
	new_folder->Add(HEPGMatchPlot.Et_ratio);
	std::cout << "DONE.\n";
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TGammaJetsInit::EndJob() {

	printf("----- end job: ---- %s\n",GetName());
	if (qMc)	std::cout << "RUN IS ON MC SAMPLE" << std::endl;
	else	std::cout << "RUN IS ON DATA SAMPLE" << std::endl;
	std::cout << "EVENTS RUN OVER ------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "EVENTS Pass Goodrun -------- = " << counter.evtsPassGoodrun << std::endl;
	std::cout << "EVENTS Pass Trigger(25,50,70)= " << counter.evtsPassTrigger << "\t(" << counter.evtsPass25Trigger
																	<<","<< counter.evtsPass50Trigger << "," << counter.evtsPass70Trigger <<")"<< std::endl;
	std::cout << "EVENTS Pass Vertex Cut ----- = " << counter.evtsPassVertexCut << std::endl;
	std::cout << "EVENTS Pass Zv Cut  -------- = " << counter.evtsPassZvCut << std::endl;
	std::cout << "EVENTS PASS THIS MODULE ---- = " << counter.evtsPassModule << std::endl;

	std::cout << " ===================== W/Z estimate from data ======================" << std::endl;
	std::cout << "Found Detector W events ---- = " << counter.dWe_detWs << std::endl;
	std::cout << "Found Detector Z events ---- = " << counter.dZe_detZs << std::endl;
	std::cout << " ===================================================================" << std::endl;

	if (qMc) {
		std::cout << " ======================= W Electron ID eff ===========================" << std::endl;
		std::cout << "Evts with a HEPG W ------------------ = " << counter.ideff_hepgWs << std::endl;
		std::cout << "HEPG W->enu ------------------------- = " << counter.ideff_hepgWenus << std::endl;
		std::cout << "W electron in central region -------- = " << counter.ideff_hepgWs_passFid << std::endl;
		std::cout << "W electron has Et>30 ---------------- = " << counter.ideff_hepgWs_passEt30 << std::endl;
		std::cout << "Det. evt passing goodrun ------------ = " << counter.ideff_passGoodrun << std::endl;
		std::cout << "passing good vertex ----------------- = " << counter.ideff_passVertex << std::endl;
		std::cout << "passing Zv < 60  -------------------- = " << counter.ideff_passZv << std::endl;
		std::cout << "passing Met > 20 -------------------- = " << counter.ideff_passMetcut << std::endl;
		std::cout << "find match to offline photon in DelR  = " << counter.ideff_match2OfflineObject << std::endl;
		std::cout << "Evts with Det W (passing Ele ID cuts) = " << counter.ideff_detWs << std::endl;
		std::cout << " ===================================================================" << std::endl;
		std::cout << " ======================= Z Electron ID eff ===========================" << std::endl;
		std::cout << "Evts with a HEPG Z ------------------ = " << counter.zideff_hepgZs << std::endl;
		std::cout << "HEPG Z->ee -------------------------- = " << counter.zideff_hepgZees << std::endl;
		std::cout << "Z electrons in central region ------- = " << counter.zideff_hepgZs_passFid << std::endl;
		std::cout << "Z electrons has Et>30 --------------- = " << counter.zideff_hepgZs_passEt30 << std::endl;
		std::cout << "Det. evt passing goodrun ------------ = " << counter.zideff_passGoodrun << std::endl;
		std::cout << "passing good vertex ----------------- = " << counter.zideff_passVertex << std::endl;
		std::cout << "passing Zv < 60  -------------------- = " << counter.zideff_passZv << std::endl;
		std::cout << "pretag Zee events(Et>20 electrons)    = " << counter.zideff_preZs << std::endl;
		std::cout << "both electrons are not conversions    = " << counter.zideff_passConv << std::endl;
		std::cout << "Evts with Det Z (both e pass EleIDcu) = " << counter.zideff_detZs << std::endl;
		std::cout << " ===================================================================" << std::endl;
	}


		// add stuff from the last run
		RunAvg_t ra;
		ra.runNumber = currRun;
		ra.events = counter.evtsRunOver;
		ra.sumWMet = sumWeventMet; 
		ra.sumZMet = sumZeventMet; 
		ra.Wcount = counter.dWe_detWs;
		ra.Zcount = counter.dZe_detZs;
		RunAvg.push_back(ra);
	
		Int_t size = RunAvg.size()-1;
		// all this must be double for the TGraph to work!
		Double_t rn[size];
		Double_t evts[size];
		Double_t swmet[size];
		Double_t szmet[size];
		Double_t wc[size];
		Double_t zc[size];
		for (int i =0 ; i < size; i++) {
			rn[i] = .0;
			evts[i] = .0;
			swmet[i] = .0;
			szmet[i] = .0;
			wc[i] = .0;
			zc[i] = .0;

		}

		std::cout << "\n============================================== "<< std::endl;

		std::cout <<"Run#   Events   avgWmet   avgZmet   Wcount   Zcount" << std::endl; 

		for (Int_t i=0; i < RunAvg.size(); ++i) {
			if (i==0) continue;
			rn[i-1] 	= RunAvg[i].runNumber;
			evts[i-1]= RunAvg[i].events - RunAvg[i-1].events;
			wc[i-1] 	= RunAvg[i].Wcount - RunAvg[i-1].Wcount;;
			zc[i-1] 	= RunAvg[i].Zcount - RunAvg[i-1].Zcount;;

			if (wc[i-1] != 0) swmet[i-1] = RunAvg[i].sumWMet / wc[i-1];

			if (zc[i-1] != 0) szmet[i-1] = RunAvg[i].sumZMet /  zc[i-1];

			std::cout << rn[i-1] << "\t";
			std::cout << evts[i-1] << "\t";
			std::cout << swmet[i-1] << "\t";
			std::cout << szmet[i-1] << "\t";
			std::cout << wc[i-1] << "\t";
			std::cout << zc[i-1] << std::endl;

			//sprintf(xlabel[i-1],"%i" ,RunAvg[i].runNumber);
		}

//		TFolder* new_folder = GetHistoFolder("Background_Plots","Background Plots");
		/*BackgroundPlot.RunAvgWMet = new TGraph(size,rn,swmet);
		BackgroundPlot.RunAvgZMet = new TGraph(size,rn,szmet);
		BackgroundPlot.RunAvgWs   = new TGraph(size,rn,wc);
		BackgroundPlot.RunAvgZs   = new TGraph(size,rn,zc);
		BackgroundPlot.RunEvents  = new TGraph(size,rn,evts);
		*/
		//		TAxis *a = BackgroundPlot.RunAvgWMet->GetXaxis();
		//		int bins = a->GetNbins();
		//		int step = bins/size;
		//		std::cout << "bins= " << bins << " step=" << step <<std::endl;

/*		for (int j=0; j < size; j++) {
				a->SetBinLabel(j*step+1,xlabel[j]);	
		}
*/
/*		BackgroundPlot.RunAvgWMet->SetNameTitle("WAvgMetPerRun","Avg. Met per run for W events");
		BackgroundPlot.RunAvgZMet->SetNameTitle("ZAvgMetPerRun"," Avg. Met per run for Z events");
		BackgroundPlot.RunAvgWs->SetNameTitle("WsPerRun"," Number of Ws pee run");
		BackgroundPlot.RunAvgZs->SetNameTitle("ZsPerRun"," Number of Zs per run");
		BackgroundPlot.RunEvents->SetNameTitle("EvtsPerRun"," Events per run");
		new_folder->Add(BackgroundPlot.RunAvgWMet);	
		new_folder->Add(BackgroundPlot.RunAvgZMet);	
		new_folder->Add(BackgroundPlot.RunAvgWs);	
		new_folder->Add(BackgroundPlot.RunAvgZs);	
		new_folder->Add(BackgroundPlot.RunEvents);
*/
	printf("---------------------------------------------------\n");
	return 0;
}
