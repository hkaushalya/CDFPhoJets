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
#include "samantha/PhoJetsMet/TWZcrossSection.hh"
#include "Stntuple/photon/TPhotonUtil.hh"


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
const static double	kL_HADEM	= 0.125;
const static double	kL_ISOETCORR_MF_LTETC20 = 0.15;	//ISOETCORR Multiplicative factor when Etc<20
const static double	kL_ISOETCORR_MF_GTETC20 = 3.0;	//when Etc>20
const static double	kL_TRKPT_MF	= 0.25;		// TRKPT Multiplicative factor (*EtCorr)
const static double	kL_TRKISO	= 5.0;

//tight photon id cut values

const static int	 	kT_CENTRAL 	= 0; 		// CENTRAL
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
bool PrePass()
bool PassGoodRun()
bool PassTrigger()
bool PassZvCut()
int Class12Vertices()
double BestVetex_Z()

void WZ_dataEstimate()
void WEleIDcutEff();
void ZEleIDcutEff();

void TagCoversions()
int FiducialLepton()
SetPhotonList

CEMPhoEleIDcut(pho,event, zvx)
TightPhotonIDcut(pho)
LoosePhotonIDcut(pho)

FillDataBlock(ientry)


FillMetPlots

TFolder* GetHistoFolder
BookGeneralHistograms
BookMetHistograms

*/



ClassImp(TWZcrossSection)

//_____________________________________________________________________________
TWZcrossSection::TWZcrossSection(const char* name, const char* title):
  TStnModule(name,title), qMc(false)
{
	std::cout << "Hello I am TWZcrossSection module" << std::endl;
}

//_____________________________________________________________________________
TWZcrossSection::~TWZcrossSection() {
}

//_____________________________________________________________________________
void TWZcrossSection::SaveHistograms() {
}

//_____________________________________________________________________________
void TWZcrossSection::BookHistograms() {

	DeleteHistograms();
	BookGeneralHistograms();
	BookMetHistograms();
	BookBackgroundHistograms();
	BookWeleIDeffHistograms();
	BookZeleIDeffHistograms();
}


//_____________________________________________________________________________
int TWZcrossSection::BeginJob()
/*{{{*/
{
	std::cout << " BEGIN JOB " << std::endl;
				// register the data block, why?  to unpack them!

	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
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
	counter.evtsPassVertexCut = 0;
	counter.evtsPassZvCut = 0;

	counter.evtsPassModule = 0;
	counter.evtsRunOver = 0;

	counter.ideff_hepgWs = 0;				//__ hepg W events
	counter.ideff_hepgWenus = 0;			//__ W that decayed to electron+ nutrino
	counter.ideff_passGoodrun = 0;		//__ from the events with a W, number of events pass good run
	counter.ideff_passVertex = 0;		//__ then pass 1 good vertex
	counter.ideff_passZv = 0;				//__ then pass Zv cut
	counter.ideff_passMetcut = 0;
	counter.ideff_detWs = 0;


	counter.zideff_hepgZs = 0;				//__ hepg Z events
	counter.zideff_hepgZees = 0;			//__ Z that decayed to electron+ nutrino
	counter.zideff_passGoodrun = 0;		//__ from the events with a Z, number of events pass good run
	counter.zideff_passVertex = 0;		//__ then pass 1 good vertex
	counter.zideff_passZv = 0;
	counter.zideff_passZv = 0;				//__ then pass Zv cut
	counter.zideff_detZs = 0;



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
int TWZcrossSection::BeginRun()
{

	currRun = fHeaderBlock->RunNumber();
	std::cout << " BEGIN RUN " << currRun << std::endl;
	
	//store previous run's summary
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
int TWZcrossSection::Event(int ientry)
{
	SetPassed(0);
	counter.evtsRunOver++;

	//_____________________________ REQUIRED Thingies 
  FillDataBlocks(ientry);
  
	if (PrePass()) {
		SetPassed(1);
		counter.evtsPassModule++;
		Cleanup();
		SetPhotonList();
		TagConversions();
		//WZ_dataEstimate();			// to get W and Z candidates from photon DATA sample
		//if (qMc) WEleIDcutEff();	// to estimate the electron id eff. for W cross section calculation
		if (qMc) ZEleIDcutEff();	// to estimate the electron id eff. for Z cross section calculation
	}

	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TWZcrossSection::Cleanup()
{
  	trig25 = false;
	trig50 = false;
	trig70 = false;

	superpho.clear();


} //Cleanup

/*-------------------------------------------------------------------*/
// Event pre-selection cuts
// REMOVE EVENTS WITHOUT PHOTONS
// do not cut on number of jets yet. do it in the TGammaJets mod
/*-------------------------------------------------------------------*/
bool TWZcrossSection::PrePass()
{
	if (!qMc) {
		if (fPhotonBlock->NPhotons() < 1) return false;
		else return (PassGoodRun() && PassTrigger() && PassZvCut());
	} else {
		//return (PassGoodRun() && PassZvCut());  // for MC runs, this is done within that module. so ignore this
		return true;
	}
} // PrePass

/*-------------------------------------------------------------------*/
// good run preselection
// remove the bad runs/sections
/*-------------------------------------------------------------------*/
bool TWZcrossSection::PassGoodRun()
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
bool TWZcrossSection::PassTrigger()
{
	if (!qMc) { //___________________ pypho22_dfc MC sample does not have trigger info
		trigbits.Event();
    	trig25 = trigbits.Pho25Iso();
	   trig50 = trigbits.Pho50();
   	trig70 = trigbits.Pho70();

	    if ( ! (trig25 || trig50 || trig70) ) {
   		return false; 
	    } else {		//_____________ use all 3 triggers. may see something at hight pt.
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
/*-------------------------------------------------------------------*/
bool TWZcrossSection::PassZvCut()
{
	int nvtx = Class12Vertices();

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
int TWZcrossSection::Class12Vertices()
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
double TWZcrossSection::BestVertex_Z()
{
  double zvx =0;
  if (fVertexBlock->GetBestVertex(12,1) != NULL) {
    zvx = fVertexBlock->GetBestVertex(12,1)->Z();
  }
  
  return zvx;
} //BestVertex_Z

/*-------------------------------------------------------------------*/
// W/Z background estimate from data
// 1. remover conversions
// 2. find photons passing EleLikeIdCUTS
// 3. calculate W/Zs
/*-------------------------------------------------------------------*/
void TWZcrossSection::WZ_dataEstimate()
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

/*-------------------------------------------------------------------*/
void TWZcrossSection::FillWplots(int index)
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


/*-------------------------------------------------------------------*/
void TWZcrossSection::FillZplots(int index1, int index2)
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
1. find hepg Z
2. apply Z window cut to remove Drell-Yan
3. find its decay e's, keep their vectors
DETECTOR LEVEL
4. apply goodrun ,vertex,Zv cut
5. look for 2 photons passing ID CUTs
*/
//--------------------------------------------------------------	
void TWZcrossSection::ZEleIDcutEff()
{
	int Nparticles = fGenpBlock->NParticles();
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
					TLorentzVector zvec;
					par->Momentum(zvec);
					ZEleIDeffPlot.hepgZpt->Fill(zvec.Pt());
					ZEleIDeffPlot.hepgZmass->Fill(zvec.M());
					
					if (zvec.M() < 66. || zvec.M() > 116.) return;  // must apply this to get the correct eff. 
					
					counter.zideff_hepgZs++;
					int p1 = par->GetFirstDaughter();
					int pl = par->GetLastDaughter();
					
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
					
					break;	//W is found, but did not decay in to enu, so quit
				} 
			} //if W
		} // if
	} //for
	
	if (ndau == 2) {
		counter.zideff_hepgZees++;
	} else return;
	

	if (PassGoodRun()) {
		counter.zideff_passGoodrun++;
	} else return;
	
	if (Class12Vertices() >= 1) {
		counter.zideff_passVertex++;
	} else return;

	if (fabs(BestVertex_Z()) < 60) {
		counter.zideff_passZv++;
	} else return;

	int ne = 0;
	for (int i=0; i < superpho.size(); i++ ) {
  		if ( (superpho[i].Conversion == 0) && (superpho[i].PhotonLikeID == 0) ) {
			ne++;
		}
	}

	if (ne ==2) counter.zideff_detZs++;


} //ZEleIDcutEff


//--------------------------------------------------------------	
// measure eff*accpt from CEM PHOTON LIKE ELE ID CUTS for Wenu
// follow steps as in ZEleIDcutEff, look for W->enu in hepg level
//  and one reconstructed electron passing electron id cuts witn Met cut
//--------------------------------------------------------------	
void TWZcrossSection::WEleIDcutEff()
{
	int Nparticles = fGenpBlock->NParticles();

	int mom_id = 0;
	TLorentzVector helevec(0,0,0,0);		// hepg stable electron
	double heleEt = 0;						// hepg electron et
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
					counter.ideff_hepgWs++;
					int p1 = par->GetFirstDaughter();
					int pl = par->GetLastDaughter();

					for (int i= p1 ; i <= pl; i++) {
						dau = fGenpBlock->Particle(i);
						int dau_id = dau->GetPdgCode();
						int dau_stat = dau->GetStatusCode();
						if ( (abs(dau_id) == 11) && (dau_stat == 1)) {
							counter.ideff_hepgWenus++;
							dau->Momentum(helevec);
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
	} else return;
	
	if (Class12Vertices() >= 1) {
		counter.ideff_passVertex++;
	} else return;

	if (fabs(BestVertex_Z()) < 60) {
		counter.ideff_passZv++;
	} else return;

	if (fMetBlock->Met(0) > 20) {
		counter.ideff_passMetcut++;
	} else return;
	
	int ne = 0;	
	for (int i=0; i< superpho.size(); i++) {	
  		if ((superpho[i].Conversion == 0) && (superpho[i].PhotonLikeID == 0)) {
			ne++;
		} 
	} 

	if (ne == 1) counter.ideff_detWs++;

} //WEleIDcutEff


/*-------------------------------------------------------------------*/
//tag conversion electrons for W calculations
/*-------------------------------------------------------------------*/
void TWZcrossSection::TagConversions()
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
bool TWZcrossSection::IsConversionElectron(TStnElectron* ele, 
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
double TWZcrossSection::ConversionSep(TStnTrack* t0, TStnTrack* t1, double & rconv)
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



/*-------------------------------------------------------------------*/
// check if the photon is fidual before matching to HEPG Level for W estimate.
/*-------------------------------------------------------------------*/
int TWZcrossSection::FiducialLepton(TLorentzVector vec_old, double z_new)
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


/*-------------------------------------------------------------------*/
//	Sorts the Photon Block according to EtCorr, from highest to lowest
// & fill the photon info to my SuperPhoton list
/*-------------------------------------------------------------------*/
void TWZcrossSection::SetPhotonList()
{
	int Npho = fPhotonBlock->NPhotons();
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

		// fill plot of electron id cumulative
		int fillto = sp.PhotonLikeID;
		if (sp.PhotonLikeID ==0) fillto = 17;
		for (int j=0; j < fillto; j++) GeneralPlot.CEMPHOEleIDeff->Fill(j);

		superpho.push_back(sp);
	}

} //SetPhotonList



/*****************  PHOTON LIKE ID CUT FUNCTION ************************/
//_____ function returns 0 if electrons passes cuts
//      requires zvx of highest Pt class12 vertex
/*---------------------------------------------------------------------*/
unsigned int TWZcrossSection::CEMPhoEleIDcut(TStnPhoton* Pho, TStnEvent* event,
															double zvx)
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
unsigned int TWZcrossSection::TightPhotonIDcut(TStnPhoton *pho)
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
unsigned int TWZcrossSection::LoosePhotonIDcut(TStnPhoton *pho)
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
void TWZcrossSection::FillDataBlocks(int ientry)
{
	fTriggerBlock->GetEntry(ientry);
	fPhotonBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
	fVertexBlock->GetEntry(ientry);
	fTrackBlock->GetEntry(ientry);
	fGenpBlock->GetEntry(ientry);
}

/*-------------------------------------------------------------------*/
// fill  Missing Et plots
/*-------------------------------------------------------------------*/
void TWZcrossSection::FillMetPlots()
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
TFolder* TWZcrossSection::GetHistoFolder(char *name, char* title)
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
void TWZcrossSection::BookGeneralHistograms()
{
	std::cout << "Booking General Plots...";
	
	GeneralPlot.Nvertices = new TH1F("Nvertices","Number of Class12 vertices per event",10,0,10);
	GeneralPlot.vertexZ 	 = new TH1F("BestVertex_Z","Best Class12 vertex Z position",100,-100,100);
	GeneralPlot.CEMPHOEleIDeff 	= new TH1F("EleIDeff","Electron ID eff for W",18,0,18);
	TAxis* x2= GeneralPlot.CEMPHOEleIDeff->GetXaxis();
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(1,"PhoCentral");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(2,"PhoEtc30");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(3,"PhoHadEm");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(4,"PhoChi2");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(5,"PhoN3d=0");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(6,"PhoN3d>2");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(7,"PhoPt2");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(8,"FindMatchEle");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(9,"EleHasTrk");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(10,"EleTrkZ");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(11,"EleEoveP");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(12,"PhoCalIso4");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(13,"PhoTrkIso");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(14,"PhoCes2Et");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(15,"PhoXCes");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(16,"PhoZCes");
	GeneralPlot.CEMPHOEleIDeff->GetXaxis()->SetBinLabel(17,"Electrons");

	TFolder* new_folder = GetHistoFolder("General_Plots","General Plots");
	new_folder->Add(GeneralPlot.Nvertices);
	new_folder->Add(GeneralPlot.vertexZ);
	std::cout << "DONE\n";

} //BookGeneralHistograms

//_____________________________________________________________________________
void TWZcrossSection::BookMetHistograms()
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
void TWZcrossSection::BookBackgroundHistograms()
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
void TWZcrossSection::BookWeleIDeffHistograms()
{
	std::cout << "Booking W ele ID Eff Plots...";
	WEleIDeffPlot.MinDelRMatch	= new TH1F("Wele2CloseOffPhoDelR","DelR between W Electron and Closest offline Photon",500,0,5);
	WEleIDeffPlot.HepgEleEt		= new TH1F("HepgEleEt","W electron Et",120,0,60);
	WEleIDeffPlot.Et_ratio		= new TH1F("Et_Offele2Hepgele","Matched Offline Electron to HEPG Object Et ratio",40,0,2);
	WEleIDeffPlot.Et_ratio->SetYTitle("Events/0.05");

	TFolder* new_folder = GetHistoFolder("WeleIDeff_Plots","Electron ID eff Plots");
	new_folder->Add(WEleIDeffPlot.HepgEleEt);
	new_folder->Add(WEleIDeffPlot.MinDelRMatch);
	new_folder->Add(WEleIDeffPlot.Et_ratio);
	std::cout << "DONE.\n";
}


//_____________________________________________________________________________
void TWZcrossSection::BookZeleIDeffHistograms()
{
	std::cout << "Booking Z ele ID Eff Plots...";
	ZEleIDeffPlot.hepgZpt  		= new TH1F("hepgZpt","HEPG Z pt",15,0,15);
	ZEleIDeffPlot.hepgZmass  	= new TH1F("hepgZmass","HEPG Z mass",60,0,120);

	TFolder* new_folder = GetHistoFolder("ZeleIDeff_Plots","Z Electrons ID eff Plots");
	new_folder->Add(ZEleIDeffPlot.hepgZpt);
	new_folder->Add(ZEleIDeffPlot.hepgZmass);
	std::cout << "DONE.\n";
}
//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TWZcrossSection::EndJob() {

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
		std::cout << " ======================= W Electron ID eff =========================" << std::endl;
		std::cout << "Evts with a HEPG W ------------------ = " << counter.ideff_hepgWs << std::endl;
		std::cout << "HEPG W->enu ------------------------- = " << counter.ideff_hepgWenus << std::endl;
		std::cout << "Det. evt passing goodrun ------------ = " << counter.ideff_passGoodrun << std::endl;
		std::cout << "passing good vertex ----------------- = " << counter.ideff_passVertex << std::endl;
		std::cout << "passing Zv < 60  -------------------- = " << counter.ideff_passZv << std::endl;
		std::cout << "passing Met > 20 -------------------- = " << counter.ideff_passMetcut << std::endl;
		std::cout << "Evts with Det W (passing Ele ID cuts) = " << counter.ideff_detWs << std::endl;
		std::cout << " ======================= Z Electron ID eff =========================" << std::endl;
		std::cout << "Evts with a HEPG Z ------------------ = " << counter.zideff_hepgZs << std::endl;
		std::cout << "HEPG Z->ee -------------------------- = " << counter.zideff_hepgZees << std::endl;
		std::cout << "Det. evt passing goodrun ------------ = " << counter.zideff_passGoodrun << std::endl;
		std::cout << "passing good vertex ----------------- = " << counter.zideff_passVertex << std::endl;
		std::cout << "passing Zv < 60  -------------------- = " << counter.zideff_passZv << std::endl;
		std::cout << "Evts with Det Z (both e pass EleIDcu) = " << counter.zideff_detZs << std::endl;
		std::cout << " ===================================================================" << std::endl;
	}


		// add stuff from the final run in the job
		RunAvg_t ra;
		ra.runNumber = currRun;
		ra.events = counter.evtsRunOver;
		ra.sumWMet = sumWeventMet; 
		ra.sumZMet = sumZeventMet; 
		ra.Wcount = counter.dWe_detWs;
		ra.Zcount = counter.dZe_detZs;
		RunAvg.push_back(ra);
	
		const Int_t size = RunAvg.size()-1;		//must be 'const'
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

		std::cout <<"   Run#   Events  avgWmet  avgZmet  Wcount   Zcount" << std::endl; 

		for (Int_t i=0; i < RunAvg.size(); ++i) {
			if (i==0) continue;
			rn[i-1] 	= RunAvg[i].runNumber;
			evts[i-1]= RunAvg[i].events - RunAvg[i-1].events;
			wc[i-1] 	= RunAvg[i].Wcount - RunAvg[i-1].Wcount;;
			zc[i-1] 	= RunAvg[i].Zcount - RunAvg[i-1].Zcount;;

			if (wc[i-1] != 0) swmet[i-1] = RunAvg[i].sumWMet / wc[i-1];
			if (zc[i-1] != 0) szmet[i-1] = RunAvg[i].sumZMet /  zc[i-1];
			
			std::printf("%7.0f  %7.0f  %6.2f  %6.2f  %7.0f  %7.0f\n",rn[i-1],evts[i-1], swmet[i-1],szmet[i-1],wc[i-1],zc[i-1]); 
		}

	printf("---------------------------------------------------\n");
	return 0;
}
