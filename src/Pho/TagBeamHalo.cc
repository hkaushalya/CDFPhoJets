///////////////////////////////////////////////////////////
// This module will tag the SuperPhotons for Beam Halo,  //
// identified by based on Tower cuts. Several scenarios  //
// to choose. Can run on different modes for now, eg.    //
// Tag only etc. I'll remove these modes later once all  //
// the testing is done. Then this should only be doing   //
// the tagging. nothing else. It will keep myBeamHalo    //
// info in a vector for other modules to be used, if     //
// needed. So no histograms should be here!              //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage  <samantha@fnal.gov>

/*
 *  $Log: TagBeamHalo.cc,v $
 *  Revision 1.14  2009/06/08 03:28:58  samantha
 *  MODIFIED: No tagging is done for MC samples.
 *
 *
 */

#include "samantha/Pho/TagBeamHalo.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnVertex.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "Stntuple/obj/TStnLinkBlock.hh"
#include <iostream>
#include "TMath.h"
#include "TLorentzVector.h"

ClassImp(TagBeamHalo)


//_____________________________________________________________________________
TagBeamHalo::TagBeamHalo(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  iMode(0),
  bRunPermit(true),
  bNoSummary(false)
  
{
	std::cout << "Hello I am TagBeamHalo module" << std::endl;
}

//_____________________________________________________________________________
TagBeamHalo::~TagBeamHalo() {
}

//_____________________________________________________________________________
void TagBeamHalo::SaveHistograms() {
}

//_____________________________________________________________________________
void TagBeamHalo::BookHistograms()
{
	DeleteHistograms();
	BookBeamHaloStudyHistograms(fNoBeamHaloHist, "NoBeamHalo"," None-Beam Halo plots");
	BookBeamHaloStudyHistograms(fWithBeamHaloHist, "WithBeamHalo","Beam Halo plots");
	BookBeamHaloStudyHistograms(fBeamHaloByCutHist, "BeamHaloByCut","Beam Halo identified bu cuts");
	BookBeamHaloStudyHistograms(fBeamHaloByTimeHist, "BeamHaloByTime","Beam Halo identified by em time");
	//BookBeamHaloStudyHistograms(fBHScene0Hist, "Scene0","BH Scenario 0 plots");
	//BookBeamHaloStudyHistograms(fBHScene1Hist, "Scene1","BH Scenario 1 plots");
	//BookBeamHaloStudyHistograms(fBHScene2Hist, "Scene2","BH Scenario 2 plots");
	//BookBeamHaloStudyHistograms(fBHScene3Hist, "Scene3","BH Scenario 3 plots");
	//BookBeamHaloStudyHistograms(fBHScene4Hist, "Scene4","BH Scenario 4 plots");
	//BookBeamHaloStudyHistograms(fBHScene5Hist, "Scene5","BH Scenario 5 plots");

	TFolder* new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		bRunPermit = false;
		return;
	} else {

		Histograms HistoMan;
		vCounterHistLabels.push_back("Events Processed");
		vCounterHistLabels.push_back("Events Passed");
		vCounterHistLabels.push_back("S-0 Halos Found");
		vCounterHistLabels.push_back("S-1 Halos Found");
		vCounterHistLabels.push_back("S-2 Halos Found");
		vCounterHistLabels.push_back("S-3 Halos Found");
		vCounterHistLabels.push_back("S-4 Halos Found");
		vCounterHistLabels.push_back("S-5 Halos Found");
		HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	}
}


//_____________________________________________________________________________
int TagBeamHalo::BeginJob()
{
	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
		bRunPermit = false;
		return 0;
	}

				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fCalBlock 		= (TCalDataBlock*)     RegisterDataBlock("CalDataBlock","TCalDataBlock");
	fPhotonBlock 	= (TStnPhotonBlock*)   RegisterDataBlock("PhotonBlock","TStnPhotonBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");

	if (!fPhotonBlock) {
		StdOut(__FILE__,__LINE__,3,"PhotonBlock not found");
		bRunPermit = false;
		return 0;
	}
	if (!fCalBlock) {
		StdOut(__FILE__,__LINE__,3,"CalDataBlock Block not found");
		bRunPermit = false;
		return 0;
	}
	if (!fMetBlock)	{
		StdOut(__FILE__,__LINE__,3,"MetBlock not found");
		bRunPermit = false;
		return 0;
	}
	if (!fVertexBlock) {
		StdOut(__FILE__,__LINE__,3,"VeretxBlock not found");
		bRunPermit = false;
		return 0;
	}
	
	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.haloCandidatesFoundByCut= 0;

	iHalosFromS0 = 0;
	iHalosFromS1 = 0;
	iHalosFromS2 = 0;
	iHalosFromS3 = 0;
	iHalosFromS4 = 0;
	iHalosFromS5 = 0;

	return 0;
}

//_____________________________________________________________________________
int TagBeamHalo::BeginRun()
{

	//int currRun =  fHeaderBlock->RunNumber();
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TagBeamHalo::Event(int ientry)
{
	SetPassed(0);

	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"Run Permit not cleared!");
		exit (1);
		return 0;
	}
	counter.evtsRunOver++;
	hCount.iCounter->Fill(0);

	if (fHeaderBlock->McFlag())
	{
		SetPassed(1);
	} else {

		Cleanup();
		FillDataBlocks(ientry);
		bool bFoundHalo = false;

		int NsuperPho = initSpMod->GetSuperPhoSize();
		for (int i=0; i< NsuperPho; i++) {
			int iHalo = IsHalo(initSpMod->GetSuperPhoton(i));
			initSpMod->GetSuperPhoton(i)->SetBeamHaloId(iHalo);
			if (iHalo > 0) {
				bFoundHalo = true;
			}
			vHaloStuff.push_back(myHaloStuff);  // the ordering will be determined by the InitSuperPhotons mod.
			// since all mods are driven based on it, no worries!
		}
		if (bFoundHalo) {
			counter.haloCandidatesFoundByCut++;
			hCount.iCounter->Fill(2);
		}

		//lets try to keep track of how many of each we see, although i'll be using only
		//one type of halo. may be i can have this check only for that and reduce computing time!
		// but then i may have to rememebr this when trying make any BH study. study modules
		// looks at all types at one pass. so lets just leave it like this. it may cost a little
		// time. but will be much safer and i don't have remember all these settings! 01-29-2008
		if (iHalosFromS0) hCount.iCounter->Fill(3);
		if (iHalosFromS1) hCount.iCounter->Fill(4);
		if (iHalosFromS2) hCount.iCounter->Fill(5);
		if (iHalosFromS3) hCount.iCounter->Fill(6);
		if (iHalosFromS4) hCount.iCounter->Fill(7);
		if (iHalosFromS5) hCount.iCounter->Fill(8);

		if (GetMode() == 0) {		// by-pass mode
			SetPassed(1);
		} else if (GetMode() == 1) {
			//if (initphomod->GetSuperPhoSize() > 0) {		// check if lead photon not BH
			//	if (!(initphomod->GetSuperPhoton(0)->IsBeamHalo()) ) {
			//		SetPassed(1);
			//	}
			//}
			SetPassed(1);
			std::cout << __FILE__ << "::" << __LINE__ << ":: Mode no longer available. Passing All." << std::endl;
		} else if (GetMode() == 2) {
			SetPassed(1);
			std::cout << __FILE__ << "::" << __LINE__ << ":: Mode no longer available. Passing All." << std::endl;
		} else if (GetMode() == 3) {
			if (initSpMod->GetSuperPhoSize() > 0) {		// check if lead photon is BH
				if (initSpMod->GetSuperPhoton(0)->IsBeamHalo()) {
					SetPassed(1);
				}
			}
		}
	}

	if (GetPassed()) {
		hCount.iCounter->Fill(1);
		counter.evtsPassModule++;
	}

	return 0;
} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TagBeamHalo::Cleanup()
{
	vHaloStuff.clear();
	iHalosFromS0 = 0;
	iHalosFromS1 = 0;
	iHalosFromS2 = 0;
	iHalosFromS3 = 0;
	iHalosFromS4 = 0;
	iHalosFromS5 = 0;
} //Cleanup

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int TagBeamHalo::IsHalo(SuperPhoton* sp)
{
	DoMyBeamHalo(sp->GetPhotonBlockIndex());
	return BeamHaloCut(&myHaloStuff);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TagBeamHalo::DoMyBeamHalo(int pho_ind)
{
	//StdOut("In DoMyBeamHalo",GetName());
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
  
  /*std::cout << _seedWedge	<< "," << 
  					_sideWedge	<< "," <<
  					_seedWedgeH	<< "," <<
  					_sideWedgeH	<< "," <<
				   _eastNHads	<< "," <<
				   _westNHads	<< "," <<
				   _emSeedE		<< "," <<
				   _emSideE		<< "," <<
				   _hadSeedE	<< "," <<
				   _hadSeedWedgeE<< std::endl;
  */
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
void TagBeamHalo::getBeamHaloInfo(int pho_ind, int& seedWedge, int& seedWedgeH,
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
// returns 1 if beam halo. else 0 
/*-------------------------------------------------------------------*/

int TagBeamHalo::BeamHaloCut(HaloStuff *beamhalo)
{
	int iIdWord= 0x0;	//this will include info on all scenarios
									//if beam halo, this will be > 0;
	//scenario == 0
	if ( (beamhalo->seedwedge) > 8 ||
		  ((beamhalo->eastNhad + beamhalo->westNhad) > 2) ) {
		iIdWord |= 0x1 << 0;
		iHalosFromS0++;
	}
  
	//scenario == 1
	if ( (beamhalo->seedwedge) > 4 && 
		  ((beamhalo->eastNhad + beamhalo->westNhad) > 1) ) {
		iIdWord |= 0x1 << 1; 
		iHalosFromS1++;
	}
  
	//scenario == 2
	if ( (beamhalo->seedwedge) > 4 && 
		  ((beamhalo->eastNhad + beamhalo->westNhad) > 2) ) {
		iIdWord |= 0x1 << 2; 
		iHalosFromS2++;
	}
	 
	//scenario == 3
	if ( (beamhalo->seedwedge) > 7 && 
		  ((beamhalo->eastNhad + beamhalo->westNhad) > 2) ) {
		iIdWord |= 0x1 << 3; 
		iHalosFromS3++;
	}

	//sam::added on 10-12-07
	//scenario == 4 
	if ( (beamhalo->seedwedge) > 8 && 
		  ((beamhalo->eastNhad + beamhalo->westNhad) > 2) ) {
		iIdWord |= 0x1 << 4; 
		iHalosFromS4++;
	}
		
	//scenario == 5
	if ( (beamhalo->seedwedge) > 8 && 
		  ((beamhalo->eastNhad + beamhalo->westNhad) > 3) ) {
		iIdWord |= 0x1 << 5; 
		iHalosFromS5++;
	}
			  
	return iIdWord;
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

int TagBeamHalo::Class12Vertices()
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
void TagBeamHalo::MatchCalorTowers(int pho_ind, CalDataArray* cdholder)
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
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TagBeamHalo::BookBeamHaloStudyHistograms(BeamHaloStudyHisto_t& Hist,
									std::string sFoldName, std::string sFoldTitle)
{
	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}

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
  sprintf(title," %s%s",sFoldName.c_str(),": CES(strip+wire) #chi^{2} of halo candidate");
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
	
  sprintf(name,"Halo_MEt");
  sprintf(title,"ME_{T} of Halo Events");
  Hist.fHaloMet = new TH1F(name,title,1500,0,1500);
  new_folder->Add(Hist.fHaloMet);

  sprintf(name,"Halo_SumEt");
  sprintf(title,"Sum E_{T} of Halo Events");
  Hist.fHaloSumet = new TH1F(name,title,1500,0,1500);
  new_folder->Add(Hist.fHaloSumet);
  
  return;
}


/*-------------------------------------------------------------------*/
//_____________________________ Filling beam halo histograms 
void  TagBeamHalo::FillBeamHaloStudyHistograms(BeamHaloStudyHisto_t& Hist, 
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

	if (DelPhi >TMath::Pi()) {
		DelPhi -= TMath::TwoPi();
	}
	if (DelPhi < (-1.* TMath::Pi())) {
		DelPhi += TMath::TwoPi();
	}

	Hist.fHalo_Met_delPhi->Fill(fabs(DelPhi));
	Hist.fSeedWedge_HadTowers->Fill(beamhalo->eastNhad+beamhalo->westNhad,
			beamhalo->seedwedge);   
  return;
}



/*-------------------------------------------------------------------*/
void TagBeamHalo::FillDataBlocks(int ientry)
{
	fPhotonBlock->GetEntry(ientry);
	fCalBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);	
	fVertexBlock->GetEntry(ientry);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TagBeamHalo::EndJob() {

	if (GetSummaryStat()) return 0;

	std::string sMsg, sMsg2, sMsg3;
	if (GetMode() == 0) sMsg2 = "Tag Only and Tag all photons.";
	if (GetMode() == 1) sMsg2 = "Filter Mode 1 (Tag all and pass only if lead Photon is NOT BH)";
	if (GetMode() == 2) sMsg2 = "Filter Mode 2 (Tag all and pass if a BH is found)";
	if (GetMode() == 3) sMsg2 = "Filter Mode 3 (pass only if lead Photon is BH)";


	std::string sModName = GetName();
	sMsg += "[TBH:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[TBH:01:] Events Processed ----------- = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	sMsg += "[TBH:02:] Events Passed -------------- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	if (fHeaderBlock->McFlag())
	{
		sMsg += "[TBH:03:] NO TAGGING DONE FOR MC **** \n";
	} else {
		sMsg += "[TBH:03:] Total Halo Events Found  --- = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
		sMsg += "[TBH:04:] Halos found by S0 ---------- = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n";
		sMsg += "[TBH:05:] Halos found by S1 ---------- = " + ToStr(hCount.iCounter->GetBinContent(5)) + "\n";
		sMsg += "[TBH:06:] Halos found by S2 ---------- = " + ToStr(hCount.iCounter->GetBinContent(6)) + "\n";
		sMsg += "[TBH:07:] Halos found by S3 ---------- = " + ToStr(hCount.iCounter->GetBinContent(7)) + "\n";
		sMsg += "[TBH:08:] Halos found by S4 ---------- = " + ToStr(hCount.iCounter->GetBinContent(8)) + "\n";
		sMsg += "[TBH:09:] Halos found by S5 ---------- = " + ToStr(hCount.iCounter->GetBinContent(9)) + "\n";
	}

	sMsg += "---------------------------------------------------\n";
	std::cout << sMsg;
	return 0;
}


//---------------------------------------------------------
// FREE FUNCTIONS
//---------------------------------------------------------
void FillHaloHistograms(Histograms::BeamHaloHists_t& Hist, TagBeamHalo::HaloStuff beamhalo,
								TStnPhoton* Pho)
{
	// common method to fill all beam halo hists

	if (!Pho) {
		StdOut(__FILE__,__LINE__,3,"Invalid pointer to a TStnPhoton. Returning without doing anything!.");
		return;
	}

	Hist.fHaloSeedWedgeEM_SeedWedgeHAD->Fill(beamhalo.seedwedgeH,beamhalo.seedwedge);
	 //______________________ Max's params
	Hist.fHaloSeedWedgeH->Fill(beamhalo.seedwedgeH);   
	Hist.fHaloSideWedgeH->Fill(beamhalo.sidewedgeH); 
	Hist.fHaloSeedWedge->Fill(beamhalo.seedwedge);   
	Hist.fHaloSideWedge->Fill(beamhalo.sidewedge); 
	Hist.fHaloEastNHad->Fill(beamhalo.eastNhad);
	Hist.fHaloWestNHad->Fill(beamhalo.westNhad);
	Hist.fHaloEmSeedE->Fill(beamhalo.emseedE);  
	Hist.fHaloEmSideE->Fill(beamhalo.emsideE);   
	Hist.fHaloHadSeedE->Fill(beamhalo.hadseedE); 
	Hist.fHaloSeedWedgeHadE->Fill(beamhalo.seedwedgeHadE);
	Hist.fHaloNHad->Fill(beamhalo.eastNhad+beamhalo.westNhad);  
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
	Hist.fSeedWedge_HadTowers->Fill(beamhalo.eastNhad+beamhalo.westNhad,
			beamhalo.seedwedge);   

}

