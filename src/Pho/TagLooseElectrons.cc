///////////////////////////////////////////////////////////
// this module will tag the loose electrons in the super //
// photon list, this will be usefull when calculating    //
// errors for the fake where we need to pick a loose     //
// electron sample. only diff from TagElectrons is loose //
// E/P cut. 02-07-2008                                   //
// I am adding std. loose electron ID cuts. Need to ID   //
// all possible em obj. for met model to work properly.  //
// intead of creating two modules for central and plug,  //
// I am doing all loose electron tagging here.           //
// I decided to keep the "detector" cut in the ID word.  //
// I was thinking of dropping it from the ID word and    //
// check for "detector" region as they are needed. But   //
// there are many modules that needs this ordering.      //
// Besides there is someone already using my code. So    //
// removing things is not a good idea.  08-23-2008       //
// Now I had to change the Et threshold for MET model.   //
// I need to remove all em objects with Et>7GeV. So all  //
// tagging will be done with 7GeV cut. The final photon  //
// selection need to do an additional check for Et>30GeV //
// photon. 09-25-2008                                    //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


/* $Id: TagLooseElectrons.cc,v 1.15 2011/05/25 20:26:29 samantha Exp $
 * $Log: TagLooseElectrons.cc,v $
 * Revision 1.15  2011/05/25 20:26:29  samantha
 * Updated the two Nvtx accessors in the Triggger Module.
 *
 * Revision 1.14  2009/11/10 19:39:47  samantha
 * ADDED:1.Counter for different types of electrons (PhoLike, std etc). 2. Auto CVS Log
 * info.
 *
 *
 */


#include <iostream>
#include <fstream>
#include <TMath.h>
#include <numeric>
#include <vector>
#include <algorithm>
#include "Stntuple/loop/TStnAna.hh"
#include "samantha/Pho/TagLooseElectrons.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/obj/TStnTrack.hh"

ClassImp(TagLooseElectrons)

//_____________________________________________________________________________
TagLooseElectrons::TagLooseElectrons(const char* name, const char* title):
  TStnModule(name,title),
  	iMode(0),
	bRunPermit(true),
  	bNoSummary(false)
{
	//std::cout << "Hello I am TagLooseElectrons module" << std::endl;

	float EtMin = 7.0;
	float EtMax = 1500.0;

//plug loose electron cuts
	fPemMinEt = EtMin;
	fPemMaxEt = EtMax;
	//fPemMinEtaPes=1.2; 
	fPemMinEtaPes=1.5; 
	//fPemMaxEtaPes=2.0;
	fPemMaxEtaPes=3.0;
	fPemMinHADEM=0.0; 
	fPemMaxHADEM=0.05;       
	fPem3x3FitCut=0; // cut on PEM 3x3 tower fit status code
	fPemMin3x3Chi2=-10.0E6;  
	fPemMax3x3Chi2=10.0;
	fPemMinPes5x9U=0.65; 
	fPemMaxPes5x9U=10.0E6; 
	fPemMinPes5x9V=0.65; 
	fPemMaxPes5x9V=10.0E6;
	fPemMinIsoFr=-1.0;  
	fPemMaxIsoFr=0.1;
	fPemMinPesdR=0.0; 
	fPemMaxPesdR=3.0;

	// central loose electron cuts
  	fCemMinEt= EtMin;    
  	fCemMaxEt= EtMax; 
	fCemMinTrkZ0=0.0; 
	fCemMaxTrkZ0=60.0;
	fTrkMinPt=3.5;
	//fTrkMinPt=10.0;		// do i need to modify this for my analysis.  I could use looser cut to ge tmore electrons identified.
  								// only thing to make sure is that, I don't stray too far away from sasha's cuts. that way
								// i can always compare my met model results quikly to his. 08-24-2008
	fTrkMaxPt=1000.0; 
	fCemMinCotAxSl=3; // cut on number of COT Axial SL with 5 or more hits 
	fCemMaxCotAxSl=8;
	fCemMinCotStSl=2; // cut on number of COT Stereo SL with 5 or more hits 
	fCemMaxCotStSl=8;
	fCemMinHADEM[0]=0.0; 
	fCemMinHADEM[1]=0.0; 
	fCemMaxHADEM[0]=0.055;                      
	fCemMaxHADEM[1]=0.00045;                      
	fCemMinIsoFr=-1.0;  
	fCemMaxIsoFr=0.1;  
 	fCemFidCut=1; // cut on CEM fiduciality. Sasha does this differently???
}

//_____________________________________________________________________________
TagLooseElectrons::~TagLooseElectrons() {
}

//_____________________________________________________________________________
void TagLooseElectrons::SaveHistograms() {
}

//_____________________________________________________________________________
void TagLooseElectrons::BookHistograms()
{
	DeleteHistograms();
	TFolder* new_folder1 = GetHistFolder(this, "B4Cuts","plots of all id varibles b4 cuts");
	if (!new_folder1) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	TFolder* new_folder2 = GetHistFolder(this, "A4Cuts","plots of all id varibles a4 cuts");
	if (!new_folder2) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	Hist.GetEmObjHistograms(new_folder1, EleB4Cut,"Before cuts");
	Hist.GetEmObjHistograms(new_folder2, EleA4Cut,"After cuts");

	TFolder* new_folder3 = GetHistFolder(this, "StdIdCuts","plots of all std. electron id varibles");
	//int nid =  iNStdLooseCenEleIdBits+1; 
	idHist.StdCen = new TH1F("stdCen","Standard Central Loose Electron ID cuts",iNStdLooseCenEleIdBits,0,iNStdLooseCenEleIdBits);
	idHist.StdPlug = new TH1F("stdPlug","Standard Plug Loose Electron ID cuts",iNStdLoosePlugEleIdBits,0,iNStdLoosePlugEleIdBits);

	new_folder3->Add(idHist.StdCen);
	new_folder3->Add(idHist.StdPlug);
	
	idHist.StdPlug->GetXaxis()->SetBinLabel(1,"Central");
	idHist.StdPlug->GetXaxis()->SetBinLabel(2,"Et>7 ");
	idHist.StdPlug->GetXaxis()->SetBinLabel(3,"PesEta");
	idHist.StdPlug->GetXaxis()->SetBinLabel(4,"HadEm");
	idHist.StdPlug->GetXaxis()->SetBinLabel(5,"Pem3x3");
	idHist.StdPlug->GetXaxis()->SetBinLabel(6,"Chi2Three");
	idHist.StdPlug->GetXaxis()->SetBinLabel(7,"Pes5x9U");
	idHist.StdPlug->GetXaxis()->SetBinLabel(8,"Pes5x9V");
	idHist.StdPlug->GetXaxis()->SetBinLabel(9,"PemIsoFr");
	idHist.StdPlug->GetXaxis()->SetBinLabel(10,"PesDelR");

	idHist.StdCen->GetXaxis()->SetBinLabel(1,"Plug");
	idHist.StdCen->GetXaxis()->SetBinLabel(2,"Et>7 ");
	idHist.StdCen->GetXaxis()->SetBinLabel(3,"Fid");
	idHist.StdCen->GetXaxis()->SetBinLabel(4,"Trk Num>0");
	idHist.StdCen->GetXaxis()->SetBinLabel(5,"TrkZ0");
	idHist.StdCen->GetXaxis()->SetBinLabel(6,"TrkPt");
	idHist.StdCen->GetXaxis()->SetBinLabel(7,"Nasl");
	idHist.StdCen->GetXaxis()->SetBinLabel(8,"Nssl");
	idHist.StdCen->GetXaxis()->SetBinLabel(9,"HadEm");
	idHist.StdCen->GetXaxis()->SetBinLabel(10,"IsoFr");
	
}

//_____________________________________________________________________________
int TagLooseElectrons::BeginJob()
{
	initSp = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initSp == NULL) {
		StdOut(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
		bRunPermit = false;
	}

  	trigMod = (TriggerModule*) ((TStnAna*) GetAna()->GetModule("TriggerModule"));
	if (!trigMod) {
		StdOut(__FILE__,__LINE__,3,"TriggerModule required!.");
		bRunPermit = false;
	}


	
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");
	fTrackBlock 	= (TStnTrackBlock*)    RegisterDataBlock("TrackBlock","TStnTrackBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 		= 0;
	counter.evtsPassModule 	= 0;
	counter.electronEvents	= 0;
	counter.PhoLikeEleEvts	= 0;
	counter.PhoLikeEles	= 0;
	counter.StdEleEvts	= 0;
	counter.StdCEMEles	= 0;
	counter.StdPLUGEles	= 0;

	return 0;
}

//_____________________________________________________________________________
int TagLooseElectrons::BeginRun()
{

 return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TagLooseElectrons::Event(int ientry)
{
	SetPassed(1);
	counter.evtsRunOver++;
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"One or more dependencies not found. pl check. Run Permit not cleared!");
		return 0;
	}
	
	
	fElectronBlock->GetEntry(ientry);
	fTrackBlock->GetEntry(ientry);
	fVertexBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);
	
	double zvx  = 0;
	if (fVertexBlock->GetBestVertex(12,1) != NULL) {
			zvx = fVertexBlock->GetBestVertex(12,1)->Z();
	}

	TStnEvent* event = GetEvent();
	int iNElectrons = 0;
	int iPhoLikeEles = 0, iStdCEMEles = 0, iStdPLUGEles = 0;
	int Nsp = initSp->GetSuperPhoSize();

	for (int i = 0; i < Nsp; i++) {
		//FillHistograms(EleB4Cut, initSp->GetSuperPhoton(i));  // lets make the plots b4 cuts. anyone who wants
																								// after cuts can do that in that module
		int eid = GetPhoLikeEleIdWord(initSp->GetSuperPhoton(i)->GetPhoton(),event, fTrackBlock, zvx);
		if (eid == 0) {
			//FillHistograms(EleA4Cut, initSp->GetSuperPhoton(i));		//need to redo these. either use SP or have struct
			counter.PhoLikeEles++;
			iPhoLikeEles++;
		}
		initSp->GetSuperPhoton(i)->SetLooseElectronId(eid);
		

		// NOW TAG LOOSE ELECTRON USING STANDARD ELECTRON ID CUTS
		int iStdLooseEleId = -1;

		TStnElectron* Ele = TPhotonUtil::MatchElectron(event,initSp->GetSuperPhoton(i)->GetPhoton()); // get matching electron
		
		if (initSp->GetSuperPhoton(i)->GetDetector() == 0 )		//central
		{
			iStdLooseEleId = GetStdCenLooseEleID(Ele, fTrackBlock, trigMod->GetN12vertex());
			initSp->GetSuperPhoton(i)->SetStdLooseElectronId(iStdLooseEleId);
			if (iStdLooseEleId == 0)
			{
				counter.StdCEMEles++;
				iStdCEMEles++;
			}
		
		} else if (initSp->GetSuperPhoton(i)->GetDetector() == 1 )		//plug
		{
			iStdLooseEleId = GetStdPlugLooseEleID(Ele,event, trigMod->GetN12vertex());
			initSp->GetSuperPhoton(i)->SetStdLooseElectronId(iStdLooseEleId);
			if (iStdLooseEleId == 0)
			{
				counter.StdPLUGEles++;
				iStdPLUGEles++;
			}
		}

		//FillStdElectronIDcutPlots(initSp->GetSuperPhoton(i));

		if (eid == 0 || iStdLooseEleId == 0)
			iNElectrons++;
		
	}
	
	if (iNElectrons >0) counter.electronEvents++;
	if (iPhoLikeEles>0) counter.PhoLikeEleEvts++;
	if (iStdCEMEles>0 || iStdPLUGEles>0) counter.StdEleEvts++;

	if (GetMode() == 0) SetPassed(1);
	else if (GetMode() == 1  && iNElectrons == 1)	SetPassed(1);
	else if (GetMode() == 2  && iNElectrons == 2)	SetPassed(1);

	if (GetPassed()) counter.evtsPassModule++;

	return 0;

} // Event

/*-------------------------------------------------------------------*/
void TagLooseElectrons::FillStdElectronIDcutPlots(SuperPhoton* sp)
{
	int stdLooseEleId = sp->GetStdLooseElectronId();
	
	if (sp->GetDetector() == 0)	//central
	{
		for (int i=0; i < iNStdLooseCenEleIdBits; ++i)
		{
			if ((stdLooseEleId >> i) && 0x1) break;
			else	idHist.StdCen->Fill(i);
		}
	}
	if (sp->GetDetector() == 1)	//plug
	{
		for (int i=0; i < iNStdLoosePlugEleIdBits; ++i)
		{
			if (!((stdLooseEleId >> i) && 0x1))
			{
				idHist.StdPlug->Fill(i);
			} else 
				break;
		}
	}
}

/*-------------------------------------------------------------------*/
void TagLooseElectrons::FillHistograms(Histograms::EmObjHists_t& hist, SuperPhoton* sp)
{
	hist.Detector->Fill(sp->GetDetector());
	hist.DetEta->Fill(sp->GetDetEta());
	hist.DetPhi->Fill(sp->GetDetPhi());
	hist.EtCorr->Fill(sp->GetEtCorr());
	hist.HadEm->Fill(sp->GetHadEm());
	hist.Chi2Mean->Fill(sp->GetChi2Mean());
	hist.N3d->Fill(sp->GetN3d());
	TStnEvent* event = GetEvent();
	TStnElectron* Ele=TPhotonUtil::MatchElectron(event,sp->GetPhoton()); // get matching electron
	if (Ele) {
		if (Ele->TrackNumber() >-1) {
			hist.TrkZ0->Fill(Ele->Z0());
	 		hist.EoverP->Fill(Ele->EOverP());
	 		hist.TrkPt->Fill(Ele->TrackPt());
	 		hist.TrkIso->Fill(Ele->TrackIso());
		}
	}
	hist.Iso4->Fill(sp->GetIso4());
	hist.Ces2Strip->Fill(sp->GetCesStripE2());
	hist.Ces2Wire->Fill(sp->GetCesWireE2());
	hist.XCes->Fill(sp->GetXCes());
	hist.ZCes->Fill(sp->GetZCes());
	hist.EmTime->Fill(sp->GetEmTime());
}

/*-------------------------------------------------------------------*/
// Loose Photon-like Electron ID CUTS
// The only difference between loose and tight is the E/P cut
// There isn't a similar set of cuts for plug-electrons!
/*-------------------------------------------------------------------*/
unsigned int TagLooseElectrons::GetPhoLikeEleIdWord(TStnPhoton* Pho, TStnEvent* event,
											TStnTrackBlock* fTrkBlock, double zvx)
{
	if (!Pho) {
		StdOut(__FILE__,__LINE__,3,"TStnPhoton is NULL. pl. check. exiting");
		exit (1);
	}
	if (!event) {
		StdOut(__FILE__,__LINE__,3,"TStnEvent is NULL. pl. check. exiting");
		exit (1);
	}
	if (!fTrkBlock) {
		StdOut(__FILE__,__LINE__,3,"TStnTrackBlock is NULL. pl. check. exiting");
		exit (1);
	}

	
	int iIdWord = 0x0;
	
	//______________________________________ only CEM photons with Et>7 GeV
	if (Pho->Detector() != 0) iIdWord |= kCentralBit;

	if (Pho->Etc() < fCemMinEt || Pho->Etc() > fCemMaxEt) iIdWord |= kEtCorrBit;

	//______________________________________ HADEM cut using 3 towers
	float cutMin_HADEM = 0.0;
	float cutMax_HADEM = 0.055 + 0.00045 * (Pho->Momentum()->E());
	if ((Pho->HadEm() < cutMin_HADEM) || (Pho->HadEm() > cutMax_HADEM)) iIdWord |= kHadEmBit;
	
	//______________________________________ CES Chi^2 cut
	if ( (Pho->Chi2() < 0.0) || (Pho->Chi2() > 20.0) ) iIdWord |= kChi2Bit;

	//______________________________________ N3D cut
	// it is not an electron if N3d=0
	if ((Pho->N3d())==0) iIdWord |= kN3d0Bit;

	//______________________________________ second N3D cut
	// no more than 2 tracks (1st--ele,2nd--to match cuts for pho)
	if ((Pho->N3d())>2)  iIdWord |= kN3d2Bit;

	//______________________________________ cut on 2nd max Pt track in cluster if N3D=2
	if ((Pho->N3d())==2) {
		float trkPtcut_max=1.0+0.005*(Pho->Etc());
		if ((Pho->Pt2())>trkPtcut_max) iIdWord |= k2ndTrkPtBit;
	}

	TStnElectron* Ele=TPhotonUtil::MatchElectron(event,Pho); // get matching electron

	//______________________________________ E/p cut for electron
	if (Ele==NULL) {
		iIdWord |= kMatchEleBit;
	} else {
		// there should be a track
		if (Ele->TrackNumber()<0) iIdWord |= kEleTrkBit;
		else {
			// only events with ele from best vertex
			if (fabs(Ele->Z0()-zvx)>3.0) iIdWord |= kEleTrkZ0Bit;

			int ele_trk = Ele->TrackNumber();
			if (!fTrkBlock) {
				StdOut(__FILE__,__LINE__,3,"You gave me a NULL pointer to TrackBlock. Pl. check");
				exit(1);
			}
			TStnTrack* Trk = fTrkBlock->Track(ele_trk);

			if ((Trk->Algorithm())==2) {
				if ((Ele->TrackBcPt())<50.0 &&
						(Ele->EOverP()>2.)) iIdWord |= kEoverPBit;
			} else {
				if ((Ele->TrackPt())<50.0 &&
						(Ele->EOverP()>2.)) iIdWord |= kEoverPBit;
			}
		}
	}


	//______________________________________ CalIso4 cut

	float cutMax_CalIso=0.0;

	if ((Pho->Etc()) < 20.0) cutMax_CalIso = 0.1 * (Pho->Etc());
	else cutMax_CalIso = 2.0 + 0.02 * ( Pho->Etc() - 20.0 );

	if ((Pho->EIso4(2)) > cutMax_CalIso) iIdWord |= kIso4Bit;

	//______________________________________ TrkIso4 cut
	float trkIso=Pho->SumPt4() - Pho->Pt();

	if (trkIso < 0.0) trkIso = Pho->SumPt4();

	if (trkIso < (trkIso > (2.0 + 0.005 * Pho->Etc())) ) iIdWord |= kTrkIso;

	//______________________________________ Energy of 2nd CES cluster (Wire & Strip)
	float _Et2ndCES = 0.0;
	if ((Pho->CesStripE2()) > (Pho->CesWireE2())) _Et2ndCES = (Pho->CesStripE2()) * (Pho->SinTheta());
	else _Et2ndCES = (Pho->CesWireE2()) * (Pho->SinTheta());

	float cutMax_2ndCes = 0.0;
	if ((Pho->Etc()) < 18.0) cutMax_2ndCes = 0.14 * (Pho->Etc());
	else cutMax_2ndCes = 2.4 + 0.01 * (Pho->Etc());

	if (_Et2ndCES > cutMax_2ndCes) iIdWord |= k2ndCesBit;

	//_______________________________________ Fiducial cuts
	if (fabs(Pho->XCes()) > 21.0) iIdWord |= kXCesBit;
	if (fabs(Pho->ZCes()) < 9.0 || fabs(Pho->ZCes()) > 230.0) iIdWord |= kZCesBit;

	return iIdWord;

} // Photon-like Electron Cuts



/*-------------------------------------------------------------------*/
// Standard Loose Central Electron ID CUTS
/*-------------------------------------------------------------------*/
unsigned int TagLooseElectrons::GetStdCenLooseEleID(TStnElectron* Ele, TStnTrackBlock* fTrkBlock, const int Nvtx)
{
	if (!Ele) {
		StdOut(__FILE__,__LINE__,3,"TStnElectron is NULL. pl. check. exiting");
		exit (1);
	}
	if (!fTrkBlock) {
		StdOut(__FILE__,__LINE__,3,"TStnTrackBlock is NULL. pl. check. exiting");
		exit (1);
	}
	
 	int iIdWord = 0;
  
	if ((Ele->DetectorCode()) != 0) iIdWord |= kSLC_DetBit;
 	if ( CorrEtEle(Ele) < fCemMinEt || (CorrEtEle(Ele)) > fCemMaxEt) iIdWord |= kSLC_EtBit;

   if((Ele->FidEleSmx()) != fCemFidCut) iIdWord |= kSLC_FidBit;  //Sasha, uses two mtds. Smx based and track based.  08-28-2008

	if (Ele->TrackNumber() < 0) iIdWord |= kSLC_TrkBit;
	else
	{
		int ele_trk = Ele->TrackNumber();
		TStnTrack* Trk = fTrkBlock->Track(ele_trk);
		if (!Trk)
		{
			StdOut(__FILE__,__LINE__,3,"Did not find a track for electron. Pl. check");
			exit(1);
		} else {
			if (fabs(Ele->Z0()) < fCemMinTrkZ0 || fabs(Ele->Z0()) > fCemMaxTrkZ0) iIdWord |= kSLC_TrkZ0Bit; // note, not beam-constrained track !!!
			if ((Trk->Algorithm()) == 2)
			{
				if ((Ele->TrackBcPt()) < fTrkMinPt || (Ele->TrackBcPt()) > fTrkMaxPt) iIdWord |= kSLC_TrkPtBit; // note, beam-constrained track !!! 
			} else {
				if ((Ele->TrackPt()) < fTrkMinPt || (Ele->TrackPt()) > fTrkMaxPt) iIdWord |= kSLC_TrkPtBit; // note, not beam-constrained track !!! 
    	}
  
		}
	}

	if ((Ele->Nasl()) < fCemMinCotAxSl || (Ele->Nasl()) > fCemMaxCotAxSl) iIdWord |= kSLC_TrkNaslBit;
	if ((Ele->Nssl()) < fCemMinCotStSl || (Ele->Nssl()) > fCemMaxCotStSl) iIdWord |= kSLC_TrkNsslBit;	
	if ((Ele->HadEm()) < (fCemMinHADEM[0]+fCemMinHADEM[1]*(CorrEnergyEle(Ele))) || 
     		(Ele->HadEm()) > (fCemMaxHADEM[0]+fCemMaxHADEM[1]*(CorrEnergyEle(Ele)))) iIdWord |= kSLC_HadEmBit;
  	if ((CorrIsoFr(Ele, Nvtx)) < fCemMinIsoFr || (CorrIsoFr(Ele, Nvtx)) > fCemMaxIsoFr) iIdWord |= kSLC_IsoFrBit;
  
  return iIdWord;
  
}


//____________________________ Sasha's routine to get Phoenix track by using Ray's methods
TStnTrack* TagLooseElectrons::GetPnxTrack(TStnEvent* event, TStnElectron* ele) {
  TStnElectron* phele = TPhotonUtil::MatchPhoenixElectron(event,ele);
  TStnTrack*    phTrk = TPhotonUtil::PhoenixTrack(event,phele);
  return phTrk;
}




/*-------------------------------------------------------------------*/
// Standard Loose Plug Electron ID CUTS
// I need to varify all these cuts with sasha or from some note
/*-------------------------------------------------------------------*/
unsigned int TagLooseElectrons::GetStdPlugLooseEleID(TStnElectron* Ele, TStnEvent* event, const int Nvtx)
{
	if (!Ele) {
		StdOut(__FILE__,__LINE__,3,"TStnPhoton is NULL. pl. check. exiting");
		exit (1);
	}
	if (!event) {
		StdOut(__FILE__,__LINE__,3,"TStnEvent is NULL. pl. check. exiting");
		exit (1);
	}
	
	int iIdWord = 0x0;

	double Et_em = CorrEtEle(Ele);
	
  if (Ele->DetectorCode() != 1) iIdWord |= kSLP_DetBit;
  if (Et_em < fPemMinEt || Et_em > fPemMaxEt) iIdWord |= kSLP_EtBit;   
  if (fabs(Ele->PesEta()) < fPemMinEtaPes || fabs(Ele->PesEta()) > fPemMaxEtaPes) iIdWord |= kSLP_PesEtaBit;
  if ((Ele->HadEm()) < fPemMinHADEM || (Ele->HadEm()) > fPemMaxHADEM) iIdWord |= kSLP_HadEmBit;
  if ((Ele->Pem3x3FitTower())==fPem3x3FitCut) iIdWord |= kSLP_Pem3x3Bit;
  if ((Ele->Chi2Three()) < fPemMin3x3Chi2 || (Ele->Chi2Three()) > fPemMax3x3Chi2) iIdWord |= kSLP_Chi2ThreeBit;
  if ((Ele->Pes5x9(0)) < fPemMinPes5x9U || (Ele->Pes5x9(0)) > fPemMaxPes5x9U) iIdWord |= kSLP_Pes5x9UBit;
  if ((Ele->Pes5x9(1)) < fPemMinPes5x9V || (Ele->Pes5x9(1)) > fPemMaxPes5x9V) iIdWord |= kSLP_Pes5x9VBit;
  if ((CorrIsoFr(Ele, Nvtx)) < fPemMinIsoFr || (CorrIsoFr(Ele, Nvtx)) > fPemMaxIsoFr) iIdWord |= kSLP_PemIsoFrBit;
  TStnPhoton* myPho=TPhotonUtil::MatchPhoton(event,Ele);
  if ((myPho->PesDeltaR()) < fPemMinPesdR || (myPho->PesDeltaR()) > fPemMaxPesdR) iIdWord |= kSLP_PesDelRBit;

  return iIdWord;
}

//____________________________________________________________________________
// returns Iso/Et(ele) corrected for leakage & nvx12
double TagLooseElectrons::CorrIsoFr(TStnElectron* ele, const int Nvtx) {
  return CorrIso(ele, Nvtx)/ele->Et();
}
//____________________________________________________________________________
// returns Iso corrected for leakage & nvx12
double TagLooseElectrons::CorrIso(TStnElectron* ele, const int NVtx)
{
	int nvx12= NVtx; //num of class12 vertices  
	double Iso=ele->Iso();
	Iso=Iso - (nvx12>1? 0.3563*(nvx12-1) : 0.0) - ele->LeakageEnergy();
	return Iso;
}



//____________returns corrected E of electrons (includes all corrections)_________ 
double TagLooseElectrons::CorrEnergyEle(TStnElectron* ele)
{
	return (ele->Et()>0.0? ele->E()*CorrEtEle(ele)/(ele->Et()) : ele->E());
}

//____________returns corrected Et of electrons (includes all corrections)_________ 
double TagLooseElectrons::CorrEtEle(TStnElectron* ele)
{
	double et=ele->Et();
	double corFactor= 1.0;
	corFactor=ele->Etcor(); // provide that fEtcor has already been corrected by Ray's CorrectElectronEnergy()
	et *= corFactor;
	return et;
}




//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TagLooseElectrons::EndJob() {

	if (GetSummaryStat()) return 0;

	std::string sMsg, sMsg2,sMsg3;
	if (GetMode() == 0) sMsg = "0 - Tag Only and Tag all photons.";
	if (GetMode() == 1) sMsg = "1 - Pass if only one electron is found";
	if (GetMode() == 2) sMsg = "2 - Pass if only two electrons are found.";
	
	printf("[TLE:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[TLE:01:] Events Processed------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "[TLE:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;
	std::cout << "[TLE:03:] Evts w/Elec (PhoLike OR STD) = " << counter.electronEvents 
						<< " [PhoLike Evts= " << counter.PhoLikeEleEvts 
						<< ", STD Evts= "<< counter.StdEleEvts << "]" << std::endl;
	std::cout << "[TLE:04:] Electrons found  ----------- = " 
					<< "[" << counter.PhoLikeEles << "-PhoLike][" 
					<< counter.StdCEMEles << "-STD CEM][" 
					<< counter.StdPLUGEles << "-STD PLUG]"<< std::endl;
	printf("---------------------------------------------------\n");
	return 0;
}
