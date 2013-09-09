///////////////////////////////////////////////////////////
//this module will tag the tight electrons the super     //
// photon list.                                          //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#include <iostream>
#include <fstream>
#include <TMath.h>
#include <numeric>
#include <vector>
#include <algorithm>
#include "Stntuple/loop/TStnAna.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/obj/TStnElectron.hh"
#include "Stntuple/obj/TStnTrack.hh"

ClassImp(TagElectrons)

//_____________________________________________________________________________
TagElectrons::TagElectrons(const char* name, const char* title):
  TStnModule(name,title),
  	iMode(0),
	bRunPermit(true),
  bNoSummary(false)
{

	float EtMin = 7.0;
	float EtMax = 1500.0;
	
	//std::cout << "Hello I am TagElectrons module" << std::endl;
	fEtMin = EtMin;
	fEtMax = EtMax;
	//fCemMinHADEM[0]=0.0; 
	//fCemMinHADEM[1]=0.0; 
	//fCemMaxHADEM[0]=0.055;                      
	//fCemMaxHADEM[1]=0.00045;                      
	//fCemMinIsoFr=-1.0;  
	//fCemMaxIsoFr=0.1;  
	fCemMinLshr=0.0;    
	fCemMaxLshr=0.2;
	fCemMinEP=0.0;    
	fCemMaxEP=2.0;
	fCemMinChi2strip=-10.0E6; 
	fCemMaxChi2strip=10.0;
	fCemMindXQ=-3.0;
	fCemMaxdXQ=1.5; 
	fCemMinCesdZ=0.0; 
	fCemMaxCesdZ=3.0;
	//fCemMinTrkZ0=0.0; 
	//fCemMaxTrkZ0=60.0;
	//fCemMinCotAxSl=3; // cut on number of COT Axial SL with 5 or more hits 
	//fCemMaxCotAxSl=8;
	//fCemMinCotStSl=2; // cut on number of COT Stereo SL with 5 or more hits 
	//fCemMaxCotStSl=8;
}

//_____________________________________________________________________________
TagElectrons::~TagElectrons() {
}

//_____________________________________________________________________________
void TagElectrons::SaveHistograms() {
}

//_____________________________________________________________________________
void TagElectrons::BookHistograms()
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
}

//_____________________________________________________________________________
int TagElectrons::BeginJob()
{
	initSp = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initSp == NULL) {
		StdOut(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
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
	counter.electronEvents	= 0;;

	return 0;
}

//_____________________________________________________________________________
int TagElectrons::BeginRun()
{

	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TagElectrons::Event(int ientry)
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
	int Nsp = initSp->GetSuperPhoSize();

	for (int i = 0; i < Nsp; i++) {
		//FillHistograms(EleB4Cut, initSp->GetSuperPhoton(i));  // lets make the plots b4 cuts. anyone who wants
																								// after cuts can do that in that module
		//get the photon-like electron id word for central photons
		int eid = GetPhoLikeEleIdWord(initSp->GetSuperPhoton(i)->GetPhoton(),event, fTrackBlock, zvx);
		if (eid == 0) {
			iNElectrons++;
			//FillHistograms(EleA4Cut, initSp->GetSuperPhoton(i));		//need to redo these. either use SP or have struct
		}
		initSp->GetSuperPhoton(i)->SetTightElectronId(eid);
	
		//now get the std. tight electron id cuts for central and plug photons
		int iStdTightEleId = -1;

		TStnElectron* Ele = TPhotonUtil::MatchElectron(event,initSp->GetSuperPhoton(i)->GetPhoton()); // get matching electron
		
		if (initSp->GetSuperPhoton(i)->GetDetector() == 0 )		//central
		{
			iStdTightEleId = GetStdCenTightEleID(Ele, fTrackBlock);
			initSp->GetSuperPhoton(i)->SetStdTightElectronId(iStdTightEleId);
		
		} else if (initSp->GetSuperPhoton(i)->GetDetector() == 1 )		//plug : has no set of tight cuts!!
		{
			//iStdTightEleId = GetStdPlugTightEleID(Ele,event);
			//initSp->GetSuperPhoton(i)->SetStdTightElectronId(iStdTightEleId);
		}

		
	}
	
	if (iNElectrons >0) counter.electronEvents++;

	if (GetMode() == 0) SetPassed(1);
	else if (GetMode() == 1  && iNElectrons == 1)	SetPassed(1);
	else if (GetMode() == 2  && iNElectrons == 2)	SetPassed(1);

	if (GetPassed()) counter.evtsPassModule++;
	

	return 0;

} // Event


/*-------------------------------------------------------------------*/
// Standard Tight Central Electron ID CUTS
/*-------------------------------------------------------------------*/
unsigned int TagElectrons::GetStdCenTightEleID(TStnElectron* Ele, TStnTrackBlock* fTrkBlock)
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
  
  	//since I'll be setting the same variables in the SuperPhoton for central and plug
	//electrons, the detector bit will be common to both. and this is kind of redundant.
	if (Ele->DetectorCode() != 0) iIdWord |= kSTC_DetectorBit;

	if (Ele->Etcor() < fEtMin || Ele->Etcor() > fEtMax) iIdWord |= kSTC_EtCorrBit;
  
  	if (fabs(Ele->Lshr()) < fCemMinLshr 
		|| fabs(Ele->Lshr()) > fCemMaxLshr) iIdWord |= kSTC_LshrBit;

	int trkNum = Ele->TrackNumber();
	
  	if (trkNum < 0) iIdWord |= kSTC_TrkBit;
	else 
	{
		TStnTrack *Trk = fTrackBlock->Track(trkNum);
		
  		if (Trk->Algorithm() == 2)
    	{
	   	if ( (Ele->TrackBcPt() < 50.0) 
				&& ( (Ele->EOverP() < fCemMinEP) || (Ele->EOverP()>fCemMaxEP) ) ) iIdWord |= kSTC_EoverPBit;
				
      	if ( (fabs(Ele->BcDelZ()) < fCemMinCesdZ) 
				|| (fabs(Ele->BcDelZ()) >fCemMaxCesdZ) ) iIdWord |= kSTC_DelZBit; // note, beam-constrained track !!!
				
      	if ( (Ele->BcDelXQ() < fCemMindXQ) 
				|| (Ele->BcDelXQ() > fCemMaxdXQ) ) iIdWord |= kSTC_DelXQBit; // note, beam-constrained track !!!
    	} else
 		{
      	if ( (Ele->TrackPt() < 50.0) 
				&& ( (Ele->EOverP() < fCemMinEP) || (Ele->EOverP() > fCemMaxEP)) ) iIdWord |= kSTC_EoverPBit;
				
      	if ( (fabs(Ele->DelZ()) < fCemMinCesdZ) 
				|| (fabs(Ele->DelZ()) > fCemMaxCesdZ) ) iIdWord |= kSTC_DelZBit; // note, NOT beam-constrained track !!!
				
      	if ( (Ele->DelXQ() < fCemMindXQ)
				|| (Ele->DelXQ() > fCemMaxdXQ) ) iIdWord |= kSTC_DelXQBit; // note, NOT beam-constrained track !!!

		} //if trackalgo == 2
	} // if trkNum > 0
  
	if ( (Ele->Chi2Strip() < fCemMinChi2strip) 
		|| (Ele->Chi2Strip() > fCemMaxChi2strip)) iIdWord |= kSTC_Chi2StripBit;

	return iIdWord;
}



/*-------------------------------------------------------------------*/
void TagElectrons::FillHistograms(Histograms::EmObjHists_t& hist, SuperPhoton* sp)
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
		if (Ele->TrackNumber() >0) {
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
// Tight Photon-like Electron ID CUTS
// The only difference between loose and tight is the E/P cut
// There isn't a similar set of cuts for plug-electrons!
/*-------------------------------------------------------------------*/
unsigned int TagElectrons::GetPhoLikeEleIdWord(TStnPhoton* Pho, TStnEvent* event,
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

	if (Pho->Etc() < fEtMin) iIdWord |= kEtCorrBit;

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
						(Ele->EOverP()<0.8 || Ele->EOverP()>1.2)) iIdWord |= kEoverPBit;
			} else {
				if ((Ele->TrackPt())<50.0 &&
						(Ele->EOverP()<0.8 || Ele->EOverP()>1.2)) iIdWord |= kEoverPBit;
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

} // Electron Cuts

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TagElectrons::EndJob() {

	if (GetSummaryStat()) return 0;

	std::string sMsg, sMsg2,sMsg3;
	if (GetMode() == 0) sMsg = "0 - Tag Only and Tag all photons.";
	if (GetMode() == 1) sMsg = "1 - Pass if only one electron is found";
	if (GetMode() == 2) sMsg = "2 - Pass if only two electrons are found.";
	
	printf("[TEL:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[TEL:01:] Events Processed------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "[TEL:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;
	std::cout << "[TEL:03:] Events With an Electron ---- = " << counter.electronEvents << std::endl;
	printf("---------------------------------------------------\n");
	return 0;
}
