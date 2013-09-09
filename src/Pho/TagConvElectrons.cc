///////////////////////////////////////////////////////////
//this module will tag the conversion electrons in super //
// photon list.                                          //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

/* $Id: TagConvElectrons.cc,v 1.7 2009/11/10 19:45:14 samantha Exp $
 * $Log: TagConvElectrons.cc,v $
 * Revision 1.7  2009/11/10 19:45:14  samantha
 * ADDED: A comment about this module and Auto CVS Log info.
 *
 *
 */

#include "samantha/Pho/TagConvElectrons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include <iostream>


ClassImp(TagConvElectrons)

//_____________________________________________________________________________
TagConvElectrons::TagConvElectrons(const char* name, const char* title):
  TStnModule(name,title), qMc(false),
  bNoSummary(false)
{
	std::cout << "Hello I am TagConvElectrons module" << std::endl;
}

//_____________________________________________________________________________
TagConvElectrons::~TagConvElectrons() {
}

//_____________________________________________________________________________
void TagConvElectrons::SaveHistograms() {
}

//_____________________________________________________________________________
void TagConvElectrons::BookHistograms()
{
	DeleteHistograms();
}


//_____________________________________________________________________________
int TagConvElectrons::BeginJob()
{
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fTrackBlock 	= (TStnTrackBlock*)    RegisterDataBlock("TrackBlock","TStnTrackBlock");
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");

	BookHistograms();

	InitSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (InitSpMod == NULL) {
		StdOut(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
		return 0;
	}
 
 	if (!fTrackBlock) {
		StdOut(__FILE__,__LINE__,3,"Could not find TrackBlock!");
		return 0;
	}
 	if (!fElectronBlock) {
		StdOut(__FILE__,__LINE__,3,"Could not find ElectronBlock!");
		return 0;
	}


	// initialize vars
	counter.evtsRunOver 		= 0;
	counter.evtsPassModule 	= 0;
	counter.convEleFound		= 0;

	return 0;
}

//_____________________________________________________________________________
int TagConvElectrons::BeginRun()
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
int TagConvElectrons::Event(int ientry)
{
	SetPassed(1);
	counter.evtsRunOver++;
	if (GetPassed()) counter.evtsPassModule++;
  	
  	FillDataBlocks(ientry);
  
	int NsuperPho = InitSpMod->GetSuperPhoSize();
	TStnEvent* event = GetEvent();
	for (int i =0; i < NsuperPho; i++) {
		TStnElectron* Ele = TPhotonUtil::MatchElectron(event,InitSpMod->GetSuperPhoton(i)->GetPhoton()); // get matching electroni
		double minsep =0, mindt=0, radius=0;
		bool bConvId = IsConversionElectron(Ele, minsep,mindt,radius);
		InitSpMod->GetSuperPhoton(i)->SetConversionId(bConvId);
		if (bConvId) {
			counter.convEleFound++;
		}
	}

	return 0;

} // Event

/*-------------------------------------------------------------------*/
// Isconversion electrons
/*-------------------------------------------------------------------*/
bool TagConvElectrons::IsConversionElectron(TStnElectron* ele, 
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
double TagConvElectrons::ConversionSep(TStnTrack* t0, TStnTrack* t1, double & rconv) {

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
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TagConvElectrons::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void TagConvElectrons::FillDataBlocks(int ientry)
{
	fElectronBlock->GetEntry(ientry);
	fTrackBlock->GetEntry(ientry);
}

/*===================================================================*/
//	Create a folder in"Hist" 
/*===================================================================*/
TFolder* TagConvElectrons::GetHistoFolder(char *name, char* title)
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
//  END JOB SUMMARY
//_____________________________________________________________________________
int TagConvElectrons::EndJob() {

	if (GetSummaryStat()) return 0;
		

	printf("[TCE:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[TCE:01:] Events Processed ----------- = " << counter.evtsRunOver << std::endl;
	std::cout << "[TCE:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;
	std::cout << "[TCE:03:] Conversion Electrons Found - = " << counter.convEleFound << std::endl;
	printf("---------------------------------------------------\n");
	return 0;
}
