//_____________________________________________________________________________
// cp2 object  
// Removed the unrealted and unused headers, varibles and comments.
// compared with a reference after cleaning up is done to confirm
// that I did not remove anything critical. - sam 05-04-2008
// Lot of changes. Doing all in one pass now.
//_____________________________________________________________________________

/*
 * $Id: cp2.cc,v 1.8 2010/07/02 17:36:51 samantha Exp $
 * 	$Log: cp2.cc,v $
 * 	Revision 1.8  2010/07/02 17:36:51  samantha
 * 	MODIFIED: To look at events from different triggers to figure out the difference
 * 	in the LERs I get from my ntuples and Calib.exe. All this trig stuff is
 * 	commented out for now to do the P30 calibration.
 * 	Also added some event counters.
 * 	
 * 	Revision 1.7  2009/08/07 21:37:33  samantha
 * 	can give name to the histogram that save all the final hists. last DB LER must be provided by the user. Still stuff hanging in here to figure out a way to clearly ID the HV off periods.
 * 	
 * 	Revision 1.6  2009/04/02 19:36:44  samantha
 * 	This is modified version of 1.4 . renamed some of the hist that was supposed to be temp that started with sam_xxxx . I am still testing a good way to figure out the HV status.
 * 	
 *
 */
//_____________________________________________________________________________
//What this in simple terms:
//1. Use L3 COT tracks and extrapolate them to CP2 pads.
//2. Find the pad closest to the track and get the ADC count for that pad (this is from Qbank).
//3. Subtract 24 ADC counts from Qbank value to recreate the RawBank ADC count.
//4. If the count is below 40 (pedestal?) reject them.
//5. Else Fill out corresponding histograms to compute LERs
//6. Then to check the CP2 HV status, loop over all the CP2 pads and derive avg. hist per pad. if this is <0.01 (or very little amount of hits, classify it as a HV off event.
//_____________________________________________________________________________

#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>
#include "Riostream.h"
#include <sstream>
#include <TLine.h>
#include <set>
#include <TCanvas.h>

#include "cp2.hh"
#include "Stntuple/obj/TL3Track.hh"
#include "TStyle.h"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TStnTriggerTable.hh"
#include "Stntuple/data/TTl3d.hh"
#include "Stntuple/obj/TStnDBManager.hh"

using namespace std;

ClassImp(cp2)
	//_____________________________________________________________________________
cp2::cp2(const char* name, const char* title):
		TStnModule(name,title),
		sNewDataFile("cp2traler-new_temp.txt"),
		sMainLogFile("cpMainLog_Default.txt"),
		sFinalHistFile("CP2CalibHists.root"),
		fLERlastDB(0)
{
	//for final RatioLER Plots	// make these values read from text file
	// database ler1 was teststand
	fLERfirst = 219.87;  // December             database ler2
	//fLERlastDB = 189.827;  //256904 Jan 29, 2008   ler8 in database
	fLERnew 	  = 0000.00;	//will be set at the end of the job

}

//_____________________________________________________________________________
cp2::~cp2() {
}

//_____________________________________________________________________________
void cp2::BookHistograms() {



	fHist.cp2hit = new TH2F("cp2hit","",48,-0.5,47.5,54,-0.5,53.5);
	fHist.cp2hit->SetTitle("CPR2 Pad Hits from Tracks");
	fHist.cp2hit->SetXTitle("Wedge Number");
	fHist.cp2hit->SetYTitle("CPR2 Pad Number");

	fHist.cp2Weight = new TH2F("cp2Weight","",48,-0.5,47.5,54,-0.5,53.5);
	fHist.cp2Weight->SetTitle("CPR2 Pad Hits Weighted by Pulse Height");
	fHist.cp2Weight->SetXTitle("Wedge Number");
	fHist.cp2Weight->SetYTitle("CPR2 Pad Number");

	fHist.cp2ave = new TH2F("cp2ave","",48,-0.5,47.5,54,-0.5,53.5);
	fHist.cp2ave->SetTitle("Average pulse height");
	fHist.cp2ave->SetXTitle("Wedge Number");
	fHist.cp2ave->SetYTitle("CPR2 Pad Number");

	fHist.ymoncp2 = new TH2F("ymoncp2","",48,-0.5,47.5,54,-0.5,53.5);
	fHist.ymoncp2->SetTitle("CPR2 Hits YMon Threshold of 60 ADC Counts");
	fHist.ymoncp2->SetXTitle("Wedge Number");
	fHist.ymoncp2->SetYTitle("CPR2 Pad Number");

	fHist.west = new TH1F("west","",100,0,1000);
	fHist.west->SetXTitle("ADC Counts");
	fHist.west->SetTitle("Pulse Height From Tracks on West Side");

	fHist.east = new TH1F("east","",100,0,1000);
	fHist.east->SetXTitle("ADC Counts");
	fHist.east->SetTitle("Pulse Height From Tracks on East Side");

	fHist.cp2tra = new TH1F("cp2tra","",100,0,1000);
	fHist.cp2tra->SetXTitle("ADC Counts");
	fHist.cp2tra->SetTitle("Pulse Height From Tracks");

	fHist.rancp2 = new TH1F("rancp2","",100,0,1000);
	fHist.rancp2->SetXTitle("ADC Counts");
	fHist.rancp2->SetTitle("Random CPR2 Hits");

	fHist.zerocp2 = new TH1F("zerocp2","",48,-0.5,47.5);
	fHist.zerocp2->SetXTitle("Wedge Number");
	fHist.zerocp2->SetTitle("Wedges with Zero Pulse Height from a Track");

	fHist.zero3d = new TH2F("zero3d","",48,-0.5,47.5,54,-0.5,53.5);
	fHist.zero3d->SetTitle("Channels with Zero Pulse Height from a Track");
	fHist.zero3d->SetXTitle("Wedge Number");
	fHist.zero3d->SetYTitle("CP2 Pad Number");

	char name[100], title[100];
	for(int i=0; i<48; i++) {		// loop over wedges
		for(int p=0; p<54; p++) {		// loop over pads
			sprintf(name,"cp2pad_%i_%i",i,p);
			sprintf(title,"cp2pad[%i][%i]",i,p);
			fHist.cp2pad[i][p]= new TH1F(name,title,100,0,1000);
			fHist.cp2pad[i][p]->SetXTitle("CP2 ADC counts"); 
			fHist.cp2pad[i][p]->SetYTitle("Number of Events"); 
		}
	}

	//define Final LER histograms
		//these are settings to make the final plots all nice and greeny
/*		gStyle.SetCanvasColor(10); 
		gStyle.SetOptStat(1111111);
		gStyle.SetHistLineWidth(2); gStyle.SetFuncWidth(3); 
		gStyle.SetHistLineColor(1);
		gStyle.SetMarkerStyle(22); gStyle.SetMarkerSize(1.0);
		gStyle.SetMarkerColor(2); gStyle.SetTitleOffset(1.3,"Y");
		gStyle.SetPadLeftMargin(0.15);
		gStyle.SetPadTickX(1); gStyle.SetPadTickY(1); 
		gStyle.SetOptFit(111); 
		gStyle->SetPalette(1,0); 
		gStyle.SetTitleTextColor(4); gStyle.SetTitleSize(0.05,"X");
		gStyle.SetTitleSize(0.05,"Y"); gStyle.SetPadBottomMargin(0.15);
		gROOT->ForceStyle();
*/

	histLER.Ratio = new TH1F("LERratio","",100,0,3);
	histLER.Ratio->SetXTitle("New / Old");
	histLER.Ratio->SetTitle("Ratio of New to Old LERs");

	histLER.OldLER = new TH1F("oldLERs","",100,0.0,3.0);
	histLER.OldLER->SetXTitle("Old LERs");
	histLER.OldLER->SetLineColor(2);

	histLER.NewLER = new TH1F("newLERs","",100,0.0,3.0);
	histLER.NewLER->SetXTitle("New LERs");

	//add a 2D plot of wedge vs pad # weighted by this ratio
	histLER.PadHitsWgted_2D = new TH2F("PadsWgted2D","",48,-0.5,47.5,54,-0.5,53.5);
	histLER.PadHitsWgted_2D->SetTitle("CPR2 Pad Hits Weighted by the ratio");
	histLER.PadHitsWgted_2D->SetXTitle("Wedge Number");
	histLER.PadHitsWgted_2D->SetYTitle("CPR2 Pad Number");
	histLER.PadHitsWgted_2D->SetMaximum(2.0); 
	histLER.PadHitsWgted_2D->SetMinimum(0.0); 
	histLER.PadHitsWgted_2D->SetStats(kFALSE);

	//add a 1D plot of wedge weighted by this ratio
	histLER.WedgeHitsWgted = new TH1F("WedgeHitsWgted","",48,-0.5,47.5);
	histLER.WedgeHitsWgted->SetTitle("CPR2 Wedge Hits Weighted by the LER Ratio");
	histLER.WedgeHitsWgted->SetXTitle("Wedge Number");
	histLER.WedgeHitsWgted->SetMaximum(1.5); 
	histLER.WedgeHitsWgted->SetMinimum(0.5); 
	histLER.WedgeHitsWgted->SetStats(kFALSE);

	//add a 1D plot of pad weighted by this ratio
	histLER.PadHitsWgted_1D = new TH1F("PadHitsWgted_1D","",54,-0.5,53.5);
	histLER.PadHitsWgted_1D->SetTitle("CPR2 Pad Hits Weighted by the LER Ratio");
	histLER.PadHitsWgted_1D->SetXTitle("Pad Number");
	histLER.PadHitsWgted_1D->SetMaximum(1.5); 
	histLER.PadHitsWgted_1D->SetMinimum(0.5); 
	histLER.PadHitsWgted_1D->SetStats(kFALSE);

	//add a 1D plot of pad weighted by this ratio for a special wedge cut
	histLER.PadHitsWgted_special = new TH1F("PadHitsWgtedSpecial","",54,-0.5,53.5);
	histLER.PadHitsWgted_special->SetTitle("CPR2 Pad Hits Weighted by the LER Ratio");
	histLER.PadHitsWgted_special->SetXTitle("Pad Number");
	histLER.PadHitsWgted_special->SetMaximum(2.0); 
	histLER.PadHitsWgted_special->SetMinimum(0.); 
	histLER.PadHitsWgted_special->SetStats(kFALSE);



	histLER.sumet_large = new TH1F("largesumet","events with large sumet",1000,0.0,1000.0);
	histLER.sumet_small = new TH1F("smallsumet","events with small sumet",100,0.0,100.0);
	histLER.met_large = new TH1F("largemet","events with large met",1000,0.0,1000.0);
	histLER.met_small = new TH1F("smallmet","events with small met",100,0.0,100.0);
	histLER.ntracks_large = new TH1F("largentrks","events with large ntracks",200,0.0,200.0);
	histLER.ntracks_small = new TH1F("smallntrks","events with small ntracks",200,0.0,200.0);
	histLER.avghitsperpad = new TH1F("avghitsperpad","avg hist per pad",10000,0.0,1.0);


	histLER.sumetVsNtrks = new TH2F("sumetVsNtrks","SumEt Vs Ntrks",100,0,100,1000,0.0,1000.0);
	histLER.metVsNtrks = new TH2F("metVsNtrks","Met Vs Ntrks",100,0,100,1000,0.0,1000.0);
	histLER.avghitsperpadVssumet = new TH2F("avghitsperpadVssumet","avg hist per pad Vs SumEt",1000,0,1000,1000,0.0,1.0);
	histLER.avghitsperpadVsNtrks = new TH2F("avghitsperpadVsNtrks","avg hist per pad Vs Ntrks",100,0,100,1000,0.0,1.0);

	histLER.avghitsperpadNtrksRatio = new TH1F("avghitsperpadNtrksRatio","avg hist per pad/Ntks",10000,0.0,1.0);
	histLER.avghitsperpadSumetRatio = new TH1F("avghitsperpadSumetRatio","avg hist per pad/Sumet",10000,0.0,1.0);
	histLER.SumetNtrksRatio = new TH1F("SumetNtrksRatio","Sumet/Ntrks",1000,0.0,300);
	histLER.avghitsperpadVsNVtx = new TH2F("avghitsperpadVsNVtx","avg hist per pad Vs NVtx",10,0,10,1000,0.0,1.0);
	histLER.avghitsperpadVsNVtx12 = new TH2F("avghitsperpadVsNVtx12","avg hist per pad Vs NVtx12",10,0,10,1000,0.0,1.0);
}
//_____________________________________________________________________________
int cp2::BeginJob() 
{
	fHeaderBlock = (TStnHeaderBlock*) 
		RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fCp2DataBlock = (TCp2DataBlock*) 
		RegisterDataBlock("Cp2DataBlock","TCp2DataBlock");
	fL3SummaryBlock =   (TL3SummaryBlock*) RegisterDataBlock("L3SummaryBlock","TL3SummaryBlock");
	fMetBlock = (TStnMetBlock*) 
		RegisterDataBlock("MetBlock","TStnMetBlock");
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");
	fTriggerBlock 	= (TStnTriggerBlock*)  RegisterDataBlock("TriggerBlock","TStnTriggerBlock");

	assert(fHeaderBlock != NULL && "HeaderBlock not found!");
	assert(fCp2DataBlock != NULL && "Cp2DataBlock not found!");
	assert(fL3SummaryBlock != NULL && "L3SummaryBlock not found!");
	assert(fMetBlock != NULL && "MetBlock not found!");
	assert(fVertexBlock != NULL && "VertexBlock not found!");
	assert(fTriggerBlock != NULL && "TriggerBlock not found!");

	BookHistograms();


	if (fLERlastDB == 0)
	{
		std::cout << __FILE__ <<"::"<< __LINE__ << ":: Global LER corresponding to last DN entry not specified.! Exiting!." << std::endl;
		exit (1);
	}

	FILE *fpp=fopen(sLastDBDataFile.c_str(),"r");
	if (! fpp)
	{
		std::cout << __FILE__ << "::" << __LINE__ << "::ERROR! default file did not open! returning without making final hists!" << std::endl;
		exit (1);
	} else {
		std::cout << "Last DB entry data file specified is " << sLastDBDataFile << std::endl;
	}

	fclose(fpp);


	
	//OPEN LOG FILES
	mainLogFile = new ofstream(sMainLogFile.c_str());
	if (mainLogFile->good()) {
		std::cout << "Mail log file " << sMainLogFile <<" opened" << std::endl;
	} else {
		std::cout << __FILE__ <<"::"<< __LINE__ << ":: MainLogFile did not open.! Exiting!." << std::endl;
		exit (1);
	}


	if (sFinalHistFile.length()<1)
	{
		std::cout << "Fatal Error! No output file name specified for the final hists to be written! exiting!" << std::endl;
		exit (1);
	}

	trigbits.SetTriggerBlock(fTriggerBlock);

	iEvtsProc[0] = 0;
	iEvtsProc[1] = 0;
	return 0;
}

//_____________________________________________________________________________
int cp2::BeginRun() 
{
	int currRun =  fHeaderBlock->RunNumber();
//	trigbits.BeginRun();	
//	TStntuple::Init(currRun);


  TStnTriggerTable* ptab = (TStnTriggerTable*) 
    TStnDBManager::Instance()->GetTable("TriggerTable");

  if(!ptab) {
	  std::cout << "TPhoTrigBits: did not find TriggeTable\n" << std::endl;
    return 0;
  }
  //ptab->Print();

  int nl3 = ptab->NL3Triggers();
  const TStnTrigger* trig;
  iTrigBit = 0;
  for(int i=0; i<nl3; i++) {
    trig = ptab->GetTrigger(3,i);

    if( (  (trig->Name()).Contains(sTrigName.c_str(),TString::kIgnoreCase) ) ) 
	 {
      iTrigBit= trig->Bit();
		//std::cout << "Found trig bit = " << iTrigBit << std::endl;
		break;
	 }
  }
	 
  
	return 0;
}


//_____________________________________________________________________________
int cp2::Event(int ientry) 
{
	++iEvtsProc[0];

	//fTriggerBlock->GetEntry(ientry);
  	//TTl3d* tl3d = fTriggerBlock->Tl3d();
	//int bit = 0;
  	//if(tl3d!=NULL && tl3d->NPaths()>0)
	//{
   // if(iTrigBit>=0)    bit = tl3d->PathPassed(iTrigBit);
	//}
	
//	if (! bit) return 0;
	/*std::cout << "bit = " << bit << std::endl;
	if (! bit) { std::cout << "trig failed " << std::endl; return 0; }
	else std::cout << "PASSED " << std::endl;
	*/

	++iEvtsProc[1];

	fHeaderBlock->GetEntry(ientry);
	fCp2DataBlock->GetEntry(ientry);
	fL3SummaryBlock->GetEntry(ientry);
	//fMetBlock->GetEntry(ientry);
	//fVertexBlock->GetEntry(ientry);

	//check for trigger. for debugging Calib.exe stntuples
	//and mine. See elog:
	//trigbits.Event();
	//if (CheckMinBiasTrigger())
	


	
	// pad position - z cdts in local cdts? - sam
	/*
		layout of cp2 pads!

		phi 
		0 x x x x x .. 18 pads along z
		z   1 x x x x x ..
		2 x x x x x ..

		3 pads along phi

		below are the local! cdts of each pad (to the center of the pad)
		- sam
	 */

	// These are in cm 
	const float fZ_MIDDLE[54] ={8.5,8.5,8.5,21.,21.,21.,33.5,33.5,33.5,
		46.,46.,46.,58.5,58.5,58.5,71.,71.,71.,
		83.5,83.5,83.5,96.,96.,96.,108.5,108.5,108.5,
		121.,121.,121.,133.5,133.5,133.5,
		146.,146.,146.,158.5,158.5,158.5,
		171.,171.,171.,183.5,183.5,183.5,
		196.,196.,196.,208.5,208.5,208.5,
		221.,221.,221.};

	// These are in degrees since a wedge is 15 degrees
	const float fPHI_MIDDLE[54] ={3.31,7.50,11.69,3.31,7.50,11.69,
		3.31,7.50,11.69,3.31,7.50,11.69,
		3.31,7.50,11.69,3.31,7.50,11.69,
		3.31,7.50,11.69,3.31,7.50,11.69,
		3.31,7.50,11.69,3.31,7.50,11.69,
		3.31,7.50,11.69,3.31,7.50,11.69,
		3.31,7.50,11.69,3.31,7.50,11.69,
		3.31,7.50,11.69,3.31,7.50,11.69,
		3.31,7.50,11.69,3.31,7.50,11.69};


	Bool_t  debug 	  = 0;
	Float_t fCP2RAD   = 170.47;			// radius cp2 is located from beam line? - sam
	Float_t xyzcp2[3];
	int   isidecpr 	  = 99;
	int   iwedgecpr   = 99;
	int   ipadcpr 	  = 99;
	float delz;
	float delx;
	float delr;
	float bestdelr;  

	int l3ntrk = fL3SummaryBlock->NCotTrack(); //cout<<l3ntrk<<endl;

	for (int j=0; j<l3ntrk; j++)
	{
		TL3Track *l3trk 	= fL3SummaryBlock->CotTrack(j);
		double pt 			= fabs(l3trk->Pt());
		double charge 	= l3trk->Charge();
		double z0 			= l3trk->Z0();
		double curv 		= (0.0029979*1.4116*charge)/(2.0*pt);

		Int_t old_ieta = 0;
		Int_t old_iphi = 0;   

		TrackExtr(fCP2RAD,curv,l3trk->Phi0(),l3trk->Lambda(),l3trk->D0(),l3trk->Z0(),
				xyzcp2,old_ieta,old_iphi,debug);

		if (pt > 1.0) {
			float zGlobal = xyzcp2[2];
			isidecpr = 0;
			if (zGlobal > 0) isidecpr = 1;
			double phicpr = atan2(xyzcp2[1],xyzcp2[0]);
			if (phicpr < 0) phicpr = phicpr+(2.0*M_PI);
			// Now convert phiGlobal from radians to degrees
			float phiGlobal = 180.0 * phicpr / M_PI; 
			int tmpphi 	= (int) phiGlobal;
			iwedgecpr 	= phiGlobal/15;
			float local_phi = (tmpphi % 15) + (phiGlobal-tmpphi);
			bestdelr = 999.;

			for (int ipad=0; ipad<54; ipad++) {
				// delta-Z is in cm
				delz = fabs(zGlobal)-fZ_MIDDLE[ipad];
				// Find delta-local-phi and convert to cm
				delx = (local_phi-fPHI_MIDDLE[ipad])*(12.5/4.13);
				// delta-R is in cm
				delr = sqrt(delz*delz + delx*delx);
				if(delr<bestdelr) {
					bestdelr = delr;
					ipadcpr = ipad;
				}
			}

			Short_t rawdatp = fCp2DataBlock->GetPadData(isidecpr, iwedgecpr, ipadcpr) - 24;   //we are making the raw data bank from Qbank data 

			// Begin "zeros" plots with raw data to find bugs or bad hardware
			if(bestdelr<4 && rawdatp<40) fHist.zerocp2->Fill(iwedgecpr+isidecpr*24,1);
			if(bestdelr<4 && rawdatp<40) fHist.zero3d->Fill(iwedgecpr+isidecpr*24,ipadcpr,1);

			// Correct datp with teststand gains and sin(theta) and make some plots
			float zdif = fabs(zGlobal-z0); float hyp = sqrt((zdif*zdif)+(fCP2RAD*fCP2RAD));
			float sth  = fCP2RAD/hyp;
			Short_t datp = (rawdatp*sth);
			if (bestdelr < 4) fHist.cp2tra->Fill(datp,1);
			if (bestdelr < 4) fHist.cp2pad[iwedgecpr+isidecpr*24][ipadcpr]->Fill(datp,1);
			if (bestdelr < 4 && isidecpr==0) fHist.west->Fill(datp,1);
			bool chimney = (isidecpr==1 && iwedgecpr==5 && ipadcpr>47);
			if (bestdelr<4 && isidecpr==1 && !chimney) fHist.east->Fill(datp,1);

			// Begin LER calibration plots
			if (bestdelr<4 && rawdatp>50 && rawdatp<1000) fHist.cp2hit->Fill(iwedgecpr+isidecpr*24,ipadcpr,1);
			if (bestdelr<4 && rawdatp>50 && rawdatp<1000) fHist.cp2Weight->Fill(iwedgecpr+isidecpr*24,ipadcpr,datp);
		}
	}


	//This part is to identify the HV status.
	//We are saying it there is no activity in the detector
	//the detector is probably off! Which is giving us false negatives.
	//So we need to tighten it up or find a better method. -- 02-20-2009 sam
	//=======updated with Steve's permission. see email on 1-5-08 - sam

	int nrun = fHeaderBlock->RunNumber();
	int nsec = fHeaderBlock->SectionNumber();
	int nevt = fHeaderBlock->EventNumber();
	std::pair<unsigned, unsigned> runsec(nrun, nsec);
	processedRunSec.insert(runsec);

	int nhits=0;
	for (int iside=0; iside<2; iside++) {		//two side, east and west barrels - sam
		for (int iwed=0; iwed<24; iwed++) {		// 24 phi wedges -sam
			for(int ipad=0; ipad<54; ipad++) {	// 54 pads for on each phi wedge -sam
				Short_t rawdatp = fCp2DataBlock->GetPadData(iside, iwed, ipad) - 24;		// what is this 24?   //we are making the raw data bank from Qbank data 
				if(rawdatp>10) nhits++;
				//rawdatap is prob. the raw ADC count, and 24 is the pedestal?
			}
		}
	}

	//float sumet = fMetBlock->Sumet(0);
	//float met = fMetBlock->Met(0);		// it seems the Met information is not calculated when I make the Stntuples!
	//new 1-5-08
	float hits=nhits/2592.;   // devide by the total number of pads(=2*24*54) to get the avg. per pad? - sam
	histLER.avghitsperpad->Fill(hits);	
	if(hits<0.01) {
		cout<<nrun<<"::" << nsec << "::"<<nevt<<"::"<<hits<<endl;		// this is where the zeros are spit out indicating HV down!
		foundHVoff.insert(runsec);
		//histLER.sumet_small->Fill(sumet);
		//histLER.met_small->Fill(met);
		//histLER.ntracks_small->Fill(l3ntrk);
	}

	//this to figure out a better HV checking 
	//histLER.sumet_large->Fill (sumet);
	//histLER.met_large->Fill (met);
	//stdcout << "sumet=" << sumet << std::endl;
	//std::cout << "met=" << met << std::endl;
	//histLER.ntracks_large->Fill(l3ntrk);
	//histLER.sumetVsNtrks->Fill(l3ntrk,sumet);
	//histLER.metVsNtrks->Fill(l3ntrk,met);
	//histLER.avghitsperpadVssumet->Fill(sumet,hits);	
	//histLER.avghitsperpadVsNtrks->Fill(met,hits);	
	
	//if (l3ntrk>0)	histLER.avghitsperpadNtrksRatio->Fill(hits/l3ntrk);	
	//if (sumet>0) histLER.avghitsperpadSumetRatio->Fill(hits/sumet);	
	//if (l3ntrk>0)histLER.SumetNtrksRatio->Fill(sumet/l3ntrk);	


	//int iNVertices = fVertexBlock->NVertices();
	//int iN12Vertices = 0;

	//for (int iVtx = 0; iVtx < iNVertices; ++iVtx) {
	//	TStnVertex* vert = fVertexBlock->Vertex(iVtx);
	//	if (vert->VClass() >= 12) ++iN12Vertices;
	//}


	//histLER.avghitsperpadVsNVtx->Fill(iNVertices,hits);	
	//histLER.avghitsperpadVsNVtx12->Fill(iN12Vertices,hits);	


	return 0;
}



//_____________________________________________________________________________
int cp2::EndJob() {
	printf("----- end of job: ---- %s\n",GetName());

	float totchan=0., globalave=0., ave, cfnew[48][54];

	for (int i=0; i<48; i++) {
		for (int p=0; p<54; p++) {
			cfnew[i][p]=1.0;
			if (fHist.cp2hit->GetBinContent(i+1,p+1)!=0) {
				ave = fHist.cp2Weight->GetBinContent(i+1,p+1)/(double)fHist.cp2hit->GetBinContent(i+1,p+1);
				fHist.cp2ave->Fill(i,p,ave);
				cfnew[i][p] = ave;
				totchan++;
				globalave = globalave+ave;
			}
		}
	}
	globalave = globalave/totchan;

	cout<<"global average = "<<globalave<<" totchan "<<totchan<<endl;
	fLERnew = globalave;

	float bgg;
	FILE *fpp=fopen(sNewDataFile.c_str(),"w");
	for (int i=0;i<2592;i++) {
		int id0 = i%54 +256*((i/54)%24);
		if (i/54>23) id0+=64;
		int ipad = (id0 & 0x0000003f); 
		int iwed = ((id0 >> 8) & 0x0000001f);
		int iside = ((id0 >> 6) & 0x00000001);
		if (cfnew[iwed+iside*24][ipad]==1.0) cfnew[iwed+iside*24][ipad]=globalave;
		bgg = globalave/cfnew[iwed + iside*24][ipad];
		fprintf(fpp," %d  %.4f .0001 \n",id0,bgg);
	}
	fclose(fpp);

	std::string::size_type dotInd = sNewDataFile.find_last_of(".");
	std::string substr = sNewDataFile.substr(0,dotInd);
	std::string rootFile = substr + "_tracks.root";
	assert (rootFile.length() > 0 && "Err :: Root filename not specified!");

	TFile *File = new TFile(rootFile.c_str(),"RECREATE","tracks");
	for(int i=0; i<48; i++) {
		for(int p=0; p<54; p++) {
			fHist.cp2pad[i][p]->Write();
		}
	}
	fHist.cp2Weight->Write(); 
	fHist.west->Write(); fHist.east->Write();
	fHist.cp2hit->Write();  fHist.cp2ave->Write();
	fHist.ymoncp2->Write(); fHist.rancp2->Write();
	fHist.cp2tra->Write(); fHist.zero3d->Write();
	fHist.zerocp2->Write();

	File->Close(); delete File;
	std::cout << "Log files created :: " << sNewDataFile << " & " << rootFile << std::endl;

	MakeRatioLERhists();
	HVStatus();

	
	std::cout << "====================================" << std::endl;
	std::cout << "Events Processed = " << iEvtsProc[0] << std::endl;
	std::cout << "Events Passed    = " << iEvtsProc[1] << std::endl;
	std::cout << "Trigger used     = " << GetTrigger() << std::endl;
	cout<<"global average = "<<globalave<<" totchan "<<totchan<<endl;
	return 0;
}
/////////////////////////////////////////////////////
// Spitout the HV zero run/sec to std::cout so if the
// HV log file is not created, you can still get this
// from job log file.
////////////////////////////////////////////////////

void cp2::HVStatus()
{
	std::cout << "HV OFF PERIODS" << std::endl;
	std::cout << "RUN\tSEC" << std::endl;
	myIt = foundHVoff.begin();
	while (myIt != foundHVoff.end())
	{
		std::cout << myIt->first << "\t" << myIt->second << std::endl;
		myIt++;
	}
}

/////////////////////////////////////////////////////
//Makes the final two plots, ratio and pad hits.
//You need to apply some styling to get them in nice
//greeny color.
////////////////////////////////////////////////////
void cp2::MakeRatioLERhists()
{

	std::cout << "ATTEMPTING TO MAKE FINAL HISTS.." <<std::endl;
	std::cout << "First and Last DB Entries provided : " << fLERfirst << ", " << fLERlastDB << std::endl;
	//read in the data corresponding to last DB entry

	Float_t fOldLers[48][54]; 
	float bgg,error;
	int err,id0,ipad,iwed,iside;

	FILE *fpp=fopen(sLastDBDataFile.c_str(),"r");
	//assert(fpp != NULL && "default file did not open!");
	if (! fpp)
	{
		std::cout << __FILE__ << "::" << __LINE__ << "::ERROR! default file did not open! returning without making final hists!" << std::endl;
		return;
	}
	std::cout << "reading default file " << sLastDBDataFile  << std::endl;

	//now check if the LER corresponding to last DB entry is specified
	if (fLERlastDB<0 || fLERlastDB<100  || fLERlastDB>300)
	{
		if (fLERlastDB<0) 
		{
			std::cout << __FILE__ << "::" << __LINE__ << "::ERROR! Specified Last DB LER cannot be <0! Returning without making final hists." << std::endl; 
			return;
		}
		if (fLERlastDB<100) std::cout << __FILE__ << "::" << __LINE__ << "::ERROR! Please double check. the LER you specified is very low!" << std::endl; 
		if (fLERlastDB>300) std::cout << __FILE__ << "::" << __LINE__ << "::ERROR! Please double check. the LER you specified is very high!" << std::endl; 
	}

	for (int i = 0; i < 2592; i++) {
		err = fscanf(fpp,"%d %f %f\n", &id0, &bgg, &error);
		if (err < 0) break;
		ipad  = (id0 & 0x0000003f);
		iwed  = ((id0 >> 8) & 0x0000001f);
		iside = ((id0 >> 6) & 0x00000001);
		fOldLers[iwed+iside*24][ipad] = bgg*fLERlastDB/fLERfirst;
	}

	fclose(fpp);

	// read in data for this run

	Float_t fNewLers[48][54];
	float bgg1,error1;
	int err1;

	FILE *fpp1=fopen(sNewDataFile.c_str(),"r");
	std::cout << "reading new file " << sNewDataFile << std::endl;
	assert(fpp1 != NULL && "new file did not open!");

	for (int i=0;i<2592;i++) {
		err1 = fscanf(fpp1,"%d %f %f\n",&id0,&bgg1,&error1);
		if (err1<0) break;
		ipad  = (id0 & 0x0000003f);
		iwed  = ((id0 >> 8) & 0x0000001f);
		iside = ((id0 >> 6) & 0x00000001);
		fNewLers[iwed+iside*24][ipad] = bgg1*fLERnew/fLERfirst;
		// we are looking a shift wrt to the last database entry - sam
	}
	fclose(fpp1);


	Float_t rat;
	for (int iside=0; iside<2; iside++) {
		for (int iwed=0; iwed<24; iwed++) {
			for (int ipad=0; ipad<54; ipad++) {
				rat = fNewLers[iwed + iside * 24][ipad] / fOldLers[iwed + iside * 24][ipad];
				if (iwed + iside * 24 ==  0 && ipad == 2) continue;
				if (iwed + iside * 24 == 31 && ipad == 5) continue;
				if (iwed + iside * 24 == 39 && ipad == 3) continue;
				histLER.Ratio->Fill(rat,1);
				histLER.OldLER->Fill(fOldLers[iwed+iside*24][ipad],1);
				histLER.NewLER->Fill(fNewLers[iwed+iside*24][ipad],1);
				//cp2pad[iwed+iside*24][ipad]->Fill(rat,1);
				histLER.PadHitsWgted_2D->Fill(iwed+iside*24,ipad,rat);
				//PadHitsWgted_2D->Fill(iwed+iside*24,ipad,fNewLers[iwed+iside*24][ipad]);
				histLER.WedgeHitsWgted->Fill(iwed+iside*24,rat/54.);
				histLER.PadHitsWgted_1D->Fill(ipad,rat/48.);
				if ((iwed+iside*24)==29) histLER.PadHitsWgted_special->Fill(ipad,rat);
			}
		}
	}


	TCanvas* c2 = new TCanvas("CP2_Ratios","CP2 Ratios",0,0,900,700);
	c2->Divide(2,2); c2->cd(1);
	histLER.OldLER->Fit("gaus");
	histLER.OldLER->Draw();
	c2->cd(2);
	histLER.PadHitsWgted_2D->Draw();   //is this is correct hist for this place? compare with cvs version
	histLER.NewLER->Fit("gaus");
	histLER.NewLER->Draw();
	c2->cd(3);
	histLER.Ratio->Fit("gaus"); 
	histLER.Ratio->Draw();
	histLER.OldLER->Draw("same");  
	c2->cd(4);
	histLER.WedgeHitsWgted->Draw();
	TLine * ln0 = new TLine (-0.5, 1., 47.5, 1.0); ln0->Draw();
	//c2.Print("ratio-222957-256904-fits.eps");   /////// change output plot name here
	std::string::size_type dotpos = sNewDataFile.find(".",0);
	std::string c1name = "ratio-fits_"+ sNewDataFile.substr(0,dotpos) + ".gif";
	assert(c1name.length()>0 && "canvas name not specified!");
	//c2.Print("ratio-266163-fits.gif");   /////// change output plot name here
	//c2->Print(c1name.c_str());   /////// change output plot name here

	TCanvas* c3 = new TCanvas("CP2_Occupancy","CP2 Occupancy",0,0,900,700);
	c3->Divide(2,2); c3->cd(1);
	histLER.PadHitsWgted_1D->Draw();
	TLine * ln1 = new TLine (-0.5, 1., 53.5, 1.0); ln1->Draw();
	c3->cd(2);
	histLER.PadHitsWgted_2D->Draw("colz");
	c3->cd(3);
	histLER.PadHitsWgted_special->SetTitle("Wedge 5East");
	histLER.PadHitsWgted_special->Draw();
	//c3.Print("ratio-222957-256904-pads.eps");    ///////// change output plot name here
	std::string c3name = "ratio-pads_"+ sNewDataFile.substr(0,dotpos) + ".gif";
	//c3.Print("ratio-266163-pads.gif");    ///////// change output plot name here
	//c3->Print(c3name.c_str());    ///////// change output plot name here

	
	TFile *f = new TFile (sFinalHistFile.c_str(),"recreate","Final CP2 calibration hists");
	if (f->IsZombie())
	{
		std::cout << "Root file "<< sFinalHistFile << " cannot be created! Final hists are not written!" << std::endl;
	}
	f->ls();
	c2->Write();
	c3->Write();
	histLER.sumet_large->Write();
	histLER.sumet_small->Write();
	histLER.met_large->Write();
	histLER.met_small->Write();
	histLER.ntracks_large->Write();
	histLER.ntracks_small->Write();
	histLER.avghitsperpad->Write();
	histLER.sumetVsNtrks->Write();
	histLER.metVsNtrks->Write();
	histLER.avghitsperpadVssumet->Write();
	histLER.avghitsperpadVsNtrks->Write();
	histLER.avghitsperpadNtrksRatio->Write();
	histLER.avghitsperpadSumetRatio->Write();
	histLER.SumetNtrksRatio->Write();

	histLER.avghitsperpadVsNVtx->Write();
	histLER.avghitsperpadVsNVtx12->Write();
	
	std::cout << "Final Hists were written to " <<  GetFinalHistOutputFileName() << std::endl;
	f->Close();


}




//_____________________________________________________________________________
void cp2::TrackExtr(float  radius, float  curvature, float  Phi0, 
		float  cotan, float  D_O, float  Z_0, float* xyz,
		int &ieta, int &iphi, bool debug) {
	/*
		Calculates the x,y,z position of a charged track passing through the
		COT, at a given radius from the origin.
		radius    is the radius from x=0,y=0 to determine the x,y,z
		curvature is the half curvature for the charged track
		Phi0      is the phi at the origin for the charged track
		cotan     is the  cot(theta) for the charged track
		D_O       is the impact parameter for the charged track
		Z_O       is the z position at the distance of closest approach
	 */
	// radii of calorimeters
	//  float const COIL_BFIELD_ENDS = 149.1;
	//  float const CEM_RMAX = 207.26;
	//  float const CHA_RMAX = 345.09;
	//  float const PEM_RMIN = 172.72;

	// boundaries of towers
	const float BOUND[11]={0.0, 24.16, 48.32, 72.48, 96.64,120.8,144.97,
		169.12,193.28,217.45,245.96};


	float epsilon = 1.0;    
	float B   = curvature* sqrt((radius*radius-D_O*D_O)/(1.0+2.0*curvature*D_O));
	float U_O = cos(Phi0);
	float V_O = sin(Phi0);
	float rho = 2.0*curvature;
	float XO  = -D_O*V_O;
	float YO  = D_O*U_O;                     

	ieta = 0;
	iphi = 0;

	if (B*B < 1.0) {
		xyz[0] = XO + (U_O*2.0*epsilon*B*sqrt(1.0-B*B)-V_O*2.0*B*B)/rho;
		xyz[1] = YO + (V_O*2.0*epsilon*B*sqrt(1.0-B*B)+U_O*2.0*B*B)/rho;

		// in principle I have to take the B-field into account only up to coil and
		// do a straight line extrapolation afterwards
		xyz[2] = Z_0 + ((cotan)*asin(B))/curvature;

		float exphi=atan2(xyz[1],xyz[0]);

		if (exphi<0) {
			exphi=exphi+2*M_PI;
		}

		float absz=float(abs(xyz[2]));

		// do central and wall for now
		if (absz<BOUND[10]) {
			iphi=int(exphi/(2*M_PI) * 24);
			for (int i=0;i<10;i++) {
				if ( absz > BOUND[i] && absz <= BOUND[i+1] ) {
					if (xyz[2]>0) {
						ieta=26+i;
					} else {
						ieta=25-i;
					}
					break;
				}
			}
		} else {
			// set them to fixed values of the next tower
			iphi=int(exphi/(2*M_PI) * 24);
			if (xyz[2]>0) {
				ieta=37;
			} else {
				ieta=14;
			}
		}

	} else {
		xyz[0] = 0.0;
		xyz[1] = 0.0;
		xyz[2] = 0.0;
	}
}

