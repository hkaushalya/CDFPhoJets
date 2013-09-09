///////////////////////////////////////////////////////////
// this module will tag the tight photons in the super   //
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
#include "samantha/Pho/TagTightPhotons.hh"
#include <cmath>

ClassImp(TagTightPhotons)

//_____________________________________________________________________________
TagTightPhotons::TagTightPhotons(const char* name, const char* title):
  TStnModule(name,title),
  bRunPermit(true),
  iNTightBits(12),
  bNoSummary(false),
  bDebug(false)
{
	//std::cout << "Hello I am TagTightPhotons module" << std::endl;
}

//_____________________________________________________________________________
TagTightPhotons::~TagTightPhotons() {
}

//_____________________________________________________________________________
void TagTightPhotons::SaveHistograms() {
}

//_____________________________________________________________________________
void TagTightPhotons::BookHistograms()
{
	DeleteHistograms();
	BookPhotonHistograms(hTightPho,"TightPhotons","Tight Photons' Variables Histograms(only tight photons");
	BookIdCutsHistograms(hIdCut,"IdCutsCumm","Id cuts cummulative Histograms of all photons");
}


//_____________________________________________________________________________
int TagTightPhotons::BeginJob()
{
	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initSpMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
		bRunPermit = false;
		return 0;
	}
				// register the data blocks

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.oneCEMtightphoton 		= 0;
	counter.onePLUGtightphoton 		= 0;
	counter.twoCEMtightphotons 	= 0;
	counter.twoPLUGtightphotons 	= 0;
	counter.morethan2CEMtightphotons = 0;
	counter.morethan2PLUGtightphotons = 0;
	counter.PhoEvts				= 0;
	return 0;
}

//_____________________________________________________________________________
int TagTightPhotons::BeginRun()
{
	return 0;
} //BeginRun

//_____________________________________________________________________________
int TagTightPhotons::Event(int ientry)
{
	if (bDebug)
	{
		std::cout << "\n\n************************** TagTightPhoton:: ";
		GetHeaderBlock()->Print();
	}
	SetPassed(1);
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"One or more dependencies not found. pl check. Run Permit not cleared!");
		exit (1);
	}
	counter.evtsRunOver++;
   
	int iNSp = 0;
	int iCEMtight = 0, iPLUGtight = 0;
	
	for (int i = 0; i < initSpMod->GetSuperPhoSize() ; i++) {
		SuperPhoton* aSp = initSpMod->GetSuperPhoton(i); // does this poniter becomes dangling afterward??
		int iIdWord = -1;

		if (aSp->GetDetector() == 0)
		{
			iIdWord = GetCenTightIdWord(aSp);
			initSpMod->GetSuperPhoton(i)->SetTightPhotonId(iIdWord);
			if (iIdWord==0) ++iCEMtight;
		} else if (aSp->GetDetector() == 1)
		{
			iIdWord = GetPlugTightIdWord(aSp);
			initSpMod->GetSuperPhoton(i)->SetTightPhotonId(iIdWord);
			if (iIdWord==0) ++iPLUGtight;
		}

		if (iIdWord == 0) {
			//FillPhotonsHist(hTightPho,aSp);  //fill tight photon variables histos
			iNSp++;
		}
		//FillIdCutHist(iIdWord);
	}

	if (iNSp>0) counter.PhoEvts++;

	if (iCEMtight==1) counter.oneCEMtightphoton++;
	else if (iCEMtight==2) counter.twoCEMtightphotons++;
	else if (iCEMtight>=3) counter.morethan2CEMtightphotons++;

	if (iPLUGtight==1) counter.onePLUGtightphoton++;
	else if (iPLUGtight==2) counter.twoPLUGtightphotons++;
	else if (iPLUGtight>=3) counter.morethan2PLUGtightphotons++;

	if (GetPassed()) {
		counter.evtsPassModule++;
	}

	return 0;

} // Event



/*-------------------------------------------------------------------*/
// Tight Photon ID CUTS
/*-------------------------------------------------------------------*/
int TagTightPhotons::GetCenTightIdWord(SuperPhoton* sp)
{
	if (!sp) {
		StdOut(__FILE__,__LINE__,3,"SuperPhoton is NULL. pl. check. exiting");
		exit (1);
	}

	unsigned int IDWord = 0x0;

	// change this. no need to copy stuff anymore!
	int   detector 	= sp->GetDetector();
	float ECorr 		= sp->GetECorr();
	float EtCorr 		= sp->GetEtCorr();
	float XCes			= sp->GetXCes();
	float ZCes 			= sp->GetZCes();
	float HadEm 		= sp->GetHadEm();
	float Chi2Mean 	= sp->GetChi2Mean();  			// (strip+wire)/2
	int   N3d			= sp->GetN3d();
	float IsoEtCorr	= sp->GetIso4();			//Corrected for leakage and MI
	float TrkPt  		= sp->GetTrkPt();     				//using max pt track in the cluster
	float TrkIso 		= sp->GetTrkIso();
	float CesWireE2	= sp->GetCesWireE2();
	float CesStripE2	= sp->GetCesStripE2();

	if (detector != 0						) {
		IDWord |= kCentralBit;
		if (bDebug) std::cout << "\tdetector failed "<< std::endl;
	}
	if (EtCorr < 7.0 						) {
		IDWord |= kEtCorrBit;
		if (bDebug) std::cout << "\tEtc failed  "<< std::endl;
	}
	if (fabs(XCes) > 21					) {
		IDWord |= kXCesBit;
		if (bDebug) std::cout << "\tXCes failed " << std::endl;
	}
	if (fabs(ZCes) < 9 ||  fabs(ZCes) > 230) {
		IDWord |= kZCesBit;
		if (bDebug) std::cout << "\tZCes failed" << std::endl;
	}
	if ( !((HadEm < 0.125) || (HadEm < (0.055 + 0.00045 * ECorr)))) {
		IDWord |= kHadEmBit;
		if (bDebug) std::cout << "\tHadEm failed" << std::endl;
	}

	if (EtCorr < 20) {
		if (IsoEtCorr > 0.1*EtCorr) {
			IDWord |= kIso4Bit;
			if (bDebug) std::cout << "\tIso failed " << std::endl;
		}
	} else {
		if (IsoEtCorr > (2.0+0.02*(EtCorr-20.0)) ) {
			IDWord |= kIso4Bit;
			if (bDebug) std::cout << "\tIso failed " << std::endl;
		}
	}

	if (Chi2Mean > 20						) {
		IDWord |= kChi2MeanBit;
		if (bDebug) std::cout << "\tChi2Mean failed" << std::endl;
	}
	if (N3d > 1								) {
		IDWord |= kN3dBit;
		if (bDebug) std::cout << "\tN3d failed" << std::endl;
	}
	if (N3d == 1) {
		if (TrkPt > (1+0.005*EtCorr)		) {
			IDWord |= kTrkPtBit;
			if (bDebug) std::cout << "\tTrkPt failed" << std::endl;
		}
	}
	if (TrkIso > (2.0+0.005*EtCorr)	) {
		IDWord |= kTrkIsoBit;
		if (bDebug) std::cout << "\tTrkIso failed"<< std::endl;
	}

	if (EtCorr<18) {
		if ( (CesWireE2 * sp->GetSinTheta()) > (0.14*EtCorr)) {
			IDWord |= kCesWireE2Bit;
			if (bDebug) std::cout << "\tCES wire failed"  << std::endl;
		}
		if ( (CesStripE2 * sp->GetSinTheta()) > (0.14*EtCorr)) {
			IDWord |= kCesStripE2Bit;
			if (bDebug) std::cout << "\tCES wire failed"  << std::endl;
		}
	} else {
		if ( (CesWireE2 * sp->GetSinTheta()) > (2.4+0.01*EtCorr)) {
			IDWord |= kCesWireE2Bit;
			if (bDebug) std::cout << "\tCES wire failed"  << std::endl;
		}
		if ( (CesStripE2 * sp->GetSinTheta()) > (2.4+0.01*EtCorr)) {
			IDWord |= kCesStripE2Bit;
			if (bDebug) std::cout << "\tCES wire failed"  << std::endl;
		}
	}

	if (bDebug) std::cout << "Tight Central Photon IdWord = " << IDWord << std::endl;

	return IDWord;

} // TightPhotonCuts



/*-------------------------------------------------------------------*/
// Std Plug Tight Photon ID CUTS
/*-------------------------------------------------------------------*/
int TagTightPhotons::GetPlugTightIdWord(SuperPhoton* sp)
{
	if (!sp) {
		StdOut(__FILE__,__LINE__,3,"SuperPhoton is NULL. pl. check. exiting");
		exit (1);
	}

	unsigned int IDWord = 0x0;

	// change this. no need to copy stuff anymore!
	int   detector 	= sp->GetDetector();
	float ECorr 		= sp->GetECorr();
	float EtCorr 		= sp->GetEtCorr();
	float Pes5by9U		= sp->GetPhoton()->Prof5by9U();
	float Pes5by9V		= sp->GetPhoton()->Prof5by9V();
	float HadEm 		= sp->GetHadEm();
	float PemChi2 		= sp->GetPhoton()->Chi3x3();  			// PEm Shape CHI2
	float IsoEtCorr	= sp->GetIso4();			//Corrected for leakage and MI
	float TrkIso 		= sp->GetTrkIso();

	if (detector != 1						) {
		IDWord |= kPlugDetBit;
		if (bDebug) std::cout << "tightPlug::detector failed:: det/DetEta"<< detector << "/ " << sp->GetDetEta() << std::endl;
	}

	if (EtCorr < 7.0 						) {
		IDWord |= kPlugEtCorrBit;
	}
	
	if (ECorr <= 100)
	{
		if (HadEm > 0.05)
		{
			IDWord |= kPlugHadEmBit;
		}
	} else {
		if (HadEm > (0.05 + 0.026 * log(ECorr/100)) )		//in C++ 'log' returns natural log. 'log10' is for base10.
		{
			IDWord |= kPlugHadEmBit;
		}
	}


	if (EtCorr < 20) {
		if (IsoEtCorr > 0.1*EtCorr) {
			IDWord |= kPlugIso4Bit;
		}
	} else {
		if (IsoEtCorr > (2.0+0.02*(EtCorr-20.0)) ) {
			IDWord |= kPlugIso4Bit;
		}
	}
	
	if (PemChi2 > 10						) {
		IDWord |= kPlugChi2Bit;
	}

	if (fabs(Pes5by9U) < 0.65					) {
		IDWord |= kPlugPes5by9UBit;
	}

	if (fabs(Pes5by9V) < 0.65					) {
		IDWord |= kPlugPes5by9VBit;
	}

	if (TrkIso > (2.0+0.005*EtCorr)	) {
		IDWord |= kPlugTrkIsoBit;
	}

	if (bDebug) std::cout << "Tight Plug Photon IdWord = " << IDWord << std::endl;
	
	return IDWord;

}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TagTightPhotons::FillIdCutHist(int& iIdWord)
{
	bool bStop = false;
	if (iIdWord != 0) {
		for (int i=1; i <= iNTightBits; i++) {		//first bin is not filled for now.
			int iBit = (iIdWord >> (i-1)) & 0x1;
			if (iBit) {
				hIdCut.Sum->Fill(i);			// fill failed cuts buy all
				bStop = true;
			} else {
				if (!bStop) hIdCut.Cummulative->Fill(i);   // stop filling at the first failed cut
			}
		}
	} else {
		for (int i=1; i <= iNTightBits; i++) hIdCut.Cummulative->Fill(i);
	}
}


/*-------------------------------------------------------------------*/
//
/*-------------------------------------------------------------------*/
void TagTightPhotons::FillPhotonsHist(PhotonPlots_t& hist, SuperPhoton* sp)
{
	hist.Detector->Fill(sp->GetDetector());
	hist.DetEta->Fill(sp->GetDetEta());
	hist.DetPhi->Fill(sp->GetDetPhi());
	hist.EtCorr->Fill(sp->GetEtCorr());
	hist.XCes->Fill(sp->GetXCes());
	hist.ZCes->Fill(sp->GetZCes());
	hist.HadEm->Fill(sp->GetHadEm());
	hist.Chi2Mean->Fill(sp->GetChi2Mean());
	hist.N3d->Fill(sp->GetN3d());
	hist.Iso4->Fill(sp->GetIso4());
	hist.TrkPt->Fill(sp->GetTrkPt());
	hist.TrkIso->Fill(sp->GetTrkIso());
	hist.Ces2Wire->Fill(sp->GetCesWireE2());
	hist.Ces2Strip->Fill(sp->GetCesStripE2());
	hist.HadTDCtime->Fill(sp->GetHadTdcTime());
}
//_____________________________________________________________________________
void TagTightPhotons::BookIdCutsHistograms(IdCutsPlots_t& hist,
							std::string sFoldName, std::string sFoldTitle)
{
	hist.Cummulative = new TH1F("TightCutsCumm","Cummulative Plot of passing Tight Photon Id Cuts",16,-0.5,15.5);
	hist.Cummulative->SetXTitle("Objects in a bin has passed all cuts to the left of that bin.");
	hist.Cummulative->GetXaxis()->SetBinLabel(1,"All Photons");
	hist.Cummulative->GetXaxis()->SetBinLabel(2,"Central");
	hist.Cummulative->GetXaxis()->SetBinLabel(3,"Etc>7");
	hist.Cummulative->GetXaxis()->SetBinLabel(4,"XCes");
	hist.Cummulative->GetXaxis()->SetBinLabel(5,"ZCes");
	hist.Cummulative->GetXaxis()->SetBinLabel(6,"HadEm");
	hist.Cummulative->GetXaxis()->SetBinLabel(7,"IsoEtCorr");
	hist.Cummulative->GetXaxis()->SetBinLabel(8,"Chi2Mean");
	hist.Cummulative->GetXaxis()->SetBinLabel(9,"N3d");
	hist.Cummulative->GetXaxis()->SetBinLabel(10,"TrkPt");
	hist.Cummulative->GetXaxis()->SetBinLabel(11,"TrkIso");
	hist.Cummulative->GetXaxis()->SetBinLabel(12,"2ndCesWire");
	hist.Cummulative->GetXaxis()->SetBinLabel(13,"2ndCesStrip");

	hist.Sum = new TH1F("TightCutsSum","Plot of #gamma failed Tight Photon Cuts",16,-0.5,15.5);
	hist.Cummulative->SetXTitle("all the cuts failed by an object.");
	hist.Sum->GetXaxis()->SetBinLabel(1,"All Photons");
	hist.Sum->GetXaxis()->SetBinLabel(2,"Central");
	hist.Sum->GetXaxis()->SetBinLabel(3,"Etc>7");
	hist.Sum->GetXaxis()->SetBinLabel(4,"XCes");
	hist.Sum->GetXaxis()->SetBinLabel(5,"ZCes");
	hist.Sum->GetXaxis()->SetBinLabel(6,"HadEm");
	hist.Sum->GetXaxis()->SetBinLabel(7,"IsoEtCorr");
	hist.Sum->GetXaxis()->SetBinLabel(8,"Chi2Mean");
	hist.Sum->GetXaxis()->SetBinLabel(9,"N3d");
	hist.Sum->GetXaxis()->SetBinLabel(10,"TrkPt");
	hist.Sum->GetXaxis()->SetBinLabel(11,"TrkIso");
	hist.Sum->GetXaxis()->SetBinLabel(12,"2ndCesWire");
	hist.Sum->GetXaxis()->SetBinLabel(13,"2ndCesStrip");
	
	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	
	new_folder->Add(hist.Cummulative);
	new_folder->Add(hist.Sum);

}
//_____________________________________________________________________________
void TagTightPhotons::BookPhotonHistograms(PhotonPlots_t& hist,
							std::string sFoldName, std::string sFoldTitle)
{
	hist.Detector 	= new TH1F("Detector","Photon Detector Region",3,0,3);
	hist.DetEta 	= new TH1F("DetEta","Photon Detector Eta",100,-5,5);
	hist.DetPhi 	= new TH1F("DetPhi","Photon Detector Phi",105,-3.5,7);
	hist.EtCorr 	= new TH1F("EtCorr","Photon Corrected Et",600,0,120);
	hist.XCes 		= new TH1F("XCes","Photon XCes",640,-32,32);
	hist.ZCes 		= new TH1F("ZCes","Photon ZCes",2000,-250,250);
	hist.HadEm 		= new TH1F("HadEm","Photon HadEm",50,0,0.5);
	hist.Chi2Mean 	= new TH1F("Chi2Mean","Photon Chi2Mean (Wire+Strip/2)",800,0,80);
	hist.N3d 		= new TH1F("N3d","Photon N3d",10,0,10);
	hist.Iso4 		= new TH1F("IsoEtCorr","Photon IsoCorr",1500,-5,10);
	hist.TrkPt 		= new TH1F("Trkpt","Photon TrkPt",1000,0,100);
	hist.TrkIso 	= new TH1F("TrkIso","Photon TrkIso",150,0,15);
	hist.Ces2Wire 	= new TH1F("Ces2Wire","Photon CES(2nd) Wire",400,0,40);
	hist.Ces2Strip = new TH1F("Ces2Strip","Photon CES(2nd) Strip",400,0,40);
	hist.HadTDCtime = new TH1F("HadTDCtime","Photon Had TDC time",200,-50,150);

	//use GetBinWidth() in here!!!! u moran
	char ytitle[30];
	sprintf(ytitle,"Events/%.2f",hist.DetEta->GetBinWidth(1));
	hist.DetEta->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.DetPhi->GetBinWidth(1));
	hist.DetPhi->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.EtCorr->GetBinWidth(1));
	hist.EtCorr->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.XCes->GetBinWidth(1));
	hist.XCes->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.ZCes->GetBinWidth(1));
	hist.ZCes->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.HadEm->GetBinWidth(1));
	hist.HadEm->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.Chi2Mean->GetBinWidth(1));
	hist.Chi2Mean->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.N3d->GetBinWidth(1));
	hist.N3d->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.Iso4->GetBinWidth(1));
	hist.Iso4->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.TrkPt->GetBinWidth(1));
	hist.TrkPt->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.TrkIso->GetBinWidth(1));
	hist.TrkIso->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.Ces2Wire->GetBinWidth(1));
	hist.Ces2Wire->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.Ces2Strip->GetBinWidth(1));
	hist.Ces2Strip->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.HadTDCtime->GetBinWidth(1));
	hist.HadTDCtime->SetYTitle(ytitle);

	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	
	new_folder->Add(hist.Detector);
	new_folder->Add(hist.DetEta);
	new_folder->Add(hist.DetPhi);
	new_folder->Add(hist.EtCorr);
	new_folder->Add(hist.XCes);
	new_folder->Add(hist.ZCes);
	new_folder->Add(hist.HadEm);
	new_folder->Add(hist.Chi2Mean);
	new_folder->Add(hist.N3d);
	new_folder->Add(hist.Iso4);
	new_folder->Add(hist.TrkPt);
	new_folder->Add(hist.TrkIso);
	new_folder->Add(hist.Ces2Wire);
	new_folder->Add(hist.Ces2Strip);
	new_folder->Add(hist.HadTDCtime);

} //BookPhotonHistograms


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TagTightPhotons::EndJob()
{
	if (GetSummaryStat()) return 0;	
	std::string sMsg, sMsg1;
	std::string sModName = GetName();

	sMsg += "[TTP:00:]----- end job: ----- " + sModName + "\n";
	sMsg += "[TTP:01:] Events Processed------------ = " + ToStr(counter.evtsRunOver) + "\n";
	sMsg += "[TTP:02:] Events Passed -------------- = " + ToStr(counter.evtsPassModule) + "\n";
	sMsg += "[TTP:03:] Evt w. >=1 T. pho (cen+plug) = " + ToStr(counter.PhoEvts) + "\n";
	sMsg += "[TTP:04:] Evt w. 1   T. pho (cen/plug) = " + ToStr(counter.oneCEMtightphoton) +"/"+ ToStr(counter.onePLUGtightphoton) + "\n";
	sMsg += "[TTP:05:] Evt w. 2   T. pho (cen/plug) = " + ToStr(counter.twoCEMtightphotons) +"/"+ ToStr(counter.twoPLUGtightphotons) +"\n";
	sMsg += "[TTP:06:] Evt w. >2  T. pho (cen/plug) = " + ToStr(counter.morethan2CEMtightphotons) +"/"+ToStr(counter.morethan2PLUGtightphotons)+ "\n";
	sMsg += "---------------------------------------------------\n";
	std::cout << sMsg;

	return 0;
}
