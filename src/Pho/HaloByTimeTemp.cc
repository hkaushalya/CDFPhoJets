///////////////////////////////////////////////////////////
// this is to study Beam Halo using EM timing.           //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "samantha/Pho/HaloByTimeTemp.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>
#include "TVector2.h"
ClassImp(HaloByTimeTemp)

//_____________________________________________________________________________
HaloByTimeTemp::HaloByTimeTemp(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bRunPermit(true),
  iHaloScenario(0)
{
	std::cout << "Hello I am HaloByTimeTemp module" << std::endl;
}

//_____________________________________________________________________________
HaloByTimeTemp::~HaloByTimeTemp() {
}

//_____________________________________________________________________________
void HaloByTimeTemp::SaveHistograms() {
}

//_____________________________________________________________________________
void HaloByTimeTemp::BookHistograms()
{
	DeleteHistograms();
  	//char name [200];
  	//char title[200];
	TFolder* new_folder = GetHistFolder(this, "NotHalo","Not Beam Halos");
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		bRunPermit = false;
		return;
	} else {
	
		Histograms HistoMan;
		HistoMan.GetHaloHistograms(new_folder,notHalo,"Not Beam Halo");
		
		//i'll assume the first check for a 'null' folder is enough. if it fails, others will too.
		// so will not check for null folder
		// for each time bin, there will be two subfolders for 0 & 23 phi wedge plots and 1 to 22 phi
		// wedge plots

		new_folder = GetHistFolder(this, "Halo13","Beam Halos - Em Time >= -13ns & <-12ns");
		HistoMan.GetHaloHistograms(new_folder,halo13,"Halos EM Time >-13ns & <-12ns");
		TFolder *sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo13pw023,sub_folder,"Halos EM Time >= -13ns & <-12ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo13pw122,sub_folder,"Halos EM Time >= -13ns & <-12ns,#Phi wedge 1-22");
		
		new_folder = GetHistFolder(this, "Halo12","Beam Halos - Em Time >= -12ns & <-11ns");
		HistoMan.GetHaloHistograms(new_folder,halo12,"Halos EM Time >= -12ns & <-11ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo12pw023,sub_folder,"Halos EM Time >= -12ns & <-11ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo12pw122,sub_folder,"Halos EM Time >= -12ns & <-11ns,#Phi wedge 1-22");
		
		new_folder = GetHistFolder(this, "Halo11","Beam Halos - Em Time >= -11ns & <-10ns");
		HistoMan.GetHaloHistograms(new_folder,halo11,"Halos EM Time >= -11ns & <-10ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo11pw023,sub_folder,"Halos EM Time >= -11ns & <-10ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo11pw122,sub_folder,"Halos EM Time >= -11ns & <-10ns,#Phi wedge 1-22");
		
		new_folder = GetHistFolder(this, "Halo10","Beam Halos - Em Time >= -10ns & <-9ns");
		HistoMan.GetHaloHistograms(new_folder,halo10,"Halos EM Time >= -10ns & <-9ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo10pw023,sub_folder,"Halos EM Time >= -10ns & <-9ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo10pw122,sub_folder,"Halos EM Time >= -10ns & <-9ns,#Phi wedge 1-22");
		
		new_folder = GetHistFolder(this, "Halo9","Beam Halos - Em Time >= -9ns & <-8ns");
		HistoMan.GetHaloHistograms(new_folder,halo9,"Halos EM Time >= -9ns & <-8ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo9pw023,sub_folder,"Halos EM Time >= -9ns & <-8ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo9pw122,sub_folder,"Halos EM Time >= -9ns & <-8ns,#Phi wedge 1-22");
		
		new_folder = GetHistFolder(this, "Halo8","Beam Halos - Em Time >= -8ns & <-7ns");
		HistoMan.GetHaloHistograms(new_folder,halo8,"Halos EM Time >= -8ns & <-7ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo8pw023,sub_folder,"Halos EM Time >= -8ns & <-7ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo8pw122,sub_folder,"Halos EM Time >= -8ns & <-7ns,#Phi wedge 1-22");
		
		new_folder = GetHistFolder(this, "Halo7","Beam Halos - Em Time >= -7ns & <-6ns");
		HistoMan.GetHaloHistograms(new_folder,halo7,"Halos EM Time >= -7ns & <-6ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo7pw023,sub_folder,"Halos EM Time >= -7ns & <-6ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo7pw122,sub_folder,"Halos EM Time >= -7ns & <-6ns,#Phi wedge 1-22");
		
		new_folder = GetHistFolder(this, "Halo6","Beam Halos - Em Time >= -6ns & <-5ns");
		HistoMan.GetHaloHistograms(new_folder,halo6,"Halos EM Time >= -6ns & <-5ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo6pw023,sub_folder,"Halos EM Time >= -6ns & <-5ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo6pw122,sub_folder,"Halos EM Time >= -6ns & <-5ns,#Phi wedge 1-22");
		
		new_folder = GetHistFolder(this, "Halo5","Beam Halos - Em Time >= -5ns & <-4ns");
		HistoMan.GetHaloHistograms(new_folder,halo5,"Halos EM Time >= -5ns & <-4ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo5pw023,sub_folder,"Halos EM Time >= -5ns & <-4ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo5pw122,sub_folder,"Halos EM Time >= -5ns & <-4ns,#Phi wedge 1-22");
		
		new_folder = GetHistFolder(this, "Halo4","Beam Halos - Em Time >= -4ns & <-3ns");
		HistoMan.GetHaloHistograms(new_folder,halo4,"Halos EM Time >= -4ns & <-3ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo4pw023,sub_folder,"Halos EM Time >= -4ns & <-3ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo4pw122,sub_folder,"Halos EM Time >= -4ns & <-3ns,#Phi wedge 1-22");

		new_folder = GetHistFolder(this, "Halo13_6","Beam Halos - Em Time >= -13ns & <-6ns");
		HistoMan.GetHaloHistograms(new_folder,halo13_6,"Halos EM Time >= -13ns & <-6ns");
		sub_folder = new_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo136pw023,sub_folder,"Halos EM Time >= -13ns & <-6ns,#Phi wedge 0,23");
		sub_folder = new_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo136pw122,sub_folder,"Halos EM Time >= -13ns & <-6ns,#Phi wedge 1-22");
	
		new_folder =0;
		sub_folder =0;

	}
	

	//BookEvtHistograms(aHaloByTime,"HaloByTime","Beam Halo: EM Time >-12ns & < -6ns ");



}


//_____________________________________________________________________________
int HaloByTimeTemp::BeginJob()
{
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}

  	jetMod = (JetFilterModuleV2*) ((TStnAna*) GetAna()->GetModule("JetFilterV2"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		bRunPermit = false;
	}

  	tightMod = (TagTightPhotons*) ((TStnAna*) GetAna()->GetModule("TagTightPhotons"));
	if (!tightMod) {
		StdOut(__FILE__,__LINE__,2,"Tight photon tag module not found.");
		bRunPermit = false;
	}

  	haloMod = (TagBeamHalo*) ((TStnAna*) GetAna()->GetModule("TagBeamHalo"));
	if (!haloMod) {
		StdOut(__FILE__,__LINE__,2,"Beam Halo tag module not found.");
		bRunPermit = false;
	}
  	setTimeMod = (SetEmTimes*) ((TStnAna*) GetAna()->GetModule("SetEmTimes"));
	if (!setTimeMod) {
		StdOut(__FILE__,__LINE__,2,"SetEmTimes module not found.");
		bRunPermit = false;
	}
	
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.haloEvtsFound 		= 0;
  

	return 0;
}

//_____________________________________________________________________________
int HaloByTimeTemp::BeginRun()
{
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int HaloByTimeTemp::Event(int ientry)
{
	SetPassed(1);
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,2,"Run Permit Failed. One or more dependencies not found.");
		exit (1);
		return 0;
	}
	counter.evtsRunOver++;
  
	//now select tight photon+2jets. assume jet filter mod will take care of 
	// number of jet requirement and pass accordingly. so I'll not check for that.
	
	int Run = fHeaderBlock->RunNumber();
	if (Run<190851) return 0; // exclude first 400pb-1 without em timing
		
  	int NsuperPho = initSpMod->GetSuperPhoSize();

	for (int i=0; i < NsuperPho; i++) {
		if (! initSpMod->GetSuperPhoton(i)->IsTightPhoton()) continue;
	
//01-04-08 now i don't remember why i wanted to pick things based on 
//topological cuts. so comment these out for now.
		//look at only for a specific halo scenario (0-5)
		//int iFoundHaloScene = (initSpMod->GetSuperPhoton(i)->GetBeamHaloId() >> GetHaloScenario()) & 0x1;
		//if (! iFoundHaloScene) continue;

		float fEmTime = initSpMod->GetSuperPhoton(i)->GetEmTime();
		int iPhiSeedIndex = initSpMod->GetSuperPhoton(i)->GetPhoton()->PhiSeedIndex();
		if (fEmTime < -13. ||  fEmTime  >-6.) continue;
		
		//this  FREE function is in TagBeamHalo. should move it to FreeFunctions.hh/cc
		FillHaloHistograms(halo13_6,haloMod->GetHaloStuff(i),
								initSpMod->GetSuperPhoton(i)->GetPhoton());

		if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
			FillPhiWedgeHistogram(halo136pw023, initSpMod->GetSuperPhoton(i));
		} else {
			FillPhiWedgeHistogram(halo136pw122, initSpMod->GetSuperPhoton(i));
		}

		//now look at each ns to see which is more efficient	
		if (fEmTime >= -13.&& fEmTime  <-12.) {
			FillHaloHistograms(halo13,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo13pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo13pw122, initSpMod->GetSuperPhoton(i));
			}
		} else if (fEmTime >= -12.&& fEmTime  <-11.) {
			FillHaloHistograms(halo12,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo12pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo12pw122, initSpMod->GetSuperPhoton(i));
			}
		} else if (fEmTime >= -11.&& fEmTime  <-10.) {
			FillHaloHistograms(halo11,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo11pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo11pw122, initSpMod->GetSuperPhoton(i));
			}
		} else if (fEmTime >= -10.&& fEmTime  <-9.) {
			FillHaloHistograms(halo10,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo10pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo10pw122, initSpMod->GetSuperPhoton(i));
			}
		} else if (fEmTime >= -9.&& fEmTime  <-8.) {
			FillHaloHistograms(halo9,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo9pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo9pw122, initSpMod->GetSuperPhoton(i));
			}
		} else if (fEmTime >= -8.&& fEmTime  <-7.) {
			FillHaloHistograms(halo8,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo8pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo8pw122, initSpMod->GetSuperPhoton(i));
			}
		} else if (fEmTime >= -7.&& fEmTime  <-6.) {
			FillHaloHistograms(halo7,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo7pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo7pw122, initSpMod->GetSuperPhoton(i));
			}
		} else if (fEmTime >= -6.&& fEmTime  <-5.) {
			FillHaloHistograms(halo6,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo6pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo6pw122, initSpMod->GetSuperPhoton(i));
			}
		} else if (fEmTime >= -5.&& fEmTime  <-4.) {
			FillHaloHistograms(halo5,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo5pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo5pw122, initSpMod->GetSuperPhoton(i));
			}
		} else if (fEmTime >= -4.&& fEmTime  <-3.) {
			FillHaloHistograms(halo4,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());
			if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
				FillPhiWedgeHistogram(halo4pw023, initSpMod->GetSuperPhoton(i));
			} else {
				FillPhiWedgeHistogram(halo4pw122, initSpMod->GetSuperPhoton(i));
			}
		}

			
		counter.haloEvtsFound++;

		//break; // don't care about rest of the photons
				//if not i'll be putting multiple entries into met.sumet plots
				//01-05-08 let it be. this mod is to study halo.		
	}

	if (GetPassed()) {
		counter.evtsPassModule++;
	}
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void HaloByTimeTemp::Cleanup()
{
	// clean and clear things for the run
}


/*-------------------------------------------------------------------*/
void HaloByTimeTemp::FillDataBlocks(int ientry)
{
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void HaloByTimeTemp::FillPhiWedgeHistogram(PhiWedgeHist_t& Hist, SuperPhoton* sp)
{
	if (!sp) {
		StdErr(__FILE__,__LINE__,3,"SuperPhoton Pointer is NULL. No histograms are filled.");
		return;
	}

	float fMetCor = jetMod->GetMyMetCorr(0,3);
	Hist.Met->Fill(fMetCor);
	if (fMetCor!=0) Hist.PhoEtMetRatio->Fill(fMetCor/sp->GetEtCorr());

	float fPhoPhi = TVector2::Phi_0_2pi(sp->GetCorVec().Phi());	// TLorentzVector return Phi from -pi to +pi
	float fMetPhi = jetMod->GetMyMetPhiCorr(0,3); 			// return Phi from 0 to 2pi
	Hist.PhoMetDelPhi->Fill(fabs(DelPhi(fPhoPhi,fMetPhi)));	// call from FreeFunctions.hh
	Hist.PhoEmTime->Fill(sp->GetEmTime());
	Hist.MetPhi->Fill(fMetPhi);
	Hist.PhoPhi->Fill(fPhoPhi);
}



/*-------------------------------------------------------------------*/
//	no longer used!
/*-------------------------------------------------------------------*/
void HaloByTimeTemp::FillHistogram(EventHist_t& Hist, SuperPhoton* sp)
{
		Hist.PhoEmTime->Fill(sp->GetEmTime());
		Hist.PhoEmTimeVsRun->Fill(fHeaderBlock->RunNumber(), sp->GetEmTime());
		Hist.Met->Fill(jetMod->GetMyMetCorr(0,3));
		Hist.Sumet->Fill(jetMod->GetMySumEtCorr(0,3));
}

/*-------------------------------------------------------------------*/
//	no longer used!
/*-------------------------------------------------------------------*/
void HaloByTimeTemp::BookEvtHistograms(EventHist_t& Hist,
										std::string sFoldName, std::string text2add) 
{
	
	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldName);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
  char name [200];
  char title[200];

  sprintf(name,"PhoTime");
  sprintf(title,"%s - Em Time Time", text2add.c_str());
  Hist.PhoEmTime=new TH1F(name,title,480,-80.0,160.0);
  new_folder->Add(Hist.PhoEmTime);

  sprintf(name,"PhoTimeVsRun");
  sprintf(title,"%s - Mean EM Times of #gamma (+2Jets>15GeV) Vs Run Number", text2add.c_str());
  Hist.PhoEmTimeVsRun = new TProfile(name,title,270000,130000,400000,-80,160);
  Hist.PhoEmTimeVsRun->GetXaxis()->SetTitle("Run Number");
  Hist.PhoEmTimeVsRun->GetYaxis()->SetTitle("Mean Em Times of Photons");
  new_folder->Add(Hist.PhoEmTimeVsRun);

  sprintf(name,"EvtMet");
  sprintf(title,"%s - Event Met", text2add.c_str());
  Hist.Met = new TH1F(name,title,1200,0,1200);
  Hist.Met->GetXaxis()->SetTitle("MEt(0) (GeV)");
  Hist.Met->GetYaxis()->SetTitle("Events/ 1 GeV");
  new_folder->Add(Hist.Met);


  sprintf(name,"EvtSumEt");
  sprintf(title,"%s - Event SumEt", text2add.c_str());
  Hist.Sumet = new TH1F(name,title,1200,0,1200);
  Hist.Sumet->GetXaxis()->SetTitle("SumEt(0) (GeV)");
  Hist.Sumet->GetYaxis()->SetTitle("Events/ 1 GeV");
  new_folder->Add(Hist.Sumet);

}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void HaloByTimeTemp::BookPhiWedgeHistograms(PhiWedgeHist_t& Hist,
										TFolder* new_folder, std::string text2add) 
{
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	char name [200];
	char title[200];
	char ytitle[200];

	sprintf(name,"EvtMet");
	sprintf(title,"%s - Event Met", text2add.c_str());
	Hist.Met = new TH1F(name,title,1200,0,1200);
	Hist.Met->GetXaxis()->SetTitle("ME_{T}^{(corr)} (GeV)");
	sprintf(ytitle,"Events / %3.3d", Hist.Met->GetBinWidth(1));
	Hist.Met->GetYaxis()->SetTitle(ytitle);
	new_folder->Add(Hist.Met);


	sprintf(name,"PhoMetDelPhi");
	sprintf(title,"%s - #Delta #Phi of #gamma & ME_{T}", text2add.c_str());
	Hist.PhoMetDelPhi = new TH1F(name,title,320,-8,8);
	Hist.PhoMetDelPhi->GetXaxis()->SetTitle("#Delta #Phi");
	sprintf(ytitle,"Events / %3.3f", Hist.PhoMetDelPhi->GetBinWidth(1));
	Hist.PhoMetDelPhi->GetYaxis()->SetTitle(ytitle);
	new_folder->Add(Hist.PhoMetDelPhi);

	sprintf(name,"PhoEtMetRatio");
	sprintf(title,"%s - ME_{T}/E_{T}^{#gamma}", text2add.c_str());
	Hist.PhoEtMetRatio = new TH1F(name,title,320,-8,8);
	Hist.PhoEtMetRatio->GetXaxis()->SetTitle("ME_{T}/E_{T}^{#gamma}");
	sprintf(ytitle,"Events / %3.3f", Hist.PhoEtMetRatio->GetBinWidth(1));
	Hist.PhoEtMetRatio->GetYaxis()->SetTitle(ytitle);
	new_folder->Add(Hist.PhoEtMetRatio);

	sprintf(name,"PhoTime");
	sprintf(title,"%s - Em Time Time", text2add.c_str());
	Hist.PhoEmTime=new TH1F(name,title,480,-80.0,160.0);
	Hist.PhoEmTime->GetXaxis()->SetTitle("EM Time (ns)");
	sprintf(ytitle,"Events / %3.3f", Hist.PhoEmTime->GetBinWidth(1));
	Hist.PhoEmTime->GetYaxis()->SetTitle(ytitle);
	new_folder->Add(Hist.PhoEmTime);
	
	sprintf(name,"MetPhi");
	sprintf(title,"%s - ME_{T} #Phi", text2add.c_str());
	Hist.MetPhi = new TH1F(name,title,320,-8,8);
	Hist.MetPhi->GetXaxis()->SetTitle("#Phi");
	sprintf(ytitle,"Events / %3.3f", Hist.MetPhi->GetBinWidth(1));
	Hist.MetPhi->GetYaxis()->SetTitle(ytitle);
	new_folder->Add(Hist.MetPhi);

	sprintf(name,"PhoPhi");
	sprintf(title,"%s - ME_{T} #Phi", text2add.c_str());
	Hist.PhoPhi = new TH1F(name,title,320,-8,8);
	Hist.PhoPhi->GetXaxis()->SetTitle("#Phi");
	sprintf(ytitle,"Events / %3.3f", Hist.PhoPhi->GetBinWidth(1));
	Hist.PhoPhi->GetYaxis()->SetTitle(ytitle);
	new_folder->Add(Hist.PhoPhi);

}


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int HaloByTimeTemp::EndJob() {
	std::string modname = GetName();
	std::string msg;
	msg  = "[HTT:00:]----- end job: ----" + modname + "\n";
	msg += "[HTT:01:] Events Processed ----------- = " + ToStr(counter.evtsRunOver) + "\n";
	msg += "[HTT:02:] Events Pass this module ---- = " + ToStr(counter.evtsPassModule) + "\n";
	msg += "[HTT:03:] Halo Events Found ---------- = " + ToStr(counter.haloEvtsFound) + "\n";
	msg += "---------------------------------------------------\n";
	std::cout << msg;
	return 0;
}
