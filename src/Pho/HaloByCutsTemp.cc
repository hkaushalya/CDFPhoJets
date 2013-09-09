///////////////////////////////////////////////////////////
//originally designed to study halos id by time and      //
// topological cuts. but then realized, better to do them//
// seperately. so there are several similar modules.     //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>
#include "samantha/Pho/HaloByCutsTemp.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>

ClassImp(HaloByCutsTemp)

//_____________________________________________________________________________
HaloByCutsTemp::HaloByCutsTemp(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bRunPermit(true)
{
	std::cout << "Hello I am HaloByCutsTemp module" << std::endl;
}

//_____________________________________________________________________________
HaloByCutsTemp::~HaloByCutsTemp() {
}

//_____________________________________________________________________________
void HaloByCutsTemp::SaveHistograms() {
}

//_____________________________________________________________________________
void HaloByCutsTemp::BookHistograms()
{
	DeleteHistograms();
  	//char name [200];
  	//char title[200];
	TFolder* new_folder = GetHistFolder(this, "1Jet","1 Jet Case");
	
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		bRunPermit = false;
		exit (0);
		return;
	} else {
		// assume the first check is enough for null older etc..
		Histograms HistoMan;
		
		TFolder *sub_folder = new_folder->AddFolder("Scene0","Scene0");
		HistoMan.GetHaloHistograms(sub_folder,halo1js0,"Halo Scenario 0, 1 Jet");
		TFolder *subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js0pw023,subsub_folder,"Halo Scenario 0, 1 Jet, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js0pw122,subsub_folder,"Halo Scenario 0, 1 Jet, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene1","Scene1");
		HistoMan.GetHaloHistograms(sub_folder,halo1js1,"Halo Scenario 1, 1 Jet");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js1pw023,subsub_folder,"Halo Scenario 1, 1 Jet, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js1pw122,subsub_folder,"Halo Scenario 1, 1 Jet, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene2","Scene2");
		HistoMan.GetHaloHistograms(sub_folder,halo1js2,"Halo Scenario 2, 1 Jet");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js2pw023,subsub_folder,"Halo Scenario 2, 1 Jet, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js2pw122,subsub_folder,"Halo Scenario 2, 1 Jet, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene3","Scene3");
		HistoMan.GetHaloHistograms(sub_folder,halo1js3,"Halo Scenario 3, 1 Jet");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js3pw023,subsub_folder,"Halo Scenario 3, 1 Jet, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js3pw122,subsub_folder,"Halo Scenario 3, 1 Jet, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene4","Scene4");
		HistoMan.GetHaloHistograms(sub_folder,halo1js4,"Halo Scenario 4, 1 Jet");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js4pw023,subsub_folder,"Halo Scenario 4, 1 Jet, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js4pw122,subsub_folder,"Halo Scenario 4, 1 Jet, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene5","Scene5");
		HistoMan.GetHaloHistograms(sub_folder,halo1js5,"Halo Scenario 5, 1 Jet");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js5pw023,subsub_folder,"Halo Scenario 5, 1 Jet, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js5pw122,subsub_folder,"Halo Scenario 5, 1 Jet, #Phi wedge 1-22");
	
		//2 Jets case 
		new_folder = GetHistFolder(this, "2Jet","2 Jet Case");

		sub_folder = new_folder->AddFolder("Scene0","Scene0");
		HistoMan.GetHaloHistograms(sub_folder,halo2js0,"Halo Scenario 0, 2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js0pw023,subsub_folder,"Halo Scenario 0, 2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js0pw122,subsub_folder,"Halo Scenario 0, 2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene1","Scene1");
		HistoMan.GetHaloHistograms(sub_folder,halo2js1,"Halo Scenario 1, 2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js1pw023,subsub_folder,"Halo Scenario 1, 2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js1pw122,subsub_folder,"Halo Scenario 1, 2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene2","Scene2");
		HistoMan.GetHaloHistograms(sub_folder,halo2js2,"Halo Scenario 2, 2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js2pw023,subsub_folder,"Halo Scenario 2, 2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js2pw122,subsub_folder,"Halo Scenario 2, 2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene3","Scene3");
		HistoMan.GetHaloHistograms(sub_folder,halo2js3,"Halo Scenario 3, 2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js3pw023,subsub_folder,"Halo Scenario 3, 2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js3pw122,subsub_folder,"Halo Scenario 3, 2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene4","Scene4");
		HistoMan.GetHaloHistograms(sub_folder,halo2js4,"Halo Scenario 4, 2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js4pw023,subsub_folder,"Halo Scenario 4, 2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js4pw122,subsub_folder,"Halo Scenario 4, 2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene5","Scene5");
		HistoMan.GetHaloHistograms(sub_folder,halo2js5,"Halo Scenario 5, 2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js5pw023,subsub_folder,"Halo Scenario 5, 2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js5pw122,subsub_folder,"Halo Scenario 5, 2 Jets, #Phi wedge 1-22");
	
	}
	
	BookEvtHistograms(aHaloByCut,"HaloByCut","Beam Halos By Id Cuts");


}


//_____________________________________________________________________________
int HaloByCutsTemp::BeginJob()
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
		StdOut(__FILE__,__LINE__,2,"Halo tag module not found.");
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
int HaloByCutsTemp::BeginRun()
{

	//int currRun =  fHeaderBlock->RunNumber();
	//std::cout << " BEGINING RUN# " << currRun << std::endl;
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int HaloByCutsTemp::Event(int ientry)
{
	SetPassed(1);
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,2,"Run Permit Failed. One or more dependencies not found.");
		return 0;
	}


	counter.evtsRunOver++;
  
  	int NsuperPho = initSpMod->GetSuperPhoSize();
	// number of jet requirement and pass accordingly. so I'll not check for that.
	
	for (int i=0; i < NsuperPho; i++) {
		if (! initSpMod->GetSuperPhoton(i)->IsTightPhoton()) continue;
		if (! initSpMod->GetSuperPhoton(i)->IsBeamHalo()) continue;		// select BH
		if (fabs(initSpMod->GetSuperPhoton(i)->GetEmTime()) > 4.8) continue;		//cut on em time

		counter.haloEvtsFound++;
		FillHistograms(aHaloByCut, initSpMod->GetSuperPhoton(i));

		int iHaloId =  initSpMod->GetSuperPhoton(i)->GetBeamHaloId();
		int iPhiSeedIndex = initSpMod->GetSuperPhoton(i)->GetPhoton()->PhiSeedIndex();
		
		//now look at each diff scenario with 1jet and 2jets
		if (jetMod->GetMyNjet(0,3) >= 1) {
			if ((iHaloId>>0) & 0x1) {
				//this  FREE function is in TagBeamHalo. should move it to FreeFunctions.hh/cc
				FillHaloHistograms(halo1js0,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo1js0pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo1js0pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			
			if ((iHaloId>>1) & 0x1) {
				FillHaloHistograms(halo1js1,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo1js1pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo1js1pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			if ((iHaloId>>2) & 0x1) {
				FillHaloHistograms(halo1js2,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo1js2pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo1js2pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			if ((iHaloId>>3) & 0x1) {
				FillHaloHistograms(halo1js3,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo1js3pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo1js3pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			if ((iHaloId>>4) & 0x1) {
				FillHaloHistograms(halo1js4,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo1js4pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo1js4pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			if ((iHaloId>>5) & 0x1) {
				FillHaloHistograms(halo1js5,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo1js5pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo1js5pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			
		
		} // end 1 jet case
		
		if (jetMod->GetMyNjet(0,3) >= 2) {
			if ((iHaloId>>0) & 0x1) {
				FillHaloHistograms(halo2js0,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo2js0pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo2js0pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			
			if ((iHaloId>>1) & 0x1) {
				FillHaloHistograms(halo2js1,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo2js1pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo2js1pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			if ((iHaloId>>2) & 0x1) {
				FillHaloHistograms(halo2js2,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo2js2pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo2js2pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			if ((iHaloId>>3) & 0x1) {
				FillHaloHistograms(halo2js3,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo2js3pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo2js3pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			if ((iHaloId>>4) & 0x1) {
				FillHaloHistograms(halo2js4,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo2js4pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo2js4pw122, initSpMod->GetSuperPhoton(i));
				}
			}
			if ((iHaloId>>5) & 0x1) {
				FillHaloHistograms(halo2js5,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo2js5pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo2js5pw122, initSpMod->GetSuperPhoton(i));
				}
			}
		
		}// di-jet case:
	
		//break; // don't care about rest of the photons
				// and to avoid mutiple entries to met/sumet plots
				//1-05-08.. 2nd thoughts, if we are looking to figure out a good eff. cuts
				// lets loook at all possible BH. there will not be not many events with 2 BHs
				
	}

	if (GetPassed()) counter.evtsPassModule++;
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void HaloByCutsTemp::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void HaloByCutsTemp::FillDataBlocks(int ientry)
{
}
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void HaloByCutsTemp::FillPhiWedgeHistogram(PhiWedgeHist_t& Hist, SuperPhoton* sp)
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
/*-------------------------------------------------------------------*/
void HaloByCutsTemp::BookPhiWedgeHistograms(PhiWedgeHist_t& Hist,
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



/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void HaloByCutsTemp::FillHistograms(EventHist_t& Hist, SuperPhoton* sp)
{
		Hist.PhoEmTime->Fill(sp->GetEmTime());
		Hist.PhoEmTimeVsRun->Fill(fHeaderBlock->RunNumber(), sp->GetEmTime());
		Hist.Met->Fill(jetMod->GetMyMetCorr(0,3));
		Hist.Sumet->Fill(jetMod->GetMySumEtCorr(0,3));
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void HaloByCutsTemp::BookEvtHistograms(EventHist_t& Hist,
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

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int HaloByCutsTemp::EndJob() {
	std::string modname = GetName();
	std::string msg;
	msg  = "[HAT:00:]----- end job: ---- " + modname + "\n";
	msg += "[HAT:01:] Events Run Over ------------ = " + ToStr(counter.evtsRunOver) + "\n";
	msg += "[HAT:02:] Events Pass this module ---- = " + ToStr(counter.evtsPassModule) + "\n";
	msg += "[HAT:03:] Halo Events Found ---------- = " + ToStr(counter.haloEvtsFound) + "\n";
	msg += "---------------------------------------------------\n";
	std::cout << msg;
	return 0;
}
