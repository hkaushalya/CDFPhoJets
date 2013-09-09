///////////////////////////////////////////////////////////
// This module is to calculate Beam Halo mis-id rate     //
// using electrons. So this module must be identical to  //
// HaloByCutsTemp_Pho.cc/hh except select electrons      //
// instead of photons. --  01-09-2008                    //
// use this to get beam halo mis-id rates                //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#include "samantha/Pho/HaloByCutsTemp_Ele.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>

ClassImp(HaloByCutsTemp_Ele)

//_____________________________________________________________________________
HaloByCutsTemp_Ele::HaloByCutsTemp_Ele(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bRunPermit(true)
{
	std::cout << "Hello I am HaloByCutsTemp_Ele module" << std::endl;
}

//_____________________________________________________________________________
HaloByCutsTemp_Ele::~HaloByCutsTemp_Ele() {
}

//_____________________________________________________________________________
void HaloByCutsTemp_Ele::SaveHistograms() {
}

//_____________________________________________________________________________
void HaloByCutsTemp_Ele::BookHistograms()
{
	DeleteHistograms();
  	//char name [200];
  	//char title[200];
	TFolder* new_folder = GetHistFolder(this, "1Jet"," >=1 Jets Case");
	
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		bRunPermit = false;
		return;
	} else {
		// assume the first check is enough for null older etc..
		Histograms HistoMan;
		
		TFolder *sub_folder = new_folder->AddFolder("Scene0","Scene0");
		HistoMan.GetHaloHistograms(sub_folder,halo1js0,"Halo Scenario 0,  >=1 Jets");
		TFolder *subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js0pw023,subsub_folder,"Halo Scenario 0,  >=1 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js0pw122,subsub_folder,"Halo Scenario 0,  >=1 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene1","Scene1");
		HistoMan.GetHaloHistograms(sub_folder,halo1js1,"Halo Scenario 1,  >=1 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js1pw023,subsub_folder,"Halo Scenario 1,  >=1 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js1pw122,subsub_folder,"Halo Scenario 1,  >=1 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene2","Scene2");
		HistoMan.GetHaloHistograms(sub_folder,halo1js2,"Halo Scenario 2,  >=1 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js2pw023,subsub_folder,"Halo Scenario 2,  >=1 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js2pw122,subsub_folder,"Halo Scenario 2,  >=1 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene3","Scene3");
		HistoMan.GetHaloHistograms(sub_folder,halo1js3,"Halo Scenario 3,  >=1 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js3pw023,subsub_folder,"Halo Scenario 3,  >=1 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js3pw122,subsub_folder,"Halo Scenario 3,  >=1 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene4","Scene4");
		HistoMan.GetHaloHistograms(sub_folder,halo1js4,"Halo Scenario 4,  >=1 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js4pw023,subsub_folder,"Halo Scenario 4,  >=1 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js4pw122,subsub_folder,"Halo Scenario 4,  >=1 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene5","Scene5");
		HistoMan.GetHaloHistograms(sub_folder,halo1js5,"Halo Scenario 5,  >=1 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo1js5pw023,subsub_folder,"Halo Scenario 5,  >=1 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo1js5pw122,subsub_folder,"Halo Scenario 5,  >=1 Jets, #Phi wedge 1-22");
	
		// >=2 Jets case 
		new_folder = GetHistFolder(this, "2Jet","2 Jet Case");

		sub_folder = new_folder->AddFolder("Scene0","Scene0");
		HistoMan.GetHaloHistograms(sub_folder,halo2js0,"Halo Scenario 0, >=2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js0pw023,subsub_folder,"Halo Scenario 0, >=2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js0pw122,subsub_folder,"Halo Scenario 0, >=2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene1","Scene1");
		HistoMan.GetHaloHistograms(sub_folder,halo2js1,"Halo Scenario 1, >=2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js1pw023,subsub_folder,"Halo Scenario 1, >=2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js1pw122,subsub_folder,"Halo Scenario 1, >=2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene2","Scene2");
		HistoMan.GetHaloHistograms(sub_folder,halo2js2,"Halo Scenario 2, >=2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js2pw023,subsub_folder,"Halo Scenario 2, >=2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js2pw122,subsub_folder,"Halo Scenario 2, >=2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene3","Scene3");
		HistoMan.GetHaloHistograms(sub_folder,halo2js3,"Halo Scenario 3, >=2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js3pw023,subsub_folder,"Halo Scenario 3, >=2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js3pw122,subsub_folder,"Halo Scenario 3, >=2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene4","Scene4");
		HistoMan.GetHaloHistograms(sub_folder,halo2js4,"Halo Scenario 4, >=2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js4pw023,subsub_folder,"Halo Scenario 4, >=2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js4pw122,subsub_folder,"Halo Scenario 4, >=2 Jets, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene5","Scene5");
		HistoMan.GetHaloHistograms(sub_folder,halo2js5,"Halo Scenario 5, >=2 Jets");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23","PhiWedge0_23");
		BookPhiWedgeHistograms(halo2js5pw023,subsub_folder,"Halo Scenario 5, >=2 Jets, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22","PhiWedge1_22");
		BookPhiWedgeHistograms(halo2js5pw122,subsub_folder,"Halo Scenario 5, >=2 Jets, #Phi wedge 1-22");


		// same as above, but for MEt>25Gev Events
		new_folder = GetHistFolder(this, "1Jet"," >=1 Jets Case");

		sub_folder = new_folder->AddFolder("Scene0Met25","Scene0Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo1js0Met25,"Halo Scenario 0, >=1 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo1js0pw023Met25,subsub_folder,"Halo Scenario 0,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo1js0pw122Met25,subsub_folder,"Halo Scenario 0,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene1Met25","Scene1Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo1js1Met25,"Halo Scenario 1,  >=1 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo1js1pw023Met25,subsub_folder,"Halo Scenario 1,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo1js1pw122Met25,subsub_folder,"Halo Scenario 1,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene2Met25","Scene2Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo1js2Met25,"Halo Scenario 2,  >=1 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo1js2pw023Met25,subsub_folder,"Halo Scenario 2,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo1js2pw122Met25,subsub_folder,"Halo Scenario 2,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene3Met25","Scene3Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo1js3Met25,"Halo Scenario 3,  >=1 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo1js3pw023Met25,subsub_folder,"Halo Scenario 3,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo1js3pw122Met25,subsub_folder,"Halo Scenario 3,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene4Met25","Scene4Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo1js4Met25,"Halo Scenario 4,  >=1 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo1js4pw023Met25,subsub_folder,"Halo Scenario 4,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo1js4pw122Met25,subsub_folder,"Halo Scenario 4,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene5Met25","Scene5Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo1js5Met25,"Halo Scenario 5,  >=1 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo1js5pw023Met25,subsub_folder,"Halo Scenario 5,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo1js5pw122Met25,subsub_folder,"Halo Scenario 5,  >=1 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		//>=2 Jets case 
		new_folder = GetHistFolder(this, "2Jet","2 Jet Case");

		sub_folder = new_folder->AddFolder("Scene0Met25","Scene0Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo2js0Met25,"Halo Scenario 0, >=2 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo2js0pw023Met25,subsub_folder,"Halo Scenario 0, >=2 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo2js0pw122Met25,subsub_folder,"Halo Scenario 0, >=2 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene1Met25","Scene1Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo2js1Met25,"Halo Scenario 1, >=2 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo2js1pw023Met25,subsub_folder,"Halo Scenario 1, >=2 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo2js1pw122Met25,subsub_folder,"Halo Scenario 1, >=2 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene2Met25","Scene2Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo2js2Met25,"Halo Scenario 2, >=2 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo2js2pw023Met25,subsub_folder,"Halo Scenario 2, >=2 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo2js2pw122Met25,subsub_folder,"Halo Scenario 2, >=2 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene3Met25","Scene3Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo2js3Met25,"Halo Scenario 3, >=2 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo2js3pw023Met25,subsub_folder,"Halo Scenario 3, >=2 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo2js3pw122Met25,subsub_folder,"Halo Scenario 3, >=2 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene4Met25","Scene4Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo2js4Met25,"Halo Scenario 4, >=2 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo2js4pw023Met25,subsub_folder,"Halo Scenario 4, >=2 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo2js4pw122Met25,subsub_folder,"Halo Scenario 4, >=2 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	
		sub_folder = new_folder->AddFolder("Scene5Met25","Scene5Met25");
		HistoMan.GetHaloHistograms(sub_folder,halo2js5Met25,"Halo Scenario 5, >=2 Jets, ME_{T}>25GeV");
		subsub_folder = sub_folder->AddFolder("PhiWedge0_23Met25","PhiWedge0_23Met25");
		BookPhiWedgeHistograms(halo2js5pw023Met25,subsub_folder,"Halo Scenario 5, >=2 Jets, ME_{T}>25GeV, #Phi wedge 0,23");
		subsub_folder = sub_folder->AddFolder("PhiWedge1_22Met25","PhiWedge1_22Met25");
		BookPhiWedgeHistograms(halo2js5pw122Met25,subsub_folder,"Halo Scenario 5, >=2 Jets, ME_{T}>25GeV, #Phi wedge 1-22");
	

		new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
		vCounterHistLabels.push_back("Events Processed");
		vCounterHistLabels.push_back("Events Passed");
		vCounterHistLabels.push_back("1 Jet Events");
		vCounterHistLabels.push_back("2 Jets Events");
		HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);

	}
	
	BookEvtHistograms(aHaloByCut,"HaloByCut","Beam Halos By Id Cuts");

	

}


//_____________________________________________________________________________
int HaloByCutsTemp_Ele::BeginJob()
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

  	eleMod = (TagElectrons*) ((TStnAna*) GetAna()->GetModule("TagElectrons"));
	if (!eleMod) {
		StdOut(__FILE__,__LINE__,2,"Electron tag module not found.");
		bRunPermit = false;
	}
  	haloMod = (TagBeamHalo*) ((TStnAna*) GetAna()->GetModule("TagBeamHalo"));
	if (!haloMod) {
		StdOut(__FILE__,__LINE__,2,"Halo tag module not found.");
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
int HaloByCutsTemp_Ele::BeginRun()
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
/*-------------------------------------------------------------------*/
int HaloByCutsTemp_Ele::Event(int ientry)
{
	// Event  loop

	SetPassed(1);
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,2,"Run Permit Failed. One or more dependencies or requirements not met.");
		exit (1);
		return 0;
	}

	hCount.iCounter->Fill(0);
	counter.evtsRunOver++;
  
  	int NsuperPho = initSpMod->GetSuperPhoSize();
	float fMetCor = jetMod->GetMyMetCorr(0,3);		//cone0.4, Et>15GeV
	// number of jet requirement and pass accordingly. so I'll not check for that.
	
	for (int i=0; i < NsuperPho; i++) {
		if (! initSpMod->GetSuperPhoton(i)->IsTightElectron()) continue;
		if (! initSpMod->GetSuperPhoton(i)->IsBeamHalo()) continue;
		if (fabs(initSpMod->GetSuperPhoton(i)->GetEmTime()) > 4.8) continue;		//cut on em time
		hCount.iCounter->Fill(4);
		counter.haloEvtsFound++;
		FillHistograms(aHaloByCut, initSpMod->GetSuperPhoton(i));

		int iHaloId =  initSpMod->GetSuperPhoton(i)->GetBeamHaloId();
		int iPhiSeedIndex = initSpMod->GetSuperPhoton(i)->GetPhoton()->PhiSeedIndex();
		
		//now look at each diff scenario with 1jet and 2jets
		if (jetMod->GetMyNjet(0,3) >= 1) {
			hCount.iCounter->Fill(2);

			if ((iHaloId>>0) & 0x1) {
				//this  FREE function is in TagBeamHalo. should move it to FreeFunctions.hh/cc
				FillHaloHistograms(halo1js0,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo1js0pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo1js0pw122, initSpMod->GetSuperPhoton(i));
				}

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo1js0Met25,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo1js0pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo1js0pw122Met25, initSpMod->GetSuperPhoton(i));
					}
				}//met>25
				
			} //halo scenario=0
			
			
			if ((iHaloId>>1) & 0x1) {
				FillHaloHistograms(halo1js1,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo1js1pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo1js1pw122, initSpMod->GetSuperPhoton(i));
				}

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo1js1Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo1js1pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo1js1pw122Met25, initSpMod->GetSuperPhoton(i));
					}
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

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo1js2Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo1js2pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo1js2pw122Met25, initSpMod->GetSuperPhoton(i));
					}
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

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo1js3Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo1js3pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo1js3pw122Met25, initSpMod->GetSuperPhoton(i));
					}
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

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo1js4Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo1js4pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo1js4pw122Met25, initSpMod->GetSuperPhoton(i));
					}
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

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo1js5Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo1js5pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo1js5pw122Met25, initSpMod->GetSuperPhoton(i));
					}
				}
			}
			
		
		} // end  >=1 Jets case
		
		if (jetMod->GetMyNjet(0,3) >= 2) {
			hCount.iCounter->Fill(3);
			
			if ((iHaloId>>0) & 0x1) {
				FillHaloHistograms(halo2js0,haloMod->GetHaloStuff(i),
										initSpMod->GetSuperPhoton(i)->GetPhoton());

				if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
					FillPhiWedgeHistogram(halo2js0pw023, initSpMod->GetSuperPhoton(i));
				} else {
					FillPhiWedgeHistogram(halo2js0pw122, initSpMod->GetSuperPhoton(i));
				}

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo2js0Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo2js0pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo2js0pw122Met25, initSpMod->GetSuperPhoton(i));
					}
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

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo2js1Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo2js1pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo2js1pw122Met25, initSpMod->GetSuperPhoton(i));
					}
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

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo2js2Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo2js2pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo2js2pw122Met25, initSpMod->GetSuperPhoton(i));
					}
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

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo2js3Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo2js3pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo2js3pw122Met25, initSpMod->GetSuperPhoton(i));
					}
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

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo2js4Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo2js4pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo2js4pw122Met25, initSpMod->GetSuperPhoton(i));
					}
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

				//for Met>25Gev
				if (fMetCor >25) {
					FillHaloHistograms(halo2js5Met25,haloMod->GetHaloStuff(i),
											initSpMod->GetSuperPhoton(i)->GetPhoton());

					if (iPhiSeedIndex == 0 || iPhiSeedIndex == 23 ) {
						FillPhiWedgeHistogram(halo2js5pw023Met25, initSpMod->GetSuperPhoton(i));
					} else {
						FillPhiWedgeHistogram(halo2js5pw122Met25, initSpMod->GetSuperPhoton(i));
					}
				}
			}
		
		}// di-jet case:
	
		break; // don't care about rest of the electrons, if any.
				// all i need is ele+1jet or ele+2jet signature. this is all fine as long as.. read more..
				// thing is that the denominator and the numerator has to counted the same way.
				// for eg. if more than one electron per event is allwoed in the denom.,  count the numerator
				// the sameway. simply put, if we allow multiple entries from one event,
				// make sure to count both deno. and numerator the same way!
				
	}

	if (GetPassed()) {
		counter.evtsPassModule++;
		hCount.iCounter->Fill(1);
	}
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void HaloByCutsTemp_Ele::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void HaloByCutsTemp_Ele::FillDataBlocks(int ientry)
{
}
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void HaloByCutsTemp_Ele::FillPhiWedgeHistogram(PhiWedgeHist_t& Hist, SuperPhoton* sp)
{
	if (!sp) {
		StdErr(__FILE__,__LINE__,3,"SuperPhoton Pointer is NULL. No histograms are filled.");
		return;
	}

	float fMetCor = jetMod->GetMyMetCorr(0,3);		//cone0.4, Et>15GeV
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
void HaloByCutsTemp_Ele::BookPhiWedgeHistograms(PhiWedgeHist_t& Hist,
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
void HaloByCutsTemp_Ele::FillHistograms(EventHist_t& Hist, SuperPhoton* sp)
{
		Hist.PhoEmTime->Fill(sp->GetEmTime());
		Hist.PhoEmTimeVsRun->Fill(fHeaderBlock->RunNumber(), sp->GetEmTime());
		Hist.Met->Fill(jetMod->GetMyMetCorr(0,3));
		Hist.Sumet->Fill(jetMod->GetMySumEtCorr(0,3));
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void HaloByCutsTemp_Ele::BookEvtHistograms(EventHist_t& Hist,
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
int HaloByCutsTemp_Ele::EndJob() {
	std::string modname = GetName();
	std::string msg;
	msg  = "[HAT:00:]----- end job: ---- " + modname + "\n";
	msg += "[HAT:01:] Events Run Over ------------ = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	msg += "[HAT:02:] Events Pass this module ---- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	msg += "[HAT:03:] Halo+ 1 Jet Events --------- = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
	msg += "[HAT:04:] Halo+ 2 Jet Events --------- = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n";
	msg += "[HAT:05:] Halo Events Found ---------- = " + ToStr(hCount.iCounter->GetBinContent(5)) + "\n";
	msg += "---------------------------------------------------\n";
	std::cout << msg;
	return 0;
}
