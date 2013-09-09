#include "TBenchmark.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TIterator.h"
#include <iostream>
#include <sstream>
#include <string>
#include "Stntuple/Stntuple/loop/TStnAna.hh"
#include "Stntuple/Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/Stntuple/loop/TStnOutputModule.hh"
#include "Stntuple/Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/base/base/TStnDataset.hh"
#include "Stntuple/ana/ana/TEmTimingModule.hh"

#include "samantha/samantha/Pho/InitSuperPhotons.hh"
#include "samantha/samantha/Pho/TagTightPhotons.hh"
#include "samantha/samantha/Pho/TagLoosePhotons.hh"
#include "samantha/samantha/Pho/TagElectrons.hh"
#include "samantha/samantha/Pho/TagLooseElectrons.hh"
#include "samantha/samantha/Pho/TagBeamHalo.hh"
#include "samantha/samantha/Pho/TriggerModule.hh"
//#include "samantha/samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/samantha/Pho/TMyJetFilterModule.hh"
#include "samantha/samantha/Pho/TagPMTSpikes.hh"
#include "samantha/samantha/Pho/SetEmTimes.hh"
#include "samantha/samantha/utils/PassFailTree.hh"
#include "samantha/samantha/utils/BadPointerChecker.hh"
#include "samantha/samantha/Pho/TagCosmicPhotons.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/FlatStupleMaker.hh"
#include "samantha/samantha/Pho/TagConvElectrons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"
#include "samantha/samantha/Pho/CosmicStudy.hh"
#include "samantha/samantha/Pho/EventProperties.hh"

void StudyCosmics(int evts = 100)
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnAna* ap;
	TChain *chain = 0;

		/* TO TEST USE A LOCAL FILE */ 
		chain = new TChain("STNTUPLE");
		//chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1.root");
		//chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2.root");
		//chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_0.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_1.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_2.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_3.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_4.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_5.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_6.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_7.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_8.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_9.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_10.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_11.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_12.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_13.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_14.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_15.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_16.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_17.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_18.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_19.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_20.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_21.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_22.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_23.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_24.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_25.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_26.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_27.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_28.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_29.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_30.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_31.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_32.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_33.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_34.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_35.root");
		
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_0.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_1.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_2.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_3.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_4.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_5.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_6.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_7.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_8.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_9.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_10.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_11.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_12.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_13.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_14.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_15.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_16.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_17.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_18.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_19.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_20.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_21.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_22.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_23.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_24.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_25.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_26.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_27.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_28.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_29.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_30.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_31.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_32.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_33.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_34.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_2_35.root");

		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_0.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_1.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_2.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_3.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_4.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_5.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_6.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_7.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_8.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_9.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_10.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_11.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_12.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_13.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_14.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_15.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_16.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_17.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_18.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_19.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_20.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_21.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_22.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_23.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_24.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_25.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_26.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_27.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_28.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_29.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_30.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_31.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_32.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_33.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_34.root");
		chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_3_35.root");

		
		std::cout << chain->GetEntries() << std::endl;

		ap = new TStnAna(chain);

	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

  /******************************************************/
  // define all module here
  /******************************************************/
  
  //PassFailTree *passFailFile = new PassFailTree;
//  BadPointerChecker *ptChecker = new BadPointerChecker;
  
	TriggerModule* trigMod;
	trigMod = new TriggerModule;
		trigMod->SetGoodRunListFile("goodrun_v19_pho_00.txt");  // to p13
		trigMod->SetGoodRunBit(1);		//1= use good run, 0=do not use it
		trigMod->SetTriggerBit(1);		// 1= require tigger, 0=NO
		trigMod->SetMinVtx(1);
		trigMod->SetMaxVtx(1000);
		trigMod->SetUseVtxZcut(1);		//1= require z<60 0=not cut on z
	
	InitSuperPhotons* myinit = new InitSuperPhotons;

	TagTightPhotons* tagTight = new TagTightPhotons;
	TagLoosePhotons* tagLoose = new TagLoosePhotons;
	
	TagElectrons* tagTightEle = new TagElectrons;
	TagLooseElectrons* tagLooseEle = new TagLooseElectrons;
											//1,2 - pass only if a electron found and remove others from list

	PhoEleFilter *peFilter = new PhoEleFilter;

																									  
	//these are common parameters for central/JES UP/JES DOWN jet collections
	int MinNjet15 = 1;
	int JTC_imode = 1;		// 0==MC, 1==data
	int UnclParamSwitch = 0;
	float MaxDeltaPhiJMet = 0.3;
	int NpointsToGenerate = 0;
	int SelectMetEvent = 0;
	int RemoveDuplicate = 0;
	int RemoveBadMet = 0;
	int DumpEvent = 0;

	TMyJetFilterModule *m5 = new TMyJetFilterModule();
		m5->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		m5->SetJTC_imode(JTC_imode); 			// 0==MC; 1==data
		m5->SetJTC_systcode(0);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		m5->SetUnclParamSwitch(UnclParamSwitch); 	// 0==photon parametrization; 1==Zee parametrization 
		m5->SetMaxDeltaPhiJMet(MaxDeltaPhiJMet); 	// controls dPhi in MyMetCleanUpCut
		m5->SetNpointsToGenerate(NpointsToGenerate); 	// number of pseudo-experiments per each event
		m5->SetSelectMetEvent(SelectMetEvent); 		// 0==do nothing; 1==select event with large MET
		m5->SetRemoveDuplicate(RemoveDuplicate); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		m5->SetRemoveBadMet(RemoveBadMet); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		///m5->SetUseMetPDFscenario(1);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		m5->SetDumpEvent(DumpEvent); 			// 0==do nothing; 1==dump event
		m5->SetDumpEventFileName("DATA_dumpEvent.dat");
		m5->SetLargeMetEventFileName("DATA_largeMetEvent.dat");
		m5->SetDatFileName("DATA_jet.dat");
		

	TagBeamHalo* tagBH = new TagBeamHalo();

	TagPMTSpikes* tagPMT = new TagPMTSpikes;
		tagPMT->SetMode(1);	//1= pass only if no spikes found
		
	TagPhoenixPhotons* phoenix = new TagPhoenixPhotons;


	TEmTimingModule* EmT = new TEmTimingModule();
	SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
		setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt");
		
	EventProperties *evt = new EventProperties;
		
	CosmicStudy* cosmic = new CosmicStudy;
	
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(peFilter,1);
	ap->AddModule(m5,1);
	ap->AddModule(tagPMT,1);
	ap->AddModule(tagBH,1);
	ap->AddModule(EmT,1);
	ap->AddModule(setEmTime,1);

	ap->AddModule(evt,1);
	ap->AddModule(cosmic,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	ap->Run(evts);
  	
	
	std::string filename = "CosmicStudy.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	//std::cout << "::::::::: Created Pass/Fail file => "<< passFailFile->GetFileName() << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
