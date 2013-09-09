#include "TBenchmark.h"
#include "TSystem.h"
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
#include "samantha/samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/samantha/Pho/TMyJetFilterModule.hh"
#include "samantha/samantha/Pho/TagPMTSpikes.hh"
#include "samantha/samantha/Pho/SetEmTimes.hh"
#include "samantha/samantha/Pho/HaloByCutsTemp_Pho.hh"
#include "samantha/samantha/utils/PassFailTree.hh"
#include "samantha/samantha/utils/BadPointerChecker.hh"
#include "samantha/samantha/Pho/HaloByCutsTemp_Ele.hh"
#include "samantha/samantha/Pho/Pho2JetsTemp.hh"
#include "samantha/samantha/Pho/EventProperties.hh"
#include "samantha/samantha/Pho/TagCosmicPhotons.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/HaloByCutsTemp_Pho_Skim.hh"
#include "samantha/samantha/Pho/HaloByCutsTemp_Pho_Skim.hh"
#include "samantha/samantha/Pho/HaloStudyAllWithMetEt.hh"
#include "samantha/samantha/Pho/FlatStupleMaker.hh"


void TestEWKMCFlatNtupleMaker(unsigned runEvts=29999)
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnAna* ap;
	TChain *chain = 0;

		/* TO TEST USE A LOCAL FILE */ 
			chain = new TChain("STNTUPLE");
				printf("\n <<<<<<<<<< RUNNING ON LOCAL EWK DATA SAMPLE >>>>>>>>>>\n");
		  		//chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Wenu/we02f713.0002e0sg");	//W->enu
		  		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Zee/zd02d3a0.0005e1s6");	//Z->ee
		  		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Zee/zd02d3f4.0003e1s6");	//Z->ee
				


		ap = new TStnAna(chain);



	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

	//******************************************//
	//for stripping evts
	//----  need this part for event preselection only
   //TStnOutputModule* out= new TStnOutputModule("HaloS4_5_Skim.root");
	//out->SetMaxFileSize(300);
	//ap->SetOutputModule(out);
	
	//******************************************//		 

  /******************************************************/
  // define all module here
  /******************************************************/
  
  PassFailTree *passFailFile = new PassFailTree;
  BadPointerChecker *ptChecker = new BadPointerChecker;
  
	TriggerModule* trigMod;
	trigMod = new TriggerModule;
		trigMod->SetGoodRunListFile("goodrun_v19_pho_00.txt");  // to p13
		trigMod->SetGoodRunBit(1);		//1= use good run, 0=do not use it
		trigMod->SetTriggerBit(0);		// 1= require tigger, 0=NO
		trigMod->SetMinVtx(0);
		trigMod->SetMaxVtx(1000);
		trigMod->SetUseVtxZcut(0);		//1= require z<60 0=not cut on z
	
	InitSuperPhotons* myinit = new InitSuperPhotons;

	TagTightPhotons* tagTight = new TagTightPhotons;
	TagLoosePhotons* tagLoose = new TagLoosePhotons;
	
	TagElectrons* tagTightEle = new TagElectrons;
	TagLooseElectrons* tagLooseEle = new TagLooseElectrons;
											//1,2 - pass only if a electron found and remove others from list

/*
	JetFilterModuleV2 *jf = new JetFilterModuleV2();  //---------- My Vertex Filter Initialization
		jf->SetMinNjet15(1);				//Jets required to pass (Et>15GeV)
		jf->SetJTC_imode(1); 			// 0==MC; 1==data
		jf->SetUnclParamSwitch(0); 	// 0==photon parametrization; 1==Zee parametrization 
		jf->SetMaxDeltaPhiJMet(0.3); 	// controls dPhi in MyMetCleanUpCut
		jf->SetNpointsToGenerate(1); 	// number of pseudo-experiments per each event
		jf->SetSelectMetEvent(0); 		// 0==do nothing; 1==select event with large MET
		jf->SetRemoveDuplicate(0); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		jf->SetRemoveBadMet(0); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		///jf->SetUseMetPDFscenario(1);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		jf->SetDumpEvent(0); 			// 0==do nothing; 1==dump event
		jf->SetDumpEventFileName("V2DATA_dumpEvent.dat");
		jf->SetLargeMetEventFileName("V2DATA_largeMetEvent.dat");
		jf->SetDatFileName("V2DATA_jet.dat");
		
		//new stuff for v2
		jf->SetSampleID(0); // 0==signal
		jf->SetMinDeltaPhiJMet(-10.0); // every event should pass this cut
		//-----------------------------
		jf->SetSelectSigMetEvent(0); // =0 do nothing; =1 select events with MetSig>metsigCut
		jf->SetMetSigCut(0);
		jf->SetAnalysisMode(1);      // =0 doesn't fill AnalysisHisto, but fills MetStudyHisto
												// =1 fills AnalysisHisto, but doesn't fill MetStudyHisto
		//-----------------------------
		jf->SetSelectExoticEvent(0);
		jf->SetDoVxSwap(0); // 0=do nothing; 1=swap vertices to minimize MET
*/
																									  

	TMyJetFilterModule *m5 = new TMyJetFilterModule("");  //---------- My Vertex Filter Initialization
		m5->SetMinNjet15(1);				//Jets required to pass (Et>15GeV)
		m5->SetJTC_imode(0); 			// 0==MC; 1==data
		m5->SetJTC_systcode(6);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		m5->SetUnclParamSwitch(0); 	// 0==photon parametrization; 1==Zee parametrization 
		m5->SetMaxDeltaPhiJMet(0.3); 	// controls dPhi in MyMetCleanUpCut
		m5->SetNpointsToGenerate(1); 	// number of pseudo-experiments per each event
		m5->SetSelectMetEvent(0); 		// 0==do nothing; 1==select event with large MET
		m5->SetRemoveDuplicate(0); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		m5->SetRemoveBadMet(0); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		///m5->SetUseMetPDFscenario(1);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		m5->SetDumpEvent(0); 			// 0==do nothing; 1==dump event
		m5->SetDumpEventFileName("DATA_dumpEvent.dat");
		m5->SetLargeMetEventFileName("DATA_largeMetEvent.dat");
		m5->SetDatFileName("DATA_jet.dat");

		
	EventProperties *evtProp = new EventProperties;

	TagBeamHalo* tagBH = new TagBeamHalo();

	//TagPMTSpikes* tagPMT = new TagPMTSpikes;
		//tagPMT->SetMode(1);	//1= pass only if no spikes found
	
	FlatStupleMaker *flatS = new FlatStupleMaker;	
	
	//ap->AddModule(passFailFile,1);
//	ap->AddModule(ptChecker,1);
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(m5,1);
	//ap->AddModule(jf,1);
	//ap->AddModule(tagPMT,1);
	ap->AddModule(tagBH,1);
	ap->AddModule(flatS,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	ap->Run(runEvts);
	//ap->Run(100000);
	//ap->ProcessEvent(207078,1150507);
	//ap->ProcessEvent(206990,247577);
  	
	
	std::string filename = "EWKMCFlatStupleTest.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	std::cout << "::::::::: Created Pass/Fail file => "<< passFailFile->GetFileName() << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
