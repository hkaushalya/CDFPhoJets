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
	
#include "samantha/samantha/Pho/TriggerModule.hh"
#include "samantha/samantha/Pho/InitSuperPhotons.hh"
#include "samantha/samantha/Pho/TagTightPhotons.hh"
#include "samantha/samantha/Pho/TagLoosePhotons.hh"
#include "samantha/samantha/Pho/TagElectrons.hh"
#include "samantha/samantha/Pho/TagLooseElectrons.hh"
//#include "samantha/samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/samantha/Pho/TMyJetFilterModule.hh"
#include "samantha/samantha/Pho/DumpStntupleEvent.hh"

void DebugEmJetMisMatch(unsigned runEvts=999999999)
{
	TStnAna* ap;
	TChain *chain = 0;

	chain = new TChain("STNTUPLE");
	chain->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/JetEmMisMatch/JetEmMismatchAndJetEta.root");

	ap = new TStnAna(chain);

	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

  /******************************************************/
  // define all module here
  /******************************************************/
  
	TriggerModule* trigMod = new TriggerModule;
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
																									  

	TMyJetFilterModule *m5 = new TMyJetFilterModule();  //---------- My Vertex Filter Initialization
		m5->SetMinNjet15(1);				//Jets required to pass (Et>15GeV)
		m5->SetJTC_imode(1); 			// 0==MC; 1==data
		m5->SetJTC_systcode(0);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
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


		DumpStntupleEvent *dump = new DumpStntupleEvent;


	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(m5,1);
	//ap->AddModule(jf,1);
	ap->AddModule(dump,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	ap->Run(runEvts);
	
//	gROOT->ProcessLine(".q");
}
