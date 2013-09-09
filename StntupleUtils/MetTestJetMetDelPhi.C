#include "TBenchmark.h"
#include "TSystem.h"
#include <iostream>
#include <sstream>
#include <string>
#include "TObjArray.h"
#include "TIterator.h"
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
#include "samantha/samantha/Pho/TagPMTSpikes.hh"
#include "samantha/samantha/Pho/SetEmTimes.hh"
#include "samantha/samantha/utils/PassFailTree.hh"
#include "samantha/samantha/utils/BadPointerChecker.hh"
#include "samantha/samantha/Pho/Pho2JetsTemp.hh"
#include "samantha/samantha/Pho/EventProperties.hh"
#include "samantha/samantha/Pho/TagCosmicPhotons.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/FlatStupleMaker.hh"
#include "samantha/samantha/Pho/TagConvElectrons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"
#include "samantha/samantha/Pho/Filter10.hh"
#include "samantha/samantha/Pho/EleJets.hh"
#include "samantha/samantha/Pho/JetFilterModuleV2.hh"


void MetTestJetMetDelPhi(int ind, int dataset, int subset, const int MinNjet15=1, const int iPhoType=0, const int SampleId = 10, const float metSigCut=0, const int debug = 0) 
{
	assert (MinNjet15 >= 0 && "MinNjet Cannot be less than 0");
	std::cout << "I got ind, dataset, subset =" << ind <<", " << dataset << ", " << subset << std::endl;
	
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnAna* ap = NULL;
	
	TChain* chain = new TChain("STNTUPLE");
	//chain->Add("ExtremeMetJetDelPhi.root");

	for (int i=1; i <=100; ++i)
	{
		std::stringstream file;
		file << "~/RESULTS/MetTest/NewStuff/02062009_NoMetSigCutAndSTNskimOfMetJetDElPhi/STNS/MetJet_DelPhiExtremes.root_" << i;
		std::cout << file.str() << std::endl;
		chain->Add(file.str().c_str());
	}

	//chain->Merge("ExtremeMetJetDelPhi.root");

		ap = new TStnAna(chain);
		ap->GetInputModule()->SetPrintLevel(1);   // print file name as they are opened
					 
	
	TH1::AddDirectory(kFALSE); //do not add these histos to memory 


  /******************************************************/
  // define all module here
  /******************************************************/
  
  //PassFailTree *passFailFile = new PassFailTree;
  BadPointerChecker *ptChecker = new BadPointerChecker;
  Filter10 *f10 = new Filter10;
  
	TriggerModule* trigMod;
	trigMod = new TriggerModule;
		trigMod->SetGoodRunListFile("goodrun_v23_pho_00.txt");  // to p13
		trigMod->SetGoodRunBit(1);		//1= use good run, 0=do not use it
		trigMod->SetTriggerBit(1);		// 1= require tigger, 0=NO (for MC)
		trigMod->SetMinVtx(1);
		trigMod->SetUseVtxZcut(1);		//1= require z<60 0=not cut on z
		trigMod->SetMaxVtx(1000);
	
	InitSuperPhotons* myinit = new InitSuperPhotons;

	TagTightPhotons* tagTight = new TagTightPhotons;
	TagLoosePhotons* tagLoose = new TagLoosePhotons;
	
	TagElectrons* tagTightEle = new TagElectrons;
	TagLooseElectrons* tagLooseEle = new TagLooseElectrons;

	PhoEleFilter *peFilter = new PhoEleFilter;


	//these are common parameters for central/JES UP/JES DOWN jet collections
	//int JTC_imode = 0;		// 0==MC, 1==data
	//if (debug) std::cout  << " WARNING! RUNNING ON LOCAL DATASET. JTC_imode is default value (0 = MC)!" << std::endl;
	//if (dataset==10)  JTC_imode = 1;		// 0==MC, 1==data
	int UnclParamSwitch = 0; 	// 0==photon parametrization; 1==Zee parametrization 
	//int SampleId = 3;		// >=0 , <=4, 3- use Pythia pho+jet parameterization
	//int SampleId = 10;		// calibration mode
	//int SampleId = 3;		// 3- use Pythia for MC OR data signal pho+jet parameterization
								// 4 = DATA sideband pho+jet parameterization

	std::cout << __FILE__ << "::" << __LINE__ << "::SAMPLE ID = " << SampleId  << " for Dataset " << dataset << std::endl;

	int NpointsToGenerate = 1;		//need at least 1 for METsig calculation
	//int NpointsToGenerate = 10;//for calibration
	//int MinNjet15 = 1;
	int SelectSigMetEvent = 0; // =0 do nothing; =1 select events with MetSig>metsigCut
	//int SelectSigMetEvent = 0; // for calibration
	//float MetSigCut = 0.0;   //for calibration
	//float MetSigCut = metSigCut;
	int AnalysisMode=1;		//0==off 1==on
	//int AnalysisMode=0;		// for calibration
	int SelectMetEvent = 0;
	int SelectExoticEvent = 0;
	int RemoveDuplicate = 1;
	int RemoveBadMet = 1;// -1=select events with bad met; 1=remove events with bad met
	int DoVxSwap = 0; // 0=do nothing; 1=swap vertices to minimize MET
	//float MaxDeltaPhiJMet = 0.3;
	int DumpEvent = 0;

	
	JetFilterModuleV2 *jf = new JetFilterModuleV2();  //---------- My Vertex Filter Initialization
   jf->SetDebug(0);
		jf->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf->SetMaxNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf->SetJTC_imode(JTC_imode); 			// 0==MC; 1==data
		jf->SetUnclParamSwitch(UnclParamSwitch); 	// 0==photon parametrization; 1==Zee parametrization 
//		jf->SetMaxDeltaPhiJMet(MaxDeltaPhiJMet); 	// controls dPhi in MyMetCleanUpCut
		jf->SetNpointsToGenerate(NpointsToGenerate); 	// number of pseudo-experiments per each event
		jf->SetSelectMetEvent(SelectMetEvent); 		// 0==do nothing; 1==select event with large MET
		jf->SetRemoveDuplicate(RemoveDuplicate); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		jf->SetRemoveBadMet(RemoveBadMet); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		///jf->SetUseMetPDFscenario(1);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		jf->SetDumpEvent(DumpEvent); 			// 0==do nothing; 1==dump event
		jf->SetDumpEventFileName("V2DATA_dumpEvent.dat");
		jf->SetLargeMetEventFileName("V2DATA_largeMetEvent.dat");
		jf->SetDatFileName("V2DATA_jet.dat");
		
		//new stuff for v2
		jf->SetSampleID(SampleId); // 0==signal
		//jf->SetMinDeltaPhiJMet(-10.0); // every event should pass this cut
		//-----------------------------
		jf->SetSelectSigMetEvent(SelectSigMetEvent); // =0 do nothing; =1 select events with MetSig>metsigCut
		jf->SetMetSigCut(metSigCut);
		jf->SetAnalysisMode(AnalysisMode);      // =0 doesn't fill AnalysisHisto, but fills MetStudyHisto
												// =1 fills AnalysisHisto, but doesn't fill MetStudyHisto
		//-----------------------------
		jf->SetSelectExoticEvent(SelectExoticEvent);		// i have commented out these parts 
		jf->SetDoVxSwap(DoVxSwap); // 0=do nothing; 1=swap vertices to minimize MET

		//jf->SetMaxJetEta(3.6);   //Sasha wanted this for Fix3_2



	TagPMTSpikes* tagPMT = new TagPMTSpikes;
		tagPMT->SetMode(1);	//1= pass only if no spikes found
	
	TagBeamHalo* tagBH = new TagBeamHalo();
	
	TEmTimingModule* EmT = new TEmTimingModule();
	
	SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
		setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt");

	TagPhoenixPhotons* tagPhoenix = new TagPhoenixPhotons;
	
	TagConvElectrons* tagConvEle = new TagConvElectrons;

		EventProperties * evtProp = new EventProperties;
		
	Pho2JetsTemp *phoSel = new Pho2JetsTemp();
		phoSel->SetPhoMinDetEta(0);
		phoSel->SetPhoMaxDetEta(0.95);
		phoSel->SetPhoMinEt(30.0);
		phoSel->SetHaloType(5);
		phoSel->SetPhotonType(iPhoType);		//0=signal, 1=sideband

	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(tagPhoenix,1);
	ap->AddModule(peFilter,1);
	ap->AddModule(tagPMT,1);
	ap->AddModule(tagBH,1);
	ap->AddModule(EmT,1);
	ap->AddModule(setEmTime,1);
	ap->AddModule(tagConvEle,1);
	ap->AddModule(evtProp,1);
	ap->AddModule(phoSel,1);
	ap->AddModule(jf,1);

	if (debug) ap->Run(debug);
	else ap->Run();
  	
	
	std::string filename = "MetTest_PhotonMC.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
