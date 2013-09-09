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

//Writing this module to quikly skim some events for cosmic/beam halo etc

void SkimSomeEvents(unsigned runEvts=100000, std::string dataset="phodata")
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnAna* ap;
	TChain *chain = 0;

		/* TO TEST USE A LOCAL FILE */ 
		chain = new TChain("STNTUPLE");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_0.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_1.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_10.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_11.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_12.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_13.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_14.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_15.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_16.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_17.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_18.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_19.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_2.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_20.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_21.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_22.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_23.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_24.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_25.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_26.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_27.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_28.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_29.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_3.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_30.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_31.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_32.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_33.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_34.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_35.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_4.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_5.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_6.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_7.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_8.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_1_9.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_2.root");
		chain->Add("/data/nbay02/a/samantha/RESULTS/Cosmics/Cosmics_phodata_3.root");


	ap = new TStnAna(chain);


	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

	//******************************************//
	//for stripping evts
	//----  need this part for event preselection only
   TStnOutputModule* out= new TStnOutputModule("Cosmic_Skim.root");
	out->SetMaxFileSize(300);
	ap->SetOutputModule(out);
	
	//******************************************//		 

  /******************************************************/
  // define all module here
  /******************************************************/
  
  //PassFailTree *passFailFile = new PassFailTree;
  //BadPointerChecker *ptChecker = new BadPointerChecker;
  
	TriggerModule* trigMod;
	trigMod = new TriggerModule;
		trigMod->SetGoodRunListFile("goodrun_v19_pho_00.txt");  // to p13
		trigMod->SetGoodRunBit(1);		//1= use good run, 0=do not use it
		trigMod->SetTriggerBit(1);		// 1= require tigger, 0=NO
		trigMod->SetMinVtx(0);
		trigMod->SetMaxVtx(1000);
		trigMod->SetUseVtxZcut(0);		//1= require z<60 0=not cut on z
	
	InitSuperPhotons* myinit = new InitSuperPhotons;

	TagTightPhotons* tagTight = new TagTightPhotons;
	TagLoosePhotons* tagLoose = new TagLoosePhotons;
	
	TagElectrons* tagTightEle = new TagElectrons;
	TagLooseElectrons* tagLooseEle = new TagLooseElectrons;
											//1,2 - pass only if a electron found and remove others from list

	PhoEleFilter *peFilter = new PhoEleFilter;

	TagBeamHalo* tagBH = new TagBeamHalo();

	TagPMTSpikes* tagPMT = new TagPMTSpikes;
		tagPMT->SetMode(1);	//1= pass only if no spikes found
		
//	TagPhoenixPhotons* phoenix = new TagPhoenixPhotons;

//	TagConvElectrons* tagConvEle = new TagConvElectrons;
	

	TEmTimingModule* EmT = new TEmTimingModule();
	SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
		setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt");
		
//	FlatStupleMaker *flatS = new FlatStupleMaker;	
//		flatS->SetStupleName("Stuple_Karen.root");
	
	//ap->AddModule(passFailFile,1);
//	ap->AddModule(ptChecker,1);
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
//	ap->AddModule(tagTight,1);
//	ap->AddModule(tagLoose,1);
//	ap->AddModule(tagTightEle,1);
//	ap->AddModule(tagLooseEle,1);
//	ap->AddModule(peFilter,1);
	ap->AddModule(tagPMT,1);
	ap->AddModule(tagBH,1);
//	ap->AddModule(phoenix,1);
//	ap->AddModule(tagConvEle,1);
	ap->AddModule(EmT,1);
	ap->AddModule(setEmTime,1);
//	ap->AddModule(flatS,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	ap->Run(runEvts);
  	
	
	std::string filename = "FlatStupleTest.root";
	//ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	//std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	//std::cout << "::::::::: Created Pass/Fail file => "<< passFailFile->GetFileName() << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
