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

void TestLooseEleTag(unsigned runEvts=100000, std::string dataset="phodata")
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnAna* ap;
	TChain *chain = 0;

		/* TO TEST USE A LOCAL FILE */ 
		chain = new TChain("STNTUPLE");

				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.000eph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0020ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0037ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0058ph1a");
/*				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.00caph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.00ecph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.00fdph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.011fph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0163ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0187ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.01baph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.01dcph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0200ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0276ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.02b7ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.02c7ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0335ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.03dfph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0406ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0413ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.053bph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.05ccph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.069cph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.06aeph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.06edph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0852ph1a");
*/


		ap = new TStnAna(chain);


	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

  /******************************************************/
  // define all module here
  /******************************************************/
  
  //PassFailTree *passFailFile = new PassFailTree;
  BadPointerChecker *ptChecker = new BadPointerChecker;
  
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

	TagBeamHalo* tagBH = new TagBeamHalo();

	TagPMTSpikes* tagPMT = new TagPMTSpikes;
		tagPMT->SetMode(1);	//1= pass only if no spikes found
		
	TagPhoenixPhotons* phoenix = new TagPhoenixPhotons;

	TagConvElectrons* tagConvEle = new TagConvElectrons;
	
	TEmTimingModule* EmT = new TEmTimingModule();
	SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
		setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt");
		
	
	//ap->AddModule(passFailFile,1);
	//ap->AddModule(ptChecker,1);
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	//ap->AddModule(tagTight,1);
	//ap->AddModule(tagLoose,1);
	//ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	//ap->AddModule(tagPMT,1);
	//ap->AddModule(tagBH,1);
	//ap->AddModule(phoenix,1);
	//ap->AddModule(tagConvEle,1);
	//ap->AddModule(EmT,1);
	//ap->AddModule(setEmTime,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	ap->Run(runEvts);
	
	std::string filename = "StdLooseEleTest.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	//std::cout << "::::::::: Created Pass/Fail file => "<< passFailFile->GetFileName() << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
