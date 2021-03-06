#pragma warning(push)		//disable warnings for these headers only
	#pragma  warning(disable:4512)
	#pragma warning(disable:4180)
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
#pragma warning(pop) 	//restore orginal warning level		
	
#include "samantha/samantha/Pho/InitSuperPhotons.hh"
#include "samantha/samantha/Pho/TagTightPhotons.hh"
#include "samantha/samantha/Pho/TagLoosePhotons.hh"
#include "samantha/samantha/Pho/TagElectrons.hh"
#include "samantha/samantha/Pho/TagLooseElectrons.hh"
#include "samantha/samantha/Pho/TagBeamHalo.hh"
#include "samantha/samantha/Pho/TriggerModule.hh"
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
#include "samantha/samantha/Pho/FlatStupleMaker_MC_EWK.hh"
#include "samantha/samantha/Pho/TagConvElectrons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"

void cafZFlatNtupleMaker(int ind, int partype, int syscode)
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");


	std::ostringstream ostr;
	ostr << ind;
	std::string base_dir("./");
	std::string str_ind = ostr.str();


	
	int datasets = 0;

		if (partype == 0) { //zee
			datasets = 8;
		} else if (partype == 1) { //zmumu
			datasets = 8;
		} else if (partype == 2) { //ztautau
			datasets = 2;
		}

	TStnAna* ap = new TStnAna();
	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	
	TStnCatalog* c = new TStnCatalog();
		
	TStnDataset *dsp[datasets];
	
		if (partype == 0) {
			dsp[0] = new TStnDataset("cdfpstn","ze1s6d");
			dsp[1] = new TStnDataset("cdfpstn","ze1sad");
			dsp[2] = new TStnDataset("cdfpstn","ze0scd");
			dsp[3] = new TStnDataset("cdfpstn","ze0sdd");
			dsp[4] = new TStnDataset("cdfpstn","ze0sed");
			dsp[5] = new TStnDataset("cdfpstn","ze0see");
			dsp[6] = new TStnDataset("cdfpstn","ze0seh");
			dsp[7] = new TStnDataset("cdfpstn","ze0sej");

			c->InitDataset(dsp[0]);
			c->InitDataset(dsp[1]);
			c->InitDataset(dsp[2]);
			c->InitDataset(dsp[3]);
			c->InitDataset(dsp[4]);
			c->InitDataset(dsp[5]);
			c->InitDataset(dsp[6]);
			c->InitDataset(dsp[7]);
			ap->AddDataset(dsp[0]);
			ap->AddDataset(dsp[1]);
			ap->AddDataset(dsp[2]);
			ap->AddDataset(dsp[3]);
			ap->AddDataset(dsp[4]);
			ap->AddDataset(dsp[5]);
			ap->AddDataset(dsp[6]);
			ap->AddDataset(dsp[7]);
		}

		if (partype == 1) {
			dsp[0] = new TStnDataset("cdfpstn","ze1s6m");
			dsp[1] = new TStnDataset("cdfpstn","ze1s9m");
			dsp[2] = new TStnDataset("cdfpstn","ze0sbm");
			dsp[3] = new TStnDataset("cdfpstn","ze0scm");
			dsp[4] = new TStnDataset("cdfpstn","ze0sdm");
			dsp[5] = new TStnDataset("cdfpstn","ze0sem");
			dsp[6] = new TStnDataset("cdfpstn","ze0sfm");
			dsp[7] = new TStnDataset("cdfpstn","ze0sgm");
			c->InitDataset(dsp[0]);
			c->InitDataset(dsp[1]);
			c->InitDataset(dsp[2]);
			c->InitDataset(dsp[3]);
			c->InitDataset(dsp[4]);
			c->InitDataset(dsp[5]);
			c->InitDataset(dsp[6]);
			c->InitDataset(dsp[7]);
			ap->AddDataset(dsp[0]);
			ap->AddDataset(dsp[1]);
			ap->AddDataset(dsp[2]);
			ap->AddDataset(dsp[3]);
			ap->AddDataset(dsp[4]);
			ap->AddDataset(dsp[5]);
			ap->AddDataset(dsp[6]);
			ap->AddDataset(dsp[7]);
		}

		if (partype == 2) {
			dsp[0] = new TStnDataset("cdfpstn","ze0s8t");
			dsp[1] = new TStnDataset("cdfpstn","ze0sat");
			c->InitDataset(dsp[0]);
			c->InitDataset(dsp[1]);
			ap->AddDataset(dsp[0]);
			ap->AddDataset(dsp[1]);
		}


		double sum=0;
		for (int i=0; i < datasets; i++) {
			std::cout << "Nfile="<< dsp[i]->GetNFiles() << std::endl;
			std::cout << "NEvts="<< dsp[i]->GetNEvents() << std::endl;
			sum+= dsp[i]->GetNEvents();
			TObjArray* myarr = dsp[i]->GetListOfFiles();
			TIterator *it = myarr->MakeIterator();
			TObject *obj;
			while (obj = (TObject*) it->Next()) {
				obj->Print();
			}
		}
		std::cout << "Total Evts="<< sum << std::endl;




	int split = 0;

	if (partype == 0 || partype == 1) split = 200;
	if (partype == 2) split = 10;	
	
		std::cout << "split, IND=" << ind << ", " << split <<std::endl;
		ap->SetSplit(ind, split);
		std::cout << "syscode, MomPdg, parType=" << syscode << ", " << 23 << ", " << partype <<std::endl;



	TH1::AddDirectory(kFALSE); //do not add these histos to memory 
	

  /******************************************************/
  // define all module here
  /******************************************************/
  
  //PassFailTree *passFailFile = new PassFailTree;
  //BadPointerChecker *ptChecker = new BadPointerChecker;
  
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

	PhoEleFilter* pef = new PhoEleFilter;

	TMyJetFilterModule *m5 = new TMyJetFilterModule();  //---------- My Vertex Filter Initialization
		m5->SetMinNjet15(1);				//Jets required to pass (Et>15GeV)
		m5->SetJTC_imode(0); 			// 0==MC; 1==data
		m5->SetJTC_systcode(syscode);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		m5->SetUnclParamSwitch(0); 	// 0==photon parametrization; 1==Zee parametrization 
		m5->SetMaxDeltaPhiJMet(0.3); 	// controls dPhi in MyMetCleanUpCut
		m5->SetNpointsToGenerate(0); 	// number of pseudo-experiments per each event
		m5->SetSelectMetEvent(0); 		// 0==do nothing; 1==select event with large MET
		m5->SetRemoveDuplicate(0); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		m5->SetRemoveBadMet(0); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		///m5->SetUseMetPDFscenario(1);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		m5->SetDumpEvent(0); 			// 0==do nothing; 1==dump event
		m5->SetDumpEventFileName("DATA_dumpEvent.dat");
		m5->SetLargeMetEventFileName("DATA_largeMetEvent.dat");
		m5->SetDatFileName("DATA_jet.dat");
		

	//TagPMTSpikes* tagPMT = new TagPMTSpikes;
		//tagPMT->SetMode(1);	//1= pass only if no spikes found
	
	TagBeamHalo* tagBH = new TagBeamHalo();
	
	//TEmTimingModule* EmT = new TEmTimingModule();
	
	//SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
		//setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt");

	TagPhoenixPhotons* tagPhoenix = new TagPhoenixPhotons;
	
	TagConvElectrons* tagConvEle = new TagConvElectrons;

	FlatStupleMaker_MC_EWK *flatS = new FlatStupleMaker_MC_EWK;	
		flatS->SetMomPDG(23);
		flatS->SetDecayType(partype);			//0=z->ee,1=z->nunu,2=z->tautau
		std::string Stuple;
			if (partype == 0) Stuple="Stuple_Zee.root"; 
			if (partype == 1) Stuple="Stuple_Zmm.root"; 
			if (partype == 2) Stuple="Stuple_Ztt.root"; 

		std::string stuplefile = base_dir + Stuple + str_ind;
		flatS->SetStupleName(stuplefile.c_str()); 
	
	//ap->AddModule(passFailFile,1);
	//ap->AddModule(ptChecker,1);
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(pef,1);
	ap->AddModule(m5,1);
	//ap->AddModule(tagPMT,1);
	ap->AddModule(tagBH,1);
	//ap->AddModule(EmT,1);
	//ap->AddModule(setEmTime,1);
	ap->AddModule(tagPhoenix,1);
	ap->AddModule(tagConvEle,1);
	ap->AddModule(flatS,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	ap->Run();
  	
	
	std::string filename = "FlatStuple.root"+str_ind;
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	//std::cout << "::::::::: Created Pass/Fail file => "<< passFailFile->GetFileName() << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
	//gROOT->ProcessLine(".q");
}
