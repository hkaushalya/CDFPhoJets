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
#include "samantha/samantha/Pho/TMyJetFilterModule.hh"
#include "samantha/samantha/Pho/TagPMTSpikes.hh"
#include "samantha/samantha/Pho/SetEmTimes.hh"
#include "samantha/samantha/utils/PassFailTree.hh"
#include "samantha/samantha/utils/BadPointerChecker.hh"
#include "samantha/samantha/Pho/EventProperties.hh"
#include "samantha/samantha/Pho/TagCosmicPhotons.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/FlatStupleMaker.hh"
#include "samantha/samantha/Pho/TagConvElectrons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"
#include "samantha/samantha/Pho/SkimCosmics.hh"


void cafCosmicSkim(int ind, int dataset, int subset) // 10=pho data,11=pho mc, 20,21,22=zee,zmm,ztt mc
															 // 30,31,32=wen,wmn,wtn mc
{
	std::cout << "I got ind, dataset, subset =" << ind <<", " << dataset << ", " << subset << std::endl;
	
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	int TOTALFILES = 0;
	int Ndatasets = 0;
	std::string suffix;

	switch (dataset) {
		case 10: //photon data
		{
			if (subset == 1) { TOTALFILES = 700; Ndatasets = 1; }
			else if (subset == 2) { TOTALFILES = 1545; Ndatasets = 1; }
			else if (subset == 3) { TOTALFILES = 1229; Ndatasets = 1; }
			else {
				TOTALFILES = 3474;
				Ndatasets = 3;
			}
			suffix = "phodata";
			break;
		}
	
		default :
		{
			std::cout << "No matching dataset found. exiting!" << std::endl;
			return;
		}
	}


	std::cout << "Total Files "<< TOTALFILES <<" from " << Ndatasets << " dataset/s to be attached." << std::endl;

	TStnDataset *dsp[Ndatasets];


	TStnAna* ap = NULL;
	
	Int_t totalfiles = 0;

	int loop = 0;
	do {
		loop++;
		if (loop > 50) {
			std::cout << __FILE__ << "::" << __LINE__ << "Failed get the correct files after 20 attempts. Exiting.." << std::endl;
			return;
		}
		if (loop >1) {
			std::cout <<"\t loop=" << loop << ":: sleeping for 10 seconds" << std::endl;
			sleep(10); 
		}
	
		ap = new TStnAna();
		ap->GetInputModule()->SetPrintLevel(1);   // print file name as they are opened
					 
		TStnCatalog* c = new TStnCatalog();
		
		if (dataset == 10) {
			std::cout<<"photon data"<<std::endl;
			if (subset == 1) {
				dsp[0] = new TStnDataset(); //p1-4
				c->InitDataset(dsp[0],"cdfpstn","cph1ah","","",190851,203779);	//dropped first 400pb-1  // 700 files
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","cph1ai"); //p5-10				//1545 files
				c->InitDataset(dsp[0]);
			} else if (subset == 3) {
				dsp[0] = new TStnDataset(); //p11-13
				c->InitDataset(dsp[0],"cdfpstn","cph1aj","","",233133,246231);  // 1229 files
			} else {
				dsp[0] = new TStnDataset(); //p1-4
				dsp[1] = new TStnDataset("cdfpstn","cph1ai"); //p5-10
				dsp[2] = new TStnDataset(); //p11-13
				c->InitDataset(dsp[0],"cdfpstn","cph1ah","","",190851,203779);	//dropped first 400pb-1
				c->InitDataset(dsp[1]);
				c->InitDataset(dsp[2],"cdfpstn","cph1aj","","",233133,246231);
			}
		}
		
		totalfiles = 0;			 
		std::cout << "total files needed and datasets="<<TOTALFILES <<","<<Ndatasets<<std::endl;
	  	for (int i=0; i < Ndatasets; i++)  {
			totalfiles += dsp[i]->GetNFiles();
			std::cout << "\t"<<i << "\t" <<dsp[i]->GetNFiles()<<std::endl;
			ap->AddDataset(dsp[i]);
		}

		std::cout << "Total Files = " << totalfiles << std::endl;
	} while (totalfiles != TOTALFILES);


	ap->SetSplit(ind,36);
	

	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

	//******************************************//
	//for stripping evts
	//----  need this part for event preselection only
	std::ostringstream sfile;
	sfile << "Cosmics_" << suffix <<"_"<< subset << "_" << ind << ".root";
	
   TStnOutputModule* out= new TStnOutputModule(sfile.str().c_str());
	out->SetMaxFileSize(300);
	ap->SetOutputModule(out);
	
	//******************************************//		 

	

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

	PhoEleFilter *peFilter = new PhoEleFilter;


	//these are common parameters for central/JES UP/JES DOWN jet collections
	int MinNjet15 = 0;
	int JTC_imode = 1;		// 0==MC, 1==data
	int UnclParamSwitch = 0;
	float MaxDeltaPhiJMet = 0.3;
	int NpointsToGenerate = 0;
	int SelectMetEvent = 0;
	int RemoveDuplicate = 0;
	int RemoveBadMet = 0;
	int DumpEvent = 0;

	if (dataset==10)  JTC_imode = 1;		// 0==MC, 1==data
	else  JTC_imode = 0;		// 0==MC, 1==data
	
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
		

	TagPMTSpikes* tagPMT = new TagPMTSpikes;
		tagPMT->SetMode(1);	//1= pass only if no spikes found
	
	TEmTimingModule* EmT = new TEmTimingModule();
	
	SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
		setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt");

	TagPhoenixPhotons* tagPhoenix = new TagPhoenixPhotons;
	
	SkimCosmics* cosmicSkim = new SkimCosmics;

	//ap->AddModule(passFailFile,1);
	ap->AddModule(ptChecker,1);
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(peFilter,1);
	ap->AddModule(m5,1);
	ap->AddModule(tagPMT,1);
	ap->AddModule(EmT,1);
	ap->AddModule(setEmTime,1);
	ap->AddModule(tagPhoenix,1);
	ap->AddModule(cosmicSkim,1);

	ap->Run();
  	
	
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
