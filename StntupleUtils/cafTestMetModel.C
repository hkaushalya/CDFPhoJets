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
#include "samantha/samantha/Pho/EventFilter.hh"


void cafTestMetModel(int ind, int dataset, int subset, const int MinNjet15=1, 
						const int iPhoType=0, const int SampleId = 0, 
						const float metSigCut=0, const int debug = 0) 
// 10=pho data,11=pho mc, 
// 20,21,22=zee,zmm,ztt mc
// 30,31,32=wen,wmn,wtn mc
// debug = use to run on local dataset when testing
{
	std::cout << "*************************************************************" << std::endl;
	std::cout << "******** Input parameters for " << __FUNCTION__ << "*********" << std::endl;
	std::cout << "******** ind, dataset, subset = " << ind << " / " << dataset << " / " << subset << std::endl;
	std::cout << "******** MinNjet15 = " << MinNjet15 << std::endl;
	std::cout << "******** iPhoType  = " << iPhoType << std::endl;
	std::cout << "******** SampleId  = " << SampleId << std::endl;
	std::cout << "******** MetSigCut = " << metSigCut << std::endl;
	std::cout << "******** debug     = " << debug << std::endl;
	std::cout << "*************************************************************" << std::endl;

	
	if (! debug)
	{
		assert ((dataset == 10 || dataset ==  11 || dataset == 20 || dataset == 21 || dataset == 22 ||
		 		dataset ==  30 || dataset == 31 || dataset == 32) && "Invalid dataset requested!");
	}
	assert (MinNjet15 >= 0 && "MinNjet Cannot be less than 0");
	
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	int TOTALFILES = 0;
	int Ndatasets = 0;
	std::string suffix;
	TStnDataset *dsp[Ndatasets];
	TStnAna* ap = NULL;
	
if (debug)
{
	TChain* chain = new TChain("STNTUPLE");
	if (dataset == 10)
	{
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/ch031a01.0247ph1a");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/ch031a53.020fph1a");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/ch031b6a.00bbph1a");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/ch031bf4.005eph1a");
	} else if (dataset == 11)
	{/*
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd028ec7.0006exo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd02ac2e.0004exo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd02b43e.0008exo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd02b473.0006exo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd02b4da.0001exo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd02b553.0001exo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd02b5dc.0001exo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd02ba61.001cexo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd02bae9.0004exo8");
		*/
		std::cout << "Using the debug file" << std::endl;
		//chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd025748.0006exo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd028c17.0001exo8");
		
	} else if (dataset == 12)
	{
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/MetTestEvt/cdfpqcd.jqcdfh.154126.4053");
	} else if (dataset == 13) // debug JetFilterV2 vector out of range problem
	{
		std::cout << "USING MET DEBUG FILES" << std::endl;
//		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/MetDebug/pd02588b.0001exo8");
//		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/MetDebug/pd02ba32.0007exo8");
//		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/MetDebug/pd02bb63.000aexo8");
		//chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/MetDebug/pd028a4f.0009exo8");
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/MetDebug/pd02baa7.000dexo8");
		
	} else {

		std::cout << "NO VALID DATASET SPECIFIED for DEBUG mode to process. returning!" << std::endl;
		return;
	}


	
	ap = new TStnAna(chain);
}
else
{

	switch (dataset) {
		case 10: //photon data
		{
			if (subset == 1) { TOTALFILES = 700; Ndatasets = 1; }
			else if (subset == 2) { TOTALFILES = 1545; Ndatasets = 1; }
			else if (subset == 3) { TOTALFILES = 1229; Ndatasets = 1; }
			else if (subset == 4) { TOTALFILES = 888; Ndatasets = 1; }
			else {
				TOTALFILES = 4362;
				Ndatasets = 4;
			}
			suffix = "phodata";
			break;
		}
	
		case 11:  //photon mc
		{
			TOTALFILES = 484;
			Ndatasets = 1;
			suffix = "phomc";
			break;
		}
		case 20:			//zee 
		{
			if (subset == 1) { TOTALFILES = 150; Ndatasets = 1; }
			else if (subset ==2)	{ TOTALFILES = 447; Ndatasets = 1; } 
			else if (subset ==3)	{ TOTALFILES = 137; Ndatasets = 1; } 
			else if (subset ==4)	{ TOTALFILES = 132; Ndatasets = 1; } 
			else if (subset ==5)	{ TOTALFILES = 228; Ndatasets = 1; } 
			else if (subset ==6)	{ TOTALFILES = 214; Ndatasets = 1; } 
			else if (subset ==7)	{ TOTALFILES = 130; Ndatasets = 1; } 
			else if (subset ==8)	{ TOTALFILES = 244; Ndatasets = 1; } 
			else { TOTALFILES = 1682; Ndatasets = 8; }
			suffix = "zeemc";
			break;
		}
		case 21: 		//zmm
		{
			if (subset == 1) { TOTALFILES = 142; Ndatasets = 1; }
			else if (subset ==2)	{ TOTALFILES = 440; Ndatasets = 1; } 
			else if (subset ==3)	{ TOTALFILES = 122; Ndatasets = 1; }
			else if (subset ==4)	{ TOTALFILES = 118; Ndatasets = 1; }
			else if (subset ==5)	{ TOTALFILES = 210; Ndatasets = 1; }
			else if (subset ==6)	{ TOTALFILES = 206; Ndatasets = 1; }
			else if (subset ==7)	{ TOTALFILES = 125; Ndatasets = 1; }
			else if (subset ==8)	{ TOTALFILES = 226; Ndatasets = 1; }
			else  { TOTALFILES = 1589; Ndatasets = 8; }
			suffix = "zmmmc";
			break;
			
		}
		case 22: 	//ztt
		{
			if (subset == 1) { TOTALFILES = 431; Ndatasets = 1; }
			else if (subset ==2)	{ TOTALFILES = 429; Ndatasets = 1; } 
			else  { TOTALFILES = 860; Ndatasets = 2; }
			suffix = "zttmc";
			break;
		}
		case 30: 	//wen
		{
			if (subset == 1) { 
				std::cout << " Dataset/subset " << dataset << "/" << subset << " is no longer exists. Returning!" << std::endl; 
				TOTALFILES = 0; Ndatasets = 1; return; }
			else if (subset ==2)	{ TOTALFILES = 646; Ndatasets = 1; } 
			else if (subset ==3)	{ TOTALFILES = 494; Ndatasets = 1; }
			else if (subset ==4)	{ TOTALFILES = 398; Ndatasets = 1; }
			else if (subset ==5)	{ TOTALFILES = 124; Ndatasets = 1; }
			else if (subset ==6)	{ TOTALFILES = 232; Ndatasets = 1; }
			else  { TOTALFILES = 1894; Ndatasets = 5; }
			suffix = "wenmc";
			break;
		}
		case 31: 	//wmn
		{
			if (subset == 1) { TOTALFILES = 127; Ndatasets = 1; }
			else if (subset ==2)	{ TOTALFILES = 378; Ndatasets = 1; } 
			else if (subset ==3)	{ TOTALFILES = 249; Ndatasets = 1; }
			else if (subset ==4)	{ TOTALFILES = 389; Ndatasets = 1; }
			else if (subset ==5)	{ TOTALFILES = 121; Ndatasets = 1; }
			else if (subset ==6)	{ TOTALFILES = 227; Ndatasets = 1; }
			else  { TOTALFILES = 1491; Ndatasets = 6; }

			suffix = "wmnmc";
			break;
		}
		case 32: 	//wtn
		{
			TOTALFILES = 420;
			Ndatasets = 1;
			suffix = "wtnmc";
			break;
		}
		default :
		{
			std::cout << "No matching dataset found. exiting!" << std::endl;
			return;
		}

		
	}


	std::cout << "Total Files "<< TOTALFILES <<" from " << suffix << " "<< Ndatasets << " dataset/s to be attached." << std::endl;


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
			} else if (subset == 4) {
				dsp[0] = new TStnDataset(); //p14-17
				c->InitDataset(dsp[0],"cdfpstn","cph1ak","","",252836,261005);
			} else {
				std::cout << "Attaching all datasets.." << std::endl;
				dsp[0] = new TStnDataset(); //p1-4
				dsp[1] = new TStnDataset("cdfpstn","cph1ai"); //p5-10
				dsp[2] = new TStnDataset(); //p11-13
				dsp[3] = new TStnDataset(); //p14-17
				c->InitDataset(dsp[0],"cdfpstn","cph1ah","","",190851,203779);	//dropped first 400pb-1
				c->InitDataset(dsp[1]);
				c->InitDataset(dsp[2],"cdfpstn","cph1aj","","",233133,246231);
				c->InitDataset(dsp[3],"cdfpstn","cph1ak","","",252836,261005);
			}

		} else if (dataset == 11) { 
			std::cout<<"DATASET photon MC pypho22_dfc"<<std::endl;
			dsp[0] = new TStnDataset();
			c->InitDataset(dsp[0],"stntuple/dev_242","pypho22_dfc","","");

		} else if (dataset == 20) { 
			std::cout<<"zee"<<std::endl;

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","ze1s6d");
				c->InitDataset(dsp[0]);
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","ze1sad");
				c->InitDataset(dsp[0]);
			} else if (subset == 3) {
				dsp[0] = new TStnDataset("cdfpstn","ze0scd");
				c->InitDataset(dsp[0]);
			} else if (subset == 4) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sdd");
				c->InitDataset(dsp[0]);
			} else if (subset == 5) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sed");
				c->InitDataset(dsp[0]);
			} else if (subset == 6) {
				dsp[0] = new TStnDataset("cdfpstn","ze0see");
				c->InitDataset(dsp[0]);
			} else if (subset == 7) {
				dsp[0] = new TStnDataset("cdfpstn","ze0seh");
				c->InitDataset(dsp[0]);
			} else if (subset ==8) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sej");
				c->InitDataset(dsp[0]);
			} else {
			
				dsp[0] = new TStnDataset("cdfpstn","ze1s6d");
				dsp[1] = new TStnDataset("cdfpstn","ze1sad");
				dsp[2] = new TStnDataset("cdfpstn","ze0scd");
				dsp[3] = new TStnDataset("cdfpstn","ze0sdd");
				dsp[4] = new TStnDataset("cdfpstn","ze0sed");
				dsp[5] = new TStnDataset("cdfpstn","ze0see");
				dsp[6] = new TStnDataset("cdfpstn","ze0seh");
				dsp[7] = new TStnDataset("cdfpstn","ze0sej");
				for (int i=0; i< Ndatasets; i++) 
						c->InitDataset(dsp[i]);
			}

		} else if (dataset == 21) { 
			std::cout<<"zmm"<<std::endl;

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","ze1s6m");
				c->InitDataset(dsp[0]);
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","ze1s9m");
				c->InitDataset(dsp[0]);
			} else if (subset == 3) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sbm");
				c->InitDataset(dsp[0]);
			} else if (subset == 4) {
				dsp[0] = new TStnDataset("cdfpstn","ze0scm");
				c->InitDataset(dsp[0]);
			} else if (subset == 5) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sdm");
				c->InitDataset(dsp[0]);
			} else if (subset == 6) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sem");
				c->InitDataset(dsp[0]);
			} else if (subset == 7) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sfm");
				c->InitDataset(dsp[0]);
			} else if (subset == 8) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sgm");
				c->InitDataset(dsp[0]);
			} else {
			
				dsp[0] = new TStnDataset("cdfpstn","ze1s6m");
				dsp[1] = new TStnDataset("cdfpstn","ze1s9m");
				dsp[2] = new TStnDataset("cdfpstn","ze0sbm");
				dsp[3] = new TStnDataset("cdfpstn","ze0scm");
				dsp[4] = new TStnDataset("cdfpstn","ze0sdm");
				dsp[5] = new TStnDataset("cdfpstn","ze0sem");
				dsp[6] = new TStnDataset("cdfpstn","ze0sfm");
				dsp[7] = new TStnDataset("cdfpstn","ze0sgm");
				for (int i=0; i< Ndatasets; i++) 
						c->InitDataset(dsp[i]);
			}

		} else if (dataset == 22) { 
			std::cout<<"ztt"<<std::endl;

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","ze0s8t");
				c->InitDataset(dsp[0]);
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sat");
				c->InitDataset(dsp[0]);
			} else {
				dsp[0] = new TStnDataset("cdfpstn","ze0s8t");
				dsp[1] = new TStnDataset("cdfpstn","ze0sat");
				for (int i=0; i< Ndatasets; i++) 
						c->InitDataset(dsp[i]);
			}

		} else if (dataset == 30) { 
			std::cout<<"wen"<<std::endl;

			if (subset == 1) {

			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","we0sge");
				c->InitDataset(dsp[0]);

			} else if (subset ==3) {
			dsp[0] = new TStnDataset("cdfpstn","we0she");
				c->InitDataset(dsp[0]);
			} else if (subset ==4) {
			dsp[0] = new TStnDataset("cdfpstn","we0sie");
				c->InitDataset(dsp[0]);
			} else if (subset ==5) {
			dsp[0] = new TStnDataset("cdfpstn","we0seh");
				c->InitDataset(dsp[0]);
			} else if (subset ==6) {
			dsp[0] = new TStnDataset("cdfpstn","we0sej");
				c->InitDataset(dsp[0]);
			} else {
			
				//dsp[0] = new TStnDataset("cdfpstn","we0sfe");
				dsp[0] = new TStnDataset("cdfpstn","we0sge");
				dsp[1] = new TStnDataset("cdfpstn","we0she");
				dsp[2] = new TStnDataset("cdfpstn","we0sie");
				dsp[3] = new TStnDataset("cdfpstn","we0seh");
				dsp[4] = new TStnDataset("cdfpstn","we0sej");
				for (int i=0; i< Ndatasets; i++) 
						c->InitDataset(dsp[i]);
			}
					
		} else if (dataset == 31) { 
			std::cout<<"wmn"<<std::endl;

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","we0s7m");
				c->InitDataset(dsp[0]);
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","we0s8m");
				c->InitDataset(dsp[0]);
			} else if (subset == 3) {
				dsp[0] = new TStnDataset("cdfpstn","we0s9m");
				c->InitDataset(dsp[0]);
			} else if (subset == 4) {
				dsp[0] = new TStnDataset("cdfpstn","we0sam");
				c->InitDataset(dsp[0]);
			} else if (subset == 5) {
				dsp[0] = new TStnDataset("cdfpstn","we0sbm");
				c->InitDataset(dsp[0]);
			} else if (subset == 6) {
				dsp[0] = new TStnDataset("cdfpstn","we0sgm");
				c->InitDataset(dsp[0]);
			} else {

				dsp[0] = new TStnDataset("cdfpstn","we0s7m");
				dsp[1] = new TStnDataset("cdfpstn","we0s8m");
				dsp[2] = new TStnDataset("cdfpstn","we0s9m");
				dsp[3] = new TStnDataset("cdfpstn","we0sam");
				dsp[4] = new TStnDataset("cdfpstn","we0sbm");
				dsp[5] = new TStnDataset("cdfpstn","we0sgm");
				for (int i=0; i< Ndatasets; i++) 
						c->InitDataset(dsp[i]);
			}
					
		} else if (dataset == 32) { 
			std::cout<<"wtn"<<std::endl;
			//dsp[0] = new TStnDataset("cdfpstn","we0s9t");
			dsp[0] = new TStnDataset("cdfpstn","we0sat");
			for (int i=0; i< Ndatasets; i++) 
					c->InitDataset(dsp[i]);
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


	if (!debug)
	{
		int split = 0;
		if (dataset == 10) split = 300;
		else split = 50;						// all MC datasets have more than 100 files.
		std::cout << "Splitting job to " << split << std::endl;

		if (TOTALFILES < split ) 
		{
			std::cout << __FILE__ << "::" << __LINE__ << "::ERROR::Splitting the job to " << split << " but datasets have only " << TOTALFILES << " files! exiting." << std::cout;
		}
		ap->SetSplit(ind,split);
	}

}
	
	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

	//******************************************//
	//for stripping evts
	//----  need this part for event preselection only
//   TStnOutputModule* out= new TStnOutputModule("MetJet_Dump.root");
//	out->SetMaxFileSize(300);
//	ap->SetOutputModule(out);
	
	//******************************************//		 

	

  /******************************************************/
  // define all module here
  /******************************************************/

	bool bSummaryStat = 0; 	// 1 = no summary
  
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
		trigMod->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
	
	InitSuperPhotons* myinit = new InitSuperPhotons;
		myinit->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	TagTightPhotons* tagTight = new TagTightPhotons;
		tagTight->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
		tagTight->SetDebug(0);

	TagLoosePhotons* tagLoose = new TagLoosePhotons;
		tagLoose->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
	
	TagElectrons* tagTightEle = new TagElectrons;
		tagTightEle->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
	TagLooseElectrons* tagLooseEle = new TagLooseElectrons;
		tagLooseEle->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	PhoEleFilter *peFilter = new PhoEleFilter;
		peFilter->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary


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

	int NpointsToGenerate = 1;		//need at least 1 for  pseudo experitments
	//int NpointsToGenerate = 10;//for calibration
	//int MinNjet15 = 1;
	int SelectSigMetEvent = 0; // =0 do nothing; =1 select events with MetSig>metsigCut
	//int SelectSigMetEvent = 0; // for calibration
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

	//jer settings
	bool oldGmean=1, oldGsigma=1, oldLmpv=1, oldLsigma=1, oldGnorm=1;
	
	
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
		jf->SetUseHadJets(0);		//0=use det jets 1=use had jets
		jf->SetDoubleSmearJets(0);		//0=smear once jets 1=smear twice
		//jf->SetMaxJetEta(3.6);   //Sasha wanted this for Fix3_2
		jf->SetZeroJERoffset(0);		//1=no offset ,0=use offset
		jf->SetMetCut(0);
		jf->SetPrintLevel(10);		//higher number => more prints
		jf->SetMEtAddMethod(2);		// how we add MEt to a jet 0=not added, 1=most mismeasure jet
											// 2= closeset jet to MEt

		jf->SetRemoveTightPhotons(1);
		//jf->SetRemoveLoosePhotons(0);
		//jf->SetRemoveSidebandPhotons(1); this does not work like this.
		// I need to move this selection to the point when the objects are
		//actually removed. not when they are stored in JetFilterV2.

		jf->UseOldJERForGaussMean(oldGmean);
		jf->UseOldJERForGaussSigma(oldGsigma);
		jf->UseOldJERForLandauMPV(oldLmpv);
		jf->UseOldJERForLandauSigma(oldLsigma);
		jf->UseOldJERForGaussNorm(oldGnorm);


	JetFilterModuleV2 *jf2 = new JetFilterModuleV2("JetFilterV2_2");  //---------- My Vertex Filter Initialization
   jf2->SetDebug(0);
		jf2->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf2->SetMaxNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf2->SetJTC_imode(JTC_imode); 			// 0==MC; 1==data
		jf2->SetUnclParamSwitch(UnclParamSwitch); 	// 0==photon parametrization; 1==Zee parametrization 
//		jf2->SetMaxDeltaPhiJMet(MaxDeltaPhiJMet); 	// controls dPhi in MyMetCleanUpCut
		jf2->SetNpointsToGenerate(NpointsToGenerate); 	// number of pseudo-experiments per each event
		jf2->SetSelectMetEvent(SelectMetEvent); 		// 0==do nothing; 1==select event with large MET
		jf2->SetRemoveDuplicate(RemoveDuplicate); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		jf2->SetRemoveBadMet(RemoveBadMet); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		///jf2->SetUseMetPDFscenario(1);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		jf2->SetDumpEvent(DumpEvent); 			// 0==do nothing; 1==dump event
		jf2->SetDumpEventFileName("V2DATA_dumpEvent.dat");
		jf2->SetLargeMetEventFileName("V2DATA_largeMetEvent.dat");
		jf2->SetDatFileName("V2DATA_jet.dat");
		
		//new stuff for v2
		jf2->SetSampleID(SampleId); // 0==signal
		//jf2->SetMinDeltaPhiJMet(-10.0); // every event should pass this cut
		//-----------------------------
		jf2->SetSelectSigMetEvent(SelectSigMetEvent); // =0 do nothing; =1 select events with MetSig>metsigCut
		jf2->SetMetSigCut(metSigCut+2);
		jf2->SetAnalysisMode(AnalysisMode);      // =0 doesn't fill AnalysisHisto, but fills MetStudyHisto
												// =1 fills AnalysisHisto, but doesn't fill MetStudyHisto
		//-----------------------------
		jf2->SetSelectExoticEvent(SelectExoticEvent);		// i have commented out these parts 
		jf2->SetDoVxSwap(DoVxSwap); // 0=do nothing; 1=swap vertices to minimize MET
		jf2->SetUseHadJets(0);		//0=use det jets 1=use had jets
		jf2->SetDoubleSmearJets(0);		//0=smear once jets 1=smear twice
		//jf2->SetMaxJetEta(3.6);   //Sasha wanted this for Fix3_2
		jf2->SetZeroJERoffset(0);		//1=no offset ,0=use offset
		jf2->SetMetCut(0);
		jf2->SetPrintLevel(10);		//higher number => more prints
		jf2->SetMEtAddMethod(2);		// how we add MEt to a jet 0=not added, 1=most mismeasure jet
											// 2= closeset jet to MEt

		jf2->SetRemoveTightPhotons(1);

		jf2->UseOldJERForGaussMean(oldGmean);
		jf2->UseOldJERForGaussSigma(oldGsigma);
		jf2->UseOldJERForLandauMPV(oldLmpv);
		jf2->UseOldJERForLandauSigma(oldLsigma);
		jf2->UseOldJERForGaussNorm(oldGnorm);




	TagPMTSpikes* tagPMT = new TagPMTSpikes;
		tagPMT->SetMode(1);	//1= pass only if no spikes found
		tagPMT->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
	
	TagBeamHalo* tagBH = new TagBeamHalo();
		tagBH->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
	
	TEmTimingModule* EmT = new TEmTimingModule();
	
	SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
		setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt");
		setEmTime->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	TagPhoenixPhotons* tagPhoenix = new TagPhoenixPhotons;
		tagPhoenix->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
	
	TagConvElectrons* tagConvEle = new TagConvElectrons;
		tagConvEle->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	EventProperties * evtProp = new EventProperties;
		evtProp->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
		
	Pho2JetsTemp *phoSel = new Pho2JetsTemp();
		phoSel->SetPhoMinDetEta(0);
		phoSel->SetPhoMaxDetEta(0.95);
		phoSel->SetPhoMinEt(30.0);
		phoSel->SetHaloType(5);
		phoSel->SetPhotonType(iPhoType);		//0=signal, 1=sideband
		phoSel->SetSummaryStat(bSummaryStat);		//1=disable 0=enable endjob summary


	//EleJets *eleSel = new EleJets();
	
	EventFilter *filter = new EventFilter;
		//filter->SetSummaryStat(1);
		filter->SetDebug(0);
		//filter->AddEventRange(153739,8876,10000);
		filter->AddEvent(153739,8876);
		//filter->AddRun(153739);


	//ap->AddModule(passFailFile,1);
	if (dataset == 10) {
		ap->AddModule(ptChecker,1);
	} else {
		std::cout << __FILE__ << "::" << __LINE__ << ":: Not adding BadPointerChecker for MC data OR debug mode." << std::endl;
	}
	//ap->AddModule(filter,1);
	//ap->AddModule(f10,1);
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
	//ap->AddModule(evtProp,1);
	ap->AddModule(phoSel,1);
	
	//ap->AddModule(eleSel,1);
	ap->AddModule(jf,1);
	ap->AddModule(jf2,1);

	if (debug) ap->Run(debug);
	else ap->Run();
	
	//TO DEBUG THE INFINITE LOOPING OF THE JER INTEGRAL
	//ap->ProcessEvent(141619, 3124);

//event with the wrong MEt prob calc for the closest jet
//	ap->ProcessEvent(153416,34644);
//	ap->ProcessEvent(167023,17793);
	
	std::stringstream filename;
	filename << "MetTest_PhotonMC.root_" << ind; 
	ap->SaveHist(filename.str().c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
