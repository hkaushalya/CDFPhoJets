#include "TSystem.h"
#include "TChain.h"
#include "Stntuple/Stntuple/loop/TStnAna.hh"
#include <iostream>
#include <sstream>
#include <string>
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
#include "samantha/samantha/Pho/PhoJetsCount.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"
#include "samantha/samantha/Pho/JetFilterModuleV3.hh"
#include "Stntuple/Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/Stntuple/loop/TStnOutputModule.hh"
#include "Stntuple/base/base/TStnDataset.hh"
#include "TSystem.h"

/////////////////////////////////////////////////////////////

void runPhoJetCount(const int nEvts = 100, const int dataset=10, const int subset=1, const int cafind=0, const int NcafSections=20) 
{
	std::cout << "Nevts to run = " << nEvts << std::endl;
	std::cout << "cafind/total sections  = " << cafind << " / " << NcafSections << std::endl;

	int Ndatasets = 10;	//just some max number. Ndsp will controll the end point to read
	std::string suffix;

	TStnDataset *dsp[Ndatasets];

	TStnAna* ap = new TStnAna();
	
		TStnCatalog* c = new TStnCatalog();
		//c->SetPrintLevel(100);
		
		if (dataset == 10) {
			std::cout<<"photon data"<<std::endl;
			suffix = "DATA_Pho";
			if (subset == 1) {
				dsp[0] = new TStnDataset(); //p1-4
				//dsp[0]->SetPrintLevel(100);
				c->InitDataset(dsp[0],"cdfpstn","cph1ah","","",190851,203779);	//dropped first 400pb-1  // 700 files
				ap->AddDataset(dsp[0]);
				suffix += "_1P4";
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","cph1ai"); //p5-10				//1545 files
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_5P10";
			} else if (subset == 3) {
				dsp[0] = new TStnDataset(); //p11-13
				c->InitDataset(dsp[0],"cdfpstn","cph1aj","","",233133,246231);  // 1229 files
				ap->AddDataset(dsp[0]);
				suffix += "_11P13";
			} else if (subset == 4) {
				dsp[0] = new TStnDataset(); //p14-17
				c->InitDataset(dsp[0],"cdfpstn","cph1ak","","",252836,261005);		//888 files
				ap->AddDataset(dsp[0]);
				suffix += "_14P17";
			} else if (subset == 5) {
				dsp[0] = new TStnDataset(); //p18-26
				c->InitDataset(dsp[0],"cdfpstn","cph1am","","",261119,284843);   //3942 files
				ap->AddDataset(dsp[0]);
				suffix += "_18P26";
			} else {
				std::cout << "Attaching all datasets.." << std::endl;
				const int Ndsp = 5; 
				dsp[0] = new TStnDataset(); //p1-4
				dsp[1] = new TStnDataset("cdfpstn","cph1ai"); //p5-10
				dsp[2] = new TStnDataset(); //p11-13
				dsp[3] = new TStnDataset(); //p14-17
				dsp[4] = new TStnDataset(); //p18-26
				c->InitDataset(dsp[0],"cdfpstn","cph1ah","","",190851,203779);	//dropped first 400pb-1
				c->InitDataset(dsp[1]);
				c->InitDataset(dsp[2],"cdfpstn","cph1aj","","",233133,246231);
				c->InitDataset(dsp[3],"cdfpstn","cph1ak","","",252836,261005);
				c->InitDataset(dsp[4],"cdfpstn","cph1am","","",261119,274055);
				c->InitDataset(dsp[5],"cdfpstn","cph1am","","",261119,284843);
				ap->AddDataset(dsp[0]);
				ap->AddDataset(dsp[1]);
				ap->AddDataset(dsp[2]);
				ap->AddDataset(dsp[3]);
				ap->AddDataset(dsp[4]);
			}

		} else if (dataset == 11) { 

			dsp[0] = new TStnDataset("cdfpstn","gq0sqd");  //PYTHIA PHO+JETS+MINBIAS
			c->InitDataset(dsp[0]);   //this had minbias
			suffix = "Pythia_PhoJetMinBiasGen6_gq0sqd";
			std::cout << "Running over " << suffix << std::endl;
			ap->AddDataset(dsp[0]);

		} 
		std::cout << "suffix = " << suffix << std::endl;

		std::cout << "******** SPLITTING JOB TO " << NcafSections << " ******** "<< std::endl;
		ap->SetSplit(cafind, NcafSections);


		TH1::AddDirectory(kFALSE); //do not add these histos to memory 



	

	bool bSummaryStat = 0; 	// 1 = no summary



	TriggerModule* trigMod;
		trigMod = new TriggerModule;
		trigMod->SetGoodRunListFile("goodrun_v31_pho_00.txt");  // to p31
		trigMod->SetNeedGoodRun(0);				//1= use good run, 0=do not use it
		trigMod->SetNeedTriggers(0);				// 1= require tigger, 0=NO (for MC)
		trigMod->SetNeedMinVtx(0);
		trigMod->SetNeedMaxVtx(1000);
		trigMod->SetNeedVtxZcut(0);				//1= require z<60 0=not cut on z
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
	peFilter->SetMinEmEt(0.);
	peFilter->SetMaxEmDetEta(4.0);

	TagPMTSpikes* tagPMT = new TagPMTSpikes;
	tagPMT->SetMode(0);	//1= pass only if no spikes found
	tagPMT->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	TagBeamHalo* tagBH = new TagBeamHalo();
	tagBH->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	TEmTimingModule* EmT = new TEmTimingModule();

	SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
	setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt");
	setEmTime->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	TagPhoenixPhotons* tagPhoenix = new TagPhoenixPhotons;
	tagPhoenix->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	int MinNjet15 = 0;
	JetFilterModuleV3 *jf = new JetFilterModuleV3();  //---------- My Vertex Filter Initialization
   jf->SetDebug(0);
		jf->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf->SetMaxNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		jf->SetMaxJetEta(3.0);
		//remove all em objects for accurate met calculations. as I do when making flat ntuples
		jf->SetRemoveTightPhotons(1);
		jf->SetRemoveLoosePhotons(1);
		jf->SetRemoveTightPhoLikeElectrons(1);
		jf->SetRemoveLoosePhoLikeElectrons(1);
		jf->SetRemoveStdLooseElectrons(1);
		jf->SetPrintLevel(1);

	
	PhoJetsCount *phoSel = new PhoJetsCount();
		phoSel->SetPhoMinDetEta(0);
		//phoSel->SetPhoMaxDetEta(1.1);
		phoSel->SetPhoMaxDetEta(3.0);
		phoSel->SetPhoMinEt(30.0);
		phoSel->SetHaloType(5);
		phoSel->SetPhotonType(0);		//0=signal, 1=sideband
		phoSel->SetSummaryStat(bSummaryStat);		//1=disable 0=enable endjob summary

	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(tagPhoenix,1);
	//ap->AddModule(peFilter,1);
	ap->AddModule(tagPMT,1);
	ap->AddModule(tagBH,1);
	ap->AddModule(EmT,1);
	ap->AddModule(setEmTime,1);
	ap->AddModule(jf,1);
	ap->AddModule(phoSel,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened

	if (nEvts>0) ap->Run(nEvts); 
	else ap->Run(); 

	std::stringstream file;
	file << "PhoJetsCount_" << cafind << "Of" << NcafSections << ".root";
	ap->SaveHist(file.str().c_str(),2);
	std::cout << "Histograms save to file " << file.str() << std::endl;

}

