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
#include "samantha/samantha/Pho/EleJets.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"
#include "samantha/samantha/Pho/JetFilterModuleV3.hh"
#include "Stntuple/Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/Stntuple/loop/TStnOutputModule.hh"
#include "Stntuple/base/base/TStnDataset.hh"
#include "TSystem.h"

/////////////////////////////////////////////////////////////

void runEleJets_caf(const int nEvts = 100, const int dataset = 10,
							const int subset = 1,const int cafind=0, 
							const int NcafSections=20) 
{
	std::cout << "Nevts to run = " << nEvts << std::endl;
	std::cout << "cafind/total sections  = " << cafind << " / " << NcafSections << std::endl;

	TStnAna*     ap = new TStnAna();
	TStnCatalog* c = new TStnCatalog();
	int Ndatasets = 10;	//just some max number. Ndsp will controll the end point to read
	std::string suffix("");
	TStnDataset *dsp[Ndatasets];

		if (dataset == 10) {
			std::cout<<"photon data"<<std::endl;
			suffix = "DATA_Pho";
			if (subset == 1) {
				dsp[0] = new TStnDataset(); //p1-4
				dsp[0]->SetPrintLevel(100);
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

			if (subset == 1)
			{
				dsp[0] = new TStnDataset();
				c->InitDataset(dsp[0],"stntuple/dev_242","pypho22_dfc","",""); //my first dataset used for all APS/ICHEP 2008 and upto March 2010
				//this is 0MB 
				//dsp[0] = new TStnDataset("cdfpstn","pexo8d");
				c->InitDataset(dsp[0]); 
				suffix = "Pythia_PhoJetNoMinBiasGen5_jqcdfh";
			} else if (subset == 2)
			{
				dsp[0] = new TStnDataset("cdfpstn","gq0sqd");  //PYTHIA PHO+JETS+MINBIAS
				c->InitDataset(dsp[0]);   //this had minbias
				suffix = "Pythia_PhoJetMinBiasGen6_gq0sqd";
			} else if (subset == 3)
			{
				dsp[0] = new TStnDataset("cdfpstn","hq0sqd");
				c->InitDataset(dsp[0]);   //HERWIG PHO+JET+MINBIAS (with JIMMY???)
				suffix = "Herwig_PhoJetMinBiasGen6_hq0sqd";
			}
			
			ap->AddDataset(dsp[0]);

		} else if (dataset == 20) { 
			std::cout<<"zee"<<std::endl;
			suffix = "MC_Zee";

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","ze1s6d");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze1s6d";
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","ze1sad");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze1sad";
			} else if (subset == 3) {
				dsp[0] = new TStnDataset("cdfpstn","ze0scd");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0scd";
			} else if (subset == 4) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sdd");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0sdd";
			} else if (subset == 5) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sed");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0sed";
			} else if (subset == 6) {
				dsp[0] = new TStnDataset("cdfpstn","ze0see");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0see";
			} else if (subset == 7) {
				dsp[0] = new TStnDataset("cdfpstn","ze0seh");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0seh";
			} else if (subset ==8) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sej");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0sej";
			} else {

				const int Ndsp = 8; 
				dsp[0] = new TStnDataset("cdfpstn","ze1s6d");
				dsp[1] = new TStnDataset("cdfpstn","ze1sad");
				dsp[2] = new TStnDataset("cdfpstn","ze0scd");
				dsp[3] = new TStnDataset("cdfpstn","ze0sdd");
				dsp[4] = new TStnDataset("cdfpstn","ze0sed");
				dsp[5] = new TStnDataset("cdfpstn","ze0see");
				dsp[6] = new TStnDataset("cdfpstn","ze0seh");
				dsp[7] = new TStnDataset("cdfpstn","ze0sej");
				for (int i=0; i< Ndsp; i++) 
				{
					c->InitDataset(dsp[i]);
					ap->AddDataset(dsp[i]);
				}
			}

		} else if (dataset == 21) { 
			std::cout<<"zmm"<<std::endl;
			suffix = "MC_Zmm";

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","ze1s6m");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze1s6m";
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","ze1s9m");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze1s9m";
			} else if (subset == 3) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sbm");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0sbm";
			} else if (subset == 4) {
				dsp[0] = new TStnDataset("cdfpstn","ze0scm");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0scm";
			} else if (subset == 5) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sdm");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0sdm";
			} else if (subset == 6) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sem");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0sem";
			} else if (subset == 7) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sfm");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0sfm";
			} else if (subset == 8) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sgm");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix += "_ze0sgm";
			} else {

				const int Ndsp = 8;
				dsp[0] = new TStnDataset("cdfpstn","ze1s6m");
				dsp[1] = new TStnDataset("cdfpstn","ze1s9m");
				dsp[2] = new TStnDataset("cdfpstn","ze0sbm");
				dsp[3] = new TStnDataset("cdfpstn","ze0scm");
				dsp[4] = new TStnDataset("cdfpstn","ze0sdm");
				dsp[5] = new TStnDataset("cdfpstn","ze0sem");
				dsp[6] = new TStnDataset("cdfpstn","ze0sfm");
				dsp[7] = new TStnDataset("cdfpstn","ze0sgm");
				for (int i=0; i< Ndsp; i++) 
				{
					c->InitDataset(dsp[i]);
					ap->AddDataset(dsp[i]);
				}
			}

		} else if (dataset == 22) { 
			std::cout<<"ztt"<<std::endl;
			suffix = "MC_Ztt";

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","ze0s8t");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "ze0s8t";
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sat");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "ze0sat";
			} else if (subset == 3) {
				dsp[0] = new TStnDataset("cdfpstn","ze0sbt");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "ze0sbt";

			} else {
				const int Ndsp = 3;
				dsp[0] = new TStnDataset("cdfpstn","ze0s8t");
				dsp[1] = new TStnDataset("cdfpstn","ze0sat");
				dsp[2] = new TStnDataset("cdfpstn","ze0sbt");
				for (int i=0; i< Ndsp; i++) 
				{
					c->InitDataset(dsp[i]);
					ap->AddDataset(dsp[i]);
				}
			}

		} else if (dataset == 30) { 
			std::cout<<"wen"<<std::endl;
			suffix = "MC_Wen";

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","we0sfe");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sfe";
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","we0sge");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sge";

			} else if (subset ==3) {
				dsp[0] = new TStnDataset("cdfpstn","we0she");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0she";
			} else if (subset ==4) {
				dsp[0] = new TStnDataset("cdfpstn","we0sie");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sie";
			} else if (subset ==5) {
				dsp[0] = new TStnDataset("cdfpstn","we0seh");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0seh";
			} else if (subset ==6) {
				dsp[0] = new TStnDataset("cdfpstn","we0sej");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sej";
			} else {

				const int Ndsp = 6;
				dsp[0] = new TStnDataset("cdfpstn","we0sfe");
				dsp[1] = new TStnDataset("cdfpstn","we0sge");
				dsp[2] = new TStnDataset("cdfpstn","we0she");
				dsp[3] = new TStnDataset("cdfpstn","we0sie");
				dsp[4] = new TStnDataset("cdfpstn","we0seh");
				dsp[5] = new TStnDataset("cdfpstn","we0sej");
				for (int i=0; i< Ndsp; i++) 
				{
					c->InitDataset(dsp[i]);
					ap->AddDataset(dsp[i]);
				}
			}

		} else if (dataset == 31) { 
			std::cout<<"wmn"<<std::endl;
			suffix = "MC_Wmn";

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","we0s7m");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0s7m";
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","we0s8m");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0s8m";
			} else if (subset == 3) {
				dsp[0] = new TStnDataset("cdfpstn","we0s9m");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0s9m";
			} else if (subset == 4) {
				dsp[0] = new TStnDataset("cdfpstn","we0sam");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sam";
			} else if (subset == 5) {
				dsp[0] = new TStnDataset("cdfpstn","we0sbm");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sbm";
			} else if (subset == 6) {
				dsp[0] = new TStnDataset("cdfpstn","we0sgm");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sgm";
			} else {

				const int Ndsp = 6;
				dsp[0] = new TStnDataset("cdfpstn","we0s7m");
				dsp[1] = new TStnDataset("cdfpstn","we0s8m");
				dsp[2] = new TStnDataset("cdfpstn","we0s9m");
				dsp[3] = new TStnDataset("cdfpstn","we0sam");
				dsp[4] = new TStnDataset("cdfpstn","we0sbm");
				dsp[5] = new TStnDataset("cdfpstn","we0sgm");
				for (int i=0; i< Ndsp; i++) 
				{
					c->InitDataset(dsp[i]);
					ap->AddDataset(dsp[i]);
				}
			}

		} else if (dataset == 32) { 
			std::cout<<"wtn"<<std::endl;
			suffix = "MC_Wtn";

			if (subset == 1) {
				dsp[0] = new TStnDataset("cdfpstn","we0s9t");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0s9t";
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","we0sat");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sat";
			} else if (subset == 3) {
				dsp[0] = new TStnDataset("cdfpstn","we0sbt");
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sbt";
			} else {
				const int Ndsp = 3;
				dsp[0] = new TStnDataset("cdfpstn","we0s9t");
				dsp[1] = new TStnDataset("cdfpstn","we0sat");
				dsp[2] = new TStnDataset("cdfpstn","we0sbt");
				for (int i=0; i< Ndsp; i++)
				{
					c->InitDataset(dsp[i]);
					ap->AddDataset(dsp[i]);
				}
			}
			
		} else if (dataset == 40) { 
			std::cout<<"dipho mc"<<std::endl;
			suffix = "MC_DiPho";
			const int Ndsp = 1;
			dsp[0] = new TStnDataset("cdfpstn","gx0s1g");
			for (int i=0; i< Ndsp; i++)
			{
				c->InitDataset(dsp[i]);
				ap->AddDataset(dsp[i]);
			}

		}

		std::cout << " >>>>>>>> suffix = " << suffix << std::endl;

	
	ap->SetSplit(cafind, NcafSections);

	bool bSummaryStat = 0; 	// 1 = no summary

	TriggerModule* trigMod;
	trigMod = new TriggerModule;
	trigMod->SetGoodRunListFile("goodrun_v31_pho_00.txt");  // to p31
	trigMod->SetGoodRunBit(1);				//1= use good run, 0=do not use it
	trigMod->SetTriggerBit(1);				// 1= require tigger, 0=NO (for MC)
	trigMod->SetMinVtx(1);
	trigMod->SetUseVtxZcut(1);				//1= require z<60 0=not cut on z
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
	peFilter->SetMinEmEt(0.);
	peFilter->SetMaxEmDetEta(4.0);

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

	int MinNjet15 = 1;
	JetFilterModuleV3 *jf = new JetFilterModuleV3();  //---------- My Vertex Filter Initialization
   jf->SetDebug(0);
		jf->SetJTC_systcode(0);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		jf->SetMaxJetEta(3.2);
		//jf->SetRemoveTightPhotons(1);
		jf->SetRemoveTightPhoLikeElectrons(1);
		jf->SetPrintLevel(0);
		jf->SetMyMetScenario(3);
		jf->SetSummaryStat(bSummaryStat);

	
	EleJets *eleSel = new EleJets();
	eleSel->SetEleMinDetEta(0);
	eleSel->SetEleMaxDetEta(1.1);
	eleSel->SetEleMinEt(30.0);
	eleSel->SetRemoveConvEles(1);
	eleSel->SetHaloType(5);
	eleSel->SetDataSample(dataset);
	eleSel->UseNvtxWeigts(1);
	eleSel->SetSummaryStat(bSummaryStat);		//1=disable 0=enable endjob summary

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
	ap->AddModule(jf,1);
	ap->AddModule(eleSel,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened

	if (nEvts>0) ap->Run(nEvts); 
	else ap->Run(); 

	std::stringstream file;
	file << "EleJets_" << suffix <<"_" << cafind << "_" << NcafSections << ".root";
	ap->SaveHist(file.str().c_str(),2);
	std::cout << "Histograms save to file " << file.str() << std::endl;

}

