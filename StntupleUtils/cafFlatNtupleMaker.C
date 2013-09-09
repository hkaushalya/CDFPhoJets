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
#include "samantha/samantha/Pho/HaloByCutsTemp_Pho.hh"
#include "samantha/samantha/utils/PassFailTree.hh"
#include "samantha/samantha/utils/BadPointerChecker.hh"
#include "samantha/samantha/Pho/HaloByCutsTemp_Ele.hh"
#include "samantha/samantha/Pho/Pho2JetsTemp.hh"
#include "samantha/samantha/Pho/EventProperties.hh"
#include "samantha/samantha/Pho/TagCosmicPhotons.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/FlatStupleMaker.hh"
#include "samantha/samantha/Pho/TagConvElectrons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"


void cafFlatNtupleMaker(int ind, int dataset, int subset, int cafsections=100, int debug=0) // 10=pho data,11=pho mc, 20,21,22=zee,zmm,ztt mc
															 // 30,31,32=wen,wmn,wtn mc
{
	std::cout << "I got ind, dataset, subset =" << ind <<", " << dataset << ", " << subset << std::endl;
	
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

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
				dsp[0] = new TStnDataset(); //p18-28
				c->InitDataset(dsp[0],"cdfpstn","cph1am","","",261119,289197);   //3942 files up to p26
				ap->AddDataset(dsp[0]);
				suffix += "_18P28";
			} else if (subset == 6) {
				dsp[0] = new TStnDataset(); //p29-30
				c->InitDataset(dsp[0],"cdfpstn","cph1ap","","",289273,293800);   //
				ap->AddDataset(dsp[0]);
				suffix += "_29P30";

			} else {
				std::cout << "Attaching all datasets.." << std::endl;
				const int Ndsp = 6; 
				dsp[0] = new TStnDataset(); //p1-4
				dsp[1] = new TStnDataset("cdfpstn","cph1ai"); //p5-10
				dsp[2] = new TStnDataset(); //p11-13
				dsp[3] = new TStnDataset(); //p14-17
				dsp[4] = new TStnDataset(); //p18-28
				dsp[5] = new TStnDataset(); //p29-30
				c->InitDataset(dsp[0],"cdfpstn","cph1ah","","",190851,203779);	//dropped first 400pb-1
				c->InitDataset(dsp[1]);
				c->InitDataset(dsp[2],"cdfpstn","cph1aj","","",233133,246231);
				c->InitDataset(dsp[3],"cdfpstn","cph1ak","","",252836,261005);
				c->InitDataset(dsp[4],"cdfpstn","cph1am","","",261119,289197);
				c->InitDataset(dsp[5],"cdfpstn","cph1ap","","",289273,293800);
				ap->AddDataset(dsp[0]);
				ap->AddDataset(dsp[1]);
				ap->AddDataset(dsp[2]);
				ap->AddDataset(dsp[3]);
				ap->AddDataset(dsp[4]);
				ap->AddDataset(dsp[5]);
			}

		} else if (dataset == 11) { 

			//if (subset == 1)
			//{
			//	dsp[0] = new TStnDataset();
			//	c->InitDataset(dsp[0],"stntuple/dev_242","pypho22_dfc","",""); //my first dataset used for all APS/ICHEP 2008 and upto March 2010
				//this is 0MB 
				//dsp[0] = new TStnDataset("cdfpstn","pexo8d");
			//	c->InitDataset(dsp[0]); 
			//	suffix = "Pythia_PhoJetNoMinBiasGen5_jqcdfh";
			//} else if (subset == 2)
			//{
				dsp[0] = new TStnDataset("cdfpstn","gq0sqd");  //PYTHIA PHO+JETS+MINBIAS
				c->InitDataset(dsp[0]);   //this had minbias
				suffix = "Pythia_PhoJetMinBiasGen6_gq0sqd";
			//} else if (subset == 3)
			//{
			//	dsp[0] = new TStnDataset("cdfpstn","hq0sqd");
			//	c->InitDataset(dsp[0]);   //HERWIG PHO+JET+MINBIAS (with JIMMY???)
			//	suffix = "Herwig_PhoJetMinBiasGen6_hq0sqd";
			//}
			
			std::cout << "Running over " << suffix << std::endl;
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
				dsp[0] = new TStnDataset("cdfpstn","we0sfe"); // 237 files   6135951 events
				c->InitDataset(dsp[0]);
				ap->AddDataset(dsp[0]);
				suffix = "we0sfe";
			} else if (subset == 2) {
				dsp[0] = new TStnDataset("cdfpstn","we0sge");  // 646 files  11664464 events
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

		} else if (dataset == 50) { 
			std::cout<<"WW diboson mc"<<std::endl; //WW->dilepton (e,u, or t) NLO S=1.27pb
			suffix = "MC_WW";
			const int Ndsp = 6;
			dsp[0] = new TStnDataset("cdfpstn","we0s5d");
			dsp[1] = new TStnDataset("cdfpstn","we0sbd");
			dsp[2] = new TStnDataset("cdfpstn","we0sgd");
			dsp[3] = new TStnDataset("cdfpstn","we0skd");
			dsp[4] = new TStnDataset("cdfpstn","we0snd");
			dsp[5] = new TStnDataset("cdfpstn","we0saf");
			for (int i=0; i< Ndsp; i++)
			{
				c->InitDataset(dsp[i]);
				ap->AddDataset(dsp[i]);
			}
		} else if (dataset == 51) { 
			std::cout<<"WZ diboson mc: WZ->dilepton (Z->dilepton, W->generic, filter on two e,u,or t with pT>1GeV. Pythia, NLO CS=0.365pb, Filter eff=0.76"<<std::endl;
			suffix = "MC_WZ";
			const int Ndsp = 6;
			dsp[0] = new TStnDataset("cdfpstn","we0s6d");
			dsp[1] = new TStnDataset("cdfpstn","we0scd");
			dsp[2] = new TStnDataset("cdfpstn","we0shd");
			dsp[3] = new TStnDataset("cdfpstn","we0sld");
			dsp[4] = new TStnDataset("cdfpstn","we0sod");
			dsp[5] = new TStnDataset("cdfpstn","we0sbf");
			for (int i=0; i< Ndsp; i++)
			{
				c->InitDataset(dsp[i]);
				ap->AddDataset(dsp[i]);
			}

		} else if (dataset == 52) { 
			std::cout<<"ZZ diboson mc: ZZ -> dilepton (Includes γ* component, both Z->generic, filter on two e or μ with pT > 1 GeV, Pythia NLO CS= 1.511 pb, Filter Efficiency = 0.23, M_{ll} > 15 GeV)"<< std::endl;
			suffix = "MC_ZZ";
			const int Ndsp = 6;
			dsp[0] = new TStnDataset("cdfpstn","we0s7d");
			dsp[1] = new TStnDataset("cdfpstn","we0sdd");
			dsp[2] = new TStnDataset("cdfpstn","we0sid");
			dsp[3] = new TStnDataset("cdfpstn","we0smd");
			dsp[4] = new TStnDataset("cdfpstn","we0spd");
			dsp[5] = new TStnDataset("cdfpstn","we0scf");
			for (int i=0; i< Ndsp; i++)
			{
				c->InitDataset(dsp[i]);
				ap->AddDataset(dsp[i]);
			}
		} else if (dataset == 60) { 
			std::cout<<"ttbar mc: Pythia NLO CS= 0.68778 pb: period 1-11"<<std::endl;
			suffix = "MC_TTbar";
			const int Ndsp = 1;
			dsp[0] = new TStnDataset("cdfpstn","te0s2z"); 
			for (int i=0; i< Ndsp; i++)
			{
				c->InitDataset(dsp[i]);
				ap->AddDataset(dsp[i]);
			}

		} else if (dataset == 70) { 
			std::cout<<"Costas W/g mc: Pythia+MAD event"<<std::endl;
			suffix = "MC_Wgamma_PytMad";
			const int Ndsp = 1;
			dsp[0] = new TStnDataset("cdfpstn","ws06jj"); 
			for (int i=0; i< Ndsp; i++)
			{
				c->InitDataset(dsp[i]);
				ap->AddDataset(dsp[i]);
			}
		} else if (dataset == 80) { 
			std::cout<<"Costas Z/g mc: Pythia+MAD event"<<std::endl;
			suffix = "MC_Zgamma_PytMad";
			const int Ndsp = 1;
			dsp[0] = new TStnDataset("cdfpstn","zs06jj"); 
			for (int i=0; i< Ndsp; i++)
			{
				c->InitDataset(dsp[i]);
				ap->AddDataset(dsp[i]);
			}
		}

		std::cout << "suffix = " << suffix << std::endl;

/*		int split = 0;
		if (dataset == 10 && subset == 0) split = 2000;
		else if (dataset == 10 && subset == 1) split = 200;
		else if (dataset == 10 && subset == 5) split = 500;
		else split = 100;						// all MC datasets have more than 100 files.
		std::cout << "******** Splitting job to " << split << std::endl;
		ap->SetSplit(ind,split);
*/
		std::cout << "******** SPLITTING JOB TO " << cafsections << " ******** "<< std::endl;
		ap->SetSplit(ind,cafsections);


		TH1::AddDirectory(kFALSE); //do not add these histos to memory 

		//******************************************//
		//for stripping evts
		//----  need this part for event preselection only
		//TStnOutputModule* out= new TStnOutputModule("CosmicsAndJetEmMismatch.root");
		//out->SetMaxFileSize(300);
		//ap->SetOutputModule(out);

		//******************************************//		 



		/******************************************************/
		// define all module here
		/******************************************************/

		//these are common parameters for central/JES UP/JES DOWN jet collections
		bool bSummaryStat = 0; 	// 1 = no summary
		const float fMaxEmEta=1.1;
		const int MinNjet15 = 0;
		const float fMaxJetEta=3.0;
		const int AnalysisMode=2;		//0==off 1==on, >1 = generate corrected jets and MET/SUMET
		//ok, this MinErThr in JetMod3 is used only to calculate Ht from Jets Et>fMinEtThr, which is set to 15GeV by default.
		//It got nothing to do with the EM objetcs I passed in. It uses all the EM objects passed in to calucalte Ht.
		//MET scenario does not matter in this stripped version. As long as I specify the scenario 3 (Et>15Gev jets)
		//when retreiving MET I am fine. This MET will include corrections to all the EM objects I passed, how low their Et is.
		//08-26-2010
		
		const float fMinEmEtThr = 25.0; //em objects with Et>, lowest should 25 if trigger is required as iso25 will make sure there is a object with 25GeV
		//do not lower this or pass EM objects wit hEt below this to be removed from jet list.
		//this will cause the Ht/met calculation to change unpredictably!
		const bool bRemAllEm = true;
		const int iMetScenario =3; // 0--raw Met; corrected for jets with: 1-- et>5, 2-- et>10, 3-- et>15, 4-- et>20
		//this does not matter in ntuple making. but use '3' which is Et>15GeV objects for Ht and MET
		//when using ht/met in final plots



		//PassFailTree *passFailFile = new PassFailTree;
		BadPointerChecker *ptChecker = new BadPointerChecker;

		TriggerModule* trigMod;
		trigMod = new TriggerModule;
		//trigMod->SetGoodRunListFile("goodrun_v19_pho_00.txt");  // to 
		//trigMod->SetGoodRunListFile("goodrun_v23_pho_00.txt");  // to 
		//trigMod->SetGoodRunListFile("goodrun_v29_pho_00.txt");  // to p25
	//	trigMod->SetGoodRunListFile("goodrun_v25_pho_00.txt");  // to p25
		trigMod->SetGoodRunListFile("goodrun_v36_pho_00.txt");  // to and including p30
		trigMod->SetNeedGoodRun(1);		//1= use good run, 0=do not use it
		trigMod->SetNeedTriggers(1);		// 1= require tigger, 0=NO (for MC) //DO NOT TURN THIS OFF!
		trigMod->SetNeedMinVtx(0);
		trigMod->SetNeedVtxZcut(0);				//1= require z<60 0=not cut on z
		trigMod->SetNeedMaxVtx(1000);
		trigMod->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

		InitSuperPhotons* myinit = new InitSuperPhotons;
		myinit->SetPrintLevel(0);
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
		peFilter->SetMinEmEt(fMinEmEtThr);
		peFilter->SetMaxEmDetEta(fMaxEmEta);

	TagBeamHalo* tagBH = new TagBeamHalo();
		tagBH->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	TagPMTSpikes* tagPMT = new TagPMTSpikes;
		tagPMT->SetMode(1);	//1= pass only if no spikes found
		tagPMT->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
		
	TagPhoenixPhotons* tagPhoenix = new TagPhoenixPhotons;
		tagPhoenix->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	TagConvElectrons* tagConvEle = new TagConvElectrons;
		tagConvEle->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
	

	TEmTimingModule* EmT = new TEmTimingModule();
	SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
		//setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt"); // use no correction. does not matter to me
		setEmTime->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary



	JetFilterModuleV3 *jf3 = new JetFilterModuleV3();  //---------- My Vertex Filter Initialization
   	jf3->SetDebug(0);
		jf3->SetJTC_systcode(0);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf3->SetMaxJetEta(fMaxJetEta);
		jf3->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf3->SetMinEtThr(fMinJetEtThr);		//Et>15.0 is the default for min Et for jets used to calculate Ht
		jf3->SetRemoveTightPhotons(bRemAllEm);
		jf3->SetRemoveLoosePhotons(bRemAllEm);
		jf3->SetRemoveTightPhoLikeElectrons(bRemAllEm);
		jf3->SetRemoveLoosePhoLikeElectrons(bRemAllEm);
		jf3->SetRemoveStdLooseElectrons(bRemAllEm);
		jf3->SetPrintLevel(0);		//higher number => more prints
		jf3->SetSummaryStat(bSummaryStat);
		jf3->SetMyMetScenario(iMetScenario);


	JetFilterModuleV3 *jf3Up = new JetFilterModuleV3("JetFilterV3JESup");  //---------- My Vertex Filter Initialization
   	jf3Up->SetDebug(0);
		jf3Up->SetJTC_systcode(1);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf3Up->SetMaxJetEta(fMaxJetEta);
		jf3Up->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf3Up->SetMinEtThr(fMinJetEtThr);		//Et>15.0 is the default for min Et for jets used to calculate Ht
		jf3Up->SetRemoveTightPhotons(bRemAllEm);
		jf3Up->SetRemoveLoosePhotons(bRemAllEm);
		jf3Up->SetRemoveTightPhoLikeElectrons(bRemAllEm);
		jf3Up->SetRemoveLoosePhoLikeElectrons(bRemAllEm);
		jf3Up->SetRemoveStdLooseElectrons(bRemAllEm);
		jf3Up->SetPrintLevel(0);		//higher number => more prints
		jf3Up->SetSummaryStat(bSummaryStat);
		jf3Up->SetMyMetScenario(iMetScenario);

	JetFilterModuleV3 *jf3Down = new JetFilterModuleV3("JetFilterV3JESdown");  //---------- My Vertex Filter Initialization
   	jf3Down->SetDebug(0);
		jf3Down->SetJTC_systcode(-1);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf3Down->SetMaxJetEta(fMaxJetEta);
		jf3Down->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf3Down->SetMinEtThr(fMinJetEtThr);		//Et>15.0 is the default for min Et for jets used to calculate Ht
		jf3Down->SetRemoveTightPhotons(bRemAllEm);
		jf3Down->SetRemoveLoosePhotons(bRemAllEm);
		jf3Down->SetRemoveTightPhoLikeElectrons(bRemAllEm);
		jf3Down->SetRemoveLoosePhoLikeElectrons(bRemAllEm);
		jf3Down->SetRemoveStdLooseElectrons(bRemAllEm);
		jf3Down->SetPrintLevel(0);		//higher number => more prints
		jf3Down->SetSummaryStat(bSummaryStat);
		jf3Down->SetMyMetScenario(iMetScenario);
		
	FlatStupleMaker *flatS = new FlatStupleMaker;
		const float fStupleVer = 8;
		std::stringstream stupleName;
		stupleName << "StupleV"<< fStupleVer << "_" << suffix << "_" << ind << ".root";
		flatS->SetStupleName(stupleName.str());
		flatS->SetSummaryStat(bSummaryStat);
		flatS->SetStupleVersion(fStupleVer);


	//ap->AddModule(passFailFile,1);
	ap->AddModule(ptChecker,1);
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(tagPhoenix,1);
	ap->AddModule(tagConvEle,1);
	ap->AddModule(peFilter,1);
	ap->AddModule(tagPMT,1);
	ap->AddModule(tagBH,1);
	ap->AddModule(EmT,1);
	ap->AddModule(setEmTime,1);
	ap->AddModule(jf3,1);
	ap->AddModule(jf3Up,1);
	ap->AddModule(jf3Down,1);
	ap->AddModule(flatS,1);

	ap->GetInputModule()->SetPrintLevel(100);	// print file name as they are opened

	if (debug>0) ap->Run(debug);
	else ap->Run();
  	
	std::stringstream filename;
	filename  << "StupleMaking" << suffix << "_" << subset << "_" << ind << ".root";
	ap->SaveHist(filename.str().c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename.str() << std::endl;
	//std::cout << "::::::::: Created Pass/Fail file => "<< passFailFile->GetFileName() << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
