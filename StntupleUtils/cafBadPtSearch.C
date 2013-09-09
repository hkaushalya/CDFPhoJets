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
#include "samantha/samantha/utils/BadPointerChecker.hh"
/*********************************************************************
 * created this to find the run/sections with bad pointers
 * Running the module to skip thiese is very expensive and this
 * happend very rarely. So I though I could find these and exclude them
 * from processing.
 *********************************************************************/

void cafBadPtSearch(int ind, int dataset, int subset, const int debug = 0) 
// 10=pho data,11=pho mc, 
// 20,21,22=zee,zmm,ztt mc
// 30,31,32=wen,wmn,wtn mc
// debug = use to run on local dataset when testing
{
	if (! debug)
	{
		assert ((dataset == 10 || dataset ==  11 || dataset == 20 || dataset == 21 || dataset == 22 ||
		 		dataset ==  30 || dataset == 31 || dataset == 32) && "Invalid dataset requested!");
	}
	std::cout << "I got ind, dataset, subset =" << ind <<", " << dataset << ", " << subset << std::endl;
	
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
	chain->Add("/data/nbay02/a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/ch031a01.0247ph1a");
	chain->Add("/data/nbay02/a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/ch031a53.020fph1a");
	chain->Add("/data/nbay02/a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/ch031b6a.00bbph1a");
	chain->Add("/data/nbay02/a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/ch031bf4.005eph1a");
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


	std::cout << "Total Files "<< TOTALFILES <<" from " << Ndatasets << " dataset/s to be attached." << std::endl;


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


	/*int split = 0;
	if (dataset == 10) split = 300;
	else split = 200;						// all MC datasets have more than 100 files.
	std::cout << "Splitting job to " << split << std::endl;
	ap->SetSplit(ind,split);
	*/

}
	
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
  
  BadPointerChecker *ptChecker = new BadPointerChecker;
  		ap->AddModule(ptChecker,1);

	if (debug) ap->Run(5000);
	else ap->Run();
  	
	
	std::string filename = "BadPointerCheck.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
