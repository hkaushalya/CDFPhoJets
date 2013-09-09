#include <iostream>
#include <sstream>
#include <string>
#include "Stntuple/Stntuple/loop/TStnAna.hh"
#include "Stntuple/Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/Stntuple/loop/TStnOutputModule.hh"
#include "Stntuple/Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/base/base/TStnDataset.hh"
#include "samantha/samantha/Pho/TriggerModule.hh"


void TestMetModel(int nevets) 
{
	TChain* chain = new TChain("STNTUPLE");

		//PHOTON MC FILES
		chain->Add("/data/nbay02/a/samantha/STNTUPLES/MC/Pho/pd028ec7.0006exo8");
		

		
	TStnAna *ap = new TStnAna(chain);
	
	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

	

  
	TriggerModule* trigMod;
	trigMod = new TriggerModule;
		trigMod->SetGoodRunListFile("goodrun_v23_pho_00.txt");  // to p13
		trigMod->SetGoodRunBit(1);				//1= use good run, 0=do not use it
		trigMod->SetTriggerBit(1);				// 1= require tigger, 0=NO (for MC)
		trigMod->SetMinVtx(1);
		trigMod->SetUseVtxZcut(1);				//1= require z<60 0=not cut on z
		trigMod->SetMaxVtx(1000);
		trigMod->SetSummaryStat(false);		//enable/disable endjob summary
	

	
	ap->AddModule(trigMod,1);
	ap->Run(nevets);

	
	std::string filename = "MetTestJay.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
