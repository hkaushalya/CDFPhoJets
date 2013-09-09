#include "TBenchmark.h"
#include "TSystem.h"
#include <iostream>
#include <sstream>
#include <string>
#include "TObjArray.h"
#include "TIterator.h"
#include "Stntuple/Stntuple/loop/TStnAna.hh"
#include "Stntuple/Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/base/base/TStnDataset.hh"
#include "samantha/samantha/Pho/JetSmearStudy.hh"


void JetSmear(const double evts,unsigned int npts=10) 
{
	TChain* chain = new TChain("STNTUPLE");
	chain->Add("/data/nbay02/a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/ch031a01.0247ph1a");
	TStnAna *ap = new TStnAna(chain);
	
	TH1::AddDirectory(kFALSE); //do not add these histos to memory 
	
	JetSmearStudy *js = new JetSmearStudy();
		js->NptsToGenerate(npts);

	ap->AddModule(js,1);
	ap->Run(evts);

	std::string filename = "MetTest_PhotonMC.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
}
