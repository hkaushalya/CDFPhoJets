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
#include "cp2.hh"


void cafCP2Study() 
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("cp2_time");

	TStnDataset *dsp;
	TStnAna* ap = new TStnAna();
	TStnCatalog* c = new TStnCatalog();
	dsp = new TStnDataset(); //p18-20
	//267616
	//267718 this seems HV off run
	//c->InitDataset(dsp,"cdfpstn","gjt4bm","","",261119,267718);
	//c->InitDataset(dsp,"cdfpstn","gjt4bm","","",267616,267616);
	c->InitDataset(dsp,"cdfpstn","gjt4bm","","",267718,267718);
	//c->InitDataset(dsp,"cdfpstn","cph1ah","","",190851,203779);
	//0m (18-20) (d,"cdfpstn","gjt4bm","","",261119,267718)     complete
	
	//ap->SetSplit(ind,split);
	ap->AddDataset(dsp);
	ap->Print();

	TH1::AddDirectory(0);
	gSystem->CompileMacro("cp2.cc","=");
	cp2 *mycp2 = new cp2;
	mycp2->SetNewDataFileName("cp2traler-test.txt");
	mycp2->SetLastDBDataFileName("cp2traler-256904_LastDBEntry.txt");

	ap->AddModule(mycp2,1);
	
	ap->Run();
  	
	
	std::string filename = "CP2Study_267718.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	gBenchmark->Show("cp2_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
