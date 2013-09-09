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
#include "TObjArray.h"
#include "TIterator.h"
#include "TFile.h"

void ListDataFiles(int dataset=11)
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnDataset *dsp;

	if (dataset == 1) dsp= new TStnDataset(); //p1-4
	if (dataset == 5)	dsp= new TStnDataset("cdfpstn","cph1ai"); //p5-10
	if (dataset	 == 11) dsp= new TStnDataset(); //p11-13

	TStnAna* ap = new TStnAna();
	
	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened

	TStnCatalog* c = new TStnCatalog();
	TObjArray* obj;
	
	
	if (dataset == 1) {
		c->InitDataset(dsp,"cdfpstn","cph1ah","","",190851,203779);
		std::cout << "files attached for cph1ah = " << dsp->GetNFiles() << std::endl;
		dsp->Print();
	}
	if (dataset == 5) {
		c->InitDataset(dsp);
		std::cout << "files attached for cph1ai = " << dsp->GetNFiles() << std::endl;
		dsp->Print();
	}

	if (dataset == 11) {
		c->InitDataset(dsp,"cdfpstn","cph1aj","","",233133,246231);
		std::cout << "files attached for cph1aj = " << dsp->GetNFiles() << std::endl;
		dsp->Print();
		//ap->AddDataset(dsp);
	}
	
	//TestModule *test = new TestModule;
	
	//ap->AddModule(test,1);

	//ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	//ap->Run();
  	
	
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
