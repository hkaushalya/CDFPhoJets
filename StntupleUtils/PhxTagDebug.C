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
#include "pdftools/TPDFReweight.hh"


void PhxTagDebug(const int MinNjet15=1, 
						const int iPhoType=0, const int debug = 100) 
{
	std::cout << "*************************************************************" << std::endl;
	std::cout << "******** Input parameters for " << __FUNCTION__ << "*********" << std::endl;
	std::cout << "******** MinNjet15 = " << MinNjet15 << std::endl;
	std::cout << "******** iPhoType  = " << iPhoType << std::endl;
	std::cout << "******** debug     = " << debug << std::endl;
	std::cout << "*************************************************************" << std::endl;

	
	assert (MinNjet15 >= 0 && "MinNjet Cannot be less than 0");
	
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnAna* ap = NULL;
	
	TChain* chain = new TChain("STNTUPLE");
	chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b4da.0001exo8");
		
	
	ap = new TStnAna(chain);
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


	int AnalysisMode=2;
	
	
	JetFilterModuleV2 *jf = new JetFilterModuleV2();  //---------- My Vertex Filter Initialization
   jf->SetDebug(0);
		jf->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		jf->SetAnalysisMode(AnalysisMode); 


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


	TPDFReweight *pdf_reweight_all = new TPDFReweight();
	pdf_reweight_all->SetMCGenerator("Pythia");
	pdf_reweight_all->SetReferencePDF(4,46);
	pdf_reweight_all->SetReweightedwithPDF("CTEQ6M");

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
	ap->AddModule(pdf_reweight_all,1);
	ap->AddModule(jf,1);
	ap->AddModule(phoSel,1);
	
	//ap->AddModule(eleSel,1);

	if (debug) ap->Run(debug);
	else ap->Run();
	
	std::stringstream filename;
	filename << "PhxDebug.root"; 
	ap->SaveHist(filename.str().c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
