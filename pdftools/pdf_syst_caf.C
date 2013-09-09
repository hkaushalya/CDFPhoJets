#include "TSystem.h"
#include "TChain.h"
#include "Stntuple/Stntuple/loop/TStnAna.hh"
#include <iostream>
#include <sstream>
#include <string>
#include "./TPDFReweight.hh"
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
#include "samantha/samantha/Pho/Pho2JetsTemp.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"
#include "samantha/samantha/Pho/JetFilterModuleV2.hh"
#include "Stntuple/Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/Stntuple/loop/TStnOutputModule.hh"
#include "Stntuple/base/base/TStnDataset.hh"
#include "TSystem.h"

/////////////////////////////////////////////////////////////
// set the environment varible to find the 'tables'
// "setenv PDFTOOLS_DIR ./pdftools/"
// then compile this:(first compile and load TPDFReweight.cc)
// -- now the above env var is set within this macro, when this
// is compiled. -01-11-2009
/////////////////////////////////////////////////////////////
/***********************************************************
 * $Id: pdf_syst_caf.C,v 1.1 2013/02/28 03:42:41 samantha Exp $
 * $Log: pdf_syst_caf.C,v $
 * Revision 1.1  2013/02/28 03:42:41  samantha
 * Final commit. no checks.! these were never commited to cvs before!
 *
 *
 *
 ***********************************************************/

void pdf_syst_caf(const int nEvts = 100, const int pdfset = 0,const int cafind=0, const int NcafSections=20) 
{
	std::cout << "Nevts to run = " << nEvts << std::endl;
	std::cout << "pdfset      = " << pdfset  << " (0=CTEQ6M, 1=MRST98LO)" << std::endl;
	std::cout << "cafind/total sections  = " << cafind << " / " << NcafSections << std::endl;

	std::stringstream sPdfDir;
	sPdfDir << gSystem->Getenv("PDFTOOLS_DIR");
	if (sPdfDir.str().length() == 0)
	{
		std::cout << __FILE__ <<":" <<  __FUNCTION__ << "::Required environment variable PDFTOOLS_DIR is not specified! Assigning the default " << std::endl;
		std::string sDefPdfDir("./samantha/pdftools/");
		gSystem->Setenv("PDFTOOLS_DIR",sDefPdfDir.c_str());	
		std::cout << "Using PDFTOOLS_DIR = " << gSystem->Getenv("PDFTOOLS_DIR") << std::endl;
	}

	gSystem->Load("./samantha/pdftools/libpdftools.so");
	gSystem->CompileMacro("./samantha/pdftools/TPDFReweight.cc","k");

	TStnAna*     ap;

	TStnCatalog* c = new TStnCatalog();
	TStnDataset *dsp = new TStnDataset();

	std::cout << "Using gen 6 pho MC SAMPLE " << std::endl;
	dsp = new TStnDataset("cdfpstn","gq0sqd");  //PYTHIA PHO+JETS+MINBIAS
	c->InitDataset(dsp);

	ap = new TStnAna();
	ap->AddDataset(dsp);
	ap->SetSplit(cafind, NcafSections);


	bool bSummaryStat = 0; 	// 1 = no summary


	TPDFReweight *pdf_reweight_all = new TPDFReweight();

	pdf_reweight_all->SetMCGenerator("Pythia");
	pdf_reweight_all->SetReferencePDF(4,46);
	if (pdfset == 0) pdf_reweight_all->SetReweightedwithPDF("CTEQ6M");
	else if (pdfset == 1) pdf_reweight_all->SetReweightedwithPDF("MRST98LO");
	//pdf_reweight_all->SetReweightedwithPDF("CTEQ5L");
	//pdf_reweight_all->SetMCScale(228.);  //for resonance production only
	pdf_reweight_all->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary


	TriggerModule* trigMod;
	trigMod = new TriggerModule;
	trigMod->SetGoodRunListFile("goodrun_v31_pho_00.txt");  // to p31
	//trigMod->SetGoodRunListFile("goodrun_v23_pho_00.txt");  // to p23
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
	JetFilterModuleV2 *jf = new JetFilterModuleV2();  //---------- My Vertex Filter Initialization
   jf->SetDebug(0);
		jf->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf->SetMaxNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		jf->SetAnalysisMode(2);      // =2 no anaysis done. only producing corrected jets
		jf->SetMaxJetEta(3.0);   //Sasha wanted this for Fix3_2
		jf->SetRemoveTightPhotons(1);
		jf->SetPrintLevel(1);

	
	Pho2JetsTemp *phoSel = new Pho2JetsTemp();
	phoSel->SetPhoMinDetEta(0);
	phoSel->SetPhoMaxDetEta(1.1);
	phoSel->SetPhoMinEt(30.0);
	phoSel->SetHaloType(5);
	phoSel->SetPhotonType(0);		//0=signal, 1=sideband
	phoSel->SetDoPDFSyst(1);
	phoSel->SetSummaryStat(bSummaryStat);		//1=disable 0=enable endjob summary

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
	ap->AddModule(pdf_reweight_all,1);
	ap->AddModule(phoSel,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened

	if (nEvts>0) ap->Run(nEvts); 
	else ap->Run(); 

	std::stringstream file;
	if (pdfset == 0) file << "PDFsyst_"<< "CTEQ6M_" << cafind << ".root";
	else if (pdfset == 1) file << "PDFsyst_"<< "MRST98LO_" << cafind << ".root";
	ap->SaveHist(file.str().c_str(),2);
	std::cout << "Histograms save to file " << file.str() << std::endl;

}

