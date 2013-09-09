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
 * $Id: pdf_syst.C,v 1.5 2011/05/26 16:47:15 samantha Exp $
 *
 * $Log: pdf_syst.C,v $
 * Revision 1.5  2011/05/26 16:47:15  samantha
 * Minor changes for new jobs. Added a new dataset, 10.
 *
 * Revision 1.4  2010/02/04 18:12:07  samantha
 * Modified to break the dataset into several job segments.
 *
 *
 ***********************************************************/

void pdf_syst(const int nEvts = 100, int fileset=1, const int local=1, const int cafind=0) 
{
	std::cout << "Nevts to run = " << nEvts << std::endl;
	std::cout << "fileset      = " << fileset << std::endl;
	std::cout << "Local?       = " << local << std::endl;

	std::stringstream sPdfDir;
	sPdfDir << gSystem->Getenv("PDFTOOLS_DIR");
	if (sPdfDir.str().length() == 0)
	{
		std::cout << __FILE__ <<":" <<  __FUNCTION__ << "::Required environment variable PDFTOOLS_DIR is not specified! Assigning the default " << std::endl;
		std::string sDefPdfDir("./samantha/pdftools/");
		gSystem->Setenv("PDFTOOLS_DIR",sDefPdfDir.c_str());	
		std::cout << "Using PDFTOOLS_DIR = " << gSystem->Getenv("PDFTOOLS_DIR") << std::endl;
	}

	if (local==0)
	{
		std::stringstream temp;
		temp << gSystem->Getenv("DATASET");
		if (temp.str().length())
		{
			std::string sd(temp.str());
			std::cout << "sd = " << sd << std::endl;
			if (sd == "1") fileset = 1;
			if (sd == "2") fileset = 2;
			if (sd == "3") fileset = 3;
			if (sd == "4") fileset = 4;
			if (sd == "5") fileset = 5;
			if (sd == "6") fileset = 6;
			if (sd == "10") fileset = 10;
			std::cout << "file set from env variable: fileset " << fileset << std::endl;
		} else 
		{
			std::cout << __FILE__ <<":" <<  __FUNCTION__ 
				<< "::Required environment variable DATASET is not specified for remote data! exiting!"
				<< std::endl;
			exit (1);
		}
	}
	
	gSystem->Load("./samantha/pdftools/libpdftools.so");
	gSystem->CompileMacro("./samantha/pdftools/TPDFReweight.cc","k");

	TStnAna*     ap;

	// All signal events //split the job to run over smaller
	if (local)
	{
		TChain*      inChain = new TChain("STNTUPLE");
		if (fileset == 1 || fileset ==0)
		{
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b4da.0001exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028ad6.000bexo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028624.0004exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02860c.0007exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd0285b1.0007exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028501.0009exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd0287f1.000cexo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02829f.0004exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027e67.000aexo8");
		} else if (fileset == 2 || fileset ==0)
		{
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027e39.0001exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027cc4.0001exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028399.0002exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027b7e.0001exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd025748.0006exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027982.0003exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027a5c.0001exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027ad7.0005exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027c04.0009exo8");
		} else if (fileset == 3 || fileset ==0)
		{
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028c17.0001exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028ec7.0006exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02ac2e.0004exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b43e.0008exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b473.0006exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b553.0001exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b5dc.0001exo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02ba61.001cexo8");
			inChain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02bae9.0004exo8");
		}

		ap = new TStnAna(inChain);
		
	} else 
	{
			TStnCatalog* c = new TStnCatalog();
			TStnDataset *dsp = new TStnDataset();
			
			//old GEN5 ZERO min bias MC sample 141544,179056
			if (fileset == 1)
			{
				std::cout << "Using gen 6 pho MC SAMPLE " << std::endl;
				c->InitDataset(dsp,"stntuple/dev_242","pypho22_dfc","","",141597,148000);
			} 
			//else if (fileset == 2) c->InitDataset(dsp,"stntuple/dev_242","pypho22_dfc","","",148001,154000);
			//else if (fileset == 3) c->InitDataset(dsp,"stntuple/dev_242","pypho22_dfc","","",154001,160000);
			//else if (fileset == 4) c->InitDataset(dsp,"stntuple/dev_242","pypho22_dfc","","",160001,166000);
			//else if (fileset == 5) c->InitDataset(dsp,"stntuple/dev_242","pypho22_dfc","","",166001,172000);
			//else if (fileset == 6) c->InitDataset(dsp,"stntuple/dev_242","pypho22_dfc","","",172001,179056);
			else if (fileset == 10) { //gen 6 MC sample with minbias
				std::cout << "Using gen 6 pho MC SAMPLE " << std::endl;
				dsp = new TStnDataset("cdfpstn","gq0sqd");  //PYTHIA PHO+JETS+MINBIAS
				c->InitDataset(dsp);
			}
			
		ap = new TStnAna();
		ap->AddDataset(dsp);
		ap->SetSplit(cafind, 6);

	}

	bool bSummaryStat = 0; 	// 1 = no summary


	TPDFReweight *pdf_reweight_all = new TPDFReweight();

	pdf_reweight_all->SetMCGenerator("Pythia");
	pdf_reweight_all->SetReferencePDF(4,46);
	pdf_reweight_all->SetReweightedwithPDF("CTEQ6M");
	//pdf_reweight_all->SetReweightedwithPDF("MRST98LO");
	//pdf_reweight_all->SetReweightedwithPDF("CTEQ5L");
	//pdf_reweight_all->SetMCScale(228.);  //for resonance production only
	pdf_reweight_all->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary


	TriggerModule* trigMod;
	trigMod = new TriggerModule;
	//trigMod->SetGoodRunListFile("goodrun_v31_pho_00.txt");  // to p31
	trigMod->SetGoodRunListFile("goodrun_v23_pho_00.txt");  // to p23
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
	file << "PDFsyst"<< fileset << ".root";
	ap->SaveHist(file.str().c_str(),2);
	std::cout << "Histograms save to file " << file.str() << std::endl;

}

