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
#include "samantha/samantha/Pho/PhoJets.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"
#include "samantha/samantha/Pho/JetFilterModuleV3.hh"
#include "Stntuple/Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/Stntuple/loop/TStnOutputModule.hh"
#include "Stntuple/base/base/TStnDataset.hh"
#include "TSystem.h"

/////////////////////////////////////////////////////////////

void runPho2Jets(const int nEvts = 100, const int cafind=0, const int NcafSections=20) 
{
	std::cout << "Nevts to run = " << nEvts << std::endl;
	std::cout << "cafind/total sections  = " << cafind << " / " << NcafSections << std::endl;

	TStnAna*     ap;

	TStnCatalog* c = new TStnCatalog();
	TStnDataset *dsp = new TStnDataset();

	//std::cout << "Using gen 6 pho MC SAMPLE " << std::endl;
	//dsp = new TStnDataset("cdfpstn","gq0sqd");  //PYTHIA PHO+JETS+MINBIAS
	//c->InitDataset(dsp);
	//ap->AddDataset(dsp);

	TChain *chain = new TChain("STNTUPLE");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd025748.0006exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027982.0003exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027a5c.0001exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027ad7.0005exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027b7e.0001exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027c04.0009exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027cc4.0001exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027e39.0001exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd027e67.000aexo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02829f.0004exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028399.0002exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028501.0009exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd0285b1.0007exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02860c.0007exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028624.0004exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd0287f1.000cexo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028ad6.000bexo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028c17.0001exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd028ec7.0006exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02ac2e.0004exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b43e.0008exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b473.0006exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b4da.0001exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b553.0001exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02b5dc.0001exo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02ba61.001cexo8");
chain->Add("/data/nbay02/b/samantha/STNTUPLES/MC/Pho/pd02bae9.0004exo8");
	
	
	ap = new TStnAna(chain);
	
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
		jf->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf->SetMaxNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		jf->SetMaxJetEta(3.2);
		jf->SetRemoveTightPhotons(1);
		jf->SetPrintLevel(1);

	
	PhoJets *phoSel = new PhoJets();
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
	ap->AddModule(phoSel,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened

	if (nEvts>0) ap->Run(nEvts); 
	else ap->Run(); 

	std::stringstream file;
	file << "PhoJets_" << cafind << "Of" << NcafSections << ".root";
	ap->SaveHist(file.str().c_str(),2);
	std::cout << "Histograms save to file " << file.str() << std::endl;

}

