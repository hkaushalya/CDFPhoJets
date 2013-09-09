#include "TBenchmark.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TIterator.h"
#include <iostream>
#include <sstream>
#include <string>
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
#include "samantha/samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/samantha/Pho/TagPMTSpikes.hh"
#include "samantha/samantha/Pho/SetEmTimes.hh"
#include "samantha/samantha/utils/PassFailTree.hh"
#include "samantha/samantha/utils/BadPointerChecker.hh"
#include "samantha/samantha/Pho/TagCosmicPhotons.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/FlatStupleMaker.hh"
#include "samantha/samantha/Pho/TagConvElectrons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"
#include "samantha/samantha/Pho/Pho2JetsTemp.hh"
#include "samantha/samantha/Pho/EventProperties.hh"
#include "samantha/samantha/Pho/JetFilterModuleV3.hh"

void TestFlatNtupleMaker(unsigned runEvts=60000, std::string dataset="phodata")
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnAna* ap;
	TChain *chain = 0;

	//TStnDataset *dsp = new TStnDataset();
	//TStnCatalog* c = new TStnCatalog();
	//		c->InitDataset(dsp,"stntuple/dev_242","pypho22_dfc","","");
	//		ap->AddDataset(dsp);

		/* TO TEST USE A LOCAL FILE */ 
		chain = new TChain("STNTUPLE");

		  		//chain->Add("/data/nbay02/b/samantha/STNTUPLES/DATA/PHOTON-cph1ah/ch031a01.0247ph1a");
		  		//chain->Add("/data/nbay02/b/samantha/STNTUPLES/DATA/PHOTON-cph1ah/ch031a53.020fph1a");
		  		//chain->Add("/data/nbay02/b/samantha/STNTUPLES/DATA/PHOTON-cph1ah/ch031b6a.00bbph1a");
		  		//chain->Add("/data/nbay02/b/samantha/STNTUPLES/DATA/PHOTON-cph1ah/ch031bf4.005eph1a");
		  		//chain->Add("/data/nbay02/b/samantha/STNTUPLES/DATA/PHOTON-cph1ah/ch02ec50.01feph1a");
		  	
				//DI-PHO MC
				//chain->Add("~/STNTUPLES/MC/DIPHO/gg038b5e.0019x0s1");

		//susy
		chain->Add("/data/nbay04/c/samantha/susy/0/stntuple.root");

		ap = new TStnAna(chain);

	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

	//******************************************//
	//for stripping evts
	//----  need this part for event preselection only
   //TStnOutputModule* out= new TStnOutputModule("MetModelTest_EMobj_gt_1.root");
	//out->SetMaxFileSize(300);
	//ap->SetOutputModule(out);
	
	//******************************************//		 

  /******************************************************/
  // define all module here
  /******************************************************/
	bool bSummaryStat = 0; 	// 1 = no summary
  
  PassFailTree *passFailFile = new PassFailTree;
  BadPointerChecker *ptChecker = new BadPointerChecker;
  
	TriggerModule* trigMod;
	trigMod = new TriggerModule;
		trigMod->SetGoodRunListFile("goodrun_v31_pho_00.txt");  // to p25
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
											//1,2 - pass only if a electron found and remove others from list

	PhoEleFilter *peFilter = new PhoEleFilter;
		peFilter->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
		peFilter->SetMinEmEt(25.);
		peFilter->SetMaxEmDetEta(1.1);

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
		//setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt"); //does not matter to me. these are small corrections
		setEmTime->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

/*	Pho2JetsTemp *phoSel = new Pho2JetsTemp();
		phoSel->SetPhoMinDetEta(0);
		phoSel->SetPhoMaxDetEta(0.95);
		phoSel->SetPhoMinEt(30.0);
		phoSel->SetHaloType(5);
		phoSel->SetPhotonType(0);		//0=signal, 1=sideband
		phoSel->SetSummaryStat(bSummaryStat);		//1=disable 0=enable endjob summary
*/
																									  
	//these are common parameters for central/JES UP/JES DOWN jet collections
	const int MinNjet15 = 0;		//use this always. 0 jets will not give you correct met/ht. as there is some hard coded Et cuts for those!
	const float fMaxJetEta=3.2;
	const int AnalysisMode=2;		//0==off 1==on, >1 = generate corrected jets and MET/SUMET
	const float fMinEmEtThr = 15.0; //em objects with Et> this will be removed and used in Ht calculation
	//do not lower this or pass EM objects wit hEt below this to be removed from jet list.
	//this will cause the Ht/met calculation to change unpredictably!
	const bool bRemAllEm = true;
	const int iMetScenario =3; // 0--raw Met; corrected for jets with: 1-- et>5, 2-- et>10, 3-- et>15, 4-- et>20
	//this does not matter in ntuple making. but use '3' which is Et>15GeV objects for Ht and MET
	//when using ht/met in final plots
	
/*	JetFilterModuleV2 *jf = new JetFilterModuleV2();  //---------- My Vertex Filter Initialization
   jf->SetDebug(0);
		jf->SetJTC_systcode(0);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf->SetMaxJetEta(fMaxJetEta);
		jf->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		jf->SetMinEtThr(fMinEmEtThr);		//Et>15.0 is the default for min Et for all em objects to be removed from jet list
		jf->SetRemoveTightPhotons(bRemAllEm);
		jf->SetRemoveLoosePhotons(bRemAllEm);
		jf->SetRemoveTightPhoLikeElectrons(bRemAllEm);
		jf->SetRemoveLoosePhoLikeElectrons(bRemAllEm);
		jf->SetRemoveStdLooseElectrons(bRemAllEm);
		jf->SetAnalysisMode(AnalysisMode);      // =0 doesn't fill AnalysisHisto, but fills MetStudyHisto
												// =1 fills AnalysisHisto, but doesn't fill MetStudyHisto
		jf->SetPrintLevel(0);		//higher number => more prints
		jf->SetSummaryStat(bSummaryStat);
		jf->SetMyMetScenario(iMetScenario);


	JetFilterModuleV2 *jf_jesup = new JetFilterModuleV2("JetFilterV2JESup","JetFilterV2JESup");
		jf_jesup->SetJTC_systcode(1);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf_jesup->SetMaxJetEta(fMaxJetEta);
		jf_jesup->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		jf_jesup->SetMinEtThr(fMinEmEtThr);		//Et>15.0 is the default for min Et for all em objects to be removed from jet list
		jf_jesup->SetRemoveTightPhotons(bRemAllEm);
		jf_jesup->SetRemoveLoosePhotons(bRemAllEm);
		jf_jesup->SetRemoveTightPhoLikeElectrons(bRemAllEm);
		jf_jesup->SetRemoveLoosePhoLikeElectrons(bRemAllEm);
		jf_jesup->SetRemoveStdLooseElectrons(bRemAllEm);
		jf_jesup->SetAnalysisMode(AnalysisMode);      // =0 doesn't fill AnalysisHisto, but fills MetStudyHisto
												// =1 fills AnalysisHisto, but doesn't fill MetStudyHisto
		jf_jesup->SetPrintLevel(0);		//higher number => more prints
		jf_jesup->SetSummaryStat(bSummaryStat);
		jf_jesup->SetMyMetScenario(iMetScenario);

	
	JetFilterModuleV2 *jf_jesdown = new JetFilterModuleV2("JetFilterV2JESdown","JetFilterV2JESdown");
		jf_jesdown->SetJTC_systcode(-1);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf_jesdown->SetMaxJetEta(fMaxJetEta);
		jf_jesdown->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		jf_jesdown->SetMinEtThr(fMinEmEtThr);		//Et>15.0 is the default for min Et for all em objects to be removed from jet list
		jf_jesdown->SetRemoveTightPhotons(bRemAllEm);
		jf_jesdown->SetRemoveLoosePhotons(bRemAllEm);
		jf_jesdown->SetRemoveTightPhoLikeElectrons(bRemAllEm);
		jf_jesdown->SetRemoveLoosePhoLikeElectrons(bRemAllEm);
		jf_jesdown->SetRemoveStdLooseElectrons(bRemAllEm);
		jf_jesdown->SetAnalysisMode(AnalysisMode);      // =0 doesn't fill AnalysisHisto, but fills MetStudyHisto
												// =1 fills AnalysisHisto, but doesn't fill MetStudyHisto
		jf_jesdown->SetPrintLevel(0);		//higher number => more prints
		jf_jesdown->SetSummaryStat(bSummaryStat);
		jf_jesdown->SetMyMetScenario(iMetScenario);
		*/

	
	JetFilterModuleV3 *jf3 = new JetFilterModuleV3();  //---------- My Vertex Filter Initialization
   	jf3->SetDebug(0);
		jf3->SetJTC_systcode(0);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf3->SetMaxJetEta(fMaxJetEta);
		jf3->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf3->SetMinEtThr(fMinEmEtThr);		//Et>15.0 is the default for min Et for all em objects to be removed from jet list
		jf3->SetRemoveTightPhotons(bRemAllEm);
		jf3->SetRemoveLoosePhotons(bRemAllEm);
		jf3->SetRemoveTightPhoLikeElectrons(bRemAllEm);
		jf3->SetRemoveLoosePhoLikeElectrons(bRemAllEm);
		jf3->SetRemoveStdLooseElectrons(bRemAllEm);
		jf3->SetPrintLevel(0);		//higher number => more prints
		jf3->SetSummaryStat(bSummaryStat);
		//jf3->SetMinEtThr(fMinEmEtThr);
		jf3->SetMyMetScenario(iMetScenario);


	JetFilterModuleV3 *jf3Up = new JetFilterModuleV3("JetFilterV3JESup");  //---------- My Vertex Filter Initialization
   	jf3Up->SetDebug(0);
		jf3Up->SetJTC_systcode(1);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf3Up->SetMaxJetEta(fMaxJetEta);
		jf3Up->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf3Up->SetMinEtThr(fMinEmEtThr);		//Et>15.0 is the default for min Et for all em objects to be removed from jet list
		jf3Up->SetRemoveTightPhotons(bRemAllEm);
		jf3Up->SetRemoveLoosePhotons(bRemAllEm);
		jf3Up->SetRemoveTightPhoLikeElectrons(bRemAllEm);
		jf3Up->SetRemoveLoosePhoLikeElectrons(bRemAllEm);
		jf3Up->SetRemoveStdLooseElectrons(bRemAllEm);
		jf3Up->SetPrintLevel(0);		//higher number => more prints
		jf3Up->SetSummaryStat(bSummaryStat);
		jf3Up->SetMinEtThr(fMinEmEtThr);
		jf3Up->SetMyMetScenario(iMetScenario);

	JetFilterModuleV3 *jf3Down = new JetFilterModuleV3("JetFilterV3JESdown");  //---------- My Vertex Filter Initialization
   	jf3Down->SetDebug(0);
		jf3Down->SetJTC_systcode(-1);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		jf3Down->SetMaxJetEta(fMaxJetEta);
		jf3Down->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf3Down->SetMinEtThr(fMinEmEtThr);		//Et>15.0 is the default for min Et for all em objects to be removed from jet list
		jf3Down->SetRemoveTightPhotons(bRemAllEm);
		jf3Down->SetRemoveLoosePhotons(bRemAllEm);
		jf3Down->SetRemoveTightPhoLikeElectrons(bRemAllEm);
		jf3Down->SetRemoveLoosePhoLikeElectrons(bRemAllEm);
		jf3Down->SetRemoveStdLooseElectrons(bRemAllEm);
		jf3Down->SetPrintLevel(0);		//higher number => more prints
		jf3Down->SetSummaryStat(bSummaryStat);
		//jf3Down->SetMinEtThr(fMinEmEtThr);
		jf3Down->SetMyMetScenario(iMetScenario);
		
		
	FlatStupleMaker *flatS = new FlatStupleMaker;
		const float fStupleVer = 8;
		std::stringstream stupleName;
		stupleName << "StupleV"<< fStupleVer << "_SUSY_0.root";
		flatS->SetStupleName(stupleName.str());
		flatS->SetSummaryStat(bSummaryStat);
		flatS->SetStupleVersion(fStupleVer);

	//EventProperties * evtProp = new EventProperties;

	//Pho2JetsTemp *pho2j = new Pho2JetsTemp();
	//	pho2j->SetHaloType(5);
	
	
	//ap->AddModule(passFailFile,1);
	//ap->AddModule(ptChecker,1);
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(tagPhoenix,1);
	ap->AddModule(tagPMT,1);
	ap->AddModule(tagBH,1);
	ap->AddModule(peFilter,1);
	ap->AddModule(tagConvEle,1);
	
	ap->AddModule(EmT,1);
	ap->AddModule(setEmTime,1);
	//ap->AddModule(evtProp,1);
	//ap->AddModule(phoSel,1);
	//ap->AddModule(jf,1);
	//ap->AddModule(jf_jesup,1);
	//ap->AddModule(jf_jesdown,1);
	
	ap->AddModule(jf3,1);
	ap->AddModule(jf3Up,1);
	ap->AddModule(jf3Down,1);
	ap->AddModule(flatS,1);
	//ap->AddModule(pho2j,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	ap->Run(runEvts);
	//ap->ProcessEvent(190925, 1049612);
	//ap->ProcessEvent(191568, 3290140);
	//ap->ProcessEvent(203265,4176808);
	//ap->ProcessEvent(203265,6292805);
	//ap->ProcessEvent(203265,6397380);


	


	std::string filename = "MetTest_My20Evts.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	//std::cout << "::::::::: Created Pass/Fail file => "<< passFailFile->GetFileName() << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
