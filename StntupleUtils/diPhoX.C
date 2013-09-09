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
#include "samantha/samantha/Pho/DiPhoX.hh"


void diPhoX(const int nevts = -1, const int index=1, const int MinNjet15=0, const int iPhoType=0, const int SampleId = 0, const float metSigCut=0) 
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnAna* ap = NULL;
	
	TChain* chain = new TChain("STNTUPLE");

	std::cout << "<<< Adding data files >>> " << std::endl;
	int iMin,iMax;
	if (index==0) { iMin=0;iMax=105;}
	if (index==1) { iMin=0;iMax=10;}
	if (index==2) { iMin=10;iMax=20;}
	if (index==3) { iMin=20;iMax=30;}
	if (index==4) { iMin=30;iMax=40;}
	if (index==5) { iMin=40;iMax=50;}
	if (index==6) { iMin=50;iMax=60;}
	if (index==7) { iMin=60;iMax=70;}
	if (index==8) { iMin=70;iMax=80;}
	if (index==9) { iMin=80;iMax=90;}
	if (index==10) { iMin=90;iMax=100;}
	if (index==11) { iMin=100;iMax=105;}

	std::cout << "Index, iMin,iMax=" << index << ", " << iMin << ", " << iMax << std::endl;
	for (int i=iMin; i<iMax; ++i)
	{
		std::stringstream file;
		std::string node;
		if (i<50)
		{
			node = "europium";
		} else {
			node = "praseodymium";
		}
		if (i<10)
		{
			//file << "/data/nbay02/b/samantha/STNTUPLES2/Pythia_dipho/Pythia_ggXsignal_GMSBtemplate_010209.root.0" << i;
			//file << "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/STNTUPLES2/Pythia_dipho/Pythia_ggXsignal_GMSBtemplate_010209.root.0" << i;
			file << "root://" << node << ".fnal.gov:5151//cdf/scratch/samantha/Pythia_dipho/Pythia_ggXsignal_GMSBtemplate_010209.root.0" << i;
		} else {
			//file << "/data/nbay02/b/samantha/STNTUPLES2/Pythia_dipho/Pythia_ggXsignal_GMSBtemplate_010209.root." << i;
			//file << "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/STNTUPLES2/Pythia_dipho/Pythia_ggXsignal_GMSBtemplate_010209.root." << i;
			file << "root://" << node << ".fnal.gov:5151//cdf/scratch/samantha/Pythia_dipho/Pythia_ggXsignal_GMSBtemplate_010209.root." << i;
		}
			
		chain->Add(file.str().c_str());
		std::cout << i << "\t" << file.str() << std::endl;
		
	}
	
	ap = new TStnAna(chain);
	
	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

  /******************************************************/
  // define all module here
  /******************************************************/

	bool bSummaryStat = 0; 	// 1 = no summary
  
  
  //PassFailTree *passFailFile = new PassFailTree;
  
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
	TagLoosePhotons* tagLoose = new TagLoosePhotons;
		tagLoose->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
	
	TagElectrons* tagTightEle = new TagElectrons;
		tagTightEle->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary
	TagLooseElectrons* tagLooseEle = new TagLooseElectrons;
		tagLooseEle->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary

	PhoEleFilter *peFilter = new PhoEleFilter;
		peFilter->SetSummaryStat(bSummaryStat);		//enable/disable endjob summary


	//these are common parameters for central/JES UP/JES DOWN jet collections
	//int JTC_imode = 0;		// 0==MC, 1==data
	//if (debug) std::cout  << " WARNING! RUNNING ON LOCAL DATASET. JTC_imode is default value (0 = MC)!" << std::endl;
	//if (dataset==10)  JTC_imode = 1;		// 0==MC, 1==data
	int UnclParamSwitch = 0; 	// 0==photon parametrization; 1==Zee parametrization 
	//int SampleId = 3;		// >=0 , <=4, 3- use Pythia pho+jet parameterization
	//int SampleId = 10;		// calibration mode
	//int SampleId = 3;		// 3- use Pythia for MC OR data signal pho+jet parameterization
								// 4 = DATA sideband pho+jet parameterization

	int NpointsToGenerate = 1;		//need at least 1 for  pseudo experitments
	int SelectSigMetEvent = 0; // =0 do nothing; =1 select events with MetSig>metsigCut
	//int SelectSigMetEvent = 0; // for calibration
	int AnalysisMode=1;		//0==off 1==on
	//int AnalysisMode=0;		// for calibration
	int SelectMetEvent = 0;
	int SelectExoticEvent = 0;
	int RemoveDuplicate = 1;
	int RemoveBadMet = 1;// -1=select events with bad met; 1=remove events with bad met
	int DoVxSwap = 0; // 0=do nothing; 1=swap vertices to minimize MET
	//float MaxDeltaPhiJMet = 0.3;
	int DumpEvent = 0;

	
	JetFilterModuleV2 *jf = new JetFilterModuleV2();  //---------- My Vertex Filter Initialization
   jf->SetDebug(0);
		//jf->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf->SetMaxNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf->SetJTC_imode(JTC_imode); 			// 0==MC; 1==data
		jf->SetUnclParamSwitch(UnclParamSwitch); 	// 0==photon parametrization; 1==Zee parametrization 
//		jf->SetMaxDeltaPhiJMet(MaxDeltaPhiJMet); 	// controls dPhi in MyMetCleanUpCut
		jf->SetNpointsToGenerate(NpointsToGenerate); 	// number of pseudo-experiments per each event
		jf->SetSelectMetEvent(SelectMetEvent); 		// 0==do nothing; 1==select event with large MET
		jf->SetRemoveDuplicate(RemoveDuplicate); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		jf->SetRemoveBadMet(RemoveBadMet); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		///jf->SetUseMetPDFscenario(1);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		jf->SetDumpEvent(DumpEvent); 			// 0==do nothing; 1==dump event
		jf->SetDumpEventFileName("V2DATA1_dumpEvent.dat");
		jf->SetLargeMetEventFileName("V2DATA1_largeMetEvent.dat");
		jf->SetDatFileName("V2DATA1_jet.dat");
		
		//new stuff for v2
		jf->SetSampleID(SampleId); // 0==signal
		//jf->SetMinDeltaPhiJMet(-10.0); // every event should pass this cut
		//-----------------------------
		jf->SetSelectSigMetEvent(SelectSigMetEvent); // =0 do nothing; =1 select events with MetSig>metsigCut
		jf->SetMetSigCut(metSigCut);
		jf->SetAnalysisMode(AnalysisMode);      // =0 doesn't fill AnalysisHisto, but fills MetStudyHisto
												// =1 fills AnalysisHisto, but doesn't fill MetStudyHisto
		//-----------------------------
		jf->SetSelectExoticEvent(SelectExoticEvent);		// i have commented out these parts 
		jf->SetDoVxSwap(DoVxSwap); // 0=do nothing; 1=swap vertices to minimize MET
		jf->SetUseHadJets(0);		//0=use det jets 1=use had jets
		jf->SetDoubleSmearJets(0);		//0=smear once jets 1=smear twice
		//jf->SetMaxJetEta(3.6);   //Sasha wanted this for Fix3_2
		jf->SetZeroJERoffset(0);		//1=no offset ,0=use offset
		jf->SetMetCut(0);
		jf->SetPrintLevel(10);		//higher number => more prints
		jf->SetMEtAddMethod(2);		// how we add MEt to a jet 0=not added, 1=most mismeasure jet
											// 2= closeset jet to MEt

		jf->SetRemoveTightPhotons(1);

	JetFilterModuleV2 *jf2 = new JetFilterModuleV2("JetFilterV2_2");  //---------- My Vertex Filter Initialization
   jf2->SetDebug(0);
		//jf2->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf2->SetMaxNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		//jf2->SetJTC_imode(JTC_imode); 			// 0==MC; 1==data
		jf2->SetUnclParamSwitch(UnclParamSwitch); 	// 0==photon parametrization; 1==Zee parametrization 
//		jf2->SetMaxDeltaPhiJMet(MaxDeltaPhiJMet); 	// controls dPhi in MyMetCleanUpCut
		jf2->SetNpointsToGenerate(NpointsToGenerate); 	// number of pseudo-experiments per each event
		jf2->SetSelectMetEvent(SelectMetEvent); 		// 0==do nothing; 1==select event with large MET
		jf2->SetRemoveDuplicate(RemoveDuplicate); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		jf2->SetRemoveBadMet(RemoveBadMet); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		///jf2->SetUseMetPDFscenario(1);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		jf2->SetDumpEvent(DumpEvent); 			// 0==do nothing; 1==dump event
		jf2->SetDumpEventFileName("V2DATA1_dumpEvent.dat");
		jf2->SetLargeMetEventFileName("V2DATA1_largeMetEvent.dat");
		jf2->SetDatFileName("V2DATA1_jet.dat");
		
		//new stuff for v2
		jf2->SetSampleID(SampleId); // 0==signal
		//jf2->SetMinDeltaPhiJMet(-10.0); // every event should pass this cut
		//-----------------------------
		jf2->SetSelectSigMetEvent(SelectSigMetEvent); // =0 do nothing; =1 select events with MetSig>metsigCut
		jf2->SetMetSigCut(metSigCut+2);
		jf2->SetAnalysisMode(AnalysisMode);      // =0 doesn't fill AnalysisHisto, but fills MetStudyHisto
												// =1 fills AnalysisHisto, but doesn't fill MetStudyHisto
		//-----------------------------
		jf2->SetSelectExoticEvent(SelectExoticEvent);		// i have commented out these parts 
		jf2->SetDoVxSwap(DoVxSwap); // 0=do nothing; 1=swap vertices to minimize MET
		jf2->SetUseHadJets(0);		//0=use det jets 1=use had jets
		jf2->SetDoubleSmearJets(0);		//0=smear once jets 1=smear twice
		//jf2->SetMaxJetEta(3.6);   //Sasha wanted this for Fix3_2
		jf2->SetZeroJERoffset(0);		//1=no offset ,0=use offset
		jf2->SetMetCut(0);
		jf2->SetPrintLevel(10);		//higher number => more prints
		jf2->SetMEtAddMethod(2);		// how we add MEt to a jet 0=not added, 1=most mismeasure jet
											// 2= closeset jet to MEt

		jf2->SetRemoveTightPhotons(1);




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
		
/*	Pho2JetsTemp *phoSel = new Pho2JetsTemp();
		phoSel->SetPhoMinDetEta(0);
		phoSel->SetPhoMaxDetEta(0.95);
		phoSel->SetPhoMinEt(30.0);
		phoSel->SetHaloType(5);
		phoSel->SetPhotonType(iPhoType);		//0=signal, 1=sideband
		phoSel->SetSummaryStat(bSummaryStat);		//1=disable 0=enable endjob summary
*/
	DiPhoX *phoSel = new DiPhoX();
		phoSel->SetPhoMinDetEta(0);
		phoSel->SetPhoMaxDetEta(0.9);
		phoSel->SetPhoMinEt(13.0);
		phoSel->SetHaloType(5);
		phoSel->SetPhotonType(iPhoType);		//0=signal, 1=sideband
		phoSel->SetSummaryStat(bSummaryStat);		//1=disable 0=enable endjob summary
		phoSel->SetManualMCflag(1);			//set the MC flag manually
														// so the PMT/BH mod requirement will be ignored

	
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(tagPhoenix,1);
	ap->AddModule(peFilter,1);
	//ap->AddModule(tagPMT,1);
	//ap->AddModule(tagBH,1);
	//ap->AddModule(EmT,1);
	//ap->AddModule(setEmTime,1);
	ap->AddModule(tagConvEle,1);
	ap->AddModule(phoSel,1);
	
	ap->AddModule(jf,1);
	ap->AddModule(jf2,1);
	
	if (nevts>0)
		ap->Run(nevts);
	else 
		ap->Run();


	
	std::stringstream filename;
	filename << "DiPhoX.root_" << index;
	//std::string filename = "DiPhoX_PhotonMC.root";
	//ap->SaveHist(filename.c_str(),2);
	ap->SaveHist(filename.str().c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename.str() << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
