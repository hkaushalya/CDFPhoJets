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
//#include "samantha/samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/samantha/Pho/TMyJetFilterModule.hh"
#include "samantha/samantha/Pho/TagPMTSpikes.hh"
#include "samantha/samantha/Pho/SetEmTimes.hh"
#include "samantha/samantha/utils/PassFailTree.hh"
#include "samantha/samantha/utils/BadPointerChecker.hh"
#include "samantha/samantha/Pho/TagCosmicPhotons.hh"
#include "samantha/samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/samantha/Pho/FlatStupleMaker.hh"
#include "samantha/samantha/Pho/TagConvElectrons.hh"
#include "samantha/samantha/Pho/PhoEleFilter.hh"
#include "samantha/samantha/Pho/CosmicStudy.hh"
#include "samantha/samantha/Pho/EventProperties.hh"

void TestMetV2(unsigned runEvts=100000, std::string dataset="phodata")
{
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");

	TStnAna* ap;
	TChain *chain = 0;

		/* TO TEST USE A LOCAL FILE */ 
		chain = new TChain("STNTUPLE");

		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_0.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_1.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_10.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_11.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_12.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_13.root");
/*		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_14.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_15.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_16.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_17.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_18.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_19.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_2.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_20.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_21.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_22.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_23.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_24.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_25.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_26.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_27.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_28.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_29.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_3.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_30.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_31.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_32.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_33.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_34.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_35.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_4.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_5.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_6.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_7.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_8.root");
		chain->Add("~/RESULTS/Cosmics/CosmicsSkim_novtxAndEmTime30to90/Cosmics_phodata_1_9.root");
		*/

/*
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.000eph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0020ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0037ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0058ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.00caph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.00ecph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.00fdph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.011fph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0163ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0187ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.01baph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.01dcph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0200ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0276ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.02b7ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.02c7ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0335ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.03dfph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0406ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0413ph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.053bph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.05ccph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.069cph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.06aeph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.06edph1a");
				chain->Add("/mnt/autofs/misc/nbay03.a/samantha/STNTUPLES/DATA/PHOTON-cph1ai/cj039630.0852ph1a");
*/

		ap = new TStnAna(chain);
		

	TH1::AddDirectory(kFALSE); //do not add these histos to memory 

  /******************************************************/
  // define all module here
  /******************************************************/
  
  //PassFailTree *passFailFile = new PassFailTree;
  BadPointerChecker *ptChecker = new BadPointerChecker;
  
	TriggerModule* trigMod;
	trigMod = new TriggerModule;
		trigMod->SetGoodRunListFile("goodrun_v19_pho_00.txt");  // to p13
		trigMod->SetGoodRunBit(1);		//1= use good run, 0=do not use it
		trigMod->SetTriggerBit(1);		// 1= require tigger, 0=NO
		trigMod->SetMinVtx(0);
		trigMod->SetMaxVtx(1000);
		trigMod->SetUseVtxZcut(0);		//1= require z<60 0=not cut on z
	
	InitSuperPhotons* myinit = new InitSuperPhotons;

	TagTightPhotons* tagTight = new TagTightPhotons;
	TagLoosePhotons* tagLoose = new TagLoosePhotons;
	
	TagElectrons* tagTightEle = new TagElectrons;
	TagLooseElectrons* tagLooseEle = new TagLooseElectrons;
											//1,2 - pass only if a electron found and remove others from list

	PhoEleFilter *peFilter = new PhoEleFilter;

	TagBeamHalo* tagBH = new TagBeamHalo();

	TagPMTSpikes* tagPMT = new TagPMTSpikes;
		tagPMT->SetMode(1);	//1= pass only if no spikes found
	

/*
	JetFilterModuleV2 *jf = new JetFilterModuleV2();  //---------- My Vertex Filter Initialization
		jf->SetMinNjet15(1);				//Jets required to pass (Et>15GeV)
		jf->SetJTC_imode(1); 			// 0==MC; 1==data
		jf->SetUnclParamSwitch(0); 	// 0==photon parametrization; 1==Zee parametrization 
		jf->SetMaxDeltaPhiJMet(0.3); 	// controls dPhi in MyMetCleanUpCut
		jf->SetNpointsToGenerate(1); 	// number of pseudo-experiments per each event
		jf->SetSelectMetEvent(0); 		// 0==do nothing; 1==select event with large MET
		jf->SetRemoveDuplicate(0); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		jf->SetRemoveBadMet(0); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		///jf->SetUseMetPDFscenario(1);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		jf->SetDumpEvent(0); 			// 0==do nothing; 1==dump event
		jf->SetDumpEventFileName("V2DATA_dumpEvent.dat");
		jf->SetLargeMetEventFileName("V2DATA_largeMetEvent.dat");
		jf->SetDatFileName("V2DATA_jet.dat");
		
		//new stuff for v2
		jf->SetSampleID(0); // 0==signal
		jf->SetMinDeltaPhiJMet(-10.0); // every event should pass this cut
		//-----------------------------
		jf->SetSelectSigMetEvent(0); // =0 do nothing; =1 select events with MetSig>metsigCut
		jf->SetMetSigCut(0);
		jf->SetAnalysisMode(1);      // =0 doesn't fill AnalysisHisto, but fills MetStudyHisto
												// =1 fills AnalysisHisto, but doesn't fill MetStudyHisto
		//-----------------------------
		jf->SetSelectExoticEvent(0);
		jf->SetDoVxSwap(0); // 0=do nothing; 1=swap vertices to minimize MET
*/
																									  
	// parameters for jet filter
	int MinNjet15 = 1;
	int JTC_imode = 1;		// 0==MC, 1==data
	int UnclParamSwitch = 0;
	float MaxDeltaPhiJMet = 0.3;
	int NpointsToGenerate = 20;
	int SelectMetEvent = 0;
	int RemoveDuplicate = 0;
	int RemoveBadMet = 0;
	int DumpEvent = 1;
	int UseMetPdf=0;  	// MetPDF scenarios: 0=no MetPDF; 
							// 1= 3 sigma (0.27%) & Npoints=370; 
							// 2= 3.29 sigma (0.1%) & Npoints=1000;
							// 3= 3.89 sigma (0.01%) & Npoints=10000;
							// 4= 4 sigma (0.0063%) & Npoints=15873;
							// 5= 2.58 sigma (1%), to calculate Mean & Sigma, Npoints=100.
										  
	//systcode shoule be 0, +6 or -6 for central, jes up and jes down collections.
	// 6 indicates the correction level

	TMyJetFilterModule *m5 = new TMyJetFilterModule();
		m5->SetMinNjet15(MinNjet15);				//Jets required to pass (Et>15GeV)
		m5->SetJTC_imode(JTC_imode); 			// 0==MC; 1==data
		m5->SetJTC_systcode(0);			// if >0 +sigma deviation, <0 -sigma deviation, all corrections are added in quadratures up to the given correction level
		m5->SetUnclParamSwitch(UnclParamSwitch); 	// 0==photon parametrization; 1==Zee parametrization 
		m5->SetMaxDeltaPhiJMet(MaxDeltaPhiJMet); 	// controls dPhi in MyMetCleanUpCut
		m5->SetNpointsToGenerate(NpointsToGenerate); 	// number of pseudo-experiments per each event
		m5->SetSelectMetEvent(SelectMetEvent); 		// 0==do nothing; 1==select event with large MET
		m5->SetRemoveDuplicate(RemoveDuplicate); 	// 0==do nothing; 1==remove "duplicate" events; -1==select "duplicate" events; "duplicate"==bad match with jets  
		m5->SetRemoveBadMet(RemoveBadMet); 		// -1=select events with bad met; 1=remove events with bad met; 0==do nothing
		//m5->SetUseMetPDFscenario(UseMetPdf);  // this over writes the Npoints generates. scenarios 1=370
												//	2=1000, 3=100000,4=15873,5=100 default=100
		m5->SetDumpEvent(DumpEvent); 			// 0==do nothing; 1==dump event
		m5->SetDumpEventFileName("cosmic_dumpEvent.dat");
		m5->SetLargeMetEventFileName("cosmic_largeMetEvent.dat");
		m5->SetDatFileName("cosmic_jet.dat");
		
	
	TagPhoenixPhotons* phoenix = new TagPhoenixPhotons;

	TagConvElectrons* tagConvEle = new TagConvElectrons;
	

	TEmTimingModule* EmT = new TEmTimingModule();
	SetEmTimes* setEmTime = new SetEmTimes;  // set EM time of super photons
		setEmTime->SetEmTimeCorrectionFile("EmTCEMtable_p13.txt");


	EventProperties* eventPropMod = new EventProperties;

	CosmicStudy* cosmicStudMod = new CosmicStudy;

	//ap->AddModule(passFailFile,1);
	ap->AddModule(ptChecker,1);
	ap->AddModule(trigMod,1);
	ap->AddModule(myinit,1);
	ap->AddModule(tagTight,1);
	ap->AddModule(tagLoose,1);
	ap->AddModule(tagTightEle,1);
	ap->AddModule(tagLooseEle,1);
	ap->AddModule(peFilter,1);
	ap->AddModule(tagPMT,1);
	ap->AddModule(tagBH,1);
	ap->AddModule(phoenix,1);
	ap->AddModule(tagConvEle,1);
	ap->AddModule(EmT,1);
	ap->AddModule(setEmTime,1);
	ap->AddModule(m5,1);
//	ap->AddModule(jf,1);
	ap->AddModule(eventPropMod,1);
	ap->AddModule(cosmicStudMod,1);

	ap->GetInputModule()->SetPrintLevel(1);	// print file name as they are opened
	ap->Run(runEvts);
	//ap->Run(100000);
	//ap->ProcessEvent(207078,1150507);
	//ap->ProcessEvent(206990,247577);
  	
	
	std::string filename = "TestMetV1_CosmicData.root";
	ap->SaveHist(filename.c_str(),2);
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	//std::cout << "::::::::: Created Pass/Fail file => "<< passFailFile->GetFileName() << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
