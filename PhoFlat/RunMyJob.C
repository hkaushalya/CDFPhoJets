//void RunMyJob (int events = -1, int dataset=0, const int tightencut=1)
void RunMyJob (int events = -1, int dataset=0, int stupleVer=7)
{
	std::cout << "events, dataset = " << events << "," << dataset << std::endl;
	gROOT->Reset();
	gSystem->Load("$ROOTSYS/lib/libPhysics.so");
	//if (gSystem->Load("../RootTools/CommonTools_cc.so") == 0)
	{
	//	gSystem->Load("../RoooTools/CommonTools_cc.so");
	}

	// gROOT->LoadMacro() method does not seeme to return 1 when the compilation is a success.
	// not sure why . but gSystem->CompileMacro does
	//int st = gROOT->LoadMacro("Stuple.cc+");
	//
	std::vector<std::string> vLibs;
	//vLibs.push_back("Stuple.cc");
	//vLibs.push_back("../RootTools/CommonTools.cc");
	vLibs.push_back("PhotonList.cc");
	vLibs.push_back("ElectronList.cc");
	vLibs.push_back("JetList.cc");
	vLibs.push_back("CommonVars.cc");
	vLibs.push_back("ReadInAll.cc");
	vLibs.push_back("FreeFunctions.cc");
	vLibs.push_back("Histograms.cc");
	vLibs.push_back("HistManager.cc");
	vLibs.push_back("HaloTypes.cc");
	vLibs.push_back("NewJetList.cc");
	vLibs.push_back("HaloStudy.cc");
	vLibs.push_back("HaloStudy_DelR.cc");
	vLibs.push_back("PhotonJets.cc");
	  vLibs.push_back("HaloJetsTemp.cc");
	  vLibs.push_back("EleFakeRate.cc");
	  vLibs.push_back("EleJetsTemp.cc");
	  vLibs.push_back("LooseEleJetsTemp.cc");
	  vLibs.push_back("CosmicJetsTemp.cc");
	  vLibs.push_back("QCDJetsTemp.cc");
	  vLibs.push_back("SidebandHalos.cc");
	  vLibs.push_back("HaloIdEff.cc");
	  vLibs.push_back("MCPhotonJets.cc");
	  vLibs.push_back("ReadInAll_4StUpgrade.cc");
	  vLibs.push_back("UpgradeStuple.cc");
	  vLibs.push_back("ewkMCPhotonJets.cc");
	vLibs.push_back("LoosePhoCuts.cc");
	  vLibs.push_back("Systematics.cc");
	 vLibs.push_back("CosmicNoVtx.cc");
	vLibs.push_back("PhotonTriggerStudy.cc");
	vLibs.push_back("QCDSystematics.cc");
	//must load this last, otherwise there'll be dlopen errors!
	vLibs.push_back("ZSel.cc");
	vLibs.push_back("WSel.cc");
	vLibs.push_back("ZJetSel.cc");
	vLibs.push_back("WJetSel.cc");
	vLibs.push_back("DoRunMyJob.cc");


	bool bReady = true;
	for (unsigned i=0; i<vLibs.size(); ++i)
	{
		if (gSystem->CompileMacro(vLibs[i].c_str(),"k") == 1) continue;
		bReady = false;
		break;
	}

	if (! bReady) 
	{
		//gSystem->ListLibraries();
		std::cout << __FILE__ << ":: NOT ALL REQUIRED LIBRARIES ARE LOADED PROPERLY!" << std::endl;
		return;
	}
	// gObjectTable->Print();
	
	gSystem->Setenv("REWGT_SIDEBAND","0");
	gSystem->Setenv("ADD_ISO_2_SIDEBAND","0");

	//change file name for testing. DON'T BE LAZY! and mess up!
	//std::string prefix("/data/nbay04/c/samantha/PhoFlatResults/ICHEP_LaptopVersion/p1_26/PhoJets_");
	std::string prefix("/data/nbay04/c/samantha/PhoFlatResults/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_final/PhoJets_");
	//std::string prefix("WJets_");
	//std::string prefix("PhoJets_Cosmics");
	//std::string prefix("~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_final_PhoEt200/PhoJets_");
	//std::string prefix("~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_Met20_final/PhoJets_");
	//std::string prefix("~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_final/NvtxWgting_dibo_ttbar/PhoJets_");
	//std::string prefix("~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_Met20_final/WJetsStudy/Incl2Jet/WJets_2InclNjet15_");
	//std::string suffix("_Met50_NvtxSigmaWgtd_MetCleaned.root");
	//std::string suffix("_Met20NvtxWgtedCleaned.root");
	//std::string suffix("_Met30NvtxWgtdCleaned.root");
	//std::string suffix("_Met50.root");
	std::string suffix("_BuckleysWgZgTest.root");
	//std::string suffix("_Test.root");

	std::string sdataset("unknown");
	if (dataset == 1 || dataset == 50) sdataset = "data";
	//else if (dataset == 2) sdataset = "phomc";
	else if (dataset == 2) sdataset = "phomcWithMB";
	else if (dataset == 3) sdataset = "zeemc";
	else if (dataset == 4) sdataset = "zmmmc";
	else if (dataset == 5) sdataset = "zttmc";
	else if (dataset == 6) sdataset = "wenmc";
	else if (dataset == 7) sdataset = "wmnmc";
	else if (dataset == 8) sdataset = "wtnmc";
	else if (dataset == 9) sdataset = "diphomc";
	else if (dataset == 10) sdataset = "wwmc";
	else if (dataset == 11) sdataset = "wzmc";
	else if (dataset == 12) sdataset = "zzmc";
	else if (dataset == 13) sdataset = "ttbarmc";
	else if (dataset == 70) sdataset = "Wgmc";
	else if (dataset == 80) sdataset = "Zgmc";



	std::string sRewgt("");
	std::string rewgt(gSystem->Getenv("REWGT_SIDEBAND"));
	if (rewgt.length()>0)
	{
		if (atoi(rewgt.c_str()) == 1) sRewgt += "_SidebandNominalWgtByPhoEtAndJetEta_";
		else if (atoi(rewgt.c_str()) == 2) sRewgt += "_SidebandSigmaWgtByPhoEtAndJetEta_";
		else if (atoi(rewgt.c_str()) == 3) sRewgt += "_SidebandNominalWgtByPhoEt_";
		else if (atoi(rewgt.c_str()) == 4) sRewgt += "_SidebandSigmaWgtByPhoEt_";
	}
		
	
	std::ostringstream outfilename;
	outfilename << prefix << sdataset << sRewgt << suffix;
	gSystem->Setenv("OUT_FILE",outfilename.str().c_str());
	//gSystem->Setenv("OUT_FILE","PhoJets_wenmc_g30_withsideband_03282010.root");
	
	//DoRunMyJob (events, dataset, tightencut);
	//DoRunMyJob (events, dataset, stupleVer);
	//gObjectTable->Print();
}
