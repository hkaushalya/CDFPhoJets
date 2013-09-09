#include <TChain.h>
#include <PhoFlat/PhotonJets.hh>
#include <PhoFlat/HaloStudy.hh>
#include <PhoFlat/HaloStudy_DelR.hh>
#include <PhoFlat/SidebandHalos.hh>
#include <PhoFlat/EleJetsTemp.hh>
#include <PhoFlat/LooseEleJetsTemp.hh>
#include <PhoFlat/HaloJetsTemp.hh>
#include <PhoFlat/CosmicJetsTemp.hh>
#include <PhoFlat/HaloIdEff.hh>
#include <PhoFlat/QCDJetsTemp.hh>
#include <PhoFlat/MCPhotonJets.hh>
#include <PhoFlat/UpgradeStuple.hh>
#include <PhoFlat/ewkMCPhotonJets.hh>
#include <PhoFlat/PhoIsoStudy.hh>
#include <PhoFlat/LoosePhoCuts.hh>
#include <PhoFlat/CosmicNoVtx.hh>
#include <PhoFlat/PhotonTriggerStudy.hh>
#include <PhoFlat/PhoenixPhoStudy.hh>
#include "TObjArray.h"
#include "TChainElement.h"
#include "TCollection.h"// for TIter
#include <sstream>
//#include "../RootTools/CommonTools.hh"
#include "../RootTools/IOColors.hh"
#include "TSystem.h"
#include <cmath>
#include <PhoFlat/QCDSystematics.hh>
#include <PhoFlat/ZSel.hh>
#include <PhoFlat/WSel.hh>
#include <PhoFlat/ZJetSel.hh>
#include <PhoFlat/WJetSel.hh>

//void DoRunMyJob (const int events, const int dataset=0, const int tightencut=0)
void DoRunMyJob2 (const int events, const int dataset=0, const int stupleVer=7)
{

	TChain *ch = new TChain("Stuple");
 	const std::string fileprefix="PhoJets";
  	std::stringstream rootfile;
	
	std::string outfile(gSystem->Getenv("OUT_FILE"));
	if (outfile.length()>0)
	{
		std::cout << "OUT_FILE = " << outfile << std::endl;
		rootfile << outfile;
	} else 
	{
		std::cout << "OUT_FILE env var not found! using default file name '" << fileprefix << "'" << std::endl;
  		rootfile << fileprefix;
	}

	

	int MCFLAG = 0;
	if (dataset != 1) MCFLAG = 1;

	if (stupleVer == 7)
	{
		std::cout << blue << std::endl;

		if (dataset == 1)
		{
			std::cout << "Using V7 PHO DATA ("<< stupleVer << ")"  << std::endl;
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_1P4.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_1P4_1.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_5P10.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_5P10_1.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_5P10_2.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_5P10_3.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_11P13.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_11P13_1.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_11P13_2.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_11P13_3.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_14P17_1.root");
			ch->Add("~/STUPLES2/StupleV7_DATA_Pho_NoJetVtxCuts_14P17_2.root");
			ch->Add("~/STUPLES3/StupleV7DATA_Pho_NoJetVtxCuts_18P26.root");
			ch->Add("~/STUPLES3/StupleV7DATA_Pho_NoJetVtxCuts_18P26_1.root");
			ch->Add("~/STUPLES3/StupleV7DATA_Pho_NoJetVtxCuts_18P26_2.root");
			ch->Add("~/STUPLES3/StupleV7DATA_Pho_NoJetVtxCuts_18P26_3.root");
			ch->Add("~/STUPLES3/StupleV7DATA_Pho_NoJetVtxCuts_18P26_4.root");
			ch->Add("~/STUPLES3/StupleV7DATA_Pho_NoJetVtxCuts_18P26_5.root");
			ch->Add("~/STUPLES3/StupleV7DATA_Pho_NoJetVtxCuts_18P26_6.root");
			ch->Add("~/STUPLES3/StupleV7DATA_Pho_NoJetVtxCuts_18P26_7.root");
			ch->Add("~/STUPLES3/StupleV7DATA_Pho_NoJetVtxCuts_18P26_8.root");
				
			
		} else if (dataset == 2)
		{
			std::cout << "Using V7 PHO MC with MINBIAS" << std::endl;
			ch->Add("~/STUPLES/StupleV7_MCPythia_PhoJetMinBiasGen6_gq0sqd_NoJetVtxCuts.root");
		} else if (dataset == 3)
		{
			std::cout << "Using V7 ZEE MC" << std::endl;
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zee_ze0scd_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zee_ze0sdd_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zee_ze0sed_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zee_ze0see_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zee_ze0seh_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zee_ze0sej_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zee_ze1s6d_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zee_ze1sad_NoJetVtxCuts_1.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zee_ze1sad_NoJetVtxCuts.root");
			

		} else if (dataset == 4)
		{
			std::cout << "Using V7 ZMM MC" << std::endl;
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zmm_ze0sbm_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zmm_ze0scm_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zmm_ze0sdm_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zmm_ze0sem_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zmm_ze0sfm_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zmm_ze0sgm_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zmm_ze1s6m_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Zmm_ze1s9m_NoJetVtxCuts.root");

		} else if (dataset == 5)
		{
			std::cout << "Using V7 ZTT MC" << std::endl;
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Ztt_ze0s8t_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Ztt_ze0sat_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Ztt_ze0sbt_NoJetVtxCuts.root");

		} else if (dataset == 6)
		{
			std::cout << "Using V7 WEN MC" << std::endl;
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wen_we0seh_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wen_we0sej_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wen_we0sfe_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wen_we0sge_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wen_we0sge_NoJetVtxCuts_1.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wen_we0she_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wen_we0sie_NoJetVtxCuts.root");

		} else if (dataset == 7)
		{
			std::cout << "Using V7 WMN MC" << std::endl;
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wmn_we0s7m_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wmn_we0s8m_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wmn_we0s9m_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wmn_we0sam_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wmn_we0sbm_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wmn_we0sgm_NoJetVtxCuts.root");

		} else if (dataset == 8)
		{
			std::cout << "Using V7 WTN MC" << std::endl;
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wtn_we0s9t_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wtn_we0sat_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MCPythia_Wtn_we0sbt_NoJetVtxCuts.root");

		} else if (dataset == 9)
		{

			std::cout << "Using V7 DIPHO MC" << std::endl;
			ch->Add("~/STUPLES2/StupleV7_MC_DiPho_NoJetVtxCuts.root");
			ch->Add("~/STUPLES2/StupleV7_MC_DiPho_NoJetVtxCuts_1.root");
		} else if (dataset == 10)
		{
			std::cout << "Using V8 WW MC" << std::endl;
			ch->Add("~/STUPLES4/StupleV8_MCPythia_WW_NoJetVtxCuts.root");
		} else if (dataset == 11)
		{
			std::cout << "Using V8 WZ MC" << std::endl;
			ch->Add("~/STUPLES4/StupleV8_MCPythia_WZ_NoJetVtxCuts.root");
			ch->Add("~/STUPLES4/StupleV8_MCPythia_WZ_NoJetVtxCuts_1.root");
		} else if (dataset == 12)
		{
			std::cout << "Using V8 ZZ MC" << std::endl;
			ch->Add("~/STUPLES4/StupleV8_MCPythia_ZZ_NoJetVtxCuts.root");
			ch->Add("~/STUPLES4/StupleV8_MCPythia_ZZ_NoJetVtxCuts_1.root");
		} else if (dataset == 13)
		{
			std::cout << "Using V8 TTBAR MC" << std::endl;
			ch->Add("~/STUPLES4/StupleV8_MCPythia_TTBAR_NoJetVtxCuts.root");


		}
	}


	TObjArray *fileElements=ch->GetListOfFiles();
	TIter next(fileElements);
	TChainElement *chEl=0;
	while (( chEl=(TChainElement*)next() )) 
	{
		TFile f(chEl->GetTitle());
		assert ( (! f.IsZombie()) && " a file is not found. pl check!");
		std::cout << "Attached data file: ";  f.Print();
	}
	std::cout << __FILE__ << ": Total Entries in chain = " << ch->GetEntries() << clearatt << std::endl;

	//std::cout << __FILE__ << " ch = " << ch << std::endl;
	TObjArray *brs = ch->GetListOfBranches();
	TIter n(brs);
	TChainElement *b=0;
	while (( b=(TChainElement*)n()))
	{
		//b->Print();
	}
	


	
//	ch->Draw("vtx_NClass12");

	
 	int iAddIsoToSidebandPhoton = atoi(gSystem->Getenv("ADD_ISO_2_SIDEBAND"));
	if (std::isnan(iAddIsoToSidebandPhoton))
	{
		std::cout << "ADD_ISO_2_SIDEBAND env var not found! returning." << std::endl;
	} else 
	{
		std::cout << "ADD_ISO_2_SIDEBAND = " << iAddIsoToSidebandPhoton << std::endl;
	}
	int iRewgtSidebandPhoton = atoi(gSystem->Getenv("REWGT_SIDEBAND"));
	if (std::isnan(iRewgtSidebandPhoton))
	{
		std::cout << "REWGT_SIDEBAND env var not found! returning." << std::endl;
	} else 
	{
		std::cout << "REWGT_SIDEBAND = " << iRewgtSidebandPhoton << std::endl;
	}

	const int iMinVtx = 1;
	const int iMaxVtx = 100;
	const float fMinPhoEt = 30.;
	const float fMaxPhoEt = 1200.;
	const float fMinJetEt = 15.;
	const float fMaxJetEt = 1200.;
	const float fMaxPhoEta = 1.1;
	const float fMinJetEta = -3.0;
	const float fMaxJetEta = 3.0;
	const int iMinNjet15 = 1;
	const int iMaxNjet15 = 100;
	const float fMinMet = 20;
	const float fMaxMet = 1500;
	const bool bDoMetCleanUp = 1;
	int iUseNvtxWgts = 1;
	if (dataset != 1) iUseNvtxWgts = 1;
	std::cout << "iUseNvtxWgts = " << iUseNvtxWgts << std::endl;

	
/*	PhotonJets *myPho = new PhotonJets;
	myPho->SetHistFileName(rootfile.str());
	myPho->SetMinPhotonEt(fMinPhoEt);
	myPho->SetMaxPhotonEt(fMaxPhoEt);
	myPho->SetMaxPhotonEta(fMaxPhoEta);
	myPho->SetMinJetEt(fMinJetEt);
	myPho->SetMaxJetEt(fMaxJetEt);
	myPho->SetMinJetEta(fMinJetEta);
	myPho->SetMaxJetEta(fMaxJetEta);
	myPho->SetSignalTimes(-4.8,4.8);
	myPho->SetCosmicTimes(30,90);
	myPho->SetMinClass12Vtx(iMinVtx);
	myPho->SetMaxClass12Vtx(iMaxVtx);
	myPho->SetMinNjet15(iMinNjet15);
	myPho->SetMaxNjet15(iMaxNjet15);
	myPho->SetHaloType(5);
	myPho->SetMinMet(fMinMet);
	myPho->SetMaxMet(fMaxMet);
	myPho->SetApplyMetCleanUpCuts(bDoMetCleanUp);
	myPho->SetMcFlag(MCFLAG);		//0=data ,1=MC
	myPho->SetAddIsoEtoSidebandPho(iAddIsoToSidebandPhoton); 
	myPho->ReweightSideband(iRewgtSidebandPhoton); 
	myPho->UseDATAMCNvtxWeights(iUseNvtxWgts); 
	//myPho->UseDATAMCNvtxWeights(0); 
	myPho->SetReportProgress(50000);
	myPho->SetDataset(dataset);
	myPho->SetPrintLevel(0);
	myPho->Main(ch,events);
*/	


/*	ZSel* zs = new ZSel();
	 zs->SetHistFileName(rootfile.str());
	 zs->SetMcFlag(MCFLAG);		//0=data ,1=MC
	 zs->SetUseNvtxWeights(1);
	 zs->Main(ch,events);
*/
/*	WSel* ws = new WSel();
	 ws->SetHistFileName(rootfile.str());
	 ws->SetMcFlag(MCFLAG);		//0=data ,1=MC
	 ws->SetUseNvtxWeights(0);
	 ws->Main(ch,events);
*/

	
	/* WJetSel* wjs = new WJetSel();
	 wjs->SetHistFileName(rootfile.str());
	 wjs->SetMcFlag(MCFLAG);		//0=data ,1=MC
	 wjs->SetUseNvtxWeights(0);
	 wjs->SetMinEleEt(20);
	 wjs->SetMinMet(20);
	 wjs->SetMinWPt(30);
	 wjs->SetMinJetEt(15.);
	 wjs->SetDataset(dataset);
	 wjs->Main(ch,events);
	 */

	
/*	ZJetSel* zjs = new ZJetSel();
	 zjs->SetHistFileName(rootfile.str());
	 zjs->SetMcFlag(MCFLAG);		//0=data ,1=MC
	 zjs->SetUseNvtxWeights(0);
	 zjs->SetMinZPt(30);
	 zjs->SetMaxZEta(1.1);
	 zjs->SetMinJetEt(15.);
	 zjs->SetMaxJetEta(3.);
	 zjs->Main(ch,events);
*/
/*	CosmicJetsTemp *myCosmic = new CosmicJetsTemp;
	myCosmic->SetHistFileName("CosmicJets.root");
	myCosmic->Main(ch,events);
*/

/*
	HaloJetsTemp *myHalo = new HaloJetsTemp;
	myHalo->SetHistFileName("HaloJets.root");
	myHalo->Main(ch,events);
*/

/*
	QCDJetsTemp *myQCD = new QCDJetsTemp;
	myQCD->SetHistFileName("QCDJets.root");
	myQCD->Main(ch,events);

*/
	//SidebandHalos *mySideHalo = new SidebandHalos;
	//mySideHalo->SetHistFileName("SidebandHalo.root");
	//mySideHalo->Main(ch,events);

/*	HaloIdEff *myHaloEff = new HaloIdEff;
	myHaloEff->SetHistFileName("HaloIdEff_g30_1njet15.root");
	myHaloEff->SetReportProgress(100000);
	myHaloEff->SetMinPhoEt(30);
	myHaloEff->SetMinPhoEta(1.1);
	myHaloEff->SetNVtx(1,1);
	myHaloEff->SetMinMet(0);
	myHaloEff->SetMinNjet15(1);
	myHaloEff->Main(ch,events);
*/

/*
	UpgradeStuple *stUp = new UpgradeStuple;
	stUp->SetHistFileName("/data/nbay02/a/samantha/STUPLESStupleV2_MC_PHO_noMinVtx_noMinJet_TLphoEleRemoved.root");
	stUp->Main(ch,events);
*/



/*	CosmicNoVtx *novtxCosmic = new CosmicNoVtx;
	novtxCosmic->SetHistFileName("NoVtxCosmics.root");
	novtxCosmic->Main(ch,events);
*/
/*
	if (dataset == 1)
	{
		for (int tightencut = 1; tightencut <=4; ++tightencut)
		{
			QCDSystematics *sys = new QCDSystematics;
			sys->SetMinPhoEt(fMinPhoEt);
			sys->SetMaxPhoEt(fMaxPhoEt);
			sys->SetMinJetEt(fMinJetEt);
			sys->SetMaxPhoEta(fMaxPhoEta);
			sys->SetMaxJetEta(fMaxJetEta);
			sys->SetMinMet(fMinMet);
			//const int tightencut = 4;
			std::stringstream sysfilename;
			if (tightencut == 1) sysfilename << "SidebandIDSystematics_TightHadEm_" << outfile;
			else if (tightencut == 2) sysfilename << "SidebandIDSystematics_TightIso_" << outfile;
			else if (tightencut == 3) sysfilename << "SidebandIDSystematics_TightTrkPt_" << outfile;
			else if (tightencut == 4) sysfilename << "SidebandIDSystematics_TightTrkIso_" << outfile;
			sys->SetTightenUpCut(tightencut);
			sys->SetHistFileName(sysfilename.str());
			sys->Loop(ch,events);
		}
	}
*/

/*
	PhoIsoStudy *phoIso = new PhoIsoStudy;
	phoIso->SetHistFileName("PhoIsoStudy.root");
	phoIso->Main(ch,events);
*/

/*	EleJetsTemp *ej = new EleJetsTemp();
	ej->SetHistFileName(rootfile.str());
	ej->SetMinEleEt(fMinPhoEt);
	ej->SetMaxEleEta(fMaxPhoEta);
	//ej->SetMinJetEt(15);
	//ej->SetMaxJetEta(3.2);
	ej->SetSignalTimes(-4.8,4.8);
	ej->SetMinClass12Vtx(1);
	ej->SetMaxClass12Vtx(100);
	ej->SetMinNjet15(1);
	ej->SetMinMet(0);
	ej->SetMcFlag(MCFLAG);		//0=data ,1=MC
	ej->SetReportProgress(10000);
	ej->Main(ch, events);
*/

	HaloStudy *haloS = new HaloStudy();
	haloS->SetHistFileName("HaloStudy.root");
	haloS->Main(ch,events);

	delete ch;
}
