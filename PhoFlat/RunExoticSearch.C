#include <iostream>
#include <TChain.h>
#include "TObjArray.h"
#include "TChainElement.h"
#include "TCollection.h"// for TIter
#include <sstream>
//#include "../RootTools/CommonTools.hh"
#include "../RootTools/IOColors.hh"
#include "TSystem.h"
#include <cmath>

void RunExoticSearch(const int events)
{

	TChain *ch = new TChain("Stuple");
 	const std::string fileprefix="PhoJets";
  	std::stringstream rootfile;
	
	int MCFLAG = 0;

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
				
			

	TObjArray *fileElements=ch->GetListOfFiles();
	TIter next(fileElements);
	TChainElement *chEl=0;
	while (( chEl=(TChainElement*)next() )) 
	{
		TFile f(chEl->GetTitle());
		assert ( (! f.IsZombie()) && " a file is not found. pl check!");
		std::cout << "Attached data file: ";  f.Print();
	}
	std::cout << __FILE__ << ": Total Entries in chain = " << ch->GetEntries() << std::endl;

	std::cout << __FILE__ << " ch = " << ch << std::endl;
	TObjArray *brs = ch->GetListOfBranches();
	TIter n(brs);
	TChainElement *b=0;
	while (( b=(TChainElement*)n()))
	{
		//b->Print();
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
	const float fMinMet = 0;
	const float fMaxMet = 1500;
	const bool bDoMetCleanUp = 1;
	const int iUseNvtxWgts = 0;

	std::vector<std::string> vLibs;
	vLibs.push_back("Stuple.cc");
	vLibs.push_back("../RootTools/CommonTools.cc");
	vLibs.push_back("PhotonList.cc");
	vLibs.push_back("ElectronList.cc");
	vLibs.push_back("JetList.cc");
	vLibs.push_back("CommonVars.cc");
	vLibs.push_back("ReadInAll.cc");
	vLibs.push_back("ExoticEvts.C");

	bool bReady = true;
	for (unsigned i=0; i<vLibs.size(); ++i)
	{
		if (gSystem->CompileMacro(vLibs[i].c_str(),"k") == 1) continue;
		bReady = false;
		break;
	}
	
	ExoticEvts(ch, events);
	
	delete ch;
}
