#include <iostream>
#include <set>
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <../PhoFlat/PhoFlat/ReadInAll.hh>
#include <../PhoFlat/PhoFlat/Stuple.hh>
#include <fstream>
#include <TSystem.h>
#include <vector>
#include <TChain.h>
#include <TH1F.h>

void DebugPhotonEt (std::string filepath=""
							, unsigned int start=0, unsigned int end=999999999)
{

	std::vector<std::string> files;  
	std::vector<std::string>::const_iterator fileIt;  
 
	if (filepath.length()) {
		files.push_back(filepath.c_str());
	} else {
	 	//files.push_back("~/stn_rel/Stuple_cph1i_allEvts.root");
	 	files.push_back("~/stn_rel/Stuple.root");
 		//files.push_back("~/RESULTS/JetEmMisMatch/Stuple.root");
	
// 	files.push_back("~/RESULTS/STUPLES/StupleV2_noMinVtx_noMinJet_TLphoEleRemoved.root");
// 	files.push_back("~/RESULTS/STUPLES/StupleV2_noMinVtx_noMinJet_TLphoEleRemoved_1.root");


	}
	
	TChain *tree = new TChain("Stuple");

  for (fileIt = files.begin(); fileIt != files.end(); ++fileIt) {
	  TFile* file = new TFile(fileIt->c_str(),"READ");
	  if (file->IsZombie()) {
			std::cout << "Cannot find file:" << *fileIt  << " pl. check." << std::endl;
			exit (1);
  		} else {
			tree->Add(fileIt->c_str());
		}
  }

	Stuple st;
	tree->SetBranchStatus  ("*", 0);
	tree->SetBranchStatus  ("evt_RunNumber", 1);
	tree->SetBranchAddress ("evt_RunNumber", &st.evt_RunNumber);
	tree->SetBranchStatus  ("evt_EventNumber", 1);
	tree->SetBranchAddress ("evt_EventNumber", &st.evt_EventNumber);
	tree->SetBranchStatus ("pho_num", 1);
	tree->SetBranchAddress ("pho_num", &st.pho_num);
	tree->SetBranchStatus ("pho_Etc", 1);
	tree->SetBranchAddress ("pho_Etc", &st.pho_Etc);
	tree->SetBranchStatus ("ele_num", 1);
	tree->SetBranchAddress ("ele_num", &st.ele_num);
	tree->SetBranchStatus ("ele_Etc", 1);
	tree->SetBranchAddress ("ele_Etc", &st.ele_Etc);
	tree->SetBranchStatus ("jet_num", 1);
	tree->SetBranchAddress ("jet_num", &st.jet_num);
	tree->SetBranchStatus ("jet_Pt", 1);
	tree->SetBranchAddress ("jet_Pt", &st.jet_Pt);
	

	unsigned event = start;
	unsigned ENTRIES = (unsigned)tree->GetEntries();
	std::cout << "Entries Found = " << ENTRIES << std::endl;

	while (tree->GetEntry (event)) {
						
//		std::cout << "======================================" <<std::endl;
//		std::cout << "run,evt=" << st.evt_RunNumber << ", " << st.evt_EventNumber << std::endl;
//		std::cout << "pho_num=" << st.pho_num << std::endl;
		std::cout << "\npho_et= ";
		for (unsigned i = 0; i < st.pho_num; ++i) {
			if (st.pho_Etc[i] < 0)
				std::cout << "\t\t" << i << "\tPho Et=" << st.pho_Etc[i] << std::endl;
			else std::cout << "\t" << st.pho_Etc[i] << "\t";
		}
		//std::cout << "ele_num=" << st.ele_num << std::endl;
		std::cout << "\nele_et= ";
		for (unsigned i = 0; i < st.ele_num; ++i) {
			if (st.ele_Etc[i] < 0)
				std::cout << "\t" << i << "\tEle Et=" << st.ele_Etc[i] << std::endl;
			else std::cout << "\t" << st.pho_Etc[i] << "\t";
		}
//		std::cout << "jet_num=" << st.jet_num << std::endl;
//		for (unsigned i = 0; i < st.jet_num; ++i) {
//			std::cout << "\t" << i << "\tEt=" << st.jet_Pt[i] << std::endl;
//		}
		
		++event;
	}


};
