////////////////////////////////////////////////////////////////
// This is to test the Stuple by reading it back.             //
////////////////////////////////////////////////////////////////

#include <iostream>
#include "PhoFlat/Stuple.hh"
#include <TChain.h>
//-------------------------------------------------------------------
void TestReadStuple()
{
	  TChain *myTree = new TChain("Stuple");

 //DATA
 //myTree->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/STUPLES/StupleV2_noMinVtx_noMinJet_TLphoEleRemoved.root");
 myTree->Add("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/STUPLES/Stuple_noMinVtx_noMinJet_TLphoEleRemoved.root");

	Stuple *myStuple = new Stuple; 

	myTree->SetBranchStatus("pho_num", 1);		// number of electrons so we know how much to read from arrays
	myTree->SetBranchStatus("ele_num", 1);
	myTree->SetBranchStatus("jet_num", 1);

	myTree->SetBranchStatus("ele_Index", 1);
	myTree->SetBranchStatus("ele_Etc", 1);
	myTree->SetBranchStatus("ele_Ntight", 1);
	myTree->SetBranchStatus("ele_Nloose", 1);
	myTree->SetBranchStatus("ele_matchJetIndex", 1);	//if a matching jet is found and removed from jet list


	myTree->SetBranchAddress("pho_num", &myStuple->pho_num);		// number of electrons so we know how much to read from arrays
	myTree->SetBranchAddress("ele_num", &myStuple->ele_num);
	myTree->SetBranchAddress("jet_num", &myStuple->jet_num);

	myTree->SetBranchAddress("ele_Ntight", &myStuple->ele_Ntight);
	myTree->SetBranchAddress("ele_Nloose", &myStuple->ele_Nloose);
	myTree->SetBranchAddress("ele_Index", &myStuple->ele_Index);
	myTree->SetBranchAddress("ele_Etc", &myStuple->ele_Etc);
	myTree->SetBranchAddress("ele_matchJetIndex", &myStuple->ele_matchJetIndex);	//if a matching jet is found and removed from jet list


	//for (unsigned i =0 ; i < myTree->GetEntries(); ++i) {
	for (unsigned i =0 ; i < 100; ++i) {
		myTree->GetEntry(i);
		std::cout << myStuple->ele_num << std::endl;
		for (unsigned j = 0 ; j < myStuple->ele_num; ++j) {
			std::cout << myStuple->ele_matchJetIndex[j] << std::endl;

		}
	}
}

