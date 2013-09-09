#include <iostream>
#include "PhoFlat/TestModule.hh"
#include "/PhoFlat/ReadInAll.hh"

ClassImp(TestModule)

TestModule::TestModule()
{
}
TestModule::~TestModule()
{
}

//---- get my Stuple ------------------------------------------------
TTree* TestModule::GetMyTree(TFile* file)
{
  file = new TFile("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/Stuples/p1_13/Stuples/MergedStuple.root", "READ");
  if (file->IsZombie()) {
	  	std::cout << __FILE__ << "::" << __LINE__ << "::"  << "file not found " << std::endl;
	    exit (1); 
	}
  TTree *tree = dynamic_cast<TTree*>(file->Get ("Stuple"));
  tree->SetBranchStatus ("*", 0);
  return tree;
}
//--- do all my stuff here ------------------------------------------
void TestModule::DoMyStuff(Stuple& stuple)
{
	for (unsigned int j = 0; j< stuple.pho_num; ++j ) {
		if (stuple.pho_PhoenixId[j] != 0) continue;
		++iEventsSelected;
	};

}

//------ main -------------------------------------------------------
void TestModule::Init()
{
	iEventsSelected = 0;
	tree->SetBranchStatus 	("pho_num", 1);
	tree->SetBranchStatus 	("ele_num", 1);
	tree->SetBranchStatus 	("evt_RunNumber", 1);
	tree->SetBranchAddress 	("pho_num", &stuple.pho_num);
	tree->SetBranchAddress 	("ele_num", &stuple.ele_num);
	tree->SetBranchAddress 	("evt_RunNumber", &stuple.evt_RunNumber);
}
//------ main -------------------------------------------------------
void TestModule::Main(int iRunEvents)
{
	tree = GetMyTree(myFile);
	Init();
	
	Double_t entries = tree->GetEntries();
	std::cout << "Entries found = " << entries << std::endl;
	if (iRunEvents > 0 && iRunEvents < entries) entries = iRunEvents;
	else std::cout << __FILE__ << "::" << __LINE__ << "::" << "Only " << entries << " found and those will be processed." << std::endl;
	
	for (Int_t i=0 ; i < entries; ++i) {
	  	tree->GetEntry(i);
		if (i !=0 && (i%iProgressBy == 0)) std::cout << i << "\t events processed." << std::endl;
		if (stuple.evt_RunNumber < 190851) continue;  //drop the first 400pb-1
	  	if (stuple.pho_num < 1) continue;
	  	if (stuple.ele_num > 0) continue; 					// reject e+pho events
		ReadInAll(tree, i, stuple);
		DoMyStuff(stuple);
	};
	
	std::cout << "Events Selected = " << iEventsSelected << std::endl;
}
