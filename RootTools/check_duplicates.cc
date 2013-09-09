#include <iostream>
#include <set>
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <sstream>
#include <TSystem.h>
#include <vector>
#include <TChain.h>
#include "../samantha/obj/Stuple.hh"
#include "CommonTools.hh"
#include "TSystem.h"

//void check_duplicates (std::string filepath, int start, int end)
void check_duplicates (std::string filepath=""
							, unsigned int start=0, unsigned int end=999999999)
{

	gSystem->Load("../PhoFlat/Stuple_cc.so");
	std::vector<std::string> files;  
	std::vector<std::string>::const_iterator fileIt;  
 
	if (filepath.length()) {
		files.push_back(filepath.c_str());
	} else {

		files.push_back("StupleV8_PhoData_5P10_NoJetVtxCuts_bTagInfo.root");
		files.push_back("StupleV8_PhoData_5P10_NoJetVtxCuts_bTagInfo_1.root");
		files.push_back("StupleV8_PhoData_5P10_NoJetVtxCuts_bTagInfo_2.root");
		files.push_back("StupleV8_PhoData_5P10_NoJetVtxCuts_bTagInfo_3.root");


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

 // TTree *tree = dynamic_cast<TTree*>(file.Get("Stuple"));
 // if (!tree) {
//		std::cout << "Cannot find a Stuple tree in the given file. pl. check." << std::endl;
//		exit (1);
//  }


  unsigned runnr, eventnr;

  tree->SetBranchStatus  ("*", 0);
  tree->SetBranchStatus  ("evt_RunNumber", 1);
  tree->SetBranchAddress ("evt_RunNumber", &runnr);
  tree->SetBranchStatus  ("evt_EventNumber", 1);
  tree->SetBranchAddress ("evt_EventNumber", &eventnr);

  unsigned event = start;
  std::set<std::pair<unsigned,unsigned> > found, dupEvts;
  std::set<std::pair<unsigned,unsigned> >::iterator myIt;
  unsigned ENTRIES = (unsigned)tree->GetEntries();
  std::cout << "Entries Found = " << ENTRIES << std::endl;

  std::vector<std::pair<unsigned,unsigned> > dups;

  while (tree->GetEvent (event))
  {
	  ++ event;
	  /*if ((event % 1000000) == 0)
	  {
		  std::cout << "Event: " << event << std::endl;
	  }
	  */
	  if ((event % 100000) == 0) std::cout << (int) (event * 100. /tree->GetEntries())<< " % complete." << std::endl; 
     
	  std::pair<unsigned,unsigned> id (runnr, eventnr);

	  myIt = found.find(id);
	  if (myIt != found.end())
	  {
		 std::cout << myIt->first << "\t" << myIt->second <<std::endl;
		 //dups.push_back(id);
	    dupEvts.insert (id);
	  }
     
	  found.insert (id);
	 
	  if (event >=end) break;
  	}
	


  
  //std::cout << found.size() << "/" << event << std::endl;
  std::cout << found.size() << "/" << event - start << std::endl;
  std::cout << "Duplicate events = " << event - start - found.size() << " found in root file " << filepath << std::endl;
  
 // for (std::set<std::pair<unsigned,unsigned> >::const_iterator it = found.begin() ; it != found.end(); ++it) {
 //		std::cout << it->first << ",\t" << it->second << std::endl;
 // }

  if (dupEvts.size()<1) return;

	//if duplicates are found, run this additional option
	//to create a new ntuple without any duplcates.
  
  TTree *oldTree=0;
  Stuple oldStuple;

  tree->SetBranchStatus  ("*", 0);
  PrepStupleTreeForReading(tree, &oldStuple);
  //oldStuple.Dump();

  //tree->SetBranchStatus  ("evt_RunNumber", 1);
  //tree->SetBranchAddress ("evt_RunNumber", &(oldStuple.evt_RunNumber));
  //tree->SetBranchStatus  ("evt_EventNumber", 1);
  //tree->SetBranchAddress ("evt_EventNumber", &(oldStuple.evt_EventNumber));


  
  std::cout << "ch = " << tree << std::endl;
  std::cout << "N = " << tree->GetEntries() << std::endl;
  

  TTree* newTree=0;
  Stuple newStuple;
  
  std::stringstream newfile;
  newfile << files.at(0) << "_NoDuplicates";
  TFile f(newfile.str().c_str(),"RECREATE");
  newTree = GetStupleTree("Stuple", newStuple);
  assert (newTree != NULL && "new tree is null");

  std::set<std::pair<unsigned,unsigned> >::iterator dupIt;
  
  unsigned iDupCount=0, iOrigCount=0;
  //for (unsigned i=1; i< 50; ++i)
  for (unsigned i=0; i< tree->GetEntries(); ++i)
  {
	  if ((i % 100000) == 0) std::cout << (int) (i * 100. /tree->GetEntries())<< " % complete." << std::endl; 
	  tree->GetEntry(i);
	  //std::cout << oldStuple.evt_RunNumber << ", "<< oldStuple.evt_EventNumber << std::endl;
	  std::pair<unsigned,unsigned> id (oldStuple.evt_RunNumber, oldStuple.evt_EventNumber);

	  dupIt = dupEvts.find(id);
	  if (dupIt == dupEvts.end())
	  {
		  newStuple.Init();
		  newStuple = oldStuple;
		  newTree->Fill();
		  iOrigCount++;
	  } else {
		  dupEvts.erase(dupIt);
		  std::cout << "dulicate found " << dupIt->first << "\t" << dupIt->second << " ------------> dupEvts.size = " << dupEvts.size() << std::endl;
		  iDupCount++;
	  }
  }

  std::cout << "orig / dulicates " << iOrigCount << " / " << iDupCount <<std::endl;
	
  f.Write();
  f.ls();
  f.Close();
  
  //delete tree;
  //delete oldTree;
  //delete newTree;
};
