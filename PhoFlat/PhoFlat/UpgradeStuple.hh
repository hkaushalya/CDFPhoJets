#ifndef UpgradeStuple_HH
#define UpgradeStuple_HH

#include "Stuple.hh"
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include <vector>
#include "TFolder.h"
#include "TDirectory.h"
#include "ReadInAll_4StUpgrade.hh"
#include "HistManager.hh"

class UpgradeStuple
{
	public:
		UpgradeStuple();
		int iProgressBy;
	
		void Main(TChain* chain, int iRunEvents = 0);
		void Init(TChain* chain);
		void CleanUp();
		std::string GetHistFileName() const { return sHistFileName; }
		void SetHistFileName(std::string name)  { sHistFileName = name; }
		
	private:
		int iEventsSelected;				// final number of events
		Stuple oldStuple, newStuple;
		TChain* myChain;
		std::string sHistFileName;

		TFile *rootFile;
		TTree *newTree;
		TDirectory* topDir;
		ReadInAll_4StUpgrade *readIn;
		HistManager myHistMan;
};
#endif
