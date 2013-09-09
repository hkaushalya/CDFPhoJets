///////////////////////////////////////////////////////////
// See sourcce file for class description                 //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#ifndef PASSFAILTREE_HH
#define PASSFAILTREE_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TTree.h"
#include <string>
#include <iostream>
#endif

class PassFailTree: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		//TStnHeaderBlock*    fHeaderBlock;


	public:
		PassFailTree(const char* name="PassFailTree", const char* title="PassFailTree");
		~PassFailTree();

		// ****** accessors
		//TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/
		struct Global_Counters_t {
			unsigned int evtsProcessed;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
		};



		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
		void Cleanup();
		std::string GetFileName() const { return sFileName; }
		void SetFileName(std::string fn_) { sFileName = fn_; }
		void CheckAndFill(); 

	private:
		bool bMcFlag;
		Global_Counters_t counter;
		TFile *rootFile;
		TTree *tree;
		std::string sFileName;
		int iCurrRun, iCurrEvent;
		int iPrevRun, iPrevEvent;
		bool bSkipThis;
  		TList *modList;
		TIterator* it;
		std::vector<TStnModule*> vTstnMods;
		std::vector<int> vPassFail;


	ClassDef(PassFailTree,1)
};

#endif
