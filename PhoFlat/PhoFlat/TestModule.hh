#ifndef TestModule_HH
#define TestModule_HH

#include "Stuple.hh"
#include <iostream>
#include "TTree.h"
#include "TFile.h"

class TestModule
{
	public:
		TestModule();
		~TestModule();
		static const int iProgressBy = 1000;
	
		TTree* GetMyTree(TFile* file);
		void DoMyStuff(Stuple& stuple);
		void Init();
		void Main(int iRunEvents = 0);
		
	private:
		int iEventsSelected;		// final number of events
		Stuple stuple;
		TFile* myFile;
		TTree* tree;


	ClassDef(TestModule,1)
};
#endif
