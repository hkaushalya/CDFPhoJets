#ifndef DEBUGPHOET_HH
#define DEBUGPHOET_HH

#include "TTree.h"
#include "Stuple.hh"
#include "TChain.h"

class DebugPhoEt
{
	public:
	
		DebugPhoEt(std::string stupleFile);

	private:
		TChain* myTree;		//pointer to the tree/chain
		Stuple* myStuple;	//pointer to the stuple to be filled
		unsigned int iENTRIES;	//total number of entries
		unsigned int iEntry;  //current entry
};
#endif
