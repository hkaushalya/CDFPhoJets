#ifndef ReadInAll_4StUpgrade_HH
#define ReadInAll_4StUpgrade_HH

#include "TTree.h"
#include "Stuple.hh"
#include "TChain.h"

class ReadInAll_4StUpgrade
{
	public:
	
		ReadInAll_4StUpgrade(TChain*, Stuple*);
		void GetEntry(const int i);

	private:
		TChain* myTree;		//pointer to the tree/chain
		Stuple* myStuple;	//pointer to the stuple to be filled
		unsigned int iENTRIES;	//total number of entries
		unsigned int iEntry;  //current entry
};
#endif
