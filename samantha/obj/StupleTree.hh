#ifndef STUPLETREE_HH
#define STUPLETREE_HH

#include "Stuple.hh"
#include "TTree.h"

class StupleTree
{
	public:
		StupleTree(const std::string sTreeName, TTree*, Stuple&);
	private:
};
#endif
