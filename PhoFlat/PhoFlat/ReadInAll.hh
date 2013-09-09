#ifndef READINALL_HH
#define READINALL_HH

#include "TTree.h"
#include "Stuple.hh"
#include "TChain.h"

class ReadInAll
{
	public:
		ReadInAll(TChain*, Stuple*);
		void Init();
		void GetEntry(const int i);
		void AllPhotons(const int i) 		{ bAllPho = i; }
		void CentralPhotons(const int i)	{ bCentralPho = i; }
		void EmUpPhotons(const int i) 	{ bEmUpPho = i; }
		void EmDownPhotons(const int i) 	{ bEmDownPho = i; }
		void AllElectrons(const int i) 	{ bAllEle = i; }
		void CentralElectrons(const int i) { bCentralEle = i; }
		void EmUpElectrons(const int i) 	{ bEmUpEle = i; }
		void EmDownElectrons(const int i){ bEmDownEle = i; }

		void AllJets(const int i) 			{ bAllJet = i; }
		void CentralJets(const int i) 	{ bCentralJet = i; }
		void JesUpJets(const int i) 		{ bJesUpJet = i; }
		void JesDownJets(const int i) 	{ bJesDownJet = i; }

		void Met(const int i) 				{ bMet = i; }
		void GenLevel(const int i) 		{ bGen = i; }

	private:
		TChain* myTree;		//pointer to the tree/chain
		Stuple* myStuple;	//pointer to the stuple to be filled
		unsigned int iENTRIES;	//total number of entries
		unsigned int iEntry;  //current entry
		bool bAllPho, bCentralPho, bEmUpPho, bEmDownPho;
		bool bAllEle, bCentralEle, bEmUpEle, bEmDownEle;
		bool bAllJet, bCentralJet, bJesUpJet, bJesDownJet;
		bool bMet, bGen;
};
#endif
