#ifndef NEWJETLIST_HH
#define NEWJETLIST_HH

#include "Stuple.hh"
#include <string>
#include <vector>
#include "FreeFunctions.hh"
#include "PhotonList.hh"
#include "ElectronList.hh"
#include "JetList.hh"


class NewJetList
{
	public:
		NewJetList();
		void AddUnused(Stuple& stuple, std::vector<int>UsedPho, std::vector<int>UsedEle);
		void AddUnusedCentral(Stuple& stuple, std::vector<int>UsedPho, std::vector<int>UsedEle);		//add unused central em objs back to central jet list
		void AddUnusedUp(Stuple& stuple, std::vector<int>UsedPho, std::vector<int>UsedEle);				//add unused em/jes up em objs back to jes up jet list
		void AddUnusedDown(Stuple& stuple, std::vector<int>UsedPho, std::vector<int>UsedEle);			//add unused em/jes down em objs back to jes down jet list
		
		NewJetList(Stuple&, int iObjType);		// ObjType == 0 (photon), 1 (electron)
		bool IsUsed(int, std::vector<int>);
		void AddEleToJetList(Stuple&, std::vector<int> EleNoMatch);
		void AddPhoToJetList(Stuple&, std::vector<int> PhoNoMatch);
		void ReorderJets(Stuple&);
	

		void AddEleToCentralJetList(Stuple&, std::vector<int> EleNoMatch);
		void AddPhoToCentralJetList(Stuple&, std::vector<int> PhoNoMatch);
		void ReorderCentralJets(Stuple&);
		void AddUpEleToUpJetList(Stuple&, std::vector<int> EleNoMatch);
		void AddUpPhoToUpJetList(Stuple&, std::vector<int> PhoNoMatch);
		void ReorderUpJets(Stuple&);
		void AddDownEleToDownJetList(Stuple&, std::vector<int> EleNoMatch);
		void AddDownPhoToDownJetList(Stuple&, std::vector<int> PhoNoMatch);
		void ReorderDownJets(Stuple&);



		void AddUnused(PhotonList& phos, ElectronList& eles, JetList& jets,
							const std::vector<int>UsedPho, const std::vector<int>UsedEle,
							const float fMinJetEt, const float fMaxJetEta);		//add unused central em objs back to central jet list
		void ReorderJets(JetList& jets);

	
};
#endif
