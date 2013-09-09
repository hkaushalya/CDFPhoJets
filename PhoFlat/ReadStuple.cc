#include "TFile.h"
#include "TTree.h"
#include <iostream>

// TO run:
// root>.x ReadStuple.C+(10)

void ReadStuple(const unsigned int nevts=10)
{
	TFile f("PHOJETS_DATA_NTUPLE.root");		//change you input file name
	if (f.IsZombie())
	{
		std::cout << __LINE__ << ":: File not found" << std::endl;
		exit(1);
	}

	TTree *myTree = dynamic_cast<TTree*> (f.Get("Stuple"));
	assert (myTree != NULL && "tree not found!") ;
	
	std::cout << "entries = " << myTree->GetEntries() << std::endl;

	int iRunNumber, iEvtNumber;
	int iNPhoNum;		//each event with have many photons and each photon
							//info is saved in an array.
	float iPhoEtc[iNPhoNum];		//just like this.
	float iPhoTightId[iNPhoNum];		//just like this.
	float iPhoLooseId[iNPhoNum];		//just like this.
	int iNtight;
	
	//enable only the braches you want to read (faster)
	myTree->SetBranchStatus("*", 0);
	myTree->SetBranchStatus("evt_EventNumber", 1);
	myTree->SetBranchStatus("evt_RunNumber", 1);
	myTree->SetBranchStatus("pho_num", 1);
	myTree->SetBranchStatus("pho_Ntight", 1);
	myTree->SetBranchStatus("pho_Etc", 1);
	myTree->SetBranchStatus("pho_TightId", 1);
	myTree->SetBranchStatus("pho_LooseId", 1);

	//allocate memory to read in from tree
	myTree->SetBranchAddress("evt_RunNumber", &iRunNumber);
	myTree->SetBranchAddress("evt_EventNumber", &iEvtNumber);
	myTree->SetBranchAddress("pho_num", &iNPhoNum);
	myTree->SetBranchAddress("pho_Ntight", &iNtight);
	myTree->SetBranchAddress("pho_Etc", iPhoEtc);
	myTree->SetBranchAddress("pho_TightId", iPhoTightId);
	myTree->SetBranchAddress("pho_LooseId", iPhoLooseId);

	unsigned int entries = (unsigned int) myTree->GetEntries();

	if (nevts<entries) entries = nevts;

	for (unsigned int i = 0; i < entries; ++i)
	{
		//to fill in the allocated memory
		myTree->GetEntry(i);
		std::cout << "==========Run,Evt:" << iRunNumber << ", "<< iEvtNumber << std::endl;
		std::cout << "\tPhoton list: Ntight= " << iNtight << std::endl;
		std::cout << "\tIndex\tEtc\tTightId\tLooseId" <<std::endl;
		//now loop over the photons in the event
		for (int j=0; j < iNPhoNum; j++)
		{
			std::cout << "\t" << j << "\t" << iPhoEtc[j]  << "\t" 
						<< iPhoTightId[j] << "\t" << iPhoLooseId[j] << std::endl;
		}

	}

	f.Close();

}	
