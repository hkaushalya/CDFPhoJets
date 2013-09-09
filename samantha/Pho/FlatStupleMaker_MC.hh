#ifndef FLATSTUPLEMAKER_MC_HH
#define FLATSTUPLEMAKER_MC_HH

///////////////////////////////////////////////////////////
// See source file for class description                //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TTree.h"
#include <string>
#include <iostream>
#include <vector>
#include "samantha/obj/Stuple.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/Pho/TriggerModule.hh"
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "Stntuple/obj/TStnElectronBlock.hh"
#include "Stntuple/obj/TStnJetBlock.hh"
#include "samantha/Pho/TagBeamHalo.hh"
#include "Stntuple/obj/TGenpBlock.hh"

class FlatStupleMaker_MC: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()
		TStnElectronBlock*  fElectronBlock;
		TStnJetBlock*       fJetBlock;
		TGenpBlock* 		  fGenpBlock;
		TStnHeaderBlock*    fHeaderBlock;

	public:
		FlatStupleMaker_MC(const char* name="FlatStupleMaker_MC", const char* title="FlatStupleMaker_MC");
		~FlatStupleMaker_MC();

		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnJetBlock*      GetJetBlock()  	 { return fJetBlock; }
		TGenpBlock* 		 GetGenpBlock() 	 { return fGenpBlock;    }
		TStnHeaderBlock* 		 GetHeaderBlock()  { return fHeaderBlock;    }


				/************** abstract objects *************************/
		struct Global_Counters_t {
			unsigned int evtsProcessed;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int evtsWritten;
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
		void SkimEvent(std::vector<int> vPho, std::vector<int> vEle);
		void GetPhotonInfo(int iArrInd, SuperPhoton sp);
		void GetElectronInfo(int iArrInd, SuperPhoton sp);
		void GetTriggerInfo();
		void GetVertexInfo();
		void GetJetInfo();
		void GetMetInfo();
		void FillDataBlocks(int ientry);
		bool MatchToGen(std::vector<int> vPho);
		void SetStupleName(std::string name) { sFileName = name; }
		std::string GetStupleName(std::string name) { return sFileName; }

	private:
		Global_Counters_t counter;
		TFile *rootFile;
		TTree *tree;
		std::string sFileName;
		Stuple stuple;
		bool bRunPermit;

		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TriggerModule*  trigMod;
		TagBeamHalo* tagHaloMod;
		TH1F* hDelR;

	ClassDef(FlatStupleMaker_MC,1)
};

#endif
