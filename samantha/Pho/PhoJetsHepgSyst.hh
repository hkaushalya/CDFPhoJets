#ifndef PHOJETSHEPGSYST_HH
#define PHOJETSHEPGSYST_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/obj/Histograms.hh"
#include "samantha/Pho/EventProperties.hh"
#include <fstream>
#include "Stntuple/obj/TGenpBlock.hh"
#include <Stntuple/obj/TStnJetBlock.hh>
#include "samantha/utils/FreeFunctions.hh"
#include "TTree.h"

class PhoJetsHepgSyst: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TGenpBlock* fGenpBlock;
		TStnJetBlock*      fJetBlockClu04;
		TStnJetBlock*      fJetBlockHadClu04;

		static const int iSignalPhoton = 0;
		static const int iSidebandPhoton = 1;


	public:
		PhoJetsHepgSyst(const char* name="PhoJetsHepgSyst", const char* title="PhoJetsHepgSyst");
		~PhoJetsHepgSyst();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TGenpBlock* GetGenpBlock() { return fGenpBlock; }


				/************** abstract objects *************************/

		// enumerators to count the photon events
			enum COUNTS {
				EVENTS_PROCESSED = 0,
				EVENTS_PASSED    = 1,
				TIGHT_PHO_EVENTS   = 2,				//inclusive all events with a tight photon
				ONE_TIGHT_PHO_ONLY = 3,
				TWO_TIGHT_PHO_ONLY = 4,
				ONE_TIGHT_ANY_LOOSE_PHO = 5,
			};

		//holds final hist filled numbers to be written to the tree
		struct FinalVal_t {
			int RunNumber;
			int EvtNumber;
			float pj1_Et_pho;
			float pj1_Et_j1;
			float pj1_InvM_pj1;
			float pj1_Ht;
			int pj1_Njet15;
			float EtRatLeadJetLeadPho;
			float DelRLeadPhoLeadJet;
			float DelPhiLeadPhoLeadJet;
			float pj2_Et_pho;
			float pj2_Et_j1;
			float pj2_Et_j2;
			float pj2_InvM_pj1;
			float pj2_InvM_pj2;
			float pj2_InvM_pj1j2;
			float pj2_InvM_j1j2;
			float pj2_Ht;
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
		void FillDataBlocks(int ientry);

		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
		float GetHEPGHt(const TLorentzVector tlPho, const std::vector<TLorentzVector>& vJet);

	private:
		bool 	qMc;  				// true if MonteCarlo data

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		float fPhoMaxEta;   //absolute upper limit on Det Eta of the photon
		float fPhoMinEta;   //absolute lower limit on Det Eta of the photon
		float fPhoMinEt;			//Min Et of the photon
		float fJetMaxEta;   //absolute upper limit on Det Eta of the lead jet
		float fJetMinEt;			//Min Et of the lead jet
		
		bool bNoSummary;			//control the end job summary print
		TTree *tree;			//holds final(histogram) values to be used when normalizing samples
									// with different pt-hat
		FinalVal_t aFinalVal;
		int iJetCone;			//jet cone to use. chooses jet block based on this. default 0.4 cone

	public:
		
		float GetPhoMaxEta() const { return fPhoMaxEta; }
		float GetPhoMinEta() const { return fPhoMinEta; }
		float GetPhoMinEt() const { return fPhoMinEt; }
		void SetPhoMaxEta(const float eta);
		void SetPhoMinEta(const float eta);
		void SetPhoMinEt(const float et);

		float GetMaxJetEta() const { return fJetMaxEta; }
		float GetMinJetEt() const { return fJetMinEt; }
		void SetMaxJetEta(const float eta);
		void SetMinJetEt(const float et);
		void SetJetCone(const int i)
		{
			if (i>2 || i<0) 
			{
				std::cout << "Invalid JET cone specified, must be 0 (0.4), "
					<< "1(0.7), 2(1.0)! Using default 0.4 cone jets." <<std::endl;
			} else iJetCone = i;
		}
		float GetJetCone()
		{
			float cone = 0.4;
			if (iJetCone == 0) cone = 0.4;
			else if (iJetCone == 1) cone = 0.7;
			else if (iJetCone == 2) cone = 1.0;
			return cone;
		}
		
	ClassDef(PhoJetsHepgSyst,1)
};

#endif
