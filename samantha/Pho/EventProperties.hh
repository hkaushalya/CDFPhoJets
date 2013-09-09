#ifndef EVENTPROPERTIES_HH
#define EVENTPROPERTIES_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include "TLorentzVector.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/JetFilterModuleV3.hh"
#include "samantha/obj/Histograms.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/Pho/TriggerModule.hh"
#include "Stntuple/obj/TStnJetBlock.hh"

class EventProperties: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnJetBlock*       fJetBlock;


	public:
		EventProperties(const char* name="EventProperties", const char* title="EventProperties");
		~EventProperties();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnJetBlock*      GetJetBlock() 	{  return fJetBlock; }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsProcessed;	//___number of evts processed
			unsigned int evtsPassed;			//___number of evts pass this module
		};


		struct EvtProps_t {
			float MetCorr;
			float SumetCorr;
			float HtCorr;
			int NVertices;
			int N12Vertices;
			int NJet15;
			int RunNumber;
			int EvtNumber;
			float BestVertexZ;
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

		void FillEventHists(Histograms::EventHists_t& hist, const float fWeight=1);
		void FillPhotonHists(Histograms::PhotonHists_t& hist, int iSuperPhoIndex, float fWght=1);
		void FillJetHists(Histograms::JetHists_t& hist, int iIndex, float fWght=1);
		void FillPhoton1JetHists(Histograms::Photon1JetHists_t& hist,
								int iSuperPhoIndex, int iJetIndex, const float fWeight=1);
		void FillPhoton2JetsHists(Histograms::Photon2JetsHists_t& hist,
								int iSuperPhoIndex, int iJet1Index, int iJet2Index, const float fWeight=1);
		void FillTwoJetsHists(Histograms::TwoJetsHists_t& hist,
								int iJet1Index, int iJet2Index, const float fWeight=1);
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
	
	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		EvtProps_t aEvtProp;
		InitSuperPhotons* initSpMod;
		JetFilterModuleV3* jetMod;
		TriggerModule* trigMod;
		
		bool bRunPermit;				//checks all req. mods are present to run successfully
		
		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;
		bool bNoSummary;			//control the end job summary print

	ClassDef(EventProperties,1)
};

#endif
