#ifndef HALOSTUDYALLWITHMETET_HH
#define HALOSTUDYALLWITHMETET_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/Pho/TagBeamHalo.hh"
#include "samantha/obj/Histograms.hh"
#include "samantha/Pho/EventProperties.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/SetEmTimes.hh"


class HaloStudyAllWithMetEt: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		HaloStudyAllWithMetEt(const char* name="HaloStudyAllWithMetEt", const char* title="HaloStudyAllWithMetEt");
		~HaloStudyAllWithMetEt();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int haloEventFound;		// halo evts id by tower cuts
		};


		// all the histograms


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

		// Setters
		void SetEmTimeWindow(float tmin_,float tmax_) {
				fEmTimeMin = tmin_;
				fEmTimeMax = tmax_;
		}
		void GetEmTimeWindow(float &tmin_,float &tmax_) {
			tmin_ = fEmTimeMin;
			tmax_ = fEmTimeMax;
		}
		float GetEmTimeWindowMin() const { return fEmTimeMin; }
			
		float GetEmTimeWindowMax() const { return fEmTimeMax; }

		// Getters

	private:
		bool 	qMc;  							// true if MonteCarlo data
		Global_Counters_t counter;

		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TagTightPhotons* tightMod;
		TagBeamHalo* tagHaloMod;
		EventProperties* evtPropMod;
		SetEmTimes* setTimeMod;

		//naming convections
		// G = Met/Et>1  L = Met/Et<1
		// GA_ = 1 jet case
		// GB_ = 2 jet case etc

		// for MEtEtRatio > 1  : halo type 0-5
		std::vector<Histograms::EventHists_t> GA_Evt, GB_Evt;
		std::vector<Histograms::PhotonHists_t> GA_Pho, GB_Pho;
		std::vector<Histograms::JetHists_t> GA_Jet, GB_Jet1, GB_Jet2;
		std::vector<Histograms::Photon1JetHists_t> GA_PhoJet, GB_PhoJet1, GB_PhoJet2;
		std::vector<Histograms::Photon2JetsHists_t> GB_PhoJets;		//jet>=2 case, pho+2jets
		std::vector<Histograms::TwoJetsHists_t> GB_TwoJets;

		// for MEtEtRatio < 1  : halo type 0-5
		std::vector<Histograms::EventHists_t> LA_Evt, LB_Evt;
		std::vector<Histograms::PhotonHists_t> LA_Pho, LB_Pho;
		std::vector<Histograms::JetHists_t> LA_Jet, LB_Jet1, LB_Jet2;
		std::vector<Histograms::Photon1JetHists_t> LA_PhoJet, LB_PhoJet1, LB_PhoJet2;
		std::vector<Histograms::Photon2JetsHists_t> LB_PhoJets;		//jet>=2 case, pho+2jets
		std::vector<Histograms::TwoJetsHists_t> LB_TwoJets;

		
		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;
		
		bool bRunPermit;	
		float fEmTimeMin, fEmTimeMax;			// em time window
		TH1F* hMetEtRatio;

	ClassDef(HaloStudyAllWithMetEt,1)
};

#endif
