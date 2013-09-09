#ifndef HALO4_5STUDY_HH
#define HALO4_5STUDY_HH

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


class Halo4_5Study: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		Halo4_5Study(const char* name="Halo4_5Study", const char* title="Halo4_5Study");
		~Halo4_5Study();

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
		void SetHaloScenario(int hsc_) { 
			if (hsc_ < 0 || hsc_ > 5) {
				StdOut(__FILE__,__LINE__,3,
				"Invalid Halo scenario picked. valid choices	are from 0-5. using default scenario.");
			} else iHaloScenario = hsc_;
		}
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
		int GetHaloScenario() const { return iHaloScenario; }

	private:
		bool 	qMc;  							// true if MonteCarlo data
		Global_Counters_t counter;

		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TagTightPhotons* tightMod;
		TagBeamHalo* tagHaloMod;
		EventProperties* evtPropMod;
		SetEmTimes* setTimeMod;

		Histograms::EventHists_t h1_Evt, h2_Evt;						// for MEtEtRatio > 1
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		Histograms::EventHists_t h1_Evt_, h2_Evt_;						// for MetEtRatio < 1
		Histograms::PhotonHists_t h1_Pho_, h2_Pho_;
		Histograms::JetHists_t h1_Jet_, h2_Jet1_, h2_Jet2_;
		Histograms::Photon1JetHists_t h1_PhoJet_, h2_PhoJet1_, h2_PhoJet2_;
		Histograms::Photon2JetsHists_t h2_PhoJets_;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets_;



		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;
		
		bool bRunPermit;	
		int iHaloScenario;				//halo scenario to be used
		float fEmTimeMin, fEmTimeMax;			// em time window

	ClassDef(Halo4_5Study,1)
};

#endif
