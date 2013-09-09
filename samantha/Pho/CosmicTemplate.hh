#ifndef COSMICTEMPLATE_HH
#define COSMICTEMPLATE_HH

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
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/obj/Histograms.hh"
#include "samantha/Pho/EventProperties.hh"

class CosmicTemplate: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		CosmicTemplate(const char* name="CosmicTemplate", const char* title="CosmicTemplate");
		~CosmicTemplate();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/



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
		void DoJetSelection(int iSpIndex);

		void SetCosmicTimeWindow(float Tmin, float Tmax) {
			fCosmicTimeMin = Tmin;
			fCosmicTimeMax = Tmax;
		}
		void  GetCosmicTimeWindow(float &tmin_, float &tmax_) const {
			tmin_ = fCosmicTimeMin;
			tmax_ = fCosmicTimeMax;
		}
		float GetCosmicTimeMin() const { return fCosmicTimeMin; }
		float GetCosmicTimeMax() const { return fCosmicTimeMax; }

	private:
		bool 	qMc;  				// true if MonteCarlo data
		bool bPermitRun;			// make sure all pre-req are met before run
		float fCosmicTimeMin;
		float fCosmicTimeMax;

		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TagTightPhotons* tightMod;
		EventProperties* evtPropMod;
		
		float fCorrection;

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;
		
	ClassDef(CosmicTemplate,1)
};

#endif
