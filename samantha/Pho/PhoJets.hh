#ifndef PHOJETS_HH
#define PHOJETS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/JetFilterModuleV3.hh"
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/obj/Histograms.hh"
#include "samantha/Pho/EventProperties.hh"
#include "samantha/Pho/TagPMTSpikes.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "samantha/Pho/TagBeamHalo.hh"
#include "samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/Pho/SetEmTimes.hh"
#include <fstream>
#include "Stntuple/obj/TGenpBlock.hh"
#include <Stntuple/obj/TStnJetBlock.hh>
#include "TTree.h"

class PhoJets: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TGenpBlock* fGenpBlock;
		TStnJetBlock*      fJetBlockClu04;


	public:
		PhoJets(const char* name="PhoJets", const char* title="PhoJets");
		~PhoJets();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TGenpBlock* GetGenpBlock() { return fGenpBlock; }


				/************** abstract objects *************************/

		// enumerators to count the photon events
			enum COUNTS {
				EVENTS_PROCESSED = 0,
				EVENTS_PASSED    = 1,
				TIGHT_PHO_EVENTS   = 2,				//inclusive all events with a tight photon
				TIGHT_PHO1JET_EVENTS   = 3,
				TIGHT_PHO2JET_EVENTS   = 4,
				ONE_TIGHT_PHO_ONLY = 5,
				TWO_TIGHT_PHO_ONLY = 6,
				ONE_TIGHT_ANY_LOOSE_PHO = 7,
			};
			enum PHOTON_TYPE {
				iSignalPhoton = 0,
				iSidebandPhoton = 1,
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


		bool DoJetSelection(int iSpIndex);

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

		void FillPhotonIDHist(const int iSpInd);
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }

	private:
		bool 	qMc;  				// true if MonteCarlo data
		InitSuperPhotons* initSpMod;
		TagTightPhotons* tightMod;
		EventProperties* evtPropMod;
		TagPMTSpikes* pmtMod;
		TagElectrons* eleMod;
		TagBeamHalo* haloMod;
		SetEmTimes* setEMTMod;
		TagPhoenixPhotons* phoenixMod;
		JetFilterModuleV3 *jetMod;

		bool bPermitRun;			// make sure all pre-req are met before run

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		float fEmTimeMin, fEmTimeMax;			// em time window
		int iPhotonType;     //0 == signal (tight photon) 1 == sideband photon
		
		float fPhoMaxDetEta;   //absolute upper limit on Det Eta of the electorn
		float fPhoMinDetEta;   //absolute lower limit on Det Eta of the electron
		float fPhoMinEt;			//Min Et of the photon
		float fJetMinEt;			//Min Et of the jet
		float fJetMaxDetEta;			//Max Eta of the jet
		int bRejPhxPhos;
		int iHaloType;
		float fMinMEt; //minimum MET
		
		SuperPhoton* thePhoton;  // this is THE photon passed all selections cuts
		bool bNoSummary;			//control the end job summary print
		TH1F* metb4;
		TH1F* meta4;
		TH1F* phometdelphi_a4;
		TH1F* jetmetdelphi_a4;


	public:
		int GetHaloType() const { return iHaloType; }
		void SetHaloType(const int iT); 
		void SetPhotonType(const int type);
		int GetPhotonType() const { return iPhotonType; }
		std::string GetPhotonTypeLabel() const
		{
			if (GetPhotonType() == 0) return "Signal Photon";
			else if (GetPhotonType() == 1) return "Sideband Photon";
			else 
			{
				return "??? YOU TELL ME! ???????";
			}
		};


		float GetPhoMaxDetEta() const { return fPhoMaxDetEta; }
		float GetPhoMinDetEta() const { return fPhoMinDetEta; }
		float GetPhoMinEt() const { return fPhoMinEt; }
		void SetPhoMaxDetEta(const float eta);
		void SetPhoMinDetEta(const float eta);
		void SetPhoMinEt(const float et);
	
		float GetJetMaxDetEta() const { return fJetMaxDetEta; }
		float GetJetMinEt() const { return fJetMinEt; }
		void SetJetMaxDetEta(const float eta) { fJetMaxDetEta = fabs(eta); }
		void SetJetMinEt(const float et) { fJetMinEt = et; }

		void SetRejPhxPhos(const bool rej) { bRejPhxPhos = rej; }
		bool GetRejPhxPhos() const { return bRejPhxPhos; }

		float GetMinMEt() const { return fMinMEt; }
		void SetMinMEt(const float met) { fMinMEt = met; }

		SuperPhoton* GetPhoton() const;
		
	ClassDef(PhoJets,1)
};

#endif
