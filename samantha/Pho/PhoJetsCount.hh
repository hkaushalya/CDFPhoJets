#ifndef PHOJETSCOUNT_HH
#define PHOJETSCOUNT_HH

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
#include <Stntuple/obj/TStnJetBlock.hh>
#include "TTree.h"
#include "samantha/Pho/TriggerModule.hh"

class PhoJetsCount: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TGenpBlock* fGenpBlock;
		TStnJetBlock*      fJetBlockClu04;


	public:
		PhoJetsCount(const char* name="PhoJetsCount", const char* title="PhoJetsCount");
		~PhoJetsCount();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TGenpBlock* GetGenpBlock() { return fGenpBlock; }


				/************** abstract objects *************************/

		// enumerators to count the photon events
			enum COUNTS {
				EVENTS_PROCESSED = 0,
				PASS_GOORUN = 1,
				PASS_TRIGGER = 2,
				PASS_NVTX = 3,
				PASS_VTXZ = 4,
				PASS_ET_TRIGMATCHED = 5,
				PASS_TIGHTPHOID = 6, //AT LEAST 1 PHOTON IN THE EVENT
				PASS_NOT_PMT = 7,
				PASS_EMTIME = 8,
				PASS_NO_MUON_STUBS = 9,
				PASS_NOT_PHOENIX = 10,
				PASS_NOT_BH5 = 11,
				//PASS_JET_ET = 12,
				//PASS_JET_ETA = 13,
				PASS_1NJET15 = 12,
				PASS_DELPHI_JETMET = 13,
				PASS_MET_20 = 14,
				PASS_MET_50 = 15,
				
				EVENTS_PASSED    = 16,
				TIGHT_PHO_EVENTS   = 17,				//inclusive all events with a tight photon
				ONE_TIGHT_PHO_ONLY = 18,
				TWO_TIGHT_PHO_ONLY = 19,
				ONE_TIGHT_ANY_LOOSE_PHO = 20,
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


		void DoJetSelection(int iSpIndex);

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

		void FillPhotonIDHist(Histograms::PhotonHists_t& hist, const int iSpInd);
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
		TriggerModule* trigMod;

		bool bPermitRun;			// make sure all pre-req are met before run

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJetsCount;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		Histograms::EventHists_t h1_Evt_a4, h2_Evt_a4;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho_a4, h2_Pho_a4;
		Histograms::JetHists_t h1_Jet_a4, h2_Jet1_a4, h2_Jet2_a4;
		Histograms::Photon1JetHists_t h1_PhoJet_a4, h2_PhoJet1_a4, h2_PhoJet2_a4;
		Histograms::Photon2JetsHists_t h2_PhoJetsCount_a4;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets_a4;


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
		
		bool bNoSummary;			//control the end job summary print
		SuperPhoton* thePhoton; //leading tight photon

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
		void SetPhoMaxDetEta(const double eta);
		void SetPhoMinDetEta(const double eta);
		void SetPhoMinEt(const float et);

		float GetJetMaxDetEta() const { return fJetMaxDetEta; }
		float GetJetMinEt() const { return fJetMinEt; }
		void SetJetMaxDetEta(const float eta) { fJetMaxDetEta = fabs(eta); }
		void SetJetMinEt(const float et) { fJetMinEt = et; }

	
		void SetRejPhxPhos(const bool rej) { bRejPhxPhos = rej; }
		bool GetRejPhxPhos() const { return bRejPhxPhos; }
		SuperPhoton* GetPhoton() const;
		
		bool PassDetEta(const SuperPhoton *sp);
		bool PassEtTrigMatched(const SuperPhoton *sp);
		bool PassXCes(const SuperPhoton* sp);
		bool PassZCes(const SuperPhoton* sp);
		bool PassHadEm(const SuperPhoton* sp);
		bool PassIso(const SuperPhoton *sp);
		bool PassChi2Mean(const SuperPhoton* sp);
		bool PassN3d(const SuperPhoton* sp);
		bool PassTrkPt(const SuperPhoton* sp);
		bool PassTrkIso(const SuperPhoton* sp);
		bool PassCes(const SuperPhoton* sp);
		void FillPhotonHists(Histograms::PhotonHists_t& hist, const SuperPhoton* sp, const float fWeight=1);


		ClassDef(PhoJetsCount,1)
};

#endif
