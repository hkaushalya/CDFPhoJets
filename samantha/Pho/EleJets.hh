#ifndef ELEJETS_HH
#define ELEJETS_HH

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
#include "samantha/Pho/TriggerModule.hh"
#include "samantha/Pho/SetEmTimes.hh"
#include <fstream>

class EleJets: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		EleJets(const char* name="EleJets", const char* title="EleJets");
		~EleJets();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/



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
		float GetEleMaxDetEta() const { return fEleMaxDetEta; }
		float GetEleMinDetEta() const { return fEleMinDetEta; }
		float GetEleMinEt() const { return fEleMinEt; }
		float GetEleMaxEt() const { return fEleMaxEt; }
		void SetEleMaxDetEta(const float eta);
		void SetEleMinDetEta(const float eta);
		void SetEleMinEt(const float et);
		void SetEleMaxEt(const float et);
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
		int GetHaloType() const { return iHaloType; }
		void SetHaloType(int iT) { 
			if (iT <0 || iT >5) {
				StdOut(__FILE__,__LINE__,3,"No such halo type. pl. pick from 0-5.");
			} else {
				iHaloType = iT;
			}
		}
		void SetRemoveConvEles(const bool b) { bRemConvEles = b; }
		bool GetRemoveConvEles() const { return bRemConvEles; }

		void FillPhotonIDHist(const int iSpInd);
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
		bool GetUseNvtxWeigts() const { return bUseNvtxWgts; }
		void UseNvtxWeigts(const bool b) { bUseNvtxWgts =  b; }
		void SetDataSample(const int s) { iDataSample =  s; }
		int GetDataSample() const { return iDataSample; }
		void SetNvtxWeights();

	private:
		bool 	qMc;  				// true if MonteCarlo data
		InitSuperPhotons* initSpMod;
		JetFilterModuleV3 *jetMod;
		TagTightPhotons* tightMod;
		EventProperties* evtPropMod;
		TagPMTSpikes* pmtMod;
		TagElectrons* eleMod;
		TagBeamHalo* haloMod;
		SetEmTimes* setEMTMod;
		TagPhoenixPhotons* phoenixMod;
		TriggerModule*  trigMod;

		bool bPermitRun;			// make sure all pre-req are met before run

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		int iHaloType;
		float fEmTimeMin, fEmTimeMax;			// em time window
		float fEleMaxDetEta;   //absolute upper limit on Det Eta of the photon
		float fEleMinDetEta;   //absolute lower limit on Det Eta of the photon
		float fEleMinEt;			//Min Et of the photon
		float fEleMaxEt;			//Min Et of the photon
		bool bRemConvEles;			// 1=Remove conversion electrons

		bool bNoSummary;			//control the end job summary print
		TH1F* metb4;
		TH1F* meta4;
		TH1F* phometdelphi_a4;
		TH1F* jetmetdelphi_a4;
		TH1F* jetEt_a4;
		TH1F* phoEt_a4;
		TH1F* phoMetTrnMass_a4;
		TH1F* phoMetTrnMass_b4;
		TH1F* phoLeadJetInvMass_a4;
		TH1F* Nvtx12_a4;
		std::vector<float> vNvtxWgts;
		int iDataSample;
		bool bUseNvtxWgts;

		
	ClassDef(EleJets,1)
};

#endif
