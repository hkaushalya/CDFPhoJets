#ifndef COSMICSTUDY_HH
#define COSMICSTUDY_HH

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
#include <Stntuple/obj/TStnEmTimingBlock.hh>
#include "Stntuple/ana/TEmTimingModule.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "CalorGeometry/CalConstants.hh"
#include "Stntuple/data/TCalTower.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnVertexBlock.hh"


class CosmicStudy: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TCalDataBlock*      fCalBlock;
		TStnEmTimingBlock*  fTimingBlock;
		TStnJetBlock*       fJetBlock;
		TStnVertexBlock*    fVertexBlock;

		typedef std::vector<TCalTower*> CalDataArray;
		typedef std::vector<TCalTower*>::iterator CalDataArrayI;


	public:
		CosmicStudy(const char* name="CosmicStudy", const char* title="CosmicStudy");
		~CosmicStudy();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TCalDataBlock* 	 GetCalBlock()  { return fCalBlock;  }
		TStnJetBlock*      GetJetBlock() 	{  return fJetBlock; }
		TStnVertexBlock*   GetVertexBlock()  const { return fVertexBlock;  }


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
		void MatchCalorTowers(int jetInd, CalDataArray* cdholder);
		double  GetEmTiming(const int jetInd);
		void FillDataBlocks(int ientry);
		int NClass12Vertices();
		float BestVertexZ();

	private:
		bool 	qMc;  				// true if MonteCarlo data
		bool bPermitRun;			// make sure all pre-req are met before run
		float fCosmicTimeMin;
		float fCosmicTimeMax;

		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TagTightPhotons* tightMod;
		EventProperties* evtPropMod;
		

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		TH1F *hPho1Jet_Jet1EmTime, *hPho2Jet_Jet2EmTime, *hPho2Jet_Jet1EmTime;
		TH1F *hPho1Jet_PhoJet1DelEmTime, *hPho2Jet_PhoJet2DelEmTime, *hPho2Jet_PhoJet1DelEmTime;
		
	ClassDef(CosmicStudy,1)
};

#endif
