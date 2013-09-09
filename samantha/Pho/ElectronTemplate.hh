#ifndef ELECTRONTEMPLATE_HH
#define ELECTRONTEMPLATE_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/utils/FreeFunctions.hh"
#include <TMath.h>
#include "samantha/Pho/TagConvElectrons.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include <fstream>
#include "TH1F.h"
#include "samantha/obj/Histograms.hh"
#include "samantha/Pho/EventProperties.hh"
#include "samantha/Pho/LooseElectronTemplate.hh"

class ElectronTemplate: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		ElectronTemplate(const char* name="ElectronTemplate", const char* title="ElectronTemplate");
		~ElectronTemplate();

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
		void FillDataBlocks(int ientry);
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

		void SetGoodRunBit(int bit) {	iGoodRunBit = bit;}
		void SetPhoenixBit(int bit) { iPhoenixBit = bit;}
		int GetGoodRunBit() const { return iGoodRunBit; } 
		int GetPhoenixBit() const { return iPhoenixBit; }
		float GetSumFakeRate_1Jet() const { return fSumFakeRate_1Jet; }
		float GetSumFakeRate_2Jet() const { return fSumFakeRate_2Jet; }


	private:
		TagElectrons* tagEleMod;
		TagConvElectrons* tagConvEleMod;
		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		EventProperties* evtPropMod;
		LooseElectronTemplate* looseEleTemp; //i am going to pass this module 100% and 
														// let loose to pass/fail. just so i don't forget this!
		
		std::string sErrs;
		
		bool bRunPermit;	

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		float fEmTimeMax, fEmTimeMin;			// em time window

		ofstream *ofFRlog_1Jet, *ofFRlog_2Jet;		//write the fake rate to files
		float fSumFakeRate_1Jet;	//total sum of all e+>=1jet events
		float fSumFakeRate_2Jet;	//total sum of all e+>=2jet events
		int iGoodRunBit;				// don't remove although i don't need it. bcos the function is complicated
											// and need careful thinking before removing this!
		int iPhoenixBit;				//if true, calcaulate fake after Phoenix rejection of photons	


	ClassDef(ElectronTemplate,1)
};

#endif
