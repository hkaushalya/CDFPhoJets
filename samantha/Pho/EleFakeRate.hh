#ifndef ELEFAKERATE_HH
#define ELEFAKERATE_HH

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
#include "samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include "Stntuple/obj/TStnMetBlock.hh"
#include <fstream>
#include "TH1F.h"
#include "samantha/obj/Histograms.hh"
#include "samantha/Pho/EventProperties.hh"

class EleFakeRate: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnMetBlock*       fMetBlock;


	public:
		EleFakeRate(const char* name="EleFakeRate", const char* title="EleFakeRate");
		~EleFakeRate();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnMetBlock*      GetMetBlock() 	{ return fMetBlock; }


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

		void SetGoodRunBit(int bit) {	iGoodRunBit = bit;}
		void SetPhoenixBit(int bit) { iPhoenixBit = bit;}
		
		float GetSumFakeRate_1Jet() const { return fSumFakeRate_1Jet; }
		float GetSumFakeRate_2Jet() const { return fSumFakeRate_2Jet; }
		int GetGoodRunBit() const { return iGoodRunBit; } 
		int GetPhoenixBit() const { return iPhoenixBit; }


		enum { FkRtParArray = 4 };
		enum { PhnxParArray = 3 };

		double fakerate_par[FkRtParArray]; // fake rate parametrization
		double fakerate_parerr[FkRtParArray]; // fake rate parametrization uncertainty
		double fakerate_corr_coeff[FkRtParArray][FkRtParArray]; // covariance matrix coefficients
		double phnx_par[PhnxParArray]; // phnx eff. parametrization
		double phnx_parerr[PhnxParArray]; // phnx eff. parametrization uncertainty
		double phnx_corr_coeff[PhnxParArray][PhnxParArray]; // covariance matrix coefficients

		void MyFakeRate(double eleEt, int goodrun_bit, int phnx_bit, 
							double &wght_d, double &wght_m, double &wght_p); // returns fake rate & syst
		double MyShiftCorrFunc(double eleEt); // correction for Et(det)/Et(genp) difference
		double MyPhnxEffFunc(double eleEt);  // phoenix efficiency function
		double MyPhnxEffFuncErr(double eleEt);  // uncertainty on phoenix efficiency 
		double MyFakeRateFunc(double eleEt); // returns fake rate
		double MyFakeRateFuncErr(double eleEt); // returns uncertainty on fake rate

	private:
		bool 	qMc;  				// true if MonteCarlo data
		float fSumFakeRate_1Jet;	//total sum of all e+>=1jet events
		float fSumFakeRate_2Jet;	//total sum of all e+>=2jet events
		int iGoodRunBit;				// i don't need this. delete this
		int iPhoenixBit;				//if true calcaulate fake after Phoenix rejection of photons	
		TagElectrons* tagEleMod;
		TagConvElectrons* tagConvEleMod;
		TagPhoenixPhotons* tagPhxMod;
		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		EventProperties* evtPropMod;
		
		std::string sErrs;
		ofstream *ofFRlog_1Jet, *ofFRlog_2Jet;		//write the fake rate to files
		
		bool bRunPermit;	

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

	ClassDef(EleFakeRate,1)
};

#endif
