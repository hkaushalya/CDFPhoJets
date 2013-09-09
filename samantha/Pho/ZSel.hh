#ifndef ZSEL_HH
#define ZSEL_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TStnElectronBlock.hh"
#include "Stntuple/obj/TStnMetBlock.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "samantha/Pho/TagConvElectrons.hh"
#include "samantha/Pho/JetFilterModuleV3.hh"
#include "samantha/Pho/TriggerModule.hh"
#endif

class ZSel: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnElectronBlock*  fElectronBlock;
		TStnMetBlock*       fMetBlock;


	public:
		ZSel(const char* name="ZSel", const char* title="ZSel");
		~ZSel();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnMetBlock*      GetMetBlock()     { return fMetBlock; }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int Zevents;
		};


		// all the histograms
		struct ZPlots_t {	//__  plots
			TH1F* ZMass;
			TH1F* ZPt;
			TH1F* ZEleEtratio;
			TH1F* ZSumet0;
			TH1F* ZHt;
			TH1F* ZMet;
			TH1F* ZeleDelPhi;
		};
		struct ElePlots_t { //__ vars from the matching electron
			TH1F* EtCorr;
			TH1F* HadEm;
			TH1F* IsoEtCorr;
			TH1F* TrkPt;
			TH1F* EoverP;
			TH1F* Chi2Mean;
			TH1F* Ces2Wire;
			TH1F* Ces2Strip;
			TH1F* DelPhiEleMet;
		};

		struct AcceptPlots_t {  // acceptance plots
			TH1F* Cumm;
			TH1F* All;
		};


		struct RunDep_Hists_t {			//run dependent histograms counters
			TH1F* Nevts;			// events in the run
			TH1F* NZs;				//z candidates with run number 
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
		void BookZHist(ZPlots_t& hist, std::string sFoldName, std::string sFoldTitle);
		void BookEleHist(ElePlots_t& hist, std::string sFoldName, std::string sFoldTitle, std::string sText);
		void FillZPlots(ZPlots_t& hist, SuperPhoton* sp1, SuperPhoton* sp2, const float fWgt=1);
		void FillElePlots(ElePlots_t& hist, SuperPhoton* sp);
		void BookAcceptanceHist(AcceptPlots_t& hist, std::string sFoldName, std::string sFoldTitle);
		void FillAcceptancePlots(AcceptPlots_t& hist, SuperPhoton* sp);
		void BookRunDependentHistograms(RunDep_Hists_t& hist,
				std::string sFoldName, std::string sFoldTitle);
		void SetMinEleEt(const float et) { fMinEleEt = et; }
		void SetMaxDetEta(const float eta) { fMaxDetEta = eta; }
		float GetMinEleEt() const { return fMinEleEt; }
		float GetMaxDetEta() const { return fMaxDetEta; }

		bool GetUseNvtxWeigts() const { return bUseNvtxWgts; }
		void UseNvtxWeigts(const bool b) { bUseNvtxWgts =  b; }
		void SetDataSample(const int s) { iDataSample =  s; }
		int GetDataSample() const { return iDataSample; }
		void SetNvtxWeights();

	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		ZPlots_t hZPlot;
		ElePlots_t hEle1Plot, hEle2Plot, hAllElePlot;
		InitSuperPhotons* initSpMod;
		TagConvElectrons* tagConvEleMod;
		TagElectrons* tagEleMod;
		TriggerModule*  trigMod;
		TH1F *metb4cut,*meta4cut;
		TH1F *nvtx12b4cut,*nvtx12a4cut;
		AcceptPlots_t  hAccept;
		RunDep_Hists_t hRDHist;
		float fMinEleEt;
		float fMaxDetEta;
		JetFilterModuleV3 *jetMod;
		std::vector<float> vNvtxWgts;
		int iDataSample;
		bool bUseNvtxWgts;
	
	ClassDef(ZSel,1)
};

#endif
