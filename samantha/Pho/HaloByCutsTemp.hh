#ifndef HALOBYCUTSTEMP_HH
#define HALOBYCUTSTEMP_HH

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
#include "samantha/Pho/TagBeamHalo.hh"
#include "samantha/obj/Histograms.hh"
#include "samantha/Pho/HaloByTimeTemp.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/SetEmTimes.hh"

class HaloByCutsTemp: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		HaloByCutsTemp(const char* name="HaloByCutsTemp", const char* title="HaloByCutsTemp");
		~HaloByCutsTemp();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int haloEvtsFound;
		};

		struct EventHist_t {
			TH1F* PhoEmTime ;
			TProfile* PhoEmTimeVsRun;
			TH1F* Met;
			TH1F* Sumet;
		};
		struct PhiWedgeHist_t {				// to look at stuff in the Halo Phi wedge plot
			TH1F* Met;
			TH1F* PhoMetDelPhi;
			TH1F* PhoEtMetRatio;
			TH1F* PhoEmTime ;
			TH1F* MetPhi;
			TH1F* PhoPhi;
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
		void BookEvtHistograms(EventHist_t& Hist,std::string sFoldName, std::string text2add); 
		void FillHistograms(EventHist_t& Hist, SuperPhoton* sp);
		void BookPhiWedgeHistograms(PhiWedgeHist_t& Hist,
										TFolder* FoldName, std::string text2add);
		void FillPhiWedgeHistogram(PhiWedgeHist_t& Hist, SuperPhoton* sp);


	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TagTightPhotons* tightMod;
		TagBeamHalo *haloMod;
		SetEmTimes *setTimeMod;

		EventHist_t aHaloByCut; 
	
		bool bRunPermit;			//make sure we have all dependencies
		
		//get and fill histos using HaloByTimeTemp mod
		Histograms::BeamHaloHists_t halo1js0,halo1js1,halo1js2,halo1js3,halo1js4,halo1js5;
		Histograms::BeamHaloHists_t halo2js0,halo2js1,halo2js2,halo2js3,halo2js4,halo2js5;
		PhiWedgeHist_t halo1js0pw023,halo1js1pw023,halo1js2pw023,halo1js3pw023,halo1js4pw023,halo1js5pw023;	//halos in Phi wedge 0  & 23
		PhiWedgeHist_t halo1js0pw122,halo1js1pw122,halo1js2pw122,halo1js3pw122,halo1js4pw122,halo1js5pw122;	//halos in Phi wedge 0  & 23
		PhiWedgeHist_t halo2js0pw023,halo2js1pw023,halo2js2pw023,halo2js3pw023,halo2js4pw023,halo2js5pw023;	//halos in Phi wedge 0  & 23
		PhiWedgeHist_t halo2js0pw122,halo2js1pw122,halo2js2pw122,halo2js3pw122,halo2js4pw122,halo2js5pw122;	//halos in Phi wedge 0  & 23


	ClassDef(HaloByCutsTemp,1)
};

#endif
