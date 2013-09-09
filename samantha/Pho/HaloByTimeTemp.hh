#ifndef HALOBYTIMETEMP_HH
#define HALOBYTIMETEMP_HH

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
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/Pho/TagBeamHalo.hh"
#include "samantha/obj/Histograms.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/SetEmTimes.hh"

class HaloByTimeTemp: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		HaloByTimeTemp(const char* name="HaloByTimeTemp", const char* title="HaloByTimeTemp");
		~HaloByTimeTemp();

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
		void BookEvtHistograms(EventHist_t& Hist,std::string sFoldName,
											std::string text2add); 
		void FillHistogram(EventHist_t& Hist, SuperPhoton* sp);
		void BookPhiWedgeHistograms(PhiWedgeHist_t& Hist,
										TFolder* FoldName, std::string text2add);
		void FillPhiWedgeHistogram(PhiWedgeHist_t& Hist, SuperPhoton* sp);

		// setters
		void SetHaloScenario(int hs_) {
			if (hs_ >=0 && hs_ <=5) iHaloScenario = hs_;
			else StdOut(__FILE__,__LINE__,3,"Invalid Halo Scenario specified. valids are 0-5");
		}


		//accessors
		int GetHaloScenario() const { return iHaloScenario; }

	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TagTightPhotons* tightMod;
		TagBeamHalo* haloMod;
		bool bRunPermit;			//make sure all dependendt mods are present
		EventHist_t aHaloByTime; 
		Histograms::BeamHaloHists_t notHalo,halo13_6;		// all that did not fall in to any halo time window
													//halos from 1ns time increments
													// halo13 = time >-13ns & < -12ns ... halo5 = >-5ns & < -4ns
		Histograms::BeamHaloHists_t halo13,halo12,halo11,halo10,halo9,halo8,halo7,halo6,halo5,halo4;
		PhiWedgeHist_t halo13pw023,halo12pw023,halo11pw023,halo10pw023,halo9pw023;	//halos in Phi wedge 0  & 23
		PhiWedgeHist_t halo8pw023,halo7pw023,halo6pw023,halo5pw023,halo4pw023;
		PhiWedgeHist_t halo13pw122,halo12pw122,halo11pw122,halo10pw122,halo9pw122;	//halos in Phi wedge >=1  & <=22
		PhiWedgeHist_t halo8pw122,halo7pw122,halo6pw122,halo5pw122,halo4pw122;
		PhiWedgeHist_t halo136pw122, halo136pw023;
		int iHaloScenario;			// what halo scenario to be picked

		SetEmTimes *setTimeMod;

	ClassDef(HaloByTimeTemp,1)
};

#endif
