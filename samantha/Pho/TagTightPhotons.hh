/*
what this should do
have all the photons id cuts
sets tight id for the superphotons in InitSuperPhotons module

todo
return id-word for a given TStnPhoton
create histos for all cut values
*/


#ifndef TAGTIGHTPHOTONS_HH
#define TAGTIGHTPHOTONS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#include <iostream>
#include <string>
#include "Stntuple/loop/TStnModule.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#endif


class TagTightPhotons: public TStnModule {

	protected:

		static const int iNCenTightPhoIdBits = 12;			// these are tight central photon ID cuts
		enum TightIDBits_t {  // bit i will be set when the cut is failed
			kCentralBit 	= 0x1 << 0,
			kEtCorrBit		= 0x1 << 1,
			kXCesBit			= 0x1 << 2,
			kZCesBit			= 0x1 << 3,
			kHadEmBit 		= 0x1 << 4,
			kIso4Bit			= 0x1 << 5,
			kChi2MeanBit	= 0x1 << 6,
			kN3dBit			= 0x1 << 7,
			kTrkPtBit		= 0x1 << 8,
			kTrkIsoBit		= 0x1 << 9,
			kCesWireE2Bit 	= 0x1 <<10,
			kCesStripE2Bit = 0x1 <<11
		};
	
		static const int iNPlugTightPhoIdBits = 8;			// these are tight plug photon ID cuts
		enum PlugTightIDBits_t {  // bit i will be set when the cut is failed
			kPlugDetBit 		= 0x1 << 0,
			kPlugEtCorrBit		= 0x1 << 1,
			kPlugHadEmBit 		= 0x1 << 2,
			kPlugIso4Bit		= 0x1 << 3,
			kPlugChi2Bit		= 0x1 << 4,
			kPlugPes5by9UBit	= 0x1 << 5,
			kPlugPes5by9VBit	= 0x1 << 6,
			kPlugTrkIsoBit		= 0x1 << 7
		};
	
		
	public:
		TagTightPhotons(const char* name="TagTightPhotons", const char* title="TagTightPhotons");
		~TagTightPhotons();


		//********** define abstract data types here **************
		struct GlobalCounters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int oneCEMtightphoton;
			unsigned int twoCEMtightphotons;
			unsigned int onePLUGtightphoton;
			unsigned int twoPLUGtightphotons;
			unsigned int morethan2CEMtightphotons;
			unsigned int morethan2PLUGtightphotons;
			unsigned int PhoEvts;	// events with at least one tight photon
		};	
		
		struct IdCutsPlots_t {
			TH1F* Cummulative;			//cummulative tight cuts histogram of all photos
			TH1F* Sum;						//plot of all id cuts passed by all photons
		};

		struct PhotonPlots_t {
			TH1F* Detector;
			TH1F* DetEta;
			TH1F* DetPhi;
			TH1F* EtCorr;
			TH1F* XCes;
			TH1F* ZCes;
			TH1F* HadEm;
			TH1F* Chi2Mean;
			TH1F* N3d;
			TH1F* Iso4;
			TH1F* TrkPt;
			TH1F* TrkIso;
			TH1F* Ces2Wire;
			TH1F* Ces2Strip;
			TH1F* HadTDCtime;
		};

		//********** END define abstract data types here **************


		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();
		
		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
		int GetCenTightIdWord(SuperPhoton*);
		void BookIdCutsHistograms(IdCutsPlots_t&, std::string folder_name,
								std::string folder_title);
		void BookPhotonHistograms(PhotonPlots_t&, std::string folder_name,
								std::string folder_title);
		void FillPhotonsHist(PhotonPlots_t&, SuperPhoton* sp);
		void FillIdCutHist(int& iIdWord);
		int GetPlugTightIdWord(SuperPhoton* sp);
		// ******* setters

		// ******** getters
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
		void SetDebug(const bool b) { bDebug = b; }
		bool Debug() const { return bDebug; }
		

	private:
		GlobalCounters_t counter;
		PhotonPlots_t hTightPho;	//plots of tight photons only
		IdCutsPlots_t hIdCut;
		bool bRunPermit;		//to make sure we have all dependecies
		InitSuperPhotons *initSpMod;
		const int iNTightBits;
		bool bNoSummary;			//control the end job summary print
		bool bDebug;				//prints more info for debugging
	ClassDef(TagTightPhotons,0)
};

#endif
