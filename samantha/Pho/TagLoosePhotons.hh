#ifndef TAGLOOSEPHOTONS_HH
#define TAGLOOSEPHOTONS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/obj/SuperPhoton.hh"
#endif

class TagLoosePhotons: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;

		static const int iNCenLoosePhoIdBits = 5;			// these are loose central photon ID cuts
		enum LooseIDBits_t {  // bit i will be set when the cut is failed
			kCentralBit 	= 0x1 << 0,
			kEtCorrBit		= 0x1 << 1,
			kXCesBit			= 0x1 << 2,
			kZCesBit			= 0x1 << 3,
			kChi2MeanBit   = 0x1 << 4,
			kHadEmBit 		= 0x1 << 5,
			kIso4Bit			= 0x1 << 6,
			kN3dBit 			= 0x1 << 7,
			kTrkPtBit		= 0x1 << 8,
			kTrkIsoBit		= 0x1 << 9
		};
	
		static const int iNPlugLoosePhoIdBits = 5;			// these are loose plug photon ID cuts
		enum PlugLooseIDBits_t {  // bit i will be set when the cut is failed
			kPlugDetBit 		= 0x1 << 0,
			kPlugEtCorrBit		= 0x1 << 1,
			kPlugHadEmBit 		= 0x1 << 2,
			kPlugIso4Bit		= 0x1 << 3,
			kPlugTrkIsoBit		= 0x1 << 4
		};
	

	public:
		TagLoosePhotons(const char* name="TagLoosePhotons", const char* title="TagLoosePhotons");
		~TagLoosePhotons();

	
		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int oneloosephoton;
			unsigned int twoloosephotons;
			unsigned int morethan2loosephotons;
			unsigned int PhoEvts;	// events with at least one loose photon
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
		int GetCenLooseIdWord(SuperPhoton*);
		int GetPlugLooseIdWord(SuperPhoton*);
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }

	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		bool bRunPermit;		//to make sure we have all dependecies
		InitSuperPhotons *initSpMod;
		const int iNLooseBits;
		bool bNoSummary;			//control the end job summary print

	
	ClassDef(TagLoosePhotons,1)
};

#endif
