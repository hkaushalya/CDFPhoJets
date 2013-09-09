#ifndef TAGELECTRONS_HH
#define TAGELECTRONS_HH

/* $Id: TagElectrons.hh,v 1.12 2009/11/10 19:29:31 samantha Exp $
 * $Log: TagElectrons.hh,v $
 * Revision 1.12  2009/11/10 19:29:31  samantha
 * Added Auto CVS Log info into the file.
 *
 *
 */

#if !defined (__CINT__) || defined (__MAKECINT__)
#include <iostream>
#include <string>
#include "Stntuple/loop/TStnModule.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "Stntuple/obj/TStnVertexBlock.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnElectronBlock.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/obj/Histograms.hh"
#include "Stntuple/obj/TStnMetBlock.hh"
#endif


class TagElectrons: public TStnModule {

	protected:

		TStnHeaderBlock*    fHeaderBlock;
		TStnVertexBlock*    fVertexBlock;
		TStnElectronBlock*  fElectronBlock;
		TStnTrackBlock*     fTrackBlock;
		TStnMetBlock*       fMetBlock;


		// these are the photon-like electron id variables
		static const int iNEleIdBits = 16;
		enum ElectronIDBits_t {  // bit i will be set when the cut is failed
			kCentralBit 	= 0x1 << 0,
			kEtCorrBit		= 0x1 << 1,
			kHadEmBit 		= 0x1 << 2,
			kChi2Bit			= 0x1 << 3,
			kN3d0Bit			= 0x1 << 4,
			kN3d2Bit			= 0x1 << 5,
			k2ndTrkPtBit	= 0x1 << 6,
			kMatchEleBit	= 0x1 << 7,
			kEleTrkBit     = 0x1 << 8,
			kEleTrkZ0Bit	= 0x1 << 9,	
			kEoverPBit 		= 0x1 <<10,
			kIso4Bit			= 0x1 <<11,
			kTrkIso        = 0x1 <<12,
			k2ndCesBit 	   = 0x1 <<13,
			kXCesBit 		= 0x1 <<14,
			kZCesBit			= 0x1 <<15 
		};

		// here are the std. electron id variables (central and plug)
		static const int iNStdTightEleIdBits = 8;
		enum StdTightCenEleIDbits {
			kSTC_DetectorBit 	= 0x1 << 0,
			kSTC_EtCorrBit 	= 0x1 << 1,
			kSTC_LshrBit 		= 0x1 << 2,
			kSTC_TrkBit     	= 0x1 << 3,
			kSTC_EoverPBit 	= 0x1 << 4,
			kSTC_DelZBit 		= 0x1 << 5,
			kSTC_DelXQBit 		= 0x1 << 6,
			kSTC_Chi2StripBit = 0x1 << 7
		};

		
	public:
		TagElectrons(const char* name="TagElectrons", const char* title="TagElectrons");
		~TagElectrons();
	
		TStnVertexBlock*   GetVertexBlock()  { return fVertexBlock;  }
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnTrackBlock*    GetTrackBlock()   { return fTrackBlock;   }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnMetBlock*      GetMetBlock() 	{ return fMetBlock; }
		
		/*************** define all abstract object here ****************/
		
		struct GlobalCounters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int electronEvents;
		};	


		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();
		
		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
		unsigned int GetPhoLikeEleIdWord(TStnPhoton* Pho, TStnEvent* event,
												TStnTrackBlock* fTrkBlock, double zvx);
		void FillHistograms(Histograms::EmObjHists_t& , SuperPhoton*);
		unsigned int GetStdCenTightEleID(TStnElectron* Ele, TStnTrackBlock* fTrkBlock);
		unsigned int GetStdPlugTightEleID(TStnElectron* Ele, TStnTrackBlock* fTrkBlock);
		
	
		// ******* setters

		void SetMode(int mode=0) {
			if (mode >=0 && mode <3) iMode = mode;
			else {
				iMode = 0;
				StdOut("Invalid mode!, valid 0,1,2. Using default(0) tagging only mode.",GetName()); 
			}
		};

		// ******** getters

		int GetMode() const { return iMode;}
		int GetElectronNIdBits() const { return iNEleIdBits; } 
		
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }

	private:
		GlobalCounters_t counter;
		int iMode;
		InitSuperPhotons* initSp;
		Histograms Hist;
		Histograms::EmObjHists_t EleB4Cut, EleA4Cut;
		bool bRunPermit;		//check for all dependencies

		// std. central electron id cuts values
		float fEtMin;
		float fEtMax;
		float fCemMinLshr;    
		float fCemMaxLshr;
		float fCemMinEP;    
		float fCemMaxEP;
		float fCemMinChi2strip; 
		float fCemMaxChi2strip;
		float fCemMindXQ;
		float fCemMaxdXQ; 
		float fCemMinCesdZ; 
		float fCemMaxCesdZ;
	
		bool bNoSummary;			//control the end job summary print



	ClassDef(TagElectrons,0)
};

#endif
