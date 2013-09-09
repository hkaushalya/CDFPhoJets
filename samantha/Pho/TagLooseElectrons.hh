#ifndef TAGLOOSEELECTRONS_HH
#define TAGLOOSEELECTRONS_HH

/* $Id: TagLooseElectrons.hh,v 1.6 2009/11/10 19:15:56 samantha Exp $
 * $Log: TagLooseElectrons.hh,v $
 * Revision 1.6  2009/11/10 19:15:56  samantha
 * ADDED: more counters to count the different types of electrons (pholike,std
 * etc).
 * RENAMED: GetElectronIdWord() -> GetPhoLikeEleIdWord()
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
#include "Stntuple/obj/TStnElectron.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include "samantha/Pho/TriggerModule.hh"
#endif


class TagLooseElectrons: public TStnModule {

	protected:

		TStnHeaderBlock*    fHeaderBlock;
		TStnVertexBlock*    fVertexBlock;
		TStnElectronBlock*  fElectronBlock;
		TStnTrackBlock*     fTrackBlock;
		TStnMetBlock*       fMetBlock;


		static const int iNEleIdBits = 16;			// these are photon-like electron ID cuts
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


		static const int iNStdLoosePlugEleIdBits = 10;		// these are std. electron ID cuts
		enum StdLoosePlugEleIDBits_t {  // bit i will be set when the cut is failed
			kSLP_DetBit 		= 0x1 << 0,
			kSLP_EtBit			= 0x1 << 1,
			kSLP_PesEtaBit		= 0x1 << 2,
			kSLP_HadEmBit 		= 0x1 << 3,
			kSLP_Pem3x3Bit		= 0x1 << 4,
			kSLP_Chi2ThreeBit	= 0x1 << 5,
			kSLP_Pes5x9UBit	= 0x1 << 6,
			kSLP_Pes5x9VBit	= 0x1 << 7,
			kSLP_PemIsoFrBit	= 0x1 << 8,
			kSLP_PesDelRBit	= 0x1 << 9
		};

		static const int iNStdLooseCenEleIdBits = 10;		// these are std. electron ID cuts
		enum StdLooseCenEleIDBits_t {  // bit i will be set when the cut is failed
			kSLC_DetBit 		= 0x1 << 0,
			kSLC_EtBit			= 0x1 << 1,
			kSLC_FidBit			= 0x1 << 2,
			kSLC_TrkBit			= 0x1 << 3,
			kSLC_TrkZ0Bit		= 0x1 << 4,
			kSLC_TrkPtBit		= 0x1 << 5,
			kSLC_TrkNaslBit	= 0x1 << 6,
			kSLC_TrkNsslBit	= 0x1 << 7,
			kSLC_HadEmBit 		= 0x1 << 8,
			kSLC_IsoFrBit		= 0x1 << 9
		};


		
	public:
		TagLooseElectrons(const char* name="TagLooseElectrons", const char* title="TagLooseElectrons");
		~TagLooseElectrons();
	
		TStnVertexBlock*   GetVertexBlock()  { return fVertexBlock;  }
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnTrackBlock*    GetTrackBlock()   { return fTrackBlock;   }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnMetBlock*      GetMetBlock() 	{ return fMetBlock; }
		
		/*************** define all abstract object here ****************/
		
		struct GlobalCounters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int electronEvents;	// PhoLike or STD 
			unsigned int PhoLikeEleEvts;	//Events that has at least one pho-like electron
			unsigned int PhoLikeEles;//total # of pho-like electrons found
			unsigned int StdEleEvts;		//Events that has at lease one std. electron
			unsigned int StdCEMEles;	//total # of central std. electrons
			unsigned int StdPLUGEles;//total # of plug std. electrons
		};	

		struct IDhists_t {
			TH1F* StdCen;
			TH1F* StdPlug;
			TH1F* CemIdCum;
			TH1F* PlugIdCum;
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
		
		unsigned int GetStdPlugLooseEleID(TStnElectron* pho, TStnEvent* event, const int Nvtx);
		unsigned int GetStdCenLooseEleID(TStnElectron* Ele, TStnTrackBlock*, const int Nvtx);
		double CorrIsoFr(TStnElectron* ele, const int Nvtx);
		double CorrEnergyEle(TStnElectron* ele);
		double CorrEtEle(TStnElectron* ele);
		double CorrIso(TStnElectron* ele, const int Nvtx);
		TStnTrack* GetPnxTrack(TStnEvent* event, TStnElectron* ele);
		void FillStdElectronIDcutPlots(SuperPhoton* sp);


		
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
		TriggerModule* trigMod;
		Histograms Hist;
		Histograms::EmObjHists_t EleB4Cut, EleA4Cut;
		bool bRunPermit;		//check for all dependencies

		float fPemMinEt;
		float fPemMaxEt;
		float fPemMinEtaPes; 
		float fPemMaxEtaPes;
		float fPemMinHADEM; 
		float fPemMaxHADEM;       
		float fPem3x3FitCut; // cut on PEM 3x3 tower fit status code
		float fPemMin3x3Chi2;  
		float fPemMax3x3Chi2;
		float fPemMinPes5x9U; 
		float fPemMaxPes5x9U; 
		float fPemMinPes5x9V; 
		float fPemMaxPes5x9V;
		float fPemMinIsoFr;  
		float fPemMaxIsoFr;
		float fPemMinPesdR; 
		float fPemMaxPesdR;

		float fCemMinEt;    
		float fCemMaxEt; 
		float fCemMinTrkZ0; 
		float fCemMaxTrkZ0;
		float fTrkMinPt;
		float fTrkMaxPt; 
		float fCemMinCotAxSl; // cut on number of COT Axial SL with 5 or more hits 
		float fCemMaxCotAxSl;
		float fCemMinCotStSl; // cut on number of COT Stereo SL with 5 or more hits 
		float fCemMaxCotStSl;
		float fCemMinHADEM[2]; 
		float fCemMaxHADEM[2];
		float fCemMinIsoFr;  
		float fCemMaxIsoFr;  
		float fCemFidCut;

		IDhists_t idHist;
		bool bNoSummary;			//control the end job summary print

	ClassDef(TagLooseElectrons,0)
};

#endif
