#ifndef FLATSTUPLEMAKER_HH
#define FLATSTUPLEMAKER_HH
/*{{{*/
/**********************************************************
 * $Id: FlatStupleMaker.hh,v 1.8 2011/05/24 00:42:25 samantha Exp $
 * $Log: FlatStupleMaker.hh,v $
 * Revision 1.8  2011/05/24 00:42:25  samantha
 * MODIFIED: To get the b-tag (SecVtx) jet info.
 * ADDED:1.  #include "BTagObjects/jetMistagJun2009.hh
 * I was using the wrong class to get b-tag info
 * #include "Stntuple/obj/TStnSecVtxTagBlock.hh. So deleted that.
 *
 * Revision 1.7  2010/08/26 15:45:31  samantha
 * ADDED: 1. Jet EM timing.
 *        2. EM obj threshold can be set.
 * MODIFIED: Replaced jetmod2 with jetmod3 which is stripdown version generating
 * only the jets/met etc I need.
 *
 **********************************************************/
/*}}}*/
///////////////////////////////////////////////////////////
// See sourcce file for class description                //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#if !defined (__CINT__) || defined (__MAKECINT__)

#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TTree.h"
#include <string>
#include <iostream>
#include <vector>
#include "samantha/obj/Stuple.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/JetFilterModuleV3.hh"
#include "samantha/Pho/TriggerModule.hh"
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "Stntuple/obj/TStnElectronBlock.hh"
#include "Stntuple/obj/TStnJetBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnVertexBlock.hh"
#include "samantha/Pho/TagBeamHalo.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/photon/TCesCpr.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "Stntuple/obj/TStnSecVtxTagBlock.hh"
#include "Stntuple/ana/TEmTimingModule.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "CalorGeometry/CalConstants.hh"
#include "Stntuple/data/TCalTower.hh"
#include "../../RootTools/IOColors.hh"
#include "BTagObjects/jetMistagJun2009.hh"

class FlatStupleMaker: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()
		TStnElectronBlock*  fElectronBlock;
		TStnJetBlock*       fJetBlock;
		TStnHeaderBlock*    fHeaderBlock;
		TStnVertexBlock*    fVertexBlock;
		TStnTrackBlock*     fTrackBlock;
		TGenpBlock* 		  fGenpBlock;
		TStnEmTimingBlock*  fTimingBlock;
		TCalDataBlock*      fCalBlock;
		TStnSecVtxTagBlock* fSecVtxTagBlock;

		typedef std::vector<TCalTower*> CalDataArray;
		typedef std::vector<TCalTower*>::iterator CalDataArrayI;
	public:
		FlatStupleMaker(const char* name="FlatStupleMaker", const char* title="FlatStupleMaker");
		~FlatStupleMaker();

		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnJetBlock*      GetJetBlock()  	 { return fJetBlock;     }
		TStnHeaderBlock* 	 GetHeaderBlock()  { return fHeaderBlock;  }
		TStnVertexBlock*   GetVertexBlock()  { return fVertexBlock;  }
		TStnSecVtxTagBlock* GetSecVtxTagBlock()  { return fSecVtxTagBlock;  }
		TStnTrackBlock*    GetTrackBlock()   { return fTrackBlock;   }
		TGenpBlock* 		 GetGenpBlock() 	 { return fGenpBlock;    }


				/************** abstract objects *************************/
		struct Global_Counters_t {
			unsigned int evtsProcessed;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int evtsWritten;
			unsigned int phoOnlyEvts;
			unsigned int eleOnlyEvts;
			unsigned int phoEleEvts;
		};

		struct Hist_t {
			TH1F* eleEt;
			TH1F* phoEt;
			TH1F* Npho;		//number of photons per event (tight ot loose)
			TH1F* Nele;		//number of electrons per event (tight ot loose)
			TH1F* Ncj;		// number of central jets
			TH1F* Nupj;		//number of jes up jets
			TH1F* Ndownj;	// number if jes down jets
			TH2F* NcjUpj;		//N central jets Vs N jes up jets
			TH2F* NcjDownj;	// N central jets Vs N jes down jets
			TH1F* cj_EperJet;	//sum of jet et / n jet;
			TH1F* upj_EperJet;
			TH1F* downj_EperJet;

			TH1F* w_met_c;
			TH1F* w_tmass_c;
			TH1F* w_met_p;
			TH1F* w_tmass_p;
			TH1F* z_met_cc;
			TH1F* z_mass_cc;
			TH1F* z_met_cp;
			TH1F* z_mass_cp;
			TH1F* z_met_pp;
			TH1F* z_mass_pp;
			TH1F* ele_stdVsPholike;
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
		void GetCommonValues();
		void GetCentralValues(std::vector<int> vPho, std::vector<int> vEle);
		void GetEmUpValues(std::vector<int> vPho, std::vector<int> vEle);
		void GetEmDownValues(std::vector<int> vPho, std::vector<int> vEle);
		void GetPhotonInfo(int iArrInd, SuperPhoton sp);
		void GetEmUpPhotonInfo(int iArrInd, SuperPhoton sp);
		void GetEmDownPhotonInfo(int iArrInd, SuperPhoton sp);
		void GetElectronInfo(int iArrInd, SuperPhoton sp);
		void GetEmUpElectronInfo(int iArrInd, SuperPhoton sp);
		void GetEmDownElectronInfo(int iArrInd, SuperPhoton sp);
		void GetTriggerInfo();
		void GetVertexInfo();
		void GetJetInfo();
		void GetRawJetInfo();
		void GetMetInfo();
		void FillDataBlocks(int ientry);
		void SetStupleName(std::string name) { sFileName = name; }
		std::string GetStupleName(std::string name) { return sFileName; }
		
		void TagForEmSys(std::vector<SuperPhoton>& vSp, int iSys);
		void GetEMSysPhotonInfo(std::vector<int> vPho);
		void GetEMSysElectronInfo(std::vector<int> vEle);
		void GetJesUpJetInfo();
		void GetJesDownJetInfo();
		void GetGenLevelInfo();
		void DumpSpAndJetBlock();
		void CheckLooseElectronTag();
		void GetHt() const;

		void MatchCalorTowers(int jet_ind, CalDataArray* cdholder);
		double  GetEmTiming(const int jetInd, const float jetDetEta);
		bool GetSummaryStat() const { return bNoSummary; }
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetStupleVersion() const { return fStupleVersion; }
		void SetStupleVersion(const float t) { fStupleVersion = t; }
		void SetEmEtThr(const float et) 
		{ 
			if (et <15.0)
			{
				std::cout << __FILE__ << ":" << __FUNCTION__ 
					<< red << "ERROR! Photon Et must be at least that was used to generate stuples (from version>7)\n"
					<< " This can cause jet list to have duplicate objects when adding additional EM objects\n"
					<< " back to jet list!\n"
					<< " Using default minimum value."
					<< clearatt << std::endl;
			} else 
			{	
				fMinEtThr = et; 
			}
		}
		float GetEmEtThr() const { return fMinEtThr;}
		

		
	private:
		Global_Counters_t counter;
		TFile *rootFile;
		TTree *tree;
		std::string sFileName;
		Stuple stuple;
		bool bRunPermit;
		bool bHasTrigInfo;   //if trigger info is available or not
		bool bMc;

		InitSuperPhotons* initSpMod;
		JetFilterModuleV3 *jetMod;
		JetFilterModuleV3 *jetModUp, *jetModDown;
		TriggerModule*  trigMod;
		TagBeamHalo* tagHaloMod;

		std::vector<SuperPhoton> spUp, spDown;  // these are to hold EM up/down SuperPhotons

		Hist_t Hist;

		bool bJetEmMismatch, bCosmicEvent;		// to save cosmic and jet/em mismatched events to debug
		TCesCpr cesCpr, cesCprUp, cesCprDown;
		bool bNoSummary;			//controls end job summary print
		float fStupleVersion;	//flat ntuple version
		float fMinEtThr;			//Et cut for all em objects
		jetMistagJun2009* matrix;
		
	ClassDef(FlatStupleMaker,1)
};

#endif
