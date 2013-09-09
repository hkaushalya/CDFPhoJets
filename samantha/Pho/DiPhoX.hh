#ifndef DIPHOTONX_HH
#define DIPHOTONX_HH

/*
 * $Log: DiPhoX.hh,v $
 * Revision 1.3  2009/06/26 19:54:17  samantha
 * ADDED: Option to manually set the MC flag. This way I can save some CPU time
 * by not running the halo/PMT tagging module when I run over MC samples.
 *
 * Revision 1.2  2009/06/16 04:30:54  samantha
 * MINOR CHANGES: 1: one more counter for triphoton events.
 *
 * Revision 1.1  2009/06/03 04:07:25  samantha
 * My own version of gg+X to check if my version of MetModel (JetFilterV2).
 * Require two photons (central and Et>13GeV). See elog:1206 for the first
 * result. For Sasha MetSig is not that good. But the difference is smaller
 * than for my pho+jets.
 *
 *
 */

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/Pho/TagPMTSpikes.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "samantha/Pho/TagBeamHalo.hh"
#include "samantha/Pho/TagPhoenixPhotons.hh"
#include <fstream>

class DiPhoX: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		DiPhoX(const char* name="DiPhoX", const char* title="DiPhoX");
		~DiPhoX();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


				/************** abstract objects *************************/

		// enumerators to count the photon events
			enum COUNTS {
				EVENTS_PROCESSED 		= 0,
				EVENTS_PASSED    		= 1,
				TIGHT_PHO_EVENTS   	= 2,				//inclusive all events with a tight photon
				ONE_TIGHT_PHO_ONLY   = 3,			
				TWO_TIGHT_PHO_ONLY   = 4,
				MORETHAN2_TIGHT_PHO  = 5
			};



		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
		void FillDataBlocks(int ientry);

		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }

		//HISTS
		struct Hist_t {
			TH1F *LeadPhoEt;
			TH1F *SubPhoEt;
			TH1F *LeadPhoEta;
			TH1F *SubPhoEta;
			TH1F *DiPhoInvMass;
			TH1F *DelPhi; 
		};
		

	private:
		bool 	qMc;  				// true if MonteCarlo data
		InitSuperPhotons* initSpMod;

		TagTightPhotons* tightMod;
		TagPMTSpikes* pmtMod;
		TagElectrons* eleMod;
		TagBeamHalo* haloMod;
		TagPhoenixPhotons* phoenixMod;

		bool bPermitRun;			// make sure all pre-req are met before run

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		int iHaloType;
		int iPhotonType; 		//type of the photon selection, 0=signal, 1=sideband
		
		float fPhoMaxDetEta;   //absolute upper limit on Det Eta of the photon
		float fPhoMinDetEta;   //absolute lower limit on Det Eta of the photon
		float fPhoMinEt;			//Min Et of the photon
		bool bRejPhxPhos;			// 1=Reject phoenix photons, 0=do not reject phx photons 
		
		bool bNoSummary;			//control the end job summary print

		Hist_t Hist;

		bool bManualMCflag;		// sets MC flag manually. this way I can avoid using PMT/BH
										// tag modules when running on MC samples

	public:
		int GetHaloType() const { return iHaloType; }
		void SetHaloType(const int iT); 
		void SetPhotonType(const int type);
		int GetPhotonType() const { return iPhotonType; }
		std::string GetPhotonTypeLabel() const
		{
			if (GetPhotonType() == 0) return "Signal Photon";
			else if (GetPhotonType() == 1) return "Sideband Photon";
			else 
			{
				return "??? YOU TELL ME! ???????";
			}
		};


		float GetPhoMaxDetEta() const { return fPhoMaxDetEta; }
		float GetPhoMinDetEta() const { return fPhoMinDetEta; }
		float GetPhoMinEt() const { return fPhoMinEt; }
		void SetPhoMaxDetEta(const float eta) { fPhoMaxDetEta = fabs(eta); }
		void SetPhoMinDetEta(const float eta) { fPhoMinDetEta = fabs(eta); }
		void SetPhoMinEt(const float et)
		{ 
			if (et < 0)
			{
				std::cout << "Photon Et must be > 0. Using default value." << std::endl;
			} else
			{
				fPhoMinEt = et;
			}
		}
	
		void SetRejPhxPhos(const bool rej) { bRejPhxPhos = rej; }
		bool GetRejPhxPhos() const { return bRejPhxPhos; }
		void FillPhotonHists(const int iPho1, const int iPho2);
 		void SetManualMCflag(const bool mc) { bManualMCflag = mc; }
 		bool ManualMCflag() const { return bManualMCflag; }
 	
		
	ClassDef(DiPhoX,1)
};

#endif
