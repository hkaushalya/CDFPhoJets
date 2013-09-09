#ifndef PHO2JETSTEMP_HH
#define PHO2JETSTEMP_HH

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
#include "samantha/obj/Histograms.hh"
#include "samantha/Pho/EventProperties.hh"
#include "samantha/Pho/TagPMTSpikes.hh"
#include "samantha/Pho/TagElectrons.hh"
#include "samantha/Pho/TagBeamHalo.hh"
#include "samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/Pho/SetEmTimes.hh"
#include <fstream>
#include "Stntuple/obj/TGenpBlock.hh"
#include <Stntuple/obj/TStnJetBlock.hh>
#include "pdftools/TPDFReweight.hh"
#include "TTree.h"

class Pho2JetsTemp: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TGenpBlock* fGenpBlock;
		TStnJetBlock*      fJetBlockClu04;
		static const int iSignalPhoton = 0;
		static const int iSidebandPhoton = 1;


	public:
		Pho2JetsTemp(const char* name="PhoJetsTemp", const char* title="PhoJetsTemp");
		~Pho2JetsTemp();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TGenpBlock* GetGenpBlock() { return fGenpBlock; }


				/************** abstract objects *************************/

		// enumerators to count the photon events
			enum COUNTS {
				EVENTS_PROCESSED = 0,
				EVENTS_PASSED    = 1,
				TIGHT_PHO_EVENTS   = 2,				//inclusive all events with a tight photon
				ONE_TIGHT_PHO_ONLY = 3,
				TWO_TIGHT_PHO_ONLY = 4,
				ONE_TIGHT_ANY_LOOSE_PHO = 5,
			};


			const static UInt_t iNPDFs=41;			//maximum variation for any PDF set
			struct PDFSyst_t  //holds pdf weights
			{
				Int_t RunNumber;
				Int_t EventNumber;
				UInt_t iNpdfs;						//variation for the current selected PDF. this is how much pdf elements will be saved in the leaf
				Float_t fPdfWeight[iNPDFs];
				Float_t pj1_Et_pho;
				Float_t pj1_Et_j1;
				Float_t pj1_InvM_pj1;
				Float_t pj1_Met;
				Float_t pj1_Ht;
				Int_t pj1_Njet15;
				Float_t EtRatLeadJetLeadPho;
				Float_t DelRLeadPhoLeadJet;
				Float_t DelPhiLeadPhoLeadJet;
				Float_t pj2_Et_pho;
				Float_t pj2_Et_j1;
				Float_t pj2_Et_j2;
				Float_t pj2_InvM_pj1;
				Float_t pj2_InvM_pj2;
				Float_t pj2_InvM_pj1j2;
				Float_t pj2_InvM_j1j2;
				Float_t pj2_Met;
				Float_t pj2_Ht;
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


		bool DoJetSelection(int iSpIndex);

		void SetEmTimeWindow(float tmin_,float tmax_) {
				fEmTimeMin = tmin_;
				fEmTimeMax = tmax_;
		}
		void GetEmTimeWindow(float &tmin_,float &tmax_) {
			tmin_ = fEmTimeMin;
			tmax_ = fEmTimeMax;
		}
		float GetEmTimeWindowMin() const { return fEmTimeMin; }
		float GetEmTimeWindowMax() const { return fEmTimeMax; }

		void FillPhotonIDHist(const int iSpInd);
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }

	private:
		bool 	qMc;  				// true if MonteCarlo data
		InitSuperPhotons* initSpMod;
		TagTightPhotons* tightMod;
		EventProperties* evtPropMod;
		TagPMTSpikes* pmtMod;
		TagElectrons* eleMod;
		TagBeamHalo* haloMod;
		SetEmTimes* setEMTMod;
		TagPhoenixPhotons* phoenixMod;
		JetFilterModuleV2 *jetMod;

		bool bPermitRun;			// make sure all pre-req are met before run

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::PhotonHists_t h1_Pho, h2_Pho;
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;
		Histograms::Photon1JetHists_t h1_PhoJet, h2_PhoJet1, h2_PhoJet2;
		Histograms::Photon2JetsHists_t h2_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t h2_TwoJets;

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		float fEmTimeMin, fEmTimeMax;			// em time window
		int iHaloType;
		int iPhotonType; 		//type of the photon selection, 0=signal, 1=sideband
		
		float fPhoMaxDetEta;   //absolute upper limit on Det Eta of the photon
		float fPhoMinDetEta;   //absolute lower limit on Det Eta of the photon
		float fPhoMinEt;			//Min Et of the photon
		bool bRejPhxPhos;			// 1=Reject phoenix photons, 0=do not reject phx photons 
		
		SuperPhoton* thePhoton;  // this is THE photon passed all selections cuts
		bool bNoSummary;			//control the end job summary print

		TH1F *hLeadJet_Et;		// lead jet after excluding the lead photon from the jet list
		TH1F *hPhoLeadJetPtRatio;
		//TH1F* hPDF41hist[41][2];		//pdf systematics for photon Et and jet et
		std::vector<TH1F*> vPDFhist;		//pdf systematics for photon Et and jet et
		TH2F* hQvsE;
		TH1F* hQERat;
		TH1F* hQSumPtRat;
		TH1F* hQPhoPtRat;
		TH1F* hQJetPtRat;
		TPDFReweight *pdfWmod;			//PDF weights module
		std::vector<float> vPDFweights; //holds PDF syst uncertainties
		int iNaltenativePDFs;		//number of alternative PDFs (or weights_ depends on the MC GENERATOR
		PDFSyst_t aPDFsyst;
		TTree *treePDFw;			//holds pdf weight for events to generate any hist after the job

		//temp stuff
		/*TH1F* hPhoEt;
		TH1F* hDetEta;
		TH1F* hLoosePhoEt;
		TH1F* hLooseDetEta;
		TH1F* hTightPhoEt;
		TH1F* hHepgPhoEt;
		TH1F* hHepgPhoEta;
      */
		bool bDoPDFSyst;
		TH1F* metb4;
		TH1F* meta4;
		TH1F* phometdelphi_a4;
		TH1F* jetmetdelphi_a4;



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
		void SetPhoMaxDetEta(const float eta);
		void SetPhoMinDetEta(const float eta);
		void SetPhoMinEt(const float et);
	
		void SetRejPhxPhos(const bool rej) { bRejPhxPhos = rej; }
		bool GetRejPhxPhos() const { return bRejPhxPhos; }
		SuperPhoton* GetPhoton() const;
		bool FillPhoPDFSystHists(const int iPhoIndex);
		void SetDoPDFSyst(const bool b) { bDoPDFSyst = b; };
		bool DoPDFSyst() const { return bDoPDFSyst;}
		
	ClassDef(Pho2JetsTemp,1)
};

#endif
