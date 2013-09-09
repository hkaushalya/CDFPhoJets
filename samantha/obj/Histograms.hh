///////////////////////////////////////////////////////////
// Defines all the common histograms. Any moddule can    //
// use them easily and quickly, 11-15-2007 - sam         //
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>

#ifndef HISTOGRAMS_HH
#define HISTOGRAMS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>
#include "samantha/utils/FreeFunctions.hh"
#include "TFolder.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include <string>
#include <vector>
#include "TH1I.h"
#endif

class Histograms {

	protected:


	public:
		Histograms();
		~Histograms();
		
		struct EmObjHists_t {	//all photon variables histos
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
			TH1F* EmTime;
			TH1F* EoverP;
			TH1F* TrkZ0;
			TProfile* EmTimeVsRun;
			TProfile* EtCorrVsRun;
		};

	struct PhotonHists_t {	//all photon variables histos
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
			TH1F* EmTime;
			TH1F* EoverP;
			TH1F* TrkZ0;
			TProfile* EmTimeVsRun;
			TProfile* EtCorrVsRun;
			TH1F* PhiWedge;
			TH1F* PhoMetDelPhi;			// phi separation from MET
			TH1F* PhoIDs;
		};


		struct JetHists_t { 		// all jet variables histos
			TH1F* DetEta;
			TH1F* DetPhi;
			TH1F* EtCorr;
			TH1F* HadEm;
			TH1F* Emfr;
			TH1F* NTowers;
			TH1F* NTracks;
			//TH1F* SeedEtCorr;
			TH1F* EvtPhi;
			TH1F* SeedIPhi;
			TH1F* SeedIEta;
			TH1F* JetMetDelPhi;			// phi separation from MET
			
		};

		struct TwoObjHists_t {	// all di-jet ,pho-jet histos
			TH1F* InvMass; 
			TH1F* DelPhi;
			TH1F* DelEta;
			TH1F* DelR;
			TH1F* EtRatio;
			TProfile* InvMassRun; 
			TProfile* DelPhiRun;
			TProfile* DelEtaRun;
			TProfile* DelRRun; 
			TProfile* EtRatioRun;
		};
	
		struct Photon1JetHists_t {	// pho-jet histos
			TH1F* InvMass; 
			TH1F* DelPhi;
			TH1F* DelEta;
			TH1F* DelEtaDet;
			TH1F* DelEtaPlus;
			TH1F* DelEtaDetPlus;
			TH1F* DelR;
			TH1F* EtRatio;
			TH2F* JetEmfrVsDelPhi;
			TH2F* JetEmfrVsPhoPhiWedge;
		};
		
		struct Photon2JetsHists_t {	// pho-jets histos
			TH1F* InvMass; 
		};
	
		struct TwoJetsHists_t {	// two lead jets histos
			TH1F* InvMass; 
			TH1F* DelPhi;
			TH1F* DelEta;
			TH1F* DelDetEta;
			TH1F* DelR;
			TH1F* EtRatio;
		};
	
	
		struct ThreeObjHists_t {	// all photon +di-jet histos
			TH1F* InvMass;
			TProfile* InvMassRun;
		};
		
		struct EventHists_t {    // stuff common to the event
			TH1F* Met;
			TH1F* Sumet;
			TH1F* Ht;
			TH1F* NVertices;		//for fVertexBlock->Vertices()
			TH1F* N12Vertices;	//for Class12Vertices()
			TH1F* NJet15;
			TProfile* MetRun;				// met with run
			TProfile* SumetRun;			// sumet with run
			TProfile* NVerticesRun;		// vertices in vertex block with run
			TProfile* N12VerticesRun;	// calss 12 vertices with run
			TProfile* NJet15Run;			// Jet15 (after removing em objects) with run
			TH2F* NJet15VsMet;
			TH2F* NJet15VsSumet;
			TH2F* N12VerticesVsMet;
			TH2F* N12VerticesVsSumet;
			TH1F* MetSumetRatio;			//MET/SUMET
			TH1F* BestVertexZ;
		};
		
		//______________________________ Beam Halo study histograms
		struct BeamHaloHists_t {

			TH1F* fHaloSeedWedgeH;   // number of seed wedge towers with Had>thresh 
			TH1F* fHaloSideWedgeH;   // number of side wedge towers with Em>thresh
			TH1F* fHaloSeedWedge;   // number of seed wedge towers with Em>thresh 
			TH1F* fHaloSideWedge;   // number of side wedge towers with Em>thresh
			TH1F* fHaloEastNHad;    // number of had towers on East side
			TH1F* fHaloWestNHad;    // number of had towers on West side
			TH1F* fHaloEmSeedE;     // EM energy of seed wedge towers
			TH1F* fHaloEmSideE;     // EM energy of side wedge towers
			TH1F* fHaloHadSeedE;    // Plug HAD energy of seed wegde towers
			TH1F* fHaloSeedWedgeHadE; // HAD energy of seed wegde towers
			TH1F* fHaloNHad;        // number of towers on East & West: EastNHad+WestNHad

			TH2F* fHaloSeedWedgeEM_SeedWedgeHAD; // number of EM seed wedge towers vs. number of HAD towers 

			TH1F* fHaloEast;        // output of Pho->HaloEast()
			TH1F* fHaloWest;        // output of Pho->HaloWest()
			TH1F* fHaloEastWest;    // Pho->HaloEast()+Pho->HaloWest()
			TH1F* fHaloSideWedgeR;  // number of towers in side wedges according to Ray's code  
			TH1F* fHaloSeedWedgeR;  // number of towers in seed wedge according to Ray's code  
			TH1F* fHaloTwrInRawR;   // number of continues towers in seed wedge according to Ray's code  

			TH1F* fHaloCesStripEoverE; // CesStripE/E(pho)
			TH1F* fHaloCesWireEoverE;  // CesWireE/E(pho)
			TH1F* fHaloHadTDC;         // Had TDC timing   

			TH1F* fHaloPhiWedge;  // wedge number of halo candidate
			TH1F* fHaloEta;       // eta of halo candidate
			TH1F* fHaloIso;       // iso of halo candidate
			TH1F* fHaloHadEm;     // Had/Em of halo candidate
			TH1F* fHaloCesChi2;   // CES Chi2 of halo candidate
			TH1F* fHaloEt;        // Et of halo candidate

			TH2F*  fSeedWedge_HadTowers; //seed towers vs N Had towers
		};
		
		struct CounterHists_t {
			TH1I* iCounter;
		};

		void GetEmObjHistograms(TFolder* new_folder,EmObjHists_t&, std::string text);
		void GetEventHistograms(TFolder* new_folder, EventHists_t&, std::string sText);
		void GetTwoObjHistograms(TFolder* new_folder, TwoObjHists_t&, std::string sText);
		void GetThreeObjHistograms(TFolder* new_folder, ThreeObjHists_t&, std::string sText);
		void GetHaloHistograms(TFolder* new_folder, BeamHaloHists_t&, std::string sText);

		
		void GetCounterHistograms(TFolder* new_folder, CounterHists_t&, std::string sText,
											std::vector<std::string> vBinLabels);
		void GetPhoton1JetHistograms(TFolder* new_folder, Photon1JetHists_t&, std::string sText);
		void GetPhoton2JetsHistograms(TFolder* new_folder, Photon2JetsHists_t&, std::string sText);
		void GetTwoJetsHistograms(TFolder* new_folder, TwoJetsHists_t&, std::string sText);
		void GetJetHistograms(TFolder* new_folder, JetHists_t&, std::string sText);
		void GetPhotonHistograms(TFolder* new_folder,PhotonHists_t&, std::string text);

		
		void Print(Option_t *opt="a");		// print its info

		//accessors
		
		
	private:

		char title[500];
		char ytitle[500];
		std::string sTextToAdd;

};

#endif
