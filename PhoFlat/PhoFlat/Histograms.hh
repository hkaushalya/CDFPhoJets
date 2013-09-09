#ifndef HISTOGRAMS_HH
#define HISTOGRAMS_HH

///////////////////////////////////////////////////////////
// Defines all the common histograms. Any moddule can    //
// use them easily and quickly, 11-15-2007 - sam         //
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>
/********************************************************************
 *	
 *	$Id: Histograms.hh,v 1.15 2011/05/26 17:51:39 samantha Exp $
 * $Log: Histograms.hh,v $
 * Revision 1.15  2011/05/26 17:51:39  samantha
 * PREVIOUS commit has wrong comments.
 * This version is same as previous commit.
 * ADDED: Hists ofr PtBal and PtBal against mass.
 *
 * Revision 1.14  2011/05/26 17:47:43  samantha
 * REMOVED: the MCPhotonJets:: from MCPhotonJets::PassPhoTightIdCuts().
 *
 * Revision 1.13  2011/04/26 23:45:35  samantha
 * ADDED: EmTime hist for Jet Histograms and DelEmTime for time separation
 * between photon and jet to photon-jet histograms.
 *
 * Revision 1.12  2010/08/25 15:29:34  samantha
 * ADDED: MET study hists for photon, jet and others. Profile histograms to see the
 * region of MET. Also Pt plots for parent objects.
 *
 * Revision 1.11  2010/06/21 19:44:33  samantha
 * ADDED: hists for phi separation of photon and leading jet with the met before
 * and after the delphi cut of 0.4.
 *
 * Revision 1.10  2010/04/13 16:10:11  samantha
 * ADDED: MetX, MetY hists to EventHists.
 *
 * Revision 1.9  2010/04/07 17:01:47  samantha
 * ADDED:1. delPhi(jet,MET) to the Jet hists collection.
 *       2. delPhi(pho,MET) to the Photon hists collection.
 *
 * Revision 1.8  2010/03/24 19:07:38  samantha
 * ADDED: VertexZ plot for highest sumPt vertex z distribution to the list of Event
 * hists.
 *
 * Revision 1.7  2010/02/09 23:56:25  samantha
 * ADDED a EventWeight Profile hist to the list of photon hist to see the weights used for sideband reweighing.
 *
 * Revision 1.6  2010/01/07 22:42:05  samantha
 * ADDED: pho-jet balance Pt and E plots to Photon1Jet hists
 *
 * Revision 1.5  2009/12/07 21:14:28  samantha
 * ADDED: 	1. CprWeight histogram to the photon hists.
 * 	2. Auto CVS Log into the file.
 *
 * 
 *******************************************************************/

#include <iostream>
#include "TFolder.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include <string>
#include <vector>
#include "TH1I.h"

class Histograms {

	public:
		Histograms();
		~Histograms();
		
	struct PhotonHists_t {	//all photon variables histos
			TH1F* Detector;
			TH1F* DetEta;
			TH1F* EvtEta;
			TH1F* DetPhi;
			TH1F* EvtPhi;
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
			TH1F* CprWeight;		//CPR weight distibution for the photon
			TProfile* EventWeights;		//weights used for the events
			TH1F* PhoMetDelPhi;
			TH1F* ClosestPhoMetDelPhi_b4;
			TH1F* ClosestPhoMetDelPhi_a4;
			TProfile* PhoEtMetProf;
			TProfile* PhoEMetProf;
			TProfile* PhoEtaMetProf;
		};


		struct JetHists_t { 		// all jet variables histos
			TH1F* DetEta;
			TH1F* EvtEta;
			TH1F* DetPhi;
			TH1F* EvtPhi;
			TH1F* EtCorr;
			TH1F* HadEm;
			TH1F* Emfr;
			TH1F* NTowers;
			TH1F* NTracks;
			//TH1F* SeedEtCorr;
			TH1F* SeedIPhi;
			TH1F* SeedIEta;
			TH1F* JetMetDelPhi;
			TH1F* ClosestJetMetDelPhi_b4;
			TH1F* ClosestJetMetDelPhi_a4;
			TProfile* JetEtMetProf;
			TProfile* JetEMetProf;
			TProfile* JetEtaMetProf;
			TH1F* EmTime;
		};

		struct Photon1JetHists_t {	// pho-jet histos
			TH1F* InvMass; 
			TH1F* DelPhi;
			TH1F* DelEta;
			TH1F* DelR;
			TH1F* EtRatio;
			TH2F* JetEmfrVsDelPhi;
			TH2F* JetEmfrVsPhoPhiWedge;
			TH1F* PhoJetPtBalance;	// JetEt/PhoEt - 1
			TH1F* PhoJetEBalance;	// JetE/PhoE - 1
			TH1F* Pt;		//Pt of the parent
			TH1F* DelEmTime;
		};
		
		struct Photon2JetsHists_t {	// pho-jets histos
			TH1F* InvMass; 
			TH1F* Pt;		//Pt of the parent
			TH1F* PtBal;	//Pt of the two jets to photon ratio
			TProfile* PtBalVsMass;
		};
	
		struct TwoJetsHists_t {	// two lead jets histos
			TH1F* InvMass; 
			TH1F* DelPhi;
			TH1F* DelEta;
			TH1F* DelR;
			TH1F* EtRatio;
			TH1F* Pt;		//Pt of the parent
		};
	
		struct EventHists_t {    // stuff common to the event
			TH1F* Met;
			TH1F* MetX;
			TH1F* MetY;
			TH1F* Met_Gen_d;
			TH1F* Sumet;
			TH1F* Ht;				//Ht without any weight
			TH1F* HtWgt;			// Ht with CPR weight
			TH1F *HtSys[8];		//for CES/CPR weight systematics
			TH1F* NVertices;		//for fVertexBlock->Vertices()
			TH1F* N12Vertices;	//for Class12Vertices()
			TH1F* VertexZ;			//highest sumpt vtx z cdt
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
			TH1F* Triggers;				//a bin for each trigger	0=pho25iso,1=pho50,2=pho70
			TProfile* N12VerticesMetProf;	// calss 12 vertices with MET
			TProfile* N12VerticesSumetProf;	// calss 12 vertices with SumEt
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

		void GetEventHistograms(EventHists_t&, std::string sText);
		void GetHaloHistograms(BeamHaloHists_t&, std::string sText);
		void GetCounterHistograms(CounterHists_t&, std::string sText,
											std::vector<std::string> vBinLabels);
		void GetPhoton1JetHistograms(Photon1JetHists_t&, std::string sText);
		void GetPhoton2JetsHistograms(Photon2JetsHists_t&, std::string sText);
		void GetTwoJetsHistograms(TwoJetsHists_t&, std::string sText);
		void GetJetHistograms(JetHists_t&, std::string sText);
		void GetPhotonHistograms(PhotonHists_t&, std::string text);

		
	private:

		char title[500];
		char ytitle[500];
		std::string sTextToAdd;

};

#endif
