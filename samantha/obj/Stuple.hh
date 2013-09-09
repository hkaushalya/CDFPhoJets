#ifndef STUPLE_HH
#define STUPLE_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#include <Rtypes.h>
#endif

///////////////////////////////////////////////////////////
// This is used skim the Stntuple and create the flat    //
// ntuple. 02-14-2008 -sam                               //
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>
/*
05-19-08:	added new jet collections: JES UP and DOWN.
				Partial photons with EM energy corrections
				up/down by 1%. For photons
06-08-08:   I created 3 complete set of varibles for
            central/ up/ down photons. Although 99% of the
				info is redundant, it is easiear to use.

*/

/*	$Id: Stuple.hh,v 1.15 2011/05/25 21:55:36 samantha Exp $
 * $Log: Stuple.hh,v $
 * Revision 1.15  2011/05/25 21:55:36  samantha
 * ADDED: SecVtx tag info to all jet collections.
 *
 * Revision 1.14  2011/04/03 22:40:35  samantha
 * ADDED: fTrkPt2 to get the Pt of the 2nd track pointing to the cluster.
 *
 * Revision 1.13  2010/04/13 16:31:25  samantha
 * ADDED: MetX,MetY, RawMet and raw jets info.
 *
 * Revision 1.12  2010/02/07 04:40:51  samantha
 * MODIFIED: ele_StdLooseEleId (also up/down versions) to ele_StdLooseId. because
 * the 'Ele' part is already in the prefix.
 *
 *
 */

#include "TLorentzVector.h"


class Stuple {

	protected:

	public:
			Stuple();
			void Init();		//______________________ reset values to some large negative values
			void Dump(int iObjType =1234); // __ dump all contents 1234=dump all 1=tight pho 2=tight ele 3= loose pho 4= loose ele
			void DumpAll();							//__ dump all
			void DumpPhotonBlock() const;
			void DumpPhoton(const int iIndex) const;
			void DumpElectronBlock() const;
			void DumpElectron(const int iIndex) const;
			void DumpJetBlock() const;
			void DumpJet(const int iIndex) const;
			void DumpMetBlock() const;
			void DumpVtxBlock() const;
			void DumpEvtBlock() const;
			void DumpGenBlock() const;
			void DumpEventSummary() const;
			
			static const UInt_t Npho = 10;
			static const UInt_t Nele = 10;
			static const UInt_t Njet = 30;
			static const UInt_t Nrawjet = 100;
			static const UInt_t Ngenpho = 100;		//generator level photons
			static const UInt_t Ngenele = 100;		//generator level electrons
			static const UInt_t Ngenjet = 100;		//generator level jets

			Float_t Version;								// Stuple version. I can keep track in the code by ref. to this

			Int_t evt_McFlag;		// =1 if MC , else 0
			Int_t evt_RunNumber;
			Int_t evt_EventNumber;
			Int_t tri_pho25iso;
			Int_t tri_pho50;
			Int_t tri_pho70;

			UInt_t  pho_num;							// number of electrons so we know how much to read from arrays
			UInt_t  pho_Ntight;						// # of tight photons
			UInt_t  pho_Nloose;						// # of loose photons
			Int_t   pho_Index[Npho];				//this is the index(order) in the Photon block
			Int_t   pho_PhoBlockIndex[Npho];  	//this is the index(order) in the Photon block
			Float_t pho_Etc[Npho];
			Float_t pho_E[Npho];
			Float_t pho_Px[Npho];
			Float_t pho_Py[Npho];
			Float_t pho_Pz[Npho];
			Int_t   pho_Detector[Npho];
			Float_t pho_DetEta[Npho];
			Float_t pho_DetPhi[Npho];
			Float_t pho_XCes[Npho];
			Float_t pho_ZCes[Npho];
			Float_t pho_HadEm[Npho];
			Float_t pho_Chi2Mean[Npho];
			Int_t   pho_N3d[Npho];
			Float_t pho_Iso4[Npho];
			Float_t pho_TrkPt[Npho];				// from photon block
			Float_t pho_TrkPt2[Npho];				// from photon block, 2nd highest pt track 
			Float_t pho_TrkIso[Npho];				// from photon block
			Float_t pho_CesWireE2[Npho];
			Float_t pho_CesStripE2[Npho];
			Int_t   pho_PhiWedge[Npho];
			Int_t   pho_NMuonStubs[Npho];			//trkless muon stubs around 30 degree cone of the photon
			Float_t pho_EmTime[Npho];
			Int_t   pho_TightId[Npho];
			Int_t   pho_LooseId[Npho];
			Int_t   pho_PhoenixId[Npho];
			Int_t   pho_Halo_seedWedge[Npho];
			Int_t   pho_Halo_eastNhad[Npho];
			Int_t   pho_Halo_westNhad[Npho];
			Int_t   pho_matchJetIndex[Npho]; 	//if a matching jet is found and removed from jet list
			Float_t pho_CprWgt[Npho];				//cpr weight
			Float_t pho_CprSys1[Npho];				//systematic errors for cpr weight (need 8 bins per photon)
			Float_t pho_CprSys2[Npho];
			Float_t pho_CprSys3[Npho];
			Float_t pho_CprSys4[Npho];
			Float_t pho_CprSys5[Npho];
			Float_t pho_CprSys6[Npho];
			Float_t pho_CprSys7[Npho];
			Float_t pho_CprSys8[Npho];
			


			// EM E uncertainty up 1%
			UInt_t  pho_up_num;							// number of electrons so we know how much to read from arrays
			UInt_t  pho_up_Ntight;						// # of tight photons
			UInt_t  pho_up_Nloose;						// # of loose photons
			Int_t   pho_up_Index[Npho];				//this is the index(order) in the Photon block
			Int_t   pho_up_PhoBlockIndex[Npho];  	//this is the index(order) in the Photon block
			Float_t pho_up_Etc[Npho];
			Float_t pho_up_E[Npho];
			Float_t pho_up_Px[Npho];
			Float_t pho_up_Py[Npho];
			Float_t pho_up_Pz[Npho];
			Int_t   pho_up_Detector[Npho];
			Float_t pho_up_DetEta[Npho];
			Float_t pho_up_DetPhi[Npho];
			Float_t pho_up_XCes[Npho];
			Float_t pho_up_ZCes[Npho];
			Float_t pho_up_HadEm[Npho];
			Float_t pho_up_Chi2Mean[Npho];
			Int_t   pho_up_N3d[Npho];
			Float_t pho_up_Iso4[Npho];
			Float_t pho_up_TrkPt[Npho];				// from photon block
			Float_t pho_up_TrkPt2[Npho];				// from photon block, 2nd highest pt track
			Float_t pho_up_TrkIso[Npho];				// from photon block
			Float_t pho_up_CesWireE2[Npho];
			Float_t pho_up_CesStripE2[Npho];
			Int_t   pho_up_PhiWedge[Npho];
			Int_t   pho_up_NMuonStubs[Npho];			//trkless muon stubs around 30 degree cone of the photon
			Float_t pho_up_EmTime[Npho];
			Int_t   pho_up_TightId[Npho];
			Int_t   pho_up_LooseId[Npho];
			Int_t   pho_up_PhoenixId[Npho];
			Int_t   pho_up_Halo_seedWedge[Npho];
			Int_t   pho_up_Halo_eastNhad[Npho];
			Int_t   pho_up_Halo_westNhad[Npho];
			Int_t   pho_up_matchJetIndex[Npho]; 	//if a matching jet is found and removed from jet list
			Float_t pho_up_CprWgt[Npho];
			Float_t pho_up_CprSys1[Npho];				//systematic errors for cpr weight (need 8 bins per photon)
			Float_t pho_up_CprSys2[Npho];
			Float_t pho_up_CprSys3[Npho];
			Float_t pho_up_CprSys4[Npho];
			Float_t pho_up_CprSys5[Npho];
			Float_t pho_up_CprSys6[Npho];
			Float_t pho_up_CprSys7[Npho];
			Float_t pho_up_CprSys8[Npho];
			

			// EM E uncertainty down 1%
			UInt_t  pho_down_num;							// number of electrons so we know how much to read from arrays
			UInt_t  pho_down_Ntight;						// # of tight photons
			UInt_t  pho_down_Nloose;						// # of loose photons
			Int_t   pho_down_Index[Npho];				//this is the index(order) in the Photon block
			Int_t   pho_down_PhoBlockIndex[Npho];  	//this is the index(order) in the Photon block
			Float_t pho_down_Etc[Npho];
			Float_t pho_down_E[Npho];
			Float_t pho_down_Px[Npho];
			Float_t pho_down_Py[Npho];
			Float_t pho_down_Pz[Npho];
			Int_t   pho_down_Detector[Npho];
			Float_t pho_down_DetEta[Npho];
			Float_t pho_down_DetPhi[Npho];
			Float_t pho_down_XCes[Npho];
			Float_t pho_down_ZCes[Npho];
			Float_t pho_down_HadEm[Npho];
			Float_t pho_down_Chi2Mean[Npho];
			Int_t   pho_down_N3d[Npho];
			Float_t pho_down_Iso4[Npho];
			Float_t pho_down_TrkPt[Npho];				// from photon block
			Float_t pho_down_TrkPt2[Npho];				// from photon block, 2nd highest pt track
			Float_t pho_down_TrkIso[Npho];				// from photon block
			Float_t pho_down_CesWireE2[Npho];
			Float_t pho_down_CesStripE2[Npho];
			Int_t   pho_down_PhiWedge[Npho];
			Int_t   pho_down_NMuonStubs[Npho];			//trkless muon stubs around 30 degree cone of the photon
			Float_t pho_down_EmTime[Npho];
			Int_t   pho_down_TightId[Npho];
			Int_t   pho_down_LooseId[Npho];
			Int_t   pho_down_PhoenixId[Npho];
			Int_t   pho_down_Halo_seedWedge[Npho];
			Int_t   pho_down_Halo_eastNhad[Npho];
			Int_t   pho_down_Halo_westNhad[Npho];
			Int_t   pho_down_matchJetIndex[Npho]; 	//if a matching jet is found and removed from jet list
			Float_t pho_down_CprWgt[Npho];
			Float_t pho_down_CprSys1[Npho];				//systematic errors for cpr weight (need 8 bins per photon)
			Float_t pho_down_CprSys2[Npho];
			Float_t pho_down_CprSys3[Npho];
			Float_t pho_down_CprSys4[Npho];
			Float_t pho_down_CprSys5[Npho];
			Float_t pho_down_CprSys6[Npho];
			Float_t pho_down_CprSys7[Npho];
			Float_t pho_down_CprSys8[Npho];
	

			//electrons collections
			UInt_t  ele_num;							// number of electrons so we know how much to read from arrays
			UInt_t  ele_Ntight;
			UInt_t  ele_Nloose;
			Int_t   ele_Index[Nele];				//this is the index(order) in the Photon block
			Int_t   ele_PhoBlockIndex[Nele];		//this is the index(order) in the Photon block
			Int_t   ele_EleBlockIndex[Nele];		//this is the index(order) in the Electron block
			Float_t ele_Etc[Nele];
			Float_t ele_E[Nele];
			Float_t ele_Px[Nele];
			Float_t ele_Py[Nele];
			Float_t ele_Pz[Nele];
			Int_t   ele_Detector[Nele];
			Float_t ele_DetEta[Nele];
			Float_t ele_DetPhi[Nele];
			Float_t ele_XCes[Nele];
			Float_t ele_ZCes[Nele];
			Float_t ele_HadEm[Nele];
			Float_t ele_Chi2Mean[Nele];
			Int_t   ele_N3d[Nele];
			Float_t ele_Iso4[Nele];
			Float_t ele_TrkIso[Nele];				//from photon block
			Float_t ele_CesWireE2[Nele];
			Float_t ele_CesStripE2[Nele];
			Int_t   ele_PhiWedge[Nele];
			Int_t   ele_NMuonStubs[Nele];			//trkless muon stubs around 30 degree cone of the photon
			Float_t ele_EmTime[Nele];
			Int_t   ele_PhoenixId[Nele];
			Int_t   ele_Halo_seedWedge[Nele];
			Int_t   ele_Halo_eastNhad[Nele];
			Int_t   ele_Halo_westNhad[Nele];
			Int_t   ele_matchJetIndex[Nele];		//if a matching jet is found and removed from jet list
			Int_t   ele_Ntracks[Nele];				//these are from TStnElectron
			Float_t ele_Emfr[Nele];
			Float_t ele_EoverP[Nele];
			Float_t ele_TrackPt[Nele];
			Float_t ele_TrackPt2[Nele];   //2nd highest pt track
			Float_t ele_TrackBcPt[Nele];
			Float_t ele_TrackPhi[Nele];
			Int_t   ele_Nssl[Nele];
			Int_t   ele_Nasl[Nele];
			Int_t   ele_TightId[Nele];				// this is using photon-like electron id cuts
			Int_t   ele_LooseId[Nele];				// this is using photon-like electron id cuts with a loose e/p cut
			Int_t   ele_ConversionId[Nele];
			Int_t   ele_StdLooseId[Nele];		// this is using the Standard electron id
			Int_t   ele_StdTightId[Nele];		// this is using the Standard electron id


			// EM E uncertainty up 1%
			UInt_t  ele_up_num;							// number of electrons so we know how much to read from arrays
			UInt_t  ele_up_Ntight;
			UInt_t  ele_up_Nloose;
			Int_t   ele_up_Index[Nele];				//this is the index(order) in the Photon block
			Int_t   ele_up_PhoBlockIndex[Nele];		//this is the index(order) in the Photon block
			Int_t   ele_up_EleBlockIndex[Nele];		//this is the index(order) in the Electron block
			Float_t ele_up_Etc[Nele];
			Float_t ele_up_E[Nele];
			Float_t ele_up_Px[Nele];
			Float_t ele_up_Py[Nele];
			Float_t ele_up_Pz[Nele];
			Int_t   ele_up_Detector[Nele];
			Float_t ele_up_DetEta[Nele];
			Float_t ele_up_DetPhi[Nele];
			Float_t ele_up_XCes[Nele];
			Float_t ele_up_ZCes[Nele];
			Float_t ele_up_HadEm[Nele];
			Float_t ele_up_Chi2Mean[Nele];
			Int_t   ele_up_N3d[Nele];
			Float_t ele_up_Iso4[Nele];
			Float_t ele_up_TrkIso[Nele];				//from photon block
			Float_t ele_up_CesWireE2[Nele];
			Float_t ele_up_CesStripE2[Nele];
			Int_t   ele_up_PhiWedge[Nele];
			Int_t   ele_up_NMuonStubs[Nele];			//trkless muon stubs around 30 degree cone of the photon
			Float_t ele_up_EmTime[Nele];
			Int_t   ele_up_PhoenixId[Nele];
			Int_t   ele_up_Halo_seedWedge[Nele];
			Int_t   ele_up_Halo_eastNhad[Nele];
			Int_t   ele_up_Halo_westNhad[Nele];
			Int_t   ele_up_matchJetIndex[Nele];		//if a matching jet is found and removed from jet list
			Int_t   ele_up_Ntracks[Nele];				//these are from TStnElectron
			Float_t ele_up_Emfr[Nele];
			Float_t ele_up_EoverP[Nele];
			Float_t ele_up_TrackPt[Nele];
			Float_t ele_up_TrackPt2[Nele];  //2nd highest pt track
			Float_t ele_up_TrackBcPt[Nele];
			Float_t ele_up_TrackPhi[Nele];
			Int_t   ele_up_Nssl[Nele];
			Int_t   ele_up_Nasl[Nele];
			Int_t   ele_up_TightId[Nele];
			Int_t   ele_up_LooseId[Nele];
			Int_t   ele_up_ConversionId[Nele];
			Int_t   ele_up_StdLooseId[Nele];
			Int_t   ele_up_StdTightId[Nele];

			// EM E uncertainty down 1%
			UInt_t  ele_down_num;							// number of electrons so we know how much to read from arrays
			UInt_t  ele_down_Ntight;
			UInt_t  ele_down_Nloose;
			Int_t   ele_down_Index[Nele];				//this is the index(order) in the Photon block
			Int_t   ele_down_PhoBlockIndex[Nele];		//this is the index(order) in the Photon block
			Int_t   ele_down_EleBlockIndex[Nele];		//this is the index(order) in the Electron block
			Float_t ele_down_Etc[Nele];
			Float_t ele_down_E[Nele];
			Float_t ele_down_Px[Nele];
			Float_t ele_down_Py[Nele];
			Float_t ele_down_Pz[Nele];
			Int_t   ele_down_Detector[Nele];
			Float_t ele_down_DetEta[Nele];
			Float_t ele_down_DetPhi[Nele];
			Float_t ele_down_XCes[Nele];
			Float_t ele_down_ZCes[Nele];
			Float_t ele_down_HadEm[Nele];
			Float_t ele_down_Chi2Mean[Nele];
			Int_t   ele_down_N3d[Nele];
			Float_t ele_down_Iso4[Nele];
			Float_t ele_down_TrkIso[Nele];				//from photon block
			Float_t ele_down_CesWireE2[Nele];
			Float_t ele_down_CesStripE2[Nele];
			Int_t   ele_down_PhiWedge[Nele];
			Int_t   ele_down_NMuonStubs[Nele];			//trkless muon stubs around 30 degree cone of the photon
			Float_t ele_down_EmTime[Nele];
			Int_t   ele_down_PhoenixId[Nele];
			Int_t   ele_down_Halo_seedWedge[Nele];
			Int_t   ele_down_Halo_eastNhad[Nele];
			Int_t   ele_down_Halo_westNhad[Nele];
			Int_t   ele_down_matchJetIndex[Nele];		//if a matching jet is found and removed from jet list
			Int_t   ele_down_Ntracks[Nele];				//these are from TStnElectron
			Float_t ele_down_Emfr[Nele];
			Float_t ele_down_EoverP[Nele];
			Float_t ele_down_TrackPt[Nele];
			Float_t ele_down_TrackPt2[Nele];    //2nd highest pt track
			Float_t ele_down_TrackBcPt[Nele];
			Float_t ele_down_TrackPhi[Nele];
			Int_t   ele_down_Nssl[Nele];
			Int_t   ele_down_Nasl[Nele];
			Int_t   ele_down_TightId[Nele];
			Int_t   ele_down_LooseId[Nele];
			Int_t   ele_down_ConversionId[Nele];
			Int_t   ele_down_StdLooseId[Nele];
			Int_t   ele_down_StdTightId[Nele];



			UInt_t  jet_num;							// number of jets so we know how much to read from arrays
			UInt_t  jet_NJet15;
			Int_t   jet_Index[Njet];				//original index in jet block
			Float_t jet_Pt[Njet];					//easy to get from JetFilter
			Float_t jet_E[Njet];
			Float_t jet_Px[Njet];
			Float_t jet_Py[Njet];
			Float_t jet_Pz[Njet];
			Float_t jet_DetEta[Njet];
			Float_t jet_DetPhi[Njet];
			Float_t jet_HadEm[Njet];
			Float_t jet_Emfr[Njet];
			Int_t   jet_Ntowers[Njet];
			Int_t   jet_Ntracks[Njet];
			Int_t   jet_SeedIPhi[Njet];
			Int_t   jet_SeedIEta[Njet];
			Float_t jet_EmTime[Njet];
			Int_t   jet_SecVtxTag[Njet];  // if it has a tag
			Float_t jet_SecVtxppb[Njet];	// positive tag prob
			Float_t jet_SecVtxnpb[Njet];  // negative tag prob
			Float_t jet_SecVtxTrkmass[Njet];
			
			//raw jet info
			UInt_t  jet_raw_num;
			Int_t   jet_raw_Index[Nrawjet];
			Float_t jet_raw_Pt[Nrawjet];
			Float_t jet_raw_E[Nrawjet];
			Float_t jet_raw_Px[Nrawjet];
			Float_t jet_raw_Py[Nrawjet];
			Float_t jet_raw_Pz[Nrawjet];

			
			// JES UP JET COLLECTION
			//some of this info overlaps with other jet collections. but it is ok. this way it is much clear
			// and I can stick to one jet collection to get all info.
			UInt_t  jet_up_num;							// number of jets so we know how much to read from arrays
			UInt_t  jet_up_NJet15;
			Int_t   jet_up_Index[Njet];				//original index in jet block
			Float_t jet_up_Pt[Njet];					//easy to get from JetFilter
			Float_t jet_up_E[Njet];
			Float_t jet_up_Px[Njet];
			Float_t jet_up_Py[Njet];
			Float_t jet_up_Pz[Njet];
			Float_t jet_up_DetEta[Njet];
			Float_t jet_up_DetPhi[Njet];
			Float_t jet_up_HadEm[Njet];
			Float_t jet_up_Emfr[Njet];
			Int_t   jet_up_Ntowers[Njet];
			Int_t   jet_up_Ntracks[Njet];
			Int_t   jet_up_SeedIPhi[Njet];
			Int_t   jet_up_SeedIEta[Njet];
			Float_t jet_up_EmTime[Njet];

			Int_t   jet_up_SecVtxTag[Njet];  // if it has a tag
			Float_t jet_up_SecVtxppb[Njet];	// positive tag prob
			Float_t jet_up_SecVtxnpb[Njet];  // negative tag prob
			Float_t jet_up_SecVtxTrkmass[Njet];
	
			// JES DOWN JET COLLECTION
			UInt_t  jet_down_num;							// number of jets so we know how much to read from arrays
			UInt_t  jet_down_NJet15;
			Int_t   jet_down_Index[Njet];				//original index in jet block
			Float_t jet_down_Pt[Njet];					//easy to get from JetFilter
			Float_t jet_down_E[Njet];
			Float_t jet_down_Px[Njet];
			Float_t jet_down_Py[Njet];
			Float_t jet_down_Pz[Njet];
			Float_t jet_down_DetEta[Njet];
			Float_t jet_down_DetPhi[Njet];
			Float_t jet_down_HadEm[Njet];
			Float_t jet_down_Emfr[Njet];
			Int_t   jet_down_Ntowers[Njet];
			Int_t   jet_down_Ntracks[Njet];
			Int_t   jet_down_SeedIPhi[Njet];
			Int_t   jet_down_SeedIEta[Njet];
			Float_t jet_down_EmTime[Njet];
			Int_t   jet_down_SecVtxTag[Njet];  // if it has a tag
			Float_t jet_down_SecVtxppb[Njet];	// positive tag prob
			Float_t jet_down_SecVtxnpb[Njet];  // negative tag prob
			Float_t jet_down_SecVtxTrkmass[Njet];
	
			
			Int_t   vtx_N;
			Int_t   vtx_NClass12;
			Float_t vtx_z;				// highest sum pt vertex Z
			Int_t   vtx_Ntracks;		// in highest sum pt vertex
			Float_t vtx_SumPt;			// highest

			Float_t met_RawMet;			//uncorrected met
			Float_t met_Met;				// corrected MET
			Float_t met_MetX;				// corrected MET X
			Float_t met_MetY;				// corrected MET Y
			Float_t met_SumEt;			// corrected
			Float_t met_Ht;				// corrected
			Float_t met_MetPhi;			// corrected
			//Float_t met_Gen_d;			//defualt generated met 
			//Float_t met_Gen_m;			//defualt generated met  - dev
			//Float_t met_Gen_p;			//defualt generated met  + dev
			//Float_t met_Gen_mUn;			//defualt generated met  - devUn
			//Float_t met_Gen_pUn;			//defualt generated met  + devUn



			//these are for MC generator level info
			Int_t gen_elenum;     // number of electrons at gen level (for Ws = 1, for Zs =2)
			Int_t gen_phonum;

			// i assume event will have only one mother and rest are its daughter
			Int_t   gen_MomIndex;		// order in the gen level
			Int_t   gen_MomPDG;      //PDG code
			Int_t   gen_MomStatus;  //stable, unstable ..
			Float_t gen_MomEtc;
			Float_t gen_MomE;
			Float_t gen_MomPx;
			Float_t gen_MomPy;
			Float_t gen_MomPz;
			Float_t gen_ProdVtxX;			//production vtx
			Float_t gen_ProdVtxY;
			Float_t gen_ProdVtxZ;
			Float_t gen_ProdVtxT;
			
			
			Int_t   gen_pho_Index[Ngenpho];		//order in the gen level;
			Int_t   gen_pho_PDG[Ngenpho];
			Int_t   gen_pho_Status[Ngenpho];
			Float_t gen_pho_Etc[Ngenpho];
			Float_t gen_pho_E[Ngenpho];
			Float_t gen_pho_Px[Ngenpho];
			Float_t gen_pho_Py[Ngenpho];
			Float_t gen_pho_Pz[Ngenpho];
			Float_t gen_pho_ProdVtxX[Ngenpho];			//production vtx
			Float_t gen_pho_ProdVtxY[Ngenpho];
			Float_t gen_pho_ProdVtxZ[Ngenpho];
			Float_t gen_pho_ProdVtxT[Ngenpho];

			Int_t   gen_ele_Index[Ngenele];		//order in the gen level;
			Int_t   gen_ele_PDG[Ngenele];
			Int_t   gen_ele_Status[Ngenele];
			Float_t gen_ele_Etc[Ngenele];
			Float_t gen_ele_E[Ngenele];
			Float_t gen_ele_Px[Ngenele];
			Float_t gen_ele_Py[Ngenele];
			Float_t gen_ele_Pz[Ngenele];
			Float_t gen_ele_ProdVtxX[Ngenele];			//production vtx
			Float_t gen_ele_ProdVtxY[Ngenele];
			Float_t gen_ele_ProdVtxZ[Ngenele];
			Float_t gen_ele_ProdVtxT[Ngenele];

			//level 3 summary Block
			//Int_t l3_em_num;
			//Float_t l3_em_Elet;
			//Float_t l3_em_Phet;
			//Float_t l3_em_Eta;
			//Float_t l3_em_Pt;
			//Float_t l3_em_;
			//Float_t l3_em_;
			//Float_t l3_em_;

			void SetVersion(const float v) { Version = v; }
			float GetVersion() const{ return Version; }

	private:

};

#endif
