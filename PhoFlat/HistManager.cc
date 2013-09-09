///////////////////////////////////////////////////////////
// This will hold common event properties like MEt, etc. //
// This will have function/s to fill a given Event Hists //
// derived from Histograms.cc/hh.                        //
// This talks to the JetFilterModule for correct MEt and //
// SumEt etc.                                            //
// All other mods can access these info and can get their//
// own histograms filled by functions in here.           //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>
/*********************************************************
 *
 * $Id: HistManager.cc,v 1.6 2010/04/07 17:34:57 samantha Exp $
 *
 * $Log: HistManager.cc,v $
 * Revision 1.6  2010/04/07 17:34:57  samantha
 * UPDATED: Electron and Eel+jets hist fill method to choose between different
 * electron/jet collections (central, em/jes up/down)
 *
 * Revision 1.5  2010/02/08 18:07:23  samantha
 * MINOR CHANGES: Commented out two print statements.
 *
 *
 *********************************************************/


#include "PhoFlat/HistManager.hh"
#include "TMath.h"
#include "TLorentzVector.h"





//___________________________________________________________________
HistManager::HistManager()
{
}

//___________________________________________________________________
HistManager::~HistManager()
{
}


//___________________________________________________________________
void HistManager::FillEventHists(const Stuple& stuple, 
						Histograms::EventHists_t& hist,
						const float fWeight)
{ 
// public method to be used by other mods to fill event properties histos
	hist.Met->Fill(stuple.met_Met, fWeight);
	//hist.Met_Gen_d->Fill(stuple.met_Gen_d, fWeight);
	hist.Sumet->Fill(stuple.met_SumEt, fWeight);
	hist.Ht->Fill(stuple.met_Ht, fWeight);
	hist.NVertices->Fill(stuple.vtx_N, fWeight);
	hist.N12Vertices->Fill(stuple.vtx_NClass12, fWeight);
	hist.VertexZ->Fill(stuple.vtx_z, fWeight);
	hist.NJet15->Fill(stuple.jet_NJet15, fWeight);
	hist.MetRun->Fill(stuple.evt_RunNumber, stuple.met_Met, fWeight);
	hist.SumetRun->Fill(stuple.evt_RunNumber, stuple.met_SumEt, fWeight);
	hist.NVerticesRun->Fill(stuple.evt_RunNumber, stuple.vtx_N, fWeight);
	hist.N12VerticesRun->Fill(stuple.evt_RunNumber, stuple.vtx_NClass12, fWeight);
	hist.NJet15Run->Fill(stuple.evt_RunNumber, stuple.jet_NJet15, fWeight);
	hist.NJet15VsMet->Fill(stuple.met_Met, stuple.jet_NJet15, fWeight);
	hist.NJet15VsSumet->Fill(stuple.met_SumEt, stuple.jet_NJet15, fWeight);
	hist.N12VerticesVsMet->Fill(stuple.met_Met, stuple.vtx_NClass12, fWeight);
	hist.N12VerticesVsSumet->Fill(stuple.met_SumEt, stuple.vtx_NClass12, fWeight);
}


//___________________________________________________________________
void HistManager::FillPhotonHists(const Stuple& stuple, const int collection,
						Histograms::PhotonHists_t& hist, 
						const int iPhoIndex, const float fWeight)
{

	int	Detector 	= 0;
	float DetEta 		= 0;
	float DetPhi 		= 0;
	float EtCorr 		= 0;
	float XCes 			= 0;
	float ZCes 			= 0;
	float HadEm 		= 0;
	float Chi2Mean 	= 0;
	int   N3d 			= 0;
	float Iso4 			= 0;
	float TrkPt 		= 0;
	float TrkIso 		= 0;
	float Ces2Wire 	= 0;
	float Ces2Strip 	= 0;
	float EmTime 		= 0;
	float PhiWedge 	= 0;
	float Run			= stuple.evt_RunNumber;

	switch (collection)
	{
		case 0:
		{
			Detector 	= stuple.pho_Detector[iPhoIndex];
			DetEta 		= stuple.pho_DetEta[iPhoIndex];
			DetPhi 		= stuple.pho_DetPhi[iPhoIndex];
			EtCorr 		= stuple.pho_Etc[iPhoIndex];
			XCes 			= stuple.pho_XCes[iPhoIndex];
			ZCes 			= stuple.pho_ZCes[iPhoIndex];
			HadEm 		= stuple.pho_HadEm[iPhoIndex];
			Chi2Mean 	= stuple.pho_Chi2Mean[iPhoIndex];
			N3d 			= stuple.pho_N3d[iPhoIndex];
			Iso4 			= stuple.pho_Iso4[iPhoIndex];
			TrkPt 		= stuple.pho_TrkPt[iPhoIndex];
			TrkIso 		= stuple.pho_TrkIso[iPhoIndex];
			Ces2Wire 	= stuple.pho_CesWireE2[iPhoIndex];
			Ces2Strip 	= stuple.pho_CesStripE2[iPhoIndex];
			EmTime 		= stuple.pho_EmTime[iPhoIndex];
			PhiWedge 	= stuple.pho_PhiWedge[iPhoIndex];
			break;
		}
		case 1:
		{
			Detector 	= stuple.pho_up_Detector[iPhoIndex];
			DetEta 		= stuple.pho_up_DetEta[iPhoIndex];
			DetPhi 		= stuple.pho_up_DetPhi[iPhoIndex];
			EtCorr 		= stuple.pho_up_Etc[iPhoIndex];
			XCes 			= stuple.pho_up_XCes[iPhoIndex];
			ZCes 			= stuple.pho_up_ZCes[iPhoIndex];
			HadEm 		= stuple.pho_up_HadEm[iPhoIndex];
			Chi2Mean 	= stuple.pho_up_Chi2Mean[iPhoIndex];
			N3d 			= stuple.pho_up_N3d[iPhoIndex];
			Iso4 			= stuple.pho_up_Iso4[iPhoIndex];
			TrkPt 		= stuple.pho_up_TrkPt[iPhoIndex];
			TrkIso 		= stuple.pho_up_TrkIso[iPhoIndex];
			Ces2Wire 	= stuple.pho_up_CesWireE2[iPhoIndex];
			Ces2Strip 	= stuple.pho_up_CesStripE2[iPhoIndex];
			EmTime 		= stuple.pho_up_EmTime[iPhoIndex];
			PhiWedge 	= stuple.pho_up_PhiWedge[iPhoIndex];
			break;
		}

		case -1:
		{
			Detector 	= stuple.pho_down_Detector[iPhoIndex];
			DetEta 		= stuple.pho_down_DetEta[iPhoIndex];
			DetPhi 		= stuple.pho_down_DetPhi[iPhoIndex];
			EtCorr 		= stuple.pho_down_Etc[iPhoIndex];
			XCes 			= stuple.pho_down_XCes[iPhoIndex];
			ZCes 			= stuple.pho_down_ZCes[iPhoIndex];
			HadEm 		= stuple.pho_down_HadEm[iPhoIndex];
			Chi2Mean 	= stuple.pho_down_Chi2Mean[iPhoIndex];
			N3d 			= stuple.pho_down_N3d[iPhoIndex];
			Iso4 			= stuple.pho_down_Iso4[iPhoIndex];
			TrkPt 		= stuple.pho_down_TrkPt[iPhoIndex];
			TrkIso 		= stuple.pho_down_TrkIso[iPhoIndex];
			Ces2Wire 	= stuple.pho_down_CesWireE2[iPhoIndex];
			Ces2Strip 	= stuple.pho_down_CesStripE2[iPhoIndex];
			EmTime 		= stuple.pho_down_EmTime[iPhoIndex];
			PhiWedge 	= stuple.pho_down_PhiWedge[iPhoIndex];
			break;
		}
	
	}


	hist.Detector->Fill(Detector, fWeight);
	hist.DetEta->Fill(DetEta, fWeight);
	hist.DetPhi->Fill(DetPhi, fWeight);
	hist.EtCorr->Fill(EtCorr, fWeight);
	hist.XCes->Fill(XCes, fWeight);
	hist.ZCes->Fill(ZCes, fWeight);
	hist.HadEm->Fill(HadEm, fWeight);
	hist.Chi2Mean->Fill(Chi2Mean, fWeight);
	hist.N3d->Fill(N3d, fWeight);
	hist.Iso4->Fill(Iso4, fWeight);
	hist.TrkPt->Fill(TrkPt, fWeight);
	hist.TrkIso->Fill(TrkIso, fWeight);
	hist.Ces2Wire->Fill(Ces2Wire, fWeight);
	hist.Ces2Strip->Fill(Ces2Strip, fWeight);
	hist.EmTime->Fill(EmTime, fWeight);
	hist.PhiWedge->Fill(PhiWedge, fWeight);
	hist.EmTimeVsRun->Fill(Run, EmTime, fWeight);
	hist.EtCorrVsRun->Fill(Run, EtCorr, fWeight);
}
//___________________________________________________________________
void HistManager::FillElectronHists(const Stuple& stuple, const int collection,
						Histograms::PhotonHists_t& hist, 
						const int iEleIndex, const float fWeight)
{
	int	Detector 	= 0;
	float DetEta 		= 0;
	float DetPhi 		= 0;
	float EtCorr 		= 0;
	float XCes 			= 0;
	float ZCes 			= 0;
	float HadEm 		= 0;
	float Chi2Mean 	= 0;
	int   N3d 			= 0;
	float Iso4 			= 0;
	float EoverP      = 0;
	float TrkPt 		= 0;
	float TrkIso 		= 0;
	float Ces2Wire 	= 0;
	float Ces2Strip 	= 0;
	float EmTime 		= 0;
	float PhiWedge 	= 0;

	switch (collection)
	{
		case 0:
		{
			Detector 	= stuple.ele_Detector[iEleIndex];
			DetEta 		= stuple.ele_DetEta[iEleIndex];
			DetPhi 		= stuple.ele_DetPhi[iEleIndex];
			EtCorr 		= stuple.ele_Etc[iEleIndex];
			XCes 			= stuple.ele_XCes[iEleIndex];
			ZCes 			= stuple.ele_ZCes[iEleIndex];
			HadEm 		= stuple.ele_HadEm[iEleIndex];
			Chi2Mean 	= stuple.ele_Chi2Mean[iEleIndex];
			N3d 			= stuple.ele_N3d[iEleIndex];
			Iso4 			= stuple.ele_Iso4[iEleIndex];
			EoverP      = stuple.ele_EoverP[iEleIndex];
			TrkPt 		= stuple.ele_TrackPt[iEleIndex];
			TrkIso 		= stuple.ele_TrkIso[iEleIndex];
			Ces2Wire 	= stuple.ele_CesWireE2[iEleIndex];
			Ces2Strip 	= stuple.ele_CesStripE2[iEleIndex];
			EmTime 		= stuple.ele_EmTime[iEleIndex];
			PhiWedge 	= stuple.ele_PhiWedge[iEleIndex];
			break;
		}
		case 1:
		{
			Detector 	= stuple.ele_up_Detector[iEleIndex];
			DetEta 		= stuple.ele_up_DetEta[iEleIndex];
			DetPhi 		= stuple.ele_up_DetPhi[iEleIndex];
			EtCorr 		= stuple.ele_up_Etc[iEleIndex];
			XCes 			= stuple.ele_up_XCes[iEleIndex];
			ZCes 			= stuple.ele_up_ZCes[iEleIndex];
			HadEm 		= stuple.ele_up_HadEm[iEleIndex];
			Chi2Mean 	= stuple.ele_up_Chi2Mean[iEleIndex];
			N3d 			= stuple.ele_up_N3d[iEleIndex];
			Iso4 			= stuple.ele_up_Iso4[iEleIndex];
			EoverP      = stuple.ele_up_EoverP[iEleIndex];
			TrkPt 		= stuple.ele_up_TrackPt[iEleIndex];
			TrkIso 		= stuple.ele_up_TrkIso[iEleIndex];
			Ces2Wire 	= stuple.ele_up_CesWireE2[iEleIndex];
			Ces2Strip 	= stuple.ele_up_CesStripE2[iEleIndex];
			EmTime 		= stuple.ele_up_EmTime[iEleIndex];
			PhiWedge 	= stuple.ele_up_PhiWedge[iEleIndex];
			break;
		}
		case -1:
		{
			Detector 	= stuple.ele_down_Detector[iEleIndex];
			DetEta 		= stuple.ele_down_DetEta[iEleIndex];
			DetPhi 		= stuple.ele_down_DetPhi[iEleIndex];
			EtCorr 		= stuple.ele_down_Etc[iEleIndex];
			XCes 			= stuple.ele_down_XCes[iEleIndex];
			ZCes 			= stuple.ele_down_ZCes[iEleIndex];
			HadEm 		= stuple.ele_down_HadEm[iEleIndex];
			Chi2Mean 	= stuple.ele_down_Chi2Mean[iEleIndex];
			N3d 			= stuple.ele_down_N3d[iEleIndex];
			Iso4 			= stuple.ele_down_Iso4[iEleIndex];
			EoverP      = stuple.ele_down_EoverP[iEleIndex];
			TrkPt 		= stuple.ele_down_TrackPt[iEleIndex];
			TrkIso 		= stuple.ele_down_TrkIso[iEleIndex];
			Ces2Wire 	= stuple.ele_down_CesWireE2[iEleIndex];
			Ces2Strip 	= stuple.ele_down_CesStripE2[iEleIndex];
			EmTime 		= stuple.ele_down_EmTime[iEleIndex];
			PhiWedge 	= stuple.ele_down_PhiWedge[iEleIndex];
			break;
		}
	}

	
	hist.Detector->Fill(Detector, fWeight);
	hist.DetEta->Fill(DetEta, fWeight);
	hist.DetPhi->Fill(DetPhi, fWeight);
	hist.EtCorr->Fill(EtCorr, fWeight);
	hist.XCes->Fill(XCes, fWeight);
	hist.ZCes->Fill(ZCes, fWeight);
	hist.HadEm->Fill(HadEm, fWeight);
	hist.Chi2Mean->Fill(Chi2Mean, fWeight);
	hist.N3d->Fill(N3d, fWeight);
	hist.Iso4->Fill(Iso4, fWeight);
	hist.EoverP->Fill(EoverP, fWeight);
	hist.TrkPt->Fill(TrkPt, fWeight);
	hist.TrkIso->Fill(TrkIso, fWeight);
	hist.Ces2Wire->Fill(Ces2Wire, fWeight);
	hist.Ces2Strip->Fill(Ces2Strip, fWeight);
	hist.EmTime->Fill(EmTime, fWeight);
	hist.PhiWedge->Fill(PhiWedge, fWeight);
	hist.EmTimeVsRun->Fill(stuple.evt_RunNumber, EmTime, fWeight);
	hist.EtCorrVsRun->Fill(stuple.evt_RunNumber, EtCorr, fWeight);
}


//___________________________________________________________________
void HistManager::FillJetHists(const Stuple& stuple, const int collection,
						Histograms::JetHists_t& hist,
						const int iJetIndex, const float fWeight)
{
	float Emfr  	= 0;
	float DetEta  	= 0;
	float Ntracks  = 0;
	float Ntowers  = 0;
	float EtCorr  	= 0;
	float EvtPhi  	= 0;
	float EvtEta  	= 0;
	int SeedIEta  	= 0;
	int SeedIPhi  	= 0;
	TLorentzVector tlVec(0,0,0,0);


	switch (collection)
	{
		case 0:
		{
			Emfr 		= stuple.jet_Emfr[iJetIndex];
			DetEta 	= stuple.jet_DetEta[iJetIndex];
			Ntracks 	= stuple.jet_Ntracks[iJetIndex];
			Ntowers 	= stuple.jet_Ntowers[iJetIndex];
			EtCorr 	= stuple.jet_Pt[iJetIndex];
			tlVec.SetPxPyPzE(stuple.jet_Px[iJetIndex], stuple.jet_Py[iJetIndex], stuple.jet_Pz[iJetIndex], stuple.jet_E[iJetIndex]);
			EvtPhi 		= fabs(tlVec.Phi());
			EvtEta 		= tlVec.Eta();
			SeedIEta 	= stuple.jet_SeedIEta[iJetIndex];
			SeedIPhi 	= stuple.jet_SeedIPhi[iJetIndex];
			break;
		}
		case 1:
		{
			Emfr 		= stuple.jet_up_Emfr[iJetIndex];
			DetEta 	= stuple.jet_up_DetEta[iJetIndex];
			Ntracks 	= stuple.jet_up_Ntracks[iJetIndex];
			Ntowers 	= stuple.jet_up_Ntowers[iJetIndex];
			EtCorr 	= stuple.jet_up_Pt[iJetIndex];
			tlVec.SetPxPyPzE(stuple.jet_up_Px[iJetIndex], stuple.jet_up_Py[iJetIndex], stuple.jet_up_Pz[iJetIndex], stuple.jet_up_E[iJetIndex]);
			EvtPhi 		= fabs(tlVec.Phi());
			EvtEta 		= tlVec.Eta();
			SeedIEta 	= stuple.jet_up_SeedIEta[iJetIndex];
			SeedIPhi 	= stuple.jet_up_SeedIPhi[iJetIndex];
			break;
		}
		case -1:
		{
			Emfr 		= stuple.jet_down_Emfr[iJetIndex];
			DetEta 	= stuple.jet_down_DetEta[iJetIndex];
			Ntracks 	= stuple.jet_down_Ntracks[iJetIndex];
			Ntowers 	= stuple.jet_down_Ntowers[iJetIndex];
			EtCorr 	= stuple.jet_down_Pt[iJetIndex];
			tlVec.SetPxPyPzE(stuple.jet_down_Px[iJetIndex], stuple.jet_down_Py[iJetIndex], stuple.jet_down_Pz[iJetIndex], stuple.jet_down_E[iJetIndex]);
			EvtPhi 		= fabs(tlVec.Phi());
			EvtEta 		= tlVec.Eta();
			SeedIEta 	= stuple.jet_down_SeedIEta[iJetIndex];
			SeedIPhi 	= stuple.jet_down_SeedIPhi[iJetIndex];
			break;
		}

	}

	
	hist.Emfr->Fill(Emfr, fWeight);
	hist.DetEta->Fill(DetEta, fWeight);
	hist.NTracks->Fill(Ntracks, fWeight);
	hist.NTowers->Fill(Ntowers, fWeight);
	hist.EtCorr->Fill(EtCorr, fWeight);
	hist.EvtPhi->Fill(EvtPhi, fWeight);
	hist.EvtEta->Fill(EvtEta, fWeight);
	hist.SeedIEta->Fill(SeedIEta, fWeight);
	hist.SeedIPhi->Fill(SeedIPhi, fWeight);
}


//___________________________________________________________________
void HistManager::FillPhoton1JetHists(const Stuple& stuple, const int collection,
						Histograms::Photon1JetHists_t& hist,
						const int iPhoIndex, const int iJetIndex,
						const float fWeight)
{

	TLorentzVector tlPho, tlJet, tlSum;
	float InvMass 	= 0;
	float DelPhi 	= 0;
	float DelEta 	= 0;
	float DelR 		= 0;
	float EtRatio 	= 0;
	float JetEmFr 	= 0;
	float PhoPhiWedge = 0; 

	switch (collection)
	{	
		case 0:
		{
			tlJet.SetPxPyPzE(stuple.jet_Px[iJetIndex], stuple.jet_Py[iJetIndex], stuple.jet_Pz[iJetIndex], stuple.jet_E[iJetIndex]);
			tlPho.SetPxPyPzE(stuple.pho_Px[iPhoIndex], stuple.pho_Py[iPhoIndex], stuple.pho_Pz[iPhoIndex], stuple.pho_E[iPhoIndex]);
			tlSum = tlPho + tlJet; 
			InvMass 	= tlSum.M();
			
			//std::cout << "stuple jet pt = " << stuple.jet_Pt[iJetIndex] << std::endl;
			//std::cout << "ph pto, jet pt ,invmass= " << tlPho.Pt() << "\t" <<  tlJet.Pt()  << "\t" << tlSum.M() << std::endl;

			DelPhi 	= fabs(tlPho.DeltaPhi(tlJet));
			DelEta 	= fabs(stuple.pho_DetEta[iPhoIndex] - stuple.jet_DetEta[iJetIndex]);
			DelR 		= tlPho.DeltaR(tlJet);
			EtRatio 	= (tlJet.E() * TMath::Sin(tlJet.Theta())) / (tlPho.E() * TMath::Sin(tlPho.Theta()));
			JetEmFr 	= stuple.jet_Emfr[iJetIndex];
			PhoPhiWedge = stuple.pho_PhiWedge[iPhoIndex];
			break;
		}
		case 1:
		{
			tlJet.SetPxPyPzE(stuple.jet_up_Px[iJetIndex], stuple.jet_up_Py[iJetIndex], stuple.jet_up_Pz[iJetIndex], stuple.jet_up_E[iJetIndex]);
			tlPho.SetPxPyPzE(stuple.pho_up_Px[iPhoIndex], stuple.pho_up_Py[iPhoIndex], stuple.pho_up_Pz[iPhoIndex], stuple.pho_up_E[iPhoIndex]);
			tlSum = tlPho + tlJet; 

			InvMass 	= tlSum.M();
			DelPhi 	= fabs(tlPho.DeltaPhi(tlJet));
			DelEta 	= fabs(stuple.pho_up_DetEta[iPhoIndex] - stuple.jet_up_DetEta[iJetIndex]);
			DelR 		= tlPho.DeltaR(tlJet);
			EtRatio 	= (tlJet.E() * TMath::Sin(tlJet.Theta())) / (tlPho.E() * TMath::Sin(tlPho.Theta()));
			JetEmFr 	= stuple.jet_up_Emfr[iJetIndex];
			PhoPhiWedge = stuple.pho_up_PhiWedge[iPhoIndex];
			break;
		}
		case -1:
		{
			tlJet.SetPxPyPzE(stuple.jet_down_Px[iJetIndex], stuple.jet_down_Py[iJetIndex], stuple.jet_down_Pz[iJetIndex], stuple.jet_down_E[iJetIndex]);
			tlPho.SetPxPyPzE(stuple.pho_down_Px[iPhoIndex], stuple.pho_down_Py[iPhoIndex], stuple.pho_down_Pz[iPhoIndex], stuple.pho_down_E[iPhoIndex]);
			tlSum = tlPho + tlJet; 

			InvMass 	= tlSum.M();
			DelPhi 	= fabs(tlPho.DeltaPhi(tlJet));
			DelEta 	= fabs(stuple.pho_down_DetEta[iPhoIndex] - stuple.jet_down_DetEta[iJetIndex]);
			DelR 		= tlPho.DeltaR(tlJet);
			EtRatio 	= (tlJet.E() * TMath::Sin(tlJet.Theta())) / (tlPho.E() * TMath::Sin(tlPho.Theta()));
			JetEmFr 	= stuple.jet_down_Emfr[iJetIndex];
			PhoPhiWedge = stuple.pho_down_PhiWedge[iPhoIndex];
			break;
		}

	}

	
	hist.InvMass->Fill(InvMass, fWeight);
	hist.DelPhi->Fill(DelPhi, fWeight);
	hist.DelEta->Fill(DelEta, fWeight);
	hist.DelR->Fill(DelR, fWeight);
	hist.EtRatio->Fill(EtRatio, fWeight);
	hist.JetEmfrVsDelPhi->Fill(DelPhi, JetEmFr, fWeight);
	hist.JetEmfrVsPhoPhiWedge->Fill(PhoPhiWedge, JetEmFr, fWeight);
	hist.PhoJetPtBalance->Fill(tlJet.Pt()/tlPho.Pt() - 1, fWeight);
}

//___________________________________________________________________
void HistManager::FillElectron1JetHists(const Stuple& stuple, const int collection,
						Histograms::Photon1JetHists_t& hist,
						const int iEleIndex, const int iJetIndex,
						const float fWeight)
{
	
	TLorentzVector tlEle, tlJet, tlSum;
	float InvMass 	= -1e10;
	float DelPhi 	= -1e10;
	float DelEta 	= -1e10;
	float DelR 		= -1e10;
	float EtRatio 	= -1e10;
	float JetEmFr 	= -1e10;
	float PhoPhiWedge = -1e10; 
	float EleJetBal = -1e10;


	switch (collection)
	{	
		case 0:
		{
			tlJet.SetPxPyPzE(stuple.jet_Px[iJetIndex], stuple.jet_Py[iJetIndex], stuple.jet_Pz[iJetIndex], stuple.jet_E[iJetIndex]);
			tlEle.SetPxPyPzE(stuple.ele_Px[iEleIndex], stuple.ele_Py[iEleIndex], stuple.ele_Pz[iEleIndex], stuple.ele_E[iEleIndex]);
			tlSum = tlEle + tlJet; 
			InvMass 	= tlSum.M();
			
			//std::cout << "stuple jet pt = " << stuple.jet_Pt[iJetIndex] << std::endl;
			//std::cout << "ph pto, jet pt ,invmass= " << tlEle.Pt() << "\t" <<  tlJet.Pt()  << "\t" << tlSum.M() << std::endl;

			DelPhi 	= fabs(tlEle.DeltaPhi(tlJet));
			DelEta 	= fabs(stuple.ele_DetEta[iEleIndex] - stuple.jet_DetEta[iJetIndex]);
			DelR 		= tlEle.DeltaR(tlJet);
			EtRatio 	= (tlJet.E() * TMath::Sin(tlJet.Theta())) / (tlEle.E() * TMath::Sin(tlEle.Theta()));
			JetEmFr 	= stuple.jet_Emfr[iJetIndex];
			PhoPhiWedge = stuple.ele_PhiWedge[iEleIndex];
			if (tlEle.Pt()>0) EleJetBal = tlJet.Pt()/tlEle.Pt() - 1;
			break;
		}
		case 1:
		{
			tlJet.SetPxPyPzE(stuple.jet_up_Px[iJetIndex], stuple.jet_up_Py[iJetIndex], stuple.jet_up_Pz[iJetIndex], stuple.jet_up_E[iJetIndex]);
			tlEle.SetPxPyPzE(stuple.ele_up_Px[iEleIndex], stuple.ele_up_Py[iEleIndex], stuple.ele_up_Pz[iEleIndex], stuple.ele_up_E[iEleIndex]);
			tlSum = tlEle + tlJet; 

			InvMass 	= tlSum.M();
			DelPhi 	= fabs(tlEle.DeltaPhi(tlJet));
			DelEta 	= fabs(stuple.ele_up_DetEta[iEleIndex] - stuple.jet_up_DetEta[iJetIndex]);
			DelR 		= tlEle.DeltaR(tlJet);
			EtRatio 	= (tlJet.E() * TMath::Sin(tlJet.Theta())) / (tlEle.E() * TMath::Sin(tlEle.Theta()));
			JetEmFr 	= stuple.jet_up_Emfr[iJetIndex];
			PhoPhiWedge = stuple.ele_up_PhiWedge[iEleIndex];
			if (tlEle.Pt()>0) EleJetBal = tlJet.Pt()/tlEle.Pt() - 1;
			break;
		}
		case -1:
		{
			tlJet.SetPxPyPzE(stuple.jet_down_Px[iJetIndex], stuple.jet_down_Py[iJetIndex], stuple.jet_down_Pz[iJetIndex], stuple.jet_down_E[iJetIndex]);
			tlEle.SetPxPyPzE(stuple.ele_down_Px[iEleIndex], stuple.ele_down_Py[iEleIndex], stuple.ele_down_Pz[iEleIndex], stuple.ele_down_E[iEleIndex]);
			tlSum = tlEle + tlJet; 

			InvMass 	= tlSum.M();
			DelPhi 	= fabs(tlEle.DeltaPhi(tlJet));
			DelEta 	= fabs(stuple.ele_down_DetEta[iEleIndex] - stuple.jet_down_DetEta[iJetIndex]);
			DelR 		= tlEle.DeltaR(tlJet);
			EtRatio 	= (tlJet.E() * TMath::Sin(tlJet.Theta())) / (tlEle.E() * TMath::Sin(tlEle.Theta()));
			JetEmFr 	= stuple.jet_down_Emfr[iJetIndex];
			PhoPhiWedge = stuple.ele_down_PhiWedge[iEleIndex];
			if (tlEle.Pt()>0) EleJetBal = tlJet.Pt()/tlEle.Pt() - 1;
			break;
		}

	}
	
	hist.InvMass->Fill(InvMass, fWeight);
	hist.DelPhi->Fill(DelPhi, fWeight);
	hist.DelEta->Fill(DelEta, fWeight);
	hist.DelR->Fill(DelR, fWeight);
	hist.EtRatio->Fill(EtRatio, fWeight);
	hist.JetEmfrVsDelPhi->Fill(DelPhi, JetEmFr, fWeight);
	hist.JetEmfrVsPhoPhiWedge->Fill(PhoPhiWedge, JetEmFr, fWeight);
	//std::cout << "balance = " << (tlJet.Pt()/tlEle.Pt() - 1) << std::endl;
	hist.PhoJetPtBalance->Fill(EleJetBal, fWeight);
	//hist.PhoJetPtBalance->Print();

}




//___________________________________________________________________
void HistManager::FillPhoton2JetsHists(const Stuple& stuple, const int collection,
						Histograms::Photon2JetsHists_t& hist,
						const int iPhoIndex, const int iJet1Index,
						const int iJet2Index, const float fWeight)
{
	TLorentzVector tlJet1, tlJet2, tlPho, tlSum;

	switch (collection)
	{
		case 0:
		{
			tlJet1.SetPxPyPzE(stuple.jet_Px[iJet1Index], stuple.jet_Py[iJet1Index], stuple.jet_Pz[iJet1Index], stuple.jet_E[iJet1Index]);
			tlJet2.SetPxPyPzE(stuple.jet_Px[iJet2Index], stuple.jet_Py[iJet2Index], stuple.jet_Pz[iJet2Index], stuple.jet_E[iJet2Index]);
			tlPho.SetPxPyPzE(stuple.pho_Px[iPhoIndex], stuple.pho_Py[iPhoIndex], stuple.pho_Pz[iPhoIndex], stuple.pho_E[iPhoIndex]);
			break;
		}
		case 1:
		{
			tlJet1.SetPxPyPzE(stuple.jet_up_Px[iJet1Index], stuple.jet_up_Py[iJet1Index], stuple.jet_up_Pz[iJet1Index], stuple.jet_up_E[iJet1Index]);
			tlJet2.SetPxPyPzE(stuple.jet_up_Px[iJet2Index], stuple.jet_up_Py[iJet2Index], stuple.jet_up_Pz[iJet2Index], stuple.jet_up_E[iJet2Index]);
			tlPho.SetPxPyPzE(stuple.pho_up_Px[iPhoIndex], stuple.pho_up_Py[iPhoIndex], stuple.pho_up_Pz[iPhoIndex], stuple.pho_up_E[iPhoIndex]);
			break;
		}
		case -1:
		{
			tlJet1.SetPxPyPzE(stuple.jet_down_Px[iJet1Index], stuple.jet_down_Py[iJet1Index], stuple.jet_down_Pz[iJet1Index], stuple.jet_down_E[iJet1Index]);
			tlJet2.SetPxPyPzE(stuple.jet_down_Px[iJet2Index], stuple.jet_down_Py[iJet2Index], stuple.jet_down_Pz[iJet2Index], stuple.jet_down_E[iJet2Index]);
			tlPho.SetPxPyPzE(stuple.pho_down_Px[iPhoIndex], stuple.pho_down_Py[iPhoIndex], stuple.pho_down_Pz[iPhoIndex], stuple.pho_down_E[iPhoIndex]);
			break;
		}
	}
	
	tlSum = tlPho + tlJet1 + tlJet2; 
	hist.InvMass->Fill(tlSum.M(), fWeight);
}

//___________________________________________________________________
void HistManager::FillElectron2JetsHists(const Stuple& stuple, const int collection,
						Histograms::Photon2JetsHists_t& hist,
						const int iEleIndex, const int iJet1Index,
						const int iJet2Index, const float fWeight)
{

	TLorentzVector tlJet1(0,0,0,0), tlJet2(0,0,0,0), tlEle(0,0,0,0), tlSum(0,0,0,0);

	switch (collection)
	{
		case 0:
		{
			tlJet1.SetPxPyPzE(stuple.jet_Px[iJet1Index], stuple.jet_Py[iJet1Index], stuple.jet_Pz[iJet1Index], stuple.jet_E[iJet1Index]);
			tlJet2.SetPxPyPzE(stuple.jet_Px[iJet2Index], stuple.jet_Py[iJet2Index], stuple.jet_Pz[iJet2Index], stuple.jet_E[iJet2Index]);
			tlEle.SetPxPyPzE(stuple.ele_Px[iEleIndex], stuple.ele_Py[iEleIndex], stuple.ele_Pz[iEleIndex], stuple.ele_E[iEleIndex]);
			break;
		}
		case 1:
		{
			tlJet1.SetPxPyPzE(stuple.jet_up_Px[iJet1Index], stuple.jet_up_Py[iJet1Index], stuple.jet_up_Pz[iJet1Index], stuple.jet_up_E[iJet1Index]);
			tlJet2.SetPxPyPzE(stuple.jet_up_Px[iJet2Index], stuple.jet_up_Py[iJet2Index], stuple.jet_up_Pz[iJet2Index], stuple.jet_up_E[iJet2Index]);
			tlEle.SetPxPyPzE(stuple.ele_up_Px[iEleIndex], stuple.ele_up_Py[iEleIndex], stuple.ele_up_Pz[iEleIndex], stuple.ele_up_E[iEleIndex]);
			break;
		}
		case -1:
		{
			tlJet1.SetPxPyPzE(stuple.jet_down_Px[iJet1Index], stuple.jet_down_Py[iJet1Index], stuple.jet_down_Pz[iJet1Index], stuple.jet_down_E[iJet1Index]);
			tlJet2.SetPxPyPzE(stuple.jet_down_Px[iJet2Index], stuple.jet_down_Py[iJet2Index], stuple.jet_down_Pz[iJet2Index], stuple.jet_down_E[iJet2Index]);
			tlEle.SetPxPyPzE(stuple.ele_down_Px[iEleIndex], stuple.ele_down_Py[iEleIndex], stuple.ele_down_Pz[iEleIndex], stuple.ele_down_E[iEleIndex]);
			break;
		}
	}
	
	tlSum = tlEle + tlJet1 + tlJet2; 
	if (tlSum.M()<0) 
	{
		std::cout << __FILE__ << ":" << __FUNCTION__ << ": Event with inv mass < 0!";
		std::cout << stuple.evt_RunNumber << ", " << stuple.evt_EventNumber << std::endl;
	}
	hist.InvMass->Fill(tlSum.M(), fWeight);
}




//___________________________________________________________________
void HistManager::FillTwoJetsHists(const Stuple& stuple, const int collection,
						Histograms::TwoJetsHists_t& hist,
						const int iJet1Index, const int iJet2Index,
						const float fWeight)
{

	TLorentzVector tlJet1, tlJet2, tlSum;
	float DelEta = 0;

	switch (collection)
	{
		case 0:
		{
			tlJet1.SetPxPyPzE(stuple.jet_Px[iJet1Index], stuple.jet_Py[iJet1Index], stuple.jet_Pz[iJet1Index], stuple.jet_E[iJet1Index]);
			tlJet2.SetPxPyPzE(stuple.jet_Px[iJet2Index], stuple.jet_Py[iJet2Index], stuple.jet_Pz[iJet2Index], stuple.jet_E[iJet2Index]);
			DelEta = fabs(stuple.jet_DetEta[iJet1Index] - stuple.jet_DetEta[iJet2Index]);
			break;
		}
		case 1:
		{
			tlJet1.SetPxPyPzE(stuple.jet_up_Px[iJet1Index], stuple.jet_up_Py[iJet1Index], stuple.jet_up_Pz[iJet1Index], stuple.jet_up_E[iJet1Index]);
			tlJet2.SetPxPyPzE(stuple.jet_up_Px[iJet2Index], stuple.jet_up_Py[iJet2Index], stuple.jet_up_Pz[iJet2Index], stuple.jet_up_E[iJet2Index]);
			DelEta = fabs(stuple.jet_up_DetEta[iJet1Index] - stuple.jet_up_DetEta[iJet2Index]);
			break;
		}
		case -1:
		{
			tlJet1.SetPxPyPzE(stuple.jet_down_Px[iJet1Index], stuple.jet_down_Py[iJet1Index], stuple.jet_down_Pz[iJet1Index], stuple.jet_down_E[iJet1Index]);
			tlJet2.SetPxPyPzE(stuple.jet_down_Px[iJet2Index], stuple.jet_down_Py[iJet2Index], stuple.jet_down_Pz[iJet2Index], stuple.jet_down_E[iJet2Index]);
			DelEta = fabs(stuple.jet_down_DetEta[iJet1Index] - stuple.jet_down_DetEta[iJet2Index]);
			break;
		}
	}

	
	tlSum = tlJet1 + tlJet2; 
	float fEtRatio = (tlJet2.E() * TMath::Sin(tlJet2.Theta())) / (tlJet1.E() * TMath::Sin(tlJet1.Theta()));
	float DelPhi = fabs(tlJet1.DeltaPhi(tlJet2));
	float DelR = tlJet1.DeltaR(tlJet2);
	
	hist.InvMass->Fill(tlSum.M(), fWeight);
	hist.DelPhi->Fill(DelPhi, fWeight);
	hist.DelEta->Fill(DelEta, fWeight);
	hist.DelR->Fill(DelR, fWeight);
	hist.EtRatio->Fill(fEtRatio, fWeight);
}

