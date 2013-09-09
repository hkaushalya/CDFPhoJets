#ifndef ZJETSEL_HH
#define ZJETSEL_HH

/**********************************************************
 * $Id: ZJetSel.hh,v 1.1 2011/05/31 17:59:51 samantha Exp $
 * $Log: ZJetSel.hh,v $
 * Revision 1.1  2011/05/31 17:59:51  samantha
 * This is to get corrections to EWK MC by comparing Z+Jets data and W MC.
 * I am simply treating Z as the photon and applying cuts similar to the photon.
 * Also standard mass windows cuts on the di-electron mass to reduce the
 * bacgrounds.
 *
 **********************************************************
 *Samantha K. Hewamanage
 *********************************************************/

#include "Stuple.hh"
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include <vector>
#include "TFolder.h"
#include "Histograms.hh"
#include "TDirectory.h"
#include "ReadInAll.hh"
#include "HistManager.hh"
#include "FreeFunctions.hh"
#include "PhotonList.hh"
#include "ElectronList.hh"
#include "JetList.hh"
#include "CommonVars.hh"
#include "TF1.h"
#include "TProfile.h"
#include "../RootTools/IOColors.hh"
#include <cmath>

class ZJetSel
{
	public:
		ZJetSel();
		int iProgressBy;
	
		struct Hist_t {
			Histograms::EventHists_t Evt;
			Histograms::PhotonHists_t Ele1, Ele2;
			Histograms::JetHists_t Jet;
			Histograms::TwoJetsHists_t Z;
			Histograms::TwoJetsHists_t ZJet;
			/*TH1F* z1j_zpt;
			TH1F* z1j_zet;
			TH1F* z1j_zmass;
			TH1F* z1j_zeta;
			TH1F* z1j_zj_mass;
			TH1F* z2j_zpt;
			TH1F* z2j_zet;
			TH1F* z2j_zmass;
			TH1F* z2j_zeta;
			TH1F* z2j_z1j_mass;
			TH1F* z2j_z2j_mass;
			TH1F* z2j_zj1j2_mass;*/

		};

		struct EvtTag_t {
			int Run;
			int Evt;
		};

		void Main(TChain* chain, int iRunEvents = 0);
		void Init(TChain* chain);
		void DoMyStuff(Stuple& stuple);
		void CleanUp();
		void BookHistograms();
		std::string GetHistFileName() const { return sHistFileName; }
		void SetHistFileName(std::string name)  { sHistFileName = name; }


	void SetMcFlag(const int ii) { bMcFlag = (bool) ii; }
	bool GetMcFlag() const { return bMcFlag; }
	void CreatePhotonLists(const Stuple st, PhotonList& pholist, const int collection);
	void CreateElectronLists(const Stuple st, ElectronList& elelist, const int collection);
	void CreateJetLists(const Stuple st, JetList& jetlist, const int collection);
	void DoMyStuff(CommonVars cVars, ElectronList cEle, JetList cJets, PhotonList cPhos);
	void GetCommonVariables(const Stuple stuple, CommonVars& commVars);
	int JetSelAndHistFill(const CommonVars Vars, const ElectronList eles, 
											const JetList jets, 
											std::vector<int> vEleIndex, 
											Hist_t* hist,
											float fWeight=1, const bool debug=false);
	
	//hist filling
		void FillEventHists(const CommonVars&, const ElectronList&, Histograms::EventHists_t* , 
								const float fWeight=1);
		void FillElectronHists(const CommonVars&, const ElectronList&, Histograms::PhotonHists_t*,
								const int iEleIndex, const float fWeight=1, const bool debug=false);									
		void FillZHists(const CommonVars&, const ElectronList&, 
								Histograms::TwoJetsHists_t* hist,
								const int iEle1Index, const int iEle2Index,
								const float fWeight=1);

		void FillZJetHists(const CommonVars& cVars, const ElectronList& eles, const JetList& jets, 
				Histograms::TwoJetsHists_t* hist, const int iEle1Index, const int iEle2Index, const int iJetIndex,
				const float fWeight);
		void FillJetHists(const CommonVars& cVars, const JetList& jets, Histograms::JetHists_t* hist, 
				const int iJetIndex,	const float fWeight);

		void PrintHeader(const CommonVars&) const;
		void SetReportProgress(const unsigned int l) { iProgressBy = l; }
		
		void SetMaxElectronEt(const float et) { fMaxEleEt = et; }
		int GetMinClass12Vtx() const { return iMinClass12Vtx; }
		int GetMaxClass12Vtx() const { return iMaxClass12Vtx; }
		
		void PrintJobSettings();
		float GetMinMet() const { return fMinMet; }
		void SetMinMet(const float met) { fMinMet = met; }
		float GetMaxMet() const { return fMaxMet; }
		void SetMaxMet(const float met) { fMaxMet = met; }
		void SetUseNvtxWeights(const bool b) { bUseNvtxWgts = b; }
		bool UseNvtxWeights() const { return bUseNvtxWgts; }
		void SetMinZPt(const float et) { fMinZpt = et; }
		float GetMinZPt() const { return fMinZpt; }
		void SetMaxZEta(const float eta) { fMaxZeta = eta; }
		float GetMaxZEta() const { return fMaxZeta; }
		float GetMinJetEt() const { return fMinJetEt; }
		void SetMinJetEt(const float et) { fMinJetEt = et; }
		float GetMaxJetEta() const { return fMaxJetEta; }
		void SetMaxJetEta(const float eta) { fMaxJetEta = eta; }

	private:
		std::string sHistFileName;		//name of the root file to save the hists
		Stuple stuple;
		TChain* myChain;
		int iCurrTree;

		TFile *histFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		
		TH1F *phoEt_iso25, *phoEt_50, *phoEt_70;

		Hist_t Hmc;		///for pho mc, ewk mc (CENTRAL/JES UP/DOWN)
		
		bool bMcFlag;			//must be set by user for different datasets, default is true, to make it safe

		std::vector<unsigned int> vMcCount;

		std::vector<EvtTag_t> EVT;
		float fMinEleEt;
		float fMaxEleEt;
		int iMinClass12Vtx;
		int iMaxClass12Vtx;
		float fMaxVtxz;
		float fMinMet, fMaxMet;
		int iExtraEle;
		std::vector<float> vNvtxWgts;
		bool bUseNvtxWgts;

		TH1F *hClosestJetMetDelPhi_b4, *hClosestJetMetDelPhi_a4;	
		TH1F *hZMetDelPhi_b4, * hZMetDelPhi_a4;
		TH1F *hClosestJetPt_a4;
		float fMinZpt, fMaxZeta;
		float fMinJetEt, fMaxJetEta;
};
#endif
