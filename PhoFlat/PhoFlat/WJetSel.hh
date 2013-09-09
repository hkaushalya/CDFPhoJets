#ifndef WJETSEL_HH
#define WJETSEL_HH
/*{{{*/
/**********************************************************
 * $Id: WJetSel.hh,v 1.2 2011/04/25 04:01:49 samantha Exp $
 * $Log: WJetSel.hh,v $
 * Revision 1.2  2011/04/25 04:01:49  samantha
 * ADDED: 1. Hists for jet2.
 * 2. Hists to see jet/met separation.
 * 2. Count number of events as cuts are applied.
 *
 * Revision 1.1  2010/08/23 03:05:18  samantha
 * W+Jet events are studied to derive additional corrections to EWK MC samples
 * to improve the MET predictions. W+jets are similar to g+jets. I use Pt and Eta
 * cuts similar to the g, on the W. I also apply a standard W transverse mass cut
 * to reduce the backgrounds.
 *
 **********************************************************
 *Samantha K. Hewamanage
 *********************************************************/
/*}}}*/

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

class WJetSel
{
	public:
		WJetSel();
		int iProgressBy;
	
		struct Hist_t {
			Histograms::EventHists_t Evt;
			Histograms::PhotonHists_t Ele1;
			Histograms::JetHists_t Jet, Jet2;
			Histograms::TwoJetsHists_t W, WJet;
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
	void CreateElectronLists(const Stuple st, ElectronList& elelist, const int collection);
	void CreatePhotonLists(const Stuple st, PhotonList& pholist, const int collection);
	void CreateJetLists(const Stuple st, JetList& jetlist, const int collection);
	void DoMyStuff(CommonVars cVars, ElectronList cEle);
	void GetCommonVariables(const Stuple stuple, CommonVars& commVars);
	int JetSelAndHistFill(const CommonVars Vars, const ElectronList eles, const JetList jets, 
											std::vector<int> vEleIndex, 
											Hist_t* hist,
											std::vector<unsigned int>& count, 
											float fWeight=1, const bool debug=false);
	
	//hist filling
		void FillEventHists(const CommonVars&, const ElectronList&, 
								const JetList&,
								Histograms::EventHists_t* , 
								const float fWeight=1);
		void FillElectronHists(const CommonVars&, const ElectronList&, Histograms::PhotonHists_t*,
								const int iEleIndex, const float fWeight=1, const bool debug=false);									
		void FillJetHists(const CommonVars& cVars, const JetList& jets, Histograms::JetHists_t* hist, 
								const int iJetIndex,	const float fWeight);
		void FillWHists(const CommonVars&, const ElectronList&, 
								Histograms::TwoJetsHists_t* hist,
								const int iEle1Index, const float fWeight=1);
		void FillWJetHists(const CommonVars&, const ElectronList& eles, const JetList& jets, 
								Histograms::TwoJetsHists_t* hist,
								const int iEle1Index, const float fWeight=1);

		void PrintHeader(const CommonVars&) const;
		void SetReportProgress(const unsigned int l) { iProgressBy = l; }
		
		int GetMinClass12Vtx() const { return iMinClass12Vtx; }
		int GetMaxClass12Vtx() const { return iMaxClass12Vtx; }
		
		void PrintJobSettings();
		float GetMinMet() const { return fMinMet; }
		void SetMinMet(const float met) { fMinMet = met; }
		float GetMaxMet() const { return fMaxMet; }
		void SetMaxMet(const float met) { fMaxMet = met; }
		void SetUseNvtxWeights(const int b) { iUseNvtxWgts = b; }
		int UseNvtxWeights() const { return iUseNvtxWgts; }
		float GetMinWPt() const { return fMinWpt; }
		void SetMinWPt(const float pt) { fMinWpt = pt; }
		float GetMinJetEt() const { return fMinJetEt; }
		void SetMinJetEt(const float et) { fMinJetEt = et; }
		float GetMaxJetEta() const { return fMaxJetEta; }
		void SetMaxJetEta(const float eta) { fMaxJetEta = eta; }
		float GetMinEleEt() const { return fMinEleEt; }
		void SetMinEleEt(const float et) { fMinEleEt = et; }
		float GetMaxEleEt() const { return fMaxEleEt; } 
		void SetMaxEleEt(const float et) { fMaxEleEt = et; }
		float GetMaxEleEta() const { return fMaxEleEta; }
		void SetMaxEleEta(const float eta) { fMaxEleEta = eta; }
		int GetDataset() const { return iDataset;}
		void  SetDataset(const int t) { iDataset = t;}
			
		enum counts {
			iMet_cut = 1,
			ielenum_cut = 2,
			ivtx_cut = 3,
			ivtxz_cut = 4,
			itele_cut = 5,
			iWpt_cut = 6,
			iWeta_cut = 7,
			iWmass_cut = 8,
			iWjet_cut = 9,
			iMetCleanup_cut = 10,
		};
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
		float fMaxEleEta;
		int iMinClass12Vtx;
		int iMaxClass12Vtx;
		float fMaxVtxz;
		float fMinMet, fMaxMet;
		int iExtraEle;
		std::vector<float> vNvtxWgts;
		int iUseNvtxWgts;
		float fMinWpt, fMaxWpt;
		float fMinJetEt, fMaxJetEta;
		unsigned njet15counterr;
		int iDataset;
		TH1F *hClosestJetMetDelPhi_b4, *hClosestJetMetDelPhi_a4;	
		TH1F *hWMetDelPhi_b4, * hWMetDelPhi_a4;
		TH1F *hClosestJetPt_a4;
		unsigned iCount[11];

};
#endif
