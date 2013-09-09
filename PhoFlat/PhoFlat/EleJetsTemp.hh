#ifndef ELEJETSTEMP_HH
#define ELEJETSTEMP_HH

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
#include "../RootTools/IOColors.hh"
#include "../RootTools/CommonTools.hh"

class EleJetsTemp
{
	public:
		EleJetsTemp();
		int iProgressBy;
	
		void Main(TChain* chain, int iRunEvents = 0);
		void Init(TChain* chain);
		void DoMyStuff(Stuple& stuple);
		void CleanUp();
		void BookHistograms();
		std::string GetHistFileName() const { return sHistFileName; }
		void SetHistFileName(std::string name)  { sHistFileName = name; }
		void DoDump(Stuple& stuple);
		void SetMaxEleDetEta(const float eta) { fMaxEleEta = fabs(eta);}
		void SetMinEleEt(const float et) { fMinEleEt = fabs(et); }
		float GetMaxEleDetEta() const { return fMaxEleEta;}
		float GetMinEleEt() const { return fMinEleEt; }
		void SetMcFlag(const int ii) { bMcFlag = (bool) ii; }
		bool GetMcFlag() const { return bMcFlag; }
		void PrintJobSettings();
		float GetMinMet() const { return fMinMet; }
		void SetMinMet(const float met) { fMinMet = met; }
		int GetMinClass12Vtx() const { return iMinClass12Vtx; }
		int GetMaxClass12Vtx() const { return iMaxClass12Vtx; }
		int GetMinNjet15() const { return iMinNjet15; }
		void SetSignalTimes(const float tmin, const float tmax) { fMinSignalEmTime = tmin; fMaxSignalEmTime = tmax; }
		void SetMinClass12Vtx(const int nvt) { iMinClass12Vtx = nvt; }
		void SetMaxClass12Vtx(const int nvt) { iMaxClass12Vtx = nvt; }
		void SetMinNjet15(const int nj) { iMinNjet15 = nj; }
		float GetSignalTimeLoLimit() const { return fMinSignalEmTime; }
		float GetSignalTimeUpLimit() const { return fMaxSignalEmTime; }

		float GetMinJetEt() const { return fMinJetEt; }
		float GetMaxEleEta() const { return fMaxEleEta; }
		float GetMaxJetEta() const { return fMaxJetEta; }
		void SetMaxEleEta(const float e) { fMaxEleEta = fabs(e); }
		void SetReportProgress(const int t) { iProgressBy = t; }
		int GetReportProgress() const { return iProgressBy; }
	private:
		int iEventsSelected;				// final number of events
		unsigned int iEle1JetEvts, iEle2JetEvts;
		std::string sHistFileName;		//name of the root file to save the hists
		Stuple stuple;
		TChain* myChain;

		Histograms::EventHists_t GA_Evt, GB_Evt;
		Histograms::PhotonHists_t GA_Pho, GB_Pho;
		Histograms::JetHists_t GA_Jet, GB_Jet1, GB_Jet2;
		Histograms::Photon1JetHists_t GA_PhoJet, GB_PhoJet1, GB_PhoJet2;
		Histograms::Photon2JetsHists_t GB_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t GB_TwoJets;

		TFile *rootFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		HistManager myHistMan;
		TH2F *delREt;
		TH1F *T1,*T2,*T3;

		float fSFR_1j, fSFR_m_1j, fSFR_p_1j;			// sum of fake rates
		float fSFR_2j, fSFR_m_2j, fSFR_p_2j;			// sum of fake rates
		
		unsigned int iEvtsDelR2;
		
		float fMaxEleEta;
		float fMinEleEt;
		bool bMcFlag;			//must be set by user for different datasets, default is true, to make it safe
		float fMinMet;
		float fMinSignalEmTime;
		int iMinClass12Vtx;
		int iMaxClass12Vtx;
		float fMaxSignalEmTime;
		float fMaxVtxz;
		int iMinNjet15;
		float fMinJetEt;
		float fMaxJetEta;
//		std::vector<std::string> vJobSummary;
		StringVector sVec;
};
#endif
