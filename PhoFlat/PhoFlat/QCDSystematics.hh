#ifndef QCDSYSTEMATICS_HH
#define QCDSYSTEMATICS_HH

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
#include "TH1F.h"

class QCDSystematics
{
	public:
		QCDSystematics();
		int iProgressBy;
	
		void Loop(TChain* chain, int iRunEvents);
		void Init(TChain* chain);
		void DoMyStuff(Stuple& stuple);
		void CleanUp();
		void BookHistograms();
		std::string GetHistFileName() const { return sHistFileName; }
		void SetHistFileName(std::string name)  { sHistFileName = name; }
		void SetTightenUpCut(int cut) { 
			iCutUsed = cut;
		}
		int GetTightenUpCut() const { return iCutUsed; }
		void PrintJobSettings();
		void SetMinPhoEt(const float et) { fMinPhoEt = et; }
		void SetMaxPhoEt(const float et) { fMaxPhoEt = et; }
		void SetMaxPhoEta(const float eta) { fMaxPhoEta = eta; }
		void SetMinJetEt(const float et) { fMinJetEt = et; }
		void SetMaxJetEta(const float eta) { fMaxJetEta = eta; }
		float GetMinPhoEt() const { return fMinPhoEt; }
		float GetMaxPhoEt() const { return fMaxPhoEt; }
		float GetMaxPhoEta() const { return fMaxPhoEta; }
		float GetMinJetEt() const { return fMinJetEt; }
		float GetMaxJetEta() const { return fMaxJetEta; }
		void SetMinMet(const float m) { fMinMet = m; }
		float GetMinMet() const { return fMinMet; }

		
	private:
		int iEventsSelected;				// final number of events
		unsigned int iPho1JetEvts, iPho2JetEvts;
		std::string sHistFileName;		//name of the root file to save the hists
		Stuple stuple;
		TChain* myChain;

		Histograms::EventHists_t GA_Evt, GB_Evt;
		Histograms::PhotonHists_t GA_Pho, GB_Pho;
		Histograms::JetHists_t GA_Jet, GB_Jet1, GB_Jet2;
		Histograms::Photon1JetHists_t GA_PhoJet, GB_PhoJet1, GB_PhoJet2;
		Histograms::Photon2JetsHists_t GB_PhoJets;		//jet>=2 case, pho+2jets
		Histograms::TwoJetsHists_t GB_TwoJets;

		TH1F* hCounter;

		TFile *rootFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		HistManager myHistMan;
		TH1F *T1,*T2,*T3;
		TH1F *haloPhiWedge_1j, *haloEmTime_1j;
		TH1F *haloPhiWedge_2j, *haloEmTime_2j;

		int iCutUsed; //set only of the loose cut to tight
		float fMinPhoEt;
		float fMaxPhoEt;
		float fMaxPhoEta;
		float fMinJetEt;
		float fMaxJetEta;
		float fMinMet;
};
#endif
