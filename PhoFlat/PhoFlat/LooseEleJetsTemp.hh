#ifndef LOOSEELEJETSTEMP_HH
#define LOOSEELEJETSTEMP_HH

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

class LooseEleJetsTemp
{
	public:
		LooseEleJetsTemp();
		int iProgressBy;
	
		void Main(TChain* chain, int iRunEvents = 0);
		void Init(TChain* chain);
		void DoMyStuff(Stuple& stuple);
		void CleanUp();
		void BookHistograms();
		std::string GetHistFileName() const { return sHistFileName; }
		void SetHistFileName(std::string name)  { sHistFileName = name; }
		
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
		TH1F *T1,*T2,*T3;

		float fSFR_1j, fSFR_m_1j, fSFR_p_1j;			// sum of fake rates
		float fSFR_2j, fSFR_m_2j, fSFR_p_2j;			// sum of fake rates
		
		
};
#endif
