#ifndef HALOSTUDY_HH
#define HALOSTUDY_HH

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

class HaloStudy
{
	public:
		HaloStudy();
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
		std::string sHistFileName;		//name of the root file to save the hists
		Stuple stuple;
		TChain* myChain;
		std::vector<int> vHalosFound;
		//naming convections
		// G = Met/Et>1  L = Met/Et<1
		// GA_ = 1 jet case
		// GB_ = 2 jet case etc
		
		std::vector<TH1F*> MetEtRatioHist;

		// for MEtEtRatio > 1  : halo type 0-5
		std::vector<Histograms::EventHists_t> GA_Evt, GB_Evt;
		std::vector<Histograms::PhotonHists_t> GA_Pho, GB_Pho;
		std::vector<Histograms::JetHists_t> GA_Jet, GB_Jet1, GB_Jet2;
		std::vector<Histograms::Photon1JetHists_t> GA_PhoJet, GB_PhoJet1, GB_PhoJet2;
		std::vector<Histograms::Photon2JetsHists_t> GB_PhoJets;		//jet>=2 case, pho+2jets
		std::vector<Histograms::TwoJetsHists_t> GB_TwoJets;

		// for MEtEtRatio < 1  : halo type 0-5
		std::vector<Histograms::EventHists_t> LA_Evt, LB_Evt;
		std::vector<Histograms::PhotonHists_t> LA_Pho, LB_Pho;
		std::vector<Histograms::JetHists_t> LA_Jet, LB_Jet1, LB_Jet2;
		std::vector<Histograms::Photon1JetHists_t> LA_PhoJet, LB_PhoJet1, LB_PhoJet2;
		std::vector<Histograms::Photon2JetsHists_t> LB_PhoJets;		//jet>=2 case, pho+2jets
		std::vector<Histograms::TwoJetsHists_t> LB_TwoJets;

		std::vector<TFolder*> halofolder;
		TFile *rootFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		HistManager myHistMan;

};
#endif
