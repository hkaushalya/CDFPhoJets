#ifndef CosmicNoVtx_HH
#define CosmicNoVtx_HH

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

class CosmicNoVtx
{
	public:
		CosmicNoVtx();
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

		// A== 1JET , B== 2JETS
		std::vector<Histograms::EventHists_t> GA_Evt, GB_Evt;
		std::vector<Histograms::PhotonHists_t> GA_Pho, GB_Pho;
		std::vector<Histograms::JetHists_t> GA_Jet, GB_Jet1, GB_Jet2;
		std::vector<Histograms::Photon1JetHists_t> GA_PhoJet, GB_PhoJet1,GB_PhoJet2;
		std::vector<Histograms::Photon2JetsHists_t> GB_PhoJets;		//jet>=2 case, pho+2jets
		std::vector<Histograms::TwoJetsHists_t> GB_TwoJets;


		std::vector<TFolder*> halofolder;
		TFile *rootFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		HistManager myHistMan;

};
#endif
