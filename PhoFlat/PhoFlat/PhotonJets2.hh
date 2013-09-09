#ifndef PhotonJets22_HH
#define PhotonJets22_HH

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

class PhotonJets2
{
	public:
		PhotonJets2();
		int iProgressBy;
	

		struct Hist_t {
			Histograms::EventHists_t p1j_Evt, p2j_Evt;				//
			Histograms::PhotonHists_t p1j_Pho, p2j_Pho;
			Histograms::JetHists_t p1j_Jet, p2j_Jet1, p2j_Jet2;
			Histograms::Photon1JetHists_t p1j_PhoJet, p2j_PhoJet1, p2j_PhoJet2;
			Histograms::Photon2JetsHists_t p2j_PhoJets;		//jet>=2 case, pho+2jets
			Histograms::TwoJetsHists_t p2j_TwoJets;
		};

		struct EvtTag_t {
			int Run;
			int Evt;
			int Halo;
			int Cosmic;
			int Signal;
			int QCD;
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
	void DoMyStuff(CommonVars cVars, PhotonList cPho, PhotonList uPho, PhotonList dPho,
									ElectronList cEle, ElectronList uEle, ElectronList dEle,
									JetList cJet, JetList uJet, JetList dJet);
	void GetCommonVariables(const Stuple stuple, CommonVars& commVars);
	int JetSelAndHistFill(const CommonVars Vars, const PhotonList phos, const ElectronList eles, 
											const JetList jets, Hist_t* hist,
											const int iPhoIndex, 
											std::vector<unsigned int>& count, 
											const bool UseCprWgtForHt=false, 
											const float fWeight=1);

	//hist filling
		void FillEventHists(const CommonVars&, const PhotonList&, const ElectronList&, const JetList&, Histograms::EventHists_t* , 
								const int iPhoIndex, const float fWeight=1, const bool UseCprWgtForHt=false);
		void FillPhotonHists(const CommonVars&, const PhotonList&, Histograms::PhotonHists_t*,
								const int iPhoIndex, const float fWeight=1);									
		void FillElectronHists(const ElectronList&, Histograms::PhotonHists_t* ,
								const int iEleIndex, const float fWeight=1);
		void FillJetHists(const CommonVars&, const JetList&, Histograms::JetHists_t* , 
								const int iJetIndex,	const float fWeight=1);
		void FillPhoton1JetHists(const CommonVars&, const PhotonList&, const JetList&, Histograms::Photon1JetHists_t* ,
								const int iPhoIndex, const int iJetIndex,
								const float fWeight=1);
		void FillPhoton2JetsHists(const CommonVars&, const PhotonList&, const JetList&, Histograms::Photon2JetsHists_t* hist,
								const int iPhoIndex, const int iJet1Index,
								const int iJet2Index, const float fWeight=1);
		void FillTwoJetsHists(const CommonVars&, const JetList&, Histograms::TwoJetsHists_t* hist,
								const int iJet1Index, const int iJet2Index,
								const float fWeight=1);

		void PrintHeader(const CommonVars&);
		void SetReportProgress(const unsigned int l) { iProgressBy = l; }
		
	private:
		std::string sHistFileName;		//name of the root file to save the hists
		Stuple stuple;
		TChain* myChain;

		TFile *histFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		//HistManager myHistMan;
		TH1F *T1,*T2,*T3;
		TH1F *haloPhiWedge_1j, *haloEmTime_1j;
		TH1F *haloPhiWedge_2j, *haloEmTime_2j;
		TH1F *phoEt_iso25, *phoEt_50, *phoEt_70;
		TH1F *phoEmTime;

		Hist_t Hsignal, Hhalo, Hcosmic, Hqcd;
		Hist_t Hmc, HmcUp, HmcDown;		///for pho mc, ewk mc (CENTRAL/JES UP/DOWN)
		
		bool bMcFlag;			//must be set by user for different datasets, default is true, to make it safe

		std::vector<unsigned int> vHaloCount, vCosmicCount, vQCDCount, vSignalCount;		// 0=total,1=1jet,2=2jets
		std::vector<unsigned int> vMcCount, vMcUpCount, vMcDownCount;


		std::vector<EvtTag_t> EVT;
		float fJetEt;			// sets the jet et for Rays W+Z peack search
		TF1* tf1_sideband;
		
};
#endif
