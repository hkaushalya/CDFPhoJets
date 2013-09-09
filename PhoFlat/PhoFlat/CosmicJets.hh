#ifndef COSMICJETS_HH
#define COSMICJETS_HH

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

class CosmicJets
{
	public:
		CosmicJets();
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
			int SidebandCosmic;
			int SidebandHalo;
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
	void DoMCStuff(CommonVars Vars, PhotonList phos, ElectronList eles, JetList jets, Hist_t* hist, std::vector<unsigned int>& count);
	void GetCommonVariables(const Stuple stuple, CommonVars& commVars);
	int JetSelAndHistFill(const CommonVars Vars, const PhotonList phos, const ElectronList eles, 
											const JetList jets, Hist_t* hist,
											std::vector<int> vPhoIndex, std::vector<int> vEleIndex, 
											std::vector<unsigned int>& count, 
											const bool UseCprWgtForHt=false, 
											const float fWeight=1, const bool debug=false);
	void SelectPhoton(const CommonVars Vars, const PhotonList phos,
									const ElectronList eles, const JetList jets, 
									std::vector<int>& UsedPho, std::vector<int>& UsedEle);
	void SelectSidebandPhoton(const CommonVars Vars, const PhotonList phos,
									const ElectronList eles, const JetList jets, 
									std::vector<int>& UsedSidebandPho, std::vector<int>& UsedSidebandEle);

	//hist filling
		void FillEventHists(const CommonVars&, const PhotonList&, const ElectronList&, const JetList&, Histograms::EventHists_t* , 
								const int iPhoIndex, const float fWeight=1, const bool UseCprWgtForHt=false);
		void FillPhotonHists(const CommonVars&, const PhotonList&, Histograms::PhotonHists_t*,
								const int iPhoIndex, const float fWeight=1, const bool debug=false);									
		void FillElectronHists(const ElectronList&, Histograms::PhotonHists_t* ,
								const int iEleIndex, const float fWeight=1);
		void FillJetHists(const CommonVars&, const JetList&, Histograms::JetHists_t* , 
								const int iJetIndex,	const float fWeight=1);
		void FillPhoton1JetHists(const CommonVars&, const PhotonList&, const JetList&, Histograms::Photon1JetHists_t* ,
								const int iPhoIndex, const int iJetIndex,
								const float fWeight=1);
		void FillElectron1JetHists(const ElectronList&, const JetList&, Histograms::Photon1JetHists_t* ,
								const int iEleIndex, const int iJetIndex,
								const float fWeight=1);
		void FillPhoton2JetsHists(const CommonVars&, const PhotonList&, const JetList&, Histograms::Photon2JetsHists_t* hist,
								const int iPhoIndex, const int iJet1Index,
								const int iJet2Index, const float fWeight=1);
		void FillElectron2JetsHists(const ElectronList&, const JetList&, Histograms::Photon2JetsHists_t* hist,
								const int iEleIndex, const int iJet1Index,
								const int iJet2Index, const float fWeight=1);
		void FillTwoJetsHists(const CommonVars&, const JetList&, Histograms::TwoJetsHists_t* hist,
								const int iJet1Index, const int iJet2Index,
								const float fWeight=1);
		int Prepass(const JetList jets);


		void PrintHeader(const CommonVars&) const;
		void SetReportProgress(const unsigned int l) { iProgressBy = l; }
		void RemoveDuplicateEMobjects(PhotonList& phos, ElectronList& eles);
		void SetAddIsoEtoSidebandPho(const bool b) { bAddIsoToSidebandPho = b; }
		bool GetAddIsoEtoSidebandPho() const { return bAddIsoToSidebandPho; }
		
		void SetMinPhotonEt(const float et) { fMinPhotonEt = et; }
		void SetMinJetEt(const float et) { fMinJetEt = et; }
		void SetMaxPhotonEta(const float eta) { fMaxPhotonEta = eta; }
		void SetMaxJetEta(const float eta) { fMaxJetEta = eta; }
		void SetCosmicTimes(const float tmin, const float tmax) { fMinCosmicEmTime = tmin; fMaxCosmicEmTime = tmax; }
		void SetSignalTimes(const float tmin, const float tmax) { fMinSignalEmTime = tmin; fMaxSignalEmTime = tmax; }
		void SetMinClass12Vtx(const int nvt) { iMinClass12Vtx = nvt; }
		void SetMaxClass12Vtx(const int nvt) { iMaxClass12Vtx = nvt; }
		void SetHaloType(const int t) { iHaloType = t; }
		void SetMinNjet15(const int nj) { iMinNjet15 = nj; }

		float GetMinPhotonEt() const { return fMinPhotonEt; }
		float GetMinJetEt() const { return fMinJetEt; }
		float GetMaxPhotonEta() const { return fMaxPhotonEta; }
		float GetMaxJetEta() const { return fMaxJetEta; }
		float GetCosmicTimeLoLimit() const { return fMinCosmicEmTime; }
		float GetCosmicTimeUpLimit() const { return fMaxCosmicEmTime; }
		float GetSignalTimeLoLimit() const { return fMinSignalEmTime; }
		float GetSignalTimeUpLimit() const { return fMaxSignalEmTime; }
		int GetMinClass12Vtx() const { return iMinClass12Vtx; }
		int GetMaxClass12Vtx() const { return iMaxClass12Vtx; }
		int GetHaloType() const { return iHaloType; }
		int GetMinNjet15() const { return iMinNjet15; }
		
		void ReweightSideband(const bool b) { fRewgtSideband = b; }
		bool GetReweightSideband() const { return fRewgtSideband; }
		void PrintJobSettings();

		
	private:
		std::string sHistFileName;		//name of the root file to save the hists
		Stuple stuple;
		TChain* myChain;

		TFile *histFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		TH1F *phoEt_iso25, *phoEt_50, *phoEt_70;

		Hist_t Hsignal, Hhalo, Hcosmic, Hqcd;
		Hist_t Hmc, HmcUp, HmcDown;		///for pho mc, ewk mc (CENTRAL/JES UP/DOWN)
		Hist_t Hsideband_ewk;			//tom reomove EWK component in Hcd
		Hist_t Hsideband_cosmic, Hsideband_halo;
		
		bool bMcFlag;			//must be set by user for different datasets, default is true, to make it safe

		std::vector<unsigned int> vHaloCount, vCosmicCount, vQCDCount, vSignalCount;		// 0=total,1=1jet,2=2jets
		std::vector<unsigned int> vMcCount, vMcUpCount, vMcDownCount, vMcSidebandEWKCount;
		std::vector<unsigned int> vSidebandCosmicCount, vSidebandHaloCount;


		std::vector<EvtTag_t> EVT;
		float fJetEt;			// sets the jet et for Rays W+Z peack search
		TF1* tf1_sideband;
		TProfile* hist_tf1_wghts;

		Clock *clk;
		bool bAddIsoToSidebandPho; // on/off addition of isolation energy to the sideband photon

		float fMinPhotonEt;
		float fMinJetEt;
		float fMaxPhotonEta;
		float fMaxJetEta;
		float fMinSignalEmTime;
		float fMaxSignalEmTime;
		float fMinCosmicEmTime;
		float fMaxCosmicEmTime;
		int iMinClass12Vtx;
		int iMaxClass12Vtx;
		float fMaxVtxz;
		int iHaloType;
		int iMinNjet15;
		bool fRewgtSideband;


};
#endif
