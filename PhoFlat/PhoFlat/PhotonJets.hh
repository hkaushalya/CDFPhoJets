#ifndef PHOTONJETS_HH
#define PHOTONJETS_HH

/*{{{*/
/**********************************************************
 * $Id: PhotonJets.hh,v 1.18 2011/05/02 17:23:10 samantha Exp $
 * $Log: PhotonJets.hh,v $
 * Revision 1.18  2011/05/02 17:23:10  samantha
 * ADDED: To check for b-tag jets. This version regquires on
 * one b-tag jets. I have made plots with 2-btag jets. see elog:1792
 *
 * Revision 1.17  2011/04/26 23:33:42  samantha
 * Many spacing asjustments.
 * ADDED: ExoticEvent() method to search for fancy events. This method
 * needs to uncommented in cc file to be effective.
 *
 * Revision 1.16  2010/07/26 21:08:11  samantha
 * MAJOR CHANGES: Added all QCD ID systematic to here.
 * ADDED: 1. fMaxMet cut of 1500.GeV default.
 * 2. Few more plots to look at Vtx distribution with MET.
 * 3. FailStripCuts() to speedup things a bit.
 * 4. Many funtions to do debugging and to understand some events.
 * CheckObjOverlap(), CheckJetTowers(), FoundPartialJet().
 *
 * Revision 1.15  2010/04/13 16:11:17  samantha
 * REPLACED: CommonTools.hh which is used to get colors with IOColors.hh. I moved
 * those color function to IOColors.hh.
 * ADDED: Max njet 15 cut Set/Get methods.
 *
 * Revision 1.14  2010/04/07 17:44:39  samantha
 * MODIFIED: 1. Dropped the fabs() from SetMin/MaxJetEta(). I think this caused
 * some problems when I was applying cuts to jets in the JetList.
 * ADDED: 1. Option to enable/disable nvtx based weighing.
 * 2. The dataset used for the job can be passed in using SetDataset(). I added
 * this to get the correct Nvtx weight function for each of the MC samples from a
 * root file (this root file has weight functions for all MC samples). This can be
 * used to automcatically execute methods speicific for certain dataset.
 *
 * Revision 1.13  2010/03/31 23:14:29  samantha
 * ADDED: 1. fMaxJetEt and fMinJetEta.
 * 2. PhotonJets::GetJetEtaWeight() to generate nominal and error weights for
 * reweighing. GetHistBin() find the corresponding bin values for the jet eta, to get the
 * weights.
 * 3. bMetCleanup option to apply met cleanup cuts to see if that will improve the
 * MET plot. Also added these methods for this and they needs to be tested with new
 * V7 stuples which has full MET vector info.BadMet()
 *
 * Revision 1.12  2010/03/24 19:01:50  samantha
 * ADDED: 	1. Set/Get max photon et
 * 	2. DoHaloRejHists(), DoHaloRejCalc()  to get the halo rejection power
 * 	when run over data.;
 * 	3. GetCosmicEstAndErr() to get comsic estimated when run over data.
 *
 * 	for above calculation I had added additional sets of folder/hists toe
 * 	the job.
 * 	---
 * 	TH1F *haloPhiWedge_1j[2], *haloEmTime_1j[2];    //0=before 1= after BH cuts
 * 	TH1F *haloPhiWedge_2j[2], *haloEmTime_2j[2];
 * 	TH1F *haloEt_1j[2], *haloEt_2j[2];
 * 	TH1F *haloNvtx[2];
 *
 * 	Histograms::EventHists_t HhaloNoVtx_Evt_b4, HhaloNoVtx_Evt_a4;
 * 	Histograms::PhotonHists_t HhaloNoVtx_Pho_b4, HhaloNoVtx_Pho_a4;
 * 	Histograms::EventHists_t HhalosidebandNoVtx_Evt_b4, HhalosidebandNoVtx_Evt_a4;
 * 	Histograms::PhotonHists_t HhalosidebandNoVtx_Pho_b4, HhalosidebandNoVtx_Pho_a4;
 *
 * 	4. GetJetEtaWeight(const TF1* func, const float x) which derives an
 * 	addition weight for the sideband reweighing wbase on lead jet detector
 * 	eta. see elog#1590.
 *
 * Revision 1.11  2010/02/13 03:43:26  samantha
 * ADDED: MEt cut setting. Deafult is set to 0 in cc file.
 *
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
#include "TRandom.h"
#include <cmath>

class PhotonJets
{
	public:
		PhotonJets();
		int iProgressBy;
	
		enum PRINTLEVELS {
			iDEBUG=6,
			iINFO=5,
			iFATAL=4,
		};

		//const int iNCUTS=10;
		//enum FAILEDCUTS {
		//	iTrigger = 0;
		//	iGoodPhoton = 1;
		//	iBadMet = ;
		//};

		struct Hist_t {
			Histograms::EventHists_t p1j_Evt, p2j_Evt;				//
			Histograms::PhotonHists_t p1j_Pho, p2j_Pho;
			Histograms::JetHists_t p1j_Jet, p2j_Jet1, p2j_Jet2;
			Histograms::Photon1JetHists_t p1j_PhoJet, p2j_PhoJet1, p2j_PhoJet2;
			Histograms::Photon2JetsHists_t p2j_PhoJets;		//jet>=2 case, pho+2jets
			Histograms::TwoJetsHists_t p2j_TwoJets;
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
											const float fWeight=1, const bool debug=false,
											const bool bUseJetEtaWeights=false);
	
	void SelectPhoton(const CommonVars Vars, const PhotonList phos,
									const ElectronList eles, const JetList jets, 
									std::vector<int>& UsedPho, std::vector<int>& UsedEle);
	void SelectMCSidebandPhoton(const CommonVars Vars, const PhotonList phos,
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
		
		void SetMinPhotonEt(const float et) 
		{ 
			if (et <15.0)
			{
				std::cout << __FILE__ << ":" << __FUNCTION__ 
					<< red << "ERROR! Photon Et must be at least that was used to generate stuples (from version>7)\n"
					<< " This can cause jet list to have duplicate objects when adding additional EM objects\n"
					<< " back to jet list!\n"
					<< " Using default minimum value."
					<< clearatt << std::endl;
			} else 
			{	
				fMinPhotonEt = et; 
			}
		}
		void SetMaxPhotonEt(const float et) { fMaxPhotonEt = et; }
		void SetMinJetEt(const float et) { fMinJetEt = et; }
		void SetMaxJetEt(const float et) { fMaxJetEt = et; }
		void SetMaxPhotonEta(const float eta) { fMaxPhotonEta = fabs(eta); }
		void SetMinJetEta(const float eta) { fMinJetEta = eta; }
		void SetMaxJetEta(const float eta) { fMaxJetEta = eta; }
		void SetCosmicTimes(const float tmin, const float tmax) { fMinCosmicEmTime = tmin; fMaxCosmicEmTime = tmax; }
		void SetSignalTimes(const float tmin, const float tmax) { fMinSignalEmTime = tmin; fMaxSignalEmTime = tmax; }
		void SetMinClass12Vtx(const int nvt) { iMinClass12Vtx = nvt; }
		void SetMaxClass12Vtx(const int nvt) { iMaxClass12Vtx = nvt; }
		void SetHaloType(const int t) { iHaloType = t; }
		void SetMinNjet15(const int nj) { iMinNjet15 = nj; }
		void SetMaxNjet15(const int nj) { iMaxNjet15 = nj; }

		float GetMinPhotonEt() const { return fMinPhotonEt; }
		float GetMaxPhotonEt() const { return fMaxPhotonEt; }
		float GetMinJetEt() const { return fMinJetEt; }
		float GetMaxJetEt() const { return fMaxJetEt; }
		float GetMaxPhotonEta() const { return fMaxPhotonEta; }
		float GetMaxJetEta() const { return fMaxJetEta; }
		float GetMinJetEta() const { return fMinJetEta; }
		float GetCosmicTimeLoLimit() const { return fMinCosmicEmTime; }
		float GetCosmicTimeUpLimit() const { return fMaxCosmicEmTime; }
		float GetSignalTimeLoLimit() const { return fMinSignalEmTime; }
		float GetSignalTimeUpLimit() const { return fMaxSignalEmTime; }
		int GetMinClass12Vtx() const { return iMinClass12Vtx; }
		int GetMaxClass12Vtx() const { return iMaxClass12Vtx; }
		int GetHaloType() const { return iHaloType; }
		int GetMinNjet15() const { return iMinNjet15; }
		int GetMaxNjet15() const { return iMaxNjet15; }
		
		void ReweightSideband(const int b) { fRewgtSideband = b; }
		int GetReweightSideband() const { return fRewgtSideband; }
		void PrintJobSettings();
		float GetMinMet() const { return fMinMet; }
		void SetMinMet(const float met) { fMinMet = met; }
		float GetMaxMet() const { return fMaxMet; }
		void SetMaxMet(const float met) { fMaxMet = met; }

		void DoHaloRejHists(const Stuple& stuple);
		void DoHaloRejCalc();
		void DoHaloRejCalc(const TH1F *haloPhiWedge_b4, const TH1F *haloPhiWedge_a4, const Double_t Nid);
		void GetCosmicEstAndErr();
		void GetCosmicEstAndErr(const TH1 *hj1, const TH1 *hj2);
		double GetJetEtaWeight(const float phoEt, const float jetEta);
		void GetHistBin(const TH1* hist_JetEtaWghts, const float jetEta, float &val, float &err);
		void SetApplyMetCleanUpCuts(const bool t) { bMetCleanup = t;}	
		bool GetApplyMetCleanUpCuts() const { return bMetCleanup;}	
		bool BadMet(const CommonVars& Vars, const PhotonList& phos, 
							const JetList& jets, const int iPhoIndex, 
							const int iLeadJetIndex, TH1* phist, TH1* jhist); 
		void UseDATAMCNvtxWeights(const int i)
		{ 
			if (i<0 || i>2)
			{
				std::cout << __FILE__ << ":" << __FUNCTION__ << ":" 
							<< "Invalid value! pl. check!" << std::endl;
				assert(false);
			}
			iUseDataMcNvtxWgts = i;
			std::cout << "Setting from hh = " << iUseDataMcNvtxWgts << std::endl;
		}
		int GetUseDATAMCNvtxWeights() const { return iUseDataMcNvtxWgts;}
		float GetDATAMCNvtxWgt(const int Nvtx12);
		int GetDataset() const { return iDataset;}
		void  SetDataset(const int t) { iDataset = t;}
		bool FailStripCuts(const Stuple&);

		bool SidebandPhoPassTightHadEmCut(const float& fEc, const float& fHadEm);
		bool SidebandPhoPassTightIsoCut(const float& fEtc, const float& fIso);
		bool SidebandPhoPassTightTrkPtCut(const float& fEtc, const float& fTrkPt);
		bool SidebandPhoPassTightTrkIsoCut(const float& fEtc, const float& fTrkIso);
		void CheckObjOverlap(const CommonVars& cvars, const PhotonList& phos, const ElectronList& eles, const JetList& jets);
		void CheckJetTowers(const CommonVars& cvars, const JetList& jets);
		void SetPrintLevel(const int p) { iPrintLevel = p; }
		int PrintLevel() const { return iPrintLevel; }
		int FoundPartialJet(const CommonVars& cvars, const PhotonList& phos, const ElectronList& eles, const JetList& jets);

		bool ExoticEvent(const Stuple& st);
		void RequireBTagJets(const bool b) { bRequireTagJets = b; }		
		bool UseBTagJets() const { return bRequireTagJets; }		
			
	private:
		std::string sHistFileName;		//name of the root file to save the hists
		Stuple stuple;
		TChain* myChain;
		int iCurrTree;

		TFile *histFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		//HistManager myHistMan;
		TH1F *haloPhiWedge_1j[2], *haloEmTime_1j[2];	//0=before 1= after BH cuts
		TH1F *haloPhiWedge_2j[2], *haloEmTime_2j[2];
		TH1F *haloEt_1j[2], *haloEt_2j[2];
		TH1F *haloNvtx[2];
		//TH1F *myMet;
		Histograms::EventHists_t HhaloNoVtx_Evt_b4, HhaloNoVtx_Evt_a4;
		Histograms::PhotonHists_t HhaloNoVtx_Pho_b4, HhaloNoVtx_Pho_a4;
		Histograms::EventHists_t HhalosidebandNoVtx_Evt_b4, HhalosidebandNoVtx_Evt_a4;
		Histograms::PhotonHists_t HhalosidebandNoVtx_Pho_b4, HhalosidebandNoVtx_Pho_a4;

		
		TH1F *phoEt_iso25, *phoEt_50, *phoEt_70;

		Hist_t Hsignal, Hhalo, Hcosmic, Hqcd;
		Hist_t Hmc, HmcUp, HmcDown;		///for pho mc, ewk mc (CENTRAL/JES UP/DOWN)
		Hist_t Hsideband_ewk;			//tom reomove EWK component in Hcd
		Hist_t Hsideband_cosmic, Hsideband_halo;
		Hist_t HsidebandSyst_HadEm, HsidebandSyst_Iso, HsidebandSyst_TrkPt, HsidebandSyst_TrkIso;
		
		bool bMcFlag;			//must be set by user for different datasets, default is true, to make it safe

		std::vector<unsigned int> vHaloCount, vCosmicCount, vQCDCount, vSignalCount;		// 0=total,1=1jet,2=2jets
		std::vector<unsigned int> vMcCount, vMcUpCount, vMcDownCount, vMcSidebandEWKCount;
		std::vector<unsigned int> vSidebandCosmicCount, vSidebandHaloCount;
		std::vector<unsigned int> vSidebandSystHadEmCount, vSidebandSystIsoCount, vSidebandSystTrkPtCount,	vSidebandSystTrkIsoCount;


		std::vector<EvtTag_t> EVT;
		TF1* tf1_sideband_phoet;
		TF1* tf1_sideband_phoet_error;
		TF1* tf1_sideband_jeteta;
		TProfile* hSidebandJetEtaWeights;
		TProfile* hSidebandPhoEtWeights;

		bool bAddIsoToSidebandPho; // on/off addition of isolation energy to the sideband photon

		float fMinPhotonEt;
		float fMaxPhotonEt;
		float fMinJetEt;
		float fMaxJetEt;
		float fMaxPhotonEta;
		float fMinJetEta;
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
		int iMaxNjet15;
		int fRewgtSideband;
		float fMinMet, fMaxMet;
		//TCanvas *c;
		std::string sReadMeText;
		int iHaloCount[2];
		TRandom ran;
		TH1* hist_JetEtaWghts;
		bool bMetCleanup;
		int iUseDataMcNvtxWgts;
		TH1* hist_NvtxWghts;
		std::vector<float> vNvtxWeights;
		int iDataset;
		int iPrintLevel;
		unsigned int iRejEvts;
		bool bRequireTagJets;

};
#endif
