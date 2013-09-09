#ifndef PHOENIXPHOSTUDY_HH
#define PHOENIXPHOSTUDY_HH

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
#include <string>

class PhoenixPhoStudy
{
	public:
		PhoenixPhoStudy();
		int iProgressBy;
	

		struct Hist_t {
			Histograms::EventHists_t p1j_Evt, p2j_Evt;				//
			Histograms::PhotonHists_t p1j_Pho, p2j_Pho;
			Histograms::JetHists_t p1j_Jet, p2j_Jet1, p2j_Jet2;
			Histograms::Photon1JetHists_t p1j_PhoJet, p2j_PhoJet1, p2j_PhoJet2;
			Histograms::Photon2JetsHists_t p2j_PhoJets;		//jet>=2 case, pho+2jets
			Histograms::TwoJetsHists_t p2j_TwoJets;
			TH1F* z1j_zpt;
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
			TH1F* z2j_zj1j2_mass;

		};

		struct EvtTag_t {
			unsigned int Run;
			unsigned int Evt;
			unsigned int Halo;
			unsigned int Cosmic;
			unsigned int Signal;
			unsigned int QCD;
			unsigned int Phx;
			
		};

		void Main(TChain* chain, int iRunEvents = 0);
		void Init(TChain* chain);
		void CleanUp();
		void BookHistograms();
		std::string GetHistFileName() const { return sHistFileName; }
		void SetHistFileName(std::string name)  { sHistFileName = name; }


	void SetMcFlag(const int ii) { bMcFlag = (bool) ii; }
	bool GetMcFlag() const { return bMcFlag; }
	void DoMCStuff(CommonVars Vars, PhotonList phos, ElectronList eles, JetList jets, Hist_t* hist, std::vector<unsigned int>& count);

	//hist filling
		void FillEventHists(Stuple&, Histograms::EventHists_t* , 
								const float fWeight=1);
		void FillPhotonHists(Stuple&, Histograms::PhotonHists_t*,
								const int iPhoIndex, const float fWeight=1);									
		void FillElectronHists(const ElectronList&, Histograms::PhotonHists_t* ,
								const int iEleIndex, const float fWeight=1);

		void DoMyStuff(Stuple st);

		bool GetPhoenixRejectStatus() const { return bRejectPhoenix; } 
		void RejectPhoenix(int ss) 
		{ 
			bRejectPhoenix = (bool) ss;
			if (bRejectPhoenix) sPhxStatus = " After Phoenix Rejection ";
			else sPhxStatus = " Before Phoenix Rejection ";
		} 
	
		bool CheckOrder(Stuple);


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

		Hist_t Hsignal,Hqcd, Hphx;
		
		bool bMcFlag;			//must be set by user for different datasets, default is true, to make it safe

		std::vector<unsigned int> vQCDCount, vSignalCount, vPhxCount;		// 0=total,1=1jet,2=2jets

		std::vector<EvtTag_t> EVT;
		bool bRejectPhoenix; 
		std::string sPhxStatus;
		
};
#endif
