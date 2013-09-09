#ifndef PhotonTriggerStudy_HH
#define PhotonTriggerStudy_HH

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

class PhotonTriggerStudy
{
	public:
		PhotonTriggerStudy();
		int iProgressBy;
	
		void Main(TChain* chain, int iRunEvents = 0);
		void Init(TChain* chain);
		void DoMyStuff(Stuple& stuple);
		void CleanUp();
		void BookHistograms();
		std::string GetHistFileName() const { return sHistFileName; }
		void SetHistFileName(std::string name)  { sHistFileName = name; }

	void SetMcFlag(const int ii) { bMcFlag = (bool) ii; }
	bool GetMcFlag() const { return bMcFlag; }

		
	private:
		std::string sHistFileName;		//name of the root file to save the hists
		Stuple stuple;
		TChain* myChain;

		TFile *histFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		//HistManager myHistMan;
		TH1F *phoEt_iso25, *phoEt_50, *phoEt_70;
		TH1F *lphoEt_iso25, *lphoEt_50, *lphoEt_70;
		TH1F *sphoEt_iso25, *sphoEt_50, *sphoEt_70;
		TH1F *pphoEt_iso25, *pphoEt_50, *pphoEt_70;
		TH1F *plphoEt_iso25, *plphoEt_50, *plphoEt_70;
		TH1F *psphoEt_iso25, *psphoEt_50, *psphoEt_70;
		TH1F *cprWgtTight, *cprWgtLoose, *cprWgtSideband;
		TH1F *cprWgtTight_a, *cprWgtLoose_a, *cprWgtSideband_a;

		bool bMcFlag;			//must be set by user for different datasets, default is true, to make it safe

};
#endif
