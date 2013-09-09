#ifndef PHOELESTUDY_HH
#define PHOELESTUDY_HH

/**********************************************************
 * $Id: PhoEleStudy.hh,v 1.1 2013/02/28 03:37:57 samantha Exp $
 * $Log: PhoEleStudy.hh,v $
 * Revision 1.1  2013/02/28 03:37:57  samantha
 * Final commit. no checks.! these were never commited to cvs before!
 *
 **********************************************************
 *Samantha K. Hewamanage
 *********************************************************/

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

class PhoEleStudy
{
	public:
		PhoEleStudy();
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
	void CreatePhotonLists(const Stuple st, PhotonList& pholist, const int collection);
	void CreateElectronLists(const Stuple st, ElectronList& elelist, const int collection);
	void CreateJetLists(const Stuple st, JetList& jetlist, const int collection);
	void GetCommonVariables(const Stuple stuple, CommonVars& commVars);

		void PrintHeader(const CommonVars&) const;
		void SetReportProgress(const unsigned int l) { iProgressBy = l; }
		void RemoveDuplicateEMobjects(PhotonList& phos, ElectronList& eles);
		
			
	private:
		std::string sHistFileName;		//name of the root file to save the hists
		Stuple stuple;
		TChain* myChain;
		int iCurrTree;

		TFile *histFile;
		TDirectory* topDir;
		ReadInAll *readIn;
		
		bool bMcFlag;			//must be set by user for different datasets, default is true, to make it safe

		std::vector<unsigned int> vSignalCount;		// 0=total,1=1jet,2=2jets

};
#endif
