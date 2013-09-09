#ifndef JERSTUDY_HH
#define JERSTUDY_HH

// $Log: JERStudy.hh,v $
// Revision 1.5  2009/03/29 04:24:36  samantha
// Another test for CVS to include comments in the source file. Let see if this also works. See Elog#1120
//

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif


#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/Pho/TagTightPhotons.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/obj/Histograms.hh"
#include "samantha/Pho/EventProperties.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/ana/TEmTimingModule.hh"
#include <Stntuple/obj/TStnEmTimingBlock.hh>
#include "TH1.h"
#include "TH2.h"


class JERStudy: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TCalDataBlock*      fCalBlock;
		TStnJetBlock*       fJetBlock;

		typedef std::vector<TCalTower*> CalDataArray;
		typedef std::vector<TCalTower*>::iterator CalDataArrayI;

	public:
		JERStudy(const char* name="JERStudy", const char* title="JERStudy");
		~JERStudy();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnJetBlock*      GetJetBlock() 	{  return fJetBlock; }
		TCalDataBlock*  GetCalDataBlock() 	{  return fCalBlock; }


				/************** abstract objects *************************/



		// all the histograms


		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
		void Cleanup();
		void DoJetSelection(int iSpIndex);

		void FillDataBlocks(int ientry);
		double  GetEmTiming(const int jetInd);
		void MatchCalorTowers(int jet_ind, CalDataArray* cdholder);

	private:
		bool 	qMc;  				// true if MonteCarlo data

		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TagTightPhotons* tightMod;
		EventProperties* evtPropMod;
		
		bool bPermitRun;			// make sure all pre-req are met before run

		Histograms::EventHists_t h1_Evt, h2_Evt;						// properties of events(met, sumet etc)
		Histograms::JetHists_t h1_Jet, h2_Jet1, h2_Jet2;

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		TH1F *hJetEmTime;
		TH1F *hLev6RawJetEtRatio; 				//
		TH2F *hL6RawJetEtRatioVsRawJetEt;
		
	ClassDef(JERStudy,1)
};

#endif
