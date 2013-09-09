#ifndef PLOTEMTIMES_HH
#define PLOTEMTIMES_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include <iostream>
#include <string>
#include "Stntuple/loop/TStnModule.hh"
#include "samantha/obj/SuperPhoton.hh"
#include <Stntuple/obj/TStnEmTimingBlock.hh>
#include "Stntuple/ana/TEmTimingModule.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "CalorGeometry/CalConstants.hh"
#include "Stntuple/data/TCalTower.hh"
#include <Stntuple/obj/TStnPhotonBlock.hh>
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include "samantha/Pho/EmTimeCorrection.hh"
#include "Stntuple/obj/TStnMetBlock.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"



class PlotEmTimes: public TStnModule {

	protected:
		TStnHeaderBlock*    fHeaderBlock;
		TStnEmTimingBlock*  fTimingBlock;
		TCalDataBlock*      fCalBlock;
		TStnPhotonBlock*    fPhotonBlock;
		TStnMetBlock*       fMetBlock;

		typedef std::vector<TCalTower*> CalDataArray;
		typedef std::vector<TCalTower*>::iterator CalDataArrayI;

	public:
	
		PlotEmTimes(const char* name="PlotEmTimes", const char* title="PlotEmTimes");
		~PlotEmTimes();
	
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnPhotonBlock*   GetPhotonBlock()  { return fPhotonBlock;  }
		TCalDataBlock* 	 GetCalBlock()  { return fCalBlock;  }
		TStnMetBlock*      GetMetBlock() 	{ return fMetBlock; }
		
		/*********** define all abstract data types here **********/
		struct GlobalCounters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
		};	

		struct Hists_t {
			TProfile* PhoEmTimeVsRun_b4;
			TProfile* TPhoEmTimeVsRun_b4;
			TProfile* TPho1JetEmTimeVsRun_b4;
			TProfile* TPho2JetEmTimeVsRun_b4;
			TProfile* PhoEmTimeVsRun_a4;
			TProfile* TPhoEmTimeVsRun_a4;
			TProfile* TPho1JetEmTimeVsRun_a4;
			TProfile* TPho2JetEmTimeVsRun_a4;
			TH1F* PhoEmTime_b4;
			TH1F* TPhoEmTime_b4;
			TH1F* TPho1JetEmTime_b4;
			TH1F* TPho2JetEmTime_b4;
			TH1F* PhoEmTime_a4;
			TH1F* TPhoEmTime_a4;
			TH1F* TPho1JetEmTime_a4;
			TH1F* TPho2JetEmTime_a4;
		};




		/*******************END abstract data types ***************/

		
		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();
		
		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
	
		double  GetEmTiming(const int phoInd);
		float  GetHadTdcTime(const int phoInd);
		void FillDataBlocks(int ientry);
		void MatchCalorTowers(int pho_ind, CalDataArray* cdholder);

		void SetEmTimeCorrectionFile(std::string file_) { sTimeCorrectionFile = file_; }
		std::string GetEmTimeCorrectionFileName() const  { return sTimeCorrectionFile; }
		
	private:
	
		GlobalCounters_t counter;
		EmTimeCorrection *cemT0Runs;
		InitSuperPhotons* initSpMod;
		float fCorrection;
		bool bRunPermit;	
		std::string sTimeCorrectionFile;
		JetFilterModuleV2* jetMod;
		Hists_t hist;

	ClassDef(PlotEmTimes,0)
};

		//bool SortTowersByEnergy(TCalTower* c1, TCalTower* c2);
#endif
