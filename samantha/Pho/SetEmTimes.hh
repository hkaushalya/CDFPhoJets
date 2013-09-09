#ifndef SETEMTIMES_HH
#define SETEMTIMES_HH

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
#include "samantha/Pho/TagPMTSpikes.hh"
#include "samantha/utils/FreeFunctions.hh"



class SetEmTimes: public TStnModule {

	protected:
		TStnHeaderBlock*    fHeaderBlock;
		TStnEmTimingBlock*  fTimingBlock;
		TCalDataBlock*      fCalBlock;
		TStnPhotonBlock*    fPhotonBlock;
		TStnMetBlock*       fMetBlock;

		typedef std::vector<TCalTower*> CalDataArray;
		typedef std::vector<TCalTower*>::iterator CalDataArrayI;

	public:
	
		SetEmTimes(const char* name="SetEmTimes", const char* title="SetEmTimes");
		~SetEmTimes();
	
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnPhotonBlock*   GetPhotonBlock()  { return fPhotonBlock;  }
		TCalDataBlock* 	 GetCalBlock()  { return fCalBlock;  }
		TStnMetBlock*      GetMetBlock() 	{ return fMetBlock; }
		
		/*********** define all abstract data types here **********/
		struct GlobalCounters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
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
	
		double  GetEmTiming(const int phoInd, const float phoDetEta);
		float  GetHadTdcTime(const int phoInd);
		void FillDataBlocks(int ientry);
		void MatchCalorTowers(int pho_ind, CalDataArray* cdholder);

		void SetEmTimeCorrectionFile(std::string file_) { sTimeCorrectionFile = file_; }
		std::string GetEmTimeCorrectionFileName() const  { return sTimeCorrectionFile; }
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
		
	private:
	
		GlobalCounters_t counter;
		EmTimeCorrection *cemT0Runs;
		InitSuperPhotons* initSpMod;
		float fCorrection;
		bool bRunPermit;	
		std::string sTimeCorrectionFile;
		bool bNoSummary;			//control the end job summary print

	ClassDef(SetEmTimes,0)
};

		//bool SortTowersByEnergy(TCalTower* c1, TCalTower* c2);
#endif
