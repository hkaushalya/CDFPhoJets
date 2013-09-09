#ifndef TAGCOSMICPHOTONS_HH
#define TAGCOSMICPHOTONS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#include <iostream>
#include <string>
#include "Stntuple/loop/TStnModule.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/SetEmTimes.hh"
#include "samantha/obj/Histograms.hh"

#endif


class TagCosmicPhotons: public TStnModule {

	protected:
		
	public:
	
		TagCosmicPhotons(const char* name="TagCosmicPhotons", const char* title="TagCosmicPhotons");
		~TagCosmicPhotons();
	
		/*********** define all abstract data types here **********/
		/*******************END abstract data types ***************/

		
		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();
		
		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
	
		void SetCosmicTimeWindow(float Tmin, float Tmax) {
			fCosmicTimeMin = Tmin;
			fCosmicTimeMax = Tmax;
		}
		
		void FillDataBlocks(int ientry);
		
		
		void  GetCosmicTimeWindow(float &tmin_, float &tmax_) const {
			tmin_ = fCosmicTimeMin;
			tmax_ = fCosmicTimeMax;
		}
		float GetCosmicTimeMin() const { return fCosmicTimeMin; }
		float GetCosmicTimeMax() const { return fCosmicTimeMax; }
	
	private:
	
		float fCosmicTimeMin;
		float fCosmicTimeMax;
	
		InitSuperPhotons* initSpMod;
		SetEmTimes* setEmTimeMod;
		bool bRunPermit;				//make sure all depndencies are met before running.

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;
		
	ClassDef(TagCosmicPhotons,0)
};

		//bool SortTowersByEnergy(TCalTower* c1, TCalTower* c2);
#endif
