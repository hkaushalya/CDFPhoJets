#ifndef SKIMCOSMICS_HH
#define SKIMCOSMICS_HH

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
#include "Stntuple/loop/TStnAna.hh"


class SkimCosmics: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;


	public:
		SkimCosmics(const char* name="SkimCosmics", const char* title="SkimCosmics");
		~SkimCosmics();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }


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

		void SetCosmicTimeWindow(float Tmin, float Tmax) {
			fCosmicTimeMin = Tmin;
			fCosmicTimeMax = Tmax;
		}
		void  GetCosmicTimeWindow(float &tmin_, float &tmax_) const {
			tmin_ = fCosmicTimeMin;
			tmax_ = fCosmicTimeMax;
		}
		float GetCosmicTimeMin() const { return fCosmicTimeMin; }
		float GetCosmicTimeMax() const { return fCosmicTimeMax; }
		void FillDataBlocks(int ientry);

	private:
		bool 	qMc;  				// true if MonteCarlo data
		bool bPermitRun;			// make sure all pre-req are met before run
		float fCosmicTimeMin;
		float fCosmicTimeMax;

		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		TagTightPhotons* tightMod;
		

		unsigned int iEvtsProcessed, iEvtsPassed, iCosmicEvts;
		
	ClassDef(SkimCosmics,1)
};

#endif
