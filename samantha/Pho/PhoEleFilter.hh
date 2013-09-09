#ifndef PhoEleFilter_HH
#define PhoEleFilter_HH

///////////////////////////////////////////////////////////
// See sourcce file for class description                //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#if !defined (__CINT__) || defined (__MAKECINT__)

#endif
#include "Stntuple/loop/TStnModule.hh"
#include <string>
#include <iostream>
#include "samantha/Pho/InitSuperPhotons.hh"

class PhoEleFilter: public TStnModule {

	protected:

	public:
		PhoEleFilter(const char* name="PhoEleFilter", const char* title="PhoEleFilter");
		~PhoEleFilter();

				/************** abstract objects *************************/
		struct Global_Counters_t {
			unsigned int evtsProcessed;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
		};




		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void SaveHistograms();

		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }

		float GetMinEmEt() const { return fMinEmEt; }
		float GetMaxEmDetEta() const { return fMaxEmDetEta; }
		void SetMinEmEt(const float et) { fMinEmEt = fabs(et); }
		void SetMaxEmDetEta(const float eta) {  fMaxEmDetEta = fabs(eta); }
				
	private:
		Global_Counters_t counter;
		bool bRunPermit;


		InitSuperPhotons* initSpMod;

		unsigned int phoOnlyEvts;
		unsigned int phoEleEvts;
		unsigned EleOnlyEvts;
		unsigned int cemEle, pemEle, cemPho, extraPho; 
		bool bNoSummary;			//control the end job summary print

		float fMinEmEt;
		float fMaxEmDetEta;

	ClassDef(PhoEleFilter,1)
};

#endif
