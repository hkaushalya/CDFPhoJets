#ifndef TAGCONVELECTRONS_HH
#define TAGCONVELECTRONS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnTrack.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnElectron.hh"
#include "Stntuple/obj/TStnElectronBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "TFolder.h"
#endif

class TagConvElectrons: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnElectronBlock*  fElectronBlock;
		TStnTrackBlock*     fTrackBlock;


	public:
		TagConvElectrons(const char* name="TagConvElectrons", const char* title="TagConvElectrons");
		~TagConvElectrons();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}
		TStnTrackBlock*    GetTrackBlock()   { return fTrackBlock;   }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int convEleFound;
		};


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
		void FillDataBlocks(int ientry);
		TFolder* GetHistoFolder(char *name, char* title);
		double 	ConversionSep(TStnTrack* t0, TStnTrack* t1, double & rconv);
		bool 	IsConversionElectron(TStnElectron* ele,
				     double& minsep, double& mindt, double& radius);

		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		InitSuperPhotons* InitSpMod;
		bool bNoSummary;			//control the end job summary print
	ClassDef(TagConvElectrons,1)
};

#endif
