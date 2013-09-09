#ifndef DUMPSTNTUPLEEVENT_HH
#define DUMPSTNTUPLEEVENT_HH

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include <iostream>
#include <Stntuple/obj/TStnJetBlock.hh>
#include <Stntuple/obj/TCalDataBlock.hh>
#include <Stntuple/obj/TStnTrackBlock.hh>
#include <Stntuple/obj/TStnPhotonBlock.hh>
#include <Stntuple/obj/TStnElectronBlock.hh>
#include <Stntuple/data/TCalTower.hh>
#include "samantha/Pho/InitSuperPhotons.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"


class DumpStntupleEvent: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnTrackBlock*     fTrackBlock;
		TStnJetBlock*		  fJetBlockClu04;
		TCalDataBlock*      fCalDataBlock;
		TStnPhotonBlock*    fPhotonBlock;
		TStnElectronBlock*  fElectronBlock;


	public:
		DumpStntupleEvent(const char* name="DumpStntupleEvent", const char* title="DumpStntupleEvent");
		~DumpStntupleEvent();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock() { return fHeaderBlock;  }
		TStnTrackBlock*    GetTrackBlock()	{ return fTrackBlock; 	}
		TStnJetBlock*		 GetJetBlock()    { return fJetBlockClu04;}
		TCalDataBlock*     GetCalDataBlock() { return fCalDataBlock;}


		typedef std::vector<TCalTower*> CalDataArray; // need for beam halo
		typedef std::vector<TCalTower*>::iterator CalDataArrayI; // need for beam halo

				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsProcessed;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
		};


		// all the histograms


		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void FillDataBlocks(int ientry);
		void DumpJetTowers(int jet_ind); // for jets
		void DumpPhotonTowers(int pho_ind); // for photons
	  	void DumpElectronTowers(int ele_ind); // for electrons

		void DumpJetMatchStuff();

	private:
		Global_Counters_t counter;
		InitSuperPhotons* initSpMod;
		JetFilterModuleV2* jetMod;
		bool bRunPermit;

	ClassDef(DumpStntupleEvent,1)
};

#endif
