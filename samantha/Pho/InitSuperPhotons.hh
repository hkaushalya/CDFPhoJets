#ifndef INITSUPERPHOTONS_HH
#define INITSUPERPHOTONS_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <vector>
#include <string>
#include "TLorentzVector.h"
#include <Stntuple/loop/TStnModule.hh>
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include <Stntuple/obj/TStnPhotonBlock.hh>
#include <Stntuple/obj/TStnPhoton.hh>
#include <Stntuple/obj/TStnVertexBlock.hh>
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "Stntuple/obj/TStnElectronBlock.hh"
#endif

class InitSuperPhotons: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnPhotonBlock*    fPhotonBlock;
		TStnVertexBlock*    fVertexBlock;
		TStnElectronBlock*  fElectronBlock;


	public:
		InitSuperPhotons(const char* name="InitSuperPhotons", const char* title="InitSuperPhotons");
		~InitSuperPhotons();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TStnPhotonBlock*   GetPhotonBlock()  { return fPhotonBlock;  }
		TStnVertexBlock*   GetVertexBlock()  { return fVertexBlock;  }
		TStnElectronBlock* GetElectronBlock(){ return fElectronBlock;}


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
		};


		// all the histograms

		struct GeneralPlots_t {
			TH1F* RunNumbers;		// run number 
			TH1F* EvtsRun;			// number of evts ran through
			TH1F* EvtsPass;		// number of events pass this module
		}; //GeneralPlots_t

		struct PhotonPlots_t {
			TH1F* Detector;
			TH1F* DetEta;
			TH1F* DetPhi;
			TH1F* EtCorr;
			TH1F* XCes;
			TH1F* ZCes;
			TH1F* HadEm;
			TH1F* Chi2Mean;
			TH1F* N3d;
			TH1F* Iso4;
			TH1F* TrkPt;
			TH1F* TrkIso;
			TH1F* Ces2Wire;
			TH1F* Ces2Strip;
			TH1F* HadTDCtime;
		};



		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void 	BookHistograms();
		void 	SaveHistograms();
		void 	Cleanup();

		//my methods
		void 	FillDataBlocks(int ientry);  // put all GetEntry statements here
		void 	CreateSuperPhotonVector();
		int   Class12Vertices();
		void  FillSuperPhotonHist();
		void 	BookPhotonHistograms(PhotonPlots_t&, std::string folder_name,
								std::string folder_title);
		void  BookGeneralHistograms(GeneralPlots_t& hist,
						std::string fold_name, std::string fold_title);

		// accessors for my 2nd mod
		int GetSuperPhoSize()  { return superpho.size(); }
		SuperPhoton* GetSuperPhoton(int i) { 
			SuperPhoton* sp = NULL;
			if ( i >= 0 && i < (int)superpho.size()) sp = &superpho[i];
			else StdOut(__FILE__,__LINE__,4,"Request for a non-existing SuperPhoton. Returning NULL!");
			return sp;
				
		}
		
		void DumpSuperPhotons() {
			StdOut("********* Event Number:",fHeaderBlock->EventNumber(),GetName());
			StdOut("NSuperPhotons",superpho.size(),GetName());
			for (itSp=superpho.begin(); itSp != superpho.end(); itSp++)
				itSp->Print();
		};
	

		int GetElectronBlockIndex(TStnPhoton* pho);		// To get the electron Block index of the electrons, needed  by JetFilterModule - 05-13-2008

		void SetPrintLevel(const int pp) { iPrintLevel = pp; }
		int GetPrintLevel() const { return iPrintLevel; }
	
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
		
	private:
		Global_Counters_t counter;
		
		//my super photons vector
		std::vector<SuperPhoton> superpho;
		std::vector<SuperPhoton>::iterator itSp;

		// histograms
		PhotonPlots_t hPhoton;
		GeneralPlots_t hGeneral;
		int iPrintLevel;
		bool bNoSummary;			//control the end job summary print

	ClassDef(InitSuperPhotons,1)
};

#endif
