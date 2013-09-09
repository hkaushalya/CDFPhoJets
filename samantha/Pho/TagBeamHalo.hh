#ifndef TAGBEAMHALO_HH
#define TAGBEAMHALO_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "CalorGeometry/CalConstants.hh"
#include "Stntuple/data/TCalTower.hh"
#include <Stntuple/obj/TStnPhotonBlock.hh>
#include "Stntuple/obj/TStnMetBlock.hh"
#include "Stntuple/obj/TStnVertexBlock.hh"
#include "samantha/obj/SuperPhoton.hh"
#include "TFolder.h"
#include <string>
#include "samantha/utils/FreeFunctions.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include "samantha/obj/Histograms.hh"
#endif

class TagBeamHalo: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TCalDataBlock*      fCalBlock;
		TStnPhotonBlock*    fPhotonBlock;
		TStnVertexBlock*    fVertexBlock;
		TStnMetBlock*       fMetBlock;

		typedef std::vector<TCalTower*> CalDataArray;  // need for beam halo
		typedef std::vector<TCalTower*>::iterator CalDataArrayI;   // need for BH

	public:
		TagBeamHalo(const char* name="TagBeamHalo", const char* title="TagBeamHalo");
		~TagBeamHalo();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock;  }
		TCalDataBlock*      GetCalDataBlock()  { return fCalBlock; }
		TStnPhotonBlock*   GetPhotonBlock()  { return fPhotonBlock;  }
		TStnVertexBlock*   GetVertexBlock()  { return fVertexBlock;  }
		TStnMetBlock*      GetMetBlock() 	{ return fMetBlock; }

		int GetMode() { return iMode; }

				/************** abstract objects *************************/


		struct GlobalCounters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int haloCandidatesFoundByCut;		//__ identified by the NHad and Seed wedge cuts
		};


	  //-------------------------- parameters for Beam Halo
		struct HaloStuff {
			int seedwedge;
			int sidewedge;
			int seedwedgeH;
			int sidewedgeH;
			int eastNhad;
			int westNhad;
			double emseedE;
			double emsideE;
			double hadseedE;
			double seedwedgeHadE;
		};

		//______________________________ Beam Halo study histograms
		struct BeamHaloStudyHisto_t {

			//_____________________ My params
			TH1F* fHaloSeedWedgeH;   // number of seed wedge towers with Had>thresh 
			TH1F* fHaloSideWedgeH;   // number of side wedge towers with Em>thresh
			//______________________ Max's params
			TH1F* fHaloSeedWedge;   // number of seed wedge towers with Em>thresh 
			TH1F* fHaloSideWedge;   // number of side wedge towers with Em>thresh
			TH1F* fHaloEastNHad;    // number of had towers on East side
			TH1F* fHaloWestNHad;    // number of had towers on West side
			TH1F* fHaloEmSeedE;     // EM energy of seed wedge towers
			TH1F* fHaloEmSideE;     // EM energy of side wedge towers
			TH1F* fHaloHadSeedE;    // Plug HAD energy of seed wegde towers
			TH1F* fHaloSeedWedgeHadE; // HAD energy of seed wegde towers
			TH1F* fHaloNHad;        // number of towers on East & West: EastNHad+WestNHad

			TH2F* fHaloSeedWedgeEM_SeedWedgeHAD; // number of EM seed wedge towers vs. number of HAD towers 

			//______________________ Ray's params
			TH1F* fHaloEast;        // output of Pho->HaloEast()
			TH1F* fHaloWest;        // output of Pho->HaloWest()
			TH1F* fHaloEastWest;    // Pho->HaloEast()+Pho->HaloWest()
			TH1F* fHaloSideWedgeR;  // number of towers in side wedges according to Ray's code  
			TH1F* fHaloSeedWedgeR;  // number of towers in seed wedge according to Ray's code  
			TH1F* fHaloTwrInRawR;   // number of continues towers in seed wedge according to Ray's code  

			//______________________ my params
			TH1F* fHaloCesStripEoverE; // CesStripE/E(pho)
			TH1F* fHaloCesWireEoverE;  // CesWireE/E(pho)
			TH1F* fHaloHadTDC;         // Had TDC timing   

			TH1F* fHaloPhiWedge;  // wedge number of halo candidate
			TH1F* fHaloEta;       // eta of halo candidate
			TH1F* fHaloIso;       // iso of halo candidate
			TH1F* fHaloHadEm;     // Had/Em of halo candidate
			TH1F* fHaloCesChi2;   // CES Chi2 of halo candidate
			TH1F* fHaloEt;        // Et of halo candidate

			//my plots
			TH1F*  fHaloEt_to_Met;		//Beam halo candidate Et to Event Met ratio
			TH1F*  fHalo_Met_delPhi;		//Phi seperation between BH photon and MEt
			TH1F*  fHalo_Met_InvMass;		//Inv. mass of BH and MEt
			TH1F*  fnHalo_to_nPho;			// # of BH to photon ratio
			TH1F*	 fHalo_nPho_to_nJets;	// number of Photons to Jets in BH events
			TH2F*  fSeedWedge_HadTowers; //seed towers vs N Had towers
			TH1F*  fHaloMet;					//Beam halo MEt
			TH1F*  fHaloSumet;					//Beam halo MEt

		};


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
		
		void  DoMyBeamHalo(int pho_ind);
		void  getBeamHaloInfo(int pho_ind, int& seedWedge, int& seedWedgeH, 
							int& sideWedge, int& sideWedgeH, int& eastNHads,
							int& westNHads, double& emSeedE, double& emSideE, 
							double& hadSeedE, double& seedWedgeHadE);
		int   BeamHaloCut(HaloStuff *beamhalo);
		void  MatchCalorTowers(int pho_ind, CalDataArray* cdholder);
		int   Class12Vertices();
		void  FillBeamHaloStudyHistograms(BeamHaloStudyHisto_t& Hist, 
							 HaloStuff *beamhalo, TStnPhoton* Pho);

		void  BookBeamHaloStudyHistograms(BeamHaloStudyHisto_t& Hist, std::string sFoldName,
								std::string sFoldTitle);
		int IsHalo(SuperPhoton*);

		//setters
		void SetMode(int mode=0) {
			if (mode >=0 && mode <=3) iMode = mode;
			else {
				iMode = 0;
				StdOut("Invalid mode!, valid 0-3. Using default(0) tagging only mode.",GetName());
			}
		}
		void SetEmTimeCutMax(float Tmax_) { fEmTimeMax = Tmax_; }
		void SetEmTimeCutMin(float Tmin_) { fEmTimeMin = Tmin_; }
		void SetUseEmTimeCut(int c_) { iUseEmTimeCut = c_; }			//use em time to select beam halo

		//accessors
		float GetEmTimeCutMax() 	const { return fEmTimeMax; 	}
		float GetEmTimeCutMin() 	const { return fEmTimeMin; 	}
		int 	GetUseEmTimeCut() 	const { return iUseEmTimeCut;	}			//use em time to select beam halo
		HaloStuff GetHaloStuff(int phoIndex) const {		// the ordering will be determined by the InitSuperPhotons mod.
																		// since all mods are driven based on it, no worries!
			if ((phoIndex < 0) && (phoIndex > (int) vHaloStuff.size())) 
			{
				std::cout << __FILE__ << "::" << __LINE__ << ":: Requested index if out of range. Returning empty struct!" << std::endl;
				HaloStuff hs;
				return hs;
			}
			else return vHaloStuff.at(phoIndex);
		}

		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
	private:
		bool 	qMc;  				// true if MonteCarlo data
		GlobalCounters_t counter;
		HaloStuff myHaloStuff;
		std::vector<HaloStuff> vHaloStuff;			// keep all halo info for other mods to be accessed
		BeamHaloStudyHisto_t fNoBeamHaloHist;   	// histograms for without beam halo
		BeamHaloStudyHisto_t fWithBeamHaloHist; 	// histograms for photons with beam halo
		BeamHaloStudyHisto_t fBeamHaloByCutHist; 	// histograms for photons with beam halo from em time
		BeamHaloStudyHisto_t fBeamHaloByTimeHist; 	// histograms for photons with beam halo from em time
		//BeamHaloStudyHisto_t fBHScene0Hist,fBHScene1Hist,fBHScene2Hist;
		//BeamHaloStudyHisto_t fBHScene3Hist,fBHScene4Hist,fBHScene5Hist;

		int iMode;			// see EndJob() in cc file for details 
		InitSuperPhotons* initSpMod;
		float fEmTimeMax, fEmTimeMin;
		int iUseEmTimeCut;
		bool bRunPermit;			//make sure all dependencies are met

		std::vector<std::string> vCounterHistLabels;		//labels for counting histogram
		Histograms::CounterHists_t hCount;

		int iHalosFromS0, iHalosFromS1, iHalosFromS2;
		int iHalosFromS3, iHalosFromS4, iHalosFromS5;
		bool bNoSummary;			//control the end job summary print
		
	ClassDef(TagBeamHalo,1)
};


void FillHaloHistograms(Histograms::BeamHaloHists_t& Hist, TagBeamHalo::HaloStuff beamhalo,
								TStnPhoton* Pho); 	// common method to fill all beam halo hists

#endif
