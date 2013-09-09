#ifndef TAGPMTSPIKES_HH
#define TAGPMTSPIKES_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "CalorGeometry/CalConstants.hh"
#include "Stntuple/data/TCalTower.hh"
#include "Stntuple/obj/TStnJetBlock.hh"
#include "Stntuple/obj/TStnPhotonBlock.hh"
#include "Stntuple/obj/TStnMetBlock.hh"
#include "samantha/utils/FreeFunctions.hh"
#include "TFolder.h"
#include <string>
#endif

class TagPMTSpikes: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TCalDataBlock*      fCalBlock;
		TStnJetBlock*       fJetBlock;
		TStnPhotonBlock*    fPhotonBlock;
		TStnMetBlock*       fMetBlock;
		typedef std::vector<TCalTower*> CalDataArray;  // need for beam halo
		typedef std::vector<TCalTower*>::iterator CalDataArrayI;   // need for BH

	public:
		TagPMTSpikes(const char* name="TagPMTSpikes", const char* title="TagPMTSpikes");
		~TagPMTSpikes();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock() { return fHeaderBlock; }
		TCalDataBlock*     GetCalDataBlock(){ return  fCalBlock; }
		TStnJetBlock*      GetJetBlock() 	{  return fJetBlock; }
		TStnPhotonBlock*   GetPhotonBlock() { return  fPhotonBlock; }
		TStnMetBlock*      GetMetBlock() 	{ return fMetBlock; }


				/************** abstract objects *************************/
	  struct CalHist_t {

		 TH2F*  fCalEMtwr;  // spike tower map for EM 
		 TH2F*  fCalCHAtwr;  // spike tower map for CHA 
		 TH2F*  fCalWHAtwr;  // spike tower map for WHA 
		 TH1F*  fCalHadTDCspike; // timing of HAD spikes;

		 TH1F*  fCalPmtEmAssm_bad;  // pmt assymetry for bad EM towers
		 TH1F*  fCalPmtHadAssm_bad;  // pmt assymetry for bad HAD towers
		 TProfile* fCalPmtAssmVsRun[19]; // pmt assymetry vs. Run num ber for bad towers
		 TH1F*  fCalBad_dPhi[19]; // dPhi=Phi(met)-Phi(bad tower) if |asymmetry|>0.55
		 TH1F*  fCalBad_MetRatio[19]; // MET/Et(bad tower) if |asymmetry|>0.55

		 TH1F*  fCalEm_dPhi; // dPhi=Phi(met)-Phi(em spike) if EmEnergy>10
		 TH1F*  fCalEm_MetRatio; // MET/Et(em spike) if EmEnergy>10 
		 TH1F*  fCalHad_dPhi; // dPhi=Phi(met)-Phi(had spike) if HadEnergy>10
		 TH1F*  fCalHad_MetRatio; // MET/Et(had spike) if HadEnergy>10 

		 TH1F*  fCalPmtEmAssm_le10;  // pmt assymetry for towers with EmEnergy<10
		 TH1F*  fCalPmtHadAssm_le10; // pmt assymetry for towers with HadEnergy<10
		 TH1F*  fCalPmtEmAssm_ge10;  // pmt assymetry for towers with EmEnergy>10
		 TH1F*  fCalPmtHadAssm_ge10; // pmt assymetry for towers with HadEnergy>10

		 TH1F*  fCalPmtCemAssm_le10;  //CEM pmt assymetry for towers with EmEnergy<10
		 TH1F*  fCalPmtChaAssm_le10; // CHA pmt assymetry for towers with HadEnergy<10, outside overlap
		 TH1F*  fCalPmtWhaAssm_le10; // WHA pmt assymetry for towers with HadEnergy<10, outside overlap
		 TH1F*  fCalPmtCemAssm_ge10;  // CEM pmt assymetry for towers with EmEnergy>10
		 TH1F*  fCalPmtChaAssm_ge10; // CHA pmt assymetry for towers with HadEnergy>10, outside overlap
		 TH1F*  fCalPmtWhaAssm_ge10; // WHA pmt assymetry for towers with HadEnergy>10, outside overlap

		 TH1F*  fCalPmtChaAssm_overlap_le10; // CHA pmt assymetry for towers with HadEnergy<10, inside overlap
		 TH1F*  fCalPmtWhaAssm_overlap_le10; // WHA pmt assymetry for towers with HadEnergy<10, inside overlap
		 TH1F*  fCalPmtChaAssm_overlap_ge10; // CHA pmt assymetry for towers with HadEnergy>10, inside overlap
		 TH1F*  fCalPmtWhaAssm_overlap_ge10; // WHA pmt assymetry for towers with HadEnergy>10, inside overlap

		 TH1F*  fCalPmtHadAssm_ge10_intime;  // pmt assymetry for towers with HadEnergy>10 and |HadTdc|<40.0 ns
		 TH1F*  fCalPmtHadAssm_ge10_outtime; // pmt assymetry for towers with HadEnergy>10 and |HadTdc|>40.0 ns
		 TH2F*  fCalEmPmtAssm_vs_PmtSum; // (pmt1-pmt2)/(pmt1+pmt2) for em energy
		 TH2F*  fCalHadPmtAssm_vs_PmtSum; // (pmt1-pmt2)/(pmt1+pmt2) for had energy
		 
		};

		struct BadTower {
			 int phiI;
			 int etaI;
			 int histoI;
			 int run1;
			 int run2;
		};



		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
			unsigned int spikeEvts;       //___number of evts with PMT spikes
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
	
		void SetUseFindSpike  (int cut)   { fUseFindSpike = cut; }
		void SetUseSpikeFilter(int cut)   { fUseSpikeFilter = cut; }
		void SetMinEmPMTasym(double cut)  { fMinEmPMTasym = cut; }
		void SetMinHadPMTasym(double cut) { fMinHadPMTasym = cut; }
		void SetMaxEmPMTasym(double cut)  { fMaxEmPMTasym = cut; }
		void SetMaxHadPMTasym(double cut) { fMaxHadPMTasym = cut; }
		void SetMinEmPmtE(double cut)     { fMinEmPmtE = cut; }
		void SetMinHadPmtE(double cut)    { fMinHadPmtE = cut; }

		void MatchCalorTowers(TStnJet* jet, TStnJetBlock* fJetBlock,
								  TCalDataBlock *fCalDataBlock,CalDataArray* cdholder);
		void MatchCalorTowers(int pho_ind, TStnPhotonBlock* fPhotonBlock,
									TCalDataBlock *fCalDataBlock,CalDataArray* cdholder);

		void CreateBadTowerList();
		int FindPmtSpikes(CalDataArray* towers, double& EmEpike, double& HadESpike, double& spikeHadTdc, CalHist_t& Hist);
		static int towerInd(int iEta) { return  iEta >= TENETA/2 ? iEta - TENETA/2 : -(TENETA/2 - 1 - iEta); }

		void 	BookCalHistograms(CalHist_t& Hist, std::string sFoldName, std::string sFoldTitle);

		int GetMode() { return iMode; }

		//setters
		void SetMode(int mode=0) {
			if (mode >=0 && mode <=2) iMode = mode;
			else {
				iMode = 0;
				StdOut("Invalid mode!, valid 0-2. Using default(0) tagging only mode.",GetName());
			}
		}

		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }

	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		int fUseFindSpike;		// status code for PMT spike studies
		int fUseSpikeFilter; 	// status code for PMT spike filter

		double fMinEmPMTasym;  	//minimum EM PMT asym
		double fMinHadPMTasym;	// minimum HAD PMT asym
		double fMaxEmPMTasym;	// maximum EM PMT asym
		double fMaxHadPMTasym;	// maximum HAD PMT asym

		double fMinEmPmtE;		//minimum EM PMT energy
		double fMinHadPmtE;		//minimum HAD PMT energy

		std::vector<BadTower> badEmTwr;  // list of bad EM towers
		std::vector<BadTower> badHadTwr; // list of bad HAD towers

		int Nspike_global;   // global number of spikes; to be used in cut
		CalHist_t CalHist;

		int iMode;			// see EndJob() in cc file for details 
		bool bNoSummary;			//control the end job summary print

	
	ClassDef(TagPMTSpikes,1)
};

#endif
