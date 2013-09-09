#ifndef CP2CALIB_HH
#define CP2CALIB_HH
//_____________________________________________________________________________
// cp2 object  
/*
 *	$Id: cp2.hh,v 1.8 2010/07/02 17:36:04 samantha Exp $
 *		$Log: cp2.hh,v $
 *		Revision 1.8  2010/07/02 17:36:04  samantha
 *		MODIFIED: To look at events from different triggers to figure out the difference
 *		in the LERs I get from my ntuples and Calib.exe. Also added some event counters.
 *		
 *		Revision 1.7  2009/08/07 21:35:33  samantha
 *		can give name to the histogram that save all the final hists.
 *		
 *		Revision 1.6  2009/04/02 22:05:12  samantha
 *		This is modified version of 1.4 . renamed some of the hist that was supposed to be temp that started with sam_xxxx . I am still testing a good way to figure out the HV status.
 *		
 */
//_____________________________________________________________________________


#include "TH1.h"
#include "TH2.h"
#include <TFile.h>
#include <set>
#include <utility>

#include <Stntuple/loop/TStnModule.hh>
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include <Stntuple/obj/TCp2DataBlock.hh>
#include "Stntuple/obj/TL3SummaryBlock.hh"
#include <Stntuple/obj/TStnMetBlock.hh>
#include "Stntuple/obj/TStnVertexBlock.hh"
#include "Stntuple/obj/TStnTriggerBlock.hh"

#include <string>
#include <TTree.h>
#include "Stntuple/photon/TPhoTrigBits.hh"
#include "Stntuple/alg/TStntuple.hh"

#if !defined (__CINT__) || defined (__MAKECINT__)
#endif


class cp2: public TStnModule {

	protected:
		// pointers to the data blocks used,
		TStnHeaderBlock*    fHeaderBlock;
		TCp2DataBlock*      fCp2DataBlock;
		TL3SummaryBlock*    fL3SummaryBlock;
		TStnMetBlock*       fMetBlock;
		TStnVertexBlock*    fVertexBlock;
		TStnTriggerBlock*   fTriggerBlock;


	public:
		cp2(const char* name="cp2", const char* title="cp2");
		~cp2();

		struct Hist_t {
			TH1F* west; TH1F* east; TH1F* zerocp2; TH1F* cp2pad[48][54];
			TH1F* rancp2; TH1F* cp2tra;
			TH2F* cp2hit; TH2F* zero3d; TH2F* cp2Weight; TH2F* cp2ave;
			TH2F* ymoncp2;


		};


		struct LERHists_t {
			TH1F* Ratio;					// Ratio of LERS New/Old
			TH1F* OldLER;					// OLD LER
			TH1F* NewLER;					// NEW LER
			TH2F* PadHitsWgted_2D;		// 2D-Pad Hits weighted by ratio
			TH1F* WedgeHitsWgted;		//Wedge Hits weighted by LER ratio			
			TH1F* PadHitsWgted_1D; 		// 1D-Pad hist weighted by the LER 
			TH1F* PadHitsWgted_special; // 1D-oad hits weighted by LER for a specail wedge
			TH1F* sumet_large;
			TH1F* sumet_small;
			TH1F* met_large;
			TH1F* met_small;
			TH1F* ntracks_large;
			TH1F* ntracks_small;
			TH1F* avghitsperpad;
			TH2F* sumetVsNtrks;
			TH2F* metVsNtrks;
			TH2F* avghitsperpadVssumet;
			TH2F* avghitsperpadVsNtrks;
			TH1F* avghitsperpadNtrksRatio;
			TH1F* avghitsperpadSumetRatio;
			TH1F* SumetNtrksRatio;
			TH2F* avghitsperpadVsNVtx;
			TH2F* avghitsperpadVsNVtx12;
		};

		
		// ****** accessors

		Hist_t*            GetHist        () { return &fHist;        }
		TStnHeaderBlock*   GetHeaderBlock()  { return fHeaderBlock; }
		TCp2DataBlock*     GetCp2DataBlock() { return fCp2DataBlock; }
		TL3SummaryBlock*   GetL3SummaryBlock()  { return fL3SummaryBlock; }
		TStnVertexBlock*   GetVertexBlock()  const { return fVertexBlock;  }
		TStnTriggerBlock*  GetTriggerBlock() const { return fTriggerBlock; }

		// ****** setters

		// overloaded methods of TStnModule
		int     BeginJob();
		int     BeginRun();
		int     Event   (int ientry);
		int     EndJob  ();
		// ****** other methods
		void    BookHistograms();
		void    TrackExtr(float   radius, float  curvature, 
				float  Phi0, float  cotan, float  D_O, float  Z_0, 
				float* xyz, int &ieta, int &iphi, bool debug);

		void SetNewDataFileName(const std::string name)
		{ 
			if (name.find(".txt",0) == std::string::npos)
			{
				sNewDataFile = name + ".txt";
			} else 
			{
				sNewDataFile = name;
			}

		}

		void SetLastDBDataFileName(const std::string name) { sLastDBDataFile = name; }
		std::string GetNewDataFileName() const { return sNewDataFile;}
		std::string GetLastDBDataFileName() const { return sLastDBDataFile;}
		void MakeRatioLERhists();
		void HVStatus();
		void SetFirstLER(const float ler)
		{
			if (ler >0) fLERfirst = ler;
			else std::cout << __LINE__<<"::" << __LINE__ <<  ":: ler must be >0. LERfirst is set to default." << std::endl;
		}
		void SetLastDBLER(const float ler)
		{
			if (ler >0) fLERlastDB = ler;
			else std::cout << __LINE__<<"::" << __LINE__ <<  ":: ler must be >0. LERlastDB is set to default." << std::endl;
		}
		float GetFirstLER() const { return fLERfirst; }
		float GetLastDBLER() const { return fLERlastDB; }

		void SetMainLogFileName(const std::string name) { sMainLogFile = name; }
		std::string GetMainLogFileName() const { return sMainLogFile; }
		void SetFinalHistOutputFileName(const std::string name) { sFinalHistFile = name; }
		std::string GetFinalHistOutputFileName() const { return sFinalHistFile; }
		void SetTrigger(const std::string s) { sTrigName = s; }
		std::string GetTrigger() const { return sTrigName; }



	private:
		std::string sNewDataFile, sLastDBDataFile, sMainLogFile;//MainLogFile will have all the setting and the results one after another.
		std::string sFinalHistFile;				//Root file where all the final calibration hists will be written to
		float fLERfirst, fLERlastDB, fLERnew;
		// histograms filled
		Hist_t fHist;
		LERHists_t histLER;
		//var  for HV check
		std::set<std::pair<unsigned,unsigned> > foundHVoff;
		std::set<std::pair<unsigned,unsigned> > processedRunSec;
		std::set<std::pair<unsigned,unsigned> >::iterator myIt;

		//keep trees for
		TTree *processed;  //list of run#,sec# processed in this job
		TTree *hvOffList;		//list of run,sec when HV was off

		ofstream *mainLogFile;
		TPhoTrigBits 	trigbits;
		int iTrigBit;
		std::string sTrigName;
		unsigned int iEvtsProc[2];

	ClassDef(cp2,0)

};
#endif
