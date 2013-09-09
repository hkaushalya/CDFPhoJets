#ifndef TRIGGERMODULE_HH
#define TRIGGERMODULE_HH


/* $Id: TriggerModule.hh,v 1.12 2011/05/23 20:30:31 samantha Exp $
 * $Log: TriggerModule.hh,v $
 * Revision 1.12  2011/05/23 20:30:31  samantha
 * MODIFIED: ,Minor spacing adjustments throughout.
 * ADDED: PassAnyPhoTriggers() to check if any trigger is passed.
 * RENAMED: 1. GetNClass12Vertices() -> GetN12vertex()
 * 2. GetNVertices() -> GetNvertex()
 *
 * Revision 1.11  2010/01/22 21:53:51  samantha
 * MOVED: All the heades out of the if( __CINT__ .. ) block.
 *
 * Revision 1.10  2009/11/10 19:12:58  samantha
 * Did a explicit conversion of unsigned to int in GetVtxZ() to remove compiler
 * warning. Added cvs log info into the file.
 *
 *
 */

#if !defined (__CINT__) || defined (__MAKECINT__)

#endif
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/loop/TStnModule.hh"
#include "TFolder.h"
#include "TString.h"
#include <string>
#include "Stntuple/obj/TStnTriggerBlock.hh"
#include "Stntuple/photon/TPhoTrigBits.hh"
#include "Stntuple/obj/TStnVertexBlock.hh"
#include "rlc/Pho/TGoodRun.hh"
#include "samantha/utils/FreeFunctions.hh"

class TriggerModule: public TStnModule {

	protected:
					// pointers to the data blocks used, header block is always
					// available via TStnModule::GetHeaderBlock()

		TStnHeaderBlock*    fHeaderBlock;
		TStnTriggerBlock*   fTriggerBlock;
		TStnVertexBlock*    fVertexBlock;


	public:
		TriggerModule(const char* name="TriggerModule", const char* title="TriggerModule");
		~TriggerModule();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  const { return fHeaderBlock;  }
		TStnTriggerBlock*  GetTriggerBlock() const { return fTriggerBlock; }
		TStnVertexBlock*   GetVertexBlock()  const { return fVertexBlock;  }


				/************** abstract objects *************************/


		struct Global_Counters_t {
			unsigned int evtsRunOver;		//___number of evts we ran over
			unsigned int evtsPassModule;	//___number of evts pass this module
		};

		struct EventParameters_t {
			bool bGoodRun; //true if good run is paased.
			int trig25;
			int trig50;
			int trig70;
			int Nverts;
			int N12verts;
			float BestVertZ;
			int BestVertNtrks;
			float BestVertSumPt;
		};

		// all the histograms
		struct GeneralHist_t {
			TH1F* Triggers;
			TH1F* Nvertices;
			TH1F* N12vertices;
			TProfile* NverticesVsRun;
			TProfile* N12verticesVsRun;
			TH1F* BestVertexZ;
		};

		struct RunDep_Hists_t {			//run dependent histograms counters
			TH1F* Nevts;			// events in the run
		};



		struct VTX {
			float Z;
			float SumPt;
			bool Class12;        // =1 if class 12, else 0
		};

		struct JetFilter_V2 {  					//info for JetFilterV2
			int VClass;								// fixed to class 12 for now
			std::vector<VTX> vAllVtx;			// all vertices sorted by SumPt highest to lowest
			int Nvtx;								// all vertices
			int NClass12vtx;								// number of class 12 vertices2vtx;								// number of class 12 vertices
			float BestClass12vtxZ;							// highest sum pt vtx z
		};




		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void SaveHistograms();
		void SetGoodRunListFile(TString filename) { goodRunListFile = filename; }
		TString GetGoodRunListFile() const { return goodRunListFile; }
	

		// setters
		void SetNeedTriggers(int tg = 1) { iNeedTrigger = tg; }
		int GetNeedTriggers() const { return iNeedTrigger; } // keep int so for diff triggers, diff ints can be used
																				// 0= for not required
																				// 1 = all three
																				// 2 = ISO_25, 3 = 50, 4 = 70
		void SetNeedGoodRun(const bool gd = true) { bNeedGoodRun = gd; }
		void SetNeedMinVtx(const int vt = 1) { iNMinVtx = vt; }
		void SetNeedMaxVtx(const int vt = 1000) { iNMaxVtx = vt; }
	
		void SetNeedVtxZcut(const bool b) { bNeedVtxZcut = b; }
		bool RequireGoodRun() const { return bNeedGoodRun; }
		bool RequireVtxZcut() const { return bNeedVtxZcut; }
		int GetMinVtx() const { return iNMinVtx; }
		int GetMaxVtx() const { return iNMaxVtx; }

		int GetTrig25IsoBit() const { return aEvtPara.trig25; }
		int GetTrig50Bit() const { return aEvtPara.trig50; }
		int GetTrig70Bit() const { return aEvtPara.trig70; }
		int GetNvertex() const { return aEvtPara.Nverts; }
		int GetN12vertex() const { return aEvtPara.N12verts; }
		float GetBestVertexZ() const { return aEvtPara.BestVertZ; }
		int GetBestVertexNTracks() const { return aEvtPara.BestVertNtrks; }
		float GetBestVertexSumPt() const { return aEvtPara.BestVertSumPt; }
		void SetJetFilterV2Stuff();

		int GetNClass12Vtx() const { return JetFilterV2Stuff.NClass12vtx; }
		int GetNVtx() const { return 	JetFilterV2Stuff.Nvtx; }
		float GetBestClass12VtxZ() const { return JetFilterV2Stuff.BestClass12vtxZ; }
		float GetVtxZ(int st, int Class ) const {		//st =1st, 2nd, 3rd highest sumpt vtx with Class12 or not
				if (st <=0) {
					StdOut(__FILE__,__LINE__,3,"Vtx request must be >1. vtx are in decreasing sum pt");
					exit(1);
				}
				//now search for the 1st,2nd .. highest sum pt vtx 
				int ST = 0;
				for (unsigned int i=0; i < JetFilterV2Stuff.vAllVtx.size(); ++i) {
					if (JetFilterV2Stuff.vAllVtx[i].Class12 != Class) continue;
					 ++ST;
					if (ST == st) return JetFilterV2Stuff.vAllVtx.at(i).Z;
				}
				return 0;	//if no vtx is found
		}

		float GetVtxZ(const int i) const {		//st =1st, 2nd, 3rd highest sumpt vtx no class check
				if (i>=0 &&  i < (int) JetFilterV2Stuff.vAllVtx.size()) return JetFilterV2Stuff.vAllVtx.at(i).Z;
				else { 
					StdOut(__FILE__,__LINE__,3,"Vtx request must be >=0. vtx are in decreasing sum pt");
					exit(1);
				}
		}

		
		float GetBestVtxZsep() const { 	//returns two best class 12 vtx z seperation
			float delz = 0;
			if (JetFilterV2Stuff.NClass12vtx > 1) {
					float z1 = GetVtxZ(1,1);
					float z2 = GetVtxZ(2,1);
					delz = fabs(z1-z2);
			}
			return delz;
			
		}
		
		void SetSummaryStat(const bool t) { bNoSummary = t; }
		bool GetSummaryStat() const { return bNoSummary; }
		bool PassGoodRun() const { return aEvtPara.bGoodRun; }
		bool PassTrigger() const;
		bool PassNVertexCut() const;
		bool PassVertexZcut() const;
		bool PassAnyPhoTriggers() const;
		bool PassIso25Trigger() const { return aEvtPara.trig25; }
		bool Pass50Trigger() const { return aEvtPara.trig50; }
		bool Pass70Trigger() const { return aEvtPara.trig70; }
	

	private:
		bool 	qMc;  				// true if MonteCarlo data
		Global_Counters_t counter;
		TPhoTrigBits 	trigbits;
		TGoodRun 		goodrun;
		TString 			goodRunListFile; // do not change to std::string, i don't own the GoodRun mod. Ray wrote it!
		GeneralHist_t  hGeneral_b4, hGeneral_a4;
		EventParameters_t aEvtPara;	
		bool bNeedGoodRun;
		int iNeedTrigger; 				// trigger can be broken to 3.
		int iNMinVtx,iNMaxVtx;			// min,max number of class 12 vertices required to pass.
		bool bNeedVtxZcut;  //switch for Best_Vtz_z < 60 used or not

		void Cleanup();
		void FillDataBlocks(int ientry) const;

		void SetVertexInfo();
		void SetTriggerInfo();
		void SetGoodRunInfo();
		void BookGeneralHistograms(GeneralHist_t& hist, std::string sFoldName,
												std::string sFoldTitle);
		void FillGeneralHist(GeneralHist_t& hist, EventParameters_t& para);
		void CleanUp();

		RunDep_Hists_t hRDHist;
		void BookRunDependentHistograms(RunDep_Hists_t& hist,
				std::string sFoldName, std::string sFoldTitle);

		bool bRunPermit;			//make sure all dependencies are met
		bool bNoSummary;			//control the end job summary print

	
		JetFilter_V2 JetFilterV2Stuff; 					//info for JetFilterV2
	
	ClassDef(TriggerModule,1)
};

#endif
