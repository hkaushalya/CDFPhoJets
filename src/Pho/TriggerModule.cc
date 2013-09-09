#include "samantha/Pho/TriggerModule.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnVertex.hh"
#include "samantha/utils/FreeFunctions.hh"
#include <iostream>
#include <assert.h>

ClassImp(TriggerModule)

//_____________________________________________________________________________
TriggerModule::TriggerModule(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bNeedGoodRun(0),
  iNeedTrigger(0),
  bNeedVtxZcut(0),
  bRunPermit(true),
  bNoSummary(false)
{
	std::cout << "Hello I am TriggerModule module" << std::endl;
}

//_____________________________________________________________________________
TriggerModule::~TriggerModule() {
}

//_____________________________________________________________________________
void TriggerModule::SaveHistograms() {
}

//_____________________________________________________________________________
void TriggerModule::BookHistograms()
{
	DeleteHistograms();
	BookGeneralHistograms(hGeneral_b4,"General_b4","General Histos before any cut");
	BookGeneralHistograms(hGeneral_a4,"General_a4","General Histos after any cut");
	BookRunDependentHistograms(hRDHist,"RunDepCounts","Run Dependent Counter plots");
}


//_____________________________________________________________________________
int TriggerModule::BeginJob()
{

				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fTriggerBlock 	= (TStnTriggerBlock*)  RegisterDataBlock("TriggerBlock","TStnTriggerBlock");
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");

	if (!fHeaderBlock || !fTriggerBlock || !fVertexBlock) {
		StdOut(__FILE__,__LINE__,3,"Required data Block not found. pl. check!.");
		bRunPermit = false;
	}
		

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;

	if (goodRunListFile.Length()<=0) {
		std::string sMsg("Good run list file not specified!");
		sMsg.append(GetGoodRunListFile());
		StdOut(__FILE__,__LINE__,2,sMsg);
		bRunPermit = false;
	} else {
		if (goodrun.Read(goodRunListFile.Data()) == 1) //failed reading goodrun file
		{
			assert (false && "TriggerModule::Error reading goodrun file from TGoodRun");
		}
	}
	trigbits.SetTriggerBlock(fTriggerBlock);

	//check for conflicting settings
	if (RequireVtxZcut() && GetMinVtx()<1)
	{
		StdOut(__FILE__, __LINE__, 3,"Conflicting setting used. Require cut on z position of vtx but min nvtx is <1?");
		assert(false);
	}
	return 0;
}

//_____________________________________________________________________________
int TriggerModule::BeginRun()
{

	int currRun =  fHeaderBlock->RunNumber();
	std::cout << " Begining Run# " << currRun << std::endl;

	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	trigbits.BeginRun();	
	TStntuple::Init(currRun);
	
	 int good = goodrun.Good(currRun);
	 std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__<< ": goodrun = " << good << " for run " << currRun << std::endl;

	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TriggerModule::Event(int ientry)
{
	SetPassed(0);
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"Run Permit not cleared!");
		exit (1);
	}

	//hRDHist.Nevts->Fill(fHeaderBlock->RunNumber()); 		// counts the number of events in each run
	counter.evtsRunOver++;
	
	CleanUp();
	FillDataBlocks(ientry);

	//this stuff is for JetFilerV2  so does not care if this module pass or fail
	// can keep it here. no hram.
	SetGoodRunInfo();
	SetVertexInfo();
	SetTriggerInfo();
	SetJetFilterV2Stuff();
	FillGeneralHist(hGeneral_b4,aEvtPara);
	
	bool bPass = true;

	if (RequireGoodRun()) {
		if (!PassGoodRun()) { bPass = false; std::cout << "good run failed !"; fHeaderBlock->Print(); }
	}

	if (!qMc) {				// Some MC data has not trig info. some only simulate L1 and L2 only. so drop this for MC data
		if (GetNeedTriggers()>0) {
			if (!PassTrigger()) { bPass = false;  /*std::cout << "trigger failed !"; fHeaderBlock->Print(); */}
		}
	}
		
	if (!PassNVertexCut()) { bPass = false;  /*std::cout << "nvtx failed !" << std::endl; fHeaderBlock->Print(); */}

	if (RequireVtxZcut() && ! PassVertexZcut()) { bPass = false;  /*std::cout << "Vtx Z failed !"; fHeaderBlock->Print();*/ }
	
	if (bPass) SetPassed(1);

	if (GetPassed()) {
		FillGeneralHist(hGeneral_a4,aEvtPara);
		counter.evtsPassModule++;
	}

	return 0;

} // Event

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TriggerModule::CleanUp()
{
	aEvtPara.bGoodRun 	= false;
	aEvtPara.trig25 		= 0;
	aEvtPara.trig50 		= 0;
	aEvtPara.trig70 		= 0;
	aEvtPara.Nverts 		= 0;
	aEvtPara.N12verts 	= 0;
	aEvtPara.BestVertZ 	= -9999.99;


	JetFilterV2Stuff.VClass = -999999;
	JetFilterV2Stuff.vAllVtx.clear();
	JetFilterV2Stuff.Nvtx   = -999999;
	JetFilterV2Stuff.NClass12vtx		= -999999;
		
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TriggerModule::FillGeneralHist(GeneralHist_t& hist, EventParameters_t& para)
{
	if (aEvtPara.trig25 || aEvtPara.trig50 || aEvtPara.trig70) hist.Triggers->Fill(1);
	else  hist.Triggers->Fill(0);
	if (aEvtPara.trig25) hist.Triggers->Fill(2);
	if (aEvtPara.trig50) hist.Triggers->Fill(3);
	if (aEvtPara.trig70) hist.Triggers->Fill(4);
	hist.Nvertices->Fill(aEvtPara.Nverts);
	hist.N12vertices->Fill(aEvtPara.N12verts);
	hist.BestVertexZ->Fill(aEvtPara.BestVertZ);
	int run = fHeaderBlock->RunNumber();
	hist.NverticesVsRun->Fill(run,aEvtPara.Nverts);
	hist.N12verticesVsRun->Fill(run,aEvtPara.N12verts);
}

//-------------------------------------------------
// this is for JetFilterModule v2. to get vertex info
// i used an external function before this to get this info
// but now i needed more than that. so need this. now i can refer to
// this and remove the external referral. but later. lets get v2 working!
//_____________________________________________________________________________
void TriggerModule::SetJetFilterV2Stuff()
{
	int VClass = 12;		// for now -sam
	JetFilterV2Stuff.VClass = VClass;
	
	TStnVertex   *v;

	JetFilterV2Stuff.Nvtx = fVertexBlock->NVertices();

	int nvtx12 = 0;	
	float bestvtx_pt = -1, bestvtx_z = -999999.9;
	for (int i=0; i< fVertexBlock->NVertices(); i++) {
   	v = fVertexBlock->Vertex(i);
		VTX vx;
		
		if (v->VClass() >= VClass) {
			if (v->SumPt() > bestvtx_pt) {
				bestvtx_pt = v->SumPt();
				bestvtx_z  = v->Z();
			}
			vx.Class12 = 1;
			nvtx12++;
		} else vx.Class12 = 0;

		vx.Z = v->Z();
      vx.SumPt = v->SumPt();
		JetFilterV2Stuff.vAllVtx.push_back(vx);
	}

	JetFilterV2Stuff.NClass12vtx = nvtx12;
	JetFilterV2Stuff.BestClass12vtxZ = bestvtx_z;
	
	int loopupto = (int)JetFilterV2Stuff.vAllVtx.size() - 1;
	while (loopupto >0) {
		for (int i=0; i<loopupto;i++) {
			if (JetFilterV2Stuff.vAllVtx[i].SumPt < JetFilterV2Stuff.vAllVtx[i+1].SumPt) {
				VTX vx = JetFilterV2Stuff.vAllVtx[i];
				JetFilterV2Stuff.vAllVtx[i] = JetFilterV2Stuff.vAllVtx[i+1];
				JetFilterV2Stuff.vAllVtx[i+1] = vx;
			}
		} // for
		loopupto--;
	} // while

	
}


/*-------------------------------------------------------------------*/
// cut on Nvtx position of the vertex,
/*-------------------------------------------------------------------*/
bool TriggerModule::PassNVertexCut() const
{
	if (aEvtPara.N12verts < GetMinVtx() || aEvtPara.N12verts > GetMaxVtx()) return false;
	else return true;
}

/*-------------------------------------------------------------------*/
// cut on Z position of the vertex,
/*-------------------------------------------------------------------*/
bool TriggerModule::PassVertexZcut() const
{
	if (fabs(GetBestVertexZ()) > 60.) return false;
	else return true;
}

/*-------------------------------------------------------------------*/
// set vertex, info
/*-------------------------------------------------------------------*/
void TriggerModule::SetVertexInfo()
{
	aEvtPara.N12verts = fVertexBlock->NVertices(12);
  if (fVertexBlock->GetBestVertex(12,1) != NULL) {
	  aEvtPara.BestVertNtrks = fVertexBlock->GetBestVertex(12,1)->NTracks();
	  aEvtPara.BestVertSumPt = fVertexBlock->GetBestVertex(12,1)->SumPt();
     aEvtPara.BestVertZ = fVertexBlock->GetBestVertex(12,1)->Z();
  }
}




/*-------------------------------------------------------------------*/
// good run preselection
// remove the bad runs/sections
/*-------------------------------------------------------------------*/
void TriggerModule::SetGoodRunInfo()
{
  int run = fHeaderBlock->RunNumber();
  int sec = fHeaderBlock->SectionNumber();
	
   aEvtPara.bGoodRun = goodrun.Good(run,sec);
	//std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__<< ": goodrun = " << aEvtPara.bGoodRun << " for run " << run << std::endl;
}


/*-------------------------------------------------------------------*/
// set trigger info
/*-------------------------------------------------------------------*/
void TriggerModule::SetTriggerInfo()
{
	trigbits.Event();
	aEvtPara.trig25 = trigbits.Pho25Iso();
	aEvtPara.trig50 = trigbits.Pho50();
	aEvtPara.trig70 = trigbits.Pho70();
}

bool TriggerModule::PassTrigger() const
{
	if (GetNeedTriggers() == 1 && PassAnyPhoTriggers()) return true;
	else if (GetNeedTriggers() == 2 && GetTrig25IsoBit()) return true;
	else if (GetNeedTriggers() == 3 && GetTrig50Bit()) return true;
	else if (GetNeedTriggers() == 4 && GetTrig70Bit()) return true;
	else return false;
}

/*-------------------------------------------------------------------*/
// get trigger info
/*-------------------------------------------------------------------*/
bool TriggerModule::PassAnyPhoTriggers() const
{
	if (aEvtPara.trig25 || aEvtPara.trig50 || aEvtPara.trig70) return true; 
	else return true;
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TriggerModule::FillDataBlocks(int ientry) const
{
	fTriggerBlock->GetEntry(ientry);
	fVertexBlock->GetEntry(ientry);
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TriggerModule::BookRunDependentHistograms(RunDep_Hists_t& hist,
				std::string sFoldName, std::string sFoldTitle)
{
	hist.Nevts = new TH1F("Nevts","Events per Run",270000,130000,400000);
	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(hist.Nevts);	
}
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TriggerModule::BookGeneralHistograms(GeneralHist_t& hist,
						std::string sFoldName, std::string sFoldTitle)
{
	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
  	char name [200];
  	char title[200];
	
	sprintf(name,"Triggers");
	sprintf(title,"Triggers used and Events passed each trigger");
	hist.Triggers = new TH1F(name,title,10,0,10);
	hist.Triggers->GetXaxis()->SetBinLabel(1,"pass none");
	hist.Triggers->GetXaxis()->SetBinLabel(2,"pass any one");
	hist.Triggers->GetXaxis()->SetBinLabel(3,"pass 25Iso");
	hist.Triggers->GetXaxis()->SetBinLabel(4,"pass 50");
	hist.Triggers->GetXaxis()->SetBinLabel(5,"pass 70");
	new_folder->Add(hist.Triggers);
	
	sprintf(name,"Nverts");
	sprintf(title,"Number of Vertices/Event");
	hist.Nvertices = new TH1F(name,title,15,0,15);
	new_folder->Add(hist.Nvertices);
	
	sprintf(name,"N12verts");
	sprintf(title,"Number of Class 12 vertices/Event");
	hist.N12vertices = new TH1F(name,title,15,0,15);
	new_folder->Add(hist.N12vertices);
	
	sprintf(name,"NvertsVsRun");
	sprintf(title,"Mean Number of Vertices/Event Vs Run");
	hist.NverticesVsRun = new TProfile(name,title,270000,130000,400000,0,15);
	new_folder->Add(hist.NverticesVsRun);
	
	sprintf(name,"N12vertsVsRun");
	sprintf(title,"Mean Number of Class 12 Vertices/Event Vs Run");
	hist.N12verticesVsRun = new TProfile(name,title,270000,130000,400000,0,15);
	new_folder->Add(hist.N12verticesVsRun);
	
	sprintf(name,"BestVertZ");
	sprintf(title,"Best Class 12 Vertex Z coordinate");
	hist.BestVertexZ = new TH1F(name,title,200,-200,200);
	new_folder->Add(hist.BestVertexZ);
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TriggerModule::EndJob() {

	if (GetSummaryStat()) return 0;
		
	std::string sTrigs, sGR;
	if (GetNeedTriggers()>=0 && !qMc)
	{
		if (GetNeedTriggers() == 0) sTrigs = "NONE";
		else if (GetNeedTriggers() == 1) sTrigs = "All 3 Photon triggers";
		else
		{
			if (GetNeedTriggers() == 2) sTrigs = "[25iso]";
			if (GetNeedTriggers() == 3) sTrigs += "[50]";
			if (GetNeedTriggers() == 4) sTrigs += "[70]";
		}
	} else {
		sTrigs = "NONE";
		if (qMc) sTrigs+=" (MC)";
	}
	
	if (RequireGoodRun()) sGR = "YES";
	else sGR = "NO";

	std::string sModName = GetName();
	std::string sMsg;
	sMsg  = "[TRM:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[TRM:01:] Events Processed ----------- = " + ToStr(counter.evtsRunOver) +"\n";
	sMsg += "[TRM:02:] Events Passed -------------- = " + ToStr(counter.evtsPassModule) +"\n";
	sMsg += "[TRM:03:] Data/Mc ? ------------------ = ";
	if (qMc)	sMsg += "MC\n";
	else	sMsg += "DATA\n";
	sMsg += "[TRM:04:] Good Run List File --------- = " + GetGoodRunListFile()+"\n";
	sMsg += "[TRM:05:] Triggers Used? ------------- = " + sTrigs +"\n";
	sMsg += "[TRM:06:] Good Run Used? ------------- = " + sGR +"\n";
	sMsg += "[TRM:07:] Min vertices required ------ = " + ToStr(GetMinVtx()) +"\n";
	sMsg += "[TRM:08:] Max vertices required ------ = " + ToStr(GetMaxVtx()) +"\n";
	sMsg += "[TRM:09:] Best Vtx z< 60 required----- = " + ToStr(RequireVtxZcut()) + " (0=NO, 1=YES)\n";

	sMsg += "---------------------------------------------------\n";
	std::cout << sMsg;
	return 0;
}
