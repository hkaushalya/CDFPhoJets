///////////////////////////////////////////////////////////
// This is same as HaloByCutsTemp_Pho.cc/hh              //
// use this to skim beam halos from scenario 4,5         //
// 01-23-2008                                            //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>


#include "samantha/Pho/HaloByCutsTemp_Pho_Skim.hh"
#include "Stntuple/loop/TStnAna.hh"
#include <iostream>

ClassImp(HaloByCutsTemp_Pho_Skim)

//_____________________________________________________________________________
HaloByCutsTemp_Pho_Skim::HaloByCutsTemp_Pho_Skim(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bRunPermit(true)
{
	std::cout << "Hello I am HaloByCutsTemp_Pho_Skim module" << std::endl;
}

//_____________________________________________________________________________
HaloByCutsTemp_Pho_Skim::~HaloByCutsTemp_Pho_Skim() {
}

//_____________________________________________________________________________
void HaloByCutsTemp_Pho_Skim::SaveHistograms() {
}

//_____________________________________________________________________________
void HaloByCutsTemp_Pho_Skim::BookHistograms()
{
	DeleteHistograms();

	Histograms HistoMan;
	TFolder* new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		bRunPermit = false;
		return;
	} else {
		vCounterHistLabels.push_back("Events Processed");
		vCounterHistLabels.push_back("Events Passed");
		vCounterHistLabels.push_back("S-4 Events found");
		vCounterHistLabels.push_back("S-5 Events found");
		HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);
	}
}


//_____________________________________________________________________________
int HaloByCutsTemp_Pho_Skim::BeginJob()
{
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}

  	jetMod = (JetFilterModuleV2*) ((TStnAna*) GetAna()->GetModule("JetFilterV2"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		bRunPermit = false;
	}

  	tightMod = (TagTightPhotons*) ((TStnAna*) GetAna()->GetModule("TagTightPhotons"));
	if (!tightMod) {
		StdOut(__FILE__,__LINE__,2,"TightPhoton tag module not found.");
		bRunPermit = false;
	}
  	haloMod = (TagBeamHalo*) ((TStnAna*) GetAna()->GetModule("TagBeamHalo"));
	if (!haloMod) {
		StdOut(__FILE__,__LINE__,2,"Halo tag module not found.");
		bRunPermit = false;
	}

  	setTimeMod = (SetEmTimes*) ((TStnAna*) GetAna()->GetModule("SetEmTimes"));
	if (!setTimeMod) {
		StdOut(__FILE__,__LINE__,2,"SetEmTime module not found.");
		bRunPermit = false;
	}


				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.haloEvtsFound 		= 0;

	return 0;
}

//_____________________________________________________________________________
int HaloByCutsTemp_Pho_Skim::BeginRun()
{

	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int HaloByCutsTemp_Pho_Skim::Event(int ientry)
{
	// Event  loop

	SetPassed(0);
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,2,"Run Permit Failed. One or more dependencies or requirements not met.");
		exit (1);
		return 0;
	}


	counter.evtsRunOver++;
	hCount.iCounter->Fill(0);
  
  	int NsuperPho = initSpMod->GetSuperPhoSize();
	
	for (int i=0; i < NsuperPho; i++) {
		if (! initSpMod->GetSuperPhoton(i)->IsTightPhoton()) continue;
		if (! initSpMod->GetSuperPhoton(i)->IsBeamHalo()) continue;
		if (fabs(initSpMod->GetSuperPhoton(i)->GetEmTime()) > 4.8) continue;		//cut on em time
		
		int iHaloId = initSpMod->GetSuperPhoton(i)->GetBeamHaloId();
	
		//now look at each diff scenario with 1jet and 2jets
		if (jetMod->GetMyNjet(0,3) >= 1) {
			if ((iHaloId>>4) & 0x1) {
				std::cout << " ******************* found 4" << std::endl;
				GetHeaderBlock()->Print();
				SetPassed(1);
				hCount.iCounter->Fill(2);
				break;
			}
			if ((iHaloId>>5) & 0x1) {
				std::cout << " ******************* found 5" << std::endl;
				GetHeaderBlock()->Print();
				SetPassed(1);
				hCount.iCounter->Fill(3);
				break;
			}
		}
	
	}

	if (GetPassed()) {
		hCount.iCounter->Fill(1);
		counter.evtsPassModule++;
	}
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void HaloByCutsTemp_Pho_Skim::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void HaloByCutsTemp_Pho_Skim::FillDataBlocks(int ientry)
{
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int HaloByCutsTemp_Pho_Skim::EndJob() {
	std::string modname = GetName();
	std::string msg;
	msg  = "[HPS:00:]----- end job: ---- " + modname + "\n";
	msg += "[HPS:01:] Events Processed ----------- = " + ToStr(counter.evtsRunOver) + "\n";
	msg += "[HPS:02:] Events Pass this module ---- = " + ToStr(counter.evtsPassModule) + "\n";
	msg += "[HPS:03:] Halo Events Found ---------- = " + ToStr(counter.haloEvtsFound) + "\n";
	msg += "[HAT:04:] Halo+ 1 Jet Events --------- = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
	msg += "---------------------------------------------------\n";
	std::cout << msg;
	return 0;
}
