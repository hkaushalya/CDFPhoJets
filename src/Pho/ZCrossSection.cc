#include "samantha/Pho/ZCrossSection.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/obj/TStnElectron.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include <iostream>
#include "TLorentzVector.h"
#include <vector>
ClassImp(ZCrossSection)

//_____________________________________________________________________________
ZCrossSection::ZCrossSection(const char* name, const char* title):
  TStnModule(name,title), qMc(false)
{
	std::cout << "Hello I am ZCrossSection module" << std::endl;
}

//_____________________________________________________________________________
ZCrossSection::~ZCrossSection() {
}

//_____________________________________________________________________________
void ZCrossSection::SaveHistograms() {
}

//_____________________________________________________________________________
void ZCrossSection::BookHistograms()
{
	DeleteHistograms();
	BookZHist(hZPlot, "Z_Plots", "Z properties");
	BookEleHist(hEle1Plot, "LeadEle_Plots","Leading Electron properties of the Z");
	BookEleHist(hEle2Plot, "Lead2Ele_Plots","2nd Leading Electron properties of the Z");
}

//_____________________________________________________________________________
int ZCrossSection::BeginJob()
{ 
  	tagEleMod = (TagElectrons*) ((TStnAna*) GetAna()->GetModule("TagElectrons"));
	if (!tagEleMod) {
		StdOut(__FILE__,__LINE__,3,"Required TagElectron module is not found!");
		return 0;
	}
  	tagConvEleMod = (TagConvElectrons*) ((TStnAna*) GetAna()->GetModule("TagConvElectrons"));
	if (!tagConvEleMod) {
		StdOut(__FILE__,__LINE__,3,"Required TagConvElectron module is not found!");
		return 0;
	}
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"Required InitSuperPhotons module is not found!");
		return 0;
	}

				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.Zevents = 0;

	return 0;
}

//_____________________________________________________________________________
int ZCrossSection::BeginRun()
{
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int ZCrossSection::Event(int ientry)
{
	SetPassed(1);
	counter.evtsRunOver++;
	if (GetPassed()) counter.evtsPassModule++;
   	
	int NsuperPho = initSpMod->GetSuperPhoSize();
	std::vector<SuperPhoton*> spvSuperPho;
	for (int i=0; i< NsuperPho; i++) {
		if (! initSpMod->GetSuperPhoton(i)->IsTightElectron()) continue;
		if (initSpMod->GetSuperPhoton(i)->IsConversion()) continue;
		spvSuperPho.push_back(initSpMod->GetSuperPhoton(i));
	}
	if (spvSuperPho.size() == 2) {
		counter.Zevents++;
		//TLorentzVector ele1_vec = spvSuperPho[0]->GetCorVec();
		//TLorentzVector ele2_vec = spvSuperPho[1]->GetCorVec();
		//TLorentzVector sum = ele1_vec + ele2_vec;
		FillZPlots(hZPlot, spvSuperPho[0], spvSuperPho[1]);
		FillElePlots(hEle1Plot, spvSuperPho[0]);
		FillElePlots(hEle2Plot, spvSuperPho[1]);
	}

	return 0;
} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void ZCrossSection::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void ZCrossSection::FillElePlots(ElePlots_t& hist, SuperPhoton* sp)
{
	hist.EtCorr->Fill(sp->GetEtCorr());
	hist.HadEm->Fill(sp->GetHadEm());
	hist.IsoEtCorr->Fill(sp->GetIso4());
	hist.TrkPt->Fill(sp->GetTrkPt());
	TStnEvent* evt = GetEvent();
	TStnElectron* Ele = TPhotonUtil::MatchElectron(evt,sp->GetPhoton()); // get matching electron
	if (Ele->TrackNumber() > 0) {
		hist.EoverP->Fill(Ele->EOverP());
	}
}

/*-------------------------------------------------------------------*/
void ZCrossSection::FillZPlots(ZPlots_t& hist, SuperPhoton* sp1, SuperPhoton* sp2)
{
	hist.ZEleEtratio->Fill(sp2->GetEtCorr()/sp1->GetEtCorr());
	TLorentzVector ele1_vec = sp1->GetCorVec();
	TLorentzVector ele2_vec = sp2->GetCorVec();
	TLorentzVector sum = ele1_vec + ele2_vec;
	hist.ZMass->Fill(sum.M());
	hist.ZPt->Fill(sum.Perp());
	hist.ZMet->Fill(fMetBlock->Met(0));
	hist.ZSumet0->Fill(fMetBlock->Sumet(0));
}
/*-------------------------------------------------------------------*/
void ZCrossSection::BookEleHist(ElePlots_t& hist, std::string sFoldName,
								std::string sFoldTitle)
{
	hist.EtCorr   	= new TH1F("EtCorr","Electron Etc",100,0,100);
	hist.HadEm		= new TH1F("HadEm","Electron HadEm",100,0,2);
	hist.IsoEtCorr	= new TH1F("Iso4","Electron Iso4",140,-4,10);
	hist.TrkPt 		= new TH1F("Trkpt","Electron TrkPt",500,0,100);
	hist.EoverP 	= new TH1F("EoverP","Electron E/p",100,0,5);

	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(hist.EtCorr);	
	new_folder->Add(hist.HadEm);	
	new_folder->Add(hist.IsoEtCorr);	
	new_folder->Add(hist.TrkPt);	
	new_folder->Add(hist.EoverP);	
}
//_____________________________________________________________________________
void ZCrossSection::BookZHist(ZPlots_t& hist, std::string sFoldName,
							std::string sFoldTitle)
{
	hist.ZMass 			= new TH1F("Z_mass","Z Invriant Mass",180,0,180);
	hist.ZPt 			= new TH1F("Z_Pt","Z Pt",100,0,200);
	hist.ZEleEtratio	= new TH1F("Z_eleEtratio","Z electrons Et ratio (2nd leading to leading) ",100,0,2);
	hist.ZSumet0 		= new TH1F("Z_Sumet0","Z Events Sumet(0)",100,0,500);
	hist.ZMet 			= new TH1F("Z_Met","Met for Z events",100,0,50);
	hist.ZeleDelPhi   = new TH1F("Z_EleDelPhi","#Delta #Phi of the two electrons of the Z events",700,-7,7);

	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(hist.ZMass);	
	new_folder->Add(hist.ZPt);
	new_folder->Add(hist.ZEleEtratio);	
	new_folder->Add(hist.ZSumet0);
	new_folder->Add(hist.ZMet);
	new_folder->Add(hist.ZeleDelPhi);
}
/*-------------------------------------------------------------------*/
void ZCrossSection::FillDataBlocks(int ientry)
{
	fMetBlock->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
}


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int ZCrossSection::EndJob() {

	printf("[ZCS:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[ZCS:01:]Events Run Over ------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "[ZCS:02:]Events Pass this module ---- = " << counter.evtsPassModule << std::endl;
	std::cout << "[ZCS:03:]Z Events ------------------- = " << counter.Zevents << std::endl;
	if (qMc)	StdOut(__FILE__,__LINE__,3,"This module is designed for data only. But the eun was on MC.");

	printf("ZCS::---------------------------------------------------\n");
	return 0;
}
