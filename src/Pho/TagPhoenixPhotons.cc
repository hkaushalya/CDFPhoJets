#include "samantha/Pho/TagPhoenixPhotons.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include <iostream>

ClassImp(TagPhoenixPhotons)

//_____________________________________________________________________________
TagPhoenixPhotons::TagPhoenixPhotons(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  iMode(0),
  bRunPermit(true),
  bNoSummary(false)
{
	std::cout << "Hello I am TagPhoenixPhotons module" << std::endl;
}

//_____________________________________________________________________________
TagPhoenixPhotons::~TagPhoenixPhotons() {
}

//_____________________________________________________________________________
void TagPhoenixPhotons::SaveHistograms() {
}

//_____________________________________________________________________________
void TagPhoenixPhotons::BookHistograms()
{
	DeleteHistograms();
	BookPhotonHistograms(fPhoenixHist,"Phoenix","Phoenix Photon Histos");
	BookPhotonHistograms(fNoPhoenixHist,"NonePhoenix","Non-Phoenix Photon histos");
	

	TFolder* new_folder = GetHistFolder(this,"Event","Event histos");
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
 	fPhxHist.Met = new TH1F("PhxEvtMet","Phoenix Event ME_{T}(0) ",1200,0,1200); 
	fPhxHist.Sumet = new TH1F("PhxSumEt","Phoenix Event SumE_{T}(0)",1200,0,1200); 
 	fNoPhxHist.Met = new TH1F("NoPhxEvtMet","Non-Phoenix Event ME_{T}(0) ",1200,0,1200); 
	fNoPhxHist.Sumet = new TH1F("NoPhxSumEt","Non-Phoenix Event SumE_{T}(0)",1200,0,1200); 

	new_folder->Add(fPhxHist.Met);
	new_folder->Add(fPhxHist.Sumet);
	new_folder->Add(fNoPhxHist.Met);
	new_folder->Add(fNoPhxHist.Sumet);

}


//_____________________________________________________________________________
int TagPhoenixPhotons::BeginJob()
{
	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initSpMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
		bRunPermit = false;
		return 0;
	}
				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	RegisterDataBlock("Phoenix_Electrons","TPhoenixElectronBlock",&fPhxEleBlock);
   RegisterDataBlock("PROD@PhoenixSI_Tracking","TStnTrackBlock",&fPhxSiTrackBlock);
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.phoenixCandidatesFound = 0;

	return 0;
}

//_____________________________________________________________________________
int TagPhoenixPhotons::BeginRun()
{
	if (fHeaderBlock->McFlag()) qMc = true;
	else 	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TagPhoenixPhotons::Event(int ientry)
{
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"One or more dependencies not found. pl check. Run Permit not cleared!");
		return 0;
	}
	SetPassed(1);
	counter.evtsRunOver++;
  
	FillDataBlocks(ientry);

  	int NsuperPho = initSpMod->GetSuperPhoSize();
	int iNphx = 0;
	for (int i=0; i< NsuperPho; i++) {
		bool bPhx = PhoenixCut(initSpMod->GetSuperPhoton(i), GetEvent());
		initSpMod->GetSuperPhoton(i)->SetPhoenixId(bPhx);
		if (bPhx) {
			iNphx++;
			//FillPhotonPlots(initSpMod->GetSuperPhoton(i),fPhoenixHist);
		} else {
			//FillPhotonPlots(initSpMod->GetSuperPhoton(i),fNoPhoenixHist);
		}
	}
	
	if (iNphx > 0) {
		counter.phoenixCandidatesFound++;
		//FillEventHist(fPhxHist);
	} else {
		//FillEventHist(fNoPhxHist);
	}
	
	if (GetPassed()) counter.evtsPassModule++;

	return 0;

} // Event

/*-------------------------------------------------------------------*/
//__ returns 1 if photon has a matched phoenix track
/*-------------------------------------------------------------------*/
int TagPhoenixPhotons::PhoenixCut(SuperPhoton* sp, TStnEvent* event)
{
	int passcode=0;
	TStnElectron* phele = TPhotonUtil::MatchPhoenixElectron(event,sp->GetPhoton());
	TStnTrack*    phTrk = TPhotonUtil::PhoenixTrack(event,phele);
	if(phTrk!=NULL) passcode=1; 
	return passcode;
}

/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TagPhoenixPhotons::Cleanup()
{

} //Cleanup

/*-------------------------------------------------------------------*/
void TagPhoenixPhotons::FillEventHist(EventHist_t& hist)
{ 
	hist.Met->Fill(fMetBlock->Met(0));
	hist.Sumet->Fill(fMetBlock->Sumet(0));
}
/*-------------------------------------------------------------------*/
void TagPhoenixPhotons::FillDataBlocks(int ientry)
{
	fPhxEleBlock->GetEntry(ientry);
	fPhxSiTrackBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);
}

/*-------------------------------------------------------------------*/
// fill photons plots
/*-------------------------------------------------------------------*/
void TagPhoenixPhotons::FillPhotonPlots(SuperPhoton* sp, PhotonPlots_t& PhotonPlot)
{
	PhotonPlot.Detector->Fill(sp->GetDetector());
	PhotonPlot.EtCorr->Fill(sp->GetEtCorr());
	PhotonPlot.XCes->Fill(sp->GetXCes());
	PhotonPlot.ZCes->Fill(sp->GetZCes());
	PhotonPlot.HadEm->Fill(sp->GetHadEm());
	PhotonPlot.IsoEtCorr->Fill(sp->GetIso4());
	PhotonPlot.Chi2Mean->Fill(sp->GetChi2Mean());
	PhotonPlot.N3d->Fill(sp->GetN3d());
	PhotonPlot.TrkPt->Fill(sp->GetTrkPt());
	PhotonPlot.TrkIso->Fill(sp->GetTrkIso());
	PhotonPlot.Ces2Wire->Fill(sp->GetCesWireE2());
	PhotonPlot.Ces2Strip->Fill(sp->GetCesStripE2());
	PhotonPlot.EvtMet->Fill(fMetBlock->Met(0));
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void TagPhoenixPhotons::BookPhotonHistograms(PhotonPlots_t& PhotonPlot,
							std::string sFoldName, std::string sFoldTitle)
{
	PhotonPlot.Detector 	= new TH1F("Detector","Photon Detector Region",3,0,3);
	PhotonPlot.EtCorr 	= new TH1F("EtCorr","Photon Corrected Et",600,0,120);
	PhotonPlot.XCes 		= new TH1F("XCes","Photon XCes",640,-32,32);
	PhotonPlot.ZCes 		= new TH1F("ZCes","Photon ZCes",2000,-250,250);
	PhotonPlot.HadEm 		= new TH1F("HadEm","Photon HadEm",50,0,0.5);
	PhotonPlot.IsoEtCorr = new TH1F("IsoEtCorr","Photon IsoCorr",1500,-5,10);
	PhotonPlot.Chi2Mean 	= new TH1F("Chi2Mean","Photon Chi2Mean (Wire+Strip/2)",800,0,80);
	PhotonPlot.N3d 		= new TH1F("N3d","Photon N3d",10,0,10);
	PhotonPlot.TrkPt 		= new TH1F("Trkpt","Photon TrkPt",1000,0,100);
	PhotonPlot.TrkIso 	= new TH1F("TrkIso","Photon TrkIso",150,0,15);
	PhotonPlot.Ces2Wire 	= new TH1F("Ces2Wire","Photon CES(2nd) Wire",400,0,40);
	PhotonPlot.Ces2Strip = new TH1F("Ces2Strip","Photon CES(2nd) Strip",400,0,40);
	PhotonPlot.HadTDCtime = new TH1F("HadTDCtime","Photon Had TDC time",200,-10,10);
	PhotonPlot.EvtMet 	= new TH1F("EvtMet","MEt(0) of the Event",150,0,150);

	PhotonPlot.EtCorr->SetYTitle("Events/0.2");
	PhotonPlot.XCes->SetYTitle("Events/0.1");
	PhotonPlot.ZCes->SetYTitle("Events/0.25");
	PhotonPlot.HadEm->SetYTitle("Events/0.01");
	PhotonPlot.IsoEtCorr->SetYTitle("Events/0.01");
	PhotonPlot.Chi2Mean->SetYTitle("Events/0.1");
	PhotonPlot.N3d->SetYTitle("Events/1.0");
	PhotonPlot.TrkPt->SetYTitle("Events/0.1");
	PhotonPlot.TrkIso->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Wire->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Strip->SetYTitle("Events/0.1");
	PhotonPlot.HadTDCtime->SetYTitle("Events/0.1");
	PhotonPlot.EvtMet->SetYTitle("Events/1GeV");

	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(PhotonPlot.Detector);
	new_folder->Add(PhotonPlot.EtCorr);
	new_folder->Add(PhotonPlot.XCes);
	new_folder->Add(PhotonPlot.ZCes);
	new_folder->Add(PhotonPlot.HadEm);
	new_folder->Add(PhotonPlot.IsoEtCorr);
	new_folder->Add(PhotonPlot.Chi2Mean);
	new_folder->Add(PhotonPlot.N3d);
	new_folder->Add(PhotonPlot.TrkPt);
	new_folder->Add(PhotonPlot.TrkIso);
	new_folder->Add(PhotonPlot.Ces2Wire);
	new_folder->Add(PhotonPlot.Ces2Strip);
	new_folder->Add(PhotonPlot.EvtMet);;


} //BookPhotonHistograms

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TagPhoenixPhotons::EndJob() {

	if (GetSummaryStat()) return 0;

	std::string sMsg;
	if (GetMode() == 0) sMsg = "Tag Only and Tag all photons.";
	if (GetMode() == 1) sMsg = "Filter Mode 1 (Remove PhxPhos and pass if a Photon is found)";
	if (GetMode() == 2) sMsg = "Filter Mode 2 (Tag all and pass if A Phoenix Photon is found)";
	if (GetMode() == 3) sMsg = "Filter Mode 3 (pass only if LEAD Photon is Phoenix)";

	printf("[TPP:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[TPP:01:] Events Processed------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "[TPP:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;
	std::cout << "[TPP:03:] Mode used ------------------ = " << GetMode() << " (" << sMsg << ")" << std::endl;
	std::cout << "[TPP:04:] Events with a Phoenix Photon = " << counter.phoenixCandidatesFound << std::endl;
	printf("---------------------------------------------------\n");
	return 0;
}
