#include "samantha/Pho/WSel.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/obj/TStnElectron.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include <iostream>
#include "TLorentzVector.h"
#include <vector>
ClassImp(WSel)

//_____________________________________________________________________________
WSel::WSel(const char* name, const char* title):
  TStnModule(name,title), qMc(false)
{
	std::cout << "Hello I am WSel module" << std::endl;
}

//_____________________________________________________________________________
WSel::~WSel() {
}

//_____________________________________________________________________________
void WSel::SaveHistograms() {
}

//_____________________________________________________________________________
void WSel::BookHistograms()
{
	DeleteHistograms();
	BookWHist(hWPlot, "W_Plots", "W properties");
	BookEleHist(hWElePlot, "LeadEle_Plots_a4","Leading Electron properties of the W","W electron");
	BookEleHist(hAllElePlot, "AllEle_Plots_b4","All Electron properties", "Before Cuts");
	BookRunDependentHistograms(hRDHist,"RunDepCounts","Run Dependent Counter plots");
	TFolder* new_folder = GetHistFolder(this, "MEt","MEt plots");
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	metb4cut = new TH1F("metb4","ME_{T} before cut",250,0,500);
	meta4cut = new TH1F("meta4","ME_{T} after cut",250,0,500);
	char ytitle[200];
	sprintf(ytitle,"Events/%.2f",metb4cut->GetBinWidth(1)); //true as long as both plots have same bin sizes
	metb4cut->GetXaxis()->SetTitle("ME_{T} (GeV)");
	metb4cut->GetYaxis()->SetTitle(ytitle);
	meta4cut->GetXaxis()->SetTitle("ME_{T} (GeV)");
	meta4cut->GetYaxis()->SetTitle(ytitle);
	
	new_folder->Add(metb4cut);	
	new_folder->Add(meta4cut);	
	BookAcceptanceHist(hAccept, "Accept_Plots","Acceptance plots");
}

//_____________________________________________________________________________
int WSel::BeginJob()
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
	counter.Wevents = 0;

	return 0;
}

//_____________________________________________________________________________
int WSel::BeginRun()
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
int WSel::Event(int ientry)
{
	SetPassed(1);
	hRDHist.Nevts->Fill(fHeaderBlock->RunNumber()); 		// counts the number of events in each run
	counter.evtsRunOver++;
	if (GetPassed()) counter.evtsPassModule++;
	
	FillDataBlocks(ientry);
	
	metb4cut->Fill(fMetBlock->Met(0));
	int NsuperPho = initSpMod->GetSuperPhoSize();
	std::vector<SuperPhoton*> spvSuperPho;
	for (int i=0; i< NsuperPho; i++) {
		FillAcceptancePlots(hAccept, initSpMod->GetSuperPhoton(i));
		FillElePlots(hAllElePlot, initSpMod->GetSuperPhoton(i));
		if (initSpMod->GetSuperPhoton(i)->GetTightElectronId() != 0) continue;
		if (initSpMod->GetSuperPhoton(i)->IsConversion()) continue;
		spvSuperPho.push_back(initSpMod->GetSuperPhoton(i));
	}
	if (spvSuperPho.size() == 1) {
			hAccept.Cumm->Fill(17);				// for 1 electron event
			hAccept.All->Fill(17);				// for 1 electron event
		if (fMetBlock->Met(0) > 20.) {
			hRDHist.NWs->Fill(fHeaderBlock->RunNumber()); 		// counts the number of W events in each run
			counter.Wevents++;
			meta4cut->Fill(fMetBlock->Met(0));
			hAccept.Cumm->Fill(18);				// pass met cut
			hAccept.All->Fill(18);				// pass met cut
			FillWPlots(hWPlot, spvSuperPho[0]);
			FillElePlots(hWElePlot, spvSuperPho[0]);
		}
	}

	return 0;
} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void WSel::Cleanup()
{

} //Cleanup
/*-------------------------------------------------------------------*/
void WSel::FillAcceptancePlots(AcceptPlots_t& hist, SuperPhoton* sp)
{
	int iNBits =  tagEleMod->GetElectronNIdBits();
	int iIdWord = sp->GetTightElectronId();

	// fill the cummulative passing plot
	for (int i=0 ; i < iNBits + 1; i++) {	// Electron id bits + conversion tag
		if (i == 0) {
			int iTemp = iIdWord & 0x1;
			if (iTemp == 0) hist.Cumm->Fill(0);
			else break;
		} else if (i == 1) {
			int iTemp = sp->IsConversion();
			if (iTemp == 0) hist.Cumm->Fill(1);	//conversion bucket
			else break;
		} else if (i >=2) {
			int iTemp = (iIdWord >> (i-1)) & 0x1;
			if (iTemp == 0) hist.Cumm->Fill(i);
			else break;
		}

	}
	//plot all passing (not cumulative)
	for (int i=0 ; i < iNBits + 1; i++) {	// Electron id bits + conversion tag
		if (i == 0) {
			int iTemp = iIdWord & 0x1;
			if (iTemp == 0) hist.All->Fill(0);
		} else if (i == 1) {
			int iTemp = sp->IsConversion();
			if (iTemp == 0) hist.All->Fill(1);	//conversion bucket
		} else if (i >=2) {
			int iTemp = (iIdWord >> (i-1)) & 0x1;
			if (iTemp == 0) hist.All->Fill(i);
		}

	}
	
}


/*-------------------------------------------------------------------*/
void WSel::FillElePlots(ElePlots_t& hist, SuperPhoton* sp)
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
	hist.Chi2Mean->Fill(sp->GetChi2Mean());
	hist.Ces2Wire->Fill(sp->GetCesWireE2());
	hist.Ces2Strip->Fill(sp->GetCesStripE2());
}

/*-------------------------------------------------------------------*/
void WSel::FillWPlots(WPlots_t& hist, SuperPhoton* sp)
{
	TLorentzVector ele_vec = sp->GetCorVec();
	TLorentzVector met_vec(0,0,0,0);
	float Met_0 = fMetBlock->Met(0);
	met_vec.SetPxPyPzE(Met_0 * TMath::Cos(fMetBlock->MetPhi(0)),Met_0 * TMath::Sin(fMetBlock->MetPhi(0)),0,Met_0);
	TLorentzVector sum = ele_vec + met_vec;
	hist.WMass->Fill(sum.M());
	hist.WPt->Fill(sum.Perp());
	hist.WMet->Fill(fMetBlock->Met(0));
	hist.WSumet0->Fill(fMetBlock->Sumet(0));
	hist.WHt->Fill(fMetBlock->Sumet(1));
	hist.WeleMetDelPhi->Fill(fabs(ele_vec.DeltaPhi(met_vec)));
	if (Met_0 !=0 ) hist.WEleEtMetRatio->Fill( (ele_vec.E() * TMath::Sin(ele_vec.Theta())) 
							/ (met_vec.E() * TMath::Sin(met_vec.Theta())) );
												//this is not needed as Sin(Theta) ==1 always in Pz=0
}
/*-------------------------------------------------------------------*/
void WSel::BookEleHist(ElePlots_t& hist, std::string sFoldName,
								std::string sFoldTitle, std::string sTextToAdd)
{
	char title[200];
	char ytitle[200];
	
	sprintf(title, "%s - E_{T}^{corr}",sTextToAdd.c_str());
	hist.EtCorr   	= new TH1F("EtCorr",title,200,0,200);
	sprintf(title, "%s - Had/Em",sTextToAdd.c_str());
	hist.HadEm		= new TH1F("HadEm",title,100,0,10);
	sprintf(title, "%s - Iso(0.4)",sTextToAdd.c_str());
	hist.IsoEtCorr	= new TH1F("Iso4",title,100,-10,90);
	sprintf(title, "%s - Track P_{T}",sTextToAdd.c_str());
	hist.TrkPt 		= new TH1F("Trkpt",title,200,0,200);
	sprintf(title, "%s - E/P",sTextToAdd.c_str());
	hist.EoverP 	= new TH1F("EoverP",title,500,0,100);
	sprintf(title, "%s - #chi^{2} Mean",sTextToAdd.c_str());
	hist.Chi2Mean = new TH1F("Chi2Mean",title,1210,-1,120);
	sprintf(title, "%s - CES(2nd) Wire",sTextToAdd.c_str());
	hist.Ces2Wire 	= new TH1F("Ces2Wire",title,400,-40,40);
	sprintf(title, "%s - CES(2nd) Strip",sTextToAdd.c_str());
	hist.Ces2Strip = new TH1F("Ces2Strip",title,400,-40,40);

	sprintf(ytitle, "Events/%.3f",hist.EtCorr->GetBinWidth(1));
	hist.EtCorr->GetYaxis()->SetTitle(ytitle);
	sprintf(ytitle, "Events/%.3f",hist.HadEm->GetBinWidth(1));
	hist.HadEm->GetYaxis()->SetTitle(ytitle);
	sprintf(ytitle, "Events/%.3f",hist.IsoEtCorr->GetBinWidth(1));
	hist.IsoEtCorr->GetYaxis()->SetTitle(ytitle);
	sprintf(ytitle, "Events/%.3f",hist.TrkPt->GetBinWidth(1));
	hist.TrkPt->GetYaxis()->SetTitle(ytitle);
	sprintf(ytitle, "Events/%.3f",hist.EoverP->GetBinWidth(1));
	hist.EoverP->GetYaxis()->SetTitle(ytitle);
	sprintf(ytitle, "Events/%.3f",hist.Chi2Mean->GetBinWidth(1));
	hist.Chi2Mean->GetYaxis()->SetTitle(ytitle);
	sprintf(ytitle, "Events/%.3f",hist.Ces2Wire->GetBinWidth(1));
	hist.Ces2Wire->GetYaxis()->SetTitle(ytitle);
	sprintf(ytitle, "Events/%.3f",hist.Ces2Strip->GetBinWidth(1));
	hist.Ces2Strip->GetYaxis()->SetTitle(ytitle);

	hist.EtCorr->GetXaxis()->SetTitle("E_{T} (GeV)");
	hist.HadEm->GetXaxis()->SetTitle("Had/Em");
	hist.IsoEtCorr->GetXaxis()->SetTitle("Isolation E_{T} (GeV)");
	hist.TrkPt->GetXaxis()->SetTitle("Track P_{T} (GeV)");
	hist.EoverP->GetXaxis()->SetTitle("E/P");
	hist.Chi2Mean->GetXaxis()->SetTitle("#Chi^{2} Mean");
	hist.Ces2Wire->GetXaxis()->SetTitle("2nd CES WIRE cluster energy");
	hist.Ces2Strip->GetXaxis()->SetTitle("2nd CES STRIP cluster energy");

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
	new_folder->Add(hist.Chi2Mean);	
	new_folder->Add(hist.Ces2Wire);	
	new_folder->Add(hist.Ces2Strip);	
}
//_____________________________________________________________________________
void WSel::BookWHist(WPlots_t& hist, std::string sFoldName,
							std::string sFoldTitle)
{
	hist.WMass 			= new TH1F("W_mass","W Transverse Invriant Mass",250,0,250);
	hist.WPt 			= new TH1F("W_Pt","W P_{T}",250,0,250);
	hist.WEleEtMetRatio	= new TH1F("W_eleEtMetRatio","Electron E_{T} to ME_{T} ratio ",100,0,10);
	hist.WSumet0 		= new TH1F("W_Sumet0","W Events #Sigma E_{T}",150,0,1500);
	hist.WHt 			= new TH1F("W_Ht","W Events #H_{T}",150,0,1500);
	hist.WMet 			= new TH1F("W_Met","Met for W events",250,0,250);
	hist.WeleMetDelPhi   = new TH1F("W_EleMetDelPhi","#Delta #Phi of electron and ME_{T}",700,-7,7);

	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	char ytitle[200];
	sprintf(ytitle,"Events/%.3f",hist.WMass->GetBinWidth(1));
	hist.WMass->GetYaxis()->SetTitle(ytitle);
	hist.WMass->GetXaxis()->SetTitle("Transverse Invariant Mass (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.WPt->GetBinWidth(1));
	hist.WPt->GetYaxis()->SetTitle(ytitle);
	hist.WPt->GetXaxis()->SetTitle("P_{T} (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.WEleEtMetRatio->GetBinWidth(1));
	hist.WEleEtMetRatio->GetYaxis()->SetTitle(ytitle);
	hist.WEleEtMetRatio->GetXaxis()->SetTitle("#frac{Ele. E_{T}}{ME_{T}}");
	sprintf(ytitle,"Events/%.3f",hist.WSumet0->GetBinWidth(1));
	hist.WSumet0->GetYaxis()->SetTitle(ytitle);
	hist.WSumet0->GetXaxis()->SetTitle("#Sigma E_{T} (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.WHt->GetBinWidth(1));
	hist.WHt->GetYaxis()->SetTitle(ytitle);
	hist.WHt->GetXaxis()->SetTitle("H_{T} (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.WMet->GetBinWidth(1));
	hist.WMet->GetYaxis()->SetTitle(ytitle);
	hist.WMet->GetXaxis()->SetTitle("ME_{T} (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.WeleMetDelPhi->GetBinWidth(1));
	hist.WeleMetDelPhi->GetYaxis()->SetTitle(ytitle);
	hist.WeleMetDelPhi->GetXaxis()->SetTitle("#delta #phi");

	
	new_folder->Add(hist.WMass);	
	new_folder->Add(hist.WPt);
	new_folder->Add(hist.WEleEtMetRatio);	
	new_folder->Add(hist.WSumet0);
	new_folder->Add(hist.WHt);
	new_folder->Add(hist.WMet);
	new_folder->Add(hist.WeleMetDelPhi);
}
/*-------------------------------------------------------------------*/
void WSel::BookAcceptanceHist(AcceptPlots_t& hist, std::string sFoldName,
								std::string sFoldTitle)
{
	hist.Cumm   = new TH1F("Cumm","Passing Cumulative. Each bin has passed all cuts to the left of it.",25,-.5,24.5);
	hist.Cumm->GetXaxis()->SetBinLabel(1,"Central");
	hist.Cumm->GetXaxis()->SetBinLabel(2,"Not Conversion");
	hist.Cumm->GetXaxis()->SetBinLabel(3,"E_{T}");
	hist.Cumm->GetXaxis()->SetBinLabel(4,"HadEm");
	hist.Cumm->GetXaxis()->SetBinLabel(5,"Chi2");
	hist.Cumm->GetXaxis()->SetBinLabel(6,"Tracks>=1");
	hist.Cumm->GetXaxis()->SetBinLabel(7,"Tracks>2");
	hist.Cumm->GetXaxis()->SetBinLabel(8,"2nd Track Pt");
	hist.Cumm->GetXaxis()->SetBinLabel(9,"Matching Ele");
	hist.Cumm->GetXaxis()->SetBinLabel(10,"Ele Z0");
	hist.Cumm->GetXaxis()->SetBinLabel(11,"Ele Trk Pt");
	hist.Cumm->GetXaxis()->SetBinLabel(12,"EoverP");
	hist.Cumm->GetXaxis()->SetBinLabel(13,"Iso4");
	hist.Cumm->GetXaxis()->SetBinLabel(14,"TrkIso");
	hist.Cumm->GetXaxis()->SetBinLabel(15,"2nd Ces");
	hist.Cumm->GetXaxis()->SetBinLabel(16,"XCes");
	hist.Cumm->GetXaxis()->SetBinLabel(17,"ZCes");
	hist.Cumm->GetXaxis()->SetBinLabel(18,"One e");
	hist.Cumm->GetXaxis()->SetBinLabel(19,"ME_{T}>20");

	hist.All		= new TH1F("All","All Passing",25,-.5,24.5);
	hist.All->GetXaxis()->SetBinLabel(1,"Central");
	hist.All->GetXaxis()->SetBinLabel(2,"Not Conversion");
	hist.All->GetXaxis()->SetBinLabel(3,"E_{T}");
	hist.All->GetXaxis()->SetBinLabel(4,"HadEm");
	hist.All->GetXaxis()->SetBinLabel(5,"Chi2");
	hist.All->GetXaxis()->SetBinLabel(6,"Tracks>=1");
	hist.All->GetXaxis()->SetBinLabel(7,"Tracks>2");
	hist.All->GetXaxis()->SetBinLabel(8,"2nd Track Pt");
	hist.All->GetXaxis()->SetBinLabel(9,"Matching Ele");
	hist.All->GetXaxis()->SetBinLabel(10,"Ele Z0");
	hist.All->GetXaxis()->SetBinLabel(11,"Ele Trk Pt");
	hist.All->GetXaxis()->SetBinLabel(12,"EoverP");
	hist.All->GetXaxis()->SetBinLabel(13,"Iso4");
	hist.All->GetXaxis()->SetBinLabel(14,"TrkIso");
	hist.All->GetXaxis()->SetBinLabel(15,"2nd Ces");
	hist.All->GetXaxis()->SetBinLabel(16,"XCes");
	hist.All->GetXaxis()->SetBinLabel(17,"ZCes");
	hist.All->GetXaxis()->SetBinLabel(18,"One e");
	hist.All->GetXaxis()->SetBinLabel(19,"ME_{T}>20");


	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(hist.Cumm);	
	new_folder->Add(hist.All);	
	
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void WSel::BookRunDependentHistograms(RunDep_Hists_t& hist,
				std::string sFoldName, std::string sFoldTitle)
{
	hist.Nevts = new TH1F("Nevts","Events per Run",270000,130000,400000);
	hist.NWs = new TH1F("NW","W Events per Run",270000,130000,400000);
	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(hist.Nevts);	
	new_folder->Add(hist.NWs);	
}



/*-------------------------------------------------------------------*/
void WSel::FillDataBlocks(int ientry)
{
	fMetBlock->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
}


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int WSel::EndJob()
{

	printf("[WSL:00:]----- end job: ---- %s\n",GetName());
	if (qMc)	StdOut(__FILE__,__LINE__,3,"This module is designed for data only. But the run is on MC.");
	std::cout << "[WSL:01:] Events Run Over ------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "[WSL:02:] Events Pass this module ---- = " << counter.evtsPassModule << std::endl;
	std::cout << "[WSL:03:] W Events ------------------- = " << counter.Wevents << std::endl;
	hAccept.Cumm->Print();
	for (int i=0;i < hAccept.Cumm->GetNbinsX() ; i++) std::cout << "Wbin"<< i+1 << "=" << hAccept.Cumm->GetBinContent(i+1) << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
