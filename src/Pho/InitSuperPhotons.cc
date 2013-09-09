/*


READ THIS EVERYDAY!
keep it simple and readable. forget about performance!
*/
/*
this module will setup photons to be used later. all it will do is sort photons, if any found,
set up 'super photons' .
NO good run/trigger/veretx or any other selection should be done here.


this will be the top level mod. all others will be built around this.
only trigger module may be run before this.

*/


#include <iostream>
#include <TMath.h>
#include <numeric>
#include <vector>
#include <algorithm>
#include "Stntuple/loop/TStnAna.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include "Stntuple/obj/TStnElectron.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/photon/TPhotonUtil.hh"

ClassImp(InitSuperPhotons)

//_____________________________________________________________________________
InitSuperPhotons::InitSuperPhotons(const char* name, const char* title):
   TStnModule(name,title),
	iPrintLevel(10),
  	bNoSummary(false)
{
	std::cout << "Hello I am InitSuperPhotons module." << std::endl;
}

//_____________________________________________________________________________
InitSuperPhotons::~InitSuperPhotons() {
}

//_____________________________________________________________________________
void InitSuperPhotons::SaveHistograms() {
}

//_____________________________________________________________________________
void InitSuperPhotons::BookHistograms()
{
	DeleteHistograms();
	BookGeneralHistograms(hGeneral,"General","General Histograms");
	BookPhotonHistograms(hPhoton,"Photons","Photon Histograms");
}


//_____________________________________________________________________________
int InitSuperPhotons::BeginJob()
{
				// register the data blocks

	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fPhotonBlock 	= (TStnPhotonBlock*)   RegisterDataBlock("PhotonBlock","TStnPhotonBlock");
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;


	return 0;
}

//_____________________________________________________________________________
int InitSuperPhotons::BeginRun()
{
	int currRun =  fHeaderBlock->RunNumber();
	//StdOut("Begining Run#",currRun,GetName());
	hGeneral.RunNumbers->Fill(currRun);
	return 0;
}

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int InitSuperPhotons::Event(int ientry)
{
	SetPassed(0);
	counter.evtsRunOver++;
  
	FillDataBlocks(ientry);
  
	if (fPhotonBlock->NPhotons() >0) {
		SetPassed(1);
		Cleanup();
		double zvx  = 0;
		if (fVertexBlock->GetBestVertex(12,1) != NULL) {
			zvx = fVertexBlock->GetBestVertex(12,1)->Z();
		}
		TStnEvent* event = GetEvent();
		TStntuple::CorrectPhotonEnergy(event,zvx);	//	sam: std minor corrections!
		TPhotonUtil::CorrectPhotonIsolation(event);	// sam,09/01/07: ray's isolation corrections
		TStntuple::CorrectElectronEnergy(event,zvx);	//	sam: 08-26-2008 I am making corrections to all electrons too. So I don't need 
																	//         worry about calling this again and again when getting matching electron
																	//         of a photon. see elog#692
		
		CreateSuperPhotonVector();
		//FillSuperPhotonHist();
	}

	if (GetPassed()) counter.evtsPassModule++;
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void InitSuperPhotons::Cleanup()
{
	superpho.clear();
} //Cleanup


/*-------------------------------------------------------------------*/
// counts the number of class 12 vertices in the event
/*-------------------------------------------------------------------*/
int InitSuperPhotons::Class12Vertices()
{
 	int nvtx = 0;
	int Nvtx = fVertexBlock->NVertices();

	for (int ivtx = 0; ivtx < Nvtx; ++ivtx) {
		TStnVertex* vert = fVertexBlock->Vertex(ivtx);
		if (vert->VClass() >= 12) ++nvtx;
	}
  
  return nvtx;

}  // Class12Vertices

/*-------------------------------------------------------------------*/
//	Sorts the Photon Block according to EtCorr, from highest to lowest
// & fill the photon info to my SuperPhoton list
/*-------------------------------------------------------------------*/
void InitSuperPhotons::CreateSuperPhotonVector()
{
	int Npho = fPhotonBlock->NPhotons();
	std::vector<TStnPhoton*> photon;
	std::vector<int> index;

	for (int i = 0; i < Npho ; i++) {
		photon.push_back(fPhotonBlock->Photon(i));
		index.push_back(i);
	}

	int loopupto = Npho-1;
	while (loopupto >0) {
		for (int i=0; i<loopupto;i++) {
			if ( (photon[i]->ECorr()*photon[i]->SinTheta()) <
					(photon[i+1]->ECorr()*photon[i+1]->SinTheta()) ) {

				TStnPhoton *temp = photon[i];
				photon[i] = photon[i+1];
				photon[i+1] = temp;

				int t = index[i];
				index[i] = index[i+1];
				index[i+1] = t;
			}
		} // for
		loopupto--;
	} // while

	for (unsigned int i=0; i < photon.size(); i++) {
		SuperPhoton sp;
		//sp.SetIndex(index[i]);   //No longer used --05-13-2008,sam
		sp.SetPhotonBlockIndex(index[i]);
		sp.SetElectronBlockIndex(GetElectronBlockIndex(photon[i]));
		sp.SetPhoton(photon[i]);
		sp.SetDetector(photon[i]->Detector());
		sp.SetDetEta(photon[i]->DetEta());
		sp.SetDetPhi(photon[i]->Phi());
		sp.SetEtRaw(photon[i]->Et());
		sp.SetECorr(photon[i]->ECorr());
		sp.SetEtCorr(photon[i]->ECorr() * photon[i]->SinTheta());
		sp.SetXCes(photon[i]->XCes());
		sp.SetZCes(photon[i]->ZCes());
		sp.SetHadEm(photon[i]->HadEm());
		sp.SetChi2Mean(photon[i]->Chi2Mean());
		sp.SetN3d(photon[i]->N3d());
		sp.SetIso4(photon[i]->EIso4(2));
		sp.SetTrkPt(photon[i]->Pt());
		sp.SetTrkPt2(photon[i]->Pt2());
		sp.SetTrkIso(photon[i]->SumPt4());
		sp.SetCesWireE2(photon[i]->CesWireE2());
		sp.SetCesStripE2(photon[i]->CesStripE2());
		sp.SetHadTdcTime(photon[i]->Time());
		sp.SetSinTheta(photon[i]->SinTheta());
		sp.SetCorVec(*(photon[i]->Momentum()));
		TLorentzVector vec;
		vec.SetPtEtaPhiM(photon[i]->Et(),photon[i]->Eta(),photon[i]->Phi(),0.0);
		sp.SetRawVec(vec);
		sp.SetNMuonStubs(photon[i]->NCosStubPho());	//N cosmic stubs in 30degree (in phi direction)
		sp.SetPrintLevel(GetPrintLevel());
		superpho.push_back(sp);
	}
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
int InitSuperPhotons::GetElectronBlockIndex(TStnPhoton* pho)
{
	if (!pho) {
		StdOut(__FILE__,__LINE__,3," pho points to NULL object. returning -1 for EleBlockIndex.");
		return -1;
	}

	int match_ele_ind = -1;

	float dmin=999.0;
	for (int i=0; i < fElectronBlock->NElectrons(); i++) {
		TStnElectron* ele = fElectronBlock->Electron(i);
		float dr = TStntuple::DeltaR(pho->TowPhi(),pho->TowEta(),
											ele->EmClusPhi(),ele->DetEta());
		if (dr < dmin) {
			match_ele_ind = i;
			dmin = dr;
		}
	}

	return match_ele_ind;
}

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void InitSuperPhotons::FillDataBlocks(int ientry)
{
	fPhotonBlock->GetEntry(ientry);
	fVertexBlock->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void InitSuperPhotons::FillSuperPhotonHist()
{
	for (int i = 0; i < GetSuperPhoSize() ; i++) {
			hPhoton.Detector->Fill(GetSuperPhoton(i)->GetDetector());
			hPhoton.DetEta->Fill(GetSuperPhoton(i)->GetDetEta());
			hPhoton.DetPhi->Fill(GetSuperPhoton(i)->GetDetPhi());
			hPhoton.EtCorr->Fill(GetSuperPhoton(i)->GetEtCorr());
			hPhoton.XCes->Fill(GetSuperPhoton(i)->GetXCes());
			hPhoton.ZCes->Fill(GetSuperPhoton(i)->GetZCes());
			hPhoton.HadEm->Fill(GetSuperPhoton(i)->GetHadEm());
			hPhoton.Chi2Mean->Fill(GetSuperPhoton(i)->GetChi2Mean());
			hPhoton.N3d->Fill(GetSuperPhoton(i)->GetN3d());
			hPhoton.Iso4->Fill(GetSuperPhoton(i)->GetIso4());
			hPhoton.TrkPt->Fill(GetSuperPhoton(i)->GetTrkPt());
			hPhoton.TrkIso->Fill(GetSuperPhoton(i)->GetTrkIso());
			hPhoton.Ces2Wire->Fill(GetSuperPhoton(i)->GetCesWireE2());
			hPhoton.Ces2Strip->Fill(GetSuperPhoton(i)->GetCesStripE2());
			hPhoton.HadTDCtime->Fill(GetSuperPhoton(i)->GetHadTdcTime());
	}
}

//_____________________________________________________________________________
void InitSuperPhotons::BookPhotonHistograms(PhotonPlots_t& hist,
							std::string sFoldName, std::string sFoldTitle)
{
	hist.Detector 	= new TH1F("Detector","Photon Detector Region",3,0,3);
	hist.DetEta 	= new TH1F("DetEta","Photon Detector Eta",100,-5,5);
	hist.DetPhi 	= new TH1F("DetPhi","Photon Detector Phi",105,-3.5,7);
	hist.EtCorr 	= new TH1F("EtCorr","Photon Corrected Et",600,0,120);
	hist.XCes 		= new TH1F("XCes","Photon XCes",640,-32,32);
	hist.ZCes 		= new TH1F("ZCes","Photon ZCes",2000,-250,250);
	hist.HadEm 		= new TH1F("HadEm","Photon HadEm",50,0,0.5);
	hist.Chi2Mean 	= new TH1F("Chi2Mean","Photon Chi2Mean (Wire+Strip/2)",800,0,80);
	hist.N3d 		= new TH1F("N3d","Photon N3d",10,0,10);
	hist.Iso4 		= new TH1F("IsoEtCorr","Photon IsoCorr",1500,-5,10);
	hist.TrkPt 		= new TH1F("Trkpt","Photon TrkPt",1000,0,100);
	hist.TrkIso 	= new TH1F("TrkIso","Photon TrkIso",150,0,15);
	hist.Ces2Wire 	= new TH1F("Ces2Wire","Photon CES(2nd) Wire",400,0,40);
	hist.Ces2Strip = new TH1F("Ces2Strip","Photon CES(2nd) Strip",400,0,40);
	hist.HadTDCtime = new TH1F("HadTDCtime","Photon Had TDC time",200,-50,150);

	//use GetBinWidth() in here!!!! u moran
	char ytitle[30];
	sprintf(ytitle,"Events/%.2f",hist.DetEta->GetBinWidth(1));
	hist.DetEta->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.DetPhi->GetBinWidth(1));
	hist.DetPhi->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.EtCorr->GetBinWidth(1));
	hist.EtCorr->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.XCes->GetBinWidth(1));
	hist.XCes->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.ZCes->GetBinWidth(1));
	hist.ZCes->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.HadEm->GetBinWidth(1));
	hist.HadEm->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.Chi2Mean->GetBinWidth(1));
	hist.Chi2Mean->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.N3d->GetBinWidth(1));
	hist.N3d->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.Iso4->GetBinWidth(1));
	hist.Iso4->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.TrkPt->GetBinWidth(1));
	hist.TrkPt->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.TrkIso->GetBinWidth(1));
	hist.TrkIso->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.Ces2Wire->GetBinWidth(1));
	hist.Ces2Wire->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.Ces2Strip->GetBinWidth(1));
	hist.Ces2Strip->SetYTitle(ytitle);
	sprintf(ytitle,"Events/%.2f",hist.HadTDCtime->GetBinWidth(1));
	hist.HadTDCtime->SetYTitle(ytitle);

	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	
	new_folder->Add(hist.Detector);
	new_folder->Add(hist.DetEta);
	new_folder->Add(hist.DetPhi);
	new_folder->Add(hist.EtCorr);
	new_folder->Add(hist.XCes);
	new_folder->Add(hist.ZCes);
	new_folder->Add(hist.HadEm);
	new_folder->Add(hist.Chi2Mean);
	new_folder->Add(hist.N3d);
	new_folder->Add(hist.Iso4);
	new_folder->Add(hist.TrkPt);
	new_folder->Add(hist.TrkIso);
	new_folder->Add(hist.Ces2Wire);
	new_folder->Add(hist.Ces2Strip);
	new_folder->Add(hist.HadTDCtime);


} //BookPhotonHistograms
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void InitSuperPhotons::BookGeneralHistograms(GeneralPlots_t& hist,
				std::string sFoldName, std::string sFoldTitle)
{
	hist.RunNumbers = new TH1F("RunNumbers","Run Numbers",270000,130000,400000);
	hist.RunNumbers->SetMarkerStyle(kPlus);
	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(hist.RunNumbers);	
} 

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int InitSuperPhotons::EndJob()
{
	if (GetSummaryStat()) return 0;

	std::cout << "[ISP:00:]----- end job: ---- " << GetName() << std::endl; 
	std::cout << "[ISP:01:] Events Processed ----------- = "  << counter.evtsRunOver << std::endl; 
	std::cout << "[ISP:02:] Events Passed -------------- = "  << counter.evtsPassModule << std::endl; 
	std::cout << "[ISP:03:] PrintLevel(def 10) --------- = "  << GetPrintLevel() << std::endl; 
	std::cout << "---------------------------------------------------" << std::endl; 
	return 0;
}
