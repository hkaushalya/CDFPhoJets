#include "samantha/Pho/ZSel.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/obj/TStnElectron.hh"
#include "Stntuple/obj/TStnEvent.hh"
#include <iostream>
#include "TLorentzVector.h"
#include <vector>
#include <assert.h>

ClassImp(ZSel)

//_____________________________________________________________________________
ZSel::ZSel(const char* name, const char* title):
  TStnModule(name,title), qMc(false)
{
	std::cout << "Hello I am ZSel module" << std::endl;
}

//_____________________________________________________________________________
ZSel::~ZSel() {
}

//_____________________________________________________________________________
void ZSel::SaveHistograms() {
}

//_____________________________________________________________________________
void ZSel::BookHistograms()
{
	DeleteHistograms();
	BookZHist(hZPlot, "Z_Plots", "Z properties");
	BookEleHist(hAllElePlot, "AllEle_Plots_b4","All Electron properties", "Before Cuts");
	BookEleHist(hEle1Plot, "LeadEle_Plots","Leading Electron properties of the Z", "Lead e");
	BookEleHist(hEle2Plot, "Lead2Ele_Plots","2nd Leading Electron properties of the Z","2nd e");
	TFolder* new_folder = GetHistFolder(this, "MEt","MEt plots");
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}

	nvtx12b4cut = new TH1F("nvtx12b4","Nvtx12 before cut",10,0,10);
	nvtx12a4cut = new TH1F("nvtx12a4","Nvtx12 after cut",10,0,10);
	
	metb4cut = new TH1F("metb4","ME_{T} before cut",300,0,300);
	meta4cut = new TH1F("meta4","ME_{T} after cut",300,0,300);
	char ytitle[200];
	sprintf(ytitle,"Events/%.2f",metb4cut->GetBinWidth(1)); //true as long as both plots have same bin sizes
	metb4cut->GetXaxis()->SetTitle("ME_{T} (GeV)");
	metb4cut->GetYaxis()->SetTitle(ytitle);
	meta4cut->GetXaxis()->SetTitle("ME_{T} (GeV)");
	meta4cut->GetYaxis()->SetTitle(ytitle);
	
	new_folder->Add(metb4cut);	
	new_folder->Add(meta4cut);	
	new_folder->Add(nvtx12a4cut);	
	new_folder->Add(nvtx12b4cut);	
	BookAcceptanceHist(hAccept, "Accept_Plots","Acceptance plots");
	BookRunDependentHistograms(hRDHist,"RunDepCounts","Run Dependent Counter plots");
}

//_____________________________________________________________________________
int ZSel::BeginJob()
{ 
  	tagEleMod = (TagElectrons*) ((TStnAna*) GetAna()->GetModule("TagElectrons"));
	if (!tagEleMod) {
		StdOut(__FILE__,__LINE__,3,"Required TagElectron module is not found!");
		assert(false);
	}
  	tagConvEleMod = (TagConvElectrons*) ((TStnAna*) GetAna()->GetModule("TagConvElectrons"));
	if (!tagConvEleMod) {
		StdOut(__FILE__,__LINE__,3,"Required TagConvElectron module is not found!");
		assert(false);
	}
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"Required InitSuperPhotons module is not found!");
		assert(false);
	}

  	jetMod = (JetFilterModuleV3*) ((TStnAna*) GetAna()->GetModule("JetFilterV3"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilterModuleV3 module required!.");
		assert(false);
	}
  	trigMod = (TriggerModule*) ((TStnAna*) GetAna()->GetModule("TriggerModule"));
	if (!trigMod) {
		StdOut(__FILE__,__LINE__,3,"TriggerModule required!.");
		assert(false);
	}

				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");

	assert (fHeaderBlock != NULL && "ZSel::BeginJob:: fHeaderBlock is null");
	assert (fElectronBlock != NULL && "ZSel::BeginJob:: fElectronBlock is null");
	assert (fMetBlock != NULL && "ZSel::BeginJob:: fMetBlock is null");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.Zevents = 0;

	if (GetUseNvtxWeigts())
	{
		SetNvtxWeights();
	}
	return 0;
}

//_____________________________________________________________________________
int ZSel::BeginRun()
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
int ZSel::Event(int ientry)
{
	SetPassed(1);
	hRDHist.Nevts->Fill(fHeaderBlock->RunNumber()); 		// counts the number of events in each run
	counter.evtsRunOver++;
	if (GetPassed()) counter.evtsPassModule++;
   	
	FillDataBlocks(ientry);

	//metb4cut->Fill(fMetBlock->Met(0));
	metb4cut->Fill(jetMod->GetMyMetCorr(0,3));
	nvtx12b4cut->Fill(trigMod->GetN12vertex());
	int NsuperPho = initSpMod->GetSuperPhoSize();
	std::vector<SuperPhoton*> spvSuperPho;
	for (int i=0; i< NsuperPho; i++) {
		FillAcceptancePlots(hAccept, initSpMod->GetSuperPhoton(i));
		FillElePlots(hAllElePlot, initSpMod->GetSuperPhoton(i));
		if (initSpMod->GetSuperPhoton(i)->GetTightElectronId() != 0) continue;  //pho-like tight electrons
		if (initSpMod->GetSuperPhoton(i)->IsConversion()) continue;
		if (initSpMod->GetSuperPhoton(i)->GetEtCorr()< GetMinEleEt()) continue;
		if (fabs(initSpMod->GetSuperPhoton(i)->GetDetEta())> GetMaxDetEta()) continue;
					
		spvSuperPho.push_back(initSpMod->GetSuperPhoton(i));
	}
	// I have modified this to study Zee EWK MC and DATA for MET plot
	// in g+jets. Since the extra em-objetcs becomes jets in my analyses,
	// this is little tricky. I'll try to keep at least one leg in the central
	// so that the radiated photon will probably be from it and allow 2nd electron
	// to be in central or plug.
	
	if (spvSuperPho.size() >= 2) {
		
		// now check if both legs are 
		const float fEta1 = spvSuperPho.at(0)->GetDetEta();
		const float fEta2 = spvSuperPho.at(1)->GetDetEta();
		//assuming pre-selection is done with upto plug region
		if ( fabs(fEta1) > 1.1 && fabs (fEta2) >1.1 ) return 0; //both z legs are in plug, so ignore

		const TLorentzVector ele1_vec = spvSuperPho.at(0)->GetCorVec();
		const TLorentzVector ele2_vec = spvSuperPho.at(1)->GetCorVec();
		const TLorentzVector sum = ele1_vec + ele2_vec;
		if (sum.M() < 76. || sum.M() > 106.) return 0 ; //std. Z mass window cut
		
		float fWgt = 1;
		if (GetUseNvtxWeigts())
		{
			const int iNvtx12 = trigMod->GetN12vertex();
			if (iNvtx12 <= vNvtxWgts.size()) fWgt = vNvtxWgts.at(iNvtx12 - 1);
		}
		counter.Zevents++;
		hRDHist.NZs->Fill(fHeaderBlock->RunNumber()); 		// counts the number of Z events in each run
		//meta4cut->Fill(fMetBlock->Met(0));
		meta4cut->Fill(jetMod->GetMyMetCorr(0,3), fWgt);
		nvtx12a4cut->Fill(trigMod->GetN12vertex(), fWgt);
		hAccept.Cumm->Fill(17);				// for 2 electron event
		hAccept.All->Fill(17);				// for 2 electron event
		FillZPlots(hZPlot, spvSuperPho[0], spvSuperPho[1], fWgt);
		FillElePlots(hEle1Plot, spvSuperPho[0]);
		FillElePlots(hEle2Plot, spvSuperPho[1]);
	}

	return 0;
} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void ZSel::Cleanup()
{

} //Cleanup

/*-------------------------------------------------------------------*/
void ZSel::FillAcceptancePlots(AcceptPlots_t& hist, SuperPhoton* sp)
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
void ZSel::FillElePlots(ElePlots_t& hist, SuperPhoton* sp)
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


	const TVector2 tv2MetVec(jetMod->GetMyMetXCorr(0,3), jetMod->GetMyMetYCorr(0,3));
	const TLorentzVector tlEleVec = sp->GetCorVec();
	const float dphi_em =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlEleVec.Phi())));

	hist.DelPhiEleMet->Fill(dphi_em);
	
}

/*-------------------------------------------------------------------*/
void ZSel::FillZPlots(ZPlots_t& hist, SuperPhoton* sp1, SuperPhoton* sp2, const float fWgt)
{
	hist.ZEleEtratio->Fill(sp2->GetEtCorr()/sp1->GetEtCorr(), fWgt);
	TLorentzVector ele1_vec = sp1->GetCorVec();
	TLorentzVector ele2_vec = sp2->GetCorVec();
	hist.ZeleDelPhi->Fill(fabs(ele1_vec.DeltaPhi(ele2_vec)), fWgt);
	TLorentzVector sum = ele1_vec + ele2_vec;
	hist.ZMass->Fill(sum.M(), fWgt);
	hist.ZPt->Fill(sum.Perp(), fWgt);
	/*hist.ZMet->Fill(fMetBlock->Met(0));
	hist.ZSumet0->Fill(fMetBlock->Sumet(0));
	hist.ZHt->Fill(fMetBlock->Sumet(1));
	*/

	hist.ZMet->Fill(jetMod->GetMyMetCorr(0,3), fWgt);
	hist.ZSumet0->Fill(jetMod->GetMySumEtCorr(0,3), fWgt);
	hist.ZHt->Fill(jetMod->GetMyHtCorr(0,3), fWgt);

}
/*-------------------------------------------------------------------*/
void ZSel::BookEleHist(ElePlots_t& hist, std::string sFoldName,
								std::string sFoldTitle, std::string sTextToAdd)
{
	char title[200];
	char ytitle[200];

	sprintf(title, "%s - E_{T}^{corr}",sTextToAdd.c_str());
	hist.EtCorr   	= new TH1F("EtCorr",title,200,0,200);
	sprintf(title, "%s - Had/Em",sTextToAdd.c_str());
	hist.HadEm		= new TH1F("HadEm",title,40,0,4);
	sprintf(title, "%s - Iso(0.4)",sTextToAdd.c_str());
	hist.IsoEtCorr	= new TH1F("Iso4",title,100,-10,90);
	sprintf(title, "%s - Track P_{T}",sTextToAdd.c_str());
	hist.TrkPt 		= new TH1F("Trkpt",title,200,0,200);
	sprintf(title, "%s - E/P",sTextToAdd.c_str());
	hist.EoverP 	= new TH1F("EoverP",title,50,0,10);
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

	sprintf(title, "%s - #Delta#Phi(e,#slash{E}_{T})",sTextToAdd.c_str());
	hist.DelPhiEleMet		= new TH1F("DelPhiEleMet",title,70,0,3.5);

	hist.EtCorr->GetXaxis()->SetTitle("E_{T} (GeV)");
	hist.HadEm->GetXaxis()->SetTitle("Had/Em");
	hist.IsoEtCorr->GetXaxis()->SetTitle("Isolation E_{T} (GeV)");
	hist.TrkPt->GetXaxis()->SetTitle("Track P_{T} (GeV)");
	hist.EoverP->GetXaxis()->SetTitle("E/P");
	hist.Chi2Mean->GetXaxis()->SetTitle("#Chi^{2} Mean");
	hist.Ces2Wire->GetXaxis()->SetTitle("2nd CES WIRE cluster energy");
	hist.Ces2Strip->GetXaxis()->SetTitle("2nd CES STRIP cluster energy");
	hist.DelPhiEleMet->GetXaxis()->SetTitle("#Delta#Phi(e,#slash{E}_{T})");

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
	new_folder->Add(hist.DelPhiEleMet);	
}
//_____________________________________________________________________________
void ZSel::BookZHist(ZPlots_t& hist, std::string sFoldName,
							std::string sFoldTitle)
{
	hist.ZMass 			= new TH1F("Z_mass","Z Invriant Mass",250,0,250);
	hist.ZPt 			= new TH1F("Z_Pt","Z Pt",250,0,250);
	hist.ZEleEtratio	= new TH1F("Z_eleEtratio","Z electrons Et ratio (2nd leading to leading) ",100,0,2);
	hist.ZSumet0 		= new TH1F("Z_Sumet0","Z Events Sumet(0)",150,0,1500);
	hist.ZHt 			= new TH1F("Z_Ht","Z Events Ht",150,0,1500);
	hist.ZMet 			= new TH1F("Z_Met","Met for Z events",250,0,250);
	hist.ZeleDelPhi   = new TH1F("Z_EleDelPhi","#Delta #Phi of the two electrons of the Z events",700,-7,7);

	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}


	char ytitle[200];
	sprintf(ytitle,"Events/%.3f",hist.ZMass->GetBinWidth(1));
	hist.ZMass->GetYaxis()->SetTitle(ytitle);
	hist.ZMass->GetXaxis()->SetTitle("Invariant Mass (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.ZPt->GetBinWidth(1));
	hist.ZPt->GetYaxis()->SetTitle(ytitle);
	hist.ZPt->GetXaxis()->SetTitle("P_{T} (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.ZEleEtratio->GetBinWidth(1));
	hist.ZEleEtratio->GetYaxis()->SetTitle(ytitle);
	hist.ZEleEtratio->GetXaxis()->SetTitle("#frac{2nd Lead Ele. E_{T}}{Lead Ele. E_{T}}");
	sprintf(ytitle,"Events/%.3f",hist.ZSumet0->GetBinWidth(1));
	hist.ZSumet0->GetYaxis()->SetTitle(ytitle);
	hist.ZSumet0->GetXaxis()->SetTitle("#Sigma E_{T} (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.ZHt->GetBinWidth(1));
	hist.ZHt->GetYaxis()->SetTitle(ytitle);
	hist.ZHt->GetXaxis()->SetTitle("H_{T} (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.ZMet->GetBinWidth(1));
	hist.ZMet->GetYaxis()->SetTitle(ytitle);
	hist.ZMet->GetXaxis()->SetTitle("ME_{T} (GeV)");
	sprintf(ytitle,"Events/%.3f",hist.ZeleDelPhi->GetBinWidth(1));
	hist.ZeleDelPhi->GetYaxis()->SetTitle(ytitle);
	hist.ZeleDelPhi->GetXaxis()->SetTitle("#delta #phi");

	
	new_folder->Add(hist.ZMass);	
	new_folder->Add(hist.ZPt);
	new_folder->Add(hist.ZEleEtratio);	
	new_folder->Add(hist.ZSumet0);
	new_folder->Add(hist.ZHt);
	new_folder->Add(hist.ZMet);
	new_folder->Add(hist.ZeleDelPhi);
}
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void ZSel::BookRunDependentHistograms(RunDep_Hists_t& hist,
				std::string sFoldName, std::string sFoldTitle)
{
	hist.Nevts = new TH1F("Nevts","Events per Run",270000,130000,400000);
	hist.NZs = new TH1F("NZ","Z Events per Run",270000,130000,400000);
	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(hist.Nevts);	
	new_folder->Add(hist.NZs);	
}
/*-------------------------------------------------------------------*/
void ZSel::BookAcceptanceHist(AcceptPlots_t& hist, std::string sFoldName,
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
	hist.Cumm->GetXaxis()->SetBinLabel(18,"Two e's");

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
	hist.All->GetXaxis()->SetBinLabel(18,"Two e's");


	TFolder* new_folder = GetHistFolder(this, sFoldName,sFoldTitle);
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Request for a new Histogram Folder returned NULL. No histograms are added.");
		return;
	}
	new_folder->Add(hist.Cumm);	
	new_folder->Add(hist.All);	
	
}



/*-------------------------------------------------------------------*/
void ZSel::FillDataBlocks(int ientry)
{
	fMetBlock->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
}

/*-------------------------------------------------------------------*/
void ZSel::SetNvtxWeights()
{

	vNvtxWgts.clear();

	if (iDataSample == 20) // Zee MC
	{
		vNvtxWgts.push_back(0.595432);
		vNvtxWgts.push_back(1.07053);
		vNvtxWgts.push_back(1.72818);
		vNvtxWgts.push_back(2.43945);
		vNvtxWgts.push_back(3.35186);
		vNvtxWgts.push_back(5.33692);
		vNvtxWgts.push_back(10.226);
		vNvtxWgts.push_back(17.4119);

	} else if (iDataSample == 21) //Zmm MC
	{
		vNvtxWgts.push_back(0.296238);
	} else if (iDataSample == 22) //Ztt MC
	{
		vNvtxWgts.push_back(0.514841);
		vNvtxWgts.push_back(1.39528);
		vNvtxWgts.push_back(1.74787);
		vNvtxWgts.push_back(2.57032);
		vNvtxWgts.push_back(2.46959);
	} else if (iDataSample == 30) //Wen MC)
	{
		vNvtxWgts.push_back(0.370298);
		vNvtxWgts.push_back(1.7441);
	} else if (iDataSample == 31) //Wmn MC
	{
		//no events
	} else if (iDataSample == 32) //Wtn MC
	{
		//no events
	} else {
		assert (false && "Unknown sample, No Nvtx weights found!");
	}

}


//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int ZSel::EndJob() {

	printf("[ZSL:00:]----- end job: ---- %s\n",GetName());
	if (qMc)	StdOut(__FILE__,__LINE__,3,"This module is designed for data only. But the run is on MC.");
	std::cout << "[ZSL:01:] Events Processed ----------- = " << counter.evtsRunOver << std::endl;
	std::cout << "[ZSL:02:] Events Passed -------------- = " << counter.evtsPassModule << std::endl;
	std::cout << "[ZSL:03:] Z Events ------------------- = " << counter.Zevents << std::endl;
	std::cout << "[ZSL:04:] Min. Ele Et ---------------- = " << GetMinEleEt() << std::endl;
	std::cout << "[ZSL:05:] Max. Ele DetEta ------------ = " << GetMaxDetEta() << std::endl;
	std::cout << "[ZSL:06:] Uset Nvtx Wegihts?  -------- = " << GetUseNvtxWeigts() << " (1==YES)" << std::endl;
	std::cout << "[ZSL:07:] iDataSample ---------------- = " << GetDataSample() << " (20=Zee,21=Zmm,22=Ztt,30=Wen,31=Wmn,32=Wtn)" << std::endl;
	hAccept.Cumm->Print();
	for (int i=0;i < hAccept.Cumm->GetNbinsX() ; i++) std::cout << "Zbin"<< i+1 << "=" << hAccept.Cumm->GetBinContent(i+1) << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
