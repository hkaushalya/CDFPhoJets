/*
READ THIS EVERYDAY!
keep it simple and readable. forget about performance!
*/

#include <iostream>
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnTriggerTable.hh"
#include "samantha/PhoJetsMet/TIDcutsAna.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "Stntuple/ana/TPrintModule.hh"
#include "Stntuple/ana/TMathModule.hh"



/*
function list in order

void Cleanup()
bool PrePass()
bool PassGoodRun()
bool PassTrigger()
bool PassZvCut()
int Class12Vertices()
double BestVetex_Z()

CEMPhoEleIDcut(pho,event, zvx)
TightPhotonIDcut(pho)
LoosePhotonIDcut(pho)

FillDataBlock(ientry)

FillPhotonPlots(SuperPhoton_t)
FillElectronPlots(SuperPhoton_t)

TFolder* GetHistoFolder
BookPhotonHistograms
BookElectronHistograms

*/



ClassImp(TIDcutsAna)

//_____________________________________________________________________________
TIDcutsAna::TIDcutsAna(const char* name, const char* title):
  TStnModule(name,title), qMc(false)
{
	std::cout << "Hello I am TIDcutsAna module" << std::endl;
}

//_____________________________________________________________________________
TIDcutsAna::~TIDcutsAna() {
}

//_____________________________________________________________________________
void TIDcutsAna::SaveHistograms() {
}

//_____________________________________________________________________________
void TIDcutsAna::BookHistograms() {

	DeleteHistograms();
	BookGeneralHistograms();
	BookElectronHistograms();
	BookPhotonHistograms();
}


//_____________________________________________________________________________
int TIDcutsAna::BeginJob()
/*{{{*/
{
	std::cout << " BEGIN JOB " << std::endl;
				// register the data block, why?  to unpack them!

	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fTrackBlock 	= (TStnTrackBlock*)    RegisterDataBlock("TrackBlock","TStnTrackBlock");
	fPhotonBlock 	= (TStnPhotonBlock*)   RegisterDataBlock("PhotonBlock","TStnPhotonBlock");
	fTriggerBlock 	= (TStnTriggerBlock*)  RegisterDataBlock("TriggerBlock","TStnTriggerBlock");
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");
	fVertexBlock 	= (TStnVertexBlock*)   RegisterDataBlock("ZVertexBlock","TStnVertexBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");

	BookHistograms();

	// initialize vars
	counter.evtsPassTrigger = 0;
	counter.evtsPass25Trigger = 0;
	counter.evtsPass50Trigger = 0;
	counter.evtsPass70Trigger = 0;
	counter.evtsPassGoodrun = 0;
	counter.evtsPassVertexCut = 0;
	counter.evtsPassZvCut = 0;

	counter.evtsPassAllTriggers = 0;
	counter.evtsPass2550Triggers = 0;
	counter.evtsPass2570Triggers = 0;
	counter.evtsPass5070Triggers = 0;

	if (goodRunListFile.Length()<=0)
		goodRunListFile = "goodrun_v17_pho_02.txt";

	goodrun.Read(goodRunListFile.Data());

	trigbits.SetTriggerBlock(fTriggerBlock);

  return 0;
}
/*}}}*/

//_____________________________________________________________________________
int TIDcutsAna::BeginRun()
{

 	int 	currRun = fHeaderBlock->RunNumber();
	std::cout << " BEGIN RUN " << currRun << std::endl;
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	trigbits.BeginRun();	
	TStntuple::Init(currRun);
	
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TIDcutsAna::Event(int ientry)
{
	SetPassed(0);
	counter.evtsRunOver++;

	//_____________________________ REQUIRED Thingies 
  FillDataBlocks(ientry);
  
	if (PrePass()) {
		SetPassed(1);
		counter.evtsPassModule++;
		Cleanup();
		PhoLikeEleIDeff();
	}

	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TIDcutsAna::Cleanup()
{
  	trig25 = false;
	trig50 = false;
	trig70 = false;

} //Cleanup


/*-------------------------------------------------------------------*/
void TIDcutsAna::PhoLikeEleIDeff()
{

	double zvx  = 0;
	if (fVertexBlock->GetBestVertex(12,1) != NULL) {
		zvx = fVertexBlock->GetBestVertex(12,1)->Z();
	}
	TStnEvent* event = GetEvent();
	TStntuple::CorrectPhotonEnergy(event,zvx);
	
	for (int i=0; i < fPhotonBlock->NPhotons();i++) {
		int code = CEMPhoEleIDcut(fPhotonBlock->Photon(i), event, zvx);
		if (code==17) std::cout << " \tfound electron" << std::endl;
		for (int j=0; j < code; j++) {
			GeneralPlot.IdCummPlot->Fill(j);
		}
	}
	


}


void TIDcutsAna::FillEleIDcummPlot(unsigned int ID)
{
	for (int i=0 ; i < 17; i++) {
		if ( ( (ID >> i) && 0x1 ) == 0) {
		}
	}


}

/*-------------------------------------------------------------------*/
// Event pre-selection cuts
// REMOVE EVENTS WITHOUT PHOTONS/JETS
// do not cut on number of jets yet. do it in the TGammaJets mod
/*-------------------------------------------------------------------*/
bool TIDcutsAna::PrePass()
{
	if (!qMc) {
		if (fPhotonBlock->NPhotons() < 1) return false;
		else return (PassGoodRun() && PassTrigger() && PassZvCut());
	} else {
		return (PassGoodRun() && PassZvCut());
	}
} // PrePass

/*-------------------------------------------------------------------*/
// good run preselection
// remove the bad runs/sections
/*-------------------------------------------------------------------*/
bool TIDcutsAna::PassGoodRun()
{
  int run = fHeaderBlock->RunNumber();
  int sec = fHeaderBlock->SectionNumber();
	
  bool good = goodrun.Good(run,sec);
	
  if (!good) {
    return false;
  } else {
    counter.evtsPassGoodrun++;
    return true;
  }
} //PassGoodRun

/*-------------------------------------------------------------------*/
// trigger preselection
/*-------------------------------------------------------------------*/
bool TIDcutsAna::PassTrigger()
{
	if (!qMc) { //___________________ pypho22_dfc MC sample does not have trigger info
		trigbits.Event();
    	trig25 = trigbits.Pho25Iso();
	   trig50 = trigbits.Pho50();
   	trig70 = trigbits.Pho70();

	    if (!(trig25 || trig50 || trig70) ) {
   		return false; //____ use all 3 triggers. may see something at hight pt.
	    } else {
      	counter.evtsPassTrigger++;
			if (trig25 && trig50 && trig70) counter.evtsPassAllTriggers++;
			if (trig25 && trig50 && !trig70) counter.evtsPass2550Triggers++;
			if (trig25 && !trig50 && trig70) counter.evtsPass2570Triggers++;
			if (!trig25 && trig50 && trig70) counter.evtsPass5070Triggers++;
			if (trig25 && !trig50 && !trig70) counter.evtsPass25Trigger++;
			if (!trig25 && trig50 && !trig70) counter.evtsPass50Trigger++;
			if (!trig25 && !trig50 && trig70) counter.evtsPass70Trigger++;
			
      	return true;
   	}
	}
}// PassTrigger


/*-------------------------------------------------------------------*/
// Z vertex preselection
// cut on Z position on the vertex,
// so we are within the acceptance of the silicon.
/*-------------------------------------------------------------------*/
bool TIDcutsAna::PassZvCut()
{
	int nvtx = Class12Vertices();
	if (nvtx < 1) {
		return false;
	} else {
		counter.evtsPassVertexCut++;
		double best_z = BestVertex_Z(); 
    	GeneralPlot.vertexZ->Fill(best_z);
		if (fabs(best_z) < 60.) {
	  		counter.evtsPassZvCut++;
			return true;
  		} else {
			return false;
		}
	}
} //PassZvCut

/*-------------------------------------------------------------------*/
// for cut on vetices: require 1 vertex now. can remove some
// background by doing so.
/*-------------------------------------------------------------------*/
int TIDcutsAna::Class12Vertices()
{
 	int nvtx = 0;
	double zvx  = 0;

	int Nvtx = fVertexBlock->NVertices();
	GeneralPlot.Nvertices->Fill(Nvtx);

	for (int ivtx = 0; ivtx < Nvtx; ++ivtx) {
		TStnVertex* vert = fVertexBlock->Vertex(ivtx);
		if (vert->VClass() >= 12) ++nvtx;
	}
  
  return nvtx;
} //Class12Vertices

/*-------------------------------------------------------------------*/
// z-position cut on vertex
/*-------------------------------------------------------------------*/
double TIDcutsAna::BestVertex_Z()
{
  //double sumpt = 0;
  //TStnVertex* bestvert;
  double zvx =0;
  if (fVertexBlock->GetBestVertex(12,1) != NULL) {
    zvx = fVertexBlock->GetBestVertex(12,1)->Z();
  }
  return zvx;
} //BestVertex_Z


/*****************  PHOTON LIKE ID CUT FUNCTION ************************/
//_____ function returns 0 if electrons passes cuts
//      requires zvx of highest Pt class12 vertex
unsigned int TIDcutsAna::CEMPhoEleIDcut(TStnPhoton* Pho, TStnEvent* event, double zvx)
{
	//______________________________________ only CEM photons with Et>7 GeV
	if (Pho->Detector() !=0 ) {
		return 1;
	}
	if (Pho->Etc()< 30.) {
		return 2;
	}

	//______________________________________ HADEM cut using 3 towers
	float cutMin_HADEM=0.0;
	float cutMax_HADEM=0.055+0.00045*(Pho->Momentum()->E());
	if ((Pho->HadEm())<cutMin_HADEM || (Pho->HadEm())>cutMax_HADEM){
		return 3;
	}

	//______________________________________ CES Chi^2 cut
	if ((Pho->Chi2())<0.0 || (Pho->Chi2())>20.0) {
		return 4;
	}

	//______________________________________ N3D cut
	if ((Pho->N3d())==0) {
		return 5; // it is not an electron if N3d=0
	}

	//______________________________________ second N3D cut
	if ((Pho->N3d())>2) {
		return 6;  // no more than 2 tracks (1st--ele,2nd--to match cuts for pho)
	}

	//______________________________________ cut on 2nd max Pt track in cluster if N3D=2
	if ((Pho->N3d())==2) {
		float trkPtcut_max=1.0+0.005*(Pho->Etc());

		if ((Pho->Pt2())>trkPtcut_max) {
			return 7;
		}
	}

	TStnElectron* Ele=TPhotonUtil::MatchElectron(event,Pho); // get matching electron

	//______________________________________ E/p cut for electron
	if (Ele==NULL) {
		return 8;
	}

	if (Ele->TrackNumber()<0) {
		return 9; // there should be a track
	}
	if (fabs(Ele->Z0()-zvx)>3.0) {
		return 10; // only events with ele from best vertex
	}

	int ele_trk = Ele->TrackNumber();
	TStnTrack* Trk = fTrackBlock->Track(ele_trk);

	if ((Trk->Algorithm())==2) {
		if ((Ele->TrackBcPt())<50.0 &&
	 			(Ele->EOverP()<0.8 || Ele->EOverP()>1.2)) return 11;
	} else {
		if ((Ele->TrackPt())<50.0 &&
	 			(Ele->EOverP()<0.8 || Ele->EOverP()>1.2)) return 11;
	}

	//______________________________________ CalIso4 cut
	float cutMax_CalIso=0.0;

	if ((Pho->Etc()) < 20.0) cutMax_CalIso = 0.1 * (Pho->Etc());
	else cutMax_CalIso = 2.0 + 0.02 * ( Pho->Etc() - 20.0 );

	if ((Pho->EIso4(2)) > cutMax_CalIso) return 12;

	//______________________________________ TrkIso4 cut
	float trkIso=Pho->SumPt4() - Pho->Pt();

	if (trkIso < 0.0) trkIso = Pho->SumPt4();

	if (trkIso < (trkIso > (2.0 + 0.005 * (Pho->Etc()) ) ) ) return 13;

	//______________________________________ Energy of 2nd CES cluster (Wire & Strip)
	float _Et2ndCES=0.0;
	if ((Pho->CesStripE2()) > (Pho->CesWireE2())) _Et2ndCES = (Pho->CesStripE2()) * (Pho->SinTheta());
	else _Et2ndCES = (Pho->CesWireE2()) * (Pho->SinTheta());

	float cutMax_2ndCes = 0.0;
	if ((Pho->Etc()) < 18.0) cutMax_2ndCes = 0.14 * (Pho->Etc());
	else cutMax_2ndCes = 2.4 + 0.01 * (Pho->Etc());

	if (_Et2ndCES > cutMax_2ndCes) return 14;

	//_______________________________________ Fiducial cuts
	if (fabs(Pho->XCes()) > 21.0) return 15;
	if (fabs(Pho->ZCes()) < 9.0 || fabs(Pho->ZCes()) > 230.0) return 16;

	return 17;

} // CEMPhoEleIDcut


/*-------------------------------------------------------------------*/
void TIDcutsAna::FillDataBlocks(int ientry)
{
	fTriggerBlock->GetEntry(ientry);
	fPhotonBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
	fVertexBlock->GetEntry(ientry);
	fTrackBlock->GetEntry(ientry);
	fGenpBlock->GetEntry(ientry);
}

/*-------------------------------------------------------------------*/
// fill photons plots
/*-------------------------------------------------------------------*/
void TIDcutsAna::FillPhotonPlots(const TSuperPhoton& sp)
{
	PhotonPlot.Detector->Fill(sp.Detector);
	PhotonPlot.EtCorr->Fill(sp.EtCorr);
	PhotonPlot.XCes->Fill(sp.XCes);
	PhotonPlot.ZCes->Fill(sp.ZCes);
	PhotonPlot.HadEm->Fill(sp.HadEm);
	PhotonPlot.IsoEtCorr->Fill(sp.IsoEtCorr);
	PhotonPlot.Chi2Mean->Fill(sp.Chi2Mean);
	PhotonPlot.N3d->Fill(sp.N3d);
	PhotonPlot.TrkPt->Fill(sp.TrkPt);
	PhotonPlot.TrkIso->Fill(sp.TrkIso);
	PhotonPlot.Ces2Wire->Fill(sp.CesWireE2);
	PhotonPlot.HadTDCtime->Fill(sp.HadTDCtime);
	PhotonPlot.Ces2Strip->Fill(sp.CesStripE2);
}
	
/*-------------------------------------------------------------------*/
// fill electron plots
/*-------------------------------------------------------------------*/
void TIDcutsAna::FillElectronPlots(const TSuperPhoton& sp)
{
	ElectronPlot.Detector->Fill(sp.Detector);
	ElectronPlot.EtCorr->Fill(sp.EtCorr);
	ElectronPlot.XCes->Fill(sp.XCes);
	ElectronPlot.ZCes->Fill(sp.ZCes);
	ElectronPlot.HadEm->Fill(sp.HadEm);
	ElectronPlot.IsoEtCorr->Fill(sp.IsoEtCorr);
	ElectronPlot.Chi2Mean->Fill(sp.Chi2Mean);
	ElectronPlot.N3d->Fill(sp.N3d);
	ElectronPlot.TrkPt->Fill(sp.TrkPt);
	ElectronPlot.TrkIso->Fill(sp.TrkIso);
	ElectronPlot.Ces2Wire->Fill(sp.CesWireE2);
	ElectronPlot.HadTDCtime->Fill(sp.HadTDCtime);
	ElectronPlot.Ces2Strip->Fill(sp.CesStripE2);
}

/*===================================================================*/
//	Create a folder in"Hist" 
/*===================================================================*/
TFolder* TIDcutsAna::GetHistoFolder(char *name, char* title)
{
	char folder_name[200];
	char folder_title[200];
	TFolder* hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
	sprintf(folder_name,name);
	sprintf(folder_title,title);
	TFolder* new_folder = (TFolder*) hist_folder->FindObject(folder_name);
	if (! new_folder) new_folder = hist_folder->AddFolder(folder_name,folder_title,NULL);
	return new_folder;
}
//_____________________________________________________________________________
void TIDcutsAna::BookGeneralHistograms()
{
	std::cout << "Booking General Plots...";
	
	GeneralPlot.IdCummPlot = new TH1F("PhoEleCumm","Cummulative Photon Like Electron ID Cuts",22,0,22);
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(1,"PhoCentral");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(2,"PhoEtc30");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(3,"PhoHadEm");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(4,"PhoChi2");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(5,"PhoN3d=0");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(6,"PhoN3d>2");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(7,"PhoPt2");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(8,"FindMatchEle");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(9,"EleHasTrk");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(10,"EleTrkZ");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(11,"EleEoveP");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(12,"PhoCalIso4");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(13,"PhoTrkIso");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(14,"PhoCes2Et");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(15,"PhoXCes");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(16,"PhoZCes");
	GeneralPlot.IdCummPlot->GetXaxis()->SetBinLabel(17,"Electrons");
	
	GeneralPlot.Nvertices = new TH1F("Nvertices","Number of Class12 vertices per event",10,0,10);
	GeneralPlot.vertexZ = new TH1F("BestVertex_Z","Best Class12 vertex Z position",100,-100,100);

	TFolder* new_folder = GetHistoFolder("General_Plots","General Plots");
	new_folder->Add(GeneralPlot.IdCummPlot);	
	new_folder->Add(GeneralPlot.pho_passIDs);
	new_folder->Add(GeneralPlot.Nvertices);
	new_folder->Add(GeneralPlot.vertexZ);
	std::cout << "DONE\n";

} //BookGeneralHistograms

//_____________________________________________________________________________
void TIDcutsAna::BookPhotonHistograms()
{
	std::cout << "Booking Photon Plots...";

	PhotonPlot.Detector 	= new TH1F("Detector","Wenu MC: Tight Photon Detector Region",3,0,3);
	PhotonPlot.EtCorr 	= new TH1F("EtCorr","Wenu MC: Tight Photon Corrected Et",600,0,120);
	PhotonPlot.XCes 		= new TH1F("XCes","Wenu MC: Tight Photon XCes",640,-32,32);
	PhotonPlot.ZCes 		= new TH1F("ZCes","Wenu MC: Tight Photon ZCes",2000,-250,250);
	PhotonPlot.HadEm 		= new TH1F("HadEm","Wenu MC: Tight Photon HadEm",50,0,0.5);
	PhotonPlot.IsoEtCorr = new TH1F("IsoEtCorr","Wenu MC: Tight Photon IsoCorr",1500,-5,10);
	PhotonPlot.Chi2Mean 	= new TH1F("Chi2Mean","Wenu MC: Tight Photon Chi2Mean (Wire+Strip/2)",800,0,80);
	PhotonPlot.N3d 		= new TH1F("N3d","Wenu MC: Tight Photon N3d",10,0,10);
	PhotonPlot.TrkPt 		= new TH1F("Trkpt","Wenu MC: Tight Photon TrkPt",1000,0,100);
	PhotonPlot.TrkIso 	= new TH1F("TrkIso","Wenu MC: Tight Photon TrkIso",150,0,15);
	PhotonPlot.Ces2Wire 	= new TH1F("Ces2Wire","Wenu MC: Tight Photon CES(2nd) Wire",400,0,40);
	PhotonPlot.Ces2Strip = new TH1F("Ces2Strip","Wenu MC: Tight Photon CES(2nd) Strip",400,0,40);
	PhotonPlot.HadTDCtime = new TH1F("HadTDCtime","Wenu MC: Tight Photon Had TDC time",200,-10,10);

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

	TFolder* new_folder = GetHistoFolder("Photon_Plots","Photon Plots");
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
	new_folder->Add(PhotonPlot.HadTDCtime);

	std::cout << "DONE.\n";

} //BookPhotonHistograms


//_____________________________________________________________________________
void TIDcutsAna::BookElectronHistograms()
{
	std::cout << "Booking Electron Plots...";

	ElectronPlot.Detector 	= new TH1F("Detector","Wenu MC:  Electron Detector Region",3,0,3);
	ElectronPlot.EtCorr 	= new TH1F("EtCorr","Wenu MC:  Electron Corrected Et",600,0,120);
	ElectronPlot.XCes 		= new TH1F("XCes","Wenu MC:  Electron XCes",640,-32,32);
	ElectronPlot.ZCes 		= new TH1F("ZCes","Wenu MC:  Electron ZCes",2000,-250,250);
	ElectronPlot.HadEm 		= new TH1F("HadEm","Wenu MC:  Electron HadEm",50,0,0.5);
	ElectronPlot.IsoEtCorr = new TH1F("IsoEtCorr","Wenu MC:  Electron IsoCorr",1500,-5,10);
	ElectronPlot.Chi2Mean 	= new TH1F("Chi2Mean","Wenu MC:  Electron Chi2Mean (Wire+Strip/2)",800,0,80);
	ElectronPlot.N3d 		= new TH1F("N3d","Wenu MC:  Electron N3d",10,0,10);
	ElectronPlot.TrkPt 		= new TH1F("Trkpt","Wenu MC:  Electron TrkPt",1000,0,100);
	ElectronPlot.TrkIso 	= new TH1F("TrkIso","Wenu MC:  Electron TrkIso",150,0,150);
	ElectronPlot.Ces2Wire 	= new TH1F("Ces2Wire","Wenu MC:  Electron CES(2nd) Wire",400,0,40);
	ElectronPlot.Ces2Strip = new TH1F("Ces2Strip","Wenu MC:  Electron CES(2nd) Strip",400,0,40);
	ElectronPlot.HadTDCtime = new TH1F("HadTDCtime","Wenu MC:  Electron Had TDC time",200,-10,10);

	ElectronPlot.EtCorr->SetYTitle("Events/0.2");
	ElectronPlot.XCes->SetYTitle("Events/0.1");
	ElectronPlot.ZCes->SetYTitle("Events/0.25");
	ElectronPlot.HadEm->SetYTitle("Events/0.01");
	ElectronPlot.IsoEtCorr->SetYTitle("Events/0.01");
	ElectronPlot.Chi2Mean->SetYTitle("Events/0.1");
	ElectronPlot.N3d->SetYTitle("Events/1.0");
	ElectronPlot.TrkPt->SetYTitle("Events/0.1");
	ElectronPlot.TrkIso->SetYTitle("Events/0.1");
	ElectronPlot.Ces2Wire->SetYTitle("Events/0.1");
	ElectronPlot.Ces2Strip->SetYTitle("Events/0.1");
	ElectronPlot.HadTDCtime->SetYTitle("Events/0.1");

	TFolder* new_folder = GetHistoFolder("Electron_Plots","Electron Plots");
	new_folder->Add(ElectronPlot.Detector);
	new_folder->Add(ElectronPlot.EtCorr);
	new_folder->Add(ElectronPlot.XCes);
	new_folder->Add(ElectronPlot.ZCes);
	new_folder->Add(ElectronPlot.HadEm);
	new_folder->Add(ElectronPlot.IsoEtCorr);
	new_folder->Add(ElectronPlot.Chi2Mean);
	new_folder->Add(ElectronPlot.N3d);
	new_folder->Add(ElectronPlot.TrkPt);
	new_folder->Add(ElectronPlot.TrkIso);
	new_folder->Add(ElectronPlot.Ces2Wire);
	new_folder->Add(ElectronPlot.Ces2Strip);
	new_folder->Add(ElectronPlot.HadTDCtime);

	std::cout << "DONE.\n";

} //BookElectronHistograms



//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TIDcutsAna::EndJob() {

	printf("----- end job: ---- %s\n",GetName());
	if (qMc)	std::cout << "RUN IS ON MC SAMPLE" << std::endl;
	else	std::cout << "RUN IS ON DATA SAMPLE" << std::endl;
	std::cout << "EVENTS RUN OVER ------------ = " << counter.evtsRunOver << std::endl;
	std::cout << "EVENTS Pass Goodrun -------- = " << counter.evtsPassGoodrun << std::endl;
	std::cout << "EVENTS Pass Trigger -------- = " << counter.evtsPassTrigger << "\t(only 25,50,70 = " << counter.evtsPass25Trigger
																	<<","<< counter.evtsPass50Trigger << "," << counter.evtsPass70Trigger 
																	<<") (25&50, 25&70, 50&70, all = "
																	<< counter.evtsPass2550Triggers << "," << counter.evtsPass2570Triggers << ","
																	<< counter.evtsPass5070Triggers << "," << counter.evtsPassAllTriggers << ")"
																	<< std::endl;
	std::cout << "EVENTS Pass Vertex Cut ----- = " << counter.evtsPassVertexCut << std::endl;
	std::cout << "EVENTS Pass Zv Cut  -------- = " << counter.evtsPassZvCut << std::endl;
	std::cout << "EVENTS PASS THIS MODULE ---- = " << counter.evtsPassModule << std::endl;

	printf("---------------------------------------------------\n");
	return 0;
}
