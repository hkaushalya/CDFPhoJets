///////////////////////////////////////////////////////////
// Defines all the common histograms. Any moddule can    //
// use them easily and quickly.                          //
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>

#include "samantha/obj/Histograms.hh"
#include <fstream>
#include <iostream>


//==========================================================
Histograms::Histograms():
sTextToAdd("")
{
	// do not make anything until requested.
		//photon hists
		//electron hists
		//jet hists
		//two obj hists
		//three obj hists

}

//==========================================================
Histograms::~Histograms()
{
}

//==========================================================
void Histograms::GetEmObjHistograms(TFolder* new_folder,
						EmObjHists_t& hEmObj, std::string sText)
{
	if (!new_folder) {
		StdOut(__FILE__,__LINE__,3,"Null Histogram Folder. No histograms are added.");
	} else {
		sTextToAdd = sText;
		sprintf(title, "%s - Detector Region",sTextToAdd.c_str());
		hEmObj.Detector 	= new TH1F("Detector",title,3,0,3);
		hEmObj.Detector->GetXaxis()->SetBinLabel(1,"Central");
		hEmObj.Detector->GetXaxis()->SetBinLabel(2,"Plug");
		hEmObj.Detector->GetXaxis()->SetBinLabel(3,"Forward");
		sprintf(title, "%s - Detector #eta",sTextToAdd.c_str());
		hEmObj.DetEta 	= new TH1F("DetEta",title,2000,-10,10);
		sprintf(title, "%s - Detector #Phi",sTextToAdd.c_str());
		hEmObj.DetPhi 	= new TH1F("DetPhi",title,140,-7,7);
		sprintf(title, "%s - E_{T}",sTextToAdd.c_str());
		hEmObj.EtCorr 	= new TH1F("EtCorr",title,500,0,200);
		sprintf(title, "%s - XCes",sTextToAdd.c_str());
		hEmObj.XCes 		= new TH1F("XCes",title,640,-32,32);
		sprintf(title, "%s - ZCes",sTextToAdd.c_str());
		hEmObj.ZCes 		= new TH1F("ZCes",title,2000,-250,250);
		sprintf(title, "%s - HadEm",sTextToAdd.c_str());
		hEmObj.HadEm 		= new TH1F("HadEm",title,400,0,2.);
		sprintf(title, "%s - #chi^{2} Mean",sTextToAdd.c_str());
		hEmObj.Chi2Mean 	= new TH1F("Chi2Mean",title,1210,-1,120);
		sprintf(title, "%s - N3d",sTextToAdd.c_str());
		hEmObj.N3d 		= new TH1F("N3d",title,50,0,50);
		sprintf(title, "%s - Iso",sTextToAdd.c_str());
		hEmObj.Iso4 		= new TH1F("IsoEtCorr",title,3000,-5,25);
		sprintf(title, "%s - TrkPt",sTextToAdd.c_str());
		hEmObj.TrkPt 		= new TH1F("Trkpt",title,1000,0,100);
		sprintf(title, "%s - TrkIso",sTextToAdd.c_str());
		hEmObj.TrkIso 	= new TH1F("TrkIso",title,150,0,15);
		sprintf(title, "%s - CES(2nd) Wire",sTextToAdd.c_str());
		hEmObj.Ces2Wire 	= new TH1F("Ces2Wire",title,400,-40,40);
		sprintf(title, "%s - CES(2nd) Strip",sTextToAdd.c_str());
		hEmObj.Ces2Strip = new TH1F("Ces2Strip",title,400,-40,40);
		sprintf(title, "%s - EM Time ",sTextToAdd.c_str());
		hEmObj.EmTime 	= new TH1F("EmTime",title,300,-100,200);
		sprintf(title, "%s - E/p",sTextToAdd.c_str());
		hEmObj.EoverP 	= new TH1F("EoverP",title,110,-1,10);
		sprintf(title, "%s - Track Z0",sTextToAdd.c_str());
		hEmObj.TrkZ0 	= new TH1F("TrkZ0",title,500,-500,500);

		sprintf(title, "%s - EM Timing of #gamma",sTextToAdd.c_str());
		hEmObj.EmTimeVsRun 	 = new TProfile("EmTimeVsRun",title,270000,130000,400000,-200,200);
		sprintf(title, "%s - E_{T}^{corr} of #gamma",sTextToAdd.c_str());
		hEmObj.EtCorrVsRun 	 = new TProfile("EtCorrVsRun",title,270000,130000,400000,-200,200);

		sprintf(ytitle,"Events/%.2f",hEmObj.DetEta->GetBinWidth(1));
		hEmObj.DetEta->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.DetPhi->GetBinWidth(1));
		hEmObj.DetPhi->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.EtCorr->GetBinWidth(1));
		hEmObj.EtCorr->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.XCes->GetBinWidth(1));
		hEmObj.XCes->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.ZCes->GetBinWidth(1));
		hEmObj.ZCes->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.HadEm->GetBinWidth(1));
		hEmObj.HadEm->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.Chi2Mean->GetBinWidth(1));
		hEmObj.Chi2Mean->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.N3d->GetBinWidth(1));
		hEmObj.N3d->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.Iso4->GetBinWidth(1));
		hEmObj.Iso4->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.TrkPt->GetBinWidth(1));
		hEmObj.TrkPt->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.TrkIso->GetBinWidth(1));
		hEmObj.TrkIso->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.Ces2Wire->GetBinWidth(1));
		hEmObj.Ces2Wire->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.Ces2Strip->GetBinWidth(1));
		hEmObj.Ces2Strip->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.EmTime->GetBinWidth(1));
		hEmObj.EmTime->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.EoverP->GetBinWidth(1));
		hEmObj.EoverP->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hEmObj.TrkZ0->GetBinWidth(1));
		hEmObj.TrkZ0->SetYTitle(ytitle);
		
		new_folder->Add(hEmObj.Detector);
		new_folder->Add(hEmObj.DetEta);
		new_folder->Add(hEmObj.DetPhi);
		new_folder->Add(hEmObj.EtCorr);
		new_folder->Add(hEmObj.XCes);
		new_folder->Add(hEmObj.ZCes);
		new_folder->Add(hEmObj.HadEm);
		new_folder->Add(hEmObj.Chi2Mean);
		new_folder->Add(hEmObj.N3d);
		new_folder->Add(hEmObj.Iso4);
		new_folder->Add(hEmObj.TrkPt);
		new_folder->Add(hEmObj.TrkIso);
		new_folder->Add(hEmObj.Ces2Wire);
		new_folder->Add(hEmObj.Ces2Strip);
		new_folder->Add(hEmObj.EmTime);
		new_folder->Add(hEmObj.EoverP);
		new_folder->Add(hEmObj.TrkZ0);
		new_folder->Add(hEmObj.EmTimeVsRun);
		new_folder->Add(hEmObj.EtCorrVsRun);
	}
	
}


//==========================================================
void Histograms::GetEventHistograms(TFolder* new_folder,
									EventHists_t& Hist, std::string sText)
{
	// Histograms of Event variables like MEt, SumEt etc.

	
	sTextToAdd = sText;
	if (!new_folder) {
		StdOut(__FILE__,__LINE__,3,"Null Histogram Folder. No histograms are added.");
		return;
	} else {
		sprintf(title, "%s - ME_{T}",sTextToAdd.c_str());
		Hist.Met 			 = new TH1F("Met",title,1500,0,1500);
		
		sprintf(title, "%s - SumE_{T}",sTextToAdd.c_str());
		Hist.Sumet 		 = new TH1F("Sumet",title,1500,0,1500);
		
		sprintf(title, "%s - H_{T}",sTextToAdd.c_str());
		Hist.Ht 		 = new TH1F("Ht",title,1500,0,1500);
		
		sprintf(title, "%s - # of Vertices",sTextToAdd.c_str());
		Hist.NVertices 	 = new TH1F("NVertices",title,25,0,25);
		
		sprintf(title, "%s - # of Class 12 Vertices",sTextToAdd.c_str());
		Hist.N12Vertices = new TH1F("N12Vertices",title,25,0,25);
		
		sprintf(title, "%s - # of Jets(>15GeV)",sTextToAdd.c_str());
		Hist.NJet15		 = new TH1F("NJet15",title,50,0,50);
		
		sprintf(title, "%s - Average ME_{T}",sTextToAdd.c_str());
		Hist.MetRun      = new TProfile("MEtVsRun",title,270000,130000,400000,0,1500);
		
		sprintf(title, "%s - Average SumE_{T}",sTextToAdd.c_str());
		Hist.SumetRun 	 = new TProfile("SumetVsRun",title,270000,130000,400000,0,1500);
		
		sprintf(title, "%s - Average # of Vertices",sTextToAdd.c_str());
		Hist.NVerticesRun 	 = new TProfile("VerticesVsRun",title,270000,130000,400000,0,25);
		
		sprintf(title, "%s - Average # of Class 12 Vertices",sTextToAdd.c_str());
		Hist.N12VerticesRun = new TProfile("N12VerticesVsRun",title,270000,130000,400000,0,25);
		
		sprintf(title, "%s - Average # of Jets (After removing EM objects)",sTextToAdd.c_str());
		Hist.NJet15Run 	 = new TProfile("NJetsVsRun",title,270000,130000,400000,0,50);

		sprintf(title, "%s - NJets(15GeV) Vs ME_{T}",sTextToAdd.c_str());
		Hist.NJet15VsMet      = new TH2F("NJet15VsMet",title,1500,0,1500,50,0,50);
		sprintf(title, "%s - NJets(15GeV) Vs SumE_{T}",sTextToAdd.c_str());
		Hist.NJet15VsSumet 	 = new TH2F("NJet15VsSumet",title,1500,0,1500,50,0,50);
		sprintf(title, "%s - # of Class 12 Vertices Vs ME_{T}",sTextToAdd.c_str());
		Hist.N12VerticesVsMet 	 = new TH2F("N12VerticesVsMet",title,1500,0,1500,50,0,50);
		sprintf(title, "%s - # of Class 12 Vertices Vs SumE_{T}",sTextToAdd.c_str());
		Hist.N12VerticesVsSumet = new TH2F("N12VerticesVsSumet",title,1500,0,1500,50,0,50);

		sprintf(title, "%s - ME_{T}/SumE_{T}",sTextToAdd.c_str());
		Hist.MetSumetRatio 	 = new TH1F("MetSumetRatio",title,50,0,5);

		sprintf(title,"Best Class 12 Vertex Z coordinate");
		Hist.BestVertexZ = new TH1F("VertexZ",title,200,-200,200);


		new_folder->Add(Hist.Met);
		new_folder->Add(Hist.Sumet);
		new_folder->Add(Hist.Ht);
		new_folder->Add(Hist.NVertices);
		new_folder->Add(Hist.N12Vertices);
		new_folder->Add(Hist.NJet15);
		new_folder->Add(Hist.MetRun);
		new_folder->Add(Hist.SumetRun);
		new_folder->Add(Hist.NVerticesRun);
		new_folder->Add(Hist.N12VerticesRun);
		new_folder->Add(Hist.NJet15Run);
		new_folder->Add(Hist.NJet15VsMet);
		new_folder->Add(Hist.NJet15VsSumet);
		new_folder->Add(Hist.N12VerticesVsMet);
		new_folder->Add(Hist.N12VerticesVsSumet);
		new_folder->Add(Hist.MetSumetRatio);
		new_folder->Add(Hist.BestVertexZ);
	}
}


//==========================================================
void Histograms::GetTwoObjHistograms(TFolder* new_folder,
							TwoObjHists_t& Hist, std::string sText)
{
	sTextToAdd = sText;
	if (!new_folder) {
		StdOut(__FILE__,__LINE__,3,"Null Histogram Folder. No histograms are added.");
		return;
	} else {
		sprintf(title, "%s - Invariant Mass",sTextToAdd.c_str());
		Hist.InvMass	 = new TH1F("InvMass",title,1000,0,1000);
		sprintf(title, "%s - #Delta#Phi",sTextToAdd.c_str());
		Hist.DelPhi		 = new TH1F("DelPhi",title,1400,-7,7);
		sprintf(title, "%s - #Delta#eta",sTextToAdd.c_str());
		Hist.DelEta		 = new TH1F("DelEta",title,100,-5,5);
		sprintf(title, "%s - #DeltaR=#sqrt{#Delta #Phi^{2}+#Delta #eta^{2}}",sTextToAdd.c_str());
		Hist.DelR		 = new TH1F("DelR",title,1000,0,10);
		sprintf(title, "%s - Average Invaraint Mass",sTextToAdd.c_str());
		Hist.InvMassRun = new TProfile("InvMassVsRun",title,270000,130000,400000,0,1000);
		sprintf(title, "%s - Average #Delta#Phi",sTextToAdd.c_str());
		Hist.DelPhiRun	 = new TProfile("DelPhiVsRun",title,270000,130000,400000,-7,7);
		sprintf(title, "%s - Average #Delta#eta",sTextToAdd.c_str());
		Hist.DelEtaRun	 = new TProfile("DelEtaVsRun",title,270000,130000,400000,-5,5);
		sprintf(title, "%s - Average #DeltaR=#sqrt{#Delta#phi^{2}+#Delta#eta^{2}",sTextToAdd.c_str());
		Hist.DelRRun	 = new TProfile("DelRVsRun",title,270000,130000,400000,0,20);

		sprintf(title, "%s - E_{T} ratio",sTextToAdd.c_str());
		Hist.EtRatio	 = new TH1F("EtRatio",title,1000,0,10);
		sprintf(title, "%s - Average E_{T} ratio",sTextToAdd.c_str());
		Hist.EtRatioRun	 = new TProfile("EtRatioVsRun",title,270000,130000,400000,0,20);
		

		new_folder->Add(Hist.InvMass);
		new_folder->Add(Hist.DelPhi);
		new_folder->Add(Hist.DelEta);
		new_folder->Add(Hist.DelR);
		new_folder->Add(Hist.InvMassRun);
		new_folder->Add(Hist.DelPhiRun);
		new_folder->Add(Hist.DelEtaRun);
		new_folder->Add(Hist.DelRRun);
		new_folder->Add(Hist.EtRatio);
		new_folder->Add(Hist.EtRatioRun);
	}
}


//==========================================================
void Histograms::GetThreeObjHistograms(TFolder* new_folder,
							ThreeObjHists_t& Hist, std::string sText)
{
	sTextToAdd = sText;
	if (!new_folder) {
		StdOut(__FILE__,__LINE__,3,"Null Histogram Folder. No histograms are added.");
		return;
	} else {
		sprintf(title, "%s - Invariant Mass",sTextToAdd.c_str());
		Hist.InvMass	 = new TH1F("InvMass",title,1000,0,1000);
		sprintf(title, "%s - Average Invaraint Mass",sTextToAdd.c_str());
		Hist.InvMassRun = new TProfile("InvMassVsRun",title,270000,130000,400000,0,1000);

		new_folder->Add(Hist.InvMass);
		new_folder->Add(Hist.InvMassRun);
	}

}


//==========================================================
void Histograms::GetHaloHistograms(TFolder* new_folder, BeamHaloHists_t& Hist, 
						std::string sText)
{
	if (!new_folder) {
		StdErr(__FILE__,__LINE__,3,"Pl. specify a valid folder. No histograms are added.");
		return;
	}

  char name [200];
  char title[200];

  sprintf(name,"HaloSeedWedgeEM_SeedWedgeHAD");
  sprintf(title,"%s - N_{EMtwr} vs. N_{HADtwr} in seed wedge", sText.c_str());
  Hist.fHaloSeedWedgeEM_SeedWedgeHAD = new TH2F(name,title,25,-0.5,24.5,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloSeedWedgeEM_SeedWedgeHAD);

  sprintf(name,"HaloSeedWedge");
  sprintf(title,"%s - number of EM towers in seed wedge, Max's version", sText.c_str());
  Hist.fHaloSeedWedge=new TH1F(name,title,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloSeedWedge);
  sprintf(name,"HaloSideWedge");
  sprintf(title,"%s - number of EM towers in side wedges, Max's version", sText.c_str());
  Hist.fHaloSideWedge=new TH1F(name,title,45,-0.5,44.5);
  new_folder->Add(Hist.fHaloSideWedge);
  sprintf(name,"HaloSeedWedgeH");
  sprintf(title,"%s - number of HAD towers in seed wedge, my version", sText.c_str());
  Hist.fHaloSeedWedgeH=new TH1F(name,title,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloSeedWedgeH);
  sprintf(name,"HaloSideWedgeH");
  sprintf(title,"%s - number of HAD towers in side wedges, my version", sText.c_str());
  Hist.fHaloSideWedgeH=new TH1F(name,title,45,-0.5,44.5);
  new_folder->Add(Hist.fHaloSideWedgeH);
  sprintf(name,"HaloEastNHad");
  sprintf(title,"%s - number of HAD towers on East side", sText.c_str());
  Hist.fHaloEastNHad=new TH1F(name,title,35,-0.5,34.5);
  new_folder->Add(Hist.fHaloEastNHad);
  sprintf(name,"HaloWestNHad");
  sprintf(title,"%s - number of HAD towers on West side", sText.c_str());
  Hist.fHaloWestNHad=new TH1F(name,title,35,-0.5,34.5);
  new_folder->Add(Hist.fHaloWestNHad);
  sprintf(name,"HaloNHad");
  sprintf(title,"%s - number of HAD towers on East+West", sText.c_str());
  Hist.fHaloNHad=new TH1F(name,title,70,-0.5,69.5);
  new_folder->Add(Hist.fHaloNHad);
  sprintf(name,"HaloEmSeedE");
  sprintf(title,"%s - EM energy of seed wedge towers", sText.c_str());
  Hist.fHaloEmSeedE=new TH1F(name,title,1000,0.0,1000.0);
  new_folder->Add(Hist.fHaloEmSeedE);
  sprintf(name,"HaloEmSideE");
  sprintf(title,"%s - EM energy of side wedge towers", sText.c_str());
  Hist.fHaloEmSideE=new TH1F(name,title,1000,0.0,1000.0);
  new_folder->Add(Hist.fHaloEmSideE);
  sprintf(name,"HaloHadSeedE");
  sprintf(title,"%s - Plug HAD energy of seed wedge towers", sText.c_str());
  Hist.fHaloHadSeedE=new TH1F(name,title,1000,0.0,1000.0);
  new_folder->Add(Hist.fHaloHadSeedE);
  sprintf(name,"HaloSeedWedgeHadE");
  sprintf(title,"%s - HAD energy of seed wedge towers", sText.c_str());
  Hist.fHaloSeedWedgeHadE=new TH1F(name,title,5000,0.0,100.0);
  new_folder->Add(Hist.fHaloSeedWedgeHadE);
  sprintf(name,"HaloEast");
  sprintf(title,"%s - HaloEast", sText.c_str());
  Hist.fHaloEast=new TH1F(name,title,10,-0.5,9.5);
  new_folder->Add(Hist.fHaloEast);
  sprintf(name,"HaloWest");
  sprintf(title,"%s - HaloWest", sText.c_str());
  Hist.fHaloWest=new TH1F(name,title,10,-0.5,9.5);
  new_folder->Add(Hist.fHaloWest);
  sprintf(name,"HaloEastWest");
  sprintf(title,"%s - HaloWest+HaloEast", sText.c_str());
  Hist.fHaloEastWest=new TH1F(name,title,10,-0.5,9.5);
  new_folder->Add(Hist.fHaloEastWest);
  sprintf(name,"HaloSeedWedgeR");
  sprintf(title,"%s - number of towers in seed wedge, Ray's version", sText.c_str());
  Hist.fHaloSeedWedgeR=new TH1F(name,title,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloSeedWedgeR);
  sprintf(name,"HaloSideWedgeR");
  sprintf(title,"%s - number of towers in side wedges, Ray's version", sText.c_str());
  Hist.fHaloSideWedgeR=new TH1F(name,title,45,-0.5,44.5);
  new_folder->Add(Hist.fHaloSideWedgeR);
  sprintf(name,"HaloTwrInRawR");
  sprintf(title,"%s - number of continues towers in seed wedge, Ray's version", sText.c_str());
  Hist.fHaloTwrInRawR=new TH1F(name,title,25,-0.5,24.5);
  new_folder->Add(Hist.fHaloTwrInRawR);

  sprintf(name,"HaloCesStripEoverE");
  sprintf(title,"%s - strip energy E_{CES}/E_{#gamma}", sText.c_str());
  Hist.fHaloCesStripEoverE=new TH1F(name,title,500,0.0,5.0);
  new_folder->Add(Hist.fHaloCesStripEoverE);
  sprintf(name,"HaloCesWireEoverE");
  sprintf(title,"%s - wire energy E_{CES}/E_{#gamma}", sText.c_str());
  Hist.fHaloCesWireEoverE=new TH1F(name,title,500,0.0,5.0);
  new_folder->Add(Hist.fHaloCesWireEoverE);
  sprintf(name,"HaloHadTDC");
  sprintf(title,"%s - photons Had TDC", sText.c_str());
  Hist.fHaloHadTDC=new TH1F(name,title,800,-100.0,100.0);
  new_folder->Add(Hist.fHaloHadTDC);

  sprintf(name,"HaloPhiWedge");
  sprintf(title,"%s - wedge number of halo candidate", sText.c_str());
  Hist.fHaloPhiWedge=new TH1F(name,title,24,-0.5,23.5);
  new_folder->Add(Hist.fHaloPhiWedge);
  sprintf(name,"HaloEta");
  sprintf(title,"%s - #eta of halo candidate", sText.c_str());
  Hist.fHaloEta=new TH1F(name,title,120,-3.0,3.0);
  new_folder->Add(Hist.fHaloEta);
  sprintf(name,"HaloIso");
  sprintf(title,"%s -  CalIso of halo candidate", sText.c_str());
  Hist.fHaloIso= new TH1F(name,title,140,-4.0,10.0);
  new_folder->Add(Hist.fHaloIso);
  sprintf(name,"HaloHadEm");
  sprintf(title,"%s - Had/Em of halo candidate", sText.c_str());
  Hist.fHaloHadEm= new TH1F(name,title,400,0,2);
  new_folder->Add(Hist.fHaloHadEm);
  sprintf(name,"HaloCesChi2");
  sprintf(title,"%s - CES(strip+wire) #chi^{2} of halo candidate", sText.c_str());
  Hist.fHaloCesChi2= new TH1F(name,title,500,0.0,100.0);
  new_folder->Add(Hist.fHaloCesChi2);
  sprintf(name,"HaloEt");
  sprintf(title,"%s -  corrected E_{T} of halo candidate", sText.c_str());
  Hist.fHaloEt= new TH1F(name,title,1000,0.0,1000.0);
  new_folder->Add(Hist.fHaloEt);


	// my stuff

  sprintf(name,"HaloEvt_SeedWedge_Nhad");
  sprintf(title,"%s - Em towers in seed wedge and N Had towers (eat+west)", sText.c_str());
  Hist.fSeedWedge_HadTowers = new TH2F(name,title,20,0,20,20,0,20);
  Hist.fSeedWedge_HadTowers->GetXaxis()->SetTitle("N Had");
  Hist.fSeedWedge_HadTowers->GetYaxis()->SetTitle("Seed Em Towers");
  new_folder->Add(Hist.fSeedWedge_HadTowers);
  
}

//==========================================================
void Histograms::GetCounterHistograms(TFolder* folder, CounterHists_t& Hist, 
						std::string sText, std::vector<std::string> vBinLabels)
{
	// this is a counter histogram. each bin is a counter of something
	// that you can define. then fill each bin as you want.
	
	if (!folder) {
		StdErr(__FILE__,__LINE__,3,"Pl. specify a valid folder. No histograms are added.");
		return;
	}

	if (vBinLabels.size() < 1) {
		StdErr(__FILE__,__LINE__,3,"Pl. specify bin labels in order to determine how many counters you want!.");
		return;
	}
  
	char name [200];
	sprintf(name,"Counter");
	sprintf(title,"%s - Event Counting Histogram", sText.c_str());
	Hist.iCounter = new TH1I(name,title,vBinLabels.size(),0,vBinLabels.size());

	for (int i=0; i < vBinLabels.size(); i++) {
		Hist.iCounter->GetXaxis()->SetBinLabel(i+1,vBinLabels[i].c_str());
	}
 	
	folder->Add(Hist.iCounter);

}


//==========================================================
void Histograms::GetPhoton1JetHistograms(TFolder* new_folder,
							Photon1JetHists_t& Hist, std::string sText)
{
	sTextToAdd = sText;
	if (!new_folder) {
		StdOut(__FILE__,__LINE__,3,"Null Histogram Folder. No histograms are added.");
		return;
	}
		
	sprintf(title, "%s - #gamma + Jet Invariant Mass",sTextToAdd.c_str());
	Hist.InvMass	 = new TH1F("InvMass",title,1000,0,1000);
	sprintf(title, "%s - #Delta#Phi ^{#gamma,Jet}",sTextToAdd.c_str());
	Hist.DelPhi		 = new TH1F("DelPhi",title,1400,-7,7);
	sprintf(title, "%s - #Delta#eta^{#gamma,Jet}",sTextToAdd.c_str());
	Hist.DelEta		 = new TH1F("DelEta",title,1000,-5,5);
	sprintf(title, "%s - Detector #Delta#eta^{#gamma,Jet}",sTextToAdd.c_str());
	Hist.DelEtaDet		 = new TH1F("DetDelEta",title,1000,-5,5);

	sprintf(title, "%s - Evt #eta^{#gamma} + Evt #eta{Jet}",sTextToAdd.c_str());
	Hist.DelEtaPlus		 = new TH1F("DelEtaPlus",title,1000,-5,5);
	sprintf(title, "%s - Det #eta^{#gamma} + Det #eta^{Jet}",sTextToAdd.c_str());
	Hist.DelEtaDetPlus		 = new TH1F("DelEtaDetPlus",title,1000,-5,5);

	
	sprintf(title, "%s - #DeltaR^{#gamma,Jet}=#sqrt{#Delta #Phi^{2}+#Delta #eta^{2}}",sTextToAdd.c_str());
	Hist.DelR		 = new TH1F("DelR",title,1000,0,10);
	sprintf(title, "%s - Jet to #gamma E_{T} ratio",sTextToAdd.c_str());
	Hist.EtRatio	 = new TH1F("EtRatio",title,1000,0,10);
		
	sprintf(title, "%s - Jet Emfr Vs #Delta#Phi^{#gamma,Jet}",sTextToAdd.c_str());
	Hist.JetEmfrVsDelPhi	= new TH2F("JetEmfrVsDelPhi",title,560,-7,7,1500,0,1.5);
	sprintf(title, "%s - Jet Emfr Vs #gamma Phi wedge",sTextToAdd.c_str());
	Hist.JetEmfrVsPhoPhiWedge = new TH2F("JetEmfrVsPhoPhiWedge",title,50,0,50,1500,0,1.5);

	
	new_folder->Add(Hist.InvMass);
	new_folder->Add(Hist.DelPhi);
	new_folder->Add(Hist.DelEta);
	new_folder->Add(Hist.DelEtaDet);
	new_folder->Add(Hist.DelEtaPlus);
	new_folder->Add(Hist.DelEtaDetPlus);
	new_folder->Add(Hist.DelR);
	new_folder->Add(Hist.EtRatio);
	new_folder->Add(Hist.JetEmfrVsDelPhi);
	new_folder->Add(Hist.JetEmfrVsPhoPhiWedge);
}


//==========================================================
void Histograms::GetPhoton2JetsHistograms(TFolder* new_folder,
							Photon2JetsHists_t& Hist, std::string sText)
{
	sTextToAdd = sText;
	if (!new_folder) {
		StdOut(__FILE__,__LINE__,3,"Null Histogram Folder. No histograms are added.");
		return;
	}
		
	sprintf(title, "%s - #gamma + 2 Jets Invariant Mass",sTextToAdd.c_str());
	Hist.InvMass	 = new TH1F("InvMass",title,1000,0,1000);
		
	new_folder->Add(Hist.InvMass);
}

//==========================================================
void Histograms::GetTwoJetsHistograms(TFolder* new_folder,
							TwoJetsHists_t& Hist, std::string sText)
{
	sTextToAdd = sText;
	if (!new_folder) {
		StdOut(__FILE__,__LINE__,3,"Null Histogram Folder. No histograms are added.");
		return;
	}
		
	sprintf(title, "%s - Two Jets Invariant Mass",sTextToAdd.c_str());
	Hist.InvMass	 = new TH1F("InvMass",title,1000,0,1000);
	sprintf(title, "%s - #Delta#Phi ^{Jet1,Jet2}",sTextToAdd.c_str());
	Hist.DelPhi		 = new TH1F("DelPhi",title,140,-7,7);
	sprintf(title, "%s - #Delta#eta",sTextToAdd.c_str());
	Hist.DelEta		 = new TH1F("DelEta",title,100,-5,5);
	sprintf(title, "%s - #DeltaR^{Jet1,Jet2}=#sqrt{#Delta #Phi^{2}+#Delta #eta^{2}}",sTextToAdd.c_str());
	Hist.DelR		 = new TH1F("DelR",title,1000,0,10);
	sprintf(title, "%s - Jet1 to Jet2 E_{T} ratio",sTextToAdd.c_str());
	Hist.EtRatio	 = new TH1F("EtRatio",title,1500,0,1.5);
		
	new_folder->Add(Hist.InvMass);
	new_folder->Add(Hist.DelPhi);
	new_folder->Add(Hist.DelEta);
	new_folder->Add(Hist.DelR);
	new_folder->Add(Hist.EtRatio);
}


//==========================================================
void Histograms::GetJetHistograms(TFolder* new_folder,
						JetHists_t& hist, std::string sText)
{
	if (!new_folder) {
		StdOut(__FILE__,__LINE__,3,"Null Histogram Folder. No histograms are added.");
		return;
	} else {
		sTextToAdd = sText;
		//sprintf(title, "%s - Detector Region",sTextToAdd.c_str());
		//hist.Detector 	= new TH1F("Detector",title,3,0,3);
		sprintf(title, "%s - Jet: Detector #eta",sTextToAdd.c_str());
		hist.DetEta 	= new TH1F("DetEta",title,2000,-10,10);
		sprintf(title, "%s - Jet: Detector #Phi",sTextToAdd.c_str());
		hist.DetPhi 	= new TH1F("DetPhi",title,140,-7,7);
		sprintf(title, "%s - Jet: E_{T}",sTextToAdd.c_str());
		hist.EtCorr 	= new TH1F("EtCorr",title,500,0,500);
		sprintf(title, "%s - Jet: Emfr",sTextToAdd.c_str());
		hist.Emfr 		= new TH1F("Emfr",title,400,0,4);
		sprintf(title, "%s - Jet: HadEm",sTextToAdd.c_str());
		hist.HadEm 		= new TH1F("HadEm",title,1000,0,10);
		sprintf(title, "%s - Jet: Number of Towers",sTextToAdd.c_str());
		hist.NTowers 	= new TH1F("NTowers",title,100,0,100);
		sprintf(title, "%s - Jet: Number of Tracks",sTextToAdd.c_str());
		hist.NTracks 	= new TH1F("NTracks",title,100,0,100);
		sprintf(title, "%s - Jet: Event #Phi",sTextToAdd.c_str());
		hist.EvtPhi 	= new TH1F("EvtPhi",title,140,-7,7);

		sprintf(title, "%s - Jet: SeedIEta",sTextToAdd.c_str());
		hist.SeedIEta 	= new TH1F("SeedIEta",title,50,0,50);
		sprintf(title, "%s - Jet: SeedIPhi",sTextToAdd.c_str());
		hist.SeedIPhi 	= new TH1F("SeedIPhi",title,46,-1,45);		//so we can see plug jets too


		sprintf(title, "%s -  #Delta #Phi (Jet, MET)",sTextToAdd.c_str());
		hist.JetMetDelPhi 	= new TH1F("JetMetDelPhi",title,140,-7,7);

		sprintf(ytitle,"Events/%.2f",hist.DetEta->GetBinWidth(1));
		hist.DetEta->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.DetPhi->GetBinWidth(1));
		hist.DetPhi->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.EtCorr->GetBinWidth(1));
		hist.EtCorr->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.Emfr->GetBinWidth(1));
		hist.Emfr->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.HadEm->GetBinWidth(1));
		hist.HadEm->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.NTowers->GetBinWidth(1));
		hist.NTowers->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.NTracks->GetBinWidth(1));
		hist.NTracks->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.EvtPhi->GetBinWidth(1));
		hist.EvtPhi->SetYTitle(ytitle);
		
		sprintf(ytitle,"Events/%.2f",hist.SeedIEta->GetBinWidth(1));
		hist.SeedIEta->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.SeedIPhi->GetBinWidth(1));
		hist.SeedIPhi->SetYTitle(ytitle);
		
		sprintf(ytitle,"Events/%.2f",hist.JetMetDelPhi->GetBinWidth(1));
		hist.JetMetDelPhi->SetYTitle(ytitle);
		
		new_folder->Add(hist.DetEta);
		new_folder->Add(hist.DetPhi);
		new_folder->Add(hist.EtCorr);
		new_folder->Add(hist.Emfr);
		new_folder->Add(hist.HadEm);
		new_folder->Add(hist.NTowers);
		new_folder->Add(hist.NTracks);
		new_folder->Add(hist.EvtPhi);
		new_folder->Add(hist.SeedIEta);
		new_folder->Add(hist.SeedIPhi);
		new_folder->Add(hist.JetMetDelPhi);
	}
	
}

//==========================================================
void Histograms::GetPhotonHistograms(TFolder* new_folder,
						PhotonHists_t& hist, std::string sText)
{
	if (!new_folder) {
		StdOut(__FILE__,__LINE__,3,"Null Histogram Folder. No histograms are added.");
	} else {
		sTextToAdd = sText;
		sprintf(title, "%s - Detector Region",sTextToAdd.c_str());
		hist.Detector 	= new TH1F("Detector",title,3,0,3);
		hist.Detector->GetXaxis()->SetBinLabel(1,"Central");
		hist.Detector->GetXaxis()->SetBinLabel(2,"Plug");
		hist.Detector->GetXaxis()->SetBinLabel(3,"Forward");
		sprintf(title, "%s - Detector #eta",sTextToAdd.c_str());
		hist.DetEta 	= new TH1F("DetEta",title,2000,-10,10);
		sprintf(title, "%s - Detector #Phi",sTextToAdd.c_str());
		hist.DetPhi 	= new TH1F("DetPhi",title,140,-7,7);
		sprintf(title, "%s - E_{T}",sTextToAdd.c_str());
		hist.EtCorr 	= new TH1F("EtCorr",title,500,0,500);
		sprintf(title, "%s - XCes",sTextToAdd.c_str());
		hist.XCes 		= new TH1F("XCes",title,640,-32,32);
		sprintf(title, "%s - ZCes",sTextToAdd.c_str());
		hist.ZCes 		= new TH1F("ZCes",title,2000,-250,250);
		sprintf(title, "%s - HadEm",sTextToAdd.c_str());
		hist.HadEm 		= new TH1F("HadEm",title,400,0,2.0);
		sprintf(title, "%s - #chi^{2} Mean",sTextToAdd.c_str());
		hist.Chi2Mean 	= new TH1F("Chi2Mean",title,1210,-1,120);
		sprintf(title, "%s - N3d",sTextToAdd.c_str());
		hist.N3d 		= new TH1F("N3d",title,50,0,50);
		sprintf(title, "%s - Iso",sTextToAdd.c_str());
		hist.Iso4 		= new TH1F("IsoEtCorr",title,3000,-5,25);
		sprintf(title, "%s - TrkPt",sTextToAdd.c_str());
		hist.TrkPt 		= new TH1F("Trkpt",title,1000,0,100);
		sprintf(title, "%s - TrkIso",sTextToAdd.c_str());
		hist.TrkIso 	= new TH1F("TrkIso",title,150,0,15);
		sprintf(title, "%s - CES(2nd) Wire",sTextToAdd.c_str());
		hist.Ces2Wire 	= new TH1F("Ces2Wire",title,400,-40,40);
		sprintf(title, "%s - CES(2nd) Strip",sTextToAdd.c_str());
		hist.Ces2Strip = new TH1F("Ces2Strip",title,400,-40,40);
		sprintf(title, "%s - EM Time ",sTextToAdd.c_str());
		hist.EmTime 	= new TH1F("EmTime",title,300,-100,200);
		sprintf(title, "%s - E/p",sTextToAdd.c_str());
		hist.EoverP 	= new TH1F("EoverP",title,110,-1,10);
		sprintf(title, "%s - Track Z0",sTextToAdd.c_str());
		hist.TrkZ0 	= new TH1F("TrkZ0",title,500,-500,500);

		sprintf(title, "%s - EM Timing of #gamma",sTextToAdd.c_str());
		hist.EmTimeVsRun 	 = new TProfile("EmTimeVsRun",title,270000,130000,400000,-200,200);
		sprintf(title, "%s - E_{T}^{corr} of #gamma",sTextToAdd.c_str());
		hist.EtCorrVsRun 	 = new TProfile("EtCorrVsRun",title,270000,130000,400000,-200,200);

		sprintf(title,"%s - #Phi wedge", sTextToAdd.c_str());
		hist.PhiWedge = new TH1F("PhiWedge",title,24,-0.5,23.5);
  
		sprintf(title, "%s -  #Delta #Phi (MET, #gamma)",sTextToAdd.c_str());
		hist.PhoMetDelPhi 	= new TH1F("PhoMetDelPhi",title,140,-7,7);

		sprintf(title, "%s -  #gamma IDs : (=1 if the set of cuts is passed.)",sTextToAdd.c_str());
		hist.PhoIDs 	= new TH1F("PhoIDs",title,10,0,10);
		hist.PhoIDs->GetXaxis()->SetBinLabel(1,"Tight #gamma");
		hist.PhoIDs->GetXaxis()->SetBinLabel(2,"Loose #gamma");
		hist.PhoIDs->GetXaxis()->SetBinLabel(3,"Sideband #gamma");
		hist.PhoIDs->GetXaxis()->SetBinLabel(4,"#gamma-like e-ID");
		hist.PhoIDs->GetXaxis()->SetBinLabel(5,"STD Tight e-ID");
		hist.PhoIDs->GetXaxis()->SetBinLabel(6,"STD Loose e-ID");
		hist.PhoIDs->GetXaxis()->SetBinLabel(7,"PHX ID ");
		hist.PhoIDs->GetXaxis()->SetBinLabel(8,"Coversion ID ");
		hist.PhoIDs->GetXaxis()->SetBinLabel(9,"BH ID");
		

  
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
		sprintf(ytitle,"Events/%.2f",hist.EmTime->GetBinWidth(1));
		hist.EmTime->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.EoverP->GetBinWidth(1));
		hist.EoverP->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.TrkZ0->GetBinWidth(1));
		hist.TrkZ0->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.PhiWedge->GetBinWidth(1));
		hist.PhiWedge->SetYTitle(ytitle);
		sprintf(ytitle,"Events/%.2f",hist.PhoMetDelPhi->GetBinWidth(1));
		hist.PhoMetDelPhi->SetYTitle(ytitle);
		
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
		new_folder->Add(hist.EmTime);
		new_folder->Add(hist.EoverP);
		new_folder->Add(hist.TrkZ0);
		new_folder->Add(hist.EmTimeVsRun);
		new_folder->Add(hist.EtCorrVsRun);
		new_folder->Add(hist.PhiWedge);
		new_folder->Add(hist.PhoMetDelPhi);
		new_folder->Add(hist.PhoIDs);
	}
	
}



//==========================================================
void Histograms::Print(Option_t *opt)
{
	printf("============= Histograms:: Values of Current Data Members =============\n");
}
