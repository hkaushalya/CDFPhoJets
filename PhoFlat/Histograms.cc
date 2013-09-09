///////////////////////////////////////////////////////////
// Defines all the common histograms. Any moddule can    //
// use them easily and quickly.                          //
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>
/*{{{*/
/********************************************************************
 *
 * $Id: Histograms.cc,v 1.15 2011/05/26 18:45:02 samantha Exp $
 * $Log: Histograms.cc,v $
 * Revision 1.15  2011/05/26 18:45:02  samantha
 * ADDED: PtBal and PtBalVsMass hists to Pho2Jet hists.
 *
 * Revision 1.14  2011/04/26 23:45:12  samantha
 * ADDED: EmTime hist for Jet Histograms and DelEmTime for time separation
 * between photon and jet to photon-jet histograms.
 *
 * Revision 1.13  2010/08/25 15:30:01  samantha
 * ADDED: MET study hists for photon, jet and others. Profile histograms to see the
 * region of MET. Also Pt plots for parent objects.
 *
 * Revision 1.12  2010/06/21 19:47:35  samantha
 * MODIFIED: MetX,MetY, InvMass, DetEta, EvtEta, HadEm and few other variables histograms
 * bin sizes.
 * ADDED: hists for phi separation of photon and leading jet with the met before
 * and after the delphi cut of 0.4
 *
 * Revision 1.11  2010/04/13 16:12:09  samantha
 * ADDED: MetX, MetY hists to events histograms.
 *
 * Revision 1.10  2010/04/07 17:01:18  samantha
 * ADDED:1. delPhi(jet,MET) to the Jet hists collection.
 *       2. delPhi(pho,MET) to the Photon hists collection.
 *
 * Revision 1.9  2010/03/24 18:40:58  samantha
 * ADDED: VertexZ plot for highest sumPt vertex z distribution to the list of Event
 * hists.
 * MODIFIED: Increased the number of bins in HADEM hists.
 *
 * Revision 1.8  2010/02/09 23:56:11  samantha
 * ADDED a EventWeight Profile hist to the list of photon hist to see the weights used for sideband reweighing.
 *
 * Revision 1.7  2010/01/07 22:40:40  samantha
 * ADDED: 1. pho-jet Pt and E balance hists to Photon1JetHists
 * MODIFIED: 1. HadEm hist in PhotonHists, is modified for smaller bin sizes.
 *
 * Revision 1.6  2009/12/07 21:16:47  samantha
 * ADDED:  1. CprWeight histogram to the photon hists.
 * 	2. CVS auto log into into the file.
 *
 *
 *******************************************************************/
/*}}}*/

#include <fstream>
#include <iostream>
#include "PhoFlat/Histograms.hh"
#include "PhoFlat/FreeFunctions.hh"
#include <sstream>


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
void Histograms::GetEventHistograms(EventHists_t& Hist, std::string sText)
{
	// Histograms of Event variables like MEt, SumEt etc.

	sTextToAdd = sText;
		sprintf(title, "%s - ME_{T}",sTextToAdd.c_str());
		Hist.Met 			 = new TH1F("Met",title,300,0,1500);
		
		sprintf(title, "%s - ME_{T}-X",sTextToAdd.c_str());
		Hist.MetX 			 = new TH1F("MetX",title,500,-500,500);
		
		sprintf(title, "%s - ME_{T}-Y",sTextToAdd.c_str());
		Hist.MetY 			 = new TH1F("MetY",title,500,-500,500);
		
		sprintf(title, "%s - Generated ME_{T}",sTextToAdd.c_str());
		Hist.Met_Gen_d 	 = new TH1F("Met_Gen_d",title,300,0,1500);
		
		sprintf(title, "%s - SumE_{T}",sTextToAdd.c_str());
		Hist.Sumet 		 = new TH1F("Sumet",title,300,0,1500);
		
		float Ht_min=0, Ht_max=1500;
		int Ht_bins=300;
		
		sprintf(title, "%s - H_{T}",sTextToAdd.c_str());
		Hist.Ht 		 = new TH1F("Ht",title,Ht_bins,Ht_min,Ht_max);


		sprintf(title, "%s - H_{T} (CPR Weighted)",sTextToAdd.c_str());
		Hist.HtWgt 		 = new TH1F("Ht_wgt",title,Ht_bins,Ht_min,Ht_max);
		for (int i=0; i < 8; i++)
		{
			std::ostringstream tt, nn;
			tt << sTextToAdd << " - H_{T} Systemactic_" << i;
			nn << "Ht_sys_" << i;
			Hist.HtSys[i] = new TH1F(nn.str().c_str(), tt.str().c_str(), Ht_bins, Ht_min, Ht_max);
		}
		
		sprintf(title, "%s - # of Vertices",sTextToAdd.c_str());
		Hist.NVertices 	 = new TH1F("NVertices",title,25,0,25);
		
		sprintf(title, "%s - # of Class 12 Vertices",sTextToAdd.c_str());
		Hist.N12Vertices = new TH1F("N12Vertices",title,25,0,25);

		sprintf(title, "%s - z position of the highest sump Pt vertex", sTextToAdd.c_str());
		Hist.VertexZ = new TH1F("VertexZ",title,400,-200,200);
		
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
		Hist.NJet15VsMet      = new TH2F("NJet15VsMet",title,300,0,1500,50,0,50);
		sprintf(title, "%s - NJets(15GeV) Vs SumE_{T}",sTextToAdd.c_str());
		Hist.NJet15VsSumet 	 = new TH2F("NJet15VsSumet",title,300,0,1500,50,0,50);
		sprintf(title, "%s - # of Class 12 Vertices Vs ME_{T}",sTextToAdd.c_str());
		Hist.N12VerticesVsMet 	 = new TH2F("N12VerticesVsMet",title,300,0,1500,50,0,50);
		sprintf(title, "%s - # of Class 12 Vertices Vs SumE_{T}",sTextToAdd.c_str());
		Hist.N12VerticesVsSumet = new TH2F("N12VerticesVsSumet",title,300,0,1500,50,0,50);


		sprintf(title, "%s - Average # of Vertices Vs #slash{E}_{T}",sTextToAdd.c_str());
		Hist.N12VerticesMetProf = new TProfile("N12VerticesMetProf",title,100,0,300,0,10);

		sprintf(title, "%s - Average # of Vertices Vs #slash{E}_{T}",sTextToAdd.c_str());
		Hist.N12VerticesSumetProf = new TProfile("N12VerticesSumetProf",title,100,0,500,0,10);
		
		sprintf(title, "%s - Triggers Passed",sTextToAdd.c_str());
		Hist.Triggers			 = new TH1F("Triggers",title,10,0,10);
		Hist.Triggers->GetXaxis()->SetBinLabel(1,"phoIso25");
		Hist.Triggers->GetXaxis()->SetBinLabel(2,"pho50");
		Hist.Triggers->GetXaxis()->SetBinLabel(3,"pho70");
}

//==========================================================
void Histograms::GetHaloHistograms(BeamHaloHists_t& Hist, 
						std::string sText)
{
  char name [200];
  char title[200];

  sprintf(name,"HaloSeedWedgeEM_SeedWedgeHAD");
  sprintf(title,"%s - N_{EMtwr} vs. N_{HADtwr} in seed wedge", sText.c_str());
  Hist.fHaloSeedWedgeEM_SeedWedgeHAD = new TH2F(name,title,25,-0.5,24.5,25,-0.5,24.5);
  

  sprintf(name,"HaloSeedWedge");
  sprintf(title,"%s - number of EM towers in seed wedge, Max's version", sText.c_str());
  Hist.fHaloSeedWedge=new TH1F(name,title,25,-0.5,24.5);
  
  sprintf(name,"HaloSideWedge");
  sprintf(title,"%s - number of EM towers in side wedges, Max's version", sText.c_str());
  Hist.fHaloSideWedge=new TH1F(name,title,45,-0.5,44.5);
  
  sprintf(name,"HaloSeedWedgeH");
  sprintf(title,"%s - number of HAD towers in seed wedge, my version", sText.c_str());
  Hist.fHaloSeedWedgeH=new TH1F(name,title,25,-0.5,24.5);
  
  sprintf(name,"HaloSideWedgeH");
  sprintf(title,"%s - number of HAD towers in side wedges, my version", sText.c_str());
  Hist.fHaloSideWedgeH=new TH1F(name,title,45,-0.5,44.5);
  
  sprintf(name,"HaloEastNHad");
  sprintf(title,"%s - number of HAD towers on East side", sText.c_str());
  Hist.fHaloEastNHad=new TH1F(name,title,35,-0.5,34.5);
  
  sprintf(name,"HaloWestNHad");
  sprintf(title,"%s - number of HAD towers on West side", sText.c_str());
  Hist.fHaloWestNHad=new TH1F(name,title,35,-0.5,34.5);
  
  sprintf(name,"HaloNHad");
  sprintf(title,"%s - number of HAD towers on East+West", sText.c_str());
  Hist.fHaloNHad=new TH1F(name,title,70,-0.5,69.5);
  
  sprintf(name,"HaloEmSeedE");
  sprintf(title,"%s - EM energy of seed wedge towers", sText.c_str());
  Hist.fHaloEmSeedE=new TH1F(name,title,1000,0.0,1000.0);
  
  sprintf(name,"HaloEmSideE");
  sprintf(title,"%s - EM energy of side wedge towers", sText.c_str());
  Hist.fHaloEmSideE=new TH1F(name,title,1000,0.0,1000.0);
  
  sprintf(name,"HaloHadSeedE");
  sprintf(title,"%s - Plug HAD energy of seed wedge towers", sText.c_str());
  Hist.fHaloHadSeedE=new TH1F(name,title,1000,0.0,1000.0);
  
  sprintf(name,"HaloSeedWedgeHadE");
  sprintf(title,"%s - HAD energy of seed wedge towers", sText.c_str());
  Hist.fHaloSeedWedgeHadE=new TH1F(name,title,5000,0.0,100.0);
  
  sprintf(name,"HaloEast");
  sprintf(title,"%s - HaloEast", sText.c_str());
  Hist.fHaloEast=new TH1F(name,title,10,-0.5,9.5);
  
  sprintf(name,"HaloWest");
  sprintf(title,"%s - HaloWest", sText.c_str());
  Hist.fHaloWest=new TH1F(name,title,10,-0.5,9.5);
  
  sprintf(name,"HaloEastWest");
  sprintf(title,"%s - HaloWest+HaloEast", sText.c_str());
  Hist.fHaloEastWest=new TH1F(name,title,10,-0.5,9.5);
  
  sprintf(name,"HaloSeedWedgeR");
  sprintf(title,"%s - number of towers in seed wedge, Ray's version", sText.c_str());
  Hist.fHaloSeedWedgeR=new TH1F(name,title,25,-0.5,24.5);
  
  sprintf(name,"HaloSideWedgeR");
  sprintf(title,"%s - number of towers in side wedges, Ray's version", sText.c_str());
  Hist.fHaloSideWedgeR=new TH1F(name,title,45,-0.5,44.5);
  
  sprintf(name,"HaloTwrInRawR");
  sprintf(title,"%s - number of continues towers in seed wedge, Ray's version", sText.c_str());
  Hist.fHaloTwrInRawR=new TH1F(name,title,25,-0.5,24.5);
  

  sprintf(name,"HaloCesStripEoverE");
  sprintf(title,"%s - strip energy E_{CES}/E_{#gamma}", sText.c_str());
  Hist.fHaloCesStripEoverE=new TH1F(name,title,500,0.0,5.0);
  
  sprintf(name,"HaloCesWireEoverE");
  sprintf(title,"%s - wire energy E_{CES}/E_{#gamma}", sText.c_str());
  Hist.fHaloCesWireEoverE=new TH1F(name,title,500,0.0,5.0);
  
  sprintf(name,"HaloHadTDC");
  sprintf(title,"%s - photons Had TDC", sText.c_str());
  Hist.fHaloHadTDC=new TH1F(name,title,800,-100.0,100.0);
  

  sprintf(name,"HaloPhiWedge");
  sprintf(title,"%s - wedge number of halo candidate", sText.c_str());
  Hist.fHaloPhiWedge=new TH1F(name,title,24,-0.5,23.5);
  
  sprintf(name,"HaloEta");
  sprintf(title,"%s - #eta of halo candidate", sText.c_str());
  Hist.fHaloEta=new TH1F(name,title,120,-3.0,3.0);
  
  sprintf(name,"HaloIso");
  sprintf(title,"%s -  CalIso of halo candidate", sText.c_str());
  Hist.fHaloIso= new TH1F(name,title,1400,-4.0,10.0);
  
  sprintf(name,"HaloHadEm");
  sprintf(title,"%s - Had/Em of halo candidate", sText.c_str());
  Hist.fHaloHadEm= new TH1F(name,title,400,-0.5,1.5);
  
  sprintf(name,"HaloCesChi2");
  sprintf(title,"%s - CES(strip+wire) #chi^{2} of halo candidate", sText.c_str());
  Hist.fHaloCesChi2= new TH1F(name,title,500,0.0,100.0);
  
  sprintf(name,"HaloEt");
  sprintf(title,"%s -  corrected E_{T} of halo candidate", sText.c_str());
  Hist.fHaloEt= new TH1F(name,title,1000,0.0,1000.0);
  


	// my stuff

  sprintf(name,"HaloEvt_SeedWedge_Nhad");
  sprintf(title,"%s - Em towers in seed wedge and N Had towers (eat+west)", sText.c_str());
  Hist.fSeedWedge_HadTowers = new TH2F(name,title,20,0,20,20,0,20);
  Hist.fSeedWedge_HadTowers->GetXaxis()->SetTitle("N Had");
  Hist.fSeedWedge_HadTowers->GetYaxis()->SetTitle("Seed Em Towers");
  
  
}

//==========================================================
void Histograms::GetCounterHistograms(CounterHists_t& Hist, 
						std::string sText, std::vector<std::string> vBinLabels)
{
	// this is a counter histogram. each bin is a counter of something
	// that you can define. then fill each bin as you want.
	
	if (vBinLabels.size() < 1) {
		StdErr(__FILE__,__LINE__,3,"Pl. specify bin labels in order to determine how many counters you want!.");
		return;
	}
  
	char name [200];
	sprintf(name,"Counter");
	sprintf(title,"%s - Event Counting Histogram", sText.c_str());
	Hist.iCounter = new TH1I(name,title,vBinLabels.size(),0,vBinLabels.size());

	for (int i=0; i < (int) vBinLabels.size(); i++) {
		Hist.iCounter->GetXaxis()->SetBinLabel(i+1,vBinLabels[i].c_str());
	}

}


//==========================================================
void Histograms::GetPhoton1JetHistograms(	Photon1JetHists_t& Hist,
								std::string sText)
{
	sTextToAdd = sText;
		
	sprintf(title, "%s - #gamma + Jet Invariant Mass",sTextToAdd.c_str());
	Hist.InvMass	 = new TH1F("InvMass",title,1000,0,1000);
	sprintf(title, "%s - #Delta#phi ^{#gamma,Jet}",sTextToAdd.c_str());
	Hist.DelPhi		 = new TH1F("DelPhi",title,140,-7,7);
	sprintf(title, "%s - #Delta#eta^{#gamma,Jet}",sTextToAdd.c_str());
	Hist.DelEta		 = new TH1F("DelEta",title,100,-5,5);
	sprintf(title, "%s - #DeltaR^{#gamma,Jet}=#sqrt{#Delta #phi^{2}+#Delta #eta^{2}}",sTextToAdd.c_str());
	Hist.DelR		 = new TH1F("DelR",title,100,0,10);
	sprintf(title, "%s - Jet to #gamma E_{T} ratio",sTextToAdd.c_str());
	Hist.EtRatio	 = new TH1F("EtRatio",title,50,0,10);
		
	sprintf(title, "%s - Jet Emfr Vs #Delta#phi^{#gamma,Jet}",sTextToAdd.c_str());
	Hist.JetEmfrVsDelPhi	= new TH2F("JetEmfrVsDelPhi",title,140,-7,7,15,0,1.5);
	sprintf(title, "%s - Jet Emfr Vs #gamma phi wedge",sTextToAdd.c_str());
	Hist.JetEmfrVsPhoPhiWedge = new TH2F("JetEmfrVsPhoPhiWedge",title,50,0,50,15,0,1.5);

	sprintf(title, "%s - Photon, Jet Balance;#frac{E_{T}^{jet}}{E_{T}^{#gamma}} - 1;Events", sTextToAdd.c_str());
	Hist.PhoJetPtBalance = new TH1F("PhoJetPtBalance",title,150,-1.,4);
	
	sprintf(title, "%s - Photon, Jet Balance;#frac{E^{jet}}{E^{#gamma}} - 1;Events", sTextToAdd.c_str());
	Hist.PhoJetEBalance = new TH1F("PhoJetEBalance",title,150,-1.,4);

	sprintf(title, "%s - Photon+Jet: P_{T}",sTextToAdd.c_str());
	Hist.Pt 	= new TH1F("Pt",title,250,0,500);

	sprintf(title, "%s - EM Time Difference",sTextToAdd.c_str());
	Hist.DelEmTime 	= new TH1F("DeltaEmTime",title,300,-100,200);
}


//==========================================================
void Histograms::GetPhoton2JetsHistograms(Photon2JetsHists_t& Hist,
								std::string sText)
{
	sTextToAdd = sText;
		
	sprintf(title, "%s - #gamma + 2 Jets Invariant Mass",sTextToAdd.c_str());
	Hist.InvMass	 = new TH1F("InvMass",title,1000,0,1000);
	sprintf(title, "%s - #gamma+Jet1+Jet2: P_{T}",sTextToAdd.c_str());
	Hist.Pt 	= new TH1F("Pt",title,250,0,500);
	sprintf(title, "%s - #gamma+Jet1+Jet2: P_{T}^{j1j2}/P_{T}^{#gamma}",sTextToAdd.c_str());
	Hist.PtBal  = new TH1F("PtBal",title,100,0,10);
	sprintf(title, "%s - #gamma+Jet1+Jet2: P_{T}^{#gamma}/P_{T}^{j1j2} Vs M{j1,j2}",sTextToAdd.c_str());
	Hist.PtBalVsMass  = new TProfile("PtBalVsMass",title,500,0,500,0,5);
		
}

//==========================================================
void Histograms::GetTwoJetsHistograms(TwoJetsHists_t& Hist,
										std::string sText)
{
	sTextToAdd = sText;
		
	sprintf(title, "%s - Two Jets Invariant Mass",sTextToAdd.c_str());
	Hist.InvMass	 = new TH1F("InvMass",title,1000,0,1000);
	sprintf(title, "%s - #Delta#phi ^{Jet1,Jet2}",sTextToAdd.c_str());
	Hist.DelPhi		 = new TH1F("DelPhi",title,140,-7,7);
	sprintf(title, "%s - #Delta#eta",sTextToAdd.c_str());
	Hist.DelEta		 = new TH1F("DelEta",title,100,-5,5);
	sprintf(title, "%s - #DeltaR^{Jet1,Jet2}=#sqrt{#Delta #phi^{2}+#Delta #eta^{2}}",sTextToAdd.c_str());
	Hist.DelR		 = new TH1F("DelR",title,100,0,10);
	sprintf(title, "%s - Jet1 to Jet2 E_{T} ratio",sTextToAdd.c_str());
	Hist.EtRatio	 = new TH1F("EtRatio",title,15,0,1.5);
	sprintf(title, "%s - Jet1+Jet2: P_{T}",sTextToAdd.c_str());
	Hist.Pt 	= new TH1F("Pt",title,250,0,500);
}


//==========================================================
void Histograms::GetJetHistograms(JetHists_t& hist, std::string sText)
{
		sTextToAdd = sText;
		//sprintf(title, "%s - Detector Region",sTextToAdd.c_str());
		//hist.Detector 	= new TH1F("Detector",title,3,0,3);
		sprintf(title, "%s - Jet: Detector #eta",sTextToAdd.c_str());
		hist.DetEta 	= new TH1F("DetEta",title,200,-5,5);
		sprintf(title, "%s - Jet: Event #eta",sTextToAdd.c_str());
		hist.EvtEta 	= new TH1F("EvtEta",title,200,-5,5);
		sprintf(title, "%s - Jet: Detector #phi",sTextToAdd.c_str());
		hist.DetPhi 	= new TH1F("DetPhi",title,140,-7,7);
		sprintf(title, "%s - Jet: E_{T}",sTextToAdd.c_str());
		hist.EtCorr 	= new TH1F("EtCorr",title,250,0,500);
		sprintf(title, "%s - Jet: Emfr",sTextToAdd.c_str());
		hist.Emfr 		= new TH1F("Emfr",title,40,0,4);
		sprintf(title, "%s - Jet: HadEm",sTextToAdd.c_str());
		hist.HadEm 		= new TH1F("HadEm",title,200,0,2);
		sprintf(title, "%s - Jet: Number of Towers",sTextToAdd.c_str());
		hist.NTowers 	= new TH1F("NTowers",title,50,0,50);
		sprintf(title, "%s - Jet: Number of Tracks",sTextToAdd.c_str());
		hist.NTracks 	= new TH1F("NTracks",title,50,0,50);
		sprintf(title, "%s - Jet: Event #phi",sTextToAdd.c_str());
		hist.EvtPhi 	= new TH1F("EvtPhi",title,140,-7,7);

		sprintf(title, "%s - Jet: SeedIEta",sTextToAdd.c_str());
		hist.SeedIEta 	= new TH1F("SeedIEta",title,50,0,50);
		sprintf(title, "%s - Jet: SeedIPhi",sTextToAdd.c_str());
		hist.SeedIPhi 	= new TH1F("SeedIPhi",title,46,-1,45);		//so we can see plug jets too

		sprintf(title, "%s - Jet: #delta#phi(jet,#slash{E}_{T})",sTextToAdd.c_str());
		hist.JetMetDelPhi	= new TH1F("JetMetDelPhi",title, 700, -3.5, 3.5);
		sprintf(ytitle,"Events/%.2f",hist.JetMetDelPhi->GetBinWidth(1));
		hist.JetMetDelPhi->SetYTitle(ytitle);

		sprintf(title, "%s - Closest Jet before cut: #delta#phi(jet,#slash{E}_{T})",sTextToAdd.c_str());
		hist.ClosestJetMetDelPhi_b4	= new TH1F("ClosesetJetMetDelPhi_b4",title, 700, -3.5, 3.5);
		sprintf(ytitle,"Events/%.2f",hist.ClosestJetMetDelPhi_b4->GetBinWidth(1));
		hist.ClosestJetMetDelPhi_b4->SetYTitle(ytitle);
		
		sprintf(title, "%s - Closest Jet after cut: #delta#phi(jet,#slash{E}_{T})",sTextToAdd.c_str());
		hist.ClosestJetMetDelPhi_a4	= new TH1F("ClosesetJetMetDelPhi_a4",title, 700, -3.5, 3.5);
		sprintf(ytitle,"Events/%.2f",hist.ClosestJetMetDelPhi_a4->GetBinWidth(1));
		hist.ClosestJetMetDelPhi_a4->SetYTitle(ytitle);

		sprintf(title, "%s - Average Jet E_{T} Vs #slash{E}_{T}",sTextToAdd.c_str());
		hist.JetEtMetProf      = new TProfile("JetEtVsMet",title,350,0,350,0,600);
		sprintf(title, "%s - Average Jet E Vs #slash{E}_{T}",sTextToAdd.c_str());
		hist.JetEMetProf      = new TProfile("JetEVsMet",title,350,0,350,0,600);
		sprintf(title, "%s - Average Jet #eta Vs #slash{E}_{T}",sTextToAdd.c_str());
		hist.JetEtaMetProf      = new TProfile("JetEtaVsMet",title,350,0,350,-4,4);
		
		sprintf(title, "%s - EM Time ",sTextToAdd.c_str());
		hist.EmTime 	= new TH1F("EmTime",title,300,-100,200);
		
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
		
}

//==========================================================
void Histograms::GetPhotonHistograms(PhotonHists_t& hist, std::string sText)
{
		sTextToAdd = sText;
		sprintf(title, "%s - Detector Region",sTextToAdd.c_str());
		hist.Detector 	= new TH1F("Detector",title,3,0,3);
		hist.Detector->GetXaxis()->SetBinLabel(1,"Central");
		hist.Detector->GetXaxis()->SetBinLabel(2,"Plug");
		hist.Detector->GetXaxis()->SetBinLabel(3,"Forward");
		sprintf(title, "%s - Detector #eta",sTextToAdd.c_str());
		hist.DetEta 	= new TH1F("DetEta",title,200,-5,5);
		sprintf(title, "%s - Detector #phi",sTextToAdd.c_str());
		hist.DetPhi 	= new TH1F("DetPhi",title,140,-7,7);
		sprintf(title, "%s - E_{T}",sTextToAdd.c_str());
		hist.EtCorr 	= new TH1F("EtCorr",title,250,0,500);
		sprintf(title, "%s - XCes",sTextToAdd.c_str());
		hist.XCes 		= new TH1F("XCes",title,320,-32,32);
		sprintf(title, "%s - ZCes",sTextToAdd.c_str());
		hist.ZCes 		= new TH1F("ZCes",title,250,-250,250);
		sprintf(title, "%s - HadEm",sTextToAdd.c_str());
		hist.HadEm 		= new TH1F("HadEm",title,400,-0.5,1.5);
		sprintf(title, "%s - #chi^{2} Mean",sTextToAdd.c_str());
		hist.Chi2Mean 	= new TH1F("Chi2Mean",title,121,-1,120);
		sprintf(title, "%s - N3d",sTextToAdd.c_str());
		hist.N3d 		= new TH1F("N3d",title,30,0,30);
		sprintf(title, "%s - Iso",sTextToAdd.c_str());
		hist.Iso4 		= new TH1F("IsoEtCorr",title,60,-5,25);
		sprintf(title, "%s - TrkPt",sTextToAdd.c_str());
		hist.TrkPt 		= new TH1F("Trkpt",title,200,0,50);
		sprintf(title, "%s - TrkIso",sTextToAdd.c_str());
		hist.TrkIso 	= new TH1F("TrkIso",title,200,0,50);
		sprintf(title, "%s - CES(2nd) Wire",sTextToAdd.c_str());
		hist.Ces2Wire 	= new TH1F("Ces2Wire",title,100,0,20);
		sprintf(title, "%s - CES(2nd) Strip",sTextToAdd.c_str());
		hist.Ces2Strip = new TH1F("Ces2Strip",title,100,0,20);
		sprintf(title, "%s - EM Time ",sTextToAdd.c_str());
		hist.EmTime 	= new TH1F("EmTime",title,300,-100,200);
		sprintf(title, "%s - E/p",sTextToAdd.c_str());
		hist.EoverP 	= new TH1F("EoverP",title,200,0,10);
		sprintf(title, "%s - Track Z0",sTextToAdd.c_str());
		hist.TrkZ0 	= new TH1F("TrkZ0",title,500,-500,500);

		sprintf(title, "%s - EM Timing of #gamma",sTextToAdd.c_str());
		hist.EmTimeVsRun 	 = new TProfile("EmTimeVsRun",title,270000,130000,400000,-200,200);
		sprintf(title, "%s - E_{T}^{corr} of #gamma",sTextToAdd.c_str());
		hist.EtCorrVsRun 	 = new TProfile("EtCorrVsRun",title,270000,130000,400000,-200,200);

		sprintf(title,"%s - #Phi wedge", sTextToAdd.c_str());
		hist.PhiWedge = new TH1F("PhiWedge",title,24,-0.5,23.5);
		sprintf(title,"%s - CES/CPR Weight", sTextToAdd.c_str());
		hist.CprWeight = new TH1F("CesCprWeight",title,200,-5,15);
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


		sprintf(title, "%s - Event Weights;E_{T}^{#gamma};Weight",sTextToAdd.c_str());
		hist.EventWeights = new TProfile("EventWeights",title,300,0,300);

	sprintf(title, "%s - Photon: #delta#phi(#gamma,#slash{E}_{T})",sTextToAdd.c_str());
	hist.PhoMetDelPhi	= new TH1F("PhoMetDelPhi",title, 700, -3.5, 3.5);
	sprintf(ytitle,"Events/%.2f",hist.PhoMetDelPhi->GetBinWidth(1));
	hist.PhoMetDelPhi->SetYTitle(ytitle);

		
		sprintf(title, "%s - before cuts: #delta#phi(#gamma, #slash{E}_{T})",sTextToAdd.c_str());
		hist.ClosestPhoMetDelPhi_b4	= new TH1F("ClosesetPhoMetDelPhi_b4",title, 700, -3.5, 3.5);
		sprintf(ytitle,"Events/%.2f",hist.ClosestPhoMetDelPhi_b4->GetBinWidth(1));
		hist.ClosestPhoMetDelPhi_b4->SetYTitle(ytitle);
		
		sprintf(title, "%s - after cut: #delta#phi(#gamma, #slash{E}_{T})",sTextToAdd.c_str());
		hist.ClosestPhoMetDelPhi_a4	= new TH1F("ClosesetPhoMetDelPhi_a4",title, 700, -3.5, 3.5);
		sprintf(ytitle,"Events/%.2f",hist.ClosestPhoMetDelPhi_a4->GetBinWidth(1));
		hist.ClosestPhoMetDelPhi_a4->SetYTitle(ytitle);

		sprintf(title, "%s - Average #gamma E_{T} Vs #slash{E}_{T}",sTextToAdd.c_str());
		hist.PhoEtMetProf      = new TProfile("PhoEtVsMet",title,350,0,350,0,400);
		sprintf(title, "%s - Average #gamma E Vs #slash{E}_{T}",sTextToAdd.c_str());
		hist.PhoEMetProf      = new TProfile("PhoEVsMet",title,350,0,350,0,400);
		
		sprintf(title, "%s - Average #gamma #eta Vs #slash{E}_{T}",sTextToAdd.c_str());
		hist.PhoEtaMetProf      = new TProfile("PhoEtaVsMet",title,350,0,350,-4,4);
	
}
