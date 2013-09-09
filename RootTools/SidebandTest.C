/*This is to test the photon sideband.
 * If it correclty represent the fakes passing
 * tight photon id cuts. I am using a di-jet
 * MC sample from Sasha.
 * See Elog#https://hep.baylor.edu/hep/samantha/1343
 *
 * To Run:
 * root>.L SidebandTest.C+
 * root>SidebandTest()
 *
 *Sep 25, 2009
 */

/* $Id: SidebandTest.C,v 1.4 2010/02/13 21:44:56 samantha Exp $
 * $Log: SidebandTest.C,v $
 * Revision 1.4  2010/02/13 21:44:56  samantha
 * MODIFIED: To get the correct values for x- values, the average Pt (As the bin
 * sizes are different). So the final hist is a TGraphError.
 *
 * Revision 1.3  2009/10/07 00:05:51  samantha
 * This version in designed to get max stats by using variable bin sizes and
 * picking the MC samples with more stats in a given region. Dropped samples
 * that had very low stats (eg 60GeV sample), and used the neighbours's stat
 * to fill that range.See elog#https://hep.baylor.edu/hep/samantha/1377
 *
 * Revision 1.2  2009/10/06 02:40:53  samantha
 * This has more accurate way of normalizing successive samples with different
 * thresholds. In the previous version the final result changed its trends with
 * the bin size of the histogram. This is because of the poor statistics in the
 * overlaping regions of adjacent samples where I calculate the normalization
 * constant.Even if I used a very finely binned histogram to derive the
 * normalization and then rebinned the results at the very end, it still showed
 * bin size dependency. So the solution was not to derive normalization from a
 * histogram at all. Instead count the needed numbers using simple counters and
 * derive scales from them. This scaling was little trikier than I initially
 * anticipated. This version still makes a histogram with constant bin size. Next
 * version will have variable bin sizes.
 * See Elog#https://hep.baylor.edu/hep/samantha/1371
 *
 * Revision 1.1  2009/09/29 21:40:55  samantha
 * This is to test the sideband using di-jet MC samples. We a testing how good is
 * the sideband sample in representing the fakes in data  (or jets passsing tight
 * photon ID cuts.) See Elog# https://hep.baylor.edu/hep/samantha/1343 .
 *
 *
 */
#include <iostream>
#include "TTree.h" 
#include "TCanvas.h"
#include "TChain.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLegend.h"
#include <sstream>
#include "TMath.h"
#include <iomanip>
#include <cmath>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"

struct PhoIdVars_
{
	Double_t E;
	Double_t Et;
	Double_t Eta;
	Double_t EtaDet;
	Double_t XCes;
	//Double_t ZCes;
	Double_t HadEm;
	Double_t Chi2;
	Double_t N3d;  // already applied
	Double_t Iso;
	Double_t TrkPt;
	Double_t TrkIso;
	//Double_t CesWireE2;
	//Double_t CesStripE2;
	Double_t EtCes2nd;
	//Double_t Nvx12;
};

const float fEtMin = 15.0;
const float fXCesMax = 21.0;
const float fEtaMax = 1.1; //my photon is limited to this region


//=================================================================//
void FindXrange(const TH1* h1, const TH1* h2, float& xmin, float& xmax)
{
	for (int bin=0; bin<=h1->GetNbinsX()+1; ++bin)
	{
		if (h1->GetBinContent(bin)>0 || h2->GetBinContent(bin)>0)
		{
			xmin = h1->GetXaxis()->GetBinLowEdge(bin);
			break;
		}
	}
	for (int bin=h1->GetNbinsX()+1; bin>=0; --bin)
	{
		if (h1->GetBinContent(bin)>0 || h2->GetBinContent(bin)>0)
		{
			xmax = h1->GetXaxis()->GetBinUpEdge(bin);
			break;
		}
	}
}

//=================================================================//
//finds the bin number corresponding to a value
//=================================================================//
int GetBin(const TH1* hist, const double val)
{
	assert (hist != NULL && "GetBin::Hist is null!");
	for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
	{
		if (val > hist->GetXaxis()->GetBinLowEdge(bin) 
			  && val <= hist->GetXaxis()->GetBinUpEdge(bin)) return bin;
	}
	std::cout << __FUNCTION__ << "::WARNING! Did not find a bin for Value =" << val
		<< ". Hist bounds are [" << hist->GetXaxis()->GetBinLowEdge(1)
		<< ", " << hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()) 
		<< "] ! returning -1" << std::endl;
	return -1;
}


//=================================================================//
// derives a scale factor for h2 based on h1 when h1,h2 are generated
// with two different Et threshold
//=================================================================//
float GetScale(const TH1* h1, const TH1* h2 , const float fEth, 
					const bool debug=false)
{
	assert (h1 != NULL && "GetScale::h1 is NULL!");
	assert (h2 != NULL && "GetScale::h2 is NULL!");
	assert ( (  fEth >= h1->GetXaxis()->GetBinLowEdge(0) 
				&& fEth <= h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX())) 
			  && "GetScale::Threshold is out of bounds of hists!");

	int bin1 = GetBin(h2,fEth) + 1;
	int binlast = h2->GetNbinsX();

	double N1 = h1->Integral(bin1,binlast);
	double N2 = h2->Integral(bin1,binlast);
	double sc1 = N1/N2;

	if (debug)
	{
		std::cout << "bin1, binlast, fEtThr=" << bin1 << "\t "
			<< binlast << "\t" << fEth << std::endl;
	}

	std::cout << "fEtThr="<< fEth << "::: N1, N2, scale = " 
					<< N1 << "\t" <<  N2 << "\t" 
					<< sc1 << std::endl;
	return sc1;
}


//=================================================================//
bool PassLoosePhoId(const PhoIdVars_& idVars, const float fEtThr)
{
	if (idVars.Et < fEtThr ) return false;
	if (fabs(idVars.EtaDet) > fEtaMax) return false;
	if (idVars.HadEm > 0.125) return false;
	if (idVars.Et<20)
	{
		if (idVars.Iso > 0.15 * idVars.Et) return false;
	} else {
		if (idVars.Iso > 3.0 + 0.02 * (idVars.Et - 20.0)) return false;
	}
	if (idVars.TrkPt > 5.0) return false;
	if (idVars.TrkIso > 5.0) return false;
	
	return true;
}
//=================================================================//
bool PassTightPhoId(const PhoIdVars_& idVars, const float fEtThr)
{
	if (idVars.Et < fEtThr) return false;
	if (fabs(idVars.EtaDet) > fEtaMax) return false;
	
	if (! (idVars.HadEm < 0.125 || idVars.HadEm < 0.055 + 0.00045 * idVars.E) ) return false;
	if (idVars.Et<20)
	{
		if (idVars.Iso > 0.1 * idVars.Et) return false;
	} else {
		if (idVars.Iso > 2.0 + 0.02 * (idVars.Et - 20.0)) return false;
	}
	if (idVars.Chi2 > 20.0) return false;
	if (idVars.TrkPt > 1.0 + 0.005 * idVars.Et) return false;
	if (idVars.TrkIso > 2.0 + 0.005 * idVars.Et) return false;
	//what is saved in the tree is the maximum of wire and strip. so this is good.
	if (idVars.Et<18)
	{
		if (idVars.EtCes2nd > 0.14 * idVars.Et) return false;
	} else {
		if (idVars.EtCes2nd > 2.4 + 0.01 * (idVars.Et - 20.0)) return false;
	}
	return true;
}

//=================================================================//
//// main function
//=================================================================//
void SidebandTest()
{
	const Int_t nThresholds = 9; //different et threholds
	//10,18,40,60,90,120,150,200,300,400,500,600
	//adding +5GeV to avoid any threshold effects
	const float Eth[] = {18.,40.,90.,120.,150.,200.,300.,400.,500.,600.};
   
	unsigned int NtightAboveThres[nThresholds][2],NlooseAboveThres[nThresholds][2];
	for (int i=0; i<nThresholds; ++i) 
	{
		NtightAboveThres[i][0] = 0;		//number of evts in sample i,  above the Et threshol of sample i
		NlooseAboveThres[i][0] = 0;	
		NtightAboveThres[i][1] = 0;		//number of evts in sample i, above the Et threhold of sample i+1
		NlooseAboveThres[i][1] = 0;
	}

	TH1F *loosePhoEt[nThresholds];
	TH1F *tightPhoEt[nThresholds];

	TH1F *avgPt[13][2];
	avgPt[0][0] = new TH1F ("looseavgPt30_45","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,30,45);
	avgPt[1][0] = new TH1F ("looseavgPt45_50","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,45,50);
	avgPt[2][0] = new TH1F ("looseavgPt50_55","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,50,55);
	avgPt[3][0] = new TH1F ("looseavgPt55_65","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,55,65);
	avgPt[4][0] = new TH1F ("looseavgPt65_80","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,65,80);
	avgPt[5][0] = new TH1F ("looseavgPt80_95","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,80,95);
	avgPt[6][0] = new TH1F ("looseavgPt95_110","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,95,110);
	avgPt[7][0] = new TH1F ("looseavgPt110_130","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,110,130);
	avgPt[8][0] = new TH1F ("looseavgPt130_160","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,130,160);
	avgPt[9][0] = new TH1F ("looseavgPt160_220","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,160,220);
	avgPt[10][0] = new TH1F ("looseavgPt220_320","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,220,320);
	avgPt[11][0] = new TH1F ("looseavgPt320_420","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,320,420);
	avgPt[12][0] = new TH1F ("looseavgPt420_550","To get Avg. Pt (loose);E^{#gamma}_{T};Events",1,420,550);

	avgPt[0][1] = new TH1F ("tightavgPt30_45","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,30,45);
	avgPt[1][1] = new TH1F ("tightavgPt45_50","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,45,50);
	avgPt[2][1] = new TH1F ("tightavgPt50_55","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,50,55);
	avgPt[3][1] = new TH1F ("tightavgPt55_65","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,55,65);
	avgPt[4][1] = new TH1F ("tightavgPt65_80","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,65,80);
	avgPt[5][1] = new TH1F ("tightavgPt80_95","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,80,95);
	avgPt[6][1] = new TH1F ("tightavgPt95_110","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,95,110);
	avgPt[7][1] = new TH1F ("tightavgPt110_130","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,110,130);
	avgPt[8][1] = new TH1F ("tightavgPt130_160","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,130,160);
	avgPt[9][1] = new TH1F ("tightavgPt160_220","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,160,220);
	avgPt[10][1] = new TH1F ("tightavgPt220_320","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,220,320);
	avgPt[11][1] = new TH1F ("tightavgPt320_420","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,320,420);
	avgPt[12][1] = new TH1F ("tightavgPt420_550","To get Avg. Pt (tight);E^{#gamma}_{T};Events",1,420,550);




	
	for (Int_t iTh=0; iTh < nThresholds; ++iTh)
	{
		std::stringstream lhistname, thistname;
		lhistname << "SidebandPhotonEt_" << iTh;
		thistname << "TightPhotonEt_" << iTh;

	 	const float  vBins [] = {30,45,50,55,65,80,95,110,130,160,220,320,420,550};
		Int_t Nbins = 13;
		loosePhoEt[iTh] = new TH1F (lhistname.str().c_str(),"Pythia di-jet MC;E_{T}^{#gamma};#gamma tight/#gamma sideband",Nbins,vBins);
		//loosePhoEt[iTh] = new TH1F (lhistname.str().c_str(),"Pythia di-jet MC;E_{T}^{#gamma};#gamma^{tight}/#gamma^{sideband}",200,0,1000);
		loosePhoEt[iTh]->SetLineColor(kRed);
		loosePhoEt[iTh]->SetMarkerColor(kRed);
		loosePhoEt[iTh]->SetMarkerStyle(21);
		tightPhoEt[iTh] = new TH1F (thistname.str().c_str(),"Pythia di-jet MC;E_{T}^{#gamma}; #gamma tight / #gamma sideband", Nbins,vBins);
		//tightPhoEt[iTh] = new TH1F (thistname.str().c_str(),"Pythia di-jet MC;E_{T}^{#gamma};#gamma^{tight}/#gamma^{sideband}",200,0,1000);
		tightPhoEt[iTh]->SetMarkerStyle(22);
		tightPhoEt[iTh]->SetLineColor(kBlue);
		tightPhoEt[iTh]->SetMarkerColor(kBlue);
		tightPhoEt[iTh]->Sumw2();
		loosePhoEt[iTh]->Sumw2();

		

		TChain *myTree = new TChain("mytree");

		if (iTh == 0) //18 GeV
		{
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet18_bq0sqc_p9p11_071709.root");
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet18_bt0sqb_p1p8_071709.root");
		} else if (iTh == 1) //40 GeV
		{
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet40_bq0src_p9p11_071709.root");
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet40_bt0srb_p1p8_071709.root");
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet40_q8is01_p0p10_071709.root");
		} else if (iTh == 2) // 90GeV
		{
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet90_bq0snd_p1p11_071709.root");
		} else if (iTh == 3) // 120GeV
		{
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet120_bt0sub_p0p8_071709.root");
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet120_bq0suc_p9p11_071709.root");
		} else if (iTh == 4) // 150GeV
		{
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet150_bt0svb_p0p8_071709.root");
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet150_bq0svc_p9p11_071709.root");
		} else if (iTh == 5) // 200GeV
		{
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet200_bq0swd_p1p11_071709.root");
		} else if (iTh == 6) // 300GeV
		{
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet300_bq0sxd_p1p11_071709.root");
		} else if (iTh == 7) // 400GeV
		{
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet400_bq0syd_p1p11_071709.root");
		} else if (iTh == 8) // 500GeV
		{
			myTree->Add("/data/nbay02/b/samantha/FlatDiJetMCsamples/myTree_phoISOstudy_PythDiJet500_bq0szd_p1p11_071709.root");
		}

		
		std::stringstream hist_title;
		hist_title << loosePhoEt[iTh]->GetTitle() <<  " (min. P_{T}^{Jet} > "<< Eth[iTh] << " GeV)";
		loosePhoEt[iTh]->SetTitle(hist_title.str().c_str());
		tightPhoEt[iTh]->SetTitle(hist_title.str().c_str());

		
		PhoIdVars_ phoId;

		myTree->SetBranchStatus("*",0);
		myTree->SetBranchStatus("phoID_Et", 1);
		myTree->SetBranchStatus("phoID_Eta", 1);
		myTree->SetBranchStatus("phoID_EtaDet", 1);
		myTree->SetBranchStatus("phoID_HadEm", 1);
		myTree->SetBranchStatus("phoID_Chi2", 1);
		myTree->SetBranchStatus("phoID_ISO", 1);
		myTree->SetBranchStatus("phoID_ISOtrk", 1);
		myTree->SetBranchStatus("phoID_TrkMaxPt", 1);
		myTree->SetBranchStatus("phoID_EtCes2nd", 1);
		myTree->SetBranchStatus("phoID_CesX", 1);
		// this does not seem to have all the variables I need. E/cesZ/ces2ndE wire and strip separately etc?
		// Need to check this?

		myTree->SetBranchAddress("phoID_Et", &phoId.Et);
		myTree->SetBranchAddress("phoID_Eta", &phoId.Eta);
		myTree->SetBranchAddress("phoID_EtaDet", &phoId.EtaDet);
		myTree->SetBranchAddress("phoID_HadEm", &phoId.HadEm);
		myTree->SetBranchAddress("phoID_Chi2", &phoId.Chi2);
		myTree->SetBranchAddress("phoID_ISO", &phoId.Iso);
		myTree->SetBranchAddress("phoID_ISOtrk", &phoId.TrkIso);
		myTree->SetBranchAddress("phoID_TrkMaxPt", &phoId.TrkPt);
		myTree->SetBranchAddress("phoID_EtCes2nd", &phoId.EtCes2nd);
		myTree->SetBranchAddress("phoID_CesX", &phoId.XCes);

		std::cout << ">>>>>>> iTh = " << iTh 
					<<" >>>>>>> Entries = " << myTree->GetEntries() << std::endl;

		double sumEt=0;
		int N = 0;
		for (unsigned int i = 1; i <= myTree->GetEntries(); ++i)
		{
			myTree->GetEntry(i);

			//E is not stored in the tree. need to calculate from Eta and Et.
			//eta = -ln (tan(theta/2))
			float theta = TMath::ATan(2*TMath::Exp(-1 * phoId.Eta));
			phoId.E = phoId.Et / TMath::Sin(theta);
			
			//add 5% to minPt to avoid threshold effects
			bool bPassLoose = PassLoosePhoId(phoId, Eth[iTh] + 5 * Eth[iTh]/100);
			bool bPassTight = PassTightPhoId(phoId, Eth[iTh] + 5 * Eth[iTh]/100);
			if (bPassLoose && !bPassTight) // sideband 
			{
				loosePhoEt[iTh]->Fill(phoId.Et);
				
				//count how many events are above the next threshold for normalization 
				// of the next sample for stitching
				++NlooseAboveThres[iTh][0];
				if (iTh != nThresholds-1)
				{
					if (phoId.Et> (Eth[iTh+1] + 5 * Eth[iTh+1]/100)) ++NlooseAboveThres[iTh][1];
				}
	

			}


				//now fill hists to get avg pt for a bin
				if (iTh == 0 && phoId.Et >= 30 && phoId.Et < 45)
				{
					avgPt[0][0]->Fill(phoId.Et);
					sumEt += phoId.Et;
					++N;

				} else if (iTh == 1)
				{
					if (phoId.Et >= 45 && phoId.Et < 50)      { avgPt[1][0]->Fill(phoId.Et);}
					else if (phoId.Et >= 50 && phoId.Et < 55) { avgPt[2][0]->Fill(phoId.Et);}
					else if (phoId.Et >= 55 && phoId.Et < 65) { avgPt[3][0]->Fill(phoId.Et);}
					else if (phoId.Et >= 65 && phoId.Et < 80) { avgPt[4][0]->Fill(phoId.Et);}
					else if (phoId.Et >= 80 && phoId.Et < 95) { avgPt[5][0]->Fill(phoId.Et);}

				} else if (iTh == 2)
				{  
					if (phoId.Et >= 95 && phoId.Et < 110)       { avgPt[6][0]->Fill(phoId.Et);}
					else if (phoId.Et >= 110 && phoId.Et < 130) { avgPt[7][0]->Fill(phoId.Et);}

				} else if (iTh == 3 && phoId.Et >= 130 && phoId.Et < 160)
				{
					avgPt[8][0]->Fill(phoId.Et);
				} else if (iTh == 4 && phoId.Et >= 160 && phoId.Et < 220)
				{
					avgPt[9][0]->Fill(phoId.Et);
				} else if (iTh == 5 && phoId.Et >= 220 && phoId.Et < 320)
				{
					avgPt[10][0]->Fill(phoId.Et);
				} else if (iTh == 6)
				{
					if (phoId.Et >= 320 && phoId.Et < 420)      { avgPt[11][0]->Fill(phoId.Et);}
					else if (phoId.Et >= 420 && phoId.Et < 550) { avgPt[12][0]->Fill(phoId.Et);}
				}

			

			if (bPassTight) //tight photons
			{
				tightPhoEt[iTh]->Fill(phoId.Et);
				++NtightAboveThres[iTh][0];
				if (iTh != nThresholds-1)
				{
					if (phoId.Et> (Eth[iTh+1] + 5 * Eth[iTh+1]/100)) ++NtightAboveThres[iTh][1];
				}

				//now fill hists to get avg pt for a bin
				if (iTh == 0 && phoId.Et >= 30 && phoId.Et < 45)
				{
					avgPt[0][1]->Fill(phoId.Et);

				} else if (iTh == 1)
				{
					if (phoId.Et >= 45 && phoId.Et < 50)      { avgPt[1][1]->Fill(phoId.Et);}
					else if (phoId.Et >= 50 && phoId.Et < 55) { avgPt[2][1]->Fill(phoId.Et);}
					else if (phoId.Et >= 55 && phoId.Et < 65) { avgPt[3][1]->Fill(phoId.Et);}
					else if (phoId.Et >= 65 && phoId.Et < 80) { avgPt[4][1]->Fill(phoId.Et);}
					else if (phoId.Et >= 80 && phoId.Et < 95) { avgPt[5][1]->Fill(phoId.Et);}

				} else if (iTh == 2)
				{  
					if (phoId.Et >= 95 && phoId.Et < 110)       { avgPt[6][1]->Fill(phoId.Et);}
					else if (phoId.Et >= 110 && phoId.Et < 130) { avgPt[7][1]->Fill(phoId.Et);}

				} else if (iTh == 3 && phoId.Et >= 130 && phoId.Et < 160)
				{
					avgPt[8][1]->Fill(phoId.Et);
				} else if (iTh == 4 && phoId.Et >= 160 && phoId.Et < 220)
				{
					avgPt[9][1]->Fill(phoId.Et);
				} else if (iTh == 5 && phoId.Et >= 220 && phoId.Et < 320)
				{
					avgPt[10][1]->Fill(phoId.Et);
				} else if (iTh == 6)
				{
					if (phoId.Et >= 320 && phoId.Et < 420)      { avgPt[11][1]->Fill(phoId.Et);}
					else if (phoId.Et >= 420 && phoId.Et < 550) { avgPt[12][1]->Fill(phoId.Et);}
				}

			}


		}
		std::cout << "iTh, Avg, N = " << iTh << ", " << sumEt/(1. * N) << "\t " << N << std::endl;
		
		delete myTree;
	}

	
	std::cout << std::setw(8) << "iTh[Et]" << std::setw(12) << "Ntight[iTh/iTh+1]"  << std::endl;
	for (int i=0; i<nThresholds; ++i) 
	{
		std::cout << std::setw(8) << "[" << Eth[i] << "]" << std::setw(12)
			<< NtightAboveThres[i][0] << "/" << NtightAboveThres[i][1] 
			<< std::endl;
		
	}

	//calculate normalization factors
	std::vector<double> vTightScales, vLooseScales;
	for (int i=1; i<nThresholds; ++i) 
	{

		double dScT = 0, dScL = 0;
		if (i==1)
		{
			dScT = (double)NtightAboveThres[i-1][1]/(double) NtightAboveThres[i][0];
			dScL = (double)NlooseAboveThres[i-1][1]/(double) NlooseAboveThres[i][0];
		} else
		{
			dScT = 1 / (double)NtightAboveThres[i][0] * ( vTightScales.at(i-2) * (double)NtightAboveThres[i][1]);
			dScL = 1 / (double)NlooseAboveThres[i][0] * ( vLooseScales.at(i-2) * (double)NlooseAboveThres[i][1]);
		}
		
		vTightScales.push_back(dScT);
		vLooseScales.push_back(dScL);
	}

		
	/// number of scale factors must be one less than the number of samples 
	assert ( (vTightScales.size() == nThresholds-1) && "Scaled factors generated is too little or too much! please check!");
	
	//now normalize each hist to the calculated scales.
	for (int i=1; i<nThresholds; ++i)
	{
		tightPhoEt[i]->Scale(vTightScales.at(i-1));
		loosePhoEt[i]->Scale(vLooseScales.at(i-1));
	}
	
	//for debug only
	std::cout << "after normalization " << std::endl;
	for (int iTh =0; iTh < nThresholds; ++iTh)
	{
		std::cout << " iTh = " << iTh << std::endl;
		for (int bin=1; bin <= tightPhoEt[iTh]->GetNbinsX(); ++bin)
		{
			if (tightPhoEt[iTh]->GetBinContent(bin) || loosePhoEt[iTh]->GetBinContent(bin) )
			{
				if (tightPhoEt[iTh]->GetBinContent(bin))
				{
					std::cout << "low edge / tight val = " << tightPhoEt[iTh]->GetBinLowEdge(bin) 
						<< " / " << tightPhoEt[iTh]->GetBinContent(bin);
				} 
				if (loosePhoEt[iTh]->GetBinContent(bin))
				{
					std::cout << "\tlow edge / loose val = " << loosePhoEt[iTh]->GetBinLowEdge(bin) 
						<< " / " << loosePhoEt[iTh]->GetBinContent(bin);
				}

				std::cout << " " << std::endl;
			}

		}
	}




	
	//now stich them up
	//use first hist and add rest to it.
	//first need to zero out all bins above the next threhold
	for (int bin = 2; bin <=tightPhoEt[0]->GetNbinsX(); ++bin)
	{	
		tightPhoEt[0]->SetBinContent(bin, 0);
		tightPhoEt[0]->SetBinError(bin, 0);
		loosePhoEt[0]->SetBinContent(bin, 0);
		loosePhoEt[0]->SetBinError(bin, 0);
	}
	

	//pick stuff from 40GeV sample
	for (int bin = 2; bin <= 6; ++bin)
	{	
		tightPhoEt[0]->SetBinContent(bin, tightPhoEt[1]->GetBinContent(bin));
		tightPhoEt[0]->SetBinError(bin, tightPhoEt[1]->GetBinError(bin));
		loosePhoEt[0]->SetBinContent(bin, loosePhoEt[1]->GetBinContent(bin) );
		loosePhoEt[0]->SetBinError(bin, loosePhoEt[1]->GetBinError(bin));
	}

	//stuff from 90 sample
	for (int bin = 7; bin <= 8; ++bin)
	{	
		tightPhoEt[0]->SetBinContent(bin, tightPhoEt[2]->GetBinContent(bin));
		tightPhoEt[0]->SetBinError(bin, tightPhoEt[2]->GetBinError(bin));
		loosePhoEt[0]->SetBinContent(bin, loosePhoEt[2]->GetBinContent(bin) );
		loosePhoEt[0]->SetBinError(bin, loosePhoEt[2]->GetBinError(bin));
	}

		//stuff from 120 sample
	for (int bin = 9; bin <= 9; ++bin)
	{	
		tightPhoEt[0]->SetBinContent(bin, tightPhoEt[3]->GetBinContent(bin));
		tightPhoEt[0]->SetBinError(bin, tightPhoEt[3]->GetBinError(bin));
		loosePhoEt[0]->SetBinContent(bin, loosePhoEt[3]->GetBinContent(bin) );
		loosePhoEt[0]->SetBinError(bin, loosePhoEt[3]->GetBinError(bin));
	}

		//stuff from 150 sample
	for (int bin =10; bin <=10; ++bin)
	{	
		tightPhoEt[0]->SetBinContent(bin, tightPhoEt[4]->GetBinContent(bin));
		tightPhoEt[0]->SetBinError(bin, tightPhoEt[4]->GetBinError(bin));
		loosePhoEt[0]->SetBinContent(bin, loosePhoEt[4]->GetBinContent(bin) );
		loosePhoEt[0]->SetBinError(bin, loosePhoEt[4]->GetBinError(bin));
	}

		//stuff from 200 sample
	for (int bin =11; bin <=11; ++bin)
	{	
		tightPhoEt[0]->SetBinContent(bin, tightPhoEt[5]->GetBinContent(bin));
		tightPhoEt[0]->SetBinError(bin, tightPhoEt[5]->GetBinError(bin));
		loosePhoEt[0]->SetBinContent(bin, loosePhoEt[5]->GetBinContent(bin) );
		loosePhoEt[0]->SetBinError(bin, loosePhoEt[5]->GetBinError(bin));
	}
		//stuff from 300 sample
	for (int bin =12; bin <=13; ++bin)
	{	
		tightPhoEt[0]->SetBinContent(bin, tightPhoEt[6]->GetBinContent(bin));
		tightPhoEt[0]->SetBinError(bin, tightPhoEt[6]->GetBinError(bin));
		loosePhoEt[0]->SetBinContent(bin, loosePhoEt[6]->GetBinContent(bin) );
		loosePhoEt[0]->SetBinError(bin, loosePhoEt[6]->GetBinError(bin));
	}

		//stuff from 400 sample
/*	for (int bin =13; bin <=13; ++bin)
	{	
		tightPhoEt[0]->SetBinContent(bin, tightPhoEt[7]->GetBinContent(bin));
		tightPhoEt[0]->SetBinError(bin, tightPhoEt[7]->GetBinError(bin));
		loosePhoEt[0]->SetBinContent(bin, loosePhoEt[7]->GetBinContent(bin) );
		loosePhoEt[0]->SetBinError(bin, loosePhoEt[7]->GetBinError(bin));
	}
	*/
	
	std::cout << "\nwhat was filled to Et[0] hist " << std::endl;
	for (int bin=1; bin <= tightPhoEt[0]->GetNbinsX(); ++bin)
	{
		if (tightPhoEt[0]->GetBinContent(bin))
		{
			std::cout << "low edge / tight val = " << tightPhoEt[0]->GetBinLowEdge(bin) 
				<< " / " << tightPhoEt[0]->GetBinContent(bin);
		}
		if (loosePhoEt[0]->GetBinContent(bin))
		{
			std::cout << "\tlow edge / loose val = " << loosePhoEt[0]->GetBinLowEdge(bin) 
				<< " / " << loosePhoEt[0]->GetBinContent(bin);
		}

		if (tightPhoEt[0]->GetBinContent(bin) || loosePhoEt[0]->GetBinContent(bin) )
		{
			std::cout << " " << std::endl;
		}
	}



	//this is a hack to remove 'nan' values in some bins. I do not know why/how it happens.
/*	for (unsigned int bin=0; bin <= tightPhoEt[0]->GetNbinsX()+1; ++bin)
	{
		if (std::isnan(tightPhoEt[0]->GetBinContent(bin)))
		{
			tightPhoEt[0]->SetBinContent(bin, 0);
			tightPhoEt[0]->SetBinError(bin, 0);
		}
		if (std::isnan(loosePhoEt[0]->GetBinContent(bin)))
		{
			loosePhoEt[0]->SetBinContent(bin, 0);
			loosePhoEt[0]->SetBinError(bin, 0);
		}

	}
*/
	//	tightPhoEt[0]->Print("all");
//	loosePhoEt[0]->Print("all");
	

	//now scale tight to loose and take the ratio

	tightPhoEt[0]->Scale(loosePhoEt[0]->Integral()/(1.0 * tightPhoEt[0]->Integral()) );
	std::cout << "Integrals After normalizing = " << loosePhoEt[0]->Integral() << ", "
		<< tightPhoEt[0]->Integral() << std::endl;
	
	tightPhoEt[0]->Divide(loosePhoEt[0]);
	

	gStyle->SetOptStat(0);
	new TCanvas();
	gPad->SetGridx();
	gPad->SetGridy();
	std::stringstream title;
	title << "Pythia di-jet MC;E_{T}^{#gamma};#gamma^{tight}/#gamma^{sideband};E_{T}^{#gamma}";
	      //<<"- bin size="<< tightPhoEt[0]->GetBinWidth(1) << "GeV";
	tightPhoEt[0]->SetTitle(title.str().c_str());
	//tightPhoEt[0]->SetMinimum(0);
	//tightPhoEt[0]->SetMaximum(2);
	tightPhoEt[0]->SetMarkerColor(kRed);
	tightPhoEt[0]->GetXaxis()->CenterTitle(1);
	tightPhoEt[0]->GetYaxis()->CenterTitle(1);
	tightPhoEt[0]->Draw("PE");
	


	//now make a 
   Int_t n = 12;
	Double_t x[n],y[n], ex[n],ey[n];
	
	for (int i=0; i <12; ++i)
	{
		x[i] = avgPt[i][0]->GetMean();
		ex[i] = avgPt[i][0]->GetRMS();
	}

	assert (tightPhoEt[0]->GetNbinsX() == 13 && "nbins did not match!");
	for (int bin=1; bin < tightPhoEt[0]->GetNbinsX(); ++bin)
	{
		std::cout << bin << "\t" << tightPhoEt[0]->GetBinLowEdge(bin) << "\t" << tightPhoEt[0]->GetBinContent(bin) << std::endl;
		y[bin-1] = tightPhoEt[0]->GetBinContent(bin);
		ey[bin-1] = tightPhoEt[0]->GetBinError(bin);
		ex[bin-1] = tightPhoEt[0]->GetBinWidth(bin)/2;
	}

	gStyle->SetOptFit(1111);
	TGraphErrors *tg = new TGraphErrors(n,x,y,ex,ey);
	tg->SetTitle("Pythia di-jet MC;E_{T}^{#gamma};#gamma^{tight}/#gamma^{sideband}");
	tg->SetMarkerColor(kRed);
	tg->SetLineColor(kBlue);
	tg->SetMarkerStyle(22);
	tg->GetXaxis()->CenterTitle(kTRUE);
	tg->GetYaxis()->CenterTitle(kTRUE);
	new TCanvas();
	tg->Fit("pol1");
	tg->Draw("AP");

	TList *funcList = tg->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1* fitfunc = (TF1*) it->Next();
	std::string func_name("tf1_sideband");
	fitfunc->SetName(func_name.c_str());

	TFile f("DiJet_Sideband_FitFuncForPhoEt.root","UPDATE");
	if ( f.Get(func_name.c_str()) != NULL)
	{
		std::stringstream old_c1;
		old_c1 << func_name << ";1";
		f.Delete(old_c1.str().c_str());
	}
	fitfunc->Write();
	f.ls();
	f.Close();
	
	return;
}
 
