/*This is to study the ISR/FSR effects.
 *A quick test using only the hepg informantion
 *Samantha K. Hewamanage, samantha@fnal.gov, Oct 1, 2009
 */

/* $Id: IsrFsrHepgSyst.C,v 1.6 2009/12/04 04:08:18 samantha Exp $
 * $Log: IsrFsrHepgSyst.C,v $
 * Revision 1.6  2009/12/04 04:08:18  samantha
 * Jay is trying to debug this code to see if I am doing this right.
 * I have added some dumps and minor code shiftting just for this.
 *
 * Revision 1.5  2009/11/15 04:19:41  samantha
 * Added DumpHist to dump debug info at each step of this process.
 *
 * Revision 1.4  2009/11/10 03:21:58  samantha
 * Added anothe sample with pt-hat=100GeV sample. Had moved some input files
 * to a new location, so had to update the paths. Deleted the stuff I had for
 * trying to use the data pho et specturm as base. Results are in
 * https://hep04.baylor.edu/hep/samantha/1455.
 *
 * Revision 1.3  2009/11/06 17:55:56  samantha
 * This version stiches two pt-hat samples, 8GeV and 50GeV samples. New ROOT files
 * have a tree with kinematic variables to derive the normalization need for
 * stiching and avoid binning effects. Added two new functions to fill out a
 * histogram for a given pt-hat sample and to cound the events above a certain
 * threshols (needed for deriving the normalization scale.) Combined result is in
 * https://hep04.baylor.edu/hep/samantha/1426. The input ROOT files are scattered
 * around the disks.
 *
 * Revision 1.2  2009/10/21 19:09:40  samantha
 * Three histograms, PhotonEt, Lead jet Et and InvMass(pho,lead jet) are made with
 * variable bin sizes (same as in the final hists). Generated results are in
 * https://hep.baylor.edu/hep/samantha/1404 with pt-hat=8GeV and photon
 * pt-min=22GeV sample. With pt-hat=50GeV and photon min-pt=22GeV results are in
 * https://hep.baylor.edu/hep/samantha/1419.
 * Only changed made to generate the second (1419) result was chaning the input
 * files.
 *
 * Revision 1.1  2009/10/14 03:04:29  samantha
 * This is to study the ISR/FSR uncertainty of the signal sample (pho+jets MC).
 * Results of this initial version is in https://hep.baylor.edu/hep/samantha/1368.
 * Still doing a contant bin size histograms. Will change in the next version.
 *
 */

#include <iostream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLegend.h"
#include <sstream>
#include "TMath.h"
#include <iomanip>
#include "TFile.h"
#include "TH1.h"
#include <vector>
#include "TTree.h"
#include "TSystem.h"
#include "CommonTools.hh"
#include "TF1.h"
#include <memory>

using namespace std;

// ToRun:
// Load ~/samantha/RootTools/CommonTools_C.so
// .L IsrFsrHepgSyst.C+

void MakeIndividualPtHatPlot(TH1F **h, const int iNHIST)
{

	TH1* hh[iNHIST];
	TCanvas *c1 = new TCanvas();
	c1->Divide(2,1);
	c1->cd(1);
	gPad->SetLogy();
	for (int ihist=0; ihist<iNHIST; ++ihist)
	{
		hh[ihist] = dynamic_cast<TH1*> (h[ihist]->Clone("copy"));
		if (ihist==0)
		{
			c1->cd(1);
			hh[ihist]->DrawClone();
		} else
		{
			c1->cd(1);
			hh[ihist]->DrawClone("SAME");
			hh[ihist]->Scale(hh[0]->Integral()/ (1. * hh[ihist]->Integral()));
			hh[ihist]->Divide(hh[0]);
			if (ihist==1)
			{
				c1->cd(2);
				hh[ihist]->Draw("PE");
			} 	else 
			{
				c1->cd(2);
				hh[ihist]->Draw("PE SAME");
			}
		}
	}
	c1->cd();

}

//-----------------------------------------------------------------------------
// effects: normalize each bin by the bin width
// side effect: normalizes histogram bins in hist
// guarantee: no-throw
// requires: hist != NULL
// requires: hist->GetDimension() == 1
//-----------------------------------------------------------------------------
void norm_bins (TH1 *hist)
{
	assert (hist != NULL && "requirement failed"); //spec
	assert (hist->GetDimension() == 1 && "requirement failed"); //spec

	for (unsigned bin = 1; bin <= unsigned (hist->GetNbinsX()); ++ bin)
	{
		const float width = hist->GetXaxis()->GetBinWidth (bin);
		hist->SetBinContent (bin, hist->GetBinContent (bin) / width);
		hist->SetBinError (bin, hist->GetBinError (bin) / width);
	};
};


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
//finds the number of events above a threhold for given cut
//=================================================================//
double NEvtsAbove(TFile *file, const std::string sBrName, TH1* hist , 
		const float fLoCut , const float fHiCut , double &count1, double &count2)
{
	assert (hist != NULL && "HIST is null");
	TTree *tree = (TTree*) file->Get("/Ana/PhoJetsHepgSyst/Hist/HistValuesTree/PhoJetsHepgSyst_FlatTree");
	assert (tree != NULL && "tree is null");
	if (fHiCut>0)
	{
		assert (fLoCut<fHiCut && "Lower cut off value must less then upper cut off");
	}

	// init to be safe
	count1 =0, count2=0;
	
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus(sBrName.c_str(),1);
	float fVal=0;
	tree->SetBranchAddress(sBrName.c_str(), &fVal);
	double count = 0;
	for (int i=0; i < tree->GetEntries(); ++i)
	{
			tree->GetEntry(i);
			//std::cout << i << "\t" << fVal << std::endl;
			if (fVal>fLoCut)
			{
				hist->Fill(fVal);
				++count1;
				//no counting for the last pt-hat sample above the next sample threshold
				if (fHiCut>0 && fVal>fHiCut) ++count2;
			}
	}
	std::cout << "Branch " << sBrName << " = " << count1  << " [cut>" 
				<< fLoCut << "]\t" << count2 << "[cut>" << fHiCut << "]" << std::endl;
	delete tree;
	return count;
}

//=================================================================//
//// main function
//=================================================================//
void IsrFsrHepgSyst(
		const std::string brname, const std::string title, 
		const float xmin, 
		const float xpoint1, const float xpoint2, 
		const float xpoint3, const float xpoint4, 
		const float width1,  const float width2, 
		const float width3,  const float width4,
		const float cut1, const float cut2, const float cut3, const float cut4, const float cut5,
		const  int tweak, const std::string eps_name=""
		)
{
	const bool bMakeIndividualPtHatPlot = false;

	if (gSystem->Load("~/samantha/RootTools/CommonTools_cc.so") == 0)
	{
		std::cout << "Required CommonTools_cc.so did not load. please check!" <<std::endl;
	}
	
	const int iNTH = 5; //stiching 3 samples pt-hat=20GeV and pt-hat=60GeV ,pt-hat=100GeV,210GeV
	const float fEth[] = {cut1,cut2,cut3,cut4, cut5};
	const int iNHIST = 5;	

	TFile* file[iNHIST][iNTH];

	file[0][0] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat20_BaseResults.root");
	file[1][0] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat20_MoreISRResults.root");
	file[2][0] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat20_LessISRResults.root");
	file[3][0] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat20_MoreFSRResults.root");
	file[4][0] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat20_LessFSRResults.root");

	file[0][1] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat60_BaseResults.root");
	file[1][1] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat60_MoreISRResults.root");
	file[2][1] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat60_LessISRResults.root");
	file[3][1] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat60_MoreFSRResults.root");
	file[4][1] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat60_LessFSRResults.root");

	file[0][2] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat100_BaseResults.root");
	file[1][2] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat100_MoreISRResults.root");
	file[2][2] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat100_LessISRResults.root");
	file[3][2] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat100_MoreFSRResults.root");
	file[4][2] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat100_LessFSRResults.root");

	file[0][3] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat160_BaseResults.root");
	file[1][3] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat160_MoreISRResults.root");
	file[2][3] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat160_LessISRResults.root");
	file[3][3] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat160_MoreFSRResults.root");
	file[4][3] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat160_LessFSRResults.root");

	file[0][4] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat210_BaseResults.root");
	file[1][4] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat210_MoreISRResults.root");
	file[2][4] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat210_LessISRResults.root");
	file[3][4] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat210_MoreFSRResults.root");
	file[4][4] = new TFile("~/RESULTS/10212009_ISRFSR/11122209_AllFinalHists/IsrFsrPythiaPtHat210_LessFSRResults.root");

	for (int i=0; i<iNHIST; ++i) //5 hists base/ISRup/ISRdown/FSRup/FSRdown
	{
		for (int iTh=0; iTh<iNTH; ++iTh)
		{
			std::stringstream fname;
			fname << "File " << file[i][iTh]->GetName() << " is not found!";
			assert (file[i][iTh] != NULL && fname.str().c_str());
		}
	}


	TH1F* h1[iNHIST];		//pt-hat 20 hists
	TH1F* h2[iNHIST];		//pt-hat 60 hists
	TH1F* h3[iNHIST];		//pt-hat 100 hists
	TH1F* h4[iNHIST];		//pt-hat 160 hists
	TH1F* h5[iNHIST];		//pt-hat 210 hists


	double nEvtsAbove[iNHIST][iNTH][2];			//[iNTH][0] = entries in iNTH sample [iNTH][1] entries in iNTH sample above the threshold of iNTH+1 sample
	for (int iNHIST=0; iNHIST>iNHIST; ++iNHIST)
	{
		for (int iNTH=0; iNTH<iNTH; ++iNTH) 
		{	
			for (int i=0; i<2; ++i)	nEvtsAbove[iNHIST][iNTH][i] = 0;
		}
	}
	
	for (int i=0; i<iNHIST; ++i)
	{
		std::stringstream name1, name2, name3, name4, name5;
		name1 << "histPt20_"<< i;
		name2 << "histPt60_"<< i;
		name3 << "histPt100_"<< i;
		name4 << "histPt160_"<< i;
		name5 << "histPt210_"<< i;
		std::stringstream title1, title2, title3, title4, title5;
		if (i==0) 
		{
			title1 << "Photon Et"<< fEth[0] << "_base";
			title2 << "Photon Et"<< fEth[1] << "_base";
			title3 << "Photon Et"<< fEth[2] << "_base";
			title4 << "Photon Et"<< fEth[3] << "_base";
			title5 << "Photon Et"<< fEth[4] << "_base";
		} else if (i==1) 
		{
			title1 << "Photon Et"<< fEth[0] << "_moreISR";
			title2<< "Photon Et"<< fEth[1] << "_moreISR";
			title3<< "Photon Et"<< fEth[2] << "_moreISR";
			title4<< "Photon Et"<< fEth[3] << "_moreISR";
			title5<< "Photon Et"<< fEth[4] << "_moreISR";
		} else if (i==2)
		{
			title1 << "Photon Et"<< fEth[0] << "_lessISR";
			title2<< "Photon Et"<< fEth[1] << "_lessISR";
			title3<< "Photon Et"<< fEth[2] << "_lessISR";
			title4<< "Photon Et"<< fEth[3] << "_lessISR";
			title5<< "Photon Et"<< fEth[4] << "_lessISR";
		} else if (i==3) 
		{
			title1 << "Photon Et"<< fEth[0] << "_moreFSR";
			title2<< "Photon Et"<< fEth[1] << "_moreFSR";
			title3<< "Photon Et"<< fEth[2] << "_moreFSR";
			title4<< "Photon Et"<< fEth[3] << "_moreFSR";
			title5<< "Photon Et"<< fEth[4] << "_moreFSR";
		} else if (i==4) 
		{
			title1 << "Photon Et"<< fEth[0] << "_lessFSR";
			title2<< "Photon Et"<< fEth[1] << "_lessFSR";
			title3<< "Photon Et"<< fEth[2] << "_lessFSR";
			title4<< "Photon Et"<< fEth[3] << "_lessFSR";
			title5<< "Photon Et"<< fEth[4] << "_lessFSR";
		}
		
  		h1[i] = MakeVariableBinHist (name1.str(), title1.str(), xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
  		h2[i] = MakeVariableBinHist (name2.str(), title2.str(), xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
  		h3[i] = MakeVariableBinHist (name3.str(), title3.str(), xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
  		h4[i] = MakeVariableBinHist (name4.str(), title4.str(), xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
  		h5[i] = MakeVariableBinHist (name5.str(), title5.str(), xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
		h1[i]->Sumw2();
		h2[i]->Sumw2();
		h3[i]->Sumw2();
		h4[i]->Sumw2();
		h5[i]->Sumw2();

		
		h1[i]->SetLineColor(i+1);
		
		NEvtsAbove(file[i][0], brname,h1[i], fEth[0], fEth[1], nEvtsAbove[i][0][0], nEvtsAbove[i][0][1]);
		NEvtsAbove(file[i][1], brname,h2[i], fEth[1], fEth[2], nEvtsAbove[i][1][0], nEvtsAbove[i][1][1]);
		NEvtsAbove(file[i][2], brname,h3[i], fEth[2], fEth[3], nEvtsAbove[i][2][0], nEvtsAbove[i][2][1]);
		NEvtsAbove(file[i][3], brname,h4[i], fEth[3], fEth[4], nEvtsAbove[i][3][0], nEvtsAbove[i][3][1]);
		NEvtsAbove(file[i][4], brname,h5[i], fEth[4],   -1.  , nEvtsAbove[i][4][0], nEvtsAbove[i][4][1]);

		norm_bins(h1[i]);
		norm_bins(h2[i]);
		norm_bins(h3[i]);
		norm_bins(h4[i]);
		norm_bins(h5[i]);
		//DrawClone(h1[i]);

	}

	if (bMakeIndividualPtHatPlot)
	{
		MakeIndividualPtHatPlot(h1,iNHIST);
		MakeIndividualPtHatPlot(h2,iNHIST);
		MakeIndividualPtHatPlot(h3,iNHIST);
		MakeIndividualPtHatPlot(h4,iNHIST);
		MakeIndividualPtHatPlot(h5,iNHIST);
	}

	//calculate normalization factors
	double vScales[iNHIST][iNTH-1];
	for (int ihist=0; ihist<iNHIST; ++ihist)
	{ 
		for (int nth=1; nth<iNTH; ++nth) 
		{
			double dScT=0;
			if (nth==1)
			{
				dScT = nEvtsAbove[ihist][nth-1][1] / nEvtsAbove[ihist][nth][0];
			} else
			{
				
				dScT = vScales[ihist][nth-2]  * ( nEvtsAbove[ihist][nth-1][1] / nEvtsAbove[ihist][nth][0]);
			}

			vScales[ihist][nth-1] = dScT;
		}
	}

	
	for (int ihist=0; ihist<iNHIST; ++ihist)
	{ 
		h2[ihist]->Scale(vScales[ihist][0]);
		h3[ihist]->Scale(vScales[ihist][1]);
		h4[ihist]->Scale(vScales[ihist][2]);
		h5[ihist]->Scale(vScales[ihist][3]);
	}
			
	//pick contributions from each pt-hat samples
	//
	int picbin2 = h2[0]->GetMaximumBin();
	int picbin3 = h3[0]->GetMaximumBin();
	int picbin4 = h4[0]->GetMaximumBin();
	int picbin5 = h5[0]->GetMaximumBin();

		for (int ihist=0; ihist <iNHIST; ++ihist)
		{
			//zero out rest of the bins before picking up contributions 
			//from higher pt-hat samples
			for (int bin=picbin2; bin <= h1[ihist]->GetNbinsX(); ++bin)
			{
				h1[ihist]->SetBinContent(bin,0);
				h1[ihist]->SetBinError(bin,0);
			}

			for (int bin=picbin2; bin <picbin3; ++bin)
			{
				h1[ihist]->SetBinContent(bin, h2[ihist]->GetBinContent(bin));
				h1[ihist]->SetBinError(bin, h2[ihist]->GetBinError(bin));
			}

			for (int bin=picbin3; bin <picbin4; ++bin)
			{
				h1[ihist]->SetBinContent(bin, h3[ihist]->GetBinContent(bin));
				h1[ihist]->SetBinError(bin, h3[ihist]->GetBinError(bin));
			}
			for (int bin=picbin4; bin < picbin5; ++bin)
			{
				h1[ihist]->SetBinContent(bin, h4[ihist]->GetBinContent(bin));
				h1[ihist]->SetBinError(bin, h4[ihist]->GetBinError(bin));
			}
			for (int bin=picbin5; bin <= h1[ihist]->GetNbinsX(); ++bin)
			{
				h1[ihist]->SetBinContent(bin, h5[ihist]->GetBinContent(bin));
				h1[ihist]->SetBinError(bin, h5[ihist]->GetBinError(bin));
			}

		}

			//now normalize all variation to base
	for (int i=1; i<iNHIST; ++i)
	{
		h1[i]->Scale(h1[0]->Integral()/ (1. * h1[i]->Integral()));
		if (i==1) 
		{
			//DumpHist(h1[1],"====== MoreISR sample ====== after normalizing to base sample bby integral");
			std::cout << "Integral of base sample     = " << h1[0]->Integral() << std::endl;
			std::cout << "Integral of more ISR sample = " << h1[1]->Integral() << std::endl;
		}
		//h1[i]->Divide(h1[0]);
	}	

	std::vector<float> base_val, base_err, lessFsr_val, lessFsr_err, final_val, final_err;
	std::vector<float> moreFsr_val, moreFsr_err;
	std::vector<float> moreIsr_val, moreIsr_err;
	std::vector<float> lessIsr_val, lessIsr_err;
	std::vector<float> maxIsr_val, maxIsr_err;
	std::vector<float> maxFsr_val, maxFsr_err;

	for (int bin=1; bin < h1[0]->GetNbinsX(); ++bin) 
	{
		base_val.push_back(h1[0]->GetBinContent(bin));
		base_err.push_back(h1[0]->GetBinError(bin));
		
		moreIsr_val.push_back(h1[1]->GetBinContent(bin));
		moreIsr_err.push_back(h1[1]->GetBinError(bin));
		lessIsr_val.push_back(h1[2]->GetBinContent(bin));
		lessIsr_err.push_back(h1[2]->GetBinError(bin));
		
		moreFsr_val.push_back(h1[3]->GetBinContent(bin));
		moreFsr_err.push_back(h1[3]->GetBinError(bin));
		lessFsr_val.push_back(h1[4]->GetBinContent(bin));
		lessFsr_err.push_back(h1[4]->GetBinError(bin));
	}


	//for debugging only. delete this for loop
	//and uncomment the division in line 451
	for (int i=1; i<iNHIST; ++i)
	{
		h1[i]->Divide(h1[0]);
	}
	
	//DumpHist(h1[1],"====== MoreISR sample ====== after dividing by base sample");
	
	h1[0]->SetLineColor(2);
	h1[1]->SetLineColor(4);
	h1[2]->SetLineColor(6);
	h1[3]->SetLineColor(8);
	h1[4]->SetLineColor(28);

	h1[0]->SetMarkerColor(2);
	h1[1]->SetMarkerColor(4);
	h1[2]->SetMarkerColor(6);
	h1[3]->SetMarkerColor(8);
	h1[4]->SetMarkerColor(28);

	h1[0]->SetMarkerStyle(20);
	h1[1]->SetMarkerStyle(22);
	h1[2]->SetMarkerStyle(23);
	h1[3]->SetMarkerStyle(22);
	h1[4]->SetMarkerStyle(23);

	
	

	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas();
	gPad->SetGridx();
	gPad->SetGridy();
	h1[1]->SetTitle(title.c_str());
	h1[1]->SetMinimum(0.75);
	h1[1]->SetMaximum(1.1);
	for (int i=1; i<5; ++i)
	{
		if (i==1) h1[i]->Draw();
		else h1[i]->Draw("same");
	}

	TLegend *leg1 = new TLegend (0.8,0.7,0.99,.99);
	leg1->SetTextFont(42);
	leg1->SetTextSize(0.03);

	//leg1->AddEntry(h1[0],"Base Sample");
	leg1->AddEntry(h1[1],"More ISR");
	leg1->AddEntry(h1[2],"Less ISR");
	leg1->AddEntry(h1[3],"More FSR");
	leg1->AddEntry(h1[4],"Less FSR");
	leg1->Draw();



	std::stringstream c1_name;
	c1_name << "Fit_" << brname;
	TCanvas *c1 = new TCanvas(c1_name.str().c_str(),c1_name.str().c_str());
  	TH1F *final = MakeVariableBinHist ("ISRFSR_final", title, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);

	for (int bin=1; bin < final->GetNbinsX(); ++bin)
	{
		float y1_p = h1[1]->GetBinContent(bin);
		float y1_m = h1[2]->GetBinContent(bin);
		float y2_p = h1[3]->GetBinContent(bin);
		float y2_m = h1[4]->GetBinContent(bin);


		float dy1_p = h1[1]->GetBinError(bin);
		float dy1_m = h1[2]->GetBinError(bin);
		float dy2_p = h1[3]->GetBinError(bin);
		float dy2_m = h1[4]->GetBinError(bin);


		float dy1 = 0, dy2 =0;
		float e_dy1 = 0, e_dy2 =0;
		
		if (y1_p>0) // must to this do avoid adding errors to bin with no values
		{
			if ( (y1_p > 1 && y1_m > 1) || (y1_p < 1 && y1_m < 1) )
			{
				dy1 = (y1_p+y1_m)/2;
				dy1 = fabs (dy1 - 1);			
				e_dy1 = 0.5 * sqrt ( pow(dy1_p,2) + pow(dy1_m,2) );
			} else 
			{
				dy1 = std::max(fabs(y1_p-1), fabs(y1_m-1));

				if (dy1 == fabs(y1_p-1)) e_dy1 = dy1_p;
				else e_dy1 = dy1_m;
			}

			if ( (y2_p > 1 && y2_m > 1) || (y2_p < 1 && y2_m < 1) )
			{
				dy2 = (y2_p+y2_m)/2;
				dy2 = fabs (dy2 - 1);			
				e_dy2 = 0.5 * sqrt ( pow(dy2_p,2) + pow(dy2_m,2) );
			} else 
			{
				dy2 = std::max(fabs(y2_p-1), fabs(y2_m-1));
				if (dy2 == fabs(y2_p-1)) e_dy2 = dy2_p;
				else e_dy2 = dy2_m;
			}

			
			//std::cout << "Y1:: " << y1_p << "\t" << y1_m << "\t" << dy1 << std::endl;	
			//std::cout << "Y2:: " << y2_p << "\t" << y2_m << "\t" << dy2 << std::endl;	
			maxIsr_val.push_back(dy1);
			maxIsr_err.push_back(e_dy1);
			maxFsr_val.push_back(dy2);
			maxFsr_err.push_back(e_dy2);
			
			float val = sqrt(pow(dy1,2) + pow(dy2,2));
			float err = sqrt( ( pow(dy1 * e_dy1, 2) + pow(dy2 * e_dy2,2) ) / val );  //using simple error propagation
			final->SetBinContent(bin,val);
			final->SetBinError(bin,err);
			final_val.push_back(val);
			final_err.push_back(err);
		 } else 
		{
			maxIsr_val.push_back(0);
			maxIsr_err.push_back(0);
			maxFsr_val.push_back(0);
			maxFsr_err.push_back(0);
			final_val.push_back(0);
			final_err.push_back(0);

		 }
		

	}
	

	std::cout << red <<std::setw(5) << "bin" << std::setw(7) << "loEdge"
			<< std::setw(10) << "base val/err" 
			<< std::setw(15) << "moreIsr val/err" 
			<< std::setw(15) << "lessIsr val/err" 
			<< std::setw(15) << "moreFsr val/err" 
			<< std::setw(15) << "lessFsr val/err" 
			<< std::setw(15) << "lessFsr final val/err" 
			<< clearatt<< std::endl;

	for (int i=0; i < base_val.size(); ++i)
	{

		//float my_err = lessFsr_val.at(i)/base_val.at(i) * sqrt ( pow(lessFsr_err.at(i)/lessFsr_val.at(i),2) + pow(base_err.at(i)/base_val.at(i),2) );

	std::cout << blue << std::setw(5) << i << std::setw(7) << h1[0]->GetBinLowEdge(i+1)
			<< std::setw(10) << base_val.at(i) << " / " << base_err.at(i) 
			<< std::setw(10) << moreIsr_val.at(i) << " / " <<  moreIsr_err.at(i) 
			<< std::setw(10) << lessIsr_val.at(i) << " / " <<  lessIsr_err.at(i) 
			<< std::setw(10) << moreFsr_val.at(i) << " / " <<  moreFsr_err.at(i) 
			<< std::setw(10) << lessFsr_val.at(i) << " / " <<  lessFsr_err.at(i) 
			<< std::setw(14) << final_val.at(i) << " / " << final_err.at(i) 
			<< clearatt << std::endl;

	

	}


	
	if (tweak)
	{
		//tweak crossing low end min pointto improve fitting
		for (int bin=1; bin <final->GetNbinsX(); ++bin)
		{
			//ignore first zero bins
			if ( ! (final->GetBinContent(bin)) ) continue;

			//now find the first minimum to tweak
			if (final->GetBinContent(bin) > final->GetBinContent(bin+1) ) continue;
			final->SetBinContent(bin, (final->GetBinContent(bin-1)+final->GetBinContent(bin+1))/2 );
			break;
		}
	}
	
	//find fit range
	int xminbin=0, xmaxbin=0;
	FindNonZeroXrange(final,xminbin,xmaxbin);
	std::cout << "minx, maxx = " << xminbin << "\t" << xmaxbin << std::endl;
	
	final->Fit("pol1","","",final->GetBinCenter(xminbin),final->GetBinCenter(xmaxbin));


	if (eps_name.length())
	{
		c1->Print(eps_name.c_str());
		std::stringstream rationame;
		int dot = eps_name.find(".");
		std::string subs= eps_name.substr(0,dot);
		rationame << subs << "_ratio.eps";
		//c->Print(rationame.str().c_str());
	}




	TList *funcList = final->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1* fitfunc = (TF1*) it->Next();
	fitfunc->SetName(brname.c_str());
	//check if there is a saved version already.
	//if so delete that before writing the new one.
	//otherwise there will be many copied of the same obj
	//with different cycles
	TFile f("ISRFSR_FitFunctions.root","UPDATE");
	if ( f.Get(brname.c_str()) != NULL)
	{
		std::stringstream old_tf;
		//old_tf << brname << ";*";  //don't do this if you want to see the final fit in canvas.
											//this deletes objects in disk AND memory.
		old_tf << brname << ";1";
			f.Delete(old_tf.str().c_str());
	}
	fitfunc->Write();
	if ( f.Get(c1_name.str().c_str()) != NULL)
	{
		std::stringstream old_c1;
		old_c1 << c1_name.str() << ";1";
		f.Delete(old_c1.str().c_str());
	}

	c1->Write();
	f.ls();
	f.Close();
	


	

		for (int ihist=0; ihist <iNHIST; ++ihist)
		{
		/*	delete h1[ihist];
			delete h2[ihist];
			delete h3[ihist];
			delete h4[ihist];
			delete h5[ihist];
		*/
			}


		for (int ihist=0; ihist <iNHIST; ++ihist)
		{
			for (int nth=0; nth <iNTH; ++nth)
			{
		//		delete file[ihist][nth];
			}
		}


	
	return;
}

void IsrFsrHepgSyst(int jets, const std::string which)
{	

	if (jets == 1)
	{
		if (which == "PhoEt")
		{
			IsrFsrHepgSyst("pj1_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};Maximum % ISR/FSR error",
					 30,150,250,300,650, 10,20,50,250,
					 //30,80,120,190,250,
					 30,31,32,33,34,
					 0,
					 "IsrFsr_pj1_Et_pho.eps"
					 );
		}

		if (which == "JetEt")
		{
			IsrFsrHepgSyst("pj1_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};Maximum % ISR/FSR error"
					,10,200,250,300,600, 10,10,50,200,
					 //15,70,120,180,240,
					 30,31,32,33,34,
					0,
					 "IsrFsr_pj1_Et_leadjet.eps"
					);
		}
		if (which == "InvMass")
		{
			IsrFsrHepgSyst("pj1_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet):; Inv. Mass(#gamma,lead jet);Maximum % ISR/FSR error"
					, 50,500, 600, 700, 1000, 25,50,75,300,
					//50,150,250,400,500,
					 30,31,32,33,34,
					0,
					 "IsrFsr_pj1_InvM_phojet.eps"
					);
		}

		if (which == "NJet")
		{
    		IsrFsrHepgSyst("pj1_Njet15",
				"Jet Multiplicity (E_{T}>15GeV)",
				1,2,3,4,15, 1,1,1,1,
				30,31,32,33,34,
				0,
					 "IsrFsr_pj1_njet15.eps"
				);

		}

		if (which == "Ht")
		{
			
			IsrFsrHepgSyst("pj1_Ht",
				"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);Maximum % ISR/FSR error",
				0,300, 600, 800, 1200, 20,50,100,300,
					//40,160,240,400,500,
					 30,31,32,33,34,
				0,
					 "IsrFsr_pj1_Ht.eps"
				);

		}


	}

		
	if (jets == 2)
	{

		if (which == "PhoEt")
		{
			IsrFsrHepgSyst("pj2_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};Maximum % ISR/FSR error",
					30,150,250,300,650, 10,20,50,250,
				//30,70,130,190,250,
					 30,31,32,33,34,
					0,
					 "IsrFsr_pj2_Et_pho.eps"
					);
		}

		if (which == "Jet1Et")
		{
			IsrFsrHepgSyst("pj2_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};Maximum % ISR/FSR error",
					10,200,260,300,600, 10,20,40,200,
				//10,60,100,160,220,
					 30,31,32,33,34,
					0,
					 "IsrFsr_pj2_Et_leadJet.eps"
					);
		}



		if (which == "InvMass_pj1j2")
		{
			IsrFsrHepgSyst("pj2_InvM_pj1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet, subleading jet):; Inv. Mass(#gamma,lead jet, subleading jet);Maximum % ISR/FSR error",
					100, 200, 300, 700, 1300, 10,20,50,200,
				//100,200,300,450,600,
					 30,31,32,33,34,
					0,
					 "IsrFsr_pj2_InvM_pj1j2.eps"
					);
		}
		if (which == "InvMass_pj1")
		{
			IsrFsrHepgSyst("pj2_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet); Inv. Mass(#gamma,lead jet);Maximum % ISR/FSR error",
					100, 200, 300, 700, 1000, 10,20,50,200,
					//100,150,240,400,500,
					 30,31,32,33,34,
					0,
					 "IsrFsr_pj2_InvM_pj1.eps"
					);
		}
		if (which == "TwoJetInvMass")
		{
			IsrFsrHepgSyst("pj2_InvM_j1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(lead jet, subleading jet):; Inv. Mass(lead jet, subleading jet);Maximum % ISR/FSR error",
				 	100, 200, 300, 700, 1200, 10,20,50,400,
				//30,50,60,70,90,
					 30,31,32,33,34,
					0,
					 "IsrFsr_pj2_InvM_twoJets.eps"
					);
		}

		if (which == "Ht")
		{

			IsrFsrHepgSyst("pj2_Ht",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);Maximum % ISR/FSR error",
					0,300, 600, 700, 1200, 20,50,100,400,
					//30,31,180,400,550,
					 30,31,32,33,34,
					0,
					 "IsrFsr_pj2_Ht.eps"
					);

		}
	}
 
}

void IsrFsrHepgSyst(int i)
{
	if (i==1) IsrFsrHepgSyst(1,"PhoEt");
	if (i==2) IsrFsrHepgSyst(1,"JetEt");
	if (i==3) IsrFsrHepgSyst(1,"InvMass");
	if (i==4) IsrFsrHepgSyst(1,"Ht");
	//if (i==4) IsrFsrHepgSyst(1, "NJet"); //need to change the branch memory address type to 'int' for this
	if (i==5) IsrFsrHepgSyst(2,"PhoEt");
	if (i==6) IsrFsrHepgSyst(2,"Jet1Et");
	if (i==7) IsrFsrHepgSyst(2,"InvMass_pj1");
	if (i==8) IsrFsrHepgSyst(2,"InvMass_pj1j2");
	if (i==9) IsrFsrHepgSyst(2,"TwoJetInvMass");
	if (i==10) IsrFsrHepgSyst(2,"Ht");
}

void IsrFsrHepgSyst()
{
	for (int i=1; i<=10; ++i) IsrFsrHepgSyst(i);
}

