/*This is to study the Q2 effects.
 *A quick test using only the hepg informantion
 *Samantha K. Hewamanage, samantha@fnal.gov, Oct 14, 2009
 */

/* $Id: MakeQ2HepgSystHists.C,v 1.4 2010/02/18 16:59:05 samantha Exp $
 * $Log: MakeQ2HepgSystHists.C,v $
 * Revision 1.4  2010/02/18 16:59:05  samantha
 * MODIFIED: 1. Mainly the histograms are binned the same way as the final plots
 * are.
 * 	2. The input root files were relocated and hence path is updated.
 * ADDED:	1.Code is modified to autmatically pick the bin with the maximum bin
 *         content to make the merging smooth. Often if I pick a bin on the left of
 * 	the peak, things become spiky.
 * 	2. This version makes the final sysytematic plots for the Q2
 * 	uncertainty. It picks the maximum deviation from the nominal for both up
 * 	and down samples and uses to derive a percentage uncertainty.
 * 	3. A fit function is derived for each distribution to represent the
 * 	systematic uncerainty and saved to a root file. The name of the function
 * 	is derived from the vaiable's branch name.
 * RENAMED: NaboveThres -> nEvtsAboveThr
 *
 * Revision 1.2  2009/11/04 19:55:50  samantha
 * This version tries to stich many pt-hat samples to improve stats at high pt.
 * But there seems to be a bug as the results does not seem to show much
 * improvement
 * in stats. Looking at the log plot, it shows I have 1% precision which translates
 * to ~2% in ratio plot. I see more than 2% stat error in the plot.
 * see elog#https://hep.baylor.edu/hep/samantha/1443. I am goind to DEBUG this now.
 *
 * Revision 1.1  2009/10/31 17:51:42  samantha
 * This is the simple first version to make a ratio plots of Q^2 systematic
 * variations wrt to a base sample. I did not have a base sample matching the
 * pt threhold so I used a base sample with a different threhold and the result
 * is in elog#https://hep.baylor.edu/hep/samantha/1407
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
#include "CommonTools.hh"
#include "TF1.h"

using namespace std;

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
//finds the number of events above a threhold for given cut
//=================================================================//
void NEvtsAbove(TFile *file, const std::string sBrName, TH1* hist , 
		const float fLoCut , const float fHiCut , 
		double &count1, double &count2,
		const bool debug = false
		)
{

	assert (hist != NULL && "HIST is null");
	assert (file != NULL && "file is null");
	TTree *tree = (TTree*) file->Get("/Ana/PhoJetsHepgSyst/Hist/HistValuesTree/PhoJetsHepgSyst_FlatTree");
	assert (tree != NULL && "tree is null");
	if (fHiCut>0)
	{
		assert (fLoCut<fHiCut && "Lower cut off value must less then upper cut off");
	}
	if (debug) std::cout << "Using file " << file->GetName() << std::endl;

	// init to be safe
	count1 =0, count2=0;
	
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus(sBrName.c_str(),1);
	float fVal=0;
	tree->SetBranchAddress(sBrName.c_str(), &fVal);
	if (debug) std::cout << "Entries found = " << tree->GetEntries() << std::endl;
		
	for (int i=0; i < tree->GetEntries(); ++i)
	{
			tree->GetEntry(i);
			//if (i<5)	std::cout << i << "\t" << fVal << std::endl;
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
}

//=================================================================//
//// main function
//=================================================================//
TF1* MakeQ2HepgSystHists(const std::string brname, const std::string title,
		const float xmin, 
		const float xpoint1, const float xpoint2, 
		const float xpoint3, const float xpoint4, 
		const float width1,  const float width2, 
		const float width3,  const float width4,
		//const int picbin2, const int picbin3, const int picbin4, const int picbin5,
		//int picbin2, int picbin3, int picbin4, int picbin5,
		const float cut1, const float cut2, const float cut3, const float cut4, const float cut5,
		const  int tweak, const std::string eps_name="",
		const std::string polyfit="pol1"
		)
{

	const int iNhists = 3; // base/Q2up/Q2down
	const int iNthr=5;	//20,60,100, 160,210
	const float fEth[iNthr]={cut1, cut2, cut3, cut4, cut5};
	TFile *file[iNhists][iNthr];


	file[0][0] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2BasePtHat20Results.root");
	file[0][1] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2BasePtHat60Results.root");
	file[0][2] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2BasePtHat100Results.root");
	file[0][3] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2BasePtHat160Results.root");
	file[0][4] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2BasePtHat210Results.root");

	file[1][0] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2UpPtHat20Results.root");
	file[1][1] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2UpPtHat60Results.root");
	file[1][2] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2UpPtHat100Results.root");
	file[1][3] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2UpPtHat160Results.root");
	file[1][4] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2UpPtHat210Results.root");

	file[2][0] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2DownPtHat20Results.root");
	file[2][1] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2DownPtHat60Results.root");
	file[2][2] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2DownPtHat100Results.root");
	file[2][3] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2DownPtHat160Results.root");
	file[2][4] = new TFile("/data/nbay02/b/samantha/RESULTS/10302009_Q2/11132209_ToMakeAllPlots/Q2DownPtHat210Results.root");

	for (int i=0; i < iNhists; ++i)
	{
		for (int j=0; j < iNthr; ++j)
		{
			if (file[i][j] == NULL)
			{
				std::cout << "A file "<< i << ", " << j << " is not found!returning!" << std::endl;
				return 0;
			}
		}
	}
	


	TH1F* hist[iNhists][iNthr];		//[base/Q2up/Q2down][20/60/110/160/210]
	double nEvtsAboveThr[iNhists][iNthr][2];		// number of events in sample i above the i+1  threshold
	for (int ihist=0; ihist<iNhists; ++ihist) 
	{
		for (int ith=0; ith<iNthr; ++ith) 
		{
			nEvtsAboveThr[ihist][ith][0] = 0.;		//number of evts in sample i, above the Et threshol of sample i
			nEvtsAboveThr[ihist][ith][1] = 0.;		//number of evts in sample i, above the Et threhold of sample i+1
		}
	}

	for (int i=0; i<iNhists; i++) //3 hists, base/Q2up/Q2down
	{
		std::string sType;
		if (i==0) sType = "base";
		else if (i==1) sType = "q2up";
		else if (i==2) sType = "q2down";

		for (int nth=0;nth<iNthr; ++nth)
		{
			std::stringstream name;
			name << "histPt_"<< sType << "_"<< fEth[nth];
			std::stringstream htitle;
			htitle << sType << "pt-cut (" << fEth[nth] << " GeV): " << brname ;

  			hist[i][nth] = MakeVariableBinHist (name.str(), htitle.str(), xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
			hist[i][nth]->Sumw2();
			hist[i][nth]->SetLineColor(nth+1);

			if (nth != iNthr -1)
			{
				NEvtsAbove(file[i][nth], brname,hist[i][nth], fEth[nth], fEth[nth+1], nEvtsAboveThr[i][nth][0], nEvtsAboveThr[i][nth][1]);
			} else 
			{
				NEvtsAbove(file[i][nth], brname,hist[i][nth], fEth[nth], -1, nEvtsAboveThr[i][nth][0], nEvtsAboveThr[i][nth][1]);
			}
			norm_bins(hist[i][nth]);
			//DrawClone(hist[i][nth]);

		} //end of Nthres loop
		//break;

	} //end of nhists loop
	//return;


	//calculate normalization factors
	double vScales[iNhists][iNthr-1];
	for (int ihist=0; ihist< iNhists; ++ihist)
	{ 
		for (int nth=1; nth<iNthr; ++nth) 
		{
			double dScT=0;
			if (nth==1)
			{
				dScT = (double)nEvtsAboveThr[ihist][nth-1][1]/(double) nEvtsAboveThr[ihist][nth][0];
			} else
			{
				if (nth==iNthr-1)
				{
				std::cout <<  nEvtsAboveThr[ihist][nth][0]  << "\t" << vScales[ihist][nth-2]  << "\t"
							<< nEvtsAboveThr[ihist][nth][1] << "\t" << nEvtsAboveThr[ihist][nth][0] <<std::endl;
				}
				
				dScT = ( 1. / nEvtsAboveThr[ihist][nth][0] ) * ( vScales[ihist][nth-2] * nEvtsAboveThr[ihist][nth-1][1]);
			}

			vScales[ihist][nth-1] = dScT;

		//std::cout << __LINE__ << "::Scale [" << ihist << "][" << nth-1 << "]=" <<  vScales[ihist][nth-1] << std::endl;
		}
	}

	
	for (int ihist=0; ihist< iNhists; ++ihist)
	{ 
		//new TCanvas();
		//gPad->SetLogy();
		//hist[ihist][0]->DrawClone();
		for (int nth=1; nth< iNthr; ++nth)
		{ 
			hist[ihist][nth]->Scale(vScales[ihist][nth-1]);
			//hist[ihist][nth]->DrawClone("same");
		}
	}


	int picbin2, picbin3, picbin4, picbin5;
	picbin2 = hist[0][1]->GetMaximumBin();
	picbin3 = hist[0][2]->GetMaximumBin();
	picbin4 = hist[0][3]->GetMaximumBin();
	picbin5 = hist[0][4]->GetMaximumBin();
		
	//now stich them up
	//use first hist and add rest to it.
	//first need to zero out all bins above the next threhold
	for (int ihist=0; ihist< iNhists; ++ihist)
	{ 
		for (int bin = picbin2; bin <=hist[ihist][0]->GetNbinsX(); ++bin)
		{	
			hist[ihist][0]->SetBinContent(bin, 0);
			hist[ihist][0]->SetBinError(bin, 0);
		}

		////use 60GeV sample for >Et>80
		for (int bin = picbin2; bin < picbin3; ++bin)
		{	
			//std::cout << "60GeV bin, cont = "<< bin << "\t" << hist[ihist][1]->GetBinContent(bin) << std::endl;
			hist[ihist][0]->SetBinContent(bin, hist[ihist][1]->GetBinContent(bin));
			hist[ihist][0]->SetBinError(bin, hist[ihist][1]->GetBinError(bin));
		}

		////use 110GeV sample for >Et>140
		//for (int bin = 10; bin <= hist[ihist][0]->GetNbinsX(); ++bin)
		for (int bin = picbin3; bin < picbin4; ++bin)
		{	
			//std::cout << "110GeV bin, cont = "<< bin << "\t" << hist[ihist][2]->GetBinContent(bin) << std::endl;
			hist[ihist][0]->SetBinContent(bin, hist[ihist][2]->GetBinContent(bin));
			hist[ihist][0]->SetBinError(bin, hist[ihist][2]->GetBinError(bin));
		}

		////use 160GeV sample for >Et>200
		//for (int bin = 8; bin <= hist[ihist][0]->GetNbinsX(); ++bin)
		for (int bin = picbin4; bin < picbin5; ++bin)
		{	
			//std::cout << "160GeV bin, cont = "<< bin << "\t" << hist[ihist][3]->GetBinContent(bin) << std::endl;
			hist[ihist][0]->SetBinContent(bin, hist[ihist][3]->GetBinContent(bin));
			hist[ihist][0]->SetBinError(bin, hist[ihist][3]->GetBinError(bin));
		}


		////use 210GeV sample for >Et>240
		for (int bin = picbin5; bin <=hist[ihist][0]->GetNbinsX(); ++bin)
		{	
			//std::cout << "210GeV bin, cont = "<< bin << "\t" << hist[ihist][4]->GetBinContent(bin) << std::endl;
			hist[ihist][0]->SetBinContent(bin, hist[ihist][4]->GetBinContent(bin));
			hist[ihist][0]->SetBinError(bin, hist[ihist][4]->GetBinError(bin));
		}
	}
	


	hist[1][0]->SetLineColor(kRed);
	hist[1][0]->SetMarkerColor(kRed);
	hist[1][0]->SetMarkerStyle(22);
	hist[2][0]->SetLineColor(kBlue);
	hist[2][0]->SetMarkerColor(kBlue);
	hist[2][0]->SetMarkerStyle(23);

	for (int ihist=1; ihist< iNhists; ++ihist)
	{ 
		hist[ihist][0]->Scale(hist[0][0]->Integral()/(1. * hist[ihist][0]->Integral()));
		hist[ihist][0]->Divide(hist[0][0]);
	}
	//gStyle->SetOptStat(0);
//	TCanvas *c = new TCanvas("HIST","HiST",1200,600);
//	c->Divide(2,1);

	TCanvas *c = new TCanvas("Ratio","Ratio");
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetGridx();
	gPad->SetGridy();
	
	hist[1][0]->SetTitle(title.c_str());
	hist[1][0]->SetStats(0);
	hist[2][0]->SetStats(0);
	hist[1][0]->GetYaxis()->SetRangeUser(0.07,1.03);
	hist[1][0]->DrawClone("PE");
	hist[2][0]->DrawClone("PE same");
	gPad->SetEditable(0);


	std::stringstream c1_name;
	c1_name << "Fit_" << brname;
	TCanvas *c1 = new TCanvas(c1_name.str().c_str(),c1_name.str().c_str());
  	TH1F *final = MakeVariableBinHist ("Q2_final", title, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);

	for (int bin=1; bin < final->GetNbinsX(); ++bin)
	{
		float y1 = hist[1][0]->GetBinContent(bin);
		float y2 = hist[2][0]->GetBinContent(bin);
		float dy1 = hist[1][0]->GetBinError(bin);
		float dy2 = hist[2][0]->GetBinError(bin);

		float dy = 0, e_dy=0;
		if (y1>0) // must to this do avoid adding errors to bin with no values
		{
			if ( (y1 > 1 && y2 > 1) || (y1 < 1 && y2 < 1) )
			{
				dy = (y1+y2)/2;
				dy = fabs (dy - 1);
				e_dy = 0.5 * sqrt ( pow(dy1,2) + pow(dy2,2) );
			} else 
			{
				dy = std::max(fabs(y1-1), fabs(y2-1));
				if (dy == fabs(y1-1)) e_dy = dy1;
				else e_dy = dy2;
			}
			//std::cout << y1 << "\t" << y2 << "\t" << dy << std::endl;	
			final->SetBinContent(bin,dy);
			final->SetBinError(bin,e_dy);
		}

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
	
	gStyle->SetOptFit(1);
	final->Fit(polyfit.c_str(),"","",final->GetBinCenter(xminbin),final->GetBinCenter(xmaxbin));


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
	TFile f("Q2_FitFunctions.root","UPDATE");
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
	

	for (int i=0; i<iNhists;++i) 
	{
		for (int nth=0;nth<iNthr; ++nth)
		{
		//delete hist[i][nth];
		//delete file[i][nth];
		}
	}
	//delete c1;
	//delete c;
	
	return fitfunc;

}


void MakeQ2HepgSystHists(int jets, const std::string which)
{	
	if (jets == 1)
	{
		if (which == "PhoEt")
		{
			MakeQ2HepgSystHists("pj1_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};Maximum % Q^{2} error",
					 30,150,250,300,650, 10,20,50,250,
					 30,60,100,160,210,
					 0,
					 "Q2_pj1_Et_pho.eps",
					 "pol1"
					 );
		}

		if (which == "JetEt")
		{
			MakeQ2HepgSystHists("pj1_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};Maximum % Q^{2} error"
					,10,200,250,300,600, 10,10,50,200,
					 15,70,120,180,240,
					1,
					 "Q2_pj1_Et_leadjet.eps",
					 "pol1"

					);
		}
		if (which == "InvMass")
		{
			MakeQ2HepgSystHists("pj1_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet):; Inv. Mass(#gamma,lead jet);Maximum % Q^{2} error"
					, 50,500, 600, 700, 1000, 25,50,75,300,
					50,150,250,400,500,
					1,
					 "Q2_pj1_InvM_phojet.eps",
					 "pol1"

					);
		}

/*		if (which == "NJet")
		{
    		MakeQ2HepgSystHists("pj1_Njet15",
				"Jet Multiplicity (E_{T}>15GeV)",
				1,2,3,4,15, 1,1,1,1,
				30,31,32,33,34,
				0,
					 "Q2_pj1_njet15.eps"
				);

		}
		*/

		if (which == "Ht")
		{
			
			MakeQ2HepgSystHists("pj1_Ht",
				"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);Maximum % Q^{2} error",
				0,300, 600, 800, 1200, 20,50,100,300,
					40,160,240,400,500,
				1,
					 "Q2_pj1_Ht.eps",
					 "pol1"

				);

		}


	}

		
	if (jets == 2)
	{

		if (which == "PhoEt")
		{
			MakeQ2HepgSystHists("pj2_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};Maximum % Q^{2} error",
					30,150,250,300,650, 10,20,50,250,
				30,70,130,190,250,
					0,
					 "Q2_pj2_Et_pho.eps",
					 "pol1"

					);
		}

		if (which == "Jet1Et")
		{
			MakeQ2HepgSystHists("pj2_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};Maximum % Q^{2} error",
					10,200,260,300,600, 10,20,40,200,
				10,60,100,160,220,
					1,
					 "Q2_pj2_Et_leadJet.eps",
					 "pol1"

					);
		}



		if (which == "InvMass_pj1j2")
		{
			MakeQ2HepgSystHists("pj2_InvM_pj1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet, subleading jet):; Inv. Mass(#gamma,lead jet, subleading jet);Maximum % Q^{2} error",
					100, 200, 300, 700, 1300, 10,20,50,200,
				100,200,300,450,600,
					0,
					 "Q2_pj2_InvM_pj1j2.eps",
					 "pol1"

					);
		}
		if (which == "InvMass_pj1")
		{
			MakeQ2HepgSystHists("pj2_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet); Inv. Mass(#gamma,lead jet);Maximum % Q^{2} error",
					100, 200, 300, 700, 1000, 10,20,50,200,
					100,150,240,400,500,
					0,
					 "Q2_pj2_InvM_pj1.eps",
					 "pol1"

					);
		}
		if (which == "TwoJetInvMass")
		{
			MakeQ2HepgSystHists("pj2_InvM_j1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(lead jet, subleading jet):; Inv. Mass(lead jet, subleading jet);Maximum % Q^{2} error",
				 	100, 200, 300, 700, 1200, 10,20,50,400,
				30,50,60,70,90,
					1,
					 "Q2_pj2_InvM_twoJets.eps",
					 "pol2"

					);
		}

		if (which == "Ht")
		{

			MakeQ2HepgSystHists("pj2_Ht",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);Maximum % Q^{2} error",
					0,300, 600, 700, 1200, 20,50,100,400,
					30,31,180,400,550,
					0,
					 "Q2_pj2_Ht.eps",
					 "pol2"

					);

		}


	}
 
}

void MakeQ2HepgSystHists(int i)
{
	if (i==1) MakeQ2HepgSystHists(1,"PhoEt");
	if (i==2) MakeQ2HepgSystHists(1,"JetEt");
	if (i==3) MakeQ2HepgSystHists(1,"InvMass");
	if (i==4) MakeQ2HepgSystHists(1,"Ht");
	//if (i==4) MakeQ2HepgSystHists(1, "NJet"); //need to change the branch memory address type to 'int' for this
	if (i==5) MakeQ2HepgSystHists(2,"PhoEt");
	if (i==6) MakeQ2HepgSystHists(2,"Jet1Et");
	if (i==7) MakeQ2HepgSystHists(2,"InvMass_pj1");
	if (i==8) MakeQ2HepgSystHists(2,"InvMass_pj1j2");
	if (i==9) MakeQ2HepgSystHists(2,"TwoJetInvMass");
	if (i==10) MakeQ2HepgSystHists(2,"Ht");
}

void MakeQ2HepgSystHists()
{
	for (int i=1; i<=10; ++i) MakeQ2HepgSystHists(i);
}

