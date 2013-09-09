/* This is to make alpha_ssystematic for final hists
 * by comparing MRST98 base to CTEQ5L nominal.
 * Samantha K. Hewamanage, samantha@fnal.gov, Oct 10, 2009
 */

/* $Id: MakeAlphaSsystFromMRST98Tree.C,v 1.3 2011/05/26 17:12:52 samantha Exp $
 * $Log: MakeAlphaSsystFromMRST98Tree.C,v $
 * Revision 1.3  2011/05/26 17:12:52  samantha
 * ADDED: Njet and MET pdf error calcuations (hists).
 * MINOR: changes for some hist bin sizes etc for jobs.
 *
 * Revision 1.2  2009/12/26 02:34:00  samantha
 * MODIFIED: 1. TF1* MakeAlphaSsystFromMRST98Tree() to accept a eps file name
 * inplace of the fit function.
 * ADDED:	1. Fit is done using a custom fit function. fitFunction() is the first
 * version I tried and it does works.
 * 	2. fitFunction2() is the 2nd attempt which is more simpler. This is the
 * 	one I have used to make the most recent fits. See
 * 	Elog#https://hep.baylor.edu/hep/samantha/1517
 *
 * Revision 1.1  2009/12/01 18:37:52  samantha
 * This derives the Alpha_s uncertainty in MC sample. I compare CTEQ5L distributions
 * to the reweighted MRST98LO PDF set. I taking ratio between CTEQ5L (my nominal
 * distributions and MRST98LO base sample. I am using a linear fit for all ratio
 * plots for now. They do not seem very good for some plots.
 * See elog#https://hep.baylor.edu/hep/samantha/1488
 *
 *
 */

#include <iostream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLegend.h"
#include <sstream>
#include <ostream> 
#include <iomanip>
#include "TFile.h"
#include "TH1.h"
#include <vector>
#include "TTree.h"
#include "TSystem.h"
#include "../RootTools/CommonTools.hh"
#include <cmath>
#include "TChain.h"
#include "TF1.h"
#include "TMath.h"


// To run:
// Load ~/samantha/RootTools/CommonTools_C.so
// .L MakeAlphaSsystFromMRST98Tree.cc+
// MakeAlphaSsystFromMRST98Tree()


Double_t fitFunction(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
		//double pol1 = par[0]+par[1]*x[0]*x[0]+par[2]/x[0]; 
		//double pol1 = par[0]  + par[1] * x[0] * x[0]; 
		//double pol1 = par[0]  + par[1] * (par[2] - x[0]) * (par[2] - x[0]); 
		double pol1 = (par[0]+par[1]*sqrt(x[0]))/x[0]+par[2];
		//double pol1 = par[0]  + par[1] * (par[2] - x[0]) * (par[2] - x[0]) /TMath::(x[0]);
		//double pol2 = par[3]  + par[4] * x[0];

		//double wgt = 1 / (1 + exp((x[0]-par[4])/par[5]));
		//double wgt = 1 / (1 + exp((x[0]-par[2] * 1.5)/par[5]));
		//val = pol1 * wgt + pol2 * (1 - wgt);
		val = pol1;
	}
	return val;
}

Double_t fitFunction2(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
			return par[0]  + par[1] * x[0] /(par[2] + par[3] * x[0]); 
	}
	return val;
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
//finds the number of events above a threhold for given cut
//=================================================================//
void FillHist(const std::string sBrName, TH1F **hist , 
		const float fLoCut)
{
	assert (hist != NULL && "HIST is null");

	//const int NFiles = 6;
	const int NFiles = 3;
	TChain *ch = new TChain("Ana/PhoJetsTemp/Hist/PJevtsWithPDFw/PhoJets");
	//const std::string path("~/RESULTS/10222009_PDFsystFullMCDataset/11242009_ForAllPlotWithPDFTree_MRST98/JOB1/");
	const std::string path("~/RESULTS/05162010_PDFsystFullPhoMCwithMBdataset/");
	//for (int inf=1; inf <= NFiles; ++inf)
	for (int inf=0; inf <= NFiles; ++inf)
	{
		std::stringstream file;
		
		//file << path << "PDFsyst" << inf << ".root";
		file << path << "PDFsyst_MRST98LO_" << inf << ".root";
		ch->Add(file.str().c_str());
	}
	assert (ch->GetEntries()>0 && "no entries found in the tree!");
	std::cout << "Total entries =  " << ch->GetEntries() << std::endl;

	ch->SetBranchStatus("*",0);
	ch->SetBranchStatus(sBrName.c_str(),1);
	ch->SetBranchStatus("PDFWeights", 1);
	float fVal=0;
	int fVal_njet=0;
	const int Nwgts = 5; //
	float w[Nwgts];

	//this hack is to switch from float to int for njet15. 	//I had defined njet15 to be Int_t
	if (sBrName == "pj1_Njet15") ch->SetBranchAddress(sBrName.c_str(), &fVal_njet);
	else ch->SetBranchAddress(sBrName.c_str(), &fVal);
	ch->SetBranchAddress("PDFWeights", w);
	for (int i=0; i < ch->GetEntries(); ++i)
	{
			ch->GetEntry(i);
			if (sBrName == "pj1_Njet15") fVal = (float) fVal_njet;
			if (fVal>fLoCut)
			{
					hist[0]->Fill(fVal);
					hist[1]->Fill(fVal,w[0]);
			}
	}
	std::cout << "Branch " << sBrName << " = " << ch->GetEntries()  << " [cut>" 
				<< fLoCut << "]" << std::endl;
	
	delete ch;
}


//=================================================================//
//// main function
//=================================================================//
TF1* MakeAlphaSsystFromMRST98Tree(
		const std::string brname, const std::string title, 
		const float xmin, 
		const float xpoint1, const float xpoint2, 
		const float xpoint3, const float xpoint4, 
		const float width1,  const float width2, 
		const float width3,  const float width4,
		const  int tweak, const std::string eps_name="",
		const std::string sFitFunc="pol1"
		)
{

	if (gSystem->Load("~/samantha/RootTools/CommonTools_cc.so") == 0)
	{
		std::cout << "Required CommonTools_cc.so did not load. please check. returning!" <<std::endl;
	}
	
	const int NPDFs = 2;	//0 = nominal, 1=MRST98LO base
	
	TH1F* hist[NPDFs];

	for (int iPDFs=0; iPDFs < NPDFs; ++iPDFs)
	{
		std::stringstream name, htitle;
		name << "TH1F_" << brname << "_" << iPDFs;
		htitle << "#alpha _{s} systematic : " << title;
		
  		hist[iPDFs] = MakeVariableBinHist (name.str(), htitle.str(), xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
		hist[iPDFs]->Sumw2();
		hist[iPDFs]->SetLineColor(iPDFs+1);
	}

	FillHist(brname,hist,0);

	norm_bins(hist[0]);
	norm_bins(hist[1]);
	hist[1]->Scale(hist[0]->Integral()/(1. * hist[1]->Integral()));
	hist[1]->Divide(hist[0]);

	//need to find first non-zero bin so I can assign the final error for 
	//the correct bin
	int iFirstBin = 1;
	for (int bin=1; bin <=hist[0]->GetNbinsX(); ++bin) 
	{
		if (hist[0]->GetBinContent(bin))
		{
			iFirstBin = bin;
			break;
		}
	}

	for (int bin=1; bin <=hist[0]->GetNbinsX(); ++bin)
	{
		float y1 = hist[1]->GetBinContent(bin);
		float dy = hist[1]->GetBinError(bin);
		if (y1>0) // must to this do avoid adding errors to bin with no values
		{
			y1 = fabs (y1 - 1);			
			hist[1]->SetBinContent(bin, y1);
			hist[1]->SetBinError(bin, dy);
		}
	}

	std::stringstream c1_name;
	c1_name << "Canvas_" << brname;
	TCanvas *c1 = new TCanvas(c1_name.str().c_str(),c1_name.str().c_str());
	gPad->SetGridx();
	gPad->SetGridy();

	std::stringstream finalname;
	finalname << "TH1F_" << brname;
  	TH1F* final = dynamic_cast<TH1F*> (hist[1]->Clone(finalname.str().c_str()));


	/*if (tweak)
	{
		int firstbin = 0;
		float firstbin_val, firstbin_err;
		for (int bin=1; bin <=final->GetNbinsX(); ++bin)
		{
			if (final->GetBinContent(bin))
			{
				firstbin = bin;
				firstbin_val = final->GetBinContent(bin);
				firstbin_err = final->GetBinError(bin);
				break;
			}
		}
	//	float binval[5],binerr[5];
		for (int bin=firstbin+1; bin<firstbin+5; ++bin)
		{
	//		binval[bin-firstbin] = final->GetBinContent(bin);
	//		binerr[bin-firstbin] = final->GetBinError(bin);
			if (final->GetBinContent(bin)<firstbin_val)
			{
				final->SetBinContent(bin, firstbin_val);
				final->SetBinError(bin, firstbin_err);
			}
		}
	}
	*/



	
	final->SetTitle(title.c_str());
	final->Sumw2();
	final->SetMarkerStyle(20);
	final->SetMarkerColor(kBlue);
	final->SetLineColor(kBlue);
	final->SetMarkerSize(1);
	final->SetLineWidth(2);

	int maxbin = final->GetMaximumBin();
	float ymax = final->GetBinContent(maxbin);
	
	final->SetMinimum(-0.1);
	final->SetMaximum(ymax+0.2);
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	final->Print();
	final->Draw();
	//find fit range
	int xminbin=0, xmaxbin=0;
	FindNonZeroXrange(final,xminbin,xmaxbin);
	std::cout << "minx, maxx = " << xminbin << "\t" << xmaxbin << std::endl;
	
	std::cout << "sFitFunc = " << sFitFunc << std::endl;

	//TF1 *f4 = new TF1("fit_func4",fitFunction,final->GetBinCenter(xminbin+3),final->GetXaxis()->GetBinUpEdge(xmaxbin),5);
	//TF1 *f4 = new TF1("fit_func4",fitFunction,final->GetBinCenter(xminbin),final->GetXaxis()->GetBinUpEdge(xmaxbin),3);
	TF1 *f4 = new TF1("fit_func4",fitFunction,final->GetBinCenter(xminbin),final->GetXaxis()->GetBinCenter(xmaxbin),3);
	f4->SetParameter(0,0.5);
	f4->SetParameter(1,0.5);
	//f4->SetParLimits(1,20,40);
	f4->SetParameter(2,30);
//	f4->SetParameter(3,1.5);
//	f4->SetParameter(4,30);
	//f4->SetParLimits(4,20,35);
	//f4->SetParameter(5,10.5); //for photon Et
	//f4->SetParLimits(5,10.5,10.5); //for photon Et

//	f4->SetParameter(5,100);
//	f4->SetParLimits(5,70,135);
	//f4->SetParameter(5,10.5); //for Ht
	//f4->SetParLimits(5,10.5,10.5);
	//f4->SetParameter(5,100); //for

	f4->SetLineColor(8);

	//final->Fit(f4,"R+");
	final->Fit("pol1","","",final->GetBinCenter(xminbin),final->GetXaxis()->GetBinCenter(xmaxbin));
	
	if (eps_name.length())
	{
		std::cout << "eps_name = " << eps_name << std::endl;
		c1->Print(eps_name.c_str());
	}


	TList *funcList = final->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1 *fitfunc = (TF1*) it->Next();
	if (fitfunc == NULL)
	{
		std::cout << "Fit function not found.! " << std::endl;
		return 0;
	}

	//std::stringstream fitfunc_name;
	//fitfunc_name << "TF1_" << brname;
	fitfunc->SetName(brname.c_str());
	assert(fitfunc != NULL && "Final fit function is null!");
	
	
	//check if there is a saved version already.
	//if so delete that before writing the new one.
	//otherwise there will be many copied of the same obj
	//with different cycles
	TFile f("AlphaS_Syst_FitFunctions.root","UPDATE");
	assert (! f.IsZombie() && "Cannot open file to write final fit info");
	if ( f.Get(brname.c_str()) != NULL)
	{
		std::stringstream old_tf;
		//old_tf << brname << ";*";  //don't do this if you want to see the final fit in canvas.
											//this deletes objects in disk AND memory.
		old_tf << brname << ";1";
		f.Delete(old_tf.str().c_str());
		std::stringstream c1_tmp;
		c1_tmp << c1_name << ";1";
		f.Delete(c1_tmp.str().c_str());

		std::stringstream old_hname;
		old_hname << final->GetName() << ";1";
		f.Delete(old_hname.str().c_str());
	}
	fitfunc->Write();
	c1->Write();
	final->Write();
//	f.ls();
	f.Close();
	
	for (int i=0; i<NPDFs;++i) delete hist[i];
	
	return fitfunc;
}

void MakeAlphaSsystFromMRST98Tree(int jets, const std::string which)
{	
	if (jets == 1)
	{
		if (which == "PhoEt")
		{
			MakeAlphaSsystFromMRST98Tree("pj1_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};% #alpha _{s} error  (MRST98LO/CTEQ5L)",
					 30,150,250,300,650, 10,20,50,250,
					 0,
					 "PDFw_pj1_Et_pho_pol1.eps",
					 "pol2"
					 );
		}

		if (which == "JetEt")
		{
			MakeAlphaSsystFromMRST98Tree("pj1_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};% #alpha _{s} error  (MRST98LO/CTEQ5L)"
					,10,200,250,300,600, 10,10,50,200,
					1,
					 "PDFw_pj1_Et_leadjet_pol1.eps",
					 "pol2"

					);
		}
		if (which == "InvMass")
		{
			MakeAlphaSsystFromMRST98Tree("pj1_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet):; Inv. Mass(#gamma,lead jet);% #alpha _{s} error  (MRST98LO/CTEQ5L)"
					, 50,500, 600, 700, 1000, 25,50,75,300,
					1,
					 "PDFw_pj1_InvM_phojet_pol1.eps"
					);
		}

		if (which == "NJet")
		{
    		MakeAlphaSsystFromMRST98Tree("pj1_Njet15",
				"Jet Multiplicity (E_{T}>15GeV)",
				1,2,3,4,15, 1,1,1,1,
				0,
					 "PDFw_pj1_njet15_pol1.eps"
				);

		}

		if (which == "Ht")
		{
			
			MakeAlphaSsystFromMRST98Tree("pj1_Ht",
				"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);% #alpha _{s} error  (MRST98LO/CTEQ5L)",
				0,300, 600, 800, 1200, 20,50,100,300,
				1,
					 "PDFw_pj1_Ht_pol1.eps",
					 "pol2"

				);

		}


		if (which == "Met")
		{
    		MakeAlphaSsystFromMRST98Tree("pj1_Met", "#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV};#slash{E}_{T} (GeV/c^{2});% #alpha _{s} error  (MRST98LO/CTEQ5L)",
    	 	 0,60,100,200, 450,10,20,50,50,
			 0,
			 "PDFw_pj1_Met_pol1.eps"
			);
		}


	}

	if (jets == 2)
	{

		if (which == "PhoEt")
		{
			MakeAlphaSsystFromMRST98Tree("pj2_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};% #alpha _{s} error  (MRST98LO/CTEQ5L)",
					30,150,250,300,650, 10,20,50,250,
					0,
					 "PDFw_pj2_Et_pho_pol1.eps",
					 "pol2"

					);
		}

		if (which == "Jet1Et")
		{
			MakeAlphaSsystFromMRST98Tree("pj2_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};% #alpha _{s} error  (MRST98LO/CTEQ5L)",
					0,200,260,300,600, 10,20,40,200,
					1,
					 "PDFw_pj2_Et_leadJet_pol1.eps",
					 "pol2"

					);
		}



		if (which == "InvMass_pj1j2")
		{
			MakeAlphaSsystFromMRST98Tree("pj2_InvM_pj1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet, subleading jet):; Inv. Mass(#gamma,lead jet, subleading jet);% #alpha _{s} error  (MRST98LO/CTEQ5L)",
					100, 200, 300, 700, 1300, 10,20,50,200,
					0,
					 "PDFw_pj2_InvM_pj1j2_pol1.eps"
					);
		}
		if (which == "InvMass_pj1")
		{
			MakeAlphaSsystFromMRST98Tree("pj2_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet); Inv. Mass(#gamma,lead jet);% #alpha _{s} error  (MRST98LO/CTEQ5L)",
					100, 200, 300, 700, 1000, 10,20,50,200,
					0,
					 "PDFw_pj2_InvM_pj1_pol1.eps"
					);
		}
		if (which == "TwoJetInvMass")
		{
			MakeAlphaSsystFromMRST98Tree("pj2_InvM_j1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(lead jet, subleading jet):; Inv. Mass(lead jet, subleading jet);% #alpha _{s} error  (MRST98LO/CTEQ5L)",
				 	100, 200, 300, 700, 1200, 10,20,50,400,
					1,
					 "PDFw_pj2_InvM_twoJets_pol1.eps"
					);
		}

		if (which == "Ht")
		{

			MakeAlphaSsystFromMRST98Tree("pj2_Ht",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);% #alpha _{s} error  (MRST98LO/CTEQ5L)",
					0,300, 600, 700, 1200, 20,50,100,400,
					0,
					 "PDFw_pj2_Ht_pol1.eps",
					 "pol2"

					);

		}

		if (which == "Met")
		{
    		MakeAlphaSsystFromMRST98Tree("pj2_Met", "#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV};#slash{E}_{T} (GeV/c^{2});% #alpha _{s} error  (MRST98LO/CTEQ5L)",
    	 	 0,60,100,250, 450,10,20,50,50,
			 0,
			 "PDFw_pj2_Met_pol1.eps"
			);
		}



	}
 
}

void MakeAlphaSsystFromMRST98Tree(int i)
{
	if (i==1) MakeAlphaSsystFromMRST98Tree(1,"PhoEt");
	if (i==2) MakeAlphaSsystFromMRST98Tree(1,"JetEt");
	if (i==3) MakeAlphaSsystFromMRST98Tree(1,"InvMass");
	if (i==4) MakeAlphaSsystFromMRST98Tree(1,"Ht");
	if (i==5) MakeAlphaSsystFromMRST98Tree(1, "NJet"); //need to change the branch memory address type to 'int' for this// added a hack to automatically do this.
	if (i==6) MakeAlphaSsystFromMRST98Tree(1,"Met");
	if (i==7) MakeAlphaSsystFromMRST98Tree(2,"PhoEt");
	if (i==8) MakeAlphaSsystFromMRST98Tree(2,"Jet1Et");
	if (i==9) MakeAlphaSsystFromMRST98Tree(2,"InvMass_pj1");
	if (i==10) MakeAlphaSsystFromMRST98Tree(2,"InvMass_pj1j2");
	if (i==11) MakeAlphaSsystFromMRST98Tree(2,"TwoJetInvMass");
	if (i==12) MakeAlphaSsystFromMRST98Tree(2,"Ht");
	if (i==13) MakeAlphaSsystFromMRST98Tree(2,"Met");
}

void MakeAlphaSsystFromMRST98Tree()
{
	for (int i=1; i<=13; ++i) MakeAlphaSsystFromMRST98Tree(i);
}
