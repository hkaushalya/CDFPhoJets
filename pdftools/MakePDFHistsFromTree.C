/*This is to make PDF systematic for final hists
 *A quick test using only the hepg informantion
 *Samantha K. Hewamanage, samantha@fnal.gov, Oct 10, 2009
 */

/* $Id: MakePDFHistsFromTree.C,v 1.2 2011/05/26 17:06:51 samantha Exp $
 * $Log: MakePDFHistsFromTree.C,v $
 * Revision 1.2  2011/05/26 17:06:51  samantha
 * ADDED: Njet and MET pdf error calcuations (hists).
 * MINOR: changes for some hist bin sizes etc for jobs.
 *
 * Revision 1.1  2009/11/19 22:49:09  samantha
 * Created this to make my final PDF systematics. Uses a tree of my final events
 * and produces a plot/fit for each of my final kinematic plots by adding 40 PDF
 * variants (up/down) in quadratures. Writes a fit function to a root file to be
 * read by the MakeHistLogAndRation.C when making final plots.
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


//TOrun:
//Load ~/samantha/RootTools/CommonTools_C.so
//.L MakePDFHistsFromTree.cc+
//MakePDFHistsFromTree()

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

	const int NFiles = 3;
	TChain *ch = new TChain("Ana/PhoJetsTemp/Hist/PJevtsWithPDFw/PhoJets");
	//const std::string path("~/RESULTS/10222009_PDFsystFullMCDataset/11122009_ForAllPlotsWithPDFTree/");
	const std::string path("~/RESULTS/05162010_PDFsystFullPhoMCwithMBdataset/");
	for (int inf=0; inf <= NFiles; ++inf)
	{
		std::stringstream file;
		
		//file << path << "PDFsyst" << inf << ".root";
		file << path << "PDFsyst_CTEQ6M_" << inf << ".root";
		ch->Add(file.str().c_str());
	}
	assert (ch->GetEntries()>0 && "no entries found in the tree!");
	std::cout << "Total entries =  " << ch->GetEntries() << std::endl;

	ch->SetBranchStatus("*",0);
	ch->SetBranchStatus(sBrName.c_str(),1);
	ch->SetBranchStatus("PDFWeights", 1);
	float fVal=0;
	int fVal_njet=0;
	const int Nwgts = 41;
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
				for (int iW=0; iW < Nwgts; ++iW)
				{
					//if (i==10) std::cout << i << "," << iW << "\t" << fVal  << "\t" << w[iW] << std::endl;
					hist[iW]->Fill(fVal,w[iW]);
				}
			}
	}
	std::cout << "Branch " << sBrName << " = " << ch->GetEntries()  << " [cut>" 
				<< fLoCut << "]" << std::endl;
	
	delete ch;
}


//=================================================================//
//// main function
//=================================================================//
TF1* MakePDFHistsFromTree(
		const std::string brname, const std::string title, 
		const float xmin, 
		const float xpoint1, const float xpoint2, 
		const float xpoint3, const float xpoint4, 
		const float width1,  const float width2, 
		const float width3,  const float width4,
		const  int tweak, const std::string eps_name=""
		)
{

	//if (gSystem->Load("~/samantha/RootTools/CommonTools_cc.so") == 0)
	if (gSystem->Load("./CommonTools_cc.so") == 0)
	{
		std::cout << "Required CommonTools_cc.so did not load. please check. returning!" <<std::endl;
	}
	

	const int NPDFs = 41;
	TH1F* hist[NPDFs];
	for (int iPDFs=0; iPDFs < NPDFs; ++iPDFs)
	{
		std::stringstream name, htitle;
		name << "pdfsyst_" << iPDFs;
		htitle << "PDF SYST (" << iPDFs << ") : " << title;
		
  		hist[iPDFs] = MakeVariableBinHist (name.str(), htitle.str(), xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
		hist[iPDFs]->Sumw2();
		hist[iPDFs]->SetLineColor(iPDFs+1);
	}

	FillHist(brname,hist,0);

//	TCanvas* c = new TCanvas("PDF Ratio", "PDF Ratios");
//	hist[0]->SetMaximum(1.30);
//	hist[0]->SetMinimum(0.97);
	for (int iW=0; iW < 41; ++iW)
	{
		norm_bins(hist[iW]);
		if (iW>0)
		{
			hist[iW]->Scale( hist[0]->Integral()/(1. * hist[iW]->Integral()));
			hist[iW]->Divide(hist[0]);
//			if (iW==1) hist[iW]->DrawClone("P");
//			else hist[iW]->DrawClone("P SAME");
		}
	}

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
	
	//now find max variation for each of the 20 PDF parameters (eigen vectors)
	std::vector<float> PDFval2;
	for (int bin=1; bin <=hist[0]->GetNbinsX(); ++bin) PDFval2.push_back(0);
	
	for (int i=1; i<NPDFs;i=i+2)
	{
		
		for (int bin=1; bin <=hist[0]->GetNbinsX(); ++bin)
		{
			float y1 = hist[i]->GetBinContent(bin);
			float y2 = hist[i+1]->GetBinContent(bin);
			float dy = 0;
			if (y1>0) // must to this do avoid adding errors to bin with no values
			{
				if ( (y1 > 1 && y2 > 1) || (y1 < 1 && y2 < 1) )
				{
					dy = (y1+y2)/2;
					dy = fabs (dy - 1);			
				} else 
				{
					dy = std::max(fabs(y1-1), fabs(y2-1));
				}
				if (bin ==1) 
				{
					DumpBin(hist[i], 1);
					std::cout << y1 << "\t" << y2 << "\t" << dy << std::endl;	
				}

				PDFval2.at(bin-1) += pow(dy,2);
			}
		}
	}
	
	std::cout << "" << std::endl;
	//for (std::vector<float>::iterator it = PDFVal2.begin(); it != PDFval2.end(); ++it)
	for (unsigned int bin=1; bin < PDFval2.size(); ++bin)
	{
			std::cout << bin << "->" << PDFval2.at(bin-1) <<  "\t" << sqrt(PDFval2.at(bin-1)) <<std::endl;  
	}
	
	std::stringstream c1_name;
	c1_name << "Fit_" << brname;
	TCanvas *c1 = new TCanvas(c1_name.str().c_str(),c1_name.str().c_str());
  	TH1F* final = MakeVariableBinHist ("final_hist", title, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
	final->Sumw2();
	final->SetMarkerStyle(20);
	final->SetMarkerColor(kBlue);
	final->SetLineColor(kBlue);
	final->SetMarkerSize(1);
	for (int bin=1; bin<final->GetNbinsX(); ++bin) //exclude the last bin
	{
		final->SetBinContent(bin,sqrt(PDFval2.at(bin-1)));
		final->SetBinError(bin,0);
	}


	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

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
	std::cout << "Fit Range = " << final->GetBinCenter(xminbin) << ", " << final->GetBinCenter(xmaxbin) << std::endl;
	
	final->DrawCopy();
	bool exclude_lastbin =false;
	if (brname == "pj1_Et_j1") exclude_lastbin = true;
	//final->Fit("pol3","","",final->GetBinCenter(xminbin),final->GetBinCenter(xmaxbin));
	if (exclude_lastbin) final->Fit("pol1","","",final->GetBinCenter(xminbin),final->GetBinCenter(xmaxbin-1));
	else  final->Fit("pol1","","",final->GetBinCenter(xminbin),final->GetBinCenter(xmaxbin));


	if (eps_name.length())
	{
		c1->Print(eps_name.c_str());
	}

	for (int i=0; i<NPDFs;++i) delete hist[i];

	TList *funcList = final->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1 *fitfunc = (TF1*) it->Next();
	fitfunc->SetName(brname.c_str());
	
	//check if there is a saved version already.
	//if so delete that before writing the new one.
	//otherwise there will be many copied of the same obj
	//with different cycles
	TFile f("PDFSyst_FitFunctions.root","UPDATE");
	if ( f.Get(brname.c_str()) != NULL)
	{
		std::stringstream old_tf;
		//old_tf << brname << ";*";  //don't do this if you want to see the final fit in canvas.
											//this deletes objects in disk AND memory.
		old_tf << brname << ";1";
		f.Delete(old_tf.str().c_str());
	}
	fitfunc->Write();
	c1->Write();
	f.ls();
	f.Close();

	
	return fitfunc;
}

void MakePDFHistsFromTree(int jets, const std::string which)
{	
	if (jets == 1)
	{
		if (which == "PhoEt")
		{
			MakePDFHistsFromTree("pj1_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};Total % PDF error",
					 30,150,250,400,550, 10,20,150,150,
					 0,
					 "PDFw_pj1_Et_pho.eps"
					 );
		}

		if (which == "JetEt")
		{
			MakePDFHistsFromTree("pj1_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};Total % PDF error"
					,10,200,250,300,400, 10,10,50,50,
					1,
					 "PDFw_pj1_Et_leadjet.eps"
					);
		}
		if (which == "InvMass")
		{
			MakePDFHistsFromTree("pj1_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet):; Inv. Mass(#gamma,lead jet);Total % PDF error"
					, 50,500, 600, 700, 1000, 25,50,50,50,
					1,
					 "PDFw_pj1_InvM_phojet.eps"
					);
		}

		if (which == "NJet")
		{
    		MakePDFHistsFromTree("pj1_Njet15",
				"Jet Multiplicity (E_{T}>15GeV)",
				1,2,3,4,15, 1,1,1,1,
				0,
					 "PDFw_pj1_njet15.eps"
				);

		}

		if (which == "Ht")
		{
			
			MakePDFHistsFromTree("pj1_Ht",
				"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);Total % PDF error",
				0,300, 600, 800, 1200, 20,50,50,50,
				1,
					 "PDFw_pj1_Ht.eps"
				);

		}

		if (which == "Met")
		{
    		MakePDFHistsFromTree("pj1_Met", "#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV};#slash{E}_{T} (GeV/c^{2});Total % PDF error",
    	 	 0,60,100,200, 450,10,20,50,50,
			 0,
			 "PDFw_pj1_Met.eps"
			);
		}
	}

	if (jets == 2)
	{

		if (which == "PhoEt")
		{
			MakePDFHistsFromTree("pj2_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};Total % PDF error",
					30,150,250,300,650, 10,20,50,50,
					0,
					 "PDFw_pj2_Et_pho.eps"
					);
		}

		if (which == "Jet1Et")
		{
			MakePDFHistsFromTree("pj2_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};Total % PDF error",
					0,200,260,300,600, 10,20,40,50,
					1,
					 "PDFw_pj2_Et_leadJet.eps"
					);
		}



		if (which == "InvMass_pj1j2")
		{
			MakePDFHistsFromTree("pj2_InvM_pj1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet, subleading jet):; Inv. Mass(#gamma,lead jet, subleading jet);Total % PDF error",
					100, 200, 300, 700, 1300, 10,20,50,50,
					0,
					 "PDFw_pj2_InvM_pj1j2.eps"
					);
		}
		if (which == "InvMass_pj1")
		{
			MakePDFHistsFromTree("pj2_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet); Inv. Mass(#gamma,lead jet);Total % PDF error",
					100, 200, 300, 700, 1000, 10,20,50,50,
					0,
					 "PDFw_pj2_InvM_pj1.eps"
					);
		}
		if (which == "TwoJetInvMass")
		{
			MakePDFHistsFromTree("pj2_InvM_j1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(lead jet, subleading jet):; Inv. Mass(lead jet, subleading jet);Total % PDF error",
				 	100, 200, 300, 700, 1200, 10,20,50,50,
					1,
					 "PDFw_pj2_InvM_twoJets.eps"
					);
		}

		if (which == "Ht")
		{

			MakePDFHistsFromTree("pj2_Ht",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);Total % PDF error",
					0,300, 600, 700, 1200, 20,50,50,50,
					0,
					 "PDFw_pj2_Ht.eps"
					);
		}

		if (which == "Met")
		{
    		MakePDFHistsFromTree("pj2_Met", "#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV};#slash{E}_{T} (GeV/c^{2});Total % PDF error",
    	 	 0,60,100,250, 450,10,20,50,50,
			 0,
			 "PDFw_pj2_Met.eps"
			);
		}


	}
 
}

void MakePDFHistsFromTree(int i)
{
	if (i==1) MakePDFHistsFromTree(1,"PhoEt");
	if (i==2) MakePDFHistsFromTree(1,"JetEt");
	if (i==3) MakePDFHistsFromTree(1,"InvMass");
	if (i==4) MakePDFHistsFromTree(1,"Ht");
	if (i==5) MakePDFHistsFromTree(1,"NJet"); //need to change the branch memory address type to 'int' for this// added a hack to automatically do this.
	if (i==6) MakePDFHistsFromTree(1,"Met");
	if (i==7) MakePDFHistsFromTree(2,"PhoEt");
	if (i==8) MakePDFHistsFromTree(2,"Jet1Et");
	if (i==9) MakePDFHistsFromTree(2,"InvMass_pj1");
	if (i==10) MakePDFHistsFromTree(2,"InvMass_pj1j2");
	if (i==11) MakePDFHistsFromTree(2,"TwoJetInvMass");
	if (i==12) MakePDFHistsFromTree(2,"Ht");
	if (i==13) MakePDFHistsFromTree(2,"Met");
}

void MakePDFHistsFromTree()
{
	for (int i=1; i<=13; ++i) MakePDFHistsFromTree(i);
}
