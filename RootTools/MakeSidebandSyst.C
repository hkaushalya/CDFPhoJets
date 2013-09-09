/*This is to make the systematic unceratinty for my photons sideband
 * representation of the fake tight photons (jets) in the signal sample
 * This methods derived a difference between two plots using a weighted
 * events and comparing them to my nominal sideband distributions.
 * The weights are derived from the di-jet MC compareing jets passing
 * tight photon id cuts and jets passing loose photon id cuts.
 * see elog#https://hep.baylor.edu/hep/samantha/1459
 * Above fit for photon Et is used to rewegiht all the sideband events 
 * to make a variation. This systematic seems to give a larger unceratinty
 * compared to ID cut vary method 
 * (see https://hep04.baylor.edu/hep/samantha/1298)
 * Samantha K. Hewamanage, samantha@fnal.gov, 11-24-2009
 */

/* $Id: MakeSidebandSyst.C,v 1.1 2013/02/28 03:27:20 samantha Exp $
 * $Log: MakeSidebandSyst.C,v $
 * Revision 1.1  2013/02/28 03:27:20  samantha
 * Final commit. no checks. These were never commited to cvs.
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
Double_t fitFunction(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
		val = par[0]+par[1]*sqrt(x[0])+par[2]*x[0];
		//val = (par[0]+par[1]*sqrt(x[0]))+ par[2]/sqrt(x[0]);
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
//// main function
//=================================================================//
TF1* MakeSidebandSyst(const std::string brname, const std::string title,
		const std::string path,
		const float xmin, 
		const float xpoint1, const float xpoint2, 
		const float xpoint3, const float xpoint4, 
		const float width1,  const float width2, 
		const float width3,  const float width4,
		const  int tweak, const std::string eps_name="",
		const std::string polyfit="pol1"
		)
{

	TFile *file_base = new TFile ("~/ICHEP_LaptopVersion/PhoJets_data.root");
	TFile *file_vary = new TFile ("~/RESULTS/11232009_SidebandSyst/SidebandSystFromDiJets.root");
	assert ((! file_base->IsZombie()) && "base file not found!");
	assert ((! file_vary->IsZombie()) && "vary file not found!");

	TH1* hist_base = (TH1*) file_base->Get(path.c_str());
	TH1* hist_vary = (TH1*) file_vary->Get(path.c_str());
	
	assert (hist_base != NULL && "base hist not found"); 
	assert (hist_vary != NULL && "vary hist not found"); 
	gStyle->SetOptStat(111111);
	/*new TCanvas();
	gPad->SetLogy();
	hist_base->DrawClone();
	new TCanvas();
	gPad->SetLogy();
	hist_vary->DrawClone();
*/

	hist_base->Print();
	hist_vary->Print();
	//rebin to ti match final hist format
  	hist_base = (TH1*) MakeVariableBinHist (hist_base,xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
  	hist_vary = (TH1*) MakeVariableBinHist (hist_vary, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
	
	hist_base->SetDirectory(0);
	hist_vary->SetDirectory(0);

	//number of entried in both base and vary must be the same.
	// integral's are different.
	hist_vary->Scale(hist_base->Integral()/(1. * hist_vary->Integral()));
	hist_vary->Divide(hist_base);
	
//	new TCanvas();
//	hist_vary->DrawClone();


	std::stringstream c1_name;
	c1_name << "Fit_" << brname;
	TCanvas *c1 = new TCanvas(c1_name.str().c_str(),c1_name.str().c_str(),1000,800);
	c1->Divide(1,2);
	c1->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetGridx();
	gPad->SetGridy();
	std::stringstream finalhist_name;
	finalhist_name << "TH1F_"<< brname;
  	TH1F *final = MakeVariableBinHist (finalhist_name.str().c_str(), title, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
  	TH1F *residual = MakeVariableBinHist ("Residual_SidebandSyst_final", "Residual Plot;;Fit -- Val / Val Error", xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);

	for (int bin=1; bin < final->GetNbinsX(); ++bin)
	{
		float y1 = hist_vary->GetBinContent(bin);
		float err = hist_vary->GetBinError(bin);

		if (y1>0) // must to this do avoid adding errors to bin with no values
		{
			float dy = fabs (y1 - 1);
			final->SetBinContent(bin,dy);
			final->SetBinError(bin,err);
		}

	}
	
	final->SetLineWidth(2);
	final->SetLineColor(kBlue);
	final->SetMarkerStyle(20);
	final->SetMarkerSize(1.2);
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
	
	TF1 *f4 = new TF1("fit_func4",fitFunction,final->GetBinCenter(xminbin),final->GetXaxis()->GetBinCenter(xmaxbin),3);
	f4->SetParameter(0,0.5);
	f4->SetParameter(1,0.5);
	f4->SetParameter(2,30);
	f4->SetLineColor(8);

	gStyle->SetOptFit(1);
	//final->Fit(polyfit.c_str(),"","",final->GetBinCenter(xminbin),final->GetBinCenter(xmaxbin));
	final->Fit(f4,"R+");


	TList *funcList = final->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1* fitfunc = (TF1*) it->Next();
	fitfunc->SetName(brname.c_str());
	//check if there is a saved version already.
	//if so delete that before writing the new one.
	//otherwise there will be many copied of the same obj
	//with different cycles
	TFile f("SidebandSyst_FitFunctions.root","UPDATE");
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

	if ( f.Get(final->GetName()) != NULL)
	{
		std::stringstream old_th1f;
		old_th1f << final->GetName() << ";1";
		f.Delete(old_th1f.str().c_str());
	}
	final->Write();

	f.ls();
	f.Close();
	
	for (int bin=1; bin < residual->GetNbinsX(); ++bin)
	{
		float val = final->GetBinContent(bin);
		float err = final->GetBinError(bin);
		float fit = fitfunc->Eval(residual->GetBinCenter(bin));
		if (err>0) residual->SetBinContent(bin,(fit-val)/err);
		residual->SetBinError(bin,0);
	}

	c1->cd(2);
	gPad->SetGridx();
	gPad->SetGridy();
	residual->SetMarkerStyle(20);
	residual->SetMarkerSize(1);
	
	residual->Draw();
	c1->cd();
	
	if (eps_name.length())
	{
		c1->Print(eps_name.c_str());
		std::stringstream rationame;
		int dot = eps_name.find(".");
		std::string subs= eps_name.substr(0,dot);
		rationame << subs << "_ratio.eps";
		//c1->Print(rationame.str().c_str());
	}

	
	return fitfunc;

}


void MakeSidebandSyst(int jets, const std::string which)
{	
	if (jets == 1)
	{
		if (which == "PhoEt")
		{
			MakeSidebandSyst("pj1_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/1Jet/Photon/EtCorr",
					 30,150,250,300,650, 10,20,50,250,
					 0,
					 "Sideband_pj1_Et_pho.eps",
					 "pol2"
					 );
		}

		if (which == "JetEt")
		{
			MakeSidebandSyst("pj1_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/1Jet/LeadJet/EtCorr",
					10,200,250,300,600, 10,10,50,200,
					1,
					 "Sideband_pj1_Et_leadjet.eps",
					 "pol2"

					);
		}
		if (which == "InvMass")
		{
			MakeSidebandSyst("pj1_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet):; Inv. Mass(#gamma,lead jet);Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/1Jet/PhotonLeadJet/InvMass",
					 50,500, 600, 700, 1000, 25,50,75,300,
					1,
					 "Sideband_pj1_InvM_phojet.eps",
					 "pol2"
					);
		}

		if (which == "NJet")
		{
    		MakeSidebandSyst("pj1_Njet15", "Jet Multiplicity (E_{T}>15GeV)",
					"Hist/SIDEBAND/1Jet/Event/NJet15",
				1,2,3,4,15, 1,1,1,1,
				0,
					 "Sideband_pj1_njet15.eps",
					 "pol1"
				);

		}
		

		if (which == "Ht")
		{
			
			MakeSidebandSyst("pj1_Ht",
				"#gamma ^{E_{T}>30GeV} + #geq 1 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/1Jet/Event/Ht",
				0,300, 600, 800, 1200, 20,50,100,300,
				1,
					 "Sideband_pj1_Ht.eps",
					 "pol2"

				);

		}


	}

		
	if (jets == 2)
	{

		if (which == "PhoEt")
		{
			MakeSidebandSyst("pj2_Et_pho",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{#gamma};E_{T}^{#gamma};Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/2Jet/Photon/EtCorr",
					30,150,250,300,650, 10,20,50,250,
					0,
					 "Sideband_pj2_Et_pho.eps",
					 "pol2"

					);
		}

		if (which == "Jet1Et")
		{
			MakeSidebandSyst("pj2_Et_j1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: E_{T}^{Lead jet}:;E_{T}^{lead jet};Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/2Jet/LeadJet/EtCorr",
					10,200,260,300,600, 10,20,40,200,
					1,
					 "Sideband_pj2_Et_leadJet.eps",
					 "pol2"

					);
		}



		if (which == "InvMass_pj1j2")
		{
			MakeSidebandSyst("pj2_InvM_pj1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet, subleading jet):; Inv. Mass(#gamma,lead jet, subleading jet);Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/2Jet/Photon2Jets/InvMass",
					100, 200, 300, 700, 1300, 10,20,50,200,
					0,
					 "Sideband_pj2_InvM_pj1j2.eps",
					 "pol1"

					);
		}
		if (which == "InvMass_pj1")
		{
			MakeSidebandSyst("pj2_InvM_pj1",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(#gamma, lead jet); Inv. Mass(#gamma,lead jet);Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/2Jet/PhotonLeadJet/InvMass",
					100, 200, 300, 700, 1000, 10,20,50,200,
					0,
					 "Sideband_pj2_InvM_pj1.eps",
					 "pol2"

					);
		}
		if (which == "TwoJetInvMass")
		{
			MakeSidebandSyst("pj2_InvM_j1j2",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: Inv. Mass(lead jet, subleading jet):; Inv. Mass(lead jet, subleading jet);Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/2Jet/2Jets/InvMass",
				 	100, 200, 300, 700, 1200, 10,20,50,400,
					1,
					 "Sideband_pj2_InvM_twoJets.eps",
					 "pol1"

					);
		}

		if (which == "Ht")
		{

			MakeSidebandSyst("pj2_Ht",
					"#gamma ^{E_{T}>30GeV} + #geq 2 Jets^{E_{T}#geq 15GeV}: H_{T} (#gamma+all jets E_{T}>15GeV);H_{T} (Gev);Maximum %  Error for #gamma sideband error",
					"Hist/SIDEBAND/2Jet/Event/Ht",
					0,300, 600, 700, 1200, 20,50,100,400,
					0,
					 "Sideband_pj2_Ht.eps",
					 "pol2"

					);

		}


	}
 
}

void MakeSidebandSyst(int i)
{
	if (i==1) MakeSidebandSyst(1,"PhoEt");
	if (i==2) MakeSidebandSyst(1,"JetEt");
	if (i==3) MakeSidebandSyst(1,"InvMass");
	if (i==4) MakeSidebandSyst(1,"Ht");
	//if (i==4) MakeSidebandSyst(1, "NJet"); //need to change the branch memory address type to 'int' for this
	if (i==5) MakeSidebandSyst(2,"PhoEt");
	if (i==6) MakeSidebandSyst(2,"Jet1Et");
	if (i==7) MakeSidebandSyst(2,"InvMass_pj1");
	if (i==8) MakeSidebandSyst(2,"InvMass_pj1j2");
	if (i==9) MakeSidebandSyst(2,"TwoJetInvMass");
	if (i==10) MakeSidebandSyst(2,"Ht");
}

void MakeSidebandSyst()
{
	for (int i=1; i<=10; ++i) MakeSidebandSyst(i);
}

