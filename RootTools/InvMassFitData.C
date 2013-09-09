#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPad.h"
#include "TH1F.h"
#include "TMath.h"
#include "TStyle.h"
#include <math.h>
#include "TSystem.h"
//#include "/cdf/home/samantha/samantha/RootTools/CommonTools.hh"
#include "/cdf/home/samantha/samantha/RootTools/IOColors.hh"
#include "TMinuit.h"
#include "TLegend.h"
#include <vector>
#include <algorithm>


using namespace std;

void MakeFreqPlot(const TH1* hist, const int plot)
{
	assert(hist != NULL && "MakeFreqPlot:: hist is null");

	TH1F* h_frq = new TH1F("FREQ",";(Observed - Expected) / #sqrt{Expected};FREQUENCY" , 60, -1.5,1.5);
	//h_frq->SetTitle(hist->GetTitle());
	h_frq->GetXaxis()->CenterTitle(1);
	h_frq->SetStats(0);
	for (int bin=1; bin<=hist->GetNbinsX(); ++bin)
	{
		h_frq->Fill(hist->GetBinContent(bin));
	}
	
	//TCanvas *c2 = new TCanvas("c2","FREQ HIST",700,500);
	TPad* pad_guass = new TPad("gaus_pad","Freq Dist",0.7,0.65,0.99,0.87);
	pad_guass->Draw();
	pad_guass->cd();
	h_frq->SetLineColor(kBlue);
	h_frq->Draw();

	/*if (plot == 1) gPad->Print("pj1_pj1Mass_Freq.eps");
	else if (plot == 2) gPad->Print("pj2_pj1Mass_Freq.eps");
	else if (plot == 3) gPad->Print("pj2_pj1j2Mass_Freq.eps");
	else if (plot == 4) gPad->Print("pj2_j1j2Mass_Freq.eps");
   */
}


void DumpHist(const TH1* hist)
{
	assert(hist != NULL && "DumpHist:: hist is null");
	for (int bin=1; bin<=hist->GetNbinsX(); ++bin)
	{
		std::cout << "bin [" << bin << "], loEg["<< hist->GetBinLowEdge(bin) 
			<< " ] ,cont,err ["  << hist->GetBinContent(bin) 
			<< ", " << hist->GetBinError(bin) << std::endl;
	}

}

int FindLastBin(TH1* hist)
{
	assert (hist != NULL && "FindLastBin:: hist null");
	for (int bin = hist->GetNbinsX(); bin >=1 ; --bin)
	{
		if (hist->GetBinContent(bin)>0) return bin;
	}

	std::cout << "No non zero bins found! returning -1" << std::endl;
	return -1;
}
//------ Gauss mean and Landau MPV: A0+A1*Pt+A3/Pt ---- default
double fitFunc0(double *x, double *par)
{
  double fitval=0.0;
  if(x[0]>0.0) fitval=par[0]+par[1]*x[0]+par[2]/x[0];
  else fitval=-1.0E6;
  return fitval;
}
//-------------- fit function: [Const*Gaus(-x/(1+x))+Landau(-x/(1+x))]/(1+Const)
double myFunc1(double *x, double *par)
{
	double fitval=0.0;
	double arg=-x[0]/(x[0]+1.0);
	double arg_L=-x[0]/(x[0]+1.0);
	double mean=par[1];
	double sigmaG=par[2];
	double mpv=par[3];
	double sigmaL=par[4];
	double normL=par[5];
	if(normL<0.0) return -1.0E6;
	double f1=normL/(1.0+normL)*TMath::Gaus(arg,mean,sigmaG);
	double f2=TMath::Landau(arg_L,mpv,sigmaL)/(1.0+normL);
	fitval=par[0]*(f1+f2);
	return fitval;
}

void GetError(TH1* h, const TF1* f, const bool de=false)
{

	for (int bin=1; bin <= h->GetNbinsX(); ++bin)
	{
		if (h->GetBinContent(bin)>0)
		{
			if (de)
			{
			//	std::cout << green << bin << " = " 
			//		<< (h->GetBinContent(bin)  - f->Eval(h->GetBinCenter(bin))) / h->GetBinError(bin)
			//		<< clearatt << std::endl;
			} else 
			{
			//	std::cout << bin << " = " 
			//		<< (h->GetBinContent(bin)  - f->Eval(h->GetBinCenter(bin)))  - 1.0 
			//		<< std::endl;

			}

			if (de)
			{
				if (h->GetBinError(bin)>0) h->SetBinContent(bin, (h->GetBinContent(bin)  - f->Eval(h->GetBinCenter(bin))) / h->GetBinError(bin));
			}
			else 
			{
				if (f->Eval(h->GetBinCenter(bin))>0)h->SetBinContent(bin, (h->GetBinContent(bin) / f->Eval(h->GetBinCenter(bin))) - 1.);
			}
			h->SetBinError(bin,0);
		}

	}
	//if (de) h->GetYaxis()->SetTitle("(DATA - Fit)/Data Stat Error");
	//else h->GetYaxis()->SetTitle("DATA/Fit - 1");
	h->GetYaxis()->SetTitle("");
	
	h->SetMarkerStyle(20);
	h->SetMarkerSize(1.0);
	
	
}

void GetObsExpt(TH1* h, const TF1* f)
{

	std::cout << blue << "bin[loEdge]\t Obsrv(O)\t Exptd(E) \t O-E/sqrt(E) "
			 << clearatt << std::endl;
	for (int bin=1; bin <= h->GetNbinsX(); ++bin)
	{

		const double dObs = h->GetBinContent(bin);
		const double dExp = f->Eval(h->GetBinCenter(bin));
		double rat = 0;
		if (dExp>0) rat = (dObs - dExp)/sqrt(dExp);
		
	//	 std::cout << blue << bin  << "[" << h->GetBinLowEdge(bin) << "] = " 
	//		 << dObs << "\t" << dExp <<"\t" << rat
	//		 << clearatt << std::endl;

		h->SetBinContent(bin,rat);
		h->SetBinError(bin,0);

	}
	h->GetYaxis()->SetTitle("");
	
	h->SetMarkerStyle(22);
	h->SetMarkerSize(1.0);
	h->SetMaximum(5.);
	h->SetMinimum(-5.);
}

//=================================================================//
// divide bins by bin width
//=================================================================//
void Local_NormalizeBinWidth(TH1* hist)
{
	assert (hist != NULL && "requirement failed"); //spec
	assert (hist->GetDimension() == 1 && "CommonTools::NormalizeBinWidth:: hist is not 1-D");

	for (unsigned bin = 1; bin <= unsigned (hist->GetNbinsX()); ++ bin)
	{
		const float width = hist->GetXaxis()->GetBinWidth (bin);
		hist->SetBinContent (bin, hist->GetBinContent (bin) / width);
		hist->SetBinError (bin, hist->GetBinError (bin) / width);
	};

}



std::auto_ptr<TH1> Local_FillVarBinHist (const std::vector<float>& bins, TH1 *input)
{
	assert (input != NULL && "requirement failed"); //spec
	assert (input->GetDimension() == 1 && "requirement failed"); //spec
	const unsigned nbin = unsigned (input->GetNbinsX ());

	std::auto_ptr<TH1> result (new TH1F (input->GetName(), input->GetTitle(),
						 bins.size() - 1, &bins[0]));

	result->SetBinContent (0, input->GetBinContent (0));
	result->SetBinError (0, input->GetBinError (0));
	result->SetBinContent (bins.size(), input->GetBinContent (nbin));
	result->SetBinError (bins.size(), input->GetBinError (nbin));
	for (unsigned bin = 1; bin != nbin; ++ bin)
	{
		const float low = input->GetXaxis()->GetBinLowEdge (bin);
		const float high = input->GetXaxis()->GetBinUpEdge (bin);
		const unsigned bin1 = result->FindBin (0.99 * low + 0.01 * high);
		const unsigned bin2 = result->FindBin (0.01 * low + 0.99 * high);
		assert (bin1 == bin2 && "requirement failed: histogram bin edges don't match"); //spec

		const float va = input->GetBinContent (bin);
		const float ea = input->GetBinError (bin);
		const float vb = result->GetBinContent (bin2);
		const float eb = result->GetBinError (bin2);
		const float vc = va + vb;
		const float ec = sqrt (ea * ea + eb * eb);
		result->SetBinContent (bin2, vc);
		result->SetBinError (bin2, ec);
	};

	assert (result.get() != NULL && "postcondition failed"); //spec
	return result;
}

//-----------------------------------------------------------------------------
// effects: make a new binning with variable bin sizes, and the given
//   cutoff points between bin sizes
// returns: the binning (by bin edges)
// guarantee: strong
//// throws: std::bad_alloc
//-----------------------------------------------------------------------------
std::vector<float> 
Local_GetVarBinVector (float down, float up5, float up10, float up20, float up50, 
					float width1, float width2, float width3, float width4)
{
	std::vector<float> result;
	float point = down;
	const unsigned nregion = 4;
	const float up [] = {up5, up10, up20, up50};
	const float step [] = {width1, width2, width3, width4};		//1j pet

	result.push_back (point);
	for (unsigned region = 0; region != nregion; ++ region)
	{
		while (point + step[region] < up[region] + 1)
		{
			point += step[region];
			result.push_back (point);
	 	};
	};
	
	if (point + 1 < up[nregion-1])
		result.push_back (up[nregion-1]);
	
	return result;
};

//-----------------------------------------------------------------------------
// make variable binned hist from a given hist. this method is overloaded
//-----------------------------------------------------------------------------
TH1* Local_MakeVariableBinHist (TH1 *hist, const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4)
{
	assert (hist != NULL && "CommonTools::MakeVariableBinHist:: hist is NULL!");
	assert (hist->GetDimension() == 1 && "CommonTools::MakeVariableBinHist:: hist is not 1-D");
	
  	std::auto_ptr<TH1> result = Local_FillVarBinHist ( 
											Local_GetVarBinVector(xmin, xpoint1, xpoint2, 
															xpoint3, xpoint4, width1, 
															width2, width3, width4)
											, hist);
  	Local_NormalizeBinWidth(result.get());
  	return result.release ();
};




void InvMassFitData(const int plot=1, const int test=0) {

	gStyle->SetOptFit(1);
	gSystem->Load("~/samantha/RootTools/CommonTools_cc.so");
	TFile *f = new TFile ("PhoJets_data.root");
	assert(f != NULL && "file is null");

	TH1* hist = 0;
	if (plot == 1) hist = dynamic_cast<TH1*> (f->Get("Hist/SIGNAL/1Jet/PhotonLeadJet/InvMass"));
	else if (plot == 2) hist = dynamic_cast<TH1*> (f->Get("Hist/SIGNAL/2Jet/PhotonLeadJet/InvMass"));
	else if (plot == 3) hist = dynamic_cast<TH1*> (f->Get("Hist/SIGNAL/2Jet/Photon2Jets/InvMass"));
	else if (plot == 4) hist = dynamic_cast<TH1*> (f->Get("Hist/SIGNAL/2Jet/2Jets/InvMass"));
	assert (hist != NULL && "hist is null");

	//std::cout << blue << "Before rebin " << clearatt << std::endl;
	//DumpHist(hist);
	
	//const 1 GeV binning, only cut off the turn-on
	/*if (plot ==1 || plot ==2) hist = Local_MakeVariableBinHist(hist,150,200,400,700,1200,1,1,1,1);
	else if (plot == 3) hist = Local_MakeVariableBinHist(hist,150,250,400,700,1100,1,1,1,1);
	else if (plot == 4) hist = Local_MakeVariableBinHist(hist,100,300,500,800,1050,1,1,1,1);
	hist->Rebin(5);
	*/
	
	if (plot ==1 || plot ==2) hist = Local_MakeVariableBinHist(hist,130,200,400,700,900,10,10,20,75);
	else if (plot == 3) hist = Local_MakeVariableBinHist(hist,150,250,400,700,1100,10,10,20,75);
	else if (plot == 4) hist = Local_MakeVariableBinHist(hist,100,300,400,700,1075,10,10,20,75);

	
	float xmin=0,xmax=1200.;
	
	if (plot ==1 ) { xmin == 50; xmax=950.;}
	if (plot ==3 ) { xmin == 100; xmax=950.;}
	if (plot ==4 ) { xmin == 100; xmax=750.;}
	
	//hist = Local_MakeVariableBinHist(hist,100,200,400,600,900,20,25,50,100);
	//if (plot ==1 || plot ==2) hist = Local_MakeVariableBinHist(hist,100,200,400,700,900,20,25,50,200);
	//else if (plot == 3) hist = Local_MakeVariableBinHist(hist,150,250,400,700,1100,10,20,50,100);
	//else if (plot == 4) hist = Local_MakeVariableBinHist(hist,100,300,500,800,1050,20,50,100,250);

	//new finer bins accord to mass resolution 
	//if (plot ==1 || plot ==2) hist = Local_MakeVariableBinHist(hist,130,200,400,700,900,10,10,20,75);
	//else if (plot == 3) hist = Local_MakeVariableBinHist(hist,150,250,400,700,1100,10,10,20,75);
	//else if (plot == 4) hist = Local_MakeVariableBinHist(hist,100,300,400,700,1075,10,10,20,75);

	//std::cout << red << "after rebin " << std::endl;
	//DumpHist(hist);
	//std::cout << clearatt << std::endl;
	
	//hist->Sumw2();
	//hist->Sumw2();
	//hist->GetXaxis()->SetRangeUser(0,700);
	//hist->Rebin(5);
	TH1* err = (TH1*) hist->Clone("hist_copy");
	TH1* err2 = (TH1*) hist->Clone("hist_copy1");
	TH1* obsexp = (TH1*) hist->Clone("hist_copy2");

	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(1,2);
	c1->cd(1);
	gPad->SetLogy();
	hist->DrawCopy();
	


	//TF1* xx = new TF1("xxf","[0] * pow(1-x/sqrt(1.96),[1])/pow(x/sqrt(1.96),[2])",100.0,700.0);
	//TF1* xx = new TF1("xxf","[0] * pow(x,[1] * x))",100.0,900.0);
	//const float xmin = hist->GetBinCenter(1);
	//const float xmin = hist->GetBinLowEdge(hist->GetMaximumBin());
	//const float xmax = hist->GetBinCenter(hist->GetNbinsX())+100.;
	//const float xmax = hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX());
	//const float xmax = hist->GetBinCenter(hist->GetNbinsX()) + hist->GetBinWidth(hist->GetNbinsX())/2.;
	//const float xmax = hist->GetBinCenter(hist->GetNbinsX());
	//const float xmax = hist->GetBinCenter(FindLastBin(hist));
	//const float xmax = 900.;
	//TF1* xx = new TF1("xxf","[0]/pow(x,[1])",xmin, xmax);
	//TF1* xx = new TF1("xxf","[0]*([2]+[3]*sqrt(x))/pow(x,[1])",xmin, xmax);
	//TF1* xx = new TF1("xxf","[0]*(1.+[2]*sqrt(x))/pow(x,[1])",xmin, xmax);
	//TF1* xx = new TF1("xxf","[0]*(1.+[2]/pow(x/[3],[1]))",xmin, xmax);
	//TF1* xx = new TF1("xxf","([0]/pow(x-50.,[1])) * pow(1.-x/2000,[2])",xmin, xmax);
//	TF1* xx = new TF1("xxf","([0]/pow(x-50.,[1])) * pow(1-x/2000+x*x/1000,[2])",xmin, xmax); //best so far
	//TF1* xx = new TF1("xxf","([0]/pow(x-50.,[1])) * pow(1-x/1960.+x*x/1960.,[2])",xmin, xmax); //best so far
	

	//works for plot 1
	//TF1* xx = new TF1("xxf","([0]/pow(x-50.,[1])) * pow(1-x/1960.+x*x/1960.,[2]) * pow((1400-x)/1960,[3])",xmin, xmax); //best so far
	//xx->SetParameters(0.9 * max, 3.5, 4,6,0.3);
	

	//TF1* xx = new TF1("xxf","[0] * exp([1]*x)/(1+pow(x,[2])/1960)",xmin, xmax);
	//xx->SetParameters(5000, -0.1, 0.1,1.0,-1.0);
	//TF1* xx = new TF1("xxf","[0] *1000.0* exp([1]*x+[4])/(1+pow(x*[3],0.01*[2]*[0]))",xmin, xmax);
	
	
	//try 1
	//TF1* xx = new TF1("xxf","[0] * ([3] * exp([1]*x) + exp([4]*x))/(1.0+pow(x,[2])/20000.)",xmin, xmax);
	//xx->SetParameters(5000, -0.1, 0.0003,1.0,-1.0);
	//xx->SetParameters(50000, -0.1, 0.0003, 5.0, -1.0,25);
	
	//TF1* xx = new TF1("xxf","[0] * exp([1]*x)/(1+pow(x*[3],[2]))",xmin, xmax);
	//xx->SetParameters(5000, -0.1, -4.2,1.0);
	//xx->SetParLimits(3,0,0.5);
	
	
	
	//xx->SetParameters(0.9 * max, -0.1, 4,6,0.3);
	//xx->SetParameters(5000, -0.1, 0.0003,1.0);
	//xx->SetParLimits(2,5,10);
	//xx->SetParLimits(2,0,4);
	//xx->SetParLimits(3,0,10);
	
	//this works for plots 3and 4
	//TF1* xx = new TF1("xxf","[0] * exp([1]*x)/([2]+pow(x,[3]))",xmin, xmax);
	//xx->SetParameters(9., -0.1, 4,6,0.3);
	//xx->SetParLimits(2,0,10);
	//xx->SetParameters(5000, -0.1, 4,2,0.3);
	
	//for plot 1:
	//TF1* xx = new TF1("xxf","[0] * exp([1]*x)/(1+[2]*x*x/1960.*1960.)",xmin, xmax);
	//TF1* xx = new TF1("xxf","[0] * exp([1]*x)/(1+[3]*pow(x,[2]))",xmin, xmax);
//	TF1* xx = new TF1("xxf","[0] / (1 + [3]*pow(x,[2])) ",xmin, xmax);
	//TF1* xx = new TF1("xxf","[0] * (exp(-[1]*x) + [3] * pow(x,[2])) ",xmin, xmax);
	//xx->SetParameters(5000, 0.1, 0.2,1,0.3,1);

	
	//TF1* xx = new TF1("xxf","([0]/pow(x-50.,[1])) * pow(1-x/1960.+x*x/1960.,[2])",xmin, xmax); //best so far
	//TF1* xx = new TF1("xxf","([0]/pow(x-50.,[1])) * pow([3]*x+x*x/(1960.*1960.) + x*x*x/1960 , [2])",xmin, xmax); //best so far
	//xx->SetParLimits(4,1000,50000);
	//TF1* xx = new TF1("xxf","([0]/pow(x,[1])) * pow(1+x*x/1960.,[2])",xmin, xmax); //best so far


//test area
	//TF1* xx = new TF1("xxf","[0] * exp(-x*[1])/(1+ [3]*[2]*[1]*pow(x/1960,[2]))",xmin, xmax); //best so far
	//xx->SetParameters(5000, 0.1, 1,1,0.3,1);

	TF1* xx = 0;
	if (plot == 1)
	{
		if (test == 1)
		{
			//xx = new TF1("xxf","[0] * exp(-x*[1])",xmin, xmax/2); //best so far
			xx = new TF1("xxf","( [0] * exp(-x*[1]) + [2] * exp(-x*[3]))",150, 750); //best so far
			xx->SetParameters(5000,0.024,20,0.0154);
			xx->SetParLimits(1,0,0.05);
			xx->SetParLimits(3,0,0.05);
		} else
		{
			xx = new TF1("xxf","[0] * (exp(-x*[1]) + [2] * exp(-[3]*x))/((1+pow(x,[4])/5000)*(1+pow(x,[5])))",xmin, xmax); //best so far
			xx->SetParLimits(1,0,5);
			xx->SetParLimits(3,0,5);
			xx->SetParLimits(4,2.6,5);
			xx->SetParLimits(5,0.5,2);
		}
		
	} else if (plot == 2)
	{
		if (test == 1)
		{
			//xx = new TF1("xxf","[0] * exp(-x*[1])",xmin, xmax/2); //best so far
			xx = new TF1("xxf","( [0] * exp(-x*[1]) + [2] * exp(-x*[3]))",150, 750); //best so far
			xx->SetParameters(5000,0.024,20,0.0154);
			xx->SetParLimits(1,0,0.05);
			xx->SetParLimits(3,0,0.05);
		} else
		{
			xx = new TF1("xxf","[0] * (exp(-x*[1]) + [2] * exp(-[3]*x))/((1+pow(x,[4])/5000)*(1+pow(x,[5])))",xmin, xmax); //best so far
			xx->SetParLimits(1,0,5);
			xx->SetParLimits(3,0,5);
			xx->SetParLimits(4,2.6,5);
			xx->SetParLimits(5,0.5,2);
		}
	} else if (plot == 3)
	{
		if (test == 1)
		{
			//xx = new TF1("xxf","[0] * exp(-x*[1])",150, 400); //best so far
			xx = new TF1("xxf","[0] * exp(-x*[1]) + [2] * exp(-x*[3])",xmin, xmax); //best so far
			//xx = new TF1("xxf","[0] exp(-x*[1]) + [2] * exp(-x*[3])",xmin, xmax);
			xx->SetParameters(5000,0.01,20,0.017);
			xx->SetParLimits(1,0,0.02);  //1.78262e-02
			xx->SetParLimits(3,0,0.0172);
			
		} else
		{
			xx = new TF1("xxf","[0] * (exp(-x*[1]))/(pow(x,[2]))",xmin, xmax); //best so far
			xx->SetParameter(0,5000);
			xx->SetParLimits(1,0,5);
			xx->SetParLimits(2,1,5);
			xx->SetParLimits(3,1,5);
		}
	} else if (plot == 4)
	{
		if (test == 1)
		{
			//xx = new TF1("xxf","[0] * exp(-x*[1])",xmin, xmax/2); //best so far
			xx = new TF1("xxf","( [0] * exp(-x*[1]) + [2] * exp(-x*[3]))",xmin, xmax); //best so far
			xx->SetParameters(5000,0.024,20,0.0154);
			xx->SetParLimits(1,0,0.05);
			xx->SetParLimits(3,0,0.05);
		} else
		{

			xx = new TF1("xxf","[0] * exp([1]*x)/(1+[2]*x*x/1960)",xmin, xmax);
			xx->SetParameters(5000, -0.1, 4,6,0.3);
		}
	}

	if (xx == 0)
	{
		std::cout << red << __LINE__ << ": xx is not defined ! check. returning." << clearatt << std::endl;
		return;
	}

	//TF1* xx = new TF1("xxf","[0] * (exp(-x*[1]) + [3]*[2] * exp(-x*[2]))",xmin, xmax); //best so far
	//TF1* xx = new TF1("xxf","[0] * (exp(-x*[1]) + [2] * exp(-[3]*x))/(1+pow(x,[4])/5000)/(1+pow(x,[6]))",xmin, xmax); //best so far
	//TF1* xx = new TF1("xxf","[0] * (exp(-x*[1]) + [2] * exp(-[3]*x))/((1+pow(x,[4])/5000)*(1+pow(x,[5])))",xmin, xmax); //best so far
	//xx->SetParLimits(1,0,5);
	//xx->SetParLimits(3,0,5);
	//xx->SetParLimits(4,2.6,5);
	//xx->SetParLimits(5,0.5,2);
	
	//*********************** good ones 06-05-2010
	//for plto 4:
	//TF1* xx = new TF1("xxf","[0] * exp([1]*x)/(1+[2]*x*x/1960)",xmin, xmax);
	//xx->SetParameters(5000, -0.1, 4,6,0.3);

	//for plot 1
	//TF1* xx = new TF1("xxf","([0]/pow(x-50.,[1])) * pow([3]*x+x*x/(1960.*1960.) + x*x*x/1960 , [2])",xmin, xmax); //best so far
	
	//best for plot 2, ok for plot1 -6-9-2010 : see elog:1662
	//TF1* xx = new TF1("xxf","[0] * (exp(-x*[1]) + [2] * exp(-[3]*x))/((1+pow(x,[4])/5000)*(1+pow(x,[5])))",xmin, xmax); //best so far
	//xx->SetParLimits(1,0,5);
	//xx->SetParLimits(3,0,5);
	//xx->SetParLimits(4,2.6,5);
	//xx->SetParLimits(5,0.5,2);


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	hist->Fit(xx,"LR+");
	//hist->Fit(xx,"R+");
	
	
	//if((gMinuit->fCstatu).Contains("SUCCESSFUL"))
	if(gMinuit->GetStatus() == 0)
	{
		TList *funcList = hist->GetListOfFunctions();
		TIterator *it = funcList->MakeIterator();
		TF1 *fitfunc = (TF1*) it->Next();

		std::cout << red << "Fit Results chi2/ndf= " 
			<< fitfunc->GetChisquare() << "/" << fitfunc->GetNDF() << " = " 
			<< fitfunc->GetChisquare() / fitfunc->GetNDF() << clearatt << std::endl; 

		GetError(err,fitfunc);
		GetError(err2,fitfunc,1);
		//this statement does not give an error!
		//err2,fitfunc->SetLineColor(kRed);
		fitfunc->SetLineColor(kRed);
		err->SetLineColor(kBlack);
		err2->SetMarkerColor(kBlue);
		//err2->SetMarkerStyle(21);
		//err2->SetMarkerSize(0.5);

		
		err2->SetLineColor(8);
		err2->SetMarkerColor(8);
		//err2->SetMarkerStyle(21);
		//err2->SetMarkerSize(0.5);
		

		GetObsExpt(obsexp,fitfunc);
		obsexp->SetLineColor(kBlue);
		obsexp->SetMarkerColor(kBlue);
		

		//set a commong y range
		float ymax = std::max(obsexp->GetBinContent(obsexp->GetMaximumBin()), std::max(err->GetBinContent(err->GetMaximumBin()) , err2->GetBinContent(err2->GetMaximumBin())));
		float ymin= std::min(obsexp->GetBinContent(obsexp->GetMinimumBin()), min(err->GetBinContent(err->GetMinimumBin()) , err2->GetBinContent(err2->GetMinimumBin())));
		int ylim = (int) (std::max(fabs(ymax),fabs(ymin)) + 0.5);
	
		//err->SetMaximum((double) (ylim));
		//err->SetMinimum( (double) (-1*ylim));

		err->SetMaximum(1.);
		err->SetMinimum(-1.);
		
		
		
		c1->cd(2);
		gStyle->SetGridStyle(7);
		gPad->SetGridx();
		gPad->SetGridy();
		err->SetStats(0);
		err2->SetStats(0);
		obsexp->SetStats(0);
		err->DrawCopy("P");
		//err2->DrawCopy("Psame");
		obsexp->DrawCopy("Psame");


		TLegend *leg = new TLegend (0.55,0.82,0.9,1.0);
		leg->SetBorderSize(1);
		leg->SetTextFont(42);
		leg->SetTextSize(0.04);

		leg->AddEntry(err,"Data/Fit - 1");
		//leg->AddEntry(err2,"(Data - Fit)/ Data Stat. Err");
		leg->AddEntry(obsexp,"(Obs. - expt) /  #sqrt{expt}");
		leg->Draw();

		c1->cd();

		if (plot == 1) c1->Print("pj1_pj1Mass.eps");
		else if (plot == 2) c1->Print("pj2_pj1Mass.eps");
		else if (plot == 3) c1->Print("pj2_pj1j2Mass.eps");
		else if (plot == 4) c1->Print("pj2_j1j2Mass.eps");
		
		MakeFreqPlot(obsexp, plot);
		
		
	} else 
	{
		std::cout << "Fit failed! Skipping residual plot generation" << std::endl;
		return;
	}
		
}
