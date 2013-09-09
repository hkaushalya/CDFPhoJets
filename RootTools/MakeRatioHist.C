#include "MakeRatioHist.h"
#include<iostream>
#include<TFile.h>
#include<TCanvas.h>
#include<TStyle.h>
#include<TLegend.h>
#include<TF1.h>
#include<sstream>
#include<cmath>
/****************************************************/
// This will return a ratio hist of give data and bckg hists
// some basic sanity check are done to make sure the two histograms
// have the same number of bins and the bin edges.
/****************************************************/


MakeRatioHist::MakeRatioHist()
{
}

float MakeRatioHist::CorrelatedErrors(const float val1, const float err1,
								const float val2, const float err2)
{
	//this will calculated the error bars if the two hists have
	// correlated errors
	// for eg: when overlaying just the background predictions from 
	// two different jobs
	// X = A/B - 1
	// delX = X * sqrt ((delA/A)^2 + X * (delB/B)^2);
	float X = val1/val2;
	return X * sqrt( pow(err1/val1,2) + pow(X * err2/val2,2) );

}

void MakeRatioHist::RatioHist(const TH1 *data, const TH1 *bg,
							const std::string sTitle)
{
	if (data == NULL)
	{
		std::cout << __FILE__ << "::" << __LINE__ << ":: data hist pointer is null! exiting!" << std::endl; 
		exit (1);
	}
	if (bg == NULL)
	{
		std::cout << __FILE__ << "::" << __LINE__ << ":: bckg hist pointer is null! exiting!" << std::endl; 
		exit (1);
	}

	TH1 *hist_data = dynamic_cast<TH1*> (data->Clone("data_copy")); 
	TH1 *hist_bg = dynamic_cast<TH1*> (bg->Clone("bg_copy")); 

	
	if (sTitle.length())
	{
		hist_data->SetTitle(sTitle.c_str());
		hist_bg->SetTitle(sTitle.c_str());
	}

	if (hist_data->GetNbinsX() != hist_bg->GetNbinsX())
	{
		std::cout << __FILE__ << "::" << __LINE__ << ":: two hists have different number of bins. cannot make ratio hist! exiting!" << std::endl; 
		exit (1);
	}


	for (unsigned i = 1; i <= (unsigned) hist_data->GetNbinsX(); ++i)
	{
		if (hist_data->GetBinCenter(i) != hist_bg->GetBinCenter(i))
		{
			std::cout << __FILE__ << "::" << __LINE__ << ":: two hists have different bins sizes. cannot make ratio hist! exiting!" << std::endl; 
		exit (1);
		}
	}
	
	// now to make the ratio plots
	//std::cout << "bin\tdata,err\tbg,err\tscale\t\t\tdata,err\tbg,err" <<std::endl;
	for (unsigned bin = 0; bin <= (unsigned) hist_data->GetNbinsX() + 1; ++ bin)
	{
		const float val = hist_bg->GetBinContent (bin);
		const float error = hist_bg->GetBinError (bin);
		const float scale = val ? 1. / val : 0;
		const float valdata = hist_data->GetBinContent (bin);
		const float errdata = hist_data->GetBinError (bin);
		hist_data->SetBinContent (bin, (hist_data->GetBinContent (bin) - val) * scale);
		hist_data->SetBinError (bin, hist_data->GetBinError (bin) * scale);
		hist_bg->SetBinError (bin, val ? error / val : 0);
		hist_bg->SetBinContent (bin, 0);
	//	if (bin%2)
	//	std::cout << bin << "\t" << valdata << "," << errdata << "\t" << val << "," << error << "\t" << scale << "\t\t\t" << hist_data->GetBinContent(bin) << "," << hist_data->GetBinError(bin) << "\t" << hist_bg->GetBinContent(bin) << "," << hist_bg->GetBinError(bin) << std::endl; 
	};


	gStyle->SetOptStat(0);
  	hist_bg->SetFillStyle(3002);
  	hist_bg->SetFillColor(kRed);
	hist_data->SetLineColor(kBlue);
	hist_data->SetMarkerColor(kBlue);
	hist_data->SetMarkerStyle (8);
	new TCanvas();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	//hist_bg->GetYaxis()->SetRangeUser(-1,1);
	//hist_data->GetYaxis()->SetRangeUser(-1,1);
	hist_bg->SetMaximum(1);
	hist_bg->SetMinimum(-1);
	//hist_bg->GetYaxis()->SetTitle("DATA - MC / MC");
	hist_bg->GetYaxis()->SetTitle("PYTHIA MC/ MEt Model - 1");
	hist_bg->Draw("E2");
	hist_data->Draw("same");
	//line along 1
	TF1 *line = new TF1("line","0",-1500,1500);
	line->SetLineWidth(2);
	line->Draw("same");
	
	

}
void MakeRatioHist::TestMethod(const int njet)
{
	//assert ( (njet>=1 && njet<=4) && "njets invalid. must be 1-4!"); 
	gStyle->SetOptStat(0);	
	//TFile *f = new TFile("AnaMode_njet4.root");
	TFile *f = new TFile("MetTest_100K.root");
	//TFile *g = new TFile("CalibMode.root");
	f->ls();
	std::string ana_path("/Ana/JetFilter2/Hist/Ana_data/");
	std::string calib_path("/Ana/JetFilter2/Hist/Ana_bckg/");
	//std::string calib_path("/Ana/JetFilterV2/Hist/METJetClu0.4_scen3/");
	
	std::stringstream histname1, histname2;
	//histname1 << "fMetSig_estimate_njet" << njet; 
	histname1 << "Met"; 
	histname2 << "Met"; 
	//histname2 << "fMetSigCalib_estimate_njet" << njet; 
	
	std::cout << "Searching for " << histname1.str() << " & " << histname2.str() << std::endl;

	f->cd();
	gDirectory->cd(ana_path.c_str());
	gDirectory->ls();
	TH1 *data = dynamic_cast<TH1*> (gDirectory->FindObjectAny(histname1.str().c_str()));
	assert (data != NULL && "data hist is null");

	//g->cd();
	f->cd();
	gDirectory->cd(calib_path.c_str());
	TH1 *bg = dynamic_cast<TH1*> (gDirectory->FindObjectAny(histname2.str().c_str()));
	assert (bg != NULL && "bg hist is null");

	//TH1 *data = dynamic_cast<TH1*> (f->Get(histname1.str().c_str()));
	//TH1 *bg = dynamic_cast<TH1*> (f->Get(histname2.str().c_str()));
	assert (data != NULL && "data hist is null");
	assert (bg != NULL && "bg hist is null");
	data->SetDirectory(0);
	bg->SetDirectory(0);
	int rebin = 4;
	data->Rebin(rebin);
	bg->Rebin(rebin);
	std::cout << data->GetBinWidth(1) << std::endl;
	std::cout << bg->GetBinWidth(1) << std::endl;
	
	std::stringstream title;
	title<< data->GetTitle();
	data->SetTitle(title.str().c_str());

	//normalize to unity
	double d1 = data->Integral();
	double b1 = bg->Integral();

	//data->Scale(1.0/(1. * d1));
	//bg->Scale(1.0/(1. * b1));

	data->Print();
	bg->Print();

	TCanvas *c = new TCanvas("c","Log And Ratio Hists",1200,600);
	c->Divide(2,1);
	c->cd(1);
	LogHist(data,bg);
	c->cd(2);
	RatioHist(data,bg);
	c->cd();

}

void MakeRatioHist::LogHist(const TH1 *data, const TH1 *bg)
{	
	TH1 *data1 = dynamic_cast<TH1*> (data->Clone("data_copy")); 
	TH1 *bg1 = dynamic_cast<TH1*> (bg->Clone("bg_copy")); 
	TH1 *bg2 = dynamic_cast<TH1*> (bg->Clone("bg_copy_copy")); 

	bg1->SetMarkerColor(kBlack);
	bg1->SetLineColor(kBlack);
  	bg1->SetFillStyle(3002);
	bg1->SetFillColor(kRed);
	data1->SetLineColor(kBlue);
	data1->SetMarkerColor(kBlue);
	data1->SetMarkerStyle (8);
	data1->SetMarkerSize(1.0);
	bg2->SetLineWidth(2);

	data1->SetMinimum(0.005);
	bg1->SetMinimum  (0.005);

	//new TCanvas();
	//gPad->SetLogy();
	bg1->Draw("E2");
	data1->Draw("samePE");
	//bg2->Draw("same");
	TLegend *leg = new TLegend (0.7,0.8,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.035);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);
	leg->AddEntry(data1,"DATA");
	leg->AddEntry(bg1, "MC");
	leg->Draw();
	
}
void MakeRatioHist::TestMethod2(const std::string name, const std::string filename)
{
	gStyle->SetOptStat(0);	
	TFile *f = new TFile("old_result.root");
	f->ls();
	TFile *g = new TFile("new_result.root");
	g->ls();
	
	f->cd();
	TH1 *data = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name.c_str()));
	assert (data != NULL && "data hist is null");

	g->cd();
	TH1 *bg = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name.c_str()));
	assert (bg != NULL && "bg hist is null");

	assert (data != NULL && "data hist is null");
	assert (bg != NULL && "bg hist is null");
	data->SetDirectory(0);
	bg->SetDirectory(0);

	std::cout << data->GetBinWidth(1) << std::endl;
	std::cout << bg->GetBinWidth(1) << std::endl;
	
	TCanvas *c = new TCanvas("c","Log And Ratio Hists",1200,600);
	c->Divide(2,1);
	c->cd(1);
	LogHist(data,bg);
	c->cd(2);
	RatioHist(data,bg);
	c->cd();
	gPad->Print(filename.c_str());

}


