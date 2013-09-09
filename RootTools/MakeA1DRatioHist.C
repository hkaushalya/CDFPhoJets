#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"
#include <sstream>
#include <iomanip>

using namespace std;

void MakeA1DRatioHist(const TH1 *data, const TH1 *bg,
				const int rebin,
				const std::vector<std::string> legend,
				const std::string title="", 
				const std::string xtitle="",
				const std::string ytitle="",
				
				const int debug=0)
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

	TH1 *hist_data = dynamic_cast<TH1*> (data->Clone()); 
	hist_data->SetDirectory(0);
	TH1 *hist_bg = dynamic_cast<TH1*> (bg->Clone()); 
	hist_bg->SetDirectory(0);

	
	// verify of the two hists have the same number of bins and bin widths 
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
	

	// rebin hists
	hist_data->Rebin(rebin);
	hist_bg->Rebin(rebin);

	// now to make the ratio plots
	if (debug)	std::cout << "bin\tdata,err\tbg,err\tscale\t\t\tdata,err\tbg,err" <<std::endl;
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
		if (debug)
		{
			std::cout << bin << "\t" << valdata << "," << errdata << "\t" 
						<< val << "," << error << "\t" << scale 
						<< "\t\t\t" << hist_data->GetBinContent(bin) 
						<< "," << hist_data->GetBinError(bin) << "\t" 
						<< hist_bg->GetBinContent(bin) << "," 
						<< hist_bg->GetBinError(bin) << std::endl; 
		}
	};
	
	
	gStyle->SetOptStat("");
	int labelfont = 10 * 4 + 2;		//10 * font ID + precision (2 = scalable)
	int titlefont = 10 * 4 + 2;		//10 * font ID + precision (2 = scalable)
	gStyle->SetLabelFont(labelfont,"X");
	gStyle->SetLabelFont(labelfont,"Y");
	gStyle->SetTitleFont(titlefont,"X");
	gStyle->SetTitleFont(titlefont,"Y");
	gStyle->SetLabelSize(0.04,"X");
	gStyle->SetLabelSize(0.027,"Y");
	//gStyle->SetLabelOffset(0.9);
	gStyle->SetTitleSize(0.03,"Y");
	gStyle->SetTitleOffset(1.8,"Y");


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
	hist_data->SetMaximum(1);
	hist_data->SetMinimum(-1);


	if (title.length()) hist_bg->SetTitle(title.c_str());
	if (xtitle.length()) hist_bg->GetXaxis()->SetTitle(xtitle.c_str());
	std::stringstream Ytitle;
	if (ytitle.length())
	{
		Ytitle << ytitle << " : Bin Width= " <<  setw(3) << hist_bg->GetBinWidth(1) << std::endl;
	} else {
		Ytitle << " Bin Width= " <<  setw(3) << hist_bg->GetBinWidth(1) << std::endl;
	}
	hist_bg->GetYaxis()->SetTitle(Ytitle.str().c_str());

	hist_bg->Draw("E2");
	hist_data->Draw("same");


	TLegend *leg = new TLegend (0.65,0.75,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);
	leg->AddEntry(hist_data, legend.at(0).c_str());
	leg->AddEntry(hist_bg, legend.at(1).c_str());

	leg->Draw();

}

