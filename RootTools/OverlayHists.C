#include <iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TStyle.h"
#include <sstream>
#include <algorithm>
#include "OverlayHists.h"
#include <iomanip>

//const float fTITLE_FONT = 102; //courier
const float fTITLE_FONT = 42; //arial
const float fLABEL_FONT = 102;
const float fAXIS_TITLE_FONT = 102;

/***************************************************/
// Created this to overlay two hists quikly.
// Not complete yet. I want to fully automate this 
// so it will figure out the type of the hist and
// do the right thing.
/***************************************************/

OverlayHists::OverlayHists()
{
	bMakeIndividual = false;
	bNormBgToData = false;
	bBinWtoYtitle = false;
	lineColor1 = kRed;
	lineColor2 = kBlue;
	lineWidth = 2;
	markerStyle1 = 20;
	markerStyle2 = 20;
	drawOptions += "L";

}

void OverlayHists::Two1DHists(const TH1* hist1, const TH1* hist2,
				const std::string legend1, const std::string legend2,
				const std::string finalTitle,
				const std::string xtitle, const std::string ytitle,
				const bool bNorm, const bool bIndividual,
				const int rebin, const double fLoX, const double fUpX)
{
	if (hist1 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! First hist is null!"<< std::endl;
		exit (1);
	}
	if (hist2 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! Second hist is null!"<< std::endl;
		exit (1);
	}

	gStyle->SetOptStat("neouirms");
	gStyle->SetOptFit(1);
	TH1* h1 = dynamic_cast<TH1*> (hist1->Clone());
	TH1* h2 = dynamic_cast<TH1*> (hist2->Clone());
	//h1->SetDirectory(0);
	//h2->SetDirectory(0);
	h1->Rebin(rebin);
	h2->Rebin(rebin);

h1->GetXaxis()->SetLabelFont(fTITLE_FONT);
h1->GetXaxis()->SetTitleFont(fTITLE_FONT);
h1->GetYaxis()->SetLabelFont(fTITLE_FONT);
h1->GetYaxis()->SetTitleFont(fTITLE_FONT);
h2->GetXaxis()->SetLabelFont(fTITLE_FONT);
h2->GetXaxis()->SetTitleFont(fTITLE_FONT);
h2->GetYaxis()->SetLabelFont(fTITLE_FONT);
h2->GetYaxis()->SetTitleFont(fTITLE_FONT);

//h1->GetTitle()->SetTextFont(fTITLE_FONT);
//h2->GetTitle()->SetTextFont(fTITLE_FONT);

	h1->SetLineColor(lineColor1);
	h1->SetMarkerColor(lineColor1);
	h1->SetMarkerSize(1);
	//h1->SetMarkerStyle(markerStyle1);
	h1->SetLineWidth(lineWidth);
	if (bIndividual)
	{
		new TCanvas;
		gPad->SetTickx();
		gPad->SetTicky();
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogy();
		h1->Draw();
		//h1->Fit("gaus","","",-0.1,0.1);
	}
	
	h2->SetLineColor(lineColor2);
	h2->SetMarkerColor(lineColor2);
	h2->SetMarkerSize(0.4);
	//h2->SetMarkerStyle(markerStyle2);
	h2->SetLineWidth(lineWidth);
	if (bIndividual)
	{
		new TCanvas;
		gPad->SetTickx();
		gPad->SetTicky();
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogy();
		h2->Draw();
		//h2->Fit("gaus","","",-0.1,0.1);
	}
	
	new TCanvas();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
	gStyle->SetTitleFont(102,"XYZ");
	gStyle->SetLabelFont(102,"XYZ");
	
	TLegend *leg = new TLegend (0.60,0.75,0.9,0.9);
	leg->SetTextFont(fTITLE_FONT);
	leg->SetTextSize(0.04);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);
	leg->AddEntry(h1, legend1.c_str());
	leg->AddEntry(h2, legend2.c_str());


	// Need to clone these to remove the stat box in the final
	// overlayed hist
	std::stringstream h1copyname;
	h1copyname << h1->GetName() << "_copy";
	TH1* h1copy = dynamic_cast<TH1*> (h1->Clone(h1copyname.str().c_str()));
	std::stringstream h2copyname;
	h2copyname << h2->GetName() << "_copy";
	TH1* h2copy = dynamic_cast<TH1*> (h2->Clone(h2copyname.str().c_str()));

	if(finalTitle.length()>0) 
	{
		h1copy->SetTitle(finalTitle.c_str());
		h2copy->SetTitle(finalTitle.c_str());
	}
	if(xtitle.length()>0) 
	{
		h1copy->GetXaxis()->SetTitle(xtitle.c_str());
		h2copy->GetXaxis()->SetTitle(xtitle.c_str());
	}
	if(ytitle.length()>0) 
	{
		h1copy->GetYaxis()->SetTitle(ytitle.c_str());
		h2copy->GetYaxis()->SetTitle(ytitle.c_str());
	}
	

	float max1 = h1copy->GetMaximum();
	float max2 = h2copy->GetMaximum();
	float maxy = 0;
	if (max1>max2) maxy = max1;
	else maxy = max2;
	
	maxy *= 100.;

	h1copy->SetMaximum(maxy);
	h2copy->SetMaximum(maxy);
	h1copy->SetMinimum(0.5); //do not make this zero if you want to keep log scale. ymin>0 for log scale
	h2copy->SetMinimum(0.5);

	if (fLoX> -9999. && fUpX>-9999.)
	{ 
		h1copy->GetXaxis()->SetRangeUser(fLoX, fUpX);
		h2copy->GetXaxis()->SetRangeUser(fLoX, fUpX);

		std::cout << "X range for at loX (=" 
			<< fLoX << ") or UpX (=" << fUpX 
			<< ")"  << std::endl;
	} else 
	{
		std::cout << "Invalid X range for at loX (=" 
			<< fLoX << ") or UpX (=" << fUpX 
			<< "). check again! X range not set." << std::endl;
	}

	h1copy->SetStats(0);
	h2copy->SetStats(0);
	
	if (bBinWtoYtitle)
	{
		std::stringstream yt;
		yt << h1copy->GetYaxis()->GetTitle() << "   " 
			<< "Events/" << h1copy->GetBinWidth(1);
		h1copy->GetYaxis()->SetTitle(yt.str().c_str());
	}

	if (bNorm)
	{
		double h1Int = h1copy->Integral();
		double h2Int = h2copy->Integral();
		std::cout << "Normalizing scale = " << h1Int/h2Int << std::endl;
		h2copy->Scale(h1Int/h2Int);
	}
	

	//h1copy->Draw("E P");
	h1copy->Draw();
	h2copy->Draw("SAME");
	leg->Draw();

}


void OverlayHists::TwoProfileHists(const TProfile* hist1, const TProfile* hist2
							, const std::vector<std::string> legend, const std::string finalTitle)
{
	if (hist1 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! First hist is null!"<< std::endl;
		exit (1);
	}
	if (hist2 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! Second hist is null!"<< std::endl;
		exit (1);
	}
	if (legend.size() != 2)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! Too many or too little Entries for Legend!"<< std::endl;
		exit (1);
	}

	gStyle->SetOptStat("neouirms");
	gStyle->SetOptFit(1);
	TProfile* h1 = dynamic_cast<TProfile*> (hist1->Clone());
	TProfile* h2 = dynamic_cast<TProfile*> (hist2->Clone());
	h1->SetDirectory(0);
	h2->SetDirectory(0);

	new TCanvas;
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	h1->SetLineColor(kBlue);
	h1->SetMarkerColor(kBlue);
	h1->SetMarkerSize(0.4);
	h1->SetMarkerStyle(20);
	gPad->SetGridy();
	h1->Draw();
	//h1->Fit("gaus","","",-0.2,0.2);
	
	std::stringstream h1file;
	h1file<< h1->GetName() << ".gif";
	//gPad->Print(h1file.str().c_str());

	new TCanvas;
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	h2->SetLineColor(kRed);
	h2->SetMarkerColor(kRed);
	h2->SetMarkerSize(0.4);
	h2->SetMarkerStyle(20);
	h2->Draw();
	//h2->Fit("gaus","","",-0.2,0.2);
	std::stringstream h2file;
	h2file<< h2->GetName() << ".gif";
	//gPad->Print(h2file.str().c_str());
	
	new TCanvas();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	
	TLegend *leg = new TLegend (0.65,0.75,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);
	leg->AddEntry(h1, legend.at(0).c_str());
	leg->AddEntry(h2, legend.at(1).c_str());


	// Need to clone these to remove the stat box in the final
	// overlayed hist
	std::stringstream h1copyname;
	h1copyname << h1->GetName() << "_copy";
	TH1* h1copy = dynamic_cast<TH1*> (h1->Clone(h1copyname.str().c_str()));
	std::stringstream h2copyname;
	h2copyname << h2->GetName() << "_copy";
	TH1* h2copy = dynamic_cast<TH1*> (h2->Clone(h2copyname.str().c_str()));

	h1copy->SetTitle(finalTitle.c_str());

	float max1 = h1copy->GetMaximum();
	float max2 = h2copy->GetMaximum();
	float maxy = 0;
	if (max1>max2) maxy = max1;
	else maxy = max2;
	
	maxy *= 1.1;

	h1copy->SetMaximum(maxy);
	h2copy->SetMaximum(maxy);
	h1copy->SetStats(0);
	h2copy->SetStats(0);
	
	h1copy->Draw();
	h2copy->Draw("same");
	leg->Draw();
	//gPad->Print("Overlayed.gif");

}


void OverlayHists::Three1DHists(const TH1* hist1, const TH1* hist2, const TH1* hist3,
				const std::vector<std::string> legend, const std::string finalTitle,const int rebin)
{
	if (hist1 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! First hist is null!"<< std::endl;
		exit (1);
	}
	if (hist2 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! Second hist is null!"<< std::endl;
		exit (1);
	}
	if (hist3 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! Third hist is null!"<< std::endl;
		exit (1);
	}
	if (legend.size() != 3)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! Too many or too little Entries for Legend!"<< std::endl;
		exit (1);
	}

	gStyle->SetOptStat("neouirms");
	gStyle->SetOptFit(1);
	TH1* h1 = dynamic_cast<TH1*> (hist1->Clone());
	TH1* h2 = dynamic_cast<TH1*> (hist2->Clone());
	TH1* h3 = dynamic_cast<TH1*> (hist3->Clone());
	h1->SetDirectory(0);
	h2->SetDirectory(0);
	h3->SetDirectory(0);
	h1->Rebin(rebin);
	h2->Rebin(rebin);
	h3->Rebin(rebin);

	new TCanvas;
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	h1->SetLineColor(kBlue);
	h1->SetMarkerColor(kBlue);
	h1->SetMarkerSize(0.4);
	h1->SetMarkerStyle(20);
	gPad->SetGridy();
	h1->Draw();
	
	new TCanvas;
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	h2->SetLineColor(kRed);
	h2->SetMarkerColor(kRed);
	h2->SetMarkerSize(0.4);
	h2->SetMarkerStyle(20);
	h2->Draw();
	
	new TCanvas;
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	h3->SetLineColor(kGreen);
	h3->SetMarkerColor(kGreen);
	h3->SetMarkerSize(0.4);
	h3->SetMarkerStyle(20);
	h3->Draw();


	new TCanvas();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	
	TLegend *leg = new TLegend (0.65,0.75,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);
	leg->AddEntry(h1, legend.at(0).c_str());
	leg->AddEntry(h2, legend.at(1).c_str());
	leg->AddEntry(h3, legend.at(2).c_str());


	// Need to clone these to remove the stat box in the final
	// overlayed hist
	std::stringstream h1copyname;
	h1copyname << h1->GetName() << "_copy";
	TH1* h1copy = dynamic_cast<TH1*> (h1->Clone(h1copyname.str().c_str()));
	std::stringstream h2copyname;
	h2copyname << h2->GetName() << "_copy";
	TH1* h2copy = dynamic_cast<TH1*> (h2->Clone(h2copyname.str().c_str()));
	std::stringstream h3copyname;
	h2copyname << h3->GetName() << "_copy";
	TH1* h3copy = dynamic_cast<TH1*> (h3->Clone(h3copyname.str().c_str()));


	h1copy->SetTitle(finalTitle.c_str());

	float max1 = h1copy->GetMaximum();
	float max2 = h2copy->GetMaximum();
	float max3 = h3copy->GetMaximum();
	float maxy = 0;
	if (max1>max2) maxy = max1;
	else maxy = max2;
	if (maxy<max3) maxy = max3;
	
	maxy *= 1.1;

	h1copy->SetMaximum(maxy);
	h2copy->SetMaximum(maxy);
	h3copy->SetMaximum(maxy);
	h1copy->SetStats(0);
	h2copy->SetStats(0);
	h3copy->SetStats(0);
	
	h1copy->Draw();
	h2copy->Draw("same");
	h3copy->Draw("same");
	leg->Draw();

}




void OverlayHists::ThreeProfileHists(const TProfile* hist1, const TProfile* hist2, const TProfile* hist3
							, const std::vector<std::string> legend, const std::string finalTitle)
{

	if (hist1 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! First hist is null!"<< std::endl;
		exit (1);
	}
	if (hist2 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! Second hist is null!"<< std::endl;
		exit (1);
	}
	if (hist3 == NULL)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! Third hist is null!"<< std::endl;
		exit (1);
	}
	if (legend.size() != 3)
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":FATAL ERROR! Too many or too little Entries for Legend!"<< std::endl;
		exit (1);
	}

	gStyle->SetOptStat("neoui");
	gStyle->SetOptFit(1);
	TProfile* h1 = dynamic_cast<TProfile*> (hist1->Clone());
	TProfile* h2 = dynamic_cast<TProfile*> (hist2->Clone());
	TProfile* h3 = dynamic_cast<TProfile*> (hist3->Clone());

	new TCanvas;
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	h1->SetLineColor(lineColor1);
	h1->SetMarkerColor(lineColor1);
	h1->SetMarkerSize(0.4);
	h1->SetMarkerStyle(markerStyle1);
	gPad->SetGridy();
	h1->Draw();
	//h1->Fit("gaus","","",-0.2,0.2);
	
	std::stringstream h1file;
	h1file<< h1->GetName() << ".gif";
	//gPad->Print(h1file.str().c_str());

	new TCanvas;
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	h2->SetLineColor(lineColor2);
	h2->SetMarkerColor(lineColor2);
	h2->SetMarkerSize(0.4);
	h2->SetMarkerStyle(markerStyle2);
	h2->Draw();
	//h2->Fit("gaus","","",-0.2,0.2);
	std::stringstream h2file;
	h2file<< h2->GetName() << ".gif";
	//gPad->Print(h2file.str().c_str());
	
	new TCanvas;
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	h3->SetLineColor(kGreen);
	h3->SetMarkerColor(kGreen);
	h3->SetMarkerSize(0.4);
	h3->SetMarkerStyle(20);
	//h3->Fit("gaus","","",-0.2,0.2);
	h3->Draw();
	std::stringstream h3file;
	h3file<< h3->GetName() << ".gif";
	//gPad->Print(h3file.str().c_str());



	new TCanvas();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	
	TLegend *leg = new TLegend (0.65,0.75,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);
	leg->AddEntry(h1, legend.at(0).c_str());
	leg->AddEntry(h2, legend.at(1).c_str());
	leg->AddEntry(h3, legend.at(2).c_str());


	// Need to clone these to remove the stat box in the final
	// overlayed hist
	std::stringstream h1copyname;
	h1copyname << h1->GetName() << "_copy";
	TH1* h1copy = dynamic_cast<TH1*> (h1->Clone(h1copyname.str().c_str()));
	std::stringstream h2copyname;
	h2copyname << h2->GetName() << "_copy";
	TH1* h2copy = dynamic_cast<TH1*> (h2->Clone(h2copyname.str().c_str()));
	std::stringstream h3copyname;
	h3copyname << h3->GetName() << "_copy";
	TH1* h3copy = dynamic_cast<TH1*> (h3->Clone(h3copyname.str().c_str()));

	h1copy->SetTitle(finalTitle.c_str());

	float max1 = h1copy->GetMaximum();
	float max2 = h2copy->GetMaximum();
	float max3 = h3copy->GetMaximum();
	float maxtemp=0,maxy = 0;
	if (max1>max2) maxtemp = max1;
	else maxtemp = max2;
	if (maxtemp>max3) maxy = maxtemp;
	else maxy = max3;
	
	maxy *= 1.1;

	h1copy->SetMaximum(maxy);
	h2copy->SetMaximum(maxy);
	h3copy->SetMaximum(maxy);
	h1copy->SetStats(0);
	h2copy->SetStats(0);
	h3copy->SetStats(0);
	
	h1copy->Draw();
	h2copy->Draw("same");
	h3copy->Draw("same");
	leg->Draw();
	//gPad->Print("Overlayed.gif");

}

