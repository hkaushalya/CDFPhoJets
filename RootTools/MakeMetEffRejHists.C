#include<iostream>
#include<sstream>
#include<TFile.h>
#include<TCanvas.h>
#include<TH1.h>
#include<TStyle.h>
#include<string>
#include<algorithm>
#include<TLegend.h>
#include<TStyle.h>
#include<TSystem.h>
#include<TDirectory.h>
#include<TF1.h>
#include<TMath.h>

using namespace std;


float GetMaxX(const TH1* htemp_MetAll)
{
	assert (htemp_MetAll != NULL && "MetAll hist is NULL!");
	
	float fX = 0;
	for (unsigned i = 1; i<= htemp_MetAll->GetNbinsX(); ++i)
	{
		double val = htemp_MetAll->GetBinContent(i);
		if (val) fX = htemp_MetAll->GetBinCenter(i);
		//std::cout << val << "\t" << fX << std::endl;
	}
	//now add few more bins of width
	//std::cout << "fX= " << fX << std::endl;
	fX = fX + 5 * htemp_MetAll->GetBinWidth(1);
	return fX;
}

//
void CalculateBinomialErrors(TH1* hRatio, TH1* hMetBefore)
{
	assert (hRatio != NULL && "Ratio hist is NULL!");
	assert (hMetBefore != NULL && "MetBefore hist is NULL!");

	std::cout << "bin\tolderr\tval\tprob\terr\tvariance" << std::endl;
	for (unsigned i = 1; i<= hRatio->GetNbinsX(); ++i)
	{
		float olderr = hRatio->GetBinError(i);
		// err = sqrt(Np(1-p))/sqrt(N)= sqrt(p(1-p))
		double val = hMetBefore->GetBinContent(i);
		if (val == 0 ) continue;
		float prob = hRatio->GetBinContent(i);
		float variance = prob * (1 - prob);
		float err = TMath::Sqrt(variance);
		hRatio->SetBinError(i,err);
		//std::cout << i << "\t" << olderr << "\t" << val << "\t" << prob << "\t" << err << "\t" << variance<< std::endl;
	}
}


//=====================================================================
//-------------- fit function: c3 + (c0 - c3)(1 + exp(-(x-c1)/c2)
double Func1(double *x, double *par)
{
	double fitval = par[3] + (par[0] - par[3]) * (1 + exp(-1.0 * (x[0] - par[1])/par[2]));
	return fitval;
}

// effects: normalize each bin by the bin width
// side effect: normalizes histogram bins in hist
// guarantee: no-throw
// requires: hist != NULL
// requires: hist->GetDimension() == 1
void norm_bins (TH1 *hist)
{
	assert (hist != NULL && "requirement failed"); //spec
	assert (hist->GetDimension() == 1 && "requirement failed"); //spec

	for (unsigned bin = 1; bin != unsigned (hist->GetNbinsX()); ++ bin)
	{
		const float width = hist->GetXaxis()->GetBinWidth (bin);
		hist->SetBinContent (bin, hist->GetBinContent (bin) / width);
		hist->SetBinError (bin, hist->GetBinError (bin) / width);
	};
};



// effects: make a new binning with variable bin sizes, and the given
//   cutoff points between bin sizes
// returns: the binning (by bin edges)
// guarantee: strong
//// throws: std::bad_alloc
std::vector<float> 
GetVarBinVector (float down, float up5, float up10, float up20, float up50, 
					float width1, float width2, float width3, float width4 ,bool debug)
{
	std::vector<float> result;
	float point = down;
	const unsigned nregion = 3;
	const float up [] = {up5, up10, up20};
	const float step [] = {width1, width2, width3};

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
// effects: turn this histogram into a histogram with variable bin
//   size
// returns: the new histogram
// guarantee: strong
// throws: std::bad_alloc
// requires: input != NULL
// requires: input->GetDimension() == 1
// requires: histogram bin edges must match
// postcondition: result != NULL
std::auto_ptr<TH1>
make_var_bin (const std::vector<float>& bins,
	      TH1 *input, bool debug)
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


	if (debug)
	{
		for (unsigned i=1;  i<=(unsigned)input->GetNbinsX (); i++)
		{
			std::cout << "input bin, lowedge = " << i << ", " << input->GetXaxis()->GetBinLowEdge(i) << "\t" << input->GetBinContent(i) << std::endl;
		}
		std::cout << "\n\n"<< std::endl;
		for (unsigned i=1; i<=(unsigned)result->GetNbinsX (); i++)
		{
			std::cout << "otput bin, lowedge = " << i << ", " << result->GetXaxis()->GetBinLowEdge(i) << "\t" << result->GetBinContent(i) << std::endl;
		}
	}

	assert (result.get() != NULL && "postcondition failed"); //spec
	return result;
};

// make a variable binned hist with given specification
TH1* MakeVariableBins (TH1 *hist, const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4 , bool use_errors=false, bool debug=false)
{
	if (debug)
	{
		std::cout << hist->GetName() << " hist bin size=" << hist->GetBinWidth(1) << std::endl;
	}

  	std::auto_ptr<TH1> result = 
		make_var_bin ( GetVarBinVector(xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, debug) , hist, debug);
	
  	norm_bins (result.get());
  	return result.release ();
};





////////////////////////////////////////////////////////
// Created this make the final
// MET model Rejection Power and Efficiency . This adds 
// up MetAll (Met before MetSig cut) and Met (Met after 
// met sig cut) for data and background.
// And make the final plots
////////////////////////////////////////////////////////

void MakeMetEffRejHists(const bool data = 0, const bool rebin = true)
//data = 0 == DATA
//data = 1 == BACKGROUND PREDICTION
{
	//data = 0 == background prediction
	//data = 1 == data
	
	TFile *f[6];
	TH1 *hist_MetRatio[6];

	
	f[0] = new TFile ("Merged_MetSig2.root");
	f[1] = new TFile ("Merged_MetSig3.root");
	f[2] = new TFile ("Merged_MetSig4.root");
	f[3] = new TFile ("Merged_MetSig5.root");
	f[4] = new TFile ("Merged_MetSig6.root");
	f[5] = new TFile ("Merged_MetSig7.root");

	float fGlobalMaxX = 0;

	for (int i=0; i<6; i++)
	{
		if (f[i]->IsZombie()) 
		{
			std::cout << "ERROR::File "<< i << " did not open! " << std::endl;
			return;
		} else {
			f[i]->Print();
			//f[i]->ls();

			std::string name = "Met";
			std::string name1 = "MetAll";

			f[i]->cd();
			gDirectory->pwd();

			TH1* htemp_Met = (TH1*) gDirectory->FindObjectAny(name.c_str());
			TH1* htemp_MetAll = (TH1*) gDirectory->FindObjectAny(name1.c_str());
			
			htemp_Met->Rebin(10);
			htemp_MetAll->Rebin(10);

			float xmin = 0;
			float xpoint1 = 20;
			float xpoint2 = 100;
			float xpoint3 = 400;
			float xpoint4 = 500;
			float width1 = 10;
			float width2 = 10;
			float width3 = 50;
			float width4 = 100;
			
			if (rebin)
			{
				htemp_Met = (TH1F*) MakeVariableBins (htemp_Met, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, true, false);
				htemp_MetAll = (TH1F*) MakeVariableBins (htemp_MetAll, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, true, false);
			}
/*			std::cout << "\tMET ALL(BIN/ERR)\tMET (BIN\ERR)" <<std::endl; 
			for (unsigned j = 1; j<=htemp_Met->GetNbinsX(); ++j)	
			{
				std::cout << j << "\t" << htemp_MetAll->GetBinContent(j) << "/" << htemp_MetAll->GetBinError(j) 
						<< "\t\t" << htemp_Met->GetBinContent(j) << "/" << htemp_Met->GetBinError(j) << std::endl; 
			}
*/
	//std::cout << "Met/MetAll Entries = " << htemp_Met->GetEntries() << " / " << htemp_MetAll->GetEntries() << std::endl; 

			assert (htemp_Met != NULL && "Met object not found!");
			assert (htemp_MetAll != NULL && "MetAll object not found!");

			//htemp_Met->Sumw2();
			//htemp_MetAll->Sumw2();

			hist_MetRatio[i] = (TH1*) htemp_Met->Clone("MetAfterCopy");
			hist_MetRatio[i]->SetDirectory(0);
			hist_MetRatio[i]->Divide(htemp_MetAll);
			if (data) // for background errors are already calculated at run-time
			{
				CalculateBinomialErrors(hist_MetRatio[i], htemp_Met);
			}
			float fMaxX = GetMaxX(hist_MetRatio[i]);
			if (fMaxX > fGlobalMaxX) fGlobalMaxX = fMaxX;
			//new TCanvas();
			//hist_MetRatio[i]->Draw();
			if (!data) // rejection power = 1 - met after /met before
			{
				int bin;
				for (unsigned j = 1; j<=hist_MetRatio[i]->GetNbinsX(); ++j)	
				{
					double val = hist_MetRatio[i]->GetBinContent(j);
					double err = hist_MetRatio[i]->GetBinError(j);
					hist_MetRatio[i]->SetBinContent(j, 1 - val);
					hist_MetRatio[i]->SetBinError(j, err);
				}
			}
			
			delete f[i];
		}
	}

//return;

	TF1 *funy1 = new TF1("y1","1",0,500);
	funy1->SetLineStyle(1);
	funy1->SetLineWidth(1);
	funy1->SetLineColor(46);
	TF1 *fitFun1=new TF1("fit1",Func1,0,500.0,5);
	fitFun1->SetLineStyle(2);
	fitFun1->SetLineWidth(1);
	fitFun1->SetLineColor(6);
	//fitFun1->SetParameter(0,1.0);				//plataeu
	fitFun1->FixParameter(0,1.0);				//plataeu
	fitFun1->SetParameter(1,10);				//turnon point
	fitFun1->SetParameter(2,30);				//turnon width
	fitFun1->SetParameter(3,2);
	fitFun1->SetParameter(4,1);
	float cut1 = 0;
	float cut2 = 500;
	
	gStyle->SetOptStat("");
	
	gStyle->SetCanvasColor (10);
	gStyle->SetCanvasBorderSize (0);
	gStyle->SetCanvasBorderMode (0);
				 
	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (0);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);
	//new TCanvas();
	TLegend *leg = new TLegend (0.85,0.6,0.99,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);
	leg->AddEntry(hist_MetRatio[0],"#slash{E}_{T}-Sig>2");
	leg->AddEntry(hist_MetRatio[1],"#slash{E}_{T}-Sig>3");
	leg->AddEntry(hist_MetRatio[2],"#slash{E}_{T}-Sig>4");
	leg->AddEntry(hist_MetRatio[3],"#slash{E}_{T}-Sig>5");
	leg->AddEntry(hist_MetRatio[4],"#slash{E}_{T}-Sig>6");
	leg->AddEntry(hist_MetRatio[5],"#slash{E}_{T}-Sig>7");

	for (int i=0; i<6; i++)
	{
		std::stringstream title, gif;
		if (data) //real met eff.
		{
			title << "MetModel Efficiency (MetSig > " << i + 2 << ")";
			//title << "MetModel Efficiency - W#rightarrowe#nu MC (all datasets)";
			gif << "Eff_MetSig" << i + 2 << ".gif";
		} else {  //qcd rejection power
			title << "MetModel Rejection Power - #gamma MC (jqdfh)";
			gif << "Rej_MetSig" << i + 2 << ".gif";
		}
		hist_MetRatio[i]->SetTitle(title.str().c_str());
		//hist_MetRatio[i]->SetTitle("");
		hist_MetRatio[i]->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");
		if (data)
		{
			hist_MetRatio[i]->GetYaxis()->SetTitle("Efficiency = #slash{E}_{T} after #slash{E}_{T}-sig cut / #slash{E}_{T}");
		} else {
			hist_MetRatio[i]->GetYaxis()->SetTitle("Rejection Power =1 - #frac{#slash{E}_{T} after #slash{E}_{T}-sig cut}{#slash{E}_{T}}");
		}
		//hist_MetRatio[i]->SetLineColor(kBlue);
		//hist_MetRatio[i]->SetMarkerColor(kBlue);
		int color;
		if (i==0) color = 2;
		else if (i==1) color = 4;
		else if (i==2) color = 30;
		else if (i==3) color = 3;
		else if (i==4) color = 6;
		else if (i==5) color = 12;
		hist_MetRatio[i]->SetLineColor(color);
		hist_MetRatio[i]->SetMarkerColor(color);
		hist_MetRatio[i]->SetMarkerSize(1);
		hist_MetRatio[i]->SetMinimum(0);
		hist_MetRatio[i]->GetXaxis()->SetLimits(0,fGlobalMaxX);
		
		new TCanvas();
		hist_MetRatio[i]->Draw("EP");
		//if (i==0) hist_MetRatio[i]->Draw("E");
		//else hist_MetRatio[i]->Draw("sameE");
		funy1->Draw("sameC");
		//hist_MetRatio[i]->Fit(fitFun1,"","",cut1,cut2);  //"E"=perform better error estimate using Minos // 
		//gPad->Print(gif.str().c_str());
		//break;
	}

	leg->Draw();
}
