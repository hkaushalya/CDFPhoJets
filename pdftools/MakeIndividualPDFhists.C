////////////////////////////////////////////////////////
//This to make individual PDF systematic plots afteri //
// generating hists from pdf_syst.C (run_pdf_syst.csh)//
// 10-29-2009- sam                                    //
////////////////////////////////////////////////////////

/* $Id: MakeIndividualPDFhists.C,v 1.1 2009/11/19 22:46:27 samantha Exp $
 * $Log: MakeIndividualPDFhists.C,v $
 * Revision 1.1  2009/11/19 22:46:27  samantha
 * Created this to make individual pairs of PDF plots.
 *
 */


#include <iostream>
#include <sstream>
#include "TStyle.h"
#include "TCanvas.h" 
#include "TH1.h" 
#include "TFile.h" 
#include "TAxis.h" 
#include "TLegend.h" 
#include <cmath>
#include <vector>


//-----------------------------------------------------------------------------
// effects: make a new binning with variable bin sizes, and the given
//   cutoff points between bin sizes
// returns: the binning (by bin edges)
// guarantee: strong
//// throws: std::bad_alloc
//-----------------------------------------------------------------------------
std::vector<float> 
GetVarBinVector (float down, float up5, float up10, float up20, float up50, 
					float width1, float width2, float width3, float width4 ,bool debug)
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
// effects: turn this histogram into a histogram with variable bin
//   size
// returns: the new histogram
// guarantee: strong
// throws: std::bad_alloc
// requires: input != NULL
// requires: input->GetDimension() == 1
// requires: histogram bin edges must match
// postcondition: result != NULL
//-----------------------------------------------------------------------------
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

	for (unsigned bin = 1; bin != unsigned (hist->GetNbinsX()); ++ bin)
	{
		const float width = hist->GetXaxis()->GetBinWidth (bin);
		hist->SetBinContent (bin, hist->GetBinContent (bin) / width);
		hist->SetBinError (bin, hist->GetBinError (bin) / width);
	};
};

//-----------------------------------------------------------------------------
// make a variable binned hist with given specification
// TESTED THIS METHOD AND ITS SUB-METHODS - sam
//-----------------------------------------------------------------------------
TH1* MakeVariableBins (TH1 *hist, const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4 , bool debug=false)
{
	if (debug)
	{
		std::cout << hist->GetName() << " hist bin size=" << hist->GetBinWidth(1) 
					<< std::endl;
	}

  	std::auto_ptr<TH1> result = make_var_bin ( 
											GetVarBinVector(xmin, xpoint1, xpoint2, 
															xpoint3, xpoint4, width1, 
															width2, width3, width4, debug) , 
											hist, debug);
  	norm_bins (result.get());
  	return result.release ();
};

void MakeIndividualPDFhists(const std::string name)
{
	const int nfiles = 6;
	
	TFile* f[nfiles];
	f[0] = new TFile("~/RESULTS/10222009_PDFsystFullMCDataset/10272009_Result3/PDFsyst1.root");
	f[1] = new TFile("~/RESULTS/10222009_PDFsystFullMCDataset/10272009_Result3/PDFsyst2.root");
	f[2] = new TFile("~/RESULTS/10222009_PDFsystFullMCDataset/10272009_Result3/PDFsyst3.root");
	f[3] = new TFile("~/RESULTS/10222009_PDFsystFullMCDataset/10272009_Result3/PDFsyst4.root");
	f[4] = new TFile("~/RESULTS/10222009_PDFsystFullMCDataset/10272009_Result3/PDFsyst5.root");
	f[5] = new TFile("~/RESULTS/10222009_PDFsystFullMCDataset/10272009_Result3/PDFsyst6.root");


	assert (f[0] !=NULL && "FILE 1 NOT FOUND!");
	assert (f[1] !=NULL && "FILE 2 NOT FOUND!");
	assert (f[2] !=NULL && "FILE 3 NOT FOUND!");
	assert (f[3] !=NULL && "FILE 4 NOT FOUND!");
	assert (f[4] !=NULL && "FILE 5 NOT FOUND!");
	assert (f[5] !=NULL && "FILE 6 NOT FOUND!");

	TH1* h[41];
	
	float xmin = 0;
	float xpoint1 = 200;
	float xpoint2 = 250;
	float xpoint3 = 300;
	float xpoint4 = 650;
	float width1 = 10;
	float width2 = 10;
	float width3 = 50;
	float width4 = 250;
	
	//0th hist is the reference. other 40 are the variants.
	for (int i=0; i < 41; ++i)
	{
		std::stringstream hist;
		hist << "Ana/PhoJetsTemp/Hist/1Jet/PDFSyst/PDFSystPara"<< i << "_" << name;
		std::cout << "looking for hist " << hist.str() << std::endl;
		h[i] = dynamic_cast<TH1*> (f[0]->Get(hist.str().c_str()));
		assert (h[i] != NULL && "hist is null");
		//Sumw2() for these hists are called at creation
		
		//add the other 2 jobs resutls
		for (int nf=1;nf< nfiles;++nf)
		{
			TH1* temp = dynamic_cast<TH1*> (f[nf]->Get(hist.str().c_str()));
			assert (temp != NULL && "hist is null");
			h[i]->Add(temp);
		}
		

		h[i] = (TH1F*) MakeVariableBins (h[i], xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		if (i>0)
		{
			h[i]->Scale(h[0]->Integral()/(1. * h[i]->Integral()));
			h[i]->Divide(h[0]);
			int marker = 22, color=2;
			if (i%2)
			{
				marker = 23;
				color = 4;
			}
				
			h[i]->SetMarkerColor(color);
			h[i]->SetLineColor(color);
			h[i]->SetMarkerStyle(marker);
			h[i]->SetMarkerSize(1);
		}
	}


			gStyle->SetOptStat(0);
		new TCanvas();
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetTickx();
		gPad->SetTicky();



	//make a plot for each PDF parameter wrt to reference
	for (int ipdf = 1; ipdf <41; ipdf=ipdf+2)
	{
		//h[ipdf]->Scale(h[0]->Integral()/(1. * h[ipdf]->Integral()));
		//h[ipdf+1]->Scale(h[0]->Integral()/(1. * h[ipdf+1]->Integral()));
		std::string xtitle, ytitle;
		if (name == "PhoEt") xtitle = "E_{T}^{#gamma}";
		ytitle = "Ratio wrt base (Parameter 0)";
		
		std::stringstream title, le1,le2, file;
		int jpdf = ipdf+1;
		title << "PDF parameter " << ipdf << " & "<< jpdf << "(normalized to base and ratio wrt to base)"
				<< ";" << xtitle << ";" << ytitle;
		le1 << "Parameter " << ipdf;
		le2 << "Parameter " << jpdf;
		file << "PDFpara_"<< ipdf << "_" << jpdf << ".eps";

		
		//h[ipdf]->SetMinimum(0.9);
		//h[ipdf]->SetMaximum(1.1);
		//h[ipdf]->SetTitle(title.str().c_str());
		h[ipdf]->GetXaxis()->CenterTitle(1);
		h[ipdf]->GetYaxis()->CenterTitle(1);

		if (ipdf==1) 
		{
			h[ipdf]->Draw();
			h[ipdf+1]->Draw("SAME");
		} else 
		{
			h[ipdf]->Draw("SAME");
			h[ipdf+1]->Draw("SAME");
		}

		TLegend *leg = new TLegend (0.7,0.8,0.9,0.9);
		leg->SetTextFont(42);
		leg->SetTextSize(0.03);


		leg->AddEntry(h[ipdf],le1.str().c_str());
		leg->AddEntry(h[ipdf+1],le2.str().c_str());
		//leg->Draw();
		//gPad->Print(file.str().c_str());
		
	}



	return;

	
	for (int i=1; i < 40; ++i)
	{
		if (i==1)
		{
			h[i]->Draw("P");
			h[i]->SetTitle("PDF: 40 variants wrt to CTEQ6;E_{T}^{#gamma};Ratio");
			//h[i]->SetTitle("PDF: 40 variants wrt to CTEQ6;E_{T}^{jet};Ratio");
		}
		else h[i]->Draw("P same");
	}

	//no make a plot with max/min values
	//exclude the 0th hist which is the reference
	
	for (int bin=0; bin<= h[1]->GetNbinsX(); ++bin)
	{
		if (h[1]->GetBinContent(bin))
		{
			float max=1,min=1;
			float maxer=0,miner=0;
			for (int nh=2; nh<41; ++nh)
			{
				if (h[nh]->GetBinContent(bin) > max )
				{
					max = h[nh]->GetBinContent(bin);
					maxer = h[nh]->GetBinError(bin);
				}
				if (h[nh]->GetBinContent(bin) < min )
				{
					min = h[nh]->GetBinContent(bin);
					miner = h[nh]->GetBinError(bin);
				}
			}
			
			h[1]->SetBinContent(bin,max);
			h[1]->SetBinError(bin,maxer);
			h[2]->SetBinContent(bin,min);
			h[2]->SetBinError(bin,miner);
		}

	}

	new TCanvas();
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetTickx(2);
	gPad->SetTicky(2);
	h[1]->SetTitle("Max/Min of 40 PDF variants");
	h[1]->Draw("PE1");
	h[2]->Draw("SAME PE1");
	
	

return;

	TH1* ph = dynamic_cast<TH1*> (f[0]->Get("Ana/PhoJetsTemp/Hist/1Jet/Photon/EtCorr"));
	assert (ph != NULL && " photon hist is null");
	for (int nf=1;nf< nfiles;++nf)
	{
		TH1* temp = dynamic_cast<TH1*> (f[nf]->Get("Ana/PhoJetsTemp/Hist/1Jet/Photon/EtCorr"));
		assert (temp != NULL && "photon jobs hist is null");
		temp->Sumw2();
		ph->Add(temp);
	}

	ph = (TH1F*) MakeVariableBins (ph, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);

	new TCanvas();
	ph->Draw();
	return;
	
	TH1* final = dynamic_cast<TH1*> (ph->Clone("copy"));
	assert (final != NULL && " photon hist is null");
	for (int bin =0; bin <= ph->GetNbinsX()+1; ++bin)
	{
		final->SetBinError(bin,0);
		final->SetBinContent(bin,0);
	}
		
	final->Print();
		
	final->SetDirectory(0);
	//ph->Sumw2();
	//final->Sumw2();
	
	//ph->Rebin(rebin);
	//final->Rebin(rebin);

	int count =0;
	for (int bin =1; bin < ph->GetNbinsX(); ++bin)
	{
		float err2=0;
		for (int i=0; i<40; ++i)
		{
			err2 += pow( (h[i]->GetBinContent(bin) - 1),2);

			if (ph->GetBinContent(bin) && count<10)
			{
				std::cout << err2 << "\t" << h[i]->GetBinContent(bin)  << "\t" 
						<< h[i]->GetBinContent(bin) - 1 << "\t" 
						<< pow( (h[i]->GetBinContent(bin) - 1),2) << std::endl;
				count++;
			}

		}
		std::cout << "err = " << sqrt(err2) << "\t" << sqrt(err2) * ph->GetBinContent(bin) << std::endl;
		//final->SetBinError(bin, sqrt(err2) * ph->GetBinContent(bin));
		final->SetBinContent(bin, sqrt(err2) * ph->GetBinContent(bin));
		//final->SetBinContent(bin,0);

	}

	final->SetLineColor(kBlue);
	final->SetMarkerColor(kBlue);
	//final->SetFillColor(kBlue);
	final->Draw();

	new TCanvas();
	ph->Draw();


}

void MakeIndividualPDFhists()
{

	MakeIndividualPDFhists("PhoEt");
	//MakeIndividualPDFhists("LeadJetEt");

}
