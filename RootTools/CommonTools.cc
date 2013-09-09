/*******************************************************************
*  This has a bunch of utilities that is used under Stntuple
*  framework and when I make the final hitsts.
*  Author: Sam Hewamanage <samantha@fnal.gov>
********************************************************************
*	$Id: CommonTools.cc,v 1.4 2012/10/29 19:41:08 samantha Exp $
*	$Log: CommonTools.cc,v $
*	Revision 1.4  2012/10/29 19:41:08  samantha
*	Final commit.
*	
*	Revision 1.3  2010/03/18 15:34:45  samantha
*	ADDED: 1. Two overloaded GetMaxBinContent() method to find the contents of the
*	maximum bin in a hist or a set of hists.
*		2. class StringVector, which can hold a set of strings to which can be
*		printed at anytime during a job. Same info can be written to a canvas
*		and saved to the root file at the end of the job. I created this so  do
*		not have to search for the log file for a job. The root file will have
*		all the info I need about the job.
*	
*	Revision 1.2  2010/02/13 22:11:41  samantha
*	MODIFIED: FindNonZeroXrange() is modified to inly search the visible
*	region of the x-axis (i.e. exclude underflow and overflow bins.)
*	
*	Revision 1.1  2009/12/15 23:23:54  samantha
*	A collection of utitlity functions, color couts, Draw hists, rebin, etc.
*	
*
********************************************************************/

#include "CommonTools.hh"
#include <iostream>
#include "TStyle.h"
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "TLatex.h"

//to compile this first load Stuple_cc.so

//=================================================================//
//dump a bin info
//=================================================================//
void DumpBin(const TH1* hist, const int iBin,const std::string sMsg)
{
	assert (hist != NULL && "DumpBin::Hist is null!");
	assert (hist->GetDimension() == 1 && "CommonTools::DumpBin:: hist is not 1-D");
	assert ( (iBin>=1 && iBin<= hist->GetNbinsX()) && "DumpBin::Bin requested is out of range");
	
	if (sMsg.length()>0) std::cout << sMsg << std::endl;
	std::cout << "\t\t" << __FUNCTION__ << "::bin[loedge] : [val/err] = "
		<< std::setw(3) << iBin << "[" << hist->GetXaxis()->GetBinLowEdge(iBin) << "]" << std::setw(20)
		<< hist->GetBinContent(iBin) << " / " << hist->GetBinError(iBin) << std::endl;
}

//=================================================================//
void DrawClone(const TH1* h, const bool log)
{
	assert (h != NULL && "hist is null");
	assert (h->GetDimension() == 1 && "CommonTools::DrawClone:: hist is not 1-D");
	
	TH1* hist = dynamic_cast<TH1*> (h->Clone("copy"));
	new TCanvas();
	gStyle->SetOptStat("neou");
	gPad->SetGridx();
	gPad->SetGridy();
	if (log) gPad->SetLogy();
	if (hist->GetLineColor() == 5) hist->SetLineColor(46);
	hist->SetStats(1);
	hist->DrawClone();
}

//=================================================================//
// divide bins by bin width
//=================================================================//
void NormalizeBinWidth(TH1* hist)
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

//=================================================================//
//dumps hist bins including underflow and overflow
//=================================================================//
void DumpHist(const TH1* hist, const std::string sMsg)
{
	assert (hist != NULL && "requirement failed"); //spec
	assert (hist->GetDimension() == 1 && "CommonTools::DumpHist:: hist is not 1-D");

	if (sMsg.length()>0) std::cout << sMsg << std::endl;
	std::cout << std::setw(4) << "Bin" << std::setw(11) << "Low Edge" 
		<< std::setw(15) << "Bin Val"
		<< std::setw(15) << "Bin Err"
		<< std::endl;

	for (int bin = 0; bin <= hist->GetNbinsX()+1; ++bin)
	{
		if (hist->GetBinContent(bin)>0 || hist->GetBinError(bin)>0)
		{
			std::cout << "[" << std::setw(3) << bin << "]" 
				<< std::setw(9) << hist->GetXaxis()->GetBinLowEdge(bin)
				<< std::setw(15) << hist->GetBinContent(bin) 
				<< std::setw(15) << hist->GetBinError(bin) << std::endl;
		}
	}

}
//=================================================================//
//finds the bin number corresponding to a value
//=================================================================//
int GetBin(const TH1* hist, const double val)
{
	assert (hist != NULL && "GetBin::Hist is null!");
	for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
	{
		if (val > hist->GetXaxis()->GetBinLowEdge(bin) 
			  && val <= hist->GetXaxis()->GetBinUpEdge(bin)) return bin;
	}
	std::cout << __FUNCTION__ << "::WARNING! Did not find a bin for Value =" << val
		<< ". Hist bounds are [" << hist->GetXaxis()->GetBinLowEdge(1)
		<< ", " << hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX()) 
		<< "] ! returning -1" << std::endl;
	return -1;
}


//=================================================================//
// find the nozero x range two two 1-D hists
//=================================================================//
void FindCommonNonZeroXrange(const TH1* h1, const TH1* h2, float& xmin, float& xmax)
{
	assert (h1 != NULL && "requirement failed"); //spec
	assert (h2 != NULL && "requirement failed"); //spec
	assert (h1->GetDimension() == 1 && "CommonTools::FindNonZeroXrange:: hist1 is not 1-D");
	assert (h2->GetDimension() == 1 && "CommonTools::FindNonZeroXrange:: hist2 is not 1-D");
	for (int bin=0; bin<=h1->GetNbinsX(); ++bin)
	{
		if (h1->GetBinContent(bin)>0 || h2->GetBinContent(bin)>0)
		{
			xmin = h1->GetXaxis()->GetBinLowEdge(bin);
			break;
		}
	}
	for (int bin=h1->GetNbinsX(); bin>=1; --bin)
	{
		if (h1->GetBinContent(bin)>0 || h2->GetBinContent(bin)>0)
		{
			xmax = h1->GetXaxis()->GetBinUpEdge(bin);
			break;
		}
	}
}
//=================================================================//
// find the nozero x range for a 1-D hist in visible region
//=================================================================//
void FindNonZeroXrange(const TH1* h1, int& xminbin, int& xmaxbin)
{
	assert (h1 != NULL && "CommonTools::FindNonZeroXrange:: hist1 is NULL!");
	assert (h1->GetDimension() == 1 && "CommonTools::FindNonZeroXrange:: hist1 is not 1-D!");
	for (int bin=1; bin<=h1->GetNbinsX(); ++bin)
	{
		if (h1->GetBinContent(bin)>0)
		{
			xminbin = bin;
			break;
		}
	}
	for (int bin=h1->GetNbinsX()+1; bin>=0; --bin)
	{
		if (h1->GetBinContent(bin)>0)
		{
			xmaxbin = bin;
			break;
		}
	}
}

//-----------------------------------------------------------------------------
// effects: make a new binning with variable bin sizes, and the given
//   cutoff points between bin sizes
// returns: the binning (by bin edges)
// guarantee: strong
//// throws: std::bad_alloc
//-----------------------------------------------------------------------------
std::vector<float> 
GetVarBinVector (float down, float up5, float up10, float up20, float up50, 
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
// make a variable binned hist with given specification (overloaded method)
//-----------------------------------------------------------------------------
TH1F* MakeVariableBinHist (const std::string name, const std::string title,
							const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4)
{

	std::vector<float> vBins = GetVarBinVector(xmin, xpoint1, xpoint2, 
															xpoint3, xpoint4, width1, 
															width2, width3, width4);
	return new TH1F(name.c_str(), title.c_str(),vBins.size()-1, &vBins[0]);
}

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
std::auto_ptr<TH1> FillVarBinHist (const std::vector<float>& bins, TH1 *input)
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
// make variable binned hist from a given hist. this method is overloaded
//-----------------------------------------------------------------------------
TH1* MakeVariableBinHist (TH1 *hist, const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4)
{
	assert (hist != NULL && "CommonTools::MakeVariableBinHist:: hist is NULL!");
	assert (hist->GetDimension() == 1 && "CommonTools::MakeVariableBinHist:: hist is not 1-D");
	
  	std::auto_ptr<TH1> result = FillVarBinHist ( 
											GetVarBinVector(xmin, xpoint1, xpoint2, 
															xpoint3, xpoint4, width1, 
															width2, width3, width4)
											, hist);
  	NormalizeBinWidth(result.get());
  	return result.release ();
};


void WilsonInterval(float k, float n, float &center, float &error)
{
	// See https://hep.baylor.edu/elog/mfrank/125 or search
	// "Wilson Score Interval" for a derivation of the below.
	// We pick z = 1 for a 68.3% confidence interval

	float p;

	if (n > 0)
	{
		p = k/n;
		center = ( p + (1/(2*n)) ) / (1 + (1/n));
		error = sqrt( (1/n) * (p*(1-p) + (1/(4*n))) ) / (1+(1/n));
	} else {
		center = 0;
		error = 0;
	}
}



void WilsonInterval(TH1* &nom, TH1* denom)
{
	for (int bin = 1; bin <= nom->GetNbinsX(); ++bin)
	{
		float k = nom->GetBinContent(bin);
		float n = denom->GetBinContent(bin);
		float center, error;

		WilsonInterval(k, n, center, error);

		nom->SetBinContent(bin, center);
		nom->SetBinError(bin, error);
	}
}

Clock::Clock()
{
	start = -1;
	end = -1;
};
Clock::Clock(const std::string n)
{
	name = n;
	Clock();
};
void Clock::Start() 
{
	std::cout << "in start " << std::endl;
	start = time_t(NULL);	//CINT seem to return 0 for this! 
	std::cout << "start time = " << start << std::endl;
}
void Clock::Stop() 
{ 
	std::cout << "in stop " << std::endl;
	end = time_t(NULL);  
	std::cout << "end  time = " << end << std::endl;
}

void Clock::ShowEndTimeHours() const 
{ 
	if (end > 0)
		std::cout << name << ": Total Time " <<  (end - start)/3600. << " (hrs)" << std::endl; 
	else 
		std::cout << "Clock has not stopped!" << start << ", " << end << std::endl;

	std::cout << name << ": Total Time " <<  (end - start)/3600. << " (hrs)" << std::endl; 
}

//find the value of the bin with the largest number of entries
//in a given hist
double GetMaxBinContent(const TH1* hist)
{
	assert (hist != NULL && "CommonTools::GetMaxBinContent()::Given hist is null!");
	assert (hist->GetDimension() == 1 && "CommonTools::GetMaxBinContent()::Hist must be 1-D!");
	return hist->GetBinContent(hist->GetMaximumBin());
}
//find the value of the bin with the largest number of entries
//in a set of given hists
double GetMaxBinContent(const std::vector<TH1*> hist)
{
	double max = 0;
	for (unsigned int i=0; i< hist.size(); ++i)
	{
		max = std::max(max, GetMaxBinContent(hist.at(i)));
	}
	return max;
}

StringVector::StringVector()
{
	c = 0;
};
StringVector::~StringVector()
{
	if (c != 0 || c != NULL) delete c;
};

TCanvas* StringVector::writeToCanvas()
{
	if (strVec.size())
	{
		if (strVec.size()>50) std::cout << __FUNCTION__ 
			<< ": Warning! Too many lines to print. Some may not print correctly. Think of creating a new StringVector..." << std::endl;
		if (c == 0 || c == NULL) c = new TCanvas();
		//const float h = c->GetWindowHeight();
	//	const float w = c->GetWindowWidth();

		float lineheight = 0.02;
		float ypos = 0.95;
		//float ypos = 0.1;
		TLatex l;
		l.SetTextSize(0.02);
		l.SetTextFont(80);
		l.SetTextAlign(12);
		
		for (std::vector<std::string>::const_iterator it = strVec.begin(); it != strVec.end(); ++it)
		{
			l.DrawLatex(0.03,ypos, (*it).c_str());
			ypos -= lineheight;
		}
		return dynamic_cast<TCanvas*> (c->Clone());
	} else 
	{
		std::cout << __FUNCTION__ << ":Vector is empty! Nothing to write.! returning null pointer!" << std::endl;
		return c;
	}


};

TTree* GetStupleTree(const std::string sTreeName, Stuple& stuple)
{

		TTree *tree = new TTree(sTreeName.c_str(), sTreeName.c_str());
		tree->Branch("evt_McFlag"   	 , &stuple.evt_McFlag   	, "evt_McFlag/I");
		tree->Branch("evt_RunNumber"   , &stuple.evt_RunNumber   , "evt_RunNumber/i");
		tree->Branch("evt_EventNumber" , &stuple.evt_EventNumber , "ect_EventNumber/i");

		tree->Branch("tri_pho25iso" , &stuple.tri_pho25iso , "tri_pho25iso/I");
		tree->Branch("tri_pho50"    , &stuple.tri_pho50    , "tri_pho50/I");
		tree->Branch("tri_pho70"    , &stuple.tri_pho70    , "tri_pho70/I");
	
		tree->Branch("pho_num"     , &stuple.pho_num      , "pho_num/i");
		tree->Branch("pho_Ntight"  , &stuple.pho_Ntight   , "pho_Ntight/i");
		tree->Branch("pho_Nloose"  , &stuple.pho_Nloose   , "pho_Nloose/i");
		tree->Branch("pho_Index"   , &stuple.pho_Index    , "pho_Index[pho_num]/I");
		tree->Branch("pho_PhoBlockIndex"   , &stuple.pho_PhoBlockIndex    , "pho_PhoBlockIndex[pho_num]/I");
		tree->Branch("pho_Etc"     , &stuple.pho_Etc      , "pho_Etc[pho_num]/F");
		tree->Branch("pho_E"       , &stuple.pho_E        , "pho_E[pho_num]/F");
		tree->Branch("pho_Px"      , &stuple.pho_Px       , "pho_Px[pho_num]/F");
		tree->Branch("pho_Py"      , &stuple.pho_Py       , "pho_Py[pho_num]/F");
		tree->Branch("pho_Pz"      , &stuple.pho_Pz       , "pho_Pz[pho_num]/F");
		tree->Branch("pho_Detector", &stuple.pho_Detector , "pho_Detector[pho_num]/I");
		tree->Branch("pho_DetEta"  , &stuple.pho_DetEta   , "pho_DetEta[pho_num]/F");
		tree->Branch("pho_DetPhi"  , &stuple.pho_DetPhi   , "pho_DetPhi[pho_num]/F");
		tree->Branch("pho_XCes"    , &stuple.pho_XCes     , "pho_XCes[pho_num]/F");
		tree->Branch("pho_ZCes"    , &stuple.pho_ZCes     , "pho_ZCes[pho_num]/F");
		tree->Branch("pho_HadEm"   , &stuple.pho_HadEm    , "pho_HadEm[pho_num]/F");
		tree->Branch("pho_Chi2Mean", &stuple.pho_Chi2Mean , "pho_Chi2Mean[pho_num]/F");
		tree->Branch("pho_N3d"     , &stuple.pho_N3d      , "pho_N3d[pho_num]/I");
		tree->Branch("pho_Iso4"    , &stuple.pho_Iso4     , "pho_Iso4[pho_num]/F");
		tree->Branch("pho_TrkPt"   , &stuple.pho_TrkPt    , "pho_TrkPt[pho_num]/F");
		tree->Branch("pho_TrkIso"  , &stuple.pho_TrkIso   , "pho_TrkIso[pho_num]/F");
		tree->Branch("pho_CesWireE2"      , &stuple.pho_CesWireE2      , "pho_CesWireE2[pho_num]/F");
		tree->Branch("pho_CesStripE2"     , &stuple.pho_CesStripE2     , "pho_CesStripE2[pho_num]/F");
		tree->Branch("pho_PhiWedge"       , &stuple.pho_PhiWedge       , "pho_PhiWedge[pho_num]/I");
		tree->Branch("pho_NMuonStubs"     , &stuple.pho_NMuonStubs     , "pho_NMuonStubs[pho_num]/I");
		tree->Branch("pho_EmTime"         , &stuple.pho_EmTime         , "pho_EmTime[pho_num]/F");
		tree->Branch("pho_TightId"        , &stuple.pho_TightId        , "pho_TightId[pho_num]/I");
		tree->Branch("pho_LooseId"        , &stuple.pho_LooseId        , "pho_LooseId[pho_num]/I");
		tree->Branch("pho_PhoenixId"      , &stuple.pho_PhoenixId      , "pho_PhoenixId[pho_num]/I");
		tree->Branch("pho_Halo_seedWedge" , &stuple.pho_Halo_seedWedge , "pho_Halo_seedWedge[pho_num]/I");
		tree->Branch("pho_Halo_eastNhad"  , &stuple.pho_Halo_eastNhad  , "pho_Halo_eastNhad[pho_num]/I");
		tree->Branch("pho_Halo_westNhad"  , &stuple.pho_Halo_westNhad  , "pho_Halo_westNhad[pho_num]/I");
		tree->Branch("pho_matchJetIndex"  , &stuple.pho_matchJetIndex  , "pho_matchJetIndex[pho_num]/I");
		tree->Branch("pho_CprWgt"  , &stuple.pho_CprWgt  , "pho_CprWgt[pho_num]/F");
		tree->Branch("pho_CprSys1"  , &stuple.pho_CprSys1  , "pho_CprSys1[pho_num]/F");
		tree->Branch("pho_CprSys2"  , &stuple.pho_CprSys2  , "pho_CprSys2[pho_num]/F");
		tree->Branch("pho_CprSys3"  , &stuple.pho_CprSys3  , "pho_CprSys3[pho_num]/F");
		tree->Branch("pho_CprSys4"  , &stuple.pho_CprSys4  , "pho_CprSys4[pho_num]/F");
		tree->Branch("pho_CprSys5"  , &stuple.pho_CprSys5  , "pho_CprSys5[pho_num]/F");
		tree->Branch("pho_CprSys6"  , &stuple.pho_CprSys6  , "pho_CprSys6[pho_num]/F");
		tree->Branch("pho_CprSys7"  , &stuple.pho_CprSys7  , "pho_CprSys7[pho_num]/F");
		tree->Branch("pho_CprSys8"  , &stuple.pho_CprSys8  , "pho_CprSys8[pho_num]/F");
		

			// Pho EM E uncertainty up 1%
		tree->Branch("pho_up_num"     , &stuple.pho_up_num      , "pho_up_num/i");
		tree->Branch("pho_up_Ntight"  , &stuple.pho_up_Ntight   , "pho_up_Ntight/i");
		tree->Branch("pho_up_Nloose"  , &stuple.pho_up_Nloose   , "pho_up_Nloose/i");
		tree->Branch("pho_up_Index"   , &stuple.pho_up_Index    , "pho_up_Index[pho_up_num]/I");
		tree->Branch("pho_up_PhoBlockIndex"   , &stuple.pho_up_PhoBlockIndex    , "pho_up_PhoBlockIndex[pho_up_num]/I");
		tree->Branch("pho_up_Etc"     , &stuple.pho_up_Etc      , "pho_up_Etc[pho_up_num]/F");
		tree->Branch("pho_up_E"       , &stuple.pho_up_E        , "pho_up_E[pho_up_num]/F");
		tree->Branch("pho_up_Px"      , &stuple.pho_up_Px       , "pho_up_Px[pho_up_num]/F");
		tree->Branch("pho_up_Py"      , &stuple.pho_up_Py       , "pho_up_Py[pho_up_num]/F");
		tree->Branch("pho_up_Pz"      , &stuple.pho_up_Pz       , "pho_up_Pz[pho_up_num]/F");
		tree->Branch("pho_up_Detector", &stuple.pho_up_Detector , "pho_up_Detector[pho_up_num]/I");
		tree->Branch("pho_up_DetEta"  , &stuple.pho_up_DetEta   , "pho_up_DetEta[pho_up_num]/F");
		tree->Branch("pho_up_DetPhi"  , &stuple.pho_up_DetPhi   , "pho_up_DetPhi[pho_up_num]/F");
		tree->Branch("pho_up_XCes"    , &stuple.pho_up_XCes     , "pho_up_XCes[pho_up_num]/F");
		tree->Branch("pho_up_ZCes"    , &stuple.pho_up_ZCes     , "pho_up_ZCes[pho_up_num]/F");
		tree->Branch("pho_up_HadEm"   , &stuple.pho_up_HadEm    , "pho_up_HadEm[pho_up_num]/F");
		tree->Branch("pho_up_Chi2Mean", &stuple.pho_up_Chi2Mean , "pho_up_Chi2Mean[pho_up_num]/F");
		tree->Branch("pho_up_N3d"     , &stuple.pho_up_N3d      , "pho_up_N3d[pho_up_num]/I");
		tree->Branch("pho_up_Iso4"    , &stuple.pho_up_Iso4     , "pho_up_Iso4[pho_up_num]/F");
		tree->Branch("pho_up_TrkPt"   , &stuple.pho_up_TrkPt    , "pho_up_TrkPt[pho_up_num]/F");
		tree->Branch("pho_up_TrkIso"  , &stuple.pho_up_TrkIso   , "pho_up_TrkIso[pho_up_num]/F");
		tree->Branch("pho_up_CesWireE2"      , &stuple.pho_up_CesWireE2      , "pho_up_CesWireE2[pho_up_num]/F");
		tree->Branch("pho_up_CesStripE2"     , &stuple.pho_up_CesStripE2     , "pho_up_CesStripE2[pho_up_num]/F");
		tree->Branch("pho_up_PhiWedge"       , &stuple.pho_up_PhiWedge       , "pho_up_PhiWedge[pho_up_num]/I");
		tree->Branch("pho_up_NMuonStubs"     , &stuple.pho_up_NMuonStubs     , "pho_up_NMuonStubs[pho_up_num]/I");
		tree->Branch("pho_up_EmTime"         , &stuple.pho_up_EmTime         , "pho_up_EmTime[pho_up_num]/F");
		tree->Branch("pho_up_TightId"        , &stuple.pho_up_TightId        , "pho_up_TightId[pho_up_num]/I");
		tree->Branch("pho_up_LooseId"        , &stuple.pho_up_LooseId        , "pho_up_LooseId[pho_up_num]/I");
		tree->Branch("pho_up_PhoenixId"      , &stuple.pho_up_PhoenixId      , "pho_up_PhoenixId[pho_up_num]/I");
		tree->Branch("pho_up_Halo_seedWedge" , &stuple.pho_up_Halo_seedWedge , "pho_up_Halo_seedWedge[pho_up_num]/I");
		tree->Branch("pho_up_Halo_eastNhad"  , &stuple.pho_up_Halo_eastNhad  , "pho_up_Halo_eastNhad[pho_up_num]/I");
		tree->Branch("pho_up_Halo_westNhad"  , &stuple.pho_up_Halo_westNhad  , "pho_up_Halo_westNhad[pho_up_num]/I");
		tree->Branch("pho_up_matchJetIndex"  , &stuple.pho_up_matchJetIndex  , "pho_up_matchJetIndex[pho_up_num]/I");
		tree->Branch("pho_up_CprWgt"  , &stuple.pho_up_CprWgt  , "pho_up_CprWgt[pho_num]/F");
		tree->Branch("pho_up_CprSys1"  , &stuple.pho_up_CprSys1  , "pho_up_CprSys1[pho_num]/F");
		tree->Branch("pho_up_CprSys2"  , &stuple.pho_up_CprSys2  , "pho_up_CprSys2[pho_num]/F");
		tree->Branch("pho_up_CprSys3"  , &stuple.pho_up_CprSys3  , "pho_up_CprSys3[pho_num]/F");
		tree->Branch("pho_up_CprSys4"  , &stuple.pho_up_CprSys4  , "pho_up_CprSys4[pho_num]/F");
		tree->Branch("pho_up_CprSys5"  , &stuple.pho_up_CprSys5  , "pho_up_CprSys5[pho_num]/F");
		tree->Branch("pho_up_CprSys6"  , &stuple.pho_up_CprSys6  , "pho_up_CprSys6[pho_num]/F");
		tree->Branch("pho_up_CprSys7"  , &stuple.pho_up_CprSys7  , "pho_up_CprSys7[pho_num]/F");
		tree->Branch("pho_up_CprSys8"  , &stuple.pho_up_CprSys8  , "pho_up_CprSys8[pho_num]/F");
	
			
			// Pho EM E uncertainty down 1%
		tree->Branch("pho_down_num"     , &stuple.pho_down_num      , "pho_down_num/i");
		tree->Branch("pho_down_Ntight"  , &stuple.pho_down_Ntight   , "pho_down_Ntight/i");
		tree->Branch("pho_down_Nloose"  , &stuple.pho_down_Nloose   , "pho_down_Nloose/i");
		tree->Branch("pho_down_Index"   , &stuple.pho_down_Index    , "pho_down_Index[pho_down_num]/I");
		tree->Branch("pho_down_PhoBlockIndex"   , &stuple.pho_down_PhoBlockIndex    , "pho_down_PhoBlockIndex[pho_down_num]/I");
		tree->Branch("pho_down_Etc"     , &stuple.pho_down_Etc      , "pho_down_Etc[pho_down_num]/F");
		tree->Branch("pho_down_E"       , &stuple.pho_down_E        , "pho_down_E[pho_down_num]/F");
		tree->Branch("pho_down_Px"      , &stuple.pho_down_Px       , "pho_down_Px[pho_down_num]/F");
		tree->Branch("pho_down_Py"      , &stuple.pho_down_Py       , "pho_down_Py[pho_down_num]/F");
		tree->Branch("pho_down_Pz"      , &stuple.pho_down_Pz       , "pho_down_Pz[pho_down_num]/F");
		tree->Branch("pho_down_Detector", &stuple.pho_down_Detector , "pho_down_Detector[pho_down_num]/I");
		tree->Branch("pho_down_DetEta"  , &stuple.pho_down_DetEta   , "pho_down_DetEta[pho_down_num]/F");
		tree->Branch("pho_down_DetPhi"  , &stuple.pho_down_DetPhi   , "pho_down_DetPhi[pho_down_num]/F");
		tree->Branch("pho_down_XCes"    , &stuple.pho_down_XCes     , "pho_down_XCes[pho_down_num]/F");
		tree->Branch("pho_down_ZCes"    , &stuple.pho_down_ZCes     , "pho_down_ZCes[pho_down_num]/F");
		tree->Branch("pho_down_HadEm"   , &stuple.pho_down_HadEm    , "pho_down_HadEm[pho_down_num]/F");
		tree->Branch("pho_down_Chi2Mean", &stuple.pho_down_Chi2Mean , "pho_down_Chi2Mean[pho_down_num]/F");
		tree->Branch("pho_down_N3d"     , &stuple.pho_down_N3d      , "pho_down_N3d[pho_down_num]/I");
		tree->Branch("pho_down_Iso4"    , &stuple.pho_down_Iso4     , "pho_down_Iso4[pho_down_num]/F");
		tree->Branch("pho_down_TrkPt"   , &stuple.pho_down_TrkPt    , "pho_down_TrkPt[pho_down_num]/F");
		tree->Branch("pho_down_TrkIso"  , &stuple.pho_down_TrkIso   , "pho_down_TrkIso[pho_down_num]/F");
		tree->Branch("pho_down_CesWireE2"      , &stuple.pho_down_CesWireE2      , "pho_down_CesWireE2[pho_down_num]/F");
		tree->Branch("pho_down_CesStripE2"     , &stuple.pho_down_CesStripE2     , "pho_down_CesStripE2[pho_down_num]/F");
		tree->Branch("pho_down_PhiWedge"       , &stuple.pho_down_PhiWedge       , "pho_down_PhiWedge[pho_down_num]/I");
		tree->Branch("pho_down_NMuonStubs"     , &stuple.pho_down_NMuonStubs     , "pho_down_NMuonStubs[pho_down_num]/I");
		tree->Branch("pho_down_EmTime"         , &stuple.pho_down_EmTime         , "pho_down_EmTime[pho_down_num]/F");
		tree->Branch("pho_down_TightId"        , &stuple.pho_down_TightId        , "pho_down_TightId[pho_down_num]/I");
		tree->Branch("pho_down_LooseId"        , &stuple.pho_down_LooseId        , "pho_down_LooseId[pho_down_num]/I");
		tree->Branch("pho_down_PhoenixId"      , &stuple.pho_down_PhoenixId      , "pho_down_PhoenixId[pho_down_num]/I");
		tree->Branch("pho_down_Halo_seedWedge" , &stuple.pho_down_Halo_seedWedge , "pho_down_Halo_seedWedge[pho_down_num]/I");
		tree->Branch("pho_down_Halo_eastNhad"  , &stuple.pho_down_Halo_eastNhad  , "pho_down_Halo_eastNhad[pho_down_num]/I");
		tree->Branch("pho_down_Halo_westNhad"  , &stuple.pho_down_Halo_westNhad  , "pho_down_Halo_westNhad[pho_down_num]/I");
		tree->Branch("pho_down_matchJetIndex"  , &stuple.pho_down_matchJetIndex  , "pho_down_matchJetIndex[pho_down_num]/I");
		tree->Branch("pho_down_CprWgt"  , &stuple.pho_down_CprWgt  , "pho_down_CprWgt[pho_num]/F");
		tree->Branch("pho_down_CprSys1"  , &stuple.pho_down_CprSys1  , "pho_down_CprSys1[pho_num]/F");
		tree->Branch("pho_down_CprSys2"  , &stuple.pho_down_CprSys2  , "pho_down_CprSys2[pho_num]/F");
		tree->Branch("pho_down_CprSys3"  , &stuple.pho_down_CprSys3  , "pho_down_CprSys3[pho_num]/F");
		tree->Branch("pho_down_CprSys4"  , &stuple.pho_down_CprSys4  , "pho_down_CprSys4[pho_num]/F");
		tree->Branch("pho_down_CprSys5"  , &stuple.pho_down_CprSys5  , "pho_down_CprSys5[pho_num]/F");
		tree->Branch("pho_down_CprSys6"  , &stuple.pho_down_CprSys6  , "pho_down_CprSys6[pho_num]/F");
		tree->Branch("pho_down_CprSys7"  , &stuple.pho_down_CprSys7  , "pho_down_CprSys7[pho_num]/F");
		tree->Branch("pho_down_CprSys8"  , &stuple.pho_down_CprSys8  , "pho_down_CprSys8[pho_num]/F");
	


		//electron collections	

		tree->Branch("ele_num"     , &stuple.ele_num      , "ele_num/i");
		tree->Branch("ele_Ntight"  , &stuple.ele_Ntight   , "ele_Ntight/i");
		tree->Branch("ele_Nloose"  , &stuple.ele_Nloose   , "ele_Nloose/i");
		tree->Branch("ele_Index"   , &stuple.ele_Index    , "ele_Index[ele_num]/I");
		tree->Branch("ele_PhoBlockIndex", &stuple.ele_PhoBlockIndex , "ele_PhoBlockIndex[ele_num]/I");
		tree->Branch("ele_EleBlockIndex", &stuple.ele_EleBlockIndex, "ele_EleBlockIndex[ele_num]/I");
		tree->Branch("ele_Etc"     , &stuple.ele_Etc      , "ele_Etc[ele_num]/F");
		tree->Branch("ele_E"       , &stuple.ele_E        , "ele_E[ele_num]/F");
		tree->Branch("ele_Px"      , &stuple.ele_Px       , "ele_Px[ele_num]/F");
		tree->Branch("ele_Py"      , &stuple.ele_Py       , "ele_Py[ele_num]/F");
		tree->Branch("ele_Pz"      , &stuple.ele_Pz       , "ele_Pz[ele_num]/F");
		tree->Branch("ele_Detector", &stuple.ele_Detector , "ele_Detector[ele_num]/I");
		tree->Branch("ele_DetEta"  , &stuple.ele_DetEta   , "ele_DetEta[ele_num]/F");
		tree->Branch("ele_DetPhi"  , &stuple.ele_DetPhi   , "ele_DetPhi[ele_num]/F");
		tree->Branch("ele_XCes"    , &stuple.ele_XCes     , "ele_XCes[ele_num]/F");
		tree->Branch("ele_ZCes"    , &stuple.ele_ZCes     , "ele_ZCes[ele_num]/F");
		tree->Branch("ele_HadEm"   , &stuple.ele_HadEm    , "ele_HadEm[ele_num]/F");
		tree->Branch("ele_Chi2Mean", &stuple.ele_Chi2Mean , "ele_Chi2Mean[ele_num]/F");
		tree->Branch("ele_N3d"     , &stuple.ele_N3d      , "ele_N3d[ele_num]/I");
		tree->Branch("ele_Iso4"    , &stuple.ele_Iso4     , "ele_Iso4[ele_num]/F");
		tree->Branch("ele_TrkIso"  , &stuple.ele_TrkIso   , "ele_TrkIso[ele_num]/F");
		tree->Branch("ele_CesWireE2"      , &stuple.ele_CesWireE2      , "ele_CesWireE2[ele_num]/F");
		tree->Branch("ele_CesStripE2"     , &stuple.ele_CesStripE2     , "ele_CesStripE2[ele_num]/F");
		tree->Branch("ele_PhiWedge"       , &stuple.ele_PhiWedge       , "ele_PhiWedge[ele_num]/I");
		tree->Branch("ele_NMuonStubs"     , &stuple.ele_NMuonStubs     , "ele_NMuonStubs[ele_num]/I");
		tree->Branch("ele_EmTime"         , &stuple.ele_EmTime         , "ele_EmTime[ele_num]/F");
		tree->Branch("ele_PhoenixId"      , &stuple.ele_PhoenixId      , "ele_PhoenixId[ele_num]/I");
		tree->Branch("ele_Halo_seedWedge" , &stuple.ele_Halo_seedWedge , "ele_Halo_seedWedge[ele_num]/I");
		tree->Branch("ele_Halo_eastNhad"  , &stuple.ele_Halo_eastNhad  , "ele_Halo_eastNhad[ele_num]/I");
		tree->Branch("ele_Halo_westNhad"  , &stuple.ele_Halo_westNhad  , "ele_Halo_westNhad[ele_num]/I");
		tree->Branch("ele_matchJetIndex"  , &stuple.ele_matchJetIndex  , "ele_matchJetIndex[ele_num]/I");
		tree->Branch("ele_Ntracks"      , &stuple.ele_Ntracks      , "ele_Ntracks[ele_num]/I");
		tree->Branch("ele_Emfr"         , &stuple.ele_Emfr         , "ele_Emfr[ele_num]/F");
		tree->Branch("ele_EoverP"       , &stuple.ele_EoverP       , "ele_EoverP[ele_num]/F");
		tree->Branch("ele_TrackPt"      , &stuple.ele_TrackPt      , "ele_TrackPt[ele_num]/F");
		tree->Branch("ele_TrackBcPt"    , &stuple.ele_TrackBcPt    , "ele_TrackBcPt[ele_num]/F");
		tree->Branch("ele_TrackPhi"     , &stuple.ele_TrackPhi     , "ele_TrackPhi[ele_num]/F");
		tree->Branch("ele_Nssl"         , &stuple.ele_Nssl         , "ele_Nssl[ele_num]/I");
		tree->Branch("ele_Nasl"         , &stuple.ele_Nasl         , "ele_Nasl[ele_num]/I");
		tree->Branch("ele_TightId"      , &stuple.ele_TightId      , "ele_TightId[ele_num]/I");
		tree->Branch("ele_LooseId"      , &stuple.ele_LooseId      , "ele_LooseId[ele_num]/I");
		tree->Branch("ele_ConversionId" , &stuple.ele_ConversionId , "ele_ConversionId[ele_num]/I");
		tree->Branch("ele_StdTightId"      , &stuple.ele_StdTightId      , "ele_StdTightId[ele_num]/I");
		tree->Branch("ele_StdLooseId"      , &stuple.ele_StdLooseId      , "ele_StdLooseId[ele_num]/I");


			// Ele EM E uncertainty up 1%
		tree->Branch("ele_up_num"     , &stuple.ele_up_num      , "ele_up_num/i");
		tree->Branch("ele_up_Ntight"  , &stuple.ele_up_Ntight   , "ele_up_Ntight/i");
		tree->Branch("ele_up_Nloose"  , &stuple.ele_up_Nloose   , "ele_up_Nloose/i");
		tree->Branch("ele_up_Index"   , &stuple.ele_up_Index    , "ele_up_Index[ele_up_num]/I");
		tree->Branch("ele_up_PhoBlockIndex", &stuple.ele_up_PhoBlockIndex , "ele_up_PhoBlockIndex[ele_up_num]/I");
		tree->Branch("ele_up_EleBlockIndex", &stuple.ele_up_EleBlockIndex, "ele_up_EleBlockIndex[ele_up_num]/I");
		tree->Branch("ele_up_Etc"     , &stuple.ele_up_Etc      , "ele_up_Etc[ele_up_num]/F");
		tree->Branch("ele_up_E"       , &stuple.ele_up_E        , "ele_up_E[ele_up_num]/F");
		tree->Branch("ele_up_Px"      , &stuple.ele_up_Px       , "ele_up_Px[ele_up_num]/F");
		tree->Branch("ele_up_Py"      , &stuple.ele_up_Py       , "ele_up_Py[ele_up_num]/F");
		tree->Branch("ele_up_Pz"      , &stuple.ele_up_Pz       , "ele_up_Pz[ele_up_num]/F");
		tree->Branch("ele_up_Detector", &stuple.ele_up_Detector , "ele_up_Detector[ele_up_num]/I");
		tree->Branch("ele_up_DetEta"  , &stuple.ele_up_DetEta   , "ele_up_DetEta[ele_up_num]/F");
		tree->Branch("ele_up_DetPhi"  , &stuple.ele_up_DetPhi   , "ele_up_DetPhi[ele_up_num]/F");
		tree->Branch("ele_up_XCes"    , &stuple.ele_up_XCes     , "ele_up_XCes[ele_up_num]/F");
		tree->Branch("ele_up_ZCes"    , &stuple.ele_up_ZCes     , "ele_up_ZCes[ele_up_num]/F");
		tree->Branch("ele_up_HadEm"   , &stuple.ele_up_HadEm    , "ele_up_HadEm[ele_up_num]/F");
		tree->Branch("ele_up_Chi2Mean", &stuple.ele_up_Chi2Mean , "ele_up_Chi2Mean[ele_up_num]/F");
		tree->Branch("ele_up_N3d"     , &stuple.ele_up_N3d      , "ele_up_N3d[ele_up_num]/I");
		tree->Branch("ele_up_Iso4"    , &stuple.ele_up_Iso4     , "ele_up_Iso4[ele_up_num]/F");
		tree->Branch("ele_up_TrkIso"  , &stuple.ele_up_TrkIso   , "ele_up_TrkIso[ele_up_num]/F");
		tree->Branch("ele_up_CesWireE2"      , &stuple.ele_up_CesWireE2      , "ele_up_CesWireE2[ele_up_num]/F");
		tree->Branch("ele_up_CesStripE2"     , &stuple.ele_up_CesStripE2     , "ele_up_CesStripE2[ele_up_num]/F");
		tree->Branch("ele_up_PhiWedge"       , &stuple.ele_up_PhiWedge       , "ele_up_PhiWedge[ele_up_num]/I");
		tree->Branch("ele_up_NMuonStubs"     , &stuple.ele_up_NMuonStubs     , "ele_up_NMuonStubs[ele_up_num]/I");
		tree->Branch("ele_up_EmTime"         , &stuple.ele_up_EmTime         , "ele_up_EmTime[ele_up_num]/F");
		tree->Branch("ele_up_PhoenixId"      , &stuple.ele_up_PhoenixId      , "ele_up_PhoenixId[ele_up_num]/I");
		tree->Branch("ele_up_Halo_seedWedge" , &stuple.ele_up_Halo_seedWedge , "ele_up_Halo_seedWedge[ele_up_num]/I");
		tree->Branch("ele_up_Halo_eastNhad"  , &stuple.ele_up_Halo_eastNhad  , "ele_up_Halo_eastNhad[ele_up_num]/I");
		tree->Branch("ele_up_Halo_westNhad"  , &stuple.ele_up_Halo_westNhad  , "ele_up_Halo_westNhad[ele_up_num]/I");
		tree->Branch("ele_up_matchJetIndex"  , &stuple.ele_up_matchJetIndex  , "ele_up_matchJetIndex[ele_up_num]/I");
		tree->Branch("ele_up_Ntracks"      , &stuple.ele_up_Ntracks      , "ele_up_Ntracks[ele_up_num]/I");
		tree->Branch("ele_up_Emfr"         , &stuple.ele_up_Emfr         , "ele_up_Emfr[ele_up_num]/F");
		tree->Branch("ele_up_EoverP"       , &stuple.ele_up_EoverP       , "ele_up_EoverP[ele_up_num]/F");
		tree->Branch("ele_up_TrackPt"      , &stuple.ele_up_TrackPt      , "ele_up_TrackPt[ele_up_num]/F");
		tree->Branch("ele_up_TrackBcPt"    , &stuple.ele_up_TrackBcPt    , "ele_up_TrackBcPt[ele_up_num]/F");
		tree->Branch("ele_up_TrackPhi"     , &stuple.ele_up_TrackPhi     , "ele_up_TrackPhi[ele_up_num]/F");
		tree->Branch("ele_up_Nssl"         , &stuple.ele_up_Nssl         , "ele_up_Nssl[ele_up_num]/I");
		tree->Branch("ele_up_Nasl"         , &stuple.ele_up_Nasl         , "ele_up_Nasl[ele_up_num]/I");
		tree->Branch("ele_up_TightId"      , &stuple.ele_up_TightId      , "ele_up_TightId[ele_up_num]/I");
		tree->Branch("ele_up_LooseId"      , &stuple.ele_up_LooseId      , "ele_up_LooseId[ele_up_num]/I");
		tree->Branch("ele_up_ConversionId" , &stuple.ele_up_ConversionId , "ele_up_ConversionId[ele_up_num]/I");
		tree->Branch("ele_up_StdTightId"      , &stuple.ele_up_StdTightId      , "ele_up_StdTightId[ele_up_num]/I");
		tree->Branch("ele_up_StdLooseId"      , &stuple.ele_up_StdLooseId      , "ele_up_StdLooseId[ele_up_num]/I");

			
			// Ele EM E uncertainty down 1%
		tree->Branch("ele_down_num"     , &stuple.ele_down_num      , "ele_down_num/i");
		tree->Branch("ele_down_Ntight"  , &stuple.ele_down_Ntight   , "ele_down_Ntight/i");
		tree->Branch("ele_down_Nloose"  , &stuple.ele_down_Nloose   , "ele_down_Nloose/i");
		tree->Branch("ele_down_Index"   , &stuple.ele_down_Index    , "ele_down_Index[ele_down_num]/I");
		tree->Branch("ele_down_PhoBlockIndex", &stuple.ele_down_PhoBlockIndex , "ele_down_PhoBlockIndex[ele_down_num]/I");
		tree->Branch("ele_down_EleBlockIndex", &stuple.ele_down_EleBlockIndex, "ele_down_EleBlockIndex[ele_down_num]/I");
		tree->Branch("ele_down_Etc"     , &stuple.ele_down_Etc      , "ele_down_Etc[ele_down_num]/F");
		tree->Branch("ele_down_E"       , &stuple.ele_down_E        , "ele_down_E[ele_down_num]/F");
		tree->Branch("ele_down_Px"      , &stuple.ele_down_Px       , "ele_down_Px[ele_down_num]/F");
		tree->Branch("ele_down_Py"      , &stuple.ele_down_Py       , "ele_down_Py[ele_down_num]/F");
		tree->Branch("ele_down_Pz"      , &stuple.ele_down_Pz       , "ele_down_Pz[ele_down_num]/F");
		tree->Branch("ele_down_Detector", &stuple.ele_down_Detector , "ele_down_Detector[ele_down_num]/I");
		tree->Branch("ele_down_DetEta"  , &stuple.ele_down_DetEta   , "ele_down_DetEta[ele_down_num]/F");
		tree->Branch("ele_down_DetPhi"  , &stuple.ele_down_DetPhi   , "ele_down_DetPhi[ele_down_num]/F");
		tree->Branch("ele_down_XCes"    , &stuple.ele_down_XCes     , "ele_down_XCes[ele_down_num]/F");
		tree->Branch("ele_down_ZCes"    , &stuple.ele_down_ZCes     , "ele_down_ZCes[ele_down_num]/F");
		tree->Branch("ele_down_HadEm"   , &stuple.ele_down_HadEm    , "ele_down_HadEm[ele_down_num]/F");
		tree->Branch("ele_down_Chi2Mean", &stuple.ele_down_Chi2Mean , "ele_down_Chi2Mean[ele_down_num]/F");
		tree->Branch("ele_down_N3d"     , &stuple.ele_down_N3d      , "ele_down_N3d[ele_down_num]/I");
		tree->Branch("ele_down_Iso4"    , &stuple.ele_down_Iso4     , "ele_down_Iso4[ele_down_num]/F");
		tree->Branch("ele_down_TrkIso"  , &stuple.ele_down_TrkIso   , "ele_down_TrkIso[ele_down_num]/F");
		tree->Branch("ele_down_CesWireE2"      , &stuple.ele_down_CesWireE2      , "ele_down_CesWireE2[ele_down_num]/F");
		tree->Branch("ele_down_CesStripE2"     , &stuple.ele_down_CesStripE2     , "ele_down_CesStripE2[ele_down_num]/F");
		tree->Branch("ele_down_PhiWedge"       , &stuple.ele_down_PhiWedge       , "ele_down_PhiWedge[ele_down_num]/I");
		tree->Branch("ele_down_NMuonStubs"     , &stuple.ele_down_NMuonStubs     , "ele_down_NMuonStubs[ele_down_num]/I");
		tree->Branch("ele_down_EmTime"         , &stuple.ele_down_EmTime         , "ele_down_EmTime[ele_down_num]/F");
		tree->Branch("ele_down_PhoenixId"      , &stuple.ele_down_PhoenixId      , "ele_down_PhoenixId[ele_down_num]/I");
		tree->Branch("ele_down_Halo_seedWedge" , &stuple.ele_down_Halo_seedWedge , "ele_down_Halo_seedWedge[ele_down_num]/I");
		tree->Branch("ele_down_Halo_eastNhad"  , &stuple.ele_down_Halo_eastNhad  , "ele_down_Halo_eastNhad[ele_down_num]/I");
		tree->Branch("ele_down_Halo_westNhad"  , &stuple.ele_down_Halo_westNhad  , "ele_down_Halo_westNhad[ele_down_num]/I");
		tree->Branch("ele_down_matchJetIndex"  , &stuple.ele_down_matchJetIndex  , "ele_down_matchJetIndex[ele_down_num]/I");
		tree->Branch("ele_down_Ntracks"      , &stuple.ele_down_Ntracks      , "ele_down_Ntracks[ele_down_num]/I");
		tree->Branch("ele_down_Emfr"         , &stuple.ele_down_Emfr         , "ele_down_Emfr[ele_down_num]/F");
		tree->Branch("ele_down_EoverP"       , &stuple.ele_down_EoverP       , "ele_down_EoverP[ele_down_num]/F");
		tree->Branch("ele_down_TrackPt"      , &stuple.ele_down_TrackPt      , "ele_down_TrackPt[ele_down_num]/F");
		tree->Branch("ele_down_TrackBcPt"    , &stuple.ele_down_TrackBcPt    , "ele_down_TrackBcPt[ele_down_num]/F");
		tree->Branch("ele_down_TrackPhi"     , &stuple.ele_down_TrackPhi     , "ele_down_TrackPhi[ele_down_num]/F");
		tree->Branch("ele_down_Nssl"         , &stuple.ele_down_Nssl         , "ele_down_Nssl[ele_down_num]/I");
		tree->Branch("ele_down_Nasl"         , &stuple.ele_down_Nasl         , "ele_down_Nasl[ele_down_num]/I");
		tree->Branch("ele_down_TightId"      , &stuple.ele_down_TightId      , "ele_down_TightId[ele_down_num]/I");
		tree->Branch("ele_down_LooseId"      , &stuple.ele_down_LooseId      , "ele_down_LooseId[ele_down_num]/I");
		tree->Branch("ele_down_ConversionId" , &stuple.ele_down_ConversionId , "ele_down_ConversionId[ele_down_num]/I");
		tree->Branch("ele_down_StdTightId"      , &stuple.ele_down_StdTightId      , "ele_down_StdTightId[ele_down_num]/I");
		tree->Branch("ele_down_StdLooseId"      , &stuple.ele_down_StdLooseId      , "ele_down_StdLooseId[ele_down_num]/I");


		// JET COLLECTIONS
		tree->Branch("jet_num"      , &stuple.jet_num      , "jet_num/i");
		tree->Branch("jet_NJet15"   , &stuple.jet_NJet15   , "jet_Njet15/I");
		tree->Branch("jet_Index"    , &stuple.jet_Index    , "jet_Index[jet_num]/I");
		tree->Branch("jet_Pt"       , &stuple.jet_Pt       , "jet_Pt[jet_num]/F");
		tree->Branch("jet_E"        , &stuple.jet_E        , "jet_E[jet_num]/F");
		tree->Branch("jet_Px"       , &stuple.jet_Px       , "jet_Px[jet_num]/F");
		tree->Branch("jet_Py"       , &stuple.jet_Py       , "jet_Py[jet_num]/F");
		tree->Branch("jet_Pz"       , &stuple.jet_Pz       , "jet_Pz[jet_num]/F");
		tree->Branch("jet_DetEta"   , &stuple.jet_DetEta   , "jet_DetEta[jet_num]/F");
		tree->Branch("jet_DetPhi"   , &stuple.jet_DetPhi   , "jet_DetPhi[jet_num]/F");
		tree->Branch("jet_HadEm"    , &stuple.jet_HadEm    , "jet_HadEm[jet_num]/F");
		tree->Branch("jet_Emfr"     , &stuple.jet_Emfr     , "jet_Emfr[jet_num]/F");
		tree->Branch("jet_Ntowers"  , &stuple.jet_Ntowers  , "jet_Ntowers[jet_num]/I");
		tree->Branch("jet_Ntracks"  , &stuple.jet_Ntracks  , "jet_Ntracks[jet_num]/I");
		tree->Branch("jet_SeedIPhi" , &stuple.jet_SeedIPhi , "jet_SeedIPhi[jet_num]/I");
		tree->Branch("jet_SeedIEta" , &stuple.jet_SeedIEta , "jet_SeedIEta[jet_num]/I");
		tree->Branch("jet_EmTime"   , &stuple.jet_EmTime   , "jet_EmTime[jet_num]/F");

		// RAW JET COLLECTIONS
		tree->Branch("jet_raw_num"   , &stuple.jet_raw_num   , "jet_raw_num/i");
		tree->Branch("jet_raw_Index" , &stuple.jet_raw_Index , "jet_raw_Index[jet_raw_num]/I");
		tree->Branch("jet_raw_Pt"    , &stuple.jet_raw_Pt    , "jet_raw_Pt[jet_raw_num]/F");
		tree->Branch("jet_raw_E"     , &stuple.jet_raw_E     , "jet_raw_E[jet_raw_num]/F");
		tree->Branch("jet_raw_Px"    , &stuple.jet_raw_Px    , "jet_raw_Px[jet_raw_num]/F");
		tree->Branch("jet_raw_Py"    , &stuple.jet_raw_Py    , "jet_raw_Py[jet_raw_num]/F");
		tree->Branch("jet_raw_Pz"    , &stuple.jet_raw_Pz    , "jet_raw_Pz[jet_raw_num]/F");


		
			// JES UP JET COLLECTION
		tree->Branch("jet_up_num"     , &stuple.jet_up_num      , "jet_up_num/i");
		tree->Branch("jet_up_NJet15"  , &stuple.jet_up_NJet15   , "jet_up_NJet15/i");
		tree->Branch("jet_up_Index"   , &stuple.jet_up_Index    , "jet_up_Index[jet_up_num]/I");
		tree->Branch("jet_up_Pt"      , &stuple.jet_up_Pt       , "jet_up_Pt[jet_up_num]/F");
		tree->Branch("jet_up_E"       , &stuple.jet_up_E        , "jet_up_E[jet_up_num]/F");
		tree->Branch("jet_up_Px"      , &stuple.jet_up_Px       , "jet_up_Px[jet_up_num]/F");
		tree->Branch("jet_up_Py"      , &stuple.jet_up_Py       , "jet_up_Py[jet_up_num]/F");
		tree->Branch("jet_up_Pz"      , &stuple.jet_up_Pz       , "jet_up_Pz[jet_up_num]/F");
		tree->Branch("jet_up_DetEta"  , &stuple.jet_up_DetEta   , "jet_up_DetEta[jet_up_num]/F");
		tree->Branch("jet_up_DetPhi"  , &stuple.jet_up_DetPhi   , "jet_up_DetPhi[jet_up_num]/F");
		tree->Branch("jet_up_HadEm"   , &stuple.jet_up_HadEm    , "jet_up_HadEm[jet_up_num]/F");
		tree->Branch("jet_up_Emfr"    , &stuple.jet_up_Emfr     , "jet_up_Emfr[jet_up_num]/F");
		tree->Branch("jet_up_Ntowers" , &stuple.jet_up_Ntowers  , "jet_up_Ntowers[jet_up_num]/I");
		tree->Branch("jet_up_Ntracks" , &stuple.jet_up_Ntracks  , "jet_up_Ntracks[jet_up_num]/I");
		tree->Branch("jet_up_SeedIPhi", &stuple.jet_up_SeedIPhi , "jet_up_SeedIPhi[jet_up_num]/I");
		tree->Branch("jet_up_SeedIEta", &stuple.jet_up_SeedIEta , "jet_up_SeedIEta[jet_up_num]/I");
		tree->Branch("jet_up_EmTime"  , &stuple.jet_up_EmTime   , "jet_up_EmTime[jet_up_num]/F");
			
			// JES DOWN JET COLLECTION
		tree->Branch("jet_down_num"     , &stuple.jet_down_num      , "jet_down_num/i");
		tree->Branch("jet_down_NJet15"  , &stuple.jet_down_NJet15   , "jet_down_NJet15/i");
		tree->Branch("jet_down_Index"   , &stuple.jet_down_Index    , "jet_down_Index[jet_down_num]/I");
		tree->Branch("jet_down_Pt"      , &stuple.jet_down_Pt       , "jet_down_Pt[jet_down_num]/F");
		tree->Branch("jet_down_E"       , &stuple.jet_down_E        , "jet_down_E[jet_down_num]/F");
		tree->Branch("jet_down_Px"      , &stuple.jet_down_Px       , "jet_down_Px[jet_down_num]/F");
		tree->Branch("jet_down_Py"      , &stuple.jet_down_Py       , "jet_down_Py[jet_down_num]/F");
		tree->Branch("jet_down_Pz"      , &stuple.jet_down_Pz       , "jet_down_Pz[jet_down_num]/F");
		tree->Branch("jet_down_DetEta"  , &stuple.jet_down_DetEta   , "jet_down_DetEta[jet_down_num]/F");
		tree->Branch("jet_down_DetPhi"  , &stuple.jet_down_DetPhi   , "jet_down_DetPhi[jet_down_num]/F");
		tree->Branch("jet_down_HadEm"   , &stuple.jet_down_HadEm    , "jet_down_HadEm[jet_down_num]/F");
		tree->Branch("jet_down_Emfr"    , &stuple.jet_down_Emfr     , "jet_down_Emfr[jet_down_num]/F");
		tree->Branch("jet_down_Ntowers" , &stuple.jet_down_Ntowers  , "jet_down_Ntowers[jet_down_num]/I");
		tree->Branch("jet_down_Ntracks" , &stuple.jet_down_Ntracks  , "jet_down_Ntracks[jet_down_num]/I");
		tree->Branch("jet_down_SeedIPhi", &stuple.jet_down_SeedIPhi , "jet_down_SeedIPhi[jet_down_num]/I");
		tree->Branch("jet_down_SeedIEta", &stuple.jet_down_SeedIEta , "jet_down_SeedIEta[jet_down_num]/I");
		tree->Branch("jet_down_EmTime"  , &stuple.jet_down_EmTime   , "jet_down_EmTime[jet_down_num]/F");
	

			// VERTEX INFO
		tree->Branch("vtx_N"        , &stuple.vtx_N        , "vtx_N/I");
		tree->Branch("vtx_NClass12" , &stuple.vtx_NClass12 , "vtx_NClass12/I");
		tree->Branch("vtx_z"        , &stuple.vtx_z        , "vtx_z/F");
		tree->Branch("vtx_Ntracks"  , &stuple.vtx_Ntracks  , "vtx_Ntracks/I");
		tree->Branch("vtx_SumPt"    , &stuple.vtx_SumPt    , "vtx_SumPt/F");

			// MET INFO
		tree->Branch("met_Met"    , &stuple.met_Met    , "met_Met/F");
		tree->Branch("met_RawMet"    , &stuple.met_RawMet    , "met_RawMet/F");
		tree->Branch("met_MetX"   , &stuple.met_MetX   , "met_MetX/F");
		tree->Branch("met_MetY"   , &stuple.met_MetY   , "met_MetY/F");
		tree->Branch("met_SumEt"  , &stuple.met_SumEt  , "met_SumEt/F");
		tree->Branch("met_Ht"     , &stuple.met_Ht     , "met_Ht/F");
		tree->Branch("met_MetPhi" , &stuple.met_MetPhi , "met_MetPhi/F");

			// GENERATED MET FROM MET MODEL
//		tree->Branch("met_Gen_d"   , &stuple.met_Gen_d   , "met_Gen_d/F");
//		tree->Branch("met_Gen_m"   , &stuple.met_Gen_m   , "met_Gen_m/F");
//		tree->Branch("met_Gen_p"   , &stuple.met_Gen_p   , "met_Gen_p/F");
//		tree->Branch("met_Gen_mUn" , &stuple.met_Gen_mUn , "met_Gen_mUn/F");
//		tree->Branch("met_Gen_pUn" , &stuple.met_Gen_pUn , "met_Gen_pUn/F");

		
			// Generator level info
		tree->Branch("gen_elenum" , &stuple.gen_elenum , "gen_elenum/i");
		tree->Branch("gen_phonum" , &stuple.gen_phonum , "gen_phonum/i");
		
		tree->Branch("gen_MomIndex" , &stuple.gen_MomIndex  , "gen_MomIndex/I");
		tree->Branch("gen_MomPDG"   , &stuple.gen_MomPDG    , "gen_MomPDG/I");
		tree->Branch("gen_MomStaus" , &stuple.gen_MomStatus , "gen_MomStatus/I");
		tree->Branch("gen_MomEtc"   , &stuple.gen_MomEtc    , "gen_MomEtc/F");
		tree->Branch("gen_MomE"     , &stuple.gen_MomE      , "gen_MomE/F");
		tree->Branch("gen_MomPx"    , &stuple.gen_MomPx     , "gen_MomPx/F");
		tree->Branch("gen_MomPy"    , &stuple.gen_MomPy     , "gen_MomPy/F");
		tree->Branch("gen_MomPz"    , &stuple.gen_MomPz     , "gen_MomPz/F");

		tree->Branch("gen_ProdVtxX" , &stuple.gen_ProdVtxX , "gen_ProdVtxX/F");
		tree->Branch("gen_ProdVtxY" , &stuple.gen_ProdVtxY , "gen_ProdVtxY/F");
		tree->Branch("gen_ProdVtxZ" , &stuple.gen_ProdVtxZ , "gen_ProdVtxZ/F");
		tree->Branch("gen_ProdVtxT" , &stuple.gen_ProdVtxT , "gen_ProdVtxT/F");

		tree->Branch("gen_pho_Index"  , &stuple.gen_pho_Index  , "gen_pho_Index[gen_phonum]/I");
		tree->Branch("gen_pho_PDG"    , &stuple.gen_pho_PDG    , "gen_pho_PDG[gen_phonum]/I");
		tree->Branch("gen_pho_Status" , &stuple.gen_pho_Status , "gen_pho_Status[gen_phonum]/I");
		tree->Branch("gen_pho_Etc"    , &stuple.gen_pho_Etc    , "gen_pho_Etc[gen_phonum]/F");
		tree->Branch("gen_pho_E"      , &stuple.gen_pho_E      , "gen_pho_E[gen_phonum]/F");
		tree->Branch("gen_pho_Px"     , &stuple.gen_pho_Px     , "gen_pho_Px[gen_phonum]/F");
		tree->Branch("gen_pho_Py"     , &stuple.gen_pho_Py     , "gen_pho_Py[gen_phonum]/F");
		tree->Branch("gen_pho_Pz"     , &stuple.gen_pho_Pz     , "gen_pho_Pz[gen_phonum]/F");

		tree->Branch("gen_ele_Index"  , &stuple.gen_ele_Index  , "gen_ele_Index[gen_elenum]/I");
		tree->Branch("gen_ele_PDG"    , &stuple.gen_ele_PDG    , "gen_ele_PDG[gen_elenum]/I");
		tree->Branch("gen_ele_Status" , &stuple.gen_ele_Status , "gen_ele_Status[gen_elenum]/I");
		tree->Branch("gen_ele_Etc"    , &stuple.gen_ele_Etc    , "gen_ele_Etc[gen_elenum]/F");
		tree->Branch("gen_ele_E"      , &stuple.gen_ele_E      , "gen_ele_E[gen_elenum]/F");
		tree->Branch("gen_ele_Px"     , &stuple.gen_ele_Px     , "gen_ele_Px[gen_elenum]/F");
		tree->Branch("gen_ele_Py"     , &stuple.gen_ele_Py     , "gen_ele_Py[gen_elenum]/F");
		tree->Branch("gen_ele_Pz"     , &stuple.gen_ele_Pz     , "gen_ele_Pz[gen_elenum]/F");

	std::cout << __FILE__ << ":: tree " << tree << std::endl;
	return tree;
}

/*
void PrepStupleTreeForReading(TChain *tree, Stuple *stuple)
{
	if (!tree) {
		std::cout << "tree not found" << std::endl;
		return;
	}
	if (!stuple) {
		std::cout << "Stuple obj not found" << std::endl;
		return;
	}

	tree->SetBranchStatus("*", 1);

		tree->SetBranchAddress("evt_McFlag"   	 , &(stuple->evt_McFlag));
		tree->SetBranchAddress("evt_RunNumber"   , &(stuple->evt_RunNumber));
		tree->SetBranchAddress("evt_EventNumber" , &(stuple->evt_EventNumber));

		tree->SetBranchAddress("tri_pho25iso" , &(stuple->tri_pho25iso));
		tree->SetBranchAddress("tri_pho50"    , &(stuple->tri_pho50));
		tree->SetBranchAddress("tri_pho70"    , &(stuple->tri_pho70));
	
		tree->SetBranchAddress("pho_num"     , &(stuple->pho_num));
		tree->SetBranchAddress("pho_Ntight"  , &(stuple->pho_Ntight));
		tree->SetBranchAddress("pho_Nloose"  , &(stuple->pho_Nloose));
		tree->SetBranchAddress("pho_Index"   , &(stuple->pho_Index));
		tree->SetBranchAddress("pho_PhoBlockIndex"   , &(stuple->pho_PhoBlockIndex));
		tree->SetBranchAddress("pho_Etc"     , &(stuple->pho_Etc));
		tree->SetBranchAddress("pho_E"       , &(stuple->pho_E));
		tree->SetBranchAddress("pho_Px"      , &(stuple->pho_Px));
		tree->SetBranchAddress("pho_Py"      , &(stuple->pho_Py));
		tree->SetBranchAddress("pho_Pz"      , &(stuple->pho_Pz));
		tree->SetBranchAddress("pho_Detector", &(stuple->pho_Detector));
		tree->SetBranchAddress("pho_DetEta"  , &(stuple->pho_DetEta));
		tree->SetBranchAddress("pho_DetPhi"  , &(stuple->pho_DetPhi));
		tree->SetBranchAddress("pho_XCes"    , &(stuple->pho_XCes));
		tree->SetBranchAddress("pho_ZCes"    , &(stuple->pho_ZCes));
		tree->SetBranchAddress("pho_HadEm"   , &(stuple->pho_HadEm));
		tree->SetBranchAddress("pho_Chi2Mean", &(stuple->pho_Chi2Mean));
		tree->SetBranchAddress("pho_N3d"     , &(stuple->pho_N3d));
		tree->SetBranchAddress("pho_Iso4"    , &(stuple->pho_Iso4));
		tree->SetBranchAddress("pho_TrkPt"   , &(stuple->pho_TrkPt));
		tree->SetBranchAddress("pho_TrkIso"  , &(stuple->pho_TrkIso));
		tree->SetBranchAddress("pho_CesWireE2"      , &(stuple->pho_CesWireE2));
		tree->SetBranchAddress("pho_CesStripE2"     , &(stuple->pho_CesStripE2));
		tree->SetBranchAddress("pho_PhiWedge"       , &(stuple->pho_PhiWedge));
		tree->SetBranchAddress("pho_NMuonStubs"     , &(stuple->pho_NMuonStubs));
		tree->SetBranchAddress("pho_EmTime"         , &(stuple->pho_EmTime));
		tree->SetBranchAddress("pho_TightId"        , &(stuple->pho_TightId));
		tree->SetBranchAddress("pho_LooseId"        , &(stuple->pho_LooseId));
		tree->SetBranchAddress("pho_PhoenixId"      , &(stuple->pho_PhoenixId));
		tree->SetBranchAddress("pho_Halo_seedWedge" , &(stuple->pho_Halo_seedWedge));
		tree->SetBranchAddress("pho_Halo_eastNhad"  , &(stuple->pho_Halo_eastNhad));
		tree->SetBranchAddress("pho_Halo_westNhad"  , &(stuple->pho_Halo_westNhad));
		tree->SetBranchAddress("pho_matchJetIndex"  , &(stuple->pho_matchJetIndex));
		tree->SetBranchAddress("pho_CprWgt"  , &(stuple->pho_CprWgt));
		tree->SetBranchAddress("pho_CprSys1"  , &(stuple->pho_CprSys1));
		tree->SetBranchAddress("pho_CprSys2"  , &(stuple->pho_CprSys2));
		tree->SetBranchAddress("pho_CprSys3"  , &(stuple->pho_CprSys3));
		tree->SetBranchAddress("pho_CprSys4"  , &(stuple->pho_CprSys4));
		tree->SetBranchAddress("pho_CprSys5"  , &(stuple->pho_CprSys5));
		tree->SetBranchAddress("pho_CprSys6"  , &(stuple->pho_CprSys6));
		tree->SetBranchAddress("pho_CprSys7"  , &(stuple->pho_CprSys7));
		tree->SetBranchAddress("pho_CprSys8"  , &(stuple->pho_CprSys8));
		

			// Pho EM E uncertainty up 1%
		tree->SetBranchAddress("pho_up_num"     , &(stuple->pho_up_num));
		tree->SetBranchAddress("pho_up_Ntight"  , &(stuple->pho_up_Ntight));
		tree->SetBranchAddress("pho_up_Nloose"  , &(stuple->pho_up_Nloose));
		tree->SetBranchAddress("pho_up_Index"   , &(stuple->pho_up_Index));
		tree->SetBranchAddress("pho_up_PhoBlockIndex"   , &(stuple->pho_up_PhoBlockIndex));
		tree->SetBranchAddress("pho_up_Etc"     , &(stuple->pho_up_Etc));
		tree->SetBranchAddress("pho_up_E"       , &(stuple->pho_up_E));
		tree->SetBranchAddress("pho_up_Px"      , &(stuple->pho_up_Px));
		tree->SetBranchAddress("pho_up_Py"      , &(stuple->pho_up_Py));
		tree->SetBranchAddress("pho_up_Pz"      , &(stuple->pho_up_Pz));
		tree->SetBranchAddress("pho_up_Detector", &(stuple->pho_up_Detector));
		tree->SetBranchAddress("pho_up_DetEta"  , &(stuple->pho_up_DetEta));
		tree->SetBranchAddress("pho_up_DetPhi"  , &(stuple->pho_up_DetPhi));
		tree->SetBranchAddress("pho_up_XCes"    , &(stuple->pho_up_XCes));
		tree->SetBranchAddress("pho_up_ZCes"    , &(stuple->pho_up_ZCes));
		tree->SetBranchAddress("pho_up_HadEm"   , &(stuple->pho_up_HadEm));
		tree->SetBranchAddress("pho_up_Chi2Mean", &(stuple->pho_up_Chi2Mean));
		tree->SetBranchAddress("pho_up_N3d"     , &(stuple->pho_up_N3d));
		tree->SetBranchAddress("pho_up_Iso4"    , &(stuple->pho_up_Iso4));
		tree->SetBranchAddress("pho_up_TrkPt"   , &(stuple->pho_up_TrkPt));
		tree->SetBranchAddress("pho_up_TrkIso"  , &(stuple->pho_up_TrkIso));
		tree->SetBranchAddress("pho_up_CesWireE2"      , &(stuple->pho_up_CesWireE2));
		tree->SetBranchAddress("pho_up_CesStripE2"     , &(stuple->pho_up_CesStripE2));
		tree->SetBranchAddress("pho_up_PhiWedge"       , &(stuple->pho_up_PhiWedge));
		tree->SetBranchAddress("pho_up_NMuonStubs"     , &(stuple->pho_up_NMuonStubs));
		tree->SetBranchAddress("pho_up_EmTime"         , &(stuple->pho_up_EmTime));
		tree->SetBranchAddress("pho_up_TightId"        , &(stuple->pho_up_TightId));
		tree->SetBranchAddress("pho_up_LooseId"        , &(stuple->pho_up_LooseId));
		tree->SetBranchAddress("pho_up_PhoenixId"      , &(stuple->pho_up_PhoenixId));
		tree->SetBranchAddress("pho_up_Halo_seedWedge" , &(stuple->pho_up_Halo_seedWedge));
		tree->SetBranchAddress("pho_up_Halo_eastNhad"  , &(stuple->pho_up_Halo_eastNhad));
		tree->SetBranchAddress("pho_up_Halo_westNhad"  , &(stuple->pho_up_Halo_westNhad));
		tree->SetBranchAddress("pho_up_matchJetIndex"  , &(stuple->pho_up_matchJetIndex));
		tree->SetBranchAddress("pho_up_CprWgt"  , &(stuple->pho_up_CprWgt));
		tree->SetBranchAddress("pho_up_CprSys1"  , &(stuple->pho_up_CprSys1));
		tree->SetBranchAddress("pho_up_CprSys2"  , &(stuple->pho_up_CprSys2));
		tree->SetBranchAddress("pho_up_CprSys3"  , &(stuple->pho_up_CprSys3));
		tree->SetBranchAddress("pho_up_CprSys4"  , &(stuple->pho_up_CprSys4));
		tree->SetBranchAddress("pho_up_CprSys5"  , &(stuple->pho_up_CprSys5));
		tree->SetBranchAddress("pho_up_CprSys6"  , &(stuple->pho_up_CprSys6));
		tree->SetBranchAddress("pho_up_CprSys7"  , &(stuple->pho_up_CprSys7));
		tree->SetBranchAddress("pho_up_CprSys8"  , &(stuple->pho_up_CprSys8));
	
			
			// Pho EM E uncertainty down 1%
		tree->SetBranchAddress("pho_down_num"     , &(stuple->pho_down_num));
		tree->SetBranchAddress("pho_down_Ntight"  , &(stuple->pho_down_Ntight));
		tree->SetBranchAddress("pho_down_Nloose"  , &(stuple->pho_down_Nloose));
		tree->SetBranchAddress("pho_down_Index"   , &(stuple->pho_down_Index));
		tree->SetBranchAddress("pho_down_PhoBlockIndex"   , &(stuple->pho_down_PhoBlockIndex));
		tree->SetBranchAddress("pho_down_Etc"     , &(stuple->pho_down_Etc));
		tree->SetBranchAddress("pho_down_E"       , &(stuple->pho_down_E));
		tree->SetBranchAddress("pho_down_Px"      , &(stuple->pho_down_Px));
		tree->SetBranchAddress("pho_down_Py"      , &(stuple->pho_down_Py));
		tree->SetBranchAddress("pho_down_Pz"      , &(stuple->pho_down_Pz));
		tree->SetBranchAddress("pho_down_Detector", &(stuple->pho_down_Detector));
		tree->SetBranchAddress("pho_down_DetEta"  , &(stuple->pho_down_DetEta));
		tree->SetBranchAddress("pho_down_DetPhi"  , &(stuple->pho_down_DetPhi));
		tree->SetBranchAddress("pho_down_XCes"    , &(stuple->pho_down_XCes));
		tree->SetBranchAddress("pho_down_ZCes"    , &(stuple->pho_down_ZCes));
		tree->SetBranchAddress("pho_down_HadEm"   , &(stuple->pho_down_HadEm));
		tree->SetBranchAddress("pho_down_Chi2Mean", &(stuple->pho_down_Chi2Mean));
		tree->SetBranchAddress("pho_down_N3d"     , &(stuple->pho_down_N3d));
		tree->SetBranchAddress("pho_down_Iso4"    , &(stuple->pho_down_Iso4));
		tree->SetBranchAddress("pho_down_TrkPt"   , &(stuple->pho_down_TrkPt));
		tree->SetBranchAddress("pho_down_TrkIso"  , &(stuple->pho_down_TrkIso));
		tree->SetBranchAddress("pho_down_CesWireE2"      , &(stuple->pho_down_CesWireE2));
		tree->SetBranchAddress("pho_down_CesStripE2"     , &(stuple->pho_down_CesStripE2));
		tree->SetBranchAddress("pho_down_PhiWedge"       , &(stuple->pho_down_PhiWedge));
		tree->SetBranchAddress("pho_down_NMuonStubs"     , &(stuple->pho_down_NMuonStubs));
		tree->SetBranchAddress("pho_down_EmTime"         , &(stuple->pho_down_EmTime));
		tree->SetBranchAddress("pho_down_TightId"        , &(stuple->pho_down_TightId));
		tree->SetBranchAddress("pho_down_LooseId"        , &(stuple->pho_down_LooseId));
		tree->SetBranchAddress("pho_down_PhoenixId"      , &(stuple->pho_down_PhoenixId));
		tree->SetBranchAddress("pho_down_Halo_seedWedge" , &(stuple->pho_down_Halo_seedWedge));
		tree->SetBranchAddress("pho_down_Halo_eastNhad"  , &(stuple->pho_down_Halo_eastNhad));
		tree->SetBranchAddress("pho_down_Halo_westNhad"  , &(stuple->pho_down_Halo_westNhad));
		tree->SetBranchAddress("pho_down_matchJetIndex"  , &(stuple->pho_down_matchJetIndex));
		tree->SetBranchAddress("pho_down_CprWgt"  , &(stuple->pho_down_CprWgt));
		tree->SetBranchAddress("pho_down_CprSys1"  , &(stuple->pho_down_CprSys1));
		tree->SetBranchAddress("pho_down_CprSys2"  , &(stuple->pho_down_CprSys2));
		tree->SetBranchAddress("pho_down_CprSys3"  , &(stuple->pho_down_CprSys3));
		tree->SetBranchAddress("pho_down_CprSys4"  , &(stuple->pho_down_CprSys4));
		tree->SetBranchAddress("pho_down_CprSys5"  , &(stuple->pho_down_CprSys5));
		tree->SetBranchAddress("pho_down_CprSys6"  , &(stuple->pho_down_CprSys6));
		tree->SetBranchAddress("pho_down_CprSys7"  , &(stuple->pho_down_CprSys7));
		tree->SetBranchAddress("pho_down_CprSys8"  , &(stuple->pho_down_CprSys8));
	


		//electron collections	

		tree->SetBranchAddress("ele_num"     , &(stuple->ele_num));
		tree->SetBranchAddress("ele_Ntight"  , &(stuple->ele_Ntight));
		tree->SetBranchAddress("ele_Nloose"  , &(stuple->ele_Nloose));
		tree->SetBranchAddress("ele_Index"   , &(stuple->ele_Index));
		tree->SetBranchAddress("ele_PhoBlockIndex", &(stuple->ele_PhoBlockIndex));
		tree->SetBranchAddress("ele_EleBlockIndex", &(stuple->ele_EleBlockIndex));
		tree->SetBranchAddress("ele_Etc"     , &(stuple->ele_Etc));
		tree->SetBranchAddress("ele_E"       , &(stuple->ele_E));
		tree->SetBranchAddress("ele_Px"      , &(stuple->ele_Px));
		tree->SetBranchAddress("ele_Py"      , &(stuple->ele_Py));
		tree->SetBranchAddress("ele_Pz"      , &(stuple->ele_Pz));
		tree->SetBranchAddress("ele_Detector", &(stuple->ele_Detector));
		tree->SetBranchAddress("ele_DetEta"  , &(stuple->ele_DetEta));
		tree->SetBranchAddress("ele_DetPhi"  , &(stuple->ele_DetPhi));
		tree->SetBranchAddress("ele_XCes"    , &(stuple->ele_XCes));
		tree->SetBranchAddress("ele_ZCes"    , &(stuple->ele_ZCes));
		tree->SetBranchAddress("ele_HadEm"   , &(stuple->ele_HadEm));
		tree->SetBranchAddress("ele_Chi2Mean", &(stuple->ele_Chi2Mean));
		tree->SetBranchAddress("ele_N3d"     , &(stuple->ele_N3d));
		tree->SetBranchAddress("ele_Iso4"    , &(stuple->ele_Iso4));
		tree->SetBranchAddress("ele_TrkIso"  , &(stuple->ele_TrkIso));
		tree->SetBranchAddress("ele_CesWireE2"      , &(stuple->ele_CesWireE2));
		tree->SetBranchAddress("ele_CesStripE2"     , &(stuple->ele_CesStripE2));
		tree->SetBranchAddress("ele_PhiWedge"       , &(stuple->ele_PhiWedge));
		tree->SetBranchAddress("ele_NMuonStubs"     , &(stuple->ele_NMuonStubs));
		tree->SetBranchAddress("ele_EmTime"         , &(stuple->ele_EmTime));
		tree->SetBranchAddress("ele_PhoenixId"      , &(stuple->ele_PhoenixId));
		tree->SetBranchAddress("ele_Halo_seedWedge" , &(stuple->ele_Halo_seedWedge));
		tree->SetBranchAddress("ele_Halo_eastNhad"  , &(stuple->ele_Halo_eastNhad));
		tree->SetBranchAddress("ele_Halo_westNhad"  , &(stuple->ele_Halo_westNhad));
		tree->SetBranchAddress("ele_matchJetIndex"  , &(stuple->ele_matchJetIndex));
		tree->SetBranchAddress("ele_Ntracks"      , &(stuple->ele_Ntracks));
		tree->SetBranchAddress("ele_Emfr"         , &(stuple->ele_Emfr));
		tree->SetBranchAddress("ele_EoverP"       , &(stuple->ele_EoverP));
		tree->SetBranchAddress("ele_TrackPt"      , &(stuple->ele_TrackPt));
		tree->SetBranchAddress("ele_TrackBcPt"    , &(stuple->ele_TrackBcPt));
		tree->SetBranchAddress("ele_TrackPhi"     , &(stuple->ele_TrackPhi));
		tree->SetBranchAddress("ele_Nssl"         , &(stuple->ele_Nssl));
		tree->SetBranchAddress("ele_Nasl"         , &(stuple->ele_Nasl));
		tree->SetBranchAddress("ele_TightId"      , &(stuple->ele_TightId));
		tree->SetBranchAddress("ele_LooseId"      , &(stuple->ele_LooseId));
		tree->SetBranchAddress("ele_ConversionId" , &(stuple->ele_ConversionId));
		tree->SetBranchAddress("ele_StdTightId"      , &(stuple->ele_StdTightId));
		tree->SetBranchAddress("ele_StdLooseId"      , &(stuple->ele_StdLooseId));


			// Ele EM E uncertainty up 1%
		tree->SetBranchAddress("ele_up_num"     , &(stuple->ele_up_num));
		tree->SetBranchAddress("ele_up_Ntight"  , &(stuple->ele_up_Ntight));
		tree->SetBranchAddress("ele_up_Nloose"  , &(stuple->ele_up_Nloose));
		tree->SetBranchAddress("ele_up_Index"   , &(stuple->ele_up_Index));
		tree->SetBranchAddress("ele_up_PhoBlockIndex", &(stuple->ele_up_PhoBlockIndex));
		tree->SetBranchAddress("ele_up_EleBlockIndex", &(stuple->ele_up_EleBlockIndex));
		tree->SetBranchAddress("ele_up_Etc"     , &(stuple->ele_up_Etc));
		tree->SetBranchAddress("ele_up_E"       , &(stuple->ele_up_E));
		tree->SetBranchAddress("ele_up_Px"      , &(stuple->ele_up_Px));
		tree->SetBranchAddress("ele_up_Py"      , &(stuple->ele_up_Py));
		tree->SetBranchAddress("ele_up_Pz"      , &(stuple->ele_up_Pz));
		tree->SetBranchAddress("ele_up_Detector", &(stuple->ele_up_Detector));
		tree->SetBranchAddress("ele_up_DetEta"  , &(stuple->ele_up_DetEta));
		tree->SetBranchAddress("ele_up_DetPhi"  , &(stuple->ele_up_DetPhi));
		tree->SetBranchAddress("ele_up_XCes"    , &(stuple->ele_up_XCes));
		tree->SetBranchAddress("ele_up_ZCes"    , &(stuple->ele_up_ZCes));
		tree->SetBranchAddress("ele_up_HadEm"   , &(stuple->ele_up_HadEm));
		tree->SetBranchAddress("ele_up_Chi2Mean", &(stuple->ele_up_Chi2Mean));
		tree->SetBranchAddress("ele_up_N3d"     , &(stuple->ele_up_N3d));
		tree->SetBranchAddress("ele_up_Iso4"    , &(stuple->ele_up_Iso4));
		tree->SetBranchAddress("ele_up_TrkIso"  , &(stuple->ele_up_TrkIso));
		tree->SetBranchAddress("ele_up_CesWireE2"      , &(stuple->ele_up_CesWireE2));
		tree->SetBranchAddress("ele_up_CesStripE2"     , &(stuple->ele_up_CesStripE2));
		tree->SetBranchAddress("ele_up_PhiWedge"       , &(stuple->ele_up_PhiWedge));
		tree->SetBranchAddress("ele_up_NMuonStubs"     , &(stuple->ele_up_NMuonStubs));
		tree->SetBranchAddress("ele_up_EmTime"         , &(stuple->ele_up_EmTime));
		tree->SetBranchAddress("ele_up_PhoenixId"      , &(stuple->ele_up_PhoenixId));
		tree->SetBranchAddress("ele_up_Halo_seedWedge" , &(stuple->ele_up_Halo_seedWedge));
		tree->SetBranchAddress("ele_up_Halo_eastNhad"  , &(stuple->ele_up_Halo_eastNhad));
		tree->SetBranchAddress("ele_up_Halo_westNhad"  , &(stuple->ele_up_Halo_westNhad));
		tree->SetBranchAddress("ele_up_matchJetIndex"  , &(stuple->ele_up_matchJetIndex));
		tree->SetBranchAddress("ele_up_Ntracks"      , &(stuple->ele_up_Ntracks));
		tree->SetBranchAddress("ele_up_Emfr"         , &(stuple->ele_up_Emfr));
		tree->SetBranchAddress("ele_up_EoverP"       , &(stuple->ele_up_EoverP));
		tree->SetBranchAddress("ele_up_TrackPt"      , &(stuple->ele_up_TrackPt));
		tree->SetBranchAddress("ele_up_TrackBcPt"    , &(stuple->ele_up_TrackBcPt));
		tree->SetBranchAddress("ele_up_TrackPhi"     , &(stuple->ele_up_TrackPhi));
		tree->SetBranchAddress("ele_up_Nssl"         , &(stuple->ele_up_Nssl));
		tree->SetBranchAddress("ele_up_Nasl"         , &(stuple->ele_up_Nasl));
		tree->SetBranchAddress("ele_up_TightId"      , &(stuple->ele_up_TightId));
		tree->SetBranchAddress("ele_up_LooseId"      , &(stuple->ele_up_LooseId));
		tree->SetBranchAddress("ele_up_ConversionId" , &(stuple->ele_up_ConversionId));
		tree->SetBranchAddress("ele_up_StdTightId"      , &(stuple->ele_up_StdTightId));
		tree->SetBranchAddress("ele_up_StdLooseId"      , &(stuple->ele_up_StdLooseId));

			
			// Ele EM E uncertainty down 1%
		tree->SetBranchAddress("ele_down_num"     , &(stuple->ele_down_num));
		tree->SetBranchAddress("ele_down_Ntight"  , &(stuple->ele_down_Ntight));
		tree->SetBranchAddress("ele_down_Nloose"  , &(stuple->ele_down_Nloose));
		tree->SetBranchAddress("ele_down_Index"   , &(stuple->ele_down_Index));
		tree->SetBranchAddress("ele_down_PhoBlockIndex", &(stuple->ele_down_PhoBlockIndex));
		tree->SetBranchAddress("ele_down_EleBlockIndex", &(stuple->ele_down_EleBlockIndex));
		tree->SetBranchAddress("ele_down_Etc"     , &(stuple->ele_down_Etc));
		tree->SetBranchAddress("ele_down_E"       , &(stuple->ele_down_E));
		tree->SetBranchAddress("ele_down_Px"      , &(stuple->ele_down_Px));
		tree->SetBranchAddress("ele_down_Py"      , &(stuple->ele_down_Py));
		tree->SetBranchAddress("ele_down_Pz"      , &(stuple->ele_down_Pz));
		tree->SetBranchAddress("ele_down_Detector", &(stuple->ele_down_Detector));
		tree->SetBranchAddress("ele_down_DetEta"  , &(stuple->ele_down_DetEta));
		tree->SetBranchAddress("ele_down_DetPhi"  , &(stuple->ele_down_DetPhi));
		tree->SetBranchAddress("ele_down_XCes"    , &(stuple->ele_down_XCes));
		tree->SetBranchAddress("ele_down_ZCes"    , &(stuple->ele_down_ZCes));
		tree->SetBranchAddress("ele_down_HadEm"   , &(stuple->ele_down_HadEm));
		tree->SetBranchAddress("ele_down_Chi2Mean", &(stuple->ele_down_Chi2Mean));
		tree->SetBranchAddress("ele_down_N3d"     , &(stuple->ele_down_N3d));
		tree->SetBranchAddress("ele_down_Iso4"    , &(stuple->ele_down_Iso4));
		tree->SetBranchAddress("ele_down_TrkIso"  , &(stuple->ele_down_TrkIso));
		tree->SetBranchAddress("ele_down_CesWireE2"      , &(stuple->ele_down_CesWireE2));
		tree->SetBranchAddress("ele_down_CesStripE2"     , &(stuple->ele_down_CesStripE2));
		tree->SetBranchAddress("ele_down_PhiWedge"       , &(stuple->ele_down_PhiWedge));
		tree->SetBranchAddress("ele_down_NMuonStubs"     , &(stuple->ele_down_NMuonStubs));
		tree->SetBranchAddress("ele_down_EmTime"         , &(stuple->ele_down_EmTime));
		tree->SetBranchAddress("ele_down_PhoenixId"      , &(stuple->ele_down_PhoenixId));
		tree->SetBranchAddress("ele_down_Halo_seedWedge" , &(stuple->ele_down_Halo_seedWedge));
		tree->SetBranchAddress("ele_down_Halo_eastNhad"  , &(stuple->ele_down_Halo_eastNhad));
		tree->SetBranchAddress("ele_down_Halo_westNhad"  , &(stuple->ele_down_Halo_westNhad));
		tree->SetBranchAddress("ele_down_matchJetIndex"  , &(stuple->ele_down_matchJetIndex));
		tree->SetBranchAddress("ele_down_Ntracks"      , &(stuple->ele_down_Ntracks));
		tree->SetBranchAddress("ele_down_Emfr"         , &(stuple->ele_down_Emfr));
		tree->SetBranchAddress("ele_down_EoverP"       , &(stuple->ele_down_EoverP));
		tree->SetBranchAddress("ele_down_TrackPt"      , &(stuple->ele_down_TrackPt));
		tree->SetBranchAddress("ele_down_TrackBcPt"    , &(stuple->ele_down_TrackBcPt));
		tree->SetBranchAddress("ele_down_TrackPhi"     , &(stuple->ele_down_TrackPhi));
		tree->SetBranchAddress("ele_down_Nssl"         , &(stuple->ele_down_Nssl));
		tree->SetBranchAddress("ele_down_Nasl"         , &(stuple->ele_down_Nasl));
		tree->SetBranchAddress("ele_down_TightId"      , &(stuple->ele_down_TightId));
		tree->SetBranchAddress("ele_down_LooseId"      , &(stuple->ele_down_LooseId));
		tree->SetBranchAddress("ele_down_ConversionId" , &(stuple->ele_down_ConversionId));
		tree->SetBranchAddress("ele_down_StdTightId"      , &(stuple->ele_down_StdTightId));
		tree->SetBranchAddress("ele_down_StdLooseId"      , &(stuple->ele_down_StdLooseId));


		// JET COLLECTIONS
		tree->SetBranchAddress("jet_num"      , &(stuple->jet_num));
		tree->SetBranchAddress("jet_NJet15"   , &(stuple->jet_NJet15));
		tree->SetBranchAddress("jet_Index"    , &(stuple->jet_Index));
		tree->SetBranchAddress("jet_Pt"       , &(stuple->jet_Pt));
		tree->SetBranchAddress("jet_E"        , &(stuple->jet_E));
		tree->SetBranchAddress("jet_Px"       , &(stuple->jet_Px));
		tree->SetBranchAddress("jet_Py"       , &(stuple->jet_Py));
		tree->SetBranchAddress("jet_Pz"       , &(stuple->jet_Pz));
		tree->SetBranchAddress("jet_DetEta"   , &(stuple->jet_DetEta));
		tree->SetBranchAddress("jet_DetPhi"   , &(stuple->jet_DetPhi));
		tree->SetBranchAddress("jet_HadEm"    , &(stuple->jet_HadEm));
		tree->SetBranchAddress("jet_Emfr"     , &(stuple->jet_Emfr));
		tree->SetBranchAddress("jet_Ntowers"  , &(stuple->jet_Ntowers));
		tree->SetBranchAddress("jet_Ntracks"  , &(stuple->jet_Ntracks));
		tree->SetBranchAddress("jet_SeedIPhi" , &(stuple->jet_SeedIPhi));
		tree->SetBranchAddress("jet_SeedIEta" , &(stuple->jet_SeedIEta));
		tree->SetBranchAddress("jet_EmTime"   , &(stuple->jet_EmTime));

		// RAW JET COLLECTIONS
		tree->SetBranchAddress("jet_raw_num"   , &(stuple->jet_raw_num));
		tree->SetBranchAddress("jet_raw_Index" , &(stuple->jet_raw_Index));
		tree->SetBranchAddress("jet_raw_Pt"    , &(stuple->jet_raw_Pt));
		tree->SetBranchAddress("jet_raw_E"     , &(stuple->jet_raw_E));
		tree->SetBranchAddress("jet_raw_Px"    , &(stuple->jet_raw_Px));
		tree->SetBranchAddress("jet_raw_Py"    , &(stuple->jet_raw_Py));
		tree->SetBranchAddress("jet_raw_Pz"    , &(stuple->jet_raw_Pz));


		
			// JES UP JET COLLECTION
		tree->SetBranchAddress("jet_up_num"     , &(stuple->jet_up_num));
		tree->SetBranchAddress("jet_up_NJet15"  , &(stuple->jet_up_NJet15));
		tree->SetBranchAddress("jet_up_Index"   , &(stuple->jet_up_Index));
		tree->SetBranchAddress("jet_up_Pt"      , &(stuple->jet_up_Pt));
		tree->SetBranchAddress("jet_up_E"       , &(stuple->jet_up_E));
		tree->SetBranchAddress("jet_up_Px"      , &(stuple->jet_up_Px));
		tree->SetBranchAddress("jet_up_Py"      , &(stuple->jet_up_Py));
		tree->SetBranchAddress("jet_up_Pz"      , &(stuple->jet_up_Pz));
		tree->SetBranchAddress("jet_up_DetEta"  , &(stuple->jet_up_DetEta));
		tree->SetBranchAddress("jet_up_DetPhi"  , &(stuple->jet_up_DetPhi));
		tree->SetBranchAddress("jet_up_HadEm"   , &(stuple->jet_up_HadEm));
		tree->SetBranchAddress("jet_up_Emfr"    , &(stuple->jet_up_Emfr));
		tree->SetBranchAddress("jet_up_Ntowers" , &(stuple->jet_up_Ntowers));
		tree->SetBranchAddress("jet_up_Ntracks" , &(stuple->jet_up_Ntracks));
		tree->SetBranchAddress("jet_up_SeedIPhi", &(stuple->jet_up_SeedIPhi));
		tree->SetBranchAddress("jet_up_SeedIEta", &(stuple->jet_up_SeedIEta));
		tree->SetBranchAddress("jet_up_EmTime"  , &(stuple->jet_up_EmTime));
			
			// JES DOWN JET COLLECTION
		tree->SetBranchAddress("jet_down_num"     , &(stuple->jet_down_num));
		tree->SetBranchAddress("jet_down_NJet15"  , &(stuple->jet_down_NJet15));
		tree->SetBranchAddress("jet_down_Index"   , &(stuple->jet_down_Index));
		tree->SetBranchAddress("jet_down_Pt"      , &(stuple->jet_down_Pt));
		tree->SetBranchAddress("jet_down_E"       , &(stuple->jet_down_E));
		tree->SetBranchAddress("jet_down_Px"      , &(stuple->jet_down_Px));
		tree->SetBranchAddress("jet_down_Py"      , &(stuple->jet_down_Py));
		tree->SetBranchAddress("jet_down_Pz"      , &(stuple->jet_down_Pz));
		tree->SetBranchAddress("jet_down_DetEta"  , &(stuple->jet_down_DetEta));
		tree->SetBranchAddress("jet_down_DetPhi"  , &(stuple->jet_down_DetPhi));
		tree->SetBranchAddress("jet_down_HadEm"   , &(stuple->jet_down_HadEm));
		tree->SetBranchAddress("jet_down_Emfr"    , &(stuple->jet_down_Emfr));
		tree->SetBranchAddress("jet_down_Ntowers" , &(stuple->jet_down_Ntowers));
		tree->SetBranchAddress("jet_down_Ntracks" , &(stuple->jet_down_Ntracks));
		tree->SetBranchAddress("jet_down_SeedIPhi", &(stuple->jet_down_SeedIPhi));
		tree->SetBranchAddress("jet_down_SeedIEta", &(stuple->jet_down_SeedIEta));
		tree->SetBranchAddress("jet_down_EmTime"  , &(stuple->jet_down_EmTime));
	

			// VERTEX INFO
		tree->SetBranchAddress("vtx_N"        , &(stuple->vtx_N));
		tree->SetBranchAddress("vtx_NClass12" , &(stuple->vtx_NClass12));
		tree->SetBranchAddress("vtx_z"        , &(stuple->vtx_z));
		tree->SetBranchAddress("vtx_Ntracks"  , &(stuple->vtx_Ntracks));
		tree->SetBranchAddress("vtx_SumPt"    , &(stuple->vtx_SumPt));

			// MET INFO
		tree->SetBranchAddress("met_Met"    , &(stuple->met_Met));
		tree->SetBranchAddress("met_RawMet"    , &(stuple->met_RawMet));
		tree->SetBranchAddress("met_MetX"   , &(stuple->met_MetX));
		tree->SetBranchAddress("met_MetY"   , &(stuple->met_MetY));
		tree->SetBranchAddress("met_SumEt"  , &(stuple->met_SumEt));
		tree->SetBranchAddress("met_Ht"     , &(stuple->met_Ht));
		tree->SetBranchAddress("met_MetPhi" , &(stuple->met_MetPhi));

			// GENERATED MET FROM MET MODEL
//		tree->SetBranchAddress("met_Gen_d"   , &(stuple->met_Gen_d));
//		tree->SetBranchAddress("met_Gen_m"   , &(stuple->met_Gen_m));
//		tree->SetBranchAddress("met_Gen_p"   , &(stuple->met_Gen_p));
//		tree->SetBranchAddress("met_Gen_mUn" , &(stuple->met_Gen_mUn));
//		tree->SetBranchAddress("met_Gen_pUn" , &(stuple->met_Gen_pUn));

		
			// Generator level info
		tree->SetBranchAddress("gen_elenum" , &(stuple->gen_elenum));
		tree->SetBranchAddress("gen_phonum" , &(stuple->gen_phonum));
		
		tree->SetBranchAddress("gen_MomIndex" , &(stuple->gen_MomIndex));
		tree->SetBranchAddress("gen_MomPDG"   , &(stuple->gen_MomPDG));
		tree->SetBranchAddress("gen_MomStaus" , &(stuple->gen_MomStatus));
		tree->SetBranchAddress("gen_MomEtc"   , &(stuple->gen_MomEtc));
		tree->SetBranchAddress("gen_MomE"     , &(stuple->gen_MomE));
		tree->SetBranchAddress("gen_MomPx"    , &(stuple->gen_MomPx));
		tree->SetBranchAddress("gen_MomPy"    , &(stuple->gen_MomPy));
		tree->SetBranchAddress("gen_MomPz"    , &(stuple->gen_MomPz));

		tree->SetBranchAddress("gen_ProdVtxX" , &(stuple->gen_ProdVtxX));
		tree->SetBranchAddress("gen_ProdVtxY" , &(stuple->gen_ProdVtxY));
		tree->SetBranchAddress("gen_ProdVtxZ" , &(stuple->gen_ProdVtxZ));
		tree->SetBranchAddress("gen_ProdVtxT" , &(stuple->gen_ProdVtxT));

		tree->SetBranchAddress("gen_pho_Index"  , &(stuple->gen_pho_Index));
		tree->SetBranchAddress("gen_pho_PDG"    , &(stuple->gen_pho_PDG));
		tree->SetBranchAddress("gen_pho_Status" , &(stuple->gen_pho_Status));
		tree->SetBranchAddress("gen_pho_Etc"    , &(stuple->gen_pho_Etc));
		tree->SetBranchAddress("gen_pho_E"      , &(stuple->gen_pho_E));
		tree->SetBranchAddress("gen_pho_Px"     , &(stuple->gen_pho_Px));
		tree->SetBranchAddress("gen_pho_Py"     , &(stuple->gen_pho_Py));
		tree->SetBranchAddress("gen_pho_Pz"     , &(stuple->gen_pho_Pz));

		tree->SetBranchAddress("gen_ele_Index"  , &(stuple->gen_ele_Index));
		tree->SetBranchAddress("gen_ele_PDG"    , &(stuple->gen_ele_PDG));
		tree->SetBranchAddress("gen_ele_Status" , &(stuple->gen_ele_Status));
		tree->SetBranchAddress("gen_ele_Etc"    , &(stuple->gen_ele_Etc));
		tree->SetBranchAddress("gen_ele_E"      , &(stuple->gen_ele_E));
		tree->SetBranchAddress("gen_ele_Px"     , &(stuple->gen_ele_Px));
		tree->SetBranchAddress("gen_ele_Py"     , &(stuple->gen_ele_Py));
		tree->SetBranchAddress("gen_ele_Pz"     , &(stuple->gen_ele_Pz));

	
}
*/