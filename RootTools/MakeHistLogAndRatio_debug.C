/* Makes my final plots
 * don't mess it up!
 *
 */
/*{{{*/
/*  $Id: MakeHistLogAndRatio_debug.C,v 1.1 2013/02/28 03:27:20 samantha Exp $
 *  $Log: MakeHistLogAndRatio_debug.C,v $
 *  Revision 1.1  2013/02/28 03:27:20  samantha
 *  Final commit. no checks. These were never commited to cvs.
 *
 *
 *   =============================================================================
 */
/*}}}*/
 
 
#include <iostream>
#include <sstream>
#include <string>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <TF1.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLine.h>
#include <cmath>
#include <memory>
#include <iomanip>
#include "CommonTools.hh"
#include "IOColors.hh"
#include "Rtypes.h" //ROOT color keyword. this is included in ROOT headers though.
#include "TMatrixTSym.h"
//#include "TFitResultPtr.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TDirectory.h"
//#include <stdio>

using namespace std;


/********************** GLOBAL VARIABLES ***************************/
//const std::string sPRINTFOLDER("/home/samantha/TMP/eps/g30_nominal/");
const bool bDEBUG = 1;
std::string sPRINTFOLDER("/home/samantha/TMP/eps/");
std::string sPRINTFORMAT("eps");			//generate eps/pdf/gif files as output
const bool iDRAWERROR_BREAKDOWN = 1;	//on/off error break down canvas
const bool bDO_VAR_BINNING = 1;
int iREBIN = 2;		//DEFINES THE CONST REBIN FOR EACH PLOTS

//different samples
const static int iG30      = 0; //default gamma30+jets case
const static int iG40      = 1; //default gamma40+jets case
const static int iG30MET25 = 2; //gamma30+jets+met25 case
const static int iG30MET40 = 3; //gamma30+jets+met40 case
const static int iG30EXCL1J= 4; // gamma30+1==jet case
const static std::string sG30 = "iG30";
const static std::string sG40 = "iG40";
const static std::string sG30MET25 = "iG30MET25";
const static std::string sG30MET40 = "iG30MET40";
const static std::string sG30EXCL1J = "iG30EXCL1J";

//fake photon fraction for different sample
const static float fG30_FAKE_PHO_FRAC = 0.319;  			//for g30+>=1 jets
const static float fG40_FAKE_PHO_FRAC = 0.21;  				//for g40+>=1 jets
const static float fG30MET25_FAKE_PHO_FRAC_1JET = 0.15; 	//for g30+>=1jets+met25
const static float fG30MET25_FAKE_PHO_FRAC_2JET = 0.14; 	//for g30+>=2jets+met25
const static float fG30MET40_FAKE_PHO_FRAC_1JET = 0.20; 	//for g30+>=1jets+met40
const static float fG30MET40_FAKE_PHO_FRAC_2JET = 0.13; 	//for g30+>=2jets+met40
const static float fG30EXCL1J_FAKE_PHO_FRAC = 0.313;  			//for g30+==1 jets

//systematic error for all sample are taken to this assuming the change
//of this value is small - 02-17-2010
const static float fFAKE_PHO_FRAC_SIGMA = 0.068;

//cosmic background estimates for different samples
//const static float fG30_COSMIC_EST_1JET = 110;  //original APS/ICHEP estimates
//const static float fG30_COSMIC_EST_2JET = 7;
const static float fG30_COSMIC_EST_1JET = 98; //new estimates 02-22-2010
const static float fG30_COSMIC_EST_2JET = 6;
const static float fG40_COSMIC_EST_1JET = 60;
const static float fG40_COSMIC_EST_2JET = 5;
const static float fG30MET25_COSMIC_EST_1JET = 82;
const static float fG30MET25_COSMIC_EST_2JET = 6;
const static float fG30MET40_COSMIC_EST_1JET = 63;
const static float fG30MET40_COSMIC_EST_2JET = 5;
const static float fG30EXCL1J_COSMIC_EST_1JET = 83;

//halo background estimates for different samples
const static float fG30_HALO_EST_1JET = 9;
const static float fG30_HALO_EST_2JET = 1;
const static float fG40_HALO_EST_1JET = 9;   		//needs to be updated, using g30 values
const static float fG40_HALO_EST_2JET = 1;   		//needs to be updated, using g30 values
const static float fG30MET25_HALO_EST_1JET = 9;   	//needs to be updated, using g30 values
const static float fG30MET25_HALO_EST_2JET = 1;   	//needs to be updated, using g30 values
const static float fG30MET40_HALO_EST_1JET = 9;   	//needs to be updated, using g30 values
const static float fG30MET40_HALO_EST_2JET = 1;   	//needs to be updated, using g30 values
const static float fG30EXCL1J_HALO_EST_1JET = 9;

const int iTITLE_FONT = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable) 
const int iWHITE = 10;		//ROOT color white

const static int iNOMINAL       = 0; //nominal sideband
const static int iISOADDED      = 1; //iso added sideband
const static int iISOADDEDREWGT = 2; //iso added reweighed sideband


const static float fLUM_DATA    = 2043.0;	//pb-1
//DIPHOMC crosssection after filter(12.7%) = 92.3 pb
//total number of events in the sample (before goodrun) 11,762,656
const static float fLUM_DIPHOMC = 127439.39;  //pb after filter


//crosssections
const float kFAC_MC = 1.4;
// SF = DATA_LUM * ( 1 / EWK_LUM ) * KFAC
//    = DATA_LUM * ( 1 / (TOT EVTS PROCESSED/CROSS SECTION) ) * KFAC
//    here are the values picked up from pp.102
//    DATASET   TOTAL EVTS   CROSSSECTION   LUM == TOTAL EVTS/CROSSSECTION
//    Zee       12,092,155    355pb          34062 pb-1
//    Zmm       13,755,133    355pb          38746 pb-1
//    Ztt       
const static float fCS_ZEE_INC 		= 355;	//pb-1
const static float fCS_ZMM_INC 		= 355;	//pb-1
const static float fCS_ZTT_INC_P0 	= 355;	//pb-1
const static float fCS_ZTT_INC_P1_7 = 238;	//pb-1
const static float fCS_WEN_INC 		= 1960;	//pb-1
const static float fCS_WMN_INC 		= 1960;	//pb-1
const static float fCS_WTN_INC 		= 1960;	//pb-1

//these are totoal nuber of events found in each one of the
//above EWM MC samples
const static float fNTOT_ZEE 		= 12092155;
const static float fNTOT_ZMM 		= 23755133;
const static float fNTOT_ZTT_P0 	= 9438815;
const static float fNTOT_ZTT_P1_7= 6916050;
const static float fNTOT_WEN 		= 18510028;
const static float fNTOT_WMN 		= 10166426;
const static float fNTOT_WTN 		= 6931813;


//sideband comsic estimate
//I did not have a root file with nominal sideband. these
//two values are derived using a iso-added sideband job
//it should not be that different.
const static float fG30_SIDEBAND_COSMIC_EST_1JET = 22; //new estimates 02-26-2010
const static float fG30_SIDEBAND_COSMIC_EST_2JET = 2;
const static float fG40_SIDEBAND_COSMIC_EST_1JET = 15; //new estimates 03-02-2010
const static float fG40_SIDEBAND_COSMIC_EST_2JET = 2;
const static float fG30MET25_SIDEBAND_COSMIC_EST_1JET = 21; //new estimates 02-26-2010
const static float fG30MET25_SIDEBAND_COSMIC_EST_2JET = 2;
const static float fG30EXCL1J_SIDEBAND_COSMIC_EST_1JET = 20;

const static float fG30MET40_SIDEBAND_COSMIC_EST_1JET = 18; //new estimates 03-08-2010
const static float fG30MET40_SIDEBAND_COSMIC_EST_2JET = 2;

//globle error storage
std::string sERRORS("");
std::string sWARNINGS("");

/*************** END OF GLOBAL VARIABLES ***************************/

std::string WhichSample(const int iSample)
{
	std::string sample("UNKNOWN");

	if (iSample == iG30) sample = sG30;
	else if (iSample == iG40) sample = sG40;
	else if (iSample == iG30MET25) sample = sG30MET25;
	else if (iSample == iG30MET40) sample = sG30MET40;
	else if (iSample == iG30EXCL1J) sample = sG30EXCL1J;

	return sample;
}


//for debugging only  functions ------------------------


std::string GetNewHistName(TH1* h, std::string prefix="")
{
	assert (h != NULL && "GetNewHistName:: hist null");
	std::stringstream newname;
	newname << prefix << "_" << h->GetName();
	//char str[newname.str().length()+1];
	//strcpy(str, newname.str().c_str());
	//std::cout << "str = " << str << 
	//return str;
	return newname.str();
}


void PrintFileAndHist(const TFile* f, const TH1* h, const TDirectory* dir)
{
	assert (f != NULL && "PrintFileAndHist:: file null");
	assert (h != NULL && "PrintFileAndHist:: hist null");
	assert (dir != NULL && "PrintFileAndHist:: dir null");

	std::cout << "Opened file "; f->Print();
	std::cout << "Retrived hist in dir "; gDirectory->pwd();
	std::cout << "Integral (with width) = " << h->Integral("width") << std::endl;
	h->Print("all");
	std::cout << " ==============================================" << std::endl;
}
void PrintAll(const TH1* h)
{
	assert (h != NULL && "PrintFileAndHist:: hist null");
	std::cout << "Info for hist " << h->GetName() << std::endl;
	std::cout << "Integral (with width) = " << h->Integral("width") << std::endl;
	h->Print("all");
	std::cout << " ==============================================" << std::endl;
}

void MakeTruePhoFracStudy(const TH1* data, const TH1* mcpho, const TH1* qcd)
{
	
	TH1* hdata = (TH1*) data->Clone("datacopy");
	TH1* hmc = (TH1*) mcpho->Clone("mccopy");
	TH1* hmc_copy = (TH1*) mcpho->Clone("mccopy2");
	TH1* hqcd = (TH1*) qcd->Clone("qcdcopy");

	hmc_copy->Add(hqcd);
	hmc->Divide(hmc_copy);
	
	//WilsonInterval(hmc, hmc_copy);
	
	std::cout << "integrals data/mc/qcd/sum " 
			<< hdata->Integral() << "\t" 
			<< hmc->Integral() << "\t" 
			<< hqcd->Integral() << "\t" 
			<< hmc->Integral()+ hqcd->Integral()
			<< std::endl;

	std::cout << std::setw(3) << "bin" << std::setw(7) << "loedge"
		<< std::setw(10) << "data"
		<< std::setw(10) << "mc"
		<< std::setw(10) << "qcd"
		<< std::setw(15) << "ratio(mc/data)"
		<< std::setw(20) << "ratio(mc+qcd/data)"
		<< std::endl;
	for (int bin=1; bin<= hdata->GetNbinsX(); ++bin)
	{
	std::cout << std::setw(3) << bin<< std::setw(7) << hdata->GetXaxis()->GetBinLowEdge(bin)
		<< std::setw(10) << hdata->GetBinContent(bin)
		<< std::setw(10) << hmc->GetBinContent(bin)
		<< std::setw(10) << hqcd->GetBinContent(bin)
		<< std::setw(15) << hmc->GetBinContent(bin)/hdata->GetBinContent(bin)
		<< std::setw(20) << (hmc->GetBinContent(bin)+hqcd->GetBinContent(bin))/hdata->GetBinContent(bin)
		<< std::endl;
			//float deno = (1 - fake_pho_frac_d) * hmc->GetBinContent(bin) + fake_pho_frac_d * hqcd->GetBinContent(bin);
			//if ((hmc->GetBinContent(bin) + hqcd->GetBinContent(bin)) > 0 )
			//hmc->SetBinContent(bin,  hmc->GetBinContent(bin) / (1. * (hmc->GetBinContent(bin) + hqcd->GetBinContent(bin) ) ) );
			
		//hmc->SetBinContent(bin, (1-fake_pho_frac_d) * hmc->GetBinContent(bin));
	}
	//hmc->Divide(hdata);
	hmc->SetTitle("#gamma^{E_{T}>30GeV, |#eta|<1.1} + #geq 1 Jet^{E_{T}>15GeV, |#eta|<3.0} : Photon Purity, p(E_{T}^{#gamma},#eta) = #frac{ #epsilon N_{S}}{ #epsilon N_{s} + (1- #epsilon)N_{B}};E_{T}^{#gamma} (GeV); #frac{dp}{dE_{T}^{#gamma} d#eta^{#gamma}}");
	hmc->GetYaxis()->CenterTitle(1);
	hmc->GetXaxis()->CenterTitle(1);

	new TCanvas();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptStat(0);
	hmc->SetMarkerStyle(20);
	hmc->SetMarkerSize(1);
	hmc->SetMarkerColor(kBlue);
	
	hmc->Draw("E1");
	gPad->SetEditable(0);



}


//-----------------------------------------------------------------------------
void DebugLogAndRatioPlots(const TH1* phojet, const  TH1* hist_err, const std::vector<float> vErr)
{

	assert (phojet != NULL && "phojet hist is null@DebugLogAndRatioPlots()");
	assert (hist_err != NULL && "hist_err hist is null@DebugLogAndRatioPlots()");
	assert (phojet->GetNbinsX() == (int) vErr.size() && "Error vector size does not match the number of bins in hist!");
	//assuming all hists are rebinned (same sizes)
	
	std::cout << "=====================================================================" << std::endl;
	std::cout << setw(5) << "bin" << setw(8) << "LoEdge"
						<<   setw(10) << "phojet"  
						<<   setw(10) << "phojet:err"  
						<<   setw(10) << "herr"  
						<<   setw(10) << "herr:err"  
						<<   setw(10) << "total:err"  
						<< std::endl;
	std::cout << "=====================================================================" << std::endl;
			

	for (int i=1; i <=phojet->GetNbinsX();++i)
	{
		if (phojet->GetBinContent(i) || hist_err->GetBinContent(i)
			||hist_err->GetBinError(i)	)
		{
			std::cout << setiosflags(ios::fixed) << setprecision(3) 
						<< setw(5) << i << setw(10) << phojet->GetBinLowEdge(i) 
						<< setw(10) << phojet->GetBinContent(i) 
						<< setw(10)	<< phojet->GetBinError(i)
						<< setw(10) << hist_err->GetBinContent(i) 
						<< setw(10)	<< hist_err->GetBinError(i)
						<< setw(10) << vErr.at(i-1)
						<< std::endl;
			
		}
	}
	std::cout << "=====================================================================" << std::endl;

}



void DebugTwoHistRatio(const TH1* hist1, const TH1* hist2)
{
	assert ((hist1 != NULL &&  hist2 != NULL) && "DebugTwoHistRatio got a NULL hist");

	TH1* h1 = (TH1*) hist1->Clone("h1copy");
	TH1* h2 = (TH1*) hist2->Clone("h2copy");
	h1->Scale(1./(1.0 * h1->Integral()));
	h2->Scale(1./(1.0 * h2->Integral()));
	
	std::cout << "bin\trat" << std::endl;
	
	for (int bin=1; bin <= h1->GetNbinsX(); ++bin)
	{
		if (h1->GetBinContent(bin))
		{
			double cval = h1->GetBinContent(bin);
			double uval = h2->GetBinContent(bin);
			//c->SetBinContent(bin, uval/cval - 1);
			//cdhistrat->SetBinContent(bin, dval/cval - 1);
			float rat = uval/cval - 1;
			std::cout << bin << "\t" << rat << std::endl;
		}

	}
	
	
	h2->Add(h1,-1);
	h2->Divide(h1);
	h2->SetTitle("#gamma sideband and #gamma MC compared (each normlaized to unity before comparison)");
	new TCanvas();
	gPad->SetGridx();
	gPad->SetGridy();
	h2->SetMarkerStyle(20);
	h2->SetMarkerColor(kBlue);
	h2->Draw("P");
	
  TLegend *leg = new TLegend (0.7,0.6,0.9,0.9);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  
  leg->AddEntry(h2,"#frac{#gamma MC}{Sideband} - 1");
  leg->Draw();
 
}


//plots the Stat err for each hist
void DebugStatErr(const TH1* zeejet, const TH1* zmmjet,
						const TH1* zttjet, const TH1* wenjet, 
						const TH1* wmnjet, const TH1* wtnjet, 
						const TH1* mcphojet,  const TH1* qcdjet, const TH1* sumbg = 0)
{


	assert ((zeejet != NULL && zmmjet != NULL && zttjet != NULL 
				&& wenjet != NULL && wmnjet != NULL && wtnjet != NULL 
				&& mcphojet != NULL)
				&& "DebugStatErr got a NULL hist");


	TH1* zee = (TH1*) zeejet->Clone("zeejet_copy");
	TH1* zmm = (TH1*) zmmjet->Clone("zmmjet_copy");
	TH1* ztt = (TH1*) zttjet->Clone("zttjet_copy");
	TH1* wen = (TH1*) wenjet->Clone("wenjet_copy");
	TH1* wmn = (TH1*) wmnjet->Clone("wmnjet_copy");
	TH1* wtn = (TH1*) wtnjet->Clone("wtnjet_copy");
	TH1* phomc = (TH1*) mcphojet->Clone("phomcjet_copy");
	TH1* qcd = (TH1*) qcdjet->Clone("qcd_copy");

	TH1* sumBG=0;
	if (sumbg)
	{
		sumBG = (TH1*) sumbg->Clone("hist_err_copy");
	}

	//assuming all hists are rebinned (same sizes)
	for (int i=1; i <=zee->GetNbinsX();++i)
	{
		if (zee->GetBinContent(i))
		{
			zee->SetBinContent(i, zee->GetBinError(i)/zee->GetBinContent(i));
			zee->SetBinError(i,0);
		}
		if (zmm->GetBinContent(i))
		{
			zmm->SetBinContent(i, zmm->GetBinError(i)/zmm->GetBinContent(i));
			zmm->SetBinError(i,0);
		}
		if (ztt->GetBinContent(i))
		{
			ztt->SetBinContent(i, ztt->GetBinError(i)/ztt->GetBinContent(i));
			ztt->SetBinError(i,0);
		}
		if (wen->GetBinContent(i))
		{
			wen->SetBinContent(i, wen->GetBinError(i)/wen->GetBinContent(i));
			wen->SetBinError(i,0);
		}
		if (wmn->GetBinContent(i))
		{
			wmn->SetBinContent(i, wmn->GetBinError(i)/wmn->GetBinContent(i));
			wmn->SetBinError(i,0);
		}
		if (wtn->GetBinContent(i))
		{
			wtn->SetBinContent(i, wtn->GetBinError(i)/wtn->GetBinContent(i));
			wtn->SetBinError(i,0);
		}
		if (phomc->GetBinContent(i))
		{
			phomc->SetBinContent(i, phomc->GetBinError(i)/phomc->GetBinContent(i));
			phomc->SetBinError(i,0);
		}
		if (qcd->GetBinContent(i))
		{
			qcd->SetBinContent(i, qcd->GetBinError(i)/qcd->GetBinContent(i));
			qcd->SetBinError(i,0);
		}
		if (sumBG)
		{
			if (sumBG->GetBinContent(i))
			{
				sumBG->SetBinContent(i, sumBG->GetBinError(i)/sumBG->GetBinContent(i));
				sumBG->SetBinError(i,0);
			}
		}
	}

	zee->SetMarkerStyle(20);
	zmm->SetMarkerStyle(22);
	ztt->SetMarkerStyle(23);
	wen->SetMarkerStyle(20);
	wmn->SetMarkerStyle(22);
	wtn->SetMarkerStyle(23);
	phomc->SetMarkerStyle(21);
	qcd->SetMarkerStyle(21);
	if (sumBG)
	{
		sumBG->SetMarkerStyle(21);
		sumBG->SetMarkerColor(kCyan);
	}

	zee->SetMarkerColor(kRed);
	zmm->SetMarkerColor(kRed);
	ztt->SetMarkerColor(kRed);
	wen->SetMarkerColor(kBlue);
	wmn->SetMarkerColor(kBlue);
	wtn->SetMarkerColor(kBlue);
	phomc->SetMarkerColor(kGreen);
	qcd->SetMarkerColor(kYellow);
	
	zee->SetTitle("Relative Statistical Error in MC samples before scaling");

  TLegend *leg = new TLegend (0.7,0.6,0.9,0.9);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  
  leg->AddEntry(zee,"Zee");
  leg->AddEntry(zmm,"Zmm");
  leg->AddEntry(ztt,"Ztt");
  leg->AddEntry(wen,"Wen");
  leg->AddEntry(wmn,"Wmn");
  leg->AddEntry(wtn,"Wtn");
  leg->AddEntry(phomc,"#gamma MC");
  leg->AddEntry(qcd,"QCD");
  if (sumBG) leg->AddEntry(sumBG,"Tot. from Sum BG");
	
	new TCanvas();
	gPad->SetGridx();
	gPad->SetGridy();
	zee->Draw("P");
	zmm->Draw("P SAME");
	ztt->Draw("P SAME");
	wen->Draw("P SAME");
	wmn->Draw("P SAME");
	wtn->Draw("P SAME");
	phomc->Draw("P SAME");
	qcd->Draw("P SAME");
  	if (sumBG) sumBG->Draw("P SAME");
	leg->Draw();	

}

//dumps only one bin content of a hist
void DumpHistBin(const TH1* hist, const int bin=0, const bool header=true
		,const int precision=1)
{
	assert (hist != NULL && "DumpHistBin:: got null pointer!");
	assert (hist->GetDimension() == 1 && "DumpHistBin:: hist is not 1D!");

	if (header)
	{
		std::cout << setw(8) << "bin"<< setw(16) << "lo/up edges"<<  setw(13) << "val" << setw(13) << "err" << std::endl; 
	}
	
	if (bin>=0 && bin <= hist->GetNbinsX())
	{
			std::cout << setiosflags(ios::fixed) << setprecision(precision) 
					<< setw(8) << bin << setw(8) << "[" <<  hist->GetBinLowEdge(bin) 
					<< "," <<setw(4) << hist->GetXaxis()->GetBinUpEdge(bin) << "]"
					<< setw(13) << hist->GetBinContent(bin) << ", " << setw(13)
					<< hist->GetBinError(bin) << "\t" << hist->GetName() << std::endl;
	}
}

//dumps the all the bin contents of a hist
void Debug_DumpHist(TH1* hist, const int iCallerLine, const int Nbins = 0)
{
	assert (hist != NULL && "Debug_DumpHist:: got null pointer!");
	assert (( hist->GetNbinsX() >= Nbins && Nbins>=0) 
				&& "Debug_DumpHist:: reuested invalid number of Nbins");


	std::cout << __FUNCTION__<<  "::Caller Line#" << iCallerLine;
	if (Nbins>=0)
	{
		std::cout << " - First " << Nbins << " non-zero bins dump for :" << std::endl;
	} else {
		std::cout << " - All non-zero bins dump for :" << std::endl;
	}
	hist->Print();
	std::cout << setw(8) << "bin"<< setw(16) << "lo/up edges"<<  setw(10) 
				 << "val" << setw(10) << "err" << std::endl; 


	for  (int bin=0; bin <= hist->GetNbinsX(); ++bin)
	{
		if (hist->GetBinContent(bin))
		{
			std::cout << setiosflags(ios::fixed) << setprecision(1) 
				<< setw(8) << bin << setw(8) << hist->GetBinLowEdge(bin) 
				<< setw(1) << ", " << setw(4) << hist->GetXaxis()->GetBinUpEdge(bin)
				<< setw(12) << hist->GetBinContent(bin) << ", " << setw(10)
				<< hist->GetBinError(bin) << std::endl;
		}
		//dump all bin if Nbins==0 (default)
		if (bin>=Nbins && Nbins>0) break;
	}
}


//-----------------------------------------------------------------------------
void QuickDebug(const TH1* halojet, const  TH1* zeejet, const TH1* zmmjet,
					const TH1* zttjet, const TH1* wenjet, const 	TH1* wmnjet,
					const TH1* wtnjet, const  TH1* cosmicjet, const  TH1* qcdjet,
					const TH1* mcphojet, const  TH1* qcdjet_100, const  TH1* phojet)
//QuickDebug( halojet,  zeejet, zmmjet, zttjet, wenjet, wmnjet,
					//wtnjet,  cosmicjet,  qcdjet,  mcphojet,  qcdjet_100,  phojet);
{


	assert (halojet != NULL && "halojet hist is null@QuickDebug()");
	//assuming all hists are rebinned (same sizes)
	for (int i=1; i <=halojet->GetNbinsX();++i)
	{
		if (phojet->GetBinContent(i))
		{
			std::cout << "bin " << i << "[" << halojet->GetBinLowEdge(i) << ", " << 
						halojet->GetXaxis()->GetBinUpEdge(i) << "]" << std::endl;
			std::cout <<   " pho jet  [" << phojet->GetBinContent(i) << " / " 
												<< phojet->GetBinError(i) << " / " 
												<< sqrt(phojet->GetBinContent(i))
						<< "]\n qcd100   [" << qcdjet_100->GetBinContent(i) << " / " 
												<< qcdjet_100->GetBinError(i)   << " / " 
												<< sqrt(qcdjet_100->GetBinContent(i))
						<< "]\n mcpho    [" << mcphojet->GetBinContent(i) << " / " 
												<< mcphojet->GetBinError(i)   << " / " 
												<< sqrt(mcphojet->GetBinContent(i))
						<< "]\n qcdjet   [" << qcdjet->GetBinContent(i) << " / " 
												<< qcdjet->GetBinError(i)   << " / " 
												<< sqrt(qcdjet->GetBinContent(i))
						<< "]\n cosmicjet[" << cosmicjet->GetBinContent(i) << " / " 
												<< cosmicjet->GetBinError(i)   << " / " 
												<< sqrt(cosmicjet->GetBinContent(i))
						<< "]\n wtnjet   [" << wtnjet->GetBinContent(i) << " / " 
												<< wtnjet->GetBinError(i)   << " / " 
												<< sqrt(wtnjet->GetBinContent(i))
						<< "]\n wmnjet   [" << wmnjet->GetBinContent(i) << " / " 
												<< wmnjet->GetBinError(i)   << " / " 
												<< sqrt(wmnjet->GetBinContent(i))
						<< "]\n wenjet   [" << wenjet->GetBinContent(i) << " / " 
												<< wenjet->GetBinError(i)   << " / " 
												<< sqrt(wenjet->GetBinContent(i))
						<< "]\n zttjet   [" << zttjet->GetBinContent(i) << " / " 
												<< zttjet->GetBinError(i)   << " / " 
												<< sqrt(zttjet->GetBinContent(i))
						<< "]\n zmmjet   [" << zmmjet->GetBinContent(i) << " / " 
												<< zmmjet->GetBinError(i)   << " / " 
												<< sqrt(zmmjet->GetBinContent(i))
						<< "]\n zeejet   [" << zeejet->GetBinContent(i) << " / " 
												<< zeejet->GetBinError(i)   << " / " 
												<< sqrt(zeejet->GetBinContent(i))
						<< "]\n halojet  [" << halojet->GetBinContent(i) << " / " 
												<< halojet->GetBinError(i)   << " / " 
												<< sqrt(halojet->GetBinContent(i))
						<< "]\n" << std::endl;  
			
			if (i == 9) break;
		}
	}


	assert (phojet != NULL && "cphohist is null");
	
	TH1* chist = (TH1*) halojet->Clone("halojet_copy");
	TH1* cphojet = (TH1*) phojet->Clone("phojet_copy");

	assert (chist != NULL && "chist null!");
	assert (cphojet != NULL && "cphohist is null");
	chist->Add(zeejet,1);
	chist->Add(zmmjet,1);
	chist->Add(zttjet,1);
	chist->Add(wenjet,1);
	chist->Add(wmnjet,1);
	chist->Add(wtnjet,1);
	chist->Add(cosmicjet,1);
	chist->Add(qcdjet,1);
	chist->Add(mcphojet,1);

	cphojet->SetLineColor(kBlue);

	new TCanvas();
	gPad->SetLogy();
	cphojet->Draw();
	chist->Draw("same");
	std::cout << "RETURNING FROM " << __FUNCTION__ << std::endl; 
}

//-----------------------------------------------------------------------------
void DebugQcdMcMix(const TH1* hist_varyMcQcdMix,const TH1* hist_err)
{
	assert (hist_varyMcQcdMix!=NULL && "qcdMcMix hist null");
	assert (hist_err!=NULL && "hist_err hist null");
	
	TH1* hist_mixcp = (TH1*) hist_varyMcQcdMix->Clone("mix_copy");
	TH1* hist_errcp = (TH1*) hist_err->Clone("err_copy");

	
	new TCanvas();
	gPad->SetLogy();
	hist_mixcp->SetTitle("QCD/MC mix DEBUG plot");	
	hist_mixcp->SetMarkerColor(kBlue);
	hist_mixcp->SetLineColor(kBlue);
	hist_errcp->SetLineColor(kRed);
	hist_errcp->SetMarkerColor(kRed);
	hist_mixcp->Draw();
	hist_errcp->Draw("same");
/*	
	for (unsigned bin = 0; bin <= hist_mixcp->GetNbinsX() + 1; ++ bin)
	{
		float value = hist_mixcp->GetBinContent (bin);
		float error = hist_mixcp->GetBinError (bin);
		hist_mixcp->SetBinError (bin, value ? error / value : 0);
		hist_mixcp->SetBinContent (bin, 0);
	};


	
	for (unsigned bin = 0; bin <= hist_errcp->GetNbinsX() + 1; ++ bin)
	{
		const float val = hist_mixcp->GetBinContent (bin);
		const float scale = val ? 1. / val : 0;
		hist_errcp->SetBinContent (bin, (hist_errcp->GetBinContent (bin) - val) * scale);
		hist_errcp->SetBinError (bin, hist_errcp->GetBinError (bin) * scale);
	};


	assert(hist_errcp!=NULL && "hist_err_copy is null");
	new TCanvas();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	hist_errcp->SetFillColor(kRed);
	hist_errcp->SetFillStyle(3002);
	hist_mixcp->SetTitle("QCD/MC mix DEBUG plot");	
	hist_errcp->Draw("same E2");	
*/
}


//-----------------------------------------------------------------------------
// assumes all three hists have a the same bin sizes
void DebugJES(const TH1* hist, const TH1* jesup, const TH1* jesdown,
				 const int iLINE=0, const std::string title="",
				 const std::string filename = "")
{

	assert(hist != NULL && "hist_err null!.");
	assert(jesup != NULL && "jesup null!.");
	assert(jesdown != NULL && "jesdown null!.");

	TH1* chist = (TH1*) hist->Clone("hist_copy");
	TH1* cuhistrat = (TH1*) hist->Clone("hist_copy_uratio");
	TH1* cdhistrat = (TH1*) hist->Clone("hist_copy_dratio");
	TH1* uhist = (TH1*) jesup->Clone("jesup_copy");
	TH1* dhist = (TH1*) jesdown->Clone("jesdown_copy");

	TCanvas *c1 = new TCanvas();
	unsigned int w = 1200;
	unsigned int h = 600;
	c1->SetCanvasSize(w,h);
	c1->SetWindowSize(w + (w - c1->GetWw()) + 20, h + (h - c1->GetWh())+30);
	c1->Divide(2,1);
	c1->cd(1);
	//gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	chist->SetLineColor(kBlue);
	chist->SetMarkerColor(kBlue);
	chist->SetMarkerStyle (20);
	chist->SetMarkerSize(0.8);
	uhist->SetLineColor(kRed);
	uhist->SetMarkerColor(kRed);
	uhist->SetMarkerStyle (22);
	dhist->SetMarkerStyle (23);
	dhist->SetMarkerColor(kGreen);
	chist->SetMinimum(0.05);
	if (title.length()) chist->SetTitle(title.c_str());
	//find x max ann zoom-in
	double x_max = 0;
	for (int bin=chist->GetNbinsX(); bin>0; --bin)
	{
		if ( chist->GetBinContent(bin) || 
			  uhist->GetBinContent(bin) || 
			  dhist->GetBinContent(bin) )
		{
			x_max = chist->GetXaxis()->GetBinUpEdge(bin);
			break;
		}
	}
	assert (x_max > 0 && "DebugJES:: Max value for X axis retuned zero!");
	
	chist->GetXaxis()->SetRangeUser(0,x_max);
	chist->Draw("P");	
	uhist->Draw("sameP");	
	dhist->Draw("sameP");	

	TLegend *leg = new TLegend (0.7,0.72,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);

	leg->AddEntry(chist,"CENTRAL");
	leg->AddEntry(uhist,"JES UP");
	leg->AddEntry(dhist,"JES DOWN");
	//leg->AddEntry(chist,"DATA");
	//leg->AddEntry(uhist,"PHO MC");
	//leg->AddEntry(dhist,"PHO SIDEBAND");
	//leg->AddEntry(chist,"Background");
	//leg->AddEntry(chist,"FakeFrac+#sigma");
	//leg->AddEntry(uhist,"Scaled Pho MC");
	//leg->AddEntry(uhist,"Scaled Sideband");
	//leg->AddEntry(dhist,"PHO SIDEBAND");
	leg->Draw();

	//now make aratio hist of the differeces
	//zero out the two ratio hist (up/down)
	for (int bin=0;bin <= cuhistrat->GetNbinsX()+1; ++bin)
	{
			cuhistrat->SetBinContent(bin, 0);
			cuhistrat->SetBinError(bin, 0);
			cdhistrat->SetBinContent(bin, 0);
			cdhistrat->SetBinError(bin, 0);
	}
	
	
	for (int bin=1; bin <= chist->GetNbinsX(); ++bin)
	{
		if (chist->GetBinContent(bin))
		{
			double cval = chist->GetBinContent(bin);
			double uval = uhist->GetBinContent(bin);
			double dval = dhist->GetBinContent(bin);
			if ( std::isnan(cval) || std::isnan(uval) ||  std::isnan(dval))
			{
				std::cout << __FUNCTION__ << ":" << __LINE__ 
							<< ": one or bin values returned NAN!" << std::endl;
			}
			
			cuhistrat->SetBinContent(bin, uval/cval - 1);
			cdhistrat->SetBinContent(bin, dval/cval - 1);
		}

	}
	
	
	//new TCanvas();
	c1->cd(2);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();

	cuhistrat->SetMarkerStyle(22);
	cdhistrat->SetMarkerStyle(23);
	cuhistrat->SetMarkerColor(kRed);
	cdhistrat->SetMarkerColor(kGreen);
	
	cuhistrat->SetMinimum (-1.0);
	cuhistrat->SetMaximum (1.0);
	if (title.length()) cuhistrat->SetTitle(title.c_str());
	cuhistrat->GetXaxis()->SetRangeUser(0,x_max);
	cuhistrat->Draw("P");
	cdhistrat->Draw("P SAME");
	leg->Draw();

	c1->cd();
	
	if (filename.length())
	{
		c1->cd();
		c1->Print(filename.c_str());
	}

	//dump first few bins
	std::cout << __FUNCTION__<<  "::Caller Line#" << iLINE<< " - First 4 non-zero bins dump for :" << std::endl;
	std::cout << "CENTRAL - ";
	hist->Print();
	std::cout << "JES UP  - ";
	jesup->Print();
	std::cout << "JES DOWN- ";
	jesdown->Print();
	std::cout << setiosflags(ios::fixed) << setprecision(1) 
		<< setw(4) << "bin " << setw(12) << "lo/up edges" <<  setw(11) << "*central*" 
		<< setw(19) << "*jes up[diff]*" <<  setw(24) << "*jes down [diff]*"
		<< std::endl;
	int dumps = 0;
	for (int bin=1;bin <= chist->GetNbinsX()+1; ++bin)
	{
		if (chist->GetBinContent(bin))
		{
			double cval = chist->GetBinContent(bin);
			double uval = uhist->GetBinContent(bin);
			double dval = dhist->GetBinContent(bin);
			double udiff = cval - uval;
			double ddiff = cval - dval;
			std::cout << setw(4) << bin << setw(6) << hist->GetBinLowEdge(bin) 
				<< setw(1) << ", " << setw(4) << hist->GetXaxis()->GetBinUpEdge(bin)
				<< setw(13) << cval
				<< setw(12) << uval << "[ " << udiff << "]"
				<< setw(12) << dval << "[ " << ddiff << "]"
				<< std::endl;
			++dumps;
		}
		if (dumps>4) break;
	}
}

//-----------------------------------------------------------------------------
//debug only: make a error (ratio hist) with a given set of values
void DebugMakeHist(std::vector<float> contents, const TH1* hist_err ,
					std::string title="DEBUG plot")
{
	assert(hist_err != NULL && "hist is null"); 
	assert((unsigned)contents.size() == (unsigned)hist_err->GetNbinsX() 
					&& "size of contetns does not match Nbins in hist");

	TH1* hist = (TH1*) hist_err->Clone("err_copy");

	for (int i=1; i<=hist->GetNbinsX(); ++i)
	{
		hist->SetBinContent(i,contents.at(i-1));
		hist->SetBinError(i,0);
	}

	new TCanvas();
	hist->SetTitle (title.c_str());
	hist->SetMarkerStyle(kPlus);
	hist->Draw("P");
}

//end of debugging only functions  ------------------------
/*****************************************************************************/

Double_t jeteta_fitFunction(Double_t *x, Double_t *par)
{
	double val = 0;
	//val = par[0] + par[1] * sqrt(x[0]) * x[0] + par[2] / x[0];

  //if(x[0]>0.0) val=(par[0]+par[1]*sqrt(x[0]))/x[0]+par[2];
  
	//if(x[0]>0.0) val=par[0]+par[1]*x[0]*x[0]* sin(x[0]);
	//val=par[0]+par[1]*x[0]*x[0]* sin(x[0]);
	val=par[0]+par[1]*x[0]*x[0];
	return val;
}


Double_t fitFunction(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
		//val = (par[0]+par[1]*sqrt(x[0]))/x[0]+par[2];
		val = (par[0]+par[1]*sqrt(x[0]))+par[2]*x[0];
	}
	return val;
}


/*******************************************************
 * This method is to check if there are statistically 
 * significant negative bins after subtractin off 
 * some components from a hist.
 * if found , issue a error!
 *******************************************************/
void CheckNegativeBin(const TH1* hist)
{
	assert (hist != NULL && "CheckNegativeBin:: hist give is null");
	assert (hist->GetDimension() == 1 && "CheckNegativeBin:: hist dimension not 1!");
	
	for (int bin=1; bin <= hist->GetNbinsX(); ++bin)
	{
		if (hist->GetBinContent(bin)>0.0) continue;
		//now check if this bin value is statistically
		//within 0

		if ((hist->GetBinContent(bin) + hist->GetBinError(bin))>=0) continue;

		std::cout <<  __FUNCTION__ << ":" << __LINE__ 
			<< ": STATISTICALLY SIGNIFICANT NEGATIVE BIN VALUE FOUND! PLEAS CHECK!" 
			 << "bin/val/err = " << bin << "\t" << hist->GetBinContent(bin)
			 << ", " << hist->GetBinError(bin) <<  std::endl;
	}
	

}

/*********************************************************
 * I am using this method to generate a reweight function
 * to match the sideband photon Et distribution to the
 * total background (pho sideband+ pho MC) photon Et
 * distribution. Call this method within the main
 * MakeHistLogAndRatio_debug() method soon after rebinning and
 * fake photon fraction is defined.
 *
 * TO BE SAFE: add a return statement soon after 
 * this method is called.  - 03-09-2010
 *********************************************************/
void MakePhoEtWeightsForSideband(const TH1* hqcd, const TH1* hmcpho, 
						const float fake_frac, const float fake_frac_sigma, 
						const std::string brname, const int njets)
{

	assert (hqcd != NULL && "MakePhoEtWeightsForSideband:: qcd hist is null!");
	assert (hqcd->GetDimension() == 1 && "MakePhoEtWeightsForSideband:: qcd hist is not 1D!");
	assert (hmcpho != NULL && "MakePhoEtWeightsForSideband:: mcpho hist is null!");
	assert (hmcpho->GetDimension() == 1 && "MakePhoEtWeightsForSideband:: mcpho hist is not 1D!");
	assert ( (fake_frac>0 && fake_frac < 1.) && "MakePhoEtWeightsForSideband:: not a valid fake photon fraction!");

	cout << __FUNCTION__ << ": Settings passed in:" << endl;
	cout << "Fake photon fraction = " << fake_frac << "+/-" << fake_frac_sigma << endl;
	
	TH1* hist_qcd = (TH1*) hqcd->Clone("hist_nominal_wghts");
	TH1* hist_qcd2 = (TH1*) hqcd->Clone("qcd_copy2");
	TH1* hist_qcd_delta = (TH1*) hqcd->Clone("hist_plusSigma_wghts");
	TH1* hist_qcd_delta2 = (TH1*) hqcd->Clone("qcd_copy4");
	TH1* residual = dynamic_cast<TH1*> (hqcd->Clone("residual"));
	TH1* hist_mcpho = (TH1*) hmcpho->Clone("mcpho_copy");
	TH1* hist_mcpho_delta = (TH1*) hmcpho->Clone("mcpho_copy2");
	hist_qcd->Print();
	hist_mcpho->Print();

	new TCanvas();
	hist_qcd->DrawCopy();
	hist_mcpho->DrawCopy("same");
	//return;

	
	//sumw2 is already called for original hists

	hist_qcd->Scale(fake_frac / (double) hist_qcd->Integral());
	hist_mcpho->Scale((1 - fake_frac) / (double) hist_mcpho->Integral());
	hist_qcd->Add(hist_mcpho);

	//normalize back to orignal
	hist_qcd->Scale(hist_qcd2->Integral()/(double)hist_qcd->Integral());
	hist_qcd->Divide(hist_qcd2);

	std::stringstream basetitle;
	//basetitle << "#gamma^{E_{T}>30GeV,|#eta|<1.1}+ #geq " << njets << " Jet^{E_{T}>15GeV, |#eta|<3.2}: E_{T}^{#gamma} : ";
	basetitle << "#gamma^{E_{T}>30GeV,|#eta|<1.1}+ == " << njets << " Jet^{E_{T}>15GeV, |#eta|<3.2}: E_{T}^{#gamma} : ";
	std::stringstream htitle;
	htitle << basetitle.str() << "True-#gamma fraction = #epsilon =" << 1-fake_frac;
	
	hist_qcd->SetTitle(htitle.str().c_str());
	hist_qcd->GetYaxis()->SetTitle("(#epsilon  #gamma-MC + (1 - #epsilon)Sideband )/Sideband");

	std::stringstream c1_name;
	//c1_name << "Fit_" << brname<<"_error";
	c1_name << "Fit_" << brname;
	TCanvas *c1 = new TCanvas(c1_name.str().c_str(),c1_name.str().c_str(),1000,800);
	//new TCanvas();
	c1->Divide(1,2);
	c1->cd(1);

	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptFit(1);
	//hist_qcd->SetStats(0);

	hist_qcd->SetLineColor(kBlue);
	//hist_qcd->SetMarkerColor(kBlue);
	hist_qcd->SetMarkerStyle (20);
	hist_qcd->SetMarkerSize(1);


	
	//find fit range
	int xminbin=0, xmaxbin=0;
	FindNonZeroXrange(hist_qcd,xminbin,xmaxbin);
	std::cout << "minx, maxx = " << xminbin << "\t" << xmaxbin << std::endl;
	TF1 *f4 = new TF1("fit_func4",fitFunction,hist_qcd->GetBinCenter(xminbin),hist_qcd->GetXaxis()->GetBinCenter(xmaxbin),3);
	f4->SetParameter(0,0.5);
	f4->SetParameter(1,0.5);
	f4->SetParameter(2,30);
	f4->SetLineColor(8);

	hist_qcd->Fit(f4,"R+");

	TList *funcList = hist_qcd->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1* fitfunc = (TF1*) it->Next();
	if (fitfunc == NULL)
	{
		std::cout <<  __FUNCTION__ << ":" << __LINE__  
			<< "::Could not retrieve fit function! returning !" <<  std::endl;
	}
	fitfunc->SetNameTitle(brname.c_str(),htitle.str().c_str());
	//fitfunc->SetNameTitle("tf1_leadjeteta_nominalweights", htitle.str().c_str());

	//residual of the fit
	residual->GetYaxis()->SetTitle("#frac{Fit -- Val}{Error}");
	for (int bin=1; bin < residual->GetNbinsX(); ++bin)
	{
		float val = hist_qcd->GetBinContent(bin);

		float err = hist_qcd->GetBinError(bin);
		float fit = fitfunc->Eval(residual->GetBinCenter(bin));
		if (err>0) residual->SetBinContent(bin,(fit-val)/err);
		residual->SetBinError(bin,0);
	}

/*	c1->cd(2);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	residual->Draw();
	c1->cd();
*/

	
/******************************************/
//uncertainty for this method
//derive another set of wgt with fake_frac+/-sigma
/******************************************/
	// this is to derive the syst error for the above see elog#1520
	// delta w = (MC/SB - 1)  * pho_frac_sigma
	//hist_mcpho_delta->Scale(fake_frac_sigma/(double)hist_mcpho_delta->Integral());
	//hist_qcd_delta->Scale(fake_frac_sigma/(double)hist_qcd_delta->Integral());
	//hist_mcpho_delta->Add(hist_qcd_delta, -1);
	//hist_mcpho_delta->Scale(hist_qcd_delta2->Integral()/(double)hist_mcpho_delta->Integral());
	//hist_mcpho_delta->Divide(hist_qcd_delta2);
	//new TCanvas();
	//hist_mcpho_delta->DrawClone();
	hist_qcd_delta->Scale( (fake_frac +fake_frac_sigma) / (double) hist_qcd_delta->Integral());
	hist_mcpho_delta->Scale((1 - (fake_frac + fake_frac_sigma)) / (double) hist_mcpho_delta->Integral());
	hist_qcd_delta->Add(hist_mcpho_delta);

	//normalize back to orignal
	hist_qcd_delta->Scale(hist_qcd_delta2->Integral()/(double)hist_qcd_delta->Integral());
	hist_qcd_delta->Divide(hist_qcd_delta2);

	
	std::stringstream htitle_delta;
	htitle_delta << basetitle.str() << " : True-#gamma fraction = #epsilon + #sigma = " << 1 - fake_frac + fake_frac_sigma;
	
	hist_qcd_delta->SetTitle(htitle_delta.str().c_str());
	hist_qcd_delta->GetYaxis()->SetTitle("( (#epsilon+#sigma) #gamma-MC + (1 - (#epsilon+#sigma)) Sideband )/Sideband");
	hist_qcd_delta->GetYaxis()->SetTitleSize(0.04);

/*	std::stringstream c2_name;
	c2_name << "Fit_" << brname << "_plusSigma";
	TCanvas *c2 = new TCanvas(c2_name.str().c_str(),c2_name.str().c_str(),1000,800);
	c2->Divide(1,2);
	c2->cd(1);

	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptFit(1);
	*/
	//hist_qcd->SetStats(0);

	hist_qcd_delta->SetLineColor(kBlue);
	hist_qcd_delta->SetMarkerStyle (20);
	hist_qcd_delta->SetMarkerSize(1);

	//find fit range
	int xminbin_delta=0, xmaxbin_delta=0;
	FindNonZeroXrange(hist_qcd_delta,xminbin_delta,xmaxbin_delta);
	std::cout << "minx_delta, maxx_delta = " << xminbin_delta << "\t" << xmaxbin_delta << std::endl;
	TF1 *f4p = new TF1("fit_func4_delta",fitFunction,hist_qcd_delta->GetBinCenter(xminbin_delta),hist_qcd_delta->GetXaxis()->GetBinCenter(xmaxbin_delta),3);
	f4p->SetParameter(0,0.5);
	f4p->SetParameter(1,0.5);
	f4p->SetParameter(2,30);
	f4p->SetLineColor(8);

	c1->cd(2);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptFit(1);

	hist_qcd_delta->Fit(f4p,"R+");


	TList *funcList_delta = hist_qcd_delta->GetListOfFunctions();
	TIterator *it_delta = funcList_delta->MakeIterator();
	TF1* fitfunc_delta = (TF1*) it_delta->Next();
	std::stringstream tf1deltaname;
	tf1deltaname << brname << "_plusSigma";
	fitfunc_delta->SetNameTitle(tf1deltaname.str().c_str(),htitle_delta.str().c_str());

	//residual of the fit
/*	residual->GetYaxis()->SetTitle("#frac{Fit -- Val}{Error}");
	for (int bin=1; bin < residual->GetNbinsX(); ++bin)
	{
		float val = hist_qcd->GetBinContent(bin);
		float err = hist_qcd->GetBinError(bin);
		float fit = fitfunc->Eval(residual->GetBinCenter(bin));
		if (err>0) residual->SetBinContent(bin,(fit-val)/err);
		residual->SetBinError(bin,0);
	}

	c1->cd(2);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	residual->Draw();
*/
	c1->cd();

/******************************************/


	
	//now save all above to a hist
	
	//check if there is a saved version already.
	//if so delete that before writing the new one.
	//otherwise there will be many copied of the same obj
	//with different cycles
	//TFile f("g40jetsmet40_SidebandReweightFunctions.root","UPDATE");
	TFile f("g30jets_SidebandJetEtaReweightFunctions.root","UPDATE");
	//TFile f("SidebandWgtFunctionFromBckgOnly_PhoFracPlusSigma.root","UPDATE");
	//TFile f("SidebandWgtFunctionFromBckgOnly_PhoFracMinusSigma.root","UPDATE");
	if ( f.Get(brname.c_str()) != NULL)
	{
		std::stringstream old_tf;
		old_tf << brname << ";1";
		f.Delete(old_tf.str().c_str());
	}
	fitfunc->Write();

	if ( f.Get(fitfunc_delta->GetName()) != NULL)
	{
		std::stringstream old_tf_delta;
		old_tf_delta << fitfunc_delta->GetName() << ";1";
		f.Delete(old_tf_delta.str().c_str());
	}
	fitfunc_delta->Write();

	
	if ( f.Get(c1_name.str().c_str()) != NULL)
	{
		std::stringstream old_c1;
		old_c1 << c1_name.str() << ";1";
		f.Delete(old_c1.str().c_str());
	}
	c1->cd();
	c1->Write();
/*
	if ( f.Get(c2->GetName()) != NULL)
	{
		std::stringstream old_c2;
		old_c2 << c2->GetName() << ";1";
		f.Delete(old_c2.str().c_str());
	}
	c2->cd();
	c2->Write();
*/

	
	if ( f.Get(hist_qcd->GetName()) != NULL)
	{
		std::stringstream old_th1f;
		old_th1f << hist_qcd->GetName() << ";1";
		f.Delete(old_th1f.str().c_str());
	}
	hist_qcd->Write();

	if ( f.Get(hist_qcd_delta->GetName()) != NULL)
	{
		std::stringstream old_th1f_delta;
		old_th1f_delta << hist_qcd_delta->GetName() << ";1";
		f.Delete(old_th1f_delta.str().c_str());
	}
	hist_qcd_delta->Write();

	
	f.ls();
	f.Close();

}

/*********************************************************
 * This s to generate additional weights to help the 
 * photo Et based sideband reweighing.
 *	This compares lead jet det eta in photon MC and 
 *	data sideband.
 * TO BE SAFE: add a return statement soon after 
 * this method is called.  - 03-19-2010
 *********************************************************/
void MakeJetEtaWeightsForSideband(const int qcdMixMtd, const TH1* hqcd, const TH1* hphodata, 
						const std::string brname, const int njets)
{

	assert (qcdMixMtd  == 1 && "MakeJetEtaWeightsForSideband:: this method works only for method 1");
	assert (hqcd != NULL && "MakeJetEtaWeightsForSideband:: qcd hist is null!");
	assert (hqcd->GetDimension() == 1 && "MakeJetEtaWeightsForSideband:: qcd hist is not 1D!");
	assert (hphodata != NULL && "MakeJetEtaWeightsForSideband:: phodata hist is null!");
	assert (hphodata->GetDimension() == 1 && "MakeJetEtaWeightsForSideband:: phodata hist is not 1D!");

	cout << __FUNCTION__ << ": Settings passed in:" << endl;
	
	TH1* hist_qcd = (TH1*) hqcd->Clone("hist_nominal_wghts");
	TH1* hist_data = (TH1*) hphodata->Clone("phodata_copy");
	hist_qcd->Print();
	hist_data->Print();

	new TCanvas();
	hist_qcd->DrawCopy();
	hist_data->DrawCopy("same");

	
	//sumw2 is already called for original hists

	hist_data->Scale(hist_qcd->Integral()/(double) hist_qcd->Integral());
	hist_data->Divide(hist_qcd);

	std::stringstream basetitle;
	basetitle << "#gamma^{E_{T}>30GeV,|#eta|<1.1}+ #geq " << njets << " Jet^{E_{T}>15GeV, |#eta|<3.2}: lead jet #eta : ";
	hist_data->SetTitle(basetitle.str().c_str());
	hist_data->GetYaxis()->SetTitle("DATA/ TOTAL BACKGROUND");

	std::stringstream c1_name;
	c1_name << "Fit_" << brname;
	TCanvas *c1 = new TCanvas(c1_name.str().c_str(),c1_name.str().c_str(),1000,800);

	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptFit(1);

	hist_data->SetLineColor(kBlue);
	hist_data->SetMarkerStyle (20);
	hist_data->SetMarkerSize(1);

	hist_data->Draw();
	
	//find fit range
	int xminbin=0, xmaxbin=0;
	FindNonZeroXrange(hist_data,xminbin,xmaxbin);
	std::cout << "minx, maxx = " << xminbin << "\t" << xmaxbin << std::endl;

	TF1 *f6 = new TF1("fit2","gaus", -1.5, 1.5);
	hist_data->Fit(f6,"R");

	TF1 *f7 = new TF1("fit3","gaus", 2.4,3.4);
	hist_data->Fit(f7,"R");

	TF1 *f8 = new TF1("fit1","gaus", -3.4,-2.4);
	hist_data->Fit(f8,"R");

	double par[9];
	f8->GetParameters(&par[0]);
	f6->GetParameters(&par[3]);
	f7->GetParameters(&par[6]);
	delete f6;
	delete f7;
	delete f8;
	TF1 *final = new TF1("final","gaus(0)+gaus(3)+gaus(6)", -3.4, 3.4);
	//TF1 *final = new TF1("final","gaus", -3.4, 3.4);
	final->SetParameters(par);
	//final->SetParameter(0,0.5);
	final->SetLineColor(kRed);
	hist_data->Fit(final,"R");

/*
 * all this is to figure out a way to get an error for this weight
 *
 * const int nPar = 9;
    double Par[20], ParMin[20], fVal, fValMin;
	 final->GetParameters(ParMin);
	 double ParErr[20];
    for(int i=0; i<nPar; i++) ParErr[i] = final->GetParError(i);
	 
    //gMinuit->Eval(nPar, 0, fValMin, ParMin,0);
	
    //for(int i=0; i<nPar; i++) {
	//	 Par[i] = ParMin[i]+ 0.1 * ParErr[i];
	 //}

  //  gMinuit->Eval(nPar, 0, fVal, Par,0);

	 
	// const float x = 1.2;
	// fValMin = final->Eval(x);
	// std::cout << " fValMin = " << fValMin << std::endl;
*/
/*	TRandom ran;
	const float step = 1./6000.;
	float x = -3.;
	
	TProfile *h1 = new TProfile("nomWgt","Nominal Weights",600,-3,3);
	TProfile *h2 = new TProfile("errWgt","Error Weights",600,-3,3);
	
	
 	while (x < 3)
	{
		final->SetParameters(ParMin);
	 	fValMin = final->Eval(x);
	 	std::cout << " fValMin[" << x << "] = ";

		fVal = 999.;
		
		for (int i=0; i<1000; ++i)
		{
			for (int npar =0 ; npar<nPar; ++npar)
			{
				Par[npar] = ParMin[npar]+(2.0*(ran.Uniform()-0.5)) * ParErr[npar];
				final->SetParameters(Par);
				fVal = final->Eval(x);
			}
			if (fVal - fValMin<1) 
			{
				std::cout << "Min val found = " << fVal  << " (" << fVal - fValMin << ")" << std::endl;
				break;
			}
		}

		if (fVal - fValMin<1)
		{
			h1->Fill(x, fValMin);
			h2->Fill(x, fVal);
		}
		x += step;
	}


	new TCanvas();
	h1->Draw();
	new TCanvas();
	h2->Draw();

*/	 

	/* this doed not work.
	
		typedef TMatrixTSym<double>	TMatrixDSym;
	TFitResultPtr r = hist_data->Fit(final,"S");
	TMatrixTDSym cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
	Double_t chi2   = r->Chi2(); // to retrieve the fit chi2
	Double_t par0   = r->Value(0); // retrieve the value for the parameter 0
	Double_t err0   = r->Error(0); // retrieve the error for the parameter 0

	std::cout << "chi2 = " << chi2 << std::endl;
	std::cout << "par0 = " << par0 << std::endl;
	std::cout << "err0 = " << err0 << std::endl;
	
	r->Print("V");     // print full information of fit including covariance matrix
	//r->Write();        // store the result in a file
 */
	
	TList *funcList = hist_data->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1* fitfunc = (TF1*) it->Next();
	assert (fitfunc != NULL && "MakeJetEtaWeightsForSideband::Cannot get Final fit function to save! it is null! pl. check");
	fitfunc->SetNameTitle("tf1_leadjeteta_nominalweights", hist_data->GetYaxis()->GetTitle());

	TFile f("g30jets_SidebandJetEtaReweightFunctions.root","UPDATE");
	if ( f.Get(brname.c_str()) != NULL)
	{
		std::stringstream old_tf;
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
	c1->cd();
	c1->Write();
	
	if ( f.Get(hist_data->GetName()) != NULL)
	{
		std::stringstream old_th1f;
		old_th1f << hist_data->GetName() << ";1";
		f.Delete(old_th1f.str().c_str());
	}
	hist_data->Write();

	f.ls();
	f.Close();

}





Double_t AlphaSsystFunc(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
		 val = (par[0]+par[1]*sqrt(x[0]))/x[0]+par[2];
	}
	return val;
}


void GetAlphaS_Systematic(const std::string brname, std::vector<float> &vErr, 
								const TH1* hist_mc)
{
	std::cout << "In " << __FUNCTION__ << std::endl;
	assert (hist_mc != NULL && "GetAlphaS_Systematic:: hist is null!");
	assert (hist_mc->GetDimension() == 1 && "GetAlphaS_Systematic:: hist is not 1D!");

	TFile alphaSyst("AlphaS_Syst_FitFunctions.root");
	assert (! alphaSyst.IsZombie() && "AlphaS_Syst file not fouund!");
	//alphaSyst.ls();


	TF1* tf1_temp = dynamic_cast<TF1*> (alphaSyst.Get(brname.c_str()));
	assert (tf1_temp != NULL && "GetAlphaS_Systematic:: Fit function not found to read the parameters!");

	const int NparNeed = 3;
	
	assert (tf1_temp->GetNdim() == 1 && "GetAlphaS_Systematic:: Read fit function is not 1D!");
	assert (tf1_temp->GetNpar() == NparNeed && "GetAlphaS_Systematic:: Read fit function NPar does not match with expected!");

	TF1 *tf1A = new TF1("AlaphaS", AlphaSsystFunc,
			hist_mc->GetBinCenter(1), hist_mc->GetBinCenter(hist_mc->GetNbinsX()), NparNeed);

	std::cout << "opened file "; alphaSyst.Print();
	double par[NparNeed];

	tf1_temp->GetParameters(par);
	for (int i =0; i <NparNeed; ++i)
	{
	//	std::cout << "par[" << i << "]=" << par[i] << std::endl;
	}
	tf1A->SetParameters(par);
	std::cout << "function "; tf1A->Print();

	std::stringstream hist_name;
	hist_name << "TH1F_"<< brname;
	TH1* hist_temp = dynamic_cast<TH1*> (alphaSyst.Get(hist_name.str().c_str()));
	assert (hist_temp != NULL && "GetAlphaS_Systematic:: Required hist not found!");
	//new TCanvas();
	//hist_temp->Draw();
	//gPad->SetEditable(0);
	

	float first_bin_val = 0;
	for (int bin = 1; bin<= hist_temp->GetNbinsX(); ++bin)
	{
		if (hist_temp->GetBinContent(bin)>0) 
		{
			first_bin_val = hist_temp->GetBinContent(bin);
			break;
		}
	}
	//delete hist_temp;
	assert (first_bin_val>0 && "GetAlphaS_Systematic:: No first_bin_val found!");
	//std::cout << "first_bin_val = " << first_bin_val << std::endl;

	if (tf1A != NULL)
	{
		vErr.clear();// removed stuffed zeros
		std::cout << std::setw(5) << "bin" << std::setw(10) << "err at bin center " << std::endl;
		for (int bin = 1; bin<= hist_mc->GetNbinsX(); ++bin)
		{
			float fAe = tf1A->Eval(hist_mc->GetBinCenter(bin)); 
			//keep minimum error at first bin val. fits are not perfect elog#1514
			if (fAe < first_bin_val) fAe = first_bin_val;
			vErr.push_back(fAe * hist_mc->GetBinContent(bin));
			//std::cout <<  bin << "[" << hist_mc->GetBinCenter(bin) << "] = " << fAe << std::endl;
			std::cout << std::setw(5) << bin << std::setw(10) << vErr.at(bin-1) << std::endl;
		}
		delete tf1A;
	} else 
	{
		std::cout << __LINE__ << ":: WARNING !  NO Alpha-S Systematic function found for "  << brname  << std::endl;
	}

}


void DrawErrorBreakDown(const int QCDerrMtd, 
		const int iSideband,
		const TH1* hist_err,
		const std::vector<float>& JESerr, const std::vector<float>&  qcdmcMixErr, 
		const std::vector<float>& vTotalEWKLumErr,
		const std::vector<float>& cosmicErr, const std::vector<float>& haloErr, 
		const std::vector<float>&  statErr, const std::vector<float>& IdErr, 
		const std::vector<float>& vPDFerror, const std::vector<float>&  vQ2error, 
		const std::vector<float>& vISRFSRerror, const std::vector<float>& vAlphaSerror, 
		const std::vector<float>& vSidebandRwtErr,
		const std::vector<float>& ERROR, const std::string title,
		const int jets, const std::string which
		)
{

	assert ( (QCDerrMtd == 1 || QCDerrMtd == 2) 
			&& "DrawErrorBreakDown:: QCDerrMtd is not valid!");
	assert ( hist_err != NULL 
			&& "DrawErrorBreakDown:: hist_err is null!");
	assert (hist_err->GetDimension() == 1 
			&& "DrawErrorBreakDown:: hist_err is is not 1D!");

	TH1::AddDirectory(kFALSE);
	TH1* hist = dynamic_cast<TH1*> (hist_err->Clone("Total Systematic"));
	
	for (int bin=0; bin <= hist->GetNbinsX()+1;++bin)
	{
		hist->SetBinContent(bin,0);
		hist->SetBinError(bin,0);
	}

	hist->SetTitle(title.c_str());

	TH1* h_jes = dynamic_cast<TH1*> (hist->Clone("JES"));
	TH1* h_phoFrac = dynamic_cast<TH1*> (hist->Clone("True Photon Fration"));
	TH1* h_EWKLum = dynamic_cast<TH1*> (hist->Clone("Total EWK Lum"));
	TH1* h_cosmic = dynamic_cast<TH1*> (hist->Clone("Cosmic"));
	TH1* h_halo = dynamic_cast<TH1*> (hist->Clone("Beam Halo"));
	TH1* h_stat = dynamic_cast<TH1*> (hist->Clone("Statistical"));
	TH1* h_qcdId = dynamic_cast<TH1*> (hist->Clone("QCD ID"));
	TH1* h_pdf = dynamic_cast<TH1*> (hist->Clone("PDF"));
	TH1* h_q2 = dynamic_cast<TH1*> (hist->Clone("Q2"));
	TH1* h_isrfsr = dynamic_cast<TH1*> (hist->Clone("ISR/FSR"));
	TH1* h_alphas = dynamic_cast<TH1*> (hist->Clone("Alpha_s"));
	TH1* h_sidebandrewgt = dynamic_cast<TH1*> (hist->Clone("Sideband Reweighing"));

	for (int bin=1; bin <= hist->GetNbinsX();++bin)
	{
		if (hist_err->GetBinContent(bin))
		{
			hist->SetBinContent(bin,ERROR.at(bin-1));
			h_jes->SetBinContent(bin,JESerr.at(bin-1));
			h_cosmic->SetBinContent(bin,cosmicErr.at(bin-1));
			h_halo->SetBinContent(bin,haloErr.at(bin-1));
			h_stat->SetBinContent(bin,statErr.at(bin-1));
			h_qcdId->SetBinContent(bin,IdErr.at(bin-1));
			h_EWKLum->SetBinContent(bin,vTotalEWKLumErr.at(bin-1));

			if (QCDerrMtd == 1 && iSideband == iISOADDEDREWGT)
			{
				h_sidebandrewgt->SetBinContent(bin,vSidebandRwtErr.at(bin-1));
			}

			//if (QCDerrMtd == 2)
			{
				h_phoFrac->SetBinContent(bin,qcdmcMixErr.at(bin-1));
				h_pdf->SetBinContent(bin,vPDFerror.at(bin-1));
				h_q2->SetBinContent(bin,vQ2error.at(bin-1));
				h_isrfsr->SetBinContent(bin,vISRFSRerror.at(bin-1));
				h_alphas->SetBinContent(bin,vAlphaSerror.at(bin-1));
			}
		}

	}



	hist->SetMarkerStyle(29); hist->SetMarkerColor(kRed);
	hist->SetMarkerSize(1.5);
	h_jes->SetMarkerStyle(20);
	h_cosmic->SetMarkerStyle(21);
	h_cosmic->SetMarkerColor(30);
	h_halo->SetMarkerStyle(22);
	h_halo->SetMarkerColor(46);

	h_stat->SetMarkerStyle(23);
	h_stat->SetMarkerColor(33);
	h_qcdId->SetMarkerStyle(20);
	h_qcdId->SetMarkerColor(38);

	h_phoFrac->SetMarkerStyle(2); h_phoFrac->SetMarkerSize(1.5);
	h_phoFrac->SetMarkerColor(7);

	h_pdf->SetMarkerStyle(20);
	h_pdf->SetMarkerColor(kBlue);
	h_q2->SetMarkerStyle(21);
	h_q2->SetMarkerColor(kGreen);
	h_isrfsr->SetMarkerStyle(3);
	h_isrfsr->SetMarkerColor(8);
	h_alphas->SetMarkerStyle(23);
	h_alphas->SetMarkerColor(6);

	h_EWKLum->SetMarkerStyle(23);
	h_EWKLum->SetMarkerColor(kRed);

	h_sidebandrewgt->SetMarkerStyle(23);
	h_sidebandrewgt->SetMarkerColor(kBlue);
	

	TLegend *leg = new TLegend (0.75,0.45,1.,1.);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->AddEntry(hist,hist->GetName());
	leg->AddEntry(h_jes,h_jes->GetName());
	leg->AddEntry(h_cosmic,h_cosmic->GetName());
	leg->AddEntry(h_halo,h_halo->GetName());
	leg->AddEntry(h_stat,h_stat->GetName());
	leg->AddEntry(h_qcdId,h_qcdId->GetName());
	leg->AddEntry(h_EWKLum,h_EWKLum->GetName());

	
	TCanvas *c2 = new TCanvas("err_log","ERRORHIST",600,800);
	c2->Divide(1,2);
	c2->cd(1);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetLogy();

	//Absolute errors, LOG plot
	hist->GetYaxis()->SetTitle("Absolute Total Error");
	hist->DrawCopy("P");
	h_jes->DrawCopy("P SAME");
	h_cosmic->DrawCopy("P SAME");
	h_halo->DrawCopy("P SAME");
	h_stat->DrawCopy("P SAME");
	h_qcdId->DrawCopy("P SAME");
	h_EWKLum->DrawCopy("P SAME");



	if (QCDerrMtd == 1 && iSideband == iISOADDEDREWGT)
	{
		h_sidebandrewgt->DrawCopy("P SAME");
		leg->AddEntry(h_sidebandrewgt,h_sidebandrewgt->GetName());
	}

	//if (QCDerrMtd == 2)
	{
		h_phoFrac->DrawCopy("P SAME");
		h_pdf->DrawCopy("P SAME");
		h_q2->DrawCopy("P SAME");
		h_isrfsr->DrawCopy("P SAME");
		h_alphas->DrawCopy("P SAME");

		
		leg->AddEntry(h_phoFrac,h_phoFrac->GetName());
		leg->AddEntry(h_pdf,h_pdf->GetName());
		leg->AddEntry(h_q2,h_q2->GetName());
		leg->AddEntry(h_isrfsr,h_isrfsr->GetName());
		leg->AddEntry(h_alphas,h_alphas->GetName());
	}

	leg->Draw();


	//TCanvas *c3 = new TCanvas("err_rat","ERROR BREAKDOWN RATIO");
	
	for (int bin=1; bin <= hist->GetNbinsX();++bin)
	{
		float val = 0;
		if (hist_err->GetBinContent(bin)>0) val = 1/hist_err->GetBinContent(bin);

		hist->SetBinContent(bin,ERROR.at(bin-1) * val);
		h_jes->SetBinContent(bin,JESerr.at(bin-1) * val);
		h_cosmic->SetBinContent(bin,cosmicErr.at(bin-1) * val);
		h_halo->SetBinContent(bin,haloErr.at(bin-1) * val);
		h_stat->SetBinContent(bin,statErr.at(bin-1) * val);
		h_qcdId->SetBinContent(bin,IdErr.at(bin-1) * val);
		h_EWKLum->SetBinContent(bin,vTotalEWKLumErr.at(bin-1) * val);

		if (QCDerrMtd == 1 && iSideband == iISOADDEDREWGT)
		{
			h_sidebandrewgt->SetBinContent(bin, vSidebandRwtErr.at(bin-1) * val);
		}

		//if (QCDerrMtd == 2)
		{
			h_phoFrac->SetBinContent(bin,qcdmcMixErr.at(bin-1) * val);
			h_pdf->SetBinContent(bin,vPDFerror.at(bin-1) * val);
			h_q2->SetBinContent(bin,vQ2error.at(bin-1) * val);
			h_isrfsr->SetBinContent(bin,vISRFSRerror.at(bin-1) * val);
			h_alphas->SetBinContent(bin,vAlphaSerror.at(bin-1) * val);
		}
	}

	
	c2->cd(2);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetGridx();
	gPad->SetGridy();

	hist->SetMaximum(1);
	hist->GetYaxis()->SetTitle("Error/Prediction");
	hist->Draw("P");
	h_jes->Draw("P SAME");
	h_cosmic->Draw("P SAME");
	h_halo->Draw("P SAME");
	h_stat->Draw("P SAME");
	h_qcdId->Draw("P SAME");
	h_EWKLum->Draw("P SAME");

	if (QCDerrMtd == 1 && iSideband == iISOADDEDREWGT)
	{
		h_sidebandrewgt->Draw("P SAME");
	}

	//if (QCDerrMtd == 2)
	{
		h_phoFrac->Draw("P SAME");
		h_pdf->Draw("P SAME");
		h_q2->Draw("P SAME");
		h_isrfsr->Draw("P SAME");
		h_alphas->Draw("P SAME");
	}

	leg->Draw();

	std::ostringstream str;
	//str << "ErrBreakdown_plot" << jets << "_" << which << ".eps";
	str << sPRINTFOLDER << "ErrBreakdown_plot" << jets << "_" << which << "." << sPRINTFORMAT.c_str();
	c2->Print (str.str().c_str());


}


void ZeroOutNegativeBins(TH1* hist)
{
	//this is hack to remove bins with val<0
	//set bin value to zero if val<0
	assert (hist != NULL && "hist is null!");
	assert (hist->GetDimension() == 1 && "not a 1D hist");

	for (int bin=1; bin <=hist->GetNbinsX(); ++bin)
	{
		if (hist->GetBinContent(bin)<0)
		{
			hist->SetBinContent(bin,0);
			//hist->SetBinError(bin,0);
		}
	}
}

//-----------------------------------------------------------------------------
void GetCosmicErr(std::vector<float>& ErrVec, const TH1 *hist , const bool debug)
{
	if (debug) std::cout << "IN COSMIC ERROR"<< std::endl;
	assert (hist!=NULL && "GetComicErr::hist is null");
	ErrVec.clear();
	for (int bin=1; bin <= hist->GetNbinsX(); bin++) {
		float val = hist->GetBinContent(bin);
		float err = 0;
		//if (val) err = 1 / sqrt(val);		//take stat error as the syst
		//if (val) err = 1 ;		//????????????????? for now 04-15-2010
		err = hist->GetBinError(bin);
		ErrVec.push_back(err);
		std::cout << "cosmic bin[" << bin << "]=" << val  << " err = " << err << "[" << hist->GetBinError(bin) << "]" << std::endl;
		if (debug) std::cout << "cosmic bin[" << bin << "]=" << val << std::endl;
	}
}

//-----------------------------------------------------------------------------
void GetHaloErr(std::vector<float>& ErrVec, const TH1 *hist, bool debug=false)
{
	assert(hist!=NULL && "hist is null!");
	ErrVec.clear();
	for (int i=1; i <= hist->GetNbinsX(); i++) {
		float bin = hist->GetBinContent(i);
		float err = bin * 0.5;		//take 50% to be the syst
		ErrVec.push_back(err);
		if (debug) std::cout << "halo err [" << i << "]="<< err <<std::endl;
	}
	if (debug)
	{
		std::cout << "halo err arr size = " << ErrVec.size() << std::endl;
		std::cout << "\n\n\n i am out" << std::endl;
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
	if (debug) std::cout << "Done creating var bin vector. " << std::endl;
	
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


//-----------------------------------------------------------------------------
//this make the ratio plot and returns the error function (TF1) used for the
//systematics band to the log plot. so the log plot errors will be based
//on this same error function.
//-----------------------------------------------------------------------------
void MakeRatioPlot(const std::string which,  const int jets, const TH1* phojet,
						const TH1* hist_err, std::vector<float> Errors, TPad *pad=0, bool debug=false)
{
	assert(Errors.size() == (unsigned int) hist_err->GetNbinsX() && "error and bin number does not match");
	assert(phojet!=NULL && "phojet hist null");
	assert(hist_err!=NULL && "hist_err hist null");

	TH1* pho = (TH1*) phojet->Clone("phojet_copy");
	TH1* err = (TH1*) hist_err->Clone("hist_err_copy");
	
	if (debug)
	{
		std::cout << "\n ======= "<< __FUNCTION__ << ":" << __LINE__ << ":" << std::endl;
		std::cout << " ======= Hists given to MakeRatioHist first 4 bins" << std::endl;
		Debug_DumpHist(pho, __LINE__,4);
		Debug_DumpHist(err, __LINE__,4);
	}
	
	for (unsigned bin = 1; bin <= (unsigned)err->GetNbinsX(); ++ bin)
	{
		const float val = err->GetBinContent (bin);
		const float scale = val ? 1. / val : 0;
		pho->SetBinContent (bin, (pho->GetBinContent (bin) - val) * scale);
		pho->SetBinError (bin, pho->GetBinError (bin) * scale);
	};


	//std::cout << "\n ======= "<< __FUNCTION__ << ":" << __LINE__ << ":" << std::endl;
	//std::cout << " ======= Err hist val and ERROR vector values" << std::endl;
	//std::cout << setw(8) << "bin" << setw(13) << "value" << setw(13) << "Error" 
	//	<< setw(10) << "rat_err"
	//	<< std::endl; 
	TH1 *hist_err_copy = NULL;
	{
		std::string myname = hist_err->GetName() + std::string ("_copy");
		hist_err_copy = dynamic_cast<TH1*>(hist_err->Clone (myname.c_str()));
		for (unsigned bin = 1; bin <= (unsigned) hist_err_copy->GetNbinsX(); ++ bin)
		{
			float value = hist_err_copy->GetBinContent (bin);
			float error = hist_err_copy->GetBinError (bin);
			hist_err_copy->SetBinError (bin, value ? error / value : 0);
			
			//hist_err_copy->SetBinError (bin, value ? Errors.at(bin-1) / value : 0);
			hist_err_copy->SetBinContent (bin, 0);

			//std::cout << setiosflags(ios::fixed) << setprecision(3) 
			//	<< setw(8) << bin << setw(13) << value << setw(13) << Errors.at(bin-1)
			//	<< setw(10) << hist_err_copy->GetBinError (bin)
			 //	<< std::endl; 
		};
	};


	
	//for (unsigned bin = 1; bin <= (unsigned)err->GetNbinsX(); ++ bin)
	
	pho->SetTitle("");
	hist_err_copy->SetTitle("");
	std::ostringstream ytitle;


	std::string xtitle, yunits;
	if (which =="Et_pho") { xtitle += "E_{T}^{ #gamma} (GeV)"; yunits= "GeV"; }
	if (jets ==1 && which =="InvMass_pj1") { xtitle += "Invariant Mass(#gamma, Lead Jet) (GeV/c^{2})"; yunits= "GeV/c^{2}"; }
	if (jets ==2 && which =="InvMass_pj1") { xtitle += "Invariant Mass(#gamma, Lead Jet) (GeV/c^{2})"; yunits= "GeV/c^{2}"; }
	if (jets ==2 && which =="InvMass_pj1j2") { xtitle += "Invariant Mass(#gamma, Two Lead Jets) (GeV/c^{2})"; yunits= "GeV/c^{2}"; }
	if (which =="Ht") xtitle += "H_{T} (GeV)";
	if (which =="InvMass_j1j2") { xtitle += "Invariant Mass(Two Lead Jets) (GeV/c^{2})"; yunits= "GeV/c^{2}"; }
	if (which =="Et_j1") xtitle += "E_{T}^{Lead Jet} (GeV)";
	if (which =="SecondLeadJetEt") xtitle += "E_{T}^{Second Lead Jet} (GeV)";
	if (which =="NJet") xtitle += "Jet Multiplicity";
	if (which == "DelR") xtitle += "#Delta R (#gamma, lead jet) = #sqrt{#delta #phi ^{2} + #delta #eta ^{2}}";
	if (which == "DelEta") xtitle += "|#Delta #eta (#gamma, lead jet)|";
	if (which == "DelPhi") xtitle += "|#Delta #phi (#gamma, lead jet)|";
	if (which == "DetEta_pho") xtitle += "#gamma detector #eta";
	if (which == "DetPhi_pho") xtitle += "#gamma detector #phi";
	if (which == "DetEta_j1") xtitle += "Lead jet detector #eta";
	if (which == "DetPhi_j1") xtitle += "Lead jet detector #phi";
	if (which == "EvtPhi_j1") xtitle += "Lead jet 4-vec #phi";
	
	//ytitle << "#frac{Data -- Background}{Background}        Events/"<< phojet->GetBinWidth(1) << " " << yunits;
	ytitle << "Data #topbar Background / Background";
	//std::cout << "title= " << title << std::endl;
	hist_err_copy->GetYaxis()->CenterTitle(true);
	hist_err_copy->GetXaxis()->CenterTitle(true);
	hist_err_copy->SetTitleOffset(0.9,"Y");
	hist_err_copy->SetTitleOffset(0.9,"X");
	hist_err_copy->GetYaxis()->SetTitle(ytitle.str().c_str());
	hist_err_copy->GetXaxis()->SetTitle(xtitle.c_str());
	hist_err_copy->SetMinimum (-1.);
	hist_err_copy->SetMaximum (1.);
	
	pho->SetMarkerStyle(8);
	

	if (pad == 0)
	{
		std::cout <<  __FUNCTION__ << ":" << __LINE__ << ":: PAD given is NULL creating new PAD"<<  std::endl;
		new TCanvas();
	} else pad->cd();

	gStyle->SetOptStat(0);
	gStyle->SetCanvasColor (10);
	//gStyle->SetCanvasBorderSize (0);
	//gStyle->SetCanvasBorderMode (0);
				 
	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (0);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);

	int labelfont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	int titlefont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	hist_err_copy->GetXaxis()->SetLabelFont(labelfont);
	hist_err_copy->GetYaxis()->SetLabelFont(labelfont);
	//hist_err_copy->GetYaxis()->SetLabelSize(0.05);
	//hist_err_copy->GetXaxis()->SetLabelSize(0.05);
	hist_err_copy->GetXaxis()->SetTitleFont(titlefont);
	hist_err_copy->GetYaxis()->SetTitleFont(titlefont);
	//hist_err_copy->GetXaxis()->SetTitleSize(0.05);
	//hist_err_copy->GetYaxis()->SetTitleSize(0.05);
	hist_err_copy->GetXaxis()->CenterTitle(true);
	hist_err_copy->GetYaxis()->CenterTitle(true);
	hist_err_copy->GetXaxis()->SetTitleOffset(0.9);
	hist_err_copy->GetYaxis()->SetTitleOffset(0.9);
	
	pho->SetMarkerColor(kBlue);
	pho->SetLineColor(kBlue);
	// phojet->GetXaxis()->SetTitle(title.c_str());
	// phojet->GetXaxis()->CenterTitle(true);
  	hist_err_copy->SetFillStyle(3002);
  	hist_err_copy->SetFillColor(kBlack);
	
  	//TPaveText *tp = new TPaveText(0.5,0.91,0.9,0.99,"NDC");
	//tp->SetLineColor(10);
	//tp->SetTextFont(titlefont);
	//tp->AddText("CDF Run II Preliminary 2.0 fb^{-1}");


	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();


	//TLegend *leg = new TLegend (0.1,0.72,0.41,0.9);
	//leg->SetTextFont(42);
	//leg->SetTextSize(0.03);
	//leg->SetBorderSize (1);
	//leg->SetFillColor (10);

	//std::string lbl;
	//if (jets == 1) lbl = "Data (#gamma +#geq 1 Jet)";
	//if (jets == 2) lbl = "Data (#gamma +#geq 2 Jets)";

	//leg->AddEntry(pho,lbl.c_str());
	//leg->AddEntry(hist_err_copy,"Systematic Uncertainty");


	//hist_err_copy->SetLineColor(kRed);

	hist_err_copy->Draw ("E2");
	pho->Draw("SAME");
	//leg->Draw();
	//tp->Draw();

/*	TCanvas *c = dynamic_cast<TCanvas*>(gPad);
	if (c)
	{
		//c->SetFillColor(kRed);
		std::ostringstream str,str1,str2;
		//str << "ratioplot" << jets << "_" << which << ".gif";
		//c->Print (str.str().c_str());
		//str1 << "ratioplot" << jets << "_" << which << ".pdf";
		//c->Print (str1.str().c_str(),"pdf");
		str2 << "ratioplot" << jets << "_" << which << ".eps";
		c->Print (str2.str().c_str(),"eps");

	};
*/
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
double GetMax(const double def, const double up, const double down)
{
	return max(fabs(def - up), fabs(def - down));
}

//-----------------------------------------------------------------------------
//return max of 4 numbers
//-----------------------------------------------------------------------------
double GetMax(const double x1, const double x2, const double x3, const double x4)
{
	return max(max(x1,x2), max(x3,x4));
}


//-----------------------------------------------------------------------------
// derive the uncertainty for the method B, for the reweighing 
// of the sideband photon Et to the total of QCD(30%)+MC PHO(70%)
//-----------------------------------------------------------------------------
void GetSidebandReweighedError(std::vector<float>& ErrVec, const TH1* qcdhist,
						const float qcdscale,
						const int iSample,
						const std::string name,
						const std::string abspath,
						const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4, 
							bool debug=false)
{

	if (debug) std::cout << " IN " << __FUNCTION__ << ":: " << __LINE__ << std::endl;
	assert (qcdhist != NULL && "GetSidebandReweighedError:: qcdhist is NULL!");
	ErrVec.clear();

	std::cout <<  __FUNCTION__ << "::abs path = " << abspath << std::endl;

	//only when using reweighed sideband 100% from Elog#1520
	//this is the error for the reweighed sideband in addition
	//to the choice of sideband to represent fake jets in
	//tight sample.

	TFile *f = 0;
	if (iSample == iG30)
	{
		//f =  new TFile("PhoJets_data_g30_wtihsidebandCMBH_SidebandPlusSigmaRewgtByPhoEtAndJetEta.root");
		//f = new TFile ("PhoJets_data_SidebandNominalWgtPhoEtAndJetEta_g30_withsidebandCMBH.root");
		f = new TFile ("PhoJets_data_SidebandSigmaWgtPhoEtAndJetEta__g30_withsidebandCMBH_take2.root");
	} else if (iSample == iG40)
	{
		f =  new TFile("PhoJets_data_gammaEt40_RewgtIsoAddedSidebandByPlusSigmawgts.root");
	} else if (iSample == iG30MET40)
	{
		f =  new TFile("PhoJets_data_g30_met40_withsideband_sidebandreweighedPlusSigmaWeights.root");

	} else if (iSample == iG30EXCL1J)
	{
		f = new TFile("PhoJets_data_SidebandSigmaWgtPhoEt__g30Exc1jet_PhoEtReweighedSideband.root");
	} else
	{
		std::cout <<  __FUNCTION__ << ": no reweight error function for the " 
			<< WhichSample(iSample) << " sample! please check " <<  std::endl;
		assert (false);
	}
	
	assert (! f->IsZombie() && "GetSidebandReweighedError: File not found!");

	//temp hack to path
	size_t first = abspath.find_first_of("/");
	std::string part1 = abspath.substr(0,first+1);
	std::string part3 = abspath.substr(first+1,abspath.length());
	std::string newpath(part1+"SIDEBAND/"+part3+"/"+name);
	//	std::cout << "newpath = " << newpath << std::endl;

	TH1F* qcdwgthist = dynamic_cast<TH1F*> (f->Get(newpath.c_str()));
	if (! qcdwgthist){
		std::cout << __FUNCTION__ << ":" << __LINE__ << ": hist not found in the dir" <<std::endl;
		return;
	}
	
	if (bDO_VAR_BINNING)
	{
		qcdwgthist  = (TH1F*) MakeVariableBins (qcdwgthist, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false); 
	} else 
	{
		qcdwgthist->Sumw2();
		qcdwgthist->Rebin(iREBIN);
	}
	qcdwgthist->SetLineColor(kBlue);

	std::stringstream title;
	title << qcdwgthist->GetTitle() << ": ERROR :Reweighed (#epsilon) sideband comapared to Reweighed (#epsilon+#sigma) sideband";
	qcdwgthist->SetTitle(title.str().c_str());

	//do not normalize the two as they only differ by weights! 
	//normazlise the wgted to nominal
	//qcdwgthist->Scale(qcdhist->Integral()/(double)qcdwgthist->Integral());

	qcdwgthist->Scale(qcdscale);
	qcdhist->Print();
	qcdwgthist->Print();
	TCanvas *c4 = 0;
	if (debug)
	{
		c4 = new TCanvas();
		c4->Divide(1,2);
		c4->cd(1);
		gPad->SetLogy();
		qcdwgthist->DrawClone();
		qcdhist->DrawClone("same");
	}

	//qcdwgthist->Divide(qcdhist);
	if (debug)
	{
		//new TCanvas();
		c4->cd(2);
		gPad->SetTickx();
		gPad->SetTicky();
		gPad->SetGridx();
		gPad->SetGridy();

		qcdwgthist->DrawClone("P");
	}

	//now find the max for each bin
	//std::cout << __FUNCTION__ << ":" << __LINE__ << std::endl;
	//std::cout << setw(5) << "bin" << setw(13) << "val" 
	//	<< setw(13) << "FracErr" << setw(10) << "Err" << std::endl;
	for (int i=1; i <= qcdhist->GetNbinsX(); i++) 
	{
		//double err = 0;
		double dy = 0;
		if (qcdhist->GetNbinsX())
		{
			//dy = fabs(qcdwgthist->GetBinContent(i) - 1);
			//err = dy * qcdhist->GetBinContent(i);
			dy = fabs(qcdwgthist->GetBinContent(i) - qcdhist->GetBinContent(i));
		}
		
		//std::cout << setw(5) << i << setw(13) << qcdhist->GetBinContent(i) 
		//	<< setw(13) << dy << setw(10) << err << std::endl;
		//ErrVec.push_back(err);
		ErrVec.push_back(dy);
	}

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GetQCD100Err(std::vector<float>& ErrVec, const TH1* qcdhist, 
						const std::string name,
						const std::string abspath,
						const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4, 
							const int method, const std::string brname="",
							bool debug=false)
{

	if (debug) std::cout << " IN " << __FUNCTION__ << ":: " << __LINE__ << std::endl;

	//check the id vary values again.!
	assert (qcdhist != NULL && "GetQCD100Err:: qcdhist is NULL!");
	ErrVec.clear();

	std::cout <<  __FUNCTION__ << "::abs path = " << abspath << std::endl;

	/*	if (method == 1) //error based on di-jet MC studies
		{
		TFile f("SidebandSyst_FitFunctions.root");
		assert (! f.IsZombie() && "GetQCD100Err:: Root file for 'method 1' is not found!");
		assert ( brname.length() > 0 && "GetQCD100Err::brname not psecified for 'method 1'");

		TF1* func = dynamic_cast<TF1*> (f.Get(brname.c_str()));
		assert ( func != NULL && "GetQCD100Err::function not found for 'method1'");
		const int iNeedPars = 3;	//see elog#1518
		assert ( func->GetNpar() == iNeedPars && "GetQCD100Err::function not found for 'method1'");
		double par[iNeedPars];
		func->GetParameters(par);

	//new TCanvas();
	//func->Draw();
	//func->Print();

	std::stringstream hist_name;
	hist_name << "TH1F_" << brname;
	TH1* hist_temp = dynamic_cast<TH1*> (f.Get(hist_name.str().c_str()));
	assert ( hist_temp != NULL && "GetQCD100Err::hist not found for 'method1'");
	assert ( hist_temp->GetDimension() == 1 && "GetQCD100Err::hist found is not 1D!");
	//new TCanvas();
	//hist_temp->Draw();

	float first_bin_val = 0;

	for (int bin=1; bin<=hist_temp->GetNbinsX(); ++bin)
	{
	if (hist_temp->GetBinContent(bin)>0)
	{
	first_bin_val = hist_temp->GetBinContent(bin);
	break;
	}
	}
	assert (first_bin_val > 0 && "GetQCD100Err:: first_bin_val==0!");

	for (int bin=1; bin<=qcdhist->GetNbinsX(); ++bin)
	{
	//if (qcdhist->GetBinContent(bin))
	{

	float err = func->Eval(qcdhist->GetBinCenter(bin));			
	if (err < first_bin_val) err = first_bin_val;
	ErrVec.push_back(err * qcdhist->GetBinContent(bin));
	//std::cout << "bin, err% ,errval = " << bin << "\t" << err << "\t"
	//	<< qcdhist->GetBinContent(bin) * err
	//	<< std::endl;
	}
	}
	delete hist_temp;
	delete func;
	} 
	*/	

	std::string hademfile("Systematics_TightHadEm.root");
	std::string isofile("Systematics_TightIso.root");
	std::string trkptfile("Systematics_TightTrkPt.root");
	std::string trkisofile("Systematics_TightTrkIso.root");

	TFile* fhadem = new TFile(hademfile.c_str());
	TFile* fiso = new TFile(isofile.c_str());
	TFile* ftrkpt = new TFile(trkptfile.c_str());
	TFile* ftrkiso = new TFile(trkisofile.c_str());

	if (fhadem->IsZombie()) { std::cout  << hademfile << "file not found" <<std::endl; return;}
	if (fiso->IsZombie()) { std::cout  << isofile << "file not found" <<std::endl; return;}
	if (ftrkpt->IsZombie()) { std::cout  << trkptfile << "file not found" <<std::endl; return;}
	if (ftrkiso->IsZombie()) { std::cout  << trkisofile << "file not found" <<std::endl; return;}


	std::cout << abspath << std::endl;

	fhadem->cd();
	if (! gDirectory->cd(abspath.c_str())) {
		fhadem->ls();
		std::cout << __FUNCTION__<< ":" << __LINE__ << ": path '" << abspath <<"' not found?:" << abspath <<std::endl;
		assert(false);
	}

	TH1F* hademhist = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (! hademhist){
		std::cout << __FUNCTION__ << ":" << __LINE__ << ": hist not found in the dir" <<std::endl;
		assert(false);
	}

	fiso->cd();
	gDirectory->cd(abspath.c_str());
	TH1F* isohist = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	ftrkpt->cd();
	gDirectory->cd(abspath.c_str());
	TH1F* trkpthist = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	ftrkiso->cd();
	gDirectory->cd(abspath.c_str());
	TH1F* trkisohist = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	assert ( isohist != NULL && "GetQCD100Err:: isohist is null");
	assert ( ftrkpt != NULL	 && "GetQCD100Err:: ftrkpt is null");
	assert ( ftrkiso != NULL && "GetQCD100Err:: ftrkiso is null");


	if (bDO_VAR_BINNING)
	{
		hademhist  = (TH1F*) MakeVariableBins (hademhist, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false); 
		isohist    = (TH1F*) MakeVariableBins (isohist, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);  
		trkpthist  = (TH1F*) MakeVariableBins (trkpthist, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);  
		trkisohist = (TH1F*) MakeVariableBins (trkisohist, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);  
	} else
	{

		hademhist->Sumw2();
		isohist->Sumw2();
		trkpthist->Sumw2();
		trkisohist->Sumw2();
		hademhist->Rebin(iREBIN);
		isohist->Rebin(iREBIN);
		trkpthist->Rebin(iREBIN);
		trkisohist->Rebin(iREBIN);
	}

	if (hademhist->Integral("width")>0)	 hademhist->Scale(qcdhist->Integral("width")/(double) hademhist->Integral("width"));
	if (isohist->Integral("width")>0)    isohist->Scale(qcdhist->Integral("width")/(double) isohist->Integral("width"));
	if (trkpthist->Integral("width")>0)  trkpthist->Scale(qcdhist->Integral("width")/(double) trkpthist->Integral("width"));
	if (trkisohist->Integral("width")>0) trkisohist->Scale(qcdhist->Integral("width")/(double) trkisohist->Integral("width"));

	hademhist->Divide(qcdhist);
	isohist->Divide(qcdhist);
	trkpthist->Divide(qcdhist);
	trkisohist->Divide(qcdhist);

	//need to compare these varitions to the nominal sideband
	//so I need the sideband without any subtraction or addition
	//or scaling. - 01-14-2010
	//THIS is true. But when you derive the absolute error
	//from these percentages, you must use the final scaled qcd hist
	//but then what if some component like ewk is already subtracted
	//from the the final! do i derive absoluet error from that?
	//but for now this is ok, as I am not subtracting anything from the
	//sideband!-01-23-2010

	//debug = true;
	if (debug)
	{
		hademhist->SetMaximum(10);
		isohist->SetLineColor(kRed);
		isohist->SetMarkerColor(kRed);
		trkpthist->SetLineColor(kBlue);
		trkpthist->SetMarkerColor(kBlue);
		trkisohist->SetLineColor(kGreen);
		trkisohist->SetMarkerColor(kGreen);
		hademhist->SetMarkerStyle(20);
		isohist->SetMarkerStyle(21);
		trkpthist->SetMarkerStyle(22);
		trkisohist->SetMarkerStyle(23);

		std::stringstream title;
		title << qcdhist->GetTitle() << " : QCD Shape Sytematic by varying #gamma ID cuts;;Scale(QCD/Varied) / QCD" << std::endl;
		hademhist->SetTitle(title.str().c_str());

		new TCanvas();
		gPad->SetGridx();
		gPad->SetGridy();
		hademhist->DrawClone("P");
		isohist->DrawClone("sameP");
		trkpthist->DrawClone("sameP");
		trkisohist->DrawClone("sameP");

		TLegend *leg = new TLegend (0.7,0.6,0.9,0.9);
		leg->SetTextFont(42);
		leg->SetTextSize(0.03);

		leg->AddEntry(isohist,"Iso cut varied");
		leg->AddEntry(hademhist,"HadEm cut varied");
		leg->AddEntry(trkpthist,"TrkPt cut varied");
		leg->AddEntry(trkisohist,"TrkIso cut varied");
		leg->Draw();
		//gPad->SetEditable(0);
	}

	//now find the max for each bin
	//std::cout << __FUNCTION__ << ":" << __LINE__ << std::endl;
	//std::cout << setw(5) << "bin" << setw(13) << "val" 
	//	<< setw(13) << "FracErr" << setw(10) << "Err" << std::endl;
	for (int i=1; i <= qcdhist->GetNbinsX(); i++) 
	{
		float y1 = fabs(hademhist->GetBinContent(i) - 1.);
		float y2 = fabs(isohist->GetBinContent(i) - 1);
		float y3	= fabs(trkpthist->GetBinContent(i) - 1);
		float y4 = fabs(trkisohist->GetBinContent(i) - 1);

		float max = GetMax(y1,y2,y3,y4);
		float err = max * qcdhist->GetBinContent(i);
		//std::cout << setw(5) << i << setw(13) << qcdNominalhist->GetBinContent(i) 
		//	<< setw(13) << max << setw(10) << err << std::endl;
		ErrVec.push_back(err);
	}

}

 //*************** SUBTRACT HALO/COSMIC/EWK FROM QCD BACKGROUND *******
void SubtractCosmicFromQCD(TH1* qcdhist,
						const int iSample,
						const int njet,
						const std::string name,
						const std::string path,
						const float xmin, const float xpoint1,
						const float xpoint2, const float xpoint3, const float xpoint4,
						const float width1, const float width2, const float width3, 
						const float width4,
						const bool debug=false
		)
{
	if (debug)
	{
		std::cout << " IN " << __FUNCTION__ << ":: " << __LINE__ << std::endl;
		std::cout <<" got qcdhist " << std::endl;
		PrintAll(qcdhist);
	}

	
	//the old root files do not have these hists
	//only the new root files have SIDEBAND cosmic hists
	
	//When all jobs are redone, I need to change this to
	//PhoJets_data.root. -02-26-2010
	TFile *file = 0;
	if (iSample == iG30)
	{
		//const std::string sFile("PhoJets_data_gammaEt30.root");
		//std::cout <<  __FUNCTION__ << __LINE__ 
		//	<< " : WARNING ! iG30 SAMPLE USING NON STANDARD COSMIC/HALO SUBTRACTION FILE! "
		//	<< sFile <<  std::endl;
		file = new TFile("PhoJets_data_gammaEt30.root");
		//file = new TFile("PhoJets_data.root");
	} else {
		file = new TFile("PhoJets_data.root"); // the new subsamples have the sideband cosmic/halo hists
															// in the same file as all data hists in main routine
	}														//03-02-2010

	if (file->IsZombie())
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ": sideband cosmic/halo file is not found!" << std::endl;
		assert (false);
	}
	std::stringstream sideband_cosmicpath;
	sideband_cosmicpath << "Hist/SIDEBAND_COSMIC/" << path << "/" << name;
	std::cout << "sideband cosmic path = " << sideband_cosmicpath.str() << std::endl;

	TH1F* sideband_cosmic = dynamic_cast<TH1F*> (file->Get(sideband_cosmicpath.str().c_str()));
	assert ( sideband_cosmic != NULL && "sideband_cosmic hist is null!");

	if (debug)
	{
		std::cout << "Opened file "; file->Print();
		std::cout << "Retrived hist "  << sideband_cosmicpath.str() << std::endl;
		PrintAll(sideband_cosmic);
	}
	
	if (bDO_VAR_BINNING)
	{
		sideband_cosmic = (TH1F*) MakeVariableBins (sideband_cosmic, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
	} else 
	{
		sideband_cosmic->Sumw2();
		sideband_cosmic->Rebin(iREBIN);
	}

	if (debug)
	{
		std::cout << " <<<< after rebin >>> " << std::endl;
		PrintAll(sideband_cosmic);
	}


	float fSidebandCosmicEst = 0;
	if (iSample == iG30)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG30_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG40)
	{
		if (njet == 1)	fSidebandCosmicEst = fG40_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG40_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG30MET25)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30MET25_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG30MET25_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG30MET40)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30MET40_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG30MET40_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG30EXCL1J)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30EXCL1J_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	
		{
			std::cout <<  "ERROR! No 2jet case for excl njet=1 sample!";
			assert (false);
		}
	} else
	{
		assert (false && "SubtractCosmicFromQCD::Unknown Sample. No cosmic estiamtes found!");
	}


	if (fSidebandCosmicEst == 0)
		cout <<  __FUNCTION__ << ": NO SIDEBAND COSMIC ESTIMATE FOUND FOR THIS SAMPLE! CHECK!" 
			<<  endl;
	if (sideband_cosmic->Integral("width")>0)
	{
		sideband_cosmic->Scale(fSidebandCosmicEst/sideband_cosmic->Integral("width"));
	}

	std::cout <<  "sideband cosmic integral = " << sideband_cosmic->Integral("width") << endl;
	
	qcdhist->Add(sideband_cosmic, -1);

	if (debug)
	{
		std::cout << "cosmic hist info " << std::endl;
		std::cout << "sideband cosmic estiamte = " << fSidebandCosmicEst << std::endl;
		std::cout << "sideband cosmic hist after scaling " << std::endl;
		PrintAll(sideband_cosmic);
		std::cout <<"qcdhist  after subtration " << std::endl;
		PrintAll(qcdhist);
	}
	
}




 //*************** SUBTRACT EWK COMPONENT IN QCD BACKGROUND *******
 //During this subtraction if there we get a bin that it statitstically 
 //significantly negative we are doing something wrong. not just for that
 //bin, but the whole subtraction.
 //statistically signifcant mean, if the bin value c an be negative as long as
 //the statistical error includes '0' in the range.
void SubtractEWKfromQCD(TH1* qcdhist,
						const std::string name,
						const std::string path,
						const float xmin, const float xpoint1,
						const float xpoint2, const float xpoint3, const float xpoint4,
						const float width1, const float width2, const float width3, 
						const float width4, const bool debug="false" 
		)
{
	std::cout << " IN " << __FUNCTION__ << ":: " << __LINE__ << std::endl;
	// this will avoid the double counting
	//I have removed halo/cosmics from both signal and qcd templates,
	// I need to subtract the templates by normalizing them to signal/qcd.
	//but since I have an estimate of the expected number of halo/cosmics 
	// events, I can just subtract the template. halo selection is orthogonal to 
	// photon or qcd selection. cosmic fraction in signal/qcd can be asumed to
	// be same, I think.

	std::cout << __LINE__ << ":: Got qcdhist info" << std::endl;
	PrintAll(qcdhist);

	
	//pick up the ewk component in the sideband to subtract off.

	TFile *file_zee = new TFile("PhoJets_zeemc_withsideband.root");	
	TFile *file_zmm = new TFile("PhoJets_zmmmc_withsideband.root");	
	TFile *file_ztt = new TFile("PhoJets_zttmc_withsideband.root");	
	TFile *file_wen = new TFile("PhoJets_wenmc_withsideband.root");	
	TFile *file_wmn = new TFile("PhoJets_wmnmc_withsideband.root");	
	TFile *file_wtn = new TFile("PhoJets_wtnmc_withsideband.root");	

/*
	TFile *file_zee = new TFile("PhoJets_zeemc.root");	
	TFile *file_zmm = new TFile("PhoJets_zmmmc.root");	
	TFile *file_ztt = new TFile("PhoJets_zttmc.root");	
	TFile *file_wen = new TFile("PhoJets_wenmc.root");	
	TFile *file_wmn = new TFile("PhoJets_wmnmc.root");	
	TFile *file_wtn = new TFile("PhoJets_wtnmc.root");	
*/
	
	if (file_zee->IsZombie() || file_zmm->IsZombie() || file_ztt->IsZombie()
			|| file_wen->IsZombie() || file_wmn->IsZombie() || file_wtn->IsZombie())
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ": sideband ewk file is not found!" << std::endl;
		assert (false);
	}
	std::stringstream sideband_ewkpath;
	sideband_ewkpath << "Hist/SIDEBAND/" << path << "/" << name;
	std::cout << "ewk path = " << sideband_ewkpath.str() << std::endl;


	TH1F* sideband_zee = dynamic_cast<TH1F*> (file_zee->Get(sideband_ewkpath.str().c_str()));
	TH1F* sideband_zmm = dynamic_cast<TH1F*> (file_zmm->Get(sideband_ewkpath.str().c_str()));
	TH1F* sideband_ztt = dynamic_cast<TH1F*> (file_ztt->Get(sideband_ewkpath.str().c_str()));
	TH1F* sideband_wen = dynamic_cast<TH1F*> (file_wen->Get(sideband_ewkpath.str().c_str()));
	TH1F* sideband_wmn = dynamic_cast<TH1F*> (file_wmn->Get(sideband_ewkpath.str().c_str()));
	TH1F* sideband_wtn = dynamic_cast<TH1F*> (file_wtn->Get(sideband_ewkpath.str().c_str()));

	assert ( (sideband_zee != NULL && sideband_zmm != NULL && sideband_ztt != NULL
				&& sideband_wen != NULL && sideband_wmn != NULL && sideband_wtn != NULL)
			&& "one or more of the sideband_EWK hists are null!");

	sideband_zee->SetName(GetNewHistName(sideband_zee,"sideband_zee").c_str());
	sideband_zmm->SetName(GetNewHistName(sideband_zmm,"sideband_zmm").c_str());
	sideband_ztt->SetName(GetNewHistName(sideband_ztt,"sideband_ztt").c_str());
	sideband_wen->SetName(GetNewHistName(sideband_wen,"sideband_wen").c_str());
	sideband_wmn->SetName(GetNewHistName(sideband_wmn,"sideband_wmn").c_str());
	sideband_wtn->SetName(GetNewHistName(sideband_wtn,"sideband_wtn").c_str());

	if (debug)
	{
		std::cout << " <<<<<<< before rebin : line" << __LINE__ << std::endl;
		PrintFileAndHist(file_zee, sideband_zee, gDirectory);
		PrintFileAndHist(file_zmm, sideband_zmm, gDirectory);
		PrintFileAndHist(file_ztt, sideband_ztt, gDirectory);
		PrintFileAndHist(file_wen, sideband_wen, gDirectory);
		PrintFileAndHist(file_wmn, sideband_wmn, gDirectory);
		PrintFileAndHist(file_wtn, sideband_wtn, gDirectory);
	}

	if (bDO_VAR_BINNING)
	{
		sideband_zee = (TH1F*) MakeVariableBins (sideband_zee, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		sideband_zmm = (TH1F*) MakeVariableBins (sideband_zmm, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		sideband_ztt = (TH1F*) MakeVariableBins (sideband_ztt, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		sideband_wen = (TH1F*) MakeVariableBins (sideband_wen, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		sideband_wmn = (TH1F*) MakeVariableBins (sideband_wmn, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		sideband_wtn = (TH1F*) MakeVariableBins (sideband_wtn, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
	} else
	{
		sideband_zee->Sumw2();
		sideband_zmm->Sumw2();
		sideband_ztt->Sumw2();
		sideband_wen->Sumw2();
		sideband_wmn->Sumw2();
		sideband_wtn->Sumw2();

		sideband_zee->Rebin(iREBIN);
		sideband_zmm->Rebin(iREBIN);
		sideband_ztt->Rebin(iREBIN);
		sideband_wen->Rebin(iREBIN);
		sideband_wmn->Rebin(iREBIN);
		sideband_wtn->Rebin(iREBIN);
	}

	if (debug)
	{
		std::cout << " <<<<<<< after rebin : line" << __LINE__ << std::endl;
		PrintAll(sideband_zee);
		PrintAll(sideband_zmm);
		PrintAll(sideband_ztt);
		PrintAll(sideband_wen);
		PrintAll(sideband_wmn);
		PrintAll(sideband_wtn);
	}

		
	//kFAC_MC is to account for all NLO contributions to the CROSS SECTION 
	//(simply a correction to everything that we do not know
	//beyond LO). 
	//Contribution to the differential crossection by the kFAC can change 
	//sign from LO to NLO to NNLO and so on. So it is not safe to assume that
	//the kFac will decrease from LO to NLO.
	//(this sign change happens only in 'differential crosssection'. In total
	// crosssection, all higher order terms gives a POSITIVE contribution)
	const float zee_Lum = fNTOT_ZEE / (fCS_ZEE_INC * kFAC_MC);
	const float zmm_Lum = fNTOT_ZMM / (fCS_ZMM_INC * kFAC_MC);
	//I am taking the average Lum for this dataset for now 
	const float ztt_Lum = (fNTOT_ZTT_P0 / (fCS_ZTT_INC_P0 * kFAC_MC) + fNTOT_ZTT_P1_7 / (fCS_ZTT_INC_P1_7 * kFAC_MC))/2;
	const float wen_Lum = fNTOT_WEN / (fCS_WEN_INC * kFAC_MC);
	const float wmn_Lum = fNTOT_WMN / (fCS_WMN_INC * kFAC_MC);
	const float wtn_Lum = fNTOT_WTN / (fCS_WTN_INC * kFAC_MC);

	const float sideband_zee_scale = fLUM_DATA / zee_Lum;
	const float sideband_zmm_scale = fLUM_DATA / zmm_Lum;
	const float sideband_ztt_scale = fLUM_DATA / ztt_Lum;
	const float sideband_wen_scale = fLUM_DATA / wen_Lum;
	const float sideband_wmn_scale = fLUM_DATA / wmn_Lum;
	const float sideband_wtn_scale = fLUM_DATA / wtn_Lum;

	sideband_zee->Scale(sideband_zee_scale);
	sideband_zmm->Scale(sideband_zmm_scale);
	sideband_ztt->Scale(sideband_ztt_scale);
	sideband_wen->Scale(sideband_wen_scale);
	sideband_wmn->Scale(sideband_wmn_scale);
	sideband_wtn->Scale(sideband_wtn_scale);

	if (debug)
	{
		std::cout << " <<<<<<< after scaling : line" << __LINE__ << std::endl;
		std::cout << "zeenorm = " << sideband_zee_scale << std::endl;
		std::cout << "zmmnorm = " << sideband_zmm_scale << std::endl;
		std::cout << "zttnorm = " << sideband_ztt_scale << std::endl;
		std::cout << "wennorm = " << sideband_wen_scale << std::endl;
		std::cout << "wmnnorm = " << sideband_wmn_scale << std::endl;
		std::cout << "wtnnorm = " << sideband_wtn_scale << std::endl;
		PrintAll(sideband_zee);
		PrintAll(sideband_zmm);
		PrintAll(sideband_ztt);
		PrintAll(sideband_wen);
		PrintAll(sideband_wmn);
		PrintAll(sideband_wtn);
	}

	

	std::cout <<  "sideband ewk stuff integrals" << std::endl;
	cout << "zee = " << sideband_zee->Integral("width") << endl;
	cout << "zmm = " << sideband_zmm->Integral("width") << endl;
	cout << "ztt = " << sideband_ztt->Integral("width") << endl;
	cout << "wen = " << sideband_wen->Integral("width") << endl;
	cout << "wmn = " << sideband_wmn->Integral("width") << endl;
	cout << "wtn = " << sideband_wtn->Integral("width") << endl;
	

	std::cout << __FUNCTION__ << ":debug  = " << debug << std::endl;
	if (debug)
	{
		new TCanvas();
		sideband_zee->Draw();
		sideband_zmm->Draw("same");
		sideband_ztt->Draw("same");
		sideband_wen->Draw("same");
		sideband_wmn->Draw("same");
		sideband_wtn->Draw("same");
		gPad->SetEditable(0);

	}
	

	
	sideband_zee->Add(sideband_zmm);
	sideband_zee->Add(sideband_ztt);
	sideband_zee->Add(sideband_wen);
	sideband_zee->Add(sideband_wmn);
	sideband_zee->Add(sideband_wtn);
	
	qcdhist->Add(sideband_zee, -1);

/*	new TCanvas();
	gPad->SetLogy();
	qcdhist->DrawCopy();
	gPad->SetEditable(0);
*/	
	//this is a hack to avoid any bin getting negative values after this subtraction
	//I am going set those bins to zero Dec 15,2009

	// need to do this after the ID error is derived. because when I do this here
	// my QCD ID error is zero for lots of bins? need to talk with someone about this!
	//ZeroOutNegativeBins(qcdjet);
	//ZeroOutNegativeBins(qcdjet_100);


	std::cout << __LINE__ << ":: qcdhist after subtracting EWK" << std::endl;
	PrintAll(qcdhist);
	
}




void MakeHistLogAndRatio_debug (std::string which, const int jets, 
					const std::string name,const std::string title,
				 	const std::string path, const int QCDerrMtd,
				 	const float xmin, const float xpoint1, 
					const float xpoint2, const float xpoint3, 
					const float xpoint4, const float width1, 
					const float width2, const float width3, 
					const float width4, 
					const std::string brname,
					const int iSample,
					const int iSideband=0,   //0=nominal, 1=iso added, 2= reweighed
					const bool debug=false
					)
{

	if (debug) std::cout << " <<<<<<<<<<<<<<<<<<<< IN " << __FUNCTION__ << ":: " << __LINE__  << " >>>>>>>>>>>>>>>>>>>>>> " << std::endl;
	
	if ( !( iSample == iG30 || iSample == iG40 || iSample == iG30MET25 
			|| iSample == iG30MET40 || iSample == iG30EXCL1J))
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ 
			<< ": Invalid sample selection! returning." << std::endl;
		return;
	}
	if ( !(iSideband == iNOMINAL || iSideband == iISOADDED 
				|| iSideband == iISOADDEDREWGT))
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ 
			<< ": Invalid sideband selection! returning." << std::endl;
		return;
	}
		
	
	TH1::AddDirectory(kFALSE);

	//the string 'which' will be used to identify hists with same name
	// for eg: Photon Et has the name EtCorr
	// and so does the Lead jet Et, and so on..
	// if which=="Photon" (with jets==1) -> Photon EtCorr of photon+1jet
	// if which=="Lead Jet" (with jets==1) -> Lead Jet EtCorr of photon+1jet
	// if which=="Photon" (with jets==2) -> Photon EtCorr of photon+2jet
	// if which=="Lead Jet" (with jets==2) -> Lead Jet EtCorr of photon+2jet
	// if which=="Second Lead Jet" (with jets==2) -> Second Lead Jet EtCorr of photon+2jet
	// same for Inv Mass
	
	
  std::string phofile("PhoJets_data.root");
  std::string mcphofile("PhoJets_phomc.root");
  std::string zeefile("PhoJets_zeemc.root");
  std::string zmmfile("PhoJets_zmmmc.root");
  std::string zttfile("PhoJets_zttmc.root");
  std::string wenfile("PhoJets_wenmc.root");
  std::string wmnfile("PhoJets_wmnmc.root");
  std::string wtnfile("PhoJets_wtnmc.root");
  std::string diphofile("PhoJets_diphomc.root");
 
 
  TFile* fpho = new TFile(phofile.c_str());
  TFile* fmcpho = new TFile(mcphofile.c_str());
  TFile* fzee = new TFile(zeefile.c_str());
  TFile* fzmm = new TFile(zmmfile.c_str());
  TFile* fztt = new TFile(zttfile.c_str());
  TFile* fwen = new TFile(wenfile.c_str());
  TFile* fwmn = new TFile(wmnfile.c_str());
  TFile* fwtn = new TFile(wtnfile.c_str());
  TFile* fdipho = 0;
  
  if (which == "Met") 
  {
	  fdipho = new TFile(diphofile.c_str());
	  if (fdipho == NULL)
	  {
		  std::cout << __FUNCTION__ << ":" << __LINE__
			  << ": DI PHO MC FILE NOT FOUND!" << std::endl;
		  return;
	  }
  }


  if (fpho->IsZombie() ||fmcpho->IsZombie() ||
      fzee->IsZombie() ||fzmm->IsZombie() ||
      fztt->IsZombie() ||fwen->IsZombie() ||
      fwmn->IsZombie() ||fwtn->IsZombie() ) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__  
			<< " a file not found. pl check" << std::endl;
		return;
	}
  if (which == "Met" && fdipho->IsZombie())
  {
		std::cout  << "dipho mc file not found. pl check" << std::endl;
		return;
  }

  // ok path :1Jet/Photon/
  std::string sHaloDir     ("Hist/HALO/"     + path);
  std::string sCosmicDir   ("Hist/COSMIC/"   + path);
  std::string sQcdDir      ("Hist/SIDEBAND/" + path);
  std::string sSigDir      ("Hist/SIGNAL/"   + path);
  std::string sMcCentralDir("Hist/CENTRAL/"  + path);
  std::string sMcUpDir     ("Hist/EMJESUP/"  + path);
  std::string sMcDownDir   ("Hist/EMJESDOWN/"+ path);



	fpho->cd();
	//std::cout << "path="<<path<< std::endl;
	//std::cout << "dir="<<sSigDir<< std::endl;

	if (! gDirectory->cd(sSigDir.c_str())) {
	 std::cout << "path not found "<< std::endl;
	 return;
	}
	
	
	TH1F* phojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (! phojet){
	 std::cout << "hist not found in the dir" <<std::endl;
	 return;
	}
	phojet->SetName(GetNewHistName(phojet,"phojet").c_str());
	
	if (debug)
	{
		std::cout << __LINE__ << ":retreiving phojet : " << name << " hist from file/dir:" << std::endl;
		fpho->Print();
		gDirectory->pwd();
		PrintAll(phojet);
	}

	//now check if this hist has been filled (or anything is in there ti begin with
	if (phojet->GetEntries()==0)
	{
		std::cout <<  __FUNCTION__ << ":" <<__LINE__ << ":" << name
			<< " hist it empty! please check!";
		phojet->Print();
		std::cout <<  __FUNCTION__ << ":" <<__LINE__ << ":" << name
			<< " hist it empty! please check! returning!" <<  std::endl;
		return;
	}
			
	//gDirectory->pwd();
	fpho->cd();
	gDirectory->cd(sHaloDir.c_str());
	gDirectory->pwd();
	TH1F* halojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (halojet == NULL ) { std::cout << __LINE__ << ": Halo hist null!" << std::endl; return; }
	halojet->SetName(GetNewHistName(halojet,"halojet").c_str());

	if (debug)
	{
		std::cout << __LINE__ << ":retreiving halojet: " << name << " hist from file/dir:" << std::endl;
		fpho->Print();
		gDirectory->pwd();
		PrintAll(halojet);
	}
	

	fpho->cd();
	gDirectory->cd(sCosmicDir.c_str());
	gDirectory->pwd();
	TH1F* cosmicjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (cosmicjet == NULL ) { std::cout << __LINE__ << ": cosmic hist null!" << std::endl; return; }
	cosmicjet->SetName(GetNewHistName(cosmicjet,"cosmicjet").c_str());

	if (debug)
	{
		std::cout << __LINE__ << ":retreiving cosmicjet: " << name << " hist from file/dir:" << std::endl;
		fpho->Print();
		gDirectory->pwd();
		PrintAll(cosmicjet);
	}
	

	TH1F* qcdjet = 0;			      //for QCD+MC combined method
	TH1F* qcdjet_100 = 0;			//100% sideband (deafult method njet plot) and all plots in reweighed mtd)
	TFile *ftemp =0;
	if (iSideband == iNOMINAL)
	{
		fpho->cd();
		gDirectory->cd(sQcdDir.c_str());
		qcdjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
		qcdjet_100 = (TH1F*) (qcdjet->Clone("qcdjet_copy"));

		if (debug)
		{
			std::cout << __LINE__ << ":retreiving qcdjet: " << name << " hist from file/dir:" << std::endl;
			fpho->Print();
			gDirectory->pwd();
			PrintAll(qcdjet);
			std::cout << __LINE__ << ":retreiving qcdjet_100: " << name << " hist from file/dir: " << std::endl;
			PrintAll(qcdjet_100);
		}
	} else if (iSideband == iISOADDED) //iso added sideband
	{
		if (iSample == iG30) ftemp = new TFile("PhoJets_data_gammaEt30_IsoAddedSideband.root");
		else if (iSample == iG40) ftemp = new TFile("PhoJets_data_gammaEt40_IsoAddedSideband_1538_again.root");

		assert (ftemp != NULL && "file not found");
		gDirectory->cd(sQcdDir.c_str());
		qcdjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());			//for QCD+MC combined method");
		qcdjet_100 = (TH1F*) gDirectory->FindObjectAny(name.c_str());
		//qcdjet_100 = dynamic_cast<TH1F*> (qcdjet->Clone("qcdjet_copy"));
	
	} else if (iSideband == iISOADDEDREWGT)
	{
		//get nomial qcd hist first from default file
		fpho->cd();
		gDirectory->cd(sQcdDir.c_str());
		qcdjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());			//for QCD+MC combined method");
 		if (iSample == iG30)
		{
			//get the reweighed sideband from the new file
			//if (iSample == iG30) ftemp = new TFile("PhoJets_data_g30_wtihsidebandCMBH_SidebandRewgtByPhoEtAndJetEta.root");
			//if (iSample == iG30) ftemp = new TFile("PhoJets_data_SidebandSigmaWgtPhoEtAndJetEta_g30_withsidebandCMBH.root");
			
			ftemp = new TFile ("PhoJets_data_SidebandNominalWgtPhoEtAndJetEta_g30_withsidebandCMBH_take2.root");

			/*ftemp = new TFile("SidebandPhoReweightFromBG_IsoENOTaddedToQCDPho_1520_1529.root");
			std::cout <<  __LINE__ << "THIS QCD REWEIGHED FILE IS NOT CHECKED!" << std::endl;
			*/
		}
		else if (iSample == iG40) ftemp = new TFile("PhoJets_data_gammaEt40_RewgtIsoAddedSidebandByNominalwgts.root");
		else if (iSample == iG30MET40) ftemp = new TFile("PhoJets_data_g30_met40_withsideband_sidebandreweighed.root");
		else if (iSample == iG30EXCL1J) ftemp = new TFile("PhoJets_data_SidebandNominalWgtPhoEt__g30Exc1jet_PhoEtReweighedSideband.root");

		assert (ftemp != NULL && "file not found");
		//now get the reweighted side
		if (gDirectory->cd(sQcdDir.c_str()) == kTRUE)
		{
			qcdjet_100 = (TH1F*) gDirectory->FindObjectAny(name.c_str());
		} else {
			std::cout <<  " rewgt hist dir not found!" << std::endl;
			std::cout << "QCD dir = " << sQcdDir <<  std::endl;
			ftemp->ls();

			assert (false);
		}
	} else
	{
		std::cout <<  __FUNCTION__ << ":" << __LINE__ << "UNKNOWN sideband selection!" <<  std::endl;
		assert (false);
	}

	if (qcdjet == NULL ) { std::cout << __LINE__ << ": qcdjet hist null!" << std::endl; return; }
	if (qcdjet_100 == NULL ) { std::cout << __LINE__ << ": qcdjet100 hist null!" << std::endl; return; }
	qcdjet->SetName(GetNewHistName(qcdjet,"qcdjet").c_str());
	qcdjet_100->SetName(GetNewHistName(qcdjet_100,"qcdjet_100").c_str());


	//MC HISTS: jesup= jesup && emup : jesdown = jesdown && emdown

	fmcpho->cd();
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* mcphojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (mcphojet == NULL ) { std::cout << __LINE__ << ": mcpho hist null!" << std::endl; return; }
	mcphojet->SetName(GetNewHistName(mcphojet,"mcphojet").c_str());
	if (debug)
	{
		PrintFileAndHist(fmcpho,mcphojet, gDirectory);
	}

	
	fmcpho->cd();
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* mcphojetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (mcphojetJESUP == NULL ) { std::cout << __LINE__ << ": mcphoJESUP hist null!" << std::endl; return; }
	mcphojetJESUP->SetName(GetNewHistName(mcphojetJESUP,"mcphojetJESUP").c_str());
	if (debug)
	{
		PrintFileAndHist(fmcpho,mcphojetJESUP, gDirectory);
	}
	fmcpho->cd();
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* mcphojetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (mcphojetJESDOWN == NULL ) { std::cout << __LINE__ << ": mcphoJESDOWN hist null!" << std::endl; return; }
	mcphojetJESDOWN->SetName(GetNewHistName(mcphojetJESDOWN,"mcphojetJESDOWN").c_str());
	if (debug)
	{
		PrintFileAndHist(fmcpho,mcphojetJESDOWN, gDirectory);
	}

			

	fzee->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* zeejet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (zeejet == NULL ) { std::cout << __LINE__ << ": zee hist null!" << std::endl; return; }
	zeejet->SetName(GetNewHistName(zeejet,"zeejet").c_str());
	if (debug)
	{
		PrintFileAndHist(fzee,zeejet, gDirectory);
	}


	
	fzee->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* zeejetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (zeejetJESUP == NULL ) { std::cout << __LINE__ << ": zeeJESUP hist null!" << std::endl; return; }
	zeejetJESUP->SetName(GetNewHistName(zeejetJESUP,"zeejetJESUP").c_str());
	if (debug)
	{
		PrintFileAndHist(fzee,zeejetJESUP, gDirectory);
	}

	
	fzee->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* zeejetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (zeejetJESDOWN == NULL ) { std::cout << __LINE__ << ":zeeJESDOWN hist null!" << std::endl; return; }
	zeejetJESDOWN->SetName(GetNewHistName(zeejetJESDOWN,"zeejetJESDOWN").c_str());
	if (debug)
	{
		PrintFileAndHist(fzee,zeejetJESDOWN, gDirectory);
	}
			

	fzmm->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* zmmjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (zmmjet == NULL ) { std::cout << __LINE__ << ": zmm hist null!" << std::endl; return; }
	zmmjet->SetName(GetNewHistName(zmmjet,"zmmjet").c_str());
	if (debug)
	{
		PrintFileAndHist(fzmm,zmmjet, gDirectory);
	}


	
	fzmm->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* zmmjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (zmmjetJESUP == NULL ) { std::cout << __LINE__ << ": zmmJESUP hist null!" << std::endl; return; }
	zmmjetJESUP->SetName(GetNewHistName(zmmjetJESUP,"zmmjetJESUP").c_str());
	if (debug)
	{
		PrintFileAndHist(fzmm,zmmjetJESUP, gDirectory);
	}



	fzmm->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* zmmjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (zmmjetJESDOWN == NULL ) { std::cout << __LINE__ << ": zmmJESDOWN hist null!" << std::endl; return; }
	zmmjetJESDOWN->SetName(GetNewHistName(zmmjetJESDOWN,"zmmjetJESDOWN").c_str());
	if (debug)
	{
		PrintFileAndHist(fzmm,zmmjetJESDOWN, gDirectory);
	}
			

	
	fztt->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* zttjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (zttjet == NULL ) { std::cout << __LINE__ << ": ztt hist null!" << std::endl; return; }
	zttjet->SetName(GetNewHistName(zttjet,"zttjet").c_str());
	if (debug)
	{
		PrintFileAndHist(fztt,zttjet, gDirectory);
	}


	
	fztt->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* zttjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (zttjetJESUP == NULL ) { std::cout << __LINE__ << ": zttJESUP hist null!" << std::endl; return; }
	zttjetJESUP->SetName(GetNewHistName(zttjetJESUP,"zttjetJESUP").c_str());
	if (debug)
	{
		PrintFileAndHist(fztt,zttjetJESUP, gDirectory);
	}


	
	fztt->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* zttjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (zttjetJESDOWN == NULL ) { std::cout << __LINE__ << ": zttJESDOWN hist null!" << std::endl; return; }
	zttjetJESDOWN->SetName(GetNewHistName(zttjetJESDOWN,"zttjetJESDOWN").c_str());
	if (debug)
	{
		PrintFileAndHist(fztt,zttjetJESDOWN, gDirectory);
	}
			

	
	fwen->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* wenjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (wenjet == NULL ) { std::cout << __LINE__ << ": wen hist null!" << std::endl; return; }
	wenjet->SetName(GetNewHistName(wenjet,"wenjet").c_str());
	if (debug)
	{
		PrintFileAndHist(fwen,wenjet, gDirectory);
	}


	
	fwen->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* wenjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (wenjetJESUP == NULL ) { std::cout << __LINE__ << ": wenJESUP hist null!" << std::endl; return; }
	wenjetJESUP->SetName(GetNewHistName(wenjetJESUP,"wenjetJESUP").c_str());
	if (debug)
	{
		PrintFileAndHist(fwen,wenjetJESUP, gDirectory);
	}


	
	fwen->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* wenjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (wenjetJESDOWN == NULL ) { std::cout << __LINE__ << ": wenJESDOWN hist null!" << std::endl; return; }
	wenjetJESDOWN->SetName(GetNewHistName(wenjetJESDOWN,"wenjetJESDOWN").c_str());
	if (debug)
	{
		PrintFileAndHist(fwen,wenjetJESDOWN, gDirectory);
	}
			

	
	fwmn->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* wmnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (wmnjet == NULL ) { std::cout << __LINE__ << ": wmn hist null!" << std::endl; return; }
	wmnjet->SetName(GetNewHistName(wmnjet,"wmnjet").c_str());
	if (debug)
	{
		PrintFileAndHist(fwmn,wmnjet, gDirectory);
	}


	
	fwmn->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* wmnjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (wmnjetJESUP == NULL ) { std::cout << __LINE__ << ": wmnJESUP hist null!" << std::endl; return; }
	wmnjetJESUP->SetName(GetNewHistName(wmnjetJESUP,"wmnjetJESUP").c_str());
	if (debug)
	{
		PrintFileAndHist(fwmn,wmnjetJESUP, gDirectory);
	}


	
	fwmn->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* wmnjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (wmnjetJESDOWN == NULL ) { std::cout << __LINE__ << ": wmnJESDOWN hist null!" << std::endl; return; }
	wmnjetJESDOWN->SetName(GetNewHistName(wmnjetJESDOWN,"wmnjetJESDOWN").c_str());
	if (debug)
	{
		PrintFileAndHist(fwmn,wmnjetJESDOWN, gDirectory);
	}
			

	
	fwtn->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* wtnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (wtnjet == NULL ) { std::cout << __LINE__ << ": wtn hist null!" << std::endl; return; }
	wtnjet->SetName(GetNewHistName(wtnjet,"wtnjet").c_str());
	if (debug)
	{
		PrintFileAndHist(fwtn,wtnjet, gDirectory);
	}


	
	fwtn->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* wtnjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (wtnjetJESUP == NULL ) { std::cout << __LINE__ << ": wtnJESUP hist null!" << std::endl; return; }
	wtnjetJESUP->SetName(GetNewHistName(wtnjetJESUP,"wtnjetJESUP").c_str());
	if (debug)
	{
		PrintFileAndHist(fwtn,wtnjetJESUP, gDirectory);
	}


	
	fwtn->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* wtnjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (wtnjetJESDOWN == NULL ) { std::cout << __LINE__ << ": wtnJESDWON hist null!" << std::endl; return; }
	wtnjetJESDOWN->SetName(GetNewHistName(wtnjetJESDOWN,"wtnjetJESDOWN").c_str());
	if (debug)
	{
		PrintFileAndHist(fwtn,wtnjetJESDOWN, gDirectory);
	}

	TH1F *diphojet = 0, *diphojetJESUP = 0, *diphojetJESDOWN = 0;
	if (which == "Met")
	{
		fdipho->cd();
		gDirectory->cd(sMcCentralDir.c_str());
		diphojet = dynamic_cast<TH1F*> (gDirectory->FindObjectAny(name.c_str()));
		assert (diphojet != NULL && "diphojet central hist not found !");

		gDirectory->cd(sMcUpDir.c_str());
		diphojetJESUP = dynamic_cast<TH1F*> (gDirectory->FindObjectAny(name.c_str()));
		assert (diphojetJESUP != NULL && "diphojet jes up hist not found !");

		gDirectory->cd(sMcDownDir.c_str());
		diphojetJESDOWN = dynamic_cast<TH1F*> (gDirectory->FindObjectAny(name.c_str()));
		assert (diphojetJESDOWN != NULL && "diphojet jes down hist not found !");
	}

	cout <<  __LINE__ << ":before rebin " << endl;
	phojet->Print();
	qcdjet->Print();
	halojet->Print();
	cosmicjet->Print();
	zeejet->Print();
	zmmjet->Print();
	zttjet->Print();
	cout <<  endl;

	
	if (bDO_VAR_BINNING)
	{
		phojet = (TH1F*) MakeVariableBins (phojet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		halojet = (TH1F*) MakeVariableBins (halojet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		zeejet = (TH1F*) MakeVariableBins (zeejet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		zmmjet = (TH1F*) MakeVariableBins (zmmjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		zttjet = (TH1F*) MakeVariableBins (zttjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		wenjet = (TH1F*) MakeVariableBins (wenjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		wmnjet = (TH1F*) MakeVariableBins (wmnjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		wtnjet = (TH1F*) MakeVariableBins (wtnjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		cosmicjet = (TH1F*) MakeVariableBins (cosmicjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		qcdjet = (TH1F*) MakeVariableBins (qcdjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		mcphojet = (TH1F*) MakeVariableBins (mcphojet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		qcdjet_100 = (TH1F*) MakeVariableBins (qcdjet_100, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);

		zeejetJESUP = (TH1F*) MakeVariableBins (zeejetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		zmmjetJESUP = (TH1F*) MakeVariableBins (zmmjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		zttjetJESUP = (TH1F*) MakeVariableBins (zttjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		wenjetJESUP = (TH1F*) MakeVariableBins (wenjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		wmnjetJESUP = (TH1F*) MakeVariableBins (wmnjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		wtnjetJESUP = (TH1F*) MakeVariableBins (wtnjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		mcphojetJESUP = (TH1F*) MakeVariableBins (mcphojetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);

		zeejetJESDOWN = (TH1F*) MakeVariableBins (zeejetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		zmmjetJESDOWN = (TH1F*) MakeVariableBins (zmmjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		zttjetJESDOWN = (TH1F*) MakeVariableBins (zttjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		wenjetJESDOWN = (TH1F*) MakeVariableBins (wenjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		wmnjetJESDOWN = (TH1F*) MakeVariableBins (wmnjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		wtnjetJESDOWN = (TH1F*) MakeVariableBins (wtnjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
		mcphojetJESDOWN = (TH1F*) MakeVariableBins (mcphojetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,false);

		if (which == "Met")
		{
			diphojet = (TH1F*) MakeVariableBins (diphojet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,false);
			diphojetJESUP = (TH1F*) MakeVariableBins (diphojetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,false);
			diphojetJESDOWN = (TH1F*) MakeVariableBins (diphojetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,false);
		}

	} else
	{
		phojet->Print("all");
		halojet->Print("all");
		zeejet->Print("all");
		zmmjet->Print("all");
		zttjet->Print("all");
		wenjet->Print("all");
		wmnjet->Print("all");
		wtnjet->Print("all");
		cosmicjet->Print("all");
		qcdjet->Print("all");
		mcphojet->Print("all");
		qcdjet_100->Print("all");

		zeejetJESUP->Print("all");
		zmmjetJESUP->Print("all");
		zttjetJESUP->Print("all");
		wenjetJESUP->Print("all");
		wmnjetJESUP->Print("all");
		wtnjetJESUP->Print("all");
		mcphojetJESUP->Print("all");

		zeejetJESDOWN->Print("all");
		zmmjetJESDOWN->Print("all");
		zttjetJESDOWN->Print("all");
		wenjetJESDOWN->Print("all");
		wmnjetJESDOWN->Print("all");
		wtnjetJESDOWN->Print("all");
		mcphojetJESDOWN->Print("all");

		
		phojet->Rebin(iREBIN);
		halojet->Rebin(iREBIN);
		zeejet->Rebin(iREBIN);
		zmmjet->Rebin(iREBIN);
		zttjet->Rebin(iREBIN);
		wenjet->Rebin(iREBIN);
		wmnjet->Rebin(iREBIN);
		wtnjet->Rebin(iREBIN);
		cosmicjet->Rebin(iREBIN);
		qcdjet->Rebin(iREBIN);
		mcphojet->Rebin(iREBIN);
		qcdjet_100->Rebin(iREBIN);

		zeejetJESUP->Rebin(iREBIN);
		zmmjetJESUP->Rebin(iREBIN);
		zttjetJESUP->Rebin(iREBIN);
		wenjetJESUP->Rebin(iREBIN);
		wmnjetJESUP->Rebin(iREBIN);
		wtnjetJESUP->Rebin(iREBIN);
		mcphojetJESUP->Rebin(iREBIN);

		zeejetJESDOWN->Rebin(iREBIN);
		zmmjetJESDOWN->Rebin(iREBIN);
		zttjetJESDOWN->Rebin(iREBIN);
		wenjetJESDOWN->Rebin(iREBIN);
		wmnjetJESDOWN->Rebin(iREBIN);
		wtnjetJESDOWN->Rebin(iREBIN);
		mcphojetJESDOWN->Rebin(iREBIN);

		if (which == "Met")
		{
			diphojet->Print("all");
			diphojetJESUP->Print("all");
			diphojetJESDOWN->Print("all");

			diphojet->Rebin(iREBIN);
			diphojetJESUP->Rebin(iREBIN);
			diphojetJESDOWN->Rebin(iREBIN);

		}

	}

	if (debug)
	{
		cout <<  "<<<<<<<<< " << __LINE__ << ": after rebin " << endl;

		PrintAll(phojet);
		PrintAll(halojet);
		PrintAll(cosmicjet);
		PrintAll(qcdjet);
		PrintAll(qcdjet_100);
		PrintAll(mcphojet);
		PrintAll(mcphojetJESUP);
		PrintAll(mcphojetJESDOWN);
		PrintAll(zeejet);
		PrintAll(zeejetJESUP);
		PrintAll(zeejetJESDOWN);
		PrintAll(zmmjet);
		PrintAll(zmmjetJESUP);
		PrintAll(zmmjetJESDOWN);
		PrintAll(zttjet);
		PrintAll(zttjetJESUP);
		PrintAll(zttjetJESDOWN);
		PrintAll(wenjet);
		PrintAll(wenjetJESUP);
		PrintAll(wenjetJESDOWN);
		PrintAll(wmnjet);
		PrintAll(wmnjetJESUP);
		PrintAll(wmnjetJESDOWN);
		PrintAll(wtnjet);
		PrintAll(wtnjetJESUP);
		PrintAll(wtnjetJESDOWN);
	}
	
/*	phojet->Print();
	qcdjet->Print();
	halojet->Print();
	cosmicjet->Print();
	zeejet->Print();
	zmmjet->Print();
	zttjet->Print();
	cout <<  endl;
*/

	
/************** FAKE PHOTON FRACTION ************************************/
	//this is used when the combination of QCD+PHO MC is used
  // this is the amount of fake photons in the signal we select(jets faking photon) which is ~30%
 	
  float fake_pho_frac_d = 0;

	if (iSample == iG30) 		fake_pho_frac_d = fG30_FAKE_PHO_FRAC;  	//for pho Et>30GeV+jets
	else if (iSample == iG40) 	fake_pho_frac_d = fG40_FAKE_PHO_FRAC; 	//for pho Et>40GeV+jets
	else if (iSample == iG30MET25)
	{
		if (jets == 1) 		fake_pho_frac_d = fG30MET25_FAKE_PHO_FRAC_1JET;
		else if (jets == 2) 	fake_pho_frac_d = fG30MET25_FAKE_PHO_FRAC_2JET;
	} else if (iSample == iG30MET40)
	{
		if (jets == 1) 		fake_pho_frac_d = fG30MET40_FAKE_PHO_FRAC_1JET;
		else if (jets == 2) 	fake_pho_frac_d = fG30MET40_FAKE_PHO_FRAC_2JET;
	} else if (iSample == iG30EXCL1J) 
	{
		if (jets == 1) fake_pho_frac_d = fG30EXCL1J_FAKE_PHO_FRAC;
		else if (jets == 2)
		{
			std::cout <<  __FUNCTION__ << __LINE__ << ": ERROR! no fake photon fraction for excl 1jet sample!" <<  std::endl;
			assert (false);
		}
	} else
	{
		std::cout <<  __FUNCTION__ << __LINE__ << ": UNKNOWN sample selction! invalid fake photon fraction!" <<  std::endl;
		assert (false);
	}

	//systematic error for all sample are taken to this assuming the change
	//of this value is small - 02-17-2010
  const float fake_pho_frac_sigma = fFAKE_PHO_FRAC_SIGMA;
  
  float fake_pho_frac_p = fake_pho_frac_d + fake_pho_frac_sigma;
  //float fake_pho_frac_m = fake_pho_frac_d - fake_pho_frac_sigma;
  
  float qcd_scale = fake_pho_frac_d;					//these background need to be scaled to these values
  float phomc_scale  = 1.0 - fake_pho_frac_d;

  std::string qcd_str, mc_str;
  std::ostringstream qcdnum, mcnum;
  qcdnum << "qcd (#gamma sideband, " << qcd_scale<< " of data)";
  mcnum << "#gamma MC (" << phomc_scale << " of data)";
  qcd_str = qcdnum.str();		// to be used in the legend of the final plot
  mc_str = mcnum.str();

   //this part is only to generate reweight functions
	//MakePhoEtWeightsForSideband(qcdjet,mcphojet,fake_pho_frac_d, fake_pho_frac_sigma, brname, jets);
	//return;


  
	/************************************************************************/
  //******************* NORMALIZING *************************************
	float haloEst 		= 0;
	float cosmicEst 	= 0;
	
	if (jets ==1 ) {
		if (iSample == iG30)
		{
			haloEst = fG30_HALO_EST_1JET;				//these estimates are for the full dataset
			cosmicEst = fG30_COSMIC_EST_1JET;
		} else if (iSample == iG40)
		{
			haloEst = fG40_HALO_EST_1JET;
			cosmicEst = fG40_COSMIC_EST_1JET;
		} else if (iSample == iG30MET25)
		{
			haloEst = fG30MET25_HALO_EST_1JET;
			cosmicEst = fG30MET25_COSMIC_EST_1JET;
		} else if (iSample == iG30MET40)
		{
			haloEst = fG30MET40_HALO_EST_1JET;
			cosmicEst = fG30MET40_COSMIC_EST_1JET;
		} else if (iSample == iG30EXCL1J)
		{
			haloEst = fG30EXCL1J_HALO_EST_1JET;				//these estimates are for the full dataset
			cosmicEst = fG30EXCL1J_COSMIC_EST_1JET;
		}
		
	} else if (jets ==2 ) {
		if (iSample == iG30)
		{
			haloEst = fG30_HALO_EST_2JET;
			cosmicEst = fG30_COSMIC_EST_2JET ;
		} else if (iSample == iG40)
		{
			haloEst = fG40_HALO_EST_2JET;
			cosmicEst = fG40_COSMIC_EST_2JET;
		} else if (iSample == iG30MET25)
		{
			haloEst = fG30MET25_HALO_EST_2JET;
			cosmicEst = fG30MET25_COSMIC_EST_2JET;
		} else if (iSample == iG30MET40)
		{
			haloEst = fG30MET40_HALO_EST_2JET;
			cosmicEst = fG30MET40_COSMIC_EST_2JET;
		} else if (iSample == iG30EXCL1J)
		{
			std::cout <<  " ERROR! no cosmic/halo estimate for excl 1jet case!" <<  std::endl;
			assert (false);
		}

	} else {
		std::cout << __FILE__ <<"::"<<__LINE__<<":: Invalid number of jets !" << std::endl;
		exit (1);
	}


	float qcdNorm    = 0; 
	float mcphoNorm  = 0;
	float qcd100Norm = 0;		//for the plots that we use 100% of QCD sideband

	//check this for constant rebin option --????????????????????????
	if (halojet->Integral("width")>0) halojet->Scale (haloEst /(double) halojet->Integral("width"));		//need to use "width" as I am rebinning these hists at the begining
	if (cosmicjet->Integral("width")>0) cosmicjet->Scale (cosmicEst/(double) cosmicjet->Integral("width"));

	
	if (debug)
	{
		std::cout << " halo estimate = " << haloEst << std::endl;
		std::cout << " cosmic estimate = " << cosmicEst << std::endl;
		std::cout << "halojet hist after scaling by integral to estimated" << std::endl;
		PrintAll(halojet);
		std::cout << "cosmicjet hist after scaling by integral to estimated" << std::endl;
		PrintAll(cosmicjet);
	}

	
	const float dataLum=2043.0;	//pb-1

	const float kFac = 1.4;
	// SF = DATA_LUM * ( 1 / EWK_LUM ) * KFAC
	//    = DATA_LUM * ( 1 / (TOT EVTS PROCESSED/CROSS SECTION) ) * KFAC
	//    here are the values picked up from pp.102
	//    DATASET   TOTAL EVTS   CROSSSECTION   LUM == TOTAL EVTS/CROSSSECTION
	//    Zee       12,092,155    355pb          34062 pb-1
	//    Zmm       13,755,133    355pb          38746 pb-1
	//    Ztt       
/*	const float zee_inc_cross_section = 355; //pb-1
	const float zmm_inc_cross_section = 355; //pb-1
	const float ztt_inc_cross_section_p0 = 355; //pb-1
	const float ztt_inc_cross_section_p1_7 = 238; //pb-1
	const float wen_inc_cross_section = 1960; //pb-1
	const float wmn_inc_cross_section = 1960; //pb-1
	const float wtn_inc_cross_section = 1960; //pb-1
*/	
	const float zeenorm = (dataLum/34056)*kFac;    // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const float zmmnorm = (dataLum/38732)*kFac;    // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const float zttnorm = (dataLum/27755)*kFac;    // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const float wennorm = (dataLum/9438)*kFac;     // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const float wmnnorm = (dataLum/5183)*kFac;     // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const float wtnnorm = (dataLum/3520)*kFac;     // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	
	if (which == "Met")
	{
		const float diphomcnorm = (dataLum/fLUM_DIPHOMC) * kFac;	
		diphojet->Scale(diphomcnorm);
		diphojetJESUP->Scale (diphomcnorm);
		diphojetJESDOWN->Scale (diphomcnorm);
		std::cout << __LINE__ << ":di pho int after scale = " << diphomcnorm 
			<< ", " << diphojet->Integral("width") << std::endl;
		/*new TCanvas();
		diphojet->Draw();
		diphojetJESUP->Draw("same");
		diphojetJESDOWN->Draw("same");
		gPad->SetEditable(0);
		*/
	}


	//makes copies for luminsoity systematic error calculations
  	TH1* zeelumsyst = (TH1*) zeejet->Clone("zeelumcopy");
  	TH1* zmmlumsyst = (TH1*) zmmjet->Clone("zmmlumcopy");
  	TH1* zttlumsyst = (TH1*) zttjet->Clone("zttlumcopy");
  	TH1* wenlumsyst = (TH1*) wenjet->Clone("wenlumcopy");
  	TH1* wmnlumsyst = (TH1*) wmnjet->Clone("wmnlumcopy");
  	TH1* wtnlumsyst = (TH1*) wtnjet->Clone("wtnlumcopy");

	
	zeejet->Scale (zeenorm);
	zmmjet->Scale (zmmnorm);
	zttjet->Scale (zttnorm);
	wenjet->Scale (wennorm);
	wmnjet->Scale (wmnnorm);
	wtnjet->Scale (wtnnorm);
	

	if (debug)
	{
		std::cout << " <<<<<< EWK NORMALIZATION >>>>>>> " << std::endl;
		std::cout << " data Lum = " << dataLum << std::endl;
		std::cout << " zeejet scale = " << zeenorm << std::endl;
		std::cout << " zmmjet scale = " << zmmnorm << std::endl;
		std::cout << " zttjet scale = " << zttnorm << std::endl;
		std::cout << " wenjet scale = " << wennorm << std::endl;
		std::cout << " wmnjet scale = " << wmnnorm << std::endl;
		std::cout << " wtnjet scale = " << wtnnorm << std::endl;
		std::cout << "zeejet hist after scaling " << std::endl;
		PrintAll(zeejet);
		PrintAll(zmmjet);
		PrintAll(zttjet);
		PrintAll(wenjet);
		PrintAll(wmnjet);
		PrintAll(wtnjet);
	}

	

	
	//need to do these steps for the QCDID samples 04-21-2010
	SubtractEWKfromQCD(qcdjet, name, path, xmin, xpoint1, xpoint2, 
					xpoint3, xpoint4,	width1, width2, width3,	width4,1); 

	//SubtractEWKfromQCD(qcdjet_100, name, path, xmin, xpoint1, xpoint2, 
	//				xpoint3, xpoint4,	width1, width2, width3,	width4,0); 

	SubtractCosmicFromQCD(qcdjet,	iSample,	jets,	name,	path,
						xmin, xpoint1,	xpoint2, xpoint3, xpoint4,
						width1, width2, width3,	width4,1); 
	//SubtractCosmicFromQCD(qcdjet_100,	iSample,	jets,	name,	path,
	//					xmin, xpoint1,	xpoint2, xpoint3, xpoint4,
	//					width1, width2, width3,	width4); 

	CheckNegativeBin(qcdjet);


	
	//now calculate the difference (or remainder) between data and halo+cosmic+ewk
	//backgrounds.
	//double leftover = phojet->Integral() - halojet->Integral() - cosmicjet->Integral() 
	//						- zeejet->Integral() - zmmjet->Integral() - zttjet->Integral()
	//						- wenjet->Integral() - wmnjet->Integral() - wtnjet->Integral();
	
	cout <<  __LINE__ << ": leftover calculation " << endl;

	double leftover = phojet->Integral("width") - halojet->Integral("width") - cosmicjet->Integral("width") 
							- zeejet->Integral("width") - zmmjet->Integral("width") - zttjet->Integral("width")
							- wenjet->Integral("width") - wmnjet->Integral("width") - wtnjet->Integral("width");


	if (which =="Met") leftover = leftover - diphojet->Integral("width");
	
	std::cout << __LINE__ <<  "leftover  = " << leftover << std::endl;
	std::cout << __LINE__ <<  "qcd_scale   = " << qcd_scale << std::endl;
	std::cout << __LINE__ <<  "phomc_scale = " << phomc_scale << std::endl;

	qcdNorm     = (leftover * qcd_scale) / qcdjet->Integral("width");
	mcphoNorm   = (leftover * phomc_scale) / mcphojet->Integral("width");

	qcd100Norm  = leftover / qcdjet_100->Integral("width"); 


	std::cout << __LINE__ << "BFORE SCALING integral of qcdjet     = " << qcdjet->Integral("width") << std::endl;
	std::cout << __LINE__ << "BFORE SCALING integral of qcdjet_100 = " << qcdjet_100->Integral("width") << std::endl;
	std::cout << __LINE__ << "BFORE SCALING integral of mcphojet   = " << mcphojet->Integral("width") << std::endl;
	std::cout << __LINE__ << "BFORE SCALING integral of mcphojetJESUP   = " << mcphojetJESUP->Integral("width") << std::endl;
	std::cout << __LINE__ << "BFORE SCALING integral of mcphojetJESDOWN = " << mcphojetJESDOWN->Integral("width") << std::endl;
	
	qcdjet_100->Scale (leftover / (double) qcdjet_100->Integral("width"));
	qcdjet->Scale ((leftover * qcd_scale) / (double) qcdjet->Integral("width"));
	
	//qcdjet_100->Scale (qcd100Norm);
	//qcdjet->Scale (qcdNorm);
	mcphojet->Scale (mcphoNorm);

	//JES PHO MC
	float mcphoJESUPNorm = (leftover * phomc_scale)/mcphojetJESUP->Integral("width");
	float mcphoJESDOWNNorm = (leftover * phomc_scale)/mcphojetJESDOWN->Integral("width");
	

	std::cout << __LINE__ << "AFTER SCALING integral of qcdjet     = " << qcdjet->Integral("width") << std::endl;
	std::cout << __LINE__ << "AFTER SCALING integral of qcdjet_100 = " << qcdjet_100->Integral("width") << std::endl;
	std::cout << __LINE__ << "AFTER SCALING integral of mcphojet   = " << mcphojet->Integral("width") << std::endl;
	std::cout << __LINE__ << "AFTER SCALING integral of mcphojetJESUP   = " << mcphojetJESUP->Integral("width") << std::endl;
	std::cout << __LINE__ << "AFTER SCALING integral of mcphojetJESDOWN = " << mcphojetJESDOWN->Integral("width") << std::endl;
	
	std::cout << " integral data / sum = " << phojet->Integral() << " / ";
	if (QCDerrMtd == 1)
	{
		std::cout << halojet->Integral() + cosmicjet->Integral()
					+ zeejet->Integral() + zmmjet->Integral() + zttjet->Integral()
					+ wenjet->Integral() + wmnjet->Integral() + wtnjet->Integral() 
					+ qcdjet_100->Integral()
					<< std::endl;
	} else if (QCDerrMtd == 2)
	{

		std::cout << halojet->Integral() + cosmicjet->Integral()
					+ zeejet->Integral() + zmmjet->Integral() + zttjet->Integral()
					+ wenjet->Integral() + wmnjet->Integral() + wtnjet->Integral() 
					+ qcdjet->Integral() + mcphojet->Integral()
					<< std::endl;
	}
	//for debug only
	//MakeTruePhoFracStudy(phojet,mcphojet,qcdjet);
	//return;
	//end of debug stuff
	

	// JES 
	mcphojetJESUP->Scale (mcphoJESUPNorm);
	zeejetJESUP->Scale (zeenorm);
	zmmjetJESUP->Scale (zmmnorm);
	zttjetJESUP->Scale (zttnorm);
	wenjetJESUP->Scale (wennorm);
	wmnjetJESUP->Scale (wmnnorm);
	wtnjetJESUP->Scale (wtnnorm);

	mcphojetJESDOWN->Scale (mcphoJESDOWNNorm);
	zeejetJESDOWN->Scale (zeenorm);
	zmmjetJESDOWN->Scale (zmmnorm);
	zttjetJESDOWN->Scale (zttnorm);
	wenjetJESDOWN->Scale (wennorm);
	wmnjetJESDOWN->Scale (wmnnorm);
	wtnjetJESDOWN->Scale (wtnnorm);


	if (debug)
	{
		
		cout <<  "<<<<<<<<< " << __LINE__ << ": after scaling all hists" << endl;

		PrintAll(phojet);
		PrintAll(halojet);
		PrintAll(cosmicjet);
		PrintAll(qcdjet);
		PrintAll(qcdjet_100);
		PrintAll(mcphojet);
		PrintAll(mcphojetJESUP);
		PrintAll(mcphojetJESDOWN);
		PrintAll(zeejet);
		PrintAll(zeejetJESUP);
		PrintAll(zeejetJESDOWN);
		PrintAll(zmmjet);
		PrintAll(zmmjetJESUP);
		PrintAll(zmmjetJESDOWN);
		PrintAll(zttjet);
		PrintAll(zttjetJESUP);
		PrintAll(zttjetJESDOWN);
		PrintAll(wenjet);
		PrintAll(wenjetJESUP);
		PrintAll(wenjetJESDOWN);
		PrintAll(wmnjet);
		PrintAll(wmnjetJESUP);
		PrintAll(wmnjetJESDOWN);
		PrintAll(wtnjet);
		PrintAll(wtnjetJESUP);
		PrintAll(wtnjetJESDOWN);
	}


	
  //***************END  NORMALIZING *************************************


	std::cout << "\n" << std::endl; 
	std::cout <<  __LINE__ << ":: ALL  HISTS " << std::endl;
	std::cout << __LINE__ << ":: phojet Integral     =  " << phojet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: halojet Integral    =  " << halojet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: cosmicjet Integral  =  " << cosmicjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: qcdjet Integral     =  " << qcdjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: qcdjet_100 Integral =  " << qcdjet_100->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: mcphojet Integral   =  " << mcphojet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: zeejet Integral     =  " << zeejet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: zmmjet Integral     =  " << zmmjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: zttjet Integral     =  " << zttjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: wenjet Integral     =  " << wenjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: wmnjet Integral     =  " << wmnjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: wtnjet Integral     =  " << wtnjet->Integral("width") << std::endl;
	if (which == "Met")
	{
		std::cout << __LINE__ << ":: diphojet Integral   =  " << diphojet->Integral("width") <<  std::endl;
	}


  
	//****************** SYSTEMATICS *************************************

  std::vector<float> cosmicErr;   //relative error for each bin is stored here
  std::vector<float> haloErr;

  GetCosmicErr(cosmicErr, cosmicjet, false);
  GetHaloErr(haloErr, halojet);

	// Luminosity Error +/-6%
	// zeeLumErr = (Ldata/ (Lmc +Lmc*0.06)) *kFac = zeeLumErr/1.06)
  //const float fact = fabs(1.0-100.0/106.0); //?????????????????????? why did I include 1???? 4-15-2010
  //
  //ok. here is what I tried to accomplish here.
  //I can change Ldata or Lmc by +/-6%. I choose to vary Lmc. 
  //i.e. Lmc' = Lmc * 1.06
  //Nominal scale = Ldata/Lmc
  //Varied scale = Ldata/Lmc' = Ldata/(Lmc *1.06) = Ldata/Lmc * (1/1.06) = zeenorm * (1/1.06)
  //so fact = 1/1.06 = 100/106;  
  //so the 1-100/106 is wrong. it should be 100/106.
  
  //const float fact = 1.06;		//Jay's fact					//this is where Jay and I stop on 04-15-2010
  const float fact = 1.0/1.06;
  const float zeeLumErr = zeenorm * fact;
  const float zmmLumErr = zmmnorm * fact;
  const float zttLumErr = zttnorm * fact;
  const float wenLumErr = wennorm * fact;
  const float wmnLumErr = wmnnorm * fact;
  const float wtnLumErr = wtnnorm * fact;

 	zeelumsyst->Scale(zeeLumErr);
 	zmmlumsyst->Scale(zmmLumErr);
 	zttlumsyst->Scale(zttLumErr);
 	wenlumsyst->Scale(wenLumErr);
 	wmnlumsyst->Scale(wmnLumErr);
 	wtnlumsyst->Scale(wtnLumErr);
	
	

	std::vector<float> vZeeLumErr,vZmmLumErr, vZttLumErr, vWenLumErr, vWmnLumErr, vWtnLumErr;

	std::cout << __LINE__ << ":: EWK LUM error calculation " << std::endl;
	std::cout << __LINE__ << ":: Scaled EWK LUM error hists " << std::endl;
	PrintAll(zeelumsyst);
	PrintAll(zmmlumsyst);
	PrintAll(zttlumsyst);
	PrintAll(wenlumsyst);
	PrintAll(wmnlumsyst);
	PrintAll(wtnlumsyst);


	std::cout << "bin" 
		<< std::setw(10) << "zee"
		<< std::setw(10) << "zmm" 
		<< std::setw(10) << "ztt" 
		<< std::setw(10) << "wen" 
		<< std::setw(10) << "wmn" 
		<< std::setw(10) << "wtn" << std::endl; 

	for (int bin=1;bin<= zeelumsyst->GetNbinsX(); ++bin)
	{
		vZeeLumErr.push_back(fabs(zeejet->GetBinContent(bin) - zeelumsyst->GetBinContent(bin)));
		vZmmLumErr.push_back(fabs(zmmjet->GetBinContent(bin) - zmmlumsyst->GetBinContent(bin)));
		vZttLumErr.push_back(fabs(zttjet->GetBinContent(bin) - zttlumsyst->GetBinContent(bin)));
		vWenLumErr.push_back(fabs(wenjet->GetBinContent(bin) - wenlumsyst->GetBinContent(bin)));
		vWmnLumErr.push_back(fabs(wmnjet->GetBinContent(bin) - wmnlumsyst->GetBinContent(bin)));
		vWtnLumErr.push_back(fabs(wtnjet->GetBinContent(bin) - wtnlumsyst->GetBinContent(bin)));
		
		std::cout << bin 
		<< std::setw(10) << vZeeLumErr.at(bin-1) 
		<< std::setw(10) << vZmmLumErr.at(bin-1) 
		<< std::setw(10) << vZttLumErr.at(bin-1) 
		<< std::setw(10) << vWenLumErr.at(bin-1) 
		<< std::setw(10) << vWmnLumErr.at(bin-1) 
		<< std::setw(10) << vWtnLumErr.at(bin-1) << std::endl; 
	}


	// JES and QCD mixture 
	// here I am trying to get two things done at once. 1. JES syst 2. fake pho fraction syste using pho mc and qcd

	TH1F *hist_err = (TH1F*) halojet->Clone("hist_err");					//central values and use for SYSTEMATIC ERROR BAND
  																							// this must be identical to the final stacked plot
	hist_err->Add(zeejet);
	hist_err->Add(zmmjet);
	hist_err->Add(zttjet);
	hist_err->Add(wmnjet);
	hist_err->Add(wenjet);
	hist_err->Add(wtnjet);
	hist_err->Add(cosmicjet);

	if (which == "Met") hist_err->Add(diphojet);

	std::vector<float> qcdmcMixErr;		//this will hold the systematics for either QCD mthd 1 or 2 for a given case
													//to get systematics when 100% QCD is used
													
	std::vector<float> IdErr, vSidebandRwtErr;
	//init this so no prob for QCDerrMtd==1
	for (int bin =1; bin <= hist_err->GetNbinsX(); ++bin) 
	{
		qcdmcMixErr.push_back(0);
		vSidebandRwtErr.push_back(0);
	}

	//choose the sideband systematic methods for
	//1. QCDerrMtd == 2 (100% sideband) -> only QCD ID vary error -> iSideband_syst_mtd = 1
	//2. QCDerrMtd == 2 (100% reweighted sideband) -> QCD ID vary error + rewght error -> iSideband_syst_mtd = 2
	//3. QCDerrMtd == 1 (phomc+sideband) -> only QCD ID vary error -> iSideband_syst_mtd = 1

	int iSideband_syst_mtd = -1;
	if (iSideband == iNOMINAL && QCDerrMtd == 2) iSideband_syst_mtd = 0;
	if (iSideband == iISOADDEDREWGT && QCDerrMtd == 1) iSideband_syst_mtd = 2;


	//QCD/MC mix error should be calcualted for both deafault method and reweighed method
	//but not for the NJET plot in default method
	std::cout <<__LINE__ <<  "Deriving error for true photon fraction " << std::endl;
	if ( ! (iSideband == iNOMINAL && QCDerrMtd == 1 && which == "NJet"))
	{
		TH1F* hist_qcdcopy = (TH1F*) qcdjet->Clone("hist_qcdcopy");
		TH1F* hist_qcdVaried = (TH1F*) qcdjet->Clone("hist_qcdVaried");
		TH1F* hist_mcphoVaried = (TH1F*) mcphojet->Clone("hist_mcphoVaried");
		std::cout << "input hists before saling " << std::endl;
		PrintAll(hist_qcdcopy);
		PrintAll(hist_qcdVaried);
		PrintAll(hist_mcphoVaried);
		
		std::cout << __LINE__ << " fake pho fraction + 1sigma= " << fake_pho_frac_p << std::endl;
		float mcphoMixVaried_Norm = (leftover * (1 - fake_pho_frac_p) )/mcphojet->Integral("width");
		hist_mcphoVaried->Scale(mcphoMixVaried_Norm);
		
		float qcdMixVaried_Norm = (leftover * fake_pho_frac_p)/qcdjet->Integral("width");
		hist_qcdVaried->Scale(qcdMixVaried_Norm);

		std::cout << __LINE__ << "hists after saling " << std::endl;
		PrintAll(hist_qcdcopy);
		PrintAll(hist_qcdVaried);
		PrintAll(hist_mcphoVaried);

		//nominal
		hist_qcdcopy->Add(mcphojet);
		//varied
		hist_qcdVaried->Add(hist_mcphoVaried);

		std::cout << "hist_qcdVaried after scaling and addding hist_qcdcopy" << std::endl;
		PrintAll(hist_qcdVaried);

		
		qcdmcMixErr.clear();
		std::cout << "qcdmcMixErr " << std::endl;
		std::cout << std::setw(5) << "bin" << std::setw(10) << "err" << std::endl;
		for (int i=1; i <= hist_qcdcopy->GetNbinsX(); i++)
		{
			float err = fabs(hist_qcdcopy->GetBinContent(i) - hist_qcdVaried->GetBinContent(i));
			qcdmcMixErr.push_back(err);
			std::cout << std::setw(5) << i << std::setw(10) << err << std::endl;
		}
	}


	if (QCDerrMtd == 1)		//_______________________ plots that use 100% qcd
  	{
		TH1F* hist_qcd100Syst = (TH1F*) qcdjet_100->Clone("hist_qcd100_clone");	//this is for default method, njet case

		hist_err->Add(qcdjet_100);
	
		std::string abspath("Hist/"+path);
		if ( ! (iSideband == iNOMINAL && which == "NJet"))
		{
			GetQCD100Err(IdErr, qcdjet, name, abspath, 
					xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, 
					iSideband_syst_mtd, "", false);

		} else {
			GetQCD100Err(IdErr, hist_qcd100Syst, name, abspath, 
					xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, 
					iSideband_syst_mtd, "", false);
		}

		if ( !(iSideband == iNOMINAL && (which =="Met" || which == "NJet")))
		{
			GetSidebandReweighedError(vSidebandRwtErr, qcdjet_100, qcd100Norm, iSample, name, abspath, 
				xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4); 
		}
 	
	} else if (QCDerrMtd == 2)  //_____________________ plots that use 30%70%
	{
  		//hist_varyMcQcdMixSyst = (TH1F*) hist_err->Clone("hist_varyMcQcdMixSyst");	//use this to get a varied mix of qcd+pho mc
		
		hist_err->Add(qcdjet);
		hist_err->Add(mcphojet);
 
		//this is new systamatic for the 
    	std::string abspath("Hist/"+path);
	   GetQCD100Err(IdErr, qcdjet, name, abspath, 
				xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, 
				iSideband_syst_mtd, brname,
				false);
  	} else 
	{
		std::cerr << "bad value: QCDerrMtd=" << QCDerrMtd << std::endl;
      return;
	};


	std::cout << __LINE__ << "hist_err = Sum of all backgrounds after scaling " << std::endl;
	PrintAll(hist_err);


	//this to reweight the sideband according to the method 2 jet eta
	//total background to data. 03-19-2010
	//MakeJetEtaWeightsForSideband(QCDerrMtd, hist_err, phojet, brname, jets);
	//return;

	// NOW GET THE SYSTEMATICS FROM JES BY comparing central values with the JES UP/DOWN
	// USE THE HIST USED TO PUT THE ERROR BAND AS THE CENTRAL VALUE

	//to get JES error for all the MC samples, they must be added separately, centra,jesup, and jesdown,
	//and then take the maximum difference. this must be done after each hist is correctly scaled.
	//
	//<1> njet plot in default method: should EXCLUDE the pho mc contribution
	//<2> met  plot in default method: should INCLUDE the dipho mc contribution
	

	std::cout << __LINE__ << " ::Deriving JES Error" << std::endl;
	std::cout << std::setw(5) << "bin" << std::setw(10) << "err " << std::endl;
	std::vector<float> JESerr;
	for (int bin=1; bin <= zeejet->GetNbinsX(); bin++)
	{
		double jes_cen = wenjet->GetBinContent(bin) + wmnjet->GetBinContent(bin) + wtnjet->GetBinContent(bin)
							+ zeejet->GetBinContent(bin) + zmmjet->GetBinContent(bin) + zttjet->GetBinContent(bin);

		double jes_up = wenjetJESUP->GetBinContent(bin) + wmnjetJESUP->GetBinContent(bin) + wtnjetJESUP->GetBinContent(bin)
							+ zeejetJESUP->GetBinContent(bin) + zmmjetJESUP->GetBinContent(bin) + zttjetJESUP->GetBinContent(bin);

		double jes_down = wenjetJESDOWN->GetBinContent(bin) + wmnjetJESDOWN->GetBinContent(bin) + wtnjetJESDOWN->GetBinContent(bin)
							+ zeejetJESDOWN->GetBinContent(bin) + zmmjetJESDOWN->GetBinContent(bin) + zttjetJESDOWN->GetBinContent(bin);

		//FOR NJET PLOT IN DEFAULT METHOD: NO PHO MC contribution
		if (!(iSideband == iNOMINAL && QCDerrMtd == 1 && which == "NJet"))
		{
			jes_cen  += mcphojet->GetBinContent(bin);
			jes_up   += mcphojetJESUP->GetBinContent(bin);
			jes_down += mcphojetJESDOWN->GetBinContent(bin);
		}

		//for MET plots in DEFAULT METHOD: ADD DI PHO MC contribution
		if (iSideband == iNOMINAL && QCDerrMtd == 2 && which == "Met")
		{
			//assert (false && "di pho mc stuff not complte yet !");
			
			jes_cen += diphojet->GetBinContent(bin);
			jes_up  += diphojetJESUP->GetBinContent(bin);
			jes_down + diphojetJESDOWN->GetBinContent(bin);
		}
		
		const float err = GetMax(jes_cen, jes_up, jes_down);
		JESerr.push_back(err);
		std::cout << std::setw(5) << bin << std::setw(10) << JESerr.at(bin-1) << std::endl;

	}

		// NOW GET THE STAT ERROR
	std::vector<float> statErr;
	std::cout << __LINE__ <<  "::Deriving stat Error from hist_err" << std::endl;
	std::cout << std::setw(5) << "bin" << std::setw(10) << "err " << std::endl;
	for (int i=1; i <= hist_err->GetNbinsX(); i++)
	{
		statErr.push_back(hist_err->GetBinError(i));
		std::cout << std::setw(5) << i << std::setw(10) << statErr.at(i-1) << std::endl;
	}

	
	//collecting PDF/ISR/FSR/Q2/AlphaS errors ---------------------------------
	std::vector<float> vPDFerror, vISRFSRerror, vQ2error, vAlphaSerror;
	//set all values to zero so no need to worry if any of these errors
	//are not evaluated for a certain plot
	for (int bin = 1; bin<= hist_err->GetNbinsX(); ++bin) 
	{
		vPDFerror.push_back(0);
		vISRFSRerror.push_back(0);
		vQ2error.push_back(0);
		vAlphaSerror.push_back(0);

	}
	
	//if (brname.length()>0 && QCDerrMtd == 2)		//_______________________ plots using PHOTON MC
	// these error should be included in all cases except 
	// in qcderr mtd 2 : 100% QCD in default method (30% qcd, 70%pho mc) for NJET plot
	if (brname.length()>0 && ! (QCDerrMtd == 1 && iSideband == iNOMINAL))	//_______________________ plots using PHOTON MC
	{
		//now get PDF error, only when we use PHO MC
		TFile pdfSyst("PDFSyst_FitFunctions.root");
		assert (! pdfSyst.IsZombie() && "pdfSyst file not fouund!");

		TF1* tf1PDF = dynamic_cast<TF1*> (pdfSyst.Get(brname.c_str()));
		if (tf1PDF != NULL)
		{
			std::cout << __LINE__ << ":: Deriving PDF errors " << std::endl;
			std::cout << "opened file "; pdfSyst.Print();
			std::cout << "function "; tf1PDF->Print();
			
			std::cout << std::setw(5) << "bin" << std::setw(10) << "err at bin center " << std::endl;
			for (int bin = 1; bin<= hist_err->GetNbinsX(); ++bin)
			{
				float fPe = tf1PDF->Eval(hist_err->GetBinCenter(bin));
				if (fPe<0.01) fPe = 0.01;  //just like ISR/FSR and Q2 keep it 1% minimal
				vPDFerror.at(bin-1) = fPe * mcphojet->GetBinContent(bin);
				std::cout << std::setw(5) << bin << std::setw(10) << vPDFerror.at(bin-1) << std::endl;
			}
			delete tf1PDF;

		} else 
		{
			std::cout << __LINE__ << ":: WARNING !  NO PDF Systematic function found for " << brname  << std::endl;
		}


		// ISR/FSR syst
		TFile isrfsrSyst("ISRFSRsyst_FitFunctions.root");
		assert (! isrfsrSyst.IsZombie() && "isrfsrSyst file not fouund!");

		TF1* tf1ISR = dynamic_cast<TF1*> (isrfsrSyst.Get(brname.c_str()));
		if (tf1ISR != NULL)
		{
			std::cout << __LINE__ << ":: Deriving ISR errors " << std::endl;
			std::cout << "opened file "; isrfsrSyst.Print();
			std::cout << "function "; tf1ISR->Print();
			std::cout << std::setw(5) << "bin" << std::setw(10) << "err at bin center " << std::endl;
			for (int bin = 1; bin<= hist_err->GetNbinsX(); ++bin)
			{
				float fRe = tf1ISR->Eval(hist_err->GetBinCenter(bin));
				if (fRe<0.01) fRe = 0.01;	//keep minimum error for this at 1%. fits are not perfect elog#1482
				vISRFSRerror.at(bin-1) =  fRe * mcphojet->GetBinContent(bin);
				std::cout << std::setw(5) << bin << std::setw(10) << vISRFSRerror.at(bin-1) << std::endl;
			}
			delete tf1ISR;
		} else 
		{
			std::cout << __LINE__ << ":: WARNING !  NO ISR/FSR Systematic function found for "  << brname  << std::endl;
		}

		// Q2 syst
		TFile q2Syst("Q2Syst_FitFunctions.root");
		assert (! q2Syst.IsZombie() && "q2Syst file not fouund!");

		TF1* tf1Q2 = dynamic_cast<TF1*> (q2Syst.Get(brname.c_str()));
		if (tf1Q2 != NULL)
		{
			std::cout << __LINE__ << ":: Deriving Q2 errors " << std::endl;
			std::cout << "opened file "; q2Syst.Print();
			std::cout << "function "; tf1Q2->Print();
			
			std::cout << std::setw(5) << "bin" << std::setw(10) << "err at bin center " << std::endl;
			for (int bin = 1; bin<= hist_err->GetNbinsX(); ++bin)
			{
				float fQe = tf1Q2->Eval(hist_err->GetBinCenter(bin)); 
				if (fQe <0.01) fQe = 0.01;	//keep minimum error for this at 1%. fits are not perfect elog#1481,1480
				vQ2error.at(bin-1) = fQe * mcphojet->GetBinContent(bin);
				std::cout << std::setw(5) << bin << std::setw(10) << vQ2error.at(bin-1) << std::endl;
			}
			delete tf1Q2;
		} else 
		{
			std::cout << __LINE__ << ":: WARNING !  NO Q2 Systematic function found for "  << brname  << std::endl;
		}


		// alpha_s syst
		std::cout << __LINE__ << ":: Deriving AlphaS errors" << std::endl;
		GetAlphaS_Systematic(brname,vAlphaSerror,mcphojet);

		
	} //end of collecting PDF/ISR/FSR/Q2/ALPHA-S errors



	
 	// NOW COLLECT ALL ERRORS AND PUT THEM TOGETHER FOR ONE FINAL NUMBER
	std::vector<float> ERRORS;

	std::cout << __FUNCTION__ << ":" << __LINE__ << "::All Errors" << std::endl;
	std::cout << setw(4) << "bin"<< setw(8) << "JES"
				 << setw(8) << "PhoFrac" << setw(8) << "EWK Lum" 
				 << setw(8) << "Cosmic" << setw(6) << "Halo"
				 << setw(8) << "Stat" << setw(10) << "QCD Id" 
				 << setw(10) << "PDF" 
				 << setw(10) << "Q^2" 
				 << setw(10) << "IsrFsr" 
				 << setw(10) << "AlphaS" 
				 << setw(10) << "REWGT" 
				 << setw(10) << "Tot Err" 
				 << setw(10) << "Sum BG" 
				 << setw(10) << "Data" 
				 << setw(10) << "D-B/B" 
				 <<  std::endl; 
	
	assert ( (cosmicErr.size() == haloErr.size() 
				&& haloErr.size()==qcdmcMixErr.size() 
				&& qcdmcMixErr.size()==JESerr.size()
				&& JESerr.size() == statErr.size() 
				&& statErr.size() == cosmicErr.size() 
				&& cosmicErr.size() == vPDFerror.size()
				&& vPDFerror.size() == vISRFSRerror.size()
				&& vISRFSRerror.size() == vQ2error.size()
				&& vQ2error.size() == vAlphaSerror.size()
				&& vAlphaSerror.size() == IdErr.size()
				&& IdErr.size() == vZeeLumErr.size()
				&& vZmmLumErr.size() == vZeeLumErr.size()
				&& vZttLumErr.size() == vZeeLumErr.size()
				&& vWenLumErr.size() == vWtnLumErr.size()
				&& vWenLumErr.size() == vWmnLumErr.size()
				&& vWmnLumErr.size() == vSidebandRwtErr.size()
				) 
			&& "MakeHistLogAndRatio_debug::Error vectors are not same in size!");


	std::vector<float> vTotalEWKLumErr;

	for (unsigned int i=0; i < vZeeLumErr.size(); i++)
	{
		vTotalEWKLumErr.push_back(sqrt (pow(vZeeLumErr[i],2) 
						+ pow(vZmmLumErr[i],2)+ pow(vZttLumErr[i],2) + pow(vWenLumErr[i],2) 
						+ pow(vWmnLumErr[i],2) + pow(vWtnLumErr[i],2)));
	}
	
	
	for (unsigned int i=0; i < cosmicErr.size(); i++)
	{
		float sum = pow(cosmicErr[i],2) 
						+ pow(haloErr[i],2)
						+ pow(vTotalEWKLumErr[i],2) 
						+ pow(qcdmcMixErr[i],2)
						+ pow(IdErr.at(i),2) 
						+ pow(JESerr[i],2) 
						+ pow(statErr[i],2)
						+ pow(vPDFerror.at(i),2)
						+ pow(vISRFSRerror.at(i),2)
						+ pow(vQ2error.at(i),2)
						+ pow(vAlphaSerror.at(i),2)
						+ pow(vSidebandRwtErr.at(i),2)
						;

		ERRORS.push_back(sqrt(sum));

		// now set the sum of background hists error to this total
		if (hist_err->GetBinContent(i+1)) //no errors if there are no background prediction!?
		{
			hist_err->SetBinError(i+1, ERRORS.at(i));
		}
		
		std::cout << setiosflags(ios::fixed) << setprecision(1) 
			<< setw(4) << i+1 << setw(8) << JESerr.at(i)
			<< setw(8) << qcdmcMixErr.at(i) << setw(8) << vTotalEWKLumErr.at(i) 
			<< setw(8) << cosmicErr.at(i) << setw(6) << haloErr.at(i)
			<< setw(8) << statErr.at(i) << setw(10) << IdErr.at(i) 
			<< setw(10) << vPDFerror.at(i) 
			<< setw(10) << vQ2error.at(i) 
			<< setw(10) << vISRFSRerror.at(i) 
			<< setw(10) << vAlphaSerror.at(i)
			<< setw(10) << vSidebandRwtErr.at(i)
			<< setw(10) << ERRORS.at(i) 
			<< setw(10) << hist_err->GetBinContent(i+1) 
			<< setw(10) << phojet->GetBinContent(i+1) 
			<< setw(10) << setiosflags(ios::fixed) << setprecision(2);
			if (hist_err->GetBinContent(i+1)>0) 
			{
				std::cout << (phojet->GetBinContent(i+1) - hist_err->GetBinContent(i+1)) /hist_err->GetBinContent(i+1);
			} else
			{
				std::cout << 0.0;
			}
			std::cout <<  std::endl; 

	}



	//ONLY THE COSMATIC CHANGES TO THE HISTOGRAMS SHOULD BE DONE BELOW THIS POINT!
	
	if (iDRAWERROR_BREAKDOWN)
	{
	DrawErrorBreakDown(QCDerrMtd, iSideband, hist_err,JESerr,qcdmcMixErr, vTotalEWKLumErr, cosmicErr,haloErr, statErr,IdErr,
			vPDFerror, vQ2error,vISRFSRerror,vAlphaSerror, vSidebandRwtErr, ERRORS, title, jets, which);
	}

	TCanvas *c1 = new TCanvas("HIST","HIST",600,800);
	TPad *pLOG = new TPad("Plog","Log Plot",0.01,0.4,0.99,0.99);   //in NDC cdts
	TPad *pRAT = new TPad("Pratio","Ratio Plot",0.01,0.01,0.99,0.4);  //in NDC cdts
	pLOG->SetBorderMode(0);
	pLOG->SetBorderSize(0);
	pLOG->SetBottomMargin(0);
	pRAT->SetBorderMode(0);
	pRAT->SetBorderSize(0);
	pRAT->SetTopMargin(0);
	pRAT->SetBottomMargin(0.1);
	//pLOG->SetFillColor(kBlue);
	//pRAT->SetFillColor(kCyan);
	//pRAT->Draw();
	//MakeRatioPlot(which, jets, phojet,hist_err, ERRORS, pRAT, false);		//drawn at the end
	//MakeRatioPlot(which, jets, phojet,hist_err, ERRORS,false);
	//pRAT->SetEditable(0);
	c1->cd();
	pLOG->Draw();
	pLOG->cd();

	


  //************* SET HIST COLORS ***************************************

	Int_t linecolor=47;

  mcphojet->SetLineColor(linecolor);
  mcphojet->SetFillColor(6);
  qcdjet->SetLineColor(linecolor);
  qcdjet->SetFillColor(kYellow);
  cosmicjet->SetLineColor(linecolor);
  cosmicjet->SetFillColor(kGreen);
  halojet->SetLineColor(linecolor);
  halojet->SetFillColor(12);
  int ewkColor = 29;
  zeejet->SetLineColor(linecolor);
  zeejet->SetFillColor(ewkColor);
  zmmjet->SetLineColor(ewkColor);
  zmmjet->SetFillColor(ewkColor);
  zttjet->SetLineColor(ewkColor);
  zttjet->SetFillColor(ewkColor);
  wenjet->SetLineColor(ewkColor);
  wenjet->SetFillColor(ewkColor);
  wmnjet->SetLineColor(ewkColor);
  wmnjet->SetFillColor(ewkColor);
  wtnjet->SetLineColor(ewkColor);
  wtnjet->SetFillColor(ewkColor);

  if (which == "Met")
  {
	  //diphojet->SetFillColor(kGreen+4);
	  //diphojet->SetLineColor(kGreen+4);
	  diphojet->SetFillColor(kCyan);
	  diphojet->SetLineColor(linecolor);
  }

  
  phojet->SetLineColor(kBlack);
  phojet->SetMarkerStyle (8);

  //this is the line histo showing 100% of QCD 0% MC pho
  qcdjet_100->SetLineColor(kRed);
  qcdjet_100->SetFillColor(kYellow);
  hist_err->SetFillStyle(3001);
  hist_err->SetFillColor(13);
  hist_err->SetLineColor(10);
 

  //********* MAKE THE PLOT ***************************************
  THStack *hs = new THStack ("hs", NULL);

  if (QCDerrMtd == 1) {//plots that use 100% qcd
    hs->Add(halojet);
    hs->Add(cosmicjet);
    hs->Add(zmmjet);
    hs->Add(wtnjet);
    hs->Add(zttjet);
    hs->Add(wmnjet);
    hs->Add(wenjet);
    hs->Add(zeejet);
    hs->Add(qcdjet_100);

	 std::cout <<  __FUNCTION__ << ":" << __LINE__ << ": QCD100 METHOD" << std::endl;

	/* int bin = 2;
	 int prec = 3;
	 DumpHistBin(halojet, bin, 1, prec);
	 DumpHistBin(cosmicjet, bin, 0, prec);
	 DumpHistBin(zeejet, bin, 0, prec);
	 DumpHistBin(zmmjet, bin, 0, prec);
	 DumpHistBin(zttjet, bin, 0, prec);
	 DumpHistBin(wenjet, bin, 0, prec);
	 DumpHistBin(wmnjet, bin, 0, prec);
	 DumpHistBin(wtnjet, bin, 0, prec);
	 DumpHistBin(qcdjet_100, bin, 0, prec);
	 DumpHistBin(phojet, bin, 0, prec);
	 DumpHistBin(hist_err, bin, 0, prec);
	 */
	 std::cout <<  "\n\n" << std::endl;

  } else if (QCDerrMtd == 2) {//plots that use 30%70%
    hs->Add(halojet);
    hs->Add(cosmicjet);
    hs->Add(zmmjet);
    hs->Add(wtnjet);
    hs->Add(zttjet);
    hs->Add(wmnjet);
    hs->Add(wenjet);
    hs->Add(zeejet);
	 if (which == "Met") hs->Add(diphojet);
    hs->Add(qcdjet);
    hs->Add(mcphojet);

	std::cout <<  __FUNCTION__ << ":" << __LINE__ << ": MIX METHOD" << std::endl;
	
/*	int bin = 11;
	int prec = 3;
	DumpHistBin(halojet, bin, 1, prec);
	DumpHistBin(cosmicjet, bin, 0, prec);
	DumpHistBin(zeejet, bin, 0, prec);
	DumpHistBin(zmmjet, bin, 0, prec);
	DumpHistBin(zttjet, bin, 0, prec);
	DumpHistBin(wenjet, bin, 0, prec);
	DumpHistBin(wmnjet, bin, 0, prec);
	DumpHistBin(wtnjet, bin, 0, prec);
	DumpHistBin(qcdjet, bin, 0, prec);
	DumpHistBin(mcphojet, bin, 0, prec);
	DumpHistBin(phojet, bin, 0, prec);
	DumpHistBin(hist_err, bin, 0, prec);
*/	

	std::cout <<  "\n\n" << std::endl;

	 
  } else {
    std::cerr << "bad value: QCDerrMtd=" << QCDerrMtd << std::endl;
    return;
  };

	
  const float fLegend_x1=0.48;
  const float fLegend_y1=0.6;
  const float fLegend_x2=0.9;
  const float fLegend_y2=0.9;
  TLegend *leg = new TLegend (fLegend_x1,fLegend_y1,fLegend_x2,fLegend_y2,"","NDC");
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  std::string str_pho,str_cosmic,str_halo,str_zee,str_zmm,str_ztt,str_wen,str_wmn,str_wtn;
  std::string str_diphomc;

  if (jets==1) {
    //str_pho 	= "Data (#gamma + #geq 1 Jet + #slash{E}_{T}>25GeV";	
    str_pho 	= "Data (#gamma + #geq 1 Jet)";	
    str_cosmic = "#gamma^{cosmic} + #geq 1 Jet";
    str_zee 	= "Z->ee MC (e + #geq 1 Jet)";
    str_wen 	= "W->e#nu MC (e + #geq 1 Jet)";
    str_ztt 	= "Z->#tau#tau MC (e + #geq 1 Jet)";
    str_wtn 	= "W->#tau#nu MC (e + #geq 1 Jet)";
    str_wmn 	= "W->#mu#nu MC (e + #geq 1 Jet)";
    str_zmm 	= "Z->#mu#mu MC (e + #geq 1 Jet)";
    str_halo 	= "#gamma^{halo} + #geq 1 Jet";
    str_diphomc= "Di-#gamma MC + #geq 1 Jet";
  }
  if (jets==2) {
    str_pho 	= "Data (#gamma + #geq 2 Jets)";	
    str_cosmic = "#gamma^{cosmic} + #geq 2 Jets";
    str_zee 	= "Z->ee MC (e + #geq 2 Jets)";
    str_wen 	= "W->e#nu MC (e + #geq 2 Jets)";
    str_ztt 	= "Z->#tau#tau MC (e + #geq 2 Jets)";
    str_wtn 	= "W->#tau#nu MC (e + #geq 2 Jets)";
    str_wmn 	= "W->#mu#nu MC (e + #geq 2 Jets)";
    str_zmm 	= "Z->#mu#mu MC (e + #geq 2 Jets)";
    str_halo 	= "#gamma^{halo} + #geq 2 Jets";
    str_diphomc= "Di-#gamma MC + #geq 2 Jets";
  }
	

  if (QCDerrMtd == 1) {//plots that use 100% qcd
    leg->AddEntry(phojet,str_pho.c_str());
    leg->AddEntry(qcdjet_100, "QCD (100% #gamma sideband)");
    leg->AddEntry(wtnjet,"EWK MC");
    leg->AddEntry(cosmicjet,str_cosmic.c_str());
    leg->AddEntry(halojet,str_halo.c_str());
    leg->AddEntry(hist_err,"Systematics Uncertainty");

  } else if (QCDerrMtd == 2) {//plots that use 30%70%
    leg->AddEntry(phojet,str_pho.c_str());
    leg->AddEntry(mcphojet,mc_str.c_str());
    leg->AddEntry(qcdjet,qcd_str.c_str());
    leg->AddEntry(wtnjet,"EWK MC");
    leg->AddEntry(cosmicjet,str_cosmic.c_str());
    leg->AddEntry(halojet,str_halo.c_str());
	 if (which == "Met") leg->AddEntry(diphojet,str_diphomc.c_str());
    leg->AddEntry(hist_err,"Systematic Uncertainty");
  } else {
    std::cerr << "bad value: QCDerrMtd=" << QCDerrMtd << std::endl;
    return;
  };

	leg->SetBorderSize (1);
	leg->SetFillColor (10);
	//gStyle->SetCanvasColor (10);
	//gStyle->SetCanvasBorderSize (0);
	//gStyle->SetCanvasBorderMode (0);

	//gStyle->SetPadColor (10);
	//gStyle->SetFillColor (10);
	//gStyle->SetTitleFillColor (10);
	//gStyle->SetTitleBorderSize (1);
	//gStyle->SetStatColor (10);
	//gStyle->SetStatBorderSize (1);


	
	//dynamically adjust the y-scale
	double dYscale_max = hist_err->GetBinContent(hist_err->GetMaximumBin()) * 10.;
	double dYscale_min = 1e10;
	for (int ibin=1; ibin <= hist_err->GetNbinsX();++ibin)
	{
		if (hist_err->GetBinContent(ibin) <= 0) continue; //some bins can be slightly <0 after subtracting other backgrounds.
		if (hist_err->GetBinContent(ibin) < dYscale_min)
		{
			dYscale_min = hist_err->GetBinContent(ibin);
		}
	}
	dYscale_min *= 1/10.;
	std::cout <<  __LINE__ << ": Y-scale min, max = " << dYscale_min << ", " << dYscale_max << std::endl;

	
	//new TCanvas();
	gStyle->SetOptStat(0);
	gStyle->SetTextFont(132);
	gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();

	//hs->SetMinimum(dYscale_min);
	hs->SetMinimum(0.005);
	hs->SetMaximum(dYscale_max);
	

	//original method. put this back after debug --11-22-2009
	hs->Draw("HIST");		//need this as I am calling sumw2 for all hists. if not it will draw all hists with error bars
	//new TCanvas();
	//hs->Draw("HIST nostack,e1p");		//need this as I am calling sumw2 for all hists. if not it will draw all hists with error bars

	hs->GetXaxis()->SetTitle (title.c_str());
	std::ostringstream ytitle;
	//if (which == "PhotonEt") ytitle << "#Delta N/ #Delta E_{T}";
	//if (which == "InvMass") ytitle << "#Delta N/ #Delta M";
	std::string xtitle, yunits;
	if (which =="PhotonEt") { xtitle += "E_{T}^{ #gamma} (GeV)"; yunits= "GeV"; }
	if (which =="InvMass") { xtitle += "Invariant Mass(#gamma, Lead Jet) (GeV/c^{2})"; yunits= "GeV/c^{2}"; }
	if (jets ==2 && which =="PhoJetsInvMass") { xtitle += "Invariant Mass(#gamma, Two Lead Jets) (GeV/c^{2})"; yunits= "GeV/c^{2}"; }
	if (which =="Ht") xtitle += "H_{T} (GeV)";
	if (which =="JetsInvMass") { xtitle += "Invariant Mass(Two Lead Jets) (GeV/c^{2})"; yunits= "GeV/c^{2}"; }
	if (which =="LeadJetEt") xtitle += "E_{T}^{Lead Jet} (GeV)";
	if (which =="SecondLeadJetEt") xtitle += "E_{T}^{Second Lead Jet} (GeV)";

	if (which == "NJet") ytitle << "Events";
	else ytitle << "Events / " << phojet->GetBinWidth(1) << " " << yunits;
	hs->GetYaxis()->SetTitle (ytitle.str().c_str());

	int labelfont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	int titlefont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	hs->GetXaxis()->SetLabelFont(labelfont);
	hs->GetYaxis()->SetLabelFont(labelfont);
	//hs->GetYaxis()->SetLabelSize(0.05);
	//hs->GetXaxis()->SetLabelSize(0.05);
	hs->GetXaxis()->SetTitleFont(titlefont);
	hs->GetYaxis()->SetTitleFont(titlefont);
	//hs->GetYaxis()->SetTitleSize(0.05);
	//hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetXaxis()->SetTitleOffset(0.9);
	hs->GetYaxis()->SetTitleOffset(0.9);
	hs->GetXaxis()->CenterTitle(true);
	hs->GetYaxis()->CenterTitle(true);


	hist_err->SetDrawOption("HIST");
	hist_err->Draw("sameE2");
	phojet->Draw("same");
	leg->Draw();

  
	TPaveText *tp = new TPaveText(0.5,0.91,0.9,0.99,"NDC");
	tp->SetLineColor(10);
	tp->SetBorderSize(0);
	tp->SetTextFont(titlefont);
	tp->AddText("CDF Run II Preliminary 2.0 fb^{-1}");
	tp->Draw();

	std::string sLegendTitle("");
	if (iSample == iG30) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+#geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2}";
	else if (iSample == iG40) sLegendTitle = "#gamma^{E_{T}>40GeV}_{|#eta|<1.1}+#geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2}";
	else if (iSample == iG30MET25) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1} + #geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2} + #slash{E}_{T}>25GeV";
	else if (iSample == iG30MET40) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1} + #geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2} + #slash{E}_{T}>40GeV";
	else if (iSample == iG30EXCL1J) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+ == 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2}";
	TPaveText *tpTitle = new TPaveText(0.1,0.91,0.5,0.99,"NDC");
	tpTitle->SetLineColor(iWHITE);
	tpTitle->SetBorderSize(0);
	tpTitle->SetTextFont(iTITLE_FONT);
	tpTitle->AddText(sLegendTitle.c_str());
	tpTitle->Draw();


	const float fTp2_h = 0.15;
	TPaveText *tp2 = new TPaveText(fLegend_x1,fLegend_y1-fTp2_h,fLegend_x2,fLegend_y1,"NDC");
	tp2->SetLineColor(kBlack);
	tp2->SetBorderSize(1);
	//tp->SetTextFont();
	if (iSideband == iNOMINAL) tp2->AddText("Nominal sideband");
	else if (iSideband == iISOADDED) tp2->AddText("Iso added sideband");
	else if (iSideband == iISOADDEDREWGT) 
	{
		//tp2->AddText("Iso added sideband");
		//tp2->AddText("reweighed to total background");
		tp2->AddText("Sideband reweighed");
		tp2->AddText("to background");
	}
	tp2->Draw();
	

	bool mark_overflow = false;
  	if (QCDerrMtd == 1) {//plots that use 100% qcd
    	if ( halojet->GetBinContent(halojet->GetNbinsX()) ||
				 cosmicjet->GetBinContent(cosmicjet->GetNbinsX()) ||
				 zmmjet->GetBinContent(zmmjet->GetNbinsX()) ||
				 wtnjet->GetBinContent(wtnjet->GetNbinsX()) ||
				 zttjet->GetBinContent(zttjet->GetNbinsX()) ||
				 wmnjet->GetBinContent(wmnjet->GetNbinsX()) ||
				 wenjet->GetBinContent(wenjet->GetNbinsX()) ||
				 zeejet->GetBinContent(zeejet->GetNbinsX()) ||
				 qcdjet_100->GetBinContent(qcdjet_100->GetNbinsX()) )
			mark_overflow = true;

  } else if (QCDerrMtd == 2) {//plots that use 30%70%
    if ( halojet->GetBinContent(halojet->GetNbinsX()) ||
				 cosmicjet->GetBinContent(cosmicjet->GetNbinsX()) ||
				 zmmjet->GetBinContent(zmmjet->GetNbinsX()) ||
				 wtnjet->GetBinContent(wtnjet->GetNbinsX()) ||
				 zttjet->GetBinContent(zttjet->GetNbinsX()) ||
				 wmnjet->GetBinContent(wmnjet->GetNbinsX()) ||
				 wenjet->GetBinContent(wenjet->GetNbinsX()) ||
				 zeejet->GetBinContent(zeejet->GetNbinsX()) ||
				 qcdjet->GetBinContent(qcdjet->GetNbinsX()) ||
				 mcphojet->GetBinContent(mcphojet->GetNbinsX()) )
			mark_overflow = true;
  }
  
  if (mark_overflow)
  {
		gPad->Range(0,0,1,1);

		TLine *line = new TLine();
		line->SetLineWidth(2);
		line->DrawLineNDC(0.895,0.218,0.93,0.218);
		
		TText *tt = new TText(0.95,0.15,"Overflow bin");
			tt->SetTextFont(titlefont);
			tt->SetTextAngle(90);
			tt->SetTextSize(0.035);
			tt->SetTextColor(kBlue);
			tt->SetNDC();

		tt->Draw();
  }


  c1->cd();
  pRAT->Draw();
  MakeRatioPlot(which, jets, phojet,hist_err, ERRORS, pRAT, false);
	c1->cd();


	//DrawErrorBreakDown(QCDerrMtd, hist_err,JESerr,qcdmcMixErr,cosmicErr,haloErr, statErr,IdErr,
	//		vPDFerror, vQ2error,vISRFSRerror,vAlphaSerror,ERRORS, title, jets, which);
	

	std::cout << " >>>>>>> Settings used for " << jets << " jet/s " << which << " plot <<<<<<" << std::endl;
	std::string sSampleName("UNKNOWN");
	if (iSample == iG30) sSampleName = sG30;
	else if (iSample == iG40) 		 sSampleName = sG40;
	else if (iSample == iG30MET25) sSampleName = sG30MET25;
	else if (iSample == iG30MET40) sSampleName = sG30MET40;
	else if (iSample == iG30EXCL1J) sSampleName = sG30EXCL1J;

	std::string sSideband("UNKNOWN");
	if (iSideband == iNOMINAL) sSideband = "NOMINAL";
	else if (iSideband == iISOADDED) sSideband = "ISO ADDED";
	else if (iSideband == iISOADDEDREWGT) sSideband = "ISO ADDED+REWEIGHED";
	
	std::string sMethod("UNKNOWN");
	if (QCDerrMtd == 1) sMethod = "100% sideband";
	else if (QCDerrMtd == 2) sMethod = "70% MC + 30% sideband";
	

	//these are the true interagral values.
	//once the hists are rebinned their intergral changes from the original -02-17-2010
	double dataint = 0, qcdint=0, mcphoint=0, ewkint=0;//, haloint=0, cosmicint=0;
	for (int bin=1; bin <= phojet->GetNbinsX(); ++bin) dataint += phojet->GetBinContent(bin) * phojet->GetBinWidth(bin);
	if (QCDerrMtd == 1)
	{
		for (int bin=1; bin <= qcdjet_100->GetNbinsX(); ++bin) qcdint += qcdjet_100->GetBinContent(bin) * qcdjet_100->GetBinWidth(bin);
	} else if (QCDerrMtd == 2)
	{
		for (int bin=1; bin <= qcdjet->GetNbinsX(); ++bin) qcdint += qcdjet->GetBinContent(bin) * qcdjet->GetBinWidth(bin);
		for (int bin=1; bin <= mcphojet->GetNbinsX(); ++bin) mcphoint += mcphojet->GetBinContent(bin) * mcphojet->GetBinWidth(bin);
	}
	
	for (int bin=1; bin <= zeejet->GetNbinsX(); ++bin)
	{
		ewkint += (zeejet->GetBinContent(bin) + 
				zmmjet->GetBinContent(bin)+
				zttjet->GetBinContent(bin)+
				wenjet->GetBinContent(bin)+
				wmnjet->GetBinContent(bin)+
				wtnjet->GetBinContent(bin)) * zeejet->GetBinWidth(bin);
	}

	//for (int bin=1; bin <= cosmicjet->GetNbinsX(); ++bin) cosmicint += cosmicjet->GetBinContent(bin) * cosmicjet->GetBinWidth(bin);
	//for (int bin=1; bin <= halojet->GetNbinsX(); ++bin) haloint += halojet->GetBinContent(bin) * halojet->GetBinWidth(bin);

	//cosmic/halo integrals are fine as they are scaled estimated values dividing by histogram integral.
	//data hist is rebinned and divided by bin width
	//ewk is "SCALED" by luminosity. so the integral of hist has no meaning now!!!
	//qcd/mcpho are "SCALED" to data integral after subtracting cosmic/halo/ewk backgrounds and integral has no meaning either!!!
	//the qcd/mcpho gets values from data hist which is already rebinned.
	
	//plot shown at some full status and preblessing talks had constant bin sizes and the first bin
	//content of photon et seems >10E6. but now <10E5 because I am using variable bin sizes that
	//require normalizing (dividing by bin width) - 2-17-2010
	
	std::cout << "iSAMPLE                = " << sSampleName << std::endl;
	std::cout << "Fake pho frac+/-sigma  = " << fake_pho_frac_d << "+/-" << fake_pho_frac_sigma << std::endl;
	std::cout << "Cosmic estimate        = " << cosmicEst << std::endl;
	std::cout << "METHOD                 = " << sMethod << std::endl;
	std::cout << "Halo estimate          = " << haloEst << std::endl;
	std::cout << "Sideband               = " << sSideband << std::endl; 
	std::cout << "DATA hist integral     = " << phojet->Integral()  << "->" << dataint << std::endl;
	std::cout << "Halojet Integral       = " << halojet->Integral("width") << std::endl;
	std::cout << "cosmicjet Integral     = " << cosmicjet->Integral("width") << std::endl;
	if (QCDerrMtd == 1)
	{
	std::cout << "Qcdjet_100 Integral    = " << qcdjet_100->Integral() << "->" << qcdint << std::endl;
	}
	else if (QCDerrMtd == 2)
	{
	std::cout << "Qcdjet Integral        = " << qcdjet->Integral() << "->" << qcdint << std::endl;
	std::cout << "MC PHO Integral        = " << mcphojet->Integral() << "->" << mcphoint << std::endl;
	}
	std::cout << "EWK Total Integral     = " << zeejet->Integral() + zmmjet->Integral()
		+ zttjet->Integral() + wenjet->Integral() + wmnjet->Integral() + wtnjet->Integral() 
		<< "->" << ewkint <<  std::endl;
	if (which == "Met") std::cout << "Di-pho Total Integral  = " << diphojet->Integral("width") << std::endl;


/*	std::cout << magenta << "hist_cosmic" << std::endl;
	double sum = 0;
	for (int bin = 1; bin <= cosmicjet->GetNbinsX(); ++bin)
	{
		if (cosmicjet->GetBinContent(bin)) //no errors if there are no background prediction!?
		{
			std::cout << bin << "\t" << cosmicjet->GetBinContent(bin) << std::endl;
			sum += cosmicjet->GetBinContent(bin) * cosmicjet->GetBinWidth(bin);
		}
	}
	std::cout << "my sum = " << sum <<  std::endl;
*/	
}

void MakeHistLogAndRatio_debug (int jets, std::string which, int icase=0)
//void MakeHistLogAndRatio_debug (int jets, std::string which)
{
	if (bDEBUG) std::cout << "Entering " << __FUNCTION__ << ":: " << __LINE__ << std::endl;
  	TVirtualPad *pad_old = gPad;	

	float xmin=0, xpoint1=0, xpoint2=0, xpoint3=0, xpoint4=0;
	float width1=0, width2=0, width3=0, width4=0;

	//default bin sizes for any sample
	
	if (jets == 1)
	{
		if (which == "Et_pho")           { xmin=0  ; xpoint1 = 150 ; xpoint2 = 250; xpoint3 = 300;xpoint4 = 650 ; width1= 10 ; width2= 20 ; width3= 50 ; width4= 250 ; iREBIN =5;} 
		else if (which == "Et_j1")       { xmin=0  ; xpoint1 = 200 ; xpoint2 = 250; xpoint3 = 300;xpoint4 = 600 ; width1= 10 ; width2= 10 ; width3= 50 ; width4= 200 ; iREBIN =5;}
		//else if (which == "InvMass_pj1") { xmin=50 ; xpoint1 = 500 ; xpoint2 = 600; xpoint3 = 700;xpoint4 = 1000; width1= 25 ; width2= 50 ; width3= 75 ; width4= 300 ; iREBIN =5;}
		else if (which == "InvMass_pj1") { xmin=0 ; xpoint1 = 500 ; xpoint2 = 600; xpoint3 = 700;xpoint4 = 1000; width1= 25 ; width2= 50 ; width3= 75 ; width4= 300 ; iREBIN =5;}
		else if (which == "Ht")          { xmin=0  ; xpoint1 = 200 ; xpoint2 = 600; xpoint3 = 800;xpoint4 = 1200; width1= 40 ; width2= 50 ; width3= 100; width4= 300 ; iREBIN =5;}
    	else if (which == "NJet")        { xmin=0  ; xpoint1 = 1   ; xpoint2 = 2  ; xpoint3 = 3  ;xpoint4 = 15  ; width1= 1  ; width2= 1  ; width3= 1  ; width4= 1   ; iREBIN =1;}
    	else if (which == "Met")         { xmin=0  ; xpoint1 = 60  ; xpoint2 = 100; xpoint3 = 200;xpoint4 = 450 ; width1= 10 ; width2= 20 ; width3= 50 ; width4= 150 ; iREBIN =1;}
    	else if (which == "DelR")        { xmin=0  ; xpoint1 = 1   ; xpoint2 = 2  ; xpoint3 = 4  ;xpoint4 = 6   ; width1= 0.2; width2= 0.2; width3= 0.2; width4= 0.2 ; iREBIN =5;}
    	else if (which == "DelEta")      { xmin=0  ; xpoint1 = 1   ; xpoint2 = 3  ; xpoint3 = 4  ;xpoint4 = 5   ; width1= 0.1; width2= 0.1; width3= 0.1; width4= 0.1 ; iREBIN =5;}
    	else if (which == "DelPhi")      { xmin=0  ; xpoint1 = 1   ; xpoint2 = 2  ; xpoint3 = 3  ;xpoint4 = 4   ; width1= 0.1; width2= 0.1; width3= 0.1; width4= 0.1 ; iREBIN =5;}
    	else if (which == "DetEta_pho")  { xmin=-1.5; xpoint1 = -0.5; xpoint2 = 0.5; xpoint3 = 1  ;xpoint4 = 1.5 ; width1= 0.1; width2= 0.1; width3= 0.1; width4= 0.1 ; iREBIN =1;}
    	else if (which == "DetPhi_pho")  { xmin=0  ; xpoint1 = 1   ; xpoint2 = 2  ; xpoint3 = 3  ;xpoint4 = 8   ; width1= 0.1; width2= 0.1; width3= 0.1; width4= 0.1 ; iREBIN =1;}
    	else if (which == "DetEta_j1")   { xmin=-3.5; xpoint1 = -1.5; xpoint2 = 0.5; xpoint3 = 1.5  ;xpoint4 = 3.5 ; width1= 0.2; width2= 0.2; width3= 0.2; width4= 0.2 ; iREBIN =1;}
    	//else if (which == "DetEta_j1")   { xmin=-3.5; xpoint1 = -1.5; xpoint2 = 0.5; xpoint3 = 1.5  ;xpoint4 = 3.5 ; width1= 0.1; width2= 0.1; width3= 0.1; width4= 0.1 ; iREBIN =1;}
    	else if (which == "DetPhi_j1")   { xmin=0  ; xpoint1 = 1   ; xpoint2 = 2  ; xpoint3 = 3  ;xpoint4 = 4   ; width1= 0.1; width2= 0.1; width3= 0.1; width4= 0.1 ; iREBIN =1;}
    	else if (which == "EvtPhi_j1")   { xmin=0  ; xpoint1 = 1   ; xpoint2 = 2  ; xpoint3 = 3  ;xpoint4 = 4   ; width1= 0.1; width2= 0.1; width3= 0.1; width4= 0.1 ; iREBIN =1;}
	}

	if (jets == 2) {
		if (which == "Et_pho") 				 	{ xmin=0  ; xpoint1= 150; xpoint2= 250; xpoint3= 300; xpoint4= 600 ; width1=10 ; width2=20 ; width3=50 ; width4=250 ; iREBIN =5;}
		else if (which == "Et_j1") 		 	{ xmin=0  ; xpoint1= 200; xpoint2= 260; xpoint3= 300; xpoint4= 600 ; width1=10 ; width2=20 ; width3=40 ; width4=200 ; iREBIN =5;}
		/*else if (which == "InvMass_pj1")  	{ xmin=100; xpoint1= 200; xpoint2= 300; xpoint3= 700; xpoint4= 1000; width1=10 ; width2=20 ; width3=50 ; width4=200 ; iREBIN =5;}
		else if (which == "InvMass_pj1j2")	{ xmin=100; xpoint1= 200; xpoint2= 300; xpoint3= 700; xpoint4= 1300; width1=10 ; width2=20 ; width3=50 ; width4=200 ; iREBIN =5;}
		else if (which == "InvMass_j1j2") 	{ xmin=100; xpoint1= 200; xpoint2= 300; xpoint3= 600; xpoint4= 1200; width1=10 ; width2=20 ; width3=50 ; width4=500 ; iREBIN =5;}
		*/
		else if (which == "InvMass_pj1")  	{ xmin=0; xpoint1= 200; xpoint2= 300; xpoint3= 700; xpoint4= 1000; width1=10 ; width2=20 ; width3=50 ; width4=50 ; iREBIN =5;}
		else if (which == "InvMass_pj1j2")	{ xmin=0; xpoint1= 200; xpoint2= 300; xpoint3= 700; xpoint4= 1300; width1=10 ; width2=20 ; width3=50 ; width4=200 ; iREBIN =5;}
		else if (which == "InvMass_j1j2") 	{ xmin=0; xpoint1= 200; xpoint2= 300; xpoint3= 600; xpoint4= 1200; width1=20 ; width2=20 ; width3=50 ; width4=50 ; iREBIN =5;}
		else if (which == "Ht")				 	{ xmin=0  ; xpoint1= 200; xpoint2= 600; xpoint3= 700; xpoint4= 1200; width1=40 ; width2=50 ; width3=100; width4=400 ; iREBIN =5;}
    	else if (which == "Met")			 	{ xmin=0  ; xpoint1= 60 ; xpoint2= 100; xpoint3= 250; xpoint4= 450 ; width1=10 ; width2=20 ; width3=50 ; width4=100 ; iREBIN =5;}
		else if (which == "SecondLeadJetEt"){ xmin=0  ; xpoint1= 100; xpoint2= 160; xpoint3= 200; xpoint4= 600 ; width1=10 ; width2=20 ; width3=40 ; width4=200 ; iREBIN =5;}
    	else if (which == "DelR") 				{ xmin=0  ; xpoint1= 1  ; xpoint2= 2  ; xpoint3= 4  ; xpoint4= 6   ; width1=0.2; width2=0.2; width3=0.2; width4=0.2 ; iREBIN =1;}
    	else if (which == "DelEta") 			{ xmin=0  ; xpoint1= 1  ; xpoint2= 3  ; xpoint3= 4  ; xpoint4= 5   ; width1=0.1; width2=0.1; width3=0.1; width4=0.1 ; iREBIN =1;}
    	else if (which == "DelPhi") 			{ xmin=0  ; xpoint1= 1  ; xpoint2= 2  ; xpoint3= 3  ; xpoint4= 4   ; width1=0.1; width2=0.1; width3=0.1; width4=0.1 ; iREBIN =1;}
	}
	




	/*  combination of settings that would work
	 * 
	 *                                             qcdMcMix      iSample 									iSideband
	 * <1> default method (30%qcd/70%MC PHO)          2        iG30/iG40/iG30MET25/iG30MET40 		iNOMIAL/iISOADDED
	 *             but for njet plot                  1           ""                  						""
	 *
	 *
	 * <2> reweighed method (100% sideband)           1           ""               						iISOADDEDREWGT
	 *
	 */

  
	//int qcdMcMix = 1; int iSample = iG30; int iSideband = iISOADDEDREWGT; sPRINTFOLDER="~/TMP/eps/g30rewgt/";
	int qcdMcMix = 2; int iSample = iG30; int iSideband = iNOMINAL; sPRINTFOLDER="./eps/";
	//int qcdMcMix = 2; int iSample = iG40; int iSideband = iNOMINAL; sPRINTFOLDER="~/TMP/eps/g40nominal/";
	//int qcdMcMix = 1; int iSample = iG40; int iSideband = iISOADDEDREWGT; sPRINTFOLDER="~/TMP/eps/g40weighed/";
	//int qcdMcMix = 2; int iSample = iG30MET25; int iSideband = iNOMINAL; sPRINTFOLDER="~/TMP/eps/g30met25nominal/";
	//int qcdMcMix = 2; int iSample = iG30MET40; int iSideband = iNOMINAL; sPRINTFOLDER="~/TMP/eps/g30met40nominal/";
	//int qcdMcMix = 1; int iSample = iG30MET40; int iSideband = iISOADDEDREWGT; sPRINTFOLDER="~/TMP/eps/g30met40weighed/";
	//int qcdMcMix = 2; int iSample = iG30EXCL1J; int iSideband = iNOMINAL; sPRINTFOLDER="./";
	//int qcdMcMix = 1; int iSample = iG30EXCL1J; int iSideband = iISOADDEDREWGT; sPRINTFOLDER="./";

	if (icase==2)
	{
		qcdMcMix = 1;
		iSideband = iISOADDEDREWGT;
	}
	
	if (which == "NJet") qcdMcMix=1;
	//if (which == "Met") qcdMcMix=1;
	
	if (iSample == iG40) sPRINTFOLDER = "~/ICHEP_LaptopVersion/PhoEt40Files/eps/";
	if (iSample == iG30MET25) sPRINTFOLDER = "~/ICHEP_LaptopVersion/PhoEt30Met25/eps/";
	if (iSample == iG30MET40) sPRINTFOLDER = "~/ICHEP_LaptopVersion/PhoEt30Met40/eps/";


	/**********************************************************
	 * modify the default bin sizes for specific samples here
	 **********************************************************/
	
	if (iSample == iG30MET40)
	{
		sPRINTFOLDER = "~/ICHEP_LaptopVersion/PhoEt30Met40/eps/";
		
		if (jets == 1)
		{
			if (which == "Et_pho")
			{
				if (qcdMcMix == 2) { xmin=0 ; xpoint1=100; xpoint2=200; xpoint3=300; xpoint4=700;  width1=10; width2=20; width3=50 ; width4=200;}
				if (qcdMcMix == 1) { xmin=30 ; xpoint1=100; xpoint2=160; xpoint3=200; xpoint4=600;  width1=10; width2=20; width3=40 ; width4=300;}
			}
		 	else if (which == "Et_j1") { xmin=10; xpoint1=90 ; xpoint2=150; xpoint3=250; xpoint4=550;  width1=10; width2=20; width3=100; width4=200;}
			else if (which == "Ht")    { xmin=0 ; xpoint1=200; xpoint2=400; xpoint3=700; xpoint4=1200; width1=40; width2=40; width3=100; width4=420;}
			else if (which == "InvMass_pj1") 
			{ 
				if (qcdMcMix == 1) { xmin=50 ; xpoint1 = 300 ; xpoint2 = 400; xpoint3 = 550;xpoint4 = 900; width1= 25 ; width2= 50 ; width3= 50 ; width4=270 ;}
			}
				
		}

		
		//2 JETS
		if (jets == 2)
		{
			if (which == "Et_pho")		{ xmin=0; xpoint1=100; xpoint2=160; xpoint3=200; xpoint4=450; width1=10; width2=20; width3=40; width4=220;}
			else if (which == "Et_j1") 		 	{ xmin=14 ; xpoint1= 104; xpoint2= 264; xpoint3= 300; xpoint4= 600 ; width1=10 ; width2=20 ; width3=36 ; width4=200 ;}
			else if (which == "InvMass_pj1")  	{ xmin=100; xpoint1= 200; xpoint2= 300; xpoint3= 500; xpoint4= 800; width1=20 ; width2=20 ; width3=100 ; width4=200 ;}
			else if (which == "InvMass_pj1j2")	{ xmin=100; xpoint1= 200; xpoint2= 400; xpoint3= 600; xpoint4= 1100; width1=20 ; width2=50 ; width3=100 ; width4=400 ;}
			else if (which == "InvMass_j1j2") 	{ xmin=100; xpoint1= 200; xpoint2= 300; xpoint3= 500; xpoint4= 900; width1=20 ; width2=50 ; width3=200 ; width4=350 ;}
			else if (which == "Ht")				 	{ xmin=100; xpoint1= 200; xpoint2= 400; xpoint3= 600; xpoint4= 1200; width1=20 ; width2=50 ; width3=100; width4=520 ;}
		}
	}

	
	if ( bDO_VAR_BINNING && ( (xmin > xpoint1) || (xpoint1 > xpoint2) || (xpoint3 > xpoint4) ))
	{
		cout << __FUNCTION__ << ":" << __LINE__ << ": xpoints must be monotonically increasing!" << std::endl;
		return;
	}
	if ( bDO_VAR_BINNING && (width1 <= 0 || width1 <= 0 || width1 <= 0 || width1 <= 0)) 
	{
		cout << __FUNCTION__ << ":" << __LINE__ << ": bin widths must be >0 !" << std::endl;
		return;
	}
	
	// to make full scale plots
	if (jets == 1) {
		if (which == "Et_pho") 				MakeHistLogAndRatio_debug (which, jets, "EtCorr","E_{T}^{#gamma} (GeV)","1Jet/Photon", 											qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "pj1_Et_pho", 	iSample, iSideband,bDEBUG);
		else if (which == "InvMass_pj1") MakeHistLogAndRatio_debug (which, jets, "InvMass","Invariant Mass (#gamma,Lead Jet) (GeV/c^{2})","1Jet/PhotonLeadJet", 	qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "pj1_InvM_pj1",iSample, iSideband);
		else if (which == "Ht") 			MakeHistLogAndRatio_debug (which, jets, "Ht","H_{T} (GeV)","1Jet/Event",																qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "pj1_Ht", 		iSample, iSideband);
		else if (which == "Et_j1") 		MakeHistLogAndRatio_debug (which, jets, "EtCorr","E_{T}^{Lead Jet} (GeV)","1Jet/LeadJet", 										qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "pj1_Et_j1", 	iSample, iSideband);
    	else if (which == "NJet") 			MakeHistLogAndRatio_debug (which, jets, "NJet15","Jet Multiplicity (E_{T}>15GeV)", "1Jet/Event",								qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
    	else if (which == "Met") 			MakeHistLogAndRatio_debug (which, jets, "Met","MEt", "1Jet/Event",																		qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
    	else if (which == "DelR") 			MakeHistLogAndRatio_debug (which, jets, "DelR","#Delta R (#gamma, lead jet)", "1Jet/PhotonLeadJet",							qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
    	else if (which == "DelEta") 		MakeHistLogAndRatio_debug (which, jets, "DelEta","#Delta #eta (#gamma, lead jet)", "1Jet/PhotonLeadJet",		 				qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
    	else if (which == "DelPhi") 		MakeHistLogAndRatio_debug (which, jets, "DelPhi","#Delta #phi (#gamma, lead jet)", "1Jet/PhotonLeadJet", 					qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
    	else if (which == "DetEta_pho") 	MakeHistLogAndRatio_debug (which, jets, "DetEta","#gamma detector #eta", "1Jet/Photon",                  					qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
    	else if (which == "DetPhi_pho") 	MakeHistLogAndRatio_debug (which, jets, "DetPhi","#gamma detector #phi", "1Jet/Photon",											qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
    	else if (which == "DetEta_j1") 	MakeHistLogAndRatio_debug (which, jets, "DetEta","Lead jet detector #eta", "1Jet/LeadJet",										qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
    	else if (which == "DetPhi_j1") 	MakeHistLogAndRatio_debug (which, jets, "DetPhi","Lead jet detector #phi", "1Jet/LeadJet",										qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
    	else if (which == "EvtPhi_j1") 	MakeHistLogAndRatio_debug (which, jets, "EvtPhi","Lead jet 4-vec #phi", "1Jet/LeadJet",											qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, "", iSample, iSideband);
	}

	if (jets == 2) {
		if (which == "Et_pho") 					MakeHistLogAndRatio_debug (which, jets, "EtCorr","E_{T}^{#gamma} (GeV)","2Jet/Photon", 												qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"pj2_Et_pho", 		iSample, iSideband);
		else if (which == "InvMass_pj1") 	MakeHistLogAndRatio_debug (which, jets, "InvMass","Invariant Mass (#gamma,Lead Jet) (GeV/c^{2})","2Jet/PhotonLeadJet", 		qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"pj2_InvM_pj1", 	iSample, iSideband);
		else if (which == "InvMass_pj1j2") 	MakeHistLogAndRatio_debug (which, jets, "InvMass","Invariant Mass (#gamma, Two Lead Jets) (GeV/c^{2})","2Jet/Photon2Jets", qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"pj2_InvM_pj1j2", 	iSample, iSideband);
		else if (which == "Ht") 				MakeHistLogAndRatio_debug (which, jets, "Ht","H_{T} (GeV)","2Jet/Event", 																	qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"pj1_Ht", 			iSample, iSideband);
		else if (which == "InvMass_j1j2") 	MakeHistLogAndRatio_debug (which, jets, "InvMass","Invariant Mass (Two Lead Jets) (GeV/c^{2})","2Jet/2Jets", 					qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"pj2_InvM_j1j2", 	iSample, iSideband);
		else if (which == "Et_j1") 			MakeHistLogAndRatio_debug (which, jets, "EtCorr","E_{T}^{Lead Jet} (GeV)","1Jet/LeadJet", 											qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"pj2_Et_j1", 		iSample, iSideband);
		else if (which == "SecondLeadJetEt")MakeHistLogAndRatio_debug (which, jets, "EtCorr","E_{T}^{Second Lead Jet} (GeV)","2Jet/SecondLeadJet", 							qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"", iSample, iSideband);
    	else if (which == "Met") 				MakeHistLogAndRatio_debug (which, jets, "Met","MEt", "2Jet/Event",																			qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"", iSample, iSideband);
    	else if (which == "DelR") 				MakeHistLogAndRatio_debug (which, jets, "DelR","#Delta R (#gamma, lead jet)", "2Jet/PhotonLeadJet",								qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"", iSample, iSideband);
    	else if (which == "DelEta") 			MakeHistLogAndRatio_debug (which, jets, "DelEta","#Delta #eta (#gamma, lead jet)", "2Jet/PhotonLeadJet",							qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"", iSample, iSideband);
    	else if (which == "DelPhi") 			MakeHistLogAndRatio_debug (which, jets, "DelPhi","#Delta #phi (#gamma, lead jet)", "2Jet/PhotonLeadJet",							qcdMcMix, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4,"", iSample, iSideband);
	}
	

	if (gPad != pad_old)
	{
		TCanvas *c = dynamic_cast<TCanvas*>(gPad);
		//c->SetFillColor(kRed);
      if (c)
		{
			std::ostringstream str,str1,str2;
			str2 << sPRINTFOLDER << "plot" << jets << "_" << which  << "." << sPRINTFORMAT.c_str();
			c->Print(str2.str().c_str());
		};
	};

};

void MakeHistLogAndRatio_debug (const int i)
{
	if (bDEBUG) std::cout << "Entering " << __FUNCTION__ << ":: " << __LINE__ << std::endl;
	//good plots that works
	if (i == 1) MakeHistLogAndRatio_debug (1, "Et_pho");
	if (i == 2) MakeHistLogAndRatio_debug (1, "Et_j1");
	if (i == 3) MakeHistLogAndRatio_debug (1, "InvMass_pj1");
	if (i == 4) MakeHistLogAndRatio_debug (1, "Ht");
	if (i == 5) MakeHistLogAndRatio_debug (2, "Et_pho");
	if (i == 6) MakeHistLogAndRatio_debug (2, "Et_j1");
	if (i == 7) MakeHistLogAndRatio_debug (2, "InvMass_pj1");
	if (i == 8) MakeHistLogAndRatio_debug (2, "InvMass_pj1j2");
	if (i == 9) MakeHistLogAndRatio_debug (2, "InvMass_j1j2");
	
	if (i == 10) MakeHistLogAndRatio_debug (1, "NJet");

	if (i== 11) MakeHistLogAndRatio_debug (2, "Ht");
	
	if (i== 12) MakeHistLogAndRatio_debug (1, "Met");
	if (i== 13) MakeHistLogAndRatio_debug (1, "DetEta_pho");
	if (i== 14) MakeHistLogAndRatio_debug (1, "DetEta_j1");
	if (i== 15) MakeHistLogAndRatio_debug (1, "DetR_pj1");


	//MakeHistLogAndRatio_debug (1, "Met");
	//MakeHistLogAndRatio_debug (2, "SecondLeadJetEt");
	//MakeHistLogAndRatio_debug (1, "InvMass");
	//MakeHistLogAndRatio_debug (2, "InvMass");
	//MakeHistLogAndRatio_debug (2, "PhoJetsInvMass");
};

void MakeHistLogAndRatio_debug ()
{
	if (bDEBUG) std::cout << "Entering " << __FUNCTION__ << ":: " << __LINE__ << std::endl;
	for (int i=1;i<=15;++i) MakeHistLogAndRatio_debug (i);
}



//list current temporary edit

