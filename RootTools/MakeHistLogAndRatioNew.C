#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <iomanip>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"

//=================================================================//	
// effects: normalize each bin by the bin width
// side effect: normalizes histogram bins in hist
// guarantee: no-throw
// requires: hist != NULL
// requires: hist->GetDimension() == 1
//=================================================================//	
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


//=================================================================//	
// effects: make a new binning with variable bin sizes, and the given
//   cutoff points between bin sizes
// returns: the binning (by bin edges)
// guarantee: strong
// throws: std::bad_alloc
//=================================================================//	
std::vector<float>  GetVarBinVector (float down, float up5, float up10, 
								float up20, float up50, float width1, float width2,
								float width3, float width4, bool debug=false)
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

//=================================================================//	
// effects: turn this histogram into a histogram with variable bin
//   size
// returns: the new histogram
// guarantee: strong
// throws: std::bad_alloc
// requires: input != NULL
// requires: input->GetDimension() == 1
// requires: histogram bin edges must match
// postcondition: result != NULL
//=================================================================//	
std::auto_ptr<TH1> make_var_bin (const std::vector<float>& bins,
	      								TH1 *input, bool debug=false)
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
};





//=================================================================//	
// make a variable binned hist with given specification
// TESTED THIS METHOD AND ITS SUB-METHODS - sam
//=================================================================//	
std::auto_ptr<TH1> MakeVariableBins (TH1 *hist, const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4 , bool use_errors=false, bool debug=false)
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
  	//return result.release ();
  	return result;
};



//=================================================================//	
void MakeHistLogAndRatioNew (std::string which, const int jets, 
					const std::string histname,const std::string title,
				 	const float xmin, const float xpoint1, 
					const float xpoint2, const float xpoint3, 
					const float xpoint4, const float width1, 
					const float width2, const float width3, 
					const float width4)
{

  std::string phofile("PhoJets_data.root");
  std::string mcphofile("PhoJets_phomc.root");
  std::string zeefile("PhoJets_zeemc.root");
  std::string zmmfile("PhoJets_zmmmc.root");
  std::string zttfile("PhoJets_zttmc.root");
  std::string wenfile("PhoJets_wenmc.root");
  std::string wmnfile("PhoJets_wmnmc.root");
  std::string wtnfile("PhoJets_wtnmc.root");
 
 
  TFile* fpho = new TFile(phofile.c_str());
  TFile* fmcpho = new TFile(mcphofile.c_str());
  TFile* fzee = new TFile(zeefile.c_str());
  TFile* fzmm = new TFile(zmmfile.c_str());
  TFile* fztt = new TFile(zttfile.c_str());
  TFile* fwen = new TFile(wenfile.c_str());
  TFile* fwmn = new TFile(wmnfile.c_str());
  TFile* fwtn = new TFile(wtnfile.c_str());


  if (fpho->IsZombie() ||fmcpho->IsZombie() ||
      fzee->IsZombie() ||fzmm->IsZombie() ||
      fztt->IsZombie() ||fwen->IsZombie() ||
      fwmn->IsZombie() ||fwtn->IsZombie() ) 
	{
		std::cout  << "a file not found. pl check" <<std::endl; return;
	}

  // ok path :1Jet/Photon/
  std::string sHaloHist     ("Hist/HALO/"     + histname);
  std::string sCosmicHist   ("Hist/COSMIC/"   + histname);
  std::string sQcdHist      ("Hist/SIDEBAND/" + histname);
  std::string sSigHist      ("Hist/SIGNAL/"   + histname);
  std::string sMcCentralHist("Hist/CENTRAL/"  + histname);
  std::string sMcUpHist     ("Hist/EMJESUP/"  + histname);
  std::string sMcDownHist   ("Hist/EMJESDOWN/"+ histname);


	std::cout << "hist    = "<< histname << std::endl;
	std::cout << "sig.dir = " <<  sSigHist << std::endl;

	
	//all datacased backgrounds are in one file but in different 
	//folders
	TH1F* phojet = dynamic_cast<TH1F*> (fpho->Get(sSigHist.c_str()));
	assert (phojet != NULL && "fpho hist not found in the dir");
			
	TH1F* halojet = dynamic_cast<TH1F*> (fpho->Get(sHaloHist.c_str()));
	assert (halojet != NULL && "halojet hist not found in the dir");

	TH1F* cosmicjet = dynamic_cast<TH1F*> (fpho->Get(sCosmicHist.c_str()));
	assert (cosmicjet != NULL && "cosmicjet hist not found in the dir");

	TH1F* qcdjet = dynamic_cast<TH1F*> (fpho->Get(sQcdHist.c_str()));;
	assert (qcdjet != NULL && "qcdjet hist not found in the dir");

	
	//MC HISTS: jesup= jesup && emup : jesdown = jesdown && emdown
	TH1F* mcphojet = dynamic_cast<TH1F*> (fmcpho->Get(sMcCentralHist.c_str()));
	TH1F* mcphojetJESUP = dynamic_cast<TH1F*> (fmcpho->Get(sMcUpHist.c_str()));
	TH1F* mcphojetJESDOWN = dynamic_cast<TH1F*> (fmcpho->Get(sMcDownHist.c_str()));
	assert ( (mcphojet != NULL && mcphojetJESUP != NULL && mcphojetJESDOWN != NULL) 
			&& "One of the mcpho MC hists (CENTRAL/JES UP/ JES DOWN) not found!");

	TH1F* zeejet = dynamic_cast<TH1F*> (fzee->Get(sMcCentralHist.c_str()));
	TH1F* zeejetJESUP =  dynamic_cast<TH1F*> (fzee->Get(sMcUpHist.c_str()));
	TH1F* zeejetJESDOWN = dynamic_cast<TH1F*> (fzee->Get(sMcDownHist.c_str()));
	assert ( (zeejet != NULL && zeejetJESUP != NULL && zeejetJESDOWN != NULL) 
			&& "One of the zee MC hists (CENTRAL/JES UP/ JES DOWN) not found!");
			
	TH1F* zmmjet = dynamic_cast<TH1F*> (fzmm->Get(sMcCentralHist.c_str()));
	TH1F* zmmjetJESUP = dynamic_cast<TH1F*> (fzmm->Get(sMcUpHist.c_str())); 
	TH1F* zmmjetJESDOWN = dynamic_cast<TH1F*> (fzmm->Get(sMcDownHist.c_str()));
	assert ( (zmmjet != NULL && zmmjetJESUP != NULL && zmmjetJESDOWN != NULL) 
			&& "One of the zmm MC hists (CENTRAL/JES UP/ JES DOWN) not found!");
			
	TH1F* zttjet = dynamic_cast<TH1F*> (fztt->Get(sMcCentralHist.c_str()));
	TH1F* zttjetJESUP = dynamic_cast<TH1F*> (fztt->Get(sMcUpHist.c_str())); 
	TH1F* zttjetJESDOWN = dynamic_cast<TH1F*> (fztt->Get(sMcDownHist.c_str()));
	assert ( (zttjet != NULL && zttjetJESUP != NULL && zttjetJESDOWN != NULL) 
			&& "One of the ztt MC hists (CENTRAL/JES UP/ JES DOWN) not found!");
			
	TH1F* wenjet = dynamic_cast<TH1F*> (fwen->Get(sMcCentralHist.c_str()));
	TH1F* wenjetJESUP = dynamic_cast<TH1F*> (fwen->Get(sMcUpHist.c_str())); 
	TH1F* wenjetJESDOWN = dynamic_cast<TH1F*> (fwen->Get(sMcDownHist.c_str()));
	assert ( (wenjet != NULL && wenjetJESUP != NULL && wenjetJESDOWN != NULL) 
			&& "One of the wmn MC hists (CENTRAL/JES UP/ JES DOWN) not found!");
			
	TH1F* wmnjet = dynamic_cast<TH1F*> (fwmn->Get(sMcCentralHist.c_str()));
	TH1F* wmnjetJESUP = dynamic_cast<TH1F*> (fwmn->Get(sMcUpHist.c_str())); 
	TH1F* wmnjetJESDOWN = dynamic_cast<TH1F*> (fwmn->Get(sMcDownHist.c_str()));
	assert ( (wmnjet != NULL && wmnjetJESUP != NULL && wmnjetJESDOWN != NULL) 
			&& "One of the wmn MC hists (CENTRAL/JES UP/ JES DOWN) not found!");
			
	TH1F* wtnjet = dynamic_cast<TH1F*> (fwtn->Get(sMcCentralHist.c_str()));
	TH1F* wtnjetJESUP = dynamic_cast<TH1F*> (fwtn->Get(sMcUpHist.c_str())); 
	TH1F* wtnjetJESDOWN = dynamic_cast<TH1F*> (fwtn->Get(sMcDownHist.c_str()));
	assert ( (wtnjet != NULL && wtnjetJESUP != NULL && wtnjetJESDOWN != NULL) 
			&& "One of the wtn MC hists (CENTRAL/JES UP/ JES DOWN) not found!");


	//rebin data and central histograms
	std::auto_ptr<TH1> hPhojet = MakeVariableBins (phojet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hHalojet = MakeVariableBins (halojet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hZeejet = MakeVariableBins (zeejet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hZmmjet = MakeVariableBins (zmmjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hZttjet = MakeVariableBins (zttjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hWenjet = MakeVariableBins (wenjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hWmnjet = MakeVariableBins (wmnjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hWtnjet = MakeVariableBins (wtnjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hCosmicjet = MakeVariableBins (cosmicjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hQcdjet = MakeVariableBins (qcdjet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hMcphojet = MakeVariableBins (mcphojet, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
	//rebin MC jes up histograms
  	std::auto_ptr<TH1> hZeejetJESUP = MakeVariableBins (zeejetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hZmmjetJESUP = MakeVariableBins (zmmjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hZttjetJESUP = MakeVariableBins (zttjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hWenjetJESUP = MakeVariableBins (wenjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hWmnjetJESUP = MakeVariableBins (wmnjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hWtnjetJESUP = MakeVariableBins (wtnjetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hMcphojetJESUP = MakeVariableBins (mcphojetJESUP, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
	//rebin MC jes down histograms
  	std::auto_ptr<TH1> hZeejetJESDOWN = MakeVariableBins (zeejetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hZmmjetJESDOWN = MakeVariableBins (zmmjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hZttjetJESDOWN = MakeVariableBins (zttjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hWenjetJESDOWN = MakeVariableBins (wenjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hWmnjetJESDOWN = MakeVariableBins (wmnjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hWtnjetJESDOWN = MakeVariableBins (wtnjetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);
  	std::auto_ptr<TH1> hMcphojetJESDOWN = MakeVariableBins (mcphojetJESDOWN, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4, false);


	//==========================important normalization constants

	const float fDataLum = 2043.0;	//pb-1
	const float kFac = 1.4;				// MC kfactor

	//halo and cosmic background hists should be normalized to these numbers
	float fHaloEst =0, fCosmicEst =0;
	if (jets == 1)
	{
		fHaloEst   = 9;
		fCosmicEst = 110;
	} else 
	{
		fHaloEst   = 1;
		fCosmicEst = 7;
	}

	//fake photon fraction decided how much of QCD and SIGNAL MC is taken
	//fake photon fraction = 0.319+/-0.068(syst) for photon Et>30GeV
	const float fFakePhoFrac       = 0.319;
	const float fFakePhoFrac_sigma = 0.068;

	const double dHalonorm   = fHaloEst / hHalojet->Integral();
	const double dCosmicnorm = fCosmicEst / hCosmicjet->Integral();
	const double dZeenorm = (fDataLum/34056)*kFac;    // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const double dZmmnorm = (fDataLum/38732)*kFac;    // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const double dZttnorm = (fDataLum/27755)*kFac;    // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const double dWennorm = (fDataLum/9438)*kFac;     // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const double dWmnnorm = (fDataLum/5183)*kFac;     // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)
	const double dWtnnorm = (fDataLum/3520)*kFac;     // for EWK mc see. log book#2 pp.72 (old) pp.102 (new)



	//========================== END of important normalization constants

	//======================== BEGIN scaling the hists
	
	hHalojet->Scale (dHalonorm);
	hCosmicjet->Scale (dCosmicnorm);

	hZeejet->Scale (dZeenorm);
	hZmmjet->Scale (dZmmnorm);
	hZttjet->Scale (dZttnorm);
	hWenjet->Scale (dWennorm);
	hWmnjet->Scale (dWmnnorm);
	hWtnjet->Scale (dWtnnorm);

	//QCD and PHO MC samples must be normalized after normalizng all other
	//backgrounds so we can subtract them from data and then split the
	//left over between QCD and PHO MC.
	const double dDataMinusOthers = hPhojet->Integral() - hHalojet->Integral() - hCosmicjet->Integral()
												- hZeejet->Integral() - hZmmjet->Integral() - hZttjet->Integral()
												- hWenjet->Integral() - hWmnjet->Integral() - hWtnjet->Integral();

	hQcdjet->Scale ( (fFakePhoFrac * dDataMinusOthers) / (1.0 * hQcdjet->Integral()));
	hMcphojet->Scale ( ( (1 - fFakePhoFrac) * dDataMinusOthers) / (1.0 * hMcphojet->Integral()));

	const double dSumBackgrouds = hHalojet->Integral() + hCosmicjet->Integral()
											+ hZeejet->Integral() + hZmmjet->Integral() 
											+ hZttjet->Integral() + hWenjet->Integral() 
											+ hWmnjet->Integral() + hWtnjet->Integral()
											+ hQcdjet->Integral() + hMcphojet->Integral();

 	std::cout << "phojet integral = " <<  hPhojet->Integral() << "\t sum backgrounds = "
				 << dSumBackgrouds  << std::endl;















}




//=================================================================//	
void MakeHistLogAndRatioNew (const int jets, const std::string which)
{

	if (jets == 1)
	{
		if (which == "PhotonEt") 	MakeHistLogAndRatioNew (which, jets,"1Jet/Photon/EtCorr" ,"E_{T}^{#gamma} (GeV)", 
																			0,200,250,300,650, 10,10,50,250);
	}

};
//=================================================================//	
void MakeHistLogAndRatioNew()
{
	MakeHistLogAndRatioNew(1, "PhotonEt");
}
