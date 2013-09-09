#include <iostream>
#include <sstream>
#include <string>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <TPaveText.h>
#include <TVectorF.h>
#include <TGraph.h>
#include <TList.h>

using namespace std;

// effects: make a new binning with variable bin sizes, and the given
//   cutoff points between bin sizes
// returns: the binning (by bin edges)
// guarantee: strong
// throws: std::bad_alloc
std::vector<float>
make_var_bin (float down, float up5, float up10, float up20, float up50, bool debug)
{
  std::vector<float> result;
  float point = down;
  const unsigned nregion = 4;
  const float up [] = {up5, up10, up20, up50};
  const float step [] = {10, 10, 50, 250};		//1j pet
  //const float step [] = {10, 20, 50, 150};		//2j pet
  //const float step [] = {25, 50, 75, 300};		//1j inv mass
  //const float step [] = {20, 50, 100, 200};		// 2j inv mass

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
make_var_bin (const std::string& name, const std::vector<float>& bins,
	      TH1 *input, bool debug)
{
  assert (input != NULL && "requirement failed"); //spec
  assert (input->GetDimension() == 1 && "requirement failed"); //spec
  const unsigned nbin = unsigned (input->GetNbinsX ());

  std::auto_ptr<TH1> result (new TH1F (name.c_str(), input->GetTitle(),
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
		for (unsigned i=1; i<=input->GetNbinsX (); i++)
		{
			std::cout << "input bin, lowedge = " << i << ", " << input->GetXaxis()->GetBinLowEdge(i) << "\t" << input->GetBinContent(i) << std::endl;
		}
		std::cout << "\n\n"<< std::endl;
		for (unsigned i=1; i<=result->GetNbinsX (); i++)
		{
			std::cout << "otput bin, lowedge = " << i << ", " << result->GetXaxis()->GetBinLowEdge(i) << "\t" << result->GetBinContent(i) << std::endl;
		}
	}

  assert (result.get() != NULL && "postcondition failed"); //spec
  return result;
};

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

TH1 *get_variableBins (TH1 *hist, int jets, std::string name, float lolimit, float hilimit, bool use_errors, bool debug=false)
{

	std::cout << name << " hist bin size=" << hist->GetBinWidth(1) << std::endl;
  std::auto_ptr<TH1> result;
	if (jets == 1)
	{
		if (name == "EtCorr") result = make_var_bin (name, make_var_bin (0, 200, 250, 300, 650, debug), hist, debug);
		if (name == "InvMass") result = make_var_bin (name, make_var_bin (50, 500, 600, 700, 1000, debug), hist, debug);
	}

	if (jets == 2)
	{
		if (name == "EtCorr") result = make_var_bin (name, make_var_bin (0, 150, 250, 300, 600, debug), hist, debug);
		if (name == "InvMass") result = make_var_bin (name, make_var_bin (100, 600, 700, 800, 1200, debug), hist, debug);
	}
	
	
  norm_bins (result.get());
  return result.release ();
};



float GetMax(const float def, const float up, const float down)
{
	//this is written for JES syst.	
	float m1 = fabs(def - up);
	float m2 = fabs(def - down);
	if (m1>m2) return m1;
	else return m2;
}



void GetJESerrors(std::vector<float>& vErr, TH1* hDef, TH1* hJesUp, TH1* hJesDown)
{
	//call this to get JES syst. make sure all hists have same bin size and normalized properly 
	// and overflow bin is moved in to visible region before calculating systs.
	//for all backgrounds, call this and get the JES error and add in quadrature

	for (int i=1; i< hDef->GetNbinsX(); i++)
		vErr.push_back(GetMax(hDef->GetBinContent(i),
									 hJesUp->GetBinContent(i),
									 hJesDown->GetBinContent(i)));   //________________ dummy value

}


TF1* get_fit_function(TH1* hist,int jets, std::string name, std::string func, int pm=1)
{
	assert (hist != NULL && "hist is a null pointer");

	TH1* test_hist = dynamic_cast<TH1*>(hist->Clone ("fit"));
	for (unsigned bin = 1; bin <= hist->GetNbinsX(); ++ bin)
	{
		test_hist->SetBinContent(bin, pm * hist->GetBinError(bin));
		test_hist->SetBinError(bin, hist->GetBinError(bin)* 0.1);
	}

	test_hist->SetDirectory(0);
	float min=0, max=1000;
	if (name == "EtCorr") 
	{
		min = 30;
		max= 550;
	}
	if (name == "InvMass") 
	{
		min = 50;
		max= 1000;
	}

	std::cout << " function fit range=" << min << "," << max << std::endl;
	test_hist->Fit(func.c_str(),"","", min, max);
	//test_hist->Fit(func.c_str(),"","",30, 450);
	
	TList *list = test_hist->GetListOfFunctions();
	TIterator *it = list->MakeIterator();
	
	float x_min = test_hist->GetBinLowEdge(1);
	float x_max = test_hist->GetXaxis()->GetXmax();
	
//	std::cout << x_min << "\t" << x_max << std::endl;
//	std::cout << test_hist->GetBinLowEdge(1) << "\t" << test_hist->GetXaxis()->GetXmax() << std::endl;

	TF1* tf = new TF1("fit_func",x_min,x_max, 100);
	
	while (tf = (TF1*) it->Next() )
	{
		//tf->Print();
		Double_t p1[5];
		tf->GetParameters(p1);
//		for (int i=0 ;i <5; i++)
//			std::cout << "p=[" << i << "]=" << p1[i] << std::endl;
		if (tf) break;
	}
	
	if (tf) return tf;
	
	std::cout << __FILE__ << "::" << __LINE__ << ":: ERROR! returning NULL pointer!"  << std::endl;
	return NULL;
	
}


TH1 *move_overflow (TH1 *hist, float lolimit, float hilimit, bool use_errors)
{
	gStyle->SetCanvasColor (10);
	gStyle->SetCanvasBorderSize (0);
	gStyle->SetCanvasBorderMode (0);
				 
	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (1);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);


	TH1 *result = new TH1F ((std::string (hist->GetName())+"_tmp").c_str(),hist->GetTitle(),(int)((hilimit - lolimit)/hist->GetBinWidth(1)),lolimit, hilimit);
	result->SetDirectory (NULL);

	if (use_errors) result->Sumw2();

	unsigned target = 1;			//first move all undeflow to first bin
	
	for (unsigned bin = 1; bin <= unsigned (hist->GetNbinsX()); ++ bin)
	{
		if (hist->GetXaxis()->GetBinCenter(bin-1) > lolimit) target++;
		if (hist->GetXaxis()->GetBinCenter(bin) > hilimit) target = result->GetNbinsX();

		
		float val0 = result->GetBinContent (target);
		float val1 = hist->GetBinContent (bin);
		
		result->SetBinContent (target, val0 + val1);
		
		if (use_errors)
		{
			float err0 = result->GetBinError (target);
			float err1 = hist->GetBinError (bin);
			result->SetBinError (target, sqrt (err0*err0 + err1*err1));
		};
	};

	return result;
};

void correct_errors (TH1 *hist, bool use_last = false)
{
	const float fudge = 0.99;
	const unsigned first = 2;
	const unsigned last = hist->GetNbinsX() - 1 - !use_last;
	bool changed = true;

	while (changed)
	{
		changed = false;
		if (hist->GetBinContent (first-1) > 0)
		{
			const float eb = hist->GetBinError (first-1);
			const float ec = hist->GetBinError (first);
			if (eb < fudge * ec)
			{
				hist->SetBinError (first-1, ec);
				changed = true;
			};
		};
		if (hist->GetBinContent (last+1) > 0 || !use_last)
		{
			const float ea = hist->GetBinError (last);
			const float eb = hist->GetBinError (last+1);
			if (eb < fudge * ea)
			{
				hist->SetBinError (last+1, ea);
				changed = true;
			};
		};
		for (unsigned bin = first; bin <= last; ++ bin)
		{
			const float ea = hist->GetBinError (bin-1);
			const float eb = hist->GetBinError (bin);
			const float ec = hist->GetBinError (bin+1);
			if (eb < ea * fudge && eb < ec * fudge)
			{
				hist->SetBinError (bin, ea < ec ? ea : ec);
				changed = true;
			};
		};
	};
};



void GetCosmicErr(std::vector<float>& ErrVec, const TH1 *hist)
{
	for (int i=1; i <= hist->GetNbinsX(); i++) {
		float bin = hist->GetBinContent(i);
		float err = 0;
		if (bin) err = 1 / sqrt(bin);		//take stat error as the syst
		ErrVec.push_back(err);
	}
}

void GetHaloErr(std::vector<float>& ErrVec, const TH1 *hist)
{
	for (int i=1; i <= hist->GetNbinsX(); i++) {
		float bin = hist->GetBinContent(i);
		float err = bin * 0.5;		//take 50% to be the syst
		ErrVec.push_back(err);
	}
}





float GetMax(const float x1, const float x2, const float x3, const float x4)
{
	float max1 = 0, max2 =0;
	
	if (x1>x2) max1 = x1;
	else max1 = x2;
	
	if (x3>x4) max2 = x3;
	else max2 = x4;

	if (max1 > max2) return max1;
	else return max2;
	
}


void GetQCD100Err(std::vector<float>& ErrVec, const TH1* qcdhist, const std::string name,
						const std::string abspath, const int rebin)
{
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
			std::cout << "path not found :" << abspath <<std::endl;
			return;
		}
	
		TH1F* hademhist = (TH1F*) gDirectory->FindObjectAny(name.c_str());
		if (! hademhist){
			std::cout << __LINE__ << " ::hist not found in the dir" <<std::endl;
			return;
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


	//make sure the bin sizes are same as the QCD hist
	//hademhist->Rebin(rebin);
	//isohist->Rebin(rebin);
	//trkpthist->Rebin(rebin);
	//trkisohist->Rebin(rebin);

	float lolimit = qcdhist->GetBinLowEdge(1);
	float hilimit = qcdhist->GetXaxis()->GetBinUpEdge(qcdhist->GetNbinsX());

	int jets = 0;// temporary, need to change this to get from caller!!

   hademhist  = (TH1F*) get_variableBins (hademhist, jets, name ,lolimit, hilimit,true, false);
   isohist  = (TH1F*) get_variableBins (isohist, jets, name ,lolimit, hilimit,true, false);
   trkpthist  = (TH1F*) get_variableBins (trkpthist, jets, name ,lolimit, hilimit,true, false);
   trkisohist  = (TH1F*) get_variableBins (trkisohist, jets, name ,lolimit, hilimit,true, false);
	//hademhist  = (TH1F*) move_overflow (hademhist, lolimit,  hilimit,false);
	//isohist    = (TH1F*) move_overflow (isohist, lolimit,  hilimit,false);
	//trkpthist  = (TH1F*) move_overflow (trkpthist, lolimit,  hilimit,false);
	//trkisohist     = (TH1F*) move_overflow (trkisohist, lolimit,  hilimit,false);

	hademhist->Scale(qcdhist->Integral()/hademhist->Integral());
	isohist->Scale(qcdhist->Integral()/isohist->Integral());
	trkpthist->Scale(qcdhist->Integral()/trkpthist->Integral());
	trkisohist->Scale(qcdhist->Integral()/trkisohist->Integral());

	hademhist->Divide(qcdhist);
	isohist->Divide(qcdhist);
	trkpthist->Divide(qcdhist);
	trkisohist->Divide(qcdhist);

	hademhist->SetMinimum(0);
	hademhist->SetMaximum(10);
	isohist->SetLineColor(kRed);
	trkpthist->SetLineColor(kBlue);
	trkisohist->SetLineColor(kGreen);


	//now find the max for each bin
	for (int i=1; i <= qcdhist->GetNbinsX(); i++) {
		// i do not need to worry about values below 1, as there will be none.
		// each additional cut will give you always equal or lesser number of events.
		// hence all hist values will be above 1.
		float max = GetMax(hademhist->GetBinContent(i),
									isohist->GetBinContent(i),
									trkpthist->GetBinContent(i),
									trkisohist->GetBinContent(i));
		float err = fabs(max -1) * qcdhist->GetBinContent(i);
		//std::cout << "i, max, err=" << i << "\t" << max << "\t" << qcdhist->GetBinContent(i)  << "\t" << err <<std::endl;
		ErrVec.push_back(err);
	}

}




TF1* MakeHistExptBackRatio (const int jets, const std::string& name, const std::string& title,
				 const float lolimit, const float hilimit, const int rebin,
				 const std::string& path, const int qcd_opt, const int QCDerrMtd)
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
		std::cout  << "a file not found. pl check" <<std::endl; exit (1);
	}

	// ok path :1Jet/Photon/
	std::string sHaloDir("Hist/HALO/"+path), sCosmicDir("Hist/COSMIC/"+path), sQcdDir("Hist/SIDEBAND/"+path), sSigDir("Hist/SIGNAL/"+path);
	std::string sMcCentralDir("Hist/CENTRAL/"+path), sMcUpDir("Hist/EMJESUP/"+path), sMcDownDir("Hist/EMJESDOWN/"+path);



	fpho->cd();
	std::cout << "path="<<path<< std::endl;
	std::cout << "dir="<<sSigDir<< std::endl;

	if (! gDirectory->cd(sSigDir.c_str())) {
	 std::cout << "path not found "<< std::endl;
	 exit(1);
	}
	
	TH1F* phojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (! phojet){
	 std::cout << "hist not found in the dir" <<std::endl;
	 exit(1);
	}
			
	gDirectory->pwd();
	fpho->cd();
	gDirectory->cd(sHaloDir.c_str());
	gDirectory->pwd();
	TH1F* halojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	fpho->cd();
	gDirectory->cd(sCosmicDir.c_str());
	gDirectory->pwd();
	TH1F* cosmicjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());

	fpho->cd();
	gDirectory->cd(sQcdDir.c_str());
	gDirectory->pwd();
	TH1F* qcdjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	TH1F* qcdjet_100 = (TH1F*) qcdjet->Clone("qcdjet_100");

	fmcpho->cd();
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* mcphojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fmcpho->cd();
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* mcphojetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fmcpho->cd();
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* mcphojetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			

	fzee->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* zeejet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fzee->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* zeejetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fzee->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* zeejetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fzmm->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* zmmjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fzmm->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* zmmjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fzmm->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* zmmjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fztt->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* zttjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fztt->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* zttjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fztt->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* zttjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fwen->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* wenjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwen->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* wenjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwen->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* wenjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fwmn->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* wmnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwmn->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* wmnjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwmn->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* wmnjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
			
	fwtn->cd();		
	gDirectory->cd(sMcCentralDir.c_str());
	TH1F* wtnjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwtn->cd();		
	gDirectory->cd(sMcUpDir.c_str());
	TH1F* wtnjetJESUP = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	fwtn->cd();		
	gDirectory->cd(sMcDownDir.c_str());
	TH1F* wtnjetJESDOWN = (TH1F*) gDirectory->FindObjectAny(name.c_str());
		
			
  gStyle->SetOptStat(0);


  //REBIN IF NEEDED ------------------------
 phojet    = (TH1F*) get_variableBins (phojet, jets, name ,lolimit, hilimit,true, false);
  halojet    = (TH1F*) get_variableBins (halojet, jets, name ,lolimit, hilimit,true, false);
  zeejet    = (TH1F*) get_variableBins (zeejet, jets, name ,lolimit, hilimit,true, false);
  zmmjet    = (TH1F*) get_variableBins (zmmjet, jets, name ,lolimit, hilimit,true, false);
  zttjet    = (TH1F*) get_variableBins (zttjet, jets, name ,lolimit, hilimit,true, false);
  wenjet    = (TH1F*) get_variableBins (wenjet, jets, name ,lolimit, hilimit,true, false);
  wmnjet    = (TH1F*) get_variableBins (wmnjet, jets, name ,lolimit, hilimit,true, false);
  wtnjet    = (TH1F*) get_variableBins (wtnjet, jets, name ,lolimit, hilimit,true, false);
  cosmicjet    = (TH1F*) get_variableBins (cosmicjet, jets, name ,lolimit, hilimit,true, false);
  qcdjet    = (TH1F*) get_variableBins (qcdjet, jets, name ,lolimit, hilimit,true, false);
  mcphojet->Print();
  mcphojet    = (TH1F*) get_variableBins (mcphojet, jets, name ,lolimit, hilimit,true, false);
  qcdjet_100    = (TH1F*) get_variableBins (qcdjet_100, jets, name ,lolimit, hilimit,true, false);

  zeejetJESUP    = (TH1F*) get_variableBins (zeejetJESUP, jets, name ,lolimit, hilimit,true, false);
  zmmjetJESUP    = (TH1F*) get_variableBins (zmmjetJESUP, jets, name ,lolimit, hilimit,true, false);
  zttjetJESUP    = (TH1F*) get_variableBins (zttjetJESUP, jets, name ,lolimit, hilimit,true, false);
  wenjetJESUP    = (TH1F*) get_variableBins (wenjetJESUP, jets, name ,lolimit, hilimit,true, false);
  wmnjetJESUP    = (TH1F*) get_variableBins (wmnjetJESUP, jets, name ,lolimit, hilimit,true, false);
  wtnjetJESUP    = (TH1F*) get_variableBins (wtnjetJESUP, jets, name ,lolimit, hilimit,true, false);
  mcphojetJESUP    = (TH1F*) get_variableBins (mcphojetJESUP, jets, name ,lolimit, hilimit,true, false);
 
  zeejetJESDOWN    = (TH1F*) get_variableBins (zeejetJESDOWN, jets, name ,lolimit, hilimit,true, false);
  zmmjetJESDOWN    = (TH1F*) get_variableBins (zmmjetJESDOWN, jets, name ,lolimit, hilimit,true, false);
  zttjetJESDOWN    = (TH1F*) get_variableBins (zttjetJESDOWN, jets, name ,lolimit, hilimit,true, false);
  wenjetJESDOWN    = (TH1F*) get_variableBins (wenjetJESDOWN, jets, name ,lolimit, hilimit,true, false);
  wmnjetJESDOWN    = (TH1F*) get_variableBins (wmnjetJESDOWN, jets, name ,lolimit, hilimit,true, false);
  wtnjetJESDOWN    = (TH1F*) get_variableBins (wtnjetJESDOWN, jets, name ,lolimit, hilimit,true, false);
  mcphojetJESDOWN    = (TH1F*) get_variableBins (mcphojetJESDOWN, jets, name ,lolimit, hilimit,true, false);
	

  /***********************************************************************/

  //QCD Scaling options ---------------------
  //fake photon fraction = 0.319+/-0.068(syst)
  float qcd_d = 0.319;
  float qcd_m = 0.251;
  float qcd_p = 0.387;
  float mc_d = 1 - qcd_d;	//.681	//real pho + qcd fake pho = 1
  float mc_p = 1 - qcd_p; //.749
  float mc_m = 1 - qcd_m; //.613

  std::string qcd_str, qcd_title, mc_str;

  float qcd_scale = qcd_d;
  float mc_scale  = mc_d;

  qcd_str = "QCD (#gamma sideband, 0.319(def) of signal)";
  //qcd_str = "QCD (#gamma sideband)";
  mc_str = "MC (0.681(def) of signal)";
  if (qcd_opt == 1) {
    qcd_scale = qcd_p;		//one goes up and other goes down
    mc_scale  = mc_m;
    qcd_str = "QCD (#gamma sideband, 0.387(def+) of signal)";
    mc_str = "#gamma MC (0.613 (def-) of signal)";
  }
  if (qcd_opt == -1) {
    qcd_scale = qcd_m;
    mc_scale  = mc_p;
    qcd_str = "QCD (#gamma sideband, 0.251(def-) of signal)";
    mc_str = "MC (0.749 (def+) of signal)";
  }

  std::ostringstream qcdnum, mcnum;
  qcdnum << "QCD (#gamma sideband, " << qcd_scale<< " of signal)";
  mcnum << "#gamma MC (" << mc_scale << " of signal)";
  qcd_str = qcdnum.str();
  mc_str = mcnum.str();




  // NORMALIZING ---------------------------
  float haloEst 		= 0;
  float cosmicEst 	= 0;

  float haloNorm    = 0;
  float cosmicNorm  = 0;
  float qcdNorm     = 0; 
  float mcphoNorm 	= 0;

  float qcd100Norm     = 1;


  if (jets ==1 ) {
    haloEst = 9;
    cosmicEst = 110;
  }
  if (jets ==2 ) {
    haloEst = 1 ;
    cosmicEst = 7 ;
  }


  haloNorm    = haloEst / halojet->Integral();
  cosmicNorm  = cosmicEst / cosmicjet->Integral();
  qcdNorm     = ((phojet->Integral()) * qcd_scale)/qcdjet->Integral();  //i am looking at only 1/10 of the signal
  qcd100Norm     = ((phojet->Integral()))/qcdjet_100->Integral();  //i am looking at only 1/10 of the signal
  mcphoNorm = ((phojet->Integral()) * mc_scale)/mcphojet->Integral();


  //JES PHO MC
  float mcphoJESUPNorm = ((phojet->Integral()) * mc_scale)/mcphojetJESUP->Integral();
  float mcphoJESDOWNNorm = ((phojet->Integral()) * mc_scale)/mcphojetJESDOWN->Integral();

	//float dataLum=2542.9;	//nb-1
	float dataLum=2043.9;	//nb-1
  float kFac = 1.4;
  //float zeenorm = (dataLum/47605)*kFac;                   // for EWK mc see. log book#2 pp.72
  //float zmmnorm = (dataLum/49577)*kFac;                   // for EWK mc see. log book#2 pp.72
  //float zttnorm = (dataLum/28991)*kFac;                   // for EWK mc see. log book#2 pp.72
  //float wennorm = (dataLum/15408)*kFac;                   // for EWK mc see. log book#2 pp.72
  //float wmnnorm = (dataLum/7704)*kFac;                   // for EWK mc see. log book#2 pp.72
  //float wtnnorm = (dataLum/7704)*kFac;                   // for EWK mc see. log book#2 pp.72

  float zeenorm = (dataLum/34056)*kFac;                   // for EWK mc see. log book#2 pp.72
  float zmmnorm = (dataLum/38732)*kFac;                   // for EWK mc see. log book#2 pp.72
  float zttnorm = (dataLum/27755)*kFac;                   // for EWK mc see. log book#2 pp.72
  float wennorm = (dataLum/9438)*kFac;                   // for EWK mc see. log book#2 pp.72
  float wmnnorm = (dataLum/5183)*kFac;                   // for EWK mc see. log book#2 pp.72
  float wtnnorm = (dataLum/3520)*kFac;                   // for EWK mc see. log book#2 pp.72



  halojet->Scale (haloNorm);
  cosmicjet->Scale (cosmicNorm);
  qcdjet->Scale (qcdNorm);
  qcdjet_100->Scale (qcd100Norm);
  mcphojet->Scale (mcphoNorm);
  zeejet->Scale (zeenorm);
  zmmjet->Scale (zmmnorm);
  zttjet->Scale (zttnorm);
  wenjet->Scale (wennorm);
  wmnjet->Scale (wmnnorm);
  wtnjet->Scale (wtnnorm);

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



	  // subtract the halo? and cosmic from the qcd background. so we are not 
  // double counting them in the MET plot
  // the expections are not very different. see elog#512
  // so did not bother to renormalize them
  qcdjet->Add(cosmicjet, -1);
  qcdjet->Add(halojet, -1);
  qcdjet_100->Add(halojet, -1);
  qcdjet_100->Add(cosmicjet, -1);

					 
	//now subtract EWK background. this is wrong need to fix this	
  qcdjet->Add(zeejet, -1);
  qcdjet->Add(zttjet, -1);
  qcdjet->Add(wenjet, -1);
  qcdjet->Add(wtnjet, -1);
  qcdjet_100->Add(zeejet, -1);
  qcdjet_100->Add(zttjet, -1);
  qcdjet_100->Add(wenjet, -1);
  qcdjet_100->Add(wtnjet, -1);

				 

					 
  int ewkColor = 29;

  mcphojet->SetLineColor(1);
  mcphojet->SetFillColor(5);
  qcdjet->SetLineColor(kRed);
  qcdjet->SetFillColor(6);
  cosmicjet->SetLineColor(1);
  cosmicjet->SetFillColor(kGreen);
  //halojet->SetLineColor(kGreen);
  //halojet->SetFillColor(8);
  halojet->SetLineColor(1);
  halojet->SetFillColor(12);

  zeejet->SetLineColor(10);
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


  phojet->SetLineColor(kBlack);
  phojet->SetMarkerStyle (8);

  //this is the line histo showing 100% of QCD 0% MC pho
  qcdjet_100->SetLineColor(kRed);
  qcdjet_100->SetFillColor(kYellow);

  std::cout << "Ewk: " << (wtnjet->Integral()+zttjet->Integral()+wenjet->Integral()+zeejet->Integral()) << std::endl;



  /************************** GET INFO FOR ERROR BAND ********************/

  std::vector<float> cosmicErr;   //relative error for each bin is stored here
  std::vector<float> haloErr;

  std::vector<float> qcdmcErr;   // use only one of these at a time.

  GetCosmicErr(cosmicErr, cosmicjet);
  GetHaloErr(haloErr, halojet);

  //clones to be used in the syst calculations. need this before normalizing them
  TH1F* qcdjet_syst = (TH1F*) qcdjet->Clone("qcdjet_syst");
  TH1F* mcphojet_syst = (TH1F*) mcphojet->Clone("mcphojet_syst");


  //sort of hist that uses the following two methods
  if (QCDerrMtd == 1) {//plots that use 100% qcd
    std::string histname = name;
    //std::string abspath = newpath.str();
    std::string abspath("Hist/"+path);
    int rebin_temp = rebin;
    GetQCD100Err(qcdmcErr, qcdjet_syst, histname, abspath, rebin_temp);
		//std::cout << "not used .. returning.." << std::endl;
		//exit(1);
  } else if (QCDerrMtd == 2) {											// plots that use the 70%/30% mixture
    //no need for new fucntion
    //just use the info to make new plot with the given mixture
    // and store in the array

    //here i am using some dummy mixture for now. change them wisely!
    float mcphojet_systNorm = ((phojet->Integral()) * 0.40)/mcphojet_syst->Integral();
    mcphojet_syst->Scale(mcphojet_systNorm);
    float qcdjet_systNorm = ((phojet->Integral()) * 0.60)/qcdjet_syst->Integral();
    qcdjet_syst->Scale(qcdjet_systNorm);

    //now create 2 temp hists.
    // one with default qvd/mc values
    //other with varied qcd/mc values them subtract them to get the error
    TH1F *hist_temp1 = (TH1F*) halojet->Clone("hist_temp1");
    TH1F *hist_temp2 = (TH1F*) halojet->Clone("hist_temp2");

    hist_temp1->Add(halojet);
    hist_temp1->Add(wtnjet);
    hist_temp1->Add(zttjet);
    hist_temp1->Add(wenjet);
    hist_temp1->Add(zeejet);
    hist_temp1->Add(cosmicjet);
    hist_temp1->Add(qcdjet_syst);
    hist_temp1->Add(mcphojet_syst);

    hist_temp1->SetLineColor(kRed);
    hist_temp1->SetMarkerColor(kRed);

    hist_temp2->Add(halojet);
    hist_temp2->Add(wtnjet);
    hist_temp2->Add(zttjet);
    hist_temp2->Add(wenjet);
    hist_temp2->Add(zeejet);
    hist_temp2->Add(cosmicjet);
    hist_temp2->Add(qcdjet);
    hist_temp2->Add(mcphojet);
    hist_temp2->SetLineColor(kBlue);
    hist_temp2->SetMarkerColor(kBlue);

    for (int i=1; i <= hist_temp1->GetNbinsX(); i++) {
      float err = fabs(hist_temp1->GetBinContent(i) - hist_temp2->GetBinContent(i));
      //std::cout << "bin i=" << i <<"\t" <<  hist_temp1->GetBinContent(i) <<"\t" << hist_temp2->GetBinContent(i)<< std::endl;
      qcdmcErr.push_back(err);
    }

    //qcdjet_syst->Delete();
    //qcdjet_syst = 0;
  }


  //now get the EWK systematics
  float LumErr = 106/100.; //move up the Lum err by 6%
  //float zeeLumUp = (dataLum/ (47605 * LumErr))*kFac;                   // for EWK mc see. log book#2 pp.72
  //float zmmLumUp = (dataLum/(49577 * LumErr))*kFac;                   // for EWK mc see. log book#2 pp.72
  //float zttLumUp = (dataLum/(28991 * LumErr))*kFac;                   // for EWK mc see. log book#2 pp.72
  //float wenLumUp = (dataLum/(15408 * LumErr))*kFac;                   // for EWK mc see. log book#2 pp.72
  //float wmnLumUp = (dataLum/(7704 * LumErr))*kFac;                   // for EWK mc see. log book#2 pp.72
  //float wtnLumUp = (dataLum/(7704 * LumErr))*kFac;                   // for EWK mc see. log book#2 pp.72

  float zeeLumUp = (dataLum/ (34056 * LumErr))*kFac;  
  float zmmLumUp = (dataLum/(38732 * LumErr))*kFac;  
  float zttLumUp = (dataLum/(27755 * LumErr))*kFac; 
  float wenLumUp = (dataLum/(9438 * LumErr))*kFac; 
  float wmnLumUp = (dataLum/(5183 * LumErr))*kFac;
  float wtnLumUp = (dataLum/(3520 * LumErr))*kFac;


  float zeeErr = fabs(zeeLumUp - zeenorm);
  float zmmErr = fabs(zmmLumUp - zmmnorm);
  float zttErr = fabs(zttLumUp - zttnorm);
  float wenErr = fabs(wenLumUp - wennorm);
  float wmnErr = fabs(wmnLumUp - wmnnorm);
  float wtnErr = fabs(wtnLumUp - wtnnorm);


	// NOW GET JES SYSTEMATICS

	std::vector<float> zeeJESerr, zmmJESerr, zttJESerr;
	std::vector<float> wenJESerr, wmnJESerr, wtnJESerr ,mcphoJESerr;
	
	GetJESerrors(zeeJESerr, zeejet, zeejetJESUP, zeejetJESDOWN);
	GetJESerrors(zmmJESerr, zmmjet, zmmjetJESUP, zmmjetJESDOWN);
	GetJESerrors(zttJESerr, zttjet, zttjetJESUP, zttjetJESDOWN);
	GetJESerrors(wenJESerr, wenjet, wenjetJESUP, wenjetJESDOWN);
	GetJESerrors(wmnJESerr, wmnjet, wmnjetJESUP, wmnjetJESDOWN);
	GetJESerrors(wtnJESerr, wtnjet, wtnjetJESUP, wtnjetJESDOWN);
	GetJESerrors(mcphoJESerr, mcphojet, mcphojetJESUP, mcphojetJESDOWN);


  TH1F *hist_err = (TH1F*) halojet->Clone("hist_err");
  if (QCDerrMtd == 1) {//plots that use 100% qcd
    hist_err->Add(zmmjet);
    hist_err->Add(wtnjet);
    hist_err->Add(zttjet);
    hist_err->Add(wmnjet);
    hist_err->Add(wenjet);
    hist_err->Add(zeejet);
    hist_err->Add(cosmicjet);
    hist_err->Add(qcdjet_100);
  } else if (QCDerrMtd == 2) {//plots that use 30%70%
    hist_err->Add(zmmjet);
    hist_err->Add(wtnjet);
    hist_err->Add(zttjet);
    hist_err->Add(wmnjet);
    hist_err->Add(wenjet);
    hist_err->Add(zeejet);
    hist_err->Add(cosmicjet);
    hist_err->Add(qcdjet);
    hist_err->Add(mcphojet);
  } else
    {
      std::cerr << "bad value: QCDerrMtd=" << QCDerrMtd << std::endl;
      exit(1);
    };


	// NOW GET THE STAT ERROR
	std::vector<float> statErr;
	for (int i=1; i <= hist_err->GetNbinsX(); i++)
	{
		float err = hist_err->GetBinError(i);
		//if (err == 0 && i > 5)
		//{
		//	int lastbin= hist_err->GetNbinsX();
		//	if (i != lastbin)
		//	{
		//		float errup = hist_err->GetBinError(i+1);
		//		float errdown = hist_err->GetBinError(i-1);
		//		err = errup+errdown/2;
		//	}
		//}
		statErr.push_back(err);
	}




  // NOW COLLECT ALL ERRORS AND PUT THEM TOGETHER FOR ONE FINAL NUMBER
  std::vector<float> ERRORS;

  
  TH1F *hist_sys_jes = (TH1F*) halojet->Clone("hist_sys_jes");
  TH1F *hist_sys_cosmic = (TH1F*) halojet->Clone("hist_sys_cosmic");
  TH1F *hist_sys_qcdmc = (TH1F*) halojet->Clone("hist_sys_qcdmc");
  TH1F *hist_sys_halo = (TH1F*) halojet->Clone("hist_sys_halo");
  TH1F *hist_sys_ewk = (TH1F*) halojet->Clone("hist_sys_ewk");
  TH1F *hist_sys_stat = (TH1F*) halojet->Clone("hist_sys_stat");
  

  for (unsigned int i=0; i < qcdmcErr.size(); i++) {
    float sum = pow(cosmicErr[i],2) + pow(haloErr[i],2) + pow(zeeErr,2) + pow(zmmErr,2)
      + pow(zttErr,2) + pow(wenErr,2) + pow(wmnErr,2) + pow(wtnErr,2) + pow(qcdmcErr[i],2)
		+ pow(zeeJESerr[i],2) + pow(zttJESerr[i],2) + pow(wenJESerr[i],2) 
		+ pow(wtnJESerr[i],2) 
		+ pow(mcphoJESerr[i],2)
		+ pow(statErr[i],2)
		+ pow(wmnJESerr[i],2) 
		+pow(zmmJESerr[i],2); 
		
		 float err = sqrt(sum);
		 ERRORS.push_back(err);

		 float sumJes = sqrt( pow(zeeJESerr[i],2)+pow(zmmJESerr[i],2)+ pow(zttJESerr[i],2) + pow(wenJESerr[i],2) 
									+ pow(wmnJESerr[i],2) + pow(wtnJESerr[i],2) + pow(mcphoJESerr[i],2) );
    	 float sumEwk = sqrt ( pow(zeeErr,2) + pow(zmmErr,2) + pow(zttErr,2) 
									+ pow(wenErr,2) + pow(wmnErr,2) + pow(wtnErr,2));

		hist_sys_jes->SetBinContent(i+1, sumJes);
		hist_sys_cosmic->SetBinContent(i+1, cosmicErr.at(i));
		hist_sys_qcdmc->SetBinContent(i+1, qcdmcErr.at(i));
		hist_sys_halo->SetBinContent(i+1, haloErr.at(i));
		hist_sys_ewk->SetBinContent(i+1, sumEwk);
		hist_sys_stat->SetBinContent(i+1, statErr.at(i));
	 
  }
 correct_errors(hist_err);

 new TCanvas();
 halojet->Draw();
 for (unsigned i=1; i <= halojet->GetNbinsX(); i++)
 {
	std::cout << "bin " << i << " = " << halojet->GetXaxis()->GetBinLowEdge(i) << std::endl;
 }

	hist_sys_jes->SetMarkerColor(kRed);
	hist_sys_cosmic->SetMarkerColor(kGreen);
	hist_sys_qcdmc->SetMarkerColor(kBlue);
	hist_sys_halo->SetMarkerColor(5);
	hist_sys_ewk->SetMarkerColor(28);
	hist_sys_stat->SetMarkerColor(kBlack);
	hist_sys_jes->SetLineColor(kRed);
	hist_sys_cosmic->SetLineColor(kGreen);
	hist_sys_qcdmc->SetLineColor(kBlue);
	hist_sys_halo->SetLineColor(kYellow);
	hist_sys_ewk->SetLineColor(28);
	hist_sys_stat->SetLineColor(kBlack);

	

	new TCanvas();
	gPad->SetLogy();
	hist_sys_jes->Draw("P");
	hist_sys_qcdmc->Draw("sameP");
	hist_sys_cosmic->Draw("sameP");
	hist_sys_halo->Draw("sameP");
	hist_sys_ewk->Draw("sameP");
	hist_sys_stat->Draw("sameP");
	new TCanvas();
	
	TLegend *leg2= new TLegend (0.5,0.7,0.9,0.9);
	leg2->SetTextFont(42);
	leg2->SetTextSize(0.03);
	leg2->SetBorderSize (1);
	leg2->SetFillColor (10);
	leg2->AddEntry(hist_sys_jes,"JES");
	leg2->AddEntry(hist_sys_qcdmc,"QCD/MC Mix");
	leg2->AddEntry(hist_sys_cosmic,"Cosmic");
	leg2->AddEntry(hist_sys_halo,"Halo");
	leg2->AddEntry(hist_sys_ewk,"EWK Lum");
	leg2->AddEntry(hist_sys_stat,"Stat");
	leg2->Draw();


  hist_err->SetFillStyle(3002);
  hist_err->SetFillColor(kRed);


  /************************* end error bars stuff **********************/



	THStack *hs = new THStack ("hs", NULL);
	//cosmicjet->SetXTitle (title.c_str());
	//halojet->SetXTitle (title.c_str());

	std::ostringstream ytitle;
	ytitle << "(Data -- Background) / Background";
//  hs->GetYaxis()->SetTitle(ytitle.str().c_str());

	std::string temp_title;
	if (jets == 1)
	{
		if (name == "EtCorr") temp_title += "E_{T}^{ #gamma} (GeV)";
		if (name == "InvMass") temp_title += "M^{ #gamma, Lead Jet} (GeV/c^{2})";
		temp_title += " (#gamma+>=1 Jet)";
	}
	if (jets == 2) 
	{
		if (name == "EtCorr") temp_title += "E_{T}^{ #gamma} (GeV)";
		if (name == "InvMass") temp_title += "M^{ #gamma, Lead Jet} (GeV/c^{2})";
		temp_title += " (#gamma+>=2 Jets)";
	}

	for (unsigned bin = 0; bin <= hist_err->GetNbinsX() + 1; ++ bin)
	{
		const float val = hist_err->GetBinContent (bin);
		const float scale = val ? 1. / val : 0;
		phojet->SetBinContent (bin, (phojet->GetBinContent (bin) - val) * scale);
		phojet->SetBinError (bin, phojet->GetBinError (bin) * scale);
	};
	
	phojet->SetTitle("");
	phojet->GetXaxis()->SetTitle(temp_title.c_str());
	std::cout << "title= " << title << std::endl;
	phojet->GetYaxis()->CenterTitle(true);
	phojet->GetYaxis()->SetTitle(ytitle.str().c_str());
	phojet->Print();
	phojet->SetMinimum (-1.5);
	phojet->SetMaximum (1.5);
	
	TH1 *hist_err_copy = NULL;
	{
		std::string myname = hist_err->GetName() + std::string ("_copy");
		hist_err_copy = dynamic_cast<TH1*>(hist_err->Clone (myname.c_str()));
		for (unsigned bin = 0; bin <= hist_err_copy->GetNbinsX() + 1; ++ bin)
		{
			float value = hist_err_copy->GetBinContent (bin);
			float error = hist_err_copy->GetBinError (bin);
			hist_err_copy->SetBinError (bin, value ? error / value : 0);
			hist_err_copy->SetBinContent (bin, 0);
		};
	};
	
	std::string func("pol2");
	if (jets == 1 && name == "EtCorr") func = "pol3";
	if (jets == 2 && name == "EtCorr") func = "pol3";
	
	TH1* tt = dynamic_cast<TH1*>(hist_err_copy->Clone("temp"));
		tt->SetDirectory(0);
	TF1 *tf = get_fit_function(tt, jets, name, func, 1);
	TF1 *tf2 = get_fit_function(tt, jets, name, func, -1);

	
	new TCanvas();
	gStyle->SetCanvasColor (10);
	gStyle->SetCanvasBorderSize (0);
	gStyle->SetCanvasBorderMode (0);
				 
	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (0);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);

	gPad->UseCurrentStyle();

  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetGridx();
  gPad->SetGridy();
  

 int labelfont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
 int titlefont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
 phojet->GetXaxis()->SetLabelFont(labelfont);
 phojet->GetYaxis()->SetLabelFont(labelfont);
 phojet->GetXaxis()->SetTitleFont(titlefont);
 phojet->GetYaxis()->SetTitleFont(titlefont);
 phojet->GetXaxis()->CenterTitle(true);
 phojet->GetYaxis()->CenterTitle(true);



  phojet->SetMarkerColor(kBlue);
  phojet->SetLineColor(kBlue);
 // phojet->GetXaxis()->SetTitle(title.c_str());
 // phojet->GetXaxis()->CenterTitle(true);
  phojet->Draw();
  hist_err_copy->Draw ("SAME E2");


	if (!tf) std::cout << "errr" << std::endl;
	else 	{
		tf->Draw("sameC");
		tf2->Draw("sameC");
	}
	
  	TPaveText *tp = new TPaveText(0.5,0.91,0.9,0.99,"NDC");
	tp->SetLineColor(10);
	tp->SetTextFont(titlefont);
	tp->AddText("CDF Run II Preliminary 2.0 fb^{-1}");
	tp->Draw();

  TLegend *leg = new TLegend (0.1,0.72,0.5,0.9);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetBorderSize (1);
  leg->SetFillColor (10);
  leg->AddEntry(phojet,"Data");
  leg->AddEntry(hist_err_copy,"Systematic Uncertainty");
  leg->AddEntry(tf,"Polynomial fit for systematic band");
  leg->Draw();

		TCanvas *c = dynamic_cast<TCanvas*>(gPad);
      if (c)
		{
	  		std::ostringstream str,str1;
			str << "ratioplot" << jets << "_" << name << ".gif";
			c->Print (str.str().c_str());
			str1 << "ratioplot" << jets << "_" << name << ".pdf";
			c->Print (str1.str().c_str(),"pdf");
		};



	return tf;

  
/*
	if ((jets==1) && (name=="InvMass")) {
		float val0 = 0.2;
		float val1 = 1.5;
		float width = 700;
		float power = 3.0;
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_1jInvMass1 = new TF1("eq_1jInvMass1",func.str().c_str(),0,700);
		TF1 *eq_1jInvMass2 = new TF1("eq_1jInvMass2",("-("+func.str()+")").c_str(),0,700);				
		//TF1 *eq_1jInvMass1 = new TF1("eq_1jInvMass1","0.1846-0.0020*x+1.0000052*x*x",0,700);
		//TF1 *eq_1jInvMass2 = new TF1("eq_1jInvMass2","-0.1846+0.0012*x-0.0000052*x*x",0,700);
		eq_1jInvMass1->Draw("sameC");
		eq_1jInvMass2->Draw("sameC");
	}
	bool pjinvmass = true;
	if (path.find("Photon2Jets") == std::string::npos) pjinvmass = false;
	
	if ((jets==2) && (name=="InvMass") && pjinvmass) {		//pho+2jets
		//float val0 = 0.05;
		//float val1 = 1.05;
		//float width = 800;
		//float power = 3.1;
		
		float val0 = 0.2;
		float val1 = 1.2;
		float width = 800;
		float power = 3.1;
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_2jInvMass1 = new TF1("eq_2jInivMass1",func.str().c_str(),0,width);
		TF1 *eq_2jInvMass2 = new TF1("eq_2jInvMass2",("-("+func.str()+")").c_str(),0,width);				
		eq_2jInvMass1->Draw("sameC");
		eq_2jInvMass2->Draw("sameC");
	}

	if ((jets==2) && (name=="InvMass") && !pjinvmass) {  // jets inv mass
		float val0 = 0.05;
		float val1 = 1.05;
		float width = 800;
		float power = 3.1;
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_2jInvMass1 = new TF1("eq_2jInivMass1",func.str().c_str(),0,800);
		TF1 *eq_2jInvMass2 = new TF1("eq_2jInvMass2",("-("+func.str()+")").c_str(),0,800);				
		//TF1 *eq_1jInvMass1 = new TF1("eq_1jInvMass1","0.1846-0.0020*x+1.0000052*x*x",0,700);
		//TF1 *eq_1jInvMass2 = new TF1("eq_1jInvMass2","-0.1846+0.0012*x-0.0000052*x*x",0,700);
		eq_2jInvMass1->Draw("sameC");
		eq_2jInvMass2->Draw("sameC");
	}

	bool phoet = true;
	if (path.find("Photon") == std::string::npos) phoet = false;
	if ((jets==1) && (name=="EtCorr" && phoet)) {
		//float val0 = 0.1;
		//float val1 = 1.1;
		//float width = 300;
		//float power = 2.0;
		
		float val0 = 0.15;
		float val1 = 1.5;
		float width = 300;
		float power = 2.5;
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_1jpet1 = new TF1("eq_1jpeta",func.str().c_str(),0,300);
		TF1 *eq_1jpet2 = new TF1("eq_1jpetb",("-("+func.str()+")").c_str(),0,300);				
		eq_1jpet1->Draw("sameC");
		eq_1jpet2->Draw("sameC");
	}

	if ((jets==2) && (name=="EtCorr") && phoet) {
		//float val0 = 0.11;
		//float val1 = 1.0;
		//float width = 300;
		//float power = 2.0;
		
		float val0 = 0.14;
		float val1 = 1.5;
		float width = 300;
		float power = 2.2;
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_2jpet1 = new TF1("eq_2jpeta",func.str().c_str(),0,300);
		TF1 *eq_2jpet2 = new TF1("eq_2jpetb",("-("+func.str()+")").c_str(),0,300);				
		eq_2jpet1->Draw("sameC");
		eq_2jpet2->Draw("sameC");
	}

	if ((jets==1) && name=="EtCorr" && !phoet) {		//lead jet et
		float val0 = 0.1;
		float val1 = 1.0;
		float width = 300;
		float power = 2;
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_1jpet1 = new TF1("eq_1jpeta",func.str().c_str(),0,300);
		TF1 *eq_1jpet2 = new TF1("eq_1jpetb",("-("+func.str()+")").c_str(),0,300);				
		eq_1jpet1->Draw("sameC");
		eq_1jpet2->Draw("sameC");
	}

	if ((jets==1) && (name=="NJet15")) {
		//float val0 = 0.05;
		//float val1 = 1.25;
		//float width = 15;
		//float power = 2;
		
		float val0 = 0.1;
		float val1 = 5;
		float width = 10;
		float power = 2.3;
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_1jpet1 = new TF1("eq_1jnjeta",func.str().c_str(),0,15);
		TF1 *eq_1jpet2 = new TF1("eq_1jnjetb",("-("+func.str()+")").c_str(),0,15);				
		eq_1jpet1->Draw("sameC");
		eq_1jpet2->Draw("sameC");
	}
*/
/*
	if ((jets==2) && (name=="NJet15")) {
		float val0 = 0.1;
		float val1 = 0.1;
		float width = 15;
		float power = 2;
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_2jpet1 = new TF1("eq_2jnjeta",func.str().c_str(),0,15);
		TF1 *eq_2jpet2 = new TF1("eq_2jnjetb",("-("+func.str()+")").c_str(),0,15);				
		eq_2jpet1->Draw("sameC");
		eq_2jpet2->Draw("sameC");
	}

	if ((jets==1) && (name=="Ht")) {		//Ht
		float val0 = 0.05;
		float val1 = 0.95;
		float width = 700;
		float power = 2;
		
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_1jpet1 = new TF1("eq_1jpeta",func.str().c_str(),0,700);
		TF1 *eq_1jpet2 = new TF1("eq_1jpetb",("-("+func.str()+")").c_str(),0,700);				
		eq_1jpet1->Draw("sameC");
		eq_1jpet2->Draw("sameC");
	}
	if ((jets==2) && (name=="Ht")) {		//Ht
		float val0 = 0.05;
		float val1 = 0.95;
		float width = 700;
		float power = 2;
		
		std::ostringstream func;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		TF1 *eq_1jpet1 = new TF1("eq_1jpeta",func.str().c_str(),0,700);
		TF1 *eq_1jpet2 = new TF1("eq_1jpetb",("-("+func.str()+")").c_str(),0,700);				
		eq_1jpet1->Draw("sameC");
		eq_1jpet2->Draw("sameC");
	}
*/

	
};

void MakeHistExptBackRatio (int jets,std::string name,  std::string opt = "")
{
  const bool opt_qcdup = (opt.find ("QCDUP") != std::string::npos);
  const bool opt_qcddown = (opt.find ("QCDDOWN") != std::string::npos);


  int qcd_scale = 0;
  TVirtualPad *pad_old = gPad;	

	
  if (opt_qcdup) qcd_scale = 1;
  if (opt_qcddown) qcd_scale = -1;
/*
  if (jets == 1) {
    if (name == "InvMass") 	MakeHistExptBackRatio (jets, "InvMass","Invariant Mass(#gamma,Lead Jet)",0,700,5,"1Jet/PhotonLeadJet", qcd_scale, 2);
    else if (name == "pet") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{#gamma}",0,300,5,"1Jet/Photon", qcd_scale,2);
    else if (name == "met")	MakeHistExptBackRatio (jets, "Met","MET",0,200,1,"1Jet/Event", qcd_scale,1);
    else if (name == "njet") 	MakeHistExptBackRatio (jets, "NJet15","NJet15",0,15,1,"1Jet/Event", qcd_scale,2);
    else if (name == "ht")  	MakeHistExptBackRatio (jets, "Ht","H_{T}",0,700,4,"1Jet/Event", qcd_scale,2);
    else if (name == "jet") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{lead Jet}",0,300,5,"1Jet/LeadJet", qcd_scale,2);
    else if (name == "jetEta") 	MakeHistExptBackRatio (jets, "DetEta","Detector #eta ^{lead Jet}",-5,5,4,"1Jet/LeadJet", qcd_scale,2);
    else if (name == "etratio") MakeHistExptBackRatio (jets, "EtRatio","E_{T} ratio (Lead Jet/ #gamma)",0,10,2,"1Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "delphi") 	MakeHistExptBackRatio (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,5,1,"1Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "delr") MakeHistExptBackRatio (jets, "DelR", "#Delta R^{#gamma, lead Jet}",0,7,1,"1Jet/PhotonLeadJet", qcd_scale,2);
  }
  if (jets == 2) {
    if (name == "InvMass") 	MakeHistExptBackRatio (jets, "InvMass","Invariant Mass(#gamma,Two Lead Jets)",0,800,5,"2Jet/Photon2Jets", qcd_scale,2);
    else if (name == "Met")	MakeHistExptBackRatio (jets, "Met","MET",0,125,1,"2Jet/Event", qcd_scale,1);
    if (name == "jetsInvMass") 	MakeHistExptBackRatio (jets, "InvMass","Invariant Mass(Lead two Lead Jets)",0,700,5,"2Jet/2Jets", qcd_scale,2);
    else if (name == "NJet15") 	MakeHistExptBackRatio (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale,1);
    else if (name == "PhoEt") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{#gamma}",0,300,5,"2Jet/Photon", qcd_scale,2);
    else if (name == "pj1DelPhi") 	MakeHistExptBackRatio (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,4,1,"2Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "pj2DelPhi") 	MakeHistExptBackRatio (jets, "DelPhi","#Delta #phi^{#gamma,2nd lead jet}",0,4,1,"2Jet/Photon2ndLeadJet", qcd_scale,2);
    else if (name == "j1j2DelPhi") 	MakeHistExptBackRatio (jets, "DelPhi","#Delta #phi^{lead jet,2nd lead jet}",0,4,1,"2Jet/2Jets", qcd_scale,2);
    else if (name == "Ht")  	MakeHistExptBackRatio (jets, "Ht","H_{T}",0,700,4,"2Jet/Event", qcd_scale,2);
    else if (name == "j1Et") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{lead Jet}",0,300,5,"2Jet/LeadJet", qcd_scale,2);
    else if (name == "j2Et") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{2nd lead Jet}",0,300,5,"2Jet/SecondLeadJet", qcd_scale,2);
    else if (name == "pj1Etratio") MakeHistExptBackRatio (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,1,"2Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "pj2Etratio") MakeHistExptBackRatio (jets, "EtRatio","E_{T} ratio (2nd lead Jet/ #gamma)",0,10,1,"2Jet/Photon2ndLeadJet", qcd_scale,2);
    else if (name == "j1j2Etratio") MakeHistExptBackRatio (jets, "EtRatio","E_{T} ratio (2nd lead Jet/lead Jet)",0,2.5,3,"2Jet/2Jets", qcd_scale,2);
    else if (name == "j1jetEta") 	MakeHistExptBackRatio(jets, "DetEta","Detector #eta ^{Lead Jet}",-5,5,4,"2Jet/LeadJet", qcd_scale,2);
    else if (name == "j2jetEta") 	MakeHistExptBackRatio (jets, "DetEta","Detector #eta ^{Second lead Jet}",-5,5,4,"2Jet/SecondLeadJet", qcd_scale,2);
  }
*/	
	// to make full scale plots
  if (jets == 1) {
    if (name == "InvMass") 	MakeHistExptBackRatio (jets, "InvMass","Invariant Mass (#gamma,Lead Jet) (GeV/c^{2})",0,1100,5,"1Jet/PhotonLeadJet", qcd_scale, 2);
    else if (name == "pet") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{#gamma} (GeV)",0,600,5,"1Jet/Photon", qcd_scale,2);
    //else if (name == "met")	MakeHistExptBackRatio (jets, "Met","MET",0,1000,1,"1Jet/Event", qcd_scale,1);
    else if (name == "met")	MakeHistExptBackRatio (jets, "Met","MET",0,1000,1,"1Jet/Event", qcd_scale,2);
    //else if (name == "njet") 	MakeHistExptBackRatio (jets, "NJet15","Jet Multiplicity (E_{T}>15GeV)",0,15,1,"1Jet/Event", qcd_scale,1);
    else if (name == "njet") 	MakeHistExptBackRatio (jets, "NJet15","Jet Multiplicity (E_{T}>15GeV)",0,15,1,"1Jet/Event", qcd_scale,1);
    else if (name == "ht")  	MakeHistExptBackRatio (jets, "Ht","H_{T} (GeV)",0,1500,4,"1Jet/Event", qcd_scale,1);
    else if (name == "jet") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{lead Jet}",0,600,5,"1Jet/LeadJet", qcd_scale,2);
    else if (name == "etratio") MakeHistExptBackRatio (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,2,"1Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "delphi") 	MakeHistExptBackRatio (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,5,1,"1Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "delr") MakeHistExptBackRatio (jets, "DelR", "#Delta R^{#gamma, lead Jet}",0,7,1,"1Jet/PhotonLeadJet", qcd_scale,2);
  }
  if (jets == 2) {
    if (name == "InvMass") 	MakeHistExptBackRatio (jets, "InvMass","Invariant Mass (#gamma,Two Lead Jets) (GeV/c^{2})",0,1100,5,"2Jet/Photon2Jets", qcd_scale,2);
    //else if (name == "Met")	MakeHistExptBackRatio (jets, "Met","MET",0,1000,1,"2Jet/Event", qcd_scale,1);
    else if (name == "Met")	MakeHistExptBackRatio (jets, "Met","MET",0,1000,1,"2Jet/Event", qcd_scale,2);
    else if (name == "jetsInvMass") 	MakeHistExptBackRatio (jets, "InvMass","Invariant Mass(Lead two Lead Jets)",0,1100,5,"2Jet/2Jets", qcd_scale,2);
    //else if (name == "NJet15") 	MakeHistExptBackRatio (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale,1);
    //else if (name == "NJet15") 	MakeHistExptBackRatio (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale,2);
    else if (name == "PhoEt") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{#gamma} (GeV)",0,500,5,"2Jet/Photon", qcd_scale,2);
    else if (name == "pj1DelPhi") 	MakeHistExptBackRatio (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,4,1,"2Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "pj2DelPhi") 	MakeHistExptBackRatio (jets, "DelPhi","#Delta #phi^{#gamma,2nd lead jet}",0,4,1,"2Jet/Photon2ndLeadJet", qcd_scale,2);
    else if (name == "j1j2DelPhi") 	MakeHistExptBackRatio (jets, "DelPhi","#Delta #phi^{lead jet,2nd lead jet}",0,4,1,"2Jet/2Jets", qcd_scale,2);
    else if (name == "Ht")  	MakeHistExptBackRatio (jets, "Ht","H_{T} (GeV)",0,1500,4,"2Jet/Event", qcd_scale,2);
    else if (name == "j1Et") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{lead Jet}",0,1000,5,"2Jet/LeadJet", qcd_scale,2);
    else if (name == "j2Et") 	MakeHistExptBackRatio (jets, "EtCorr","E_{T}^{2nd lead Jet}",0,1000,5,"2Jet/SecondLeadJet", qcd_scale,2);
    else if (name == "pj1Etratio") MakeHistExptBackRatio (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,1,"2Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "pj2Etratio") MakeHistExptBackRatio (jets, "EtRatio","E_{T} ratio (2nd lead Jet/ #gamma)",0,10,1,"2Jet/Photon2ndLeadJet", qcd_scale,2);
    else if (name == "j1j2Etratio") MakeHistExptBackRatio (jets, "EtRatio","E_{T} ratio (2nd lead Jet/lead Jet)",0,2.5,3,"2Jet/2Jets", qcd_scale,2);
  }
	



	if (gPad != pad_old)
	{
		TCanvas *c = dynamic_cast<TCanvas*>(gPad);
      if (c)
		{
	  		std::ostringstream str,str1;
			str << "ratioplot" << jets << "_" << name << ".gif";
			c->Print (str.str().c_str());
			//str1 << "ratioplot" << jets << "_" << name << ".pdf";
			//c->Print (str1.str().c_str(),"pdf");
		};
	};

	
};




void MakeHistExptBackRatio ()
{
	//MakeHistExptBackRatio (1, "InvMass");
	MakeHistExptBackRatio (1, "pet");
	//MakeHistExptBackRatio (1, "met");
	//MakeHistExptBackRatio (1, "njet");
	//MakeHistExptBackRatio (1, "ht");
	//MakeHistExptBackRatio (1, "jet");
	//MakeHistExptBackRatio (1, "jetEta");
	//MakeHistExptBackRatio (1, "etratio");
	//MakeHistExptBackRatio (1, "delphi");
	//MakeHistExptBackRatio (1, "delr");
	//MakeHistExptBackRatio (2, "InvMass");
	//MakeHistExptBackRatio (2, "Met");
	//MakeHistExptBackRatio (2, "jetsInvMass");		//need to fix name. this uses the function for pho+jets inav mass to get errors
	//MakeHistExptBackRatio (2, "NJet15");
	//MakeHistExptBackRatio (2, "PhoEt");
	//MakeHistExptBackRatio (2, "pj1DelPhi");
	//MakeHistExptBackRatio (2, "pj2DelPhi");
	//MakeHistExptBackRatio (2, "j1j2DelPhi");
	//MakeHistExptBackRatio (2, "Ht");
	//MakeHistExptBackRatio (2, "j1Et");
	//MakeHistExptBackRatio (2, "j2Et");
	//MakeHistExptBackRatio (2, "pj1Etratio");
	//MakeHistExptBackRatio (2, "pj2Etratio");
	//MakeHistExptBackRatio (2, "j1j2Etratio");
	//MakeHistExptBackRatio (2, "j1jetEta");
	//MakeHistExptBackRatio (2, "j2jetEta");
};
