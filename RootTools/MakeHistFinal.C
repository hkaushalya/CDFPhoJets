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

using namespace std;


TF1* MakeHistExptBackRatio (const int jets, const std::string& name, const std::string& title,
				 const float lolimit, const float hilimit, const int rebin,
				 const std::string& path, const int qcd_opt, const int QCDerrMtd);



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

float GetMax(const float def, const float up, const float down)
{
	float m1 = fabs(def - up);
	float m2 = fabs(def - down);
	if (m1>m2) return m1;
	else return m2;
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
TH1* get_variableBins (TH1 *hist, int jets, std::string name, float lolimit, float hilimit, bool use_errors, bool debug=false)
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




TH1 *move_overflow (TH1 *hist, float lolimit, float hilimit, bool use_errors, bool debug=false)
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


	TH1 *result = new TH1F ((std::string (hist->GetName())+"_tmp").c_str(),hist->GetTitle(),(int)(fabs(hilimit - lolimit)/hist->GetBinWidth(1)),lolimit, hilimit);
	result->SetDirectory (NULL);

	if (use_errors) result->Sumw2();

	unsigned target = 1;			//first move all undeflow to first bin
//	
	for (unsigned bin = 1; bin <= unsigned (hist->GetNbinsX()); ++ bin)
	{
		if (hist->GetXaxis()->GetBinCenter(bin-1) > lolimit) target++;
		if (hist->GetXaxis()->GetBinCenter(bin) > hilimit) target = result->GetNbinsX();
		
		//if ( hist->GetXaxis()->GetBinCenter(bin) > result->GetXaxis()->GetBinCenter(1)  && 
		//if ( hist->GetXaxis()->GetBinCenter(bin) > lolimit  && 
		//	( hist->GetXaxis()->GetBinCenter(bin) <= result->GetXaxis()->GetBinCenter(result->GetNbinsX())) )
		 //target++;																		// new hist bin increment only when > 1 st bin and < last bin
		
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
	std::cout << "\n\n\n i am in" << std::endl;
	assert(hist!=NULL && "hist is null!");
	for (int i=1; i <= hist->GetNbinsX(); i++) {
		float bin = hist->GetBinContent(i);
		float err = bin * 0.5;		//take 50% to be the syst
		ErrVec.push_back(err);
		std::cout << "halo err [" << i << "]="<< err <<std::endl;
	}
	std::cout << "halo err arr size = " << ErrVec.size() << std::endl;
	std::cout << "\n\n\n i am out" << std::endl;
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
	//isohist->SetLineColor(kRed);
	//trkpthist->SetLineColor(kBlue);
	//trkisohist->SetLineColor(kGreen);


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
		ErrVec.push_back(err);
	}

}

std::auto_ptr<TF1> get_error_function (unsigned jets, const std::string& name, const std::string& sign = "+")
{
	float val0 = 0, val1 = 0, width = 0, power = 2;

	if ((jets==1) && (name=="InvMass")) {
		val0 = 0.2;
		val1 = 1.5;
		width = 700;
		power = 3.0;
	}
	if ((jets==2) && (name=="InvMass")) {
		val0 = 0.2;
		val1 = 1.2;
		width = 800;
		power = 3.1;
	}

	if ((jets==1) && (name=="EtCorr")) {
		val0 = 0.15;
		val1 = 1.5;
		width = 300;
		power = 2.5;
		//val0 = 0.1;
		//val1 = 1.0;
		//width = 300;
		//power = 2;
	}

	if ((jets==2) && (name=="EtCorr")) {
		val0 = 0.14;
		val1 = 1.5;
		width = 300;
		power = 2.2;
	}


	if ((jets==1) && (name=="NJet15")) {
		val0 = 0.05;
		val1 = 1.25;
		width = 15;
		power = 2;
	}

	if ((jets==2) && (name=="NJet15")) {
		val0 = 0.1;
		val1 = 5;
		width = 10;
		power = 2.3;
	}
	if ((jets==1) && (name=="Ht")) {
		/*val0 = 0.05;
		val1 = 0.95;
		width = 700;
		power = 2;
		*/
		val0 = 0.4;
		val1 = 3.0;
		width = 700;
		power = 2;
	}
	if ((jets==2) && (name=="Ht")) {
		/*val0 = 0.05;
		val1 = 0.95;
		width = 700;
		power = 2;
		*/
		val0 = 0.03;
		val1 = 1.5;
		width = 700;
		power = 2;
	}
	if ((jets==1) && (name=="DetEta")) {
		val0 = 0.2;
		val1 = 1.5;
		width = 10;
		power = 3.0;
	}


	if (width > 0)
	{
		std::ostringstream myname, func;
		myname << "eq_" << jets << name;
		func << val0 << " + " << (val1 - val0) / pow (width, power) << " *pow(x," << power << ")";
		return std::auto_ptr<TF1>(new TF1 (myname.str().c_str(), func.str().c_str(), 0, width));
	};
	std::cout << "couldn't find errors for " << name << " " << jets << " jets" << std::endl;
	assert (false);
	return std::auto_ptr<TF1>();
};

void MakeHistFinal (const int jets, const std::string name, const std::string title,
				 const float lolimit, const float hilimit, const int rebin,
				 const std::string path, const int qcd_opt, const int QCDerrMtd)
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
	std::string sHaloDir("Hist/HALO/"+path), sCosmicDir("Hist/COSMIC/"+path), sQcdDir("Hist/SIDEBAND/"+path), sSigDir("Hist/SIGNAL/"+path);
	std::string sMcCentralDir("Hist/CENTRAL/"+path), sMcUpDir("Hist/EMJESUP/"+path), sMcDownDir("Hist/EMJESDOWN/"+path);



	fpho->cd();
	std::cout << "path="<<path<< std::endl;
	std::cout << "dir="<<sSigDir<< std::endl;

	if (! gDirectory->cd(sSigDir.c_str())) {
	 std::cout << "path not found "<< std::endl;
	 return;
	}
	
	TH1F* phojet = (TH1F*) gDirectory->FindObjectAny(name.c_str());
	if (! phojet){
	 std::cout << "hist not found in the dir" <<std::endl;
	 return;
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
	TH1F* qcdjet = (TH1F*) gDirectory->FindObjectAny(name.c_str());			//for QCD+MC combined method
	TH1F* qcdjet_100 = (TH1F*) qcdjet->Clone("qcdjet_100");						// for 100% QCD method


	//MC HISTS: jesup= jesup && emup : jesdown = jesdown && emdown

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
  cosmicjet    = (TH1F*) get_variableBins (cosmicjet, jets, name ,lolimit, hilimit,true, true);
  qcdjet    = (TH1F*) get_variableBins (qcdjet, jets, name ,lolimit, hilimit,true, false);
  mcphojet    = (TH1F*) get_variableBins (mcphojet, jets, name ,lolimit, hilimit,true, true);
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
	
	/************** FAKE PHOTON FRACTION ************************************/
	//this is used when the combination of QCD+PHO MC is used
  //fake photon fraction = 0.319+/-0.068(syst)
  // this is the amount of fake photons in the signal we select(jets faking photon) which is ~30%
  float fake_pho_frac_d = 0.319;
  float fake_pho_frac_sigma = 0.068;
  float fake_pho_frac_p = fake_pho_frac_d + fake_pho_frac_sigma;
  float fake_pho_frac_m = fake_pho_frac_d - fake_pho_frac_sigma;
  
  float qcd_d = 0.319;			//mean
  float qcd_m = 0.251;			//mean+sigma
  float qcd_p = 0.387;			//mean-sigma

  float qcd_scale = qcd_d;					//these background need to be scaled to these values
  float phomc_scale  = 1 - qcd_d;

  std::string qcd_str, mc_str;
  std::ostringstream qcdnum, mcnum;
  qcdnum << "QCD (#gamma sideband, " << qcd_scale<< " of data)";
  mcnum << "#gamma MC (" << phomc_scale << " of data)";
  qcd_str = qcdnum.str();		// to be used in the legend of the final plot
  mc_str = mcnum.str();

	/************************************************************************/
  //******************* NORMALIZING *************************************
  float haloEst 		= 0;
  float cosmicEst 	= 0;
  if (jets ==1 ) {
    haloEst = 9;				//these estimates are for the full dataset
    cosmicEst = 110;
  } else if (jets ==2 ) {
    haloEst = 1 ;
    cosmicEst = 7 ;
  } else {
	  std::cout << __FILE__ <<"::"<<__LINE__<<":: Invalid number of jets !" << std::endl;
	  exit (1);
  }


  float haloNorm   = 0;
  float cosmicNorm = 0;
  float qcdNorm    = 0; 
  float mcphoNorm	 = 0;
  float qcd100Norm = 0;		//for the plots that we use 100% of QCD sideband

  haloNorm    = haloEst / halojet->Integral();
  cosmicNorm  = cosmicEst / cosmicjet->Integral();
  qcdNorm     = ((phojet->Integral()) * qcd_scale)/qcdjet->Integral();
  qcd100Norm  = ((phojet->Integral()))/qcdjet_100->Integral(); 
  mcphoNorm   = ((phojet->Integral()) * phomc_scale)/mcphojet->Integral();


  //JES PHO MC
  float mcphoJESUPNorm = ((phojet->Integral()) * phomc_scale)/mcphojetJESUP->Integral();
  float mcphoJESDOWNNorm = ((phojet->Integral()) * phomc_scale)/mcphojetJESDOWN->Integral();


	//float dataLum=2542.9;	//nb-1
	float dataLum=2043.0;	//nb-1
	
  float kFac = 1.4;
//  float zeenorm = (dataLum/47605)*kFac;                   // for EWK mc see. log book#2 pp.72
//  float zmmnorm = (dataLum/49577)*kFac;                   // for EWK mc see. log book#2 pp.72
//  float zttnorm = (dataLum/28991)*kFac;                   // for EWK mc see. log book#2 pp.72
//  float wennorm = (dataLum/15408)*kFac;                   // for EWK mc see. log book#2 pp.72
//  float wmnnorm = (dataLum/7704)*kFac;                   // for EWK mc see. log book#2 pp.72
//  float wtnnorm = (dataLum/7704)*kFac;                   // for EWK mc see. log book#2 pp.72

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

  //***************END  NORMALIZING *************************************

  //*************** SUBTRACT HALO AND COSMIC FROM QCD BACKGROUND *******
  // this will avoid the double counting
  //I have removed halo/cosmics from both signal and qcd templates,
  // I need to subtract the templates by normalizing them to signal/qcd.
  //but since I have an estimate of the expected number of halo/cosmics 
  // events, I can just subtract the template. halo selection is orthogonal to 
  // photon or qcd selection. cosmic fraction ib signal/qcd can be asumed to
  // be same, I think.
  
  qcdjet->Add(cosmicjet, -1);
  qcdjet->Add(halojet, -1);
  qcdjet_100->Add(halojet, -1);
  qcdjet_100->Add(cosmicjet, -1);

  //*********************************************************************
	//****************** SYSTEMATICS *************************************

  std::vector<float> cosmicErr;   //relative error for each bin is stored here
  std::vector<float> haloErr;

  GetCosmicErr(cosmicErr, cosmicjet);
  GetHaloErr(haloErr, halojet);


	// Luminosity Error +/-6%
  float fact = fabs(1.0-100.0/106.0);
  float zeeLumErr = zeenorm * fact;
  float zmmLumErr = zmmnorm * fact;
  float zttLumErr = zttnorm * fact;
  float wenLumErr = wennorm * fact;
  float wmnLumErr = wmnnorm * fact;
  float wtnLumErr = wtnnorm * fact;

	// JES and QCD mixture 
	// here I am trying to get two things done at once. 1. JES syst 2. fake pho fraction syste using pho mc and qcd

  TH1F *hist_err = (TH1F*) halojet->Clone("hist_err");					//central values and use for SYSTEMATIC ERROR BAND
  																							// this must be identical to the final stacked plot
  TH1F *hist_jesup = (TH1F*) halojet->Clone("hist_err_jesup");			//sum of em/jes up plots, just like the sum of central plots
  TH1F *hist_jesdown = (TH1F*) halojet->Clone("hist_err_jesdown"); 	//sum of em/jes down plots
  
    hist_err->Add(cosmicjet);
    hist_err->Add(zmmjet);
    hist_err->Add(wtnjet);
    hist_err->Add(zttjet);
    hist_err->Add(wmnjet);
    hist_err->Add(wenjet);
    hist_err->Add(zeejet);
	  
    hist_jesup->Add(cosmicjet);
    hist_jesup->Add(zmmjetJESUP);
    hist_jesup->Add(wtnjetJESUP);
    hist_jesup->Add(zttjetJESUP);
    hist_jesup->Add(wmnjetJESUP);
    hist_jesup->Add(wenjetJESUP);
    hist_jesup->Add(zeejetJESUP);
 
    hist_jesdown->Add(cosmicjet);
    hist_jesdown->Add(zmmjetJESDOWN);
    hist_jesdown->Add(wtnjetJESDOWN);
    hist_jesdown->Add(zttjetJESDOWN);
    hist_jesdown->Add(wmnjetJESDOWN);
    hist_jesdown->Add(wenjetJESDOWN);
    hist_jesdown->Add(zeejetJESDOWN);
 
 	std::vector<float> qcdmcMixErr;		//this will hold the systematics for either QCD mthd 1 or 2 for a given case
	TH1F* hist_qcd100Syst;					//to get systematics when 100% QCD is used
 	TH1F* hist_varyMcQcdMixSyst;			// this will be used to get the systematic error in the fake photon fraction
												// when I use the QCD+MC combination
												// I need a copy of the central hists before the qcd+mc (central) is added

	if (QCDerrMtd == 1)		//_______________________ plots that use 100% qcd
  	{
  		hist_qcd100Syst = (TH1F*) hist_err->Clone("hist_qcd100");

		hist_err->Add(qcdjet_100);
		hist_jesup->Add(qcdjet_100);
		hist_jesdown->Add(qcdjet_100);
	
		//systematics for fake photon fraction
    	std::string abspath("Hist/"+path);
	   GetQCD100Err(qcdmcMixErr, hist_qcd100Syst, name, abspath, rebin);

		
 	
	} else if (QCDerrMtd == 2)  //_____________________ plots that use 30%70%
	{
  		hist_varyMcQcdMixSyst = (TH1F*) hist_err->Clone("hist_varyMcQcdMixSyst");	//use this to get a varied mix of qcd+pho mc
		
		hist_err->Add(qcdjet);
		hist_err->Add(mcphojet);
 
		hist_jesup->Add(qcdjet);
		hist_jesup->Add(mcphojetJESUP);
 
		hist_jesdown->Add(qcdjet);
		hist_jesdown->Add(mcphojetJESDOWN);


		//get the systematics for the mixture
		//here i am using some dummy mixture for now. change them wisely!

		//clones to be used in the syst calculations. need this before normalizing them
		TH1F* hist_qcdVaried = (TH1F*) qcdjet->Clone("hist_qcdVaried");
		TH1F* hist_mcphoVaried = (TH1F*) mcphojet->Clone("hist_mcphoVaried");

		//not sure what i am doing in the middle of the night???	
		float mcphoMixVaried_Norm = ((phojet->Integral()) * (1 - fake_pho_frac_p) )/mcphojet->Integral();
		hist_mcphoVaried->Scale(mcphoMixVaried_Norm);
		
		float qcdMixVaried_Norm = ((phojet->Integral()) * fake_pho_frac_p)/qcdjet->Integral();
		hist_qcdVaried->Scale(qcdMixVaried_Norm);

		hist_varyMcQcdMixSyst->Add(hist_qcdVaried);
		hist_varyMcQcdMixSyst->Add(hist_mcphoVaried);

		for (int i=1; i <= hist_varyMcQcdMixSyst->GetNbinsX(); i++)
		{
			float err = fabs(hist_varyMcQcdMixSyst->GetBinContent(i) - hist_err->GetBinContent(i));
			qcdmcMixErr.push_back(err);
		}
	 
  	} else 
	{
		std::cerr << "bad value: QCDerrMtd=" << QCDerrMtd << std::endl;
      return;
	};

	// NOW GET THE SYSTEMATICS FROM JES BY comparing central values with the JES UP/DOWN
	// USE THE HIST USED TO PUT THE ERROR BAND AS THE CENTRAL VALUE
 	std::vector<float> JESerr;
	for (int i=1; i <= hist_err->GetNbinsX(); i++)
	{
			JESerr.push_back(GetMax(hist_err->GetBinContent(i),
											hist_jesup->GetBinContent(i),
											hist_jesdown->GetBinContent(i)));
	}

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
	std::cout << "hist_err nbins = " << hist_err->GetNbinsX() << std::endl;
	std::cout << "stat err arr size= " << statErr.size() << std::endl;
	std::cout << "cosmic err arr size= " << cosmicErr.size() << std::endl;
	std::cout << "halo err arr size= " << haloErr.size() << std::endl;
	std::cout << "zeeJes err arr size= " << JESerr.size() << std::endl;
	std::cout << "QCDMC err arr size= " << qcdmcMixErr.size() << std::endl;


  // NOW COLLECT ALL ERRORS AND PUT THEM TOGETHER FOR ONE FINAL NUMBER
  std::vector<float> ERRORS;

	float jesmax=0, jesmaxbin = 0, statmax=0, statmaxbin, qcdmcmax=0 , qcdmcmaxbin=0;
	float halomax=0, halomaxbin=0;
	float cosmicmax=0, cosmicmaxbin=0;
  for (unsigned int i=0; i < cosmicErr.size(); i++) {
    float sum = pow(cosmicErr[i],2) + pow(haloErr[i],2) + pow(zeeLumErr,2) + pow(zmmLumErr,2)
      + pow(zttLumErr,2) + pow(wenLumErr,2) + pow(wmnLumErr,2) + pow(wtnLumErr,2) + pow(qcdmcMixErr[i],2)
		+ pow(JESerr[i],2) + pow(statErr[i],2);
		
	if (cosmicErr.size() - 2 == i) {
    if (jesmax < JESerr.at(i)) {
		 jesmax =JESerr.at(i);
		 jesmaxbin = hist_err->GetBinLowEdge(i+1);
	 }

    if (statmax < statErr.at(i)) {
	 	statmax =statErr.at(i);
		statmaxbin = hist_err->GetBinLowEdge(i+1);
	 }

    if (halomax < haloErr.at(i)) {
	 	halomax =haloErr.at(i);
		halomaxbin = hist_err->GetBinLowEdge(i+1);
	 }

    if (cosmicmax < cosmicErr.at(i)) {
	 	cosmicmax = cosmicErr.at(i);
		cosmicmaxbin = hist_err->GetBinLowEdge(i+1);
	 }

	}
    float err = sqrt(sum);
    ERRORS.push_back(err);
    hist_err->SetBinError(i+1, err);
  }

	std::cout << "jes max, bin = " << jesmax << ", " << jesmaxbin << std::endl;
	std::cout << "stat max,bin = " << statmax << ", " << statmaxbin  << std::endl;
	std::cout << "halo max, bin= " << halomax << ", " << halomaxbin  << std::endl;
	std::cout << "cosmic max, bin= " << cosmicmax << ", " << cosmicmaxbin  << std::endl;

 //correct_errors(hist_err);
  //std::auto_ptr<TF1> func = get_error_function (jets, name);
  TF1* func = MakeHistExptBackRatio (jets, name, title, lolimit, hilimit, rebin, path, qcd_opt, QCDerrMtd);

  for (int bin = 1; bin <= hist_err->GetNbinsX(); ++ bin)
  {
	 hist_err->SetBinError (bin, func->Eval (hist_err->GetXaxis()->GetBinCenter (bin)) * hist_err->GetBinContent (bin));
  };


 
  std::cout << "Ewk: " << (wtnjet->Integral()+zttjet->Integral()+wenjet->Integral()+zeejet->Integral()) << std::endl;
  std::cout << "Ewk JESUP: " << (wtnjetJESUP->Integral()+zttjetJESUP->Integral()+wenjetJESUP->Integral()+zeejetJESUP->Integral()) << std::endl;
  std::cout << "Ewk JESDOWN: " << (wtnjetJESDOWN->Integral()+zttjetJESDOWN->Integral()+wenjetJESDOWN->Integral()+zeejetJESDOWN->Integral()) << std::endl;
  std::cout << "QCD frac =" << qcdjet->Integral() << std::endl;
  std::cout << "Pho frac =" << mcphojet->Integral() << std::endl;

	

	//********************************************************************

  //************* SET HIST COLORS ***************************************

	Int_t linecolor=47;

  mcphojet->SetLineColor(linecolor);
  mcphojet->SetFillColor(6);
  //qcdjet->SetLineColor(6);
  //qcdjet->SetFillColor(6);
  //qcdjet->SetLineColor(42);
  //qcdjet->SetFillColor(42);
  qcdjet->SetLineColor(linecolor);
  qcdjet->SetFillColor(kYellow);
  cosmicjet->SetLineColor(linecolor);
  cosmicjet->SetFillColor(kGreen);
  halojet->SetLineColor(linecolor);
  halojet->SetFillColor(12);
  /*
    zeejet->SetLineColor(1);
    zeejet->SetFillColor(50);
    zmmjet->SetLineColor(1);
    zmmjet->SetFillColor(16);
    zttjet->SetLineColor(1);
    zttjet->SetFillColor(19);
    wenjet->SetLineColor(1);
    wenjet->SetFillColor(41);
    wmnjet->SetLineColor(1);
    wmnjet->SetFillColor(42);
    wtnjet->SetLineColor(1);
    wtnjet->SetFillColor(38);
  */
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

  phojet->SetLineColor(kBlack);
  phojet->SetMarkerStyle (8);

  //this is the line histo showing 100% of QCD 0% MC pho
  qcdjet_100->SetLineColor(kRed);
  qcdjet_100->SetFillColor(kYellow);
  hist_err->SetFillStyle(3001);
  hist_err->SetFillColor(13);
  hist_err->SetLineColor(10);
 

  
  //********* END SET HIST COLORS ***************************************



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

  } else if (QCDerrMtd == 2) {//plots that use 30%70%
    hs->Add(halojet);
    hs->Add(cosmicjet);
    hs->Add(zmmjet);
    hs->Add(wtnjet);
    hs->Add(zttjet);
    hs->Add(wmnjet);
    hs->Add(wenjet);
    hs->Add(zeejet);
    hs->Add(qcdjet);
    hs->Add(mcphojet);
  } else {
    std::cerr << "bad value: QCDerrMtd=" << QCDerrMtd << std::endl;
    return;
  };

  TLegend *leg = new TLegend (0.5,0.6,0.9,0.9);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  std::string str_pho,str_cosmic,str_halo,str_zee,str_zmm,str_ztt,str_wen,str_wmn,str_wtn;

  if (jets==1) {
    str_pho 	= "Data (#gamma + >=1 Jet)";	
    str_cosmic = "#gamma^{cosmic} + >=1 Jet";
    str_zee 	= "Z->ee MC (e + >=1 Jet)";
    str_wen 	= "W->e#nu MC (e + >=1 Jet)";
    str_ztt 	= "Z->#tau#tau MC (e + >=1 Jet)";
    str_wtn 	= "W->#tau#nu MC (e + >=1 Jet)";
    str_wmn 	= "W->#mu#nu MC (e + >=1 Jet)";
    str_zmm 	= "Z->#mu#mu MC (e + >=1 Jet)";
    str_halo 	= "#gamma^{halo} + >=1 Jet";
  }
  if (jets==2) {
    str_pho 	= "Data (#gamma + >=2 Jets)";	
    str_cosmic = "#gamma^{cosmic} + >=2 Jets";
    str_zee 	= "Z->ee MC (e + >=2 Jets)";
    str_wen 	= "W->e#nu MC (e + >=2 Jets)";
    str_ztt 	= "Z->#tau#tau MC (e + >=2 Jets)";
    str_wtn 	= "W->#tau#nu MC (e + >=2 Jets)";
    str_wmn 	= "W->#mu#nu MC (e + >=2 Jets)";
    str_zmm 	= "Z->#mu#mu MC (e + >=2 Jets)";
    str_halo 	= "#gamma^{halo} + >=2 Jets";
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
    //leg->AddEntry(qcdjet,"QCD (#gamma sideband+MC)");
    //leg->AddEntry(qcdjet,"QCD (#gamma sideband)");
    //leg->AddEntry(zeejet,str_zee.c_str());
    //leg->AddEntry(wenjet,str_wen.c_str());
    //leg->AddEntry(zttjet,str_ztt.c_str());
    //leg->AddEntry(wtnjet,str_wtn.c_str());
    leg->AddEntry(wtnjet,"EWK MC");
    leg->AddEntry(cosmicjet,str_cosmic.c_str());
    leg->AddEntry(halojet,str_halo.c_str());
    leg->AddEntry(hist_err,"Systematic Uncertainty");
  } else {
    std::cerr << "bad value: QCDerrMtd=" << QCDerrMtd << std::endl;
    return;
  };

  leg->SetBorderSize (1);
  leg->SetFillColor (10);
	gStyle->SetCanvasColor (10);
	gStyle->SetCanvasBorderSize (0);
	gStyle->SetCanvasBorderMode (0);
				 
	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (1);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);


  new TCanvas();

  //gPad->Print();
  gStyle->SetTextFont(132);
  hs->SetMinimum(0.05);
  if (jets == 1) hs->SetMaximum(0.5e6);
  if (jets == 2) hs->SetMaximum(0.5e5);

	
  
 // hs->SetTitle("CDF Run II Preliminary 2.0 fb^{-1}");
  gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  hs->Draw("HIST");		//need this as I am calling sumw2 for all hists. if not it will draw all hists with error bars

  hs->GetXaxis()->SetTitle (title.c_str());
  std::ostringstream ytitle;
  if (name == "EtCorr") ytitle << "#Delta N/ #Delta E_{T}";
  if (name == "InvMass") ytitle << "#Delta N/ #Delta M";
  ytitle << "    Bin Size = " << phojet->GetBinWidth(1);
  hs->GetYaxis()->SetTitle (ytitle.str().c_str());

 int labelfont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
 int titlefont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
 hs->GetXaxis()->SetLabelFont(labelfont);
 hs->GetYaxis()->SetLabelFont(labelfont);
 hs->GetXaxis()->SetTitleFont(titlefont);
 hs->GetYaxis()->SetTitleFont(titlefont);
 hs->GetXaxis()->CenterTitle(true);
 hs->GetYaxis()->CenterTitle(true);

  
  hist_err->SetDrawOption("HIST");
  hist_err->Draw("sameE2");
  phojet->Draw("same");
  leg->Draw ();

	TPaveText *tp = new TPaveText(0.5,0.91,0.9,0.99,"NDC");
	tp->SetLineColor(10);
	tp->SetTextFont(titlefont);
	tp->AddText("CDF Run II Preliminary 2.0 fb^{-1}");
	tp->Draw();

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
	 

};

void MakeHistFinal (int jets,std::string name,  std::string opt = "")
{
  const bool opt_qcdup = (opt.find ("QCDUP") != std::string::npos);
  const bool opt_qcddown = (opt.find ("QCDDOWN") != std::string::npos);


  int qcd_scale = 0;
  TVirtualPad *pad_old = gPad;	

	
  if (opt_qcdup) qcd_scale = 1;
  if (opt_qcddown) qcd_scale = -1;
/*
  if (jets == 1) {
    if (name == "InvMass") 	MakeHistFinal (jets, "InvMass","Invariant Mass (#gamma,Lead Jet) (GeV/c^{2})",60,640,4,"1Jet/PhotonLeadJet", qcd_scale, 2);
    else if (name == "pet") 	MakeHistFinal (jets, "EtCorr","E_{T}^{#gamma} (GeV)",0,300,5,"1Jet/Photon", qcd_scale,2);
    else if (name == "met")	MakeHistFinal (jets, "Met","MET (GeV)",0,200,1,"1Jet/Event", qcd_scale,1);
    //else if (name == "met")	MakeHistFinal (jets, "Met","MET",0,200,1,"1Jet/Event", qcd_scale,2);
    else if (name == "njet") 	MakeHistFinal (jets, "NJet15","Jet Multiplicity (E_{T}>15GeV)",0,15,1,"1Jet/Event", qcd_scale,1);
    //else if (name == "njet") 	MakeHistFinal (jets, "NJet15","Jet Multiplicity (E_{T}>15GeV)",0,15,1,"1Jet/Event", qcd_scale,2);
    else if (name == "ht")  	MakeHistFinal (jets, "Ht","H_{T} (GeV)",0,700,4,"1Jet/Event", qcd_scale,2);
    else if (name == "jet") 	MakeHistFinal (jets, "EtCorr","E_{T}^{lead Jet} (GeV)",0,300,5,"1Jet/LeadJet", qcd_scale,2);
    else if (name == "jetEta") 	MakeHistFinal (jets, "DetEta","Detector #eta ^{lead Jet}",-5,5,4,"1Jet/LeadJet", qcd_scale,2);
    else if (name == "etratio") MakeHistFinal (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,2,"1Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "delphi") 	MakeHistFinal (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,5,1,"1Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "delr") MakeHistFinal (jets, "DelR", "#Delta R^{#gamma, lead Jet}",0,7,1,"1Jet/PhotonLeadJet", qcd_scale,2);
  }
  if (jets == 2) {
    if (name == "InvMass") 	MakeHistFinal (jets, "InvMass","Invariant Mass (#gamma,Two Lead Jets) (GeV/c^{2})",80,780,4,"2Jet/Photon2Jets", qcd_scale,2);
    else if (name == "Met")	MakeHistFinal (jets, "Met","MET",0,125,1,"2Jet/Event", qcd_scale,1);
    //else if (name == "Met")	MakeHistFinal (jets, "Met","MET",0,125,1,"2Jet/Event", qcd_scale,2);
    else if (name == "jetsInvMass") 	MakeHistFinal (jets, "InvMass","Invariant Mass(Lead two Lead Jets) (GeV/c^{2})",0,650,5,"2Jet/2Jets", qcd_scale,2);
    //else if (name == "NJet15") 	MakeHistFinal (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale,1);
    //else if (name == "NJet15") 	MakeHistFinal (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale,2);
    else if (name == "PhoEt") 	MakeHistFinal (jets, "EtCorr","E_{T}^{#gamma} (GeV)",0,300,5,"2Jet/Photon", qcd_scale,2);
    else if (name == "pj1DelPhi") 	MakeHistFinal (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,4,1,"2Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "pj2DelPhi") 	MakeHistFinal (jets, "DelPhi","#Delta #phi^{#gamma,2nd lead jet}",0,4,1,"2Jet/Photon2ndLeadJet", qcd_scale,2);
    else if (name == "j1j2DelPhi") 	MakeHistFinal (jets, "DelPhi","#Delta #phi^{lead jet,2nd lead jet} (GeV)",0,4,1,"2Jet/2Jets", qcd_scale,2);
    else if (name == "Ht")  	MakeHistFinal (jets, "Ht","H_{T} (GeV)",0,700,4,"2Jet/Event", qcd_scale,2);
    //else if (name == "Ht")  	MakeHistFinal (jets, "Ht","H_{T}",0,700,4,"2Jet/Event", qcd_scale,2);
    else if (name == "j1Et") 	MakeHistFinal (jets, "EtCorr","E_{T}^{lead Jet} (GeV)",0,300,5,"2Jet/LeadJet", qcd_scale,2);
    else if (name == "j2Et") 	MakeHistFinal (jets, "EtCorr","E_{T}^{2nd lead Jet} (GeV)",0,300,5,"2Jet/SecondLeadJet", qcd_scale,2);
    else if (name == "jetEta") 	MakeHistFinal (jets, "DetEta","Detector #eta ^{lead Jet}",-5,5,4,"2Jet/LeadJet", qcd_scale,2);
    else if (name == "pj1Etratio") MakeHistFinal (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,1,"2Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "pj2Etratio") MakeHistFinal (jets, "EtRatio","E_{T} ratio (2nd lead Jet/ #gamma)",0,10,1,"2Jet/Photon2ndLeadJet", qcd_scale,2);
    else if (name == "j1j2Etratio") MakeHistFinal (jets, "EtRatio","E_{T} ratio (2nd lead Jet/lead Jet)",0,2.5,3,"2Jet/2Jets", qcd_scale,2);
  }
*/	
	// to make full scale plots
  if (jets == 1) {
    if (name == "InvMass") 	MakeHistFinal (jets, "InvMass","Invariant Mass (#gamma,Lead Jet) (GeV/c^{2})",60,660,4,"1Jet/PhotonLeadJet", qcd_scale, 2);
    //if (name == "InvMass") 	MakeHistFinal (jets, "InvMass","Invariant Mass (#gamma,Lead Jet) (GeV/c^{2})",0,1100,5,"1Jet/PhotonLeadJet", qcd_scale, 2);
    else if (name == "pet") 	MakeHistFinal (jets, "EtCorr","E_{T}^{#gamma} (GeV)",0,600,5,"1Jet/Photon", qcd_scale,2);
    //else if (name == "met")	MakeHistFinal (jets, "Met","MET",0,1000,1,"1Jet/Event", qcd_scale,1);
    else if (name == "met")	MakeHistFinal (jets, "Met","MET",0,1000,1,"1Jet/Event", qcd_scale,2);
    //else if (name == "njet") 	MakeHistFinal (jets, "NJet15","Jet Multiplicity (E_{T}>15GeV)",0,15,1,"1Jet/Event", qcd_scale,1);
    else if (name == "njet") 	MakeHistFinal (jets, "NJet15","Jet Multiplicity (E_{T}>15GeV)",0,15,1,"1Jet/Event", qcd_scale,1);
    else if (name == "ht")  	MakeHistFinal (jets, "Ht","H_{T} (GeV)",0,1500,4,"1Jet/Event", qcd_scale,1);
    else if (name == "jet") 	MakeHistFinal (jets, "EtCorr","E_{T}^{lead Jet}",0,600,5,"1Jet/LeadJet", qcd_scale,2);
    else if (name == "etratio") MakeHistFinal (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,2,"1Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "delphi") 	MakeHistFinal (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,5,1,"1Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "delr") MakeHistFinal (jets, "DelR", "#Delta R^{#gamma, lead Jet}",0,7,1,"1Jet/PhotonLeadJet", qcd_scale,2);
  }
  if (jets == 2) {
    if (name == "InvMass") 	MakeHistFinal (jets, "InvMass","Invariant Mass (#gamma,Two Lead Jets) (GeV/c^{2})",0,1100,5,"2Jet/Photon2Jets", qcd_scale,2);
    //else if (name == "Met")	MakeHistFinal (jets, "Met","MET",0,1000,1,"2Jet/Event", qcd_scale,1);
    else if (name == "Met")	MakeHistFinal (jets, "Met","MET",0,1000,1,"2Jet/Event", qcd_scale,2);
    else if (name == "jetsInvMass") 	MakeHistFinal (jets, "InvMass","Invariant Mass(Lead two Lead Jets)",0,1100,5,"2Jet/2Jets", qcd_scale,2);
    //else if (name == "NJet15") 	MakeHistFinal (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale,1);
    //else if (name == "NJet15") 	MakeHistFinal (jets, "NJet15","NJet15",0,15,1,"2Jet/Event", qcd_scale,2);
    else if (name == "PhoEt") 	MakeHistFinal (jets, "EtCorr","E_{T}^{#gamma} (GeV)",0,500,5,"2Jet/Photon", qcd_scale,2);
    else if (name == "pj1DelPhi") 	MakeHistFinal (jets, "DelPhi","#Delta #phi^{#gamma,lead jet}",0,4,1,"2Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "pj2DelPhi") 	MakeHistFinal (jets, "DelPhi","#Delta #phi^{#gamma,2nd lead jet}",0,4,1,"2Jet/Photon2ndLeadJet", qcd_scale,2);
    else if (name == "j1j2DelPhi") 	MakeHistFinal (jets, "DelPhi","#Delta #phi^{lead jet,2nd lead jet}",0,4,1,"2Jet/2Jets", qcd_scale,2);
    else if (name == "Ht")  	MakeHistFinal (jets, "Ht","H_{T} (GeV)",0,1500,4,"2Jet/Event", qcd_scale,2);
    else if (name == "j1Et") 	MakeHistFinal (jets, "EtCorr","E_{T}^{lead Jet}",0,1000,5,"2Jet/LeadJet", qcd_scale,2);
    else if (name == "j2Et") 	MakeHistFinal (jets, "EtCorr","E_{T}^{2nd lead Jet}",0,1000,5,"2Jet/SecondLeadJet", qcd_scale,2);
    else if (name == "pj1Etratio") MakeHistFinal (jets, "EtRatio","E_{T} ratio (lead Jet/ #gamma)",0,10,1,"2Jet/PhotonLeadJet", qcd_scale,2);
    else if (name == "pj2Etratio") MakeHistFinal (jets, "EtRatio","E_{T} ratio (2nd lead Jet/ #gamma)",0,10,1,"2Jet/Photon2ndLeadJet", qcd_scale,2);
    else if (name == "j1j2Etratio") MakeHistFinal (jets, "EtRatio","E_{T} ratio (2nd lead Jet/lead Jet)",0,2.5,3,"2Jet/2Jets", qcd_scale,2);
  }
	

	if (gPad != pad_old)
	{
		TCanvas *c = dynamic_cast<TCanvas*>(gPad);
      if (c)
		{
			std::ostringstream str,str1;
			str << "plot" << jets << "_" << name << ".gif";
			c->Print (str.str().c_str());
			str1 << "plot" << jets << "_" << name << ".pdf";
			c->Print (str1.str().c_str(),"pdf");
		};
	};

	
};




void MakeHistFinal ()
{
//	MakeHistFinal (1, "InvMass");
	MakeHistFinal (1, "pet");
	//MakeHistFinal (1, "met");
	//MakeHistFinal (1, "njet");
	//MakeHistFinal (1, "ht");
	//MakeHistFinal (1, "jet");
	//MakeHistFinal (1, "jetEta");
	//MakeHistFinal (1, "etratio");
	//MakeHistFinal (1, "delphi");
	//MakeHistFinal (1, "delr");
	//MakeHistFinal (2, "InvMass");
	//MakeHistFinal (2, "Met");
	//MakeHistFinal (2, "jetsInvMass");
	//MakeHistFinal (2, "NJet15");
	//MakeHistFinal (2, "PhoEt");
	//MakeHistFinal (2, "pj1DelPhi");
	//MakeHistFinal (2, "pj2DelPhi");
	//MakeHistFinal (2, "j1j2DelPhi");
	//MakeHistFinal (2, "Ht");
	//MakeHistFinal (2, "j1Et");
	//MakeHistFinal (2, "j2Et");
	//MakeHistFinal (2, "pj1Etratio");
	//MakeHistFinal (2, "pj2Etratio");
	//MakeHistFinal (2, "j1j2Etratio");
};
