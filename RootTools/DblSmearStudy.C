/********************************************************************
 * This is ROOT macro I created to study the double smearing effect
 * on the MET model. This simple runs over the CAF output files,
 * merge the relavent histograms, and overlay them.
 * this version looks at different region of the DelPhi(Jet,Met)
 * distribution, low,middle, and high regions. 
 * ******************************************************************
 *  $Id: DblSmearStudy.C,v 1.1 2010/02/18 16:47:11 samantha Exp $ 
 *
 *  $Log: DblSmearStudy.C,v $
 *  Revision 1.1  2010/02/18 16:47:11  samantha
 *  Created this to study the double smear effect on Met Model. This version looks
 *  different regions of the DelPhi(jet,Met) distributions separately by overlaying
 *  histograms.
 *
 * *****************************************************************/

{
	std::vector<std::string> files,hists,path;
	
	gSystem->CompileMacro("~/samantha/RootTools/MergeA1DPlot.C","f");
	//gSystem->CompileMacro("~/samantha/RootTools/MergeA2DPlot.C","f");
	gSystem->CompileMacro("~/samantha/RootTools/OverlayHists.C","f");
	//gSystem->CompileMacro("~/samantha/RootTools/MakeA1DRatioHist.C","f");
	//for (unsigned int i=1; i <=50; ++i)
	for (unsigned int i=1; i <=100; ++i)
	{
		std::stringstream name;
		name << "MetTest_PhotonMC.root_" << i;
		files.push_back(name.str());
	}

	hists.push_back("debug_PhoPtRatio_DelPhiJetMetTop");
	path.push_back("Ana/JetFilterV2/Hist/Ana_data");
	hists.push_back("debug_PhoPtRatio_DelPhiJetMetMid");
	path.push_back("Ana/JetFilterV2/Hist/Ana_data");
	hists.push_back("debug_PhoPtRatio_DelPhiJetMetLow");
	path.push_back("Ana/JetFilterV2/Hist/Ana_data");
	MergeA1DPlot(files, hists, path,"DATA_MEtSig0.root");	

	hists.clear();
	path.clear();
	
	hists.push_back("debug_PhoPtRatio_DelPhiJetMetTop");
	path.push_back("Ana/JetFilterV2/Hist/Ana_def");
	hists.push_back("debug_PhoPtRatio_DelPhiJetMetMid");
	path.push_back("Ana/JetFilterV2/Hist/Ana_def");
	hists.push_back("debug_PhoPtRatio_DelPhiJetMetLow");
	path.push_back("Ana/JetFilterV2/Hist/Ana_def");
	MergeA1DPlot(files, hists, path,"DEF_MEtSig0.root");	



	std::vector<std::string> legends;
	TFile *f = new TFile("DATA_MEtSig0.root");
	TH1 *h1 = dynamic_cast<TH1*> (gDirectory->FindObjectAny("debug_PhoPtRatio_DelPhiJetMetTop")); //data
	TH1 *h10 = dynamic_cast<TH1*> (gDirectory->FindObjectAny("debug_PhoPtRatio_DelPhiJetMetMid")); //data
	TH1 *h1000 = dynamic_cast<TH1*> (gDirectory->FindObjectAny("debug_PhoPtRatio_DelPhiJetMetLow")); //data
	assert (h1 != NULL);
	h1->Print();
	assert (h10 != NULL);
	assert (h100 != NULL);
	legends.push_back("PYTHIA MC");


	TFile *g = new TFile("DEF_MEtSig0.root");
	TH1 *h2 = dynamic_cast<TH1*> (gDirectory->FindObjectAny("debug_PhoPtRatio_DelPhiJetMetTop")); //bg
	TH1 *h20 = dynamic_cast<TH1*> (gDirectory->FindObjectAny("debug_PhoPtRatio_DelPhiJetMetMid")); //
	TH1 *h200 = dynamic_cast<TH1*> (gDirectory->FindObjectAny("debug_PhoPtRatio_DelPhiJetMetLow")); //
	assert (h2 != NULL);
	h2->Print();
	assert (h20 != NULL);
	assert (h200 != NULL);
	legends.push_back("Fake MEt");

	
	//MakeA1DRatioHist(h1,h2,4,legends,
	std::string title("#gamma^{E_{T}>30GeV,#eta <0.95}+>=Jets^{Et>15GeV,#eta<3.0}, #slash{E}_{T}Sig>0");
	std::string ytitle("#frac{With JER Offset}{No JER Offset} - 1");
	MakeA1DRatioHist(h1,h2,4,legends,title,"ME_{T}", ytitle);
	MakeA1DRatioHist(h10,h20,4,legends,title,"ME_{T} All", ytitle);
	MakeA1DRatioHist(h100,h200,4,legends,title,"ME_{T}-Sig", ytitle);
	MakeA1DRatioHist(h1000,h2000,4,legends,title,"Lead Jet E_{T}", ytitle);

}
