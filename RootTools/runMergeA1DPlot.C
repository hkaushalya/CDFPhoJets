/********************************************************************
 * *****************************************************************/

{
	std::vector<std::string> files,hists,path;
	
	gSystem->CompileMacro("~/samantha/RootTools/MergeA1DPlot.C","f");
	//gSystem->CompileMacro("~/samantha/RootTools/MergeA2DPlot.C","f");
	//gSystem->CompileMacro("~/samantha/RootTools/OverlayHists.C","f");
	//gSystem->CompileMacro("~/samantha/RootTools/MakeA1DRatioHist.C","f");
	int nfiles = 0;
	for (unsigned int i=0; i <100; ++i)
	{
		std::stringstream name;
		name << "PhoJets_" << i << "Of100.root";
		files.push_back(name.str());
		nfiles++;
	}
	std::cout << "given files = " << nfiles << std::endl;
	hists.push_back("JetMetDelPhi");
	path.push_back("Ana/PhoJets/Hist/1Jet/LeadJet/");
	MergeA1DPlot(files, hists, path,"MergeResult.root");	


}
