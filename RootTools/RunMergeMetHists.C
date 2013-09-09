{

	//gSystem->CompileMacro("/cdf/home/samantha/samantha/RootTools/MergeMetHists.C");
	//gSystem->CompileMacro("/cdf/home/samantha/samantha/RootTools/MergeA2DPlot.C","f");
	gSystem->CompileMacro("/cdf/home/samantha/samantha/RootTools/MergeA1DPlot.C","f");
	

	std::vector<std::string> files;
	for (int i =1; i<=100; ++i)
	{
		std::stringstream file;
		file << "MetTest_PhotonMC.root_" << i;
		files.push_back(file.str());
	}

	std::vector<std::string> hists, paths;
	hists.push_back("debug_PhoRes");
	paths.push_back("Ana/JetFilterV2/Hist/Ana_data/");
//	hists.push_back("debug_PhoRes");
//	paths.push_back("Ana/JetFilterV2/Hist/Ana_def/");
	
	MergeA1DPlot(files,hists,paths,"MEtSig0.root",1);
	

	hists.clear();
	paths.clear();
	hists.push_back("debug_PhoRes");
	paths.push_back("Ana/JetFilterV2_2/Hist/Ana_data/");
//	hists.push_back("debug_PhoRes");
//	paths.push_back("Ana/JetFilterV2_2/Hist/Ana_def/");
	MergeA1DPlot(files,hists,paths,"MEtSig2.root",1);

	return;

	//void MergeA2DPlot(const std::vector<std::string> file, const std::vector<std::string> hist, const std::vector<std::string> path, const std::string outputfile="MergedHists.root")

//MergeMetHists (100,"MetTest_PhotonMC.root_");
//MergeMetHists ();
///MergeMetHists (100,"MetTest_PhotonMC.root_","PYTHIA #gamma MC: Uncl=1, MetSig=6", "Met",-1,-1,4, "PhoMC_Uncl1MetSig6_Met.gif");
//MergeMetHists (100,"MetTest_PhotonMC.root_","PYTHIA #gamma MC: Uncl=1, MetSig=6", "MetSig",-1,-1,4, "PhoMC_Uncl1MetSig6_MetSig.gif");

}
