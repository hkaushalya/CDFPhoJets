#include<TFile.h>
#include<iostream>
#include<TDirectory.h>
#include<vector>
#include<TSystem.h>

void myMetMerge(){
	gSystem->CompileMacro("~/samantha/RootTools/MergeA1DPlot.C","f");
	//gSystem->CompileMacro("~/samantha/RootTools/MakeARootFile.C","f");
	//gSystem->CompileMacro("~/samantha/RootTools/MakeARootFile.C");
	gSystem->CompileMacro("~/samantha/RootTools/OverlayHists.C","f");
	
	std::vector<std::string> files,hists,vPath;

	for (unsigned int i=1; i <=49; ++i)
	{
		std::stringstream name;
		name << "MetTest_PhotonMC.root_" << i;
		files.push_back(name.str());
	}


	TFile f("TestFile.root","RECREATE");
	MakeARootFile(f);

	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_data/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_def/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_def/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_ue1/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_ue2/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_ue3/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_ue4/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_ue5/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_ue6/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_ue7/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_ue8/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_ue9/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer1/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer2/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer3/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer4/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer5/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer6/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer7/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer8/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer9/");
	vPath.push_back("/Ana/JetFilterV2/Hist/Ana_jer10/");


	
	hists.push_back("Etjet");
	hists.push_back("Etjet2");
	hists.push_back("MetAll");
	hists.push_back("MetSig");
	hists.push_back("Met");
	hists.push_back("dPhi3");
	hists.push_back("Njet15");
	hists.push_back("Et1");
	hists.push_back("dPhi1");
	hists.push_back("Ht");
	

	std::cout << files.size() << "\t" << hists.size() << "\t" << vPath.size() << std::endl;
	

	//for (unsigned i=0; i<vPath.size(); ++i)
	for (unsigned i=0; i<1; ++i)
	{
		//for (unsigned int j=0; j<hists.size(); j++)
		for (unsigned int j=0; j<1; j++)
		{
//			std::cout << " vPath, hist = " << vPath[i] << ", " << hists[j] << std::endl;
			f.cd(vPath[i]);
			MergeA1DPlot(files, vPath[i], hists[j],gDirectory);
		}
	}

	f.Close();

};
