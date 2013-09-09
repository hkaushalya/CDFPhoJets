#include "mergehist.h"
#include<TSystem.h>
#include<iostream>
#include<sstream>
#include<TFile.h>

//add the files. then compile and run.


void mergeHists(const int Nfiles){
	assert (Nfiles>0 && "Number of files must be >0.");
	//gSystem->CompileMacro("~/samantha/RootTools/mergehist.C","f");

	MergeHist merge;

	//specify all input files
	//merge.input("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/MetTest/WenuToEleJets_WenMc/MetTest_PhotonMC.root_1");
	std::string dir("/mnt/autofs/misc/nbay03.a/samantha/RESULTS/MetTest/NewStuff/04152009_UsingHadJets/InclNjet15_MetSig2/");

	//for (int i=1; i<=Nfiles; ++i) 
/*	for (int i=41; i<=50; ++i) 
	{
		std::stringstream file;
		file << dir << "MetTest_PhotonMC.root_" << i;
		TFile f(file.str().c_str());
		if (f.IsZombie())
		{
			std::cout << file.str() << " Not found!" << std::endl;
			exit (1);
		} else {
			std::cout << "Adding file " << file.str() << std::endl;
		}
		
		merge.input(file.str());
	}
*/

	std::stringstream f1,f2,f3,f4,f5;
	f1 << dir << "MergedHists.root";
	f2 << dir << "MergedHists2.root";
	f3 << dir << "MergedHists3.root";
	f4 << dir << "MergedHists4.root";
	f5 << dir << "MergedHists5.root";

	//merge.input(file.str());
	merge.input(f1.str());
	merge.input(f2.str());
	merge.input(f3.str());
	merge.input(f4.str());
	merge.input(f5.str());
	std::stringstream outpath;
	//outpath << dir << "MergedHists5.root";
	outpath << dir << "MERGEDHISTS.root";
	merge.output(outpath.str());
	std::cout << "Output written to " << outpath.str() << std::endl;
}
