#include <TH1.h>
#include <iostream>
#include <string>
#include <TFile.h>
/*
GIVEN THE ROOT FILE, PATH AND THE HIST(1D) THIS WILL DUMP ALL BINS
WITH NON ZERO VALUES
*/

void DumpHistBins(std::string file, std::string path, std::string histname, bool showoverflowbins = false)

{

	TFile f(file.c_str());
	if (f.IsZombie()) {
		std::cout << " ERROR: " << file << " cannot be found." << std::endl;
		exit (1);
	}

	if (f.cd(path.c_str())) {
		TH1* hist;
		gDirectory->GetObject(histname.c_str(), hist);
		if (hist) {
			std::cout << "Dump for hist " << histname << " in dir " 
						<< path << " in file " << file << std::endl;
			for (int i=1 ; i <= hist->GetNbinsX(); ++i) {
				float n = hist->GetBinContent(i);
				if (n) std::cout << "Bin \t" << i << "\t LowEdge=\t" << hist->GetXaxis()->GetBinLowEdge(i) << "\t" << n << std::endl; 
			}

		} else {
			std::cout << " ERROR: " << histname << " histogram not found." << std::endl;
			exit (1);
		}

	} else {
		std::cout << " ERROR: " << path << " path not found." << std::endl;
		exit (1);
	}

}

