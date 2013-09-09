#include<iostream>
#include<sstream>
#include<TFile.h>
#include<TCanvas.h>
#include<TH2.h>
#include<TStyle.h>
#include<string>
#include<algorithm>
#include<TLegend.h>
#include<TStyle.h>
#include<TSystem.h>
#include<TDirectory.h>

using namespace std;

////////////////////////////////////////////////////////
//Adds up a hist 2D split among many files.
////////////////////////////////////////////////////////

void MergeA2DPlot(const std::vector<std::string> file, const std::vector<std::string> hist, const std::vector<std::string> path, const std::string outputfile="MergedHists.root")
{
	assert (file.size()>0 && "Specify number of files to process");
	assert (hist.size()>0 && "File name required!");
	assert ((path.size()>0) && "path required!");
	assert ((hist.size() == path.size()) && "number of hists and paths do not match!");
	
	
	TFile out(outputfile.c_str(),"UPDATE");
	for (unsigned int i=0; i<hist.size(); ++i)
	{
		TH2* hist_data =0;
		double sum=0;
		std::cout << "Merging hist " << hist.at(i)  << " path = " << path.at(i) << std::endl;
		
		for (unsigned int j = 0; j< file.size(); ++j) 
		{
			TFile f(file.at(j).c_str());

			if (f.IsZombie()) 
			{
				std::cout << "ERROR::File " << file.at(i) << " did not open! " << std::endl;
			}

			f.cd(path.at(i).c_str());
			gDirectory->pwd();
			TH2* temp_data = (TH2*) gDirectory->FindObjectAny(hist.at(i).c_str());
			if (temp_data == NULL)
			{
				std::cout << "Hist "<< hist.at(i) << " not found in file " << file.at(j) << std::endl;
				exit (1);
			}
			assert (temp_data != NULL && "hist not found!");
			//temp_data->Print();
			//new TCanvas();

			sum += temp_data->GetEntries();

			if (!hist_data)
			{
				std::string name;
				//name = temp_data->GetName() + std::string("_copy");
				name = temp_data->GetName();
				hist_data = dynamic_cast<TH2*> (temp_data->Clone(name.c_str()));
				hist_data->SetDirectory(0);
			} else {
				if (temp_data->GetEntries())
				{
					hist_data->Add(temp_data);
				}
			}
		}
		assert (hist_data != NULL && "HIST_DATA null");
		out.cd();
		hist_data->Write();

		if (sum != hist_data->GetEntries()) std::cout << "sum=" << sum << " hits entries = " << hist_data->GetEntries() << " did not match!" << std::endl;

		//TCanvas *c = new TCanvas("c1","Merged Histogram",200,10,700,500);
		new TCanvas();
		hist_data->Draw();

	}

	out.Close();

}
