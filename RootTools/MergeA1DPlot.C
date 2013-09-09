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
//Adds up a hist split among many files.
////////////////////////////////////////////////////////

void MergeA1DPlot(const std::vector<std::string> file, 
						const std::vector<std::string> hist, 
						const std::vector<std::string> path, 
						const std::string outputfile="MergedHists.root" , 
						const bool bDrawFinalHist=false)
{
	assert (file.size()>0 && "Specify number of files to process");
	assert (hist.size()>0 && "File name required!");
	assert ((path.size()>0) && "path required!");
	
	
	TFile out(outputfile.c_str(),"UPDATE");
	for (unsigned int i=0; i<hist.size(); ++i)
	{
		std::cout << "Merging the hist " << hist.at(i);
		TH1* hist_data =0;
		double sum=0;

		for (unsigned int j = 0; j< file.size(); ++j) 
		{
			TFile f(file.at(j).c_str());

			if (f.IsZombie()) 
			{
				std::cout << "ERROR::File " << file.at(i) << " did not open! " 
					<< std::endl;
			}

			f.cd(path.at(i).c_str());
			//f.ls();
			TH1* temp_data = (TH1*) gDirectory->FindObjectAny(hist.at(i).c_str());
			if (temp_data == NULL)
			{
				std::cout << "Hist "<< hist.at(i) << " not found in file " 
					<< file.at(j) << std::endl;
				continue;
			} else temp_data->Print();
			assert (temp_data != NULL && "hist not found!");

			sum += temp_data->GetEntries();

			if (!hist_data)
			{
				std::string name;
				name = temp_data->GetName();
				hist_data = dynamic_cast<TH1*> (temp_data->Clone());
				hist_data->SetDirectory(0);
			} else {
				if (temp_data->GetEntries())
				{
					hist_data->Add(temp_data);
				}
			}
		}
		assert (hist_data != NULL && "HIST_DATA null");
		if (sum != hist_data->GetEntries()) std::cout << "sum=" << sum 
			<< " hits entries = " << hist_data->GetEntries() 
			<< " did not match!" << std::endl;
		out.cd();
		hist_data->Write();
		std::cout << " written to output file " << outputfile << std::endl; 

		if (bDrawFinalHist) {
			new TCanvas();
			hist_data->Draw("E");
		}

	}

	out.Close();

}



TH1* MergeA1DPlot(const std::vector<std::string> file, 
						const std::string path, const std::string hist_name,
						TDirectory* dir)
{
	assert (file.size()>0 && "Specify number of files to process");
	assert (hist_name.length()>0 && "File name required!");
	assert ((path.length()>0) && "path required!");
	
	
		TH1* hist_data =0;
		double sum=0;

		for (unsigned int j = 0; j< file.size(); ++j) 
		{
			TFile f(file.at(j).c_str());

			if (f.IsZombie()) 
			{
				std::cout << "ERROR::File " << file.at(j) << " did not open! " 
					<< std::endl;
			}

			f.cd(path.c_str());
			TH1* temp_data = (TH1*) gDirectory->FindObjectAny(hist_name.c_str());
			if (temp_data == NULL)
			{
				std::cout << "Hist "<< hist_name << " not found in file " 
					<< file.at(j) << std::endl;
				continue;
			}
			assert (temp_data != NULL && "hist not found!");

			sum += temp_data->GetEntries();

			if (!hist_data)
			{
				std::string name;
				name = temp_data->GetName();
				//hist_data = dynamic_cast<TH1*> (temp_data->Clone(name.c_str()));
				//outfile.cd(path.c_str());
				hist_data = dynamic_cast<TH1*> (temp_data->Clone());
				//hist_data->SetDirectory(gDirectory);
				dir->Add(hist_data);
				dir->pwd();
			} else {
				if (temp_data->GetEntries())
				{
					hist_data->Add(temp_data);
				}
			}
		}
		assert (hist_data != NULL && "HIST_DATA null");
		hist_data->Print();
		if (sum != hist_data->GetEntries()) std::cout << "sum=" << sum 
				<< " hits entries = " << hist_data->GetEntries() 
				<< " did not match!" << std::endl;


	assert (hist_data != NULL);

//	outfile.cd(path.c_str());
//	outfile.ls();
//	gDirectory->Add(hist_data);
//	outfile.ls();
	
	return hist_data;


}
