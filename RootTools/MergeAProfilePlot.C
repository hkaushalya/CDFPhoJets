#include<iostream>
#include<sstream>
#include<TFile.h>
#include<TCanvas.h>
#include<TProfile.h>
#include<TStyle.h>
#include<string>
#include<algorithm>
#include<TLegend.h>
#include<TStyle.h>
#include<TSystem.h>
#include<TDirectory.h>

using namespace std;

////////////////////////////////////////////////////////
//Adds up a Profile hist split among many files.
////////////////////////////////////////////////////////

void MergeAProfilePlot(const int iNfiles, const std::string filenamebase,
		const std::string path, const std::string name, const std::string outfile="",
		const unsigned int rebin=1, const bool logy = false,
		const std::string title="", const std::string xtitle="") 
{
	assert (iNfiles>0 && "Specify number of files to process");
	assert ((filenamebase.length()>0) && "File name required!");
	assert ((path.length()>0) && "path required!");
	assert ((name.length()>0) && "name of the hist required!");
	std::cout << "Reading from " << std::endl;
	std::string data_path(path);
	
	TFile *f=0;
	TProfile* hist_data =0;
	double sum=0;
	
	for (int i = 1; i<= iNfiles; ++i) 
	{
		std::stringstream file;
		file << filenamebase << i;
	
		f = new TFile (file.str().c_str());

		if (f->IsZombie()) 
		{
			std::cout << "ERROR::File " << file.str() << " did not open! " << std::endl;
			exit (1);
		} else {
				//std::cout << "File Added::";
				//f->Print();
		}

			f->cd(data_path.c_str());
			gDirectory->pwd();
			TProfile* temp_data = dynamic_cast<TProfile*> (gDirectory->FindObjectAny(name.c_str()));
			assert (temp_data != NULL && "hist not found!");
			temp_data->Print();

			sum += temp_data->GetEntries();

			if (!hist_data)
			{
					std::string name;
					//name = temp_data->GetName() + std::string("_copy");
					name = temp_data->GetName();
					hist_data = dynamic_cast<TProfile*> (temp_data->Clone(name.c_str()));
					hist_data->SetDirectory(0);
			} else {
				if (temp_data->GetEntries())
				{
					hist_data->Add(temp_data);
				}
				delete f;
			}
			hist_data->Print();
	}
	assert (hist_data != NULL && "HIST_DATA null");

	//TFile *file = new TFile ("Merged.root","RECREATE");
	//hist_data->Write();

	if (sum != hist_data->GetEntries()) std::cout << "sum=" << sum << " hits entries = " << hist_data->GetEntries() << " did not match!" << std::endl;
	
	//file->Print();
	//file->Close();
	//hist_data->RebinAxis(20,"X");
	//hist_data->RebinAxis(20,"Y");
	//hist_data->SetMarkerColor(kBlue);
	if (title.length())
	{
		std::stringstream ti;
		ti << hist_data->GetTitle() << title;
		//hist_data->SetTitle(ti.str().c_str());
		hist_data->SetTitle(title.c_str());
	}
		


	hist_data->Rebin(rebin);
	hist_data->SetLineColor(kBlue);
	TCanvas *c = new TCanvas();
	gPad->SetGridx();
	gPad->SetGridy();
	if (logy) {
		gPad->SetLogy();
		hist_data->SetMinimum(0.005);
	} else hist_data->SetMinimum(0);


	//search for max X
	unsigned lastbin = 0;
	for (unsigned bin = hist_data->GetNbinsX(); bin > 0; --bin)
	{
		if (hist_data->GetBinContent(bin))
		{
			lastbin = bin;
			break;
		}
	}

	if (lastbin > 0 )
	{
		hist_data->GetXaxis()->SetRangeUser(hist_data->GetBinLowEdge(1), hist_data->GetXaxis()->GetBinUpEdge(lastbin + 2));
	}

	
	hist_data->Draw();
	//hist_data->DrawNormalized();


	std::stringstream prntfile;
	//prntfile << name << ".gif";
	//gPad->Print(prntfile.str().c_str());
	std::string outputfile("MergedHists.root");
	if (outfile.length())
	{
		outputfile = outfile;
		TFile out(outputfile.c_str(),"UPDATE");
		hist_data->Write();
		out.ls();
		out.Close();
	
	}
}



void MergeAProfilePlot(const int iNfiles) 
{
	std::string filename("MetTest_PhotonMC.root_");
	//std::string datapath("/Ana/JetFilterV2/Hist/Ana_data/");
	std::string datapath("/Ana/JetFilterV2/Hist/");
	std::string datafile("data.root");
	std::string bgpath("/Ana/JetFilterV2/Hist/Ana_bckg/");
	std::string bgfile("bckg.root");

	/*MergeAProfilePlot(iNfiles, filename.c_str(), datapath.c_str(), "JetPhoEtRatioVsDelPhiJetPho",datafile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), datapath.c_str(), "newJetPhoEtRatioVsDelPhiJetPho",datafile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), datapath.c_str(), "JetMetRatioVsDelPhiJetMet",datafile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), datapath.c_str(), "newJetMetRatioVsDelPhiJetMet",datafile.c_str());
*/
	MergeAProfilePlot(iNfiles, filename.c_str(), datapath.c_str(), "hDbl_JetPhoRatioVsDelPhiJetMet",datafile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), datapath.c_str(), "hDbl_JetMetRatioVsDelPhiJetMet",datafile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), datapath.c_str(), "hDbl_JetPhoEtRatioVsDelPhiJetPho",datafile.c_str());

	
/*	MergeAProfilePlot(iNfiles, filename.c_str(), bgpath.c_str(), "JetPhoEtRatioVsDelPhiJetPho",bgfile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), bgpath.c_str(), "newJetPhoEtRatioVsDelPhiJetPho",bgfile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), bgpath.c_str(), "JetMetRatioVsDelPhiJetMet",bgfile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), bgpath.c_str(), "newJetMetRatioVsDelPhiJetMet",bgfile.c_str());
	*/
	MergeAProfilePlot(iNfiles, filename.c_str(), bgpath.c_str(), "JetPhoRatioVsDelPhiJetMet",bgfile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), bgpath.c_str(), "JetMetRatioVsDelPhiJetMet",bgfile.c_str());
	MergeAProfilePlot(iNfiles, filename.c_str(), bgpath.c_str(), "JetPhoEtRatioVsDelPhiJetPho",bgfile.c_str());
/*	MergeAProfilePlot(iNfiles, filename.c_str(), bgpath.c_str(), "newJetPhoRatioVsDelPhiJetMet",bgfile.c_str());
*/

//	MergeAProfilePlot(iNfiles, filename.c_str(), datapath.c_str(), "newJetPhoRatioVsDelPhiJetMet",datafile.c_str());
//	MergeAProfilePlot(iNfiles, filename.c_str(), bgpath.c_str(), "newJetPhoRatioVsDelPhiJetMet",bgfile.c_str());
}
