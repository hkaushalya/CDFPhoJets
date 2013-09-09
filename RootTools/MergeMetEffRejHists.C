#include<iostream>
#include<sstream>
#include<TFile.h>
#include<TCanvas.h>
#include<TH1.h>
#include<TStyle.h>
#include<string>
#include <algorithm>
#include <TLegend.h>
#include<TStyle.h>
#include<TSystem.h>
#include<TDirectory.h>

using namespace std;

////////////////////////////////////////////////////////
// Created this to add up the histograms to make the
// MET model Rejection Power and Efficiency . This adds 
// up MetAll (Met before MetSig cut) and Met (Met after 
// met sig cut) for data and background.
// And make the final plots
////////////////////////////////////////////////////////

void MergeMetEffRejHists(const int iNfiles=0, const bool data = 0, const std::string sOutFileName="") 
{
	//data = 0 == background prediction
	//data = 1 == data
	
	if (iNfiles <=0)
	{ 
		std::cout << "Specify number of files to process!" << std::endl;
		return;
	}
	if (sOutFileName.length() == 0)
	{
		std::cout << "Output file name required!" << std::endl;
		return;
	}

	std::string path;
	if (data)
	{
		path = "/Ana/JetFilterV2/Hist/Ana_data/";
	} else {
		path = "/Ana/JetFilterV2/Hist/Ana_bckg/";
	}
	std::cout << "Hist search path = " << path << std::endl;
	
	TFile *f;
	TH1* hist_MetAll;
	TH1* hist_MetAfter;
	double sum=0;
	double sum0=0;
	
	for (int i = 1; i<= iNfiles; ++i) 
	{
		std::stringstream file;
		file << "MetTest_PhotonMC.root_"<< i;
	
		f = new TFile (file.str().c_str());

		if (f->IsZombie()) 
		{
			std::cout << "ERROR::File " << file.str() << "did not open! " << std::endl;
		} else {
				//std::cout << "File Added::";
				//f->Print();
		}

			std::string name = "Met";
			std::string name1 = "MetAll";

			f->cd(path.c_str());
			gDirectory->pwd();
			//gDirectory->ls();

			TH1* htemp_Met = (TH1*) gDirectory->FindObjectAny(name.c_str());
			TH1* htemp_MetAll = (TH1*) gDirectory->FindObjectAny(name1.c_str());

			assert (htemp_Met != NULL && "Met object not found!");
			assert (htemp_MetAll != NULL && "MetAll object not found!");

			htemp_Met->Print();
			htemp_MetAll->Print();

			sum += htemp_Met->GetEntries();
			sum0 += htemp_MetAll->GetEntries();

			if (i==1)
			{
				hist_MetAll = (TH1*) htemp_Met->Clone();
				hist_MetAfter = (TH1*) htemp_MetAll->Clone();
				hist_MetAll->SetDirectory(0);
				hist_MetAfter->SetDirectory(0);
			} else {
				hist_MetAll->Add(htemp_Met);
				hist_MetAfter->Add(htemp_MetAll);
			}

			delete f;
	}
	assert (hist_MetAll != NULL && "hist_MetAll null");
	assert (hist_MetAfter != NULL && "HIST_DATA null");

	TFile *outfile = new TFile (sOutFileName.c_str(),"RECREATE");
	hist_MetAll->Write();
	hist_MetAfter->Write();


	if (sum != hist_MetAll->GetEntries()) std::cout << "sum=" << sum << " hits entries = " << hist_MetAll->GetEntries() << " did not match!" << std::endl;
	if (sum0 != hist_MetAfter->GetEntries()) std::cout << "sum0=" << sum0 << " hits entries = " << hist_MetAfter->GetEntries() << " did not match!" << std::endl;

	std::cout << "Met/MetAll Entries = " << sum << " / " << sum0 << std::endl;
	
	outfile->Close();
	std::cout << "Wrote hist/s to file: ";
	outfile->Print();

}
