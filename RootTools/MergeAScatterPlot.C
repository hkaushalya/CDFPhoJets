#include<iostream>
#include<sstream>
#include<TFile.h>
#include<TCanvas.h>
#include<TH2.h>
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
// MET model predictions of MET. This adds up All METSig
// Calibration hists and write to a root file. Then Run
// MetSigFitter on it to make the final METSig Calibration
// plots. This DOES NOT DRAW ANY HISTOGRAMS.
////////////////////////////////////////////////////////
//fMetSigCalib_estimate_njet0
//			if (j==1) name = "fMetSigCalib_estimate_njet0";
//			//fMetSigCalib_estimate_njet1
//			if (j==2) name = "fMetSigCalib_estimate_njet1";
//			if (j==3) name = "fMetSigCalib_estimate_njet2";
//			if (j==4) name = "fMetSigCalib_estimate_njet3";
//			if (j==5) name = "fMetSigCalib_estimate_njet4";



void MergeAScatterPlot(const int iNfiles=0) 
{
	assert (iNfiles>0 && "Specify number of files to process");
	std::cout << "Reading from " << std::endl;
	gSystem->pwd();
	std::cout << "" << std::endl;
	std::string data_path("/Ana/JetFilterV2/Hist/Ana_data/");
	
	TFile *f;
	TH2* hist_data;
	double sum=0;
	
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

			std::string name = "MetAllVsMetSig";
			f->cd(data_path.c_str());
			gDirectory->pwd();

			TH2* temp_data = (TH2*) gDirectory->FindObjectAny(name.c_str());

			assert (temp_data != NULL && "0 object not found!");

			sum += temp_data->GetEntries();

			if (i==1)
			{
					hist_data = (TH2*) temp_data->Clone("copy");
					hist_data->SetDirectory(0);
			} else {
				hist_data->Add(temp_data);
				delete f;
			}
	}
	assert (hist_data != NULL && "HIST_DATA null");

	TFile *file = new TFile ("Merged.root","RECREATE");
	hist_data->Write();

	if (sum != hist_data->GetEntries()) std::cout << "sum=" << sum << " hits entries = " << hist_data->GetEntries() << " did not match!" << std::endl;
	
	file->Print();
	file->Close();
	hist_data->RebinAxis(20,"X");
	hist_data->RebinAxis(20,"Y");
	hist_data->SetMarkerColor(kBlue);

	hist_data->Draw();

}
