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
#include"MakeRatioHist.h"


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


float GetMaximum(const TH1* hist1, const TH1* hist2, const int bin)
{
	assert (hist1 != NULL && "hist1 is null!");
	assert (hist2 != NULL && "hist2 is null!");
	float max1=0, max2=0;
	if (bin > 0 && bin < hist1->GetNbinsX()) max1= hist1->GetBinContent(bin);
	if (bin > 0 && bin < hist2->GetNbinsX()) max2= hist2->GetBinContent(bin);

	return (max1>max2)? max1: max2;
}


void MergeCalibMetHists(const int iNfiles=0) 
{
	assert (iNfiles>0 && "Specify number of files to process");
	gSystem->CompileMacro("~/samantha/RootTools/MakeRatioHist.C","=");
	std::string data_path("/Ana/JetFilterV2/Hist/METJetClu0.4_scen3/");
	
	TFile *f = 0;
	TH1* hist_mc = 0;
	TH1* hist_mc1 = 0;
	TH1* hist_mc2 = 0;
	TH1* hist_mc3 = 0;
	TH1* hist_mc4 = 0;
	TH1* hist_mc5 = 0;
	TH1* hist_data = 0;
	TH1* hist_data1 = 0;
	TH1* hist_data2 = 0;
	TH1* hist_data3 = 0;
	TH1* hist_data4 = 0;
	TH1* hist_data5 = 0;
	double sum=0;
	double sum0=0;
	double sum1=0;
	double sum2=0;
	double sum3=0;
	double sum4=0;
	double sumdata=0;
	double sumdata1=0;
	double sumdata2=0;
	double sumdata3=0;
	double sumdata4=0;
	double sumdata5=0;
	
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

			std::string name = "fMetSigCalib_estimate";
			std::string name0 = "fMetSigCalib_estimate_njet0";
			std::string name1 = "fMetSigCalib_estimate_njet1";
			std::string name2 = "fMetSigCalib_estimate_njet2";
			std::string name3 = "fMetSigCalib_estimate_njet3";
			std::string name4 = "fMetSigCalib_estimate_njet4";
			
			std::string name_ = "fMetSig_estimate";
			std::string name0_ = "fMetSig_estimate_njet0";
			std::string name1_ = "fMetSig_estimate_njet1";
			std::string name2_ = "fMetSig_estimate_njet2";
			std::string name3_ = "fMetSig_estimate_njet3";
			std::string name4_ = "fMetSig_estimate_njet4";
	
			//std::string name5 = "";
			/*std::cout << "looking for " << name << std::endl;
			std::cout << "looking for " << name0 << std::endl;
			std::cout << "looking for " << name1 << std::endl;
			std::cout << "looking for " << name2 << std::endl;
			std::cout << "looking for " << name3 << std::endl;
			std::cout << "looking for " << name4 << std::endl;
			*/
			f->cd(data_path.c_str());
			gDirectory->pwd();

			TH1* temp_mc = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name.c_str()));
			TH1* temp_mc1 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name0.c_str()));
			TH1* temp_mc2 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name1.c_str()));
			TH1* temp_mc3 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name2.c_str()));
			TH1* temp_mc4 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name3.c_str()));
			TH1* temp_mc5 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name4.c_str()));
			TH1* temp_data = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name_.c_str()));
			TH1* temp_data1 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name0_.c_str()));
			TH1* temp_data2 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name1_.c_str()));
			TH1* temp_data3 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name2_.c_str()));
			TH1* temp_data4 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name3_.c_str()));
			TH1* temp_data5 = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name4_.c_str()));
			temp_data2->Print();
			temp_mc2->Print();

			assert (temp_mc != NULL && "0 object not found!");
			assert (temp_mc1 != NULL && "1 object not found!");
			assert (temp_mc2 != NULL && "2 object not found!");
			assert (temp_mc3 != NULL && "3 object not found!");
			assert (temp_mc4 != NULL && "4 object not found!");
			assert (temp_mc5 != NULL && "5 object not found!");
			assert (temp_data != NULL && "data object not found!");
			assert (temp_data1 != NULL && "data1 object not found!");
			assert (temp_data2 != NULL && "data2 object not found!");
			assert (temp_data3 != NULL && "data3 object not found!");
			assert (temp_data4 != NULL && "data4 object not found!");
			assert (temp_data5 != NULL && "data5 object not found!");

			sum += temp_mc->GetEntries();
			sum0 += temp_mc1->GetEntries();
			sum1 += temp_mc2->GetEntries();
			sum2 += temp_mc3->GetEntries();
			sum3 += temp_mc4->GetEntries();
			sum4 += temp_mc5->GetEntries();
			sumdata += temp_data->GetEntries();
			sumdata1 += temp_data1->GetEntries();
			sumdata2 += temp_data2->GetEntries();
			sumdata3 += temp_data3->GetEntries();
			sumdata4 += temp_data4->GetEntries();
			sumdata5 += temp_data5->GetEntries();
			

			if (i==1)
			{
					hist_mc = dynamic_cast<TH1*> (temp_mc->Clone());
					hist_mc1 = dynamic_cast<TH1*> (temp_mc1->Clone());
					hist_mc2 = dynamic_cast<TH1*> (temp_mc2->Clone());
					hist_mc3 = dynamic_cast<TH1*> (temp_mc3->Clone());
					hist_mc4 = dynamic_cast<TH1*> (temp_mc4->Clone());
					hist_mc5 = dynamic_cast<TH1*> (temp_mc5->Clone());
					hist_data = dynamic_cast<TH1*> (temp_data->Clone());
					hist_data1 = dynamic_cast<TH1*> (temp_data1->Clone());
					hist_data2 = dynamic_cast<TH1*> (temp_data2->Clone());
					hist_data3 = dynamic_cast<TH1*> (temp_data3->Clone());
					hist_data4 = dynamic_cast<TH1*> (temp_data4->Clone());
					hist_data5 = dynamic_cast<TH1*> (temp_data5->Clone());
					hist_mc->SetDirectory(0);
					hist_mc1->SetDirectory(0);
					hist_mc2->SetDirectory(0);
					hist_mc3->SetDirectory(0);
					hist_mc4->SetDirectory(0);
					hist_mc5->SetDirectory(0);
					hist_data->SetDirectory(0);
					hist_data1->SetDirectory(0);
					hist_data2->SetDirectory(0);
					hist_data3->SetDirectory(0);
					hist_data4->SetDirectory(0);
					hist_data5->SetDirectory(0);
			} else {
				hist_mc->Add(temp_mc);
				hist_mc1->Add(temp_mc1);
				hist_mc2->Add(temp_mc2);
				hist_mc3->Add(temp_mc3);
				hist_mc4->Add(temp_mc4);
				hist_mc5->Add(temp_mc5);
				hist_data->Add(temp_data);
				hist_data1->Add(temp_data1);
				hist_data2->Add(temp_data2);
				hist_data3->Add(temp_data3);
				hist_data4->Add(temp_data4);
				hist_data5->Add(temp_data5);
				delete f;
			}
	}
	
	assert (hist_mc != NULL && "hist_mc null");
	assert (hist_mc1 != NULL && "hist_mc null");
	assert (hist_mc2 != NULL && "hist_mc null");
	assert (hist_mc3 != NULL && "hist_mc null");
	assert (hist_mc4 != NULL && "hist_mc null");
	assert (hist_mc5 != NULL && "hist_mc null");
	assert (hist_data != NULL && "HIST_DATA null");
	assert (hist_data1 != NULL && "HIST_DATA null");
	assert (hist_data2 != NULL && "HIST_DATA null");
	assert (hist_data3 != NULL && "HIST_DATA null");
	assert (hist_data4 != NULL && "HIST_DATA null");
	assert (hist_data5 != NULL && "HIST_DATA null");

	//TFile *file = new TFile ("Merged_MetSigCalibHists.root","RECREATE");
	//TFile *file = new TFile ("Merged.root","RECREATE");
	TFile *file = new TFile ("Merged.root","UPDATE");
	hist_mc->Write();
	hist_mc1->Write();
	hist_mc2->Write();
	hist_mc3->Write();
	hist_mc4->Write();
	hist_mc5->Write();
	hist_data->Write();
	hist_data1->Write();
	hist_data2->Write();
	hist_data3->Write();
	hist_data4->Write();
	hist_data5->Write();

	if (sum != hist_mc->GetEntries()) std::cout << "sum=" << sum << " hits entries = " << hist_mc->GetEntries() << " did not match!" << std::endl;
	
	if (sum0 != hist_mc1->GetEntries()) std::cout << "sum0=" << sum0 << " hits entries = " << hist_mc1->GetEntries() << " did not match!" << std::endl;
	if (sum1 != hist_mc2->GetEntries()) std::cout << "sum1=" << sum1 << " hits entries = " << hist_mc2->GetEntries() << " did not match!" << std::endl;
	if (sum2 != hist_mc3->GetEntries()) std::cout << "sum2=" << sum2 << " hits entries = " << hist_mc3->GetEntries() << " did not match!" << std::endl;
	if (sum3 != hist_mc4->GetEntries()) std::cout << "sum3=" << sum3 << " hits entries = " << hist_mc4->GetEntries() << " did not match!" << std::endl;
	if (sum4 != hist_mc5->GetEntries()) std::cout << "sum4=" << sum4 << " hits entries = " << hist_mc5->GetEntries() << " did not match!" << std::endl;
	file->Close();
	std::cout << "Wrote hist/s to file: ";
	file->Print();

	unsigned rebin = 20;

	TLegend *leg = new TLegend (0.5,0.72,0.41,0.9);
	std::cout << "leg = " << leg << std::endl;
	gStyle->SetOptStat(1111111);
/*
new TCanvas();
	gPad->SetLogy();
	hist_mc->Rebin(rebin);
	//leg->AddEntry(hist_mc,"DATA");
	hist_mc->SetLineColor(kBlue);
	hist_mc->DrawNormalized();
	hist_data->Rebin(rebin);
	hist_data->DrawNormalized("same");
	//leg->AddEntry(hist_data,"MC");
	//leg->Draw();

*/

//hist_mc1->Print();
//hist_data1->Print();
/*
	new TCanvas();
	gPad->SetLogy();
	hist_mc1->Sumw2();
	hist_mc1->Rebin(rebin);
	hist_mc1->SetLineColor(kBlue);
	hist_mc1->DrawNormalized();
	hist_data1->Sumw2();
	hist_data1->Rebin(rebin);
	hist_data1->DrawNormalized("same");
//	leg->AddEntry(hist_mc1,"DATA");
//	leg->AddEntry(hist_data1,"MC");
//	leg->Draw();
	
*/
/*
	std::cout << "hist_data2 = " << hist_data2 << std::endl;
	std::cout << "hist_mc2 = " << hist_mc2 << std::endl;
	hist_data2->Print();
	hist_mc2->Print();

	hist_mc2->Rebin(rebin);
	hist_data2->Rebin(rebin);
	double intgl_mc2 = hist_mc2->Integral();
	double intgl_data2 = hist_data2->Integral();
	hist_mc2->Scale(1.0/intgl_mc2);
	hist_data2->Scale(1.0/intgl_data2);
	hist_data2->SetLineColor(kBlue);
	leg->AddEntry(hist_data2,"DATA");;
	leg->AddEntry(hist_mc2,"MC");;
	//new TCanvas();
//	hist_data2->SetMinimum(0);
//:	hist_data2->SetMaximum(GetMaximum(hist_data2, hist_mc2, 1)+0.1);
	//hist_data2->Draw();
	//hist_mc2->Draw("same");


	for (int bin=1; bin <= hist_data2->GetNbinsX(); ++bin)
	{
			std::cout << "cont data= " << hist_data2->GetBinContent(bin) << ", " << hist_data2->GetBinError(bin) << std::endl;
			std::cout << "cont bg  = " << hist_mc2->GetBinContent(bin) << ", " << hist_mc2->GetBinError(bin) <<  std::endl;
			float val = hist_mc2->GetBinContent(bin);
			float err = hist_mc2->GetBinError(bin);
			float scale = val ? 1.0/val : 0;
			std::cout << "scale = 1/bg cont = " << scale << std::endl;
			hist_data2->SetBinContent(bin,(hist_data2->GetBinContent(bin) - val) * scale);
			hist_data2->SetBinError(bin, hist_data2->GetBinError(bin) * scale);

			hist_mc2->SetBinError(bin, val ? err * scale : 0);
			hist_mc2->SetBinContent(bin,0);

			std::cout << "after cont data= " << hist_data2->GetBinContent(bin) << ", " << hist_data2->GetBinError(bin) << std::endl;
			std::cout << "after cont bg  = " << hist_mc2->GetBinContent(bin) << ", " << hist_mc2->GetBinError(bin) << std::endl;
			break;

	}

	
//	new TCanvas();
//	hist_data2->Draw("P");
//	hist_mc2->Draw("sameE");


return;
*/

/*
new TCanvas();
	gPad->SetLogy();
	hist_mc3->Rebin(rebin);
	hist_mc3->SetLineColor(kBlue);
	hist_mc3->DrawNormalized();
	hist_data3->Rebin(rebin);
	hist_data3->DrawNormalized("same");
	//leg->Draw();
gPad->Print("Njet2.gif");
	new TCanvas();
	gPad->SetLogy();
	hist_mc4->Rebin(rebin);
	hist_mc4->SetLineColor(kBlue);
	hist_mc4->DrawNormalized();
	hist_data4->Rebin(rebin);
	hist_data4->DrawNormalized("same");
//	leg->Draw();
gPad->Print("Njet3.gif");
	new TCanvas();
	gPad->SetLogy();
	hist_mc5->Rebin(rebin);
	hist_mc5->SetLineColor(kBlue);
	hist_mc5->DrawNormalized();
	hist_data5->Rebin(rebin);
	hist_data5->DrawNormalized("same");
//	leg->Draw();
gPad->Print("Njet4.gif");
*/
}
