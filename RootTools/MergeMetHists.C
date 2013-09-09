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
#include <TPaveText.h>
#include<iomanip>
#include <cmath>
#include "TGaxis.h"
#include "TDirectory.h"
#include "TFolder.h"

using namespace std;

////////////////////////////////////////////////////////
// Created this to add up the histograms to make the
// MET model predictions of MET. This adds up MET All/
// MET after cuts/Met Sig hists.
// added  stuff to get the ratio plots of the above
// and put it side-by-side - 10-18-2008
////////////////////////////////////////////////////////
/*******************************************************
 * 
 * $Id: MergeMetHists.C,v 1.10 2010/02/13 21:55:10 samantha Exp $
 *
 * $Log: MergeMetHists.C,v $
 * Revision 1.10  2010/02/13 21:55:10  samantha
 * ADDED: 1. more stuff. Exclusive case, DiPhoX case.
 *
 *
 * ****************************************************/

//this is to assert systematics calculations are done correctly.
void DebugSystError(const TH1* hist_data,const TH1* hist_bg, const TH1* hist_data_copy, const TH1* hist_err_copy)
{
	//assert all these hists have the same bin widths
	unsigned Nbin1 = hist_data->GetNbinsX(); 
	unsigned Nbin2 = hist_bg->GetNbinsX(); 
	unsigned Nbin3 = hist_data_copy->GetNbinsX(); 
	unsigned Nbin4 = hist_err_copy->GetNbinsX(); 
	assert ( (Nbin1 == Nbin2 && Nbin2 == Nbin3 && Nbin3 == Nbin4 && Nbin4 == Nbin1) && "Number of bins did not match!");
	for (unsigned bin = 1; bin <= (unsigned) hist_data->GetNbinsX(); ++bin)
	{
		float loedge1 = hist_data->GetBinLowEdge(bin);
		float loedge2 = hist_bg->GetBinLowEdge(bin);
		float loedge3 = hist_data_copy->GetBinLowEdge(bin);
		float loedge4 = hist_err_copy->GetBinLowEdge(bin);
		assert ( (loedge1 == loedge2 && loedge2 == loedge3 && loedge3 == loedge4 && loedge4 == loedge1) && "Bin edges did not match!");

		//need just first non zero bin
		if (hist_data->GetBinContent(bin) && hist_bg->GetBinContent(bin))
		{
			std::cout << "DATA : " << hist_data->GetName() << " Integral = " << hist_data->Integral() << std::endl;
			std::cout << "BCKG : " << hist_bg->GetName() << " Integral = " << hist_bg->Integral() << std::endl;
			double diff = hist_data->Integral()-hist_bg->Integral();
			std::cout << "DIFF : " << diff  << std::endl;
			std::cout << std::endl;
			//std::cout << "bin = " << bin << std::endl;
			//std::cout << "data +/-err = " << hist_data->GetBinContent(bin) << ", " << hist_data->GetBinError(bin) << std::endl;
			//std::cout << "bckg +/-err = " << hist_bg->GetBinContent(bin) << ", " << hist_bg->GetBinError(bin) << std::endl;
			//std::cout << "ratio plot data +/-err = " << hist_data_copy->GetBinContent(bin) << ", " << hist_data_copy->GetBinError(bin) << std::endl;
			//std::cout << "ratio plot bckg +/-err = " << hist_err_copy->GetBinContent(bin) << ", " << hist_err_copy->GetBinError(bin) << std::endl;

			break;
		}
	}

}


//return the maximum of the difference
float MaxDiff(const TH1 *def, const TH1 *up, const TH1 *down, const int bin)
{
	//if (bin == def->GetMaximumBin())
/*	if (def->GetBinContent(bin))
	{
		std::cout << " def, up, down, max = " 
					<< def->GetBinContent(bin) << "\t"
					<< up->GetBinContent(bin) << "\t"
					<< down->GetBinContent(bin) << "\t | ";


		std::cout << "max = " << max(def->GetBinContent(bin) - up->GetBinContent(bin),
					def->GetBinContent(bin) - down->GetBinContent(bin)) << std::endl;
					
	}
*/
	return max(def->GetBinContent(bin) - up->GetBinContent(bin),
					def->GetBinContent(bin) - down->GetBinContent(bin));
}



void DoFinalStep(std::vector<TH1*> vHist)
{
//call this function only for MEt plot to get final prediction

		double N_data = vHist.at(0)->Integral();
		double N_qcd  = vHist.at(2)->Integral();
		double N_qcd_stat = sqrt(vHist.at(2)->Integral());
		double N_qcd_syst = 0;
		double dN_syst[10];

		dN_syst[0]=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(3)->Integral()));

			double d1;
			double d2;
			d1=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(4)->Integral()));
			d2=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(5)->Integral()));
			dN_syst[1]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(6)->Integral()));
			d2=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(7)->Integral()));
			dN_syst[2]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(8)->Integral()));
			d2=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(9)->Integral()));
			dN_syst[3]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(10)->Integral()));
			d2=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(11)->Integral()));
			dN_syst[4]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(12)->Integral()));
			d2=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(13)->Integral()));
			dN_syst[5]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(14)->Integral()));
			d2=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(15)->Integral()));
			dN_syst[6]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(16)->Integral()));
			d2=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(17)->Integral()));
			dN_syst[7]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(18)->Integral()));
			d2=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(19)->Integral()));
			dN_syst[8]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(20)->Integral()));
			d2=fabs(1.0*(vHist.at(2)->Integral()-vHist.at(21)->Integral()));
			dN_syst[9]= d1>d2 ? d1 : d2;
			
			for(int i=0; i<10; i++)
			{
				N_qcd_syst=N_qcd_syst+dN_syst[i]*dN_syst[i];
			}
			N_qcd_syst=sqrt(N_qcd_syst);

	std::cout << "DATA = " << N_data << std::endl;
	std::cout << "BCKG = " << N_qcd << "+/-" << N_qcd_stat << "+/-" << N_qcd_syst << std::endl;

}
void DoSystematics(std::vector<TH1*> vHist)
{
/*
	syst1 = def - ue1;
	syst2 = max(def-ue2,def-ue3)
	syst3 = max(def-ue4,def-ue5)
	syst4 = max(def-ue6,def-ue7)
	syst5 = max(def-ue8,def-ue9)
	syst6 = max(def-jer1,def-jer2)
	syst7 = max(def-jer3,def-jer4)
	syst8 = max(def-jer5,def-jer6)
	syst9 = max(def-jer7,def-jer8)
	syst = max(def-jer9,def-jer10)
	
*/
	//do the errors bin by bin
/*	std::cout << setw(6) << "Bin" << setw(6) << "DATA" << setw(6) << "BCKG" 
				<< setw(10) << "Tot.SYS" << setw(8) << "STAT"
				<< setw(8) << "sys1" << setw(8) << "sys2" << setw(8) << "sys3"
				<< setw(8) << "sys4" << setw(8) << "sys5" << setw(8) << "sys6" 
				<< setw(8) << "sys7" << setw(8) << "sys8" << setw(8) << "sys9" 
				<< setw(8) << "sys10" << std::endl;
*/
	for (int bin = 1; bin <= vHist.at(0)->GetNbinsX(); ++bin)
	{

	// I STILL NEED TO DIVIDE BY THE NUMBER OF THE PSEUDOEXPERIMENTS 
	// I ALSO ASSUME THAT SumW2() has been called for all the hists
	// in the EndJob of JetFilterModuleV2.cc!
		
 		
		float stat = sqrt(vHist.at(2)->GetBinContent(bin));
		//std::cout << "bin, content, stat, err = " << bin << "\t" << vHist.at(2)->GetBinContent(bin)
		//				<< "\t" << stat << "\t" << vHist.at(2)->GetBinError(bin) << std::endl; 
		float sys1 = fabs(vHist.at(2)->GetBinContent(bin) - vHist.at(3)->GetBinContent(bin));  //ue1
		float sys2 = fabs(MaxDiff(vHist.at(2), vHist.at(4), vHist.at(5), bin));						//ue2,ue3
		float sys3 = fabs(MaxDiff(vHist.at(2), vHist.at(6), vHist.at(7), bin));						//ue4, ue5
		float sys4 = fabs(MaxDiff(vHist.at(2), vHist.at(8), vHist.at(9), bin));						//ue6, ue7
		float sys5 = fabs(MaxDiff(vHist.at(2), vHist.at(10), vHist.at(11), bin));					//ue8, ue9 
		float sys6 = fabs(MaxDiff(vHist.at(2), vHist.at(12), vHist.at(13), bin));					//jer1, jer2
		float sys7 = fabs(MaxDiff(vHist.at(2), vHist.at(14), vHist.at(15), bin));					//jer3, jer4
		float sys8 = fabs(MaxDiff(vHist.at(2), vHist.at(16), vHist.at(17), bin));					//jer5, jer6
		float sys9 = fabs(MaxDiff(vHist.at(2), vHist.at(18), vHist.at(19), bin));					//jer7, jer8 
		float sys10 = fabs(MaxDiff(vHist.at(2), vHist.at(20), vHist.at(21), bin));					//jer9, jer10 

		float sysTotal = sqrt( pow(stat,2) + pow(sys1,2) + pow(sys2,2) + pow(sys3,2)
										+ pow(sys4,2) + pow(sys5,2) + pow(sys6,2) + pow(sys7,2)
										+ pow(sys8,2) + pow(sys9,2) + pow(sys10,2));

		if (vHist.at(0)->GetBinContent(bin) 
			|| vHist.at(2)->GetBinContent(bin)
			|| sysTotal>0)
		{
			//std::cout << "---------------------Info for bin =" <<  bin << std::endl;
			//std::cout << "\nMax Bin        = " << bin << std::endl;
			//std::cout << "Max Bin stat   = " << stat << std::endl;
			//std::cout << "Max Bin sys1-10= " << sys1 << "\t" << sys2 << "\t" << sys3 << "\t" << sys4 << "\t" << sys5
			//			<< "\t" << sys6 << "\t" << sys7 << "\t" << sys8 << "\t" << sys9 << "\t" << sys10 << std::endl;
			//std::cout << "Max Bin SysTot= " << sysTotal << std::endl;
/*			std::cout << setw(6) << bin << setw(6) << vHist.at(0)->GetBinContent(bin) << setw(6) << vHist.at(2)->GetBinContent(bin)
						<< setw(10) << sysTotal << setw(8) << stat
						<< setw(8) << sys1 << setw(8) << sys2 << setw(8) << sys3
						<< setw(8) << sys4 << setw(8) << sys5 << setw(8) << sys6 
						<< setw(8) << sys7 << setw(8) << sys8 << setw(8) << sys9 
						<< setw(8) << sys10 << std::endl;
*/		}

		//setup the background hist using def hist and sysTotal
		vHist.at(1)->SetBinContent(bin, vHist.at(2)->GetBinContent(bin));
		vHist.at(1)->SetBinError(bin, sysTotal);
	}



	
	//vHist.at(1)->Draw("E");
	
}



// Returns last non-zero bin upper axis value
float FindUpLimit(const TH1* hist, const std::string sAxis="X")
{
	assert (hist != NULL && "hist is NULL!");
	assert ((sAxis=="X" || sAxis=="Y") && "Axis must be X or Y?");

	float uplimit=-10E10;

	for (int i=1; i <= hist->GetNbinsX(); ++i)
	{
		if (sAxis == "X")
		{
			if (hist->GetBinContent(i) > 0)
			{
				uplimit = hist->GetXaxis()->GetBinUpEdge(i);
			}
		} 
		if (sAxis == "Y")
		{
			if (hist->GetBinContent(i) > uplimit) uplimit = hist->GetBinContent(i);
		}
	}
	return uplimit;
}


//void MergeMetHists(const int Nfiles=0, 
//							const std::string filenamebase="",
//							const std::string title="", const std::string name="",
//							const float LoX=-1, const float HiX=-1, const int rebin=1,
//							const std::string printfile="", bool debug=false)
void MergeMetHists(const int Nfiles, 
		const std::string filenamebase,
		const std::string title, const std::string name,
		const float LoX=-1, const float HiX=-1, const int rebin=1,
		const std::string printfile="", const int MetSig=4, const int iPathSelect=0,
		bool debug=false, const bool logy=false, const bool bWriteToFiles = false)
{
	assert (Nfiles>0 && "Number of files must be > 0");
	assert (name.length()>0 && "Require hist name!");
	assert (filenamebase.length()>0 && "Require base of the file name!");

	std::cout << "MetSig=" << MetSig << std::endl;
	std::cout << "Searching for obj = " << name <<std::endl;

	std::string data_path, bg_path, def_path;
	std::string ue1_path, ue2_path, ue3_path, ue4_path, ue5_path, ue6_path;
	std::string ue7_path, ue8_path, ue9_path;
	std::string jer1_path, jer2_path, jer3_path, jer4_path, jer5_path;
	std::string jer6_path, jer7_path, jer8_path, jer9_path, jer10_path;
	
	if (iPathSelect ==0)
	{
		data_path="/Ana/JetFilterV2/Hist/Ana_data/";
		std::cout << "Path Selected " << data_path << std::endl;
		bg_path="/Ana/JetFilterV2/Hist/Ana_bckg/";
		def_path="/Ana/JetFilterV2/Hist/Ana_def/";
		ue1_path="/Ana/JetFilterV2/Hist/Ana_ue1/";
		ue2_path="/Ana/JetFilterV2/Hist/Ana_ue2/";
		ue3_path="/Ana/JetFilterV2/Hist/Ana_ue3/";
		ue4_path="/Ana/JetFilterV2/Hist/Ana_ue4/";
		ue5_path="/Ana/JetFilterV2/Hist/Ana_ue5/";
		ue6_path="/Ana/JetFilterV2/Hist/Ana_ue6/";
		ue7_path="/Ana/JetFilterV2/Hist/Ana_ue7/";
		ue8_path="/Ana/JetFilterV2/Hist/Ana_ue8/";
		ue9_path="/Ana/JetFilterV2/Hist/Ana_ue9/";
		jer1_path="/Ana/JetFilterV2/Hist/Ana_jer1/";
		jer2_path="/Ana/JetFilterV2/Hist/Ana_jer2/";
		jer3_path="/Ana/JetFilterV2/Hist/Ana_jer3/";
		jer4_path="/Ana/JetFilterV2/Hist/Ana_jer4/";
		jer5_path="/Ana/JetFilterV2/Hist/Ana_jer5/";
		jer6_path="/Ana/JetFilterV2/Hist/Ana_jer6/";
		jer7_path="/Ana/JetFilterV2/Hist/Ana_jer7/";
		jer8_path="/Ana/JetFilterV2/Hist/Ana_jer8/";
		jer9_path="/Ana/JetFilterV2/Hist/Ana_jer9/";
		jer10_path="/Ana/JetFilterV2/Hist/Ana_jer10/";
	} else if (iPathSelect == 1)
	{
		data_path="/Ana/JetFilterV2_2/Hist/Ana_data/";
		std::cout << "Path Selected " << data_path << std::endl;
		bg_path="/Ana/JetFilterV2_2/Hist/Ana_bckg/";
		def_path="/Ana/JetFilterV2_2/Hist/Ana_def/";
		ue1_path="/Ana/JetFilterV2_2/Hist/Ana_ue1/";
		ue2_path="/Ana/JetFilterV2_2/Hist/Ana_ue2/";
		ue3_path="/Ana/JetFilterV2_2/Hist/Ana_ue3/";
		ue4_path="/Ana/JetFilterV2_2/Hist/Ana_ue4/";
		ue5_path="/Ana/JetFilterV2_2/Hist/Ana_ue5/";
		ue6_path="/Ana/JetFilterV2_2/Hist/Ana_ue6/";
		ue7_path="/Ana/JetFilterV2_2/Hist/Ana_ue7/";
		ue8_path="/Ana/JetFilterV2_2/Hist/Ana_ue8/";
		ue9_path="/Ana/JetFilterV2_2/Hist/Ana_ue9/";
		jer1_path="/Ana/JetFilterV2_2/Hist/Ana_jer1/";
		jer2_path="/Ana/JetFilterV2_2/Hist/Ana_jer2/";
		jer3_path="/Ana/JetFilterV2_2/Hist/Ana_jer3/";
		jer4_path="/Ana/JetFilterV2_2/Hist/Ana_jer4/";
		jer5_path="/Ana/JetFilterV2_2/Hist/Ana_jer5/";
		jer6_path="/Ana/JetFilterV2_2/Hist/Ana_jer6/";
		jer7_path="/Ana/JetFilterV2_2/Hist/Ana_jer7/";
		jer8_path="/Ana/JetFilterV2_2/Hist/Ana_jer8/";
		jer9_path="/Ana/JetFilterV2_2/Hist/Ana_jer9/";
		jer10_path="/Ana/JetFilterV2_2/Hist/Ana_jer10/";
	} else {
		std::cout << __FUNCTION__ <<": Invalid path selected! exiting !" << std::endl;
		exit (1);
	}

	std::vector<std::string> vPaths;
	vPaths.push_back(data_path);
	vPaths.push_back(bg_path);
	vPaths.push_back(def_path);
	vPaths.push_back(ue1_path);
	vPaths.push_back(ue2_path);
	vPaths.push_back(ue3_path);
	vPaths.push_back(ue4_path);
	vPaths.push_back(ue5_path);
	vPaths.push_back(ue6_path);
	vPaths.push_back(ue7_path);
	vPaths.push_back(ue8_path);
	vPaths.push_back(ue9_path);
	vPaths.push_back(jer1_path);
	vPaths.push_back(jer2_path);
	vPaths.push_back(jer3_path);
	vPaths.push_back(jer4_path);
	vPaths.push_back(jer5_path);
	vPaths.push_back(jer6_path);
	vPaths.push_back(jer7_path);
	vPaths.push_back(jer8_path);
	vPaths.push_back(jer9_path);
	vPaths.push_back(jer10_path);

	TFile *f = 0;
	TH1 *hist_data = 0, *hist_bg = 0;

	int iNHists = vPaths.size();

	std::vector<TH1*> vHist;
	for (int n= 0; n < iNHists; ++n)
	{
		TH1 *temp=0;
		vHist.push_back(temp);
	}

	int NfilesOpened = 0;

	//for (int i=1; i<=1; ++i) 
	//for (int i=1; i<=Nfiles; ++i) 
	for (int i=0; i<Nfiles; ++i) 
	{
		std::stringstream file;
		file << filenamebase << i;

		f = new TFile (file.str().c_str());
		if (f->IsZombie())
		{
			std::cout << "ERROR::File " << file.str() << " did not open! Exiting." << std::endl;
			exit (1);
		} else {
			NfilesOpened++;
			if (debug)
			{
				std::cout << "File Added::";
				f->Print();
			}
		}

		for (unsigned int iPath = 0; iPath < vPaths.size(); ++iPath)
		{
			//std::cout << "iPath = " << iPath << std::endl;
			//f->cd();
			//gDirectory->pwd();
			f->cd(vPaths.at(iPath).c_str());
			//gDirectory->pwd();
			//if (iPath<2) f->ls();
			TH1 *hTemp = dynamic_cast<TH1*> (gDirectory->FindObjectAny(name.c_str()));

			if (hTemp->GetEntries())	// this has to be done to avoid crashes when adding hists which some how have 'sum=nan' instead of 'sum=0' when they do not have any entries.
			{
				if (! vHist.at(iPath))
				{
					//std::string localname = hTemp->GetName() + std::string ("_Clone");
					vHist.at(iPath) = dynamic_cast<TH1*>(hTemp->Clone (name.c_str()));
					assert(vHist.at(iPath) != NULL &&  "ERROR! hist cast failed");
					vHist.at(iPath)->SetDirectory(0);
				} else
				{
					vHist.at(iPath)->Add(hTemp);
				}
			}
		}

		delete f;

	}

//	assert(vHist.size() == vPaths.size());
//	for (int k=0; k < vHist.size(); ++k)
//	{
//		vHist.at(k)->Print();
//	}

	DoSystematics(vHist);
	if (name == "Met") DoFinalStep(vHist);
	
	hist_data = vHist.at(0);
	hist_bg = vHist.at(1);

	// I added this so I can comapre different background predictions
	// when I need, like the double smeared background and the singly smeared
	if (bWriteToFiles)
	{
		TFile f("DATA_HISTS.root","UPDATE");
		if (f.FindObjectAny(hist_data->GetName()) != NULL)
		{
			std::cout << "Warnning! Hist with name "<< hist_data->GetName() << " already exists in file " <<f.GetName() << std::endl;
		}
		hist_data->Write();
		std::cout << "Wrote DATA hist " << hist_data->GetName() << " to file " << f.GetName() << std::endl;
		f.Close();
		TFile g("BCKG_HISTS.root","UPDATE");
		if (g.FindObjectAny(hist_bg->GetName()) != NULL)
		{
			std::cout << "Warnning! Hist with name "<< hist_bg->GetName() << " already exists in file " << g.GetName() << std::endl;
		}
		hist_bg->Write();
		std::cout << "Wrote BCKG hist " << hist_bg->GetName() << " to file " << g.GetName() << std::endl;
		g.Close();
	}


/*	std::cout << "NORMALIZING BG TO DATA " << std::endl;
	double data_int = hist_data->Integral();
	double bg_int = hist_bg->Integral();
	hist_bg->Scale(data_int/bg_int);
	std::cout << "SCALE = " << data_int/bg_int << std::endl;
*/

	if (debug) 
	{
		hist_data->Print("all");
		hist_bg->Print("all");
	}

	std::cout << "Total file added = " << NfilesOpened << std::endl;
	gStyle->SetOptStat("");

	if (hist_data->GetEntries()) 	hist_data->Rebin(rebin);
	if (hist_bg->GetEntries()) hist_bg->Rebin(rebin);

	TH1 *hist_err_copy = NULL;
	std::string bgname = hist_bg->GetName() + std::string ("err_copy");
	hist_err_copy = dynamic_cast<TH1*>(hist_bg->Clone (bgname.c_str()));
	hist_err_copy->SetDirectory(0);

	TH1 *hist_data_copy = NULL;
	std::string dataname = hist_data->GetName() + std::string ("data_copy");
	hist_data_copy = dynamic_cast<TH1*>(hist_data->Clone (dataname.c_str()));
	hist_data_copy->SetDirectory(0);



	float x_loLim, x_hiLim;
	if (LoX >0) x_loLim = LoX;
	else x_loLim = hist_data->GetBinLowEdge(1);

	if (HiX >0) x_hiLim = HiX;
	else x_hiLim = max(FindUpLimit(hist_data), FindUpLimit(hist_bg)) + hist_data->GetBinWidth(1) * 2;
	if (debug)
	{
		std::cout << "min, max = " << x_loLim << ", " << x_hiLim << std::endl;
	}

	float y_hiLim = max(FindUpLimit(hist_data,"Y"), FindUpLimit(hist_bg,"Y"));
	if (logy) y_hiLim *= 10;
	else y_hiLim += y_hiLim * 0.1;



	gStyle->SetCanvasColor (10);
	gStyle->SetCanvasBorderSize (0);
	gStyle->SetCanvasBorderMode (0);
	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (0);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);
	gStyle->SetCanvasDefW(1200);
	gStyle->SetCanvasDefH(600);
	gStyle->SetAxisColor(8,"XY");
	
	int labelfont = 10 * 4 + 2;		//10 * font ID + precision (2 = scalable)
	int titlefont = 10 * 4 + 2;		//10 * font ID + precision (2 = scalable)
	gStyle->SetLabelFont(labelfont,"X");
	gStyle->SetLabelFont(labelfont,"Y");
	gStyle->SetTitleFont(titlefont,"X");
	gStyle->SetTitleFont(titlefont,"Y");
	gStyle->SetLabelSize(0.04,"X");
	gStyle->SetLabelSize(0.027,"Y");
	//gStyle->SetLabelOffset(0.9);
	gStyle->SetTitleSize(0.03,"Y");
	gStyle->SetTitleOffset(1.7,"Y");
	//TGaxis::SetMaxDigits(3);



	hist_data->UseCurrentStyle();
	hist_bg->UseCurrentStyle();

	TCanvas *c1= new TCanvas;
	c1->Divide(2,1);
	c1->cd(1);
	if (logy)
	{
		if (hist_data->GetEntries() > 0 && hist_bg->GetEntries() > 0) gPad->SetLogy();
	}
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();

	TPaveText *tp = new TPaveText(0.02,0.92,0.98,0.99,"NDC");
	tp->SetLineColor(10);
	tp->SetTextFont(titlefont);


	std::string tt(hist_data->GetTitle());
	hist_data->SetTitle("");
	hist_bg->SetTitle("");
	if (title.length()>0)
	{
		//tt += " - ";
		tt += title;
		//tt = title;
		if (debug)
		{
			std::cout << tt << std::endl;
		}
		std::cout << "title = " << title << std::endl;
		tp->AddText(tt.c_str());
		//tp->AddText("#Delta#phi(#slash{E}_{T}-Jet1) after cuts: gg-signal Cor #slash{E}_{T}-Sig>0");
	}

	std::stringstream ytitle, xtitle;
	ytitle << "Events / " << setprecision(3) << hist_data->GetXaxis()->GetBinWidth(1) << " GeV";
	if (debug)
	{
		std::cout << hist_data->GetBinWidth(1) <<std::endl;
	}
	if (name == "MetAll") xtitle << "#slash{E}_{T} (for all events) (GeV)";
	else if (name == "Met") xtitle << "#slash{E}_{T} (after cuts) (GeV)";
	else if (name == "MetSig") xtitle << "#slash{E}_{T} Significance (for all events)";
	else if (name == "Njet15") xtitle << "Njets^{E_{T}>15GeV} (After #slash{E}_{T}-Sig cut)";
	else if (name == "Njet20") xtitle << "Njets^{E_{T}>20GeV} (After #slash{E}_{T}-Sig cut)";
	else if (name == "Njet25") xtitle << "Njets^{E_{T}>25GeV} (After #slash{E}_{T}-Sig cut)";
	else if (name == "Njet30") xtitle << "Njets^{E_{T}>30GeV} (After #slash{E}_{T}-Sig cut)";
	else if (name == "Njet35") xtitle << "Njets^{E_{T}>35GeV} (After #slash{E}_{T}-Sig cut)";

	if (debug)
	{
		std::cout << "xtitle=" << xtitle.str() << std::endl;
	}

	if (debug)
	{
		hist_data->Print();
		std::cout << "bin#\tLoEdge\tdata\tdata_err\tbg_cont\tbg_err"<<std::endl;
		for (unsigned bin = 0; bin <= (unsigned) hist_data_copy->GetNbinsX() + 1; ++ bin)
		{
			float val = hist_data->GetBinContent (bin);
			float err = hist_data->GetBinError(bin);
			float val2 = hist_bg->GetBinContent (bin);
			float err2 = hist_bg->GetBinError(bin);
			float loEdge = hist_data->GetBinLowEdge(bin);

			if (val>0 || err>0 || val2>0 || err2>0)
				std::cout << bin << "\t" << loEdge <<"\t" << val << "\t" << err << "\t\t" << val2 << "\t" << err2 << std::endl;
		}
	}

	hist_data->GetXaxis()->SetTitle(xtitle.str().c_str());
	if (name.find("Njet") == std::string::npos) hist_data->GetYaxis()->SetTitle(ytitle.str().c_str());
	hist_data->GetXaxis()->CenterTitle(true);
	hist_data->GetYaxis()->CenterTitle(true);
	hist_bg->GetXaxis()->SetTitle(xtitle.str().c_str());
	if (name.find("Njet") == std::string::npos) hist_bg->GetYaxis()->SetTitle(ytitle.str().c_str());
	hist_bg->GetXaxis()->CenterTitle(true);
	hist_bg->GetYaxis()->CenterTitle(true);

	//temp
	x_loLim = 0; x_hiLim = 200;
	hist_data->GetXaxis()->SetRangeUser(x_loLim, x_hiLim);
	hist_bg->GetXaxis()->SetRangeUser(x_loLim, x_hiLim);

	hist_data->SetMinimum(0.1);
	hist_bg->SetMinimum(0.1);
	hist_data->SetMaximum(y_hiLim);
	hist_bg->SetMaximum(y_hiLim);


	hist_data->SetMarkerStyle(8);
	hist_data->SetMarkerSize(1.0);
	hist_data->SetMarkerColor(kBlue);
	hist_data->SetLineColor(kBlue);
	hist_bg->SetMarkerColor(kRed);
	hist_bg->SetLineColor(kRed);
	hist_bg->SetFillColor(kRed);
	hist_data->SetTitleSize(0.04);
	hist_bg->SetTitleSize(0.04);
	//hist_bg->SetLineWidth(2);
	//hist_data->SetLineWidth(2);
	TH1* hist_bg_copy = dynamic_cast<TH1*>(hist_bg->Clone ("hist_bg_BlkLine"));
	hist_bg_copy->SetLineColor(kBlack);
	hist_bg->SetLineColor(kBlack);
	hist_bg_copy->SetLineWidth(2);
	hist_bg_copy->SetFillColor(0);
	//hist_bg_copy->Draw("L");
	//hist_bg->Draw("sameE2");
//	hist_bg->SetFillColor(10);
//	hist_bg->SetLineColor(10);
	hist_bg->Draw("E2");
	hist_bg_copy->Draw("SAME HIST");
	hist_data->Draw("sameP");
	tp->Draw();

	TLegend *leg = new TLegend (0.5,0.8,0.90,0.90);
	leg->SetTextFont(42);
	leg->SetTextSize(0.025);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);

	//std::stringstream leg_data, leg_bg;
	//leg_data << "Data (E=" << hist_data->GetEntries() << " M=" << hist_data->GetMean()
	//		<< " R=" << hist_data->GetRMS() << ")";
	//leg_bg << "Bkg (E=" << hist_bg->GetEntries() << " M=" << hist_bg->GetMean()
	//		<< " R=" << hist_bg->GetRMS() << ")";

	//leg->AddEntry(hist_data,leg_data.str().c_str());
	//leg->AddEntry(hist_bg,leg_bg.str().c_str());
	//leg->AddEntry(hist_data,"Data (Measured) (DET Jets) ");
	//leg->AddEntry(hist_bg, "MC Prediction (HAD Jets, Norm to Data)");
	leg->AddEntry(hist_data,"PYTHIA #gamma MC");
	//leg->AddEntry(hist_data,"PYTHIA DiPho MC");
	leg->AddEntry(hist_bg, "Fake ME_{T}");
	leg->Draw();

	// now to make the ratio plots
	for (unsigned bin = 0; bin <= (unsigned) hist_data_copy->GetNbinsX() + 1; ++ bin)
	{
		const float val = hist_err_copy->GetBinContent (bin);
		const float scale = val ? 1. / val : 0;
		hist_data_copy->SetBinContent (bin, (hist_data_copy->GetBinContent (bin) - val) * scale);
		hist_data_copy->SetBinError (bin, hist_data_copy->GetBinError (bin) * scale);
	};
	for (unsigned bin = 0; bin <= (unsigned) hist_err_copy->GetNbinsX() + 1; ++ bin)
	{
		float value = hist_err_copy->GetBinContent (bin);
		float error = hist_err_copy->GetBinError (bin);
		hist_err_copy->SetBinError (bin, value ? error / value : 0);
		hist_err_copy->SetBinContent (bin, 0);
	};


	/*
		TH1 *hist_ratio = NULL;
		std::string myname = hist_data->GetName() + std::string ("_copy");
		hist_ratio = dynamic_cast<TH1*>(hist_data->Clone (myname.c_str()));
		hist_ratio->Divide(hist_bg);
	 */


	hist_data_copy->UseCurrentStyle();
	hist_err_copy->UseCurrentStyle();
	//new TCanvas();
	c1->cd(2);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	hist_data_copy->SetTitle("");
	hist_err_copy->SetTitle("");
	//hist_data_copy->SetTitle(tt.c_str());
	//hist_err_copy->SetTitle(tt.c_str());
	hist_data_copy->GetXaxis()->SetTitle(xtitle.str().c_str());
	hist_err_copy->GetXaxis()->SetTitle(xtitle.str().c_str());
	hist_data_copy->GetXaxis()->SetRangeUser(x_loLim, x_hiLim);
	hist_err_copy->GetXaxis()->SetRangeUser(x_loLim, x_hiLim);
	float fRatioHist_ymax = 1.0;
	float fRatioHist_ymin = -1.0;
	hist_data_copy->SetMinimum(fRatioHist_ymin);
	hist_data_copy->SetMaximum(fRatioHist_ymax);
	hist_err_copy->SetMinimum(fRatioHist_ymin);
	hist_err_copy->SetMaximum(fRatioHist_ymax);
	std::stringstream ratio_ytitle;
	ratio_ytitle << "#gamma MC/Fake ME_{T} - 1";
	hist_data_copy->GetYaxis()->SetTitle(ratio_ytitle.str().c_str());
	hist_err_copy->GetYaxis()->SetTitle(ratio_ytitle.str().c_str());
	//hist_data_copy->SetTitle(ratio_ytitle.str().c_str());
	//hist_err_copy->SetTitle(ratio_ytitle.str().c_str());

	hist_data_copy->GetXaxis()->CenterTitle(true);
	hist_data_copy->GetYaxis()->CenterTitle(true);
	hist_err_copy->GetXaxis()->CenterTitle(true);
	hist_err_copy->GetYaxis()->CenterTitle(true);
	//	hist_data->GetYaxis()->SetRangeUser(;
	////	hist_bg->SetMinimum(0.1);


	hist_data_copy->SetLineColor(kBlue);
	hist_data_copy->SetMarkerColor(kBlue);
	hist_data_copy->SetMarkerStyle (8);
	hist_data_copy->SetMarkerSize(1.0);
	hist_err_copy->SetFillColor(kRed);
	hist_err_copy->SetFillStyle(3002);
	hist_data_copy->SetLineWidth(2);
	hist_data_copy->Draw("P");	
	hist_err_copy->Draw("same E2");	
	tp->Draw();


	c1->cd();
	if (printfile.length()>0)
	{
		std::stringstream epsfile, giffile;
		epsfile << printfile << ".eps";
		giffile << printfile << ".gif";
		c1->Print(giffile.str().c_str());
		c1->Print(epsfile.str().c_str());
	}

	DebugSystError(hist_data,hist_bg, hist_data_copy, hist_err_copy);
	

}


void MergeMetHists(const int metsig, const int iNfiles, const int Incl, const float jetcut, 
		const int iPathSelect = 0, const int debug=0)
{
	assert (iNfiles>0 && "Number of files must be > 0");
	std::stringstream title1, title2, title3, file1, file2, file3, file;
	std::stringstream title4, title5, title6, title7, title8;
	std::stringstream file4, file5, file6, file7, file8,file9;


	if (Incl == 1)  //incl njet15
	{
		std::cout << "Using njet>=" << jetcut << "GeV settings." << std::endl;
		//title1 << "PYTHIA #gamma^{Det#Eta<0.95}_{signal}+>=1 jet^{E_{T}>="<< jetcut << "GeV, Det#Eta<3.0} : #slash{E}_{T}-Sig>" << metsig << ", Dbl.Smear+JER offset";
	 //title2 << "PYTHIA #gamma^{Det#Eta<0.95}_{signal}+>=1 jet^{E_{T}>="<< jetcut << "GeV, Det#Eta<3.0}: Dbl.Smear+JER offset";
		title1 << "#gamma^{DetEta<0.95}_{E_{T}>30GeV}+>=1 jet_{E_{T}>="<< jetcut << "GeV}^{Det#Eta<3.0} : Raw#slash{E}_{T}-Sig>" << metsig;
		title2 << "#gamma^{DetEta<0.95}_{E_{T}>30GeV}+>=1 jet_{E_{T}>="<< jetcut << "GeV}^{Det#Eta<3.0} : Raw#slash{E}_{T}-Sig>" << metsig;
		file  << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_NJet15";
		file1 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_MetAll";
		file2 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_Met";
		file3 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_MetSig";
		file4 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_EtJet1";
		file5 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_PhoEt";
		file6 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_DelPhiPhoMet";
		file7 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_DelPhiJetMet";
		file8 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_Ht";
		file9 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_EtJet2";

		title4 << " :Njets>=1 Raw#slash{E}_{T}-Sig>"<< metsig;
	}

	if (Incl == 0)  //excl njet15
	{
		std::cout << "Using njet==15GeV settings." << std::endl;
		title1 << "#gamma^{DetEta<0.95}_{signal}+=1 jet^{E_{T}>="<< jetcut<< "GeV, DetEta<3.0} : Pseudo Expts=1, #slash{E}_{T}-Sig>" << metsig;
		title2 << "#gamma^{DetEta<0.95}_{signal}+=1 jet^{E_{T}>="<< jetcut << "GeV, DetEta<3.0} : Pseudo Expts=1";
		file  << "PhoMCSignal_ExclNJet" << jetcut<<"_MetSig"<< metsig << "_NJet15";
		file1 << "PhoMCSignal_ExclNJet" << jetcut<<"_MetSig"<< metsig << "_MetAll";
		file2 << "PhoMCSignal_ExclNJet" << jetcut<<"_MetSig"<< metsig << "_Met";
		file3 << "PhoMCSignal_ExclNJet" << jetcut<<"_MetSig"<< metsig << "_MetSig";
		file4 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_EtJet1";
		file5 << "PhoMCSignal_ExclNJet" << jetcut<<"_MetSig"<< metsig << "_PhoEt";
		file6 << "PhoMCSignal_ExclNJet" << jetcut<<"_MetSig"<< metsig << "_DelPhiPhoMet";
		file7 << "PhoMCSignal_ExclNJet" << jetcut<<"_MetSig"<< metsig << "_DelPhiJetMet";
		file8 << "PhoMCSignal_ExclNJet" << jetcut<<"_MetSig"<< metsig << "_Ht";
		file9 << "PhoMCSignal_InclNJet" << jetcut<<"_MetSig"<< metsig << "_EtJet2";

		title4 << " :Njets==1 #slash{E}_{T}-Sig>"<< metsig;
	}

	if (Incl == 2)  //no jet cut
	{
		std::cout << "Using njet>=" << jetcut << "GeV settings." << std::endl;
		title1 << "#gamma^{DetEta<0.95}_{E_{T}>30GeV}+ NO Jet cut : Pseudo Expts=1, Raw #slash{E}_{T}-Sig>" << metsig;
		title2 << "#gamma^{DetEta<0.95}_{E_{T}>30GeV}+ NO jet cut : Pseudo Expts=1";
		file  << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_NJet15";
		file1 << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_MetAll";
		file2 << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_Met";
		file3 << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_MetSig";
		file4 << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_EtJet1";
		file5 << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_PhoEt";
		file6 << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_DelPhiPhoMet";
		file7 << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_DelPhiJetMet";
		file8 << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_Ht";
		file9 << "PhoMCSignal_noJetCut" << jetcut<<"_MetSig"<< metsig << "_EtJet2";

		title4 << " :NO jet cut, Raw #slash{E}_{T}-Sig>"<< metsig;
	}

if (Incl == 3)  //diPhoXt no cut
	{
		std::cout << "Using njet>=" << jetcut << "GeV settings." << std::endl;
		title1 << "DiPhoX: #gamma^{DetEta<0.95}_{E_{T}>13GeV}+ NO Jet cut : gg-signal Corr #slash{E}_{T}-Sig>" << metsig;
		title2 << "DiPhoX: #gamma^{DetEta<0.95}_{E_{T}>13GeV}+ NO jet cut";
		file  << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_NJet15";
		file1 << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_MetAll";
		file2 << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_Met";
		file3 << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_MetSig";
		file4 << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_EtJet1";
		file5 << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_PhoEt";
		file6 << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_DelPhiPhoMet";
		file7 << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_DelPhiJetMet";
		file8 << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_Ht";
		file9 << "DiPhoX" << jetcut<<"_MetSig"<< metsig << "_EtJet2";

		title4 << " :NO jet cut, gg-signal Corr #slash{E}_{T}-Sig>"<< metsig;
	}

	std::string sRootFileBase("MetTest_PhotonMC.root_");
	//std::string sRootFileBase("DiPhoX.root_");

	if (metsig == 0)
	{
		MergeMetHists (iNfiles,sRootFileBase,title2.str(), "MetSig",-1,-1,4, file3.str(), metsig, debug,1,1);
		MergeMetHists (iNfiles,sRootFileBase,title1.str(), "MetAll",-1,-1,4, file1.str(), metsig, debug,1,1);
	}

	if (metsig != 0)
	{
		MergeMetHists (iNfiles, sRootFileBase, title1.str(), "Met",-1,-1,4, file2.str(), metsig, iPathSelect,debug,1,1);
	}
	MergeMetHists (iNfiles, sRootFileBase, title4.str(), "Etjet",-1,-1,4, file4.str(), metsig, iPathSelect, debug,1,1);
	MergeMetHists (iNfiles, sRootFileBase, title4.str(), "Etjet2",-1,-1,4, file9.str(), metsig, iPathSelect, debug);
	MergeMetHists (iNfiles, sRootFileBase, title1.str(), "Njet15",-1,-1,1, file.str(), metsig, iPathSelect, debug);
	MergeMetHists (iNfiles, sRootFileBase, title4.str(), "dPhi3",-1,-1,4, file7.str(), metsig, iPathSelect, debug);
	MergeMetHists (iNfiles, sRootFileBase, title4.str(), "dPhi1",-1,-1,4, file6.str(), metsig, iPathSelect, debug);
	MergeMetHists (iNfiles, sRootFileBase, title4.str(), "Ht",-1,-1,4, file8.str(), metsig, iPathSelect, debug);
	MergeMetHists (iNfiles, sRootFileBase, title4.str(), "Et1",-1,-1,4, file5.str(), metsig, iPathSelect, debug);
	//MergeMetHists (iNfiles, sRootFileBase, title1.str(), "Et_a[0]",-1,-1,5, file1.str(), metsig, iPathSelect, debug);
	//MergeMetHists (iNfiles, sRootFileBase, title1.str(), "Et_a[1]",-1,-1,5, file2.str(), metsig, iPathSelect, debug);

}

void MergeMetHists(const int iNfiles, const int debug = 0)
{
	//	MergeMetHists(0,iNfiles, debug);
	//MergeMetHists(3,iNfiles, debug);
	//MergeMetHists(4,iNfiles, debug);
	//MergeMetHists(5,iNfiles, debug);
	//MergeMetHists(6,iNfiles, debug);
	//MergeMetHists(7,iNfiles, debug);
}
