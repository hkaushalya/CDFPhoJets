#include<TFile.h>
#include<TCanvas.h>
#include<iostream>
#include<sstream>
#include<TH1.h>
#include<vector>
#include<cmath>
#include <exception>
#include <stdexcept> // out_of_range exception

//This is to spit out and debug my final MergeMetHists.C which makes the
//log/ratio plots pf Met/MetSig/MetAll.
//Somehow the systematics are very small. Infact they seems to be <1% at low values in all hits. . MetModel precision is in on the order of ~20% according to Sasha.
// This will given path, hist and a bin number will spit out
//data, the central value or def and the systematic hist values
// so I can see what is going on.

TH1* GetHist(const std::vector<std::string>& file, const std::string path, const std::string histname, const int debug)
{
	TH1 *hist = 0;
	//std::cout << "hist=" << hist << std::endl;
	
	//first add the hists
	for (unsigned i=0; i < file.size(); ++i)
	{
		TFile f(file.at(i).c_str());
		if (f.IsZombie())
		{
			std::cout << "File " << file.at(i) << " not found!" << std::endl;
			exit (1);
		} else 
		{
			f.cd(path.c_str());
		 	if (debug)
			{
				f.Print();
				gDirectory->pwd();
				std::cout << "Searching for " << histname << std::endl;
			}
			
			if (!hist)
			{
				TH1 *temp_hist = (TH1*) gDirectory->FindObjectAny(histname.c_str());
				//if (debug && temp_hist->GetEntries()) temp_hist->Print();
				if (debug) temp_hist->Print();
				assert (temp_hist != NULL && "temp_hist is null!");
				std::string name = temp_hist->GetName() + std::string ("Clone");
				hist = (TH1*) temp_hist->Clone(name.c_str());
				assert (hist != NULL && "temp_hist is null!");
				hist->SetDirectory(0);
			} else
			{
				TH1* temp_hist = (TH1*) gDirectory->FindObjectAny(histname.c_str());
				//if (debug && temp_hist->GetEntries()) temp_hist->Print();
				if (debug) temp_hist->Print();
				assert (temp_hist != NULL && "temp_hist is null!");
				if (debug) std::cout << temp_hist->GetTitle() << " = " << temp_hist->GetEntries() << std::endl;
				if (temp_hist->GetEntries())
				{
					hist->Add(temp_hist);
				}
			} // if (hist)

		} // if (f.IsZombie)
	} //for 

	assert(hist != NULL && "hist is null");
	return hist;
}


//return the diff between two hist 
double HistoBinDiff(const TH1* h, const TH1* h1, const TH1* h2,int bin_ind)
{
	assert(h != NULL && "HistoBinDiff:: h is null!");
	assert(h1 != NULL && "HistoBinDiff:: h1 is null!");
	assert(h2 != NULL && "HistoBinDiff:: h2 is null!");
	double value=0.0;
	double v0=h->GetBinContent(bin_ind);
	double err=h->GetBinError(bin_ind);
	double v1=h1->GetBinContent(bin_ind);
	double v2=h2->GetBinContent(bin_ind);
	
	value=(fabs(v0-v1) > fabs(v0-v2)) ? sqrt(err*err+(v0-v1)*(v0-v1)) : sqrt(err*err+(v0-v2)*(v0-v2));
	return value;
}


//dump systematic errors for a bin
void DumpSystErr(const TH1 *data_hist, const TH1 *bckg_hist, const TH1 *def_hist, const std::vector<TH1*> syst_ue_hist, const std::vector<TH1*> syst_jer_hist, const int bin)
{

	assert (data_hist != NULL && "data hist null!");
	assert (bckg_hist != NULL && "bckg hist null!");
	assert (def_hist != NULL && "def hist null!");
	for (unsigned i=0; i < syst_ue_hist.size(); ++i)
	{
		assert (syst_ue_hist.at(i) != NULL && "DumpSystErr:: syst_ue_hist is null!");
	}
	for (unsigned i=0; i < syst_jer_hist.size(); ++i)
	{
		assert (syst_jer_hist.at(i) != NULL && "DumpSystErr:: syst_jer_hist is null!");
	}


	
	if (bin>=0 && bin <= data_hist->GetNbinsX())
	{
		std::cout << data_hist->GetName() << " for bin = " << bin <<std::endl;
		std::cout << "data,err\t= " << data_hist->GetBinContent(bin) << ", " << data_hist->GetBinError(bin) << std::endl;
		std::cout << "def,err \t= " << def_hist->GetBinContent(bin) << ", " << def_hist->GetBinError(bin) << std::endl;
		std::cout << "bckg,err\t= " << bckg_hist->GetBinContent(bin) << ", " << bckg_hist->GetBinError(bin) << std::endl;
		for (unsigned i=0; i< syst_ue_hist.size(); ++i)
		{
			std::cout << "ue" << i <<" ,err\t= " << syst_ue_hist.at(i)->GetBinContent(bin) << ", " << syst_ue_hist.at(i)->GetBinError(bin) << std::endl;
		}
		for (unsigned i=0; i< syst_jer_hist.size(); ++i)
		{
			std::cout << "jer" << i <<" ,err\t= " << syst_jer_hist.at(i)->GetBinContent(bin) << ", " << syst_jer_hist.at(i)->GetBinError(bin) << std::endl;
		}

		double ueerr_sq=0, ueerr=0;
		double jererr_sq=0, jererr=0;
		try
		{
		double sys_ue1d = HistoBinDiff(bckg_hist, syst_ue_hist.at(0), def_hist, bin);
		double sys_ue23 = HistoBinDiff(bckg_hist, syst_ue_hist.at(1), syst_ue_hist.at(2), bin);
		double sys_ue45 = HistoBinDiff(bckg_hist, syst_ue_hist.at(3), syst_ue_hist.at(4), bin);
		double sys_ue67 = HistoBinDiff(bckg_hist, syst_ue_hist.at(5), syst_ue_hist.at(6), bin);
		double sys_ue89 = HistoBinDiff(bckg_hist, syst_ue_hist.at(7), syst_ue_hist.at(8), bin);
		
		double sys_jer12 = HistoBinDiff(bckg_hist, syst_jer_hist.at(0), syst_jer_hist.at(1), bin);
		double sys_jer34 = HistoBinDiff(bckg_hist, syst_jer_hist.at(2), syst_jer_hist.at(3), bin);
		double sys_jer56 = HistoBinDiff(bckg_hist, syst_jer_hist.at(4), syst_jer_hist.at(5), bin);
		double sys_jer78 = HistoBinDiff(bckg_hist, syst_jer_hist.at(6), syst_jer_hist.at(7), bin);
		double sys_jer910= HistoBinDiff(bckg_hist, syst_jer_hist.at(8), syst_jer_hist.at(9), bin);
		
		double sys_def = def_hist->GetBinError(bin);
		std::cout << "sys_ue1d   = " << sys_ue1d << std::endl; 
		std::cout << "sys_ue23   = " << sys_ue23 << std::endl; 
		std::cout << "sys_ue45   = " << sys_ue45 << std::endl; 
		std::cout << "sys_ue67   = " << sys_ue67 << std::endl; 
		std::cout << "sys_ue89   = " << sys_ue89 << std::endl; 
		std::cout << "sys_jer12  = " << sys_jer12 << std::endl; 
		std::cout << "sys_jer34  = " << sys_jer34 << std::endl; 
		std::cout << "sys_jer56  = " << sys_jer56 << std::endl; 
		std::cout << "sys_jer78  = " << sys_jer78 << std::endl; 
		std::cout << "sys_jer910 = " << sys_jer910 << std::endl; 

		
		ueerr_sq = pow(sys_ue1d,2) + pow(sys_ue23,2) + pow(sys_ue45,2) + pow(sys_ue67,2) + pow(sys_ue89,2);
		jererr_sq = pow(sys_jer12,2) + pow(sys_jer34,2) + pow(sys_jer56,2) + pow(sys_jer78,2) + pow(sys_jer910,2);
		
		double syst_tot = sqrt(ueerr_sq + jererr_sq + pow(sys_def,2));

		std::cout << "total syst = " << syst_tot << std::endl;
		} catch (std::out_of_range& e)
		{
			std::cout << __LINE__ << "::HistoBinDiff::" << e.what() << std::endl;
		}
	} else
	{
		std::cout << "Bin number requested is out of range!" << std::endl;
	}

}

void DebugMetModelSyst(const std::vector<std::string> file, const std::string histname, const int bin, const int rebin, const int debug, TFile *rootFile = NULL)
{

	assert (file.size() >0 && "No files specified.");
	assert (histname.length() > 0 && "hist name not specified");
	assert (bin >=0 && "Bin must be >=0"); //do i need the underflow bin and may be overflow bin?

//JetFilterFolders
//JetClu-0.4
//MatchPhoJet0.4
// Ana_data
// Ana_bckg
// Ana_def
// Ana_ue1,Ana_ue2,Ana_ue3,Ana_ue4,Ana_ue5,Ana_ue6,Ana_ue7,Ana_ue8,Ana_ue9
// Ana_jer1,Ana_jer2,Ana_jer3,Ana_jer4,Ana_jer5,Ana_jer6,Ana_jer7,Ana_jer8,Ana_jer9,Ana_jer10

	std::string data_path("/Ana/JetFilterV2/Hist/Ana_data");
	std::string bckg_path("/Ana/JetFilterV2/Hist/Ana_bckg");
	std::string def_path("/Ana/JetFilterV2/Hist/Ana_def");
	std::vector<std::string> syst_ue_path, syst_jer_path;

	TH1 *data_hist;
	TH1 *bckg_hist;
	TH1 *def_hist;
	std::vector<TH1*> syst_ue_hist;
	std::vector<TH1*> syst_jer_hist;


	data_hist = GetHist(file,data_path, histname, debug);
	std::stringstream n1; n1<< data_hist->GetTitle() << "(Ana_data)";
	data_hist->SetTitle(n1.str().c_str());
	data_hist->SetName(n1.str().c_str());
	new TCanvas("data hist");
	data_hist->Draw();
	
	bckg_hist = GetHist(file,bckg_path, histname, debug);
	std::stringstream n2; n2<< bckg_hist->GetTitle() << "(Ana_bckg)";
	bckg_hist->SetTitle(n2.str().c_str());
	bckg_hist->SetName(n2.str().c_str());
	new TCanvas("background hist");
	bckg_hist->Draw();
	
	def_hist = GetHist(file,def_path, histname, debug);
	std::stringstream n3; n3<< def_hist->GetTitle() << "(Ana_def)";
	def_hist->SetTitle(n3.str().c_str());
	def_hist->SetName(n3.str().c_str());
	new TCanvas("default background hist");
	def_hist->Draw();

	for (int i=0; i<9; ++i)
	{
		std::stringstream path;
		path << "/Ana/JetFilterV2/Hist/Ana_ue"<< (i+1);
		syst_ue_path.push_back(path.str());
		TH1* temp = GetHist(file, path.str(),histname, debug);
		assert ( temp != NULL && "ue temp hist is null!");
		syst_ue_hist.push_back(temp);
		std::stringstream n; n<< syst_ue_hist.at(i)->GetTitle() << "(Ana_ue" << (i+1) << ")";
		syst_ue_hist.at(i)->SetTitle(n.str().c_str());
		syst_ue_hist.at(i)->SetName(n.str().c_str());
		if (i==0)
		{
			new TCanvas(path.str().c_str());
			syst_ue_hist.at(i)->Draw();
		}
	}
	for (int i=0; i<10; ++i)
	{
		std::stringstream path;
		path << "/Ana/JetFilterV2/Hist/Ana_jer"<< (i+1);
		syst_jer_path.push_back(path.str());
		TH1* temp = GetHist(file, path.str(), histname, debug);
		assert ( temp != NULL && "jer temp hist is null!");
		syst_jer_hist.push_back(temp);
		std::stringstream n; n<< syst_jer_hist.at(i)->GetTitle() << "(Ana_jer" << (i+1) << ")";
		syst_jer_hist.at(i)->SetTitle(n.str().c_str());
		syst_jer_hist.at(i)->SetName(n.str().c_str());
		if (i==0)
		{
			new TCanvas(path.str().c_str());
			syst_jer_hist.at(i)->Draw();
		}
	}
	
	if (rebin>1)
	{
		data_hist->Rebin(rebin);
		def_hist->Rebin(rebin);
		bckg_hist->Rebin(rebin);
		for (unsigned i=0; i<syst_ue_hist.size(); ++i)
		{
			try {
				syst_ue_hist.at(i)->Rebin(rebin);
			} catch (std::out_of_range& e)
			{
				std::cout << __LINE__ << "::ue rebin hist i = " << i << "::" << e.what() << std::endl;
			}
		}
		for (unsigned i=0; i<syst_jer_hist.size(); ++i)
		{
			syst_jer_hist.at(i)->Rebin(rebin);
		}	
	}

	assert((syst_ue_hist.size() == 9) && " syst_ue size != 9"); 
	for (int i=0; i<9; ++i)
	{
		assert ((syst_ue_hist.at(i) != NULL) && "syst_ue null assertion failed");
	}
	assert((syst_jer_hist.size() == 10) && " syst_ue size != 10"); 
	for (int i=0; i<10; ++i)
	{
		assert( (syst_jer_hist.at(i) != NULL) && "syst_jer null assertion failed");;
	}


	if (rootFile->IsOpen())
	{
		//rootFile->Write(); //this does not work as I have used SetDirectory(0) for all hists
		//I must bring the rootFile to focus. otherwise the gDirectory points to Rint:/ wihtin this  function scope.
		//gDirectory->pwd();
		rootFile->cd();
		data_hist->Write();
		def_hist->Write();
		bckg_hist->Write();
		for (int i=0; i<8; ++i)
		{
			syst_ue_hist.at(i)->Write();
			if (debug) syst_ue_hist.at(i)->Print();
		}
		for (int i=0; i<10; ++i)
		{
			syst_jer_hist.at(i)->Write();
			if (debug) syst_jer_hist.at(i)->Print();
		}
		
		std::cout << "Hists written to ";
		rootFile->Print();
		if (debug)  rootFile->ls();
	}


	DumpSystErr(data_hist,bckg_hist,def_hist,syst_ue_hist,syst_jer_hist,bin);
		
	
}

void DebugMetModelSyst(const std::string filebase = "MetTest_PhotonMC.root_", const int Nfiles=100, const std::string histname="Njet15", const int bin=2, const int rebin=1, std::string rootfile="HistFiles.root", const int debug = false)
{
	std::vector<std::string> files;

	for (int i=1; i<= Nfiles; ++i)
	{
		std::stringstream filename;
		filename << filebase << i;
		files.push_back(filename.str());
	}
	TFile *f = new TFile(rootfile.c_str(),"RECREATE");

	if (f->IsZombie())
	{
		std::cout << "TFile " << rootfile << " open failed! exiting!" << std::endl; 
	}
	
	DebugMetModelSyst(files, histname, bin, rebin, debug, f);
	f->Close();

}

void DebugMetModelSyst(const int Nfiles, const int bin)
{
	DebugMetModelSyst("MetTest_PhotonMC.root_", Nfiles, "Njet30", bin);
}
