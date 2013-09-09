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
#include<TIterator.h>
#include<TKey.h>

using namespace std;

////////////////////////////////////////////////////////
// Try to replicate the dir structure of a root file
////////////////////////////////////////////////////////


void GetDirs(const std::string path, std::vector<std::string> vDir)
{
	size_t found;
	size_t length = path.length();
	std::vector<size_t> vPts;

	size_t begin = 0, end = 0;
/*
	int i =0;
	while (i<5)
	{
		end = path.find_first_of("/");
		
		found = path.find("/",pos);
		
		if (found == std::string::npos) break;
		else {
			std::cout << "post = " << found << std::endl;
			std::string subs = path.substr(pos,found);
			std::cout << "substr = " << subs << std::endl;
			pos = found+1;
		}
		vPts.push_back(found);
		i++;
	}

*/
	

}

//void MakeARootFile(const std::vector<std::vector> vPaths, TFile &file)
//void MakeARootFile(TFile &file)
void MakeARootFile()
{

//	fordd (unsigned int iPath=0; iPath < vPath.size(); iPath++)
//	{
//		std::vector<std::string> vDir;
//		GetDirs(vPath.at(iPath), vDir);
		
//	}

	std::vector<std::string> vDir;
	TFile *outfile = new TFile ("newRoot.root","RECREATE");
	
	if (outfile->IsZombie()) { std::cout << "File not created! " << std::endl; assert (false); }

	outfile->mkdir("Hist");
	gDirectory->cd("Hist");
	gDirectory->mkdir("CENTRAL");
	gDirectory->cd("CENTRAL");
	gDirectory->mkdir("1Jet");
	gDirectory->cd("1Jet");
	gDirectory->mkdir("Event");
	gDirectory->mkdir("Photon");
	gDirectory->mkdir("LeadJet");
	gDirectory->mkdir("PhotonLeadJet");
	outfile->cd();
	outfile->cd("Hist/CENTRAL");
	gDirectory->mkdir("2Jet");
	gDirectory->cd("2Jet");
	gDirectory->mkdir("Event");
	gDirectory->mkdir("Photon");
	gDirectory->mkdir("LeadJet");
	gDirectory->mkdir("SecondLeadJet");
	gDirectory->mkdir("PhotonLeadJet");
	gDirectory->mkdir("Photon2ndLeadJet");
	gDirectory->mkdir("Photon2Jets");
	gDirectory->mkdir("2Jets");

	outfile->cd();
	outfile->cd("Hist");
	gDirectory->mkdir("EMJESUP");
	gDirectory->cd("EMJESUP");
	gDirectory->mkdir("1Jet");
	gDirectory->cd("1Jet");
	gDirectory->mkdir("Event");
	gDirectory->mkdir("Photon");
	gDirectory->mkdir("LeadJet");
	gDirectory->mkdir("PhotonLeadJet");
	outfile->cd();
	outfile->cd("Hist/EMJESUP");
	gDirectory->mkdir("2Jet");
	gDirectory->cd("2Jet");
	gDirectory->mkdir("Event");
	gDirectory->mkdir("Photon");
	gDirectory->mkdir("LeadJet");
	gDirectory->mkdir("SecondLeadJet");
	gDirectory->mkdir("PhotonLeadJet");
	gDirectory->mkdir("Photon2ndLeadJet");
	gDirectory->mkdir("Photon2Jets");
	gDirectory->mkdir("2Jets");



	outfile->cd();
	outfile->cd("Hist");
	gDirectory->mkdir("EMJESDOWN");
	gDirectory->cd("EMJESDOWN");
	gDirectory->mkdir("1Jet");
	gDirectory->cd("1Jet");
	gDirectory->mkdir("Event");
	gDirectory->mkdir("Photon");
	gDirectory->mkdir("LeadJet");
	gDirectory->mkdir("PhotonLeadJet");
	outfile->cd();
	outfile->cd("Hist/EMJESDOWN");
	gDirectory->mkdir("2Jet");
	gDirectory->cd("2Jet");
	gDirectory->mkdir("Event");
	gDirectory->mkdir("Photon");
	gDirectory->mkdir("LeadJet");
	gDirectory->mkdir("SecondLeadJet");
	gDirectory->mkdir("PhotonLeadJet");
	gDirectory->mkdir("Photon2ndLeadJet");
	gDirectory->mkdir("Photon2Jets");
	gDirectory->mkdir("2Jets");

	outfile->cd();
	outfile->cd("Hist");
	gDirectory->mkdir("SIDEBAND");
	gDirectory->cd("SIDEBAND");
	gDirectory->mkdir("1Jet");
	gDirectory->cd("1Jet");
	gDirectory->mkdir("Event");
	gDirectory->mkdir("Photon");
	gDirectory->mkdir("LeadJet");
	gDirectory->mkdir("PhotonLeadJet");
	outfile->cd();
	outfile->cd("Hist/SIDEBAND");
	gDirectory->mkdir("2Jet");
	gDirectory->cd("2Jet");
	gDirectory->mkdir("Event");
	gDirectory->mkdir("Photon");
	gDirectory->mkdir("LeadJet");
	gDirectory->mkdir("SecondLeadJet");
	gDirectory->mkdir("PhotonLeadJet");
	gDirectory->mkdir("Photon2ndLeadJet");
	gDirectory->mkdir("Photon2Jets");
	gDirectory->mkdir("2Jets");


	std::vector<std::string> vHist;
 	vHist.push_back("Ana/PhoJets/Hist/1Jet/Event/NJet15");
 	vHist.push_back("Ana/PhoJets/Hist/1Jet/Event/Ht");
 	vHist.push_back("Ana/PhoJets/Hist/1Jet/Event/Met");
 	vHist.push_back("Ana/PhoJets/Hist/1Jet/Photon/EtCorr");
 	vHist.push_back("Ana/PhoJets/Hist/1Jet/LeadJet/EtCorr");
 	vHist.push_back("Ana/PhoJets/Hist/1Jet/PhotonLeadJet/InvMass");
 	vHist.push_back("Ana/PhoJets/Hist/2Jet/Event/NJet15");
 	vHist.push_back("Ana/PhoJets/Hist/2Jet/Event/Ht");
 	vHist.push_back("Ana/PhoJets/Hist/2Jet/Event/Met");
 	vHist.push_back("Ana/PhoJets/Hist/2Jet/Photon/EtCorr");
 	vHist.push_back("Ana/PhoJets/Hist/2Jet/LeadJet/EtCorr");
 	vHist.push_back("Ana/PhoJets/Hist/2Jet/PhotonLeadJet/InvMass");
 	vHist.push_back("Ana/PhoJets/Hist/2Jet/Photon2Jets/InvMass");
 	vHist.push_back("Ana/PhoJets/Hist/2Jet/2Jets/InvMass");

/*	
 vHist.push_back("Hist/CENTRAL/1Jet/Event/NJet15");
   vHist.push_back("Hist/CENTRAL/1Jet/Event/Ht");
	vHist.push_back("Hist/CENTRAL/1Jet/Event/Met");
	vHist.push_back("Hist/CENTRAL/1Jet/Photon/EtCorr");
	vHist.push_back("Hist/CENTRAL/1Jet/LeadJet/EtCorr");
	vHist.push_back("Hist/CENTRAL/1Jet/PhotonLeadJet/InvMass");
	vHist.push_back("Hist/CENTRAL/2Jet/Event/NJet15");
	vHist.push_back("Hist/CENTRAL/2Jet/Event/Ht");
	vHist.push_back("Hist/CENTRAL/2Jet/Event/Met");
	vHist.push_back("Hist/CENTRAL/2Jet/Photon/EtCorr");
	vHist.push_back("Hist/CENTRAL/2Jet/LeadJet/EtCorr");
	vHist.push_back("Hist/CENTRAL/2Jet/PhotonLeadJet/InvMass");
	vHist.push_back("Hist/CENTRAL/2Jet/Photon2Jets/InvMass");
	vHist.push_back("Hist/CENTRAL/2Jet/2Jets/InvMass");
*/

	std::vector<std::string> vFileName;
	std::vector<TFile*> vFile;

	const int Nfiles = 200;
	for (int i=0; i<Nfiles ; ++i)
	{
		std::stringstream fname;
		fname << "PhoJets_MC_WW_" << i << "_" << Nfiles <<".root";
		vFileName.push_back(fname.str());
		TFile *f = new TFile(fname.str().c_str());
		if (f->IsZombie()) {std::cout << "file " << fname.str() << " is not found!" << std::endl; assert(false); }
		vFile.push_back(f);
	}

	
	for (int iHist=0; iHist < vHist.size(); ++iHist)
	{
		TH1* hist = 0;
		for (int iFile=0; iFile < vFile.size(); ++iFile)
		{
			vFile.at(iFile)->cd();
			//vFile.at(iFile)->cd("Hist/CENTRAL/1Jet/Event/");
			//gDirectory->ls();
			TH1* htemp = dynamic_cast<TH1*> (vFile.at(iFile)->Get(vHist.at(iHist).c_str()));
			if (htemp == NULL) {std::cout << "Hist " << vHist.at(iHist) << " not found in file " << vFile.at(iFile)->GetName() << std::endl; assert(false);}
			htemp->Print();
			if (hist == 0) { hist = htemp; hist->Sumw2();}
			else hist->Add(htemp);
		}

		
		if (vHist.at(iHist).find("1Jet") != std::string::npos 
			&& vHist.at(iHist).find("Event") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/1Jet/Event");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/1Jet/Event");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/1Jet/Event");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/1Jet/Event");
			hist->Write();
		}
		if (vHist.at(iHist).find("1Jet") != std::string::npos 
			&& vHist.at(iHist).find("Photon") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/1Jet/Photon");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/1Jet/Photon");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/1Jet/Photon");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/1Jet/Photon");
			hist->Write();
		}
		if (vHist.at(iHist).find("1Jet") != std::string::npos 
			&& vHist.at(iHist).find("LeadJet") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/1Jet/LeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/1Jet/LeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/1Jet/LeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/1Jet/LeadJet");
			hist->Write();
		}
		if (vHist.at(iHist).find("1Jet") != std::string::npos 
			&& vHist.at(iHist).find("PhotonLeadJet") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/1Jet/PhotonLeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/1Jet/PhotonLeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/1Jet/PhotonLeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/1Jet/PhotonLeadJet");
			hist->Write();
		}



		if (vHist.at(iHist).find("2Jet") != std::string::npos 
			&& vHist.at(iHist).find("Event") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/2Jet/Event");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/2Jet/Event");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/2Jet/Event");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/2Jet/Event");
			hist->Write();
		}
		if (vHist.at(iHist).find("2Jet") != std::string::npos 
			&& vHist.at(iHist).find("Photon") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/2Jet/Photon");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/2Jet/Photon");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/2Jet/Photon");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/2Jet/Photon");
			hist->Write();
		}
		if (vHist.at(iHist).find("2Jet") != std::string::npos 
			&& vHist.at(iHist).find("LeadJet") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/2Jet/LeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/2Jet/LeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/2Jet/LeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/2Jet/LeadJet");
			hist->Write();
		}
		if (vHist.at(iHist).find("2Jet") != std::string::npos 
			&& vHist.at(iHist).find("PhotonLeadJet") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/2Jet/PhotonLeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/2Jet/PhotonLeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/2Jet/PhotonLeadJet");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/2Jet/PhotonLeadJet");
			hist->Write();
		}


		if (vHist.at(iHist).find("2Jet") != std::string::npos 
			&& vHist.at(iHist).find("Photon2Jets") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/2Jet/Photon2Jets");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/2Jet/Photon2Jets");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/2Jet/Photon2Jets");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/2Jet/Photon2Jets");
			hist->Write();
		}
		if (vHist.at(iHist).find("2Jet") != std::string::npos 
			&& vHist.at(iHist).find("2Jets") != std::string::npos) 
		{
			outfile->cd();
			gDirectory->cd("Hist/CENTRAL/2Jet/2Jets");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESUP/2Jet/2Jets");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/EMJESDOWN/2Jet/2Jets");
			hist->Write();
			outfile->cd();
			gDirectory->cd("Hist/SIDEBAND/2Jet/2Jets");
			hist->Write();
		}
		new TCanvas();
		hist->DrawCopy();


	}


	outfile->ls();
	//outfile->Close();
}
