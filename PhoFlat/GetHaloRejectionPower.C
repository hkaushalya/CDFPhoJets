#include<iostream>
#include<iomanip>
void GetHaloRejectionPower()
{

	//derive beam halo photon rejection power
	//halos tend to occupy phi wedges 0-23.
	//also has the flat cosmic background.
	//so subtract off the cosmic background from
	//wedges 0-23 by averaging the contents
	//in wedges 1-22.
	//repeat this after applying halo id cuts
	//and take the ratio to estimate the
	//halo rejection power.
	//
	//to do the systematics, assume one of the wedges
	//has no cosmics and repeat the procedure.


	
	//TFile *f = new TFile("HaloIdEff_g30.root");
	TFile *f = new TFile("HaloIdEff_g30_1njet15.root");
	assert (! f->IsZombie() && "root file not found!");

	TH1F *phiw_b4cuts = dynamic_cast<TH1F*> (f->Get("/Hist/HaloAndCosmic_b4/Photon/PhiWedge"));
	assert (phiw_b4cuts != NULL && "phi wedge before cut hist not found!");
	TH1F *phiw_a4cuts = dynamic_cast<TH1F*> (f->Get("/Hist/HaloAndCosmic_a4/Photon/PhiWedge"));
	assert (phiw_a4cuts != NULL && "phi wedge after cut hist not found!");

	new TCanvas();
	phiw_b4cuts->Draw();
	new TCanvas();
	phiw_a4cuts->Draw();

	double w0_total_b4 = phiw_b4cuts->GetBinContent(1);
	double w23_total_b4 = phiw_b4cuts->GetBinContent(24);

	double w0_total_a4 = phiw_a4cuts->GetBinContent(1);
	double w23_total_a4 = phiw_a4cuts->GetBinContent(24);
	

	double w1_22_cosmic_avg_b4=0;
	double w1_22_cosmic_avg_a4=0;

	const int iNPhiWedges = 24;
	assert (phiw_b4cuts->GetNbinsX() == iNPhiWedges && "number of bins is not 24!");
	assert (phiw_a4cuts->GetNbinsX() == iNPhiWedges && "number of bins is not 24!");
	
	for (int bin=2; bin < phiw_b4cuts->GetNbinsX(); ++bin)
	{
			w1_22_cosmic_avg_b4 += phiw_b4cuts->GetBinContent(bin);
			w1_22_cosmic_avg_a4 += phiw_a4cuts->GetBinContent(bin);
			std::cout << "bin / val = " << bin << ", " << phiw_b4cuts->GetBinContent(bin) << std::endl;
	}
	
	w1_22_cosmic_avg_b4 = w1_22_cosmic_avg_b4/ (double) (iNPhiWedges-2);
	w1_22_cosmic_avg_a4 = w1_22_cosmic_avg_a4/ (double) (iNPhiWedges-2);

	std::cout << std::setw(15) << "Before cuts" << std::setw(10) << "after cuts" << std::endl;
	std::cout << std::setw(10) << "W0 total  " << std::setw(5) << w0_total_b4 << std::setw(10) << w0_total_a4 << std::endl;
	std::cout << std::setw(10) << "W23 total " << std::setw(5) << w23_total_b4 << std::setw(10) << w23_total_a4 << std::endl;

	std::cout << std::setw(10) << "W1-22 Avg " << std::setw(5) << w1_22_cosmic_avg_b4 << std::setw(10) << w1_22_cosmic_avg_a4 << std::endl;

	const double halo_avg_b4 = fabs(w0_total_b4 + w23_total_b4 - 2 * w1_22_cosmic_avg_b4); 
	const double halo_avg_a4 = fabs(w0_total_a4 + w23_total_a4 - 2 * w1_22_cosmic_avg_a4);

	const double rej_pow = (halo_avg_a4/halo_avg_b4) * 100;

	std::cout << "halo rejection power (%) = " << rej_pow << std::endl;


	//systematics
	//assume one of the w0 or w23 has no cosmics but only halos.

	const double halo_avg_b4_syst = (w0_total_b4 + w23_total_b4 -  w1_22_cosmic_avg_b4); 
	const double halo_avg_a4_syst = (w0_total_a4 + w23_total_a4 -  w1_22_cosmic_avg_a4);

	const double rej_pow_syst = (halo_avg_a4_syst/halo_avg_b4_syst) * 100;

	std::cout << "halo rejection power systematic (%) = " << rej_pow_syst << std::endl;
	
}
