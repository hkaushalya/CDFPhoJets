#include <iostream>
#include "TH1.h"
#include "TH1F.h"
#include "../RootTools/IOColors.hh"
#include <cmath>

void Sumw2()
{

	TH1F *den = new TH1F("den","denom hist", 20,0,20);
	den->Sumw2();
	TH1F *num = new TH1F("num","numer hist", 20,0,20);
	num->Sumw2();

	den->Fill(1);
	den->Fill(1);
	den->Fill(10);
	den->Fill(10);
	den->Fill(10);
	den->Fill(10);
	den->Fill(10);
	
	num->Fill(1);
	num->Fill(10);
	num->Fill(10);
	num->Fill(10);
	num->Fill(10);

	num->Print("all");
	//return;
	

	float x,y,dx, dy, dp;
 	 for (int bin = 1; bin < num->GetNbinsX(); ++bin)
	 {
		if (num->GetBinContent(bin))
		{
			std::cout << green << "num bin=" << bin << "\t" << num->GetBinContent(bin) << "/" << num->GetBinError(bin) << clearatt << std::endl;
			x = num->GetBinContent(bin);
			dx = num->GetBinError(bin);
		}
	 }

 	 for (int bin = 1; bin < den->GetNbinsX(); ++bin)
	 {
		if (den->GetBinContent(bin))
		{
			std::cout << red << "den bin=" << bin << "\t" << den->GetBinContent(bin) << "/" << den->GetBinError(bin) << clearatt << std::endl;
			y = den->GetBinContent(bin);
			dy = den->GetBinError(bin);
		}
	 }

	 num->Divide(den);


 	 for (int bin = 1; bin < num->GetNbinsX(); ++bin)
	 {
		if (num->GetBinContent(bin))
			std::cout << blue << "num bin=" << bin << "\t" << num->GetBinContent(bin) << "/" << num->GetBinError(bin) << clearatt << std::endl;
	 }

	 dp = x/y * sqrt ( pow(dy/y,2) + pow(dx/x,2) );
	 std::cout << magenta << "my dp = " <<  dp << clearatt << std::endl;
}
