#include <iostream>
#include "TH1F.h"
#include "TF1.h"

void LandauSigma()
{
	TF1* f1 = new TF1("lan","sqrt(25 + 100/x)",-100,100);	
	f1->SetLineColor(kRed);
	f1->SetLineWidth(2);
	TF1* f2 = new TF1("lan2","sqrt(25 - 100/x)",-100,100);	
	f2->SetLineColor(kBlue);
	f2->SetLineWidth(2);
	f1->Draw();
	f2->Draw("SAME");
}
