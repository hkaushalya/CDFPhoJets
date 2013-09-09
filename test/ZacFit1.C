#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"

void ZacFit1(const unsigned int pts =1000)
{
	TF1 *f1 = new TF1("f1","expo",-20,20);
	
	f1->SetParameter(0, 0.0001);
	f1->SetParameter(1, 0.02);
	f1->SetLineColor(5);
	f1->Draw();
}
