#include<iostream>
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"

// to make a nice gaus hist for my talk
void MakeGausHist()
{

	new TCanvas();
	gPad->SetGridx();
	gPad->SetGridy();
	TF1 *gaus = new TF1("F1","gaus",-2,2);
	gaus->SetTitle("Jet Energy Resolution: #frac{Measured}{Generated} - 1");
	gaus->SetParameters(2.7,0,0.3);
	gaus->Draw();


}
