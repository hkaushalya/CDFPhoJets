#include "TPad.h"
#include <iostream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TFrame.h"
//study custom pad creation in a canvas
void Pad()
{

	TCanvas *c1 = new TCanvas("c1","TCANVAS",600,800);
	c1->SetLeftMargin(0);
	c1->SetRightMargin(0);
	c1->SetTopMargin(0);
	c1->SetBottomMargin(0);


	TH1F *rat1 = new TH1F("rat1","rat1",80,-4,4);
	rat1->FillRandom("gaus",100);

	//TFrame *fr = new TFrame (0.01,0.4,0.99,0.99);
	//fr->Draw();
	

	TPad *p1 = new TPad("log","Log Plot",0.01,0.4,0.99,0.99); 	//in NDC cdts
	//p1->SetFillColor(kYellow);
	p1->SetBorderMode(0);
	p1->SetLineColor(kBlack);
	p1->SetBottomMargin(0);
	p1->Draw();
	p1->cd();
	
	rat1->Draw();
	TFrame *fr = c1->GetFrame();
	fr->Print();
	
	c1->cd();

	TPad *p2 = new TPad("ratio","Ratio Plot",0.01,0.01,0.99,0.4); 	//in NDC cdts
	//p2->SetFillColor(kBlue);
	p2->SetBorderMode(0);
	//p2->SetLineColor(kBlue);
	//p2->SetFrameLineColor(kBlack);
	p2->SetTopMargin(0);
	p2->Draw();
	p2->cd();

	TH1F *rat2 = new TH1F("rat2","rat2",80,-4,4);
	rat2->FillRandom("gaus",100);
	rat2->Draw();

	c1->cd();
};
