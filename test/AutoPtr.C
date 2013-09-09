#include<iostream>
#include"TH1F.h"
#include<memory>


TH1F* FillSome(std::auto_ptr<TH1F> hist, const int rebin)
{

	hist->FillRandom("pol3",100);
	hist->Print();
	return hist.release();
}

void AutoPtr()
{
	//study the auto_ptr 
	std::auto_ptr<TH1F> h1 (new TH1F("H1","TEST HIST",100,0,100));

	TH1*  h2 = FillSome(h1,10);
	h2->Draw();
}
