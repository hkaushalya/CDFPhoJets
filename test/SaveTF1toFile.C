#include<iostream>
#include "TFile.h"
#include "TF1.h"
void SaveTF1toFile()
{
	TFile *f = new TFile("SaveTF1.root","RECREATE");
	TF1 *tf = new TF1("f1","x+x^2");
	tf->Write();
	f->Close();

	std::cout << "reading " << std::endl;
	f = new TFile("SaveTF1.root","UPDATE");
	f->ls();
	TF1 *tf2 = (TF1*) f->Get("f1");
	assert (tf2 != NULL);
	tf2->Print();
	//f->Delete(tf2->GetName());  //here you need to specify the name cycle. else no use
	f->Delete("f1;*"); //this one works
	f->ls();
	f->Close();

	std::cout << "reading after deletion" << std::endl;
	f = new TFile("SaveTF1.root");
	f->ls();
	f->Close();


}
