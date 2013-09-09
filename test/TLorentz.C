#include <iostream>
#include "TLorentzVector.h"

void Print(const TLorentzVector& tlV)
{
	std::cout << "Pt= "<< tlV.Pt() << std::endl;
	std::cout << "V.V= "<< tlV.Dot(tlV) << std::endl;
	std::cout << "V.M()= "<< tlV.M() << std::endl;
	std::cout << "V.E()= "<< tlV.E() << "\n" << std::endl;

}
void TLorentz()
{

	//TLorentzVector tlV(-0.81,-22.42,-143.45,145.37);
	TLorentzVector tlV(1.44,35.38,-26.30,44.11);
	Print(tlV);
	tlV.SetE(47);
	Print(tlV);
	tlV.Print();

}
