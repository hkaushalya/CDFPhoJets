#ifndef JERFUNCTION_HH
#define JERFUNCTION_HH

#include <vector>
#include <memory>
#include "TF1.h"

class JERFunction
{
	
	public:
		enum JERCLASS_t
		{
			DEFAULT= 1,
			SYST   = 2
		};
		
		JERFunction(JERCLASS_t, const double fJetE, const double fDetEta, const int syst=0);
		JERFunction();

		int EtaBin(const double fEta) const;				//returns eta bin
		double RandomPoint() const;
		double Integral(const double fXmin, const double fXmin) const; 
		double JERParameter(const int iEtaBin, const int iPar, const int iSystCode); //sets of parameters for different JER classes


	private:
		std::vector< std::vector<double> > vParam;			//etabin,parameter
		std::auto_ptr<TF1> FuncTf1;
		JERCLASS_t JERCLASS;										//which JER function selection is used


};
#endif
