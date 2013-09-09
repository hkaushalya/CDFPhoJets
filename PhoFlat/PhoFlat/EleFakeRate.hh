#ifndef ELEFAKERATE_HH
#define ELEFAKERATE_HH

#include "Stuple.hh"

class EleFakeRate 
{
		public:
			EleFakeRate();

			enum { FkRtParArray = 4 };
			enum { PhnxParArray = 3 };

			double fakerate_par[FkRtParArray]; // fake rate parametrization
			double fakerate_parerr[FkRtParArray]; // fake rate parametrization uncertainty
			double fakerate_corr_coeff[FkRtParArray][FkRtParArray]; // covariance matrix coefficients
			double phnx_par[PhnxParArray]; // phnx eff. parametrization
			double phnx_parerr[PhnxParArray]; // phnx eff. parametrization uncertainty
			double phnx_corr_coeff[PhnxParArray][PhnxParArray]; // covariance matrix coefficients

			void MyFakeRate(double eleEt, int goodrun_bit, int phnx_bit, 
							double &wght_d, double &wght_m, double &wght_p); // returns fake rate & syst
			double MyShiftCorrFunc(double eleEt); // correction for Et(det)/Et(genp) difference
			double MyPhnxEffFunc(double eleEt);  // phoenix efficiency function
			double MyPhnxEffFuncErr(double eleEt);  // uncertainty on phoenix efficiency 
			double MyFakeRateFunc(double eleEt); // returns fake rate
			double MyFakeRateFuncErr(double eleEt); // returns uncertainty on fake rate

};
#endif
