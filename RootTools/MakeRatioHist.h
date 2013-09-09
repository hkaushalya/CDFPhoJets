#ifndef MAKERATIOHIST_HH
#define MAKERATIOHIST_HH
#include<TH1.h>


class MakeRatioHist
{
	private:

	public:
		MakeRatioHist();
		void RatioHist(const TH1 *data, const TH1 *bckg, const std::string sTitle="");
		void LogHist(const TH1 *data, const TH1 *bckg);
		void TestMethod(const int njet);
		void TestMethod2(const std::string,const std::string);
		float CorrelatedErrors(const float val1, const float err1,
										const float val2, const float err2);

};
#endif
