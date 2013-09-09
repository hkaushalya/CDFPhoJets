#ifndef OVERLAYHISTS_H
#define OVERLAYHISTS_H

#include <TProfile.h>
#include <string>
#include "Rtypes.h"
#include "TColor.h"

class OverlayHists
{

	public:
		OverlayHists();
		void TwoProfileHists(const TProfile* h1, const TProfile* h2,
							const std::vector<std::string> legend, const std::string finalHistTitle);
		void ThreeProfileHists(const TProfile* h1, const TProfile* h2, const TProfile* h3, 
							const std::vector<std::string> legend, const std::string finalHistTitle);

		void Two1DHists(const TH1* hist1, const TH1* hist2,
				const std::string legend1, const std::string legend2,
				const std::string finalTitle="",
				const std::string xtitle="", const std::string ytitle="",
				const bool norm=0, const bool individuals=1, const int rebin=1, 
				const double fLoX=-9999.0, const double fUpX=-9999.0);
		
		void Three1DHists(const TH1* hist1, const TH1* hist2, const TH1* hist3,
				const std::vector<std::string> legend, const std::string finalTitle, const int rebin=1);
		
		void SetMakeIndividual(const bool i) { bMakeIndividual = i; }
		bool MakeIndividual() const { return bMakeIndividual; }
		void SetNormBgToData(const bool norm) { bNormBgToData = norm; }
		bool NormBgToData() const { return bNormBgToData; }
		void SetAddBinWToYtitle(const bool i) { bBinWtoYtitle = i; }
		bool AddBinWToYtitle() const { return bBinWtoYtitle; }
		

	private:
		bool bMakeIndividual;		// make individual plots in addition to the overlayed
		bool bNormBgToData;			//normalize background to data 
		bool bBinWtoYtitle;			// adds the bin width to y title
		Color_t lineColor1;
		Color_t lineColor2;
		Width_t lineWidth;
		Style_t markerStyle1, markerStyle2;
		std::string drawOptions;
};

#endif
