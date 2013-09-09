#ifndef LoosePhoCuts_HH
#define LoosePhoCuts_HH

#include "Stuple.hh"
#include <iostream>

class LoosePhoCuts
{
	public:
		LoosePhoCuts(const Stuple&);

		enum LooseIDBits_t {  // bit i will be set when the cut is failed
			kCentralBit 	= 0x1 << 0,
			kEtCorrBit		= 0x1 << 1,
			kXCesBit			= 0x1 << 2,
			kZCesBit			= 0x1 << 3,
			kHadEmBit 		= 0x1 << 4,
			kIso4Bit			= 0x1 << 5,
			kN3dBit			= 0x1 << 6,
			kTrkPtBit		= 0x1 << 7,
			kTrkIsoBit		= 0x1 << 8
		};

		unsigned int GetIdWord() const { return IDWord; }
		int IDword(const int iPhoInd);
		void SetTightenUpCut(const int tt) {
			TightenUpCut  = tt;	//1=hadem,2=iso,3=trkpt,4=trkiso
		};

	private:
		Stuple myStuple;
		unsigned int IDWord;
		int TightenUpCut;  // which cut should match to tight pho cuts

};
#endif
