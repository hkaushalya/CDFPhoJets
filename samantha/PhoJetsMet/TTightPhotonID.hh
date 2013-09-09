/*
what this should do
have all the photons id cuts

return id-word for a given photon

*/


#ifndef TTIGHTPHOTONID_HH
#define TTIGHTPHOTONID_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <Stntuple/loop/TStnModule.hh>
#include <Stntuple/obj/TStnPhoton.hh>
#include "Stntuple/ana/TSuperPhoton.hh"
#endif


class TTightPhotonID: public TStnModule {

	protected:

		enum TightIDBits_t {  // bit i will be set when a cut is failed
			kCentralBit 	= 0x1 << 0,
			kEtCorrBit		= 0x1 << 1,
			kXCesBit			= 0x1 << 2,
			kZCesBit			= 0x1 << 3,
			kHadEmBit 		= 0x1 << 4,
			kChi2MeanBit	= 0x1 << 5,
			kN3dBit			= 0x1 << 6,
			kTrkPtBit		= 0x1 << 7,
			kTrkIsoBit		= 0x1 << 8,
			kCalIsoBit		= 0x1 << 9,
			kCesWireE2Bit 	= 0x1 <<10,
			kCesStripE2Bit = 0x1 <<11
		};
		
		
	public:
		TTightPhotonID(const char* name="TightPhotonID", const char* title="Tight Photon IDs");
		TTightPhotonID(TStnPhoton* pho);
		TTightPhotonID(TSuperPhoton& superpho);
		~TTightPhotonID();
	
		int GetIdWord();
		int GetIdWord(TStnPhoton* pho);
		int GetIdWord(TSuperPhoton& superpho);
		void PrintIDWord();
		void PrintCutValues();
		void PrintIDBits();


	private:
		int iIDWORD;
		//holds the current object's valuse
		int   iDetector;
		float fECorr;
		float fEtCorr;
		float fXCes;
		float fZCes;
		float fHadEm;
		float fChi2Mean;
		int   iN3d;
		float fCalIso;
		float fTrkPt;
		float fTrkIso;
		float fCesWireEt2;
		float fCesStripEt2;

		void Init(TStnPhoton*);
		void Init(TSuperPhoton&);
		void SetCutValues();
		//holds tight photon cuts values
		int 	kCENTRAL;
		float	kETCORR;
		float	kXCES_MAX;
		float	kZCES_MIN;
		float	kZCES_MAX;
		float	kHADEM_1	;
		float	kHADEM_2;
		float	kCHI2MEAN_MAX;
		int 	kN3DTRACKS_MAX;
		float	kTRKPT_MAX;
		float	kTRKISO_MAX;
		float kCALISO_MAX;
		float kCES2ET_MAX;
		
	ClassDef(TTightPhotonID,0)
};

#endif
