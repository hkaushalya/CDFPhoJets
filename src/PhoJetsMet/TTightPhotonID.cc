#include <iostream>
#include "samantha/PhoJetsMet/TTightPhotonID.hh"


ClassImp(TTightPhotonID)
//==========================================================
TTightPhotonID::TTightPhotonID(const char* name, const char* title):
	TStnModule(name,title),
	iIDWORD(-1)
{
}

//==========================================================
TTightPhotonID::TTightPhotonID(TStnPhoton* _pho)
{
	if ( _pho == NULL) std::cout << "TTightPhotonID:: WARNING! TStnPhoton pointer is NULL." << std::endl;
	else Init(_pho);
}
//==========================================================
void TTightPhotonID::Init(TStnPhoton* _pho)
{
		iDetector 	= _pho->Detector();
		fECorr		= _pho->ECorr();
		fEtCorr 		= _pho->ECorr() * _pho->SinTheta();
		fXCes			= _pho->XCes();
		fZCes 		= _pho->ZCes();
		fHadEm 		= _pho->HadEm();
		fChi2Mean 	= _pho->Chi2Mean();  			// (strip+wire)/2
		iN3d			= _pho->N3d();
		fCalIso		= _pho->EIso4(2);				//Corrected for leakage and MI
		fTrkPt  		= _pho->Pt();     				//using max pt track in the cluster
		fTrkIso 		= _pho->SumPt4();
		fCesWireEt2	= _pho->CesWireE2() * _pho->SinTheta();
		fCesStripEt2= _pho->CesStripE2() * _pho->SinTheta();
		iIDWORD = -1;
	
}	

//==========================================================
TTightPhotonID::TTightPhotonID(TSuperPhoton& _superpho)
{
	// need a assignment operator for super photon class!
	//if ( _superpho == NULL) std::cout << "TTightPhotonID:: ERROR! You passed a NULL reference!" << std::endl;
	//else {
		Init(_superpho);
}

//==========================================================
void TTightPhotonID::Init(TSuperPhoton& _superpho)
{
		std::cout << "CREATED TTIGHT with superpho" <<std::endl;
		iDetector 	= _superpho.Detector;
		fECorr 		= _superpho.ECorr;
		fEtCorr 		= _superpho.EtCorr;
		fXCes			= _superpho.XCes;
		fZCes			= _superpho.ZCes;
		fHadEm		= _superpho.HadEm;
		fChi2Mean	= _superpho.Chi2Mean;
		iN3d			= _superpho.N3d;
		fCalIso		= _superpho.IsoEtCorr;
		fTrkPt		= _superpho.TrkPt;
		fTrkIso		= _superpho.TrkIso;
		fCesWireEt2	= _superpho.CesWireE2 * _superpho.SinTheta;
		fCesStripEt2= _superpho.CesStripE2 * _superpho.SinTheta;
		iIDWORD = -1;
	
}	

//==========================================================
TTightPhotonID::~TTightPhotonID()
{
}

/*-------------------------------------------------------------------*/
void TTightPhotonID::SetCutValues()
{
	//define tight photon cuts here
	kCENTRAL			= 0;
	kETCORR			= 30.0;
	kXCES_MAX		= 21;
	kZCES_MIN   	= 9.0;
	kZCES_MAX		= 230.0;
	kHADEM_1			= 0.125;
	kHADEM_2			= 0.055 + 0.00045 * fECorr;
	kCHI2MEAN_MAX 	= 20.0;
	kN3DTRACKS_MAX = 1;
	kTRKPT_MAX   	= 1.0 + 0.005 * fEtCorr;
	kTRKISO_MAX  	= 2.0 + 0.005 * fEtCorr;

	if (fEtCorr < 20) kCALISO_MAX = 0.1 *fEtCorr;
	else kCALISO_MAX =  2.0 + 0.02 * (fEtCorr * 20.0);
	
	if (fEtCorr < 18) kCES2ET_MAX = 0.14 * fEtCorr;
	else kCES2ET_MAX = 2.4 + 0.01 * fEtCorr;
}

/*-------------------------------------------------------------------*/
int TTightPhotonID::GetIdWord(TSuperPhoton& _sp)
{
		Init(_sp);
		return GetIdWord();
}
/*-------------------------------------------------------------------*/
int TTightPhotonID::GetIdWord(TStnPhoton* _pho)
{
		Init(_pho);
		return GetIdWord();
}
/*-------------------------------------------------------------------*/
int TTightPhotonID::GetIdWord()
{
	SetCutValues();
	iIDWORD = 0x0;


	if (iDetector != kCENTRAL			) iIDWORD |= kCentralBit;
	if (fEtCorr < kETCORR 				) iIDWORD |= kEtCorrBit;
	if (fabs(fXCes) > kXCES_MAX		) iIDWORD |= kXCesBit;

	if ( ( fabs(fZCes) < kZCES_MIN ) ||
		  ( fabs(fZCes) > kZCES_MAX ) ) iIDWORD |= kZCesBit;

	if ( ! ( (fHadEm < kHADEM_1) ||
		      (fHadEm < kHADEM_2) )   ) iIDWORD |= kHadEmBit;

	if (fChi2Mean > kCHI2MEAN_MAX		) iIDWORD |= kChi2MeanBit;
	if (iN3d > kN3DTRACKS_MAX			) iIDWORD |= kN3dBit;
	if (fTrkPt > kTRKPT_MAX				) iIDWORD |= kTrkPtBit;
	if (fTrkIso > kTRKISO_MAX			) iIDWORD |= kTrkIsoBit;
	if (fCalIso > kCALISO_MAX			) iIDWORD |= kCalIsoBit;
	if (fCesWireEt2 > kCES2ET_MAX		) iIDWORD |= kCesWireE2Bit;
	if (fCesStripEt2 > kCES2ET_MAX	) iIDWORD |= kCesStripE2Bit;


	return iIDWORD;

} // GetIDword

/*-------------------------------------------------------------------*/
void TTightPhotonID::PrintIDWord()
{
	printf("ID WORD=%8i\n", iIDWORD);
}

/*-------------------------------------------------------------------*/
void TTightPhotonID::PrintIDBits()
{
	printf("kCentralBit 		= 0x1 << 0\n");
	printf("kEtCorrBit 		= 0x1 << 1\n");
	printf("kXCesBit 			= 0x1 << 2\n");
	printf("kZCesBit 			= 0x1 << 3\n");
	printf("kHadEmBit 		= 0x1 << 4\n");
	printf("kChi2MeanBit 	= 0x1 << 5\n");
	printf("kN3dBit 			= 0x1 << 6\n");
	printf("kTrkPtBit 		= 0x1 << 7\n");
	printf("kTrkIsoBit 		= 0x1 << 8\n");
	printf("kCalIsoBit 		= 0x1 << 9\n");
	printf("kCesWireE2Bit 	= 0x1 << 10\n");
	printf("kCesStripE2Bit 	= 0x1 << 11\n");
	
}

/*-------------------------------------------------------------------*/
void TTightPhotonID::PrintCutValues()
{
	printf("kCENRAL 	= %i\n", kCENTRAL);
	printf("kETCORR	= %6.2f\n",kETCORR);
	printf("kXCES_MAX	= %5.2f\n",kXCES_MAX);
	printf("kZCES_MAX	= %5.2f\n",kZCES_MAX);
	printf("kZCES_MIN	= %5.2f\n",kZCES_MIN);
	printf("kHADEM_1, kHADEM_2 = %6.4f, %6.4f\n",kHADEM_1, kHADEM_2);
	printf("kCHI2MEAN_MAX = %5.2f\n",kCHI2MEAN_MAX);
	printf("kN3DTRACKS	 = %2i\n",kN3DTRACKS_MAX);
	printf("kTRKPT_MAX	 = %6.2f\n",kTRKPT_MAX);
	printf("kTRKISO_MAX	 = %5.3f\n",kTRKISO_MAX);
	printf("kCALISO_MAX	 = %5.3f\n", kCALISO_MAX);
	printf("kCES2ET_MAX   = %5.3f\n",kCES2ET_MAX);
}
