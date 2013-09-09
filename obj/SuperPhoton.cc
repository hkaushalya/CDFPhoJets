///////////////////////////////////////////////////////////
// This will hold all the information about the photon.  //
// Almost all the information in the TStnPhoton class    //
// and some more stuff like Tight Photon, Cosmic, Beam   //
// Halo, ... tags.                                       //
///////////////////////////////////////////////////////////
// Author: Samantha Hewamanage <samantha@fnal.gov>

#include "samantha/obj/SuperPhoton.hh"
#include <fstream>
#include <iostream>
#include <iomanip>

//______________________________________________________________
SuperPhoton::SuperPhoton() :
	iIndex(-1),
	aPhoton(0),
	tlCorVec(0,0,0,0),
	tlRawVec(0,0,0,0),
	iDetector(-1),
	fDetEta(-99999.99),
	fDetPhi(-99999.99),
	fEtRaw(-99999.99),
	fECorr(-99999.99),
	fEtCorr(-99999.99),
	fXCes(-99999.99),
	fZCes(-99999.99),
	fHadEm(-99999.99),
	fChi2Mean(-99999.99),
	iN3d(-99999),
	fIsoEtCorr(-99999.99),
	fTrkPt(-99999.99),
	fTrkPt2(-99999.99),
	fTrkIso(-99999.99),
	fCesWireE2(-99999.99),
	fCesStripE2(-99999.99),
	iLoosePhotonId(-1),
	iTightPhotonId(-1),
	iLooseElectronId(-1),
	iTightElectronId(-1),
	iIsBeamHalo(-1),
	iIsCosmic(-1),
	iIsConversion(-1),
	iIsPhoenix(-1),
	iIsPMTSpike(-1),
	fHadTdcTime(-99999.99),
	fSinTheta(-99999.99),
	fEmTime(-99999.99),
	iTimeWindow(-1),
	iStdLooseElectronId(-1),
	iStdTightElectronId(-1),
	iPrintLevel(10)
{
	// Constructor
}

//______________________________________________________________
SuperPhoton::~SuperPhoton()
{
	// Destructor
}

//______________________________________________________________
void SuperPhoton::Print(Option_t *opt, const int iPrtLevel)
{
	// Print all data memeber values
	printf("============= SuperPhoton:: Values of Current Data Members =============\n");
	printf("Index Pho Blk, Ele Blk = %3i , %3i\n"  , GetPhotonBlockIndex(), GetElectronBlockIndex());  
	printf("Detector ------------ = %9i\n"  , GetDetector());  
	std::cout << "Detector ------------ = " << std::setw(8) << GetDetector() << std::endl;  
//	printf("Corr E -------------- = %5.3f\n", GetECorr());  
	printf("Corr Et ------------- = %5.3f\n", GetEtCorr());  
	printf("Loose Photon IDword   = %9i\n"  , GetLoosePhotonId());  
	printf("Tight Photon IDword   = %9i\n"  , GetTightPhotonId());  
	printf("PhoLike Loose Ele ID  = %9i\n"  , GetLooseElectronId());  
	printf("PhoLike Tight Ele ID  = %9i\n"  , GetTightElectronId());  
	printf("Std. Loose Ele ID --- = %9i\n"  , GetLooseElectronId());  
	printf("Std. Tight Ele ID --- = %9i\n"  , GetTightElectronId());  
	printf("IsBeamHalo? --------- = %9i\n"  , GetBeamHaloId());  
	printf("IsCosmic? ----------- = %9i\n"  , GetCosmicId());  
	printf("IsConversion? ------- = %9i\n"  , GetConversionId());  
	printf("IsPhoenix? ---------- = %9i\n"  , GetPhoenixId());  
//	printf("IsPMTSpike? --------- = %9i\n"  , GetPMTSpikeId());  
//	printf("Raw Et -------------- = %5.3f\n", GetEtRaw());  
/*	printf("Det. Eta ------------ = %5.3f\n", GetDetEta());  
	printf("Det. Phi ------------ = %5.3f\n", GetDetPhi());  
	printf("XCes ---------------- = %5.3f\n", GetXCes());  
	printf("ZCes ---------------- = %7.1f\n", GetZCes());  
	printf("Had/Em -------------- = %5.3f\n", GetHadEm());  
	printf("Mean Chi^2 ---------- = %5.3f\n", GetChi2Mean());  
	printf("N3d tracks ---------- = %9i\n"  , GetN3d());  
	printf("Corr. Iso(0.4) ------ = %6.2f\n", GetIso4());  
	printf("Trk. Pt ------------- = %6.2f\n", GetTrkPt());  
	printf("Trk. Pt2 ------------ = %6.2f\n", GetTrkPt2());  
	printf("Trk. Iso ------------ = %6.2f\n", GetTrkIso()); 
	printf("2nd CES(Wire) E ----- = %6.2f\n", GetCesWireE2());  
	printf("2nd CES(Strip) E ---- = %6.2f\n", GetCesStripE2());  
	printf("SinTheta ------------ = %6.2f\n", GetSinTheta());  
	printf("EM time ------------- = %6.2f\n", GetEmTime());  
*/
	std::string sWin("N/A");
	/*
	if (IsEarly()) sWin = saTimeWindows[0];
	if (IsTimeBeamHalo()) sWin = saTimeWindows[1];
	if (IsInTime()) sWin = saTimeWindows[2];
	if (IsLate()) sWin = saTimeWindows[3];
	*/
	if (IsEarly()) sWin = "Early Window";
	if (IsTimeBeamHalo()) sWin = "BeamHalo Window";
	if (IsInTime()) sWin = "InTime Window";
	if (IsLate()) sWin = "Late Window";
/*	printf("Time Window --------- = %s\n", sWin.c_str());  
	printf("Had TDC time -------- = %6.2f\n", GetHadTdcTime());  
	printf("4-Vec Raw(E,Px,Py,Pz) =  %6.2f, %6.2f, %6.2f, %6.2f\n",GetRawVec().E(), GetRawVec().Px(), GetRawVec().Py(), GetRawVec().Pz() );  
	printf("4-Vec Cor(E,Px,Py,Pz) =  %6.2f, %6.2f, %6.2f, %6.2f\n\n",GetCorVec().E(), GetCorVec().Px(), GetCorVec().Py(), GetCorVec().Pz() );  
*/
}

