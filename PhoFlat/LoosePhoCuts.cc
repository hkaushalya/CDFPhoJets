////////////////////////////////////////////////////////////////
// can check for a photon passing Loose photon ID cuts.       //
////////////////////////////////////////////////////////////////


#include "PhoFlat/LoosePhoCuts.hh"
#include <iostream>



//-------------------------------------------------------------------
LoosePhoCuts::LoosePhoCuts(const Stuple& st) 
{
	myStuple = st;
}

//-------------------------------------------------------------------
int LoosePhoCuts::IDword(const int iPhoIndex)
{
	IDWord = 0x0;			//the id word will hold the last checked obj info.
	
	if (myStuple.pho_Detector[iPhoIndex] != 0						) {
		IDWord |= kCentralBit;
	}

	if (myStuple.pho_Etc[iPhoIndex] < 7. 						) {
		IDWord |= kEtCorrBit;
	}
	
	if (fabs(myStuple.pho_XCes[iPhoIndex]) > 21					) {
		IDWord |= kXCesBit;
	}
	
	if (fabs(myStuple.pho_ZCes[iPhoIndex]) < 9 ||  fabs(myStuple.pho_ZCes[iPhoIndex]) > 230) {
		IDWord |= kZCesBit;
	}
	
	if (myStuple.pho_N3d[iPhoIndex] > 1								) {
		IDWord |= kN3dBit;
	}


	
	if (TightenUpCut == 1) {  // hadem
		//std::cout << "\t st E, Et, HadEm=" << myStuple.pho_E[iPhoIndex] << "," << myStuple.pho_Etc[iPhoIndex] << "," << myStuple.pho_HadEm[iPhoIndex]<< std::endl;
		if ( !((myStuple.pho_HadEm[iPhoIndex] < 0.125) || (myStuple.pho_HadEm[iPhoIndex] < (0.055 + 0.00045 * myStuple.pho_E[iPhoIndex])))) {
			IDWord |= kHadEmBit;
		} //else std::cout << "I FAILED" <<std::endl;
		
	} else {
		if (myStuple.pho_HadEm[iPhoIndex] > 0.125) {
			IDWord |= kHadEmBit;
		}
	}

	if (TightenUpCut == 2) { // iso
		if (myStuple.pho_Etc[iPhoIndex] < 20) {
			if (myStuple.pho_Iso4[iPhoIndex] > 0.1*myStuple.pho_Etc[iPhoIndex]) {
				IDWord |= kIso4Bit;
			}
		} else {
			if (myStuple.pho_Iso4[iPhoIndex] > (2.0+0.02*(myStuple.pho_Etc[iPhoIndex]-20.0)) ) {
				IDWord |= kIso4Bit;
			}
		}

	} else {
		if (myStuple.pho_Etc[iPhoIndex] < 20) {
			if (myStuple.pho_Iso4[iPhoIndex] > 0.15*myStuple.pho_Etc[iPhoIndex]) {
				IDWord |= kIso4Bit;
			}
		} else {
			if (myStuple.pho_Iso4[iPhoIndex] > 3.0 + 0.02 * (myStuple.pho_Etc[iPhoIndex] - 20.0)) {
				IDWord |= kIso4Bit;
			}
		}
	}
	
	if (TightenUpCut == 3) {	//trakpt
		if (myStuple.pho_TrkPt[iPhoIndex] > (1+0.005*myStuple.pho_Etc[iPhoIndex])		) {
			IDWord |= kTrkPtBit;
		}
	} else {
		if (myStuple.pho_TrkPt[iPhoIndex] > (0.25*myStuple.pho_Etc[iPhoIndex])	) {
			IDWord |= kTrkPtBit;
		}
	}

	if (TightenUpCut == 4) {	//trackiso
		if (myStuple.pho_TrkIso[iPhoIndex] > (2.0+0.005*myStuple.pho_Etc[iPhoIndex])	) {
			IDWord |= kTrkIsoBit;
		}
		
	} else {
		if (myStuple.pho_TrkIso[iPhoIndex] >5.0 ) {
			IDWord |= kTrkIsoBit;
		}
	}
	
	return IDWord;

}
