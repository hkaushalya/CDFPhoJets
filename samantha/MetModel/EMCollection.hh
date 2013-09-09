#ifndef EMCOLLECTION_HH
#define EMCOLLECTION_HH

#include<vector>
#include "TLorentzVector.h"

// This is a common object to hold all em objects
// ele/pho/muon/tau etc
// create a collection for each different type of EM Objects
// you have for MetModel
// for eg: EMcollection tightElectrons; (or Electrons including both tight and loose)


class EMCollection
{
	public:
		EMCollection();
		EMCollection(const int iCollType);


		enum {
			iELECTRON = 0,
			iPHOTON   = 1,
			iMUON     = 2,
			iTAU      = 3,
		};

	private:
			int iType;													// 0=electron,1=photon,2=muon,3=tau
			std::vector<TLorentzVector> tlvRawMomentum;    	// raw momentum
			std::vector<TLorentzVector> tlvCorrMomentum;   	// corrected momentum (leakage, em scale,etc.)
			std::vector<int> ivBlockIndex; 						// original index in photon/electron .. block
			std::vector<float> fvEmFr;      				// EM energy fraction
			std::vector<float> fvEtaDet;    						// detector eta for photons
			std::vector<float> fvXces;    						// X_ces for photons
			std::vector<float> fvZces;    						// Z_ces for photons

};
#endif
