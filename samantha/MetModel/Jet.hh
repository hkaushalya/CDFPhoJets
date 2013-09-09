#ifndef JET_HH
#define JET_HH

#include"TLorentzVector.h"
#include "Stntuple/obj/TStnJet.hh"
#include "samantha/MetModel/MEtObject.hh"
// this is the Jet object for MetModel
//  holds common info for a jet
// JetCollection-s can be created 
// with different jet correction levels
//  (rawJetColl, lev1JetColl etc)

class Jet : public MEtObject
{
	public:
		Jet();

	private:
		float fEtaDetCorr;   			// jet detector eta (after removing matched EM objects)
		float fEmFrCorr; 				// EM fraction (after removing matched EM objects)
		int   iNobjMatch; 				// total number of matched objects (pho,ele,mu,tau,btag) for this jet collection
		int   iNphoMatch; 				// number of matched photons for this jet collection
		int   iNeleMatch; 				// number of matched electrons for this jet collection

};

#endif
