#ifndef METOBJECT_HH
#define METOBJECT_HH

#include"TLorentzVector.h"
// this is the abstract base object for MetModel
//  holds common info for both jets and em objects

class MEtObject
{
	public:
		MEtObject();

		virtual void SetIndex(const unsigned int ind) { iBlockInd = ind; }
		virtual void SetRawMomentum(const TLorentzVector mom) { tlRawMomentum = mom; }
		virtual void SetCorMomentum(const TLorentzVector mom) { tlCorMomentum = mom; }
		virtual void SetNtraks(const unsigned int trks) { iNtrks = trks; }
		virtual void SetNtowers(const unsigned int twrs) { iNtwrs = twrs; }
		virtual void SetDetEta(const float eta) { fEtaDet = eta; }
		virtual void SetEmFration(const float emf) { fEmFr = emf; }


		virtual unsigned int Index() const { return iBlockInd; }
		virtual TLorentzVector RawMomentum() const { return tlRawMomentum; }
		virtual TLorentzVector CorMomentum() const { return tlCorMomentum; }
		virtual unsigned int Ntracks() const { return iNtrks; }
		virtual unsigned int Ntowers() const { return iNtwrs; }
		virtual float DetEta() const { return fEtaDet; }
		virtual float EmFraction() const { return fEmFr; }

	private:
		int iBlockInd; 				// original index in the data Block (this is assumed unique) 
		TLorentzVector tlRawMomentum;	// 4vec without corrections
		TLorentzVector tlCorMomentum;	// 4-vector after corrections
		int   iNtrks;  				// number of tracks in cluster (Ntrk as defined in Stntuple)
		int   iNtwrs;	 				// number of towers in cluster (Ntwr as defined in Stntuple)
		float fEtaDet;   				// detector eta 
		float fEmFr;  					// EM fraction 

};

#endif
