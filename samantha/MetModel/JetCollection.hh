#ifndef JETCOLLECTION_HH
#define JETCOLLECTION_HH

#include <vector>
#include "TLorentzVector.h"
#include "Jet.hh"


//this is container for different jet collections
// for MetModel. jets with different correction levels.
// for eg: JetCollection JetsLev6NoEmObj; 


class JetCollection
{
	public:
		JetCollection();

		struct NJets_t {
			int iNjet; 		// number of clusters ("raw jets")
			int iNjet_5; 	// number of level-6 jets with Et>5, after matching with EM object
			int iNjet_10; 	// number of level-6 jets with Et>10, after matching with EM object
			int iNjet_15; 	// number of level-6 jets with Et>15, after matching with EM object
			int iNjet_20; 	// number of level-6 jets with Et>20, after matching with EM object
			int iNjet_25; 	// number of level-6 jets with Et>25, after matching with EM object
			int iNjet_30; 	// number of level-6 jets with Et>30, after matching with EM object
			int iNjet_35; 	// number of level-6 jets with Et>35, after matching with EM object
		};

		//------------------- parameters for Jet Correction initialization
		struct JetCorrectPara {
			int iLevel;					//jet correction level
			int iNvx;					//number of vertices
			int iJetCone;				//cone size of the jet collection //0.4 for now
			int iVersion;				//version of the corrections
			int iSys;					// systematic code
			int iRunNumber;			// Run number for run dependent corrections
			int iImode;					// this is similar to name used in JER intructions page for jet corrections
		};

		struct Sumet {
			float fRaw;					//SumEt from MEt block
			float fTh5;					//SumEt calcuated with Jets with Et>5GeV
			float fTh10;				//SumEt calcuated with Jets with Et>10GeV
			float fTh15;				//SumEt calcuated with Jets with Et>15GeV
			float fTh20;				//SumEt calcuated with Jets with Et>20GeV
			float fTh25;				//SumEt calcuated with Jets with Et>25GeV
			float fTh30;				//SumEt calcuated with Jets with Et>30GeV
			float fTh35;				//SumEt calcuated with Jets with Et>35GeV
		};

		struct Met {
			TVector2 fRaw;					//MEt from MEt block
			TVector2 fTh5;					//MEt calcuated with Jets with Et>5GeV
			TVector2 fTh10;				//MEt calcuated with Jets with Et>10GeV
			TVector2 fTh15;				//MEt calcuated with Jets with Et>15GeV
			TVector2 fTh20;				//MEt calcuated with Jets with Et>20GeV
			TVector2 fTh25;				//MEt calcuated with Jets with Et>25GeV
			TVector2 fTh30;				//MEt calcuated with Jets with Et>30GeV
			TVector2 fTh35;				//MEt calcuated with Jets with Et>35GeV
		};

		float SumEtJets(const float JetEtThreshold) const;		//return SumEt of Jets above the given threhold
		void ClearAll();								//clear all vectors and reset all variables to defaults

	private:
		NJets_t Njet;									//all njets at different thresholds
		JetCorrectPara JetCorrPara;  			// general parameters for this jet collection
		//std::vector<Jet> vJets;					// All the jets
		std::vector<float> vfSmearFactors;	// generated smear factors for each jet

};

#endif
