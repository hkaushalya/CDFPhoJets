#ifndef HISTMANAGER_HH
#define HISTMANAGER_HH

#include "Stuple.hh"
#include "Histograms.hh"
#include "FreeFunctions.hh"
#include "PhotonList.hh"
#include "ElectronList.hh"
#include "JetList.hh"
#include "CommonVars.hh"

class HistManager
{
	public:
		HistManager();
		~HistManager();
		void FillEventHists(const Stuple&, Histograms::EventHists_t& , 
								const float fWeight=1);
		//here the Photon and  the Jet Index is the order that they are in the branch!
		//collection=0(central),1(em/jes up),2(em/jes down)
		void FillPhotonHists(const Stuple&, const int collection, Histograms::PhotonHists_t&,
								const int iPhoIndex, const float fWeight=1);									
		void FillElectronHists(const Stuple&, const int collection, Histograms::PhotonHists_t& ,
								const int iEleIndex, const float fWeight=1);
		void FillJetHists(const Stuple&, const int collections, Histograms::JetHists_t& , 
								const int iJetIndex,	const float fWeight=1);
		void FillPhoton1JetHists(const Stuple&, const int collection, Histograms::Photon1JetHists_t& ,
								const int iPhoIndex, const int iJetIndex,
								const float fWeight=1);
		void FillElectron1JetHists(const Stuple&, const int collection, Histograms::Photon1JetHists_t& ,
								const int iEleIndex, const int iJetIndex,
								const float fWeight=1);
		void FillPhoton2JetsHists(const Stuple&, const int collection, Histograms::Photon2JetsHists_t& hist,
								const int iPhoIndex, const int iJet1Index,
								const int iJet2Index, const float fWeight=1);
		void FillElectron2JetsHists(const Stuple&, const int collection, Histograms::Photon2JetsHists_t& hist,
								const int iEleIndex, const int iJet1Index,
								const int iJet2Index, const float fWeight=1);
		void FillTwoJetsHists(const Stuple&, const int collection, Histograms::TwoJetsHists_t& hist,
								const int iJet1Index, const int iJet2Index,
								const float fWeight=1);


	
};
#endif
