#ifndef METMODEL_HH
#define METMODEL_HH

#include "JetCollection.hh"
#include "EMCollection.hh"
#include "TLorentzVector.h"
#include <Stntuple/loop/TStnModule.hh>
//------declaration of blocks to be read
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include <Stntuple/obj/TStnJetBlock.hh>
#include <Stntuple/obj/TStnMetBlock.hh>
#include "JetUser/JetEnergyCorrections.hh"
//----- need this for jet-EM matching
#include <Stntuple/obj/TCalDataBlock.hh>
#include <CalorGeometry/CalConstants.hh>
#include <Stntuple/obj/TStnPhotonBlock.hh>
#include <Stntuple/obj/TStnElectronBlock.hh>



class MetModel: public TStnModule
{
	protected:
		TStnJetBlock*      fJetBlockClu04;
		TStnMetBlock*      fMetBlock;		// sam added
		TStnHeaderBlock*   fHeaderBlock;	// sam added

		//----------------------- need this for jet-EM matching
		TCalDataBlock* fCalData;
		TStnPhotonBlock* fPhotonBlock;
		TStnElectronBlock* fElectronBlock;
	//	typedef std::vector<TCalTower*> CalDataArray; 				// need for beam halo
	//	typedef std::vector<TCalTower*>::iterator CalDataArrayI; // need for beam halo

	
	
	public:
		MetModel(const char* name="MetModel", const char* title="MetModel");
		~MetModel();

		// ****** accessors
		TStnHeaderBlock*   GetHeaderBlock()  const { return fHeaderBlock;  }

		// ****** over ridden methods of TStnModule
		int BeginJob();
		int BeginRun();
		int Event   (int ientry);
		int EndJob  ();

		// ****** other methods
		void BookHistograms();
		void FillDataBlocks(int ientry);


		struct Jets_t {
			JetCollection rawJets;
			JetCollection lev1Jets;
			JetCollection lev4Jets;
			JetCollection lev5Jets;
			JetCollection lev6Jets;
			JetCollection lev7Jets;
			JetCollection metAddedJets;
			JetCollection smearedJets;
		};
		
		struct EMObjects_t {
			EMCollection Photons;
			EMCollection Electrons;
		};


		void ClearAll(); 		//clear all for next event
		void GetRawJets();
		void GetElectrons(EMObjects_t&);
		void GetEMobjects(EMObjects_t& emobjs);


	private:
		EMObjects_t EMobjects;
		Jets_t Jets;


	ClassDef(MetModel,1);
};
#endif
