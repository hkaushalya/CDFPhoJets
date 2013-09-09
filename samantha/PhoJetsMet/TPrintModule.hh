#ifndef TPRINTMODULE_HH
#define TPRINTMODULE_HH

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <iostream>
#include <Stntuple/loop/TStnModule.hh>
#include <Stntuple/loop/TStnAna.hh>
#include <Stntuple/obj/TStnHeaderBlock.hh>
#include <Stntuple/obj/TStnPhotonBlock.hh>
#include <Stntuple/obj/TStnJetBlock.hh>
#include <Stntuple/obj/TStnElectronBlock.hh>
#include <Stntuple/obj/TGenpBlock.hh>

#endif

class TPrintModule {

	protected:

	
	public:
		TPrintModule(const char* name="samsPrintModule", const char* title="samsPrintModule", TStnHeaderBlock* fHeaderBlock=NULL);
		TPrintModule(TStnHeaderBlock* fHeaderBlock);
		TPrintModule(TStnHeaderBlock*, TStnPhoton*);
		TPrintModule(TStnHeaderBlock*, TStnElectron*);
		TPrintModule(TStnHeaderBlock*, TStnJet*);
		TPrintModule(TStnHeaderBlock*, TGenParticle*);
		TPrintModule(TStnHeaderBlock*, TStnPhoton*,TStnJet*);
		~TPrintModule();

		void Reset();		//____ clear everything and re-init
		
		//setters
		void AddThisJet(TStnJet*);
		void AddThisPhoton(TStnPhoton*);
		void AddThisElectron(TStnElectron*);
		void AddThisHEPGparticle(TGenParticle*);
	
		//print methods
		void PhotonInfo(TStnPhoton*);
		void JetInfo(TStnJet*);
		void ElectronInfo(TStnElectron*);
		void HepgInfo(TGenParticle*);
		void Print(Option_t* opt="a");
		void Print(TStnPhotonBlock*);
		void Print(TStnElectronBlock*);
		void Print(TStnJetBlock*);
		void Print(TGenpBlock*);

	private:
		TObject *obj;
		bool photonheader, jetheader, electronheader,hepgheader;
		Int_t photoncounter,jetcounter,electroncounter,hepgcounter;
		void Init(TStnHeaderBlock*);
		int eventNumber, runNumber;
		std::vector<TStnJet*> jet_list;
		std::vector<TStnPhoton*> photon_list;
		std::vector<TStnElectron*> electron_list;
		std::vector<TGenParticle*> hepgpar_list;
	
};

#endif
