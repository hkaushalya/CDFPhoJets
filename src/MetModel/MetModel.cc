#include "samantha/MetModel/MetModel.hh"
#include "Stntuple/loop/TStnAna.hh"

ClassImp(MetModel)

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
MetModel::MetModel(const char* name, const char* title):
   TStnModule(name,title)
{

}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
int MetModel::BeginJob()
{
	//_____________________________________________________ register the data block
	RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlockClu04);

	//---------- need this for jet-EM matching
	RegisterDataBlock("CalDataBlock","TCalDataBlock",&fCalData);
	RegisterDataBlock("PhotonBlock","TStnPhotonBlock",&fPhotonBlock);
	RegisterDataBlock("ElectronBlock","TStnElectronBlock",&fElectronBlock);
	RegisterDataBlock("MetBlock","TStnMetBlock",&fMetBlock);		//added to replace Sasha's EventFiltermod 03-14-2008
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	//_____________________________________________________ book histograms

	BookHistograms();

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
int MetModel::BeginRun()
{
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
int MetModel::Event(int ientry)
{
	ClearAll();
	GetRawJets();
	GetEMobjects(EMobjects);
	
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
int MetModel::EndJob()
{
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void MetModel::ClearAll()
{
	Jets.rawJets.ClearAll();
	Jets.lev1Jets.ClearAll();
	Jets.lev4Jets.ClearAll();
	Jets.lev5Jets.ClearAll();
	Jets.lev6Jets.ClearAll();
	Jets.lev7Jets.ClearAll();
	Jets.metAddedJets.ClearAll();
	Jets.smearedJets.ClearAll();

}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void MetModel::BookHistograms()
{

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void MetModel::GetEMobjects(EMObjects_t& emobjs)
{

	/*photons.iType = EMCollection::iPHOTON;		//photons
	
	int Npho=initSpMod->GetSuperPhoSize();

	for (int i=0; i<Npho; i++)
	{
		if ( initSpMod->GetSuperPhoton(i)->IsTightPhoton() ||
				initSpMod->GetSuperPhoton(i)->IsLoosePhoton())
		{
			TLorentzVector _pho_raw = initSpMod->GetSuperPhoton(i)->GetRawVec();
			TLorentzVector _pho_cor = initSpMod->GetSuperPhoton(i)->GetCorVec();
			photons.ivBlockIndex.push_back(initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex());
			photons.tlvRawMomentum.push_back(_pho_raw);
			photons.tlvCorrMomentum.push_back(_pho_cor);
			double hadem=initSpMod->GetSuperPhoton(i)->GetHadEm();
			photons.fvEmFr.push_back(1.0/(hadem+1.0));
			photons.fvEtaDet.push_back(initSpMod->GetSuperPhoton(i)->GetDetEta());
			photons.fvXces.push_back(initSpMod->GetSuperPhoton(i)->GetXCes());
			photons.fvZces.push_back(initSpMod->GetSuperPhoton(i)->GetZCes());
		}
	}
*/
}



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void MetModel::GetElectrons(EMObjects_t& electrons)
{
/*
	electrons.iType = EMCollections::iELECTRON;		//electrons

	int _myNele=initSpMod->GetSuperPhoSize();
	for (int i=0; i<_myNele; i++)
	{
		if (initSpMod->GetSuperPhoton(i)->IsTightElectron()
				|| initSpMod->GetSuperPhoton(i)->IsLooseElectron()
				|| initSpMod->GetSuperPhoton(i)->IsStdLooseElectron())
		{

			TLorentzVector _ele_raw = initSpMod->GetSuperPhoton(i)->GetRawVec();
			TLorentzVector _ele_cor = initSpMod->GetSuperPhoton(i)->GetCorVec();
			electrons.ivBlockIndex.push_back(initSpMod->GetSuperPhoton(i)->GetElectronBlockIndex());
			electrons.tlvRawMomentum.push_back(_ele_raw);
			electrons.tlvCorrMomentum.push_back(_ele_cor);
			double hadem=initSpMod->GetSuperPhoton(i)->GetHadEm();
			electrons.fvEmFr.push_back(1.0/(hadem+1.0));
			electrons.fvEtaDet.push_back(initSpMod->GetSuperPhoton(i)->GetDetEta());
			electrons.fvEtaDet.push_back(initSpMod->GetSuperPhoton(i)->GetDetPhi());
			electrons.fvXces.push_back(initSpMod->GetSuperPhoton(i)->GetXCes());
			electrons.fvZces.push_back(initSpMod->GetSuperPhoton(i)->GetZCes());
		}
	}
*/
}




