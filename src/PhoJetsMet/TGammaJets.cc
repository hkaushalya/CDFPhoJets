/*
when I say "pho_", generally it means the leading photon (most of the time)

READ THIS EVERYDAY!
keep it simple and readable. forget about performance!
*/

#include <iostream>
#include <fstream>
#include <TMath.h>
#include <algorithm>
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnTriggerTable.hh"
#include "samantha/PhoJetsMet/TGammaJets.hh"
#include "samantha/PhoJetsMet/TGammaJetsInit.hh"
#include "samantha/PhoJetsMet/TMyJetFilterModule.hh"
#include "Stntuple/obj/TL3Em.hh"
//#include "Stntuple/obj/TL3Met.hh"
#include "Stntuple/photon/TPhotonUtil.hh"
#include "samantha/PhoJetsMet/TPrintModule.hh"
#include "samantha/PhoJetsMet/TMathModule.hh"
#include "Stntuple/obj/TGenParticle.hh"
//for jay
#include "Stntuple/data/TCalTower.hh"


double const kPI = TMath::Pi();
double const kTWOPI = 2.*kPI;
const static unsigned int kPHO_ID_BITS = 19; //total number photon ID bits (loose + tight)
const static unsigned int kPHO_LOOSE_ID_BITS = 8;    //number of Loose ID variables
const static double kPHO_ET_CUT = 30; 	//my leading photon. so i make the cut on 30GeV!

ClassImp(TGammaJets)

//_____________________________________________________________________________
TGammaJets::TGammaJets(const char* name, const char* title):
  TStnModule(name,title), qMc(false)
{
	std::cout << " CREATED A TGAMMAJETS MODULE." << std::endl;
}

//_____________________________________________________________________________
TGammaJets::~TGammaJets() {
}

//_____________________________________________________________________________
void TGammaJets::SaveHistograms() {
}

//_____________________________________________________________________________
void TGammaJets::BookHistograms() {

	DeleteHistograms();

	BookGeneralHistograms();
	BookInvariantMassHistograms();
	BookPhotonHistograms();
	BookPhotonIDcutsHistograms();
	//BookL3SummaryHistograms();
	//BookGenpStudyHistograms();
	//BookHEPGPhotonMatchingHistograms();
	BookMetHistograms();
	BookBackgroundHistograms();
	BookEmObjectRemovedFromJetsHistograms();
	BookJetHistograms();


}


//_____________________________________________________________________________
int TGammaJets::BeginJob()
{
				// register the data block, why?  to unpack them!

	fHeaderBlock 	= (TStnHeaderBlock*) RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fCprDataBlock 	= (TCprDataBlock*) 	RegisterDataBlock("CprDataBlock","TCprDataBlock");
	fJetBlock 		= (TStnJetBlock*)    RegisterDataBlock("JetBlock","TStnJetBlock");
	//fJetBlock7 		= (TStnJetBlock*)    RegisterDataBlock("PROD@JetCluModule-cone0.7","TStnJetBlock");
	fTrackBlock 	= (TStnTrackBlock*) 	RegisterDataBlock("TrackBlock","TStnTrackBlock");
	fPhotonBlock 	= (TStnPhotonBlock*) RegisterDataBlock("PhotonBlock","TStnPhotonBlock");
	fTriggerBlock 	= (TStnTriggerBlock*) RegisterDataBlock("TriggerBlock","TStnTriggerBlock");
	fElectronBlock = (TStnElectronBlock*) RegisterDataBlock("ElectronBlock","TStnElectronBlock");
	fVertexBlock 	= (TStnVertexBlock*) RegisterDataBlock("VertexBlock","TStnVertexBlock");
	fMetBlock 		= (TStnMetBlock*)  	RegisterDataBlock("MetBlock","TStnMetBlock");
	fL3SummaryBlock= (TL3SummaryBlock*) RegisterDataBlock("L3SummaryBlock","TL3SummaryBlock");
	fGenpBlock     = (TGenpBlock*) RegisterDataBlock("GenpBlock","TGenpBlock");

	fCalDataBlock = (TCalDataBlock*) RegisterDataBlock("CalDataBlock","TCalDataBlock");

	BookHistograms();

	// initialize vars
	counter.evtsWithoutAphoton = 0;
	counter.evtsWithoutAjet    = 0;
	counter.evtsPassTrigger    = 0;
	counter.evtsPass25Trigger    = 0;
	counter.evtsPass50Trigger    = 0;
	counter.evtsPass70Trigger    = 0;
	counter.evtsWithPhotonAnd2Jets = 0;
	counter.evtsPassLooseCuts  = 0;
	counter.evtsPassTightCuts  = 0;
	counter.sidebandevts       = 0;
	counter.evtsPass2JetsCut   = 0;
	counter.total_evts_b4_strip = 0;
	counter.total_evts_a4_strip = 0;
	counter.phos_match2_genp = 0;
	counter.phos_nomatch2_genp = 0;
	counter.HEPGmomPi0 = 0;
	counter.HEPGmomChargedPi = 0;
	counter.HEPGmomEta = 0;
	counter.HEPGmomOmega = 0;
	counter.HEPGmomOther = 0;
	counter.Wevents = 0;
	counter.Zevents = 0;
	counter.pho_2j_events = 0;
	counter.goodrunevts = 0;

	fakeRateSum = 0.;

  return 0;
}


//_____________________________________________________________________________
int TGammaJets::BeginRun()

{
	run = fHeaderBlock->RunNumber();
//	bool good = goodrun.Good(run);

//	if (!good) return 0;

	if(fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int TGammaJets::Event(int ientry)
{
	//std::cout << "\n\n>>>>>>>>>>>>> entry:"<< ientry << " IN TG EVENT LOOP\n";
	SetPassed(0);
	/*	
	TStnTriggerTable* ptab = (TStnTriggerTable*);
	DBman.Instance()->GetTable("TriggerTable");

	int nl3 = ptab->NL3Triggers();
	const TStnTrigger* trig;

	for (int i =0; i<nl3;i++) {
		trig = ptab->GetTrigger(3,i);
		std::cout << trig->Name() <<std::endl;
	}
	*/

	//_______________________________ remove the bad runs/sections
	int run = fHeaderBlock->RunNumber();
	evt = fHeaderBlock->EventNumber();
	sec = fHeaderBlock->SectionNumber();
	float instlum = fHeaderBlock->InstLum();
	
	//______________________________ REQUIRED METHODS
	Cleanup();

	if (IsWhatIWant()) {
		SetPassed(1);
		Background_ElectronFakeRate(ientry);
	}
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// do my event selection
// 1. loop to find photon (pass tight id with Etc > 30)
// 2. check i have two jets > 15 gev
// 3. do things with them
/*-------------------------------------------------------------------*/
bool TGammaJets::IsWhatIWant()
{
	bool foundPho = false;
	bool foundEle = false;
	bool found2Jets = false;
		//_talk to TGammaJetsInit and get my photon info	
	TGammaJetsInit* MyPhoton = (TGammaJetsInit*) ((TStnAna*) GetAna()->GetModule("GammaJetsInit")); // "GammaJetsInit" is a default name

	int Npho = MyPhoton->GetSuperPhotonSize();

	for (int i = 0; i < Npho; i++ ) {
		int TightID = MyPhoton->GetSuperPhotonTightID(i);
		int EleID = MyPhoton->GetSuperPhotonPhoLikeID(i);
		std::cout << "==i=" << i << "   TID= " << TightID << "\t  EID=" << EleID << std::endl;
		PhoTightID.push_back(TightID);
		PhoEleID.push_back(EleID);
		PhoVec.push_back(MyPhoton->GetSuperPhotonCorVec(i));
		if (!foundPho && (TightID == 0) ) {
			foundPho = true;
			leadPhoIndex = i;			// this is MY leading photon
			leadPhoVec = MyPhoton->GetSuperPhotonCorVec(i);
		}
	}
	
	TMyJetFilterModule* MyJet = (TMyJetFilterModule*) ((TStnAna*) GetAna()->GetModule("MyJetFilter")); // "GammaJetsInit" is a default name

	int nj = MyJet->GetMyJetNoEMobj_lev7size();
	MyJet->GetMyJetNoEMobj_lev7size();
	if ( nj > 2) {
		leadJet1Vec = *(MyJet->GetMyJetNoEMobj_lev7(0,0));   	// cone 0.4, leading jet
		leadJet2Vec = *(MyJet->GetMyJetNoEMobj_lev7(0,1));		// cone 0.4, 2nd leading jet

		double et1 = leadJet1Vec.Energy() * TMath::Sin(leadJet1Vec.Theta());
		double et2 = leadJet2Vec.Energy() * TMath::Sin(leadJet2Vec.Theta());
		assert( et1 > et2);
		if (et2 > 15) {
			found2Jets = true;
			counter.pho_2j_events++;
		}
		
	}
	return (foundPho && found2Jets);
}


/*-------------------------------------------------------------------*/
// for JAY
/*-------------------------------------------------------------------*/
/*{{{*/
void TGammaJets::JayTower(int ientry)
{
	if (fJetBlock->NJets() <4) return;
	fCalDataBlock->GetEntry(ientry);
	cout << "=============================== " << ientry   << " ========================= "<< endl;
	int Ntows = fCalDataBlock->NTowers();
	cout << ">>>>>>>>>>>>>>>> Cal Tower Summary. <<<<<<<<<<<<<<<<<<" << endl;
	cout << "\tNumber of HitTowers= " << Ntows <<"\n"<< endl;
	
	for (int i=0; i< Ntows; i++) {
		TCalTower *t = fCalDataBlock->Tower(i);
		double HadE = t->HadEnergy();
		double EmE 	= t->EmEnergy();
		double E    = t->Energy();
		double Et   = t->Et();
		double HadEt = HadE * fabs(sin(t->Theta()));
		double EmEt = EmE * fabs(sin(t->Theta()));
		int index_Phi = t->IPhi();
		int index_Eta = t->IEta();

		printf("%3i  %3i  %9.3f  %9.3f  %9.3f  %9.3f %9.3f  %9.3f\n", index_Eta,index_Phi, E, EmE, HadE, Et, EmEt, HadEt); 
/*		cout << i <<"  #Eta_i,Phi_i: " << index_Eta << "," << index_Phi << endl;
		cout << i <<"  #HadE,EmE   : " << HadE << "," << EmE << endl;
		cout << "\tHadEt = " << HadEt  << endl;
		cout << "\tEmEt  = " << EmEt  << endl;
		cout << "\tEt    = " << Et << endl;
*/
	}

	
	cout << "\n>>>>>>>>>>>>>>>> Offline Jet Summary. <<<<<<<<<<<<<<<<<<" << endl;
	cout << "\tNumber of Jets= " << fJetBlock->NJets() <<"\n"<< endl;
	
	for (int i=0; i< fJetBlock->NJets(); i++) {
		TStnJet* j= fJetBlock->Jet(i);
		double jeta = j->DetEta();
		double jphi = j->Phi();
		TLorentzVector lp = *(j->Momentum());
		double je   = lp.E();
		double jet = j->Et();
		double emfr = j->Emfr();

		printf("%9.3f  %9.3f  %9.3f  %9.3f  %9.3f \n", jeta,jphi, je, jet, emfr); 
		

	}
}
/*}}}*/


/*-------------------------------------------------------------------*/
// setup things common to each event  
/*-------------------------------------------------------------------*/
void TGammaJets::DoCommonThings(int ientry)
{
	
	SetPhotonList(ientry); 	//________ fill all the photon info to my SuperPhoton list
	SetJetList(ientry);		//________ need to do the same thing as photons 
	RemovePhotonsFromJetBlock();
	
	// set what i mostly used
	leading_pho  = GetLeadingPhoton();
	leading_jet1 = GetLeadingJet1();	// want to change these like the photon
	leading_jet2 = GetLeadingJet2();

	bool l  = (bool)(superpho_vec[0].LooseID);
	bool t  = (bool)(superpho_vec[0].TightID);
	bool pe = (bool)(superpho_vec[0].PhotonLikeID);

	if (l == 0) counter.evtsPassLooseCuts++;
	if (t == 0) counter.evtsPassTightCuts++;
	
	if (pe == 0) { 		//____ assume W, if a photon passed PhotonLikeID cuts

		//_______________ note:: my W,Z count will be lower since Photon like id cuts uses central cut.
		//_______________ so imay want to remove it!
		if (fMetBlock->Met(0) > 30) {
			counter.Wevents++;
		}

		if (superpho_vec.size() > 1) {	//______________ assume Z, if two photons, both pass Photon like IDs
			if (superpho_vec[1].PhotonLikeID == 1) {
				counter.Zevents++;
			}
		}

		//do the fake rate calculation
	}

	//________________ fill ID bin plot
	// bit = 0 if it pass the ID , bit == 1 if failed
	if ( l &&  t &&  pe ) GeneralPlot.pho_IDbins_pass->Fill(0); 	// failed all
	if (!l &&  t &&  pe ) GeneralPlot.pho_IDbins_pass->Fill(1); 	// pass l, fail t,pe
	if ( l && !t &&  pe ) GeneralPlot.pho_IDbins_pass->Fill(2); 	// fail 1, pass t, fail pe
	if ( l &&  t && !pe ) GeneralPlot.pho_IDbins_pass->Fill(3);		// fail 1, fail t, pass pe
	if (!l && !t &&  pe ) GeneralPlot.pho_IDbins_pass->Fill(4);  	// pass l,t, fail pe
	if (!l &&  t && !pe ) GeneralPlot.pho_IDbins_pass->Fill(5);		// pass l, fail t, pass pe
	if ( l && !t && !pe ) GeneralPlot.pho_IDbins_pass->Fill(6);		// fail l, pass t, pass pe
	if (!l && !t && !pe ) GeneralPlot.pho_IDbins_pass->Fill(7);		// pass all

	if ( l &&  t &&  pe ) GeneralPlot.pho_IDbins_fail->Fill(0); 	// failed all
	if ( l && !t && !pe ) GeneralPlot.pho_IDbins_fail->Fill(1); 	// fail l, pass t,pe
	if (!l &&  t && !pe ) GeneralPlot.pho_IDbins_fail->Fill(2); 	// pass 1, fail t, pass pe
	if (!l && !t &&  pe ) GeneralPlot.pho_IDbins_fail->Fill(3);		// pass 1, pass t, fail pe
	if ( l &&  t && !pe ) GeneralPlot.pho_IDbins_fail->Fill(4);  	// fail l,t, pass pe
	if ( l && !t &&  pe ) GeneralPlot.pho_IDbins_fail->Fill(5);		// fail l, pass t, fail pe
	if (!l &&  t &&  pe ) GeneralPlot.pho_IDbins_fail->Fill(6);		// pass l, fail t, fail pe
	if (!l && !t && !pe ) GeneralPlot.pho_IDbins_fail->Fill(7);		// pass all

	
}


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void TGammaJets::Cleanup()
{
	leading_pho = NULL;
	leading_jet1 = NULL;
	leading_jet2 = NULL;

	trig25 = false;
	trig50 = false;
	trig70 = false;

	superpho_vec.clear();
	superjet_vec.clear();

	leadPhoIndex = -99999;
	leadPhoVec.SetPxPyPzE(0,0,0,0);
	leadJet1Vec.SetPxPyPzE(0,0,0,0);
	leadJet2Vec.SetPxPyPzE(0,0,0,0);
	PhoVec.clear();
} //Cleanup


/*-------------------------------------------------------------------*/
// generate stuff needed by sasha's code
/*-------------------------------------------------------------------*/
bool TGammaJets::InitForSashaCode(int ientry)
{
/*	int Npho = fPhotonBlock->NPhotons();

	for (int i =0 ; i < Npho; i++) {
		TStnPhoton* pho = fPhotonBlock->Photon(i);
	
		if	(LoosePhotonCuts(pho)) {
			//Nphotons++;  //done
	
			TLorentzVector corvec;
			//corvec.SetPxPyPzE( pho->Momentum()->Px(),pho->Momentum()->Py(),
			//								pho->Momentum()->Pz(),pho->Momentum()->E());
			//phoCorrected.push_back(corvec);
			
			TLorentzVector rawvec;
			rawvec.SetPtEtaPhiM(pho->Et(),pho->Eta(),pho->Phi(),0.0);
			phoUncorrected.push_back(rawvec);
			
			//phoEtaDet.push_back(pho->DetEta());
			//phoHadEm.push_back(pho->HadEm());
			//phoCesX.push_back(pho->XCes());
			//phoCesZ.push_back(pho->ZCes());
			//phoIndex.push_back(i);

			std::cout << "Corr.E,PX,PY,PZ=" << corvec.E()<<"," << corvec.Px()<< "," << corvec.Py() <<","<< corvec.Pz() << std::endl;
			//std::cout << "RAW.E,PX,PY,PZ=" << rawvec.E() <<"," << rawvec.Px() <<","<< rawvec.Py() <<","<< rawvec.Pz() << std::endl;

		} // if
	} // for

	if (Nphotons>0) return true;
	else return false;
*/
} //InitForSashaCode


/*-------------------------------------------------------------------*/
// generator level study
/*-------------------------------------------------------------------*/
void TGammaJets::GenpBlockStudy(int ientry)
{
	if (fPhotonBlock->NPhotons() < 1) {

		int Nparticles = fGenpBlock->NParticles();
		float Et = 0;
		TGenParticle *hEtphoton;

		for (int i=0; i<Nparticles; i++) {
			TGenParticle *par = fGenpBlock->Particle(i);
		
			if (( par->GetStatusCode()==1) && par->IsPhoton()) {
			
				if ((par->Energy()*TMath::Sin(par->Theta())) > Et) {
					Et = par->Energy() * TMath::Sin(par->Theta());
					hEtphoton = par;
				}
			}
		} // for

		GenpPlot.Photon_Et->Fill(hEtphoton->Energy()*TMath::Sin(hEtphoton->Theta()));
		GenpPlot.PhotonsLostInTheDetector_Eta->Fill(hEtphoton->Eta());
		GenpPlot.PhotonsLostInTheDetector_Phi->Fill(hEtphoton->Phi());

	} //if

} //GenpBlockStudy


/*-------------------------------------------------------------------*/
// ID-cut analyses (individual and N-1 cuts)
/*-------------------------------------------------------------------*/
void TGammaJets::PhotonIDcutsAna(int ientry)
{
	int Detector   = leading_pho->Detector();
	float EtCorr 	= leading_pho->ECorr() * leading_pho->SinTheta();
	float XCes 		= leading_pho->XCes();
	float ZCes 		= leading_pho->ZCes();
	float HadEm 	= leading_pho->HadEm();
	float IsoEtCorr= leading_pho->EIso4(2);	//___ corrected for leakage and MI
	float TrkPt  	= leading_pho->Pt();     	//___ using max pt track in the cluster
	float TrkIso 	= leading_pho->SumPt4();

	// fill histos before cut/s
	PhoIDcutPlot.Detector_b->Fill(Detector);
	PhoIDcutPlot.EtCorr_b->Fill(EtCorr);
	PhoIDcutPlot.XCes_b->Fill(XCes);
	PhoIDcutPlot.ZCes_b->Fill(ZCes);
	PhoIDcutPlot.HadEm_b->Fill(HadEm);
	PhoIDcutPlot.IsoEtCorr_b->Fill(IsoEtCorr);
	PhoIDcutPlot.TrkPt_b->Fill(TrkPt);
	PhoIDcutPlot.TrkIso_b->Fill(TrkIso);

	bool passDetector = false;
	bool passEtCorr 	= false;
	bool passXCes 		= false;
	bool passZCes 		= false;
	bool passHadEm 	= false;
	bool passTrkPt 	= false;
	bool passTrkIso 	= false;
	bool passIsoEtCorr= false;

	for (int i = 0 ; i < kPHO_LOOSE_ID_BITS; ++i) {
	
		bool pass = true;
		
		if (i != 0) {
			if (Detector != 0) pass = false;
			else passDetector = true;
		}
		if (i != 1) {
			if (EtCorr < kPHO_ET_CUT) pass = false;
			else passEtCorr = true;
		}
		if (i != 2) {
			if (fabs(XCes) > 21) pass = false;
			else passXCes = true;
		}
		if (i != 3) {
			if ((fabs(ZCes) < 9) || (fabs(ZCes) > 230)) pass = false;
			else passZCes = true;
		}
		if (i != 4) {
			if (HadEm > 0.125) pass = false;
			else passHadEm = true;
		}
		if (i != 5) {
			if (TrkPt > 0.25*EtCorr) pass = false;
			else passTrkPt = true;
		}
		if (i != 6) {
			if (TrkIso > 5.0) pass = false;
			else passTrkIso = true;
		}
		if (i != 7) {
			if (EtCorr < 20) {
				if (IsoEtCorr > 0.15*EtCorr) pass = false;
			else passIsoEtCorr = true;
			} else {
				if (IsoEtCorr > 3.0*EtCorr) pass = false;
			else passIsoEtCorr = true;
			}
		}

		if (pass) {	//________________________________ fill N-1 cut plots
			if (i == 0)	{ PhoIDcutPlot.Detector_a_n_1->Fill(Detector); }
			else if (i == 1)	{ PhoIDcutPlot.EtCorr_a_n_1->Fill(EtCorr); }
			else if (i == 2)	{ PhoIDcutPlot.XCes_a_n_1->Fill(XCes);}
			else if (i == 3)	{ PhoIDcutPlot.ZCes_a_n_1->Fill(ZCes);}
			else if (i == 4)	{ PhoIDcutPlot.HadEm_a_n_1->Fill(HadEm);}
			else if (i == 5)	{ PhoIDcutPlot.TrkPt_a_n_1->Fill(TrkPt);}
			else if (i == 6)	{ PhoIDcutPlot.TrkIso_a_n_1->Fill(TrkIso);}
			else if (i == 7)	{ PhoIDcutPlot.IsoEtCorr_a_n_1->Fill(IsoEtCorr);}
			
		} //if
	} //for
	
		//std::cout << "\n" << ientry << "#\t leading photon\n\tEta=" << leading_pho->DetEta() << "\tPhi="
		//		<< leading_pho->Phi()  << "\tIsoEtCorr=" << IsoEtCorr << "\n";
	//if (TrkIso < 6.0) {
/*		if (LoosePhotonCuts(leading_pho) == 0) {
			if (MatchPho2GenLevel(leading_pho) > 0) {
				std::cout << "\t========== START ==========================\n";
				TPrintModule* mm = new TPrintModule(fHeaderBlock);
				mm->Print(fGenpBlock);
				std::cout << "\t=========== END ===========================\n";
				
				counter.phos_match2_genp++;
			} else {
				counter.phos_nomatch2_genp++;
			}
		}
*/	
	//}
	
	//_______________________________________________ fill individual id cuts
	if (passDetector) PhoIDcutPlot.Detector_a->Fill(Detector);
	if (passEtCorr) PhoIDcutPlot.EtCorr_a->Fill(EtCorr);
	if (passXCes) PhoIDcutPlot.XCes_a->Fill(XCes);
	if (passZCes) PhoIDcutPlot.ZCes_a->Fill(ZCes);
	if (passHadEm) PhoIDcutPlot.HadEm_a->Fill(HadEm);
	if (passTrkPt) PhoIDcutPlot.TrkPt_a->Fill(TrkPt);
	if (passTrkIso) PhoIDcutPlot.TrkIso_a->Fill(TrkIso);
	if (passIsoEtCorr) PhoIDcutPlot.IsoEtCorr_a->Fill(IsoEtCorr);

} // PhotonIDcutsAna



/*-------------------------------------------------------------------*/
/* match leading offline photon to a generator level                 */
/* return: number of macthes found.                                  */
/*-------------------------------------------------------------------*/
unsigned int TGammaJets::MatchPho2GenLevel(TStnPhoton* pho)
{
	if (!qMc) {
		return 0;   // false
	}
	unsigned int found_matches = 0;
	bool found_match = false;
	double pho_eta = pho->Eta();		//_ use event eta. not detEta
	double pho_phi = pho->Phi();
	double pho_EtCorr = pho->ECorr() * pho->SinTheta();
	
	int Nparticles = fGenpBlock->NParticles();
	if (Nparticles < 1) {

	}
	
	TPrintModule* tp = new TPrintModule(fHeaderBlock);
	tp->AddThisPhoton(pho);
	
	TGenParticle* closest_par;
	double closest_par_delR = 999.0;
	double closest_par_Et_ratio = 0.;
	
	for (int i = 0 ; i < Nparticles; i++) {
		TGenParticle* par = fGenpBlock->Particle(i);
		int im = par->GetFirstMother();

		if(im>=0) {	//_________________________________im=-1 means, no mother for imcoming particles
			TGenParticle* mom = fGenpBlock->Particle(im);
			if (mom != 0) {
				int mom_id = mom->GetPdgCode();
				int mom_stable = mom->GetStatusCode();
				int mom_number = mom->Number();
				//double mom_charge = mom->Charge();  // give an error in TGenParticle.hh line 56??
				int par_id = par->GetPdgCode();
				int par_stable = par->GetStatusCode();
				int par_number = par->Number();
				//double par_charge = par->Charge();
				
				if ( par->GetStatusCode()==1 ) {
					
						float gen_eta = par->Eta();
						float gen_phi = par->Phi();
						TMathModule math;
						//found_match = math.MatchEtaPhi(pho_phi,pho_eta,gen_phi,gen_eta,0.1); 		 // match with this precision
						
						//if (found_match) {
							double delEta = math.GetDelEta(pho_eta,gen_eta);
							double delPhi = math.GetDelPhi(pho_phi,gen_phi);
							double delR   = math.GetDelR(pho_phi,pho_eta,gen_phi,gen_eta);
							double par_Et = par->Energy() * TMath::Sin(par->Theta());
							//HEPGPhoMatchPlot.DelEta->Fill(delEta);	
							//HEPGPhoMatchPlot.DelPhi->Fill(delPhi);	
							//HEPGPhoMatchPlot.DelR->Fill(delR);
							//HEPGPhoMatchPlot.HEPGObj_Et->Fill(par->Energy() * TMath::Sin(par->Theta()));
							//HEPGPhoMatchPlot.Et_ratio->Fill(pho_EtCorr / (par->Energy() * TMath::Sin(par->Theta())) );
							//HEPGPhoMatchPlot.HEPGObj_Eta->Fill(par->Eta());
							double match_et_ratio = par_Et/pho_EtCorr;
							
						//if (delR < closest_par_delR && (match_et_ratio > 0.4 )) {
						if (delR < closest_par_delR) {
							tp->AddThisHEPGparticle(par);
							closest_par = par;
							closest_par_delR = delR;
							closest_par_Et_ratio = match_et_ratio;
							++found_matches;
							//cout << "\t" << i <<"# closest par_delR = "<< closest_par_delR <<endl;

							if (par_Et < 15) {		
								if (mom->IsPi0()) { 
									counter.HEPGmomPi0++;
									HEPGPhoMatchPlot.matched_Object->Fill(0);
								} else if (mom->IsElectron()) { 
									counter.HEPGmomElectron++;
									HEPGPhoMatchPlot.matched_Object->Fill(1);
								} else if (mom->IsMuon()) { 
									counter.HEPGmomMuon++;
									HEPGPhoMatchPlot.matched_Object->Fill(2);
								} else if (mom->IsTau()) { 
									counter.HEPGmomTau++;
									HEPGPhoMatchPlot.matched_Object->Fill(3);
								} else if ((mom_id == 211) || (mom_id == -211)) {	//_______________________ Pi+/-
									counter.HEPGmomChargedPi++;
									HEPGPhoMatchPlot.matched_Object->Fill(4);
								} else if (mom_id == 221) {	//_______________________ eta particle
									counter.HEPGmomEta++;
									HEPGPhoMatchPlot.matched_Object->Fill(5);
								} else if (mom_id == 223) {	//_______________________ omega particle
									counter.HEPGmomOmega++;
									HEPGPhoMatchPlot.matched_Object->Fill(6);
								} else {
									counter.HEPGmomOther++;		//_______________________ some other
									//cout << "\t" << i<<"#"<< par->GetPdgCode() << endl;
									HEPGPhoMatchPlot.matched_Object->Fill(7);
/*{{{*/
/*								cout << "---> matching mom info\n";
								cout << i <<"#\tPDGCode="<<mom_id << "\tStable?="<< mom_stable << "\tUniqNumber="<<mom_number << endl;
								cout << "--------> Genp dump <-------------------------------\n";
								for (int i = 0 ; i < Nparticles; i++) {
									TGenParticle* par_temp = fGenpBlock->Particle(i);
									int imtemp 				= par_temp->GetFirstMother();
									if (imtemp>=0) {
										TGenParticle* mom_temp= fGenpBlock->Particle(imtemp);
										if (mom_temp !=0) {
											int mom_id_temp 		= mom_temp->GetPdgCode();
											int mom_stable_temp 	= mom_temp->GetStatusCode();
											int mom_number_temp 	= mom_temp->Number();
											int par_id_temp 		= par_temp->GetPdgCode();
											int par_stable_temp 	= par_temp->GetStatusCode();
											int par_number_temp 	= par_temp->Number();
											cout << "---> mom info\n";
											cout << i <<"#\tPDGCode="<<mom_id_temp << "\tStable?="<< mom_stable_temp << "\tUniqNumber="<< mom_number_temp << endl;
											cout << "\t---> particle info\n";
											cout << i <<"#\tPDGCode="<<par_id_temp << "\tStable?="<< par_stable_temp << "\tUniqNumber="<< par_number_temp << endl;
										}
									}
								}
*/
							//}
							/*
							cout << "---> mom info\n";
							cout << "\tPDGCode="<<mom_id << "\tStable?="<< mom_stable << "\tUniqNumber="<<mom_number << "\tCharge="<<mom_charge<<endl;
							cout << "-----> daughter info\n";
							cout << "\tPDGCode="<<par_id << "\tStable?="<< par_stable << "\tUniqNumber="<<par_number << "\tCharge="<<par_charge<<endl;
							std::cout << " Generator level" << i <<"th photon: Eta=" << gen_eta << "\tPhi="
								<< gen_phi << "\tmatched? " << found_match <<"\n";
								math.Print();
							cout <<"\t///////////////////////////////////////////////////////////\n";
							*/
/*}}}*/
							}
						} // if 
						
					}

				} // if

			} // if
		} // if
			
	} // for

	if (found_matches > 0 && closest_par_delR > 0.4) {
		std::cout << "\t========== START ==========================" << "delR="<<closest_par_delR << "\n";
		tp->AddThisHEPGparticle(closest_par);
		tp->Print("ph");
		TPrintModule* pp = new TPrintModule(fHeaderBlock);
		pp->Print(fGenpBlock);
		std::cout << "\t========== END   ==========================\n";
	}
	//if (closest_par != NULL && found_matches > 0) {
	//if (closest_par_delR < 999.0) {
	if (found_matches>0) {
		HEPGPhoMatchPlot.DelR->Fill(closest_par_delR);
		HEPGPhoMatchPlot.HEPGObj_Et->Fill(closest_par->Energy() * TMath::Sin(closest_par->Theta()));
		HEPGPhoMatchPlot.Et_ratio->Fill(closest_par_Et_ratio);
	}
	
	//cout << "\t matches =" << found_matches <<endl;
	HEPGPhoMatchPlot.Nmatches->Fill(found_matches);
	return found_matches;
}


/*-------------------------------------------------------------------*/
// Gamma+2 Jets Selection
/*-------------------------------------------------------------------*/
void TGammaJets::GammaJetsSelection(int ientry)
{

	if (superpho_vec[0].TightID == 0) {
		if ( (leading_jet1->Et() > 15) && (leading_jet2->Et() > 15) ) {
		}
	} 

} // GammJetsSelection

/*-------------------------------------------------------------------*/
// Calculate 3 body (photon+2jets) invaraint mass 
/*-------------------------------------------------------------------*/
double TGammaJets::Get3BodyMass(TLorentzVector phovec, TLorentzVector j1vec,
						TLorentzVector j2vec)
{

	TLorentzVector sum3 = phovec + j1vec + j2vec;
	return sum3.M();

}


/*-------------------------------------------------------------------*/
//  photon+2jets  propeties
/*-------------------------------------------------------------------*/

	//MassPlot.j1_j2->Fill(j1j2sum.M());

/*
	//fraction of energy
	GeneralPlot.j1_to_photon_energy_ratio->Fill(j1vec.E()/phovec.E());
	GeneralPlot.j2_to_photon_energy_ratio->Fill(j2vec.E()/phovec.E());
	GeneralPlot.j2_to_j1_energy_ratio->Fill(j2vec.E()/j1vec.E());	
	//fraction of momentum to energy
	GeneralPlot.j1_momentum_to_energy_ratio->Fill(j1vec.P()/j1vec.E());
	GeneralPlot.j2_momentum_to_energy_ratio->Fill(j2vec.P()/j2vec.E());

	//fraction of momentum w.r.t photon
	GeneralPlot.j1_to_photon_momentum_ratio->Fill(j1vec.P()/phovec.P());
	GeneralPlot.j2_to_photon_momentum_ratio->Fill(j1vec.P()/phovec.P());
	GeneralPlot.j2_to_j1_momentum_ratio->Fill(j2vec.P()/j1vec.P());

	//Et plots
	PhotonPlot.EtCorr->Fill(plist[0]->ECorr() * plist[0]->SinTheta());
	Jet1Plot.EtRaw->Fill(jlist[0]->Et());
	Jet2Plot.EtRaw->Fill(jlist[1]->Et());
*/
/*
	//Phi plots
	float pPhi = pho->Phi();
	float j1Phi = j1->Phi();
	float j2Phi = j2->Phi();
	float pj1Phi = pj1sum.Phi();
	float j2vecPhi = j2vec.Phi();

	if (fabs(pPhi - j1Phi) > kPI) {
		if (pPhi> kPI) pPhi = kTWOPI - pPhi;
		if (j1Phi> kPI) j1Phi = kTWOPI - j1Phi;
		//fHist.fH[35]->Fill(fabs(pPhi+j1Phi));
	} else {
		//fHist.fH[35]->Fill(fabs(pPhi-j1Phi));
	}
	

	if (fabs(pPhi - j2Phi) > kPI) {
		if (pPhi> kPI) pPhi = kTWOPI - pPhi;
		if (j2Phi> kPI) j2Phi = kTWOPI - j2Phi;
		//fHist.fH[36]->Fill(fabs(pPhi+j2Phi));
	} else {
		//fHist.fH[36]->Fill(fabs(pPhi-j2Phi));
	}


	if (fabs(j1Phi - j2Phi) > kPI) {
		if (j1Phi> kPI) j1Phi = kTWOPI - j1Phi;
		if (j2Phi> kPI) j2Phi = kTWOPI - j2Phi;
		//fHist.fH[37]->Fill(fabs(j1Phi+j2Phi));
	} else {
		//fHist.fH[37]->Fill(fabs(j1Phi-j2Phi));
	}

	
	if (fabs(pj1Phi - j2vecPhi) > kPI) {
		if (pj1Phi> kPI) pj1Phi = kTWOPI - pj1Phi;
		if (j2vecPhi> kPI) j2vecPhi = kTWOPI - j2vecPhi;
		//fHist.fH[38]->Fill(fabs(pj1Phi+j2vecPhi));
	} else {
		//fHist.fH[38]->Fill(fabs(pj1Phi-j2vecPhi));
	}
*/




/*-------------------------------------------------------------------*/
// plot Missing Et
/*-------------------------------------------------------------------*/
void TGammaJets::MissingEtPlots()
{
	MetPlot.Met_0->Fill(fMetBlock->Met(0));
	MetPlot.MetPhi_0->Fill(fMetBlock->MetPhi(0));
	MetPlot.MetX_0->Fill(fMetBlock->MetX(0));
	MetPlot.MetY_0->Fill(fMetBlock->MetY(0));
	MetPlot.Sumet_0->Fill(fMetBlock->Sumet(0));
	MetPlot.Sumet_1->Fill(fMetBlock->Sumet(1));
	MetPlot.Sumet_2->Fill(fMetBlock->Sumet(2));
	MetPlot.Metsig->Fill(fMetBlock->MetSig());
	MetPlot.Z0->Fill(fMetBlock->Z0());
}


/*-------------------------------------------------------------------*/
// plot leading photon and leading 2 jets Phi and Eta Seperation
/*-------------------------------------------------------------------*/
void TGammaJets::PlotPhotonJetsPhi(int ientry)
{

	if (leading_pho == NULL) return;
	if (leading_jet1 == NULL) return;
	if (leading_jet2 == NULL) return;

	//leading photon
	float pPhi = leading_pho->Phi();
	float pEta = leading_pho->DetEta();
	//leading jets > 15 GeV
	if ( leading_jet1->Et() < 20 || leading_jet2->Et() < 20 ) return;
	float j1Phi = leading_jet1->Phi();
	float j1Eta = leading_jet1->DetEta();
	float j2Phi = leading_jet2->Phi();
	float j2Eta = leading_jet2->DetEta();

	float delPhi1 =0;   //photon and leading jet Phi seperation
	float delPhi2 =0;   //photon and 2nd leading jet Phi seperation
	float delPhi3 =0;   // 2 leading jets Phi seperation

	if ( fabs(pPhi - j1Phi) >  kPI) {
		float j1Phi_temp = 0, pPhi_temp =0;
		if (pPhi > kPI) pPhi_temp = kTWOPI - pPhi;
		if (j1Phi > kPI) j1Phi_temp = kTWOPI - j1Phi;
		delPhi1 = fabs(pPhi_temp + j1Phi_temp);
	} else {
		delPhi1 = fabs(pPhi - j1Phi);
	}

	if ( fabs(pPhi - j2Phi) >  kPI) {
		float j2Phi_temp = 0, pPhi_temp =0;
		if (pPhi > kPI) pPhi_temp = kTWOPI - pPhi;
		if (j2Phi > kPI) j2Phi_temp  = kTWOPI - j2Phi;
		delPhi2 = fabs(pPhi_temp + j2Phi_temp);
	} else {
		delPhi2 = fabs(pPhi - j2Phi);
	}

	if ( fabs(j1Phi - j2Phi) >  kPI) {
		float j2Phi_temp = 0, j1Phi_temp =0;
		if (j1Phi > kPI) j1Phi_temp = kTWOPI - j1Phi;
		if (j2Phi > kPI) j2Phi_temp = kTWOPI - j2Phi;
		delPhi3 = fabs(j1Phi_temp + j2Phi_temp);
	} else {
		delPhi3 = fabs(j1Phi - j2Phi);
	}

	assert(delPhi1 < kPI);
	assert(delPhi2 < kPI);
	assert(delPhi3 < kPI);

	////fHist.fH[40]->Fill(delPhi1);
	////fHist.fH[41]->Fill(delPhi2);
	////fHist.fH[42]->Fill(delPhi3);
//	std::cout << "************* " << ientry <<std::endl;
//	std::cout << "p-j1=" << delPhi1 << "  p-j2=" < delPhi2 << "  j1-j2=" << delPhi3  << " SUM=" << delPhi1+delPhi2+delPhi3 <<std::endl;
}

/*}}}*/


/*-------------------------------------------------------------------*/
// to find a cut value to remove my leading photon from the jet Block
// I set it to 0.15
/*-------------------------------------------------------------------*/

void TGammaJets::PlotJetPhotonSeperation(int ientry)
/*{{{*/
{

	if (superpho_vec.size() < 1) return;
	int Njets = superjet_vec.size();
	if (Njets< 1) return;

	std::vector<float> delRlist (Njets);
	std::vector<TStnJet*> jlist (Njets);
	jlist = GetSortedJetList();
	
	for (int i=0; i < Njets; i++) {
		TMathModule* mymathmod = new TMathModule;
		delRlist[i] = mymathmod->GetDelR(leading_pho,superjet_vec[i].jet);
	}

	//std::sort(delRlist.begin(),delRlist.end());

	std::cout << "******** " << ientry << "\n";
	for (int i = 0; i < delRlist.size() ; i++) {
		std::cout << i << " #  " << delRlist[i] ;
	}

	int loopupto = Njets-1;	
	while (loopupto !=0) {
		for (int i=0; i<loopupto;i++) {
			if (delRlist[i] > delRlist[i+1]) {
				float delR = delRlist[i+1];
				delRlist[i+1] = delRlist[i];
				delRlist[i] = delR;
				
				TStnJet* temp = jlist[i+1];
				jlist[i+1] = jlist[i];
				jlist[i] = temp;
			}
		}
		loopupto--;
	}

	//std::cout << ">>>>>>>>>> " << ientry << "\tdelR[0] = " << delRlist[0] <<"\t"<< " delR[1] = " << delRlist[1] <<"\n";
	if (delRlist[1] < 0.15) {
		std::cout << "\n====================== START SUMMARY ======================\n";
		std::cout << "    delR \t Et\n";
		for (int i = 0; i < delRlist.size() ; i++) {
			printf("%2i  %6.3f  %6.2f\n", i, delRlist[i], jlist[i]->Et());
		}
		std::cout << "\n====================== END   SUMMARY ======================\n";
	}
}
/*}}}*/


/*-------------------------------------------------------------------*/
//	METHOD 1:
// Removes the leading Photon from the JetBlock by
// looking at Eta/Phi Space
// make sure to pass in the sorted objects lists!
/*-------------------------------------------------------------------*/

void TGammaJets::RemovePhotonsFromJetBlock()
{
	TPrintModule* tp = new TPrintModule(fHeaderBlock);

	int np = superpho_vec.size();
	int nj = superjet_vec.size();
	bool foundmatch = false;
	for (int i = 0; i < superpho_vec.size() ; i++) {
		if (foundmatch) break;
		TMathModule mymathmod;
		if (superjet_vec.size() < 1) break;
		for (int j =0; j < superjet_vec.size(); j++) {
			float delR = mymathmod.GetDelR(superpho_vec[i].pho, superjet_vec[j].jet);
			float delEta = mymathmod.GetDelEta();
			float delPhi = mymathmod.GetDelPhi();
			if (delR < 0.1) {
				foundmatch = true;
				//cout << "\t,matched=" << i <<"," << j <<endl;
				//mymathmod->Print();
				EmJetMatchPlot.DelR->Fill(delR);	
				EmJetMatchPlot.DelEta->Fill(delEta);	
				EmJetMatchPlot.DelPhi->Fill(delPhi);
				EmJetMatchPlot.Et_ratio->Fill(superpho_vec[i].Etraw/superjet_vec[j].Et);
				
				superjet_vec.erase(superjet_vec.begin()+j);
				break;
			}
		}
	}
	int nja = superjet_vec.size();

	int sum = (nja+np -nj);
	EmJetMatchPlot.NOmatches->Fill(foundmatch);	
/*	if (sum !=0) {
		tp->Print(fPhotonBlock);
		tp->Print(fJetBlock);
		tp->Reset();
		for (int i=0; i < superjet_vec.size() ; i++) {
			tp->AddThisJet(superjet_vec[i].jet);
		}
		tp->Print("j");
	}
*/
}

/*-------------------------------------------------------------------*/
// plot photon loose/tight cuts cumulative plot 
/*-------------------------------------------------------------------*/
void TGammaJets::PhotonCutsCumulativePlot()
{
	GeneralPlot.photon_loose_and_tight_cuts_cumm->Fill(0); 	//_____ for passing the trigger

	for (int i = 0 ; i < kPHO_ID_BITS ;i++) {
		int bit = -9;
		if (i < kPHO_LOOSE_ID_BITS) {
			bit = (superpho_vec[0].LooseID >> i) & 0x1;
		} else {
			bit = (superpho_vec[0].TightID >> i) & 0x1;
		}
		if ( !bit ) { 	//_________________________________________ bit was set if a cut was failed
			GeneralPlot.photon_loose_and_tight_cuts_cumm->Fill(i+1);
		} else break;	//_________________________________________ since i am making cumulative plot
	}
}



/*-------------------------------------------------------------------*/
// plot photon variables before and after cuts 
/*-------------------------------------------------------------------*/
void TGammaJets::PhotonIDPlots(int ientry)
{
	TStnPhoton* pho = GetLeadingPhoton();
	if (pho == NULL) return;

	int 	Detector 	= pho->Detector();
	float EtCorr 		= pho->ECorr() * pho->SinTheta(); 	//corrected energy this is what i want to cut on
	float XCes 			= pho->XCes();
	float ZCes 			= pho->ZCes();
	float HadEm 		= pho->HadEm();
	float IsoEtCorr	= pho->EIso4(2);  						// corrected for leakage and MI
	float Chi2Mean		= pho->Chi2Mean();  						// (strip+wire)/2
	int   N3d			= pho->N3d();
	float TrkPt  		= pho->Pt();     							//using max pt track in the cluster
	float TrkIso 		= pho->SumPt4();  					
	float CesWireE2	= pho->CesWireE2();
	float CesStripE2	= pho->CesStripE2();

	PhotonPlot.Detector_b->Fill(Detector);
	PhotonPlot.EtCorr_b->Fill(EtCorr);
	PhotonPlot.XCes_b->Fill(XCes);
	PhotonPlot.ZCes_b->Fill(ZCes);
	PhotonPlot.HadEm_b->Fill(HadEm);
	PhotonPlot.IsoEtCorr_b->Fill(IsoEtCorr);
	PhotonPlot.Chi2Mean_b->Fill(Chi2Mean);
	PhotonPlot.N3d_b->Fill(N3d);
	PhotonPlot.TrkPt_b->Fill(TrkPt);
	PhotonPlot.TrkIso_b->Fill(TrkIso);
	PhotonPlot.Ces2Wire_b->Fill(CesWireE2);
	PhotonPlot.Ces2Strip_b->Fill(CesStripE2);

	if (!superpho_vec[0].LooseID) {
		PhotonPlot.Detector_aL->Fill(Detector);
		PhotonPlot.EtCorr_aL->Fill(EtCorr);
		PhotonPlot.XCes_aL->Fill(XCes);
		PhotonPlot.ZCes_aL->Fill(ZCes);
		PhotonPlot.HadEm_aL->Fill(HadEm);
		PhotonPlot.IsoEtCorr_aL->Fill(IsoEtCorr);
		PhotonPlot.Chi2Mean_aL->Fill(Chi2Mean);
		PhotonPlot.N3d_aL->Fill(N3d);
		PhotonPlot.TrkPt_aL->Fill(TrkPt);
		PhotonPlot.TrkIso_aL->Fill(TrkIso);
		PhotonPlot.Ces2Wire_aL->Fill(CesWireE2);
		PhotonPlot.Ces2Strip_aL->Fill(CesStripE2);

		if (!superpho_vec[0].TightID) {
			PhotonPlot.Detector_aT->Fill(Detector);
			PhotonPlot.EtCorr_aT->Fill(EtCorr);
			PhotonPlot.XCes_aT->Fill(XCes);
			PhotonPlot.ZCes_aT->Fill(ZCes);
			PhotonPlot.HadEm_aT->Fill(HadEm);
			PhotonPlot.IsoEtCorr_aT->Fill(IsoEtCorr);
			PhotonPlot.Chi2Mean_aT->Fill(Chi2Mean);
			PhotonPlot.N3d_aT->Fill(N3d);
			PhotonPlot.TrkPt_aT->Fill(TrkPt);
			PhotonPlot.TrkIso_aT->Fill(TrkIso);
			PhotonPlot.Ces2Wire_aT->Fill(CesWireE2);
			PhotonPlot.Ces2Strip_aT->Fill(CesStripE2);
		}
	}
	
}

/*-------------------------------------------------------------------*/
// Plot variables of two jets	
/*-------------------------------------------------------------------*/
void TGammaJets::JetPlots(TStnJet *j1, TStnJet *j2)
{
		float	DetEta1 	= j1->DetEta();	
		float Phi1   	= j1->Phi();
		float Et1		= j1->Et();
		float Emfr1		= j1->Emfr();
		float EOverP1	= j1->EOverP();
	
		Jet1Plot.DetEta->Fill(DetEta1);
		Jet1Plot.Phi->Fill(Phi1);
		Jet1Plot.Etraw->Fill(Et1);
		Jet1Plot.Emfr->Fill(Emfr1);
		Jet1Plot.EOverP->Fill(EOverP1);

		float	DetEta2 	= j2->DetEta();
		float Phi2   	= j2->Phi();
		float Et2		= j2->Et();
		float Emfr2		= j2->Emfr();
		float EOverP2	= j2->EOverP();
	
		Jet2Plot.DetEta->Fill(DetEta2);
		Jet2Plot.Phi->Fill(Phi2);
		Jet2Plot.Etraw->Fill(Et2);
		Jet2Plot.Emfr->Fill(Emfr2);
		Jet2Plot.EOverP->Fill(EOverP2);

		TMathModule mm;
		double delPhi= mm.GetDelPhi(j1,j2);
		double delEta= mm.GetDelEta(j1,j2);
		double delR  = mm.GetDelR(j1,j2);
		JetPlot.DelPhi->Fill(delPhi);
		JetPlot.DelEta->Fill(delEta);
		JetPlot.DelR->Fill(delR);
		JetPlot.EtRatio->Fill(Et2/Et1);
}



/*-------------------------------------------------------------------*/
//
/*-------------------------------------------------------------------*/
void TGammaJets::TriggerEfficiency(int ientry)
{/*{{{*/
	int   index = 0;
	float HighestEtCorr = 0; 
	TStnPhoton* pho;
	TStnPhoton* mypho;

	int Npho = fPhotonBlock->NPhotons();
	if (Npho < 1) return;

	for (int j=0; j < Npho; j++) {
	
		pho = fPhotonBlock->Photon(j);
	   float EtCorr = pho->ECorr() * pho->SinTheta();	
		if (EtCorr > HighestEtCorr) {
			HighestEtCorr = EtCorr;
			mypho = pho;
		}
	}
	
	if (((mypho->ECorr()*mypho->SinTheta()) < 25) && (mypho->Detector()==0)) {
		////fHist.fH[1]->Fill(mypho->EtCorr());
		////fHist.fH[3]->Fill(mypho->VertexZ());
	}
}


/*}}}*/



/*-------------------------------------------------------------------*/
//	to study L3 SUMMARY
/*-------------------------------------------------------------------*/
void TGammaJets::L3Summary(int ientry)
{/*{{{*/


	int NEm 			= fL3SummaryBlock->NEm();
	int NMuon 		= fL3SummaryBlock->NMuon();
	int NJet4 		= fL3SummaryBlock->NJet4();  // number of L3 jets with cone size 0.4
	int NJet7 		= fL3SummaryBlock->NJet7();  // number of L3 jets with cone size 0.6
	int NCotTrack	= fL3SummaryBlock->NCotTrack();
	int NSvtTrack 	= fL3SummaryBlock->NSvtTrack();
	int NMet 		= fL3SummaryBlock->NMet();   //??? what is NMet
	int NTau 		= fL3SummaryBlock->NTau();
   //std::cout << NEm <<","<<NMuon<<","<< NCotTrack<<","<<NSvtTrack<<"\n";
	
	L3SummaryPlot.NEm->Fill(NEm);
	L3SummaryPlot.NMuon->Fill(NMuon);
	L3SummaryPlot.NJet4->Fill(NJet4);
	L3SummaryPlot.NJet7->Fill(NJet7);
	L3SummaryPlot.NCotTrack->Fill(NCotTrack);
	L3SummaryPlot.NSvtTrack->Fill(NSvtTrack);
	L3SummaryPlot.NMet->Fill(NMet);
	L3SummaryPlot.NTau->Fill(NTau);
	if (NEm>0) {
		L3SummaryPlot.NTau_to_NEm->Fill(NTau/(NEm*1.0));
		L3SummaryPlot.NMuon_to_NEm->Fill(NMuon/(NEm*1.0));
	}
	if (NJet7>0) {
		L3SummaryPlot.NJet4_to_NJet7->Fill(NJet4/(NJet7*1.0));
	}
	if (NJet4>0) {
		L3SummaryPlot.NJet7_to_NJet4->Fill(NJet7/(NJet4*1.0));
	}
	if (NCotTrack>0) {
		L3SummaryPlot.NSvtTrack_to_NCotTrack->Fill(NSvtTrack/(NCotTrack*1.0));
		L3SummaryPlot.NCotTrack_vs_NEm->Fill(NEm,NCotTrack);
	}
	if (NSvtTrack>0) {
		L3SummaryPlot.NCotTrack_to_NSvtTrack->Fill(NCotTrack/(NSvtTrack*1.0));
		L3SummaryPlot.NSvtTrack_vs_NEm->Fill(NEm,NSvtTrack);
	}
	//std::cout << " OUT\n";
/*	if (NEm<1) return;

	float EmObjEt = 0.0;
	int index = 0;
	for (int i=0;i<NEm;i++) {
		TL3Em* L3EmObject = fL3SummaryBlock->Em(i);
		if (L3EmObject->Phet() > EmObjEt) {
			EmObjEt = L3EmObject->Phet();
			index = i;
		}
	}	
	TL3Em* L3Em = fL3SummaryBlock->Em(index);
*/
}

/*}}}*/



/*-------------------------------------------------------------------*/
// Remove the overlap of Photons and Electrons	
/*-------------------------------------------------------------------*/
void TGammaJets::RemoveOverlap(int ientry)
{/*{{{*/


	int phoMask = 0xffff;  
	int ind = 0xff;  //??????
	
	int Nele = fElectronBlock->NElectrons();
	int Npho = fPhotonBlock->NPhotons();
	

}/*}}}*/



/*-------------------------------------------------------------------*/
//	Sorts the Photon Block according to EtCorr, from highest to lowest
// & fill the photon info to my SuperPhoton list
/*-------------------------------------------------------------------*/
void TGammaJets::SetPhotonList(int ientry)
{
	int Npho = fPhotonBlock->NPhotons();
	std::vector<TStnPhoton*> photon;
	std::vector<int> index;

	for (int i = 0; i < Npho ; i++) {
		photon.push_back(fPhotonBlock->Photon(i));
		index.push_back(i);
	}

	if (Npho>1) {   //________________________ dont really need this since i check for this at the beginnig!
		int loopupto = Npho-1;	
		while (loopupto !=0) {
			for (int i=0; i<loopupto;i++) {
				if ( (photon[i]->ECorr()*photon[i]->SinTheta()) < 
						(photon[i+1]->ECorr()*photon[i+1]->SinTheta()) ) {

					TStnPhoton *temp = photon[i];
					photon[i] = photon[i+1];
					photon[i+1] = temp;

					int t = index[i];
					index[i] = index[i+1];
					index[i+1] = t;
					
				} //_____________________________ if
			} //________________________________ for
			loopupto--;
		} //___________________________________ while
	} //______________________________________ if
	
	double zvx  = 0;
	if (fVertexBlock->GetBestVertex(12,1) != NULL) {
		zvx = fVertexBlock->GetBestVertex(12,1)->Z();
	}
	
	TStnEvent* event = GetEvent();

	//TPrintModule* pp = new TPrintModule(fHeaderBlock);
	//std::cout << " xxxxxxxxxxxxxxxxxxx Start dump \n";
	//pp->Print(fPhotonBlock);
	//pp->Print(fElectronBlock);
	
	for (int i=0; i < photon.size(); i++) {
		SuperPhoton_t sp;
		sp.index    = index[i];
		sp.pho 		= photon[i];
		sp.Detector = photon[i]->Detector();
		sp.DetEta 	= photon[i]->DetEta();
		sp.Phi 		= photon[i]->Phi();
		sp.Etraw 	= photon[i]->Et();
		sp.ECorr 	= photon[i]->ECorr();
		sp.EtCorr 	= photon[i]->ECorr() * photon[i]->SinTheta();
		sp.XCes 		= photon[i]->XCes();
		sp.ZCes 		= photon[i]->ZCes();
		sp.HadEm 	= photon[i]->HadEm();
		sp.Chi2Mean = photon[i]->Chi2Mean();
		sp.N3d 		= photon[i]->N3d();
		sp.IsoEtCorr= photon[i]->EIso4(2);
		sp.TrkPt 	= photon[i]->Pt();
		sp.TrkIso 	= photon[i]->SumPt4();
		sp.CesWireE2= photon[i]->CesWireE2();
		sp.CesStripE2= photon[i]->CesStripE2();
		sp.LooseID	= LoosePhotonCuts(photon[i]);
		sp.TightID	= TightPhotonCuts(photon[i]);
		sp.PhotonLikeID	= !(CEMPhoEleIDcut(photon[i], event, zvx));  // i am using the negative logic
//		std::cout << "\n\t\t:photonloop:  " <<i<<"# PhotonLikeID  =" <<sp.PhotonLikeID <<std::endl;
//		std::cout << "\t\t:photonloop:  " <<i<<"# PhotonLooseID =" <<sp.LooseID <<std::endl;
//		std::cout << "\t\t:photonloop:  " <<i<<"# PhotonTightID =" <<sp.TightID <<std::endl;
		sp.BeamHalo = false;				// for now????
		sp.Cosmic   = false;				// for now???
		
		sp.corvec.SetPxPyPzE( photon[i]->Momentum()->Px(),photon[i]->Momentum()->Py(),
											photon[i]->Momentum()->Pz(),photon[i]->Momentum()->E());
		sp.rawvec.SetPtEtaPhiM(photon[i]->Et(),photon[i]->Eta(),photon[i]->Phi(),0.0);

		superpho_vec.push_back(sp);
	}
//	std::cout << " xxxxxxxxxxxxxxxxxxx End dump \n";

	if (photon.size() > 1)	assert(photon[0]->ECorr()*photon[0]->SinTheta() 
										> photon[1]->ECorr()*photon[1]->SinTheta());

}



/*===================================================================*/
//	Sorts the Jet Block according to Et, from Higest to Lowest,
// & fill the jet list
/*===================================================================*/
void TGammaJets::SetJetList(int ientry)
{
	int Njets = fJetBlock->NJets();
	std::vector<TStnJet*> jlist;
	std::vector<int> index;
	
	for (int i = 0; i < Njets ; i++) {
		jlist.push_back(fJetBlock->Jet(i));
		index.push_back(i);
	}

	if (Njets>1) {
		int loopupto = Njets-1;	
		while (loopupto !=0) {
			for (int i=0; i<loopupto;i++) {
				if (jlist[i]->Et() < jlist[i+1]->Et()) {
					TStnJet *temp = jlist[i];
					jlist[i] = jlist[i+1];
					jlist[i+1] = temp;

					int t = index[i];
					index[i] = index[i+1];
					index[i+1] = t;
				}
			}
			loopupto--;
		}
	}
	
	for (int i=0; i< jlist.size(); i++){
		SuperJet_t sj;
		sj.jet = jlist[i];
		sj.index = index[i];
		sj.Et = jlist[i]->Et();
		superjet_vec.push_back(sj);
	}
}


/*-------------------------------------------------------------------*/
/* for cut on vetices: require 1 vertex now. can remove lots of      */
/* background by doing so.															*/
/*-------------------------------------------------------------------*/
int TGammaJets::Class12Vertices()
{
	int nvtx = 0;
	double zvx  = 0;
	
	for (int ivtx = 0; ivtx < fVertexBlock->NVertices(); ++ivtx) {
			TStnVertex* vert = fVertexBlock->Vertex(ivtx);
			if (vert->VClass() >= 12) ++nvtx;
	}
	return nvtx;
}

/*-------------------------------------------------------------------*/
/* z-position cut on vertex 														*/
/*-------------------------------------------------------------------*/
double TGammaJets::BestVertex_Z()
{
	//double sumpt = 0;
	//TStnVertex* bestvert;
	double zvx =0;
	if (fVertexBlock->GetBestVertex(12,1) != NULL) {
		zvx = fVertexBlock->GetBestVertex(12,1)->Z();
	}


	/*for ( int ivtx = 0; ivtx < fVertexBlock->NVertices(); ++ivtx) {
		TStnVertex* ver = fVertexBlock->Vertex(ivtx);
		if (ver->SumPt() > sumpt) {
			bestvert = ver;
		}
	}
	if (bestvert) { std::cout << "found vertex\n"; }
	else { std::cout <<"no vertex\n"; }
	if (bestvert != NULL) return bestvert->Z();
	else return 0.0;
	*/
	return zvx;
}


/*-------------------------------------------------------------------*/
//	Tight Photon ID CUTS
/*-------------------------------------------------------------------*/
unsigned int TGammaJets::TightPhotonCuts(TStnPhoton *pho)
{
	unsigned int IDWord = 0x0;  // now, it does not matter if this is called many times
									// but, the counting of events passing/failing cuts will be wrong
									// if called many times.

	int 	detector = pho->Detector();
	float ECorr 		= pho->ECorr(); //total corrected energy // this is what i want to cut on
	float EtCorr 		= ECorr * pho->SinTheta(); //total corrected energy // this is waht i want to cut on
	float XCes 			= pho->XCes();
	float ZCes 			= pho->ZCes();
	float HadEm 		= pho->HadEm();
	float Chi2Mean 	= pho->Chi2Mean();  			// (strip+wire)/2
	int   N3d			= pho->N3d();
	float IsoEtCorr	= pho->EIso4(2);				//Corrected for leakage and MI
	float TrkPt  		= pho->Pt();     				//using max pt track in the cluster
	float TrkIso 		= pho->SumPt4();  			
	float CesWireE2	= pho->CesWireE2();
	float CesStripE2	= pho->CesStripE2();

	if (detector != 0						) 	IDWord |= kCentral_T;
	if (EtCorr < 30 						) IDWord |= kEtCorr7_T; // has no effect! looser cut
	if (fabs(XCes) > 21					) IDWord |= kXCes_T;
	if (fabs(ZCes) < 9 ||  fabs(ZCes) > 230) IDWord |= kZCes_T;  // same as loose cuts
	if ( !((HadEm < 0.125) || (HadEm < (0.055 + 0.00045 * ECorr))) )	IDWord |= kHadEm_T;
	if (Chi2Mean > 20						) IDWord |= kChi2Mean_T;
	if (N3d > 1							   ) IDWord |= kN3d_T;
	if (TrkPt > (1+0.005*EtCorr)		) IDWord |= kTrkPt_T;
	if (TrkIso > (2.0+0.005*EtCorr)	) IDWord |= kTrkIso_T;

	if (EtCorr < 20) {
		if (IsoEtCorr > 0.1*EtCorr) IDWord |= kIsoEtCorr_T;
	} else {
		if (IsoEtCorr > (2.0+0.02*(EtCorr-20.0)) ) IDWord |= kIsoEtCorr_T;
	}

	if (EtCorr<18) {
		if ( (CesWireE2 * pho->SinTheta()) > (0.14*EtCorr)) IDWord |= kCesWireE2_T;
		if ( (CesStripE2 * pho->SinTheta()) > (0.14*EtCorr)) IDWord |= kCesStripE2_T;
	} else {
		if ( (CesWireE2 * pho->SinTheta()) > (2.4+0.01*EtCorr)) IDWord |= kCesWireE2_T;
		if ( (CesStripE2 * pho->SinTheta()) > (2.4+0.01*EtCorr)) IDWord |= kCesStripE2_T;
	}

	return IDWord;

} // TightPhotonCuts


/*-------------------------------------------------------------------*/
//	Loose Photon ID CUTS
/*-------------------------------------------------------------------*/
unsigned int TGammaJets::LoosePhotonCuts(TStnPhoton *pho)
{
	unsigned int IDWord = 0x0;
	int 	detector = pho->Detector();
	float EtCorr 	= pho->ECorr() * pho->SinTheta();	// corrected energy this is what i want to cut on
	float XCes 		= pho->XCes();
	float ZCes 		= pho->ZCes();
	float HadEm 	= pho->HadEm();
	float IsoEtCorr= pho->EIso4(2);	//__________________ corrected for leakage and MI
	float TrkPt  	= pho->Pt();    	//__________________ using max pt track in the cluster
	float TrkIso 	= pho->SumPt4();
			
	if (detector != 0				) 	IDWord |= kCentral_L;
	if (EtCorr < 30				) 	IDWord |= kEtCorr30_L;
	if (fabs(XCes) > 21			) 	IDWord |= kXCes_L;
	if ((fabs(ZCes) < 9) || (fabs(ZCes) > 230)) IDWord |= kZCes_L;
	if (HadEm > 0.125				) 	IDWord |= kHadEm_L;
	if (TrkPt > 0.25*EtCorr		) 	IDWord |= kTrkPt_L;
	if (TrkIso > 5.0				) 	IDWord |= kTrkIso_L;

	if (EtCorr < 20) {
		if (IsoEtCorr > 0.15*EtCorr) IDWord |= kIsoEtCorr_L;
	} else {
		if (IsoEtCorr > 3.0*EtCorr) IDWord |= kIsoEtCorr_L;
	}

	return IDWord;

} //LoosePhotonCuts


//	ALL SETTERS GOES BELOW THIS POINT
/*-------------------------------------------------------------------*/


/*===================================================================*/
//	ALL BACKGROUND STUDIES GOES HERE
/*===================================================================*/

/*-------------------------------------------------------------------*/
void TGammaJets::Background_W()
{
	if (superpho_vec[0].PhotonLikeID == 1) { 		//____ assume W, if a photon passed PhotonLikeID cuts
		if (fMetBlock->Met(0) > 30) {
				//found a W. now do something!
			
			double Met_0 = fMetBlock->Met(0);
			BackgroundPlot.Wele_EtCorr->Fill(superpho_vec[0].EtCorr);
			BackgroundPlot.Wele_HadEm->Fill(superpho_vec[0].HadEm);
			BackgroundPlot.Wele_IsoEtCorr->Fill(superpho_vec[0].IsoEtCorr);
			BackgroundPlot.Wele_TrkPt->Fill(superpho_vec[0].TrkPt);
			BackgroundPlot.W_Met->Fill(Met_0);
			BackgroundPlot.W_eleEt_to_Met_ratio->Fill(superpho_vec[0].EtCorr/Met_0);
			BackgroundPlot.W_Sumet_0->Fill(fMetBlock->Sumet(0));
			BackgroundPlot.W_Sumet_1_Htc->Fill(fMetBlock->Sumet(1));
			BackgroundPlot.W_Sumet_2_Sumetjet->Fill(fMetBlock->Sumet(2));
			
			TLorentzVector pho_vec = *(superpho_vec[0].pho->Momentum());
			TLorentzVector phoT_vec;	//transverse component
			phoT_vec.SetPxPyPzE(pho_vec.Px(),pho_vec.Py(),0,pho_vec.E() * superpho_vec[0].pho->SinTheta());
			TLorentzVector met_vec;
			met_vec.SetPxPyPzE(Met_0 * TMath::Cos(fMetBlock->MetPhi(0)),Met_0 * TMath::Sin(fMetBlock->MetPhi(0)),0,Met_0);

			TLorentzVector sum = phoT_vec + met_vec;
			BackgroundPlot.W_Mass->Fill(sum.M());
			BackgroundPlot.W_Pt->Fill(sum.Perp());
			
			

		}
	}
}


/*-------------------------------------------------------------------*/
void TGammaJets::Background_Z()
{

	if (superpho_vec.size() > 1) {	//______________ Z, if two photons, both pass Photon like IDs
		if ( (superpho_vec[0].PhotonLikeID && superpho_vec[1].PhotonLikeID) == 1 ) {
			//found a Z. now do something!	
			
			double Met_0 = fMetBlock->Met(0);
			BackgroundPlot.Zele1_EtCorr->Fill(superpho_vec[0].EtCorr);
			BackgroundPlot.Zele2_EtCorr->Fill(superpho_vec[1].EtCorr);
			BackgroundPlot.Zele1_HadEm->Fill(superpho_vec[0].HadEm);
			BackgroundPlot.Zele1_IsoEtCorr->Fill(superpho_vec[0].IsoEtCorr);
			BackgroundPlot.Zele1_TrkPt->Fill(superpho_vec[0].TrkPt);
			BackgroundPlot.Zele2_TrkPt->Fill(superpho_vec[1].TrkPt);
			BackgroundPlot.Z_Met->Fill(Met_0);
			BackgroundPlot.Z_Sumet_0->Fill(fMetBlock->Sumet(0));
			BackgroundPlot.Z_Sumet_1_Htc->Fill(fMetBlock->Sumet(1));
			BackgroundPlot.Z_Sumet_2_Sumetjet->Fill(fMetBlock->Sumet(2));
			
			TLorentzVector pho1_vec = *(superpho_vec[0].pho->Momentum());
			TLorentzVector pho2_vec = *(superpho_vec[1].pho->Momentum());
			TLorentzVector sum = pho1_vec + pho2_vec;
			BackgroundPlot.Z_Mass->Fill(sum.M());
			BackgroundPlot.Z_Pt->Fill(sum.Perp());
			
			TMathModule math;
			BackgroundPlot.Z_eleDelPhi->Fill(math.GetDelPhi(superpho_vec[0].Phi,superpho_vec[1].Phi));
			BackgroundPlot.Z_eleDelEta->Fill(fabs(superpho_vec[0].DetEta - superpho_vec[1].DetEta));
			BackgroundPlot.Z_eleEtratio->Fill(superpho_vec[1].EtCorr/superpho_vec[0].EtCorr);

		}
	}

}



/*####################################################################*/

/*------- ElectronFakeRate: do f.r. and put on top of Mass plot ------*/
void TGammaJets::Background_ElectronFakeRate(int ientry)
{
	//2 methods
	//1: apply photon-like cuts
	//2: apply std. electron cuts
	//eiko: for photon Et>40, phoenix rejection is better

	//what i want to do
	//1. select my events with photon Et>30, 2 jets >15 GeV
	//2. photon must pass std. photon cuts
	//3. electron must pass photon-like id cuts
	//   if only one electron, do nothing
	//4. remove conversion electrons using std. conversion filter
	//5. (optional) phoenix-to-photon rejections
	//6. get fake rates for each electron and fill hito with that weight
	
	//1. select my events 
	MassPlot.photon_j1_j2->Fill(Get3BodyMass(leadPhoVec, leadJet1Vec, leadJet2Vec));
	if (PhoVec.size() < 2) return;
	
	double fr = 0; // fake rate
	double fr_minus_sigma=0, fr_plus_sigma=0;
	int goodrun_bit = 0;
	int phnx_bit = 1;  // phoenix rejection for now
	
	for (int i=0; i < PhoVec.size(); i++) {
		std::cout << "\t i=" << i <<  "\t" << PhoEleID[i] << std::endl;	
		if (PhoEleID[i] == 0) {	
		
			fr = 0; // fake rate
			fr_minus_sigma=0, fr_plus_sigma=0;
		
			double EtCorr = PhoVec[i].Energy() * TMath::Sin(PhoVec[i].Theta());
			FakeRate(EtCorr, goodrun_bit, phnx_bit, fr,fr_minus_sigma, fr_plus_sigma);
			std::cout << fr << "\t"<< fr_minus_sigma << "\t" <<fr_plus_sigma << "\n";
			fakeRateSum += fr;
			double mass = Get3BodyMass(PhoVec[i], leadJet1Vec, leadJet2Vec);
			MassPlot.photon_j1_j2_fake->Fill(mass,fr);
		}
	}

} // Background_W


/*****************  PHOTON LIKE ID CUT FUNCTION ************************/
//_____ function returns 1 if electrons passes cuts
//      requires zvx of highest Pt class12 vertex
int TGammaJets::CEMPhoEleIDcut(TStnPhoton* Pho, TStnEvent* event, double zvx)
{
	int passcode=1;

	//______________________________________ only CEM photons with Et>7 GeV
	//if ((Pho->Detector())!=0 || (Pho->Etc())<7.0) return 0;
	if ((Pho->Detector())!=0 || (Pho->Etc())<30.0) {
		//std::cout << "PHO-LIKE: detector||Et<30 failed." <<std::endl;
		return 0;
	}

	//______________________________________ HADEM cut using 3 towers
	float cutMin_HADEM=0.0;
	float cutMax_HADEM=0.055+0.00045*(Pho->Momentum()->E());
	if ((Pho->HadEm())<cutMin_HADEM || (Pho->HadEm())>cutMax_HADEM){
		//std::cout << "PHO-LIKE: hadem failed." <<std::endl;
		return 0;
	}

	//______________________________________ CES Chi^2 cut 
	if ((Pho->Chi2())<0.0 || (Pho->Chi2())>20.0) {
		//std::cout << "PHO-LIKE: chi2 failed" <<std::endl;
		return 0;
	}

	//______________________________________ N3D cut
	if ((Pho->N3d())==0) {
		//std::cout << "PHO-LIKE: N3d==0 failed." <<std::endl;
		return 0; // it is not an electron if N3d=0
	}

	//______________________________________ second N3D cut
	if ((Pho->N3d())>2) {
		//std::cout << "PHO-LIKE: N3d > 2 failed." <<std::endl;
		return 0;  // no more than 2 tracks (1st--ele,2nd--to match cuts for pho)
	}

	//______________________________________ cut on 2nd max Pt track in cluster if N3D=2
	if ((Pho->N3d())==2) {
		float trkPtcut_max=1.0+0.005*(Pho->Etc());

		if ((Pho->Pt2())>trkPtcut_max) {
		//std::cout << "PHO-LIKE: if (N3d==2) " <<std::endl;
			return 0;
		}
	}

	TStnElectron* Ele=TPhotonUtil::MatchElectron(event,Pho); // get matching electron

	//______________________________________ E/p cut for electron
	if (Ele==NULL) {
		//std::cout << "PHO-LIKE: NO ELECTRON MATCHED!." <<std::endl;
		return 0;
	} //else std:cout << "\t\t FOUND A MATCHING ELECTRON: with tracks="<< Ele->TrackNumber() << "\tZ0=" <<Ele->Z0()<< "\n";
	
	if (Ele->TrackNumber()<0) {
		//std::cout << "PHO-LIKE: no track." <<std::endl;
		return 0; // there should be a track
	}
	if(fabs(Ele->Z0()-zvx)>3.0) {
		//std::cout << "PHO-LIKE: Z seperation faile." <<std::endl;
		return 0; // only events with ele from best vertex
	}
	
	int ele_trk = Ele->TrackNumber();
	TStnTrack* Trk = fTrackBlock->Track(ele_trk);
	
	if ((Trk->Algorithm())==2) {
		if ((Ele->TrackBcPt())<50.0 && 
	 			(Ele->EOverP()<0.8 || Ele->EOverP()>1.2)) return 0;
	} else {
		if ((Ele->TrackPt())<50.0 && 
	 			(Ele->EOverP()<0.8 || Ele->EOverP()>1.2)) return 0;
	}

	//______________________________________ CalIso4 cut
	float cutMax_CalIso=0.0;

	if ((Pho->Etc()) < 20.0) cutMax_CalIso = 0.1 * (Pho->Etc());
	else cutMax_CalIso = 2.0 + 0.02 * ( Pho->Etc() - 20.0 );      

	if ((Pho->EIso4(2)) > cutMax_CalIso) return 0; 

	//______________________________________ TrkIso4 cut
	float trkIso=Pho->SumPt4() - Pho->Pt();

	if (trkIso < 0.0) trkIso = Pho->SumPt4(); 

	if (trkIso < (trkIso > (2.0 + 0.005 * (Pho->Etc()) ) ) ) return 0;

	//______________________________________ Energy of 2nd CES cluster (Wire & Strip)
	float _Et2ndCES=0.0;
	if ((Pho->CesStripE2()) > (Pho->CesWireE2())) _Et2ndCES = (Pho->CesStripE2()) * (Pho->SinTheta());
	else _Et2ndCES = (Pho->CesWireE2()) * (Pho->SinTheta());

	float cutMax_2ndCes = 0.0;
	if ((Pho->Etc()) < 18.0) cutMax_2ndCes = 0.14 * (Pho->Etc());	  
	else cutMax_2ndCes = 2.4 + 0.01 * (Pho->Etc());	  	  

	if (_Et2ndCES > cutMax_2ndCes) return 0;

	//_______________________________________ Fiducial cuts
	if (fabs(Pho->XCes()) > 21.0) return 0;
	if (fabs(Pho->ZCes()) < 9.0 || fabs(Pho->ZCes()) > 230.0) return 0;

	return passcode;
} // CEMPhoEleIDcut



/*********** FAKE RATE FUNCTION **************************************/

//- main routine to calculate fake rate
//_________________________________ returns fake rate & systematics
//
// 
//_________________________________ Main function _______________________
//__________ input:  eleEt= Pho->Etc(); (where Pho must be a TStnPhoton object)
//		     goodrun_bit=0 --- for goodrun_v13_pho_00.txt good run list
//		     goodrun_bit=1 --- for goodrun_v13_pho_lepphx_00.txt good run list		     
//		     phnx_bit=0 --- if you DON'T require photon-to-phoenix rejection
//		     phnx_bit=1 --- if you DO require photon-to-phoenix rejection
//____________________________________________________________________________________________________
//__________ output: wght_d=fake rate; wght_m=(fake rate) - (1 sigma); wght_p=(fake rate) + (1 sigma) 

void TGammaJets::FakeRate(double eleEt, int goodrun_bit, int phnx_bit, double &wght_d, double &wght_m, double &wght_p) {

  wght_d=0.0;
  wght_m=0.0;
  wght_p=0.0;

  if(goodrun_bit==0 && phnx_bit==0) // with phnx, all runs
    {
      phnx_par[0]=0.0; phnx_par[1]=0.0; phnx_par[2]=0.0; 
      phnx_parerr[0]=0.0; phnx_parerr[1]=0.0; phnx_parerr[2]=0.0; 
      phnx_corr_coeff[0][1]=phnx_corr_coeff[1][0]=0.0;
      phnx_corr_coeff[0][2]=phnx_corr_coeff[2][0]=0.0;
      phnx_corr_coeff[1][2]=phnx_corr_coeff[2][1]=0.0;
      phnx_corr_coeff[0][0]=phnx_corr_coeff[0][0]=1.0;
      phnx_corr_coeff[1][1]=phnx_corr_coeff[1][1]=1.0;
      phnx_corr_coeff[2][2]=phnx_corr_coeff[2][2]=1.0;
      fakerate_par[0]=-3.089; fakerate_par[1]=-0.0304; fakerate_par[2]=0.0059; 
      fakerate_parerr[0]=0.067; fakerate_parerr[1]=0.0055; fakerate_parerr[2]=0.0025;
      fakerate_par[3]=1.09;      // data/MC scale factor (using 40<Et<50)
      fakerate_parerr[3]=0.07; 
      fakerate_corr_coeff[0][0]=1.0; // diagonal covariance matrix correlation coefficients
      fakerate_corr_coeff[1][1]=1.0; 
      fakerate_corr_coeff[2][2]=1.0; 
      fakerate_corr_coeff[3][3]=1.0; 
      fakerate_corr_coeff[0][1]=fakerate_corr_coeff[1][0]=-0.514;
      fakerate_corr_coeff[0][2]=fakerate_corr_coeff[2][0]=0.237;
      fakerate_corr_coeff[1][2]=fakerate_corr_coeff[2][1]=-0.951;
      fakerate_corr_coeff[0][3]=fakerate_corr_coeff[3][0]=0.0;
      fakerate_corr_coeff[1][3]=fakerate_corr_coeff[3][1]=0.0;
      fakerate_corr_coeff[2][3]=fakerate_corr_coeff[3][2]=0.0;
    }
  if(goodrun_bit==1 && phnx_bit==0) // with phnx, good Silicon
    {
      phnx_par[0]=0.0; phnx_par[1]=0.0; phnx_par[2]=0.0; 
      phnx_parerr[0]=0.0; phnx_parerr[1]=0.0; phnx_parerr[2]=0.0; 
      phnx_corr_coeff[0][1]=phnx_corr_coeff[1][0]=0.0;
      phnx_corr_coeff[0][2]=phnx_corr_coeff[2][0]=0.0;
      phnx_corr_coeff[1][2]=phnx_corr_coeff[2][1]=0.0;
      phnx_corr_coeff[0][0]=phnx_corr_coeff[0][0]=1.0;
      phnx_corr_coeff[1][1]=phnx_corr_coeff[1][1]=1.0;
      phnx_corr_coeff[2][2]=phnx_corr_coeff[2][2]=1.0;
      fakerate_par[0]=-3.086; fakerate_par[1]=-0.0270; fakerate_par[2]=0.0040; // if v11_lepphx
      fakerate_parerr[0]=0.063; fakerate_parerr[1]=0.0054; fakerate_parerr[2]=0.0031;
      fakerate_par[3]=1.03; 
      fakerate_parerr[3]=0.08; 
      fakerate_corr_coeff[0][0]=1.0; // diagonal covariance matrix correlation coefficients
      fakerate_corr_coeff[1][1]=1.0; 
      fakerate_corr_coeff[2][2]=1.0; 
      fakerate_corr_coeff[3][3]=1.0; 
      fakerate_corr_coeff[0][1]=fakerate_corr_coeff[1][0]=-0.286;
      fakerate_corr_coeff[0][2]=fakerate_corr_coeff[2][0]=0.0;
      fakerate_corr_coeff[1][2]=fakerate_corr_coeff[2][1]=-0.955;
      fakerate_corr_coeff[0][3]=fakerate_corr_coeff[3][0]=0.0;
      fakerate_corr_coeff[1][3]=fakerate_corr_coeff[3][1]=0.0;
      fakerate_corr_coeff[2][3]=fakerate_corr_coeff[3][2]=0.0;
    }
  if(goodrun_bit==1 && phnx_bit==0) // no phnx, all runs
    {
      phnx_par[0]=0.394; phnx_par[1]=0.106; phnx_par[2]=24.54; 
      phnx_parerr[0]=0.004; phnx_parerr[1]=0.006; phnx_parerr[2]=0.32; 
      phnx_corr_coeff[0][1]=phnx_corr_coeff[1][0]=-0.482;
      phnx_corr_coeff[0][2]=phnx_corr_coeff[2][0]=0.450;
      phnx_corr_coeff[1][2]=phnx_corr_coeff[2][1]=-0.265;
      phnx_corr_coeff[0][0]=phnx_corr_coeff[0][0]=1.0;
      phnx_corr_coeff[1][1]=phnx_corr_coeff[1][1]=1.0;
      phnx_corr_coeff[2][2]=phnx_corr_coeff[2][2]=1.0;
      fakerate_par[0]=-3.089; fakerate_par[1]=-0.0304; fakerate_par[2]=0.0059; 
      fakerate_parerr[0]=0.067; fakerate_parerr[1]=0.0055; fakerate_parerr[2]=0.0025;
      fakerate_par[3]=1.42;      // data/MC scale factor (using 40<Et<50)
      fakerate_parerr[3]=0.15; 
      fakerate_corr_coeff[0][0]=1.0; // diagonal covariance matrix correlation coefficients
      fakerate_corr_coeff[1][1]=1.0; 
      fakerate_corr_coeff[2][2]=1.0; 
      fakerate_corr_coeff[3][3]=1.0; 
      fakerate_corr_coeff[0][1]=fakerate_corr_coeff[1][0]=-0.514;
      fakerate_corr_coeff[0][2]=fakerate_corr_coeff[2][0]=0.237;
      fakerate_corr_coeff[1][2]=fakerate_corr_coeff[2][1]=-0.951;
      fakerate_corr_coeff[0][3]=fakerate_corr_coeff[3][0]=0.0;
      fakerate_corr_coeff[1][3]=fakerate_corr_coeff[3][1]=0.0;
      fakerate_corr_coeff[2][3]=fakerate_corr_coeff[3][2]=0.0;
    }
  if(goodrun_bit==1 && phnx_bit==1) // no phnx, good Silicon
    {
      phnx_par[0]=0.421; phnx_par[1]=0.104; phnx_par[2]=24.61; 
      phnx_parerr[0]=0.004; phnx_parerr[1]=0.006; phnx_parerr[2]=0.32; 
      phnx_corr_coeff[0][1]=phnx_corr_coeff[1][0]=-0.487;
      phnx_corr_coeff[0][2]=phnx_corr_coeff[2][0]=0.426;
      phnx_corr_coeff[1][2]=phnx_corr_coeff[2][1]=-0.185;
      phnx_corr_coeff[0][0]=phnx_corr_coeff[0][0]=1.0;
      phnx_corr_coeff[1][1]=phnx_corr_coeff[1][1]=1.0;
      phnx_corr_coeff[2][2]=phnx_corr_coeff[2][2]=1.0;
      fakerate_par[0]=-3.086; fakerate_par[1]=-0.0270; fakerate_par[2]=0.0040; // if v11_lepphx
      fakerate_parerr[0]=0.063; fakerate_parerr[1]=0.0054; fakerate_parerr[2]=0.0031;
      fakerate_par[3]=1.19; 
      fakerate_parerr[3]=0.18; 
      fakerate_corr_coeff[0][0]=1.0; // diagonal covariance matrix correlation coefficients
      fakerate_corr_coeff[1][1]=1.0; 
      fakerate_corr_coeff[2][2]=1.0; 
      fakerate_corr_coeff[3][3]=1.0; 
      fakerate_corr_coeff[0][1]=fakerate_corr_coeff[1][0]=-0.286;
      fakerate_corr_coeff[0][2]=fakerate_corr_coeff[2][0]=0.0;
      fakerate_corr_coeff[1][2]=fakerate_corr_coeff[2][1]=-0.955;
      fakerate_corr_coeff[0][3]=fakerate_corr_coeff[3][0]=0.0;
      fakerate_corr_coeff[1][3]=fakerate_corr_coeff[3][1]=0.0;
      fakerate_corr_coeff[2][3]=fakerate_corr_coeff[3][2]=0.0;
    }

  if(goodrun_bit<0 || goodrun_bit>1 || phnx_bit<0 || phnx_bit>1) return;

  if(eleEt>=13.0)
    {
      double term1=FakeRateFuncErr(eleEt);
      double multp=(1.0-PhnxEffFunc(eleEt));
      double term2=PhnxEffFuncErr(eleEt)*FakeRateFunc(eleEt);
      double my_err=ShiftCorrFunc(eleEt)*sqrt(term1*term1*multp*multp+term2*term2);
      wght_d=ShiftCorrFunc(eleEt)*FakeRateFunc(eleEt)*(1.0-PhnxEffFunc(eleEt));
      wght_m=wght_d-my_err;
      wght_p=wght_d+my_err;
      if(wght_d<0.0) wght_d=0.0;
      if(wght_m<0.0) wght_m=0.0;
      if(wght_p<0.0) wght_p=0.0;
    }
  return;
  
} // FakeRate

//____________ correction for Et(det)/Et(genp) difference
double TGammaJets::ShiftCorrFunc(double eleEt) {

  double arg1=0.0048+exp(-3.031-0.027*eleEt);
  double arg2=0.0056+exp(-3.015-0.030*eleEt);
  double value=arg1/arg2;
  return value;

}

//__________________________________ phoenix efficiency function
double TGammaJets::PhnxEffFunc(double eleEt) {

double arg=phnx_par[1]*(phnx_par[2]-eleEt);
  double value=phnx_par[0]*TMath::Erfc(arg);
  if(value<0.0 || value>1.0) value=0.0;

return value;

}

//__________________________________ uncertainty on phoenix efficiency
double TGammaJets::PhnxEffFuncErr(double eleEt) {

double arg=phnx_par[1]*(phnx_par[2]-eleEt);
  double erfc_drv=2.0*exp(-1.0*arg*arg)/TMath::Pi();
  double term1=phnx_parerr[0]*TMath::Erfc(arg);
  double term2=phnx_parerr[1]*phnx_par[0]*(phnx_par[2]-eleEt)*erfc_drv;
  double term3=phnx_parerr[2]*phnx_par[0]*phnx_par[1]*erfc_drv;
  double term4=2.0*phnx_corr_coeff[0][1]*phnx_parerr[0]*phnx_parerr[1]*TMath::Erfc(arg)
    *erfc_drv*phnx_par[0]*(phnx_par[2]-eleEt);
  double term5=2.0*phnx_corr_coeff[0][2]*phnx_parerr[0]*phnx_parerr[2]*TMath::Erfc(arg)
    *erfc_drv*phnx_par[0]*phnx_par[1];
  double term6=2.0*phnx_corr_coeff[1][2]*phnx_parerr[1]*phnx_parerr[2]
    *arg*erfc_drv*erfc_drv*phnx_par[0]*phnx_par[0];
  double term=term1*term1+term2*term2+term3*term3+term4+term5+term6;
  double err=0.0;
  if(term>0.0) err=sqrt(term);
  return err;

}


//__________________ returns fake rate
double TGammaJets::FakeRateFunc(double eleEt) {

double fkrt=0.0;
  if(eleEt>=13.0)
    {
      double arg1=fakerate_par[0]+fakerate_par[1]*eleEt;
      double lin1=fakerate_par[2];
      double scale1=fakerate_par[3];
      fkrt=scale1*(exp(arg1)+lin1);
    }  
  return fkrt;

}

//________________________________________ returns fake rates
double TGammaJets::FakeRateFuncErr(double eleEt)
{

  double err=0.0;
  double _expo=0.0;
  if(fakerate_par[3]>0.0) _expo=FakeRateFunc(eleEt)/fakerate_par[3];
  double mult1=_expo-fakerate_par[2];
  double term1=_expo*_expo*fakerate_parerr[3]*fakerate_parerr[3];
  double sub_term1=mult1*mult1*(fakerate_parerr[0]*fakerate_parerr[0]
				+2.0*fakerate_corr_coeff[1][0]*fakerate_parerr[0]*fakerate_parerr[1]*eleEt
				+fakerate_parerr[1]*fakerate_parerr[1]*eleEt*eleEt);
  double sub_term2=2.0*fakerate_parerr[2]*mult1*(fakerate_corr_coeff[2][0]*fakerate_parerr[0]
						 +fakerate_corr_coeff[2][1]*fakerate_parerr[1]*eleEt);
  double sub_term3=fakerate_parerr[2]*fakerate_parerr[2];
  double term2=fakerate_par[3]*fakerate_par[3]*(sub_term1+sub_term2+sub_term3);
  if((term1+term2)>0.0) err=sqrt(term1+term2);
  return err;

}




/*############## END FAKE RATE FUNCTION STUFF ########################*/

void TGammaJets::FindConversion()
{
	//1. match the leading photon to the electron
	//2. see if that is a conversion
	//3. if so throw that event away. do not count it as W
	


}
bool TGammaJets::IsConversionElectron(TStnElectron* ele, 
				     double& minsep, double& mindt, double& radius) 
{

	const double dlamCut = 0.04;
	const double sepCut  = 0.2;
	const double ptCut   = 0.0;
	const bool reqCot    = false;

	int nconv = 0, ntrident = 0;
	bool hasCot = false;

	minsep =  999;
	mindt  =  999;
	radius = -999;

	int ntrk = fTrackBlock->NTracks();

	// this should not be happening
	int iEleTrk = ele->TrackNumber();
	if(iEleTrk<0 || iEleTrk>=ntrk) return false;
	TStnTrack* eleTrk = fTrackBlock->Track(iEleTrk);

	double  sep=-999, dt=-999;
	float conv_sign=-999;
	bool cand = false;

	for(int j=0; j<ntrk; j++) {
		TStnTrack* oTrk = fTrackBlock->Track(j);
		
		// only ask for a track from the electron track
		if ( j==iEleTrk ) continue;

		// require number of COT axial and stereo segments > 1
		hasCot = ((oTrk->NAxSeg()> 1) && (oTrk->NStSeg() > 1));
		if (!hasCot && reqCot)continue;
		 
		dt        = eleTrk->Lam0() - oTrk->Lam0();
		conv_sign = eleTrk->Charge() + oTrk->Charge();

    	bool eleCand = false;
    	
		for (int i=0; i < fElectronBlock->NElectrons(); i++) {
      	if (j == fElectronBlock->Electron(i)->TrackNumber()) eleCand = true;
    	}

  
		// opposite charge
		if (conv_sign != 0)continue;

		// delta cottheta cut
		if (fabs(dt) > dlamCut)continue;

		double rconv;
		sep = ConversionSep(eleTrk,oTrk,rconv);

		//separation cuts
		if (fabs(sep)> sepCut) continue;
    	
		nconv++;
    
		if (fabs(sep) <= minsep && fabs(dt) <= mindt) {
				minsep = fabs(sep);
				mindt  = fabs(dt);
				radius = rconv;
		}

		 // now check for trident
		 // by looking for a partner of second conversion 
		for (int k=0; k<ntrk; k++) {
			if ( (k != iEleTrk) && (k != j) ) {
				TStnTrack* kTrk = fTrackBlock->Track(k);
				dt        = oTrk->Lam0() - kTrk->Lam0();
				conv_sign = oTrk->Charge() + kTrk->Charge();
				sep       = ConversionSep(oTrk,kTrk,rconv);
				hasCot 	 = ( (kTrk->NAxSeg() > 1) && (kTrk->NStSeg() > 1) );

				if ( (kTrk->Pt() > ptCut) && (hasCot || (!reqCot)) &&
			   	conv_sign == 0 && fabs(dt) < dlamCut && fabs(sep) < sepCut) {
			  		ntrident++;
				}
			}
		} // end trident loop
	}

  
  if (ntrident > 0) cand = false;
  else if (nconv > 0) cand = true;
  else cand = false;
  return cand;
}


//____________________________________________________________________________
double TGammaJets::ConversionSep(TStnTrack* t0, TStnTrack* t1, double & rconv) {

  double  r1, r2, x1, x2, y1, y2, ff, dd, sep;

  r1 = 1.0/fabs(2*t0->C0());
  ff = t0->Charge()*r1+t0->D0();
  x1 = -ff*sin(t0->Phi0());
  y1 =  ff*cos(t0->Phi0());

  r2 = 1.0/fabs(2*t1->C0());
  ff = t1->Charge()*r2+t1->D0();
  x2 = -ff*sin(t1->Phi0());
  y2 =  ff*cos(t1->Phi0());

  dd  = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  
  if(t0->Charge()!=t1->Charge()) {
    sep = dd-r1-r2;
  } else {
    if (r1<r2) sep = r2 - dd - r1;
    else       sep = r1 - dd - r2;
  }

  double dx = x2 - x1;
  double dy = y2 - y1;
  double dr = sqrt( dx*dx + dy*dy );
    
  double xconv = x1 + dx/dr * (r1 + sep/2);
  double yconv = y1 + dy/dr * (r1 + sep/2);
  rconv = sqrt( xconv*xconv + yconv*yconv );	

  return sep;
}


/*############## END W/Z BACKGROUND ###################################*/


/*-------------- CALCULATE INVARIANT MASS ----------*/
double TGammaJets::GetInvariantMass(TLorentzVector* v1, TLorentzVector*v2,
												TLorentzVector* v3)
{
	TLorentzVector sum = *v1+*v2+*v3;
	return sum.M();
	
}

//_____________________________________________________________________________
double TGammaJets::GetInvariantMass(TLorentzVector* v1, TLorentzVector*v2)
{
	TLorentzVector sum = *v1 + *v2;
	return sum.M();
}


//_____________________________________________________________________________
//--- return the Et sorted photons. filter out the photons from the super photon vector.
//_____________________________________________________________________________
std::vector<TStnPhoton*> TGammaJets::GetSortedPhotonList()
{
	std::vector<TStnPhoton*> plist;

	for (int i =0 ; i< superpho_vec.size() ;++i)
		plist.push_back(superpho_vec[i].pho);
	
	return plist;
}
//_____________________________________________________________________________
//--- return the jets. filter out the jets from the super jet vector.
//_____________________________________________________________________________
std::vector<TStnJet*> TGammaJets::GetSortedJetList()
{
	std::vector<TStnJet*> jlist;

	for (int i =0 ; i< superjet_vec.size() ;++i)
		jlist.push_back(superjet_vec[i].jet);
	
	return jlist;
}
/*===================================================================*/
//	Fill data Block 
/*===================================================================*/
void TGammaJets::FillDataBlocks(int ientry)
{

	fPhotonBlock->GetEntry(ientry);
	fJetBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);
	fGenpBlock->GetEntry(ientry);
	fL3SummaryBlock->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
	fVertexBlock->GetEntry(ientry);
	fTrackBlock->GetEntry(ientry);
}

/*===================================================================*/
//	ALL BOOKINGS FOR HISTOGRAMS SHOULD GO DOWN HERE
/*===================================================================*/
void TGammaJets::BookInvariantMassHistograms()
{
	std::cout << "BOOKING MASS PLOTS...";

	float mp_min = 0, mp_max = 800, mp_bin_size = 0.2;
	int mp_bins = (int) ((mp_max - mp_min) / mp_bin_size);

	MassPlot.photon_j1    = new TH1F("photon_j1_mass","Invariant Mass of the Photon and the Leading Jet",mp_bins,mp_min,mp_max);
	MassPlot.photon_j2    = new TH1F("photon_j2_mass","Invariant Mass of the Photon and the 2nd Leading Jet",mp_bins,mp_min,mp_max);
	MassPlot.photon_j1_j2 = new TH1F("photon_j1_j2_mass","Invariant Mass of the Photon + 2 Leading Jets",mp_bins,mp_min,mp_max);
	MassPlot.j1_j2        = new TH1F("j1_j2_mass","Invariant Mass of the 2 leading Jets",mp_bins,mp_min,mp_max);
	MassPlot.photon_j1_j2_fake = new TH1F("photon_j1_j2_mass_fake","Invariant Mass of the Electron + 2 Leading Jets",mp_bins,mp_min,mp_max);

	TFolder* new_folder = GetHistoFolder("Inv_Mass_Plots","Invriant Mass Plots");
	new_folder->Add(MassPlot.photon_j1);
	new_folder->Add(MassPlot.photon_j2);
	new_folder->Add(MassPlot.photon_j1_j2);
	new_folder->Add(MassPlot.j1_j2);
	new_folder->Add(MassPlot.photon_j1_j2_fake);

	std::cout << "DONE.\n";

}




/*===================================================================*/
//	Create a folder in"Hist" 
/*===================================================================*/
TFolder* TGammaJets::GetHistoFolder(char *name, char* title)
{
	char folder_name[200];
	char folder_title[200];
	TFolder* hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
	sprintf(folder_name,name);
	sprintf(folder_title,title);
	TFolder* new_folder = (TFolder*) hist_folder->FindObject(folder_name);
	if (! new_folder) new_folder = hist_folder->AddFolder(folder_name,folder_title,NULL);
	return new_folder;
}

//_____________________________________________________________________________
void TGammaJets::BookGeneralHistograms()
{
	std::cout << "Booking General Plots...";
	
	GeneralPlot.photon_loose_and_tight_cuts_cumm = new TH1F("loose_tight_cuts_cummulative","Cummulative Photon Cuts (Loose and Tight)",22,0,22);
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(1,"Trigger");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(2,"L:Central");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(3,"L:Etc>30");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(4,"L:XCes");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(5,"L:ZCes");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(6,"L:HadEm");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(7,"L:TrkPt");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(8,"L:TrkIso");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(9,"L:IsoEtCorr");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(10,"T:EtCorr");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(11,"T:XCes");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(12,"T:ZCes");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(13,"T:HadEm");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(14,"T:Chi2Mean");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(15,"T:N3d");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(16,"T:TrkPt");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(17,"T:TrkIso");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(18,"T:IsoEt4Corr");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(19,"T:2ndCesWire");
	GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(20,"T:2ndCesStrip");
	//GeneralPlot.photon_loose_and_tight_cuts_cumm->GetXaxis()->SetBinLabel(21,"2 Jets>15GeV");
	
	
	GeneralPlot.pho_IDbins_pass = new TH1F("pho_IDbins_pass","The Leading Photon Pass Binning",8,0,8);
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(1,"Fail all");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(2,"Pass Loose Pho");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(3,"Pass Tight Pho");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(4,"Pass Pho-Like Ele");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(5,"Pass Loose && Tight");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(6,"Pass Loose && Pho-Like");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(7,"Pass Tight && Pho-Like");
	GeneralPlot.pho_IDbins_pass->GetXaxis()->SetBinLabel(8,"Pass all");

	GeneralPlot.pho_IDbins_fail = new TH1F("pho_IDbins_fail","The Leading Photon Fail Binning",8,0,8);
	GeneralPlot.pho_IDbins_fail->GetXaxis()->SetBinLabel(1,"Fail all");
	GeneralPlot.pho_IDbins_fail->GetXaxis()->SetBinLabel(2,"Fail Loose Pho");
	GeneralPlot.pho_IDbins_fail->GetXaxis()->SetBinLabel(3,"Fail Tight Pho");
	GeneralPlot.pho_IDbins_fail->GetXaxis()->SetBinLabel(4,"Fail Pho-Like Ele");
	GeneralPlot.pho_IDbins_fail->GetXaxis()->SetBinLabel(5,"Fail Loose && Tight");
	GeneralPlot.pho_IDbins_fail->GetXaxis()->SetBinLabel(6,"Fail Loose && Pho-Like");
	GeneralPlot.pho_IDbins_fail->GetXaxis()->SetBinLabel(7,"Fail Tight && Pho-Like");
	GeneralPlot.pho_IDbins_fail->GetXaxis()->SetBinLabel(8,"Pass all");

	GeneralPlot.Nvertices = new TH1F("Nvertices","Number of Class12 vertices per event",10,0,10);
	GeneralPlot.vertexZ = new TH1F("BestVertex_Z","Best Class12 vertex Z position(events with vertices>=1)",100,-100,100);

	TFolder* new_folder = GetHistoFolder("General_Plots","General Plots");
	new_folder->Add(GeneralPlot.photon_loose_and_tight_cuts_cumm);	
	new_folder->Add(GeneralPlot.pho_IDbins_pass);
	new_folder->Add(GeneralPlot.pho_IDbins_fail);
	new_folder->Add(GeneralPlot.Nvertices);
	new_folder->Add(GeneralPlot.vertexZ);
	std::cout << "DONE\n";

} //BookGeneralHistograms


//_____________________________________________________________________________
void TGammaJets::BookMetHistograms()
{
	std::cout << "BOOKING Met Plots...";
	
	MetPlot.Met_0   	= new TH1F("Met_0","Met(0) of events pass Iso25 trigger",100,0,100);
	MetPlot.MetPhi_0  = new TH1F("MetPhi_0","Met(0) Phi of events pass Iso25 trigger",70,0,7);
	MetPlot.MetX_0   	= new TH1F("MetX_0","MetX(0) of events pass Iso25 trigger",50,0,50);
	MetPlot.MetY_0   	= new TH1F("MetY_0","MetY(0) of events pass Iso25 trigger",50,0,50);
	MetPlot.Sumet_0   = new TH1F("Sumet_0","Sumet_0 of events pass Iso25 trigger",100,0,500);
	MetPlot.Sumet_1   = new TH1F("Sumet_1","Sumet_1:Htc of events pass Iso25 trigger",100,0,500);
	MetPlot.Sumet_2   = new TH1F("Sumet_2","Sumet_2:Sumetjet of events pass Iso25 trigger",100,0,500);
	MetPlot.Metsig   	= new TH1F("Metsig","MetSig of events pass Iso25 trigger",50,0,50);
	MetPlot.Z0   		= new TH1F("Z0","Z0 position of events pass Iso25 trigger",50,0,50);

	TFolder* new_folder = GetHistoFolder("Met_Plots"," Missing Et Plots");
	new_folder->Add(MetPlot.Met_0);	
	new_folder->Add(MetPlot.MetPhi_0);	
	new_folder->Add(MetPlot.MetX_0);	
	new_folder->Add(MetPlot.MetY_0);	
	new_folder->Add(MetPlot.Sumet_0);	
	new_folder->Add(MetPlot.Sumet_1);	
	new_folder->Add(MetPlot.Sumet_2);	
	new_folder->Add(MetPlot.Metsig);	
	new_folder->Add(MetPlot.Z0);	
	std::cout << "DONE.\n";
}
//_____________________________________________________________________________
void TGammaJets::BookBackgroundHistograms()
{
	std::cout << "BOOKING Background Plots...";

	BackgroundPlot.Wele_EtCorr   	= new TH1F("Wele_EtCorr","W Electron Et",100,0,100);
	BackgroundPlot.Wele_HadEm		= new TH1F("Wele_HadEm","W Electron HadEm",50,0,1);
	BackgroundPlot.Wele_IsoEtCorr	= new TH1F("Wele_IsoEtCorr","W Electron IsoEtCorr",140,-4,10);
	BackgroundPlot.Wele_TrkPt 		= new TH1F("Wele_Trkpt","W Electron TrkPt",1000,0,100);
	BackgroundPlot.W_Met 			= new TH1F("W_Met","Met for W events",100,0,100);
	BackgroundPlot.W_eleEt_to_Met_ratio		= new TH1F("W_eleEt2Met","W Eelectron Et to Met Ratio",40,0,2);
	BackgroundPlot.W_Sumet_0 		= new TH1F("W_Sumet_0","W Events Sumet(0)",100,0,500);
	BackgroundPlot.W_Sumet_1_Htc 	= new TH1F("W_Sumet_1_Htc","W Events Sumet(1):Htc",100,0,500);
	BackgroundPlot.W_Sumet_2_Sumetjet	= new TH1F("W_Sumet_2_Sumetjet","W Events Sumet(2):Sumetjet",100,0,500);
	BackgroundPlot.W_Mass 			= new TH1F("W_mass","W Transverse Invriant Mass",100,20,120);
	BackgroundPlot.W_Pt 				= new TH1F("W_Pt","W Pt",100,0,100);

	BackgroundPlot.Zele1_EtCorr   	= new TH1F("Zele1_EtCorr","Z Leading Electron Et",100,0,100);
	BackgroundPlot.Zele2_EtCorr   	= new TH1F("Zele2_EtCorr","Z 2nd Leading Electron Et",100,0,100);
	BackgroundPlot.Zele1_HadEm			= new TH1F("Zele1_HadEm","Z Leading Electron HadEm",50,0,1);
	BackgroundPlot.Zele1_IsoEtCorr	= new TH1F("Zele1_IsoEtCorr","Z Leading Electron IsoEtCorr",140,-4,10);
	BackgroundPlot.Zele1_TrkPt 		= new TH1F("Zele1_Trkpt","Z Leading Electron TrkPt",1000,0,100);
	BackgroundPlot.Zele2_TrkPt 		= new TH1F("Zele2_Trkpt","Z 2nd Leading Electron TrkPt",1000,0,100);
	BackgroundPlot.Z_Met 				= new TH1F("Z_Met","Met for Z events",100,0,50);
	BackgroundPlot.Z_Sumet_0 			= new TH1F("Z_Sumet_0","Z Events Sumet(0)",100,0,500);
	BackgroundPlot.Z_Sumet_1_Htc 		= new TH1F("Z_Sumet_1_Htc","Z Events Sumet(1):Htc",100,0,500);
	BackgroundPlot.Z_Sumet_2_Sumetjet= new TH1F("Z_Sumet_2_Sumetjet","Z Events Sumet(2):Sumetjet",100,0,500);
	BackgroundPlot.Z_Mass 				= new TH1F("Z_mass","Z Invriant Mass",100,20,120);
	BackgroundPlot.Z_Pt 					= new TH1F("Z_Pt","Z Pt",100,0,100);
	BackgroundPlot.Z_eleDelPhi 		= new TH1F("Z_eleDelPhi","Z: electrons Phi seperation",320,0,3.2);
	BackgroundPlot.Z_eleDelEta 		= new TH1F("Z_eleDelEta","Z: electrons Eta seperation",40,0,4); // cal. resolution,order of 0.1
	BackgroundPlot.Z_eleEtratio		= new TH1F("Z_eleEtratio","Z: electrons Et ratio (leading to 2nd leading) ",100,0,2);

	TFolder* new_folder = GetHistoFolder("Background_Plots","Background Plots");
	new_folder->Add(BackgroundPlot.Wele_EtCorr);	
	new_folder->Add(BackgroundPlot.Wele_HadEm);	
	new_folder->Add(BackgroundPlot.Wele_IsoEtCorr);	
	new_folder->Add(BackgroundPlot.Wele_TrkPt);	
	new_folder->Add(BackgroundPlot.W_Met);	
	new_folder->Add(BackgroundPlot.W_eleEt_to_Met_ratio);	
	new_folder->Add(BackgroundPlot.W_Sumet_0);	
	new_folder->Add(BackgroundPlot.W_Sumet_1_Htc);	
	new_folder->Add(BackgroundPlot.W_Sumet_2_Sumetjet);	
	new_folder->Add(BackgroundPlot.W_Mass);
	new_folder->Add(BackgroundPlot.W_Pt);

	new_folder->Add(BackgroundPlot.Zele1_EtCorr);	
	new_folder->Add(BackgroundPlot.Zele2_EtCorr);	
	new_folder->Add(BackgroundPlot.Zele1_HadEm);	
	new_folder->Add(BackgroundPlot.Zele1_IsoEtCorr);	
	new_folder->Add(BackgroundPlot.Zele1_TrkPt);	
	new_folder->Add(BackgroundPlot.Zele2_TrkPt);	
	new_folder->Add(BackgroundPlot.Z_Met);	
	new_folder->Add(BackgroundPlot.Z_Sumet_0);	
	new_folder->Add(BackgroundPlot.Z_Sumet_1_Htc);	
	new_folder->Add(BackgroundPlot.Z_Sumet_2_Sumetjet);	
	new_folder->Add(BackgroundPlot.Z_Mass);	
	new_folder->Add(BackgroundPlot.Z_Pt);
	new_folder->Add(BackgroundPlot.Z_eleDelPhi);	
	new_folder->Add(BackgroundPlot.Z_eleDelEta);	
	new_folder->Add(BackgroundPlot.Z_eleEtratio);	
	
	std::cout << "DONE.\n";
}



//_____________________________________________________________________________
void TGammaJets::BookL3SummaryHistograms()
{
	std::cout << "BOOKING L3Summary Plots...";
/*	char folder_name[200];
	TFolder* hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
	sprintf(folder_name,"L3Summary_plots");
	TFolder* new_folder = (TFolder*) hist_folder->FindObject(folder_name);
	if (! new_folder) new_folder = hist_folder->AddFolder(folder_name,"L3Summary Plots");
*/
	L3SummaryPlot.NEm   = new TH1F("NEm","L3Summary: Number of Em Objects/event",15,0,15);
	L3SummaryPlot.NMuon = new TH1F("NMuon","L3Summary: Number of Muons/event",15,0,15);
	L3SummaryPlot.NTau  = new TH1F("NTau","L3Summary: Number of Taus/event",15,0,15);
	L3SummaryPlot.NJet4 = new TH1F("NJet4","L3Summary: Number of Jets(0.4)/event",15,0,15);
	L3SummaryPlot.NJet7 = new TH1F("NJet7","L3Summary: Number of Jets(0.7)/event",15,0,15);
	L3SummaryPlot.NMet = new TH1F("NMet","L3Summary: Number of Met Objects/event",15,0,15);
	L3SummaryPlot.NCotTrack = new TH1F("NCotTrack","L3Summary: Number of COT Tracks/event",50,0,50);
	L3SummaryPlot.NSvtTrack = new TH1F("NSvtTrack","L3Summary: Number of Silicon Tracks/event",50,0,50);
	L3SummaryPlot.NTau_to_NEm = new TH1F("NTau_to_NEm","L3Summary: Taus/Em  per event",50,0,5);
	L3SummaryPlot.NMuon_to_NEm = new TH1F("NMuon_to_NEm","L3Summary: Muons/Em  per event",50,0,5);
	L3SummaryPlot.NJet4_to_NJet7 = new TH1F("NJet4_to_NJet7","L3Summary: Jets(0.4)/Jets(0.7) ",50,0,5);
	L3SummaryPlot.NJet7_to_NJet4 = new TH1F("NJet7_to_NJet4","L3Summary: Jets(0.7)/Jets(0.4) ",50,0,5);
	L3SummaryPlot.NSvtTrack_to_NCotTrack = new TH1F("NSvtTrack_to_NCotTrack","L3Summary: SvtTracks/CotTracks",50,0,5);
	L3SummaryPlot.NCotTrack_to_NSvtTrack = new TH1F("NCotTrack_to_NSvtTrack","L3Summary: CotTracks/SvtTracks",50,0,5);
	L3SummaryPlot.NCotTrack_vs_NEm = new TH2F("NCotTrack_vs_NEm","L3Summary: Number of COT Tracks vs Em Objects",20,0,20,50,0,50);
	L3SummaryPlot.NSvtTrack_vs_NEm = new TH2F("NSvtTrack_vs_NEm","L3Summary: Number of Silicon Tracks vs Em Objects",20,0,20,50,0,50);
	
	TFolder* new_folder = GetHistoFolder("L3Summary_Plots","L3 Summary Plots");
	new_folder->Add(L3SummaryPlot.NEm);	
	new_folder->Add(L3SummaryPlot.NMuon);	
	new_folder->Add(L3SummaryPlot.NTau);	
	new_folder->Add(L3SummaryPlot.NJet4);	
	new_folder->Add(L3SummaryPlot.NJet7);	
	new_folder->Add(L3SummaryPlot.NMet);	
	new_folder->Add(L3SummaryPlot.NCotTrack);	
	new_folder->Add(L3SummaryPlot.NSvtTrack);	
	new_folder->Add(L3SummaryPlot.NTau_to_NEm);	
	new_folder->Add(L3SummaryPlot.NMuon_to_NEm);	
	new_folder->Add(L3SummaryPlot.NJet4_to_NJet7);	
	new_folder->Add(L3SummaryPlot.NJet7_to_NJet4);	
	new_folder->Add(L3SummaryPlot.NSvtTrack_to_NCotTrack);	
	new_folder->Add(L3SummaryPlot.NCotTrack_to_NSvtTrack);	
	new_folder->Add(L3SummaryPlot.NCotTrack_vs_NEm);	
	new_folder->Add(L3SummaryPlot.NSvtTrack_vs_NEm);

	std::cout << "DONE\n";

} //BookL3SummaryHistograms

//_____________________________________________________________________________
void TGammaJets::BookHEPGPhotonMatchingHistograms()
{
	std::cout << "Booking HEPG Photon Matching Plots...";
	HEPGPhoMatchPlot.DelEta	= new TH1F("DelEta","Delta Eta of all HEPG Matches",80,-4,4);
	HEPGPhoMatchPlot.DelEta->SetYTitle("Events/0.1");
	HEPGPhoMatchPlot.DelPhi	= new TH1F("DelPhi","Delta Phi of all HEPG Matches",175,0,3.5);
	HEPGPhoMatchPlot.DelPhi->SetYTitle("Events/0.02");
	HEPGPhoMatchPlot.DelR	= new TH1F("DelR","Delta R of all HEPG Matches",100,0,1.0);
	HEPGPhoMatchPlot.DelR->SetYTitle("Events/0.002");
	HEPGPhoMatchPlot.Nmatches	= new TH1F("Nmatches","Number of matches of per event",10,0,10);
	HEPGPhoMatchPlot.matched_Object	= new TH1F("MacthedObject","Matched photon's Mother",8,0,8);
	HEPGPhoMatchPlot.matched_Object->GetXaxis()->SetBinLabel(1,"Pi0");
	HEPGPhoMatchPlot.matched_Object->GetXaxis()->SetBinLabel(2,"Electron");
	HEPGPhoMatchPlot.matched_Object->GetXaxis()->SetBinLabel(3,"Muon");
	HEPGPhoMatchPlot.matched_Object->GetXaxis()->SetBinLabel(4,"Tau");
	HEPGPhoMatchPlot.matched_Object->GetXaxis()->SetBinLabel(5,"Pi+/-");
	HEPGPhoMatchPlot.matched_Object->GetXaxis()->SetBinLabel(6,"Eta");
	HEPGPhoMatchPlot.matched_Object->GetXaxis()->SetBinLabel(7,"Omega");
	HEPGPhoMatchPlot.matched_Object->GetXaxis()->SetBinLabel(8,"Other");
	HEPGPhoMatchPlot.HEPGObj_Et	= new TH1F("HepgObjEt","All matching HEPG Objects Et",50,0,100);
	HEPGPhoMatchPlot.HEPGObj_Et->SetYTitle("Events/2.0");
	HEPGPhoMatchPlot.Et_ratio		= new TH1F("OffPhoEt_to_HepgPhoEt","Offline leading Photon Et to Matched HEPG Object Et",40,0,2);
	HEPGPhoMatchPlot.Et_ratio->SetYTitle("Events/0.05");
	HEPGPhoMatchPlot.HEPGObj_Eta	= new TH1F("HepgPhoEta","All matching HEPG Objects Eta",80,-4,4);
	HEPGPhoMatchPlot.HEPGObj_Eta->SetYTitle("Events/0.1");

	TFolder* new_folder = GetHistoFolder("HEPGMatching_Plots","Photon HEPG Matching Plots");
	new_folder->Add(HEPGPhoMatchPlot.DelEta);
	new_folder->Add(HEPGPhoMatchPlot.DelPhi);
	new_folder->Add(HEPGPhoMatchPlot.DelR);
	new_folder->Add(HEPGPhoMatchPlot.Nmatches);
	new_folder->Add(HEPGPhoMatchPlot.matched_Object);
	new_folder->Add(HEPGPhoMatchPlot.HEPGObj_Et);
	new_folder->Add(HEPGPhoMatchPlot.Et_ratio);
	new_folder->Add(HEPGPhoMatchPlot.HEPGObj_Eta);

	std::cout << "DONE.\n";
}

//_____________________________________________________________________________
void TGammaJets::BookEmObjectRemovedFromJetsHistograms()
{
	std::cout << "Booking EmObjects Removed from Jet Block Plots...";
	EmJetMatchPlot.DelEta	= new TH1F("DelEta","Delta Eta of EMObject and Matched Jet",80,0,4);
	EmJetMatchPlot.DelEta->SetYTitle("Events/0.05");
	EmJetMatchPlot.DelPhi	= new TH1F("DelPhi","Delta Phi of EMObject and Matched Jet",175,0,3.5);
	EmJetMatchPlot.DelPhi->SetYTitle("Events/0.02");
	EmJetMatchPlot.DelR	= new TH1F("DelR","Delta R of EMObject and Matched Jet",100,0,1.0);
	EmJetMatchPlot.DelR->SetYTitle("Events/0.002");
	EmJetMatchPlot.NOmatches	= new TH1F("Nmatches","Number of Em objects did not match /event",10,0,10);
	EmJetMatchPlot.Et_ratio		= new TH1F("EmJetEtRatio","EmObject Et to Matched Jet Et",40,0,2);
	EmJetMatchPlot.Et_ratio->SetYTitle("Events/0.05");

	TFolder* new_folder = GetHistoFolder("EM2JetMatching_Plots","Photons removed from Jet Block Plots");
	new_folder->Add(EmJetMatchPlot.DelEta);
	new_folder->Add(EmJetMatchPlot.DelPhi);
	new_folder->Add(EmJetMatchPlot.DelR);
	new_folder->Add(EmJetMatchPlot.NOmatches);
	new_folder->Add(EmJetMatchPlot.Et_ratio);

	std::cout << "DONE.\n";
}

//_____________________________________________________________________________
void TGammaJets::BookJetHistograms()
{
	std::cout << "Booking Leading Jets Plots...";

	Jet1Plot.Etraw	= new TH1F("j1Etraw","Jet1 raw Et",50,0,100);
	Jet1Plot.Etraw->SetYTitle("Events/2.0");
	Jet1Plot.DetEta= new TH1F("j1Eta","Eta of Jet1 ",80,0,4);
	Jet1Plot.DetEta->SetYTitle("Events/0.05");
	Jet1Plot.Phi	= new TH1F("j1Phi","Phi of Jet1",175,0,3.5);
	Jet1Plot.Phi->SetYTitle("Events/0.02");
	Jet1Plot.Emfr	= new TH1F("j1Emfr","Emfr of Jet1",20,0,2);
	Jet1Plot.Emfr->SetYTitle("Events/0.1");
	Jet1Plot.EOverP= new TH1F("j1EP","E/P of Jet1",20,0,2);
	Jet1Plot.EOverP->SetYTitle("Events/0.1");

	Jet2Plot.Etraw	= new TH1F("j2Etraw","Jet2 raw Et",50,0,100);
	Jet2Plot.Etraw->SetYTitle("Events/2.0");
	Jet2Plot.DetEta= new TH1F("j2Eta","Eta of Jet2 ",80,0,4);
	Jet2Plot.DetEta->SetYTitle("Events/0.05");
	Jet2Plot.Phi	= new TH1F("j2Phi","Phi of Jet2",175,0,3.5);
	Jet2Plot.Phi->SetYTitle("Events/0.02");
	Jet2Plot.Emfr	= new TH1F("j2Emfr","Emfr of Jet2",20,0,2);
	Jet2Plot.Emfr->SetYTitle("Events/0.1");
	Jet2Plot.EOverP= new TH1F("j2EP","E/P of Jet2",20,0,2);


	JetPlot.DelEta	= new TH1F("jetsDelEta","Delta Eta of two Jets",80,0,4);
	JetPlot.DelEta->SetYTitle("Events/0.05");
	JetPlot.DelPhi	= new TH1F("jetsDelPhi","Delta Phi of two Jet",175,0,3.5);
	JetPlot.DelPhi->SetYTitle("Events/0.02");
	JetPlot.DelR	= new TH1F("jetsDelR","Delta R of of two Jet",100,0,1.0);
	JetPlot.DelR->SetYTitle("Events/0.002");
	JetPlot.EtRatio= new TH1F("jetsEtRatio","Et ration of two Jets",40,0,2);
	JetPlot.EtRatio->SetYTitle("Events/0.05");

	TFolder* new_folder = GetHistoFolder("TwoJets_Plots","Leading Two Jets Plots");

	new_folder->Add(Jet1Plot.Etraw);
	new_folder->Add(Jet1Plot.DetEta);
	new_folder->Add(Jet1Plot.Phi);
	new_folder->Add(Jet1Plot.Emfr);
	new_folder->Add(Jet1Plot.EOverP);

	new_folder->Add(Jet2Plot.Etraw);
	new_folder->Add(Jet2Plot.DetEta);
	new_folder->Add(Jet2Plot.Phi);
	new_folder->Add(Jet2Plot.Emfr);
	new_folder->Add(Jet2Plot.EOverP);
	
	new_folder->Add(JetPlot.DelEta);
	new_folder->Add(JetPlot.DelPhi);
	new_folder->Add(JetPlot.DelR);
	new_folder->Add(JetPlot.EtRatio);

	std::cout << "DONE.\n";
}
//_____________________________________________________________________________
void TGammaJets::BookPhotonHistograms()
{
	std::cout << "Booking Photon Plots...";
	PhotonPlot.Detector_b 	= new TH1F("Detector_b","Photon Detector Region",3,0,3);
	PhotonPlot.EtCorr_b 		= new TH1F("EtCorr_b","Photon Corrected Et",600,0,120);
	PhotonPlot.XCes_b 		= new TH1F("XCes_b","Photon XCes",640,-32,32);
	PhotonPlot.ZCes_b 		= new TH1F("ZCes_b","Photon ZCes",2000,-250,250);
	PhotonPlot.HadEm_b 		= new TH1F("HadEm_b","Photon HadEm",50,0,0.5);
	PhotonPlot.IsoEtCorr_b 	= new TH1F("IsoEtCorr_b","Photon IsoCorr",1500,-5,10);
	PhotonPlot.Chi2Mean_b 	= new TH1F("Chi2Mean_b","Photon Chi2Mean (Wire+Strip/2)",800,0,80);
	PhotonPlot.N3d_b 			= new TH1F("N3d_b","Photon N3d",10,0,10);
	PhotonPlot.TrkPt_b 		= new TH1F("Trkpt_b","Photon TrkPt",1000,0,100);
	PhotonPlot.TrkIso_b 		= new TH1F("TrkIso_b","Photon TrkIso",150,0,15);
	PhotonPlot.Ces2Wire_b 	= new TH1F("Ces2Wire_b","Photon CES(2nd) Wire",400,0,40);
	PhotonPlot.Ces2Strip_b 	= new TH1F("Ces2Strip_b","Photon CES(2nd) Strip",400,0,40);

	PhotonPlot.Detector_aL 	= new TH1F("Detector_aL","Photon Detector Region",3,0,3);
	PhotonPlot.EtCorr_aL		= new TH1F("EtCorr_aL","Photon Corrected Et",600,0,120);
	PhotonPlot.XCes_aL 		= new TH1F("XCes_aL","Photon XCes",640,-32,32);
	PhotonPlot.ZCes_aL 		= new TH1F("ZCes_aL","Photon ZCes",2000,-250,250);
	PhotonPlot.HadEm_aL 		= new TH1F("HadEm_aL","Photon HadEm",50,0,0.5);
	PhotonPlot.IsoEtCorr_aL	= new TH1F("IsoCorr_aL","Photon IsoCorr",1500,-5,10);
	PhotonPlot.Chi2Mean_aL 	= new TH1F("Chi2Mean_aL","Photon Chi2Mean (Wire+Strip/2)",800,0,80);
	PhotonPlot.N3d_aL 		= new TH1F("N3d_aL","Photon N3d",10,0,10);
	PhotonPlot.TrkPt_aL 		= new TH1F("Trkpt_aL","Photon TrkPt",1000,0,100);
	PhotonPlot.TrkIso_aL 	= new TH1F("TrkIso_aL","Photon TrkIso",150,0,15);
	PhotonPlot.Ces2Wire_aL 	= new TH1F("Ces2Wire_aL","Photon CES(2nd) Wire",400,0,40);
	PhotonPlot.Ces2Strip_aL	= new TH1F("Ces2Strip_aL","Photon CES(2nd) Strip",400,0,40);
	
	PhotonPlot.Detector_aT 	= new TH1F("Detector_aT","Photon Detector Region",3,0,3);
	PhotonPlot.EtCorr_aT		= new TH1F("EtCorr_aT","Photon Corrected Et",600,0,120);
	PhotonPlot.XCes_aT 		= new TH1F("XCes_aT","Photon XCes",640,-32,32);
	PhotonPlot.ZCes_aT 		= new TH1F("ZCes_aT","Photon ZCes",2000,-250,250);
	PhotonPlot.HadEm_aT 		= new TH1F("HadEm_aT","Photon HadEm",50,0,0.5);
	PhotonPlot.IsoEtCorr_aT	= new TH1F("IsoCorr_aT","Photon IsoCorr",1500,-5,10);
	PhotonPlot.Chi2Mean_aT 	= new TH1F("Chi2Mean_aT","Photon Chi2Mean (Wire+Strip/2)",800,0,80);
	PhotonPlot.N3d_aT 		= new TH1F("N3d_aT","Photon N3d",10,0,10);
	PhotonPlot.TrkPt_aT 		= new TH1F("Trkpt_aT","Photon TrkPt",1000,0,100);
	PhotonPlot.TrkIso_aT 	= new TH1F("TrkIso_aT","Photon TrkIso",150,0,15);
	PhotonPlot.Ces2Wire_aT 	= new TH1F("Ces2Wire_aT","Photon CES(2nd) Wire",400,0,40);
	PhotonPlot.Ces2Strip_aT	= new TH1F("Ces2Strip_aT","Photon CES(2nd) Strip",400,0,40);

	PhotonPlot.EtCorr_b->SetYTitle("Events/0.2");
	PhotonPlot.XCes_b->SetYTitle("Events/0.1");
	PhotonPlot.ZCes_b->SetYTitle("Events/0.25");
	PhotonPlot.HadEm_b->SetYTitle("Events/0.01");
	PhotonPlot.IsoEtCorr_b->SetYTitle("Events/0.01");
	PhotonPlot.Chi2Mean_b->SetYTitle("Events/0.1");
	PhotonPlot.N3d_b->SetYTitle("Events/1.0");
	PhotonPlot.TrkPt_b->SetYTitle("Events/0.1");
	PhotonPlot.TrkIso_b->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Wire_b->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Strip_b->SetYTitle("Events/0.1");

	PhotonPlot.EtCorr_aL->SetYTitle("Events/0.2");
	PhotonPlot.XCes_aL->SetYTitle("Events/0.1");
	PhotonPlot.ZCes_aL->SetYTitle("Events/0.25");
	PhotonPlot.HadEm_aL->SetYTitle("Events/0.01");
	PhotonPlot.IsoEtCorr_aL->SetYTitle("Events/0.01");
	PhotonPlot.Chi2Mean_aL->SetYTitle("Events/0.1");
	PhotonPlot.N3d_aL->SetYTitle("Events/1.0");
	PhotonPlot.TrkPt_aL->SetYTitle("Events/0.1");
	PhotonPlot.TrkIso_aL->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Wire_aL->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Strip_aL->SetYTitle("Events/0.1");

	PhotonPlot.EtCorr_aT->SetYTitle("Events/0.2");
	PhotonPlot.XCes_aT->SetYTitle("Events/0.1");
	PhotonPlot.ZCes_aT->SetYTitle("Events/0.25");
	PhotonPlot.HadEm_aT->SetYTitle("Events/0.01");
	PhotonPlot.IsoEtCorr_aT->SetYTitle("Events/0.01");
	PhotonPlot.Chi2Mean_aT->SetYTitle("Events/0.1");
	PhotonPlot.N3d_aT->SetYTitle("Events/1.0");
	PhotonPlot.TrkPt_aT->SetYTitle("Events/0.1");
	PhotonPlot.TrkIso_aT->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Wire_aT->SetYTitle("Events/0.1");
	PhotonPlot.Ces2Strip_aT->SetYTitle("Events/0.1");

	TFolder* new_folder = GetHistoFolder("Photon_Plots","Photon Plots");
	new_folder->Add(PhotonPlot.Detector_b);
	new_folder->Add(PhotonPlot.EtCorr_b);
	new_folder->Add(PhotonPlot.XCes_b);
	new_folder->Add(PhotonPlot.ZCes_b);
	new_folder->Add(PhotonPlot.HadEm_b);
	new_folder->Add(PhotonPlot.IsoEtCorr_b);
	new_folder->Add(PhotonPlot.Chi2Mean_b);
	new_folder->Add(PhotonPlot.N3d_b);
	new_folder->Add(PhotonPlot.TrkPt_b);
	new_folder->Add(PhotonPlot.TrkIso_b);
	new_folder->Add(PhotonPlot.Ces2Wire_b);
	new_folder->Add(PhotonPlot.Ces2Strip_b);

	new_folder->Add(PhotonPlot.Detector_aL);
	new_folder->Add(PhotonPlot.EtCorr_aL);
	new_folder->Add(PhotonPlot.XCes_aL);
	new_folder->Add(PhotonPlot.ZCes_aL);
	new_folder->Add(PhotonPlot.HadEm_aL);
	new_folder->Add(PhotonPlot.IsoEtCorr_aL);
	new_folder->Add(PhotonPlot.Chi2Mean_aL);
	new_folder->Add(PhotonPlot.N3d_aL);
	new_folder->Add(PhotonPlot.TrkPt_aL);
	new_folder->Add(PhotonPlot.TrkIso_aL);
	new_folder->Add(PhotonPlot.Ces2Wire_aL);
	new_folder->Add(PhotonPlot.Ces2Strip_aL);
	
	new_folder->Add(PhotonPlot.Detector_aT);
	new_folder->Add(PhotonPlot.EtCorr_aT);
	new_folder->Add(PhotonPlot.XCes_aT);
	new_folder->Add(PhotonPlot.ZCes_aT);
	new_folder->Add(PhotonPlot.HadEm_aT);
	new_folder->Add(PhotonPlot.IsoEtCorr_aT);
	new_folder->Add(PhotonPlot.Chi2Mean_aT);
	new_folder->Add(PhotonPlot.N3d_aT);
	new_folder->Add(PhotonPlot.TrkPt_aT);
	new_folder->Add(PhotonPlot.TrkIso_aT);
	new_folder->Add(PhotonPlot.Ces2Wire_aT);
	new_folder->Add(PhotonPlot.Ces2Strip_aT);
	std::cout << "DONE.\n";

} //BookPhotonHistograms

//_____________________________________________________________________________
void TGammaJets::BookPhotonIDcutsHistograms()
{
	std::cout << "Booking Photon \"N-1\" ID cuts Plots...";

	PhoIDcutPlot.Detector_b 	= new TH1F("Detector_b","Photon Detector Region :before any cuts",3,0,3);
	PhoIDcutPlot.EtCorr_b 		= new TH1F("EtCorr_b","Photon Corrected Et :before any cuts",600,0,120);
	PhoIDcutPlot.XCes_b 			= new TH1F("XCes_b","Photon XCes :before any cuts",640,-32,32);
	PhoIDcutPlot.ZCes_b 			= new TH1F("ZCes_b","Photon ZCes :before any cuts",2000,-250,250);
	PhoIDcutPlot.HadEm_b 		= new TH1F("HadEm_b","Photon HadEm :before any cuts",50,0,0.5);
	PhoIDcutPlot.IsoEtCorr_b 	= new TH1F("IsoEtCorr_b","Photon IsoEtCorr :before any cuts",1500,-5,10);
	PhoIDcutPlot.TrkPt_b 		= new TH1F("Trkpt_b","Photon TrkPt :before any cuts",1000,0,100);
	PhoIDcutPlot.TrkIso_b 		= new TH1F("TrkIso_b","Photon TrkIso :before any cuts",150,0,15);

	PhoIDcutPlot.Detector_a_n_1= new TH1F("Detector_a_n_1","Photon Detector Region :after N-1 cuts",3,0,3);
	PhoIDcutPlot.EtCorr_a_n_1	= new TH1F("EtCorr_a_n_1","Photon Corrected Et :after N-1 cuts",600,0,120);
	PhoIDcutPlot.XCes_a_n_1	 	= new TH1F("XCes_a_n_1","Photon XCes :after N-1 cuts",640,-32,32);
	PhoIDcutPlot.ZCes_a_n_1	 	= new TH1F("ZCes_a_n_1","Photon ZCes :after N-1 cuts",2000,-250,250);
	PhoIDcutPlot.HadEm_a_n_1 	= new TH1F("HadEm_a_n_1","Photon HadEm :after N-1 cuts",50,0,0.5);
	PhoIDcutPlot.IsoEtCorr_a_n_1	= new TH1F("IsoEtCorr_a_n_1","Photon IsoEtCorr :after N-1 cuts",1500,-5,10);
	PhoIDcutPlot.TrkPt_a_n_1 	= new TH1F("Trkpt_a_n_1","Photon TrkPt :after N-1 cuts",1000,0,100);
	PhoIDcutPlot.TrkIso_a_n_1	= new TH1F("TrkIso_a_n_1","Photon TrkIso :after N-1 cuts",150,0,15);

	PhoIDcutPlot.Detector_a 	= new TH1F("Detector_a","Photon Detector Region :after Central cut",3,0,3);
	PhoIDcutPlot.EtCorr_a		= new TH1F("EtCorr_a","Photon Corrected Et :after Et>30 cut",600,0,120);
	PhoIDcutPlot.XCes_a	 		= new TH1F("XCes_a","Photon XCes :after XCes cut",640,-32,32);
	PhoIDcutPlot.ZCes_a	 		= new TH1F("ZCes_a","Photon ZCes :after ZCes cut",2000,-250,250);
	PhoIDcutPlot.HadEm_a 		= new TH1F("HadEm_a","Photon HadEm :after HadEm cut",50,0,0.5);
	PhoIDcutPlot.IsoEtCorr_a	= new TH1F("IsoEtCorr_a","Photon IsoEtCorr :after IsoEtCorr cut",1500,-5,10);
	PhoIDcutPlot.TrkPt_a 		= new TH1F("Trkpt_a","Photon TrkPt :after TrkPt cut",1000,0,100);
	PhoIDcutPlot.TrkIso_a	 	= new TH1F("TrkIso_a","Photon TrkIso :after TrkIso cut",150,0,15);

	PhoIDcutPlot.EtCorr_b->SetYTitle("Events/0.2");
	PhoIDcutPlot.XCes_b->SetYTitle("Events/0.1");
	PhoIDcutPlot.ZCes_b->SetYTitle("Events/0.25");
	PhoIDcutPlot.HadEm_b->SetYTitle("Events/0.01");
	PhoIDcutPlot.IsoEtCorr_b->SetYTitle("Events/0.01");
	PhoIDcutPlot.TrkPt_b->SetYTitle("Events/0.1");
	PhoIDcutPlot.TrkIso_b->SetYTitle("Events/0.1");

	PhoIDcutPlot.EtCorr_a_n_1->SetYTitle("Events/0.2");
	PhoIDcutPlot.XCes_a_n_1->SetYTitle("Events/0.1");
	PhoIDcutPlot.ZCes_a_n_1->SetYTitle("Events/0.25");
	PhoIDcutPlot.HadEm_a_n_1->SetYTitle("Events/0.01");
	PhoIDcutPlot.IsoEtCorr_a_n_1->SetYTitle("Events/0.01");
	PhoIDcutPlot.TrkPt_a_n_1->SetYTitle("Events/0.1");
	PhoIDcutPlot.TrkIso_a_n_1->SetYTitle("Events/0.1");

	PhoIDcutPlot.EtCorr_a->SetYTitle("Events/0.2");
	PhoIDcutPlot.XCes_a->SetYTitle("Events/0.1");
	PhoIDcutPlot.ZCes_a->SetYTitle("Events/0.25");
	PhoIDcutPlot.HadEm_a->SetYTitle("Events/0.01");
	PhoIDcutPlot.IsoEtCorr_a->SetYTitle("Events/0.01");
	PhoIDcutPlot.TrkPt_a->SetYTitle("Events/0.1");
	PhoIDcutPlot.TrkIso_a->SetYTitle("Events/0.1");

	TFolder* new_folder = GetHistoFolder("Photon_IDcutPlots","Photon N-1 ID cut Plots");
	new_folder->Add(PhoIDcutPlot.Detector_b);
	new_folder->Add(PhoIDcutPlot.EtCorr_b);
	new_folder->Add(PhoIDcutPlot.XCes_b);
	new_folder->Add(PhoIDcutPlot.ZCes_b);
	new_folder->Add(PhoIDcutPlot.HadEm_b);
	new_folder->Add(PhoIDcutPlot.IsoEtCorr_b);
	new_folder->Add(PhoIDcutPlot.TrkPt_b);
	new_folder->Add(PhoIDcutPlot.TrkIso_b);

	new_folder->Add(PhoIDcutPlot.Detector_a_n_1);
	new_folder->Add(PhoIDcutPlot.EtCorr_a_n_1);
	new_folder->Add(PhoIDcutPlot.XCes_a_n_1);
	new_folder->Add(PhoIDcutPlot.ZCes_a_n_1);
	new_folder->Add(PhoIDcutPlot.HadEm_a_n_1);
	new_folder->Add(PhoIDcutPlot.IsoEtCorr_a_n_1);
	new_folder->Add(PhoIDcutPlot.TrkPt_a_n_1);
	new_folder->Add(PhoIDcutPlot.TrkIso_a_n_1);

	new_folder->Add(PhoIDcutPlot.Detector_a);
	new_folder->Add(PhoIDcutPlot.EtCorr_a);
	new_folder->Add(PhoIDcutPlot.XCes_a);
	new_folder->Add(PhoIDcutPlot.ZCes_a);
	new_folder->Add(PhoIDcutPlot.HadEm_a);
	new_folder->Add(PhoIDcutPlot.IsoEtCorr_a);
	new_folder->Add(PhoIDcutPlot.TrkPt_a);
	new_folder->Add(PhoIDcutPlot.TrkIso_a);
	std::cout << "DONE.\n";
}




//_____________________________________________________________________________
void TGammaJets::BookGenpStudyHistograms()
/*{{{*/
{
	std::cout << "BOOKING Genp Study Plots...";

	GenpPlot.Photon_Et       = new TH1F("PhotonEt","Genp. Highest Et Photon's Et",750,0,150);
	GenpPlot.Photon_Eta	    = new TH1F("PhotonEta","Genp. Highest Et Photon's Eta",400,-4,4);
	GenpPlot.Photon_Phi	    = new TH1F("PhotonEta","Genp. Highest Et Photon's Phi",350,0,3.5);
	GenpPlot.JetsLostInTheDetector_Eta	= new TH1F("JetsLost_Eta","Jets Lost in the Detector-Eta",400,-4,4);
	GenpPlot.JetsLostInTheDetector_Phi	= new TH1F("JetsLost_Phi","Jets Lost in the Detector-Phi",350,0,3.5);
	GenpPlot.PhotonsLostInTheDetector_Eta	= new TH1F("PhotonsLost_Eta","Photons Lost in the Detector-Eta",400,-4,4);
	GenpPlot.PhotonsLostInTheDetector_Phi	= new TH1F("PhotonsLost_Phi","Photons Lost in the Detector-Phi",350,0,3.5);

	TFolder* new_folder = GetHistoFolder("GenpStudy_Plots","Genp Study Plots");
	new_folder->Add(GenpPlot.Photon_Et);
	new_folder->Add(GenpPlot.Photon_Eta);
	new_folder->Add(GenpPlot.Photon_Phi);
	new_folder->Add(GenpPlot.JetsLostInTheDetector_Eta);
	new_folder->Add(GenpPlot.JetsLostInTheDetector_Phi);
	new_folder->Add(GenpPlot.PhotonsLostInTheDetector_Eta);
	new_folder->Add(GenpPlot.PhotonsLostInTheDetector_Phi);
		
	std::cout << "DONE.\n";

} //BookGenpStudyHistograms

/*}}}*/
//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int TGammaJets::EndJob() {

	printf("----- end job: ---- %s\n",GetName());

	printf("RUN SUMMARY\n");
	printf("TRIGGERS: Photon_ISO25,50,70\n");
	printf("Reject all events without a jet or a photon\n");
	printf("\n");


	std::cout << "                           MC ? = " << qMc ;
	if (qMc) std::cout << " YES\n" << std::endl; 
	else std::cout << " NO\n" << std::endl; 
	
	std::cout.precision(3);
	float loose_eff = 0, tight_eff=0, w_2_trig_perc=0,w_2_pho_perc=0,z_2_trig_perc=0,z_2_pho_perc=0;
	if (counter.evtsPassTrigger) {
		loose_eff = 100*(counter.evtsPassLooseCuts/(counter.evtsPassTrigger*1.));
		tight_eff = 100 *(counter.evtsPassTightCuts/(counter.evtsPassTrigger*1.));
		w_2_trig_perc = (counter.Wevents/(counter.evtsPassTrigger*1.0))*100;
		z_2_trig_perc = (counter.Zevents/(counter.evtsPassTrigger*1.0))*100;
		z_2_trig_perc = (counter.Zevents/(counter.evtsPassTrigger*1.0))*100;
	}
	if (counter.evtsPassTightCuts) {
		w_2_pho_perc = (counter.Wevents/(counter.evtsPassTightCuts*1.0))*100;
		z_2_pho_perc = (counter.Zevents/(counter.evtsPassTightCuts*1.0))*100;
	}

	std::cout << "Evts. without a Jet ---------------- : "<< counter.evtsWithoutAjet  << std::endl;
	std::cout << "Evts. without a Photon               : "<< counter.evtsWithoutAphoton  << std::endl;
	std::cout << "Evts. pass good run selection  ----- : "<< counter.goodrunevts << std::endl;
	std::cout << "Evts. Passed Trigger(25+50+70) ----- : "<< counter.evtsPassTrigger << "\t(25=" << counter.evtsPass25Trigger << ", 50=" << counter.evtsPass50Trigger << ", 70=" << counter.evtsPass70Trigger << ")" <<std::endl;
	std::cout << "Evts. Passed Loose Photon Cuts       : "<< counter.evtsPassLooseCuts;
		std::cout << "\t(" << loose_eff <<"\% trig)" << std::endl;
	std::cout << "Evts. Passed Tight Photon Cuts       : "<< counter.evtsPassTightCuts;
		std::cout << "\t("<< tight_eff << "\% trig)" << std::endl;
//	std::cout << "Total evts before stripping          : "<< counter.total_evts_b4_strip << std::endl;
//	std::cout << "Total evts after stripping --------- : "<< counter.total_evts_a4_strip << std::endl;
	if (qMc) {
	std::cout << "Photons matched to Genp              : "<< counter.phos_match2_genp << std::endl;
	std::cout << "Photons did not match to Genp        : "<< counter.phos_nomatch2_genp << std::endl;
	std::cout << "Photons matched to a HEPG Pi0 ------ : "<< counter.HEPGmomPi0 << std::endl;
	std::cout << "Photons matched to a HEPG Pi+/-      : "<< counter.HEPGmomChargedPi << std::endl;
	std::cout << "Photons matched to a HEPG Eta        : "<< counter.HEPGmomEta << std::endl;
	std::cout << "Photons matched to a HEPG Omega ---- : "<< counter.HEPGmomOmega << std::endl;
	std::cout << "Photons matched to a HEPG Other      : "<< counter.HEPGmomOther << std::endl;
	}
	std::cout << "W candidate events ----------------- : "<< counter.Wevents;
		std::cout << "\t(" << w_2_trig_perc << "\%trig, " << w_2_pho_perc << "\%tight pho)" << std::endl;
	std::cout << "Z candidate events ----------------- : "<< counter.Zevents;
		std::cout << "\t(" << z_2_trig_perc << "\%trig, " << z_2_pho_perc << "\%tight pho)" << std::endl;
	std::cout << "Photon(std)+2Jets(>15Gev) events     : "<< counter.pho_2j_events << std::endl;
	

	std::cout << "Sum of fakerate ---------------------: "<< fakeRateSum << std::endl;
	
	std::cout << ""<< std::endl;
	return 0;
}
