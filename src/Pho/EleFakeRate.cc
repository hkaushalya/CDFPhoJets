//_____________________________________________________________________________
// this will calculate the fake rate for an event with e+2jets before phoenix rejection of photons
// it can be aceesed via GetNoPhxFakeRate()
// 1. electron tagging module must exitst
// 2. conversion removal module must exist
// 3. if two electrons are found sum of the two fake rates will be retuned.



#include "samantha/Pho/EleFakeRate.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>

ClassImp(EleFakeRate)
//_____________________________________________________________________________
EleFakeRate::EleFakeRate(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  iGoodRunBit(0),			//get rid of this option. i don't need this
  iPhoenixBit(1),
  sErrs(""),
  bRunPermit(true)
{
	//constructor
	std::cout << "Hello I am EleFakeRate module" << std::endl;
}

//_____________________________________________________________________________
EleFakeRate::~EleFakeRate() {
	//destructor
}

//_____________________________________________________________________________
void EleFakeRate::SaveHistograms() {
}

//_____________________________________________________________________________
void EleFakeRate::BookHistograms()
{
	// book all histograms here
	DeleteHistograms();

  	//char name [200];
  	//char title[200];
	TFolder *new_folder, *sub_folder;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");
	vCounterHistLabels.push_back("Events Passed");
	vCounterHistLabels.push_back("1 Jet Events");
	vCounterHistLabels.push_back("2 Jets Events");
	HistoMan.GetCounterHistograms(new_folder,hCount,this->GetName(),vCounterHistLabels);

	new_folder = GetHistFolder(this, "1Jet","#gamma + >=1 Jet case");
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h1_Evt,"#gamma + >=1 Jets(15GeV)");
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h1_Pho,"#gamma + >=1 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h1_Jet,"#gamma + >=1 Jets(15GeV)");
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h1_PhoJet,"#gamma + >=1 Jets(15GeV)");



	new_folder = GetHistFolder(this, "2Jet","#gamma + >=2 Jets case");
	
	sub_folder = new_folder->AddFolder("Event","Event Histograms");
	HistoMan.GetEventHistograms(sub_folder,h2_Evt,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
	HistoMan.GetPhotonHistograms(sub_folder,h2_Pho,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet1,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
	HistoMan.GetJetHistograms(sub_folder,h2_Jet2,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet1,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
	HistoMan.GetPhoton1JetHistograms(sub_folder,h2_PhoJet2,"#gamma + >=2 Jets(15GeV)");
	
	sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
	HistoMan.GetPhoton2JetsHistograms(sub_folder,h2_PhoJets,"#gamma + >=2 Jets(15GeV)");

	sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
	HistoMan.GetTwoJetsHistograms(sub_folder,h2_TwoJets,"#gamma + >=2 Jets(15GeV)");




}


//_____________________________________________________________________________
int EleFakeRate::BeginJob()
{
	//begin job
	// check for all module dependencies and report errors.
	// initialize data members
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}
	
	tagEleMod = (TagElectrons*) ((TStnAna*) GetAna()->GetModule("TagElectrons"));
	if (!tagEleMod) {
		StdOut(__FILE__,__LINE__,3,"Electron tagging module required!.");
		sErrs.append("Electron Tagging mod not found. ");
		bRunPermit = false;
	}
	
	tagConvEleMod = (TagConvElectrons*) ((TStnAna*) GetAna()->GetModule("TagConvElectrons"));
	if (!tagConvEleMod) {
		StdOut(__FILE__,__LINE__,3,"Conversion Electrons tagging module required!.");
		sErrs.append("Conversion Electron Tagging mod not found. ");
		bRunPermit = false;
	}
	
	if (GetPhoenixBit()) {
		tagPhxMod = (TagPhoenixPhotons*) ((TStnAna*) GetAna()->GetModule("TagPhoenixPhotons"));
		if (!tagPhxMod) {
			StdOut(__FILE__,__LINE__,3,"Phoenix Photons tagging module required!.");
			sErrs.append("PhoenixBit set, but Phoenix Tagging mod not found. ");
			bRunPermit = false;
		}
	}

  	jetMod = (JetFilterModuleV2*) ((TStnAna*) GetAna()->GetModule("JetFilterV2"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		sErrs.append("JetFilter mod not found. ");
		bRunPermit = false;
	}

  	evtPropMod = (EventProperties*) ((TStnAna*) GetAna()->GetModule("EventProperties"));
	if (!evtPropMod) {
		StdOut(__FILE__,__LINE__,2,"EventProperties tag module not found.");
		sErrs.append("Event Properties mod not found. ");
		bRunPermit = false;
	}
	
 	// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fMetBlock 		= (TStnMetBlock*)      RegisterDataBlock("MetBlock","TStnMetBlock");

	BookHistograms();


	// initialize vars
/*	
	fSumFakeRate_1Jet = 0;
	fSumFakeRate_2Jet = 0;

	std::string file1("FakeRate_1Jet.sxt");
	ofFRlog_1Jet = new ofstream(file1.c_str());
	if (ofFRlog_1Jet->good()) {
		std::cout << "1 Jet Electron Fake Rate file opened" << std::endl;
		*ofFRlog_1Jet << "Fake rate\tFr-1sigma\tFr+1sigma\tsigma" <<std::endl;
	} else {
		StdOut(__FILE__,__LINE__,3,"1 Jet Fake Rate file open Error.");
		bRunPermit = false;
	}

	std::string file2("FakeRate_2Jet.sxt");
	ofFRlog_2Jet = new ofstream(file2.c_str());
	if (ofFRlog_2Jet->good()) {
		std::cout << "2 Jets Electron Fake Rate file opened" << std::endl;
		*ofFRlog_2Jet << "Fake rate\tFr-1sigma\tFr+1sigma\tsigma" <<std::endl;
	} else {
		std::cout << "2 Jet Fake Rate file open error" << std::endl;
		StdOut(__FILE__,__LINE__,3,"2 Jets Fake Rate file open Error.");
		bRunPermit = false;
	}
*/
	return 0;
}

//_____________________________________________________________________________
int EleFakeRate::BeginRun()
{

	//int currRun =  fHeaderBlock->RunNumber();
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int EleFakeRate::Event(int ientry)
{
	SetPassed(1);
	if (!bRunPermit) {
		StdErr(__FILE__,__LINE__,3,"One or more dependencies missing/ OR log file failed to open. pl. check!");
		exit (1);
		return 0;
	}
	
	hCount.iCounter->Fill(0);
	FillDataBlocks(ientry);

	//first look for a good electron. then look at number of jets
	//if there are several electrons, reject those events.
	/*
	std::vector<int> vElectronIndex;
	std::vector<double> vElectronFakeRate_d;
	std::vector<double> vElectronFakeRate_m;
	std::vector<double> vElectronFakeRate_p;
	std::vector<double> vElectronFakeRate_sigma;

	int NsuperPho = initSpMod->GetSuperPhoSize();
	for (int i=0; i < NsuperPho; i++) {
		if (! initSpMod->GetSuperPhoton(i)->IsTightElectron())  continue; // select electrons
		if (initSpMod->GetSuperPhoton(i)->IsConversion())  continue; // remove conversion
		//should be included later once the BH stuff is finalized
		//if (initSpMod->GetSuperPhoton(i)->IsBeamHalo())  continue; // remove conversion
		if (iPhoenixBit) {	
			if (initSpMod->GetSuperPhoton(i)->GetPhoenixId() != 0) continue;   // reject phoenix photons
		}
		
		float fEleId = initSpMod->GetSuperPhoton(i)->GetTightElectronId();
			
		double eleEt = initSpMod->GetSuperPhoton(i)->GetEtCorr();
		double wght_d = 0, wght_m =0, wght_p=0;
		MyFakeRate(eleEt, iGoodRunBit, iPhoenixBit, wght_d, wght_m, wght_p);
		if (wght_d > 0) {
			vElectronIndex.push_back(i);
			vElectronFakeRate_d.push_back(wght_d);
			vElectronFakeRate_m.push_back(wght_m);
			vElectronFakeRate_p.push_back(wght_p);
			vElectronFakeRate_sigma.push_back(fabs(wght_d - wght_m));
		}
	}
	//reject events with more tha 1 electron for now		
	if (vElectronIndex.size() == 1) {
		// >=1 jet case
		int iLeadJetIndex = 0;		
		int i2ndLeadJetIndex = 1;
		
		if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
			SetPassed(0);			//making exclusive templates
			hCount.iCounter->Fill(2);
			
			float  fEventFakeRate = 0;
			for (int i=0; i < vElectronIndex.size(); i++) {
				*ofFRlog_1Jet  << vElectronFakeRate_d[i] << "\t" 
									<< vElectronFakeRate_m[i] << "\t"
									<< vElectronFakeRate_p[i] << "\t"
									<< vElectronFakeRate_sigma[i] << std::endl;

				fEventFakeRate += vElectronFakeRate_d[i];
				// use individual fake rate of electrons to weight the plots like
				// invaraint mass
				// use sum of fake rate of electron to weight the plots like
				// MEt, Ht
				evtPropMod->FillPhotonHists(h1_Pho,vElectronIndex[i],vElectronFakeRate_d[i]);
				evtPropMod->FillPhoton1JetHists(h1_PhoJet,vElectronIndex[i],iLeadJetIndex,vElectronFakeRate_d[i]);
			}
			
			evtPropMod->FillEventHists(h1_Evt,fEventFakeRate);
			evtPropMod->FillJetHists(h1_Jet,iLeadJetIndex,fEventFakeRate);		//lead jet
	
			fSumFakeRate_1Jet += fEventFakeRate;
		} // 1 JET CASE
		
		
		// >=2 jet case
		if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
			SetPassed(0);			//making exclusive templates
			hCount.iCounter->Fill(3);
			float  fEventFakeRate = 0;
			for (int i=0; i < vElectronIndex.size(); i++) {
				*ofFRlog_2Jet  << vElectronFakeRate_d[i] << "\t" 
									<< vElectronFakeRate_m[i] << "\t"
									<< vElectronFakeRate_p[i] << "\t"
									<< vElectronFakeRate_sigma[i] << std::endl;

				fEventFakeRate += vElectronFakeRate_d[i];
				// use individual fake rate of electrons to weight the plots like
				// invaraint mass
				// use sum of fake rate of electron to weight the plots like
				// MEt, Ht

				evtPropMod->FillPhotonHists(h2_Pho,vElectronIndex[i],vElectronFakeRate_d[i]);
				evtPropMod->FillPhoton1JetHists(h2_PhoJet1,vElectronIndex[i],iLeadJetIndex,vElectronFakeRate_d[i]);
				evtPropMod->FillPhoton1JetHists(h2_PhoJet2,vElectronIndex[i],i2ndLeadJetIndex,vElectronFakeRate_d[i]);
				evtPropMod->FillPhoton2JetsHists(h2_PhoJets,i,iLeadJetIndex,i2ndLeadJetIndex,vElectronFakeRate_d[i]);
			}
			evtPropMod->FillEventHists(h2_Evt,fEventFakeRate);
			evtPropMod->FillJetHists(h2_Jet1,iLeadJetIndex,fEventFakeRate);		//lead jet
			evtPropMod->FillJetHists(h2_Jet2,i2ndLeadJetIndex,fEventFakeRate);		//2nd lead jet
			evtPropMod->FillTwoJetsHists(h2_TwoJets,iLeadJetIndex,i2ndLeadJetIndex,fEventFakeRate);

			fSumFakeRate_2Jet += fEventFakeRate;
			
		} // 2 JET CASE

	}
  
  */

	if (GetPassed()) {
		hCount.iCounter->Fill(1);
	}
	
	return 0;

} // Event


//_________________________________ returns fake rate & systematics
//_________________________________ Main function _______________________
//__________ input:  eleEt= Pho->Etc(); (where Pho must be a TStnPhoton object)
//		     goodrun_bit=0 --- for goodrun_v13_pho_00.txt good run list
//		     goodrun_bit=1 --- for goodrun_v13_pho_lepphx_00.txt good run list		     
//		     phnx_bit=0 --- if you DON'T require photon-to-phoenix rejection
//		     phnx_bit=1 --- if you DO require photon-to-phoenix rejection
//____________________________________________________________________________________________________
//__________ output: wght_d=fake rate; wght_m=(fake rate) - (1 sigma); wght_p=(fake rate) + (1 sigma) 

void EleFakeRate::MyFakeRate(double eleEt, int goodrun_bit, int phnx_bit,
											double &wght_d, double &wght_m, double &wght_p)
{
	wght_d=0.0;
	wght_m=0.0;
	wght_p=0.0;

	if (goodrun_bit==0 && phnx_bit==0) { // with phnx, all runs
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
	
	if (goodrun_bit==1 && phnx_bit==0) { // with phnx, good Silicon //sam: isnt this and the following case same
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
	
	//if (goodrun_bit==1 && phnx_bit==0) { // no phnx, all runs //sam::this shouldn't this be goodrun=0 phx=1 
	if (goodrun_bit==0 && phnx_bit==1) { // no phnx, all runs //sam::this shouldn't this be goodrun=0 phx=1 
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
	
	if (goodrun_bit==1 && phnx_bit==1) { // no phnx, good Silicon
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
	
	if (goodrun_bit<0 || goodrun_bit>1 || phnx_bit<0 || phnx_bit>1) return;

	if (eleEt>=13.0) {
		double term1 = MyFakeRateFuncErr(eleEt);
		double multp = (1.0 - MyPhnxEffFunc(eleEt));
		double term2 = MyPhnxEffFuncErr(eleEt) * MyFakeRateFunc(eleEt);
		double my_err = MyShiftCorrFunc(eleEt) * sqrt(term1 * term1 * multp * multp + term2 * term2);
		wght_d = MyShiftCorrFunc(eleEt) * MyFakeRateFunc(eleEt) * (1.0 - MyPhnxEffFunc(eleEt));
		wght_m = wght_d - my_err;
		wght_p = wght_d + my_err;
		if (wght_d < 0.0) wght_d = 0.0;
		if (wght_m < 0.0) wght_m = 0.0;
		if (wght_p < 0.0) wght_p = 0.0;
	}

	return;
}


//____________ correction for Et(det)/Et(genp) difference
double EleFakeRate::MyShiftCorrFunc(double eleEt) 
{
	double arg1 = 0.0048 + exp(-3.031-0.027 * eleEt);
	double arg2 = 0.0056 + exp(-3.015-0.030 * eleEt);
	if (arg2 == 0) StdOut(__FILE__,__LINE__,3,"Division by zero");
	double value = arg1 / arg2;
	return value;
}

//__________________________________ phoenix efficiency function
double EleFakeRate::MyPhnxEffFunc(double eleEt)
{
	double arg = phnx_par[1] * (phnx_par[2] - eleEt);
	double value = phnx_par[0] * TMath::Erfc(arg);
	if(value<0.0 || value>1.0) value = 0.0;
	return value;
}

//__________________________________ uncertainty on phoenix efficiency
double EleFakeRate::MyPhnxEffFuncErr(double eleEt)
{
	double arg = phnx_par[1] * (phnx_par[2]-eleEt);
	double erfc_drv = 2.0 * exp(-1.0 * arg * arg) / TMath::Pi();
	double term1 = phnx_parerr[0] * TMath::Erfc(arg);
	double term2 = phnx_parerr[1] * phnx_par[0] * (phnx_par[2]-eleEt) * erfc_drv;
	double term3 = phnx_parerr[2] * phnx_par[0] * phnx_par[1] * erfc_drv;
	double term4 = 2.0 * phnx_corr_coeff[0][1] * phnx_parerr[0] * phnx_parerr[1] * TMath::Erfc(arg)
						* erfc_drv * phnx_par[0] * (phnx_par[2]-eleEt);
	double term5 = 2.0 * phnx_corr_coeff[0][2] * phnx_parerr[0] * phnx_parerr[2] * TMath::Erfc(arg)
						* erfc_drv * phnx_par[0] * phnx_par[1];
	double term6 = 2.0 * phnx_corr_coeff[1][2] * phnx_parerr[1] * phnx_parerr[2]
						* arg * erfc_drv * erfc_drv * phnx_par[0] * phnx_par[0];
	double term = term1 * term1 + term2 * term2 + term3 * term3 + term4 + term5 + term6;
	double err = 0.0;
	if(term>0.0) err = sqrt(term);
	return err;
}


//__________________ returns fake rate
double EleFakeRate::MyFakeRateFunc(double eleEt)
{
	double fkrt = 0.0;
	if (eleEt >= 13.0) {
		double arg1 = fakerate_par[0] + fakerate_par[1] * eleEt;
		double lin1 = fakerate_par[2];
		double scale1 = fakerate_par[3];
		fkrt = scale1 * (exp(arg1) + lin1);
	}  
	return fkrt;
}

//________________________________________ returns fake rates
double EleFakeRate::MyFakeRateFuncErr(double eleEt)
{
	double err=0.0;
	double _expo=0.0;
	
	if (fakerate_par[3] > 0.0) {
		_expo=MyFakeRateFunc(eleEt)/fakerate_par[3];
		if (fakerate_par[3] == 0) StdOut(__FILE__,__LINE__,3,"Division by zero");
	}
	double mult1 = _expo - fakerate_par[2];
	double term1 = _expo * _expo * fakerate_parerr[3] * fakerate_parerr[3];
	double sub_term1 = mult1 * mult1 * (fakerate_parerr[0] * fakerate_parerr[0]
							+ 2.0 * fakerate_corr_coeff[1][0] * fakerate_parerr[0] * fakerate_parerr[1] * eleEt
							+ fakerate_parerr[1] * fakerate_parerr[1] * eleEt * eleEt);
	double sub_term2 = 2.0 * fakerate_parerr[2] * mult1 * (fakerate_corr_coeff[2][0] * fakerate_parerr[0]
							+ fakerate_corr_coeff[2][1] * fakerate_parerr[1] * eleEt);
	double sub_term3 = fakerate_parerr[2] * fakerate_parerr[2];
	double term2 = fakerate_par[3] * fakerate_par[3] * (sub_term1 + sub_term2 + sub_term3);
	if ((term1 + term2)>0.0) err = sqrt(term1 + term2);
	return err;
}


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void EleFakeRate::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void EleFakeRate::FillDataBlocks(int ientry)
{
	fMetBlock->GetEntry(ientry);	
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int EleFakeRate::EndJob() {
	
	ofFRlog_1Jet->close();
	ofFRlog_2Jet->close();
	

	std::string modname = GetName();
	std::string msg;
	msg  = "[EFR:00:]----- end job: ---- " + modname + "\n";
	msg += "[EFR:01:] Events Run Over ------------ = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	msg += "[EFR:02:] Events Pass this module ---- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	msg += "[EFR:03:] ele + 1 Jet Events --------- = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
	msg += "[EFR:04:] ele + 2 Jet Events --------- = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n";
	msg += "[EFR:05:] Phoenix Bit ---------------- = " + ToStr(GetPhoenixBit()) + "\n";
	msg += "[EFR:06:] Good run Bit --------------- = " + ToStr(GetGoodRunBit()) + "\n";
	msg += "[EFR:07:] Sum of Fake Rates (>=1 Jet)  = " + ToStr(GetSumFakeRate_1Jet()) + "\n";
	msg += "[EFR:08:] Sum of Fake Rates (>=2 Jet)  = " + ToStr(GetSumFakeRate_2Jet()) +"\n";
	msg += "---------------------------------------------------\n";
	
	std::cout << msg;
	if (sErrs.length() > 0) std::cout << sErrs << std::endl;
	
	return 0;
}
