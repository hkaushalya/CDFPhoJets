///////////////////////////////////////////////////////////
// this is my e + jets template same as other temps      //
// except select electrons.                              //
// 1. electron tagging module must exitst                //
// 2. conversion removal module must exist               //
// 3. if two electrons are found sum of the two fake     //
//     rates will be retuned.                            //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "samantha/Pho/ElectronTemplate.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include <iostream>
#include "samantha/Pho/EleFakeRate.hh"

ClassImp(ElectronTemplate)
//_____________________________________________________________________________
ElectronTemplate::ElectronTemplate(const char* name, const char* title):
  TStnModule(name,title),
  sErrs(""),
  bRunPermit(true),
  fEmTimeMax(4.8),
  fEmTimeMin(-4.8)
{
	//constructor
	std::cout << "Hello I am ElectronTemplate module" << std::endl;
}

//_____________________________________________________________________________
ElectronTemplate::~ElectronTemplate() {
	//destructor
}

//_____________________________________________________________________________
void ElectronTemplate::SaveHistograms() {
}

//_____________________________________________________________________________
void ElectronTemplate::BookHistograms()
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
int ElectronTemplate::BeginJob()
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
		StdOut(__FILE__,__LINE__,3,"Tight Electron tagging module required!.");
		sErrs.append("Electron Tagging mod not found. ");
		bRunPermit = false;
	}
	
	tagConvEleMod = (TagConvElectrons*) ((TStnAna*) GetAna()->GetModule("TagConvElectrons"));
	if (!tagConvEleMod) {
		StdOut(__FILE__,__LINE__,3,"Conversion Electrons tagging module required!.");
		sErrs.append("Conversion Electron Tagging mod not found. ");
		bRunPermit = false;
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
	
  	looseEleTemp = (LooseElectronTemplate*) ((TStnAna*) GetAna()->GetModule("LooseElectronTemplate"));
	if (looseEleTemp) {
	   TList* modList = GetAna()->GetModuleList();
		if (modList->GetSize() > 0) {
			TIterator* it = modList->MakeIterator();
			TStnModule* tstnMod;
			std::string lEleTempName(looseEleTemp->GetName());
			std::string tEleTempName(GetName());
			//bool bFoundTEleTemp = false;
			
			while (tstnMod = (TStnModule*) it->Next()) {
				std::string sModName(tstnMod->GetName());
				//if () {
				//}
			}
		tstnMod =0;
		}
																			  
	} else {
		StdOut(__FILE__,__LINE__,2,"LooseElectronTemplate not found. it controls pass/fail for this mod!");
		sErrs.append("LooseElectronTemplate not found. ");
		bRunPermit = false;
	}
	
 	// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();


	// initialize vars
	fSumFakeRate_1Jet = 0;
	fSumFakeRate_2Jet = 0;

	std::string file1("TightEleFakeRate_1Jet.sxt");
	ofFRlog_1Jet = new ofstream(file1.c_str());
	if (ofFRlog_1Jet->good()) {
		std::cout << "1 Jet Electron Fake Rate file opened" << std::endl;
		*ofFRlog_1Jet << "Fake rate\tFr-1sigma\tFr+1sigma\tsigma" <<std::endl;
	} else {
		StdOut(__FILE__,__LINE__,3,"1 Jet Fake Rate file open Error.");
		bRunPermit = false;
	}

	std::string file2("TightEleFakeRate_2Jet.sxt");
	ofFRlog_2Jet = new ofstream(file2.c_str());
	if (ofFRlog_2Jet->good()) {
		std::cout << "2 Jets Electron Fake Rate file opened" << std::endl;
		*ofFRlog_2Jet << "Fake rate\tFr-1sigma\tFr+1sigma\tsigma" <<std::endl;
	} else {
		std::cout << "2 Jet Fake Rate file open error" << std::endl;
		StdOut(__FILE__,__LINE__,3,"2 Jets Fake Rate file open Error.");
		bRunPermit = false;
	}

	return 0;
}

//_____________________________________________________________________________
int ElectronTemplate::BeginRun()
{
	return 0;
} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int ElectronTemplate::Event(int ientry)
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
		//if (initSpMod->GetSuperPhoton(i)->IsBeamHalo())  continue;
		float fEmTime = initSpMod->GetSuperPhoton(i)->GetEmTime();
		if (fEmTime > fEmTimeMax || fEmTime < fEmTimeMin) continue;		//cut on em time, same as real photon template
		vElectronIndex.push_back(i);
		double eleEt = initSpMod->GetSuperPhoton(i)->GetEtCorr();
		double wght_d = 0, wght_m =0, wght_p=0;
		EleFakeRate aEleFR;
		aEleFR.MyFakeRate(eleEt, GetGoodRunBit(), GetPhoenixBit(), wght_d, wght_m, wght_p);
		vElectronFakeRate_d.push_back(wght_d);
		vElectronFakeRate_m.push_back(wght_m);
		vElectronFakeRate_p.push_back(wght_p);
		vElectronFakeRate_sigma.push_back(fabs(wght_d - wght_m));
	}
	
	//reject events with more tha 1 electron for now		
	//if (vElectronIndex.size() == 1) { //this could be wrong. as we would be looking at W electrons only
												// we want a template that includes all e's from Zs etc???
	for (unsigned int i=0; i < vElectronIndex.size(); i++) {
		int iLeadJetIndex = 0;		
		int i2ndLeadJetIndex = 1;
		
		// >=1 jet case
		if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
			SetPassed(0);			//making exclusive templates
			hCount.iCounter->Fill(2);
			*ofFRlog_1Jet  << vElectronFakeRate_d[i] << "\t" 
								<< vElectronFakeRate_m[i] << "\t"
								<< vElectronFakeRate_p[i] << "\t"
								<< vElectronFakeRate_sigma[i] << std::endl;
			fSumFakeRate_1Jet += vElectronFakeRate_d[i];
			
			evtPropMod->FillPhotonHists(h1_Pho,vElectronIndex[0],1);
			evtPropMod->FillPhoton1JetHists(h1_PhoJet,vElectronIndex[0],iLeadJetIndex,1);
			evtPropMod->FillEventHists(h1_Evt,1);
			evtPropMod->FillJetHists(h1_Jet,iLeadJetIndex,1);		//lead jet
	
		} // 1 JET CASE
		
		// >=2 jet case
		if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
			SetPassed(0);			//making exclusive templates
			hCount.iCounter->Fill(3);
			*ofFRlog_2Jet  << vElectronFakeRate_d[i] << "\t" 
								<< vElectronFakeRate_m[i] << "\t"
								<< vElectronFakeRate_p[i] << "\t"
								<< vElectronFakeRate_sigma[i] << std::endl;
			fSumFakeRate_2Jet += vElectronFakeRate_d[i];

			evtPropMod->FillPhotonHists(h2_Pho,vElectronIndex[0],1);
			evtPropMod->FillPhoton1JetHists(h2_PhoJet1,vElectronIndex[0],iLeadJetIndex,1);
			evtPropMod->FillPhoton1JetHists(h2_PhoJet2,vElectronIndex[0],i2ndLeadJetIndex,1);
			evtPropMod->FillPhoton2JetsHists(h2_PhoJets,vElectronIndex[0],iLeadJetIndex,i2ndLeadJetIndex,1);
			evtPropMod->FillEventHists(h2_Evt,1);
			evtPropMod->FillJetHists(h2_Jet1,iLeadJetIndex,1);		//lead jet
			evtPropMod->FillJetHists(h2_Jet2,i2ndLeadJetIndex,1);		//2nd lead jet
			evtPropMod->FillTwoJetsHists(h2_TwoJets,iLeadJetIndex,i2ndLeadJetIndex,1);

		} // 2 JET CASE

	}

	if (GetPassed()) {
		hCount.iCounter->Fill(1);
	}

	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void ElectronTemplate::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void ElectronTemplate::FillDataBlocks(int ientry)
{
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int ElectronTemplate::EndJob() {
	ofFRlog_1Jet->close();
	ofFRlog_2Jet->close();
	
	std::string modname = GetName();
	std::string msg;
	msg  = "[ETP:00:]----- end job: ---- " + modname + "\n";
	msg += "[ETP:01:] Events Run Over ------------ = " + ToStr(hCount.iCounter->GetBinContent(1)) + "\n";
	msg += "[ETP:02:] Events Pass this module ---- = " + ToStr(hCount.iCounter->GetBinContent(2)) + "\n";
	msg += "[ETP:03:] ele + 1 Jet Events --------- = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
	msg += "[ETP:04:] ele + 2 Jet Events --------- = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n";
	msg += "[ETP:05:] EM Time window used -------- = " + ToStr(GetEmTimeWindowMin()) + "," + ToStr(GetEmTimeWindowMax())+ "\n";
	msg += "[ETP:06:] Phoenix Bit ---------------- = " + ToStr(GetPhoenixBit()) + "\n";
	msg += "[ETP:07:] Good run Bit --------------- = " + ToStr(GetGoodRunBit()) + "\n";
	msg += "[ETP:08:] Sum of Fake Rates (>=1 Jet)  = " + ToStr(GetSumFakeRate_1Jet()) + "\n";
	msg += "[ETP:09:] Sum of Fake Rates (>=2 Jet)  = " + ToStr(GetSumFakeRate_2Jet()) +"\n";
	msg += "---------------------------------------------------\n";
	
	std::cout << msg;
	if (sErrs.length() > 0) std::cout << sErrs << std::endl;
	
	return 0;
}
