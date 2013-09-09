///////////////////////////////////////////////////////////
// this is to same as Halo4_5Study except this look at   //
// all halo types with                                   //
// Met/PhoEt > 1 and <1 to figure out if a halo can      //
// produce a jet not only a photons.                     //
///////////////////////////////////////////////////////////
// Author: Samantha K. Hewamanage <samantha@fnal.gov>

#include "samantha/Pho/HaloStudyAllWithMetEt.hh"
#include "samantha/Pho/InitSuperPhotons.hh"
#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnPhoton.hh"
#include "samantha/utils/FreeFunctions.hh"
#include <iostream>

ClassImp(HaloStudyAllWithMetEt)

//_____________________________________________________________________________
HaloStudyAllWithMetEt::HaloStudyAllWithMetEt(const char* name, const char* title):
  TStnModule(name,title),
  qMc(false),
  bRunPermit(true),
  fEmTimeMin(-4.8),
  fEmTimeMax(4.8)
{
	// Default constructor
	std::cout << "Hello I am HaloStudyAllWithMetEt module" << std::endl;
}

//_____________________________________________________________________________
HaloStudyAllWithMetEt::~HaloStudyAllWithMetEt()
{
	// Default destructor
}

//_____________________________________________________________________________
void HaloStudyAllWithMetEt::SaveHistograms()
{
}

//_____________________________________________________________________________
void HaloStudyAllWithMetEt::BookHistograms()
{
	// Book all the histograms
	
	DeleteHistograms();
	
  	//char name [200];
  	//char title[200];
	TFolder *halofolder, *folder, *new_folder, *sub_folder;
	Histograms HistoMan;
	
	folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	hMetEtRatio = new TH1F("MetEtRatio","MEt/Pho Et", 100,0,5);
	folder->Add(hMetEtRatio);
	vCounterHistLabels.push_back("Events Processed");			//0
	vCounterHistLabels.push_back("Events Passed");				//1
	vCounterHistLabels.push_back("Type 0, Met/Pho Et>1");		//2
	vCounterHistLabels.push_back("Type 0, Met/Pho Et<1");		//3
	vCounterHistLabels.push_back("Type 1, Met/Pho Et>1");		//4
	vCounterHistLabels.push_back("Type 1, Met/Pho Et<1");		//5
	vCounterHistLabels.push_back("Type 2, Met/Pho Et>1");		//6
	vCounterHistLabels.push_back("Type 2, Met/Pho Et<1");		//7
	vCounterHistLabels.push_back("Type 3, Met/Pho Et>1");		//8	
	vCounterHistLabels.push_back("Type 3, Met/Pho Et<1");		//9
	vCounterHistLabels.push_back("Type 4, Met/Pho Et>1");		//10
	vCounterHistLabels.push_back("Type 4, Met/Pho Et<1");		//11
	vCounterHistLabels.push_back("Type 5, Met/Pho Et>1");		//12
	vCounterHistLabels.push_back("Type 5, Met/Pho Et<1");		//13
	HistoMan.GetCounterHistograms(folder,hCount,this->GetName(),vCounterHistLabels);

	Histograms::EventHists_t hEvt;
	Histograms::PhotonHists_t hPho;
	Histograms::JetHists_t hJet;
	Histograms::Photon1JetHists_t hPhoJet;
	Histograms::Photon2JetsHists_t hPhoJets;
	Histograms::TwoJetsHists_t hTwoJets;

	// for MEtEtRatio < 1  : halo type 0-5
	//folder structure
	// HaloType/MetEtRatio_gt1/1Jet/...

	for (int i =0; i <6; i++) {
		GA_Evt.push_back(hEvt);
		GB_Evt.push_back(hEvt);
		GA_Pho.push_back(hPho);
		GB_Pho.push_back(hPho);
		GA_Jet.push_back(hJet);
		GB_Jet1.push_back(hJet);
		GB_Jet2.push_back(hJet);
		GA_PhoJet.push_back(hPhoJet);
		GB_PhoJet1.push_back(hPhoJet);
		GB_PhoJet2.push_back(hPhoJet);
		GB_PhoJets.push_back(hPhoJets);
		GB_TwoJets.push_back(hTwoJets);

		LA_Evt.push_back(hEvt);
		LB_Evt.push_back(hEvt);
		LA_Pho.push_back(hPho);
		LB_Pho.push_back(hPho);
		LA_Jet.push_back(hJet);
		LB_Jet1.push_back(hJet);
		LB_Jet2.push_back(hJet);
		LA_PhoJet.push_back(hPhoJet);
		LB_PhoJet1.push_back(hPhoJet);
		LB_PhoJet2.push_back(hPhoJet);
		LB_PhoJets.push_back(hPhoJets);
		LB_TwoJets.push_back(hTwoJets);


		std::ostringstream foldername, text2title1, text2title2;
		foldername << "Type"<< i;
		text2title1 << "Halo("<< i << ")+>=2Jets";
		text2title2 << "Halo("<< i << ")+>=1Jets";
			
		halofolder = GetHistFolder(this,foldername.str().c_str(),foldername.str().c_str());
		folder = halofolder->AddFolder("MetEtRatio_gt1","Met/Photon Et Ratio > 1");
			new_folder = folder->AddFolder("1Jet","1 Jet case");
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
		
				HistoMan.GetEventHistograms(sub_folder,GA_Evt[i],text2title1.str().c_str());
				sub_folder = new_folder->AddFolder("Photon",text2title1.str().c_str());
				HistoMan.GetPhotonHistograms(sub_folder,GA_Pho[i],text2title1.str().c_str());

				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,GA_Jet[i],text2title1.str().c_str());
				sub_folder = new_folder->AddFolder("PhotonLeadJet",text2title1.str().c_str());
				HistoMan.GetPhoton1JetHistograms(sub_folder,GA_PhoJet[i],text2title1.str().c_str());


			new_folder = folder->AddFolder("2Jet","#gamma + >=2 Jet case");
	
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
				HistoMan.GetEventHistograms(sub_folder,GB_Evt[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
				HistoMan.GetPhotonHistograms(sub_folder,GB_Pho[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,GB_Jet1[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,GB_Jet2[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(sub_folder,GB_PhoJet1[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(sub_folder,GB_PhoJet2[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
				HistoMan.GetPhoton2JetsHistograms(sub_folder,GB_PhoJets[i],text2title2.str().c_str());

				sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
				HistoMan.GetTwoJetsHistograms(sub_folder,GB_TwoJets[i],text2title2.str().c_str());


		folder = halofolder->AddFolder("MetEtRatio_lt1","Met/Photon Et Ratio < 1");
			new_folder = folder->AddFolder("1Jet","1 Jet case");
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
		
				HistoMan.GetEventHistograms(sub_folder,LA_Evt[i],text2title1.str().c_str());
				sub_folder = new_folder->AddFolder("Photon",text2title1.str().c_str());
				HistoMan.GetPhotonHistograms(sub_folder,LA_Pho[i],text2title1.str().c_str());

				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,LA_Jet[i],text2title1.str().c_str());
				sub_folder = new_folder->AddFolder("PhotonLeadJet",text2title1.str().c_str());
				HistoMan.GetPhoton1JetHistograms(sub_folder,LA_PhoJet[i],text2title1.str().c_str());


			new_folder = folder->AddFolder("2Jet","#gamma + >=2 Jet case");
	
				sub_folder = new_folder->AddFolder("Event","Event Histograms");
				HistoMan.GetEventHistograms(sub_folder,LB_Evt[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon","#gamma Histograms");
				HistoMan.GetPhotonHistograms(sub_folder,LB_Pho[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("LeadJet","Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,LB_Jet1[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("SecondLeadJet","2^{nd} Lead Jet Histograms");
				HistoMan.GetJetHistograms(sub_folder,LB_Jet2[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("PhotonLeadJet","#gamma and Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(sub_folder,LB_PhoJet1[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
				HistoMan.GetPhoton1JetHistograms(sub_folder,LB_PhoJet2[i],text2title2.str().c_str());
	
				sub_folder = new_folder->AddFolder("Photon2Jets","#gamma + 2Jets Histograms");
				HistoMan.GetPhoton2JetsHistograms(sub_folder,LB_PhoJets[i],text2title2.str().c_str());

				sub_folder = new_folder->AddFolder("2Jets","2Jets Histograms");
				HistoMan.GetTwoJetsHistograms(sub_folder,LB_TwoJets[i],text2title2.str().c_str());
		
		}


}


//_____________________________________________________________________________
int HaloStudyAllWithMetEt::BeginJob()
{
	// Initialize and check for all dependecies
	
  	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (!initSpMod) {
		StdOut(__FILE__,__LINE__,3,"InitSuperPhotons module required!.");
		bRunPermit = false;
	}

  	jetMod = (JetFilterModuleV2*) ((TStnAna*) GetAna()->GetModule("JetFilterV2"));
	if (!jetMod) {
		StdOut(__FILE__,__LINE__,3,"JetFilter module required!.");
		bRunPermit = false;
	}

  	tightMod = (TagTightPhotons*) ((TStnAna*) GetAna()->GetModule("TagTightPhotons"));
	if (!tightMod) {
		StdOut(__FILE__,__LINE__,2,"Tight photon tag module not found.");
		bRunPermit = false;
	}

  	tagHaloMod = (TagBeamHalo*) ((TStnAna*) GetAna()->GetModule("TagBeamHalo"));
	if (!tagHaloMod) {
		StdOut(__FILE__,__LINE__,2,"Halo tag module not found.");
		bRunPermit = false;
	}
	
  	evtPropMod = (EventProperties*) ((TStnAna*) GetAna()->GetModule("EventProperties"));
	if (!evtPropMod) {
		StdOut(__FILE__,__LINE__,2,"EventProperties tag module not found.");
		bRunPermit = false;
	}

  	setTimeMod = (SetEmTimes*) ((TStnAna*) GetAna()->GetModule("SetEmTimes"));
	if (!setTimeMod) {
		StdOut(__FILE__,__LINE__,2,"SetEmTimes tag module not found.");
		bRunPermit = false;
	}


				// register the data blocks
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");

	BookHistograms();

	// initialize vars
	counter.evtsRunOver 			= 0;
	counter.evtsPassModule 		= 0;
	counter.haloEventFound = 0;

	return 0;
}

//_____________________________________________________________________________
int HaloStudyAllWithMetEt::BeginRun()
{

	//int currRun =  fHeaderBlock->RunNumber();
	//std::cout << " BEGINING RUN# " << currRun << std::endl;
	
	if (fHeaderBlock->McFlag()) 
   	qMc = true;
	else
   	qMc = false;
  
	return 0;

} //BeginRun

/*-------------------------------------------------------------------*/
// EVENT LOOP
/*-------------------------------------------------------------------*/
int HaloStudyAllWithMetEt::Event(int ientry)
{
	// Event loop
	
	SetPassed(1);
	
	// this makes sure that all required mods exists
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,1,"Permit Failed.");
		exit (1);
		return 0;
	}
		
	counter.evtsRunOver++;
	hCount.iCounter->Fill(0);
  	int NsuperPho = initSpMod->GetSuperPhoSize();
	
	for (int i=0; i < NsuperPho; i++) {
		if (! initSpMod->GetSuperPhoton(i)->IsTightPhoton()) continue;
		if (! initSpMod->GetSuperPhoton(i)->IsBeamHalo()) continue;		// select BH
		float fEmTime = initSpMod->GetSuperPhoton(i)->GetEmTime();
		if (fEmTime > fEmTimeMax || fEmTime < fEmTimeMin) continue;		//cut on em time, same as real photon template
		int iHaloId = initSpMod->GetSuperPhoton(i)->GetBeamHaloId();
			
		float MetEtRatio =  jetMod->GetMyMetCorr(0,3)/ initSpMod->GetSuperPhoton(i)->GetEtCorr();
		std::cout << "ratio=" << MetEtRatio << std::endl;
		hMetEtRatio->Fill(MetEtRatio);

		if (MetEtRatio >= 1.0) {
			for (int j=0; j <6 ; j++ ) {
				if ( ! ((iHaloId >> j) & 0x1) ) continue;   //ignore non halo

				int iBin = 2 + 2*j;
				hCount.iCounter->Fill(iBin);		//2 already in use counters + halo

				if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
					evtPropMod->FillPhotonHists(GA_Pho[j],i,1);
					evtPropMod->FillEventHists(GA_Evt[j],1);
					evtPropMod->FillJetHists(GA_Jet[j],0,1);		//lead jet
					evtPropMod->FillPhoton1JetHists(GA_PhoJet[j],i,0,1);
				}
				
				if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
					evtPropMod->FillPhotonHists(GB_Pho[j],i,1);
					evtPropMod->FillEventHists(GB_Evt[j],1);
					evtPropMod->FillJetHists(GB_Jet1[j],0,1);		//lead jet
					evtPropMod->FillJetHists(GB_Jet2[j],1,1);		//2nd lead jet
					evtPropMod->FillPhoton1JetHists(GB_PhoJet1[j],i,0,1);
					evtPropMod->FillPhoton1JetHists(GB_PhoJet2[j],i,1,1);
					evtPropMod->FillPhoton2JetsHists(GB_PhoJets[j],i,0,1,1);
					evtPropMod->FillTwoJetsHists(GB_TwoJets[j],0,1,1);
				}
			}

		} else {
			
			PrintInBinaryForm(iHaloId);
			for (int j=0; j <6 ; j++ ) {
				if ( ! ((iHaloId >> j) & 0x1) ) continue;   //ignore non halo
				int iBin = 2 + 2*j + 1;
				hCount.iCounter->Fill(iBin);		//2 already in use counters + halo

				if (jetMod->GetMyNjet(0,3) >= 1) {	// mono-jet case
					evtPropMod->FillPhotonHists(LA_Pho[j],i,1);
					evtPropMod->FillEventHists(LA_Evt[j],1);
					evtPropMod->FillJetHists(LA_Jet[j],0,1);		//lead jet
					evtPropMod->FillPhoton1JetHists(LA_PhoJet[j],i,0,1);
				}
				
				if (jetMod->GetMyNjet(0,3) >= 2) {	// di-jet case
					evtPropMod->FillPhotonHists(LB_Pho[j],i,1);
					evtPropMod->FillEventHists(LB_Evt[j],1);
					evtPropMod->FillJetHists(LB_Jet1[j],0,1);		//lead jet
					evtPropMod->FillJetHists(LB_Jet2[j],1,1);		//2nd lead jet
					evtPropMod->FillPhoton1JetHists(LB_PhoJet1[j],i,0,1);
					evtPropMod->FillPhoton1JetHists(LB_PhoJet2[j],i,1,1);
					evtPropMod->FillPhoton2JetsHists(LB_PhoJets[j],i,0,1,1);
					evtPropMod->FillTwoJetsHists(LB_TwoJets[j],0,1,1);
				}//if 
			}//for
			
		}//if  
		break; // don't care about rest of the photons
					//if not i'll be putting multiple entries into met.sumet plots
	} //for

	if (GetPassed()) {
		counter.evtsPassModule++;
		hCount.iCounter->Fill(1);
	}
	
	return 0;

} // Event


/*-------------------------------------------------------------------*/
// clean and clear things for the run
/*-------------------------------------------------------------------*/
void HaloStudyAllWithMetEt::Cleanup()
{

} //Cleanup


/*-------------------------------------------------------------------*/
void HaloStudyAllWithMetEt::FillDataBlocks(int ientry)
{
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int HaloStudyAllWithMetEt::EndJob() {

	std::string sModName(GetName());
	std::string sMsg;
	
	sMsg  = "[HAL:00:]----- end job: ---- " + sModName + "\n";
	sMsg += "[HAL:01:] Events Processed ----------- = " + ToStr(counter.evtsRunOver) + "\n";
	sMsg += "[HAL:02:] Events Pass this module ---- = " + ToStr(counter.evtsPassModule) + "\n";
	sMsg += "[HAL:03:] ==  MEt/PhoEt > 1 ===\n";
	sMsg += "[HAL:04:] Halo Type 0 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(3)) + "\n";
	sMsg += "[HAL:05:] Halo Type 1 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(5)) + "\n";
	sMsg += "[HAL:06:] Halo Type 2 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(7)) + "\n";
	sMsg += "[HAL:07:] Halo Type 3 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(9)) + "\n";
	sMsg += "[HAL:08:] Halo Type 4 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(11)) + "\n";
	sMsg += "[HAL:09:] Halo Type 5 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(13)) + "\n";
	sMsg += "[HAL:10:] ==  MEt/PhoEt < 1 ===\n";
	sMsg += "[HAL:11:] Halo Type 0 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(4)) + "\n";
	sMsg += "[HAL:12:] Halo Type 1 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(6)) + "\n";
	sMsg += "[HAL:13:] Halo Type 2 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(8)) + "\n";
	sMsg += "[HAL:14:] Halo Type 3 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(10)) + "\n";
	sMsg += "[HAL:15:] Halo Type 4 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(12)) + "\n";
	sMsg += "[HAL:16:] Halo Type 5 ---------------- = " + ToStr(hCount.iCounter->GetBinContent(14)) + "\n";

	sMsg += "---------------------------------------------------\n";
	std::cout << sMsg;
	return 0;
}
