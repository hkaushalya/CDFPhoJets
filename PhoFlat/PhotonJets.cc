/*{{{*/
/*  $Id: PhotonJets.cc,v 1.26 2011/05/26 19:40:08 samantha Exp $
 *  $Log: PhotonJets.cc,v $
 *  Revision 1.26  2011/05/26 19:40:08  samantha
 *  MODIFIED: Almost back to regular g+jets study. I commented out 1jet hist fills
 *  and b-tag jet info to get the estimate for Wg and Zg events for Buckley.
 *
 *  Revision 1.25  2011/05/02 17:24:32  samantha
 *  ADDED: To check for b-tag jets. This version regquires on
 *  one b-tag jets. I have made plots with 2-btag jets. see elog:1792
 *
 *  Revision 1.24  2011/04/26 23:42:09  samantha
 *  MODIFIED:  1. Changed bUseDataMcNvtxWgts -> iUseDataMcNvtxWgts to get
 *  uncertainties for the nvtx weighting process.
 *  ADDED:  1. new Z+Jets weights
 *  2. WW/WZ/ZZ/TTbar nvtx weights
 *  3. ExoticEvent() method to search for fancy stuff.
 *  MINOR CHANGES: 1. Cleaned up some unused/old pieces of code lying around.
 *  2. Deleted some commented out stuff in BadMet()
 *
 *  Revision 1.23  2010/08/25 15:25:43  samantha
 *  BUG FIX: I have beeing ommiting the data beyong P13. Fixed this. However this
 *  does not affect my final results ebcause of the way things are normalized.
 *
 *  MAJOR CHANGES: 1. Nvtx andsideband weight function/values changed for plots with
 *  met clean up cuts. MEt is required to be away from any jet Et>15GeV. Do not care
 *  about the photon. Moved the met cut to the strip cuts.
 *  2. EWK  nvtx weights come from W+jets and Z+jets DATA/MC comparison.
 *
 *  Revision 1.22  2010/07/26 21:05:20  samantha
 *  MAJOR CHANGES: Added all QCD ID systematic to here.
 *  ADDED: 1. fMaxMet cut of 1500.GeV default.
 *         2. Few more plots to look at Vtx distribution with MET.
 *         3. FailStripCuts() to speedup things a bit.
 *         4. Many funtions to do debugging and to understand some events.
 *         CheckObjOverlap(), CheckJetTowers(), FoundPartialJet().
 *
 *  Revision 1.21  2010/04/13 16:29:12  samantha
 *  I was not able to read the new MetX,MetY, RawMet and CprWgt branches. ROOT
 *  complained that it can't find the branch to esxecute "SetBranchStatus" to 1.
 *  This method was caclled within the ReadInAll. After many tests I gave up.
 *  Following day, eveything worked just fine!
 *
 *  ADDED:1. more options to choose different reweighing methods. This includes
 *  photon Et only reweighing now (3 and 4)
 *  2. MetX and MetY info to the job. This info is needed to apply the met clean up
 *  cuts. Also these hists are added to the EventHists and DelPhi plots with Met/Pho
 *  and Met/jet is added to PhotonHist and JetHist.
 *  3.Modified the BadMet() (Met clean up cut) to be included in the job.
 *
 *  Revision 1.20  2010/04/07 17:57:21  samantha
 *  MODIFIED: 1. Dropped the fabs() from SetMin/MaxJetEta() in hh file. I think this caused
 *  some problems when I was applying cuts to jets in the JetList.
 *  ADDED: 1. nvtx based weighing, GetDATAMCNvtxWgt().
 *  2. The dataset used for the job can be passed in using SetDataset(). I added
 *  this to get the correct Nvtx weight function for each of the MC samples from a
 *  root file (this root file has weight functions for all MC samples). This can be
 *  used to automcatically execute methods speicific for certain dataset.
 *  3. Enabled the ApplyCuts method in JetList in the JetSelAndHistFill().
 *  4. Fixed a minor bug in parameter validation in the FillPhoton1JetHists()
 *  method. The upper bound of the photon index should be checked as '>' not with
 *  '>='. Same is true with the jetIndex.
 *  5. FillPhoton1JetHists() dumps events with InvMass of pho and jet <5GeV. This
 *  seems like I am adding the same obj. Which seems like I am not removing some of the
 *  objects from the jet list correctly OR I make mistake when additing extra em
 *  objects back to the jet list.! need to check this.
 *
 *  Revision 1.19  2010/03/31 23:12:42  samantha
 *  MAJOR CHANGES: Modified to reweight the sideband based on photon Et and jet eta
 *  weights. I tried to use a fit function to derive jet eta weights. For some
 *  reason it showed up nasty. see elog#1599 for details. So I decided to use the
 *  exact bin values instead. And this time it cam up right! elog#1602. I have
 *  commented out the par that I did the first test in GetJetEtaWeight() and
 *  included the stuff to use bin values.
 *
 *  BUG FIX: I found out that sideband EWK/PHO MC hists did not had photon Et cut
 *  applied. Also I was passing in the wrong sideband photon index to JetSelAndHistFill for
 *  MC jobs. I was using central photon index instead of pho_up or down in
 *  respective cases. No major changes in the final plots though.
 *
 *  ADDED: 1. Aditional libraries for exception handling. For debugging only.
 *  2. fMaxJetEt and fMinJetEta.
 *  3. PhotonJets::GetJetEtaWeight() to generate nominal and error weights for
 *  reweighing. GetHistBin() find the corresponding bin values for the jet eta, to get the
 *  weights.
 *  4. bMetCleanup option to apply met cleanup cuts to see if that will improve the
 *  MET plot. Also added these methods for this and they needs to be tested with new
 *  V7 stuples which has full MET vector info.BadMet()
 *  5. I have some commented out stuff in FillEventHists() where I tried to figure
 *  out why my Ht/Met calculation differs from what I have saved from MET model. I
 *  figured out the difference and see elog#1596 for my MET definition.
 *  6. The weight function are cleanly derived and read from the files. Delete all
 *  leftovers from previous jobs.
 *
 *  Revision 1.18  2010/03/24 19:05:42  samantha
 *  Modified for photon Et>30+jets sideband reweighing using photon Et and jet  eta
 *  weights.
 *  ADDED: 	1. DoHaloRejHists(), DoHaloRejCalc()  to get the halo rejection power
 *  	when run over data.;
 *  	2. GetCosmicEstAndErr() to get comsic estimated when run over data.
 *
 *  	for above calculation I had added additional sets of folder/hists toe
 *  	the job.
 *  	---
 *  	TH1F *haloPhiWedge_1j[2], *haloEmTime_1j[2];    //0=before 1= after BH cuts
 *  	TH1F *haloPhiWedge_2j[2], *haloEmTime_2j[2];
 *  	TH1F *haloEt_1j[2], *haloEt_2j[2];
 *  	TH1F *haloNvtx[2];
 *
 *  	Histograms::EventHists_t HhaloNoVtx_Evt_b4, HhaloNoVtx_Evt_a4;
 *  	Histograms::PhotonHists_t HhaloNoVtx_Pho_b4, HhaloNoVtx_Pho_a4;
 *  	Histograms::EventHists_t HhalosidebandNoVtx_Evt_b4, HhalosidebandNoVtx_Evt_a4;
 *  	Histograms::PhotonHists_t HhalosidebandNoVtx_Pho_b4, HhalosidebandNoVtx_Pho_a4;
 *
 *  	4. GetJetEtaWeight(const TF1* func, const float x) which derives an
 *  	addition weight for the sideband reweighing wbase on lead jet detector
 *  	eta. see elog#1590. For this I had to modify the JetSelAndHist fill
 *  	method.
 *
 *  Revision 1.17  2010/03/13 19:28:19  samantha
 *  MODIFIED: tf1_sideband function uses g40jetsmet40 function for the g30jetsmet40
 *  sample. I have modified the function name to include the function format and
 *  the parameters are printed out in the job summary.
 *  DELETED: retrieving of reweighing function based on di-jet sideband is removed
 *  as I am not using it.
 *
 *  Revision 1.16  2010/03/09 17:07:31  samantha
 *  Comments for the previous commit: v.1.15
 *  1. Added DEBUG() method to print out messages during debug process.
 *  2. Modified: DoMyStuff() vtx cut. No I can specify the min/max vertices needed.
 *  3. Modified: JetSelAndHistFill(): I have commented out the the photon-jet
 *  balance stuff.
 *
 *  Revision 1.15  2010/03/09 16:43:39  samantha
 *  PhoFlat/QCDSystematics.hh
 *
 *  Revision 1.14  2010/02/13 04:33:13  samantha
 *  ADDED: 1. Reweight function fitFunction3() for the iso added sideband.
 *         2. Events with Ht<fMinPhoEt+fMinJetEt are thrown away. this is to
 *         get rid of some events which show Ht below the sum of the objects.
 *         3.Events weights arerecorded in hist in photon list now.
 *  MODIFIED: 1. phoEt_iso25, phoEt_50 and phoEt_70 trigger photon Et hists
 *  bin size is reduced to 1GeV.
 *            2. The AddUnused() method in NewJetList method was upadated to
 *  	  give a jet et/eta cut. so those methods in here are updated.
 *  RENAMED: SelectSidebandPhoton() -> SelectMCSidebandPhoton()
 *
 *  DELETED: 1. DoZJetsStuff() that is no longer in need.
 *           2. pho-jet balance plot stuff.
 *  	 3. GetHt() method remnants.
 *
 *  Revision 1.13  2010/01/12 21:17:38  samantha
 *  MAJOR CHANGES: Modified sideband selection to study isolation addition to
 *  sideband photon and reweigh it to match the MC. See elog#1533,1536
 *  Jet selection is modified to for this to veto extra jets and require
 *  DelPhi(p,j)>2.7 .
 *  ADDED: Lots of changes with all the testing.
 *  	1. Vtx cut and Et cuts are reapplied for each type of event selection.
 *  	2. All the settings, Photon Et, MinNjet15 etc are stored in variables
 *  	and used in the code. Even the sideband reweighing, Isolation addition
 *  	can be controlled externally.
 *  	3. Those settings will be printed out by PrintJobSettings() at the
 *  	beginning and end of the job.
 *  	4. Histograms for beam halo and cosmics in the sideband. Also added
 *  	hists to subtract EWK component in the photon sideband.
 *
 *  Revision 1.12  2009/12/29 00:15:54  samantha
 *  MODIFIED: 1. Event selection is modified as follows. If there is a tight photon,
 *  	it will be first tested for CosmicTemplate, then Halo template. If it did not
 *  	fit either, will go into signal bin. If the event has a tight photon and a
 *  	sideband photon, it will go into signal. Sideband bin and signal bin are now
 *  	mutually exclusive.
 *  	2. Init() method: code is modifed to reweight the SIDEBAND (see
 *  	elog#1520, sideband photon Et is matched to total background photon Et
 *  	and used 100% to make the final plots, see elog#1526)
 *  	The weight function is read in from a root file. I get the parameters from
 *  	this read-in and creates a new one (tf1_sideband) and sets those values to
 *  	reconstruct the funtion. I am doing this to avoid any errorS in forming the
 *  	function (see Elog#1507). This reweighed sideband is used to make
 *  DELETED: The previous weight function in Init()
 *  ADDED  : 1. fitFunction2() and fitFunction() are the two set of sideband reweighing
 *  	function I have used.
 *  	2. Main() will check to make sure the MC flag set matches to the data
 *  	set run over. (MC flag required to set manually for each job)
 *  	3. Working on a timer which can tell time to finish. But it seems the
 *  	ROOT CINT does not quite like the C++ clock processing. I testes thie
 *  	clock funcntion C++ code and it worked. But same code does not work in
 *  	here.
 *  	4. RemoveDuplicateEMobjects() method is added to over come a special
 *  	case where a sideband photon also passed loose electron ID cuts and
 *  	added to both PhotonList and Electron list. This electron may end up as
 *  	the leading jet when the new jet list  created. So we are adding the
 *  	same object at the end.
 *
 *  Revision 1.11  2009/12/07 21:06:05  samantha
 *  MODIFIED: To derive the true photon fraction using CPR weights. Commented out
 *  parts of the photon hist filling routing for this.
 *  See elog#https://hep.baylor.edu/hep/samantha/1498
 *  ADDED: Few comments to avoid making mistakes of rearranging or commenting things
 *  that could lead to error, like appylying the vertex cut at the top will
 *  contradict with Beam Halo selection (meaning leave no events for BH)
 *  COMMENTED OUT: some of the variables not used in Hist fillings, like JetEmFr
 *  etc.
 *
 *  Revision 1.10  2009/11/28 00:07:03  samantha
 *  MODIFIED: To make a systematic calculation for photon sideband based on a di-jet
 *  MC study (elog#https://hep.baylor.edu/hep/samantha/1459). All the sideband
 *  kinematic plots are weghted based on the photon Et distribution in #1459 (using
 *  the linear fit) and compared to the nominal. Results are in
 *  https://hep.baylor.edu/hep/samantha/1485.
 *  Commented out the CPR weights that I thought may no longer need.
 *  Added auto CVS log into the file.
 *  PROBLEM: I noticed that I am selected fewer events (thousands or more) I used
 *  to. Do not remember making such a big change. Need to investigate.
 *
 *
 */

/*}}}*/

#include <iostream>
#include <sstream>
#include "PhoFlat/PhotonJets.hh"
#include "PhoFlat/ReadInAll.hh"
#include "PhoFlat/HaloTypes.hh"
#include "PhoFlat/NewJetList.hh"
#include "TLorentzVector.h"
#include "PhoFlat/HistManager.hh"
#include "TBenchmark.h"
#include "TCanvas.h"
#include <iomanip>
#include "TROOT.h"
//for debugging
#include <exception>
#include <stdexcept> // out_of_range exception

void DEBUG(const std::string func, const double line, const std::string msg)
{
	std::cout << func << ":" << line << ": " << 	msg << std::endl;
}
 
Double_t fitFunction2(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
		val = (par[0]+par[1]*sqrt(x[0]))+par[2]*x[0];
	}
	return val;
}

Double_t fitFunction(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
		double pol1 = par[0]  + par[1] * x[0] * x[0]; 
		double pol2 = par[2]  + par[3] * x[0];

		double wgt = 1 / (1 + exp((x[0]-par[4])/10.5));
		val = pol1 * wgt + pol2 * (1 - wgt);
	}
	return val;
}

//reweight iso-added sideband to total background
Double_t fitFunction3(Double_t *x, Double_t *par)
{
   double val = 0;
   if (x[0]>0)
   {
      val = (par[0]+par[1]*sqrt(x[0]))+par[2]*x[0];
   }
   return val;
}



//-------------------------------------------------------------------
PhotonJets::PhotonJets():
iProgressBy(50000),
sHistFileName("PhotonJets.root"),
bMcFlag(true),
bAddIsoToSidebandPho(false),
bRequireTagJets(false)
{

	fMinPhotonEt = 15.; // GeV  //do not lower this any more. See header file Set method for more details!
	fMaxPhotonEt = 1500.; // GeV
	fMinJetEt = 15.;
	fMaxJetEt = 1500.;
	fMaxPhotonEta = 1.1;
	fMinJetEta    = 0.0;
	fMaxJetEta    = 3.0;
	fMinSignalEmTime = -4.8; //ns
	fMaxSignalEmTime = 4.8;
	fMinCosmicEmTime = 30.;
	fMaxCosmicEmTime = 90.;
	iMinClass12Vtx = 1;
	iMaxClass12Vtx = 100;
	fMaxVtxz = 60.; //cm
	iHaloType = 5;
	iMinNjet15 = 1;
	iMaxNjet15 = 100;
	fRewgtSideband = false;
	fMinMet = 0.0;	//minimum Missing transverse energy
	fMaxMet = 1500.0;	//maximum Missing transverse energy
	bMetCleanup = false;		// no met clean up cuts by default
	iUseDataMcNvtxWgts = 0;  //0=no weghting, 1=central,2=central+1sigma
	
}

//-------------------------------------------------------------
void PhotonJets::Init(TChain *tree)
{
	readIn = new ReadInAll(tree, &stuple);
	readIn->CentralPhotons(1);
	readIn->CentralElectrons(1);
	readIn->CentralJets(1);
	readIn->Met(1);

	if (bMcFlag) {
		readIn->EmUpPhotons(1);
		readIn->EmUpElectrons(1);
		readIn->EmDownPhotons(1);
		readIn->EmDownElectrons(1);
		readIn->JesUpJets(1);
		readIn->JesDownJets(1);
		readIn->GenLevel(0);
		
		for (int i=0; i<4; i++) {
			vMcCount.push_back(0);
			vMcUpCount.push_back(0);
			vMcDownCount.push_back(0);
			vMcSidebandEWKCount.push_back(0);
		}
	} else {
			//counter assigninig
			// 0 = all events with at least a qualifying photon
			// 1 = event with a qualifying photon and >=1 jet 
			// 2 = events with a qualifying photon and >=2 jets
			// 3 = fails met clean up

		
		for (int i=0; i<4; i++) {
			vHaloCount.push_back(0);
			vCosmicCount.push_back(0);
			vQCDCount.push_back(0);
			vSignalCount.push_back(0);
			vSidebandHaloCount.push_back(0);
			vSidebandCosmicCount.push_back(0);
			vSidebandSystHadEmCount.push_back(0);
			vSidebandSystIsoCount.push_back(0);
			vSidebandSystTrkPtCount.push_back(0);
			vSidebandSystTrkIsoCount.push_back(0);
		}

	}
	
	iRejEvts = 0;
	
	readIn->Init();

	if (GetReweightSideband())
	{
		/* *************************************************************
		 * I tried to go high tech by storing and reading in the stored
		 * TF1 from a root file but it did not work. For some reason the
		 * function returns 0 for Et< some value.
		 * So now I am only reading the function to get the final fit
		 * values. I create a new function inside this module and assign
		 * those read parameters to the new function. -02-09-2010
		 ***************************************************************/

		std::cout << " ####### SIDEBAND REWEIGHING INFO ########## " << magenta << std::endl;

		if (GetReweightSideband() == 3 || GetReweightSideband() == 4) //photon Et weights for sideband reweight
		{
			//TFile *fSidebandSyst = new TFile("RootFiles/g30Exc1jet_SidebandPhoEtReweightFunctions.root");
			//TFile *fSidebandSyst = new TFile("~/ICHEP_LaptopVersion/p1_26/g30jets_p26_SidebandPhoEtReweightFunctions.root");
			//TFile *fSidebandSyst = new TFile("~/ICHEP_LaptopVersion/p1_26/NvtxWgtMCs/g30jets_SidebandJetEtaReweightFunctions.root");
			//TFile *fSidebandSyst = new TFile("~/ICHEP_LaptopVersion/p1_26/g30jets_p26_SidebandPhoEtReweightFunctions.root");
			//TFile *fSidebandSyst = new TFile("g30jets_p26_SidebandPhoEtReweightFunctions.root");
			//TFile *fSidebandSyst = new TFile("~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned/g30jetsNvtxWgtMC_MetCleaned_SidebandPhoEtReweightFunctions.root");
			//TFile *fSidebandSyst = new TFile("~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleanJetMetOnly_gEta9/g30jets_SidebandPhoEtReweightFunctions.root");
			//TFile *fSidebandSyst = new TFile("~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_final/g30jets_SidebandPhoEtReweightFunctions.root");
			TFile *fSidebandSyst = new TFile("~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_Met20_final/g30jets_SidebandPhoEtReweightFunctions.root");
			//TFile *fSidebandSyst = new TFile("~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_Met30_final/g30jetsMet30_PhoMCnvtxWgtd_SidebandPhoEtReweightFunctions.root");
			assert (!fSidebandSyst->IsZombie() && "Sideband Wgt Function file is not found!");
			std::cout << "Reading parameters from file "; fSidebandSyst->Print();

			double par1[3], par2[3];
			//double par1[2], par2[3];
			TF1 *tf1 = 0;
			TF1 *tf2 = 0;

			tf1 = dynamic_cast<TF1*> (fSidebandSyst->Get("pj1_Et_pho"));
			tf2 = dynamic_cast<TF1*> (fSidebandSyst->Get("pj1_Et_pho_plusSigma"));

			assert (tf1 != NULL && "sideband pho et  reweight function not found in file!");
			assert (tf2 != NULL && "sideband pho et error  reweight function not found in file!");

			tf1->GetParameters(par1);
			tf1_sideband_phoet = new TF1("sidebandPhoEtRewgtFunc_fitFunction3_p0+p1sqrt(x)+p2x",fitFunction3,0,1000,3);
			//tf1_sideband_phoet = new TF1("sidebandPhoEtRewgtFunc_p0+p1x","[0]+[1]*x",0,1000);
			assert (tf1_sideband_phoet != NULL && "tf1_sideband_phoet reweight function is not reconstructed porperly!");
			tf1_sideband_phoet->SetParameters(par1);

			tf2->GetParameters(par2);
			//tf1_sideband_phoet_error = new TF1("sidebandPhoEtRewgtErrorFunc_fitFunction3_p0+p1sqrt(x)+p2x",fitFunction3,0,1000,3);
			tf1_sideband_phoet_error = new TF1("sidebandPhoEtRewgtErrorFunc_pol2","[0]+[1]*x+[2]*x*x",0,1000);
			assert (tf1_sideband_phoet_error != NULL && "tf1_sideband_phoet_error reweight function is not reconstructed porperly!");
			tf1_sideband_phoet_error->SetParameters(par2);

			std::cout << "Photon Et Reweight function parameters " << std::endl;
			tf1_sideband_phoet->Print();
			std::cout << green << "Photon Et Reweight Error function parameters " << std::endl;
			tf1_sideband_phoet_error->Print();

		} else if (GetReweightSideband() == 1 || GetReweightSideband() == 2) //photon Et and Jet Eta weights for sideband reweight
		{
			TFile *fSidebandSyst = new TFile("RootFiles/g30jets_SidebandPhoEtReweightFunctions.root");
			assert (!fSidebandSyst->IsZombie() && "Sideband Wgt Function file is not found!");

			TFile *fSidebandSyst2 = new TFile("RootFiles/g30jets_SidebandJetEtaReweightFunctions.root");
			assert (!fSidebandSyst2->IsZombie() && "Sideband jet eta Wgt Function file is not found!");


			std::cout << red <<"Using sideband weight function/s in file/s " << std::endl;;
			fSidebandSyst->Print();
			fSidebandSyst2->Print();
			std::cout << clearatt << std::endl;

			double par1[3], par2[3], par3[9];
			TF1 *tf1 = 0;
			TF1 *tf2 = 0;
			TF1 *tf3 = 0;

			tf1 = dynamic_cast<TF1*> (fSidebandSyst->Get("pj1_Et_pho"));
			tf2 = dynamic_cast<TF1*> (fSidebandSyst->Get("pj1_Et_pho_plusSigma"));
			tf3 = dynamic_cast<TF1*> (fSidebandSyst2->Get("tf1_leadjeteta_nominalweights")); //see note below
			hist_JetEtaWghts = dynamic_cast<TH1*> (fSidebandSyst2->Get("phodata_copy"));

			assert (tf1 != NULL && "sideband pho et  reweight function not found in file!");
			assert (tf2 != NULL && "sideband pho et error  reweight function not found in file!");
			assert (tf3 != NULL && "sideband jet eta reweight function not found in file!");
			assert (hist_JetEtaWghts != NULL && "sideband jet eta reweight hist not found in file!");

			tf1->GetParameters(par1);
			tf1_sideband_phoet = new TF1("sidebandPhoEtRewgtFunc_fitFunction3_p0+p1sqrt(x)+p2x",fitFunction3,0,1000,3);
			assert (tf1_sideband_phoet != NULL && "tf1_sideband_phoet reweight function is not reconstructed porperly!");
			tf1_sideband_phoet->SetParameters(par1);

			tf2->GetParameters(par2);
			tf1_sideband_phoet_error = new TF1("sidebandPhoEtRewgtErrorFunc_fitFunction3_p0+p1sqrt(x)+p2x",fitFunction3,0,1000,3);
			assert (tf1_sideband_phoet_error != NULL && "tf1_sideband_phoet_error reweight function is not reconstructed porperly!");
			tf1_sideband_phoet_error->SetParameters(par2);

			std::cout << magenta << "Photon Et Reweight function parameters " << std::endl;
			tf1_sideband_phoet->Print();
			std::cout << green << "Photon Et Reweight Error function parameters " << std::endl;
			tf1_sideband_phoet_error->Print();


			tf3->GetParameters(par3);
			tf1_sideband_jeteta = new TF1("sidebandJetEtaRewgtFunc","gaus(0)+gaus(3)+gaus(6)", -3.4, 3.4);
			assert (tf1_sideband_jeteta != NULL && "tf1_sideband_jeteta reweight function is not reconstructed porperly!");
			tf1_sideband_jeteta->SetParameters(par3);

			/*	double perr[tf1_sideband_jeteta->GetNpar()];
				perr[0] = 0.028;
				perr[1] = 0.14;
				perr[2] = 0.202;
				perr[3] = 0.053;
				perr[4] = 0.1055;
				perr[5] = 0.0;
				perr[6] = 0.0250;
				perr[7] = 0.076;
				perr[8] = 0.0638;
				tf1_sideband_jeteta->SetParErrors(perr);
				*/

			std::cout << cyan << "Jet eta reweight function parameters " << std::endl;
			tf1_sideband_jeteta->Print();


			std::cout << clearatt << std::endl;
			std::cout << " ### END SIDEBAND REWEIGHING INFO ########## " << std::endl;
		} else {
			std::cout << __FUNCTION__ << ": " << __LINE__ << " UNKOWN SIDEBAND Rewight option given!" << std::endl;
			assert (false);
		}
	}
	
	if (GetUseDATAMCNvtxWeights()>0)
	{
		vNvtxWeights.clear();

		TH1 *hist =0;
		if (GetDataset() == 2)   //weights for photon MC 
		{
			//new weights from g+jets with MET cleanup
			vNvtxWeights.push_back(0.552832);
			vNvtxWeights.push_back(1.0657);
			vNvtxWeights.push_back(1.86598);
			vNvtxWeights.push_back(2.90396);
			vNvtxWeights.push_back(4.02392);
			vNvtxWeights.push_back(6.31669);
			vNvtxWeights.push_back(15.4261);
			vNvtxWeights.push_back(30.6712);

		} else if (GetDataset() == 9)  //diphoton mc rewgt 
		{
			//new weights from g+jets with MET cleanup
			vNvtxWeights.push_back(1.02879);
			vNvtxWeights.push_back(1.00456);
			vNvtxWeights.push_back(0.989921);
			vNvtxWeights.push_back(0.944556);
			vNvtxWeights.push_back(0.892121);
			vNvtxWeights.push_back(1.01257);
			vNvtxWeights.push_back(1.63495);
			vNvtxWeights.push_back(3.53568);
			vNvtxWeights.push_back(7.1378);
			vNvtxWeights.push_back(8.13572);

		} else if (GetDataset() == 3) //Zee MC
		{
/*			//central+stat err vals
			vNvtxWeights.push_back(0.600929);
			vNvtxWeights.push_back(1.07984);
			vNvtxWeights.push_back(1.74758);
			vNvtxWeights.push_back(2.48331);
			vNvtxWeights.push_back(3.46219);
			vNvtxWeights.push_back(5.69621);
			vNvtxWeights.push_back(11.6082);
			vNvtxWeights.push_back(22.8686);
*/

			//z+jets with MET cleanup cuts 8-25-2010
			vNvtxWeights.push_back(0.611806);
			vNvtxWeights.push_back(1.0525);
			vNvtxWeights.push_back(1.58434);
			vNvtxWeights.push_back(2.72301);
			vNvtxWeights.push_back(3.81061);
			vNvtxWeights.push_back(3.89529);
			vNvtxWeights.push_back(3.24607);
			if (GetUseDATAMCNvtxWeights() == 2) //do the central+sigma weighting
			{
				vNvtxWeights.at(0) += 0.0485179;
				vNvtxWeights.at(1) += 0.0771926;
				vNvtxWeights.at(2) += 0.158618;
				vNvtxWeights.at(3) += 0.401018;
				vNvtxWeights.at(4) += 1.05942;
				vNvtxWeights.at(5) += 2.30448;
				vNvtxWeights.at(6) += 3.74824;
			}

			
		} else if (GetDataset() == 4) //Zmm MC
		{
			//I am using these z->ee weights as I get no events 
			//from z+jets with clean up cuts! 9-2-2010
			vNvtxWeights.push_back(0.296238);
			//central+stat err vals
			if (GetUseDATAMCNvtxWeights() == 2) //do the central+sigma weighting
			{
				vNvtxWeights.at(0) += 0.29625;
			}
			
			//z+jets with MET cleanup cuts 8-25-2010
			//no events
		} else if (GetDataset() == 5) //Ztt MC
		{
			//I am using these z->ee weights as I get no events 
			//from z+jets with clean up cuts! 9-2-2010

			vNvtxWeights.push_back(0.514841);
			vNvtxWeights.push_back(1.39528);
			vNvtxWeights.push_back(1.74787);
			vNvtxWeights.push_back(2.57032);
			vNvtxWeights.push_back(2.46959);

			//central+stat err vals
			if (GetUseDATAMCNvtxWeights() == 2) //do the central+sigma weighting
			{
				vNvtxWeights.at(0) += 0.043004;
				vNvtxWeights.at(1) += 0.17616;
				vNvtxWeights.at(2) += 0.30952;
				vNvtxWeights.at(3) += 0.85777;
				vNvtxWeights.at(4) += 1.42755;
			}
			

			//z+jets with MET cleanup cuts 8-25-2010
			//no events

		} else if (GetDataset() == 6) //Wen MC)
		{
			/*
			//08-12-2010 : this is from w+jets see elog:
			vNvtxWeights.push_back(0.692896);
			vNvtxWeights.push_back(1.04487);
			vNvtxWeights.push_back(1.52813);
			vNvtxWeights.push_back(2.22268);
			vNvtxWeights.push_back(3.94164);
			vNvtxWeights.push_back(10.6812);
			vNvtxWeights.push_back(32.7236);
			vNvtxWeights.push_back(185.434);
			vNvtxWeights.push_back(1.0);
			vNvtxWeights.push_back(48.0754);
			// use central + stat err weights for systematics
			vNvtxWeights.push_back(0.692896 + 0.00590148);
			  vNvtxWeights.push_back(1.04487 + 0.00946897);
			  vNvtxWeights.push_back(1.52813 + 0.0193913);
			  vNvtxWeights.push_back(2.22268 + 0.0455044);
			  vNvtxWeights.push_back(3.94164 + 0.139866);
			  vNvtxWeights.push_back(10.6812 + 0.738173);
			  vNvtxWeights.push_back(32.7236 + 5.04018);
			  vNvtxWeights.push_back(185.434 + 109.025);
			  vNvtxWeights.push_back(1.0 + 1.0); //100% error );
			  vNvtxWeights.push_back(48.0754 + 51.3947);
			  */

			//08-25-2010 : this is from w+jets with met clean up cuts(delPhi(jet15,MET)>0.4 see elog:

			vNvtxWeights.push_back(0.531763);
			vNvtxWeights.push_back(1.09089);
			vNvtxWeights.push_back(2.01056);
			vNvtxWeights.push_back(3.42637);
			vNvtxWeights.push_back(4.65062);
			vNvtxWeights.push_back(8.09207);
			vNvtxWeights.push_back(32.3683);

			if (GetUseDATAMCNvtxWeights() == 2) //do the central+sigma weighting
			{
				vNvtxWeights.at(0) += 0.0154345;
				vNvtxWeights.at(1) += 0.0312086;
				vNvtxWeights.at(2) += 0.0790176;
				vNvtxWeights.at(3) += 0.227053;
				vNvtxWeights.at(4) += 0.626729;
				vNvtxWeights.at(5) += 2.49727;
				vNvtxWeights.at(6) += 33.69;
			}

		} else if (GetDataset() == 7) //Wmn MC
		{
			//I am using these w->ll weights as I get no events 
			//from w+jets with clean up cuts! 9-2-2010
			vNvtxWeights.push_back(0.569584);
			vNvtxWeights.push_back(0.823644);
			vNvtxWeights.push_back(2.2337);
			vNvtxWeights.push_back(2.53752);
			//central+stat weights
			if (GetUseDATAMCNvtxWeights() == 2) //do the central+sigma weighting
			{
				vNvtxWeights.at(0) += 0.105787;
				vNvtxWeights.at(1) += 0.16155;
				vNvtxWeights.at(2) += 0.84429;
				vNvtxWeights.at(3) += 1.46509;
			}

			//08-25-2010 : this is from w+jets with met clean up cuts(delPhi(jet15,MET)>0.4 see elog:
			//no events were found in MC

		} else if (GetDataset() == 8) //Wtn MC
		{
			/*vNvtxWeights.push_back(0.51911);
			  vNvtxWeights.push_back(1.04513);
			  vNvtxWeights.push_back(1.83233);
			  vNvtxWeights.push_back(2.5778);
			  vNvtxWeights.push_back(2.84779);
			  vNvtxWeights.push_back(3.36053);
			  vNvtxWeights.push_back(7.3911);
			  vNvtxWeights.push_back(7.64155);
			  //central+stat weights
			vNvtxWeights.push_back(0.526618);
			vNvtxWeights.push_back(1.06369);
			vNvtxWeights.push_back(1.88243);
			vNvtxWeights.push_back(2.69731);
			vNvtxWeights.push_back(3.08479);
			vNvtxWeights.push_back(3.90827);
			vNvtxWeights.push_back(10.415);
			vNvtxWeights.push_back(13.0562);
			*/

			//08-25-2010 : this is from w+jets with met clean up cuts(delPhi(jet15,MET)>0.4 see elog:
			//got only 80 events in MC.
			vNvtxWeights.push_back(0.541767);
			vNvtxWeights.push_back(1.09912);
			vNvtxWeights.push_back(2.26181);
			vNvtxWeights.push_back(2.75954);
			vNvtxWeights.push_back(1.20096);

			if (GetUseDATAMCNvtxWeights() == 2)
			{
				vNvtxWeights.at(0) += 0.0848229;
				vNvtxWeights.at(1) += 0.221419;
				vNvtxWeights.at(2) += 0.802496;
				vNvtxWeights.at(3) += 1.59784;
				vNvtxWeights.at(4) += 0.85485;
			}


		} else if (GetDataset() == 10) //WW  MC
		{
			vNvtxWeights.push_back(0.771982);
			vNvtxWeights.push_back(1.03336);
			vNvtxWeights.push_back(1.21094);
			vNvtxWeights.push_back(1.65976);
			vNvtxWeights.push_back(2.35735);
			vNvtxWeights.push_back(3.61096);
			vNvtxWeights.push_back(31.1623);
			if (GetUseDATAMCNvtxWeights() == 2)
			{
				vNvtxWeights.at(0) += 0.0166032;
				vNvtxWeights.at(1) += 0.0262575;
				vNvtxWeights.at(2) += 0.0449111;
				vNvtxWeights.at(3) += 0.110284;
				vNvtxWeights.at(4) += 0.293543;
				vNvtxWeights.at(5) += 0.842994;
				vNvtxWeights.at(6) += 31.2263;
			}

		} else if (GetDataset() == 11) //WZ  MC
		{
			vNvtxWeights.push_back(0.75121);
			vNvtxWeights.push_back(1.01193);
			vNvtxWeights.push_back(1.31046);
			vNvtxWeights.push_back(1.85815);
			vNvtxWeights.push_back(2.71961);
			vNvtxWeights.push_back(7.15026);

			if (GetUseDATAMCNvtxWeights() == 2)
			{
				vNvtxWeights.at(0) += 0.0183239;
				vNvtxWeights.at(1) += 0.0293686;
				vNvtxWeights.at(2) += 0.0582636;
				vNvtxWeights.at(3) += 0.151025;
				vNvtxWeights.at(4) += 0.421811;
				vNvtxWeights.at(5) += 2.72017;
			}

		} else if (GetDataset() == 12) //ZZ MC
		{
			vNvtxWeights.push_back(0.835025);
			vNvtxWeights.push_back(1.00522);
			vNvtxWeights.push_back(1.15153);
			vNvtxWeights.push_back(1.32211);
			vNvtxWeights.push_back(3.33063);
			vNvtxWeights.push_back(2.69623);
			vNvtxWeights.push_back(14.2741);
			if (GetUseDATAMCNvtxWeights() == 2)
			{
				vNvtxWeights.at(0) += 0.0413295;
				vNvtxWeights.at(1) += 0.0535735;
				vNvtxWeights.at(2) += 0.0912061;
				vNvtxWeights.at(3) += 0.172617;
				vNvtxWeights.at(4) += 0.868697;
				vNvtxWeights.at(5) += 1.11147;
				vNvtxWeights.at(6) += 14.9709;
			}

		} else if (GetDataset() == 13) //TTbar MC
		{
			vNvtxWeights.push_back(0.664237);
			vNvtxWeights.push_back(1.12296);
			vNvtxWeights.push_back(1.51586);
			vNvtxWeights.push_back(2.43992);
			vNvtxWeights.push_back(3.70212);
			vNvtxWeights.push_back(6.33803);
			vNvtxWeights.push_back(8.63631);
			if (GetUseDATAMCNvtxWeights() == 2)
			{
				vNvtxWeights.at(0) += 0.0239256;
				vNvtxWeights.at(1) += 0.0541302;
				vNvtxWeights.at(2) += 0.114942;
				vNvtxWeights.at(3) += 0.362358;
				vNvtxWeights.at(4) += 1.07383;
				vNvtxWeights.at(5) += 3.66951;
				vNvtxWeights.at(6) += 8.65406;
			}

		}

		std::cout << green << "======= Nvtx weights ===========" << std::endl; 
		for (int i=0;i<vNvtxWeights.size();++i)
			std::cout << "Nvtx " << i+1 << " = " << vNvtxWeights.at(i) << std::endl;
		std::cout << "================================" << clearatt << std::endl; 
	}

	
	//cannot add folder without creating a dir first!
	histFile = new TFile(sHistFileName.c_str(),"RECREATE");
	topDir = histFile->mkdir("Hist","All Histograms");	//create top level dir.
	topDir->cd();
	gDirectory->pwd();

	
} // Init

//-------------------------------------------------------------------
void PhotonJets::CleanUp()
{
	vMcCount.clear();
	vMcUpCount.clear();
	vMcDownCount.clear();
	vMcSidebandEWKCount.clear();
	vHaloCount.clear();
	vCosmicCount.clear();
	vSignalCount.clear();
	vSidebandHaloCount.clear();
	vSidebandCosmicCount.clear();
	vSidebandSystHadEmCount.clear();
	vSidebandSystIsoCount.clear();
	vSidebandSystTrkPtCount.clear();
	vSidebandSystTrkIsoCount.clear();

	histFile->Write();
	histFile->Close();
} 


//------ main -------------------------------------------------------
void PhotonJets::Main(TChain *ch, int iRunEvents)
{
	TBenchmark timer;
	timer.Start("phoana_time");
	
	
	if (ch) {
		myChain = ch;
		//ch->Print();	
	} else {
		std::cout << "NULL chain!" << std::endl;
		assert(false);
	}
	if (ch->GetEntries()<1)
	{
		std::cout << "No entries found!" << std::endl;
		return;
	}

	gROOT->Reset();
	Init(myChain);
	std::cout << yellow << "******* Job Settings ***********" << std::endl;
	PrintJobSettings();
	std::cout << "********************************" << clearatt << std::endl;
	
	BookHistograms();
	
	unsigned int iENTRIES = (unsigned int) myChain->GetEntries();
	unsigned int iEvt2Process = 0; 	//_____________________________________________ events to be processed
	if (iRunEvents < 0 ) {
		iRunEvents = iENTRIES;			//______________________________________________ requested all entries
		iEvt2Process = iENTRIES;
	} else if (iRunEvents > 0 && iRunEvents <= iENTRIES) iEvt2Process = iRunEvents;
	else iEvt2Process = iENTRIES;		//______________________________________________ default

	std::cout << " ==> Entries found :: " << iENTRIES << std::endl;

	//do a check for the MC. this is a common mistake I make
	
	  readIn->GetEntry(1);
	  if (stuple.evt_McFlag != GetMcFlag()) 
	  {
		  std::cout << red << "MC Flag," << GetMcFlag()
			  << ", setting does not match what is in Stuple, " 
			  << stuple.evt_McFlag << ". pleas check." 
			  << " returning!" << clearatt << std::endl;
		  //std::cout << "Changing the MCFlag to Stuple setting McFlag = " << stuple.evt_McFlag << clearatt << std::endl;
		  //SetMcFlag(stuple.evt_McFlag); //for some reason this seems to overwrite counting arrays? 12-07-2009
		  //return; // simple return causes trouble when rerun. so drastic measure needed, even exit(1) would not work
		  assert (false);
	  }

	unsigned int iEvtProc = 0;			//______________________________________________ number of events processed
	unsigned int iEleEvts = 0, iPhoEvts = 0, iPhoEleEvts =0, iOther=0;
	unsigned int iFirst400 =0;
	std::cout << "Init done. Ready, Set , GO... " << std::endl; 
	
	timer.Start("looptimer");
	
	iCurrTree = -1;
	for ( ; iEvtProc < iEvt2Process; ++iEvtProc) {
		if (iCurrTree != myChain->GetTreeNumber())
		{
			std::cout << green << "Opening file " << myChain->GetFile()->GetName() << clearatt << std::endl;
			iCurrTree = myChain->GetTreeNumber();
		}
	  	readIn->GetEntry(iEvtProc);
		
		//if (! (stuple.evt_RunNumber == 190925 &&  stuple.evt_EventNumber == 753747) ) continue;
	 	//if ( ! (stuple.evt_RunNumber == 190863 && stuple.evt_EventNumber == 109590)) continue;
		//if (stuple.vtx_NClass12 !=1) continue;
		//stuple.Dump();
		//std::cout << cyan << "****************************************************************************" << clearatt << std::endl;
		//std::cout << "entry, evt, run,Mc " << iEvtProc << "\t"<< stuple.evt_RunNumber << "\t" << stuple.evt_EventNumber << "\t" << stuple.evt_McFlag <<  std::endl;
		//stuple.DumpJetBlock();
		//stuple.DumpPhotonBlock();
		//stuple.DumpElectronBlock();
		//std::cout << cyan << "****************************************************************************" << clearatt << std::endl;
		//stuple.Dump(1);
		//std::cout << cyan << "****************************************************************************" << clearatt << std::endl;

		if (iEvtProc > 0 && (iEvtProc%iProgressBy == 0)) {
			std::cout << iEvtProc << "\t events processed [" << (int)( iEvtProc/(double)iEvt2Process * 100)<< "%] ";
			timer.Show("looptimer");
			timer.Start("looptimer");
		}

		  //___________________________________________________________________________ drop the first 400pb-1
		if (stuple.evt_McFlag == 0) {  //______________________________________________ if data
			if (stuple.evt_RunNumber < 190851) {
				iFirst400++;
				continue;
			}
		}
	
		//for halo id eff only 3-13-2010
		if (stuple.vtx_NClass12 == 0) 
		{
			DoHaloRejHists(stuple);
			//continue;
		} //else continue;
		
		
		//std::cout << "pho_num=" << stuple.pho_num << std::endl;
		//std::cout << "pho_num status=" << myChain->GetBranchStatus("pho_num") << std::endl;
	  	if (stuple.pho_num > 0 && stuple.ele_num > 0) iPhoEleEvts++;
	  	else if (stuple.ele_num > 0 && stuple.pho_num < 1) iEleEvts++;
	  	else if (stuple.ele_num < 1 && stuple.pho_num > 0) iPhoEvts++;
		else iOther++;

		if (FailStripCuts(stuple)) continue;
		//if (ExoticEvent(stuple)) 
		//{
		//	stuple.Dump();
		//} else continue;
		
		CommonVars commVars;
		PhotonList centralPhos, emupPhos, emdownPhos;
		ElectronList centralEles, emupEles, emdownEles;	//do not comment these out. I need these to recreate the jet list with all unused EM objects
		JetList centralJets, jesupJets, jesdownJets;

		GetCommonVariables(stuple, commVars);
		CreatePhotonLists(stuple, centralPhos,0);
		CreateElectronLists(stuple, centralEles,0);
		RemoveDuplicateEMobjects(centralPhos,centralEles); 
		CreateJetLists(stuple, centralJets,0);
		//centralJets.Dump();
		//centralPhos.Dump();
		//centralEles.Dump();
		
		if (bMcFlag)		//__________________________ need EM/JES UP/DOWN collections
		{
			CreatePhotonLists(stuple, emupPhos,1);
			CreatePhotonLists(stuple, emdownPhos,-1);
			CreateElectronLists(stuple, emupEles,1);
			CreateElectronLists(stuple, emdownEles,-1);

			RemoveDuplicateEMobjects(emupPhos,emupEles); 
			RemoveDuplicateEMobjects(emdownPhos,emdownEles); 
			CreateJetLists(stuple, jesupJets,1);
			CreateJetLists(stuple, jesdownJets,-1);
		}

		//centralJets.DumpRawJets();
	//	if (FoundPartialJet(commVars,centralPhos,centralEles,centralJets)) continue;
	//	if (bMcFlag)
	//	{
	//		if (FoundPartialJet(commVars, emupPhos,emupEles,jesupJets) || 
	//				FoundPartialJet(commVars, emdownPhos,emdownEles,jesdownJets))
	//			continue;
	//	}
		//check if the EM objects and jets overlap
		//CheckObjOverlap(commVars, centralPhos,centralEles,centralJets);
		//CheckJetTowers(commVars, centralJets);

		DoMyStuff(commVars, centralPhos, emupPhos, emdownPhos, centralEles, emupEles,
						emdownEles, centralJets, jesupJets, jesdownJets);


	}
	
	std::cout << green << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << clearatt << std::endl;
	std::cout << green << ">>>>>>>>>>> Requested events processed. Doing final calculations.." << clearatt << std::endl;
	if (!bMcFlag)
	{
		GetCosmicEstAndErr();
		DoHaloRejCalc();
		std::cout << "CPR Weights mean: for 1 jet = " << Hsignal.p1j_Pho.CprWeight->GetMean() << std::endl;
		std::cout << "CPR Weights mean: for 2 jets= " << Hsignal.p2j_Pho.CprWeight->GetMean() << std::endl;
	}
	std::cout << green << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << clearatt << std::endl;
	
	//print out the job summary
	std::cout << magenta << "======== PhotonJets =======================" << std::endl;
	std::cout << "Entries found    = " << iENTRIES << std::endl;
	std::cout << "Entries request  = " << iRunEvents << " (" << iRunEvents/(double) iENTRIES * 100 << "%)"<< std::endl;
	std::cout << "Events Processed = " << iEvtProc << " (" << iEvtProc/(double) iENTRIES * 100 << "%)"<< std::endl;
	std::cout << "Ele./Pho Events  = " << iPhoEleEvts << std::endl;
	std::cout << "Ele Events       = " << iEleEvts << std::endl;
	std::cout << "Pho Events       = " << iPhoEvts << std::endl;
	std::cout << "Other Events     = " << iOther << std::endl;
	std::cout << "First400 Events  = " << iFirst400 << std::endl;
	std::cout << "sum Events       = " << iPhoEvts + iEleEvts + iPhoEleEvts + iOther+iFirst400 << std::endl;
	std::cout << red << "__________________________________________" << std::endl;

	PrintJobSettings();
	std::cout << "Rejected Events        = " << std::setw(6) << iRejEvts << " (" << (int) iRejEvts*100/iEvtProc << "%)" << std::endl;
	
	const int itextWidth = 25;
	if (bMcFlag) {
		std::cout << "-------------- CENTRAL --------------------- "<< std::endl;
		std::cout << "Pho Evts               = " << vMcCount[0] << std::endl;
		std::cout << "Pho+>=1Jet Evts        = " << vMcCount[1] << std::endl; 
		std::cout << "Pho+>=2Jet Evts        = " << vMcCount[2] << std::endl; 
		std::cout << "-------------- EM/JES UP ------------------- "<< std::endl;
		std::cout << "Pho Evts               = " << vMcUpCount[0] << std::endl;
		std::cout << "Pho+>=1Jet Evts        = " << vMcUpCount[1] << std::endl; 
		std::cout << "Pho+>=2Jet Evts        = " << vMcUpCount[2] << std::endl; 
		std::cout << "-------------- EM/JES DOWN ----------------- "<< std::endl;
		std::cout << "Pho Evts               = " << vMcDownCount[0] << std::endl;
		std::cout << "Pho+>=1Jet Evts        = " << vMcDownCount[1] << std::endl; 
		std::cout << "Pho+>=2Jet Evts        = " << vMcDownCount[2] << std::endl; 
		std::cout << "-------------- SIDEBAND -------------------- "<< std::endl;
		std::cout << "Pho Evts               = " << vMcSidebandEWKCount[0] << std::endl;
		std::cout << "Pho+>=1Jet Evts        = " << vMcSidebandEWKCount[1] << std::endl; 
		std::cout << "Pho+>=2Jet Evts        = " << vMcSidebandEWKCount[2] << std::endl; 

	} else {
		if (GetApplyMetCleanUpCuts())
		{
		std::cout << std::setw(itextWidth) << "Fail MET Cleanup, Halo = " << std::setw(10) << vHaloCount[3];
		if (vHaloCount[0]>0) std::cout << " [" <<  vHaloCount[3]*100/vHaloCount[0] << "% of halo photons]"<< std::endl;
		else std::cout << std::endl;
		std::cout << std::setw(itextWidth) << "Fail MET Cleanup,Cosmic= " << std::setw(10) << vCosmicCount[3];
		if (vCosmicCount[0]>0) std::cout << " [" <<  vCosmicCount[3]*100/vCosmicCount[0] << " % of cosmic photons]"<< std::endl;
		else std::cout << std::endl;
		std::cout << std::setw(itextWidth) << "Fail MET Cleanup, QCD  = " << std::setw(10) << vQCDCount[3];
		if (vQCDCount[0]>0) std::cout << " [" <<  vQCDCount[3]*100/vQCDCount[0] << " % of sideband photons]"<< std::endl;
		else std::cout << std::endl;
		std::cout << std::setw(itextWidth) << "Fail MET Cleanup, Pho  = " << std::setw(10) << vSignalCount[3];
		if (vSignalCount[0]>0) std::cout << " [" <<  vSignalCount[3]*100/vSignalCount[0] << " % of signal photons]"<< std::endl;
		else std::cout << std::endl;
		}

		std::cout << red << "-------------- QCD SYETEMATIC TIGHT HADEM -- "<< std::endl;
		std::cout << std::setw(itextWidth) << "QCD Evts               = " << vSidebandSystHadEmCount[0] << std::endl;
		std::cout << std::setw(itextWidth) << "QCD+>=1Jet Evts        = " << vSidebandSystHadEmCount[1] << std::endl; 
		std::cout << std::setw(itextWidth) << "QCD+>=2Jet Evts        = " << vSidebandSystHadEmCount[2] << std::endl; 
		std::cout << red << "-------------- QCD SYETEMATIC TIGHT ISO ---- "<< std::endl;
		std::cout << std::setw(itextWidth) << "QCD Evts               = " << vSidebandSystIsoCount[0] << std::endl;
		std::cout << std::setw(itextWidth) << "QCD+>=1Jet Evts        = " << vSidebandSystIsoCount[1] << std::endl; 
		std::cout << std::setw(itextWidth) << "QCD+>=2Jet Evts        = " << vSidebandSystIsoCount[2] << std::endl; 
		std::cout << red << "-------------- QCD SYETEMATIC TIGHT TRKISO - "<< std::endl;
		std::cout << std::setw(itextWidth) << "QCD Evts               = " << vSidebandSystTrkIsoCount[0] << std::endl;
		std::cout << std::setw(itextWidth) << "QCD+>=1Jet Evts        = " << vSidebandSystTrkIsoCount[1] << std::endl; 
		std::cout << std::setw(itextWidth) << "QCD+>=2Jet Evts        = " << vSidebandSystTrkIsoCount[2] << std::endl; 
		std::cout << red << "-------------- QCD SYETEMATIC TIGHT TRKPT -- "<< std::endl;
		std::cout << std::setw(itextWidth) << "QCD Evts               = " << vSidebandSystTrkPtCount[0] << std::endl;
		std::cout << std::setw(itextWidth) << "QCD+>=1Jet Evts        = " << vSidebandSystTrkPtCount[1] << std::endl; 
		std::cout << std::setw(itextWidth) << "QCD+>=2Jet Evts        = " << vSidebandSystTrkPtCount[2] << std::endl; 
		
		std::cout << std::setw(itextWidth) <<"-------------- BEAM HALO ----tight ----- sideband" << std::endl;
		std::cout << std::setw(itextWidth) << "Halo Evts              = " << std::setw(10) << vHaloCount[0] << std::setw(10) << vSidebandHaloCount.at(0) << std::endl;
		std::cout << std::setw(itextWidth) << "Halo+>=1Jet Evts       = " << std::setw(10) << vHaloCount[1] << std::setw(10) << vSidebandHaloCount.at(1) << std::endl; 
		std::cout << std::setw(itextWidth) << "Halo+>=2Jet Evts       = " << std::setw(10) << vHaloCount[2] << std::setw(10) << vSidebandHaloCount.at(2) << std::endl; 
		std::cout << green << "-------------- COSMIC ---------------------- "<< std::endl;
		std::cout << std::setw(itextWidth) << "Cosmic Evts            = " << std::setw(10) << vCosmicCount[0] << std::setw(10) << vSidebandCosmicCount.at(0) << std::endl;
		std::cout << std::setw(itextWidth) << "Cosmic+>=1Jet Evts     = " << std::setw(10) << vCosmicCount[1] << std::setw(10) << vSidebandCosmicCount.at(1) << std::endl; 
		std::cout << std::setw(itextWidth) << "Cosmic+>=2Jet Evts     = " << std::setw(10) << vCosmicCount[2] << std::setw(10) << vSidebandCosmicCount.at(2) << std::endl; 
		std::cout << red << "-------------- QCD ------------------------- "<< std::endl;
		std::cout << std::setw(itextWidth) << "QCD Evts               = " << vQCDCount[0] << std::endl;
		std::cout << std::setw(itextWidth) << "QCD+>=1Jet Evts        = " << vQCDCount[1] << std::endl; 
		std::cout << std::setw(itextWidth) << "QCD+>=2Jet Evts        = " << vQCDCount[2] << std::endl; 
		std::cout << blue << "-------------- SIGNAL----------------------- "<< std::endl;
		std::cout << std::setw(itextWidth) << "Pho Evts               = " << vSignalCount[0] << std::endl;
		std::cout << std::setw(itextWidth) << "Pho+>=1Jet Evts        = " << vSignalCount[1] << std::endl; 
		std::cout << std::setw(itextWidth) << "Pho+>=2Jet Evts        = " << vSignalCount[2] << std::endl; 
	}
	std::cout << cyan << "Data Written to        = " << GetHistFileName() << clearatt << std::endl;
	std::cout << "======= END PhotonJets ====================" << std::endl;

/*	std::cout << "Events in mutiple bins::" << std::endl;
	std::cout << "Run" << "\t" << "Evt" << "\t" << "Halo" << "\t" << "Cosmic" << "\t" << "Signal" << "\t" << "QCD" << std::endl; 
	for (int i=0; i < EVT.size(); i++)
	{
		if (EVT[i].Halo + EVT[i].Cosmic + EVT[i].Signal + EVT[i].QCD >= 2)
			std::cout << EVT[i].Run << "\t" << EVT[i].Evt << "\t" << EVT[i].Halo << "\t" << EVT[i].Cosmic << "\t" << EVT[i].Signal << "\t" << EVT[i].QCD << std::endl; 
	}
*/
	EVT.clear();
	
	CleanUp();
	timer.Show("phoana_time");
} // Main


//-------------------------------------------------------------------
void PhotonJets::BookHistograms()
{
	//for halo estimate
	haloPhiWedge_1j[0] = new TH1F("haloPhiWedge_1j_b4","Halo(5)+>=1 Jet: before cuts: #phi^{#gamma}",24,0,24);
	haloPhiWedge_1j[1] = new TH1F("haloPhiWedge_1j_a4","Halo(5)+>=1 Jet: after cuts: #phi^{#gamma}",24,0,24);
	haloEmTime_1j[0] = new TH1F("haloEmTime_1j_b4","Halo(5)+>=1 Jet: before cuts: #gamma EmTime",30,-10,10);
	haloEmTime_1j[1] = new TH1F("haloEmTime_1j_a4","Halo(5)+>=1 Jet: after cuts: #gamma EmTime",30,-10,10);
	haloPhiWedge_2j[0] = new TH1F("haloPhiWedge_2j_b4","Halo(5)+>=2 Jet: before cuts: #phi^{#gamma}",24,0,24);
	haloPhiWedge_2j[1] = new TH1F("haloPhiWedge_2j_a4","Halo(5)+>=2 Jet: after cuts: #phi^{#gamma}",24,0,24);
	haloEmTime_2j[0] = new TH1F("haloEmTime_2j_b4","Halo(5)+>=2 Jet: before cuts: #gamma EmTime",30,-10,10);
	haloEmTime_2j[1] = new TH1F("haloEmTime_2j_a4","Halo(5)+>=2 Jet: after cuts: #gamma EmTime",30,-10,10);
	haloNvtx[0] = new TH1F("haloNvtx_b4","Halo(5): before cuts: NVtx (class 12)",10,0,10);
	haloNvtx[1] = new TH1F("haloNvtx_a4","Halo(5): after cuts: NVtx (class 12)",10,0,10);
	

	phoEt_iso25 = new TH1F("phoEt_Iso25","E_{T}^{#gamma} (pho_iso25 trigger)",1200,0,1200);
	phoEt_50 = new TH1F("phoEt_50","E_{T}^{#gamma} (pho_50 trigger)",1200,0,1200);
	phoEt_70 = new TH1F("phoEt_70","E_{T}^{#gamma} (pho_70 trigger)",1200,0,1200);

	
	//myMet = new TH1F("myMet","my #slash{E}_{T}",1200,0,1200);

	/*Hsignal.z1j_zpt	= new TH1F("z1j_zpt","#gamma+>=1Jet :- p_{T}^{Z} ",200,-500,500);
	Hsignal.z1j_zet	= new TH1F("z1j_zet","#gamma+>=1Jet :- E_{T}^{Z} ",250,0,500);
	Hsignal.z1j_zmass	= new TH1F("z1j_zmass","#gamma+>=1Jet :- Z mass ",200,0,1000);
	Hsignal.z1j_zeta	= new TH1F("z1j_zeta","#gamma+>=1Jet :- Event #Eta^{Z} ",50,-5,5);
	Hsignal.z1j_zj_mass	= new TH1F("z1j_zj1mass","#gamma+>=1Jet :- Invariant Mass (Z, Lead Jet)",200,0,1000);
	
	Hsignal.z2j_zpt	= new TH1F("z2j_zpt","#gamma+>=2Jets :- p_{T}^{Z} ",200,-500,500);
	Hsignal.z2j_zet	= new TH1F("z2j_zet","#gamma+>=2Jets :- E_{T}^{Z} ",250,0,500);
	Hsignal.z2j_zmass	= new TH1F("z2j_zmass","#gamma+>=2Jets :- Z mass ",200,0,1000);
	Hsignal.z2j_zeta	= new TH1F("z2j_zeta","#gamma+>=2Jets :- Event #Eta^{Z} ",50,-5,5);
	Hsignal.z2j_z1j_mass	= new TH1F("z2j_zj1mass","#gamma+>=2Jets :- Invariant Mass (Z, Lead Jet)",200,0,1000);
	Hsignal.z2j_z2j_mass	= new TH1F("z2j_zj2mass","#gamma+>=2Jets :- Invariant Mass (Z, Second Lead Jet)",200,0,1000);
	Hsignal.z2j_zj1j2_mass	= new TH1F("z2j_zj1j2mass","#gamma+>=2Jets :- Invariant Mass (Z, Two Lead Jets)",200,0,1000);
*/
	
	


	Histograms HistoMan;
	std::string text2title1, text2title2;

	TDirectory *l2dir, *l3dir, *l4dir;

	if (bMcFlag) {			//plots for MC

		//CENTRAL
		text2title1 = "CENTRAL #gamma_{MC}+>=1Jet :";
		text2title2 = "CENTRAL #gamma_{MC}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("CENTRAL","Central values");
		l2dir->cd();

		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hmc.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hmc.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hmc.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hmc.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hmc.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hmc.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hmc.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hmc.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hmc.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hmc.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hmc.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hmc.p2j_TwoJets,text2title2.c_str());



		//EM AND JES UP
		text2title1 = "EM/JES UP #gamma_{MC}+>=1Jet :";
		text2title2 = "EM/JES UP #gamma_{MC}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("EMJESUP","EM and JES up values");
		l2dir->cd();

		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HmcUp.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HmcUp.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcUp.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(HmcUp.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HmcUp.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HmcUp.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcUp.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcUp.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HmcUp.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HmcUp.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(HmcUp.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(HmcUp.p2j_TwoJets,text2title2.c_str());



		//EM AND JES DOWN
		text2title1 = "EM/JES DOWN #gamma_{MC}+>=1Jet :";
		text2title2 = "EM/JES DOWN #gamma_{MC}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("EMJESDOWN","EM and JES down values");
		l2dir->cd();

		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HmcDown.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HmcDown.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcDown.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(HmcDown.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HmcDown.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HmcDown.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcDown.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HmcDown.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HmcDown.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HmcDown.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(HmcDown.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(HmcDown.p2j_TwoJets,text2title2.c_str());


		//SIDEBAND MC EVENTS to be subtracted from DATA SIDEBAND TEMPLATE
		text2title1 = "SIDEBAND #gamma_{MC}+>=1Jet :";
		text2title2 = "SIDEBAND #gamma_{MC}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("SIDEBAND","SIDEBAND values");
		l2dir->cd();

		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsideband_ewk.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsideband_ewk.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsideband_ewk.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hsideband_ewk.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsideband_ewk.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsideband_ewk.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsideband_ewk.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsideband_ewk.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsideband_ewk.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsideband_ewk.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hsideband_ewk.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hsideband_ewk.p2j_TwoJets,text2title2.c_str());



	} else {  //DATA PLOTS

			//plots for signal
		text2title1 = "#gamma_{signal}+>=1Jet :";
		text2title2 = "#gamma_{signal}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("SIGNAL","Signal");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsignal.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsignal.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsignal.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hsignal.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsignal.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsignal.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsignal.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsignal.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsignal.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsignal.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hsignal.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hsignal.p2j_TwoJets,text2title2.c_str());



			//plots for halo
		text2title1 = "#gamma_{halo}+>=1Jet :";
		text2title2 = "#gamma_{halo}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("HALO","Beam Halo");
		l2dir->cd();
	
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hhalo.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hhalo.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hhalo.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hhalo.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hhalo.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hhalo.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hhalo.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hhalo.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hhalo.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hhalo.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hhalo.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hhalo.p2j_TwoJets,text2title2.c_str());



			//plots for cosmic
		text2title1 = "#gamma_{cosmic}+>=1Jet :";
		text2title2 = "#gamma_{cosmic}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("COSMIC","Cosmic");
		l2dir->cd();

		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hcosmic.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hcosmic.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hcosmic.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hcosmic.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hcosmic.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hcosmic.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hcosmic.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hcosmic.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hcosmic.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hcosmic.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hcosmic.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hcosmic.p2j_TwoJets,text2title2.c_str());



			//qcd sideband
		text2title1 = "#gamma_{sideband}+>=1Jet :";
		text2title2 = "#gamma_{sideband}+>=2Jet :";
		
		topDir->cd();

		l2dir = gDirectory->mkdir("SIDEBAND","Photon sideband to model QCD");
		l2dir->cd();
		
		hSidebandJetEtaWeights = new TProfile("sidebandJetEtaWgts","Sideband Weights based on leading jet #eta",800, -4,4);
		hSidebandPhoEtWeights = new TProfile("sidebandPhoEtWgts","Sideband Weights based on photon Et",1000, 0,500);
		
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hqcd.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hqcd.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hqcd.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hqcd.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hqcd.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hqcd.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hqcd.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hqcd.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hqcd.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hqcd.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hqcd.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hqcd.p2j_TwoJets,text2title2.c_str());


		//sideband cosmic
		text2title1 = "#gamma_{sideband cosmic}+>=1Jet :";
		text2title2 = "#gamma_{sideband cosmic}+>=2Jet :";
		
		topDir->cd();

		l2dir = gDirectory->mkdir("SIDEBAND_COSMIC","COSMICS in QCD sideband");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsideband_cosmic.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsideband_cosmic.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsideband_cosmic.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hsideband_cosmic.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsideband_cosmic.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsideband_cosmic.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsideband_cosmic.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsideband_cosmic.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsideband_cosmic.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsideband_cosmic.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hsideband_cosmic.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hsideband_cosmic.p2j_TwoJets,text2title2.c_str());


				//sideband halo
		text2title1 = "#gamma_{sideband halo}+>=1Jet :";
		text2title2 = "#gamma_{sideband halo}+>=2Jet :";
		
		topDir->cd();

		l2dir = gDirectory->mkdir("SIDEBAND_HALO","Sideband Halo hists");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsideband_halo.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsideband_halo.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsideband_halo.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(Hsideband_halo.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(Hsideband_halo.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(Hsideband_halo.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsideband_halo.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(Hsideband_halo.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsideband_halo.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(Hsideband_halo.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(Hsideband_halo.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(Hsideband_halo.p2j_TwoJets,text2title2.c_str());

			
				//plost for beam halo rejection power
		topDir->cd();

		l2dir = gDirectory->mkdir("HALOREJPOW","Plots for Halo rejection power");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("SIGNALHALO_B4","For halo in tight photon sample");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HhaloNoVtx_Evt_b4, "No Vtx beam halo (signal) before BH cuts");

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HhaloNoVtx_Pho_b4, "No Vtx beam halo (signal) before BH cuts");
		
		l2dir->cd();
		l3dir = gDirectory->mkdir("SIGNALHALO_A4","For halo in tight photon sample");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HhaloNoVtx_Evt_a4, "No Vtx beam halo (signal) after BH cuts:");

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HhaloNoVtx_Pho_a4, "No vtx beam halo (signal) after BH cuts");

		l2dir->cd();
		l3dir = gDirectory->mkdir("SIDEBANDHALO_B4","For halo in sideband photon sample");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HhalosidebandNoVtx_Evt_b4, "No Vtx beam halo (sideband) before BH cuts");

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HhalosidebandNoVtx_Pho_b4, "No Vtx beam halo (sideband) before BH cuts");
		
		l2dir->cd();
		l3dir = gDirectory->mkdir("SIDEBANDHALO_A4","For halo in tight photon sample");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HhalosidebandNoVtx_Evt_a4, "No Vtx beam halo (sideband) after BH cuts:");

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HhalosidebandNoVtx_Pho_a4, "No vtx beam halo (sideband) after BH cuts");



		//plos for sideband systematics
		// TIGHT HADEM
		text2title1 = "SIDEBAND SYSTEMATIC: Tight HadEm : #gamma_{sideband}+>=1Jet :";
		text2title2 = "SIDEBAND SYSTEMATIC: Tight HadEm : #gamma_{sideband}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("SIDEBANDSYST_HADEM","Photon Sideband Systematic with Tight HadEm cut");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HsidebandSyst_HadEm.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HsidebandSyst_HadEm.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_HadEm.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_HadEm.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HsidebandSyst_HadEm.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HsidebandSyst_HadEm.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_HadEm.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_HadEm.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_HadEm.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_HadEm.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(HsidebandSyst_HadEm.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(HsidebandSyst_HadEm.p2j_TwoJets,text2title2.c_str());



		// TIGHT ISO
		text2title1 = "SIDEBAND SYSTEMATIC: Tight Iso : #gamma_{sideband}+>=1Jet :";
		text2title2 = "SIDEBAND SYSTEMATIC: Tight Iso : #gamma_{sideband}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("SIDEBANDSYST_ISO","Photon Sideband Systematic with Tight Iso cut");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HsidebandSyst_Iso.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HsidebandSyst_Iso.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_Iso.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_Iso.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HsidebandSyst_Iso.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HsidebandSyst_Iso.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_Iso.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_Iso.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_Iso.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_Iso.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(HsidebandSyst_Iso.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(HsidebandSyst_Iso.p2j_TwoJets,text2title2.c_str());


		// TIGHT TRKPT
		text2title1 = "SIDEBAND SYSTEMATIC: Tight TrkPt : #gamma_{sideband}+>=1Jet :";
		text2title2 = "SIDEBAND SYSTEMATIC: Tight TrkPt : #gamma_{sideband}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("SIDEBANDSYST_TRKPT","Photon Sideband Systematic with Tight TrkPt cut");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HsidebandSyst_TrkPt.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HsidebandSyst_TrkPt.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_TrkPt.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_TrkPt.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HsidebandSyst_TrkPt.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HsidebandSyst_TrkPt.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_TrkPt.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_TrkPt.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_TrkPt.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_TrkPt.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(HsidebandSyst_TrkPt.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(HsidebandSyst_TrkPt.p2j_TwoJets,text2title2.c_str());

		// TIGHT TRKISO
		text2title1 = "SIDEBAND SYSTEMATIC: Tight TrkIso : #gamma_{sideband}+>=1Jet :";
		text2title2 = "SIDEBAND SYSTEMATIC: Tight TrkIso : #gamma_{sideband}+>=2Jet :";
		
		topDir->cd();
		l2dir = gDirectory->mkdir("SIDEBANDSYST_TRKISO","Photon Sideband Systematic with Tight TrkIso cut");
		l2dir->cd();
		
		l3dir = gDirectory->mkdir("1Jet","1 Jet case");
		l3dir->cd();
		
		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HsidebandSyst_TrkIso.p1j_Evt,text2title1.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HsidebandSyst_TrkIso.p1j_Pho,text2title1.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_TrkIso.p1j_Jet,text2title1.c_str());
				
		l3dir->cd();		
		l4dir = gDirectory->mkdir("PhotonLeadJet","Photon and Lead Jet Histograms");
		l4dir->cd();
		gDirectory->Print();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_TrkIso.p1j_PhoJet,text2title1.c_str());

	
		l2dir->cd();
		l3dir = gDirectory->mkdir("2Jet","#gamma + >=2 Jet case");
		l3dir->cd();

		l4dir = gDirectory->mkdir("Event","Event Histograms");
		l4dir->cd();
				HistoMan.GetEventHistograms(HsidebandSyst_TrkIso.p2j_Evt,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon","Photon Histograms");
		l4dir->cd();
				HistoMan.GetPhotonHistograms(HsidebandSyst_TrkIso.p2j_Pho,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("LeadJet","Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_TrkIso.p2j_Jet1,text2title2.c_str());

		l3dir->cd();
		l4dir = gDirectory->mkdir("SecondLeadJet","2^{nd} Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetJetHistograms(HsidebandSyst_TrkIso.p2j_Jet2,text2title2.c_str());
	
		l3dir->cd();
		l4dir = gDirectory->mkdir("PhotonLeadJet","#gamma and Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_TrkIso.p2j_PhoJet1,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2ndLeadJet","#gamma and 2nd Lead Jet Histograms");
		l4dir->cd();
				HistoMan.GetPhoton1JetHistograms(HsidebandSyst_TrkIso.p2j_PhoJet2,text2title2.c_str());
					
		l3dir->cd();
		l4dir = gDirectory->mkdir("Photon2Jets","#gamma + 2Jets Histograms");
		l4dir->cd();
				HistoMan.GetPhoton2JetsHistograms(HsidebandSyst_TrkIso.p2j_PhoJets,text2title2.c_str());
				
		l3dir->cd();
		l4dir = gDirectory->mkdir("2Jets","2Jets Histograms");
		l4dir->cd();
				HistoMan.GetTwoJetsHistograms(HsidebandSyst_TrkIso.p2j_TwoJets,text2title2.c_str());

				
		
	}
	
} //BookHistograms



void PhotonJets::CreatePhotonLists(const Stuple stuple, PhotonList& pholist, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
			pholist.pho_num 	 =  stuple.pho_num;
			pholist.pho_Ntight = stuple.pho_Ntight;
			pholist.pho_Nloose = stuple.pho_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_E[i]);
				pholist.pho_Px.push_back(stuple.pho_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_matchJetIndex[i]);
				
				pholist.pho_CprWgt.push_back(stuple.pho_CprWgt[i]);
				pholist.pho_CprSys1.push_back(stuple.pho_CprSys1[i]);
				pholist.pho_CprSys2.push_back(stuple.pho_CprSys2[i]);
				pholist.pho_CprSys3.push_back(stuple.pho_CprSys3[i]);
				pholist.pho_CprSys4.push_back(stuple.pho_CprSys4[i]);
				pholist.pho_CprSys5.push_back(stuple.pho_CprSys5[i]);
				pholist.pho_CprSys6.push_back(stuple.pho_CprSys6[i]);
				pholist.pho_CprSys7.push_back(stuple.pho_CprSys7[i]);
				pholist.pho_CprSys8.push_back(stuple.pho_CprSys8[i]);
				
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		case 1:
		{
			pholist.pho_num 	 =  stuple.pho_up_num;
			pholist.pho_Ntight = stuple.pho_up_Ntight;
			pholist.pho_Nloose = stuple.pho_up_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_up_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_up_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_up_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_up_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_up_E[i]);
				pholist.pho_Px.push_back(stuple.pho_up_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_up_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_up_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_up_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_up_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_up_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_up_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_up_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_up_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_up_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_up_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_up_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_up_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_up_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_up_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_up_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_up_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_up_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_up_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_up_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_up_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_up_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_up_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_up_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_up_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_up_matchJetIndex[i]);
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		case -1:
		{
			pholist.pho_num 	 =  stuple.pho_down_num;
			pholist.pho_Ntight = stuple.pho_down_Ntight;
			pholist.pho_Nloose = stuple.pho_down_Nloose;
			
			for (unsigned int i=0; i < stuple.pho_down_num; i++)
			{
				pholist.pho_Index.push_back(stuple.pho_down_Index[i]);
				pholist.pho_PhoBlockIndex.push_back(stuple.pho_down_PhoBlockIndex[i]);
				pholist.pho_Etc.push_back(stuple.pho_down_Etc[i]);
				pholist.pho_E.push_back(stuple.pho_down_E[i]);
				pholist.pho_Px.push_back(stuple.pho_down_Px[i]);
				pholist.pho_Py.push_back(stuple.pho_down_Py[i]);
				pholist.pho_Pz.push_back(stuple.pho_down_Pz[i]);
				pholist.pho_Detector.push_back(stuple.pho_down_Detector[i]);
				pholist.pho_DetEta.push_back(stuple.pho_down_DetEta[i]);
				pholist.pho_DetPhi.push_back(stuple.pho_down_DetPhi[i]);
				pholist.pho_XCes.push_back(stuple.pho_down_XCes[i]);
				pholist.pho_ZCes.push_back(stuple.pho_down_ZCes[i]);
				pholist.pho_HadEm.push_back(stuple.pho_down_HadEm[i]);
				pholist.pho_Chi2Mean.push_back(stuple.pho_down_Chi2Mean[i]);
				pholist.pho_N3d.push_back(stuple.pho_down_N3d[i]);
				pholist.pho_Iso4.push_back(stuple.pho_down_Iso4[i]);
				pholist.pho_TrkPt.push_back(stuple.pho_down_TrkPt[i]);	
				pholist.pho_TrkIso.push_back(stuple.pho_down_TrkIso[i]);
				pholist.pho_CesWireE2.push_back(stuple.pho_down_CesWireE2[i]);
				pholist.pho_CesStripE2.push_back(stuple.pho_down_CesStripE2[i]);
				pholist.pho_PhiWedge.push_back(stuple.pho_down_PhiWedge[i]);
				pholist.pho_NMuonStubs.push_back(stuple.pho_down_NMuonStubs[i]);
				pholist.pho_EmTime.push_back(stuple.pho_down_EmTime[i]);
				pholist.pho_TightId.push_back(stuple.pho_down_TightId[i]);
				pholist.pho_LooseId.push_back(stuple.pho_down_LooseId[i]);
				pholist.pho_PhoenixId.push_back(stuple.pho_down_PhoenixId[i]);
				pholist.pho_Halo_seedWedge.push_back(stuple.pho_down_Halo_seedWedge[i]);
				pholist.pho_Halo_eastNhad.push_back(stuple.pho_down_Halo_eastNhad[i]);
				pholist.pho_Halo_westNhad.push_back(stuple.pho_down_Halo_westNhad[i]);
				pholist.pho_matchJetIndex.push_back(stuple.pho_down_matchJetIndex[i]);
			}

			assert(pholist.pho_num == pholist.pho_Index.size());
			break;
		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid photon collection requested! can make central/em up/em down photon lists only!");
			assert(false);

		}
		
	}	//switch

/*}}}*/
}	//CreatePhotonLists



void PhotonJets::CreateElectronLists(const Stuple stuple, ElectronList& list, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
			list.ele_num    = stuple.ele_num;
			list.ele_Ntight = stuple.ele_Ntight;
			list.ele_Nloose = stuple.ele_Nloose;

			for (unsigned int i=0; i < stuple.ele_num; i++)
			{
				list.ele_Index.push_back(stuple.ele_Index[i]);
				list.ele_PhoBlockIndex.push_back(stuple.ele_PhoBlockIndex[i]);
				list.ele_EleBlockIndex.push_back(stuple.ele_EleBlockIndex[i]);
				list.ele_Etc.push_back(stuple.ele_Etc[i]);
				list.ele_E.push_back(stuple.ele_E[i]);
				list.ele_Px.push_back(stuple.ele_Px[i]);
				list.ele_Py.push_back(stuple.ele_Py[i]);
				list.ele_Pz.push_back(stuple.ele_Pz[i]);
				list.ele_Detector.push_back(stuple.ele_Detector[i]);
				list.ele_DetEta.push_back(stuple.ele_DetEta[i]);
				list.ele_DetPhi.push_back(stuple.ele_DetPhi[i]);
				list.ele_XCes.push_back(stuple.ele_XCes[i]);
				list.ele_ZCes.push_back(stuple.ele_ZCes[i]);
				list.ele_HadEm.push_back(stuple.ele_HadEm[i]);
				list.ele_Chi2Mean.push_back(stuple.ele_Chi2Mean[i]);
				list.ele_N3d.push_back(stuple.ele_N3d[i]);
				list.ele_Iso4.push_back(stuple.ele_Iso4[i]);
				list.ele_TrkIso.push_back(stuple.ele_TrkIso[i]);		
				list.ele_CesWireE2.push_back(stuple.ele_CesWireE2[i]);
				list.ele_CesStripE2.push_back(stuple.ele_CesStripE2[i]);
				list.ele_PhiWedge.push_back(stuple.ele_PhiWedge[i]);
				list.ele_NMuonStubs.push_back(stuple.ele_NMuonStubs[i]);
				list.ele_EmTime.push_back(stuple.ele_EmTime[i]);
				list.ele_PhoenixId.push_back(stuple.ele_PhoenixId[i]);
				list.ele_Halo_seedWedge.push_back(stuple.ele_Halo_seedWedge[i]);
				list.ele_Halo_eastNhad.push_back(stuple.ele_Halo_eastNhad[i]);
				list.ele_Halo_westNhad.push_back(stuple.ele_Halo_westNhad[i]);
				list.ele_matchJetIndex.push_back(stuple.ele_matchJetIndex[i]);
				list.ele_Ntracks.push_back(stuple.ele_Ntracks[i]);	
				list.ele_Emfr.push_back(stuple.ele_Emfr[i]);
				list.ele_EoverP.push_back(stuple.ele_EoverP[i]);
				list.ele_TrackPt.push_back(stuple.ele_TrackPt[i]);
				list.ele_TrackBcPt.push_back(stuple.ele_TrackBcPt[i]);
				list.ele_TrackPhi.push_back(stuple.ele_TrackPhi[i]);
				list.ele_Nssl.push_back(stuple.ele_Nssl[i]);
				list.ele_Nasl.push_back(stuple.ele_Nasl[i]);
				list.ele_TightId.push_back(stuple.ele_TightId[i]);
				list.ele_LooseId.push_back(stuple.ele_LooseId[i]);
				list.ele_ConversionId.push_back(stuple.ele_ConversionId[i]);
			}

			assert( (int)list.ele_num == (int)list.ele_Index.size());
			break;
		}

		case 1:
		{
			list.ele_num    = stuple.ele_up_num;
			list.ele_Ntight = stuple.ele_up_Ntight;
			list.ele_Nloose = stuple.ele_up_Nloose;

			for (unsigned int i=0; i < stuple.ele_up_num; i++)
			{
				list.ele_Index.push_back(stuple.ele_up_Index[i]);
				list.ele_PhoBlockIndex.push_back(stuple.ele_up_PhoBlockIndex[i]);
				list.ele_EleBlockIndex.push_back(stuple.ele_up_EleBlockIndex[i]);
				list.ele_Etc.push_back(stuple.ele_up_Etc[i]);
				list.ele_E.push_back(stuple.ele_up_E[i]);
				list.ele_Px.push_back(stuple.ele_up_Px[i]);
				list.ele_Py.push_back(stuple.ele_up_Py[i]);
				list.ele_Pz.push_back(stuple.ele_up_Pz[i]);
				list.ele_Detector.push_back(stuple.ele_up_Detector[i]);
				list.ele_DetEta.push_back(stuple.ele_up_DetEta[i]);
				list.ele_DetPhi.push_back(stuple.ele_up_DetPhi[i]);
				list.ele_XCes.push_back(stuple.ele_up_XCes[i]);
				list.ele_ZCes.push_back(stuple.ele_up_ZCes[i]);
				list.ele_HadEm.push_back(stuple.ele_up_HadEm[i]);
				list.ele_Chi2Mean.push_back(stuple.ele_up_Chi2Mean[i]);
				list.ele_N3d.push_back(stuple.ele_up_N3d[i]);
				list.ele_Iso4.push_back(stuple.ele_up_Iso4[i]);
				list.ele_TrkIso.push_back(stuple.ele_up_TrkIso[i]);		
				list.ele_CesWireE2.push_back(stuple.ele_up_CesWireE2[i]);
				list.ele_CesStripE2.push_back(stuple.ele_up_CesStripE2[i]);
				list.ele_PhiWedge.push_back(stuple.ele_up_PhiWedge[i]);
				list.ele_NMuonStubs.push_back(stuple.ele_up_NMuonStubs[i]);
				list.ele_EmTime.push_back(stuple.ele_up_EmTime[i]);
				list.ele_PhoenixId.push_back(stuple.ele_up_PhoenixId[i]);
				list.ele_Halo_seedWedge.push_back(stuple.ele_up_Halo_seedWedge[i]);
				list.ele_Halo_eastNhad.push_back(stuple.ele_up_Halo_eastNhad[i]);
				list.ele_Halo_westNhad.push_back(stuple.ele_up_Halo_westNhad[i]);
				list.ele_matchJetIndex.push_back(stuple.ele_up_matchJetIndex[i]);
				list.ele_Ntracks.push_back(stuple.ele_up_Ntracks[i]);	
				list.ele_Emfr.push_back(stuple.ele_up_Emfr[i]);
				list.ele_EoverP.push_back(stuple.ele_up_EoverP[i]);
				list.ele_TrackPt.push_back(stuple.ele_up_TrackPt[i]);
				list.ele_TrackBcPt.push_back(stuple.ele_up_TrackBcPt[i]);
				list.ele_TrackPhi.push_back(stuple.ele_up_TrackPhi[i]);
				list.ele_Nssl.push_back(stuple.ele_up_Nssl[i]);
				list.ele_Nasl.push_back(stuple.ele_up_Nasl[i]);
				list.ele_TightId.push_back(stuple.ele_up_TightId[i]);
				list.ele_LooseId.push_back(stuple.ele_up_LooseId[i]);
				list.ele_ConversionId.push_back(stuple.ele_up_ConversionId[i]);
			}

			assert(list.ele_num == list.ele_Index.size());
			break;

		}

		case -1:
		{
			list.ele_num    = stuple.ele_down_num;
			list.ele_Ntight = stuple.ele_down_Ntight;
			list.ele_Nloose = stuple.ele_down_Nloose;

			for (unsigned int i=0; i < stuple.ele_down_num; i++)
			{
				list.ele_Index.push_back(stuple.ele_down_Index[i]);
				list.ele_PhoBlockIndex.push_back(stuple.ele_down_PhoBlockIndex[i]);
				list.ele_EleBlockIndex.push_back(stuple.ele_down_EleBlockIndex[i]);
				list.ele_Etc.push_back(stuple.ele_down_Etc[i]);
				list.ele_E.push_back(stuple.ele_down_E[i]);
				list.ele_Px.push_back(stuple.ele_down_Px[i]);
				list.ele_Py.push_back(stuple.ele_down_Py[i]);
				list.ele_Pz.push_back(stuple.ele_down_Pz[i]);
				list.ele_Detector.push_back(stuple.ele_down_Detector[i]);
				list.ele_DetEta.push_back(stuple.ele_down_DetEta[i]);
				list.ele_DetPhi.push_back(stuple.ele_down_DetPhi[i]);
				list.ele_XCes.push_back(stuple.ele_down_XCes[i]);
				list.ele_ZCes.push_back(stuple.ele_down_ZCes[i]);
				list.ele_HadEm.push_back(stuple.ele_down_HadEm[i]);
				list.ele_Chi2Mean.push_back(stuple.ele_down_Chi2Mean[i]);
				list.ele_N3d.push_back(stuple.ele_down_N3d[i]);
				list.ele_Iso4.push_back(stuple.ele_down_Iso4[i]);
				list.ele_TrkIso.push_back(stuple.ele_down_TrkIso[i]);		
				list.ele_CesWireE2.push_back(stuple.ele_down_CesWireE2[i]);
				list.ele_CesStripE2.push_back(stuple.ele_down_CesStripE2[i]);
				list.ele_PhiWedge.push_back(stuple.ele_down_PhiWedge[i]);
				list.ele_NMuonStubs.push_back(stuple.ele_down_NMuonStubs[i]);
				list.ele_EmTime.push_back(stuple.ele_down_EmTime[i]);
				list.ele_PhoenixId.push_back(stuple.ele_down_PhoenixId[i]);
				list.ele_Halo_seedWedge.push_back(stuple.ele_down_Halo_seedWedge[i]);
				list.ele_Halo_eastNhad.push_back(stuple.ele_down_Halo_eastNhad[i]);
				list.ele_Halo_westNhad.push_back(stuple.ele_down_Halo_westNhad[i]);
				list.ele_matchJetIndex.push_back(stuple.ele_down_matchJetIndex[i]);
				list.ele_Ntracks.push_back(stuple.ele_down_Ntracks[i]);	
				list.ele_Emfr.push_back(stuple.ele_down_Emfr[i]);
				list.ele_EoverP.push_back(stuple.ele_down_EoverP[i]);
				list.ele_TrackPt.push_back(stuple.ele_down_TrackPt[i]);
				list.ele_TrackBcPt.push_back(stuple.ele_down_TrackBcPt[i]);
				list.ele_TrackPhi.push_back(stuple.ele_down_TrackPhi[i]);
				list.ele_Nssl.push_back(stuple.ele_down_Nssl[i]);
				list.ele_Nasl.push_back(stuple.ele_down_Nasl[i]);
				list.ele_TightId.push_back(stuple.ele_down_TightId[i]);
				list.ele_LooseId.push_back(stuple.ele_down_LooseId[i]);
				list.ele_ConversionId.push_back(stuple.ele_down_ConversionId[i]);
			}

			assert(list.ele_num == list.ele_Index.size());
			break;
		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid collection requested! can make central/em up/em down  lists only!");
			assert(false);
		}
	}
/*}}}*/
}	//CreateElectronLists

void PhotonJets::CreateJetLists(const Stuple stuple, JetList& list, const int collection)
{
/*{{{*/
	switch (collection)
	{
		case 0:
		{
		  list.jet_num  = stuple.jet_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_NJet15;
			for (unsigned int i=0; i < stuple.jet_num; i++)
			{
				if ( fabs(stuple.jet_Px[i]) > 1500 || 
					  fabs(stuple.jet_Py[i]) > 1500 || 
					  fabs(stuple.jet_Pz[i]) > 1500 || 
					  fabs(stuple.jet_E[i]) > 1500)
				{
					stuple.DumpJetBlock();
					assert(false);
				}

				list.jet_Index.push_back(stuple.jet_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_E[i]);
				list.jet_Px.push_back(stuple.jet_Px[i]);
				list.jet_Py.push_back(stuple.jet_Py[i]);
				list.jet_Pz.push_back(stuple.jet_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_SeedIEta[i]);

				//if (UseBTagJets()) // Old Stuples has no btag info
			//	{
					list.jet_EmTime.push_back(stuple.jet_EmTime[i]);
			//	}
			}

			assert(list.jet_num == list.jet_Index.size());

			list.jet_raw_num  = stuple.jet_raw_num;							// number of raw jets
			for (unsigned int i=0; i < stuple.jet_raw_num; i++)
			{
				list.jet_raw_Index.push_back(stuple.jet_raw_Index[i]);				//original index in jet block
				list.jet_raw_Pt.push_back(stuple.jet_raw_Pt[i]);
				list.jet_raw_E.push_back(stuple.jet_raw_E[i]);
				list.jet_raw_Px.push_back(stuple.jet_raw_Px[i]);
				list.jet_raw_Py.push_back(stuple.jet_raw_Py[i]);
				list.jet_raw_Pz.push_back(stuple.jet_raw_Pz[i]);
			}
			
			break;
		}

		case 1:
		{
		  list.jet_num  = stuple.jet_up_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_up_NJet15;
			for (unsigned int i=0; i < stuple.jet_up_num; i++)
			{
				list.jet_Index.push_back(stuple.jet_up_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_up_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_up_E[i]);
				list.jet_Px.push_back(stuple.jet_up_Px[i]);
				list.jet_Py.push_back(stuple.jet_up_Py[i]);
				list.jet_Pz.push_back(stuple.jet_up_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_up_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_up_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_up_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_up_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_up_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_up_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_up_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_up_SeedIEta[i]);
				
				//if (UseBTagJets()) // Old Stuples has no btag info
				{
					list.jet_EmTime.push_back(stuple.jet_up_EmTime[i]);
				}
			}

			assert(list.jet_num == list.jet_Index.size());
			break;

			
		}

		case -1:
		{
		  list.jet_num  = stuple.jet_down_num;							// number of jets so we know how much to read from arrays
		  list.jet_NJet15 = stuple.jet_down_NJet15;
			for (unsigned int i=0; i < stuple.jet_down_num; i++)
			{
				list.jet_Index.push_back(stuple.jet_down_Index[i]);				//original index in jet block
				list.jet_Pt.push_back(stuple.jet_down_Pt[i]);					//easy to get from JetFilter
				list.jet_E.push_back(stuple.jet_down_E[i]);
				list.jet_Px.push_back(stuple.jet_down_Px[i]);
				list.jet_Py.push_back(stuple.jet_down_Py[i]);
				list.jet_Pz.push_back(stuple.jet_down_Pz[i]);
				list.jet_DetEta.push_back(stuple.jet_down_DetEta[i]);
				list.jet_DetPhi.push_back(stuple.jet_down_DetPhi[i]);
				list.jet_HadEm.push_back(stuple.jet_down_HadEm[i]);
				list.jet_Emfr.push_back(stuple.jet_down_Emfr[i]);
				list.jet_Ntowers.push_back(stuple.jet_down_Ntowers[i]);
				list.jet_Ntracks.push_back(stuple.jet_down_Ntracks[i]);
				list.jet_SeedIPhi.push_back(stuple.jet_down_SeedIPhi[i]);
				list.jet_SeedIEta.push_back(stuple.jet_down_SeedIEta[i]);
				
				//if (UseBTagJets()) // Old Stuples has no btag info
				{
					list.jet_EmTime.push_back(stuple.jet_down_EmTime[i]);
				}
			}

			assert(list.jet_num == list.jet_Index.size());
			break;

		}

		default:
		{
			StdOut(__FILE__,__LINE__,3,"Invalid collection requested! can make central and JES up/down lists only!");
			assert(false);

		}

	}
/*}}}*/
} // CreateJetLists

void PhotonJets::GetCommonVariables(const Stuple stuple, CommonVars& cVars)
{
/*{{{*/
	cVars.evt_McFlag 		= stuple.evt_McFlag;
	cVars.evt_RunNumber 	= stuple.evt_RunNumber;
	cVars.evt_EventNumber= stuple.evt_EventNumber;
	cVars.tri_pho25iso 	= stuple.tri_pho25iso;
	cVars.tri_pho50 		= stuple.tri_pho50;
	cVars.tri_pho70 		= stuple.tri_pho70;

	cVars.vtx_N 			= stuple.vtx_N;
	cVars.vtx_NClass12 	= stuple.vtx_NClass12;
	cVars.vtx_z 			= stuple.vtx_z;
	cVars.vtx_Ntracks 	= stuple.vtx_Ntracks;
	cVars.vtx_SumPt 		= stuple.vtx_SumPt;

	cVars.met_Met 			= stuple.met_Met;	
	cVars.met_MetX 		= stuple.met_MetX;	
	cVars.met_MetY 		= stuple.met_MetY;	
	cVars.met_SumEt 		= stuple.met_SumEt;
	cVars.met_Ht 			= stuple.met_Ht;	
	cVars.met_MetPhi 		= stuple.met_MetPhi;
//	cVars.met_Gen_d 		= stuple.met_Gen_d;	
//	cVars.met_Gen_m 		= stuple.met_Gen_m;
//	cVars.met_Gen_p 		= stuple.met_Gen_p;
//	cVars.met_Gen_mUn 	= stuple.met_Gen_mUn;
//	cVars.met_Gen_pUn 	= stuple.met_Gen_pUn;
/*}}}*/
}


void PhotonJets::DoMyStuff(CommonVars cVars, PhotonList cPhos, PhotonList uPhos, PhotonList dPhos,
									ElectronList cEles, ElectronList uEles, ElectronList dEles,
									JetList cJets, JetList uJets, JetList dJets)
{	
	//if (cPhos.pho_num + cEles.ele_num <2) return;
	//std::cout << "___begin DOMyStuff " << std::endl;
		//cPhos.Dump();
		//cJets.Dump();

	//std::cout << "_____________________ " << std::endl;
	
	//don't apply this cut here! BH require no vtx events!
	//_________________________________________________________________ require a one good vertex 		  
	//if (cVars.vtx_NClass12 <= 0) return;
	//if (cVars.vtx_z > 60.)  return;

	if (!bMcFlag)   		//_________________ 0=data, 1=mc
	{
		// DATA path
		
		//______________________________________________________________ trigger
		if (  (cVars.tri_pho25iso != 1) && (cVars.tri_pho50 != 1) 	//_______________ only data has L1-3 trig. info
				&& (cVars.tri_pho70 != 1) ) 										//_______________ mc has L1-2.
			return;												
			
		//require at least one tight photon
		//remove e+pho events
		//if (stuple.pho_Ntight < 1) return;			//this will screw up the qcd template. so drop this
		//if (cEles.ele_Ntight > 0) return;
		
		std::vector<int> UsedPho, UsedEle, QcdPho, QcdEle;
	
		//______________________________________________________________ select the photon
		for (unsigned int i = 0; i < cPhos.pho_num; ++i )
		{
			if (cPhos.pho_PhoenixId[i] != 0) continue;  //____________________________ reject phoenix photons
			//________________________________________________________________________ this is a hack to find the trigger photon, as I don't have the L3 info to match reconstructed photon to trigger photon.
			if (cVars.tri_pho70 == 1 && cPhos.pho_Etc[i] < 70.0) continue;
			if (cVars.tri_pho50 == 1 && cPhos.pho_Etc[i] < 50.0) continue;
			if (cVars.tri_pho25iso == 1 && cPhos.pho_Etc[i] < 25.0) continue;

			if (fabs(cPhos.pho_DetEta.at(i))>fMaxPhotonEta) continue;
			//if (cPhos.pho_Etc[i]<fMinPhotonEt > cPhos.pho_Etc[i]>fMaxPhotonEt) continue;		//need this as I am tagging Et>7GeV em objects. 
			if (cPhos.pho_TightId[i] == 0)
			{
				if (UsedPho.size() == 0) UsedPho.push_back(i);	//________ need only one. else problem
			} 
			if (cPhos.pho_TightId[i] > 0 && cPhos.pho_LooseId[i] == 0)			// >0 means that this has failed tight cuts. must do this as I init them to large negative values and != could mean <0!!
			{
				if (QcdPho.size() == 0) QcdPho.push_back(i);	//________ need only one. else problem
			}
		}	
	
		//if (UsedPho.size() > 0) 
		//	std::cout << "picked pho ind = " << UsedPho[0] << std::endl;
		//if (QcdPho.size() > 0) 
		//	std::cout << "picked qcd ind = " << QcdPho[0] << std::endl;
	
		// need copy for QCD stuff. because the PHOTON may be different and hence a different jet list
		// need this because QcdPho is different. If I try to add unused jet list,
		// this will cause error! as some variables are not set correctly for many reasons!
		// for. eg. when I add an electron to jet list, jet_num is increased but ele_num is not
		// decreased and the electron is not removed from the array.!
		PhotonList qcdcPhos = cPhos;
		//PhotonList qcduPhos = uPhos;
		//PhotonList qcddPhos = dPhos;
		ElectronList qcdcEles = cEles;
		//ElectronList qcduEles = uEles;
		//ElectronList qcddEles = dEles;
		JetList qcdcJets = cJets;
		//JetList qcduJets = uJets;
		//JetList qcddJets = dJets;
		
		int used = 0;				//____ to make sure evts are not fallen into multiple bins
	
		EvtTag_t evt;
		evt.Run 		= cVars.evt_RunNumber;
		evt.Evt 		= cVars.evt_EventNumber;
		evt.Halo 	= 0;
		evt.Cosmic 	= 0;
		evt.Signal 	= 0;
		evt.QCD 		= 0;
		evt.SidebandCosmic = 0;

	
		if (UsedPho.size() > 0 && cPhos.pho_Etc[UsedPho.at(0)]>fMinPhotonEt			//need this as I am tagging Et>7GeV em objects. Need to update rest of this code
						&& cPhos.pho_Etc[UsedPho.at(0)]<fMaxPhotonEt)
																								// when looking for electrons. 11-19-2009
		{
			//cJets.Dump();
			//________________________________________________________________ now put back the un-used em objs to the jet list and reorder the jets
			NewJetList nj;		// this need to be worked out to fit the new changes and new jet lists
			nj.AddUnused(cPhos, cEles, cJets, UsedPho, UsedEle, fMinJetEt, fMaxJetEta);
			
			// ______________________________________________________________ USE em timing of THE photon to select cosmic
			if (cPhos.pho_EmTime[UsedPho[0]] > fMinCosmicEmTime &&
				 cPhos.pho_EmTime[UsedPho[0]] < fMaxCosmicEmTime)
			{
				//_________________________________________________________________ require a one good vertex 		  
				if (cVars.vtx_NClass12 >= iMinClass12Vtx 
						&& cVars.vtx_NClass12 <= iMaxClass12Vtx 
						&& fabs(cVars.vtx_z) < fMaxVtxz)
				{
					used = JetSelAndHistFill(cVars, cPhos, cEles, cJets, &Hcosmic, UsedPho, 
							UsedEle, vCosmicCount, false);
					if (used) evt.Cosmic = 1;
					if (used && PrintLevel()>iINFO) { std::cout << __LINE__ << ":: COSMIC EVENT "; cVars.PrintHeader();}
				}
			}

			if (!used)
			{
				if (cPhos.pho_EmTime[UsedPho[0]] > fMinSignalEmTime &&
					 cPhos.pho_EmTime[UsedPho[0]] < fMaxSignalEmTime) 
				{
					//_________________________________________ pick events from no vtx for halo 
					if (cVars.vtx_NClass12 == 0)
					{
						//these two lines are for halo rejection power in signal
						if (cVars.met_Ht > fMinPhotonEt+fMinJetEt) //temp hack to some crazy events! throw them away
						{
							FillPhotonHists(cVars, cPhos, &(HhaloNoVtx_Pho_b4), UsedPho[0]);
							FillEventHists(cVars, cPhos, cEles, cJets, &(HhaloNoVtx_Evt_b4), UsedPho[0]);
						}
						//_______________________________________________________________ select for Halo template
						HaloTypes ht(cPhos.pho_Halo_seedWedge[UsedPho[0]], 
							cPhos.pho_Halo_eastNhad[UsedPho[0]] + cPhos.pho_Halo_westNhad[UsedPho[0]]);

						if (ht.IsType(iHaloType))
						{
							//these two lines are for halo rejection power in signal
							if (cVars.met_Ht > fMinPhotonEt+fMinJetEt) //temp hack to some crazy events! throw them away
							{
								FillPhotonHists(cVars, cPhos, &(HhaloNoVtx_Pho_a4), UsedPho[0]);
								FillEventHists(cVars, cPhos, cEles, cJets, &(HhaloNoVtx_Evt_a4), UsedPho[0]);
							}

							if (cPhos.pho_PhiWedge[0] == 0 || cPhos.pho_PhiWedge[0] == 23)
							{
								used = JetSelAndHistFill(cVars, cPhos, cEles, cJets,
														&Hhalo, UsedPho, UsedEle, vHaloCount, false); 
								if (used) evt.Halo = 1;
								if (used && PrintLevel()>iINFO) { std::cout << __LINE__ << ":: HALO EVENT "; cVars.PrintHeader(); }
							}
						}
					}

					
					//________________________________________________________________ what ever left is my signal
					if (!used)
					{
						//_________________________________________________________________ require a one good vertex 		  
						if (cVars.vtx_NClass12 >= iMinClass12Vtx 
							&& cVars.vtx_NClass12 <= iMaxClass12Vtx
							&& fabs(cVars.vtx_z) < fMaxVtxz
				 			&& cPhos.pho_NMuonStubs[UsedPho[0]] == 0)						//________ there is 6% in-eff here. but we know ~60% of cosmics have muton stubs
						{
							if (PrintLevel()>iDEBUG)
							{
								std::cout << cyanbg << "SIGNAL JETS" << clearatt << std::endl;
								cVars.Dump();
								cPhos.Dump();
								cEles.Dump();
								cJets.DumpAll();
							}
							used = JetSelAndHistFill(cVars, cPhos, cEles, cJets,
													&Hsignal, UsedPho, UsedEle, vSignalCount, false);
							//do not set 'used'. this event could go into QCD bin and it need to know if this event is not
							// halo or cosmic!
							if (used) evt.Signal = 1;
							if (used && PrintLevel()>iINFO) { std::cout << __LINE__ << ":: SIGNAL EVENT "; cVars.PrintHeader();}
							//if (used) { std::cout << __LINE__ << ":: SIGNAL EVENT "; cVars.PrintHeader();}

							if (cVars.tri_pho70 == 1) phoEt_70->Fill(cPhos.pho_Etc.at(UsedPho.at(0))); 
							else if (cVars.tri_pho50 == 1) phoEt_50->Fill(cPhos.pho_Etc.at(UsedPho.at(0))); 
							else if (cVars.tri_pho25iso == 1) phoEt_iso25->Fill(cPhos.pho_Etc.at(UsedPho.at(0))); 
						}
					}
				}
			}
	
		} // if UsedPho.size() > 0
		
			
		// if this event is NOT halo or NOT cosmic:
		// now where we have a tight photon and a sideband
		// photon, this event must go into both sideband and signal bin.
		// if we put it only in the sideband bin, it will be filled with
		// events without any tight photons! which is not the reality!
		//_______________________________________________________________ select for QCD template
		if (!used)		// _______________________________________________ if this event is already in the cosmic or a halo bin
		{														
			if (QcdPho.size() > 0) 						// ____________________ we must have sideband photon
			{
				if (qcdcPhos.pho_EmTime[QcdPho[0]] > fMinSignalEmTime &&
					 qcdcPhos.pho_EmTime[QcdPho[0]] < fMaxSignalEmTime) 
				{
					//_________________________________________________________________ require a one good vertex 		  
					if (cVars.vtx_NClass12 >= iMinClass12Vtx 
							&& cVars.vtx_NClass12 <= iMaxClass12Vtx 
							&& fabs(cVars.vtx_z) < fMaxVtxz)
					{

						//reject beam halo in the sideband
						HaloTypes ht(qcdcPhos.pho_Halo_seedWedge[QcdPho[0]], 
							qcdcPhos.pho_Halo_eastNhad[QcdPho[0]] + qcdcPhos.pho_Halo_westNhad[QcdPho[0]]);

						if (ht.IsType(iHaloType)) return;
						
			
						//trying to add back the isolation energy taken out back to this jet
						//this should improve the results in elog#1526
						const int phoInd = QcdPho.at(0);
						if (GetAddIsoEtoSidebandPho())
						{
							const float iso = qcdcPhos.pho_Iso4.at(phoInd);
							
							//but Etc has to be set manually
							TLorentzVector tlPho(qcdcPhos.pho_Px.at(phoInd), qcdcPhos.pho_Py.at(phoInd), 
									qcdcPhos.pho_Pz.at(phoInd), qcdcPhos.pho_E.at(phoInd));
							const float pt = qcdcPhos.pho_Etc.at(phoInd) = tlPho.Pt();
							const float scale = (pt + iso)/pt;
							tlPho *= scale;
													
							qcdcPhos.pho_E.at(phoInd) = tlPho.E();  	//this will take care of the correct Pt derivation 
							qcdcPhos.pho_Px.at(phoInd) = tlPho.Px();  // using TLorentzVector during hist filling
							qcdcPhos.pho_Py.at(phoInd) = tlPho.Py();  
							qcdcPhos.pho_Pz.at(phoInd) = tlPho.Pz();  
							qcdcPhos.pho_Etc.at(phoInd) = tlPho.Pt();

							//now update the Ht too
							cVars.met_Ht += iso;

						}
						if (qcdcPhos.pho_Etc.at(phoInd)<fMinPhotonEt || qcdcPhos.pho_Etc.at(phoInd)>fMaxPhotonEt) return;
						
						const float fPhoEt = qcdcPhos.pho_Etc.at(QcdPho.at(0));
						float fWgt = 1;
						//if (GetReweightSideband()) fWgt = tf1_sideband_phoet->Eval(fPhoEt);

						//________________________________________________________________ now put back the un-used em objs to the jet list and reorder the jets
						NewJetList nj2;		// this need to be worked out to fit the new changes and new jet lists
						nj2.AddUnused(qcdcPhos, qcdcEles, qcdcJets, QcdPho, QcdEle, fMinJetEt, fMaxJetEta);

						
						used = JetSelAndHistFill(cVars, qcdcPhos, qcdcEles, qcdcJets,
													&Hqcd, QcdPho, QcdEle, vQCDCount, false, fWgt, false , GetReweightSideband());
						if (used) evt.QCD = 1;
						if (used && PrintLevel()>iINFO) { std::cout << __LINE__ << ":: SIDEBAND EVENT "; cVars.PrintHeader(); }

						
						//now do the systematic error for the choice
						//assuming the jet list remained unchanged.
						//when reweightin, the variations should be reweithted and comapred to reweigthed-central value -- 05-24-2010

						if (SidebandPhoPassTightHadEmCut(qcdcPhos.pho_E.at(QcdPho.at(0)), qcdcPhos.pho_HadEm.at(QcdPho.at(0))))
						{
							JetSelAndHistFill(cVars, qcdcPhos, qcdcEles, qcdcJets,
									&HsidebandSyst_HadEm, QcdPho, QcdEle, vSidebandSystHadEmCount, false, 1, false ,GetReweightSideband());
						}
						if (SidebandPhoPassTightIsoCut(qcdcPhos.pho_Etc.at(QcdPho.at(0)), qcdcPhos.pho_Iso4.at(QcdPho.at(0))))
						{
							JetSelAndHistFill(cVars, qcdcPhos, qcdcEles, qcdcJets,
									&HsidebandSyst_Iso, QcdPho, QcdEle, vSidebandSystIsoCount, false, 1, false ,GetReweightSideband());
						}
						if (SidebandPhoPassTightTrkPtCut(qcdcPhos.pho_Etc.at(QcdPho.at(0)), qcdcPhos.pho_TrkPt.at(QcdPho.at(0))))
						{
							JetSelAndHistFill(cVars, qcdcPhos, qcdcEles, qcdcJets,
									&HsidebandSyst_TrkPt, QcdPho, QcdEle, vSidebandSystTrkPtCount, false, 1, false ,GetReweightSideband());
						}
						if (SidebandPhoPassTightTrkIsoCut(qcdcPhos.pho_Etc.at(QcdPho.at(0)), qcdcPhos.pho_TrkIso.at(QcdPho.at(0))))
						{
							JetSelAndHistFill(cVars, qcdcPhos, qcdcEles, qcdcJets,
									&HsidebandSyst_TrkIso, QcdPho, QcdEle, vSidebandSystTrkIsoCount, false, 1, false ,GetReweightSideband());
						}

						
					} else if (cVars.vtx_NClass12 == 0) // sideband beam halo
					{

						if (qcdcPhos.pho_Etc.at(QcdPho[0])<fMinPhotonEt || qcdcPhos.pho_Etc.at(QcdPho[0]) > fMaxPhotonEt) return;
						//_______________________________________________________________ select for Halo template
						HaloTypes ht(qcdcPhos.pho_Halo_seedWedge[QcdPho[0]], 
							qcdcPhos.pho_Halo_eastNhad[QcdPho[0]] + qcdcPhos.pho_Halo_westNhad[QcdPho[0]]);

						if (ht.IsType(iHaloType))
						{
							if (qcdcPhos.pho_PhiWedge[QcdPho[0]] == 0 || qcdcPhos.pho_PhiWedge[QcdPho[0]] == 23)
							{
								used = JetSelAndHistFill(cVars, qcdcPhos, cEles, cJets,
														&Hsideband_halo, QcdPho, QcdEle, vSidebandHaloCount, false, 1, false, GetReweightSideband()); 
								if (used) evt.SidebandHalo = 1;
							}
						}
					}

				
				} else if (qcdcPhos.pho_EmTime[QcdPho.at(0)] > fMinCosmicEmTime &&  //cosmics in the sideband
				 				qcdcPhos.pho_EmTime[QcdPho.at(0)] < fMaxCosmicEmTime)
				{
					if (qcdcPhos.pho_Etc.at(QcdPho.at(0))<fMinPhotonEt || qcdcPhos.pho_Etc.at(QcdPho.at(0))>fMaxPhotonEt) return;
					//_________________________________________________________________ require a one good vertex 		  
					if (cVars.vtx_NClass12 >= iMinClass12Vtx 
							&& cVars.vtx_NClass12<= iMaxClass12Vtx 
							&& fabs(cVars.vtx_z) < fMaxVtxz)
					{
						int used = JetSelAndHistFill(cVars, qcdcPhos, qcdcEles, qcdcJets,
													&Hsideband_cosmic, QcdPho, QcdEle, vSidebandCosmicCount, false, 1, false, GetReweightSideband());
						if (used) evt.SidebandCosmic = 1;
					}
				}

				
			}
		}

		// now store the summary of this event
		EVT.push_back(evt);
		
	} else {

		// I had this at beginnig of this method. this causes NO BH events as they are
		// choosen from nvtx==0 events! So had to apply these cuts here.
		//12-02-2009
		//_________________________________________________________________ require a one good vertex 		  
		if (cVars.vtx_NClass12 < iMinClass12Vtx || cVars.vtx_NClass12> iMaxClass12Vtx) return;
		else if (fabs(cVars.vtx_z) > fMaxVtxz)  return;

		// MC path
		std::vector<int> UsedPhoc, UsedElec, UsedSidebandPho, UsedSidebandEle;
		std::vector<int> UsedPhou, UsedEleu;
		std::vector<int> UsedPhod, UsedEled;

		
		// CENTRAL
		//first select the photon, then add un-used em obj to jet list, then fill the hists
		SelectPhoton(cVars, cPhos, cEles, cJets, UsedPhoc, UsedElec);
		
		if (UsedPhoc.size() > 0)
		{
			const int phoInd = UsedPhoc.at(0);
			if (GetAddIsoEtoSidebandPho())
			{
				const float iso = cPhos.pho_Iso4.at(phoInd);

				//but Etc has to be set manually
				TLorentzVector tlPho(cPhos.pho_Px.at(phoInd), cPhos.pho_Py.at(phoInd), 
						cPhos.pho_Pz.at(phoInd), cPhos.pho_E.at(phoInd));
				const float pt = cPhos.pho_Etc.at(phoInd) = tlPho.Pt();
				const float scale = (pt + iso)/pt;
				//std::cout << "tlPho M= " << tlPho.E();
				tlPho *= scale;
				//std::cout << " -> " << tlPho.E()  << "(iso, scale = " << iso << ", " << scale<< std::endl;

				cPhos.pho_E.at(phoInd) = tlPho.E();  	//this will take care of the correct Pt derivation 
				cPhos.pho_Px.at(phoInd) = tlPho.Px();  // using TLorentzVector during hist filling
				cPhos.pho_Py.at(phoInd) = tlPho.Py();  
				cPhos.pho_Pz.at(phoInd) = tlPho.Pz();  
				cPhos.pho_Etc.at(phoInd) = tlPho.Pt();
			}
			//do not return, let the JES/EM UP DOWN methods run to see if they have a photon with min Et! 03-27-2010
			//if (cPhos.pho_Etc.at(phoInd)<fMinPhotonEt || cPhos.pho_Etc.at(phoInd)>fMaxPhotonEt ) return;
			
			if (cPhos.pho_Etc.at(phoInd)>=fMinPhotonEt 
				&& cPhos.pho_Etc.at(phoInd)<=fMaxPhotonEt )
			{
				//________________________________________________________________ now put back the un-used em objs to the jet list and reorder the jets
				NewJetList nj;
				nj.AddUnused(cPhos, cEles, cJets, UsedPhoc, UsedElec, fMinJetEt, fMaxJetEta);
				//________________________________________________________________ Select Jets and fill hists
				JetSelAndHistFill(cVars, cPhos, cEles, cJets, &Hmc, UsedPhoc, UsedElec, vMcCount, false);
			}

		} else {
			
			// SIDEBAND // central and sideband are made exclusive. see SelectPhoton() routine 01-02-2010
			//first select the photon, then add un-used em obj to jet list, then fill the hists
			
			SelectMCSidebandPhoton(cVars, cPhos, cEles, cJets, UsedSidebandPho, UsedSidebandEle);
			if (UsedSidebandPho.size()>0)
			{
				const int phoInd = UsedSidebandPho.at(0);
				if (cPhos.pho_Etc.at(phoInd)>=fMinPhotonEt 
						&& cPhos.pho_Etc.at(phoInd)<=fMaxPhotonEt )
				{
					//________________________________________________________________ now put back the un-used em objs to the jet list and reorder the jets
					NewJetList nj;
					nj.AddUnused(cPhos, cEles, cJets, UsedSidebandPho, UsedSidebandEle, fMinJetEt, fMaxJetEta);
					//________________________________________________________________ Select Jets and fill hists
					JetSelAndHistFill(cVars, cPhos, cEles, cJets, &Hsideband_ewk, UsedSidebandPho, UsedSidebandEle, vMcSidebandEWKCount, false);
				}
			}
			
		}


		// JES AND EM UP
		//first select the photon, then add un-used em obj to jet list, then fill the hists
		SelectPhoton(cVars, uPhos, uEles, uJets, UsedPhou, UsedEleu);
		if (UsedPhou.size() > 0)
		{
			const int phoInd = UsedPhou.at(0);
			if (uPhos.pho_Etc.at(phoInd)>=fMinPhotonEt 
					&& uPhos.pho_Etc.at(phoInd)<=fMaxPhotonEt )
			{
				//________________________________________________________________ now put back the un-used em objs to the jet list and reorder the jets
				NewJetList nj;
				nj.AddUnused(uPhos, uEles, uJets, UsedPhou, UsedEleu, fMinJetEt, fMaxJetEta);

				//________________________________________________________________ Select Jets and fill hists
				JetSelAndHistFill(cVars, uPhos, uEles, uJets, &HmcUp, UsedPhou, UsedEleu, vMcUpCount, false);
			}
		}

		// JES AND EM DOWN
		//first select the photon, then add un-used em obj to jet list, then fill the hists
		SelectPhoton(cVars, dPhos, dEles, dJets, UsedPhod, UsedEled);
		if (UsedPhod.size() > 0)
		{
			const int phoInd = UsedPhod.at(0);
			if (dPhos.pho_Etc.at(phoInd)>=fMinPhotonEt 
					&& dPhos.pho_Etc.at(phoInd)<=fMaxPhotonEt )
			{
				//________________________________________________________________ now put back the un-used em objs to the jet list and reorder the jets
				NewJetList nj;
				nj.AddUnused(dPhos, dEles, dJets, UsedPhod, UsedEled, fMinJetEt, fMaxJetEta);
				//________________________________________________________________ Select Jets and fill hists
				JetSelAndHistFill(cVars, dPhos, dEles, dJets, &HmcDown, UsedPhod, 
						UsedEled, vMcDownCount, false);
			}
		}


	}	// if data/mc
		

}


void PhotonJets::SelectPhoton(const CommonVars Vars, const PhotonList phos,
									const ElectronList eles, const JetList jets, 
									std::vector<int>& UsedPho, std::vector<int>& UsedEle)
{

		//after general evt selection (run, trigger, vtx) everything else should be same as signal selection
		//require at least one tight photon
		//remove e+pho events
		//if (phos.pho_Ntight < 1) return; 
	
		UsedPho.clear();
		UsedEle.clear();
		
		
		//______________________________________________________________ select the photon
		for (unsigned int i = 0; i < phos.pho_num; ++i ) {
			if (phos.pho_PhoenixId[i] != 0) continue;  //____________________________ reject phoenix photons
			if (fabs(phos.pho_DetEta.at(i)) > fMaxPhotonEta) continue;
			if (phos.pho_TightId[i] == 0 && UsedPho.size() == 0) UsedPho.push_back(i); // _____ tight photons
			//break;		//______________________________ one is enough. and this is the highest tight Et photon in the event as I have sorted the SuperPhoton list.
		}	
	
		if (UsedPho.size() > 0) {					//______ any photons were found
			//_______________________________________________________________ reject if for Halo template
			HaloTypes ht(phos.pho_Halo_seedWedge[UsedPho[0]], 
					phos.pho_Halo_eastNhad[UsedPho[0]] + phos.pho_Halo_westNhad[UsedPho[0]]);

			if (ht.IsType(iHaloType)) {
				UsedPho.clear();
			}
		}
}

void PhotonJets::SelectMCSidebandPhoton(const CommonVars Vars, const PhotonList phos,
									const ElectronList eles, const JetList jets, 
									std::vector<int>& UsedSidebandPho, std::vector<int>& UsedSidebandEle)
{

		//after general evt selection (run, trigger, vtx) everything else should be same as signal selection
		//require at least one tight photon
		//remove e+pho events
	
		UsedSidebandPho.clear();
		UsedSidebandEle.clear();
	
		//______________________________________________________________ select the photon
		for (unsigned int i = 0; i < phos.pho_num; ++i ) {
			if (phos.pho_PhoenixId[i] != 0) continue;  //____________________________ reject phoenix photons
			if (phos.pho_TightId[i] != 0 && phos.pho_LooseId[i] == 0
					&& UsedSidebandPho.size() == 0) UsedSidebandPho.push_back(i); // _____ sideband photons
			//break;		//______________________________ one is enough. and this is the highest tight Et photon in the event as I have sorted the SuperPhoton list.
		}	
	
		if (UsedSidebandPho.size() > 0) {					//______ any photons were found
			//_______________________________________________________________ reject if for Halo template
			HaloTypes ht(phos.pho_Halo_seedWedge[UsedSidebandPho[0]], 
					phos.pho_Halo_eastNhad[UsedSidebandPho[0]] + phos.pho_Halo_westNhad[UsedSidebandPho[0]]);

			if (ht.IsType(iHaloType)) {
				UsedSidebandPho.clear();
			}
		}

}

int PhotonJets::JetSelAndHistFill(const CommonVars Vars, const PhotonList phos, 
											const ElectronList eles, 
											JetList jets, 
											Hist_t* hist,
											std::vector<int> vPhoIndex, 
											std::vector<int> vEleIndex, 
											std::vector<unsigned int>& count, 
											const bool UseCprWgtForHt, 
											const float fWgt,
											const bool debug,
											const bool bUseJetEtaWeights
											)
{
	int used = 0;
	
	if (vPhoIndex.size() < 1) {
		StdOut(__FILE__,__LINE__,3,"::No photons were selected to make plots!");
		assert(false);
	}
	if (count.size() != 4) {
		StdOut(__FILE__,__LINE__,3,"::counter vector size must be 4!");
		assert(false);
	}
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"::Hist pointer is NULL!");
		assert(false);

	}
	
	if (Vars.met_Ht < fMinPhotonEt+fMinJetEt) return 1; //temp hack to some crazy events! throw them away
	
	int iPhoIndex = vPhoIndex[0];	//assuming the lead photon index is stored first

	++count[0];

	//double check the effect of this cut before uncommenting!
	jets.ApplyEtEtaCuts(fMinJetEt, fMaxJetEt, fMinJetEta, fMaxJetEta);
	
	//for (int i=0; i< jets.jet_num; ++i)
	//{
	//	if (jets.jet_Pt.at(i)<fMinJetEt) assert(false);
	//	if (jets.jet_DetEta.at(i)>fMaxJetEta) assert(false);
	//}

	float fWeight = fWgt; //temp hack for sideband reweighing using both photon Et and jet eta -03-19-2010
	
	if (jets.jet_NJet15 >= (unsigned) iMinNjet15 &&
			jets.jet_NJet15 <= (unsigned) iMaxNjet15)
	{
		int iLeadJetIndex = 0;
		int iSecondLeadJetIndex = 1;

		if (GetApplyMetCleanUpCuts()) 
		{
			if (BadMet(Vars, phos, jets, iPhoIndex, iLeadJetIndex, hist->p1j_Pho.ClosestPhoMetDelPhi_b4, hist->p1j_Jet.ClosestJetMetDelPhi_b4))
			{
				//++iFailCutCount.at(iBadMet);
				++count[3];
				return 1; //discard this event
			}
		}

		if (GetUseDATAMCNvtxWeights())
		{
			const float vtxW = GetDATAMCNvtxWgt(Vars.vtx_NClass12);
			fWeight *= vtxW;
			//std::cout << "Nvtx wgt = " << fWgt << "fWgt -> newWeight =  " << fWeight << std::endl;
		}
		

		++count[1];
		used = 1;
		if (debug) std::cout << __FUNCTION__ << ":" << __LINE__ << " PASS NJET1 = " << jets.jet_NJet15 << std::endl;


		//need to think!
		//I need nominal and delta w to caulculate total delta w
		//and it should be done only in sideband reweighted method == 2?
		
		if (bUseJetEtaWeights)
		{
			const float wtemp = fWeight;
			const float fWeta = GetJetEtaWeight(phos.pho_Etc.at(iPhoIndex), jets.jet_DetEta.at(0));
			fWeight *= fWeta;
			//std::cout << "phoEt/JetEta fWeight -> newWeight =  " << wtemp << " -> " << fWeight << std::endl;
		}

		//temp hack to get Wmc MET weights.
		/*if (GetDataset() == 6 || GetDataset() == 7 || GetDataset() == 8)
		{
			const float x = Vars.met_Met;
			float w = 1;
			if (x>=20 && x<25) w = 0.989131;
			else if (x>=25 && x<30) w = 1.02107;
			else if (x>=30 && x<35) w = 1.03474;
			else if (x>=35 && x<40) w = 0.98105;
			else if (x>=40 && x<45) w = 0.915306;
			else if (x>=45 && x<50) w = 0.948611;
			else if (x>=50 && x<55) w = 0.996498;
			else if (x>=55 && x<60) w = 0.778355;
			//std::cout << "using met weight before  w = " << fWeight << std::endl;
			fWeight *= w;
			//std::cout << "\tusing met weight after   w = " << fWeight << std::endl;
		}
		*/

		

		//FillPhotonHists(Vars, phos, &(hist->p1j_Pho), iPhoIndex, fWeight);
		//FillEventHists(Vars, phos, eles, jets, &(hist->p1j_Evt), iPhoIndex, fWeight, UseCprWgtForHt);
		//FillJetHists(Vars, jets, &(hist->p1j_Jet), iLeadJetIndex, fWeight);		//lead jet
		//FillPhoton1JetHists(Vars, phos, jets, &(hist->p1j_PhoJet),iPhoIndex, iLeadJetIndex, fWeight);


		if (jets.jet_NJet15 >= (unsigned)(iMinNjet15+1)) {		// di-jet case
		/*	int njet30 = 0;
			for (int ij =0; ij < jets.jet_NJet15; ++ij)
			{
				if (jets.jet_Pt.at(ij)<40. || fabs(jets.jet_DetEta.at(ij))>2.4) continue;
				njet30++;
				if (njet30==1) iLeadJetIndex = ij;
				if (njet30==2) iSecondLeadJetIndex = ij;
			}
			if (njet30 != 2) return 1;

			TLorentzVector tlJet1, tlJet2, tlSum;
			tlJet1.SetPxPyPzE(jets.jet_Px[iLeadJetIndex], jets.jet_Py[iLeadJetIndex], jets.jet_Pz[iLeadJetIndex], jets.jet_E[iLeadJetIndex]);
			tlJet2.SetPxPyPzE(jets.jet_Px[iSecondLeadJetIndex], jets.jet_Py[iSecondLeadJetIndex], jets.jet_Pz[iSecondLeadJetIndex], jets.jet_E[iSecondLeadJetIndex]);
			tlSum = tlJet1 + tlJet2; 
			if (tlSum.Pt()<40) return 1;
*/
			//std::cout << "njet30 = " << njet30 << std::endl;
			//std::cout << "lead, second = " << iLeadJetIndex << " [" << jets.jet_Pt.at(iLeadJetIndex) 
			//			<< "], " << iSecondLeadJetIndex << "[" 
			//			<< jets.jet_Pt.at(iSecondLeadJetIndex) << "]" << std::endl;
			++count[2];
			used = 1;

			FillPhotonHists(Vars, phos, &(hist->p2j_Pho), iPhoIndex, fWeight);
			FillEventHists(Vars, phos, eles, jets, &(hist->p2j_Evt), iPhoIndex, fWeight, UseCprWgtForHt);
			FillJetHists(Vars, jets, &(hist->p2j_Jet1), iLeadJetIndex, fWeight);		//lead jet
			FillJetHists(Vars, jets, &(hist->p2j_Jet2),iSecondLeadJetIndex, fWeight);		//2nd lead jet
			FillPhoton1JetHists(Vars, phos, jets, &(hist->p2j_PhoJet1),iPhoIndex, iLeadJetIndex, fWeight);
			//std::cout << __LINE__ << ": nj15=" << jets.jet_NJet15  << "\t iMinNjet15+1 = " << iMinNjet15+1 << std::endl;
			FillPhoton1JetHists(Vars, phos, jets, &(hist->p2j_PhoJet2),iPhoIndex, iSecondLeadJetIndex, fWeight);
			FillPhoton2JetsHists(Vars, phos, jets, &(hist->p2j_PhoJets),iPhoIndex, iLeadJetIndex, iSecondLeadJetIndex, fWeight);
			FillTwoJetsHists(Vars, jets, &(hist->p2j_TwoJets), iLeadJetIndex, iSecondLeadJetIndex, fWeight);
		}
	}

//	std::cout << __LINE__ << ":" << __FUNCTION__ << ": OUT "; PrintHeader(Vars); 
	return used;		//1=this event has been used!
}


void PhotonJets::FillEventHists(const CommonVars& vars, const PhotonList& phos,
								const ElectronList& eles, const JetList& jets, 
								Histograms::EventHists_t* hist, const int iPhoIndex,
								const float fWeight, const bool UseCprWgtForHt)
{ 
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}

	//redo MET AND HT/sumet to see if my calcuation 
	//checks out with  JetFilter 02-22-2010
/*	CommonVars vars = cvars;
	
	//use all objects Et>15
	float sumEt = 0, Ht=0;
	TLorentzVector tlSum(0,0,0,0);
	//only the selected photon should be the taken from photon list.
	//all unused photoons are in the jet list!
		TLorentzVector tlPho(0,0,0,0);
		tlPho.SetPxPyPzE(phos.pho_Px[iPhoIndex], phos.pho_Py[iPhoIndex], phos.pho_Pz[iPhoIndex], phos.pho_E[iPhoIndex]);
		tlSum += tlPho;
		sumEt += tlPho.Pt();
		
	for (int i=0; i<jets.jet_num; ++i)
	{
		TLorentzVector tlJet(0,0,0,0);
		tlJet.SetPxPyPzE(jets.jet_Px[i], jets.jet_Py[i], jets.jet_Pz[i], jets.jet_E[i]);
		if (tlJet.Pt()<15.0) continue;
		if (fabs(jets.jet_DetEta.at(i)) > fMaxJetEta) continue;
		tlSum += tlJet;
		sumEt += tlJet.Pt();
	}

	//no need to worry about electrons. As they are not used in this
	//modules, they are already added to jet list
	vars.met_Met = tlSum.Pt();
	vars.met_Ht = sumEt + tlSum.Pt();
	vars.met_SumEt = sumEt;
*/	
	
	hist->Met->Fill(vars.met_Met, fWeight);
	hist->MetX->Fill(vars.met_MetX, fWeight);
	hist->MetY->Fill(vars.met_MetY, fWeight);
	//hist->Met_Gen_d->Fill(vars.met_Gen_d, fWeight);
	hist->Sumet->Fill(vars.met_SumEt, fWeight);
	hist->NVertices->Fill(vars.vtx_N, fWeight);
	hist->N12Vertices->Fill(vars.vtx_NClass12, fWeight);
	hist->VertexZ->Fill(vars.vtx_z, fWeight);
	hist->NJet15->Fill(jets.jet_NJet15, fWeight);
	hist->MetRun->Fill(vars.evt_RunNumber, vars.met_Met, fWeight);
	hist->SumetRun->Fill(vars.evt_RunNumber, vars.met_SumEt, fWeight);
	hist->NVerticesRun->Fill(vars.evt_RunNumber, vars.vtx_N, fWeight);
	hist->N12VerticesRun->Fill(vars.evt_RunNumber, vars.vtx_NClass12, fWeight);
	hist->NJet15Run->Fill(vars.evt_RunNumber, jets.jet_NJet15, fWeight);
	hist->NJet15VsMet->Fill(vars.met_Met, jets.jet_NJet15, fWeight);
	hist->NJet15VsSumet->Fill(vars.met_SumEt, jets.jet_NJet15, fWeight);
	hist->N12VerticesVsMet->Fill(vars.met_Met, vars.vtx_NClass12, fWeight);
	hist->N12VerticesVsSumet->Fill(vars.met_SumEt, vars.vtx_NClass12, fWeight);
	hist->N12VerticesMetProf->Fill(vars.met_Met, vars.vtx_NClass12, fWeight);
	hist->N12VerticesSumetProf->Fill(vars.met_SumEt, vars.vtx_NClass12, fWeight);
	if (vars.tri_pho25iso == 1) hist->Triggers->Fill(0);
	if (vars.tri_pho50 == 1) hist->Triggers->Fill(1);
	if (vars.tri_pho70 == 1) hist->Triggers->Fill(2);

	//if (vars.met_Ht <45.)
	if (vars.met_Ht < fMinPhotonEt+fMinJetEt) //temp hack to some crazy events! throw them away
	{
		std::cout << __FILE__ <<":" << __FUNCTION__ << ":" << __LINE__ <<	"::Event with Ht = " << vars.met_Ht << "<" << fMinPhotonEt+fMinJetEt << std::endl;
		PrintHeader(vars); 
	}
	
	hist->Ht->Fill(vars.met_Ht, fWeight);
	if (UseCprWgtForHt)			// these will be filled only if CPR weight is required
	{
		
//		hist->HtWgt->Fill(vars.met_Ht, phos.pho_CprWgt.at(iPhoIndex));	
//		hist->HtSys[0]->Fill(vars.met_Ht, phos.pho_CprSys1.at(iPhoIndex));
//		hist->HtSys[1]->Fill(vars.met_Ht, phos.pho_CprSys2.at(iPhoIndex));
//		hist->HtSys[2]->Fill(vars.met_Ht, phos.pho_CprSys3.at(iPhoIndex));
//		hist->HtSys[3]->Fill(vars.met_Ht, phos.pho_CprSys4.at(iPhoIndex));
//		hist->HtSys[4]->Fill(vars.met_Ht, phos.pho_CprSys5.at(iPhoIndex));
//		hist->HtSys[5]->Fill(vars.met_Ht, phos.pho_CprSys6.at(iPhoIndex));
//		hist->HtSys[6]->Fill(vars.met_Ht, phos.pho_CprSys7.at(iPhoIndex));
//		hist->HtSys[7]->Fill(vars.met_Ht, phos.pho_CprSys8.at(iPhoIndex));
	}




}

void PhotonJets::FillPhotonHists(const CommonVars& cVars, const PhotonList& phos, Histograms::PhotonHists_t* hist,
								const int iPhoIndex, const float fWeight, const bool debug)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}

	if (iPhoIndex<0 || iPhoIndex >= phos.pho_Etc.size())
	{
		StdOut(__FILE__,__LINE__,3,"Photon index out bounds!");
		assert(false);
	}

	hist->Detector->Fill(phos.pho_Detector[iPhoIndex], fWeight);
	hist->DetEta->Fill(phos.pho_DetEta[iPhoIndex], fWeight);
	hist->DetPhi->Fill(phos.pho_DetPhi[iPhoIndex], fWeight);
	hist->EtCorr->Fill(phos.pho_Etc[iPhoIndex], fWeight);
	hist->EventWeights->Fill(phos.pho_Etc[iPhoIndex], fWeight);
	//if (debug && phos.pho_Etc.at(iPhoIndex)>70)
	if (debug)
	{
		std::cout << __FUNCTION__ << " wgt = " << fWeight << std::endl;
		if (fWeight == 0) std::cout << " ************** Wgt == 0 " << std::endl;
			std::cout << "pho Et -> fWgt = " << phos.pho_Etc.at(iPhoIndex)<< "-> " << fWeight  << std::endl;
			hist->EtCorr->Print();
	}
	hist->XCes->Fill(phos.pho_XCes[iPhoIndex], fWeight);
	hist->ZCes->Fill(phos.pho_ZCes[iPhoIndex], fWeight);
	hist->HadEm->Fill(phos.pho_HadEm[iPhoIndex], fWeight);
	hist->Chi2Mean->Fill(phos.pho_Chi2Mean[iPhoIndex], fWeight);
	hist->N3d->Fill(phos.pho_N3d[iPhoIndex], fWeight);
	hist->Iso4->Fill(phos.pho_Iso4[iPhoIndex], fWeight);
	hist->TrkPt->Fill(phos.pho_TrkPt[iPhoIndex], fWeight);
	hist->TrkIso->Fill(phos.pho_TrkIso[iPhoIndex], fWeight);
	hist->Ces2Wire->Fill(phos.pho_CesWireE2[iPhoIndex], fWeight);
	hist->Ces2Strip->Fill(phos.pho_CesStripE2[iPhoIndex], fWeight);
	hist->EmTime->Fill(phos.pho_EmTime[iPhoIndex], fWeight);
	hist->PhiWedge->Fill(phos.pho_PhiWedge[iPhoIndex], fWeight);
	hist->EmTimeVsRun->Fill(cVars.evt_RunNumber, phos.pho_EmTime[iPhoIndex], fWeight);
	hist->EtCorrVsRun->Fill(cVars.evt_RunNumber, phos.pho_Etc[iPhoIndex], fWeight);
	if (phos.pho_CprWgt[iPhoIndex] != 0) //0 weight means the photon failed the 
	{												// CES/CPR fiducial cuts. see elog#1491 and #1495
													// these photon need to be excluded to get the
													// correct photon fraction. 
		hist->CprWeight->Fill(phos.pho_CprWgt[iPhoIndex]);	
	}
	TLorentzVector tlPhoVec(0,0,0,0);
	TVector2 tv2MetVec(cVars.met_MetX, cVars.met_MetY);
	tlPhoVec.SetPxPyPzE(phos.pho_Px[iPhoIndex], phos.pho_Py[iPhoIndex], phos.pho_Pz[iPhoIndex], phos.pho_E[iPhoIndex]);

	//from V7 stuples and above has complete MET vector to do this
	const float dphi_pm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlPhoVec.Phi())));
	hist->PhoMetDelPhi->Fill(dphi_pm);
	hist->ClosestPhoMetDelPhi_a4->Fill(dphi_pm);
	
	hist->PhoEtMetProf->Fill(cVars.met_Met, phos.pho_Etc[iPhoIndex], fWeight);
	hist->PhoEMetProf->Fill(cVars.met_Met, tlPhoVec.E(), fWeight);
	hist->PhoEtaMetProf->Fill(cVars.met_Met, tlPhoVec.Eta(), fWeight);
}
	

void PhotonJets::FillJetHists(const CommonVars& cVars, const JetList& jets, Histograms::JetHists_t* hist, 
								const int iJetIndex,	const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	hist->Emfr->Fill(jets.jet_Emfr[iJetIndex], fWeight);
	hist->DetEta->Fill(jets.jet_DetEta[iJetIndex], fWeight);
	hist->NTracks->Fill(jets.jet_Ntracks[iJetIndex], fWeight);
	if (jets.jet_Ntowers[iJetIndex]==0)
	{
		std::cout << __FUNCTION__ << "::" << __LINE__ << "A jet with ntower==0 found!!! ";
		//cVars.PrintHeader();
		//jets.DumpAll();
	}
	hist->NTowers->Fill(jets.jet_Ntowers[iJetIndex], fWeight);
	hist->EtCorr->Fill(jets.jet_Pt[iJetIndex], fWeight);

	TLorentzVector tlVec(jets.jet_Px[iJetIndex], jets.jet_Py[iJetIndex], jets.jet_Pz[iJetIndex], jets.jet_E[iJetIndex]);
	hist->EvtPhi->Fill(tlVec.Phi(), fWeight);
	hist->EvtEta->Fill(tlVec.Eta(), fWeight);
	//hist->SeedIEta->Fill(jets.jet_SeedIEta[iJetIndex], fWeight);
	//hist->SeedIPhi->Fill(jets.jet_SeedIPhi[iJetIndex], fWeight);
	
	TVector2 tv2MetVec(cVars.met_MetX, cVars.met_MetY);

	//from V7 stuples and above has complete MET vector to do this
	const float dphi_jm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlVec.Phi())));
	hist->JetMetDelPhi->Fill(dphi_jm);
	hist->ClosestJetMetDelPhi_a4->Fill(dphi_jm);
	
	hist->JetEtMetProf->Fill(cVars.met_Met, tlVec.Pt(), fWeight);
	hist->JetEMetProf->Fill(cVars.met_Met, tlVec.E(), fWeight);
	hist->JetEtaMetProf->Fill(cVars.met_Met, tlVec.Eta(), fWeight);
}

void PhotonJets::FillPhoton1JetHists(const CommonVars& cVars, const PhotonList& phos, const JetList& jets, 
							Histograms::Photon1JetHists_t* hist, const int iPhoIndex, const int iJetIndex,
							const float fWeight)
{
	//std::cout << "\tfillpho1jethist:: npho=" << phos.pho_num << std::endl;
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	if (iPhoIndex < 0 || iPhoIndex > (int)phos.pho_num)
	{
		StdOut(__FILE__,__LINE__,3,"Invalid Photon index!");
		assert(false);
	}
	if (iJetIndex < 0 || iJetIndex > (int)jets.jet_NJet15)
	{
		StdOut(__FILE__,__LINE__,3,"Invalid Jet index!");
		assert(false);
	}
	TLorentzVector tlPho, tlJet, tlSum;
	float InvMass 	= 0;
	float DelPhi 	= 0;
	float DelEta 	= 0;
	float DelR 		= 0;
	float EtRatio 	= 0;
	float JetEmFr 	= 0;
	float PhoPhiWedge = 0; 

//	std::cout << "-----------jet dump in hist fill " << std::endl;
//	jets.Dump();

	tlJet.SetPxPyPzE(jets.jet_Px[iJetIndex], jets.jet_Py[iJetIndex], jets.jet_Pz[iJetIndex], jets.jet_E[iJetIndex]);
	tlPho.SetPxPyPzE(phos.pho_Px[iPhoIndex], phos.pho_Py[iPhoIndex], phos.pho_Pz[iPhoIndex], phos.pho_E[iPhoIndex]);

	tlSum = tlPho + tlJet; 
	InvMass 	= tlSum.M();
	if (InvMass< 1) 
	{
		std::cout << __FILE__ <<":" << __FUNCTION__ << ":" << __LINE__ <<	"::Event with InvMass = " << InvMass << " <5" << std::endl;
		std::cout << "Run,Event Number::" << cVars.evt_RunNumber << ", " << cVars.evt_EventNumber << std::endl;
		std::cout << "Jet Pt/E = " << tlJet.Pt() << ", " << tlJet.E() << std::endl;
		std::cout << "Pho Pt/E = " << tlPho.Pt() << ", " << tlPho.E() << std::endl;
		std::cout << "Sum Pt/E = " << tlSum.Pt() << ", " << tlSum.E() << std::endl;
		phos.Dump();
		jets.Dump();
		PrintHeader(cVars);
		stuple.Dump(0);
	}
	
	DelPhi 	= fabs(tlPho.DeltaPhi(tlJet));
	DelEta 	= fabs(phos.pho_DetEta[iPhoIndex] - jets.jet_DetEta[iJetIndex]);
	DelR 		= tlPho.DeltaR(tlJet);
	EtRatio 	= (tlJet.Energy() * TMath::Sin(tlJet.Theta())) / (tlPho.Energy() * TMath::Sin(tlPho.Theta()));
	//JetEmFr 	= jets.jet_Emfr[iJetIndex];
	PhoPhiWedge = phos.pho_PhiWedge[iPhoIndex];
	
	hist->InvMass->Fill(InvMass, fWeight);
	hist->DelPhi->Fill(DelPhi, fWeight);
	hist->DelEta->Fill(DelEta, fWeight);
	hist->DelR->Fill(DelR, fWeight);
	hist->EtRatio->Fill(EtRatio, fWeight);
	//hist->JetEmfrVsDelPhi->Fill(DelPhi, JetEmFr, fWeight);
	//hist->JetEmfrVsPhoPhiWedge->Fill(PhoPhiWedge, JetEmFr, fWeight);
	
	//photon jet balance
	//std::cout << "Pho tlPt, Etc = " << tlPho.Pt() << ", " << phos.pho_Etc.at(iPhoIndex) << std::endl;
	//std::cout << "Jet tlPt, Etc = " << tlJet.Pt() << ", " << jets.jet_Pt.at(iJetIndex) << std::endl;
	hist->PhoJetPtBalance->Fill(tlJet.Pt()/tlPho.Pt() - 1, fWeight);
	hist->PhoJetEBalance->Fill(tlJet.E()/tlPho.E() - 1, fWeight);

}


void PhotonJets::FillPhoton2JetsHists(const CommonVars& cVars, const PhotonList& phos,
						const JetList& jets, Histograms::Photon2JetsHists_t* hist,
						const int iPhoIndex, const int iJet1Index,
						const int iJet2Index, const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	TLorentzVector tlJet1, tlJet2, tlPho, tlSum;

	tlJet1.SetPxPyPzE(jets.jet_Px[iJet1Index], jets.jet_Py[iJet1Index], jets.jet_Pz[iJet1Index], jets.jet_E[iJet1Index]);
	tlJet2.SetPxPyPzE(jets.jet_Px[iJet2Index], jets.jet_Py[iJet2Index], jets.jet_Pz[iJet2Index], jets.jet_E[iJet2Index]);
	tlPho.SetPxPyPzE(phos.pho_Px[iPhoIndex], phos.pho_Py[iPhoIndex], phos.pho_Pz[iPhoIndex], phos.pho_E[iPhoIndex]);
	
	tlSum = tlPho + tlJet1 + tlJet2; 
	hist->InvMass->Fill(tlSum.M(), fWeight);
	hist->Pt->Fill(tlSum.Pt(), fWeight);
	hist->PtBal->Fill((tlJet1+tlJet2).Pt()/tlPho.Pt(), fWeight);
	hist->PtBalVsMass->Fill((tlJet1+tlJet2).M(), tlPho.Pt()/(tlJet1+tlJet2).Pt());

}

void PhotonJets::FillTwoJetsHists(const CommonVars& cVars, const JetList& jets,
						Histograms::TwoJetsHists_t* hist,
						const int iJet1Index, const int iJet2Index,
						const float fWeight)
{
	if (!hist)
	{
		StdOut(__FILE__,__LINE__,3,"Hist pointer is null!");
		assert(false);
	}
	TLorentzVector tlJet1, tlJet2, tlSum;
	float DelEta = 0;

	tlJet1.SetPxPyPzE(jets.jet_Px[iJet1Index], jets.jet_Py[iJet1Index], jets.jet_Pz[iJet1Index], jets.jet_E[iJet1Index]);
	tlJet2.SetPxPyPzE(jets.jet_Px[iJet2Index], jets.jet_Py[iJet2Index], jets.jet_Pz[iJet2Index], jets.jet_E[iJet2Index]);
	DelEta = fabs(jets.jet_DetEta[iJet1Index] - jets.jet_DetEta[iJet2Index]);
	
	tlSum = tlJet1 + tlJet2; 
	float fEtRatio = (tlJet2.Energy() * TMath::Sin(tlJet2.Theta())) / (tlJet1.Energy() * TMath::Sin(tlJet1.Theta()));
	float DelPhi = fabs(tlJet1.DeltaPhi(tlJet2));
	float DelR = tlJet1.DeltaR(tlJet2);
	
	hist->InvMass->Fill(tlSum.M(), fWeight);
	hist->DelPhi->Fill(DelPhi, fWeight);
	hist->DelEta->Fill(DelEta, fWeight);
	hist->DelR->Fill(DelR, fWeight);
	hist->EtRatio->Fill(fEtRatio, fWeight);
	hist->Pt->Fill(tlSum.Pt(), fWeight);

}



int PhotonJets::Prepass(const JetList jets)
{
	if (jets.jet_NJet15<2) return 0;
	
	int iJet1Index=0, iJet2Index=1;
	TLorentzVector tlJet1,tlJet2, tlSum;
	tlJet1.SetPxPyPzE(jets.jet_Px[iJet1Index], jets.jet_Py[iJet1Index], jets.jet_Pz[iJet1Index], jets.jet_E[iJet1Index]);
	tlJet2.SetPxPyPzE(jets.jet_Px[iJet2Index], jets.jet_Py[iJet2Index], jets.jet_Pz[iJet2Index], jets.jet_E[iJet2Index]);
	
	tlSum = tlJet1 + tlJet2; 

	if (tlSum.M() > 50.0 && tlSum.M() < 200.0) return 1;
	
	return 0;	

}

void PhotonJets::PrintHeader(const CommonVars& cVars) const
{
	std::cout << " ============== " << cVars.evt_RunNumber << ", " << cVars.evt_EventNumber << std::endl;
}

void PhotonJets::RemoveDuplicateEMobjects(PhotonList& phos, ElectronList& eles)
{
	//this is to account for the special case of an EM Object passing 
	//both loose photon id cuts and loose electron id cuts  (tight id
	//cuts are failed in both cases it seems)so it ends up in both
	//PhotonList and in ElectronList. When this rooutine is called
	//the loose electron will be added to jet list and if it becomes
	//the leading jet I'll be adding the same object to derive the
	//invariant mass, which becomes <0? Dec 12-2009
	//So I'll try to delete the copy form electron list

	std::vector<int> vEleToDelete;
	if (phos.pho_num>=1 && eles.ele_num>=1)
	{
		for (unsigned int ipho = 0; ipho<phos.pho_num; ++ipho)
		{
			for (unsigned int iele = 0; iele<eles.ele_num; ++iele)
			{

				if (phos.pho_Etc.at(ipho) == eles.ele_Etc.at(iele))
				{
					vEleToDelete.push_back(iele);
				}

			}

		}

	}
	
	//now removed the duplicate em objects
	for (unsigned int i=0; i< vEleToDelete.size(); ++i)
	{
		//stuple.Dump();
		eles.Delete(vEleToDelete.at(i));
	}

}

void PhotonJets::PrintJobSettings()
{	

	std::cout << "Data Written to        = " << GetHistFileName() << std::endl;
	std::cout << "Add Iso4(2) to QCD Pho = ";
	if (GetAddIsoEtoSidebandPho())
	{
		std::cout << GetAddIsoEtoSidebandPho() << " (YES)" << std::endl; 
	} else 
	{
		std::cout << GetAddIsoEtoSidebandPho() << " (NO)" << std::endl; 
	}
	std::cout << "Reweight Sideband Pho  = ";
	if (GetReweightSideband())
	{
		if (GetReweightSideband() == 1)
		{
			std::cout << GetReweightSideband() << " (NOMINAL WEIGHTS) : ";
			tf1_sideband_phoet->Print(); 
			std::cout << std::endl;
			tf1_sideband_jeteta->Print(); 
		} else if (GetReweightSideband() == 2)
		{
			std::cout << GetReweightSideband() << " (+SIGMA WEIGHTS) : " ; 
			tf1_sideband_phoet->Print(); 
			std::cout << std::endl;
			tf1_sideband_jeteta->Print(); 
		} else if (GetReweightSideband() == 3)
		{
			std::cout << GetReweightSideband() << " (PHO ET ONLY WEIGHTS) : ";
			tf1_sideband_phoet->Print(); 
			std::cout << std::endl;
		} else if (GetReweightSideband() == 4)
		{
			std::cout << GetReweightSideband() << " (+SIGMA PHO ET ONLY WEIGHTS) : " ; 
			tf1_sideband_phoet_error->Print(); 
			std::cout << std::endl;

		} else
		{
			std::cout << GetReweightSideband() << " (ERROR! UNKNOWN WEIGHTS!)" << std::endl; 
		}
	} else 
	{
		std::cout << GetReweightSideband() << " (NO)" << std::endl; 
	}

	
	std::cout << "MC Flag setting        = ";
	if (bMcFlag) std::cout << "MC" << std::endl;
	else std::cout << "DATA" << std::endl;
	std::cout << "min,max Vtx for Signal = " << std::setw(4) << GetMinClass12Vtx() << ", " << std::setw(4) << GetMaxClass12Vtx() << std::endl;
	std::cout << "Min,max njet15         = " << std::setw(4) << GetMinNjet15() << ", " << std::setw(4) << GetMaxNjet15() << std::endl;
	std::cout << "Photon Et(min,max)     = " << std::setw(4) << GetMinPhotonEt() << ", " << std::setw(5) << GetMaxPhotonEt() << std::endl;
	std::cout << "Photon Eta(max)        = " << std::setw(4) << GetMaxPhotonEta() << std::endl;
	std::cout << "Jet Et(min,max)        = " << std::setw(4) << GetMinJetEt() << ", " << std::setw(4)<< GetMaxJetEt() << std::endl;
	std::cout << "Jet Eta (min,max)      = " << std::setw(4) << GetMinJetEta() << ", " << std::setw(4)<< GetMaxJetEta() << std::endl;
	std::cout << "Signal Time Window     = " << std::setw(4) << GetSignalTimeLoLimit() << ", " << std::setw(4) << GetSignalTimeUpLimit() << std::endl;
	std::cout << "Cosmic Time Window     = " << std::setw(4) << GetCosmicTimeLoLimit() << ", " << std::setw(4) << GetCosmicTimeUpLimit() << std::endl;
	std::cout << "Halo Type              = " << std::setw(4) << GetHaloType() << std::endl;
	std::cout << "Min/max MEt            = " << std::setw(4) << GetMinMet() << ", " << std::setw(4) << GetMaxMet()<< std::endl;
	std::cout << "MEt cleanup cuts?      = " << std::setw(4) << GetApplyMetCleanUpCuts()  << " (1=YES)"<< std::endl;
	std::string sNvtxSetting("UNKNOWN");
	std::cout << "Nvtx setting           = " << GetUseDATAMCNvtxWeights() << std::endl;
	if (GetUseDATAMCNvtxWeights() == 1) sNvtxSetting = "CENTRAL WEIGHTS";
	else if (GetUseDATAMCNvtxWeights() == 2) sNvtxSetting = "SIGMA WEIGHTS";
	std::cout << "Nvtx Weighting?        = " << sNvtxSetting << std::endl;
	if (UseBTagJets()) std::cout << "b-Tag jets required?   = YES" << std::endl;
	else std::cout << "b-Tag jets required?   = NO" << std::endl;
	std::cout << "PrintLevel             = " << std::setw(4) << PrintLevel() << std::endl;
}


void PhotonJets::DoHaloRejHists(const Stuple& stuple)
{

	if (  (stuple.tri_pho25iso != 1) && (stuple.tri_pho50 != 1) 	//_______________ only data has L1-3 trig. info
				&& (stuple.tri_pho70 != 1) ) 										//_______________ mc has L1-2.
			return;												
			
	
		//______________________________________________________________ select the photon
		for (unsigned int i = 0; i < stuple.pho_num; ++i )
		{
			if (stuple.pho_PhoenixId[i] != 0) continue;  //____________________________ reject phoenix photons (this is redundant as nvtx==0 guarantees this)
			//________________________________________________________________________ this is a hack to find the trigger photon, as I don't have the L3 info to match reconstructed photon to trigger photon.
			if (stuple.tri_pho70 == 1 && stuple.pho_Etc[i] < 70.0) continue;
			if (stuple.tri_pho50 == 1 && stuple.pho_Etc[i] < 50.0) continue;
			if (stuple.tri_pho25iso == 1 && stuple.pho_Etc[i] < 25.0) continue;

			if (stuple.pho_NMuonStubs[i] != 0) continue;
			if (fabs(stuple.pho_DetEta[i])>fMaxPhotonEta) continue;
			if (stuple.pho_Etc[i]<fMinPhotonEt || stuple.pho_Etc[i]>fMaxPhotonEt) continue;		//need this as I am tagging Et>7GeV em objects. 
			if (stuple.pho_EmTime[i] < fMinSignalEmTime || stuple.pho_EmTime[i] > fMaxSignalEmTime) continue;
			if (stuple.pho_TightId[i] == 0)
			{
				//do the signal halo rejection power
				haloPhiWedge_1j[0]->Fill(stuple.pho_PhiWedge[i]);
				haloEmTime_1j[0]->Fill(stuple.pho_EmTime[i]);
				haloNvtx[0]->Fill(stuple.vtx_NClass12);
				
				//FillPhotonHists(cVars, cPhos, HhaloNoVtx_Evt_b4, i);
				//FillEventHists(cVars, cPhos, cEles, cJets, HhaloNoVtx_Pho_b4, i);

				
				HaloTypes ht(stuple.pho_Halo_seedWedge[i], stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);

				if (ht.IsType(iHaloType))
				{
					haloPhiWedge_1j[1]->Fill(stuple.pho_PhiWedge[i]);
					haloEmTime_1j[1]->Fill(stuple.pho_EmTime[i]);
					haloNvtx[1]->Fill(stuple.vtx_NClass12);
					//FillPhotonHists(cVars, cPhos, HhaloNoVtx_Evt_a4, i);
					//FillEventHists(cVars, cPhos, cEles, cJets, HhaloNoVtx_Pho_a4, i);
				}

				break;  //no more porcessing of this events
			} else if (stuple.pho_TightId[i] > 0 && stuple.pho_LooseId[i] == 0)			// >0 means that this has failed tight cuts. must do this as I init them to large negative values and != could mean <0!!
			{
				//do the sideband halo rejections power
				haloPhiWedge_2j[0]->Fill(stuple.pho_PhiWedge[i]);
				haloEmTime_2j[0]->Fill(stuple.pho_EmTime[i]);
				
				HaloTypes ht(stuple.pho_Halo_seedWedge[i], stuple.pho_Halo_eastNhad[i] + stuple.pho_Halo_westNhad[i]);

				if (ht.IsType(iHaloType))
				{
					haloPhiWedge_2j[1]->Fill(stuple.pho_PhiWedge[i]);
					haloEmTime_2j[1]->Fill(stuple.pho_EmTime[i]);
				}

				break;  //no more porcessing of this events
			}
		}	

}

/*w0 is very peculiar to use it byself to do calucations.
 *
 *
 *
 */
void PhotonJets::DoHaloRejCalc()
{	
	std::cout << magenta << "========= halo rejection power and estimates ==== " << std::endl;
	std::cout << "Halo in tight photon sample ... " << std::endl;
	if (Hhalo.p1j_Pho.EtCorr->GetEntries()>0)
	{
		std::cout << " Identified tight halo +1jet events   = " << Hhalo.p1j_Pho.EtCorr->Integral() << std::endl;
		//DoHaloRejCalc(haloPhiWedge_1j[0], haloPhiWedge_1j[1], Hhalo.p1j_Pho.EtCorr->GetIntegral());
	} else std::cout << "NO DATA" << std::endl;
	std::cout << "Halo insideband photon sample ... " << std::endl;
	if (Hsideband_halo.p1j_Pho.EtCorr->GetEntries()>0)
	{
		std::cout << " Identified sideband halo+1jet events = " << Hsideband_halo.p1j_Pho.EtCorr->GetIntegral() << std::endl;
		//DoHaloRejCalc(haloPhiWedge_2j[0], haloPhiWedge_2j[1], Hsideband_halo.p1j_Pho.EtCorr->GetIntegral());
	} else std::cout << "NO DATA" << std::endl;
	std::cout << "=== end of halo rejection power and estimates ==== " << clearatt << std::endl;

}

void PhotonJets::DoHaloRejCalc(const TH1F *haloPhiWedge_b4, const TH1F *haloPhiWedge_a4, const Double_t Nid)
{	
	double Nw0_b = haloPhiWedge_b4->GetBinContent(1);
	double Nw23_b = haloPhiWedge_b4->GetBinContent(24);
	double Nw0_a = haloPhiWedge_a4->GetBinContent(1);
	double Nw23_a = haloPhiWedge_a4->GetBinContent(24);
	double Nw1_22_b = 0;
	double Nw1_22_a = 0;

	for (int bin=2; bin < 24; ++bin)
	{
		Nw1_22_b += haloPhiWedge_b4->GetBinContent(bin);
		Nw1_22_a += haloPhiWedge_a4->GetBinContent(bin);
	}

	std::cout << "before " << std::endl;
	std::cout << "w1 = " << Nw0_b << std::endl;
	std::cout << "w23 = " << Nw23_b << std::endl;
	std::cout << "w1_22 = " << Nw1_22_b  << " (" << Nw1_22_b/22. << ")"<< std::endl;
	double Nh_b = Nw0_b + Nw23_b - 2*Nw1_22_b/22.;
	std::cout << "total before = " << Nh_b << std::endl;

	std::cout << "after " << std::endl;
	std::cout << "w1 = " << Nw0_a << std::endl;
	std::cout << "w23 = " << Nw23_a << std::endl;
	std::cout << "w1_22 = " << Nw1_22_a  << " (" << Nw1_22_a/22. << ")"<< std::endl;
	double Nh_a = Nw0_a + Nw23_a - 2*Nw1_22_a/22.;
	std::cout << "total after = " << Nh_a << std::endl;
	double rejpow = 0;
	if (Nh_b>0) rejpow = Nh_a/Nh_b * 100;
	std::cout << " Halo rej power = " << rejpow << std::endl;
	std::cout << " Halo est       = " << Nid * (100-rejpow)/100 << std::endl;

/*	//now do the systeamtic
	// I assume wedge0 is 100% halo
	const double Nh_b_s1 = Nw0_b + Nw23_b - Nw1_22_b/22.;
	const double Nh_a_s1 = Nw0_a + Nw23_a - Nw1_22_a/22.;
	double rejpow_s1 = 0;
	if (Nh_b_s1>0) rejpow_s1 = Nh_a_s1/Nh_b_s1 * 100;
	// Now assume wedge23 is 100% halo
	const double Nh_b_s2 = Nw0_b + Nw23_b - Nw1_22_b/22.;
	const double Nh_a_s2 = Nw0_a + Nw23_a - Nw1_22_a/22.;
	double rejpow_s1 = 0;
	if (Nh_b_s1>0) rejpow_s1 = Nh_a_s1/Nh_b_s1 * 100;
*/	


}

/*****************************************************
 *this calculated the cosmic photon background
 * at the end of job running over data
 ****************************************************/

void PhotonJets::GetCosmicEstAndErr()
{
	std::cout << blue << __FUNCTION__ << " :: Cosmic Photon Estimate Summary ****" << std::endl;
	GetCosmicEstAndErr(Hcosmic.p1j_Pho.EmTime, Hcosmic.p2j_Pho.EmTime);
	std::cout << __FUNCTION__ << " :: Sideband Cosmic Photon Estimate Summary ****" << std::endl;
	GetCosmicEstAndErr(Hsideband_cosmic.p1j_Pho.EmTime, Hsideband_cosmic.p2j_Pho.EmTime);
	std::cout << __FUNCTION__ << " :: END Cosmic Photon Estimate Summary " << clearatt << std::endl;
}
void PhotonJets::GetCosmicEstAndErr(const TH1 *hj1, const TH1 *hj2)
{
	if  (hj1 == NULL || hj2 == NULL)
	{ 
		std::cout << __FUNCTION__ << ": One of ths input hists is null. returning!" << std::endl;
		return;
	}
	

	const float min_time= GetCosmicTimeLoLimit();
	const float max_time= GetCosmicTimeUpLimit();
	const float mid_time = fabs(min_time - max_time)/2.0;
	const float sig_time_window = fabs(GetSignalTimeUpLimit() - GetSignalTimeLoLimit());
	
	double tot_1j = 0, tot_2j = 0;
	double half_1j = 0, half_2j = 0;
	
	for (int bin=1; bin <= hj1->GetNbinsX(); ++bin)
	{
		if (hj1->GetBinContent(bin))
		{
			tot_1j += hj1->GetBinContent(bin);
			//std::cout << hj1->GetBinLowEdge(bin) << "\t" << (min_time + mid_time) << std::endl;
			if (hj1->GetBinLowEdge(bin)< (min_time + mid_time))
			{
				//std::cout << "\t" << hj1->GetBinContent(bin) << std::endl;
				half_1j += hj1->GetBinContent(bin);
			}
		}
		if (hj2->GetBinContent(bin))
		{
			tot_2j += hj2->GetBinContent(bin);
			if (hj2->GetBinLowEdge(bin)< (min_time + mid_time))
			{
				half_2j += hj2->GetBinContent(bin);
			}
		}

	}
	
	const float est_1j = sig_time_window * ( tot_1j / (max_time - min_time));
	const float est_2j = sig_time_window * ( tot_2j / (max_time - min_time));
	const float estsyst_1j = sig_time_window * ( half_1j / mid_time);
	const float estsyst_2j = sig_time_window * ( half_2j / mid_time);

	std::cout << "sig widow, cosmic time window, half width = " 
		<< sig_time_window << ", " << max_time - min_time << ", "
		<< mid_time << std::endl;
	std::cout << "cosmic est 1j = " << est_1j << "(" << estsyst_1j <<  ") +/- " << fabs(est_1j -estsyst_1j) << std::endl;
	std::cout << "cosmic est 2j = " << est_2j << "(" << estsyst_2j << ") +/- "  << fabs(est_2j -estsyst_2j )<< std::endl;
}

double PhotonJets::GetJetEtaWeight(const float phoEt, const float jetEta)
{
		//std::cout << "In " << __FUNCTION__ << std::endl;
		const float x = phoEt;
		const float y = jetEta;
		

	if (GetReweightSideband() == 1)
	{
		//const float W = tf1_sideband_phoet->Eval(x) * tf1_sideband_jeteta->Eval(y);
		//hSidebandJetEtaWeights->Fill(jetEta, tf1_sideband_jeteta->Eval(jetEta));
		//hSidebandPhoEtWeights->Fill(phoEt, tf1_sideband_phoet->Eval(phoEt));
		

				
	/*	std::cout << "nominal weight for phoEt("<< phoEt << ") and jetEta("<<jetEta 
			<< ") = " << W 
			<< "    wEt = " << tf1_sideband_phoet->Eval(x)
			<< ", wEta = " << tf1_sideband_jeteta->Eval(y)
			<< std::endl;
	*/	

		float w_eta = 0, w_eta_err=0;
		GetHistBin(hist_JetEtaWghts, jetEta, w_eta, w_eta_err);
		const float W = tf1_sideband_phoet->Eval(x) * w_eta ;
		hSidebandJetEtaWeights->Fill(jetEta, w_eta);
		hSidebandPhoEtWeights->Fill(phoEt, tf1_sideband_phoet->Eval(phoEt));
		return W;

	} else if (GetReweightSideband() == 2)
	{

		float w_eta = 0, w_eta_err=0;
		GetHistBin(hist_JetEtaWghts, jetEta, w_eta, w_eta_err);
		hSidebandJetEtaWeights->Fill(jetEta, w_eta);
		hSidebandPhoEtWeights->Fill(phoEt, tf1_sideband_phoet->Eval(phoEt));

		const float w_pet     = tf1_sideband_phoet->Eval(phoEt);
		const float w_pet_alt = tf1_sideband_phoet_error->Eval(phoEt);
		const float w_pet_delta = fabs(w_pet - w_pet_alt);

		const float W = w_pet * w_eta;
		const float delW = sqrt( pow(w_pet * w_eta_err,2) + pow(w_eta * w_pet_delta,2));
		//std::cout << __LINE__ << "W+delW = " << W << "+/- " << delW << std::endl;
		return W+delW;

		
/*		const double funcVal = tf1_sideband_jeteta->Eval(jetEta);
		//std::cout << __LINE__<< ": func val = " << funcVal << std::endl;
		const float p0 = tf1_sideband_jeteta->GetParameter(0);
		const float p1 = tf1_sideband_jeteta->GetParameter(1);
		const float p2 = tf1_sideband_jeteta->GetParameter(2);
		const float p3 = tf1_sideband_jeteta->GetParameter(3);
		const float p4 = tf1_sideband_jeteta->GetParameter(4);
		const float p5 = tf1_sideband_jeteta->GetParameter(5);
		const float p6 = tf1_sideband_jeteta->GetParameter(6);
		const float p7 = tf1_sideband_jeteta->GetParameter(7);
		const float p8 = tf1_sideband_jeteta->GetParameter(8);
		
		//g is the jet eta reweight function
		const float g1 = p0 * exp( -0.5 * pow((y-p1)/p2,2) );
		const float g2 = p3 * exp( -0.5 * pow((y-p4)/p5,2) );
		const float g3 = p6 * exp( -0.5 * pow((y-p7)/p8,2) );
		const double g = g1 + g2 + g3;
		
		//std::cout << __LINE__<< ": my   val = " << g << std::endl;
		if ((g - funcVal) > 0.000001)
		{
			std::cout << red << __LINE__ <<  ": my func and ROOT func vals disagree!" << std::endl;
		}

		//1st derivative of g

		const float g1prime = g1 * ( (p1-y)/pow(p2,2) );
		const float g2prime = g2 * ( (p4-y)/pow(p5,2) );
		const float g3prime = g3 * ( (p7-y)/pow(p8,2) );
		const float gprime  =  g1prime + g2prime + g3prime;

		const float deltag = fabs(g - gprime);

		
		//f is the photon et weight function

		const float f = tf1_sideband_phoet->Eval(phoEt);
		const float falt = tf1_sideband_phoet_error->Eval(phoEt);
		const float deltaf = fabs(f-falt);

		//total weight
		const float W = f * g;

		//now the error weight
		const float deltaW = falt * gprime * deltag;
		const float Walt   =  W + deltaW;
		

		hSidebandJetEtaWeights->Fill(jetEta, funcVal);
		hSidebandPhoEtWeights->Fill(phoEt, falt);
		
		return Walt;
		
*/		
/*		
		TF1 *final = dynamic_cast<TF1*> (func->Clone("func_copy"));																					;

		const int nPar = final->GetNpar();
		const double fValMin = final->Eval(x);
		const double fValMin2 = func->Eval(x);
		assert (fValMin == fValMin2 && "two function val does not match!");
		double fVal = 999999.;
		double minPar[nPar], minParErr[nPar], Par[nPar];

		for (int i=0; i < nPar; ++i)
		{
			minPar[i] = func->GetParameter(i);
			minParErr[i] = func->GetParError(i);
		}

		//TRandom ran; // if I do this, I'll get the same set of random
		//numbers everytime. Which means I'll be shifting every
		//nominal val by some constant!
		for (int i=0; i<1000; ++i)
		{
			for (int npar =0 ; npar<nPar; ++npar)
			{
				const double dRan = ran.Uniform();
				//Par[npar] = minPar[npar]+(2.0*(ran.Uniform()-0.5)) * minParErr[npar];
				Par[npar] = minPar[npar]+(2.0*(dRan-0.5)) * minParErr[npar];
				//std::cout << blue << " par list = " << npar << "min = " << minPar[npar] << ", gen= " << Par[npar] << " ** minParErr = " << minParErr[npar] << clearatt << std::endl;
				//std::cout << blue << "npar, ran = " << npar << " - > " << dRan << clearatt << std::endl;
				final->SetParameters(Par);
			}
			//final->Print();
			fVal = final->Eval(x);
			if (fVal - fValMin<1) 
			{
				//std::cout << cyan << "Minval/val found = " << fValMin << ", " << fVal  << " (" << fVal - fValMin << ")"  << clearatt << std::endl;
				return fVal;
			}
		}
		*/
	} else	if (GetReweightSideband() == 3)
	{
		const float W = tf1_sideband_phoet->Eval(phoEt);
		hSidebandPhoEtWeights->Fill(phoEt, W);
		return W;

	} else if (GetReweightSideband() == 4)
	{
		const float W = tf1_sideband_phoet_error->Eval(phoEt);
		hSidebandPhoEtWeights->Fill(phoEt, W);
		return W;
	}

	
	return 1;

}
void ApplyJetCuts(JetList& jets)
{




}
void PhotonJets::GetHistBin(const TH1* hist_JetEtaWghts, const float jetEta, float &val, float &err)
{
	val =1, err=1;
	for (int i=1; i<=hist_JetEtaWghts->GetNbinsX(); ++i)
	{
		if (jetEta<hist_JetEtaWghts->GetBinLowEdge(i) 
			|| jetEta >= hist_JetEtaWghts->GetXaxis()->GetBinUpEdge(i)) 
			continue;
		val = hist_JetEtaWghts->GetBinContent(i);	
		err = hist_JetEtaWghts->GetBinError(i);	
	}

	if (val==1 && err ==1)
		std::cout << __FUNCTION__ <<": val and err are set to 1!" << std::endl;
}

bool PhotonJets::BadMet(const CommonVars& Vars, const PhotonList& phos, 
							const JetList& jets, const int iPhoIndex, 
							const int iJetIndex, TH1* phist, TH1* jhist) 
{
/*	std::cout << __FUNCTION__ << ": NOT COMPLETE YET! Works only on version >= V7 stuples." 
		 << " returning false!" << std::endl;
	return false;
*/
	//find the jet closest to met and see if it is withitin <0.4
	
	TVector2 tv2MetVec(Vars.met_MetX, Vars.met_MetY);
	TLorentzVector tlPhoVec(phos.pho_Px[iPhoIndex], phos.pho_Py[iPhoIndex], phos.pho_Pz[iPhoIndex], phos.pho_E[iPhoIndex]);
	const float dphi_pm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlPhoVec.Phi())));

	const float fDelPhiCut = 0.4;
	bool pbad = false;

	if (!pbad) phist->Fill(dphi_pm);
	
	bool bad = false;
	float dphi_jm_closest = 10.0;
	for (int i=0; i <jets.jet_num ; ++i)
	{
		TLorentzVector tlJetVec(jets.jet_Px[i], jets.jet_Py[i], jets.jet_Pz[i], jets.jet_E[i]);
		if (tlJetVec.Pt()<GetMinJetEt())
		{
			//std::cout << __FUNCTION__ << ": skipping jet Et<MinJetEt : " << tlJetVec.Pt() << " < " << GetMinJetEt() << std::endl;
			continue;
		}
		
		const float dphi_jm =fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(tv2MetVec.Phi())-TVector2::Phi_0_2pi(tlJetVec.Phi())));
		//std::cout << " ijet, dphi = " << i << ", " << dphi_jm  << "(TMath::Pi()-0.4 = " << TMath::Pi()-0.4 << ")"<< std::endl;
		if (dphi_jm <0.4 || dphi_jm>(TMath::Pi()-0.4))
		{
			//std::cout << red << "Failed MET CLEAN UP : jet Pt = " << tlJetVec.Pt() << clearatt << std::endl;
			bad = true;
		}
		if (dphi_jm < dphi_jm_closest) dphi_jm_closest = dphi_jm;
		
	}
	if (!bad) jhist->Fill(dphi_jm_closest);
	//return (bad || pbad);
	return (bad);  //looking at jet-met only for now
	
}

float PhotonJets::GetDATAMCNvtxWgt(const int Nvtx12)
{
	//this method only works when requiring nvtx>=1
	
	float w = 1;

/*	if (Nvtx12>=1 && Nvtx12<=hist_NvtxWghts->GetNbinsX())
	{
		w = hist_NvtxWghts->GetBinContent(Nvtx12);
	}
	*/
	if (Nvtx12>0 && Nvtx12<=vNvtxWeights.size()) w = vNvtxWeights.at(Nvtx12 -1);
		
	//std::cout << "nvtx, w = " << Nvtx12 << ", " << w << std::endl;
	return w;
}

bool PhotonJets::FailStripCuts(const Stuple& st)
{
	//apply the met cut here to speed things up
	//met is not calculated again in this so this should be fine
	//2-9-2010
	if (stuple.met_Met < fMinMet || stuple.met_Met > fMaxMet) return true;

	
	if (!bMcFlag)   		//_________________ 0=data, 1=mc
	{
		//______________________________________________________________ trigger
		if (  (st.tri_pho25iso != 1) && (st.tri_pho50 != 1)
				&& (st.tri_pho70 != 1) )
			return true;												
	}

	//______________________________________________________________ select the photon
	for (unsigned int i = 0; i < st.pho_num; ++i )
	{
		if (fabs(st.pho_DetEta[i])<fMaxPhotonEta &&
				(st.pho_Etc[i]>fMinPhotonEt && st.pho_Etc[i]<fMaxPhotonEt)) return false;
	}


		
	return true;
}

bool PhotonJets::SidebandPhoPassTightHadEmCut(const float& fEc, const float& fHadEm)
{
	if ( fHadEm< 0.125 || fHadEm < (0.055 + 0.00045 * fEc)) return true;
	else return false;
}

bool PhotonJets::SidebandPhoPassTightIsoCut(const float& fEtc, const float& fIso)
{
	if (fEtc < 20) {
		if (fIso > 0.1 * fEtc) return false;
	} else {
		if (fIso > (2.0+0.02 * (fEtc - 20.0)) ) return false;
	}
	return true;

}
bool PhotonJets::SidebandPhoPassTightTrkPtCut(const float& fEtc, const float& fTrkPt)
{
	if (fTrkPt> (1.0 + 0.005 * fEtc) ) return false;
	else return true;
}
bool PhotonJets::SidebandPhoPassTightTrkIsoCut(const float& fEtc, const float& fTrkIso)
{
	if (fTrkIso > (2.0 + 0.005 * fEtc)) return false;
	else return true;
}
void PhotonJets::CheckObjOverlap(const CommonVars& cvars, const PhotonList& phos, const ElectronList& eles, const JetList& jets)
{
	for (int ij=0; ij < jets.jet_num; ++ij)
	{
		TLorentzVector tlJetVec(0,0,0,0);
		tlJetVec.SetPxPyPzE(jets.jet_Px.at(ij), jets.jet_Py.at(ij), jets.jet_Pz.at(ij), jets.jet_E.at(ij));
		for (int ip=0; ip< phos.pho_num; ++ip)
		{
			TLorentzVector tlPhoVec(0,0,0,0);
			tlPhoVec.SetPxPyPzE(phos.pho_Px.at(ip), phos.pho_Py.at(ip), phos.pho_Pz.at(ip), phos.pho_E.at(ip));
			if (tlJetVec.DeltaR(tlPhoVec) < 0.4 && phos.pho_Ntight>0)
			{
				std::cout << blue << "OVER LAPPING PHOTON FOUND! ij, ip = " << ij << ", " << ip  << clearatt << std::endl;
				cvars.PrintHeader();
				phos.Dump();
				jets.Dump();
				jets.DumpRawJets();
				std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n" << std::endl;
			}
		}
		for (int ie=0; ie< eles.ele_num; ++ie)
		{
			TLorentzVector tlEleVec(0,0,0,0);
			tlEleVec.SetPxPyPzE(eles.ele_Px.at(ie), eles.ele_Py.at(ie), eles.ele_Pz.at(ie), eles.ele_E.at(ie));
			if (tlJetVec.DeltaR(tlEleVec) < 0.1)
			{
				std::cout << blue << "OVER LAPPING ELECTRON FOUND" << clearatt << std::endl;
			}
		}
	}
}

void PhotonJets::CheckJetTowers(const CommonVars& cvars, const JetList& jets)
{
	for (int ij=0; ij < jets.jet_num; ++ij)
	{
		if (jets.jet_Ntowers.at(ij) < 1 )
		{
			std::cout << blue << "JET WITH 0 TOWERS FOUND! ij = " << ij << clearatt << std::endl;
			cvars.PrintHeader();
			jets.Dump();
			jets.DumpRawJets();
			std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n" << std::endl;

		};
	}

}

int PhotonJets::FoundPartialJet(const CommonVars& cvars, const PhotonList& phos, const ElectronList& eles, const JetList& jets)
{
	for (int ij=0; ij < jets.jet_num; ++ij)
	{
		for (int ip=0; ip< phos.pho_num; ++ip)
		{
			if (phos.pho_matchJetIndex.at(ip) == jets.jet_Index.at(ij))
			{
				//std::cout << blue << "OVER LAPPING PHOTON FOUND! ij, ip = " << ij << ", " << ip  << clearatt << std::endl;
				//cvars.PrintHeader();
				++iRejEvts;
				return 1;
			}
		}
		for (int ie=0; ie< eles.ele_num; ++ie)
		{
			if (eles.ele_matchJetIndex.at(ie) == jets.jet_Index.at(ij))
			{
				//std::cout << blue << "OVER LAPPING ELECTRON FOUND" << clearatt << std::endl;
				++iRejEvts;
				return 1;
			}
		}
	}
	return 0;
}


bool PhotonJets::ExoticEvent(const Stuple& st)
{
	if (st.pho_Ntight > 0 && st.ele_num > 0)
	{
		for (int i=0; i< st.pho_num ; ++i)
		{
			if (st.pho_TightId[i] != 0) continue;
			if (st.pho_Etc[i] > 30.) continue;
			std::cout << "Run, Evt = " << st.evt_RunNumber << ", " << st.evt_EventNumber << std::endl;
			return true;
		}
	}
	return false;
}
