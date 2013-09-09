/*****************************************************************************
 * This is a newere version of MakeHistLogAndRatio.C. I have tried to use
 * globle variables to simply overall setting in the code.
 * I tested results of this and previous version using g30 and nominal
 * sideband method, and they are a match. I'll try to clean this up further
 * and declare file names as global varibales too. 05-20-2010
 *****************************************************************************/
/*{{{*/
/********************************************************************
 *  $Id: MakeHistLogAndRatioV2.C,v 1.13 2012/10/29 19:41:08 samantha Exp $
 *  $Log: MakeHistLogAndRatioV2.C,v $
 *  Revision 1.13  2012/10/29 19:41:08  samantha
 *  Final commit.
 *
 *  Revision 1.12  2011/04/21 17:13:53  samantha
 *  This is the version I used to make the blessing plot and for my dissertation.
 *  MAJOR CHANGES: 1. I have added lots of vectors to get the breakdown of different
 *  systeatic errors. Also to calculate the final event count with total errors.
 *  2. DiPho MC lum is added to the global lum variable list and is used in the
 *  error calculations.
 *  3. KS test results are added too. Though I never showed it to anyone, I have put
 *  the plots in an elog entry.
 *
 *  Minor Changes: Some cosmetic and spacing changes and added more comments
 *  for different parts of the code.
 *
 *  Revision 1.11  2010/12/01 17:49:47  samantha
 *  This is the version I used to make the blessing plots.
 *  Now the stopry is simple. I have only two methods.
 *  Method A: for all plots which uses nominal sideband and photon MC, even for Njet
 *  plot.
 *  Method B: for all plots which uses weighted sideband.
 *
 *  MAJOR CHANGES: I am using inclusive g30 fake photon fraction for all Method A
 *  plots. Ray said there is no justification for changing this for different
 *  samples. Especially there is an 10% increase in photon fraction in MET>20 sample
 *  which has no real physics explanation.
 *  DELETED: MakeJetEtaWeightsForSideband() method which is not used.
 *  ADDED: some stuff to get the final events count with systematics. This will be
 *  further broken down into individual backgrounds in the next version.
 *
 *  Revision 1.10  2010/11/13 22:54:48  samantha
 *  ADDED: g30jetsmet30 plots as jay requested. See elog:1748.
 *  Put back the NJET making method as before. I edited this to get njet
 *  plot using pho mc +qcd.
 *
 *  Revision 1.8  2010/10/21 03:29:37  samantha
 *  Previous Comitted revision and this is identical. This is a force
 *  commit because I could not add comments to previous commit.
 *
 *  Major Changes:
 *  1. This version has correct Lum up to p25 which is 4.8fb-1. I forgot to subtract
 *  the first 400pb-1.
 *  2. Has W+jets, correction for MET in W MC samples. Need to turn on the switch to
 *  use.
 *  3. Method B gets all the errors just like Method A, and the difference in
 *  background predictions. (Also reweighting error). But I need to run it 3 times
 *  to get Mtd B correctly. First run on Mtb B without A-B background errors. Then
 *  switch to Mtd A and make the plots using Mtd A-B errors. Then with'rerun' switch
 *  on run on Mtd B to make the plots with A-B errors.
 *  4. Added two more samples. g30jetsMet20 and g30jetsMet50.
 *  5. Also when I run out of statics in photon MC, as in g30JetsMet50 sample, and
 *  total EWK+DIPHO sum is more than data integral forces me to scale pho mc and QCD
 *  to zero.
 *  5. I named all the histograms so it's easy to id them.
 *  6. halo and cosmic are collectively called non-collision.
 *  7. Titles and labels have been modified according to Jay and Ray's suggestions.
 *  8. I have some stuff to manually smoth out some hist's bins with abnormally
 *  large systematic errors.
 *
 *  I need to check if the photon MC is correclty scaled when ISR/FSR/PDF/Q2 is
 *  calculated in mtd B. Because was scaling them to 0 ealier to be safe in mtd B.
 *  Now we are considering keeping both plots, this need to be checked.
 *
 *  Revision 1.7  2010/10/21 02:50:43  samantha
 *  *** empty log message ***
 *
 *  Revision 1.6  2010/09/17 17:57:40  samantha
 *  This version has some modification I made to make a nice MET plot
 *  for TTU talk.
 *  CORRECTIONS for MET:
 *  1. I have applied w+jets MET correction for MET in Wmc.
 *  2. Has made changes to the bin sizes in mass and met plots.
 *
 *  Revision 1.5  2010/07/06 19:48:43  samantha
 *  BUG FIX: Some EWK sample values used in luminosity were wrong. At some point
 *  some of the EWK samples were not accessible. So I dropped them from the
 *  calculation. But now they are availbale again. And I have been using them
 *  without updating the total number if event in each sample. So I fixed that.
 *
 *  ADDED:
 *  1. CorrectEWKMC() for EWK MC sample correction for MET plot. For now this method is temp fix
 *  and is npt valied for any other distributions. See elog:1681
 *  2. SubtractDiPhofromPhoMC() to remove the di-pho contribution in pho+jets MC
 *  sample to remove double counting.
 *
 *  MODIFIED: MET distribution binning.
 *
 *  Revision 1.4  2010/06/22 00:10:37  samantha
 *  Modified to include generate the method B errors. For this I need to run it
 *  under WEIGHTED_SIDEBAND once and save background prediction to a file Then
 *  run with NOMINAL_SIDEBAND. I do not draw the ErrorBreackdown any more in the
 *  WEIGHTED _SIDEBAND method.Also prints the mtdBerr error info.
 *
 *  Revision 1.3  2010/06/13 16:15:10  samantha
 *  MAJOR CHANGES: xpoints and widdth for variable bins are made global now. this
 *  reduces the clutter in the number of parameters to many functions. cleaned up
 *  all debug methods.
 *  So far things seems to hold together in this new version.
 *  Now this needs to be run in 2 modes:
 *  1. is to derive unceratinties from method 2 to be used in method 1. I need to
 *  save the background prediction from method 2 in file.
 *  2. then recall them in method 1 to derive uncertainties by comparing to method 1
 *  prediction.
 *
 *  Revision 1.2  2010/05/23 19:53:02  samantha
 *  MAJOR CHANGES:
 *  1. Lot of cleaninig up done. Sorted estimates and fake frations
 *  setting for each sample into functions etc.
 *  2. Testing on METHOD B predictions. Ran into trouble when calculating systematics
 *  for pho fraction. I removed the qcdjet_100 hist altogether and tried to use
 *  qcdjet for method B pointing to weighted hist. But to calcualate pho.frac syst I
 *  need nominal sideband too! But for now I disabled all syst calcualtions for
 *  method B as I am going to use this to derive errors for method A. So i'll plot
 *  the log/ratio plot for method B without syst errors.
 *
 *  Revision 1.1  2010/05/20 15:07:43  samantha
 *  This is a newere version of MakeHistLogAndRatio.C. I have tried to use
 *  globle variables to simply overall setting in the code.
 *  I tested results of this and previous version using g30 and nominal
 *  sideband method, and they are a match. I'll try to clean this up further
 *  and declare file names as global varibales too. 05-20-2010
 *
 ********************************************************************/
/*}}}*/
 
#include <iostream>
#include <sstream>
#include <string>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <TF1.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLine.h>
#include <cmath>
#include <memory>
#include <iomanip>
#include "CommonTools.hh"
#include "IOColors.hh"
#include "Rtypes.h" //ROOT color keyword. this is included in ROOT headers though.
#include "TMatrixTSym.h"
//#include "TFitResultPtr.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TProfile.h"

using namespace std;


/********************** GLOBAL VARIABLES ***************************/
//const std::string sPRINTFOLDER("/home/samantha/TMP/eps/g30_nominal/");
std::string sPRINTFOLDER("/home/samantha/TMP/eps/");
const std::string sPRINTFORMAT("eps");			//generate eps/pdf/gif files as output
const int iPRL_MODE = 1;				//Enable/disable addition info on the canvas
const bool iDRAWERROR_BREAKDOWN = 0;	//on/off error break down canvas
const bool bDO_VAR_BINNING = 1;
int iREBIN = 5;		//DEFINES THE CONST REBIN FOR EACH PLOTS
const bool bCorrectWMC_Met = true;  //apply W+jets met corrections to W mc samples
int iHAVE_MTDB_ERRORS = 1;  //some samples have no mtd b errors

//for debugging and printing additiona info
const static int iFUNCTION_NAMES = 5;
const static int iDEBUG = 0;

//differnt path
const std::string sDATA_DIR("Hist/SIGNAL/");
const std::string sHALO_DIR("Hist/HALO/");
const std::string sCOSMIC_DIR("Hist/COSMIC/");
const std::string sSIDEBAND_DIR("Hist/SIDEBAND/");
const std::string sSIDEBANDSYST_HADEM_DIR("Hist/SIDEBANDSYST_HADEM/");
const std::string sSIDEBANDSYST_ISO_DIR("Hist/SIDEBANDSYST_ISO/");
const std::string sSIDEBANDSYST_TRKPT_DIR("Hist/SIDEBANDSYST_TRKPT/");
const std::string sSIDEBANDSYST_TRKISO_DIR("Hist/SIDEBANDSYST_TRKISO/");
const std::string sSIDEBAND_COSMIC_DIR("Hist/SIDEBAND_COSMIC/");
const std::string sSIDEBAND_HALO_DIR("Hist/SIDEBAND_HALO/");
const std::string sMCCENTRAL_DIR("Hist/CENTRAL/");
const std::string sMCEMJESUP_DIR("Hist/EMJESUP/");
const std::string sMCEMJESDOWN_DIR("Hist/EMJESDOWN/");
const std::string sMCSIDEBAND_DIR("Hist/SIDEBAND/");

//these will be set for each hist based in njet
std::string sDATA_HIST_PATH("");  
std::string sHALO_HIST_PATH("");  
std::string sCOSMIC_HIST_PATH("");  
std::string sSIDEBAND_HIST_PATH("");  
std::string sSIDEBANDSYST_HADEM_HIST_PATH("");  
std::string sSIDEBANDSYST_ISO_HIST_PATH("");
std::string sSIDEBANDSYST_TRKPT_HIST_PATH("");
std::string sSIDEBANDSYST_TRKISO_HIST_PATH("");
std::string sSIDEBAND_COSMIC_HIST_PATH("");  
std::string sSIDEBAND_HALO_HIST_PATH("");  
std::string sMCCENTRAL_HIST_PATH("");  
std::string sMCEMJESUP_HIST_PATH("");  
std::string sMCEMJESDOWN_HIST_PATH("");  
std::string sMCSIDEBAND_HIST_PATH("");  

//all the required files need to run the job
//these are generic names. create softlinks with 
//these names to point to the appropriate files
const std::string sDATA_FILE_NOMINAL("PhoJets_data.root");
const std::string sDATA_FILE_WGTSIDEBAND("PhoJets_data_wgtdSideband.root"); //in rewegihted methods this will be used as the file_DATA_NOMINAL
const std::string sDATA_FILE_WGTSIDEBAND_ERR("PhoJets_data_wgtdSideband_err.root");
const std::string sMCPHO_FILE("PhoJets_phomc.root");
const std::string sMCZEE_FILE("PhoJets_zeemc.root");
const std::string sMCZMM_FILE("PhoJets_zmmmc.root");
const std::string sMCZTT_FILE("PhoJets_zttmc.root");
const std::string sMCWEN_FILE("PhoJets_wenmc.root");
const std::string sMCWMN_FILE("PhoJets_wmnmc.root");
const std::string sMCWTN_FILE("PhoJets_wtnmc.root");
const std::string sMCDIPHO_FILE("PhoJets_diphomc.root");

TFile *file_DATA_NOMINAL = 0;
TFile *file_DATA_WGTSIDEBAND = 0;
TFile *file_DATA_WGTSIDEBAND_ERR = 0; //in reweighted method, file_DATA_NOMINAL point to the sideband wighted with nominal weights.
TFile *file_MCPHO = 0;
TFile *file_ZEEMC = 0;
TFile *file_ZMMMC = 0;
TFile *file_ZTTMC = 0;
TFile *file_WENMC = 0;
TFile *file_WMNMC = 0;
TFile *file_WTNMC = 0;
TFile *file_DIPHOMC = 0;


//different samples
const static int iG30      = 0; //default gamma30+jets case
const static int iG40      = 1; //default gamma40+jets case
const static int iG30MET20 = 2; //gamma30+jets+met20 case
const static int iG30MET25 = 3; //gamma30+jets+met25 case
const static int iG30MET30 = 4; //gamma30+jets+met20 case
const static int iG30MET40 = 5; //gamma30+jets+met40 case
const static int iG30EXCL1J= 6; // gamma30+1==jet case
const static int iG30MET50 = 7; //gamma30+jets+met40 case

std::string sSAMPLE_NAME("UNKNOWN"); //this will be set at runtime
const static std::string sG30 = "G30Jets";
const static std::string sG40 = "G40";
const static std::string sG30MET20 = "G30JetsMet20";
const static std::string sG30MET25 = "G30Met25";
const static std::string sG30MET30 = "G30Met30";
const static std::string sG30MET40 = "G30Met40";
const static std::string sG30EXCL1J = "G30EXCL1J";
const static std::string sG30MET50 = "G30MET50";

//different sidebands
const static int iNOMINAL_SIDEBAND = 0;
const static int iWEIGHTED_SIDEBAND = 1;
int iMTDB_RERUN = 0;   //need this to include A-B diff in mtd B.
const static std::string sNOMINAL_SIDEBAND("NOMINAL_SIDEBAND");
const static std::string sWEIGHTED_SIDEBAND("WEIGHTED_SIDEBAND");

//fake photon fraction for different sample
float fFAKE_PHO_FRAC_DEF = 0;   //default val, will be assigned based on the sample selection dynamically
float fFAKE_PHO_FRAC_PSIG = 0;   //+sigma val, will be assigned based on the sample selection dynamically
const static float fG30_FAKE_PHO_FRAC_1JET = 0.303;  			//for g30+>=1 jets
const static float fG30_FAKE_PHO_FRAC_2JET = 0.273;  			//for g30+>=2 jets   ::elog#1648
const static float fG40_FAKE_PHO_FRAC = 0.21;  				//for g40+>=1 jets
const static float fG30MET20_FAKE_PHO_FRAC_1JET = 0.21; 	//for g30+>=1jets+met20
const static float fG30MET20_FAKE_PHO_FRAC_2JET = 0.19; 	//for g30+>=2jets+met20
const static float fG30MET25_FAKE_PHO_FRAC_1JET = 0.15; 	//for g30+>=1jets+met25
const static float fG30MET25_FAKE_PHO_FRAC_2JET = 0.14; 	//for g30+>=2jets+met25
const static float fG30MET30_FAKE_PHO_FRAC_1JET = 0.235; 	//for g30+>=1jets+met30
const static float fG30MET30_FAKE_PHO_FRAC_2JET = 0.116; 	//for g30+>=2jets+met30
const static float fG30MET40_FAKE_PHO_FRAC_1JET = 0.20; 	//for g30+>=1jets+met40
const static float fG30MET40_FAKE_PHO_FRAC_2JET = 0.13; 	//for g30+>=2jets+met40
const static float fG30EXCL1J_FAKE_PHO_FRAC = 0.313;  			//for g30+==1 jets
const static float fG30MET50_FAKE_PHO_FRAC_1JET = 0.324; 	//for g30+>=1jets+met50
const static float fG30MET50_FAKE_PHO_FRAC_2JET = 0.269; 	//for g30+>=2jets+met50

//systematic error for all sample are taken to this assuming the change
//of this value is small - 02-17-2010
const static float fFAKE_PHO_FRAC_SIGMA = 0.068;

//cosmic background estimates for different samples
//const static float fG30_COSMIC_EST_1JET = 110;  //original APS/ICHEP estimates
//const static float fG30_COSMIC_EST_2JET = 7;
const static float fG30_COSMIC_EST_1JET = 131; //for p26
const static float fG30_COSMIC_EST_2JET = 7;
const static float fG40_COSMIC_EST_1JET = 60;
const static float fG40_COSMIC_EST_2JET = 5;
const static float fG30MET20_COSMIC_EST_1JET = 124;
const static float fG30MET20_COSMIC_EST_2JET = 8;
const static float fG30MET25_COSMIC_EST_1JET = 82;
const static float fG30MET25_COSMIC_EST_2JET = 6;
const static float fG30MET30_COSMIC_EST_1JET = 108;
const static float fG30MET30_COSMIC_EST_2JET = 8;
const static float fG30MET40_COSMIC_EST_1JET = 63;
const static float fG30MET40_COSMIC_EST_2JET = 5;
const static float fG30EXCL1J_COSMIC_EST_1JET = 83;
const static float fG30MET50_COSMIC_EST_1JET = 65;
const static float fG30MET50_COSMIC_EST_2JET = 5;
//estimates for sideband (to be subtracted.
const static float fG30_SIDEBAND_COSMIC_EST_1JET = 32; //for p26
const static float fG30_SIDEBAND_COSMIC_EST_2JET = 2;  //for p26

const static float fG40_SIDEBAND_COSMIC_EST_1JET = 15; //new estimates 03-02-2010
const static float fG40_SIDEBAND_COSMIC_EST_2JET = 2;
const static float fG30MET20_SIDEBAND_COSMIC_EST_1JET = 30;
const static float fG30MET20_SIDEBAND_COSMIC_EST_2JET = 2;
const static float fG30MET25_SIDEBAND_COSMIC_EST_1JET = 21; //new estimates 02-26-2010
const static float fG30MET25_SIDEBAND_COSMIC_EST_2JET = 2;
const static float fG30EXCL1J_SIDEBAND_COSMIC_EST_1JET = 20;
const static float fG30MET30_SIDEBAND_COSMIC_EST_1JET = 28;
const static float fG30MET30_SIDEBAND_COSMIC_EST_2JET = 2;
const static float fG30MET40_SIDEBAND_COSMIC_EST_1JET = 18; //new estimates 03-08-2010
const static float fG30MET40_SIDEBAND_COSMIC_EST_2JET = 2;
const static float fG30MET50_SIDEBAND_COSMIC_EST_1JET = 17; //new estimates 03-08-2010
const static float fG30MET50_SIDEBAND_COSMIC_EST_2JET = 1;


//halo background estimates for different samples
//const static float fG30_HALO_EST_1JET = 9;
//const static float fG30_HALO_EST_2JET = 1;
const static float fG30_HALO_EST_1JET = 1; //<1 event with delphi(jet-met) cut
const static float fG30_HALO_EST_2JET = 1;
const static float fG40_HALO_EST_1JET = 9;   //needs to be updated, using g30 values
const static float fG40_HALO_EST_2JET = 1;   //needs to be updated, using g30 values
const static float fG30MET20_HALO_EST_1JET = 1;
const static float fG30MET20_HALO_EST_2JET = 1;
const static float fG30MET25_HALO_EST_1JET = 9;   //needs to be updated, using g30 values
const static float fG30MET25_HALO_EST_2JET = 1;   //needs to be updated, using g30 values
const static float fG30MET30_HALO_EST_1JET = 1;
const static float fG30MET30_HALO_EST_2JET = 1;
const static float fG30MET40_HALO_EST_1JET = 9;   //needs to be updated, using g30 values
const static float fG30MET40_HALO_EST_2JET = 1;   //needs to be updated, using g30 values
const static float fG30EXCL1J_HALO_EST_1JET = 9;
const static float fG30MET50_HALO_EST_1JET = 1;
const static float fG30MET50_HALO_EST_2JET = 1;

const int iTITLE_FONT = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable) 
const float iTITLE_FONT_SIZE = 0.04;
const float iTITLE_TEXT_ALIGN = 10*1+3; //10*horizontal+vertical  :horizontal: 1=left 2=center 3= right verticle 1=bottom 2=center 3=top
const int iWHITE = 10;		//ROOT color whiTE
//const int iAXIS_LABEL_FONT = 102;
const int iAXIS_LABEL_FONT = iTITLE_FONT;

//const static float fLUM_DATA    = 5360.11;	//pb wrong this incluedes p26 which goodrun31 does not cover!
const static float fLUM_DATA    = 4832.6;	//pb up to and including p25 from goodrun v31.

//NLO crosssections correction
const float kFAC_MC = 1.4;

//DIPHOMC crosssection after filter(12.7%) = 92.3 pb
// 727 pb times a filter ratio of 12.7% = 92.3pb
//total number of events in the sample (before goodrun) 11,762,656
const static float fLUM_DIPHOMC = 127439.39;  // L= n/ (e * cs) = n/cs' = 11,762,656/92.3
const float fDIPHOMC_SCALE = (fLUM_DATA/fLUM_DIPHOMC) * kFAC_MC;	

// SF = DATA_LUM * ( 1 / EWK_LUM ) * KFAC
//    = DATA_LUM * ( 1 / (TOT EVTS PROCESSED/CROSS SECTION) ) * KFAC
const static float fCS_ZEE_INC 		= 355;	//pb-1
const static float fCS_ZMM_INC 		= 355;	//pb-1
const static float fCS_ZTT_INC_P0 	= 355;	//pb-1
const static float fCS_ZTT_INC_P1_7 = 238;	//pb-1
const static float fCS_ZTT_INC_P8_13= 355;	//pb-1
const static float fCS_WEN_INC 		= 1960;	//pb-1
const static float fCS_WMN_INC 		= 1960;	//pb-1
const static float fCS_WTN_INC 		= 1960;	//pb-1

//these are the new numbers which includes all the EWK dataset.
//some of them was not accessible at some point and I excluded
//them (above numbers are for those exclude datasets).
//But they are available now.
const static float fNTOT_ZEE 		= 23323254;
const static float fNTOT_ZMM 		= 23279458 ;
const static float fNTOT_ZTT_P0 	= 9438815;
const static float fNTOT_ZTT_P1_7= 6916050;
const static float fNTOT_ZTT_P8_13= 6575236;
const static float fNTOT_WEN 		= 34134323;
const static float fNTOT_WMN 		= 23252069;
const static float fNTOT_WTN 		= 29471407;


const float fLUM_ZEE = fNTOT_ZEE / (fCS_ZEE_INC * kFAC_MC);
const float fLUM_ZMM = fNTOT_ZMM / (fCS_ZMM_INC * kFAC_MC);
//I am taking the average Lum for this dataset for now 
const float fLUM_ZTT = (fNTOT_ZTT_P0 / (fCS_ZTT_INC_P0 * kFAC_MC) + fNTOT_ZTT_P1_7 / (fCS_ZTT_INC_P1_7 * kFAC_MC) + fNTOT_ZTT_P8_13/ (fCS_ZTT_INC_P8_13 * kFAC_MC))/3.;
const float fLUM_WEN = fNTOT_WEN / (fCS_WEN_INC * kFAC_MC);
const float fLUM_WMN = fNTOT_WMN / (fCS_WMN_INC * kFAC_MC);
const float fLUM_WTN = fNTOT_WTN / (fCS_WTN_INC * kFAC_MC);

const float fZEE_SCALE = fLUM_DATA / fLUM_ZEE;
const float fZMM_SCALE = fLUM_DATA / fLUM_ZMM;
const float fZTT_SCALE = fLUM_DATA / fLUM_ZTT;
const float fWEN_SCALE = fLUM_DATA / fLUM_WEN;
const float fWMN_SCALE = fLUM_DATA / fLUM_WMN;
const float fWTN_SCALE = fLUM_DATA / fLUM_WTN;





//globle error storage
std::string sERRORS("");
std::string sWARNINGS("");

//hist bin widths
float fXMIN=0, fXPOINT1=0, fXPOINT2=0, fXPOINT3=0, fXPOINT4=0; 
float fWIDTH1=0, fWIDTH2=0, fWIDTH3=0, fWIDTH4=0;


//resets global histpaths
void ClearHistPaths()
{
	sDATA_HIST_PATH = "";  
	sHALO_HIST_PATH = "";  
	sCOSMIC_HIST_PATH = "";  
	sSIDEBAND_HIST_PATH = "";  
	sSIDEBANDSYST_HADEM_HIST_PATH = "";  
	sSIDEBANDSYST_ISO_HIST_PATH = "";
	sSIDEBANDSYST_TRKPT_HIST_PATH = "";
	sSIDEBANDSYST_TRKISO_HIST_PATH = "";
	sSIDEBAND_COSMIC_HIST_PATH = "";  
	sSIDEBAND_HALO_HIST_PATH = "";  
	sMCCENTRAL_HIST_PATH = "";  
	sMCEMJESUP_HIST_PATH = "";  
	sMCEMJESDOWN_HIST_PATH = "";  
	sMCSIDEBAND_HIST_PATH = "";  
}

/*************** END OF GLOBAL VARIABLES ***************************/

/*************** PROTOTYPES ***************************************/
void SubtractEWKfromQCD(TH1* qcdhist, const bool debug=false);
void SubtractCosmicFromQCD(TH1* qcdhist, const int iSample, const int njet);
void SubtractDiPhofromPhoMC(TH1* mcpho, const bool debug=false);
void GetSidebandReweighedError(std::vector<float>& ErrVec, const TH1* qcdhist,
						const float qcdscale, const int iSample, const int njets,
						bool debug=false);
void CheckUnderflowOverflow(const TH1* hist);

void CheckUnderflowOverflow(const TH1* hist)
{
	assert (hist != NULL && "CheckUnderflowOverflow:: Given Hist is null!");
	if (hist->GetBinContent(0) != 0)
	{
		std::cout << red << "This hist " << hist->GetName() << " has underflow of = " << hist->GetBinContent(0)  << clearatt << std::endl;
	}
	if (hist->GetBinContent(hist->GetNbinsX()+1) != 0)
	{
		std::cout << red << "This hist " << hist->GetName() << " has overflow of = " << hist->GetBinContent(hist->GetNbinsX()+1)  << clearatt << std::endl;
	}
}
void CheckIntegral(const TH1* hist)
{
	assert (hist != NULL && "CheckIntegral:: Given Hist is null!");
	double ww = 0;
	for (int bin=1; bin <= hist->GetNbinsX(); bin++)
	{
		ww += hist->GetBinContent(bin) * hist->GetBinWidth(bin);
	}
	std::cout << red << "This hist " << hist->GetName() << " integral with width from ROOT/calculated= " << hist->Integral("width") << "/ " << ww  << clearatt << std::endl;
}



/* Correct EWK MC to match EWK data elog:
 * this should be only used for MET plot for now.
 * I need to run over the whole EWK MC sample
 * with these weights. Then there is no need for this method. 06-28-2010
 */
void CorrectEWKMC(TH1* hist, const int met_type)
{
	//met_type == 0 for QCD region (Zll)
	//met_type == 1 for HiMET region (Wll)
	
	assert (hist != NULL && "CorrectEWKMC:: hist ptr given is null!");
	assert (( met_type == 0 || met_type == 1) && "CorrectEWKMC:: unknown MET type! must be 0 or 1!");
	TF1 *fw= 0;
	
	//fit for Zll
	 //p0                        =     1.05609         +/-     0.0141337
	// p1                        =     -0.00816382     +/-     0.00169421
	 if (met_type == 0) fw = new TF1("fw","1.05609 - 0.00816382 * x");
	 else fw = new TF1("fw","1");

	std::cout << red << "bin\tval*\tw = newval" << std::endl;
	for (int bin=1; bin <= hist->GetNbinsX(); ++bin)
	{
		const float w = fw->Eval(hist->GetBinCenter(bin));
		std::cout << bin << "\t" << hist->GetBinContent(bin) <<", w= " << w; 
		hist->SetBinContent(bin, hist->GetBinContent(bin) * w);
		hist->SetBinError(bin, hist->GetBinError(bin) * w);
		std::cout << "->"<< hist->GetBinContent(bin)	<< std::endl;
	}
	std::cout << clearatt << std::endl;

}

void CorrectWMC(TH1* hist)
{
//this is W+JETS with met clean up MET weights from elog:1720
//
//bin loEdge val err
//5, 20, 0.989131, 0.0343441
//6, 25, 1.02107, 0.0315942
//7, 30, 1.03474, 0.0361743
//8, 35, 0.98105, 0.0461473
//9, 40, 0.915306, 0.0613873
//10, 45, 0.948611, 0.098749
//11, 50, 0.996498, 0.193267
//12, 55, 0.778355, 0.40225
//SCALE = DATA/MC = 0.00153665

	assert (hist != NULL && "CorrectWMC:: hist ptr given is null!");

	std::cout << green << "bin\tval*\tw = newval" << std::endl;
	for (int bin=1; bin <= hist->GetNbinsX(); ++bin)
	{
		const float x = hist->GetBinCenter(bin);
		if (x<20. || x>55.) continue; 
		float w = 1;

		if (x>=20 && x<25) w = 0.989131;
		else if (x>=25 && x<30) w = 1.02107;
		else if (x>=30 && x<35) w = 1.03474;
		else if (x>=35 && x<40) w = 0.98105;
		else if (x>=40 && x<45) w = 0.915306;
		else if (x>=45 && x<50) w = 0.948611;
		else if (x>=50 && x<55) w = 0.996498;
		else if (x>=55 && x<60) w = 0.778355;

		hist->SetBinContent(bin, hist->GetBinContent(bin) * w);
		hist->SetBinError(bin, hist->GetBinError(bin) * w);
		std::cout << "->"<< hist->GetBinContent(bin)	<< std::endl;
		
	}
	std::cout << clearatt << std::endl;

}






std::string WhichSample(const int iSample)
{
	std::string sample("UNKNOWN");

	if (iSample == iG30) sample = sG30;
	else if (iSample == iG40) sample = sG40;
	else if (iSample == iG30MET20) sample = sG30MET20;
	else if (iSample == iG30MET25) sample = sG30MET25;
	else if (iSample == iG30MET30) sample = sG30MET30;
	else if (iSample == iG30MET40) sample = sG30MET40;
	else if (iSample == iG30MET50) sample = sG30MET50;
	else if (iSample == iG30EXCL1J) sample = sG30EXCL1J;

	return sample;
}
std::string WhichSideband(const int iSideband)
{
	std::string sideband("UNKNOWN");
	if (iSideband == iNOMINAL_SIDEBAND) sideband = sNOMINAL_SIDEBAND;
	else if (iSideband == iWEIGHTED_SIDEBAND) sideband = sWEIGHTED_SIDEBAND;

	return sideband;
}
void GetHaloCosmicEstimates(const int isample, const int jets, float &haloEst, float &cosmicEst)
{
	if (jets ==1 ) {
		if (isample == iG30)
		{
			haloEst = fG30_HALO_EST_1JET;				//these estimates are for the full dataset
			cosmicEst = fG30_COSMIC_EST_1JET;
		} else if (isample == iG40)
		{
			haloEst = fG40_HALO_EST_1JET;
			cosmicEst = fG40_COSMIC_EST_1JET;
		} else if (isample == iG30MET20)
		{
			haloEst = fG30MET20_HALO_EST_1JET;
			cosmicEst = fG30MET20_COSMIC_EST_1JET;
		} else if (isample == iG30MET25)
		{
			haloEst = fG30MET25_HALO_EST_1JET;
			cosmicEst = fG30MET25_COSMIC_EST_1JET;
		} else if (isample == iG30MET30)
		{
			haloEst = fG30MET30_HALO_EST_1JET;
			cosmicEst = fG30MET30_COSMIC_EST_1JET;
		} else if (isample == iG30MET40)
		{
			haloEst = fG30MET40_HALO_EST_1JET;
			cosmicEst = fG30MET40_COSMIC_EST_1JET;
		} else if (isample == iG30MET50)
		{
			haloEst = fG30MET50_HALO_EST_1JET;
			cosmicEst = fG30MET50_COSMIC_EST_1JET;
		} else if (isample == iG30EXCL1J)
		{
			haloEst = fG30EXCL1J_HALO_EST_1JET;				//these estimates are for the full dataset
			cosmicEst = fG30EXCL1J_COSMIC_EST_1JET;
		}
		
	} else if (jets ==2 ) {
		if (isample == iG30)
		{
			haloEst = fG30_HALO_EST_2JET;
			cosmicEst = fG30_COSMIC_EST_2JET ;
		} else if (isample == iG40)
		{
			haloEst = fG40_HALO_EST_2JET;
			cosmicEst = fG40_COSMIC_EST_2JET;
		} else if (isample == iG30MET20)
		{
			haloEst = fG30MET20_HALO_EST_2JET;
			cosmicEst = fG30MET20_COSMIC_EST_2JET;
		} else if (isample == iG30MET25)
		{
			haloEst = fG30MET25_HALO_EST_2JET;
			cosmicEst = fG30MET25_COSMIC_EST_2JET;
		} else if (isample == iG30MET30)
		{
			haloEst = fG30MET30_HALO_EST_2JET;
			cosmicEst = fG30MET30_COSMIC_EST_2JET;
		} else if (isample == iG30MET40)
		{
			haloEst = fG30MET40_HALO_EST_2JET;
			cosmicEst = fG30MET40_COSMIC_EST_2JET;
		} else if (isample == iG30MET50)
		{
			haloEst = fG30MET50_HALO_EST_2JET;
			cosmicEst = fG30MET50_COSMIC_EST_2JET;
		} else if (isample == iG30EXCL1J)
		{
			std::cout << red << " ERROR! no cosmic/halo estimate for excl 1jet case!" << clearatt << std::endl;
			assert (false);
		}

	} else {
		std::cout << __FILE__ <<"::"<<__LINE__<<":: Invalid number of jets !" << std::endl;
		exit (1);
	}
}

void SetFakePhotonFraction(const int isample, const int jets)
{

	fFAKE_PHO_FRAC_DEF = fG30_FAKE_PHO_FRAC_1JET;  	//for pho Et>30GeV+jets
   fFAKE_PHO_FRAC_PSIG = fFAKE_PHO_FRAC_DEF + fFAKE_PHO_FRAC_SIGMA;
	return;

	if (isample == iG30)
	{
		if (jets == 1) fFAKE_PHO_FRAC_DEF = fG30_FAKE_PHO_FRAC_1JET;  	//for pho Et>30GeV+jets
		else if (jets == 2) fFAKE_PHO_FRAC_DEF = fG30_FAKE_PHO_FRAC_2JET;  	//for pho Et>30GeV+jets

	}
	else if (isample == iG40) 	fFAKE_PHO_FRAC_DEF = fG40_FAKE_PHO_FRAC; 	//for pho Et>40GeV+jets
	else if (isample == iG30MET20)
	{
		if (jets == 1) 		fFAKE_PHO_FRAC_DEF = fG30MET20_FAKE_PHO_FRAC_1JET;
		else if (jets == 2) 	fFAKE_PHO_FRAC_DEF = fG30MET20_FAKE_PHO_FRAC_2JET;
	} else if (isample == iG30MET25)
	{
		if (jets == 1) 		fFAKE_PHO_FRAC_DEF = fG30MET25_FAKE_PHO_FRAC_1JET;
		else if (jets == 2) 	fFAKE_PHO_FRAC_DEF = fG30MET25_FAKE_PHO_FRAC_2JET;
	} else if (isample == iG30MET30)
	{
		if (jets == 1) 		fFAKE_PHO_FRAC_DEF = fG30MET30_FAKE_PHO_FRAC_1JET;
		else if (jets == 2) 	fFAKE_PHO_FRAC_DEF = fG30MET30_FAKE_PHO_FRAC_2JET;
	} else if (isample == iG30MET40)
	{
		if (jets == 1) 		fFAKE_PHO_FRAC_DEF = fG30MET40_FAKE_PHO_FRAC_1JET;
		else if (jets == 2) 	fFAKE_PHO_FRAC_DEF = fG30MET40_FAKE_PHO_FRAC_2JET;
	} else if (isample == iG30MET50)
	{
		if (jets == 1) 		fFAKE_PHO_FRAC_DEF = fG30MET50_FAKE_PHO_FRAC_1JET;
		else if (jets == 2) 	fFAKE_PHO_FRAC_DEF = fG30MET50_FAKE_PHO_FRAC_2JET;
	} else if (isample == iG30EXCL1J) 
	{
		if (jets == 1) fFAKE_PHO_FRAC_DEF = fG30EXCL1J_FAKE_PHO_FRAC;
		else if (jets == 2)
		{
			std::cout << red << __FUNCTION__ << __LINE__ << ": ERROR! no fake photon fraction for excl 1jet sample!" << clearatt << std::endl;
			assert (false);
		}
	} else
	{
		std::cout << red << __FUNCTION__ << __LINE__ << ": UNKNOWN sample selction! invalid fake photon fraction!" << clearatt << std::endl;
		assert (false);
	}

  fFAKE_PHO_FRAC_PSIG = fFAKE_PHO_FRAC_DEF + fFAKE_PHO_FRAC_SIGMA;

  assert ( (fFAKE_PHO_FRAC_DEF > 0. && fFAKE_PHO_FRAC_PSIG >0 ) && "SetFakePhotonFraction:: fFAKE_PHO_FRAC_DEF is 0!");
  assert ( (fFAKE_PHO_FRAC_DEF != fFAKE_PHO_FRAC_PSIG) && "SetFakePhotonFraction:: fFAKE_PHO_FRAC_DEF cannot be equal to fFAKE_PHO_FRAC_PSIG!");
  
}

void In(const std::string str) { std::cout << redbg << "In  " << str << clearatt << std::endl; }
void Out(const std::string str){ std::cout << red << "Out " << str << clearatt << std::endl; }
void AtPoint(const std::string fun, const int line) { std::cout << blue << "AtPoint:"<< fun << ":" << line << clearatt << std::endl; }
bool Msg(const std::string func, const int line, const std::string msg="")
{
	std::cout << red << func << ":" << line << ":" << msg << clearatt << std::endl;
	return true;
}


//for debugging only  functions ------------------------

void MakeTruePhoFracStudy(const TH1* data, const TH1* mcpho, const TH1* qcd)
{
	
	TH1* hdata = (TH1*) data->Clone("datacopy");
	TH1* hmc = (TH1*) mcpho->Clone("mccopy");
	TH1* hmc_copy = (TH1*) mcpho->Clone("mccopy2");
	TH1* hqcd = (TH1*) qcd->Clone("qcdcopy");

	hmc_copy->Add(hqcd);
	hmc->Divide(hmc_copy);
	
	//WilsonInterval(hmc, hmc_copy);
	
	std::cout << "integrals data/mc/qcd/sum " 
			<< hdata->Integral() << "\t" 
			<< hmc->Integral() << "\t" 
			<< hqcd->Integral() << "\t" 
			<< hmc->Integral()+ hqcd->Integral()
			<< std::endl;

	std::cout << std::setw(3) << "bin" << std::setw(7) << "loedge"
		<< std::setw(10) << "data"
		<< std::setw(10) << "mc"
		<< std::setw(10) << "qcd"
		<< std::setw(15) << "ratio(mc/data)"
		<< std::setw(20) << "ratio(mc+qcd/data)"
		<< std::endl;
	for (int bin=1; bin<= hdata->GetNbinsX(); ++bin)
	{
	std::cout << std::setw(3) << bin<< std::setw(7) << hdata->GetXaxis()->GetBinLowEdge(bin)
		<< std::setw(10) << hdata->GetBinContent(bin)
		<< std::setw(10) << hmc->GetBinContent(bin)
		<< std::setw(10) << hqcd->GetBinContent(bin)
		<< std::setw(15) << hmc->GetBinContent(bin)/hdata->GetBinContent(bin)
		<< std::setw(20) << (hmc->GetBinContent(bin)+hqcd->GetBinContent(bin))/hdata->GetBinContent(bin)
		<< std::endl;
			//float deno = (1 - fFAKE_PHO_FRAC_DEF) * hmc->GetBinContent(bin) + fFAKE_PHO_FRAC_DEF * hqcd->GetBinContent(bin);
			//if ((hmc->GetBinContent(bin) + hqcd->GetBinContent(bin)) > 0 )
			//hmc->SetBinContent(bin,  hmc->GetBinContent(bin) / (1. * (hmc->GetBinContent(bin) + hqcd->GetBinContent(bin) ) ) );
			
		//hmc->SetBinContent(bin, (1-fFAKE_PHO_FRAC_DEF) * hmc->GetBinContent(bin));
	}
	//hmc->Divide(hdata);
	hmc->SetTitle("#gamma^{E_{T}>30GeV, |#eta|<1.1} + #geq 1 Jet^{E_{T}>15GeV, |#eta|<3.0} : Photon Purity, p(E_{T}^{#gamma},#eta) = #frac{ #epsilon N_{S}}{ #epsilon N_{s} + (1- #epsilon)N_{B}};E_{T}^{#gamma} (GeV); #frac{dp}{dE_{T}^{#gamma} d#eta^{#gamma}}");
	hmc->GetYaxis()->CenterTitle(1);
	hmc->GetXaxis()->CenterTitle(1);

	new TCanvas();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptStat(0);
	hmc->SetMarkerStyle(20);
	hmc->SetMarkerSize(1);
	hmc->SetMarkerColor(kBlue);
	
	hmc->Draw("E1");
	gPad->SetEditable(0);



}


//-----------------------------------------------------------------------------
// assumes all three hists have a the same bin sizes
void DebugJES(const TH1* hist, const TH1* jesup, const TH1* jesdown,
				 const int iLINE=0, const std::string title="",
				 const std::string filename = "")
{

	assert(hist != NULL && "hist_err null!.");
	assert(jesup != NULL && "jesup null!.");
	assert(jesdown != NULL && "jesdown null!.");

	TH1* chist = (TH1*) hist->Clone("hist_copy");
	TH1* cuhistrat = (TH1*) hist->Clone("hist_copy_uratio");
	TH1* cdhistrat = (TH1*) hist->Clone("hist_copy_dratio");
	TH1* uhist = (TH1*) jesup->Clone("jesup_copy");
	TH1* dhist = (TH1*) jesdown->Clone("jesdown_copy");

	TCanvas *c1 = new TCanvas();
	unsigned int w = 1200;
	unsigned int h = 600;
	c1->SetCanvasSize(w,h);
	c1->SetWindowSize(w + (w - c1->GetWw()) + 20, h + (h - c1->GetWh())+30);
	c1->Divide(2,1);
	c1->cd(1);
	//gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	chist->SetLineColor(kBlue);
	chist->SetMarkerColor(kBlue);
	chist->SetMarkerStyle (20);
	chist->SetMarkerSize(0.8);
	uhist->SetLineColor(kRed);
	uhist->SetMarkerColor(kRed);
	uhist->SetMarkerStyle (22);
	dhist->SetMarkerStyle (23);
	dhist->SetMarkerColor(kGreen);
	chist->SetMinimum(0.05);
	if (title.length()) chist->SetTitle(title.c_str());
	//find x max ann zoom-in
	double x_max = 0;
	for (int bin=chist->GetNbinsX(); bin>0; --bin)
	{
		if ( chist->GetBinContent(bin) || 
			  uhist->GetBinContent(bin) || 
			  dhist->GetBinContent(bin) )
		{
			x_max = chist->GetXaxis()->GetBinUpEdge(bin);
			break;
		}
	}
	assert (x_max > 0 && "DebugJES:: Max value for X axis retuned zero!");
	
	chist->GetXaxis()->SetRangeUser(0,x_max);
	chist->Draw("P");	
	uhist->Draw("sameP");	
	dhist->Draw("sameP");	

	TLegend *leg = new TLegend (0.7,0.72,0.9,0.9);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);

	leg->AddEntry(chist,"CENTRAL");
	leg->AddEntry(uhist,"JES UP");
	leg->AddEntry(dhist,"JES DOWN");
	//leg->AddEntry(chist,"DATA");
	//leg->AddEntry(uhist,"PHO MC");
	//leg->AddEntry(dhist,"PHO SIDEBAND");
	//leg->AddEntry(chist,"Background");
	//leg->AddEntry(chist,"FakeFrac+#sigma");
	//leg->AddEntry(uhist,"Scaled Pho MC");
	//leg->AddEntry(uhist,"Scaled Sideband");
	//leg->AddEntry(dhist,"PHO SIDEBAND");
	leg->Draw();

	//now make aratio hist of the differeces
	//zero out the two ratio hist (up/down)
	for (int bin=0;bin <= cuhistrat->GetNbinsX()+1; ++bin)
	{
			cuhistrat->SetBinContent(bin, 0);
			cuhistrat->SetBinError(bin, 0);
			cdhistrat->SetBinContent(bin, 0);
			cdhistrat->SetBinError(bin, 0);
	}
	
	
	for (int bin=1; bin <= chist->GetNbinsX(); ++bin)
	{
		if (chist->GetBinContent(bin))
		{
			double cval = chist->GetBinContent(bin);
			double uval = uhist->GetBinContent(bin);
			double dval = dhist->GetBinContent(bin);
			if ( std::isnan(cval) || std::isnan(uval) ||  std::isnan(dval))
			{
				std::cout << __FUNCTION__ << ":" << __LINE__ 
							<< ": one or bin values returned NAN!" << std::endl;
			}
			
			cuhistrat->SetBinContent(bin, uval/cval - 1);
			cdhistrat->SetBinContent(bin, dval/cval - 1);
		}

	}
	
	
	//new TCanvas();
	c1->cd(2);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();

	cuhistrat->SetMarkerStyle(22);
	cdhistrat->SetMarkerStyle(23);
	cuhistrat->SetMarkerColor(kRed);
	cdhistrat->SetMarkerColor(kGreen);
	
	cuhistrat->SetMinimum (-1.0);
	cuhistrat->SetMaximum (1.0);
	if (title.length()) cuhistrat->SetTitle(title.c_str());
	cuhistrat->GetXaxis()->SetRangeUser(0,x_max);
	cuhistrat->Draw("P");
	cdhistrat->Draw("P SAME");
	leg->Draw();

	c1->cd();
	
	if (filename.length())
	{
		c1->cd();
		c1->Print(filename.c_str());
	}

	//dump first few bins
	std::cout << __FUNCTION__<<  "::Caller Line#" << iLINE<< " - First 4 non-zero bins dump for :" << std::endl;
	std::cout << "CENTRAL - ";
	hist->Print();
	std::cout << "JES UP  - ";
	jesup->Print();
	std::cout << "JES DOWN- ";
	jesdown->Print();
	std::cout << setiosflags(ios::fixed) << setprecision(1) 
		<< setw(4) << "bin " << setw(12) << "lo/up edges" <<  setw(11) << "*central*" 
		<< setw(19) << "*jes up[diff]*" <<  setw(24) << "*jes down [diff]*"
		<< std::endl;
	int dumps = 0;
	for (int bin=1;bin <= chist->GetNbinsX()+1; ++bin)
	{
		if (chist->GetBinContent(bin))
		{
			double cval = chist->GetBinContent(bin);
			double uval = uhist->GetBinContent(bin);
			double dval = dhist->GetBinContent(bin);
			double udiff = cval - uval;
			double ddiff = cval - dval;
			std::cout << setw(4) << bin << setw(6) << hist->GetBinLowEdge(bin) 
				<< setw(1) << ", " << setw(4) << hist->GetXaxis()->GetBinUpEdge(bin)
				<< setw(13) << cval
				<< setw(12) << uval << "[ " << udiff << "]"
				<< setw(12) << dval << "[ " << ddiff << "]"
				<< std::endl;
			++dumps;
		}
		if (dumps>4) break;
	}
}

//end of debugging only functions  ------------------------
/*****************************************************************************/

Double_t fitFunction(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
		//val = (par[0]+par[1]*sqrt(x[0]))/x[0]+par[2];
		val = (par[0]+par[1]*sqrt(x[0]))+par[2]*x[0];
		//val = par[0]+par[1]*sqrt(x[0]);
	}
	return val;
}


/*******************************************************
 * This method is to check if there are statistically 
 * significant negative bins after subtractin off 
 * some components from a hist.
 * if found , issue a error!
 *******************************************************/
void CheckNegativeBin(TH1* hist)
{
	assert (hist != NULL && "CheckNegativeBin:: hist give is null");
	assert (hist->GetDimension() == 1 && "CheckNegativeBin:: hist dimension not 1!");
	
	for (int bin=1; bin <= hist->GetNbinsX(); ++bin)
	{
		if (hist->GetBinContent(bin)>0.0) continue;
		//now check if this bin value is statistically
		//within 0

		if ((hist->GetBinContent(bin) + hist->GetBinError(bin))>=0) continue;

		std::cout << red << __FUNCTION__ << ":" << __LINE__ 
			<< ": STATISTICALLY SIGNIFICANT NEGATIVE BIN VALUE FOUND! PLEAS CHECK! Setting val/err to zero" 
			 << "bin/val/err = " << bin << "\t" << hist->GetBinContent(bin)
			 << ", " << hist->GetBinError(bin) << clearatt << std::endl;
		hist->SetBinContent(bin,0);
		hist->SetBinError(bin,0);
	}
	

}

/*********************************************************
 * I am using this method to generate a reweight function
 * to match the sideband photon Et distribution to the
 * total background (pho sideband+ pho MC) photon Et
 * distribution. Call this method within the main
 * MakeHistLogAndRatioV2() method soon after rebinning and
 * fake photon fraction is defined.
 *
 * TO BE SAFE: add a return statement soon after 
 * this method is called.  - 03-09-2010
 *********************************************************/
void MakePhoEtWeightsForSideband( TH1* hqcd, TH1* hmcpho, 
						const float fake_frac, const float fake_frac_sigma, 
						const std::string brname, const int njets)
{

	assert (hqcd != NULL && "MakePhoEtWeightsForSideband:: qcd hist is null!");
	assert (hqcd->GetDimension() == 1 && "MakePhoEtWeightsForSideband:: qcd hist is not 1D!");
	assert (hmcpho != NULL && "MakePhoEtWeightsForSideband:: mcpho hist is null!");
	assert (hmcpho->GetDimension() == 1 && "MakePhoEtWeightsForSideband:: mcpho hist is not 1D!");
	assert ( (fake_frac>0 && fake_frac < 1.) && "MakePhoEtWeightsForSideband:: not a valid fake photon fraction!");

	cout << __FUNCTION__ << ": Settings passed in:" << endl;
	cout << "Fake photon fraction = " << fake_frac << "+/-" << fake_frac_sigma << endl;

	/*hqcd->Print("all");
	hmcpho->Print("all");
	new TCanvas();
	gPad->SetLogy();
	hqcd->DrawCopy();
	hmcpho->DrawCopy("same");
	return;
	*/
	
/*	new TCanvas();
	THStack *hs3 = new THStack ("hs3", NULL);
	hs3->Add(hqcd);
	hs3->Add(hmcpho);
	gPad->SetLogy();
	hs3->Draw("HIST");
	return;
*/
	
	TH1* hist_qcd = (TH1*) hqcd->Clone("hist_nominal_wghts");
	hist_qcd->SetLineColor(kYellow);
	hist_qcd->SetFillColor(kYellow);
	hist_qcd->SetMarkerColor(kYellow);
	TH1* hist_qcd2 = (TH1*) hqcd->Clone("qcd_copy2");
	TH1* hist_qcd_delta = (TH1*) hqcd->Clone("hist_plusSigma_wghts");
	TH1* hist_qcd_delta2 = (TH1*) hqcd->Clone("qcd_copy4");
	TH1* residual = dynamic_cast<TH1*> (hqcd->Clone("residual"));
	TH1* hist_mcpho = (TH1*) hmcpho->Clone("mcpho_copy");
	hist_mcpho->SetLineColor(6);
	hist_mcpho->SetFillColor(6);
	hist_mcpho->SetMarkerColor(6);
	TH1* hist_mcpho_delta = (TH1*) hmcpho->Clone("mcpho_copy2");
	hist_qcd->Print();
	hist_mcpho->Print();



/*	new TCanvas();
	hist_qcd->DrawCopy("L");
	hist_mcpho->DrawCopy("same L");
	hist_mcpho->Add(hist_qcd);
	hist_mcpho->SetLineColor(kRed);
	hist_mcpho->DrawCopy("same L");
	//return;
*/


	
	//sumw2 is already called for original hists

	hist_qcd->Scale(fake_frac / (double) hist_qcd->Integral());
	hist_mcpho->Scale((1 - fake_frac) / (double) hist_mcpho->Integral());
	hist_qcd->Add(hist_mcpho);

	//normalize back to orignal
	hist_qcd->Scale(hist_qcd2->Integral()/(double)hist_qcd->Integral());
	hist_qcd->Divide(hist_qcd2);

	std::stringstream basetitle;
	basetitle << "#gamma^{E_{T}>30GeV,|#eta|<1.1}+ #geq " << njets << " Jet^{E_{T}>15 GeV,|#eta|<3.0}+#slash{E}_{T}>20 GeV: E_{T}^{#gamma} : ";
	//basetitle << "#gamma^{E_{T}>30GeV,|#eta|<1.1}+ == " << njets << " Jet^{E_{T}>15GeV, |#eta|<3.2}: E_{T}^{#gamma} : ";
	std::stringstream htitle;
	htitle << basetitle.str() << "True-#gamma fraction = #epsilon =" << 1-fake_frac;
	
	hist_qcd->SetTitle(htitle.str().c_str());
	hist_qcd->GetYaxis()->SetTitle("(#epsilon  #gamma-MC + (1 - #epsilon)Sideband )/Sideband");

	std::stringstream c1_name;
	//c1_name << "Fit_" << brname<<"_error";
	c1_name << "Fit_" << brname;
	//TCanvas *c1 = new TCanvas(c1_name.str().c_str(),c1_name.str().c_str(),1000,800);
	new TCanvas();
	//new TCanvas();
	//c1->Divide(1,2);
	//c1->cd(1);

	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptFit(1);
	//hist_qcd->SetStats(0);

	hist_qcd->SetLineColor(kBlue);
	//hist_qcd->SetMarkerColor(kBlue);
	hist_qcd->SetMarkerStyle (20);
	hist_qcd->SetMarkerSize(1);


	
	//find fit range
	int xminbin=0, xmaxbin=0;
	FindNonZeroXrange(hist_qcd,xminbin,xmaxbin);
	xmaxbin = xmaxbin;
	std::cout << "minx, maxx = " << xminbin << "\t" << xmaxbin << std::endl;
	TF1 *f4 = new TF1("fit_func4",fitFunction,hist_qcd->GetBinCenter(xminbin),hist_qcd->GetXaxis()->GetBinCenter(xmaxbin),3);
	f4->SetParameter(0,0.5);
	f4->SetParameter(1,0.5);
	f4->SetParameter(2,30);
	f4->SetLineColor(8);

	hist_qcd->Fit(f4,"R+");
	//hist_qcd->Fit("pol1","R","",hist_qcd->GetBinCenter(xminbin),hist_qcd->GetXaxis()->GetBinCenter(xmaxbin));

	TList *funcList = hist_qcd->GetListOfFunctions();
	TIterator *it = funcList->MakeIterator();
	TF1* fitfunc = (TF1*) it->Next();
	if (fitfunc == NULL)
	{
		std::cout << red << __FUNCTION__ << ":" << __LINE__  
			<< "::Could not retrieve fit function! returning !" << clearatt << std::endl;
	}
	fitfunc->SetNameTitle(brname.c_str(),htitle.str().c_str());
	//fitfunc->SetNameTitle("tf1_leadjeteta_nominalweights", htitle.str().c_str());

	//residual of the fit
	residual->GetYaxis()->SetTitle("#frac{Fit -- Val}{Error}");
	for (int bin=1; bin < residual->GetNbinsX(); ++bin)
	{
		float val = hist_qcd->GetBinContent(bin);

		float err = hist_qcd->GetBinError(bin);
		float fit = fitfunc->Eval(residual->GetBinCenter(bin));
		if (err>0) residual->SetBinContent(bin,(fit-val)/err);
		residual->SetBinError(bin,0);
	}

/*	c1->cd(2);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	residual->Draw();
	c1->cd();
*/

	
/******************************************/
//uncertainty for this method
//derive another set of wgt with fake_frac+/-sigma
/******************************************/
	// this is to derive the syst error for the above see elog#1520
	// delta w = (MC/SB - 1)  * pho_frac_sigma
	//hist_mcpho_delta->Scale(fake_frac_sigma/(double)hist_mcpho_delta->Integral());
	//hist_qcd_delta->Scale(fake_frac_sigma/(double)hist_qcd_delta->Integral());
	//hist_mcpho_delta->Add(hist_qcd_delta, -1);
	//hist_mcpho_delta->Scale(hist_qcd_delta2->Integral()/(double)hist_mcpho_delta->Integral());
	//hist_mcpho_delta->Divide(hist_qcd_delta2);
	//new TCanvas();
	//hist_mcpho_delta->DrawClone();
	hist_qcd_delta->Scale( (fake_frac +fake_frac_sigma) / (double) hist_qcd_delta->Integral());
	hist_mcpho_delta->Scale((1 - (fake_frac + fake_frac_sigma)) / (double) hist_mcpho_delta->Integral());
	hist_qcd_delta->Add(hist_mcpho_delta);

	//normalize back to orignal
	hist_qcd_delta->Scale(hist_qcd_delta2->Integral()/(double)hist_qcd_delta->Integral());
	hist_qcd_delta->Divide(hist_qcd_delta2);

	
	std::stringstream htitle_delta;
	htitle_delta << basetitle.str() << " : True-#gamma fraction = #epsilon - #sigma = " << 1 - (fake_frac + fake_frac_sigma);
	
	hist_qcd_delta->SetTitle(htitle_delta.str().c_str());
	hist_qcd_delta->GetYaxis()->SetTitle("( (#epsilon+#sigma) #gamma-MC + (1 - (#epsilon+#sigma)) Sideband )/Sideband");
	hist_qcd_delta->GetYaxis()->SetTitleSize(0.04);

/*	std::stringstream c2_name;
	c2_name << "Fit_" << brname << "_plusSigma";
	TCanvas *c2 = new TCanvas(c2_name.str().c_str(),c2_name.str().c_str(),1000,800);
	c2->Divide(1,2);
	c2->cd(1);

	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptFit(1);
	*/
	//hist_qcd->SetStats(0);

	hist_qcd_delta->SetLineColor(kBlue);
	hist_qcd_delta->SetMarkerStyle (20);
	hist_qcd_delta->SetMarkerSize(1);

	//find fit range
	int xminbin_delta=0, xmaxbin_delta=0;
	FindNonZeroXrange(hist_qcd_delta,xminbin_delta,xmaxbin_delta);
	std::cout << "minx_delta, maxx_delta = " << xminbin_delta << "\t" << xmaxbin_delta << std::endl;
	TF1 *f4p = new TF1("fit_func4_delta",fitFunction,hist_qcd_delta->GetBinCenter(xminbin_delta),hist_qcd_delta->GetXaxis()->GetBinCenter(xmaxbin_delta),3);
	f4p->SetParameter(0,0.5);
	f4p->SetParameter(1,0.5);
	f4p->SetParameter(2,30);
	f4p->SetLineColor(8);

	//c1->cd(2);
	new TCanvas();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetOptFit(1);

	//hist_qcd_delta->Fit(f4p,"R+");
	hist_qcd_delta->Fit("pol2","R","",hist_qcd->GetBinCenter(xminbin),hist_qcd->GetXaxis()->GetBinCenter(xmaxbin));


return;
	TList *funcList_delta = hist_qcd_delta->GetListOfFunctions();
	TIterator *it_delta = funcList_delta->MakeIterator();
	TF1* fitfunc_delta = (TF1*) it_delta->Next();
	std::stringstream tf1deltaname;
	tf1deltaname << brname << "_plusSigma";
	fitfunc_delta->SetNameTitle(tf1deltaname.str().c_str(),htitle_delta.str().c_str());

	//residual of the fit
/*	residual->GetYaxis()->SetTitle("#frac{Fit -- Val}{Error}");
	for (int bin=1; bin < residual->GetNbinsX(); ++bin)
	{
		float val = hist_qcd->GetBinContent(bin);
		float err = hist_qcd->GetBinError(bin);
		float fit = fitfunc->Eval(residual->GetBinCenter(bin));
		if (err>0) residual->SetBinContent(bin,(fit-val)/err);
		residual->SetBinError(bin,0);
	}

	c1->cd(2);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();
	residual->Draw();
*/
	//c1->cd();

/******************************************/


	
	//now save all above to a hist
	
	//check if there is a saved version already.
	//if so delete that before writing the new one.
	//otherwise there will be many copied of the same obj
	//with different cycles
	//TFile f("g40jetsmet40_SidebandReweightFunctions.root","UPDATE");
	//TFile f("g30jets_SidebandPhoEtReweightFunctions.root","UPDATE");
	//TFile f("SidebandWgtFunctionFromBckgOnly_PhoFracPlusSigma.root","UPDATE");
	//TFile f("SidebandWgtFunctionFromBckgOnly_PhoFracMinusSigma.root","UPDATE");
	TFile f("temp.root","UPDATE");
	if ( f.Get(brname.c_str()) != NULL)
	{
		std::stringstream old_tf;
		old_tf << brname << ";1";
		f.Delete(old_tf.str().c_str());
	}
	fitfunc->Write();

	if ( f.Get(fitfunc_delta->GetName()) != NULL)
	{
		std::stringstream old_tf_delta;
		old_tf_delta << fitfunc_delta->GetName() << ";1";
		f.Delete(old_tf_delta.str().c_str());
	}
	fitfunc_delta->Write();

	
	if ( f.Get(c1_name.str().c_str()) != NULL)
	{
		std::stringstream old_c1;
		old_c1 << c1_name.str() << ";1";
		f.Delete(old_c1.str().c_str());
	}
	//c1->cd();
	//c1->Write();
/*
	if ( f.Get(c2->GetName()) != NULL)
	{
		std::stringstream old_c2;
		old_c2 << c2->GetName() << ";1";
		f.Delete(old_c2.str().c_str());
	}
	c2->cd();
	c2->Write();
*/

	
	if ( f.Get(hist_qcd->GetName()) != NULL)
	{
		std::stringstream old_th1f;
		old_th1f << hist_qcd->GetName() << ";1";
		f.Delete(old_th1f.str().c_str());
	}
	hist_qcd->Write();

	if ( f.Get(hist_qcd_delta->GetName()) != NULL)
	{
		std::stringstream old_th1f_delta;
		old_th1f_delta << hist_qcd_delta->GetName() << ";1";
		f.Delete(old_th1f_delta.str().c_str());
	}
	hist_qcd_delta->Write();

	
	f.ls();
	f.Close();

}


Double_t AlphaSsystFunc(Double_t *x, Double_t *par)
{
	double val = 0;
	if (x[0]>0)
	{
		 val = (par[0]+par[1]*sqrt(x[0]))/x[0]+par[2];
	}
	return val;
}


void GetAlphaS_Systematic(const std::string brname, std::vector<float> &vErr, 
								const TH1* hist_mc)
{
	assert (hist_mc != NULL && "GetAlphaS_Systematic:: hist is null!");
	assert (hist_mc->GetDimension() == 1 && "GetAlphaS_Systematic:: hist is not 1D!");

	TFile alphaSyst("AlphaS_Syst_FitFunctions.root");
	assert (! alphaSyst.IsZombie() && "AlphaS_Syst file not fouund!");
	//alphaSyst.ls();


	TF1* tf1_temp = dynamic_cast<TF1*> (alphaSyst.Get(brname.c_str()));
	assert (tf1_temp != NULL && "GetAlphaS_Systematic:: Fit function not found to read the parameters!");

	const int NparNeed = 3;
	
	assert (tf1_temp->GetNdim() == 1 && "GetAlphaS_Systematic:: Read fit function is not 1D!");
	assert (tf1_temp->GetNpar() == NparNeed && "GetAlphaS_Systematic:: Read fit function NPar does not match with expected!");

	TF1 *tf1A = new TF1("AlaphaS", AlphaSsystFunc,
			hist_mc->GetBinCenter(1), hist_mc->GetBinCenter(hist_mc->GetNbinsX()), NparNeed);

	double par[NparNeed];

	tf1_temp->GetParameters(par);
	for (int i =0; i <NparNeed; ++i)
	{
	//	std::cout << "par[" << i << "]=" << par[i] << std::endl;
	}
	tf1A->SetParameters(par);
	//tf1A->Print();

	std::stringstream hist_name;
	hist_name << "TH1F_"<< brname;
	TH1* hist_temp = dynamic_cast<TH1*> (alphaSyst.Get(hist_name.str().c_str()));
	assert (hist_temp != NULL && "GetAlphaS_Systematic:: Required hist not found!");
	//new TCanvas();
	//hist_temp->Draw();
	//gPad->SetEditable(0);
	

	float first_bin_val = 0;
	for (int bin = 1; bin<= hist_temp->GetNbinsX(); ++bin)
	{
		if (hist_temp->GetBinContent(bin)>0) 
		{
			first_bin_val = hist_temp->GetBinContent(bin);
			break;
		}
	}
	//delete hist_temp;
	assert (first_bin_val>0 && "GetAlphaS_Systematic:: No first_bin_val found!");
	//std::cout << "first_bin_val = " << first_bin_val << std::endl;

	if (tf1A != NULL)
	{
		vErr.clear();// removed stuffed zeros
		for (int bin = 1; bin<= hist_mc->GetNbinsX(); ++bin)
		{
			float fAe = tf1A->Eval(hist_mc->GetBinCenter(bin)); 
			//keep minimum error at first bin val. fits are not perfect elog#1514
			if (fAe < first_bin_val) fAe = first_bin_val;
			vErr.push_back(fAe * hist_mc->GetBinContent(bin));
			//std::cout << green << bin << "[" << hist_mc->GetBinCenter(bin) << "] = " << fAe << std::endl;
		}
		delete tf1A;
	} else 
	{
		std::cout << __LINE__ << ":: WARNING !  NO Alpha-S Systematic function found for "  << brname  << std::endl;
	}

}


void DrawErrorBreakDown(const int QCDerrMtd, 
		const int iSideband,
		const TH1* hist_err,
		const std::vector<float>& JESerr, const std::vector<float>&  qcdmcMixErr, 
		const std::vector<float>& vTotalEWKLumErr,
		const std::vector<float>& cosmicErr, const std::vector<float>& haloErr, 
		const std::vector<float>&  statErr, const std::vector<float>& IdErr, 
		const std::vector<float>& vPDFerror, const std::vector<float>&  vQ2error, 
		const std::vector<float>& vISRFSRerror, const std::vector<float>& vAlphaSerror, 
		const std::vector<float>& vSidebandRwtErr, const std::vector<float> vMtdBerr,
		const std::vector<float>& ERROR, const std::string title,
		const int jets, const std::string which
		)
{
	if (iDEBUG > iFUNCTION_NAMES) In(__FUNCTION__);

	assert ( (QCDerrMtd == 1 || QCDerrMtd == 2) 
			&& "DrawErrorBreakDown:: QCDerrMtd is not valid!");
	assert ( hist_err != NULL 
			&& "DrawErrorBreakDown:: hist_err is null!");
	assert (hist_err->GetDimension() == 1 
			&& "DrawErrorBreakDown:: hist_err is is not 1D!");

	TH1::AddDirectory(kFALSE);
	TH1* hist = dynamic_cast<TH1*> (hist_err->Clone("Total Systematic"));
	
	for (int bin=0; bin <= hist->GetNbinsX()+1;++bin)
	{
		hist->SetBinContent(bin,0);
		hist->SetBinError(bin,0);
	}

	hist->SetTitle(title.c_str());

	TH1* h_jes = dynamic_cast<TH1*> (hist->Clone("JES"));
	TH1* h_phoFrac = dynamic_cast<TH1*> (hist->Clone("True Photon Fration"));
	TH1* h_EWKLum = dynamic_cast<TH1*> (hist->Clone("Total EWK Lum"));
	TH1* h_cosmic = dynamic_cast<TH1*> (hist->Clone("Cosmic"));
	TH1* h_halo = dynamic_cast<TH1*> (hist->Clone("Beam Halo"));
	TH1* h_stat = dynamic_cast<TH1*> (hist->Clone("Statistical"));
	TH1* h_qcdId = dynamic_cast<TH1*> (hist->Clone("QCD ID"));
	TH1* h_pdf = dynamic_cast<TH1*> (hist->Clone("PDF"));
	TH1* h_q2 = dynamic_cast<TH1*> (hist->Clone("Q2"));
	TH1* h_isrfsr = dynamic_cast<TH1*> (hist->Clone("ISR/FSR"));
	TH1* h_alphas = dynamic_cast<TH1*> (hist->Clone("Alpha_s"));
	TH1* h_sidebandrewgt = dynamic_cast<TH1*> (hist->Clone("Sideband Reweighing"));
	TH1* h_mtdBerr = dynamic_cast<TH1*> (hist->Clone("Mtd B errors"));

	for (int bin=1; bin <= hist->GetNbinsX();++bin)
	{
		if (hist_err->GetBinContent(bin))
		{
			hist->SetBinContent(bin,ERROR.at(bin-1));
			h_jes->SetBinContent(bin,JESerr.at(bin-1));
			h_cosmic->SetBinContent(bin,cosmicErr.at(bin-1));
			h_halo->SetBinContent(bin,haloErr.at(bin-1));
			h_stat->SetBinContent(bin,statErr.at(bin-1));
			h_qcdId->SetBinContent(bin,IdErr.at(bin-1));
			h_EWKLum->SetBinContent(bin,vTotalEWKLumErr.at(bin-1));

			if (iSideband == iWEIGHTED_SIDEBAND)
			{
				h_sidebandrewgt->SetBinContent(bin,vSidebandRwtErr.at(bin-1));
			} else //for nominal sideband usage to make final plots
			{
				h_mtdBerr->SetBinContent(bin,vMtdBerr.at(bin-1));
			}

			//if (QCDerrMtd == 2)
			{
				h_phoFrac->SetBinContent(bin,qcdmcMixErr.at(bin-1));
				h_pdf->SetBinContent(bin,vPDFerror.at(bin-1));
				h_q2->SetBinContent(bin,vQ2error.at(bin-1));
				h_isrfsr->SetBinContent(bin,vISRFSRerror.at(bin-1));
				h_alphas->SetBinContent(bin,vAlphaSerror.at(bin-1));
			}
		}

	}



	hist->SetMarkerStyle(29); hist->SetMarkerColor(kRed);
	hist->SetMarkerSize(1.5);
	h_jes->SetMarkerStyle(20);
	h_cosmic->SetMarkerStyle(21);
	h_cosmic->SetMarkerColor(30);
	h_halo->SetMarkerStyle(22);
	h_halo->SetMarkerColor(46);

	h_stat->SetMarkerStyle(23);
	h_stat->SetMarkerColor(33);
	h_qcdId->SetMarkerStyle(20);
	h_qcdId->SetMarkerColor(38);

	h_phoFrac->SetMarkerStyle(2); h_phoFrac->SetMarkerSize(1.5);
	h_phoFrac->SetMarkerColor(7);

	h_pdf->SetMarkerStyle(20);
	h_pdf->SetMarkerColor(kBlue);
	h_q2->SetMarkerStyle(21);
	h_q2->SetMarkerColor(kGreen);
	h_isrfsr->SetMarkerStyle(3);
	h_isrfsr->SetMarkerColor(8);
	h_alphas->SetMarkerStyle(23);
	h_alphas->SetMarkerColor(6);

	h_EWKLum->SetMarkerStyle(23);
	h_EWKLum->SetMarkerColor(kRed);

	if (iSideband == iWEIGHTED_SIDEBAND)
	{
		h_sidebandrewgt->SetMarkerStyle(23);
		h_sidebandrewgt->SetMarkerColor(kBlue);
	} else 
	{
		h_mtdBerr->SetMarkerStyle(21);
		h_mtdBerr->SetMarkerColor(28);
	}
	

	TLegend *leg = new TLegend (0.74,0.30,0.99,1.);
	leg->SetTextFont(42);
	leg->SetTextSize(0.03);
	leg->AddEntry(hist,hist->GetName());
	leg->AddEntry(h_jes,h_jes->GetName());
	leg->AddEntry(h_cosmic,h_cosmic->GetName());
	leg->AddEntry(h_halo,h_halo->GetName());
	leg->AddEntry(h_stat,h_stat->GetName());
	leg->AddEntry(h_qcdId,h_qcdId->GetName());
	leg->AddEntry(h_EWKLum,h_EWKLum->GetName());

	
	TCanvas *c2 = new TCanvas("err_log","ERRORHIST",600,800);
	c2->Divide(1,2);
	c2->cd(1);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetLogy();

	//Absolute errors, LOG plot
	hist->GetYaxis()->SetTitle("Absolute Total Error");
	hist->DrawCopy("P");
	h_jes->DrawCopy("P SAME");
	h_cosmic->DrawCopy("P SAME");
	h_halo->DrawCopy("P SAME");
	h_stat->DrawCopy("P SAME");
	h_qcdId->DrawCopy("P SAME");
	h_EWKLum->DrawCopy("P SAME");
	if (iSideband == iWEIGHTED_SIDEBAND)
	{
		h_sidebandrewgt->DrawCopy("P SAME");
		leg->AddEntry(h_sidebandrewgt,h_sidebandrewgt->GetName());
	} else
	{
		h_mtdBerr->DrawCopy("P SAME");
		leg->AddEntry(h_mtdBerr,h_mtdBerr->GetName());
	}
	h_phoFrac->DrawCopy("P SAME");
	h_pdf->DrawCopy("P SAME");
	h_q2->DrawCopy("P SAME");
	h_isrfsr->DrawCopy("P SAME");
	h_alphas->DrawCopy("P SAME");


	leg->AddEntry(h_phoFrac,h_phoFrac->GetName());
	leg->AddEntry(h_pdf,h_pdf->GetName());
	leg->AddEntry(h_q2,h_q2->GetName());
	leg->AddEntry(h_isrfsr,h_isrfsr->GetName());
	leg->AddEntry(h_alphas,h_alphas->GetName());

	leg->Draw();


	//TCanvas *c3 = new TCanvas("err_rat","ERROR BREAKDOWN RATIO");
	
	for (int bin=1; bin <= hist->GetNbinsX();++bin)
	{
		float val = 0;
		if (hist_err->GetBinContent(bin)>0) val = 1/hist_err->GetBinContent(bin);

		hist->SetBinContent(bin,ERROR.at(bin-1) * val);
		h_jes->SetBinContent(bin,JESerr.at(bin-1) * val);
		h_cosmic->SetBinContent(bin,cosmicErr.at(bin-1) * val);
		h_halo->SetBinContent(bin,haloErr.at(bin-1) * val);
		h_stat->SetBinContent(bin,statErr.at(bin-1) * val);
		h_qcdId->SetBinContent(bin,IdErr.at(bin-1) * val);
		h_EWKLum->SetBinContent(bin,vTotalEWKLumErr.at(bin-1) * val);

		h_sidebandrewgt->SetBinContent(bin, vSidebandRwtErr.at(bin-1) * val); //this vector is stuffed with zeros when this error is not calcualted. and will not be drawn
		h_mtdBerr->SetBinContent(bin,vMtdBerr.at(bin-1) * val);

		h_phoFrac->SetBinContent(bin,qcdmcMixErr.at(bin-1) * val);
		h_pdf->SetBinContent(bin,vPDFerror.at(bin-1) * val);
		h_q2->SetBinContent(bin,vQ2error.at(bin-1) * val);
		h_isrfsr->SetBinContent(bin,vISRFSRerror.at(bin-1) * val);
		h_alphas->SetBinContent(bin,vAlphaSerror.at(bin-1) * val);
	}

	
	c2->cd(2);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetGridx();
	gPad->SetGridy();

	hist->SetMaximum(1);
	hist->GetYaxis()->SetTitle("Error/Prediction");
	hist->Draw("P");
	h_jes->Draw("P SAME");
	h_cosmic->Draw("P SAME");
	h_halo->Draw("P SAME");
	h_stat->Draw("P SAME");
	h_qcdId->Draw("P SAME");
	h_EWKLum->Draw("P SAME");

	if (iSideband == iWEIGHTED_SIDEBAND)
	{
		h_sidebandrewgt->Draw("P SAME");
	} else
	{
		h_mtdBerr->Draw("P SAME");
	}

	//if (QCDerrMtd == 2)
	{
		h_phoFrac->Draw("P SAME");
		h_pdf->Draw("P SAME");
		h_q2->Draw("P SAME");
		h_isrfsr->Draw("P SAME");
		h_alphas->Draw("P SAME");
	}

	leg->Draw();

	std::ostringstream str;
	//str << "ErrBreakdown_plot" << jets << "_" << which << ".eps";
	if (iSideband == iWEIGHTED_SIDEBAND)
	{
		str << sPRINTFOLDER << sSAMPLE_NAME << "_Errs_MtdB_plot" << jets << "_" << which << "." << sPRINTFORMAT.c_str();
	} else {
		str << sPRINTFOLDER << sSAMPLE_NAME << "_Errs_MtdA_plot" << jets << "_" << which << "." << sPRINTFORMAT.c_str();
	}

	c2->Print (str.str().c_str());

	if (iDEBUG > iFUNCTION_NAMES) Out(__FUNCTION__);

}


void ZeroOutNegativeBins(TH1* hist)
{
	//this is hack to remove bins with val<0
	//set bin value to zero if val<0
	assert (hist != NULL && "hist is null!");
	assert (hist->GetDimension() == 1 && "not a 1D hist");

	for (int bin=1; bin <=hist->GetNbinsX(); ++bin)
	{
		if (hist->GetBinContent(bin)<0)
		{
			hist->SetBinContent(bin,0);
			//hist->SetBinError(bin,0);
		}
	}
}

//-----------------------------------------------------------------------------
void GetCosmicErr(std::vector<float>& ErrVec, const TH1 *hist , const bool debug)
{
	if (debug) std::cout << "IN COSMIC ERROR"<< std::endl;
	assert (hist!=NULL && "GetComicErr::hist is null");
	ErrVec.clear();
	for (int bin=1; bin <= hist->GetNbinsX(); bin++) {
		float val = hist->GetBinContent(bin);
		float err = 0;
		//if (val) err = 1 / sqrt(val);		//take stat error as the syst
		//if (val) err = 1 ;		//????????????????? for now 04-15-2010
		if (val) err = sqrt(val);		//????????????????? for now 04-15-2010
		ErrVec.push_back(err);
		//std::cout << "cosmic bin[" << bin << "]=" << val  << " err = " << err << "[" << hist->GetBinError(bin) << "]" << std::endl;
		if (debug) std::cout << "cosmic bin[" << bin << "]=" << val << std::endl;
	}
}

//-----------------------------------------------------------------------------
void GetHaloErr(std::vector<float>& ErrVec, const TH1 *hist, bool debug=false)
{
	assert(hist!=NULL && "hist is null!");
	ErrVec.clear();
	for (int i=1; i <= hist->GetNbinsX(); i++) {
		float bin = hist->GetBinContent(i);
		float err = bin * 0.5;		//take 50% to be the syst
		ErrVec.push_back(err);
		if (debug) std::cout << "halo err [" << i << "]="<< err <<std::endl;
	}
	if (debug)
	{
		std::cout << "halo err arr size = " << ErrVec.size() << std::endl;
		std::cout << "\n\n\n i am out" << std::endl;
	}
}




//-----------------------------------------------------------------------------
// effects: make a new binning with variable bin sizes, and the given
//   cutoff points between bin sizes
// returns: the binning (by bin edges)
// guarantee: strong
//// throws: std::bad_alloc
//-----------------------------------------------------------------------------
std::vector<float> 
GetVarBinVector (float down, float up5, float up10, float up20, float up50, 
					float width1, float width2, float width3, float width4 ,bool debug)
{
	std::vector<float> result;
	float point = down;
	const unsigned nregion = 4;
	const float up [] = {up5, up10, up20, up50};
	const float step [] = {width1, width2, width3, width4};		//1j pet

	result.push_back (point);
	for (unsigned region = 0; region != nregion; ++ region)
	{
		while (point + step[region] < up[region] + 1)
		{
			point += step[region];
			result.push_back (point);
	 	};
	};
	
	if (point + 1 < up[nregion-1])
		result.push_back (up[nregion-1]);
	if (debug) std::cout << "Done creating var bin vector. " << std::endl;
	
	return result;
};
//-----------------------------------------------------------------------------
// effects: turn this histogram into a histogram with variable bin
//   size
// returns: the new histogram
// guarantee: strong
// throws: std::bad_alloc
// requires: input != NULL
// requires: input->GetDimension() == 1
// requires: histogram bin edges must match
// postcondition: result != NULL
//-----------------------------------------------------------------------------
std::auto_ptr<TH1>
make_var_bin (const std::vector<float>& bins,
	      TH1 *input, bool debug)
{
	assert (input != NULL && "requirement failed"); //spec
	assert (input->GetDimension() == 1 && "requirement failed"); //spec
	const unsigned nbin = unsigned (input->GetNbinsX ());

	std::auto_ptr<TH1> result (new TH1F (input->GetName(), input->GetTitle(),
						 bins.size() - 1, &bins[0]));

	result->SetBinContent (0, input->GetBinContent (0));
	result->SetBinError (0, input->GetBinError (0));
	result->SetBinContent (bins.size(), input->GetBinContent (nbin));
	result->SetBinError (bins.size(), input->GetBinError (nbin));
	for (unsigned bin = 1; bin != nbin; ++ bin)
	{
		const float low = input->GetXaxis()->GetBinLowEdge (bin);
		const float high = input->GetXaxis()->GetBinUpEdge (bin);
		const unsigned bin1 = result->FindBin (0.99 * low + 0.01 * high);
		const unsigned bin2 = result->FindBin (0.01 * low + 0.99 * high);
		assert (bin1 == bin2 && "requirement failed: histogram bin edges don't match"); //spec

		const float va = input->GetBinContent (bin);
		const float ea = input->GetBinError (bin);
		const float vb = result->GetBinContent (bin2);
		const float eb = result->GetBinError (bin2);
		const float vc = va + vb;
		const float ec = sqrt (ea * ea + eb * eb);
		result->SetBinContent (bin2, vc);
		result->SetBinError (bin2, ec);
	};


	if (debug)
	{
		for (unsigned i=1;  i<=(unsigned)input->GetNbinsX (); i++)
		{
			std::cout << "input bin, lowedge = " << i << ", " << input->GetXaxis()->GetBinLowEdge(i) << "\t" << input->GetBinContent(i) << std::endl;
		}
		std::cout << "\n\n"<< std::endl;
		for (unsigned i=1; i<=(unsigned)result->GetNbinsX (); i++)
		{
			std::cout << "otput bin, lowedge = " << i << ", " << result->GetXaxis()->GetBinLowEdge(i) << "\t" << result->GetBinContent(i) << std::endl;
		}
	}

	assert (result.get() != NULL && "postcondition failed"); //spec
	return result;
};

//-----------------------------------------------------------------------------
// effects: normalize each bin by the bin width
// side effect: normalizes histogram bins in hist
// guarantee: no-throw
// requires: hist != NULL
// requires: hist->GetDimension() == 1
//-----------------------------------------------------------------------------
void norm_bins (TH1 *hist, const float fOldBinWidth)
{
	assert (hist != NULL && "requirement failed"); //spec
	assert (hist->GetDimension() == 1 && "requirement failed"); //spec
	assert ( fOldBinWidth>0 && "norm_bins:: fOldBinWidth must be >0!");

	for (unsigned bin = 1; bin != unsigned (hist->GetNbinsX()); ++ bin)
	{
		const float width = hist->GetXaxis()->GetBinWidth (bin);
		//hist->SetBinContent (bin, hist->GetBinContent (bin) / width);
		//hist->SetBinError (bin, hist->GetBinError (bin) / width);


		//hist->SetBinContent (bin, hist->GetBinContent (bin) * fOldBinWidth / width);
		//hist->SetBinError (bin, hist->GetBinError (bin) *fOldBinWidth/ width);
		hist->SetBinContent (bin, hist->GetBinContent (bin) / width);
		hist->SetBinError (bin, hist->GetBinError (bin) / width);
	};
};

//-----------------------------------------------------------------------------
// make a variable binned hist with given specification
//-----------------------------------------------------------------------------
//TH1* MakeVariableBins (TH1 *hist, const float xmin, const float xpoint1,
//							const float xpoint2, const float xpoint3, const float xpoint4,
//				 			const float width1, const float width2, const float width3, 
//							const float width4 , bool debug=false)
TH1* MakeVariableBins (TH1 *hist, bool debug=false)

{
	if (debug)
	{
		std::cout << hist->GetName() << " hist bin size=" << hist->GetBinWidth(1) 
					<< std::endl;
	}

  	std::auto_ptr<TH1> result = make_var_bin ( 
											GetVarBinVector(fXMIN, fXPOINT1, fXPOINT2, 
														fXPOINT3, fXPOINT4, fWIDTH1, 
														fWIDTH2, fWIDTH3, fWIDTH4, debug) 
											, hist, debug);
  	norm_bins (result.get(), hist->GetBinWidth(1));
  	return result.release ();
};


//-----------------------------------------------------------------------------
//this make the ratio plot
//Assumes that hist_err has total syst+stat errors
////-----------------------------------------------------------------------------
void MakeRatioPlot(const TH1* phojet,
						const TH1* hist_err, const std::string xtitle, 
						TPad *pad=0, bool debug=false)
{
	assert(phojet!=NULL && "phojet hist null");
	assert(hist_err!=NULL && "hist_err hist null");

	TH1* pho = (TH1*) phojet->Clone("phojet_copy");
	TH1* err = (TH1*) hist_err->Clone("hist_err_copy");
	
	for (unsigned bin = 1; bin <= (unsigned)err->GetNbinsX(); ++ bin)
	{
		const float val = err->GetBinContent (bin);
		const float scale = val ? 1. / val : 0;
		pho->SetBinContent (bin, (pho->GetBinContent (bin) - val) * scale);
		pho->SetBinError (bin, pho->GetBinError (bin) * scale);
	};


	//std::cout << "\n ======= "<< __FUNCTION__ << ":" << __LINE__ << ":" << std::endl;
	//std::cout << " ======= Err hist val and ERROR vector values" << std::endl;
	//std::cout << setw(8) << "bin" << setw(13) << "value" << setw(13) << "Error" 
	//	<< setw(10) << "rat_err"
	//	<< std::endl; 
	TH1 *hist_err_copy = NULL;
	{
		std::string myname = hist_err->GetName() + std::string ("_copy");
		hist_err_copy = dynamic_cast<TH1*>(hist_err->Clone (myname.c_str()));
		for (unsigned bin = 1; bin <= (unsigned) hist_err_copy->GetNbinsX(); ++ bin)
		{
			float value = hist_err_copy->GetBinContent (bin);
			float error = hist_err_copy->GetBinError (bin);
			hist_err_copy->SetBinError (bin, value ? error / value : 0);
			
			//hist_err_copy->SetBinError (bin, value ? Errors.at(bin-1) / value : 0);
			hist_err_copy->SetBinContent (bin, 0);

			//std::cout << setiosflags(ios::fixed) << setprecision(3) 
			//	<< setw(8) << bin << setw(13) << value << setw(13) << Errors.at(bin-1)
			//	<< setw(10) << hist_err_copy->GetBinError (bin)
			 //	<< std::endl; 
		};
	};


	
	//for (unsigned bin = 1; bin <= (unsigned)err->GetNbinsX(); ++ bin)
	
	pho->SetTitle("");
	hist_err_copy->SetTitle("");

	//ytitle << "#frac{Data -- Background}{Background}        Events/"<< phojet->GetBinWidth(1) << " " << yunits;
	const std::string ytitle("(Data #topbar Background) / Background");
	//std::cout << "title= " << title << std::endl;
	hist_err_copy->GetYaxis()->CenterTitle(true);
	hist_err_copy->GetXaxis()->CenterTitle(true);
	hist_err_copy->SetTitleOffset(0.95,"Y");
	hist_err_copy->SetTitleOffset(0.98,"X");
	hist_err_copy->GetYaxis()->SetTitle(ytitle.c_str());
	hist_err_copy->GetXaxis()->SetTitle(xtitle.c_str());
	hist_err_copy->GetYaxis()->SetTitleColor(kBlack);
	hist_err_copy->GetXaxis()->SetTitleColor(kBlack);
	hist_err_copy->SetMinimum (-1.);
	hist_err_copy->SetMaximum (1.);
	
	pho->SetMarkerStyle(8);
	

	if (pad == 0)
	{
		std::cout << red << __FUNCTION__ << ":" << __LINE__ << ":: PAD given is NULL creating new PAD"<< clearatt << std::endl;
		new TCanvas();
	} else pad->cd();

	gStyle->SetOptStat(0);
	gStyle->SetCanvasColor (10);
	//gStyle->SetCanvasBorderSize (0);
	//gStyle->SetCanvasBorderMode (0);
				 
	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (0);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);

	int labelfont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	int titlefont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	hist_err_copy->GetXaxis()->SetLabelFont(iAXIS_LABEL_FONT);
	hist_err_copy->GetYaxis()->SetLabelFont(iAXIS_LABEL_FONT);
	hist_err_copy->GetYaxis()->SetLabelSize(0.055);
	hist_err_copy->GetXaxis()->SetLabelSize(0.055);
	hist_err_copy->GetXaxis()->SetTitleFont(iAXIS_LABEL_FONT);
	hist_err_copy->GetYaxis()->SetTitleFont(iAXIS_LABEL_FONT);
	hist_err_copy->GetXaxis()->SetTitleSize(0.055);
	hist_err_copy->GetYaxis()->SetTitleSize(0.055);
	hist_err_copy->GetXaxis()->CenterTitle(true);
	hist_err_copy->GetYaxis()->CenterTitle(true);
	hist_err_copy->GetXaxis()->SetTitleOffset(0.95);
	hist_err_copy->GetYaxis()->SetTitleOffset(0.95);
	
	//pho->SetMarkerColor(kBlue);
	//pho->SetLineColor(kBlue);
	pho->SetMarkerColor(kBlack);
	pho->SetLineColor(kBlack);
	// phojet->GetXaxis()->SetTitle(title.c_str());
	// phojet->GetXaxis()->CenterTitle(true);
  	//hist_err_copy->SetFillStyle(3002);
  	hist_err_copy->SetFillStyle(3344);
  	//hist_err_copy->SetFillColor(kBlack);
  	hist_err_copy->SetFillColor(9);
	hist_err_copy->SetLineColor(kWhite);
	
  	//TPaveText *tp = new TPaveText(0.4,0.7,0.9,0.95,"NDC");
	//tp->SetLineColor(10);
	//tp->SetTextFont(titlefont);
	//tp->AddText("Shaded Region = Syst + Stat Err on Prediction");


	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridx();
	gPad->SetGridy();


	TLegend *leg = new TLegend (0.55,0.75,0.9,1.0);
	leg->SetTextFont(iAXIS_LABEL_FONT);
	//leg->SetTextFont(42);

	leg->SetTextSize(0.05);
	leg->SetBorderSize (1);
	leg->SetFillColor (10);

	//std::string lbl;
	//if (jets == 1) lbl = "Data (#gamma +#geq 1 Jet)";
	//if (jets == 2) lbl = "Data (#gamma +#geq 2 Jets)";

	leg->AddEntry(pho,"Only Stat Err on Data");
	leg->AddEntry(hist_err_copy,"Syst + Stat on Bg.");
	//leg->AddEntry(pho,"Data");
	//leg->AddEntry(hist_err_copy,"Total Bkgd. Uncertainty");

	hist_err_copy->Draw ("E2");
	pho->Draw("SAME");
	leg->Draw();
	//tp->Draw();

/*	TCanvas *c = dynamic_cast<TCanvas*>(gPad);
	if (c)
	{
		//c->SetFillColor(kRed);
		std::ostringstream str,str1,str2;
		//str << "ratioplot" << jets << "_" << which << ".gif";
		//c->Print (str.str().c_str());
		//str1 << "ratioplot" << jets << "_" << which << ".pdf";
		//c->Print (str1.str().c_str(),"pdf");
		str2 << "ratioplot" << jets << "_" << which << ".eps";
		c->Print (str2.str().c_str(),"eps");

	};
*/
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
double GetMax(const double def, const double up, const double down)
{
	return max(fabs(def - up), fabs(def - down));
}

//-----------------------------------------------------------------------------
//return max of 4 numbers
//-----------------------------------------------------------------------------
double GetMax(const double x1, const double x2, const double x3, const double x4)
{
	return max(max(x1,x2), max(x3,x4));
}


//-----------------------------------------------------------------------------
// derive the uncertainty for the method B, for the reweighing 
// of the sideband photon Et to the total of QCD(30%)+MC PHO(70%)
//-----------------------------------------------------------------------------
void GetSidebandReweighedError(std::vector<float>& ErrVec, const TH1* qcdhist,
						const float qcdscale, const int iSample, const int njets,
						bool debug)
{

	if (iDEBUG > iFUNCTION_NAMES) In(__FUNCTION__);
	assert (qcdhist != NULL && "GetSidebandReweighedError:: qcdhist is NULL!");
	ErrVec.clear();

	//only when using reweighed sideband 100% from Elog#1520
	//this is the error for the reweighed sideband in addition
	//to the choice of sideband to represent fake jets in
	//tight sample.

	
	AtPoint(__FUNCTION__,__LINE__);
	std::cout << "sSIDEBAND_HIST_PATH = " << sSIDEBAND_HIST_PATH << std::endl;
	assert (file_DATA_WGTSIDEBAND_ERR != 0 && "file not found!");
	TH1F* qcdwgthist = dynamic_cast<TH1F*> (file_DATA_WGTSIDEBAND_ERR->Get(sSIDEBAND_HIST_PATH.c_str()));
	AtPoint(__FUNCTION__,__LINE__);
	if (! qcdwgthist){
		std::cout << __FUNCTION__ << ":" << __LINE__ << ": hist not found in the dir" <<std::endl;
		return;
	}
	
	if (bDO_VAR_BINNING)
	{
		qcdwgthist  = (TH1F*) MakeVariableBins (qcdwgthist, false); 
	} else 
	{
		qcdwgthist->Sumw2();
		qcdwgthist->Rebin(iREBIN);
	}



	SubtractEWKfromQCD(qcdwgthist); 
	SubtractCosmicFromQCD(qcdwgthist,	iSample,	njets);
	CheckNegativeBin(qcdwgthist);
	
	AtPoint(__FUNCTION__,__LINE__);
	qcdwgthist->SetLineColor(kBlue);

	std::stringstream title;
	title << qcdwgthist->GetTitle() << ": ERROR :Reweighed (#epsilon) sideband comapared to Reweighed (#epsilon+#sigma) sideband";
	qcdwgthist->SetTitle(title.str().c_str());
	AtPoint(__FUNCTION__,__LINE__);

	//do not normalize the two as they only differ by weights! 
	//normazlise the wgted to nominal
	//qcdwgthist->Scale(qcdhist->Integral()/(double)qcdwgthist->Integral());

	qcdwgthist->Scale(qcdscale);
	AtPoint(__FUNCTION__,__LINE__);
	qcdhist->Print();
	qcdwgthist->Print();
	AtPoint(__FUNCTION__,__LINE__);
	TCanvas *c4 = 0;
	if (debug)
	{
		c4 = new TCanvas();
		c4->Divide(1,2);
		c4->cd(1);
		gPad->SetLogy();
		qcdwgthist->DrawClone();
		qcdhist->DrawClone("same");
	}

	//qcdwgthist->Divide(qcdhist);
	if (debug)
	{
		//new TCanvas();
		c4->cd(2);
		gPad->SetTickx();
		gPad->SetTicky();
		gPad->SetGridx();
		gPad->SetGridy();

		qcdwgthist->DrawClone("P");
	}
	AtPoint(__FUNCTION__,__LINE__);

	//now find the max for each bin
	//std::cout << __FUNCTION__ << ":" << __LINE__ << std::endl;
	//std::cout << setw(5) << "bin" << setw(13) << "val" 
	//	<< setw(13) << "FracErr" << setw(10) << "Err" << std::endl;
	for (int i=1; i <= qcdhist->GetNbinsX(); i++) 
	{
		//double err = 0;
		double dy = 0;
		if (qcdhist->GetNbinsX())
		{
			//dy = fabs(qcdwgthist->GetBinContent(i) - 1);
			//err = dy * qcdhist->GetBinContent(i);
			dy = fabs(qcdwgthist->GetBinContent(i) - qcdhist->GetBinContent(i));
		}
		
		//std::cout << setw(5) << i << setw(13) << qcdhist->GetBinContent(i) 
		//	<< setw(13) << dy << setw(10) << err << std::endl;
		//ErrVec.push_back(err);
		ErrVec.push_back(dy);
	}
	if (iDEBUG > iFUNCTION_NAMES) Out(__FUNCTION__);
}



//*************** Derive uncertainty in normalization of QCD and photon MC
void GetPhoFracError(const double leftover, const float fQCDScale, const float fQCDscale_psigma,
							const float fQCDScale_msigma, 
							std::vector<float>& vErr) 
{
	if (iDEBUG > iFUNCTION_NAMES) In(__FUNCTION__);
	assert (leftover>0 && "GetPhoFracError:: leftover must be >0!");
	
	//in either cases (sidebas is nominal or wgted, I will pick the nominal sideband
	//and photon mc and vary thier relative normalization by fake photon fraction for this
	//error calculation
	TH1F* nominal_qcdjet = (TH1F*) file_DATA_NOMINAL->Get(sSIDEBAND_HIST_PATH.c_str());
	TH1F* nominal_mcphojet = (TH1F*) file_MCPHO->Get(sMCCENTRAL_HIST_PATH.c_str());
	if (nominal_qcdjet == NULL ) { std::cout << __LINE__ << ": qcdjet hist null!" << std::endl; return; }
	if (nominal_mcphojet == NULL ) { std::cout << __LINE__ << ": mcphojet hist null!" << std::endl; return; }
	if (bDO_VAR_BINNING)
	{
		nominal_qcdjet = (TH1F*) MakeVariableBins (nominal_qcdjet, false);
		nominal_mcphojet = (TH1F*) MakeVariableBins (nominal_mcphojet, false);
	} else
	{
		nominal_qcdjet->Sumw2();
		nominal_mcphojet->Sumw2();
		nominal_qcdjet->Rebin(iREBIN);
		nominal_mcphojet->Rebin(iREBIN);
	}

	TH1F* psigma_qcdjet = (TH1F*) nominal_qcdjet->Clone("psigma_qcdjetcopy");	//to scale by f+1sigma
	TH1F* msigma_mcphojet = (TH1F*) nominal_mcphojet->Clone("a_mcphoVaried");      //to scale by f-1sigma

	
	//setup nominal distribution
	//nominal_qcdjet->Scale(fQCDScale);
	if (nominal_qcdjet->Integral("width")>0) nominal_qcdjet->Scale((leftover * fFAKE_PHO_FRAC_DEF) / (double) nominal_qcdjet->Integral("width"));
	if (nominal_mcphojet->Integral("width")>0) nominal_mcphojet->Scale((leftover * (1. - fFAKE_PHO_FRAC_DEF)) / (double) nominal_mcphojet->Integral("width"));
	std::cout << __FUNCTION__ << "::qcd/mc nominal = " << nominal_qcdjet->Integral() << ", " << nominal_mcphojet->Integral() << std::endl;
	nominal_qcdjet->Add(nominal_mcphojet);
		
	//setup fake_frac+sigma distribution
	if (psigma_qcdjet->Integral("width")>0) psigma_qcdjet->Scale((leftover * (fFAKE_PHO_FRAC_PSIG)) / (double) psigma_qcdjet->Integral("width"));
	if (msigma_mcphojet->Integral("width")>0) msigma_mcphojet->Scale((leftover * (1 - fFAKE_PHO_FRAC_PSIG)) / (double) msigma_mcphojet->Integral("width"));
	std::cout << __FUNCTION__ << "::qcd/mc psigma  = " << psigma_qcdjet->Integral() << ", " << msigma_mcphojet->Integral() << std::endl;
	psigma_qcdjet->Add(msigma_mcphojet);

	vErr.clear();
	for (int i=1; i <= nominal_qcdjet->GetNbinsX(); i++)
	{
		float err = fabs(nominal_qcdjet->GetBinContent(i) - psigma_qcdjet->GetBinContent(i));
		//std::cout << "i , err = " << i << ", " << err << std::endl;
		vErr.push_back(err);
	}
	//std::cout << " left over / varied = " << leftover << " / " << nominal_qcdjet->Integral("width") << " / "  << psigma_qcdjet->Integral("width") << std::endl;
//	for (int i=1; i <= vErr.size(); i++)
//	std::cout << "i , err = " << i << ", " << vErr.at(i-1) << std::endl;
	
/*	new TCanvas();
	nominal_qcdjet->DrawCopy();
	psigma_qcdjet->DrawCopy("same");

	gPad->SetEditable(0);
*/
/*		delete nominal_qcdjet;
		delete psigma_qcdjet;
		delete msigma_mcphojet;
		*/

	if (iDEBUG > iFUNCTION_NAMES) Out(__FUNCTION__);
}

 //*************** SUBTRACT HALO/COSMIC/EWK FROM QCD BACKGROUND *******
void SubtractCosmicFromQCD(TH1* qcdhist,
						const int iSample, const int njet)
{
	//the old root files do not have these hists
	//only the new root files have SIDEBAND cosmic hists
	
	TH1F* sideband_cosmic = dynamic_cast<TH1F*> (file_DATA_NOMINAL->Get(sSIDEBAND_COSMIC_HIST_PATH.c_str()));
	assert ( sideband_cosmic != NULL && "sideband_cosmic hist is null!");
	sideband_cosmic->Print();


	if (bDO_VAR_BINNING)
	{
		sideband_cosmic = (TH1F*) MakeVariableBins (sideband_cosmic, false);
	} else 
	{
		sideband_cosmic->Sumw2();
		sideband_cosmic->Rebin(iREBIN);
	}


	float fSidebandCosmicEst = 0;
	if (iSample == iG30)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG30_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG40)
	{
		if (njet == 1)	fSidebandCosmicEst = fG40_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG40_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG30MET20)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30MET20_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG30MET20_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG30MET25)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30MET25_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG30MET25_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG30MET30)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30MET30_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG30MET30_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG30MET40)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30MET40_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG30MET40_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG30MET50)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30MET50_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	fSidebandCosmicEst = fG30MET50_SIDEBAND_COSMIC_EST_2JET;
	} else if (iSample == iG30EXCL1J)
	{
		if (njet == 1)	fSidebandCosmicEst = fG30EXCL1J_SIDEBAND_COSMIC_EST_1JET;
		else if (njet == 2)	
		{
			std::cout << red << "ERROR! No 2jet case for excl njet=1 sample!";
			assert (false);
		}
	} else
	{
		assert (false && "SubtractCosmicFromQCD::Unknown Sample. No cosmic estiamtes found!");
	}


	if (fSidebandCosmicEst == 0)
		cout << red << __FUNCTION__ << ": NO SIDEBAND COSMIC ESTIMATE FOUND FOR THIS SAMPLE! CHECK!" 
			<< clearatt << endl;
	if (sideband_cosmic->Integral("width")>0)
	{
		sideband_cosmic->Scale(fSidebandCosmicEst/sideband_cosmic->Integral("width"));
	}

	std::cout << cyan << "sideband cosmic integral = " << sideband_cosmic->Integral("width") << endl;
	
	qcdhist->Add(sideband_cosmic, -1);

}

 //*************** SUBTRACT DIPHO COMPONENT IN INCLUSIVE PHOTON MC *******
 //I am passing in the phojet mc hist before scaling. I create diphojet
 //hist and scale to photon mc by lum and subtract it from phojet mc.
 //I had to do this separately as diphojet in mainrotine is scaled to
 //data and that is not what I want to use here.
void SubtractDiPhofromPhoMC(TH1* mcpho, const bool debug)
{
	assert (mcpho != NULL && "SubtractDiPhofromPhoMC:: phomc hist ptr is null!");
	TH1F* mcdipho = dynamic_cast<TH1F*> (file_DIPHOMC->Get(sMCCENTRAL_HIST_PATH.c_str()));

	if (bDO_VAR_BINNING)
	{
		mcdipho = (TH1F*) MakeVariableBins(mcdipho, false);
	} else
	{
		mcdipho->Sumw2();
		mcdipho->Rebin(iREBIN);
	}

	
	//now need the proper nomalization to get the correct number of 
	//dipho events in the pho mc sample
	// for this i would scale dipho hist by lum to pho mc lum
	
	//photon mc weighted cross section .0000624357197484 (mb), total events 3444250 (~3.4M)
	//from qcd web page for gen 6 mc samples
	const float fPhoMcLum = 3444250/62435.7;  //(pb)
	const float fDiPhoScale = fPhoMcLum/fLUM_DIPHOMC;
	mcdipho->Scale(fDiPhoScale);
	std::cout << __FUNCTION__ << ":fDiPhoScale = " << fDiPhoScale << " ,dipho integral = " << mcdipho->Integral("width") << std::endl;
	mcpho->Add(mcdipho,-1);
	
	CheckNegativeBin(mcpho);
	
}

 //*************** SUBTRACT EWK COMPONENT IN QCD BACKGROUND *******
 //During this subtraction if there we get a bin that it statitstically 
 //significantly negative we are doing something wrong. not just for that
 //bin, but the whole subtraction.
 //statistically signifcant mean, if the bin value c an be negative as long as
 //the statistical error includes '0' in the range.
void SubtractEWKfromQCD(TH1* qcdhist, const bool debug)
{
	//std::cout << __FUNCTION__ << ":debug  = " << debug << std::endl;
	// this will avoid the double counting
	//I have removed halo/cosmics from both signal and qcd templates,
	// I need to subtract the templates by normalizing them to signal/qcd.
	//but since I have an estimate of the expected number of halo/cosmics 
	// events, I can just subtract the template. halo selection is orthogonal to 
	// photon or qcd selection. cosmic fraction in signal/qcd can be asumed to
	// be same, I think.

	TH1F* sideband_zee = dynamic_cast<TH1F*> (file_ZEEMC->Get(sMCSIDEBAND_HIST_PATH.c_str()));
	TH1F* sideband_zmm = dynamic_cast<TH1F*> (file_ZMMMC->Get(sMCSIDEBAND_HIST_PATH.c_str()));
	TH1F* sideband_ztt = dynamic_cast<TH1F*> (file_ZTTMC->Get(sMCSIDEBAND_HIST_PATH.c_str()));
	TH1F* sideband_wen = dynamic_cast<TH1F*> (file_WENMC->Get(sMCSIDEBAND_HIST_PATH.c_str()));
	TH1F* sideband_wmn = dynamic_cast<TH1F*> (file_WMNMC->Get(sMCSIDEBAND_HIST_PATH.c_str()));
	TH1F* sideband_wtn = dynamic_cast<TH1F*> (file_WTNMC->Get(sMCSIDEBAND_HIST_PATH.c_str()));

	assert ( (sideband_zee != NULL && sideband_zmm != NULL && sideband_ztt != NULL
				&& sideband_wen != NULL && sideband_wmn != NULL && sideband_wtn != NULL)
			&& "one or more of the sideband_EWK hists are null!");
	sideband_zee->Print();

	if (bDO_VAR_BINNING)
	{
		sideband_zee = (TH1F*) MakeVariableBins (sideband_zee, false);
		sideband_zmm = (TH1F*) MakeVariableBins (sideband_zmm, false);
		sideband_ztt = (TH1F*) MakeVariableBins (sideband_ztt, false);
		sideband_wen = (TH1F*) MakeVariableBins (sideband_wen, false);
		sideband_wmn = (TH1F*) MakeVariableBins (sideband_wmn, false);
		sideband_wtn = (TH1F*) MakeVariableBins (sideband_wtn, false);
	} else
	{
		sideband_zee->Sumw2();
		sideband_zmm->Sumw2();
		sideband_ztt->Sumw2();
		sideband_wen->Sumw2();
		sideband_wmn->Sumw2();
		sideband_wtn->Sumw2();

		sideband_zee->Rebin(iREBIN);
		sideband_zmm->Rebin(iREBIN);
		sideband_ztt->Rebin(iREBIN);
		sideband_wen->Rebin(iREBIN);
		sideband_wmn->Rebin(iREBIN);
		sideband_wtn->Rebin(iREBIN);
	}


	//kFAC_MC is to account for all NLO contributions to the CROSS SECTION 
	//(simply a correction to everything that we do not know
	//beyound LO). 
	//Contribution to the differential crossection by the kFAC can change 
	//sign from LO to NLO to NNLO and so on. So it is not safe to assume that
	//the kFac will decrease from LO to NLO.
	//(this sign change happens only in 'differential crosssection'. In total
	// crosssection, all higher order terms gives a POSITIVE contribution)

	if (iDEBUG)
	{
		std::cout << __FUNCTION__ << "::" << __LINE__ << ":: Sideband EWK subtration scales"<< std::endl;
		std::cout << "zeenorm = " << fZEE_SCALE << std::endl;
		std::cout << "ZMMnorm = " << fZMM_SCALE << std::endl;
		std::cout << "ZTTnorm = " << fZTT_SCALE << std::endl;
		std::cout << "wennorm = " << fWEN_SCALE << std::endl;
		std::cout << "WMNnorm = " << fWMN_SCALE << std::endl;
		std::cout << "WTNnorm = " << fWTN_SCALE << std::endl;
	}

	
	sideband_zee->Scale(fZEE_SCALE);
	sideband_zmm->Scale(fZMM_SCALE);
	sideband_ztt->Scale(fZTT_SCALE);
	sideband_wen->Scale(fWEN_SCALE);
	sideband_wmn->Scale(fWMN_SCALE);
	sideband_wtn->Scale(fWTN_SCALE);


	if (iDEBUG)
	{
		std::cout << cyan << "sideband ewk stuff" << std::endl;
		cout << "zee = " << sideband_zee->Integral("width") << endl;
		cout << "zmm = " << sideband_zmm->Integral("width") << endl;
		cout << "ztt = " << sideband_ztt->Integral("width") << endl;
		cout << "wen = " << sideband_wen->Integral("width") << endl;
		cout << "wmn = " << sideband_wmn->Integral("width") << endl;
		cout << "wtn = " << sideband_wtn->Integral("width") << endl;
	}

	if (debug)
	{
		new TCanvas();
		sideband_zee->Draw();
		sideband_zmm->Draw("same");
		sideband_ztt->Draw("same");
		sideband_wen->Draw("same");
		sideband_wmn->Draw("same");
		sideband_wtn->Draw("same");
		gPad->SetEditable(0);
	}
	
	/*CorrectEWKMC(sideband_zee, 0);
	CorrectEWKMC(sideband_zmm, 0);
	CorrectEWKMC(sideband_ztt, 0);
	CorrectEWKMC(sideband_wen, 1);
	CorrectEWKMC(sideband_wmn, 1);
	CorrectEWKMC(sideband_wtn, 1);
	*/
	sideband_zee->Add(sideband_zmm);
	sideband_zee->Add(sideband_ztt);
	sideband_zee->Add(sideband_wen);
	sideband_zee->Add(sideband_wmn);
	sideband_zee->Add(sideband_wtn);

	qcdhist->Add(sideband_zee, -1);

/*	new TCanvas();
	gPad->SetLogy();
	qcdhist->DrawCopy();
	gPad->SetEditable(0);
*/	
	//this is a hack to avoid any bin getting negative values after this subtraction
	//I am going set those bins to zero Dec 15,2009

	// need to do this after the ID error is derived. because when I do this here
	// my QCD ID error is zero for lots of bins? need to talk with someone about this!
	//ZeroOutNegativeBins(qcdjet);
	//ZeroOutNegativeBins(qcdjet_100);

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GetQCD100Err(const int isideband, std::vector<float>& ErrVec, const TH1* qcdhist, 
						const int njets, const int iSample,
						bool debug=false)
{
	if (iDEBUG > iFUNCTION_NAMES) In(__FUNCTION__);
	assert (qcdhist != NULL && "GetQCD100Err:: qcdhist is NULL!");
	ErrVec.clear();
	std::cout << __FUNCTION__ << "qcdhist int = " << qcdhist->Integral("width") << std::endl;

	TH1F* hademhist  = 0; 
	TH1F* isohist    = 0;
	TH1F* trkpthist  = 0;
	TH1F* trkisohist = 0;
	
	if (isideband == iNOMINAL_SIDEBAND)
	{
		hademhist  = (TH1F*) file_DATA_NOMINAL->Get(sSIDEBANDSYST_HADEM_HIST_PATH.c_str());
		isohist    = (TH1F*) file_DATA_NOMINAL->Get(sSIDEBANDSYST_ISO_HIST_PATH.c_str());
		trkpthist  = (TH1F*) file_DATA_NOMINAL->Get(sSIDEBANDSYST_TRKPT_HIST_PATH.c_str());
		trkisohist = (TH1F*) file_DATA_NOMINAL->Get(sSIDEBANDSYST_TRKISO_HIST_PATH.c_str());
   } else if (isideband == iWEIGHTED_SIDEBAND)
	{
		std::cout << "using reweighted stuff " << std::endl;
		hademhist  = (TH1F*) file_DATA_WGTSIDEBAND->Get(sSIDEBANDSYST_HADEM_HIST_PATH.c_str());
		isohist    = (TH1F*) file_DATA_WGTSIDEBAND->Get(sSIDEBANDSYST_ISO_HIST_PATH.c_str());
		trkpthist  = (TH1F*) file_DATA_WGTSIDEBAND->Get(sSIDEBANDSYST_TRKPT_HIST_PATH.c_str());
		trkisohist = (TH1F*) file_DATA_WGTSIDEBAND->Get(sSIDEBANDSYST_TRKISO_HIST_PATH.c_str());
	}

	std::cout << __FUNCTION__ << "hadem int = " << hademhist->Integral("width") << std::endl;
	std::cout << __FUNCTION__ << "iso   int = " << isohist->Integral("width") << std::endl;
	std::cout << __FUNCTION__ << "trkpt int = " << trkpthist->Integral("width") << std::endl;
	std::cout << __FUNCTION__ << "trkiso int= " << trkisohist->Integral("width") << std::endl;
	if (hademhist->Integral("width") < 1 
		 || isohist->Integral("width") < 1
		 || trkpthist->Integral("width") < 1
		 || trkisohist->Integral("width") < 1)
	{
		assert (false && "GetQCD100Err:: Id varied hist/s is empty!");
	}
	
	//hademhist->Print("all");
	assert ( hademhist != 0 && "GetQCD100Err:: hademhist is null");
	assert ( isohist   != 0 && "GetQCD100Err:: isohist is null");
	assert ( trkpthist != 0	 && "GetQCD100Err:: ftrkpt is null");
	assert ( trkisohist!= 0 && "GetQCD100Err:: ftrkiso is null");


	if (bDO_VAR_BINNING)
	{
		hademhist  = (TH1F*) MakeVariableBins (hademhist, false); 
		isohist    = (TH1F*) MakeVariableBins (isohist, false);  
		trkpthist  = (TH1F*) MakeVariableBins (trkpthist, false);  
		trkisohist = (TH1F*) MakeVariableBins (trkisohist, false);  
	} else
	{

		hademhist->Sumw2();
		isohist->Sumw2();
		trkpthist->Sumw2();
		trkisohist->Sumw2();
		hademhist->Rebin(iREBIN);
		isohist->Rebin(iREBIN);
		trkpthist->Rebin(iREBIN);
		trkisohist->Rebin(iREBIN);
	}


	SubtractEWKfromQCD(hademhist); 
	SubtractCosmicFromQCD(hademhist,	iSample,	njets);
	CheckNegativeBin(hademhist);

	SubtractEWKfromQCD(isohist); 
	SubtractCosmicFromQCD(isohist,	iSample,	njets);
	CheckNegativeBin(isohist);

	SubtractEWKfromQCD(trkpthist); 
	SubtractCosmicFromQCD(trkpthist,	iSample,	njets);
	CheckNegativeBin(trkpthist);

	SubtractEWKfromQCD(trkisohist); 
	SubtractCosmicFromQCD(trkisohist,	iSample,	njets);
	CheckNegativeBin(trkisohist);

	
	if (hademhist->Integral("width")>0)	 hademhist->Scale(qcdhist->Integral("width")/(double) hademhist->Integral("width"));
	if (isohist->Integral("width")>0)    isohist->Scale(qcdhist->Integral("width")/(double) isohist->Integral("width"));
	if (trkpthist->Integral("width")>0)  trkpthist->Scale(qcdhist->Integral("width")/(double) trkpthist->Integral("width"));
	if (trkisohist->Integral("width")>0) trkisohist->Scale(qcdhist->Integral("width")/(double) trkisohist->Integral("width"));

	hademhist->Divide(qcdhist);
	isohist->Divide(qcdhist);
	trkpthist->Divide(qcdhist);
	trkisohist->Divide(qcdhist);

	//need to compare these varitions to the nominal sideband
	//so I need the sideband without any subtraction or addition
	//or scaling. - 01-14-2010
	//THIS is true. But when you derive the absolute error
	//from these percentages, you must use the final scaled qcd hist
	//but then what if some component like ewk is already subtracted
	//from the the final! do i derive absoluet error from that?
	//but for now this is ok, as I am not subtracting anything from the
	//sideband!-01-23-2010

	//debug = true;
	if (debug)
	{
		hademhist->SetMaximum(10);
		isohist->SetLineColor(kRed);
		isohist->SetMarkerColor(kRed);
		trkpthist->SetLineColor(kBlue);
		trkpthist->SetMarkerColor(kBlue);
		trkisohist->SetLineColor(kGreen);
		trkisohist->SetMarkerColor(kGreen);
		hademhist->SetMarkerStyle(20);
		isohist->SetMarkerStyle(21);
		trkpthist->SetMarkerStyle(22);
		trkisohist->SetMarkerStyle(23);

		std::stringstream title;
		title << qcdhist->GetTitle() << " : QCD Shape Sytematic by varying #gamma ID cuts;;Scale(QCD/Varied) / QCD" << std::endl;
		hademhist->SetTitle(title.str().c_str());

		new TCanvas();
		gPad->SetGridx();
		gPad->SetGridy();
		hademhist->DrawClone("P");
		isohist->DrawClone("sameP");
		trkpthist->DrawClone("sameP");
		trkisohist->DrawClone("sameP");

		TLegend *leg = new TLegend (0.7,0.6,0.9,0.9);
		leg->SetTextFont(42);
		leg->SetTextSize(0.03);

		leg->AddEntry(isohist,"Iso cut varied");
		leg->AddEntry(hademhist,"HadEm cut varied");
		leg->AddEntry(trkpthist,"TrkPt cut varied");
		leg->AddEntry(trkisohist,"TrkIso cut varied");
		leg->Draw();
		//gPad->SetEditable(0);
	}

	//now find the max for each bin
	//std::cout << __FUNCTION__ << ":" << __LINE__ << std::endl;
	//std::cout << setw(5) << "bin" << setw(13) << "val" 
	//	<< setw(13) << "FracErr" << setw(10) << "Err" << std::endl;
	for (int i=1; i <= qcdhist->GetNbinsX(); i++) 
	{
		float y1 = fabs(hademhist->GetBinContent(i) - 1.);
		float y2 = fabs(isohist->GetBinContent(i) - 1);
		float y3	= fabs(trkpthist->GetBinContent(i) - 1);
		float y4 = fabs(trkisohist->GetBinContent(i) - 1);

		float max = GetMax(y1,y2,y3,y4);
		float err = max * qcdhist->GetBinContent(i);
	//	std::cout << setw(5) << i << setw(13) << qcdhist->GetBinContent(i) 
	//		<< setw(13) << max << setw(10) << err << std::endl;
		ErrVec.push_back(err);
	}

	if (iDEBUG > iFUNCTION_NAMES) Out(__FUNCTION__);
}



void MakeHistLogAndRatioV2 (std::string which, const int jets, 
					const std::string name,const std::string xtitle,
				 	const std::string path, const int QCDerrMtd,
				 	//const float xmin, const float xpoint1, 
					//const float xpoint2, const float xpoint3, 
					//const float xpoint4, const float width1, 
					//const float width2, const float width3, 
					//const float width4, 
					const std::string brname,
					const int iSample,
					const int iSideband=0,   //0=nominal, 1= reweighed
					const std::string sMtdBerr=""
					)
{


	if ( !( iSample == iG30 || iSample == iG40 || iSample == iG30MET20 || iSample == iG30MET25 
				|| iSample == iG30MET30 || iSample == iG30MET40 || iSample == iG30MET50 || iSample == iG30EXCL1J))
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ 
			<< ": Invalid sample selection! returning." << std::endl;
		return;
	}
	if ( !(iSideband == iNOMINAL_SIDEBAND || iSideband == iWEIGHTED_SIDEBAND))
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ 
			<< ": Invalid sideband selection! returning." << std::endl;
		return;
	}

	std::cout << __FUNCTION__ << ": iSample, iSideband = "  << WhichSample(iSample) << ", " << WhichSideband(iSideband) << std::endl;

	//vector to hold systematic unceratinties for individual backgrounds.
	//QCD, PHOTON MC, DIPHO, EWK, COSMIC, HALO
	
	std::vector <double> vJES_phomc ,vJES_diphomc, vJES_ewkmc, vJES_zee, vJES_zmm, vJES_ztt, vJES_wen, vJES_wmn, vJES_wtn;
	std::vector <double> vPhoFrac_qcd, vPhoFrac_phomc;  //fake photon fraction
	std::vector <double> vSideband_qcd; 	//for the choice of sideband to represent QCD
	std::vector <double> vISR_phomc, vFSR_phomc, vPDF_phomc, vQ2_phomc, vAlphaS_phomc;  //I have systematic samples ont  pho mc only
	std::vector <double> vCosmicEst, vHaloEst;
	std::vector <double> vStat_phomc, vStat_diphomc, vStat_ewk, vStat_comic, vStat_halo;





	//initialize some of the global variables
	if (iSample == iG30) sSAMPLE_NAME = sG30;
	else if (iSample == iG40) 		 sSAMPLE_NAME = sG40;
	else if (iSample == iG30MET20) sSAMPLE_NAME = sG30MET20;
	else if (iSample == iG30MET30) sSAMPLE_NAME = sG30MET30;
	else if (iSample == iG30MET50) sSAMPLE_NAME = sG30MET50;
	else if (iSample == iG30MET25) sSAMPLE_NAME = sG30MET25;
	else if (iSample == iG30MET40) sSAMPLE_NAME = sG30MET40;
	else if (iSample == iG30EXCL1J) sSAMPLE_NAME = sG30EXCL1J;

	std::string sSideband("UNKNOWN");
	if (iSideband == iNOMINAL_SIDEBAND) sSideband = "NOMINAL";
	else if (iSideband == iWEIGHTED_SIDEBAND) sSideband = "REWEIGHED";

	std::string sMethod("UNKNOWN");
	if (QCDerrMtd == 1) sMethod = "100% sideband";
	else if (QCDerrMtd == 2) sMethod = "70% MC + 30% sideband";

	
	TH1::AddDirectory(kFALSE);

	//the string 'which' will be used to identify hists with same name
	// for eg: Photon Et has the name EtCorr
	// and so does the Lead jet Et, and so on..
	// if which=="Photon" (with jets==1) -> Photon EtCorr of photon+1jet
	// if which=="Lead Jet" (with jets==1) -> Lead Jet EtCorr of photon+1jet
	// if which=="Photon" (with jets==2) -> Photon EtCorr of photon+2jet
	// if which=="Lead Jet" (with jets==2) -> Lead Jet EtCorr of photon+2jet
	// if which=="Second Lead Jet" (with jets==2) -> Second Lead Jet EtCorr of photon+2jet
	// same for Inv Mass

	//get all hists for nominal sideband from this file. also in the case of reweighted sideband, except the wgted sideband.
	//wgted sideband should be picked form the file_DATA_WGTSIDEBAND
	file_DATA_NOMINAL = new TFile(sDATA_FILE_NOMINAL.c_str());
	if (iSideband == iWEIGHTED_SIDEBAND) file_DATA_WGTSIDEBAND = new TFile(sDATA_FILE_WGTSIDEBAND.c_str());
	file_MCPHO = new TFile(sMCPHO_FILE.c_str());
	file_ZEEMC = new TFile(sMCZEE_FILE.c_str());
	file_ZMMMC = new TFile(sMCZMM_FILE.c_str());
	file_ZTTMC = new TFile(sMCZTT_FILE.c_str());
	file_WENMC = new TFile(sMCWEN_FILE.c_str());
	file_WMNMC = new TFile(sMCWMN_FILE.c_str());
	file_WTNMC = new TFile(sMCWTN_FILE.c_str());
	file_DIPHOMC = new TFile(sMCDIPHO_FILE.c_str());

	if (file_DATA_NOMINAL->IsZombie() ||file_MCPHO->IsZombie() ||
			file_ZEEMC->IsZombie() ||file_ZMMMC->IsZombie() ||
			file_ZTTMC->IsZombie() ||file_WENMC->IsZombie() ||
			file_WMNMC->IsZombie() ||file_WTNMC->IsZombie() || file_DIPHOMC->IsZombie()) 
	{
		std::cout << __FUNCTION__ << ":" << __LINE__  
			<< " a file not found. pl check" << std::endl;
		return;
	}

	if (iSideband == iWEIGHTED_SIDEBAND)  //check if the file for wieghting error exists
	{
		if (file_DATA_WGTSIDEBAND->IsZombie())
		{
			std::cout << __FUNCTION__ << ":" << __LINE__  
				<< " file_DATA_WGTSIDEBAND file not found! Given file name is, " << sDATA_FILE_WGTSIDEBAND << ".  returning!" << std::endl;
			return;
		} else std::cout << "rewgt file found!"<< std::endl;


		file_DATA_WGTSIDEBAND_ERR = new TFile(sDATA_FILE_WGTSIDEBAND_ERR.c_str());
		if (file_DATA_WGTSIDEBAND_ERR->IsZombie())
		{
			std::cout << __FUNCTION__ << ":" << __LINE__  
				<< " file_DATA_WGTSIDEBAND_ERR file not found! Given file name is, " << sDATA_FILE_WGTSIDEBAND_ERR << ".  returning!" << std::endl;
			return;
		} else std::cout << "rewgt file found!"<< std::endl;
	}

	TH1F* phojet = (TH1F*) file_DATA_NOMINAL->Get(sDATA_HIST_PATH.c_str());
	if (! phojet){
		std::cout << "hist not found in the dir" <<std::endl;
		return;
	} else phojet->SetName("phojet");
	//now check if this hist has been filled (or anything is in there ti begin with
	if (phojet->GetEntries()==0)
	{
		std::cout << red << __FUNCTION__ << ":" <<__LINE__ << ":" << name
			<< " hist it empty! please check!";
		phojet->Print();
		std::cout << red << __FUNCTION__ << ":" <<__LINE__ << ":" << name
			<< " hist it empty! please check! returning!" << clearatt << std::endl;
		return;
	}

	TH1F* halojet = (TH1F*) file_DATA_NOMINAL->Get(sHALO_HIST_PATH.c_str());
	if (halojet == NULL ) { std::cout << __LINE__ << ": Halo hist null!" << std::endl; return; }
	else halojet->SetName("halojet");

	TH1F* cosmicjet = (TH1F*) file_DATA_NOMINAL->Get(sCOSMIC_HIST_PATH.c_str());
	if (cosmicjet == NULL ) { std::cout << __LINE__ << ": cosmic hist null!" << std::endl; return; }
	else cosmicjet->SetName("cosmicjet");

	//for QCD+MC combined method and holds the nominal sideband in 100% sideband (NJET and REWEIGHTED mtd)
	//100% sideband (deafult method njet plot) and all plots in reweighed mtd)
	TH1F* qcdjet = 0;
	if (iSideband == iNOMINAL_SIDEBAND)
	{
		std::cout << green << "using nominal sideband ..." << clearatt << std::endl;
		qcdjet = (TH1F*) file_DATA_NOMINAL->Get(sSIDEBAND_HIST_PATH.c_str());
	} else if (iSideband == iWEIGHTED_SIDEBAND) 
	{
		std::cout << green << "using weighted sideband ..." << clearatt << std::endl;
		qcdjet = (TH1F*) file_DATA_WGTSIDEBAND->Get(sSIDEBAND_HIST_PATH.c_str()); 
	}
	if (qcdjet == NULL ) { std::cout << __LINE__ << ": qcdjet hist null!" << std::endl; return; }
	else qcdjet->SetName("qcdjet");


	//MC HISTS: jesup= jesup && emup : jesdown = jesdown && emdown

	TH1F* mcphojet = (TH1F*) file_MCPHO->Get(sMCCENTRAL_HIST_PATH.c_str());
	if (mcphojet == NULL ) { std::cout << __LINE__ << ": mcpho hist null!" << std::endl; return; }
	else mcphojet->SetName("mcphojet");
	TH1F* mcphojetJESUP = (TH1F*) file_MCPHO->Get(sMCEMJESUP_HIST_PATH.c_str());
	if (mcphojetJESUP == NULL ) { std::cout << __LINE__ << ": mcphoJESUP hist null!" << std::endl; return; }
	else mcphojetJESUP->SetName("mcphojetJESUP");
	TH1F* mcphojetJESDOWN = (TH1F*) file_MCPHO->Get(sMCEMJESDOWN_HIST_PATH.c_str());
	if (mcphojetJESDOWN == NULL ) { std::cout << __LINE__ << ": mcphoJESDOWN hist null!" << std::endl; return; }
	else mcphojetJESDOWN->SetName("mcphojetJESDOWN");

	TH1F* zeejet = (TH1F*) file_ZEEMC->Get(sMCCENTRAL_HIST_PATH.c_str());
	if (zeejet == NULL ) { std::cout << __LINE__ << ": zee hist null!" << std::endl; return; }
	else zeejet->SetName("zeejet");
	TH1F* zeejetJESUP = (TH1F*) file_ZEEMC->Get(sMCEMJESUP_HIST_PATH.c_str());
	if (zeejetJESUP == NULL ) { std::cout << __LINE__ << ": zeeJESUP hist null!" << std::endl; return; }
	TH1F* zeejetJESDOWN = (TH1F*)  file_ZEEMC->Get(sMCEMJESDOWN_HIST_PATH.c_str());
	if (zeejetJESDOWN == NULL ) { std::cout << __LINE__ << ":zeeJESDOWN hist null!" << std::endl; return; }


	TH1F* zmmjet = (TH1F*) file_ZMMMC->Get(sMCCENTRAL_HIST_PATH.c_str());
	if (zmmjet == NULL ) { std::cout << __LINE__ << ": zmm hist null!" << std::endl; return; }
	else zmmjet->SetName("zmmjet");
	TH1F* zmmjetJESUP = (TH1F*) file_ZMMMC->Get(sMCEMJESUP_HIST_PATH.c_str());
	if (zmmjetJESUP == NULL ) { std::cout << __LINE__ << ": zmmJESUP hist null!" << std::endl; return; }
	TH1F* zmmjetJESDOWN = (TH1F*) file_ZMMMC->Get(sMCEMJESDOWN_HIST_PATH.c_str());
	if (zmmjetJESDOWN == NULL ) { std::cout << __LINE__ << ": zmmJESDOWN hist null!" << std::endl; return; }

	TH1F* zttjet = (TH1F*) file_ZTTMC->Get(sMCCENTRAL_HIST_PATH.c_str());
	if (zttjet == NULL ) { std::cout << __LINE__ << ": ztt hist null!" << std::endl; return; }
	else zttjet->SetName("zttjet");
	TH1F* zttjetJESUP = (TH1F*) file_ZTTMC->Get(sMCEMJESUP_HIST_PATH.c_str());
	if (zttjetJESUP == NULL ) { std::cout << __LINE__ << ": zttJESUP hist null!" << std::endl; return; }
	TH1F* zttjetJESDOWN = (TH1F*) file_ZTTMC->Get(sMCEMJESDOWN_HIST_PATH.c_str());
	if (zttjetJESDOWN == NULL ) { std::cout << __LINE__ << ": zttJESDOWN hist null!" << std::endl; return; }

	TH1F* wenjet = (TH1F*) file_WENMC->Get(sMCCENTRAL_HIST_PATH.c_str());
	if (wenjet == NULL ) { std::cout << __LINE__ << ": wen hist null!" << std::endl; return; }
	else wenjet->SetName("wenjet");
	TH1F* wenjetJESUP = (TH1F*) file_WENMC->Get(sMCEMJESUP_HIST_PATH.c_str());
	if (wenjetJESUP == NULL ) { std::cout << __LINE__ << ": wenJESUP hist null!" << std::endl; return; }
	TH1F* wenjetJESDOWN = (TH1F*) file_WENMC->Get(sMCEMJESDOWN_HIST_PATH.c_str());
	if (wenjetJESDOWN == NULL ) { std::cout << __LINE__ << ": wenJESDOWN hist null!" << std::endl; return; }

	TH1F* wmnjet = (TH1F*) file_WMNMC->Get(sMCCENTRAL_HIST_PATH.c_str());
	if (wmnjet == NULL ) { std::cout << __LINE__ << ": wmn hist null!" << std::endl; return; }
	else wmnjet->SetName("wmnjet");
	TH1F* wmnjetJESUP = (TH1F*) file_WMNMC->Get(sMCEMJESUP_HIST_PATH.c_str());
	if (wmnjetJESUP == NULL ) { std::cout << __LINE__ << ": wmnJESUP hist null!" << std::endl; return; }
	TH1F* wmnjetJESDOWN = (TH1F*) file_WMNMC->Get(sMCEMJESDOWN_HIST_PATH.c_str());
	if (wmnjetJESDOWN == NULL ) { std::cout << __LINE__ << ": wmnJESDOWN hist null!" << std::endl; return; }

	TH1F* wtnjet = (TH1F*) file_WTNMC->Get(sMCCENTRAL_HIST_PATH.c_str());
	if (wtnjet == NULL ) { std::cout << __LINE__ << ": wtn hist null!" << std::endl; return; }
	else wtnjet->SetName("wtnjet");
	TH1F* wtnjetJESUP = (TH1F*) file_WTNMC->Get(sMCEMJESUP_HIST_PATH.c_str());
	if (wtnjetJESUP == NULL ) { std::cout << __LINE__ << ": wtnJESUP hist null!" << std::endl; return; }
	TH1F* wtnjetJESDOWN = (TH1F*) file_WTNMC->Get(sMCEMJESDOWN_HIST_PATH.c_str());
	if (wtnjetJESDOWN == NULL ) { std::cout << __LINE__ << ": wtnJESDWON hist null!" << std::endl; return; }

	TH1F* diphojet = dynamic_cast<TH1F*> (file_DIPHOMC->Get(sMCCENTRAL_HIST_PATH.c_str()));
	assert (diphojet != NULL && "diphojet central hist not found !");
	diphojet->SetName("diphojet");
	TH1F* diphojetJESUP = dynamic_cast<TH1F*> (file_DIPHOMC->Get(sMCEMJESUP_HIST_PATH.c_str()));
	assert (diphojetJESUP != NULL && "diphojet jesup hist not found !");
	TH1F* diphojetJESDOWN = dynamic_cast<TH1F*> (file_DIPHOMC->Get(sMCEMJESDOWN_HIST_PATH.c_str()));
	assert (diphojetJESDOWN != NULL && "diphojet jesdown hist not found !");

	cout << red << __LINE__ << ":before rebin " << endl;
	phojet->Print();
	double dInt1 =0;
	for (int bin=1; bin<phojet->GetNbinsX(); ++bin) dInt1 += phojet->GetBinContent(bin)*phojet->GetBinWidth(bin);
	std::cout << "\t Integral = " << phojet->Integral("width")  << ", measured = " << dInt1 << std::endl;
	CheckIntegral(phojet);
	CheckIntegral(mcphojet);
	CheckIntegral(qcdjet);
	CheckIntegral(diphojet);
	CheckIntegral(zeejet);
	CheckIntegral(zmmjet);
	CheckIntegral(zttjet);
	CheckIntegral(wenjet);
	CheckIntegral(wmnjet);
	CheckIntegral(wtnjet);
	CheckIntegral(cosmicjet);
	CheckIntegral(halojet);

	mcphojet->Print();
	qcdjet->Print();
	halojet->Print();
	cosmicjet->Print();
	zeejet->Print();
	zmmjet->Print();
	zttjet->Print();
	diphojet->Print();
	cout << clearatt << endl;

	//this is to get the final predictions with systematics
	const float fOrigBinWidth = phojet->GetBinWidth(1);
	std::cout << blue << "fOrigBinWidth = " << fOrigBinWidth << std::endl;
	double fFirst5binVal = 0;
	int nbins = 0;
	for (int bin = 1; bin <=50; ++bin)
	{
		//if (phojet->GetBinContent(bin)>0)
		{
		fFirst5binVal += phojet->GetBinContent(bin);
		std::cout << "bin[lo edge] , sum = " << bin << "[" 
				<< phojet->GetBinLowEdge(bin) << "]"
				<< std::setprecision(10)<< fFirst5binVal << std::endl;
		++nbins;
		//if (nbins>=5) break;
		}
	}
	

	if (bDO_VAR_BINNING)
	{
		phojet = (TH1F*) MakeVariableBins (phojet, false);
		halojet = (TH1F*) MakeVariableBins (halojet, false);
		zeejet = (TH1F*) MakeVariableBins (zeejet, false);
		zmmjet = (TH1F*) MakeVariableBins (zmmjet, false);
		zttjet = (TH1F*) MakeVariableBins (zttjet, false);
		wenjet = (TH1F*) MakeVariableBins (wenjet, false);
		wmnjet = (TH1F*) MakeVariableBins (wmnjet, false);
		wtnjet = (TH1F*) MakeVariableBins (wtnjet, false);
		cosmicjet = (TH1F*) MakeVariableBins (cosmicjet, false);
		qcdjet = (TH1F*) MakeVariableBins (qcdjet, false);
		mcphojet = (TH1F*) MakeVariableBins (mcphojet, false);
		diphojet = (TH1F*) MakeVariableBins (diphojet,false);

		zeejetJESUP = (TH1F*) MakeVariableBins (zeejetJESUP, false);
		zmmjetJESUP = (TH1F*) MakeVariableBins (zmmjetJESUP, false);
		zttjetJESUP = (TH1F*) MakeVariableBins (zttjetJESUP, false);
		wenjetJESUP = (TH1F*) MakeVariableBins (wenjetJESUP, false);
		wmnjetJESUP = (TH1F*) MakeVariableBins (wmnjetJESUP, false);
		wtnjetJESUP = (TH1F*) MakeVariableBins (wtnjetJESUP, false);
		mcphojetJESUP = (TH1F*) MakeVariableBins (mcphojetJESUP, false);
		diphojetJESUP = (TH1F*) MakeVariableBins (diphojetJESUP,false);

		zeejetJESDOWN = (TH1F*) MakeVariableBins (zeejetJESDOWN, false);
		zmmjetJESDOWN = (TH1F*) MakeVariableBins (zmmjetJESDOWN, false);
		zttjetJESDOWN = (TH1F*) MakeVariableBins (zttjetJESDOWN, false);
		wenjetJESDOWN = (TH1F*) MakeVariableBins (wenjetJESDOWN, false);
		wmnjetJESDOWN = (TH1F*) MakeVariableBins (wmnjetJESDOWN, false);
		wtnjetJESDOWN = (TH1F*) MakeVariableBins (wtnjetJESDOWN, false);
		mcphojetJESDOWN = (TH1F*) MakeVariableBins (mcphojetJESDOWN,false);
		diphojetJESDOWN = (TH1F*) MakeVariableBins (diphojetJESDOWN,false);

	} else
	{
		phojet->Sumw2();
		halojet->Sumw2();
		zeejet->Sumw2();
		zmmjet->Sumw2();
		zttjet->Sumw2();
		wenjet->Sumw2();
		wmnjet->Sumw2();
		wtnjet->Sumw2();
		cosmicjet->Sumw2();
		qcdjet->Sumw2();
		mcphojet->Sumw2();

		zeejetJESUP->Sumw2();
		zmmjetJESUP->Sumw2();
		zttjetJESUP->Sumw2();
		wenjetJESUP->Sumw2();
		wmnjetJESUP->Sumw2();
		wtnjetJESUP->Sumw2();
		mcphojetJESUP->Sumw2();

		zeejetJESDOWN->Sumw2();
		zmmjetJESDOWN->Sumw2();
		zttjetJESDOWN->Sumw2();
		wenjetJESDOWN->Sumw2();
		wmnjetJESDOWN->Sumw2();
		wtnjetJESDOWN->Sumw2();
		mcphojetJESDOWN->Sumw2();


		phojet->Rebin(iREBIN);
		halojet->Rebin(iREBIN);
		zeejet->Rebin(iREBIN);
		zmmjet->Rebin(iREBIN);
		zttjet->Rebin(iREBIN);
		wenjet->Rebin(iREBIN);
		wmnjet->Rebin(iREBIN);
		wtnjet->Rebin(iREBIN);
		cosmicjet->Rebin(iREBIN);
		qcdjet->Rebin(iREBIN);
		mcphojet->Rebin(iREBIN);

		zeejetJESUP->Rebin(iREBIN);
		zmmjetJESUP->Rebin(iREBIN);
		zttjetJESUP->Rebin(iREBIN);
		wenjetJESUP->Rebin(iREBIN);
		wmnjetJESUP->Rebin(iREBIN);
		wtnjetJESUP->Rebin(iREBIN);
		mcphojetJESUP->Rebin(iREBIN);

		zeejetJESDOWN->Rebin(iREBIN);
		zmmjetJESDOWN->Rebin(iREBIN);
		zttjetJESDOWN->Rebin(iREBIN);
		wenjetJESDOWN->Rebin(iREBIN);
		wmnjetJESDOWN->Rebin(iREBIN);
		wtnjetJESDOWN->Rebin(iREBIN);
		mcphojetJESDOWN->Rebin(iREBIN);

		diphojet->Sumw2();
		diphojetJESUP->Sumw2();
		diphojetJESDOWN->Sumw2();

		diphojet->Rebin(iREBIN);
		diphojetJESUP->Rebin(iREBIN);
		diphojetJESDOWN->Rebin(iREBIN);
	}

	cout << green << __LINE__ << ": after rebin " << endl;
	double sum =0;
	for (int bin = 1; bin <=5; ++bin)
	{
		//if (phojet->GetBinContent(bin)>0)
		{
			sum += phojet->GetBinContent(bin);
		std::cout << "bin[lo edge] , val, sum = " << bin << "[" 
				<< phojet->GetBinLowEdge(bin) << "]\t " << phojet->GetBinContent(bin)  
				<< std::setprecision(10)<< sum << std::endl;
		}
	}


	phojet->Print();
	double dInt2 =0;
	for (int bin=1; bin<phojet->GetNbinsX(); ++bin) dInt2 += phojet->GetBinContent(bin)*phojet->GetBinWidth(bin);
	std::cout << "\t Integral = " << phojet->Integral("width")  << ", measured = " << dInt2 << std::endl;
	mcphojet->Print();
	qcdjet->Print();
	halojet->Print();
	cosmicjet->Print();
	zeejet->Print();
	zmmjet->Print();
	zttjet->Print();
	cout << clearatt << endl;


	/************** FAKE PHOTON FRACTION ************************************/
	//this is used when the combination of QCD+PHO MC is used
	// this is the amount of fake photons in the signal we select(jets faking photon) which is ~30%

	SetFakePhotonFraction(iSample, jets);

	std::string qcd_str, mc_str;
	std::ostringstream qcdnum, mcnum;
	qcdnum << "QCD (#gamma sideband, " << fFAKE_PHO_FRAC_DEF<< " of data)";
	//mcnum << "#gamma MC (+MB)" << (1 - fFAKE_PHO_FRAC_DEF) << " of data)";
	mcnum << "#gamma MC (" << (1 - fFAKE_PHO_FRAC_DEF) << " of data)";
	qcd_str = qcdnum.str();		// to be used in the legend of the final plot
	mc_str = mcnum.str();

	//this part is only to generate reweight functions
   //MakePhoEtWeightsForSideband(qcdjet,mcphojet,fFAKE_PHO_FRAC_DEF, fFAKE_PHO_FRAC_SIGMA, brname, jets);
	//return;



	/************************************************************************/
	//******************* NORMALIZING *************************************
	float haloEst = 0, cosmicEst = 0;
	GetHaloCosmicEstimates(iSample, jets, haloEst, cosmicEst);

	if (halojet->Integral("width")>0) halojet->Scale (haloEst /(double) halojet->Integral("width"));		//need to use "width" as I am rebinning these hists at the begining
	if (cosmicjet->Integral("width")>0) cosmicjet->Scale (cosmicEst/(double) cosmicjet->Integral("width"));

   TH1* dipholumsyst = (TH1*) diphojet->Clone("dipholumcopy");
	diphojet->Scale(fDIPHOMC_SCALE);
	diphojetJESUP->Scale (fDIPHOMC_SCALE);
	diphojetJESDOWN->Scale (fDIPHOMC_SCALE);


	//makes copies for luminsoity systematic error calculations
	TH1* zeelumsyst = (TH1*) zeejet->Clone("zeelumcopy");
	TH1* zmmlumsyst = (TH1*) zmmjet->Clone("zmmlumcopy");
	TH1* zttlumsyst = (TH1*) zttjet->Clone("zttlumcopy");
	TH1* wenlumsyst = (TH1*) wenjet->Clone("wenlumcopy");
	TH1* wmnlumsyst = (TH1*) wmnjet->Clone("wmnlumcopy");
	TH1* wtnlumsyst = (TH1*) wtnjet->Clone("wtnlumcopy");


	zeejet->Scale (fZEE_SCALE);
	zmmjet->Scale (fZMM_SCALE);
	zttjet->Scale (fZTT_SCALE);
	wenjet->Scale (fWEN_SCALE);
	wmnjet->Scale (fWMN_SCALE);
	wtnjet->Scale (fWTN_SCALE);


	SubtractEWKfromQCD(qcdjet); 
	SubtractCosmicFromQCD(qcdjet,	iSample,	jets);
	CheckNegativeBin(qcdjet);
	
	std::cout << green << __FUNCTION__ << ":dipho int. = " << diphojet->Integral("width") << std::endl;
	
	SubtractDiPhofromPhoMC(mcphojet);
	CheckNegativeBin(mcphojet);
	
	/*cout << red << __LINE__ << ": before leftover calculation " << endl;
	phojet->Print();
	std::cout << "\t integral = " << phojet->Integral("width")<< std::endl;
	mcphojet->Print();
	std::cout << "\t integral = " << mcphojet->Integral("width")<< std::endl;
	qcdjet->Print();
	std::cout << "\t integral = " << qcdjet->Integral("width")<< std::endl;
	halojet->Print();
	std::cout << "\t integral = " << halojet->Integral("width")<< std::endl;
	cosmicjet->Print();
	std::cout << "\t integral = " << cosmicjet->Integral("width")<< std::endl;
	zeejet->Print();
	std::cout << "\t integral = " << zeejet->Integral("width")<< std::endl;
	zmmjet->Print();
	std::cout << "\t integral = " << zmmjet->Integral("width")<< std::endl;
	zttjet->Print();
	std::cout << "\t integral = " << zttjet->Integral("width")<< std::endl;
	wenjet->Print();
	std::cout << "\t integral = " << wenjet->Integral("width")<< std::endl;
	wmnjet->Print();
	std::cout << "\t integral = " << wmnjet->Integral("width")<< std::endl;
	wtnjet->Print();
	std::cout << "\t integral = " << wtnjet->Integral("width")<< std::endl;
	cout << clearatt << endl;
  */
	double leftover = phojet->Integral("width") - halojet->Integral("width") - cosmicjet->Integral("width") 
		- zeejet->Integral("width") - zmmjet->Integral("width") - zttjet->Integral("width")
		- wenjet->Integral("width") - wmnjet->Integral("width") - wtnjet->Integral("width")
		- diphojet->Integral("width");

	std::cout << red << "leftover = " << leftover  << clearatt << std::endl;
	//this is needed in MET>50 plots where the Total bakcground is more than data
	//even before pho mc and QCD is added.
	if (leftover < 0) 
	{
		leftover =0;
		std::cout << red << "Forced leftover to be zero!" << clearatt << std::endl;
	}

	
	float fQCDScale = 0;
	float fQCDScale_psigma=0, fQCDScale_msigma=0;  //needed for PhoFraction error derivation only
	//use 100% sideband for reweighted and method A njet
	//if (iSideband == iNOMINAL_SIDEBAND && which == "NJet")
	//{
	//	fQCDScale = leftover / qcdjet->Integral("width");
	//	qcdjet->Scale(fQCDScale);
	//	mcphojet->Scale(0); // scale this to zer as there is no photon mc contribution in this case.
	//	mcphojetJESUP->Scale(0); // scale this to zer as there is no photon mc contribution in this case.
	//	mcphojetJESDOWN->Scale(0); // scale this to zer as there is no photon mc contribution in this case.
//
//	} else if (iSideband == iWEIGHTED_SIDEBAND) 
	if (iSideband == iWEIGHTED_SIDEBAND) 
	{

		fQCDScale = leftover / qcdjet->Integral("width");
		qcdjet->Scale(fQCDScale);
		//photon MC needs to be properly scaled to get correct ISR/FSR/PDF/Q2/ALPHA_S errors 
		//for both nominal and weighted sideband methods. But do not add MC PHO as a background 
		//in weighted sideband method.
		if (mcphojet->Integral("width")>0)        mcphojet->Scale((leftover * (1 - fFAKE_PHO_FRAC_DEF)) / mcphojet->Integral("width"));
		if (mcphojetJESUP->Integral("width")>0)   mcphojetJESUP->Scale((leftover * (1 - fFAKE_PHO_FRAC_DEF))/mcphojetJESUP->Integral("width"));
		if (mcphojetJESDOWN->Integral("width")>0) mcphojetJESDOWN->Scale((leftover * (1 - fFAKE_PHO_FRAC_DEF))/mcphojetJESDOWN->Integral("width"));
	} else {  //all plots except njet in Method A

		fQCDScale = (leftover * fFAKE_PHO_FRAC_DEF) / qcdjet->Integral("width");
		fQCDScale_psigma = (leftover * (fFAKE_PHO_FRAC_DEF + fFAKE_PHO_FRAC_SIGMA)) / qcdjet->Integral("width");
		fQCDScale_msigma = (leftover * (fFAKE_PHO_FRAC_DEF - fFAKE_PHO_FRAC_SIGMA)) / qcdjet->Integral("width");
		qcdjet->Scale(fQCDScale);
		//photon MC needs to be properly scaled to get correct ISR/FSR/PDF/Q2/ALPHA_S errors 
		//for both nominal and weighted sideband methods. But do not add MC PHO as a background 
		//in weighted sideband method.
		if (mcphojet->Integral("width")>0)        mcphojet->Scale((leftover * (1 - fFAKE_PHO_FRAC_DEF)) / mcphojet->Integral("width"));
		if (mcphojetJESUP->Integral("width")>0)   mcphojetJESUP->Scale((leftover * (1 - fFAKE_PHO_FRAC_DEF))/mcphojetJESUP->Integral("width"));
		if (mcphojetJESDOWN->Integral("width")>0) mcphojetJESDOWN->Scale((leftover * (1 - fFAKE_PHO_FRAC_DEF))/mcphojetJESDOWN->Integral("width"));
	}

	std::cout << __LINE__ << "AFTER SCALING qcdjet     = " << qcdjet->Integral("width") << std::endl;

	std::cout << " integral data / background = " << phojet->Integral("width") << " / ";

	double dTotBckg =  halojet->Integral() + cosmicjet->Integral()
		+ zeejet->Integral() + zmmjet->Integral() + zttjet->Integral()
		+ wenjet->Integral() + wmnjet->Integral() + wtnjet->Integral() 
		+ diphojet->Integral()
		+ qcdjet->Integral();
	// pho mc is only for method A.
	// But not for njet in method A.
	if (iSideband == iNOMINAL_SIDEBAND && which != "NJet")
	{
		dTotBckg += mcphojet->Integral();
	}
	std::cout << dTotBckg << std::endl;

	// JES 
	
	zeejetJESUP->Scale (fZEE_SCALE);
	zmmjetJESUP->Scale (fZMM_SCALE);
	zttjetJESUP->Scale (fZTT_SCALE);
	wenjetJESUP->Scale (fWEN_SCALE);
	wmnjetJESUP->Scale (fWMN_SCALE);
	wtnjetJESUP->Scale (fWTN_SCALE);

	zeejetJESDOWN->Scale (fZEE_SCALE);
	zmmjetJESDOWN->Scale (fZMM_SCALE);
	zttjetJESDOWN->Scale (fZTT_SCALE);
	wenjetJESDOWN->Scale (fWEN_SCALE);
	wmnjetJESDOWN->Scale (fWMN_SCALE);
	wtnjetJESDOWN->Scale (fWTN_SCALE);

	

	//***************END  NORMALIZING *************************************


	std::cout << "\n" << std::endl; 
	std::cout << green << __LINE__ << ":: ALL  HISTS " << std::endl;
	std::cout << __LINE__ << ":: phojet Integral     =  " << phojet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: halojet Integral    =  " << halojet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: cosmicjet Integral  =  " << cosmicjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: qcdjet Integral     =  " << qcdjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: mcphojet Integral   =  " << mcphojet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: zeejet Integral     =  " << zeejet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: zmmjet Integral     =  " << zmmjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: zttjet Integral     =  " << zttjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: wenjet Integral     =  " << wenjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: wmnjet Integral     =  " << wmnjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: wtnjet Integral     =  " << wtnjet->Integral("width") << std::endl;
	std::cout << __LINE__ << ":: diphojet Integral   =  " << diphojet->Integral("width") << clearatt << std::endl;



	//****************** SYSTEMATICS *************************************

	std::vector<float> cosmicErr;   //relative error for each bin is stored here
	std::vector<float> haloErr;

	GetCosmicErr(cosmicErr, cosmicjet, false);
	GetHaloErr(haloErr, halojet);

	const float fact = 1.0/1.06;

	const float diPhoLumErr = fDIPHOMC_SCALE * fact;
	dipholumsyst->Scale(diPhoLumErr);
	std::vector <float> vDiPhoLumErr;

	for (int bin=1;bin<= dipholumsyst->GetNbinsX(); ++bin)
	{
		vDiPhoLumErr.push_back(fabs(diphojet->GetBinContent(bin) - dipholumsyst->GetBinContent(bin)));
	}

	const float zeeLumErr = fZEE_SCALE * fact;
	const float zmmLumErr = fZMM_SCALE * fact;
	const float zttLumErr = fZTT_SCALE * fact;
	const float wenLumErr = fWEN_SCALE * fact;
	const float wmnLumErr = fWMN_SCALE * fact;
	const float wtnLumErr = fWTN_SCALE* fact;

	zeelumsyst->Scale(zeeLumErr);
	zmmlumsyst->Scale(zmmLumErr);
	zttlumsyst->Scale(zttLumErr);
	wenlumsyst->Scale(wenLumErr);
	wmnlumsyst->Scale(wmnLumErr);
	wtnlumsyst->Scale(wtnLumErr);

	std::vector<float> vZeeLumErr,vZmmLumErr, vZttLumErr, vWenLumErr, vWmnLumErr, vWtnLumErr;

	for (int bin=1;bin<= zeelumsyst->GetNbinsX(); ++bin)
	{
		vZeeLumErr.push_back(fabs(zeejet->GetBinContent(bin) - zeelumsyst->GetBinContent(bin)));
		vZmmLumErr.push_back(fabs(zmmjet->GetBinContent(bin) - zmmlumsyst->GetBinContent(bin)));
		vZttLumErr.push_back(fabs(zttjet->GetBinContent(bin) - zttlumsyst->GetBinContent(bin)));
		vWenLumErr.push_back(fabs(wenjet->GetBinContent(bin) - wenlumsyst->GetBinContent(bin)));
		vWmnLumErr.push_back(fabs(wmnjet->GetBinContent(bin) - wmnlumsyst->GetBinContent(bin)));
		vWtnLumErr.push_back(fabs(wtnjet->GetBinContent(bin) - wtnlumsyst->GetBinContent(bin)));
	}




	// JES and QCD mixture 
	// here I am trying to get two things done at once. 1. JES syst 2. fake pho fraction syste using pho mc and qcd

	TH1F *hist_err = (TH1F*) halojet->Clone("hist_err");					//central values and use for SYSTEMATIC ERROR BAND
	// this must be identical to the final stacked plot
	hist_err->Add(zeejet);
	hist_err->Add(zmmjet);
	hist_err->Add(zttjet);
	hist_err->Add(wmnjet);
	hist_err->Add(wenjet);
	hist_err->Add(wtnjet);
	hist_err->Add(cosmicjet);
	hist_err->Add(diphojet);

	std::vector<float> qcdmcMixErr;		//this will hold the systematics for either QCD mthd 1 or 2 for a given case
	//to get systematics when 100% QCD is used
	std::vector<float> IdErr, vSidebandRwtErr;
	//init this so no prob for QCDerrMtd==1
	for (int bin =1; bin <= hist_err->GetNbinsX(); ++bin) 
	{
		qcdmcMixErr.push_back(0);
		vSidebandRwtErr.push_back(0);
	}

	//QCD/MC mix error should be calcualted for both deafault method and reweighed method
		GetPhoFracError(leftover, fQCDScale, fQCDScale_psigma, fQCDScale_msigma,qcdmcMixErr);



	hist_err->Add(qcdjet);
	if (QCDerrMtd == 2)  //_____________________ plots that use 30%70%
	{
		hist_err->Add(mcphojet);
	}

	std::string abspath("Hist/"+path);

	GetQCD100Err(iSideband, IdErr, qcdjet, jets, iSample);
	if (iSideband == iWEIGHTED_SIDEBAND)
	{
		GetSidebandReweighedError(vSidebandRwtErr, qcdjet, fQCDScale, iSample, jets);
	}

	// NOW GET THE SYSTEMATICS FROM JES BY comparing central values with the JES UP/DOWN
	// USE THE HIST USED TO PUT THE ERROR BAND AS THE CENTRAL VALUE

	//to get JES error for all the MC samples, they must be added separately, centra,jesup, and jesdown,
	//and then take the maximum difference. this must be done after each hist is correctly scaled.
	//
	//<1> njet plot in method A: should EXCLUDE the pho mc contribution


	std::vector<float> JESerr;
	for (int bin=1; bin <= zeejet->GetNbinsX(); bin++)
	{
		double jes_cen_ewk = wenjet->GetBinContent(bin) + wmnjet->GetBinContent(bin) + wtnjet->GetBinContent(bin)
			+ zeejet->GetBinContent(bin) + zmmjet->GetBinContent(bin) + zttjet->GetBinContent(bin);

		double jes_up_ewk = wenjetJESUP->GetBinContent(bin) + wmnjetJESUP->GetBinContent(bin) + wtnjetJESUP->GetBinContent(bin)
			+ zeejetJESUP->GetBinContent(bin) + zmmjetJESUP->GetBinContent(bin) + zttjetJESUP->GetBinContent(bin);

		double jes_down_ewk = wenjetJESDOWN->GetBinContent(bin) + wmnjetJESDOWN->GetBinContent(bin) + wtnjetJESDOWN->GetBinContent(bin)
			+ zeejetJESDOWN->GetBinContent(bin) + zmmjetJESDOWN->GetBinContent(bin) + zttjetJESDOWN->GetBinContent(bin);


		double jes_cen = jes_cen_ewk + diphojet->GetBinContent(bin) + mcphojet->GetBinContent(bin);
		double jes_up = jes_up_ewk + diphojetJESUP->GetBinContent(bin) + mcphojetJESUP->GetBinContent(bin);
		double jes_down = jes_down_ewk + diphojetJESDOWN->GetBinContent(bin) + mcphojetJESDOWN->GetBinContent(bin);

		const float err = GetMax(jes_cen, jes_up, jes_down);
		//std::cout << cyan <<" JES ERROR FOR BIN :" << bin << "[" << zeejet->GetBinLowEdge(bin) << "] cen["  
		//	<< jes_cen << "]\t up[" << jes_up << "]\t down[" << jes_down << "] max = " << err << clearatt << std::endl;
		JESerr.push_back(err);

		//now to keep track of indivual background contributions
		vJES_phomc.push_back(GetMax(mcphojet->GetBinContent(bin), mcphojetJESUP->GetBinContent(bin), mcphojetJESDOWN->GetBinContent(bin)));
		vJES_diphomc.push_back(GetMax(diphojet->GetBinContent(bin), diphojetJESUP->GetBinContent(bin), diphojetJESDOWN->GetBinContent(bin)));
		vJES_ewkmc.push_back(GetMax(jes_cen_ewk, jes_up_ewk, jes_down_ewk));


	}

	// NOW GET THE STAT ERROR
	std::vector<float> statErr;
	for (int i=1; i <= hist_err->GetNbinsX(); i++)
	{
		statErr.push_back(hist_err->GetBinError(i));
	}


	//collecting PDF/ISR/FSR/Q2/AlphaS errors ---------------------------------
	std::vector<float> vPDFerror, vISRFSRerror, vQ2error, vAlphaSerror;
	//set all values to zero so no need to worry if any of these errors
	//are not evaluated for a certain plot
	for (int bin = 1; bin<= hist_err->GetNbinsX(); ++bin) 
	{
		vPDFerror.push_back(0);
		vISRFSRerror.push_back(0);
		vQ2error.push_back(0);
		vAlphaSerror.push_back(0);

	}

	//if (brname.length()>0 && QCDerrMtd == 2)		//_______________________ plots using PHOTON MC
	// these error should be included in all cases except 
	// in qcderr mtd 2 : 100% QCD in default method (30% qcd, 70%pho mc) for NJET plot
	//if (brname.length()>0 && QCDerrMtd == 2)	//_______________________ plots using PHOTON MC
	if (brname.length()>0)	//_______________________ plots using PHOTON MC
	{
		//std::cout << red << "IN HERE" << clearatt << std::endl;
		//now get PDF error, only when we use PHO MC
		TFile pdfSyst("PDFSyst_FitFunctions.root");
		assert (! pdfSyst.IsZombie() && "pdfSyst file not fouund!");

		TF1* tf1PDF = dynamic_cast<TF1*> (pdfSyst.Get(brname.c_str()));
		if (tf1PDF != NULL)
		{
			for (int bin = 1; bin<= hist_err->GetNbinsX(); ++bin)
			{
				float fPe = tf1PDF->Eval(hist_err->GetBinCenter(bin));
				if (fPe<0.01) fPe = 0.01;  //just like ISR/FSR and Q2 keep it 1% minimal
				vPDFerror.at(bin-1) = fPe * mcphojet->GetBinContent(bin);
			}
			delete tf1PDF;

		} else 
		{
			std::cout << __LINE__ << ":: WARNING !  NO PDF Systematic function found for " << brname  << std::endl;
		}


		// ISR/FSR syst
		TFile isrfsrSyst("ISRFSRsyst_FitFunctions.root");
		assert (! isrfsrSyst.IsZombie() && "isrfsrSyst file not fouund!");

		TF1* tf1ISR = dynamic_cast<TF1*> (isrfsrSyst.Get(brname.c_str()));
		if (tf1ISR != NULL)
		{
			for (int bin = 1; bin<= hist_err->GetNbinsX(); ++bin)
			{
				float fRe = tf1ISR->Eval(hist_err->GetBinCenter(bin));
				if (fRe<0.01) fRe = 0.01;	//keep minimum error for this at 1%. fits are not perfect elog#1482
				vISRFSRerror.at(bin-1) =  fRe * mcphojet->GetBinContent(bin);
			}
			delete tf1ISR;
		} else 
		{
			std::cout << __LINE__ << ":: WARNING !  NO ISR/FSR Systematic function found for "  << brname  << std::endl;
		}

		// Q2 syst
		TFile q2Syst("Q2Syst_FitFunctions.root");
		assert (! q2Syst.IsZombie() && "q2Syst file not fouund!");

		TF1* tf1Q2 = dynamic_cast<TF1*> (q2Syst.Get(brname.c_str()));
		if (tf1Q2 != NULL)
		{
			for (int bin = 1; bin<= hist_err->GetNbinsX(); ++bin)
			{
				float fQe = tf1Q2->Eval(hist_err->GetBinCenter(bin)); 
				if (fQe <0.01) fQe = 0.01;	//keep minimum error for this at 1%. fits are not perfect elog#1481,1480
				vQ2error.at(bin-1) = fQe * mcphojet->GetBinContent(bin);
			}
			delete tf1Q2;
		} else 
		{
			std::cout << __LINE__ << ":: WARNING !  NO Q2 Systematic function found for "  << brname  << std::endl;
		}


		// alpha_s syst
		GetAlphaS_Systematic(brname,vAlphaSerror,mcphojet);


	} //end of collecting PDF/ISR/FSR/Q2/ALPHA-S errors


	//new systematic for METHOD A, from METHOD B predictions
	std::vector<float> vMtdBerr;
	for (int bin = 1; bin<= hist_err->GetNbinsX(); ++bin) vMtdBerr.push_back(0);

	if (iHAVE_MTDB_ERRORS == 1)
	{
		if ((iSideband == iNOMINAL_SIDEBAND && sMtdBerr.length()>0) || (iSideband == iWEIGHTED_SIDEBAND && sMtdBerr.length()>0 && iMTDB_RERUN == 1))
		{

			std::stringstream mtdBfilename;
			if (iSideband == iNOMINAL_SIDEBAND)	mtdBfilename << "MethodB_ErrHists_" << sSAMPLE_NAME << ".root"; 
			else if (iSideband == iWEIGHTED_SIDEBAND && iMTDB_RERUN == 1) mtdBfilename << "MethodA_ErrHists_" << sSAMPLE_NAME << ".root";

			TFile f(mtdBfilename.str().c_str());
			if (f.IsZombie())
			{
				std::cout << red << __FUNCTION__ << ":" << __LINE__ << "method B error file not found!!!!" << clearatt << std::endl;
				return;
			}
			//data hist is same in both METHOD A and B.
			//so I'll use what I have in METHOD A. compare it Method B prediction

			std::stringstream bgh_name;
			bgh_name << "bg_" << sMtdBerr;

			TH1* bghist_mtdB = dynamic_cast<TH1*> (f.Get(bgh_name.str().c_str()));
			assert(bghist_mtdB != NULL && "method B err hist is null");
			//subtract the two hists and take the abosolute val as err
			//bghist_mtdB = MakeVariableBins(bghist_mtdB, false);
			cout << "mtd b bins = " << bghist_mtdB->GetNbinsX() << endl;
			cout << "hist_err bins = " << hist_err->GetNbinsX() << endl;

			//use that hist test method I found to check if these two have the sam bin boundries
			assert (bghist_mtdB->GetNbinsX() == hist_err->GetNbinsX() && "Metd B err hist has wrong number of bins");

			bghist_mtdB->Add(hist_err, -1);

			vMtdBerr.clear();
			for (int bin = 1; bin<= bghist_mtdB->GetNbinsX(); ++bin) 
			{
				vMtdBerr.push_back(fabs(bghist_mtdB->GetBinContent(bin)));


				//these are manual tweaks to total systematic to smooth out any abnormally
				//large errors in some bins.

				//for iG30 plots
				if (iSample == iG30  && iSideband == iNOMINAL_SIDEBAND)
				{
					if (jets == 2 && which == "InvMass_pj1")
					{
						if (bin == 9) vMtdBerr.at(bin-1) *= 0.5;

					}
				}




				//hack for nice plots, smooth out some abnormally large values.
				/*if (jets== 1 && which == "InvMass_pj1")
				{
					std::cout << "bin, loedge [val] = " << bin << ", " << bghist_mtdB->GetBinLowEdge(bin)
						<< "["<< bghist_mtdB->GetBinError(bin) << std::endl;
					if (bin == 15) vMtdBerr.at(bin-1) *= 0.5;
					if (bin == 17) vMtdBerr.at(bin-1) *= 0.8;
				}

				if (jets == 2 && which == "InvMass_pj1")
				{
					std::cout << "bin, loedge [val] = " << bin << ", " << bghist_mtdB->GetBinLowEdge(bin)
						<< "["<< bghist_mtdB->GetBinError(bin) << std::endl;
					//if (bin == 15) vMtdBerr.at(bin-1) *= 0.5;
					if (bin == 18) vMtdBerr.at(bin-1) *= 0.5;
				}

				if (jets== 1 && which == "Met")
				{
					std::cout << "bin, loedge [val] = " << bin << ", " << bghist_mtdB->GetBinLowEdge(bin)
						<< "["<< bghist_mtdB->GetBinError(bin) << std::endl;
					if (bin == 15) vMtdBerr.at(bin-1) *= 0.5;
					if (bin == 17) vMtdBerr.at(bin-1) *= 0.8;
				}
				*/


			}

		}
	} //if (iHAVE_MTDB_ERRORS)
	

	// NOW COLLECT ALL ERRORS AND PUT THEM TOGETHER FOR ONE FINAL NUMBER
	std::vector<float> ERRORS;

	std::cout << __FUNCTION__ << ":" << __LINE__ << "::All Errors" << std::endl;
	std::cout << red <<setw(4) << "bin"<< setw(8) << "[LoEdge]"<< setw(8) << "JES"
		<< setw(8) << "PhoFrac" << setw(8) << "EWK Lum" 
		<< setw(8) << "Cosmic" << setw(6) << "Halo"
		<< setw(8) << "Stat" << setw(10) << "QCD Id" 
		<< setw(10) << "PDF" 
		<< setw(10) << "Q^2" 
		<< setw(10) << "IsrFsr" 
		<< setw(10) << "AlphaS" 
		<< setw(10) << "REWGT" 
		<< setw(10) << "MTB B" 
		<< setw(10) << "Tot Err" 
		<< setw(10) << "Sum BG" 
		<< setw(10) << "Data" 
		<< setw(10) << "D-B/B" 
		<< clearatt << std::endl; 

	assert ( (cosmicErr.size() == haloErr.size() 
				&& haloErr.size()==qcdmcMixErr.size() 
				&& qcdmcMixErr.size()==JESerr.size()
				&& JESerr.size() == statErr.size() 
				&& statErr.size() == cosmicErr.size() 
				&& cosmicErr.size() == vPDFerror.size()
				&& vPDFerror.size() == vISRFSRerror.size()
				&& vISRFSRerror.size() == vQ2error.size()
				&& vQ2error.size() == vAlphaSerror.size()
				&& vAlphaSerror.size() == IdErr.size()
				&& IdErr.size() == vDiPhoLumErr.size()
				&& vDiPhoLumErr.size() == vZeeLumErr.size()
				&& vZmmLumErr.size() == vZeeLumErr.size()
				&& vZttLumErr.size() == vZeeLumErr.size()
				&& vWenLumErr.size() == vWtnLumErr.size()
				&& vWenLumErr.size() == vWmnLumErr.size()
				&& vWmnLumErr.size() == vSidebandRwtErr.size()
				&& vSidebandRwtErr.size() == vMtdBerr.size()
				) 
			&& "MakeHistLogAndRatioV2::Error vectors are not same in size!");


	std::vector<float> vTotalEWKLumErr;

	for (unsigned int i=0; i < vZeeLumErr.size(); i++)
	{
		vTotalEWKLumErr.push_back(sqrt (pow(vZeeLumErr[i],2) 
					+ pow(vZmmLumErr[i],2)+ pow(vZttLumErr[i],2) + pow(vWenLumErr[i],2) 
					+ pow(vWmnLumErr[i],2) + pow(vWtnLumErr[i],2)));
	}


	for (unsigned int i=0; i < cosmicErr.size(); i++)
	{
		float sum = pow(cosmicErr[i],2) 
			+ pow(haloErr[i],2)
			+ pow(vDiPhoLumErr[i],2)
			+ pow(vTotalEWKLumErr[i],2) 
			+ pow(qcdmcMixErr[i],2)
			+ pow(IdErr.at(i),2) 
			+ pow(JESerr[i],2) 
			+ pow(statErr[i],2)
			+ pow(vPDFerror.at(i),2)
			+ pow(vISRFSRerror.at(i),2)
			+ pow(vQ2error.at(i),2)
			+ pow(vAlphaSerror.at(i),2)
			+ pow(vSidebandRwtErr.at(i),2)
			+ pow(vMtdBerr.at(i),2)
			;

		ERRORS.push_back(sqrt(sum));

		// now set the sum of background hists error to this total
		if (hist_err->GetBinContent(i+1)) //no errors if there are no background prediction!?
		{
			hist_err->SetBinError(i+1, ERRORS.at(i));
		}
		
		std::cout << blue  << setiosflags(ios::fixed) << setprecision(1) 
			<< setw(4) << i+1 << "[" << setw(6) << hist_err->GetBinLowEdge(i+1) << "]" << setw(8) << JESerr.at(i)
			<< setw(8) << qcdmcMixErr.at(i) << setw(8) << vTotalEWKLumErr.at(i) 
			<< setw(8) << cosmicErr.at(i) << setw(6) << haloErr.at(i)
			<< setw(8) << statErr.at(i) << setw(10) << IdErr.at(i) 
			<< setw(10) << vPDFerror.at(i) 
			<< setw(10) << vQ2error.at(i) 
			<< setw(10) << vISRFSRerror.at(i) 
			<< setw(10) << vAlphaSerror.at(i)
			<< setw(10) << vSidebandRwtErr.at(i)
			<< setw(10) << vMtdBerr.at(i)
			<< setw(10) << ERRORS.at(i) 
			<< setw(10) << hist_err->GetBinContent(i+1) 
			<< setw(10) << phojet->GetBinContent(i+1) 
			<< setw(10) << setiosflags(ios::fixed) << setprecision(2);
		if (hist_err->GetBinContent(i+1)>0) 
		{
			std::cout << (phojet->GetBinContent(i+1) - hist_err->GetBinContent(i+1)) /hist_err->GetBinContent(i+1);
		} else
		{
			std::cout << 0.0;
		}
		std::cout << clearatt << std::endl; 

	}

	//create a copy for writing to a file in METHOD B errors
	TH1F* hist_err_copy = dynamic_cast<TH1F*> (hist_err->Clone(sMtdBerr.c_str()));
	hist_err_copy->SetDirectory(0);
//	hist_err_copy->Print("all");

	//ONLY THE COSMETIC CHANGES TO THE HISTOGRAMS SHOULD BE DONE BELOW THIS POINT!
	

	//if (iDRAWERROR_BREAKDOWN && ! iSideband == iWEIGHTED_SIDEBAND)
	if (iDRAWERROR_BREAKDOWN)
	{
		DrawErrorBreakDown(QCDerrMtd, iSideband, hist_err,JESerr,qcdmcMixErr, vTotalEWKLumErr, cosmicErr,haloErr, statErr,IdErr,
				vPDFerror, vQ2error,vISRFSRerror,vAlphaSerror, vSidebandRwtErr, vMtdBerr, ERRORS, xtitle, jets, which);
	}

	TCanvas *c1 = new TCanvas("HIST","HIST",600,800);
	TPad *pLOG = new TPad("Plog","Log Plot",0.01,0.4,0.99,0.99);   //in NDC cdts
	TPad *pRAT = new TPad("Pratio","Ratio Plot",0.01,0.01,0.99,0.4);  //in NDC cdts
	pLOG->SetBorderMode(0);
	pLOG->SetBorderSize(0);
	pLOG->SetBottomMargin(0);
	pRAT->SetBorderMode(0);
	pRAT->SetBorderSize(0);
	pRAT->SetTopMargin(0);
	pRAT->SetBottomMargin(0.1);
	c1->cd();
	pLOG->Draw();
	pLOG->cd();



	/**************************************************************
	 *               SET HIST COLORS
	 *************************************************************/

	Int_t linecolor=47;

	Int_t cosmic_color = 8;
	Int_t halo_color   = 8;

	mcphojet->SetLineColor(6);
	mcphojet->SetFillColor(6);
	mcphojet->SetMarkerColor(6);
	//mcphojet->SetFillColor(29);
	qcdjet->SetLineColor(kYellow);
	qcdjet->SetFillColor(kYellow);
	qcdjet->SetMarkerColor(kYellow);
	cosmicjet->SetLineColor(cosmic_color);
	cosmicjet->SetFillColor(cosmic_color);
	cosmicjet->SetMarkerColor(cosmic_color);
	halojet->SetLineColor(halo_color);
	halojet->SetFillColor(halo_color);
	halojet->SetMarkerColor(halo_color);
	//int ewkColor = 29;
	int ewkColor = kOrange-3;
	zeejet->SetLineColor(ewkColor);
	zeejet->SetFillColor(ewkColor);
	zeejet->SetMarkerColor(ewkColor);
	zmmjet->SetLineColor(ewkColor);
	zmmjet->SetFillColor(ewkColor);
	zttjet->SetLineColor(ewkColor);
	zttjet->SetFillColor(ewkColor);
	wenjet->SetLineColor(ewkColor);
	wenjet->SetFillColor(ewkColor);
	wmnjet->SetLineColor(ewkColor);
	wmnjet->SetFillColor(ewkColor);
	wtnjet->SetLineColor(ewkColor);
	wtnjet->SetFillColor(ewkColor);

	//diphojet->SetFillColor(kGreen+4);
	//diphojet->SetLineColor(kGreen+4);
	diphojet->SetFillColor(kCyan);
	diphojet->SetLineColor(kCyan);
	diphojet->SetMarkerColor(kCyan);


	phojet->SetLineColor(kBlack);
	phojet->SetMarkerStyle (8);

	//this is the line histo showing 100% of QCD 0% MC pho
	//hist_err->SetFillStyle(3001);
	hist_err->SetFillStyle(3344);
	//hist_err->SetFillColor(13);
	hist_err->SetFillColor(kBlue);
	hist_err->SetLineColor(kWhite);
	hist_err->GetXaxis()->SetTitleColor(kBlack);
	hist_err->GetYaxis()->SetTitleColor(kBlack);


	/*************************************************************
	 *   MAKE THE LOG PLOT
	 ************************************************************/
	THStack *hs = new THStack ("hs", NULL);

	hs->Add(halojet);
	hs->Add(cosmicjet);
	hs->Add(zmmjet);
	hs->Add(wtnjet);
	hs->Add(zttjet);
	hs->Add(wmnjet);
	hs->Add(wenjet);
	hs->Add(zeejet);
	hs->Add(diphojet);
	hs->Add(qcdjet);

	//if (QCDerrMtd == 2) {//plots that use 30%70%
	if (iSideband == iNOMINAL_SIDEBAND) {//plots that use 30%70%
		hs->Add(mcphojet);
		std::cout << cyan << __FUNCTION__ << ":" << __LINE__ << ": MIX METHOD" << std::endl;
	}



	/*************************************************************
	 *  Set up legend/s
	 *************************************************************/

	const float fLegend_x1=0.56;
	const float fLegend_y1=0.6;
	const float fLegend_x2=0.9;
	const float fLegend_y2=0.9;
	TLegend *leg = new TLegend (fLegend_x1,fLegend_y1,fLegend_x2,fLegend_y2,"","NDC");
	leg->SetTextFont(iAXIS_LABEL_FONT);
	leg->SetTextSize(0.03);
	std::string str_pho,str_cosmic,str_halo,str_zee,str_zmm,str_ztt,str_wen,str_wmn,str_wtn;
	std::string str_diphomc;

	if (jets==1) {
		//str_pho 	= "Data (#gamma + #geq 1 Jet + #slash{E}_{T}>25GeV";	
		str_pho 	= "Data";	
		str_cosmic = "#gamma^{cosmic} + #geq 1 Jet";
		str_zee 	= "Z->ee MC (e + #geq 1 Jet)";
		str_wen 	= "W->e#nu MC (e + #geq 1 Jet)";
		str_ztt 	= "Z->#tau#tau MC (e + #geq 1 Jet)";
		str_wtn 	= "W->#tau#nu MC (e + #geq 1 Jet)";
		str_wmn 	= "W->#mu#nu MC (e + #geq 1 Jet)";
		str_zmm 	= "Z->#mu#mu MC (e + #geq 1 Jet)";
		str_halo 	= "#gamma^{halo} + #geq 1 Jet";
		str_diphomc= "Di-#gamma MC + #geq 1 Jet";
	}
	if (jets==2) {
		str_pho 	= "Data";	
		str_cosmic = "#gamma^{cosmic} + #geq 2 Jets";
		str_zee 	= "Z->ee MC (e + #geq 2 Jets)";
		str_wen 	= "W->e#nu MC (e + #geq 2 Jets)";
		str_ztt 	= "Z->#tau#tau MC (e + #geq 2 Jets)";
		str_wtn 	= "W->#tau#nu MC (e + #geq 2 Jets)";
		str_wmn 	= "W->#mu#nu MC (e + #geq 2 Jets)";
		str_zmm 	= "Z->#mu#mu MC (e + #geq 2 Jets)";
		str_halo 	= "#gamma^{halo} + #geq 2 Jets";
		str_diphomc= "Di-#gamma MC + #geq 2 Jets";
	}


	if (QCDerrMtd == 1) {//plots that use 100% qcd
		std::cout << "USING METHOD 1 " << std::endl;
		leg->AddEntry(phojet,str_pho.c_str());
		if (iSideband == iNOMINAL_SIDEBAND)
		{
			leg->AddEntry(qcdjet, "QCD");
		} else {
			leg->AddEntry(qcdjet, "QCD (weighted)");
		}
		leg->AddEntry(diphojet,"Di-#gamma");
		leg->AddEntry(zeejet,"EWK");
		//leg->AddEntry(diphojet,str_diphomc.c_str());
		//leg->AddEntry(cosmicjet,str_cosmic.c_str());
		//leg->AddEntry(halojet,str_halo.c_str());
		leg->AddEntry(cosmicjet,"Non-collision (Data)");
		leg->AddEntry(hist_err,"Systematic Uncertainty");

	} else if (QCDerrMtd == 2) {//plots that use 30%70%
		//leg->AddEntry(phojet,str_pho.c_str());
		//leg->AddEntry(mcphojet,mc_str.c_str());
		//leg->AddEntry(qcdjet,qcd_str.c_str());
		leg->AddEntry(mcphojet, "#gamma MC");
		leg->AddEntry(qcdjet, "QCD");
		leg->AddEntry(diphojet,"Di-#gamma");
		leg->AddEntry(zeejet,"EWK");
		//leg->AddEntry(diphojet,str_diphomc.c_str());
		//leg->AddEntry(cosmicjet,str_cosmic.c_str());
		//leg->AddEntry(halojet,str_halo.c_str());
		leg->AddEntry(cosmicjet,"Non-collision");
		//leg->AddEntry(hist_err,"Systematic Uncertainty");
	} else {
		std::cerr << "bad value: QCDerrMtd=" << QCDerrMtd << std::endl;
		return;
	};

	leg->SetBorderSize (1);
	leg->SetFillColor (10);

	//dynamically adjust the y-scale
	double dYscale_max = hist_err->GetBinContent(hist_err->GetMaximumBin()) * 10.;
	double dYscale_min = 1e10;
	for (int ibin=1; ibin <= hist_err->GetNbinsX();++ibin)
	{
		if (hist_err->GetBinContent(ibin) <= 0) continue; //some bins can be slightly <0 after subtracting other backgrounds.
		if (hist_err->GetBinContent(ibin) < dYscale_min)
		{
			dYscale_min = hist_err->GetBinContent(ibin);
		}
	}
	dYscale_min *= 1./2.;
	std::cout << blue << __LINE__ << ": Y-scale min, max = " << dYscale_min << ", " << dYscale_max << std::endl;


	//new TCanvas();
	gStyle->SetOptStat(0);
	//gStyle->SetTextFont(132);
	gStyle->SetTextFont(iAXIS_LABEL_FONT);
	gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();

	if (bDO_VAR_BINNING) 
	{
		hs->SetMinimum(dYscale_min);
		//hs->SetMinimum(0.005);
	} else hs->SetMinimum(0.5);
	hs->SetMaximum(dYscale_max);


	//original method. put this back after debug --11-22-2009
	hs->Draw("HIST");		//need this as I am calling sumw2 for all hists. if not it will draw all hists with error bars
	//new TCanvas();
	//hs->Draw("HIST nostack,e1p");		//need this as I am calling sumw2 for all hists. if not it will draw all hists with error bars

	hs->GetXaxis()->SetTitle (xtitle.c_str());
	std::ostringstream ytitle;
	if (which == "NJet") ytitle << "Events";
	else 
	{
		/*ytitle << "Events / " << std::setprecision(2) << phojet->GetBinWidth(1);
		if (which == "Et_pho" || which == "Et_j1" || which == "Ht" || which == "Met")  ytitle << " GeV/c";
		else if (which.find("Mass") != string::npos)  ytitle << " GeV/c^{2}";
		*/
		if (which == "Et_pho" || which == "Et_j1")  ytitle << "dN / dE_{T}";
		else if (which == "Ht") ytitle << "dN / dH_{T}";
		else if (which == "Met")  ytitle << "dN / d#slash{E}_{T}";
		else if (which.find("Mass") != string::npos)  ytitle << "dN / dM";

	}
	hs->GetYaxis()->SetTitle (ytitle.str().c_str());

	int labelfont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	int titlefont = 10 * 4 + 2;		//10 * fond ID + precision (2 = scalable)
	hs->GetXaxis()->SetLabelFont(iAXIS_LABEL_FONT);
	hs->GetYaxis()->SetLabelFont(iAXIS_LABEL_FONT);
	hs->GetYaxis()->SetLabelSize(0.04);
	//hs->GetXaxis()->SetLabelSize(0.05);
	hs->GetXaxis()->SetTitleFont(iAXIS_LABEL_FONT);
	hs->GetYaxis()->SetTitleFont(iAXIS_LABEL_FONT);
	hs->GetYaxis()->SetTitleSize(0.04);
	//hs->GetXaxis()->SetTitleSize(0.05);
	hs->GetXaxis()->SetTitleOffset(0.9);
	//hs->GetYaxis()->SetTitleOffset(0.965);
	hs->GetYaxis()->SetTitleOffset(0.99);
	hs->GetXaxis()->CenterTitle(true);
	hs->GetYaxis()->CenterTitle(true);
	hs->GetYaxis()->SetTitleColor(kBlack);
	hs->GetXaxis()->SetTitleColor(kBlack);


	hist_err->SetDrawOption("HIST");
	hist_err->Draw("sameE2");
	phojet->Draw("same");
	leg->Draw();


	TPaveText *tp = new TPaveText(0.4,0.91,0.9,0.96,"NDC");
	tp->SetLineColor(10);
	tp->SetBorderSize(0);
	tp->SetTextFont(iTITLE_FONT);
	tp->SetTextSize(iTITLE_FONT_SIZE);
	tp->SetTextAlign(10*3);
	tp->AddText("CDF Run II Preliminary 4.8 fb^{-1}");
	tp->Draw();

	/**********************************************************
	 * Set Titles with less info for publications
	 **********************************************************/
	std::string sLegendTitle("");
	if (iPRL_MODE)
	{
		if (iSample == iG30 && jets == 1) sLegendTitle = "#gamma + #geq 1 Jet";
		else if (iSample == iG30 && jets == 2) sLegendTitle = "#gamma + #geq 2 Jet";
		//else if (iSample == iG30 && jets == 2) sLegendTitle = "#gamma + == 2 Jet";
		else if (iSample == iG40) sLegendTitle = "#gamma^{E_{T}>40GeV}_{|#eta|<1.1}+#geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2}";
		else if (iSample == iG30MET20 && jets == 1) sLegendTitle = "#gamma + #geq 1 Jet + #slash{E}_{T} > 20GeV";
		else if (iSample == iG30MET20 && jets == 2) sLegendTitle = "#gamma + #geq 2 Jet + #slash{E}_{T} > 20GeV";
		else if (iSample == iG30MET30 && jets == 1) sLegendTitle = "#gamma + #geq 1 Jet + #slash{E}_{T} > 30GeV";
		else if (iSample == iG30MET30 && jets == 2) sLegendTitle = "#gamma + #geq 2 Jet + #slash{E}_{T} > 30GeV";
		else if (iSample == iG30MET25) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1} + #geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2} + #slash{E}_{T}>25GeV";
		else if (iSample == iG30MET40) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1} + #geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2} + #slash{E}_{T}>40GeV";
		else if (iSample == iG30EXCL1J) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+ == 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2}";
		else if (iSample == iG30MET50 && jets == 1) sLegendTitle = "#gamma + #geq 1 Jet + #slash{E}_{T} > 50 GeV";
		else if (iSample == iG30MET50 && jets == 2) sLegendTitle = "#gamma + #geq 2 Jet + #slash{E}_{T} > 50 GeV";

	} else {	
		if (iSample == iG30 && jets == 1) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+#geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.0}, #Delta#phi(#slash{E}_{T},jet)>0.4";
		else if (iSample == iG30 && jets == 2) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+#geq 2 Jet^{E_{T}>15GeV}_{|#eta|<3.0}, #Delta#phi(#slash{E}_{T},jet)>0.4";
		else if (iSample == iG40) sLegendTitle = "#gamma^{E_{T}>40GeV}_{|#eta|<1.1}+#geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2}";
		else if (iSample == iG30MET20 && jets == 1) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+#geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2}+#slash{E}_{T}>20GeV, #Delta#phi(#slash{E}_{T},jet)>0.4";
		else if (iSample == iG30MET20 && jets == 2) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+#geq 2 Jet^{E_{T}>15GeV}_{|#eta|<3.2}+#slash{E}_{T}>20GeV, #Delta#phi(#slash{E}_{T},jet)>0.4";
		else if (iSample == iG30MET25) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1} + #geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2} + #slash{E}_{T}>25GeV";
		else if (iSample == iG30MET30 && jets == 1) sLegendTitle = "#gamma + #geq 1 Jet + #slash{E}_{T} > 30GeV";
		else if (iSample == iG30MET30 && jets == 2) sLegendTitle = "#gamma + #geq 2 Jet + #slash{E}_{T} > 30GeV";
		else if (iSample == iG30MET40) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1} + #geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2} + #slash{E}_{T}>40GeV";
		else if (iSample == iG30EXCL1J) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+ == 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2}";
		else if (iSample == iG30MET50 && jets == 1) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+#geq 1 Jet^{E_{T}>15GeV}_{|#eta|<3.2}+#slash{E}_{T}>50GeV, #Delta#phi(#slash{E}_{T},jet)>0.4";
		else if (iSample == iG30MET50 && jets == 2) sLegendTitle = "#gamma^{E_{T}>30GeV}_{|#eta|<1.1}+#geq 2 Jet^{E_{T}>15GeV}_{|#eta|<3.2}+#slash{E}_{T}>50GeV, #Delta#phi(#slash{E}_{T},jet)>0.4";
	}

	
	
	TPaveText *tpTitle = new TPaveText(0.1,0.91,0.4,0.96,"NDC");
	tpTitle->SetLineColor(iWHITE);
	tpTitle->SetBorderSize(0);
	tpTitle->SetTextFont(iTITLE_FONT);
	tpTitle->SetTextSize(iTITLE_FONT_SIZE);
	tpTitle->SetTextAlign(10*1);
	tpTitle->AddText(sLegendTitle.c_str());
	tpTitle->Draw();


	const float fTp2_h = 0.045;
	TPaveText *tp2 = new TPaveText(fLegend_x1,fLegend_y1-fTp2_h,fLegend_x2,fLegend_y1,"NDC");
	tp2->SetLineColor(kBlack);
	tp2->SetBorderSize(1);
	//tp->SetTextFont();
	//if (iSideband == iNOMINAL_SIDEBAND) tp2->AddText("Nominal sideband");
	if (iSideband == iNOMINAL_SIDEBAND)
	{
		tp2->AddText("Method A");
		tp2->SetFillColor(6);
	} else if (iSideband == iWEIGHTED_SIDEBAND) 
	{
		//tp2->AddText("Iso added sideband");
		//tp2->AddText("reweighed to total background");
		//tp2->AddText("Weighted Sideband");
		//tp2->AddText("to background");
		tp2->AddText("Method B");
		tp2->SetFillColor(kYellow);
	}

	tp2->Draw();
	//if (! iPRL_MODE) tp2->Draw(); // make this addition info only when I am testing things to make sure 
											// I do the correct thing


	/******************************************************************
	 *     Marks the overflow bin (if any)
	 ******************************************************************/
	bool mark_overflow = false;
	if (QCDerrMtd == 1) {//plots that use 100% qcd
		if ( halojet->GetBinContent(halojet->GetNbinsX()) ||
				cosmicjet->GetBinContent(cosmicjet->GetNbinsX()) ||
				zmmjet->GetBinContent(zmmjet->GetNbinsX()) ||
				wtnjet->GetBinContent(wtnjet->GetNbinsX()) ||
				zttjet->GetBinContent(zttjet->GetNbinsX()) ||
				wmnjet->GetBinContent(wmnjet->GetNbinsX()) ||
				wenjet->GetBinContent(wenjet->GetNbinsX()) ||
				zeejet->GetBinContent(zeejet->GetNbinsX()) ||
				qcdjet->GetBinContent(qcdjet->GetNbinsX()) ||
				diphojet->GetBinContent(diphojet->GetNbinsX()) )
			mark_overflow = true;

	} else if (QCDerrMtd == 2) {//plots that use 30%70%
		if ( halojet->GetBinContent(halojet->GetNbinsX()) ||
				cosmicjet->GetBinContent(cosmicjet->GetNbinsX()) ||
				zmmjet->GetBinContent(zmmjet->GetNbinsX()) ||
				wtnjet->GetBinContent(wtnjet->GetNbinsX()) ||
				zttjet->GetBinContent(zttjet->GetNbinsX()) ||
				wmnjet->GetBinContent(wmnjet->GetNbinsX()) ||
				wenjet->GetBinContent(wenjet->GetNbinsX()) ||
				zeejet->GetBinContent(zeejet->GetNbinsX()) ||
				qcdjet->GetBinContent(qcdjet->GetNbinsX()) ||
				diphojet->GetBinContent(diphojet->GetNbinsX()) ||
				mcphojet->GetBinContent(mcphojet->GetNbinsX()) )
			mark_overflow = true;
	}

	if (mark_overflow)
	{
		gPad->Range(0,0,1,1);

		TLine *line = new TLine();
		line->SetLineWidth(2);
		line->DrawLineNDC(0.895,0.218,0.93,0.218);

		TText *tt = new TText(0.95,0.15,"Overflow bin");
		tt->SetTextFont(titlefont);
		tt->SetTextAngle(90);
		tt->SetTextSize(0.035);
		//tt->SetTextColor(kBlue);
		tt->SetNDC();

		tt->Draw();
	}

	//KS TEST RESULTS
	const float fKSProb = hist_err->KolmogorovTest(phojet);
	std::stringstream ssKStest;
	ssKStest << "KS Prob. " << std::setprecision(2) << fKSProb;
	//new TCanvas();
	//phojet->Draw();
	//hist_err->Draw("same");
	std::cout << red <<  " *****************************" << std::endl;
	std::cout << "KS Prob = " << fKSProb << std::endl;
	std::cout << " *****************************"  << clearatt << std::endl;
	TPaveText *tpKStest = new TPaveText(fLegend_x1,fLegend_y1-2.*fTp2_h, fLegend_x2,fLegend_y1-fTp2_h,"NDC");

	tpKStest->SetLineColor(kBlack);
	tpKStest->SetBorderSize(1);
	tpKStest->AddText(ssKStest.str().c_str());
	//tpKStest->Draw();


	//std::stringstream myname;
	//myname << "logplot" << jets << "_" << which << ".eps";
	//gPad->Print(myname.str().c_str());

	c1->cd();
	pRAT->Draw();
	//hist_err should have all the errors set by now
	MakeRatioPlot(phojet, hist_err, xtitle, pRAT, false);
	c1->cd();


	//DrawErrorBreakDown(QCDerrMtd, hist_err,JESerr,qcdmcMixErr,cosmicErr,haloErr, statErr,IdErr,
	//		vPDFerror, vQ2error,vISRFSRerror,vAlphaSerror,ERRORS, xtitle, jets, which);


	std::cout << " >>>>>>> Settings used for " << jets << " jet/s " << which << " plot <<<<<<" << std::endl;
	std::cout << "DATA Lum               = " << fLUM_DATA << std::endl;	
	std::cout << "iSAMPLE                = " << sSAMPLE_NAME << std::endl;
	std::cout << "Fake pho frac+/-sigma  = " << fFAKE_PHO_FRAC_DEF << "+/-" << fFAKE_PHO_FRAC_SIGMA << std::endl;
	std::cout << "METHOD                 = " << sMethod << std::endl;
	std::cout << "Apply W+jets MET corr to W MC? = " <<  bCorrectWMC_Met << " (0=NO)" << std::endl;
	std::cout << "Sideband               = " << sSideband << std::endl; 
	std::cout << "DATA integral          = " << phojet->Integral("width") << std::endl;
	std::cout << "Halo integral          = " << halojet->Integral("width") << std::endl;
	std::cout << "Cosmic integral        = " << cosmicjet->Integral("width") << std::endl;
	std::cout << "Qcdjet Integral        = " << qcdjet->Integral("width") << std::endl;
	std::cout << "MC PHO Integral        = " << mcphojet->Integral("width") << std::endl;
	std::cout << "EWK Total Integral     = " << zeejet->Integral("width") + zmmjet->Integral("width")
		+ zttjet->Integral("width") + wenjet->Integral("width") 
		+ wmnjet->Integral("width") + wtnjet->Integral("width") 
		<< std::endl;
	std::cout << "Di-pho Integral        = " << diphojet->Integral("width") << clearatt << std::endl;
	std::cout << red << "DiPHO MC Lum err is not included in the error breakdown list!"  << clearatt << std::endl;



	//write out the final data/err and total syst to a file for scanning
	if (iSideband == iWEIGHTED_SIDEBAND && iMTDB_RERUN == 0) //write to a error file, METHOD B HISTS. do not write in 2nd run
	{
		std::stringstream mtdBfilename;
		mtdBfilename << "MethodB_ErrHists_" << sSAMPLE_NAME << ".root"; 
		TFile f(mtdBfilename.str().c_str(),"UPDATE");
		
		std::stringstream dh_name, bgh_name;
		dh_name << "data_" << sMtdBerr;
		bgh_name << "bg_" << sMtdBerr;
		phojet->Write(dh_name.str().c_str());
		hist_err_copy->Write(bgh_name.str().c_str());
	} else if (iSideband == iNOMINAL_SIDEBAND) //write to a error file //METHOD A HISTST
	{
		std::stringstream mtdBfilename;
		mtdBfilename << "MethodA_ErrHists_" << sSAMPLE_NAME << ".root"; 
		TFile f(mtdBfilename.str().c_str(),"UPDATE");
		
		std::stringstream dh_name, bgh_name;
		dh_name << "data_" << sMtdBerr;
		bgh_name << "bg_" << sMtdBerr;
		phojet->Write(dh_name.str().c_str());
		hist_err_copy->Write(bgh_name.str().c_str());
	}
	
	//for p-value scan method.
	/*TFile f("OutPut.root","RECREATE");
	phojet->Write();
	hist_err->Write();
	TH1F *syst_plus = dynamic_cast<TH1F*> (phojet->Clone("SystPlus"));
	for (int bin=1; bin<=syst_plus->GetNbinsX(); ++bin)
	{
		syst_plus->SetBinError(bin,0);
		syst_plus->SetBinContent(bin, ERRORS.at(bin-1));
	}
	TH1F *syst_minus = dynamic_cast<TH1F*> (syst_plus->Clone("SystMinus"));
	syst_plus->Write();
	syst_minus->Write();
	*/

	//generate number for tables
	
	double dTotData = 0, dTotBg = 0, dTotSys = 0;

	
	std::cout << "bin[loEdge] \t data \t bg \t syst" << std::endl;
	for (int bin = 1; bin <= phojet->GetNbinsX(); ++bin)
	//for (int bin = 4; bin <=4; ++bin)
	{
		//float fBinWidth = phojet->GetBinWidth(bin);
		std::cout << bin << "[" << phojet->GetBinLowEdge(bin) << "] \t " << phojet->GetBinContent(bin)   
				<< " \t "<< hist_err->GetBinContent(bin)  <<" \t " 
				<< ERRORS.at(bin-1) << std::endl;
		dTotData += phojet->GetBinContent(bin) ;
		dTotBg += hist_err->GetBinContent(bin) ;
		dTotSys += ERRORS.at(bin-1); 
	}
		
	std::cout << "Default const bin width = " << fOrigBinWidth << std::endl;
	std::cout << "DATA / BG+/-Sys = " << dTotData 
			<< " / " << dTotBg
			<< " / " << dTotSys
			<< std::endl;
	
	//photon MC errors
	//JES+EM/pho frac+isr+fsr+pdf+alphas+q2
	double dPhoMcSyst = 0, dQCDSyst =0, dDiphoSyst =0, dEWKSyst =0, dCosmicSyst=0, dHaloSyst=0;
	double dPhoMcStat = 0 , dQCDStat = 0 , dDiPhoStat =0, dEWKStat = 0, dCosmicStat = 0, dHaloStat=0;
	double tJES=0, tqcdMix =0, tPdf = 0, tIsr = 0, tQ2 = 0, tAlphas=0, tQcdId =0;
	double tLum=0, tStat=0;
	//for (int bin = 4; bin <=4; ++bin)
	for (int bin = 1; bin <= phojet->GetNbinsX(); ++bin)
	{
		//const double dBW = phojet->GetBinWidth(bin);
		const double dBW = 1.;
		
		tJES =  (vJES_phomc.at(bin-1), 2) + pow(IdErr.at(bin-1), 2)
					+ pow(vJES_ewkmc.at(bin-1), 2); 
		tqcdMix =  pow(qcdmcMixErr.at(bin-1), 2);


		dPhoMcSyst += (pow(vJES_phomc.at(bin-1) * dBW, 2) + pow(qcdmcMixErr.at(bin-1) * dBW, 2) 
						+ pow(vPDFerror.at(bin-1) * dBW, 2) + pow(vISRFSRerror.at(bin-1) * dBW, 2) 
						+ pow(vQ2error.at(bin-1) * dBW, 2) + pow(vAlphaSerror.at(bin-1) * dBW,2));
		//dQCDSyst += (pow(qcdmcMixErr.at(bin-1) * dBW, 2) + pow(IdErr.at(bin-1) * dBW, 2));
		dQCDSyst += (pow(IdErr.at(bin-1) * dBW, 2));  //pho frac err included only in pho mc , otherwise disagreement with the histogram values
		if (iSideband==iWEIGHTED_SIDEBAND) dQCDSyst += (pow(vSidebandRwtErr.at(bin-1) * dBW, 2)); //add MC error to this for method B.
		dDiphoSyst += (pow(vJES_diphomc.at(bin-1) * dBW, 2) + pow(vDiPhoLumErr.at(bin-1) * dBW,2));
		dEWKSyst += pow(vTotalEWKLumErr.at(bin-1) * dBW, 2) + pow(vJES_ewkmc.at(bin-1) * dBW, 2);
		dCosmicSyst += pow(cosmicErr.at(bin-1)  * dBW, 2);
		dHaloSyst += pow(haloErr.at(bin-1) * dBW, 2);
		
		dPhoMcStat += mcphojet->GetBinContent(bin-1)* dBW;
		dQCDStat += qcdjet->GetBinContent(bin-1)* dBW;
		dDiPhoStat += diphojet->GetBinContent(bin-1)* dBW;
		dEWKStat += (zeejet->GetBinContent(bin-1) + zmmjet->GetBinContent(bin-1)+ zttjet->GetBinContent(bin-1)
							+ wenjet->GetBinContent(bin-1) + wmnjet->GetBinContent(bin-1) + wtnjet->GetBinContent(bin-1))* dBW;
		dCosmicStat += cosmicjet->GetBinContent(bin-1) * dBW;
		dHaloStat += halojet->GetBinContent(bin-1) * dBW;
		
		tStat = dPhoMcStat+dQCDStat+dDiPhoStat+dEWKStat+dCosmicStat+dHaloStat;

	}

	std::cout << "My JES/ QCDMIX /STAT = " <<sqrt(tJES) << " / "<< sqrt(tqcdMix) << " / " << sqrt(tStat) << std::endl;
		dPhoMcSyst = sqrt(dPhoMcSyst);
		dQCDSyst = sqrt(dQCDSyst);
		dDiphoSyst = sqrt(dDiphoSyst);
		dEWKSyst = sqrt(dEWKSyst);
		dCosmicSyst = sqrt(dCosmicSyst);
		dHaloSyst = sqrt(dHaloSyst);

		/*dPhoMcStat = sqrt(dPhoMcStat);
		dQCDStat = sqrt(dQCDStat);
		dDiPhoStat = sqrt(dDiPhoStat);
		dEWKStat = sqrt(dEWKStat);
		dCosmicStat = sqrt(dCosmicStat);
		dHaloStat = sqrt(dHaloStat);
		*/
		double totStat = 0;
		if (iSideband == iWEIGHTED_SIDEBAND) totStat = sqrt(dQCDStat + dDiPhoStat+ dEWKStat + dCosmicStat +dHaloStat);
		else  totStat = sqrt(dPhoMcStat + dQCDStat + dDiPhoStat+ dEWKStat + dCosmicStat +dHaloStat);

		const double dTotSyst = sqrt(pow(dPhoMcSyst,2)  + pow(dQCDSyst,2)
											 + pow(dDiphoSyst,2) + pow(dEWKSyst,2)
											 + pow(dCosmicSyst,2)+ pow(dHaloSyst,2)
											 + pow(totStat,2));
		std::cout << "tot systs = " << dTotSyst << " / " << dTotSyst/fOrigBinWidth << std::endl;

		std::cout << "datat  / stat      = " << phojet->Integral() << " / " << sqrt(phojet->Integral()) << std::endl;
		std::cout << "Pho MC syst / stat = " << mcphojet->Integral() << " / " << dPhoMcSyst << " / " << sqrt (dPhoMcStat) << std::endl;
		std::cout << "QCD    syst / stat = " << qcdjet->Integral() << " / " << dQCDSyst << " / " << sqrt (dQCDStat) << std::endl;
		std::cout << "DIpho  syst / stat = " << diphojet->Integral() << " / " <<  dDiphoSyst << " / " << sqrt (dDiPhoStat) << std::endl;

		std::cout << "EWK    syst / stat = " << zeejet->Integral()+zmmjet->Integral()
															+ zttjet->Integral() + wenjet->Integral()
															+wmnjet->Integral() + wtnjet->Integral()
															<< " / " << dEWKSyst << " / " << sqrt (dEWKStat) << std::endl;
		std::cout << "Cosmic syst / stat = " << dCosmicSyst << " / " << sqrt (dCosmicStat) << std::endl;
		std::cout << "Halo   syst / stat = " << dHaloSyst << " / " << sqrt (dHaloStat) << std::endl;

		
		std::cout << "==================test 2 ====================" << std::endl;
		std::cout << "datat  / stat      = " << phojet->Integral("width") << " / " << sqrt(phojet->Integral("width")) << std::endl;
		std::cout << "Pho MC syst / stat = " << mcphojet->Integral("width") << " / " << dPhoMcSyst << " / " << sqrt (dPhoMcStat) << std::endl;
		std::cout << "QCD    syst / stat = " << qcdjet->Integral("width") << " / " << dQCDSyst << " / " << sqrt (dQCDStat) << std::endl;
		std::cout << "DIpho  syst / stat = " << diphojet->Integral("width") << " / " <<  dDiphoSyst << " / " << sqrt (dDiPhoStat) << std::endl;

		std::cout << "EWK    syst / stat = " << zeejet->Integral("width")+zmmjet->Integral("width")
															+ zttjet->Integral("width") + wenjet->Integral("width")
															+wmnjet->Integral("width") + wtnjet->Integral("width")
															<< " / " << dEWKSyst << " / " << sqrt (dEWKStat) << std::endl;
		std::cout << "Cosmic syst / stat = " << dCosmicSyst << " / " << sqrt (dCosmicStat) << std::endl;
		std::cout << "Halo   syst / stat = " << dHaloSyst << " / " << sqrt (dHaloStat) << std::endl;

/*	new TCanvas();

	TH1* hist_data_clone = dynamic_cast<TH1*> (phojet->Clone("phojet_data_clone"));
	TH1* hist_err_clone = dynamic_cast<TH1*> (hist_err->Clone("bg__clone"));
		
	hist_data_clone->Add(hist_err_clone, -1);
	hist_data_clone->SetTitle(sLegendTitle.c_str());
	hist_data_clone->GetYaxis()->SetTitleFont(iTITLE_FONT);
	hist_data_clone->GetYaxis()->SetLabelFont(iTITLE_FONT);
	hist_data_clone->GetYaxis()->CenterTitle(1);
	hist_data_clone->GetYaxis()->SetTitle("Data -- Background");
	hist_data_clone->GetXaxis()->CenterTitle(1);
	hist_data_clone->GetXaxis()->SetTitle(xtitle.c_str());

	hist_data_clone->Draw("PE");
	gPad->SetGridx();
	gPad->SetGridy();
		
*/

}

void MakeHistLogAndRatioV2 (int jets, std::string which, int icase=0)
//void MakeHistLogAndRatioV2 (int jets, std::string which)
{
  TVirtualPad *pad_old = gPad;	


	//SETUP ABSOLTE HIST PATHS FOR ALL HIST (INCLUDING SYSTEMATIC HISTS)
	ClearHistPaths();
	std::string sRelPath("");
	if (jets == 1)
	{
		if (which == "Et_pho") 				sRelPath +=  "1Jet/Photon/EtCorr";
		else if (which == "InvMass_pj1") sRelPath +=  "1Jet/PhotonLeadJet/InvMass";
		else if (which == "Ht") 			sRelPath +=  "1Jet/Event/Ht";
		else if (which == "Et_j1") 		sRelPath +=  "1Jet/LeadJet/EtCorr";
    	else if (which == "NJet") 			sRelPath +=  "1Jet/Event/NJet15";
    	else if (which == "Met") 			sRelPath +=  "1Jet/Event/Met";
    	else if (which == "DelR") 			sRelPath +=  "1Jet/PhotonLeadJet/DelR";
    	else if (which == "DelEta") 		sRelPath +=  "1Jet/PhotonLeadJet/DelEta";
    	else if (which == "DelPhi") 		sRelPath +=  "1Jet/PhotonLeadJet/DelPhi";
    	else if (which == "DetEta_pho") 	sRelPath +=  "1Jet/Photon/DelEta";
    	else if (which == "DetPhi_pho") 	sRelPath +=  "1Jet/Photon/DetPhi";
    	else if (which == "DetEta_j1") 	sRelPath +=  "1Jet/LeadJet/DetEta";
    	else if (which == "DetPhi_j1") 	sRelPath +=  "1Jet/LeadJet/DetPhi";
    	else if (which == "EvtPhi_j1") 	sRelPath +=  "1Jet/LeadJet/EvtPhi";

	} else if (jets == 2) {
		if (which == "Et_pho") 					sRelPath +=  "2Jet/Photon/EtCorr";
		else if (which == "InvMass_pj1") 	sRelPath +=  "2Jet/PhotonLeadJet/InvMass";
		else if (which == "InvMass_pj1j2") 	sRelPath +=  "2Jet/Photon2Jets/InvMass";
		else if (which == "Ht") 				sRelPath +=  "2Jet/Event/Ht";
		else if (which == "InvMass_j1j2") 	sRelPath +=  "2Jet/2Jets/InvMass";
		else if (which == "Et_j1") 			sRelPath +=  "2Jet/LeadJet/EtCorr";
		else if (which == "SecondLeadJetEt")sRelPath +=  "2Jet/SecondLeadJet/EtCorr";
    	else if (which == "Met") 				sRelPath +=  "2Jet/Event/Met";
    	else if (which == "DelR") 				sRelPath +=  "2Jet/PhotonLeadJet/DelR";
    	else if (which == "DelEta") 			sRelPath +=  "2Jet/PhotonLeadJet/DelEta";
    	else if (which == "DelPhi") 			sRelPath +=  "2Jet/PhotonLeadJet/DelPhi";
	}

	sDATA_HIST_PATH 	+= sDATA_DIR + sRelPath;  
	sHALO_HIST_PATH 	+= sHALO_DIR + sRelPath;  
	sCOSMIC_HIST_PATH 	+= sCOSMIC_DIR + sRelPath;  
	sSIDEBAND_HIST_PATH 	+= sSIDEBAND_DIR + sRelPath;  
	sSIDEBANDSYST_HADEM_HIST_PATH 	+= sSIDEBANDSYST_HADEM_DIR + sRelPath;  
	sSIDEBANDSYST_ISO_HIST_PATH 		+= sSIDEBANDSYST_ISO_DIR + sRelPath;
	sSIDEBANDSYST_TRKPT_HIST_PATH 	+= sSIDEBANDSYST_TRKPT_DIR + sRelPath;
	sSIDEBANDSYST_TRKISO_HIST_PATH 	+= sSIDEBANDSYST_TRKISO_DIR + sRelPath;
	sSIDEBAND_COSMIC_HIST_PATH 		+= sSIDEBAND_COSMIC_DIR + sRelPath;  
	sSIDEBAND_HALO_HIST_PATH 	+= sSIDEBAND_HALO_DIR + sRelPath;  
	sMCCENTRAL_HIST_PATH 		+= sMCCENTRAL_DIR + sRelPath;  
	sMCEMJESUP_HIST_PATH 		+= sMCEMJESUP_DIR + sRelPath;  
	sMCEMJESDOWN_HIST_PATH 		+= sMCEMJESDOWN_DIR + sRelPath;  
	sMCSIDEBAND_HIST_PATH 		+= sMCSIDEBAND_DIR + sRelPath;  



	//default bin sizes for any sample
	
	if (jets == 1)
	{
		if (which == "Et_pho")           { fXMIN=0  ; fXPOINT1 = 100 ; fXPOINT2 = 160; fXPOINT3 = 200;fXPOINT4 = 650 ; fWIDTH1= 10 ; fWIDTH2= 20 ; fWIDTH3= 40 ; fWIDTH4=250 ; iREBIN =5;} 
		//if (which == "Et_pho")           { fXMIN=0  ; fXPOINT1 = 150 ; fXPOINT2 = 200; fXPOINT3 = 300;fXPOINT4 = 650 ; fWIDTH1= 10 ; fWIDTH2= 10 ; fWIDTH3= 100 ; fWIDTH4= 250 ; iREBIN =5;} 
		else if (which == "Et_j1")       { fXMIN=0  ; fXPOINT1 = 100 ; fXPOINT2 = 200; fXPOINT3 = 250;fXPOINT4 = 600 ; fWIDTH1= 20 ; fWIDTH2= 50; fWIDTH3= 50 ; fWIDTH4= 200; iREBIN =5;}
		else if (which == "InvMass_pj1") { fXMIN=50 ; fXPOINT1 = 250 ; fXPOINT2 = 350; fXPOINT3 = 550;fXPOINT4 = 1200; fWIDTH1= 25 ; fWIDTH2= 50 ; fWIDTH3= 100; fWIDTH4= 350; iREBIN =10;}
		//else if (which == "InvMass_pj1") { fXMIN=0 ; fXPOINT1 = 500 ; fXPOINT2 = 600; fXPOINT3 = 700;fXPOINT4 = 1000; fWIDTH1= 25 ; fWIDTH2= 50 ; fWIDTH3= 75 ; fWIDTH4= 300 ; iREBIN =5;}
		//else if (which == "InvMass_pj1") { fXMIN=0 ; fXPOINT1 = 500 ; fXPOINT2 = 600; fXPOINT3 = 700;fXPOINT4 = 1000; fWIDTH1= 5 ; fWIDTH2= 5 ; fWIDTH3= 5 ; fWIDTH4= 5 ; iREBIN =5;}
		else if (which == "Ht")          { fXMIN=0  ; fXPOINT1 = 200 ; fXPOINT2 = 400; fXPOINT3 = 600;fXPOINT4 = 1200; fWIDTH1= 40 ; fWIDTH2= 50 ; fWIDTH3=200; fWIDTH4= 450; iREBIN =5;}
    	else if (which == "NJet")        { fXMIN=0  ; fXPOINT1 = 1   ; fXPOINT2 = 2  ; fXPOINT3 = 3  ;fXPOINT4 = 15  ; fWIDTH1= 1  ; fWIDTH2= 1  ; fWIDTH3= 1  ; fWIDTH4= 1   ; iREBIN =1;}
    	else if (which == "Met")         { fXMIN=0  ; fXPOINT1 = 60  ; fXPOINT2 = 100; fXPOINT3 = 150;fXPOINT4 = 350 ; fWIDTH1= 10 ; fWIDTH2= 20 ; fWIDTH3= 50 ; fWIDTH4=150 ; iREBIN =1;}
    	else if (which == "DelR")        { fXMIN=0  ; fXPOINT1 = 1   ; fXPOINT2 = 2  ; fXPOINT3 = 4  ;fXPOINT4 = 6   ; fWIDTH1= 0.2; fWIDTH2= 0.2; fWIDTH3= 0.2; fWIDTH4= 0.2 ; iREBIN =5;}
    	else if (which == "DelEta")      { fXMIN=0  ; fXPOINT1 = 1   ; fXPOINT2 = 3  ; fXPOINT3 = 4  ;fXPOINT4 = 5   ; fWIDTH1= 0.1; fWIDTH2= 0.1; fWIDTH3= 0.1; fWIDTH4= 0.1 ; iREBIN =5;}
    	else if (which == "DelPhi")      { fXMIN=0  ; fXPOINT1 = 1   ; fXPOINT2 = 2  ; fXPOINT3 = 3  ;fXPOINT4 = 4   ; fWIDTH1= 0.1; fWIDTH2= 0.1; fWIDTH3= 0.1; fWIDTH4= 0.1 ; iREBIN =5;}
    	else if (which == "DetEta_pho")  { fXMIN=-1.5; fXPOINT1 = -0.5; fXPOINT2 = 0.5; fXPOINT3 = 1  ;fXPOINT4 = 1.5 ; fWIDTH1= 0.1; fWIDTH2= 0.1; fWIDTH3= 0.1; fWIDTH4= 0.1 ; iREBIN =1;}
    	else if (which == "DetPhi_pho")  { fXMIN=0  ; fXPOINT1 = 1   ; fXPOINT2 = 2  ; fXPOINT3 = 3  ;fXPOINT4 = 8   ; fWIDTH1= 0.1; fWIDTH2= 0.1; fWIDTH3= 0.1; fWIDTH4= 0.1 ; iREBIN =1;}
    	else if (which == "DetEta_j1")   { fXMIN=-3.5; fXPOINT1 = -1.5; fXPOINT2 = 0.5; fXPOINT3 = 1.5  ;fXPOINT4 = 3.5 ; fWIDTH1= 0.2; fWIDTH2= 0.2; fWIDTH3= 0.2; fWIDTH4= 0.2 ; iREBIN =1;}
    	//else if (which == "DetEta_j1")   { fXMIN=-3.5; fXPOINT1 = -1.5; fXPOINT2 = 0.5; fXPOINT3 = 1.5  ;fXPOINT4 = 3.5 ; fWIDTH1= 0.1; fWIDTH2= 0.1; fWIDTH3= 0.1; fWIDTH4= 0.1 ; iREBIN =1;}
    	else if (which == "DetPhi_j1")   { fXMIN=0  ; fXPOINT1 = 1   ; fXPOINT2 = 2  ; fXPOINT3 = 3  ;fXPOINT4 = 4   ; fWIDTH1= 0.1; fWIDTH2= 0.1; fWIDTH3= 0.1; fWIDTH4= 0.1 ; iREBIN =1;}
    	else if (which == "EvtPhi_j1")   { fXMIN=0  ; fXPOINT1 = 1   ; fXPOINT2 = 2  ; fXPOINT3 = 3  ;fXPOINT4 = 4   ; fWIDTH1= 0.1; fWIDTH2= 0.1; fWIDTH3= 0.1; fWIDTH4= 0.1 ; iREBIN =1;}
	} else if (jets == 2) {
		if (which == "Et_pho") 				 	{ fXMIN=0  ; fXPOINT1= 100; fXPOINT2= 200; fXPOINT3= 450; fXPOINT4= 600 ; fWIDTH1=10 ; fWIDTH2=50 ; fWIDTH3=150; fWIDTH4=150; iREBIN =5;}
		else if (which == "Et_j1") 		 	{ fXMIN=0  ; fXPOINT1= 100; fXPOINT2= 200; fXPOINT3= 300; fXPOINT4= 600 ; fWIDTH1=10 ; fWIDTH2=50 ; fWIDTH3=100; fWIDTH4=200; iREBIN =5;}
		//else if (which == "InvMass_pj1")  	{ fXMIN=100; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 500; fXPOINT4= 1100; fWIDTH1=20 ; fWIDTH2=50 ; fWIDTH3=100 ; fWIDTH4=500 ; iREBIN =10;}
		else if (which == "InvMass_pj1")  	{ fXMIN=50; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 500; fXPOINT4= 1100; fWIDTH1=20 ; fWIDTH2=50 ; fWIDTH3=100 ; fWIDTH4=500 ; iREBIN =10;}
		//else if (which == "InvMass_pj1j2")	{ fXMIN=100; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 600; fXPOINT4= 1300; fWIDTH1=10 ; fWIDTH2=20 ; fWIDTH3=50 ; fWIDTH4=500 ; iREBIN =10;}
		else if (which == "InvMass_pj1j2")	{ fXMIN=80; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 600; fXPOINT4= 1300; fWIDTH1=10 ; fWIDTH2=20 ; fWIDTH3=50 ; fWIDTH4=500 ; iREBIN =10;}
		else if (which == "InvMass_j1j2") 	{ fXMIN=60; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 600; fXPOINT4= 1200; fWIDTH1=20 ; fWIDTH2=50 ; fWIDTH3=100 ; fWIDTH4=500 ; iREBIN =10;}
		//else if (which == "InvMass_j1j2") 	{ fXMIN=0; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 600; fXPOINT4= 1200; fWIDTH1=20 ; fWIDTH2=50 ; fWIDTH3=100 ; fWIDTH4=500; iREBIN=1;}
		
		/*else if (which == "InvMass_pj1")  	{ fXMIN=0; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 700; fXPOINT4= 1000; fWIDTH1=10 ; fWIDTH2=20 ; fWIDTH3=50 ; fWIDTH4=50 ; iREBIN =5;}
		else if (which == "InvMass_pj1j2")	{ fXMIN=0; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 700; fXPOINT4= 1300; fWIDTH1=10 ; fWIDTH2=20 ; fWIDTH3=50 ; fWIDTH4=200 ; iREBIN =5;}
		else if (which == "InvMass_j1j2") 	{ fXMIN=0; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 600; fXPOINT4= 1200; fWIDTH1=20 ; fWIDTH2=20 ; fWIDTH3=50 ; fWIDTH4=50 ; iREBIN =5;}
		*/
		else if (which == "Ht")				 	{ fXMIN=0  ; fXPOINT1= 50; fXPOINT2= 400; fXPOINT3= 600; fXPOINT4= 1200; fWIDTH1=40 ; fWIDTH2=50 ; fWIDTH3=200; fWIDTH4=450 ; iREBIN =5;}
    	else if (which == "Met")			 	{ fXMIN=0  ; fXPOINT1= 50 ; fXPOINT2= 150; fXPOINT3= 250; fXPOINT4= 300 ; fWIDTH1=10 ; fWIDTH2=100 ; fWIDTH3=100; fWIDTH4=50; iREBIN =1;}
		else if (which == "SecondLeadJetEt"){ fXMIN=0  ; fXPOINT1= 100; fXPOINT2= 160; fXPOINT3= 200; fXPOINT4= 600 ; fWIDTH1=10 ; fWIDTH2=20 ; fWIDTH3=40 ; fWIDTH4=200 ; iREBIN =5;}
    	else if (which == "DelR") 				{ fXMIN=0  ; fXPOINT1= 1  ; fXPOINT2= 2  ; fXPOINT3= 4  ; fXPOINT4= 6   ; fWIDTH1=0.2; fWIDTH2=0.2; fWIDTH3=0.2; fWIDTH4=0.2 ; iREBIN =1;}
    	else if (which == "DelEta") 			{ fXMIN=0  ; fXPOINT1= 1  ; fXPOINT2= 3  ; fXPOINT3= 4  ; fXPOINT4= 5   ; fWIDTH1=0.1; fWIDTH2=0.1; fWIDTH3=0.1; fWIDTH4=0.1 ; iREBIN =1;}
    	else if (which == "DelPhi") 			{ fXMIN=0  ; fXPOINT1= 1  ; fXPOINT2= 2  ; fXPOINT3= 3  ; fXPOINT4= 4   ; fWIDTH1=0.1; fWIDTH2=0.1; fWIDTH3=0.1; fWIDTH4=0.1 ; iREBIN =1;}
	} else 
	{
		std::cout << __FUNCTION__ << ":: Invalid Njets request! only 1 or 2 jets can be give." << std::endl;
	}
	




	/*  combination of settings that would work
	 * 
	 *                                             qcdMcMix      iSample 									iSideband
	 * <1> default method (30%qcd/70%MC PHO)          2        iG30/iG40/iG30MET25/iG30MET40 		     iNOMINAL_SIDEBAND
	 *             but for njet plot                  1           ""                  						""
	 *
	 *
	 * <2> reweighed method (100% sideband)           1           ""               						iWEIGHTED_SIDEBAND
	 *
	 */

  
 	//int qcdMcMix = 1; int iSample = iG30; int iSideband = iWEIGHTED_SIDEBAND; iMTDB_RERUN=0;sPRINTFOLDER="./eps/";
	int qcdMcMix = 2; int iSample = iG30; int iSideband = iNOMINAL_SIDEBAND; iMTDB_RERUN=0; iHAVE_MTDB_ERRORS=0;sPRINTFOLDER="./eps/";
	//int qcdMcMix = 2; int iSample = iG40; int iSideband = iNOMINAL_SIDEBAND; sPRINTFOLDER="~/TMP/eps/g40nominal/";
	//int qcdMcMix = 1; int iSample = iG40; int iSideband = iWEIGHTED_SIDEBAND; sPRINTFOLDER="~/TMP/eps/g40weighed/";
	//int qcdMcMix = 2; int iSample = iG30MET25; int iSideband = iNOMINAL_SIDEBAND; sPRINTFOLDER="~/TMP/eps/g30met25nominal/";
	//int qcdMcMix = 2; int iSample = iG30MET20; int iSideband = iNOMINAL_SIDEBAND; iMTDB_RERUN=0; iHAVE_MTDB_ERRORS=0;sPRINTFOLDER="./eps/";
	//int qcdMcMix = 1; int iSample = iG30MET20; int iSideband = iWEIGHTED_SIDEBAND; iMTDB_RERUN=0; sPRINTFOLDER="./eps/";
	//int qcdMcMix = 2; int iSample = iG30MET30; int iSideband = iNOMINAL_SIDEBAND; iMTDB_RERUN=0; sPRINTFOLDER="~/TMP/eps/g30met30nominal/";
	//int qcdMcMix = 1; int iSample = iG30MET30; int iSideband = iWEIGHTED_SIDEBAND; iMTDB_RERUN=1; sPRINTFOLDER="~/TMP/eps/g30met30nominal/";
	//int qcdMcMix = 2; int iSample = iG30MET40; int iSideband = iNOMINAL_SIDEBAND; sPRINTFOLDER="~/TMP/eps/g30met40nominal/";
	//int qcdMcMix = 1; int iSample = iG30MET40; int iSideband = iWEIGHTED_SIDEBAND; sPRINTFOLDER="~/TMP/eps/g30met40weighed/";
	//int qcdMcMix = 2; int iSample = iG30EXCL1J; int iSideband = iNOMINAL_SIDEBAND; sPRINTFOLDER="./";
	//int qcdMcMix = 1; int iSample = iG30EXCL1J; int iSideband = iWEIGHTED_SIDEBAND; sPRINTFOLDER="./";
	//int qcdMcMix = 2; int iSample = iG30MET50; int iSideband = iNOMINAL_SIDEBAND; iMTDB_RERUN=0; iHAVE_MTDB_ERRORS=0; sPRINTFOLDER="~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_Met50_final/eps/";
	//int qcdMcMix = 1; int iSample = iG30MET50; int iSideband = iWEIGHTED_SIDEBAND; iMTDB_RERUN=1; sPRINTFOLDER="~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_Met50_final/eps/";

	if (icase==2)
	{
		qcdMcMix = 1;
		iSideband = iWEIGHTED_SIDEBAND;
	}
	
	//if (which == "NJet") qcdMcMix=1;
	//if (which == "Met") qcdMcMix=1;
	
	if (iSample == iG40) sPRINTFOLDER = "~/ICHEP_LaptopVersion/PhoEt40Files/eps/";
	if (iSample == iG30MET25) sPRINTFOLDER = "~/ICHEP_LaptopVersion/PhoEt30Met25/eps/";
//	if (iSample == iG30MET20) sPRINTFOLDER = "~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_Met20_final/eps/";
	if (iSample == iG30MET30) sPRINTFOLDER = "~/ICHEP_LaptopVersion/p1_26/NvtxWgtAndMetCleaned_Met30_final/eps/";
	if (iSample == iG30MET40) sPRINTFOLDER = "~/ICHEP_LaptopVersion/PhoEt30Met40/eps/";


	/**********************************************************
	 * modify the default bin sizes for specific samples here
	 **********************************************************/
	
	if (iSample == iG30MET20)
	{
		if (jets == 1)
		{
			if (which == "Et_pho")           { fXMIN=0  ; fXPOINT1 = 100 ; fXPOINT2 = 150; fXPOINT3 = 200;fXPOINT4 = 650 ; fWIDTH1= 10 ; fWIDTH2= 20 ; fWIDTH3= 50 ; fWIDTH4=250 ; iREBIN =5;} 
			if (which == "Et_j1")            { fXMIN=0  ; fXPOINT1 = 100 ; fXPOINT2 = 150; fXPOINT3 = 200;fXPOINT4 = 600 ; fWIDTH1= 20 ; fWIDTH2= 50; fWIDTH3= 50 ; fWIDTH4= 250; iREBIN =5;}
			if (which == "InvMass_pj1")      { fXMIN=50 ; fXPOINT1 = 250 ; fXPOINT2 = 350; fXPOINT3 = 450;fXPOINT4 = 1200; fWIDTH1= 25 ; fWIDTH2= 50 ; fWIDTH3= 100; fWIDTH4= 450; iREBIN =5;}
		}

		if (jets == 2) 
		{
			if (which == "Et_pho") 		 	{ fXMIN=0  ; fXPOINT1= 100; fXPOINT2= 200; fXPOINT3= 450; fXPOINT4= 500 ; fWIDTH1=10 ; fWIDTH2=50 ; fWIDTH3=250; fWIDTH4=50; iREBIN =5;}
			if (which == "Et_j1") 		 	{ fXMIN=0  ; fXPOINT1= 100; fXPOINT2= 200; fXPOINT3= 400; fXPOINT4= 450 ; fWIDTH1=10 ; fWIDTH2=50 ; fWIDTH3=200; fWIDTH4=50; iREBIN =5;}
			if (which == "InvMass_pj1j2")	{ fXMIN=100; fXPOINT1= 200; fXPOINT2= 600; fXPOINT3= 1000; fXPOINT4= 1100; fWIDTH1=50 ; fWIDTH2=100 ; fWIDTH3=400 ; fWIDTH4=100 ; iREBIN =5;}
			if (which == "InvMass_j1j2") 	{ fXMIN=100; fXPOINT1= 200; fXPOINT2= 600; fXPOINT3= 800; fXPOINT4= 850; fWIDTH1=50 ; fWIDTH2=100 ; fWIDTH3=200 ; fWIDTH4=50 ; iREBIN =5;}
		}
	}


	if (iSample == iG30MET40)
	{
		sPRINTFOLDER = "~/ICHEP_LaptopVersion/PhoEt30Met40/eps/";
		
		if (jets == 1)
		{
			if (which == "Et_pho")
			{
				if (qcdMcMix == 2) { fXMIN=0 ; fXPOINT1=100; fXPOINT2=200; fXPOINT3=300; fXPOINT4=700;  fWIDTH1=10; fWIDTH2=20; fWIDTH3=50 ; fWIDTH4=200;}
				if (qcdMcMix == 1) { fXMIN=30 ; fXPOINT1=100; fXPOINT2=160; fXPOINT3=200; fXPOINT4=600;  fWIDTH1=10; fWIDTH2=20; fWIDTH3=40 ; fWIDTH4=300;}
			}
		 	else if (which == "Et_j1") { fXMIN=10; fXPOINT1=90 ; fXPOINT2=150; fXPOINT3=250; fXPOINT4=550;  fWIDTH1=10; fWIDTH2=20; fWIDTH3=100; fWIDTH4=200;}
			else if (which == "Ht")    { fXMIN=0 ; fXPOINT1=200; fXPOINT2=400; fXPOINT3=700; fXPOINT4=1200; fWIDTH1=40; fWIDTH2=40; fWIDTH3=100; fWIDTH4=420;}
			else if (which == "InvMass_pj1") 
			{ 
				if (qcdMcMix == 1) { fXMIN=50 ; fXPOINT1 = 300 ; fXPOINT2 = 400; fXPOINT3 = 550;fXPOINT4 = 900; fWIDTH1= 25 ; fWIDTH2= 50 ; fWIDTH3= 50 ; fWIDTH4=270 ;}
			}
				
		}

		
		//2 JETS
		if (jets == 2)
		{
			if (which == "Et_pho")		{ fXMIN=0; fXPOINT1=100; fXPOINT2=160; fXPOINT3=200; fXPOINT4=450; fWIDTH1=10; fWIDTH2=20; fWIDTH3=40; fWIDTH4=220;}
			else if (which == "Et_j1") 		 	{ fXMIN=14 ; fXPOINT1= 104; fXPOINT2= 264; fXPOINT3= 300; fXPOINT4= 600 ; fWIDTH1=10 ; fWIDTH2=20 ; fWIDTH3=36 ; fWIDTH4=200 ;}
			else if (which == "InvMass_pj1")  	{ fXMIN=100; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 500; fXPOINT4= 800; fWIDTH1=20 ; fWIDTH2=20 ; fWIDTH3=100 ; fWIDTH4=200 ;}
			else if (which == "InvMass_pj1j2")	{ fXMIN=100; fXPOINT1= 200; fXPOINT2= 400; fXPOINT3= 600; fXPOINT4= 1100; fWIDTH1=20 ; fWIDTH2=50 ; fWIDTH3=100 ; fWIDTH4=400 ;}
			else if (which == "InvMass_j1j2") 	{ fXMIN=100; fXPOINT1= 200; fXPOINT2= 300; fXPOINT3= 500; fXPOINT4= 900; fWIDTH1=20 ; fWIDTH2=50 ; fWIDTH3=200 ; fWIDTH4=350 ;}
			else if (which == "Ht")				 	{ fXMIN=100; fXPOINT1= 200; fXPOINT2= 400; fXPOINT3= 600; fXPOINT4= 1200; fWIDTH1=20 ; fWIDTH2=50 ; fWIDTH3=100; fWIDTH4=520 ;}
		}
	}

	
	if ( bDO_VAR_BINNING )
	{
		if ( (fXMIN > fXPOINT1) || (fXPOINT1 > fXPOINT2) || (fXPOINT3 > fXPOINT4) )
		{
			cout << __FUNCTION__ << ":" << __LINE__ << ": xpoints must be monotonically increasing!" << std::endl;
			return;
		}
		if (fWIDTH1 <= 0 || fWIDTH2 <= 0 || fWIDTH3 <= 0 || fWIDTH4 <= 0) 
		{
			cout << __FUNCTION__ << ":" << __LINE__ << ": bin widths must be >0 !" << std::endl;
			return;
		}
	}
	
	// to make full scale plots
	if (jets == 1) {
		if (which == "Et_pho") 				MakeHistLogAndRatioV2 (which, jets, "EtCorr","Photon E_{T} [GeV/c]","1Jet/Photon", 											qcdMcMix, "pj1_Et_pho", 	iSample, iSideband, "mBerr_pj1_Et_pho");
		else if (which == "InvMass_pj1") MakeHistLogAndRatioV2 (which, jets, "InvMass","Invariant Mass (#gamma,Lead Jet) [GeV/c^{2}]","1Jet/PhotonLeadJet", 	qcdMcMix, "pj1_InvM_pj1",iSample, iSideband, "mBerr_pj1_InvM_pj1");
		else if (which == "Ht") 			MakeHistLogAndRatioV2 (which, jets, "Ht","H_{T} [GeV/c]","1Jet/Event",																qcdMcMix, "pj1_Ht", 		iSample, iSideband, "mBerr_pj1_Ht");
		else if (which == "Et_j1") 		MakeHistLogAndRatioV2 (which, jets, "EtCorr","E_{T}^{Lead Jet} [GeV/c]","1Jet/LeadJet", 										qcdMcMix, "pj1_Et_j1", 	iSample, iSideband, "mBerr_pj1_Et_j1");
    	else if (which == "NJet") 			MakeHistLogAndRatioV2 (which, jets, "NJet15","Jet Multiplicity (E_{T}>15 GeV)", "1Jet/Event",								qcdMcMix, "", iSample, iSideband, "mBerr_pj1_Njet");
    	else if (which == "Met") 			MakeHistLogAndRatioV2 (which, jets, "Met","#slash{E}_{T} [GeV/c]", "1Jet/Event",										qcdMcMix, "", iSample, iSideband ,"mBerr_pj1_Met");
    	else if (which == "DelR") 			MakeHistLogAndRatioV2 (which, jets, "DelR","#Delta R (#gamma, lead jet)", "1Jet/PhotonLeadJet",							qcdMcMix, "", iSample, iSideband);
    	else if (which == "DelEta") 		MakeHistLogAndRatioV2 (which, jets, "DelEta","#Delta #eta (#gamma, lead jet)", "1Jet/PhotonLeadJet",		 			qcdMcMix, "", iSample, iSideband);
    	else if (which == "DelPhi") 		MakeHistLogAndRatioV2 (which, jets, "DelPhi","#Delta #phi (#gamma, lead jet)", "1Jet/PhotonLeadJet", 					qcdMcMix, "", iSample, iSideband);
    	else if (which == "DetEta_pho") 	MakeHistLogAndRatioV2 (which, jets, "DetEta","#gamma detector #eta", "1Jet/Photon",                  					qcdMcMix, "", iSample, iSideband);
    	else if (which == "DetPhi_pho") 	MakeHistLogAndRatioV2 (which, jets, "DetPhi","#gamma detector #phi", "1Jet/Photon",											qcdMcMix, "", iSample, iSideband);
    	else if (which == "DetEta_j1") 	MakeHistLogAndRatioV2 (which, jets, "DetEta","Lead jet detector #eta", "1Jet/LeadJet",										qcdMcMix, "", iSample, iSideband);
    	else if (which == "DetPhi_j1") 	MakeHistLogAndRatioV2 (which, jets, "DetPhi","Lead jet detector #phi", "1Jet/LeadJet",										qcdMcMix, "", iSample, iSideband);
    	else if (which == "EvtPhi_j1") 	MakeHistLogAndRatioV2 (which, jets, "EvtPhi","Lead jet 4-vec #phi", "1Jet/LeadJet",											qcdMcMix, "", iSample, iSideband);
	}

	if (jets == 2) {
		if (which == "Et_pho") 					MakeHistLogAndRatioV2 (which, jets, "EtCorr","Photon E_{T} [GeV/c]","2Jet/Photon", 												qcdMcMix,"pj2_Et_pho", 		iSample, iSideband, "mBerr_pj2_Et_pho");
		else if (which == "InvMass_pj1") 	MakeHistLogAndRatioV2 (which, jets, "InvMass","Invariant Mass (#gamma,Lead Jet) [GeV/c^{2}]","2Jet/PhotonLeadJet", 	qcdMcMix,"pj2_InvM_pj1", 	iSample, iSideband, "mBerr_pj2_InvM_pj1");
		else if (which == "InvMass_pj1j2") 	MakeHistLogAndRatioV2 (which, jets, "InvMass","Invariant Mass (#gamma, Two Lead Jets) [GeV/c^{2}]","2Jet/Photon2Jets", qcdMcMix,"pj2_InvM_pj1j2", 	iSample, iSideband, "mBerr_pj2_InvM_pj1j2");
		else if (which == "Ht") 				MakeHistLogAndRatioV2 (which, jets, "Ht","H_{T} (GeV/c)","2Jet/Event", 																qcdMcMix,"pj1_Ht", 			iSample, iSideband, "mBerr_pj2_Ht");
		else if (which == "InvMass_j1j2") 	MakeHistLogAndRatioV2 (which, jets, "InvMass","Invariant Mass (Two Lead Jets) [GeV/c^{2}]","2Jet/2Jets", 				qcdMcMix,"pj2_InvM_j1j2", 	iSample, iSideband, "mBerr_pj2_InvM_j1j2");
		else if (which == "Et_j1") 			MakeHistLogAndRatioV2 (which, jets, "EtCorr","E_{T}^{Lead Jet} [GeV/c]","2Jet/LeadJet", 											qcdMcMix,"pj2_Et_j1", 		iSample, iSideband, "mBerr_pj2_Et_j1");
		else if (which == "SecondLeadJetEt")MakeHistLogAndRatioV2 (which, jets, "EtCorr","E_{T}^{Second Lead Jet} [GeV/c]","2Jet/SecondLeadJet", 						qcdMcMix,"", iSample, iSideband);
    	else if (which == "Met") 				MakeHistLogAndRatioV2 (which, jets, "Met","#slash{E}_{T} [GeV/c]", "2Jet/Event",											qcdMcMix,"", iSample, iSideband, "mBerr_pj2_Met");
    	else if (which == "DelR") 				MakeHistLogAndRatioV2 (which, jets, "DelR","#Delta R (#gamma, lead jet)", "2Jet/PhotonLeadJet",								qcdMcMix,"", iSample, iSideband);
    	else if (which == "DelEta") 			MakeHistLogAndRatioV2 (which, jets, "DelEta","#Delta #eta (#gamma, lead jet)", "2Jet/PhotonLeadJet",						qcdMcMix,"", iSample, iSideband);
    	else if (which == "DelPhi") 			MakeHistLogAndRatioV2 (which, jets, "DelPhi","#Delta #phi (#gamma, lead jet)", "2Jet/PhotonLeadJet",						qcdMcMix,"", iSample, iSideband);
	}
	

	if (gPad != pad_old)
	{
		TCanvas *c = dynamic_cast<TCanvas*>(gPad);
		//c->SetFillColor(kRed);
      if (c)
		{
			std::ostringstream str,str1,str2;
			if (iSideband == iWEIGHTED_SIDEBAND)
			{
				str2 << sPRINTFOLDER << sSAMPLE_NAME << "_MtdB_plot" << jets << "_" << which  << "." << sPRINTFORMAT.c_str();
			} else {
				str2 << sPRINTFOLDER << sSAMPLE_NAME << "_MtdA_plot" << jets << "_" << which  << "." << sPRINTFORMAT.c_str();
			}
				
				std::cout << "Print Folder = " << sPRINTFOLDER << std::endl;
			c->Print(str2.str().c_str());
		};
	};

};

void MakeHistLogAndRatioV2 (const int i)
{
	//good plots that works
	if (i == 1) MakeHistLogAndRatioV2 (1, "Et_pho");
	if (i == 2) MakeHistLogAndRatioV2 (1, "Et_j1");
	if (i == 3) MakeHistLogAndRatioV2 (1, "InvMass_pj1");
	if (i == 4) MakeHistLogAndRatioV2 (1, "Ht");
	if (i == 5) MakeHistLogAndRatioV2 (2, "Et_pho");
	if (i == 6) MakeHistLogAndRatioV2 (2, "Et_j1");
	if (i == 7) MakeHistLogAndRatioV2 (2, "InvMass_pj1");
	if (i == 8) MakeHistLogAndRatioV2 (2, "InvMass_pj1j2");
	if (i == 9) MakeHistLogAndRatioV2 (2, "InvMass_j1j2");
	
	if (i== 10) MakeHistLogAndRatioV2 (1, "NJet");
	if (i== 11) MakeHistLogAndRatioV2 (2, "Ht");
	if (i== 12) MakeHistLogAndRatioV2 (1, "Met");
	if (i== 13) MakeHistLogAndRatioV2 (2, "Met");
	if (i== 14) MakeHistLogAndRatioV2 (1, "DetEta_pho");
	if (i== 15) MakeHistLogAndRatioV2 (1, "DetEta_j1");
	if (i== 16) MakeHistLogAndRatioV2 (1, "DetR_pj1");


	//MakeHistLogAndRatioV2 (1, "Met");
	//MakeHistLogAndRatioV2 (2, "SecondLeadJetEt");
	//MakeHistLogAndRatioV2 (1, "InvMass");
	//MakeHistLogAndRatioV2 (2, "InvMass");
	//MakeHistLogAndRatioV2 (2, "PhoJetsInvMass");
};

void MakeHistLogAndRatioV2 ()
{
	for (int i=1;i<=16;++i) MakeHistLogAndRatioV2 (i);
}

//list current temporary edit

