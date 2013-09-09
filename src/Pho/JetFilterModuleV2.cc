/*****************************************************************************/
/*	This is a modified version of TMyJetFilterModule. It is to be used in     */
/*                                 	g+Jets                                   */
/*****************************************************************************/


/*
 * $Id: JetFilterModuleV2.cc,v 1.65 2011/05/25 20:27:24 samantha Exp $
 * $Log: JetFilterModuleV2.cc,v $
 * Revision 1.65  2011/05/25 20:27:24  samantha
 * Included #include <assert.h>
 *
 * Revision 1.64  2009/11/10 19:55:00  samantha
 * DELETED:1. Unused iSavedEvents variable.
 * 	2. Unused SmearedHtAll(), GetGenLevelInfo(),DebugInfo()
 * UNCOMMENTED: 1. MyHtALL() as I am no longer using MET model, just need correct
 * HT and MET and the jets.
 * ADDED  :1  Another mode that produces only the corrected jets and MET/HT for me.
 * 	2. Accessor to get lev7 no EM jets.
 * 	3. moved GetMatch_JetIndex() from hh
 * Minor modificcation trhough out the code, commenting out debug info. Also to the
 * end job summary with the new mode.
 *
 * Revision 1.63  2009/08/25 16:10:42  samantha
 * ADDED: 1: But not complete, several multi-dimentional std::vectors to store JER
 *        parameters and its errors. I have made them static. The usage is
 *        somewhat unorganzied. The parameter values are under a funtion that does
 *        not belong to the class, and will be used to fill the vector at the
 *        beginning of the job,which is fine. Then when we try to reconstruct the
 *        JER from these, the JER function which is not a member of the class, uses
 *        public access to the parameters. I was able test and 2D vectors. But ran
 *        into trouble with 3D vector, where I wanted to store the exact JER bin
 *        values for E<60GeV. So all this is not used for now. Compilor complains
 *        the 3D vector is not defined! But I have defined it exactly the same way
 *        as the 2D vectors. I doubt there is limit on the dimensionality of
 *        std::vector. For all it see is, just some object, in this case another
 *        vector! So all these vectors a commented still.
 * MODIFIED:1. JERparam() is modified by commenting out the old parameterization for
 * 	the full range, and adding the parameterization only for E>60GeV. For
 * 	E<60GeV region, we use exact bin values. See reults in Elog#1287. JER
 * 	fits and things are in #1286.
 * 	2. I have replaced the way we access std::vector elements. Instead of
 * 	using array notation I use the bound safe method 'at()'.
 * 	3. Lots of prints in InitMetProbStuff() to understand the problem I had
 * 	with JER function collapsing and MetSig getting NAN values. They are
 * 	commented out for now.
 *
 * Revision 1.62  2009/08/07 21:21:48  samantha
 * ADDED: 1: new2JERparam() and new3JERparam to hold 2nd and 3d(final) round of fit parameters.
 * Though I can quickly switch between old and new fit functions, I have to manually
 * edit the file to use between new2Jer and new3Jer fits. But this will be not necessay once we
 * decide which of these we'll keep to do our systematics.
 * Note on new3JERparam():: The top list of parameters are the most recent once and comes from
 * elog#1275 plots. The once below it (commmented out) are previous to this are in elog#1274.
 * The difference between the two is that i have limited the floating of the crossing point of
 * the fit function in this most recent fits. This is in an attempt to overcome the predicament
 * mentioned in elog#1270. Almost worked, but needed another tweak. See MODIFIED: 1.
 *
 *         2: FitFuncForMeanMPVNorm() and FitFuncForSigmas() with this newJER3.
 * FitFuncForMeanMPVNorm() is for Gauss Mean, Landau MPV and Gauss Norm.
 *         3: fitFunc3_new() is using fitFunc3() but stitches it between two region using the
 * two step functions. Used for 1st (or 2nd) roound of fitting for the sigams.
 * 	4:Filling the new Photons energy resolution plot by finding a matching HEPG photon.
 * 	Method is executed only in MC mode.
 * 	5: Few prints statements hangs in the InitMetProbStuff to understand and remedy the
 * 	case when the JER integral returns 0! see MODFIED:1 for more details.
 *
 *
 * FitFuncForSigmas() is for Landau and Gauss sigmas.
 * MODIFIED: 1: MyJER() so that in case of very small Landau Sigma, its contribution
 * is set to zero. This is part of the solution to the problem I encountered as
 * explianed in elog#1270. But now occasionally the integral returns 0 and the 'test_prob'
 * in InitMetProbStuff(), gets an undefined value.This causes  stuff.max_MetProbInteg
 * to get un undefined values. Though this happens in about 10 events, I am need to fix
 * it. Other question I have is, even if the Landau contribution is zero shouldn't there
 * be non-zero Gaussian contribution and hence the integral should not be zero!
 *            2: Minor spacing adjustments. Tryin to keep things within 80 columns.
 *
 * Revision 1.61  2009/07/24 23:00:41  samantha
 * ADDED: 1: New JER fit method. I nicely stictched the code access the new JER
 * fits
 * with minimal modifications. The new fit function is implemented in very compact
 * way. Now I can actually mix old and new JER fits for the 5 parameters (Gauss
 * Mean,
 * Gauss Sigma, Landau MPV, Landau Sigma, and Gauss Norm).
 *
 * Revision 1.60  2009/07/21 15:45:38  samantha
 * MODIFIED: 1: JERparam() is modified with new JER parameterization. So far only
 * Gauss Mean, Landau Mean, and Gauss Norm is fitted with the new function,
 * dJayFunc1. Gauss Sigma and Landau Sigma still uses the old fit function (Sasha
 * said these two must have physics motivated fit functions, a function that
 * represents calorimeter behaviour. The new function has no such thing! But Jay
 * disagrees.) See elog#1243(New JER fits), 1246 (Final PHOJETC MC hists from this new fit
 * parameters), 1249(same results for DiPhoX)
 *           2: JER generation (or getting parameter values for JER) points is
 * 	  modified to use new dJayFunc1 for Gauss Mean, Landau Mean, and Gauss
 * 	  Norm.
 * DELETED: 1: Xsection functions/hists and VTX swapping stuff.
 *
 * Revision 1.59  2009/06/26 22:57:09  samantha
 * MAJOR CHANGES: 1: In GenerateMyTotalMet() reintroduced the njet cut. So now the
 * MEt is generated using smeared jet only if there is atleast one 15GeV jet. So is
 * the SumEt. In the cases of events without a 15GeV jet, only the uncluster MEt is
 * generated as total generated MEt. I had to move the njet from the original
 * postion so the smear factors are generated for all the events.
 * ADDED: 	1: Option to remove only sideband photons from jet list. Not complete yet.
 * When I do this, code will fail as it require a tight photon in the event to make
 * the final plots. I got to update the code to give the lead photon info
 * separtely.
 * 	2: Partially complete the lead jet plot after removing only the lead
 * photon from the jet list to debug my problems as Ray requested. I need to add
 * the hist and fill it up.
 * BUGFIX: tightId vector which holds the tight photon info was not cleared after
 * the event. Fixed it.
 *
 * Revision 1.58  2009/06/18 22:30:24  samantha
 * MINOR CHANGES: Initializes all the EM objects removal parameters to zero(they
 * are not to be removed).
 *
 * Revision 1.57  2009/06/18 22:04:38  samantha
 * MAJOR CHANGES: 1:Controls what types of EM Objects get removed from the jet list.
 * 2:When adding met to closest jet I look for a jet with Pt>3GeV now.
 * MINOR CHANGES: Few Pt>0.001 statements to avoid the annoying Warning I get from
 * TLorenzVector::Eta() when Pt=0.
 *
 * Revision 1.56  2009/06/16 04:22:52  samantha
 * MAJOR CHANGES: 1: Swicthed back to the original (version 1) of the MEtSig
 * calculation by replacing InitMetProbStuff(). Sasha email me that part again.
 * We are trying to figure out if the MEt model is broken while we try to fix it
 * by making all the changes. Still I do not see a good agreement in MEtSig plot.
 * See elog#1216.
 *
 * Revision 1.55  2009/06/10 23:30:24  samantha
 * ADDED: 	1: Option to decide the method I add MEt to a jet. To the jet closest to
 * MEt or to the most mismeasured jet.
 * MODIFIED1: FindMatchingHEPGPar() scans only the first 30 or so particles. This
 * helps in the DiPhoX test for some reason I forgot I does not find mathches to
 * HEPG partciles.!
 *
 * Revision 1.54  2009/06/01 23:41:50  samantha
 * ADDED:	1: iPrintLvl to control the amount of debug print statements. Higher
 * 	number means more prints.
 * 	2: Few new hists for debugging. Compared the lead det pho and lead det
 * 	jet to hepg photon in the three diffrent regions of the
 * 	DelPhi(Lead Jet, MEt) plot. See elog#1203,1204,1207,1208.
 * 	3: Upadeted the FillAnalysisHistograms() to get the needed info
 * 	to fill the above new hists.
 * 	4: Unpacking HEPG info for photon energy resolution study.
 * 	5: FindMatchingHEPGPar() will find a matching HEPG particle for a
 * 	given 4-vec and particle type and return the 4-vec of the matched.
 * 	6: FindMatchingHADJet() will find the matching HAD jet for a given
 * 	4-vec within a given DeltaR.
 * MODIFIED: 1: In EndJob summary, changed the njet cut info to reflect
 * 	'no cut on jets'.
 * 	2: Formatted the output of the DumpEMobjects()
 * DELETED:1: Previoud debugging hists. Profile hists of Pho/Jet ratio etc.
 * 	2: Print statements in the EndJob about JetEt,Njet hists.
 *
 * Revision 1.53  2009/05/25 18:31:14  samantha
 * ADDED: Introduced a hard MetCut to see if it will fix the DelPhi(Jet,MEt) plot.
 * I placed the cut after MetSig cut. MetAll and MetSig is filled before this cut.
 *
 * Revision 1.52  2009/05/21 22:15:24  samantha
 * DELETE: All the hDblSmear hists.
 *
 * Revision 1.51  2009/05/21 18:50:45  samantha
 * MAJOR CHANGES: Replaced TRandom with TRandom3 as per Nils suggestion in Double
 * smearing Jets. MEt addition is based on the most mismeasured jet. DELETED: All
 * the ratio profile hists that was lying around.
 *
 * Revision 1.50  2009/05/13 22:25:05  samantha
 * MAJOR CHNAGES: 1. Added stuff to smear the jets for pseudo experiments twice.
 * 2. Modified things to use the same hist fill routine FillAnalysisHistograms()
 * for both DATA and pseudo experiments. 3. MyHtAll is modified to be used in both
 * DATA and BG.
 * ADDED:  1. Switch to enable/disbale additional smearing and function Smear_newJetL6FirstTime().
 * 2. Added 2nd Jet Et hist and Eta hists for two lead jets.  3. GetNjets() returns number of
 * jet above a given Pt threshold for given set of jets (a vector with 4-momenta) 4. Expanded the
 * DumpJets() to dump more info about other jets.
 * DELETED: 1. MyMaddCut() which is not used.  2. FillBackgroundAnalysisHistograms()
 * which is now reaplced by FillAnalysisHistograms(). 3. Some of the debug hists and stuff in
 * FillAnalysisHistograms() and DoMyAnalysis()  4. SmearedHtAll() is deleted and MyHtAll
 * is modified to use in all occasions.
 *
 * Revision 1.49  2009/05/07 02:16:27  samantha
 * MAJOR CHANGES: Option to use HAD jets for pseudo experiments or detector jets. MEt will be added to the jet closest to the MEt when using detector jets for pseudo experiments. Checked the effect of artificial 'off_set' in JER peak position when generating pseudo experiments (see elog #1168). Fixed the bug in EM obj Et hist (I was requiring >=2 photons when filling this) MINOR CHANGES: lot of spacing and comments adjustments. Also trying align operators for easy read like 'if (blah && blah || blah..) etc.
 *
 * Revision 1.48  2009/04/17 19:14:33  samantha
 * MAJOR CHANGES: This version I am using as many HADron jets as possible for the pseudo-experiments. See Elog#1155 for more description. Identified detector EM Objects are macthed to HAD jets and removed from HADjet list. Remaining HAD jets are stuffed in for pseudo experiments jet vectors. I do not add MEt to jets though! DELETED:  ReplaceLev6JetWithHadJets(). CHANGED: 1] MyNewJetWithMet() to pass in EM object info.  2] Setting the hDblSmear.hDbl_MetAddedJetMisMeasureProb_ByJetEt hists X-axis title. 3] Commented out the wrongly filled histograms in final hist filling routines, FillAnalysisHistograms() and FillBackgroundHistograms().  4] Lots of print statements in the JetMismeasureProb() to figure out the problem in mismeasure calculation (Elog#1152). 5] Modified the way of summarizing Min/Max NJets requirement in the endjob summary. Now it will show only njets>=1 cases required from a particular Et.   ADDED: Dump_newJetLev6() to dump MetAdded jets.
 *
 * Revision 1.47  2009/04/14 20:48:29  samantha
 * ADDED: Lot of print statements to JetMismeasureProb() to study the problem in MEtAddition. See Elog#1144 attachment 1, the very first event. This shows that MEt should have been added to jet 0 not jet 1 as it is closer to MEt. So this bug has to be fixed. I did a test without MET addition for inclusive njet15 (Elog#1150) but I forgot to fix the final scaling of histograms correctly. I should use the hist->Integral() instead of hist->GetEntries() CHANGED: Removed all the compiler warnings about implicit type conversion, like unsigned to int etc. DELETED: some hists that are not used like the generated JerRandom point distribution (complete list - hDbl_JerRandom[6], hDbl_JerRandom_JetE[6], hDbl_JerRandom_JetDetEta[6], smearfacVsGenJer,smearfacVsJetE,smearfacVsGenJerRatio,smearfac_withOffset,offsets,JerRandom.)
 *
 * Revision 1.46  2009/04/10 21:21:34  samantha
 * ADDED: hDbl_MetAddedJetMisMeasureProb_ByJetEt to plots the mismeasure probability of the MEt added jet by Jet Et bins. See Elog#1146. Also dupmping information of MEt addition to low Pt jets.  see Elog#1144. I also added the jetstuff.smeared_Njet_th0,5... initialization to the standard clearing/reset to function.  RENAMED: MostMisMeasuredJet() to JetMismeasureProb(). OTHER: Commented out the wrong/unwatnted hists from FillAnalysisHitograms() and FillBackgroundHistograms().
 *
 * Revision 1.45  2009/04/07 17:52:54  samantha
 * MAJOR CHANGES: Now I am adding MEt to the least mismeasured jet. Results are in elog#1125.  ADDED: :MostMisMeasuredJet() to measure the mismeasure probability of a jet. Information is stored for later used under JetStuff. There are statements to dump the MEt addition.  DELETED: Removed the unwanted Sasha's header files and SortTowersByEnergy() function. SortTowersByEnergy() is under utils/FreeFunctions.cc .
 *
 *
 */


#include <iomanip>
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Stntuple/loop/TStnAna.hh"
#include "samantha/Pho/JetFilterModuleV2.hh"
#include "samantha/utils/FreeFunctions.hh"
#include <set>
#include <cmath>
#include "TRandom3.h"
#include <assert.h>


std::vector< std::vector<double> > JetFilterModuleV2::vJerParam;
std::vector< std::vector<double> > JetFilterModuleV2::vJerParamErr;
std::vector< std::vector<double> > JetFilterModuleV2::vJerSystParam;
std::vector< std::vector<double> > JetFilterModuleV2::vJerSystParamErr;

//std::vector<std::vector<std::vector<double> > > vJerBinValue;

//------- Function for MET resolution due to unclustered energy
double UnclMetResolution(double* x, double* par){
	double value=0.0;
	double mean=par[0];
	double sigma1=par[1];
	double sigma2=par[3]*sigma1;
	double arg=x[0];
	//   value=(sigma1*TMath::Gaus(arg,mean,sigma1)+par[2]*sigma2*TMath::Gaus(arg,mean,sigma2))/(sigma1+sigma2*par[2]);
	value=(TMath::Gaus(arg,mean,sigma1)+par[2]*TMath::Gaus(arg,mean,sigma2))/(1.0+par[2]);
	return value;
}

//-------- Jet Energy Resolution: centered at zero
// //-------------- fit function: [Const*Gaus(-x/(1+x))+Landau(-x/(1+x))]/(1+Const)
double MyJER(double* x, double* par) {
	double value = 0.0;
	double arg    = -x[0]/(x[0] + 1.0);
	double arg_L  = -x[0]/(x[0] + 1.0);
	double mean   = par[0];
	double sigmaG = par[1];
	double mpv    = par[2];
	double sigmaL = par[3];
	double normL  = par[4];
	if (normL < 0.0) normL = 0.0;
	double f1 = normL / (1.0 + normL) * TMath::Gaus(arg, mean, sigmaG);
	// I had to do this trick to avoid ROOT getting into infinite
	// loop while attempting evaluate f2 for Landau with a very very narrow
	// width! See elog#1270,1271,1273
	// So in cases when Landau width is negligible I avoid evaluating it
	// and set its contribution to zero.
	// 0.001 is just some values I picked. No particular reason.
	double f2 = 0;
	if (sigmaL>0.001) f2 = TMath::Landau(arg_L, mpv, sigmaL) / (1.0 + normL);

	value = f1 + f2;
	if (value<0.0) value = 0.0;

//	std::cout << setw(10) << "Index" << setw(10) << "Pt" << setw(10) << "Eta" << setw(10) << "Phi" << std::endl;
	
//	std::cout << "$$$"<< __FUNCTION__ << ":: f1, f2, val = " << f1 << ", " 
//				<< f2 << ", " << value << std::endl;
	return value;
}


//________________ returns correlation coefficients for JER params
double JERCorrCoeff(int i, int j)
{
	double corr_param[15][11]; // [eta][param]
	
	corr_param[0][0]=-0.953359; // Gauss mean: c01
	corr_param[0][1]=-0.883424; // Gauss mean: c02
	corr_param[0][2]=0.783976; // Gauss mean: c12
	corr_param[0][3]=-0.769484; // Gauss sigma: c01
	corr_param[0][4]=-0.916005; // Landau mpv: c01
	corr_param[0][5]=-0.920646; // Landau mpv: c02
	corr_param[0][6]=0.774913; // Landau mpv: c12
	corr_param[0][7]=-0.798127; // Landau sigma: c01
	corr_param[0][8]=-0.979048; // norm: c01
	corr_param[0][9]=0.914083; // norm: c02
	corr_param[0][10]=-0.972488; // norm: c12
	corr_param[1][0]=-0.930448; // Gauss mean: c01
	corr_param[1][1]=-0.893613; // Gauss mean: c02
	corr_param[1][2]=0.743464; // Gauss mean: c12
	corr_param[1][3]=-0.855253; // Gauss sigma: c01
	corr_param[1][4]=-0.93506; // Landau mpv: c01
	corr_param[1][5]=-0.790079; // Landau mpv: c02
	corr_param[1][6]=0.633928; // Landau mpv: c12
	corr_param[1][7]=-0.791488; // Landau sigma: c01
	corr_param[1][8]=-0.979059; // norm: c01
	corr_param[1][9]=0.935759; // norm: c02
	corr_param[1][10]=-0.978165; // norm: c12
	corr_param[2][0]=-0.933467; // Gauss mean: c01
	corr_param[2][1]=-0.903664; // Gauss mean: c02
	corr_param[2][2]=0.769636; // Gauss mean: c12
	corr_param[2][3]=-0.834006; // Gauss sigma: c01
	corr_param[2][4]=-0.941456; // Landau mpv: c01
	corr_param[2][5]=-0.863316; // Landau mpv: c02
	corr_param[2][6]=0.730447; // Landau mpv: c12
	corr_param[2][7]=-0.795485; // Landau sigma: c01
	corr_param[2][8]=-0.973561; // norm: c01
	corr_param[2][9]=0.922559; // norm: c02
	corr_param[2][10]=-0.979269; // norm: c12
	corr_param[3][0]=-0.92146; // Gauss mean: c01
	corr_param[3][1]=-0.875047; // Gauss mean: c02
	corr_param[3][2]=0.720247; // Gauss mean: c12
	corr_param[3][3]=-0.820492; // Gauss sigma: c01
	corr_param[3][4]=-0.932934; // Landau mpv: c01
	corr_param[3][5]=-0.806167; // Landau mpv: c02
	corr_param[3][6]=0.666462; // Landau mpv: c12
	corr_param[3][7]=-0.752656; // Landau sigma: c01
	corr_param[3][8]=-0.97604; // norm: c01
	corr_param[3][9]=0.916931; // norm: c02
	corr_param[3][10]=-0.970166; // norm: c12
	corr_param[4][0]=-0.921373; // Gauss mean: c01
	corr_param[4][1]=-0.865997; // Gauss mean: c02
	corr_param[4][2]=0.713135; // Gauss mean: c12
	corr_param[4][3]=-0.812212; // Gauss sigma: c01
	corr_param[4][4]=-0.933552; // Landau mpv: c01
	corr_param[4][5]=-0.826206; // Landau mpv: c02
	corr_param[4][6]=0.696959; // Landau mpv: c12
	corr_param[4][7]=-0.731748; // Landau sigma: c01
	corr_param[4][8]=-0.972365; // norm: c01
	corr_param[4][9]=0.900603; // norm: c02
	corr_param[4][10]=-0.96338; // norm: c12
	corr_param[5][0]=-0.945056; // Gauss mean: c01
	corr_param[5][1]=-0.918674; // Gauss mean: c02
	corr_param[5][2]=0.808282; // Gauss mean: c12
	corr_param[5][3]=-0.811252; // Gauss sigma: c01
	corr_param[5][4]=-0.941345; // Landau mpv: c01
	corr_param[5][5]=-0.908828; // Landau mpv: c02
	corr_param[5][6]=0.805962; // Landau mpv: c12
	corr_param[5][7]=-0.737379; // Landau sigma: c01
	corr_param[5][8]=-0.974342; // norm: c01
	corr_param[5][9]=0.878683; // norm: c02
	corr_param[5][10]=-0.952568; // norm: c12
	corr_param[6][0]=-0.863692; // Gauss mean: c01
	corr_param[6][1]=-0.898828; // Gauss mean: c02
	corr_param[6][2]=0.683443; // Gauss mean: c12
	corr_param[6][3]=-0.864674; // Gauss sigma: c01
	corr_param[6][4]=-0.910948; // Landau mpv: c01
	corr_param[6][5]=-0.878636; // Landau mpv: c02
	corr_param[6][6]=0.709797; // Landau mpv: c12
	corr_param[6][7]=-0.790837; // Landau sigma: c01
	corr_param[6][8]=-0.976028; // norm: c01
	corr_param[6][9]=0.920123; // norm: c02
	corr_param[6][10]=-0.976345; // norm: c12
	corr_param[7][0]=-0.901329; // Gauss mean: c01
	corr_param[7][1]=-0.868814; // Gauss mean: c02
	corr_param[7][2]=0.664341; // Gauss mean: c12
	corr_param[7][3]=-0.879818; // Gauss sigma: c01
	corr_param[7][4]=-0.910909; // Landau mpv: c01
	corr_param[7][5]=-0.848758; // Landau mpv: c02
	corr_param[7][6]=0.654596; // Landau mpv: c12
	corr_param[7][7]=-0.854805; // Landau sigma: c01
	corr_param[7][8]=-0.98475; // norm: c01
	corr_param[7][9]=0.944593; // norm: c02
	corr_param[7][10]=-0.983111; // norm: c12
	corr_param[8][0]=-0.93546; // Gauss mean: c01
	corr_param[8][1]=-0.919387; // Gauss mean: c02
	corr_param[8][2]=0.775089; // Gauss mean: c12
	corr_param[8][3]=-0.885525; // Gauss sigma: c01
	corr_param[8][4]=-0.902932; // Landau mpv: c01
	corr_param[8][5]=-0.858334; // Landau mpv: c02
	corr_param[8][6]=0.664962; // Landau mpv: c12
	corr_param[8][7]=-0.83488; // Landau sigma: c01
	corr_param[8][8]=-0.989491; // norm: c01
	corr_param[8][9]=0.946333; // norm: c02
	corr_param[8][10]=-0.980041; // norm: c12
	corr_param[9][0]=-0.958622; // Gauss mean: c01
	corr_param[9][1]=-0.92783; // Gauss mean: c02
	corr_param[9][2]=0.813899; // Gauss mean: c12
	corr_param[9][3]=-0.918357; // Gauss sigma: c01
	corr_param[9][4]=-0.946212; // Landau mpv: c01
	corr_param[9][5]=-0.874383; // Landau mpv: c02
	corr_param[9][6]=0.732471; // Landau mpv: c12
	corr_param[9][7]=-0.890866; // Landau sigma: c01
	corr_param[9][8]=-0.993186; // norm: c01
	corr_param[9][9]=0.965285; // norm: c02
	corr_param[9][10]=-0.987298; // norm: c12
	corr_param[10][0]=-0.959116; // Gauss mean: c01
	corr_param[10][1]=-0.934235; // Gauss mean: c02
	corr_param[10][2]=0.822224; // Gauss mean: c12
	corr_param[10][3]=-0.929578; // Gauss sigma: c01
	corr_param[10][4]=-0.945125; // Landau mpv: c01
	corr_param[10][5]=-0.889286; // Landau mpv: c02
	corr_param[10][6]=0.741957; // Landau mpv: c12
	corr_param[10][7]=-0.900695; // Landau sigma: c01
	corr_param[10][8]=-0.992615; // norm: c01
	corr_param[10][9]=0.967849; // norm: c02
	corr_param[10][10]=-0.989424; // norm: c12
	corr_param[11][0]=-0.961204; // Gauss mean: c01
	corr_param[11][1]=-0.933701; // Gauss mean: c02
	corr_param[11][2]=0.822536; // Gauss mean: c12
	corr_param[11][3]=-0.938163; // Gauss sigma: c01
	corr_param[11][4]=-0.95; // Landau mpv: c01
	corr_param[11][5]=-0.875959; // Landau mpv: c02
	corr_param[11][6]=0.737455; // Landau mpv: c12
	corr_param[11][7]=-0.895806; // Landau sigma: c01
	corr_param[11][8]=-0.993612; // norm: c01
	corr_param[11][9]=0.968624; // norm: c02
	corr_param[11][10]=-0.988782; // norm: c12
	corr_param[12][0]=-0.961494; // Gauss mean: c01
	corr_param[12][1]=-0.888003; // Gauss mean: c02
	corr_param[12][2]=0.756614; // Gauss mean: c12
	corr_param[12][3]=-0.940026; // Gauss sigma: c01
	corr_param[12][4]=-0.971454; // Landau mpv: c01
	corr_param[12][5]=-0.907369; // Landau mpv: c02
	corr_param[12][6]=0.819359; // Landau mpv: c12
	corr_param[12][7]=-0.890781; // Landau sigma: c01
	corr_param[12][8]=-0.995606; // norm: c01
	corr_param[12][9]=0.976755; // norm: c02
	corr_param[12][10]=-0.991439; // norm: c12
	corr_param[13][0]=-0.979353; // Gauss mean: c01
	corr_param[13][1]=-0.901815; // Gauss mean: c02
	corr_param[13][2]=0.819153; // Gauss mean: c12
	corr_param[13][3]=-0.906452; // Gauss sigma: c01
	corr_param[13][4]=-0.977947; // Landau mpv: c01
	corr_param[13][5]=-0.902621; // Landau mpv: c02
	corr_param[13][6]=0.821221; // Landau mpv: c12
	corr_param[13][7]=-0.907155; // Landau sigma: c01
	corr_param[13][8]=-0.996171; // norm: c01
	corr_param[13][9]=0.98616; // norm: c02
	corr_param[13][10]=-0.996056; // norm: c12
	corr_param[14][0]=0; // Gauss mean: c01
	corr_param[14][1]=0; // Gauss mean: c02
	corr_param[14][2]=0; // Gauss mean: c12
	corr_param[14][3]=0; // Gauss sigma: c01
	corr_param[14][4]=-0.993129; // Landau mpv: c01
	corr_param[14][5]=-0.979168; // Landau mpv: c02
	corr_param[14][6]=0.957156; // Landau mpv: c12
	corr_param[14][7]=-0.895607; // Landau sigma: c01
	corr_param[14][8]=0; // norm: c01
	corr_param[14][9]=0; // norm: c02
	corr_param[14][10]=0; // norm: c12


	double value=0.0;
	if(i<15 && j<11) value=corr_param[i][j];

	return value;
}

//________________ returns 2nd NEW JER param or its uncertainty
double new2JERparam(int code, int i, int j)
{
	//these are derived from stiching old sigma function.
	//split the fit range to 2 and fit either side with
	//the sigma function and merge them via 2 step functions.
	double new2jer_param[15][25]; // [eta][param]
	double new2jer_param_er[15][25]; // [eta][param]
	double value=0.0;

	// NEW parameters 07-28-2009
	new2jer_param[0][0]=19.8678;  new2jer_param_er[0][0]=4.66326;
	new2jer_param[0][1]=-4.49957e-05;  new2jer_param_er[0][1]=4.59938e-06;
	new2jer_param[0][2]=0.153883;  new2jer_param_er[0][2]=0.360298;
	new2jer_param[0][3]=35.1405;  new2jer_param_er[0][3]=1.04008;
	new2jer_param[0][4]=-271.177;  new2jer_param_er[0][4]=8.25098;
	new2jer_param[0][5]=-0.000330151;  new2jer_param_er[0][5]=0.000661854;
	new2jer_param[0][6]=5.91603e-05;  new2jer_param_er[0][6]=3.19889e-06;
	new2jer_param[0][7]=12.3922;  new2jer_param_er[0][7]=0.246686;
	new2jer_param[0][8]=29.8272;  new2jer_param_er[0][8]=4.64062;
	new2jer_param[0][9]=37.8957;  new2jer_param_er[0][9]=4.36968;
	new2jer_param[0][10]=0.00147467;  new2jer_param_er[0][10]=0.000896253;
	new2jer_param[0][11]=0.000225551;  new2jer_param_er[0][11]=1.86132e-05;
	new2jer_param[0][12]=-1.99797;  new2jer_param_er[0][12]=1.1247;
	new2jer_param[0][13]=46.7961;  new2jer_param_er[0][13]=5.89435;
	new2jer_param[0][14]=51.9756;  new2jer_param_er[0][14]=49.5571;
	new2jer_param[0][15]=1.94363;  new2jer_param_er[0][15]=0.228898;
	new2jer_param[0][16]=0.00021461;  new2jer_param_er[0][16]=7.31778e-06;
	new2jer_param[0][17]=2.13276;  new2jer_param_er[0][17]=0.0496794;
	new2jer_param[0][18]=90.3924;  new2jer_param_er[0][18]=3.23837;
	new2jer_param[0][19]=-639.947;  new2jer_param_er[0][19]=10.6542;
	new2jer_param[0][20]=3.3013;  new2jer_param_er[0][20]=0.912313;
	new2jer_param[0][21]=0.000536978;  new2jer_param_er[0][21]=2.50028e-05;
	new2jer_param[0][22]=-2.05617;  new2jer_param_er[0][22]=0.429164;
	new2jer_param[0][23]=31.908;  new2jer_param_er[0][23]=1.6553;
	new2jer_param[0][24]=-174.311;  new2jer_param_er[0][24]=8.83889;
	new2jer_param[1][0]=1.82786e-05;  new2jer_param_er[1][0]=8.44509e-06;
	new2jer_param[1][1]=112.775;  new2jer_param_er[1][1]=36.1985;
	new2jer_param[1][2]=103535;  new2jer_param_er[1][2]=31076.9;
	new2jer_param[1][3]=-93.83;  new2jer_param_er[1][3]=3.90891;
	new2jer_param[1][4]=-1027.05;  new2jer_param_er[1][4]=30.4902;
	new2jer_param[1][5]=0.000145558;  new2jer_param_er[1][5]=4.2197e-06;
	new2jer_param[1][6]=24.2616;  new2jer_param_er[1][6]=1.2647;
	new2jer_param[1][7]=40984;  new2jer_param_er[1][7]=2127.99;
	new2jer_param[1][8]=-88.6887;  new2jer_param_er[1][8]=1.72304;
	new2jer_param[1][9]=-820.441;  new2jer_param_er[1][9]=4.45765;
	new2jer_param[1][10]=0.000133326;  new2jer_param_er[1][10]=2.19806e-05;
	new2jer_param[1][11]=0.000247177;  new2jer_param_er[1][11]=0.000170501;
	new2jer_param[1][12]=-73.1561;  new2jer_param_er[1][12]=48.0352;
	new2jer_param[1][13]=102.206;  new2jer_param_er[1][13]=7.95808;
	new2jer_param[1][14]=439.84;  new2jer_param_er[1][14]=71.728;
	new2jer_param[1][15]=-8.66329e-06;  new2jer_param_er[1][15]=9.66547e-05;
	new2jer_param[1][16]=7.98535e-05;  new2jer_param_er[1][16]=2.49972e-06;
	new2jer_param[1][17]=7.88279;  new2jer_param_er[1][17]=0.220199;
	new2jer_param[1][18]=45.2114;  new2jer_param_er[1][18]=3.95301;
	new2jer_param[1][19]=68.3827;  new2jer_param_er[1][19]=5.57978;
	new2jer_param[1][20]=4.19301;  new2jer_param_er[1][20]=0.841744;
	new2jer_param[1][21]=0.00018651;  new2jer_param_er[1][21]=1.49985e-05;
	new2jer_param[1][22]=-0.348334;  new2jer_param_er[1][22]=0.147113;
	new2jer_param[1][23]=67.2871;  new2jer_param_er[1][23]=1.29538;
	new2jer_param[1][24]=-406.969;  new2jer_param_er[1][24]=13.5397;
	new2jer_param[2][0]=-2.63618e-07;  new2jer_param_er[2][0]=6.25021e-06;
	new2jer_param[2][1]=37.2498;  new2jer_param_er[2][1]=19.9259;
	new2jer_param[2][2]=8644.6;  new2jer_param_er[2][2]=4212.62;
	new2jer_param[2][3]=-82.4164;  new2jer_param_er[2][3]=2.53039;
	new2jer_param[2][4]=-795.108;  new2jer_param_er[2][4]=34.6097;
	new2jer_param[2][5]=0.000136734;  new2jer_param_er[2][5]=4.10764e-06;
	new2jer_param[2][6]=28.7246;  new2jer_param_er[2][6]=5.67768;
	new2jer_param[2][7]=47787;  new2jer_param_er[2][7]=9197.38;
	new2jer_param[2][8]=-88.4299;  new2jer_param_er[2][8]=1.6018;
	new2jer_param[2][9]=-830.843;  new2jer_param_er[2][9]=16.1476;
	new2jer_param[2][10]=-4.50036e-05;  new2jer_param_er[2][10]=6.48786e-05;
	new2jer_param[2][11]=0.00662261;  new2jer_param_er[2][11]=0.00297834;
	new2jer_param[2][12]=-168.688;  new2jer_param_er[2][12]=74.1123;
	new2jer_param[2][13]=-392.793;  new2jer_param_er[2][13]=41.2839;
	new2jer_param[2][14]=-1804.87;  new2jer_param_er[2][14]=196.051;
	new2jer_param[2][15]=9.50783e-05;  new2jer_param_er[2][15]=3.14005e-06;
	new2jer_param[2][16]=13.32;  new2jer_param_er[2][16]=2.09324;
	new2jer_param[2][17]=30657.7;  new2jer_param_er[2][17]=4485.86;
	new2jer_param[2][18]=-114.822;  new2jer_param_er[2][18]=3.32755;
	new2jer_param[2][19]=-1091.61;  new2jer_param_er[2][19]=20.3077;
	new2jer_param[2][20]=0.000722077;  new2jer_param_er[2][20]=0.00120837;
	new2jer_param[2][21]=7.33606e-05;  new2jer_param_er[2][21]=1.59523e-05;
	new2jer_param[2][22]=23.827;  new2jer_param_er[2][22]=0.893251;
	new2jer_param[2][23]=17.9333;  new2jer_param_er[2][23]=2.16877;
	new2jer_param[2][24]=54.0246;  new2jer_param_er[2][24]=6.41162;
	new2jer_param[3][0]=2.6076e-05;  new2jer_param_er[3][0]=1.00655e-05;
	new2jer_param[3][1]=5.62049;  new2jer_param_er[3][1]=0.743139;
	new2jer_param[3][2]=5079.88;  new2jer_param_er[3][2]=710.214;
	new2jer_param[3][3]=-96.9768;  new2jer_param_er[3][3]=4.19675;
	new2jer_param[3][4]=-768.2;  new2jer_param_er[3][4]=12.295;
	new2jer_param[3][5]=0.000160302;  new2jer_param_er[3][5]=5.47767e-06;
	new2jer_param[3][6]=36.6061;  new2jer_param_er[3][6]=6.21387;
	new2jer_param[3][7]=56614.5;  new2jer_param_er[3][7]=8842.1;
	new2jer_param[3][8]=-85.6586;  new2jer_param_er[3][8]=1.90631;
	new2jer_param[3][9]=-818.654;  new2jer_param_er[3][9]=14.0542;
	new2jer_param[3][10]=102695;  new2jer_param_er[3][10]=28701.5;
	new2jer_param[3][11]=-3.21476e-05;  new2jer_param_er[3][11]=1.36107e-05;
	new2jer_param[3][12]=-1.78761;  new2jer_param_er[3][12]=0.0608095;
	new2jer_param[3][13]=184.409;  new2jer_param_er[3][13]=2.52695;
	new2jer_param[3][14]=-3607.21;  new2jer_param_er[3][14]=51.601;
	new2jer_param[3][15]=9.15004e-06;  new2jer_param_er[3][15]=8.04309e-05;
	new2jer_param[3][16]=7.17268e-05;  new2jer_param_er[3][16]=2.6516e-06;
	new2jer_param[3][17]=9.42404;  new2jer_param_er[3][17]=0.217615;
	new2jer_param[3][18]=41.5946;  new2jer_param_er[3][18]=2.44982;
	new2jer_param[3][19]=71.7862;  new2jer_param_er[3][19]=4.03082;
	new2jer_param[3][20]=5.19853;  new2jer_param_er[3][20]=1.14223;
	new2jer_param[3][21]=0.000164887;  new2jer_param_er[3][21]=1.652e-05;
	new2jer_param[3][22]=0.140838;  new2jer_param_er[3][22]=0.14477;
	new2jer_param[3][23]=71.745;  new2jer_param_er[3][23]=2.39114;
	new2jer_param[3][24]=-469.562;  new2jer_param_er[3][24]=15.7655;
	new2jer_param[4][0]=47.5497;  new2jer_param_er[4][0]=7.02596;
	new2jer_param[4][1]=0.000156976;  new2jer_param_er[4][1]=7.72677e-06;
	new2jer_param[4][2]=0.87687;  new2jer_param_er[4][2]=0.13449;
	new2jer_param[4][3]=73.551;  new2jer_param_er[4][3]=2.54369;
	new2jer_param[4][4]=-709.306;  new2jer_param_er[4][4]=10.8735;
	new2jer_param[4][5]=-0.000161065;  new2jer_param_er[4][5]=0.000257453;
	new2jer_param[4][6]=0.00017144;  new2jer_param_er[4][6]=4.38609e-06;
	new2jer_param[4][7]=12.6641;  new2jer_param_er[4][7]=0.254432;
	new2jer_param[4][8]=28.5781;  new2jer_param_er[4][8]=1.79537;
	new2jer_param[4][9]=42.3101;  new2jer_param_er[4][9]=2.35822;
	new2jer_param[4][10]=0.17238;  new2jer_param_er[4][10]=0.0419682;
	new2jer_param[4][11]=-5.02959e-05;  new2jer_param_er[4][11]=2.47804e-05;
	new2jer_param[4][12]=-1.58198;  new2jer_param_er[4][12]=0.0653919;
	new2jer_param[4][13]=211.993;  new2jer_param_er[4][13]=42.7189;
	new2jer_param[4][14]=-1272.34;  new2jer_param_er[4][14]=463.252;
	new2jer_param[4][15]=4.14315e-05;  new2jer_param_er[4][15]=8.95393e-05;
	new2jer_param[4][16]=8.26805e-05;  new2jer_param_er[4][16]=3.53552e-06;
	new2jer_param[4][17]=11.439;  new2jer_param_er[4][17]=0.323368;
	new2jer_param[4][18]=42.9396;  new2jer_param_er[4][18]=2.6298;
	new2jer_param[4][19]=77.5841;  new2jer_param_er[4][19]=4.61069;
	new2jer_param[4][20]=0.00115133;  new2jer_param_er[4][20]=0.00136645;
	new2jer_param[4][21]=7.65547e-05;  new2jer_param_er[4][21]=1.78474e-05;
	new2jer_param[4][22]=28.0857;  new2jer_param_er[4][22]=1.26509;
	new2jer_param[4][23]=18.7907;  new2jer_param_er[4][23]=2.85836;
	new2jer_param[4][24]=68.1781;  new2jer_param_er[4][24]=11.5539;
	new2jer_param[5][0]=0.000141175;  new2jer_param_er[5][0]=0.000208054;
	new2jer_param[5][1]=6.8018e-05;  new2jer_param_er[5][1]=6.47675e-06;
	new2jer_param[5][2]=12.2424;  new2jer_param_er[5][2]=0.66104;
	new2jer_param[5][3]=41.7701;  new2jer_param_er[5][3]=6.1064;
	new2jer_param[5][4]=91.7719;  new2jer_param_er[5][4]=11.738;
	new2jer_param[5][5]=-0.000167167;  new2jer_param_er[5][5]=0.000159012;
	new2jer_param[5][6]=0.000126344;  new2jer_param_er[5][6]=5.1593e-06;
	new2jer_param[5][7]=17.0935;  new2jer_param_er[5][7]=0.569867;
	new2jer_param[5][8]=49.6915;  new2jer_param_er[5][8]=3.73185;
	new2jer_param[5][9]=65.245;  new2jer_param_er[5][9]=5.05287;
	new2jer_param[5][10]=0.00047374;  new2jer_param_er[5][10]=5.62231e-05;
	new2jer_param[5][11]=6.4775e-05;  new2jer_param_er[5][11]=4.93433e-05;
	new2jer_param[5][12]=-30.6499;  new2jer_param_er[5][12]=14.3355;
	new2jer_param[5][13]=97.842;  new2jer_param_er[5][13]=5.39819;
	new2jer_param[5][14]=237.931;  new2jer_param_er[5][14]=46.9107;
	new2jer_param[5][15]=-6.04337e-05;  new2jer_param_er[5][15]=0.000236725;
	new2jer_param[5][16]=6.19245e-05;  new2jer_param_er[5][16]=5.1413e-06;
	new2jer_param[5][17]=12.3866;  new2jer_param_er[5][17]=0.647418;
	new2jer_param[5][18]=50.6366;  new2jer_param_er[5][18]=9.65717;
	new2jer_param[5][19]=90.581;  new2jer_param_er[5][19]=16.3769;
	new2jer_param[5][20]=0.00614505;  new2jer_param_er[5][20]=0.00112699;
	new2jer_param[5][21]=0.000360885;  new2jer_param_er[5][21]=5.19201e-05;
	new2jer_param[5][22]=1.91428;  new2jer_param_er[5][22]=2.35163;
	new2jer_param[5][23]=76.9499;  new2jer_param_er[5][23]=9.95366;
	new2jer_param[5][24]=117.18;  new2jer_param_er[5][24]=27.9691;
	new2jer_param[6][0]=-0.000101154;  new2jer_param_er[6][0]=1.66316e-05;
	new2jer_param[6][1]=0.000806924;  new2jer_param_er[6][1]=0.000129755;
	new2jer_param[6][2]=0.919908;  new2jer_param_er[6][2]=0.136703;
	new2jer_param[6][3]=-70.5321;  new2jer_param_er[6][3]=10.7575;
	new2jer_param[6][4]=184.316;  new2jer_param_er[6][4]=21.6584;
	new2jer_param[6][5]=9.87787e-05;  new2jer_param_er[6][5]=8.31655e-06;
	new2jer_param[6][6]=5.2306;  new2jer_param_er[6][6]=0.348162;
	new2jer_param[6][7]=8317.97;  new2jer_param_er[6][7]=547.541;
	new2jer_param[6][8]=-104.162;  new2jer_param_er[6][8]=3.42991;
	new2jer_param[6][9]=-795.571;  new2jer_param_er[6][9]=6.35185;
	new2jer_param[6][10]=1.96618e-05;  new2jer_param_er[6][10]=4.99713e-06;
	new2jer_param[6][11]=-0.00503536;  new2jer_param_er[6][11]=0.00788996;
	new2jer_param[6][12]=-4.45447;  new2jer_param_er[6][12]=3.41862;
	new2jer_param[6][13]=-28.2124;  new2jer_param_er[6][13]=3.28918;
	new2jer_param[6][14]=14.2267;  new2jer_param_er[6][14]=52.3373;
	new2jer_param[6][15]=0.00018074;  new2jer_param_er[6][15]=4.86509e-06;
	new2jer_param[6][16]=0.0537848;  new2jer_param_er[6][16]=0.027314;
	new2jer_param[6][17]=73.3283;  new2jer_param_er[6][17]=37.5512;
	new2jer_param[6][18]=-87.5383;  new2jer_param_er[6][18]=2.73876;
	new2jer_param[6][19]=-316.586;  new2jer_param_er[6][19]=47.93;
	new2jer_param[6][20]=0.00395055;  new2jer_param_er[6][20]=0.000711368;
	new2jer_param[6][21]=8.67745e-05;  new2jer_param_er[6][21]=4.37582e-05;
	new2jer_param[6][22]=9.32493;  new2jer_param_er[6][22]=3.66678;
	new2jer_param[6][23]=28.3516;  new2jer_param_er[6][23]=4.36477;
	new2jer_param[6][24]=75.6393;  new2jer_param_er[6][24]=4.54423;
	new2jer_param[7][0]=-5.34974e-05;  new2jer_param_er[7][0]=1.48981e-05;
	new2jer_param[7][1]=0.00268381;  new2jer_param_er[7][1]=0.00346928;
	new2jer_param[7][2]=-0.268993;  new2jer_param_er[7][2]=0.616655;
	new2jer_param[7][3]=-90.1315;  new2jer_param_er[7][3]=18.0139;
	new2jer_param[7][4]=-52.5316;  new2jer_param_er[7][4]=162.211;
	new2jer_param[7][5]=-30.1315;  new2jer_param_er[7][5]=5.99022e-05;
	new2jer_param[7][6]=31.6336;  new2jer_param_er[7][6]=6.28897e-05;
	new2jer_param[7][7]=10.1294;  new2jer_param_er[7][7]=0.0891712;
	new2jer_param[7][8]=-1.09792e+07;  new2jer_param_er[7][8]=448.731;
	new2jer_param[7][9]=-533706;  new2jer_param_er[7][9]=21.8262;
	new2jer_param[7][10]=-5.01939e-05;  new2jer_param_er[7][10]=8.09824e-06;
	new2jer_param[7][11]=-0.000471913;  new2jer_param_er[7][11]=9.69877e-05;
	new2jer_param[7][12]=-2.81359;  new2jer_param_er[7][12]=0.099287;
	new2jer_param[7][13]=-19.2502;  new2jer_param_er[7][13]=3.1731;
	new2jer_param[7][14]=106.724;  new2jer_param_er[7][14]=6.08989;
	new2jer_param[7][15]=0.000109567;  new2jer_param_er[7][15]=1.02712e-05;
	new2jer_param[7][16]=0.964028;  new2jer_param_er[7][16]=0.134618;
	new2jer_param[7][17]=1949.65;  new2jer_param_er[7][17]=305.007;
	new2jer_param[7][18]=-87.8154;  new2jer_param_er[7][18]=6.91029;
	new2jer_param[7][19]=-594.478;  new2jer_param_er[7][19]=11.1255;
	new2jer_param[7][20]=-0.00970539;  new2jer_param_er[7][20]=0.0083618;
	new2jer_param[7][21]=0.000292597;  new2jer_param_er[7][21]=7.15134e-05;
	new2jer_param[7][22]=24.6347;  new2jer_param_er[7][22]=1.92959;
	new2jer_param[7][23]=15.9126;  new2jer_param_er[7][23]=1.68159;
	new2jer_param[7][24]=35.6568;  new2jer_param_er[7][24]=9.1284;
	new2jer_param[8][0]=5.81928e-06;  new2jer_param_er[8][0]=9.51715e-06;
	new2jer_param[8][1]=0.000345541;  new2jer_param_er[8][1]=4.18947e-05;
	new2jer_param[8][2]=-1.03657;  new2jer_param_er[8][2]=0.187075;
	new2jer_param[8][3]=-79.0016;  new2jer_param_er[8][3]=11.8852;
	new2jer_param[8][4]=193.348;  new2jer_param_er[8][4]=18.31;
	new2jer_param[8][5]=0.000104963;  new2jer_param_er[8][5]=3.49292e-06;
	new2jer_param[8][6]=23.1058;  new2jer_param_er[8][6]=0.939441;
	new2jer_param[8][7]=55786.7;  new2jer_param_er[8][7]=2694.09;
	new2jer_param[8][8]=-94.9157;  new2jer_param_er[8][8]=1.89123;
	new2jer_param[8][9]=-886.55;  new2jer_param_er[8][9]=4.87381;
	new2jer_param[8][10]=-5.51873;  new2jer_param_er[8][10]=29.9464;
	new2jer_param[8][11]=-3.00207e-05;  new2jer_param_er[8][11]=6.59182e-06;
	new2jer_param[8][12]=-2.53429;  new2jer_param_er[8][12]=0.250205;
	new2jer_param[8][13]=46.6411;  new2jer_param_er[8][13]=3.50745;
	new2jer_param[8][14]=-359.44;  new2jer_param_er[8][14]=256.978;
	new2jer_param[8][15]=-0.00025644;  new2jer_param_er[8][15]=0.00024549;
	new2jer_param[8][16]=1.16327e-05;  new2jer_param_er[8][16]=1.00571e-05;
	new2jer_param[8][17]=7.13343;  new2jer_param_er[8][17]=1.49364;
	new2jer_param[8][18]=57.045;  new2jer_param_er[8][18]=18.3435;
	new2jer_param[8][19]=50.9732;  new2jer_param_er[8][19]=28.241;
	new2jer_param[8][20]=-0.00644595;  new2jer_param_er[8][20]=0.00551788;
	new2jer_param[8][21]=0.000691224;  new2jer_param_er[8][21]=0.000137453;
	new2jer_param[8][22]=53.8604;  new2jer_param_er[8][22]=7.50788;
	new2jer_param[8][23]=21.4632;  new2jer_param_er[8][23]=2.14259;
	new2jer_param[8][24]=62.5323;  new2jer_param_er[8][24]=14.8349;
	new2jer_param[9][0]=-0.000313161;  new2jer_param_er[9][0]=8.7212e-05;
	new2jer_param[9][1]=3.42748e-05;  new2jer_param_er[9][1]=6.85578e-06;
	new2jer_param[9][2]=4.18957;  new2jer_param_er[9][2]=0.33638;
	new2jer_param[9][3]=10.9825;  new2jer_param_er[9][3]=1.33956;
	new2jer_param[9][4]=103.769;  new2jer_param_er[9][4]=3.06172;
	new2jer_param[9][5]=-0.00350841;  new2jer_param_er[9][5]=0.000973922;
	new2jer_param[9][6]=6.37045e-05;  new2jer_param_er[9][6]=7.02804e-06;
	new2jer_param[9][7]=13.6276;  new2jer_param_er[9][7]=0.470483;
	new2jer_param[9][8]=30.3997;  new2jer_param_er[9][8]=1.79715;
	new2jer_param[9][9]=28.4588;  new2jer_param_er[9][9]=4.48326;
	new2jer_param[9][10]=49.0892;  new2jer_param_er[9][10]=85.1444;
	new2jer_param[9][11]=8.90324e-06;  new2jer_param_er[9][11]=8.72579e-06;
	new2jer_param[9][12]=-5.70279;  new2jer_param_er[9][12]=0.252636;
	new2jer_param[9][13]=13.4923;  new2jer_param_er[9][13]=3.00577;
	new2jer_param[9][14]=-112.163;  new2jer_param_er[9][14]=23.4186;
	new2jer_param[9][15]=-0.000332528;  new2jer_param_er[9][15]=0.000426399;
	new2jer_param[9][16]=1.58786e-05;  new2jer_param_er[9][16]=2.04611e-05;
	new2jer_param[9][17]=6.39391;  new2jer_param_er[9][17]=1.91296;
	new2jer_param[9][18]=47.7182;  new2jer_param_er[9][18]=22.525;
	new2jer_param[9][19]=36.3698;  new2jer_param_er[9][19]=30.5127;
	new2jer_param[9][20]=-0.151388;  new2jer_param_er[9][20]=0.316865;
	new2jer_param[9][21]=0.000851038;  new2jer_param_er[9][21]=0.000131773;
	new2jer_param[9][22]=34.6581;  new2jer_param_er[9][22]=3.9375;
	new2jer_param[9][23]=19.6211;  new2jer_param_er[9][23]=1.0444;
	new2jer_param[9][24]=0.11371;  new2jer_param_er[9][24]=39.8918;
	new2jer_param[10][0]=-0.00024438;  new2jer_param_er[10][0]=0.000141817;
	new2jer_param[10][1]=-2.11148e-05;  new2jer_param_er[10][1]=7.59845e-06;
	new2jer_param[10][2]=5.8104;  new2jer_param_er[10][2]=0.509704;
	new2jer_param[10][3]=14.3354;  new2jer_param_er[10][3]=2.70173;
	new2jer_param[10][4]=122.236;  new2jer_param_er[10][4]=7.29174;
	new2jer_param[10][5]=-0.00189466;  new2jer_param_er[10][5]=0.000800394;
	new2jer_param[10][6]=3.8643e-05;  new2jer_param_er[10][6]=8.9635e-06;
	new2jer_param[10][7]=16.5577;  new2jer_param_er[10][7]=0.83526;
	new2jer_param[10][8]=38.3546;  new2jer_param_er[10][8]=3.58156;
	new2jer_param[10][9]=43.8401;  new2jer_param_er[10][9]=8.09963;
	new2jer_param[10][10]=-0.00312013;  new2jer_param_er[10][10]=0.00221981;
	new2jer_param[10][11]=-5.73549e-05;  new2jer_param_er[10][11]=2.09888e-05;
	new2jer_param[10][12]=-4.24355;  new2jer_param_er[10][12]=1.69442;
	new2jer_param[10][13]=33.5526;  new2jer_param_er[10][13]=10.2252;
	new2jer_param[10][14]=17.9554;  new2jer_param_er[10][14]=32.5678;
	new2jer_param[10][15]=-0.00125281;  new2jer_param_er[10][15]=0.00123924;
	new2jer_param[10][16]=3.1922e-05;  new2jer_param_er[10][16]=1.00836e-05;
	new2jer_param[10][17]=6.27554;  new2jer_param_er[10][17]=0.907545;
	new2jer_param[10][18]=38.7817;  new2jer_param_er[10][18]=11.2734;
	new2jer_param[10][19]=10.8881;  new2jer_param_er[10][19]=15.3541;
	new2jer_param[10][20]=0.00238845;  new2jer_param_er[10][20]=0.000105573;
	new2jer_param[10][21]=0.000860129;  new2jer_param_er[10][21]=0.000138294;
	new2jer_param[10][22]=21.1862;  new2jer_param_er[10][22]=6.44095;
	new2jer_param[10][23]=0.335593;  new2jer_param_er[10][23]=4.40729;
	new2jer_param[10][24]=157.256;  new2jer_param_er[10][24]=3.21485;
	new2jer_param[11][0]=-0.000677771;  new2jer_param_er[11][0]=0.000483256;
	new2jer_param[11][1]=1.71819e-06;  new2jer_param_er[11][1]=1.13147e-05;
	new2jer_param[11][2]=2.91113;  new2jer_param_er[11][2]=1.03978;
	new2jer_param[11][3]=24.7033;  new2jer_param_er[11][3]=6.22734;
	new2jer_param[11][4]=119.195;  new2jer_param_er[11][4]=25.4557;
	new2jer_param[11][5]=-0.00107938;  new2jer_param_er[11][5]=0.000614445;
	new2jer_param[11][6]=1.36526e-05;  new2jer_param_er[11][6]=1.80198e-05;
	new2jer_param[11][7]=20.3567;  new2jer_param_er[11][7]=2.55611;
	new2jer_param[11][8]=53.6875;  new2jer_param_er[11][8]=9.02942;
	new2jer_param[11][9]=60.8061;  new2jer_param_er[11][9]=19.107;
	new2jer_param[11][10]=-0.987088;  new2jer_param_er[11][10]=0.556797;
	new2jer_param[11][11]=-6.63354e-05;  new2jer_param_er[11][11]=1.24183e-05;
	new2jer_param[11][12]=-3.25544;  new2jer_param_er[11][12]=0.665516;
	new2jer_param[11][13]=55.2693;  new2jer_param_er[11][13]=4.78842;
	new2jer_param[11][14]=-316.759;  new2jer_param_er[11][14]=31.3587;
	new2jer_param[11][15]=-0.000508183;  new2jer_param_er[11][15]=0.00038616;
	new2jer_param[11][16]=-7.37601e-06;  new2jer_param_er[11][16]=2.18165e-05;
	new2jer_param[11][17]=11.2976;  new2jer_param_er[11][17]=3.40538;
	new2jer_param[11][18]=58.7582;  new2jer_param_er[11][18]=16.7415;
	new2jer_param[11][19]=73.3543;  new2jer_param_er[11][19]=41.2621;
	new2jer_param[11][20]=13.7691;  new2jer_param_er[11][20]=0.00206641;
	new2jer_param[11][21]=-7.30096;  new2jer_param_er[11][21]=0.00109591;
	new2jer_param[11][22]=1.25341;  new2jer_param_er[11][22]=0.963188;
	new2jer_param[11][23]=1.3521e+06;  new2jer_param_er[11][23]=319.963;
	new2jer_param[11][24]=-857213;  new2jer_param_er[11][24]=202.917;
	new2jer_param[12][0]=-4.44425;  new2jer_param_er[12][0]=1.7212;
	new2jer_param[12][1]=-1.45874e-05;  new2jer_param_er[12][1]=1.54716e-05;
	new2jer_param[12][2]=1.12838;  new2jer_param_er[12][2]=1.54617;
	new2jer_param[12][3]=49.3183;  new2jer_param_er[12][3]=1.71588;
	new2jer_param[12][4]=-316.532;  new2jer_param_er[12][4]=19.166;
	new2jer_param[12][5]=-0.000241368;  new2jer_param_er[12][5]=0.000318295;
	new2jer_param[12][6]=-8.56154e-05;  new2jer_param_er[12][6]=8.23535e-05;
	new2jer_param[12][7]=40.712;  new2jer_param_er[12][7]=17.1509;
	new2jer_param[12][8]=93.7904;  new2jer_param_er[12][8]=34.2986;
	new2jer_param[12][9]=183.419;  new2jer_param_er[12][9]=106.471;
	new2jer_param[12][10]=0.00100941;  new2jer_param_er[12][10]=0.000744039;
	new2jer_param[12][11]=0.000164103;  new2jer_param_er[12][11]=0.000108537;
	new2jer_param[12][12]=-39.4559;  new2jer_param_er[12][12]=25.0679;
	new2jer_param[12][13]=78.7969;  new2jer_param_er[12][13]=29.8539;
	new2jer_param[12][14]=111.055;  new2jer_param_er[12][14]=95.1736;
	new2jer_param[12][15]=-0.000462314;  new2jer_param_er[12][15]=0.00065739;
	new2jer_param[12][16]=-9.39592e-06;  new2jer_param_er[12][16]=2.6614e-05;
	new2jer_param[12][17]=11.5005;  new2jer_param_er[12][17]=4.42764;
	new2jer_param[12][18]=56.4271;  new2jer_param_er[12][18]=21.1917;
	new2jer_param[12][19]=80.3647;  new2jer_param_er[12][19]=82.9741;
	new2jer_param[12][20]=-0.0121976;  new2jer_param_er[12][20]=0.0534953;
	new2jer_param[12][21]=4.3661e-05;  new2jer_param_er[12][21]=0.000224644;
	new2jer_param[12][22]=69.6908;  new2jer_param_er[12][22]=19.6412;
	new2jer_param[12][23]=35.8245;  new2jer_param_er[12][23]=4.25445;
	new2jer_param[12][24]=68.53;  new2jer_param_er[12][24]=59.303;
	new2jer_param[13][0]=-2708.59;  new2jer_param_er[13][0]=2954.52;
	new2jer_param[13][1]=2.99885e-06;  new2jer_param_er[13][1]=3.84965e-05;
	new2jer_param[13][2]=0.86867;  new2jer_param_er[13][2]=3.99404;
	new2jer_param[13][3]=73.7363;  new2jer_param_er[13][3]=22.1834;
	new2jer_param[13][4]=-1017.54;  new2jer_param_er[13][4]=81.2896;
	new2jer_param[13][5]=-1.01156;  new2jer_param_er[13][5]=1.46718;
	new2jer_param[13][6]=0.000168164;  new2jer_param_er[13][6]=1.7852e-05;
	new2jer_param[13][7]=12.753;  new2jer_param_er[13][7]=1.26857;
	new2jer_param[13][8]=11.1167;  new2jer_param_er[13][8]=4.7917;
	new2jer_param[13][9]=-15.6119;  new2jer_param_er[13][9]=16.3482;
	new2jer_param[13][10]=-0.000165482;  new2jer_param_er[13][10]=1.84846e-05;
	new2jer_param[13][11]=-0.000379693;  new2jer_param_er[13][11]=0.000301507;
	new2jer_param[13][12]=-9.50198;  new2jer_param_er[13][12]=0.594789;
	new2jer_param[13][13]=-44.6713;  new2jer_param_er[13][13]=17.052;
	new2jer_param[13][14]=243.04;  new2jer_param_er[13][14]=55.4226;
	new2jer_param[13][15]=-0.533868;  new2jer_param_er[13][15]=7.46394e-05;
	new2jer_param[13][16]=0.528462;  new2jer_param_er[13][16]=7.38634e-05;
	new2jer_param[13][17]=8.89243;  new2jer_param_er[13][17]=0.456141;
	new2jer_param[13][18]=-766514;  new2jer_param_er[13][18]=113912;
	new2jer_param[13][19]=8380.14;  new2jer_param_er[13][19]=102.893;
	new2jer_param[13][20]=-0.00181166;  new2jer_param_er[13][20]=0.00211858;
	new2jer_param[13][21]=-0.000243397;  new2jer_param_er[13][21]=0.000218896;
	new2jer_param[13][22]=69.9992;  new2jer_param_er[13][22]=31.0877;
	new2jer_param[13][23]=56.4935;  new2jer_param_er[13][23]=9.24641;
	new2jer_param[13][24]=150.92;  new2jer_param_er[13][24]=60.4588;
	new2jer_param[14][0]=0;  new2jer_param_er[14][0]=0;
	new2jer_param[14][1]=0;  new2jer_param_er[14][1]=0;
	new2jer_param[14][2]=0;  new2jer_param_er[14][2]=0;
	new2jer_param[14][3]=0;  new2jer_param_er[14][3]=0;
	new2jer_param[14][4]=0;  new2jer_param_er[14][4]=0;
	new2jer_param[14][5]=0;  new2jer_param_er[14][5]=0;
	new2jer_param[14][6]=0;  new2jer_param_er[14][6]=0;
	new2jer_param[14][7]=0;  new2jer_param_er[14][7]=0;
	new2jer_param[14][8]=0;  new2jer_param_er[14][8]=0;
	new2jer_param[14][9]=0;  new2jer_param_er[14][9]=0;
	new2jer_param[14][10]=7.47433e-05;  new2jer_param_er[14][10]=0.000317255;
	new2jer_param[14][11]=0.00159928;  new2jer_param_er[14][11]=0.00193466;
	new2jer_param[14][12]=-781.099;  new2jer_param_er[14][12]=876.806;
	new2jer_param[14][13]=184.541;  new2jer_param_er[14][13]=70.1105;
	new2jer_param[14][14]=786.746;  new2jer_param_er[14][14]=327.327;
	new2jer_param[14][15]=0.00106029;  new2jer_param_er[14][15]=0.00797798;
	new2jer_param[14][16]=5.75598e-05;  new2jer_param_er[14][16]=2.09668e-05;
	new2jer_param[14][17]=8.22004;  new2jer_param_er[14][17]=5.32011;
	new2jer_param[14][18]=64.8387;  new2jer_param_er[14][18]=39.2822;
	new2jer_param[14][19]=43.711;  new2jer_param_er[14][19]=106.089;
	new2jer_param[14][20]=0;  new2jer_param_er[14][20]=0;
	new2jer_param[14][21]=0;  new2jer_param_er[14][21]=0;
	new2jer_param[14][22]=0;  new2jer_param_er[14][22]=0;
	new2jer_param[14][23]=0;  new2jer_param_er[14][23]=0;
	new2jer_param[14][24]=0;  new2jer_param_er[14][24]=0;


	if(i<15 && j<25)
	{
		if(code==0) value=new2jer_param[i][j];
		if(code==1) value=new2jer_param_er[i][j];
	}
	return value;
}

//________________ returns NEW JER param or its uncertainty
double newJERparam(int code, int i, int j)
{
	double newjer_param[15][25]; // [eta][param]
	double newjer_param_er[15][25]; // [eta][param]
	double value=0.0;

	// NEW parameters 07-14-2009
	newjer_param[0][0]=0.119343;  newjer_param_er[0][0]=0.00314468;
	newjer_param[0][1]=-0.00105953;  newjer_param_er[0][1]=6.95625e-05;
	newjer_param[0][2]=3.45175e-05;  newjer_param_er[0][2]=1.36251e-06;
	newjer_param[0][3]=0.0306632;  newjer_param_er[0][3]=0.00256322;
	newjer_param[0][4]=-0.000131089;  newjer_param_er[0][4]=8.70229e-06;
	newjer_param[0][5]=0.311833;  newjer_param_er[0][5]=0.00150468;
	newjer_param[0][6]=-0.00675844;  newjer_param_er[0][6]=3.37299e-05;
	newjer_param[0][7]=7.65782e-05;  newjer_param_er[0][7]=6.69574e-07;
	newjer_param[0][8]=0.105356;  newjer_param_er[0][8]=0.00142981;
	newjer_param[0][9]=-0.000136616;  newjer_param_er[0][9]=5.03254e-06;
	newjer_param[0][10]=-0.0361686;  newjer_param_er[0][10]=0.0185456;
	newjer_param[0][11]=0.00128243;  newjer_param_er[0][11]=4.24584e-05;
	newjer_param[0][12]=-1.52357e-06;  newjer_param_er[0][12]=1.00586e-05;
	newjer_param[0][13]=0.0215701;  newjer_param_er[0][13]=0.00392601;
	newjer_param[0][14]=0.000136695;  newjer_param_er[0][14]=2.02636e-05;
	newjer_param[0][15]=0.191215;  newjer_param_er[0][15]=0.00918844;
	newjer_param[0][16]=-0.00364137;  newjer_param_er[0][16]=2.24066e-05;
	newjer_param[0][17]=3.83763e-05;  newjer_param_er[0][17]=4.91331e-06;
	newjer_param[0][18]=0.0896749;  newjer_param_er[0][18]=0.00165894;
	newjer_param[0][19]=7.85796e-06;  newjer_param_er[0][19]=7.98429e-06;
	newjer_param[0][20]=-0.0272588;  newjer_param_er[0][20]=0.0342609;
	newjer_param[0][21]=0.00449702;  newjer_param_er[0][21]=0.000111858;
	newjer_param[0][22]=-8.08661e-06;  newjer_param_er[0][22]=2.2715e-05;
	newjer_param[0][23]=0.0368499;  newjer_param_er[0][23]=0.0109694;
	newjer_param[0][24]=0.000367415;  newjer_param_er[0][24]=5.70766e-05;
	newjer_param[1][0]=0.0743446;  newjer_param_er[1][0]=0.0183135;
	newjer_param[1][1]=0.000247436;  newjer_param_er[1][1]=4.80135e-05;
	newjer_param[1][2]=2.87667e-06;  newjer_param_er[1][2]=8.85104e-06;
	newjer_param[1][3]=0.0904459;  newjer_param_er[1][3]=0.00177181;
	newjer_param[1][4]=-0.000188996;  newjer_param_er[1][4]=8.124e-06;
	newjer_param[1][5]=0.285438;  newjer_param_er[1][5]=0.00973598;
	newjer_param[1][6]=-0.00591245;  newjer_param_er[1][6]=2.57752e-05;
	newjer_param[1][7]=5.83689e-05;  newjer_param_er[1][7]=4.73984e-06;
	newjer_param[1][8]=0.104292;  newjer_param_er[1][8]=0.000972824;
	newjer_param[1][9]=-0.000109701;  newjer_param_er[1][9]=4.46752e-06;
	newjer_param[1][10]=-0.147575;  newjer_param_er[1][10]=0.0106425;
	newjer_param[1][11]=0.00551904;  newjer_param_er[1][11]=0.000567946;
	newjer_param[1][12]=-6.02653e-05;  newjer_param_er[1][12]=6.89239e-06;
	newjer_param[1][13]=-0.00799666;  newjer_param_er[1][13]=0.0012872;
	newjer_param[1][14]=1.88714e-05;  newjer_param_er[1][14]=5.32563e-06;
	newjer_param[1][15]=0.168263;  newjer_param_er[1][15]=0.00552407;
	newjer_param[1][16]=-0.00381363;  newjer_param_er[1][16]=0.000301486;
	newjer_param[1][17]=3.72841e-05;  newjer_param_er[1][17]=3.69896e-06;
	newjer_param[1][18]=0.0552846;  newjer_param_er[1][18]=0.000664559;
	newjer_param[1][19]=-1.34569e-05;  newjer_param_er[1][19]=2.6725e-06;
	newjer_param[1][20]=0.0494562;  newjer_param_er[1][20]=0.0278234;
	newjer_param[1][21]=0.00238938;  newjer_param_er[1][21]=0.00206893;
	newjer_param[1][22]=3.70568e-05;  newjer_param_er[1][22]=2.8803e-05;
	newjer_param[1][23]=0.235407;  newjer_param_er[1][23]=0.00688744;
	newjer_param[1][24]=-0.00045021;  newjer_param_er[1][24]=2.72024e-05;
	newjer_param[2][0]=-0.0275628;  newjer_param_er[2][0]=0.0193417;
	newjer_param[2][1]=0.00348816;  newjer_param_er[2][1]=4.0658e-05;
	newjer_param[2][2]=-2.58846e-05;  newjer_param_er[2][2]=9.08771e-06;
	newjer_param[2][3]=0.088725;  newjer_param_er[2][3]=0.00178948;
	newjer_param[2][4]=-0.000219068;  newjer_param_er[2][4]=7.83811e-06;
	newjer_param[2][5]=0.24067;  newjer_param_er[2][5]=0.0104227;
	newjer_param[2][6]=-0.0044073;  newjer_param_er[2][6]=2.16975e-05;
	newjer_param[2][7]=4.59682e-05;  newjer_param_er[2][7]=4.88015e-06;
	newjer_param[2][8]=0.105543;  newjer_param_er[2][8]=0.00101879;
	newjer_param[2][9]=-0.000120634;  newjer_param_er[2][9]=4.579e-06;
	newjer_param[2][10]=-0.175947;  newjer_param_er[2][10]=0.0147411;
	newjer_param[2][11]=0.00545286;  newjer_param_er[2][11]=4.14128e-05;
	newjer_param[2][12]=-5.21046e-05;  newjer_param_er[2][12]=7.58811e-06;
	newjer_param[2][13]=-0.00656557;  newjer_param_er[2][13]=0.00133665;
	newjer_param[2][14]=-2.77298e-06;  newjer_param_er[2][14]=5.36955e-06;
	newjer_param[2][15]=0.157946;  newjer_param_er[2][15]=0.00780927;
	newjer_param[2][16]=-0.00335858;  newjer_param_er[2][16]=2.26226e-05;
	newjer_param[2][17]=3.29693e-05;  newjer_param_er[2][17]=4.06204e-06;
	newjer_param[2][18]=0.0598815;  newjer_param_er[2][18]=0.000688747;
	newjer_param[2][19]=-2.43654e-05;  newjer_param_er[2][19]=2.70315e-06;
	newjer_param[2][20]=0.105951;  newjer_param_er[2][20]=0.0507535;
	newjer_param[2][21]=0.00108184;  newjer_param_er[2][21]=0.000209741;
	newjer_param[2][22]=4.83183e-05;  newjer_param_er[2][22]=3.42105e-05;
	newjer_param[2][23]=0.222155;  newjer_param_er[2][23]=0.00708712;
	newjer_param[2][24]=-0.000384935;  newjer_param_er[2][24]=2.76368e-05;
	newjer_param[3][0]=0.170694;  newjer_param_er[3][0]=0.0184655;
	newjer_param[3][1]=-0.00400434;  newjer_param_er[3][1]=5.12281e-05;
	newjer_param[3][2]=4.51275e-05;  newjer_param_er[3][2]=9.53991e-06;
	newjer_param[3][3]=0.094768;  newjer_param_er[3][3]=0.00222356;
	newjer_param[3][4]=-0.000184981;  newjer_param_er[3][4]=9.72646e-06;
	newjer_param[3][5]=0.340462;  newjer_param_er[3][5]=0.00977226;
	newjer_param[3][6]=-0.00796849;  newjer_param_er[3][6]=2.6721e-05;
	newjer_param[3][7]=7.85343e-05;  newjer_param_er[3][7]=5.09179e-06;
	newjer_param[3][8]=0.109677;  newjer_param_er[3][8]=0.0012721;
	newjer_param[3][9]=-0.000116944;  newjer_param_er[3][9]=5.67218e-06;
	newjer_param[3][10]=-0.203812;  newjer_param_er[3][10]=0.010309;
	newjer_param[3][11]=0.00677219;  newjer_param_er[3][11]=0.000524366;
	newjer_param[3][12]=-6.90694e-05;  newjer_param_er[3][12]=6.44028e-06;
	newjer_param[3][13]=0.00968912;  newjer_param_er[3][13]=0.00153634;
	newjer_param[3][14]=-2.1932e-05;  newjer_param_er[3][14]=6.07106e-06;
	newjer_param[3][15]=0.157946;  newjer_param_er[3][15]=0.0051451;
	newjer_param[3][16]=-0.00295073;  newjer_param_er[3][16]=0.000268408;
	newjer_param[3][17]=2.76299e-05;  newjer_param_er[3][17]=3.34027e-06;
	newjer_param[3][18]=0.0673421;  newjer_param_er[3][18]=0.000780854;
	newjer_param[3][19]=-4.24974e-05;  newjer_param_er[3][19]=3.03686e-06;
	newjer_param[3][20]=0.0861196;  newjer_param_er[3][20]=0.0238194;
	newjer_param[3][21]=2.85701e-05;  newjer_param_er[3][21]=0.00158696;
	newjer_param[3][22]=5.64572e-05;  newjer_param_er[3][22]=2.28046e-05;
	newjer_param[3][23]=0.200478;  newjer_param_er[3][23]=0.00701648;
	newjer_param[3][24]=-0.000353516;  newjer_param_er[3][24]=2.73206e-05;
	newjer_param[4][0]=0.11804;  newjer_param_er[4][0]=0.0184238;
	newjer_param[4][1]=-0.00232059;  newjer_param_er[4][1]=5.23777e-05;
	newjer_param[4][2]=4.23751e-05;  newjer_param_er[4][2]=1.05065e-05;
	newjer_param[4][3]=0.104218;  newjer_param_er[4][3]=0.00235565;
	newjer_param[4][4]=-0.000118992;  newjer_param_er[4][4]=9.97604e-06;
	newjer_param[4][5]=0.306236;  newjer_param_er[4][5]=0.00986915;
	newjer_param[4][6]=-0.00605338;  newjer_param_er[4][6]=2.76995e-05;
	newjer_param[4][7]=6.20639e-05;  newjer_param_er[4][7]=5.77336e-06;
	newjer_param[4][8]=0.125262;  newjer_param_er[4][8]=0.00135793;
	newjer_param[4][9]=-0.000106561;  newjer_param_er[4][9]=5.88961e-06;
	newjer_param[4][10]=-0.192637;  newjer_param_er[4][10]=0.013221;
	newjer_param[4][11]=0.00701698;  newjer_param_er[4][11]=4.63014e-05;
	newjer_param[4][12]=-7.554e-05;  newjer_param_er[4][12]=8.76423e-06;
	newjer_param[4][13]=0.0207201;  newjer_param_er[4][13]=0.00200845;
	newjer_param[4][14]=-3.80897e-05;  newjer_param_er[4][14]=7.4085e-06;
	newjer_param[4][15]=0.167815;  newjer_param_er[4][15]=0.00677706;
	newjer_param[4][16]=-0.00263561;  newjer_param_er[4][16]=2.53907e-05;
	newjer_param[4][17]=2.13312e-05;  newjer_param_er[4][17]=4.63341e-06;
	newjer_param[4][18]=0.0801541;  newjer_param_er[4][18]=0.00105189;
	newjer_param[4][19]=-5.20055e-05;  newjer_param_er[4][19]=3.84686e-06;
	newjer_param[4][20]=0.138705;  newjer_param_er[4][20]=0.03147;
	newjer_param[4][21]=-0.00334811;  newjer_param_er[4][21]=0.0020746;
	newjer_param[4][22]=9.2358e-05;  newjer_param_er[4][22]=2.96098e-05;
	newjer_param[4][23]=0.255835;  newjer_param_er[4][23]=0.00874165;
	newjer_param[4][24]=-0.000434759;  newjer_param_er[4][24]=3.08506e-05;
	newjer_param[5][0]=0.0744578;  newjer_param_er[5][0]=0.0246458;
	newjer_param[5][1]=0.000218351;  newjer_param_er[5][1]=4.40605e-05;
	newjer_param[5][2]=-6.74408e-06;  newjer_param_er[5][2]=7.54392e-07;
	newjer_param[5][3]=0.0866891;  newjer_param_er[5][3]=0.00163009;
	newjer_param[5][4]=-7.76345e-05;  newjer_param_er[5][4]=6.14488e-06;
	newjer_param[5][5]=0.297965;  newjer_param_er[5][5]=0.0142519;
	newjer_param[5][6]=-0.00507016;  newjer_param_er[5][6]=2.36136e-05;
	newjer_param[5][7]=4.56093e-05;  newjer_param_er[5][7]=4.06007e-07;
	newjer_param[5][8]=0.118217;  newjer_param_er[5][8]=0.000975398;
	newjer_param[5][9]=-7.34549e-05;  newjer_param_er[5][9]=3.74462e-06;
	newjer_param[5][10]=-0.230767;  newjer_param_er[5][10]=0.0214145;
	newjer_param[5][11]=0.00789996;  newjer_param_er[5][11]=4.30248e-05;
	newjer_param[5][12]=-9.88319e-05;  newjer_param_er[5][12]=7.74325e-07;
	newjer_param[5][13]=0.00343542;  newjer_param_er[5][13]=0.00250523;
	newjer_param[5][14]=-3.3107e-05;  newjer_param_er[5][14]=8.44783e-06;
	newjer_param[5][15]=0.169462;  newjer_param_er[5][15]=0.0109132;
	newjer_param[5][16]=-0.00262482;  newjer_param_er[5][16]=2.41364e-05;
	newjer_param[5][17]=1.59232e-05;  newjer_param_er[5][17]=4.34128e-07;
	newjer_param[5][18]=0.0742118;  newjer_param_er[5][18]=0.00131377;
	newjer_param[5][19]=-4.78425e-05;  newjer_param_er[5][19]=4.40894e-06;
	newjer_param[5][20]=0.176869;  newjer_param_er[5][20]=0.0678367;
	newjer_param[5][21]=-0.00259215;  newjer_param_er[5][21]=0.000196127;
	newjer_param[5][22]=2.45926e-05;  newjer_param_er[5][22]=3.83147e-06;
	newjer_param[5][23]=0.479562;  newjer_param_er[5][23]=0.0175298;
	newjer_param[5][24]=-0.000678321;  newjer_param_er[5][24]=5.60546e-05;
	newjer_param[6][0]=0.00938799;  newjer_param_er[6][0]=0.0214178;
	newjer_param[6][1]=0.00286763;  newjer_param_er[6][1]=5.49692e-05;
	newjer_param[6][2]=-4.64671e-05;  newjer_param_er[6][2]=9.12845e-07;
	newjer_param[6][3]=0.110373;  newjer_param_er[6][3]=0.00285928;
	newjer_param[6][4]=-0.000322779;  newjer_param_er[6][4]=1.1545e-05;
	newjer_param[6][5]=0.281354;  newjer_param_er[6][5]=0.011295;
	newjer_param[6][6]=-0.00413676;  newjer_param_er[6][6]=2.89717e-05;
	newjer_param[6][7]=2.84304e-05;  newjer_param_er[6][7]=4.77439e-07;
	newjer_param[6][8]=0.137881;  newjer_param_er[6][8]=0.00162782;
	newjer_param[6][9]=-0.00019431;  newjer_param_er[6][9]=7.11215e-06;
	newjer_param[6][10]=-0.261284;  newjer_param_er[6][10]=0.0164498;
	newjer_param[6][11]=0.00726371;  newjer_param_er[6][11]=4.72184e-05;
	newjer_param[6][12]=-7.91647e-05;  newjer_param_er[6][12]=1.01505e-05;
	newjer_param[6][13]=-0.0219306;  newjer_param_er[6][13]=0.00243926;
	newjer_param[6][14]=9.96534e-05;  newjer_param_er[6][14]=1.15067e-05;
	newjer_param[6][15]=0.15611;  newjer_param_er[6][15]=0.0087093;
	newjer_param[6][16]=-0.00301021;  newjer_param_er[6][16]=2.57324e-05;
	newjer_param[6][17]=2.65328e-05;  newjer_param_er[6][17]=5.41176e-06;
	newjer_param[6][18]=0.0667241;  newjer_param_er[6][18]=0.00122799;
	newjer_param[6][19]=9.27018e-06;  newjer_param_er[6][19]=5.558e-06;
	newjer_param[6][20]=0.332473;  newjer_param_er[6][20]=0.0524506;
	newjer_param[6][21]=-0.0139085;  newjer_param_er[6][21]=0.000198343;
	newjer_param[6][22]=0.000222125;  newjer_param_er[6][22]=3.63124e-05;
	newjer_param[6][23]=0.168506;  newjer_param_er[6][23]=0.00873082;
	newjer_param[6][24]=-0.00037344;  newjer_param_er[6][24]=3.87933e-05;
	newjer_param[7][0]=-0.021659;  newjer_param_er[7][0]=0.0183932;
	newjer_param[7][1]=0.00176846;  newjer_param_er[7][1]=6.69075e-05;
	newjer_param[7][2]=-1.70584e-05;  newjer_param_er[7][2]=1.08547e-05;
	newjer_param[7][3]=0.0550046;  newjer_param_er[7][3]=0.00234551;
	newjer_param[7][4]=-0.000172672;  newjer_param_er[7][4]=1.18965e-05;
	newjer_param[7][5]=0.282372;  newjer_param_er[7][5]=0.00986023;
	newjer_param[7][6]=-0.00502482;  newjer_param_er[7][6]=3.67484e-05;
	newjer_param[7][7]=4.8236e-05;  newjer_param_er[7][7]=6.04829e-06;
	newjer_param[7][8]=0.115812;  newjer_param_er[7][8]=0.00133564;
	newjer_param[7][9]=-0.000179829;  newjer_param_er[7][9]=6.89609e-06;
	newjer_param[7][10]=-0.355076;  newjer_param_er[7][10]=0.0175825;
	newjer_param[7][11]=0.0126345;  newjer_param_er[7][11]=6.849e-05;
	newjer_param[7][12]=-0.000141083;  newjer_param_er[7][12]=1.11936e-05;
	newjer_param[7][13]=-0.0399519;  newjer_param_er[7][13]=0.00299535;
	newjer_param[7][14]=0.000117128;  newjer_param_er[7][14]=1.67435e-05;
	newjer_param[7][15]=0.1295;  newjer_param_er[7][15]=0.00989812;
	newjer_param[7][16]=-0.00132418;  newjer_param_er[7][16]=3.78495e-05;
	newjer_param[7][17]=8.74715e-06;  newjer_param_er[7][17]=6.22161e-06;
	newjer_param[7][18]=0.0528986;  newjer_param_er[7][18]=0.00157671;
	newjer_param[7][19]=-2.67672e-05;  newjer_param_er[7][19]=8.62368e-06;
	newjer_param[7][20]=0.362403;  newjer_param_er[7][20]=0.067433;
	newjer_param[7][21]=-0.0112584;  newjer_param_er[7][21]=0.00400909;
	newjer_param[7][22]=0.000155376;  newjer_param_er[7][22]=5.09827e-05;
	newjer_param[7][23]=0.27005;  newjer_param_er[7][23]=0.0196289;
	newjer_param[7][24]=-0.000376317;  newjer_param_er[7][24]=0.00011391;
	newjer_param[8][0]=-0.147359;  newjer_param_er[8][0]=0.030064;
	newjer_param[8][1]=0.00860543;  newjer_param_er[8][1]=9.03519e-05;
	newjer_param[8][2]=-0.00012039;  newjer_param_er[8][2]=1.40453e-05;
	newjer_param[8][3]=0.0405336;  newjer_param_er[8][3]=0.00134233;
	newjer_param[8][4]=-6.27006e-05;  newjer_param_er[8][4]=5.70287e-06;
	newjer_param[8][5]=0.357916;  newjer_param_er[8][5]=0.0167814;
	newjer_param[8][6]=-0.00855328;  newjer_param_er[8][6]=5.21371e-05;
	newjer_param[8][7]=8.19137e-05;  newjer_param_er[8][7]=8.08121e-06;
	newjer_param[8][8]=0.104416;  newjer_param_er[8][8]=0.000827841;
	newjer_param[8][9]=-0.000140373;  newjer_param_er[8][9]=3.57486e-06;
	newjer_param[8][10]=-0.340401;  newjer_param_er[8][10]=0.0332533;
	newjer_param[8][11]=0.0106851;  newjer_param_er[8][11]=0.000115048;
	newjer_param[8][12]=-0.000118701;  newjer_param_er[8][12]=1.85101e-05;
	newjer_param[8][13]=-0.0537098;  newjer_param_er[8][13]=0.00254475;
	newjer_param[8][14]=0.000105504;  newjer_param_er[8][14]=1.11912e-05;
	newjer_param[8][15]=0.165392;  newjer_param_er[8][15]=0.0173911;
	newjer_param[8][16]=-0.0031747;  newjer_param_er[8][16]=6.03268e-05;
	newjer_param[8][17]=3.15669e-05;  newjer_param_er[8][17]=9.72956e-06;
	newjer_param[8][18]=0.0450839;  newjer_param_er[8][18]=0.00128918;
	newjer_param[8][19]=-5.9899e-05;  newjer_param_er[8][19]=5.55408e-06;
	newjer_param[8][20]=0.291315;  newjer_param_er[8][20]=0.145622;
	newjer_param[8][21]=-0.00487384;  newjer_param_er[8][21]=0.000420369;
	newjer_param[8][22]=1.11175e-05;  newjer_param_er[8][22]=0.0001095;
	newjer_param[8][23]=0.452202;  newjer_param_er[8][23]=0.0306732;
	newjer_param[8][24]=-0.00016021;  newjer_param_er[8][24]=0.000159679;
	newjer_param[9][0]=-0.299768;  newjer_param_er[9][0]=0.0536478;
	newjer_param[9][1]=0.0158593;  newjer_param_er[9][1]=0.000136572;
	newjer_param[9][2]=-0.00021132;  newjer_param_er[9][2]=2.47158e-05;
	newjer_param[9][3]=0.0301886;  newjer_param_er[9][3]=0.00206139;
	newjer_param[9][4]=-1.59033e-05;  newjer_param_er[9][4]=9.69668e-06;
	newjer_param[9][5]=0.229547;  newjer_param_er[9][5]=0.0311322;
	newjer_param[9][6]=-0.0014667;  newjer_param_er[9][6]=8.05009e-05;
	newjer_param[9][7]=-2.875e-07;  newjer_param_er[9][7]=1.45449e-05;
	newjer_param[9][8]=0.117019;  newjer_param_er[9][8]=0.0012277;
	newjer_param[9][9]=-0.000184489;  newjer_param_er[9][9]=5.82147e-06;
	newjer_param[9][10]=-0.476516;  newjer_param_er[9][10]=0.0607404;
	newjer_param[9][11]=0.0153612;  newjer_param_er[9][11]=0.000125725;
	newjer_param[9][12]=-0.000159452;  newjer_param_er[9][12]=3.31762e-05;
	newjer_param[9][13]=-0.0590107;  newjer_param_er[9][13]=0.0034089;
	newjer_param[9][14]=0.000149436;  newjer_param_er[9][14]=1.73205e-05;
	newjer_param[9][15]=0.116423;  newjer_param_er[9][15]=0.0312739;
	newjer_param[9][16]=-0.00106857;  newjer_param_er[9][16]=6.67136e-05;
	newjer_param[9][17]=1.01815e-05;  newjer_param_er[9][17]=1.70455e-05;
	newjer_param[9][18]=0.0502626;  newjer_param_er[9][18]=0.00177539;
	newjer_param[9][19]=-8.62682e-05;  newjer_param_er[9][19]=8.83492e-06;
	newjer_param[9][20]=0.819459;  newjer_param_er[9][20]=0.0137594;
	newjer_param[9][21]=-0.0314006;  newjer_param_er[9][21]=0.000374014;
	newjer_param[9][22]=0.000308394;  newjer_param_er[9][22]=9.04288e-06;
	newjer_param[9][23]=0.325153;  newjer_param_er[9][23]=0.0360778;
	newjer_param[9][24]=0.000131873;  newjer_param_er[9][24]=0.000206635;
	newjer_param[10][0]=-0.953732;  newjer_param_er[10][0]=0.00869531;
	newjer_param[10][1]=0.0410471;  newjer_param_er[10][1]=0.000196105;
	newjer_param[10][2]=-0.000444301;  newjer_param_er[10][2]=3.91545e-06;
	newjer_param[10][3]=0.0307947;  newjer_param_er[10][3]=0.00219417;
	newjer_param[10][4]=-5.3727e-05;  newjer_param_er[10][4]=9.27327e-06;
	newjer_param[10][5]=0.125666;  newjer_param_er[10][5]=0.00508167;
	newjer_param[10][6]=0.00192099;  newjer_param_er[10][6]=0.00011747;
	newjer_param[10][7]=-1.90869e-05;  newjer_param_er[10][7]=2.39803e-06;
	newjer_param[10][8]=0.117313;  newjer_param_er[10][8]=0.00136483;
	newjer_param[10][9]=-0.000169817;  newjer_param_er[10][9]=5.93895e-06;
	newjer_param[10][10]=-0.268869;  newjer_param_er[10][10]=0.00658315;
	newjer_param[10][11]=0.00743125;  newjer_param_er[10][11]=0.000153269;
	newjer_param[10][12]=-9.31827e-05;  newjer_param_er[10][12]=3.2684e-06;
	newjer_param[10][13]=-0.0614147;  newjer_param_er[10][13]=0.00327184;
	newjer_param[10][14]=0.000104405;  newjer_param_er[10][14]=1.53567e-05;
	newjer_param[10][15]=0.22591;  newjer_param_er[10][15]=0.00342316;
	newjer_param[10][16]=-0.00542276;  newjer_param_er[10][16]=7.92462e-05;
	newjer_param[10][17]=5.56294e-05;  newjer_param_er[10][17]=1.6729e-06;
	newjer_param[10][18]=0.0491936;  newjer_param_er[10][18]=0.00167363;
	newjer_param[10][19]=-6.25096e-05;  newjer_param_er[10][19]=7.68511e-06;
	newjer_param[10][20]=-0.403697;  newjer_param_er[10][20]=0.0180041;
	newjer_param[10][21]=0.0298899;  newjer_param_er[10][21]=0.000429542;
	newjer_param[10][22]=-0.000414844;  newjer_param_er[10][22]=9.54328e-06;
	newjer_param[10][23]=0.28128;  newjer_param_er[10][23]=0.0329415;
	newjer_param[10][24]=0.000116748;  newjer_param_er[10][24]=0.000171607;
	newjer_param[11][0]=-0.721504;  newjer_param_er[11][0]=0.0123063;
	newjer_param[11][1]=0.0285797;  newjer_param_er[11][1]=0.000254757;
	newjer_param[11][2]=-0.000301013;  newjer_param_er[11][2]=4.87484e-06;
	newjer_param[11][3]=-0.0105836;  newjer_param_er[11][3]=0.00269672;
	newjer_param[11][4]=6.65018e-05;  newjer_param_er[11][4]=1.06714e-05;
	newjer_param[11][5]=0.383707;  newjer_param_er[11][5]=0.00728806;
	newjer_param[11][6]=-0.00946731;  newjer_param_er[11][6]=0.000154518;
	newjer_param[11][7]=0.000108082;  newjer_param_er[11][7]=3.01887e-06;
	newjer_param[11][8]=0.110545;  newjer_param_er[11][8]=0.00170017;
	newjer_param[11][9]=-0.000142146;  newjer_param_er[11][9]=6.9127e-06;
	newjer_param[11][10]=-0.610956;  newjer_param_er[11][10]=0.00768133;
	newjer_param[11][11]=0.0191861;  newjer_param_er[11][11]=0.000172439;
	newjer_param[11][12]=-0.000196723;  newjer_param_er[11][12]=3.61353e-06;
	newjer_param[11][13]=-0.0809709;  newjer_param_er[11][13]=0.00406365;
	newjer_param[11][14]=0.00014617;  newjer_param_er[11][14]=1.73409e-05;
	newjer_param[11][15]=0.180601;  newjer_param_er[11][15]=0.00386761;
	newjer_param[11][16]=-0.00436925;  newjer_param_er[11][16]=8.70863e-05;
	newjer_param[11][17]=5.30058e-05;  newjer_param_er[11][17]=1.82574e-06;
	newjer_param[11][18]=0.0516735;  newjer_param_er[11][18]=0.00205294;
	newjer_param[11][19]=-6.66337e-05;  newjer_param_er[11][19]=8.61556e-06;
	newjer_param[11][20]=-0.882516;  newjer_param_er[11][20]=0.0179587;
	newjer_param[11][21]=0.0451445;  newjer_param_er[11][21]=0.000415322;
	newjer_param[11][22]=-0.000517426;  newjer_param_er[11][22]=9.05856e-06;
	newjer_param[11][23]=0.196388;  newjer_param_er[11][23]=0.0366708;
	newjer_param[11][24]=0.00033582;  newjer_param_er[11][24]=0.0001734;
	newjer_param[12][0]=-0.0417805;  newjer_param_er[12][0]=0.010778;
	newjer_param[12][1]=0.00107649;  newjer_param_er[12][1]=0.00021101;
	newjer_param[12][2]=-6.06845e-05;  newjer_param_er[12][2]=3.99126e-06;
	newjer_param[12][3]=-0.0556253;  newjer_param_er[12][3]=0.00408475;
	newjer_param[12][4]=0.000160627;  newjer_param_er[12][4]=1.62227e-05;
	newjer_param[12][5]=1.39369;  newjer_param_er[12][5]=0.00590008;
	newjer_param[12][6]=-0.0473597;  newjer_param_er[12][6]=0.000117664;
	newjer_param[12][7]=0.000450922;  newjer_param_er[12][7]=2.28152e-06;
	newjer_param[12][8]=0.102307;  newjer_param_er[12][8]=0.00257808;
	newjer_param[12][9]=-0.000113637;  newjer_param_er[12][9]=1.03826e-05;
	newjer_param[12][10]=-1.81198;  newjer_param_er[12][10]=0.00882027;
	newjer_param[12][11]=0.0585738;  newjer_param_er[12][11]=0.000165491;
	newjer_param[12][12]=-0.000541543;  newjer_param_er[12][12]=2.99558e-06;
	newjer_param[12][13]=-0.124696;  newjer_param_er[12][13]=0.006;
	newjer_param[12][14]=0.000249786;  newjer_param_er[12][14]=2.38289e-05;
	newjer_param[12][15]=-0.0477552;  newjer_param_er[12][15]=0.00463068;
	newjer_param[12][16]=0.0048045;  newjer_param_er[12][16]=8.71889e-05;
	newjer_param[12][17]=-4.57489e-05;  newjer_param_er[12][17]=1.59042e-06;
	newjer_param[12][18]=0.0556088;  newjer_param_er[12][18]=0.00299674;
	newjer_param[12][19]=-7.86181e-05;  newjer_param_er[12][19]=1.1584e-05;
	newjer_param[12][20]=-1.49496;  newjer_param_er[12][20]=0.0268586;
	newjer_param[12][21]=0.0732949;  newjer_param_er[12][21]=0.000502876;
	newjer_param[12][22]=-0.000831473;  newjer_param_er[12][22]=9.13321e-06;
	newjer_param[12][23]=0.289883;  newjer_param_er[12][23]=0.0557916;
	newjer_param[12][24]=-0.000132964;  newjer_param_er[12][24]=0.000224823;
	newjer_param[13][0]=-2.94925;  newjer_param_er[13][0]=0.0262199;
	newjer_param[13][1]=0.085559;  newjer_param_er[13][1]=0.000430568;
	newjer_param[13][2]=-0.000628295;  newjer_param_er[13][2]=6.9836e-06;
	newjer_param[13][3]=-0.0540865;  newjer_param_er[13][3]=0.0136904;
	newjer_param[13][4]=0.000139313;  newjer_param_er[13][4]=4.936e-05;
	newjer_param[13][5]=0.867206;  newjer_param_er[13][5]=0.0151517;
	newjer_param[13][6]=-0.0311954;  newjer_param_er[13][6]=0.000248982;
	newjer_param[13][7]=0.000338698;  newjer_param_er[13][7]=4.04128e-06;
	newjer_param[13][8]=0.094809;  newjer_param_er[13][8]=0.00915288;
	newjer_param[13][9]=-3.07479e-06;  newjer_param_er[13][9]=3.34433e-05;
	newjer_param[13][10]=-0.909461;  newjer_param_er[13][10]=0.0136892;
	newjer_param[13][11]=0.0269072;  newjer_param_er[13][11]=0.00022884;
	newjer_param[13][12]=-0.000263558;  newjer_param_er[13][12]=3.78923e-06;
	newjer_param[13][13]=-0.120874;  newjer_param_er[13][13]=0.00759237;
	newjer_param[13][14]=0.000149562;  newjer_param_er[13][14]=2.68079e-05;
	newjer_param[13][15]=1.39308;  newjer_param_er[13][15]=0.00713046;
	newjer_param[13][16]=-0.0432826;  newjer_param_er[13][16]=0.000118901;
	newjer_param[13][17]=0.000363841;  newjer_param_er[13][17]=1.96356e-06;
	newjer_param[13][18]=0.0615054;  newjer_param_er[13][18]=0.00404197;
	newjer_param[13][19]=-6.12759e-05;  newjer_param_er[13][19]=1.42769e-05;
	newjer_param[13][20]=-2.45903;  newjer_param_er[13][20]=0.0337559;
	newjer_param[13][21]=0.080548;  newjer_param_er[13][21]=0.000562604;
	newjer_param[13][22]=-0.000635591;  newjer_param_er[13][22]=9.28515e-06;
	newjer_param[13][23]=0.123775;  newjer_param_er[13][23]=0.0329547;
	newjer_param[13][24]=-6.06731e-05;  newjer_param_er[13][24]=0.00011711;
	newjer_param[14][0]=0;  newjer_param_er[14][0]=0;
	newjer_param[14][1]=0;  newjer_param_er[14][1]=0;
	newjer_param[14][2]=0;  newjer_param_er[14][2]=0;
	newjer_param[14][3]=0;  newjer_param_er[14][3]=0;
	newjer_param[14][4]=0;  newjer_param_er[14][4]=0;
	newjer_param[14][5]=0;  newjer_param_er[14][5]=0;
	newjer_param[14][6]=0;  newjer_param_er[14][6]=0;
	newjer_param[14][7]=0;  newjer_param_er[14][7]=0;
	newjer_param[14][8]=0;  newjer_param_er[14][8]=0;
	newjer_param[14][9]=0;  newjer_param_er[14][9]=0;
	newjer_param[14][10]=-430.789;  newjer_param_er[14][10]=0.377943;
	newjer_param[14][11]=8.10679;  newjer_param_er[14][11]=0.00358588;
	newjer_param[14][12]=-0.0382566;  newjer_param_er[14][12]=3.38665e-05;
	newjer_param[14][13]=-0.143086;  newjer_param_er[14][13]=0.00730127;
	newjer_param[14][14]=6.88638e-05;  newjer_param_er[14][14]=1.81481e-05;
	newjer_param[14][15]=134.235;  newjer_param_er[14][15]=0.218195;
	newjer_param[14][16]=-2.62041;  newjer_param_er[14][16]=0.0020727;
	newjer_param[14][17]=0.0129575;  newjer_param_er[14][17]=1.96062e-05;
	newjer_param[14][18]=0.0746147;  newjer_param_er[14][18]=0.00451432;
	newjer_param[14][19]=-6.99466e-05;  newjer_param_er[14][19]=1.1196e-05;
	newjer_param[14][20]=0;  newjer_param_er[14][20]=0;
	newjer_param[14][21]=0;  newjer_param_er[14][21]=0;
	newjer_param[14][22]=0;  newjer_param_er[14][22]=0;
	newjer_param[14][23]=0;  newjer_param_er[14][23]=0;
	newjer_param[14][24]=0;  newjer_param_er[14][24]=0;

	if(i<15 && j<25)
	{
		if(code==0) value=newjer_param[i][j];
		if(code==1) value=newjer_param_er[i][j];
	}

	return value;
}



//________________ returns old JER param or its uncertainty
double JERparam(int code, int i, int j)
{
	double jer_param[15][13]; // [eta][param]
	double jer_param_er[15][13]; // [eta][param]

	double value=0.0;

	//below is the original JER parameters
	//I am going to replace this with a fit parameters
	//for E>60GeV region for a test. E<60GeV is going to
	//use the exact bin values.
	
/*	jer_param[0][0]=-0.00505533;  jer_param_er[0][0]=0.00366313;
	jer_param[0][1]=-7.37332E-5;  jer_param_er[0][1]=1.02366E-5;
	jer_param[0][2]=5.06548;  jer_param_er[0][2]=0.198112;
	jer_param[0][3]=1.03604;  jer_param_er[0][3]=0.0162871;
	jer_param[0][4]=0.000363984;  jer_param_er[0][4]=9.99763E-5;
	jer_param[0][5]=0.0404354;  jer_param_er[0][5]=0.00394954;
	jer_param[0][6]=7.70406E-5;  jer_param_er[0][6]=1.80066E-5;
	jer_param[0][7]=-1.05943;  jer_param_er[0][7]=0.148324;
	jer_param[0][8]=0.157408;  jer_param_er[0][8]=0.00948233;
	jer_param[0][9]=0.00745406;  jer_param_er[0][9]=0.000132492;
	jer_param[0][10]=-4.48426;  jer_param_er[0][10]=0.997815;
	jer_param[0][11]=1.14401;  jer_param_er[0][11]=0.269919;
	jer_param[0][12]=0.0459092;  jer_param_er[0][12]=0.0168848;
	jer_param[1][0]=0.0928237;  jer_param_er[1][0]=0.00248223;
	jer_param[1][1]=-0.000201847;  jer_param_er[1][1]=8.82907E-6;
	jer_param[1][2]=0.177168;  jer_param_er[1][2]=0.137301;
	jer_param[1][3]=0.670063;  jer_param_er[1][3]=0.0119732;
	jer_param[1][4]=0.00308811;  jer_param_er[1][4]=9.34083E-5;
	jer_param[1][5]=0.000336138;  jer_param_er[1][5]=0.00152807;
	jer_param[1][6]=5.34059E-6;  jer_param_er[1][6]=5.37503E-6;
	jer_param[1][7]=-1.15979;  jer_param_er[1][7]=0.0761452;
	jer_param[1][8]=0.137285;  jer_param_er[1][8]=0.00491652;
	jer_param[1][9]=0.0020658;  jer_param_er[1][9]=3.46802E-5;
	jer_param[1][10]=-14.7515;  jer_param_er[1][10]=0.612517;
	jer_param[1][11]=4.80581;  jer_param_er[1][11]=0.195418;
	jer_param[1][12]=-0.127518;  jer_param_er[1][12]=0.0109911;
	jer_param[2][0]=0.100634;  jer_param_er[2][0]=0.00241696;
	jer_param[2][1]=-0.000253429;  jer_param_er[2][1]=8.43763E-6;
	jer_param[2][2]=-0.688927;  jer_param_er[2][2]=0.125514;
	jer_param[2][3]=0.655545;  jer_param_er[2][3]=0.0103474;
	jer_param[2][4]=0.00300742;  jer_param_er[2][4]=8.73698E-5;
	jer_param[2][5]=0.0105033;  jer_param_er[2][5]=0.00178783;
	jer_param[2][6]=-3.76786E-5;  jer_param_er[2][6]=5.83178E-6;
	jer_param[2][7]=-1.85956;  jer_param_er[2][7]=0.096736;
	jer_param[2][8]=0.139571;  jer_param_er[2][8]=0.00510727;
	jer_param[2][9]=0.00227237;  jer_param_er[2][9]=3.75247E-5;
	jer_param[2][10]=-19.2885;  jer_param_er[2][10]=1.05958;
	jer_param[2][11]=5.66883;  jer_param_er[2][11]=0.262952;
	jer_param[2][12]=-0.161499;  jer_param_er[2][12]=0.013508;
	jer_param[3][0]=0.0952545;  jer_param_er[3][0]=0.00263113;
	jer_param[3][1]=-0.000189574;  jer_param_er[3][1]=9.63887E-6;
	jer_param[3][2]=0.236287;  jer_param_er[3][2]=0.12846;
	jer_param[3][3]=0.731283;  jer_param_er[3][3]=0.0128143;
	jer_param[3][4]=0.00341218;  jer_param_er[3][4]=0.00011154;
	jer_param[3][5]=0.0258263;  jer_param_er[3][5]=0.0016674;
	jer_param[3][6]=-4.83423E-5;  jer_param_er[3][6]=5.80565E-6;
	jer_param[3][7]=-2.22975;  jer_param_er[3][7]=0.0737835;
	jer_param[3][8]=0.179229;  jer_param_er[3][8]=0.00515981;
	jer_param[3][9]=0.0024204;  jer_param_er[3][9]=4.113E-5;
	jer_param[3][10]=-10.4553;  jer_param_er[3][10]=0.561997;
	jer_param[3][11]=3.36745;  jer_param_er[3][11]=0.171045;
	jer_param[3][12]=-0.0618917;  jer_param_er[3][12]=0.00980421;
	jer_param[4][0]=0.115229;  jer_param_er[4][0]=0.0027871;
	jer_param[4][1]=-0.000157746;  jer_param_er[4][1]=9.95943E-6;
	jer_param[4][2]=-0.0905143;  jer_param_er[4][2]=0.134539;
	jer_param[4][3]=0.795373;  jer_param_er[4][3]=0.016243;
	jer_param[4][4]=0.00635758;  jer_param_er[4][4]=0.000146956;
	jer_param[4][5]=0.0345064;  jer_param_er[4][5]=0.00216488;
	jer_param[4][6]=-5.72569E-5;  jer_param_er[4][6]=7.07075E-6;
	jer_param[4][7]=-2.21206;  jer_param_er[4][7]=0.0958281;
	jer_param[4][8]=0.212192;  jer_param_er[4][8]=0.00725209;
	jer_param[4][9]=0.00349351;  jer_param_er[4][9]=6.06418E-5;
	jer_param[4][10]=-10.6021;  jer_param_er[4][10]=0.660483;
	jer_param[4][11]=3.23948;  jer_param_er[4][11]=0.191871;
	jer_param[4][12]=-0.0236045;  jer_param_er[4][12]=0.0106671;
	jer_param[5][0]=0.0900639;  jer_param_er[5][0]=0.00267652;
	jer_param[5][1]=-8.17704E-5;  jer_param_er[5][1]=7.7302E-6;
	jer_param[5][2]=-0.567021;  jer_param_er[5][2]=0.160166;
	jer_param[5][3]=0.706037;  jer_param_er[5][3]=0.0137402;
	jer_param[5][4]=0.0066424;  jer_param_er[5][4]=0.000110101;
	jer_param[5][5]=0.0161473;  jer_param_er[5][5]=0.00313166;
	jer_param[5][6]=-3.72664E-5;  jer_param_er[5][6]=8.95092E-6;
	jer_param[5][7]=-3.40086;  jer_param_er[5][7]=0.154954;
	jer_param[5][8]=0.163174;  jer_param_er[5][8]=0.00755835;
	jer_param[5][9]=0.00291679;  jer_param_er[5][9]=6.87186E-5;
	jer_param[5][10]=-10.8209;  jer_param_er[5][10]=1.39492;
	jer_param[5][11]=2.05925;  jer_param_er[5][11]=0.352187;
	jer_param[5][12]=0.184539;  jer_param_er[5][12]=0.0192903;
	jer_param[6][0]=0.110669;  jer_param_er[6][0]=0.00336557;
	jer_param[6][1]=-0.000302301;  jer_param_er[6][1]=1.15245E-5;
	jer_param[6][2]=-1.70381;  jer_param_er[6][2]=0.17128;
	jer_param[6][3]=0.84678;  jer_param_er[6][3]=0.0199186;
	jer_param[6][4]=0.00418531;  jer_param_er[6][4]=0.00022283;
	jer_param[6][5]=-0.00347175;  jer_param_er[6][5]=0.0027646;
	jer_param[6][6]=8.50977E-5;  jer_param_er[6][6]=1.10633E-5;
	jer_param[6][7]=-3.2283;  jer_param_er[6][7]=0.126136;
	jer_param[6][8]=0.0758951;  jer_param_er[6][8]=0.00708642;
	jer_param[6][9]=0.00415232;  jer_param_er[6][9]=8.11982E-5;
	jer_param[6][10]=-16.6687;  jer_param_er[6][10]=1.09606;
	jer_param[6][11]=5.18449;  jer_param_er[6][11]=0.288507;
	jer_param[6][12]=-0.190928;  jer_param_er[6][12]=0.016249;
	jer_param[7][0]=0.0658253;  jer_param_er[7][0]=0.00272349;
	jer_param[7][1]=-0.000196113;  jer_param_er[7][1]=1.14286E-5;
	jer_param[7][2]=-1.07207;  jer_param_er[7][2]=0.137955;
	jer_param[7][3]=0.889453;  jer_param_er[7][3]=0.015712;
	jer_param[7][4]=0.00152891;  jer_param_er[7][4]=0.000140114;
	jer_param[7][5]=-0.0287057;  jer_param_er[7][5]=0.00306348;
	jer_param[7][6]=0.000128678;  jer_param_er[7][6]=1.50299E-5;
	jer_param[7][7]=-2.67266;  jer_param_er[7][7]=0.129521;
	jer_param[7][8]=0.170338;  jer_param_er[7][8]=0.00800203;
	jer_param[7][9]=0.00119121;  jer_param_er[7][9]=7.89355E-5;
	jer_param[7][10]=-8.77829;  jer_param_er[7][10]=1.94651;
	jer_param[7][11]=2.5475;  jer_param_er[7][11]=0.581462;
	jer_param[7][12]=0.0673843;  jer_param_er[7][12]=0.0379263;
	jer_param[8][0]=0.058632;  jer_param_er[8][0]=0.00246133;
	jer_param[8][1]=-8.6882E-5;  jer_param_er[8][1]=7.27552E-6;
	jer_param[8][2]=-2.81222;  jer_param_er[8][2]=0.177198;
	jer_param[8][3]=0.87612;  jer_param_er[8][3]=0.0135608;
	jer_param[8][4]=0.000862571;  jer_param_er[8][4]=8.72723E-5;
	jer_param[8][5]=-0.0272999;  jer_param_er[8][5]=0.00336864;
	jer_param[8][6]=5.7274E-5;  jer_param_er[8][6]=1.19158E-5;
	jer_param[8][7]=-3.16373;  jer_param_er[8][7]=0.187887;
	jer_param[8][8]=0.218132;  jer_param_er[8][8]=0.00804389;
	jer_param[8][9]=-0.000115495;  jer_param_er[8][9]=5.26047E-5;
	jer_param[8][10]=6.59552;  jer_param_er[8][10]=4.14997;
	jer_param[8][11]=-3.91415;  jer_param_er[8][11]=1.16124;
	jer_param[8][12]=0.689584;  jer_param_er[8][12]=0.0728198;
	jer_param[9][0]=0.047587;  jer_param_er[9][0]=0.00405247;
	jer_param[9][1]=-3.67703E-5;  jer_param_er[9][1]=1.32232E-5;
	jer_param[9][2]=-2.82013;  jer_param_er[9][2]=0.277916;
	jer_param[9][3]=0.90659;  jer_param_er[9][3]=0.0207463;
	jer_param[9][4]=0.00138775;  jer_param_er[9][4]=0.000128634;
	jer_param[9][5]=-0.0150813;  jer_param_er[9][5]=0.00500665;
	jer_param[9][6]=4.34831E-5;  jer_param_er[9][6]=1.98595E-5;
	jer_param[9][7]=-4.24538;  jer_param_er[9][7]=0.258383;
	jer_param[9][8]=0.235095;  jer_param_er[9][8]=0.0112121;
	jer_param[9][9]=-0.000176939;  jer_param_er[9][9]=7.11021E-5;
	jer_param[9][10]=17.7377;  jer_param_er[9][10]=6.10388;
	jer_param[9][11]=-6.60975;  jer_param_er[9][11]=1.55793;
	jer_param[9][12]=0.757554;  jer_param_er[9][12]=0.0890909;
	jer_param[10][0]=0.0702999;  jer_param_er[10][0]=0.00524409;
	jer_param[10][1]=-0.0001322;  jer_param_er[10][1]=1.44554E-5;
	jer_param[10][2]=-4.64174;  jer_param_er[10][2]=0.433494;
	jer_param[10][3]=1.06671;  jer_param_er[10][3]=0.0316977;
	jer_param[10][4]=0.00105544;  jer_param_er[10][4]=0.000168547;
	jer_param[10][5]=-0.02437;  jer_param_er[10][5]=0.00587413;
	jer_param[10][6]=2.26747E-5;  jer_param_er[10][6]=1.97589E-5;
	jer_param[10][7]=-3.8046;  jer_param_er[10][7]=0.382791;
	jer_param[10][8]=0.268162;  jer_param_er[10][8]=0.0148949;
	jer_param[10][9]=-6.26511E-5;  jer_param_er[10][9]=8.37528E-5;
	jer_param[10][10]=10.2852;  jer_param_er[10][10]=7.77639;
	jer_param[10][11]=-4.34279;  jer_param_er[10][11]=1.80774;
	jer_param[10][12]=0.564928;  jer_param_er[10][12]=0.0947709;
	jer_param[11][0]=0.0276176;  jer_param_er[11][0]=0.00700865;
	jer_param[11][1]=-1.08182E-5;  jer_param_er[11][1]=1.77066E-5;
	jer_param[11][2]=-4.30341;  jer_param_er[11][2]=0.642783;
	jer_param[11][3]=1.15381;  jer_param_er[11][3]=0.0466366;
	jer_param[11][4]=0.000660203;  jer_param_er[11][4]=0.000219101;
	jer_param[11][5]=-0.0365288;  jer_param_er[11][5]=0.007012;
	jer_param[11][6]=5.34815E-5;  jer_param_er[11][6]=2.21572E-5;
	jer_param[11][7]=-4.83064;  jer_param_er[11][7]=0.460712;
	jer_param[11][8]=0.28416;  jer_param_er[11][8]=0.0195007;
	jer_param[11][9]=-4.88533E-7;  jer_param_er[11][9]=0.000101491;
	jer_param[11][10]=14.821;  jer_param_er[11][10]=8.92697;
	jer_param[11][11]=-5.32396;  jer_param_er[11][11]=2.01108;
	jer_param[11][12]=0.568491;  jer_param_er[11][12]=0.10153;
	jer_param[12][0]=-0.00456382;  jer_param_er[12][0]=0.00797029;
	jer_param[12][1]=5.93327E-5;  jer_param_er[12][1]=2.2512E-5;
	jer_param[12][2]=-6.05263;  jer_param_er[12][2]=0.641518;
	jer_param[12][3]=1.01168;  jer_param_er[12][3]=0.0580327;
	jer_param[12][4]=0.00124107;  jer_param_er[12][4]=0.000266261;
	jer_param[12][5]=-0.0672398;  jer_param_er[12][5]=0.0103374;
	jer_param[12][6]=0.000144969;  jer_param_er[12][6]=3.13729E-5;
	jer_param[12][7]=-7.44458;  jer_param_er[12][7]=0.643695;
	jer_param[12][8]=0.236166;  jer_param_er[12][8]=0.0233694;
	jer_param[12][9]=0.000272798;  jer_param_er[12][9]=0.000115376;
	jer_param[12][10]=15.9149;  jer_param_er[12][10]=16.3677;
	jer_param[12][11]=-5.06773;  jer_param_er[12][11]=3.26237;
	jer_param[12][12]=0.515062;  jer_param_er[12][12]=0.14614;
	jer_param[13][0]=-0.0364159;  jer_param_er[13][0]=0.0230243;
	jer_param[13][1]=0.000103308;  jer_param_er[13][1]=6.3942E-5;
	jer_param[13][2]=-2.04414;  jer_param_er[13][2]=1.61387;
	jer_param[13][3]=1.16241;  jer_param_er[13][3]=0.154151;
	jer_param[13][4]=0.00449765;  jer_param_er[13][4]=0.000680983;
	jer_param[13][5]=-0.0713694;  jer_param_er[13][5]=0.0128;
	jer_param[13][6]=6.17634E-5;  jer_param_er[13][6]=3.50159E-5;
	jer_param[13][7]=-6.69936;  jer_param_er[13][7]=0.883717;
	jer_param[13][8]=0.401773;  jer_param_er[13][8]=0.0401819;
	jer_param[13][9]=0.00048766;  jer_param_er[13][9]=0.000171866;
	jer_param[13][10]=-6.84164;  jer_param_er[13][10]=15.4788;
	jer_param[13][11]=1.14003;  jer_param_er[13][11]=2.85757;
	jer_param[13][12]=0.063226;  jer_param_er[13][12]=0.116727;
	jer_param[14][0]=0;  jer_param_er[14][0]=0;
	jer_param[14][1]=0;  jer_param_er[14][1]=0;
	jer_param[14][2]=0;  jer_param_er[14][2]=0;
	jer_param[14][3]=0;  jer_param_er[14][3]=0;
	jer_param[14][4]=0;  jer_param_er[14][4]=0;
	jer_param[14][5]=-0.127208;  jer_param_er[14][5]=0.0191958;
	jer_param[14][6]=5.65717E-5;  jer_param_er[14][6]=3.57586E-5;
	jer_param[14][7]=-4.42422;  jer_param_er[14][7]=1.87829;
	jer_param[14][8]=0.9097;  jer_param_er[14][8]=0.0515384;
	jer_param[14][9]=-0.000149165;  jer_param_er[14][9]=0.000158924;
	jer_param[14][10]=0;  jer_param_er[14][10]=0;
	jer_param[14][11]=0;  jer_param_er[14][11]=0;
	jer_param[14][12]=0;  jer_param_er[14][12]=0;
*/
	
	//========================== OLD FIT FUNCTION PARAMETERS for E>60GeV only
  jer_param[0][0]=-0.0421892;  jer_param_er[0][0]=0.00582121;
  jer_param[0][1]=9.47639e-06;  jer_param_er[0][1]=1.40057e-05;
  jer_param[0][2]=8.07454;  jer_param_er[0][2]=0.465693;
  jer_param[0][3]=1.14135;  jer_param_er[0][3]=0.0254592;
  jer_param[0][4]=-5.71421e-05;  jer_param_er[0][4]=0.000123532;
  jer_param[0][5]=0.00673628;  jer_param_er[0][5]=0.00933377;
  jer_param[0][6]=0.000176218;  jer_param_er[0][6]=3.13972e-05;
  jer_param[0][7]=1.20457;  jer_param_er[0][7]=0.567125;
  jer_param[0][8]=0.15466;  jer_param_er[0][8]=0.0222841;
  jer_param[0][9]=0.00753183;  jer_param_er[0][9]=0.000195714;
  jer_param[0][10]=51.2133;  jer_param_er[0][10]=7.39744;
  jer_param[0][11]=-9.6038;  jer_param_er[0][11]=1.38217;
  jer_param[0][12]=0.531334;  jer_param_er[0][12]=0.061669;
  jer_param[1][0]=0.0868524;  jer_param_er[1][0]=0.00437809;
  jer_param[1][1]=-0.000184468;  jer_param_er[1][1]=1.29375e-05;
  jer_param[1][2]=0.549579;  jer_param_er[1][2]=0.309168;
  jer_param[1][3]=0.692304;  jer_param_er[1][3]=0.0157638;
  jer_param[1][4]=0.00295614;  jer_param_er[1][4]=0.000109943;
  jer_param[1][5]=0.013038;  jer_param_er[1][5]=0.00315384;
  jer_param[1][6]=-2.44657e-05;  jer_param_er[1][6]=8.59298e-06;
  jer_param[1][7]=-2.30164;  jer_param_er[1][7]=0.247848;
  jer_param[1][8]=0.104279;  jer_param_er[1][8]=0.00791589;
  jer_param[1][9]=0.00223899;  jer_param_er[1][9]=4.71293e-05;
  jer_param[1][10]=-16.5429;  jer_param_er[1][10]=5.24013;
  jer_param[1][11]=5.60005;  jer_param_er[1][11]=0.871842;
  jer_param[1][12]=-0.172535;  jer_param_er[1][12]=0.0345462;
  jer_param[2][0]=0.0821827;  jer_param_er[2][0]=0.00418742;
  jer_param[2][1]=-0.000203719;  jer_param_er[2][1]=1.2152e-05;
  jer_param[2][2]=0.592087;  jer_param_er[2][2]=0.2922;
  jer_param[2][3]=0.735988;  jer_param_er[2][3]=0.015049;
  jer_param[2][4]=0.00252903;  jer_param_er[2][4]=0.000106018;
  jer_param[2][5]=0.0128072;  jer_param_er[2][5]=0.00314417;
  jer_param[2][6]=-4.25571e-05;  jer_param_er[2][6]=8.44537e-06;
  jer_param[2][7]=-2.09644;  jer_param_er[2][7]=0.244244;
  jer_param[2][8]=0.132923;  jer_param_er[2][8]=0.00812173;
  jer_param[2][9]=0.00230918;  jer_param_er[2][9]=4.89787e-05;
  jer_param[2][10]=-8.33188;  jer_param_er[2][10]=5.50062;
  jer_param[2][11]=4.03066;  jer_param_er[2][11]=0.916829;
  jer_param[2][12]=-0.103124;  jer_param_er[2][12]=0.0362279;
  jer_param[3][0]=0.0939778;  jer_param_er[3][0]=0.00580377;
  jer_param[3][1]=-0.000188265;  jer_param_er[3][1]=1.65369e-05;
  jer_param[3][2]=0.418804;  jer_param_er[3][2]=0.41026;
  jer_param[3][3]=0.781879;  jer_param_er[3][3]=0.0200264;
  jer_param[3][4]=0.00313047;  jer_param_er[3][4]=0.00014161;
  jer_param[3][5]=0.0395259;  jer_param_er[3][5]=0.00377854;
  jer_param[3][6]=-8.32808e-05;  jer_param_er[3][6]=1.00488e-05;
  jer_param[3][7]=-3.24558;  jer_param_er[3][7]=0.291717;
  jer_param[3][8]=0.19488;  jer_param_er[3][8]=0.00993259;
  jer_param[3][9]=0.002342;  jer_param_er[3][9]=5.85965e-05;
  jer_param[3][10]=-11.3482;  jer_param_er[3][10]=5.70448;
  jer_param[3][11]=4.18195;  jer_param_er[3][11]=0.949661;
  jer_param[3][12]=-0.113168;  jer_param_er[3][12]=0.0373629;
  jer_param[4][0]=0.0769937;  jer_param_er[4][0]=0.00648534;
  jer_param[4][1]=-5.93526e-05;  jer_param_er[4][1]=1.80139e-05;
  jer_param[4][2]=2.73755;  jer_param_er[4][2]=0.462881;
  jer_param[4][3]=0.90762;  jer_param_er[4][3]=0.0268729;
  jer_param[4][4]=0.00569015;  jer_param_er[4][4]=0.000187464;
  jer_param[4][5]=0.0531015;  jer_param_er[4][5]=0.00499601;
  jer_param[4][6]=-0.000103791;  jer_param_er[4][6]=1.27635e-05;
  jer_param[4][7]=-3.51752;  jer_param_er[4][7]=0.382673;
  jer_param[4][8]=0.23916;  jer_param_er[4][8]=0.0148451;
  jer_param[4][9]=0.0033711;  jer_param_er[4][9]=8.53063e-05;
  jer_param[4][10]=-30.4048;  jer_param_er[4][10]=6.66339;
  jer_param[4][11]=7.51952;  jer_param_er[4][11]=1.09774;
  jer_param[4][12]=-0.212285;  jer_param_er[4][12]=0.042077;
  jer_param[5][0]=0.0973182;  jer_param_er[5][0]=0.00406679;
  jer_param[5][1]=-9.88762e-05;  jer_param_er[5][1]=1.04671e-05;
  jer_param[5][2]=-1.12998;  jer_param_er[5][2]=0.297544;
  jer_param[5][3]=0.648102;  jer_param_er[5][3]=0.0184293;
  jer_param[5][4]=0.00692417;  jer_param_er[5][4]=0.000126653;
  jer_param[5][5]=0.0588128;  jer_param_er[5][5]=0.00563585;
  jer_param[5][6]=-0.000138796;  jer_param_er[5][6]=1.38978e-05;
  jer_param[5][7]=-6.22123;  jer_param_er[5][7]=0.397175;
  jer_param[5][8]=0.168724;  jer_param_er[5][8]=0.0136361;
  jer_param[5][9]=0.00293732;  jer_param_er[5][9]=8.61362e-05;
  jer_param[5][10]=-98.751;  jer_param_er[5][10]=10.4943;
  jer_param[5][11]=18.773;  jer_param_er[5][11]=1.8084;
  jer_param[5][12]=-0.48381;  jer_param_er[5][12]=0.069637;
  jer_param[6][0]=0.167261;  jer_param_er[6][0]=0.00705154;
  jer_param[6][1]=-0.000437154;  jer_param_er[6][1]=1.84479e-05;
  jer_param[6][2]=-5.40154;  jer_param_er[6][2]=0.461716;
  jer_param[6][3]=0.963706;  jer_param_er[6][3]=0.0268102;
  jer_param[6][4]=0.00333274;  jer_param_er[6][4]=0.000255135;
  jer_param[6][5]=0.0369433;  jer_param_er[6][5]=0.00603786;
  jer_param[6][6]=-3.27303e-05;  jer_param_er[6][6]=1.87188e-05;
  jer_param[6][7]=-5.74418;  jer_param_er[6][7]=0.383447;
  jer_param[6][8]=0.0165817;  jer_param_er[6][8]=0.0124165;
  jer_param[6][9]=0.00463673;  jer_param_er[6][9]=0.000112166;
  jer_param[6][10]=2.46588;  jer_param_er[6][10]=5.8467;
  jer_param[6][11]=2.19141;  jer_param_er[6][11]=1.03661;
  jer_param[6][12]=-0.0807845;  jer_param_er[6][12]=0.0435705;
  jer_param[7][0]=0.0663829;  jer_param_er[7][0]=0.00595593;
  jer_param[7][1]=-0.000199426;  jer_param_er[7][1]=1.8693e-05;
  jer_param[7][2]=-1.05385;  jer_param_er[7][2]=0.397169;
  jer_param[7][3]=0.96766;  jer_param_er[7][3]=0.0196309;
  jer_param[7][4]=0.000956836;  jer_param_er[7][4]=0.00016031;
  jer_param[7][5]=0.0315858;  jer_param_er[7][5]=0.00790363;
  jer_param[7][6]=-5.74461e-05;  jer_param_er[7][6]=2.7121e-05;
  jer_param[7][7]=-6.65697;  jer_param_er[7][7]=0.497394;
  jer_param[7][8]=0.117688;  jer_param_er[7][8]=0.0121169;
  jer_param[7][9]=0.00158927;  jer_param_er[7][9]=0.000108159;
  jer_param[7][10]=-22.4305;  jer_param_er[7][10]=14.0853;
  jer_param[7][11]=5.74351;  jer_param_er[7][11]=2.63617;
  jer_param[7][12]=-0.102966;  jer_param_er[7][12]=0.119883;
  jer_param[8][0]=0.0958991;  jer_param_er[8][0]=0.00379364;
  jer_param[8][1]=-0.000174509;  jer_param_er[8][1]=9.92417e-06;
  jer_param[8][2]=-5.98343;  jer_param_er[8][2]=0.304702;
  jer_param[8][3]=0.875918;  jer_param_er[8][3]=0.0148333;
  jer_param[8][4]=0.000862512;  jer_param_er[8][4]=9.29965e-05;
  jer_param[8][5]=0.0161104;  jer_param_er[8][5]=0.00786032;
  jer_param[8][6]=-4.38698e-05;  jer_param_er[8][6]=2.0328e-05;
  jer_param[8][7]=-6.75709;  jer_param_er[8][7]=0.621678;
  jer_param[8][8]=0.181482;  jer_param_er[8][8]=0.0104139;
  jer_param[8][9]=7.30734e-05;  jer_param_er[8][9]=6.63151e-05;
  jer_param[8][10]=-40.9875;  jer_param_er[8][10]=30.2176;
  jer_param[8][11]=6.70382;  jer_param_er[8][11]=5.16006;
  jer_param[8][12]=0.158348;  jer_param_er[8][12]=0.214499;
  jer_param[9][0]=0.121962;  jer_param_er[9][0]=0.00836305;
  jer_param[9][1]=-0.000235837;  jer_param_er[9][1]=2.35369e-05;
  jer_param[9][2]=-8.88778;  jer_param_er[9][2]=0.662939;
  jer_param[9][3]=0.897486;  jer_param_er[9][3]=0.0236274;
  jer_param[9][4]=0.00143759;  jer_param_er[9][4]=0.000142724;
  jer_param[9][5]=0.0257533;  jer_param_er[9][5]=0.0141059;
  jer_param[9][6]=-6.64452e-05;  jer_param_er[9][6]=4.06557e-05;
  jer_param[9][7]=-7.58012;  jer_param_er[9][7]=1.10958;
  jer_param[9][8]=0.200961;  jer_param_er[9][8]=0.0156056;
  jer_param[9][9]=4.88112e-06;  jer_param_er[9][9]=9.3964e-05;
  jer_param[9][10]=-43.4458;  jer_param_er[9][10]=28.5424;
  jer_param[9][11]=6.24739;  jer_param_er[9][11]=5.11566;
  jer_param[9][12]=0.133833;  jer_param_er[9][12]=0.225589;
  jer_param[10][0]=0.101777;  jer_param_er[10][0]=0.00962828;
  jer_param[10][1]=-0.000203577;  jer_param_er[10][1]=2.33106e-05;
  jer_param[10][2]=-7.70192;  jer_param_er[10][2]=0.898148;
  jer_param[10][3]=1.09267;  jer_param_er[10][3]=0.0351777;
  jer_param[10][4]=0.000931535;  jer_param_er[10][4]=0.00018284;
  jer_param[10][5]=0.0122738;  jer_param_er[10][5]=0.0137083;
  jer_param[10][6]=-6.43889e-05;  jer_param_er[10][6]=3.508e-05;
  jer_param[10][7]=-7.20028;  jer_param_er[10][7]=1.23086;
  jer_param[10][8]=0.210035;  jer_param_er[10][8]=0.0196102;
  jer_param[10][9]=0.000219703;  jer_param_er[10][9]=0.000107601;
  jer_param[10][10]=-54.9674;  jer_param_er[10][10]=27.9653;
  jer_param[10][11]=6.99112;  jer_param_er[10][11]=4.82626;
  jer_param[10][12]=0.0923327;  jer_param_er[10][12]=0.204234;
  jer_param[11][0]=0.0613393;  jer_param_er[11][0]=0.00860218;
  jer_param[11][1]=-8.09778e-05;  jer_param_er[11][1]=2.0718e-05;
  jer_param[11][2]=-7.93629;  jer_param_er[11][2]=0.819867;
  jer_param[11][3]=1.08434;  jer_param_er[11][3]=0.0488564;
  jer_param[11][4]=0.000952762;  jer_param_er[11][4]=0.000229649;
  jer_param[11][5]=-0.0381364;  jer_param_er[11][5]=0.055397;
  jer_param[11][6]=8.24029e-05;  jer_param_er[11][6]=0.000132147;
  jer_param[11][7]=-5.8031;  jer_param_er[11][7]=5.40652;
  jer_param[11][8]=0.142164;  jer_param_er[11][8]=0.117338;
  jer_param[11][9]=0.000693611;  jer_param_er[11][9]=0.000571238;
  jer_param[11][10]=-19.0979;  jer_param_er[11][10]=25.4207;
  jer_param[11][11]=0.3485;  jer_param_er[11][11]=4.50045;
  jer_param[11][12]=0.341987;  jer_param_er[11][12]=0.191721;
  jer_param[12][0]=0.0248507;  jer_param_er[12][0]=0.0130463;
  jer_param[12][1]=-3.14293e-07;  jer_param_er[12][1]=3.05648e-05;
  jer_param[12][2]=-9.44294;  jer_param_er[12][2]=1.36865;
  jer_param[12][3]=1.12047;  jer_param_er[12][3]=0.0815018;
  jer_param[12][4]=0.000790941;  jer_param_er[12][4]=0.000354797;
  jer_param[12][5]=-0.057686;  jer_param_er[12][5]=0.0131712;
  jer_param[12][6]=0.000124849;  jer_param_er[12][6]=3.6104e-05;
  jer_param[12][7]=-8.49394;  jer_param_er[12][7]=1.06807;
  jer_param[12][8]=0.214141;  jer_param_er[12][8]=0.0341813;
  jer_param[12][9]=0.000361616;  jer_param_er[12][9]=0.000154593;
  jer_param[12][10]=-47.2874;  jer_param_er[12][10]=28.4302;
  jer_param[12][11]=6.08229;  jer_param_er[12][11]=5.20522;
  jer_param[12][12]=0.0640882;  jer_param_er[12][12]=0.218879;
  jer_param[13][0]=-0.0437888;  jer_param_er[13][0]=0.0247287;
  jer_param[13][1]=0.000117632;  jer_param_er[13][1]=6.63009e-05;
  jer_param[13][2]=-1.13873;  jer_param_er[13][2]=1.95754;
  jer_param[13][3]=1.26067;  jer_param_er[13][3]=0.192074;
  jer_param[13][4]=0.00413018;  jer_param_er[13][4]=0.000803293;
  jer_param[13][5]=-0.0676323;  jer_param_er[13][5]=0.0142384;
  jer_param[13][6]=5.46473e-05;  jer_param_er[13][6]=3.69747e-05;
  jer_param[13][7]=-7.1633;  jer_param_er[13][7]=1.17489;
  jer_param[13][8]=0.417182;  jer_param_er[13][8]=0.0520734;
  jer_param[13][9]=0.000431625;  jer_param_er[13][9]=0.000209578;
  jer_param[13][10]=-7.62579;  jer_param_er[13][10]=17.9923;
  jer_param[13][11]=1.26336;  jer_param_er[13][11]=3.20107;
  jer_param[13][12]=0.058685;  jer_param_er[13][12]=0.128244;
  jer_param[14][0]=0;  jer_param_er[14][0]=0;
  jer_param[14][1]=0;  jer_param_er[14][1]=0;
  jer_param[14][2]=0;  jer_param_er[14][2]=0;
  jer_param[14][3]=0;  jer_param_er[14][3]=0;
  jer_param[14][4]=0;  jer_param_er[14][4]=0;
  jer_param[14][5]=-0.0735557;  jer_param_er[14][5]=0.00925058;
  jer_param[14][6]=-3.92455e-05;  jer_param_er[14][6]=1.94005e-05;
  jer_param[14][7]=-10.0728;  jer_param_er[14][7]=0.626313;
  jer_param[14][8]=0.487154;  jer_param_er[14][8]=0.0249419;
  jer_param[14][9]=0.000992211;  jer_param_er[14][9]=0.000111403;
  jer_param[14][10]=0;  jer_param_er[14][10]=0;
  jer_param[14][11]=0;  jer_param_er[14][11]=0;
  jer_param[14][12]=0;  jer_param_er[14][12]=0;

	if(i<15 && j<13)
	{
		if(code==0) value=jer_param[i][j];
		if(code==1) value=jer_param_er[i][j];
	}
	return value;
}

// returns the exact bin values for JERs E<60GeV.
// >60GeV uses the fit values (from old fit function) 
double oldJERBinValues(const int iEtaBin, const int iPar, const int iBin)
{
	//====== Exact bin values for E<60GeV (zero bins filled with extrapolated
	//(nearest neighbor) values
	//======= par0-4 == Gaus Mean/Sigma, Landau MPV/Sigma and Relative fraction 

	if ( (iEtaBin >= 0 && iEtaBin < 15)
			&& (iPar >= 0 && iPar < 5)
			&& (iBin >= 0 && iBin<60) )
	{
		double vBinVals[15][5][60];		//etabin/parameter/energy bins(1GeV)
		vBinVals[0][0][0] = -0.163304;
		vBinVals[0][0][1] = -0.163304;
		vBinVals[0][0][2] = -0.163304;
		vBinVals[0][0][3] = -0.163304;
		vBinVals[0][0][4] = -0.163304;
		vBinVals[0][0][5] = -0.163304;
		vBinVals[0][0][6] = -0.163304;
		vBinVals[0][0][7] = 0.111493;
		vBinVals[0][0][8] = 0.111493;
		vBinVals[0][0][9] = 0.111493;
		vBinVals[0][0][10] = 0.111493;
		vBinVals[0][0][11] = 0.111493;
		vBinVals[0][0][12] = 0.313686;
		vBinVals[0][0][13] = 0.313686;
		vBinVals[0][0][14] = 0.313686;
		vBinVals[0][0][15] = 0.313686;
		vBinVals[0][0][16] = 0.313686;
		vBinVals[0][0][17] = 0.05472;
		vBinVals[0][0][18] = 0.05472;
		vBinVals[0][0][19] = 0.05472;
		vBinVals[0][0][20] = 0.05472;
		vBinVals[0][0][21] = 0.05472;
		vBinVals[0][0][22] = 0.0202207;
		vBinVals[0][0][23] = 0.0202207;
		vBinVals[0][0][24] = 0.0202207;
		vBinVals[0][0][25] = 0.0202207;
		vBinVals[0][0][26] = 0.0202207;
		vBinVals[0][0][27] = 0.12236;
		vBinVals[0][0][28] = 0.12236;
		vBinVals[0][0][29] = 0.12236;
		vBinVals[0][0][30] = 0.12236;
		vBinVals[0][0][31] = 0.12236;
		vBinVals[0][0][32] = 0.107328;
		vBinVals[0][0][33] = 0.107328;
		vBinVals[0][0][34] = 0.107328;
		vBinVals[0][0][35] = 0.107328;
		vBinVals[0][0][36] = 0.107328;
		vBinVals[0][0][37] = 0.128613;
		vBinVals[0][0][38] = 0.128613;
		vBinVals[0][0][39] = 0.128613;
		vBinVals[0][0][40] = 0.128613;
		vBinVals[0][0][41] = 0.128613;
		vBinVals[0][0][42] = 0.111842;
		vBinVals[0][0][43] = 0.111842;
		vBinVals[0][0][44] = 0.111842;
		vBinVals[0][0][45] = 0.111842;
		vBinVals[0][0][46] = 0.111842;
		vBinVals[0][0][47] = 0.127156;
		vBinVals[0][0][48] = 0.127156;
		vBinVals[0][0][49] = 0.127156;
		vBinVals[0][0][50] = 0.127156;
		vBinVals[0][0][51] = 0.127156;
		vBinVals[0][0][52] = 0.12639;
		vBinVals[0][0][53] = 0.12639;
		vBinVals[0][0][54] = 0.12639;
		vBinVals[0][0][55] = 0.12639;
		vBinVals[0][0][56] = 0.12639;
		vBinVals[0][0][57] = 0.107631;
		vBinVals[0][0][58] = 0.107631;
		vBinVals[0][0][59] = 0.107631;
		vBinVals[0][1][0] = 0.18161;
		vBinVals[0][1][1] = 0.18161;
		vBinVals[0][1][2] = 0.18161;
		vBinVals[0][1][3] = 0.18161;
		vBinVals[0][1][4] = 0.18161;
		vBinVals[0][1][5] = 0.18161;
		vBinVals[0][1][6] = 0.18161;
		vBinVals[0][1][7] = 0.246076;
		vBinVals[0][1][8] = 0.246076;
		vBinVals[0][1][9] = 0.246076;
		vBinVals[0][1][10] = 0.246076;
		vBinVals[0][1][11] = 0.246076;
		vBinVals[0][1][12] = 0.32651;
		vBinVals[0][1][13] = 0.32651;
		vBinVals[0][1][14] = 0.32651;
		vBinVals[0][1][15] = 0.32651;
		vBinVals[0][1][16] = 0.32651;
		vBinVals[0][1][17] = 0.197646;
		vBinVals[0][1][18] = 0.197646;
		vBinVals[0][1][19] = 0.197646;
		vBinVals[0][1][20] = 0.197646;
		vBinVals[0][1][21] = 0.197646;
		vBinVals[0][1][22] = 0.162135;
		vBinVals[0][1][23] = 0.162135;
		vBinVals[0][1][24] = 0.162135;
		vBinVals[0][1][25] = 0.162135;
		vBinVals[0][1][26] = 0.162135;
		vBinVals[0][1][27] = 0.182794;
		vBinVals[0][1][28] = 0.182794;
		vBinVals[0][1][29] = 0.182794;
		vBinVals[0][1][30] = 0.182794;
		vBinVals[0][1][31] = 0.182794;
		vBinVals[0][1][32] = 0.162474;
		vBinVals[0][1][33] = 0.162474;
		vBinVals[0][1][34] = 0.162474;
		vBinVals[0][1][35] = 0.162474;
		vBinVals[0][1][36] = 0.162474;
		vBinVals[0][1][37] = 0.160287;
		vBinVals[0][1][38] = 0.160287;
		vBinVals[0][1][39] = 0.160287;
		vBinVals[0][1][40] = 0.160287;
		vBinVals[0][1][41] = 0.160287;
		vBinVals[0][1][42] = 0.155634;
		vBinVals[0][1][43] = 0.155634;
		vBinVals[0][1][44] = 0.155634;
		vBinVals[0][1][45] = 0.155634;
		vBinVals[0][1][46] = 0.155634;
		vBinVals[0][1][47] = 0.152787;
		vBinVals[0][1][48] = 0.152787;
		vBinVals[0][1][49] = 0.152787;
		vBinVals[0][1][50] = 0.152787;
		vBinVals[0][1][51] = 0.152787;
		vBinVals[0][1][52] = 0.148066;
		vBinVals[0][1][53] = 0.148066;
		vBinVals[0][1][54] = 0.148066;
		vBinVals[0][1][55] = 0.148066;
		vBinVals[0][1][56] = 0.148066;
		vBinVals[0][1][57] = 0.135741;
		vBinVals[0][1][58] = 0.135741;
		vBinVals[0][1][59] = 0.135741;
		vBinVals[0][2][0] = -0.306665;
		vBinVals[0][2][1] = -0.306665;
		vBinVals[0][2][2] = -0.306665;
		vBinVals[0][2][3] = -0.306665;
		vBinVals[0][2][4] = -0.306665;
		vBinVals[0][2][5] = -0.306665;
		vBinVals[0][2][6] = -0.306665;
		vBinVals[0][2][7] = -0.194985;
		vBinVals[0][2][8] = -0.194985;
		vBinVals[0][2][9] = -0.194985;
		vBinVals[0][2][10] = -0.194985;
		vBinVals[0][2][11] = -0.194985;
		vBinVals[0][2][12] = -0.0608578;
		vBinVals[0][2][13] = -0.0608578;
		vBinVals[0][2][14] = -0.0608578;
		vBinVals[0][2][15] = -0.0608578;
		vBinVals[0][2][16] = -0.0608578;
		vBinVals[0][2][17] = -0.0185246;
		vBinVals[0][2][18] = -0.0185246;
		vBinVals[0][2][19] = -0.0185246;
		vBinVals[0][2][20] = -0.0185246;
		vBinVals[0][2][21] = -0.0185246;
		vBinVals[0][2][22] = 0.0162521;
		vBinVals[0][2][23] = 0.0162521;
		vBinVals[0][2][24] = 0.0162521;
		vBinVals[0][2][25] = 0.0162521;
		vBinVals[0][2][26] = 0.0162521;
		vBinVals[0][2][27] = -0.0177118;
		vBinVals[0][2][28] = -0.0177118;
		vBinVals[0][2][29] = -0.0177118;
		vBinVals[0][2][30] = -0.0177118;
		vBinVals[0][2][31] = -0.0177118;
		vBinVals[0][2][32] = 0.00101698;
		vBinVals[0][2][33] = 0.00101698;
		vBinVals[0][2][34] = 0.00101698;
		vBinVals[0][2][35] = 0.00101698;
		vBinVals[0][2][36] = 0.00101698;
		vBinVals[0][2][37] = 0.0108564;
		vBinVals[0][2][38] = 0.0108564;
		vBinVals[0][2][39] = 0.0108564;
		vBinVals[0][2][40] = 0.0108564;
		vBinVals[0][2][41] = 0.0108564;
		vBinVals[0][2][42] = 0.0256493;
		vBinVals[0][2][43] = 0.0256493;
		vBinVals[0][2][44] = 0.0256493;
		vBinVals[0][2][45] = 0.0256493;
		vBinVals[0][2][46] = 0.0256493;
		vBinVals[0][2][47] = 0.0170728;
		vBinVals[0][2][48] = 0.0170728;
		vBinVals[0][2][49] = 0.0170728;
		vBinVals[0][2][50] = 0.0170728;
		vBinVals[0][2][51] = 0.0170728;
		vBinVals[0][2][52] = 0.0227435;
		vBinVals[0][2][53] = 0.0227435;
		vBinVals[0][2][54] = 0.0227435;
		vBinVals[0][2][55] = 0.0227435;
		vBinVals[0][2][56] = 0.0227435;
		vBinVals[0][2][57] = 0.0284676;
		vBinVals[0][2][58] = 0.0284676;
		vBinVals[0][2][59] = 0.0284676;
		vBinVals[0][3][0] = 0.0953565;
		vBinVals[0][3][1] = 0.0953565;
		vBinVals[0][3][2] = 0.0953565;
		vBinVals[0][3][3] = 0.0953565;
		vBinVals[0][3][4] = 0.0953565;
		vBinVals[0][3][5] = 0.0953565;
		vBinVals[0][3][6] = 0.0953565;
		vBinVals[0][3][7] = 0.121164;
		vBinVals[0][3][8] = 0.121164;
		vBinVals[0][3][9] = 0.121164;
		vBinVals[0][3][10] = 0.121164;
		vBinVals[0][3][11] = 0.121164;
		vBinVals[0][3][12] = 0.143362;
		vBinVals[0][3][13] = 0.143362;
		vBinVals[0][3][14] = 0.143362;
		vBinVals[0][3][15] = 0.143362;
		vBinVals[0][3][16] = 0.143362;
		vBinVals[0][3][17] = 0.141051;
		vBinVals[0][3][18] = 0.141051;
		vBinVals[0][3][19] = 0.141051;
		vBinVals[0][3][20] = 0.141051;
		vBinVals[0][3][21] = 0.141051;
		vBinVals[0][3][22] = 0.133079;
		vBinVals[0][3][23] = 0.133079;
		vBinVals[0][3][24] = 0.133079;
		vBinVals[0][3][25] = 0.133079;
		vBinVals[0][3][26] = 0.133079;
		vBinVals[0][3][27] = 0.110676;
		vBinVals[0][3][28] = 0.110676;
		vBinVals[0][3][29] = 0.110676;
		vBinVals[0][3][30] = 0.110676;
		vBinVals[0][3][31] = 0.110676;
		vBinVals[0][3][32] = 0.108428;
		vBinVals[0][3][33] = 0.108428;
		vBinVals[0][3][34] = 0.108428;
		vBinVals[0][3][35] = 0.108428;
		vBinVals[0][3][36] = 0.108428;
		vBinVals[0][3][37] = 0.105472;
		vBinVals[0][3][38] = 0.105472;
		vBinVals[0][3][39] = 0.105472;
		vBinVals[0][3][40] = 0.105472;
		vBinVals[0][3][41] = 0.105472;
		vBinVals[0][3][42] = 0.106455;
		vBinVals[0][3][43] = 0.106455;
		vBinVals[0][3][44] = 0.106455;
		vBinVals[0][3][45] = 0.106455;
		vBinVals[0][3][46] = 0.106455;
		vBinVals[0][3][47] = 0.100473;
		vBinVals[0][3][48] = 0.100473;
		vBinVals[0][3][49] = 0.100473;
		vBinVals[0][3][50] = 0.100473;
		vBinVals[0][3][51] = 0.100473;
		vBinVals[0][3][52] = 0.100282;
		vBinVals[0][3][53] = 0.100282;
		vBinVals[0][3][54] = 0.100282;
		vBinVals[0][3][55] = 0.100282;
		vBinVals[0][3][56] = 0.100282;
		vBinVals[0][3][57] = 0.0978895;
		vBinVals[0][3][58] = 0.0978895;
		vBinVals[0][3][59] = 0.0978895;
		vBinVals[0][4][0] = 0.276402;
		vBinVals[0][4][1] = 0.276402;
		vBinVals[0][4][2] = 0.276402;
		vBinVals[0][4][3] = 0.276402;
		vBinVals[0][4][4] = 0.276402;
		vBinVals[0][4][5] = 0.276402;
		vBinVals[0][4][6] = 0.276402;
		vBinVals[0][4][7] = 0.188836;
		vBinVals[0][4][8] = 0.188836;
		vBinVals[0][4][9] = 0.188836;
		vBinVals[0][4][10] = 0.188836;
		vBinVals[0][4][11] = 0.188836;
		vBinVals[0][4][12] = 0.0474377;
		vBinVals[0][4][13] = 0.0474377;
		vBinVals[0][4][14] = 0.0474377;
		vBinVals[0][4][15] = 0.0474377;
		vBinVals[0][4][16] = 0.0474377;
		vBinVals[0][4][17] = 0.0317284;
		vBinVals[0][4][18] = 0.0317284;
		vBinVals[0][4][19] = 0.0317284;
		vBinVals[0][4][20] = 0.0317284;
		vBinVals[0][4][21] = 0.0317284;
		vBinVals[0][4][22] = 0.0776162;
		vBinVals[0][4][23] = 0.0776162;
		vBinVals[0][4][24] = 0.0776162;
		vBinVals[0][4][25] = 0.0776162;
		vBinVals[0][4][26] = 0.0776162;
		vBinVals[0][4][27] = 0.123312;
		vBinVals[0][4][28] = 0.123312;
		vBinVals[0][4][29] = 0.123312;
		vBinVals[0][4][30] = 0.123312;
		vBinVals[0][4][31] = 0.123312;
		vBinVals[0][4][32] = 0.109676;
		vBinVals[0][4][33] = 0.109676;
		vBinVals[0][4][34] = 0.109676;
		vBinVals[0][4][35] = 0.109676;
		vBinVals[0][4][36] = 0.109676;
		vBinVals[0][4][37] = 0.123622;
		vBinVals[0][4][38] = 0.123622;
		vBinVals[0][4][39] = 0.123622;
		vBinVals[0][4][40] = 0.123622;
		vBinVals[0][4][41] = 0.123622;
		vBinVals[0][4][42] = 0.133726;
		vBinVals[0][4][43] = 0.133726;
		vBinVals[0][4][44] = 0.133726;
		vBinVals[0][4][45] = 0.133726;
		vBinVals[0][4][46] = 0.133726;
		vBinVals[0][4][47] = 0.13902;
		vBinVals[0][4][48] = 0.13902;
		vBinVals[0][4][49] = 0.13902;
		vBinVals[0][4][50] = 0.13902;
		vBinVals[0][4][51] = 0.13902;
		vBinVals[0][4][52] = 0.126203;
		vBinVals[0][4][53] = 0.126203;
		vBinVals[0][4][54] = 0.126203;
		vBinVals[0][4][55] = 0.126203;
		vBinVals[0][4][56] = 0.126203;
		vBinVals[0][4][57] = 0.160674;
		vBinVals[0][4][58] = 0.160674;
		vBinVals[0][4][59] = 0.160674;
		vBinVals[1][0][0] = -0.173962;
		vBinVals[1][0][1] = -0.173962;
		vBinVals[1][0][2] = -0.173962;
		vBinVals[1][0][3] = -0.173962;
		vBinVals[1][0][4] = -0.173962;
		vBinVals[1][0][5] = -0.173962;
		vBinVals[1][0][6] = -0.173962;
		vBinVals[1][0][7] = 0.0696978;
		vBinVals[1][0][8] = 0.0696978;
		vBinVals[1][0][9] = 0.0696978;
		vBinVals[1][0][10] = 0.0696978;
		vBinVals[1][0][11] = 0.0696978;
		vBinVals[1][0][12] = 0.134907;
		vBinVals[1][0][13] = 0.134907;
		vBinVals[1][0][14] = 0.134907;
		vBinVals[1][0][15] = 0.134907;
		vBinVals[1][0][16] = 0.134907;
		vBinVals[1][0][17] = 0.0270979;
		vBinVals[1][0][18] = 0.0270979;
		vBinVals[1][0][19] = 0.0270979;
		vBinVals[1][0][20] = 0.0270979;
		vBinVals[1][0][21] = 0.0270979;
		vBinVals[1][0][22] = 0.0710903;
		vBinVals[1][0][23] = 0.0710903;
		vBinVals[1][0][24] = 0.0710903;
		vBinVals[1][0][25] = 0.0710903;
		vBinVals[1][0][26] = 0.0710903;
		vBinVals[1][0][27] = 0.0796466;
		vBinVals[1][0][28] = 0.0796466;
		vBinVals[1][0][29] = 0.0796466;
		vBinVals[1][0][30] = 0.0796466;
		vBinVals[1][0][31] = 0.0796466;
		vBinVals[1][0][32] = 0.0796466;
		vBinVals[1][0][33] = 0.0808828;
		vBinVals[1][0][34] = 0.0808828;
		vBinVals[1][0][35] = 0.0808828;
		vBinVals[1][0][36] = 0.0808828;
		vBinVals[1][0][37] = 0.097228;
		vBinVals[1][0][38] = 0.097228;
		vBinVals[1][0][39] = 0.097228;
		vBinVals[1][0][40] = 0.097228;
		vBinVals[1][0][41] = 0.097228;
		vBinVals[1][0][42] = 0.0940544;
		vBinVals[1][0][43] = 0.0940544;
		vBinVals[1][0][44] = 0.0940544;
		vBinVals[1][0][45] = 0.0940544;
		vBinVals[1][0][46] = 0.0940544;
		vBinVals[1][0][47] = 0.0975833;
		vBinVals[1][0][48] = 0.0975833;
		vBinVals[1][0][49] = 0.0975833;
		vBinVals[1][0][50] = 0.0975833;
		vBinVals[1][0][51] = 0.0975833;
		vBinVals[1][0][52] = 0.0842198;
		vBinVals[1][0][53] = 0.0842198;
		vBinVals[1][0][54] = 0.0842198;
		vBinVals[1][0][55] = 0.0842198;
		vBinVals[1][0][56] = 0.0842198;
		vBinVals[1][0][57] = 0.0936287;
		vBinVals[1][0][58] = 0.0936287;
		vBinVals[1][0][59] = 0.0936287;
		vBinVals[1][1][0] = 0.170472;
		vBinVals[1][1][1] = 0.170472;
		vBinVals[1][1][2] = 0.170472;
		vBinVals[1][1][3] = 0.170472;
		vBinVals[1][1][4] = 0.170472;
		vBinVals[1][1][5] = 0.170472;
		vBinVals[1][1][6] = 0.170472;
		vBinVals[1][1][7] = 0.226716;
		vBinVals[1][1][8] = 0.226716;
		vBinVals[1][1][9] = 0.226716;
		vBinVals[1][1][10] = 0.226716;
		vBinVals[1][1][11] = 0.226716;
		vBinVals[1][1][12] = 0.245448;
		vBinVals[1][1][13] = 0.245448;
		vBinVals[1][1][14] = 0.245448;
		vBinVals[1][1][15] = 0.245448;
		vBinVals[1][1][16] = 0.245448;
		vBinVals[1][1][17] = 0.175182;
		vBinVals[1][1][18] = 0.175182;
		vBinVals[1][1][19] = 0.175182;
		vBinVals[1][1][20] = 0.175182;
		vBinVals[1][1][21] = 0.175182;
		vBinVals[1][1][22] = 0.173569;
		vBinVals[1][1][23] = 0.173569;
		vBinVals[1][1][24] = 0.173569;
		vBinVals[1][1][25] = 0.173569;
		vBinVals[1][1][26] = 0.173569;
		vBinVals[1][1][27] = 0.161311;
		vBinVals[1][1][28] = 0.161311;
		vBinVals[1][1][29] = 0.161311;
		vBinVals[1][1][30] = 0.161311;
		vBinVals[1][1][31] = 0.161311;
		vBinVals[1][1][32] = 0.161311;
		vBinVals[1][1][33] = 0.140974;
		vBinVals[1][1][34] = 0.140974;
		vBinVals[1][1][35] = 0.140974;
		vBinVals[1][1][36] = 0.140974;
		vBinVals[1][1][37] = 0.136968;
		vBinVals[1][1][38] = 0.136968;
		vBinVals[1][1][39] = 0.136968;
		vBinVals[1][1][40] = 0.136968;
		vBinVals[1][1][41] = 0.136968;
		vBinVals[1][1][42] = 0.130826;
		vBinVals[1][1][43] = 0.130826;
		vBinVals[1][1][44] = 0.130826;
		vBinVals[1][1][45] = 0.130826;
		vBinVals[1][1][46] = 0.130826;
		vBinVals[1][1][47] = 0.130119;
		vBinVals[1][1][48] = 0.130119;
		vBinVals[1][1][49] = 0.130119;
		vBinVals[1][1][50] = 0.130119;
		vBinVals[1][1][51] = 0.130119;
		vBinVals[1][1][52] = 0.124709;
		vBinVals[1][1][53] = 0.124709;
		vBinVals[1][1][54] = 0.124709;
		vBinVals[1][1][55] = 0.124709;
		vBinVals[1][1][56] = 0.124709;
		vBinVals[1][1][57] = 0.123287;
		vBinVals[1][1][58] = 0.123287;
		vBinVals[1][1][59] = 0.123287;
		vBinVals[1][2][0] = -0.317176;
		vBinVals[1][2][1] = -0.317176;
		vBinVals[1][2][2] = -0.317176;
		vBinVals[1][2][3] = -0.317176;
		vBinVals[1][2][4] = -0.317176;
		vBinVals[1][2][5] = -0.317176;
		vBinVals[1][2][6] = -0.317176;
		vBinVals[1][2][7] = -0.202306;
		vBinVals[1][2][8] = -0.202306;
		vBinVals[1][2][9] = -0.202306;
		vBinVals[1][2][10] = -0.202306;
		vBinVals[1][2][11] = -0.202306;
		vBinVals[1][2][12] = -0.0851608;
		vBinVals[1][2][13] = -0.0851608;
		vBinVals[1][2][14] = -0.0851608;
		vBinVals[1][2][15] = -0.0851608;
		vBinVals[1][2][16] = -0.0851608;
		vBinVals[1][2][17] = -0.0575838;
		vBinVals[1][2][18] = -0.0575838;
		vBinVals[1][2][19] = -0.0575838;
		vBinVals[1][2][20] = -0.0575838;
		vBinVals[1][2][21] = -0.0575838;
		vBinVals[1][2][22] = -0.056881;
		vBinVals[1][2][23] = -0.056881;
		vBinVals[1][2][24] = -0.056881;
		vBinVals[1][2][25] = -0.056881;
		vBinVals[1][2][26] = -0.056881;
		vBinVals[1][2][27] = -0.0458422;
		vBinVals[1][2][28] = -0.0458422;
		vBinVals[1][2][29] = -0.0458422;
		vBinVals[1][2][30] = -0.0458422;
		vBinVals[1][2][31] = -0.0458422;
		vBinVals[1][2][32] = -0.0458422;
		vBinVals[1][2][33] = -0.0505789;
		vBinVals[1][2][34] = -0.0505789;
		vBinVals[1][2][35] = -0.0505789;
		vBinVals[1][2][36] = -0.0505789;
		vBinVals[1][2][37] = -0.0351266;
		vBinVals[1][2][38] = -0.0351266;
		vBinVals[1][2][39] = -0.0351266;
		vBinVals[1][2][40] = -0.0351266;
		vBinVals[1][2][41] = -0.0351266;
		vBinVals[1][2][42] = -0.0294399;
		vBinVals[1][2][43] = -0.0294399;
		vBinVals[1][2][44] = -0.0294399;
		vBinVals[1][2][45] = -0.0294399;
		vBinVals[1][2][46] = -0.0294399;
		vBinVals[1][2][47] = -0.0214843;
		vBinVals[1][2][48] = -0.0214843;
		vBinVals[1][2][49] = -0.0214843;
		vBinVals[1][2][50] = -0.0214843;
		vBinVals[1][2][51] = -0.0214843;
		vBinVals[1][2][52] = -0.000856292;
		vBinVals[1][2][53] = -0.000856292;
		vBinVals[1][2][54] = -0.000856292;
		vBinVals[1][2][55] = -0.000856292;
		vBinVals[1][2][56] = -0.000856292;
		vBinVals[1][2][57] = -0.0176067;
		vBinVals[1][2][58] = -0.0176067;
		vBinVals[1][2][59] = -0.0176067;
		vBinVals[1][3][0] = 0.0859657;
		vBinVals[1][3][1] = 0.0859657;
		vBinVals[1][3][2] = 0.0859657;
		vBinVals[1][3][3] = 0.0859657;
		vBinVals[1][3][4] = 0.0859657;
		vBinVals[1][3][5] = 0.0859657;
		vBinVals[1][3][6] = 0.0859657;
		vBinVals[1][3][7] = 0.117386;
		vBinVals[1][3][8] = 0.117386;
		vBinVals[1][3][9] = 0.117386;
		vBinVals[1][3][10] = 0.117386;
		vBinVals[1][3][11] = 0.117386;
		vBinVals[1][3][12] = 0.129728;
		vBinVals[1][3][13] = 0.129728;
		vBinVals[1][3][14] = 0.129728;
		vBinVals[1][3][15] = 0.129728;
		vBinVals[1][3][16] = 0.129728;
		vBinVals[1][3][17] = 0.11665;
		vBinVals[1][3][18] = 0.11665;
		vBinVals[1][3][19] = 0.11665;
		vBinVals[1][3][20] = 0.11665;
		vBinVals[1][3][21] = 0.11665;
		vBinVals[1][3][22] = 0.0935906;
		vBinVals[1][3][23] = 0.0935906;
		vBinVals[1][3][24] = 0.0935906;
		vBinVals[1][3][25] = 0.0935906;
		vBinVals[1][3][26] = 0.0935906;
		vBinVals[1][3][27] = 0.0847934;
		vBinVals[1][3][28] = 0.0847934;
		vBinVals[1][3][29] = 0.0847934;
		vBinVals[1][3][30] = 0.0847934;
		vBinVals[1][3][31] = 0.0847934;
		vBinVals[1][3][32] = 0.0847934;
		vBinVals[1][3][33] = 0.0757152;
		vBinVals[1][3][34] = 0.0757152;
		vBinVals[1][3][35] = 0.0757152;
		vBinVals[1][3][36] = 0.0757152;
		vBinVals[1][3][37] = 0.0764909;
		vBinVals[1][3][38] = 0.0764909;
		vBinVals[1][3][39] = 0.0764909;
		vBinVals[1][3][40] = 0.0764909;
		vBinVals[1][3][41] = 0.0764909;
		vBinVals[1][3][42] = 0.072766;
		vBinVals[1][3][43] = 0.072766;
		vBinVals[1][3][44] = 0.072766;
		vBinVals[1][3][45] = 0.072766;
		vBinVals[1][3][46] = 0.072766;
		vBinVals[1][3][47] = 0.0669271;
		vBinVals[1][3][48] = 0.0669271;
		vBinVals[1][3][49] = 0.0669271;
		vBinVals[1][3][50] = 0.0669271;
		vBinVals[1][3][51] = 0.0669271;
		vBinVals[1][3][52] = 0.0736983;
		vBinVals[1][3][53] = 0.0736983;
		vBinVals[1][3][54] = 0.0736983;
		vBinVals[1][3][55] = 0.0736983;
		vBinVals[1][3][56] = 0.0736983;
		vBinVals[1][3][57] = 0.061058;
		vBinVals[1][3][58] = 0.061058;
		vBinVals[1][3][59] = 0.061058;
		vBinVals[1][4][0] = 0.316622;
		vBinVals[1][4][1] = 0.316622;
		vBinVals[1][4][2] = 0.316622;
		vBinVals[1][4][3] = 0.316622;
		vBinVals[1][4][4] = 0.316622;
		vBinVals[1][4][5] = 0.316622;
		vBinVals[1][4][6] = 0.316622;
		vBinVals[1][4][7] = 0.177416;
		vBinVals[1][4][8] = 0.177416;
		vBinVals[1][4][9] = 0.177416;
		vBinVals[1][4][10] = 0.177416;
		vBinVals[1][4][11] = 0.177416;
		vBinVals[1][4][12] = 0.0780594;
		vBinVals[1][4][13] = 0.0780594;
		vBinVals[1][4][14] = 0.0780594;
		vBinVals[1][4][15] = 0.0780594;
		vBinVals[1][4][16] = 0.0780594;
		vBinVals[1][4][17] = 0.10498;
		vBinVals[1][4][18] = 0.10498;
		vBinVals[1][4][19] = 0.10498;
		vBinVals[1][4][20] = 0.10498;
		vBinVals[1][4][21] = 0.10498;
		vBinVals[1][4][22] = 0.165017;
		vBinVals[1][4][23] = 0.165017;
		vBinVals[1][4][24] = 0.165017;
		vBinVals[1][4][25] = 0.165017;
		vBinVals[1][4][26] = 0.165017;
		vBinVals[1][4][27] = 0.155853;
		vBinVals[1][4][28] = 0.155853;
		vBinVals[1][4][29] = 0.155853;
		vBinVals[1][4][30] = 0.155853;
		vBinVals[1][4][31] = 0.155853;
		vBinVals[1][4][32] = 0.155853;
		vBinVals[1][4][33] = 0.263704;
		vBinVals[1][4][34] = 0.263704;
		vBinVals[1][4][35] = 0.263704;
		vBinVals[1][4][36] = 0.263704;
		vBinVals[1][4][37] = 0.261954;
		vBinVals[1][4][38] = 0.261954;
		vBinVals[1][4][39] = 0.261954;
		vBinVals[1][4][40] = 0.261954;
		vBinVals[1][4][41] = 0.261954;
		vBinVals[1][4][42] = 0.282686;
		vBinVals[1][4][43] = 0.282686;
		vBinVals[1][4][44] = 0.282686;
		vBinVals[1][4][45] = 0.282686;
		vBinVals[1][4][46] = 0.282686;
		vBinVals[1][4][47] = 0.242429;
		vBinVals[1][4][48] = 0.242429;
		vBinVals[1][4][49] = 0.242429;
		vBinVals[1][4][50] = 0.242429;
		vBinVals[1][4][51] = 0.242429;
		vBinVals[1][4][52] = 0.264995;
		vBinVals[1][4][53] = 0.264995;
		vBinVals[1][4][54] = 0.264995;
		vBinVals[1][4][55] = 0.264995;
		vBinVals[1][4][56] = 0.264995;
		vBinVals[1][4][57] = 0.234468;
		vBinVals[1][4][58] = 0.234468;
		vBinVals[1][4][59] = 0.234468;
		vBinVals[2][0][0] = -0.274551;
		vBinVals[2][0][1] = -0.274551;
		vBinVals[2][0][2] = -0.274551;
		vBinVals[2][0][3] = -0.274551;
		vBinVals[2][0][4] = -0.274551;
		vBinVals[2][0][5] = -0.274551;
		vBinVals[2][0][6] = -0.274551;
		vBinVals[2][0][7] = 0.0260218;
		vBinVals[2][0][8] = 0.0260218;
		vBinVals[2][0][9] = 0.0260218;
		vBinVals[2][0][10] = 0.0260218;
		vBinVals[2][0][11] = 0.0260218;
		vBinVals[2][0][12] = 0.134474;
		vBinVals[2][0][13] = 0.134474;
		vBinVals[2][0][14] = 0.134474;
		vBinVals[2][0][15] = 0.134474;
		vBinVals[2][0][16] = 0.134474;
		vBinVals[2][0][17] = 0.0274719;
		vBinVals[2][0][18] = 0.0274719;
		vBinVals[2][0][19] = 0.0274719;
		vBinVals[2][0][20] = 0.0274719;
		vBinVals[2][0][21] = 0.0274719;
		vBinVals[2][0][22] = 0.0490225;
		vBinVals[2][0][23] = 0.0490225;
		vBinVals[2][0][24] = 0.0490225;
		vBinVals[2][0][25] = 0.0490225;
		vBinVals[2][0][26] = 0.0490225;
		vBinVals[2][0][27] = 0.0692719;
		vBinVals[2][0][28] = 0.0692719;
		vBinVals[2][0][29] = 0.0692719;
		vBinVals[2][0][30] = 0.0692719;
		vBinVals[2][0][31] = 0.0692719;
		vBinVals[2][0][32] = 0.0692719;
		vBinVals[2][0][33] = 0.0508256;
		vBinVals[2][0][34] = 0.0508256;
		vBinVals[2][0][35] = 0.0508256;
		vBinVals[2][0][36] = 0.0508256;
		vBinVals[2][0][37] = 0.0819087;
		vBinVals[2][0][38] = 0.0819087;
		vBinVals[2][0][39] = 0.0819087;
		vBinVals[2][0][40] = 0.0819087;
		vBinVals[2][0][41] = 0.0819087;
		vBinVals[2][0][42] = 0.0888502;
		vBinVals[2][0][43] = 0.0888502;
		vBinVals[2][0][44] = 0.0888502;
		vBinVals[2][0][45] = 0.0888502;
		vBinVals[2][0][46] = 0.0888502;
		vBinVals[2][0][47] = 0.0875144;
		vBinVals[2][0][48] = 0.0875144;
		vBinVals[2][0][49] = 0.0875144;
		vBinVals[2][0][50] = 0.0875144;
		vBinVals[2][0][51] = 0.0875144;
		vBinVals[2][0][52] = 0.0900286;
		vBinVals[2][0][53] = 0.0900286;
		vBinVals[2][0][54] = 0.0900286;
		vBinVals[2][0][55] = 0.0900286;
		vBinVals[2][0][56] = 0.0900286;
		vBinVals[2][0][57] = 0.0731152;
		vBinVals[2][0][58] = 0.0731152;
		vBinVals[2][0][59] = 0.0731152;
		vBinVals[2][1][0] = -0.126545;
		vBinVals[2][1][1] = -0.126545;
		vBinVals[2][1][2] = -0.126545;
		vBinVals[2][1][3] = -0.126545;
		vBinVals[2][1][4] = -0.126545;
		vBinVals[2][1][5] = -0.126545;
		vBinVals[2][1][6] = -0.126545;
		vBinVals[2][1][7] = 0.225621;
		vBinVals[2][1][8] = 0.225621;
		vBinVals[2][1][9] = 0.225621;
		vBinVals[2][1][10] = 0.225621;
		vBinVals[2][1][11] = 0.225621;
		vBinVals[2][1][12] = 0.240498;
		vBinVals[2][1][13] = 0.240498;
		vBinVals[2][1][14] = 0.240498;
		vBinVals[2][1][15] = 0.240498;
		vBinVals[2][1][16] = 0.240498;
		vBinVals[2][1][17] = 0.174617;
		vBinVals[2][1][18] = 0.174617;
		vBinVals[2][1][19] = 0.174617;
		vBinVals[2][1][20] = 0.174617;
		vBinVals[2][1][21] = 0.174617;
		vBinVals[2][1][22] = 0.162862;
		vBinVals[2][1][23] = 0.162862;
		vBinVals[2][1][24] = 0.162862;
		vBinVals[2][1][25] = 0.162862;
		vBinVals[2][1][26] = 0.162862;
		vBinVals[2][1][27] = 0.160188;
		vBinVals[2][1][28] = 0.160188;
		vBinVals[2][1][29] = 0.160188;
		vBinVals[2][1][30] = 0.160188;
		vBinVals[2][1][31] = 0.160188;
		vBinVals[2][1][32] = 0.160188;
		vBinVals[2][1][33] = 0.137745;
		vBinVals[2][1][34] = 0.137745;
		vBinVals[2][1][35] = 0.137745;
		vBinVals[2][1][36] = 0.137745;
		vBinVals[2][1][37] = 0.135599;
		vBinVals[2][1][38] = 0.135599;
		vBinVals[2][1][39] = 0.135599;
		vBinVals[2][1][40] = 0.135599;
		vBinVals[2][1][41] = 0.135599;
		vBinVals[2][1][42] = 0.131798;
		vBinVals[2][1][43] = 0.131798;
		vBinVals[2][1][44] = 0.131798;
		vBinVals[2][1][45] = 0.131798;
		vBinVals[2][1][46] = 0.131798;
		vBinVals[2][1][47] = 0.130183;
		vBinVals[2][1][48] = 0.130183;
		vBinVals[2][1][49] = 0.130183;
		vBinVals[2][1][50] = 0.130183;
		vBinVals[2][1][51] = 0.130183;
		vBinVals[2][1][52] = 0.12652;
		vBinVals[2][1][53] = 0.12652;
		vBinVals[2][1][54] = 0.12652;
		vBinVals[2][1][55] = 0.12652;
		vBinVals[2][1][56] = 0.12652;
		vBinVals[2][1][57] = 0.118323;
		vBinVals[2][1][58] = 0.118323;
		vBinVals[2][1][59] = 0.118323;
		vBinVals[2][2][0] = -0.0669441;
		vBinVals[2][2][1] = -0.0669441;
		vBinVals[2][2][2] = -0.0669441;
		vBinVals[2][2][3] = -0.0669441;
		vBinVals[2][2][4] = -0.0669441;
		vBinVals[2][2][5] = -0.0669441;
		vBinVals[2][2][6] = -0.0669441;
		vBinVals[2][2][7] = -0.211174;
		vBinVals[2][2][8] = -0.211174;
		vBinVals[2][2][9] = -0.211174;
		vBinVals[2][2][10] = -0.211174;
		vBinVals[2][2][11] = -0.211174;
		vBinVals[2][2][12] = -0.128739;
		vBinVals[2][2][13] = -0.128739;
		vBinVals[2][2][14] = -0.128739;
		vBinVals[2][2][15] = -0.128739;
		vBinVals[2][2][16] = -0.128739;
		vBinVals[2][2][17] = -0.0814983;
		vBinVals[2][2][18] = -0.0814983;
		vBinVals[2][2][19] = -0.0814983;
		vBinVals[2][2][20] = -0.0814983;
		vBinVals[2][2][21] = -0.0814983;
		vBinVals[2][2][22] = -0.0746794;
		vBinVals[2][2][23] = -0.0746794;
		vBinVals[2][2][24] = -0.0746794;
		vBinVals[2][2][25] = -0.0746794;
		vBinVals[2][2][26] = -0.0746794;
		vBinVals[2][2][27] = -0.0501763;
		vBinVals[2][2][28] = -0.0501763;
		vBinVals[2][2][29] = -0.0501763;
		vBinVals[2][2][30] = -0.0501763;
		vBinVals[2][2][31] = -0.0501763;
		vBinVals[2][2][32] = -0.0501763;
		vBinVals[2][2][33] = -0.0609723;
		vBinVals[2][2][34] = -0.0609723;
		vBinVals[2][2][35] = -0.0609723;
		vBinVals[2][2][36] = -0.0609723;
		vBinVals[2][2][37] = -0.0499201;
		vBinVals[2][2][38] = -0.0499201;
		vBinVals[2][2][39] = -0.0499201;
		vBinVals[2][2][40] = -0.0499201;
		vBinVals[2][2][41] = -0.0499201;
		vBinVals[2][2][42] = -0.037919;
		vBinVals[2][2][43] = -0.037919;
		vBinVals[2][2][44] = -0.037919;
		vBinVals[2][2][45] = -0.037919;
		vBinVals[2][2][46] = -0.037919;
		vBinVals[2][2][47] = -0.0268103;
		vBinVals[2][2][48] = -0.0268103;
		vBinVals[2][2][49] = -0.0268103;
		vBinVals[2][2][50] = -0.0268103;
		vBinVals[2][2][51] = -0.0268103;
		vBinVals[2][2][52] = -0.0224025;
		vBinVals[2][2][53] = -0.0224025;
		vBinVals[2][2][54] = -0.0224025;
		vBinVals[2][2][55] = -0.0224025;
		vBinVals[2][2][56] = -0.0224025;
		vBinVals[2][2][57] = -0.0064964;
		vBinVals[2][2][58] = -0.0064964;
		vBinVals[2][2][59] = -0.0064964;
		vBinVals[2][3][0] = 0.0448832;
		vBinVals[2][3][1] = 0.0448832;
		vBinVals[2][3][2] = 0.0448832;
		vBinVals[2][3][3] = 0.0448832;
		vBinVals[2][3][4] = 0.0448832;
		vBinVals[2][3][5] = 0.0448832;
		vBinVals[2][3][6] = 0.0448832;
		vBinVals[2][3][7] = 0.113789;
		vBinVals[2][3][8] = 0.113789;
		vBinVals[2][3][9] = 0.113789;
		vBinVals[2][3][10] = 0.113789;
		vBinVals[2][3][11] = 0.113789;
		vBinVals[2][3][12] = 0.115631;
		vBinVals[2][3][13] = 0.115631;
		vBinVals[2][3][14] = 0.115631;
		vBinVals[2][3][15] = 0.115631;
		vBinVals[2][3][16] = 0.115631;
		vBinVals[2][3][17] = 0.113547;
		vBinVals[2][3][18] = 0.113547;
		vBinVals[2][3][19] = 0.113547;
		vBinVals[2][3][20] = 0.113547;
		vBinVals[2][3][21] = 0.113547;
		vBinVals[2][3][22] = 0.0995586;
		vBinVals[2][3][23] = 0.0995586;
		vBinVals[2][3][24] = 0.0995586;
		vBinVals[2][3][25] = 0.0995586;
		vBinVals[2][3][26] = 0.0995586;
		vBinVals[2][3][27] = 0.0949353;
		vBinVals[2][3][28] = 0.0949353;
		vBinVals[2][3][29] = 0.0949353;
		vBinVals[2][3][30] = 0.0949353;
		vBinVals[2][3][31] = 0.0949353;
		vBinVals[2][3][32] = 0.0949353;
		vBinVals[2][3][33] = 0.074715;
		vBinVals[2][3][34] = 0.074715;
		vBinVals[2][3][35] = 0.074715;
		vBinVals[2][3][36] = 0.074715;
		vBinVals[2][3][37] = 0.0743902;
		vBinVals[2][3][38] = 0.0743902;
		vBinVals[2][3][39] = 0.0743902;
		vBinVals[2][3][40] = 0.0743902;
		vBinVals[2][3][41] = 0.0743902;
		vBinVals[2][3][42] = 0.0717557;
		vBinVals[2][3][43] = 0.0717557;
		vBinVals[2][3][44] = 0.0717557;
		vBinVals[2][3][45] = 0.0717557;
		vBinVals[2][3][46] = 0.0717557;
		vBinVals[2][3][47] = 0.0713462;
		vBinVals[2][3][48] = 0.0713462;
		vBinVals[2][3][49] = 0.0713462;
		vBinVals[2][3][50] = 0.0713462;
		vBinVals[2][3][51] = 0.0713462;
		vBinVals[2][3][52] = 0.0676397;
		vBinVals[2][3][53] = 0.0676397;
		vBinVals[2][3][54] = 0.0676397;
		vBinVals[2][3][55] = 0.0676397;
		vBinVals[2][3][56] = 0.0676397;
		vBinVals[2][3][57] = 0.0728309;
		vBinVals[2][3][58] = 0.0728309;
		vBinVals[2][3][59] = 0.0728309;
		vBinVals[2][4][0] = 0.352932;
		vBinVals[2][4][1] = 0.352932;
		vBinVals[2][4][2] = 0.352932;
		vBinVals[2][4][3] = 0.352932;
		vBinVals[2][4][4] = 0.352932;
		vBinVals[2][4][5] = 0.352932;
		vBinVals[2][4][6] = 0.352932;
		vBinVals[2][4][7] = 0.16383;
		vBinVals[2][4][8] = 0.16383;
		vBinVals[2][4][9] = 0.16383;
		vBinVals[2][4][10] = 0.16383;
		vBinVals[2][4][11] = 0.16383;
		vBinVals[2][4][12] = 0.125062;
		vBinVals[2][4][13] = 0.125062;
		vBinVals[2][4][14] = 0.125062;
		vBinVals[2][4][15] = 0.125062;
		vBinVals[2][4][16] = 0.125062;
		vBinVals[2][4][17] = 0.120055;
		vBinVals[2][4][18] = 0.120055;
		vBinVals[2][4][19] = 0.120055;
		vBinVals[2][4][20] = 0.120055;
		vBinVals[2][4][21] = 0.120055;
		vBinVals[2][4][22] = 0.154611;
		vBinVals[2][4][23] = 0.154611;
		vBinVals[2][4][24] = 0.154611;
		vBinVals[2][4][25] = 0.154611;
		vBinVals[2][4][26] = 0.154611;
		vBinVals[2][4][27] = 0.165258;
		vBinVals[2][4][28] = 0.165258;
		vBinVals[2][4][29] = 0.165258;
		vBinVals[2][4][30] = 0.165258;
		vBinVals[2][4][31] = 0.165258;
		vBinVals[2][4][32] = 0.165258;
		vBinVals[2][4][33] = 0.252561;
		vBinVals[2][4][34] = 0.252561;
		vBinVals[2][4][35] = 0.252561;
		vBinVals[2][4][36] = 0.252561;
		vBinVals[2][4][37] = 0.224165;
		vBinVals[2][4][38] = 0.224165;
		vBinVals[2][4][39] = 0.224165;
		vBinVals[2][4][40] = 0.224165;
		vBinVals[2][4][41] = 0.224165;
		vBinVals[2][4][42] = 0.253403;
		vBinVals[2][4][43] = 0.253403;
		vBinVals[2][4][44] = 0.253403;
		vBinVals[2][4][45] = 0.253403;
		vBinVals[2][4][46] = 0.253403;
		vBinVals[2][4][47] = 0.245191;
		vBinVals[2][4][48] = 0.245191;
		vBinVals[2][4][49] = 0.245191;
		vBinVals[2][4][50] = 0.245191;
		vBinVals[2][4][51] = 0.245191;
		vBinVals[2][4][52] = 0.236232;
		vBinVals[2][4][53] = 0.236232;
		vBinVals[2][4][54] = 0.236232;
		vBinVals[2][4][55] = 0.236232;
		vBinVals[2][4][56] = 0.236232;
		vBinVals[2][4][57] = 0.245509;
		vBinVals[2][4][58] = 0.245509;
		vBinVals[2][4][59] = 0.245509;
		vBinVals[3][0][0] = -0.222055;
		vBinVals[3][0][1] = -0.222055;
		vBinVals[3][0][2] = -0.222055;
		vBinVals[3][0][3] = -0.222055;
		vBinVals[3][0][4] = -0.222055;
		vBinVals[3][0][5] = -0.222055;
		vBinVals[3][0][6] = -0.222055;
		vBinVals[3][0][7] = -0.0251376;
		vBinVals[3][0][8] = -0.0251376;
		vBinVals[3][0][9] = -0.0251376;
		vBinVals[3][0][10] = -0.0251376;
		vBinVals[3][0][11] = -0.0251376;
		vBinVals[3][0][12] = 0.160462;
		vBinVals[3][0][13] = 0.160462;
		vBinVals[3][0][14] = 0.160462;
		vBinVals[3][0][15] = 0.160462;
		vBinVals[3][0][16] = 0.160462;
		vBinVals[3][0][17] = 0.104165;
		vBinVals[3][0][18] = 0.104165;
		vBinVals[3][0][19] = 0.104165;
		vBinVals[3][0][20] = 0.104165;
		vBinVals[3][0][21] = 0.104165;
		vBinVals[3][0][22] = 0.0437674;
		vBinVals[3][0][23] = 0.0437674;
		vBinVals[3][0][24] = 0.0437674;
		vBinVals[3][0][25] = 0.0437674;
		vBinVals[3][0][26] = 0.0437674;
		vBinVals[3][0][27] = 0.0968675;
		vBinVals[3][0][28] = 0.0968675;
		vBinVals[3][0][29] = 0.0968675;
		vBinVals[3][0][30] = 0.0968675;
		vBinVals[3][0][31] = 0.0968675;
		vBinVals[3][0][32] = 0.0768937;
		vBinVals[3][0][33] = 0.0768937;
		vBinVals[3][0][34] = 0.0768937;
		vBinVals[3][0][35] = 0.0768937;
		vBinVals[3][0][36] = 0.0768937;
		vBinVals[3][0][37] = 0.0768937;
		vBinVals[3][0][38] = 0.0759483;
		vBinVals[3][0][39] = 0.0759483;
		vBinVals[3][0][40] = 0.0759483;
		vBinVals[3][0][41] = 0.0759483;
		vBinVals[3][0][42] = 0.0885227;
		vBinVals[3][0][43] = 0.0885227;
		vBinVals[3][0][44] = 0.0885227;
		vBinVals[3][0][45] = 0.0885227;
		vBinVals[3][0][46] = 0.0885227;
		vBinVals[3][0][47] = 0.0936193;
		vBinVals[3][0][48] = 0.0936193;
		vBinVals[3][0][49] = 0.0936193;
		vBinVals[3][0][50] = 0.0936193;
		vBinVals[3][0][51] = 0.0936193;
		vBinVals[3][0][52] = 0.097238;
		vBinVals[3][0][53] = 0.097238;
		vBinVals[3][0][54] = 0.097238;
		vBinVals[3][0][55] = 0.097238;
		vBinVals[3][0][56] = 0.097238;
		vBinVals[3][0][57] = 0.093037;
		vBinVals[3][0][58] = 0.093037;
		vBinVals[3][0][59] = 0.093037;
		vBinVals[3][1][0] = 0.153614;
		vBinVals[3][1][1] = 0.153614;
		vBinVals[3][1][2] = 0.153614;
		vBinVals[3][1][3] = 0.153614;
		vBinVals[3][1][4] = 0.153614;
		vBinVals[3][1][5] = 0.153614;
		vBinVals[3][1][6] = 0.153614;
		vBinVals[3][1][7] = 0.211381;
		vBinVals[3][1][8] = 0.211381;
		vBinVals[3][1][9] = 0.211381;
		vBinVals[3][1][10] = 0.211381;
		vBinVals[3][1][11] = 0.211381;
		vBinVals[3][1][12] = 0.257737;
		vBinVals[3][1][13] = 0.257737;
		vBinVals[3][1][14] = 0.257737;
		vBinVals[3][1][15] = 0.257737;
		vBinVals[3][1][16] = 0.257737;
		vBinVals[3][1][17] = 0.215258;
		vBinVals[3][1][18] = 0.215258;
		vBinVals[3][1][19] = 0.215258;
		vBinVals[3][1][20] = 0.215258;
		vBinVals[3][1][21] = 0.215258;
		vBinVals[3][1][22] = 0.184073;
		vBinVals[3][1][23] = 0.184073;
		vBinVals[3][1][24] = 0.184073;
		vBinVals[3][1][25] = 0.184073;
		vBinVals[3][1][26] = 0.184073;
		vBinVals[3][1][27] = 0.180221;
		vBinVals[3][1][28] = 0.180221;
		vBinVals[3][1][29] = 0.180221;
		vBinVals[3][1][30] = 0.180221;
		vBinVals[3][1][31] = 0.180221;
		vBinVals[3][1][32] = 0.161382;
		vBinVals[3][1][33] = 0.161382;
		vBinVals[3][1][34] = 0.161382;
		vBinVals[3][1][35] = 0.161382;
		vBinVals[3][1][36] = 0.161382;
		vBinVals[3][1][37] = 0.161382;
		vBinVals[3][1][38] = 0.142435;
		vBinVals[3][1][39] = 0.142435;
		vBinVals[3][1][40] = 0.142435;
		vBinVals[3][1][41] = 0.142435;
		vBinVals[3][1][42] = 0.13811;
		vBinVals[3][1][43] = 0.13811;
		vBinVals[3][1][44] = 0.13811;
		vBinVals[3][1][45] = 0.13811;
		vBinVals[3][1][46] = 0.13811;
		vBinVals[3][1][47] = 0.13401;
		vBinVals[3][1][48] = 0.13401;
		vBinVals[3][1][49] = 0.13401;
		vBinVals[3][1][50] = 0.13401;
		vBinVals[3][1][51] = 0.13401;
		vBinVals[3][1][52] = 0.130081;
		vBinVals[3][1][53] = 0.130081;
		vBinVals[3][1][54] = 0.130081;
		vBinVals[3][1][55] = 0.130081;
		vBinVals[3][1][56] = 0.130081;
		vBinVals[3][1][57] = 0.127016;
		vBinVals[3][1][58] = 0.127016;
		vBinVals[3][1][59] = 0.127016;
		vBinVals[3][2][0] = -0.395829;
		vBinVals[3][2][1] = -0.395829;
		vBinVals[3][2][2] = -0.395829;
		vBinVals[3][2][3] = -0.395829;
		vBinVals[3][2][4] = -0.395829;
		vBinVals[3][2][5] = -0.395829;
		vBinVals[3][2][6] = -0.395829;
		vBinVals[3][2][7] = -0.235363;
		vBinVals[3][2][8] = -0.235363;
		vBinVals[3][2][9] = -0.235363;
		vBinVals[3][2][10] = -0.235363;
		vBinVals[3][2][11] = -0.235363;
		vBinVals[3][2][12] = -0.136434;
		vBinVals[3][2][13] = -0.136434;
		vBinVals[3][2][14] = -0.136434;
		vBinVals[3][2][15] = -0.136434;
		vBinVals[3][2][16] = -0.136434;
		vBinVals[3][2][17] = -0.0771568;
		vBinVals[3][2][18] = -0.0771568;
		vBinVals[3][2][19] = -0.0771568;
		vBinVals[3][2][20] = -0.0771568;
		vBinVals[3][2][21] = -0.0771568;
		vBinVals[3][2][22] = -0.0541099;
		vBinVals[3][2][23] = -0.0541099;
		vBinVals[3][2][24] = -0.0541099;
		vBinVals[3][2][25] = -0.0541099;
		vBinVals[3][2][26] = -0.0541099;
		vBinVals[3][2][27] = -0.0779678;
		vBinVals[3][2][28] = -0.0779678;
		vBinVals[3][2][29] = -0.0779678;
		vBinVals[3][2][30] = -0.0779678;
		vBinVals[3][2][31] = -0.0779678;
		vBinVals[3][2][32] = -0.0474549;
		vBinVals[3][2][33] = -0.0474549;
		vBinVals[3][2][34] = -0.0474549;
		vBinVals[3][2][35] = -0.0474549;
		vBinVals[3][2][36] = -0.0474549;
		vBinVals[3][2][37] = -0.0474549;
		vBinVals[3][2][38] = -0.0484017;
		vBinVals[3][2][39] = -0.0484017;
		vBinVals[3][2][40] = -0.0484017;
		vBinVals[3][2][41] = -0.0484017;
		vBinVals[3][2][42] = -0.0376311;
		vBinVals[3][2][43] = -0.0376311;
		vBinVals[3][2][44] = -0.0376311;
		vBinVals[3][2][45] = -0.0376311;
		vBinVals[3][2][46] = -0.0376311;
		vBinVals[3][2][47] = -0.0221899;
		vBinVals[3][2][48] = -0.0221899;
		vBinVals[3][2][49] = -0.0221899;
		vBinVals[3][2][50] = -0.0221899;
		vBinVals[3][2][51] = -0.0221899;
		vBinVals[3][2][52] = -0.0179116;
		vBinVals[3][2][53] = -0.0179116;
		vBinVals[3][2][54] = -0.0179116;
		vBinVals[3][2][55] = -0.0179116;
		vBinVals[3][2][56] = -0.0179116;
		vBinVals[3][2][57] = -0.0118887;
		vBinVals[3][2][58] = -0.0118887;
		vBinVals[3][2][59] = -0.0118887;
		vBinVals[3][3][0] = 0.0473717;
		vBinVals[3][3][1] = 0.0473717;
		vBinVals[3][3][2] = 0.0473717;
		vBinVals[3][3][3] = 0.0473717;
		vBinVals[3][3][4] = 0.0473717;
		vBinVals[3][3][5] = 0.0473717;
		vBinVals[3][3][6] = 0.0473717;
		vBinVals[3][3][7] = 0.110146;
		vBinVals[3][3][8] = 0.110146;
		vBinVals[3][3][9] = 0.110146;
		vBinVals[3][3][10] = 0.110146;
		vBinVals[3][3][11] = 0.110146;
		vBinVals[3][3][12] = 0.120544;
		vBinVals[3][3][13] = 0.120544;
		vBinVals[3][3][14] = 0.120544;
		vBinVals[3][3][15] = 0.120544;
		vBinVals[3][3][16] = 0.120544;
		vBinVals[3][3][17] = 0.124711;
		vBinVals[3][3][18] = 0.124711;
		vBinVals[3][3][19] = 0.124711;
		vBinVals[3][3][20] = 0.124711;
		vBinVals[3][3][21] = 0.124711;
		vBinVals[3][3][22] = 0.117797;
		vBinVals[3][3][23] = 0.117797;
		vBinVals[3][3][24] = 0.117797;
		vBinVals[3][3][25] = 0.117797;
		vBinVals[3][3][26] = 0.117797;
		vBinVals[3][3][27] = 0.0869056;
		vBinVals[3][3][28] = 0.0869056;
		vBinVals[3][3][29] = 0.0869056;
		vBinVals[3][3][30] = 0.0869056;
		vBinVals[3][3][31] = 0.0869056;
		vBinVals[3][3][32] = 0.0939158;
		vBinVals[3][3][33] = 0.0939158;
		vBinVals[3][3][34] = 0.0939158;
		vBinVals[3][3][35] = 0.0939158;
		vBinVals[3][3][36] = 0.0939158;
		vBinVals[3][3][37] = 0.0939158;
		vBinVals[3][3][38] = 0.080454;
		vBinVals[3][3][39] = 0.080454;
		vBinVals[3][3][40] = 0.080454;
		vBinVals[3][3][41] = 0.080454;
		vBinVals[3][3][42] = 0.0796709;
		vBinVals[3][3][43] = 0.0796709;
		vBinVals[3][3][44] = 0.0796709;
		vBinVals[3][3][45] = 0.0796709;
		vBinVals[3][3][46] = 0.0796709;
		vBinVals[3][3][47] = 0.0812341;
		vBinVals[3][3][48] = 0.0812341;
		vBinVals[3][3][49] = 0.0812341;
		vBinVals[3][3][50] = 0.0812341;
		vBinVals[3][3][51] = 0.0812341;
		vBinVals[3][3][52] = 0.0768822;
		vBinVals[3][3][53] = 0.0768822;
		vBinVals[3][3][54] = 0.0768822;
		vBinVals[3][3][55] = 0.0768822;
		vBinVals[3][3][56] = 0.0768822;
		vBinVals[3][3][57] = 0.0761275;
		vBinVals[3][3][58] = 0.0761275;
		vBinVals[3][3][59] = 0.0761275;
		vBinVals[3][4][0] = 0.573881;
		vBinVals[3][4][1] = 0.573881;
		vBinVals[3][4][2] = 0.573881;
		vBinVals[3][4][3] = 0.573881;
		vBinVals[3][4][4] = 0.573881;
		vBinVals[3][4][5] = 0.573881;
		vBinVals[3][4][6] = 0.573881;
		vBinVals[3][4][7] = 0.239895;
		vBinVals[3][4][8] = 0.239895;
		vBinVals[3][4][9] = 0.239895;
		vBinVals[3][4][10] = 0.239895;
		vBinVals[3][4][11] = 0.239895;
		vBinVals[3][4][12] = 0.127262;
		vBinVals[3][4][13] = 0.127262;
		vBinVals[3][4][14] = 0.127262;
		vBinVals[3][4][15] = 0.127262;
		vBinVals[3][4][16] = 0.127262;
		vBinVals[3][4][17] = 0.0567916;
		vBinVals[3][4][18] = 0.0567916;
		vBinVals[3][4][19] = 0.0567916;
		vBinVals[3][4][20] = 0.0567916;
		vBinVals[3][4][21] = 0.0567916;
		vBinVals[3][4][22] = 0.0924312;
		vBinVals[3][4][23] = 0.0924312;
		vBinVals[3][4][24] = 0.0924312;
		vBinVals[3][4][25] = 0.0924312;
		vBinVals[3][4][26] = 0.0924312;
		vBinVals[3][4][27] = 0.194894;
		vBinVals[3][4][28] = 0.194894;
		vBinVals[3][4][29] = 0.194894;
		vBinVals[3][4][30] = 0.194894;
		vBinVals[3][4][31] = 0.194894;
		vBinVals[3][4][32] = 0.1529;
		vBinVals[3][4][33] = 0.1529;
		vBinVals[3][4][34] = 0.1529;
		vBinVals[3][4][35] = 0.1529;
		vBinVals[3][4][36] = 0.1529;
		vBinVals[3][4][37] = 0.1529;
		vBinVals[3][4][38] = 0.191036;
		vBinVals[3][4][39] = 0.191036;
		vBinVals[3][4][40] = 0.191036;
		vBinVals[3][4][41] = 0.191036;
		vBinVals[3][4][42] = 0.212157;
		vBinVals[3][4][43] = 0.212157;
		vBinVals[3][4][44] = 0.212157;
		vBinVals[3][4][45] = 0.212157;
		vBinVals[3][4][46] = 0.212157;
		vBinVals[3][4][47] = 0.217158;
		vBinVals[3][4][48] = 0.217158;
		vBinVals[3][4][49] = 0.217158;
		vBinVals[3][4][50] = 0.217158;
		vBinVals[3][4][51] = 0.217158;
		vBinVals[3][4][52] = 0.212936;
		vBinVals[3][4][53] = 0.212936;
		vBinVals[3][4][54] = 0.212936;
		vBinVals[3][4][55] = 0.212936;
		vBinVals[3][4][56] = 0.212936;
		vBinVals[3][4][57] = 0.21782;
		vBinVals[3][4][58] = 0.21782;
		vBinVals[3][4][59] = 0.21782;
		vBinVals[4][0][0] = -0.267635;
		vBinVals[4][0][1] = -0.267635;
		vBinVals[4][0][2] = -0.267635;
		vBinVals[4][0][3] = -0.267635;
		vBinVals[4][0][4] = -0.267635;
		vBinVals[4][0][5] = -0.267635;
		vBinVals[4][0][6] = -0.267635;
		vBinVals[4][0][7] = -0.122486;
		vBinVals[4][0][8] = -0.122486;
		vBinVals[4][0][9] = -0.122486;
		vBinVals[4][0][10] = -0.122486;
		vBinVals[4][0][11] = -0.122486;
		vBinVals[4][0][12] = 0.0950279;
		vBinVals[4][0][13] = 0.0950279;
		vBinVals[4][0][14] = 0.0950279;
		vBinVals[4][0][15] = 0.0950279;
		vBinVals[4][0][16] = 0.0950279;
		vBinVals[4][0][17] = 0.119284;
		vBinVals[4][0][18] = 0.119284;
		vBinVals[4][0][19] = 0.119284;
		vBinVals[4][0][20] = 0.119284;
		vBinVals[4][0][21] = 0.119284;
		vBinVals[4][0][22] = 0.0783432;
		vBinVals[4][0][23] = 0.0783432;
		vBinVals[4][0][24] = 0.0783432;
		vBinVals[4][0][25] = 0.0783432;
		vBinVals[4][0][26] = 0.0783432;
		vBinVals[4][0][27] = 0.10727;
		vBinVals[4][0][28] = 0.10727;
		vBinVals[4][0][29] = 0.10727;
		vBinVals[4][0][30] = 0.10727;
		vBinVals[4][0][31] = 0.10727;
		vBinVals[4][0][32] = 0.0933506;
		vBinVals[4][0][33] = 0.0933506;
		vBinVals[4][0][34] = 0.0933506;
		vBinVals[4][0][35] = 0.0933506;
		vBinVals[4][0][36] = 0.0933506;
		vBinVals[4][0][37] = 0.068682;
		vBinVals[4][0][38] = 0.068682;
		vBinVals[4][0][39] = 0.068682;
		vBinVals[4][0][40] = 0.068682;
		vBinVals[4][0][41] = 0.068682;
		vBinVals[4][0][42] = 0.0938318;
		vBinVals[4][0][43] = 0.0938318;
		vBinVals[4][0][44] = 0.0938318;
		vBinVals[4][0][45] = 0.0938318;
		vBinVals[4][0][46] = 0.0938318;
		vBinVals[4][0][47] = 0.114268;
		vBinVals[4][0][48] = 0.114268;
		vBinVals[4][0][49] = 0.114268;
		vBinVals[4][0][50] = 0.114268;
		vBinVals[4][0][51] = 0.114268;
		vBinVals[4][0][52] = 0.111775;
		vBinVals[4][0][53] = 0.111775;
		vBinVals[4][0][54] = 0.111775;
		vBinVals[4][0][55] = 0.111775;
		vBinVals[4][0][56] = 0.111775;
		vBinVals[4][0][57] = 0.114152;
		vBinVals[4][0][58] = 0.114152;
		vBinVals[4][0][59] = 0.114152;
		vBinVals[4][1][0] = 0.134449;
		vBinVals[4][1][1] = 0.134449;
		vBinVals[4][1][2] = 0.134449;
		vBinVals[4][1][3] = 0.134449;
		vBinVals[4][1][4] = 0.134449;
		vBinVals[4][1][5] = 0.134449;
		vBinVals[4][1][6] = 0.134449;
		vBinVals[4][1][7] = 0.180477;
		vBinVals[4][1][8] = 0.180477;
		vBinVals[4][1][9] = 0.180477;
		vBinVals[4][1][10] = 0.180477;
		vBinVals[4][1][11] = 0.180477;
		vBinVals[4][1][12] = 0.234928;
		vBinVals[4][1][13] = 0.234928;
		vBinVals[4][1][14] = 0.234928;
		vBinVals[4][1][15] = 0.234928;
		vBinVals[4][1][16] = 0.234928;
		vBinVals[4][1][17] = 0.250829;
		vBinVals[4][1][18] = 0.250829;
		vBinVals[4][1][19] = 0.250829;
		vBinVals[4][1][20] = 0.250829;
		vBinVals[4][1][21] = 0.250829;
		vBinVals[4][1][22] = 0.186688;
		vBinVals[4][1][23] = 0.186688;
		vBinVals[4][1][24] = 0.186688;
		vBinVals[4][1][25] = 0.186688;
		vBinVals[4][1][26] = 0.186688;
		vBinVals[4][1][27] = 0.19967;
		vBinVals[4][1][28] = 0.19967;
		vBinVals[4][1][29] = 0.19967;
		vBinVals[4][1][30] = 0.19967;
		vBinVals[4][1][31] = 0.19967;
		vBinVals[4][1][32] = 0.172104;
		vBinVals[4][1][33] = 0.172104;
		vBinVals[4][1][34] = 0.172104;
		vBinVals[4][1][35] = 0.172104;
		vBinVals[4][1][36] = 0.172104;
		vBinVals[4][1][37] = 0.153347;
		vBinVals[4][1][38] = 0.153347;
		vBinVals[4][1][39] = 0.153347;
		vBinVals[4][1][40] = 0.153347;
		vBinVals[4][1][41] = 0.153347;
		vBinVals[4][1][42] = 0.157006;
		vBinVals[4][1][43] = 0.157006;
		vBinVals[4][1][44] = 0.157006;
		vBinVals[4][1][45] = 0.157006;
		vBinVals[4][1][46] = 0.157006;
		vBinVals[4][1][47] = 0.148261;
		vBinVals[4][1][48] = 0.148261;
		vBinVals[4][1][49] = 0.148261;
		vBinVals[4][1][50] = 0.148261;
		vBinVals[4][1][51] = 0.148261;
		vBinVals[4][1][52] = 0.146812;
		vBinVals[4][1][53] = 0.146812;
		vBinVals[4][1][54] = 0.146812;
		vBinVals[4][1][55] = 0.146812;
		vBinVals[4][1][56] = 0.146812;
		vBinVals[4][1][57] = 0.141695;
		vBinVals[4][1][58] = 0.141695;
		vBinVals[4][1][59] = 0.141695;
		vBinVals[4][2][0] = 0.143556;
		vBinVals[4][2][1] = 0.143556;
		vBinVals[4][2][2] = 0.143556;
		vBinVals[4][2][3] = 0.143556;
		vBinVals[4][2][4] = 0.143556;
		vBinVals[4][2][5] = 0.143556;
		vBinVals[4][2][6] = 0.143556;
		vBinVals[4][2][7] = 0.143556;
		vBinVals[4][2][8] = 0.143556;
		vBinVals[4][2][9] = 0.143556;
		vBinVals[4][2][10] = 0.143556;
		vBinVals[4][2][11] = 0.143556;
		vBinVals[4][2][12] = -0.129198;
		vBinVals[4][2][13] = -0.129198;
		vBinVals[4][2][14] = -0.129198;
		vBinVals[4][2][15] = -0.129198;
		vBinVals[4][2][16] = -0.129198;
		vBinVals[4][2][17] = -0.0602385;
		vBinVals[4][2][18] = -0.0602385;
		vBinVals[4][2][19] = -0.0602385;
		vBinVals[4][2][20] = -0.0602385;
		vBinVals[4][2][21] = -0.0602385;
		vBinVals[4][2][22] = -0.0570728;
		vBinVals[4][2][23] = -0.0570728;
		vBinVals[4][2][24] = -0.0570728;
		vBinVals[4][2][25] = -0.0570728;
		vBinVals[4][2][26] = -0.0570728;
		vBinVals[4][2][27] = -0.0489995;
		vBinVals[4][2][28] = -0.0489995;
		vBinVals[4][2][29] = -0.0489995;
		vBinVals[4][2][30] = -0.0489995;
		vBinVals[4][2][31] = -0.0489995;
		vBinVals[4][2][32] = -0.0322852;
		vBinVals[4][2][33] = -0.0322852;
		vBinVals[4][2][34] = -0.0322852;
		vBinVals[4][2][35] = -0.0322852;
		vBinVals[4][2][36] = -0.0322852;
		vBinVals[4][2][37] = -0.0255665;
		vBinVals[4][2][38] = -0.0255665;
		vBinVals[4][2][39] = -0.0255665;
		vBinVals[4][2][40] = -0.0255665;
		vBinVals[4][2][41] = -0.0255665;
		vBinVals[4][2][42] = -0.0325868;
		vBinVals[4][2][43] = -0.0325868;
		vBinVals[4][2][44] = -0.0325868;
		vBinVals[4][2][45] = -0.0325868;
		vBinVals[4][2][46] = -0.0325868;
		vBinVals[4][2][47] = -0.0287027;
		vBinVals[4][2][48] = -0.0287027;
		vBinVals[4][2][49] = -0.0287027;
		vBinVals[4][2][50] = -0.0287027;
		vBinVals[4][2][51] = -0.0287027;
		vBinVals[4][2][52] = -0.0148054;
		vBinVals[4][2][53] = -0.0148054;
		vBinVals[4][2][54] = -0.0148054;
		vBinVals[4][2][55] = -0.0148054;
		vBinVals[4][2][56] = -0.0148054;
		vBinVals[4][2][57] = -0.0105999;
		vBinVals[4][2][58] = -0.0105999;
		vBinVals[4][2][59] = -0.0105999;
		vBinVals[4][3][0] = 0.0298375;
		vBinVals[4][3][1] = 0.0298375;
		vBinVals[4][3][2] = 0.0298375;
		vBinVals[4][3][3] = 0.0298375;
		vBinVals[4][3][4] = 0.0298375;
		vBinVals[4][3][5] = 0.0298375;
		vBinVals[4][3][6] = 0.0298375;
		vBinVals[4][3][7] = -0.00420713;
		vBinVals[4][3][8] = -0.00420713;
		vBinVals[4][3][9] = -0.00420713;
		vBinVals[4][3][10] = -0.00420713;
		vBinVals[4][3][11] = -0.00420713;
		vBinVals[4][3][12] = 0.131328;
		vBinVals[4][3][13] = 0.131328;
		vBinVals[4][3][14] = 0.131328;
		vBinVals[4][3][15] = 0.131328;
		vBinVals[4][3][16] = 0.131328;
		vBinVals[4][3][17] = 0.137275;
		vBinVals[4][3][18] = 0.137275;
		vBinVals[4][3][19] = 0.137275;
		vBinVals[4][3][20] = 0.137275;
		vBinVals[4][3][21] = 0.137275;
		vBinVals[4][3][22] = 0.124784;
		vBinVals[4][3][23] = 0.124784;
		vBinVals[4][3][24] = 0.124784;
		vBinVals[4][3][25] = 0.124784;
		vBinVals[4][3][26] = 0.124784;
		vBinVals[4][3][27] = 0.109967;
		vBinVals[4][3][28] = 0.109967;
		vBinVals[4][3][29] = 0.109967;
		vBinVals[4][3][30] = 0.109967;
		vBinVals[4][3][31] = 0.109967;
		vBinVals[4][3][32] = 0.109538;
		vBinVals[4][3][33] = 0.109538;
		vBinVals[4][3][34] = 0.109538;
		vBinVals[4][3][35] = 0.109538;
		vBinVals[4][3][36] = 0.109538;
		vBinVals[4][3][37] = 0.101421;
		vBinVals[4][3][38] = 0.101421;
		vBinVals[4][3][39] = 0.101421;
		vBinVals[4][3][40] = 0.101421;
		vBinVals[4][3][41] = 0.101421;
		vBinVals[4][3][42] = 0.0867476;
		vBinVals[4][3][43] = 0.0867476;
		vBinVals[4][3][44] = 0.0867476;
		vBinVals[4][3][45] = 0.0867476;
		vBinVals[4][3][46] = 0.0867476;
		vBinVals[4][3][47] = 0.0860074;
		vBinVals[4][3][48] = 0.0860074;
		vBinVals[4][3][49] = 0.0860074;
		vBinVals[4][3][50] = 0.0860074;
		vBinVals[4][3][51] = 0.0860074;
		vBinVals[4][3][52] = 0.0860586;
		vBinVals[4][3][53] = 0.0860586;
		vBinVals[4][3][54] = 0.0860586;
		vBinVals[4][3][55] = 0.0860586;
		vBinVals[4][3][56] = 0.0860586;
		vBinVals[4][3][57] = 0.0835267;
		vBinVals[4][3][58] = 0.0835267;
		vBinVals[4][3][59] = 0.0835267;
		vBinVals[4][4][0] = 0.463549;
		vBinVals[4][4][1] = 0.463549;
		vBinVals[4][4][2] = 0.463549;
		vBinVals[4][4][3] = 0.463549;
		vBinVals[4][4][4] = 0.463549;
		vBinVals[4][4][5] = 0.463549;
		vBinVals[4][4][6] = 0.463549;
		vBinVals[4][4][7] = 0.630603;
		vBinVals[4][4][8] = 0.630603;
		vBinVals[4][4][9] = 0.630603;
		vBinVals[4][4][10] = 0.630603;
		vBinVals[4][4][11] = 0.630603;
		vBinVals[4][4][12] = 0.152941;
		vBinVals[4][4][13] = 0.152941;
		vBinVals[4][4][14] = 0.152941;
		vBinVals[4][4][15] = 0.152941;
		vBinVals[4][4][16] = 0.152941;
		vBinVals[4][4][17] = 0.070523;
		vBinVals[4][4][18] = 0.070523;
		vBinVals[4][4][19] = 0.070523;
		vBinVals[4][4][20] = 0.070523;
		vBinVals[4][4][21] = 0.070523;
		vBinVals[4][4][22] = 0.0930897;
		vBinVals[4][4][23] = 0.0930897;
		vBinVals[4][4][24] = 0.0930897;
		vBinVals[4][4][25] = 0.0930897;
		vBinVals[4][4][26] = 0.0930897;
		vBinVals[4][4][27] = 0.108725;
		vBinVals[4][4][28] = 0.108725;
		vBinVals[4][4][29] = 0.108725;
		vBinVals[4][4][30] = 0.108725;
		vBinVals[4][4][31] = 0.108725;
		vBinVals[4][4][32] = 0.231195;
		vBinVals[4][4][33] = 0.231195;
		vBinVals[4][4][34] = 0.231195;
		vBinVals[4][4][35] = 0.231195;
		vBinVals[4][4][36] = 0.231195;
		vBinVals[4][4][37] = 0.173485;
		vBinVals[4][4][38] = 0.173485;
		vBinVals[4][4][39] = 0.173485;
		vBinVals[4][4][40] = 0.173485;
		vBinVals[4][4][41] = 0.173485;
		vBinVals[4][4][42] = 0.183687;
		vBinVals[4][4][43] = 0.183687;
		vBinVals[4][4][44] = 0.183687;
		vBinVals[4][4][45] = 0.183687;
		vBinVals[4][4][46] = 0.183687;
		vBinVals[4][4][47] = 0.215623;
		vBinVals[4][4][48] = 0.215623;
		vBinVals[4][4][49] = 0.215623;
		vBinVals[4][4][50] = 0.215623;
		vBinVals[4][4][51] = 0.215623;
		vBinVals[4][4][52] = 0.217431;
		vBinVals[4][4][53] = 0.217431;
		vBinVals[4][4][54] = 0.217431;
		vBinVals[4][4][55] = 0.217431;
		vBinVals[4][4][56] = 0.217431;
		vBinVals[4][4][57] = 0.247143;
		vBinVals[4][4][58] = 0.247143;
		vBinVals[4][4][59] = 0.247143;
		vBinVals[5][0][0] = 4.7246;
		vBinVals[5][0][1] = 4.7246;
		vBinVals[5][0][2] = 4.7246;
		vBinVals[5][0][3] = 4.7246;
		vBinVals[5][0][4] = 4.7246;
		vBinVals[5][0][5] = 4.7246;
		vBinVals[5][0][6] = 4.7246;
		vBinVals[5][0][7] = -0.122116;
		vBinVals[5][0][8] = -0.122116;
		vBinVals[5][0][9] = -0.122116;
		vBinVals[5][0][10] = -0.122116;
		vBinVals[5][0][11] = -0.122116;
		vBinVals[5][0][12] = 0.00356392;
		vBinVals[5][0][13] = 0.00356392;
		vBinVals[5][0][14] = 0.00356392;
		vBinVals[5][0][15] = 0.00356392;
		vBinVals[5][0][16] = 0.00356392;
		vBinVals[5][0][17] = 0.0911727;
		vBinVals[5][0][18] = 0.0911727;
		vBinVals[5][0][19] = 0.0911727;
		vBinVals[5][0][20] = 0.0911727;
		vBinVals[5][0][21] = 0.0911727;
		vBinVals[5][0][22] = 0.0725482;
		vBinVals[5][0][23] = 0.0725482;
		vBinVals[5][0][24] = 0.0725482;
		vBinVals[5][0][25] = 0.0725482;
		vBinVals[5][0][26] = 0.0725482;
		vBinVals[5][0][27] = 0.0339492;
		vBinVals[5][0][28] = 0.0339492;
		vBinVals[5][0][29] = 0.0339492;
		vBinVals[5][0][30] = 0.0339492;
		vBinVals[5][0][31] = 0.0339492;
		vBinVals[5][0][32] = 0.0753558;
		vBinVals[5][0][33] = 0.0753558;
		vBinVals[5][0][34] = 0.0753558;
		vBinVals[5][0][35] = 0.0753558;
		vBinVals[5][0][36] = 0.0753558;
		vBinVals[5][0][37] = 0.066168;
		vBinVals[5][0][38] = 0.066168;
		vBinVals[5][0][39] = 0.066168;
		vBinVals[5][0][40] = 0.066168;
		vBinVals[5][0][41] = 0.066168;
		vBinVals[5][0][42] = 0.0574576;
		vBinVals[5][0][43] = 0.0574576;
		vBinVals[5][0][44] = 0.0574576;
		vBinVals[5][0][45] = 0.0574576;
		vBinVals[5][0][46] = 0.0574576;
		vBinVals[5][0][47] = 0.0676833;
		vBinVals[5][0][48] = 0.0676833;
		vBinVals[5][0][49] = 0.0676833;
		vBinVals[5][0][50] = 0.0676833;
		vBinVals[5][0][51] = 0.0676833;
		vBinVals[5][0][52] = 0.0742555;
		vBinVals[5][0][53] = 0.0742555;
		vBinVals[5][0][54] = 0.0742555;
		vBinVals[5][0][55] = 0.0742555;
		vBinVals[5][0][56] = 0.0742555;
		vBinVals[5][0][57] = 0.0790778;
		vBinVals[5][0][58] = 0.0790778;
		vBinVals[5][0][59] = 0.0790778;
		vBinVals[5][1][0] = 1.23809;
		vBinVals[5][1][1] = 1.23809;
		vBinVals[5][1][2] = 1.23809;
		vBinVals[5][1][3] = 1.23809;
		vBinVals[5][1][4] = 1.23809;
		vBinVals[5][1][5] = 1.23809;
		vBinVals[5][1][6] = 1.23809;
		vBinVals[5][1][7] = 0.183209;
		vBinVals[5][1][8] = 0.183209;
		vBinVals[5][1][9] = 0.183209;
		vBinVals[5][1][10] = 0.183209;
		vBinVals[5][1][11] = 0.183209;
		vBinVals[5][1][12] = 0.20588;
		vBinVals[5][1][13] = 0.20588;
		vBinVals[5][1][14] = 0.20588;
		vBinVals[5][1][15] = 0.20588;
		vBinVals[5][1][16] = 0.20588;
		vBinVals[5][1][17] = 0.218332;
		vBinVals[5][1][18] = 0.218332;
		vBinVals[5][1][19] = 0.218332;
		vBinVals[5][1][20] = 0.218332;
		vBinVals[5][1][21] = 0.218332;
		vBinVals[5][1][22] = 0.209701;
		vBinVals[5][1][23] = 0.209701;
		vBinVals[5][1][24] = 0.209701;
		vBinVals[5][1][25] = 0.209701;
		vBinVals[5][1][26] = 0.209701;
		vBinVals[5][1][27] = 0.178356;
		vBinVals[5][1][28] = 0.178356;
		vBinVals[5][1][29] = 0.178356;
		vBinVals[5][1][30] = 0.178356;
		vBinVals[5][1][31] = 0.178356;
		vBinVals[5][1][32] = 0.197652;
		vBinVals[5][1][33] = 0.197652;
		vBinVals[5][1][34] = 0.197652;
		vBinVals[5][1][35] = 0.197652;
		vBinVals[5][1][36] = 0.197652;
		vBinVals[5][1][37] = 0.182954;
		vBinVals[5][1][38] = 0.182954;
		vBinVals[5][1][39] = 0.182954;
		vBinVals[5][1][40] = 0.182954;
		vBinVals[5][1][41] = 0.182954;
		vBinVals[5][1][42] = 0.149236;
		vBinVals[5][1][43] = 0.149236;
		vBinVals[5][1][44] = 0.149236;
		vBinVals[5][1][45] = 0.149236;
		vBinVals[5][1][46] = 0.149236;
		vBinVals[5][1][47] = 0.14571;
		vBinVals[5][1][48] = 0.14571;
		vBinVals[5][1][49] = 0.14571;
		vBinVals[5][1][50] = 0.14571;
		vBinVals[5][1][51] = 0.14571;
		vBinVals[5][1][52] = 0.143421;
		vBinVals[5][1][53] = 0.143421;
		vBinVals[5][1][54] = 0.143421;
		vBinVals[5][1][55] = 0.143421;
		vBinVals[5][1][56] = 0.143421;
		vBinVals[5][1][57] = 0.144054;
		vBinVals[5][1][58] = 0.144054;
		vBinVals[5][1][59] = 0.144054;
		vBinVals[5][2][0] = -0.249178;
		vBinVals[5][2][1] = -0.249178;
		vBinVals[5][2][2] = -0.249178;
		vBinVals[5][2][3] = -0.249178;
		vBinVals[5][2][4] = -0.249178;
		vBinVals[5][2][5] = -0.249178;
		vBinVals[5][2][6] = -0.249178;
		vBinVals[5][2][7] = -0.299278;
		vBinVals[5][2][8] = -0.299278;
		vBinVals[5][2][9] = -0.299278;
		vBinVals[5][2][10] = -0.299278;
		vBinVals[5][2][11] = -0.299278;
		vBinVals[5][2][12] = -0.1758;
		vBinVals[5][2][13] = -0.1758;
		vBinVals[5][2][14] = -0.1758;
		vBinVals[5][2][15] = -0.1758;
		vBinVals[5][2][16] = -0.1758;
		vBinVals[5][2][17] = -0.134005;
		vBinVals[5][2][18] = -0.134005;
		vBinVals[5][2][19] = -0.134005;
		vBinVals[5][2][20] = -0.134005;
		vBinVals[5][2][21] = -0.134005;
		vBinVals[5][2][22] = -0.0898541;
		vBinVals[5][2][23] = -0.0898541;
		vBinVals[5][2][24] = -0.0898541;
		vBinVals[5][2][25] = -0.0898541;
		vBinVals[5][2][26] = -0.0898541;
		vBinVals[5][2][27] = -0.0756562;
		vBinVals[5][2][28] = -0.0756562;
		vBinVals[5][2][29] = -0.0756562;
		vBinVals[5][2][30] = -0.0756562;
		vBinVals[5][2][31] = -0.0756562;
		vBinVals[5][2][32] = -0.0725964;
		vBinVals[5][2][33] = -0.0725964;
		vBinVals[5][2][34] = -0.0725964;
		vBinVals[5][2][35] = -0.0725964;
		vBinVals[5][2][36] = -0.0725964;
		vBinVals[5][2][37] = -0.0478089;
		vBinVals[5][2][38] = -0.0478089;
		vBinVals[5][2][39] = -0.0478089;
		vBinVals[5][2][40] = -0.0478089;
		vBinVals[5][2][41] = -0.0478089;
		vBinVals[5][2][42] = -0.0525034;
		vBinVals[5][2][43] = -0.0525034;
		vBinVals[5][2][44] = -0.0525034;
		vBinVals[5][2][45] = -0.0525034;
		vBinVals[5][2][46] = -0.0525034;
		vBinVals[5][2][47] = -0.054793;
		vBinVals[5][2][48] = -0.054793;
		vBinVals[5][2][49] = -0.054793;
		vBinVals[5][2][50] = -0.054793;
		vBinVals[5][2][51] = -0.054793;
		vBinVals[5][2][52] = -0.0677279;
		vBinVals[5][2][53] = -0.0677279;
		vBinVals[5][2][54] = -0.0677279;
		vBinVals[5][2][55] = -0.0677279;
		vBinVals[5][2][56] = -0.0677279;
		vBinVals[5][2][57] = -0.0541095;
		vBinVals[5][2][58] = -0.0541095;
		vBinVals[5][2][59] = -0.0541095;
		vBinVals[5][3][0] = 0.00157587;
		vBinVals[5][3][1] = 0.00157587;
		vBinVals[5][3][2] = 0.00157587;
		vBinVals[5][3][3] = 0.00157587;
		vBinVals[5][3][4] = 0.00157587;
		vBinVals[5][3][5] = 0.00157587;
		vBinVals[5][3][6] = 0.00157587;
		vBinVals[5][3][7] = 0.0979949;
		vBinVals[5][3][8] = 0.0979949;
		vBinVals[5][3][9] = 0.0979949;
		vBinVals[5][3][10] = 0.0979949;
		vBinVals[5][3][11] = 0.0979949;
		vBinVals[5][3][12] = 0.118912;
		vBinVals[5][3][13] = 0.118912;
		vBinVals[5][3][14] = 0.118912;
		vBinVals[5][3][15] = 0.118912;
		vBinVals[5][3][16] = 0.118912;
		vBinVals[5][3][17] = 0.118345;
		vBinVals[5][3][18] = 0.118345;
		vBinVals[5][3][19] = 0.118345;
		vBinVals[5][3][20] = 0.118345;
		vBinVals[5][3][21] = 0.118345;
		vBinVals[5][3][22] = 0.121969;
		vBinVals[5][3][23] = 0.121969;
		vBinVals[5][3][24] = 0.121969;
		vBinVals[5][3][25] = 0.121969;
		vBinVals[5][3][26] = 0.121969;
		vBinVals[5][3][27] = 0.118968;
		vBinVals[5][3][28] = 0.118968;
		vBinVals[5][3][29] = 0.118968;
		vBinVals[5][3][30] = 0.118968;
		vBinVals[5][3][31] = 0.118968;
		vBinVals[5][3][32] = 0.0958429;
		vBinVals[5][3][33] = 0.0958429;
		vBinVals[5][3][34] = 0.0958429;
		vBinVals[5][3][35] = 0.0958429;
		vBinVals[5][3][36] = 0.0958429;
		vBinVals[5][3][37] = 0.0998076;
		vBinVals[5][3][38] = 0.0998076;
		vBinVals[5][3][39] = 0.0998076;
		vBinVals[5][3][40] = 0.0998076;
		vBinVals[5][3][41] = 0.0998076;
		vBinVals[5][3][42] = 0.0981446;
		vBinVals[5][3][43] = 0.0981446;
		vBinVals[5][3][44] = 0.0981446;
		vBinVals[5][3][45] = 0.0981446;
		vBinVals[5][3][46] = 0.0981446;
		vBinVals[5][3][47] = 0.0806393;
		vBinVals[5][3][48] = 0.0806393;
		vBinVals[5][3][49] = 0.0806393;
		vBinVals[5][3][50] = 0.0806393;
		vBinVals[5][3][51] = 0.0806393;
		vBinVals[5][3][52] = 0.0714543;
		vBinVals[5][3][53] = 0.0714543;
		vBinVals[5][3][54] = 0.0714543;
		vBinVals[5][3][55] = 0.0714543;
		vBinVals[5][3][56] = 0.0714543;
		vBinVals[5][3][57] = 0.0690284;
		vBinVals[5][3][58] = 0.0690284;
		vBinVals[5][3][59] = 0.0690284;
		vBinVals[5][4][0] = 29.5814;
		vBinVals[5][4][1] = 29.5814;
		vBinVals[5][4][2] = 29.5814;
		vBinVals[5][4][3] = 29.5814;
		vBinVals[5][4][4] = 29.5814;
		vBinVals[5][4][5] = 29.5814;
		vBinVals[5][4][6] = 29.5814;
		vBinVals[5][4][7] = 0.230575;
		vBinVals[5][4][8] = 0.230575;
		vBinVals[5][4][9] = 0.230575;
		vBinVals[5][4][10] = 0.230575;
		vBinVals[5][4][11] = 0.230575;
		vBinVals[5][4][12] = 0.178568;
		vBinVals[5][4][13] = 0.178568;
		vBinVals[5][4][14] = 0.178568;
		vBinVals[5][4][15] = 0.178568;
		vBinVals[5][4][16] = 0.178568;
		vBinVals[5][4][17] = 0.185771;
		vBinVals[5][4][18] = 0.185771;
		vBinVals[5][4][19] = 0.185771;
		vBinVals[5][4][20] = 0.185771;
		vBinVals[5][4][21] = 0.185771;
		vBinVals[5][4][22] = 0.0995015;
		vBinVals[5][4][23] = 0.0995015;
		vBinVals[5][4][24] = 0.0995015;
		vBinVals[5][4][25] = 0.0995015;
		vBinVals[5][4][26] = 0.0995015;
		vBinVals[5][4][27] = 0.132787;
		vBinVals[5][4][28] = 0.132787;
		vBinVals[5][4][29] = 0.132787;
		vBinVals[5][4][30] = 0.132787;
		vBinVals[5][4][31] = 0.132787;
		vBinVals[5][4][32] = 0.139186;
		vBinVals[5][4][33] = 0.139186;
		vBinVals[5][4][34] = 0.139186;
		vBinVals[5][4][35] = 0.139186;
		vBinVals[5][4][36] = 0.139186;
		vBinVals[5][4][37] = 0.123576;
		vBinVals[5][4][38] = 0.123576;
		vBinVals[5][4][39] = 0.123576;
		vBinVals[5][4][40] = 0.123576;
		vBinVals[5][4][41] = 0.123576;
		vBinVals[5][4][42] = 0.252613;
		vBinVals[5][4][43] = 0.252613;
		vBinVals[5][4][44] = 0.252613;
		vBinVals[5][4][45] = 0.252613;
		vBinVals[5][4][46] = 0.252613;
		vBinVals[5][4][47] = 0.200182;
		vBinVals[5][4][48] = 0.200182;
		vBinVals[5][4][49] = 0.200182;
		vBinVals[5][4][50] = 0.200182;
		vBinVals[5][4][51] = 0.200182;
		vBinVals[5][4][52] = 0.252383;
		vBinVals[5][4][53] = 0.252383;
		vBinVals[5][4][54] = 0.252383;
		vBinVals[5][4][55] = 0.252383;
		vBinVals[5][4][56] = 0.252383;
		vBinVals[5][4][57] = 0.219761;
		vBinVals[5][4][58] = 0.219761;
		vBinVals[5][4][59] = 0.219761;
		vBinVals[6][0][0] = -0.183619;
		vBinVals[6][0][1] = -0.183619;
		vBinVals[6][0][2] = -0.183619;
		vBinVals[6][0][3] = -0.183619;
		vBinVals[6][0][4] = -0.183619;
		vBinVals[6][0][5] = -0.183619;
		vBinVals[6][0][6] = -0.183619;
		vBinVals[6][0][7] = -0.183619;
		vBinVals[6][0][8] = -0.183619;
		vBinVals[6][0][9] = -0.183619;
		vBinVals[6][0][10] = -0.183619;
		vBinVals[6][0][11] = -0.183619;
		vBinVals[6][0][12] = -0.153494;
		vBinVals[6][0][13] = -0.153494;
		vBinVals[6][0][14] = -0.153494;
		vBinVals[6][0][15] = -0.153494;
		vBinVals[6][0][16] = -0.153494;
		vBinVals[6][0][17] = 0.0527443;
		vBinVals[6][0][18] = 0.0527443;
		vBinVals[6][0][19] = 0.0527443;
		vBinVals[6][0][20] = 0.0527443;
		vBinVals[6][0][21] = 0.0527443;
		vBinVals[6][0][22] = 0.0694586;
		vBinVals[6][0][23] = 0.0694586;
		vBinVals[6][0][24] = 0.0694586;
		vBinVals[6][0][25] = 0.0694586;
		vBinVals[6][0][26] = 0.0694586;
		vBinVals[6][0][27] = 0.0652551;
		vBinVals[6][0][28] = 0.0652551;
		vBinVals[6][0][29] = 0.0652551;
		vBinVals[6][0][30] = 0.0652551;
		vBinVals[6][0][31] = 0.0652551;
		vBinVals[6][0][32] = 0.0144285;
		vBinVals[6][0][33] = 0.0144285;
		vBinVals[6][0][34] = 0.0144285;
		vBinVals[6][0][35] = 0.0144285;
		vBinVals[6][0][36] = 0.0144285;
		vBinVals[6][0][37] = 0.0295133;
		vBinVals[6][0][38] = 0.0295133;
		vBinVals[6][0][39] = 0.0295133;
		vBinVals[6][0][40] = 0.0295133;
		vBinVals[6][0][41] = 0.0295133;
		vBinVals[6][0][42] = 0.0572422;
		vBinVals[6][0][43] = 0.0572422;
		vBinVals[6][0][44] = 0.0572422;
		vBinVals[6][0][45] = 0.0572422;
		vBinVals[6][0][46] = 0.0572422;
		vBinVals[6][0][47] = 0.0550676;
		vBinVals[6][0][48] = 0.0550676;
		vBinVals[6][0][49] = 0.0550676;
		vBinVals[6][0][50] = 0.0550676;
		vBinVals[6][0][51] = 0.0550676;
		vBinVals[6][0][52] = 0.0413465;
		vBinVals[6][0][53] = 0.0413465;
		vBinVals[6][0][54] = 0.0413465;
		vBinVals[6][0][55] = 0.0413465;
		vBinVals[6][0][56] = 0.0413465;
		vBinVals[6][0][57] = 0.056139;
		vBinVals[6][0][58] = 0.056139;
		vBinVals[6][0][59] = 0.056139;
		vBinVals[6][1][0] = -0.156272;
		vBinVals[6][1][1] = -0.156272;
		vBinVals[6][1][2] = -0.156272;
		vBinVals[6][1][3] = -0.156272;
		vBinVals[6][1][4] = -0.156272;
		vBinVals[6][1][5] = -0.156272;
		vBinVals[6][1][6] = -0.156272;
		vBinVals[6][1][7] = -0.156272;
		vBinVals[6][1][8] = -0.156272;
		vBinVals[6][1][9] = -0.156272;
		vBinVals[6][1][10] = -0.156272;
		vBinVals[6][1][11] = -0.156272;
		vBinVals[6][1][12] = 0.145192;
		vBinVals[6][1][13] = 0.145192;
		vBinVals[6][1][14] = 0.145192;
		vBinVals[6][1][15] = 0.145192;
		vBinVals[6][1][16] = 0.145192;
		vBinVals[6][1][17] = 0.21055;
		vBinVals[6][1][18] = 0.21055;
		vBinVals[6][1][19] = 0.21055;
		vBinVals[6][1][20] = 0.21055;
		vBinVals[6][1][21] = 0.21055;
		vBinVals[6][1][22] = 0.22866;
		vBinVals[6][1][23] = 0.22866;
		vBinVals[6][1][24] = 0.22866;
		vBinVals[6][1][25] = 0.22866;
		vBinVals[6][1][26] = 0.22866;
		vBinVals[6][1][27] = 0.192746;
		vBinVals[6][1][28] = 0.192746;
		vBinVals[6][1][29] = 0.192746;
		vBinVals[6][1][30] = 0.192746;
		vBinVals[6][1][31] = 0.192746;
		vBinVals[6][1][32] = 0.175927;
		vBinVals[6][1][33] = 0.175927;
		vBinVals[6][1][34] = 0.175927;
		vBinVals[6][1][35] = 0.175927;
		vBinVals[6][1][36] = 0.175927;
		vBinVals[6][1][37] = 0.15912;
		vBinVals[6][1][38] = 0.15912;
		vBinVals[6][1][39] = 0.15912;
		vBinVals[6][1][40] = 0.15912;
		vBinVals[6][1][41] = 0.15912;
		vBinVals[6][1][42] = 0.180591;
		vBinVals[6][1][43] = 0.180591;
		vBinVals[6][1][44] = 0.180591;
		vBinVals[6][1][45] = 0.180591;
		vBinVals[6][1][46] = 0.180591;
		vBinVals[6][1][47] = 0.152734;
		vBinVals[6][1][48] = 0.152734;
		vBinVals[6][1][49] = 0.152734;
		vBinVals[6][1][50] = 0.152734;
		vBinVals[6][1][51] = 0.152734;
		vBinVals[6][1][52] = 0.130461;
		vBinVals[6][1][53] = 0.130461;
		vBinVals[6][1][54] = 0.130461;
		vBinVals[6][1][55] = 0.130461;
		vBinVals[6][1][56] = 0.130461;
		vBinVals[6][1][57] = 0.131361;
		vBinVals[6][1][58] = 0.131361;
		vBinVals[6][1][59] = 0.131361;
		vBinVals[6][2][0] = -0.292264;
		vBinVals[6][2][1] = -0.292264;
		vBinVals[6][2][2] = -0.292264;
		vBinVals[6][2][3] = -0.292264;
		vBinVals[6][2][4] = -0.292264;
		vBinVals[6][2][5] = -0.292264;
		vBinVals[6][2][6] = -0.292264;
		vBinVals[6][2][7] = -0.292264;
		vBinVals[6][2][8] = -0.292264;
		vBinVals[6][2][9] = -0.292264;
		vBinVals[6][2][10] = -0.292264;
		vBinVals[6][2][11] = -0.292264;
		vBinVals[6][2][12] = 0.118832;
		vBinVals[6][2][13] = 0.118832;
		vBinVals[6][2][14] = 0.118832;
		vBinVals[6][2][15] = 0.118832;
		vBinVals[6][2][16] = 0.118832;
		vBinVals[6][2][17] = -0.171072;
		vBinVals[6][2][18] = -0.171072;
		vBinVals[6][2][19] = -0.171072;
		vBinVals[6][2][20] = -0.171072;
		vBinVals[6][2][21] = -0.171072;
		vBinVals[6][2][22] = -0.11419;
		vBinVals[6][2][23] = -0.11419;
		vBinVals[6][2][24] = -0.11419;
		vBinVals[6][2][25] = -0.11419;
		vBinVals[6][2][26] = -0.11419;
		vBinVals[6][2][27] = -0.109441;
		vBinVals[6][2][28] = -0.109441;
		vBinVals[6][2][29] = -0.109441;
		vBinVals[6][2][30] = -0.109441;
		vBinVals[6][2][31] = -0.109441;
		vBinVals[6][2][32] = -0.0962436;
		vBinVals[6][2][33] = -0.0962436;
		vBinVals[6][2][34] = -0.0962436;
		vBinVals[6][2][35] = -0.0962436;
		vBinVals[6][2][36] = -0.0962436;
		vBinVals[6][2][37] = -0.0766795;
		vBinVals[6][2][38] = -0.0766795;
		vBinVals[6][2][39] = -0.0766795;
		vBinVals[6][2][40] = -0.0766795;
		vBinVals[6][2][41] = -0.0766795;
		vBinVals[6][2][42] = -0.0793521;
		vBinVals[6][2][43] = -0.0793521;
		vBinVals[6][2][44] = -0.0793521;
		vBinVals[6][2][45] = -0.0793521;
		vBinVals[6][2][46] = -0.0793521;
		vBinVals[6][2][47] = -0.0822503;
		vBinVals[6][2][48] = -0.0822503;
		vBinVals[6][2][49] = -0.0822503;
		vBinVals[6][2][50] = -0.0822503;
		vBinVals[6][2][51] = -0.0822503;
		vBinVals[6][2][52] = -0.0758047;
		vBinVals[6][2][53] = -0.0758047;
		vBinVals[6][2][54] = -0.0758047;
		vBinVals[6][2][55] = -0.0758047;
		vBinVals[6][2][56] = -0.0758047;
		vBinVals[6][2][57] = -0.0690533;
		vBinVals[6][2][58] = -0.0690533;
		vBinVals[6][2][59] = -0.0690533;
		vBinVals[6][3][0] = 0.106257;
		vBinVals[6][3][1] = 0.106257;
		vBinVals[6][3][2] = 0.106257;
		vBinVals[6][3][3] = 0.106257;
		vBinVals[6][3][4] = 0.106257;
		vBinVals[6][3][5] = 0.106257;
		vBinVals[6][3][6] = 0.106257;
		vBinVals[6][3][7] = 0.106257;
		vBinVals[6][3][8] = 0.106257;
		vBinVals[6][3][9] = 0.106257;
		vBinVals[6][3][10] = 0.106257;
		vBinVals[6][3][11] = 0.106257;
		vBinVals[6][3][12] = 0.0705678;
		vBinVals[6][3][13] = 0.0705678;
		vBinVals[6][3][14] = 0.0705678;
		vBinVals[6][3][15] = 0.0705678;
		vBinVals[6][3][16] = 0.0705678;
		vBinVals[6][3][17] = 0.102789;
		vBinVals[6][3][18] = 0.102789;
		vBinVals[6][3][19] = 0.102789;
		vBinVals[6][3][20] = 0.102789;
		vBinVals[6][3][21] = 0.102789;
		vBinVals[6][3][22] = 0.107833;
		vBinVals[6][3][23] = 0.107833;
		vBinVals[6][3][24] = 0.107833;
		vBinVals[6][3][25] = 0.107833;
		vBinVals[6][3][26] = 0.107833;
		vBinVals[6][3][27] = 0.0978463;
		vBinVals[6][3][28] = 0.0978463;
		vBinVals[6][3][29] = 0.0978463;
		vBinVals[6][3][30] = 0.0978463;
		vBinVals[6][3][31] = 0.0978463;
		vBinVals[6][3][32] = 0.0937328;
		vBinVals[6][3][33] = 0.0937328;
		vBinVals[6][3][34] = 0.0937328;
		vBinVals[6][3][35] = 0.0937328;
		vBinVals[6][3][36] = 0.0937328;
		vBinVals[6][3][37] = 0.0923728;
		vBinVals[6][3][38] = 0.0923728;
		vBinVals[6][3][39] = 0.0923728;
		vBinVals[6][3][40] = 0.0923728;
		vBinVals[6][3][41] = 0.0923728;
		vBinVals[6][3][42] = 0.0741421;
		vBinVals[6][3][43] = 0.0741421;
		vBinVals[6][3][44] = 0.0741421;
		vBinVals[6][3][45] = 0.0741421;
		vBinVals[6][3][46] = 0.0741421;
		vBinVals[6][3][47] = 0.0717059;
		vBinVals[6][3][48] = 0.0717059;
		vBinVals[6][3][49] = 0.0717059;
		vBinVals[6][3][50] = 0.0717059;
		vBinVals[6][3][51] = 0.0717059;
		vBinVals[6][3][52] = 0.0689247;
		vBinVals[6][3][53] = 0.0689247;
		vBinVals[6][3][54] = 0.0689247;
		vBinVals[6][3][55] = 0.0689247;
		vBinVals[6][3][56] = 0.0689247;
		vBinVals[6][3][57] = 0.0657103;
		vBinVals[6][3][58] = 0.0657103;
		vBinVals[6][3][59] = 0.0657103;
		vBinVals[6][4][0] = 0.491252;
		vBinVals[6][4][1] = 0.491252;
		vBinVals[6][4][2] = 0.491252;
		vBinVals[6][4][3] = 0.491252;
		vBinVals[6][4][4] = 0.491252;
		vBinVals[6][4][5] = 0.491252;
		vBinVals[6][4][6] = 0.491252;
		vBinVals[6][4][7] = 0.491252;
		vBinVals[6][4][8] = 0.491252;
		vBinVals[6][4][9] = 0.491252;
		vBinVals[6][4][10] = 0.491252;
		vBinVals[6][4][11] = 0.491252;
		vBinVals[6][4][12] = 0.450422;
		vBinVals[6][4][13] = 0.450422;
		vBinVals[6][4][14] = 0.450422;
		vBinVals[6][4][15] = 0.450422;
		vBinVals[6][4][16] = 0.450422;
		vBinVals[6][4][17] = 0.192092;
		vBinVals[6][4][18] = 0.192092;
		vBinVals[6][4][19] = 0.192092;
		vBinVals[6][4][20] = 0.192092;
		vBinVals[6][4][21] = 0.192092;
		vBinVals[6][4][22] = 0.0992901;
		vBinVals[6][4][23] = 0.0992901;
		vBinVals[6][4][24] = 0.0992901;
		vBinVals[6][4][25] = 0.0992901;
		vBinVals[6][4][26] = 0.0992901;
		vBinVals[6][4][27] = 0.104535;
		vBinVals[6][4][28] = 0.104535;
		vBinVals[6][4][29] = 0.104535;
		vBinVals[6][4][30] = 0.104535;
		vBinVals[6][4][31] = 0.104535;
		vBinVals[6][4][32] = 0.125538;
		vBinVals[6][4][33] = 0.125538;
		vBinVals[6][4][34] = 0.125538;
		vBinVals[6][4][35] = 0.125538;
		vBinVals[6][4][36] = 0.125538;
		vBinVals[6][4][37] = 0.255666;
		vBinVals[6][4][38] = 0.255666;
		vBinVals[6][4][39] = 0.255666;
		vBinVals[6][4][40] = 0.255666;
		vBinVals[6][4][41] = 0.255666;
		vBinVals[6][4][42] = 0.142676;
		vBinVals[6][4][43] = 0.142676;
		vBinVals[6][4][44] = 0.142676;
		vBinVals[6][4][45] = 0.142676;
		vBinVals[6][4][46] = 0.142676;
		vBinVals[6][4][47] = 0.186984;
		vBinVals[6][4][48] = 0.186984;
		vBinVals[6][4][49] = 0.186984;
		vBinVals[6][4][50] = 0.186984;
		vBinVals[6][4][51] = 0.186984;
		vBinVals[6][4][52] = 0.248957;
		vBinVals[6][4][53] = 0.248957;
		vBinVals[6][4][54] = 0.248957;
		vBinVals[6][4][55] = 0.248957;
		vBinVals[6][4][56] = 0.248957;
		vBinVals[6][4][57] = 0.217694;
		vBinVals[6][4][58] = 0.217694;
		vBinVals[6][4][59] = 0.217694;
		vBinVals[7][0][0] = -0.201208;
		vBinVals[7][0][1] = -0.201208;
		vBinVals[7][0][2] = -0.201208;
		vBinVals[7][0][3] = -0.201208;
		vBinVals[7][0][4] = -0.201208;
		vBinVals[7][0][5] = -0.201208;
		vBinVals[7][0][6] = -0.201208;
		vBinVals[7][0][7] = -0.201208;
		vBinVals[7][0][8] = -0.201208;
		vBinVals[7][0][9] = -0.201208;
		vBinVals[7][0][10] = -0.201208;
		vBinVals[7][0][11] = -0.201208;
		vBinVals[7][0][12] = -0.0901321;
		vBinVals[7][0][13] = -0.0901321;
		vBinVals[7][0][14] = -0.0901321;
		vBinVals[7][0][15] = -0.0901321;
		vBinVals[7][0][16] = -0.0901321;
		vBinVals[7][0][17] = 0.0162618;
		vBinVals[7][0][18] = 0.0162618;
		vBinVals[7][0][19] = 0.0162618;
		vBinVals[7][0][20] = 0.0162618;
		vBinVals[7][0][21] = 0.0162618;
		vBinVals[7][0][22] = 0.00509443;
		vBinVals[7][0][23] = 0.00509443;
		vBinVals[7][0][24] = 0.00509443;
		vBinVals[7][0][25] = 0.00509443;
		vBinVals[7][0][26] = 0.00509443;
		vBinVals[7][0][27] = 0.000401153;
		vBinVals[7][0][28] = 0.000401153;
		vBinVals[7][0][29] = 0.000401153;
		vBinVals[7][0][30] = 0.000401153;
		vBinVals[7][0][31] = 0.000401153;
		vBinVals[7][0][32] = 0.0235141;
		vBinVals[7][0][33] = 0.0235141;
		vBinVals[7][0][34] = 0.0235141;
		vBinVals[7][0][35] = 0.0235141;
		vBinVals[7][0][36] = 0.0235141;
		vBinVals[7][0][37] = 0.00113626;
		vBinVals[7][0][38] = 0.00113626;
		vBinVals[7][0][39] = 0.00113626;
		vBinVals[7][0][40] = 0.00113626;
		vBinVals[7][0][41] = 0.00113626;
		vBinVals[7][0][42] = -0.0176691;
		vBinVals[7][0][43] = -0.0176691;
		vBinVals[7][0][44] = -0.0176691;
		vBinVals[7][0][45] = -0.0176691;
		vBinVals[7][0][46] = -0.0176691;
		vBinVals[7][0][47] = 0.026015;
		vBinVals[7][0][48] = 0.026015;
		vBinVals[7][0][49] = 0.026015;
		vBinVals[7][0][50] = 0.026015;
		vBinVals[7][0][51] = 0.026015;
		vBinVals[7][0][52] = 0.0321351;
		vBinVals[7][0][53] = 0.0321351;
		vBinVals[7][0][54] = 0.0321351;
		vBinVals[7][0][55] = 0.0321351;
		vBinVals[7][0][56] = 0.0321351;
		vBinVals[7][0][57] = 0.0357106;
		vBinVals[7][0][58] = 0.0357106;
		vBinVals[7][0][59] = 0.0357106;
		vBinVals[7][1][0] = 0.154692;
		vBinVals[7][1][1] = 0.154692;
		vBinVals[7][1][2] = 0.154692;
		vBinVals[7][1][3] = 0.154692;
		vBinVals[7][1][4] = 0.154692;
		vBinVals[7][1][5] = 0.154692;
		vBinVals[7][1][6] = 0.154692;
		vBinVals[7][1][7] = 0.154692;
		vBinVals[7][1][8] = 0.154692;
		vBinVals[7][1][9] = 0.154692;
		vBinVals[7][1][10] = 0.154692;
		vBinVals[7][1][11] = 0.154692;
		vBinVals[7][1][12] = 0.180484;
		vBinVals[7][1][13] = 0.180484;
		vBinVals[7][1][14] = 0.180484;
		vBinVals[7][1][15] = 0.180484;
		vBinVals[7][1][16] = 0.180484;
		vBinVals[7][1][17] = 0.20722;
		vBinVals[7][1][18] = 0.20722;
		vBinVals[7][1][19] = 0.20722;
		vBinVals[7][1][20] = 0.20722;
		vBinVals[7][1][21] = 0.20722;
		vBinVals[7][1][22] = 0.185284;
		vBinVals[7][1][23] = 0.185284;
		vBinVals[7][1][24] = 0.185284;
		vBinVals[7][1][25] = 0.185284;
		vBinVals[7][1][26] = 0.185284;
		vBinVals[7][1][27] = 0.18126;
		vBinVals[7][1][28] = 0.18126;
		vBinVals[7][1][29] = 0.18126;
		vBinVals[7][1][30] = 0.18126;
		vBinVals[7][1][31] = 0.18126;
		vBinVals[7][1][32] = 0.174376;
		vBinVals[7][1][33] = 0.174376;
		vBinVals[7][1][34] = 0.174376;
		vBinVals[7][1][35] = 0.174376;
		vBinVals[7][1][36] = 0.174376;
		vBinVals[7][1][37] = 0.165648;
		vBinVals[7][1][38] = 0.165648;
		vBinVals[7][1][39] = 0.165648;
		vBinVals[7][1][40] = 0.165648;
		vBinVals[7][1][41] = 0.165648;
		vBinVals[7][1][42] = 0.136938;
		vBinVals[7][1][43] = 0.136938;
		vBinVals[7][1][44] = 0.136938;
		vBinVals[7][1][45] = 0.136938;
		vBinVals[7][1][46] = 0.136938;
		vBinVals[7][1][47] = 0.143576;
		vBinVals[7][1][48] = 0.143576;
		vBinVals[7][1][49] = 0.143576;
		vBinVals[7][1][50] = 0.143576;
		vBinVals[7][1][51] = 0.143576;
		vBinVals[7][1][52] = 0.136395;
		vBinVals[7][1][53] = 0.136395;
		vBinVals[7][1][54] = 0.136395;
		vBinVals[7][1][55] = 0.136395;
		vBinVals[7][1][56] = 0.136395;
		vBinVals[7][1][57] = 0.135764;
		vBinVals[7][1][58] = 0.135764;
		vBinVals[7][1][59] = 0.135764;
		vBinVals[7][2][0] = -0.313959;
		vBinVals[7][2][1] = -0.313959;
		vBinVals[7][2][2] = -0.313959;
		vBinVals[7][2][3] = -0.313959;
		vBinVals[7][2][4] = -0.313959;
		vBinVals[7][2][5] = -0.313959;
		vBinVals[7][2][6] = -0.313959;
		vBinVals[7][2][7] = -0.313959;
		vBinVals[7][2][8] = -0.313959;
		vBinVals[7][2][9] = -0.313959;
		vBinVals[7][2][10] = -0.313959;
		vBinVals[7][2][11] = -0.313959;
		vBinVals[7][2][12] = -0.367048;
		vBinVals[7][2][13] = -0.367048;
		vBinVals[7][2][14] = -0.367048;
		vBinVals[7][2][15] = -0.367048;
		vBinVals[7][2][16] = -0.367048;
		vBinVals[7][2][17] = -0.174408;
		vBinVals[7][2][18] = -0.174408;
		vBinVals[7][2][19] = -0.174408;
		vBinVals[7][2][20] = -0.174408;
		vBinVals[7][2][21] = -0.174408;
		vBinVals[7][2][22] = -0.118018;
		vBinVals[7][2][23] = -0.118018;
		vBinVals[7][2][24] = -0.118018;
		vBinVals[7][2][25] = -0.118018;
		vBinVals[7][2][26] = -0.118018;
		vBinVals[7][2][27] = -0.101365;
		vBinVals[7][2][28] = -0.101365;
		vBinVals[7][2][29] = -0.101365;
		vBinVals[7][2][30] = -0.101365;
		vBinVals[7][2][31] = -0.101365;
		vBinVals[7][2][32] = -0.117318;
		vBinVals[7][2][33] = -0.117318;
		vBinVals[7][2][34] = -0.117318;
		vBinVals[7][2][35] = -0.117318;
		vBinVals[7][2][36] = -0.117318;
		vBinVals[7][2][37] = -0.100033;
		vBinVals[7][2][38] = -0.100033;
		vBinVals[7][2][39] = -0.100033;
		vBinVals[7][2][40] = -0.100033;
		vBinVals[7][2][41] = -0.100033;
		vBinVals[7][2][42] = -0.0782654;
		vBinVals[7][2][43] = -0.0782654;
		vBinVals[7][2][44] = -0.0782654;
		vBinVals[7][2][45] = -0.0782654;
		vBinVals[7][2][46] = -0.0782654;
		vBinVals[7][2][47] = -0.0794711;
		vBinVals[7][2][48] = -0.0794711;
		vBinVals[7][2][49] = -0.0794711;
		vBinVals[7][2][50] = -0.0794711;
		vBinVals[7][2][51] = -0.0794711;
		vBinVals[7][2][52] = -0.0922871;
		vBinVals[7][2][53] = -0.0922871;
		vBinVals[7][2][54] = -0.0922871;
		vBinVals[7][2][55] = -0.0922871;
		vBinVals[7][2][56] = -0.0922871;
		vBinVals[7][2][57] = -0.0453136;
		vBinVals[7][2][58] = -0.0453136;
		vBinVals[7][2][59] = -0.0453136;
		vBinVals[7][3][0] = 0.131532;
		vBinVals[7][3][1] = 0.131532;
		vBinVals[7][3][2] = 0.131532;
		vBinVals[7][3][3] = 0.131532;
		vBinVals[7][3][4] = 0.131532;
		vBinVals[7][3][5] = 0.131532;
		vBinVals[7][3][6] = 0.131532;
		vBinVals[7][3][7] = 0.131532;
		vBinVals[7][3][8] = 0.131532;
		vBinVals[7][3][9] = 0.131532;
		vBinVals[7][3][10] = 0.131532;
		vBinVals[7][3][11] = 0.131532;
		vBinVals[7][3][12] = -0.0216764;
		vBinVals[7][3][13] = -0.0216764;
		vBinVals[7][3][14] = -0.0216764;
		vBinVals[7][3][15] = -0.0216764;
		vBinVals[7][3][16] = -0.0216764;
		vBinVals[7][3][17] = 0.101255;
		vBinVals[7][3][18] = 0.101255;
		vBinVals[7][3][19] = 0.101255;
		vBinVals[7][3][20] = 0.101255;
		vBinVals[7][3][21] = 0.101255;
		vBinVals[7][3][22] = 0.114489;
		vBinVals[7][3][23] = 0.114489;
		vBinVals[7][3][24] = 0.114489;
		vBinVals[7][3][25] = 0.114489;
		vBinVals[7][3][26] = 0.114489;
		vBinVals[7][3][27] = 0.112235;
		vBinVals[7][3][28] = 0.112235;
		vBinVals[7][3][29] = 0.112235;
		vBinVals[7][3][30] = 0.112235;
		vBinVals[7][3][31] = 0.112235;
		vBinVals[7][3][32] = 0.0911959;
		vBinVals[7][3][33] = 0.0911959;
		vBinVals[7][3][34] = 0.0911959;
		vBinVals[7][3][35] = 0.0911959;
		vBinVals[7][3][36] = 0.0911959;
		vBinVals[7][3][37] = 0.0740797;
		vBinVals[7][3][38] = 0.0740797;
		vBinVals[7][3][39] = 0.0740797;
		vBinVals[7][3][40] = 0.0740797;
		vBinVals[7][3][41] = 0.0740797;
		vBinVals[7][3][42] = 0.0815081;
		vBinVals[7][3][43] = 0.0815081;
		vBinVals[7][3][44] = 0.0815081;
		vBinVals[7][3][45] = 0.0815081;
		vBinVals[7][3][46] = 0.0815081;
		vBinVals[7][3][47] = 0.0747259;
		vBinVals[7][3][48] = 0.0747259;
		vBinVals[7][3][49] = 0.0747259;
		vBinVals[7][3][50] = 0.0747259;
		vBinVals[7][3][51] = 0.0747259;
		vBinVals[7][3][52] = 0.068659;
		vBinVals[7][3][53] = 0.068659;
		vBinVals[7][3][54] = 0.068659;
		vBinVals[7][3][55] = 0.068659;
		vBinVals[7][3][56] = 0.068659;
		vBinVals[7][3][57] = 0.0724492;
		vBinVals[7][3][58] = 0.0724492;
		vBinVals[7][3][59] = 0.0724492;
		vBinVals[7][4][0] = 2.29338;
		vBinVals[7][4][1] = 2.29338;
		vBinVals[7][4][2] = 2.29338;
		vBinVals[7][4][3] = 2.29338;
		vBinVals[7][4][4] = 2.29338;
		vBinVals[7][4][5] = 2.29338;
		vBinVals[7][4][6] = 2.29338;
		vBinVals[7][4][7] = 2.29338;
		vBinVals[7][4][8] = 2.29338;
		vBinVals[7][4][9] = 2.29338;
		vBinVals[7][4][10] = 2.29338;
		vBinVals[7][4][11] = 2.29338;
		vBinVals[7][4][12] = 0.707677;
		vBinVals[7][4][13] = 0.707677;
		vBinVals[7][4][14] = 0.707677;
		vBinVals[7][4][15] = 0.707677;
		vBinVals[7][4][16] = 0.707677;
		vBinVals[7][4][17] = 0.21677;
		vBinVals[7][4][18] = 0.21677;
		vBinVals[7][4][19] = 0.21677;
		vBinVals[7][4][20] = 0.21677;
		vBinVals[7][4][21] = 0.21677;
		vBinVals[7][4][22] = 0.171769;
		vBinVals[7][4][23] = 0.171769;
		vBinVals[7][4][24] = 0.171769;
		vBinVals[7][4][25] = 0.171769;
		vBinVals[7][4][26] = 0.171769;
		vBinVals[7][4][27] = 0.187562;
		vBinVals[7][4][28] = 0.187562;
		vBinVals[7][4][29] = 0.187562;
		vBinVals[7][4][30] = 0.187562;
		vBinVals[7][4][31] = 0.187562;
		vBinVals[7][4][32] = 0.171954;
		vBinVals[7][4][33] = 0.171954;
		vBinVals[7][4][34] = 0.171954;
		vBinVals[7][4][35] = 0.171954;
		vBinVals[7][4][36] = 0.171954;
		vBinVals[7][4][37] = 0.191595;
		vBinVals[7][4][38] = 0.191595;
		vBinVals[7][4][39] = 0.191595;
		vBinVals[7][4][40] = 0.191595;
		vBinVals[7][4][41] = 0.191595;
		vBinVals[7][4][42] = 0.200133;
		vBinVals[7][4][43] = 0.200133;
		vBinVals[7][4][44] = 0.200133;
		vBinVals[7][4][45] = 0.200133;
		vBinVals[7][4][46] = 0.200133;
		vBinVals[7][4][47] = 0.241419;
		vBinVals[7][4][48] = 0.241419;
		vBinVals[7][4][49] = 0.241419;
		vBinVals[7][4][50] = 0.241419;
		vBinVals[7][4][51] = 0.241419;
		vBinVals[7][4][52] = 0.21665;
		vBinVals[7][4][53] = 0.21665;
		vBinVals[7][4][54] = 0.21665;
		vBinVals[7][4][55] = 0.21665;
		vBinVals[7][4][56] = 0.21665;
		vBinVals[7][4][57] = 0.237653;
		vBinVals[7][4][58] = 0.237653;
		vBinVals[7][4][59] = 0.237653;
		vBinVals[8][0][0] = -0.217531;
		vBinVals[8][0][1] = -0.217531;
		vBinVals[8][0][2] = -0.217531;
		vBinVals[8][0][3] = -0.217531;
		vBinVals[8][0][4] = -0.217531;
		vBinVals[8][0][5] = -0.217531;
		vBinVals[8][0][6] = -0.217531;
		vBinVals[8][0][7] = -0.217531;
		vBinVals[8][0][8] = -0.217531;
		vBinVals[8][0][9] = -0.217531;
		vBinVals[8][0][10] = -0.217531;
		vBinVals[8][0][11] = -0.217531;
		vBinVals[8][0][12] = -0.118869;
		vBinVals[8][0][13] = -0.118869;
		vBinVals[8][0][14] = -0.118869;
		vBinVals[8][0][15] = -0.118869;
		vBinVals[8][0][16] = -0.118869;
		vBinVals[8][0][17] = -0.0244903;
		vBinVals[8][0][18] = -0.0244903;
		vBinVals[8][0][19] = -0.0244903;
		vBinVals[8][0][20] = -0.0244903;
		vBinVals[8][0][21] = -0.0244903;
		vBinVals[8][0][22] = 0.00469728;
		vBinVals[8][0][23] = 0.00469728;
		vBinVals[8][0][24] = 0.00469728;
		vBinVals[8][0][25] = 0.00469728;
		vBinVals[8][0][26] = 0.00469728;
		vBinVals[8][0][27] = -0.000912776;
		vBinVals[8][0][28] = -0.000912776;
		vBinVals[8][0][29] = -0.000912776;
		vBinVals[8][0][30] = -0.000912776;
		vBinVals[8][0][31] = -0.000912776;
		vBinVals[8][0][32] = -0.0113414;
		vBinVals[8][0][33] = -0.0113414;
		vBinVals[8][0][34] = -0.0113414;
		vBinVals[8][0][35] = -0.0113414;
		vBinVals[8][0][36] = -0.0113414;
		vBinVals[8][0][37] = -0.0331641;
		vBinVals[8][0][38] = -0.0331641;
		vBinVals[8][0][39] = -0.0331641;
		vBinVals[8][0][40] = -0.0331641;
		vBinVals[8][0][41] = -0.0331641;
		vBinVals[8][0][42] = -0.00313909;
		vBinVals[8][0][43] = -0.00313909;
		vBinVals[8][0][44] = -0.00313909;
		vBinVals[8][0][45] = -0.00313909;
		vBinVals[8][0][46] = -0.00313909;
		vBinVals[8][0][47] = -0.0112428;
		vBinVals[8][0][48] = -0.0112428;
		vBinVals[8][0][49] = -0.0112428;
		vBinVals[8][0][50] = -0.0112428;
		vBinVals[8][0][51] = -0.0112428;
		vBinVals[8][0][52] = -0.0352037;
		vBinVals[8][0][53] = -0.0352037;
		vBinVals[8][0][54] = -0.0352037;
		vBinVals[8][0][55] = -0.0352037;
		vBinVals[8][0][56] = -0.0352037;
		vBinVals[8][0][57] = 0.011293;
		vBinVals[8][0][58] = 0.011293;
		vBinVals[8][0][59] = 0.011293;
		vBinVals[8][1][0] = 0.131491;
		vBinVals[8][1][1] = 0.131491;
		vBinVals[8][1][2] = 0.131491;
		vBinVals[8][1][3] = 0.131491;
		vBinVals[8][1][4] = 0.131491;
		vBinVals[8][1][5] = 0.131491;
		vBinVals[8][1][6] = 0.131491;
		vBinVals[8][1][7] = 0.131491;
		vBinVals[8][1][8] = 0.131491;
		vBinVals[8][1][9] = 0.131491;
		vBinVals[8][1][10] = 0.131491;
		vBinVals[8][1][11] = 0.131491;
		vBinVals[8][1][12] = 0.167883;
		vBinVals[8][1][13] = 0.167883;
		vBinVals[8][1][14] = 0.167883;
		vBinVals[8][1][15] = 0.167883;
		vBinVals[8][1][16] = 0.167883;
		vBinVals[8][1][17] = 0.189814;
		vBinVals[8][1][18] = 0.189814;
		vBinVals[8][1][19] = 0.189814;
		vBinVals[8][1][20] = 0.189814;
		vBinVals[8][1][21] = 0.189814;
		vBinVals[8][1][22] = 0.192968;
		vBinVals[8][1][23] = 0.192968;
		vBinVals[8][1][24] = 0.192968;
		vBinVals[8][1][25] = 0.192968;
		vBinVals[8][1][26] = 0.192968;
		vBinVals[8][1][27] = 0.180559;
		vBinVals[8][1][28] = 0.180559;
		vBinVals[8][1][29] = 0.180559;
		vBinVals[8][1][30] = 0.180559;
		vBinVals[8][1][31] = 0.180559;
		vBinVals[8][1][32] = 0.176323;
		vBinVals[8][1][33] = 0.176323;
		vBinVals[8][1][34] = 0.176323;
		vBinVals[8][1][35] = 0.176323;
		vBinVals[8][1][36] = 0.176323;
		vBinVals[8][1][37] = 0.154606;
		vBinVals[8][1][38] = 0.154606;
		vBinVals[8][1][39] = 0.154606;
		vBinVals[8][1][40] = 0.154606;
		vBinVals[8][1][41] = 0.154606;
		vBinVals[8][1][42] = 0.164586;
		vBinVals[8][1][43] = 0.164586;
		vBinVals[8][1][44] = 0.164586;
		vBinVals[8][1][45] = 0.164586;
		vBinVals[8][1][46] = 0.164586;
		vBinVals[8][1][47] = 0.149646;
		vBinVals[8][1][48] = 0.149646;
		vBinVals[8][1][49] = 0.149646;
		vBinVals[8][1][50] = 0.149646;
		vBinVals[8][1][51] = 0.149646;
		vBinVals[8][1][52] = 0.120553;
		vBinVals[8][1][53] = 0.120553;
		vBinVals[8][1][54] = 0.120553;
		vBinVals[8][1][55] = 0.120553;
		vBinVals[8][1][56] = 0.120553;
		vBinVals[8][1][57] = 0.137443;
		vBinVals[8][1][58] = 0.137443;
		vBinVals[8][1][59] = 0.137443;
		vBinVals[8][2][0] = -0.434468;
		vBinVals[8][2][1] = -0.434468;
		vBinVals[8][2][2] = -0.434468;
		vBinVals[8][2][3] = -0.434468;
		vBinVals[8][2][4] = -0.434468;
		vBinVals[8][2][5] = -0.434468;
		vBinVals[8][2][6] = -0.434468;
		vBinVals[8][2][7] = -0.434468;
		vBinVals[8][2][8] = -0.434468;
		vBinVals[8][2][9] = -0.434468;
		vBinVals[8][2][10] = -0.434468;
		vBinVals[8][2][11] = -0.434468;
		vBinVals[8][2][12] = -0.213752;
		vBinVals[8][2][13] = -0.213752;
		vBinVals[8][2][14] = -0.213752;
		vBinVals[8][2][15] = -0.213752;
		vBinVals[8][2][16] = -0.213752;
		vBinVals[8][2][17] = -0.183936;
		vBinVals[8][2][18] = -0.183936;
		vBinVals[8][2][19] = -0.183936;
		vBinVals[8][2][20] = -0.183936;
		vBinVals[8][2][21] = -0.183936;
		vBinVals[8][2][22] = -0.148666;
		vBinVals[8][2][23] = -0.148666;
		vBinVals[8][2][24] = -0.148666;
		vBinVals[8][2][25] = -0.148666;
		vBinVals[8][2][26] = -0.148666;
		vBinVals[8][2][27] = -0.134934;
		vBinVals[8][2][28] = -0.134934;
		vBinVals[8][2][29] = -0.134934;
		vBinVals[8][2][30] = -0.134934;
		vBinVals[8][2][31] = -0.134934;
		vBinVals[8][2][32] = -0.131002;
		vBinVals[8][2][33] = -0.131002;
		vBinVals[8][2][34] = -0.131002;
		vBinVals[8][2][35] = -0.131002;
		vBinVals[8][2][36] = -0.131002;
		vBinVals[8][2][37] = -0.0934355;
		vBinVals[8][2][38] = -0.0934355;
		vBinVals[8][2][39] = -0.0934355;
		vBinVals[8][2][40] = -0.0934355;
		vBinVals[8][2][41] = -0.0934355;
		vBinVals[8][2][42] = -0.0944738;
		vBinVals[8][2][43] = -0.0944738;
		vBinVals[8][2][44] = -0.0944738;
		vBinVals[8][2][45] = -0.0944738;
		vBinVals[8][2][46] = -0.0944738;
		vBinVals[8][2][47] = -0.0862545;
		vBinVals[8][2][48] = -0.0862545;
		vBinVals[8][2][49] = -0.0862545;
		vBinVals[8][2][50] = -0.0862545;
		vBinVals[8][2][51] = -0.0862545;
		vBinVals[8][2][52] = -0.0834524;
		vBinVals[8][2][53] = -0.0834524;
		vBinVals[8][2][54] = -0.0834524;
		vBinVals[8][2][55] = -0.0834524;
		vBinVals[8][2][56] = -0.0834524;
		vBinVals[8][2][57] = -0.0933933;
		vBinVals[8][2][58] = -0.0933933;
		vBinVals[8][2][59] = -0.0933933;
		vBinVals[8][3][0] = 0.110948;
		vBinVals[8][3][1] = 0.110948;
		vBinVals[8][3][2] = 0.110948;
		vBinVals[8][3][3] = 0.110948;
		vBinVals[8][3][4] = 0.110948;
		vBinVals[8][3][5] = 0.110948;
		vBinVals[8][3][6] = 0.110948;
		vBinVals[8][3][7] = 0.110948;
		vBinVals[8][3][8] = 0.110948;
		vBinVals[8][3][9] = 0.110948;
		vBinVals[8][3][10] = 0.110948;
		vBinVals[8][3][11] = 0.110948;
		vBinVals[8][3][12] = 0.131853;
		vBinVals[8][3][13] = 0.131853;
		vBinVals[8][3][14] = 0.131853;
		vBinVals[8][3][15] = 0.131853;
		vBinVals[8][3][16] = 0.131853;
		vBinVals[8][3][17] = 0.109638;
		vBinVals[8][3][18] = 0.109638;
		vBinVals[8][3][19] = 0.109638;
		vBinVals[8][3][20] = 0.109638;
		vBinVals[8][3][21] = 0.109638;
		vBinVals[8][3][22] = 0.107964;
		vBinVals[8][3][23] = 0.107964;
		vBinVals[8][3][24] = 0.107964;
		vBinVals[8][3][25] = 0.107964;
		vBinVals[8][3][26] = 0.107964;
		vBinVals[8][3][27] = 0.102849;
		vBinVals[8][3][28] = 0.102849;
		vBinVals[8][3][29] = 0.102849;
		vBinVals[8][3][30] = 0.102849;
		vBinVals[8][3][31] = 0.102849;
		vBinVals[8][3][32] = 0.0822329;
		vBinVals[8][3][33] = 0.0822329;
		vBinVals[8][3][34] = 0.0822329;
		vBinVals[8][3][35] = 0.0822329;
		vBinVals[8][3][36] = 0.0822329;
		vBinVals[8][3][37] = 0.0864036;
		vBinVals[8][3][38] = 0.0864036;
		vBinVals[8][3][39] = 0.0864036;
		vBinVals[8][3][40] = 0.0864036;
		vBinVals[8][3][41] = 0.0864036;
		vBinVals[8][3][42] = 0.088691;
		vBinVals[8][3][43] = 0.088691;
		vBinVals[8][3][44] = 0.088691;
		vBinVals[8][3][45] = 0.088691;
		vBinVals[8][3][46] = 0.088691;
		vBinVals[8][3][47] = 0.0790669;
		vBinVals[8][3][48] = 0.0790669;
		vBinVals[8][3][49] = 0.0790669;
		vBinVals[8][3][50] = 0.0790669;
		vBinVals[8][3][51] = 0.0790669;
		vBinVals[8][3][52] = 0.076533;
		vBinVals[8][3][53] = 0.076533;
		vBinVals[8][3][54] = 0.076533;
		vBinVals[8][3][55] = 0.076533;
		vBinVals[8][3][56] = 0.076533;
		vBinVals[8][3][57] = 0.0634286;
		vBinVals[8][3][58] = 0.0634286;
		vBinVals[8][3][59] = 0.0634286;
		vBinVals[8][4][0] = 3.42137;
		vBinVals[8][4][1] = 3.42137;
		vBinVals[8][4][2] = 3.42137;
		vBinVals[8][4][3] = 3.42137;
		vBinVals[8][4][4] = 3.42137;
		vBinVals[8][4][5] = 3.42137;
		vBinVals[8][4][6] = 3.42137;
		vBinVals[8][4][7] = 3.42137;
		vBinVals[8][4][8] = 3.42137;
		vBinVals[8][4][9] = 3.42137;
		vBinVals[8][4][10] = 3.42137;
		vBinVals[8][4][11] = 3.42137;
		vBinVals[8][4][12] = 0.901789;
		vBinVals[8][4][13] = 0.901789;
		vBinVals[8][4][14] = 0.901789;
		vBinVals[8][4][15] = 0.901789;
		vBinVals[8][4][16] = 0.901789;
		vBinVals[8][4][17] = 0.330516;
		vBinVals[8][4][18] = 0.330516;
		vBinVals[8][4][19] = 0.330516;
		vBinVals[8][4][20] = 0.330516;
		vBinVals[8][4][21] = 0.330516;
		vBinVals[8][4][22] = 0.18878;
		vBinVals[8][4][23] = 0.18878;
		vBinVals[8][4][24] = 0.18878;
		vBinVals[8][4][25] = 0.18878;
		vBinVals[8][4][26] = 0.18878;
		vBinVals[8][4][27] = 0.192068;
		vBinVals[8][4][28] = 0.192068;
		vBinVals[8][4][29] = 0.192068;
		vBinVals[8][4][30] = 0.192068;
		vBinVals[8][4][31] = 0.192068;
		vBinVals[8][4][32] = 0.158671;
		vBinVals[8][4][33] = 0.158671;
		vBinVals[8][4][34] = 0.158671;
		vBinVals[8][4][35] = 0.158671;
		vBinVals[8][4][36] = 0.158671;
		vBinVals[8][4][37] = 0.204279;
		vBinVals[8][4][38] = 0.204279;
		vBinVals[8][4][39] = 0.204279;
		vBinVals[8][4][40] = 0.204279;
		vBinVals[8][4][41] = 0.204279;
		vBinVals[8][4][42] = 0.13829;
		vBinVals[8][4][43] = 0.13829;
		vBinVals[8][4][44] = 0.13829;
		vBinVals[8][4][45] = 0.13829;
		vBinVals[8][4][46] = 0.13829;
		vBinVals[8][4][47] = 0.214855;
		vBinVals[8][4][48] = 0.214855;
		vBinVals[8][4][49] = 0.214855;
		vBinVals[8][4][50] = 0.214855;
		vBinVals[8][4][51] = 0.214855;
		vBinVals[8][4][52] = 0.269977;
		vBinVals[8][4][53] = 0.269977;
		vBinVals[8][4][54] = 0.269977;
		vBinVals[8][4][55] = 0.269977;
		vBinVals[8][4][56] = 0.269977;
		vBinVals[8][4][57] = 0.162784;
		vBinVals[8][4][58] = 0.162784;
		vBinVals[8][4][59] = 0.162784;
		vBinVals[9][0][0] = -0.248872;
		vBinVals[9][0][1] = -0.248872;
		vBinVals[9][0][2] = -0.248872;
		vBinVals[9][0][3] = -0.248872;
		vBinVals[9][0][4] = -0.248872;
		vBinVals[9][0][5] = -0.248872;
		vBinVals[9][0][6] = -0.248872;
		vBinVals[9][0][7] = -0.248872;
		vBinVals[9][0][8] = -0.248872;
		vBinVals[9][0][9] = -0.248872;
		vBinVals[9][0][10] = -0.248872;
		vBinVals[9][0][11] = -0.248872;
		vBinVals[9][0][12] = -3024.11;
		vBinVals[9][0][13] = -3024.11;
		vBinVals[9][0][14] = -3024.11;
		vBinVals[9][0][15] = -3024.11;
		vBinVals[9][0][16] = -3024.11;
		vBinVals[9][0][17] = -0.0662623;
		vBinVals[9][0][18] = -0.0662623;
		vBinVals[9][0][19] = -0.0662623;
		vBinVals[9][0][20] = -0.0662623;
		vBinVals[9][0][21] = -0.0662623;
		vBinVals[9][0][22] = -0.108331;
		vBinVals[9][0][23] = -0.108331;
		vBinVals[9][0][24] = -0.108331;
		vBinVals[9][0][25] = -0.108331;
		vBinVals[9][0][26] = -0.108331;
		vBinVals[9][0][27] = -0.00671763;
		vBinVals[9][0][28] = -0.00671763;
		vBinVals[9][0][29] = -0.00671763;
		vBinVals[9][0][30] = -0.00671763;
		vBinVals[9][0][31] = -0.00671763;
		vBinVals[9][0][32] = 0.0153384;
		vBinVals[9][0][33] = 0.0153384;
		vBinVals[9][0][34] = 0.0153384;
		vBinVals[9][0][35] = 0.0153384;
		vBinVals[9][0][36] = 0.0153384;
		vBinVals[9][0][37] = -0.0381201;
		vBinVals[9][0][38] = -0.0381201;
		vBinVals[9][0][39] = -0.0381201;
		vBinVals[9][0][40] = -0.0381201;
		vBinVals[9][0][41] = -0.0381201;
		vBinVals[9][0][42] = -0.0476759;
		vBinVals[9][0][43] = -0.0476759;
		vBinVals[9][0][44] = -0.0476759;
		vBinVals[9][0][45] = -0.0476759;
		vBinVals[9][0][46] = -0.0476759;
		vBinVals[9][0][47] = -0.0432301;
		vBinVals[9][0][48] = -0.0432301;
		vBinVals[9][0][49] = -0.0432301;
		vBinVals[9][0][50] = -0.0432301;
		vBinVals[9][0][51] = -0.0432301;
		vBinVals[9][0][52] = -0.0181688;
		vBinVals[9][0][53] = -0.0181688;
		vBinVals[9][0][54] = -0.0181688;
		vBinVals[9][0][55] = -0.0181688;
		vBinVals[9][0][56] = -0.0181688;
		vBinVals[9][0][57] = -0.039806;
		vBinVals[9][0][58] = -0.039806;
		vBinVals[9][0][59] = -0.039806;
		vBinVals[9][1][0] = 0.11955;
		vBinVals[9][1][1] = 0.11955;
		vBinVals[9][1][2] = 0.11955;
		vBinVals[9][1][3] = 0.11955;
		vBinVals[9][1][4] = 0.11955;
		vBinVals[9][1][5] = 0.11955;
		vBinVals[9][1][6] = 0.11955;
		vBinVals[9][1][7] = 0.11955;
		vBinVals[9][1][8] = 0.11955;
		vBinVals[9][1][9] = 0.11955;
		vBinVals[9][1][10] = 0.11955;
		vBinVals[9][1][11] = 0.11955;
		vBinVals[9][1][12] = 123481;
		vBinVals[9][1][13] = 123481;
		vBinVals[9][1][14] = 123481;
		vBinVals[9][1][15] = 123481;
		vBinVals[9][1][16] = 123481;
		vBinVals[9][1][17] = 0.17631;
		vBinVals[9][1][18] = 0.17631;
		vBinVals[9][1][19] = 0.17631;
		vBinVals[9][1][20] = 0.17631;
		vBinVals[9][1][21] = 0.17631;
		vBinVals[9][1][22] = 0.147023;
		vBinVals[9][1][23] = 0.147023;
		vBinVals[9][1][24] = 0.147023;
		vBinVals[9][1][25] = 0.147023;
		vBinVals[9][1][26] = 0.147023;
		vBinVals[9][1][27] = 0.186368;
		vBinVals[9][1][28] = 0.186368;
		vBinVals[9][1][29] = 0.186368;
		vBinVals[9][1][30] = 0.186368;
		vBinVals[9][1][31] = 0.186368;
		vBinVals[9][1][32] = 0.193102;
		vBinVals[9][1][33] = 0.193102;
		vBinVals[9][1][34] = 0.193102;
		vBinVals[9][1][35] = 0.193102;
		vBinVals[9][1][36] = 0.193102;
		vBinVals[9][1][37] = 0.150848;
		vBinVals[9][1][38] = 0.150848;
		vBinVals[9][1][39] = 0.150848;
		vBinVals[9][1][40] = 0.150848;
		vBinVals[9][1][41] = 0.150848;
		vBinVals[9][1][42] = 0.144115;
		vBinVals[9][1][43] = 0.144115;
		vBinVals[9][1][44] = 0.144115;
		vBinVals[9][1][45] = 0.144115;
		vBinVals[9][1][46] = 0.144115;
		vBinVals[9][1][47] = 0.170984;
		vBinVals[9][1][48] = 0.170984;
		vBinVals[9][1][49] = 0.170984;
		vBinVals[9][1][50] = 0.170984;
		vBinVals[9][1][51] = 0.170984;
		vBinVals[9][1][52] = 0.144327;
		vBinVals[9][1][53] = 0.144327;
		vBinVals[9][1][54] = 0.144327;
		vBinVals[9][1][55] = 0.144327;
		vBinVals[9][1][56] = 0.144327;
		vBinVals[9][1][57] = 0.121807;
		vBinVals[9][1][58] = 0.121807;
		vBinVals[9][1][59] = 0.121807;
		vBinVals[9][2][0] = -0.59061;
		vBinVals[9][2][1] = -0.59061;
		vBinVals[9][2][2] = -0.59061;
		vBinVals[9][2][3] = -0.59061;
		vBinVals[9][2][4] = -0.59061;
		vBinVals[9][2][5] = -0.59061;
		vBinVals[9][2][6] = -0.59061;
		vBinVals[9][2][7] = -0.59061;
		vBinVals[9][2][8] = -0.59061;
		vBinVals[9][2][9] = -0.59061;
		vBinVals[9][2][10] = -0.59061;
		vBinVals[9][2][11] = -0.59061;
		vBinVals[9][2][12] = -0.252286;
		vBinVals[9][2][13] = -0.252286;
		vBinVals[9][2][14] = -0.252286;
		vBinVals[9][2][15] = -0.252286;
		vBinVals[9][2][16] = -0.252286;
		vBinVals[9][2][17] = -0.218736;
		vBinVals[9][2][18] = -0.218736;
		vBinVals[9][2][19] = -0.218736;
		vBinVals[9][2][20] = -0.218736;
		vBinVals[9][2][21] = -0.218736;
		vBinVals[9][2][22] = 0.152886;
		vBinVals[9][2][23] = 0.152886;
		vBinVals[9][2][24] = 0.152886;
		vBinVals[9][2][25] = 0.152886;
		vBinVals[9][2][26] = 0.152886;
		vBinVals[9][2][27] = -0.164649;
		vBinVals[9][2][28] = -0.164649;
		vBinVals[9][2][29] = -0.164649;
		vBinVals[9][2][30] = -0.164649;
		vBinVals[9][2][31] = -0.164649;
		vBinVals[9][2][32] = -0.14354;
		vBinVals[9][2][33] = -0.14354;
		vBinVals[9][2][34] = -0.14354;
		vBinVals[9][2][35] = -0.14354;
		vBinVals[9][2][36] = -0.14354;
		vBinVals[9][2][37] = -0.111503;
		vBinVals[9][2][38] = -0.111503;
		vBinVals[9][2][39] = -0.111503;
		vBinVals[9][2][40] = -0.111503;
		vBinVals[9][2][41] = -0.111503;
		vBinVals[9][2][42] = -0.103231;
		vBinVals[9][2][43] = -0.103231;
		vBinVals[9][2][44] = -0.103231;
		vBinVals[9][2][45] = -0.103231;
		vBinVals[9][2][46] = -0.103231;
		vBinVals[9][2][47] = -0.0957714;
		vBinVals[9][2][48] = -0.0957714;
		vBinVals[9][2][49] = -0.0957714;
		vBinVals[9][2][50] = -0.0957714;
		vBinVals[9][2][51] = -0.0957714;
		vBinVals[9][2][52] = -0.114921;
		vBinVals[9][2][53] = -0.114921;
		vBinVals[9][2][54] = -0.114921;
		vBinVals[9][2][55] = -0.114921;
		vBinVals[9][2][56] = -0.114921;
		vBinVals[9][2][57] = -0.0867108;
		vBinVals[9][2][58] = -0.0867108;
		vBinVals[9][2][59] = -0.0867108;
		vBinVals[9][3][0] = 0.0961599;
		vBinVals[9][3][1] = 0.0961599;
		vBinVals[9][3][2] = 0.0961599;
		vBinVals[9][3][3] = 0.0961599;
		vBinVals[9][3][4] = 0.0961599;
		vBinVals[9][3][5] = 0.0961599;
		vBinVals[9][3][6] = 0.0961599;
		vBinVals[9][3][7] = 0.0961599;
		vBinVals[9][3][8] = 0.0961599;
		vBinVals[9][3][9] = 0.0961599;
		vBinVals[9][3][10] = 0.0961599;
		vBinVals[9][3][11] = 0.0961599;
		vBinVals[9][3][12] = 0.0882695;
		vBinVals[9][3][13] = 0.0882695;
		vBinVals[9][3][14] = 0.0882695;
		vBinVals[9][3][15] = 0.0882695;
		vBinVals[9][3][16] = 0.0882695;
		vBinVals[9][3][17] = 0.0994334;
		vBinVals[9][3][18] = 0.0994334;
		vBinVals[9][3][19] = 0.0994334;
		vBinVals[9][3][20] = 0.0994334;
		vBinVals[9][3][21] = 0.0994334;
		vBinVals[9][3][22] = 0.0598369;
		vBinVals[9][3][23] = 0.0598369;
		vBinVals[9][3][24] = 0.0598369;
		vBinVals[9][3][25] = 0.0598369;
		vBinVals[9][3][26] = 0.0598369;
		vBinVals[9][3][27] = 0.0929396;
		vBinVals[9][3][28] = 0.0929396;
		vBinVals[9][3][29] = 0.0929396;
		vBinVals[9][3][30] = 0.0929396;
		vBinVals[9][3][31] = 0.0929396;
		vBinVals[9][3][32] = 0.0849565;
		vBinVals[9][3][33] = 0.0849565;
		vBinVals[9][3][34] = 0.0849565;
		vBinVals[9][3][35] = 0.0849565;
		vBinVals[9][3][36] = 0.0849565;
		vBinVals[9][3][37] = 0.103973;
		vBinVals[9][3][38] = 0.103973;
		vBinVals[9][3][39] = 0.103973;
		vBinVals[9][3][40] = 0.103973;
		vBinVals[9][3][41] = 0.103973;
		vBinVals[9][3][42] = 0.0842764;
		vBinVals[9][3][43] = 0.0842764;
		vBinVals[9][3][44] = 0.0842764;
		vBinVals[9][3][45] = 0.0842764;
		vBinVals[9][3][46] = 0.0842764;
		vBinVals[9][3][47] = 0.0805004;
		vBinVals[9][3][48] = 0.0805004;
		vBinVals[9][3][49] = 0.0805004;
		vBinVals[9][3][50] = 0.0805004;
		vBinVals[9][3][51] = 0.0805004;
		vBinVals[9][3][52] = 0.0630821;
		vBinVals[9][3][53] = 0.0630821;
		vBinVals[9][3][54] = 0.0630821;
		vBinVals[9][3][55] = 0.0630821;
		vBinVals[9][3][56] = 0.0630821;
		vBinVals[9][3][57] = 0.0651845;
		vBinVals[9][3][58] = 0.0651845;
		vBinVals[9][3][59] = 0.0651845;
		vBinVals[9][4][0] = 3.53851;
		vBinVals[9][4][1] = 3.53851;
		vBinVals[9][4][2] = 3.53851;
		vBinVals[9][4][3] = 3.53851;
		vBinVals[9][4][4] = 3.53851;
		vBinVals[9][4][5] = 3.53851;
		vBinVals[9][4][6] = 3.53851;
		vBinVals[9][4][7] = 3.53851;
		vBinVals[9][4][8] = 3.53851;
		vBinVals[9][4][9] = 3.53851;
		vBinVals[9][4][10] = 3.53851;
		vBinVals[9][4][11] = 3.53851;
		vBinVals[9][4][12] = 0.0016409;
		vBinVals[9][4][13] = 0.0016409;
		vBinVals[9][4][14] = 0.0016409;
		vBinVals[9][4][15] = 0.0016409;
		vBinVals[9][4][16] = 0.0016409;
		vBinVals[9][4][17] = 0.450261;
		vBinVals[9][4][18] = 0.450261;
		vBinVals[9][4][19] = 0.450261;
		vBinVals[9][4][20] = 0.450261;
		vBinVals[9][4][21] = 0.450261;
		vBinVals[9][4][22] = 0.555343;
		vBinVals[9][4][23] = 0.555343;
		vBinVals[9][4][24] = 0.555343;
		vBinVals[9][4][25] = 0.555343;
		vBinVals[9][4][26] = 0.555343;
		vBinVals[9][4][27] = 0.204421;
		vBinVals[9][4][28] = 0.204421;
		vBinVals[9][4][29] = 0.204421;
		vBinVals[9][4][30] = 0.204421;
		vBinVals[9][4][31] = 0.204421;
		vBinVals[9][4][32] = 0.118837;
		vBinVals[9][4][33] = 0.118837;
		vBinVals[9][4][34] = 0.118837;
		vBinVals[9][4][35] = 0.118837;
		vBinVals[9][4][36] = 0.118837;
		vBinVals[9][4][37] = 0.251713;
		vBinVals[9][4][38] = 0.251713;
		vBinVals[9][4][39] = 0.251713;
		vBinVals[9][4][40] = 0.251713;
		vBinVals[9][4][41] = 0.251713;
		vBinVals[9][4][42] = 0.146899;
		vBinVals[9][4][43] = 0.146899;
		vBinVals[9][4][44] = 0.146899;
		vBinVals[9][4][45] = 0.146899;
		vBinVals[9][4][46] = 0.146899;
		vBinVals[9][4][47] = 0.0640807;
		vBinVals[9][4][48] = 0.0640807;
		vBinVals[9][4][49] = 0.0640807;
		vBinVals[9][4][50] = 0.0640807;
		vBinVals[9][4][51] = 0.0640807;
		vBinVals[9][4][52] = 0.234438;
		vBinVals[9][4][53] = 0.234438;
		vBinVals[9][4][54] = 0.234438;
		vBinVals[9][4][55] = 0.234438;
		vBinVals[9][4][56] = 0.234438;
		vBinVals[9][4][57] = 0.227191;
		vBinVals[9][4][58] = 0.227191;
		vBinVals[9][4][59] = 0.227191;
		vBinVals[10][0][0] = -0.214372;
		vBinVals[10][0][1] = -0.214372;
		vBinVals[10][0][2] = -0.214372;
		vBinVals[10][0][3] = -0.214372;
		vBinVals[10][0][4] = -0.214372;
		vBinVals[10][0][5] = -0.214372;
		vBinVals[10][0][6] = -0.214372;
		vBinVals[10][0][7] = -0.214372;
		vBinVals[10][0][8] = -0.214372;
		vBinVals[10][0][9] = -0.214372;
		vBinVals[10][0][10] = -0.214372;
		vBinVals[10][0][11] = -0.214372;
		vBinVals[10][0][12] = -0.214372;
		vBinVals[10][0][13] = -0.214372;
		vBinVals[10][0][14] = -0.214372;
		vBinVals[10][0][15] = -0.214372;
		vBinVals[10][0][16] = -0.214372;
		vBinVals[10][0][17] = -0.130108;
		vBinVals[10][0][18] = -0.130108;
		vBinVals[10][0][19] = -0.130108;
		vBinVals[10][0][20] = -0.130108;
		vBinVals[10][0][21] = -0.130108;
		vBinVals[10][0][22] = -0.0781485;
		vBinVals[10][0][23] = -0.0781485;
		vBinVals[10][0][24] = -0.0781485;
		vBinVals[10][0][25] = -0.0781485;
		vBinVals[10][0][26] = -0.0781485;
		vBinVals[10][0][27] = -0.0196848;
		vBinVals[10][0][28] = -0.0196848;
		vBinVals[10][0][29] = -0.0196848;
		vBinVals[10][0][30] = -0.0196848;
		vBinVals[10][0][31] = -0.0196848;
		vBinVals[10][0][32] = -0.0667353;
		vBinVals[10][0][33] = -0.0667353;
		vBinVals[10][0][34] = -0.0667353;
		vBinVals[10][0][35] = -0.0667353;
		vBinVals[10][0][36] = -0.0667353;
		vBinVals[10][0][37] = -0.0381784;
		vBinVals[10][0][38] = -0.0381784;
		vBinVals[10][0][39] = -0.0381784;
		vBinVals[10][0][40] = -0.0381784;
		vBinVals[10][0][41] = -0.0381784;
		vBinVals[10][0][42] = -0.0100892;
		vBinVals[10][0][43] = -0.0100892;
		vBinVals[10][0][44] = -0.0100892;
		vBinVals[10][0][45] = -0.0100892;
		vBinVals[10][0][46] = -0.0100892;
		vBinVals[10][0][47] = -0.0189216;
		vBinVals[10][0][48] = -0.0189216;
		vBinVals[10][0][49] = -0.0189216;
		vBinVals[10][0][50] = -0.0189216;
		vBinVals[10][0][51] = -0.0189216;
		vBinVals[10][0][52] = -0.0149831;
		vBinVals[10][0][53] = -0.0149831;
		vBinVals[10][0][54] = -0.0149831;
		vBinVals[10][0][55] = -0.0149831;
		vBinVals[10][0][56] = -0.0149831;
		vBinVals[10][0][57] = -0.0372589;
		vBinVals[10][0][58] = -0.0372589;
		vBinVals[10][0][59] = -0.0372589;
		vBinVals[10][1][0] = 0.130527;
		vBinVals[10][1][1] = 0.130527;
		vBinVals[10][1][2] = 0.130527;
		vBinVals[10][1][3] = 0.130527;
		vBinVals[10][1][4] = 0.130527;
		vBinVals[10][1][5] = 0.130527;
		vBinVals[10][1][6] = 0.130527;
		vBinVals[10][1][7] = 0.130527;
		vBinVals[10][1][8] = 0.130527;
		vBinVals[10][1][9] = 0.130527;
		vBinVals[10][1][10] = 0.130527;
		vBinVals[10][1][11] = 0.130527;
		vBinVals[10][1][12] = 0.130527;
		vBinVals[10][1][13] = 0.130527;
		vBinVals[10][1][14] = 0.130527;
		vBinVals[10][1][15] = 0.130527;
		vBinVals[10][1][16] = 0.130527;
		vBinVals[10][1][17] = 0.155795;
		vBinVals[10][1][18] = 0.155795;
		vBinVals[10][1][19] = 0.155795;
		vBinVals[10][1][20] = 0.155795;
		vBinVals[10][1][21] = 0.155795;
		vBinVals[10][1][22] = 0.160258;
		vBinVals[10][1][23] = 0.160258;
		vBinVals[10][1][24] = 0.160258;
		vBinVals[10][1][25] = 0.160258;
		vBinVals[10][1][26] = 0.160258;
		vBinVals[10][1][27] = 0.182785;
		vBinVals[10][1][28] = 0.182785;
		vBinVals[10][1][29] = 0.182785;
		vBinVals[10][1][30] = 0.182785;
		vBinVals[10][1][31] = 0.182785;
		vBinVals[10][1][32] = 0.165356;
		vBinVals[10][1][33] = 0.165356;
		vBinVals[10][1][34] = 0.165356;
		vBinVals[10][1][35] = 0.165356;
		vBinVals[10][1][36] = 0.165356;
		vBinVals[10][1][37] = 0.160206;
		vBinVals[10][1][38] = 0.160206;
		vBinVals[10][1][39] = 0.160206;
		vBinVals[10][1][40] = 0.160206;
		vBinVals[10][1][41] = 0.160206;
		vBinVals[10][1][42] = 0.168024;
		vBinVals[10][1][43] = 0.168024;
		vBinVals[10][1][44] = 0.168024;
		vBinVals[10][1][45] = 0.168024;
		vBinVals[10][1][46] = 0.168024;
		vBinVals[10][1][47] = 0.173575;
		vBinVals[10][1][48] = 0.173575;
		vBinVals[10][1][49] = 0.173575;
		vBinVals[10][1][50] = 0.173575;
		vBinVals[10][1][51] = 0.173575;
		vBinVals[10][1][52] = 0.149928;
		vBinVals[10][1][53] = 0.149928;
		vBinVals[10][1][54] = 0.149928;
		vBinVals[10][1][55] = 0.149928;
		vBinVals[10][1][56] = 0.149928;
		vBinVals[10][1][57] = 0.132898;
		vBinVals[10][1][58] = 0.132898;
		vBinVals[10][1][59] = 0.132898;
		vBinVals[10][2][0] = 9.36426;
		vBinVals[10][2][1] = 9.36426;
		vBinVals[10][2][2] = 9.36426;
		vBinVals[10][2][3] = 9.36426;
		vBinVals[10][2][4] = 9.36426;
		vBinVals[10][2][5] = 9.36426;
		vBinVals[10][2][6] = 9.36426;
		vBinVals[10][2][7] = 9.36426;
		vBinVals[10][2][8] = 9.36426;
		vBinVals[10][2][9] = 9.36426;
		vBinVals[10][2][10] = 9.36426;
		vBinVals[10][2][11] = 9.36426;
		vBinVals[10][2][12] = 9.36426;
		vBinVals[10][2][13] = 9.36426;
		vBinVals[10][2][14] = 9.36426;
		vBinVals[10][2][15] = 9.36426;
		vBinVals[10][2][16] = 9.36426;
		vBinVals[10][2][17] = -0.563002;
		vBinVals[10][2][18] = -0.563002;
		vBinVals[10][2][19] = -0.563002;
		vBinVals[10][2][20] = -0.563002;
		vBinVals[10][2][21] = -0.563002;
		vBinVals[10][2][22] = -3.29963;
		vBinVals[10][2][23] = -3.29963;
		vBinVals[10][2][24] = -3.29963;
		vBinVals[10][2][25] = -3.29963;
		vBinVals[10][2][26] = -3.29963;
		vBinVals[10][2][27] = -0.175652;
		vBinVals[10][2][28] = -0.175652;
		vBinVals[10][2][29] = -0.175652;
		vBinVals[10][2][30] = -0.175652;
		vBinVals[10][2][31] = -0.175652;
		vBinVals[10][2][32] = -0.11556;
		vBinVals[10][2][33] = -0.11556;
		vBinVals[10][2][34] = -0.11556;
		vBinVals[10][2][35] = -0.11556;
		vBinVals[10][2][36] = -0.11556;
		vBinVals[10][2][37] = -0.116908;
		vBinVals[10][2][38] = -0.116908;
		vBinVals[10][2][39] = -0.116908;
		vBinVals[10][2][40] = -0.116908;
		vBinVals[10][2][41] = -0.116908;
		vBinVals[10][2][42] = -0.114027;
		vBinVals[10][2][43] = -0.114027;
		vBinVals[10][2][44] = -0.114027;
		vBinVals[10][2][45] = -0.114027;
		vBinVals[10][2][46] = -0.114027;
		vBinVals[10][2][47] = -0.108643;
		vBinVals[10][2][48] = -0.108643;
		vBinVals[10][2][49] = -0.108643;
		vBinVals[10][2][50] = -0.108643;
		vBinVals[10][2][51] = -0.108643;
		vBinVals[10][2][52] = -0.121751;
		vBinVals[10][2][53] = -0.121751;
		vBinVals[10][2][54] = -0.121751;
		vBinVals[10][2][55] = -0.121751;
		vBinVals[10][2][56] = -0.121751;
		vBinVals[10][2][57] = -0.0918439;
		vBinVals[10][2][58] = -0.0918439;
		vBinVals[10][2][59] = -0.0918439;
		vBinVals[10][3][0] = 3.491;
		vBinVals[10][3][1] = 3.491;
		vBinVals[10][3][2] = 3.491;
		vBinVals[10][3][3] = 3.491;
		vBinVals[10][3][4] = 3.491;
		vBinVals[10][3][5] = 3.491;
		vBinVals[10][3][6] = 3.491;
		vBinVals[10][3][7] = 3.491;
		vBinVals[10][3][8] = 3.491;
		vBinVals[10][3][9] = 3.491;
		vBinVals[10][3][10] = 3.491;
		vBinVals[10][3][11] = 3.491;
		vBinVals[10][3][12] = 3.491;
		vBinVals[10][3][13] = 3.491;
		vBinVals[10][3][14] = 3.491;
		vBinVals[10][3][15] = 3.491;
		vBinVals[10][3][16] = 3.491;
		vBinVals[10][3][17] = -8511.14;
		vBinVals[10][3][18] = -8511.14;
		vBinVals[10][3][19] = -8511.14;
		vBinVals[10][3][20] = -8511.14;
		vBinVals[10][3][21] = -8511.14;
		vBinVals[10][3][22] = -0.390276;
		vBinVals[10][3][23] = -0.390276;
		vBinVals[10][3][24] = -0.390276;
		vBinVals[10][3][25] = -0.390276;
		vBinVals[10][3][26] = -0.390276;
		vBinVals[10][3][27] = 0.0878062;
		vBinVals[10][3][28] = 0.0878062;
		vBinVals[10][3][29] = 0.0878062;
		vBinVals[10][3][30] = 0.0878062;
		vBinVals[10][3][31] = 0.0878062;
		vBinVals[10][3][32] = 0.104709;
		vBinVals[10][3][33] = 0.104709;
		vBinVals[10][3][34] = 0.104709;
		vBinVals[10][3][35] = 0.104709;
		vBinVals[10][3][36] = 0.104709;
		vBinVals[10][3][37] = 0.0929202;
		vBinVals[10][3][38] = 0.0929202;
		vBinVals[10][3][39] = 0.0929202;
		vBinVals[10][3][40] = 0.0929202;
		vBinVals[10][3][41] = 0.0929202;
		vBinVals[10][3][42] = 0.091115;
		vBinVals[10][3][43] = 0.091115;
		vBinVals[10][3][44] = 0.091115;
		vBinVals[10][3][45] = 0.091115;
		vBinVals[10][3][46] = 0.091115;
		vBinVals[10][3][47] = 0.0815337;
		vBinVals[10][3][48] = 0.0815337;
		vBinVals[10][3][49] = 0.0815337;
		vBinVals[10][3][50] = 0.0815337;
		vBinVals[10][3][51] = 0.0815337;
		vBinVals[10][3][52] = 0.0766127;
		vBinVals[10][3][53] = 0.0766127;
		vBinVals[10][3][54] = 0.0766127;
		vBinVals[10][3][55] = 0.0766127;
		vBinVals[10][3][56] = 0.0766127;
		vBinVals[10][3][57] = 0.0831976;
		vBinVals[10][3][58] = 0.0831976;
		vBinVals[10][3][59] = 0.0831976;
		vBinVals[10][4][0] = 0.449133;
		vBinVals[10][4][1] = 0.449133;
		vBinVals[10][4][2] = 0.449133;
		vBinVals[10][4][3] = 0.449133;
		vBinVals[10][4][4] = 0.449133;
		vBinVals[10][4][5] = 0.449133;
		vBinVals[10][4][6] = 0.449133;
		vBinVals[10][4][7] = 0.449133;
		vBinVals[10][4][8] = 0.449133;
		vBinVals[10][4][9] = 0.449133;
		vBinVals[10][4][10] = 0.449133;
		vBinVals[10][4][11] = 0.449133;
		vBinVals[10][4][12] = 0.449133;
		vBinVals[10][4][13] = 0.449133;
		vBinVals[10][4][14] = 0.449133;
		vBinVals[10][4][15] = 0.449133;
		vBinVals[10][4][16] = 0.449133;
		vBinVals[10][4][17] = 1.04377;
		vBinVals[10][4][18] = 1.04377;
		vBinVals[10][4][19] = 1.04377;
		vBinVals[10][4][20] = 1.04377;
		vBinVals[10][4][21] = 1.04377;
		vBinVals[10][4][22] = 1.20667;
		vBinVals[10][4][23] = 1.20667;
		vBinVals[10][4][24] = 1.20667;
		vBinVals[10][4][25] = 1.20667;
		vBinVals[10][4][26] = 1.20667;
		vBinVals[10][4][27] = 0.214821;
		vBinVals[10][4][28] = 0.214821;
		vBinVals[10][4][29] = 0.214821;
		vBinVals[10][4][30] = 0.214821;
		vBinVals[10][4][31] = 0.214821;
		vBinVals[10][4][32] = 0.14552;
		vBinVals[10][4][33] = 0.14552;
		vBinVals[10][4][34] = 0.14552;
		vBinVals[10][4][35] = 0.14552;
		vBinVals[10][4][36] = 0.14552;
		vBinVals[10][4][37] = 0.158247;
		vBinVals[10][4][38] = 0.158247;
		vBinVals[10][4][39] = 0.158247;
		vBinVals[10][4][40] = 0.158247;
		vBinVals[10][4][41] = 0.158247;
		vBinVals[10][4][42] = 0.144003;
		vBinVals[10][4][43] = 0.144003;
		vBinVals[10][4][44] = 0.144003;
		vBinVals[10][4][45] = 0.144003;
		vBinVals[10][4][46] = 0.144003;
		vBinVals[10][4][47] = 0.0864884;
		vBinVals[10][4][48] = 0.0864884;
		vBinVals[10][4][49] = 0.0864884;
		vBinVals[10][4][50] = 0.0864884;
		vBinVals[10][4][51] = 0.0864884;
		vBinVals[10][4][52] = 0.157069;
		vBinVals[10][4][53] = 0.157069;
		vBinVals[10][4][54] = 0.157069;
		vBinVals[10][4][55] = 0.157069;
		vBinVals[10][4][56] = 0.157069;
		vBinVals[10][4][57] = 0.257177;
		vBinVals[10][4][58] = 0.257177;
		vBinVals[10][4][59] = 0.257177;
		vBinVals[11][0][0] = -0.26084;
		vBinVals[11][0][1] = -0.26084;
		vBinVals[11][0][2] = -0.26084;
		vBinVals[11][0][3] = -0.26084;
		vBinVals[11][0][4] = -0.26084;
		vBinVals[11][0][5] = -0.26084;
		vBinVals[11][0][6] = -0.26084;
		vBinVals[11][0][7] = -0.26084;
		vBinVals[11][0][8] = -0.26084;
		vBinVals[11][0][9] = -0.26084;
		vBinVals[11][0][10] = -0.26084;
		vBinVals[11][0][11] = -0.26084;
		vBinVals[11][0][12] = -0.26084;
		vBinVals[11][0][13] = -0.26084;
		vBinVals[11][0][14] = -0.26084;
		vBinVals[11][0][15] = -0.26084;
		vBinVals[11][0][16] = -0.26084;
		vBinVals[11][0][17] = -0.205177;
		vBinVals[11][0][18] = -0.205177;
		vBinVals[11][0][19] = -0.205177;
		vBinVals[11][0][20] = -0.205177;
		vBinVals[11][0][21] = -0.205177;
		vBinVals[11][0][22] = -0.137856;
		vBinVals[11][0][23] = -0.137856;
		vBinVals[11][0][24] = -0.137856;
		vBinVals[11][0][25] = -0.137856;
		vBinVals[11][0][26] = -0.137856;
		vBinVals[11][0][27] = -0.108698;
		vBinVals[11][0][28] = -0.108698;
		vBinVals[11][0][29] = -0.108698;
		vBinVals[11][0][30] = -0.108698;
		vBinVals[11][0][31] = -0.108698;
		vBinVals[11][0][32] = -0.0540388;
		vBinVals[11][0][33] = -0.0540388;
		vBinVals[11][0][34] = -0.0540388;
		vBinVals[11][0][35] = -0.0540388;
		vBinVals[11][0][36] = -0.0540388;
		vBinVals[11][0][37] = -0.0662536;
		vBinVals[11][0][38] = -0.0662536;
		vBinVals[11][0][39] = -0.0662536;
		vBinVals[11][0][40] = -0.0662536;
		vBinVals[11][0][41] = -0.0662536;
		vBinVals[11][0][42] = -0.0640657;
		vBinVals[11][0][43] = -0.0640657;
		vBinVals[11][0][44] = -0.0640657;
		vBinVals[11][0][45] = -0.0640657;
		vBinVals[11][0][46] = -0.0640657;
		vBinVals[11][0][47] = -0.0280409;
		vBinVals[11][0][48] = -0.0280409;
		vBinVals[11][0][49] = -0.0280409;
		vBinVals[11][0][50] = -0.0280409;
		vBinVals[11][0][51] = -0.0280409;
		vBinVals[11][0][52] = -0.0382251;
		vBinVals[11][0][53] = -0.0382251;
		vBinVals[11][0][54] = -0.0382251;
		vBinVals[11][0][55] = -0.0382251;
		vBinVals[11][0][56] = -0.0382251;
		vBinVals[11][0][57] = -0.055633;
		vBinVals[11][0][58] = -0.055633;
		vBinVals[11][0][59] = -0.055633;
		vBinVals[11][1][0] = 0.103497;
		vBinVals[11][1][1] = 0.103497;
		vBinVals[11][1][2] = 0.103497;
		vBinVals[11][1][3] = 0.103497;
		vBinVals[11][1][4] = 0.103497;
		vBinVals[11][1][5] = 0.103497;
		vBinVals[11][1][6] = 0.103497;
		vBinVals[11][1][7] = 0.103497;
		vBinVals[11][1][8] = 0.103497;
		vBinVals[11][1][9] = 0.103497;
		vBinVals[11][1][10] = 0.103497;
		vBinVals[11][1][11] = 0.103497;
		vBinVals[11][1][12] = 0.103497;
		vBinVals[11][1][13] = 0.103497;
		vBinVals[11][1][14] = 0.103497;
		vBinVals[11][1][15] = 0.103497;
		vBinVals[11][1][16] = 0.103497;
		vBinVals[11][1][17] = 0.129093;
		vBinVals[11][1][18] = 0.129093;
		vBinVals[11][1][19] = 0.129093;
		vBinVals[11][1][20] = 0.129093;
		vBinVals[11][1][21] = 0.129093;
		vBinVals[11][1][22] = 0.146381;
		vBinVals[11][1][23] = 0.146381;
		vBinVals[11][1][24] = 0.146381;
		vBinVals[11][1][25] = 0.146381;
		vBinVals[11][1][26] = 0.146381;
		vBinVals[11][1][27] = 0.14494;
		vBinVals[11][1][28] = 0.14494;
		vBinVals[11][1][29] = 0.14494;
		vBinVals[11][1][30] = 0.14494;
		vBinVals[11][1][31] = 0.14494;
		vBinVals[11][1][32] = 0.162809;
		vBinVals[11][1][33] = 0.162809;
		vBinVals[11][1][34] = 0.162809;
		vBinVals[11][1][35] = 0.162809;
		vBinVals[11][1][36] = 0.162809;
		vBinVals[11][1][37] = 0.167117;
		vBinVals[11][1][38] = 0.167117;
		vBinVals[11][1][39] = 0.167117;
		vBinVals[11][1][40] = 0.167117;
		vBinVals[11][1][41] = 0.167117;
		vBinVals[11][1][42] = 0.149112;
		vBinVals[11][1][43] = 0.149112;
		vBinVals[11][1][44] = 0.149112;
		vBinVals[11][1][45] = 0.149112;
		vBinVals[11][1][46] = 0.149112;
		vBinVals[11][1][47] = 0.175347;
		vBinVals[11][1][48] = 0.175347;
		vBinVals[11][1][49] = 0.175347;
		vBinVals[11][1][50] = 0.175347;
		vBinVals[11][1][51] = 0.175347;
		vBinVals[11][1][52] = 0.177025;
		vBinVals[11][1][53] = 0.177025;
		vBinVals[11][1][54] = 0.177025;
		vBinVals[11][1][55] = 0.177025;
		vBinVals[11][1][56] = 0.177025;
		vBinVals[11][1][57] = 0.141439;
		vBinVals[11][1][58] = 0.141439;
		vBinVals[11][1][59] = 0.141439;
		vBinVals[11][2][0] = -0.58675;
		vBinVals[11][2][1] = -0.58675;
		vBinVals[11][2][2] = -0.58675;
		vBinVals[11][2][3] = -0.58675;
		vBinVals[11][2][4] = -0.58675;
		vBinVals[11][2][5] = -0.58675;
		vBinVals[11][2][6] = -0.58675;
		vBinVals[11][2][7] = -0.58675;
		vBinVals[11][2][8] = -0.58675;
		vBinVals[11][2][9] = -0.58675;
		vBinVals[11][2][10] = -0.58675;
		vBinVals[11][2][11] = -0.58675;
		vBinVals[11][2][12] = -0.58675;
		vBinVals[11][2][13] = -0.58675;
		vBinVals[11][2][14] = -0.58675;
		vBinVals[11][2][15] = -0.58675;
		vBinVals[11][2][16] = -0.58675;
		vBinVals[11][2][17] = -0.285943;
		vBinVals[11][2][18] = -0.285943;
		vBinVals[11][2][19] = -0.285943;
		vBinVals[11][2][20] = -0.285943;
		vBinVals[11][2][21] = -0.285943;
		vBinVals[11][2][22] = -0.258856;
		vBinVals[11][2][23] = -0.258856;
		vBinVals[11][2][24] = -0.258856;
		vBinVals[11][2][25] = -0.258856;
		vBinVals[11][2][26] = -0.258856;
		vBinVals[11][2][27] = 0.224876;
		vBinVals[11][2][28] = 0.224876;
		vBinVals[11][2][29] = 0.224876;
		vBinVals[11][2][30] = 0.224876;
		vBinVals[11][2][31] = 0.224876;
		vBinVals[11][2][32] = -0.198798;
		vBinVals[11][2][33] = -0.198798;
		vBinVals[11][2][34] = -0.198798;
		vBinVals[11][2][35] = -0.198798;
		vBinVals[11][2][36] = -0.198798;
		vBinVals[11][2][37] = -0.154446;
		vBinVals[11][2][38] = -0.154446;
		vBinVals[11][2][39] = -0.154446;
		vBinVals[11][2][40] = -0.154446;
		vBinVals[11][2][41] = -0.154446;
		vBinVals[11][2][42] = -429812;
		vBinVals[11][2][43] = -429812;
		vBinVals[11][2][44] = -429812;
		vBinVals[11][2][45] = -429812;
		vBinVals[11][2][46] = -429812;
		vBinVals[11][2][47] = -0.140152;
		vBinVals[11][2][48] = -0.140152;
		vBinVals[11][2][49] = -0.140152;
		vBinVals[11][2][50] = -0.140152;
		vBinVals[11][2][51] = -0.140152;
		vBinVals[11][2][52] = -0.11954;
		vBinVals[11][2][53] = -0.11954;
		vBinVals[11][2][54] = -0.11954;
		vBinVals[11][2][55] = -0.11954;
		vBinVals[11][2][56] = -0.11954;
		vBinVals[11][2][57] = -0.108259;
		vBinVals[11][2][58] = -0.108259;
		vBinVals[11][2][59] = -0.108259;
		vBinVals[11][3][0] = 0.0506118;
		vBinVals[11][3][1] = 0.0506118;
		vBinVals[11][3][2] = 0.0506118;
		vBinVals[11][3][3] = 0.0506118;
		vBinVals[11][3][4] = 0.0506118;
		vBinVals[11][3][5] = 0.0506118;
		vBinVals[11][3][6] = 0.0506118;
		vBinVals[11][3][7] = 0.0506118;
		vBinVals[11][3][8] = 0.0506118;
		vBinVals[11][3][9] = 0.0506118;
		vBinVals[11][3][10] = 0.0506118;
		vBinVals[11][3][11] = 0.0506118;
		vBinVals[11][3][12] = 0.0506118;
		vBinVals[11][3][13] = 0.0506118;
		vBinVals[11][3][14] = 0.0506118;
		vBinVals[11][3][15] = 0.0506118;
		vBinVals[11][3][16] = 0.0506118;
		vBinVals[11][3][17] = 0.0834638;
		vBinVals[11][3][18] = 0.0834638;
		vBinVals[11][3][19] = 0.0834638;
		vBinVals[11][3][20] = 0.0834638;
		vBinVals[11][3][21] = 0.0834638;
		vBinVals[11][3][22] = 0.0673625;
		vBinVals[11][3][23] = 0.0673625;
		vBinVals[11][3][24] = 0.0673625;
		vBinVals[11][3][25] = 0.0673625;
		vBinVals[11][3][26] = 0.0673625;
		vBinVals[11][3][27] = -5.03395;
		vBinVals[11][3][28] = -5.03395;
		vBinVals[11][3][29] = -5.03395;
		vBinVals[11][3][30] = -5.03395;
		vBinVals[11][3][31] = -5.03395;
		vBinVals[11][3][32] = 0.0761566;
		vBinVals[11][3][33] = 0.0761566;
		vBinVals[11][3][34] = 0.0761566;
		vBinVals[11][3][35] = 0.0761566;
		vBinVals[11][3][36] = 0.0761566;
		vBinVals[11][3][37] = 0.0879282;
		vBinVals[11][3][38] = 0.0879282;
		vBinVals[11][3][39] = 0.0879282;
		vBinVals[11][3][40] = 0.0879282;
		vBinVals[11][3][41] = 0.0879282;
		vBinVals[11][3][42] = 13919.3;
		vBinVals[11][3][43] = 13919.3;
		vBinVals[11][3][44] = 13919.3;
		vBinVals[11][3][45] = 13919.3;
		vBinVals[11][3][46] = 13919.3;
		vBinVals[11][3][47] = 0.0768985;
		vBinVals[11][3][48] = 0.0768985;
		vBinVals[11][3][49] = 0.0768985;
		vBinVals[11][3][50] = 0.0768985;
		vBinVals[11][3][51] = 0.0768985;
		vBinVals[11][3][52] = 0.0842906;
		vBinVals[11][3][53] = 0.0842906;
		vBinVals[11][3][54] = 0.0842906;
		vBinVals[11][3][55] = 0.0842906;
		vBinVals[11][3][56] = 0.0842906;
		vBinVals[11][3][57] = 0.0897124;
		vBinVals[11][3][58] = 0.0897124;
		vBinVals[11][3][59] = 0.0897124;
		vBinVals[11][4][0] = 7.45845;
		vBinVals[11][4][1] = 7.45845;
		vBinVals[11][4][2] = 7.45845;
		vBinVals[11][4][3] = 7.45845;
		vBinVals[11][4][4] = 7.45845;
		vBinVals[11][4][5] = 7.45845;
		vBinVals[11][4][6] = 7.45845;
		vBinVals[11][4][7] = 7.45845;
		vBinVals[11][4][8] = 7.45845;
		vBinVals[11][4][9] = 7.45845;
		vBinVals[11][4][10] = 7.45845;
		vBinVals[11][4][11] = 7.45845;
		vBinVals[11][4][12] = 7.45845;
		vBinVals[11][4][13] = 7.45845;
		vBinVals[11][4][14] = 7.45845;
		vBinVals[11][4][15] = 7.45845;
		vBinVals[11][4][16] = 7.45845;
		vBinVals[11][4][17] = 1.18386;
		vBinVals[11][4][18] = 1.18386;
		vBinVals[11][4][19] = 1.18386;
		vBinVals[11][4][20] = 1.18386;
		vBinVals[11][4][21] = 1.18386;
		vBinVals[11][4][22] = 0.583168;
		vBinVals[11][4][23] = 0.583168;
		vBinVals[11][4][24] = 0.583168;
		vBinVals[11][4][25] = 0.583168;
		vBinVals[11][4][26] = 0.583168;
		vBinVals[11][4][27] = 0.66301;
		vBinVals[11][4][28] = 0.66301;
		vBinVals[11][4][29] = 0.66301;
		vBinVals[11][4][30] = 0.66301;
		vBinVals[11][4][31] = 0.66301;
		vBinVals[11][4][32] = 0.190734;
		vBinVals[11][4][33] = 0.190734;
		vBinVals[11][4][34] = 0.190734;
		vBinVals[11][4][35] = 0.190734;
		vBinVals[11][4][36] = 0.190734;
		vBinVals[11][4][37] = 0.0982949;
		vBinVals[11][4][38] = 0.0982949;
		vBinVals[11][4][39] = 0.0982949;
		vBinVals[11][4][40] = 0.0982949;
		vBinVals[11][4][41] = 0.0982949;
		vBinVals[11][4][42] = 0.169488;
		vBinVals[11][4][43] = 0.169488;
		vBinVals[11][4][44] = 0.169488;
		vBinVals[11][4][45] = 0.169488;
		vBinVals[11][4][46] = 0.169488;
		vBinVals[11][4][47] = 0.124994;
		vBinVals[11][4][48] = 0.124994;
		vBinVals[11][4][49] = 0.124994;
		vBinVals[11][4][50] = 0.124994;
		vBinVals[11][4][51] = 0.124994;
		vBinVals[11][4][52] = 0.0868426;
		vBinVals[11][4][53] = 0.0868426;
		vBinVals[11][4][54] = 0.0868426;
		vBinVals[11][4][55] = 0.0868426;
		vBinVals[11][4][56] = 0.0868426;
		vBinVals[11][4][57] = 0.264306;
		vBinVals[11][4][58] = 0.264306;
		vBinVals[11][4][59] = 0.264306;
		vBinVals[12][0][0] = -0.255203;
		vBinVals[12][0][1] = -0.255203;
		vBinVals[12][0][2] = -0.255203;
		vBinVals[12][0][3] = -0.255203;
		vBinVals[12][0][4] = -0.255203;
		vBinVals[12][0][5] = -0.255203;
		vBinVals[12][0][6] = -0.255203;
		vBinVals[12][0][7] = -0.255203;
		vBinVals[12][0][8] = -0.255203;
		vBinVals[12][0][9] = -0.255203;
		vBinVals[12][0][10] = -0.255203;
		vBinVals[12][0][11] = -0.255203;
		vBinVals[12][0][12] = -0.255203;
		vBinVals[12][0][13] = -0.255203;
		vBinVals[12][0][14] = -0.255203;
		vBinVals[12][0][15] = -0.255203;
		vBinVals[12][0][16] = -0.255203;
		vBinVals[12][0][17] = -0.255203;
		vBinVals[12][0][18] = -0.255203;
		vBinVals[12][0][19] = -0.255203;
		vBinVals[12][0][20] = -0.255203;
		vBinVals[12][0][21] = -0.255203;
		vBinVals[12][0][22] = -0.218723;
		vBinVals[12][0][23] = -0.218723;
		vBinVals[12][0][24] = -0.218723;
		vBinVals[12][0][25] = -0.218723;
		vBinVals[12][0][26] = -0.218723;
		vBinVals[12][0][27] = -0.177037;
		vBinVals[12][0][28] = -0.177037;
		vBinVals[12][0][29] = -0.177037;
		vBinVals[12][0][30] = -0.177037;
		vBinVals[12][0][31] = -0.177037;
		vBinVals[12][0][32] = -0.161735;
		vBinVals[12][0][33] = -0.161735;
		vBinVals[12][0][34] = -0.161735;
		vBinVals[12][0][35] = -0.161735;
		vBinVals[12][0][36] = -0.161735;
		vBinVals[12][0][37] = -0.139786;
		vBinVals[12][0][38] = -0.139786;
		vBinVals[12][0][39] = -0.139786;
		vBinVals[12][0][40] = -0.139786;
		vBinVals[12][0][41] = -0.139786;
		vBinVals[12][0][42] = -0.131019;
		vBinVals[12][0][43] = -0.131019;
		vBinVals[12][0][44] = -0.131019;
		vBinVals[12][0][45] = -0.131019;
		vBinVals[12][0][46] = -0.131019;
		vBinVals[12][0][47] = -0.105494;
		vBinVals[12][0][48] = -0.105494;
		vBinVals[12][0][49] = -0.105494;
		vBinVals[12][0][50] = -0.105494;
		vBinVals[12][0][51] = -0.105494;
		vBinVals[12][0][52] = -0.140275;
		vBinVals[12][0][53] = -0.140275;
		vBinVals[12][0][54] = -0.140275;
		vBinVals[12][0][55] = -0.140275;
		vBinVals[12][0][56] = -0.140275;
		vBinVals[12][0][57] = -0.145035;
		vBinVals[12][0][58] = -0.145035;
		vBinVals[12][0][59] = -0.145035;
		vBinVals[12][1][0] = 0.11461;
		vBinVals[12][1][1] = 0.11461;
		vBinVals[12][1][2] = 0.11461;
		vBinVals[12][1][3] = 0.11461;
		vBinVals[12][1][4] = 0.11461;
		vBinVals[12][1][5] = 0.11461;
		vBinVals[12][1][6] = 0.11461;
		vBinVals[12][1][7] = 0.11461;
		vBinVals[12][1][8] = 0.11461;
		vBinVals[12][1][9] = 0.11461;
		vBinVals[12][1][10] = 0.11461;
		vBinVals[12][1][11] = 0.11461;
		vBinVals[12][1][12] = 0.11461;
		vBinVals[12][1][13] = 0.11461;
		vBinVals[12][1][14] = 0.11461;
		vBinVals[12][1][15] = 0.11461;
		vBinVals[12][1][16] = 0.11461;
		vBinVals[12][1][17] = 0.11461;
		vBinVals[12][1][18] = 0.11461;
		vBinVals[12][1][19] = 0.11461;
		vBinVals[12][1][20] = 0.11461;
		vBinVals[12][1][21] = 0.11461;
		vBinVals[12][1][22] = 0.127392;
		vBinVals[12][1][23] = 0.127392;
		vBinVals[12][1][24] = 0.127392;
		vBinVals[12][1][25] = 0.127392;
		vBinVals[12][1][26] = 0.127392;
		vBinVals[12][1][27] = 0.136441;
		vBinVals[12][1][28] = 0.136441;
		vBinVals[12][1][29] = 0.136441;
		vBinVals[12][1][30] = 0.136441;
		vBinVals[12][1][31] = 0.136441;
		vBinVals[12][1][32] = -0.134707;
		vBinVals[12][1][33] = -0.134707;
		vBinVals[12][1][34] = -0.134707;
		vBinVals[12][1][35] = -0.134707;
		vBinVals[12][1][36] = -0.134707;
		vBinVals[12][1][37] = 0.145651;
		vBinVals[12][1][38] = 0.145651;
		vBinVals[12][1][39] = 0.145651;
		vBinVals[12][1][40] = 0.145651;
		vBinVals[12][1][41] = 0.145651;
		vBinVals[12][1][42] = 0.14431;
		vBinVals[12][1][43] = 0.14431;
		vBinVals[12][1][44] = 0.14431;
		vBinVals[12][1][45] = 0.14431;
		vBinVals[12][1][46] = 0.14431;
		vBinVals[12][1][47] = 0.146702;
		vBinVals[12][1][48] = 0.146702;
		vBinVals[12][1][49] = 0.146702;
		vBinVals[12][1][50] = 0.146702;
		vBinVals[12][1][51] = 0.146702;
		vBinVals[12][1][52] = 0.14198;
		vBinVals[12][1][53] = 0.14198;
		vBinVals[12][1][54] = 0.14198;
		vBinVals[12][1][55] = 0.14198;
		vBinVals[12][1][56] = 0.14198;
		vBinVals[12][1][57] = 0.114556;
		vBinVals[12][1][58] = 0.114556;
		vBinVals[12][1][59] = 0.114556;
		vBinVals[12][2][0] = 1721.2;
		vBinVals[12][2][1] = 1721.2;
		vBinVals[12][2][2] = 1721.2;
		vBinVals[12][2][3] = 1721.2;
		vBinVals[12][2][4] = 1721.2;
		vBinVals[12][2][5] = 1721.2;
		vBinVals[12][2][6] = 1721.2;
		vBinVals[12][2][7] = 1721.2;
		vBinVals[12][2][8] = 1721.2;
		vBinVals[12][2][9] = 1721.2;
		vBinVals[12][2][10] = 1721.2;
		vBinVals[12][2][11] = 1721.2;
		vBinVals[12][2][12] = 1721.2;
		vBinVals[12][2][13] = 1721.2;
		vBinVals[12][2][14] = 1721.2;
		vBinVals[12][2][15] = 1721.2;
		vBinVals[12][2][16] = 1721.2;
		vBinVals[12][2][17] = 1721.2;
		vBinVals[12][2][18] = 1721.2;
		vBinVals[12][2][19] = 1721.2;
		vBinVals[12][2][20] = 1721.2;
		vBinVals[12][2][21] = 1721.2;
		vBinVals[12][2][22] = 0.748406;
		vBinVals[12][2][23] = 0.748406;
		vBinVals[12][2][24] = 0.748406;
		vBinVals[12][2][25] = 0.748406;
		vBinVals[12][2][26] = 0.748406;
		vBinVals[12][2][27] = -0.260764;
		vBinVals[12][2][28] = -0.260764;
		vBinVals[12][2][29] = -0.260764;
		vBinVals[12][2][30] = -0.260764;
		vBinVals[12][2][31] = -0.260764;
		vBinVals[12][2][32] = 0.29789;
		vBinVals[12][2][33] = 0.29789;
		vBinVals[12][2][34] = 0.29789;
		vBinVals[12][2][35] = 0.29789;
		vBinVals[12][2][36] = 0.29789;
		vBinVals[12][2][37] = -0.216649;
		vBinVals[12][2][38] = -0.216649;
		vBinVals[12][2][39] = -0.216649;
		vBinVals[12][2][40] = -0.216649;
		vBinVals[12][2][41] = -0.216649;
		vBinVals[12][2][42] = -0.20323;
		vBinVals[12][2][43] = -0.20323;
		vBinVals[12][2][44] = -0.20323;
		vBinVals[12][2][45] = -0.20323;
		vBinVals[12][2][46] = -0.20323;
		vBinVals[12][2][47] = -0.22298;
		vBinVals[12][2][48] = -0.22298;
		vBinVals[12][2][49] = -0.22298;
		vBinVals[12][2][50] = -0.22298;
		vBinVals[12][2][51] = -0.22298;
		vBinVals[12][2][52] = -0.176723;
		vBinVals[12][2][53] = -0.176723;
		vBinVals[12][2][54] = -0.176723;
		vBinVals[12][2][55] = -0.176723;
		vBinVals[12][2][56] = -0.176723;
		vBinVals[12][2][57] = -0.181457;
		vBinVals[12][2][58] = -0.181457;
		vBinVals[12][2][59] = -0.181457;
		vBinVals[12][3][0] = 625.533;
		vBinVals[12][3][1] = 625.533;
		vBinVals[12][3][2] = 625.533;
		vBinVals[12][3][3] = 625.533;
		vBinVals[12][3][4] = 625.533;
		vBinVals[12][3][5] = 625.533;
		vBinVals[12][3][6] = 625.533;
		vBinVals[12][3][7] = 625.533;
		vBinVals[12][3][8] = 625.533;
		vBinVals[12][3][9] = 625.533;
		vBinVals[12][3][10] = 625.533;
		vBinVals[12][3][11] = 625.533;
		vBinVals[12][3][12] = 625.533;
		vBinVals[12][3][13] = 625.533;
		vBinVals[12][3][14] = 625.533;
		vBinVals[12][3][15] = 625.533;
		vBinVals[12][3][16] = 625.533;
		vBinVals[12][3][17] = 625.533;
		vBinVals[12][3][18] = 625.533;
		vBinVals[12][3][19] = 625.533;
		vBinVals[12][3][20] = 625.533;
		vBinVals[12][3][21] = 625.533;
		vBinVals[12][3][22] = -32.7443;
		vBinVals[12][3][23] = -32.7443;
		vBinVals[12][3][24] = -32.7443;
		vBinVals[12][3][25] = -32.7443;
		vBinVals[12][3][26] = -32.7443;
		vBinVals[12][3][27] = 0.0815744;
		vBinVals[12][3][28] = 0.0815744;
		vBinVals[12][3][29] = 0.0815744;
		vBinVals[12][3][30] = 0.0815744;
		vBinVals[12][3][31] = 0.0815744;
		vBinVals[12][3][32] = -0.365661;
		vBinVals[12][3][33] = -0.365661;
		vBinVals[12][3][34] = -0.365661;
		vBinVals[12][3][35] = -0.365661;
		vBinVals[12][3][36] = -0.365661;
		vBinVals[12][3][37] = 0.0816271;
		vBinVals[12][3][38] = 0.0816271;
		vBinVals[12][3][39] = 0.0816271;
		vBinVals[12][3][40] = 0.0816271;
		vBinVals[12][3][41] = 0.0816271;
		vBinVals[12][3][42] = 0.0812139;
		vBinVals[12][3][43] = 0.0812139;
		vBinVals[12][3][44] = 0.0812139;
		vBinVals[12][3][45] = 0.0812139;
		vBinVals[12][3][46] = 0.0812139;
		vBinVals[12][3][47] = 0.0675373;
		vBinVals[12][3][48] = 0.0675373;
		vBinVals[12][3][49] = 0.0675373;
		vBinVals[12][3][50] = 0.0675373;
		vBinVals[12][3][51] = 0.0675373;
		vBinVals[12][3][52] = 0.0751231;
		vBinVals[12][3][53] = 0.0751231;
		vBinVals[12][3][54] = 0.0751231;
		vBinVals[12][3][55] = 0.0751231;
		vBinVals[12][3][56] = 0.0751231;
		vBinVals[12][3][57] = 0.0759548;
		vBinVals[12][3][58] = 0.0759548;
		vBinVals[12][3][59] = 0.0759548;
		vBinVals[12][4][0] = 1.58398;
		vBinVals[12][4][1] = 1.58398;
		vBinVals[12][4][2] = 1.58398;
		vBinVals[12][4][3] = 1.58398;
		vBinVals[12][4][4] = 1.58398;
		vBinVals[12][4][5] = 1.58398;
		vBinVals[12][4][6] = 1.58398;
		vBinVals[12][4][7] = 1.58398;
		vBinVals[12][4][8] = 1.58398;
		vBinVals[12][4][9] = 1.58398;
		vBinVals[12][4][10] = 1.58398;
		vBinVals[12][4][11] = 1.58398;
		vBinVals[12][4][12] = 1.58398;
		vBinVals[12][4][13] = 1.58398;
		vBinVals[12][4][14] = 1.58398;
		vBinVals[12][4][15] = 1.58398;
		vBinVals[12][4][16] = 1.58398;
		vBinVals[12][4][17] = 1.58398;
		vBinVals[12][4][18] = 1.58398;
		vBinVals[12][4][19] = 1.58398;
		vBinVals[12][4][20] = 1.58398;
		vBinVals[12][4][21] = 1.58398;
		vBinVals[12][4][22] = 4.66259;
		vBinVals[12][4][23] = 4.66259;
		vBinVals[12][4][24] = 4.66259;
		vBinVals[12][4][25] = 4.66259;
		vBinVals[12][4][26] = 4.66259;
		vBinVals[12][4][27] = 0.397344;
		vBinVals[12][4][28] = 0.397344;
		vBinVals[12][4][29] = 0.397344;
		vBinVals[12][4][30] = 0.397344;
		vBinVals[12][4][31] = 0.397344;
		vBinVals[12][4][32] = 0.000198344;
		vBinVals[12][4][33] = 0.000198344;
		vBinVals[12][4][34] = 0.000198344;
		vBinVals[12][4][35] = 0.000198344;
		vBinVals[12][4][36] = 0.000198344;
		vBinVals[12][4][37] = 0.168074;
		vBinVals[12][4][38] = 0.168074;
		vBinVals[12][4][39] = 0.168074;
		vBinVals[12][4][40] = 0.168074;
		vBinVals[12][4][41] = 0.168074;
		vBinVals[12][4][42] = 0.162918;
		vBinVals[12][4][43] = 0.162918;
		vBinVals[12][4][44] = 0.162918;
		vBinVals[12][4][45] = 0.162918;
		vBinVals[12][4][46] = 0.162918;
		vBinVals[12][4][47] = 0.180791;
		vBinVals[12][4][48] = 0.180791;
		vBinVals[12][4][49] = 0.180791;
		vBinVals[12][4][50] = 0.180791;
		vBinVals[12][4][51] = 0.180791;
		vBinVals[12][4][52] = 0.0846316;
		vBinVals[12][4][53] = 0.0846316;
		vBinVals[12][4][54] = 0.0846316;
		vBinVals[12][4][55] = 0.0846316;
		vBinVals[12][4][56] = 0.0846316;
		vBinVals[12][4][57] = 0.128283;
		vBinVals[12][4][58] = 0.128283;
		vBinVals[12][4][59] = 0.128283;
		vBinVals[13][0][0] = -0.295968;
		vBinVals[13][0][1] = -0.295968;
		vBinVals[13][0][2] = -0.295968;
		vBinVals[13][0][3] = -0.295968;
		vBinVals[13][0][4] = -0.295968;
		vBinVals[13][0][5] = -0.295968;
		vBinVals[13][0][6] = -0.295968;
		vBinVals[13][0][7] = -0.295968;
		vBinVals[13][0][8] = -0.295968;
		vBinVals[13][0][9] = -0.295968;
		vBinVals[13][0][10] = -0.295968;
		vBinVals[13][0][11] = -0.295968;
		vBinVals[13][0][12] = -0.295968;
		vBinVals[13][0][13] = -0.295968;
		vBinVals[13][0][14] = -0.295968;
		vBinVals[13][0][15] = -0.295968;
		vBinVals[13][0][16] = -0.295968;
		vBinVals[13][0][17] = -0.295968;
		vBinVals[13][0][18] = -0.295968;
		vBinVals[13][0][19] = -0.295968;
		vBinVals[13][0][20] = -0.295968;
		vBinVals[13][0][21] = -0.295968;
		vBinVals[13][0][22] = -0.295968;
		vBinVals[13][0][23] = -0.2643;
		vBinVals[13][0][24] = -0.2643;
		vBinVals[13][0][25] = -0.2643;
		vBinVals[13][0][26] = -0.2643;
		vBinVals[13][0][27] = -0.244189;
		vBinVals[13][0][28] = -0.244189;
		vBinVals[13][0][29] = -0.244189;
		vBinVals[13][0][30] = -0.244189;
		vBinVals[13][0][31] = -0.244189;
		vBinVals[13][0][32] = -0.169571;
		vBinVals[13][0][33] = -0.169571;
		vBinVals[13][0][34] = -0.169571;
		vBinVals[13][0][35] = -0.169571;
		vBinVals[13][0][36] = -0.169571;
		vBinVals[13][0][37] = -0.126272;
		vBinVals[13][0][38] = -0.126272;
		vBinVals[13][0][39] = -0.126272;
		vBinVals[13][0][40] = -0.126272;
		vBinVals[13][0][41] = -0.126272;
		vBinVals[13][0][42] = -0.101372;
		vBinVals[13][0][43] = -0.101372;
		vBinVals[13][0][44] = -0.101372;
		vBinVals[13][0][45] = -0.101372;
		vBinVals[13][0][46] = -0.101372;
		vBinVals[13][0][47] = -0.0894905;
		vBinVals[13][0][48] = -0.0894905;
		vBinVals[13][0][49] = -0.0894905;
		vBinVals[13][0][50] = -0.0894905;
		vBinVals[13][0][51] = -0.0894905;
		vBinVals[13][0][52] = -0.195719;
		vBinVals[13][0][53] = -0.195719;
		vBinVals[13][0][54] = -0.195719;
		vBinVals[13][0][55] = -0.195719;
		vBinVals[13][0][56] = -0.195719;
		vBinVals[13][0][57] = -0.0783068;
		vBinVals[13][0][58] = -0.0783068;
		vBinVals[13][0][59] = -0.0783068;
		vBinVals[13][1][0] = 0.118665;
		vBinVals[13][1][1] = 0.118665;
		vBinVals[13][1][2] = 0.118665;
		vBinVals[13][1][3] = 0.118665;
		vBinVals[13][1][4] = 0.118665;
		vBinVals[13][1][5] = 0.118665;
		vBinVals[13][1][6] = 0.118665;
		vBinVals[13][1][7] = 0.118665;
		vBinVals[13][1][8] = 0.118665;
		vBinVals[13][1][9] = 0.118665;
		vBinVals[13][1][10] = 0.118665;
		vBinVals[13][1][11] = 0.118665;
		vBinVals[13][1][12] = 0.118665;
		vBinVals[13][1][13] = 0.118665;
		vBinVals[13][1][14] = 0.118665;
		vBinVals[13][1][15] = 0.118665;
		vBinVals[13][1][16] = 0.118665;
		vBinVals[13][1][17] = 0.118665;
		vBinVals[13][1][18] = 0.118665;
		vBinVals[13][1][19] = 0.118665;
		vBinVals[13][1][20] = 0.118665;
		vBinVals[13][1][21] = 0.118665;
		vBinVals[13][1][22] = 0.118665;
		vBinVals[13][1][23] = 0.106856;
		vBinVals[13][1][24] = 0.106856;
		vBinVals[13][1][25] = 0.106856;
		vBinVals[13][1][26] = 0.106856;
		vBinVals[13][1][27] = 0.107202;
		vBinVals[13][1][28] = 0.107202;
		vBinVals[13][1][29] = 0.107202;
		vBinVals[13][1][30] = 0.107202;
		vBinVals[13][1][31] = 0.107202;
		vBinVals[13][1][32] = 0.133565;
		vBinVals[13][1][33] = 0.133565;
		vBinVals[13][1][34] = 0.133565;
		vBinVals[13][1][35] = 0.133565;
		vBinVals[13][1][36] = 0.133565;
		vBinVals[13][1][37] = -0.146474;
		vBinVals[13][1][38] = -0.146474;
		vBinVals[13][1][39] = -0.146474;
		vBinVals[13][1][40] = -0.146474;
		vBinVals[13][1][41] = -0.146474;
		vBinVals[13][1][42] = 0.151377;
		vBinVals[13][1][43] = 0.151377;
		vBinVals[13][1][44] = 0.151377;
		vBinVals[13][1][45] = 0.151377;
		vBinVals[13][1][46] = 0.151377;
		vBinVals[13][1][47] = 0.149468;
		vBinVals[13][1][48] = 0.149468;
		vBinVals[13][1][49] = 0.149468;
		vBinVals[13][1][50] = 0.149468;
		vBinVals[13][1][51] = 0.149468;
		vBinVals[13][1][52] = 0.0919117;
		vBinVals[13][1][53] = 0.0919117;
		vBinVals[13][1][54] = 0.0919117;
		vBinVals[13][1][55] = 0.0919117;
		vBinVals[13][1][56] = 0.0919117;
		vBinVals[13][1][57] = 0.148887;
		vBinVals[13][1][58] = 0.148887;
		vBinVals[13][1][59] = 0.148887;
		vBinVals[13][2][0] = -0.340376;
		vBinVals[13][2][1] = -0.340376;
		vBinVals[13][2][2] = -0.340376;
		vBinVals[13][2][3] = -0.340376;
		vBinVals[13][2][4] = -0.340376;
		vBinVals[13][2][5] = -0.340376;
		vBinVals[13][2][6] = -0.340376;
		vBinVals[13][2][7] = -0.340376;
		vBinVals[13][2][8] = -0.340376;
		vBinVals[13][2][9] = -0.340376;
		vBinVals[13][2][10] = -0.340376;
		vBinVals[13][2][11] = -0.340376;
		vBinVals[13][2][12] = -0.340376;
		vBinVals[13][2][13] = -0.340376;
		vBinVals[13][2][14] = -0.340376;
		vBinVals[13][2][15] = -0.340376;
		vBinVals[13][2][16] = -0.340376;
		vBinVals[13][2][17] = -0.340376;
		vBinVals[13][2][18] = -0.340376;
		vBinVals[13][2][19] = -0.340376;
		vBinVals[13][2][20] = -0.340376;
		vBinVals[13][2][21] = -0.340376;
		vBinVals[13][2][22] = -0.340376;
		vBinVals[13][2][23] = 12.4233;
		vBinVals[13][2][24] = 12.4233;
		vBinVals[13][2][25] = 12.4233;
		vBinVals[13][2][26] = 12.4233;
		vBinVals[13][2][27] = -0.0691615;
		vBinVals[13][2][28] = -0.0691615;
		vBinVals[13][2][29] = -0.0691615;
		vBinVals[13][2][30] = -0.0691615;
		vBinVals[13][2][31] = -0.0691615;
		vBinVals[13][2][32] = -0.295977;
		vBinVals[13][2][33] = -0.295977;
		vBinVals[13][2][34] = -0.295977;
		vBinVals[13][2][35] = -0.295977;
		vBinVals[13][2][36] = -0.295977;
		vBinVals[13][2][37] = -0.277854;
		vBinVals[13][2][38] = -0.277854;
		vBinVals[13][2][39] = -0.277854;
		vBinVals[13][2][40] = -0.277854;
		vBinVals[13][2][41] = -0.277854;
		vBinVals[13][2][42] = -0.239583;
		vBinVals[13][2][43] = -0.239583;
		vBinVals[13][2][44] = -0.239583;
		vBinVals[13][2][45] = -0.239583;
		vBinVals[13][2][46] = -0.239583;
		vBinVals[13][2][47] = -0.231759;
		vBinVals[13][2][48] = -0.231759;
		vBinVals[13][2][49] = -0.231759;
		vBinVals[13][2][50] = -0.231759;
		vBinVals[13][2][51] = -0.231759;
		vBinVals[13][2][52] = 0.00276873;
		vBinVals[13][2][53] = 0.00276873;
		vBinVals[13][2][54] = 0.00276873;
		vBinVals[13][2][55] = 0.00276873;
		vBinVals[13][2][56] = 0.00276873;
		vBinVals[13][2][57] = -0.180588;
		vBinVals[13][2][58] = -0.180588;
		vBinVals[13][2][59] = -0.180588;
		vBinVals[13][3][0] = 0.035741;
		vBinVals[13][3][1] = 0.035741;
		vBinVals[13][3][2] = 0.035741;
		vBinVals[13][3][3] = 0.035741;
		vBinVals[13][3][4] = 0.035741;
		vBinVals[13][3][5] = 0.035741;
		vBinVals[13][3][6] = 0.035741;
		vBinVals[13][3][7] = 0.035741;
		vBinVals[13][3][8] = 0.035741;
		vBinVals[13][3][9] = 0.035741;
		vBinVals[13][3][10] = 0.035741;
		vBinVals[13][3][11] = 0.035741;
		vBinVals[13][3][12] = 0.035741;
		vBinVals[13][3][13] = 0.035741;
		vBinVals[13][3][14] = 0.035741;
		vBinVals[13][3][15] = 0.035741;
		vBinVals[13][3][16] = 0.035741;
		vBinVals[13][3][17] = 0.035741;
		vBinVals[13][3][18] = 0.035741;
		vBinVals[13][3][19] = 0.035741;
		vBinVals[13][3][20] = 0.035741;
		vBinVals[13][3][21] = 0.035741;
		vBinVals[13][3][22] = 0.035741;
		vBinVals[13][3][23] = 3.97298;
		vBinVals[13][3][24] = 3.97298;
		vBinVals[13][3][25] = 3.97298;
		vBinVals[13][3][26] = 3.97298;
		vBinVals[13][3][27] = 0.0318418;
		vBinVals[13][3][28] = 0.0318418;
		vBinVals[13][3][29] = 0.0318418;
		vBinVals[13][3][30] = 0.0318418;
		vBinVals[13][3][31] = 0.0318418;
		vBinVals[13][3][32] = 0.0600584;
		vBinVals[13][3][33] = 0.0600584;
		vBinVals[13][3][34] = 0.0600584;
		vBinVals[13][3][35] = 0.0600584;
		vBinVals[13][3][36] = 0.0600584;
		vBinVals[13][3][37] = 0.0608645;
		vBinVals[13][3][38] = 0.0608645;
		vBinVals[13][3][39] = 0.0608645;
		vBinVals[13][3][40] = 0.0608645;
		vBinVals[13][3][41] = 0.0608645;
		vBinVals[13][3][42] = 0.0738441;
		vBinVals[13][3][43] = 0.0738441;
		vBinVals[13][3][44] = 0.0738441;
		vBinVals[13][3][45] = 0.0738441;
		vBinVals[13][3][46] = 0.0738441;
		vBinVals[13][3][47] = 0.0723748;
		vBinVals[13][3][48] = 0.0723748;
		vBinVals[13][3][49] = 0.0723748;
		vBinVals[13][3][50] = 0.0723748;
		vBinVals[13][3][51] = 0.0723748;
		vBinVals[13][3][52] = 0.0647442;
		vBinVals[13][3][53] = 0.0647442;
		vBinVals[13][3][54] = 0.0647442;
		vBinVals[13][3][55] = 0.0647442;
		vBinVals[13][3][56] = 0.0647442;
		vBinVals[13][3][57] = 0.0846228;
		vBinVals[13][3][58] = 0.0846228;
		vBinVals[13][3][59] = 0.0846228;
		vBinVals[13][4][0] = 0.2072;
		vBinVals[13][4][1] = 0.2072;
		vBinVals[13][4][2] = 0.2072;
		vBinVals[13][4][3] = 0.2072;
		vBinVals[13][4][4] = 0.2072;
		vBinVals[13][4][5] = 0.2072;
		vBinVals[13][4][6] = 0.2072;
		vBinVals[13][4][7] = 0.2072;
		vBinVals[13][4][8] = 0.2072;
		vBinVals[13][4][9] = 0.2072;
		vBinVals[13][4][10] = 0.2072;
		vBinVals[13][4][11] = 0.2072;
		vBinVals[13][4][12] = 0.2072;
		vBinVals[13][4][13] = 0.2072;
		vBinVals[13][4][14] = 0.2072;
		vBinVals[13][4][15] = 0.2072;
		vBinVals[13][4][16] = 0.2072;
		vBinVals[13][4][17] = 0.2072;
		vBinVals[13][4][18] = 0.2072;
		vBinVals[13][4][19] = 0.2072;
		vBinVals[13][4][20] = 0.2072;
		vBinVals[13][4][21] = 0.2072;
		vBinVals[13][4][22] = 0.2072;
		vBinVals[13][4][23] = 0.0181796;
		vBinVals[13][4][24] = 0.0181796;
		vBinVals[13][4][25] = 0.0181796;
		vBinVals[13][4][26] = 0.0181796;
		vBinVals[13][4][27] = 0.904677;
		vBinVals[13][4][28] = 0.904677;
		vBinVals[13][4][29] = 0.904677;
		vBinVals[13][4][30] = 0.904677;
		vBinVals[13][4][31] = 0.904677;
		vBinVals[13][4][32] = 0.374327;
		vBinVals[13][4][33] = 0.374327;
		vBinVals[13][4][34] = 0.374327;
		vBinVals[13][4][35] = 0.374327;
		vBinVals[13][4][36] = 0.374327;
		vBinVals[13][4][37] = 0.240616;
		vBinVals[13][4][38] = 0.240616;
		vBinVals[13][4][39] = 0.240616;
		vBinVals[13][4][40] = 0.240616;
		vBinVals[13][4][41] = 0.240616;
		vBinVals[13][4][42] = 0.166359;
		vBinVals[13][4][43] = 0.166359;
		vBinVals[13][4][44] = 0.166359;
		vBinVals[13][4][45] = 0.166359;
		vBinVals[13][4][46] = 0.166359;
		vBinVals[13][4][47] = 0.180618;
		vBinVals[13][4][48] = 0.180618;
		vBinVals[13][4][49] = 0.180618;
		vBinVals[13][4][50] = 0.180618;
		vBinVals[13][4][51] = 0.180618;
		vBinVals[13][4][52] = 0.34357;
		vBinVals[13][4][53] = 0.34357;
		vBinVals[13][4][54] = 0.34357;
		vBinVals[13][4][55] = 0.34357;
		vBinVals[13][4][56] = 0.34357;
		vBinVals[13][4][57] = 0.0958685;
		vBinVals[13][4][58] = 0.0958685;
		vBinVals[13][4][59] = 0.0958685;
		vBinVals[14][0][0] = 0;
		vBinVals[14][0][1] = 0;
		vBinVals[14][0][2] = 0;
		vBinVals[14][0][3] = 0;
		vBinVals[14][0][4] = 0;
		vBinVals[14][0][5] = 0;
		vBinVals[14][0][6] = 0;
		vBinVals[14][0][7] = 0;
		vBinVals[14][0][8] = 0;
		vBinVals[14][0][9] = 0;
		vBinVals[14][0][10] = 0;
		vBinVals[14][0][11] = 0;
		vBinVals[14][0][12] = 0;
		vBinVals[14][0][13] = 0;
		vBinVals[14][0][14] = 0;
		vBinVals[14][0][15] = 0;
		vBinVals[14][0][16] = 0;
		vBinVals[14][0][17] = 0;
		vBinVals[14][0][18] = 0;
		vBinVals[14][0][19] = 0;
		vBinVals[14][0][20] = 0;
		vBinVals[14][0][21] = 0;
		vBinVals[14][0][22] = 0;
		vBinVals[14][0][23] = 0;
		vBinVals[14][0][24] = 0;
		vBinVals[14][0][25] = 0;
		vBinVals[14][0][26] = 0;
		vBinVals[14][0][27] = 0;
		vBinVals[14][0][28] = 0;
		vBinVals[14][0][29] = 0;
		vBinVals[14][0][30] = 0;
		vBinVals[14][0][31] = 0;
		vBinVals[14][0][32] = 0;
		vBinVals[14][0][33] = 0;
		vBinVals[14][0][34] = 0;
		vBinVals[14][0][35] = 0;
		vBinVals[14][0][36] = 0;
		vBinVals[14][0][37] = 0;
		vBinVals[14][0][38] = 0;
		vBinVals[14][0][39] = 0;
		vBinVals[14][0][40] = 0;
		vBinVals[14][0][41] = 0;
		vBinVals[14][0][42] = 0;
		vBinVals[14][0][43] = 0;
		vBinVals[14][0][44] = 0;
		vBinVals[14][0][45] = 0;
		vBinVals[14][0][46] = 0;
		vBinVals[14][0][47] = 0;
		vBinVals[14][0][48] = 0;
		vBinVals[14][0][49] = 0;
		vBinVals[14][0][50] = 0;
		vBinVals[14][0][51] = 0;
		vBinVals[14][0][52] = 0;
		vBinVals[14][0][53] = 0;
		vBinVals[14][0][54] = 0;
		vBinVals[14][0][55] = 0;
		vBinVals[14][0][56] = 0;
		vBinVals[14][0][57] = 0;
		vBinVals[14][0][58] = 0;
		vBinVals[14][0][59] = 0;
		vBinVals[14][1][0] = 0;
		vBinVals[14][1][1] = 0;
		vBinVals[14][1][2] = 0;
		vBinVals[14][1][3] = 0;
		vBinVals[14][1][4] = 0;
		vBinVals[14][1][5] = 0;
		vBinVals[14][1][6] = 0;
		vBinVals[14][1][7] = 0;
		vBinVals[14][1][8] = 0;
		vBinVals[14][1][9] = 0;
		vBinVals[14][1][10] = 0;
		vBinVals[14][1][11] = 0;
		vBinVals[14][1][12] = 0;
		vBinVals[14][1][13] = 0;
		vBinVals[14][1][14] = 0;
		vBinVals[14][1][15] = 0;
		vBinVals[14][1][16] = 0;
		vBinVals[14][1][17] = 0;
		vBinVals[14][1][18] = 0;
		vBinVals[14][1][19] = 0;
		vBinVals[14][1][20] = 0;
		vBinVals[14][1][21] = 0;
		vBinVals[14][1][22] = 0;
		vBinVals[14][1][23] = 0;
		vBinVals[14][1][24] = 0;
		vBinVals[14][1][25] = 0;
		vBinVals[14][1][26] = 0;
		vBinVals[14][1][27] = 0;
		vBinVals[14][1][28] = 0;
		vBinVals[14][1][29] = 0;
		vBinVals[14][1][30] = 0;
		vBinVals[14][1][31] = 0;
		vBinVals[14][1][32] = 0;
		vBinVals[14][1][33] = 0;
		vBinVals[14][1][34] = 0;
		vBinVals[14][1][35] = 0;
		vBinVals[14][1][36] = 0;
		vBinVals[14][1][37] = 0;
		vBinVals[14][1][38] = 0;
		vBinVals[14][1][39] = 0;
		vBinVals[14][1][40] = 0;
		vBinVals[14][1][41] = 0;
		vBinVals[14][1][42] = 0;
		vBinVals[14][1][43] = 0;
		vBinVals[14][1][44] = 0;
		vBinVals[14][1][45] = 0;
		vBinVals[14][1][46] = 0;
		vBinVals[14][1][47] = 0;
		vBinVals[14][1][48] = 0;
		vBinVals[14][1][49] = 0;
		vBinVals[14][1][50] = 0;
		vBinVals[14][1][51] = 0;
		vBinVals[14][1][52] = 0;
		vBinVals[14][1][53] = 0;
		vBinVals[14][1][54] = 0;
		vBinVals[14][1][55] = 0;
		vBinVals[14][1][56] = 0;
		vBinVals[14][1][57] = 0;
		vBinVals[14][1][58] = 0;
		vBinVals[14][1][59] = 0;
		vBinVals[14][2][0] = -0.268104;
		vBinVals[14][2][1] = -0.268104;
		vBinVals[14][2][2] = -0.268104;
		vBinVals[14][2][3] = -0.268104;
		vBinVals[14][2][4] = -0.268104;
		vBinVals[14][2][5] = -0.268104;
		vBinVals[14][2][6] = -0.268104;
		vBinVals[14][2][7] = -0.268104;
		vBinVals[14][2][8] = -0.268104;
		vBinVals[14][2][9] = -0.268104;
		vBinVals[14][2][10] = -0.268104;
		vBinVals[14][2][11] = -0.268104;
		vBinVals[14][2][12] = -0.268104;
		vBinVals[14][2][13] = -0.268104;
		vBinVals[14][2][14] = -0.268104;
		vBinVals[14][2][15] = -0.268104;
		vBinVals[14][2][16] = -0.268104;
		vBinVals[14][2][17] = -0.268104;
		vBinVals[14][2][18] = -0.268104;
		vBinVals[14][2][19] = -0.268104;
		vBinVals[14][2][20] = -0.268104;
		vBinVals[14][2][21] = -0.268104;
		vBinVals[14][2][22] = -0.268104;
		vBinVals[14][2][23] = -0.268104;
		vBinVals[14][2][24] = -0.268104;
		vBinVals[14][2][25] = -0.268104;
		vBinVals[14][2][26] = -0.268104;
		vBinVals[14][2][27] = -0.268104;
		vBinVals[14][2][28] = -0.268104;
		vBinVals[14][2][29] = -0.268104;
		vBinVals[14][2][30] = -0.268104;
		vBinVals[14][2][31] = -0.268104;
		vBinVals[14][2][32] = -0.268104;
		vBinVals[14][2][33] = -0.268104;
		vBinVals[14][2][34] = -0.268104;
		vBinVals[14][2][35] = -0.268104;
		vBinVals[14][2][36] = -0.268104;
		vBinVals[14][2][37] = -0.323294;
		vBinVals[14][2][38] = -0.323294;
		vBinVals[14][2][39] = -0.323294;
		vBinVals[14][2][40] = -0.323294;
		vBinVals[14][2][41] = -0.323294;
		vBinVals[14][2][42] = -0.316823;
		vBinVals[14][2][43] = -0.316823;
		vBinVals[14][2][44] = -0.316823;
		vBinVals[14][2][45] = -0.316823;
		vBinVals[14][2][46] = -0.316823;
		vBinVals[14][2][47] = -0.292118;
		vBinVals[14][2][48] = -0.292118;
		vBinVals[14][2][49] = -0.292118;
		vBinVals[14][2][50] = -0.292118;
		vBinVals[14][2][51] = -0.292118;
		vBinVals[14][2][52] = -0.268128;
		vBinVals[14][2][53] = -0.268128;
		vBinVals[14][2][54] = -0.268128;
		vBinVals[14][2][55] = -0.268128;
		vBinVals[14][2][56] = -0.268128;
		vBinVals[14][2][57] = -0.247029;
		vBinVals[14][2][58] = -0.247029;
		vBinVals[14][2][59] = -0.247029;
		vBinVals[14][3][0] = 0.0912038;
		vBinVals[14][3][1] = 0.0912038;
		vBinVals[14][3][2] = 0.0912038;
		vBinVals[14][3][3] = 0.0912038;
		vBinVals[14][3][4] = 0.0912038;
		vBinVals[14][3][5] = 0.0912038;
		vBinVals[14][3][6] = 0.0912038;
		vBinVals[14][3][7] = 0.0912038;
		vBinVals[14][3][8] = 0.0912038;
		vBinVals[14][3][9] = 0.0912038;
		vBinVals[14][3][10] = 0.0912038;
		vBinVals[14][3][11] = 0.0912038;
		vBinVals[14][3][12] = 0.0912038;
		vBinVals[14][3][13] = 0.0912038;
		vBinVals[14][3][14] = 0.0912038;
		vBinVals[14][3][15] = 0.0912038;
		vBinVals[14][3][16] = 0.0912038;
		vBinVals[14][3][17] = 0.0912038;
		vBinVals[14][3][18] = 0.0912038;
		vBinVals[14][3][19] = 0.0912038;
		vBinVals[14][3][20] = 0.0912038;
		vBinVals[14][3][21] = 0.0912038;
		vBinVals[14][3][22] = 0.0912038;
		vBinVals[14][3][23] = 0.0912038;
		vBinVals[14][3][24] = 0.0912038;
		vBinVals[14][3][25] = 0.0912038;
		vBinVals[14][3][26] = 0.0912038;
		vBinVals[14][3][27] = 0.0912038;
		vBinVals[14][3][28] = 0.0912038;
		vBinVals[14][3][29] = 0.0912038;
		vBinVals[14][3][30] = 0.0912038;
		vBinVals[14][3][31] = 0.0912038;
		vBinVals[14][3][32] = 0.0912038;
		vBinVals[14][3][33] = 0.0912038;
		vBinVals[14][3][34] = 0.0912038;
		vBinVals[14][3][35] = 0.0912038;
		vBinVals[14][3][36] = 0.0912038;
		vBinVals[14][3][37] = 0.0545128;
		vBinVals[14][3][38] = 0.0545128;
		vBinVals[14][3][39] = 0.0545128;
		vBinVals[14][3][40] = 0.0545128;
		vBinVals[14][3][41] = 0.0545128;
		vBinVals[14][3][42] = 0.0566287;
		vBinVals[14][3][43] = 0.0566287;
		vBinVals[14][3][44] = 0.0566287;
		vBinVals[14][3][45] = 0.0566287;
		vBinVals[14][3][46] = 0.0566287;
		vBinVals[14][3][47] = 0.0661524;
		vBinVals[14][3][48] = 0.0661524;
		vBinVals[14][3][49] = 0.0661524;
		vBinVals[14][3][50] = 0.0661524;
		vBinVals[14][3][51] = 0.0661524;
		vBinVals[14][3][52] = 0.0798734;
		vBinVals[14][3][53] = 0.0798734;
		vBinVals[14][3][54] = 0.0798734;
		vBinVals[14][3][55] = 0.0798734;
		vBinVals[14][3][56] = 0.0798734;
		vBinVals[14][3][57] = 0.0840601;
		vBinVals[14][3][58] = 0.0840601;
		vBinVals[14][3][59] = 0.0840601;
		vBinVals[14][4][0] = 0;
		vBinVals[14][4][1] = 0;
		vBinVals[14][4][2] = 0;
		vBinVals[14][4][3] = 0;
		vBinVals[14][4][4] = 0;
		vBinVals[14][4][5] = 0;
		vBinVals[14][4][6] = 0;
		vBinVals[14][4][7] = 0;
		vBinVals[14][4][8] = 0;
		vBinVals[14][4][9] = 0;
		vBinVals[14][4][10] = 0;
		vBinVals[14][4][11] = 0;
		vBinVals[14][4][12] = 0;
		vBinVals[14][4][13] = 0;
		vBinVals[14][4][14] = 0;
		vBinVals[14][4][15] = 0;
		vBinVals[14][4][16] = 0;
		vBinVals[14][4][17] = 0;
		vBinVals[14][4][18] = 0;
		vBinVals[14][4][19] = 0;
		vBinVals[14][4][20] = 0;
		vBinVals[14][4][21] = 0;
		vBinVals[14][4][22] = 0;
		vBinVals[14][4][23] = 0;
		vBinVals[14][4][24] = 0;
		vBinVals[14][4][25] = 0;
		vBinVals[14][4][26] = 0;
		vBinVals[14][4][27] = 0;
		vBinVals[14][4][28] = 0;
		vBinVals[14][4][29] = 0;
		vBinVals[14][4][30] = 0;
		vBinVals[14][4][31] = 0;
		vBinVals[14][4][32] = 0;
		vBinVals[14][4][33] = 0;
		vBinVals[14][4][34] = 0;
		vBinVals[14][4][35] = 0;
		vBinVals[14][4][36] = 0;
		vBinVals[14][4][37] = 0;
		vBinVals[14][4][38] = 0;
		vBinVals[14][4][39] = 0;
		vBinVals[14][4][40] = 0;
		vBinVals[14][4][41] = 0;
		vBinVals[14][4][42] = 0;
		vBinVals[14][4][43] = 0;
		vBinVals[14][4][44] = 0;
		vBinVals[14][4][45] = 0;
		vBinVals[14][4][46] = 0;
		vBinVals[14][4][47] = 0;
		vBinVals[14][4][48] = 0;
		vBinVals[14][4][49] = 0;
		vBinVals[14][4][50] = 0;
		vBinVals[14][4][51] = 0;
		vBinVals[14][4][52] = 0;
		vBinVals[14][4][53] = 0;
		vBinVals[14][4][54] = 0;
		vBinVals[14][4][55] = 0;
		vBinVals[14][4][56] = 0;
		vBinVals[14][4][57] = 0;
		vBinVals[14][4][58] = 0;
		vBinVals[14][4][59] = 0;
		return vBinVals[iEtaBin][iPar][iBin];
	} else {
		std::cout << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__
			<< ": Undefined value requested. please check. exiting!" << std::endl;
		exit (1);

	}
	return 0;
}






//________________ returns index for eta_bin in JER
int WhatJEReta(double eta_det)
{
	int eta_bin;
	for(int i=0; i<14; i++)
	{
		if(fabs(eta_det)>=i*0.2 && fabs(eta_det)<(i+1)*0.2)
		{
			eta_bin=i;
			break;
		}
	}
	if(fabs(eta_det)>=2.8) eta_bin=14;
	return eta_bin;
}












//----- 2nd round fit functions ---------------------------------------------------------------
double new3JERparam(int code, int i, int j)
{
	double new3jer_param[15][31]; // [eta][param]
	double new3jer_param_er[15][31]; // [eta][param]
	double value=0.0;

	//========================== NEW FIT FUNCTION PARAMETERS ========= 08-05-2009
	// I limit the step function crossing point in an attempt to overcome the problem
	// mentioned in elog#1270

  new3jer_param[0][0]=-0.107537;  new3jer_param_er[0][0]=0.0566252;
  new3jer_param[0][1]=0.0101998;  new3jer_param_er[0][1]=0.00267699;
  new3jer_param[0][2]=-9.99535e-05;  new3jer_param_er[0][2]=3.24013e-05;
  new3jer_param[0][3]=0.116372;  new3jer_param_er[0][3]=0.00743799;
  new3jer_param[0][4]=-0.000824633;  new3jer_param_er[0][4]=5.77058e-05;
  new3jer_param[0][5]=1.23834e-06;  new3jer_param_er[0][5]=1.02365e-07;
  new3jer_param[0][6]=50;  new3jer_param_er[0][6]=1.57999;
  new3jer_param[0][7]=-1.10793;  new3jer_param_er[0][7]=1.0448;
  new3jer_param[0][8]=0.0121459;  new3jer_param_er[0][8]=0.0347651;
  new3jer_param[0][9]=7.9828;  new3jer_param_er[0][9]=0.284825;
  new3jer_param[0][10]=0.0339306;  new3jer_param_er[0][10]=0.0013341;
  new3jer_param[0][11]=20;  new3jer_param_er[0][11]=15.6888;
  new3jer_param[0][12]=-0.0570722;  new3jer_param_er[0][12]=0.0169371;
  new3jer_param[0][13]=0.00255033;  new3jer_param_er[0][13]=0.000750679;
  new3jer_param[0][14]=-1.95879e-05;  new3jer_param_er[0][14]=1.37464e-05;
  new3jer_param[0][15]=0.0518506;  new3jer_param_er[0][15]=0.0110479;
  new3jer_param[0][16]=-0.00015948;  new3jer_param_er[0][16]=0.000100111;
  new3jer_param[0][17]=6.38365e-07;  new3jer_param_er[0][17]=2.069e-07;
  new3jer_param[0][18]=70;  new3jer_param_er[0][18]=16.8986;
  new3jer_param[0][19]=1.35107;  new3jer_param_er[0][19]=0.259805;
  new3jer_param[0][20]=0.0686707;  new3jer_param_er[0][20]=0.0165325;
  new3jer_param[0][21]=0.793687;  new3jer_param_er[0][21]=0.279586;
  new3jer_param[0][22]=0.0874214;  new3jer_param_er[0][22]=0.00172322;
  new3jer_param[0][23]=20;  new3jer_param_er[0][23]=17.4647;
  new3jer_param[0][24]=-0.107059;  new3jer_param_er[0][24]=0.0477939;
  new3jer_param[0][25]=0.00973563;  new3jer_param_er[0][25]=0.00298643;
  new3jer_param[0][26]=-9.48027e-05;  new3jer_param_er[0][26]=4.79747e-05;
  new3jer_param[0][27]=0.165775;  new3jer_param_er[0][27]=0.028031;
  new3jer_param[0][28]=-0.000941991;  new3jer_param_er[0][28]=0.000259297;
  new3jer_param[0][29]=2.93026e-06;  new3jer_param_er[0][29]=5.5319e-07;
  new3jer_param[0][30]=50;  new3jer_param_er[0][30]=12.3864;
  new3jer_param[1][0]=0.0884825;  new3jer_param_er[1][0]=0.0165723;
  new3jer_param[1][1]=-0.000396748;  new3jer_param_er[1][1]=0.000668837;
  new3jer_param[1][2]=9.00949e-06;  new3jer_param_er[1][2]=6.82975e-06;
  new3jer_param[1][3]=0.0721112;  new3jer_param_er[1][3]=0.00530366;
  new3jer_param[1][4]=-2.41626e-05;  new3jer_param_er[1][4]=4.58352e-05;
  new3jer_param[1][5]=-3.27678e-07;  new3jer_param_er[1][5]=9.01269e-08;
  new3jer_param[1][6]=70;  new3jer_param_er[1][6]=2.42223;
  new3jer_param[1][7]=1.01944;  new3jer_param_er[1][7]=0.426821;
  new3jer_param[1][8]=0.0256876;  new3jer_param_er[1][8]=0.0244477;
  new3jer_param[1][9]=4.63705;  new3jer_param_er[1][9]=0.188783;
  new3jer_param[1][10]=0.0567771;  new3jer_param_er[1][10]=0.00106819;
  new3jer_param[1][11]=20;  new3jer_param_er[1][11]=18.9918;
  new3jer_param[1][12]=-0.136141;  new3jer_param_er[1][12]=0.0147021;
  new3jer_param[1][13]=0.0046606;  new3jer_param_er[1][13]=0.000921656;
  new3jer_param[1][14]=-4.67532e-05;  new3jer_param_er[1][14]=1.28932e-05;
  new3jer_param[1][15]=-0.0192809;  new3jer_param_er[1][15]=0.00423984;
  new3jer_param[1][16]=0.000115645;  new3jer_param_er[1][16]=3.40824e-05;
  new3jer_param[1][17]=-1.86755e-07;  new3jer_param_er[1][17]=6.37078e-08;
  new3jer_param[1][18]=60.2651;  new3jer_param_er[1][18]=7.51143;
  new3jer_param[1][19]=1.10845;  new3jer_param_er[1][19]=0.119022;
  new3jer_param[1][20]=0.0566205;  new3jer_param_er[1][20]=0.0105459;
  new3jer_param[1][21]=0.812854;  new3jer_param_er[1][21]=0.139593;
  new3jer_param[1][22]=0.0483472;  new3jer_param_er[1][22]=0.00066212;
  new3jer_param[1][23]=20;  new3jer_param_er[1][23]=18.4786;
  new3jer_param[1][24]=-0.0442165;  new3jer_param_er[1][24]=0.0385695;
  new3jer_param[1][25]=0.0106688;  new3jer_param_er[1][25]=0.00305846;
  new3jer_param[1][26]=-0.000127049;  new3jer_param_er[1][26]=5.68601e-05;
  new3jer_param[1][27]=0.387361;  new3jer_param_er[1][27]=0.0194439;
  new3jer_param[1][28]=-0.00172015;  new3jer_param_er[1][28]=0.000150375;
  new3jer_param[1][29]=2.39865e-06;  new3jer_param_er[1][29]=2.89995e-07;
  new3jer_param[1][30]=50;  new3jer_param_er[1][30]=15.8042;
  new3jer_param[2][0]=-0.0558391;  new3jer_param_er[2][0]=0.0286116;
  new3jer_param[2][1]=0.0047738;  new3jer_param_er[2][1]=0.00145382;
  new3jer_param[2][2]=-4.14672e-05;  new3jer_param_er[2][2]=1.83055e-05;
  new3jer_param[2][3]=0.0985458;  new3jer_param_er[2][3]=0.00479701;
  new3jer_param[2][4]=-0.000310602;  new3jer_param_er[2][4]=3.95796e-05;
  new3jer_param[2][5]=1.84186e-07;  new3jer_param_er[2][5]=7.44813e-08;
  new3jer_param[2][6]=50;  new3jer_param_er[2][6]=19.7341;
  new3jer_param[2][7]=0.512322;  new3jer_param_er[2][7]=0.546257;
  new3jer_param[2][8]=0.00652577;  new3jer_param_er[2][8]=0.0312757;
  new3jer_param[2][9]=5.22907;  new3jer_param_er[2][9]=0.1937;
  new3jer_param[2][10]=0.0525518;  new3jer_param_er[2][10]=0.00107728;
  new3jer_param[2][11]=20;  new3jer_param_er[2][11]=18.2862;
  new3jer_param[2][12]=-0.137934;  new3jer_param_er[2][12]=0.025341;
  new3jer_param[2][13]=0.00307188;  new3jer_param_er[2][13]=0.00140956;
  new3jer_param[2][14]=-2.02804e-05;  new3jer_param_er[2][14]=1.73514e-05;
  new3jer_param[2][15]=-0.0246058;  new3jer_param_er[2][15]=0.00573657;
  new3jer_param[2][16]=0.000142928;  new3jer_param_er[2][16]=4.31133e-05;
  new3jer_param[2][17]=-2.64319e-07;  new3jer_param_er[2][17]=7.49754e-08;
  new3jer_param[2][18]=50;  new3jer_param_er[2][18]=17.9769;
  new3jer_param[2][19]=1.48165;  new3jer_param_er[2][19]=0.258778;
  new3jer_param[2][20]=0.0167996;  new3jer_param_er[2][20]=0.0164251;
  new3jer_param[2][21]=1.42852;  new3jer_param_er[2][21]=0.142274;
  new3jer_param[2][22]=0.047546;  new3jer_param_er[2][22]=0.000682538;
  new3jer_param[2][23]=20;  new3jer_param_er[2][23]=19.3318;
  new3jer_param[2][24]=-0.0605025;  new3jer_param_er[2][24]=0.0736183;
  new3jer_param[2][25]=0.0119383;  new3jer_param_er[2][25]=0.0045531;
  new3jer_param[2][26]=-0.00012213;  new3jer_param_er[2][26]=6.53851e-05;
  new3jer_param[2][27]=0.342351;  new3jer_param_er[2][27]=0.0173133;
  new3jer_param[2][28]=-0.00135649;  new3jer_param_er[2][28]=0.000133935;
  new3jer_param[2][29]=1.76527e-06;  new3jer_param_er[2][29]=2.42139e-07;
  new3jer_param[2][30]=50;  new3jer_param_er[2][30]=15.9222;
  new3jer_param[3][0]=0.220049;  new3jer_param_er[3][0]=0.0229101;
  new3jer_param[3][1]=-0.00716605;  new3jer_param_er[3][1]=0.00121414;
  new3jer_param[3][2]=9.03346e-05;  new3jer_param_er[3][2]=1.63624e-05;
  new3jer_param[3][3]=0.0913929;  new3jer_param_er[3][3]=0.0045358;
  new3jer_param[3][4]=-0.000150115;  new3jer_param_er[3][4]=3.99071e-05;
  new3jer_param[3][5]=-7.54507e-08;  new3jer_param_er[3][5]=8.07302e-08;
  new3jer_param[3][6]=50;  new3jer_param_er[3][6]=0.912723;
  new3jer_param[3][7]=1.84997;  new3jer_param_er[3][7]=0.271656;
  new3jer_param[3][8]=0.0535898;  new3jer_param_er[3][8]=0.0167578;
  new3jer_param[3][9]=4.99511;  new3jer_param_er[3][9]=0.219681;
  new3jer_param[3][10]=0.0588264;  new3jer_param_er[3][10]=0.00124432;
  new3jer_param[3][11]=40;  new3jer_param_er[3][11]=15.8004;
  new3jer_param[3][12]=-0.184246;  new3jer_param_er[3][12]=0.00930048;
  new3jer_param[3][13]=0.00538507;  new3jer_param_er[3][13]=0.000441317;
  new3jer_param[3][14]=-4.49893e-05;  new3jer_param_er[3][14]=5.23713e-06;
  new3jer_param[3][15]=-0.000272302;  new3jer_param_er[3][15]=0.00402049;
  new3jer_param[3][16]=6.88301e-05;  new3jer_param_er[3][16]=3.04064e-05;
  new3jer_param[3][17]=-1.80491e-07;  new3jer_param_er[3][17]=5.44626e-08;
  new3jer_param[3][18]=70;  new3jer_param_er[3][18]=12.891;
  new3jer_param[3][19]=0.743505;  new3jer_param_er[3][19]=0.192455;
  new3jer_param[3][20]=0.0469858;  new3jer_param_er[3][20]=0.00935696;
  new3jer_param[3][21]=1.92447;  new3jer_param_er[3][21]=0.14965;
  new3jer_param[3][22]=0.048089;  new3jer_param_er[3][22]=0.00073469;
  new3jer_param[3][23]=40;  new3jer_param_er[3][23]=17.2114;
  new3jer_param[3][24]=0.0703947;  new3jer_param_er[3][24]=0.0316273;
  new3jer_param[3][25]=0.00104987;  new3jer_param_er[3][25]=0.002421;
  new3jer_param[3][26]=3.07433e-05;  new3jer_param_er[3][26]=4.11166e-05;
  new3jer_param[3][27]=0.275658;  new3jer_param_er[3][27]=0.0150253;
  new3jer_param[3][28]=-0.000975134;  new3jer_param_er[3][28]=0.000116931;
  new3jer_param[3][29]=1.15296e-06;  new3jer_param_er[3][29]=2.15012e-07;
  new3jer_param[3][30]=50;  new3jer_param_er[3][30]=3.06942;
  new3jer_param[4][0]=0.130408;  new3jer_param_er[4][0]=0.0222449;
  new3jer_param[4][1]=-0.00340213;  new3jer_param_er[4][1]=0.00129552;
  new3jer_param[4][2]=5.74043e-05;  new3jer_param_er[4][2]=1.81035e-05;
  new3jer_param[4][3]=0.122659;  new3jer_param_er[4][3]=0.00503038;
  new3jer_param[4][4]=-0.000292842;  new3jer_param_er[4][4]=4.61727e-05;
  new3jer_param[4][5]=3.60169e-07;  new3jer_param_er[4][5]=9.6853e-08;
  new3jer_param[4][6]=50;  new3jer_param_er[4][6]=12.3584;
  new3jer_param[4][7]=0.0994227;  new3jer_param_er[4][7]=0.502704;
  new3jer_param[4][8]=0.0864482;  new3jer_param_er[4][8]=0.0177311;
  new3jer_param[4][9]=4.44844;  new3jer_param_er[4][9]=0.221593;
  new3jer_param[4][10]=0.0786389;  new3jer_param_er[4][10]=0.00128124;
  new3jer_param[4][11]=20;  new3jer_param_er[4][11]=18.4469;
  new3jer_param[4][12]=-0.165712;  new3jer_param_er[4][12]=0.0240303;
  new3jer_param[4][13]=0.00507425;  new3jer_param_er[4][13]=0.00153825;
  new3jer_param[4][14]=-4.33962e-05;  new3jer_param_er[4][14]=2.58543e-05;
  new3jer_param[4][15]=0.00596175;  new3jer_param_er[4][15]=0.00760369;
  new3jer_param[4][16]=9.50617e-05;  new3jer_param_er[4][16]=5.59348e-05;
  new3jer_param[4][17]=-2.61393e-07;  new3jer_param_er[4][17]=9.77872e-08;
  new3jer_param[4][18]=70;  new3jer_param_er[4][18]=19.9803;
  new3jer_param[4][19]=0.807587;  new3jer_param_er[4][19]=0.108953;
  new3jer_param[4][20]=0.0527689;  new3jer_param_er[4][20]=0.00652034;
  new3jer_param[4][21]=2.06944;  new3jer_param_er[4][21]=0.191477;
  new3jer_param[4][22]=0.05748;  new3jer_param_er[4][22]=0.00089954;
  new3jer_param[4][23]=40;  new3jer_param_er[4][23]=13.6997;
  new3jer_param[4][24]=0.143285;  new3jer_param_er[4][24]=0.0409534;
  new3jer_param[4][25]=-0.0037696;  new3jer_param_er[4][25]=0.00300876;
  new3jer_param[4][26]=7.85254e-05;  new3jer_param_er[4][26]=4.86225e-05;
  new3jer_param[4][27]=0.321488;  new3jer_param_er[4][27]=0.0186633;
  new3jer_param[4][28]=-0.000971977;  new3jer_param_er[4][28]=0.000144944;
  new3jer_param[4][29]=9.86957e-07;  new3jer_param_er[4][29]=2.66791e-07;
  new3jer_param[4][30]=50;  new3jer_param_er[4][30]=2.62261;
  new3jer_param[5][0]=0.0735202;  new3jer_param_er[5][0]=0.0219619;
  new3jer_param[5][1]=0.00018131;  new3jer_param_er[5][1]=0.000842598;
  new3jer_param[5][2]=-4.52704e-06;  new3jer_param_er[5][2]=8.93316e-06;
  new3jer_param[5][3]=0.0908279;  new3jer_param_er[5][3]=0.00405447;
  new3jer_param[5][4]=-0.000112972;  new3jer_param_er[5][4]=3.28994e-05;
  new3jer_param[5][5]=6.59606e-08;  new3jer_param_er[5][5]=6.13004e-08;
  new3jer_param[5][6]=70;  new3jer_param_er[5][6]=14.2521;
  new3jer_param[5][7]=1.52183;  new3jer_param_er[5][7]=0.389597;
  new3jer_param[5][8]=0.102812;  new3jer_param_er[5][8]=0.0156864;
  new3jer_param[5][9]=3.17222;  new3jer_param_er[5][9]=0.181142;
  new3jer_param[5][10]=0.084634;  new3jer_param_er[5][10]=0.000934544;
  new3jer_param[5][11]=20;  new3jer_param_er[5][11]=15.5553;
  new3jer_param[5][12]=-0.301711;  new3jer_param_er[5][12]=0.0292409;
  new3jer_param[5][13]=0.0128286;  new3jer_param_er[5][13]=0.00169054;
  new3jer_param[5][14]=-0.000181522;  new3jer_param_er[5][14]=2.39601e-05;
  new3jer_param[5][15]=-0.0092045;  new3jer_param_er[5][15]=0.0053185;
  new3jer_param[5][16]=6.96301e-05;  new3jer_param_er[5][16]=4.17049e-05;
  new3jer_param[5][17]=-1.83794e-07;  new3jer_param_er[5][17]=7.56229e-08;
  new3jer_param[5][18]=50;  new3jer_param_er[5][18]=1.90508;
  new3jer_param[5][19]=1.95663;  new3jer_param_er[5][19]=0.182142;
  new3jer_param[5][20]=0.00437304;  new3jer_param_er[5][20]=0.00930281;
  new3jer_param[5][21]=2.23237;  new3jer_param_er[5][21]=0.246507;
  new3jer_param[5][22]=0.0514385;  new3jer_param_er[5][22]=0.00109073;
  new3jer_param[5][23]=40;  new3jer_param_er[5][23]=2.13171;
  new3jer_param[5][24]=0.178338;  new3jer_param_er[5][24]=0.0640275;
  new3jer_param[5][25]=-0.00319376;  new3jer_param_er[5][25]=0.0038986;
  new3jer_param[5][26]=5.59603e-05;  new3jer_param_er[5][26]=6.59878e-05;
  new3jer_param[5][27]=0.508882;  new3jer_param_er[5][27]=0.0433818;
  new3jer_param[5][28]=-0.000933784;  new3jer_param_er[5][28]=0.000312936;
  new3jer_param[5][29]=4.68644e-07;  new3jer_param_er[5][29]=5.33297e-07;
  new3jer_param[5][30]=68.3458;  new3jer_param_er[5][30]=11.8792;
  new3jer_param[6][0]=0.0347005;  new3jer_param_er[6][0]=0.0184864;
  new3jer_param[6][1]=0.000937425;  new3jer_param_er[6][1]=0.000830865;
  new3jer_param[6][2]=-1.30709e-05;  new3jer_param_er[6][2]=9.21559e-06;
  new3jer_param[6][3]=0.0875522;  new3jer_param_er[6][3]=0.00888453;
  new3jer_param[6][4]=-8.55151e-05;  new3jer_param_er[6][4]=8.18485e-05;
  new3jer_param[6][5]=-4.6158e-07;  new3jer_param_er[6][5]=1.52946e-07;
  new3jer_param[6][6]=70;  new3jer_param_er[6][6]=4.03686;
  new3jer_param[6][7]=2.34387;  new3jer_param_er[6][7]=0.186198;
  new3jer_param[6][8]=-0.0335922;  new3jer_param_er[6][8]=0.0115086;
  new3jer_param[6][9]=8.35159;  new3jer_param_er[6][9]=0.348084;
  new3jer_param[6][10]=0.0485484;  new3jer_param_er[6][10]=0.00217528;
  new3jer_param[6][11]=40;  new3jer_param_er[6][11]=1.03671;
  new3jer_param[6][12]=-0.344393;  new3jer_param_er[6][12]=0.0237024;
  new3jer_param[6][13]=0.0131182;  new3jer_param_er[6][13]=0.00148083;
  new3jer_param[6][14]=-0.000180064;  new3jer_param_er[6][14]=2.18422e-05;
  new3jer_param[6][15]=-0.0105614;  new3jer_param_er[6][15]=0.00699589;
  new3jer_param[6][16]=-2.36259e-05;  new3jer_param_er[6][16]=7.01264e-05;
  new3jer_param[6][17]=2.83426e-07;  new3jer_param_er[6][17]=1.56341e-07;
  new3jer_param[6][18]=50;  new3jer_param_er[6][18]=1.79359;
  new3jer_param[6][19]=1.31123;  new3jer_param_er[6][19]=0.307076;
  new3jer_param[6][20]=0.0532645;  new3jer_param_er[6][20]=0.0152332;
  new3jer_param[6][21]=-0.0297589;  new3jer_param_er[6][21]=0.275428;
  new3jer_param[6][22]=0.0689466;  new3jer_param_er[6][22]=0.0014889;
  new3jer_param[6][23]=40;  new3jer_param_er[6][23]=19.9687;
  new3jer_param[6][24]=0.325285;  new3jer_param_er[6][24]=0.0729539;
  new3jer_param[6][25]=-0.0133924;  new3jer_param_er[6][25]=0.00454512;
  new3jer_param[6][26]=0.000190133;  new3jer_param_er[6][26]=6.429e-05;
  new3jer_param[6][27]=0.329918;  new3jer_param_er[6][27]=0.0247353;
  new3jer_param[6][28]=-0.00194938;  new3jer_param_er[6][28]=0.00023475;
  new3jer_param[6][29]=3.39452e-06;  new3jer_param_er[6][29]=5.05375e-07;
  new3jer_param[6][30]=50;  new3jer_param_er[6][30]=1.63864;
  new3jer_param[7][0]=-0.0363087;  new3jer_param_er[7][0]=0.0159478;
  new3jer_param[7][1]=0.00256729;  new3jer_param_er[7][1]=0.000771292;
  new3jer_param[7][2]=-2.37768e-05;  new3jer_param_er[7][2]=8.48443e-06;
  new3jer_param[7][3]=0.0582484;  new3jer_param_er[7][3]=0.00712301;
  new3jer_param[7][4]=-0.000189063;  new3jer_param_er[7][4]=6.41661e-05;
  new3jer_param[7][5]=1.15837e-08;  new3jer_param_er[7][5]=1.27971e-07;
  new3jer_param[7][6]=70;  new3jer_param_er[7][6]=1.51805;
  new3jer_param[7][7]=0.509533;  new3jer_param_er[7][7]=0.621444;
  new3jer_param[7][8]=-0.00558115;  new3jer_param_er[7][8]=0.0336057;
  new3jer_param[7][9]=6.59333;  new3jer_param_er[7][9]=0.27519;
  new3jer_param[7][10]=0.0426632;  new3jer_param_er[7][10]=0.00171458;
  new3jer_param[7][11]=20;  new3jer_param_er[7][11]=16.3604;
  new3jer_param[7][12]=-0.296756;  new3jer_param_er[7][12]=0.0150918;
  new3jer_param[7][13]=0.00867224;  new3jer_param_er[7][13]=0.000771177;
  new3jer_param[7][14]=-8.12718e-05;  new3jer_param_er[7][14]=8.72039e-06;
  new3jer_param[7][15]=-0.0690383;  new3jer_param_er[7][15]=0.00880247;
  new3jer_param[7][16]=0.000460047;  new3jer_param_er[7][16]=8.56913e-05;
  new3jer_param[7][17]=-8.53995e-07;  new3jer_param_er[7][17]=1.89517e-07;
  new3jer_param[7][18]=70;  new3jer_param_er[7][18]=0.583276;
  new3jer_param[7][19]=0.111636;  new3jer_param_er[7][19]=0.272265;
  new3jer_param[7][20]=0.169657;  new3jer_param_er[7][20]=0.0221918;
  new3jer_param[7][21]=-0.268551;  new3jer_param_er[7][21]=0.301645;
  new3jer_param[7][22]=0.0486536;  new3jer_param_er[7][22]=0.00194;
  new3jer_param[7][23]=20;  new3jer_param_er[7][23]=3.18675;
  new3jer_param[7][24]=0.345809;  new3jer_param_er[7][24]=0.059104;
  new3jer_param[7][25]=-0.00997728;  new3jer_param_er[7][25]=0.00338212;
  new3jer_param[7][26]=0.000137981;  new3jer_param_er[7][26]=4.14079e-05;
  new3jer_param[7][27]=0.214122;  new3jer_param_er[7][27]=0.0613989;
  new3jer_param[7][28]=0.000120621;  new3jer_param_er[7][28]=0.00062674;
  new3jer_param[7][29]=-9.90705e-07;  new3jer_param_er[7][29]=1.43959e-06;
  new3jer_param[7][30]=70;  new3jer_param_er[7][30]=1.79527;
  new3jer_param[8][0]=-0.045249;  new3jer_param_er[8][0]=0.0240541;
  new3jer_param[8][1]=0.00248962;  new3jer_param_er[8][1]=0.00102856;
  new3jer_param[8][2]=-3.55084e-05;  new3jer_param_er[8][2]=1.00571e-05;
  new3jer_param[8][3]=0.0274564;  new3jer_param_er[8][3]=0.00381056;
  new3jer_param[8][4]=5.67099e-05;  new3jer_param_er[8][4]=3.09473e-05;
  new3jer_param[8][5]=-2.32842e-07;  new3jer_param_er[8][5]=5.68585e-08;
  new3jer_param[8][6]=70;  new3jer_param_er[8][6]=1.60276;
  new3jer_param[8][7]=3.74538;  new3jer_param_er[8][7]=0.31126;
  new3jer_param[8][8]=-0.0350332;  new3jer_param_er[8][8]=0.012865;
  new3jer_param[8][9]=7.09742;  new3jer_param_er[8][9]=0.206276;
  new3jer_param[8][10]=0.0369685;  new3jer_param_er[8][10]=0.00106271;
  new3jer_param[8][11]=40;  new3jer_param_er[8][11]=2.16757;
  new3jer_param[8][12]=-0.224642;  new3jer_param_er[8][12]=0.0481272;
  new3jer_param[8][13]=0.00392954;  new3jer_param_er[8][13]=0.00252884;
  new3jer_param[8][14]=-2.57769e-05;  new3jer_param_er[8][14]=3.28976e-05;
  new3jer_param[8][15]=-0.10305;  new3jer_param_er[8][15]=0.00701373;
  new3jer_param[8][16]=0.000563536;  new3jer_param_er[8][16]=5.94671e-05;
  new3jer_param[8][17]=-8.82646e-07;  new3jer_param_er[8][17]=1.10424e-07;
  new3jer_param[8][18]=70;  new3jer_param_er[8][18]=17.5701;
  new3jer_param[8][19]=0.567367;  new3jer_param_er[8][19]=0.656913;
  new3jer_param[8][20]=0.0704125;  new3jer_param_er[8][20]=0.0305535;
  new3jer_param[8][21]=2.35039;  new3jer_param_er[8][21]=0.274621;
  new3jer_param[8][22]=0.0187985;  new3jer_param_er[8][22]=0.00156133;
  new3jer_param[8][23]=20;  new3jer_param_er[8][23]=16.6277;
  new3jer_param[8][24]=0.422792;  new3jer_param_er[8][24]=0.166925;
  new3jer_param[8][25]=-0.013864;  new3jer_param_er[8][25]=0.0109505;
  new3jer_param[8][26]=0.000171851;  new3jer_param_er[8][26]=0.000191956;
  new3jer_param[8][27]=0.272589;  new3jer_param_er[8][27]=0.0977811;
  new3jer_param[8][28]=0.00155807;  new3jer_param_er[8][28]=0.000883326;
  new3jer_param[8][29]=-3.36749e-06;  new3jer_param_er[8][29]=1.69713e-06;
  new3jer_param[8][30]=70;  new3jer_param_er[8][30]=15.6354;
  new3jer_param[9][0]=-0.174942;  new3jer_param_er[9][0]=0.0410119;
  new3jer_param[9][1]=0.00849332;  new3jer_param_er[9][1]=0.0017247;
  new3jer_param[9][2]=-0.000106572;  new3jer_param_er[9][2]=1.66441e-05;
  new3jer_param[9][3]=0.0340075;  new3jer_param_er[9][3]=0.00655096;
  new3jer_param[9][4]=-3.63072e-05;  new3jer_param_er[9][4]=6.33251e-05;
  new3jer_param[9][5]=1.85558e-08;  new3jer_param_er[9][5]=1.46601e-07;
  new3jer_param[9][6]=70;  new3jer_param_er[9][6]=0.479253;
  new3jer_param[9][7]=6.61446;  new3jer_param_er[9][7]=0.634631;
  new3jer_param[9][8]=-0.105715;  new3jer_param_er[9][8]=0.0230839;
  new3jer_param[9][9]=7.25843;  new3jer_param_er[9][9]=0.286012;
  new3jer_param[9][10]=0.0412497;  new3jer_param_er[9][10]=0.00145571;
  new3jer_param[9][11]=40;  new3jer_param_er[9][11]=1.6096;
  new3jer_param[9][12]=-0.337628;  new3jer_param_er[9][12]=0.0457211;
  new3jer_param[9][13]=0.00796001;  new3jer_param_er[9][13]=0.00208963;
  new3jer_param[9][14]=-6.18984e-05;  new3jer_param_er[9][14]=2.24107e-05;
  new3jer_param[9][15]=-0.108523;  new3jer_param_er[9][15]=0.0114774;
  new3jer_param[9][16]=0.000689633;  new3jer_param_er[9][16]=0.000114841;
  new3jer_param[9][17]=-1.32862e-06;  new3jer_param_er[9][17]=2.70826e-07;
  new3jer_param[9][18]=70;  new3jer_param_er[9][18]=4.01764;
  new3jer_param[9][19]=-0.813436;  new3jer_param_er[9][19]=1.17622;
  new3jer_param[9][20]=0.103229;  new3jer_param_er[9][20]=0.0439382;
  new3jer_param[9][21]=2.49217;  new3jer_param_er[9][21]=0.35577;
  new3jer_param[9][22]=0.0188682;  new3jer_param_er[9][22]=0.00200871;
  new3jer_param[9][23]=20;  new3jer_param_er[9][23]=16.2629;
  new3jer_param[9][24]=0.838697;  new3jer_param_er[9][24]=0.195885;
  new3jer_param[9][25]=-0.0333986;  new3jer_param_er[9][25]=0.00970986;
  new3jer_param[9][26]=0.000361833;  new3jer_param_er[9][26]=0.000116387;
  new3jer_param[9][27]=0.267283;  new3jer_param_er[9][27]=0.118979;
  new3jer_param[9][28]=0.000616332;  new3jer_param_er[9][28]=0.00124207;
  new3jer_param[9][29]=-9.39781e-07;  new3jer_param_er[9][29]=2.97127e-06;
  new3jer_param[9][30]=70;  new3jer_param_er[9][30]=2.85094;
  new3jer_param[10][0]=-0.574168;  new3jer_param_er[10][0]=0.0854421;
  new3jer_param[10][1]=0.0228021;  new3jer_param_er[10][1]=0.00331966;
  new3jer_param[10][2]=-0.000230131;  new3jer_param_er[10][2]=3.14273e-05;
  new3jer_param[10][3]=0.0324899;  new3jer_param_er[10][3]=0.00590435;
  new3jer_param[10][4]=-5.99141e-05;  new3jer_param_er[10][4]=5.00397e-05;
  new3jer_param[10][5]=2.3651e-10;  new3jer_param_er[10][5]=1.00775e-07;
  new3jer_param[10][6]=70;  new3jer_param_er[10][6]=1.22064;
  new3jer_param[10][7]=7.09343;  new3jer_param_er[10][7]=1.46866;
  new3jer_param[10][8]=-0.104924;  new3jer_param_er[10][8]=0.0421128;
  new3jer_param[10][9]=7.73291;  new3jer_param_er[10][9]=0.341977;
  new3jer_param[10][10]=0.041686;  new3jer_param_er[10][10]=0.00162816;
  new3jer_param[10][11]=40;  new3jer_param_er[10][11]=3.64618;
  new3jer_param[10][12]=-0.133786;  new3jer_param_er[10][12]=0.196081;
  new3jer_param[10][13]=0.0015014;  new3jer_param_er[10][13]=0.00912009;
  new3jer_param[10][14]=-3.44832e-05;  new3jer_param_er[10][14]=0.000104836;
  new3jer_param[10][15]=-0.0996914;  new3jer_param_er[10][15]=0.00913848;
  new3jer_param[10][16]=0.000469939;  new3jer_param_er[10][16]=8.31301e-05;
  new3jer_param[10][17]=-7.81904e-07;  new3jer_param_er[10][17]=1.74619e-07;
  new3jer_param[10][18]=50;  new3jer_param_er[10][18]=13.5254;
  new3jer_param[10][19]=1.87346;  new3jer_param_er[10][19]=2.09245;
  new3jer_param[10][20]=0.0496746;  new3jer_param_er[10][20]=0.0647441;
  new3jer_param[10][21]=2.70178;  new3jer_param_er[10][21]=0.385632;
  new3jer_param[10][22]=0.0217253;  new3jer_param_er[10][22]=0.00202821;
  new3jer_param[10][23]=20;  new3jer_param_er[10][23]=10.0071;
  new3jer_param[10][24]=-0.159945;  new3jer_param_er[10][24]=0.279156;
  new3jer_param[10][25]=0.0165115;  new3jer_param_er[10][25]=0.0125563;
  new3jer_param[10][26]=-0.000229064;  new3jer_param_er[10][26]=0.000143604;
  new3jer_param[10][27]=0.364856;  new3jer_param_er[10][27]=0.106014;
  new3jer_param[10][28]=-0.000667039;  new3jer_param_er[10][28]=0.00101174;
  new3jer_param[10][29]=1.66469e-06;  new3jer_param_er[10][29]=2.21821e-06;
  new3jer_param[10][30]=70;  new3jer_param_er[10][30]=12.5585;
  new3jer_param[11][0]=-0.90196;  new3jer_param_er[11][0]=0.515815;
  new3jer_param[11][1]=0.0369416;  new3jer_param_er[11][1]=0.0222684;
  new3jer_param[11][2]=-0.000395616;  new3jer_param_er[11][2]=0.000237496;
  new3jer_param[11][3]=-0.0477924;  new3jer_param_er[11][3]=0.00998877;
  new3jer_param[11][4]=0.000390147;  new3jer_param_er[11][4]=8.47609e-05;
  new3jer_param[11][5]=-6.39862e-07;  new3jer_param_er[11][5]=1.66749e-07;
  new3jer_param[11][6]=50;  new3jer_param_er[11][6]=16.3431;
  new3jer_param[11][7]=2.32859;  new3jer_param_er[11][7]=4.81521;
  new3jer_param[11][8]=-0.00697212;  new3jer_param_er[11][8]=0.122267;
  new3jer_param[11][9]=7.32018;  new3jer_param_er[11][9]=0.42129;
  new3jer_param[11][10]=0.0429903;  new3jer_param_er[11][10]=0.0019489;
  new3jer_param[11][11]=20;  new3jer_param_er[11][11]=19.9772;
  new3jer_param[11][12]=-0.294175;  new3jer_param_er[11][12]=0.323618;
  new3jer_param[11][13]=0.00464302;  new3jer_param_er[11][13]=0.0142533;
  new3jer_param[11][14]=-2.68264e-05;  new3jer_param_er[11][14]=0.000152949;
  new3jer_param[11][15]=-0.155451;  new3jer_param_er[11][15]=0.0148115;
  new3jer_param[11][16]=0.000801954;  new3jer_param_er[11][16]=0.000125203;
  new3jer_param[11][17]=-1.30184e-06;  new3jer_param_er[11][17]=2.44338e-07;
  new3jer_param[11][18]=70;  new3jer_param_er[11][18]=19.0245;
  new3jer_param[11][19]=-1.19025;  new3jer_param_er[11][19]=3.10545;
  new3jer_param[11][20]=0.0875681;  new3jer_param_er[11][20]=0.0887295;
  new3jer_param[11][21]=3.13428;  new3jer_param_er[11][21]=0.499745;
  new3jer_param[11][22]=0.0212547;  new3jer_param_er[11][22]=0.0024322;
  new3jer_param[11][23]=20;  new3jer_param_er[11][23]=19.5753;
  new3jer_param[11][24]=-0.266033;  new3jer_param_er[11][24]=0.798582;
  new3jer_param[11][25]=0.0156233;  new3jer_param_er[11][25]=0.0365041;
  new3jer_param[11][26]=-0.000159787;  new3jer_param_er[11][26]=0.000179184;
  new3jer_param[11][27]=0.0482362;  new3jer_param_er[11][27]=0.109694;
  new3jer_param[11][28]=0.00176255;  new3jer_param_er[11][28]=0.000969561;
  new3jer_param[11][29]=-3.05466e-06;  new3jer_param_er[11][29]=1.98342e-06;
  new3jer_param[11][30]=70;  new3jer_param_er[11][30]=15.213;
  new3jer_param[12][0]=2.2418;  new3jer_param_er[12][0]=1.42044;
  new3jer_param[12][1]=-0.083243;  new3jer_param_er[12][1]=0.0512155;
  new3jer_param[12][2]=0.000714626;  new3jer_param_er[12][2]=0.000439433;
  new3jer_param[12][3]=-0.158378;  new3jer_param_er[12][3]=0.0153529;
  new3jer_param[12][4]=0.000980118;  new3jer_param_er[12][4]=0.000119128;
  new3jer_param[12][5]=-1.53723e-06;  new3jer_param_er[12][5]=2.21422e-07;
  new3jer_param[12][6]=50;  new3jer_param_er[12][6]=15.9727;
  new3jer_param[12][7]=-10.8934;  new3jer_param_er[12][7]=10.3289;
  new3jer_param[12][8]=0.265435;  new3jer_param_er[12][8]=0.206305;
  new3jer_param[12][9]=6.40372;  new3jer_param_er[12][9]=0.702999;
  new3jer_param[12][10]=0.0468618;  new3jer_param_er[12][10]=0.00304968;
  new3jer_param[12][11]=20;  new3jer_param_er[12][11]=16.5398;
  new3jer_param[12][12]=-3.86112;  new3jer_param_er[12][12]=1.53192;
  new3jer_param[12][13]=0.137826;  new3jer_param_er[12][13]=0.0581745;
  new3jer_param[12][14]=-0.00132888;  new3jer_param_er[12][14]=0.000563287;
  new3jer_param[12][15]=-0.112904;  new3jer_param_er[12][15]=0.0219059;
  new3jer_param[12][16]=0.000155259;  new3jer_param_er[12][16]=0.000168814;
  new3jer_param[12][17]=1.781e-07;  new3jer_param_er[12][17]=3.1241e-07;
  new3jer_param[12][18]=50;  new3jer_param_er[12][18]=16.07;
  new3jer_param[12][19]=7.5103;  new3jer_param_er[12][19]=7.06497;
  new3jer_param[12][20]=-0.136865;  new3jer_param_er[12][20]=0.187479;
  new3jer_param[12][21]=4.5667;  new3jer_param_er[12][21]=0.796348;
  new3jer_param[12][22]=0.0166231;  new3jer_param_er[12][22]=0.00338051;
  new3jer_param[12][23]=40;  new3jer_param_er[12][23]=19.9831;
  new3jer_param[12][24]=1.48987;  new3jer_param_er[12][24]=3.79174;
  new3jer_param[12][25]=-0.0374031;  new3jer_param_er[12][25]=0.144307;
  new3jer_param[12][26]=0.000233242;  new3jer_param_er[12][26]=0.00140381;
  new3jer_param[12][27]=-0.0442138;  new3jer_param_er[12][27]=0.155284;
  new3jer_param[12][28]=0.00262747;  new3jer_param_er[12][28]=0.00124186;
  new3jer_param[12][29]=-5.34113e-06;  new3jer_param_er[12][29]=2.39867e-06;
  new3jer_param[12][30]=50;  new3jer_param_er[12][30]=19.049;
  new3jer_param[13][0]=1.38325;  new3jer_param_er[13][0]=7.79836;
  new3jer_param[13][1]=-0.0745589;  new3jer_param_er[13][1]=0.247215;
  new3jer_param[13][2]=0.000891211;  new3jer_param_er[13][2]=0.00196207;
  new3jer_param[13][3]=-0.195559;  new3jer_param_er[13][3]=0.0524396;
  new3jer_param[13][4]=0.00110466;  new3jer_param_er[13][4]=0.000349004;
  new3jer_param[13][5]=-1.58679e-06;  new3jer_param_er[13][5]=5.6821e-07;
  new3jer_param[13][6]=50;  new3jer_param_er[13][6]=11.26;
  new3jer_param[13][7]=-45.8246;  new3jer_param_er[13][7]=52.7681;
  new3jer_param[13][8]=1.0905;  new3jer_param_er[13][8]=1.05066;
  new3jer_param[13][9]=-1.19985;  new3jer_param_er[13][9]=2.763;
  new3jer_param[13][10]=0.0982496;  new3jer_param_er[13][10]=0.0105469;
  new3jer_param[13][11]=40;  new3jer_param_er[13][11]=11.0014;
  new3jer_param[13][12]=0.564383;  new3jer_param_er[13][12]=3.15506;
  new3jer_param[13][13]=-0.024468;  new3jer_param_er[13][13]=0.0311928;
  new3jer_param[13][14]=0.000200516;  new3jer_param_er[13][14]=0.000802675;
  new3jer_param[13][15]=-0.221741;  new3jer_param_er[13][15]=0.0274687;
  new3jer_param[13][16]=0.000845349;  new3jer_param_er[13][16]=0.000183962;
  new3jer_param[13][17]=-1.15123e-06;  new3jer_param_er[13][17]=3.01128e-07;
  new3jer_param[13][18]=70;  new3jer_param_er[13][18]=19.738;
  new3jer_param[13][19]=-3.81075;  new3jer_param_er[13][19]=18.2324;
  new3jer_param[13][20]=0.115733;  new3jer_param_er[13][20]=0.347851;
  new3jer_param[13][21]=4.1508;  new3jer_param_er[13][21]=1.18806;
  new3jer_param[13][22]=0.028934;  new3jer_param_er[13][22]=0.00448828;
  new3jer_param[13][23]=40;  new3jer_param_er[13][23]=12.0145;
  new3jer_param[13][24]=-0.0561819;  new3jer_param_er[13][24]=2.53112;
  new3jer_param[13][25]=0.000710781;  new3jer_param_er[13][25]=0.0787036;
  new3jer_param[13][26]=4.00942e-05;  new3jer_param_er[13][26]=0.00275037;
  new3jer_param[13][27]=-0.0436665;  new3jer_param_er[13][27]=0.115274;
  new3jer_param[13][28]=0.001184;  new3jer_param_er[13][28]=0.000907485;
  new3jer_param[13][29]=-2.21109e-06;  new3jer_param_er[13][29]=1.44027e-06;
  new3jer_param[13][30]=70;  new3jer_param_er[13][30]=15.366;
  new3jer_param[14][0]=0;  new3jer_param_er[14][0]=0;
  new3jer_param[14][1]=0;  new3jer_param_er[14][1]=0;
  new3jer_param[14][2]=0;  new3jer_param_er[14][2]=0;
  new3jer_param[14][3]=0;  new3jer_param_er[14][3]=0;
  new3jer_param[14][4]=0;  new3jer_param_er[14][4]=0;
  new3jer_param[14][5]=0;  new3jer_param_er[14][5]=0;
  new3jer_param[14][6]=0;  new3jer_param_er[14][6]=0;
  new3jer_param[14][7]=0;  new3jer_param_er[14][7]=0;
  new3jer_param[14][8]=0;  new3jer_param_er[14][8]=0;
  new3jer_param[14][9]=0;  new3jer_param_er[14][9]=0;
  new3jer_param[14][10]=0;  new3jer_param_er[14][10]=0;
  new3jer_param[14][11]=0;  new3jer_param_er[14][11]=0;
  new3jer_param[14][12]=-581.332;  new3jer_param_er[14][12]=4.59018;
  new3jer_param[14][13]=10.9849;  new3jer_param_er[14][13]=0.0506322;
  new3jer_param[14][14]=-0.0520908;  new3jer_param_er[14][14]=0.000426188;
  new3jer_param[14][15]=-0.120105;  new3jer_param_er[14][15]=0.0264841;
  new3jer_param[14][16]=-6.1382e-05;  new3jer_param_er[14][16]=0.000142561;
  new3jer_param[14][17]=1.77937e-07;  new3jer_param_er[14][17]=1.92265e-07;
  new3jer_param[14][18]=59.1772;  new3jer_param_er[14][18]=2.19426;
  new3jer_param[14][19]=120.977;  new3jer_param_er[14][19]=240.795;
  new3jer_param[14][20]=-1.36305;  new3jer_param_er[14][20]=2.60413;
  new3jer_param[14][21]=8.27713;  new3jer_param_er[14][21]=1.78637;
  new3jer_param[14][22]=0.0251648;  new3jer_param_er[14][22]=0.00458752;
  new3jer_param[14][23]=20;  new3jer_param_er[14][23]=14.066;
  new3jer_param[14][24]=0;  new3jer_param_er[14][24]=0;
  new3jer_param[14][25]=0;  new3jer_param_er[14][25]=0;
  new3jer_param[14][26]=0;  new3jer_param_er[14][26]=0;
  new3jer_param[14][27]=0;  new3jer_param_er[14][27]=0;
  new3jer_param[14][28]=0;  new3jer_param_er[14][28]=0;
  new3jer_param[14][29]=0;  new3jer_param_er[14][29]=0;
  new3jer_param[14][30]=0;  new3jer_param_er[14][30]=0;

	//========================== NEW FIT FUNCTION PARAMETERS ========= 07-30-2009
/*	new3jer_param[0][0]=-0.107537;  new3jer_param_er[0][0]=0.0566252;
	new3jer_param[0][1]=0.0101998;  new3jer_param_er[0][1]=0.00267699;
	new3jer_param[0][2]=-9.99535e-05;  new3jer_param_er[0][2]=3.24013e-05;
	new3jer_param[0][3]=0.116372;  new3jer_param_er[0][3]=0.00743799;
	new3jer_param[0][4]=-0.000824633;  new3jer_param_er[0][4]=5.77058e-05;
	new3jer_param[0][5]=1.23834e-06;  new3jer_param_er[0][5]=1.02365e-07;
	new3jer_param[0][6]=50;  new3jer_param_er[0][6]=1.57999;
	new3jer_param[0][7]=-7089.48;  new3jer_param_er[0][7]=1734.39;
	new3jer_param[0][8]=-234.953;  new3jer_param_er[0][8]=59.5457;
	new3jer_param[0][9]=8.05359;  new3jer_param_er[0][9]=0.257165;
	new3jer_param[0][10]=0.0337035;  new3jer_param_er[0][10]=0.00125313;
	new3jer_param[0][11]=-174.28;  new3jer_param_er[0][11]=5.0391;
	new3jer_param[0][12]=-0.0570722;  new3jer_param_er[0][12]=0.0169371;
	new3jer_param[0][13]=0.00255033;  new3jer_param_er[0][13]=0.000750679;
	new3jer_param[0][14]=-1.95879e-05;  new3jer_param_er[0][14]=1.37464e-05;
	new3jer_param[0][15]=0.0518506;  new3jer_param_er[0][15]=0.0110479;
	new3jer_param[0][16]=-0.00015948;  new3jer_param_er[0][16]=0.000100111;
	new3jer_param[0][17]=6.38365e-07;  new3jer_param_er[0][17]=2.069e-07;
	new3jer_param[0][18]=70;  new3jer_param_er[0][18]=16.8986;
	new3jer_param[0][19]=25335.1;  new3jer_param_er[0][19]=41852;
	new3jer_param[0][20]=-851.295;  new3jer_param_er[0][20]=1415.52;
	new3jer_param[0][21]=0.776435;  new3jer_param_er[0][21]=0.247753;
	new3jer_param[0][22]=0.087487;  new3jer_param_er[0][22]=0.00160073;
	new3jer_param[0][23]=-259.129;  new3jer_param_er[0][23]=152.486;
	new3jer_param[0][24]=-0.107059;  new3jer_param_er[0][24]=0.0477939;
	new3jer_param[0][25]=0.00973563;  new3jer_param_er[0][25]=0.00298643;
	new3jer_param[0][26]=-9.48027e-05;  new3jer_param_er[0][26]=4.79747e-05;
	new3jer_param[0][27]=0.165775;  new3jer_param_er[0][27]=0.028031;
	new3jer_param[0][28]=-0.000941991;  new3jer_param_er[0][28]=0.000259297;
	new3jer_param[0][29]=2.93026e-06;  new3jer_param_er[0][29]=5.5319e-07;
	new3jer_param[0][30]=50;  new3jer_param_er[0][30]=12.3864;
	new3jer_param[1][0]=0.0884825;  new3jer_param_er[1][0]=0.0165723;
	new3jer_param[1][1]=-0.000396748;  new3jer_param_er[1][1]=0.000668837;
	new3jer_param[1][2]=9.00949e-06;  new3jer_param_er[1][2]=6.82975e-06;
	new3jer_param[1][3]=0.0721112;  new3jer_param_er[1][3]=0.00530366;
	new3jer_param[1][4]=-2.41626e-05;  new3jer_param_er[1][4]=4.58352e-05;
	new3jer_param[1][5]=-3.27678e-07;  new3jer_param_er[1][5]=9.01269e-08;
	new3jer_param[1][6]=70;  new3jer_param_er[1][6]=2.42223;
	new3jer_param[1][7]=-802.717;  new3jer_param_er[1][7]=2835.19;
	new3jer_param[1][8]=-47.1343;  new3jer_param_er[1][8]=166.355;
	new3jer_param[1][9]=4.62173;  new3jer_param_er[1][9]=0.166643;
	new3jer_param[1][10]=0.0568763;  new3jer_param_er[1][10]=0.000983391;
	new3jer_param[1][11]=-147.733;  new3jer_param_er[1][11]=103.925;
	new3jer_param[1][12]=-0.136141;  new3jer_param_er[1][12]=0.0147021;
	new3jer_param[1][13]=0.0046606;  new3jer_param_er[1][13]=0.000921656;
	new3jer_param[1][14]=-4.67532e-05;  new3jer_param_er[1][14]=1.28932e-05;
	new3jer_param[1][15]=-0.0192809;  new3jer_param_er[1][15]=0.00423984;
	new3jer_param[1][16]=0.000115645;  new3jer_param_er[1][16]=3.40824e-05;
	new3jer_param[1][17]=-1.86755e-07;  new3jer_param_er[1][17]=6.37078e-08;
	new3jer_param[1][18]=60.2651;  new3jer_param_er[1][18]=7.51143;
	new3jer_param[1][19]=1.04425;  new3jer_param_er[1][19]=0.0314378;
	new3jer_param[1][20]=0.0470521;  new3jer_param_er[1][20]=0.000308997;
	new3jer_param[1][21]=6.71153e+06;  new3jer_param_er[1][21]=9.61129e+06;
	new3jer_param[1][22]=-15365;  new3jer_param_er[1][22]=22045.2;
	new3jer_param[1][23]=743.232;  new3jer_param_er[1][23]=35.6213;
	new3jer_param[1][24]=-0.0442165;  new3jer_param_er[1][24]=0.0385695;
	new3jer_param[1][25]=0.0106688;  new3jer_param_er[1][25]=0.00305846;
	new3jer_param[1][26]=-0.000127049;  new3jer_param_er[1][26]=5.68601e-05;
	new3jer_param[1][27]=0.387361;  new3jer_param_er[1][27]=0.0194439;
	new3jer_param[1][28]=-0.00172015;  new3jer_param_er[1][28]=0.000150375;
	new3jer_param[1][29]=2.39865e-06;  new3jer_param_er[1][29]=2.89995e-07;
	new3jer_param[1][30]=50;  new3jer_param_er[1][30]=15.8042;
	new3jer_param[2][0]=-0.0558391;  new3jer_param_er[2][0]=0.0286116;
	new3jer_param[2][1]=0.0047738;  new3jer_param_er[2][1]=0.00145382;
	new3jer_param[2][2]=-4.14672e-05;  new3jer_param_er[2][2]=1.83055e-05;
	new3jer_param[2][3]=0.0985458;  new3jer_param_er[2][3]=0.00479701;
	new3jer_param[2][4]=-0.000310602;  new3jer_param_er[2][4]=3.95796e-05;
	new3jer_param[2][5]=1.84186e-07;  new3jer_param_er[2][5]=7.44813e-08;
	new3jer_param[2][6]=50;  new3jer_param_er[2][6]=19.7341;
	new3jer_param[2][7]=-3.8718e+06;  new3jer_param_er[2][7]=1.03329e+06;
	new3jer_param[2][8]=-255644;  new3jer_param_er[2][8]=58030.3;
	new3jer_param[2][9]=5.20742;  new3jer_param_er[2][9]=0.168824;
	new3jer_param[2][10]=0.052684;  new3jer_param_er[2][10]=0.000985285;
	new3jer_param[2][11]=-354.211;  new3jer_param_er[2][11]=5.1901;
	new3jer_param[2][12]=-0.137934;  new3jer_param_er[2][12]=0.025341;
	new3jer_param[2][13]=0.00307188;  new3jer_param_er[2][13]=0.00140956;
	new3jer_param[2][14]=-2.02804e-05;  new3jer_param_er[2][14]=1.73514e-05;
	new3jer_param[2][15]=-0.0246058;  new3jer_param_er[2][15]=0.00573657;
	new3jer_param[2][16]=0.000142928;  new3jer_param_er[2][16]=4.31133e-05;
	new3jer_param[2][17]=-2.64319e-07;  new3jer_param_er[2][17]=7.49754e-08;
	new3jer_param[2][18]=50;  new3jer_param_er[2][18]=17.9769;
	new3jer_param[2][19]=1.04798;  new3jer_param_er[2][19]=0.0359102;
	new3jer_param[2][20]=0.0494956;  new3jer_param_er[2][20]=0.000334066;
	new3jer_param[2][21]=-6.42951e+14;  new3jer_param_er[2][21]=3.00248e+14;
	new3jer_param[2][22]=1.31988e+12;  new3jer_param_er[2][22]=6.44048e+11;
	new3jer_param[2][23]=1222.92;  new3jer_param_er[2][23]=0.015738;
	new3jer_param[2][24]=-0.0605025;  new3jer_param_er[2][24]=0.0736183;
	new3jer_param[2][25]=0.0119383;  new3jer_param_er[2][25]=0.0045531;
	new3jer_param[2][26]=-0.00012213;  new3jer_param_er[2][26]=6.53851e-05;
	new3jer_param[2][27]=0.342351;  new3jer_param_er[2][27]=0.0173133;
	new3jer_param[2][28]=-0.00135649;  new3jer_param_er[2][28]=0.000133935;
	new3jer_param[2][29]=1.76527e-06;  new3jer_param_er[2][29]=2.42139e-07;
	new3jer_param[2][30]=50;  new3jer_param_er[2][30]=15.9222;
	new3jer_param[3][0]=0.220049;  new3jer_param_er[3][0]=0.0229101;
	new3jer_param[3][1]=-0.00716605;  new3jer_param_er[3][1]=0.00121414;
	new3jer_param[3][2]=9.03346e-05;  new3jer_param_er[3][2]=1.63624e-05;
	new3jer_param[3][3]=0.0913929;  new3jer_param_er[3][3]=0.0045358;
	new3jer_param[3][4]=-0.000150115;  new3jer_param_er[3][4]=3.99071e-05;
	new3jer_param[3][5]=-7.54507e-08;  new3jer_param_er[3][5]=8.07302e-08;
	new3jer_param[3][6]=50;  new3jer_param_er[3][6]=0.912723;
	new3jer_param[3][7]=2.29098;  new3jer_param_er[3][7]=0.0871666;
	new3jer_param[3][8]=0.0817857;  new3jer_param_er[3][8]=0.00334469;
	new3jer_param[3][9]=5.22692;  new3jer_param_er[3][9]=0.465349;
	new3jer_param[3][10]=0.0574838;  new3jer_param_er[3][10]=0.00206276;
	new3jer_param[3][11]=112.199;  new3jer_param_er[3][11]=23.9142;
	new3jer_param[3][12]=-0.184246;  new3jer_param_er[3][12]=0.00930048;
	new3jer_param[3][13]=0.00538507;  new3jer_param_er[3][13]=0.000441317;
	new3jer_param[3][14]=-4.49893e-05;  new3jer_param_er[3][14]=5.23713e-06;
	new3jer_param[3][15]=-0.000272302;  new3jer_param_er[3][15]=0.00402049;
	new3jer_param[3][16]=6.88301e-05;  new3jer_param_er[3][16]=3.04064e-05;
	new3jer_param[3][17]=-1.80491e-07;  new3jer_param_er[3][17]=5.44626e-08;
	new3jer_param[3][18]=70;  new3jer_param_er[3][18]=12.891;
	new3jer_param[3][19]=0.959735;  new3jer_param_er[3][19]=0.0376819;
	new3jer_param[3][20]=0.0558552;  new3jer_param_er[3][20]=0.000844841;
	new3jer_param[3][21]=2.75095;  new3jer_param_er[3][21]=0.371093;
	new3jer_param[3][22]=0.044675;  new3jer_param_er[3][22]=0.00142469;
	new3jer_param[3][23]=165.644;  new3jer_param_er[3][23]=27.3238;
	new3jer_param[3][24]=0.0703947;  new3jer_param_er[3][24]=0.0316273;
	new3jer_param[3][25]=0.00104987;  new3jer_param_er[3][25]=0.002421;
	new3jer_param[3][26]=3.07433e-05;  new3jer_param_er[3][26]=4.11166e-05;
	new3jer_param[3][27]=0.275658;  new3jer_param_er[3][27]=0.0150253;
	new3jer_param[3][28]=-0.000975134;  new3jer_param_er[3][28]=0.000116931;
	new3jer_param[3][29]=1.15296e-06;  new3jer_param_er[3][29]=2.15012e-07;
	new3jer_param[3][30]=50;  new3jer_param_er[3][30]=3.06942;
	new3jer_param[4][0]=0.130408;  new3jer_param_er[4][0]=0.0222449;
	new3jer_param[4][1]=-0.00340213;  new3jer_param_er[4][1]=0.00129552;
	new3jer_param[4][2]=5.74043e-05;  new3jer_param_er[4][2]=1.81035e-05;
	new3jer_param[4][3]=0.122659;  new3jer_param_er[4][3]=0.00503038;
	new3jer_param[4][4]=-0.000292842;  new3jer_param_er[4][4]=4.61727e-05;
	new3jer_param[4][5]=3.60169e-07;  new3jer_param_er[4][5]=9.6853e-08;
	new3jer_param[4][6]=50;  new3jer_param_er[4][6]=12.3584;
	new3jer_param[4][7]=-1100.64;  new3jer_param_er[4][7]=340.335;
	new3jer_param[4][8]=-28.0531;  new3jer_param_er[4][8]=11.3309;
	new3jer_param[4][9]=4.54038;  new3jer_param_er[4][9]=0.201987;
	new3jer_param[4][10]=0.0782949;  new3jer_param_er[4][10]=0.00121771;
	new3jer_param[4][11]=-146.2;  new3jer_param_er[4][11]=7.608;
	new3jer_param[4][12]=-0.165712;  new3jer_param_er[4][12]=0.0240303;
	new3jer_param[4][13]=0.00507425;  new3jer_param_er[4][13]=0.00153825;
	new3jer_param[4][14]=-4.33962e-05;  new3jer_param_er[4][14]=2.58543e-05;
	new3jer_param[4][15]=0.00596175;  new3jer_param_er[4][15]=0.00760369;
	new3jer_param[4][16]=9.50617e-05;  new3jer_param_er[4][16]=5.59348e-05;
	new3jer_param[4][17]=-2.61393e-07;  new3jer_param_er[4][17]=9.77872e-08;
	new3jer_param[4][18]=70;  new3jer_param_er[4][18]=19.9803;
	new3jer_param[4][19]=1.06979;  new3jer_param_er[4][19]=0.039438;
	new3jer_param[4][20]=0.0641176;  new3jer_param_er[4][20]=0.000568913;
	new3jer_param[4][21]=5.83718;  new3jer_param_er[4][21]=3.0095;
	new3jer_param[4][22]=0.0434838;  new3jer_param_er[4][22]=0.00753371;
	new3jer_param[4][23]=314.545;  new3jer_param_er[4][23]=22.8683;
	new3jer_param[4][24]=0.143285;  new3jer_param_er[4][24]=0.0409534;
	new3jer_param[4][25]=-0.0037696;  new3jer_param_er[4][25]=0.00300876;
	new3jer_param[4][26]=7.85254e-05;  new3jer_param_er[4][26]=4.86225e-05;
	new3jer_param[4][27]=0.321488;  new3jer_param_er[4][27]=0.0186633;
	new3jer_param[4][28]=-0.000971977;  new3jer_param_er[4][28]=0.000144944;
	new3jer_param[4][29]=9.86957e-07;  new3jer_param_er[4][29]=2.66791e-07;
	new3jer_param[4][30]=50;  new3jer_param_er[4][30]=2.62261;
	new3jer_param[5][0]=0.0735202;  new3jer_param_er[5][0]=0.0219619;
	new3jer_param[5][1]=0.00018131;  new3jer_param_er[5][1]=0.000842598;
	new3jer_param[5][2]=-4.52704e-06;  new3jer_param_er[5][2]=8.93316e-06;
	new3jer_param[5][3]=0.0908279;  new3jer_param_er[5][3]=0.00405447;
	new3jer_param[5][4]=-0.000112972;  new3jer_param_er[5][4]=3.28994e-05;
	new3jer_param[5][5]=6.59606e-08;  new3jer_param_er[5][5]=6.13004e-08;
	new3jer_param[5][6]=70;  new3jer_param_er[5][6]=14.2521;
	new3jer_param[5][7]=-1.90016;  new3jer_param_er[5][7]=1.41421;
	new3jer_param[5][8]=0.554732;  new3jer_param_er[5][8]=1.41421;
	new3jer_param[5][9]=2.9286;  new3jer_param_er[5][9]=0.0532931;
	new3jer_param[5][10]=0.0861254;  new3jer_param_er[5][10]=0.000518283;
	new3jer_param[5][11]=-702.403;  new3jer_param_er[5][11]=1.41421;
	new3jer_param[5][12]=-0.301711;  new3jer_param_er[5][12]=0.0292409;
	new3jer_param[5][13]=0.0128286;  new3jer_param_er[5][13]=0.00169054;
	new3jer_param[5][14]=-0.000181522;  new3jer_param_er[5][14]=2.39601e-05;
	new3jer_param[5][15]=-0.0092045;  new3jer_param_er[5][15]=0.0053185;
	new3jer_param[5][16]=6.96301e-05;  new3jer_param_er[5][16]=4.17049e-05;
	new3jer_param[5][17]=-1.83794e-07;  new3jer_param_er[5][17]=7.56229e-08;
	new3jer_param[5][18]=50;  new3jer_param_er[5][18]=1.90508;
	new3jer_param[5][19]=1.5709;  new3jer_param_er[5][19]=0.130092;
	new3jer_param[5][20]=0.0371641;  new3jer_param_er[5][20]=0.00534915;
	new3jer_param[5][21]=3.30615;  new3jer_param_er[5][21]=0.446914;
	new3jer_param[5][22]=0.0478183;  new3jer_param_er[5][22]=0.00163622;
	new3jer_param[5][23]=92.0985;  new3jer_param_er[5][23]=12.0986;
	new3jer_param[5][24]=0.178338;  new3jer_param_er[5][24]=0.0640275;
	new3jer_param[5][25]=-0.00319376;  new3jer_param_er[5][25]=0.0038986;
	new3jer_param[5][26]=5.59603e-05;  new3jer_param_er[5][26]=6.59878e-05;
	new3jer_param[5][27]=0.508882;  new3jer_param_er[5][27]=0.0433818;
	new3jer_param[5][28]=-0.000933784;  new3jer_param_er[5][28]=0.000312936;
	new3jer_param[5][29]=4.68644e-07;  new3jer_param_er[5][29]=5.33297e-07;
	new3jer_param[5][30]=68.3458;  new3jer_param_er[5][30]=11.8792;
	new3jer_param[6][0]=0.0347005;  new3jer_param_er[6][0]=0.0184864;
	new3jer_param[6][1]=0.000937425;  new3jer_param_er[6][1]=0.000830865;
	new3jer_param[6][2]=-1.30709e-05;  new3jer_param_er[6][2]=9.21559e-06;
	new3jer_param[6][3]=0.0875522;  new3jer_param_er[6][3]=0.00888453;
	new3jer_param[6][4]=-8.55151e-05;  new3jer_param_er[6][4]=8.18485e-05;
	new3jer_param[6][5]=-4.6158e-07;  new3jer_param_er[6][5]=1.52946e-07;
	new3jer_param[6][6]=70;  new3jer_param_er[6][6]=4.03686;
	new3jer_param[6][7]=-3.86139e+06;  new3jer_param_er[6][7]=1.67035e+06;
	new3jer_param[6][8]=-1.7334e+06;  new3jer_param_er[6][8]=259031;
	new3jer_param[6][9]=7.2893;  new3jer_param_er[6][9]=0.301931;
	new3jer_param[6][10]=0.0535885;  new3jer_param_er[6][10]=0.00203214;
	new3jer_param[6][11]=-379.358;  new3jer_param_er[6][11]=3.67917;
	new3jer_param[6][12]=-0.344393;  new3jer_param_er[6][12]=0.0237024;
	new3jer_param[6][13]=0.0131182;  new3jer_param_er[6][13]=0.00148083;
	new3jer_param[6][14]=-0.000180064;  new3jer_param_er[6][14]=2.18422e-05;
	new3jer_param[6][15]=-0.0105614;  new3jer_param_er[6][15]=0.00699589;
	new3jer_param[6][16]=-2.36259e-05;  new3jer_param_er[6][16]=7.01264e-05;
	new3jer_param[6][17]=2.83426e-07;  new3jer_param_er[6][17]=1.56341e-07;
	new3jer_param[6][18]=50;  new3jer_param_er[6][18]=1.79359;
	new3jer_param[6][19]=1.22181;  new3jer_param_er[6][19]=0.705988;
	new3jer_param[6][20]=0.0518528;  new3jer_param_er[6][20]=0.0275788;
	new3jer_param[6][21]=0.0433576;  new3jer_param_er[6][21]=0.35759;
	new3jer_param[6][22]=0.0686263;  new3jer_param_er[6][22]=0.00177367;
	new3jer_param[6][23]=49.1294;  new3jer_param_er[6][23]=38.8008;
	new3jer_param[6][24]=0.325285;  new3jer_param_er[6][24]=0.0729539;
	new3jer_param[6][25]=-0.0133924;  new3jer_param_er[6][25]=0.00454512;
	new3jer_param[6][26]=0.000190133;  new3jer_param_er[6][26]=6.429e-05;
	new3jer_param[6][27]=0.329918;  new3jer_param_er[6][27]=0.0247353;
	new3jer_param[6][28]=-0.00194938;  new3jer_param_er[6][28]=0.00023475;
	new3jer_param[6][29]=3.39452e-06;  new3jer_param_er[6][29]=5.05375e-07;
	new3jer_param[6][30]=50;  new3jer_param_er[6][30]=1.63864;
	new3jer_param[7][0]=-0.0363087;  new3jer_param_er[7][0]=0.0159478;
	new3jer_param[7][1]=0.00256729;  new3jer_param_er[7][1]=0.000771292;
	new3jer_param[7][2]=-2.37768e-05;  new3jer_param_er[7][2]=8.48443e-06;
	new3jer_param[7][3]=0.0582484;  new3jer_param_er[7][3]=0.00712301;
	new3jer_param[7][4]=-0.000189063;  new3jer_param_er[7][4]=6.41661e-05;
	new3jer_param[7][5]=1.15837e-08;  new3jer_param_er[7][5]=1.27971e-07;
	new3jer_param[7][6]=70;  new3jer_param_er[7][6]=1.51805;
	new3jer_param[7][7]=-911.67;  new3jer_param_er[7][7]=201.335;
	new3jer_param[7][8]=-43.8908;  new3jer_param_er[7][8]=10.7679;
	new3jer_param[7][9]=6.48176;  new3jer_param_er[7][9]=0.239973;
	new3jer_param[7][10]=0.0433151;  new3jer_param_er[7][10]=0.00157171;
	new3jer_param[7][11]=-136.428;  new3jer_param_er[7][11]=4.89411;
	new3jer_param[7][12]=-0.296756;  new3jer_param_er[7][12]=0.0150918;
	new3jer_param[7][13]=0.00867224;  new3jer_param_er[7][13]=0.000771177;
	new3jer_param[7][14]=-8.12718e-05;  new3jer_param_er[7][14]=8.72039e-06;
	new3jer_param[7][15]=-0.0690383;  new3jer_param_er[7][15]=0.00880247;
	new3jer_param[7][16]=0.000460047;  new3jer_param_er[7][16]=8.56913e-05;
	new3jer_param[7][17]=-8.53995e-07;  new3jer_param_er[7][17]=1.89517e-07;
	new3jer_param[7][18]=70;  new3jer_param_er[7][18]=0.583276;
	new3jer_param[7][19]=-4798.84;  new3jer_param_er[7][19]=1399.64;
	new3jer_param[7][20]=567.747;  new3jer_param_er[7][20]=152.571;
	new3jer_param[7][21]=-0.0761321;  new3jer_param_er[7][21]=0.260785;
	new3jer_param[7][22]=0.0476726;  new3jer_param_er[7][22]=0.00176623;
	new3jer_param[7][23]=-194.097;  new3jer_param_er[7][23]=7.08988;
	new3jer_param[7][24]=0.345809;  new3jer_param_er[7][24]=0.059104;
	new3jer_param[7][25]=-0.00997728;  new3jer_param_er[7][25]=0.00338212;
	new3jer_param[7][26]=0.000137981;  new3jer_param_er[7][26]=4.14079e-05;
	new3jer_param[7][27]=0.214122;  new3jer_param_er[7][27]=0.0613989;
	new3jer_param[7][28]=0.000120621;  new3jer_param_er[7][28]=0.00062674;
	new3jer_param[7][29]=-9.90705e-07;  new3jer_param_er[7][29]=1.43959e-06;
	new3jer_param[7][30]=70;  new3jer_param_er[7][30]=1.79527;
	new3jer_param[8][0]=-0.045249;  new3jer_param_er[8][0]=0.0240541;
	new3jer_param[8][1]=0.00248962;  new3jer_param_er[8][1]=0.00102856;
	new3jer_param[8][2]=-3.55084e-05;  new3jer_param_er[8][2]=1.00571e-05;
	new3jer_param[8][3]=0.0274564;  new3jer_param_er[8][3]=0.00381056;
	new3jer_param[8][4]=5.67099e-05;  new3jer_param_er[8][4]=3.09473e-05;
	new3jer_param[8][5]=-2.32842e-07;  new3jer_param_er[8][5]=5.68585e-08;
	new3jer_param[8][6]=70;  new3jer_param_er[8][6]=1.60276;
	new3jer_param[8][7]=947883;  new3jer_param_er[8][7]=4.32227e+06;
	new3jer_param[8][8]=-343241;  new3jer_param_er[8][8]=1.54554e+06;
	new3jer_param[8][9]=6.70206;  new3jer_param_er[8][9]=0.177496;
	new3jer_param[8][10]=0.0387028;  new3jer_param_er[8][10]=0.000965875;
	new3jer_param[8][11]=-343.303;  new3jer_param_er[8][11]=280.287;
	new3jer_param[8][12]=-0.224642;  new3jer_param_er[8][12]=0.0481272;
	new3jer_param[8][13]=0.00392954;  new3jer_param_er[8][13]=0.00252884;
	new3jer_param[8][14]=-2.57769e-05;  new3jer_param_er[8][14]=3.28976e-05;
	new3jer_param[8][15]=-0.10305;  new3jer_param_er[8][15]=0.00701373;
	new3jer_param[8][16]=0.000563536;  new3jer_param_er[8][16]=5.94671e-05;
	new3jer_param[8][17]=-8.82646e-07;  new3jer_param_er[8][17]=1.10424e-07;
	new3jer_param[8][18]=70;  new3jer_param_er[8][18]=17.5701;
	new3jer_param[8][19]=-123643;  new3jer_param_er[8][19]=121284;
	new3jer_param[8][20]=3546.45;  new3jer_param_er[8][20]=3568.02;
	new3jer_param[8][21]=2.41266;  new3jer_param_er[8][21]=0.248617;
	new3jer_param[8][22]=0.0185443;  new3jer_param_er[8][22]=0.00148184;
	new3jer_param[8][23]=-268.821;  new3jer_param_er[8][23]=152.679;
	new3jer_param[8][24]=0.422792;  new3jer_param_er[8][24]=0.166925;
	new3jer_param[8][25]=-0.013864;  new3jer_param_er[8][25]=0.0109505;
	new3jer_param[8][26]=0.000171851;  new3jer_param_er[8][26]=0.000191956;
	new3jer_param[8][27]=0.272589;  new3jer_param_er[8][27]=0.0977811;
	new3jer_param[8][28]=0.00155807;  new3jer_param_er[8][28]=0.000883326;
	new3jer_param[8][29]=-3.36749e-06;  new3jer_param_er[8][29]=1.69713e-06;
	new3jer_param[8][30]=70;  new3jer_param_er[8][30]=15.6354;
	new3jer_param[9][0]=-0.174942;  new3jer_param_er[9][0]=0.0410119;
	new3jer_param[9][1]=0.00849332;  new3jer_param_er[9][1]=0.0017247;
	new3jer_param[9][2]=-0.000106572;  new3jer_param_er[9][2]=1.66441e-05;
	new3jer_param[9][3]=0.0340075;  new3jer_param_er[9][3]=0.00655096;
	new3jer_param[9][4]=-3.63072e-05;  new3jer_param_er[9][4]=6.33251e-05;
	new3jer_param[9][5]=1.85558e-08;  new3jer_param_er[9][5]=1.46601e-07;
	new3jer_param[9][6]=70;  new3jer_param_er[9][6]=0.479253;
	new3jer_param[9][7]=-4.26486e+09;  new3jer_param_er[9][7]=1.41421;
	new3jer_param[9][8]=-3.13846e+12;  new3jer_param_er[9][8]=1.41421;
	new3jer_param[9][9]=6.02024;  new3jer_param_er[9][9]=0.186015;
	new3jer_param[9][10]=0.0466621;  new3jer_param_er[9][10]=0.00110913;
	new3jer_param[9][11]=-755.283;  new3jer_param_er[9][11]=2.45506;
	new3jer_param[9][12]=-0.337628;  new3jer_param_er[9][12]=0.0457211;
	new3jer_param[9][13]=0.00796001;  new3jer_param_er[9][13]=0.00208963;
	new3jer_param[9][14]=-6.18984e-05;  new3jer_param_er[9][14]=2.24107e-05;
	new3jer_param[9][15]=-0.108523;  new3jer_param_er[9][15]=0.0114774;
	new3jer_param[9][16]=0.000689633;  new3jer_param_er[9][16]=0.000114841;
	new3jer_param[9][17]=-1.32862e-06;  new3jer_param_er[9][17]=2.70826e-07;
	new3jer_param[9][18]=70;  new3jer_param_er[9][18]=4.01764;
	new3jer_param[9][19]=-52006;  new3jer_param_er[9][19]=39149.8;
	new3jer_param[9][20]=1318.69;  new3jer_param_er[9][20]=1014.82;
	new3jer_param[9][21]=2.56323;  new3jer_param_er[9][21]=0.332443;
	new3jer_param[9][22]=0.0185645;  new3jer_param_er[9][22]=0.00192597;
	new3jer_param[9][23]=-230.54;  new3jer_param_er[9][23]=129.566;
	new3jer_param[9][24]=0.838697;  new3jer_param_er[9][24]=0.195885;
	new3jer_param[9][25]=-0.0333986;  new3jer_param_er[9][25]=0.00970986;
	new3jer_param[9][26]=0.000361833;  new3jer_param_er[9][26]=0.000116387;
	new3jer_param[9][27]=0.267283;  new3jer_param_er[9][27]=0.118979;
	new3jer_param[9][28]=0.000616332;  new3jer_param_er[9][28]=0.00124207;
	new3jer_param[9][29]=-9.39781e-07;  new3jer_param_er[9][29]=2.97127e-06;
	new3jer_param[9][30]=70;  new3jer_param_er[9][30]=2.85094;
	new3jer_param[10][0]=-0.574168;  new3jer_param_er[10][0]=0.0854421;
	new3jer_param[10][1]=0.0228021;  new3jer_param_er[10][1]=0.00331966;
	new3jer_param[10][2]=-0.000230131;  new3jer_param_er[10][2]=3.14273e-05;
	new3jer_param[10][3]=0.0324899;  new3jer_param_er[10][3]=0.00590435;
	new3jer_param[10][4]=-5.99141e-05;  new3jer_param_er[10][4]=5.00397e-05;
	new3jer_param[10][5]=2.3651e-10;  new3jer_param_er[10][5]=1.00775e-07;
	new3jer_param[10][6]=70;  new3jer_param_er[10][6]=1.22064;
	new3jer_param[10][7]=-1.08159;  new3jer_param_er[10][7]=19.6405;
	new3jer_param[10][8]=0.284454;  new3jer_param_er[10][8]=11.4296;
	new3jer_param[10][9]=5.20702;  new3jer_param_er[10][9]=1.7609e-14;
	new3jer_param[10][10]=0.0534543;  new3jer_param_er[10][10]=1.08778e-16;
	new3jer_param[10][11]=-803.165;  new3jer_param_er[10][11]=1.06178e-10;
	new3jer_param[10][12]=-0.133786;  new3jer_param_er[10][12]=0.196081;
	new3jer_param[10][13]=0.0015014;  new3jer_param_er[10][13]=0.00912009;
	new3jer_param[10][14]=-3.44832e-05;  new3jer_param_er[10][14]=0.000104836;
	new3jer_param[10][15]=-0.0996914;  new3jer_param_er[10][15]=0.00913848;
	new3jer_param[10][16]=0.000469939;  new3jer_param_er[10][16]=8.31301e-05;
	new3jer_param[10][17]=-7.81904e-07;  new3jer_param_er[10][17]=1.74619e-07;
	new3jer_param[10][18]=50;  new3jer_param_er[10][18]=13.5254;
	new3jer_param[10][19]=-27620.7;  new3jer_param_er[10][19]=80517;
	new3jer_param[10][20]=863.932;  new3jer_param_er[10][20]=2513.35;
	new3jer_param[10][21]=2.69726;  new3jer_param_er[10][21]=0.186123;
	new3jer_param[10][22]=0.021745;  new3jer_param_er[10][22]=0.00115899;
	new3jer_param[10][23]=-237.954;  new3jer_param_er[10][23]=83.4023;
	new3jer_param[10][24]=-0.159945;  new3jer_param_er[10][24]=0.279156;
	new3jer_param[10][25]=0.0165115;  new3jer_param_er[10][25]=0.0125563;
	new3jer_param[10][26]=-0.000229064;  new3jer_param_er[10][26]=0.000143604;
	new3jer_param[10][27]=0.364856;  new3jer_param_er[10][27]=0.106014;
	new3jer_param[10][28]=-0.000667039;  new3jer_param_er[10][28]=0.00101174;
	new3jer_param[10][29]=1.66469e-06;  new3jer_param_er[10][29]=2.21821e-06;
	new3jer_param[10][30]=70;  new3jer_param_er[10][30]=12.5585;
	new3jer_param[11][0]=-0.90196;  new3jer_param_er[11][0]=0.515815;
	new3jer_param[11][1]=0.0369416;  new3jer_param_er[11][1]=0.0222684;
	new3jer_param[11][2]=-0.000395616;  new3jer_param_er[11][2]=0.000237496;
	new3jer_param[11][3]=-0.0477924;  new3jer_param_er[11][3]=0.00998877;
	new3jer_param[11][4]=0.000390147;  new3jer_param_er[11][4]=8.47609e-05;
	new3jer_param[11][5]=-6.39862e-07;  new3jer_param_er[11][5]=1.66749e-07;
	new3jer_param[11][6]=50;  new3jer_param_er[11][6]=16.3431;
	new3jer_param[11][7]=0.733495;  new3jer_param_er[11][7]=7.84629;
	new3jer_param[11][8]=0.190938;  new3jer_param_er[11][8]=3.97043;
	new3jer_param[11][9]=5.95296;  new3jer_param_er[11][9]=1.72789e-13;
	new3jer_param[11][10]=0.0491914;  new3jer_param_er[11][10]=9.34959e-16;
	new3jer_param[11][11]=-717.163;  new3jer_param_er[11][11]=3.59469e-09;
	new3jer_param[11][12]=-0.294175;  new3jer_param_er[11][12]=0.323618;
	new3jer_param[11][13]=0.00464302;  new3jer_param_er[11][13]=0.0142533;
	new3jer_param[11][14]=-2.68264e-05;  new3jer_param_er[11][14]=0.000152949;
	new3jer_param[11][15]=-0.155451;  new3jer_param_er[11][15]=0.0148115;
	new3jer_param[11][16]=0.000801954;  new3jer_param_er[11][16]=0.000125203;
	new3jer_param[11][17]=-1.30184e-06;  new3jer_param_er[11][17]=2.44338e-07;
	new3jer_param[11][18]=70;  new3jer_param_er[11][18]=19.0245;
	new3jer_param[11][19]=-55450.4;  new3jer_param_er[11][19]=78220.1;
	new3jer_param[11][20]=836.067;  new3jer_param_er[11][20]=1298.14;
	new3jer_param[11][21]=3.15804;  new3jer_param_er[11][21]=0.486865;
	new3jer_param[11][22]=0.0211629;  new3jer_param_er[11][22]=0.00239084;
	new3jer_param[11][23]=-225.611;  new3jer_param_er[11][23]=153.858;
	new3jer_param[11][24]=-0.266033;  new3jer_param_er[11][24]=0.798582;
	new3jer_param[11][25]=0.0156233;  new3jer_param_er[11][25]=0.0365041;
	new3jer_param[11][26]=-0.000159787;  new3jer_param_er[11][26]=0.000179184;
	new3jer_param[11][27]=0.0482362;  new3jer_param_er[11][27]=0.109694;
	new3jer_param[11][28]=0.00176255;  new3jer_param_er[11][28]=0.000969561;
	new3jer_param[11][29]=-3.05466e-06;  new3jer_param_er[11][29]=1.98342e-06;
	new3jer_param[11][30]=70;  new3jer_param_er[11][30]=15.213;
	new3jer_param[12][0]=2.2418;  new3jer_param_er[12][0]=1.42044;
	new3jer_param[12][1]=-0.083243;  new3jer_param_er[12][1]=0.0512155;
	new3jer_param[12][2]=0.000714626;  new3jer_param_er[12][2]=0.000439433;
	new3jer_param[12][3]=-0.158378;  new3jer_param_er[12][3]=0.0153529;
	new3jer_param[12][4]=0.000980118;  new3jer_param_er[12][4]=0.000119128;
	new3jer_param[12][5]=-1.53723e-06;  new3jer_param_er[12][5]=2.21422e-07;
	new3jer_param[12][6]=50;  new3jer_param_er[12][6]=15.9727;
	new3jer_param[12][7]=4.35732;  new3jer_param_er[12][7]=1.41421;
	new3jer_param[12][8]=-333223;  new3jer_param_er[12][8]=1.41421;
	new3jer_param[12][9]=4.88449;  new3jer_param_er[12][9]=0.240715;
	new3jer_param[12][10]=0.0535502;  new3jer_param_er[12][10]=0.00129543;
	new3jer_param[12][11]=-764.866;  new3jer_param_er[12][11]=1.536;
	new3jer_param[12][12]=-3.86112;  new3jer_param_er[12][12]=1.53192;
	new3jer_param[12][13]=0.137826;  new3jer_param_er[12][13]=0.0581745;
	new3jer_param[12][14]=-0.00132888;  new3jer_param_er[12][14]=0.000563287;
	new3jer_param[12][15]=-0.112904;  new3jer_param_er[12][15]=0.0219059;
	new3jer_param[12][16]=0.000155259;  new3jer_param_er[12][16]=0.000168814;
	new3jer_param[12][17]=1.781e-07;  new3jer_param_er[12][17]=3.1241e-07;
	new3jer_param[12][18]=50;  new3jer_param_er[12][18]=16.07;
	new3jer_param[12][19]=1.87515;  new3jer_param_er[12][19]=0.219822;
	new3jer_param[12][20]=0.0310785;  new3jer_param_er[12][20]=0.00184407;
	new3jer_param[12][21]=-23.8332;  new3jer_param_er[12][21]=21.1476;
	new3jer_param[12][22]=0.0885073;  new3jer_param_er[12][22]=0.0547578;
	new3jer_param[12][23]=327.565;  new3jer_param_er[12][23]=22.8698;
	new3jer_param[12][24]=1.48987;  new3jer_param_er[12][24]=3.79174;
	new3jer_param[12][25]=-0.0374031;  new3jer_param_er[12][25]=0.144307;
	new3jer_param[12][26]=0.000233242;  new3jer_param_er[12][26]=0.00140381;
	new3jer_param[12][27]=-0.0442138;  new3jer_param_er[12][27]=0.155284;
	new3jer_param[12][28]=0.00262747;  new3jer_param_er[12][28]=0.00124186;
	new3jer_param[12][29]=-5.34113e-06;  new3jer_param_er[12][29]=2.39867e-06;
	new3jer_param[12][30]=50;  new3jer_param_er[12][30]=19.049;
	new3jer_param[13][0]=1.38325;  new3jer_param_er[13][0]=7.79836;
	new3jer_param[13][1]=-0.0745589;  new3jer_param_er[13][1]=0.247215;
	new3jer_param[13][2]=0.000891211;  new3jer_param_er[13][2]=0.00196207;
	new3jer_param[13][3]=-0.195559;  new3jer_param_er[13][3]=0.0524396;
	new3jer_param[13][4]=0.00110466;  new3jer_param_er[13][4]=0.000349004;
	new3jer_param[13][5]=-1.58679e-06;  new3jer_param_er[13][5]=5.6821e-07;
	new3jer_param[13][6]=50;  new3jer_param_er[13][6]=11.26;
	new3jer_param[13][7]=-3.18129;  new3jer_param_er[13][7]=6.58953;
	new3jer_param[13][8]=0.216038;  new3jer_param_er[13][8]=0.117927;
	new3jer_param[13][9]=-3.1157;  new3jer_param_er[13][9]=3.6218;
	new3jer_param[13][10]=0.104868;  new3jer_param_er[13][10]=0.01315;
	new3jer_param[13][11]=122.534;  new3jer_param_er[13][11]=44.2217;
	new3jer_param[13][12]=0.564383;  new3jer_param_er[13][12]=3.15506;
	new3jer_param[13][13]=-0.024468;  new3jer_param_er[13][13]=0.0311928;
	new3jer_param[13][14]=0.000200516;  new3jer_param_er[13][14]=0.000802675;
	new3jer_param[13][15]=-0.221741;  new3jer_param_er[13][15]=0.0274687;
	new3jer_param[13][16]=0.000845349;  new3jer_param_er[13][16]=0.000183962;
	new3jer_param[13][17]=-1.15123e-06;  new3jer_param_er[13][17]=3.01128e-07;
	new3jer_param[13][18]=70;  new3jer_param_er[13][18]=19.738;
	new3jer_param[13][19]=0.875323;  new3jer_param_er[13][19]=307.001;
	new3jer_param[13][20]=0.123781;  new3jer_param_er[13][20]=26.7399;
	new3jer_param[13][21]=3.1944;  new3jer_param_er[13][21]=8.66561e-13;
	new3jer_param[13][22]=0.0325387;  new3jer_param_er[13][22]=4.72083e-15;
	new3jer_param[13][23]=-732.744;  new3jer_param_er[13][23]=2.01237e-08;
	new3jer_param[13][24]=-0.0561819;  new3jer_param_er[13][24]=2.53112;
	new3jer_param[13][25]=0.000710781;  new3jer_param_er[13][25]=0.0787036;
	new3jer_param[13][26]=4.00942e-05;  new3jer_param_er[13][26]=0.00275037;
	new3jer_param[13][27]=-0.0436665;  new3jer_param_er[13][27]=0.115274;
	new3jer_param[13][28]=0.001184;  new3jer_param_er[13][28]=0.000907485;
	new3jer_param[13][29]=-2.21109e-06;  new3jer_param_er[13][29]=1.44027e-06;
	new3jer_param[13][30]=70;  new3jer_param_er[13][30]=15.366;
	new3jer_param[14][0]=0;  new3jer_param_er[14][0]=0;
	new3jer_param[14][1]=0;  new3jer_param_er[14][1]=0;
	new3jer_param[14][2]=0;  new3jer_param_er[14][2]=0;
	new3jer_param[14][3]=0;  new3jer_param_er[14][3]=0;
	new3jer_param[14][4]=0;  new3jer_param_er[14][4]=0;
	new3jer_param[14][5]=0;  new3jer_param_er[14][5]=0;
	new3jer_param[14][6]=0;  new3jer_param_er[14][6]=0;
	new3jer_param[14][7]=0;  new3jer_param_er[14][7]=0;
	new3jer_param[14][8]=0;  new3jer_param_er[14][8]=0;
	new3jer_param[14][9]=0;  new3jer_param_er[14][9]=0;
	new3jer_param[14][10]=0;  new3jer_param_er[14][10]=0;
	new3jer_param[14][11]=0;  new3jer_param_er[14][11]=0;
	new3jer_param[14][12]=-581.332;  new3jer_param_er[14][12]=4.59018;
	new3jer_param[14][13]=10.9849;  new3jer_param_er[14][13]=0.0506322;
	new3jer_param[14][14]=-0.0520908;  new3jer_param_er[14][14]=0.000426188;
	new3jer_param[14][15]=-0.120105;  new3jer_param_er[14][15]=0.0264841;
	new3jer_param[14][16]=-6.1382e-05;  new3jer_param_er[14][16]=0.000142561;
	new3jer_param[14][17]=1.77937e-07;  new3jer_param_er[14][17]=1.92265e-07;
	new3jer_param[14][18]=59.1772;  new3jer_param_er[14][18]=2.19426;
	new3jer_param[14][19]=1.67807e+12;  new3jer_param_er[14][19]=1.28967e+13;
	new3jer_param[14][20]=-2.03984e+10;  new3jer_param_er[14][20]=1.61321e+11;
	new3jer_param[14][21]=8.28304;  new3jer_param_er[14][21]=2.26809;
	new3jer_param[14][22]=0.0251499;  new3jer_param_er[14][22]=0.00571798;
	new3jer_param[14][23]=-564.859;  new3jer_param_er[14][23]=243.287;
	new3jer_param[14][24]=0;  new3jer_param_er[14][24]=0;
	new3jer_param[14][25]=0;  new3jer_param_er[14][25]=0;
	new3jer_param[14][26]=0;  new3jer_param_er[14][26]=0;
	new3jer_param[14][27]=0;  new3jer_param_er[14][27]=0;
	new3jer_param[14][28]=0;  new3jer_param_er[14][28]=0;
	new3jer_param[14][29]=0;  new3jer_param_er[14][29]=0;
	new3jer_param[14][30]=0;  new3jer_param_er[14][30]=0;
*/

	if(i<15 && j<31)
	{
		if(code==0) value = new3jer_param[i][j];
		if(code==1) value = new3jer_param_er[i][j];
	}

	return value;
}


//----- Fit function for Gauss Sigma and Landau Sigma : sqrt(A0/Pt+A2) * step function combined
double FitFuncForSigmas(const float fJetE, const float fDetEta, const int iSystCode, 
					const JetFilterModuleV2::JERFuncTypes_t iSelFunc)
{
	if (! (iSelFunc == JetFilterModuleV2::JER_GAUSSSIGMA 
			|| iSelFunc == JetFilterModuleV2::JER_LANDAUSIGMA) )
	{
		std::string sMsg;
		sMsg += "Undefined function parameters request. please check!\n";
		sMsg += "You can access Gauss Sigma(1), Landau Sigma(3). Exiting now.";
		StdOut(__FILE__,__LINE__,3,sMsg);
		exit (1);
	}

	int iEtaBin = WhatJEReta(fDetEta);
	int iIndex  = -1;
	
	if (iSelFunc == JetFilterModuleV2::JER_GAUSSSIGMA)  //pick up Gaus Sigma parameterization
	{
		iIndex = 7;		
	} else if (iSelFunc == JetFilterModuleV2::JER_LANDAUSIGMA)  //pick up Landau Sigma parameterization
	{
		iIndex = 19;		
	}

	double par0=0,par1=0,par2=0,par3=0,par4=0;
	par0 = new3JERparam(0,iEtaBin, iIndex);
	par1 = new3JERparam(0,iEtaBin, iIndex+1);
	par2 = new3JERparam(0,iEtaBin, iIndex+2);
	par3 = new3JERparam(0,iEtaBin, iIndex+3);
	par4 = new3JERparam(0,iEtaBin, iIndex+4);
	
	double fitval = 0.0;
	double p1 = par0 / fJetE + par1;
	double p2 = par2 / fJetE + par3;

	double ws = 1 / (1 + exp((fJetE - par4)/25.));
	fitval = p2 * (1 - ws) + p1 * ws;

	return fitval;
}


//----- 2nd round fit func for Gauss Mean/Landau MPV/Gauss Norm
double FitFuncForMeanMPVNorm(const float fJetE, const float fDetEta, const int iSystCode, 
					const JetFilterModuleV2::JERFuncTypes_t iSelFunc)

{
	if ( ! (iSelFunc == JetFilterModuleV2::JER_GAUSSMEAN 
				|| iSelFunc == JetFilterModuleV2::JER_LANDAUMPV
				|| iSelFunc == JetFilterModuleV2::JER_GAUSSNORM) )
	{
		std::string sMsg;
		sMsg += "Undefined function parameters request. please check!\n";
		sMsg += "You can access Gauss Mean(0), Landau MPV(2),";
		sMsg += "and Gauss Norm(4). Exiting now.";
		StdOut(__FILE__,__LINE__,3,sMsg);
		exit (1);
	}

	int iEtaBin = WhatJEReta(fDetEta);
	int iIndex  = -1;
	
	if (iSelFunc == JetFilterModuleV2::JER_GAUSSMEAN)  //pick up Gaus Mean parameters
	{
		iIndex = 0;		
	} else if (iSelFunc == JetFilterModuleV2::JER_LANDAUMPV)  //pick up Landau MPV parameterization
	{
		iIndex = 12;		
	} else if (iSelFunc == JetFilterModuleV2::JER_GAUSSNORM)	//pick up Gaus Norm Parameterization
	{
		iIndex = 24;		
	}

	double par0=0,par1=0,par2=0,par3=0,par4=0,par5=0,par6=0;
	par0 = new3JERparam(0,iEtaBin, iIndex);
	par1 = new3JERparam(0,iEtaBin, iIndex+1);
	par2 = new3JERparam(0,iEtaBin, iIndex+2);
	par3 = new3JERparam(0,iEtaBin, iIndex+3);
	par4 = new3JERparam(0,iEtaBin, iIndex+4);
	par5 = new3JERparam(0,iEtaBin, iIndex+5);
	par6 = new3JERparam(0,iEtaBin, iIndex+6);
	
	double p1 = par0  + par1 * fJetE + par2 * fJetE * fJetE;
	double p2 = par3  + par4 * fJetE + par5 * fJetE * fJetE;

	double ws = 1 / (1 + exp((fJetE-par6)/10.5));
	double fitval = p1 * ws + p2 * (1 - ws);
	
	return fitval;
}


//----- end of 2nd round fit functions ---------------------------------------------------------------






//-------- 1st round of new fit functions ---------------

//----- Another new fit function: sqrt(A0/Pt+A2) * step function combined
double fitFunc3_new(const float fJetE, const float fDetEta, const int iSystCode, 
					const JetFilterModuleV2::JERFuncTypes_t iSelFunc)

{
	if (iSelFunc < JetFilterModuleV2::JER_GAUSSMEAN 
			|| iSelFunc > JetFilterModuleV2::JER_GAUSSNORM)
	{
		std::string sMsg;
		sMsg += "Undefined function parameters request. please check!\n";
		sMsg += "You can access Gauss Mean(0), Gauss Sigma(1), Landau MPV(2),";
		sMsg += "Landau Sigma(3), and Gauss Norm(4). Exiting now.";
		StdOut(__FILE__,__LINE__,3,sMsg);
		exit (1);
	}


	int iEtaBin = WhatJEReta(fDetEta);
	int iIndex=0;
	
	if (iSelFunc == JetFilterModuleV2::JER_GAUSSMEAN)  //pick up Gaus Mean parameters
	{
		iIndex = 0;		
	} else if (iSelFunc == JetFilterModuleV2::JER_GAUSSSIGMA)  //pick up Gaus Sigma parameterization
	{
		iIndex = 5;		
	} else if (iSelFunc == JetFilterModuleV2::JER_LANDAUMPV)  //pick up Landau MPV parameterization
	{
		iIndex = 10;		
	} else if (iSelFunc == JetFilterModuleV2::JER_LANDAUSIGMA)  //pick up Landau Sigma parameterization
	{
		iIndex = 15;		
	} else if (iSelFunc == JetFilterModuleV2::JER_GAUSSNORM)	//pick up Gaus Norm Parameterization
	{
		iIndex = 20;		
	}

	double p0=0,p1=0,p2=0,p3=0,p4=0;
	p0 = new2JERparam(0,iEtaBin, iIndex);
	p1 = new2JERparam(0,iEtaBin, iIndex+1);
	p2 = new2JERparam(0,iEtaBin, iIndex+2);
	p3 = new2JERparam(0,iEtaBin, iIndex+3);
	p4 = new2JERparam(0,iEtaBin, iIndex+4);
	
	double fitval = 0.0;
	double f1 = p0 * fJetE;
	double f2 = p1 * fJetE + p2 / fJetE;

	double ws = 1 / (1 + exp((fJetE-p4)/p3));
	fitval = f2 * (1 - ws) + f1 * ws;

	return fitval;
}

//-------------- fit function: poly3 * step func + poly2 * (step func)
// this is the new fit function used for 
// Gauss Mean, Landau MPV, and Gauss Norm
// iSelFunc = 0(Gauss Mean),1(Gaus Sigma), 2(Landau MPV), 3(Landau Sigma)
//            4(Gauss Norm or relative fraction)
double dJayFunc1(const float fJetE, const float fDetEta, const int iSystCode, 
					const JetFilterModuleV2::JERFuncTypes_t iSelFunc)
{

	if (iSelFunc < JetFilterModuleV2::JER_GAUSSMEAN 
		 || iSelFunc > JetFilterModuleV2::JER_GAUSSNORM)
	{
		std::string sMsg;
		sMsg += "Undefined function parameters request. please check!\n";
		sMsg += "You can access Gauss Mean(0), Gauss Sigma(1), Landau MPV(2),";
		sMsg += "Landau Sigma(3), and Gauss Norm(4). Exiting now.";
		StdOut(__FILE__,__LINE__,3,sMsg);
		exit (1);
	}

	double p0=0,p1=0,p2=0,p3=0,p4=0;

	int iEtaBin = WhatJEReta(fDetEta);
	int iIndex=0;
	
	if (iSelFunc == JetFilterModuleV2::JER_GAUSSMEAN)  //pick up Gaus Mean parameters
	{
		iIndex = 0;		
	} else if (iSelFunc == JetFilterModuleV2::JER_GAUSSSIGMA)  //pick up Gaus Sigma parameterization
	{
		iIndex = 5;		
	} else if (iSelFunc == JetFilterModuleV2::JER_LANDAUMPV)  //pick up Landau MPV parameterization
	{
		iIndex = 10;		
	} else if (iSelFunc == JetFilterModuleV2::JER_LANDAUSIGMA)  //pick up Landau Sigma parameterization
	{
		iIndex = 15;		
	} else if (iSelFunc == JetFilterModuleV2::JER_GAUSSNORM)	//pick up Gaus Norm Parameterization
	{
		iIndex = 20;		
	}

	p0 = newJERparam(0,iEtaBin, iIndex);
	p1 = newJERparam(0,iEtaBin, iIndex+1);
	p2 = newJERparam(0,iEtaBin, iIndex+2);
	p3 = newJERparam(0,iEtaBin, iIndex+3);
	p4 = newJERparam(0,iEtaBin, iIndex+4);

	double f1 = p0  + p1 * fJetE + p2 * fJetE * fJetE;
	double f2 = p3  + p4 * fJetE;

	double ws = 1 / (1 + exp((fJetE - 60.0)/10.5));
	double fitval = f1 * ws + f2 * (1 - ws);
	return fitval;
}

//-------- END of 1st round of new fit functions ---------------

//_________________ returns mean of Gauss in JER: A0+A1*Pt+A2/E
double MyJer_meanG(double jet_E, double eta_det, int stat_code, int syst_code, const int iUseNewFit=0)
{
	//iUseNewFit=0 (use old fit)
	//				=1 (use new for Gaus Mean)
	
	double val=0.0;
	int bin=WhatJEReta(eta_det);

	// choose between old and new parametization
	if (! iUseNewFit)
	{
		if (jet_E<60.0) //return exact bin values in this case
		{
			int iEBin = (int) jet_E;
			int iEtaBin = bin;
			//return JetFilterModuleV2::vJerBinValue.at(iEtaBin).at(0).at(iEBin);
			return oldJERBinValues(iEtaBin,0, iEBin);
		}
		
		double p1=JERparam(0,bin,0);
		double p2=JERparam(0,bin,1);
		double p3=JERparam(0,bin,2);

		double val_err=0.0;
		if(jet_E>0.0)
		{
			if(abs(syst_code)==1)
			{
				double t1=JERparam(1,bin,0);
				double t2=JERparam(1,bin,1)*jet_E;
				double t3=JERparam(1,bin,2)/jet_E;
				double t12=2.0*t1*t2*JERCorrCoeff(bin,0);
				double t13=2.0*t1*t3*JERCorrCoeff(bin,1);
				double t23=2.0*t2*t3*JERCorrCoeff(bin,2);
				double t=t1*t1+t2*t2+t3*t3+t12+t13+t23;
				if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
			}
			val = p1 + p2 * jet_E + p3 / jet_E + val_err;
			if(val<-1.0 || val>5.0) val=0.0;
		}
		
	} else {
		assert (iUseNewFit == 1 && "Wrong function choosen to get Gaus Mean from FitFuncForMeanMPVNorm(). Pl. Check!");
		val = FitFuncForMeanMPVNorm(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_GAUSSMEAN); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 1 && "Wrong function choosen to get Gaus Mean from dJayFunc1. Pl. Check!");
		//val = dJayFunc1(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_GAUSSMEAN); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 1 && "Wrong function choosen to get Gaus Mean from fitFunc3_new. Pl. Check!");
		//val = fitFunc3_new(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_GAUSSMEAN); // for now systcode is used as jer_stat_code
	}
	return val;
	
}

//_________________ returns sigma of Gauss in JER: sqrt(A3/E+A4)
double MyJer_sigmaG(double jet_E, double eta_det, int stat_code, int syst_code, const int iUseNewFit=0)
{
	//iUseNewFit=0 (use old fit)
	//				=2 (use new for Gaus Sigma)
	
	int bin=WhatJEReta(eta_det);
	double val=0.0;

	if (! iUseNewFit)
	{

		if (jet_E<60.0) //return exact bin values in this case
		{
			int iEBin = (int) jet_E;
			int iEtaBin = bin;
			//return JetFilterModuleV2::vJerBinValue.at(iEtaBin).at(1).at(iEBin);
			return oldJERBinValues(iEtaBin,1, iEBin);
		}
		
		double p1=JERparam(0,bin,3);
		double p2=JERparam(0,bin,4);
		double val_err=0.0;
		if(jet_E>0.0)
		{
			val=p1/jet_E+p2;
			if(val>0.0)
			{
				if(abs(syst_code)==2)
				{
					double t1=JERparam(1,bin,3)/jet_E;
					double t2=JERparam(1,bin,4);
					double t12=2.0*t1*t2*JERCorrCoeff(bin,3);
					double t=0.25*(t1*t1+t2*t2+t12)/val;
					if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
				}
				val=sqrt(val)+val_err;
			}
			if(val<0.0) val=0.0;
		}

	} else {
		assert (iUseNewFit == 2 && "Wrong function choosen to get Gaus Sigma from FitFuncForSigmas(). Pl. Check!");
		val = FitFuncForSigmas(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_GAUSSSIGMA); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 2 && "Wrong function choosen to get Gaus Sigma from dJayFunc1. Pl. Check!");
		//val = dJayFunc1(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_GAUSSSIGMA); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 2 && "Wrong function choosen to get Gaus Sigma from fitFunc3_new. Pl. Check!");
		//val = fitFunc3_new(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_GAUSSSIGMA); // for now systcode is used as jer_stat_code
	}
	return val;
}

//_________________ returns mpv of Landau in JER: A5+A6*Pt+A7/E
double MyJer_mpvL(double jet_E, double eta_det, int stat_code, int syst_code, const int iUseNewFit=0)
{

	//iUseNewFit=0 (use old fit)
	//				=3 (use new for Gaus Sigma)
	int bin=WhatJEReta(eta_det);
	double val=0.0;

	if (! iUseNewFit)
	{
		if (jet_E<60.0) //return exact bin values in this case
		{
			int iEBin = (int) jet_E;
			int iEtaBin = bin;
			//return JetFilterModuleV2::vJerBinValue.at(iEtaBin).at(2).at(iEBin);
			return oldJERBinValues(iEtaBin,2, iEBin);
		}
		
		double p1=JERparam(0,bin,5);
		double p2=JERparam(0,bin,6);
		double p3=JERparam(0,bin,7);
		double val_err=0.0;
		if(jet_E>0.0)
		{
			if(abs(syst_code)==3)
			{
				double t1=JERparam(1,bin,5);
				double t2=JERparam(1,bin,6)*jet_E;
				double t3=JERparam(1,bin,7)/jet_E;
				double t12=2.0*t1*t2*JERCorrCoeff(bin,4);
				double t13=2.0*t1*t3*JERCorrCoeff(bin,5);
				double t23=2.0*t2*t3*JERCorrCoeff(bin,6);
				double t=t1*t1+t2*t2+t3*t3+t12+t13+t23;
				if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
			}
			val=p1+p2*jet_E+p3/jet_E+val_err;
			if(val<-1.0 || val>5.0) val=0.0;
		}

	} else {
		assert (iUseNewFit == 3 && "Wrong function choosen to get  Landau MPV from FitFuncForMeanMPVNorm(). Pl. Check!");
		val = FitFuncForMeanMPVNorm(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_LANDAUMPV); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 3 && "Wrong function choosen to get  Landau MPV from dJayFunc1. Pl. Check!");
		//val = dJayFunc1(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_LANDAUMPV); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 3 && "Wrong function choosen to get  Landau MPV from fitFunc3_new. Pl. Check!");
		//val = fitFunc3_new(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_LANDAUMPV); // for now systcode is used as jer_stat_code
	}
	
	return val;
}

//_________________ returns sigma of Landau in JER: sqrt(A8/E+A9)
double MyJer_sigmaL(double jet_E, double eta_det, int stat_code, int syst_code, const int iUseNewFit=0)
{
	//iUseNewFit=0 (use old fit)
	//				=4 (use new for Gaus Sigma)
	int bin=WhatJEReta(eta_det);
	double val=0.0;

	if (! iUseNewFit)
	{
		if (jet_E<60.0) //return exact bin values in this case
		{
			int iEBin = (int) jet_E;
			int iEtaBin = bin;
			//return JetFilterModuleV2::vJerBinValue.at(iEtaBin).at(3).at(iEBin);
			return oldJERBinValues(iEtaBin,3, iEBin);
		}
		
		double p1=JERparam(0,bin,8);
		double p2=JERparam(0,bin,9);
		double val_err=0.0;
		if(jet_E>0.0)
		{
			val=p1/jet_E+p2;
			if(val>0.0)
			{
				if(abs(syst_code)==4)
				{
					double t1=JERparam(1,bin,8)/jet_E;
					double t2=JERparam(1,bin,9);
					double t12=2.0*t1*t2*JERCorrCoeff(bin,7);
					double t=0.25*(t1*t1+t2*t2+t12)/val;
					if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
				}
				val=sqrt(val)+val_err;
			}
			if(val<0.0) val=0.0;
		}
	} else {
		assert (iUseNewFit == 4 && "Wrong function choosen to get  Landau Sigma from FitFuncForSigmas(). Pl. Check!");
		val = FitFuncForSigmas(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_LANDAUSIGMA); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 4 && "Wrong function choosen to get  Landau Sigma from dJayFunc1. Pl. Check!");
		//val = dJayFunc1(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_LANDAUSIGMA); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 4 && "Wrong function choosen to get  Landau Sigma from fitFunc3_new. Pl. Check!");
		//val = fitFunc3_new(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_LANDAUSIGMA); // for now systcode is used as jer_stat_code
	}
	return val;
}

//_________________ returns normalization of Gauss: [A10+A11*sqrt(E)]/E+A12
double MyJer_normG(double jet_E, double eta_det, int stat_code, int syst_code, const int iUseNewFit=0)
{
	//iUseNewFit=0 (use old fit)
	//				=5 (use new for Gaus Sigma)
	int bin=WhatJEReta(eta_det);
	double val=0.0;

	if (! iUseNewFit)
	{
		if (jet_E<60.0) //return exact bin values in this case
		{
			int iEBin = (int) jet_E;
			int iEtaBin = bin;
			//return JetFilterModuleV2::vJerBinValue.at(iEtaBin).at(4).at(iEBin);
			return oldJERBinValues(iEtaBin,4, iEBin);
		}
		
		double p1=JERparam(0,bin,10);
		double p1_1=JetFilterModuleV2::vJerParam.at(bin).at(10);
		double p2=JERparam(0,bin,11);
		double p2_1=JetFilterModuleV2::vJerParam.at(bin).at(11);
		double p3=JERparam(0,bin,12);
		double p3_1=JetFilterModuleV2::vJerParam.at(bin).at(12);
		assert ((p1==p1_1 && p2==p2_1 && p3==p3_1) && "JER para error!");
		double val_err=0.0;
		if(jet_E>0.0)
		{
			if(abs(syst_code)==5)
			{
				double t1=JERparam(1,bin,10)/jet_E;
				double t2=JERparam(1,bin,11)/sqrt(jet_E);
				double t3=JERparam(1,bin,12);
				double t12=2.0*t1*t2*JERCorrCoeff(bin,8);
				double t13=2.0*t1*t3*JERCorrCoeff(bin,9);
				double t23=2.0*t2*t3*JERCorrCoeff(bin,10);
				double t=t1*t1+t2*t2+t3*t3+t12+t13+t23;
				if(t>0.0) val_err=(syst_code>0) ? sqrt(t) : -1.0*sqrt(t);
			}
			val=(p1+p2*sqrt(jet_E))/jet_E+p3+val_err;
		}
		if(val<0.0) val=0.0;

	} else {
		assert (iUseNewFit == 5 && "Wrong function choosen to get Gauss Norm from FitFuncForMeanMPVNorm(). Pl. Check!");
		val = FitFuncForMeanMPVNorm(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_GAUSSNORM); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 5 && "Wrong function choosen to get Gauss Norm from dJayFunc1. Pl. Check!");
		//val = dJayFunc1(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_GAUSSNORM); // for now systcode is used as jer_stat_code
		//assert (iUseNewFit == 5 && "Wrong function choosen to get Gauss Norm from fitFunc3_new. Pl. Check!");
		//val = fitFunc3_new(jet_E,eta_det,syst_code, JetFilterModuleV2::JER_GAUSSNORM); // for now systcode is used as jer_stat_code
	}
	return val;
}


//___ generates Sigma and Mean for "unclustered" Met
//___ takes into account correlations between params
//___ systcode=0 -default param
//___ systcode=1 -choice of parameterization: Z->ee vs. dipho-sideband
//___ systcode=2(+/-) -Met Model Uncl. En. systematics: mean+/-G
//___ systcode=3(+/-) -Met Model Uncl. En. systematics: sigma+/-G
//___ systcode=4(+/-) -Met Model Uncl. En. systematics: scale+/-G
//___ systcode=5(+/-) -Met Model Uncl. En. systematics: norm+/-G
void GetMyUnclMetResolution(double sumEt, int isData, int ParamSwitch, int systcode,
		double &sigmaMx, double &sigmaMy, double &meanMx, double &meanMy,
		double &normX, double &normY, double &scaleX, double &scaleY) {

	double sigmaX_err=0.0;
	double meanX_err=0.0;
	double sigmaY_err=0.0;
	double meanY_err=0.0;
	double normX_err=0.0;
	double normY_err=0.0;
	double scaleX_err=0.0;
	double scaleY_err=0.0;
	double dummypar=0.0;

	if(sumEt<0.0) sumEt=0.0;

	//--------------------------------------------------------------------------------
	//_______________ ***** met model parametrization.
	//                for MetX & MetY: sigma=p0+p1*sqrt(SumEt)
	//                for MetX & MetY: mean=p0+p1*SumEt
	//                first index is for di-pho sideband ([0]) & central Z->ee ([1]) parametrizations
	double metmodel_sigmaX[2][2];    // sigmaX=p[0]+p[1]*sqrt(corrSumEt)
	double metmodel_sigmaX_er[2][2];
	double metmodel_sigmaY[2][2];    // sigmaY=p[0]+p[1]*sqrt(corrSumEt)
	double metmodel_sigmaY_er[2][2];
	double metmodel_meanX[2][2];     // meanX=p[0]+p[1]*corrSumEt
	double metmodel_meanX_er[2][2];
	double metmodel_meanY[2][2];     // meanY=p[0]+p[1]*corrSumEt
	double metmodel_meanY_er[2][2];

	double metmodel_normX[2]; // normX=p[0]
	double metmodel_normX_er[2];
	double metmodel_normY[2]; // normY=p[0]
	double metmodel_normY_er[2];
	double metmodel_sigmaScaleX[2][2];    // sigmaScaleX=p[0]+p[1]*sqrt(corrSumEt)
	double metmodel_sigmaScaleX_er[2][2];
	double metmodel_sigmaScaleY[2][2];    // sigmaScaleY=p[0]+p[1]*sqrt(corrSumEt)
	double metmodel_sigmaScaleY_er[2][2];

	//-------------- For unclustered MetModel with Njet(cut)=0
	double metmodel_sigmaXcorr[2]; // correlation coefficient for sigmaX
	double metmodel_sigmaYcorr[2]; // correlation coefficient for sigmaY
	double metmodel_meanXcorr[2]; // correlation coefficient for meanX
	double metmodel_meanYcorr[2]; // correlation coefficient for meanY
	double metmodel_sigmaScaleXcorr[2]; // correlation coefficient for sigmaScaleX
	double metmodel_sigmaScaleYcorr[2]; // correlation coefficient for sigmaScaleY

	if(isData==1)   //---------- Data
	{
		//____________________________ Data di-pho sideband parameterization for events with Njet15 = 0
		//------------- parametrization after 09/12/07
		metmodel_sigmaX[0][0]=0.82;
		metmodel_sigmaX[0][1]=0.372;
		metmodel_sigmaX_er[0][0]=0.25;
		metmodel_sigmaX_er[0][1]=0.031;
		metmodel_sigmaY[0][0]=0.60;
		metmodel_sigmaY[0][1]=0.387;
		metmodel_sigmaY_er[0][0]=0.25;
		metmodel_sigmaY_er[0][1]=0.022;
		metmodel_meanX[0][0]=-0.022;
		metmodel_meanX[0][1]=0.00647;
		metmodel_meanX_er[0][0]=0.057;
		metmodel_meanX_er[0][1]=0.00076;
		metmodel_meanY[0][0]=-0.016;
		metmodel_meanY[0][1]=-0.00337;
		metmodel_meanY_er[0][0]=0.038;
		metmodel_meanY_er[0][1]=0.00033;

		metmodel_sigmaXcorr[0]=-0.911;
		metmodel_sigmaYcorr[0]=-0.906;
		metmodel_meanXcorr[0]=-0.814;
		metmodel_meanYcorr[0]=-0.764;

		metmodel_normX[0]=0.080;
		metmodel_normX_er[0]=0.022;
		metmodel_normY[0]=0.147;
		metmodel_normY_er[0]=0.036;
		metmodel_sigmaScaleX[0][0]=2.16;
		metmodel_sigmaScaleX[0][1]=-0.064;
		metmodel_sigmaScaleX_er[0][0]=0.17;
		metmodel_sigmaScaleX_er[0][1]=0.020;
		metmodel_sigmaScaleY[0][0]=1.99;
		metmodel_sigmaScaleY[0][1]=-0.046;
		metmodel_sigmaScaleY_er[0][0]=0.17;
		metmodel_sigmaScaleY_er[0][1]=0.020;
		metmodel_sigmaScaleXcorr[0]=-0.973;
		metmodel_sigmaScaleYcorr[0]=-0.976;

		//------------- Data Z->ee (CC) parameterization (as of 09/12/07)
		metmodel_sigmaX[1][0]=1.03;
		metmodel_sigmaX[1][1]=0.371;
		metmodel_sigmaX_er[1][0]=0.36;
		metmodel_sigmaX_er[1][1]=0.042;
		metmodel_sigmaY[1][0]=1.04;
		metmodel_sigmaY[1][1]=0.389;
		metmodel_sigmaY_er[1][0]=0.32;
		metmodel_sigmaY_er[1][1]=0.034;
		metmodel_meanX[1][0]=-0.048;
		metmodel_meanX[1][1]=0.00674;
		metmodel_meanX_er[1][0]=0.057;
		metmodel_meanX_er[1][1]=0.00073;
		metmodel_meanY[1][0]=-0.117;
		metmodel_meanY[1][1]=-0.00323;
		metmodel_meanY_er[1][0]=0.053;
		metmodel_meanY_er[1][1]=0.00073;

		metmodel_sigmaXcorr[1]=-0.921;
		metmodel_sigmaYcorr[1]=-0.903;
		metmodel_meanXcorr[1]=-0.774;
		metmodel_meanYcorr[1]=-0.827;

		metmodel_normX[1]=0.281;
		metmodel_normX_er[1]=0.076;
		metmodel_normY[1]=0.235;
		metmodel_normY_er[1]=0.078;
		metmodel_sigmaScaleX[1][0]=1.94;
		metmodel_sigmaScaleX[1][1]=-0.051;
		metmodel_sigmaScaleX_er[1][0]=0.11;
		metmodel_sigmaScaleX_er[1][1]=0.013;
		metmodel_sigmaScaleY[1][0]=2.19;
		metmodel_sigmaScaleY[1][1]=-0.079;
		metmodel_sigmaScaleY_er[1][0]=0.16;
		metmodel_sigmaScaleY_er[1][1]=0.019;
		metmodel_sigmaScaleXcorr[1]=-0.944;
		metmodel_sigmaScaleYcorr[1]=-0.965;
	}
	else //--------- MC
	{
		if(ParamSwitch==0) // 0=dipho parametrization
		{
			//____________________________ Pythia di-pho signal parametrization for events with Njet15 = 0
			//------------- parametrization after 06/19/08 (new Pythia with Min Bias)
			metmodel_sigmaX[0][0]=0.67;
			metmodel_sigmaX[0][1]=0.385;
			metmodel_sigmaX_er[0][0]=0.24;
			metmodel_sigmaX_er[0][1]=0.028;
			metmodel_sigmaY[0][0]=0.80;
			metmodel_sigmaY[0][1]=0.364;
			metmodel_sigmaY_er[0][0]=0.21;
			metmodel_sigmaY_er[0][1]=0.024;
			metmodel_meanX[0][0]=0.020;
			metmodel_meanX[0][1]=0.00166;
			metmodel_meanX_er[0][0]=0.023;
			metmodel_meanX_er[0][1]=0.00032;
			metmodel_meanY[0][0]=-0.012;
			metmodel_meanY[0][1]=0.00004;
			metmodel_meanY_er[0][0]=0.024;
			metmodel_meanY_er[0][1]=0.00031;

			metmodel_sigmaXcorr[0]=-0.916;
			metmodel_sigmaYcorr[0]=-0.912;
			metmodel_meanXcorr[0]=-0.845;
			metmodel_meanYcorr[0]=-0.870;

			metmodel_normX[0]=0.15;
			metmodel_normX_er[0]=0.02;
			metmodel_normY[0]=0.122;
			metmodel_normY_er[0]=0.018;
			metmodel_sigmaScaleX[0][0]=1.969;
			metmodel_sigmaScaleX[0][1]=-0.04819;
			metmodel_sigmaScaleX_er[0][0]=0.072;
			metmodel_sigmaScaleX_er[0][1]=0.0078;
			metmodel_sigmaScaleY[0][0]=2.013;
			metmodel_sigmaScaleY[0][1]=-0.0515;
			metmodel_sigmaScaleY_er[0][0]=0.075;
			metmodel_sigmaScaleY_er[0][1]=0.0075;
			metmodel_sigmaScaleXcorr[0]=-0.954;
			metmodel_sigmaScaleYcorr[0]=-0.962;

			//------------- Pythia Z->ee (CC) parametrization is the same as Pythia di-pho (for now)
			metmodel_sigmaX[1][0]=0.67;
			metmodel_sigmaX[1][1]=0.385;
			metmodel_sigmaX_er[1][0]=0.24;
			metmodel_sigmaX_er[1][1]=0.028;
			metmodel_sigmaY[1][0]=0.80;
			metmodel_sigmaY[1][1]=0.364;
			metmodel_sigmaY_er[1][0]=0.21;
			metmodel_sigmaY_er[1][1]=0.024;
			metmodel_meanX[1][0]=0.020;
			metmodel_meanX[1][1]=0.00166;
			metmodel_meanX_er[1][0]=0.023;
			metmodel_meanX_er[1][1]=0.00032;
			metmodel_meanY[1][0]=-0.012;
			metmodel_meanY[1][1]=0.00004;
			metmodel_meanY_er[1][0]=0.024;
			metmodel_meanY_er[1][1]=0.00031;

			metmodel_sigmaXcorr[1]=-0.916;
			metmodel_sigmaYcorr[1]=-0.912;
			metmodel_meanXcorr[1]=-0.845;
			metmodel_meanYcorr[1]=-0.870;

			metmodel_normX[1]=0.15;
			metmodel_normX_er[1]=0.02;
			metmodel_normY[1]=0.122;
			metmodel_normY_er[1]=0.018;
			metmodel_sigmaScaleX[1][0]=1.969;
			metmodel_sigmaScaleX[1][1]=-0.04819;
			metmodel_sigmaScaleX_er[1][0]=0.072;
			metmodel_sigmaScaleX_er[1][1]=0.0078;
			metmodel_sigmaScaleY[1][0]=2.013;
			metmodel_sigmaScaleY[1][1]=-0.0515;
			metmodel_sigmaScaleY_er[1][0]=0.075;
			metmodel_sigmaScaleY_er[1][1]=0.0075;
			metmodel_sigmaScaleXcorr[1]=-0.954;
			metmodel_sigmaScaleYcorr[1]=-0.962;
		}
		else // 1= cem-cem Z->ee parametrization
		{
			//------------- Pythia Z->ee (CC) parametrization (as of 09/11/07: double Gaussian)
			metmodel_sigmaX[1][0]=0.98;
			metmodel_sigmaX[1][1]=0.457;
			metmodel_sigmaX_er[1][0]=0.25;
			metmodel_sigmaX_er[1][1]=0.027;
			metmodel_sigmaY[1][0]=1.66;
			metmodel_sigmaY[1][1]=0.367;
			metmodel_sigmaY_er[1][0]=0.17;
			metmodel_sigmaY_er[1][1]=0.027;
			metmodel_meanX[1][0]=0.155;
			metmodel_meanX[1][1]=0.00137;
			metmodel_meanX_er[1][0]=0.020;
			metmodel_meanX_er[1][1]=0.00035;
			metmodel_meanY[1][0]=-0.179;
			metmodel_meanY[1][1]=-0.00333;
			metmodel_meanY_er[1][0]=0.026;
			metmodel_meanY_er[1][1]=0.00044;

			metmodel_sigmaXcorr[1]=-0.893;
			metmodel_sigmaYcorr[1]=-0.857;
			metmodel_meanXcorr[1]=-0.857;
			metmodel_meanYcorr[1]=-0.807;

			metmodel_normX[1]=0.110;
			metmodel_normX_er[1]=0.025;
			metmodel_normY[1]=0.134;
			metmodel_normY_er[1]=0.029;
			metmodel_sigmaScaleX[1][0]=1.755;
			metmodel_sigmaScaleX[1][1]=-0.046;
			metmodel_sigmaScaleX_er[1][0]=0.091;
			metmodel_sigmaScaleX_er[1][1]=0.011;
			metmodel_sigmaScaleY[1][0]=1.697;
			metmodel_sigmaScaleY[1][1]=-0.035;
			metmodel_sigmaScaleY_er[1][0]=0.083;
			metmodel_sigmaScaleY_er[1][1]=0.011;
			metmodel_sigmaScaleXcorr[1]=-0.966;
			metmodel_sigmaScaleYcorr[1]=-0.961;

			//----------------------- di-pho papamertization for studies with Pythia Zee events
			metmodel_sigmaX[0][0]=0.98;
			metmodel_sigmaX[0][1]=0.457;
			metmodel_sigmaX_er[0][0]=0.25;
			metmodel_sigmaX_er[0][1]=0.027;
			metmodel_sigmaY[0][0]=1.66;
			metmodel_sigmaY[0][1]=0.367;
			metmodel_sigmaY_er[0][0]=0.17;
			metmodel_sigmaY_er[0][1]=0.027;
			metmodel_meanX[0][0]=0.155;
			metmodel_meanX[0][1]=0.00137;
			metmodel_meanX_er[0][0]=0.020;
			metmodel_meanX_er[0][1]=0.00035;
			metmodel_meanY[0][0]=-0.179;
			metmodel_meanY[0][1]=-0.00333;
			metmodel_meanY_er[0][0]=0.026;
			metmodel_meanY_er[0][1]=0.00044;

			metmodel_sigmaXcorr[0]=-0.893;
			metmodel_sigmaYcorr[0]=-0.857;
			metmodel_meanXcorr[0]=-0.857;
			metmodel_meanYcorr[0]=-0.807;

			metmodel_normX[0]=0.110;
			metmodel_normX_er[0]=0.025;
			metmodel_normY[0]=0.134;
			metmodel_normY_er[0]=0.029;
			metmodel_sigmaScaleX[0][0]=1.755;
			metmodel_sigmaScaleX[0][1]=-0.046;
			metmodel_sigmaScaleX_er[0][0]=0.091;
			metmodel_sigmaScaleX_er[0][1]=0.011;
			metmodel_sigmaScaleY[0][0]=1.697;
			metmodel_sigmaScaleY[0][1]=-0.035;
			metmodel_sigmaScaleY_er[0][0]=0.083;
			metmodel_sigmaScaleY_er[0][1]=0.011;
			metmodel_sigmaScaleXcorr[0]=-0.966;
			metmodel_sigmaScaleYcorr[0]=-0.961;
		}
	}

	//___ systcode=4(+/-) -Met Model Uncl. En. systematics: scale+/-G
	//___ systcode=5(+/-) -Met Model Uncl. En. systematics: norm+/-G
	//_____________________ systcode=0 -default param
	if(systcode==0)
	{
		//___________ X-axis
		//-sigma
		sigmaMx=metmodel_sigmaX[ParamSwitch][0]+metmodel_sigmaX[ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMx<0.0) sigmaMx=0.0;
		//-mean
		meanMx=metmodel_meanX[ParamSwitch][0]+metmodel_meanX[ParamSwitch][1]*sumEt;
		//-norm
		normX=metmodel_normX[ParamSwitch];
		if(normX<0.0) normX=0.0;
		//-scale
		scaleX=metmodel_sigmaScaleX[ParamSwitch][0]+metmodel_sigmaScaleX[ParamSwitch][1]*sqrt(sumEt);
		if(scaleX<1.0) scaleX=1.0;
		//___________ Y-axis
		//-sigma
		sigmaMy=metmodel_sigmaY[ParamSwitch][0]+metmodel_sigmaY[ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMy<0.0) sigmaMy=0.0;
		//-mean
		meanMy=metmodel_meanY[ParamSwitch][0]+metmodel_meanY[ParamSwitch][1]*sumEt;
		//-norm
		normY=metmodel_normY[ParamSwitch];
		if(normY<0.0) normY=0.0;
		//-scale
		scaleY=metmodel_sigmaScaleY[ParamSwitch][0]+metmodel_sigmaScaleY[ParamSwitch][1]*sqrt(sumEt);
		if(scaleY<1.0) scaleY=1.0;
	}
	//________________________ systcode=1 -choice of parameterization: Z->ee vs. dipho-sideband
	if(abs(systcode)==1)
	{
		//___________ X-axis
		//-sigma
		sigmaMx=metmodel_sigmaX[1-ParamSwitch][0]+metmodel_sigmaX[1-ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMx<0.0) sigmaMx=0.0;
		//-mean
		meanMx=metmodel_meanX[1-ParamSwitch][0]+metmodel_meanX[1-ParamSwitch][1]*sumEt;
		//-norm
		normX=metmodel_normX[1-ParamSwitch];
		if(normX<0.0) normX=0.0;
		//-scale
		scaleX=metmodel_sigmaScaleX[1-ParamSwitch][0]+metmodel_sigmaScaleX[1-ParamSwitch][1]*sqrt(sumEt);
		if(scaleX<1.0) scaleX=1.0;
		//___________ Y-axis
		//-sigma
		sigmaMy=metmodel_sigmaY[1-ParamSwitch][0]+metmodel_sigmaY[1-ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMy<0.0) sigmaMy=0.0;
		//-mean
		meanMy=metmodel_meanY[1-ParamSwitch][0]+metmodel_meanY[1-ParamSwitch][1]*sumEt;
		//-norm
		normY=metmodel_normY[1-ParamSwitch];
		if(normY<0.0) normY=0.0;
		//-scale
		scaleY=metmodel_sigmaScaleY[1-ParamSwitch][0]+metmodel_sigmaScaleY[1-ParamSwitch][1]*sqrt(sumEt);
		if(scaleY<1.0) scaleY=1.0;
	}

	//___________________ calculating systematics
	if(abs(systcode)>1)
	{
		//___________ X-axis
		//-sigma
		dummypar=metmodel_sigmaX_er[ParamSwitch][0]*metmodel_sigmaX_er[ParamSwitch][0]
			+ metmodel_sigmaX_er[ParamSwitch][1]*metmodel_sigmaX_er[ParamSwitch][1]*sumEt
			+ 2.0*metmodel_sigmaXcorr[ParamSwitch]*metmodel_sigmaX_er[ParamSwitch][0]*metmodel_sigmaX_er[ParamSwitch][1]*sqrt(sumEt);
		if(dummypar<0.0) dummypar=0.0;
		sigmaX_err=sqrt(dummypar);
		//-mean
		dummypar=metmodel_meanX_er[ParamSwitch][0]*metmodel_meanX_er[ParamSwitch][0]
			+ metmodel_meanX_er[ParamSwitch][1]*metmodel_meanX_er[ParamSwitch][1]*sumEt*sumEt
			+ 2.0*metmodel_meanXcorr[ParamSwitch]*metmodel_meanX_er[ParamSwitch][0]*metmodel_meanX_er[ParamSwitch][1]*sumEt;
		if(dummypar<0.0) dummypar=0.0;
		meanX_err=sqrt(dummypar);
		//-norm
		dummypar=metmodel_normX_er[ParamSwitch]*metmodel_normX_er[ParamSwitch];
		normX_err=sqrt(dummypar);
		//-scale
		dummypar=metmodel_sigmaScaleX_er[ParamSwitch][0]*metmodel_sigmaScaleX_er[ParamSwitch][0]
			+ metmodel_sigmaScaleX_er[ParamSwitch][1]*metmodel_sigmaScaleX_er[ParamSwitch][1]*sumEt
			+ 2.0*metmodel_sigmaScaleXcorr[ParamSwitch]*metmodel_sigmaScaleX_er[ParamSwitch][0]*metmodel_sigmaScaleX_er[ParamSwitch][1]*sqrt(sumEt);
		if(dummypar<0.0) dummypar=0.0;
		scaleX_err=sqrt(dummypar);
		//___________ Y-axis
		//-sigma
		dummypar=metmodel_sigmaY_er[ParamSwitch][0]*metmodel_sigmaY_er[ParamSwitch][0]
			+ metmodel_sigmaY_er[ParamSwitch][1]*metmodel_sigmaY_er[ParamSwitch][1]*sumEt
			+ 2.0*metmodel_sigmaYcorr[ParamSwitch]*metmodel_sigmaY_er[ParamSwitch][0]*metmodel_sigmaY_er[ParamSwitch][1]*sqrt(sumEt);
		if(dummypar<0.0) dummypar=0.0;
		sigmaY_err=sqrt(dummypar);
		//-mean
		dummypar=metmodel_meanY_er[ParamSwitch][0]*metmodel_meanY_er[ParamSwitch][0]
			+ metmodel_meanY_er[ParamSwitch][1]*metmodel_meanY_er[ParamSwitch][1]*sumEt*sumEt
			+ 2.0*metmodel_meanYcorr[ParamSwitch]*metmodel_meanY_er[ParamSwitch][0]*metmodel_meanY_er[ParamSwitch][1]*sumEt;
		if(dummypar<0.0) dummypar=0.0;
		meanY_err=sqrt(dummypar);
		//-norm
		dummypar=metmodel_normY_er[ParamSwitch]*metmodel_normY_er[ParamSwitch];
		normY_err=sqrt(dummypar);
		//-scale
		dummypar=metmodel_sigmaScaleY_er[ParamSwitch][0]*metmodel_sigmaScaleY_er[ParamSwitch][0]
			+ metmodel_sigmaScaleY_er[ParamSwitch][1]*metmodel_sigmaScaleY_er[ParamSwitch][1]*sumEt
			+ 2.0*metmodel_sigmaScaleYcorr[ParamSwitch]*metmodel_sigmaScaleY_er[ParamSwitch][0]*metmodel_sigmaScaleY_er[ParamSwitch][1]*sqrt(sumEt);
		if(dummypar<0.0) dummypar=0.0;
		scaleY_err=sqrt(dummypar);
	}
	//____________________ systcode=2(+/-) -Met Model Uncl. En. systematics: mean+/-G
	if(abs(systcode)==2)
	{
		int int_systcode=(systcode>0) ? 1 : -1 ;
		//___________ X-axis
		//-sigma
		sigmaMx=metmodel_sigmaX[ParamSwitch][0]+metmodel_sigmaX[ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMx<0.0) sigmaMx=0.0;
		//-mean
		meanMx=metmodel_meanX[ParamSwitch][0]+metmodel_meanX[ParamSwitch][1]*sumEt+int_systcode*meanX_err;
		//-norm
		normX=metmodel_normX[ParamSwitch];
		if(normX<0.0) normX=0.0;
		//-scale
		scaleX=metmodel_sigmaScaleX[ParamSwitch][0]+metmodel_sigmaScaleX[ParamSwitch][1]*sqrt(sumEt);
		if(scaleX<1.0) scaleX=1.0;
		//___________ Y-axis
		//-sigma
		sigmaMy=metmodel_sigmaY[ParamSwitch][0]+metmodel_sigmaY[ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMy<0.0) sigmaMy=0.0;
		//-mean
		meanMy=metmodel_meanY[ParamSwitch][0]+metmodel_meanY[ParamSwitch][1]*sumEt+int_systcode*meanY_err;
		//-norm
		normY=metmodel_normY[ParamSwitch];
		if(normY<0.0) normY=0.0;
		//-scale
		scaleY=metmodel_sigmaScaleY[ParamSwitch][0]+metmodel_sigmaScaleY[ParamSwitch][1]*sqrt(sumEt);
		if(scaleY<1.0) scaleY=1.0;
	}
	//____________________ systcode=3(+/-) -Met Model Uncl. En. systematics: sigma+/-G
	if(abs(systcode)==3)
	{
		int int_systcode=(systcode>0) ? 1 : -1 ;
		//___________ X-axis
		//-sigma
		sigmaMx=metmodel_sigmaX[ParamSwitch][0]+metmodel_sigmaX[ParamSwitch][1]*sqrt(sumEt)+int_systcode*sigmaX_err;
		if(sigmaMx<0.0) sigmaMx=0.0;
		//-mean
		meanMx=metmodel_meanX[ParamSwitch][0]+metmodel_meanX[ParamSwitch][1]*sumEt;
		//-norm
		normX=metmodel_normX[ParamSwitch];
		if(normX<0.0) normX=0.0;
		//-scale
		scaleX=metmodel_sigmaScaleX[ParamSwitch][0]+metmodel_sigmaScaleX[ParamSwitch][1]*sqrt(sumEt);
		if(scaleX<1.0) scaleX=1.0;
		//___________ Y-axis
		//-sigma
		sigmaMy=metmodel_sigmaY[ParamSwitch][0]+metmodel_sigmaY[ParamSwitch][1]*sqrt(sumEt)+int_systcode*sigmaY_err;
		if(sigmaMy<0.0) sigmaMy=0.0;
		//-mean
		meanMy=metmodel_meanY[ParamSwitch][0]+metmodel_meanY[ParamSwitch][1]*sumEt;
		//-norm
		normY=metmodel_normY[ParamSwitch];
		if(normY<0.0) normY=0.0;
		//-scale
		scaleY=metmodel_sigmaScaleY[ParamSwitch][0]+metmodel_sigmaScaleY[ParamSwitch][1]*sqrt(sumEt);
		if(scaleY<1.0) scaleY=1.0;
	}
	//____________________ systcode=4(+/-) -Met Model Uncl. En. systematics: scale+/-G
	if(abs(systcode)==4)
	{
		int int_systcode=(systcode>0) ? 1 : -1 ;
		//___________ X-axis
		//-sigma
		sigmaMx=metmodel_sigmaX[ParamSwitch][0]+metmodel_sigmaX[ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMx<0.0) sigmaMx=0.0;
		//-mean
		meanMx=metmodel_meanX[ParamSwitch][0]+metmodel_meanX[ParamSwitch][1]*sumEt;
		//-norm
		normX=metmodel_normX[ParamSwitch];
		if(normX<0.0) normX=0.0;
		//-scale
		scaleX=metmodel_sigmaScaleX[ParamSwitch][0]+metmodel_sigmaScaleX[ParamSwitch][1]*sqrt(sumEt)+int_systcode*scaleX_err;
		if(scaleX<1.0) scaleX=1.0;
		//___________ Y-axis
		//-sigma
		sigmaMy=metmodel_sigmaY[ParamSwitch][0]+metmodel_sigmaY[ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMy<0.0) sigmaMy=0.0;
		//-mean
		meanMy=metmodel_meanY[ParamSwitch][0]+metmodel_meanY[ParamSwitch][1]*sumEt;
		//-norm
		normY=metmodel_normY[ParamSwitch];
		if(normY<0.0) normY=0.0;
		//-scale
		scaleY=metmodel_sigmaScaleY[ParamSwitch][0]+metmodel_sigmaScaleY[ParamSwitch][1]*sqrt(sumEt)+int_systcode*scaleY_err;
		if(scaleY<1.0) scaleY=1.0;
	}
	//____________________ systcode=5(+/-) -Met Model Uncl. En. systematics: scale+/-G
	if(abs(systcode)==5)
	{
		int int_systcode=(systcode>0) ? 1 : -1 ;
		//___________ X-axis
		//-sigma
		sigmaMx=metmodel_sigmaX[ParamSwitch][0]+metmodel_sigmaX[ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMx<0.0) sigmaMx=0.0;
		//-mean
		meanMx=metmodel_meanX[ParamSwitch][0]+metmodel_meanX[ParamSwitch][1]*sumEt;
		//-norm
		normX=metmodel_normX[ParamSwitch]+int_systcode*normX_err;
		if(normX<0.0) normX=0.0;
		//-scale
		scaleX=metmodel_sigmaScaleX[ParamSwitch][0]+metmodel_sigmaScaleX[ParamSwitch][1]*sqrt(sumEt);
		if(scaleX<1.0) scaleX=1.0;
		//___________ Y-axis
		//-sigma
		sigmaMy=metmodel_sigmaY[ParamSwitch][0]+metmodel_sigmaY[ParamSwitch][1]*sqrt(sumEt);
		if(sigmaMy<0.0) sigmaMy=0.0;
		//-mean
		meanMy=metmodel_meanY[ParamSwitch][0]+metmodel_meanY[ParamSwitch][1]*sumEt;
		//-norm
		normY=metmodel_normY[ParamSwitch]+int_systcode*normY_err;
		if(normY<0.0) normY=0.0;
		//-scale
		scaleY=metmodel_sigmaScaleY[ParamSwitch][0]+metmodel_sigmaScaleY[ParamSwitch][1]*sqrt(sumEt);
		if(scaleY<1.0) scaleY=1.0;
	}
	return;
}

//=====================================================================
//-------------- raw MetSig fit function: p0*exp(-p1*x+p2*x*x)+p3*exp(-p4*x)
double rawMetSigFunc(double *x, double *par)
{
	double fitval=0.0;
	if(x[0]>0.0)
	{
		fitval=par[0]*exp(-1.0*par[1]*x[0]+par[2]*x[0]*x[0])+par[3]*exp(-1.0*par[4]*x[0]);
	}
	else fitval=-1.0E6;
	return fitval;
}


ClassImp(JetFilterModuleV2)
//_____________________________________________________________________________
JetFilterModuleV2::JetFilterModuleV2(const char* name, const char* title):
		TStnModule(name,title),
		bRunPermit(true),
		debug(false),
		bUseHadJets(0),
	   bNoSummary(false),
		bDoubleSmearJets(false),
		fMetCut(0),
		iPrintLvl(100),
		iMetAddMethod(0),
		fMinisculeEta(-1e11),
		bRem_TightPho(0), bRem_LoosePho(0), bRem_SidebandPho(0),
		bRem_TightPhoLikeEle(0), bRem_LoosePhoLikeEle(0), bRem_StdLooseEle(0)
{

	//initilize the selection of JER parameterization old or new
	abUseJERnew.GaussMean_New   = JER_OLDFITS;
	abUseJERnew.GaussSigma_New  = JER_OLDFITS;
	abUseJERnew.LandauMPV_New   = JER_OLDFITS;
	abUseJERnew.LandauSigma_New = JER_OLDFITS;
	abUseJERnew.GaussNorm_New   = JER_OLDFITS;

		
	fUnclParamSwitch=0; // 0=sideband (default); 1=Z->ee unclustered energy parametrizations

	Npoints=1; // number of random points to generate

	//______________________________________ setting default cut values
	fEventWeight=1.0; // event weight for background samples
	fUseEventWeight=0; // code to turn ON/OFF event weight
	fAnalysisMode=0; // 0==OFF, 1==ON
	fUseVerbose=0;
	fDumpEvent=0;
	fSampleID=0; // need to set it from *.C script
	fJetAlgo=0;      // by default it's JetClu(0.4)
	fMyMetScenario=3; // by default met is corrected for jets with Et>15
	MetToy_max=200.0; // max value of toy Met
	MetToy_min=0.0; // min value of toy Met
	fSelectSigMetEvent=0; // 0= do nothing; -1=reject event; 1=select event
	fMetSig_cut=4.0; // MetSig=4.0 by default (i.e., 10E-4 events with fake MET to be selected for true MetSig)

	MinEtThr=15.0;
	MaxEtThr=1000.0;
	MinNjetThr=0;
	MaxNjetThr=100;
	MinJetEta=0.0;
	MaxJetEta=3.0;
	MinEt1st=0.0;
	MaxEt1st=1000.0;
	MinEt2nd=0.0;
	MaxEt2nd=1000.0;
	MinDeltaRgjMin=0.0;
	MaxDeltaRgjMin=1.0E6;
	MinDeltaEtagjMin=0.0;
	MaxDeltaEtagjMin=1.0E6;
	MinDeltaPhigjMin=0.0;
	MaxDeltaPhigjMin=1.0E6;
	MinDeltaRgj=0.0;
	MaxDeltaRgj=1.0E6;
	MinDeltaEtagj=0.0;
	MaxDeltaEtagj=1.0E6;
	MinDeltaPhigj=0.0;
	MaxDeltaPhigj=1.0E6;
	//   MinDeltaPhiJMet=0.056*2.0*asin(1.0); // ~10 degrees
	//   MaxDeltaPhiJMet=0.944*2.0*asin(1.0); // ~170 degrees
	MinDeltaPhiJMet=-10.0;
	MaxDeltaPhiJMet=0.3;
	MinMet2Et=0.0;
	MaxMet2Et=1.0E6;
	MinMet2Etlev6=0.0;
	MaxMet2Etlev6=1.0E6;
	MindZ=0.0;
	MaxdZ=1.0E6;
	MinMjj=0.0;
	MaxMjj=1.0E6;
	MinMj1g1=0.0;
	MaxMj1g1=1.0E6;
	MinMj1g2=0.0;
	MaxMj1g2=1.0E6;
	MinMj2g1=0.0;
	MaxMj2g1=1.0E6;
	MinMj2g2=0.0;
	MaxMj2g2=1.0E6;
	MinMjg=0.0;
	MaxMjg=1.0E6;
	MinQtJet=0.0;
	MaxQtJet=1.0E6;
	MinHt=0.0;
	MaxHt=1.0E6;
	MinNjet=0;
	MaxNjet=100;
	MinNjet5=0;
	MaxNjet5=100;
	MinNjet10=0;
	MaxNjet10=100;
	MinNjet15=0;
	MaxNjet15=100;
	MinNjet20=0;
	MaxNjet20=100;
	MinNjet25=0;
	MaxNjet25=100;
	MinNjet30=0;
	MaxNjet30=100;
	MinNjet35=0;
	MaxNjet35=100;
	fSelectMetEvent=0; // 0=no selection; 1=select events
	fSelectExoticEvent=0; // 0=no selection; 1=select events
	fRemoveDuplicate=0; // 0=no selection; 1=remove events; -1=select event
	fRemoveBadMet=0; // 0=do nothing; 1=remove event; -1=select event
	fDoVxSwap=0; // 0=do nothing; 1= swap vertices to minimize MET
	//________________________
	//________________________
	fDatFileName="MyJetFilter_output.dat";
	fRearEventFileName="exoticEventList.dat";
	fExoticEventFileName="veryExoticEventList.dat";
	fLargeMetEventFileName="largeMetEventList.dat";
	fDumpEventFileName="dumpEvent.dat";
	//________________________
	fJTC_coneSize=1;  // JetClue R=0.7
	fJTC_version=5;   // version of correction, suggested by Anwar for 5.3.1, may change later
	fJTC_level=7;     // full correction
	fJTC_systcode=0;  // default correction
	//fJTC_imode=1;     // DATA is default  //this will be automatically set at run time - sam 12-09-2008
	fJTC_imode=-1;     // DATA is default  //this will be automatically set at run time - sam 12-09-2008
	std::cout<<"Hello I am JetFilterModuleV2"<<std::endl;
}

//_____________________________________________________________________________
JetFilterModuleV2::~JetFilterModuleV2() {
}

//_____________________________________________________________________________
void JetFilterModuleV2::BookHistograms() {
	char folder_name[200];
	TFolder* fol;
	TFolder* hist_folder;

	//-----------------------------------------------------------------------------
	//  clear the histogram list
	//-----------------------------------------------------------------------------
	DeleteHistograms();
	hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

	//------ booking generic histograms


	sprintf(folder_name,"JetClu-0.4"); //----- Cone 0.4
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	BookJetHistograms(fHistJet04,Form("Hist/%s",folder_name),folder_name);

	//------ booking matching histograms
	sprintf(folder_name,"MatchPhoJet0.4"); //----- Cone 0.4
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	BookMatchStudyHistograms(fMatchPhoJet04,Form("Hist/%s",folder_name),folder_name);

	
	if(fAnalysisMode==1)
	{
		//________________ data
		sprintf(folder_name,"Ana_data");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_data,Form("Hist/%s",folder_name));

		//________________ Met Model total prediction (stat+syst)
		sprintf(folder_name,"Ana_bckg");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_bckg,Form("Hist/%s",folder_name));
		//________________ Met Model default prediction
		sprintf(folder_name,"Ana_def");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_def,Form("Hist/%s",folder_name));

		//________________ Met Model systematics: z->ee vs. dipho sideband parameterization
		sprintf(folder_name,"Ana_ue1");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_ue1,Form("Hist/%s",folder_name));
		//________________ Met Model Uncl. En. systematics: mean-G
		sprintf(folder_name,"Ana_ue2");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_ue2,Form("Hist/%s",folder_name));
		//________________ Met Model Uncl. En. systematics: mean+G
		sprintf(folder_name,"Ana_ue3");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_ue3,Form("Hist/%s",folder_name));
		//________________ Met Model Uncl. En. systematics: sigma-G
		sprintf(folder_name,"Ana_ue4");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_ue4,Form("Hist/%s",folder_name));
		//________________ Met Model Uncl. En. systematics: sigma+G
		sprintf(folder_name,"Ana_ue5");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_ue5,Form("Hist/%s",folder_name));
		//________________ Met Model Uncl. En. systematics: scale-G
		sprintf(folder_name,"Ana_ue6");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_ue6,Form("Hist/%s",folder_name));
		//________________ Met Model Uncl. En. systematics: scale+G
		sprintf(folder_name,"Ana_ue7");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_ue7,Form("Hist/%s",folder_name));
		//________________ Met Model Uncl. En. systematics: norm-G
		sprintf(folder_name,"Ana_ue8");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_ue8,Form("Hist/%s",folder_name));
		//________________ Met Model Uncl. En. systematics: norm+G
		sprintf(folder_name,"Ana_ue9");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_ue9,Form("Hist/%s",folder_name));


		//________________ Met Model JER systematics: meanG-G
		sprintf(folder_name,"Ana_jer1");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer1,Form("Hist/%s",folder_name));
		//________________ Met Model JER systematics: meanG+G
		sprintf(folder_name,"Ana_jer2");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer2,Form("Hist/%s",folder_name));
		//________________ Met Model JER systematics: sigmaG-G
		sprintf(folder_name,"Ana_jer3");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer3,Form("Hist/%s",folder_name));
		//________________ Met Model JER systematics: sigmaG+G
		sprintf(folder_name,"Ana_jer4");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer4,Form("Hist/%s",folder_name));
		//________________ Met Model JER systematics: mpvG-G
		sprintf(folder_name,"Ana_jer5");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer5,Form("Hist/%s",folder_name));
		//________________ Met Model JER systematics: mpvG+G
		sprintf(folder_name,"Ana_jer6");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer6,Form("Hist/%s",folder_name));
		//________________ Met Model JER systematics: sigmaL-G
		sprintf(folder_name,"Ana_jer7");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer7,Form("Hist/%s",folder_name));
		//________________ Met Model JER systematics: sigmaL+G
		sprintf(folder_name,"Ana_jer8");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer8,Form("Hist/%s",folder_name));
		//________________ Met Model JER systematics: norm-G
		sprintf(folder_name,"Ana_jer9");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer9,Form("Hist/%s",folder_name));
		//________________ Met Model JER systematics: norm+G
		sprintf(folder_name,"Ana_jer10");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookAnalysisHistograms(fAna_jer10,Form("Hist/%s",folder_name));

		if(fUseEventWeight==1)
		{
			FinalAnalysisHistoStep1(fAna_data);
			FinalAnalysisHistoStep1(fAna_def);
		}
	}
	else if (fAnalysisMode == 0)
	{
		//_______________________ booking met study histograms
		//----- Cone 0.4, Et>15
		sprintf(folder_name,"METJetClu0.4_scen3");
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookMetStudyHistograms(fMetStudyJet04sc3,Form("Hist/%s",folder_name),folder_name);
		
		//------- booking met cleanup histograms
		sprintf(folder_name,"MetCleanup0.4"); //----- Cone 0.4
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookMetCleanupHistograms(fCleanup,Form("Hist/%s",folder_name),folder_name);

		//------- booking generated met cleanup histograms
		sprintf(folder_name,"GenMetCleanup0.4"); //----- Cone 0.4, default MetModel
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookMetCleanupHistograms(fGenCleanup,Form("Hist/%s",folder_name),folder_name);

		sprintf(folder_name,"PhoMetStudy0.4"); //----- Cone 0.4, correlation between Photon and MET
		fol = (TFolder*) hist_folder->FindObject(folder_name);
		if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
		BookPhoMetStudyHistograms(fPhoMetStudy,Form("Hist/%s",folder_name));
	}

	/*
	//this is entirely Sam's stuff.
	char name [200];
	char title[200];
	TFolder *new_folder, *sub_folder;
	Histograms HistoMan;

	new_folder = GetHistFolder(this, "CountingRoom","Counting Histograms");
	vCounterHistLabels.push_back("Events Processed");				// 0
	vCounterHistLabels.push_back("Events With a Photon");			// 1
	vCounterHistLabels.push_back("Events Passed JetFilter");		// 2
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
	 */


	return;
}

//______________________ Booking final analysis histograms
void JetFilterModuleV2::BookAnalysisHistograms(AnalysisHisto_t& Hist, const char* Folder)
{
	char name [200];
	char title[200];

	sprintf(name,"MetAll");
	sprintf(title,"#slash{E}_{T}, all events");
	Hist.fAna_MetAll=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_MetAll,Folder);
	sprintf(name,"MetSig");
	sprintf(title,"#slash{E}_{T}-significance, all events");
	Hist.fAna_MetSig=new TH1F(name,title,200,0.0,20.0);
	AddHistogram(Hist.fAna_MetSig,Folder);

	sprintf(name,"MetAllVsMetSig");
	sprintf(title,"#slash{E}_{T} vs #slash{E}_{T}-significance all events");
	Hist.fAna_MetAllVsMetSig=new TH2F(name,title,200,0,20,5000,0,500);
	Hist.fAna_MetAllVsMetSig->GetXaxis()->SetTitle("#slash{E}_{T}-significance");
	Hist.fAna_MetAllVsMetSig->GetYaxis()->SetTitle("#slash{E}_{T} (GeV)");
	AddHistogram(Hist.fAna_MetAllVsMetSig,Folder);

	// the following histrograms are only for events that pass MetSig cut
	sprintf(name,"Met");
	sprintf(title,"#slash{E}_{T}, after cuts");
	Hist.fAna_Met=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Met,Folder);
	sprintf(name,"Njet15");
	sprintf(title,"N_{jet}(E_{T}^{jet}>15 GeV), after cuts");
	Hist.fAna_Njet15= new TH1F(name,title,21,-0.5,20.5);
	AddHistogram(Hist.fAna_Njet15,Folder);

	sprintf(name,"M");
	sprintf(title,"Mass of two leading EM objects, after cuts");
	Hist.fAna_M=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_M,Folder);
	sprintf(name,"dPhi");
	sprintf(title,"#Delta#phi of two leading EM objects, after cuts");
	Hist.fAna_dPhi= new TH1F(name,title,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi,Folder);
	sprintf(name,"Njet20");
	sprintf(title,"N_{jet}(E_{T}^{jet}>20 GeV), after cuts");
	Hist.fAna_Njet20= new TH1F(name,title,21,-0.5,20.5);
	AddHistogram(Hist.fAna_Njet20,Folder);
	sprintf(name,"Njet25");
	sprintf(title,"N_{jet}(E_{T}^{jet}>25 GeV), after cuts");
	Hist.fAna_Njet25= new TH1F(name,title,21,-0.5,20.5);
	AddHistogram(Hist.fAna_Njet25,Folder);
	sprintf(name,"Njet30");
	sprintf(title,"N_{jet}(E_{T}^{jet}>30 GeV), after cuts");
	Hist.fAna_Njet30= new TH1F(name,title,21,-0.5,20.5);
	AddHistogram(Hist.fAna_Njet30,Folder);
	sprintf(name,"Njet35");
	sprintf(title,"N_{jet}(E_{T}^{jet}>35 GeV), after cuts");
	Hist.fAna_Njet35= new TH1F(name,title,21,-0.5,20.5);
	AddHistogram(Hist.fAna_Njet35,Folder);
	sprintf(name,"Qt");
	sprintf(title,"Q_{T} of two leading EM objects, after cuts");
	Hist.fAna_Qt=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Qt,Folder);
	sprintf(name,"Et1");
	sprintf(title,"E_{T} of 1^{st} EM object, after cuts");
	Hist.fAna_Et1=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Et1,Folder);
	sprintf(name,"Et2");
	sprintf(title,"E_{T} of 2^{nd} EM object, after cuts");
	Hist.fAna_Et2=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Et2,Folder);
	sprintf(name,"Etjet");
	sprintf(title,"E_{T} of 1^{st} jet, after cuts");
	Hist.fAna_Etjet=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Etjet,Folder);
	sprintf(name,"Etjet2");
	sprintf(title,"E_{T} of 2^{nd} jet, after cuts");
	Hist.fAna_Etjet2=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Etjet2,Folder);

	sprintf(name,"Ht");
	sprintf(title,"H_{T}, after cuts");
	Hist.fAna_Ht=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Ht,Folder);
	sprintf(name,"Mjj");
	sprintf(title,"M_{jj}-- mass of two leading jets, after cuts");
	Hist.fAna_Mjj=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Mjj,Folder);
	sprintf(name,"Nem");
	sprintf(title,"N_{em}-- number of EM objects (e/#gamma), after cuts");
	Hist.fAna_Nem= new TH1F(name,title,11,-0.5,10.5);
	AddHistogram(Hist.fAna_Nem,Folder);
	sprintf(name,"Mej");
	sprintf(title,"M_{em-jet}, after cuts");
	Hist.fAna_Mej=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Mej,Folder);
	sprintf(name,"Mextra");
	sprintf(title,"M_{emL-emE}-- mass of leadingEM-extraEM if N_{em}>2, after cuts");
	Hist.fAna_Mextra=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Mextra,Folder);
	sprintf(name,"Etem");
	sprintf(title,"E_{T} of extra EM objects, after cuts");
	Hist.fAna_Etem=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Etem,Folder);
	sprintf(name,"dPhi1");
	sprintf(title,"#Delta#phi(#slash{E}_{T}-em1), after cuts");
	Hist.fAna_dPhi1= new TH1F(name,title,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi1,Folder);
	sprintf(name,"dPhi2");
	sprintf(title,"#Delta#phi(#slash{E}_{T}-em2), after cuts");
	Hist.fAna_dPhi2= new TH1F(name,title,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi2,Folder);
	sprintf(name,"dPhi3");
	sprintf(title,"#Delta#phi(#slash{E}_{T}-jet1) if E_{T}^{jet}>15 GeV, after cuts");
	Hist.fAna_dPhi3= new TH1F(name,title,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi3,Folder);
	sprintf(name,"dPhi4");
	sprintf(title,"#Delta#phi(#slash{E}_{T}-jet1) if E_{T}^{jet}>20 GeV, after cuts");
	Hist.fAna_dPhi4= new TH1F(name,title,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi4,Folder);
	sprintf(name,"dPhi5");
	sprintf(title,"#Delta#phi(#slash{E}_{T}-jet1) if E_{T}^{jet}>25 GeV, after cuts");
	Hist.fAna_dPhi5= new TH1F(name,title,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi5,Folder);
	sprintf(name,"dPhi6");
	sprintf(title,"#Delta#phi(#slash{E}_{T}-jet1) if E_{T}^{jet}>30 GeV, after cuts");
	Hist.fAna_dPhi6= new TH1F(name,title,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi6,Folder);
	sprintf(name,"dPhi7");
	sprintf(title,"#Delta#phi(#slash{E}_{T}-jet1) if E_{T}^{jet}>35 GeV, after cuts");
	Hist.fAna_dPhi7= new TH1F(name,title,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi7,Folder);

	sprintf(name,"jet1Eta");
	sprintf(title,"Lead Jet 4-vec #eta, after cuts");
	Hist.fAna_Jet1Eta= new TH1F(name,title,160,-4.0,4.0);
	AddHistogram(Hist.fAna_Jet1Eta,Folder);
	sprintf(name,"jet2Eta");
	sprintf(title,"Subleading Jet 4-vec #eta, after cuts");
	Hist.fAna_Jet2Eta= new TH1F(name,title,160,-4.0,4.0);
	AddHistogram(Hist.fAna_Jet2Eta,Folder);

	sprintf(name,"pho1Eta");
	sprintf(title,"Lead #gamma 4-vec #eta, after cuts");
	Hist.fAna_Pho1Eta= new TH1F(name,title,160,-4.0,4.0);
	AddHistogram(Hist.fAna_Pho1Eta,Folder);

	sprintf(name,"2dDalitz");
	sprintf(title,"Dalitz plot, after cuts");
	Hist.fAna_Dalitz= new TH2F(name,title,500,0,1000,500,0,1000);
	Hist.fAna_Dalitz->GetXaxis()->SetTitle("Inv. Mass(#gamma,Lead Jet)");
	Hist.fAna_Dalitz->GetYaxis()->SetTitle("Inv. Mass(#gamma,Sublead Jet)");
	AddHistogram(Hist.fAna_Dalitz,Folder);


	//for debugging only
	

	sprintf(name,"debug_PhoPtRatio_DelPhiJetMetLow");
	sprintf(title,"#gamma P_{T} Ratio for #Delta#Phi^{(Lead Jet,#slash{E}_{T})}<0.4, after cuts");
	Hist.hDebug_PtRatio_JMetdelphiLow= new TH1F(name,title,1000,-1.0,1.0);
	Hist.hDebug_PtRatio_JMetdelphiLow->GetXaxis()->SetTitle("Lead DET. #gamma/HEPG #gamma - 1");
	AddHistogram(Hist.hDebug_PtRatio_JMetdelphiLow,Folder);

	sprintf(name,"debug_PhoJetPtRatio_DelPhiJetMetLow");
	sprintf(title,"Lead Det.Jet/Hepg #gamma P_{T} Ratio for #Delta#Phi^{(Lead Jet,#slash{E}_{T})}<0.4, after cuts");
	Hist.hDebug_PhoJetPtRatio_JMetdelphiLow= new TH1F(name,title,1000,-2.0,2.0);
	Hist.hDebug_PhoJetPtRatio_JMetdelphiLow->GetXaxis()->SetTitle("Lead DET Jet/HEPG #gamma - 1");
	AddHistogram(Hist.hDebug_PhoJetPtRatio_JMetdelphiLow,Folder);

	sprintf(name,"debug_PhoPtRatio_DelPhiJetMetMid");
	sprintf(title,"#gamma P_{T} Ratio for 1.3<#Delta#Phi^{(Lead Jet,#slash{E}_{T})}<1.7, after cuts");
	Hist.hDebug_PtRatio_JMetdelphiMid= new TH1F(name,title,1000,-1.0,1.0);
	Hist.hDebug_PtRatio_JMetdelphiMid->GetXaxis()->SetTitle("Lead DET. #gamma/HEPG #gamma - 1");
	AddHistogram(Hist.hDebug_PtRatio_JMetdelphiMid,Folder);

	sprintf(name,"debug_PhoJetPtRatio_DelPhiJetMetMid");
	sprintf(title,"Lead Det.Jet/Hepg #gamma P_{T} Ratio for 1.3<#Delta#Phi^{(Lead Jet,#slash{E}_{T})}<1.7, after cuts");
	Hist.hDebug_PhoJetPtRatio_JMetdelphiMid= new TH1F(name,title,1000,-2.0,2.0);
	Hist.hDebug_PhoJetPtRatio_JMetdelphiMid->GetXaxis()->SetTitle("Lead DET Jet/HEPG #gamma - 1");
	AddHistogram(Hist.hDebug_PhoJetPtRatio_JMetdelphiMid,Folder);

	sprintf(name,"debug_PhoPtRatio_DelPhiJetMetTop");
	sprintf(title,"#gamma P_{T} Ratio for #Delta#Phi^{(Lead Jet,#slash{E}_{T})}>2.8, after cuts");
	Hist.hDebug_PtRatio_JMetdelphiTop= new TH1F(name,title,1000,-1.0,1.0);
	Hist.hDebug_PtRatio_JMetdelphiTop->GetXaxis()->SetTitle("Lead DET. #gamma/HEPG #gamma - 1");
	AddHistogram(Hist.hDebug_PtRatio_JMetdelphiTop,Folder);

	sprintf(name,"debug_PhoJetPtRatio_DelPhiJetMetTop");
	sprintf(title,"Lead Det.Jet/Hepg #gamma P_{T} Ratio for #Delta#Phi^{(Lead Jet,#slash{E}_{T})}>2.8, after cuts");
	Hist.hDebug_PhoJetPtRatio_JMetdelphiTop= new TH1F(name,title,1000,-2.0,2.0);
	Hist.hDebug_PhoJetPtRatio_JMetdelphiTop->GetXaxis()->SetTitle("Lead DET Jet/HEPG #gamma - 1");
	AddHistogram(Hist.hDebug_PhoJetPtRatio_JMetdelphiTop,Folder);


	sprintf(name,"debug_Jet1PtRatio_dPhi4");
	sprintf(title,"Lead Det.Jet/Had Jet P_{T} Ratio for #Delta#Phi^{(Lead Jet,#slash{E}_{T})}>2.8 & #Delta#Phi(Lead Jet, Sublead Jet) <0.4, after cuts");
	Hist.hDebug_Jet1PtRatio_dphi4= new TH1F(name,title,1000,-2.0,2.0);
	Hist.hDebug_Jet1PtRatio_dphi4->GetXaxis()->SetTitle("Lead DET Jet/HAD Jet - 1");
	AddHistogram(Hist.hDebug_Jet1PtRatio_dphi4,Folder);

	sprintf(name,"debug_Jet2PtRatio_dPhi4");
	sprintf(title,"Sublead Det.Jet/Had Jet P_{T} Ratio for #Delta#Phi^{(Lead Jet,#slash{E}_{T})}>2.8 & #Delta#Phi(Lead Jet, Sublead Jet)<0.4, after cuts");
	Hist.hDebug_Jet2PtRatio_dphi4= new TH1F(name,title,1000,-2.0,2.0);
	Hist.hDebug_Jet2PtRatio_dphi4->GetXaxis()->SetTitle("Sublead DET Jet/HAD Jet - 1");
	AddHistogram(Hist.hDebug_Jet2PtRatio_dphi4,Folder);

	sprintf(name,"debug_LeadJetEt");
	sprintf(title,"E_{T} of 1^{st} jet after removing only the lead photon, after cuts");
	Hist.hDebug_LeadJet = new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.hDebug_LeadJet,Folder);


	sprintf(name,"debug_PhoRes");
	sprintf(title,"Photon energy resolution, after cuts");
	Hist.hDebug_PhoRes= new TH1F(name,title,1000,-2.0,2.0);
	Hist.hDebug_PhoRes->GetXaxis()->SetTitle("Photon E Det. /HAD Photon E - 1");
	AddHistogram(Hist.hDebug_PhoRes,Folder);


	// 2D plots with MEtSig
	sprintf(name,"Met_MetSig");
	sprintf(title,"#slash{E}_{T} Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_Met_MetSig=new TH2F(name,title,200,0.0,20.0,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Met_MetSig,Folder);

	sprintf(name,"Njet15_MetSig");
	sprintf(title,"N_{jet}(E_{T}^{jet}>15 GeV) Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_Njet15_MetSig= new TH2F(name,title,200,0,20,21,-0.5,20.5);
	AddHistogram(Hist.fAna_Njet15_MetSig,Folder);

	sprintf(name,"Et1_MetSig");
	sprintf(title,"E_{T} of 1^{st} EM object Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_Et1_MetSig=new TH2F(name,title,200,0,20.0,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Et1_MetSig,Folder);

	
	sprintf(name,"Etjet_MetSig");
	sprintf(title,"E_{T} of 1^{st} jet Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_Etjet_MetSig=new TH2F(name,title,200,0,20,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Etjet_MetSig,Folder);

	sprintf(name,"Etjet2_MetSig");
	sprintf(title,"E_{T} of 2^{nd} jet Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_Etjet2_MetSig=new TH2F(name,title,200,0,20,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Etjet2_MetSig,Folder);

	sprintf(name,"Ht_MetSig");
	sprintf(title,"H_{T} Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_Ht_MetSig=new TH2F(name,title,200,0,20,1000,0.0,1000.0);
	AddHistogram(Hist.fAna_Ht_MetSig,Folder);

	sprintf(name,"dPhi1_MetSig");
	sprintf(title,"#Delta#phi(#slash{E}_{T}-em1) Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_dPhi1_MetSig= new TH2F(name,title,200,0,20,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi1_MetSig,Folder);

	sprintf(name,"dPhi3_MetSig");
	sprintf(title,"#Delta#phi(#slash{E}_{T}-jet1) if E_{T}^{jet}>15 GeV Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_dPhi3_MetSig= new TH2F(name,title,200,0,20,100,0.0,TMath::Pi());
	AddHistogram(Hist.fAna_dPhi3_MetSig,Folder);

	sprintf(name,"jet1Eta_MetSig");
	sprintf(title,"Lead Jet 4-vec #eta Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_Jet1Eta_MetSig= new TH2F(name,title,200,0,20,160,-4.0,4.0);
	AddHistogram(Hist.fAna_Jet1Eta_MetSig,Folder);

	sprintf(name,"jet2Eta_MetSig");
	sprintf(title,"Subleading Jet 4-vec #eta Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_Jet2Eta_MetSig= new TH2F(name,title,200,0,20,160,-4.0,4.0);
	AddHistogram(Hist.fAna_Jet2Eta_MetSig,Folder);

	sprintf(name,"pho1Eta_MetSig");
	sprintf(title,"Lead #gamma 4-vec #eta Vs #slash{E}_{T}Sig, after cuts");
	Hist.fAna_Pho1Eta_MetSig= new TH2F(name,title,200,0,20,160,-4.0,4.0);
	AddHistogram(Hist.fAna_Pho1Eta_MetSig,Folder);


	sprintf(name,"debug_offsetJetE");
	sprintf(title,"jer off_set ==0 cases by Jet E, after cuts");
	Hist.hDebug_offsetVsJetE =new TH3F(name,title,100,0.0,500.0, 25,0,5, 1000,0.0,1000.0);
	AddHistogram(Hist.hDebug_offsetVsJetE ,Folder);


	
	// end debugging info
	return;
}


//______________________ Booking general jet histograms
void JetFilterModuleV2::BookJetHistograms(JetGeneral_t& Hist, const char* Folder, const char* algoname) {
	char name [200];
	char title[200];
	// book histograms

	sprintf(name,"toyMET_all");
	sprintf(title,"#slash{E}_{T}^{toy}, all events");
	Hist.fEvnt_toyMET_all=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fEvnt_toyMET_all,Folder);
	sprintf(name,"toyMET_cut");
	sprintf(title,"#slash{E}_{T}^{toy}, events after Cleanup cuts");
	Hist.fEvnt_toyMET_cut=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fEvnt_toyMET_cut,Folder);
	sprintf(name,"corMET_all");
	sprintf(title,"#slash{E}_{T}^{cor}, all events before Cleanup cuts");
	Hist.fEvnt_corMET_all=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fEvnt_corMET_all,Folder);

	//__________________________________ before cuts

	sprintf(title,"%s%s",algoname,": Number of Jets");
	Hist.fEvntNjet_b= new TH1F("Njet_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>5 GeV");
	Hist.fEvntNjet5_b= new TH1F("Njet5_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet5_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>10 GeV");
	Hist.fEvntNjet10_b= new TH1F("Njet10_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet10_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>15 GeV");
	Hist.fEvntNjet15_b= new TH1F("Njet15_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet15_b,Folder);
	Hist.fEvntNjet20_b= new TH1F("Njet20_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet20_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>25 GeV");
	Hist.fEvntNjet25_b= new TH1F("Njet25_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet25_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>30 GeV");
	Hist.fEvntNjet30_b= new TH1F("Njet30_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet30_b,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>35 GeV");
	Hist.fEvntNjet35_b= new TH1F("Njet35_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet35_b,Folder);

	sprintf(title," %s%s",algoname,": dZ, jet & best vertex");
	Hist.fEvntdZ_b= new TH1F("dZ_b",title,400,-200.0,200.0);
	AddHistogram(Hist.fEvntdZ_b,Folder);
	sprintf(title," %s%s",algoname,": raw Et, 1st jet");
	Hist.fEvntEt0_b[0]= new TH1F("Et0_b[0]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_b[0],Folder);
	sprintf(title," %s%s",algoname,": raw Et, 2nd jet");
	Hist.fEvntEt0_b[1]= new TH1F("Et0_b[1]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_b[1],Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), 1st jet");
	Hist.fEvntEt_b[0]= new TH1F("Et_b[0]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_b[0],Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), 2nd jet");
	Hist.fEvntEt_b[1]= new TH1F("Et_b[1]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_b[1],Folder);
	sprintf(title," %s%s",algoname,": detector Eta, 1st jet");
	Hist.fEvntEtaDet_b[0]= new TH1F("EtaDet_b[0]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDet_b[0],Folder);
	sprintf(title," %s%s",algoname,": detector Eta, 2nd jet");
	Hist.fEvntEtaDet_b[1]= new TH1F("EtaDet_b[1]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDet_b[1],Folder);
	sprintf(title," %s%s",algoname,": Eta, 1st jet");
	Hist.fEvntEta_b[0]= new TH1F("Eta_b[0]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEta_b[0],Folder);
	sprintf(title," %s%s",algoname,": Eta, 2nd jet");
	Hist.fEvntEta_b[1]= new TH1F("Eta_b[1]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEta_b[1],Folder);
	sprintf(title," %s%s",algoname,": Phi, 1st jet");
	Hist.fEvntPhi_b[0]= new TH1F("Phi_b[0]",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhi_b[0],Folder);
	sprintf(title," %s%s",algoname,": Phi, 2nd jet");
	Hist.fEvntPhi_b[1]= new TH1F("Phi_b[1]",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhi_b[1],Folder);
	sprintf(title," %s%s",algoname,": Theta, 1st jet");
	Hist.fEvntTheta_b[0]= new TH1F("Theta_b[0]",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntTheta_b[0],Folder);
	sprintf(title," %s%s",algoname,": Theta, 2nd jet");
	Hist.fEvntTheta_b[1]= new TH1F("Theta_b[1]",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntTheta_b[1],Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, 1st jet");
	Hist.fEvntEmFr_b[0]= new TH1F("EmFr_b[0]",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFr_b[0],Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, 2nd jet");
	Hist.fEvntEmFr_b[1]= new TH1F("EmFr_b[1]",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFr_b[1],Folder);
	sprintf(title," %s%s",algoname,": Number of towers, 1st jet");
	Hist.fEvntNTowers_b[0]= new TH1F("NTowers_b[0]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowers_b[0],Folder);
	sprintf(title," %s%s",algoname,": Number of towers, 2nd jet");
	Hist.fEvntNTowers_b[1]= new TH1F("NTowers_b[1]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowers_b[1],Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, 1st jet");
	Hist.fEvntNTracks_b[0]= new TH1F("NTracks_b[0]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracks_b[0],Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, 2nd jet");
	Hist.fEvntNTracks_b[1]= new TH1F("NTracks_b[1]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracks_b[1],Folder);
	sprintf(title," %s%s",algoname,": raw Et, extra jets");
	Hist.fEvntEt0X_b= new TH1F("Et0X_b",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0X_b,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), extra jets");
	Hist.fEvntEtX_b= new TH1F("EtX_b",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEtX_b,Folder);
	sprintf(title," %s%s",algoname,": detector Eta, extra jets");
	Hist.fEvntEtaDetX_b= new TH1F("EtaDetX_b",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDetX_b,Folder);
	sprintf(title," %s%s",algoname,": Eta, extra jets");
	Hist.fEvntEtaX_b= new TH1F("EtaX_b",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaX_b,Folder);
	sprintf(title," %s%s",algoname,": Phi, extra jets");
	Hist.fEvntPhiX_b= new TH1F("PhiX_b",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhiX_b,Folder);
	sprintf(title," %s%s",algoname,": Theta, extra jets");
	Hist.fEvntThetaX_b= new TH1F("ThetaX_b",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntThetaX_b,Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, extra jets");
	Hist.fEvntEmFrX_b= new TH1F("EmFrX_b",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFrX_b,Folder);
	sprintf(title," %s%s",algoname,": Number of towers, extra jets");
	Hist.fEvntNTowersX_b= new TH1F("NTowersX_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowersX_b,Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, extra jets");
	Hist.fEvntNTracksX_b= new TH1F("NTracksX_b",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracksX_b,Folder);
	sprintf(title," %s%s",algoname,": My Invariant Mass, 2 leading jets");
	Hist.fEvntMjj_b= new TH1F("Mjj_b",title,2000,0.0,2000.0);
	AddHistogram(Hist.fEvntMjj_b,Folder);
	sprintf(title," %s%s",algoname,": Theta* of 2 leading jets");
	Hist.fEvntThetaStar_b= new TH1F("ThetaStar_b",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntThetaStar_b,Folder);
	sprintf(title," %s%s",algoname,": Delta Phi of 2 leading jets");
	Hist.fEvntDeltaPhi_b= new TH1F("DeltaPhi_b",title,350,0.0,3.5);
	AddHistogram(Hist.fEvntDeltaPhi_b,Folder);
	sprintf(title," %s%s",algoname,": Delta Eta of 2 leading jets");
	Hist.fEvntDeltaEta_b= new TH1F("DeltaEta_b",title,800,-8.0,8.0);
	AddHistogram(Hist.fEvntDeltaEta_b,Folder);
	sprintf(title," %s%s",algoname,": Delta R of 2 leading jets");
	Hist.fEvntDeltaR_b= new TH1F("DeltaR_b",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR_b,Folder);
	sprintf(title," %s%s",algoname,": dR, 1st jet & any extra jet");
	Hist.fEvntDeltaR1x_b= new TH1F("DeltaR1x_b",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR1x_b,Folder);
	sprintf(title," %s%s",algoname,": dR, 2nd jet & any extra jet");
	Hist.fEvntDeltaR2x_b= new TH1F("DeltaR2x_b",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR2x_b,Folder);
	sprintf(title," %s%s",algoname,": Kt of two best jets");
	Hist.fEvntKt2jet_b= new TH1F("Kt2jet_b",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntKt2jet_b,Folder);
	sprintf(title," %s%s",algoname,": Kt of all jets");
	Hist.fEvntKtAll_b= new TH1F("KtAll_b",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntKtAll_b,Folder);

	sprintf(title," %s%s",algoname,": Njet vs. Nvx(class>=12)");
	Hist.fEvntNJet_Nvx_b= new TProfile("NJet_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet5 vs. Nvx(class>=12)");
	Hist.fEvntNJet5_Nvx_b= new TProfile("NJet5_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet5_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet10 vs. Nvx(class>=12)");
	Hist.fEvntNJet10_Nvx_b= new TProfile("NJet10_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet10_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet15 vs. Nvx(class>=12)");
	Hist.fEvntNJet15_Nvx_b= new TProfile("NJet15_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet15_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet20 vs. Nvx(class>=12)");
	Hist.fEvntNJet20_Nvx_b= new TProfile("NJet20_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet20_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet25 vs. Nvx(class>=12)");
	Hist.fEvntNJet25_Nvx_b= new TProfile("NJet25_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet25_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet30 vs. Nvx(class>=12)");
	Hist.fEvntNJet30_Nvx_b= new TProfile("NJet30_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet30_Nvx_b,Folder);
	sprintf(title," %s%s",algoname,": Njet35 vs. Nvx(class>=12)");
	Hist.fEvntNJet35_Nvx_b= new TProfile("NJet35_Nvx_b",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet35_Nvx_b,Folder);

	sprintf(title," %s%s",algoname,": 2 leading jets, Kt vs. Mjj");
	Hist.fEvntKt_Mjj_b= new TProfile("Kt_Mjj_b",title,400,0.0,1000.0,0.0,1000.0);
	AddHistogram(Hist.fEvntKt_Mjj_b,Folder);
	sprintf(title," %s%s",algoname,": KtAll vs. Mjj");
	Hist.fEvntKtAll_Mjj_b= new TProfile("KtAll_Mjj_b",title,400,0.0,1000.0,0.0,1000.0);
	AddHistogram(Hist.fEvntKtAll_Mjj_b,Folder);
	sprintf(title," %s%s",algoname,": raw Et vs. Njet");
	Hist.fEvntEt0_Njet_b= new TProfile("Et0_Njet_b",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_Njet_b,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6) vs. Njet");
	Hist.fEvntEt_Njet_b= new TProfile("Et_Njet_b",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_Njet_b,Folder);
	sprintf(title," %s%s",algoname,": raw Et vs. Nvx(class>=12)");
	Hist.fEvntEt0_Nvx12_b= new TProfile("Et0_Nvx12_b",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_Nvx12_b,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6) vs. Nvx(class>=12)");
	Hist.fEvntEt_Nvx12_b= new TProfile("Et_Nvx12_b",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_Nvx12_b,Folder);
	sprintf(title," %s%s",algoname,": Njet vs. Lum for events with Nvx12=1");
	Hist.fEvntNjet_Lum_b= new TProfile("Njet_Lum_b",title,50,0.0,100.0,0.0,1000.0);
	AddHistogram(Hist.fEvntNjet_Lum_b,Folder);


	//________________________________________________________________ after cuts
	sprintf(title,"%s%s",algoname,": Number of Jets");
	Hist.fEvntNjet_a= new TH1F("Njet_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>5 GeV");
	Hist.fEvntNjet5_a= new TH1F("Njet5_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet5_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>10 GeV");
	Hist.fEvntNjet10_a= new TH1F("Njet10_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet10_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>15 GeV");
	Hist.fEvntNjet15_a= new TH1F("Njet15_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet15_a,Folder);
	Hist.fEvntNjet20_a= new TH1F("Njet20_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet20_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>25 GeV");
	Hist.fEvntNjet25_a= new TH1F("Njet25_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet25_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>30 GeV");
	Hist.fEvntNjet30_a= new TH1F("Njet30_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet30_a,Folder);
	sprintf(title,"%s%s",algoname,": Number of Jets with E_{T}>35 GeV");
	Hist.fEvntNjet35_a= new TH1F("Njet35_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNjet35_a,Folder);

	sprintf(title," %s%s",algoname,": dZ, jet & best vertex");
	Hist.fEvntdZ_a= new TH1F("dZ_a",title,400,-200.0,200.0);
	AddHistogram(Hist.fEvntdZ_a,Folder);
	sprintf(title," %s%s",algoname,": raw Et, 1st jet");
	Hist.fEvntEt0_a[0]= new TH1F("Et0_a[0]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_a[0],Folder);
	sprintf(title," %s%s",algoname,": raw Et, 2nd jet");
	Hist.fEvntEt0_a[1]= new TH1F("Et0_a[1]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_a[1],Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), 1st jet");
	Hist.fEvntEt_a[0]= new TH1F("Et_a[0]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_a[0],Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), 2nd jet");
	Hist.fEvntEt_a[1]= new TH1F("Et_a[1]",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_a[1],Folder);
	sprintf(title," %s%s",algoname,": detector Eta, 1st jet");
	Hist.fEvntEtaDet_a[0]= new TH1F("EtaDet_a[0]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDet_a[0],Folder);
	sprintf(title," %s%s",algoname,": detector Eta, 2nd jet");
	Hist.fEvntEtaDet_a[1]= new TH1F("EtaDet_a[1]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDet_a[1],Folder);
	sprintf(title," %s%s",algoname,": Eta, 1st jet");
	Hist.fEvntEta_a[0]= new TH1F("Eta_a[0]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEta_a[0],Folder);
	sprintf(title," %s%s",algoname,": Eta, 2nd jet");
	Hist.fEvntEta_a[1]= new TH1F("Eta_a[1]",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEta_a[1],Folder);
	sprintf(title," %s%s",algoname,": Phi, 1st jet");
	Hist.fEvntPhi_a[0]= new TH1F("Phi_a[0]",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhi_a[0],Folder);
	sprintf(title," %s%s",algoname,": Phi, 2nd jet");
	Hist.fEvntPhi_a[1]= new TH1F("Phi_a[1]",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhi_a[1],Folder);
	sprintf(title," %s%s",algoname,": Theta, 1st jet");
	Hist.fEvntTheta_a[0]= new TH1F("Theta_a[0]",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntTheta_a[0],Folder);
	sprintf(title," %s%s",algoname,": Theta, 2nd jet");
	Hist.fEvntTheta_a[1]= new TH1F("Theta_a[1]",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntTheta_a[1],Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, 1st jet");
	Hist.fEvntEmFr_a[0]= new TH1F("EmFr_a[0]",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFr_a[0],Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, 2nd jet");
	Hist.fEvntEmFr_a[1]= new TH1F("EmFr_a[1]",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFr_a[1],Folder);
	sprintf(title," %s%s",algoname,": Number of towers, 1st jet");
	Hist.fEvntNTowers_a[0]= new TH1F("NTowers_a[0]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowers_a[0],Folder);
	sprintf(title," %s%s",algoname,": Number of towers, 2nd jet");
	Hist.fEvntNTowers_a[1]= new TH1F("NTowers_a[1]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowers_a[1],Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, 1st jet");
	Hist.fEvntNTracks_a[0]= new TH1F("NTracks_a[0]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracks_a[0],Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, 2nd jet");
	Hist.fEvntNTracks_a[1]= new TH1F("NTracks_a[1]",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracks_a[1],Folder);
	sprintf(title," %s%s",algoname,": raw Et, extra jets");
	Hist.fEvntEt0X_a= new TH1F("Et0X_a",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0X_a,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6), extra jets");
	Hist.fEvntEtX_a= new TH1F("EtX_a",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntEtX_a,Folder);
	sprintf(title," %s%s",algoname,": detector Eta, extra jets");
	Hist.fEvntEtaDetX_a= new TH1F("EtaDetX_a",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaDetX_a,Folder);
	sprintf(title," %s%s",algoname,": Eta, extra jets");
	Hist.fEvntEtaX_a= new TH1F("EtaX_a",title,200,-5.0,5.0);
	AddHistogram(Hist.fEvntEtaX_a,Folder);
	sprintf(title," %s%s",algoname,": Phi, extra jets");
	Hist.fEvntPhiX_a= new TH1F("PhiX_a",title,700,-0.5,6.5);
	AddHistogram(Hist.fEvntPhiX_a,Folder);
	sprintf(title," %s%s",algoname,": Theta, extra jets");
	Hist.fEvntThetaX_a= new TH1F("ThetaX_a",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntThetaX_a,Folder);
	sprintf(title," %s%s",algoname,": Em Fraction, extra jets");
	Hist.fEvntEmFrX_a= new TH1F("EmFrX_a",title,200,0.0,2.0);
	AddHistogram(Hist.fEvntEmFrX_a,Folder);
	sprintf(title," %s%s",algoname,": Number of towers, extra jets");
	Hist.fEvntNTowersX_a= new TH1F("NTowersX_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTowersX_a,Folder);
	sprintf(title," %s%s",algoname,": Number of tracks, extra jets");
	Hist.fEvntNTracksX_a= new TH1F("NTracksX_a",title,100,-0.5,99.5);
	AddHistogram(Hist.fEvntNTracksX_a,Folder);
	sprintf(title," %s%s",algoname,": My Invariant Mass, 2 leading jets");
	Hist.fEvntMjj_a= new TH1F("Mjj_a",title,2000,0.0,2000.0);
	AddHistogram(Hist.fEvntMjj_a,Folder);
	sprintf(title," %s%s",algoname,": Theta* of 2 leading jets");
	Hist.fEvntThetaStar_a= new TH1F("ThetaStar_a",title,400,-0.5,3.5);
	AddHistogram(Hist.fEvntThetaStar_a,Folder);
	sprintf(title," %s%s",algoname,": Delta Phi of 2 leading jets");
	Hist.fEvntDeltaPhi_a= new TH1F("DeltaPhi_a",title,350,0.0,3.5);
	AddHistogram(Hist.fEvntDeltaPhi_a,Folder);
	sprintf(title," %s%s",algoname,": Delta Eta of 2 leading jets");
	Hist.fEvntDeltaEta_a= new TH1F("DeltaEta_a",title,800,-8.0,8.0);
	AddHistogram(Hist.fEvntDeltaEta_a,Folder);
	sprintf(title," %s%s",algoname,": Delta R of 2 leading jets");
	Hist.fEvntDeltaR_a= new TH1F("DeltaR_a",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR_a,Folder);
	sprintf(title," %s%s",algoname,": dR, 1st jet & any extra jet");
	Hist.fEvntDeltaR1x_a= new TH1F("DeltaR1x_a",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR1x_a,Folder);
	sprintf(title," %s%s",algoname,": dR, 2nd jet & any extra jet");
	Hist.fEvntDeltaR2x_a= new TH1F("DeltaR2x_a",title,800,0.0,8.0);
	AddHistogram(Hist.fEvntDeltaR2x_a,Folder);
	sprintf(title," %s%s",algoname,": Kt of two best jets");
	Hist.fEvntKt2jet_a= new TH1F("Kt2jet_a",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntKt2jet_a,Folder);
	sprintf(title," %s%s",algoname,": Kt of all jets");
	Hist.fEvntKtAll_a= new TH1F("KtAll_a",title,1000,0.0,1000.0);
	AddHistogram(Hist.fEvntKtAll_a,Folder);

	sprintf(title," %s%s",algoname,": Njet vs. Nvx(class>=12)");
	Hist.fEvntNJet_Nvx_a= new TProfile("NJet_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet5 vs. Nvx(class>=12)");
	Hist.fEvntNJet5_Nvx_a= new TProfile("NJet5_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet5_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet10 vs. Nvx(class>=12)");
	Hist.fEvntNJet10_Nvx_a= new TProfile("NJet10_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet10_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet15 vs. Nvx(class>=12)");
	Hist.fEvntNJet15_Nvx_a= new TProfile("NJet15_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet15_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet20 vs. Nvx(class>=12)");
	Hist.fEvntNJet20_Nvx_a= new TProfile("NJet20_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet20_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet25 vs. Nvx(class>=12)");
	Hist.fEvntNJet25_Nvx_a= new TProfile("NJet25_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet25_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet30 vs. Nvx(class>=12)");
	Hist.fEvntNJet30_Nvx_a= new TProfile("NJet30_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet30_Nvx_a,Folder);
	sprintf(title," %s%s",algoname,": Njet35 vs. Nvx(class>=12)");
	Hist.fEvntNJet35_Nvx_a= new TProfile("NJet35_Nvx_a",title,50,-0.5,49.5,-0.5,50.0);
	AddHistogram(Hist.fEvntNJet35_Nvx_a,Folder);

	sprintf(title," %s%s",algoname,": 2 leading jets, Kt vs. Mjj");
	Hist.fEvntKt_Mjj_a= new TProfile("Kt_Mjj_a",title,400,0.0,1000.0,0.0,1000.0);
	AddHistogram(Hist.fEvntKt_Mjj_a,Folder);
	sprintf(title," %s%s",algoname,": KtAll vs. Mjj");
	Hist.fEvntKtAll_Mjj_a= new TProfile("KtAll_Mjj_a",title,400,0.0,1000.0,0.0,1000.0);
	AddHistogram(Hist.fEvntKtAll_Mjj_a,Folder);
	sprintf(title," %s%s",algoname,": raw Et vs. Njet");
	Hist.fEvntEt0_Njet_a= new TProfile("Et0_Njet_a",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_Njet_a,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6) vs. Njet");
	Hist.fEvntEt_Njet_a= new TProfile("Et_Njet_a",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_Njet_a,Folder);
	sprintf(title," %s%s",algoname,": raw Et vs. Nvx(class>=12)");
	Hist.fEvntEt0_Nvx12_a= new TProfile("Et0_Nvx12_a",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt0_Nvx12_a,Folder);
	sprintf(title," %s%s",algoname,": Et(lev6) vs. Nvx(class>=12)");
	Hist.fEvntEt_Nvx12_a= new TProfile("Et_Nvx12_a",title,50,-0.5,49.5,0.0,1000.0);
	AddHistogram(Hist.fEvntEt_Nvx12_a,Folder);
	sprintf(title," %s%s",algoname,": Njet vs. Lum for events with Nvx12=1");
	Hist.fEvntNjet_Lum_a= new TProfile("Njet_Lum_a",title,50,0.0,100.0,0.0,1000.0);
	AddHistogram(Hist.fEvntNjet_Lum_a,Folder);

	return;
}


//_________________________________________________ booking met histo
void JetFilterModuleV2::BookMetStudyHistograms(MetStudyHisto_t& Hist, const char* Folder, const char* algoname) {

	char name [200];
	char title[200];

	//------ upper limit on MetSignificance for real MET
	for(int i=0; i<5; i++)
	{
		char range[200];
		if(i==0) sprintf(range,"closest object");
		if(i==1) sprintf(range,"closest EM object");
		if(i==2) sprintf(range,"closest jet");
		if(i==3) sprintf(range,"1^{st} EM object");
		if(i==4) sprintf(range,"1^{st} jet");
		sprintf(name,"fMetSig_vs_dPhi[%i]",i);
		sprintf(title,"#slash{E}_{T}-significance vs. #Delta#phi, %s",range);
		Hist.fMetSig_vs_dPhi[i]= new TH2F(name,title,100,0.0,20,60,0.0,TMath::Pi());
		AddHistogram(Hist.fMetSig_vs_dPhi[i],Folder);
		sprintf(name,"fMetSigCalib_vs_dPhi[%i]",i);
		sprintf(title,"Pseudo-experiments: #slash{E}_{T}-significance vs. #Delta#phi, %s",range);
		Hist.fMetSigCalib_vs_dPhi[i]= new TH2F(name,title,100,0.0,20,60,0.0,TMath::Pi());
		AddHistogram(Hist.fMetSigCalib_vs_dPhi[i],Folder);
	}

	sprintf(name,"fMetSig_vs_sqrtHt");
	sprintf(title,"#slash{E}_{T}-significance vs. #sqrt{H_{T}}");
	Hist.fMetSig_vs_sqrtHt= new TH2F(name,title,100,0.0,20,100,0.0,30.0);
	AddHistogram(Hist.fMetSig_vs_sqrtHt,Folder);
	sprintf(name,"fMetSig_vs_sqrtMet");
	sprintf(title,"#slash{E}_{T}-significance vs. #sqrt{#slash{E}_{T}}");
	Hist.fMetSig_vs_sqrtMet= new TH2F(name,title,100,0.0,20,100,0.0,30.0);
	AddHistogram(Hist.fMetSig_vs_sqrtMet,Folder);
	sprintf(name,"fMetSig_vs_sqrtSumEt");
	sprintf(title,"#slash{E}_{T}-significance vs. #sqrt{#SigmaE_{T}}");
	Hist.fMetSig_vs_sqrtSumEt= new TH2F(name,title,100,0.0,20,100,0.0,30.0);
	AddHistogram(Hist.fMetSig_vs_sqrtSumEt,Folder);

	sprintf(name,"fMetSig_vs_Xces");
	sprintf(title,"#slash{E}_{T}-significance vs. X_{CES} for #gamma/e if #Delta#phi_{#slash{E}_{T}-#gamma/e}<0.4");
	Hist.fMetSig_vs_Xces= new TH2F(name,title,100,0.0,20,100,0.0,25.0);
	AddHistogram(Hist.fMetSig_vs_Xces,Folder);
	sprintf(name,"fMetSig_vs_etaEM");
	sprintf(title,"#slash{E}_{T}-significance vs. |#eta_{det}| for #gamma/e if #Delta#phi_{#slash{E}_{T}-#gamma/e}<0.4");
	Hist.fMetSig_vs_etaEM= new TH2F(name,title,100,0.0,20,150,0.0,3.0);
	AddHistogram(Hist.fMetSig_vs_etaEM,Folder);
	sprintf(name,"fMetSig_Xces_vs_etaEM");
	sprintf(title,"|X_{ces}| vs. |#eta_{rel}| for #gamma/e if siginificance>4.0 and #Delta#phi_{#slash{E}_{T}-#gamma/e}<0.4");
	Hist.fMetSig_Xces_vs_etaEM= new TH2F(name,title,100,0.0,25.0,100,0.0,1.0);
	AddHistogram(Hist.fMetSig_Xces_vs_etaEM,Folder);
	sprintf(name,"fMetSig_vs_etaJET");
	sprintf(title,"#slash{E}_{T}-significance vs. |#eta_{det}| for jet if #Delta#phi_{#slash{E}_{T}-jet}<0.4");
	Hist.fMetSig_vs_etaJET= new TH2F(name,title,100,0.0,20,150,0.0,3.0);
	AddHistogram(Hist.fMetSig_vs_etaJET,Folder);
	sprintf(name,"ObjType");
	sprintf(title,"Type of closest to #slash{E}_{T} object (0=e/#gamma, 1=jet, 2=b-tag, 3=#mu, 4=#tau, 5=trk) for events MetSig>4.0");
	Hist.fObjType= new TH1F(name,title,6,-0.5,5.5);
	AddHistogram(Hist.fObjType,Folder);

	sprintf(name,"fMetSig_corrFactor");
	sprintf(title,"Correction factor for estimated #slash{E}_{T}-significance");
	Hist.fMetSig_corrFactor=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSig_corrFactor,Folder);

	sprintf(name,"fMetSig_estimate");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance");
	Hist.fMetSig_estimate=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSig_estimate,Folder);
	sprintf(name,"fMetSig_estimate_njet0");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance if N_{jet15}=0");
	Hist.fMetSig_estimate_njet0=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSig_estimate_njet0,Folder);
	sprintf(name,"fMetSig_estimate_njet1");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance if N_{jet15}=1");
	Hist.fMetSig_estimate_njet1=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSig_estimate_njet1,Folder);

	sprintf(name,"fMetSig_estimate_njet2");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance if N_{jet15}=2");
	Hist.fMetSig_estimate_njet2=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSig_estimate_njet2,Folder);

	sprintf(name,"fMetSig_estimate_njet3");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance if N_{jet15}=3");
	Hist.fMetSig_estimate_njet3=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSig_estimate_njet3,Folder);

	sprintf(name,"fMetSig_estimate_njet4");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance if N_{jet15}>=4");
	Hist.fMetSig_estimate_njet4=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSig_estimate_njet4,Folder);

	//------ upper limit on MetSignificance for MET from pseudo-experiments
	sprintf(name,"fMetSigCalib_estimate");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance for pseudo-experiments");
	Hist.fMetSigCalib_estimate=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSigCalib_estimate,Folder);
	sprintf(name,"fMetSigCalib_estimate_njet0");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance for pseudo-experiments with N_{jet15}=0");
	Hist.fMetSigCalib_estimate_njet0=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSigCalib_estimate_njet0,Folder);
	sprintf(name,"fMetSigCalib_estimate_njet1");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance for pseudo-experiments with N_{jet15}=1");
	Hist.fMetSigCalib_estimate_njet1=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSigCalib_estimate_njet1,Folder);

	sprintf(name,"fMetSigCalib_estimate_njet2");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance for pseudo-experiments with N_{jet15}=2");
	Hist.fMetSigCalib_estimate_njet2=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSigCalib_estimate_njet2,Folder);

	sprintf(name,"fMetSigCalib_estimate_njet3");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance for pseudo-experiments with N_{jet15}=3");
	Hist.fMetSigCalib_estimate_njet3=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSigCalib_estimate_njet3,Folder);

	sprintf(name,"fMetSigCalib_estimate_njet4");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance for pseudo-experiments with N_{jet15}>=4");
	Hist.fMetSigCalib_estimate_njet4=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSigCalib_estimate_njet4,Folder);

	//------ upper limit on MetSignificance for toy MET
	sprintf(name,"fMetSigToy_estimate");
	sprintf(title,"Estimated upper limit on #slash{E}_{T}-significance for toy MC");
	Hist.fMetSigToy_estimate=new TH1F(name,title,1000,0.0,20.0);
	AddHistogram(Hist.fMetSigToy_estimate,Folder);

	for(int i=0; i<5; i++)
	{
		sprintf(name,"fMetSig_Met[%i]",i);
		sprintf(title,"#slash{E}_{T} for events with Significance>%i",i);
		Hist.fMetSig_Met[i]=new TH1F(name,title,1000,0.0,1000.0);
		AddHistogram(Hist.fMetSig_Met[i],Folder);

		sprintf(name,"fMetSigCalib_Met[%i]",i);
		sprintf(title,"pseudo-experiments: #slash{E}_{T} for events with Significance>%i",i);
		Hist.fMetSigCalib_Met[i]=new TH1F(name,title,1000,0.0,1000.0);
		AddHistogram(Hist.fMetSigCalib_Met[i],Folder);

		sprintf(name,"fMetSigToy_Met[%i]",i);
		if(i==0) sprintf(title,"toy MC: #slash{E}_{T} for events with Significance>%i",i);
		else sprintf(title,"toy MC: #slash{E}_{T} for events with Significance>%i",i+2);
		Hist.fMetSigToy_Met[i]=new TH1F(name,title,1000,0.0,1000.0);
		AddHistogram(Hist.fMetSigToy_Met[i],Folder);

		if(i<4)
		{
			sprintf(name,"fMetSigToy_eff[%i]",i);
			sprintf(title,"toy MC: efficiency of MetSig>%i cut",i+3);
			Hist.fMetSigToy_eff[i]=new TH1F(name,title,1000,0.0,1000.0);
			AddHistogram(Hist.fMetSigToy_eff[i],Folder);
			Hist.fMetSigToy_eff[i]->Sumw2();
		}
	}

	sprintf(name,"SumEtCorr");
	sprintf(title," %s%s",algoname,": corrected #sum E_{T}");
	Hist.fSumEtCorr=new TH1F(name,title,2000,0.0,2000.0);
	AddHistogram(Hist.fSumEtCorr,Folder);
	sprintf(name,"SqrtSumEtCorr");
	sprintf(title," %s%s",algoname,": corrected #sqrt{#sum E_{T}}");
	Hist.fSqrtSumEtCorr=new TH1F(name,title,180,0.0,45.0);
	AddHistogram(Hist.fSqrtSumEtCorr,Folder);

	sprintf(name,"SumEtCorr_noJet_vx1");
	sprintf(title," %s%s",algoname,": corrected #sum E_{T}, N_{jet15}=0, N_{vx}=1");
	Hist.fSumEtCorr_noJet_vx1=new TH1F(name,title,2000,0.0,2000.0);
	AddHistogram(Hist.fSumEtCorr_noJet_vx1,Folder);
	sprintf(name,"SumEtCorr_noJet_vx2");
	sprintf(title," %s%s",algoname,": corrected #sum E_{T}, N_{jet15}=0, N_{vx}=2");
	Hist.fSumEtCorr_noJet_vx2=new TH1F(name,title,2000,0.0,2000.0);
	AddHistogram(Hist.fSumEtCorr_noJet_vx2,Folder);
	sprintf(name,"SumEtCorr_noJet_vx3");
	sprintf(title," %s%s",algoname,": corrected #sum E_{T}, N_{jet15}=0, N_{vx}=3");
	Hist.fSumEtCorr_noJet_vx3=new TH1F(name,title,2000,0.0,2000.0);
	AddHistogram(Hist.fSumEtCorr_noJet_vx3,Folder);
	sprintf(name,"SumEtCorrNoJet_vs_Nvx");
	sprintf(title," %s%s",algoname,": corrected #sum E_{T} .vs. N_{vx} for N_{jet15}=0 events");
	Hist.fSumEtCorrNoJet_vs_Nvx=new TProfile(name,title,11,-0.5,10.5,0.0,2000.0);
	AddHistogram(Hist.fSumEtCorrNoJet_vs_Nvx,Folder);


	sprintf(name,"SumEtCorr_withJet");
	sprintf(title," %s%s",algoname,": corrected #sum E_{T} for events with N_{jet}(E_{T}>Threshold)>0");
	Hist.fSumEtCorr_withJet=new TH1F(name,title,2000,0.0,2000.0);
	AddHistogram(Hist.fSumEtCorr_withJet,Folder);
	sprintf(name,"SqrtSumEtCorr_withJet");
	sprintf(title," %s%s",algoname,": corrected #sqrt{#sum E_{T}} for events with N_{jet}(E_{T}>Threshold)>0");
	Hist.fSqrtSumEtCorr_withJet=new TH1F(name,title,180,0.0,45.0);
	AddHistogram(Hist.fSqrtSumEtCorr_withJet,Folder);
	sprintf(name,"SumEtCorr_noJet");
	sprintf(title," %s%s",algoname,": corrected #sum E_{T} for events with N_{jet}(E_{T}>Threshold)=0");
	Hist.fSumEtCorr_noJet=new TH1F(name,title,2000,0.0,2000.0);
	AddHistogram(Hist.fSumEtCorr_noJet,Folder);
	sprintf(name,"SqrtSumEtCorr_noJet");
	sprintf(title," %s%s",algoname,": corrected #sqrt{#sum E_{T}} for events with N_{jet}(E_{T}>Threshold)=0");
	Hist.fSqrtSumEtCorr_noJet=new TH1F(name,title,180,0.0,45.0);
	AddHistogram(Hist.fSqrtSumEtCorr_noJet,Folder);

	sprintf(name,"MetRaw");
	sprintf(title," %s%s",algoname,": MEt_4 before corr. for jets (corrected events)");
	Hist.fMetRaw=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fMetRaw,Folder);
	sprintf(name,"MetCorr");
	sprintf(title," %s%s",algoname,": MEt_4 after corr. for jets (corrected events)");
	Hist.fMetCorr=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fMetCorr,Folder);
	sprintf(name,"MetPhiRaw");
	sprintf(title," %s%s",algoname,": #phi of MEt_4 before corr. for jets (corrected events)");
	Hist.fMetPhiRaw=new TH1F(name,title,128,0.0,6.4);
	AddHistogram(Hist.fMetPhiRaw,Folder);
	sprintf(name,"MetPhiCorr");
	sprintf(title," %s%s",algoname,": #phi of MEt_4 after corr. for jets (corrected events)");
	Hist.fMetPhiCorr=new TH1F(name,title,128,0.0,6.4);
	AddHistogram(Hist.fMetPhiCorr,Folder);
	sprintf(name,"MetRawAll");
	sprintf(title," %s%s",algoname,": MEt_4 before corr. for jets (all events)" );
	Hist.fMetRawAll=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fMetRawAll,Folder);
	sprintf(name,"MetCorrAll");
	sprintf(title," %s%s",algoname,": MEt_4 after corr. for jets (all events)");
	Hist.fMetCorrAll=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fMetCorrAll,Folder);
	sprintf(name,"MetPhiRawAll");
	sprintf(title," %s%s",algoname,": #phi of MEt_4 before corr. for jets (all events)");
	Hist.fMetPhiRawAll=new TH1F(name,title,128,0.0,6.4);
	AddHistogram(Hist.fMetPhiRawAll,Folder);
	sprintf(name,"MetPhiCorrAll");
	sprintf(title," %s%s",algoname,": #phi of MEt_4 after corr. for jets (all events)");
	Hist.fMetPhiCorrAll=new TH1F(name,title,128,0.0,6.4);
	AddHistogram(Hist.fMetPhiCorrAll,Folder);
	sprintf(name,"MetRawAll_X");
	sprintf(title," %s%s",algoname,": MEtX_4 before corr. for jets (all events)" );
	Hist.fMetRawAll_X=new TH1F(name,title,400,-200.0,200.0);
	AddHistogram(Hist.fMetRawAll_X,Folder);
	sprintf(name,"MetCorrAll_X");
	sprintf(title," %s%s",algoname,": MEtX_4 after corr. for jets (all events)");
	Hist.fMetCorrAll_X=new TH1F(name,title,400,-200.0,200.0);
	AddHistogram(Hist.fMetCorrAll_X,Folder);
	sprintf(name,"MetRawAll_Y");
	sprintf(title," %s%s",algoname,": MEtY_4 before corr. for jets (all events)" );
	Hist.fMetRawAll_Y=new TH1F(name,title,400,-200.0,200.0);
	AddHistogram(Hist.fMetRawAll_Y,Folder);
	sprintf(name,"MetCorrAll_Y");
	sprintf(title," %s%s",algoname,": MEtY_4 after corr. for jets (all events)");
	Hist.fMetCorrAll_Y=new TH1F(name,title,400,-200.0,200.0);
	AddHistogram(Hist.fMetCorrAll_Y,Folder);

	sprintf(name,"Met4VsNjet");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs N_{jet}");
	Hist.fMet4VsNjet=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsNjet,Folder);
	sprintf(name,"Met4VsNjet_th5");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs N_{jet} with E_{T}>5 GeV");
	Hist.fMet4VsNjet_th5=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsNjet_th5,Folder);
	sprintf(name,"Met4VsNjet_th10");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs N_{jet} with E_{T}>10 GeV");
	Hist.fMet4VsNjet_th10=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsNjet_th10,Folder);
	sprintf(name,"Met4VsNjet_th15");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs N_{jet} with E_{T}>15 GeV");
	Hist.fMet4VsNjet_th15=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsNjet_th15,Folder);
	sprintf(name,"Met4VsNjet_th20");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs N_{jet} with E_{T}>20 GeV");
	Hist.fMet4VsNjet_th20=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsNjet_th20,Folder);
	sprintf(name,"Met4VsNjet_th25");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs N_{jet} with E_{T}>25 GeV");
	Hist.fMet4VsNjet_th25=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsNjet_th25,Folder);
	sprintf(name,"Met4VsNjet_th30");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs N_{jet} with E_{T}>30 GeV");
	Hist.fMet4VsNjet_th30=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsNjet_th30,Folder);
	sprintf(name,"Met4VsNjet_th35");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs N_{jet} with E_{T}>35 GeV");
	Hist.fMet4VsNjet_th35=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsNjet_th35,Folder);

	sprintf(name,"Met4VsRun");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} Vs. Run for N_{vx12}=1");
	Hist.fMet4VsRun=new TProfile(name,title,100,130000.0,230000.0,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsRun,Folder);
	sprintf(name,"Nvx12VsRun");
	sprintf(title,"N_{vx12} Vs. Run");
	Hist.fNvx12VsRun=new TProfile(name,title,200,130000.0,230000.0,-1.0,100.0);
	AddHistogram(Hist.fNvx12VsRun,Folder);
	sprintf(name,"Met4VsNvx12");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs N_{vx12}");
	Hist.fMet4VsNvx12=new TProfile(name,title,10,-0.5,9.5,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsNvx12,Folder);

	sprintf(name,"SumEtJetFrac");
	sprintf(title," %s%s",algoname,": #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T}");
	Hist.fSumEtJetFrac=new TH1F(name,title,100,0.0,1.0);
	AddHistogram(Hist.fSumEtJetFrac,Folder);
	sprintf(name,"SumEtJet");
	sprintf(title," %s%s",algoname,": corrected #sum E^{jets}_{T}");
	Hist.fSumEtJet=new TH1F(name,title,2000,0.0,2000.0);
	AddHistogram(Hist.fSumEtJet,Folder);
	sprintf(name,"SqrtSumEtJet");
	sprintf(title," %s%s",algoname,": corrected #sqrt{#sum E^{jets}_{T}}");
	Hist.fSqrtSumEtJet=new TH1F(name,title,180,0.0,45.0);
	AddHistogram(Hist.fSqrtSumEtJet,Folder);
	sprintf(name,"Met4VsJetFr");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T}");
	Hist.fMet4VsJetFr=new TProfile(name,title,20,0.0,1.0,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsJetFr,Folder);

	sprintf(name,"Met4VsSqrtSumEtJet");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T} vs #sqrt{#Sigma E^{jets}_{T}}");
	Hist.fMet4VsSqrtSumEtJet=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
	AddHistogram(Hist.fMet4VsSqrtSumEtJet,Folder);
	sprintf(name,"JetFrVsSumEt");
	sprintf(title," %s%s",algoname,": #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T} vs #Sigma E^{tot}_{T}}");
	Hist.fJetFrVsSumEt=new TProfile(name,title,40,0.0,400.0,-1.0,10.0);
	AddHistogram(Hist.fJetFrVsSumEt,Folder);

	sprintf(name,"Met4X_vs_Met4Y");
	sprintf(title," %s%s",algoname,": meth4 #slash{E}_{T}^{X} vs. #slash{E}_{T}^{Y}");
	Hist.fMet4X_vs_Met4Y=new TH2F(name,title,40,-40.0,40.0,40,-40.0,40.0);
	AddHistogram(Hist.fMet4X_vs_Met4Y,Folder);

	for(int i=0; i<10; i++)
	{

		sprintf(name,"Met4[%i]",i);
		sprintf(title,"%s: meth4 #slash{E}_{T}, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
		Hist.fMet4[i]=new TH1F(name,title,500,0.0,500.0);
		AddHistogram(Hist.fMet4[i],Folder);
		sprintf(name,"Met4Phi[%i]",i);
		sprintf(title,"%s: #phi of meth4 #slash{E}_{T}, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
		Hist.fMet4Phi[i]=new TH1F(name,title,128,0.0,6.4);
		AddHistogram(Hist.fMet4Phi[i],Folder);

		//---- any number of vertices
		sprintf(name,"Met4X[%i]",i);
		sprintf(title,"%s: meth4 MEtX, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
		Hist.fMet4X[i]=new TH1F(name,title,400,-200.0,200.0);
		AddHistogram(Hist.fMet4X[i],Folder);
		sprintf(name,"Met4Y[%i]",i);
		sprintf(title,"%s: meth4 MEtY, %i<#sqrt{#sum E_{T}^{corr}}<%i",algoname,i*2,(i+1)*2);
		Hist.fMet4Y[i]=new TH1F(name,title,400,-200.0,200.0);
		AddHistogram(Hist.fMet4Y[i],Folder);

		//---- Njet_thr*=0
		sprintf(name,"Met4X_noJet[%i]",i);
		sprintf(title,"#slash{E}_{T}^{X}: N_{jet15}=0 and %i<#sqrt{#sum E_{T}^{uncl}}<%i",i*2,(i+1)*2);
		Hist.fMet4X_noJet[i]=new TH1F(name,title,400,-200.0,200.0);
		AddHistogram(Hist.fMet4X_noJet[i],Folder);
		sprintf(name,"Met4Y_noJet[%i]",i);
		sprintf(title,"#slash{E}_{T}^{Y}: N_{jet15}=0 and %i<#sqrt{#sum E_{T}^{uncl}}<%i",i*2,(i+1)*2);
		Hist.fMet4Y_noJet[i]=new TH1F(name,title,400,-200.0,200.0);
		AddHistogram(Hist.fMet4Y_noJet[i],Folder);
		//---- Njet_thr*>0
		sprintf(name,"Met4X_withJet[%i]",i);
		sprintf(title,"#slash{E}_{T}^{X}: N_{jet15}>0 and %i<#sqrt{#sum E_{T}^{uncl}}<%i",i*2,(i+1)*2);
		Hist.fMet4X_withJet[i]=new TH1F(name,title,400,-200.0,200.0);
		AddHistogram(Hist.fMet4X_withJet[i],Folder);
		sprintf(name,"Met4Y_withJet[%i]",i);
		sprintf(title,"#slash{E}_{T}^{Y}: N_{jet15}>0 and %i<#sqrt{#sum E_{T}^{uncl}}<%i",i*2,(i+1)*2);
		Hist.fMet4Y_withJet[i]=new TH1F(name,title,400,-200.0,200.0);
		AddHistogram(Hist.fMet4Y_withJet[i],Folder);
	}

	//______________________________________________________________________ histograms for MetModel-2
	sprintf(name,"GenV2Met_def_proj");
	sprintf(title,"projection of generated #slash{E}_{T} (default parametrization) on #slash{E}_{T}^{det} direction");
	Hist.fGenV2Met_def_proj=new TH1F(name,title,1000,-500.0,500.0);
	AddHistogram(Hist.fGenV2Met_def_proj,Folder);
	sprintf(name,"GenV2MetXY_def_proj");
	sprintf(title,"2D projection of generated #slash{E}_{T} (default parametrization) on #slash{E}_{T}^{det} direction");
	Hist.fGenV2MetXY_def_proj=new TH2F(name,title,100,-100.0,100.0,100,-100.0,100.0);
	AddHistogram(Hist.fGenV2MetXY_def_proj,Folder);
	sprintf(name,"GenV2Met_def_projNjet15");
	sprintf(title,"projection of generated #slash{E}_{T} (default parametrization) on #slash{E}_{T}^{det} direction if N_{jet15}>0");
	Hist.fGenV2Met_def_projNjet15=new TH1F(name,title,1000,-500.0,500.0);
	AddHistogram(Hist.fGenV2Met_def_projNjet15,Folder);
	sprintf(name,"GenV2MetXY_def_projNjet15");
	sprintf(title,"2D projection of generated #slash{E}_{T} (default parametrization) on #slash{E}_{T}^{det} direction if N_{jet15}>0");
	Hist.fGenV2MetXY_def_projNjet15=new TH2F(name,title,100,-100.0,100.0,100,-100.0,100.0);
	AddHistogram(Hist.fGenV2MetXY_def_projNjet15,Folder);
	sprintf(name,"GenV2_dPhiMet_def");
	sprintf(title,"#Delta#phi(#slash{E}_{T}^{gen}-#slash{E}_{T}^{det}) (default parametrization)");
	Hist.fGenV2_dPhiMet_def=new TH1F(name,title,157,0.0,TMath::Pi());
	AddHistogram(Hist.fGenV2_dPhiMet_def,Folder);
	sprintf(name,"GenV2_dPhiMet_def_Njet15");
	sprintf(title,"#Delta#phi(#slash{E}_{T}^{gen}-#slash{E}_{T}^{det}) (default parametrization) if N_{jet15}>0");
	Hist.fGenV2_dPhiMet_def_Njet15=new TH1F(name,title,157,0.0,TMath::Pi());
	AddHistogram(Hist.fGenV2_dPhiMet_def_Njet15,Folder);


	sprintf(name,"GenV2Met_max_def_proj");
	sprintf(title,"projection of maximum #slash{E}_{T}^{gen} (default parametrization) on #slash{E}_{T}^{det} direction");
	Hist.fGenV2Met_max_def_proj=new TH1F(name,title,1000,-500.0,500.0);
	AddHistogram(Hist.fGenV2Met_max_def_proj,Folder);
	sprintf(name,"GenV2MetXY_max_def_proj");
	sprintf(title,"2D projection of maximum #slash{E}_{T}^{gen} (default parametrization) on #slash{E}_{T}^{det} direction");
	Hist.fGenV2MetXY_max_def_proj=new TH2F(name,title,100,-100.0,100.0,100,-100.0,100.0);
	AddHistogram(Hist.fGenV2MetXY_max_def_proj,Folder);
	sprintf(name,"GenV2Met_max_def_projNjet15");
	sprintf(title,"projection of maximum #slash{E}_{T}^{gen} (default parametrization) on #slash{E}_{T}^{det} direction if N_{jet15}>0");
	Hist.fGenV2Met_max_def_projNjet15=new TH1F(name,title,1000,-500.0,500.0);
	AddHistogram(Hist.fGenV2Met_max_def_projNjet15,Folder);
	sprintf(name,"GenV2MetXY_max_def_projNjet15");
	sprintf(title,"2D projection of maximum #slash{E}_{T}^{gen} (default parametrization) on #slash{E}_{T}^{det} direction if N_{jet15}>0");
	Hist.fGenV2MetXY_max_def_projNjet15=new TH2F(name,title,100,-100.0,100.0,100,-100.0,100.0);
	AddHistogram(Hist.fGenV2MetXY_max_def_projNjet15,Folder);
	sprintf(name,"GenV2_dPhiMet_max_def");
	sprintf(title,"#Delta#phi(#slash{E}_{T}^{gen}-#slash{E}_{T}^{det}) (default parametrization), using maximum #slash{E}_{T}^{gen}");
	Hist.fGenV2_dPhiMet_max_def=new TH1F(name,title,157,0.0,TMath::Pi());
	AddHistogram(Hist.fGenV2_dPhiMet_max_def,Folder);
	sprintf(name,"GenV2_dPhiMet_max_def_Njet15");
	sprintf(title,"#Delta#phi(#slash{E}_{T}^{gen}-#slash{E}_{T}^{det}) (default parametrization) if N_{jet15}>0, using maximum #slash{E}_{T}^{gen}");
	Hist.fGenV2_dPhiMet_max_def_Njet15=new TH1F(name,title,157,0.0,TMath::Pi());
	AddHistogram(Hist.fGenV2_dPhiMet_max_def_Njet15,Folder);


	sprintf(name,"GenV2Met_max_def_projMet10");
	sprintf(title,"projection of maximum #slash{E}_{T}^{gen} (default parametrization) on #slash{E}_{T}^{det}>10 direction");
	Hist.fGenV2Met_max_def_projMet10=new TH1F(name,title,1000,-500.0,500.0);
	AddHistogram(Hist.fGenV2Met_max_def_projMet10,Folder);
	sprintf(name,"GenV2MetXY_max_def_projMet10");
	sprintf(title,"2D projection of maximum #slash{E}_{T}^{gen} (default parametrization) on #slash{E}_{T}^{det}>10 direction");
	Hist.fGenV2MetXY_max_def_projMet10=new TH2F(name,title,100,-100.0,100.0,100,-100.0,100.0);
	AddHistogram(Hist.fGenV2MetXY_max_def_projMet10,Folder);
	sprintf(name,"GenV2Met_max_def_projMet10Njet15");
	sprintf(title,"projection of maximum #slash{E}_{T}^{gen} (default parametrization) on #slash{E}_{T}^{det}>10 direction if N_{jet15}>0");
	Hist.fGenV2Met_max_def_projMet10Njet15=new TH1F(name,title,1000,-500.0,500.0);
	AddHistogram(Hist.fGenV2Met_max_def_projMet10Njet15,Folder);
	sprintf(name,"GenV2MetXY_max_def_projMet10Njet15");
	sprintf(title,"2D projection of maximum #slash{E}_{T}^{gen} (default parametrization) on #slash{E}_{T}^{det}>10 direction if N_{jet15}>0");
	Hist.fGenV2MetXY_max_def_projMet10Njet15=new TH2F(name,title,100,-100.0,100.0,100,-100.0,100.0);
	AddHistogram(Hist.fGenV2MetXY_max_def_projMet10Njet15,Folder);
	sprintf(name,"GenV2_dPhiMet_max_def_Met10");
	sprintf(title,"#Delta#phi(#slash{E}_{T}^{gen}-#slash{E}_{T}^{det}) (default parametrization), if #slash{E}_{T}^{det}>10 and using max #slash{E}_{T}^{gen}");
	Hist.fGenV2_dPhiMet_max_def_Met10=new TH1F(name,title,157,0.0,TMath::Pi());
	AddHistogram(Hist.fGenV2_dPhiMet_max_def_Met10,Folder);
	sprintf(name,"GenV2_dPhiMet_max_def_Met10Njet15");
	sprintf(title,"#Delta#phi(#slash{E}_{T}^{gen}-#slash{E}_{T}^{det}) (default parametrization) if #slash{E}_{T}^{det}>10 N_{jet15}>0, using max #slash{E}_{T}^{gen}");
	Hist.fGenV2_dPhiMet_max_def_Met10Njet15=new TH1F(name,title,157,0.0,TMath::Pi());
	AddHistogram(Hist.fGenV2_dPhiMet_max_def_Met10Njet15,Folder);


	sprintf(name,"GenV2Met_def");
	sprintf(title," %s%s",algoname,": generated Model-2 #slash{E}_{T} (default parametrization)");
	Hist.fGenV2Met_def=new TH1F(name,title,500,0.0,500.0);
	AddHistogram(Hist.fGenV2Met_def,Folder);
	sprintf(name,"GenV2MetPhi_def");
	sprintf(title," %s%s",algoname,": #phi of generated Model-2 #slash{E}_{T} (default parametrization)");
	Hist.fGenV2MetPhi_def=new TH1F(name,title,128,0.0,6.4);
	AddHistogram(Hist.fGenV2MetPhi_def,Folder);
	sprintf(name,"GenV2MetX_def");
	sprintf(title," %s%s",algoname,": generated Model-2 #slash{E}_{T}^{X} (default parametrization)");
	Hist.fGenV2MetX_def=new TH1F(name,title,400,-200.0,200.0);
	AddHistogram(Hist.fGenV2MetX_def,Folder);
	sprintf(name,"GenV2MetY_def");
	sprintf(title," %s%s",algoname,": generated Model-2 #slash{E}_{T}^{Y} (default parametrization)");
	Hist.fGenV2MetY_def=new TH1F(name,title,400,-200.0,200.0);
	AddHistogram(Hist.fGenV2MetY_def,Folder);

	sprintf(name,"MetGenV2VsNjet");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs N_{jet}");
	Hist.fMetGenV2VsNjet=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsNjet,Folder);
	sprintf(name,"MetGenV2VsNjet_th5");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs N_{jet} with E_{T}>5 GeV");
	Hist.fMetGenV2VsNjet_th5=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsNjet_th5,Folder);
	sprintf(name,"MetGenV2VsNjet_th10");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs N_{jet} with E_{T}>10 GeV");
	Hist.fMetGenV2VsNjet_th10=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsNjet_th10,Folder);
	sprintf(name,"MetGenV2VsNjet_th15");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs N_{jet} with E_{T}>15 GeV");
	Hist.fMetGenV2VsNjet_th15=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsNjet_th15,Folder);
	sprintf(name,"MetGenV2VsNjet_th20");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs N_{jet} with E_{T}>20 GeV");
	Hist.fMetGenV2VsNjet_th20=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsNjet_th20,Folder);
	sprintf(name,"MetGenV2VsNjet_th25");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs N_{jet} with E_{T}>25 GeV");
	Hist.fMetGenV2VsNjet_th25=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsNjet_th25,Folder);
	sprintf(name,"MetGenV2VsNjet_th30");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs N_{jet} with E_{T}>30 GeV");
	Hist.fMetGenV2VsNjet_th30=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsNjet_th30,Folder);
	sprintf(name,"MetGenV2VsNjet_th35");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs N_{jet} with E_{T}>35 GeV");
	Hist.fMetGenV2VsNjet_th35=new TProfile(name,title,40,-0.5,39.5,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsNjet_th35,Folder);


	sprintf(name,"MetGenV2VsRun");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} Vs. Run for N_{vx12}=1");
	Hist.fMetGenV2VsRun=new TProfile(name,title,100,130000.0,230000.0,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsRun,Folder);
	sprintf(name,"MetGenV2VsNvx12");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs N_{vx12}");
	Hist.fMetGenV2VsNvx12=new TProfile(name,title,10,-0.5,9.5,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsNvx12,Folder);

	sprintf(name,"MetGenV2VsJetFr");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs #Sigma E^{jets}_{T}/#Sigma E^{tot}_{T}");
	Hist.fMetGenV2VsJetFr=new TProfile(name,title,20,0.0,1.0,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsJetFr,Folder);
	sprintf(name,"MetGenV2VsSqrtSumEtJet");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T} vs #sqrt{#Sigma E^{jets}_{T}}");
	Hist.fMetGenV2VsSqrtSumEtJet=new TProfile(name,title,45,0.0,45.0,-1.0,1000.0);
	AddHistogram(Hist.fMetGenV2VsSqrtSumEtJet,Folder);
	sprintf(name,"GenV2MetX_vs_MetY");
	sprintf(title," %s%s",algoname,": generated Model-2 (def parametrization) #slash{E}_{T}^{X} vs. #slash{E}_{T}^{Y}");
	Hist.fGenV2MetX_vs_MetY=new TH2F(name,title,40,-40.0,40.0,40,-40.0,40.0);
	AddHistogram(Hist.fGenV2MetX_vs_MetY,Folder);

	sprintf(name,"Met4GenV2Met");
	sprintf(title," %s%s",algoname,": #slash{E}_{T}-#slash{E}_{T}^{gen}");
	Hist.fMet4GenV2Met=new TH1F(name,title,400,-200.0,200.0);
	AddHistogram(Hist.fMet4GenV2Met,Folder);
	sprintf(name,"Met4_vs_Met4GenV2Met");
	sprintf(title," %s%s",algoname,": #slash{E}_{T} .vs. #slash{E}_{T}-#slash{E}_{T}^{gen}");
	Hist.fMet4_vs_Met4GenV2Met=new TH2F(name,title,60,-150.0,150.0,30,0.0,150.0);
	AddHistogram(Hist.fMet4_vs_Met4GenV2Met,Folder);
	sprintf(name,"Met4_vs_GenV2Met");
	sprintf(title," %s%s",algoname,": #slash{E}_{T} .vs. #slash{E}_{T}^{gen}");
	Hist.fMet4_vs_GenV2Met=new TH2F(name,title,50,0.0,250.0,50,0.0,250.0);
	AddHistogram(Hist.fMet4_vs_GenV2Met,Folder);

	sprintf(name,"MetCorr_vs_dZvx");
	sprintf(title," %s%s",algoname,": corrected #slash{E}_{T} .vs. #DeltaZ_{vx} for N_{vx12}>1");
	Hist.fMetCorr_vs_dZvx=new TH2F(name,title,40,0.0,200.0,100,0.0,250.0);
	AddHistogram(Hist.fMetCorr_vs_dZvx,Folder);
	sprintf(name,"MetCorr_vs_dZvxWorse");
	sprintf(title," %s%s",algoname,": corrected #slash{E}_{T} .vs. worse #DeltaZ_{vx} for N_{vx}>1");
	Hist.fMetCorr_vs_dZvxWorse=new TH2F(name,title,40,0.0,200.0,100,0.0,250.0);
	AddHistogram(Hist.fMetCorr_vs_dZvxWorse,Folder);

	sprintf(name,"MetCorr_vs_ZvxBest");
	sprintf(title," %s%s",algoname,": corrected #slash{E}_{T} .vs. Z_{vx} of best vertex");
	Hist.fMetCorr_vs_ZvxBest=new TH2F(name,title,30,0.0,150.0,100,0.0,250.0);
	AddHistogram(Hist.fMetCorr_vs_ZvxBest,Folder);
	sprintf(name,"MetCorr_vs_Zvx2ndBest");
	sprintf(title," %s%s",algoname,": corrected #slash{E}_{T} .vs. Z_{vx} of 2nd best vertex");
	Hist.fMetCorr_vs_Zvx2ndBest=new TH2F(name,title,30,0.0,150.0,100,0.0,250.0);
	AddHistogram(Hist.fMetCorr_vs_Zvx2ndBest,Folder);
	sprintf(name,"MetCorr_vs_ZvxWorse");
	sprintf(title," %s%s",algoname,": corrected #slash{E}_{T} .vs. largest Z_{vx}");
	Hist.fMetCorr_vs_ZvxWorse=new TH2F(name,title,30,0.0,150.0,100,0.0,250.0);
	AddHistogram(Hist.fMetCorr_vs_ZvxWorse,Folder);
	sprintf(name,"Met0Met4_vs_Zvx");
	sprintf(title," %s%s",algoname,": #slash{E}_{T}^{meth0}-#slash{E}_{T}^{meth4} .vs. Z_{vx}");
	Hist.fMet0Met4_vs_Zvx=new TH2F(name,title,20,0.0,60.0,200,-250.0,250.0);
	AddHistogram(Hist.fMet0Met4_vs_Zvx,Folder);
	sprintf(name,"Met0_vs_Zvx");
	sprintf(title," %s%s",algoname,": #slash{E}_{T}^{meth0} .vs. Z_{vx}");
	Hist.fMet0_vs_Zvx=new TH2F(name,title,20,0.0,60.0,100,0.0,250.0);
	AddHistogram(Hist.fMet0_vs_Zvx,Folder);
	sprintf(name,"MetCorrVsZvxBest");
	sprintf(title," %s%s",algoname,": corrected #slash{E}_{T} .vs. Z_{vx} of best vertex");
	Hist.fMetCorrVsZvxBest=new TProfile(name,title,30,0.0,150.0,-1.0,250.0);
	AddHistogram(Hist.fMetCorrVsZvxBest,Folder);
	sprintf(name,"MetCorrVsZvx2ndBest");
	sprintf(title," %s%s",algoname,": corrected #slash{E}_{T} .vs. Z_{vx} of 2nd best vertex");
	Hist.fMetCorrVsZvx2ndBest=new TProfile(name,title,30,0.0,150.0,-1.0,250.0);
	AddHistogram(Hist.fMetCorrVsZvx2ndBest,Folder);
	sprintf(name,"MetCorrVsZvxWorse");
	sprintf(title," %s%s",algoname,": corrected #slash{E}_{T} .vs. largest Z_{vx}");
	Hist.fMetCorrVsZvxWorse=new TProfile(name,title,30,0.0,150.0,-1.0,250.0);
	AddHistogram(Hist.fMetCorrVsZvxWorse,Folder);

	return;
}

//________________________________________________booking pho-jet match histo
void JetFilterModuleV2::BookMatchStudyHistograms(MatchStudyHisto_t& Hist, const char* Folder, const char* algoname) {

	char name [200];
	char title[200];

	sprintf(name,"MatchNtwr");
	sprintf(title,"%s: Number of towers in matched jet",algoname);
	Hist.fMatchNtwr=new TH1F(name,title,50,-0.5,49.5);
	AddHistogram(Hist.fMatchNtwr,Folder);
	sprintf(name,"MatchDelR");
	sprintf(title,"%s: #Delta R between jet and matched EM object",algoname);
	Hist.fMatchDelR=new TH1F(name,title,500,0.0,10.0);
	AddHistogram(Hist.fMatchDelR,Folder);
	sprintf(name,"MatchDelPhi");
	sprintf(title,"%s: #Delta#phi between jet and matched EM object",algoname);
	Hist.fMatchDelPhi=new TH1F(name,title,200,0.0,4.0);
	AddHistogram(Hist.fMatchDelPhi,Folder);
	sprintf(name,"MatchDelEta");
	sprintf(title,"%s: #Delta#eta between jet and matched EM object",algoname);
	Hist.fMatchDelEta=new TH1F(name,title,200,0.0,4.0);
	AddHistogram(Hist.fMatchDelEta,Folder);
	sprintf(name,"MatchDelEtaDet");
	sprintf(title,"%s: #Delta#eta_{det} between jet and matched EM object",algoname);
	Hist.fMatchDelEtaDet=new TH1F(name,title,200,0.0,4.0);
	AddHistogram(Hist.fMatchDelEtaDet,Folder);
	sprintf(name,"MatchNmatch");
	sprintf(title," %s%s",algoname,": Number of matched EM objects");
	Hist.fMatchNmatch=new TH1F(name,title,10,-0.5,9.5);
	AddHistogram(Hist.fMatchNmatch,Folder);
	sprintf(name,"MatchEt_raw_b");
	sprintf(title," %s%s",algoname,": raw E_{T} of jet before matching");
	Hist.fMatchEt_raw_b=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fMatchEt_raw_b,Folder);
	sprintf(name,"MatchEt_raw_a");
	sprintf(title," %s%s",algoname,": raw E_{T} of jet after matching");
	Hist.fMatchEt_raw_a=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fMatchEt_raw_a,Folder);
	sprintf(name,"MatchEt_lev6_b");
	sprintf(title," %s%s",algoname,": lev6 E_{T} of jet before matching");
	Hist.fMatchEt_lev6_b=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fMatchEt_lev6_b,Folder);
	sprintf(name,"MatchEt_lev6_a");
	sprintf(title," %s%s",algoname,": lev6 E_{T} of jet after matching");
	Hist.fMatchEt_lev6_a=new TH1F(name,title,1000,0.0,1000.0);
	AddHistogram(Hist.fMatchEt_lev6_a,Folder);
	sprintf(name,"MatchEtJet2EtPho_b");
	sprintf(title," %s%s",algoname,": raw E^{jet}_{T}/E^{EMobj}_{T} before removing EM object");
	Hist.fMatchEtJet2EtPho_b=new TH1F(name,title,1000,0.0,10.0);
	AddHistogram(Hist.fMatchEtJet2EtPho_b,Folder);
	sprintf(name,"MatchEtJet2EtPho_a");
	sprintf(title," %s%s",algoname,": raw E^{jet}_{T}/E^{EMobj}_{T} after removing EM object");
	Hist.fMatchEtJet2EtPho_a=new TH1F(name,title,1000,0.0,10.0);
	AddHistogram(Hist.fMatchEtJet2EtPho_a,Folder);
	sprintf(name,"MatchDelRoldnew");
	sprintf(title,"%s: #Delta R between new and old jet",algoname);
	Hist.fMatchDelRoldnew=new TH1F(name,title,500,0.0,10.0);
	AddHistogram(Hist.fMatchDelRoldnew,Folder);
	sprintf(name,"MatchDelEtaDetoldnew");
	sprintf(title,"%s: #Delta#eta_{det} between new and old jet",algoname);
	Hist.fMatchDelEtaDetoldnew=new TH1F(name,title,200,-4.0,4.0);
	AddHistogram(Hist.fMatchDelEtaDetoldnew,Folder);

	return;
}


//_________________________________________________ filling met cleanup histo (for studies)
void JetFilterModuleV2::BookMetCleanupHistograms(MetCleanupHisto_t& Hist, const char* Folder, const char* algoname) {

	char name [200];
	char title[200];

	sprintf(name,"CleanupDelPhiVsEtaDet_jet");
	sprintf(title," %s%s",algoname,": d#phi=#phi_{jet}-#phi_{#slash{E}_{T}} vs. #eta^{det}_{jet}");
	Hist.fCleanupDelPhiVsEtaDet_jet= new TH2F(name,title,37,-3.7,3.7,16,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhiVsEtaDet_jet,Folder);
	sprintf(name,"CleanupDelPhiVsEtaDet_em");
	sprintf(title," %s%s",algoname,": d#phi=#phi_{e,#gamma}-#phi_{#slash{E}_{T}} vs. #eta^{det}_{e,#gamma}");
	Hist.fCleanupDelPhiVsEtaDet_em= new TH2F(name,title,37,-3.7,3.7,16,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhiVsEtaDet_em,Folder);

	sprintf(name,"CleanupMetEtemVsDelPhi");
	sprintf(title," %s%s",algoname,": #slash{E}_{T}/E_{T}^{e,#gamma} vs. d#phi=#phi_{e,#gamma}-#phi_{#slash{E}_{T}}");
	Hist.fCleanupMetEtemVsDelPhi= new TH2F(name,title,16,0.0,TMath::Pi(),50,0.0,10.0);
	AddHistogram(Hist.fCleanupMetEtemVsDelPhi,Folder);
	sprintf(name,"CleanupDelPhi_metem");
	sprintf(title," %s%s",algoname,": d#phi=#phi_{e,#gamma}-#phi_{#slash{E}_{T}}");
	Hist.fCleanupDelPhi_metem= new TH1F(name,title,128,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_metem,Folder);
	sprintf(name,"CleanupMetEtem");
	sprintf(title," %s%s",algoname,": #slash{E}_{T}/E_{T}^{e,#gamma}");
	Hist.fCleanupMetEtem= new TH1F(name,title,800,0.0,20.0);
	AddHistogram(Hist.fCleanupMetEtem,Folder);

	for(int i=0; i<3; i++)
	{
		char range[200];
		if(i==0) sprintf(range,", jet in crack");
		if(i==1) sprintf(range,", jet away from crack");
		if(i==2) sprintf(range,", jet in |#eta_{det}|>2.6");

		sprintf(name,"CleanupMetEtjetVsDelPhi[%i]",i);
		sprintf(title," %s: #slash{E}_{T}/E_{T}^{jet} vs. d#phi=#phi_{jet}-#phi_{#slash{E}_{T}}%s",algoname,range);
		Hist.fCleanupMetEtjetVsDelPhi[i]= new TH2F(name,title,16,0.0,TMath::Pi(),50,0.0,10.0);
		AddHistogram(Hist.fCleanupMetEtjetVsDelPhi[i],Folder);
		sprintf(name,"CleanupDelPhi_metjet[%i]",i);
		sprintf(title," %s: d#phi=#phi_{jet}-#phi_{#slash{E}_{T}}%s",algoname,range);
		Hist.fCleanupDelPhi_metjet[i]= new TH1F(name,title,128,0.0,TMath::Pi());
		AddHistogram(Hist.fCleanupDelPhi_metjet[i],Folder);
		sprintf(name,"CleanupMetEtjet[%i]",i);
		sprintf(title," %s: #slash{E}_{T}/E_{T}^{jet}%s",algoname,range);
		Hist.fCleanupMetEtjet[i]= new TH1F(name,title,800,0.0,20.0);
		AddHistogram(Hist.fCleanupMetEtjet[i],Folder);
	}

	sprintf(name,"CleanupDelPhi_RawMetRawJet3");
	sprintf(title," %s%s",algoname,": #Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}^{raw}}, closest jet with E_{T}^{raw}>3 GeV");
	Hist.fCleanupDelPhi_RawMetRawJet3= new TH1F(name,title,126,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_RawMetRawJet3,Folder);
	sprintf(name,"CleanupDelPhi_RawMetRawJet5");
	sprintf(title," %s%s",algoname,": #Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}^{raw}}, closest jet with E_{T}^{raw}>5 GeV");
	Hist.fCleanupDelPhi_RawMetRawJet5= new TH1F(name,title,126,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_RawMetRawJet5,Folder);
	sprintf(name,"CleanupDelPhi_RawMetRawJet10");
	sprintf(title," %s%s",algoname,": #Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}^{raw}}, closest jet with E_{T}^{raw}>10 GeV");
	Hist.fCleanupDelPhi_RawMetRawJet10= new TH1F(name,title,126,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_RawMetRawJet10,Folder);
	sprintf(name,"CleanupDelPhi_RawMet1stRawJet");
	sprintf(title," %s%s",algoname,": #Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}^{raw}}, 1^{st} jet with E_{T}^{raw}>3 GeV");
	Hist.fCleanupDelPhi_RawMet1stRawJet= new TH1F(name,title,126,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_RawMet1stRawJet,Folder);

	//----------------------------------------- histograms added on 04/03/07
	sprintf(name,"CleanupDelPhi_metjet5");
	sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>5 GeV");
	Hist.fCleanupDelPhi_metjet5= new TH1F(name,title,128,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_metjet5,Folder);
	sprintf(name,"CleanupDelPhi_metjet15");
	sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>15 GeV");
	Hist.fCleanupDelPhi_metjet15= new TH1F(name,title,128,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_metjet15,Folder);
	sprintf(name,"CleanupDelPhi_metjet25");
	sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>25 GeV");
	Hist.fCleanupDelPhi_metjet25= new TH1F(name,title,128,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_metjet25,Folder);
	sprintf(name,"CleanupDelPhi_met10jet15");
	sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>15 GeV and #slash{E}_{T}>10 GeV");
	Hist.fCleanupDelPhi_met10jet15= new TH1F(name,title,128,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_met10jet15,Folder);
	sprintf(name,"CleanupDelPhi_metjet1st15");
	sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, 1^{st} jet with E_{T}^{lev6}>15 GeV");
	Hist.fCleanupDelPhi_metjet1st15= new TH1F(name,title,128,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_metjet1st15,Folder);
	sprintf(name,"CleanupDelPhi_met10jet1st15");
	sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, 1^{st} jet with E_{T}^{lev6}>15 GeV and #slash{E}_{T}>10 GeV");
	Hist.fCleanupDelPhi_met10jet1st15= new TH1F(name,title,128,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_met10jet1st15,Folder);
	sprintf(name,"CleanupDelPhi_metjet15_dPhiMin");
	sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>15 GeV and min(#Delta#phi)");
	Hist.fCleanupDelPhi_metjet15_dPhiMin= new TH1F(name,title,128,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_metjet15_dPhiMin,Folder);
	sprintf(name,"CleanupDelPhi_met10jet15_dPhiMin");
	sprintf(title,"#Delta#phi=#phi_{jet}-#phi_{#slash{E}_{T}}, E_{T}^{lev6}>15 GeV, min(#Delta#phi), and #slash{E}_{T}>10 GeV");
	Hist.fCleanupDelPhi_met10jet15_dPhiMin= new TH1F(name,title,128,0.0,TMath::Pi());
	AddHistogram(Hist.fCleanupDelPhi_met10jet15_dPhiMin,Folder);

	return;
}


//______________________________________________________________ booking Pho-Met study histograms
void JetFilterModuleV2::BookPhoMetStudyHistograms(PhoMetStudyHisto_t& Hist, const char* Folder) {
	char name [200];
	char title[200];

	sprintf(name,"eta");
	sprintf(title,"|#eta_{det}^{#gamma}|");
	Hist.fPhoMet_eta=new TH1F(name,title,120,0.0,1.2);
	AddHistogram(Hist.fPhoMet_eta,Folder);
	sprintf(name,"xces");
	sprintf(title,"|X_{CES}|");
	Hist.fPhoMet_xces=new TH1F(name,title,86,0.0,23.0);
	AddHistogram(Hist.fPhoMet_xces,Folder);
	sprintf(name,"zces");
	sprintf(title,"|Z_{CES}|");
	Hist.fPhoMet_zces=new TH1F(name,title,240,0.0,240.0);
	AddHistogram(Hist.fPhoMet_zces,Folder);

	sprintf(name,"signeta");
	sprintf(title,"sign(Z_{vx})*#eta_{det}^{#gamma}");
	Hist.fPhoMet_signeta=new TH1F(name,title,240,-1.2,1.2);
	AddHistogram(Hist.fPhoMet_signeta,Folder);
	sprintf(name,"signzces");
	sprintf(title,"sign(Z_{vx})*Z_{CES}");
	Hist.fPhoMet_signzces=new TH1F(name,title,480,-240.0,240.0);
	AddHistogram(Hist.fPhoMet_signzces,Folder);
	sprintf(name,"MetSig5dPhi02_signeta");
	sprintf(title,"sign(Z_{vx})*#eta_{det}^{#gamma} if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_MetSig5dPhi02_signeta=new TH1F(name,title,240,-1.2,1.2);
	AddHistogram(Hist.fPhoMet_MetSig5dPhi02_signeta,Folder);
	sprintf(name,"MetSig5dPhi02_signzces");
	sprintf(title,"sign(Z_{vx})*Z_{CES} if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_MetSig5dPhi02_signzces=new TH1F(name,title,480,-240.0,240.0);
	AddHistogram(Hist.fPhoMet_MetSig5dPhi02_signzces,Folder);
	sprintf(name,"Fr_MetSig5dPhi02_signeta");
	sprintf(title,"Fraction of photons as a function of sign(Z_{vx})*#eta_{det}^{#gamma} if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_Fr_MetSig5dPhi02_signeta=new TH1F(name,title,240,-1.2,1.2);
	AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_signeta,Folder);
	sprintf(name,"Fr_MetSig5dPhi02_signzces");
	sprintf(title,"Fraction of photons as a function of sign(Z_{vx})*Z_{CES} if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_Fr_MetSig5dPhi02_signzces=new TH1F(name,title,480,-240.0,240.0);
	AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_signzces,Folder);

	sprintf(name,"dPhi02_eta");
	sprintf(title,"|#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
	Hist.fPhoMet_dPhi02_eta=new TH1F(name,title,120,0.0,1.2);
	AddHistogram(Hist.fPhoMet_dPhi02_eta,Folder);
	sprintf(name,"dPhi02_xces");
	sprintf(title,"|X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
	Hist.fPhoMet_dPhi02_xces=new TH1F(name,title,86,0.0,23);
	AddHistogram(Hist.fPhoMet_dPhi02_xces,Folder);
	sprintf(name,"dPhi02_zces");
	sprintf(title,"|Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
	Hist.fPhoMet_dPhi02_zces=new TH1F(name,title,240,0.0,240);
	AddHistogram(Hist.fPhoMet_dPhi02_zces,Folder);
	sprintf(name,"dPhi01_eta");
	sprintf(title,"|#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
	Hist.fPhoMet_dPhi01_eta=new TH1F(name,title,120,0.0,1.2);
	AddHistogram(Hist.fPhoMet_dPhi01_eta,Folder);
	sprintf(name,"dPhi01_xces");
	sprintf(title,"|X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
	Hist.fPhoMet_dPhi01_xces=new TH1F(name,title,86,0.0,23);
	AddHistogram(Hist.fPhoMet_dPhi01_xces,Folder);
	sprintf(name,"dPhi01_zces");
	sprintf(title,"|Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
	Hist.fPhoMet_dPhi01_zces=new TH1F(name,title,240,0.0,240);
	AddHistogram(Hist.fPhoMet_dPhi01_zces,Folder);

	sprintf(name,"Fr_dPhi02_eta");
	sprintf(title,"Fraction of photons as a function of |#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
	Hist.fPhoMet_Fr_dPhi02_eta=new TH1F(name,title,120,0.0,1.2);
	AddHistogram(Hist.fPhoMet_Fr_dPhi02_eta,Folder);
	sprintf(name,"Fr_dPhi02_xces");
	sprintf(title,"Fraction of photons as a function of |X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
	Hist.fPhoMet_Fr_dPhi02_xces=new TH1F(name,title,86,0.0,23);
	AddHistogram(Hist.fPhoMet_Fr_dPhi02_xces,Folder);
	sprintf(name,"Fr_dPhi02_zces");
	sprintf(title,"Fraction of photons as a function of |Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2");
	Hist.fPhoMet_Fr_dPhi02_zces=new TH1F(name,title,240,0.0,240);
	AddHistogram(Hist.fPhoMet_Fr_dPhi02_zces,Folder);
	sprintf(name,"Fr_dPhi01_eta");
	sprintf(title,"Fraction of photons as a function of |#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
	Hist.fPhoMet_Fr_dPhi01_eta=new TH1F(name,title,120,0.0,1.2);
	AddHistogram(Hist.fPhoMet_Fr_dPhi01_eta,Folder);
	sprintf(name,"Fr_dPhi01_xces");
	sprintf(title,"Fraction of photons as a function of |X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
	Hist.fPhoMet_Fr_dPhi01_xces=new TH1F(name,title,86,0.0,23);
	AddHistogram(Hist.fPhoMet_Fr_dPhi01_xces,Folder);
	sprintf(name,"Fr_dPhi01_zces");
	sprintf(title,"Fraction of photons as a function of |Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1");
	Hist.fPhoMet_Fr_dPhi01_zces=new TH1F(name,title,240,0.0,240);
	AddHistogram(Hist.fPhoMet_Fr_dPhi01_zces,Folder);

	sprintf(name,"MetSig5dPhi02_eta");
	sprintf(title,"|#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_MetSig5dPhi02_eta=new TH1F(name,title,120,0.0,1.2);
	AddHistogram(Hist.fPhoMet_MetSig5dPhi02_eta,Folder);
	sprintf(name,"MetSig5dPhi02_xces");
	sprintf(title,"|X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_MetSig5dPhi02_xces=new TH1F(name,title,86,0.0,23);
	AddHistogram(Hist.fPhoMet_MetSig5dPhi02_xces,Folder);
	sprintf(name,"MetSig5dPhi02_zces");
	sprintf(title,"|Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_MetSig5dPhi02_zces=new TH1F(name,title,240,0.0,240);
	AddHistogram(Hist.fPhoMet_MetSig5dPhi02_zces,Folder);
	sprintf(name,"MetSig5dPhi01_eta");
	sprintf(title,"|#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_MetSig5dPhi01_eta=new TH1F(name,title,120,0.0,1.2);
	AddHistogram(Hist.fPhoMet_MetSig5dPhi01_eta,Folder);
	sprintf(name,"MetSig5dPhi01_xces");
	sprintf(title,"|X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_MetSig5dPhi01_xces=new TH1F(name,title,86,0.0,23);
	AddHistogram(Hist.fPhoMet_MetSig5dPhi01_xces,Folder);
	sprintf(name,"MetSig5dPhi01_zces");
	sprintf(title,"|Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_MetSig5dPhi01_zces=new TH1F(name,title,240,0.0,240);
	AddHistogram(Hist.fPhoMet_MetSig5dPhi01_zces,Folder);

	sprintf(name,"Fr_MetSig5dPhi02_eta");
	sprintf(title,"Fraction of photons as a function of |#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_Fr_MetSig5dPhi02_eta=new TH1F(name,title,120,0.0,1.2);
	AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_eta,Folder);
	sprintf(name,"Fr_MetSig5dPhi02_xces");
	sprintf(title,"Fraction of photons as a function of |X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_Fr_MetSig5dPhi02_xces=new TH1F(name,title,86,0.0,23);
	AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_xces,Folder);
	sprintf(name,"Fr_MetSig5dPhi02_zces");
	sprintf(title,"Fraction of photons as a function of |Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_Fr_MetSig5dPhi02_zces=new TH1F(name,title,240,0.0,240);
	AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi02_zces,Folder);
	sprintf(name,"Fr_MetSig5dPhi01_eta");
	sprintf(title,"Fraction of photons as a function of |#eta_{det}^{#gamma}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_Fr_MetSig5dPhi01_eta=new TH1F(name,title,120,0.0,1.2);
	AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi01_eta,Folder);
	sprintf(name,"Fr_MetSig5dPhi01_xces");
	sprintf(title,"Fraction of photons as a function of |X_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_Fr_MetSig5dPhi01_xces=new TH1F(name,title,86,0.0,23);
	AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi01_xces,Folder);
	sprintf(name,"Fr_MetSig5dPhi01_zces");
	sprintf(title,"Fraction of photons as a function of |Z_{CES}| if #Delta#phi_{#gamma#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{#gamma}}>5");
	Hist.fPhoMet_Fr_MetSig5dPhi01_zces=new TH1F(name,title,240,0.0,240);
	AddHistogram(Hist.fPhoMet_Fr_MetSig5dPhi01_zces,Folder);


	//--------------------------- correlation between electrons and met
	sprintf(name,"Ele_eta");
	sprintf(title,"|#eta_{det}^{ele}|");
	Hist.fEleMet_eta=new TH1F(name,title,300,0.0,3.0);
	AddHistogram(Hist.fEleMet_eta,Folder);
	sprintf(name,"Ele_phi");
	sprintf(title,"#phi_{det}");
	Hist.fEleMet_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
	AddHistogram(Hist.fEleMet_phi,Folder);

	sprintf(name,"Ele_dPhi02_eta");
	sprintf(title,"|#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.2");
	Hist.fEleMet_dPhi02_eta=new TH1F(name,title,300,0.0,3.0);
	AddHistogram(Hist.fEleMet_dPhi02_eta,Folder);
	sprintf(name,"Ele_dPhi02_phi");
	sprintf(title,"#phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.2");
	Hist.fEleMet_dPhi02_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
	AddHistogram(Hist.fEleMet_dPhi02_phi,Folder);
	sprintf(name,"Ele_dPhi01_eta");
	sprintf(title,"|#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.1");
	Hist.fEleMet_dPhi01_eta=new TH1F(name,title,300,0.0,3.0);
	AddHistogram(Hist.fEleMet_dPhi01_eta,Folder);
	sprintf(name,"Ele_dPhi01_phi");
	sprintf(title,"#phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.1");
	Hist.fEleMet_dPhi01_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
	AddHistogram(Hist.fEleMet_dPhi01_phi,Folder);

	sprintf(name,"Ele_Fr_dPhi02_eta");
	sprintf(title,"Fraction of electrons as a function of |#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.2");
	Hist.fEleMet_Fr_dPhi02_eta=new TH1F(name,title,300,0.0,3.0);
	AddHistogram(Hist.fEleMet_Fr_dPhi02_eta,Folder);
	sprintf(name,"Ele_Fr_dPhi02_phi");
	sprintf(title,"Fraction of electrons as a function of #phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.2");
	Hist.fEleMet_Fr_dPhi02_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
	AddHistogram(Hist.fEleMet_Fr_dPhi02_phi,Folder);
	sprintf(name,"Ele_Fr_dPhi01_eta");
	sprintf(title,"Fraction of electrons as a function of |#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.1");
	Hist.fEleMet_Fr_dPhi01_eta=new TH1F(name,title,300,0.0,3.0);
	AddHistogram(Hist.fEleMet_Fr_dPhi01_eta,Folder);
	sprintf(name,"Ele_Fr_dPhi01_phi");
	sprintf(title,"Fraction of electrons as a function of #phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.1");
	Hist.fEleMet_Fr_dPhi01_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
	AddHistogram(Hist.fEleMet_Fr_dPhi01_phi,Folder);

	sprintf(name,"Ele_MetSig5dPhi02_eta");
	sprintf(title,"|#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
	Hist.fEleMet_MetSig5dPhi02_eta=new TH1F(name,title,300,0.0,3.0);
	AddHistogram(Hist.fEleMet_MetSig5dPhi02_eta,Folder);
	sprintf(name,"Ele_MetSig5dPhi02_phi");
	sprintf(title,"#phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
	Hist.fEleMet_MetSig5dPhi02_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
	AddHistogram(Hist.fEleMet_MetSig5dPhi02_phi,Folder);
	sprintf(name,"Ele_MetSig5dPhi01_eta");
	sprintf(title,"|#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
	Hist.fEleMet_MetSig5dPhi01_eta=new TH1F(name,title,300,0.0,3.0);
	AddHistogram(Hist.fEleMet_MetSig5dPhi01_eta,Folder);
	sprintf(name,"Ele_MetSig5dPhi01_phi");
	sprintf(title,"#phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
	Hist.fEleMet_MetSig5dPhi01_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
	AddHistogram(Hist.fEleMet_MetSig5dPhi01_phi,Folder);

	sprintf(name,"Ele_Fr_MetSig5dPhi02_eta");
	sprintf(title,"Fraction of electrons as a function of |#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
	Hist.fEleMet_Fr_MetSig5dPhi02_eta=new TH1F(name,title,300,0.0,3.0);
	AddHistogram(Hist.fEleMet_Fr_MetSig5dPhi02_eta,Folder);
	sprintf(name,"Ele_Fr_MetSig5dPhi02_phi");
	sprintf(title,"Fraction of electrons as a function of #phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.2 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
	Hist.fEleMet_Fr_MetSig5dPhi02_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
	AddHistogram(Hist.fEleMet_Fr_MetSig5dPhi02_phi,Folder);
	sprintf(name,"Ele_Fr_MetSig5dPhi01_eta");
	sprintf(title,"Fraction of electrons as a function of |#eta_{det}^{ele}| if #Delta#phi_{e-#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
	Hist.fEleMet_Fr_MetSig5dPhi01_eta=new TH1F(name,title,300,0.0,3.0);
	AddHistogram(Hist.fEleMet_Fr_MetSig5dPhi01_eta,Folder);
	sprintf(name,"Ele_Fr_MetSig5dPhi01_phi");
	sprintf(title,"Fraction of electrons as a function of #phi_{det} if #Delta#phi_{e-#slash{E}_{T}}<0.1 and #slash{E}_{T}/#sigma_{E_{T}^{ele}}>5");
	Hist.fEleMet_Fr_MetSig5dPhi01_phi=new TH1F(name,title,630,0.0,2.0*TMath::Pi());
	AddHistogram(Hist.fEleMet_Fr_MetSig5dPhi01_phi,Folder);

	return;
}

//___________________ filling analysis histograms
void JetFilterModuleV2::FillAnalysisHistograms(AnalysisHisto_t& Hist,
							const std::vector<TLorentzVector> vFinalJets,
							JetStuff jetstuff,
							const CommonStuff miscstuff,
							const TVector2 Met, const float metsig)
{
	
	std::vector<TLorentzVector> vJet	= vFinalJets;
	const int iNjet15 = GetNjets(vJet,15.0); 
	
	assert (vJet.size()>0 && "ERROR! Analysis hists are given an empty jet vector!");
	if (vJet.size()>1)
	{
		//this additional condition is reuqired in the case of running on dipho+X sample.
		//in cases when there are only 2 phos in the event, after removing them from jet list
		// the jet vec will have two empty TLorentz vectors. If you do not require any jets
		// both will return 1e-20 for Pt and fail this assertion.
		if (vJet.at(0).Pt()>0.0001 && vJet.at(1).Pt()>0.0001)
		{
			assert(vJet.at(0).Pt() >= vJet.at(1).Pt() && "ERROR! Analysis hists are given unsorted jets");
		}
	}
	// now check if the photons are sorted in Pt
	if (miscstuff.myCorrPhoton.size()>1)
	{
		assert(miscstuff.myCorrPhoton.at(0).Pt()>=miscstuff.myCorrPhoton.at(1).Pt()
		 && "ERROR! Analsis  hists are given unsorted photons");
	}


	
	Hist.fAna_MetAll->Fill(Met.Mod()); 		// MET for all events
	Hist.fAna_MetSig->Fill(metsig); 			// MET-significance for all events
	Hist.fAna_MetAllVsMetSig->Fill(metsig,Met.Mod());
	
	// the following histrograms are only for events that pass MetSig cut
	if(metsig>=fMetSig_cut)
	{
		if (Met.Mod() < MetCut()) return;			//applies a hard MEt cut

		Hist.fAna_Met->Fill(Met.Mod()); // MET
		Hist.fAna_Njet15->Fill(iNjet15); // Njet(Et>15)
		Hist.fAna_Met_MetSig->Fill(metsig, Met.Mod()); // MET
		Hist.fAna_Njet15_MetSig->Fill(metsig, iNjet15); // Njet(Et>15)

		std::vector<TLorentzVector> vJets = vFinalJets;
		// since I am using MyHtAll for both DATA/BG routines,
		// I am subtracting the MEt added to Ht by default and adding the generated MEt
		float ht=MyHtAll(vJets, jetstuff, miscstuff,3)-jetstuff.myMETcorr_th15.Mod()+Met.Mod();
		Hist.fAna_Ht->Fill(ht); // Ht--sum Et of all objects
		Hist.fAna_Ht_MetSig->Fill(metsig, ht); // Ht--sum Et of all objects
		Hist.fAna_Njet20->Fill(jetstuff.myNjet_th20); // Njet(Et>20)
		Hist.fAna_Njet25->Fill(jetstuff.myNjet_th25); // Njet(Et>25)
		Hist.fAna_Njet30->Fill(jetstuff.myNjet_th30); // Njet(Et>30)
		Hist.fAna_Njet35->Fill(jetstuff.myNjet_th35); // Njet(Et>35)

		float _dphi=fabs(TVector2::Phi_mpi_pi(Met.Phi()-vJet.at(0).Phi()));
		Hist.fAna_Etjet->Fill(vJet.at(0).Pt()); // Et of 1st jet 
		Hist.fAna_Etjet_MetSig->Fill(metsig, vJet.at(0).Pt()); // Et of 1st jet 
		//I am doing this trick to get rid the warning of from Eta() when Pt=0
		if (vJet.at(0).Pt() !=0)
		{
			Hist.fAna_Jet1Eta->Fill(vJet.at(0).Eta()); // Eta of 1st jet with Et>15
			Hist.fAna_Jet1Eta_MetSig->Fill(metsig, vJet.at(0).Eta()); // Eta of 1st jet with Et>15
		} else
		{
			Hist.fAna_Jet1Eta->Fill(fMinisculeEta);
			Hist.fAna_Jet1Eta_MetSig->Fill(metsig, fMinisculeEta);
		}

		//sub leading jet
		if (vJet.size() >1)
		{
			if (vJet.at(1).Pt() > 0.001)  // a jet that did not match to EM obj and does not have a neglegible value
			{
				float fJet2Eta = fabs(vJet.at(1).Eta());
				if (fJet2Eta > MinJetEta
					 && fJet2Eta < MaxJetEta)
				{
					Hist.fAna_Etjet2->Fill(vJet.at(1).Pt()); // Et of 2nd jet with Et>0
					Hist.fAna_Jet2Eta->Fill(vJet.at(1).Eta()); // Eta of 2nd jet with Et>0
					Hist.fAna_Etjet2_MetSig->Fill(metsig, vJet.at(1).Pt()); // Et of 2nd jet with Et>0
					Hist.fAna_Jet2Eta_MetSig->Fill(metsig, vJet.at(1).Eta()); // Eta of 2nd jet with Et>0
				}
			}
		}

		Hist.fAna_dPhi3->Fill(_dphi); // dPhi(met-jet1)
		Hist.fAna_dPhi3_MetSig->Fill(metsig, _dphi); // dPhi(met-jet1)
		if(jetstuff.myNjet_th20>0) Hist.fAna_dPhi4->Fill(_dphi); // dPhi(met-jet1)
		if(jetstuff.myNjet_th25>0) Hist.fAna_dPhi5->Fill(_dphi); // dPhi(met-jet1)
		if(jetstuff.myNjet_th30>0) Hist.fAna_dPhi6->Fill(_dphi); // dPhi(met-jet1)
		if(jetstuff.myNjet_th35>0) Hist.fAna_dPhi7->Fill(_dphi); // dPhi(met-jet1)


		if(iNjet15>1) Hist.fAna_Mjj->Fill((vJet.at(0)+vJet.at(1)).M()); // Mjj (of two highest Et jets) if Njet15>=2

		Hist.fAna_Nem->Fill(miscstuff.myCorrElectron.size()+miscstuff.myCorrPhoton.size()); // number of EM objects: ele+pho

		if(miscstuff.myCorrPhoton.size()>=1) // assuming these are either loose or tight photons only and in Et sorted
		{
			//find lead tight photon
			int iLeadPhoIndex = -1;
			for (unsigned int i = 0; i < miscstuff.tightId.size(); ++i)
			{
				if (miscstuff.tightId.at(i) == 0)
				{
					iLeadPhoIndex = i;		//need only the leading tight photon
					break;
				}
			}

			assert (iLeadPhoIndex>=0 && "A tight photon is not found!");
			
			Hist.fAna_Et1->Fill(miscstuff.myCorrPhoton.at(iLeadPhoIndex).Pt()); 	// Et of first EM object
			Hist.fAna_Et1_MetSig->Fill(metsig, miscstuff.myCorrPhoton.at(iLeadPhoIndex).Pt()); 	// Et of first EM object
			Hist.fAna_Pho1Eta->Fill(miscstuff.myCorrPhoton.at(iLeadPhoIndex).Eta()); 	// Eta of first EM object
			float _dphi_MetEm=fabs(TVector2::Phi_mpi_pi(miscstuff.myCorrPhoton.at(iLeadPhoIndex).Phi()-Met.Phi()));
			Hist.fAna_dPhi1->Fill(_dphi_MetEm); // dPhi(met-em1)

			if (vJet.size()>=2)
			{
				if(vJet.at(0).Pt()>15.0 && vJet.at(1).Pt()>5.0 && iLeadPhoIndex>=0)
				{
					float j1phomass = (miscstuff.myCorrPhoton.at(iLeadPhoIndex)+vJet.at(0)).M();	
					float j2phomass = (miscstuff.myCorrPhoton.at(iLeadPhoIndex)+vJet.at(1)).M();	
					Hist.fAna_Dalitz->Fill(j1phomass,j2phomass);
				}
			}

			//for debuggin only
			if (fHeaderBlock->McFlag())
			{

				TLorentzVector hepgPhoVec = FindMatchingHEPGPar(miscstuff.myCorrPhoton.at(iLeadPhoIndex),22,3,0.1);
				float PtRatio =  miscstuff.myCorrPhoton.at(iLeadPhoIndex).Pt()/ hepgPhoVec.Pt() - 1;
				float JetPhoPtRatio = vJet.at(0).Pt()/hepgPhoVec.Pt() - 1;

				//photon energy resolution
				float fERatio =  miscstuff.myCorrPhoton.at(iLeadPhoIndex).E()/ hepgPhoVec.E() - 1;
				Hist.hDebug_PhoRes->Fill(fERatio);

				if (_dphi<0.4)
				{
					Hist.hDebug_PtRatio_JMetdelphiLow->Fill(PtRatio);
					Hist.hDebug_PhoJetPtRatio_JMetdelphiLow->Fill(JetPhoPtRatio);
				}
				if (_dphi>1.3 && _dphi<1.7)
				{
					Hist.hDebug_PtRatio_JMetdelphiMid->Fill(PtRatio);
					Hist.hDebug_PhoJetPtRatio_JMetdelphiMid->Fill(JetPhoPtRatio);
				}
				if (_dphi>2.8)
				{
					Hist.hDebug_PtRatio_JMetdelphiTop->Fill(PtRatio);
					Hist.hDebug_PhoJetPtRatio_JMetdelphiTop->Fill(JetPhoPtRatio);
					if(jetstuff.myNjet_th5>=2)
					{
						TLorentzVector j1vec = vJet.at(0);
						TLorentzVector j2vec = vJet.at(1);
						float dPhij1j2 = j1vec.DeltaPhi(j2vec); 
						if (fabs(dPhij1j2) <0.4)
						{
							if (fHeaderBlock->McFlag())
							{
								TLorentzVector j1HadJet = FindMatchingHADJet(j1vec,0.4);
								TLorentzVector j2HadJet = FindMatchingHADJet(j2vec,0.4);
								if (j1HadJet.Pt() == j2HadJet.Pt())
								{
									if (PrintLevel() > PRN_WARN)
									{
										std::cout << __FUNCTION__ << "::" << __LINE__ 
											<< ":: Two det jets matched to the same had jet! Run, Event="
											<< fHeaderBlock->RunNumber() << ", "
											<< fHeaderBlock->EventNumber() << std::endl;
									}
								}

								float j1ratio = j1vec.Pt()/j1HadJet.Pt() - 1;
								float j2ratio = j2vec.Pt()/j2HadJet.Pt() - 1;
								Hist.hDebug_Jet1PtRatio_dphi4->Fill(j1ratio);
								Hist.hDebug_Jet2PtRatio_dphi4->Fill(j2ratio);
							}
						}
					}
				}

			}
			//this is plot the leading jet after removing only the leading photon

			std::vector<TLorentzVector> tlvObj;
			
			assert(miscstuff.myCorrPhoton.size() == miscstuff.tightId.size()
						&& "JetFilterV2::FillAnalysisHists:: vector size mismatch! corPhoSize!= tightPho.size()");

			//find tight photons
			for (unsigned int i = 0; i < miscstuff.tightId.size(); ++i)
			{
				if (miscstuff.tightId.at(i) == 0)
				{
					tlvObj.push_back(miscstuff.myCorrPhoton.at(i));
				}
			}

			float jetPt = vFinalJets.at(0).Pt(); //lead jet
			if (tlvObj.size()>1) //photons >1
			{
				if (tlvObj.at(1).Pt()> jetPt) jetPt = tlvObj.at(1).Pt();
			}

			Hist.hDebug_LeadJet->Fill(jetPt);

			if (PrintLevel()>10)
			{
				if (tlvObj.size()>1)
				{
					std::cout << "======= " << __FUNCTION__ << " ==== DEBUG INFO " << std::endl;
					std::cout << "DiPhoton event: g1, g2 Pt =" 
						<< miscstuff.myCorrPhoton.at(0).Pt() << " , " 
						<< miscstuff.myCorrPhoton.at(1).Pt()  << std::endl;
					std::cout << "Lead jet Pt= " <<  vFinalJets.at(0).Pt() << std::endl;
					std::cout << "choosen jet Pt= " <<  jetPt << std::endl;

				}
			}


			//end debugging stuff


			if(jetstuff.myNjet_th15>0)
			{
				Hist.fAna_Mej->Fill((miscstuff.myCorrPhoton.at(0)+vJet.at(0)).M()); // M(em-jet1) if Njet15>=1
			}
			for(unsigned int i=2; i<miscstuff.myCorrPhoton.size(); i++)
			{
				Hist.fAna_Mextra->Fill((miscstuff.myCorrPhoton.at(0)+miscstuff.myCorrPhoton.at(i)).M()); // M(em1-emN) & M(em2-emN) if Nem>=3
				Hist.fAna_Etem->Fill(miscstuff.myCorrPhoton.at(i).Pt()); // Et of extra EM objetcs
			}
			for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
			{
				Hist.fAna_Mextra->Fill((miscstuff.myCorrPhoton.at(0)+miscstuff.myCorrElectron.at(i)).M()); // M(em1-emN) & M(em2-emN) if Nem>=3
				Hist.fAna_Etem->Fill(miscstuff.myCorrElectron.at(i).Pt()); // Et of extra EM objetcs
			}
		} //photons>0

	}	// metsigcut

	return;
}



//____________________________________________________________________ filling histograms for pho-met study
void JetFilterModuleV2::FillPhoMetStudyHistograms(PhoMetStudyHisto_t& Hist, TVector2 Met, CommonStuff miscstuff, double zvx) {
	for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
	{
		double _dphi=fabs(TVector2::Phi_mpi_pi(Met.Phi()-miscstuff.myCorrPhoton.at(i).Phi()));
		double p0, p1;
		if(fabs(miscstuff.myPhoEtaDet.at(i))<1.1)
		{
			p0=0.135*0.135;
			p1=0.015*0.015;
		}
		else
		{
			p0=0.16*0.16;
			p1=0.01*0.01;
		}
		double metsig=Met.Mod()*cos(_dphi)/sqrt(miscstuff.myCorrPhoton.at(i).Pt()*p0+miscstuff.myCorrPhoton.at(i).Pt()*miscstuff.myCorrPhoton.at(i).Pt()*p1);
		Hist.fPhoMet_eta->Fill(fabs(miscstuff.myPhoEtaDet.at(i)));
		Hist.fPhoMet_xces->Fill(fabs(miscstuff.myPhoXces.at(i)));
		Hist.fPhoMet_zces->Fill(fabs(miscstuff.myPhoZces.at(i)));

		if(zvx*miscstuff.myPhoEtaDet.at(i)>0.0) Hist.fPhoMet_signeta->Fill(fabs(miscstuff.myPhoEtaDet.at(i)));
		else Hist.fPhoMet_signeta->Fill(-1.0*fabs(miscstuff.myPhoEtaDet.at(i)));
		if(zvx*miscstuff.myPhoZces.at(i)>0.0) Hist.fPhoMet_signzces->Fill(fabs(miscstuff.myPhoZces.at(i)));
		else Hist.fPhoMet_signzces->Fill(-1.0*fabs(miscstuff.myPhoZces.at(i)));

		if(Met.Mod()>0.0 && _dphi<0.2)
		{
			Hist.fPhoMet_dPhi02_eta->Fill(fabs(miscstuff.myPhoEtaDet.at(i)));
			Hist.fPhoMet_dPhi02_xces->Fill(fabs(miscstuff.myPhoXces.at(i)));
			Hist.fPhoMet_dPhi02_zces->Fill(fabs(miscstuff.myPhoZces.at(i)));
			if(metsig>5.0)
			{
				Hist.fPhoMet_MetSig5dPhi02_eta->Fill(fabs(miscstuff.myPhoEtaDet.at(i)));
				Hist.fPhoMet_MetSig5dPhi02_xces->Fill(fabs(miscstuff.myPhoXces.at(i)));
				Hist.fPhoMet_MetSig5dPhi02_zces->Fill(fabs(miscstuff.myPhoZces.at(i)));
				if(zvx*miscstuff.myPhoEtaDet.at(i)>0.0) Hist.fPhoMet_MetSig5dPhi02_signeta->Fill(fabs(miscstuff.myPhoEtaDet.at(i)));
				else Hist.fPhoMet_MetSig5dPhi02_signeta->Fill(-1.0*fabs(miscstuff.myPhoEtaDet.at(i)));
				if(zvx*miscstuff.myPhoZces.at(i)>0.0) Hist.fPhoMet_MetSig5dPhi02_signzces->Fill(fabs(miscstuff.myPhoZces.at(i)));
				else Hist.fPhoMet_MetSig5dPhi02_signzces->Fill(-1.0*fabs(miscstuff.myPhoZces.at(i)));
			}
			if(_dphi<0.1)
			{
				Hist.fPhoMet_dPhi01_eta->Fill(fabs(miscstuff.myPhoEtaDet.at(i)));
				Hist.fPhoMet_dPhi01_xces->Fill(fabs(miscstuff.myPhoXces.at(i)));
				Hist.fPhoMet_dPhi01_zces->Fill(fabs(miscstuff.myPhoZces.at(i)));
				if(metsig>5.0)
				{
					Hist.fPhoMet_MetSig5dPhi01_eta->Fill(fabs(miscstuff.myPhoEtaDet.at(i)));
					Hist.fPhoMet_MetSig5dPhi01_xces->Fill(fabs(miscstuff.myPhoXces.at(i)));
					Hist.fPhoMet_MetSig5dPhi01_zces->Fill(fabs(miscstuff.myPhoZces.at(i)));
				}
			}
		}
	}

	for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
	{
		double _dphi=fabs(TVector2::Phi_mpi_pi(Met.Phi()-miscstuff.myCorrElectron.at(i).Phi()));
		double p0, p1;
		if(fabs(miscstuff.myEleEtaDet.at(i))<1.1)
		{
			p0=0.135*0.135;
			p1=0.015*0.015;
		}
		else
		{
			p0=0.16*0.16;
			p1=0.01*0.01;
		}
		double metsig=Met.Mod()*cos(_dphi)/sqrt(miscstuff.myCorrElectron.at(i).Pt()*p0+miscstuff.myCorrElectron.at(i).Pt()*miscstuff.myCorrElectron.at(i).Pt()*p1);
		Hist.fEleMet_eta->Fill(fabs(miscstuff.myEleEtaDet.at(i)));
		Hist.fEleMet_phi->Fill(fabs(miscstuff.myElePhiDet.at(i)));

		if(Met.Mod()>0.0 && _dphi<0.2)
		{
			Hist.fEleMet_dPhi02_eta->Fill(fabs(miscstuff.myEleEtaDet.at(i)));
			Hist.fEleMet_dPhi02_phi->Fill(fabs(miscstuff.myElePhiDet.at(i)));
			if(metsig>5.0)
			{
				Hist.fEleMet_MetSig5dPhi02_eta->Fill(fabs(miscstuff.myEleEtaDet.at(i)));
				Hist.fEleMet_MetSig5dPhi02_phi->Fill(fabs(miscstuff.myElePhiDet.at(i)));
			}
			if(_dphi<0.1)
			{
				Hist.fEleMet_dPhi01_eta->Fill(fabs(miscstuff.myEleEtaDet.at(i)));
				Hist.fEleMet_dPhi01_phi->Fill(fabs(miscstuff.myElePhiDet.at(i)));
				if(metsig>5.0)
				{
					Hist.fEleMet_MetSig5dPhi01_eta->Fill(fabs(miscstuff.myEleEtaDet.at(i)));
					Hist.fEleMet_MetSig5dPhi01_phi->Fill(fabs(miscstuff.myElePhiDet.at(i)));
				}
			}
		}
	}

	return;
}

//_____________________________________________ finalize histograms for pho-met study
void JetFilterModuleV2::DoFinalPhoMetHisto(PhoMetStudyHisto_t& Hist) {

	Hist.fPhoMet_Fr_dPhi02_eta->Sumw2();
	Hist.fPhoMet_Fr_dPhi02_xces->Sumw2();
	Hist.fPhoMet_Fr_dPhi02_zces->Sumw2();
	Hist.fPhoMet_Fr_dPhi01_eta->Sumw2();
	Hist.fPhoMet_Fr_dPhi01_xces->Sumw2();
	Hist.fPhoMet_Fr_dPhi01_zces->Sumw2();

	Hist.fPhoMet_Fr_MetSig5dPhi02_eta->Sumw2();
	Hist.fPhoMet_Fr_MetSig5dPhi02_xces->Sumw2();
	Hist.fPhoMet_Fr_MetSig5dPhi02_zces->Sumw2();
	Hist.fPhoMet_Fr_MetSig5dPhi01_eta->Sumw2();
	Hist.fPhoMet_Fr_MetSig5dPhi01_xces->Sumw2();
	Hist.fPhoMet_Fr_MetSig5dPhi01_zces->Sumw2();

	Hist.fPhoMet_Fr_MetSig5dPhi02_signeta->Sumw2();
	Hist.fPhoMet_Fr_MetSig5dPhi02_signzces->Sumw2();

	Hist.fPhoMet_Fr_dPhi02_eta->Divide(Hist.fPhoMet_dPhi02_eta,Hist.fPhoMet_eta);
	Hist.fPhoMet_Fr_dPhi02_xces->Divide(Hist.fPhoMet_dPhi02_xces,Hist.fPhoMet_xces);
	Hist.fPhoMet_Fr_dPhi02_zces->Divide(Hist.fPhoMet_dPhi02_zces,Hist.fPhoMet_zces);
	Hist.fPhoMet_Fr_dPhi01_eta->Divide(Hist.fPhoMet_dPhi01_eta,Hist.fPhoMet_eta);
	Hist.fPhoMet_Fr_dPhi01_xces->Divide(Hist.fPhoMet_dPhi01_xces,Hist.fPhoMet_xces);
	Hist.fPhoMet_Fr_dPhi01_zces->Divide(Hist.fPhoMet_dPhi01_zces,Hist.fPhoMet_zces);

	Hist.fPhoMet_Fr_MetSig5dPhi02_eta->Divide(Hist.fPhoMet_MetSig5dPhi02_eta,Hist.fPhoMet_eta);
	Hist.fPhoMet_Fr_MetSig5dPhi02_xces->Divide(Hist.fPhoMet_MetSig5dPhi02_xces,Hist.fPhoMet_xces);
	Hist.fPhoMet_Fr_MetSig5dPhi02_zces->Divide(Hist.fPhoMet_MetSig5dPhi02_zces,Hist.fPhoMet_zces);
	Hist.fPhoMet_Fr_MetSig5dPhi01_eta->Divide(Hist.fPhoMet_MetSig5dPhi01_eta,Hist.fPhoMet_eta);
	Hist.fPhoMet_Fr_MetSig5dPhi01_xces->Divide(Hist.fPhoMet_MetSig5dPhi01_xces,Hist.fPhoMet_xces);
	Hist.fPhoMet_Fr_MetSig5dPhi01_zces->Divide(Hist.fPhoMet_MetSig5dPhi01_zces,Hist.fPhoMet_zces);

	Hist.fPhoMet_Fr_MetSig5dPhi02_signeta->Divide(Hist.fPhoMet_MetSig5dPhi02_signeta,Hist.fPhoMet_signeta);
	Hist.fPhoMet_Fr_MetSig5dPhi02_signzces->Divide(Hist.fPhoMet_MetSig5dPhi02_signzces,Hist.fPhoMet_signzces);


	//-------------------- correlation between electron and met ---------------------

	Hist.fEleMet_Fr_dPhi02_eta->Sumw2();
	Hist.fEleMet_Fr_dPhi02_phi->Sumw2();
	Hist.fEleMet_Fr_dPhi01_eta->Sumw2();
	Hist.fEleMet_Fr_dPhi01_phi->Sumw2();

	Hist.fEleMet_Fr_MetSig5dPhi02_eta->Sumw2();
	Hist.fEleMet_Fr_MetSig5dPhi02_phi->Sumw2();
	Hist.fEleMet_Fr_MetSig5dPhi01_eta->Sumw2();
	Hist.fEleMet_Fr_MetSig5dPhi01_phi->Sumw2();

	Hist.fEleMet_Fr_dPhi02_eta->Divide(Hist.fEleMet_dPhi02_eta,Hist.fEleMet_eta);
	Hist.fEleMet_Fr_dPhi02_phi->Divide(Hist.fEleMet_dPhi02_phi,Hist.fEleMet_phi);
	Hist.fEleMet_Fr_dPhi01_eta->Divide(Hist.fEleMet_dPhi01_eta,Hist.fEleMet_eta);
	Hist.fEleMet_Fr_dPhi01_phi->Divide(Hist.fEleMet_dPhi01_phi,Hist.fEleMet_phi);

	Hist.fEleMet_Fr_MetSig5dPhi02_eta->Divide(Hist.fEleMet_MetSig5dPhi02_eta,Hist.fEleMet_eta);
	Hist.fEleMet_Fr_MetSig5dPhi02_phi->Divide(Hist.fEleMet_MetSig5dPhi02_phi,Hist.fEleMet_phi);
	Hist.fEleMet_Fr_MetSig5dPhi01_eta->Divide(Hist.fEleMet_MetSig5dPhi01_eta,Hist.fEleMet_eta);
	Hist.fEleMet_Fr_MetSig5dPhi01_phi->Divide(Hist.fEleMet_MetSig5dPhi01_phi,Hist.fEleMet_phi);

	return;
}

//__________________________________ Filling histograms for Met cleanup studies
void JetFilterModuleV2::FillMetCleanupHistograms(MetCleanupHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, int metcode) {

	double dphi;
	double ratio;
	double etadet;
	double ptcut=0.0;
	TVector2 metvec(0.0,0.0);
	TVector2 metvec_raw(0.0,0.0);

	double dphi_3=10.0;
	double dphi_5=10.0;
	double dphi_10=10.0;
	double dphi_1st=10.0;

	if(metcode==1)
	{
		metvec.Set(jetstuff.myMETcorr_th5.Px(),jetstuff.myMETcorr_th5.Py());
		ptcut=5.0;
	}
	if(metcode==2)
	{
		metvec.Set(jetstuff.myMETcorr_th10.Px(),jetstuff.myMETcorr_th10.Py());
		ptcut=10.0;
	}
	if(metcode==3)
	{
		metvec.Set(jetstuff.myMETcorr_th15.Px(),jetstuff.myMETcorr_th15.Py());
		ptcut=15.0;
	}
	if(metcode==4)
	{
		metvec.Set(jetstuff.myMETcorr_th20.Px(),jetstuff.myMETcorr_th20.Py());
		ptcut=20.0;
	}
	if(metcode<1 || metcode>4)
	{
		metvec.Set(miscstuff.myMET_raw.Px(),miscstuff.myMET_raw.Py());
		ptcut=0.0;
	}
	metvec_raw.Set(miscstuff.myMET_raw.Px(),miscstuff.myMET_raw.Py());

	if(metvec.Mod()>0.0)
	{
		for(unsigned int i=0; i<miscstuff.myRawPhoton.size(); i++)
		{
			dphi=fabs(TVector2::Phi_mpi_pi(metvec.Phi()-miscstuff.myRawPhoton.at(i).Phi()));
			ratio=metvec.Mod()/miscstuff.myRawPhoton.at(i).Pt();
			etadet=miscstuff.myPhoEtaDet.at(i);
			Hist.fCleanupDelPhiVsEtaDet_em->Fill(etadet,dphi);
			Hist.fCleanupDelPhi_metem->Fill(dphi);
			Hist.fCleanupMetEtem->Fill(ratio);
			Hist.fCleanupMetEtemVsDelPhi->Fill(dphi,ratio);
		}
		for(unsigned int i=0; i<miscstuff.myRawElectron.size(); i++)
		{
			dphi=fabs(TVector2::Phi_mpi_pi(metvec.Phi()-miscstuff.myRawElectron.at(i).Phi()));
			ratio=metvec.Mod()/miscstuff.myRawElectron.at(i).Pt();
			etadet=miscstuff.myEleEtaDet.at(i);
			Hist.fCleanupDelPhiVsEtaDet_em->Fill(etadet,dphi);
			Hist.fCleanupDelPhi_metem->Fill(dphi);
			Hist.fCleanupMetEtem->Fill(ratio);
			Hist.fCleanupMetEtemVsDelPhi->Fill(dphi,ratio);
		}

		double maxEt=15.0; // new part as of 04/03/07
		double minDelPhi=10.0;  // new part as of 04/03/07
		double minDelPhi_1st=10.0;  // new part as of 04/03/07
		for(unsigned int i=0; i<jetstuff.Jet_lev6_noEMobj.size(); i++)
		{
			dphi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(metvec.Phi())-TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj.at(i).Phi())));
			etadet=jetstuff.EtaDetCorr.at(i);
			if((fabs(etadet)>0.2 && fabs(etadet)<0.9)
					|| (fabs(etadet)>1.3 && fabs(etadet)<2.6))
			{
				if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>maxEt && fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<3.0)  // new part as of 04/03/07
				{
					maxEt=jetstuff.Jet_lev6_noEMobj.at(i).Pt();
					minDelPhi_1st=dphi;
				}
				if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>15.0
						&& fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<3.0
						&& dphi<minDelPhi) minDelPhi=dphi;  // new part as of 04/03/07
				if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>5.0
						&& fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<3.0) Hist.fCleanupDelPhi_metjet5->Fill(dphi);  // new part as of 04/03/07
				if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>15.0
						&& fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<3.0)   // new part as of 04/03/07
				{
					Hist.fCleanupDelPhi_metjet15->Fill(dphi);
					if(metvec.Mod()>10.0) Hist.fCleanupDelPhi_met10jet15->Fill(dphi);
				}
				if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>25
						&& fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<3.0) Hist.fCleanupDelPhi_metjet25->Fill(dphi);  // new part as of 04/03/07
			}
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>ptcut)
			{
				ratio=metvec.Mod()/jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				Hist.fCleanupDelPhiVsEtaDet_jet->Fill(etadet,dphi);
				int j=1;
				if(fabs(etadet)>2.6) j=2;
				if(fabs(etadet)<0.2 || (fabs(etadet)>0.9 && fabs(etadet)<1.3)) j=0;
				Hist.fCleanupDelPhi_metjet[j]->Fill(dphi);
				Hist.fCleanupMetEtjet[j]->Fill(ratio);
				Hist.fCleanupMetEtjetVsDelPhi[j]->Fill(dphi,ratio);
			}
		}
		if(minDelPhi_1st<=TMath::Pi())  // new part as of 04/03/07
		{
			Hist.fCleanupDelPhi_metjet1st15->Fill(minDelPhi_1st);
			if(metvec.Mod()>10.0) Hist.fCleanupDelPhi_met10jet1st15->Fill(minDelPhi_1st);
		}
		if(minDelPhi<=TMath::Pi())  // new part as of 04/03/07
		{
			Hist.fCleanupDelPhi_metjet15_dPhiMin->Fill(minDelPhi);
			if(metvec.Mod()>10.0) Hist.fCleanupDelPhi_met10jet15_dPhiMin->Fill(minDelPhi);
		}
	}

	if(metvec_raw.Mod()>0.0)
	{
		double rawEtmax=0.0;
		for(unsigned int i=0; i<jetstuff.Jet_raw_noEMobj.size(); i++)
		{
			double _dphi=15.0; // dummy value
			if(jetstuff.Jet_raw_noEMobj.at(i).Pt()>3.0)
			{
				_dphi=fabs(TVector2::Phi_mpi_pi(metvec_raw.Phi()-jetstuff.Jet_raw_noEMobj.at(i).Phi()));
				if(jetstuff.Jet_raw_noEMobj.at(i).Pt()>rawEtmax)
				{
					rawEtmax=jetstuff.Jet_raw_noEMobj.at(i).Pt();
					dphi_1st=_dphi;
				}
				if(_dphi<dphi_3) dphi_3=_dphi;
				if(jetstuff.Jet_raw_noEMobj.at(i).Pt()>5.0 && _dphi<dphi_5) dphi_5=_dphi;
				if(jetstuff.Jet_raw_noEMobj.at(i).Pt()>10.0 && _dphi<dphi_10) dphi_10=_dphi;
			}
		}
		if(dphi_3<3.2) Hist.fCleanupDelPhi_RawMetRawJet3->Fill(dphi_3);
		if(dphi_5<3.2) Hist.fCleanupDelPhi_RawMetRawJet5->Fill(dphi_5);
		if(dphi_10<3.2) Hist.fCleanupDelPhi_RawMetRawJet10->Fill(dphi_10);
		if(dphi_1st<3.2) Hist.fCleanupDelPhi_RawMet1stRawJet->Fill(dphi_1st);
	}

	return;
}

//__________________________________ Filling histograms for GenMet cleanup studies
void JetFilterModuleV2::FillGenMetCleanupHistograms(MetCleanupHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff,
		TVector2 genMet, std::vector<TLorentzVector> vec) {

	double dphi;
	double ratio;
	double etadet;
	double ptcut=15.0;
	if(genMet.Mod()>0.0 && genMet.Mod()<2000.0)
	{
		for(unsigned int i=0; i<miscstuff.myRawPhoton.size(); i++)
		{
			dphi=fabs(TVector2::Phi_mpi_pi(genMet.Phi()-miscstuff.myRawPhoton.at(i).Phi()));
			ratio=genMet.Mod()/miscstuff.myRawPhoton.at(i).Pt();
			etadet=miscstuff.myPhoEtaDet.at(i);
			Hist.fCleanupDelPhiVsEtaDet_em->Fill(etadet,dphi);
			Hist.fCleanupDelPhi_metem->Fill(dphi);
			Hist.fCleanupMetEtem->Fill(ratio);
			Hist.fCleanupMetEtemVsDelPhi->Fill(dphi,ratio);
		}
		for(unsigned int i=0; i<miscstuff.myRawElectron.size(); i++)
		{
			dphi=fabs(TVector2::Phi_mpi_pi(genMet.Phi()-miscstuff.myRawElectron.at(i).Phi()));
			ratio=genMet.Mod()/miscstuff.myRawElectron.at(i).Pt();
			etadet=miscstuff.myEleEtaDet.at(i);
			Hist.fCleanupDelPhiVsEtaDet_em->Fill(etadet,dphi);
			Hist.fCleanupDelPhi_metem->Fill(dphi);
			Hist.fCleanupMetEtem->Fill(ratio);
			Hist.fCleanupMetEtemVsDelPhi->Fill(dphi,ratio);
		}

		double maxEt=15.0; // new part as of 04/03/07
		double minDelPhi=10.0;  // new part as of 04/03/07
		double minDelPhi_1st=10.0;  // new part as of 04/03/07
		for(unsigned int i=0; i<vec.size(); i++)
		{
			dphi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(genMet.Phi())-TVector2::Phi_0_2pi(vec.at(i).Phi())));
			if(vec.at(i).Pt()>0.0) ratio=genMet.Mod()/vec.at(i).Pt();
			else ratio=9999.0; // dummy value
			etadet=jetstuff.EtaDetCorr.at(i);
			if((fabs(etadet)>0.2 && fabs(etadet)<0.9)
					|| (fabs(etadet)>1.3 && fabs(etadet)<2.6))
			{
				if(vec.at(i).Pt()>maxEt && fabs(vec.at(i).Eta())<3.0)  // new part as of 04/03/07
				{
					maxEt=vec.at(i).Pt();
					minDelPhi_1st=dphi;
				}
				if(vec.at(i).Pt()>15.0
						&& fabs(vec.at(i).Eta())<3.0
						&& dphi<minDelPhi) minDelPhi=dphi;  // new part as of 04/03/07
				if(vec.at(i).Pt()>5.0 && fabs(vec.at(i).Eta())<3.0) Hist.fCleanupDelPhi_metjet5->Fill(dphi);  // new part as of 04/03/07
				if(vec.at(i).Pt()>15.0 && fabs(vec.at(i).Eta())<3.0)   // new part as of 04/03/07
				{
					Hist.fCleanupDelPhi_metjet15->Fill(dphi);
					if(genMet.Mod()>10.0) Hist.fCleanupDelPhi_met10jet15->Fill(dphi);
				}
				if(vec.at(i).Pt()>25 && fabs(vec.at(i).Eta())<3.0) Hist.fCleanupDelPhi_metjet25->Fill(dphi);  // new part as of 04/03/07
			}
			if(vec.at(i).Pt()>ptcut)
			{
				Hist.fCleanupDelPhiVsEtaDet_jet->Fill(etadet,dphi);
				int j=1;
				if(fabs(etadet)>2.6) j=2;
				if(fabs(etadet)<0.2 || (fabs(etadet)>0.9 && fabs(etadet)<1.3)) j=0;
				Hist.fCleanupDelPhi_metjet[j]->Fill(dphi);
				Hist.fCleanupMetEtjet[j]->Fill(ratio);
				Hist.fCleanupMetEtjetVsDelPhi[j]->Fill(dphi,ratio);
			}
		}
		if(minDelPhi_1st<=TMath::Pi())  // new part as of 04/03/07
		{
			Hist.fCleanupDelPhi_metjet1st15->Fill(minDelPhi_1st);
			if(genMet.Mod()>10.0) Hist.fCleanupDelPhi_met10jet1st15->Fill(minDelPhi_1st);
		}
		if(minDelPhi<=TMath::Pi())  // new part as of 04/03/07
		{
			Hist.fCleanupDelPhi_metjet15_dPhiMin->Fill(minDelPhi);
			if(genMet.Mod()>10.0) Hist.fCleanupDelPhi_met10jet15_dPhiMin->Fill(minDelPhi);
		}
	}
	return;
}

//_____________________________ calculates dPhi for metsig.vs.dPhi studies
void JetFilterModuleV2::CalculateMetDelPhi(MetSigDelPhiStuff &msdp_stuff,JetStuff jetstuff,CommonStuff miscstuff,TVector2 MetVec) {
	msdp_stuff.eta_det_em=10.0; // detEta of closest EM object
	msdp_stuff.eta_det_jet=10.0; // detEta of closest jet
	msdp_stuff.xces_em=25.0; // X_ces of closest EM
	for(int ij=0; ij<5; ij++)
	{
		msdp_stuff.dPhi_det[ij]=10.0; // initial dummy value
	}

	for(unsigned int ij=0; ij<miscstuff.myCorrPhoton.size(); ij++)
	{
		double _dphi_d=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(miscstuff.myCorrPhoton[ij].Phi())-TVector2::Phi_0_2pi(MetVec.Phi())));
		if(ij==0) msdp_stuff.dPhi_det[3]=_dphi_d;
		if(_dphi_d<msdp_stuff.dPhi_det[1])
		{
			msdp_stuff.dPhi_det[1]=_dphi_d;
			msdp_stuff.eta_det_em=miscstuff.myPhoEtaDet[ij];
			msdp_stuff.xces_em=miscstuff.myPhoXces[ij];
		}
	}
	for(unsigned int ij=0; ij<miscstuff.myCorrElectron.size(); ij++)
	{
		double _dphi_d=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(miscstuff.myCorrElectron[ij].Phi())-TVector2::Phi_0_2pi(MetVec.Phi())));
		if(ij==0)
		{
			if(miscstuff.myCorrPhoton.size()>0) // make sure there are photons before next step
			{
				if(miscstuff.myCorrElectron[0].Pt()>miscstuff.myCorrPhoton[0].Pt()) msdp_stuff.dPhi_det[3]=_dphi_d;
			}
			else msdp_stuff.dPhi_det[3]=_dphi_d;
		}
		if(_dphi_d<msdp_stuff.dPhi_det[1])
		{
			msdp_stuff.dPhi_det[1]=_dphi_d;
			msdp_stuff.eta_det_em=miscstuff.myEleEtaDet[ij];
			msdp_stuff.xces_em=miscstuff.myEleXces[ij];
		}
	}
	double max_jet_Pt=0.0;
	for(unsigned int ij=0; ij<jetstuff.Jet_lev6_noEMobj.size(); ij++)
	{
		if(jetstuff.Jet_lev6_noEMobj[ij].Pt()>3.0 && (jetstuff.Npho_match[ij]+jetstuff.Nele_match[ij])==0)
		{
			double _dphi_d=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj[ij].Phi())-TVector2::Phi_0_2pi(MetVec.Phi())));
			if(jetstuff.Jet_lev6_noEMobj[ij].Pt()>max_jet_Pt)
			{
				max_jet_Pt=jetstuff.Jet_lev6_noEMobj[ij].Pt();
				msdp_stuff.dPhi_det[4]=_dphi_d;
			}
			if(_dphi_d<msdp_stuff.dPhi_det[2])
			{
				msdp_stuff.dPhi_det[2]=_dphi_d;
				msdp_stuff.eta_det_jet=jetstuff.EtaDetCorr[ij];
			}
		}
	}
	if(msdp_stuff.dPhi_det[2]<msdp_stuff.dPhi_det[1]) msdp_stuff.dPhi_det[0]=msdp_stuff.dPhi_det[2];
	else msdp_stuff.dPhi_det[0]=msdp_stuff.dPhi_det[1];

	return;
}

//_______________________ returns difference between two histograms in bin=bin_ind
double JetFilterModuleV2::HistoBinDiff(TH1F* h,TH1F* h1,TH1F* h2,int bin_ind)
{
	double value=0.0;
	double v0=h->GetBinContent(bin_ind+1);
	double err=h->GetBinError(bin_ind+1);
	double v1=h1->GetBinContent(bin_ind+1);
	double v2=h2->GetBinContent(bin_ind+1);
	value=(fabs(v0-v1) > fabs(v0-v2)) ? sqrt(err*err+(v0-v1)*(v0-v1)) : sqrt(err*err+(v0-v2)*(v0-v2));
	return value;
}

//_______________________ finalize analysis histograms
//This calculates and sets the errors for the final prediction hists in fAna_bckg - sam
//_________________________________________________________________________________
void JetFilterModuleV2::DoFinalAnalysisHisto(AnalysisHisto_t& HistF,AnalysisHisto_t Hist1,AnalysisHisto_t Hist2) {

	for(int i=0; i<27; i++)
	{
		int nbins=0;
		if(i==0) nbins=HistF.fAna_MetAll->GetNbinsX();
		if(i==1) nbins=HistF.fAna_MetSig->GetNbinsX();
		if(i==2) nbins=HistF.fAna_Met->GetNbinsX();
		if(i==3) nbins=HistF.fAna_M->GetNbinsX();
		if(i==4) nbins=HistF.fAna_dPhi->GetNbinsX();
		if(i==5) nbins=HistF.fAna_Njet15->GetNbinsX();
		if(i==6) nbins=HistF.fAna_Qt->GetNbinsX();
		if(i==7) nbins=HistF.fAna_Et1->GetNbinsX();
		if(i==8) nbins=HistF.fAna_Et2->GetNbinsX();
		if(i==9) nbins=HistF.fAna_Etjet->GetNbinsX();
		if(i==10) nbins=HistF.fAna_Ht->GetNbinsX();
		if(i==11) nbins=HistF.fAna_Mjj->GetNbinsX();
		if(i==12) nbins=HistF.fAna_Nem->GetNbinsX();
		if(i==13) nbins=HistF.fAna_Mej->GetNbinsX();
		if(i==14) nbins=HistF.fAna_Mextra->GetNbinsX();
		if(i==15) nbins=HistF.fAna_Etem->GetNbinsX();
		if(i==16) nbins=HistF.fAna_dPhi1->GetNbinsX();
		if(i==17) nbins=HistF.fAna_dPhi2->GetNbinsX();
		if(i==18) nbins=HistF.fAna_dPhi3->GetNbinsX();
		if(i==19) nbins=HistF.fAna_dPhi4->GetNbinsX();
		if(i==20) nbins=HistF.fAna_dPhi5->GetNbinsX();
		if(i==21) nbins=HistF.fAna_dPhi6->GetNbinsX();
		if(i==22) nbins=HistF.fAna_dPhi7->GetNbinsX();
		if(i==23) nbins=HistF.fAna_Njet20->GetNbinsX();
		if(i==24) nbins=HistF.fAna_Njet25->GetNbinsX();
		if(i==25) nbins=HistF.fAna_Njet30->GetNbinsX();
		if(i==26) nbins=HistF.fAna_Njet35->GetNbinsX();

		for(int j=0; j<nbins; j++)
		{
			if(i==0) HistF.fAna_MetAll->SetBinError(j+1,HistoBinDiff(HistF.fAna_MetAll,Hist1.fAna_MetAll,Hist2.fAna_MetAll,j));
			if(i==1) HistF.fAna_MetSig->SetBinError(j+1,HistoBinDiff(HistF.fAna_MetSig,Hist1.fAna_MetSig,Hist2.fAna_MetSig,j));
			if(i==2) HistF.fAna_Met->SetBinError(j+1,HistoBinDiff(HistF.fAna_Met,Hist1.fAna_Met,Hist2.fAna_Met,j));
			if(i==3) HistF.fAna_M->SetBinError(j+1,HistoBinDiff(HistF.fAna_M,Hist1.fAna_M,Hist2.fAna_M,j));
			if(i==4) HistF.fAna_dPhi->SetBinError(j+1,HistoBinDiff(HistF.fAna_dPhi,Hist1.fAna_dPhi,Hist2.fAna_dPhi,j));
			if(i==5) HistF.fAna_Njet15->SetBinError(j+1,HistoBinDiff(HistF.fAna_Njet15,Hist1.fAna_Njet15,Hist2.fAna_Njet15,j));
			if(i==6) HistF.fAna_Qt->SetBinError(j+1,HistoBinDiff(HistF.fAna_Qt,Hist1.fAna_Qt,Hist2.fAna_Qt,j));
			if(i==7) HistF.fAna_Et1->SetBinError(j+1,HistoBinDiff(HistF.fAna_Et1,Hist1.fAna_Et1,Hist2.fAna_Et1,j));
			if(i==8) HistF.fAna_Et2->SetBinError(j+1,HistoBinDiff(HistF.fAna_Et2,Hist1.fAna_Et2,Hist2.fAna_Et2,j));
			if(i==9) HistF.fAna_Etjet->SetBinError(j+1,HistoBinDiff(HistF.fAna_Etjet,Hist1.fAna_Etjet,Hist2.fAna_Etjet,j));
			if(i==10) HistF.fAna_Ht->SetBinError(j+1,HistoBinDiff(HistF.fAna_Ht,Hist1.fAna_Ht,Hist2.fAna_Ht,j));
			if(i==11) HistF.fAna_Mjj->SetBinError(j+1,HistoBinDiff(HistF.fAna_Mjj,Hist1.fAna_Mjj,Hist2.fAna_Mjj,j));
			if(i==12) HistF.fAna_Nem->SetBinError(j+1,HistoBinDiff(HistF.fAna_Nem,Hist1.fAna_Nem,Hist2.fAna_Nem,j));
			if(i==13) HistF.fAna_Mej->SetBinError(j+1,HistoBinDiff(HistF.fAna_Mej,Hist1.fAna_Mej,Hist2.fAna_Mej,j));
			if(i==14) HistF.fAna_Mextra->SetBinError(j+1,HistoBinDiff(HistF.fAna_Mextra,Hist1.fAna_Mextra,Hist2.fAna_Mextra,j));
			if(i==15) HistF.fAna_Etem->SetBinError(j+1,HistoBinDiff(HistF.fAna_Etem,Hist1.fAna_Etem,Hist2.fAna_Etem,j));
			if(i==16) HistF.fAna_dPhi1->SetBinError(j+1,HistoBinDiff(HistF.fAna_dPhi1,Hist1.fAna_dPhi1,Hist2.fAna_dPhi1,j));
			if(i==17) HistF.fAna_dPhi2->SetBinError(j+1,HistoBinDiff(HistF.fAna_dPhi2,Hist1.fAna_dPhi2,Hist2.fAna_dPhi2,j));
			if(i==18) HistF.fAna_dPhi3->SetBinError(j+1,HistoBinDiff(HistF.fAna_dPhi3,Hist1.fAna_dPhi3,Hist2.fAna_dPhi3,j));
			if(i==19) HistF.fAna_dPhi4->SetBinError(j+1,HistoBinDiff(HistF.fAna_dPhi4,Hist1.fAna_dPhi4,Hist2.fAna_dPhi4,j));
			if(i==20) HistF.fAna_dPhi5->SetBinError(j+1,HistoBinDiff(HistF.fAna_dPhi5,Hist1.fAna_dPhi5,Hist2.fAna_dPhi5,j));
			if(i==21) HistF.fAna_dPhi6->SetBinError(j+1,HistoBinDiff(HistF.fAna_dPhi6,Hist1.fAna_dPhi6,Hist2.fAna_dPhi6,j));
			if(i==22) HistF.fAna_dPhi7->SetBinError(j+1,HistoBinDiff(HistF.fAna_dPhi7,Hist1.fAna_dPhi7,Hist2.fAna_dPhi7,j));
			if(i==23) HistF.fAna_Njet20->SetBinError(j+1,HistoBinDiff(HistF.fAna_Njet20,Hist1.fAna_Njet20,Hist2.fAna_Njet20,j));
			if(i==24) HistF.fAna_Njet25->SetBinError(j+1,HistoBinDiff(HistF.fAna_Njet25,Hist1.fAna_Njet25,Hist2.fAna_Njet25,j));
			if(i==25) HistF.fAna_Njet30->SetBinError(j+1,HistoBinDiff(HistF.fAna_Njet30,Hist1.fAna_Njet30,Hist2.fAna_Njet30,j));
			if(i==26) HistF.fAna_Njet35->SetBinError(j+1,HistoBinDiff(HistF.fAna_Njet35,Hist1.fAna_Njet35,Hist2.fAna_Njet35,j));

		}
	}
	return;
}

//_____ counts data & background(metmodel) events with MetSig>cut
//----- systcode=-1 is for data, systcode=0,..19 for background
void JetFilterModuleV2::MyEventCount(double met,double metsig,int systcode,MetResults &metstuff)
{
	if(metsig>=fMetSig_cut)
	{
		int njet15=0;
		int njet20=0;
		int njet25=0;
		njet15= (jet04stuff.myNjet_th15 < 10) ? jet04stuff.myNjet_th15 : 9;
		njet20= (jet04stuff.myNjet_th20 < 10) ? jet04stuff.myNjet_th20 : 9;
		njet25= (jet04stuff.myNjet_th25 < 10) ? jet04stuff.myNjet_th25 : 9;

		if(met>20)
		{
			if(systcode==-1)
			{
				metstuff.ana20Njet_dt[0][njet15]++;
				metstuff.ana20Njet_dt[1][njet20]++;
				metstuff.ana20Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana20Njet_bg[systcode][0][njet15]++;
				metstuff.ana20Njet_bg[systcode][1][njet20]++;
				metstuff.ana20Njet_bg[systcode][2][njet25]++;
			}
		}
		if(met>25)
		{
			if(systcode==-1)
			{
				metstuff.ana25Njet_dt[0][njet15]++;
				metstuff.ana25Njet_dt[1][njet20]++;
				metstuff.ana25Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana25Njet_bg[systcode][0][njet15]++;
				metstuff.ana25Njet_bg[systcode][1][njet20]++;
				metstuff.ana25Njet_bg[systcode][2][njet25]++;
			}
		}
		if(met>30)
		{
			if(systcode==-1)
			{
				metstuff.ana30Njet_dt[0][njet15]++;
				metstuff.ana30Njet_dt[1][njet20]++;
				metstuff.ana30Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana30Njet_bg[systcode][0][njet15]++;
				metstuff.ana30Njet_bg[systcode][1][njet20]++;
				metstuff.ana30Njet_bg[systcode][2][njet25]++;
			}
		}
		if(met>35)
		{
			if(systcode==-1)
			{
				metstuff.ana35Njet_dt[0][njet15]++;
				metstuff.ana35Njet_dt[1][njet20]++;
				metstuff.ana35Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana35Njet_bg[systcode][0][njet15]++;
				metstuff.ana35Njet_bg[systcode][1][njet20]++;
				metstuff.ana35Njet_bg[systcode][2][njet25]++;
			}
		}
		if(met>40)
		{
			if(systcode==-1)
			{
				metstuff.ana40Njet_dt[0][njet15]++;
				metstuff.ana40Njet_dt[1][njet20]++;
				metstuff.ana40Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana40Njet_bg[systcode][0][njet15]++;
				metstuff.ana40Njet_bg[systcode][1][njet20]++;
				metstuff.ana40Njet_bg[systcode][2][njet25]++;
			}
		}
		if(met>45)
		{
			if(systcode==-1)
			{
				metstuff.ana45Njet_dt[0][njet15]++;
				metstuff.ana45Njet_dt[1][njet20]++;
				metstuff.ana45Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana45Njet_bg[systcode][0][njet15]++;
				metstuff.ana45Njet_bg[systcode][1][njet20]++;
				metstuff.ana45Njet_bg[systcode][2][njet25]++;
			}
		}
		if(met>50)
		{
			if(systcode==-1)
			{
				metstuff.ana50Njet_dt[0][njet15]++;
				metstuff.ana50Njet_dt[1][njet20]++;
				metstuff.ana50Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana50Njet_bg[systcode][0][njet15]++;
				metstuff.ana50Njet_bg[systcode][1][njet20]++;
				metstuff.ana50Njet_bg[systcode][2][njet25]++;
			}
		}
		if(met>75)
		{
			if(systcode==-1)
			{
				metstuff.ana75Njet_dt[0][njet15]++;
				metstuff.ana75Njet_dt[1][njet20]++;
				metstuff.ana75Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana75Njet_bg[systcode][0][njet15]++;
				metstuff.ana75Njet_bg[systcode][1][njet20]++;
				metstuff.ana75Njet_bg[systcode][2][njet25]++;
			}
		}
		if(met>100)
		{
			if(systcode==-1)
			{
				metstuff.ana100Njet_dt[0][njet15]++;
				metstuff.ana100Njet_dt[1][njet20]++;
				metstuff.ana100Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana100Njet_bg[systcode][0][njet15]++;
				metstuff.ana100Njet_bg[systcode][1][njet20]++;
				metstuff.ana100Njet_bg[systcode][2][njet25]++;
			}
		}
		if(met>150)
		{
			if(systcode==-1)
			{
				metstuff.ana150Njet_dt[0][njet15]++;
				metstuff.ana150Njet_dt[1][njet20]++;
				metstuff.ana150Njet_dt[2][njet25]++;
			}
			else
			{
				metstuff.ana150Njet_bg[systcode][0][njet15]++;
				metstuff.ana150Njet_bg[systcode][1][njet20]++;
				metstuff.ana150Njet_bg[systcode][2][njet25]++;
			}
		}
	}
	return;
}
//_____________________________________ generates Met predictions and fills out analysis histograms
void JetFilterModuleV2::DoMyAnalysis(JetStuff jetstuff,CommonStuff miscstuff)
{
	if(fAnalysisMode==1)  //don't need this anymore. I moved this check to the Event loop. 01-29-2009 sam
	{
		//_____________________________ check the njet requirement for DATA

		jetstuff.smear_factor.clear();
		int metcode=3;
		//______________________ data
		double metSig_est=MyTotalMetSignificance(jetstuff,jetstuff.myMETcorr_th15,0,0);
		myMetSig=metSig_est; // filling my Met Significance
		int Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,jetstuff.myMETcorr_th15); // met cleanup cut
		if(Nbadmet==0)
		{
			if (MyNjetCut(jet04stuff)==1) {
				FillAnalysisHistograms(fAna_data,jetstuff.Jet_lev6_noEMobj, jetstuff, miscstuff,jetstuff.myMETcorr_th15,metSig_est);// filling analysis histograms
				MyEventCount(jetstuff.myMETcorr_th15.Mod(),metSig_est,-1,met_results);
			}
		} // if NBadmet == 0



		//______________________ MetModel part: pseudo-experiments
		for(int k=0; k<Npoints; k++)
		{
			//______________________ Met Model default prediction
			TVector2 myGenMet(0.0,0.0);
			int jer_code = 0; // JER systematics code
			int ue_code  = 0; // unclustered energy parameterization systematics
			
			TRandom3 *rn = new TRandom3();		//should not use gRandom->GetSeed() twice as it returns the same seed.
			int rnd_seed_1 = rn->GetSeed();
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	

			int rnd_seed = gRandom->GetSeed();
			//std::cout << "seed , seed1 = " << rnd_seed << "\t" << rnd_seed_1 << std::endl;
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);		//need to stuff in smeared jets once per loop

			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			myMetSig_gen=metSig_est;
			myGenMetVec.Set(myGenMet.Px(),myGenMet.Py());
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut

		
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff)==1)
				{
					FillAnalysisHistograms(fAna_def,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					FillAnalysisHistograms(fAna_bckg,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est); // Met Model total prediction (stat+syst)
					MyEventCount(myGenMet.Mod(),metSig_est,0,met_results);
				}
			}

			//______________________ Met Model systematics: z->ee vs. dipho sideband parameterization
			myGenMet.Set(0.0,0.0);
			jer_code=0;
			ue_code=1;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_ue1,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,1,met_results);
				}
			}
			//______________________ Met Model Uncl. En. systematics: mean-G
			myGenMet.Set(0.0,0.0);
			jer_code=0;
			ue_code=-2;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_ue2,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,2,met_results);
				}
			}
			//______________________ Met Model Uncl. En. systematics: mean+G
			myGenMet.Set(0.0,0.0);
			jer_code=0;
			ue_code=2;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_ue3,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,3,met_results);
				}
			}
			//______________________ Met Model Uncl. En. systematics: sigma-G
			myGenMet.Set(0.0,0.0);
			jer_code=0;
			ue_code=-3;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_ue4,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,4,met_results);
				}
			}
			//______________________ Met Model Uncl. En. systematics: sigma+G
			myGenMet.Set(0.0,0.0);
			jer_code=0;
			ue_code=3;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_ue5,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,5,met_results);
				}
			}
			//______________________ Met Model Uncl. En. systematics: scale-G
			myGenMet.Set(0.0,0.0);
			jer_code=0;
			ue_code=-4;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_ue6,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,6,met_results);
				}
			}
			//______________________ Met Model Uncl. En. systematics: scale+G
			myGenMet.Set(0.0,0.0);
			jer_code=0;
			ue_code=4;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_ue7,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,7,met_results);
				}
			}
			//______________________ Met Model Uncl. En. systematics: norm-G
			myGenMet.Set(0.0,0.0);
			jer_code=0;
			ue_code=-5;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_ue8,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,8,met_results);
				}
			}
			//______________________ Met Model Uncl. En. systematics: norm+G
			myGenMet.Set(0.0,0.0);
			jer_code=0;
			ue_code=5;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_ue9,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,9,met_results);
				}
			}

			//______________________ Met Model JER systematics: meanG-G
			myGenMet.Set(0.0,0.0);
			jer_code=-1;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer1,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,10,met_results);
				}
			}
			//______________________ Met Model JER systematics: meanG+G
			myGenMet.Set(0.0,0.0);
			jer_code=1;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer2,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,11,met_results);
				}
			}
			//______________________ Met Model JER systematics: sigmaG-G
			myGenMet.Set(0.0,0.0);
			jer_code=-2;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer3,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,12,met_results);
				}
			}
			//______________________ Met Model JER systematics: sigmaG+G
			myGenMet.Set(0.0,0.0);
			jer_code=2;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer4,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,13,met_results);
				}
			}
			//______________________ Met Model JER systematics: mpvL-G
			myGenMet.Set(0.0,0.0);
			jer_code=-3;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer5,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,14,met_results);
				}
			}
			//______________________ Met Model JER systematics: mpvL+G
			myGenMet.Set(0.0,0.0);
			jer_code=3;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer6,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,15,met_results);
				}
			}
			//______________________ Met Model JER systematics: sigmaL-G
			myGenMet.Set(0.0,0.0);
			jer_code=-4;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer7,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,16,met_results);
				}
			}
			//______________________ Met Model JER systematics: sigmaL+G
			myGenMet.Set(0.0,0.0);
			jer_code=4;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer8,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,17,met_results);
				}
			}
			//______________________ Met Model JER systematics: norm-G
			myGenMet.Set(0.0,0.0);
			jer_code=-5;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer9,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,18,met_results);
				}
			}
			//______________________ Met Model JER systematics: norm+G
			myGenMet.Set(0.0,0.0);
			jer_code=5;
			ue_code=0;
			Smear_newJetL6FirstTime(jetstuff, jer_code,metcode,rnd_seed_1);	
			GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,jer_code,ue_code,metcode,rnd_seed,myGenMet);
			GenerateSmearedJets(jetstuff);
			metSig_est=MyTotalMetSignificance(jetstuff,myGenMet,jer_code,ue_code);
			Nbadmet=MyMetCleanUpCut(jetstuff,miscstuff,jetstuff.Jet_lev6_noEMobj,myGenMet); // met cleanup cut
			if(Nbadmet==0)
			{
				if (MySmeared_NjetCut(jetstuff))
				{
					FillAnalysisHistograms(fAna_jer10,jetstuff.smeared_newJetLev6_noEMObj,jetstuff,miscstuff,myGenMet,metSig_est);
					MyEventCount(myGenMet.Mod(),metSig_est,19,met_results);
				}
			}
		} // for loop over  k points to generate

		//std::cout << "D"<< std::endl;
	}  //if Anamode == 1

	//std::cout << "\t\tRetrun point\tmyNjet_th15=" << jetstuff.myNjet_th15<< std::endl;
	return;
} //DoMyAnalysis



// this is where the calibration plots are made
//_________________________________________________ filling met histo
void JetFilterModuleV2::FillMetStudyHistograms(MetStudyHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff,
		int metcode, int nvx, double dzvx, double dzvx_worse, double zvx[3]) {

	double mgg=0.0;
	if(miscstuff.myCorrPhoton.size()>1) mgg=(miscstuff.myCorrPhoton[0]+miscstuff.myCorrPhoton[1]).M();
	double ht=MyHtAll(jetstuff.Jet_lev6_noEMobj ,jetstuff,miscstuff,metcode);

	double sqrtSumEt;
	double SumEt;
	double sqrtSumEt_jet;
	double SumEt_jet;
	double JetFr=0.0;
	double met;
	double metR=miscstuff.myMET_raw.Mod();
	double metR0=miscstuff.myMET0_raw.Mod();
	double metx;
	double mety;
	double metRx=miscstuff.myMET_raw.Px();
	double metRy=miscstuff.myMET_raw.Py();
	double metphiR=miscstuff.myMET_raw.Phi();
	double metphiC;
	int njet=0;

	if(metcode==1)
	{
		SumEt=jetstuff.mySumEtCorr_th5;
		sqrtSumEt=sqrt(jetstuff.mySumEtCorr_th5);
		SumEt_jet=jetstuff.mySumEtJet_th5;
		sqrtSumEt_jet=sqrt(jetstuff.mySumEtJet_th5);
		met=jetstuff.myMETcorr_th5.Mod();
		metx=jetstuff.myMETcorr_th5.Px();
		mety=jetstuff.myMETcorr_th5.Py();
		metphiC=jetstuff.myMETcorr_th5.Phi();
		njet=jetstuff.myNjet_th5;
	}
	if(metcode==2)
	{
		SumEt=jetstuff.mySumEtCorr_th10;
		sqrtSumEt=sqrt(jetstuff.mySumEtCorr_th10);
		SumEt_jet=jetstuff.mySumEtJet_th10;
		sqrtSumEt_jet=sqrt(jetstuff.mySumEtJet_th10);
		met=jetstuff.myMETcorr_th10.Mod();
		metx=jetstuff.myMETcorr_th10.Px();
		mety=jetstuff.myMETcorr_th10.Py();
		metphiC=jetstuff.myMETcorr_th10.Phi();
		njet=jetstuff.myNjet_th10;
	}
	if(metcode==3)
	{
		SumEt=jetstuff.mySumEtCorr_th15;
		sqrtSumEt=sqrt(jetstuff.mySumEtCorr_th15);
		SumEt_jet=jetstuff.mySumEtJet_th15;
		sqrtSumEt_jet=sqrt(jetstuff.mySumEtJet_th15);
		met=jetstuff.myMETcorr_th15.Mod();
		metx=jetstuff.myMETcorr_th15.Px();
		mety=jetstuff.myMETcorr_th15.Py();
		metphiC=jetstuff.myMETcorr_th15.Phi();
		njet=jetstuff.myNjet_th15;
	}
	if(metcode==4)
	{
		SumEt=jetstuff.mySumEtCorr_th20;
		sqrtSumEt=sqrt(jetstuff.mySumEtCorr_th20);
		SumEt_jet=jetstuff.mySumEtJet_th20;
		sqrtSumEt_jet=sqrt(jetstuff.mySumEtJet_th20);
		met=jetstuff.myMETcorr_th20.Mod();
		metx=jetstuff.myMETcorr_th20.Px();
		mety=jetstuff.myMETcorr_th20.Py();
		metphiC=jetstuff.myMETcorr_th20.Phi();
		njet=jetstuff.myNjet_th20;
	}

	if(metcode<1 || metcode>4)
	{
		SumEt=miscstuff.mySumEt_raw;
		sqrtSumEt=sqrt(miscstuff.mySumEt_raw);
		SumEt_jet=0.0; // not defined for this case
		sqrtSumEt_jet=0.0; // not defined for this case
		met=miscstuff.myMET_raw.Mod();
		metx=miscstuff.myMET_raw.Px();
		mety=miscstuff.myMET_raw.Py();
		metphiC=miscstuff.myMET_raw.Phi();
		njet=jetstuff.myNjet;
	}
	if(SumEt>0.0) JetFr=SumEt_jet/SumEt;

	if(fabs(met-metR)>0.0)
	{
		Hist.fMetRaw->Fill(metR);
		Hist.fMetCorr->Fill(met);
		Hist.fMetPhiRaw->Fill(metphiR);
		Hist.fMetPhiCorr->Fill(metphiC);
	}

	Hist.fSumEtCorr->Fill(SumEt);
	Hist.fSqrtSumEtCorr->Fill(sqrtSumEt);

	if(njet>0)
	{
		Hist.fSumEtCorr_withJet->Fill(SumEt);
		Hist.fSqrtSumEtCorr_withJet->Fill(sqrtSumEt);
	}
	else
	{
		Hist.fSumEtCorr_noJet->Fill(SumEt);
		Hist.fSqrtSumEtCorr_noJet->Fill(sqrtSumEt);
		if(nvx==1) Hist.fSumEtCorr_noJet_vx1->Fill(SumEt);
		if(nvx==2) Hist.fSumEtCorr_noJet_vx2->Fill(SumEt);
		if(nvx==3) Hist.fSumEtCorr_noJet_vx3->Fill(SumEt);
		Hist.fSumEtCorrNoJet_vs_Nvx->Fill(nvx,SumEt);
	}

	Hist.fMetRawAll->Fill(metR);
	Hist.fMetCorrAll->Fill(met);
	Hist.fMetRawAll_X->Fill(metRx);
	Hist.fMetCorrAll_X->Fill(metx);
	Hist.fMetRawAll_Y->Fill(metRy);
	Hist.fMetCorrAll_Y->Fill(mety);
	Hist.fMetPhiRawAll->Fill(metphiR);
	Hist.fMetPhiCorrAll->Fill(metphiC);

	//   if(nvx>1) Hist.fMetCorr_vs_dZvx->Fill(fabs(dzvx),met);
	if(dzvx_worse>0.0)
	{
		Hist.fMetCorr_vs_dZvx->Fill(fabs(dzvx),met);
		Hist.fMetCorr_vs_dZvxWorse->Fill(fabs(dzvx_worse),met);
		Hist.fMetCorr_vs_Zvx2ndBest->Fill(zvx[1],met);
		Hist.fMetCorr_vs_ZvxWorse->Fill(zvx[2],met);
		Hist.fMetCorrVsZvx2ndBest->Fill(zvx[1],met);
		Hist.fMetCorrVsZvxWorse->Fill(zvx[2],met);
	}
	Hist.fMetCorr_vs_ZvxBest->Fill(fabs(zvx[0]),met);
	Hist.fMetCorrVsZvxBest->Fill(fabs(zvx[0]),met);
	if(nvx=1)
	{
		Hist.fMet0Met4_vs_Zvx->Fill(zvx[0],metR0-metR);
		Hist.fMet0_vs_Zvx->Fill(zvx[0],metR0);
	}

	//---- Below this point, all histograms for met model will be filled using the last generated value of met

	Hist.fMet4VsNjet->Fill(jetstuff.myNjet,met);
	Hist.fMet4VsNjet_th5->Fill(jetstuff.myNjet_th5,met);
	Hist.fMet4VsNjet_th10->Fill(jetstuff.myNjet_th10,met);
	Hist.fMet4VsNjet_th15->Fill(jetstuff.myNjet_th15,met);
	Hist.fMet4VsNjet_th20->Fill(jetstuff.myNjet_th20,met);
	Hist.fMet4VsNjet_th25->Fill(jetstuff.myNjet_th25,met);
	Hist.fMet4VsNjet_th30->Fill(jetstuff.myNjet_th30,met);
	Hist.fMet4VsNjet_th35->Fill(jetstuff.myNjet_th35,met);

	Hist.fMet4VsNvx12->Fill(nvx,met);
	Hist.fMet4X_vs_Met4Y->Fill(metx,mety);

	if(nvx==1) Hist.fMet4VsRun->Fill(GetHeaderBlock()->RunNumber(),met);
	Hist.fNvx12VsRun->Fill(GetHeaderBlock()->RunNumber(),nvx);

	if(SumEt_jet>0.0)
	{
		Hist.fSumEtJetFrac->Fill(JetFr);
		Hist.fSumEtJet->Fill(SumEt_jet);
		Hist.fSqrtSumEtJet->Fill(sqrtSumEt_jet);
	}
	Hist.fMet4VsJetFr->Fill(JetFr,met);
	Hist.fMet4VsSqrtSumEtJet->Fill(sqrtSumEt_jet,met);
	Hist.fJetFrVsSumEt->Fill(SumEt,JetFr);

	for(int i=0; i<10; i++)
	{
		if(sqrtSumEt>i*2.0 && sqrtSumEt<=(i+1)*2.0)
		{

			Hist.fMet4[i]->Fill(met);
			Hist.fMet4Phi[i]->Fill(metphiC);

			if(njet==0)
			{
				Hist.fMet4X_noJet[i]->Fill(metx);
				Hist.fMet4Y_noJet[i]->Fill(mety);
			}
			else
			{
				Hist.fMet4X_withJet[i]->Fill(metx);
				Hist.fMet4Y_withJet[i]->Fill(mety);
			}

			Hist.fMet4X[i]->Fill(metx);
			Hist.fMet4Y[i]->Fill(mety);
			break; // just in case if nvx>5
		}
	}

	//__________________________ filling histograms for Met Model-2
	TVector2 myGenV2Met_def(0.0,0.0);
	TVector2 myGenV2Met_def_max(0.0,0.0);

	double metSig_est=MyTotalMetSignificance(jetstuff,jetstuff.myMETcorr_th15,0,0);
	myMetSig=metSig_est; // filling my Met Significance

	MetSigDelPhiStuff msdpStuff_dt;
	CalculateMetDelPhi(msdpStuff_dt,jetstuff,miscstuff,jetstuff.myMETcorr_th15);

	//______________ calculating dPhi(met-obj)
	int object_type=-1; // 0=e/pho, 1=jet, 2=b-tag, 3=mu, 4=tau, 5=trk
	if(msdpStuff_dt.dPhi_det[2]<msdpStuff_dt.dPhi_det[1]) object_type=1;
	else object_type=0;

	for(int ij=0; ij<5; ij++)
	{
		if(msdpStuff_dt.dPhi_det[ij]<7.0) Hist.fMetSig_vs_dPhi[ij]->Fill(metSig_est,msdpStuff_dt.dPhi_det[ij]);
	}
	if(ht>=0.0) Hist.fMetSig_vs_sqrtHt->Fill(metSig_est,sqrt(ht));
	if(met>=0.0) Hist.fMetSig_vs_sqrtMet->Fill(metSig_est,sqrt(met));
	if(SumEt>=0.0) Hist.fMetSig_vs_sqrtSumEt->Fill(metSig_est,sqrt(SumEt));

	if(msdpStuff_dt.dPhi_det[1]<0.4 && metSig_est>4.0)
	{
		Hist.fMetSig_vs_Xces->Fill(metSig_est,fabs(msdpStuff_dt.xces_em));
		Hist.fMetSig_vs_etaEM->Fill(metSig_est,fabs(msdpStuff_dt.eta_det_em));
		Hist.fMetSig_Xces_vs_etaEM->Fill(fabs(msdpStuff_dt.xces_em),MyInCemTowerEta(msdpStuff_dt.eta_det_em));
	}
	if(msdpStuff_dt.dPhi_det[2]<0.4 && metSig_est>4.0) Hist.fMetSig_vs_etaJET->Fill(metSig_est,fabs(msdpStuff_dt.eta_det_jet));
	if(object_type>-1) Hist.fObjType->Fill(object_type);

	Hist.fMetSig_estimate->Fill(metSig_est);
	if(njet==0) Hist.fMetSig_estimate_njet0->Fill(metSig_est);
	if(njet==1) Hist.fMetSig_estimate_njet1->Fill(metSig_est);
	if(njet==2) Hist.fMetSig_estimate_njet2->Fill(metSig_est);// for sasha: he has njet>1
	if(njet==3) Hist.fMetSig_estimate_njet3->Fill(metSig_est);//for me need njet>3  - 10-16-2008
	if(njet>=4) Hist.fMetSig_estimate_njet4->Fill(metSig_est);
	int sig_bin=0;
	for(int i_sig=0; i_sig<5; i_sig++)
	{
		if(metSig_est>=1.0*i_sig)
		{
			sig_bin=i_sig;
			Hist.fMetSig_Met[sig_bin]->Fill(met);
		}
	}
	//__ 09/18/07: end of test

	for(int k=0; k<Npoints; k++)
	{
		//------------------------ default parametrization
		int rnd_seed=gRandom->GetSeed();
		GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,0,0,metcode,rnd_seed,myGenV2Met_def);
		GenerateSmearedJets(jetstuff);
		//_________________________________________________________________________________________
		//------------- MetSignificance for pseudo-experiments ------------------------------------
		metSig_est=MyTotalMetSignificance(jetstuff,myGenV2Met_def,0,0);
		myMetSig_gen=metSig_est;
		myGenMetVec.Set(myGenV2Met_def.Px(),myGenV2Met_def.Py());
		MetSigDelPhiStuff msdpStuff_gn;
		CalculateMetDelPhi(msdpStuff_gn,jetstuff,miscstuff,myGenV2Met_def); // calculating dPhi(met-obj)
		for(int ij=0; ij<5; ij++)
		{
			if(msdpStuff_gn.dPhi_det[ij]<7.0) Hist.fMetSigCalib_vs_dPhi[ij]->Fill(metSig_est,msdpStuff_gn.dPhi_det[ij]);
		}

		Hist.fMetSigCalib_estimate->Fill(metSig_est);
		int smearedNjets=0;
		if (metcode==1) smearedNjets = jetstuff.smeared_Njet_th5;
		if (metcode==2) smearedNjets = jetstuff.smeared_Njet_th10;
		if (metcode==3) smearedNjets = jetstuff.smeared_Njet_th15;
		if (metcode==4) smearedNjets = jetstuff.smeared_Njet_th20;
		if (metcode<1 || metcode>4) smearedNjets = jetstuff.smeared_Njet_th0;
		
		if(smearedNjets==0) Hist.fMetSigCalib_estimate_njet0->Fill(metSig_est);
		if(smearedNjets==1) Hist.fMetSigCalib_estimate_njet1->Fill(metSig_est);
		if(smearedNjets==2) Hist.fMetSigCalib_estimate_njet2->Fill(metSig_est);
		if(smearedNjets==3) Hist.fMetSigCalib_estimate_njet3->Fill(metSig_est);
		if(smearedNjets>=4) Hist.fMetSigCalib_estimate_njet4->Fill(metSig_est);
		sig_bin=0;
		for(int i_sig=0; i_sig<5; i_sig++)
		{
			if(metSig_est>=1.0*i_sig)
			{
				sig_bin=i_sig;
				Hist.fMetSigCalib_Met[sig_bin]->Fill(myGenV2Met_def.Mod());
			}
		}
		//____________ End of MetSignificance for pseudo-experiments _______________________________

		//_________________________________________________________________________________________
		//------------- Study of MetSig efficiency: pseudo-experiments and toy MET ----------------
		//------------------------ default parametrization
		double toyMet=MetToy_min+(MetToy_max-MetToy_min)*(gRandom->Rndm());
		double toyMetPhi=TMath::TwoPi()*(gRandom->Rndm());
		double toyMet_x=toyMet*cos(toyMetPhi);
		double toyMet_y=toyMet*sin(toyMetPhi);
		TVector2 MetToy(toyMet_x,toyMet_y);
		TVector2 myGenV2Met_gen(0.0,0.0);
		GenerateMyTotalMet(jetstuff,miscstuff,fGenCleanup,0,0,metcode,rnd_seed,myGenV2Met_gen);
		TVector2 myMetSum(myGenV2Met_gen.Px()+MetToy.Px(),myGenV2Met_gen.Py()+MetToy.Py());
		metSig_est=MyTotalMetSignificance(jetstuff,myMetSum,0,0);
		Hist.fMetSigToy_estimate->Fill(metSig_est);
		sig_bin=0;
		Hist.fMetSigToy_Met[0]->Fill(MetToy.Mod());
		for(int i_sig=1; i_sig<5; i_sig++)
		{
			if(metSig_est>=(1.0*i_sig+2.0))
			{
				sig_bin=i_sig;
				Hist.fMetSigToy_Met[sig_bin]->Fill(MetToy.Mod());
			}
		}
		//____________ End of Study of MetSig efficiency _______________________________

		if(myGenV2Met_def.Mod()>myGenV2Met_def_max.Mod()) myGenV2Met_def_max.Set(myGenV2Met_def.X(),myGenV2Met_def.Y());

		double _met_dphi=TVector2::Phi_0_2pi(myGenV2Met_def.Phi()-metphiC);
		Hist.fGenV2Met_def_proj->Fill(myGenV2Met_def.Mod()*cos(_met_dphi));
		Hist.fGenV2MetXY_def_proj->Fill(myGenV2Met_def.Mod()*cos(_met_dphi),myGenV2Met_def.Mod()*sin(_met_dphi));
		Hist.fGenV2_dPhiMet_def->Fill(fabs(TVector2::Phi_mpi_pi(myGenV2Met_def.Phi()-metphiC)));
		if(jetstuff.myNjet_th15>0)
		{
			Hist.fGenV2Met_def_projNjet15->Fill(myGenV2Met_def.Mod()*cos(_met_dphi));
			Hist.fGenV2MetXY_def_projNjet15->Fill(myGenV2Met_def.Mod()*cos(_met_dphi),myGenV2Met_def.Mod()*sin(_met_dphi));
			Hist.fGenV2_dPhiMet_def_Njet15->Fill(fabs(TVector2::Phi_mpi_pi(myGenV2Met_def.Phi()-metphiC)));
		}
		if(myGenV2Met_def.Mod()<2000.0) // new as of 04/04/07
		{
			Hist.fMet4GenV2Met->Fill(met-myGenV2Met_def.Mod()); // Met4-MetGen
			Hist.fMet4_vs_Met4GenV2Met->Fill(met-myGenV2Met_def.Mod(),met); // Met4 .vs. Met4-MetGen
			Hist.fMet4_vs_GenV2Met->Fill(myGenV2Met_def.Mod(),met); // Met4 .vs. MetGen

			Hist.fGenV2Met_def->Fill(myGenV2Met_def.Mod());
			Hist.fGenV2MetX_def->Fill(myGenV2Met_def.Px());
			Hist.fGenV2MetY_def->Fill(myGenV2Met_def.Py());
			Hist.fGenV2MetPhi_def->Fill(myGenV2Met_def.Phi());

			//---- profile histograms
			Hist.fMetGenV2VsNjet->Fill(jetstuff.myNjet,myGenV2Met_def.Mod());
			Hist.fMetGenV2VsNjet_th5->Fill(jetstuff.myNjet_th5,myGenV2Met_def.Mod());
			Hist.fMetGenV2VsNjet_th10->Fill(jetstuff.myNjet_th10,myGenV2Met_def.Mod());
			Hist.fMetGenV2VsNjet_th15->Fill(jetstuff.myNjet_th15,myGenV2Met_def.Mod());
			Hist.fMetGenV2VsNjet_th20->Fill(jetstuff.myNjet_th20,myGenV2Met_def.Mod());
			Hist.fMetGenV2VsNjet_th25->Fill(jetstuff.myNjet_th25,myGenV2Met_def.Mod());
			Hist.fMetGenV2VsNjet_th30->Fill(jetstuff.myNjet_th30,myGenV2Met_def.Mod());
			Hist.fMetGenV2VsNjet_th35->Fill(jetstuff.myNjet_th35,myGenV2Met_def.Mod());
			Hist.fMetGenV2VsNvx12->Fill(nvx,myGenV2Met_def.Mod());
			Hist.fGenV2MetX_vs_MetY->Fill(myGenV2Met_def.Px(),myGenV2Met_def.Py());
			if(nvx==1) Hist.fMetGenV2VsRun->Fill(GetHeaderBlock()->RunNumber(),myGenV2Met_def.Mod());
			Hist.fMetGenV2VsJetFr->Fill(JetFr,myGenV2Met_def.Mod());
			Hist.fMetGenV2VsSqrtSumEtJet->Fill(sqrtSumEt_jet,myGenV2Met_def.Mod());
		}
		//______ Warning!!! this has to be changed in future to be independent of "if(metcode==3)" condition
		if(metcode==3) MyMetModelEventCount(jetstuff,myGenV2Met_def,met_results); // count generated events with Met>cut
	}

	double _metMax_dphi=TVector2::Phi_0_2pi(myGenV2Met_def_max.Phi()-metphiC);
	Hist.fGenV2Met_max_def_proj->Fill(myGenV2Met_def_max.Mod()*cos(_metMax_dphi));
	Hist.fGenV2MetXY_max_def_proj->Fill(myGenV2Met_def_max.Mod()*cos(_metMax_dphi),myGenV2Met_def_max.Mod()*sin(_metMax_dphi));
	Hist.fGenV2_dPhiMet_max_def->Fill(fabs(TVector2::Phi_mpi_pi(myGenV2Met_def_max.Phi()-metphiC)));
	if(jetstuff.myNjet_th15>0)
	{
		Hist.fGenV2Met_max_def_projNjet15->Fill(myGenV2Met_def_max.Mod()*cos(_metMax_dphi));
		Hist.fGenV2MetXY_max_def_projNjet15->Fill(myGenV2Met_def_max.Mod()*cos(_metMax_dphi),myGenV2Met_def_max.Mod()*sin(_metMax_dphi));
		Hist.fGenV2_dPhiMet_max_def_Njet15->Fill(fabs(TVector2::Phi_mpi_pi(myGenV2Met_def_max.Phi()-metphiC)));
	}
	if(met>10.0)
	{
		Hist.fGenV2Met_max_def_projMet10->Fill(myGenV2Met_def_max.Mod()*cos(_metMax_dphi));
		Hist.fGenV2MetXY_max_def_projMet10->Fill(myGenV2Met_def_max.Mod()*cos(_metMax_dphi),myGenV2Met_def_max.Mod()*sin(_metMax_dphi));
		Hist.fGenV2_dPhiMet_max_def_Met10->Fill(fabs(TVector2::Phi_mpi_pi(myGenV2Met_def_max.Phi()-metphiC)));
		if(jetstuff.myNjet_th15>0)
		{
			Hist.fGenV2Met_max_def_projMet10Njet15->Fill(myGenV2Met_def_max.Mod()*cos(_metMax_dphi));
			Hist.fGenV2MetXY_max_def_projMet10Njet15->Fill(myGenV2Met_def_max.Mod()*cos(_metMax_dphi),myGenV2Met_def_max.Mod()*sin(_metMax_dphi));
			Hist.fGenV2_dPhiMet_max_def_Met10Njet15->Fill(fabs(TVector2::Phi_mpi_pi(myGenV2Met_def_max.Phi()-metphiC)));
		}
	}

	return;
}

//_________________________________________ filling general histo for particular Jet Cone
void JetFilterModuleV2::FillJetHistogramsB(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12) {


	//   double toyMet=MetToy_min+(MetToy_max-MetToy_min)*(gRandom->Rndm());
	//   double toyMetPhi=TMath::TwoPi()*(gRandom->Rndm());
	//   double toyMet_x=toyMet*cos(toyMetPhi);
	//   double toyMet_y=toyMet*sin(toyMetPhi);
	//   TVector2 MetToy(toyMet_x,toyMet_y);
	//   TVector2 MetToyEvnt(stuff.myMETcorr_th15.Px()+toyMet_x,stuff.myMETcorr_th15.Py()+toyMet_y);
	//   if(MyMetCleanUpCut(stuff,allstuff,stuff.Jet_lev6_noEMobj,MetToyEvnt)==0) Hist.fEvnt_toyMET_cut->Fill(MetToy.Mod());
	//   Hist.fEvnt_toyMET_all->Fill(MetToy.Mod());
	Hist.fEvnt_corMET_all->Fill(stuff.myMETcorr_th15.Mod());

	TLorentzVector jetsum(0.0,0.0,0.0,0.0);
	double ave_et=0.0;
	double ave_et0=0.0;
	Hist.fEvntNjet_b->Fill(stuff.myNjet);
	Hist.fEvntNjet5_b->Fill(stuff.myNjet_th5);
	Hist.fEvntNjet10_b->Fill(stuff.myNjet_th10);
	Hist.fEvntNjet15_b->Fill(stuff.myNjet_th15);
	Hist.fEvntNjet20_b->Fill(stuff.myNjet_th20);
	Hist.fEvntNjet25_b->Fill(stuff.myNjet_th25);
	Hist.fEvntNjet30_b->Fill(stuff.myNjet_th30);
	Hist.fEvntNjet35_b->Fill(stuff.myNjet_th35);
	Hist.fEvntdZ_b->Fill(dz);

	for(unsigned int i=0; i<stuff.Jet_raw_noEMobj.size(); i++)
	{
		if((stuff.Npho_match.at(i)+stuff.Nele_match.at(i))==0)
		{
			jetsum=jetsum+stuff.Jet_lev6_noEMobj.at(i);
			ave_et=ave_et+stuff.Jet_lev6_noEMobj.at(i).Pt();
			ave_et0=ave_et0+stuff.Jet_raw_noEMobj.at(i).Pt();
			Hist.fEvntEt0_Njet_b->Fill(stuff.myNjet,stuff.Jet_raw_noEMobj.at(i).Pt());
			Hist.fEvntEt_Njet_b->Fill(stuff.myNjet,stuff.Jet_lev6_noEMobj.at(i).Pt());
			if(i<2)
			{
				Hist.fEvntEt0_b[i]->Fill(stuff.Jet_raw_noEMobj.at(i).Pt());
				Hist.fEvntEt_b[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Pt());
				Hist.fEvntEtaDet_b[i]->Fill(stuff.EtaDet.at(i));
				Hist.fEvntEta_b[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Eta());
				Hist.fEvntPhi_b[i]->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj.at(i).Phi()));
				Hist.fEvntTheta_b[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Theta());
				Hist.fEvntEmFr_b[i]->Fill(stuff.EmFrRaw.at(i));
				Hist.fEvntNTowers_b[i]->Fill(stuff.JetNtwr.at(i));
				Hist.fEvntNTracks_b[i]->Fill(stuff.JetNtrk.at(i));
			}
			else
			{
				Hist.fEvntEt0X_b->Fill(stuff.Jet_raw_noEMobj.at(i).Pt());
				Hist.fEvntEtX_b->Fill(stuff.Jet_lev6_noEMobj.at(i).Pt());
				Hist.fEvntEtaDetX_b->Fill(stuff.EtaDet.at(i));
				Hist.fEvntEtaX_b->Fill(stuff.Jet_lev6_noEMobj.at(i).Eta());
				Hist.fEvntPhiX_b->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj.at(i).Phi()));
				Hist.fEvntThetaX_b->Fill(stuff.Jet_lev6_noEMobj.at(i).Theta());
				Hist.fEvntEmFrX_b->Fill(stuff.EmFrRaw.at(i));
				Hist.fEvntNTowersX_b->Fill(stuff.JetNtwr.at(i));
				Hist.fEvntNTracksX_b->Fill(stuff.JetNtrk.at(i));
				Hist.fEvntDeltaR1x_b->Fill(stuff.Jet_lev6_noEMobj.at(i).DeltaR(stuff.Jet_lev6[0]));
				Hist.fEvntDeltaR2x_b->Fill(stuff.Jet_lev6_noEMobj.at(i).DeltaR(stuff.Jet_lev6[1]));
			}
		}
	}
	if(stuff.Jet_raw_noEMobj.size()>1
			&& (stuff.Npho_match[0]+stuff.Nele_match[0])==0
			&& (stuff.Npho_match[1]+stuff.Nele_match[1])==0)
	{
		Hist.fEvntThetaStar_b->Fill(GetThetaStar(stuff.Jet_lev6_noEMobj));
		Hist.fEvntDeltaPhi_b->Fill(fabs(TVector2::Phi_mpi_pi(stuff.Jet_lev6_noEMobj[0].DeltaPhi(stuff.Jet_lev6_noEMobj[1]))));
		Hist.fEvntDeltaEta_b->Fill(stuff.Jet_lev6_noEMobj[0].Eta()-stuff.Jet_lev6_noEMobj[1].Eta());
		Hist.fEvntDeltaR_b->Fill(stuff.Jet_lev6_noEMobj[1].DeltaR(stuff.Jet_lev6_noEMobj[0]));
		Hist.fEvntMjj_b->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M());
		Hist.fEvntKt2jet_b->Fill(GetMyKtKick(stuff.Jet_lev6_noEMobj));
		Hist.fEvntKtAll_b->Fill(jetsum.Pt());
		Hist.fEvntKt_Mjj_b->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),GetMyKtKick(stuff.Jet_lev6_noEMobj));
		Hist.fEvntKtAll_Mjj_b->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),jetsum.Pt());
	}
	Hist.fEvntNJet_Nvx_b->Fill(nvx12,stuff.myNjet);
	Hist.fEvntNJet5_Nvx_b->Fill(nvx12,stuff.myNjet_th5);
	Hist.fEvntNJet10_Nvx_b->Fill(nvx12,stuff.myNjet_th10);
	Hist.fEvntNJet15_Nvx_b->Fill(nvx12,stuff.myNjet_th15);
	Hist.fEvntNJet20_Nvx_b->Fill(nvx12,stuff.myNjet_th20);
	Hist.fEvntNJet25_Nvx_b->Fill(nvx12,stuff.myNjet_th25);
	Hist.fEvntNJet30_Nvx_b->Fill(nvx12,stuff.myNjet_th30);
	Hist.fEvntNJet35_Nvx_b->Fill(nvx12,stuff.myNjet_th35);
	if(stuff.myNjet>0)
	{
		Hist.fEvntEt0_Nvx12_b->Fill(nvx12,ave_et0/(stuff.myNjet));
		Hist.fEvntEt_Nvx12_b->Fill(nvx12,ave_et/(stuff.myNjet));
	}

	if(nvx12==1) Hist.fEvntNjet_Lum_b->Fill((GetHeaderBlock()->InstLum())*1.0E-30,stuff.myNjet);
	return;
}
//_________________________________________ filling general histo for particular Jet Cone
void JetFilterModuleV2::FillJetHistogramsA(JetGeneral_t& Hist, JetStuff stuff, double dz, int nvx12) {
	TLorentzVector jetsum(0.0,0.0,0.0,0.0);
	double ave_et=0.0;
	double ave_et0=0.0;
	Hist.fEvntNjet_a->Fill(stuff.myNjet);
	Hist.fEvntNjet5_a->Fill(stuff.myNjet_th5);
	Hist.fEvntNjet10_a->Fill(stuff.myNjet_th10);
	Hist.fEvntNjet15_a->Fill(stuff.myNjet_th15);
	Hist.fEvntNjet20_a->Fill(stuff.myNjet_th20);
	Hist.fEvntNjet25_a->Fill(stuff.myNjet_th25);
	Hist.fEvntNjet30_a->Fill(stuff.myNjet_th30);
	Hist.fEvntNjet35_a->Fill(stuff.myNjet_th35);
	Hist.fEvntdZ_a->Fill(dz);

	for(unsigned int i=0; i<stuff.Jet_raw_noEMobj.size(); i++)
	{
		if((stuff.Npho_match.at(i)+stuff.Nele_match.at(i))==0)
		{
			jetsum=jetsum+stuff.Jet_lev6_noEMobj.at(i);
			ave_et=ave_et+stuff.Jet_lev6_noEMobj.at(i).Pt();
			ave_et0=ave_et0+stuff.Jet_raw_noEMobj.at(i).Pt();
			Hist.fEvntEt0_Njet_a->Fill(stuff.myNjet,stuff.Jet_raw_noEMobj.at(i).Pt());
			Hist.fEvntEt_Njet_a->Fill(stuff.myNjet,stuff.Jet_lev6_noEMobj.at(i).Pt());
			if(i<2)
			{
				Hist.fEvntEt0_a[i]->Fill(stuff.Jet_raw_noEMobj.at(i).Pt());
				Hist.fEvntEt_a[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Pt());
				Hist.fEvntEtaDet_a[i]->Fill(stuff.EtaDet.at(i));
				Hist.fEvntEta_a[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Eta());
				Hist.fEvntPhi_a[i]->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj.at(i).Phi()));
				Hist.fEvntTheta_a[i]->Fill(stuff.Jet_lev6_noEMobj.at(i).Theta());
				Hist.fEvntEmFr_a[i]->Fill(stuff.EmFrRaw.at(i));
				Hist.fEvntNTowers_a[i]->Fill(stuff.JetNtwr.at(i));
				Hist.fEvntNTracks_a[i]->Fill(stuff.JetNtrk.at(i));
			}
			else
			{
				Hist.fEvntEt0X_a->Fill(stuff.Jet_raw_noEMobj.at(i).Pt());
				Hist.fEvntEtX_a->Fill(stuff.Jet_lev6_noEMobj.at(i).Pt());
				Hist.fEvntEtaDetX_a->Fill(stuff.EtaDet.at(i));
				Hist.fEvntEtaX_a->Fill(stuff.Jet_lev6_noEMobj.at(i).Eta());
				Hist.fEvntPhiX_a->Fill(TVector2::Phi_0_2pi(stuff.Jet_lev6_noEMobj.at(i).Phi()));
				Hist.fEvntThetaX_a->Fill(stuff.Jet_lev6_noEMobj.at(i).Theta());
				Hist.fEvntEmFrX_a->Fill(stuff.EmFrRaw.at(i));
				Hist.fEvntNTowersX_a->Fill(stuff.JetNtwr.at(i));
				Hist.fEvntNTracksX_a->Fill(stuff.JetNtrk.at(i));
				Hist.fEvntDeltaR1x_a->Fill(stuff.Jet_lev6_noEMobj.at(i).DeltaR(stuff.Jet_lev6_noEMobj[0]));
				Hist.fEvntDeltaR2x_a->Fill(stuff.Jet_lev6_noEMobj.at(i).DeltaR(stuff.Jet_lev6_noEMobj[1]));
			}
		}
	}
	if(stuff.Jet_raw_noEMobj.size()>1
			&& (stuff.Npho_match[0]+stuff.Nele_match[0])==0
			&& (stuff.Npho_match[1]+stuff.Nele_match[1])==0)
	{
		Hist.fEvntThetaStar_a->Fill(GetThetaStar(stuff.Jet_lev6_noEMobj));
		Hist.fEvntDeltaPhi_a->Fill(fabs(TVector2::Phi_mpi_pi(stuff.Jet_lev6_noEMobj[0].DeltaPhi(stuff.Jet_lev6_noEMobj[1]))));
		Hist.fEvntDeltaEta_a->Fill(stuff.Jet_lev6_noEMobj[0].Eta()-stuff.Jet_lev6_noEMobj[1].Eta());
		Hist.fEvntDeltaR_a->Fill(stuff.Jet_lev6_noEMobj[1].DeltaR(stuff.Jet_lev6_noEMobj[0]));
		Hist.fEvntMjj_a->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M());
		Hist.fEvntKt2jet_a->Fill(GetMyKtKick(stuff.Jet_lev6_noEMobj));
		Hist.fEvntKtAll_a->Fill(jetsum.Pt());
		Hist.fEvntKt_Mjj_a->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),GetMyKtKick(stuff.Jet_lev6_noEMobj));
		Hist.fEvntKtAll_Mjj_a->Fill((stuff.Jet_lev6_noEMobj[0]+stuff.Jet_lev6_noEMobj[1]).M(),jetsum.Pt());
	}
	Hist.fEvntNJet_Nvx_a->Fill(nvx12,stuff.myNjet);
	Hist.fEvntNJet5_Nvx_a->Fill(nvx12,stuff.myNjet_th5);
	Hist.fEvntNJet10_Nvx_a->Fill(nvx12,stuff.myNjet_th10);
	Hist.fEvntNJet15_Nvx_a->Fill(nvx12,stuff.myNjet_th15);
	Hist.fEvntNJet20_Nvx_a->Fill(nvx12,stuff.myNjet_th20);
	Hist.fEvntNJet25_Nvx_a->Fill(nvx12,stuff.myNjet_th25);
	Hist.fEvntNJet30_Nvx_a->Fill(nvx12,stuff.myNjet_th30);
	Hist.fEvntNJet35_Nvx_a->Fill(nvx12,stuff.myNjet_th35);

	if(stuff.myNjet>0)
	{
		Hist.fEvntEt0_Nvx12_a->Fill(nvx12,ave_et0/(stuff.myNjet));
		Hist.fEvntEt_Nvx12_a->Fill(nvx12,ave_et/(stuff.myNjet));
	}
	if(nvx12==1) Hist.fEvntNjet_Lum_a->Fill((GetHeaderBlock()->InstLum())*1.0E-30,stuff.myNjet);
	return;
}

//______________________________________________________ filling match histo
void JetFilterModuleV2::FillMatchingHistograms(MatchStudyHisto_t& Hist, JetStuff jetstuff, CommonStuff miscstuff, MatchStuff match)
{
	//hist of EM objects matched and removed from jet list
	
	Hist.fMatchNtwr->Fill(jetstuff.JetNtwr[match.JetInd_match]);
	
	if (match.EmObjType_match==0)
	{
		Hist.fMatchDelR->Fill(jetstuff.Jet_raw[match.JetInd_match].DeltaR(miscstuff.myRawPhoton[match.EmInd_match]));
		Hist.fMatchDelPhi->Fill(fabs(TVector2::Phi_mpi_pi(jetstuff.Jet_raw[match.JetInd_match].Phi()-
						miscstuff.myRawPhoton[match.EmInd_match].Phi())));
		Hist.fMatchDelEta->Fill(fabs(jetstuff.Jet_raw[match.JetInd_match].Eta()-
					miscstuff.myRawPhoton[match.EmInd_match].Eta()));
		Hist.fMatchDelEtaDet->Fill(fabs(jetstuff.EtaDet[match.JetInd_match]-
					miscstuff.myPhoEtaDet[match.EmInd_match]));

		Hist.fMatchEtJet2EtPho_b->Fill(jetstuff.Jet_raw[match.JetInd_match].Pt()/miscstuff.myRawPhoton[match.EmInd_match].Pt());
		Hist.fMatchEtJet2EtPho_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].Pt()/miscstuff.myRawPhoton[match.EmInd_match].Pt());
		//       Hist.fMatchEtJet2EtPho_b->Fill(jetstuff.Jet_raw[match.JetInd_match].E()/miscstuff.myRawPhoton[match.EmInd_match].E());
		//       Hist.fMatchEtJet2EtPho_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].E()/miscstuff.myRawPhoton[match.EmInd_match].E());
	}
	if(match.EmObjType_match==1)
	{
		Hist.fMatchDelR->Fill(jetstuff.Jet_raw[match.JetInd_match].DeltaR(miscstuff.myRawElectron[match.EmInd_match]));
		Hist.fMatchDelPhi->Fill(fabs(TVector2::Phi_mpi_pi(jetstuff.Jet_raw[match.JetInd_match].Phi()-
						miscstuff.myRawElectron[match.EmInd_match].Phi())));
		Hist.fMatchDelEta->Fill(fabs(jetstuff.Jet_raw[match.JetInd_match].Eta()-
					miscstuff.myRawElectron[match.EmInd_match].Eta()));
		Hist.fMatchDelEtaDet->Fill(fabs(jetstuff.EtaDet[match.JetInd_match]-
					miscstuff.myEleEtaDet[match.EmInd_match]));

		Hist.fMatchEtJet2EtPho_b->Fill(jetstuff.Jet_raw[match.JetInd_match].Pt()/miscstuff.myRawElectron[match.EmInd_match].Pt());
		Hist.fMatchEtJet2EtPho_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].Pt()/miscstuff.myRawElectron[match.EmInd_match].Pt());
		//       Hist.fMatchEtJet2EtPho_b->Fill(jetstuff.Jet_raw[match.JetInd_match].E()/miscstuff.myRawElectron[match.EmInd_match].E());
		//       Hist.fMatchEtJet2EtPho_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].E()/miscstuff.myRawElectron[match.EmInd_match].E());
	}
	Hist.fMatchNmatch->Fill(jetstuff.Npho_match[match.JetInd_match]+jetstuff.Nele_match[match.JetInd_match]);
	Hist.fMatchEt_raw_b->Fill(jetstuff.Jet_raw[match.JetInd_match].Pt());
	Hist.fMatchEt_raw_a->Fill(jetstuff.Jet_raw_noEMobj[match.JetInd_match].Pt());
	Hist.fMatchEt_lev6_b->Fill(jetstuff.Jet_lev6[match.JetInd_match].Pt());
	Hist.fMatchEt_lev6_a->Fill(jetstuff.Jet_lev6_noEMobj[match.JetInd_match].Pt());
	if(jetstuff.Jet_raw_noEMobj[match.JetInd_match].Pt()>0.0)
	{
		Hist.fMatchDelRoldnew->Fill(jetstuff.Jet_raw[match.JetInd_match].DeltaR(jetstuff.Jet_raw_noEMobj[match.JetInd_match]));
		Hist.fMatchDelEtaDetoldnew->Fill(jetstuff.EtaDet[match.JetInd_match]-jetstuff.EtaDetCorr[match.JetInd_match]);
	}
	return;
}


//_________________________________ returns normalized kT_perp. for resolution studies
float JetFilterModuleV2::GetMyKtPerp(std::vector<TLorentzVector> vec)
{
	float kT_perp=0.0;
	if(vec.size()>=2)
	{
		float phi12=vec[0].DeltaPhi(vec[1]);
		kT_perp=2*cos(phi12/2.0);
	}
	return kT_perp;
}
//_________________________________ returns normalized kT_parl. for resolution studies
float JetFilterModuleV2::GetMyKtParl(std::vector<TLorentzVector> vec)
{
	float kT_parl=0.0;
	if(vec.size()>=2)
	{
		float pt_1=vec[0].Pt();
		float pt_2=vec[1].Pt();
		float phi12=vec[0].DeltaPhi(vec[1]);
		kT_parl=2.0*(pt_1-pt_2)*sin(phi12/2.0)/(pt_1+pt_2);
	}
	return kT_parl;
}


//______________________________________ reads raw Met, pho, ele info and fills CommonStuff
void JetFilterModuleV2::DoCommonStuff(CommonStuff &miscstuff) {

	//_______________________________________talk to InitSuperPhotons and get Sam's photon info

	int _myNpho=initSpMod->GetSuperPhoSize();

	for(int i=0; i<_myNpho; i++)
	{

		//SetRemoveSidebandPhotons(1) or RemoveSidebandPhoton()
		// this does not work like this.
		// I need to move this selection to the point when the objects are
		//actually removed. not when they are stored in JetFilterV2.

		if ( (RemoveTightPho() && initSpMod->GetSuperPhoton(i)->IsTightPhoton() )
			 || (RemoveLoosePho() && initSpMod->GetSuperPhoton(i)->IsLoosePhoton()) 
			 //|| (RemoveSidebandPho() && ! initSpMod->GetSuperPhoton(i)->IsTightPhoton()  -->not complete yet
			 //    && initSpMod->GetSuperPhoton(i)->IsLoosePhoton())
			 )
		{
			if (PrintLevel()>10)
			{
				std::cout << __FUNCTION__ << "::phos tight, loose = " 
							<< "\t" << initSpMod->GetSuperPhoton(i)->IsTightPhoton()
							<< "\t"<< initSpMod->GetSuperPhoton(i)->IsLoosePhoton() 
							<< std::endl;
			}

			TLorentzVector _pho_raw = initSpMod->GetSuperPhoton(i)->GetRawVec();
			TLorentzVector _pho_cor = initSpMod->GetSuperPhoton(i)->GetCorVec();
			miscstuff.myPhoInd.push_back(initSpMod->GetSuperPhoton(i)->GetPhotonBlockIndex());
			miscstuff.myRawPhoton.push_back(_pho_raw);
			miscstuff.myCorrPhoton.push_back(_pho_cor);
			double hadem=initSpMod->GetSuperPhoton(i)->GetHadEm();
			miscstuff.myPhoEmFr.push_back(1.0/(hadem+1.0));
			miscstuff.myPhoEtaDet.push_back(initSpMod->GetSuperPhoton(i)->GetDetEta());
			miscstuff.myPhoXces.push_back(initSpMod->GetSuperPhoton(i)->GetXCes());
			miscstuff.myPhoZces.push_back(initSpMod->GetSuperPhoton(i)->GetZCes());
			miscstuff.tightId.push_back(initSpMod->GetSuperPhoton(i)->GetTightPhotonId());
		}	

	}


	int _myNele=initSpMod->GetSuperPhoSize();
	for(int i=0; i<_myNele; i++) {
		if ( (RemoveTightPhoLikeEle() && initSpMod->GetSuperPhoton(i)->IsTightElectron())
				|| (RemoveLoosePhoLikeEle() && initSpMod->GetSuperPhoton(i)->IsLooseElectron())
				// || (RemoveStdTightEle() && initSpMod->GetSuperPhoton(i)->IsStdTightElectron())		//not imp. yet
				|| (RemoveStdLooseEle() && initSpMod->GetSuperPhoton(i)->IsStdLooseElectron()) )
		{
			//std::cout << "index - eles tight, loose - Std.Tight, Std. Loose = " 
			//			<< i << " -\t" 
			//			<< initSpMod->GetSuperPhoton(i)->IsTightElectron()
			//			<< "\t"<< initSpMod->GetSuperPhoton(i)->IsLooseElectron() 
			//			<< "\t"<< initSpMod->GetSuperPhoton(i)->IsStdTightElectron() 
			//			<< "\t"<< initSpMod->GetSuperPhoton(i)->IsStdLooseElectron() 
			//			<< std::endl;

			TLorentzVector _ele_raw = initSpMod->GetSuperPhoton(i)->GetRawVec();
			TLorentzVector _ele_cor = initSpMod->GetSuperPhoton(i)->GetCorVec();
			miscstuff.myEleInd.push_back(initSpMod->GetSuperPhoton(i)->GetElectronBlockIndex());
			miscstuff.myRawElectron.push_back(_ele_raw);
			miscstuff.myCorrElectron.push_back(_ele_cor);
			double hadem=initSpMod->GetSuperPhoton(i)->GetHadEm();
			miscstuff.myEleEmFr.push_back(1.0/(hadem+1.0));
			miscstuff.myEleEtaDet.push_back(initSpMod->GetSuperPhoton(i)->GetDetEta());
			miscstuff.myElePhiDet.push_back(initSpMod->GetSuperPhoton(i)->GetDetPhi());
			miscstuff.myEleXces.push_back(initSpMod->GetSuperPhoton(i)->GetXCes());
			miscstuff.myEleZces.push_back(initSpMod->GetSuperPhoton(i)->GetZCes());

			double factor=1.0; // added on 11/07/05
			double extra_PemEleEt=0.0; // extra PEM ele Et due to difference in energy scale, added on 11/07/04
			if (fabs(initSpMod->GetSuperPhoton(i)->GetDetEta())>1.1) // added on 11/07/05
			{
				if (fabs(initSpMod->GetSuperPhoton(i)->GetDetEta())<1.78) factor *= (initSpMod->GetSuperPhoton(i)->GetDetEta()>0.0? 1.020 : 1.015);
				else factor *= (initSpMod->GetSuperPhoton(i)->GetDetEta()>0.0? 1.010 : 1.007);
				extra_PemEleEt=(factor-1.0)*(_ele_raw.Pt());
			}
			else factor=0.0;
			// need to verify with sasha: as I use photon variables and he used electron variables.
			//double pprEt=extra_PemEleEt+factor*(MyZee->GetmyEleCprPpr(i))*fabs((MyZee->GetmyEleFid(i))*sin(_ele_raw->Theta())); // added on 11/04/05
			float Ele_Det = initSpMod->GetSuperPhoton(i)->GetDetector();
			double pprEt= -1.0;
			if (Ele_Det == 0)
			{
				extra_PemEleEt+factor*(fabs(initSpMod->GetSuperPhoton(i)->GetPhoton()->CprEnergy())*sin(_ele_raw.Theta())); // added on 11/04/05
			} else if (Ele_Det == 1)
			{
				extra_PemEleEt+factor*(fabs(initSpMod->GetSuperPhoton(i)->GetPhoton()->PprE())*sin(_ele_raw.Theta())); // added on 11/04/05

			}
			if(pprEt<0.0) pprEt=0.0; // added on 11/04/05
			miscstuff.myElePprEt.push_back(pprEt); // added on 11/04/05
		}
	}



	//_______________________________________ get my raw MET & SumET info

	miscstuff.mySumEt_raw=fMetBlock->Sumet(2);
	miscstuff.myMET_raw.Set(fMetBlock->MetX(4),fMetBlock->MetY(4));
	miscstuff.myMET0_raw.Set(fMetBlock->MetX(0),fMetBlock->MetY(0));

	//____________________________ these are to be added later ______________________________
	//
	//     std::vector<TLorentzVector> myRawMuon;      // raw muons
	//     std::vector<TLorentzVector> myCorrMuon;     // corrected muons (corrections???)
	//     std::vector<TLorentzVector> myRawTau;       // raw taus
	//     std::vector<TLorentzVector> myCorrTau;      // corrected taus  (corrections???, for consistency)
	//     std::vector<TLorentzVector> myRawBjet;      // raw b-jets
	//     std::vector<TLorentzVector> myCorrBjet;     // corrected b-jets  (b-specific corrections???, for consistency)

	return;
}
//___________________________________ corrects Met & SumEt for jets
void JetFilterModuleV2::DoMyMet(CommonStuff miscstuff, JetStuff &jetstuff) {
	//-----correcting SumET
	jetstuff.mySumEtCorr_th5=miscstuff.mySumEt_raw;  // make it the same for a moment
	jetstuff.mySumEtCorr_th10=miscstuff.mySumEt_raw;
	jetstuff.mySumEtCorr_th15=miscstuff.mySumEt_raw;
	jetstuff.mySumEtCorr_th20=miscstuff.mySumEt_raw;
	jetstuff.mySumEtJet_th5=0.0;  // make it the same for a moment
	jetstuff.mySumEtJet_th10=0.0;
	jetstuff.mySumEtJet_th15=0.0;
	jetstuff.mySumEtJet_th20=0.0;
	//_________ correcting all SumEt for photons
	for(unsigned int i=0; i<miscstuff.myRawPhoton.size(); i++)
	{
		jetstuff.mySumEtCorr_th5=jetstuff.mySumEtCorr_th5-miscstuff.myRawPhoton.at(i).Pt();
		jetstuff.mySumEtCorr_th10=jetstuff.mySumEtCorr_th10-miscstuff.myRawPhoton.at(i).Pt();
		jetstuff.mySumEtCorr_th15=jetstuff.mySumEtCorr_th15-miscstuff.myRawPhoton.at(i).Pt();
		jetstuff.mySumEtCorr_th20=jetstuff.mySumEtCorr_th20-miscstuff.myRawPhoton.at(i).Pt();
	}
	//_________ correcting all SumEt for electrons
	for(unsigned int i=0; i<miscstuff.myRawElectron.size(); i++)
	{
		jetstuff.mySumEtCorr_th5=jetstuff.mySumEtCorr_th5-miscstuff.myRawElectron.at(i).Pt();
		jetstuff.mySumEtCorr_th10=jetstuff.mySumEtCorr_th10-miscstuff.myRawElectron.at(i).Pt();
		jetstuff.mySumEtCorr_th15=jetstuff.mySumEtCorr_th15-miscstuff.myRawElectron.at(i).Pt();
		jetstuff.mySumEtCorr_th20=jetstuff.mySumEtCorr_th20-miscstuff.myRawElectron.at(i).Pt();
	}
	//_________ correcting Met for photons and electrons, new as of 09/20/06
	TLorentzVector MetPho(0.0,0.0,0.0,0.0);
	TLorentzVector MetEle(0.0,0.0,0.0,0.0);
	for(unsigned int i=0; i<miscstuff.myRawPhoton.size(); i++)
	{
		MetPho=MetPho+miscstuff.myCorrPhoton.at(i)-miscstuff.myRawPhoton.at(i);
	}
	for(unsigned int i=0; i<miscstuff.myRawElectron.size(); i++)
	{
		MetEle=MetEle+miscstuff.myCorrElectron.at(i)-miscstuff.myRawElectron.at(i);
	}
	if(MetPho.E()<MetPho.P() || MetPho.E()<0.0) MetPho.SetPxPyPzE(0.0,0.0,0.0,0.0);
	if(MetEle.E()<MetEle.P() || MetEle.E()<0.0) MetEle.SetPxPyPzE(0.0,0.0,0.0,0.0);

	//_________ correcting all SumEt and calculating Met correction for jets
	TLorentzVector MetJet5(0.0,0.0,0.0,0.0);
	TLorentzVector MetJet10(0.0,0.0,0.0,0.0);
	TLorentzVector MetJet15(0.0,0.0,0.0,0.0);
	TLorentzVector MetJet20(0.0,0.0,0.0,0.0);
	for(unsigned int i=0; i<jetstuff.Jet_raw_noEMobj.size(); i++)
	{
		//--------- may need an eta cut here // this is ok as we treat everything (or leftover) as unclustered energy- sam
		if(fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())>=MinJetEta && fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<=MaxJetEta)
		{
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>5.0)
			{
				// corr=U.E.+M.I-raw=lev5-lev6-raw+lev1-lev4,
				jetstuff.mySumEtCorr_th5=jetstuff.mySumEtCorr_th5-jetstuff.Jet_raw_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				jetstuff.mySumEtJet_th5=jetstuff.mySumEtJet_th5+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt();
				MetJet5=MetJet5+jetstuff.Jet_lev5_noEMobj.at(i)
					+jetstuff.Jet_lev1_noEMobj.at(i)
					-jetstuff.Jet_lev4_noEMobj.at(i)
					-jetstuff.Jet_raw_noEMobj.at(i);
			}
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>10.0)
			{
				// corr=lev5-raw+M.I.=lev5-raw+lev1-lev4,
				jetstuff.mySumEtCorr_th10=jetstuff.mySumEtCorr_th10-jetstuff.Jet_raw_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				jetstuff.mySumEtJet_th10=jetstuff.mySumEtJet_th10+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt();
				MetJet10=MetJet10+jetstuff.Jet_lev5_noEMobj.at(i)
					+jetstuff.Jet_lev1_noEMobj.at(i)
					-jetstuff.Jet_lev4_noEMobj.at(i)
					-jetstuff.Jet_raw_noEMobj.at(i);
			}
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>15.0)
			{
				// corr=lev5-raw+M.I.=lev5-raw+lev1-lev4,
				jetstuff.mySumEtCorr_th15=jetstuff.mySumEtCorr_th15-jetstuff.Jet_raw_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				jetstuff.mySumEtJet_th15=jetstuff.mySumEtJet_th15+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt();
				MetJet15=MetJet15+jetstuff.Jet_lev5_noEMobj.at(i)
					+jetstuff.Jet_lev1_noEMobj.at(i)
					-jetstuff.Jet_lev4_noEMobj.at(i)
					-jetstuff.Jet_raw_noEMobj.at(i);
			}
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>20.0)
			{
				// corr=lev5-raw+M.I.=lev5-raw+lev1-lev4,
				jetstuff.mySumEtCorr_th20=jetstuff.mySumEtCorr_th20-jetstuff.Jet_raw_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				jetstuff.mySumEtJet_th20=jetstuff.mySumEtJet_th20+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
					+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
					-jetstuff.Jet_lev4_noEMobj.at(i).Pt();
				MetJet20=MetJet20+jetstuff.Jet_lev5_noEMobj.at(i)
					+jetstuff.Jet_lev1_noEMobj.at(i)
					-jetstuff.Jet_lev4_noEMobj.at(i)
					-jetstuff.Jet_raw_noEMobj.at(i);
			}
		}
	}
	//---- making sure all MetJet makes sense:
	if(MetJet5.E()<MetJet5.P()) MetJet5.SetPxPyPzE(0.0,0.0,0.0,0.0);
	if(MetJet10.E()<MetJet10.P()) MetJet10.SetPxPyPzE(0.0,0.0,0.0,0.0);
	if(MetJet15.E()<MetJet15.P()) MetJet15.SetPxPyPzE(0.0,0.0,0.0,0.0);
	if(MetJet20.E()<MetJet20.P()) MetJet20.SetPxPyPzE(0.0,0.0,0.0,0.0);

	//---- this part was modified on 09/20/06
	MetJet5=MetJet5+MetEle+MetPho; // modified on 09/20/06
	MetJet10=MetJet10+MetEle+MetPho; // modified on 09/20/06
	MetJet15=MetJet15+MetEle+MetPho; // modified on 09/20/06
	MetJet20=MetJet20+MetEle+MetPho; // modified on 09/20/06

	jetstuff.myMETcorr_th5.Set(MetJet5.Px(),MetJet5.Py());
	jetstuff.myMETcorr_th10.Set(MetJet10.Px(),MetJet10.Py());
	jetstuff.myMETcorr_th15.Set(MetJet15.Px(),MetJet15.Py());
	jetstuff.myMETcorr_th20.Set(MetJet20.Px(),MetJet20.Py());
	//---- finally, calculating corrected MET
	jetstuff.myMETcorr_th5=miscstuff.myMET_raw-jetstuff.myMETcorr_th5;
	jetstuff.myMETcorr_th10=miscstuff.myMET_raw-jetstuff.myMETcorr_th10;
	jetstuff.myMETcorr_th15=miscstuff.myMET_raw-jetstuff.myMETcorr_th15;
	jetstuff.myMETcorr_th20=miscstuff.myMET_raw-jetstuff.myMETcorr_th20;

	return;
}
//__________________________________ does my jets; main routine
void JetFilterModuleV2::DoMyJet(TStnJetBlock* fJetBlock, CommonStuff miscstuff, JetStuff &jetstuff, MatchStudyHisto_t& Hist) {
	DoMyJetNoMatch(fJetBlock,jetstuff);
	DoMyJetWithMatch(fJetBlock,miscstuff,jetstuff,Hist);
	ReorderMyJets(jetstuff);
	return;
}



//--the following three functions could have been replaced by one which works with generic data type (to be done)
void JetFilterModuleV2::myExchange_tlv(TLorentzVector& val1, TLorentzVector& val2) { //exchanges val1 and val2
	TLorentzVector dummy=val1;
	val1=val2;
	val2=dummy;
	return;
}
void JetFilterModuleV2::myExchange_dbl(double& val1, double& val2) { //exchanges val1 and val2
	double dummy=val1;
	val1=val2;
	val2=dummy;
	return;
}
void JetFilterModuleV2::myExchange_int(int& val1, int& val2) { //exchanges val1 and val2
	int dummy=val1;
	val1=val2;
	val2=dummy;
	return;
}
//_________________________________________ reorders jets after removing EM objects
void JetFilterModuleV2::ReorderMyJets(JetStuff &jetstuff) {
	int reorder_word=0;
	int reorder_word_last=0;
	for(int j=0; j<jetstuff.myNjet-1; j++)
	{
		reorder_word_last=reorder_word;
		for(int i=1; i<jetstuff.myNjet; i++)
		{
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>jetstuff.Jet_lev6_noEMobj.at(i-1).Pt())
			{
				myExchange_tlv(jetstuff.Jet_raw.at(i-1),jetstuff.Jet_raw.at(i));
				myExchange_tlv(jetstuff.Jet_lev1.at(i-1),jetstuff.Jet_lev1.at(i));
				myExchange_tlv(jetstuff.Jet_lev4.at(i-1),jetstuff.Jet_lev4.at(i));
				myExchange_tlv(jetstuff.Jet_lev5.at(i-1),jetstuff.Jet_lev5.at(i));
				myExchange_tlv(jetstuff.Jet_lev6.at(i-1),jetstuff.Jet_lev6.at(i));
				myExchange_tlv(jetstuff.Jet_lev7.at(i-1),jetstuff.Jet_lev7.at(i));
				myExchange_tlv(jetstuff.Jet_raw_noEMobj.at(i-1),jetstuff.Jet_raw_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev1_noEMobj.at(i-1),jetstuff.Jet_lev1_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev4_noEMobj.at(i-1),jetstuff.Jet_lev4_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev5_noEMobj.at(i-1),jetstuff.Jet_lev5_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev6_noEMobj.at(i-1),jetstuff.Jet_lev6_noEMobj.at(i));
				myExchange_tlv(jetstuff.Jet_lev7_noEMobj.at(i-1),jetstuff.Jet_lev7_noEMobj.at(i));

				myExchange_dbl(jetstuff.EtaDet.at(i-1),jetstuff.EtaDet.at(i));
				myExchange_dbl(jetstuff.EtaDetCorr.at(i-1),jetstuff.EtaDetCorr.at(i));
				myExchange_dbl(jetstuff.EmFrRaw.at(i-1),jetstuff.EmFrRaw.at(i));
				myExchange_dbl(jetstuff.EmFrCorr.at(i-1),jetstuff.EmFrCorr.at(i));

				myExchange_int(jetstuff.JetNtrk.at(i-1),jetstuff.JetNtrk.at(i));
				myExchange_int(jetstuff.JetNtwr.at(i-1),jetstuff.JetNtwr.at(i));
				myExchange_int(jetstuff.Nobj_match.at(i-1),jetstuff.Nobj_match.at(i));
				myExchange_int(jetstuff.Npho_match.at(i-1),jetstuff.Npho_match.at(i));
				myExchange_int(jetstuff.Nele_match.at(i-1),jetstuff.Nele_match.at(i));
				myExchange_int(jetstuff.Nmu_match.at(i-1),jetstuff.Nmu_match.at(i));
				myExchange_int(jetstuff.Ntau_match.at(i-1),jetstuff.Ntau_match.at(i));
				myExchange_int(jetstuff.Nbtag_match.at(i-1),jetstuff.Nbtag_match.at(i));
				myExchange_int(jetstuff.JetBlockInd.at(i-1),jetstuff.JetBlockInd.at(i));

				reorder_word++;
			}
		}
		if(reorder_word_last==reorder_word) break;
	}

	//dump jets if the sorting failed (check only leading two jets for now)
	if (jetstuff.Jet_lev6_noEMobj.size()>1) 
	{

		if (jetstuff.Jet_lev6_noEMobj.at(0).Pt() < jetstuff.Jet_lev6_noEMobj.at(1).Pt())
		{
			std::cout << "===================";
			GetHeaderBlock()->Print();
			for (int i=1; i<jetstuff.myNjet; i++)
			{
				std::cout << "i = " << i << " jet= " << jetstuff.Jet_lev6_noEMobj.at(i).Pt() << std::endl;
			}
			std::cout << "===================" << std::endl;
		}
		assert (jetstuff.Jet_lev6_noEMobj.at(0).Pt() >= jetstuff.Jet_lev6_noEMobj.at(1).Pt() &&
					"JetFilterModule::ReorderMyJets::Sorting of Jets has failed!");
	}
	return;
}

//-----------------------------------------------------------------
// This is to sort the smeared jets (from newJetLev6 jets) in Pt
// for pseudo experiments
// The 'map' will automatically sort the inserted in ascending order
// according to the first element. Then read back the corresponding
// indices and refill the original vector whic should be in decending
// order of Pt.
//-----------------------------------------------------------------
void JetFilterModuleV2::ReorderSmearedJets(std::vector<TLorentzVector>& vJets)
{
	if (vJets.size())
	{
		std::vector<TLorentzVector> vjets = vJets;
		std::set<std::pair<double,unsigned> > PtInd;
		std::set<std::pair<double,unsigned> >::iterator myIt;

		for (unsigned int i=0; i < vjets.size(); ++i)
		{
			std::pair<double, unsigned> jetPtInd(vjets.at(i).Pt(), i);
			PtInd.insert(jetPtInd);
		}

		vJets.clear();
		myIt = PtInd.end();
		while (myIt != PtInd.begin())
		{
			myIt--;
			if (debug) std::cout << __FUNCTION__ << "::" << __LINE__ << "Smeared Jets::" << myIt->first << "\t" << myIt->second << std::endl;
			vJets.push_back(vjets.at(myIt->second));
		}


		//I am changing this assertion a little because of following case whic causes
		//(vJets.at(0).Pt() >= vJets.at(1).Pt() to fail due some limitation!
		//ReorderSmearedJets::5837Smeared Jets::1e-20     1
		//ReorderSmearedJets::5837Smeared Jets::1e-20     0
		//ReorderSmearedJets::5841Smeared Jets 0,1 Pt::1e-20      1e-20

		if (vJets.size()>1 && vJets.at(0).Pt()>0.1)
		{
			assert (vJets.at(0).Pt() >= vJets.at(1).Pt() && "ReorderSmearJets::Sorting is wrong!");
		}
	}

} // ReorderSmearedJets	


//____________________________________ re-calculates jet EmFr after removing EM object
double JetFilterModuleV2::CorrectEmFr(double emfr_old, double emfr_obj, double e_old, double e_obj)
{
	double emfr_new = emfr_old;
	
	if(e_old > e_obj)
	{
		double term1 = emfr_obj * e_obj / e_old;
		double denom = 1.0-e_obj / e_old;
		emfr_new = (emfr_old - term1) / denom;
	}
	return emfr_new;
}
//____________________________________ re-calculates jet EtaDet after removing EM object
//_____________ This function has been updated on 06/27/08.
//_____________ Implemented fix for large unphysical DelEta=etaOld-etaNew
double JetFilterModuleV2::CorrectEtaDet(TLorentzVector *vec_pho, TLorentzVector *vec_old,
										double jetetadet_old, double pho_etadet) 
{
	double E_jet 	= vec_old->E();
	double Pt_jet 	= vec_old->Pt();
	double Phi_jet = vec_old->Phi();
	
	TLorentzVector dummy_jet;
	dummy_jet.SetPtEtaPhiE(Pt_jet,jetetadet_old,Phi_jet,E_jet);
	
	double E_pho 	= vec_pho->E();
	double Pt_pho 	= vec_pho->Pt();
	double Phi_pho = vec_pho->Phi();
	
	TLorentzVector dummy_pho;
	dummy_pho.SetPtEtaPhiE(Pt_pho,pho_etadet,Phi_pho,E_pho);
	dummy_jet = dummy_jet - dummy_pho;
	
	double etadet_new = jetetadet_old;
	
	if (dummy_jet.E() > 0.0 && dummy_jet.Pt() > 0.0 && Pt_jet > Pt_pho) etadet_new = dummy_jet.Eta();
	if (fabs(etadet_new-jetetadet_old) > 0.8) etadet_new = jetetadet_old;

	return etadet_new;
}

//__________________________________ does my jets after removing EM objetcs, to be called after DoMyJetNoMatch
void JetFilterModuleV2::DoMyJetWithMatch(TStnJetBlock* fJetBlock, CommonStuff miscstuff, JetStuff &jetstuff, MatchStudyHisto_t& Hist) {

	int _jetcone=0;

	if(fabs(fJetBlock->ConeSize()-0.4)<0.01) _jetcone=0;
	if(fabs(fJetBlock->ConeSize()-0.7)<0.01) _jetcone=1;
	if(fabs(fJetBlock->ConeSize()-1.0)<0.01) _jetcone=2;
	fJTC_coneSize=_jetcone;
	//_________________________ Initializing Jet corrections
	CorInit lev1;
	CorInit lev4;
	CorInit lev5;
	CorInit lev6;
	CorInit lev7;
	lev1.level=1;
	lev1.nvx=myNvx_class12;
	lev1.cone=fJTC_coneSize;
	lev1.version=fJTC_version;
	lev1.sys=fJTC_systcode;
	lev1.imode=fJTC_imode;
	lev1.Nrun=myJetRun;
	lev4.level=4;
	lev4.nvx=myNvx_class12;
	lev4.cone=fJTC_coneSize;
	lev4.version=fJTC_version;
	lev4.sys=fJTC_systcode;
	lev4.imode=fJTC_imode;
	lev4.Nrun=myJetRun;
	lev5.level=5;
	lev5.nvx=myNvx_class12;
	lev5.cone=fJTC_coneSize;
	lev5.version=fJTC_version;
	lev5.sys=fJTC_systcode;
	lev5.imode=fJTC_imode;
	lev5.Nrun=myJetRun;
	lev6.level=6;
	lev6.nvx=myNvx_class12;
	lev6.cone=fJTC_coneSize;
	lev6.version=fJTC_version;
	lev6.sys=fJTC_systcode;
	lev6.imode=fJTC_imode;
	lev6.Nrun=myJetRun;
	lev7.level=7;
	lev7.nvx=myNvx_class12;
	lev7.cone=fJTC_coneSize;
	lev7.version=fJTC_version;
	lev7.sys=fJTC_systcode;
	lev7.imode=fJTC_imode;
	lev7.Nrun=myJetRun;

	MatchStuff dummymatch, sam_match;
	ClearMatchStuff();

	int _Npho_match=0;
	int _Nele_match=0;
	double epsilon=1.0E-10;

	if (PrintLevel()>10)
	{
		std::cout << __FUNCTION__ << "::" << __LINE__ << ":: njets=" << jetstuff.myNjet << std::endl;
	}

	for (int i=0; i<jetstuff.myNjet; i++)
	{
		float _etadet_tmp=jetstuff.EtaDet.at(i); // new line, 10/20/05
		TLorentzVector *_jet=&jetstuff.Jet_raw.at(i);
		//_____ removing photons
		for (unsigned int j=0; j<miscstuff.myRawPhoton.size() && _Npho_match<(int)miscstuff.myRawPhoton.size(); j++)
		{
			TLorentzVector *_pho=&miscstuff.myRawPhoton[j];
			double _match_dR;
			double _match_dEta;
			double _match_dPhi;
			int match_stat=MyMatchPhoJet(jetstuff.JetBlockInd.at(i),miscstuff.myPhoInd[j],-1,
					fJetBlock,_pho,_jet,_match_dR,_match_dPhi,_match_dEta);

			if (PrintLevel()>10)
			{
				std::cout << __FUNCTION__ << "::" << __LINE__ 
								<< ":: match_stat = " << match_stat << "\tjet i, pho j=" 
								<< i << ", " << j << std::endl;
			}

			if (match_stat==1)
			{
				_Npho_match++;
				jetstuff.Nobj_match.at(i)=jetstuff.Nobj_match.at(i)+1;
				jetstuff.Npho_match.at(i)=jetstuff.Npho_match.at(i)+1;

				//i added this to get matching info correctly, above stores a dummy value for EmInd_match
				// sam -- 03/15/2008
				sam_match.JetInd_match=jetstuff.JetBlockInd.at(i) ;
				sam_match.EmInd_match= miscstuff.myPhoInd[j];
				sam_match.EmObjType_match=0; // photons
				sam_matchstuff.push_back(sam_match);

				if ((jetstuff.Npho_match.at(i)+jetstuff.Nele_match.at(i))==1)
				{ // filling only first match
					dummymatch.JetInd_match=i;
					dummymatch.EmInd_match=j;
					dummymatch.EmObjType_match=0; // photons
					matchstuff.push_back(dummymatch);
				}

				// here, I re-calculate jet EmFr
				jetstuff.EmFrCorr.at(i) = CorrectEmFr(jetstuff.EmFrCorr.at(i),miscstuff.myPhoEmFr[j],
						jetstuff.Jet_raw_noEMobj.at(i).E(),miscstuff.myRawPhoton[j].E());
				
				// corrected EtaDet is to be calculated before removing EM object
				_etadet_tmp = CorrectEtaDet(_pho,&jetstuff.Jet_raw_noEMobj.at(i),
						jetstuff.EtaDetCorr.at(i),miscstuff.myPhoEtaDet[j]); // new line, 10/20/05

				jetstuff.EtaDetCorr.at(i) = _etadet_tmp;  // new line, 10/20/05
				double _jet_phi = jetstuff.Jet_raw.at(i).Phi();
				double _jet_eta = jetstuff.Jet_raw.at(i).Eta();
				if (jetstuff.Jet_raw_noEMobj.at(i).E()>miscstuff.myRawPhoton[j].E())
					jetstuff.Jet_raw_noEMobj.at(i)=jetstuff.Jet_raw_noEMobj.at(i)-miscstuff.myRawPhoton[j]; // removing raw photon from raw jet
				else jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);

				// updated from V3
				/*	else jetstuff.Jet_raw_noEMobj.at(i).SetPxPyPzE(0.0,0.0,0.0,0.0);

					if (jetstuff.Jet_raw_noEMobj.at(i).E()<jetstuff.Jet_raw_noEMobj.at(i).P()
					|| jetstuff.Jet_raw_noEMobj.at(i).E()<=0.0
					|| (jetstuff.Jet_raw.at(i).DeltaR(jetstuff.Jet_raw_noEMobj.at(i)))>2.0*(fJetBlock->ConeSize())) {

					jetstuff.Jet_raw_noEMobj.at(i).SetPxPyPzE(0.0,0.0,0.0,0.0); // in case removing EM object leads to M()<0.0
					}
				 */
				if (jetstuff.Jet_raw_noEMobj.at(i).Pt()<epsilon) jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon); // new, 05/28/08

				if (jetstuff.Jet_raw_noEMobj.at(i).E()<jetstuff.Jet_raw_noEMobj.at(i).P()
						|| jetstuff.Jet_raw_noEMobj.at(i).E()<epsilon)
					jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);

				if (jetstuff.Jet_raw.at(i).DeltaR(jetstuff.Jet_raw_noEMobj.at(i))>2.0*(fJetBlock->ConeSize()))
				{
					jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);
				}

			} // matchstat

		}	//photon loop

		//_____ removing electrons
		for (int j=0; j< (int) miscstuff.myRawElectron.size() && _Nele_match<(int) miscstuff.myRawElectron.size(); j++) {
			TLorentzVector *_ele=&miscstuff.myRawElectron[j];
			double _match_dR;
			double _match_dEta;
			double _match_dPhi;

			//std::cout << " processing ele i, ind" << i << "\t" << miscstuff.myEleInd[j] << std::endl;
			int match_stat=MyMatchPhoJet(jetstuff.JetBlockInd.at(i),-1,miscstuff.myEleInd[j],
					fJetBlock,_ele,_jet,_match_dR,_match_dPhi,_match_dEta);

			//std::cout << "jet i, delR, dEta, dPhi= " << i <<  "\t" << _match_dR << "\t" << _match_dEta << "," << _match_dPhi << std::endl;

			if (match_stat==1)
			{
				_Nele_match++;
				jetstuff.Nobj_match.at(i)=jetstuff.Nobj_match.at(i)+1;
				jetstuff.Nele_match.at(i)=jetstuff.Nele_match.at(i)+1;

				//i added this to get matching info correctly, above stores a dummy value for EmInd_match
				// sam -- 03/02/2008
				sam_match.JetInd_match=jetstuff.JetBlockInd.at(i);
				sam_match.EmInd_match= miscstuff.myEleInd[j];
				sam_match.EmObjType_match=1; // electrons
				sam_matchstuff.push_back(sam_match);

				if ((jetstuff.Nele_match.at(i)+jetstuff.Npho_match.at(i))==1)
				{ // filling only first match
					dummymatch.JetInd_match=i;
					dummymatch.EmInd_match=j;
					dummymatch.EmObjType_match=1; // electrons
					matchstuff.push_back(dummymatch);
				}

				// here, I re-calculate jet EmFr
				jetstuff.EmFrCorr.at(i)=CorrectEmFr(jetstuff.EmFrCorr.at(i),miscstuff.myEleEmFr[j],
						jetstuff.Jet_raw_noEMobj.at(i).E(),miscstuff.myRawElectron[j].E());
				// corrected EtaDet is to be calculated before removing EM object
				_etadet_tmp=CorrectEtaDet(_ele,&jetstuff.Jet_raw_noEMobj.at(i),
						jetstuff.EtaDetCorr.at(i),miscstuff.myEleEtaDet[j]); // new line, 10/20/05
				jetstuff.EtaDetCorr.at(i)=_etadet_tmp;  // new line, 10/20/05
				double _jet_phi=jetstuff.Jet_raw.at(i).Phi();
				double _jet_eta=jetstuff.Jet_raw.at(i).Eta();

				if (jetstuff.Jet_raw_noEMobj.at(i).E()>miscstuff.myRawElectron[j].E())
					jetstuff.Jet_raw_noEMobj.at(i)=jetstuff.Jet_raw_noEMobj.at(i)-miscstuff.myRawElectron[j]; // removing raw electron from raw jet
				else jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);

				if (jetstuff.Jet_raw_noEMobj.at(i).Pt()<epsilon) jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon); // new, 05/28/08
				if (jetstuff.Jet_raw_noEMobj.at(i).E()<jetstuff.Jet_raw_noEMobj.at(i).P()
						|| jetstuff.Jet_raw_noEMobj.at(i).E()<epsilon) jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);
				if(jetstuff.Jet_raw.at(i).DeltaR(jetstuff.Jet_raw_noEMobj.at(i))>2.0*(fJetBlock->ConeSize()))
				{
					jetstuff.Jet_raw_noEMobj.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);
				}

			}	//match stat
		} // electron loop


		if ((jetstuff.Nele_match.at(i)+jetstuff.Npho_match.at(i))>0) // correcting jet energy after removing EM objetcs
		{
			TLorentzVector _rawJet=jetstuff.Jet_raw_noEMobj.at(i);
			TLorentzVector _rawJet_tmp;
			_rawJet_tmp=_rawJet;
			float _myEmFr_tmp=jetstuff.EmFrCorr.at(i);
			double corr7=GetCorrection(lev7,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			_myEmFr_tmp=jetstuff.EmFrCorr.at(i); // do this again because EmFr can be potentially changed in GetCorrection
			_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
			double corr6=GetCorrection(lev6,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			_myEmFr_tmp=jetstuff.EmFrCorr.at(i); // do this again because EmFr can be potentially changed in GetCorrection
			_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
			double corr5=GetCorrection(lev5,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			_myEmFr_tmp=jetstuff.EmFrCorr.at(i); // do this again because EmFr can be potentially changed in GetCorrection
			_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
			double corr4=GetCorrection(lev4,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			_myEmFr_tmp=jetstuff.EmFrCorr.at(i); // do this again because EmFr can be potentially changed in GetCorrection
			_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
			double corr1=GetCorrection(lev1,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
			if(corr1<0.0) corr1=epsilon;
			if(corr4<0.0) corr4=epsilon;
			if(corr5<0.0) corr5=epsilon;
			if(corr6<0.0) corr6=epsilon;
			if(corr7<0.0) corr7=epsilon;
			jetstuff.Jet_lev1_noEMobj.at(i)=_rawJet*corr1;
			jetstuff.Jet_lev4_noEMobj.at(i)=_rawJet*corr4;
			jetstuff.Jet_lev5_noEMobj.at(i)=_rawJet*corr5;
			jetstuff.Jet_lev6_noEMobj.at(i)=_rawJet*corr6;
			jetstuff.Jet_lev7_noEMobj.at(i)=_rawJet*corr7;
		}

		if(fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())>=MinJetEta && fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<=MaxJetEta)
		{
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>5.0) jetstuff.myNjet_th5++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>10.0) jetstuff.myNjet_th10++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>15.0) jetstuff.myNjet_th15++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>20.0) jetstuff.myNjet_th20++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>25.0) jetstuff.myNjet_th25++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>30.0) jetstuff.myNjet_th30++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>35.0) jetstuff.myNjet_th35++;
		}
	} // jet loop

	if(_Nele_match!=(int)miscstuff.myRawElectron.size() || _Npho_match!=(int)miscstuff.myRawPhoton.size()) bad_EMjet_match_flag=0;
	for(unsigned int i=0; i<matchstuff.size(); i++)
	{
		FillMatchingHistograms(Hist,jetstuff,miscstuff,matchstuff[i]);
	}

	if (_Npho_match != (int)miscstuff.myRawPhoton.size())
	{
		std::cout << __FILE__ << ":" << __LINE__ << "::JET:: PHO NO MATCH: Npho, Nmatched ="
			<< miscstuff.myRawPhoton.size() << ", " << _Npho_match << "\t- ";
		GetHeaderBlock()->Print();
	}
	if (_Nele_match != (int) miscstuff.myRawElectron.size())
	{
		std::cout << __FILE__ << ":" << __LINE__ << "::JET:: ELE NO MATCH: Nele, Nmatched ="
			<< miscstuff.myRawElectron.size() << ", " << _Nele_match << "\t- ";
		GetHeaderBlock()->Print();
	}


	return;
}
//__________________________________ does my jets, but doesn't remove EM objetcs
void JetFilterModuleV2::DoMyJetNoMatch(TStnJetBlock* fJetBlock, JetStuff &jetstuff) {

	int _jetcone=0;
	if(fabs(fJetBlock->ConeSize()-0.4)<0.01) _jetcone=0;
	if(fabs(fJetBlock->ConeSize()-0.7)<0.01) _jetcone=1;
	if(fabs(fJetBlock->ConeSize()-1.0)<0.01) _jetcone=2;

	fJTC_coneSize=_jetcone;
	//_________________________ Initializing Jet corrections
	CorInit lev1;
	CorInit lev4;
	CorInit lev5;
	CorInit lev6;
	CorInit lev7;
	lev1.level=1;
	lev1.nvx=myNvx_class12;
	lev1.cone=fJTC_coneSize;
	lev1.version=fJTC_version;
	lev1.sys=fJTC_systcode;
	lev1.imode=fJTC_imode;
	lev1.Nrun=myJetRun;
	lev4.level=4;
	lev4.nvx=myNvx_class12;
	lev4.cone=fJTC_coneSize;
	lev4.version=fJTC_version;
	lev4.sys=fJTC_systcode;
	lev4.imode=fJTC_imode;
	lev4.Nrun=myJetRun;
	lev5.level=5;
	lev5.nvx=myNvx_class12;
	lev5.cone=fJTC_coneSize;
	lev5.version=fJTC_version;
	lev5.sys=fJTC_systcode;
	lev5.imode=fJTC_imode;
	lev5.Nrun=myJetRun;
	lev6.level=6;
	lev6.nvx=myNvx_class12;
	lev6.cone=fJTC_coneSize;
	lev6.version=fJTC_version;
	lev6.sys=fJTC_systcode;
	lev6.imode=fJTC_imode;
	lev6.Nrun=myJetRun;
	lev7.level=7;
	lev7.nvx=myNvx_class12;
	lev7.cone=fJTC_coneSize;
	lev7.version=fJTC_version;
	lev7.sys=fJTC_systcode;
	lev7.imode=fJTC_imode;
	lev7.Nrun=myJetRun;

	jetstuff.myNjet=fJetBlock->NJets();
	for(int i=0; i<jetstuff.myNjet; i++)
	{
		TStnJet* jet = fJetBlock->Jet(i);
		jetstuff.JetBlockInd.push_back(i);
		jetstuff.JetNtrk.push_back(jet->NTracks());
		jetstuff.JetNtwr.push_back(jet->NTowers());
		jetstuff.EtaDet.push_back(jet->DetEta());
		jetstuff.EtaDetCorr.push_back(jet->DetEta()); // for a moment, make it the same as "EtaDet"
		jetstuff.EmFrRaw.push_back(jet->Emfr());
		jetstuff.EmFrCorr.push_back(jet->Emfr()); // for a moment, make it the same as "EmFrRaw"

		jetstuff.Nobj_match.push_back(0); // these are packed with zero's for a moment
		jetstuff.Npho_match.push_back(0);
		jetstuff.Nele_match.push_back(0);
		jetstuff.Nmu_match.push_back(0);
		jetstuff.Ntau_match.push_back(0);
		jetstuff.Nbtag_match.push_back(0);

		//__________________ getting raw jets
		TLorentzVector _rawJet;
		TLorentzVector _rawJet_tmp;
		_rawJet.SetPx(jet->Momentum()->Px());
		_rawJet.SetPy(jet->Momentum()->Py());
		_rawJet.SetPz(jet->Momentum()->Pz());
		_rawJet.SetE(jet->Momentum()->E());
		jetstuff.Jet_raw.push_back(_rawJet);
		jetstuff.Jet_raw_noEMobj.push_back(_rawJet); // for a moment, make it the same as "Jet_raw"

		float _myEmFr_tmp=jet->Emfr();
		float _etadet_tmp=jet->DetEta();
		_rawJet_tmp=_rawJet;
		double corr7=GetCorrection(lev7,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
		_myEmFr_tmp=jet->Emfr(); // do this again because EmFr can be potentially changed in GetCorrection
		_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
		double corr6=GetCorrection(lev6,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
		_myEmFr_tmp=jet->Emfr(); // do this again because EmFr can be potentially changed in GetCorrection
		_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
		double corr5=GetCorrection(lev5,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
		_myEmFr_tmp=jet->Emfr(); // do this again because EmFr can be potentially changed in GetCorrection
		_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
		double corr4=GetCorrection(lev4,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);
		_myEmFr_tmp=jet->Emfr(); // do this again because EmFr can be potentially changed in GetCorrection
		_rawJet_tmp=_rawJet; // do this again because _rawJet_tmp can be potentially changed in GetCorrection
		double corr1=GetCorrection(lev1,_rawJet_tmp,_myEmFr_tmp,_etadet_tmp);

		if(corr1<0.0) corr1=0.0;
		if(corr4<0.0) corr4=0.0;
		if(corr5<0.0) corr5=0.0;
		if(corr6<0.0) corr6=0.0;
		if(corr7<0.0) corr7=0.0;
		jetstuff.Jet_lev1.push_back(_rawJet*corr1);
		jetstuff.Jet_lev1_noEMobj.push_back(_rawJet*corr1); // for a moment, make it the same as "Jet_lev1"
		jetstuff.Jet_lev4.push_back(_rawJet*corr4);
		jetstuff.Jet_lev4_noEMobj.push_back(_rawJet*corr4); // for a moment, make it the same as "Jet_lev4"
		jetstuff.Jet_lev5.push_back(_rawJet*corr5);
		jetstuff.Jet_lev5_noEMobj.push_back(_rawJet*corr5); // for a moment, make it the same as "Jet_lev5"
		jetstuff.Jet_lev6.push_back(_rawJet*corr6);
		jetstuff.Jet_lev6_noEMobj.push_back(_rawJet*corr6); // for a moment, make it the same as "Jet_lev6"
		jetstuff.Jet_lev7.push_back(_rawJet*corr7);
		jetstuff.Jet_lev7_noEMobj.push_back(_rawJet*corr7); // for a moment, make it the same as "Jet_lev7"

	}
	return;
}


//--------------------------- 09/20/06 ----------------------------------
//____ This is new function to generate Met due to EM object energy resolution
void JetFilterModuleV2::GenerateMyEMobjMet(CommonStuff miscstuff, int systcode, TVector2 &myEMGenMet) {
	// Resolution: G/Et=sqrt(p0*p0/Et+p1*p1)
	double p0[2]; // resolution params: [0]=CEM; [1]=PEM
	double p1[2];
	double p0_er[2];
	double p1_er[2];
	//______________ CEM parametrization (temporary)
	p0[0]=0.135; // according to hep-ex/0510047 (jes nim)
	p1[0]=0.015;  // according to hep-ex/0510047 (jes nim)
	p0_er[0]=0.01*p0[0]; // ? guesstimate
	p1_er[0]=0.01*p1[0]; // ? guesstimate
	//______________ PEM parametrization (temporary)
	p0[1]=0.16; // according to hep-ex/0510047 (jes nim)
	p1[1]=0.01;  // according to hep-ex/0510047 (jes nim)
	p0_er[1]=0.01*p0[1]; // ? guesstimate
	p1_er[1]=0.01*p1[1]; // ? guesstimate
	//--------------------------------------
	int det_code=0;
	TLorentzVector emmetsum(0.0,0.0,0.0,0.0);
	TLorentzVector emmetsum_devm(0.0,0.0,0.0,0.0);
	TLorentzVector emmetsum_devp(0.0,0.0,0.0,0.0);
	for(unsigned int i=0; i<miscstuff.myRawPhoton.size(); i++)
	{
		if(fabs(miscstuff.myPhoEtaDet.at(i))<=1.1) det_code=0;
		else det_code=1;
		double obj_Pt=miscstuff.myCorrPhoton.at(i).Pt();
		double P_0=(p0[det_code]+systcode*p0_er[det_code])*(p0[det_code]+systcode*p0_er[det_code]);
		double P_1=(p1[det_code]+systcode*p1_er[det_code])*(p1[det_code]+systcode*p1_er[det_code]);
		double sigmaEM=sqrt(P_0/obj_Pt+P_1);
		double scale_factor=gRandom->Gaus(1.0,sigmaEM); // get jet energy scale factor
		if(miscstuff.myRawElectron.size()<2)
		{
			if(i<2) while(scale_factor*miscstuff.myCorrPhoton.at(i).Pt()<13.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
			if(i>=2) while(scale_factor*miscstuff.myCorrPhoton.at(i).Pt()<7.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
		}
		if(miscstuff.myRawElectron.size()>=2) while(scale_factor*miscstuff.myCorrPhoton.at(i).Pt()<7.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
		emmetsum=emmetsum+(1.0-scale_factor)*miscstuff.myCorrPhoton.at(i);
	}
	for(unsigned int i=0; i<miscstuff.myRawElectron.size(); i++)
	{
		if(fabs(miscstuff.myEleEtaDet.at(i))<=1.1) det_code=0;
		else det_code=1;
		double obj_Pt=miscstuff.myCorrElectron.at(i).Pt();
		double P_0=(p0[det_code]+systcode*p0_er[det_code])*(p0[det_code]+systcode*p0_er[det_code]);
		double P_1=(p1[det_code]+systcode*p1_er[det_code])*(p1[det_code]+systcode*p1_er[det_code]);
		double sigmaEM=sqrt(P_0/obj_Pt+P_1);
		double scale_factor=gRandom->Gaus(1.0,sigmaEM); // get jet energy scale factor
		if(miscstuff.myRawElectron.size()<2) while(scale_factor*miscstuff.myCorrElectron.at(i).Pt()<13.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
		if(miscstuff.myRawElectron.size()>=2)
		{
			if(i==0) while(scale_factor*miscstuff.myCorrElectron.at(i).Pt()<20.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
			if(i>0) while(scale_factor*miscstuff.myCorrElectron.at(i).Pt()<10.0) scale_factor=gRandom->Gaus(1.0,sigmaEM);
		}
		emmetsum=emmetsum+(1.0-scale_factor)*miscstuff.myCorrElectron.at(i);
	}
	myEMGenMet.Set(emmetsum.Px(),emmetsum.Py());

	return;
}

//_________ This function adds Met to closest jet, the output is "new" JetLev6
void JetFilterModuleV2::MyNewJetWithMet(JetStuff &jetstuff, const CommonStuff& commstuff)
{

	assert (jetstuff.newJetLev6.size() == 0 && "newJetLev6 is not 0 to begin with");
	//________________________________________temporaly replaces detector jets with HADron jets -- 02-11-2009
	if (bUseHadJets)
	{
		if (GetHeaderBlock()->McFlag())
		{
			// first stuff newJetLev6 with L6Detjet to get the correct size of the
			// vector 
			for (unsigned int i=0; i < jetstuff.Jet_lev6_noEMobj.size(); i++) 
			{	
				// packing "jetstuff.newJetLev6" with det jets temporily to get the right size
				jetstuff.newJetLev6.push_back(jetstuff.Jet_lev6_noEMobj.at(i));
			}


			std::vector<TLorentzVector> HadJets;

			// need to remove the photon/electron from the had jet list now!
			// as I am starting from HAD jets
			for (int i=0; i < fHadJetBlockClu04->NJets(); ++i)  		//loop over HAD jets to find match
			{
				TStnJet *hadjet = fHadJetBlockClu04->Jet(i); 
				TLorentzVector hjetvec(*hadjet->Momentum());
				bool emobj = false;

				for (unsigned int ipho = 0; ipho < commstuff.myCorrPhoton.size(); ++ipho)
				{
					if (hjetvec.DeltaR(commstuff.myCorrPhoton.at(ipho)) < 0.1)
					{
						emobj = true;
						break;
					}
				}
				for (unsigned int iele = 0; iele < commstuff.myCorrElectron.size(); ++iele)
				{
					if (hjetvec.DeltaR(commstuff.myCorrElectron.at(iele)) < 0.1)
					{
						emobj = true;
						break;
					}
				}

				if (emobj) continue;
				HadJets.push_back(hjetvec);

			} // hadron jet loop

			if (HadJets.size() > jetstuff.Jet_lev6_noEMobj.size())
			{
				std::cout << "Event with more HAD jets than DET jets! ";
				fHeaderBlock->Print();
				mycounter.evtsMoreHad++;
			}

			// if there are more had jets than the det jets, use the most energetic had jets only
			// order of the jets does not matter, as they will be resorted after smearing
			if (HadJets.size() > jetstuff.Jet_lev6_noEMobj.size())
			{
				if (HadJets.size())
				{
					std::vector<TLorentzVector> vjets = HadJets;
					std::set<std::pair<float,unsigned> > PtInd;
					std::set<std::pair<float,unsigned> >::iterator myIt;

					//sort the had jet in Pt
					for (unsigned int i=0; i < vjets.size(); ++i)
					{
						std::pair<float, unsigned> jetPtInd(vjets.at(i).Pt(), i);
						PtInd.insert(jetPtInd);
					}

					HadJets.clear();
					myIt = PtInd.end();
					while (myIt != PtInd.begin())
					{
						myIt--;
						HadJets.push_back(vjets.at(myIt->second));
					}

				}
			}


			// now replace the newJetLev6 by HAD jets and zero out remaining, if any
			int replaced = 0;
			for (unsigned int i=0; i < jetstuff.newJetLev6.size(); i++)
			{
				if (i < HadJets.size())
				{
					jetstuff.newJetLev6.at(i) = HadJets.at(replaced);
					++replaced;
				} else
				{
					//now stuff in zeros for unused(no had jets left to fill) slots in newJetLev6
					if (jetstuff.newJetLev6.at(i).Pt() > 0.1) // if this is already assigned the 1e-20 no need to redo it!
					{
						float epsilon = 1.0E-20;
						float _jet_phi = jetstuff.newJetLev6.at(i).Phi();
						float _jet_eta = jetstuff.newJetLev6.at(i).Eta();
						jetstuff.newJetLev6.at(i).SetPtEtaPhiM(epsilon,_jet_eta,_jet_phi,0.3*epsilon);
					}
				}
			}


		} else
		{
			StdOut(__FILE__,__LINE__,3,"Request to Hadron jets ignored. This is not a MC dataset. Resetting to false.");
		}


	} else {
		//use detector level jets with met addition

		for (unsigned int i=0; i < jetstuff.Jet_lev6_noEMobj.size(); i++) // packing "jetstuff.newJetLev6" with "old" jets
		{
			jetstuff.newJetLev6.push_back(jetstuff.Jet_lev6_noEMobj.at(i));
		}
		
		int iMetAddJetIndex = -1;
		
		if (MEtAddMethod()==DO_NOT_ADD_MET) { /*does nothing*/ }
		else if (MEtAddMethod()==TOMOST_MISMEASURED_JET)
		{
			iMetAddJetIndex = LeastMismeasuredJet(jetstuff.newJetLev6,jetstuff.myMETcorr_th15);
		}
		else if (MEtAddMethod()==TO_JET_CLOSEST_TO_MET)
		{
			iMetAddJetIndex = JetClosestToMet(jetstuff.newJetLev6,jetstuff.myMETcorr_th15);
		}
		else 
		{
			StdOut(__FILE__,__LINE__,3,"Method of MET addition not specified. exiting!");
			exit(1);
		}

		// if a suitable jet is found MEt will be added.
		if (iMetAddJetIndex>=0)
		{
			double phi_met    = jetstuff.myMETcorr_th15.Phi();
			double phi_jet    = jetstuff.Jet_lev6_noEMobj[iMetAddJetIndex].Phi();
			double eta_jet    = jetstuff.Jet_lev6_noEMobj[iMetAddJetIndex].Eta();
			double pt_jet     = jetstuff.Jet_lev6_noEMobj[iMetAddJetIndex].Pt();
			double raw_pt_jet = jetstuff.Jet_raw_noEMobj[iMetAddJetIndex].Pt();
			double e_jet      = jetstuff.Jet_lev6_noEMobj[iMetAddJetIndex].E();
			double dphi       = fabs(TVector2::Phi_mpi_pi(phi_jet-phi_met));

			// here we project the MEt on to the Jet axis.
			double pt_scale   = (jetstuff.myMETcorr_th15.Mod()*cos(dphi)+pt_jet)/pt_jet;
			if (pt_jet<15.0)  pt_scale = pt_scale-1.0+raw_pt_jet/pt_jet;  // I forgot what this is for
			if (pt_scale<0.0) pt_scale = 1.0E-6;
			
			jetstuff.newJetLev6[iMetAddJetIndex].SetPtEtaPhiE(pt_jet*pt_scale,eta_jet,phi_jet,e_jet*pt_scale);
			//std::cout << jetstuff.newJetLev6[iMetAddJetIndex].Pt() << std::endl;
			
			assert (jetstuff.newJetLev6.at(iMetAddJetIndex).E()>0 && "MetAdded jet energy is now <0!");
		}

	}

}

//___ This function returns upper limit of allowed jet energy fluctuation.
double JetFilterModuleV2::MyIntegralUpLimit(const double jet_E, double eta_det, int stat_code, int syst_code)
{
	double value = 5.0; // maximum allowed: E(det)/E(true)-1<=5.0
	double parJER[5];
	parJER[0] = MyJer_meanG (jet_E, eta_det, 0, syst_code, JERForGaussMean()); // for now systcode is used as jer_stat_code
	parJER[1] = MyJer_sigmaG(jet_E, eta_det, 0, syst_code, JERForGaussSigma());
	parJER[2] = MyJer_mpvL  (jet_E, eta_det, 0, syst_code, JERForLandauMPV());
	parJER[3] = MyJer_sigmaL(jet_E, eta_det, 0, syst_code, JERForLandauSigma());
	parJER[4] = MyJer_normG (jet_E, eta_det, 0, syst_code, JERForGaussNorm());

	double lim_L = 5.0;
	double lim_G = 5.0;
	//------ calculating limit for Landau part
	//--- This limit is roughly equivalent to 1.0E-6 significance level of JER Landau part if it is found in range [-1.0;10.0]
	if (fabs(1+parJER[2]-3.7*parJER[3])>0.0) lim_L=(3.7*parJER[3]-parJER[2])/(1+parJER[2]-3.7*parJER[3]);
	if (lim_L>5.0 || lim_L<parJER[3]) lim_L=5.0;
	//------ calculating limit for Gaussian part taking into account relative normalization of Gaussian
	if(parJER[4]>0.0 && parJER[4]<=1.0)
	{
		double sig_tmp=-log10(parJER[4]/(1+parJER[4]));
		double sig=4.2-sig_tmp; // 4.2 roughly corresponds to 4*Sigma significance level (~1.0E-5)
		if(sig<0.0) sig=0.0;
		double n_sigma=0.6807+1.035*sig-0.05574*sig*sig; // crude parametrization of N_sigma vs. significance
		if(fabs(1+parJER[0]-n_sigma*parJER[1])>0.0) lim_G=(n_sigma*parJER[1]-parJER[0])/(1+parJER[0]-n_sigma*parJER[1]);
		else // just to make a "smooth" transition between Landau and Gauss
		{
			if(parJER[0]<1.0 && parJER[0]>0.0) lim_G=1.0/parJER[0]-1.0;
			else lim_G=-1.0;
		}
		if(lim_G>5.0 || lim_G<parJER[1]) lim_G=5.0;
	}
	//------ deciding what limit should be used: Landau or Gauss
	if(parJER[4]>0.0) // gauss is non-zero
	{
		if(lim_L<lim_G) value=lim_G; // Landau lim smaller
		else value=lim_L;// Landau lim larger
	}
	else value=lim_L; // gauss is zero
	if(value>5.0) value=5.0;
	return value;
}

//___ This function clears params for MetProb
void JetFilterModuleV2::ClearMetProbStuff(MetProbStuff& stuff) {
	stuff.max_MetProbInteg=1.0;
	return;
}

//___ returns relative position of EM shower inside CEM tower
double JetFilterModuleV2::MyInCemTowerEta(double eta_det) {
	double rel_eta=-10.0;
	double theta_boundary[12]={90.0,82.526,75.297,68.516,62.310,56.735,51.790,47.436,43.614,40.261,36.822,33.524};
	double eta_l, eta_r, d_eta;
	double _pi=TMath::Pi();
	for(int i=1; i<12; i++)
	{
		eta_l=fabs(log(tan(_pi*theta_boundary[i-1]/(2.0*180.0))));
		eta_r=fabs(log(tan(_pi*theta_boundary[i]/(2.0*180.0))));
		d_eta=eta_r-eta_l;
		if(fabs(eta_det)>eta_l && fabs(eta_det)<=eta_r)
		{
			rel_eta=fabs(fabs(eta_det)-eta_l)/d_eta;
			break;
		}
	}
	return rel_eta;
}


/////////////////////////////////////////////////////////////////////////////////
//This is the original basic MetSig calc version that Sasha has. I am 
// using it again to verify that I have not broken the MEtModel
/////////////////////////////////////////////////////////////////////////////////
void JetFilterModuleV2::InitMetProbStuff(MetProbStuff& stuff,JetStuff jetstuff,
							TVector2 myMetVec, int systcode, int systcodeUncl) 
{

	//std::cout << "_____________________________________________________________\n "
	//			<< __LINE__ << "::In Function " << __FUNCTION__ 
	//			<< "[SYSTCODE = " << systcode << "][SYSTCODEUNCL = " << systcodeUncl  
	//			<< "] stuff.max_MetProbInteg = " <<  stuff.max_MetProbInteg
	//			<< std::endl;
	double jet_EtCut = 15.0;
	double _met = myMetVec.Mod();

	//-------------------------- Getting unclustered energy parameters
	double _meanX  = 0.0;
	double _sigmaX = 0.0;
	double _meanY  = 0.0;
	double _sigmaY = 0.0;
	double _normX  = 0.0;
	double _scaleX = 0.0;
	double _normY  = 0.0;
	double _scaleY = 0.0;
	GetMyUnclMetResolution(jetstuff.mySumEtCorr_th15, fJTC_imode, fUnclParamSwitch,
								systcodeUncl, _sigmaX, _sigmaY, _meanX, _meanY, _normX, 
								_normY, _scaleX, _scaleY);

	//------------------------------------------------------------
	// calculating precision for multi-D integration
	//______ uncl. energy part
	TF1* ueFun = new TF1("ue",UnclMetResolution,-1000.0,1000.0,4);
	ueFun->SetParameter(0,_meanX);
	ueFun->SetParameter(1,_sigmaX);
	ueFun->SetParameter(2,_normX);
	ueFun->SetParameter(3,_scaleX);
	
	double ue_lim = 5.0 * _scaleX * _sigmaX > _met ? 5.0 * _scaleX * _sigmaX : 2.0*_met;
	double scaleF = ueFun->Integral(-1.0 * ue_lim, 1.0 * ue_lim);

	double test_prob;
	if (scaleF > 0.0) test_prob = 1.0 - (ueFun->Integral(-1.0*_met, _met)) / scaleF;
	else test_prob = 1.0;

	if (test_prob <= 0.0) test_prob = 1.0E-19;
	if (test_prob > 1.0) test_prob = 1.0; 

	stuff.max_MetProbInteg = stuff.max_MetProbInteg * (1.0 - test_prob); // calculating upper limit on MetProb
	ueFun->SetParameter(0,_meanY);
	ueFun->SetParameter(1,_sigmaY);
	ueFun->SetParameter(2,_normY);
	ueFun->SetParameter(3,_scaleY);
	ue_lim=5.0*_scaleY*_sigmaY>_met ? 5.0*_scaleY*_sigmaY : 2.0*_met;
	scaleF=ueFun->Integral(-1.0*ue_lim,1.0*ue_lim);
	if(scaleF>0.0) test_prob=1.0-(ueFun->Integral(-1.0*_met,_met))/scaleF;
	else test_prob=1.0;
	if(test_prob<=0.0) test_prob=1.0E-19;
	if(test_prob>1.0) test_prob=1.0; 
	stuff.max_MetProbInteg=stuff.max_MetProbInteg*(1.0-test_prob); // calculating upper limit on MetProb
	delete ueFun;


	//std::cout << "\t Uncl MetSig :: stuff.max_MetProbInteg = " <<  stuff.max_MetProbInteg
	//			<< std::endl;
	//------------------------------------------------------------
	//______ jet part
	
	for (unsigned int i=0; i<jetstuff.newJetLev6.size(); i++)
	{
		if (jetstuff.newJetLev6.at(i).Pt()>3.0 && fabs(jetstuff.newJetLev6.at(i).Eta())<=MaxJetEta)
		{
			double parJER[5];
			parJER[0] = MyJer_meanG(jetstuff.newJetLev6.at(i).E(), jetstuff.EtaDetCorr.at(i), 
											0, systcode, JERForGaussMean()); // for now systcode is used as jer_stat_code
			parJER[1] = MyJer_sigmaG(jetstuff.newJetLev6.at(i).E(), jetstuff.EtaDetCorr.at(i), 
											0, systcode, JERForGaussSigma());
			parJER[2] = MyJer_mpvL(jetstuff.newJetLev6.at(i).E(), jetstuff.EtaDetCorr.at(i), 
											0, systcode, JERForLandauMPV());
			parJER[3] = MyJer_sigmaL(jetstuff.newJetLev6.at(i).E(), jetstuff.EtaDetCorr.at(i), 
											0, systcode, JERForLandauSigma());
			parJER[4] = MyJer_normG(jetstuff.newJetLev6.at(i).E(), jetstuff.EtaDetCorr.at(i), 
											0, systcode, JERForGaussNorm());

			double jer_limit = MyIntegralUpLimit(jetstuff.newJetLev6.at(i).E(),
											jetstuff.EtaDetCorr.at(i), 0, systcode);

			TF1* JerFun = new TF1("jer",MyJER,-1.0,jer_limit,5);
			for (int par_ind = 0; par_ind < 5; par_ind++)
			{
				JerFun->SetParameter(par_ind,parJER[par_ind]);
			}

	//		std::cout << "\tFor Jet " << i << "] E, Eta = " 
	//					<< jetstuff.newJetLev6.at(i).E() << ", "
	//					<< jetstuff.EtaDetCorr.at(i) << std::endl;
			
			// JER total integral
			double tot_jetInt = JerFun->Integral(-1.0, jer_limit);
	//		std::cout << "\t\t Total JER Integral = " << tot_jetInt << " (-1.0, " << jer_limit 
	//					<< ")" << std::endl;
			
			double off_set = JerFun->GetMaximumX(-1.0,1.0); // returns mpv of JER function;
	//		std::cout << "\t\t off_set = " << off_set << std::endl;

			// I forgot do this for this function when I was first looking at the effect of
			// this off_set effect 08-08-2009  //check with Sasha
			if (ZeroJERoffset()) off_set = 0;

			//i do not understand this????
			if (jetstuff.newJetLev6.at(i).Pt() - _met > jet_EtCut)
			{
	//			std::cout << "\tCASE 1" << std::endl;
			//i do not understand this????
				double x1 = -_met / jetstuff.newJetLev6.at(i).Pt() + off_set;
				double x2 = _met / jetstuff.newJetLev6.at(i).Pt() + off_set;

				if (x1 < -1.0) x1 = -1.0;
				if (x2 > jer_limit) x2=jer_limit;  
				
				test_prob = 1.0 - (JerFun->Integral(x1, x2)) / tot_jetInt;
	//			std::cout << "\t\t\t test_prob = " << test_prob << " (x1=" << x1 << ", x2=" << x1 << ")" << std::endl;

				if (test_prob <= 0.0) test_prob = 1.0E-19;
				if (test_prob > 1.0)  test_prob = 1.0;

				stuff.max_MetProbInteg = stuff.max_MetProbInteg * (1.0 - test_prob); // calculating upper limit on MetProb

	//			std::cout << "\t\t\tTotProb, testPro = " <<  stuff.max_MetProbInteg 
	//						<< ", " << 1 - test_prob << std::endl;
			}
			else
			{
	//			std::cout << "\tCASE 2" << std::endl;
				double x1 = jet_EtCut / jetstuff.newJetLev6.at(i).Pt() - 1.0 + off_set;
				double x2 = _met / jetstuff.newJetLev6.at(i).Pt() + off_set;
				
				if (x2 < x1) x2 = x1;
				if (x2 > jer_limit) x2 = jer_limit;  
				
				//if (tot_jetInt == 0)
				//{
				//	test_prob = 0;
				//	std::cout << "\t\t EST PROB == 0 JER Evals" << std::endl;
				//	JerFun->Print();
				//	for (int j = 0; j < 10; ++j)
				//	{
				//		std::cout << "\t\t\t Eval " << j << " = " << JerFun->Eval(i) << std::endl;  
				//	}
				//} else {
					test_prob = 1.0 - (JerFun->Integral(-1.0, x2)) / tot_jetInt;
				//}

	//			if (std::isnan(test_prob))
	//			{
	//				std::cout << "\t\t test_prob is NAN "<< std::endl; 
	//				JerFun->Print();
	//			}
				
	//			std::cout << "\t\t\t test_prob = " << test_prob << " (x1=" << x1 << ", x2=" << x1 << ")" << std::endl;
				
				if (test_prob <= 0.0) test_prob = 1.0E-19;
				if (test_prob > 1.0)  test_prob = 1.0;
				
				if (jetstuff.newJetLev6.at(i).Pt() + _met > jet_EtCut) 
				{
					stuff.max_MetProbInteg = stuff.max_MetProbInteg * (1.0 - test_prob); // upper lim on MetProb
				}
				
	//			std::cout << "\t\t\tTotProb, 1-testPro = " <<  stuff.max_MetProbInteg 
	//						<< ", " << 1 - test_prob << std::endl;
			}
			delete JerFun;
		} else {
	//		std::cout << "\t\t Skipping jet " << i << std::endl;
		}
	}

	
	if (std::isnan(stuff.max_MetProbInteg))
	{
		std::cout << __FUNCTION__ << ":" << __LINE__ << ":MetSig Calculation is NAN! (stuff.max_MetProbInteg=NAN). Exiting!" << std::endl; 
		exit (1);
	}

	///std::cout << "____________________________ Returning from " << __FUNCTION__ 
		//		<< ":: stuff.max_MetProbInteg = " <<  stuff.max_MetProbInteg
		//		<< std::endl;
	
	return;
}



//___ This function returns input parameters for MetSig correction
double JetFilterModuleV2::MetSigCorrPAR(int ind, int njet15, int isMC, int sample_id) {
	double param=0.0;
	double par_arr[5][6]; // fixed array length on 06/20/08
	int jetind=0;
	if(njet15==0) {
		jetind=0;
		// I must protect myself from using njet=0 as I do not have calibration for it.
		if (sample_id == 3 | sample_id ==4)
		{
			if (PrintLevel()>10)
			{
				std::cout << __FILE__ << ":"<< __FUNCTION__ <<":" << __LINE__ << ":: Sam, you did not do calibration for njet==0. Why are you trying to use it? I am giving you absolute zeros!" << std::endl;
			}
		}
	}
	if(njet15==1) jetind=1;
	if(njet15==2) jetind=2;
	if(njet15==3) jetind=3;
	if(njet15>=4) jetind=4;


	if(ind<0 || ind>5) return 0.0; // fixed on 06/23/08
	if(isMC==0) // for MC
	{
		if(sample_id==0) // Pythia ggX signal
		{
			// I added zeros for njet>=3 bins for these to be safe! sam-10-27-2008
			if (PrintLevel()>PRN_WARN)
			{
				std::cout << __FILE__ << "::" << __LINE__ << ":: WARNING. You have not updated calibration for this sample!" << std::endl;
			}
			// Pythia calibration params are the same as in the e-log entry on 06/23/08
			par_arr[0][0]=0.008;  par_arr[1][0]=0.1586;  par_arr[2][0]=0.1438;   par_arr[3][0]=0.0; par_arr[4][0]=0.0;
			par_arr[0][1]=-3.088; par_arr[1][1]=1.515;   par_arr[2][1]=1.4; 	  par_arr[3][1]=0.0; par_arr[4][1]=0.0;
			par_arr[0][2]=-2.625;	par_arr[1][2]=-0.104;  par_arr[2][2]=-0.081;   par_arr[3][2]=0.0; par_arr[4][2]=0.0;
			par_arr[0][3]=0.1807;	par_arr[1][3]=0.00076; par_arr[2][3]=0.000179; par_arr[3][3]=0.0; par_arr[4][3]=0.0;
			par_arr[0][4]=2.247;  par_arr[1][4]=0.507;   par_arr[2][4]=0.143;    par_arr[3][4]=0.0; par_arr[4][4]=0.0;
			// 	  par_arr[0][5]=6.17;   par_arr[1][5]=15.11;   par_arr[2][5]=12.57; par_arr[3][5]=0.0; par_arr[4][5]=0.0;
			par_arr[0][5]=6.17;   par_arr[1][5]=12.0;   par_arr[2][5]=9.0;       par_arr[3][6]=0.0; par_arr[4][6]=0.0;
		}
		if(sample_id==2) // Pythia cem-cem Z->ee
		{
			// I added zeros for njet>=3 bins for these to be safe! sam-10-27-2008
			std::cout << __FILE__ << "::" << __LINE__ << ":: WARNING. You have not updated clibration for this sample!" << std::endl;
			// fixed some entires on 06/20/08 (made the same as in 01/20/08 -- "Re-calibration of MetSig")
			par_arr[0][0]=0.00986; par_arr[1][0]=0.6469;  par_arr[2][0]=0.5971;  par_arr[3][0]=0.0; par_arr[4][0]=0.0;
			par_arr[0][1]=-2.19; 	 par_arr[1][1]=1.507;   par_arr[2][1]=1.382;   par_arr[3][1]=0.0; par_arr[4][1]=0.0;
			par_arr[0][2]=-2.10;	 par_arr[1][2]=-0.1073; par_arr[2][2]=-0.0857; par_arr[3][2]=0.0; par_arr[4][2]=0.0;
			par_arr[0][3]=0.1774;	 par_arr[1][3]=0.00564; par_arr[2][3]=0.00121; par_arr[3][3]=0.0; par_arr[4][3]=0.0;
			par_arr[0][4]=2.168;   par_arr[1][4]=0.6112;  par_arr[2][4]=0.359;   par_arr[3][4]=0.0; par_arr[4][4]=0.0;
			// 	  par_arr[0][5]=5.97;    par_arr[1][5]=10.65;   par_arr[2][5]=15.01; par_arr[3][5]=0.0; par_arr[4][5]=0.0;
			par_arr[0][5]=5.97;    par_arr[1][5]=12.0;   par_arr[2][5]=9.0;      par_arr[3][5]=0.0; par_arr[4][5]=0.0;
		}
		if(sample_id==3) // Pythia pho-jet(Pt_hat>22) sam:updated 10-24-2008
		{
			par_arr[0][0]=0; par_arr[1][0]=0.2225;	par_arr[2][0]=0.177;	par_arr[3][0]=0.01424;	par_arr[4][0]=0.02735;
			par_arr[0][1]=0; par_arr[1][1]=0.9034;	par_arr[2][1]=0.527;	par_arr[3][1]=0.06331;	par_arr[4][1]=-2.589;
			par_arr[0][2]=0; par_arr[1][2]=-0.1795;	par_arr[2][2]=-0.2684;	par_arr[3][2]=-0.3897;	par_arr[4][2]=-1.183;
			par_arr[0][3]=0; par_arr[1][3]=0.004741;	par_arr[2][3]=0.008691;	par_arr[3][3]=0.00006814;	par_arr[4][3]=0.3082;
			par_arr[0][4]=0; par_arr[1][4]=0.5934;	par_arr[2][4]=0.6857;	par_arr[3][4]=0.2228;	par_arr[4][4]=1.364;
			par_arr[0][5]=0; par_arr[1][5]=14;	par_arr[2][5]=13;	par_arr[3][5]=6.5;	par_arr[4][5]=6;
		}
	}
	if(isMC==1) // for DATA
	{
		if(sample_id==0) // ggX signal
		{
			// I added zeros for njet>=3 bins for these to be safe! sam-10-27-2008
			std::cout << __FILE__ << "::" << __LINE__ << ":: WARNING. You have not updated clibration for this sample!" << std::endl;
			// fixed some entires on 06/20/08 (made the same as in 01/20/08 -- "Re-calibration of MetSig")
			par_arr[0][0]=0.00693; par_arr[1][0]=0.6237; par_arr[2][0]=0.587; par_arr[3][0]=0.0; par_arr[4][0]=0.0;
			par_arr[0][1]=-3.37; 	 par_arr[1][1]=1.387;  par_arr[2][1]=1.316; par_arr[3][1]=0.0; par_arr[4][1]=0.0;
			par_arr[0][2]=-2.79;	 par_arr[1][2]=-0.180; par_arr[2][2]=-0.127; par_arr[3][2]=0.0; par_arr[4][2]=0.0;
			par_arr[0][3]=0.1838;	 par_arr[1][3]=0.0083; par_arr[2][3]=0.0022; par_arr[3][3]=0.0; par_arr[4][3]=0.0;
			par_arr[0][4]=2.26;    par_arr[1][4]=0.624;  par_arr[2][4]=0.36;   par_arr[3][4]=0.0; par_arr[4][4]=0.0;
			// 	  par_arr[0][5]=5.59;    par_arr[1][5]=11.81;  par_arr[2][5]=14.75;  par_arr[3][5]=0.0; par_arr[4][5]=0.0;
			par_arr[0][5]=5.59;    par_arr[1][5]=12.0;  par_arr[2][5]=9.0;   par_arr[3][5]=0.0; par_arr[4][5]=0.0;
		}
		if(sample_id==1) // ggX sideband
		{
			// I added zeros for njet>=3 bins for these to be safe! sam-10-27-2008
			std::cout << __FILE__ << "::" << __LINE__ << ":: WARNING. You have not updated clibration for this sample!" << std::endl;
			par_arr[0][0]=0.00849; par_arr[1][0]=0.5939; par_arr[2][0]=0.570; par_arr[3][0]=0.0; par_arr[4][0]=0.0;
			par_arr[0][1]=-3.04; 	 par_arr[1][1]=1.295;  par_arr[2][1]=1.273; par_arr[3][1]=0.0; par_arr[4][1]=0.0;
			par_arr[0][2]=-2.61;	 par_arr[1][2]=-0.227; par_arr[2][2]=-0.141; par_arr[3][2]=0.0; par_arr[4][2]=0.0;
			par_arr[0][3]=0.1813;	 par_arr[1][3]=0.0161; par_arr[2][3]=0.0051; par_arr[3][3]=0.0; par_arr[4][3]=0.0;
			par_arr[0][4]=2.277;   par_arr[1][4]=0.724;  par_arr[2][4]=0.50; par_arr[3][4]=0.0; par_arr[4][4]=0.0;
			// 	  par_arr[0][5]=5.39;    par_arr[1][5]=8.69;   par_arr[2][5]=11.65; par_arr[3][5]=0.0; par_arr[4][5]=0.0;
			par_arr[0][5]=5.39;    par_arr[1][5]=12.0;   par_arr[2][5]=9.0; par_arr[3][5]=0.0; par_arr[4][5]=0.0;
		}
		if(sample_id==2) // cem-cem Z->ee
		{
			// I added zeros for njet>=3 bins for these to be safe! sam-10-27-2008
			std::cout << __FILE__ << "::" << __LINE__ << ":: WARNING. You have not updated clibration for this sample!" << std::endl;
			// fixed some entires on 06/20/08 (made the same as in 01/20/08 -- "Re-calibration of MetSig")
			par_arr[0][0]=0.0065;  par_arr[1][0]=0.6395; par_arr[2][0]=0.611;   par_arr[3][0]=0.0; par_arr[4][0]=0.0;
			par_arr[0][1]=-4.03; 	 par_arr[1][1]=1.504;  par_arr[2][1]=1.466;   par_arr[3][1]=0.0; par_arr[4][1]=0.0;
			par_arr[0][2]=-3.70;	 par_arr[1][2]=-0.140; par_arr[2][2]=-0.041;  par_arr[3][2]=0.0; par_arr[4][2]=0.0;
			par_arr[0][3]=0.1794;	 par_arr[1][3]=0.0142; par_arr[2][3]=0.00056; par_arr[3][3]=0.0; par_arr[4][3]=0.0;
			par_arr[0][4]=2.151;   par_arr[1][4]=0.675;  par_arr[2][4]=0.17;     par_arr[3][4]=0.0; par_arr[4][4]=0.0;
			// 	  par_arr[0][5]=6.11;    par_arr[1][5]=12.65;  par_arr[2][5]=14.35; par_arr[3][5]=0.0; par_arr[4][5]=0.0;
			par_arr[0][5]=6.11;    par_arr[1][5]=12.0;   par_arr[2][5]=9.0;     par_arr[3][5]=0.0;  par_arr[4][5]=0.0;
		}
		if(sample_id==3) // Photon Data signal pho+jet  sam:updated 10-24-2008
		{
			par_arr[0][0]=0; par_arr[1][0]=0.5047;	par_arr[2][0]=0.4173;	par_arr[3][0]=0.1766;	par_arr[4][0]=0.464;
			par_arr[0][1]=0; par_arr[1][1]=1.149;	par_arr[2][1]=0.7897;	par_arr[3][1]=0.4634;	par_arr[4][1]=-0.4958;
			par_arr[0][2]=0; par_arr[1][2]=-0.0883;	par_arr[2][2]=-0.1734;	par_arr[3][2]=-0.2588;	par_arr[4][2]=-0.54;
			par_arr[0][3]=0; par_arr[1][3]=0.005567;	par_arr[2][3]=0.006247;	par_arr[3][3]=0.000398;	par_arr[4][3]=0.1854;
			par_arr[0][4]=0; par_arr[1][4]=0.7214;	par_arr[2][4]=0.7715;	par_arr[3][4]=0.585;	par_arr[4][4]=1.305;
			par_arr[0][5]=0; par_arr[1][5]=14;	par_arr[2][5]=14;	par_arr[3][5]=12;	par_arr[4][5]=8.5;

		}
		if(sample_id==4) // Photon Data sidband  pho+jet  sam:added 10-24-2008
		{
			par_arr[0][0]=0; par_arr[1][0]=0.6314;	par_arr[2][0]=0.545;	par_arr[3][0]=0.4539;	par_arr[4][0]=0.5585;
			par_arr[0][1]=0; par_arr[1][1]=1.131;	par_arr[2][1]=0.8589;	par_arr[3][1]=0.4924;	par_arr[4][1]=-0.01166;
			par_arr[0][2]=0; par_arr[1][2]=-0.09311;	par_arr[2][2]=-0.1503;	par_arr[3][2]=-0.2475;	par_arr[4][2]=-0.3923;
			par_arr[0][3]=0; par_arr[1][3]=0.005976;	par_arr[2][3]=0.004662;	par_arr[3][3]=0.0001962;	par_arr[4][3]=0.0002891;
			par_arr[0][4]=0; par_arr[1][4]=0.7618;	par_arr[2][4]=0.7923;	par_arr[3][4]=0.4541;	par_arr[4][4]=0.4259;
			par_arr[0][5]=0; par_arr[1][5]=15;	par_arr[2][5]=12;	par_arr[3][5]=8.5;	par_arr[4][5]=8.5;

		}
	}

	param=par_arr[jetind][ind];
	return param;
}

//______________ This functio returns corrected (true) MetSig
//______________ No correction is done for fSampleId==10 (when RawMetSig is requested)
double JetFilterModuleV2::MetSigCorrection(double rawSig) {
	double corr=1.0;
	if(fSampleID>=0 && fSampleID<=4)
	{
		float cut1=0.0;
		float cut2=20.0;
		TF1 *fitFun=new TF1("fit",rawMetSigFunc,cut1,cut2,5);
		int njet=jet04stuff.myNjet_th15;
		fitFun->SetParameter(0,MetSigCorrPAR(0,njet,fJTC_imode,fSampleID));
		fitFun->SetParameter(1,MetSigCorrPAR(1,njet,fJTC_imode,fSampleID));
		fitFun->SetParameter(2,MetSigCorrPAR(2,njet,fJTC_imode,fSampleID));
		fitFun->SetParameter(3,MetSigCorrPAR(3,njet,fJTC_imode,fSampleID));
		fitFun->SetParameter(4,MetSigCorrPAR(4,njet,fJTC_imode,fSampleID));

		double bin_int=0.0;
		double bin_int0=0.0;
		double last_bin_slope=1.0;
		double tot_int=fitFun->Integral(0.0,20.0);
		double par5=MetSigCorrPAR(5,njet,fJTC_imode,fSampleID); // fixed par5 on 06/20/08

		bin_int=fitFun->Integral(0.0,rawSig);
		if(tot_int>0.0 && par5>0.0) bin_int=bin_int/tot_int;
		else
		{
			delete fitFun;
			return rawSig;
		}
		if(bin_int>0.0 && bin_int<1.0 && rawSig<par5 && rawSig>0.0) corr=-log10(1-bin_int)/rawSig;
		else
		{
			bin_int0=fitFun->Integral(0.0,par5);
			bin_int0=bin_int0/tot_int;
			last_bin_slope=-log10(1-bin_int0)/par5;
			corr=last_bin_slope;
		}
		if(corr<0.0) corr=1.0;
		delete fitFun;
	}
	return rawSig*corr;
}


//___ This function returns my Met significance for all events
double JetFilterModuleV2::MyTotalMetSignificance(JetStuff jetstuff, TVector2 myMetVec, int systcode, int systcodeUncl) {

	ClearMetProbStuff(metsigstuff); //clears metprob stuff (reset MetProb to 1)
	InitMetProbStuff(metsigstuff,jetstuff,myMetVec,systcode,systcodeUncl); //inits params for MetProb

	double metSig_est = 0.0;
	if (metsigstuff.max_MetProbInteg > 0.0 && metsigstuff.max_MetProbInteg < 1.0)
	{
		metSig_est = -1.0 * log10(1.0 - metsigstuff.max_MetProbInteg);
	}
	//std::cout << __FUNCTION__ << "::metSig_est = " << metSig_est << std::endl;
	if (metsigstuff.max_MetProbInteg<=0.0) metSig_est = 0.0;
	if (metsigstuff.max_MetProbInteg>=1.0) metSig_est = 19.0;
	if ((19.0 - metSig_est) >0.1) metSig_est = MetSigCorrection(metSig_est);
	//std::cout << __FUNCTION__ << "::Corr metSig_est = " << metSig_est << std::endl;
	return metSig_est;
}

//___ filling final histogram for MetSig correction
void JetFilterModuleV2::FinalMetSigCorrHisto(MetStudyHisto_t& Hist) {
	int nbins=Hist.fMetSigCalib_estimate->GetNbinsX();
	//int Nevnt=(int) Hist.fMetSigCalib_estimate->GetEntries();	//i added the cast to remove the constant compiler warning! --sam-03-15-2008
	double Nevnt= Hist.fMetSigCalib_estimate->GetEntries();	//this is probably the right type that has enough space to hold a large number! --sam-11-03-2008
	double bin_val=0.0;
	double sum=0.0;
	double corr=0.0;
	double corr_err=0.0;
	for(int i=0; i<nbins; i++)
	{
		bin_val=Hist.fMetSigCalib_estimate->GetBinContent(i+1);
		if(bin_val>0.0) sum=sum+bin_val/(1.0*Nevnt);
		if(sum<1.0)
		{
			corr=-1.0*log10(1.0-sum);
			if(bin_val>0.0 && Nevnt>0) corr_err=sqrt(bin_val)/(log(10.0)*(1.0-sum)*Nevnt);
			else corr_err=0.0;
			if(corr_err>0.0)
			{
				Hist.fMetSig_corrFactor->SetBinContent(i+1,corr);
				Hist.fMetSig_corrFactor->SetBinError(i+1,corr_err);
			}
		}
		else
		{
			Hist.fMetSig_corrFactor->SetBinContent(i+1,log10(1.0*Nevnt));
			Hist.fMetSig_corrFactor->SetBinError(i+1,1.0/log(10.0));
		}
	}

	//---------------------- calculating efficiency of MetSig cut for toy MC
	for(int i=1; i<5; i++)
	{
		Hist.fMetSigToy_eff[i-1]->Divide(Hist.fMetSigToy_Met[i],Hist.fMetSigToy_Met[0]);
		nbins=Hist.fMetSigToy_eff[i-1]->GetNbinsX();
		for(int j=0; j<nbins; j++)
		{
			double val1=Hist.fMetSigToy_Met[i]->GetBinContent(j+1);
			double val2=Hist.fMetSigToy_Met[0]->GetBinContent(j+1);
			double err=0.0;
			if(val2>0.0) err=sqrt(val1*(1.0-val1/val2))/val2; // this efficiency is wrong. need a precise formula for efficiency
			Hist.fMetSigToy_eff[i-1]->SetBinError(j+1,err);
		}
	}
	return;
}


//___ This is new function to generate Met; it takes care of unclustered & jet components of Met.
//___ cleanup histograms for generated met are also filled here.
void JetFilterModuleV2::GenerateMyTotalMet(JetStuff& jetstuff, CommonStuff miscstuff, MetCleanupHisto_t& Hist,
		int systcode, int systcodeUncl, int metcode, int rnd_seed, TVector2 &myGenMet)
{

	if(systcode!=0 || systcodeUncl!=0)
	{
		//std::cout << __FUNCTION__ << ":seed=" << rnd_seed << std::endl;
		gRandom->SetSeed(rnd_seed);
	}

	std::vector<TLorentzVector> jet_smeared; // new as of 04/03/07
	jet_smeared.clear(); // new as of 04/03/07
	jetstuff.smear_factor.clear();
	myGenMet.Set(0.0,0.0);

	//___________ generating met due to jets
	double SumEtJet_smear = 0.0; 		// U.E.+M.I. for jets above threshold after smearing
	double SumEtJet       = 0.0;     // U.E.+M.I. for jets above threshold before smearing
	double SumEtRawJet    = 0.0;
	double SumEtRawJet_smear = 0.0;
	double Pt_cut[5] = {2000.0,5.0,10.0,15.0,20.0};
	std::vector<int> jet_ind; 			// index of smeared jets fluctuated above Pt-cut
	double cut;
	double scale_factor = 1.0;
	TLorentzVector jetmetsum(0.0,0.0,0.0,0.0);
	TVector2 myJetMet(0.0,0.0);

	int NjetCut=0;
	if (metcode==0) NjetCut=jetstuff.myNjet;
	if (metcode==1) NjetCut=jetstuff.myNjet_th5;
	if (metcode==2) NjetCut=jetstuff.myNjet_th10;
	if (metcode==3) NjetCut=jetstuff.myNjet_th15;
	if (metcode==4) NjetCut=jetstuff.myNjet_th20;
	
	if (PrintLevel()>10)
	{
		std::cout << "--------- " << __FUNCTION__ << " ---- DEBUG INFO" << std::endl;
		fHeaderBlock->Print();
		std::cout << "metcode = " << metcode << std::endl;
		std::cout << "NjetCut, njet, njet5, njet10, njet15, mjet20 = " << NjetCut << "," << jetstuff.myNjet 
			<< "," <<  jetstuff.myNjet_th5 << ","<< jetstuff.myNjet_th10 << "," << jetstuff.myNjet_th15
			<< "," << jetstuff.myNjet_th20 << std::endl;
	}
	
	if (metcode>0 && metcode<5)    // moved the old '&& njetcut>0' cut further down (inside for loop). so the smeared jets will be
 	{										// generated for every event! not just for event with njetcut>0
		cut = Pt_cut[metcode];
		for(unsigned int i=0; i<jetstuff.newJetLev6.size(); i++) // "jetstuff.Jet_lev6_noEMobj" is replaced by "jetstuff.newJetLev6" (11/02/06)
		{
			jetstuff.smear_factor.push_back(1.0);
			//on line below I use exactly the same eta cut as in MET correction
			if (jetstuff.newJetLev6.at(i).Pt()>3.0 
				 && fabs(jetstuff.newJetLev6.at(i).Eta())<=MaxJetEta) // "jetstuff.Jet_lev6_noEMobj" is replaced by "jetstuff.newJetLev6" (11/02/06)
			{
				double parJER[5];
				// for now systcode is used as jer_stat_code
				parJER[0] = MyJer_meanG(jetstuff.newJetLev6.at(i).E(),
									jetstuff.EtaDetCorr.at(i), 0, systcode, JERForGaussMean()); 
				parJER[1] = MyJer_sigmaG(jetstuff.newJetLev6.at(i).E(),
									jetstuff.EtaDetCorr.at(i), 0, systcode, JERForGaussSigma());
				parJER[2] = MyJer_mpvL(jetstuff.newJetLev6.at(i).E(),
									jetstuff.EtaDetCorr.at(i), 0, systcode, JERForLandauMPV());
				parJER[3] = MyJer_sigmaL(jetstuff.newJetLev6.at(i).E(),
									jetstuff.EtaDetCorr.at(i), 0, systcode, JERForLandauSigma());
				parJER[4] = MyJer_normG(jetstuff.newJetLev6.at(i).E(),
									jetstuff.EtaDetCorr.at(i), 0, systcode, JERForGaussNorm());

				double jer_limit = MyIntegralUpLimit(jetstuff.newJetLev6.at(i).E(), 
												jetstuff.EtaDetCorr.at(i),0,systcode);

				TF1* JerFun = new TF1("jer", MyJER, -1.0, jer_limit, 5);
				
				//------ setting Gaussian
				JerFun->SetParameter(0,parJER[0]);
				JerFun->SetParameter(1,parJER[1]);
				//------ setting Landau
				JerFun->SetParameter(2,parJER[2]);
				JerFun->SetParameter(3,parJER[3]);
				//------ setting Gaussian normalization
				JerFun->SetParameter(4,parJER[4]);
				// for data we do not know the mean position. but for MC we know. So we do not need
				// this offset. In data this is about 3% shift in mean which according to Sasha
				// was not a problem for his analysis. It could go to a max of ~5% for low Et jets.
				// see his preblessing webtalk on 10/23/2008 slide 19
				double off_set=JerFun->GetMaximumX(-1.0,1.0); // returns mpv of JER function;
				//scale_factor=1.0-off_set+JerFun->GetRandom(); // get jet energy scale factor; default=assume jet scale is 1.0
				//std::cout << __FUNCTION__ << ":" << __LINE__ << std::endl;
				double jer_rand = JerFun->GetRandom(); // get jet energy scale factor; default=assume jet scale is 1.0
				//std::cout << "\tFor Jet " << i << "] E, Eta = " 
				//		<< jetstuff.newJetLev6.at(i).E() << ", "
				//		<< jetstuff.EtaDetCorr.at(i) << std::endl;
				//std::cout << "\tjer_rand = " << jer_rand << std::endl;

				// set off_set to zero temporily to study the effect of it - 05-04-2009
				if (ZeroJERoffset()) off_set = 0;
				//std::cout << "\t off_set = " << off_set << std::endl;

				scale_factor = 1.0 - off_set + jer_rand; // get jet energy scale factor; default=assume jet scale is 1.0
				
				if (debug) std::cout << __LINE__ << "::" << __FUNCTION__ << ":: scale factor initial=" << scale_factor;

				if(scale_factor < 0.0) {
					scale_factor = 0.0;
				}
				//std::cout << "\tscale_factor = " << scale_factor << std::endl;
				//std::cout << __FUNCTION__ << ":" << __LINE__ << std::endl;
				delete JerFun;
				jetstuff.smear_factor[i] = scale_factor;
				jet_smeared.push_back(scale_factor * jetstuff.newJetLev6.at(i)); // new as of 04/03/07

				// Had to move this njet cut here from above, ' if(metcode>0 && metcode<5) '
				// I need to have smear factors generated when there is no jet requirement.
				// (when there is no jet cut, like in gg+X analysis)
				// If not, I am not able to generate smeared jets correctly. sam-06-22-2009
				if (NjetCut>0)
				{
					// "jetstuff.Jet_lev6_noEMobj" is replaced by "jetstuff.newJetLev6" (11/02/06)
					if((scale_factor*jetstuff.newJetLev6.at(i)).Pt()>cut) 
					{
						jet_ind.push_back(i); // temporary; commented out on 06/27/06 for studies
						jetmetsum = jetmetsum + (1.0 - scale_factor)*jetstuff.newJetLev6.at(i); // "jetstuff.Jet_lev6_noEMobj" is replaced by "jetstuff.newJetLev6" (11/02/06)
						// SumEtJet_smear= Sum(U.E.+M.I.)
						SumEtJet_smear=SumEtJet_smear+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
							+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
							-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
							-jetstuff.Jet_lev6_noEMobj.at(i).Pt(); // Jet_lev6_noEMobj added on 06/26/06
						SumEtRawJet_smear=SumEtRawJet_smear+jetstuff.Jet_raw_noEMobj.at(i).Pt();
					}
					if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>cut)
					{
						SumEtJet=SumEtJet+jetstuff.Jet_lev5_noEMobj.at(i).Pt()
							+jetstuff.Jet_lev1_noEMobj.at(i).Pt()
							-jetstuff.Jet_lev4_noEMobj.at(i).Pt()
							-jetstuff.Jet_lev6_noEMobj.at(i).Pt(); // Jet_lev6_noEMobj added on 06/26/06
						SumEtRawJet=SumEtRawJet+jetstuff.Jet_raw_noEMobj.at(i).Pt();
					}
				} //njet cut
			}

		}
		myJetMet.Set(jetmetsum.Px(),jetmetsum.Py());
	}


	//_____debug :: sam

	if (debug)
	{
		assert (jetstuff.newJetLev6.size() == jetstuff.smear_factor.size() && "smear factors do not match to njets.");
		std::cout << "============" << __FUNCTION__ << "==================" << std::endl;
		for(unsigned int i=0; i<jetstuff.newJetLev6.size(); i++)
		{
			std::cout << "\t" << i << "\tjet= " << jetstuff.newJetLev6.at(i).Pt() << "\tsmearfac= " << jetstuff.smear_factor.at(i) << std::endl;
		}
	}	
	//-------

	
	//_____________ generating unclustered component of Met
	double x;
	double y;
	double _meanX=0.0;
	double _sigmaX=0.0;
	double _meanY=0.0;
	double _sigmaY=0.0;
	double _normX=0.0;
	double _scaleX=0.0;
	double _normY=0.0;
	double _scaleY=0.0;
	double sumEt=0.0;
	double sumEt_raw=0.0;
	//--------------------- this chunk of code is modified on 11/09/06
	if(metcode==1) sumEt_raw=jetstuff.mySumEtCorr_th5-SumEtJet+SumEtRawJet;
	if(metcode==2) sumEt_raw=jetstuff.mySumEtCorr_th10-SumEtJet+SumEtRawJet;
	if(metcode==3) sumEt_raw=jetstuff.mySumEtCorr_th15-SumEtJet+SumEtRawJet;
	if(metcode==4) sumEt_raw=jetstuff.mySumEtCorr_th20-SumEtJet+SumEtRawJet;
	if(metcode<1 || metcode>4) sumEt_raw=miscstuff.mySumEt_raw;
	sumEt=sumEt_raw+SumEtJet_smear-SumEtRawJet_smear;
	//------------------------------------------------ end of 11/09/06 part
	if(sumEt<0.0) sumEt=0.0;
	GetMyUnclMetResolution(sumEt,fJTC_imode,fUnclParamSwitch,systcodeUncl,_sigmaX,_sigmaY,_meanX,_meanY,_normX,_normY,_scaleX,_scaleY);
	//------- new part (05/24/07)
	TF1* UnclMetFun=new TF1("fit",UnclMetResolution,-200.0,200.0,4);
	//------ setting 1st Gaussian
	UnclMetFun->SetParameter(0,_meanX);
	UnclMetFun->SetParameter(1,_sigmaX);
	//------ setting 2nd Gaussian
	UnclMetFun->SetParameter(2,_normX);
	UnclMetFun->SetParameter(3,_scaleX);
	x=UnclMetFun->GetRandom();
	UnclMetFun->SetParameter(0,_meanY);
	UnclMetFun->SetParameter(1,_sigmaY);
	//------ setting 2nd Gaussian
	UnclMetFun->SetParameter(2,_normY);
	UnclMetFun->SetParameter(3,_scaleY);
	y=UnclMetFun->GetRandom();
	delete UnclMetFun;
	//------- end of new part (05/24/07)
	myGenMet.Set(x,y);
	myGenMet=myGenMet+myJetMet;
	//   //_______________ new part, 09/20/06
	//   TVector2 myEMGenMet(0.0,0.0);
	//   GenerateMyEMobjMet(miscstuff,systcode,myEMGenMet);
	//   myGenMet=myGenMet+myEMGenMet;
	//___________________________________ filling GenMet cleanup histograms
	if(fAnalysisMode==0) FillGenMetCleanupHistograms(Hist,jetstuff,miscstuff,myGenMet,jetstuff.Jet_lev6_noEMobj); // histo for generated met cleanup studies

	return;
}

//______________________________________________________________________________________
void JetFilterModuleV2::Smear_newJetL6FirstTime(JetStuff& jetstuff,
							const int systcode, const int metcode,
							const int rnd_seed)
{
// To study the effect of double smearing I am smearing the jets for pseudo experiment
// twice. This is the first smearing. Then these jets will go through the
// normal routine I'll have to redo the njet count too! No the smeared Njet count is done
// after 2nd smearing.


//	DumpJets(__FUNCTION__,__LINE__,jetstuff,20);

	
	if (! DoubleSmearJets()) return;

	if (systcode!=0)
	{
		//std::cout << __FUNCTION__ << ":seed=" << rnd_seed << std::endl;
		gRandom->SetSeed(rnd_seed);
	}

	double Pt_cut[5] = {2000.0,5.0,10.0,15.0,20.0};
	double cut;

	if(metcode>0 && metcode<5)
	{
		cut=Pt_cut[metcode];
		for(unsigned int i=0; i<jetstuff.newJetLev6.size(); i++) // "jetstuff.Jet_lev6_noEMobj" is replaced by "jetstuff.newJetLev6" (11/02/06)
		{
			//on line below I use exactly the same eta cut as in MET correction
			if(jetstuff.newJetLev6.at(i).Pt()>3.0 && fabs(jetstuff.newJetLev6.at(i).Eta())<=MaxJetEta) // "jetstuff.Jet_lev6_noEMobj" is replaced by "jetstuff.newJetLev6" (11/02/06)
			{
				std::cout << __FUNCTION__ << "::" << jetstuff.newJetLev6.at(i).E() << std::endl;
				double parJER[5];
				parJER[0]=MyJer_meanG(jetstuff.newJetLev6.at(i).E(),jetstuff.EtaDetCorr.at(i),0,systcode, JERForGaussMean()); // for now systcode is used as jer_stat_code
				parJER[1]=MyJer_sigmaG(jetstuff.newJetLev6.at(i).E(),jetstuff.EtaDetCorr.at(i),0,systcode, JERForGaussSigma());
				parJER[2]=MyJer_mpvL(jetstuff.newJetLev6.at(i).E(),jetstuff.EtaDetCorr.at(i),0,systcode, JERForLandauMPV());
				parJER[3]=MyJer_sigmaL(jetstuff.newJetLev6.at(i).E(),jetstuff.EtaDetCorr.at(i),0,systcode, JERForLandauSigma());
				parJER[4]=MyJer_normG(jetstuff.newJetLev6.at(i).E(),jetstuff.EtaDetCorr.at(i),0,systcode, JERForGaussNorm());
				double jer_limit=MyIntegralUpLimit(jetstuff.newJetLev6.at(i).E(),jetstuff.EtaDetCorr.at(i),0,systcode);
				TF1* JerFun=new TF1("jer",MyJER,-1.0,jer_limit,5);
				//------ setting Gaussian
				JerFun->SetParameter(0,parJER[0]);
				JerFun->SetParameter(1,parJER[1]);
				//------ setting Landau
				JerFun->SetParameter(2,parJER[2]);
				JerFun->SetParameter(3,parJER[3]);
				//------ setting Gaussian normalization
				JerFun->SetParameter(4,parJER[4]);
				double off_set=JerFun->GetMaximumX(-1.0,1.0); // returns mpv of JER function;
				//scale_factor=1.0-off_set+JerFun->GetRandom(); // get jet energy scale factor; default=assume jet scale is 1.0
				if (ZeroJERoffset()) off_set = 0;
				std::cout << __FUNCTION__ << ":" << __LINE__ << std::endl;
				double jer_rand = JerFun->GetRandom(); // get jet energy scale factor; default=assume jet scale is 1.0
				std::cout << __FUNCTION__ << ":" << __LINE__ << std::endl;
				double scale_factor=1.0-off_set+jer_rand; // get jet energy scale factor; default=assume jet scale is 1.0

				if (scale_factor<0.0)
				{
					scale_factor=0.0;
				}
				delete JerFun;
				jetstuff.newJetLev6.at(i) *= scale_factor;
			}
		}
	}


//	std::cout << "----------- AFTER  " <<std::endl;
//	DumpJets(__FUNCTION__,__LINE__,jetstuff,20);
	
}


//---------------------------------------------------
//Hacked by Sam - 28-12-2008
//now check and set if there are any smeared jets that is above the given threshold
//must be called everytime after GenerateMyTotalMet
//---------------------------------------------------
void JetFilterModuleV2::GenerateSmearedJets(JetStuff& jetstuff)
{

	if (debug)
	{
		std::cout << "============" << __FUNCTION__ 
					<< " DEBUG INFO BEGIN ==================" << std::endl;
		std::cout << "newJetLev.size , smear_factors.size = " 
					<< jetstuff.newJetLev6.size() << "\t"
					<< jetstuff.smear_factor.size() << std::endl;

		for(unsigned int i=0; i<jetstuff.newJetLev6.size(); i++)
		{
			std::cout << "\t" << i << "\tjet= " << jetstuff.newJetLev6.at(i).Pt() << "\tsmearfac= " << jetstuff.smear_factor.at(i) << std::endl;
		}
	}
	
	assert (jetstuff.newJetLev6.size() == jetstuff.smear_factor.size() && "number of smear factors does not match the number of jets!");
	jetstuff.smeared_newJetLev6_noEMObj.clear();
	for(unsigned int i=0; i<jetstuff.newJetLev6.size(); i++)
	{
		//is the smear factor and newJetLev ordered, if not this may be wrong
		//if (jetstuff.newJetLev6.at(i).Pt()>3.0 && fabs(jetstuff.newJetLev6.at(i).Eta())<=MaxJetEta)
		jetstuff.smeared_newJetLev6_noEMObj.push_back(jetstuff.smear_factor.at(i) *jetstuff.newJetLev6.at(i));
	}
	ReorderSmearedJets(jetstuff.smeared_newJetLev6_noEMObj);

	std::vector<int> smeared_Njets;
	GetNjets(jetstuff.smeared_newJetLev6_noEMObj, smeared_Njets);
	jetstuff.smeared_Njet_th0 = smeared_Njets.at(0);
	jetstuff.smeared_Njet_th5 = smeared_Njets.at(1);
	jetstuff.smeared_Njet_th10 = smeared_Njets.at(2);
	jetstuff.smeared_Njet_th15 = smeared_Njets.at(3);
	jetstuff.smeared_Njet_th20 = smeared_Njets.at(4);
	jetstuff.smeared_Njet_th25 = smeared_Njets.at(5);
	jetstuff.smeared_Njet_th30 = smeared_Njets.at(6);
	jetstuff.smeared_Njet_th35 = smeared_Njets.at(7);

	if (debug)
	{
		std::cout <<  "smeared Njets, original njet = " << jetstuff.smeared_Njet_th0 << "\t" << jetstuff.myNjet << std::endl;
		if (jetstuff.myNjet_th15 != jetstuff.smeared_Njet_th15) std::cout << "\t njet, sm_njet=" << jetstuff.myNjet_th15 << "," << jetstuff.smeared_Njet_th15<< std::endl;

		{
			std::cout << "============" << __FUNCTION__ << "==================" << std::endl;
			for(unsigned int i=0; i<jetstuff.smeared_newJetLev6_noEMObj.size(); i++)
			{
				std::cout << "\t" << i << "\tjet= " << jetstuff.smeared_newJetLev6_noEMObj.at(i).Pt() << "\tsmearfac= " << jetstuff.smear_factor.at(i) << std::endl;
			}
		}

		std::cout << "============" << __FUNCTION__ 
					<< " DEBUG INFO ENDS ==================" << std::endl;

	}


}

//_______________ ZEROes met_results arrays; to be called in BeginJob
void JetFilterModuleV2::CleanMetResults(MetResults &metstuff) {

	for(int i=0; i<3; i++)
	{
		for(int j=0; j<10; j++)
		{
			metstuff.ana20Njet_dt[i][j]=0;
			metstuff.ana25Njet_dt[i][j]=0;
			metstuff.ana30Njet_dt[i][j]=0;
			metstuff.ana35Njet_dt[i][j]=0;
			metstuff.ana40Njet_dt[i][j]=0;
			metstuff.ana45Njet_dt[i][j]=0;
			metstuff.ana50Njet_dt[i][j]=0;
			metstuff.ana75Njet_dt[i][j]=0;
			metstuff.ana100Njet_dt[i][j]=0;
			metstuff.ana150Njet_dt[i][j]=0;
			for(int l=0; l<20; l++)
			{
				metstuff.ana20Njet_bg[l][i][j]=0;
				metstuff.ana25Njet_bg[l][i][j]=0;
				metstuff.ana30Njet_bg[l][i][j]=0;
				metstuff.ana35Njet_bg[l][i][j]=0;
				metstuff.ana40Njet_bg[l][i][j]=0;
				metstuff.ana45Njet_bg[l][i][j]=0;
				metstuff.ana50Njet_bg[l][i][j]=0;
				metstuff.ana75Njet_bg[l][i][j]=0;
				metstuff.ana100Njet_bg[l][i][j]=0;
				metstuff.ana150Njet_bg[l][i][j]=0;
			}

			metstuff.Met20Njet_dt[i][j]=0;
			metstuff.Met25Njet_dt[i][j]=0;
			metstuff.Met30Njet_dt[i][j]=0;
			metstuff.Met35Njet_dt[i][j]=0;
			metstuff.Met40Njet_dt[i][j]=0;
			metstuff.Met45Njet_dt[i][j]=0;
			metstuff.Met50Njet_dt[i][j]=0;
			metstuff.Met75Njet_dt[i][j]=0;
			metstuff.Met100Njet_dt[i][j]=0;
			metstuff.Met150Njet_dt[i][j]=0;
			metstuff.Met20Njet_def[i][j]=0;
			metstuff.Met25Njet_def[i][j]=0;
			metstuff.Met30Njet_def[i][j]=0;
			metstuff.Met35Njet_def[i][j]=0;
			metstuff.Met40Njet_def[i][j]=0;
			metstuff.Met45Njet_def[i][j]=0;
			metstuff.Met50Njet_def[i][j]=0;
			metstuff.Met75Njet_def[i][j]=0;
			metstuff.Met100Njet_def[i][j]=0;
			metstuff.Met150Njet_def[i][j]=0;

			metstuff.Met20Njet_gen[i][j]=0.0;
			metstuff.Met25Njet_gen[i][j]=0.0;
			metstuff.Met30Njet_gen[i][j]=0.0;
			metstuff.Met35Njet_gen[i][j]=0.0;
			metstuff.Met40Njet_gen[i][j]=0.0;
			metstuff.Met45Njet_gen[i][j]=0.0;
			metstuff.Met50Njet_gen[i][j]=0.0;
			metstuff.Met75Njet_gen[i][j]=0.0;
			metstuff.Met100Njet_gen[i][j]=0.0;
			metstuff.Met150Njet_gen[i][j]=0.0;
			metstuff.Met20Njet_genstat[i][j]=0.0;
			metstuff.Met25Njet_genstat[i][j]=0.0;
			metstuff.Met30Njet_genstat[i][j]=0.0;
			metstuff.Met35Njet_genstat[i][j]=0.0;
			metstuff.Met40Njet_genstat[i][j]=0.0;
			metstuff.Met45Njet_genstat[i][j]=0.0;
			metstuff.Met50Njet_genstat[i][j]=0.0;
			metstuff.Met75Njet_genstat[i][j]=0.0;
			metstuff.Met100Njet_genstat[i][j]=0.0;
			metstuff.Met150Njet_genstat[i][j]=0.0;
			metstuff.Met20Njet_gensyst[i][j]=0.0;
			metstuff.Met25Njet_gensyst[i][j]=0.0;
			metstuff.Met30Njet_gensyst[i][j]=0.0;
			metstuff.Met35Njet_gensyst[i][j]=0.0;
			metstuff.Met40Njet_gensyst[i][j]=0.0;
			metstuff.Met45Njet_gensyst[i][j]=0.0;
			metstuff.Met50Njet_gensyst[i][j]=0.0;
			metstuff.Met75Njet_gensyst[i][j]=0.0;
			metstuff.Met100Njet_gensyst[i][j]=0.0;
			metstuff.Met150Njet_gensyst[i][j]=0.0;
		}
	}
	return;
}

void JetFilterModuleV2::MyMetEventCount(JetStuff jetstuff, CommonStuff miscstuff, int metcode, MetResults &metstuff) {
	double met=0.0;
	int njet15=0;
	int njet20=0;
	int njet25=0;
	//int njet30=0; // I think I need to complete this section for these two cases
	//int njet35=0;

	if(metcode==1) met=jetstuff.myMETcorr_th5.Mod();
	if(metcode==2) met=jetstuff.myMETcorr_th10.Mod();
	if(metcode==3) met=jetstuff.myMETcorr_th15.Mod();
	if(metcode==4) met=jetstuff.myMETcorr_th20.Mod();
	if(metcode<1 || metcode>4) met=miscstuff.myMET_raw.Mod();
	njet15= (jetstuff.myNjet_th15 < 10) ? jetstuff.myNjet_th15 : 9;
	njet20= (jetstuff.myNjet_th20 < 10) ? jetstuff.myNjet_th20 : 9;
	njet25= (jetstuff.myNjet_th25 < 10) ? jetstuff.myNjet_th25 : 9;

	if(met>20)
	{
		metstuff.Met20Njet_dt[0][njet15]++;
		metstuff.Met20Njet_dt[1][njet20]++;
		metstuff.Met20Njet_dt[2][njet25]++;
	}
	if(met>25)
	{
		metstuff.Met25Njet_dt[0][njet15]++;
		metstuff.Met25Njet_dt[1][njet20]++;
		metstuff.Met25Njet_dt[2][njet25]++;
	}
	if(met>30)
	{
		metstuff.Met30Njet_dt[0][njet15]++;
		metstuff.Met30Njet_dt[1][njet20]++;
		metstuff.Met30Njet_dt[2][njet25]++;
	}
	if(met>35)
	{
		metstuff.Met35Njet_dt[0][njet15]++;
		metstuff.Met35Njet_dt[1][njet20]++;
		metstuff.Met35Njet_dt[2][njet25]++;
	}
	if(met>40)
	{
		metstuff.Met40Njet_dt[0][njet15]++;
		metstuff.Met40Njet_dt[1][njet20]++;
		metstuff.Met40Njet_dt[2][njet25]++;
	}
	if(met>45)
	{
		metstuff.Met45Njet_dt[0][njet15]++;
		metstuff.Met45Njet_dt[1][njet20]++;
		metstuff.Met45Njet_dt[2][njet25]++;
	}
	if(met>50)
	{
		metstuff.Met50Njet_dt[0][njet15]++;
		metstuff.Met50Njet_dt[1][njet20]++;
		metstuff.Met50Njet_dt[2][njet25]++;
	}
	if(met>75)
	{
		metstuff.Met75Njet_dt[0][njet15]++;
		metstuff.Met75Njet_dt[1][njet20]++;
		metstuff.Met75Njet_dt[2][njet25]++;
	}
	if(met>100)
	{
		metstuff.Met100Njet_dt[0][njet15]++;
		metstuff.Met100Njet_dt[1][njet20]++;
		metstuff.Met100Njet_dt[2][njet25]++;
	}
	if(met>150)
	{
		metstuff.Met150Njet_dt[0][njet15]++;
		metstuff.Met150Njet_dt[1][njet20]++;
		metstuff.Met150Njet_dt[2][njet25]++;
	}
	return;
}

//______________ counts generated events above Met cut
void JetFilterModuleV2::MyMetModelEventCount(JetStuff jetstuff, TVector2 myGenMet_def,MetResults &metstuff) {
	int njet15=0;
	int njet20=0;
	int njet25=0;

	double met_def=myGenMet_def.Mod();

	njet15= (jetstuff.myNjet_th15 < 10) ? jetstuff.myNjet_th15 : 9;
	njet20= (jetstuff.myNjet_th20 < 10) ? jetstuff.myNjet_th20 : 9;
	njet25= (jetstuff.myNjet_th20 < 10) ? jetstuff.myNjet_th25 : 9;

	//_____________ default parametrization
	if(met_def>20.0 && met_def<2000.0)
	{
		metstuff.Met20Njet_def[0][njet15]++;
		metstuff.Met20Njet_def[1][njet20]++;
		metstuff.Met20Njet_def[2][njet25]++;
	}
	if(met_def>25.0 && met_def<2000.0)
	{
		metstuff.Met25Njet_def[0][njet15]++;
		metstuff.Met25Njet_def[1][njet20]++;
		metstuff.Met25Njet_def[2][njet25]++;
	}
	if(met_def>30.0 && met_def<2000.0)
	{
		metstuff.Met30Njet_def[0][njet15]++;
		metstuff.Met30Njet_def[1][njet20]++;
		metstuff.Met30Njet_def[2][njet25]++;
	}
	if(met_def>35.0 && met_def<2000.0)
	{
		metstuff.Met35Njet_def[0][njet15]++;
		metstuff.Met35Njet_def[1][njet20]++;
		metstuff.Met35Njet_def[2][njet25]++;
	}
	if(met_def>40.0 && met_def<2000.0)
	{
		metstuff.Met40Njet_def[0][njet15]++;
		metstuff.Met40Njet_def[1][njet20]++;
		metstuff.Met40Njet_def[2][njet25]++;
	}
	if(met_def>45.0 && met_def<2000.0)
	{
		metstuff.Met45Njet_def[0][njet15]++;
		metstuff.Met45Njet_def[1][njet20]++;
		metstuff.Met45Njet_def[2][njet25]++;
	}
	if(met_def>50.0 && met_def<2000.0)
	{
		metstuff.Met50Njet_def[0][njet15]++;
		metstuff.Met50Njet_def[1][njet20]++;
		metstuff.Met50Njet_def[2][njet25]++;
	}
	if(met_def>75.0 && met_def<2000.0)
	{
		metstuff.Met75Njet_def[0][njet15]++;
		metstuff.Met75Njet_def[1][njet20]++;
		metstuff.Met75Njet_def[2][njet25]++;
	}
	if(met_def>100.0 && met_def<2000.0)
	{
		metstuff.Met100Njet_def[0][njet15]++;
		metstuff.Met100Njet_def[1][njet20]++;
		metstuff.Met100Njet_def[2][njet25]++;
	}
	if(met_def>150.0 && met_def<2000.0)
	{
		metstuff.Met150Njet_def[0][njet15]++;
		metstuff.Met150Njet_def[1][njet20]++;
		metstuff.Met150Njet_def[2][njet25]++;
	}
	return;
}


//______________ prints out a comparison of Met in Data to Met Model predictions
void JetFilterModuleV2::MyMetResults(MetResults &metstuff) {

	if(Npoints>0)
	{
		double dummy=0.0;
		double dummy1=0.0;
		double dummy2=0.0;
		//double syst_J1=0.0; //commented out bcos not used
		//double syst_J2=0.0;
		//double syst_U1=0.0;
		//double syst_U2=0.0;
		double syst_tmp20[3][10];
		double syst_tmp25[3][10];
		double syst_tmp30[3][10];
		double syst_tmp35[3][10];
		double syst_tmp40[3][10];
		double syst_tmp45[3][10];
		double syst_tmp50[3][10];
		double syst_tmp75[3][10];
		double syst_tmp100[3][10];
		double syst_tmp150[3][10];
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<10; j++)
			{
				if(fAnalysisMode==1)
				{
					//--------- new event count
					metstuff.Met20Njet_dt[i][j]=metstuff.ana20Njet_dt[i][j];
					metstuff.Met25Njet_dt[i][j]=metstuff.ana25Njet_dt[i][j];
					metstuff.Met30Njet_dt[i][j]=metstuff.ana30Njet_dt[i][j];
					metstuff.Met35Njet_dt[i][j]=metstuff.ana35Njet_dt[i][j];
					metstuff.Met40Njet_dt[i][j]=metstuff.ana40Njet_dt[i][j];
					metstuff.Met45Njet_dt[i][j]=metstuff.ana45Njet_dt[i][j];
					metstuff.Met50Njet_dt[i][j]=metstuff.ana50Njet_dt[i][j];
					metstuff.Met75Njet_dt[i][j]=metstuff.ana75Njet_dt[i][j];
					metstuff.Met100Njet_dt[i][j]=metstuff.ana100Njet_dt[i][j];
					metstuff.Met150Njet_dt[i][j]=metstuff.ana150Njet_dt[i][j];

					metstuff.Met20Njet_gen[i][j]=(1.0*metstuff.ana20Njet_bg[0][i][j])/Npoints;
					metstuff.Met25Njet_gen[i][j]=(1.0*metstuff.ana25Njet_bg[0][i][j])/Npoints;
					metstuff.Met30Njet_gen[i][j]=(1.0*metstuff.ana30Njet_bg[0][i][j])/Npoints;
					metstuff.Met35Njet_gen[i][j]=(1.0*metstuff.ana35Njet_bg[0][i][j])/Npoints;
					metstuff.Met40Njet_gen[i][j]=(1.0*metstuff.ana40Njet_bg[0][i][j])/Npoints;
					metstuff.Met45Njet_gen[i][j]=(1.0*metstuff.ana45Njet_bg[0][i][j])/Npoints;
					metstuff.Met50Njet_gen[i][j]=(1.0*metstuff.ana50Njet_bg[0][i][j])/Npoints;
					metstuff.Met75Njet_gen[i][j]=(1.0*metstuff.ana75Njet_bg[0][i][j])/Npoints;
					metstuff.Met100Njet_gen[i][j]=(1.0*metstuff.ana100Njet_bg[0][i][j])/Npoints;
					metstuff.Met150Njet_gen[i][j]=(1.0*metstuff.ana150Njet_bg[0][i][j])/Npoints;

					metstuff.Met20Njet_genstat[i][j]=sqrt(1.0*metstuff.ana20Njet_bg[0][i][j])/Npoints;
					metstuff.Met25Njet_genstat[i][j]=sqrt(1.0*metstuff.ana25Njet_bg[0][i][j])/Npoints;
					metstuff.Met30Njet_genstat[i][j]=sqrt(1.0*metstuff.ana30Njet_bg[0][i][j])/Npoints;
					metstuff.Met35Njet_genstat[i][j]=sqrt(1.0*metstuff.ana35Njet_bg[0][i][j])/Npoints;
					metstuff.Met40Njet_genstat[i][j]=sqrt(1.0*metstuff.ana40Njet_bg[0][i][j])/Npoints;
					metstuff.Met45Njet_genstat[i][j]=sqrt(1.0*metstuff.ana45Njet_bg[0][i][j])/Npoints;
					metstuff.Met50Njet_genstat[i][j]=sqrt(1.0*metstuff.ana50Njet_bg[0][i][j])/Npoints;
					metstuff.Met75Njet_genstat[i][j]=sqrt(1.0*metstuff.ana75Njet_bg[0][i][j])/Npoints;
					metstuff.Met100Njet_genstat[i][j]=sqrt(1.0*metstuff.ana100Njet_bg[0][i][j])/Npoints;
					metstuff.Met150Njet_genstat[i][j]=sqrt(1.0*metstuff.ana150Njet_bg[0][i][j])/Npoints;

					syst_tmp20[i][j]=0.0;
					syst_tmp25[i][j]=0.0;
					syst_tmp30[i][j]=0.0;
					syst_tmp35[i][j]=0.0;
					syst_tmp40[i][j]=0.0;
					syst_tmp45[i][j]=0.0;
					syst_tmp50[i][j]=0.0;
					syst_tmp75[i][j]=0.0;
					syst_tmp100[i][j]=0.0;
					syst_tmp150[i][j]=0.0;

					for(int l=0; l<10; l++)
					{
						if(l==0)
						{
							dummy=fabs(metstuff.ana20Njet_bg[0][i][j]-metstuff.ana20Njet_bg[1][i][j])*1.0;
							syst_tmp20[i][j]=syst_tmp20[i][j]+dummy*dummy;
							dummy=fabs(metstuff.ana25Njet_bg[0][i][j]-metstuff.ana25Njet_bg[1][i][j])*1.0;
							syst_tmp25[i][j]=syst_tmp25[i][j]+dummy*dummy;
							dummy=fabs(metstuff.ana30Njet_bg[0][i][j]-metstuff.ana30Njet_bg[1][i][j])*1.0;
							syst_tmp30[i][j]=syst_tmp30[i][j]+dummy*dummy;
							dummy=fabs(metstuff.ana35Njet_bg[0][i][j]-metstuff.ana35Njet_bg[1][i][j])*1.0;
							syst_tmp35[i][j]=syst_tmp35[i][j]+dummy*dummy;
							dummy=fabs(metstuff.ana40Njet_bg[0][i][j]-metstuff.ana40Njet_bg[1][i][j])*1.0;
							syst_tmp40[i][j]=syst_tmp40[i][j]+dummy*dummy;
							dummy=fabs(metstuff.ana45Njet_bg[0][i][j]-metstuff.ana45Njet_bg[1][i][j])*1.0;
							syst_tmp45[i][j]=syst_tmp45[i][j]+dummy*dummy;
							dummy=fabs(metstuff.ana50Njet_bg[0][i][j]-metstuff.ana50Njet_bg[1][i][j])*1.0;
							syst_tmp50[i][j]=syst_tmp50[i][j]+dummy*dummy;
							dummy=fabs(metstuff.ana75Njet_bg[0][i][j]-metstuff.ana75Njet_bg[1][i][j])*1.0;
							syst_tmp75[i][j]=syst_tmp75[i][j]+dummy*dummy;
							dummy=fabs(metstuff.ana100Njet_bg[0][i][j]-metstuff.ana100Njet_bg[1][i][j])*1.0;
							syst_tmp100[i][j]=syst_tmp100[i][j]+dummy*dummy;
							dummy=fabs(metstuff.ana150Njet_bg[0][i][j]-metstuff.ana150Njet_bg[1][i][j])*1.0;
							syst_tmp150[i][j]=syst_tmp150[i][j]+dummy*dummy;
						}
						else
						{
							dummy1=fabs(metstuff.ana20Njet_bg[0][i][j]-metstuff.ana20Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana20Njet_bg[0][i][j]-metstuff.ana20Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp20[i][j]=syst_tmp20[i][j]+dummy*dummy;
							dummy1=fabs(metstuff.ana25Njet_bg[0][i][j]-metstuff.ana25Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana25Njet_bg[0][i][j]-metstuff.ana25Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp25[i][j]=syst_tmp25[i][j]+dummy*dummy;
							dummy1=fabs(metstuff.ana30Njet_bg[0][i][j]-metstuff.ana30Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana30Njet_bg[0][i][j]-metstuff.ana30Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp30[i][j]=syst_tmp30[i][j]+dummy*dummy;
							dummy1=fabs(metstuff.ana35Njet_bg[0][i][j]-metstuff.ana35Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana35Njet_bg[0][i][j]-metstuff.ana35Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp35[i][j]=syst_tmp35[i][j]+dummy*dummy;
							dummy1=fabs(metstuff.ana40Njet_bg[0][i][j]-metstuff.ana40Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana40Njet_bg[0][i][j]-metstuff.ana40Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp40[i][j]=syst_tmp40[i][j]+dummy*dummy;
							dummy1=fabs(metstuff.ana45Njet_bg[0][i][j]-metstuff.ana45Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana45Njet_bg[0][i][j]-metstuff.ana45Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp45[i][j]=syst_tmp45[i][j]+dummy*dummy;
							dummy1=fabs(metstuff.ana50Njet_bg[0][i][j]-metstuff.ana50Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana50Njet_bg[0][i][j]-metstuff.ana50Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp50[i][j]=syst_tmp50[i][j]+dummy*dummy;
							dummy1=fabs(metstuff.ana75Njet_bg[0][i][j]-metstuff.ana75Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana75Njet_bg[0][i][j]-metstuff.ana75Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp75[i][j]=syst_tmp75[i][j]+dummy*dummy;
							dummy1=fabs(metstuff.ana100Njet_bg[0][i][j]-metstuff.ana100Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana100Njet_bg[0][i][j]-metstuff.ana100Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp100[i][j]=syst_tmp100[i][j]+dummy*dummy;
							dummy1=fabs(metstuff.ana150Njet_bg[0][i][j]-metstuff.ana150Njet_bg[l*2][i][j])*1.0;
							dummy2=fabs(metstuff.ana150Njet_bg[0][i][j]-metstuff.ana150Njet_bg[l*2+1][i][j])*1.0;
							dummy=(dummy1>dummy2) ? dummy1 : dummy2;
							syst_tmp150[i][j]=syst_tmp150[i][j]+dummy*dummy;
						}
					}
					metstuff.Met20Njet_gensyst[i][j]=sqrt(1.0*syst_tmp20[i][j])/Npoints;
					metstuff.Met25Njet_gensyst[i][j]=sqrt(1.0*syst_tmp25[i][j])/Npoints;
					metstuff.Met30Njet_gensyst[i][j]=sqrt(1.0*syst_tmp30[i][j])/Npoints;
					metstuff.Met35Njet_gensyst[i][j]=sqrt(1.0*syst_tmp35[i][j])/Npoints;
					metstuff.Met40Njet_gensyst[i][j]=sqrt(1.0*syst_tmp40[i][j])/Npoints;
					metstuff.Met45Njet_gensyst[i][j]=sqrt(1.0*syst_tmp45[i][j])/Npoints;
					metstuff.Met50Njet_gensyst[i][j]=sqrt(1.0*syst_tmp50[i][j])/Npoints;
					metstuff.Met75Njet_gensyst[i][j]=sqrt(1.0*syst_tmp75[i][j])/Npoints;
					metstuff.Met100Njet_gensyst[i][j]=sqrt(1.0*syst_tmp100[i][j])/Npoints;
					metstuff.Met150Njet_gensyst[i][j]=sqrt(1.0*syst_tmp150[i][j])/Npoints;
				}
				else
				{
					//--------- old events count
					metstuff.Met20Njet_genstat[i][j]=sqrt(1.0*metstuff.Met20Njet_def[i][j])/Npoints;
					metstuff.Met25Njet_genstat[i][j]=sqrt(1.0*metstuff.Met25Njet_def[i][j])/Npoints;
					metstuff.Met30Njet_genstat[i][j]=sqrt(1.0*metstuff.Met30Njet_def[i][j])/Npoints;
					metstuff.Met35Njet_genstat[i][j]=sqrt(1.0*metstuff.Met35Njet_def[i][j])/Npoints;
					metstuff.Met40Njet_genstat[i][j]=sqrt(1.0*metstuff.Met40Njet_def[i][j])/Npoints;
					metstuff.Met45Njet_genstat[i][j]=sqrt(1.0*metstuff.Met45Njet_def[i][j])/Npoints;
					metstuff.Met50Njet_genstat[i][j]=sqrt(1.0*metstuff.Met50Njet_def[i][j])/Npoints;
					metstuff.Met75Njet_genstat[i][j]=sqrt(1.0*metstuff.Met75Njet_def[i][j])/Npoints;
					metstuff.Met100Njet_genstat[i][j]=sqrt(1.0*metstuff.Met100Njet_def[i][j])/Npoints;
					metstuff.Met150Njet_genstat[i][j]=sqrt(1.0*metstuff.Met150Njet_def[i][j])/Npoints;

					metstuff.Met20Njet_gen[i][j]=1.0*metstuff.Met20Njet_def[i][j]/Npoints;
					metstuff.Met25Njet_gen[i][j]=1.0*metstuff.Met25Njet_def[i][j]/Npoints;
					metstuff.Met30Njet_gen[i][j]=1.0*metstuff.Met30Njet_def[i][j]/Npoints;
					metstuff.Met35Njet_gen[i][j]=1.0*metstuff.Met35Njet_def[i][j]/Npoints;
					metstuff.Met40Njet_gen[i][j]=1.0*metstuff.Met40Njet_def[i][j]/Npoints;
					metstuff.Met45Njet_gen[i][j]=1.0*metstuff.Met45Njet_def[i][j]/Npoints;
					metstuff.Met50Njet_gen[i][j]=1.0*metstuff.Met50Njet_def[i][j]/Npoints;
					metstuff.Met75Njet_gen[i][j]=1.0*metstuff.Met75Njet_def[i][j]/Npoints;
					metstuff.Met100Njet_gen[i][j]=1.0*metstuff.Met100Njet_def[i][j]/Npoints;
					metstuff.Met150Njet_gen[i][j]=1.0*metstuff.Met150Njet_def[i][j]/Npoints;
				}
			}
		}
	}

	ofstream outfile(fDatFileName, std::ios::out);
	if (! outfile)
	{
		std::cerr<<"Error in openning of .DAT file\n";
	}
	if(outfile)
	{
		outfile<<"---------------------------- Results ---------------------------"<<"\n";
		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 20 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met20Njet_dt[0][i]<<"            "
					<<metstuff.Met20Njet_gen[0][i]<<" +/- "
					<<metstuff.Met20Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met20Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met20Njet_dt[0][i]<<"            "
					<<metstuff.Met20Njet_gen[0][i]<<" +/- "
					<<metstuff.Met20Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met20Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met20Njet_dt[1][i]<<"            "
					<<metstuff.Met20Njet_gen[1][i]<<" +/- "
					<<metstuff.Met20Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met20Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met20Njet_dt[1][i]<<"            "
					<<metstuff.Met20Njet_gen[1][i]<<" +/- "
					<<metstuff.Met20Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met20Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met20Njet_dt[2][i]<<"            "
					<<metstuff.Met20Njet_gen[2][i]<<" +/- "
					<<metstuff.Met20Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met20Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met20Njet_dt[2][i]<<"            "
					<<metstuff.Met20Njet_gen[2][i]<<" +/- "
					<<metstuff.Met20Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met20Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 25 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met25Njet_dt[0][i]<<"            "
					<<metstuff.Met25Njet_gen[0][i]<<" +/- "
					<<metstuff.Met25Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met25Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met25Njet_dt[0][i]<<"            "
					<<metstuff.Met25Njet_gen[0][i]<<" +/- "
					<<metstuff.Met25Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met25Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met25Njet_dt[1][i]<<"            "
					<<metstuff.Met25Njet_gen[1][i]<<" +/- "
					<<metstuff.Met25Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met25Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met25Njet_dt[1][i]<<"            "
					<<metstuff.Met25Njet_gen[1][i]<<" +/- "
					<<metstuff.Met25Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met25Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met25Njet_dt[2][i]<<"            "
					<<metstuff.Met25Njet_gen[2][i]<<" +/- "
					<<metstuff.Met25Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met25Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met25Njet_dt[2][i]<<"            "
					<<metstuff.Met25Njet_gen[2][i]<<" +/- "
					<<metstuff.Met25Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met25Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 30 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met30Njet_dt[0][i]<<"            "
					<<metstuff.Met30Njet_gen[0][i]<<" +/- "
					<<metstuff.Met30Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met30Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met30Njet_dt[0][i]<<"            "
					<<metstuff.Met30Njet_gen[0][i]<<" +/- "
					<<metstuff.Met30Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met30Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met30Njet_dt[1][i]<<"            "
					<<metstuff.Met30Njet_gen[1][i]<<" +/- "
					<<metstuff.Met30Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met30Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met30Njet_dt[1][i]<<"            "
					<<metstuff.Met30Njet_gen[1][i]<<" +/- "
					<<metstuff.Met30Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met30Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met30Njet_dt[2][i]<<"            "
					<<metstuff.Met30Njet_gen[2][i]<<" +/- "
					<<metstuff.Met30Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met30Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met30Njet_dt[2][i]<<"            "
					<<metstuff.Met30Njet_gen[2][i]<<" +/- "
					<<metstuff.Met30Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met30Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 35 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met35Njet_dt[0][i]<<"            "
					<<metstuff.Met35Njet_gen[0][i]<<" +/- "
					<<metstuff.Met35Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met35Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met35Njet_dt[0][i]<<"            "
					<<metstuff.Met35Njet_gen[0][i]<<" +/- "
					<<metstuff.Met35Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met35Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met35Njet_dt[1][i]<<"            "
					<<metstuff.Met35Njet_gen[1][i]<<" +/- "
					<<metstuff.Met35Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met35Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met35Njet_dt[1][i]<<"            "
					<<metstuff.Met35Njet_gen[1][i]<<" +/- "
					<<metstuff.Met35Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met35Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met35Njet_dt[2][i]<<"            "
					<<metstuff.Met35Njet_gen[2][i]<<" +/- "
					<<metstuff.Met35Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met35Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met35Njet_dt[2][i]<<"            "
					<<metstuff.Met35Njet_gen[2][i]<<" +/- "
					<<metstuff.Met35Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met35Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 40 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met40Njet_dt[0][i]<<"            "
					<<metstuff.Met40Njet_gen[0][i]<<" +/- "
					<<metstuff.Met40Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met40Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met40Njet_dt[0][i]<<"            "
					<<metstuff.Met40Njet_gen[0][i]<<" +/- "
					<<metstuff.Met40Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met40Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met40Njet_dt[1][i]<<"            "
					<<metstuff.Met40Njet_gen[1][i]<<" +/- "
					<<metstuff.Met40Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met40Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met40Njet_dt[1][i]<<"            "
					<<metstuff.Met40Njet_gen[1][i]<<" +/- "
					<<metstuff.Met40Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met40Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met40Njet_dt[2][i]<<"            "
					<<metstuff.Met40Njet_gen[2][i]<<" +/- "
					<<metstuff.Met40Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met40Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met40Njet_dt[2][i]<<"            "
					<<metstuff.Met40Njet_gen[2][i]<<" +/- "
					<<metstuff.Met40Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met40Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 45 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met45Njet_dt[0][i]<<"            "
					<<metstuff.Met45Njet_gen[0][i]<<" +/- "
					<<metstuff.Met45Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met45Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met45Njet_dt[0][i]<<"            "
					<<metstuff.Met45Njet_gen[0][i]<<" +/- "
					<<metstuff.Met45Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met45Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met45Njet_dt[1][i]<<"            "
					<<metstuff.Met45Njet_gen[1][i]<<" +/- "
					<<metstuff.Met45Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met45Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met45Njet_dt[1][i]<<"            "
					<<metstuff.Met45Njet_gen[1][i]<<" +/- "
					<<metstuff.Met45Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met45Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met45Njet_dt[2][i]<<"            "
					<<metstuff.Met45Njet_gen[2][i]<<" +/- "
					<<metstuff.Met45Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met45Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met45Njet_dt[2][i]<<"            "
					<<metstuff.Met45Njet_gen[2][i]<<" +/- "
					<<metstuff.Met45Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met45Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 50 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met50Njet_dt[0][i]<<"            "
					<<metstuff.Met50Njet_gen[0][i]<<" +/- "
					<<metstuff.Met50Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met50Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met50Njet_dt[0][i]<<"            "
					<<metstuff.Met50Njet_gen[0][i]<<" +/- "
					<<metstuff.Met50Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met50Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met50Njet_dt[1][i]<<"            "
					<<metstuff.Met50Njet_gen[1][i]<<" +/- "
					<<metstuff.Met50Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met50Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met50Njet_dt[1][i]<<"            "
					<<metstuff.Met50Njet_gen[1][i]<<" +/- "
					<<metstuff.Met50Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met50Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met50Njet_dt[2][i]<<"            "
					<<metstuff.Met50Njet_gen[2][i]<<" +/- "
					<<metstuff.Met50Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met50Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met50Njet_dt[2][i]<<"            "
					<<metstuff.Met50Njet_gen[2][i]<<" +/- "
					<<metstuff.Met50Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met50Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 75 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met75Njet_dt[0][i]<<"            "
					<<metstuff.Met75Njet_gen[0][i]<<" +/- "
					<<metstuff.Met75Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met75Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met75Njet_dt[0][i]<<"            "
					<<metstuff.Met75Njet_gen[0][i]<<" +/- "
					<<metstuff.Met75Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met75Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met75Njet_dt[1][i]<<"            "
					<<metstuff.Met75Njet_gen[1][i]<<" +/- "
					<<metstuff.Met75Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met75Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met75Njet_dt[1][i]<<"            "
					<<metstuff.Met75Njet_gen[1][i]<<" +/- "
					<<metstuff.Met75Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met75Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met75Njet_dt[2][i]<<"            "
					<<metstuff.Met75Njet_gen[2][i]<<" +/- "
					<<metstuff.Met75Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met75Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met75Njet_dt[2][i]<<"            "
					<<metstuff.Met75Njet_gen[2][i]<<" +/- "
					<<metstuff.Met75Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met75Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 100 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met100Njet_dt[0][i]<<"            "
					<<metstuff.Met100Njet_gen[0][i]<<" +/- "
					<<metstuff.Met100Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met100Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met100Njet_dt[0][i]<<"            "
					<<metstuff.Met100Njet_gen[0][i]<<" +/- "
					<<metstuff.Met100Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met100Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met100Njet_dt[1][i]<<"            "
					<<metstuff.Met100Njet_gen[1][i]<<" +/- "
					<<metstuff.Met100Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met100Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met100Njet_dt[1][i]<<"            "
					<<metstuff.Met100Njet_gen[1][i]<<" +/- "
					<<metstuff.Met100Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met100Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met100Njet_dt[2][i]<<"            "
					<<metstuff.Met100Njet_gen[2][i]<<" +/- "
					<<metstuff.Met100Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met100Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met100Njet_dt[2][i]<<"            "
					<<metstuff.Met100Njet_gen[2][i]<<" +/- "
					<<metstuff.Met100Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met100Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"________________________________________________________________"<<"\n";
		outfile<<"                                                                "<<"\n";
		outfile<<" Met > 150 GeV         observed         predicted"<<"\n";
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet15="<<i<<"            "
				<<metstuff.Met150Njet_dt[0][i]<<"            "
					<<metstuff.Met150Njet_gen[0][i]<<" +/- "
					<<metstuff.Met150Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met150Njet_gensyst[0][i]<<"\n";
			else outfile<<" Njet15>="<<i<<"            "
				<<metstuff.Met150Njet_dt[0][i]<<"            "
					<<metstuff.Met150Njet_gen[0][i]<<" +/- "
					<<metstuff.Met150Njet_genstat[0][i]<<" +/- "
					<<metstuff.Met150Njet_gensyst[0][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet20="<<i<<"            "
				<<metstuff.Met150Njet_dt[1][i]<<"            "
					<<metstuff.Met150Njet_gen[1][i]<<" +/- "
					<<metstuff.Met150Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met150Njet_gensyst[1][i]<<"\n";
			else outfile<<" Njet20>="<<i<<"            "
				<<metstuff.Met150Njet_dt[1][i]<<"            "
					<<metstuff.Met150Njet_gen[1][i]<<" +/- "
					<<metstuff.Met150Njet_genstat[1][i]<<" +/- "
					<<metstuff.Met150Njet_gensyst[1][i]<<std::endl;
		}
		outfile<<"                                                                "<<"\n";
		for(int i=0; i<10; i++)
		{
			if(i<9) outfile<<"  Njet25="<<i<<"            "
				<<metstuff.Met150Njet_dt[2][i]<<"            "
					<<metstuff.Met150Njet_gen[2][i]<<" +/- "
					<<metstuff.Met150Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met150Njet_gensyst[2][i]<<"\n";
			else outfile<<" Njet25>="<<i<<"            "
				<<metstuff.Met150Njet_dt[2][i]<<"            "
					<<metstuff.Met150Njet_gen[2][i]<<" +/- "
					<<metstuff.Met150Njet_genstat[2][i]<<" +/- "
					<<metstuff.Met150Njet_gensyst[2][i]<<std::endl;
		}

		outfile<<"---- The End!------"<<std::endl;
	}
	return;
}

//_______________ selects large MET events; writes them to file
int JetFilterModuleV2::MyLargeMetEvent(JetStuff jetstuff, CommonStuff miscstuff, int metscenario, int Nvx12, int runN, int eventN) {
	int largemet_code=0;
	if(fSelectMetEvent!=0)
	{
		int Nexo_obj=0;
		//__________________cuts
		//   double cut_MEt=30.0;
		double cut_MEt=20.0;

		//_________________ contribution from MET
		if(metscenario==0 && miscstuff.myMET_raw.Mod() > cut_MEt) Nexo_obj++;
		if(metscenario==1 && jetstuff.myMETcorr_th5.Mod() > cut_MEt) Nexo_obj++;
		if(metscenario==2 && jetstuff.myMETcorr_th10.Mod() > cut_MEt) Nexo_obj++;
		if(metscenario==3 && jetstuff.myMETcorr_th15.Mod() > cut_MEt) Nexo_obj++;
		if(metscenario==4 && jetstuff.myMETcorr_th20.Mod() > cut_MEt) Nexo_obj++;
		if((metscenario<0 || metscenario>4) && miscstuff.myMET_raw.Mod() > cut_MEt) Nexo_obj++;

		//__________________________ printing out interesting event
		if(Nexo_obj>0)
		{
			largemet_code=1;
			//____________________ opening existing file
			ofstream outfile(fLargeMetEventFileName, std::ios::app);
			if (! outfile)
			{
				std::cerr<<"Error in openning of .DAT file\n";
			}
			if(outfile)
			{
				outfile<<"_____________________________________________________________________________"<<"\n";
				outfile<<"Run "<<runN<<"  Event "<<eventN<<"\n";
				outfile<<"Nvx12 "<<Nvx12<<"\n";
				outfile<<"_______ Npho="<<miscstuff.myCorrPhoton.size()<<"\n";
				//_________________ printing photons
				for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
				{
					outfile<<"...... pho-"<<i+1
						<<" Et="<<miscstuff.myCorrPhoton.at(i).Pt()
						<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrPhoton.at(i).Phi())
						<<" eta="<<miscstuff.myCorrPhoton.at(i).Eta()
						<<" eta_det="<<miscstuff.myPhoEtaDet.at(i)<<"\n";
				}
				outfile<<"_______ Nele="<<miscstuff.myCorrElectron.size()<<"\n";
				//_________________ printing electrons
				for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
				{
					outfile<<"...... ele-"<<i+1
						<<" Et="<<miscstuff.myCorrElectron.at(i).Pt()
						<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrElectron.at(i).Phi())
						<<" eta="<<miscstuff.myCorrElectron.at(i).Eta()
						<<" eta_det="<<miscstuff.myEleEtaDet.at(i)<<"\n";
				}
				// 	  outfile<<"_______ Njet15="<<jetstuff.myNjet_th15<<"\n";
				outfile<<"_______ Njet5="<<jetstuff.myNjet_th5<<"\n";
				//_________________ printing jets
				for(int i=0; i<jetstuff.myNjet; i++)
				{
					// 	      if(jetstuff.Jet_lev6_noEMobj.at(i).Pt() >=15.0) outfile<<"...... jet-"<<i+1
					if(jetstuff.Jet_lev6_noEMobj.at(i).Pt() >=5.0) outfile<<"...... jet-"<<i+1
						<<" Et="<<jetstuff.Jet_lev6_noEMobj.at(i).Pt()
							<<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj.at(i).Phi())
							<<" eta="<<jetstuff.Jet_lev6_noEMobj.at(i).Eta()
							<<" eta_det="<<jetstuff.EtaDet.at(i)
							<<" raw_em_fr="<<jetstuff.EmFrRaw.at(i)<<"\n";

				}
				//_________________ printing MET
				if(metscenario==0) outfile<<"_______ MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";
				if(metscenario==1) outfile<<"_______ MET="<<jetstuff.myMETcorr_th5.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th5.Phi())<<"\n";
				if(metscenario==2) outfile<<"_______ MET="<<jetstuff.myMETcorr_th10.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th10.Phi())<<"\n";
				if(metscenario==3) outfile<<"_______ MET="<<jetstuff.myMETcorr_th15.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th15.Phi())<<"\n";
				if(metscenario==4) outfile<<"_______ MET="<<jetstuff.myMETcorr_th20.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th20.Phi())<<"\n";
				if(metscenario<0 || metscenario>4) outfile<<"_______ MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";
				outfile<<"............................................................................."<<std::endl;
			}
		}
	}
	return largemet_code;
}

//_______________ event print-out
void JetFilterModuleV2::MyDumpEvent(JetStuff jetstuff, CommonStuff miscstuff, int Nvx12, int runN, int eventN) {

	//__________________________ printing out interesting event
	if(fDumpEvent==1)
	{
		//____________________ opening existing file
		ofstream outfile(fDumpEventFileName, std::ios::app);
		if (! outfile)
		{
			std::cerr<<"Error in openning of .DAT file\n";
		}
		if(outfile)
		{
			outfile<<"_____________________________________________________________________________"<<"\n";
			outfile<<"Run "<<runN<<"  Event "<<eventN<<"\n";
			outfile<<"Nvx12 "<<Nvx12<<"\n";
			outfile<<"_______ Npho="<<miscstuff.myCorrPhoton.size()<<"\n";
			//_________________ printing photons
			outfile<<"___ printing corrected photons ____"<<"\n";
			for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
			{
				outfile<<"...... pho-"<<i+1
					<<" Et="<<miscstuff.myCorrPhoton.at(i).Pt()
					<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrPhoton.at(i).Phi())
					<<" eta="<<miscstuff.myCorrPhoton.at(i).Eta()
					<<" eta_det="<<miscstuff.myPhoEtaDet.at(i)<<"\n";
			}
			outfile<<"___ printing raw photons ____"<<"\n";
			for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
			{
				outfile<<"...... pho-"<<i+1
					<<" raw Et="<<miscstuff.myRawPhoton.at(i).Pt()
					<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myRawPhoton.at(i).Phi())
					<<" eta="<<miscstuff.myRawPhoton.at(i).Eta()
					<<" index="<<miscstuff.myPhoInd.at(i)<<"\n";
			}

			outfile<<"_______ Nele="<<miscstuff.myCorrElectron.size()<<"\n";
			//_________________ printing electrons
			outfile<<"___ printing corrected electrons ____"<<"\n";
			for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
			{
				outfile<<"...... ele-"<<i+1
					<<" Et="<<miscstuff.myCorrElectron.at(i).Pt()
					<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrElectron.at(i).Phi())
					<<" eta="<<miscstuff.myCorrElectron.at(i).Eta()
					<<" eta_det="<<miscstuff.myEleEtaDet.at(i)<<"\n";
			}
			outfile<<"___ printing raw electrons ____"<<"\n";
			for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
			{
				outfile<<"...... ele-"<<i+1
					<<" raw Et="<<miscstuff.myRawElectron.at(i).Pt()
					<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myRawElectron.at(i).Phi())
					<<" eta="<<miscstuff.myRawElectron.at(i).Eta()
					<<" index="<<miscstuff.myEleInd.at(i)<<"\n";
			}

			outfile<<"_______ Njet="<<jetstuff.myNjet<<"\n";
			//_________________ printing jets
			outfile<<"___ printing corrected jets after removing EM object ____"<<"\n";
			for(int i=0; i<jetstuff.myNjet; i++)
			{
				outfile<<"...... jet-"<<i+1
					<<" Et="<<jetstuff.Jet_lev6_noEMobj.at(i).Pt()
					<<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj.at(i).Phi())
					<<" eta="<<jetstuff.Jet_lev6_noEMobj.at(i).Eta()
					<<" cor_eta_det="<<jetstuff.EtaDetCorr.at(i)
					<<" cor_em_fr="<<jetstuff.EmFrCorr.at(i)<<"\n";
			}
			outfile<<"___ printing raw jets after removing EM object ____"<<"\n";
			for(int i=0; i<jetstuff.myNjet; i++)
			{
				outfile<<"...... jet-"<<i+1
					<<" Et="<<jetstuff.Jet_raw_noEMobj.at(i).Pt()
					<<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_raw_noEMobj.at(i).Phi())
					<<" eta="<<jetstuff.Jet_raw_noEMobj.at(i).Eta()
					<<" cor_eta_det="<<jetstuff.EtaDetCorr.at(i)
					<<" cor_em_fr="<<jetstuff.EmFrCorr.at(i)<<"\n";
			}

			outfile<<"___ printing corrected jets before removing EM object ____"<<"\n";
			for(int i=0; i<jetstuff.myNjet; i++)
			{
				outfile<<"...... jet-"<<i+1
					<<" Et="<<jetstuff.Jet_lev6.at(i).Pt()
					<<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_lev6.at(i).Phi())
					<<" eta="<<jetstuff.Jet_lev6.at(i).Eta()
					<<" raw_eta_det="<<jetstuff.EtaDet.at(i)
					<<" raw_em_fr="<<jetstuff.EmFrRaw.at(i)<<"\n";
			}
			outfile<<"___ printing raw jets before removing EM object ____"<<"\n";
			for(int i=0; i<jetstuff.myNjet; i++)
			{
				outfile<<"...... jet-"<<i+1
					<<" Et="<<jetstuff.Jet_raw.at(i).Pt()
					<<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_raw.at(i).Phi())
					<<" eta="<<jetstuff.Jet_raw.at(i).Eta()
					<<" raw_eta_det="<<jetstuff.EtaDet.at(i)
					<<" raw_em_fr="<<jetstuff.EmFrRaw.at(i)<<"\n";
			}

			//_________________ printing SumEt
			outfile<<"_______ raw SumEt="<<miscstuff.mySumEt_raw<<"\n";
			outfile<<"_______ corr SumEt="<<jetstuff.mySumEtCorr_th15<<"\n";
			outfile<<"............................................................................."<<std::endl;

			//_________________ printing MET
			outfile<<"_______ raw MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";
			outfile<<"_______ corr MET="<<jetstuff.myMETcorr_th15.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th15.Phi())<<"\n";
			outfile<<"............................................................................."<<std::endl;
		}
	}
	return;
}

//_______________ selects my exotic events; writes them to file
int JetFilterModuleV2::MyExoticEvent(JetStuff jetstuff, CommonStuff miscstuff, int metscenario, int runN, int eventN) {
	int exo_code=0;
	if(fSelectExoticEvent!=0)
	{
		int Nexo_obj=0;
		//__________________cuts
		int cut_Nexo_obj=3;
		double cut_Et_pho[3]={20.0,20.0,13.0};
		double cut_Et_ele=13.0;
		double cut_Et_jet=50.0; // default = 15.0
		double cut_MEt=30.0;

		//_________________ counting jets
		for(int i=0; i<jetstuff.myNjet; i++)
		{
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt() > cut_Et_jet) Nexo_obj++;
		}
		//_________________ counting photons
		for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
		{
			if(i<3 && miscstuff.myCorrPhoton.at(i).Pt() > cut_Et_pho[i]) Nexo_obj++;
			if(i>2 && miscstuff.myCorrPhoton.at(i).Pt() > cut_Et_pho[2]) Nexo_obj++;
		}
		//_________________ counting electrons
		for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
		{
			if(miscstuff.myCorrElectron.at(i).Pt() > cut_Et_ele) Nexo_obj++;
		}
		//_________________ contribution from MET
		if(metscenario==0 && miscstuff.myMET_raw.Mod() > cut_MEt) Nexo_obj++;
		if(metscenario==1 && jetstuff.myMETcorr_th5.Mod() > cut_MEt) Nexo_obj++;
		if(metscenario==2 && jetstuff.myMETcorr_th10.Mod() > cut_MEt) Nexo_obj++;
		if(metscenario==3 && jetstuff.myMETcorr_th15.Mod() > cut_MEt) Nexo_obj++;
		if(metscenario==4 && jetstuff.myMETcorr_th20.Mod() > cut_MEt) Nexo_obj++;
		if((metscenario<0 || metscenario>4) && miscstuff.myMET_raw.Mod() > cut_MEt) Nexo_obj++;

		//__________________________ printing out interesting event
		if(Nexo_obj>=cut_Nexo_obj)
		{
			exo_code=1;
			//____________________ opening existing file
			ofstream outfile(fRearEventFileName, std::ios::app);
			if (! outfile)
			{
				std::cerr<<"Error in openning of .DAT file\n";
			}
			if(outfile)
			{
				outfile<<"_____________________________________________________________________________"<<"\n";
				outfile<<"Run "<<runN<<"  Event "<<eventN<<"\n";
				outfile<<"_______ Npho="<<miscstuff.myCorrPhoton.size()<<"\n";
				//_________________ printing photons
				for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
				{
					outfile<<"...... pho-"<<i+1
						<<" Et="<<miscstuff.myCorrPhoton.at(i).Pt()
						<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrPhoton.at(i).Phi())
						<<" eta="<<miscstuff.myCorrPhoton.at(i).Eta()<<"\n";
				}
				outfile<<"_______ Nele="<<miscstuff.myCorrElectron.size()<<"\n";
				//_________________ printing electrons
				for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
				{
					outfile<<"...... ele-"<<i+1
						<<" Et="<<miscstuff.myCorrElectron.at(i).Pt()
						<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrElectron.at(i).Phi())
						<<" eta="<<miscstuff.myCorrElectron.at(i).Eta()<<"\n";
				}
				outfile<<"_______ Njet15="<<jetstuff.myNjet_th15<<"\n";
				//_________________ printing jets
				for(int i=0; i<jetstuff.myNjet; i++)
				{
					if(jetstuff.Jet_lev6_noEMobj.at(i).Pt() >=15.0) outfile<<"...... jet-"<<i+1
						<<" Et="<<jetstuff.Jet_lev6_noEMobj.at(i).Pt()
							<<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj.at(i).Phi())
							<<" eta="<<jetstuff.Jet_lev6_noEMobj.at(i).Eta()<<"\n";
				}
				//_________________ printing MET
				if(metscenario==0) outfile<<"_______ MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";
				if(metscenario==1) outfile<<"_______ MET="<<jetstuff.myMETcorr_th5.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th5.Phi())<<"\n";
				if(metscenario==2) outfile<<"_______ MET="<<jetstuff.myMETcorr_th10.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th10.Phi())<<"\n";
				if(metscenario==3) outfile<<"_______ MET="<<jetstuff.myMETcorr_th15.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th15.Phi())<<"\n";
				if(metscenario==4) outfile<<"_______ MET="<<jetstuff.myMETcorr_th20.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th20.Phi())<<"\n";
				if(metscenario<0 || metscenario>4) outfile<<"_______ MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";
				outfile<<"............................................................................."<<std::endl;
			}
		}
	}
	return exo_code;
}


//_______________ selects my very exotic events; writes them to file
int JetFilterModuleV2::MyVeryExoticEvent(JetStuff jetstuff, CommonStuff miscstuff, int metscenario,
		int runN, int eventN, int nvx12, std::vector<double> zvx_vec) {
	int exo_code=0;
	if(fSelectExoticEvent!=0)
	{

		//__________________cuts

		double cut_Et_pho=200.0;
		double cut_Et_jet=150.0;
		double cut_MEt=40.0;
		double cut_Ht=450.0;
		double cut_Mall=450.0;
		double cut_Mgg=250.0;
		double cut_Mggj=400.0;
		double cut_Mggjj=500.0;
		double cut_Mgj=300.0;
		double cut_Mjj=200.0;

		double _Ht=0.0;
		double _Mall=0.0;
		double _Mgg=0.0;
		double _Mggj1=0.0;
		double _Mggjj=0.0;
		double _Mjj=0.0;
		double _Mg1j1=0.0;
		double _Mg2j2=0.0;
		double _Mg1j2=0.0;
		double _Mg2j1=0.0;

		TLorentzVector all_vec(0.0,0.0,0.0,0.0);
		//_________________ counting jets
		for(int i=0; i<jetstuff.myNjet; i++)
		{
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt() > cut_Et_jet) exo_code++;
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt() > 15.0)
			{
				_Ht=_Ht+jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				all_vec=all_vec+jetstuff.Jet_lev6_noEMobj.at(i);
			}
		}
		//_________________ counting photons
		for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
		{
			if(miscstuff.myCorrPhoton.at(i).Pt() > cut_Et_pho) exo_code++;
			_Ht=_Ht+miscstuff.myCorrPhoton.at(i).Pt();
			all_vec=all_vec+miscstuff.myCorrPhoton.at(i);
		}
		//_________________ contribution from MET

		TVector2 met_vec(0.0,0.0);
		if(metscenario==0) met_vec=miscstuff.myMET_raw;
		if(metscenario==1) met_vec=jetstuff.myMETcorr_th5;
		if(metscenario==2) met_vec=jetstuff.myMETcorr_th10;
		if(metscenario==3) met_vec=jetstuff.myMETcorr_th15;
		if(metscenario==4) met_vec=jetstuff.myMETcorr_th20;
		if(metscenario<0 || metscenario>4) met_vec=miscstuff.myMET_raw;
		if(met_vec.Mod() > cut_MEt) exo_code++;
		_Ht=_Ht+met_vec.Mod();

		//_________________ counting electrons
		for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
		{
			_Ht=_Ht+miscstuff.myCorrElectron.at(i).Pt();
			all_vec=all_vec+miscstuff.myCorrElectron.at(i);
		}

		_Mall=all_vec.M();
		_Mgg=(miscstuff.myCorrPhoton[0]+miscstuff.myCorrPhoton[1]).M();
		if(jetstuff.myNjet_th15>0) _Mggj1=(miscstuff.myCorrPhoton[0]+miscstuff.myCorrPhoton[1]+jetstuff.Jet_lev6_noEMobj[0]).M();
		if(jetstuff.myNjet_th15>1)
		{
			_Mggjj=(miscstuff.myCorrPhoton[0]+miscstuff.myCorrPhoton[1]+jetstuff.Jet_lev6_noEMobj[0]+jetstuff.Jet_lev6_noEMobj[1]).M();
			_Mjj=(jetstuff.Jet_lev6_noEMobj[0]+jetstuff.Jet_lev6_noEMobj[1]).M();
			_Mg1j1=(miscstuff.myCorrPhoton[0]+jetstuff.Jet_lev6_noEMobj[0]).M();
			_Mg2j2=(miscstuff.myCorrPhoton[1]+jetstuff.Jet_lev6_noEMobj[1]).M();
			_Mg1j2=(miscstuff.myCorrPhoton[0]+jetstuff.Jet_lev6_noEMobj[1]).M();
			_Mg2j1=(miscstuff.myCorrPhoton[1]+jetstuff.Jet_lev6_noEMobj[0]).M();
		}

		if(_Ht > cut_Ht) exo_code++;
		if(_Mall > cut_Mall) exo_code++;
		if(_Mgg > cut_Mgg) exo_code++;
		if(_Mggj1 > cut_Mggj) exo_code++;
		if(_Mggjj > cut_Mggjj) exo_code++;
		if(_Mjj > cut_Mjj) exo_code++;
		if(_Mg1j1 > cut_Mgj) exo_code++;
		if(_Mg2j2 > cut_Mgj) exo_code++;
		if(_Mg1j2 > cut_Mgj) exo_code++;
		if(_Mg2j1 > cut_Mgj) exo_code++;


		//__________________________ printing out interesting event
		if(exo_code>0)
		{
			//____________________ opening existing file
			ofstream outfile(fExoticEventFileName, std::ios::app);
			if (! outfile)
			{
				std::cerr<<"Error in openning of .DAT file\n";
			}
			if(outfile)
			{
				outfile<<"_____________________________________________________________________________"<<"\n";
				outfile<<"Run "<<runN<<"  Event "<<eventN<<"\n";
				outfile<<"_______ Nvx12="<<nvx12<<"\n";
				for(unsigned int i=0; i<zvx_vec.size(); i++)
				{
					outfile<<"...... vx-"<<i<<" z="<<zvx_vec[i]<<"\n";
				}
				//_________________ printing photons
				for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
				{
					outfile<<"...... pho-"<<i+1
						<<" Et="<<miscstuff.myCorrPhoton.at(i).Pt()
						<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrPhoton.at(i).Phi())
						<<" eta="<<miscstuff.myCorrPhoton.at(i).Eta()
						<<" eta_det="<<miscstuff.myPhoEtaDet.at(i)<<"\n";
				}
				outfile<<"_______ Nele="<<miscstuff.myCorrElectron.size()<<"\n";
				//_________________ printing electrons
				for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
				{
					outfile<<"...... ele-"<<i+1
						<<" Et="<<miscstuff.myCorrElectron.at(i).Pt()
						<<" Phi="<<TVector2::Phi_0_2pi(miscstuff.myCorrElectron.at(i).Phi())
						<<" eta="<<miscstuff.myCorrElectron.at(i).Eta()<<"\n";
				}
				outfile<<"_______ Njet15="<<jetstuff.myNjet_th15<<"\n";
				//_________________ printing jets
				for(int i=0; i<jetstuff.myNjet; i++)
				{
					if(jetstuff.Jet_lev6_noEMobj.at(i).Pt() >=15.0) outfile<<"...... jet-"<<i+1
						<<" Et="<<jetstuff.Jet_lev6_noEMobj.at(i).Pt()
							<<" Phi="<<TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj.at(i).Phi())
							<<" eta="<<jetstuff.Jet_lev6_noEMobj.at(i).Eta()
							<<" eta_det="<<jetstuff.EtaDetCorr.at(i)
							<<" emFr="<<jetstuff.EmFrCorr.at(i)<<"\n";
				}
				//_________________ printing MET
				if(metscenario==0) outfile<<"_______ MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";
				if(metscenario==1) outfile<<"_______ MET="<<jetstuff.myMETcorr_th5.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th5.Phi())<<"\n";
				if(metscenario==2) outfile<<"_______ MET="<<jetstuff.myMETcorr_th10.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th10.Phi())<<"\n";
				if(metscenario==3) outfile<<"_______ MET="<<jetstuff.myMETcorr_th15.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th15.Phi())<<"\n";
				if(metscenario==4) outfile<<"_______ MET="<<jetstuff.myMETcorr_th20.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(jetstuff.myMETcorr_th20.Phi())<<"\n";
				if(metscenario<0 || metscenario>4) outfile<<"_______ MET="<<miscstuff.myMET_raw.Mod()<<" MET_phi="<<TVector2::Phi_0_2pi(miscstuff.myMET_raw.Phi())<<"\n";

				outfile<<"_______ Ht="<<_Ht<<"\n";
				outfile<<"_______ M(all)="<<_Mall<<"\n";
				outfile<<"_______ M(dipho)="<<_Mgg<<"\n";
				outfile<<"_______ M(dipho+jet)="<<_Mggj1<<"\n";
				outfile<<"_______ M(dipho+dijet)="<<_Mggjj<<"\n";
				outfile<<"_______ M(dijet)="<<_Mjj<<"\n";
				outfile<<"_______ M(pho1+jet1)="<<_Mg1j1<<"\n";
				outfile<<"_______ M(pho2+jet2)="<<_Mg2j2<<"\n";
				outfile<<"_______ M(pho1+jet2)="<<_Mg1j2<<"\n";
				outfile<<"_______ M(pho2+jet1)="<<_Mg2j1<<"\n";

				outfile<<"............................................................................."<<std::endl;
			}
		}
	}
	return exo_code;
}


//________ creates a list of towers for all jets
void JetFilterModuleV2::MatchCalorTowers(int jet_ind,
		TStnJetBlock* fJetBlock,
		TCalDataBlock *fCalDataBlock,
		CalDataArray* cdholder) {
	cdholder->clear();
	TStnLinkBlock* links = fJetBlock->TowerLinkList();
	int nptow = links->NLinks(jet_ind);

	TCalTower* ctower = NULL;
	for(int j=0; j<nptow; j++)
	{
		int iptow = links->Index(jet_ind,j);
		int iphi = iptow & 0x3F;
		int ieta = (iptow>>8) & 0x3F;

		ctower = fCalDataBlock->Tower(ieta,iphi);
		cdholder->push_back(ctower);
	}

	std::sort(cdholder->begin(), cdholder->end(),SortTowersByEnergy);

	return;
}

//________ creates a list of towers for all EM objects
void JetFilterModuleV2::MatchCalorTowers(int pho_ind,
		TStnPhotonBlock* fPhotonBlock,
		TCalDataBlock *fCalDataBlock,
		CalDataArray* cdholder) {
	cdholder->clear();
	TStnLinkBlock* links = fPhotonBlock->TowerLinkList();
	int nptow = links->NLinks(pho_ind);
	TCalTower* ctower = NULL;
	for(int j=0; j<nptow; j++)
	{
		int iptow = links->Index(pho_ind,j);
		int iphi = iptow & 0x3F;
		int ieta = (iptow>>8) & 0x3F;

		ctower = fCalDataBlock->Tower(ieta,iphi);
		cdholder->push_back(ctower);
	}

	std::sort(cdholder->begin(),cdholder->end(),SortTowersByEnergy);

	return;
}

//________ creates a list of towers for electrons
void JetFilterModuleV2::MatchCalorTowers(int ele_ind,
		TStnElectronBlock* fElectronBlock,
		TCalDataBlock *fCalDataBlock,
		CalDataArray* cdholder) {
	cdholder->clear();
	TStnLinkBlock* links = fElectronBlock->TowerLinkList();
	int nptow = links->NLinks(ele_ind);
	TCalTower* ctower = NULL;
	for(int j=0; j<nptow; j++)
	{
		int iptow = links->Index(ele_ind,j);
		int iphi = iptow & 0x3F;
		int ieta = (iptow>>8) & 0x3F;

		ctower = fCalDataBlock->Tower(ieta,iphi);
		cdholder->push_back(ctower);
	}

	std::sort(cdholder->begin(),cdholder->end(),SortTowersByEnergy);

	return;
}

//_____________________________________________________________________________
int JetFilterModuleV2::BeginJob() {


	//_____________________________________________________ register the data block
	RegisterDataBlock("JetBlock","TStnJetBlock",&fJetBlockClu04);
	RegisterDataBlock("PROD@JetCluModule-had-cone0.4","TStnJetBlock",&fHadJetBlockClu04);

	//---------- need this for jet-EM matching
	RegisterDataBlock("CalDataBlock","TCalDataBlock",&fCalData);
	RegisterDataBlock("PhotonBlock","TStnPhotonBlock",&fPhotonBlock);
	RegisterDataBlock("ElectronBlock","TStnElectronBlock",&fElectronBlock);
	RegisterDataBlock("MetBlock","TStnMetBlock",&fMetBlock);		//added to replace Sasha's EventFiltermod 03-14-2008
	fHeaderBlock 	= (TStnHeaderBlock*)   RegisterDataBlock("HeaderBlock","TStnHeaderBlock");
	fGenpBlock     = (TGenpBlock*)        RegisterDataBlock("GenpBlock","TGenpBlock");

	//_____________________________________________________ book histograms

	BookHistograms();
	CleanMetResults(met_results);

	//_________ opening a file to write my exotic events (to be filled later)
	ofstream outfile(fRearEventFileName, std::ios::out); // re-creating file
	if (! outfile)
	{
		std::cerr<<"Error in openning of .DAT file\n";
	}
	if(outfile)
	{
		outfile<<"--------------- My exotic event list ---------------------------"<<"\n";
		outfile<<"________________________________________________________________"<<std::endl;
	}
	ofstream outfile1(fLargeMetEventFileName, std::ios::out); // re-creating file
	if (! outfile1)
	{
		std::cerr<<"Error in openning of .DAT file\n";
	}
	ofstream outfile0(fDumpEventFileName, std::ios::out); // re-creating file
	if (! outfile0)
	{
		std::cerr<<"Error in openning of .DAT file\n";
	}
	if(outfile1)
	{
		outfile1<<"--------------- My large MET event list ------------------------"<<"\n";
		outfile1<<"________________________________________________________________"<<std::endl;
	}
	ofstream outfile2(fExoticEventFileName, std::ios::out); // re-creating file
	if (! outfile2)
	{
		std::cerr<<"Error in openning of .DAT file\n";
	}
	if(outfile2)
	{
		outfile2<<"--------------- My very exotic event list ------------------------"<<"\n";
		outfile2<<"________________________________________________________________"<<std::endl;
	}

	EventCount_b=0;
	EventCount_a=0;



	// sam's stuff -- 03-14-2008

	mycounter.evtsRunOver = 0;
	mycounter.evtsPassModule = 0;
	mycounter.evtsMoreHad = 0;

	initSpMod = (InitSuperPhotons*) ((TStnAna*) GetAna()->GetModule("InitSuperPhotons"));
	if (initSpMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module InitSuperPhotons!");
		bRunPermit = false;
	}
	tightMod = (TagTightPhotons*) ((TStnAna*) GetAna()->GetModule("TagTightPhotons"));
	if (tightMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TagTightPhotons!");
		bRunPermit = false;
	}
	looseMod = (TagLoosePhotons*) ((TStnAna*) GetAna()->GetModule("TagLoosePhotons"));
	if (looseMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TagLoosePhotons!");
		bRunPermit = false;
	}
	tightEleMod = (TagElectrons*) ((TStnAna*) GetAna()->GetModule("TagElectrons"));
	if (tightEleMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TagElectrons!");
		bRunPermit = false;
	}
	looseEleMod = (TagLooseElectrons*) ((TStnAna*) GetAna()->GetModule("TagLooseElectrons"));
	if (looseEleMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TagElectrons!");
		bRunPermit = false;
	}
	trigMod = (TriggerModule*) ((TStnAna*) GetAna()->GetModule("TriggerModule"));
	if (trigMod == NULL) {
		StdErr(__FILE__,__LINE__,3,"Could not find module TriggerModule!");
		bRunPermit = false;
	}


	//removed these from Event loop
	if (!fPhotonBlock) {
		StdOut(__FILE__,__LINE__,3,"PhotonBlock not found!");
		bRunPermit = false;
		return 0;
	}
	if (!fElectronBlock) {
		StdOut(__FILE__,__LINE__,3,"ElectronBlock not found!");
		bRunPermit = false;
		return 0;
	}

	if (!fCalData) {
		StdOut(__FILE__,__LINE__,3,"CalDataBlock Block not found!");
		bRunPermit = false;
		return 0;
	}
	if (!fMetBlock)	{
		StdOut(__FILE__,__LINE__,3,"MetBlock not found!");
		bRunPermit = false;
		return 0;
	}
	if (!fJetBlockClu04)	{
		StdOut(__FILE__,__LINE__,3,"JetBlock not found!");
		bRunPermit = false;
		return 0;
	}

	vJerParam.clear();
	vJerParamErr.clear();
	vJerSystParam.clear();
	vJerSystParamErr.clear();

	for (int iEtaBin = 0; iEtaBin<15; ++iEtaBin)
	{
		std::vector<double> rowOfPar, rowOfParErr;
		for (int iPar =0; iPar<13; ++iPar)
		{
			rowOfPar.push_back(JERparam(0,iEtaBin, iPar));
			rowOfParErr.push_back(JERparam(1,iEtaBin, iPar));
		}
		vJerParam.push_back(rowOfPar);
		vJerParamErr.push_back(rowOfParErr);

		std::vector<double> rowOfSystPar, rowOfSystParErr;
		for (int iPar =0; iPar<31; ++iPar)
		{
			rowOfSystPar.push_back(new3JERparam(0,iEtaBin, iPar));
			rowOfSystParErr.push_back(new3JERparam(1,iEtaBin, iPar));
		}
		vJerSystParam.push_back(rowOfSystPar);
		vJerSystParamErr.push_back(rowOfSystParErr);
	}

	/*std::cout << "JER VECTOR: " << std::endl;
	for (int iEtaBin = 0; iEtaBin<15; ++iEtaBin)
	{
		for (int iPar =0; iPar<13; ++iPar)
		{
			double val = vJerParam[iEtaBin][iPar];
			//std::cout << "[" << iEtaBin << "]" << "[" << val  << "]" <<  std::endl;
		}
	}
	*/
	//std::cout << "JER syst VECTOR: " << std::endl;
	/*for (int iEtaBin = 0; iEtaBin<15; ++iEtaBin)
	{
		for (int iPar =0; iPar<31; ++iPar)
		{
			double val = vJerSystParam[iEtaBin][iPar];
			//std::cout << "[" << iEtaBin << ", " << iPar << "]" << "[" << val  << "]" <<  std::endl;
		}
	}
	*/
	
	// there are exact bin values for the E<60GeV
	// for old JER fits.
	// I need a 3D vector [eta-15][par-5][ebin-60] 
	// Fit histogram was binned in 1GeV hence 60 bin values
	// 
	// i cannot figure out why i get 'undefined symbol: _ZN17JetFilterModuleV212vJerBinValueE'
	// this is no different from other multi-d vector definitions!
	/*
	vJerBinValue.clear();
	for (int iEtaBin = 0; iEtaBin<15; ++iEtaBin)
	{
		std::vector<std::vector<double> > v2D;
		for (int iPar =0; iPar<5; ++iPar)
		{
			
			std::vector<double> rowOfBinValues;
			for (int iBin =0; iBin<60; ++iBin)
			{
				rowOfBinValues.push_back(oldJERBinValues(iEtaBin, iPar, iBin));
			}
			
			v2D.push_back(rowOfBinValues);	//[parameter][binvalue]
		}
		vJerBinValue.push_back(v2D);
	}
	*/

	return 0;
}


//_____________________________________________________________________________
int JetFilterModuleV2::BeginRun() {

	fJTC_imode = ! (fHeaderBlock->McFlag());		// fJTC_imode = 1(DATA), 0 (MC)
	std::cout << __FILE__ << "::" << __LINE__ << ":AUTOMATIC SETTING OF MC_FLAG = " << fHeaderBlock->McFlag() << " FOR JTCmodei = " << fJTC_imode << "(1=DATA, 0=MC)"  << std::endl;
	return 0;
}

//_____________________________________________________________________________
int JetFilterModuleV2::Event(int ientry)
{
	SetPassed(0); 
	mycounter.evtsRunOver++;
	if (!bRunPermit) {
		StdOut(__FILE__,__LINE__,3,"One or more dependencies not found. pl check. Run Permit not cleared!");
		SetPassed(0);
		exit (1);
		return 0;
	}

	//GetGenLevelInfo();

	bSaveThisEvent = false;
	EventCount_b++;
	ClearModuleOutput();
	SigMetEvent_status=0; // re-setting status code for Significan Met event
	bad_EMjet_match_flag=1; // by default no events flagged "bad" in the beginning
	bool pass   = true;

	myJetRun=GetHeaderBlock()->RunNumber(); // obtaining run number, myJetRun is globaly defined
	if (debug)
	{
		std::cout << "================" << myJetRun << ", " << GetHeaderBlock()->EventNumber() << std::endl; 
	}
	//std::cout << "================" << myJetRun << ", " << GetHeaderBlock()->EventNumber() << std::endl; 

	//_________________________________________________________________________
	//---------  Accessing Jet Blocks
	fJetBlockClu04->GetEntry(ientry);
	fHadJetBlockClu04->GetEntry(ientry);
	fPhotonBlock->GetEntry(ientry);
	fCalData->GetEntry(ientry);
	fElectronBlock->GetEntry(ientry);
	fMetBlock->GetEntry(ientry);	// sam 03-14-2008
	fGenpBlock->GetEntry(ientry);

	//------------------------------------- collect vertex info for later use
	std::vector<double> zvxvec;
	double zvx_best = trigMod->GetVtxZ(1,1);		//1st class 12 vertex-sam
	double zvx_2nd  = -1000.0;
	myNvx_class12     = trigMod->GetNClass12Vtx();
	double dzvx_2best = trigMod->GetBestVtxZsep(); // Zsep of two best class vertices, return 0 if NClass12Vtx <2 -- sam
	int myNvx_all = trigMod->GetNVtx();
	if(myNvx_all>1) zvx_2nd=trigMod->GetVtxZ(2,1);		//2nd best class 12 vtx
	double dzvx_best2worse = -1000.0; // worse Zsep between best any other vertex
	double zvx_worse = -1000.0; // largest Zvx
	for (int i=0; i<myNvx_all; i++)
	{
		zvxvec.push_back(trigMod->GetVtxZ(i));
		double _dz = zvx_best - (trigMod->GetVtxZ(i));
		if(fabs(_dz)>dzvx_best2worse) dzvx_best2worse=fabs(_dz);
		if(fabs(trigMod->GetVtxZ(i)) > zvx_worse) zvx_worse = fabs(trigMod->GetVtxZ(i));
	}
	double _zvx[3];
	_zvx[0]=zvx_best;
	_zvx[1]=zvx_2nd;
	_zvx[2]=zvx_worse;

	DoCommonStuff(allstuff); // reads raw Met, pho, ele info and fills CommonStuff
	DoMyJet(fJetBlockClu04,allstuff,jet04stuff,fMatchPhoJet04); // does my JetClu04 jets
	DoMyMet(allstuff,jet04stuff); // corrects Met & SumEt for JetClu04 jets

	// I am adding this condition here to avoid unneccessary calculation when making
	// flat ntuples. All I need is list of corrected Jets and MET/SUMET- Sam, Sep 07,2009
	if (fAnalysisMode ==0 || fAnalysisMode == 1)
	{


		MyNewJetWithMet(jet04stuff, allstuff); // adding jet and Met

		//-------------------------------------------------------------------------------
		// Filling histograms before event cuts >>>>>>>>
		//_______________________________________________________________________________
		FillJetHistogramsB(fHistJet04,jet04stuff,zvx_best-fJetBlockClu04->ZVertex(),myNvx_class12); // filling general histo for JetClu04

		//-------------------------------------------------------------------------------
		// Applying event Cuts
		//_______________________________________________________________________________

		//....... temporary cuts
		pass=false;

		if(MyAngularCut(jet04stuff,allstuff)==1) pass=true; // added on 03/08/06 (temporary)
		else
		{
			StdOut(__FILE__,__LINE__,0," Angular cut failed. ");
			SetPassed(0);
			return 0;
		}

		int Nduplicate=MyDuplicateCut(allstuff);
		if(fRemoveDuplicate==1 && Nduplicate>0) // removing duplicates
		{
			StdOut(__FILE__,__LINE__,0," Duplicate cut failed. ");
			SetPassed(0);
			return 0;
		}
		else pass=true;
		if((fRemoveDuplicate+Nduplicate)<0) // selecting duplicates
		{
			SetPassed(0);
			return 0;
		}
		else pass=true;

		if(pass==true && fAnalysisMode==1) DoMyAnalysis(jet04stuff,allstuff); // filling analysis histo; needs to be before MetCleanUpCut  //delete fAnalysisMode check in DoMyAnalysis() 01-29-2009 sam

		//   int Nbadmet=MyMetCleanUpCut(jet04stuff,allstuff);
		int Nbadmet=MyMetCleanUpCut(jet04stuff,allstuff,jet04stuff.Jet_lev6_noEMobj,jet04stuff.myMETcorr_th15); // new as of 04/04/07
		if(fRemoveBadMet==1 && Nbadmet>0) // removing events with bad MET
		{
			SetPassed(0);
			return 0;
		}
		else pass=true;
		if((fRemoveBadMet+Nbadmet)<0) // selecting events with bad MET
		{
			SetPassed(0);
			return 0;
		}
		else pass=true;

		if(pass)
		{
			SetPassed(1);

			//-------------------------------------------------------------------------------
			// Filling histos for events which pass the cuts
			//_______________________________________________________________________________
			FillJetHistogramsA(fHistJet04,jet04stuff,zvx_best-fJetBlockClu04->ZVertex(),myNvx_class12); // filling general histo for JetClu04
			if(fAnalysisMode==0) // tmp, introduced on 01/02/08
			{
				FillMetCleanupHistograms(fCleanup,jet04stuff,allstuff,3); // histo for met cleanup studies
				FillMetStudyHistograms(fMetStudyJet04sc3,jet04stuff,allstuff,3,myNvx_class12,dzvx_2best,dzvx_best2worse,_zvx); // met study histo for JetClu 0.4 and Et>15 GeV
				MyMetEventCount(jet04stuff,allstuff,3,met_results); // counting number of data events with Met>cut
				FillPhoMetStudyHistograms(fPhoMetStudy,jet04stuff.myMETcorr_th15,allstuff,zvx_best);
			}

			/*
			//sam//
			int exo_event=MyExoticEvent(jet04stuff,allstuff,3,GetHeaderBlock()->RunNumber(),GetHeaderBlock()->EventNumber());
			int veryexo_event=MyVeryExoticEvent(jet04stuff,allstuff,3,GetHeaderBlock()->RunNumber(),GetHeaderBlock()->EventNumber(),myNvx_class12,zvxvec);
			*/
			int met_event=MyLargeMetEvent(jet04stuff,allstuff,3,myNvx_class12,GetHeaderBlock()->RunNumber(),GetHeaderBlock()->EventNumber());

			//       MyDumpEvent(jet04stuff,allstuff,myNvx_class12,GetHeaderBlock()->RunNumber(),GetHeaderBlock()->EventNumber());

			//sam//
			if(fSelectMetEvent==1 && met_event==0)
			{
				pass=false;
				SetPassed(0);
			}
			//    if(fSelectExoticEvent==1 && exo_event==0)
			//{
			//  pass=false;
			//  SetPassed(0);
			//	}

			//sam//
			//_____ selecting events according to MetSignificance
			if(myMetSig>=fMetSig_cut) SigMetEvent_status=1; // new as of 12/09/07
			if(fSelectSigMetEvent==1 && SigMetEvent_status==0)
			{
				pass=false;
				SetPassed(0); // selecting event with Significant Met
			}
			if(fSelectSigMetEvent==-1 && SigMetEvent_status==1)
			{
				pass=false;
				SetPassed(0); // rejecting event with Significant Met
			}
		}
		else SetPassed(0);
		
	} else   
	{
		// I want to pass many events as possible when making flat ntuples
		// I'll pass events only if there is a jet with
		// Et>MinEt and Eta<MaxEta -- Sam, Sep 7,2009
		
		if (MyNjetCut(jet04stuff)==1)
		{

			for (unsigned int i=0; i<jet04stuff.Jet_lev7_noEMobj.size(); i++)
			{
				if (fabs(jet04stuff.Jet_lev7_noEMobj.at(i).Eta())<=MaxJetEta)
				{
					SetPassed(1);
					break;
				}
			}
		}
	} // if (fAnalysisMode=0 or 1)

	if(pass==true)
	{
		EventCount_a++;
		MyDumpEvent(jet04stuff,allstuff,myNvx_class12,GetHeaderBlock()->RunNumber(),GetHeaderBlock()->EventNumber());
	}

	if (GetPassed()) 	mycounter.evtsPassModule++;


	//std::cout << "======= END EVENT - ";
	//GetHeaderBlock()->Print();
	return 0;
}

//_____________________________________________________________________________
void JetFilterModuleV2::Display() {

	return;
}

//________ calls Sumw2 for final analysis histograms
void JetFilterModuleV2::FinalAnalysisHistoStep1(AnalysisHisto_t& Hist) {
	Hist.fAna_MetAll->Sumw2();
	Hist.fAna_MetSig->Sumw2();
	Hist.fAna_Met->Sumw2();
	Hist.fAna_Njet15->Sumw2();
	Hist.fAna_M->Sumw2();
	Hist.fAna_dPhi->Sumw2();
	Hist.fAna_Njet20->Sumw2();
	Hist.fAna_Njet25->Sumw2();
	Hist.fAna_Njet30->Sumw2();
	Hist.fAna_Njet35->Sumw2();
	Hist.fAna_Qt->Sumw2();
	Hist.fAna_Et1->Sumw2();
	Hist.fAna_Et2->Sumw2();
	Hist.fAna_Etjet->Sumw2();
	Hist.fAna_Ht->Sumw2();
	Hist.fAna_Mjj->Sumw2();
	Hist.fAna_Nem->Sumw2();
	Hist.fAna_Mej->Sumw2();
	Hist.fAna_Mextra->Sumw2();
	Hist.fAna_Etem->Sumw2();
	Hist.fAna_dPhi1->Sumw2();
	Hist.fAna_dPhi2->Sumw2();
	Hist.fAna_dPhi3->Sumw2();
	Hist.fAna_dPhi4->Sumw2();
	Hist.fAna_dPhi5->Sumw2();
	Hist.fAna_dPhi6->Sumw2();
	Hist.fAna_dPhi7->Sumw2();
	return;
}

//________ normalizes final background analysis histograms to Nevents
//ok now dividing by npoints is wrong. Not all the expts will have entries(if jets flutuated down)
//I am normizlaing them to data before met sig cut. - sam
void JetFilterModuleV2::FinalAnalysisHistoStep2(AnalysisHisto_t& bckg, const AnalysisHisto_t& data)
{

	if (Npoints>0)
	{
		//First get the scale factor for all other hists using either MetAll or MetSig before MetSig cut
		//These are the only two hists that are filled before the MetSig cut
		//All bckg hists must be scaled by this irrespective of the MetSig cut
		double dDataIntegral = data.fAna_MetAll->Integral();
		double dBckgIntegral = bckg.fAna_MetAll->Integral();
		float fScale = (float) dDataIntegral/dBckgIntegral;

		bckg.fAna_MetAll->Scale(fScale);
		bckg.fAna_MetSig->Scale(fScale);
		bckg.fAna_Met->Scale(fScale);
		bckg.fAna_Njet15->Scale(fScale);

		//std::cout << "Integrals a4 data, bckg, scale = " << data.fAna_MetAll->Integral() << "," << bckg.fAna_MetAll->Integral() << ", " << fScale << std::endl;
		bckg.fAna_M->Scale(fScale);
		bckg.fAna_dPhi->Scale(fScale);
		bckg.fAna_Njet20->Scale(fScale);
		bckg.fAna_Njet25->Scale(fScale);
		bckg.fAna_Njet30->Scale(fScale);
		bckg.fAna_Njet35->Scale(fScale);
		bckg.fAna_Qt->Scale(fScale);
		bckg.fAna_Et1->Scale(fScale);
		bckg.fAna_Et2->Scale(fScale);
		bckg.fAna_Etjet->Scale(fScale);
		bckg.fAna_Ht->Scale(fScale);
		bckg.fAna_Mjj->Scale(fScale);
		bckg.fAna_Nem->Scale(fScale);
		bckg.fAna_Mej->Scale(fScale);
		bckg.fAna_Mextra->Scale(fScale);
		bckg.fAna_Etem->Scale(fScale);
		bckg.fAna_dPhi1->Scale(fScale);
		bckg.fAna_dPhi2->Scale(fScale);
		bckg.fAna_dPhi3->Scale(fScale);
		bckg.fAna_dPhi4->Scale(fScale);
		bckg.fAna_dPhi5->Scale(fScale);
		bckg.fAna_dPhi6->Scale(fScale);
		bckg.fAna_dPhi7->Scale(fScale);
	}

	return;
}
//________ multiplies final analysis histograms by event weight
void JetFilterModuleV2::FinalAnalysisHistoStep3(AnalysisHisto_t& Hist) {

	TF1* h_norm=new TF1("f1","1.0",-1000.0,10000.0);
	Hist.fAna_MetAll->Multiply(h_norm,fEventWeight);
	Hist.fAna_MetSig->Multiply(h_norm,fEventWeight);
	Hist.fAna_Met->Multiply(h_norm,fEventWeight);
	Hist.fAna_Njet15->Multiply(h_norm,fEventWeight);
	Hist.fAna_M->Multiply(h_norm,fEventWeight);
	Hist.fAna_dPhi->Multiply(h_norm,fEventWeight);
	Hist.fAna_Njet20->Multiply(h_norm,fEventWeight);
	Hist.fAna_Njet25->Multiply(h_norm,fEventWeight);
	Hist.fAna_Njet30->Multiply(h_norm,fEventWeight);
	Hist.fAna_Njet35->Multiply(h_norm,fEventWeight);
	Hist.fAna_Qt->Multiply(h_norm,fEventWeight);
	Hist.fAna_Et1->Multiply(h_norm,fEventWeight);
	Hist.fAna_Et2->Multiply(h_norm,fEventWeight);
	Hist.fAna_Etjet->Multiply(h_norm,fEventWeight);
	Hist.fAna_Ht->Multiply(h_norm,fEventWeight);
	Hist.fAna_Mjj->Multiply(h_norm,fEventWeight);
	Hist.fAna_Nem->Multiply(h_norm,fEventWeight);
	Hist.fAna_Mej->Multiply(h_norm,fEventWeight);
	Hist.fAna_Mextra->Multiply(h_norm,fEventWeight);
	Hist.fAna_Etem->Multiply(h_norm,fEventWeight);
	Hist.fAna_dPhi1->Multiply(h_norm,fEventWeight);
	Hist.fAna_dPhi2->Multiply(h_norm,fEventWeight);
	Hist.fAna_dPhi3->Multiply(h_norm,fEventWeight);
	Hist.fAna_dPhi4->Multiply(h_norm,fEventWeight);
	Hist.fAna_dPhi5->Multiply(h_norm,fEventWeight);
	Hist.fAna_dPhi6->Multiply(h_norm,fEventWeight);
	Hist.fAna_dPhi7->Multiply(h_norm,fEventWeight);
	delete h_norm;
	return;
}

//_____________________________________________________________________________
//  END JOB SUMMARY
//_____________________________________________________________________________
int JetFilterModuleV2::EndJob() {

	if (GetSummaryStat()) return 0;

	std::string sMsg;
	if (GetJTC_imode() == 0) sMsg = "MC";
	else if (GetJTC_imode() == 1) sMsg = "DATA";
	else sMsg = "UNKNOWN";

	std::string sUnclParam;
	if (GetUnclParamSwitch() == 0) sUnclParam = "gg sideband";
	else if (GetUnclParamSwitch() == 1) sUnclParam = "Zee";
	else sUnclParam = "UNKNOWN!";

	std::string sSampleId;
	if (GetSampleID() == 0) sSampleId = "ggX signal";
	else if (GetSampleID() == 1) sSampleId = "ggX sideband";
	else if (GetSampleID() == 2) sSampleId = "Z->ee";
	else if (GetSampleID() == 3)
	{	if (GetJTC_imode() == 0) sSampleId = "PYTHIA pho+jet";
		else if (GetJTC_imode() == 1) sSampleId = "DATA pho+jet";
		else sSampleId = "Not Labeled!";
	}
	else if (GetSampleID() == 4) sSampleId = "DATA sideband pho+jet";
	else if (GetSampleID() == 10) sSampleId = "Calibration";
	else sSampleId = "Not Labeled!";

	printf("[TMJ:00:]----- end job: ---- %s\n",GetName());
	std::cout << "[TMJ:01:] Events Processed ----------- = " << mycounter.evtsRunOver << std::endl;
	std::cout << "[TMJ:02:] Events Passed -------------- = " << mycounter.evtsPassModule << std::endl;
	
	std::string njetcut("");
	if (GetMinNjet5()>0) njetcut += "[5 GeV, >=" + ToStr(GetMinNjet5()) + " & =<" +  ToStr(GetMaxNjet5()) + "] ";
	if (GetMinNjet10()>0) njetcut += "[10 GeV, >=" + ToStr(GetMinNjet10()) + " & =<" +  ToStr(GetMaxNjet10()) + "] ";
	if (GetMinNjet15()>0) njetcut += "[15 GeV, >=" + ToStr(GetMinNjet15()) + " & =<" +  ToStr(GetMaxNjet15()) + "] ";
	if (GetMinNjet20()>0) njetcut += "[20 GeV, >=" + ToStr(GetMinNjet20()) + " & =<" +  ToStr(GetMaxNjet20()) + "] ";
	if (GetMinNjet25()>0) njetcut += "[25 GeV, >=" + ToStr(GetMinNjet25()) + " & =<" +  ToStr(GetMaxNjet25()) + "] ";
	if (GetMinNjet30()>0) njetcut += "[30 GeV, >=" + ToStr(GetMinNjet30()) + " & =<" +  ToStr(GetMaxNjet30()) + "] ";
	if (GetMinNjet35()>0) njetcut += "[35 GeV, >=" + ToStr(GetMinNjet35()) + " & =<" +  ToStr(GetMaxNjet35()) + "] ";
	if (GetMinNjet5()==0 && GetMinNjet10()==0 && GetMinNjet15()==0 && GetMinNjet20()==0 && GetMinNjet25()==0
			&& GetMinNjet30()==0 && GetMinNjet35() ==0) njetcut += "No Jet Cut";
	
	std::cout << "[TMJ:03:] (Et> - Min/Max) NJets ------ = " << njetcut << std::endl;
	std::cout << "[TMJ:04:] Max Jet Detector Eta         = " << GetMaxJetEta() <<std::endl;
	if (GetAnalysisMode() ==0 || GetAnalysisMode() == 1)
	{
	std::cout << "[TMJ:05:] DatFile Name --------------- = " << GetDatFileName() << std::endl;
	std::cout << "[TMJ:06:] LargeMetEventFile Name ----- = " << GetLargeMetEventFileName() << std::endl;
	std::cout << "[TMJ:07:] DumpEventFile Name --------- = " << GetDumpEventFileName() << std::endl;
	}
	std::cout << "[TMJ:08:] JTC_imode / Syst Code ------ = " << GetJTC_imode() << " (" << sMsg << ") / " << GetJTC_systcode()  << std::endl;
	if (GetAnalysisMode() ==0 || GetAnalysisMode() == 1)
	{
	std::cout << "[TMJ:09:] Met Sig Cut ---------------- = " << GetMetSigCut() << std::endl;
	std::cout << "[TMJ:10:] Select Sig Met Evt --------- = " << GetSelectSigMetEvent() << std::endl;
	std::cout << "[TMJ:11:] Uncl Para Switch ----------- = " << GetUnclParamSwitch() << " (" << sUnclParam << ")" << std::endl;
	std::cout << "[TMJ:12:] Sample ID ------------------ = " << GetSampleID()  << " (" << sSampleId << ")" << std::endl;
	std::cout << "[TMJ:13:] N points generated  -------- = " << GetNpointsToGenerate() << std::endl;
	}
	std::string sAnaMode;
	if (GetAnalysisMode() == 0) sAnaMode = "OFF - Calibration Mode";
	else if (GetAnalysisMode() == 1) sAnaMode = "ON - Analysis Mode";
	else sAnaMode = "Generate Corr Jets/MEt/SumEt only";
	std::cout << "[TMJ:14:] Analysis mode -------------- = " << GetAnalysisMode() << " - " << sAnaMode << std::endl;

	if (GetAnalysisMode() ==0 || GetAnalysisMode() == 1)
	{
	std::cout << "[TMJ:15:] Remove Duplicate ----------- = " << GetRemoveDuplicate() << std::endl;
	std::cout << "[TMJ:16:] Remove Bad Met  ------------ = " <<  GetRemoveBadMet()<< std::endl;
	}

	if (fAnalysisMode==1)
	{
		int N_data=(int) fAna_data.fAna_Met->GetEntries();
		double N_qcd=0.0;
		double N_qcd_stat=0.0;
		double N_qcd_syst=0.0;
		double dN_syst[10];

		if (Npoints>0)
		{
			//N_qcd=1.0*fAna_def.fAna_Met->GetEntries()/(1.0*Npoints);		// change this to Integral
			//this is temporily until I come up with a new method to renormalize all
			//hists (inclduing syst hists) to match data. 04-10-2009
			double dDataIntegral = fAna_data.fAna_MetAll->Integral();
			double dBckgIntegral = fAna_bckg.fAna_MetAll->Integral();
			float fScale = (float) dDataIntegral/dBckgIntegral;

			N_qcd=(1.0*fAna_def.fAna_Met->Integral()/(1.0*Npoints))*fScale;		// change this to Integral
			N_qcd_stat=sqrt(1.0*fAna_def.fAna_Met->GetEntries())/(1.0*Npoints);
			dN_syst[0]=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_ue1.fAna_Met->GetEntries()));

			double d1;
			double d2;
			d1=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_ue2.fAna_Met->GetEntries()));
			d2=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_ue3.fAna_Met->GetEntries()));
			dN_syst[1]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_ue4.fAna_Met->GetEntries()));
			d2=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_ue5.fAna_Met->GetEntries()));
			dN_syst[2]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_ue6.fAna_Met->GetEntries()));
			d2=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_ue7.fAna_Met->GetEntries()));
			dN_syst[3]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_ue8.fAna_Met->GetEntries()));
			d2=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_ue9.fAna_Met->GetEntries()));
			dN_syst[4]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer1.fAna_Met->GetEntries()));
			d2=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer2.fAna_Met->GetEntries()));
			dN_syst[5]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer3.fAna_Met->GetEntries()));
			d2=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer4.fAna_Met->GetEntries()));
			dN_syst[6]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer5.fAna_Met->GetEntries()));
			d2=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer6.fAna_Met->GetEntries()));
			dN_syst[7]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer7.fAna_Met->GetEntries()));
			d2=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer8.fAna_Met->GetEntries()));
			dN_syst[8]= d1>d2 ? d1 : d2;
			
			d1=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer9.fAna_Met->GetEntries()));
			d2=fabs(1.0*(fAna_def.fAna_Met->GetEntries()-fAna_jer10.fAna_Met->GetEntries()));
			dN_syst[9]= d1>d2 ? d1 : d2;
			
			for(int i=0; i<10; i++)
			{
				N_qcd_syst=N_qcd_syst+dN_syst[i]*dN_syst[i];
			}
			N_qcd_syst=sqrt(N_qcd_syst)/(1.0*Npoints);
		}

		std::cout << "[TMJ:--:]---------- Analysis Summary  --------------"<<std::endl;
		std::cout << "[TMJ:17:] Data observed -------------- = " << N_data<<std::endl;
		std::cout << "[TMJ:18:] MetModel prediction -------- = " << N_qcd<<" +- "<<N_qcd_stat<<" +- "<<N_qcd_syst<<std::endl;
		std::cout << "[TMJ:19:] Number of events before cuts = " << EventCount_b<<std::endl;
		std::cout << "[TMJ:20:] Number of events after cuts  = " << EventCount_a<<std::endl;
		
		FinalAnalysisHistoStep1(fAna_bckg);				//calls SumW2() for all hists
		// I need to double check this systematic calcuation with Sasha. --02-10-2009
		// I moved all this to the final MetHistMaking macro MergeMetHists.C
		// This will solve the problem of weighting different job section differently when
		// I split the job to run on CAF.

		//DoFinalAnalysisHisto(fAna_bckg,fAna_ue1,fAna_def); // ok. here u r picking relative uncertainty
		//DoFinalAnalysisHisto(fAna_bckg,fAna_ue2,fAna_ue3); // starting from def
		//DoFinalAnalysisHisto(fAna_bckg,fAna_ue4,fAna_ue5);
		//DoFinalAnalysisHisto(fAna_bckg,fAna_ue6,fAna_ue7);
		//DoFinalAnalysisHisto(fAna_bckg,fAna_ue8,fAna_ue9);
		//DoFinalAnalysisHisto(fAna_bckg,fAna_jer1,fAna_jer2);  // but not here why not start from def??
		//DoFinalAnalysisHisto(fAna_bckg,fAna_jer3,fAna_jer4);
		//DoFinalAnalysisHisto(fAna_bckg,fAna_jer5,fAna_jer6);
		//DoFinalAnalysisHisto(fAna_bckg,fAna_jer7,fAna_jer8);
		//DoFinalAnalysisHisto(fAna_bckg,fAna_jer9,fAna_jer10);
		//FinalAnalysisHistoStep2(fAna_bckg, fAna_data);
		
		if(fUseEventWeight==1)
		{
			FinalAnalysisHistoStep3(fAna_data);
			FinalAnalysisHistoStep2(fAna_def, fAna_data); // Do I need to do this for def?
			FinalAnalysisHistoStep3(fAna_def);
			FinalAnalysisHistoStep3(fAna_bckg);
		}


		if (debug)
		{
			std::cout << __FUNCTION__ << "::" << __LINE__ << "::data hist fAna_njet15" << std::endl;

			for (int i = 0; i <= fAna_data.fAna_Njet15->GetNbinsX() + 1; ++i)
			{
				if (fAna_data.fAna_Njet15->GetBinContent(i))
				{
					std::cout << "bin=" << i << "\tloedge=" << fAna_data.fAna_Njet15->GetBinCenter(i) << "\tcontent= "  << fAna_data.fAna_Njet15->GetBinContent(i) << std::endl;

				}
			}
			std::cout << __FUNCTION__ << "::" << __LINE__ << "::bckg hist fAna_njet15" << std::endl;
			for (int i = 0; i <= fAna_bckg.fAna_Njet15->GetNbinsX() + 1; ++i)
			{
				if (fAna_bckg.fAna_Njet15->GetBinContent(i))
				{
					std::cout << "bin=" << i << "\tloedge=" << fAna_bckg.fAna_Njet15->GetBinCenter(i) << "\tcontent= "  << fAna_bckg.fAna_Njet15->GetBinContent(i) << std::endl;

				}
			}
		}

		
	} else if (fAnalysisMode == 0)
	{
		DoFinalPhoMetHisto(fPhoMetStudy);
		FinalMetSigCorrHisto(fMetStudyJet04sc3);
	}
	
	
	if (GetAnalysisMode() ==0 || GetAnalysisMode() == 1)
	{
		std::cout << "[TMJ:23:] Using Hadron Jets?           = " << ToStr(GetUseHadJets()) << " (1=YES, 0=NO)" <<std::endl;
		std::cout << "[TMJ:24:] More HAD jets than DET jets  = " << ToStr(mycounter.evtsMoreHad) << std::endl;
		std::cout << "[TMJ:24:] Additional Smearing?         = " << ToStr(DoubleSmearJets()) << " (1=YES, 0=NO)" << std::endl;
		std::cout << "[TMJ:25:] JER Offset ?                 = " << ToStr(ZeroJERoffset()) << " (1=NO offset, 0=With offset)" << std::endl;
		std::cout << "[TMJ:26:] MEt cut >                    = " << ToStr(MetCut())  << " GeV"<< std::endl;
	}
		std::cout << "[TMJ:27:] Print Level                  = " << PrintLevel() << std::endl;

	if (GetAnalysisMode() ==0 || GetAnalysisMode() == 1)
	{
		std::string sMtd("UNKNOWN!");
		if (MEtAddMethod() == DO_NOT_ADD_MET) sMtd = "NOT ADDED";
		else if (MEtAddMethod() == TOMOST_MISMEASURED_JET) sMtd = "Most Mismeasured jet";
		else if (MEtAddMethod() == TO_JET_CLOSEST_TO_MET) sMtd = "Jet Closest to MEt in DelPhi(Jet,MEt)<0.4";
		
		std::cout << "[TMJ:28:] MEt Added Method             = " << MEtAddMethod() << " (" << sMtd << ")"<< std::endl;
	}
		
		std::string sEmObj;
		if (RemoveTightPho()) sEmObj += "[Tight Pho]";
		if (RemoveLoosePho()) sEmObj += "[Loose Pho]";
		//if (RemoveSidebandPho()) sEmObj += "[Sideband Pho]";
		if (RemoveTightPhoLikeEle()) sEmObj += "[Tight PhoLikeEle]";
		if (RemoveLoosePhoLikeEle()) sEmObj += "[Loose PhoLikeEle]";
		if (RemoveStdLooseEle()) sEmObj += "[Std. Loose Ele]";
		if (! (RemoveTightPho() || RemoveLoosePho() || RemoveTightPhoLikeEle()
				|| RemoveLoosePhoLikeEle() || RemoveStdLooseEle()) ) sEmObj += "None";
		std::cout << "[TMJ:30:] Types of EM objs removed     = " << sEmObj << std::endl;

	if (GetAnalysisMode() ==0 || GetAnalysisMode() == 1)
	{
		//info about the JER functions used new or old
		std::string sJerMsg;
		if (JERForGaussMean() == JER_GAUSSMEAN) sJerMsg += "[Gmean-NEW]";
		else if (JERForGaussMean() == JER_OLDFITS) sJerMsg +="[Gmean-OLD]";
		else sJerMsg +="[Gmean-UNKNOWN]";

		if (JERForGaussSigma() == JER_GAUSSSIGMA) sJerMsg += "[Gsigma-NEW]";
		else if (JERForGaussSigma() == JER_OLDFITS) sJerMsg += "[Gsigma-OLD]";
		else sJerMsg +="[Gsigma-UNKNOWN]";

		if (JERForLandauMPV() == JER_LANDAUMPV) sJerMsg += "[Lmean-NEW]";
		else if (JERForLandauMPV() == JER_OLDFITS) sJerMsg += "[Lmpv-OLD]";
		else sJerMsg +="[Lmpv-UNKNOWN]";

		if (JERForLandauSigma() == JER_LANDAUSIGMA) sJerMsg += "[Lsigma-NEW]";
		else if (JERForLandauSigma() == JER_OLDFITS) sJerMsg += "[Lsigma-OLD]";
		else sJerMsg +="[Gnorm-UNKNOWN]";

		if (JERForGaussNorm() == JER_GAUSSNORM) sJerMsg += "[Gnorm-NEW]";
		else if (JERForGaussNorm() == JER_OLDFITS) sJerMsg += "[Gnorm-OLD]";
		else sJerMsg +="[Gnorm-UNKNOWN]";

		std::cout << "[TMJ:31:] JER fits used                = " << sJerMsg << std::endl;
	}		

	if (GetAnalysisMode() == 1)
	{
		MyMetResults(met_results);
	}

	printf("---------------------------------------------------\n");

	return 0;
}

double JetFilterModuleV2::GetCorrection(CorInit settings, TLorentzVector vec,
		float& emf, float etad) {
	int nrun = settings.Nrun;
	int nVertex = settings.nvx;
	int coneSize=settings.cone;
	int version=settings.version;
	int syscode=settings.sys;
	int level=settings.level;
	int imode=settings.imode;

	JetEnergyCorrections myJetEnergyCorrections=
		JetEnergyCorrections("JetCorrections","JetCorrections",
				level,nVertex,coneSize,version,syscode,nrun,imode);

	HepLorentzVector P4Jet;
	P4Jet.setPx(vec.Px());
	P4Jet.setPy(vec.Py());
	P4Jet.setPz(vec.Pz());
	P4Jet.setE(vec.E());
	if(abs(syscode)>0)
	{
		int n_sigma= (syscode>0) ? 1 : -1 ;
		myJetEnergyCorrections.setTotalSysUncertainties(n_sigma);
	}
	double scaleFactor = myJetEnergyCorrections.doEnergyCorrections(P4Jet,emf,etad);
	return scaleFactor;
}

double JetFilterModuleV2::GetCorrection(TLorentzVector vec, float& emf, float etad) {

	int nrun = myJetRun;
	int nVertex = myNvx_class12;
	int coneSize=fJTC_coneSize;
	int version=fJTC_version;
	int syscode=fJTC_systcode;
	int level=fJTC_level;
	int imode=fJTC_imode;

	JetEnergyCorrections myJetEnergyCorrections=
		JetEnergyCorrections("JetCorrections","JetCorrections",
				level,nVertex,coneSize,version,syscode,nrun,imode);

	HepLorentzVector P4Jet;
	P4Jet.setPx(vec.Px());
	P4Jet.setPy(vec.Py());
	P4Jet.setPz(vec.Pz());
	P4Jet.setE(vec.E());
	if(abs(syscode)>0)
	{
		int n_sigma= (syscode>0) ? 1 : -1 ;
		myJetEnergyCorrections.setTotalSysUncertainties(n_sigma);
	}
	double scaleFactor = myJetEnergyCorrections.doEnergyCorrections(P4Jet,emf,etad);
	return scaleFactor;
}

//_____________________________________________________________ Balance of two leading jets
// !!!!!!!!      OLD stuff, but I still need it  !!!!!!!!!
double JetFilterModuleV2::GetBalance2(double pt1, double pt2, double phi1, double phi2) {
	double balX;
	double balY;
	double bal=9999.0;
	balX=pt1*TMath::Cos(phi1)+pt2*TMath::Cos(phi2);
	balY=pt1*TMath::Sin(phi1)+pt2*TMath::Sin(phi2);
	if(fabs(pt1+pt2)>0.0) bal=TMath::Sqrt(balX*balX+balY*balY)/(pt1+pt2);
	return bal;
}
double JetFilterModuleV2::GetBalance2(std::vector<TLorentzVector> vec){
	double bal=9999.0;
	if(vec.size()>=2)
	{
		double pt_1=vec[0].Pt();
		double pt_2=vec[1].Pt();
		double Phi_1=vec[0].Phi();
		double Phi_2=vec[1].Phi();
		bal=GetBalance2(pt_1,pt_2,Phi_1,Phi_2);
	}
	return bal;
}

//________________________________________________ Theta* or dijet scattering angle
//------------------- Input params are eta's or true rapidities
double JetFilterModuleV2::GetThetaStar(double y1, double y2) {
	double tcos;
	double y=(y1-y2)/2.0;
	double t=0.0;
	tcos=(TMath::Exp(y)-TMath::Exp(-y))/(TMath::Exp(y)+TMath::Exp(-y));
	if(fabs(tcos)<=1.0) t=TMath::ACos(tcos);
	return t;
}
double JetFilterModuleV2::GetThetaStar(std::vector<TLorentzVector> vec) {
	double t=0.0;
	if(vec.size()>=2)
	{
		double Y1=vec[0].Rapidity();
		double Y2=vec[1].Rapidity();
		t=GetThetaStar(Y1, Y2);
	}
	return t;
}


//__________________________________________________ returns Kt of two leading jets
double JetFilterModuleV2::GetMyKtKick(std::vector<TLorentzVector> vec) {
	double kt_kick=0.0;
	if(vec.size()>=2)
	{
		double pt_1=vec[0].Pt();
		double pt_2=vec[1].Pt();
		kt_kick=(pt_1+pt_2)*GetBalance2(vec);
	}
	return kt_kick;
}
//__________________________________________________ returns Phi of two leading jets
double JetFilterModuleV2::GetMyPhiKt2L(std::vector<TLorentzVector> vec) {
	double phi_kick=-1.0E6;
	if(vec.size()>=2)
	{
		double pt=(vec[0]+vec[1]).Pt();
		if(pt>0.0) phi_kick=TVector2::Phi_0_2pi((vec[0]+vec[1]).Phi());
	}
	return phi_kick;
}
//__________________________________________________ returns Phi of sumKt of all jets
double JetFilterModuleV2::GetMyPhiKtAll(std::vector<TLorentzVector> vec) {
	double phi_extra=-1.0E6;
	if(vec.size()>0)
	{
		TLorentzVector extra_jets(0.0,0.0,0.0,0.0);
		for(unsigned int i=0; i<vec.size(); i++)
		{
			extra_jets=extra_jets+vec[i];
		}
		double pt=extra_jets.Pt();
		if(pt>0.0) phi_extra=TVector2::Phi_0_2pi(extra_jets.Phi());
	}
	return phi_extra;
}

//__________________________________________________ returns 1 if jet is matched
//                                                   calculates dR, dPhi, dEta between pho & jet
int JetFilterModuleV2::MyMatchPhoJet(int jet_ind, int pho_ind, int ele_ind,TStnJetBlock* fJetBlock,
		TLorentzVector *mypho,TLorentzVector *myjet,
		double &dR_fj, double &dPhi_fj, double &dEta_fj)
{
	int match_code=0;

	if(jet_ind<0) return match_code;
	if(pho_ind<0 && ele_ind<0) return match_code;

	if(pho_ind>=0)
	{
		CalDataArray jet_towers;
		CalDataArray pho_towers;
		MatchCalorTowers(jet_ind,fJetBlock,fCalData,&jet_towers);
		MatchCalorTowers(pho_ind,fPhotonBlock,fCalData,&pho_towers);
		CalDataArrayI tt_em;
		CalDataArrayI tt_jet;
		for(tt_em = pho_towers.begin(); tt_em != pho_towers.end(); tt_em++)
		{
			TCalTower* calt_em = *tt_em;
			int iEta_em = calt_em->IEta();
			int iPhi_em = calt_em->IPhi();
			for(tt_jet = jet_towers.begin(); tt_jet != jet_towers.end(); tt_jet++)
			{
				TCalTower* calt_jet = *tt_jet;
				int iEta_jet = calt_jet->IEta();
				int iPhi_jet = calt_jet->IPhi();
				if(iEta_em==iEta_jet && iPhi_em==iPhi_jet) match_code=1;
				if(match_code==1) break;
			}
			if(match_code==1) break;
		}
	}

	if(ele_ind>=0)
	{
		CalDataArray jet_towers;
		CalDataArray ele_towers;
		MatchCalorTowers(jet_ind,fJetBlock,fCalData,&jet_towers);
		MatchCalorTowers(ele_ind,fElectronBlock,fCalData,&ele_towers);
		CalDataArrayI tt_em;
		CalDataArrayI tt_jet;
		for(tt_em = ele_towers.begin(); tt_em != ele_towers.end(); tt_em++)
		{
			TCalTower* calt_em = *tt_em;
			int iEta_em = calt_em->IEta();
			int iPhi_em = calt_em->IPhi();
			for(tt_jet = jet_towers.begin(); tt_jet != jet_towers.end(); tt_jet++)
			{
				TCalTower* calt_jet = *tt_jet;
				int iEta_jet = calt_jet->IEta();
				int iPhi_jet = calt_jet->IPhi();
				if(iEta_em==iEta_jet && iPhi_em==iPhi_jet) match_code=1;
				if(match_code==1) break;
			}
			if(match_code==1) break;
		}
	}

	TLorentzVector vec=*myjet;
	double dR=mypho->DeltaR(vec);
	dR_fj=dR;
	dPhi_fj=mypho->DeltaPhi(vec);
	dEta_fj=mypho->Eta()-myjet->Eta();
	return match_code;
}

void JetFilterModuleV2::ClearModuleOutput()       // clears the Module output from
{                              // previous event.

	//-------- my new global output parameters
	myMetSig=0.0;
	myMetSig_gen=0.0;
	myGenMetVec.Set(0.0,0.0);
	ClearJetStuff(jet04stuff); // clears JetClu-0.4 stuff
	ClearCommonStuff(allstuff); // clears common stuff
	//----- end of my new output params
	//.............................................................
	//  ClearMatchStuff(matchstuff);

	return;
}

void JetFilterModuleV2::ClearJetStuff(JetStuff &stuff)
{ 
	// clears the JetStuff structures
	stuff.myMETcorr_th5.Set(-1.0E6,-1.0E6);
	stuff.mySumEtCorr_th5=-1.0E6;
	stuff.mySumEtJet_th5=-1.0E6;
	stuff.myMETcorr_th10.Set(-1.0E6,-1.0E6);
	stuff.mySumEtCorr_th10=-1.0E6;
	stuff.mySumEtJet_th10=-1.0E6;
	stuff.myMETcorr_th15.Set(-1.0E6,-1.0E6);
	stuff.mySumEtCorr_th15=-1.0E6;
	stuff.mySumEtJet_th15=-1.0E6;
	stuff.myMETcorr_th20.Set(-1.0E6,-1.0E6);
	stuff.mySumEtCorr_th20=-1.0E6;
	stuff.mySumEtJet_th20=-1.0E6;
	stuff.myNjet=0;
	stuff.myNjet_th5=0;
	stuff.myNjet_th10=0;
	stuff.myNjet_th15=0;
	stuff.myNjet_th20=0;
	stuff.myNjet_th25=0;
	stuff.myNjet_th30=0;
	stuff.myNjet_th35=0;
	stuff.smeared_Njet_th0=0;
	stuff.smeared_Njet_th5=0;
	stuff.smeared_Njet_th10=0;
	stuff.smeared_Njet_th15=0;
	stuff.smeared_Njet_th20=0;
	stuff.smeared_Njet_th25=0;
	stuff.smeared_Njet_th30=0;
	stuff.smeared_Njet_th35=0;
	stuff.Jet_raw.clear();
	stuff.Jet_lev1.clear();
	stuff.Jet_lev4.clear();
	stuff.Jet_lev5.clear();
	stuff.Jet_lev6.clear();
	stuff.Jet_lev7.clear();
	stuff.Jet_raw_noEMobj.clear();
	stuff.Jet_lev1_noEMobj.clear();
	stuff.Jet_lev4_noEMobj.clear();
	stuff.Jet_lev5_noEMobj.clear();
	stuff.Jet_lev6_noEMobj.clear();
	stuff.Jet_lev7_noEMobj.clear();
	stuff.newJetLev6.clear();
	stuff.smeared_newJetLev6_noEMObj.clear();
	stuff.bvMetAdded.clear();
	stuff.iL6MatchedHadJetIndex.clear();
	stuff.fL6NoEm_MisMeasureProb.clear();
	stuff.JetNtrk.clear();
	stuff.JetNtwr.clear();
	stuff.EtaDet.clear();
	stuff.EtaDetCorr.clear();
	stuff.EmFrRaw.clear();
	stuff.EmFrCorr.clear();
	stuff.Nobj_match.clear();
	stuff.Npho_match.clear();
	stuff.Nele_match.clear();
	stuff.Nmu_match.clear();
	stuff.Ntau_match.clear();
	stuff.Nbtag_match.clear();
	stuff.JetBlockInd.clear();
	stuff.smear_factor.clear();
	return;
}
//_________________________ clears parameters which do not depend on Jet Cone size
void JetFilterModuleV2::ClearCommonStuff(CommonStuff &stuff)
{ 
	// clears CommonStuff structures
	stuff.myRawPhoton.clear();
	stuff.myCorrPhoton.clear();
	stuff.myPhoInd.clear();
	stuff.tightId.clear();
	stuff.myPhoEmFr.clear();
	stuff.myPhoEtaDet.clear();
	stuff.myPhoXces.clear();
	stuff.myPhoZces.clear();
	stuff.myRawElectron.clear();
	stuff.myEleInd.clear();
	stuff.myCorrElectron.clear();
	stuff.myEleEmFr.clear();
	stuff.myEleEtaDet.clear();
	stuff.myElePhiDet.clear();
	stuff.myElePprEt.clear();
	stuff.myEleXces.clear();
	stuff.myEleZces.clear();
	stuff.myRawMuon.clear();
	stuff.myCorrMuon.clear();
	stuff.myRawTau.clear();
	stuff.myCorrTau.clear();
	stuff.myRawBjet.clear();
	stuff.myCorrBjet.clear();
	stuff.myMET_raw.Set(-1.0E6,-1.0E6);
	stuff.myMET0_raw.Set(-1.0E6,-1.0E6);
	stuff.mySumEt_raw=-1.0E6;
	stuff.newVertexZ=-1000.0;

	return;
}

void JetFilterModuleV2::ClearMatchStuff() {
	matchstuff.clear();
	return;
}


//__________________________________________________________ accessors to the
//__________________________________________________________ output params
//_________ need to specify which jets to return (variable "cone"): 0--cone 0.4; 1--cone 0.7; 2--cone 1.0

TLorentzVector* JetFilterModuleV2::GetMyJet_raw(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_raw.size()) return &_stuff.Jet_raw.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJet_lev1(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev1.size()) return &_stuff.Jet_lev1.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJet_lev4(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev4.size()) return &_stuff.Jet_lev4.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJet_lev5(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev5.size()) return &_stuff.Jet_lev5.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJet_lev6(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev6.size()) return &_stuff.Jet_lev6.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJetNoEMobj_raw(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_raw_noEMobj.size()) return &_stuff.Jet_raw_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJetNoEMobj_lev1(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev1_noEMobj.size()) return &_stuff.Jet_lev1_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJetNoEMobj_lev4(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev4_noEMobj.size()) return &_stuff.Jet_lev4_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJetNoEMobj_lev5(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev5_noEMobj.size()) return &_stuff.Jet_lev5_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJetNoEMobj_lev6(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev6_noEMobj.size()) return &_stuff.Jet_lev6_noEMobj.at(i);
	else return NULL;
}
TLorentzVector* JetFilterModuleV2::GetMyJetNoEMobj_lev7(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev7_noEMobj.size()) return &_stuff.Jet_lev7_noEMobj.at(i);
	else return NULL;
}

double JetFilterModuleV2::GetMyJetEtaDet(int cone, int i) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EtaDet.size()) return _stuff.EtaDet.at(i);
	else return -1.0E6;
}
double JetFilterModuleV2::GetMyJetEtaDetCorr(int cone, int i) { // after removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EtaDetCorr.size()) return _stuff.EtaDetCorr.at(i);
	else return -1.0E6;
}
double JetFilterModuleV2::GetMyJetEmFrRaw(int cone, int i) {  // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EmFrRaw.size()) return _stuff.EmFrRaw.at(i);
	else return -1.0E6;
}
double JetFilterModuleV2::GetMyJetEmFrCorr(int cone, int i) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EmFrCorr.size()) return _stuff.EmFrCorr.at(i);
	else return -1.0E6;
}
int JetFilterModuleV2::GetMyJetNtrk(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetNtrk.size()) return _stuff.JetNtrk.at(i);
	else return -1;
}
int JetFilterModuleV2::GetMyJetNtwr(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetNtwr.size()) return _stuff.JetNtwr.at(i);
	else return -1;
}
int JetFilterModuleV2::GetMyJetNobjMatch(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nobj_match.size()) return _stuff.Nobj_match.at(i);
	else return -1;
}
int JetFilterModuleV2::GetMyJetNphoMatch(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Npho_match.size()) return _stuff.Npho_match.at(i);
	else return -1;
}
int JetFilterModuleV2::GetMyJetNeleMatch(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nele_match.size()) return _stuff.Nele_match.at(i);
	else return -1;
}
int JetFilterModuleV2::GetMyJetNmuMatch(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nmu_match.size()) return _stuff.Nmu_match.at(i);
	else return -1;
}
int JetFilterModuleV2::GetMyJetNtauMatch(int cone, int i) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Ntau_match.size()) return _stuff.Ntau_match.at(i);
	else return -1;
}
int JetFilterModuleV2::GetMyJetNbtagMatch(int cone, int i)
{
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nbtag_match.size()) return _stuff.Nbtag_match.at(i);
	else return -1;
}
int JetFilterModuleV2::GetMyJetBlockInd(int cone, int i)
{ 
	// original jet index in JetBlock, need this after jet reordering
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetBlockInd.size()) return _stuff.JetBlockInd.at(i);
	else return -1;
}
int JetFilterModuleV2::GetMyNjet(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25,
	// 5--et>30, 5--et>35
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return _stuff.myNjet;
	if(threshold==1) return _stuff.myNjet_th5;
	if(threshold==2) return _stuff.myNjet_th10;
	if(threshold==3) return _stuff.myNjet_th15;
	if(threshold==4) return _stuff.myNjet_th20;
	if(threshold==5) return _stuff.myNjet_th25;
	if(threshold==6) return _stuff.myNjet_th30;
	if(threshold==7) return _stuff.myNjet_th35;
	//if(threshold<0 || threshold>7) return -1;
	return -1;
}
double JetFilterModuleV2::GetMySumEtCorr(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.mySumEt_raw; // for consistency
	if(threshold==1) return _stuff.mySumEtCorr_th5;
	if(threshold==2) return _stuff.mySumEtCorr_th10;
	if(threshold==3) return _stuff.mySumEtCorr_th15;
	if(threshold==4) return _stuff.mySumEtCorr_th20;
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}
double JetFilterModuleV2::GetMyHtCorr(int cone, int threshold, const int iJetCorrLev)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	// CorrLev decides the lev of corrected jet to be used to calcualte Ht.
	// default would be Lev6 as used through out MEtModel.
	// I may need use lev7 corrected jets(BE CAREFULL! need to change how SumEt 
	// is calculated too). Besides I had to do this chage since
	// I am using the MyHtAll() for both DATA/BG calcualtions and in each time 
	// a different set of jets will be passed in. Sep7,2009-Sam
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;

	double Ht= -1.0E6;
	if (iJetCorrLev==6) Ht = MyHtAll(_stuff.Jet_lev6_noEMobj, _stuff,allstuff,threshold);
	else if (iJetCorrLev==7) 
	{
		//Ht = MyHtAll(_stuff.Jet_lev7_noEMobj, _stuff,allstuff,threshold);
		std::cout << __FILE__ << ":" << __LINE__ << 
			":: Not implemented yet. You need chage MEt/SumEt calculations to use Lev7 jets too!" 
			<< std::endl;
		exit (1);
	} else
	{
		StdOut(__FILE__,__LINE__,3,
		" Not yet defined level of corrected jets requested in GetMyHtCorr()!");
		exit (1);
	}

	return Ht;
}

double JetFilterModuleV2::GetMyMetCorr(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.myMET_raw.Mod(); // for consistency
	if(threshold==1) return _stuff.myMETcorr_th5.Mod();
	if(threshold==2) return _stuff.myMETcorr_th10.Mod();
	if(threshold==3) return _stuff.myMETcorr_th15.Mod();
	if(threshold==4) return _stuff.myMETcorr_th20.Mod();
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}
double JetFilterModuleV2::GetMyMetXCorr(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.myMET_raw.Px(); // for consistency
	if(threshold==1) return _stuff.myMETcorr_th5.Px();
	if(threshold==2) return _stuff.myMETcorr_th10.Px();
	if(threshold==3) return _stuff.myMETcorr_th15.Px();
	if(threshold==4) return _stuff.myMETcorr_th20.Px();
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}
//___________________________________________________________________
double JetFilterModuleV2::GetMyMetYCorr(int cone, int threshold) 
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.myMET_raw.Py(); // for consistency
	if(threshold==1) return _stuff.myMETcorr_th5.Py();
	if(threshold==2) return _stuff.myMETcorr_th10.Py();
	if(threshold==3) return _stuff.myMETcorr_th15.Py();
	if(threshold==4) return _stuff.myMETcorr_th20.Py();
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}
//___________________________________________________________________
double JetFilterModuleV2::GetMyMetPhiCorr(int cone, int threshold)
{
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) return allstuff.myMET_raw.Phi(); // for consistency
	if(threshold==1) return _stuff.myMETcorr_th5.Phi();
	if(threshold==2) return _stuff.myMETcorr_th10.Phi();
	if(threshold==3) return _stuff.myMETcorr_th15.Phi();
	if(threshold==4) return _stuff.myMETcorr_th20.Phi();
	//if(threshold<0 || threshold>4) return -1.0E6;
	return -1.0E6;
}


//__________________________________________________________ setters of the
//__________________________________________________________ output params

//_________ need to specify which jets to return (variable "cone"): 
// 0--cone 0.4; 1--cone 0.7; 2--cone 1.0
void JetFilterModuleV2::SetMyJet_raw(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_raw.size()) _stuff.Jet_raw.at(i)=vec;
	else _stuff.Jet_raw.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJet_lev1(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev1.size()) _stuff.Jet_lev1.at(i)=vec;
	else _stuff.Jet_lev1.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJet_lev4(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev4.size()) _stuff.Jet_lev4.at(i)=vec;
	else _stuff.Jet_lev4.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJet_lev5(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev5.size()) _stuff.Jet_lev5.at(i)=vec;
	else _stuff.Jet_lev5.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJet_lev6(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev6.size()) _stuff.Jet_lev6.at(i)=vec;
	else _stuff.Jet_lev6.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNoEMobj_raw(int cone, int i, TLorentzVector vec)  {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_raw_noEMobj.size()) _stuff.Jet_raw_noEMobj.at(i)=vec;
	else _stuff.Jet_raw_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNoEMobj_lev1(int cone, int i, TLorentzVector vec)  {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev1_noEMobj.size()) _stuff.Jet_lev1_noEMobj.at(i)=vec;
	else _stuff.Jet_lev1_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNoEMobj_lev4(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev4_noEMobj.size()) _stuff.Jet_lev4_noEMobj.at(i)=vec;
	else _stuff.Jet_lev4_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNoEMobj_lev5(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev5_noEMobj.size()) _stuff.Jet_lev5_noEMobj.at(i)=vec;
	else _stuff.Jet_lev5_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNoEMobj_lev6(int cone, int i, TLorentzVector vec) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Jet_lev6_noEMobj.size()) _stuff.Jet_lev6_noEMobj.at(i)=vec;
	else _stuff.Jet_lev6_noEMobj.push_back(vec);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetEtaDet(int cone, int i, double param) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EtaDet.size()) _stuff.EtaDet.at(i)=param;
	else _stuff.EtaDet.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetEtaDetCorr(int cone, int i, double param) { // after removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EtaDetCorr.size()) _stuff.EtaDetCorr.at(i)=param;
	else _stuff.EtaDetCorr.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetEmFrRaw(int cone, int i, double param) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EmFrRaw.size()) _stuff.EmFrRaw.at(i)=param;
	else _stuff.EmFrRaw.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetEmFrCorr(int cone, int i, double param) { // before removing EM object
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.EmFrCorr.size()) _stuff.EmFrCorr.at(i)=param;
	else _stuff.EmFrCorr.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNtrk(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetNtrk.size()) _stuff.JetNtrk.at(i)=param;
	else _stuff.JetNtrk.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNtwr(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetNtwr.size()) _stuff.JetNtwr.at(i)=param;
	else _stuff.JetNtwr.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNobjMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nobj_match.size()) _stuff.Nobj_match.at(i)=param;
	else _stuff.Nobj_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNphoMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Npho_match.size()) _stuff.Npho_match.at(i)=param;
	else _stuff.Npho_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNeleMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nele_match.size()) _stuff.Nele_match.at(i)=param;
	else _stuff.Nele_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNmuMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nmu_match.size()) _stuff.Nmu_match.at(i)=param;
	else _stuff.Nmu_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNtauMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Ntau_match.size()) _stuff.Ntau_match.at(i)=param;
	else _stuff.Ntau_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetNbtagMatch(int cone, int i, int param) {
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.Nbtag_match.size()) _stuff.Nbtag_match.at(i)=param;
	else _stuff.Nbtag_match.push_back(param);
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyJetBlockInd(int cone, int i, int param) {
	// original jet index in JetBlock, need this after jet reordering
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(i>=0 && i<(int)_stuff.JetBlockInd.size()) _stuff.JetBlockInd.at(i)=param;
	else _stuff.JetBlockInd.push_back(param);
	return;
}

//___________________________________________________________________
//this is not used???28-12-2008,sam
void JetFilterModuleV2::SetMyNjet(int cone, int threshold, int param) {
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20, 5--et>25, 6--et>30, 7--et>35
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) _stuff.myNjet=param;
	if(threshold==1) _stuff.myNjet_th5=param;
	if(threshold==2) _stuff.myNjet_th10=param;
	if(threshold==3) _stuff.myNjet_th15=param;
	if(threshold==4) _stuff.myNjet_th20=param;
	if(threshold==5) _stuff.myNjet_th25=param;
	if(threshold==6) _stuff.myNjet_th30=param;
	if(threshold==7) _stuff.myNjet_th35=param;
	if(threshold<0 || threshold>7) return;
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMySumEtCorr(int cone, int threshold, double param) {
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) allstuff.mySumEt_raw=param; // for consistency
	if(threshold==1) _stuff.mySumEtCorr_th5=param;
	if(threshold==2) _stuff.mySumEtCorr_th10=param;
	if(threshold==3) _stuff.mySumEtCorr_th15=param;
	if(threshold==4) _stuff.mySumEtCorr_th20=param;
	if(threshold<0 || threshold>4) return;
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyMetCorr(int cone, int threshold, double px, double py) {
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) allstuff.myMET_raw.Set(px,py); // for consistency
	if(threshold==1) _stuff.myMETcorr_th5.Set(px,py);
	if(threshold==2) _stuff.myMETcorr_th10.Set(px,py);
	if(threshold==3) _stuff.myMETcorr_th15.Set(px,py);
	if(threshold==4) _stuff.myMETcorr_th20.Set(px,py);
	if(threshold<0 || threshold>4) return;
	return;
}
//___________________________________________________________________
void JetFilterModuleV2::SetMyMetCorr(int cone, int threshold, const TVector2& vec) {
	// thresholdcode= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
	JetStuff _stuff;
	if(cone==0) _stuff=jet04stuff;
	if(threshold==0) allstuff.myMET_raw.Set(vec); // for consistency
	if(threshold==1) _stuff.myMETcorr_th5.Set(vec);
	if(threshold==2) _stuff.myMETcorr_th10.Set(vec);
	if(threshold==3) _stuff.myMETcorr_th15.Set(vec);
	if(threshold==4) _stuff.myMETcorr_th20.Set(vec);
	if(threshold<0 || threshold>4) return;
	return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   My Cuts are defined here
//___________________________________________________________________
int JetFilterModuleV2::MyNjetCut(JetStuff jetstuff)
{
	if(jetstuff.myNjet<MinNjet || jetstuff.myNjet>MaxNjet) return 0;
	if(jetstuff.myNjet_th5<MinNjet5 || jetstuff.myNjet_th5>MaxNjet5) return 0;
	if(jetstuff.myNjet_th10<MinNjet10 || jetstuff.myNjet_th10>MaxNjet10) return 0;
	if(jetstuff.myNjet_th15<MinNjet15 || jetstuff.myNjet_th15>MaxNjet15) return 0;
	if(jetstuff.myNjet_th20<MinNjet20 || jetstuff.myNjet_th20>MaxNjet20) return 0;
	if(jetstuff.myNjet_th25<MinNjet25 || jetstuff.myNjet_th25>MaxNjet25) return 0;
	if(jetstuff.myNjet_th30<MinNjet30 || jetstuff.myNjet_th30>MaxNjet30) return 0;
	if(jetstuff.myNjet_th35<MinNjet35 || jetstuff.myNjet_th35>MaxNjet35) return 0;
	return 1;
}


//___________________________________________________________________
// cut on the number of smeared jets for pseudo experiments. We must cut on jets after we smear them.
// 12-28-2008, sam
//___________________________________________________________________
int JetFilterModuleV2::MySmeared_NjetCut(JetStuff jetstuff)
{
	if(jetstuff.smeared_Njet_th0<MinNjet || jetstuff.smeared_Njet_th0>MaxNjet) return 0;
	if(jetstuff.smeared_Njet_th5<MinNjet5 || jetstuff.smeared_Njet_th5>MaxNjet5) return 0;
	if(jetstuff.smeared_Njet_th10<MinNjet10 || jetstuff.smeared_Njet_th10>MaxNjet10) return 0;
	if(jetstuff.smeared_Njet_th15<MinNjet15 || jetstuff.smeared_Njet_th15>MaxNjet15) return 0;
	if(jetstuff.smeared_Njet_th20<MinNjet20 || jetstuff.smeared_Njet_th20>MaxNjet20) return 0;
	if(jetstuff.smeared_Njet_th25<MinNjet25 || jetstuff.smeared_Njet_th25>MaxNjet25) return 0;
	if(jetstuff.smeared_Njet_th30<MinNjet30 || jetstuff.smeared_Njet_th30>MaxNjet30) return 0;
	if(jetstuff.smeared_Njet_th35<MinNjet35 || jetstuff.smeared_Njet_th35>MaxNjet35) return 0;
	return 1;
}
//___________________________________________________________________
// sets the number of jets. We must cut on jets
// after we smear them. 12-28-2008, sam
//___________________________________________________________________
void JetFilterModuleV2::GetNjets(
			const std::vector<TLorentzVector>& vJet,	std::vector<int>& Njets)
{

	int Njet_th0=0, Njet_th5=0, Njet_th10=0, Njet_th15=0, Njet_th20=0, 
		 Njet_th25=0, Njet_th30=0, Njet_th35=0;

	for (unsigned int i=0; i < vJet.size(); ++i)
	{
		float fJetPt = vJet.at(i).Pt();
		if (fJetPt>0.0)  Njet_th0++;
		if (fJetPt>5.0)  Njet_th5++;
		if (fJetPt>10.0) Njet_th10++;
		if (fJetPt>15.0) Njet_th15++;
		if (fJetPt>20.0) Njet_th20++;
		if (fJetPt>25.0) Njet_th25++;
		if (fJetPt>30.0) Njet_th30++;
		if (fJetPt>35.0) Njet_th35++;
	}
	Njets.push_back(Njet_th0);
	Njets.push_back(Njet_th5);
	Njets.push_back(Njet_th10);
	Njets.push_back(Njet_th15);
	Njets.push_back(Njet_th20);
	Njets.push_back(Njet_th25);
	Njets.push_back(Njet_th30);
	Njets.push_back(Njet_th35);

	if (debug)
	{
		std::cout << __FUNCTION__ << "::" << __LINE__ 
			<< "::GetNjets njets 0,5,10,15,20,25,30,35=" << Njet_th0 << "," 
			<< Njet_th5 << "," << Njet_th10 << "," << Njet_th15 << "," << Njet_th20 << "," 
			<< Njet_th25 << "," << Njet_th30 << ", " << Njet_th35 << std::endl;
	}
}
//___________________________________________________________________
int JetFilterModuleV2::GetNjets(
			const std::vector<TLorentzVector>& vJet,	const float fEt)
{
// for a given set of jets and Et threshold
	int iNjets = 0; 

	for (unsigned int i=0; i < vJet.size(); ++i)
	{
		if (vJet.at(i).Pt()>=fEt) ++iNjets;
	}

	return iNjets;
}

//--------------------------------------------------------------------------
int JetFilterModuleV2::MyGenericCut(JetStuff jetstuff, double dz) {
	//__________________ cut on dZ=Zvx-Zjet
	if(fabs(dz)<MindZ || fabs(dz)>MaxdZ) return 0;
	//__________________ cut on 1st jet Et
	if(jetstuff.myNjet>0)
	{
		if(fabs(jetstuff.Jet_lev6_noEMobj[0].Eta())>MinJetEta
			&& fabs(jetstuff.Jet_lev6_noEMobj[0].Eta())<MaxJetEta
			&& (jetstuff.Jet_lev6_noEMobj[0].Pt()<MinEt1st 
			    || jetstuff.Jet_lev6_noEMobj[0].Pt()>MaxEt1st))
		{
			return 0;
		}
	}
	//__________________ cut on 2nd jet Et
	if(jetstuff.myNjet>1)
	{
		if(fabs(jetstuff.Jet_lev6_noEMobj[1].Eta())>MinJetEta
			&& fabs(jetstuff.Jet_lev6_noEMobj[1].Eta())<MaxJetEta
			&& (jetstuff.Jet_lev6_noEMobj[1].Pt()<MinEt2nd 
			    || jetstuff.Jet_lev6_noEMobj[1].Pt()>MaxEt2nd))
		{
			return 0;
		}
	}
	return 1;
}

//--------------------------------------------------------------------------
int JetFilterModuleV2::MyAngularCut(JetStuff jetstuff, CommonStuff miscstuff) {
	double dPhiMin=100.0;
	double dEtaMin=100.0;
	double dRMin=100.0;
	for(int i=0; i<jetstuff.myNjet; i++)
	{
		if(fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())>MinJetEta &&
				fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<MaxJetEta &&
				fabs(jetstuff.Jet_lev6_noEMobj.at(i).Pt())>MinEtThr &&
				fabs(jetstuff.Jet_lev6_noEMobj.at(i).Pt())<MaxEtThr)
		{
			// looping over all photons
			for(unsigned int j=0; j<miscstuff.myCorrPhoton.size(); j++)
			{
				double dPhi=fabs(TVector2::Phi_mpi_pi(jetstuff.Jet_lev6_noEMobj.at(i).Phi()-miscstuff.myCorrPhoton[j].Phi()));
				double dEta=fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta()-miscstuff.myCorrPhoton[j].Eta());
				double dR=sqrt(dPhi*dPhi+dEta*dEta);
				if(dPhi<dPhiMin) dPhiMin=dPhi;
				if(dEta<dEtaMin) dEtaMin=dEta;
				if(dR<dRMin) dRMin=dR;
				if(dPhi<MinDeltaPhigj || dPhi>MaxDeltaPhigj) return 0;
				if(dEta<MinDeltaEtagj || dEta>MaxDeltaEtagj) return 0;
				if(dR<MinDeltaRgj || dR>MaxDeltaRgj) return 0;
			}
		}
	}
	if(dPhiMin<MinDeltaPhigjMin || dPhiMin>MaxDeltaPhigjMin) return 0;
	if(dEtaMin<MinDeltaEtagjMin || dEtaMin>MaxDeltaEtagjMin) return 0;
	if(dRMin<MinDeltaRgjMin || dRMin>MaxDeltaRgjMin) return 0;
	return 1;
}

//check to make sure that there is no overlap of photons and electrons in the event
//--------------------------------------------------------------------------
int JetFilterModuleV2::MyDuplicateCut(CommonStuff miscstuff) {

	double dRMin=0.2;
	int passcode=0;
	if(fRemoveDuplicate!=0)
	{
		if(bad_EMjet_match_flag==0) return 1;
		for(unsigned int j=0; j<miscstuff.myCorrPhoton.size(); j++)
		{
			for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
			{
				if(fabs(miscstuff.myEleEtaDet.at(i))<1.2)
				{
					double dPhi=fabs(TVector2::Phi_mpi_pi(miscstuff.myElePhiDet.at(i)-miscstuff.myCorrPhoton[j].Phi()));
					double dEta=fabs(miscstuff.myEleEtaDet.at(i)-miscstuff.myPhoEtaDet[j]);
					double dR=sqrt(dPhi*dPhi+dEta*dEta);
					if(dR<dRMin) return 1;
				}
			}
		}
	}
	return passcode;
}



//----- Test version as of 04/04/07
//--------------------------------------------------------------------------
// ATTN!!! corrected jet & met are used in this cut
//
// This version of cuts can be used for generated MET in MetModel studies.
//
// This is where the large difference in MC and DATA in eta>0.95 and eta<1.1 etc
// is taken careof. Such events are thrown out.
int JetFilterModuleV2::MyMetCleanUpCut(JetStuff jetstuff, CommonStuff miscstuff, std::vector<TLorentzVector> jetvec, TVector2 metvec) {
	int Nbadmetjet=0;
	if(abs(fRemoveBadMet)>0) // new as of 12/14/07
	{
		double MetEM_dPhi_max=0.4;
		double MetPho_dPhi_max=0.3;
		double emXces_max_pho=18.5;
		double emXces_max_ele=13.0;
		for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++) // don't look at 1st photon for twr-9 effect
		{
			//__________________________ photon/MET phi separation
			double dPhi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(miscstuff.myCorrPhoton.at(i).Phi())-TVector2::Phi_0_2pi(metvec.Phi())));
			if(miscstuff.myCorrPhoton.size()>1)
			{
				//_____________________ if 2nd, 3rd.. photon is beyond eta>1.0  or eta<0.1  & too close to m
				if(i>0 && (fabs(miscstuff.myPhoEtaDet.at(i))>1.0 || fabs(miscstuff.myPhoEtaDet.at(i))<0.1) && dPhi<MetPho_dPhi_max) Nbadmetjet++;
				//______________________ if 2nd, 3rd photon
				if(i>0 && fabs(miscstuff.myPhoXces.at(i))>emXces_max_pho && dPhi<MetEM_dPhi_max) Nbadmetjet++;
			}
			else
			{
				if(fabs((miscstuff.myPhoEtaDet.at(i))>1.0 || fabs(miscstuff.myPhoEtaDet.at(i))<0.1) && dPhi<MetPho_dPhi_max) Nbadmetjet++;
				if(fabs(miscstuff.myPhoXces.at(i))>emXces_max_pho && dPhi<MetEM_dPhi_max) Nbadmetjet++;
			}
		}
		for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
		{
			double dPhi1=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(miscstuff.myCorrElectron.at(i).Phi())-TVector2::Phi_0_2pi(metvec.Phi())));
			double dPhi2=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(miscstuff.myElePhiDet.at(i))-TVector2::Phi_0_2pi(metvec.Phi())));
			double dPhi=(dPhi1<dPhi2) ? dPhi1 : dPhi2;
			if(miscstuff.myCorrElectron.size()>1)
			{
				if(i>0 && fabs(miscstuff.myEleXces.at(i))>emXces_max_ele && dPhi<MetEM_dPhi_max) Nbadmetjet++;
				if(i>0 && fabs(fabs(miscstuff.myEleEtaDet.at(i))-1.05)<0.05 && dPhi<MetEM_dPhi_max) Nbadmetjet++; // new as of 06/25/08
			}
			else
			{
				if(fabs(miscstuff.myEleXces.at(i))>emXces_max_ele && dPhi<MetEM_dPhi_max) Nbadmetjet++;
				if(fabs(fabs(miscstuff.myEleEtaDet.at(i))-1.05)<0.05 && dPhi<MetEM_dPhi_max) Nbadmetjet++; // new as of 06/25/08
			}
		}

		//       //---------------------- This part was added on 07/03/08, to remove events where a photon got lost in CEM cracks
		//       double jet_minPt=5.0;
		//       double jet_maxEta=1.3;
		//       double jet_maxEmFr=0.875;
		//       double MetJet_dPhi_max=0.4;
		//       for(int i=0; i<jetstuff.Jet_lev6_noEMobj.size(); i++)
		// 	{
		// 	  if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>jet_minPt &&
		// 	     fabs(jetstuff.EtaDet.at(i))<jet_maxEta &&
		// 	     jetstuff.Nobj_match.at(i)==0 &&
		// 	     jetstuff.EmFrRaw.at(i)>jet_maxEmFr)
		// 	    {
		// 	      double dPhi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj.at(i).Phi())-TVector2::Phi_0_2pi(metvec.Phi())));
		// 	      if(dPhi<MetJet_dPhi_max) Nbadmetjet++;
		// 	    }
		// 	}

		//---------------------- This part was added on 07/03/08, to remove events where a photon got lost in CEM cracks
		double jet_minPt=5.0;
		//double jet_maxEta=1.3;
		double jet_maxEmFr1=0.875;
		double jet_maxEmFr2=0.3;
		double MetJet_dPhi_max=0.3;
		int jet_NtwrMax=10;
		int jet_NtrkMax=5;

		double dPhi_min=TMath::Pi();
		double eta_det=10.0;
		double emfr=0.5;
		int ntwr=100;
		int ntrk=100;
		for(unsigned int i=0; i<jetstuff.Jet_lev6_noEMobj.size(); i++)
		{
			if(jetstuff.Jet_lev6_noEMobj.at(i).Pt()>jet_minPt &&
					fabs(jetstuff.EtaDet.at(i))<1.3 &&
					jetstuff.Nobj_match.at(i)==0)
			{
				double dPhi=fabs(TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(jetstuff.Jet_lev6_noEMobj.at(i).Phi())-TVector2::Phi_0_2pi(metvec.Phi())));
				if(dPhi<dPhi_min)
				{
					dPhi_min=dPhi;
					eta_det=jetstuff.EtaDet.at(i);
					emfr=jetstuff.EmFrRaw.at(i);
					ntwr=jetstuff.JetNtwr.at(i);
					ntrk=jetstuff.JetNtrk.at(i);
				}
			}
		}
		if(dPhi_min<MetJet_dPhi_max
			&& emfr>jet_maxEmFr1 
			&& ntwr<jet_NtwrMax 
			&& ntrk<jet_NtrkMax 
			&& fabs(eta_det)<1.1) 
		{
			Nbadmetjet++;
		}
		if(dPhi_min<MetJet_dPhi_max
			&& emfr<jet_maxEmFr2 
			&& ntwr<jet_NtwrMax 
			&& ntrk<jet_NtrkMax
			&& (fabs(eta_det)<0.1 || (fabs(eta_det)>1.1 && fabs(eta_det)<1.2)))
		{
			Nbadmetjet++;
		}
	}
	return Nbadmetjet;
}


//________________________________ returns Qt of all jets above Eta & Et thresholds
double JetFilterModuleV2::MyQtJet(JetStuff jetstuff) {
	double qt=-1.0E6;
	TLorentzVector vec(0.0,0.0,0.0,0.0);
	for(int i=0; i<jetstuff.myNjet; i++)
	{
		if(fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())>MinJetEta &&
				fabs(jetstuff.Jet_lev6_noEMobj.at(i).Eta())<MaxJetEta &&
				fabs(jetstuff.Jet_lev6_noEMobj.at(i).Pt())>MinEtThr &&
				fabs(jetstuff.Jet_lev6_noEMobj.at(i).Pt())<MaxEtThr)
		{
			vec=vec+jetstuff.Jet_lev6_noEMobj.at(i);
		}
	}
	if(vec.E()>0.0) qt=vec.Pt();
	return qt;
}
//________________________________ returns SumEt=jets+photons+leptons+taus+met (jets above Eta & Et thresholds)
//                                 metscenario= 0--raw, 1--et>5, 2--et>10, 3--et>15, 4--et>20
double JetFilterModuleV2::MyHtAll(const std::vector<TLorentzVector> vJets, 
									JetStuff jetstuff, CommonStuff miscstuff, 
									int metscenario)
{
	double ht=0.0;
	//_________________ contribution from jets
	for(int i=0; i<jetstuff.myNjet; i++)
	{
		if(fabs(vJets.at(i).Eta())>MinJetEta
			&& fabs(vJets.at(i).Eta())<MaxJetEta
			&& fabs(vJets.at(i).Pt())>MinEtThr
			&& fabs(vJets.at(i).Pt())<MaxEtThr)
		{
			ht=ht+vJets.at(i).Pt();
		}
	}
	//_________________ contribution from photons
	for(unsigned int i=0; i<miscstuff.myCorrPhoton.size(); i++)
	{
		ht=ht+miscstuff.myCorrPhoton.at(i).Pt();
	}
	//_________________ contribution from electrons
	for(unsigned int i=0; i<miscstuff.myCorrElectron.size(); i++)
	{
		ht=ht+miscstuff.myCorrElectron.at(i).Pt();
	}
	//_________________ contribution from muons
	for(unsigned int i=0; i<miscstuff.myCorrMuon.size(); i++)
	{
		ht=ht+miscstuff.myCorrMuon.at(i).Pt();
	}
	//_________________ contribution from taus
	for(unsigned int i=0; i<miscstuff.myCorrTau.size(); i++)
	{
		ht=ht+miscstuff.myCorrTau.at(i).Pt();
	}
	//_________________ contribution from MET
	if(metscenario==0) ht=ht+miscstuff.myMET_raw.Mod();
	if(metscenario==1) ht=ht+jetstuff.myMETcorr_th5.Mod();
	if(metscenario==2) ht=ht+jetstuff.myMETcorr_th10.Mod();
	if(metscenario==3) ht=ht+jetstuff.myMETcorr_th15.Mod();
	if(metscenario==4) ht=ht+jetstuff.myMETcorr_th20.Mod();
	if(metscenario<0 || metscenario>4) ht=ht+miscstuff.myMET_raw.Mod();
	return ht;
}

//_________________________________________ dump jet info for debugging
void JetFilterModuleV2::DumpJets(const std::string func, const unsigned int line,
								const JetStuff &jetstuff, const int lev)
{
	// lev: 0=raw, 1=lev1, 4=lev4, 5=lev5, 6=lev6, 7=lev7
	// lev: 10=raw noEM, 11=lev1 noEM, 14=lev4 noEM, 15=lev5 noEM, 16=lev6 noEM, 17=lev7 noEM
	// lev 20 = newJetLev6, 21 = smeared_newJetLev6_noEMObj
	
	if (jetstuff.myNjet>0)
	{
		std::string sLev;
		switch (lev)
		{
			case 0: { sLev = "RAW"; break;}
			case 1: { sLev = "LEV1"; break;}
			case 4: { sLev = "LEV4"; break;}
			case 5: { sLev = "LEV5"; break;}
			case 6: { sLev = "LEV6"; break;}
			case 10: { sLev = "RAW NO EM"; break;}
			case 11: { sLev = "LEV1 NO EM"; break;}
			case 14: { sLev = "LEV4 NO EM"; break;}
			case 15: { sLev = "LEV5 NO EM"; break;}
			case 16: { sLev = "LEV6 NO EM"; break;}
			case 20: { sLev = "newJetL6 NO EM"; break;}
			case 21: { sLev = "smeared_newJetL6 NO EM"; break;}
			default: { sLev = "UNKOWN!"; }
		};
		
		std::cout << "=== Caller,line::" <<func << ":" << line <<": Jets For of " << lev  << " - " << sLev << std::endl;

		for (int i=0; i<jetstuff.myNjet; i++)
		{
			float jetEt = 0.0;
			float eta = -99.99;
			switch (lev)
			{
				case 6:
					{
						jetEt = jetstuff.Jet_lev6.at(i).Pt();
						eta   = jetstuff.Jet_lev6.at(i).Eta();
						break;
					}
				
				case 16:
					{
						jetEt = jetstuff.Jet_lev6_noEMobj.at(i).Pt();
						eta   = jetstuff.Jet_lev6_noEMobj.at(i).Eta();
						break;
					}
				case 20:
					{
						jetEt = jetstuff.newJetLev6.at(i).Pt();
						eta   = jetstuff.newJetLev6.at(i).Eta();
						break;
					}
				case 21:
					{
						jetEt = jetstuff.smeared_newJetLev6_noEMObj.at(i).Pt();
						eta   = jetstuff.smeared_newJetLev6_noEMObj.at(i).Eta();
						break;
					}
				default:
					{
						std::cout << "function not defined for lev= " << lev << " jets!\n";
						return;
					}

			}
			std::cout << "i = " << i << " jetPt= " << jetEt << "\t eta= " << eta << std::endl;
		}
		
	} else {
		std::cout << "=== Caller,line::" <<func << ":" << line <<": NO Jets For lev " << lev << std::endl;
	}
	
}
//_________________________________________ dump smeared jet info for debugging
void JetFilterModuleV2::DumpSmearedJets(const std::string func, const unsigned int line,
								const JetStuff &jetstuff, const int lev)
{
	std::cout << "=== Caller,line::" <<func << ":" << line <<": SMEARED Jets For lev " << lev << std::endl;
	unsigned int njet =  jetstuff.smeared_newJetLev6_noEMObj.size();
	for (unsigned int i=0; i<njet; i++)
	{
		float jetEt = 0;
		float eta = -99;

		switch (lev)
		{
			case 6:
			{
				jetEt = jetstuff.smeared_newJetLev6_noEMObj.at(i).Pt();
				eta   = GetMyJetEtaDetCorr(0,i);
				break;
			}
			default:
			{
				std::cout << "function not defined for lev= " << lev << " jets!\n";
				return;
			}

		}
		std::cout << "i = " << i << " jetPt= " << jetEt << "\t eta= " << eta << std::endl;
	}
	std::cout << "===================" << std::endl;
}
//_________________________________________ dump EM objects info for debugging
// prints all EM Objects by default
//------------------------------------------------------------------------------------
void JetFilterModuleV2::DumpEMobjects(const std::string func, const unsigned int line,
								const CommonStuff& commonStuff, const int type)
{
	//type 0=all (pho+ele+++), 1 = pho only, 2= ele only
	if (type == 0 || type == 1)
	{
		if (commonStuff.myCorrPhoton.size())
		{
			std::cout << __FUNCTION__ << " Caller,line::" <<func << ":" << line <<": Photons" << std::endl;
			std::cout << setw(10) << "Index" << setw(10) << "Pt" << setw(10) << "Eta" << setw(10) << "Phi" << std::endl;
			for (unsigned int i=0; i<commonStuff.myCorrPhoton.size(); ++i)
			{
				std::cout << setw(10) << i 
							<< setw(10) << commonStuff.myCorrPhoton.at(i).Pt() 
							<< setw(10) << commonStuff.myPhoEtaDet.at(i) 
							<< setw(10) << commonStuff.myCorrPhoton.at(i).Phi() 
							<< std::endl;
			}
		}
	}
	
	if (type == 0 || type == 2)
	{
		if (commonStuff.myCorrElectron.size())
		{
			std::cout << "=== Caller,line::" <<func << ":" << line <<": Electrons" << std::endl;
			std::cout << setw(10) << "Index" << setw(10) << "Pt" << setw(10) << "Eta" << setw(10) << "Phi" << std::endl;
			for (unsigned int i=0; i<commonStuff.myCorrElectron.size(); i++)
			{
				std::cout << setw(10) << i 
					<< setw(10) << commonStuff.myCorrElectron.at(i).Pt() 
					<< setw(10) << commonStuff.myEleEtaDet.at(i) 
					<< setw(10) << commonStuff.myCorrElectron.at(i).Phi() 
					<< std::endl;
			}
		}
	}
	
}


//-----------------------------------------------------------------------------
void JetFilterModuleV2::JetMismeasureProb(JetStuff& jetstuff, const TVector2& myMetVec, const int systcode)
{
	// This is part of the MEt addition to the closest jet for the pseudo experiments.
	// This will use the L6NoEmObj jets and calculate the probability for each jet 
	// for mismeasurement of its energy.

	jetstuff.fL6NoEm_MisMeasureProb.clear();

	if (Npoints>0) // case of MET from pseudo-experiments
	{
		for(unsigned int i=0; i<jetstuff.Jet_lev6_noEMobj.size(); i++)
		{
			float test_prob = 0;
			jetstuff.fL6NoEm_MisMeasureProb.push_back(test_prob);  // stuff with zeros for now
			std::cout << __FUNCTION__ << ":: Jet Index = " <<  i << std::endl;
			std::cout 	<< "\t L6JetPt    = jetstuff.Jet_lev6_noEMobj.at(i).Pt() = " 
							<< jetstuff.Jet_lev6_noEMobj.at(i).Pt() 
							<< "\t\nrawJetPt  = jetstuff.Jet_raw_noEMobj.at(i).Pt()  = " 
							<< jetstuff.Jet_raw_noEMobj.at(i).Pt()
							<< "\t\nJetEtaDet = jetstuff.EtaDet.at(i)                = " 
							<< jetstuff.EtaDet.at(i) << std::endl;

			if (jetstuff.Jet_lev6_noEMobj.at(i).Pt() > 3.0 
				&& jetstuff.Jet_raw_noEMobj.at(i).Pt() > 3.0 
				&& fabs(jetstuff.EtaDet.at(i)) <= MaxJetEta)					//why aren't we using the CorrectedEtaDetector???
			{
				std::cout << "\t passed L6JetPt>3 && rawJetPt>3 and fabs(DetEta)<MaxJetEta(==3)" << std::endl;
				double dPhi_jetmet = fabs(TVector2::Phi_mpi_pi(jetstuff.Jet_lev6_noEMobj.at(i).Phi()-myMetVec.Phi()));
				double trueEjet = jetstuff.Jet_lev6_noEMobj.at(i).E();
				double scaleE   = jetstuff.Jet_lev6_noEMobj.at(i).E()/jetstuff.Jet_lev6_noEMobj.at(i).Pt();
				double cos_dphi = cos(dPhi_jetmet);
				std::cout<< "\t dPhi_jetmet = " << dPhi_jetmet
							<< " trueEjet = " << trueEjet
							<< " scaleE = " << scaleE
							<< " cos_dphi " << cos_dphi << std::endl;

				if (fabs(cos_dphi)>0.0) // only consider these cases (no contribution is dPhi_jetmet=pi/2)
				{
					std::cout << "\t passed fabs(cos_dphi)>0" << std::endl;
					trueEjet = trueEjet + scaleE * myMetVec.Mod() / cos_dphi;
					std::cout << "\t trueEjet = trueEjet + scaleE * myMetVec.Mod() / cos_dphi =" 
								 << trueEjet << " + " << " + " << scaleE << " + " 
								 << myMetVec.Mod() / cos_dphi << " = "  
								 << trueEjet << std::endl;

					if (trueEjet>3.0 && trueEjet<1000.0)
					{
						std::cout << "\t passed trueEjet>3 && trueEjet<1000 " << std::endl;
						double parJER[5];
						parJER[0]=MyJer_meanG(trueEjet,jetstuff.EtaDetCorr.at(i),0,systcode,JERForGaussMean()); // for now systcode is used as jer_stat_code
						parJER[1]=MyJer_sigmaG(trueEjet,jetstuff.EtaDetCorr.at(i),0,systcode,JERForGaussSigma());
						parJER[2]=MyJer_mpvL(trueEjet,jetstuff.EtaDetCorr.at(i),0,systcode,JERForLandauMPV());
						parJER[3]=MyJer_sigmaL(trueEjet,jetstuff.EtaDetCorr.at(i),0,systcode,JERForLandauSigma());
						parJER[4]=MyJer_normG(trueEjet,jetstuff.EtaDetCorr.at(i),0,systcode,JERForGaussNorm());
						std::cout << "\t JER Par0,1,2,3,4 = " << parJER[0] << ", " << parJER[1] 
							<< ", "<< parJER[2] << ", "<< parJER[3] << ", "<< parJER[4] << std::endl;

						std::cout << "\t trueEjet = " << trueEjet 
							<< " EtaDetCorr = " << jetstuff.EtaDetCorr.at(i)
							<< " systcode = " << systcode << std::endl;
						double jer_limit=MyIntegralUpLimit(trueEjet,jetstuff.EtaDetCorr.at(i),0,systcode);
						std::cout << "\t jet_limit = MyIntegralUpLimit(trueEjet,jetstuff.EtaDetCorr.at(i),0,systcode) = " 
							<< jer_limit << std::endl;
						TF1* JerFun=new TF1("jer",MyJER,-1.0,jer_limit,5);
						for(int par_ind=0; par_ind<5; par_ind++)
						{
							JerFun->SetParameter(par_ind,parJER[par_ind]);
						}
						double off_set=JerFun->GetMaximumX(-1.0,1.0); // returns mpv of JER function;		  
						if (ZeroJERoffset()) off_set = 0;
						std::cout << "\t trueEjet " << trueEjet << " L6JetE = "  << jetstuff.Jet_lev6_noEMobj.at(i).E() << std::endl;
						if(trueEjet>(jetstuff.Jet_lev6_noEMobj.at(i).E())) // case-1: MET points along the jet
						{
							std::cout << "\tif(trueEjet>(jetstuff.Jet_lev6_noEMobj.at(i).E())) :: MEt along jet direction" << std::endl;	
							double x1=(jetstuff.Jet_lev6_noEMobj.at(i).E())/trueEjet-1.0+off_set;
							std::cout << "\t\t x1=(jetstuff.Jet_lev6_noEMobj.at(i).E())/trueEjet-1.0+off_set = "  
								<< x1 << std::endl;
							if(x1<-1.0) x1=-1.0;
							std::cout << "\t\tif(x1<-1.0) x1=-1.0 --> x1=" << x1 << std::endl;
							if(x1>jer_limit) x1=jer_limit;
							std::cout << "\t\t if(x1>jer_limit) x1=jer_limit --> x1=" << x1 << std::endl;
							test_prob=JerFun->Integral(-1.0,x1);
							std::cout << "\t\ttest_prob = JerFun->Integral(-1.0,x1) = " << test_prob << std::endl; 
						}
						else // case-2: MET points opposite to the jet
						{
							std::cout << "\tif(trueEjet<(jetstuff.Jet_lev6_noEMobj.at(i).E())) :: MEt opposite to jet direction"  << std::endl;	
							double x1=(jetstuff.Jet_lev6_noEMobj.at(i).E())/trueEjet-1.0+off_set;
							std::cout << "\t\t x1=(jetstuff.Jet_lev6_noEMobj.at(i).E())/trueEjet-1.0+off_set = "  
								<< x1 << std::endl;
							if(x1<-1.0) x1=-1.0;
							std::cout << "\t\tif(x1<-1.0) x1=-1.0 --> x1=" << x1 << std::endl;
							if(x1>jer_limit) x1=jer_limit;
							std::cout << "\t\t if(x1>jer_limit) x1=jer_limit --> x1=" << x1 << std::endl;
							test_prob=JerFun->Integral(x1,jer_limit);
							std::cout << "\t\ttest_prob = JerFun->Integral(x1,jer_limit) = " << test_prob << std::endl; 
						}
						delete JerFun;
					} else { 
						std::cout << "\t FAILED trueEjet>3 && trueEjet<1000 " << std::endl;
					}
				} else {
					std::cout << "\t FAILED fabs(cos_dphi)>0" << std::endl;

				}
			} else {
				std::cout << "\t FAILED L6JetPt>3 && rawJetPt>3 and fabs(DetEta)<MaxJetEta(==3)" << std::endl;
			}


			std::cout << __FUNCTION__ << ":: i, prob = "  << i << "\t" <<  test_prob << std::endl;
			jetstuff.fL6NoEm_MisMeasureProb.at(i) = test_prob;
		}
	}

	//std::cout << __FUNCTION__ << ":: mismeasure size = " << jetstuff.fL6NoEm_MisMeasureProb.size() <<std::endl;

}

//-----------------------------------------------------------------------------
void JetFilterModuleV2::Dump_newJetLev6(const std::string funcname, const int line, const JetStuff& jetstuff)
{
	// Dumps the newJetLev6 vector in JetStuff 
	std::cout << __FUNCTION__ << ":: called By " << funcname << " Line " << ToStr(line) << std::endl;
	
	if (jetstuff.newJetLev6.size())
	{
		std::cout << "======================= NEWJETLEV6  " << std::endl;
		std::cout << "Index  JetPt   JetEta  JetPhi" << std::endl;
		for(unsigned int j=0; j<jetstuff.newJetLev6.size(); j++)
		{
			std::cout <<  setw(3) << j << setw(10) << jetstuff.newJetLev6.at(j).Pt() 
						<<  setw(10) << jetstuff.newJetLev6.at(j).Eta() 
						<< setw(10) << jetstuff.newJetLev6.at(j).Phi() << std::endl;
		}
	}
}

// find the jet closest to the MEt
/*int JetFilterModuleV2::ClosestJet(const std::vector<TLorentzVector> vJets, const TVector2 vMet)
{

		double phi_met    = jetstuff.myMETcorr_th15.Phi();
		for (unsigned int i=0; i < jetstuff.Jet_lev6_noEMobj.size(); i++) // packing "jetstuff.newJetLev6" with "old" jets
		{
			jetstuff.newJetLev6.push_back(jetstuff.Jet_lev6_noEMobj.at(i));

			if (jetstuff.myMETcorr_th15.Mod() > 0.0)
			{
				double phi_jet = jetstuff.Jet_lev6_noEMobj.at(i).Phi();
				double dphi    = fabs(TVector2::Phi_mpi_pi(phi_jet-phi_met));
				if (dphi< closejet_delphi)
				{
					closejet_indJ = i;
					closejet_delphi = dphi;
				}
			}

		}
}
*/


// find the least mismeasure jet for any given set of jets 
int JetFilterModuleV2::LeastMismeasuredJet(const std::vector<TLorentzVector> vJets, const TVector2 MetVec)
{
	
	const int systcode = 0;			//for now this should be good
	int iJetIndex = -1;				//default, if no suitable jet is found
	float fJetProb = 0; 

	if (Npoints>0) // case of MET from pseudo-experiments
	{
		for(unsigned int i=0; i<vJets.size(); i++)
		{
			float test_prob = 0;

			//sasha cut on det eta and uses det eta when choosing JER. I am using 4-vec Eta. 
			// I made this routine more generalized and hence I may be choosing slighly 
			// different JERs need to double check as this is only for the MEt addition 10-05-2009
			if (vJets.at(i).Pt() > 3.0 
				&& fabs(vJets.at(i).Eta()) <= MaxJetEta)					//why aren't we using the CorrectedEtaDetector???
			{
				double dPhi_jetmet = fabs(TVector2::Phi_mpi_pi(vJets.at(i).Phi()-MetVec.Phi()));
				double trueEjet = vJets.at(i).E();
				double scaleE   = vJets.at(i).E()/vJets.at(i).Pt();
				double cos_dphi = cos(dPhi_jetmet);

				if (fabs(cos_dphi)>0.0) // only consider these cases (no contribution is dPhi_jetmet=pi/2)
				{
					trueEjet = trueEjet + scaleE * MetVec.Mod() / cos_dphi;

					if (trueEjet>3.0 && trueEjet<1000.0)
					{
						double parJER[5];
						parJER[0]=MyJer_meanG(trueEjet,vJets.at(i).Eta(),0,systcode, JERForGaussMean()); // for now systcode is used as jer_stat_code
						parJER[1]=MyJer_sigmaG(trueEjet,vJets.at(i).Eta(),0,systcode, JERForGaussSigma());
						parJER[2]=MyJer_mpvL(trueEjet,vJets.at(i).Eta(),0,systcode, JERForLandauMPV());
						parJER[3]=MyJer_sigmaL(trueEjet,vJets.at(i).Eta(),0,systcode, JERForLandauSigma());
						parJER[4]=MyJer_normG(trueEjet,vJets.at(i).Eta(),0,systcode, JERForGaussNorm());

						double jer_limit=MyIntegralUpLimit(trueEjet,vJets.at(i).Eta(),0,systcode);
						TF1* JerFun=new TF1("jer",MyJER,-1.0,jer_limit,5);
						for(int par_ind=0; par_ind<5; par_ind++)
						{
							JerFun->SetParameter(par_ind,parJER[par_ind]);
						}

						double off_set=JerFun->GetMaximumX(-1.0,1.0); // returns mpv of JER function;		  
						if (ZeroJERoffset()) off_set = 0;
						
						if(trueEjet>(vJets.at(i).E())) // case-1: MET points along the jet
						{
							double x1=(vJets.at(i).E())/trueEjet-1.0+off_set;
							if(x1<-1.0) x1=-1.0;
							if(x1>jer_limit) x1=jer_limit;
							test_prob=JerFun->Integral(-1.0,x1);
						}
						else // case-2: MET points opposite to the jet
						{
							double x1=(vJets.at(i).E())/trueEjet-1.0+off_set;
							if(x1<-1.0) x1=-1.0;
							if(x1>jer_limit) x1=jer_limit;
							test_prob=JerFun->Integral(x1,jer_limit);
						}
						delete JerFun;

						if (test_prob> fJetProb)
						{
							fJetProb = test_prob;
							iJetIndex = i;
						}
						
					}
				}
			}

		}
	}

	return iJetIndex;
}

//=============================================================================
TLorentzVector JetFilterModuleV2::FindMatchingHEPGPar(const TLorentzVector tlObj,
									const int iPDGcode, const int iStatus,const float fDelR)
{
	// For a given detector obj (4-vec) this will find the matching HEPG particle
	// and returns its 4-vec
	// iStatus must be 1,2, or 3. see elog #1201

	if (iStatus<1 || iStatus>3)
	{
		StdOut(__FILE__,__LINE__,2,
		"WARNING! Requiring a unknown HEPG particle status!");
	}
	
	TLorentzVector hepgVec(0,0,0,0);
	
	if (fHeaderBlock->McFlag())
	{
		int Nparticles = fGenpBlock->NParticles();

		TGenParticle *par, *mom;

		// I'll search only the first 50 particles
		// That should be enough to find the leading decay particles
		for (int i = 0 ; (i < Nparticles || i<25) ; i++)	//_____ not much diff from running over all. my match is mostly within the fist 10 objects
		{
			par = fGenpBlock->Particle(i);
			int im = par->GetFirstMother();

			if (im >= 0) {	//_________________________________im=-1 means no mother for imcoming particles
				mom = fGenpBlock->Particle(im);

				if (mom != 0) {
					int par_id   = par->GetPdgCode();
					int par_stat = par->GetStatusCode();
					TLorentzVector parvec;
					par->Momentum(parvec);

					if (par_stat != iStatus) continue; 			// find a stable particle first. see elog#1201
					if (par_id != iPDGcode) continue;
						
					if (tlObj.DeltaR(parvec) < fDelR) {
						//std::cout << "FOUND MATCH "<< std::endl;	
						//std::cout << "i,par_id,pat_stat,delR,PtRat = " << i << "\t" << par_id << "\t" << par_stat;
						//std::cout << "\t" << tlObj.DeltaR(parvec) << "\t" << tlObj.Pt()/parvec.Pt() << std::endl;
						hepgVec = parvec;
						break;
					}
					
				}
			} 
		}

	} // if(qMc)
	
	if (hepgVec.Pt() == 0 && PrintLevel() > PRN_WARN)
	{
		std::cout << __FILE__ << ":" << __FUNCTION__ << ":" << __LINE__ 
					<< ": NO MATCH FOR PHOTON FOUND! Run,Event="
					<< fHeaderBlock->RunNumber() << ", " 
					<< fHeaderBlock->EventNumber() << std::endl;
	}

	return hepgVec;
}

//-----------------------------------------------------------------------------
TLorentzVector JetFilterModuleV2::FindMatchingHADJet(const TLorentzVector tlJetVec,
		const float fDelR) 
{

	TLorentzVector tlHadJetVec(0,0,0,0);

	if (fHeaderBlock->McFlag())
	{
		for (int i=0; i < fHadJetBlockClu04->NJets(); ++i)  		//loop over HAD jets to find match
		{
			TStnJet *hadjet = fHadJetBlockClu04->Jet(i); 
			TLorentzVector hjetvec(*hadjet->Momentum());

			if (hjetvec.DeltaR(tlJetVec) < fDelR)
			{
				tlHadJetVec = hjetvec;
				break;
			}
		}
	} else {
		StdOut(__FUNCTION__,__LINE__,3,"This is not MC. No HAD info. returning 0!"); 
	}
	return tlHadJetVec;
}


//-------------------------------------------------------------------
void JetFilterModuleV2::DumpTLvec(const std::string func, const unsigned int line,
								const std::vector<TLorentzVector> vObj , const std::string label)
{
	// Dumps the info about the objets(4-vectors) in a vector.

	std::cout << __FUNCTION__ <<" Caller-> " << func << " line " << line << std::endl;
	if (label.length())
	{
		std::cout << "\t=== " << label << " === " << std::endl;
	}

	if (vObj.size())
	{
		int width = 15;
		std::cout << setw(width) << "Index" << setw(width) << "Pt" << setw(width) << "Eta" << setw(width) << "Phi" << std::endl; 
		
		for (unsigned int i=0; i < vObj.size() ; ++i)
		{
			std::cout << setprecision(10) << setw(width) << i << setw(width) << vObj.at(i).Pt()
					<< setw(width) << vObj.at(i).Eta() << setw(width) << vObj.at(i).Phi() << std::endl;

		}

	}


}  // DumpTLvec

//-------------------------------------------------------------------
int JetFilterModuleV2::JetClosestToMet(std::vector<TLorentzVector> vJet,
				const TVector2 Met)
{
	//Returns the index of the Jet closest to the MEt direction if 
	// DelPhi(Jet,MEt)<0.4
	// returns >=0 if a jet is found
	// returns -1 else

	int iJetIndex = -1;
	int iNjets=0;		//number of jets within DelPhi<0.4
	if (Met.Mod() > 0.0)
	{
		float fPhiMet = Met.Phi();
		float fDphi = 999.99;
		
		//	std::cout << setprecision(4) << setw(10) << "fDphi" 
		//				<< setw(8) << "iJetInd" << setw(8) << "fdphi" 
		//				<< setw(8) << "Pt" <<std::endl;
		
		for (unsigned i=0; i<vJet.size(); ++i)
		{
			// if this is a leftover after removing EM objects ignore it.
			if (vJet.at(i).Pt()<3) continue;	//Will add MEt to a jet with Pt>3GeV

			float fPhiJet = vJet.at(i).Phi();
			float fdphi    = fabs(TVector2::Phi_mpi_pi(fPhiJet-fPhiMet));
			
			//std::cout << setw(2) << i << setw(8) << fDphi << setw(8) << iJetIndex 
			//				<< setw(8) << fdphi << setw(8) << vJet.at(i).Pt() <<std::endl;
			
			if(fdphi<0.4 || (TMath::Pi()-fdphi)<0.4) 
			{
				fDphi = fdphi;
				iJetIndex = i;
				iNjets++;
			}
		}
	}

	//if more than one jet is within DelPhi<0.4
	// do not add MEt to any of them.
	if (iNjets>1) iJetIndex = -1;
	//std::cout <<  "->>>>>>>>>>>>>> choosen = " << iJetIndex <<std::endl;

	return iJetIndex;
}

//---------------------------------------------------------------------------
// Returns original jet block index of the jet matched to an EM object
//---------------------------------------------------------------------------
int JetFilterModuleV2::GetMatch_JetIndex(int ObjType, int EmInd) const
{
	if (sam_matchstuff.size() >0) {
		for (unsigned int j=0; j < sam_matchstuff.size(); j++) {
			if (sam_matchstuff[j].EmObjType_match == ObjType) {
				if (sam_matchstuff[j].EmInd_match == EmInd ) {
					return sam_matchstuff[j].JetInd_match;
				}
			}
		}
		std::cerr << __FILE__ << "::" << __LINE__ << "::" << __FUNCTION__
			<< ":: ERROR::: EM obj did not match" <<std::endl;
		std::cout << "\tERROR::Run,Evt:: " << fHeaderBlock->RunNumber() << ", " << fHeaderBlock->EventNumber() << std::endl;
		//	exit (1);
	}

	std::cerr << __FILE__ << "::" << __LINE__ << "::" << __FUNCTION__
		<< ":: ERROR:: NO MATCHING INFO FOUND. Returining -1 ! " <<std::endl;
	std::cout << "\tERROR::Run,Evt:: " << fHeaderBlock->RunNumber() << ", " << fHeaderBlock->EventNumber() << std::endl;
	return -1;
}

