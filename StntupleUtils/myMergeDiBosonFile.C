#include <iomanip>
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TFolder.h"
#include <Stntuple/Stntuple/val/stntuple_val_functions.hh>

/* I think I got this from Sasha
 * this is suppose to merge hists split
 * among many files when run on CAF.
 * I tried to get it to work once
 * without a success.
 */
void myMergeDiBosonFile()
{
  gROOT->Reset();
  char fileList[500];

	char	f1[100]="/mnt/autofs/misc/nbay02.a/samantha/stn_rel/MetTest_PhotonMC_1.root";
	char	f2[100]="/mnt/autofs/misc/nbay02.a/samantha/stn_rel/MetTest_PhotonMC_1.root";

  //sprintf(fileList,"%s %s %s",f1,f2,f3);
  sprintf(fileList,"%s %s",f1,f2);

  merge_stn_hist(fileList,"Merged.root");
  //merge_stn_hist(fileList,"RESULTS/MetTest/NewStuff/04132009_MetAddedJetMisMeasureProb/InclNjet15_MetSig2/Merged.root");

  return;
}
